//! Chebyshev binomial-differential integration and fractional-power
//! rationalization.
//!
//! Two closely related substitution mechanisms for algebraic integrands:
//!
//! - **Binomial differential** `x^m · (a + b·x^n)^p` with `m, n, p ∈ ℚ` and
//!   `a, b` constant w.r.t. the integration variable. Chebyshev's theorem
//!   says the integral is elementary exactly when one of three conditions
//!   holds (`s` = denominator of `p` in lowest terms):
//!   1. `p ∈ ℤ`: with `s = lcm(den m, den n)`, `x = t^s` rationalizes.
//!   2. `(m+1)/n ∈ ℤ`, `s > 1`: `t^s = a + b·x^n` rationalizes.
//!   3. `(m+1)/n + p ∈ ℤ`, `s > 1`: `t^s = a·x^(−n) + b` rationalizes.
//! - **Fractional-power rationalization**: integrands that are rational in
//!   `x` and in fractional powers `g^(k/s)` of a single linear base
//!   `g = d·x + c` (including `sqrt(g)`). With `L = lcm` of the involved
//!   denominators, `g = t^L` (i.e. `x = (t^L − c)/d`) rationalizes.
//!
//! Both mechanisms substitute, reintegrate the resulting rational t-form via
//! [`integrate_raw`], and back-substitute. The entry returns `None` when no
//! integrability condition holds, when the t-form is not genuinely rational
//! in `t` (structural check), when the t-integral keeps an `Integral`
//! residue, or when a budget is exceeded — a wrong answer is never returned.
//! The t-forms are rational, so they cannot re-match this module's own
//! patterns (the integer-exponent binomial shape is declined), keeping the
//! mechanism idempotent under re-entry.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    contains_integral, gcd_i64, int_pow, integrate_raw, inv, is_constant, lcm_i64, linear_form,
    node_count, pick_subst_symbol, rat_atom, replace_symbol,
};

/// Largest denominator accepted for `p` in Chebyshev cases 2/3, and for the
/// case-1 common denominator of `m` and `n` (the substitution degree `s`).
const MAX_DENOM: i64 = 6;
/// Cap on the rationalization degree `L` (lcm of fractional-power
/// denominators) in `x = t^L`.
const MAX_RATIONALIZE_L: i64 = 4;
/// Cap on integer exponents applied to non-monomial t-factors: `|p|` in
/// case 1, `|k−1|` in case 2, `|k+1|` in case 3.
const MAX_INT_EXP: i64 = 8;
/// Node budget for the input integrand and for the substituted t-integrand.
const MAX_SUBST_NODES: usize = 200;
/// Cap on the number of fractional-power sites in the rationalization scan.
const MAX_FRAC_SITES: usize = 16;

/// A rational number in lowest terms (`d > 0`).
#[derive(Clone, Copy)]
struct Rat {
    n: i64,
    d: i64,
}

impl Rat {
    fn new(n: i64, d: i64) -> Option<Rat> {
        if d == 0 {
            return None;
        }
        let (n, d) = if d < 0 {
            (n.checked_neg()?, d.checked_neg()?)
        } else {
            (n, d)
        };
        let g = gcd_i64(n, d);
        Some(Rat { n: n / g, d: d / g })
    }

    fn int(n: i64) -> Rat {
        Rat { n, d: 1 }
    }

    fn is_int(self) -> bool {
        self.d == 1
    }

    fn add(self, o: Rat) -> Option<Rat> {
        let n = self
            .n
            .checked_mul(o.d)?
            .checked_add(o.n.checked_mul(self.d)?)?;
        Rat::new(n, self.d.checked_mul(o.d)?)
    }

    fn mul(self, o: Rat) -> Option<Rat> {
        Rat::new(self.n.checked_mul(o.n)?, self.d.checked_mul(o.d)?)
    }

    fn div(self, o: Rat) -> Option<Rat> {
        if o.n == 0 {
            return None;
        }
        Rat::new(self.n.checked_mul(o.d)?, self.d.checked_mul(o.n)?)
    }
}

/// The result of matching `expr` against `C · x^m · (a + b·x^n)^p`.
struct BinomMatch<'a> {
    /// Constant (w.r.t. `var`) leftover factors.
    rest: Vec<Atom<'a>>,
    /// The binomial base `a + b·x^n` as found (used for back-substitution).
    base: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    m: Rat,
    n: Rat,
    p: Rat,
}

/// Integrate `expr` via Chebyshev binomial substitution or fractional-power
/// rationalization. Returns `None` when neither applies (see module docs).
pub(crate) fn integrate_binomial<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) {
        return None;
    }
    // Top-level sums are distributed over terms by the caller chain.
    if matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    try_chebyshev(ctx, expr, var).or_else(|| try_rationalize(ctx, expr, var))
}

// =========================================================================
// Chebyshev binomial differential
// =========================================================================

fn try_chebyshev<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if node_count(expr) > MAX_SUBST_NODES {
        return None;
    }
    let m0 = match_binom(ctx, expr, var)?;
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let x = ctx.var(var.as_str());

    let (factors, t_back): (Vec<Atom<'a>>, Atom<'a>) = if m0.p.is_int() {
        // Case 1: p ∈ ℤ, at least one of m, n non-integer; x = t^s with
        // s = lcm(den m, den n) makes every exponent integral.
        if m0.m.is_int() && m0.n.is_int() {
            return None;
        }
        let s = lcm_i64(m0.m.d, m0.n.d)?;
        if !(2..=MAX_DENOM).contains(&s) {
            return None;
        }
        let pe = m0.p.n;
        if pe == 0 || pe.abs() > MAX_INT_EXP {
            return None;
        }
        let ms = m0.m.mul(Rat::int(s))?;
        let ns = m0.n.mul(Rat::int(s))?;
        // integrand_t = C · s · t^(m·s + s − 1) · (a + b·t^(n·s))^p
        let te = ms.n.checked_add(s)?.checked_sub(1)?;
        let binom_base_t = ctx.add(&[m0.a, ctx.mul(&[m0.b, int_pow(ctx, t, ns.n)])]);
        let mut factors = m0.rest.clone();
        factors.push(ctx.num(s));
        factors.push(int_pow(ctx, t, te));
        factors.push(int_pow(ctx, binom_base_t, pe));
        let back = ctx.pow(x, rat_atom(ctx, 1, s));
        (factors, back)
    } else {
        let s = m0.p.d; // > 1, p = r/s in lowest terms
        if s > MAX_DENOM {
            return None;
        }
        let r = m0.p.n;
        let k2 = m0.m.add(Rat::int(1))?.div(m0.n)?;
        if k2.is_int() {
            // Case 2: (m+1)/n = k ∈ ℤ, t^s = a + b·x^n.
            // integrand_t = C · s/(b·n) · t^(r+s−1) · ((t^s − a)/b)^(k−1)
            let e = k2.n.checked_sub(1)?;
            if e.abs() > MAX_INT_EXP {
                return None;
            }
            let coeff = ctx.mul(&[
                rat_atom(ctx, s.checked_mul(m0.n.d)?, m0.n.n),
                inv(ctx, m0.b),
            ]);
            let te = r.checked_add(s)?.checked_sub(1)?;
            let base_t = ctx.mul(&[
                ctx.add(&[int_pow(ctx, t, s), ctx.mul(&[ctx.num(-1), m0.a])]),
                inv(ctx, m0.b),
            ]);
            let mut factors = m0.rest.clone();
            factors.push(coeff);
            factors.push(int_pow(ctx, t, te));
            factors.push(int_pow(ctx, base_t, e));
            let back = ctx.pow(m0.base, rat_atom(ctx, 1, s));
            (factors, back)
        } else {
            // Case 3: (m+1)/n + p = k ∈ ℤ, t^s = a·x^(−n) + b.
            // integrand_t = C · (−s)/(a·n) · t^(r+s−1) · (a/(t^s − b))^(k+1)
            let k3 = k2.add(m0.p)?;
            if !k3.is_int() {
                return None;
            }
            let e = k3.n.checked_add(1)?;
            if e.abs() > MAX_INT_EXP {
                return None;
            }
            let coeff = ctx.mul(&[
                rat_atom(ctx, s.checked_mul(m0.n.d)?.checked_neg()?, m0.n.n),
                inv(ctx, m0.a),
            ]);
            let te = r.checked_add(s)?.checked_sub(1)?;
            let frac = ctx.mul(&[
                m0.a,
                ctx.pow(
                    ctx.add(&[int_pow(ctx, t, s), ctx.mul(&[ctx.num(-1), m0.b])]),
                    ctx.num(-1),
                ),
            ]);
            let mut factors = m0.rest.clone();
            factors.push(coeff);
            factors.push(int_pow(ctx, t, te));
            factors.push(int_pow(ctx, frac, e));
            let g = ctx.add(&[
                ctx.mul(&[
                    m0.a,
                    ctx.pow(x, rat_atom(ctx, m0.n.n.checked_neg()?, m0.n.d)),
                ]),
                m0.b,
            ]);
            let back = ctx.pow(g, rat_atom(ctx, 1, s));
            (factors, back)
        }
    };
    finish_substitution(ctx, factors, t_sym, var, t_back)
}

/// Match `expr` as `C · x^m · (a + b·x^n)^p`: pull out `x`-power factors
/// (exponent summed into `m`), exactly one binomial-power factor, and
/// constant leftover factors. Any other non-constant factor → `None`.
fn match_binom<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<BinomMatch<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut m = Rat::int(0);
    let mut rest = Vec::new();
    let mut binom: Option<(Atom<'a>, Atom<'a>, Atom<'a>, Rat, Rat)> = None;
    for f in factors {
        match f.node() {
            AtomNode::Var(v) if *v == var => {
                m = m.add(Rat::int(1))?;
            }
            AtomNode::Pow(b, e) => {
                if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                    let r = rat_of(*e)?;
                    m = m.add(Rat::new(r.0, r.1)?)?;
                } else if let Some((p, q)) = rat_of(*e)
                    && let Some((a, bb, n)) = match_base(ctx, *b, var)
                {
                    if binom.is_some() {
                        return None;
                    }
                    binom = Some((*b, a, bb, n, Rat::new(p, q)?));
                } else if is_constant(f, var) {
                    rest.push(f);
                } else {
                    return None;
                }
            }
            _ => {
                if is_constant(f, var) {
                    rest.push(f);
                } else {
                    return None;
                }
            }
        }
    }
    let (base, a, b, n, p) = binom?;
    Some(BinomMatch {
        rest,
        base,
        a,
        b,
        m,
        n,
        p,
    })
}

/// Match `base` as `a + b·x^n`: exactly one non-constant term of the shape
/// `c·x^k` (`k` rational, nonzero), the remaining terms constant with a
/// nonzero sum.
fn match_base<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Rat)> {
    let AtomNode::Add(args) = base.node() else {
        return None;
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut xterm: Option<(Atom<'a>, Rat)> = None;
    for arg in args.iter() {
        if is_constant(*arg, var) {
            consts.push(*arg);
            continue;
        }
        if xterm.is_some() {
            return None;
        }
        xterm = Some(x_power_term(ctx, *arg, var)?);
    }
    if consts.is_empty() {
        return None;
    }
    let a = normalize(ctx, ctx.add(&consts));
    if matches!(a.node(), AtomNode::Num(0)) {
        return None;
    }
    let (b, n) = xterm?;
    if n.n == 0 || matches!(b.node(), AtomNode::Num(0)) {
        return None;
    }
    Some((a, b, n))
}

/// Match `term` as `c·x^k` with `c` constant w.r.t. `var`; returns `(c, k)`.
fn x_power_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Rat)> {
    match term.node() {
        AtomNode::Var(v) if *v == var => Some((ctx.num(1), Rat::int(1))),
        AtomNode::Pow(b, e) if matches!(b.node(), AtomNode::Var(v) if *v == var) => {
            let (p, q) = rat_of(*e)?;
            Some((ctx.num(1), Rat::new(p, q)?))
        }
        AtomNode::Mul(args) => {
            let mut coeff: Vec<Atom<'a>> = Vec::new();
            let mut n: Option<Rat> = None;
            for f in args.iter() {
                match f.node() {
                    AtomNode::Var(v) if *v == var => {
                        if n.is_some() {
                            return None;
                        }
                        n = Some(Rat::int(1));
                    }
                    AtomNode::Pow(b, e) if matches!(b.node(), AtomNode::Var(v) if *v == var) => {
                        if n.is_some() {
                            return None;
                        }
                        let (p, q) = rat_of(*e)?;
                        n = Some(Rat::new(p, q)?);
                    }
                    _ => {
                        if is_constant(*f, var) {
                            coeff.push(*f);
                        } else {
                            return None;
                        }
                    }
                }
            }
            let n = n?;
            let c = if coeff.is_empty() {
                ctx.num(1)
            } else {
                normalize(ctx, ctx.mul(&coeff))
            };
            Some((c, n))
        }
        _ => None,
    }
}

// =========================================================================
// Fractional-power rationalization
// =========================================================================

fn try_rationalize<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if node_count(expr) > MAX_SUBST_NODES {
        return None;
    }
    let mut base: Option<Atom<'a>> = None;
    let mut dens: Vec<i64> = Vec::new();
    let mut sites = 0usize;
    collect_frac_sites(ctx, expr, var, &mut base, &mut dens, &mut sites)?;
    let g = base?;
    // g = d·x + c, non-constant in x.
    let (d, c) = linear_form(ctx, g, var)?;
    if matches!(d.node(), AtomNode::Num(0)) {
        return None;
    }
    let l = dens.iter().try_fold(1i64, |acc, s| lcm_i64(acc, *s))?;
    if !(2..=MAX_RATIONALIZE_L).contains(&l) {
        return None;
    }
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    // x = (t^L − c)/d, dx/dt = L·t^(L−1)/d.
    let x_t = normalize(
        ctx,
        ctx.mul(&[
            ctx.add(&[int_pow(ctx, t, l), ctx.mul(&[ctx.num(-1), c])]),
            inv(ctx, d),
        ]),
    );
    let substituted = rationalize_subst(ctx, expr, var, g, l, t, x_t)?;
    let dx = ctx.mul(&[ctx.num(l), inv(ctx, d), int_pow(ctx, t, l - 1)]);
    let g_back = ctx.pow(g, rat_atom(ctx, 1, l));
    finish_substitution(ctx, vec![substituted, dx], t_sym, var, g_back)
}

/// Collect the fractional-power sites of `expr`: `Pow(g, k/s)` with `s > 1`
/// and `sqrt(g)` nodes whose base is non-constant in `var`. All sites must
/// share a single (normalized) base, recorded in `base`; a second distinct
/// base or too many sites → `None`.
fn collect_frac_sites<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    base: &mut Option<Atom<'a>>,
    dens: &mut Vec<i64>,
    sites: &mut usize,
) -> Option<()> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if let Some((_, s)) = rat_of(*e)
                && s > 1
                && !is_constant(*b, var)
            {
                let nb = normalize(ctx, *b);
                match *base {
                    Some(g0) if nb != g0 => return None,
                    None => *base = Some(nb),
                    _ => {}
                }
                dens.push(s);
                *sites += 1;
                if *sites > MAX_FRAC_SITES {
                    return None;
                }
            }
            collect_frac_sites(ctx, *b, var, base, dens, sites)?;
            collect_frac_sites(ctx, *e, var, base, dens, sites)
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            if !is_constant(args[0], var) {
                let nb = normalize(ctx, args[0]);
                match *base {
                    Some(g0) if nb != g0 => return None,
                    None => *base = Some(nb),
                    _ => {}
                }
                dens.push(2);
                *sites += 1;
                if *sites > MAX_FRAC_SITES {
                    return None;
                }
            }
            collect_frac_sites(ctx, args[0], var, base, dens, sites)
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_frac_sites(ctx, *a, var, base, dens, sites)?;
            }
            Some(())
        }
        AtomNode::Num(_) | AtomNode::Var(_) => Some(()),
    }
}

/// Substitute `x → x_t` and every fractional-power site `g^(k/s) → t^(k·L/s)`
/// (integer exponent, since `s | L`), `sqrt(g) → t^(L/2)`. Returns `None` on
/// arithmetic overflow only; structural inapplicability is caught earlier.
fn rationalize_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    g: Atom<'a>,
    l: i64,
    t: Atom<'a>,
    x_t: Atom<'a>,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if let Some((k, s)) = rat_of(*e)
                && s > 1
                && !is_constant(*b, var)
                && normalize(ctx, *b) == g
            {
                return Some(int_pow(ctx, t, k.checked_mul(l / s)?));
            }
            let nb = rationalize_subst(ctx, *b, var, g, l, t, x_t)?;
            let ne = rationalize_subst(ctx, *e, var, g, l, t, x_t)?;
            Some(ctx.pow(nb, ne))
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            if !is_constant(args[0], var) && normalize(ctx, args[0]) == g {
                return Some(int_pow(ctx, t, l / 2));
            }
            let na = rationalize_subst(ctx, args[0], var, g, l, t, x_t)?;
            Some(ctx.fun("sqrt", &[na]))
        }
        AtomNode::Var(v) => {
            if *v == var {
                Some(x_t)
            } else {
                Some(expr)
            }
        }
        AtomNode::Num(_) => Some(expr),
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(rationalize_subst(ctx, *a, var, g, l, t, x_t)?);
            }
            Some(if matches!(expr.node(), AtomNode::Add(_)) {
                ctx.add(&rebuilt)
            } else {
                ctx.mul(&rebuilt)
            })
        }
        AtomNode::Fun(name, args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(rationalize_subst(ctx, *a, var, g, l, t, x_t)?);
            }
            Some(ctx.fun(name.as_str(), &rebuilt))
        }
    }
}

// =========================================================================
// Shared substitution plumbing
// =========================================================================

/// Reintegrate the substituted t-form and map `t` back to an x-expression.
/// Declines (`None`) when the t-form is not genuinely rational in `t`, when
/// the node budget is exceeded, or when the t-integral has a residue.
fn finish_substitution<'a>(
    ctx: &'a AtomArena<'a>,
    factors: Vec<Atom<'a>>,
    t_sym: Symbol,
    var: Symbol,
    t_back: Atom<'a>,
) -> Option<Atom<'a>> {
    let integrand_t = normalize(ctx, ctx.mul(&factors));
    if node_count(integrand_t) > MAX_SUBST_NODES {
        return None;
    }
    if !is_rational_in(integrand_t, var) {
        return None;
    }
    let result_t = integrate_raw(ctx, integrand_t, t_sym, 0, true, 0, 0);
    if contains_integral(result_t) {
        return None;
    }
    let back = replace_symbol(ctx, result_t, t_sym, t_back);
    Some(normalize(ctx, back))
}

/// Structural check: `expr` is a rational function of the substitution
/// variable — only integer powers, no function nodes, and no leftover
/// occurrence of the original integration variable (which would silently
/// integrate as a constant of the t-integral).
fn is_rational_in(expr: Atom<'_>, orig_var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        // Constant symbols are fine; the original variable must be fully
        // substituted away.
        AtomNode::Var(v) => *v != orig_var,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| is_rational_in(*a, orig_var))
        }
        AtomNode::Pow(b, e) => matches!(e.node(), AtomNode::Num(_)) && is_rational_in(*b, orig_var),
        AtomNode::Fun(_, _) => false,
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification (handles the operators
    /// the substitution results produce).
    fn eval_f64(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
        match expr.node() {
            AtomNode::Num(n) => Some(*n as f64),
            AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, val)| *val),
            AtomNode::Add(args) => args
                .iter()
                .try_fold(0.0, |acc, a| Some(acc + eval_f64(*a, env)?)),
            AtomNode::Mul(args) => args
                .iter()
                .try_fold(1.0, |acc, a| Some(acc * eval_f64(*a, env)?)),
            AtomNode::Pow(b, e) => Some(eval_f64(*b, env)?.powf(eval_f64(*e, env)?)),
            AtomNode::Fun(name, args) => {
                let v = eval_f64(*args.first()?, env)?;
                Some(match name.as_str() {
                    "sin" => v.sin(),
                    "cos" => v.cos(),
                    "tan" => v.tan(),
                    "exp" => v.exp(),
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    _ => return None,
                })
            }
        }
    }

    /// Run the mechanism directly, require `Some` without residue, and check
    /// `diff(result) == integrand` numerically at the given sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_binomial(ctx, integrand, var).expect("mechanism declined");
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(integrand, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
    }

    /// `base^(p/q)` with the codebase's rational-exponent atom shape.
    fn powq<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, p: i64, q: i64) -> Atom<'a> {
        ctx.pow(base, rat_atom(ctx, p, q))
    }

    #[test]
    fn chebyshev_case2_numeric() {
        // ∫ x^3·(1 + x^2)^(1/2) dx: m=3, n=2, p=1/2, (m+1)/n = 2 ∈ ℤ.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let base = ctx.add(&[ctx.num(1), ctx.pow(x, ctx.num(2))]);
        let expr = ctx.mul(&[ctx.pow(x, ctx.num(3)), powq(&ctx, base, 1, 2)]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn chebyshev_case2_symbolic_consts() {
        // ∫ x^2·(a + b·x^3)^(1/2) dx: m=2, n=3, p=1/2, (m+1)/n = 1 ∈ ℤ;
        // answer 2(a + b·x^3)^(3/2)/(9b).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let base = ctx.add(&[a, ctx.mul(&[b, ctx.pow(x, ctx.num(3))])]);
        let expr = ctx.mul(&[ctx.pow(x, ctx.num(2)), powq(&ctx, base, 1, 2)]);
        let env = [(Symbol::new("a"), 0.7), (Symbol::new("b"), 1.3)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.4, 0.9, 1.4]);
    }

    #[test]
    fn chebyshev_case3_negative_p() {
        // ∫ (1 + x^2)^(1/2)/x^2 dx: m=−2, n=2, p=−1/2;
        // (m+1)/n = −1/2 ∉ ℤ, (m+1)/n + p = −1 ∈ ℤ → case 3.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let base = ctx.add(&[ctx.num(1), ctx.pow(x, ctx.num(2))]);
        let expr = ctx.mul(&[powq(&ctx, base, 1, 2), ctx.pow(x, ctx.num(-2))]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.6, 1.1, 2.0]);
    }

    #[test]
    fn chebyshev_case1_fractional_mn() {
        // ∫ x^(1/2)·(1 + x^(1/3))^2 dx: p=2 ∈ ℤ, m=1/2, n=1/3 → x = t^6.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let base = ctx.add(&[ctx.num(1), powq(&ctx, x, 1, 3)]);
        let expr = ctx.mul(&[powq(&ctx, x, 1, 2), ctx.pow(base, ctx.num(2))]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.4, 1.0, 1.9]);
    }

    #[test]
    fn rationalize_corpus_symbolic() {
        // Corpus failure: ∫ x^(3/2)·(A + B·x)/(a + b·x)^3 dx → x = t^2 gives
        // 2t^4·(A + B·t^2)/(a + b·t^2)^3, rational over ℚ(A,B,a,b).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let big_a = ctx.var("A");
        let big_b = ctx.var("B");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let num = ctx.add(&[big_a, ctx.mul(&[big_b, x])]);
        let den = ctx.pow(ctx.add(&[a, ctx.mul(&[b, x])]), ctx.num(-3));
        let expr = ctx.mul(&[powq(&ctx, x, 3, 2), num, den]);
        let env = [
            (Symbol::new("A"), 1.2),
            (Symbol::new("B"), 0.8),
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.5, 1.1, 2.3]);
    }

    #[test]
    fn rationalize_fourth_root_plus_sqrt() {
        // Corpus failure: ∫ dx/((1+x)^(1/4) + sqrt(1+x)) → 1+x = t^4 gives
        // 4t^3/(t + t^2) = 4t − 4 + 4/(t + 1).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let one_px = ctx.add(&[ctx.num(1), x]);
        let fourth = powq(&ctx, one_px, 1, 4);
        let sq = ctx.fun("sqrt", &[one_px]);
        let expr = ctx.pow(ctx.add(&[fourth, sq]), ctx.num(-1));
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn decline_nonintegrable_chebyshev() {
        // Corpus failure (honest decline): (a + b·x^3)^(3/2)/x^5 has
        // m=−5, n=3, p=3/2: (m+1)/n = −4/3 ∉ ℤ, (m+1)/n + p = 1/6 ∉ ℤ.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let base = ctx.add(&[a, ctx.mul(&[b, ctx.pow(x, ctx.num(3))])]);
        let expr = ctx.mul(&[powq(&ctx, base, 3, 2), ctx.pow(x, ctx.num(-5))]);
        assert!(integrate_binomial(&ctx, expr, Symbol::new("x")).is_none());

        // ∫ (1 + x^3)^(1/2) dx: m=0, n=3, p=1/2 — no case applies.
        let expr2 = powq(&ctx, ctx.add(&[ctx.num(1), ctx.pow(x, ctx.num(3))]), 1, 2);
        assert!(integrate_binomial(&ctx, expr2, Symbol::new("x")).is_none());
    }

    #[test]
    fn decline_nonlinear_radical_base() {
        // Corpus failure (honest decline): (b·x + c·x^2)^(2/3) factors as
        // x^(2/3)·(b + c·x)^(2/3) with m=2/3, n=1, p=2/3 — neither
        // Chebyshev condition holds, and the base is not a single x-power.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let base = ctx.add(&[ctx.mul(&[b, x]), ctx.mul(&[c, ctx.pow(x, ctx.num(2))])]);
        let expr = powq(&ctx, base, 2, 3);
        assert!(integrate_binomial(&ctx, expr, Symbol::new("x")).is_none());
    }

    #[test]
    fn decline_subst_symbol_collision() {
        // Both "t" and "u" taken by the integrand → no substitution variable.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let t = ctx.var("t");
        let u = ctx.var("u");
        let base = ctx.add(&[t, ctx.mul(&[u, ctx.pow(x, ctx.num(2))])]);
        let expr = ctx.mul(&[ctx.pow(x, ctx.num(3)), powq(&ctx, base, 1, 2)]);
        assert!(integrate_binomial(&ctx, expr, Symbol::new("x")).is_none());
    }
}
