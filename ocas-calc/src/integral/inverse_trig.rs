//! Inverse-trig/hyperbolic integration mechanisms (0.27.1 Phase 1E).
//!
//! Three mechanisms for integrands built from the six inverse functions
//! `asin/acos/atan/asinh/acosh/atanh`:
//!
//! - **M1 kernel-derivative power rule**: `(a + b·F(u))^k · C·R^e` with
//!   `u = σ·x` (zero intercept) and `R^e` a constant-normalized copy of the
//!   algebraic kernel of `F'(u)`. Per family, with `R = s·E(u)` and the
//!   sign `s_F` from `F'(u) = s_F·E(u)^e`:
//!
//!   | F     | E(u)   | e     | s_F |
//!   |-------|--------|-------|-----|
//!   | asin  | 1−u²   | −1/2  | +1  |
//!   | acos  | 1−u²   | −1/2  | −1  |
//!   | atan  | 1+u²   | −1    | +1  |
//!   | atanh | 1−u²   | −1    | +1  |
//!   | asinh | 1+u²   | −1/2  | +1  |
//!   | acosh | u²−1   | −1/2  | +1  |
//!
//!   `∫ = C·s^e·s_F·(a + b·F(u))^(k+1) / (b·σ·(k+1))` for integer
//!   `1 ≤ k ≤ 6` (`k = 0` is the bare form — M3/rule-family-E territory).
//! - **M2 inv-hyp substitution**: a power `(a + b·F(arg))^k`
//!   (`F ∈ {asinh, acosh}`, `arg` linear with nonzero slope,
//!   `1 ≤ k ≤ 4`) times polynomial factors (total degree ≤ 4) and at most
//!   one radical `(s·E)^(m/2)` (`m ∈ {1, 3}`, zero-intercept `arg` only)
//!   maps under `t = F(arg)` to a `poly(t)·cosh(t)^A·sinh(t)^B` form
//!   (`x = (G(t)−ρ)/σ`, `G = sinh`/`cosh`, `dx = G'(t)/σ·dt`,
//!   `(s·E)^(m/2) = s^(m/2)·G'(t)^m`). The hyperbolic powers are reduced
//!   to multiple angles in-module (`cosh^A·sinh^B → Σ c_m·H(m·t)` via the
//!   exponential binomial expansion), the sum is re-integrated through the
//!   chain (`integrate_raw`) and declined on any `Integral` residue. When
//!   the chain leaves a residue on `t^j·H(m·t)` terms with `j ≥ 3` (its
//!   parts budget caps at depth 2), a direct tabular closed form — the
//!   hyperbolic analogue of trig_reduction's T3 — finishes the same terms.
//! - **M3 linear-argument bare forms**: `C·F(a + b·x)` integrates by
//!   `u = a + b·x` to the textbook closed forms (the rule-family-E shapes
//!   divided by the slope).
//!
//! Sign convention: a *numeric* radicand normalization `s` under a
//! square-root power must be positive; a *symbolic* `s` is assumed
//! positive (the same convention rule family G uses for `asin(x/a)`),
//! which is what lets corpus radicands like `(d + c²·d·x²)` through.
//!
//! All chain re-entries go through `integrate_raw` and are declined on any
//! `Integral` residue. Outputs keep the inverse function inside
//! `sinh`/`cosh`/power arguments, so they never re-match this module's own
//! top-level factor patterns (idempotent under chain re-entry).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    contains_integral, int_pow, integrate_raw, inv, is_constant, linear_form, node_count,
    pick_subst_symbol, rat_atom, replace_symbol,
};

/// Node budget for the input integrand.
const MAX_NODES: usize = 200;
/// Node budget for the substituted t-form before multiple-angle reduction.
const MAX_SUBST_NODES: usize = 300;
/// Node budget for the reduced t-form that is re-integrated.
const MAX_REDUCED_NODES: usize = 600;
/// M1 cap on the `(a + b·F(u))^k` exponent.
const MAX_POWER_K: i64 = 6;
/// M2 cap on the `(a + b·F(arg))^k` exponent.
const MAX_SUBST_K: i64 = 4;
/// M2 cap on the total degree of polynomial factors.
const MAX_POLY_DEG: usize = 4;
/// Cap on `A + B` in the `cosh(t)^A·sinh(t)^B` reduction.
const MAX_HYP_EXP: usize = 8;

/// The six inverse functions covered by this module.
#[derive(Clone, Copy, PartialEq, Eq)]
enum InvFun {
    Asin,
    Acos,
    Atan,
    Asinh,
    Acosh,
    Atanh,
}

impl InvFun {
    fn name(self) -> &'static str {
        match self {
            InvFun::Asin => "asin",
            InvFun::Acos => "acos",
            InvFun::Atan => "atan",
            InvFun::Asinh => "asinh",
            InvFun::Acosh => "acosh",
            InvFun::Atanh => "atanh",
        }
    }
}

fn inv_fun(name: &str) -> Option<InvFun> {
    Some(match name {
        "asin" => InvFun::Asin,
        "acos" => InvFun::Acos,
        "atan" => InvFun::Atan,
        "asinh" => InvFun::Asinh,
        "acosh" => InvFun::Acosh,
        "atanh" => InvFun::Atanh,
        _ => return None,
    })
}

/// Kernel exponent shape: `R^(−1/2)` (also `1/sqrt(R)`) or `R^(−1)`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum KernelExp {
    InvSqrt,
    Recip,
}

/// The kernel exponent shape expected per family.
fn family_kernel_exp(f: InvFun) -> KernelExp {
    match f {
        InvFun::Asin | InvFun::Acos | InvFun::Asinh | InvFun::Acosh => KernelExp::InvSqrt,
        InvFun::Atan | InvFun::Atanh => KernelExp::Recip,
    }
}

/// `s_F` from `F'(u) = s_F·E(u)^e`: only `acos` differentiates with a
/// minus sign.
fn family_sign(f: InvFun) -> i64 {
    match f {
        InvFun::Acos => -1,
        _ => 1,
    }
}

/// The kernel sum `E` per family in terms of the slope `σ` of `u = σ·x`:
/// `1−σ²x²` (asin/acos/atanh), `1+σ²x²` (atan/asinh), `σ²x²−1` (acosh) —
/// built as `σ²·x²` so it matches the parser's `c^2*x^2` shape — together
/// with the constant term `e0 ∈ {1, −1}` (used to read the normalization
/// `s` off `R(0) = s·e0`).
fn kernel_sum<'a>(
    ctx: &'a AtomArena<'a>,
    f: InvFun,
    sigma: Atom<'a>,
    var: Symbol,
) -> (Atom<'a>, i64) {
    let x = ctx.var(var.as_str());
    let s2x2 = ctx.mul(&[int_pow(ctx, sigma, 2), int_pow(ctx, x, 2)]);
    match f {
        InvFun::Asin | InvFun::Acos | InvFun::Atanh => (
            normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), s2x2])])),
            1,
        ),
        InvFun::Atan | InvFun::Asinh => (normalize(ctx, ctx.add(&[ctx.num(1), s2x2])), 1),
        InvFun::Acosh => (normalize(ctx, ctx.add(&[ctx.num(-1), s2x2])), -1),
    }
}

/// True when `e` folds to numeric zero (structurally after `normalize`,
/// then via `collect_terms` for mixed symbolic cancellations).
pub(crate) fn folds_to_zero<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> bool {
    let n = normalize(ctx, e);
    if matches!(n.node(), AtomNode::Num(0)) {
        return true;
    }
    matches!(
        crate::ode::util::collect_terms(ctx, n).node(),
        AtomNode::Num(0)
    )
}

fn is_zero_atom<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> bool {
    folds_to_zero(ctx, e)
}

/// A matched `(a + b·F(arg))^k` factor; `base` is the original
/// `(a + b·F(arg))` atom so the answer reuses the input's exact shape.
struct InvPow<'a> {
    f: InvFun,
    arg: Atom<'a>,
    #[allow(dead_code)] // kept for M2's t-form rebuild readability
    a: Atom<'a>,
    b: Atom<'a>,
    k: i64,
    base: Atom<'a>,
}

/// Match `(a + b·F(arg))`: a bare `F(arg)`, a constant multiple
/// `b·F(arg)`, or a sum of constants plus exactly one such term.
fn match_inv_base<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    var: Symbol,
) -> Option<(InvFun, Atom<'a>, Atom<'a>, Atom<'a>)> {
    match base.node() {
        AtomNode::Fun(name, args) if args.len() == 1 => {
            let f = inv_fun(name.as_str())?;
            Some((f, args[0], ctx.num(0), ctx.num(1)))
        }
        AtomNode::Mul(args) => {
            let mut coeff: Vec<Atom<'a>> = Vec::new();
            let mut found: Option<(InvFun, Atom<'a>)> = None;
            for fa in args.iter() {
                match fa.node() {
                    AtomNode::Fun(name, fargs) if fargs.len() == 1 => {
                        let f = inv_fun(name.as_str())?;
                        if found.is_some() {
                            return None;
                        }
                        found = Some((f, fargs[0]));
                    }
                    _ if is_constant(*fa, var) => coeff.push(*fa),
                    _ => return None,
                }
            }
            let (f, arg) = found?;
            let b = if coeff.is_empty() {
                ctx.num(1)
            } else {
                normalize(ctx, ctx.mul(&coeff))
            };
            Some((f, arg, ctx.num(0), b))
        }
        AtomNode::Add(args) => {
            let mut consts: Vec<Atom<'a>> = Vec::new();
            let mut found: Option<(InvFun, Atom<'a>, Atom<'a>, Atom<'a>)> = None;
            for t in args.iter() {
                if is_constant(*t, var) {
                    consts.push(*t);
                    continue;
                }
                let term = match_inv_base(ctx, *t, var)?;
                if found.is_some() {
                    return None;
                }
                found = Some(term);
            }
            let (f, arg, _zero, b) = found?;
            let a = if consts.is_empty() {
                ctx.num(0)
            } else {
                normalize(ctx, ctx.add(&consts))
            };
            Some((f, arg, a, b))
        }
        _ => None,
    }
}

/// Match `(a + b·F(arg))^k` with integer `1 ≤ k ≤ k_cap`, or the bare base
/// as `k = 1`.
fn match_inv_power<'a>(
    ctx: &'a AtomArena<'a>,
    factor: Atom<'a>,
    var: Symbol,
    k_cap: i64,
) -> Option<InvPow<'a>> {
    if let AtomNode::Pow(b, e) = factor.node()
        && let AtomNode::Num(k) = e.node()
    {
        if !(1..=k_cap).contains(k) {
            return None;
        }
        let (f, arg, a, bb) = match_inv_base(ctx, *b, var)?;
        return Some(InvPow {
            f,
            arg,
            a,
            b: bb,
            k: *k,
            base: *b,
        });
    }
    let (f, arg, a, bb) = match_inv_base(ctx, factor, var)?;
    Some(InvPow {
        f,
        arg,
        a,
        b: bb,
        k: 1,
        base: factor,
    })
}

/// Match the M1 kernel factor: `R^(−1/2)`, `1/sqrt(R)`, or `R^(−1)`.
fn match_kernel(factor: Atom<'_>) -> Option<(Atom<'_>, KernelExp)> {
    let AtomNode::Pow(b, e) = factor.node() else {
        return None;
    };
    if let AtomNode::Fun(name, args) = b.node()
        && name.as_str() == "sqrt"
        && args.len() == 1
        && matches!(e.node(), AtomNode::Num(-1))
    {
        return Some((args[0], KernelExp::InvSqrt));
    }
    match rat_of(*e)? {
        (-1, 2) => Some((*b, KernelExp::InvSqrt)),
        (-1, 1) => Some((*b, KernelExp::Recip)),
        _ => None,
    }
}

/// Read the radicand normalization `s` from `R = s·E` via `R(0) = s·e0`,
/// then verify the full identity `R ≡ s·E` symbolically. Returns `None`
/// when `s` is zero or the identity does not hold.
fn radicand_normalization<'a>(
    ctx: &'a AtomArena<'a>,
    r: Atom<'a>,
    e_sum: Atom<'a>,
    e0: i64,
    var: Symbol,
) -> Option<Atom<'a>> {
    let r0 = normalize(ctx, replace_symbol(ctx, r, var, ctx.num(0)));
    let s_atom = if e0 == 1 {
        r0
    } else {
        normalize(ctx, ctx.mul(&[ctx.num(-1), r0]))
    };
    if is_zero_atom(ctx, s_atom) {
        return None;
    }
    let diff = ctx.add(&[r, ctx.mul(&[ctx.num(-1), s_atom, e_sum])]);
    if !folds_to_zero(ctx, diff) {
        return None;
    }
    Some(s_atom)
}

/// Gate on the square-root normalization: a numeric `s` must be positive
/// (real square root); a symbolic `s` is assumed positive (module docs).
fn sqrt_norm_positive(s_atom: Atom<'_>) -> bool {
    match rat_of(s_atom) {
        Some((p, q)) => (p > 0) == (q > 0),
        None => true,
    }
}

// =========================================================================
// M1: kernel-derivative power rule
// =========================================================================

/// `∫ (a + b·F(u))^k · C·R^e dx` with `R^e` the kernel of `F'(u)` — see the
/// module table for the per-family constants.
pub(crate) fn integrate_kernel_power<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut inv_f: Option<InvPow<'a>> = None;
    let mut kern: Option<(Atom<'a>, KernelExp)> = None;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if inv_f.is_none()
            && let Some(m) = match_inv_power(ctx, f, var, MAX_POWER_K)
        {
            inv_f = Some(m);
            continue;
        }
        if kern.is_none()
            && let Some(kd) = match_kernel(f)
        {
            kern = Some(kd);
            continue;
        }
        return None;
    }
    let m = inv_f?;
    let (r, ke) = kern?;
    if ke != family_kernel_exp(m.f) {
        return None;
    }
    let (sigma, intercept) = linear_form(ctx, m.arg, var)?;
    if !is_zero_atom(ctx, intercept) || is_zero_atom(ctx, sigma) {
        return None;
    }
    let (e_sum, e0) = kernel_sum(ctx, m.f, sigma, var);
    let s_atom = radicand_normalization(ctx, r, e_sum, e0, var)?;
    if ke == KernelExp::InvSqrt && !sqrt_norm_positive(s_atom) {
        return None;
    }
    let s_factor = match ke {
        KernelExp::Recip => inv(ctx, s_atom),
        KernelExp::InvSqrt => ctx.pow(s_atom, rat_atom(ctx, -1, 2)),
    };
    let mut out = rest;
    if family_sign(m.f) < 0 {
        out.push(ctx.num(-1));
    }
    out.push(s_factor);
    out.push(int_pow(ctx, m.base, m.k + 1));
    let denom = normalize(ctx, ctx.mul(&[m.b, sigma, ctx.num(m.k + 1)]));
    out.push(inv(ctx, denom));
    Some(normalize(ctx, ctx.mul(&out)))
}

// =========================================================================
// M2: inv-hyp substitution with multiple-angle reduction
// =========================================================================

/// Degree of `expr` as a polynomial in `var` (non-negative integer powers
/// only); `None` when any non-polynomial part occurs.
fn poly_deg(expr: Atom<'_>, var: Symbol) -> Option<usize> {
    match expr.node() {
        AtomNode::Num(_) => Some(0),
        AtomNode::Var(v) => Some(if *v == var { 1 } else { 0 }),
        AtomNode::Add(args) => args
            .iter()
            .try_fold(0usize, |d, a| Some(d.max(poly_deg(*a, var)?))),
        AtomNode::Mul(args) => args
            .iter()
            .try_fold(0usize, |d, a| Some(d + poly_deg(*a, var)?)),
        AtomNode::Pow(b, e) => {
            let AtomNode::Num(n) = e.node() else {
                return None;
            };
            if *n < 0 {
                return None;
            }
            Some(poly_deg(*b, var)?.checked_mul(*n as usize)?)
        }
        AtomNode::Fun(_, _) => None,
    }
}

/// Match `R^(m/2)` with odd `m ∈ {1, 3}` (also the `sqrt(R)` / `sqrt(R)^m`
/// spellings); negative exponents produce negative hyperbolic powers the
/// reduction cannot handle, so they are out of scope.
fn match_sqrt_pow(factor: Atom<'_>) -> Option<(Atom<'_>, i64)> {
    match factor.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            Some((args[0], 1))
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && let AtomNode::Num(n) = e.node()
            {
                return match n {
                    1 | 3 => Some((args[0], *n)),
                    _ => None,
                };
            }
            let (p, q) = rat_of(*e)?;
            if q == 2 && matches!(p, 1 | 3) {
                Some((*b, p))
            } else {
                None
            }
        }
        _ => None,
    }
}

/// `∫ C·(a + b·F(arg))^k · poly(x) · (s·E)^(m/2) dx` for
/// `F ∈ {asinh, acosh}` via `t = F(arg)` — see the module docs.
pub(crate) fn integrate_invhyp_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut inv_f: Option<InvPow<'a>> = None;
    let mut polys: Vec<Atom<'a>> = Vec::new();
    let mut sqrt_f: Option<(Atom<'a>, i64)> = None;
    let mut total_deg = 0usize;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if inv_f.is_none()
            && let Some(m) = match_inv_power(ctx, f, var, MAX_SUBST_K)
        {
            if !matches!(m.f, InvFun::Asinh | InvFun::Acosh) {
                return None;
            }
            inv_f = Some(m);
            continue;
        }
        if sqrt_f.is_none()
            && let Some(sp) = match_sqrt_pow(f)
        {
            sqrt_f = Some(sp);
            continue;
        }
        if let Some(d) = poly_deg(f, var)
            && d >= 1
            && total_deg + d <= MAX_POLY_DEG
        {
            total_deg += d;
            polys.push(f);
            continue;
        }
        return None;
    }
    let m = inv_f?;
    let (sigma, rho) = linear_form(ctx, m.arg, var)?;
    if is_zero_atom(ctx, sigma) {
        return None;
    }
    if sqrt_f.is_some() && !is_zero_atom(ctx, rho) {
        return None;
    }
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let (g_name, gp_name) = match m.f {
        InvFun::Asinh => ("sinh", "cosh"),
        InvFun::Acosh => ("cosh", "sinh"),
        _ => return None,
    };
    let g_t = ctx.fun(g_name, &[t]);
    let gp_t = ctx.fun(gp_name, &[t]);
    // x = (G(t) − ρ)/σ
    let x_of_t = {
        let num_t = if is_zero_atom(ctx, rho) {
            g_t
        } else {
            ctx.add(&[g_t, ctx.mul(&[ctx.num(-1), rho])])
        };
        ctx.mul(&[num_t, inv(ctx, sigma)])
    };
    let mut t_factors = rest;
    let base_t = normalize(ctx, ctx.add(&[m.a, ctx.mul(&[m.b, t])]));
    t_factors.push(int_pow(ctx, base_t, m.k));
    for p in &polys {
        t_factors.push(replace_symbol(ctx, *p, var, x_of_t));
    }
    // dx = G'(t)/σ·dt
    t_factors.push(ctx.mul(&[gp_t, inv(ctx, sigma)]));
    if let Some((r, mpow)) = sqrt_f {
        let (e_sum, e0) = kernel_sum(ctx, m.f, sigma, var);
        let s_atom = radicand_normalization(ctx, r, e_sum, e0, var)?;
        if !sqrt_norm_positive(s_atom) {
            return None;
        }
        // (s·E)^(m/2) = s^(m/2)·G'(t)^m, since E = G'(t)² exactly.
        t_factors.push(ctx.pow(s_atom, rat_atom(ctx, mpow, 2)));
        t_factors.push(int_pow(ctx, gp_t, mpow));
    }
    let product = normalize(ctx, ctx.mul(&t_factors));
    let expanded = crate::expand::expand_bounded(ctx, product).unwrap_or(product);
    let folded = crate::ode::util::collect_terms(ctx, normalize(ctx, expanded));
    if node_count(folded) > MAX_SUBST_NODES {
        return None;
    }
    let terms: Vec<Atom<'a>> = match folded.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![folded],
    };
    let mut reduced: Vec<Atom<'a>> = Vec::new();
    for term in terms {
        let (coeff, j, na, nb) = scan_hyp_term(ctx, term, t_sym)?;
        reduce_hyp_term(ctx, &coeff, j, na, nb, t, &mut reduced)?;
    }
    let t_form = normalize(ctx, ctx.add(&reduced));
    if node_count(t_form) > MAX_REDUCED_NODES {
        return None;
    }
    let back = ctx.fun(m.f.name(), &[m.arg]);
    // Spec path: re-integrate the hyperbolic×polynomial form through the
    // chain (hyperbolic power rules D + parts), decline on residue.
    let chain = integrate_raw(ctx, t_form, t_sym, 0, true, 0, 0);
    if !contains_integral(chain) {
        return Some(normalize(ctx, replace_symbol(ctx, chain, t_sym, back)));
    }
    // Fallback: the chain's parts budget (depth 2) cannot finish
    // `t^j·H(m·t)` for `j ≥ 3`; integrate the reduced terms by the
    // tabular closed form instead (T3 analogue).
    let rterms: Vec<Atom<'a>> = match t_form.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![t_form],
    };
    let mut out = Vec::with_capacity(rterms.len());
    for term in rterms {
        out.push(integrate_hyp_term(ctx, term, t_sym)?);
    }
    let res_t = normalize(ctx, ctx.add(&out));
    Some(normalize(ctx, replace_symbol(ctx, res_t, t_sym, back)))
}

/// `C(n, k)` for small non-negative integers (0 outside the range).
fn binom(n: i64, mut k: i64) -> i64 {
    if k < 0 || k > n {
        return 0;
    }
    if k > n - k {
        k = n - k;
    }
    let mut r: i64 = 1;
    for i in 0..k {
        r = r.saturating_mul(n - i) / (i + 1);
    }
    r
}

/// True when `e` is exactly the substitution variable.
fn is_t_var(e: Atom<'_>, t_sym: Symbol) -> bool {
    matches!(e.node(), AtomNode::Var(v) if *v == t_sym)
}

/// Decompose an expanded t-term into `(constants, j, A, B)` for
/// `C·t^j·cosh(t)^A·sinh(t)^B`; `None` on any other shape.
fn scan_hyp_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    t_sym: Symbol,
) -> Option<(Vec<Atom<'a>>, i64, usize, usize)> {
    let _ = ctx;
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut j = 0i64;
    let mut na = 0usize;
    let mut nb = 0usize;
    for f in factors {
        if is_constant(f, t_sym) {
            coeff.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Var(v) if *v == t_sym => j += 1,
            AtomNode::Pow(b, e) if matches!(b.node(), AtomNode::Var(v) if *v == t_sym) => {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                if *n < 0 {
                    return None;
                }
                j += n;
            }
            AtomNode::Fun(name, args) if args.len() == 1 && is_t_var(args[0], t_sym) => {
                match name.as_str() {
                    "cosh" => na += 1,
                    "sinh" => nb += 1,
                    _ => return None,
                }
            }
            AtomNode::Pow(b, e) => {
                let AtomNode::Fun(name, args) = b.node() else {
                    return None;
                };
                if args.len() != 1 || !is_t_var(args[0], t_sym) {
                    return None;
                }
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                if *n <= 0 {
                    return None;
                }
                match name.as_str() {
                    "cosh" => na += *n as usize,
                    "sinh" => nb += *n as usize,
                    _ => return None,
                }
            }
            _ => return None,
        }
    }
    if j > MAX_SUBST_K || na + nb > MAX_HYP_EXP {
        return None;
    }
    Some((coeff, j, na, nb))
}

/// Expand `cosh(t)^A·sinh(t)^B` into a multiple-angle sum and push the
/// resulting `coeff·t^j·H(m·t)` terms. Via `(e^t+e^{−t})^A·(e^t−e^{−t})^B`:
/// the `e^{m·t}` coefficient pairs with the `e^{−m·t}` one (sign
/// `(−1)^B`), giving `2·c·cosh(m·t)` (B even) or `2·c·sinh(m·t)` (B odd)
/// for `m > 0`, and the plain constant `c(0)` for `m = 0` (which computes
/// to zero for odd B, matching `sinh^odd`'s zero constant term).
fn reduce_hyp_term<'a>(
    ctx: &'a AtomArena<'a>,
    coeff: &[Atom<'a>],
    j: i64,
    na: usize,
    nb: usize,
    t: Atom<'a>,
    out: &mut Vec<Atom<'a>>,
) -> Option<()> {
    if na == 0 && nb == 0 {
        let mut fs = coeff.to_vec();
        if j > 0 {
            fs.push(int_pow(ctx, t, j));
        }
        out.push(normalize(ctx, ctx.mul(&fs)));
        return Some(());
    }
    let nsum = na + nb;
    if nsum > MAX_HYP_EXP {
        return None;
    }
    let mut m = nsum % 2;
    while m <= nsum {
        // Coefficient of e^{m·t}: m = 2(i + j') − (A + B).
        let target = ((m + nsum) / 2) as i64;
        let mut sum: i64 = 0;
        for i in 0..=(na as i64) {
            let jj = target - i;
            if jj < 0 || jj > nb as i64 {
                continue;
            }
            let c = binom(na as i64, i) * binom(nb as i64, jj);
            sum += if (nb as i64 - jj) % 2 == 0 { c } else { -c };
        }
        if sum != 0 {
            let num = if m == 0 { sum } else { 2 * sum };
            let r = rat_atom(ctx, num, 1i64 << nsum);
            let mut fs = coeff.to_vec();
            fs.push(r);
            if j > 0 {
                fs.push(int_pow(ctx, t, j));
            }
            if m > 0 {
                let marg = if m == 1 {
                    t
                } else {
                    ctx.mul(&[ctx.num(m as i64), t])
                };
                fs.push(ctx.fun(if nb.is_multiple_of(2) { "cosh" } else { "sinh" }, &[marg]));
            }
            out.push(normalize(ctx, ctx.mul(&fs)));
        }
        m += 2;
    }
    Some(())
}

/// The multiple `m` of a reduced term's hyperbolic argument: `t` → 1,
/// `m·t` with numeric `m > 0` → `m`.
fn hyp_multiple(arg: Atom<'_>, t_sym: Symbol) -> Option<i64> {
    match arg.node() {
        AtomNode::Var(v) if *v == t_sym => Some(1),
        AtomNode::Mul(args) if args.len() == 2 => {
            let (a, b) = (args[0], args[1]);
            match (a.node(), b.node()) {
                (AtomNode::Num(m), AtomNode::Var(v)) if *v == t_sym && *m > 0 => Some(*m),
                (AtomNode::Var(v), AtomNode::Num(m)) if *v == t_sym && *m > 0 => Some(*m),
                _ => None,
            }
        }
        _ => None,
    }
}

/// Closed form of one reduced term `C·t^j·H(m·t)` (`H ∈ {cosh, sinh}`,
/// `m ≥ 1`), or `C·t^j` when no hyperbolic factor is present: tabular
/// repeated parts, `∫ t^j·cosh(m·t) dt = Σ_{i=0..j} (−1)^i·j^ underline_i·
/// t^(j−i)·H_i(m·t)/m^(i+1)` with `H_i` alternating `sinh, cosh, …`
/// (sinh first when integrating cosh, and conversely). No chain re-entry.
fn integrate_hyp_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    t_sym: Symbol,
) -> Option<Atom<'a>> {
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut j = 0i64;
    let mut hyp: Option<(bool, i64)> = None;
    for f in factors {
        if is_constant(f, t_sym) {
            coeff.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Var(v) if *v == t_sym => j += 1,
            AtomNode::Pow(b, e) if matches!(b.node(), AtomNode::Var(v) if *v == t_sym) => {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                if *n < 0 {
                    return None;
                }
                j += n;
            }
            AtomNode::Fun(name, args) if args.len() == 1 => {
                let is_sinh = match name.as_str() {
                    "sinh" => true,
                    "cosh" => false,
                    _ => return None,
                };
                if hyp.is_some() {
                    return None;
                }
                hyp = Some((is_sinh, hyp_multiple(args[0], t_sym)?));
            }
            _ => return None,
        }
    }
    if j > MAX_SUBST_K {
        return None;
    }
    let t = ctx.var(t_sym.as_str());
    let (is_sinh, m) = match hyp {
        None => {
            // C·t^j → C·t^(j+1)/(j+1)
            let mut fs = coeff;
            fs.push(int_pow(ctx, t, j + 1));
            fs.push(inv(ctx, ctx.num(j + 1)));
            return Some(normalize(ctx, ctx.mul(&fs)));
        }
        Some(h) => h,
    };
    let mut terms = Vec::new();
    let mut fall: i64 = 1; // falling factorial j!/(j−i)!
    for i in 0..=j {
        if i > 0 {
            fall = fall.checked_mul(j - i + 1)?;
        }
        let mpow = m.checked_pow((i + 1) as u32)?;
        let out_sinh = if is_sinh { i % 2 == 1 } else { i % 2 == 0 };
        let marg = if m == 1 { t } else { ctx.mul(&[ctx.num(m), t]) };
        let mut fs = coeff.clone();
        if i % 2 == 1 {
            fs.push(ctx.num(-1));
        }
        if fall != 1 {
            fs.push(ctx.num(fall));
        }
        let deg = j - i;
        if deg > 0 {
            fs.push(int_pow(ctx, t, deg));
        }
        fs.push(ctx.pow(ctx.num(mpow), ctx.num(-1)));
        fs.push(ctx.fun(if out_sinh { "sinh" } else { "cosh" }, &[marg]));
        terms.push(normalize(ctx, ctx.mul(&fs)));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

// =========================================================================
// M3: linear-argument bare forms
// =========================================================================

/// `∫ C·F(a + b·x) dx` for the six inverse functions, via `u = a + b·x`
/// and the textbook bare forms (rule family E) divided by the slope.
pub(crate) fn integrate_bare_linear<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut bare: Option<(InvFun, Atom<'a>)> = None;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if bare.is_none()
            && let AtomNode::Fun(name, args) = f.node()
            && args.len() == 1
            && inv_fun(name.as_str()).is_some()
        {
            bare = Some((inv_fun(name.as_str())?, args[0]));
            continue;
        }
        return None;
    }
    let (f, u) = bare?;
    let (slope, _intercept) = linear_form(ctx, u, var)?;
    if is_zero_atom(ctx, slope) {
        return None;
    }
    let u2 = int_pow(ctx, u, 2);
    let one_minus_u2 = || ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), u2])]);
    let fu = ctx.fun(f.name(), &[u]);
    let core = match f {
        InvFun::Asin => ctx.add(&[ctx.mul(&[u, fu]), ctx.fun("sqrt", &[one_minus_u2()])]),
        InvFun::Acos => ctx.add(&[
            ctx.mul(&[u, fu]),
            ctx.mul(&[ctx.num(-1), ctx.fun("sqrt", &[one_minus_u2()])]),
        ]),
        InvFun::Atan => ctx.add(&[
            ctx.mul(&[u, fu]),
            ctx.mul(&[
                rat_atom(ctx, -1, 2),
                ctx.fun("log", &[ctx.add(&[ctx.num(1), u2])]),
            ]),
        ]),
        InvFun::Asinh => ctx.add(&[
            ctx.mul(&[u, fu]),
            ctx.mul(&[ctx.num(-1), ctx.fun("sqrt", &[ctx.add(&[ctx.num(1), u2])])]),
        ]),
        InvFun::Acosh => ctx.add(&[
            ctx.mul(&[u, fu]),
            ctx.mul(&[ctx.num(-1), ctx.fun("sqrt", &[ctx.add(&[ctx.num(-1), u2])])]),
        ]),
        InvFun::Atanh => ctx.add(&[
            ctx.mul(&[u, fu]),
            ctx.mul(&[rat_atom(ctx, 1, 2), ctx.fun("log", &[one_minus_u2()])]),
        ]),
    };
    let mut out = rest;
    out.push(inv(ctx, slope));
    out.push(core);
    Some(normalize(ctx, ctx.mul(&out)))
}

// =========================================================================
// Combined entry
// =========================================================================

/// Try the inverse-trig/hyperbolic mechanisms in cost order: bare linear
/// forms (M3), kernel-derivative powers (M1), then the inv-hyp
/// substitution (M2).
pub(crate) fn integrate_inverse_trig<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    integrate_bare_linear(ctx, expr, var)
        .or_else(|| integrate_kernel_power(ctx, expr, var))
        .or_else(|| integrate_invhyp_subst(ctx, expr, var))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Parse and normalize (mirroring what the chain feeds the mechanism).
    fn parse_norm<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        normalize(ctx, ocas_parse::parse(ctx, s).unwrap())
    }

    /// Numeric f64 evaluator for test verification (handles every operator
    /// the M1–M3 results produce).
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
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "asin" => v.asin(),
                    "acos" => v.acos(),
                    "atan" => v.atan(),
                    "asinh" => v.asinh(),
                    "acosh" => v.acosh(),
                    "atanh" => v.atanh(),
                    _ => return None,
                })
            }
        }
    }

    /// Run the combined entry, require `Some` without residue, and check
    /// `diff(result) == integrand` numerically at the sample points.
    /// Samples must keep radicands positive and `|arg|` inside the domain
    /// of the asin/acos/atanh families.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_inverse_trig(ctx, integrand, var).expect("mechanism declined");
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

    fn assert_declined<'a>(ctx: &'a AtomArena<'a>, input: &str) {
        let expr = parse_norm(ctx, input);
        assert!(
            integrate_inverse_trig(ctx, expr, Symbol::new("x")).is_none(),
            "expected None for {input}"
        );
    }

    // ------------------------- M1 -------------------------

    #[test]
    fn m1_asinh_kernel_sqrt_fun_form() {
        // ∫ (a + b·asinh(c·x))²/√(1 + c²x²) dx — 1/sqrt spelling.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*asinh(c*x))^2/sqrt(1 + c^2*x^2)");
        let env = [
            (Symbol::new("a"), 1.2),
            (Symbol::new("b"), -0.7),
            (Symbol::new("c"), 1.5),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.5, 0.9]);
    }

    #[test]
    fn m1_asin_kernel_pow_form() {
        // ∫ (a + b·asin(c·x))³·(1 − c²x²)^(−1/2) dx — Pow spelling; |c·x| < 1.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*asin(c*x))^3*(1 - c^2*x^2)^(-1/2)");
        let env = [
            (Symbol::new("a"), 0.8),
            (Symbol::new("b"), 1.1),
            (Symbol::new("c"), 2.0),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.3, 0.45]);
    }

    #[test]
    fn m1_acos_kernel_negative_sign() {
        // k = 1 bare base; acos carries s_F = −1.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*acos(c*x))/sqrt(1 - c^2*x^2)");
        let env = [
            (Symbol::new("a"), 0.9),
            (Symbol::new("b"), -1.2),
            (Symbol::new("c"), 1.6),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.3, 0.55]);
    }

    #[test]
    fn m1_atanh_kernel_reciprocal() {
        // atanh kernel exponent is −1 (no square root); |c·x| < 1.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*atanh(c*x))^2/(1 - c^2*x^2)");
        let env = [
            (Symbol::new("a"), 1.3),
            (Symbol::new("b"), 0.7),
            (Symbol::new("c"), 1.9),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.25, 0.5]);
    }

    #[test]
    fn m1_atan_kernel_symbolic_norm() {
        // R = d + c²·d·x² = d·(1 + c²x²): symbolic normalization s = d is
        // exact under the reciprocal kernel (no sign gate).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*atan(c*x))^4/(d + c^2*d*x^2)");
        let env = [
            (Symbol::new("a"), 1.1),
            (Symbol::new("b"), -0.8),
            (Symbol::new("c"), 1.7),
            (Symbol::new("d"), 2.3),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.7, 1.3]);
    }

    #[test]
    fn m1_acosh_kernel() {
        // E = c²x² − 1 (constant term e0 = −1); samples keep c·x > 1.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*acosh(c*x))^2/sqrt(c^2*x^2 - 1)");
        let env = [
            (Symbol::new("a"), 0.6),
            (Symbol::new("b"), 1.4),
            (Symbol::new("c"), 0.8),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[1.5, 2.0, 3.0]);
    }

    #[test]
    fn m1_numeric_norm_factor() {
        // s = 4 (numeric, positive): the answer picks up 4^(−1/2).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*asinh(c*x))^2/sqrt(4 + 4*c^2*x^2)");
        let env = [
            (Symbol::new("a"), 1.2),
            (Symbol::new("b"), -0.7),
            (Symbol::new("c"), 1.5),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.5, 0.9]);
    }

    // ------------------------- M2 -------------------------

    #[test]
    fn m2_corpus_acosh_poly() {
        // rubi-00110: (d + e·x)·(a + b·acosh(c·x))² — poly × inv-hyp power.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(d + e*x)*(a + b*acosh(c*x))^2");
        let env = [
            (Symbol::new("a"), 0.9),
            (Symbol::new("b"), -1.3),
            (Symbol::new("c"), 0.7),
            (Symbol::new("d"), 1.1),
            (Symbol::new("e"), 0.6),
        ];
        // acosh domain: c·x > 1.
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[1.6, 2.2, 3.0]);
    }

    #[test]
    fn m2_corpus_asinh_sqrt_pow() {
        // rubi-00049: (d + c²·d·x²)^(3/2)·(a + b·asinh(c·x)) — the
        // sqrt-quadratic-power path (symbolic s = d assumed positive).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(d + c^2*d*x^2)^(3/2)*(a + b*asinh(c*x))");
        let env = [
            (Symbol::new("a"), 0.8),
            (Symbol::new("b"), 1.1),
            (Symbol::new("c"), 1.3),
            (Symbol::new("d"), 1.7),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.6, 1.1]);
    }

    #[test]
    fn m2_corpus_acosh_linear_arg_k4() {
        // rubi-00129: (c·e + d·e·x)·(a + b·acosh(c + d·x))^4 — linear arg
        // with intercept; k = 4 exceeds the chain's parts budget, so the
        // closed-form fallback finishes the t³/t⁴ terms.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(c*e + d*e*x)*(a + b*acosh(c + d*x))^4");
        let env = [
            (Symbol::new("a"), 0.7),
            (Symbol::new("b"), -1.1),
            (Symbol::new("c"), 0.4),
            (Symbol::new("d"), 0.9),
            (Symbol::new("e"), 1.4),
        ];
        // acosh domain: c + d·x > 1.
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.8, 1.2, 1.8]);
    }

    #[test]
    fn m2_asinh_linear_arg_bare_power() {
        // Corpus shape: (a + b·asinh(c + d·x))² with no extra factor.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(a + b*asinh(c + d*x))^2");
        let env = [
            (Symbol::new("a"), 1.1),
            (Symbol::new("b"), 0.7),
            (Symbol::new("c"), 0.4),
            (Symbol::new("d"), 1.2),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.8, 1.5]);
    }

    // ------------------------- M3 -------------------------

    #[test]
    fn m3_trig_family_linear_arg() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // u = a + b·x stays in (−1, 1) at the samples (asin/acos domain).
        let env = [(Symbol::new("a"), 0.3), (Symbol::new("b"), 1.2)];
        for input in ["asin(a + b*x)", "acos(a + b*x)", "atan(a + b*x)"] {
            let expr = parse_norm(&ctx, input);
            assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.3, 0.5]);
        }
    }

    #[test]
    fn m3_hyper_family_linear_arg() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Same env keeps |u| < 1 for atanh; asinh is entire.
        let env = [(Symbol::new("a"), 0.3), (Symbol::new("b"), 1.2)];
        for input in ["asinh(a + b*x)", "atanh(a + b*x)"] {
            let expr = parse_norm(&ctx, input);
            assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.3, 0.5]);
        }
    }

    #[test]
    fn m3_acosh_linear_arg_domain() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // u = a + b·x > 1 at the samples (acosh domain).
        let env = [(Symbol::new("a"), 1.4), (Symbol::new("b"), 0.9)];
        let expr = parse_norm(&ctx, "acosh(a + b*x)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.1, 0.4, 0.8]);
    }

    #[test]
    fn m3_constant_multiple() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 0.3), (Symbol::new("b"), 1.2)];
        let expr = parse_norm(&ctx, "q*asin(a + b*x)");
        let mut env2 = env.to_vec();
        env2.push((Symbol::new("q"), 2.5));
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env2, &[0.1, 0.3, 0.5]);
    }

    // ------------------------- declines -------------------------

    #[test]
    fn declines_out_of_scope() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Corpus: k = −4 (negative inv-hyp power — not the kernel form).
        assert_declined(&ctx, "x/acosh(a*x)^4");
        // Corpus: atanh power × matching-ish rational factor — the rational
        // factor (d + c·d·x)^(−1) is not the atanh kernel (1 − c²x²)^(−1).
        assert_declined(&ctx, "(a + b*atanh(c*x))^2/(d + c*d*x)");
        // Corpus: kernel exponent −5/2 and k = −2 — neither M1 nor M2.
        assert_declined(&ctx, "x^2/((c + a^2*c*x^2)^(5/2)*atan(a*x)^2)");
        // Corpus: asin of a quadratic — no mechanism applies.
        assert_declined(&ctx, "a + b*asin(1 + d*x^2)");
        assert_declined(&ctx, "asin(1 + d*x^2)");
        // k = 6 exceeds the M2 cap.
        assert_declined(&ctx, "(d + e*x)*(a + b*acosh(c*x))^6");
        // Numeric negative normalization under an odd square-root power.
        assert_declined(&ctx, "(-4 - 4*c^2*x^2)^(3/2)*(a + b*asinh(c*x))");
        // asinh needs E = 1 + u²; the 1 − c²x² radicand is the wrong kernel.
        assert_declined(&ctx, "(a + b*asinh(c*x))^2/sqrt(1 - c^2*x^2)");
        // atanh is outside M2's asinh/acosh substitution scope.
        assert_declined(&ctx, "(1 + c*x)^3*(a + b*atanh(c*x))^3");
    }
}
