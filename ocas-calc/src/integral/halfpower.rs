//! Half-power front-end for trig radical integrands (0.27.2 Wave C6; affine
//! arguments added in the 0.27.3 "Phase E1" wave).
//!
//! Detects integrands of the shape
//!
//! ```text
//!     C · S(cos u)^{p/2} dx,     u = c + d·x,  p odd,  deg_`cos` S ≤ 2,
//! ```
//!
//! with `C`, `c` and `d` free of `x` and `d ≢ 0`, and rewrites them into the
//! canonical algebraic form
//! `∫ W(z)/√(Q(z)) dz` with `deg Q = 4` and `Q` already in Legendre normal
//! form, then hands the result to [`super::elliptic::integrate_elliptic`] and
//! substitutes back. When the rewrite does not produce a 3rd/4th-degree
//! radicand the module declines, so elementary engines keep first claim.
//!
//! # Substitutions and their derivation
//!
//! **Weierstrass (`t = tan(u/2)`).** With `sin u = 2t/(1+t²)`,
//! `cos u = (1−t²)/(1+t²)` and `du = 2 dt/(1+t²)`, a base that is linear in
//! `cos u` becomes `S = [(c₀+c₁) + (c₀−c₁)t²]/(1+t²)`, so
//! `S^{p/2} du = 2 P(t)^{p/2} dt/(1+t²)^{(p+2)/2}`; clearing denominators
//! gives the quartic-radicand integral `2C·√(P(t)(1+t²))-shaped` forms. This
//! is the classical route and the one the wave brief specifies.
//!
//! **Half-angle refinement.** The Weierstrass quartic `(1+t²)P(t)` is *not*
//! in Legendre normal form, and its reduction needs the complementary-modulus
//! integrals for the `t²` remainder — the `u = t/√(1+t²)` composition that the
//! elliptic engine knows how to perform anyway. The engine therefore applies
//! the same substitution `z = t/√(1+t²) = sin(u/2)` first, which turns
//!
//! - `S = c₀ + c₁cos u = (c₀+c₁) − 2c₁sin²(u/2) = A(1 − M z²)` with
//!   `A = c₀+c₁`, `M = 2c₁/A`, and `du = 2 dz/√(1−z²)`;
//! - `S = c₀ + c₂cos²u = (c₀+c₂) − c₂sin²u = A(1 − M z²)` with
//!   `A = c₀+c₂`, `M = c₂/A`, and `du = dz/√(1−z²)` (here `z = sin u`),
//!
//! into the *Legendre* quartic `Q = A(1−z²)(1−M z²)`. Since
//! `sin²(u/2) = t²/(1+t²)` and `sin²u = t²/(1+t²)` hold identically, the
//! emitted amplitude `asin(z)` has derivative exactly `1/2` (resp. `1`) times
//! the required chain factor on every interval where the substitution is
//! invertible. This module emits that refined form; it is the same
//! substitution, just composed in the order that keeps the result in the
//! first/second kind.
//!
//! **Affine argument.** Everything above is written for the symbolic atom
//! `u = c + d·x` rather than for the bare integration variable: the module
//! recognises `u` from the `cos` kernel it finds in the radicand, requires
//! `u = c + d·x` with `c`, `d` free of `x` (numeric *or* symbolic) and
//! `d ≢ 0`, and runs the identical reduction with that atom. The only other
//! change is the chain factor
//!
//! ```text
//!     ∫ C·S(cos(c+d·x))^{p/2} dx = (1/d)·∫ C·S(cos u)^{p/2} du,
//! ```
//!
//! so the emitted antiderivative is scaled by `1/d`; a slope of literally `1`
//! (the bare variable, `x + c`, …) emits no chain factor at all, so the
//! reduction itself is left untouched.
//!
//! **One half power per base.** The non-constant part is folded to a single
//! `S^{p/2}` before anything else. `sqrt` is a head rather than a `Pow`, so
//! `normalize` cannot combine `S^{a}·(√S)^{b}` although its value is
//! `S^{a+b/2}`; corpus integrands reach the front-end in exactly that shape
//! (`cos u·(√(cos u))^{−1} = √(cos u)`, `cos u·√(cos u) = cos^{3/2}u`), and
//! the fold is what lets them take the same reduction as a literal `S^{p/2}`.
//! Products of integer/half-integer powers of two *different* bases still
//! decline.
//!
//! Both branches require `S` to depend on `cos u` only (no `sin u` term and,
//! for the quadratic branch, no `cos u` term): a genuine `sin u` component
//! makes the radicand a general quartic whose reduction needs a Möbius
//! normalisation this wave does not implement, and `cos u` together with
//! `cos²u` makes the Legendre form unreachable by either half-angle.
//!
//! # Known gaps (declined, never guessed)
//!
//! - Kernels that are not a polynomial in `cos u`: `sec`, `tan`, `cot`,
//!   `tanh`, `coth`, `sinh` families, and any `S` mixing `sin u`/`cos u`.
//! - Prefactors that depend on `x` (`exp(x)·√(cos x)`, `sin(x)·√(S)`).
//! - Products of two *different* radical bases (`√(1+cos u)·√(cos u)`) and
//!   `|p| > 5`.
//! - Arguments that are not affine in the integration variable
//!   (`cos(x²)`, `cos(sin x)`, `cos(exp x)`), and radicands that are a
//!   polynomial in `cos u` only after a multiple-angle expansion
//!   (`cos 2x = 2cos²x − 1`) — the module does not expand.
//! - `S` whose substituted leading coefficient `A` vanishes identically
//!   (a degenerate, non-squarefree radicand).
//!
//! # Sheet factor (why the emitted form differentiates back)
//!
//! The engine's radicand is `√(1 − z²)` and `z` is inverted by `asin` on its
//! principal branch. Because `dz/du = ½cos(u/2)` (linear-`S` branch) and
//! `dz/du = cos u` (quadratic-`S` branch) are *signed*, the engine's answer
//! differentiates to `sign(cos(u/2))` (resp. `sign(cos u)`) times the
//! integrand. Multiplying the whole antiderivative by that sign — written
//! `cos(u/2)·(1 − sin²(u/2))^{−1/2}` resp. `cos u·(1 − sin²u)^{−1/2}` so no
//! `abs` head is needed — restores the correct sheet on every interval.
//!
//! The 0.27.2 module emitted that factor only on the `z = sin u` branch and
//! claimed the other one solely on `|u| < π`, where the factor is identically
//! `1`. An affine argument puts `u = c + d·x` on an arbitrary interval —
//! whose width depends on symbols whose sign is unknown — so the factor is now
//! emitted on **both** branches. On `|u| < π` it evaluates to exactly `1`, so
//! the bare-variable output is numerically unchanged; beyond that it is what
//! makes the emitted form an antiderivative at all.
//!
//! # Branch / parameter caveat (important)
//!
//! The emitted closed form is real — and therefore the real antiderivative —
//! **only when the radicand's leading factor `A = c₀ + c₁` (respectively
//! `c₀ + c₂`) is positive**, because it carries `1/√A`. If a parameter
//! specialization makes `A < 0`, the same expression is still a formal
//! complex antiderivative, but the branch matching the real integrand is the
//! *complementary* one and would need the imaginary-modulus transformation
//! (`F(φ, m) ↦ F(φ', 1/m)/√m` and the matching `E`/`Π` identities), which is
//! not implemented. The sign of a symbolic `A` is not decidable here, so the
//! form is emitted regardless; the corpus oracle's deterministic dummy values
//! for `a` and `b` (see the `harness_parity_for_cos_family` test) happen to
//! give `a + b < 0`, where both this module's evaluator and the harness's
//! report the point as out of domain rather than as a disagreement.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::is_constant;

/// Integrate `C·S(cos(c + d·x))^{p/2}` by rewriting to a quartic radical
/// integrand.
///
/// Returns `None` when the rewrite does not produce a quartic radicand in
/// Legendre normal form.
pub(crate) fn integrate_half_power<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if super::node_count(expr) > 200 {
        return None;
    }
    // Canonicalise first: `normalize` folds `(S^{a})^{n}` for integer `n`, so
    // `1/(S)^{3/2}` reaches this front-end as the single half power `S^{-3/2}`
    // instead of a radical wrapped in a reciprocal.
    let expr = ocas_atom::normalize::normalize(ctx, expr);
    let (coeff, core) = split_constant(ctx, expr, var)?;
    // The non-constant part must be a single half power of one base.
    let (base, p) = single_half_power(core)?;
    if p % 2 == 0 || p.unsigned_abs() > 5 {
        return None;
    }
    // `S` must be a polynomial of degree ≤ 2 in one `cos(u)` kernel whose
    // argument is affine in the integration variable; `u` is that argument
    // (the bare variable is the special case `u = 0 + 1·var`) and `inv_d` the
    // chain factor `1/d` that the affine substitution adds.
    let (base_atom, inv_d, s) = find_cos_kernel(ctx, base, var)?;
    let c0 = s.first().copied().unwrap_or_else(|| ctx.num(0));
    let c1 = s.get(1).copied().unwrap_or_else(|| ctx.num(0));
    let c2 = s.get(2).copied().unwrap_or_else(|| ctx.num(0));
    let zero = |a: Atom<'a>| is_zero(ctx, a);
    // Pick the half-angle branch: linear in `cos u`, or pure `cos²u`.
    //
    // The base is written as `A − β z²` (a polynomial: no compound
    // denominator), because `S = A(1 − M z²)` with `M = β/A` is the same
    // value; the engine then sees a polynomial radicand and the Legendre
    // normal form drops out directly.
    //
    // `sign_arg` is the angle whose cosine gives the sheet factor: `u/2` on
    // the `z = sin(u/2)` branch and `u` on the `z = sin u` branch, i.e. the
    // argument of `z_repl` itself.
    let (a_coef, beta, scale, z_repl, sign_arg) = if zero(c2) {
        if zero(c1) {
            // Constant base: no radical, elementary engines own this.
            return None;
        }
        let a = cz(ctx, ctx.add(&[c0, c1]));
        if zero(a) {
            return None;
        }
        let beta = cz(ctx, ctx.mul(&[ctx.num(2), c1]));
        let half_var = ctx.mul(&[base_atom, ctx.pow(ctx.num(2), ctx.num(-1))]);
        (a, beta, ctx.num(2), ctx.fun("sin", &[half_var]), half_var)
    } else {
        if !zero(c1) {
            // A `cos u` term together with `cos²u`: neither half-angle reaches
            // Legendre normal form.
            return None;
        }
        let a = cz(ctx, ctx.add(&[c0, c2]));
        if zero(a) {
            return None;
        }
        (a, c2, ctx.num(1), ctx.fun("sin", &[base_atom]), base_atom)
    };
    let z = fresh_symbol(expr, var)?;
    let zvar = ctx.var(z.as_str());
    let radicand = ctx.add(&[
        a_coef,
        ctx.mul(&[ctx.num(-1), beta, ctx.pow(zvar, ctx.num(2))]),
    ]);
    let one_minus_z2 = ctx.add(&[
        ctx.num(1),
        ctx.mul(&[ctx.num(-1), ctx.pow(zvar, ctx.num(2))]),
    ]);
    let rewritten = ctx.mul(&[
        coeff,
        scale,
        ctx.pow(radicand, half_exp(ctx, p)),
        ctx.pow(one_minus_z2, half_exp(ctx, -1)),
    ]);
    let g = super::elliptic::integrate_elliptic(ctx, rewritten, z)?;
    let back = super::replace_symbol(ctx, g, z, z_repl);
    // Sheet factor, emitted on *both* branches (see the module docs).
    //
    // `dz/du = ½cos(u/2)` (resp. `cos u`) is *signed*, so the engine's
    // `√(1 − z²) = |cos(u/2)|` (resp. `|cos u|`) describes the wrong sheet
    // wherever that cosine is negative; multiplying the whole antiderivative
    // by `sign(cos(sign_arg)) = cos(sign_arg)/|cos(sign_arg)|` — written here
    // as `cos(sign_arg)·(1 − sin²(sign_arg))^{-1/2}`, and `1 − sin²(sign_arg)`
    // is exactly `1 − z_repl²`, so no `abs` head is needed — restores the
    // correct sheet on every interval.
    let one_minus_sin2 = ctx.add(&[
        ctx.num(1),
        ctx.mul(&[ctx.num(-1), ctx.pow(z_repl, ctx.num(2))]),
    ]);
    let back = ctx.mul(&[
        ctx.fun("cos", &[sign_arg]),
        ctx.pow(one_minus_sin2, half_exp(ctx, -1)),
        back,
    ]);
    // Chain factor of the affine substitution: `du = d·dx`.
    let back = match inv_d {
        Some(inv_d) => ctx.mul(&[inv_d, back]),
        None => back,
    };
    Some(ocas_atom::normalize::normalize(ctx, back))
}

// =========================================================================
// Helpers
// =========================================================================

/// `2^{-1}`.
fn half<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    ctx.pow(ctx.num(2), ctx.num(-1))
}

/// The atom `p/2`.
fn half_exp<'a>(ctx: &'a AtomArena<'a>, p: i64) -> Atom<'a> {
    ctx.mul(&[ctx.num(p), half(ctx)])
}

fn cz<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    crate::ode::util::collect_terms(ctx, a)
}

fn is_zero<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> bool {
    matches!(cz(ctx, a).node(), AtomNode::Num(0))
}

/// Locate the `cos(u)` kernel of a radicand and return everything the
/// reduction needs: the argument atom `u`, the chain factor `1/d` (or `None`
/// when `d ≡ 1`, the bare-variable case), and the coefficients `[c₀, c₁, c₂]`
/// of `radicand` as a polynomial in `cos(u)`.
///
/// Requires `u = c + d·var` to be **affine** in the integration variable with
/// `c`, `d` free of `var` and `d ≢ 0`. Candidates are `cos` kernels of the
/// radicand, outermost first; a candidate is accepted only once
/// [`poly_in_base`] confirms that the whole radicand really is a polynomial of
/// degree ≤ 2 in it (with `var`-free coefficients), so the accepted shape is
/// always a genuine one — never a guess.
fn find_cos_kernel<'a>(
    ctx: &'a AtomArena<'a>,
    radicand: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Option<Atom<'a>>, Vec<Atom<'a>>)> {
    let x = ctx.var(var.as_str());
    let mut candidates: Vec<Atom<'a>> = Vec::new();
    collect_cos(radicand, &mut candidates);
    for candidate in candidates {
        let AtomNode::Fun(_, args) = candidate.node() else {
            continue;
        };
        if args.len() != 1 {
            continue;
        }
        let arg = args[0];
        // `arg` must be `c + d·var`: degree ≤ 1 in the bare variable with
        // `var`-free coefficients.
        let Some(slope) = poly_in_base(ctx, arg, x, var, 1) else {
            continue;
        };
        let Some(d) = slope.get(1).copied() else {
            continue; // Constant argument: `cos` is not a kernel of `x`.
        };
        if is_zero(ctx, d) {
            continue;
        }
        let Some(coeffs) = poly_in_base(ctx, radicand, candidate, var, 2) else {
            continue;
        };
        let inv_d = if matches!(d.node(), AtomNode::Num(1)) {
            None
        } else {
            Some(ctx.pow(d, ctx.num(-1)))
        };
        return Some((arg, inv_d, coeffs));
    }
    None
}

/// Collect every distinct `cos(...)` atom in `expr`, outermost first.
///
/// The order makes the search deterministic; the first candidate that passes
/// the affine test *and* the polynomial test wins.
fn collect_cos<'a>(expr: Atom<'a>, out: &mut Vec<Atom<'a>>) {
    if let AtomNode::Fun(name, args) = expr.node() {
        if name.as_str() == "cos" && !out.contains(&expr) {
            out.push(expr);
        }
        for a in args.iter() {
            collect_cos(*a, out);
        }
        return;
    }
    match expr.node() {
        AtomNode::Pow(b, e) => {
            collect_cos(*b, out);
            collect_cos(*e, out);
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_cos(*a, out);
            }
        }
        AtomNode::Num(_) | AtomNode::Var(_) | AtomNode::Fun(_, _) => {}
    }
}

/// Split a product into `(constant part, non-constant part)`.
fn split_constant<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let (c, nc): (Vec<Atom<'a>>, Vec<Atom<'a>>) =
        factors.into_iter().partition(|f| is_constant(*f, var));
    if nc.is_empty() {
        return None;
    }
    let core = if nc.len() == 1 { nc[0] } else { ctx.mul(&nc) };
    let coeff = match c.len() {
        0 => ctx.num(1),
        1 => c[0],
        _ => ctx.mul(&c),
    };
    Some((coeff, core))
}

/// Parse a rational exponent atom `p/q` (same shapes as the elliptic engine).
fn exp_fraction<'a>(exp: Atom<'a>) -> Option<(i64, i64)> {
    match exp.node() {
        AtomNode::Num(n) => Some((*n, 1)),
        AtomNode::Pow(b, e) => {
            if let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
                && *ee == -1
                && *bb > 0
            {
                return Some((1, *bb));
            }
            None
        }
        AtomNode::Mul(args) => {
            let mut num: Option<i64> = None;
            let mut den: Option<i64> = None;
            for a in args.iter() {
                match a.node() {
                    AtomNode::Num(n) => {
                        if num.is_some() {
                            return None;
                        }
                        num = Some(*n);
                    }
                    AtomNode::Pow(b, e) => {
                        if let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
                            && *ee == -1
                            && *bb > 0
                        {
                            if den.is_some() {
                                return None;
                            }
                            den = Some(*bb);
                        } else {
                            return None;
                        }
                    }
                    _ => return None,
                }
            }
            match (num, den) {
                (Some(p), Some(q)) => Some((p, q)),
                (Some(p), None) => Some((p, 1)),
                (None, Some(q)) => Some((1, q)),
                (None, None) => None,
            }
        }
        _ => None,
    }
}

/// Structural occurrence test.
fn contains_atom<'a>(expr: Atom<'a>, needle: Atom<'a>) -> bool {
    if expr == needle {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Pow(b, e) => contains_atom(*b, needle) || contains_atom(*e, needle),
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|a| contains_atom(*a, needle))
        }
    }
}

/// A substitution variable that collides neither with `var` nor with any
/// symbol already present in `expr`.
fn fresh_symbol<'a>(expr: Atom<'a>, var: Symbol) -> Option<Symbol> {
    for name in ["_z", "_s", "_w", "_v"] {
        let s = Symbol::new(name);
        if s != var && !super::contains_symbol(expr, s) {
            return Some(s);
        }
    }
    None
}

/// Fold the non-constant part into the single half power `S^{p/2}` it is,
/// or return `None` when it is not a product of integer/half-integer powers of
/// one common base.
///
/// `normalize` folds `(S^{a})^{n}` for integer `n` but `sqrt` is a *head*, not
/// a `Pow`, so nothing folds `S^{a}·(√S)^{b}` into `S^{a+b/2}`. Corpus shapes
/// reach this front-end exactly that way — `cos u·(√(cos u))^{−1}` (the
/// distributed form of `B·cos u/√(cos u)`) is `√(cos u)`, and
/// `cos u·√(cos u)` is `cos^{3/2}u` — so the fold is what lets them take the
/// same reduction as a literal `S^{p/2}`. The three shapes the 0.27.2 module
/// accepted (`sqrt(S)`, `S^{p/2}`, `sqrt(S)^{m}`) map to exactly the same
/// `(S, p)` as before, so nothing that solved then changes.
///
/// The fold is the principal branch (`√S^n·S^m = S^{n/2+m}`), which is the
/// same branch the module's emitted form already claims: where the integrand
/// is real (`S ≥ 0`) it is exact, and where `S < 0` both the original and the
/// folded form are formal complex expressions (see the module-level caveat).
fn single_half_power<'a>(core: Atom<'a>) -> Option<(Atom<'a>, i64)> {
    let factors: Vec<Atom<'a>> = match core.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![core],
    };
    let mut base: Option<Atom<'a>> = None;
    let mut p: i64 = 0;
    for f in factors {
        let (b, e2) = half_power_factor(f)?;
        match base {
            None => base = Some(b),
            Some(prev) if prev == b => {}
            // Two different bases: no single `S^{p/2}` describes the product.
            Some(_) => return None,
        }
        p = p.checked_add(e2)?;
    }
    let base = base?;
    // An even numerator is an integer power: there is no radical to reduce.
    if p % 2 == 0 {
        return None;
    }
    Some((base, p))
}

/// `(base, 2·exponent)` for one multiplicative factor, when its exponent is an
/// integer or a half-integer; `sqrt(S)` counts as `S^{1/2}`.
fn half_power_factor<'a>(f: Atom<'a>) -> Option<(Atom<'a>, i64)> {
    match f.node() {
        // A bare number cannot carry the base; coefficients were already
        // split off by `split_constant`.
        AtomNode::Num(_) => None,
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
            {
                let (m, q) = exp_fraction(*e)?;
                if q != 1 {
                    return None;
                }
                return Some((args[0], m));
            }
            let (num, q) = exp_fraction(*e)?;
            match q {
                1 => Some((*b, num.checked_mul(2)?)),
                2 => Some((*b, num)),
                _ => None,
            }
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            Some((args[0], 1))
        }
        // Anything else carries exponent 1, i.e. `2·exponent = 2`.
        _ => Some((f, 2)),
    }
}

/// Coefficients `[c₀, …, c_n]` of `expr` as a polynomial in `base`, or `None`
/// when `expr` is not such a polynomial.
///
/// Coefficients must be free of `var`; this is what rejects `√(sin u)` (whose
/// "coefficient" would be `sin u` itself) and any other base that merely
/// contains `var` without being built from `base`.
fn poly_in_base<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    base: Atom<'a>,
    var: Symbol,
    max_deg: usize,
) -> Option<Vec<Atom<'a>>> {
    let expanded = cz(ctx, expr);
    let terms: Vec<Atom<'a>> = match expanded.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expanded],
    };
    let mut out: Vec<Atom<'a>> = Vec::new();
    for t in terms {
        let (deg, coeff) = monomial(ctx, t, base, max_deg)?;
        // Coefficients must not depend on the integration variable: this is
        // what rejects a base that merely *contains* `var` (e.g. `sin x`
        // against the base `cos x`).
        if !is_constant(coeff, var) {
            return None;
        }
        if out.len() <= deg {
            out.resize(deg + 1, ctx.num(0));
        }
        out[deg] = cz(ctx, ctx.add(&[out[deg], coeff]));
    }
    while let Some(last) = out.last().copied() {
        if is_zero(ctx, last) {
            out.pop();
        } else {
            break;
        }
    }
    Some(out)
}

/// `(degree in base, coefficient)` for one additive term.
fn monomial<'a>(
    ctx: &'a AtomArena<'a>,
    t: Atom<'a>,
    base: Atom<'a>,
    max_deg: usize,
) -> Option<(usize, Atom<'a>)> {
    if t == base {
        return Some((1, ctx.num(1)));
    }
    match t.node() {
        AtomNode::Num(_) => Some((0, t)),
        AtomNode::Var(_) => {
            if contains_atom(t, base) {
                None
            } else {
                Some((0, t))
            }
        }
        AtomNode::Fun(_, _) => {
            if contains_atom(t, base) {
                None
            } else {
                Some((0, t))
            }
        }
        AtomNode::Add(_) => None,
        AtomNode::Pow(b, e) => {
            if *b == base {
                if let AtomNode::Num(n) = e.node() {
                    let n = usize::try_from(*n).ok()?;
                    if n <= max_deg {
                        return Some((n, ctx.num(1)));
                    }
                }
                return None;
            }
            if contains_atom(t, base) {
                None
            } else {
                Some((0, t))
            }
        }
        AtomNode::Mul(args) => {
            let mut deg = 0usize;
            let mut coeffs: Vec<Atom<'a>> = Vec::with_capacity(args.len());
            for a in args.iter() {
                let (d, c) = monomial(ctx, *a, base, max_deg)?;
                deg = deg.checked_add(d)?;
                if deg > max_deg {
                    return None;
                }
                coeffs.push(c);
            }
            Some((deg, ctx.mul(&coeffs)))
        }
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::super::elliptic::testnum::{diff_local, eval_f64};
    use super::*;
    use ocas_core::arena::Arena;

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, s).unwrap_or_else(|e| panic!("parse {s}: {e:?}"))
    }

    /// Run the front-end, require a radical-free result, and check
    /// `d/dx result == integrand` numerically at the sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        src: &str,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let var = Symbol::new("x");
        let integrand = parse(ctx, src);
        let result =
            integrate_half_power(ctx, integrand, var).unwrap_or_else(|| panic!("declined: {src}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue for {src}: {result}"
        );
        let d = diff_local(ctx, result, var);
        let mut checked = 0usize;
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = match eval_f64(d, &env) {
                Some(v) => v,
                None => continue,
            };
            let rhs = eval_f64(integrand, &env).expect("eval integrand");
            let tol = 1e-9 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at x={xv}: diff={lhs} integrand={rhs} (src: {src}, result: {result})"
            );
            checked += 1;
        }
        assert!(checked >= 2, "only {checked} usable samples for {src}");
    }

    fn declines<'a>(ctx: &'a AtomArena<'a>, src: &str) {
        let e = parse(ctx, src);
        let var = Symbol::new("x");
        assert!(
            integrate_half_power(ctx, e, var).is_none(),
            "expected decline for {src}"
        );
    }

    #[test]
    fn linear_in_cos_first_kind() {
        // ∫dx/√(a+b·cos x) = (2/√(a+b))·F(x/2, 2b/(a+b))
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(a+b*cos(x))",
            &env,
            &[-2.0, -1.0, 0.4, 1.2, 2.4],
        );
    }

    #[test]
    fn linear_in_cos_second_kind() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        // ∫√(a+b·cos x) dx = 2√(a+b)·E(x/2, 2b/(a+b))
        assert_antiderivative_num(&ctx, "sqrt(a+b*cos(x))", &env, &[-2.0, -0.8, 0.5, 1.5, 2.5]);
    }

    #[test]
    fn linear_in_cos_hermite() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        assert_antiderivative_num(
            &ctx,
            "1/(a+b*cos(x))^(3/2)",
            &env,
            &[-2.2, -0.6, 0.3, 1.4, 2.6],
        );
    }

    #[test]
    fn quadratic_in_cos_first_kind() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        // ∫dx/√(a+b·cos²x) = F(x, b/(a+b))/√(a+b)
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(a+b*cos(x)^2)",
            &env,
            &[-2.4, -1.2, -0.4, 0.6, 1.9, 2.6],
        );
    }

    #[test]
    fn quadratic_in_cos_second_kind() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        // ∫√(a+b·cos²x) dx = √(a+b)·E(x, b/(a+b))
        assert_antiderivative_num(
            &ctx,
            "sqrt(a+b*cos(x)^2)",
            &env,
            &[-2.4, -1.1, -0.3, 0.7, 1.8, 2.7],
        );
    }

    #[test]
    fn quadratic_in_cos_hermite() {
        // The brief's headline case: ∫dx/(a+b·cos²x)^{3/2}
        //   = E(x, M)/(a√(a+b)) − (b/a)·sin x cos x / ((a+b)√(a+b·cos²x)),
        //   M = b/(a+b).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        assert_antiderivative_num(
            &ctx,
            "1/(a+b*cos(x)^2)^(3/2)",
            &env,
            &[-2.4, -1.3, -0.3, 0.8, 1.7, 2.8],
        );
    }

    #[test]
    fn numeric_coefficients_and_extra_powers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_antiderivative_num(&ctx, "1/sqrt(2+cos(x))", &[], &[-2.0, -0.5, 0.6, 2.0]);
        assert_antiderivative_num(&ctx, "1/(2+cos(x)^2)^(3/2)", &[], &[-2.2, -0.7, 0.5, 2.2]);
        assert_antiderivative_num(&ctx, "3/sqrt(2+cos(x))", &[], &[-1.8, -0.4, 0.9, 2.1]);
        // Symbolic coefficient on the second-kind branch.
        let env = [(Symbol::new("a"), 3.0), (Symbol::new("b"), 1.0)];
        assert_antiderivative_num(&ctx, "sqrt(a+b*cos(x)^2)", &env, &[-2.1, -0.8, 0.6, 2.2]);
    }

    // ---------------------------------------------------------------------
    // Affine arguments (0.27.3 "Phase E1"): `u = c + d·x`
    // ---------------------------------------------------------------------

    /// The sample grid is chosen so that `u = c + d·x` crosses **both** `±π`:
    /// that is exactly where `dz/du = ½cos(u/2)` changes sign, i.e. the regime
    /// the bare-variable tests cannot reach (`|x| < π` there). Without the
    /// sheet factor these cases differentiate to `±` the integrand and fail.
    const AFFINE_ENV: [f64; 4] = [2.0, 1.0, 0.3, 1.6];
    const AFFINE_SAMPLES: [f64; 6] = [-2.4, -1.2, -0.3, 0.5, 1.3, 2.5];

    fn affine_env() -> [(Symbol, f64); 4] {
        [
            (Symbol::new("a"), AFFINE_ENV[0]),
            (Symbol::new("b"), AFFINE_ENV[1]),
            (Symbol::new("c"), AFFINE_ENV[2]),
            (Symbol::new("d"), AFFINE_ENV[3]),
        ]
    }

    #[test]
    fn affine_linear_in_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = affine_env();
        // ∫dx/√(a+b·cos(c+d·x)) = (2/(d√(a+b)))·F(u/2, 2b/(a+b))
        assert_antiderivative_num(&ctx, "1/sqrt(a+b*cos(c+d*x))", &env, &AFFINE_SAMPLES);
        assert_antiderivative_num(&ctx, "sqrt(a+b*cos(c+d*x))", &env, &AFFINE_SAMPLES);
        assert_antiderivative_num(&ctx, "1/(a+b*cos(c+d*x))^(3/2)", &env, &AFFINE_SAMPLES);
        // The `z = sin(u/2)` sheet only flips past `|u| = π`, so a run that
        // straddles it is the guard for the linear-branch factor.
        let u: Vec<f64> = AFFINE_SAMPLES
            .iter()
            .map(|x| AFFINE_ENV[2] + AFFINE_ENV[3] * x)
            .collect();
        assert!(u[0] < -std::f64::consts::PI && u[u.len() - 1] > std::f64::consts::PI);
    }

    #[test]
    fn affine_quadratic_in_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = affine_env();
        assert_antiderivative_num(&ctx, "1/sqrt(a+b*cos(c+d*x)^2)", &env, &AFFINE_SAMPLES);
        assert_antiderivative_num(&ctx, "sqrt(a+b*cos(c+d*x)^2)", &env, &AFFINE_SAMPLES);
        assert_antiderivative_num(&ctx, "1/(a+b*cos(c+d*x)^2)^(3/2)", &env, &AFFINE_SAMPLES);
    }

    /// Numeric slopes, unit slope with a shift, a half slope and a negative
    /// slope.
    #[test]
    fn affine_numeric_and_unit_slopes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Slope 2: the emitted chain factor is the number `1/2`.
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(3+cos(2*x))",
            &[],
            &[-2.6, -1.1, 0.0, 0.7, 1.4],
        );
        // `d ≡ 1` with a constant shift: no chain factor is emitted.
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(3+cos(x+1))",
            &[],
            &[-2.6, -1.0, 0.2, 1.1, 2.3],
        );
        // Half slope (`x/2`): the chain factor is the number `2`.
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(3+cos(1+x/2))",
            &[],
            &[-4.0, -2.0, 0.0, 2.0, 4.0, 6.0],
        );
        // Negative slope.
        assert_antiderivative_num(
            &ctx,
            "1/sqrt(3+cos(1-2*x))",
            &[],
            &[-1.5, -0.4, 0.6, 1.9, 2.4],
        );
    }

    /// The two corpus shapes with `S = cos u` exactly (rubi-00377 and its
    /// `p = +1` twin). `A = 1`, `β = 2` fixes the emitted `m = 2`, but the
    /// evaluator's `q = 1 − m·sin²φ = cos u`, so every sample where the
    /// integrand is real is also evaluable and the derivative check is real.
    #[test]
    fn affine_pure_cos_both_kinds() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 0.3), (Symbol::new("b"), 1.6)];
        // `u = 0.3 + 1.6x` stays inside `(−π/2, π/2)` on these samples, which
        // is the integrand's own reality domain for `√(cos u)`.
        let samples = [-1.0, -0.6, -0.2, 0.2, 0.6];
        assert_antiderivative_num(&ctx, "1/cos(a+b*x)^(1/2)", &env, &samples);
        assert_antiderivative_num(&ctx, "cos(a+b*x)^(1/2)", &env, &samples);
        // Scaling the kernel keeps the same branch.
        let scaled = [
            (Symbol::new("a"), 0.3),
            (Symbol::new("b"), 1.6),
            (Symbol::new("c"), 4.0),
        ];
        assert_antiderivative_num(&ctx, "(c*cos(a+b*x))^(1/2)", &scaled, &samples);
    }

    /// Same-base power folding. `sqrt` is a head, so `normalize` leaves
    /// `S^{a}·(√S)^{b}` alone even though its value is `S^{a+b/2}`; the module
    /// folds it, which is what turns the distributed corpus shape
    /// `B·cos u/√(cos u)` into the `√(cos u)` it already knows.
    #[test]
    fn same_base_powers_fold_to_one_half_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 0.3), (Symbol::new("b"), 1.6)];
        let samples = [-1.0, -0.6, -0.2, 0.2, 0.6];
        assert_antiderivative_num(&ctx, "cos(a+b*x)/sqrt(cos(a+b*x))", &env, &samples);
        assert_antiderivative_num(&ctx, "cos(a+b*x)*(cos(a+b*x))^(-1/2)", &env, &samples);
        // Exponents add: `√S·S^{1/2}·S^{-1/2}` is `√S`.
        assert_antiderivative_num(
            &ctx,
            "sqrt(cos(a+b*x))*(cos(a+b*x))^(1/2)*(cos(a+b*x))^(-1/2)",
            &env,
            &samples,
        );
        // `cos u·√(cos u) = cos^{3/2}u` is folded to `p = 3` and then declined
        // downstream, exactly like a literal `(a+b·cos u)^{3/2}`.
        declines(&ctx, "cos(a+b*x)*sqrt(cos(a+b*x))");
        // A product that folds to an *integer* power is not a radical at all.
        declines(&ctx, "cos(a+b*x)*cos(a+b*x)^(-1)");
        // Two different bases still decline.
        declines(&ctx, "sqrt(cos(a+b*x))*sqrt(2+cos(a+b*x))");
        declines(&ctx, "sqrt(cos(a+b*x))*sqrt(sin(a+b*x))");
    }

    /// End-to-end (corpus rubi-00291). The pipeline distributes
    /// `(A + B·cos u)/√(cos u)` into `A·cos^{-1/2}u + B·cos u·(√(cos u))^{-1}`;
    /// the first term was already reducible, the second only becomes the
    /// `√(cos u)` kernel once the same-base fold fires. Before the fold this
    /// case ends in a residual `Integral(...)` and must be counted a fallback.
    #[test]
    fn affine_distributed_corpus_shape_solves() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for src in [
            "(A + B*cos(c+d*x))/sqrt(cos(c+d*x))",
            "cos(c+d*x)/sqrt(cos(c+d*x))",
        ] {
            let integrand = parse(&ctx, src);
            let result = crate::integrate(&ctx, integrand, var);
            let text = result.to_string();
            assert!(!text.contains("Integral("), "{src} fell back: {text}");
        }
    }

    /// Harness-parity check.
    ///
    /// `ocas-tests/src/integral_eval.rs` verifies a solved case with a 5-point
    /// central difference of the antiderivative at the fixed abscissae below,
    /// assigning each free parameter a deterministic dummy value keyed by its
    /// name (hash → `TABLE`, random sign). This test replicates that oracle
    /// exactly over the *whole* pipeline (`ocas_calc::integrate`) and asserts
    /// the two properties that matter:
    ///
    /// 1. With the harness's own dummy values for `a`/`b` the antiderivative
    ///    may land outside the evaluated branch (`a + b < 0` makes the emitted
    ///    `(a+b)^{-1/2}` complex, and the harness's `m = b/(a+b)` then exceeds
    ///    its `|m| ≤ 0.9` window), but it must **never** disagree where both
    ///    sides are evaluable.
    /// 2. On the same abscissae — which include `1.91` and `2.53`, i.e. past
    ///    `π/2` — with a parameter regime where the antiderivative is real,
    ///    every usable sample must agree to `1e-5` relative. This is the
    ///    regression guard for the `z = sin u` branch correction.
    #[test]
    fn harness_parity_for_cos_family() {
        const TABLE: [f64; 12] = [
            2.0,
            3.0,
            5.0,
            7.0,
            0.5,
            1.5,
            0.25,
            11.0,
            1.0 / 3.0,
            4.0,
            6.0,
            0.75,
        ];
        const SAMPLES: [f64; 8] = [-1.7, -0.9, -0.37, 0.31, 0.77, 1.23, 1.91, 2.53];
        fn param_value(name: &str) -> f64 {
            let mut h: u64 = 0xcbf2_9ce4_8422_2325;
            for b in name.bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
            let v = TABLE[(h % 12) as usize];
            if (h >> 8) & 1 == 0 { v } else { -v }
        }
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let cases = [
            "1/(a + b*cos(x)^2)^(3/2)",
            "1/sqrt(a + b*cos(x))",
            "1/sqrt(a + b*cos(x)^2)",
            "sqrt(a + b*cos(x))",
            "sqrt(a + b*cos(x)^2)",
        ];
        for mode in ["harness dummy parameters", "real regime a=2 b=1"] {
            for src in cases {
                let env = if mode.starts_with("harness") {
                    [
                        (Symbol::new("a"), param_value("a")),
                        (Symbol::new("b"), param_value("b")),
                    ]
                } else {
                    [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)]
                };
                let integrand = parse(&ctx, src);
                let result = crate::integrate(&ctx, integrand, var);
                let text = result.to_string();
                assert!(!text.contains("Integral("), "{src} fell back: {text}");
                let (mut checked, mut worst) = (0usize, 0.0f64);
                for &x in SAMPLES.iter() {
                    let h = 1e-4 * x.abs().max(1.0);
                    let eval_at = |t: f64| {
                        let mut e = env.to_vec();
                        e.push((var, t));
                        eval_f64(result, &e)
                    };
                    let (fm2, fm1, fp1, fp2) = (
                        eval_at(x - 2.0 * h),
                        eval_at(x - h),
                        eval_at(x + h),
                        eval_at(x + 2.0 * h),
                    );
                    let mut e = env.to_vec();
                    e.push((var, x));
                    let rhs = eval_f64(integrand, &e);
                    if let (Some(a), Some(b), Some(c), Some(d), Some(f)) = (fm2, fm1, fp1, fp2, rhs)
                    {
                        let deriv = (-d + 8.0 * c - 8.0 * b + a) / (12.0 * h);
                        if deriv.is_finite() {
                            checked += 1;
                            worst = worst.max((deriv - f).abs() / f.abs().max(1.0));
                        }
                    }
                }
                if mode.starts_with("harness") {
                    // Outside the evaluated branch the harness's oracle reports
                    // Domain for every sample (Indeterminate), never a
                    // mismatch; assert only the "no disagreement" half.
                    assert!(worst <= 1e-5, "{src} [{mode}]: worst rel {worst:e}");
                } else {
                    assert!(
                        checked >= 2,
                        "{src} [{mode}]: only {checked} usable samples"
                    );
                    assert!(
                        worst <= 1e-5,
                        "{src} [{mode}]: worst rel {worst:e} ({text})"
                    );
                    println!("{src}: checked={checked} worst_rel={worst:e}");
                }
            }
        }
    }

    #[test]
    fn declines_unsupported_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // `sin` component: the radicand is a general quartic.
        declines(&ctx, "sqrt(sin(x))");
        declines(&ctx, "1/sqrt(a+b*sin(x))");
        // `cos` and `cos²` together: neither half-angle reaches Legendre form.
        declines(&ctx, "1/sqrt(a+b*cos(x)+c*cos(x)^2)");
        // Non-polynomial kernels.
        declines(&ctx, "1/(a+b*sec(x))^(5/2)");
        declines(&ctx, "sqrt(tan(x))");
        declines(&ctx, "sqrt(a+b*sinh(x))");
        // Variable-dependent prefactor.
        declines(&ctx, "exp(x)*sqrt(cos(x))");
        declines(&ctx, "sin(x)*sqrt(a+b*cos(x))");
        declines(&ctx, "x/sqrt(a+b*cos(x))");
        // Radical-free / plain algebraic inputs.
        declines(&ctx, "1/sqrt(x^5+1)");
        declines(&ctx, "1/(1+cos(x))");
        declines(&ctx, "cos(x)");
        declines(&ctx, "sqrt(1+x^2)");
        // Higher half powers.
        declines(&ctx, "1/(a+b*cos(x))^(7/2)");
        // Corpus rubi-00305: same `|p| > 5` budget, affine argument.
        declines(&ctx, "1/(b*cos(c+d*x))^(7/2)");
        // Positive third powers leave a non-constant `S(u)√Q` remainder, which
        // the elliptic engine deliberately declines (see its module docs).
        declines(&ctx, "(a+b*cos(x))^(3/2)");
        declines(&ctx, "(a+b*cos(x)^2)^(5/2)");
        // Degenerate base (A ≡ 0): 1 − cos x = 2sin²(x/2) is not our family.
        declines(&ctx, "1/sqrt(1-cos(x))");
        // Arguments that are not affine in the integration variable.
        declines(&ctx, "1/sqrt(a+b*cos(x^2))");
        declines(&ctx, "1/sqrt(a+b*cos(x^2+x))");
        declines(&ctx, "1/sqrt(a+b*cos(sin(x)))");
        declines(&ctx, "1/sqrt(a+b*cos(exp(x)))");
        // Affine radicand that is not a polynomial in `cos u`.
        declines(&ctx, "1/sqrt(a+b*sec(c+d*x))");
        declines(&ctx, "1/sqrt(a+b*cos(c+d*x)+e*cos(c+d*x)^2)");
        declines(&ctx, "1/sqrt(a+b*cos(c+d*x)^(-1))");
        // Multiple different cos kernels: no single polynomial base exists.
        declines(&ctx, "1/sqrt(a+cos(c+d*x)+cos(2*(c+d*x)))");
        // The rule-table engine's resonant shapes: this module must keep
        // declining them, so it can never poach a case another stage owns.
        declines(&ctx, "cos(c+d*x)^2");
        declines(&ctx, "sin(c+d*x)^2");
        declines(&ctx, "sin(c+d*x)*cos(c+d*x)");
        declines(&ctx, "(a*cos(c+d*x)+b*sin(c+d*x))^2");
    }

    #[test]
    fn stress_stable_outcomes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let solid = parse(&ctx, "1/sqrt(a+b*cos(x))");
        let mut first: Option<String> = None;
        for i in 0..300 {
            let r = integrate_half_power(&ctx, solid, var);
            let s = r.map(|a| a.to_string());
            if i == 0 {
                first = s.clone();
            }
            assert_eq!(s, first, "non-deterministic outcome at iteration {i}");
        }
        assert!(first.is_some(), "stable outcome must be a solve");
        let none = parse(&ctx, "sqrt(sin(x))");
        for i in 0..300 {
            assert!(
                integrate_half_power(&ctx, none, var).is_none(),
                "decline not stable at iteration {i}"
            );
        }
    }
}
