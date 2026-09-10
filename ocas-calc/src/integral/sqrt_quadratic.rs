//! General quadratic-radical integration engine for `√(a + b·x + c·x²)`
//! (0.27.1, Phase 1D).
//!
//! Rule-table family G covers the special cases (`sqrt(a+b·x)`,
//! `sqrt(a²±x²)`, `1/sqrt(x²±a²)`); this module is the general engine for
//! an arbitrary quadratic `q = a + b·x + c·x²` with coefficients constant
//! w.r.t. the integration variable (symbolic coefficients allowed where
//! noted). All sub-mechanisms return `None` outside their scope — a wrong
//! answer is never emitted.
//!
//! - **S1 direct forms** (`c ≠ 0`): `∫dx/√q` (log form
//!   `log(2√c·√q + 2c·x + b)/√c` for `c > 0` or symbolic `c`; asin form
//!   `asin((−2c·x − b)/√Δ)/√(−c)` only with numeric evidence `c < 0` and
//!   `Δ = b² − 4ac > 0`), `∫√q dx = (2c·x+b)√q/(4c) − (Δ/8c)·∫dx/√q`,
//!   `∫dx/q^(3/2) = (4c·x + 2b)/((4ac − b²)·√q)`, and the
//!   numerator-polynomial ladder `(p₀ + p₁x + p₂x²)·q^(±1/2)` via
//!   `J_m = ∫x^m/√q`, `I_m = ∫x^m·√q` for `m ≤ 2`.
//! - **S2 reciprocal forms**: `∫dx/(x·√q)` (log form for rational `a > 0`,
//!   asin form for rational `a < 0 ∧ Δ > 0`; symbolic `a` declines) and
//!   `∫dx/(x²·√q) = −√q/(a·x) − (b/2a)·∫dx/(x·√q)`.
//! - **S4 linear denominator**: `∫dx/((d+e·x)·√q)` and `∫√q/(d+e·x) dx`
//!   via `u = d + e·x`, mapping to S2/S1 on the *unscaled*
//!   `q̃(u) = q((u−d)/e)`. Because `q̃(d+e·x) ≡ q`, back-substitution
//!   replaces the `√q̃` node by `√q` directly (no `|e|` factors; symbolic
//!   `d, e` are supported wherever the S2 sign gate on the constant term
//!   of `q̃` is decidable).
//! - **S3 Euler III**: `q` with two distinct rational roots (numeric
//!   evidence only): `t = √(c·(x−r₁)/(x−r₂))` rationalizes any
//!   rational-in-`(x, √q)` integrand; the t-form is reintegrated via
//!   [`integrate_raw`] and declined on any `Integral` residue.
//!
//! Outputs keep their radicals inside `log`/`asin` arguments or under a
//! quotient (S3), neither of which re-matches this module's patterns, so
//! the mechanism is idempotent under chain re-entry. As a bonus the S2/S4
//! formulas remain valid for a linear base `q` (`c = 0`), which the
//! matcher accepts; the plain S1 linear shapes are declined (rule-table
//! territory, family G).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::{quadratic_coeffs, rat_of, rational_sqrt};
use super::{
    contains_integral, gcd_i64, int_pow, integrate_raw, inv, is_constant, linear_form, node_count,
    pick_subst_symbol, rat_atom, replace_symbol,
};

/// Node budget for the input integrand and for substituted t-forms (S3).
const MAX_NODES: usize = 200;

/// A quadratic `q = a + b·x + c·x²` with coefficients constant w.r.t. the
/// integration variable, plus the variable atom and the canonical `q` /
/// `√q` atoms used when emitting results.
struct Quad<'a> {
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
    q: Atom<'a>,
    sqrt_q: Atom<'a>,
}

impl<'a> Quad<'a> {
    /// Parse a radical's inner expression as a quadratic in `var`.
    fn parse(ctx: &'a AtomArena<'a>, inner: Atom<'a>, var: Symbol) -> Option<Quad<'a>> {
        // quadratic_coeffs unwraps sqrt/pow shapes; wrap the bare inner.
        let (c, b, a) = quadratic_coeffs(ctx, ctx.fun("sqrt", &[inner]), var)?;
        Some(Quad::build(ctx, a, b, c, ctx.var(var.as_str())))
    }

    /// Assemble the canonical atoms for known coefficients.
    fn build(
        ctx: &'a AtomArena<'a>,
        a: Atom<'a>,
        b: Atom<'a>,
        c: Atom<'a>,
        x: Atom<'a>,
    ) -> Quad<'a> {
        let q = normalize(
            ctx,
            ctx.add(&[a, ctx.mul(&[b, x]), ctx.mul(&[c, int_pow(ctx, x, 2)])]),
        );
        let sqrt_q = ctx.fun("sqrt", &[q]);
        Quad {
            a,
            b,
            c,
            x,
            q,
            sqrt_q,
        }
    }
}

fn is_zero(e: Atom<'_>) -> bool {
    matches!(e.node(), AtomNode::Num(0))
}

/// `(−1)·a`, normalized so numeric atoms fold (pattern-match ready).
fn neg<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    normalize(ctx, ctx.mul(&[ctx.num(-1), a]))
}

// -------------------------------------------------------------------------
// Exact rational folding: `normalize` folds numeric `Mul` but not a
// mixed-shape `Add` like `1 + 4·9⁻¹`, and the `rat_of` sign gates below
// need canonical atoms.
// -------------------------------------------------------------------------

/// Reduced rational `(p, q)` with `q > 0`; `None` on zero denominator or
/// arithmetic overflow.
fn r_new(n: i64, d: i64) -> Option<(i64, i64)> {
    if d == 0 {
        return None;
    }
    let (n, d) = if d < 0 {
        (n.checked_neg()?, d.checked_neg()?)
    } else {
        (n, d)
    };
    let g = gcd_i64(n, d);
    Some((n / g, d / g))
}

fn r_add(a: (i64, i64), b: (i64, i64)) -> Option<(i64, i64)> {
    let n = a.0.checked_mul(b.1)?.checked_add(b.0.checked_mul(a.1)?)?;
    r_new(n, a.1.checked_mul(b.1)?)
}

fn r_mul(a: (i64, i64), b: (i64, i64)) -> Option<(i64, i64)> {
    r_new(a.0.checked_mul(b.0)?, a.1.checked_mul(b.1)?)
}

fn r_pow(a: (i64, i64), n: i64) -> Option<(i64, i64)> {
    let (p, q) = if n < 0 {
        if a.0 == 0 {
            return None;
        }
        (a.1, a.0)
    } else {
        a
    };
    let k = u32::try_from(n.unsigned_abs()).ok()?;
    r_new(p.checked_pow(k)?, q.checked_pow(k)?)
}

/// Evaluate `atom` as an exact rational (ints, rational products/sums,
/// integer powers); `None` when any non-rational structure occurs.
fn eval_rational(atom: Atom<'_>) -> Option<(i64, i64)> {
    match atom.node() {
        AtomNode::Num(n) => Some((*n, 1)),
        AtomNode::Add(args) => args
            .iter()
            .try_fold((0i64, 1i64), |acc, a| r_add(acc, eval_rational(*a)?)),
        AtomNode::Mul(args) => args
            .iter()
            .try_fold((1i64, 1i64), |acc, a| r_mul(acc, eval_rational(*a)?)),
        AtomNode::Pow(b, e) => {
            let AtomNode::Num(n) = e.node() else {
                return None;
            };
            r_pow(eval_rational(*b)?, *n)
        }
        AtomNode::Var(_) | AtomNode::Fun(_, _) => None,
    }
}

/// Fold `atom` to the canonical rational `p/q` atom when it is fully
/// rational; otherwise return it unchanged (symbolic case — the downstream
/// numeric-evidence gates decline on their own).
fn fold_rational<'a>(ctx: &'a AtomArena<'a>, atom: Atom<'a>) -> Atom<'a> {
    match eval_rational(atom) {
        Some((p, q)) => rat_atom(ctx, p, q),
        None => atom,
    }
}

/// Discriminant `Δ = b² − 4ac`, normalized. Squares are built with `mul`
/// (not `pow`): `normalize` folds numeric `Mul` but not numeric `Pow`,
/// and the `rat_of` gates downstream need folded numerics.
fn delta<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>) -> Atom<'a> {
    fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.add(&[
                ctx.mul(&[qd.b, qd.b]),
                neg(ctx, ctx.mul(&[ctx.num(4), qd.a, qd.c])),
            ]),
        ),
    )
}

/// `√a` as an atom: folded to a rational when `a` is a positive rational
/// square, left as a `sqrt` node otherwise. `None` for negative rationals
/// (imaginary — the caller's branch guard should have caught it).
fn sqrt_emit<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Option<Atom<'a>> {
    if let Some((p, q)) = rat_of(a) {
        if p < 0 {
            return None;
        }
        if let Some((r, s)) = rational_sqrt(p, q) {
            return Some(rat_atom(ctx, r, s));
        }
    }
    Some(ctx.fun("sqrt", &[a]))
}

// =========================================================================
// S1.1: ∫ dx/√q
// =========================================================================

/// `∫ dx/√q`. For `c = 0` (linear `q`, used only by the S2/S4 helpers)
/// returns `2√q/b`. For `c ≠ 0`: asin form with numeric evidence
/// `c < 0 ∧ Δ > 0`, otherwise the log form (standard CAS answer for
/// `c > 0`, also emitted for symbolic `c`).
fn int_inv_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        if is_zero(qd.b) {
            return None;
        }
        return Some(normalize(
            ctx,
            ctx.mul(&[ctx.num(2), qd.sqrt_q, inv(ctx, qd.b)]),
        ));
    }
    if let Some((p, _)) = rat_of(qd.c)
        && p < 0
    {
        let d = delta(ctx, qd);
        let (pd, _) = rat_of(d)?;
        if pd <= 0 {
            return None;
        }
        let sq_d = sqrt_emit(ctx, d)?;
        let sq_nc = sqrt_emit(ctx, neg(ctx, qd.c))?;
        // asin((−2c·x − b)/√Δ) / √(−c)
        let num = normalize(
            ctx,
            ctx.add(&[ctx.mul(&[ctx.num(-2), qd.c, qd.x]), neg(ctx, qd.b)]),
        );
        let arg = normalize(ctx, ctx.mul(&[num, inv(ctx, sq_d)]));
        return Some(normalize(
            ctx,
            ctx.mul(&[inv(ctx, sq_nc), ctx.fun("asin", &[arg])]),
        ));
    }
    let sq_c = sqrt_emit(ctx, qd.c)?;
    // log(2√c·√q + 2c·x + b) / √c
    let arg = normalize(
        ctx,
        ctx.add(&[
            ctx.mul(&[ctx.num(2), sq_c, qd.sqrt_q]),
            ctx.mul(&[ctx.num(2), qd.c, qd.x]),
            qd.b,
        ]),
    );
    Some(normalize(
        ctx,
        ctx.mul(&[inv(ctx, sq_c), ctx.fun("log", &[arg])]),
    ))
}

// =========================================================================
// S1.2/S1.5: the I_m = ∫ x^m·√q ladder (m ≤ 2)
// =========================================================================

/// `I₀ = ∫ √q dx = (2c·x+b)√q/(4c) − (Δ/8c)·J₀`. The `J₀` term is dropped
/// when `Δ = 0` (perfect square), where the first term alone is exact.
fn int_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let t1 = normalize(
        ctx,
        ctx.mul(&[
            normalize(ctx, ctx.add(&[ctx.mul(&[ctx.num(2), qd.c, qd.x]), qd.b])),
            qd.sqrt_q,
            inv(ctx, ctx.mul(&[ctx.num(4), qd.c])),
        ]),
    );
    let d = delta(ctx, qd);
    if is_zero(d) {
        return Some(t1);
    }
    let j0 = int_inv_sqrt(ctx, qd)?;
    let t2 = normalize(
        ctx,
        ctx.mul(&[d, inv(ctx, ctx.mul(&[ctx.num(8), qd.c])), j0]),
    );
    Some(normalize(ctx, ctx.add(&[t1, neg(ctx, t2)])))
}

/// `I₁ = ∫ x·√q dx = q·√q/(3c) − (b/2c)·I₀`.
fn int_x_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>, i0: Option<Atom<'a>>) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let t1 = normalize(
        ctx,
        ctx.mul(&[qd.q, qd.sqrt_q, inv(ctx, ctx.mul(&[ctx.num(3), qd.c]))]),
    );
    if is_zero(qd.b) {
        return Some(t1);
    }
    let t2 = normalize(
        ctx,
        ctx.mul(&[qd.b, inv(ctx, ctx.mul(&[ctx.num(2), qd.c])), i0?]),
    );
    Some(normalize(ctx, ctx.add(&[t1, neg(ctx, t2)])))
}

/// `I₂ = ∫ x²·√q dx = x·q·√q/(4c) − (5b/8c)·I₁ − (a/4c)·I₀`.
fn int_x2_sqrt<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    i0: Option<Atom<'a>>,
) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let t1 = normalize(
        ctx,
        ctx.mul(&[
            qd.x,
            qd.q,
            qd.sqrt_q,
            inv(ctx, ctx.mul(&[ctx.num(4), qd.c])),
        ]),
    );
    let mut terms = vec![t1];
    if !is_zero(qd.b) {
        let i1 = int_x_sqrt(ctx, qd, i0)?;
        terms.push(neg(
            ctx,
            normalize(
                ctx,
                ctx.mul(&[ctx.num(5), qd.b, inv(ctx, ctx.mul(&[ctx.num(8), qd.c])), i1]),
            ),
        ));
    }
    if !is_zero(qd.a) {
        terms.push(neg(
            ctx,
            normalize(
                ctx,
                ctx.mul(&[qd.a, inv(ctx, ctx.mul(&[ctx.num(4), qd.c])), i0?]),
            ),
        ));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// `∫ (p₀ + p₁x + p₂x²)·√q dx` via the `I_m` ladder. Coefficient paths
/// multiply by a zero coefficient are skipped, so e.g. `b = 0` never
/// requires `I₀` (and its `J₀` branch evidence).
fn poly_times_sqrt<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    p: [Atom<'a>; 3],
) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let i0 = int_sqrt(ctx, qd);
    let mut terms: Vec<Atom<'a>> = Vec::new();
    if !is_zero(p[0]) {
        terms.push(ctx.mul(&[p[0], i0?]));
    }
    if !is_zero(p[1]) {
        terms.push(ctx.mul(&[p[1], int_x_sqrt(ctx, qd, i0)?]));
    }
    if !is_zero(p[2]) {
        terms.push(ctx.mul(&[p[2], int_x2_sqrt(ctx, qd, i0)?]));
    }
    if terms.is_empty() {
        return None;
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

// =========================================================================
// S1.4/S1.5: the J_m = ∫ x^m/√q ladder (m ≤ 2)
// =========================================================================

/// `J₁ = ∫ x/√q dx = √q/c − (b/2c)·J₀`.
fn int_x_inv_sqrt<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    j0: Option<Atom<'a>>,
) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let t1 = normalize(ctx, ctx.mul(&[qd.sqrt_q, inv(ctx, qd.c)]));
    if is_zero(qd.b) {
        return Some(t1);
    }
    let t2 = normalize(
        ctx,
        ctx.mul(&[qd.b, inv(ctx, ctx.mul(&[ctx.num(2), qd.c])), j0?]),
    );
    Some(normalize(ctx, ctx.add(&[t1, neg(ctx, t2)])))
}

/// `J₂ = ∫ x²/√q dx = x·√q/(2c) − (3b/4c)·J₁ − (a/2c)·J₀`.
fn int_x2_inv_sqrt<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    j0: Option<Atom<'a>>,
) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let t1 = normalize(
        ctx,
        ctx.mul(&[qd.x, qd.sqrt_q, inv(ctx, ctx.mul(&[ctx.num(2), qd.c]))]),
    );
    let mut terms = vec![t1];
    if !is_zero(qd.b) {
        let j1 = int_x_inv_sqrt(ctx, qd, j0)?;
        terms.push(neg(
            ctx,
            normalize(
                ctx,
                ctx.mul(&[ctx.num(3), qd.b, inv(ctx, ctx.mul(&[ctx.num(4), qd.c])), j1]),
            ),
        ));
    }
    if !is_zero(qd.a) {
        terms.push(neg(
            ctx,
            normalize(
                ctx,
                ctx.mul(&[qd.a, inv(ctx, ctx.mul(&[ctx.num(2), qd.c])), j0?]),
            ),
        ));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// `∫ (p₀ + p₁x + p₂x²)/√q dx` via the `J_m` ladder (same lazy-evidence
/// policy as [`poly_times_sqrt`]).
fn poly_over_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>, p: [Atom<'a>; 3]) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let j0 = int_inv_sqrt(ctx, qd);
    let mut terms: Vec<Atom<'a>> = Vec::new();
    if !is_zero(p[0]) {
        terms.push(ctx.mul(&[p[0], j0?]));
    }
    if !is_zero(p[1]) {
        terms.push(ctx.mul(&[p[1], int_x_inv_sqrt(ctx, qd, j0)?]));
    }
    if !is_zero(p[2]) {
        terms.push(ctx.mul(&[p[2], int_x2_inv_sqrt(ctx, qd, j0)?]));
    }
    if terms.is_empty() {
        return None;
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

// =========================================================================
// S1.3: ∫ dx/q^(3/2)
// =========================================================================

/// `∫ p₀·q^(−3/2) dx = p₀·(4c·x + 2b)/((4ac − b²)·√q)`. A pure algebraic
/// identity — valid for both sign branches, no `J₀` evidence needed.
fn int_inv_q32<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>, p0: Atom<'a>) -> Option<Atom<'a>> {
    if is_zero(qd.c) {
        return None;
    }
    let d = delta(ctx, qd);
    if is_zero(d) {
        return None;
    }
    let four_ac_minus_b2 = neg(ctx, d);
    let num = normalize(
        ctx,
        ctx.add(&[
            ctx.mul(&[ctx.num(4), qd.c, qd.x]),
            ctx.mul(&[ctx.num(2), qd.b]),
        ]),
    );
    Some(normalize(
        ctx,
        ctx.mul(&[p0, num, inv(ctx, four_ac_minus_b2), inv(ctx, qd.sqrt_q)]),
    ))
}

// =========================================================================
// S2: reciprocal forms
// =========================================================================

/// `K₀ = ∫ dx/(x·√q)`, `a ≠ 0` (also valid for linear `q`, i.e. `c = 0`).
/// Rational `a > 0`: `(1/√a)·log(x / (2√a·√q + b·x + 2a))`. Rational
/// `a < 0` with `Δ > 0`: `(1/√(−a))·asin((b·x + 2a)/(x·√Δ))`. Symbolic
/// `a` declines (sign undecidable).
fn int_inv_x_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>) -> Option<Atom<'a>> {
    if is_zero(qd.a) {
        return None;
    }
    let (pa, _) = rat_of(qd.a)?;
    if pa > 0 {
        let sq_a = sqrt_emit(ctx, qd.a)?;
        let den = normalize(
            ctx,
            ctx.add(&[
                ctx.mul(&[ctx.num(2), sq_a, qd.sqrt_q]),
                ctx.mul(&[qd.b, qd.x]),
                ctx.mul(&[ctx.num(2), qd.a]),
            ]),
        );
        let arg = normalize(ctx, ctx.mul(&[qd.x, inv(ctx, den)]));
        return Some(normalize(
            ctx,
            ctx.mul(&[inv(ctx, sq_a), ctx.fun("log", &[arg])]),
        ));
    }
    let d = delta(ctx, qd);
    let (pd, _) = rat_of(d)?;
    if pd <= 0 {
        return None;
    }
    let sq_d = sqrt_emit(ctx, d)?;
    let sq_na = sqrt_emit(ctx, neg(ctx, qd.a))?;
    let num = normalize(
        ctx,
        ctx.add(&[ctx.mul(&[qd.b, qd.x]), ctx.mul(&[ctx.num(2), qd.a])]),
    );
    let arg = normalize(
        ctx,
        ctx.mul(&[num, inv(ctx, normalize(ctx, ctx.mul(&[qd.x, sq_d])))]),
    );
    Some(normalize(
        ctx,
        ctx.mul(&[inv(ctx, sq_na), ctx.fun("asin", &[arg])]),
    ))
}

/// `∫ dx/(x²·√q) = −√q/(a·x) − (b/2a)·K₀`. The `K₀` term is dropped when
/// `b = 0`.
fn int_inv_x2_sqrt<'a>(ctx: &'a AtomArena<'a>, qd: &Quad<'a>) -> Option<Atom<'a>> {
    if is_zero(qd.a) {
        return None;
    }
    let t1 = neg(
        ctx,
        normalize(
            ctx,
            ctx.mul(&[qd.sqrt_q, inv(ctx, normalize(ctx, ctx.mul(&[qd.a, qd.x])))]),
        ),
    );
    if is_zero(qd.b) {
        return Some(t1);
    }
    let k0 = int_inv_x_sqrt(ctx, qd)?;
    let t2 = normalize(
        ctx,
        ctx.mul(&[qd.b, inv(ctx, ctx.mul(&[ctx.num(2), qd.a])), k0]),
    );
    Some(normalize(ctx, ctx.add(&[t1, neg(ctx, t2)])))
}

// =========================================================================
// S4: √q over a linear denominator, u = d + e·x
// =========================================================================

/// Replace every `sqrt(arg)` node with `arg == q_arg` (structural equality
/// with the u-quadratic's canonical atom) by `sqrt_x`. Used to map
/// `√q̃(u)` back to `√q(x)` *before* substituting `u = d + e·x`, since
/// `q̃(d + e·x)` does not syntactically fold to `q` under `normalize`.
fn replace_sqrt_atom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    q_arg: Atom<'a>,
    sqrt_x: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args)
            if name.as_str() == "sqrt" && args.len() == 1 && args[0] == q_arg =>
        {
            sqrt_x
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_atom(ctx, *a, q_arg, sqrt_x))
                .collect();
            ctx.fun(name.as_str(), &rebuilt)
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_atom(ctx, *a, q_arg, sqrt_x))
                .collect();
            ctx.add(&rebuilt)
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_atom(ctx, *a, q_arg, sqrt_x))
                .collect();
            ctx.mul(&rebuilt)
        }
        AtomNode::Pow(b, e) => {
            let nb = replace_sqrt_atom(ctx, *b, q_arg, sqrt_x);
            let ne = replace_sqrt_atom(ctx, *e, q_arg, sqrt_x);
            ctx.pow(nb, ne)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// `∫ p₀·q^(k/2)/(d + e·x) dx` for `k = ±1` via `u = d + e·x`. With
/// `q̃(u) = q((u−d)/e) = ã + b̃·u + c̃·u²` (unscaled — no `|e|` factors):
/// `k = −1` gives `(p₀/e)·K₀(q̃)`; `k = +1` gives
/// `(p₀/e)·[ã·K₀(q̃) + b̃·J₀(q̃) + c̃·J₁(q̃)]` (zero coefficients skip their
/// sub-integral and its evidence requirements).
fn try_linear_den<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    qd: &Quad<'a>,
    p0: Atom<'a>,
    k: i64,
    lin: (Atom<'a>, Atom<'a>),
    var: Symbol,
) -> Option<Atom<'a>> {
    let (e, d) = lin;
    let t_sym = pick_subst_symbol(expr, var)?;
    let u = ctx.var(t_sym.as_str());
    let e2 = fold_rational(ctx, normalize(ctx, ctx.mul(&[e, e])));
    let d2 = fold_rational(ctx, normalize(ctx, ctx.mul(&[d, d])));
    // q̃ coefficients: q((u−d)/e) expanded in powers of u.
    let c_t = fold_rational(ctx, normalize(ctx, ctx.mul(&[qd.c, inv(ctx, e2)])));
    let b_t = fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.add(&[
                ctx.mul(&[qd.b, inv(ctx, e)]),
                neg(
                    ctx,
                    normalize(ctx, ctx.mul(&[ctx.num(2), qd.c, d, inv(ctx, e2)])),
                ),
            ]),
        ),
    );
    let a_t = fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.add(&[
                qd.a,
                neg(ctx, normalize(ctx, ctx.mul(&[qd.b, d, inv(ctx, e)]))),
                normalize(ctx, ctx.mul(&[qd.c, d2, inv(ctx, e2)])),
            ]),
        ),
    );
    let qt = Quad::build(ctx, a_t, b_t, c_t, u);
    let result_u = if k == -1 {
        let k0 = int_inv_x_sqrt(ctx, &qt)?;
        normalize(ctx, ctx.mul(&[p0, inv(ctx, e), k0]))
    } else {
        let j0 = int_inv_sqrt(ctx, &qt);
        let mut terms: Vec<Atom<'a>> = Vec::new();
        if !is_zero(a_t) {
            terms.push(ctx.mul(&[a_t, int_inv_x_sqrt(ctx, &qt)?]));
        }
        if !is_zero(b_t) {
            terms.push(ctx.mul(&[b_t, j0?]));
        }
        if !is_zero(c_t) {
            terms.push(ctx.mul(&[c_t, int_x_inv_sqrt(ctx, &qt, j0)?]));
        }
        if terms.is_empty() {
            return None;
        }
        normalize(
            ctx,
            ctx.mul(&[p0, inv(ctx, e), normalize(ctx, ctx.add(&terms))]),
        )
    };
    // Map √q̃(u) back to √q(x) before substituting u = d + e·x.
    let with_x_sqrt = replace_sqrt_atom(ctx, result_u, qt.q, qd.sqrt_q);
    let u_back = normalize(ctx, ctx.add(&[d, ctx.mul(&[e, ctx.var(var.as_str())])]));
    let back = replace_symbol(ctx, with_x_sqrt, t_sym, u_back);
    Some(normalize(ctx, back))
}

// =========================================================================
// Factor matching
// =========================================================================

/// The matched shape `N(x)·q^(k/2) / (x^m · (d + e·x))`: numerator
/// coefficients `p₀ + p₁x + p₂x²` (constant w.r.t. `var`), denominator
/// `x^m` (`m ∈ {0,1,2}`) xor one linear factor, and one radical factor
/// `q^(k/2)` with `k ∈ {−3,−1,1,3}`.
struct RadMatch<'a> {
    num: [Atom<'a>; 3],
    x_den: i64,
    lin: Option<(Atom<'a>, Atom<'a>)>,
    inner: Atom<'a>,
    k: i64,
}

/// Parse a numerator/denominator polynomial factor of degree ≤ 2 in `var`
/// as `[p₀, p₁, p₂]` (implicit coefficients made explicit).
fn poly_coeffs<'a>(ctx: &'a AtomArena<'a>, f: Atom<'a>, var: Symbol) -> Option<[Atom<'a>; 3]> {
    let (c2, c1, c0) = quadratic_coeffs(ctx, ctx.fun("sqrt", &[f]), var)?;
    Some([c0, c1, c2])
}

/// Multiply the numerator polynomial by `[g₀, g₁, g₂]`; `None` when the
/// product would exceed degree 2.
fn convolve<'a>(
    ctx: &'a AtomArena<'a>,
    num: [Atom<'a>; 3],
    g: [Atom<'a>; 3],
) -> Option<[Atom<'a>; 3]> {
    let mut out = [ctx.num(0), ctx.num(0), ctx.num(0)];
    for (i, ni) in num.iter().enumerate() {
        for (j, gj) in g.iter().enumerate() {
            let t = normalize(ctx, ctx.mul(&[*ni, *gj]));
            if is_zero(t) {
                continue;
            }
            if i + j > 2 {
                return None;
            }
            out[i + j] = normalize(ctx, ctx.add(&[out[i + j], t]));
        }
    }
    Some(out)
}

/// Decompose `expr` into the [`RadMatch`] shape. Any factor outside the
/// shape grammar (non-sqrt functions of `var`, higher-degree polynomials,
/// a second radical, fractional powers of anything but `q`, …) → `None`.
fn match_radical<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<RadMatch<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut num = [ctx.num(1), ctx.num(0), ctx.num(0)];
    let mut x_den = 0i64;
    let mut lin: Option<(Atom<'a>, Atom<'a>)> = None;
    let mut rad: Option<(Atom<'a>, i64)> = None;
    for f in factors {
        if is_constant(f, var) {
            for c in num.iter_mut() {
                *c = normalize(ctx, ctx.mul(&[*c, f]));
            }
            continue;
        }
        match f.node() {
            AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
                if rad.is_some() {
                    return None;
                }
                rad = Some((args[0], 1));
            }
            AtomNode::Pow(base, ex) => {
                // sqrt(q)^n with odd |n| ≤ 3.
                if let AtomNode::Fun(name, args) = base.node()
                    && name.as_str() == "sqrt"
                    && args.len() == 1
                {
                    if rad.is_some() {
                        return None;
                    }
                    match ex.node() {
                        AtomNode::Num(n) if matches!(n, -3 | -1 | 1 | 3) => {
                            rad = Some((args[0], *n));
                        }
                        _ => return None,
                    }
                    continue;
                }
                // x^±m powers.
                if matches!(base.node(), AtomNode::Var(v) if *v == var) {
                    match ex.node() {
                        AtomNode::Num(1) => {
                            num = convolve(ctx, num, [ctx.num(0), ctx.num(1), ctx.num(0)])?;
                        }
                        AtomNode::Num(2) => {
                            num = convolve(ctx, num, [ctx.num(0), ctx.num(0), ctx.num(1)])?;
                        }
                        AtomNode::Num(-1) | AtomNode::Num(-2) => {
                            let AtomNode::Num(n) = ex.node() else {
                                return None;
                            };
                            x_den += n.checked_neg()?;
                            if x_den > 2 || lin.is_some() {
                                return None;
                            }
                        }
                        _ => return None,
                    }
                    continue;
                }
                // q^(p/2) with odd |p| ≤ 3 (normalized fractional-exponent form).
                if let Some((p, 2)) = rat_of(*ex)
                    && matches!(p, -3 | -1 | 1 | 3)
                    && !is_constant(*base, var)
                {
                    if rad.is_some() {
                        return None;
                    }
                    rad = Some((*base, p));
                    continue;
                }
                // One linear denominator factor (d + e·x)^-1, d ≠ 0.
                if matches!(ex.node(), AtomNode::Num(-1)) {
                    let (e, d) = linear_form(ctx, *base, var)?;
                    if is_zero(e) || is_zero(d) || lin.is_some() || x_den > 0 {
                        return None;
                    }
                    lin = Some((e, d));
                    continue;
                }
                return None;
            }
            _ => {
                // Remaining factor must be a polynomial of degree ≤ 2.
                let g = poly_coeffs(ctx, f, var)?;
                num = convolve(ctx, num, g)?;
            }
        }
    }
    let (inner, k) = rad?;
    Some(RadMatch {
        num,
        x_den,
        lin,
        inner,
        k,
    })
}

// =========================================================================
// S1/S2/S4 dispatch
// =========================================================================

fn try_direct<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let m = match_radical(ctx, expr, var)?;
    let qd = Quad::parse(ctx, m.inner, var)?;
    let p = m.num;
    let num_const = is_zero(p[1]) && is_zero(p[2]);
    match (m.k, m.x_den, m.lin) {
        (-3, 0, None) if num_const => int_inv_q32(ctx, &qd, p[0]),
        // q^(3/2) = q·√q: fold the constant numerator into q's coefficients.
        (3, 0, None) if num_const => {
            let pn = [
                normalize(ctx, ctx.mul(&[p[0], qd.a])),
                normalize(ctx, ctx.mul(&[p[0], qd.b])),
                normalize(ctx, ctx.mul(&[p[0], qd.c])),
            ];
            poly_times_sqrt(ctx, &qd, pn)
        }
        (-1, 0, None) => poly_over_sqrt(ctx, &qd, p),
        (1, 0, None) => poly_times_sqrt(ctx, &qd, p),
        // (p₀ + p₁x)/(x·√q) = p₀·K₀ + p₁·J₀.
        (-1, 1, None) if is_zero(p[2]) => {
            let mut terms: Vec<Atom<'a>> = Vec::new();
            if !is_zero(p[0]) {
                terms.push(ctx.mul(&[p[0], int_inv_x_sqrt(ctx, &qd)?]));
            }
            if !is_zero(p[1]) {
                terms.push(ctx.mul(&[p[1], int_inv_sqrt(ctx, &qd)?]));
            }
            if terms.is_empty() {
                return None;
            }
            Some(normalize(ctx, ctx.add(&terms)))
        }
        // (p₀ + p₁x)/(x²·√q) = p₀·(−√q/(a·x) − (b/2a)·K₀) + p₁·K₀.
        (-1, 2, None) if is_zero(p[2]) => {
            let mut terms: Vec<Atom<'a>> = Vec::new();
            if !is_zero(p[0]) {
                terms.push(ctx.mul(&[p[0], int_inv_x2_sqrt(ctx, &qd)?]));
            }
            if !is_zero(p[1]) {
                terms.push(ctx.mul(&[p[1], int_inv_x_sqrt(ctx, &qd)?]));
            }
            if terms.is_empty() {
                return None;
            }
            Some(normalize(ctx, ctx.add(&terms)))
        }
        (-1, 0, Some((e, d))) if num_const => try_linear_den(ctx, expr, &qd, p[0], -1, (e, d), var),
        (1, 0, Some((e, d))) if num_const => try_linear_den(ctx, expr, &qd, p[0], 1, (e, d), var),
        _ => None,
    }
}

// =========================================================================
// S3: Euler III substitution (two distinct rational roots)
// =========================================================================

/// Replace `sqrt(q_t)` and half-integer powers of `q_t` by (powers of) the
/// rational expression `sqrt_t`; everything else is rebuilt unchanged.
/// Leftover non-rational nodes are caught by the rationality gate.
fn replace_sqrt_forms<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    q_t: Atom<'a>,
    sqrt_t: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            if normalize(ctx, args[0]) == q_t {
                sqrt_t
            } else {
                expr
            }
        }
        AtomNode::Pow(base, ex) => {
            // sqrt(q)^(±n) and q^(p/2) forms over the substituted quadratic.
            if let AtomNode::Fun(name, args) = base.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && normalize(ctx, args[0]) == q_t
                && let AtomNode::Num(n) = ex.node()
            {
                return int_pow(ctx, sqrt_t, *n);
            }
            if let Some((p, 2)) = rat_of(*ex)
                && normalize(ctx, *base) == q_t
            {
                return int_pow(ctx, sqrt_t, p);
            }
            let nb = replace_sqrt_forms(ctx, *base, q_t, sqrt_t);
            let ne = replace_sqrt_forms(ctx, *ex, q_t, sqrt_t);
            ctx.pow(nb, ne)
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_forms(ctx, *a, q_t, sqrt_t))
                .collect();
            ctx.add(&rebuilt)
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_forms(ctx, *a, q_t, sqrt_t))
                .collect();
            ctx.mul(&rebuilt)
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt_forms(ctx, *a, q_t, sqrt_t))
                .collect();
            ctx.fun(name.as_str(), &rebuilt)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Structural check: `expr` is a rational function of the substitution
/// variable — integer powers only, no function nodes, and no leftover
/// occurrence of the original integration variable.
fn is_rational_t(expr: Atom<'_>, orig_var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v != orig_var,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| is_rational_t(*a, orig_var))
        }
        AtomNode::Pow(b, e) => matches!(e.node(), AtomNode::Num(_)) && is_rational_t(*b, orig_var),
        AtomNode::Fun(_, _) => false,
    }
}

/// Euler III: `q = a + b·x + c·x²` with rational coefficients and
/// `Δ = b² − 4ac` a positive rational square, i.e. two distinct rational
/// roots `r₁, r₂`. With `t² = c·(x−r₁)/(x−r₂)`:
/// `x = (t²·r₂ − c·r₁)/(t² − c)`, `√q = t·c·(r₂−r₁)/(t² − c)`,
/// `dx/dt = 2t·c·(r₁−r₂)/(t² − c)²` — any rational-in-`(x, √q)` integrand
/// becomes rational in `t`.
fn try_euler3<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let (c2, c1, c0) = quadratic_coeffs(ctx, expr, var)?;
    let (p2, _) = rat_of(c2)?;
    if p2 == 0 {
        return None;
    }
    // Δ = c1² − 4·c0·c2 must be a positive rational square. (c1² via `mul`:
    // `normalize` folds numeric `Mul` but not numeric `Pow`.)
    let d = fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.add(&[ctx.mul(&[c1, c1]), neg(ctx, ctx.mul(&[ctx.num(4), c0, c2]))]),
        ),
    );
    let (pd, qd_) = rat_of(d)?;
    if pd <= 0 {
        return None;
    }
    let (sp, sq) = rational_sqrt(pd, qd_)?;
    let s = rat_atom(ctx, sp, sq);
    let den = normalize(ctx, ctx.mul(&[ctx.num(2), c2]));
    let r1 = fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.mul(&[
                normalize(ctx, ctx.add(&[neg(ctx, c1), neg(ctx, s)])),
                inv(ctx, den),
            ]),
        ),
    );
    let r2 = fold_rational(
        ctx,
        normalize(
            ctx,
            ctx.mul(&[normalize(ctx, ctx.add(&[neg(ctx, c1), s])), inv(ctx, den)]),
        ),
    );

    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let t2 = int_pow(ctx, t, 2);
    let den_t = normalize(ctx, ctx.add(&[t2, neg(ctx, c2)]));
    let x_t = normalize(
        ctx,
        ctx.mul(&[
            normalize(
                ctx,
                ctx.add(&[ctx.mul(&[t2, r2]), neg(ctx, ctx.mul(&[c2, r1]))]),
            ),
            inv(ctx, den_t),
        ]),
    );
    let sqrt_t = normalize(
        ctx,
        ctx.mul(&[
            t,
            c2,
            normalize(ctx, ctx.add(&[r2, neg(ctx, r1)])),
            inv(ctx, den_t),
        ]),
    );
    let dx_dt = normalize(
        ctx,
        ctx.mul(&[
            ctx.num(2),
            t,
            c2,
            normalize(ctx, ctx.add(&[r1, neg(ctx, r2)])),
            inv(ctx, int_pow(ctx, den_t, 2)),
        ]),
    );

    // Structural match target: q(x(t)), normalized.
    let x = ctx.var(var.as_str());
    let q_atom = normalize(
        ctx,
        ctx.add(&[c0, ctx.mul(&[c1, x]), ctx.mul(&[c2, int_pow(ctx, x, 2)])]),
    );
    let q_t = normalize(ctx, replace_symbol(ctx, q_atom, var, x_t));

    let subbed = replace_symbol(ctx, expr, var, x_t);
    let subbed = replace_sqrt_forms(ctx, subbed, q_t, sqrt_t);
    let integrand_t = normalize(ctx, ctx.mul(&[subbed, dx_dt]));
    if node_count(integrand_t) > MAX_NODES || !is_rational_t(integrand_t, var) {
        return None;
    }
    let result_t = integrate_raw(ctx, integrand_t, t_sym, 0, true, 0, 0);
    if contains_integral(result_t) {
        return None;
    }
    // Back-substitute t = √(c·(x−r₁)/(x−r₂)).
    let t_back = ctx.fun(
        "sqrt",
        &[normalize(
            ctx,
            ctx.mul(&[
                c2,
                normalize(ctx, ctx.add(&[x, neg(ctx, r1)])),
                inv(ctx, normalize(ctx, ctx.add(&[x, neg(ctx, r2)]))),
            ]),
        )],
    );
    let back = replace_symbol(ctx, result_t, t_sym, t_back);
    Some(normalize(ctx, back))
}

// =========================================================================
// Entry
// =========================================================================

/// Integrate `expr` via the quadratic-radical mechanisms (see module
/// docs). Returns `None` outside the supported shapes, on undecidable
/// branch evidence, on a residue from the S3 re-entry, or when a budget is
/// exceeded.
pub(crate) fn integrate_sqrt_quadratic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    if node_count(expr) > MAX_NODES {
        return None;
    }
    let r = try_direct(ctx, expr, var).or_else(|| try_euler3(ctx, expr, var))?;
    Some(fold_numeric_rationals(ctx, r))
}

/// Fold `Num(p) · Num(q)^-1` factor pairs into reduced `p/q` rationals
/// throughout the tree — normalize does not merge a Num with a numeric
/// reciprocal power, so emission points would otherwise print `2·x·2^-1`
/// for `x`.
fn fold_numeric_rationals<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Mul(args) => {
            let mut num: i128 = 1;
            let mut den: i128 = 1;
            let mut rest: Vec<Atom<'a>> = Vec::new();
            let mut folded = false;
            for a in args.iter() {
                let a = fold_numeric_rationals(ctx, *a);
                match a.node() {
                    AtomNode::Num(n) => {
                        num *= i128::from(*n);
                        folded = true;
                    }
                    AtomNode::Pow(b, e)
                        if matches!(b.node(), AtomNode::Num(_))
                            && matches!(e.node(), AtomNode::Num(-1)) =>
                    {
                        let AtomNode::Num(bb) = b.node() else {
                            unreachable!()
                        };
                        den *= i128::from(*bb);
                        folded = true;
                    }
                    _ => rest.push(a),
                }
            }
            if !folded {
                return expr;
            }
            let coeff = match (i64::try_from(num), i64::try_from(den)) {
                (Ok(p), Ok(q)) => rat_atom(ctx, p, q),
                _ => return expr,
            };
            if matches!(coeff.node(), AtomNode::Num(0)) {
                return coeff;
            }
            if !matches!(coeff.node(), AtomNode::Num(1)) {
                rest.insert(0, coeff);
            }
            match rest.len() {
                0 => coeff,
                1 => {
                    if matches!(coeff.node(), AtomNode::Num(1)) {
                        rest[0]
                    } else {
                        normalize(ctx, ctx.mul(&rest))
                    }
                }
                _ => normalize(ctx, ctx.mul(&rest)),
            }
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| fold_numeric_rationals(ctx, *a))
                .collect();
            let r = ctx.add(&rebuilt);
            if r == expr { expr } else { r }
        }
        AtomNode::Pow(b, e) => {
            let nb = fold_numeric_rationals(ctx, *b);
            let ne = fold_numeric_rationals(ctx, *e);
            if nb == *b && ne == *e {
                expr
            } else {
                ctx.pow(nb, ne)
            }
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| fold_numeric_rationals(ctx, *a))
                .collect();
            let r = ctx.fun(name.as_str(), &rebuilt);
            if r == expr { expr } else { r }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification.
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
                    "asin" => v.asin(),
                    "atan" => v.atan(),
                    _ => return None,
                })
            }
        }
    }

    /// Run the mechanism directly, require `Some` without residue, and
    /// check `diff(result) == integrand` numerically at the sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) -> Atom<'a> {
        let result = integrate_sqrt_quadratic(ctx, integrand, var).expect("mechanism declined");
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
        result
    }

    /// `a + b·x + c·x²` with integer coefficients.
    fn mk_quad<'a>(ctx: &'a AtomArena<'a>, a: i64, b: i64, c: i64) -> Atom<'a> {
        let x = ctx.var("x");
        ctx.add(&[
            ctx.num(a),
            ctx.mul(&[ctx.num(b), x]),
            ctx.mul(&[ctx.num(c), ctx.pow(x, ctx.num(2))]),
        ])
    }

    #[test]
    fn inv_sqrt_c_pos_delta_neg_both_shapes() {
        // ∫ dx/√(4 + 3x + 2x²): c = 2 > 0 (non-square), Δ = 9 − 32 < 0.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 4, 3, 2);
        let f1 = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        assert_antiderivative_num(&ctx, f1, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
        // Pow(q, -1/2) parse shape.
        let f2 = ctx.pow(q, rat_atom(&ctx, -1, 2));
        assert_antiderivative_num(&ctx, f2, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn inv_sqrt_c_pos_delta_pos() {
        // ∫ dx/√(1 + 5x + 2x²): c > 0, Δ = 25 − 8 = 17 > 0 (log form covers it).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 1, 5, 2);
        let f = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.2, 1.0, 2.0]);
    }

    #[test]
    fn inv_sqrt_c_neg_asin() {
        // ∫ dx/√(3 + x − 2x²): c = −2 < 0, Δ = 25 > 0 → asin branch.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 3, 1, -2);
        let f = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[-0.5, 0.0, 1.0]);
    }

    #[test]
    fn sqrt_q_c_pos() {
        // ∫ √(4 + 3x + 2x²) dx.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 4, 3, 2);
        let f = ctx.fun("sqrt", &[q]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn sqrt_q_c_neg() {
        // ∫ √(3 + x − 2x²) dx (asin branch inside I₀).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 3, 1, -2);
        let f = ctx.fun("sqrt", &[q]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[-0.5, 0.0, 1.0]);
    }

    #[test]
    fn inv_q32_both_shapes() {
        // ∫ dx/(4 + 3x + 2x²)^(3/2) = (8x + 6)/(23·√q).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 4, 3, 2);
        let f1 = ctx.pow(q, rat_atom(&ctx, -3, 2));
        assert_antiderivative_num(&ctx, f1, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
        let f2 = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-3));
        assert_antiderivative_num(&ctx, f2, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn linear_over_sqrt() {
        // ∫ (3 + 5x)/√(4 + 3x + 2x²) dx — J ladder with p₁ ≠ 0.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 4, 3, 2);
        let num = ctx.add(&[ctx.num(3), ctx.mul(&[ctx.num(5), x])]);
        let f = ctx.mul(&[num, ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1))]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn quad_times_sqrt() {
        // ∫ (1 + x²)·√(4 + 3x + 2x²) dx — I ladder through I₂.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 4, 3, 2);
        let num = ctx.add(&[ctx.num(1), ctx.pow(x, ctx.num(2))]);
        let f = ctx.mul(&[num, ctx.fun("sqrt", &[q])]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn q32_as_q_times_sqrt() {
        // ∫ (4 + 3x + 2x²)^(3/2) dx — k = 3 folds q into the I ladder.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q = mk_quad(&ctx, 4, 3, 2);
        let f = ctx.pow(q, rat_atom(&ctx, 3, 2));
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn x_over_sqrt_b_zero_c_neg() {
        // ∫ x/√(4 − x²) dx = −√(4 − x²): b = 0 skips the J₀ evidence path.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 4, 0, -1);
        let f = ctx.mul(&[x, ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1))]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn inv_x_sqrt_a_pos() {
        // ∫ dx/(x·√(4 + 3x + 2x²)): a = 4 > 0 → log form.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 4, 3, 2);
        let f = ctx.mul(&[
            ctx.pow(x, ctx.num(-1)),
            ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1)),
        ]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn inv_x_sqrt_a_neg() {
        // ∫ dx/(x·√(−1 + 2x + x²)): a = −1 < 0, Δ = 8 > 0 → asin form.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, -1, 2, 1);
        let f = ctx.mul(&[
            ctx.pow(x, ctx.num(-1)),
            ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1)),
        ]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.6, 1.0, 2.0]);
    }

    #[test]
    fn inv_x2_sqrt() {
        // ∫ dx/(x²·√(4 + 3x + 2x²)) — S2.2 with K₀ sub-term.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 4, 3, 2);
        let f = ctx.mul(&[
            ctx.pow(x, ctx.num(-2)),
            ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1)),
        ]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn sqrt_over_linear() {
        // ∫ √(1 + x²)/(2 + 3x) dx — S4.1 (ã̃ = 13/9 > 0 → K₀ log form).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 1, 0, 1);
        let lin = ctx.add(&[ctx.num(2), ctx.mul(&[ctx.num(3), x])]);
        let f = ctx.mul(&[ctx.fun("sqrt", &[q]), ctx.pow(lin, ctx.num(-1))]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn inv_linear_sqrt() {
        // ∫ dx/((2 + 3x)·√(1 + x²)) — S4.2.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, 1, 0, 1);
        let lin = ctx.add(&[ctx.num(2), ctx.mul(&[ctx.num(3), x])]);
        let f = ctx.mul(&[
            ctx.pow(lin, ctx.num(-1)),
            ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1)),
        ]);
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn euler3_linear_denom() {
        // ∫ dx/((x + 1)·√(x² − 1)): S4.2 declines (ã̃ = 0), Euler III
        // rationalizes to ∫ −t⁻² dt = 1/t = √((x−1)/(x+1)).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = mk_quad(&ctx, -1, 0, 1);
        let lin = ctx.add(&[ctx.num(1), x]);
        let f = ctx.mul(&[
            ctx.pow(lin, ctx.num(-1)),
            ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1)),
        ]);
        let result = assert_antiderivative_num(&ctx, f, Symbol::new("x"), &[], &[1.5, 2.0, 3.0]);
        // Idempotency: the S3 output shape (sqrt of a quotient) must not
        // re-match this module.
        assert!(integrate_sqrt_quadratic(&ctx, result, Symbol::new("x")).is_none());
    }

    #[test]
    fn symbolic_coeffs_log_form() {
        // ∫ dx/√(a + b·x + c·x²) with symbolic coefficients → log form.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let q = ctx.add(&[a, ctx.mul(&[b, x]), ctx.mul(&[c, ctx.pow(x, ctx.num(2))])]);
        let f = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        let env = [
            (Symbol::new("a"), 4.0),
            (Symbol::new("b"), 3.0),
            (Symbol::new("c"), 2.0),
        ];
        assert_antiderivative_num(&ctx, f, Symbol::new("x"), &env, &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn declines_var_dependent_coeff() {
        // 1/√(1 + x·sin(x) + x²): coefficient of x is not constant → None.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = ctx.add(&[
            ctx.num(1),
            ctx.mul(&[x, ctx.fun("sin", &[x])]),
            ctx.pow(x, ctx.num(2)),
        ]);
        let f = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        assert!(integrate_sqrt_quadratic(&ctx, f, Symbol::new("x")).is_none());
    }

    #[test]
    fn declines_degree3_under_radical() {
        // 1/√(1 + x³): degree > 2 → None.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let q = ctx.add(&[ctx.num(1), ctx.pow(x, ctx.num(3))]);
        let f = ctx.pow(ctx.fun("sqrt", &[q]), ctx.num(-1));
        assert!(integrate_sqrt_quadratic(&ctx, f, Symbol::new("x")).is_none());
    }

    #[test]
    fn declines_two_distinct_radicals() {
        // √(x² + 3x + 1)·√(x² + 1): two different quadratics → None
        // (direct matcher rejects the second radical; Euler III finds no
        // perfect-square Δ on either).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let q1 = mk_quad(&ctx, 1, 3, 1);
        let q2 = mk_quad(&ctx, 1, 0, 1);
        let f = ctx.mul(&[ctx.fun("sqrt", &[q1]), ctx.fun("sqrt", &[q2])]);
        assert!(integrate_sqrt_quadratic(&ctx, f, Symbol::new("x")).is_none());
    }
}
