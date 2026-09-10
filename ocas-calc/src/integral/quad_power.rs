//! Quadratic/linear denominator power reductions (0.27.1 Phase 2).
//!
//! Closed-form recurrences for `P(x)/(a + b·x + c·x²)^n` and
//! `P(x)/(d + e·x)^n` with `P` a polynomial of degree ≤ 4 whose
//! coefficients are constant w.r.t. `x` (symbolic allowed). These corpus
//! families (`(d+e·x)^4/(a+b·x+c·x²)^4`, …) stall the general symbolic
//! rational backend: its multivariate coefficient-field gcd explodes on
//! 5+ generator shapes. The recurrences below are pure atom arithmetic —
//! no gcd, no factorization.
//!
//! - **L1** `P(x)/(d + e·x)^n`: shift `y = d + e·x`, expand into
//!   `Σ c_k·y^(k−n)` and integrate each power directly (`k−n = −1` → log).
//! - **L2** `P(x)/q^n` (q quadratic, `c ≠ 0`): `Σ p_m·J(m,n)` with
//!   `J(0,n) = I_n`, `J(1,n) = (q^(1−n)/(1−n) − b·I_n)/(2c)` and
//!
//!   ```text
//!   I_n = (2c·x+b)/((n−1)·D·q^(n−1)) + 2c(2n−3)/((n−1)·D)·I_{n−1},
//!       D = 4ac − b² ≠ 0
//!   J(m,n) = [ x^(m−1)·q^(1−n) − (m−1)a·J(m−2,n) − b(m−n)·J(m−1,n) ]
//!            / [c(m−2n+1)],     c(m−2n+1) ≠ 0
//!   ```
//!
//!   (Both derived from `d/dx[x^k·q^(1−n)]` expansions; verified
//!   numerically in tests.) The `I_1 = ∫dx/q` base re-enters the chain
//!   (the symbolic rational backend owns quadratic denominators).
//!
//! The singular cases (`D = 0`, `c(m−2n+1) = 0` when needed) decline
//! honestly. All recursion is on decreasing `n`/`m` with tight budgets, so
//! termination is structural.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::{
    contains_integral, int_pow, integrate_raw, inv, is_constant, is_fallback, node_count, rat_atom,
};

/// Cap on the numerator polynomial degree for L2 (J-recursion depth).
const MAX_NUM_DEG: usize = 4;
/// Cap on the numerator polynomial degree for L1/L3 (direct power
/// integration is cheap at higher degrees).
const MAX_NUM_DEG_WIDE: usize = 8;
/// Cap on the denominator power.
const MAX_DENOM_POW: i64 = 20;
/// Node budget for the input integrand.
const MAX_NODES: usize = 200;
/// Cap on the combined linear-factor powers in L3 partial fractions.
const MAX_L3_POW: i64 = 8;

/// Integrate `expr` via the denominator-power recurrences (module docs).
pub(crate) fn integrate_quad_power<'a>(
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
    try_linear_denom(ctx, expr, var)
        .or_else(|| try_quad_denom(ctx, expr, var))
        .or_else(|| try_two_linear(ctx, expr, var))
}

// =========================================================================
// Shared matching
// =========================================================================

/// A matched `P(x) / base^n` shape.
struct DenomPowMatch<'a> {
    /// Numerator coefficients in ascending x-degree (sparse; constants
    /// w.r.t. `var`).
    numer: Vec<(usize, Atom<'a>)>,
    /// The denominator base atom (linear or quadratic in `var`).
    base: Atom<'a>,
    /// The denominator power `n ≥ 2`.
    n: i64,
    /// Constant factors split off the whole product.
    consts: Vec<Atom<'a>>,
}

/// Decompose `expr` as `consts · P(x) / base^n`. The numerator `P` may be
/// an explicit `Add`, a positive binomial/trinomial power (expanded), or a
/// monomial. Any non-polynomial-in-x factor → None.
fn match_denom_pow<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<DenomPowMatch<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut numer: Vec<(usize, Atom<'a>)> = Vec::new();
    // Multiplicative numerator monomial degree (x^k factors shift every
    // additive numerator term's degree).
    let mut mono_deg: usize = 0;
    let mut denom: Option<(Atom<'a>, i64)> = None;
    for f in factors {
        if is_constant(f, var) {
            consts.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Pow(b, e) => {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                let n = *n;
                if n == 0 || n.abs() > MAX_DENOM_POW {
                    return None;
                }
                // `x^k` is a multiplicative numerator monomial.
                if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                    if n < 0 {
                        return None;
                    }
                    mono_deg = mono_deg.checked_add(n as usize)?;
                    continue;
                }
                if n < 0 {
                    // Denominator power; the base must be linear/quadratic.
                    if denom.is_some() {
                        return None;
                    }
                    if !matches!(b.node(), AtomNode::Add(_)) {
                        return None;
                    }
                    denom = Some((*b, -n));
                } else {
                    // Numerator binomial power: expand into P terms (at
                    // most one additive numerator source).
                    if !numer.is_empty() {
                        return None;
                    }
                    numer = expand_pow_add(ctx, *b, n, var)?;
                }
            }
            AtomNode::Add(add_args) => {
                if !numer.is_empty() {
                    return None;
                }
                let mut terms = Vec::new();
                for t in add_args.iter() {
                    terms.push(term_as_monomial(ctx, *t, var)?);
                }
                numer = terms;
            }
            _ => {
                // A bare `x` monomial factor shifts degrees; anything else
                // must be a single-monomial numerator.
                if matches!(f.node(), AtomNode::Var(v) if *v == var) {
                    mono_deg = mono_deg.checked_add(1)?;
                    continue;
                }
                let (d, c) = term_as_monomial(ctx, f, var)?;
                if d != 0 {
                    mono_deg = mono_deg.checked_add(d)?;
                    if !matches!(normalize(ctx, c).node(), AtomNode::Num(1)) {
                        consts.push(c);
                    }
                    continue;
                }
                if !numer.is_empty() {
                    return None;
                }
                numer = vec![(0, c)];
            }
        }
    }
    let (base, n) = denom?;
    if n < 2 {
        // n = 1 belongs to the plain rational backends.
        return None;
    }
    // Pure denominator power: numerator is the constant 1.
    if numer.is_empty() {
        numer.push((0, ctx.num(1)));
    }
    // The multiplicative monomial factors shift every numerator degree.
    if mono_deg > 0 {
        for (d, _) in numer.iter_mut() {
            *d = d.checked_add(mono_deg)?;
        }
    }
    // Degree budget (wide; L2 narrows to MAX_NUM_DEG itself).
    if numer.iter().any(|(d, _)| *d > MAX_NUM_DEG_WIDE) {
        return None;
    }
    Some(DenomPowMatch {
        numer,
        base,
        n,
        consts,
    })
}

/// Match `term` as `c·x^d` with `c` constant: returns `(d, c)`.
fn term_as_monomial<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
) -> Option<(usize, Atom<'a>)> {
    if is_constant(term, var) {
        return Some((0, term));
    }
    let sub: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut deg: Option<i64> = None;
    for f in sub {
        match f.node() {
            AtomNode::Var(v) if *v == var => {
                if deg.is_some() {
                    return None;
                }
                deg = Some(1);
            }
            AtomNode::Pow(b, e)
                if matches!(b.node(), AtomNode::Var(v) if *v == var)
                    && matches!(e.node(), AtomNode::Num(_)) =>
            {
                if deg.is_some() {
                    return None;
                }
                let AtomNode::Num(k) = e.node() else {
                    return None;
                };
                if *k < 0 {
                    return None;
                }
                deg = Some(*k);
            }
            _ => {
                if is_constant(f, var) {
                    coeff.push(f);
                } else {
                    return None;
                }
            }
        }
    }
    let deg = deg?;
    if deg > MAX_NUM_DEG_WIDE as i64 {
        return None;
    }
    let c = if coeff.is_empty() {
        ctx.num(1)
    } else {
        normalize(ctx, ctx.mul(&coeff))
    };
    Some((deg as usize, c))
}

/// Expand `(add)^n` (small positive integer n) into monomial terms by
/// distributing the n-fold product of the term list; declines when any
/// term is not a monomial in `var` or the result exceeds the degree cap.
fn expand_pow_add<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    n: i64,
    var: Symbol,
) -> Option<Vec<(usize, Atom<'a>)>> {
    let AtomNode::Add(args) = base.node() else {
        return None;
    };
    if !(2..=4).contains(&args.len()) {
        return None;
    }
    let mut mono: Vec<(usize, Atom<'a>)> = Vec::with_capacity(args.len());
    for t in args.iter() {
        mono.push(term_as_monomial(ctx, *t, var)?);
    }
    expand_terms(ctx, &mono, n)
}

/// Distribute `Π_i (sum of monomials_i)` over the sums and collect
/// monomial terms: the direct n-fold-free expansion of `(Σ t)^n` is
/// computed by expanding the product of `n` copies of the term list.
fn expand_terms<'a>(
    ctx: &'a AtomArena<'a>,
    mono: &[(usize, Atom<'a>)],
    n: i64,
) -> Option<Vec<(usize, Atom<'a>)>> {
    // Start from the empty product: {(0, 1)}.
    let mut cur: Vec<(usize, Atom<'a>)> = vec![(0, ctx.num(1))];
    for _ in 0..n {
        let mut next: Vec<(usize, Atom<'a>)> = Vec::new();
        for (d1, c1) in &cur {
            for (d2, c2) in mono {
                let d = d1 + d2;
                if d > MAX_NUM_DEG_WIDE {
                    return None;
                }
                next.push((d, normalize(ctx, ctx.mul(&[*c1, *c2]))));
            }
        }
        if next.len() > 400 {
            return None;
        }
        // Merge equal degrees each round to keep the list short.
        let mut merged: std::collections::BTreeMap<usize, Vec<Atom<'a>>> =
            std::collections::BTreeMap::new();
        for (d, c) in next {
            merged.entry(d).or_default().push(c);
        }
        cur = merged
            .into_iter()
            .map(|(d, cs)| (d, normalize(ctx, ctx.add(&cs))))
            .collect();
    }
    Some(cur)
}

fn binom(n: i64, k: i64) -> i64 {
    let mut r = 1i64;
    for i in 0..k {
        r = r.saturating_mul(n - i) / (i + 1);
    }
    r
}

fn is_zero<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> bool {
    if matches!(normalize(ctx, a).node(), AtomNode::Num(0)) {
        return true;
    }
    // Mixed symbolic cancellations (e.g. `b·d³·(b·d³)⁻¹ − 1`) need the
    // collect_terms fold.
    super::inverse_trig::folds_to_zero(ctx, a)
}

/// `q(x)` coefficients `(a, b, c)` for a quadratic `Add` (missing terms
/// are zero), or `(a, b)` for a linear one. Any non-constant,
/// non-monomial term → None.
fn poly_coeffs_small<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    var: Symbol,
) -> Option<Vec<Atom<'a>>> {
    let AtomNode::Add(args) = base.node() else {
        return None;
    };
    let mut coeffs: Vec<Atom<'a>> = vec![ctx.num(0), ctx.num(0), ctx.num(0)];
    for t in args.iter() {
        let (d, c) = term_as_monomial(ctx, *t, var)?;
        if d > 2 {
            return None;
        }
        coeffs[d] = c;
    }
    Some(coeffs)
}

// =========================================================================
// L1: P(x)/(d + e·x)^n
// =========================================================================

fn try_linear_denom<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let m = match_denom_pow(ctx, expr, var)?;
    let coeffs = poly_coeffs_small(ctx, m.base, var)?;
    let (d0, e0) = (coeffs[0], coeffs[1]);
    if !is_zero(ctx, coeffs[2]) || is_zero(ctx, e0) {
        return None;
    }
    // y = d + e·x ⇒ x = (y − d)/e, dx = dy/e. P(x) → Σ_k c_k·y^k where
    // c_k = Σ_i p_i·C(i,k)·(−d)^(i−k)·e^(−i−1) (the extra e^(−1) is dx/dy).
    let x = ctx.var(var.as_str());
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for (i, p) in &m.numer {
        let i = *i as i64;
        for k in 0..=i {
            // (−d)^(i−k): sign when (i−k) odd.
            let mut cf = vec![
                *p,
                ctx.num(binom(i, k)),
                int_pow(ctx, d0, i - k),
                ctx.pow(e0, ctx.num(-(i + 1))),
            ];
            if (i - k) % 2 == 1 {
                cf.push(ctx.num(-1));
            }
            let coeff = normalize(ctx, ctx.mul(&cf));
            let ypow = k - m.n;
            // y^(k−n) integrates to y^(k−n+1)/(k−n+1) or log(y).
            let integrated = if ypow == -1 {
                ctx.fun("log", &[m.base])
            } else {
                normalize(
                    ctx,
                    ctx.mul(&[int_pow(ctx, m.base, ypow + 1), rat_atom(ctx, 1, ypow + 1)]),
                )
            };
            terms.push(normalize(ctx, ctx.mul(&[coeff, integrated])));
        }
    }
    let mut all = m.consts.clone();
    all.push(ctx.add(&terms));
    let _ = x;
    Some(normalize(ctx, ctx.mul(&all)))
}

// =========================================================================
// L2: P(x)/q^n, q quadratic
// =========================================================================

/// `(a, b, c)` of the quadratic, with `c` nonzero.
struct Quad<'a> {
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    /// The quadratic atom `q` itself.
    q: Atom<'a>,
    /// `D = 4ac − b²` (normalized).
    disc: Atom<'a>,
}

impl<'a> Quad<'a> {
    fn parse(ctx: &'a AtomArena<'a>, base: Atom<'a>, var: Symbol) -> Option<Quad<'a>> {
        let coeffs = poly_coeffs_small(ctx, base, var)?;
        let c = coeffs[2];
        if is_zero(ctx, c) {
            return None;
        }
        let (a, b) = (coeffs[0], coeffs[1]);
        let disc = normalize(
            ctx,
            ctx.add(&[ctx.mul(&[ctx.num(4), a, c]), ctx.mul(&[ctx.num(-1), b, b])]),
        );
        if is_zero(ctx, disc) {
            return None;
        }
        Some(Quad {
            a,
            b,
            c,
            q: base,
            disc,
        })
    }
}

/// `I_n = ∫ dx/q^n` by the Q1 recurrence (module docs). `n ≥ 2` here; the
/// `n = 1` base is handled by the caller via the chain.
fn int_inv_quad_pow<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    n: i64,
    var: Symbol,
) -> Option<Atom<'a>> {
    debug_assert!(n >= 2);
    let x = ctx.var(var.as_str());
    // I_1 via the chain (symbolic rational owns the quadratic denominator).
    let q_inv = ctx.pow(qd.q, ctx.num(-1));
    let i1 = integrate_raw(ctx, q_inv, var, 0, true, 0, 0);
    if contains_integral(i1) {
        return None;
    }
    let mut result = i1;
    let mut k = 2;
    while k <= n {
        // I_k = (2c·x+b)/((k−1)D·q^(k−1)) + 2c(2k−3)/((k−1)D)·I_{k−1}.
        let kd = normalize(ctx, ctx.mul(&[ctx.num(k - 1), qd.disc]));
        let boundary = normalize(
            ctx,
            ctx.mul(&[
                ctx.add(&[ctx.mul(&[ctx.num(2), qd.c, x]), qd.b]),
                inv(ctx, kd),
                ctx.pow(qd.q, ctx.num(-(k - 1))),
            ]),
        );
        let coeff = normalize(
            ctx,
            ctx.mul(&[ctx.num(2), qd.c, ctx.num(2 * k - 3), inv(ctx, kd)]),
        );
        result = normalize(ctx, ctx.add(&[boundary, ctx.mul(&[coeff, result])]));
        k += 1;
    }
    Some(result)
}

/// `J(m, n)` by the Q2 recurrence (module docs); `memo[(m, n)]` caches.
fn int_x_pow_over_quad<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    m: usize,
    n: i64,
    var: Symbol,
    memo: &mut std::collections::HashMap<(usize, i64), Option<Atom<'a>>>,
) -> Option<Atom<'a>> {
    if let Some(cached) = memo.get(&(m, n)) {
        return *cached;
    }
    let result = int_x_pow_over_quad_inner(ctx, qd, m, n, var, memo);
    memo.insert((m, n), result);
    result
}

fn int_x_pow_over_quad_inner<'a>(
    ctx: &'a AtomArena<'a>,
    qd: &Quad<'a>,
    m: usize,
    n: i64,
    var: Symbol,
    memo: &mut std::collections::HashMap<(usize, i64), Option<Atom<'a>>>,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    if m == 0 {
        if n == 1 {
            let q_inv = ctx.pow(qd.q, ctx.num(-1));
            let r = integrate_raw(ctx, q_inv, var, 0, true, 0, 0);
            if contains_integral(r) {
                return None;
            }
            return Some(r);
        }
        return int_inv_quad_pow(ctx, qd, n, var);
    }
    if n == 0 {
        // Plain monomial: x^(m+1)/(m+1).
        let p1 = m as i64 + 1;
        return Some(normalize(
            ctx,
            ctx.mul(&[int_pow(ctx, x, p1), rat_atom(ctx, 1, p1)]),
        ));
    }
    if m == 1 {
        if n == 1 {
            // J(1,1) = (log q − b·I_1)/(2c).
            let i1 = int_x_pow_over_quad(ctx, qd, 0, 1, var, memo)?;
            let logq = ctx.fun("log", &[qd.q]);
            return Some(normalize(
                ctx,
                ctx.mul(&[
                    inv(ctx, ctx.mul(&[ctx.num(2), qd.c])),
                    ctx.add(&[logq, ctx.mul(&[ctx.num(-1), qd.b, i1])]),
                ]),
            ));
        }
        // J(1,n) = (q^(1−n)/(1−n) − b·I_n)/(2c).
        let i_n = int_x_pow_over_quad(ctx, qd, 0, n, var, memo)?;
        let qterm = normalize(
            ctx,
            ctx.mul(&[ctx.pow(qd.q, ctx.num(1 - n)), rat_atom(ctx, 1, 1 - n)]),
        );
        return Some(normalize(
            ctx,
            ctx.mul(&[
                inv(ctx, ctx.mul(&[ctx.num(2), qd.c])),
                ctx.add(&[qterm, ctx.mul(&[ctx.num(-1), qd.b, i_n])]),
            ]),
        ));
    }
    // J(m,n) = [x^(m−1)q^(1−n) − (m−1)a·J(m−2,n) − b(m−n)·J(m−1,n)]
    //          / [c(m−2n+1)].
    let denom_factor = (m as i64).checked_sub(2 * n)?.checked_add(1)?;
    if denom_factor == 0 {
        return None; // singular case
    }
    let jm2 = int_x_pow_over_quad(ctx, qd, m - 2, n, var, memo)?;
    let jm1 = int_x_pow_over_quad(ctx, qd, m - 1, n, var, memo)?;
    let boundary = normalize(
        ctx,
        ctx.mul(&[int_pow(ctx, x, m as i64 - 1), ctx.pow(qd.q, ctx.num(1 - n))]),
    );
    let num = normalize(
        ctx,
        ctx.add(&[
            boundary,
            ctx.mul(&[ctx.num(-(m as i64 - 1)), qd.a, jm2]),
            ctx.mul(&[ctx.num(-1), qd.b, ctx.num(m as i64 - n), jm1]),
        ]),
    );
    let denom = normalize(ctx, ctx.mul(&[qd.c, ctx.num(denom_factor)]));
    Some(normalize(ctx, ctx.mul(&[num, inv(ctx, denom)])))
}

fn try_quad_denom<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let m = match_denom_pow(ctx, expr, var)?;
    // L2's recursion budget is tighter than the wide match cap.
    if m.numer.iter().any(|(d, _)| *d > MAX_NUM_DEG) {
        return None;
    }
    let qd = Quad::parse(ctx, m.base, var)?;
    let mut memo = std::collections::HashMap::new();
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for (deg, coeff) in &m.numer {
        let j = int_x_pow_over_quad(ctx, &qd, *deg, m.n, var, &mut memo)?;
        terms.push(normalize(ctx, ctx.mul(&[*coeff, j])));
    }
    let mut all = m.consts.clone();
    all.push(ctx.add(&terms));
    let sum = normalize(ctx, ctx.mul(&all));
    if sum == expr {
        return None;
    }
    // Sanity: the closed form must be residue-free by construction (the
    // only chain re-entry is the I_1 base, checked inside).
    if is_fallback(&sum) || contains_integral(sum) {
        return None;
    }
    Some(sum)
}

// =========================================================================
// L3: P(x) / ((d1 + e1·x)^m · (d2 + e2·x)^n) via partial fractions
// =========================================================================

/// One linear denominator power `1/(d + e·x)^k`.
struct LinDenom<'a> {
    d: Atom<'a>,
    e: Atom<'a>,
    k: i64,
    /// The base atom `d + e·x`.
    base: Atom<'a>,
}

/// Coefficient of `1/(d + e·x)^k` in the partial-fraction decomposition of
/// `P(x)·(d + e·x)^(−m)·rest(x)`, evaluated by the standard rule
/// `(d/dx)^(m−k)[P·rest](root) / ((m−k)!·e^(m−k))` at the root `x = −d/e`.
/// The `e^(m−k)` factor accounts for the slope each `d/dx` pulls out of
/// `l = d + e·x` (0.27.1 fix).
fn pf_coeff<'a>(
    ctx: &'a AtomArena<'a>,
    body: Atom<'a>,
    derivs: i64,
    root: Atom<'a>,
    slope: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let mut d = body;
    for _ in 0..derivs {
        d = normalize(ctx, crate::diff(ctx, d, var));
    }
    let at = normalize(ctx, super::replace_symbol(ctx, d, var, root));
    let mut fact = 1i64;
    for i in 2..=derivs {
        fact = fact.checked_mul(i)?;
    }
    let scale = normalize(
        ctx,
        ctx.mul(&[rat_atom(ctx, 1, fact), int_pow(ctx, slope, -derivs)]),
    );
    Some(normalize(ctx, ctx.mul(&[at, scale])))
}

/// `P(x) / (l1^m · l2^n)` with two distinct linear denominators.
fn try_two_linear<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let factors0: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    // Distribute `(f·g)^k → f^k·g^k` for integer k < 0: the parser keeps a
    // product denominator as one `Pow(Mul, -k)` node.
    let mut factors: Vec<Atom<'a>> = Vec::new();
    for f in factors0 {
        if let AtomNode::Pow(b, e) = f.node()
            && let (AtomNode::Mul(margs), AtomNode::Num(k)) = (b.node(), e.node())
            && *k < 0
        {
            for a in margs.iter() {
                // normalize folds nested powers like `(2+x)^2` to `^-2`.
                factors.push(normalize(ctx, ctx.pow(*a, *e)));
            }
            continue;
        }
        factors.push(f);
    }
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut dens: Vec<LinDenom<'a>> = Vec::new();
    let mut numer: Vec<(usize, Atom<'a>)> = Vec::new();
    let mut mono_deg: usize = 0;
    for f in factors {
        if is_constant(f, var) {
            consts.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Pow(b, e) => {
                let AtomNode::Num(k) = e.node() else {
                    return None;
                };
                let k = *k;
                if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                    if k < 0 {
                        return None;
                    }
                    mono_deg = mono_deg.checked_add(k as usize)?;
                    continue;
                }
                if k >= 0 || !matches!(b.node(), AtomNode::Add(_)) {
                    return None;
                }
                let coeffs = poly_coeffs_small(ctx, *b, var)?;
                if !is_zero(ctx, coeffs[2]) || is_zero(ctx, coeffs[1]) {
                    return None;
                }
                dens.push(LinDenom {
                    d: coeffs[0],
                    e: coeffs[1],
                    k: -k,
                    base: *b,
                });
            }
            AtomNode::Add(add_args) => {
                if !numer.is_empty() {
                    return None;
                }
                for t in add_args.iter() {
                    numer.push(term_as_monomial(ctx, *t, var)?);
                }
            }
            _ => {
                if matches!(f.node(), AtomNode::Var(v) if *v == var) {
                    mono_deg = mono_deg.checked_add(1)?;
                    continue;
                }
                let (d, c) = term_as_monomial(ctx, f, var)?;
                if d != 0 {
                    mono_deg = mono_deg.checked_add(d)?;
                    consts.push(c);
                    continue;
                }
                if !numer.is_empty() {
                    return None;
                }
                numer = vec![(0, c)];
            }
        }
    }
    if dens.len() != 2 {
        return None;
    }
    let (l1, l2) = (&dens[0], &dens[1]);
    if l1.k + l2.k > MAX_L3_POW {
        return None;
    }
    // Distinct, non-proportional linear factors required.
    let cross = normalize(
        ctx,
        ctx.add(&[ctx.mul(&[l1.d, l2.e]), ctx.mul(&[ctx.num(-1), l2.d, l1.e])]),
    );
    if is_zero(ctx, cross) {
        return None;
    }
    // Numerator must be proper: deg P < m + n. The equal-degree edge takes
    // the constant quotient out first (P = C·l1^m·l2^n + R, deg R < m + n).
    let x = ctx.var(var.as_str());
    if numer.is_empty() {
        numer.push((0, ctx.num(1)));
    }
    if mono_deg > 0 {
        for (d, _) in numer.iter_mut() {
            *d = d.checked_add(mono_deg)?;
        }
    }
    let total = (l1.k + l2.k) as usize;
    let mut quotient: Option<Atom<'a>> = None;
    let mut deg_p = numer.iter().map(|(d, _)| *d).max()?;
    if deg_p == total {
        // C = lc(P) / (e1^m · e2^n); remainder R = P − C·l1^m·l2^n.
        let lc = numer.iter().find(|(d, _)| *d == deg_p)?.1;
        let denom_lc = normalize(
            ctx,
            ctx.mul(&[int_pow(ctx, l1.e, l1.k), int_pow(ctx, l2.e, l2.k)]),
        );
        let c_quot = normalize(ctx, ctx.mul(&[lc, inv(ctx, denom_lc)]));
        // Expand l1^m · l2^n into monomial terms and subtract.
        let mono1 = [(0usize, l1.d), (1, l1.e)];
        let mono2 = [(0usize, l2.d), (1, l2.e)];
        let p1 = expand_terms(ctx, &mono1, l1.k)?;
        let p2 = expand_terms(ctx, &mono2, l2.k)?;
        let mut rem: std::collections::BTreeMap<usize, Vec<Atom<'a>>> =
            std::collections::BTreeMap::new();
        for (d1, c1) in &p1 {
            for (d2, c2) in &p2 {
                rem.entry(d1 + d2)
                    .or_default()
                    .push(normalize(ctx, ctx.mul(&[c_quot, *c1, *c2])));
            }
        }
        for (d, c) in &numer {
            rem.entry(*d)
                .or_default()
                .push(normalize(ctx, ctx.mul(&[ctx.num(-1), *c])));
        }
        // The top degree cancels by construction of C; symbolic inverse
        // pairs like `(b·d³)·(b·d³)⁻¹` don't fold, so drop it analytically.
        rem.remove(&total);
        numer = rem
            .into_iter()
            .map(|(d, cs)| (d, normalize(ctx, ctx.mul(&[ctx.num(-1), ctx.add(&cs)]))))
            .filter(|(_, c)| !is_zero(ctx, *c))
            .collect();
        if numer.is_empty() {
            // P = C·l1^m·l2^n exactly: the integral is C·x.
            let mut all = consts.clone();
            all.push(ctx.mul(&[c_quot, x]));
            return Some(normalize(ctx, ctx.mul(&all)));
        }
        deg_p = numer.iter().map(|(d, _)| *d).max()?;
        if deg_p >= total {
            return None;
        }
        quotient = Some(c_quot);
    } else if deg_p > total {
        return None;
    }
    let p_atom = normalize(
        ctx,
        ctx.add(
            &numer
                .iter()
                .map(|(d, c)| normalize(ctx, ctx.mul(&[*c, int_pow(ctx, x, *d as i64)])))
                .collect::<Vec<_>>(),
        ),
    );
    // A_i from the l1 side: body = P·l2^(−n), root r1 = −d1/e1.
    let r1 = normalize(ctx, ctx.mul(&[ctx.num(-1), l1.d, inv(ctx, l1.e)]));
    let r2 = normalize(ctx, ctx.mul(&[ctx.num(-1), l2.d, inv(ctx, l2.e)]));
    let body1 = normalize(ctx, ctx.mul(&[p_atom, ctx.pow(l2.base, ctx.num(-l2.k))]));
    let body2 = normalize(ctx, ctx.mul(&[p_atom, ctx.pow(l1.base, ctx.num(-l1.k))]));
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for i in 1..=l1.k {
        let a_i = pf_coeff(ctx, body1, l1.k - i, r1, l1.e, var)?;
        terms.push(normalize(
            ctx,
            ctx.mul(&[a_i, ctx.pow(l1.base, ctx.num(-i))]),
        ));
    }
    for j in 1..=l2.k {
        let b_j = pf_coeff(ctx, body2, l2.k - j, r2, l2.e, var)?;
        terms.push(normalize(
            ctx,
            ctx.mul(&[b_j, ctx.pow(l2.base, ctx.num(-j))]),
        ));
    }
    // The equal-degree edge contributes a constant quotient term.
    if let Some(c_quot) = quotient {
        terms.push(c_quot);
    }
    let mut all = consts.clone();
    all.push(ctx.add(&terms));
    let sum = normalize(ctx, ctx.mul(&all));
    if sum == expr {
        return None;
    }
    // Each term is a linear-denominator power (plus the quotient term) —
    // integrate through the chain and decline on residue.
    let r = integrate_raw(ctx, sum, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    Some(r)
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integral::integrate;
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
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    "atanh" => v.atanh(),
                    _ => return None,
                })
            }
        }
    }

    fn parse_norm<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        let e = ocas_parse::parse(ctx, s).unwrap();
        normalize(ctx, e)
    }

    fn assert_module_solves(input: &str, consts: &[(Symbol, f64)], samples: &[f64]) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, input);
        let result = integrate_quad_power(&ctx, expr, Symbol::new("x"))
            .unwrap_or_else(|| panic!("declined: {input}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue for {input}: {result}"
        );
        let d = crate::diff(&ctx, result, Symbol::new("x"));
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((Symbol::new("x"), xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(expr, &env).expect("eval integrand");
            let tol = 1e-5 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{input} at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
    }

    fn assert_module_declined(input: &str) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, input);
        assert!(
            integrate_quad_power(&ctx, expr, Symbol::new("x")).is_none(),
            "expected None for {input}"
        );
    }

    // ------------------------- L1 -------------------------

    #[test]
    fn l1_linear_denom_cube() {
        // ∫ (x² + 1)/(x + 2)³ dx — split into y = x+2 powers.
        assert_module_solves("(x^2 + 1)/(x + 2)^3", &[], &[0.3, 1.0, 2.5]);
    }

    #[test]
    fn l1_symbolic_binomial_numerator() {
        // Corpus shape: (d + e·x)^4/(a + b·x)^3 family — symbolic coeffs.
        let env = [
            (Symbol::new("d"), 1.2),
            (Symbol::new("e"), 0.6),
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
        ];
        assert_module_solves("(d + e*x)^4/(a + b*x)^3", &env, &[0.4, 1.0, 1.9]);
    }

    #[test]
    fn l1_n1_declines() {
        // n = 1 is the plain rational path, not ours.
        assert_module_declined("(x^2 + 1)/(x + 2)");
    }

    #[test]
    fn l1_trinomial_numerator_high_pow() {
        // Corpus shape: trinomial^3 over linear^9 (n > 6 cap raised).
        assert_module_solves("(2 + 3*x + x^2)^3/(1 + x)^9", &[], &[0.2, 0.7, 1.4]);
    }

    // ------------------------- L3 -------------------------

    #[test]
    fn l3_two_distinct_linear_denoms() {
        // ∫ x/((x+1)·(x+2)²) dx — partial fractions over two linear powers.
        assert_module_solves("x/((1 + x)*(2 + x)^2)", &[], &[0.3, 0.9, 1.8]);
    }

    #[test]
    fn l3_equal_degree_numeric() {
        // x⁴/((1+2x)(3+4x)³): the equal-degree edge (quotient constant).
        assert_module_solves("x^4/((1 + 2*x)*(3 + 4*x)^3)", &[], &[0.2, 0.6, 1.1]);
    }

    #[test]
    fn l3_symbolic_two_linear() {
        // Corpus shape: x⁴/((a + b·x)·(c + d·x)³).
        let env = [
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
            (Symbol::new("c"), 1.5),
            (Symbol::new("d"), 0.8),
        ];
        assert_module_solves("x^4/((a + b*x)*(c + d*x)^3)", &env, &[0.2, 0.6, 1.1]);
    }

    #[test]
    fn l3_proportional_declines() {
        // (2 + 4x) = 2·(1 + 2x): proportional factors are one factor.
        assert_module_declined("x/((1 + 2*x)*(2 + 4*x)^2)");
    }

    #[test]
    fn l3_corpus_linear_pair_powers() {
        // Corpus shape: (a + b·x)^10·(A + B·x)/(d + e·x)^19 → the l2 side is
        // (d+ex)^19, over the L3 cap — decline honestly.
        assert_module_declined("(a + b*x)^10*(A + B*x)/(d + e*x)^19");
    }

    // ------------------------- L2 -------------------------

    #[test]
    fn l2_inv_quad_squared() {
        // ∫ dx/(x² + 2x + 5)² — Q1 recurrence off the atan base.
        assert_module_solves("1/(x^2 + 2*x + 5)^2", &[], &[0.3, 1.0, 2.2]);
    }

    #[test]
    fn l2_inv_quad_cubed_symbolic() {
        // ∫ dx/(a + b·x + c·x²)³ — symbolic coefficients.
        let env = [
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
            (Symbol::new("c"), 1.5),
        ];
        assert_module_solves("1/(a + b*x + c*x^2)^3", &env, &[0.3, 1.0, 2.0]);
    }

    #[test]
    fn l2_x_over_quad_squared() {
        // ∫ x/(x² + 3x + 7)² dx — J(1,2) branch.
        assert_module_solves("x/(x^2 + 3*x + 7)^2", &[], &[0.4, 1.1, 2.0]);
    }

    #[test]
    fn l2_x3_over_quad_cubed() {
        // ∫ x³/(2 + x + x²)³ dx — full J(3,3) recursion.
        assert_module_solves("x^3/(2 + x + x^2)^3", &[], &[0.4, 1.0, 1.8]);
    }

    #[test]
    fn l2_corpus_quartic_over_quad_fourth() {
        // Corpus shape: (d + e·x)^4/(a + b·x + c·x²)^4.
        let env = [
            (Symbol::new("d"), 1.2),
            (Symbol::new("e"), 0.6),
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
            (Symbol::new("c"), 1.5),
        ];
        assert_module_solves("(d + e*x)^4/(a + b*x + c*x^2)^4", &env, &[0.3, 0.9, 1.7]);
    }

    #[test]
    fn l2_repeated_root_declines() {
        // q = (x+1)² has D = 0 — declines (the rational backend owns it).
        assert_module_declined("1/(1 + 2*x + x^2)^3");
    }

    #[test]
    fn chain_end_to_end() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "x^3/(2 + x + x^2)^3");
        let r = integrate(&ctx, expr, Symbol::new("x"));
        assert!(
            !r.to_string().contains("Integral"),
            "chain left residue: {r}"
        );
    }

    /// `rubi-01458`, `rubi-00597`, `rubi-01308` and `rubi-00367` were listed
    /// against this module; it is not their source. All four are *outside* its
    /// L1/L2/L3 shape classes and must decline here, so the module can never
    /// be the origin of a wrong answer for them:
    ///
    /// - `csc(e+f·x)^5·(a + b·sec(e+f·x)²)²` is transcendental in `x` (the
    ///   module only matches polynomial numerators over polynomial bases), and
    ///   the pipeline reports it as an honest `Integral(...)` residue;
    /// - the two `(1−2x)^k/((2+3x)^m·(3+5x)^n)` rationals exceed the L3
    ///   combined-power and numerator-degree caps, so they fall to
    ///   `rational` (which solves them exactly — the corpus oracle's
    ///   `Mismatch` for those two is float64 cancellation in the *oracle's*
    ///   5-point stencil, not a wrong antiderivative);
    /// - `1/((a+b·x)·(a²−b²·x²))` has a repeated (non-squarefree) quadratic
    ///   factor, which L2 deliberately leaves to `rational`.
    #[test]
    fn non_owned_corpus_shapes_decline() {
        assert_module_declined("csc(e + f*x)^5*(a + b*sec(e + f*x)^2)^2");
        assert_module_declined("(1 - 2*x)^3/((2 + 3*x)^7*(3 + 5*x)^2)");
        assert_module_declined("(1 - 2*x)^2/((2 + 3*x)^6*(3 + 5*x)^3)");
        assert_module_declined("1/((a + b*x)*(a^2 - b^2*x^2))");
    }
}
