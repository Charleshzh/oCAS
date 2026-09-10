//! Trig-denominator power reduction, linear-numerator decomposition, and
//! polynomial×trig closed forms.
//!
//! Three mechanisms, each declining honestly (`None`) when inapplicable:
//!
//! - **T1 [`integrate_trig_denom`]** — `∫ dx/(a + b·T(u))^n` for
//!   `T ∈ {sin, cos}`, `u = c + d·x` linear (symbolic coefficients allowed),
//!   `a, b` constant w.r.t. `x`, integer `2 ≤ n ≤ 6`. The classical
//!   recurrence, derived by differentiating `S(u)/D^(n−1)` with
//!   `D = a + b·T(u)` and `S` the complementary trig function:
//!
//!   ```text
//!   J_n = [ ε·b·S(u)/(d·D^(n−1)) + a(2n−3)·J_{n−1} − (n−2)·J_{n−2} ] / [ (n−1)(a²−b²) ]
//!   with  ε = −1, S = sin  when T = cos;   ε = +1, S = cos  when T = sin.
//!   ```
//!
//!   The `n = 1` base re-enters the chain via [`integrate_raw`] (Weierstrass
//!   and the symbolic-rational backend solve it); any `Integral` residue
//!   there declines the whole reduction. The singular case `a² = b²`
//!   declines — exactly when `a, b` are numeric, syntactically via
//!   `collect_terms` when symbolic.
//!
//! - **T2 [`integrate_trig_num_linear`]** — `(A + B·T(u))/(a + b·T(u))^n`
//!   splits into `A·J_n + B·K_n` with `K_n = ∫ T/D^n = (J_{n−1} − a·J_n)/b`,
//!   reusing T1. The numerator's trig function and argument must match the
//!   denominator's exactly; a product of two denominator powers declines.
//!
//! - **T3 [`integrate_poly_trig`]** — `∫ x^m·T(a·x + b) dx` (integer
//!   `1 ≤ m ≤ 6`) and `∫ (c + d·x)^k·T(a·x + b) dx` (`k ≤ 4`, binomially
//!   expanded into monomial terms) as direct closed forms by tabular
//!   reduction:
//!
//!   ```text
//!   ∫ x^m sin(ax+b) dx = Σ_{k=0..m} fall(m,k)·x^(m−k)·a^(−k−1)·s_k,
//!       s_k = −cos, +sin, +cos, −sin, …  (k mod 4 = 0,1,2,3)
//!   ∫ x^m cos(ax+b) dx = Σ_{k=0..m} fall(m,k)·x^(m−k)·a^(−k−1)·c_k,
//!       c_k = +sin, +cos, −sin, −cos, …
//!   ```
//!
//!   The closed form never re-enters the chain and contains no `Integral`,
//!   so it cannot ping-pong with the parts heuristic (the 0.27 lesson).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::{contains_integral, integrate_raw, is_constant, linear_form};

/// Cap on the denominator power `n` in T1/T2.
const MAX_DENOM_POW: i64 = 6;
/// Cap on the monomial degree `m` in T3.
const MAX_POLY_DEG: i64 = 6;
/// Cap on `k` for `(c + d·x)^k` expansion in T3.
const MAX_LIN_POW: i64 = 4;
/// Node budget for the input integrand.
const MAX_NODES: usize = 200;

/// Combined entry: T1 (denominator power) → T2 (linear numerator) → T3
/// (polynomial×trig closed form).
pub(crate) fn integrate_trig_reduction<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    integrate_trig_denom(ctx, expr, var)
        .or_else(|| integrate_trig_num_linear(ctx, expr, var))
        .or_else(|| integrate_poly_trig(ctx, expr, var))
}

// =========================================================================
// Shared matcher: one `(a + b·T(u))^(−n)` factor in a product
// =========================================================================

/// A product split around a single trig-denominator power.
struct DenomMatch<'a> {
    /// Constant (w.r.t. `var`) leftover factors.
    rest: Vec<Atom<'a>>,
    /// Non-constant factors that are not the denominator power.
    other: Vec<Atom<'a>>,
    /// `a + b·T(u)` as found.
    base: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    sin: bool,
    u: Atom<'a>,
    /// `du/dx` (nonzero).
    du: Atom<'a>,
    /// The denominator power `n ≥ 1`.
    n: i64,
    /// `a² − b²`, verified nonzero.
    disc: Atom<'a>,
}

fn match_trig_denom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<DenomMatch<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest = Vec::new();
    let mut other = Vec::new();
    let mut found: Option<(Atom<'a>, Atom<'a>, Atom<'a>, bool, Atom<'a>, Atom<'a>, i64)> = None;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if let Some((base, n)) = as_recip_pow(f)
            && let Some((a, b, sin, u, du)) = split_trig_base(ctx, base, var)
        {
            // A second trig-denominator power is out of scope (decline).
            if found.is_some() {
                return None;
            }
            found = Some((base, a, b, sin, u, du, n));
            continue;
        }
        other.push(f);
    }
    let (base, a, b, sin, u, du, n) = found?;
    let disc = discriminant(ctx, a, b)?;
    Some(DenomMatch {
        rest,
        other,
        base,
        a,
        b,
        sin,
        u,
        du,
        n,
        disc,
    })
}

/// `(base, n)` for `(base^n)^(−1)` or `base^(−n)` with integer `n ≥ 1` —
/// the two shapes the parser/normalizer produce for `1/base^n`. The nested
/// form is tried first: `(base^2)^(−1)` must read as `n = 2`, not as
/// `Pow(base^2, −1)` with `n = 1`.
fn as_recip_pow(f: Atom<'_>) -> Option<(Atom<'_>, i64)> {
    let AtomNode::Pow(b, e) = f.node() else {
        return None;
    };
    if let (AtomNode::Pow(bb, ee), AtomNode::Num(-1)) = (b.node(), e.node())
        && let AtomNode::Num(n) = ee.node()
        && *n >= 1
    {
        return Some((*bb, *n));
    }
    if let AtomNode::Num(n) = e.node()
        && *n <= -1
    {
        return Some((*b, n.checked_neg()?));
    }
    None
}

/// Decomposition of one additive term.
enum TermSplit<'a> {
    Constant,
    /// `coeff · T(u)` with `coeff` constant and `T ∈ {sin, cos}`.
    Trig {
        coeff: Atom<'a>,
        sin: bool,
        u: Atom<'a>,
    },
    Other,
}

/// Decompose one additive term: constant, or constants × a single bare
/// `sin(u)`/`cos(u)` factor.
fn split_trig_term<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>, var: Symbol) -> TermSplit<'a> {
    if is_constant(term, var) {
        return TermSplit::Constant;
    }
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut trig: Option<(bool, Atom<'a>)> = None;
    for f in factors {
        if is_constant(f, var) {
            coeff.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Fun(name, args)
                if args.len() == 1 && (name.as_str() == "sin" || name.as_str() == "cos") =>
            {
                if trig.is_some() {
                    return TermSplit::Other;
                }
                trig = Some((name.as_str() == "sin", args[0]));
            }
            _ => return TermSplit::Other,
        }
    }
    let Some((sin, u)) = trig else {
        return TermSplit::Other;
    };
    let coeff = if coeff.is_empty() {
        ctx.num(1)
    } else {
        normalize(ctx, ctx.mul(&coeff))
    };
    TermSplit::Trig { coeff, sin, u }
}

/// Match `base` as `a + b·T(u)`: constant terms plus exactly one trig term
/// whose argument is linear in `var` with nonzero slope. Returns
/// `(a, b, sin, u, du/dx)`.
fn split_trig_base<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, bool, Atom<'a>, Atom<'a>)> {
    let AtomNode::Add(args) = base.node() else {
        return None;
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut trig: Option<(Atom<'a>, bool, Atom<'a>)> = None;
    for t in args.iter() {
        match split_trig_term(ctx, *t, var) {
            TermSplit::Constant => consts.push(*t),
            TermSplit::Trig { coeff, sin, u } => {
                if trig.is_some() {
                    return None;
                }
                trig = Some((coeff, sin, u));
            }
            TermSplit::Other => return None,
        }
    }
    let (b, sin, u) = trig?;
    let a = if consts.is_empty() {
        ctx.num(0)
    } else {
        normalize(ctx, ctx.add(&consts))
    };
    let (du, _phase) = linear_form(ctx, u, var)?;
    if matches!(du.node(), AtomNode::Num(0)) {
        return None;
    }
    Some((a, b, sin, u, du))
}

/// `a² − b²`, folded for numeric `a`/`b`; `None` when it is zero (exactly
/// for numeric values, syntactically via `collect_terms` for symbolic
/// ones) — the singular case the recurrence cannot divide by.
fn discriminant<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>, b: Atom<'a>) -> Option<Atom<'a>> {
    if let (Some(av), Some(bv)) = (numeric(a), numeric(b)) {
        let d = (av as i128)
            .checked_mul(av as i128)?
            .checked_sub((bv as i128).checked_mul(bv as i128)?)?;
        if d == 0 {
            return None;
        }
        return i64::try_from(d).ok().map(|v| ctx.num(v));
    }
    let sq_a = square(ctx, a)?;
    let sq_b = square(ctx, b)?;
    let d = crate::ode::util::collect_terms(ctx, ctx.add(&[sq_a, ctx.mul(&[ctx.num(-1), sq_b])]));
    if matches!(d.node(), AtomNode::Num(0)) {
        return None;
    }
    Some(d)
}

fn square<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> Option<Atom<'a>> {
    match e.node() {
        AtomNode::Num(v) => Some(ctx.num(v.checked_mul(*v)?)),
        _ => Some(ctx.pow(e, ctx.num(2))),
    }
}

// =========================================================================
// T1: ∫ dx/(a + b·T(u))^n
// =========================================================================

/// T1 entry: `∫ C·dx/(a + b·T(u))^n` for `2 ≤ n ≤ 6` (see module docs).
pub(crate) fn integrate_trig_denom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    let m = match_trig_denom(ctx, expr, var)?;
    if !(2..=MAX_DENOM_POW).contains(&m.n) || !m.other.is_empty() {
        return None;
    }
    let core = denom_integral(ctx, &m, m.n, var)?;
    let mut factors = m.rest.clone();
    factors.push(core);
    Some(normalize(ctx, ctx.mul(&factors)))
}

/// `J_n = ∫ dx/(a + b·T(u))^n` by the recurrence in the module docs.
/// `n = 1` re-enters the chain (Weierstrass) and declines on residue.
fn denom_integral<'a>(
    ctx: &'a AtomArena<'a>,
    m: &DenomMatch<'a>,
    n: i64,
    var: Symbol,
) -> Option<Atom<'a>> {
    debug_assert!((1..=MAX_DENOM_POW).contains(&n));
    if n == 1 {
        let g = ctx.pow(m.base, ctx.num(-1));
        let r = integrate_raw(ctx, g, var, 0, true, 0, 0);
        if contains_integral(r) {
            return None;
        }
        return Some(r);
    }
    let j1 = denom_integral(ctx, m, n - 1, var)?;
    let j2 = if n == 2 {
        None
    } else {
        Some(denom_integral(ctx, m, n - 2, var)?)
    };
    // Boundary term ε·b·S(u)/(du·D^(n−1)): ε = −1, S = sin for T = cos;
    // ε = +1, S = cos for T = sin.
    let (sign, other_sin) = if m.sin { (1, false) } else { (-1, true) };
    let boundary = ctx.mul(&[
        ctx.num(sign),
        m.b,
        ctx.fun(if other_sin { "sin" } else { "cos" }, &[m.u]),
        ctx.pow(m.base, ctx.num(-(n - 1))),
        ctx.pow(m.du, ctx.num(-1)),
    ]);
    let mut terms = vec![boundary, ctx.mul(&[m.a, ctx.num(2 * n - 3), j1])];
    if let Some(j2) = j2 {
        terms.push(ctx.mul(&[ctx.num(-(n - 2)), j2]));
    }
    let bracket = normalize(ctx, ctx.add(&terms));
    let denom = normalize(ctx, ctx.mul(&[ctx.num(n - 1), m.disc]));
    Some(normalize(
        ctx,
        ctx.mul(&[bracket, ctx.pow(denom, ctx.num(-1))]),
    ))
}

// =========================================================================
// T2: (A + B·T(u))/(a + b·T(u))^n
// =========================================================================

/// T2 entry: `(A + B·T(u))/(a + b·T(u))^n` for `2 ≤ n ≤ 6`, with the
/// numerator's trig function and argument matching the denominator's.
pub(crate) fn integrate_trig_num_linear<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    let m = match_trig_denom(ctx, expr, var)?;
    if !(2..=MAX_DENOM_POW).contains(&m.n) || m.other.len() != 1 {
        return None;
    }
    let u_norm = normalize(ctx, m.u);
    let (num_a, num_b) = split_numerator(ctx, m.other[0], var, m.sin, u_norm)?;
    if matches!(num_b.node(), AtomNode::Num(0)) {
        return None; // constant numerator: T1's shape
    }
    let jn = denom_integral(ctx, &m, m.n, var)?;
    let jn1 = denom_integral(ctx, &m, m.n - 1, var)?;
    // (A + B·T)/D^n = A·J_n + (B/b)·(J_{n−1} − a·J_n), from T = (D − a)/b.
    let term_a = ctx.mul(&[num_a, jn]);
    let inner = normalize(ctx, ctx.add(&[jn1, ctx.mul(&[ctx.num(-1), m.a, jn])]));
    let term_b = ctx.mul(&[num_b, inv(ctx, m.b), inner]);
    let core = normalize(ctx, ctx.add(&[term_a, term_b]));
    let mut factors = m.rest.clone();
    factors.push(core);
    Some(normalize(ctx, ctx.mul(&factors)))
}

/// Split the numerator as `A + B·T(u)`: constant terms vs. coefficients of
/// trig factors matching (`sin`, `u_norm`). Every term must classify.
fn split_numerator<'a>(
    ctx: &'a AtomArena<'a>,
    num: Atom<'a>,
    var: Symbol,
    sin: bool,
    u_norm: Atom<'a>,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let terms: Vec<Atom<'a>> = match num.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![num],
    };
    let mut a_terms = Vec::new();
    let mut b_terms = Vec::new();
    for t in terms {
        match split_trig_term(ctx, t, var) {
            TermSplit::Constant => a_terms.push(t),
            TermSplit::Trig { coeff, sin: s2, u } => {
                if s2 != sin || normalize(ctx, u) != u_norm {
                    return None;
                }
                b_terms.push(coeff);
            }
            TermSplit::Other => return None,
        }
    }
    let a = if a_terms.is_empty() {
        ctx.num(0)
    } else {
        normalize(ctx, ctx.add(&a_terms))
    };
    let b = if b_terms.is_empty() {
        ctx.num(0)
    } else {
        normalize(ctx, ctx.add(&b_terms))
    };
    Some((a, b))
}

// =========================================================================
// T3: polynomial × trig closed form
// =========================================================================

/// T3 entry: `∫ C·x^m·T(a·x + b) dx` (`1 ≤ m ≤ 6`) and
/// `∫ C·(c + d·x)^k·T(a·x + b) dx` (`1 ≤ k ≤ 4`) as direct closed forms.
pub(crate) fn integrate_poly_trig<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest = Vec::new();
    let mut poly: Option<Vec<(Atom<'a>, i64)>> = None;
    let mut trig: Option<(bool, Atom<'a>, Atom<'a>)> = None;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        match f.node() {
            AtomNode::Fun(name, args)
                if args.len() == 1 && (name.as_str() == "sin" || name.as_str() == "cos") =>
            {
                if trig.is_some() {
                    return None;
                }
                let (slope, _phase) = linear_form(ctx, args[0], var)?;
                if matches!(slope.node(), AtomNode::Num(0)) {
                    return None;
                }
                trig = Some((name.as_str() == "sin", args[0], slope));
            }
            AtomNode::Var(v) if *v == var => {
                if poly.is_some() {
                    return None;
                }
                poly = Some(vec![(ctx.num(1), 1)]);
            }
            AtomNode::Pow(b, e) => {
                if poly.is_some() {
                    return None;
                }
                let AtomNode::Num(k) = e.node() else {
                    return None;
                };
                if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                    if !(1..=MAX_POLY_DEG).contains(k) {
                        return None;
                    }
                    poly = Some(vec![(ctx.num(1), *k)]);
                } else {
                    if !(1..=MAX_LIN_POW).contains(k) {
                        return None;
                    }
                    poly = Some(expand_linear_pow(ctx, *b, *k, var)?);
                }
            }
            _ => return None,
        }
    }
    let poly = poly?;
    let (sin, u, slope) = trig?;
    let mut terms = Vec::new();
    for (coef, deg) in poly {
        let cf = closed_form(ctx, deg, sin, u, slope, var)?;
        match cf.node() {
            AtomNode::Add(args) => {
                for t in args.iter() {
                    terms.push(ctx.mul(&[coef, *t]));
                }
            }
            _ => terms.push(ctx.mul(&[coef, cf])),
        }
    }
    let core = normalize(ctx, ctx.add(&terms));
    let mut factors = rest;
    factors.push(core);
    Some(normalize(ctx, ctx.mul(&factors)))
}

/// `(c + d·x)^k` as `(coeff, degree)` monomial pairs via the binomial
/// theorem (zero coefficients dropped).
fn expand_linear_pow<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    k: i64,
    var: Symbol,
) -> Option<Vec<(Atom<'a>, i64)>> {
    let (d, c) = linear_form(ctx, base, var)?;
    if matches!(d.node(), AtomNode::Num(0)) {
        return None;
    }
    let mut out = Vec::new();
    for j in 0..=k {
        let coef = normalize(
            ctx,
            ctx.mul(&[
                ctx.num(binom(k, j)),
                int_pow(ctx, c, k - j)?,
                int_pow(ctx, d, j)?,
            ]),
        );
        if matches!(coef.node(), AtomNode::Num(0)) {
            continue;
        }
        out.push((coef, j));
    }
    if out.is_empty() { None } else { Some(out) }
}

/// Closed form of `∫ x^m·T(u) dx` (`T` = sin when `sin`, else cos) with
/// `u = slope·x + phase`, by m-fold tabular reduction (see module docs).
/// Emits the full sum directly — no chain re-entry.
fn closed_form<'a>(
    ctx: &'a AtomArena<'a>,
    m: i64,
    sin: bool,
    u: Atom<'a>,
    slope: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    debug_assert!((0..=MAX_POLY_DEG).contains(&m));
    let x = ctx.var(var.as_str());
    let mut terms = Vec::new();
    let mut fall: i64 = 1; // falling factorial m!/(m−k)!
    for k in 0..=m {
        if k > 0 {
            fall = fall.checked_mul(m - k + 1)?;
        }
        let out_sin = if sin { k % 2 == 1 } else { k % 2 == 0 };
        let neg = if sin {
            matches!(k % 4, 0 | 3)
        } else {
            k % 4 >= 2
        };
        let acoef = match slope.node() {
            AtomNode::Num(av) => {
                let p = av.checked_pow((k + 1) as u32)?;
                ctx.pow(ctx.num(p), ctx.num(-1))
            }
            _ => ctx.pow(slope, ctx.num(-(k + 1))),
        };
        let mut fs: Vec<Atom<'a>> = Vec::new();
        if neg {
            fs.push(ctx.num(-1));
        }
        if fall != 1 {
            fs.push(ctx.num(fall));
        }
        let deg = m - k;
        if deg > 0 {
            fs.push(int_pow(ctx, x, deg)?);
        }
        fs.push(acoef);
        fs.push(ctx.fun(if out_sin { "sin" } else { "cos" }, &[u]));
        terms.push(normalize(ctx, ctx.mul(&fs)));
    }
    Some(normalize(ctx, ctx.add(&terms)))
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

// =========================================================================
// Helpers
// =========================================================================

/// Integer value of a `Num` atom.
fn numeric(e: Atom<'_>) -> Option<i64> {
    if let AtomNode::Num(n) = e.node() {
        Some(*n)
    } else {
        None
    }
}

/// `a^(−1)`, folding the trivial units.
fn inv<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    match a.node() {
        AtomNode::Num(1) => ctx.num(1),
        AtomNode::Num(-1) => ctx.num(-1),
        _ => ctx.pow(a, ctx.num(-1)),
    }
}

/// `base^e` for integer `e ≥ 0`, folding numeric bases and `e ∈ {0, 1}`.
/// Negative exponents (never produced here) fall back to a `Pow` node.
fn int_pow<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, e: i64) -> Option<Atom<'a>> {
    if e < 0 {
        return Some(ctx.pow(base, ctx.num(e)));
    }
    match e {
        0 => Some(ctx.num(1)),
        1 => Some(base),
        _ => {
            if let AtomNode::Num(v) = base.node() {
                Some(ctx.num(v.checked_pow(e as u32)?))
            } else {
                Some(ctx.pow(base, ctx.num(e)))
            }
        }
    }
}

/// Total node count of the expression tree (saturating).
fn node_count(expr: Atom<'_>) -> usize {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => 1,
        AtomNode::Pow(b, e) => node_count(*b)
            .saturating_add(node_count(*e))
            .saturating_add(1),
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => args
            .iter()
            .fold(1usize, |acc, a| acc.saturating_add(node_count(*a))),
    }
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

    /// Numeric f64 evaluator for test verification (handles the operators
    /// the reduction results produce, incl. the Weierstrass atan/log base).
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
                    "sec" => v.cos().recip(),
                    "csc" => v.sin().recip(),
                    "cot" => v.tan().recip(),
                    "exp" => v.exp(),
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    _ => return None,
                })
            }
        }
    }

    /// Run the mechanism, require `Some` without residue, and check
    /// `diff(result) == integrand` numerically at the given sample points.
    /// Sample points must keep `a + b·T(u)` away from zero and `u/2` away
    /// from the `tan` poles (the Weierstrass base answer is a sawtooth).
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_trig_reduction(ctx, integrand, var).expect("mechanism declined");
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
            integrate_trig_reduction(ctx, expr, Symbol::new("x")).is_none(),
            "expected None for {input}"
        );
    }

    // ------------------------- T1 -------------------------

    #[test]
    fn t1_cos_n2_numeric() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "1/(2 + cos(x))^2");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t1_cos_n3_corpus() {
        // Corpus sample: symbolic linear argument, numeric a/b.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "1/(-5 + 3*cos(c + d*x))^3");
        let env = [(Symbol::new("c"), 0.4), (Symbol::new("d"), 1.3)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.5, 0.9]);
    }

    #[test]
    fn t1_sin_n2_numeric() {
        // a = 2, b = 1: the Weierstrass base 1/(2+sin x) is chain-correct.
        // (a ∤ 2b bases like 1/(3+2 sin x) hit the rational.rs
        // sqrt_positive_rational missing-/q bug — main-line fix.)
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "1/(2 + sin(x))^2");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t1_cos_n4_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "1/(a + b*cos(c + d*x))^4");
        let env = [
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
            (Symbol::new("c"), 0.3),
            (Symbol::new("d"), 0.8),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.6, 1.0]);
    }

    #[test]
    fn t1_singular_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // a² = b² exactly (numeric) and syntactically (symbolic a = b).
        assert_declined(&ctx, "1/(2 + 2*cos(x))^2");
        assert_declined(&ctx, "1/(a + a*sin(c + d*x))^2");
    }

    #[test]
    fn t1_out_of_scope_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // n = 1 is the chain's (Weierstrass); n = 7 exceeds the budget;
        // a Mul base (no additive constant term) is not the T1 shape.
        assert_declined(&ctx, "1/(2 + cos(x))");
        assert_declined(&ctx, "1/(2 + cos(x))^7");
        assert_declined(&ctx, "1/(2*cos(x))^2");
    }

    // ------------------------- T2 -------------------------

    #[test]
    fn t2_cos_linear_num() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(1 + 3*cos(x))/(2 + cos(x))^2");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t2_sin_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(A + B*cos(c + d*x))/(a + b*cos(c + d*x))^3");
        let env = [
            (Symbol::new("A"), 1.5),
            (Symbol::new("B"), -0.5),
            (Symbol::new("a"), 2.0),
            (Symbol::new("b"), 0.5),
            (Symbol::new("c"), 0.3),
            (Symbol::new("d"), 1.1),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn t2_bare_trig_numerator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "sin(x)/(2 + sin(x))^2");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t2_two_denominator_powers_decline() {
        // Corpus sample: a product of two denominator powers is out of
        // scope (and a = b in the first factor is singular anyway).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_declined(
            &ctx,
            "(A + B*sin(e + f*x))/((a + a*sin(e + f*x))^3*(c - c*sin(e + f*x))^2)",
        );
    }

    #[test]
    fn t2_mismatched_trig_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // sin numerator over a cos denominator is a u-substitution shape,
        // not this decomposition.
        assert_declined(&ctx, "sin(x)/(2 + cos(x))^2");
    }

    // ------------------------- T3 -------------------------

    #[test]
    fn t3_x2_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "x^2*cos(x)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t3_x4_sin_symbolic_slope() {
        // Beyond the parts budget (PARTS_MAX_DEPTH = 2): the 0.27.1 gap.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "x^4*sin(a + b*x)");
        let env = [(Symbol::new("a"), 0.4), (Symbol::new("b"), 1.3)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t3_x6_cos_numeric() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "x^6*cos(2*x - 1)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t3_linpow4_cos() {
        // The post-trig_reduce corpus term shape: `(c+d*x)^k·T(a+b*x)`.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(c + d*x)^4*cos(a + b*x)");
        let env = [
            (Symbol::new("a"), 1.1),
            (Symbol::new("b"), 0.6),
            (Symbol::new("c"), 0.7),
            (Symbol::new("d"), -0.4),
        ];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t3_linpow3_sin_numeric() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(1 + 2*x)^3*sin(3*x)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn t3_out_of_scope_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Two trig factors: trig_reduce's product-to-sum owns this shape.
        assert_declined(&ctx, "x*sin(x)*cos(x)");
        // A trig power: rules C12/C14 and trig_reduce own this shape.
        assert_declined(&ctx, "x*sin(x)^2");
        // Degree beyond the budget.
        assert_declined(&ctx, "x^7*sin(x)");
    }
}
