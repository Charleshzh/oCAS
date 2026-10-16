//! Inverse-trig/hyperbolic integration mechanisms (0.27.1 Phase 1E).
//!
//! Three mechanisms for integrands built from the eleven inverse functions
//! `asin/acos/atan/asinh/acosh/atanh` (0.27.1) plus
//! `acsc/asec/acoth/asech/acsch` (0.30.0 Wave 5c):
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
//! - **M4 inverse-function substitution** (0.30.0 Wave 5,
//!   `docs/planning/INV_SUBST_SPEC_CN.md`): `(a + b·F(w))^k` times a
//!   polynomial weight, `w = ρ + σ·x` affine, integer `1 ≤ k ≤ 6`, the six
//!   main heads plus (W5c) the five reciprocal heads, plus (W5b) at most one
//!   companion radical
//!   `(s·E(w))^(m/2)`, integer `|m| ≤ 5`, absorbed as
//!   `s^(m/2)·kernel_head^m·σ` with the sheet factor σ = +1 on the principal
//!   band. `u = F(w)` maps the integrand to
//!   `P(u)·T(u)`; trig products reduce to multiple angles and re-enter the
//!   chain, hyperbolic families use the M2-style tabular closed form
//!   directly, and the ratio families (atan/atanh and, in the reciprocal
//!   pair, acsc/asec/acoth/asech/acsch) use a dedicated parts/
//!   `R`-reduction tabular. The back-substitution is a whitelist rewrite
//!   (`u → F(w)`, `G(u) → w`, `G′(u) → +√E(w)` on the principal band,
//!   multiple angles via Chebyshev); any residue, leftover `u`, or
//!   `F(G(·))` composite
//!   declines the whole shape. See the M4 section docs for the real-branch
//!   semantics, the Wave 5c capability record, and the deviations from the
//!   spec.
//!
//! Three extra passes run around that core:
//!
//! - **P1 inverse-composition cancellation**: a bounded, idempotent bottom-up
//!   pre-pass rewriting `atanh(tanh(u))`, `asinh(sinh(u))`, `acoth(coth(u))`,
//!   `log(exp(u))` (unconditional on the reals) and `atan(tan(u))`,
//!   `acosh(cosh(u))` (only on a certified principal band) to `u`, then
//!   re-integrating the expanded rewrite (see `integrate_expanded_terms`).
//!   This is what makes `atanh(tanh(a + b·x))^k·P(x)` shapes elementary.
//! - **P2 branch-shifted composition kernels**: `acoth(tanh(u))` equals
//!   `u ± iπ/2` and `atan(tan(u))` equals `u − kπ` on the k-th band, so
//!   cancelling them to `u` would be wrong for `k > 1`. Treating the
//!   composition as an opaque kernel `K` whose derivative is the constant
//!   slope `σ` of `u`, substituting `K = σ·x + w` (fresh constant `w`),
//!   integrating in `x` and mapping `w` back to `K − σ·x` is branch-free and
//!   exact.
//! - **P3 integral representation / companion radicands**: `M1` restricted to
//!   zero-intercept arguments; `P3` accepts an intercept and a *rational*
//!   inverse-function exponent, but only after verifying that the radical is
//!   the derivative-shaped companion `s·E(u)` of the family at the *same*
//!   argument `u` (so `atan(u)/sqrt(1 + u²)`, whose companion exponent would
//!   have to be `−1`, is declined rather than mis-integrated).
//!
//! Sign convention: a *numeric* radicand normalization `s` under a
//! square-root power must be positive; a *symbolic* `s` is assumed
//! positive (the same convention rule family G uses for `asin(x/a)`),
//! which is what lets corpus radicands like `(d + c²·d·x²)` through.
//!
//! The remaining chain re-entries go through `integrate_raw` and are declined
//! on any `Integral` residue; the P1/P2 rewrites are expanded first so that
//! `integrate_expanded_terms` can finish them with the module's own power rule
//! and never spend the pipeline's global chain-entry budget. Outputs keep the
//! inverse function inside `sinh`/`cosh`/power arguments, so they never
//! re-match this module's own top-level factor patterns (idempotent under
//! chain re-entry).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    contains_integral, contains_symbol, int_pow, integrate_raw, inv, is_constant, linear_form,
    node_count, pick_subst_symbol, rat_atom, replace_symbol,
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
/// P1 cap on the nodes visited by the cancellation pre-pass.
const MAX_CANCEL_BUDGET: usize = 800;
/// P2 cap on the `K^m` kernel exponent (`|m|`), and on the polynomial
/// cofactor degree and `x^(−j)` reciprocal exponent.
const MAX_KERNEL_POW: i64 = 6;
const MAX_COFACTOR_DEG: usize = 6;
const MAX_RECIP_EXP: i64 = 8;
/// P2 cap on the node count of the proxy integrand (before and after the
/// bounded expansion).
const MAX_KERNEL_NODES: usize = 400;
/// P3 cap on the denominator `q` of a rational inverse-function exponent.
const MAX_RAT_Q: i64 = 4;

/// The eleven inverse functions covered by this module: the six main heads of
/// M1–M3/P1–P3 (0.27.1–0.30.0 W5a/W5b) and the five reciprocal heads of
/// 0.30.0 Wave 5c (`acsc`/`asec`/`acoth`/`asech`/`acsch`, spec
/// `INV_SUBST_SPEC_CN.md` §2.1/§2.6).
///
/// The reciprocal heads are owned by M4 alone: M1/P3 decline them
/// ([`family_kernel_exp`]) because their derivatives carry the family's
/// *reciprocal* forward head (`acsc′ = −1/(w·√(w²−1))` versus
/// `asec′ = +1/(w·√(w²−1))`), which the kernel-derivative power rule of M1
/// does not model.
#[derive(Clone, Copy, PartialEq, Eq)]
enum InvFun {
    Asin,
    Acos,
    Atan,
    Asinh,
    Acosh,
    Atanh,
    Acsc,
    Asec,
    Acoth,
    Asech,
    Acsch,
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
            InvFun::Acsc => "acsc",
            InvFun::Asec => "asec",
            InvFun::Acoth => "acoth",
            InvFun::Asech => "asech",
            InvFun::Acsch => "acsch",
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
        "acsc" => InvFun::Acsc,
        "asec" => InvFun::Asec,
        "acoth" => InvFun::Acoth,
        "asech" => InvFun::Asech,
        "acsch" => InvFun::Acsch,
        _ => return None,
    })
}

/// Kernel exponent shape: `R^(−1/2)` (also `1/sqrt(R)`) or `R^(−1)`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum KernelExp {
    InvSqrt,
    Recip,
}

/// The kernel exponent shape expected per family, or `None` for a family the
/// kernel-derivative power rule (M1) and the companion-radicand pass (P3) must
/// not claim.
///
/// The Wave 5c reciprocal heads decline here: `acsc`/`asec`/`acsch` carry an
/// extra `|w|` in `F′(w)` (`d/dw acsc(w) = −1/(|w|√(w²−1))`), so a kernel
/// factor `E(w)^e` alone is *not* the derivative companion M1's chain rule
/// needs; `acoth`'s exponent shape does match (`1/(1−w²) = −E(w)^(−1)`) but
/// M4 owns the whole reciprocal-head shape (its substitution covers the same
/// integrands with a polynomial weight, so one owner keeps the answer shape
/// single-valued).
fn family_kernel_exp(f: InvFun) -> Option<KernelExp> {
    Some(match f {
        InvFun::Asin | InvFun::Acos | InvFun::Asinh | InvFun::Acosh => KernelExp::InvSqrt,
        InvFun::Atan | InvFun::Atanh => KernelExp::Recip,
        InvFun::Acsc | InvFun::Asec | InvFun::Acoth | InvFun::Asech | InvFun::Acsch => {
            return None;
        }
    })
}

/// `s_F` from `F′(u) = s_F·E(u)^e`, or `None` where that shape does not exist
/// (the Wave 5c reciprocal heads; see [`family_kernel_exp`]). Only `acos`
/// differentiates with a minus sign.
fn family_sign(f: InvFun) -> Option<i64> {
    Some(match f {
        InvFun::Acos => -1,
        InvFun::Asin | InvFun::Atan | InvFun::Asinh | InvFun::Acosh | InvFun::Atanh => 1,
        InvFun::Acsc | InvFun::Asec | InvFun::Acoth | InvFun::Asech | InvFun::Acsch => {
            return None;
        }
    })
}

/// The kernel sum `E` per family in terms of the slope `σ` of `u = σ·x`:
/// `1−σ²x²` (asin/acos/atanh/asech), `1+σ²x²` (atan/asinh/acsch),
/// `σ²x²−1` (acosh/acsc/asec/acoth) — built as `σ²·x²` so it matches the
/// parser's `c^2*x^2` shape — together with the constant term `e0 ∈ {1, −1}`
/// (used to read the normalization `s` off `R(0) = s·e0`).
fn kernel_sum<'a>(
    ctx: &'a AtomArena<'a>,
    f: InvFun,
    sigma: Atom<'a>,
    var: Symbol,
) -> (Atom<'a>, i64) {
    let x = ctx.var(var.as_str());
    let s2x2 = ctx.mul(&[int_pow(ctx, sigma, 2), int_pow(ctx, x, 2)]);
    match f {
        InvFun::Asin | InvFun::Acos | InvFun::Atanh | InvFun::Asech => (
            normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), s2x2])])),
            1,
        ),
        InvFun::Atan | InvFun::Asinh | InvFun::Acsch => {
            (normalize(ctx, ctx.add(&[ctx.num(1), s2x2])), 1)
        }
        InvFun::Acosh | InvFun::Acsc | InvFun::Asec | InvFun::Acoth => {
            (normalize(ctx, ctx.add(&[ctx.num(-1), s2x2])), -1)
        }
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
/// An `s` carrying the imaginary unit is neither (W5b: M4's companion
/// absorption reuses this gate; the rubi-00408 family declines here).
fn sqrt_norm_positive(s_atom: Atom<'_>) -> bool {
    if mentions_imaginary(s_atom) {
        return false;
    }
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
    if Some(ke) != family_kernel_exp(m.f) {
        return None;
    }
    let (sigma, intercept) = linear_form(ctx, m.arg, var)?;
    if !is_zero_atom(ctx, intercept) || is_zero_atom(ctx, sigma) {
        return None;
    }
    let sign = family_sign(m.f)?;
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
    if sign < 0 {
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
        let (coeff, j, na, nb) = scan_hyp_term(ctx, term, t_sym, MAX_SUBST_K)?;
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
        out.push(integrate_hyp_term(ctx, term, t_sym, MAX_SUBST_K)?);
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
/// `C·t^j·cosh(t)^A·sinh(t)^B`; `None` on any other shape. `j_cap` is the
/// caller's cap on the polynomial power (M2: `MAX_SUBST_K`, M4: `MAX_IS_K`).
fn scan_hyp_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    t_sym: Symbol,
    j_cap: i64,
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
    if j > j_cap || na + nb > MAX_HYP_EXP {
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
    j_cap: i64,
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
    if j > j_cap {
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

/// `∫ C·F(a + b·x) dx` for the six main inverse functions, via `u = a + b·x`
/// and the textbook bare forms (rule family E) divided by the slope. The
/// Wave 5c reciprocal heads decline here — M4 owns them and emits the same
/// textbook forms for the same shapes (one owner per shape).
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
        // Wave 5c reciprocal heads: M4 owns them (its substitution emits the
        // branch-consistent `log(w ± √E)` forms for the same bare shapes), so
        // M3 declines here and the two mechanisms never disagree on shape.
        InvFun::Acsc | InvFun::Asec | InvFun::Acoth | InvFun::Asech | InvFun::Acsch => {
            return None;
        }
    };
    let mut out = rest;
    out.push(inv(ctx, slope));
    out.push(core);
    Some(normalize(ctx, ctx.mul(&out)))
}

// =========================================================================
// P1: inverse-composition cancellation
// =========================================================================

/// Real-domain guard under which an `outer(inner(u)) → u` cancellation is an
/// identity. See the table on [`cancel_guard`].
#[derive(Clone, Copy, PartialEq, Eq)]
enum CancelGuard {
    /// The cancellation holds for every real `u` (the inner function may have
    /// its own poles there; the integrand carries them either way).
    Unconditional,
    /// `|u| < π/2`: `atan(tan(u)) = u` only on the principal band.
    AtanBand,
    /// `u ≥ 0`: `acosh(cosh(u)) = |u|`, so it equals `u` only on the
    /// non-negative half-line.
    AcoshHalf,
}

/// The cancellation table, with the real-domain condition each identity needs:
///
/// | pair           | identity                     | condition                        |
/// |----------------|------------------------------|----------------------------------|
/// | `atanh(tanh u)`| `u`                          | none: `tanh` maps ℝ into (−1, 1), where `atanh` inverts it |
/// | `asinh(sinh u)`| `u`                          | none: `sinh: ℝ → ℝ` is a bijection |
/// | `acoth(coth u)`| `u`                          | none on the domain: the only real pole, `u = 0`, is a pole of `coth` itself |
/// | `log(exp u)`   | `u`                          | none: `exp u > 0` for every real `u` |
/// | `atan(tan u)`  | `u`                          | `u ∈ (−π/2, π/2)` — outside it `atan(tan u) = u − kπ` |
/// | `acosh(cosh u)`| `u`                          | `u ≥ 0` — for `u < 0` the value is `−u` |
///
/// The branch conditions (`atan`, `acosh`) can only be certified for a
/// constant rational argument strictly inside the band, so the pre-pass
/// declines them for symbolic arguments instead of guessing a branch; `P2`
/// still covers `atan(tan(u))` exactly (and branch-independently) whenever the
/// integrand is a kernel power times a polynomial cofactor.
fn cancel_guard(outer: &str, inner: &str) -> Option<CancelGuard> {
    Some(match (outer, inner) {
        ("atanh", "tanh") | ("asinh", "sinh") | ("acoth", "coth") | ("log", "exp") => {
            CancelGuard::Unconditional
        }
        ("atan", "tan") => CancelGuard::AtanBand,
        ("acosh", "cosh") => CancelGuard::AcoshHalf,
        _ => return None,
    })
}

/// Is the guard *certified* for the argument `u`? A symbolic `u` (one that
/// contains `var`) can leave the principal band, so only a rational constant
/// inside the band passes — the honest answer when nothing more is known.
fn guard_certified(g: CancelGuard, u: Atom<'_>, var: Symbol) -> bool {
    match g {
        CancelGuard::Unconditional => true,
        CancelGuard::AtanBand => {
            // `|p/q| ≤ 3/2 < π/2` (3/2 is a safe rational lower bound).
            is_constant(u, var)
                && rat_of(u).is_some_and(|(p, q)| {
                    q > 0 && p.saturating_abs().saturating_mul(2) <= q.saturating_mul(3)
                })
        }
        CancelGuard::AcoshHalf => {
            is_constant(u, var) && rat_of(u).is_some_and(|(p, q)| q > 0 && p >= 0)
        }
    }
}

/// Rewrite every child, reporting whether any of them changed.
fn rewrite_children<'a>(
    ctx: &'a AtomArena<'a>,
    args: &[Atom<'a>],
    var: Symbol,
    budget: &mut usize,
) -> (bool, Vec<Atom<'a>>) {
    let mut changed = false;
    let mut out: Vec<Atom<'a>> = Vec::with_capacity(args.len());
    for a in args {
        match cancel_compositions(ctx, *a, var, budget) {
            Some(r) => {
                changed = true;
                out.push(r);
            }
            None => out.push(*a),
        }
    }
    (changed, out)
}

/// Bottom-up rewrite replacing every cancellable composition by its argument.
///
/// Returns `Some(_)` only when at least one replacement happened, so an
/// out-of-family integrand costs one O(nodes) scan and nothing else. Children
/// are rewritten before their parent is re-matched, so one pass reaches the
/// fixpoint for arbitrarily nested pairs; because a rewritten node contains no
/// cancellable pair by construction, a second application is a no-op
/// (idempotent). `budget` is a deterministic node counter — when it runs out
/// the whole pre-pass declines.
fn cancel_compositions<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    budget: &mut usize,
) -> Option<Atom<'a>> {
    if *budget == 0 {
        return None;
    }
    *budget -= 1;
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => None,
        AtomNode::Add(args) => {
            let (changed, out) = rewrite_children(ctx, args, var, budget);
            if changed { Some(ctx.add(&out)) } else { None }
        }
        AtomNode::Mul(args) => {
            let (changed, out) = rewrite_children(ctx, args, var, budget);
            if changed { Some(ctx.mul(&out)) } else { None }
        }
        AtomNode::Pow(b, e) => {
            let nb = cancel_compositions(ctx, *b, var, budget);
            let ne = cancel_compositions(ctx, *e, var, budget);
            match (nb, ne) {
                (None, None) => None,
                (nb, ne) => Some(ctx.pow(nb.unwrap_or(*b), ne.unwrap_or(*e))),
            }
        }
        AtomNode::Fun(name, args) => {
            let (changed, new_args) = rewrite_children(ctx, args, var, budget);
            if new_args.len() == 1
                && let AtomNode::Fun(inner, iargs) = new_args[0].node()
                && iargs.len() == 1
                && let Some(g) = cancel_guard(name.as_str(), inner.as_str())
                && guard_certified(g, iargs[0], var)
            {
                return Some(iargs[0]);
            }
            if changed {
                Some(ctx.fun(name.as_str(), &new_args))
            } else {
                None
            }
        }
    }
}

/// Rational `var`-power of one factor: `x → 1/1`, `x^n → n/1`,
/// `sqrt(x) → 1/2`, `sqrt(x)^n → n/2`, and recursively for powers of those.
/// `None` for anything that is not a `var` power (`normalize` spells `sqrt(x)`
/// as a function head, which is why the plain `Pow` check is not enough).
fn factor_var_power(factor: Atom<'_>, var: Symbol) -> Option<(i64, i64)> {
    match factor.node() {
        AtomNode::Var(v) if *v == var => Some((1, 1)),
        AtomNode::Pow(b, e) => {
            let (p, q) = rat_of(*e)?;
            if q <= 0 {
                return None;
            }
            let (bp, bq) = factor_var_power(*b, var)?;
            Some((bp.checked_mul(p)?, bq.checked_mul(q)?))
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            let (p, q) = factor_var_power(args[0], var)?;
            Some((p, q.checked_mul(2)?))
        }
        _ => None,
    }
}

/// Accumulate `add` into the running `var` power, or `None` when the common
/// denominator would grow past the module's small budget (the caller then
/// hands the term to the ordinary pipeline instead).
fn add_var_power(acc: &mut Option<(i64, i64)>, add: (i64, i64)) -> Option<()> {
    match acc {
        None => {
            *acc = Some(add);
            Some(())
        }
        Some((p, q)) => {
            let nq = q.checked_mul(add.1)?;
            if nq <= 0 || nq > MAX_RAT_Q * MAX_RAT_Q {
                return None;
            }
            let np = p.checked_mul(add.1)?.checked_add(add.0.checked_mul(*q)?)?;
            *acc = Some((np, nq));
            Some(())
        }
    }
}

/// Split a term `C·x^q` into its constant factors and the rational exponent
/// `q = p/q` (`q > 0`); `None` when a non-constant factor is not a `var`
/// power.
fn split_monomial<'a>(term: Atom<'a>, var: Symbol) -> Option<(Vec<Atom<'a>>, i64, i64)> {
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut power: Option<(i64, i64)> = None;
    for f in factors {
        if let Some(pq) = factor_var_power(f, var) {
            add_var_power(&mut power, pq)?;
            continue;
        }
        if is_constant(f, var) {
            coeff.push(f);
            continue;
        }
        return None;
    }
    let (p, q) = power?;
    if q <= 0 {
        return None;
    }
    Some((coeff, p, q))
}

/// Integrate an already expanded and collected expression term by term.
///
/// The cheap paths matter for cost, not correctness: `integrate_raw` on a
/// constant *symbolic* term (expanded shapes leave terms like `a²` behind), or
/// on a `sqrt(x)`-style factor, falls through `integrate_power` /
/// `integrate_function` into `try_risch_or_fallback`, spending entries of the
/// pipeline's global chain budget — while `C·x` and the monomial power rule
/// are exact and immediate. Keeping the re-entry budget-free is what makes the
/// mechanisms safe to repeat and lets repeated calls stay deterministic.
fn integrate_expanded_terms<'a>(
    ctx: &'a AtomArena<'a>,
    folded: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let terms: Vec<Atom<'a>> = match folded.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![folded],
    };
    let x = ctx.var(var.as_str());
    let mut out: Vec<Atom<'a>> = Vec::with_capacity(terms.len());
    for term in terms {
        if is_constant(term, var) {
            out.push(ctx.mul(&[term, x]));
            continue;
        }
        if let Some((coeff, p, q)) = split_monomial(term, var) {
            let shifted = p.checked_add(q)?;
            let mut fs = coeff;
            if shifted == 0 {
                fs.push(ctx.fun("log", &[x]));
            } else {
                let r = rat_atom(ctx, shifted, q);
                fs.push(ctx.pow(x, r));
                fs.push(inv(ctx, r));
            }
            out.push(normalize(ctx, ctx.mul(&fs)));
            continue;
        }
        let r = integrate_raw(ctx, term, var, 0, true, 0, 0);
        if contains_integral(r) {
            return None;
        }
        out.push(r);
    }
    Some(normalize(ctx, ctx.add(&out)))
}

// =========================================================================
// P2: branch-shifted composition kernels with constant slope
// =========================================================================

/// Composition pairs `F∘G` whose derivative with respect to the argument is
/// identically 1, so `d/dx F(G(u)) = σ` for a linear `u = ρ + σ·x`.
///
/// `acoth(tanh(u))` and `atan(tan(u))` are the interesting members: their
/// values are `u ± iπ/2` and `u − kπ` (band-dependent), so they are *not*
/// cancellable to `u` for powers — but their derivatives are exactly the
/// constant slope, which is all the kernel algebra below needs.
fn affine_kernel_pair(outer: &str, inner: &str) -> bool {
    matches!(
        (outer, inner),
        ("atanh", "tanh")
            | ("asinh", "sinh")
            | ("acoth", "coth")
            | ("acoth", "tanh")
            | ("atan", "tan")
            | ("log", "exp")
    )
}

/// A matched composition kernel `K = F(G(u))` with `u = ρ + σ·x`.
struct AffineKernel<'a> {
    /// The `F(G(u))` atom itself, reused verbatim in the answer.
    base: Atom<'a>,
    /// The constant slope `σ` of `u`.
    sigma: Atom<'a>,
}

/// Match a bare `F(G(u))` kernel against [`affine_kernel_pair`].
fn match_affine_kernel<'a>(
    ctx: &'a AtomArena<'a>,
    factor: Atom<'a>,
    var: Symbol,
) -> Option<AffineKernel<'a>> {
    let AtomNode::Fun(name, args) = factor.node() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let AtomNode::Fun(inner, iargs) = args[0].node() else {
        return None;
    };
    if iargs.len() != 1 || !affine_kernel_pair(name.as_str(), inner.as_str()) {
        return None;
    }
    let (sigma, _rho) = linear_form(ctx, iargs[0], var)?;
    if is_zero_atom(ctx, sigma) {
        return None;
    }
    Some(AffineKernel {
        base: factor,
        sigma,
    })
}

/// `var^(−j)` for `1 ≤ j ≤ MAX_RECIP_EXP` — the `1/x^3`, `1/x^6` cofactors of
/// the `acoth(tanh(a + b·x))^3/x^k` corpus shapes.
fn reciprocal_power(factor: Atom<'_>, var: Symbol) -> Option<i64> {
    let AtomNode::Pow(b, e) = factor.node() else {
        return None;
    };
    if !matches!(b.node(), AtomNode::Var(v) if *v == var) {
        return None;
    }
    let AtomNode::Num(j) = e.node() else {
        return None;
    };
    if *j <= -1 && *j >= -MAX_RECIP_EXP {
        Some(-*j)
    } else {
        None
    }
}

/// A fresh constant symbol for the affine proxy, or `None` when all canonical
/// candidates already occur in the integrand.
fn pick_proxy_symbol(expr: Atom<'_>, var: Symbol) -> Option<Symbol> {
    for name in ["w", "k", "r", "z"] {
        let s = Symbol::new(name);
        if s != var && !contains_symbol(expr, s) {
            return Some(s);
        }
    }
    None
}

/// `∫ K^m · P(x) · x^(−j) dx` for a branch-shifted composition kernel `K`
/// with `K' = σ` (see [`affine_kernel_pair`]).
///
/// The substitution `K = σ·x + w` turns the integrand into an elementary
/// expression in `x` and the fresh constant `w`; after integrating in `x` the
/// map `w ↦ K − σ·x` restores the kernel. Because `K' = σ` is constant, the
/// chain rule gives `d/dx F(x, K − σx) = (∂_x F + ∂_w F·(K' − σ))` evaluated
/// at `w = K − σx`, i.e. exactly the original integrand.
///
/// Shape family (everything outside it declines): exactly one kernel power
/// `K^m` (`|m| ≤ MAX_KERNEL_POW`), polynomial cofactors of total degree
/// `≤ MAX_COFACTOR_DEG`, at most one `x^(−j)` (`j ≤ MAX_RECIP_EXP`) and
/// arbitrary constants.
pub(crate) fn integrate_affine_kernel<'a>(
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
    let mut kern: Option<(AffineKernel<'a>, i64)> = None;
    let mut cofactor: Vec<Atom<'a>> = Vec::new();
    let mut cofactor_deg = 0usize;
    let mut recip = 0i64;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if kern.is_none() {
            if let Some(k) = match_affine_kernel(ctx, f, var) {
                kern = Some((k, 1));
                continue;
            }
            if let AtomNode::Pow(b, e) = f.node()
                && let AtomNode::Num(m) = e.node()
                && *m != 0
                && (-MAX_KERNEL_POW..=MAX_KERNEL_POW).contains(m)
                && let Some(k) = match_affine_kernel(ctx, *b, var)
            {
                kern = Some((k, *m));
                continue;
            }
        }
        if recip == 0
            && let Some(j) = reciprocal_power(f, var)
        {
            recip = j;
            cofactor.push(f);
            continue;
        }
        if let Some(d) = poly_deg(f, var)
            && d >= 1
            && cofactor_deg + d <= MAX_COFACTOR_DEG
        {
            cofactor_deg += d;
            cofactor.push(f);
            continue;
        }
        return None;
    }
    let (k, m) = kern?;
    let w_sym = pick_proxy_symbol(expr, var)?;
    let x = ctx.var(var.as_str());
    let w = ctx.var(w_sym.as_str());
    let lin = normalize(ctx, ctx.add(&[ctx.mul(&[k.sigma, x]), w]));
    let mut rewritten: Vec<Atom<'a>> = rest;
    rewritten.push(int_pow(ctx, lin, m));
    rewritten.extend(cofactor);
    let product = normalize(ctx, ctx.mul(&rewritten));
    if node_count(product) > MAX_KERNEL_NODES {
        return None;
    }
    let expanded = crate::expand::expand_bounded(ctx, product).unwrap_or(product);
    let folded = crate::ode::util::collect_terms(ctx, normalize(ctx, expanded));
    if node_count(folded) > MAX_KERNEL_NODES {
        return None;
    }
    let chain = integrate_expanded_terms(ctx, folded, var)?;
    if contains_integral(chain) {
        return None;
    }
    let back = normalize(ctx, ctx.add(&[k.base, ctx.mul(&[ctx.num(-1), k.sigma, x])]));
    let res = normalize(ctx, replace_symbol(ctx, chain, w_sym, back));
    if contains_symbol(res, w_sym) {
        return None;
    }
    Some(res)
}

// =========================================================================
// P3: integral representation — companion radicands of the inverse families
// =========================================================================

/// The kernel sum `E(u)` expressed in the *argument* `u` itself, so the
/// companion identity can be checked for arguments that carry an intercept
/// (M1 reads the normalization off `E(σ·x)`, which forces `ρ = 0`).
///
/// The Wave 5c reciprocal heads use the same three radicands as their main
/// counterparts — `u²−1` for acsc/asec/acoth, `1−u²` for asech, `1+u²` for
/// acsch. For `acoth` this is the *sign-corrected* choice: the spec table's
/// `1−σ²x²` is negative on `acoth`'s whole real domain (`|w| > 1`), and its
/// principal square root would carry `i` — `w²−1 = csch²u` is the radicand
/// whose square root is the family's own forward head (the no-`I` discipline
/// of spec §1.3-7).
fn family_kernel_of_u<'a>(ctx: &'a AtomArena<'a>, f: InvFun, u: Atom<'a>) -> Atom<'a> {
    let u2 = int_pow(ctx, u, 2);
    match f {
        InvFun::Asin | InvFun::Acos | InvFun::Atanh | InvFun::Asech => {
            normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), u2])]))
        }
        InvFun::Atan | InvFun::Asinh | InvFun::Acsch => normalize(ctx, ctx.add(&[ctx.num(1), u2])),
        InvFun::Acosh | InvFun::Acsc | InvFun::Asec | InvFun::Acoth => {
            normalize(ctx, ctx.add(&[ctx.num(-1), u2]))
        }
    }
}

/// Match `(a + b·F(arg))^k` for a rational `k = p/q` with `1 ≤ p ≤ 6q`,
/// `q ≤ 4` (also the `sqrt(F(arg))` spelling as `k = 1/2`); returns
/// `(F, arg, a, b, p, q)`. Negative and zero exponents stay out of scope,
/// matching M1's positive-integer family.
fn match_inv_power_rational<'a>(
    ctx: &'a AtomArena<'a>,
    factor: Atom<'a>,
    var: Symbol,
) -> Option<(InvFun, Atom<'a>, Atom<'a>, Atom<'a>, i64, i64)> {
    if let AtomNode::Fun(name, args) = factor.node()
        && name.as_str() == "sqrt"
        && args.len() == 1
        && let Some((f, arg, a, b)) = match_inv_base(ctx, args[0], var)
    {
        return Some((f, arg, a, b, 1, 2));
    }
    if let AtomNode::Pow(b, e) = factor.node()
        && let Some((p, q)) = rat_of(*e)
        && q > 0
        && q <= MAX_RAT_Q
        && p >= 1
        && p <= MAX_POWER_K * q
        && let Some((f, arg, a, bb)) = match_inv_base(ctx, *b, var)
    {
        return Some((f, arg, a, bb, p, q));
    }
    let (f, arg, a, b) = match_inv_base(ctx, factor, var)?;
    Some((f, arg, a, b, 1, 1))
}

/// Distribute integer powers over products: `(u₁·u₂·…)^n → u₁^n·u₂^n·…`.
///
/// `normalize` keeps `Pow(Mul([a, x]), 2)` intact and `collect_terms` cannot
/// see through it, so a radicand like `1 − (a·x)²` would never fold against
/// its expanded twin without this pass. Deterministic and budget-bounded.
fn flatten_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, budget: &mut usize) -> Atom<'a> {
    if *budget == 0 {
        return expr;
    }
    *budget -= 1;
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| flatten_powers(ctx, *a, budget))
                .collect();
            ctx.add(&out)
        }
        AtomNode::Mul(args) => {
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| flatten_powers(ctx, *a, budget))
                .collect();
            ctx.mul(&out)
        }
        AtomNode::Pow(b, e) => {
            let nb = flatten_powers(ctx, *b, budget);
            let ne = flatten_powers(ctx, *e, budget);
            if let AtomNode::Num(n) = ne.node()
                && *n >= 2
                && let AtomNode::Mul(fs) = nb.node()
            {
                let out: Vec<Atom<'a>> = fs
                    .iter()
                    .map(|f| ctx.pow(flatten_powers(ctx, *f, budget), ne))
                    .collect();
                return ctx.mul(&out);
            }
            ctx.pow(nb, ne)
        }
        AtomNode::Fun(name, args) => {
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| flatten_powers(ctx, *a, budget))
                .collect();
            ctx.fun(name.as_str(), &out)
        }
    }
}

/// The canonical expanded form the companion identity is checked against:
/// bounded expansion, product-power flattening, then like-term collection.
fn companion_normal_form<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    let expanded = crate::expand::expand_bounded(ctx, expr).unwrap_or(expr);
    let mut budget = MAX_CANCEL_BUDGET;
    let flat = flatten_powers(ctx, expanded, &mut budget);
    crate::ode::util::collect_terms(ctx, normalize(ctx, flat))
}

/// Scale `s` with `R ≡ s·E(u)`, confirmed by a bounded expansion fold.
///
/// The candidate is *read off* at an abscissa where `E(u)` is a number, which
/// keeps the extraction division-free in the symbolic parameters: the root
/// `x = −ρ/σ` of the argument makes `E(u) = e0 ∈ {±1}`, and `x = 0` does the
/// same for intercept-free arguments. Reading anywhere else would leave a
/// `P(x)/Q(x)` ratio that `normalize` cannot cancel (it keeps `A·A⁻¹` for
/// compound `A`), so those candidates are skipped rather than trusted.
fn radicand_scale<'a>(
    ctx: &'a AtomArena<'a>,
    r: Atom<'a>,
    e_u: Atom<'a>,
    x_root: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let re = companion_normal_form(ctx, r);
    let ee = companion_normal_form(ctx, e_u);
    for point in [x_root, ctx.num(0)] {
        let rp = companion_normal_form(ctx, replace_symbol(ctx, re, var, point));
        let ep = companion_normal_form(ctx, replace_symbol(ctx, ee, var, point));
        if !is_constant(rp, var) || is_zero_atom(ctx, ep) {
            continue;
        }
        let s = normalize(ctx, ctx.mul(&[rp, inv(ctx, ep)]));
        if is_zero_atom(ctx, s) {
            continue;
        }
        let diff = ctx.add(&[re, ctx.mul(&[ctx.num(-1), s, ee])]);
        let folded = companion_normal_form(ctx, diff);
        if folds_to_zero(ctx, folded) {
            return Some(s);
        }
    }
    None
}

/// `∫ C·(a + b·F(u))^k · R^e dx` where `R` *is* the derivative-shaped
/// companion `s·E(u)` of `F` at the same argument `u` (with `e` the family's
/// derivative exponent), by the chain rule
/// `d/dx (a + b·F(u))^(k+1) = b·σ·(k+1)·(a + b·F(u))^k·s_F·E(u)^e`.
///
/// The companion check is the whole point: `atan(a + b·x)/sqrt(1 + (a+b·x)²)`
/// has a `−1/2` radical but `atan`'s derivative exponent is `−1`, so its
/// radical is *not* a derivative companion and the antiderivative is not
/// elementary (it needs the inverse tangent integral `Ti₂`); the pass detects
/// the mismatch and declines. `k` may be rational (`sqrt(asin(a·x))` is
/// `k = 1/2`), and the argument may carry an intercept, which M1 rejects.
pub(crate) fn integrate_companion_form<'a>(
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
    let mut inv_f: Option<(InvFun, Atom<'a>, Atom<'a>, Atom<'a>, i64, i64)> = None;
    let mut kern: Option<(Atom<'a>, KernelExp)> = None;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if inv_f.is_none()
            && let Some(m) = match_inv_power_rational(ctx, f, var)
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
    let (f, arg, a, b, p, q) = inv_f?;
    let (r, ke) = kern?;
    if Some(ke) != family_kernel_exp(f) {
        return None;
    }
    let (sigma, rho) = linear_form(ctx, arg, var)?;
    if is_zero_atom(ctx, sigma) {
        return None;
    }
    let e_u = family_kernel_of_u(ctx, f, arg);
    // `x = −ρ/σ`, the abscissa where the argument vanishes and `E(u) = e0`.
    let x_root = normalize(ctx, ctx.mul(&[ctx.num(-1), rho, inv(ctx, sigma)]));
    let s_atom = radicand_scale(ctx, r, e_u, x_root, var)?;
    if ke == KernelExp::InvSqrt && !sqrt_norm_positive(s_atom) {
        return None;
    }
    let s_factor = match ke {
        KernelExp::Recip => inv(ctx, s_atom),
        KernelExp::InvSqrt => ctx.pow(s_atom, rat_atom(ctx, -1, 2)),
    };
    let base = normalize(ctx, ctx.add(&[a, ctx.mul(&[b, ctx.fun(f.name(), &[arg])])]));
    let mut out = rest;
    if family_sign(f)? < 0 {
        out.push(ctx.num(-1));
    }
    out.push(s_factor);
    out.push(ctx.pow(base, rat_atom(ctx, p + q, q)));
    let denom = normalize(ctx, ctx.mul(&[b, sigma, rat_atom(ctx, p + q, q)]));
    out.push(inv(ctx, denom));
    Some(normalize(ctx, ctx.mul(&out)))
}

// =========================================================================
// M4: inverse-function substitution (0.30.0 Waves 5a/5b/5c)
// =========================================================================
//
// Implementation of `docs/planning/INV_SUBST_SPEC_CN.md` §2.1/2.3/2.4/2.6:
// `∫ C·(a + b·F(w))^k · P(x) dx` with `F` one of the eleven inverse heads
// (the six main heads of W5a/W5b plus the five reciprocal heads of W5c),
// `w = ρ + σ·x` affine (`σ ≢ 0`), integer `1 ≤ k ≤ MAX_IS_K` and
// `P` polynomial of total degree ≤ `MAX_IS_POLY_DEG`. The substitution
// `u = F(w)` maps the integrand to `P(u)·T(u)` with `T` a trig/hyperbolic
// monomial (`x = (G(u) − ρ)/σ`, `dx = G′(u)/σ·du`, per the spec §2.1
// table); elementarity is **not** prejudged — a re-integration residue
// declines the whole shape (§2.5).
//
// W5a/W5b/W5c scope limits (deviations from the full spec are all narrowings
// or measured reorderings):
//
// - **Companion radicals** (W5b, §2.2): at most one radical factor
//   `(s·E(w))^(m/2)` with integer `0 < |m| ≤ MAX_IS_RAD_M` is accepted, and
//   only when the radicand is the family's own kernel at the same argument
//   (`radicand_scale` reuse, the P3 mechanism). It absorbs as
//   `s^(m/2)·G′(u)^m·σ`, merged with the `dx` factor's `G′(u)` power into a
//   single head power. Non-companion radicands, a second radical base, and
//   imaginary-carrying radicands all decline (§2.5-4). Narrowings: a
//   *negative* combined head power is only integrated where the existing
//   tabulars see it — the ratio families convert `R^A·S^B` (`B < 0`) to the
//   exact `sin^A·cos^{−A−B}` / `sinh^A·cosh^{−A−B}` form and re-use the
//   trig/hyperbolic machinery; the trig/hyperbolic families (asin/acos/
//   asinh/acosh) decline negative combined powers (e.g. the elementary but
//   unimplemented `u·sin u·sec²u` family).
// - **Sheet factor σ emission record (§4.2)**: for odd `m` the absorption
//   carries `σ = G′(u)·E(G(u))^(−1/2)`, a per-interval `±1` constant of zero
//   formal derivative. On the principal band of every one of the six
//   families the unsigned `G′` head is non-negative (asin `cos u ≥ 0` on
//   `[−π/2, π/2]`, acos `sin u ≥ 0` on `[0, π]`, atan `sec² u > 0`, asinh
//   `cosh u > 0`, acosh `sinh u ≥ 0` on `[0, ∞)`, atanh `sech²u > 0`), so
//   the sign evidence is definite and σ folds to `+1`; the emitted answer
//   is exact on that band (the same standing semantics as the
//   back-substitution's `G′(u) → +√E(w)`). No σ node is emitted.
// - **The t-form re-integration is fully in-module, no chain re-entry**
//   (§2.3 allows `integrate_raw` re-entry; the tabular form is its
//   sanctioned T3-analogue). Measured reason: the chain-entry budget
//   (`MAX_CHAIN_ENTRIES = 256`) is shared per top-level call and is
//   typically spent before the `inverse_trig` stage runs — traced on
//   `x^2*acos(a*x)^2` (rubi-00554), where a chain-based t-form
//   re-integration died on the budget and the pipeline fell back. The
//   tabular paths are also what keeps `exp(μu)`-shaped chain answers (the
//   hyperbolic Weierstrass routes) out of the whitelist back-substitution:
//   - trig families (asin/acos): `trig_reduce_products` multiple-angle
//     reduction, then `integrate_trig_term` (repeated-parts table for
//     `C·u^j·sin/cos(μu)`);
//   - hyperbolic families (asinh/acosh): M2's `scan_hyp_term` /
//     `reduce_hyp_term` expansion plus `integrate_hyp_term` (with the
//     j-cap widened to `MAX_IS_K`);
//   - ratio families (atan/atanh **and the Wave 5c reciprocals**):
//     `ratio_i`/`ratio_j`, parts on the exact differential `dv = R^a·S² du`
//     plus the `R^c → R^{c−2}` peel (`tan² = sec² − 1` / `tanh² = 1 − sech²`
//     / `cot² = csc² − 1` / `coth² = 1 + csch²`); the `J(j≥1, 1)` base
//     case is the polylog territory of §2.5 and declines.
// - **tan/sec even-power pre-reduction** (§2.3 second clause) is
//   subsumed by the ratio tabular, which handles `tan^A·sec²` directly.
//
// ## Wave 5c: the five reciprocal heads (0.30.0)
//
// `acsc`/`asec`/`acoth`/`asech`/`acsch` run through the same pipeline with
// the spec §2.1 reciprocal rows: `acsc → G = csc`, `asec → sec`,
// `acoth → coth`, `asech → sech`, `acsch → csch`, and
// `G′` = `−csc·cot`, `sec·tan`, `−csch²`, `−sech·tanh`, `−csch·coth`
// (the only families whose `G′` is a *product* of two heads — the `dx`
// factor therefore contributes two head factors instead of one).
//
// The t-form `u^j·P(G(u))·G′(u)/σ` is again a `C·u^j·R^A·S^B` term, but in
// the reciprocal pair: `(cot, csc)` for acsc, `(tan, sec)` for asec,
// `(coth, csch)` for acoth/acsch, `(tanh, sech)` for asech. [`RatioFam`]
// generalises the pre-existing `(tan, sec)`/`(tanh, sech)` tabular to all
// four through the three signs `R′ = ε·S²`, `S′ = δ·R·S`,
// `S² = η·(1 + γ·R²)`; the two new families are `(ε, δ) = (−1, −1)` with
// `(η, γ) = (1, 1)` (cot/csc) and `(−1, −1)` (coth/csch).
//
// Kernel radicands (spec §2.1, with one documented deviation): `w² − 1` for
// acsc/asec/acoth, `1 − w²` for asech, `1 + w²` for acsch. The deviation is
// `acoth`: the spec table lists `1 − σ²x²`, which is *negative* on
// `acoth`'s whole real domain (`|w| > 1`) and whose principal square root
// would carry `i`; `w² − 1 = csch²u` is the sign-corrected radicand whose
// square root is the family's own forward `csch` (the no-`I` discipline of
// §1.3-7). Companion radicals absorb exactly as in W5b, into the *kernel
// head* of the family (`cot`/`tan`/`csch`/`tanh`/`coth` — see
// [`family_kernel_head`]).
//
// Ownership boundaries: the reciprocal heads are M4's alone. M1/P3 decline
// them ([`family_kernel_exp`] returns `None`) because their `F′` carries an
// extra reciprocal forward head (`d/dw acsc(w) = −1/(|w|√(w²−1))`) that the
// kernel-derivative power rule does not model, and M3 declines them so the
// bare `∫F(w) dw` shape has a single owner (M4 emits the same textbook
// `w·F(w) ± log(w + √E)` / `w·F(w) ± acos(w)` forms).
//
// Capability record (measured, the numeric tabular self-check in the tests
// covers the whole `(j, A, B)` grid):
//
// - the bare and affine shapes, every `k = 1` kernel power with a
//   polynomial weight, `k = 2` with a linear weight, and the
//   companion-radical shapes all close;
// - `k ≥ 3` kernel powers, bare `(a + b·F)^k` with `k ≥ 2`, and any term
//   whose reduction reaches `∫u^j·S(u) du` with `j ≥ 1` decline (polylog
//   territory, §2.5) — including rubi-01132/01218/01347/01467 and
//   rubi-00372;
// - symbolic exponents, negative powers, non-affine arguments and
//   non-companion radicals decline through the unchanged W5a/W5b gates;
// - the two bucket rows inside the caps, rubi-00511 (`x³·acsc(a+bx)`) and
//   rubi-01708 (`x⁴·asec(a+bx)`), solve and certify.
//
// Re-entry (§1.3-5): the answers contain `log(w + √E(w))` and `acos(w)`
// pieces that the rest of the pipeline still declines to integrate
// (`∫log(w + √E) dw`, `∫√((a+bx)²−1) dx`), so — as for the W5a T3 case —
// the idempotence test asserts termination plus completeness of M4's *own*
// re-claims rather than a residue-free full-chain re-entry.
//
// Real-branch semantics (spec §4): every answer is exact on the principal
// band of `F` (asin: `u ∈ [−π/2, π/2]`; acos: `[0, π]`; atan:
// `(−π/2, π/2)`; asinh/atanh: ℝ; acosh: `[0, ∞)`; Wave 5c: `u > 0` for
// acsc/asec (where `cot u`/`tan u` is the positive `√E(w)`), `u > 0` for
// asech/acsch (`tanh u`/`coth u` positive) and `u > 0` for acoth
// (`csch u > 0`); the `w < 0` halves are deliberately out of scope, the same
// real-value discipline as acosh's `w < −1` half, §4.1). On the band `G′(u)`
// has a definite sign for all eleven families — the
// back-substitution (`cos u → +√E(w)`) already assumes it. The radical
// absorption extends the same convention: for odd `m`,
// `(s·E)^(m/2) = s^(m/2)·G′^m·σ` folds the sheet factor σ to `+1` with
// the principal-band sign evidence (§4.2), and the `√E(w)` answers
// re-enter the certificate's radical layer through the `r² = E(w)`
// relation.
// The in-module tabulars never produce `F(G(u))` composites; the terminal
// structural check (`contains_inv_composite`) is the belt-and-braces
// enforcement of §2.4's ban. No `I` is introduced anywhere.
//
// Wave 5c's self-check: the parameterised ratio tabular is verified against
// numerical quadrature over the whole `(j, A, B)` grid of all four families
// (`m4_w5c_ratio_tabular_numeric`), and every W5c answer is additionally
// checked by numeric differentiation inside its own real domain
// (`m4_w5c_*` tests). That self-check is what exposed a latent wrong closed
// form in `ratio_s_pow` for `j ≥ 1, b ≥ 3` (see its docs) — the fix declines
// those non-elementary reductions.

/// M4 cap on the `(a + b·F(w))^k` exponent (spec §3.2, aligned with
/// `MAX_POWER_K`).
const MAX_IS_K: i64 = 6;
/// M4 cap on the total degree of the polynomial weight (spec §3.2).
const MAX_IS_POLY_DEG: usize = 6;
/// M4 node budget for the substituted t-form (spec §3.2).
const MAX_IS_SUBST_NODES: usize = 400;
/// M4 cap on the multiple `μ` of a reduced trig/hyperbolic angle and on the
/// Chebyshev back-substitution order (spec §3.2 `MAX_IS_TRIG_EXP`; aligned
/// with `MAX_HYP_EXP` and trig_reduce's factor cap).
const MAX_IS_MU: i64 = 8;
/// M4 cap on the companion-radical half-power: `(s·E(w))^(m/2)` is accepted
/// for integer `0 < |m| ≤ MAX_IS_RAD_M` (spec §3.2, W5b).
const MAX_IS_RAD_M: i64 = 5;
/// M4 total work budget (spec §3.2): node counts are charged at every
/// stage so a single call's work is deterministically bounded.
const MAX_IS_WORK: u64 = 100_000;

/// Charge `nodes` against the M4 work budget; `true` when exhausted.
fn is_work_exceeded(work: &mut u64, nodes: usize) -> bool {
    *work = work.saturating_add(nodes as u64);
    *work > MAX_IS_WORK
}

/// `G` head of the family (`x = G(u)` at zero intercept), per spec §2.1 —
/// including the Wave 5c reciprocal rows (`acsc → csc`, `asec → sec`,
/// `acoth → coth`, `asech → sech`, `acsch → csch`).
fn family_g_head(f: InvFun) -> &'static str {
    match f {
        InvFun::Asin => "sin",
        InvFun::Acos => "cos",
        InvFun::Atan => "tan",
        InvFun::Asinh => "sinh",
        InvFun::Acosh => "cosh",
        InvFun::Atanh => "tanh",
        InvFun::Acsc => "csc",
        InvFun::Asec => "sec",
        InvFun::Acoth => "coth",
        InvFun::Asech => "sech",
        InvFun::Acsch => "csch",
    }
}

/// `G′(u)` per family as `(sign, head factors)`: asin `cos`, acos `−sin`,
/// atan `sec²`, asinh `cosh`, acosh `sinh`, atanh `sech²`; Wave 5c:
/// acsc `−csc·cot`, asec `sec·tan`, acoth `−csch²`, asech `−sech·tanh`,
/// acsch `−csch·coth` (spec §2.1 table).
fn family_gp_terms(f: InvFun) -> (i64, &'static [(&'static str, i64)]) {
    match f {
        InvFun::Asin => (1, &[("cos", 1)]),
        InvFun::Acos => (-1, &[("sin", 1)]),
        InvFun::Atan => (1, &[("sec", 2)]),
        InvFun::Asinh => (1, &[("cosh", 1)]),
        InvFun::Acosh => (1, &[("sinh", 1)]),
        InvFun::Atanh => (1, &[("sech", 2)]),
        InvFun::Acsc => (-1, &[("csc", 1), ("cot", 1)]),
        InvFun::Asec => (1, &[("sec", 1), ("tan", 1)]),
        InvFun::Acoth => (-1, &[("csch", 2)]),
        InvFun::Asech => (-1, &[("sech", 1), ("tanh", 1)]),
        InvFun::Acsch => (-1, &[("csch", 1), ("coth", 1)]),
    }
}

/// The forward head whose square is the family kernel at the substitution
/// value: `E(G(u)) = head²` for asin/acos/asinh/acosh and for the Wave 5c
/// reciprocals, and `E(G(u)) = head` for atan/atanh (where `G′` already
/// carries the square). A companion radical `(s·E(w))^(m/2)` therefore
/// contributes `s^(m/2)·head^m` in every family, so it merges into the same
/// head power the `dx` factor's `G′`/σ supplies (spec §2.2/§4.2).
fn family_kernel_head(f: InvFun) -> &'static str {
    match f {
        InvFun::Asin => "cos",
        InvFun::Acos => "sin",
        InvFun::Atan => "sec",
        InvFun::Asinh => "cosh",
        InvFun::Acosh => "sinh",
        InvFun::Atanh => "sech",
        InvFun::Acsc => "cot",
        InvFun::Asec => "tan",
        InvFun::Acoth => "csch",
        InvFun::Asech => "tanh",
        InvFun::Acsch => "coth",
    }
}

/// Everything the back-substitution needs: the substitution variable, the
/// family, `F(w)`, the argument `w`, `E(w)`, and `√E(w)`.
struct IsBackSub<'a> {
    t_sym: Symbol,
    f: InvFun,
    fw: Atom<'a>,
    w: Atom<'a>,
    e_w: Atom<'a>,
    sqrt_e: Atom<'a>,
}

/// `c^k` in the back-substitution, with the `√E` fold: even powers of
/// `√E(w)` land on `E(w)^(k/2)` and odd powers on `E(w)^((k−1)/2)·√E(w)`,
/// so the answer never carries `sqrt(E)·sqrt(E)` products (they neither
/// normalize nor re-enter the chain cleanly).
fn bs_power<'a>(ctx: &'a AtomArena<'a>, bs: &IsBackSub<'a>, c: Atom<'a>, k: i64) -> Atom<'a> {
    if c == bs.sqrt_e {
        if k % 2 == 0 {
            int_pow(ctx, bs.e_w, k / 2)
        } else {
            normalize(
                ctx,
                ctx.mul(&[int_pow(ctx, bs.e_w, (k - 1) / 2), bs.sqrt_e]),
            )
        }
    } else {
        int_pow(ctx, c, k)
    }
}

/// The principal-band algebraic value of one trig/hyperbolic head at `u`
/// (spec §2.1/§4.2): `G(u) = w`, `G′(u) = +√E(w)` (the sheet factor is
/// `+1` on the band), and the remaining heads of the class by the algebraic
/// identities. `None` for heads of the other class.
///
/// Wave 5c reciprocal families: the *identity* head is the one `G` names
/// (`csc u = w` for acsc, `sech u = w` for asech, …) and `√E(w)` is the
/// family's kernel head (`cot u` for acsc, `tanh u` for asech, …); the two
/// remaining heads of the pair follow from `head·identity = √E`-style
/// identities. Both reciprocal families have `E = w² − 1` for
/// acsc/asec/acoth (`u > 0` band, `√E ≥ 0`) and `E = 1 − w²`/`1 + w²` for
/// asech/acsch.
fn is_head_map<'a>(ctx: &'a AtomArena<'a>, bs: &IsBackSub<'a>, name: &str) -> Option<Atom<'a>> {
    let w = bs.w;
    let r = bs.sqrt_e;
    let w_over_r = |ctx: &'a AtomArena<'a>| normalize(ctx, ctx.mul(&[w, inv(ctx, r)]));
    let r_over_w = |ctx: &'a AtomArena<'a>| normalize(ctx, ctx.mul(&[r, inv(ctx, w)]));
    Some(match bs.f {
        InvFun::Asin => match name {
            "sin" => w,
            "cos" => r,
            "tan" => w_over_r(ctx),
            "sec" => inv(ctx, r),
            "cot" => r_over_w(ctx),
            "csc" => inv(ctx, w),
            _ => return None,
        },
        InvFun::Acos => match name {
            "cos" => w,
            "sin" => r,
            "tan" => r_over_w(ctx),
            "sec" => inv(ctx, w),
            "cot" => w_over_r(ctx),
            "csc" => inv(ctx, r),
            _ => return None,
        },
        InvFun::Atan => match name {
            "tan" => w,
            "sec" => r,
            "sin" => w_over_r(ctx),
            "cos" => inv(ctx, r),
            "cot" => inv(ctx, w),
            "csc" => r_over_w(ctx),
            _ => return None,
        },
        InvFun::Asinh => match name {
            "sinh" => w,
            "cosh" => r,
            "tanh" => w_over_r(ctx),
            "sech" => inv(ctx, r),
            "coth" => r_over_w(ctx),
            "csch" => inv(ctx, w),
            _ => return None,
        },
        InvFun::Acosh => match name {
            "cosh" => w,
            "sinh" => r,
            "tanh" => r_over_w(ctx),
            "sech" => inv(ctx, w),
            "coth" => w_over_r(ctx),
            "csch" => inv(ctx, r),
            _ => return None,
        },
        InvFun::Atanh => match name {
            "tanh" => w,
            "sech" => r,
            "sinh" => w_over_r(ctx),
            "cosh" => inv(ctx, r),
            "coth" => inv(ctx, w),
            "csch" => r_over_w(ctx),
            _ => return None,
        },
        // `w = csc u`, `√E(w) = cot u` (`E = w² − 1`).
        InvFun::Acsc => match name {
            "csc" => w,
            "cot" => r,
            "sin" => inv(ctx, w),
            "cos" => r_over_w(ctx),
            "tan" => inv(ctx, r),
            "sec" => w_over_r(ctx),
            _ => return None,
        },
        // `w = sec u`, `√E(w) = tan u`.
        InvFun::Asec => match name {
            "sec" => w,
            "tan" => r,
            "cos" => inv(ctx, w),
            "sin" => r_over_w(ctx),
            "cot" => inv(ctx, r),
            "csc" => w_over_r(ctx),
            _ => return None,
        },
        // `w = coth u`, `√E(w) = csch u` (`E = w² − 1`, `u > 0` band).
        InvFun::Acoth => match name {
            "coth" => w,
            "csch" => r,
            "sinh" => inv(ctx, r),
            "cosh" => w_over_r(ctx),
            "tanh" => inv(ctx, w),
            "sech" => r_over_w(ctx),
            _ => return None,
        },
        // `w = sech u`, `√E(w) = tanh u` (`E = 1 − w²`).
        InvFun::Asech => match name {
            "sech" => w,
            "tanh" => r,
            "cosh" => inv(ctx, w),
            "sinh" => r_over_w(ctx),
            "coth" => inv(ctx, r),
            "csch" => w_over_r(ctx),
            _ => return None,
        },
        // `w = csch u`, `√E(w) = coth u` (`E = 1 + w²`).
        InvFun::Acsch => match name {
            "csch" => w,
            "coth" => r,
            "sinh" => inv(ctx, w),
            "cosh" => r_over_w(ctx),
            "tanh" => inv(ctx, r),
            "sech" => w_over_r(ctx),
            _ => return None,
        },
    })
}

/// `Some(μ)` when `arg = μ·t` with `μ` a non-negative integer constant
/// (`t` itself → 1); `None` for any other shape.
fn u_multiple(arg: Atom<'_>, t_sym: Symbol) -> Option<i64> {
    match arg.node() {
        AtomNode::Num(n) if *n >= 0 => Some(*n),
        AtomNode::Var(v) if *v == t_sym => Some(1),
        AtomNode::Mul(args) => {
            let mut mu: Option<i64> = None;
            let mut seen_t = false;
            for a in args.iter() {
                match a.node() {
                    AtomNode::Num(n) if *n >= 0 && mu.is_none() => mu = Some(*n),
                    AtomNode::Var(v) if *v == t_sym && !seen_t => seen_t = true,
                    _ => return None,
                }
            }
            if seen_t { Some(mu.unwrap_or(1)) } else { None }
        }
        _ => None,
    }
}

/// Coefficients of `T_n` (low-to-high): `cos(n·u) = T_n(cos u)` (and the
/// same polynomial for cosh). Recurrence `T_{n+1} = 2c·T_n − T_{n−1}`.
fn chebyshev_t_coeffs(n: i64) -> Vec<i64> {
    if n == 0 {
        return vec![1];
    }
    let mut prev = vec![1];
    let mut cur = vec![0, 1];
    for _ in 1..n {
        let mut next = vec![0i64; cur.len() + 1];
        for (k, a) in cur.iter().enumerate() {
            next[k + 1] += 2 * a;
        }
        for (k, a) in prev.iter().enumerate() {
            next[k] -= a;
        }
        prev = cur;
        cur = next;
    }
    cur
}

/// Coefficients of `U_n` (low-to-high): `sin((n+1)·u) = sin u·U_n(cos u)`
/// (and the same polynomial for sinh/cosh).
fn chebyshev_u_coeffs(n: i64) -> Vec<i64> {
    if n == 0 {
        return vec![1];
    }
    let mut prev = vec![1];
    let mut cur = vec![0, 2];
    for _ in 1..n {
        let mut next = vec![0i64; cur.len() + 1];
        for (k, a) in cur.iter().enumerate() {
            next[k + 1] += 2 * a;
        }
        for (k, a) in prev.iter().enumerate() {
            next[k] -= a;
        }
        prev = cur;
        cur = next;
    }
    cur
}

/// Evaluate a coefficient polynomial at `c` with the `√E` power fold.
fn chebyshev_eval<'a>(
    ctx: &'a AtomArena<'a>,
    bs: &IsBackSub<'a>,
    coeffs: &[i64],
    c: Atom<'a>,
) -> Atom<'a> {
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for (k, a) in coeffs.iter().enumerate() {
        if *a == 0 {
            continue;
        }
        terms.push(normalize(
            ctx,
            ctx.mul(&[ctx.num(*a), bs_power(ctx, bs, c, k as i64)]),
        ));
    }
    match terms.len() {
        0 => ctx.num(0),
        1 => terms[0],
        _ => normalize(ctx, ctx.add(&terms)),
    }
}

/// Map one trig/hyperbolic head at `μ·u` to its algebraic value in `w`
/// (spec §2.4). Multiple angles go through Chebyshev polynomials
/// (`cos(μu) = T_μ(cos u)`, `sin(μu) = sin u·U_{μ−1}(cos u)`); the
/// tan/sec/tanh/sech heads only exist at `μ = 1` in the accepted fragment.
fn is_trig_back_map<'a>(
    ctx: &'a AtomArena<'a>,
    bs: &IsBackSub<'a>,
    name: &str,
    mu: i64,
) -> Option<Atom<'a>> {
    if mu > MAX_IS_MU {
        return None;
    }
    if mu == 0 {
        return Some(match name {
            "cos" | "cosh" | "sec" | "sech" => ctx.num(1),
            _ => ctx.num(0),
        });
    }
    if mu == 1 {
        return is_head_map(ctx, bs, name);
    }
    let (s_name, c_name) = match name {
        "sin" | "cos" => ("sin", "cos"),
        "sinh" | "cosh" => ("sinh", "cosh"),
        _ => return None,
    };
    let c = is_head_map(ctx, bs, c_name)?;
    if matches!(name, "cos" | "cosh") {
        return Some(chebyshev_eval(ctx, bs, &chebyshev_t_coeffs(mu), c));
    }
    // sin(μu) = sin u·U_{μ−1}(cos u): in every family the two head maps
    // land on different atoms (`w` and `√E`), so the plain product needs
    // no further folding.
    let s = is_head_map(ctx, bs, s_name)?;
    let u_poly = chebyshev_eval(ctx, bs, &chebyshev_u_coeffs(mu - 1), c);
    Some(normalize(ctx, ctx.mul(&[s, u_poly])))
}

/// The inverse/forward head pairs whose composites are banned (§2.4) —
/// the same eleven pairs [`contains_inv_composite`] enforces on the answer.
fn is_inv_pair(f: &str, g: &str) -> bool {
    matches!(
        (f, g),
        ("asin", "sin")
            | ("acos", "cos")
            | ("atan", "tan")
            | ("asinh", "sinh")
            | ("acosh", "cosh")
            | ("atanh", "tanh")
            | ("acsc", "csc")
            | ("asec", "sec")
            | ("acoth", "coth")
            | ("asech", "sech")
            | ("acsch", "csch")
    )
}

/// Back-substitution `u → F(w)` (spec §2.4): a whitelist rewrite of the
/// integrated t-form. `Var(t)` → `F(w)`; the eight trig/hyperbolic heads at
/// `μ·u` → their algebraic values; any other `t`-dependent node (inverse
/// heads, unrecognized arguments, leftover composites) declines the whole
/// substitution.
fn is_back_sub<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, bs: &IsBackSub<'a>) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        AtomNode::Var(v) => Some(if *v == bs.t_sym { bs.fw } else { expr }),
        AtomNode::Add(args) => {
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| is_back_sub(ctx, *a, bs))
                .collect::<Option<_>>()?;
            Some(ctx.add(&out))
        }
        AtomNode::Mul(args) => {
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| is_back_sub(ctx, *a, bs))
                .collect::<Option<_>>()?;
            Some(ctx.mul(&out))
        }
        AtomNode::Pow(b, e) => {
            let nb = is_back_sub(ctx, *b, bs)?;
            let ne = is_back_sub(ctx, *e, bs)?;
            // `(√E)^k` folds to `E^(k/2)` resp. `E^((k−1)/2)·√E` (W5b): the
            // answer never carries a powered `√E` literal.
            if let AtomNode::Num(k) = ne.node() {
                return Some(bs_power(ctx, bs, nb, *k));
            }
            Some(ctx.pow(nb, ne))
        }
        AtomNode::Fun(name, args) => {
            let n = name.as_str();
            if args.len() == 1 {
                if matches!(
                    n,
                    "sin"
                        | "cos"
                        | "tan"
                        | "sec"
                        | "cot"
                        | "csc"
                        | "sinh"
                        | "cosh"
                        | "tanh"
                        | "sech"
                        | "coth"
                        | "csch"
                ) {
                    return match u_multiple(args[0], bs.t_sym) {
                        Some(mu) => is_trig_back_map(ctx, bs, n, mu),
                        None if contains_symbol(args[0], bs.t_sym) => None,
                        None => {
                            let a0 = is_back_sub(ctx, args[0], bs)?;
                            Some(ctx.fun(n, &[a0]))
                        }
                    };
                }
                if inv_fun(n).is_some() && contains_symbol(args[0], bs.t_sym) {
                    // W5b cross-family algebraic value: the ratio tabular's
                    // `∫sech u du = asin(tanh u)` base maps to `asin(w)`.
                    // Only a *non-inverse-pair* trig/hyperbolic head under
                    // the inverse head maps this way — a same-pair
                    // `F(G(·))` composite stays a decline (§2.4 ban).
                    //
                    // Wave 5c: the asech family's `∫sech u du` base is the
                    // same `asin(tanh u)`, but there `w = sech u` is the
                    // identity head and `tanh u = √(1−w²)` is a radical, so
                    // the generic map would emit `asin(√(1−w²))` — a nested
                    // radical the certificate cannot relate back. The
                    // band identity `asin(tanh u) = acos(sech u) = acos(w)`
                    // states the same function without the nesting, so the
                    // asech family takes the identity-head form instead.
                    if bs.f == InvFun::Asech && n == "asin" {
                        if let AtomNode::Fun(inner, iargs) = args[0].node()
                            && iargs.len() == 1
                            && inner.as_str() == "tanh"
                            && let Some(mu) = u_multiple(iargs[0], bs.t_sym)
                            && mu == 1
                        {
                            return Some(ctx.fun("acos", &[bs.w]));
                        }
                    }
                    if let AtomNode::Fun(inner, iargs) = args[0].node()
                        && iargs.len() == 1
                        && !is_inv_pair(n, inner.as_str())
                        && matches!(
                            inner.as_str(),
                            "sin" | "cos" | "tan" | "sec" | "sinh" | "cosh" | "tanh" | "sech"
                        )
                        && let Some(mu) = u_multiple(iargs[0], bs.t_sym)
                    {
                        let v = is_trig_back_map(ctx, bs, inner.as_str(), mu)?;
                        return Some(ctx.fun(n, &[v]));
                    }
                    // An inverse head over any other `t`-dependent argument
                    // is a leftover composite — no whitelist value for it.
                    return None;
                }
            }
            let out: Vec<Atom<'a>> = args
                .iter()
                .map(|a| is_back_sub(ctx, *a, bs))
                .collect::<Option<_>>()?;
            Some(ctx.fun(n, &out))
        }
    }
}

/// Terminal structural check (spec §2.4/§8.2): the answer must not contain
/// an `F(G(·))` inverse-of-forward composite for any of the eleven pairs.
fn contains_inv_composite(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_inv_composite(*a))
        }
        AtomNode::Pow(b, e) => contains_inv_composite(*b) || contains_inv_composite(*e),
        AtomNode::Fun(name, args) => {
            if args.len() == 1
                && inv_fun(name.as_str()).is_some()
                && let AtomNode::Fun(inner, _) = args[0].node()
                && is_inv_pair(name.as_str(), inner.as_str())
            {
                return true;
            }
            args.iter().any(|a| contains_inv_composite(*a))
        }
    }
}

/// The four `(R, S)` families the ratio tabular integrates, parameterised by
/// the three signs that drive every reduction step:
///
/// | family     | `R`  | `S`  | `R′`     | `S′`      | `S²`     |
/// |------------|------|------|----------|-----------|----------|
/// | `TanSec`   | tan  | sec  | `+S²`    | `+R·S`    | `1 + R²` |
/// | `TanhSech` | tanh | sech | `+S²`    | `−R·S`    | `1 − R²` |
/// | `CotCsc`   | cot  | csc  | `−S²`    | `−R·S`    | `1 + R²` |
/// | `CothCsch` | coth | csch | `−S²`    | `−R·S`    | `R² − 1` |
///
/// `TanSec`/`TanhSech` are the pre-existing atan/atanh tabular; `CotCsc`/
/// `CothCsch` are the Wave 5c mirrors behind acsc/asech and acsch/acoth.
/// Everything below is expressed through
/// `R′ = ε·S²`, `S′ = δ·R·S`, `S² = η·(1 + γ·R²)` (with `β = η·γ`, so
/// `R² = β·(S² − η)`), which reproduces the two original families at
/// `(ε, δ, η, γ) = (1, 1, 1, 1)` and `(1, −1, 1, −1)` and covers the two new
/// ones at `(−1, −1, 1, 1)` and `(−1, −1, −1, −1)`. The sign algebra is what
/// makes a reciprocal head's answer come out with the right `±` throughout
/// (verified against SymPy's closed forms in the module tests).
#[derive(Clone, Copy, PartialEq, Eq)]
enum RatioFam {
    TanSec,
    TanhSech,
    CotCsc,
    CothCsch,
}

impl RatioFam {
    /// The `(R, S)` head names.
    fn heads(self) -> (&'static str, &'static str) {
        match self {
            RatioFam::TanSec => ("tan", "sec"),
            RatioFam::TanhSech => ("tanh", "sech"),
            RatioFam::CotCsc => ("cot", "csc"),
            RatioFam::CothCsch => ("coth", "csch"),
        }
    }

    /// `R′ = eps·S²`.
    fn eps(self) -> i64 {
        match self {
            RatioFam::TanSec | RatioFam::TanhSech => 1,
            RatioFam::CotCsc | RatioFam::CothCsch => -1,
        }
    }

    /// `S′ = delta·R·S`.
    fn delta(self) -> i64 {
        match self {
            RatioFam::TanSec => 1,
            RatioFam::TanhSech | RatioFam::CotCsc | RatioFam::CothCsch => -1,
        }
    }

    /// `S² = eta·(1 + gamma·R²)`.
    fn eta(self) -> i64 {
        match self {
            RatioFam::TanSec | RatioFam::TanhSech | RatioFam::CotCsc => 1,
            RatioFam::CothCsch => -1,
        }
    }

    fn gamma(self) -> i64 {
        match self {
            RatioFam::TanSec | RatioFam::CotCsc => 1,
            RatioFam::TanhSech | RatioFam::CothCsch => -1,
        }
    }

    /// `beta = eta·gamma`.
    fn beta(self) -> i64 {
        self.eta() * self.gamma()
    }

    /// Whether the class is the hyperbolic one (`sinh`/`cosh` rewrites).
    fn hyper(self) -> bool {
        matches!(self, RatioFam::TanhSech | RatioFam::CothCsch)
    }

    /// Coefficients (low-to-high in `R²`) of `S^(2k)`.
    fn s2_coeffs(self, k: i64) -> Vec<i64> {
        let (eta, beta) = (self.eta(), self.beta());
        (0..=k)
            .map(|i| binom(k, i) * eta.pow((k - i) as u32) * beta.pow(i as u32))
            .collect()
    }

    /// Coefficients (low-to-high in `S²`) of `R^(2k)`.
    fn r2_coeffs(self, k: i64) -> Vec<i64> {
        let (eta, beta) = (self.eta(), self.beta());
        (0..=k)
            .map(|i| binom(k, i) * beta.pow(k as u32) * (-eta).pow((k - i) as u32))
            .collect()
    }
}

/// `s·a` as an atom, folding `±1` into a bare sign (the tabular's ubiquitous
/// `±1` coefficients).
fn sc_mul<'a>(ctx: &'a AtomArena<'a>, s: i64, a: Atom<'a>) -> Atom<'a> {
    match s {
        1 => a,
        -1 => normalize(ctx, ctx.mul(&[ctx.num(-1), a])),
        _ => normalize(ctx, ctx.mul(&[ctx.num(s), a])),
    }
}

/// Decompose an expanded ratio-family t-term into `(constants, j, A, B)`
/// for `C·u^j·R(u)^A·S(u)^B`. A companion radical `(s·E)^(m/2)` contributes
/// `s^(m/2)·head^m·σ` to the same head power the `dx` factor's `G′` supplies,
/// so `B` is any integer here (negative `B` is an exact reciprocal power of
/// `cos`/`cosh` — see [`integrate_ratio_term`]).
fn scan_ratio_term<'a>(
    term: Atom<'a>,
    t_sym: Symbol,
    fam: RatioFam,
) -> Option<(Vec<Atom<'a>>, i64, i64, i64)> {
    let (tan, sec) = fam.heads();
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut j = 0i64;
    let mut tan_pow = 0i64;
    let mut sec_pow = 0i64;
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
            AtomNode::Fun(name, args)
                if args.len() == 1 && matches!(args[0].node(), AtomNode::Var(v) if *v == t_sym) =>
            {
                if name.as_str() == tan {
                    tan_pow += 1;
                } else if name.as_str() == sec {
                    sec_pow += 1;
                } else {
                    return None;
                }
            }
            AtomNode::Pow(b, e) => {
                let AtomNode::Fun(name, args) = b.node() else {
                    return None;
                };
                if args.len() != 1 || !matches!(args[0].node(), AtomNode::Var(v) if *v == t_sym) {
                    return None;
                }
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                if name.as_str() == tan && *n > 0 {
                    tan_pow += n;
                } else if name.as_str() == sec && *n != 0 {
                    sec_pow += n;
                } else {
                    return None;
                }
            }
            _ => return None,
        }
    }
    if j > MAX_IS_K || tan_pow > MAX_IS_POLY_DEG as i64 || sec_pow.abs() > 2 + MAX_IS_RAD_M {
        return None;
    }
    Some((coeff, j, tan_pow, sec_pow))
}

/// `I(j, a) = ∫ u^j·R(u)^a·S(u)² du`: parts on the exact differential
/// `dv = R^a·S² du` (`v = ε·R(u)^(a+1)/(a+1)`, since `R′ = ε·S²`), reducing
/// to [`ratio_j`]. `None` when the recursion reaches a polylog term (§2.5).
fn ratio_i<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    a: i64,
    t: Atom<'a>,
    fam: RatioFam,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(0..=2 * MAX_IS_MU).contains(&a) {
        return None;
    }
    let (tan, _) = fam.heads();
    let eps = fam.eps();
    let r_t = ctx.fun(tan, &[t]);
    let boundary = normalize(
        ctx,
        ctx.mul(&[
            ctx.num(eps),
            int_pow(ctx, t, j),
            int_pow(ctx, r_t, a + 1),
            inv(ctx, ctx.num(a + 1)),
        ]),
    );
    if j == 0 {
        return Some(boundary);
    }
    let inner = ratio_j(ctx, j - 1, a + 1, t, fam)?;
    let coeff = rat_atom(ctx, -j * eps, a + 1);
    Some(normalize(
        ctx,
        ctx.add(&[boundary, ctx.mul(&[coeff, inner])]),
    ))
}

/// `J(j, c) = ∫ u^j·R(u)^c du`. With `R² = β·S² − β·η`:
/// `J(c) = β·I(c−2) − β·η·J(c−2)` — `I − J` for `(tan, sec)` and
/// `(cot, csc)`, `J − I` for `(tanh, sech)`, `I + J` for `(coth, csch)`.
/// Bases: `J(j, 0) = u^{j+1}/(j+1)` and `J(0, 1) = −log(cos u)` (tan),
/// `log(cosh u)` (tanh), `log(sin u)` (cot), `log(sinh u)` (coth);
/// `J(j ≥ 1, 1)` is polylog territory and declines (§2.5 — elementarity
/// emerges from the reduction base case, the residue-check analogue).
fn ratio_j<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    c: i64,
    t: Atom<'a>,
    fam: RatioFam,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(0..=2 * MAX_IS_MU).contains(&c) {
        return None;
    }
    match c {
        0 => Some(normalize(
            ctx,
            ctx.mul(&[int_pow(ctx, t, j + 1), inv(ctx, ctx.num(j + 1))]),
        )),
        1 => {
            if j >= 1 {
                return None;
            }
            let (head, sign) = match fam {
                RatioFam::TanSec => ("cos", -1),
                RatioFam::TanhSech => ("cosh", 1),
                RatioFam::CotCsc => ("sin", 1),
                RatioFam::CothCsch => ("sinh", 1),
            };
            let log = ctx.fun("log", &[ctx.fun(head, &[t])]);
            Some(if sign < 0 {
                normalize(ctx, ctx.mul(&[ctx.num(-1), log]))
            } else {
                log
            })
        }
        _ => {
            let jm = ratio_j(ctx, j, c - 2, t, fam)?;
            let im = ratio_i(ctx, j, c - 2, t, fam)?;
            let (beta, be) = (fam.beta(), fam.beta() * fam.eta());
            Some(normalize(
                ctx,
                ctx.add(&[sc_mul(ctx, beta, im), sc_mul(ctx, -be, jm)]),
            ))
        }
    }
}

/// `∫ C·u^j·R(u)^A·S(u)^B du` (W5b: any integer `B`). `B < 0` converts to
/// the exact `sin/cosh`-family form and re-uses the trig/hyperbolic
/// machinery; `B ≥ 0` goes through [`ratio_rs`].
fn integrate_ratio_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    t_sym: Symbol,
    fam: RatioFam,
) -> Option<Atom<'a>> {
    let (coeff, j, a, b) = scan_ratio_term(term, t_sym, fam)?;
    let t = ctx.var(t_sym.as_str());
    if b < 0 {
        // `sec = 1/cos` / `sech = 1/cosh` / `csc = 1/sin` / `csch = 1/sinh`
        // hold exactly (no branch), so `R^A·S^B = p(u)^A·q(u)^(−A−B)` with
        // `p` the numerator head of the family: a pure `sin/cos`
        // (`sinh/cosh`) term, integrable by the trig/hyperbolic tabulars when
        // `−A−B ≥ 0`. Otherwise the reciprocal power survives — decline
        // (§2.5). The reciprocal families never take this branch (their `S`
        // power is at least 1 by construction), but the conversion is exact
        // for all four, so it is kept total.
        let c_pow = -a - b;
        if c_pow < 0 {
            return None;
        }
        let (s_name, c_name) = match fam {
            RatioFam::TanSec => ("sin", "cos"),
            RatioFam::TanhSech => ("sinh", "cosh"),
            RatioFam::CotCsc => ("cos", "sin"),
            RatioFam::CothCsch => ("cosh", "sinh"),
        };
        let mut fs = coeff;
        if j > 0 {
            fs.push(int_pow(ctx, t, j));
        }
        if a > 0 {
            fs.push(int_pow(ctx, ctx.fun(s_name, &[t]), a));
        }
        if c_pow > 0 {
            fs.push(int_pow(ctx, ctx.fun(c_name, &[t]), c_pow));
        }
        let conv = normalize(ctx, ctx.mul(&fs));
        return integrate_sincos_term(ctx, conv, t_sym, fam.hyper());
    }
    let core = ratio_rs(ctx, j, a, b, t, fam)?;
    let mut fs = coeff;
    fs.push(core);
    Some(normalize(ctx, ctx.mul(&fs)))
}

/// Integrate a `C·u^j·s(u)^A·c(u)^B` term (`sin/cos` or `sinh/cosh` at the
/// bare angle `u`): the same multiple-angle reduction + tabular the
/// Asin/Acos and Asinh/Acosh dispatch branches run, factored out for the
/// ratio family's negative-`S`-power conversions (W5b).
fn integrate_sincos_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    t_sym: Symbol,
    hyper: bool,
) -> Option<Atom<'a>> {
    if hyper {
        let t = ctx.var(t_sym.as_str());
        let (coeff, j, na, nb) = scan_hyp_term(ctx, term, t_sym, MAX_IS_K)?;
        let mut reduced: Vec<Atom<'a>> = Vec::new();
        reduce_hyp_term(ctx, &coeff, j, na, nb, t, &mut reduced)?;
        let t_form = normalize(ctx, ctx.add(&reduced));
        if node_count(t_form) > MAX_REDUCED_NODES {
            return None;
        }
        let rterms: Vec<Atom<'a>> = match t_form.node() {
            AtomNode::Add(args) => args.to_vec(),
            _ => vec![t_form],
        };
        let mut out = Vec::with_capacity(rterms.len());
        for rt in rterms {
            out.push(integrate_hyp_term(ctx, rt, t_sym, MAX_IS_K)?);
        }
        return Some(normalize(ctx, ctx.add(&out)));
    }
    let reduced = match super::trig_reduce::trig_reduce_products(ctx, term, t_sym) {
        Some(r) => crate::ode::util::collect_terms(ctx, normalize(ctx, r)),
        None => term,
    };
    if node_count(reduced) > MAX_REDUCED_NODES {
        return None;
    }
    let rterms: Vec<Atom<'a>> = match reduced.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![reduced],
    };
    let mut out = Vec::with_capacity(rterms.len());
    for rt in rterms {
        out.push(integrate_trig_term(ctx, rt, t_sym)?);
    }
    Some(normalize(ctx, ctx.add(&out)))
}

/// `∫ u^j·R(u)^A·S(u)^B du` with `B ≥ 0` (any [`RatioFam`]; the two
/// pre-existing families have `R′ = S²`, `S′ = ±R·S` with `−` hyperbolic, the
/// two Wave 5c mirrors `R′ = −S²`, `S′ = −R·S`).
///
/// - even `B ≥ 2`: `S^B = S²·(S²)^((B−2)/2)` expands to `R^a·S²` terms for
///   the [`ratio_i`] tabular (`fam.s2_coeffs`);
/// - odd `B`, odd `A ≥ 1`: parts on `d(S^B) = δ·B·R·S^B du` with
///   `R^(A−1) = (R²)^((A−1)/2)` a polynomial in `S²`, so the primitive of
///   `dv = R^A·S^B du` is a polynomial in `S` and the reduced terms are
///   pure `S` powers ([`ratio_s_pow`]);
/// - odd `B`, even `A`: `R^A = (R²)^(A/2)` expands likewise;
/// - `B = 0`: [`ratio_j`].
fn ratio_rs<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    a: i64,
    b: i64,
    t: Atom<'a>,
    fam: RatioFam,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(0..=2 * MAX_IS_MU).contains(&a) {
        return None;
    }
    if b == 0 {
        return ratio_j(ctx, j, a, t, fam);
    }
    if !(1..=2 * MAX_IS_MU).contains(&b) {
        return None;
    }
    let (_, sec) = fam.heads();
    let s_t = ctx.fun(sec, &[t]);
    if b % 2 == 0 {
        // S^(2k) = Σ_i c_i·R^(2i) with `c = fam.s2_coeffs(k)`, k = (B−2)/2.
        let k = (b - 2) / 2;
        let mut terms = Vec::new();
        for (i, c) in fam.s2_coeffs(k).iter().enumerate() {
            let core = ratio_i(ctx, j, a + 2 * i as i64, t, fam)?;
            terms.push(sc_mul(ctx, *c, core));
        }
        return Some(normalize(ctx, ctx.add(&terms)));
    }
    if a == 0 {
        return ratio_s_pow(ctx, j, b, t, fam);
    }
    let sgn = fam.delta();
    if a % 2 == 1 {
        // dv = R^A·S^B du: with `R^(A−1) = Σ_i c_i·S^(2i)` the primitive is
        // `v = Σ_i c_i·S^(B+2i)/(δ·(B+2i))` (since `d(S^p) = δ·p·R·S^p du`).
        let k = (a - 1) / 2;
        let v_terms: Vec<(i64, i64)> = fam
            .r2_coeffs(k)
            .iter()
            .enumerate()
            .map(|(i, c)| (*c, b + 2 * i as i64))
            .collect();
        let mut v_atoms = Vec::with_capacity(v_terms.len());
        let mut inner = Vec::with_capacity(v_terms.len());
        for (c, p) in v_terms {
            let coeff = rat_atom(ctx, c, sgn * p);
            v_atoms.push(normalize(ctx, ctx.mul(&[coeff, int_pow(ctx, s_t, p)])));
            if j > 0 {
                let red = ratio_s_pow(ctx, j - 1, p, t, fam)?;
                inner.push(normalize(ctx, ctx.mul(&[coeff, red])));
            }
        }
        let v = normalize(ctx, ctx.add(&v_atoms));
        let boundary = normalize(ctx, ctx.mul(&[int_pow(ctx, t, j), v]));
        if j == 0 {
            return Some(boundary);
        }
        let inner = normalize(ctx, ctx.add(&inner));
        return Some(normalize(
            ctx,
            ctx.add(&[boundary, ctx.mul(&[ctx.num(-j), inner])]),
        ));
    }
    // Even `A ≥ 2`, odd `B`: `R^A = (R²)^(A/2)` → `Σ d_i·S^(B+2i)`.
    let k = a / 2;
    let mut terms = Vec::new();
    for (i, c) in fam.r2_coeffs(k).iter().enumerate() {
        let red = ratio_s_pow(ctx, j, b + 2 * i as i64, t, fam)?;
        terms.push(sc_mul(ctx, *c, red));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// `∫ u^j·S(u)^b du`, odd `b ≥ 1`. With `d(R·S^(b−2)) = ε·S^b + (b−2)·δ·β·S^b
/// − (b−2)·δ·β·η·S^(b−2)` the `j = 0` reduction is
/// `∫S^b = (R·S^(b−2) + (b−2)·δ·β·η·∫S^(b−2))/(ε + (b−2)·δ·β)`, which
/// specialises to the pre-existing `(R·S^(b−2) + (b−2)∫S^(b−2))/(b−1)` for
/// both `(tan, sec)` and `(tanh, sech)`. The bases are
/// `∫sec u du = log(sec u + tan u)`, `∫sech u du = asin(tanh u)`,
/// `∫csc u du = −log(csc u + cot u)` and
/// `∫csch u du = −log(csch u + coth u)` (all exact on the principal band; the
/// sech base back-substitutes through the whitelist's cross-family head map).
///
/// `j ≥ 1` declines: since `ratio_rs` only routes *odd* `S` powers here, the
/// `b`-decreasing reduction bottoms out at `∫u^j·S(u) du`, which is polylog
/// territory in all four families (§2.5) — the same base case that
/// [`ratio_j`] refuses at `c = 1`. (An earlier revision of this helper
/// instead substituted `∫u^(j−1)·S^(b−2)` for `∫u^(j−1)·H_(b−2) du` and
/// emitted a *wrong* closed form there; the module's numeric tabular
/// self-check in the tests now covers the whole `(j, a, b)` grid of all four
/// families.)
fn ratio_s_pow<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    b: i64,
    t: Atom<'a>,
    fam: RatioFam,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(1..=2 * MAX_IS_MU).contains(&b) {
        return None;
    }
    let (tan, sec) = fam.heads();
    let r_t = ctx.fun(tan, &[t]);
    let s_t = ctx.fun(sec, &[t]);
    if b == 1 {
        if j >= 1 {
            return None;
        }
        return Some(match fam {
            // d/du asin(tanh u) = sech²u/√(1 − tanh²u) = sech u (sech > 0).
            RatioFam::TanhSech => ctx.fun("asin", &[r_t]),
            RatioFam::TanSec => ctx.fun("log", &[normalize(ctx, ctx.add(&[s_t, r_t]))]),
            // d/du(−log(csc u + cot u)) = csc u; d/du(−log(csch u + coth u))
            // = csch u (both arguments positive on the `u > 0` band).
            RatioFam::CotCsc | RatioFam::CothCsch => normalize(
                ctx,
                ctx.mul(&[
                    ctx.num(-1),
                    ctx.fun("log", &[normalize(ctx, ctx.add(&[s_t, r_t]))]),
                ]),
            ),
        });
    }
    if j >= 1 {
        return None;
    }
    // `H_b = (R·S^(b−2) + (b−2)·δ·β·η·H_{b−2})/(ε + (b−2)·δ·β)`.
    let db = fam.delta() * fam.beta();
    let den = fam.eps() + (b - 2) * db;
    if den == 0 {
        return None;
    }
    let num2 = (b - 2) * db * fam.eta();
    let rec = ratio_s_pow(ctx, 0, b - 2, t, fam)?;
    Some(normalize(
        ctx,
        ctx.mul(&[
            inv(ctx, ctx.num(den)),
            ctx.add(&[
                ctx.mul(&[r_t, int_pow(ctx, s_t, b - 2)]),
                sc_mul(ctx, num2, rec),
            ]),
        ]),
    ))
}

/// `Some(μ)` when `arg = μ·t` with `μ` an integer constant of either sign
/// (`t` itself → 1); `None` for any other shape.
fn signed_u_multiple(arg: Atom<'_>, t_sym: Symbol) -> Option<i64> {
    match arg.node() {
        AtomNode::Num(n) => Some(*n),
        AtomNode::Var(v) if *v == t_sym => Some(1),
        AtomNode::Mul(args) => {
            let mut mu: Option<i64> = None;
            let mut seen_t = false;
            for a in args.iter() {
                match a.node() {
                    AtomNode::Num(n) if mu.is_none() => mu = Some(*n),
                    AtomNode::Var(v) if *v == t_sym && !seen_t => seen_t = true,
                    _ => return None,
                }
            }
            if seen_t { Some(mu.unwrap_or(1)) } else { None }
        }
        _ => None,
    }
}

/// Closed form of one (reduced) trig t-term `C·u^j·sin/cos(μ·u)`
/// (`0 ≤ j ≤ MAX_IS_K`, `1 ≤ |μ| ≤ MAX_IS_MU`) or of the pure polynomial
/// `C·u^j`: repeated parts, the trig analogue of [`integrate_hyp_term`]
/// (spec §2.3's T3 analogue):
/// `∫u^j·H(μu) du = Σ_{i=0..j} (−1)^i·j^underline_i·u^{j−i}·G_i(μu)/μ^{i+1}`
/// with `G` cycling `[sin, −cos, −sin, cos]` for `H = cos` and
/// `[−cos, −sin, cos, sin]` for `H = sin`. Any other shape declines.
fn integrate_trig_term<'a>(
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
    let mut trig: Option<(bool, i64)> = None;
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
                let is_sin = match name.as_str() {
                    "sin" => true,
                    "cos" => false,
                    _ => return None,
                };
                if trig.is_some() {
                    return None;
                }
                trig = Some((is_sin, signed_u_multiple(args[0], t_sym)?));
            }
            _ => return None,
        }
    }
    if j > MAX_IS_K {
        return None;
    }
    let t = ctx.var(t_sym.as_str());
    let (is_sin, mu) = match trig {
        // C·u^j → C·u^(j+1)/(j+1)
        None => {
            let mut fs = coeff;
            fs.push(int_pow(ctx, t, j + 1));
            fs.push(inv(ctx, ctx.num(j + 1)));
            return Some(normalize(ctx, ctx.mul(&fs)));
        }
        Some(h) => h,
    };
    // sin(−μu) = −sin(μu), cos(−μu) = cos(μu): fold a negative multiple
    // into the coefficient.
    let (neg, mu) = if mu < 0 { (is_sin, -mu) } else { (false, mu) };
    if mu > MAX_IS_MU {
        return None;
    }
    if mu == 0 {
        // sin(0) = 0 kills the term; cos(0) = 1 leaves a polynomial.
        if is_sin {
            return Some(ctx.num(0));
        }
        let mut fs = coeff;
        if neg {
            fs.push(ctx.num(-1));
        }
        fs.push(int_pow(ctx, t, j + 1));
        fs.push(inv(ctx, ctx.num(j + 1)));
        return Some(normalize(ctx, ctx.mul(&fs)));
    }
    let mut terms = Vec::new();
    let mut fall: i64 = 1; // falling factorial j!/(j−i)!
    for i in 0..=j {
        if i > 0 {
            fall = fall.checked_mul(j - i + 1)?;
        }
        let mpow = mu.checked_pow((i + 1) as u32)?;
        let (h_name, h_sign) = match (is_sin, i % 4) {
            (true, 0) => ("cos", -1),
            (true, 1) => ("sin", -1),
            (true, 2) => ("cos", 1),
            (true, _) => ("sin", 1),
            (false, 0) => ("sin", 1),
            (false, 1) => ("cos", -1),
            (false, 2) => ("sin", -1),
            (false, _) => ("cos", 1),
        };
        let marg = if mu == 1 {
            t
        } else {
            ctx.mul(&[ctx.num(mu), t])
        };
        let mut fs = coeff.clone();
        if ((i % 2 == 1) != (h_sign < 0)) != neg {
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
        fs.push(ctx.fun(h_name, &[marg]));
        terms.push(normalize(ctx, ctx.mul(&fs)));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// Match a radical factor `R^(m/2)` with integer `0 < |m| ≤ MAX_IS_RAD_M`
/// (W5b): the `sqrt(R)` / `sqrt(R)^n` spellings (`m = 1` / `m = n`),
/// half-power exponents `R^(p/2)` (`m = p`), and negative integer powers
/// `R^p` (`m = 2p`, the squared-kernel reciprocals — `|2p| ≤ 5`). Positive
/// integer powers are polynomial weights and never reach here (the caller
/// tries [`poly_deg`] first). Companionship of `R` is *not* checked here;
/// the caller runs [`radicand_scale`] against the family kernel.
fn match_radical_pow(factor: Atom<'_>) -> Option<(Atom<'_>, i64)> {
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
                if *n != 0 && n.unsigned_abs() <= MAX_IS_RAD_M as u64 {
                    return Some((args[0], *n));
                }
                return None;
            }
            let (p, q) = rat_of(*e)?;
            let m = match (p, q) {
                (p, 2) if p != 0 && p.unsigned_abs() <= MAX_IS_RAD_M as u64 => p,
                (p, 1) if (-2..=-1).contains(&p) => 2 * p,
                _ => return None,
            };
            Some((*b, m))
        }
        _ => None,
    }
}

/// True when the imaginary unit (the corpus spelling `i` or the internal
/// `I`) occurs anywhere in `expr` — M4 emits no complex answers, so an
/// `i`-carrying radicand or normalization declines (rubi-00408).
fn mentions_imaginary(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => matches!(v.as_str(), "i" | "I"),
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|a| mentions_imaginary(*a))
        }
        AtomNode::Pow(b, e) => mentions_imaginary(*b) || mentions_imaginary(*e),
    }
}

/// M4 entry (spec §3.1): `∫ C·(a + b·F(w))^k · P(x) dx` by the
/// inverse-function substitution `u = F(w)`. See the section docs above
/// for the scope, the branch semantics, and the deviations from the spec.
pub(crate) fn integrate_inv_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // Cheap structural gates first (spec §3.1 and §8 risk 7): the decline
    // paths stay O(nodes); no substitution work runs before the shape is
    // confirmed. Inside a parts scope the substitution is suppressed
    // outright (spec §3.1 — all or nothing there).
    if super::chain::parts_active() {
        return None;
    }
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    if node_count(expr) > MAX_NODES {
        return None;
    }
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut inv_f: Option<InvPow<'a>> = None;
    let mut polys: Vec<Atom<'a>> = Vec::new();
    let mut rad: Option<(Atom<'a>, i64)> = None;
    let mut total_deg = 0usize;
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if inv_f.is_none()
            && let Some(m) = match_inv_power(ctx, f, var, MAX_IS_K)
        {
            inv_f = Some(m);
            continue;
        }
        if let Some(d) = poly_deg(f, var)
            && d >= 1
            && total_deg + d <= MAX_IS_POLY_DEG
        {
            total_deg += d;
            polys.push(f);
            continue;
        }
        // W5b: at most one radical factor; companionship (the family's own
        // kernel at the same argument) is verified after the affine gate.
        if rad.is_none()
            && let Some(rm) = match_radical_pow(f)
        {
            rad = Some(rm);
            continue;
        }
        // F-in-denom, extra radical bases, negative or symbolic powers of
        // anything else → decline.
        return None;
    }
    let m = inv_f?;
    let (sigma, rho) = linear_form(ctx, m.arg, var)?;
    if is_zero_atom(ctx, sigma) {
        return None;
    }
    let mut work: u64 = node_count(expr) as u64;
    // W5b companion-radical absorption (§2.2): `(s·E(w))^(m/2)` maps to
    // `s^(m/2)·G′(u)^m·σ`. Companionship is P3's `radicand_scale` (same
    // argument `w`, same family kernel, constant normalization quotient).
    let mut rad_data: Option<(i64, Atom<'a>)> = None;
    if let Some((r, rm)) = rad {
        let e_u = family_kernel_of_u(ctx, m.f, m.arg);
        // `x = −ρ/σ`, the abscissa where the argument vanishes and
        // `E(u) = e0 ∈ {±1}` (the P3 candidate point).
        let x_root = normalize(ctx, ctx.mul(&[ctx.num(-1), rho, inv(ctx, sigma)]));
        let s_atom = radicand_scale(ctx, r, e_u, x_root, var)?;
        // Odd `m` takes the half power of `s`: a numeric normalization must
        // be positive (real square root), a symbolic one is assumed
        // positive (the module's standing sign convention).
        if rm % 2 != 0 && !sqrt_norm_positive(s_atom) {
            return None;
        }
        // No `I` (spec §1.3-7): an imaginary radicand normalization would
        // put `i^(m/2)` into the answer — decline (rubi-00408 family).
        if mentions_imaginary(r) || mentions_imaginary(s_atom) {
            return None;
        }
        // Charge the companion check against the M4 work budget.
        if is_work_exceeded(&mut work, 2 * node_count(r)) {
            return None;
        }
        rad_data = Some((rm, s_atom));
    }
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let g_t = ctx.fun(family_g_head(m.f), &[t]);
    // x = (G(u) − ρ)/σ
    let x_of_t = {
        let num_t = if is_zero_atom(ctx, rho) {
            g_t
        } else {
            normalize(ctx, ctx.add(&[g_t, ctx.mul(&[ctx.num(-1), rho])]))
        };
        normalize(ctx, ctx.mul(&[num_t, inv(ctx, sigma)]))
    };
    let mut t_factors = rest;
    let base_t = normalize(ctx, ctx.add(&[m.a, ctx.mul(&[m.b, t])]));
    t_factors.push(int_pow(ctx, base_t, m.k));
    for p in &polys {
        t_factors.push(normalize(ctx, replace_symbol(ctx, *p, var, x_of_t)));
    }
    // dx = gp_sign·G′(u)^…/σ·du, with the Wave 5c heads contributing *two*
    // head factors (`acsc`: `−csc·cot`, `asec`: `sec·tan`, `asech`:
    // `−sech·tanh`, `acsch`: `−csch·coth`). A companion radical contributes
    // `s^(rm/2)·head^rm·σ` where `head` is the family's kernel head
    // (`E(G(u)) = head²`), which merges additively into that head's power
    // (§2.2). The sheet factor σ (odd `rm`) is `+1` on the principal band of
    // every family — see the section docs for the §4.2 emission record.
    let (gp_sign, gp_terms) = family_gp_terms(m.f);
    if gp_sign < 0 {
        t_factors.push(ctx.num(-1));
    }
    let mut head_powers: Vec<(&'static str, i64)> = gp_terms.to_vec();
    if let Some((rm, s_atom)) = &rad_data {
        let kernel_head = family_kernel_head(m.f);
        match head_powers.iter_mut().find(|(h, _)| *h == kernel_head) {
            Some(slot) => slot.1 += rm,
            None => head_powers.push((kernel_head, *rm)),
        }
        t_factors.push(ctx.pow(*s_atom, rat_atom(ctx, *rm, 2)));
        if std::env::var_os("OCAS_M4_DEBUG").is_some() {
            eprintln!(
                "[m4] companion radical (s·E)^({rm}/2) absorbed; sheet factor σ folded to +1 (principal band)"
            );
        }
    }
    {
        let mut gp_factors: Vec<Atom<'a>> = Vec::with_capacity(head_powers.len() + 1);
        for (head, pow) in &head_powers {
            if *pow != 0 {
                gp_factors.push(int_pow(ctx, ctx.fun(head, &[t]), *pow));
            }
        }
        gp_factors.push(inv(ctx, sigma));
        t_factors.push(normalize(ctx, ctx.mul(&gp_factors)));
    }
    let product = normalize(ctx, ctx.mul(&t_factors));
    if is_work_exceeded(&mut work, node_count(product)) {
        return None;
    }
    let expanded = crate::expand::expand_bounded(ctx, product).unwrap_or(product);
    // `x^m` maps to `(G(u)/σ)^m`, a power of a *product* — distribute those
    // integer powers (`flatten_powers`, the P3 helper) so every term is a
    // flat `C·u^j·H(u)^…` product for the per-term scanners.
    let mut flat_budget = MAX_CANCEL_BUDGET;
    let flattened = flatten_powers(ctx, expanded, &mut flat_budget);
    let folded = crate::ode::util::collect_terms(ctx, normalize(ctx, flattened));
    if node_count(folded) > MAX_IS_SUBST_NODES || is_work_exceeded(&mut work, node_count(folded)) {
        return None;
    }
    let terms: Vec<Atom<'a>> = match folded.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![folded],
    };
    let t_answer = match m.f {
        // Hyperbolic families: M2's multiple-angle reduction plus the
        // tabular closed form directly (no chain-first attempt — see the
        // section docs).
        InvFun::Asinh | InvFun::Acosh => {
            let mut reduced: Vec<Atom<'a>> = Vec::new();
            for term in terms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                let (coeff, j, na, nb) = scan_hyp_term(ctx, term, t_sym, MAX_IS_K)?;
                reduce_hyp_term(ctx, &coeff, j, na, nb, t, &mut reduced)?;
            }
            let t_form = normalize(ctx, ctx.add(&reduced));
            if node_count(t_form) > MAX_REDUCED_NODES
                || is_work_exceeded(&mut work, node_count(t_form))
            {
                return None;
            }
            let rterms: Vec<Atom<'a>> = match t_form.node() {
                AtomNode::Add(args) => args.to_vec(),
                _ => vec![t_form],
            };
            let mut out = Vec::with_capacity(rterms.len());
            for term in rterms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                out.push(integrate_hyp_term(ctx, term, t_sym, MAX_IS_K)?);
            }
            normalize(ctx, ctx.add(&out))
        }
        // Ratio families (atan/atanh and the Wave 5c reciprocals): the
        // t-terms are `C·u^j·R(u)^A·S(u)^B` in the family's own `(R, S)`
        // pair (`(tan, sec)`, `(tanh, sech)`, `(cot, csc)`, `(coth, csch)`),
        // integrated by the parts/R-reduction tabular (a polylog base case
        // declines). The reciprocal rows are the sign-flipped mirrors of the
        // two originals — see [`RatioFam`].
        InvFun::Atan
        | InvFun::Atanh
        | InvFun::Acsc
        | InvFun::Asec
        | InvFun::Acoth
        | InvFun::Asech
        | InvFun::Acsch => {
            let fam = match m.f {
                InvFun::Atan | InvFun::Asec => RatioFam::TanSec,
                InvFun::Atanh | InvFun::Asech => RatioFam::TanhSech,
                InvFun::Acsc => RatioFam::CotCsc,
                InvFun::Acoth | InvFun::Acsch => RatioFam::CothCsch,
                _ => return None,
            };
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                out.push(integrate_ratio_term(ctx, term, t_sym, fam)?);
            }
            normalize(ctx, ctx.add(&out))
        }
        // Trig families (asin/acos): sin/cos products reduce to multiple
        // angles in-module (trig_reduce_products), then the T3-analogue
        // tabular (`integrate_trig_term`) finishes each
        // `C·u^j·sin/cos(μu)` term without touching the chain (§2.3).
        InvFun::Asin | InvFun::Acos => {
            let mut out: Vec<Atom<'a>> = Vec::new();
            for term in terms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                let reduced = match super::trig_reduce::trig_reduce_products(ctx, term, t_sym) {
                    Some(r) => crate::ode::util::collect_terms(ctx, normalize(ctx, r)),
                    None => term,
                };
                if node_count(reduced) > MAX_REDUCED_NODES
                    || is_work_exceeded(&mut work, node_count(reduced))
                {
                    return None;
                }
                let rterms: Vec<Atom<'a>> = match reduced.node() {
                    AtomNode::Add(args) => args.to_vec(),
                    _ => vec![reduced],
                };
                for rt in rterms {
                    match integrate_trig_term(ctx, rt, t_sym) {
                        Some(v) => out.push(v),
                        None => {
                            if std::env::var_os("OCAS_M4_DEBUG").is_some() {
                                eprintln!("[m4] trig term declined: {rt} (from {term})");
                            }
                            return None;
                        }
                    }
                }
            }
            normalize(ctx, ctx.add(&out))
        }
    };
    // Back-substitution (§2.4) plus the terminal checks.
    let w = m.arg;
    let e_w = family_kernel_of_u(ctx, m.f, w);
    let sqrt_e = ctx.fun("sqrt", &[e_w]);
    let fw = ctx.fun(m.f.name(), &[w]);
    let bs = IsBackSub {
        t_sym,
        f: m.f,
        fw,
        w,
        e_w,
        sqrt_e,
    };
    let ans = is_back_sub(ctx, t_answer, &bs)?;
    if is_work_exceeded(&mut work, node_count(ans)) {
        return None;
    }
    let ans = normalize(ctx, ans);
    if contains_symbol(ans, t_sym) || contains_integral(ans) || contains_inv_composite(ans) {
        return None;
    }
    Some(ans)
}

// =========================================================================
// Combined entry
// =========================================================================

/// Try the inverse-trig/hyperbolic mechanisms in cost order:
///
/// 1. **P1** — the inverse-composition cancellation pre-pass. When it fires,
///    the rewritten integrand (expanded and folded) is re-integrated through
///    `integrate_raw`; a residue there just falls through to the stages below
///    (which see the *original* integrand, since only a correct answer may be
///    returned). A second application is a no-op, so the chain re-entry
///    terminates.
/// 2. bare linear forms (M3), kernel-derivative powers (M1), companion
///    radicands (P3), branch-shifted composition kernels (P2), the inv-hyp
///    substitution (M2), and finally the inverse-function substitution
///    (M4 — the most general and most expensive, so it runs last).
pub(crate) fn integrate_inverse_trig<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    let mut budget = MAX_CANCEL_BUDGET;
    if let Some(rewritten) = cancel_compositions(ctx, expr, var, &mut budget) {
        let rewritten = normalize(ctx, rewritten);
        if node_count(rewritten) <= MAX_NODES {
            let expanded = crate::expand::expand_bounded(ctx, rewritten).unwrap_or(rewritten);
            let folded = crate::ode::util::collect_terms(ctx, normalize(ctx, expanded));
            if node_count(folded) <= MAX_KERNEL_NODES
                && let Some(r) = integrate_expanded_terms(ctx, folded, var)
            {
                return Some(normalize(ctx, r));
            }
        }
    }
    integrate_bare_linear(ctx, expr, var)
        .or_else(|| integrate_kernel_power(ctx, expr, var))
        .or_else(|| integrate_companion_form(ctx, expr, var))
        .or_else(|| integrate_affine_kernel(ctx, expr, var))
        .or_else(|| integrate_invhyp_subst(ctx, expr, var))
        .or_else(|| integrate_inv_subst(ctx, expr, var))
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
                    "cot" => 1.0 / v.tan(),
                    "sec" => 1.0 / v.cos(),
                    "csc" => 1.0 / v.sin(),
                    "exp" => v.exp(),
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "tanh" => v.tanh(),
                    "coth" => 1.0 / v.tanh(),
                    "sech" => 1.0 / v.cosh(),
                    "csch" => 1.0 / v.sinh(),
                    "asin" => v.asin(),
                    "acos" => v.acos(),
                    "atan" => v.atan(),
                    "acot" => std::f64::consts::FRAC_PI_2 - v.atan(),
                    "asinh" => v.asinh(),
                    "acosh" => v.acosh(),
                    "atanh" => v.atanh(),
                    "acoth" => (1.0 / v).atanh(),
                    // Wave 5c reciprocal heads: `acsc(w) = asin(1/w)` etc.
                    // (the samples always sit on the positive band).
                    "acsc" => (1.0 / v).asin(),
                    "asec" => (1.0 / v).acos(),
                    "asech" => (1.0 / v).acosh(),
                    "acsch" => (1.0 / v).asinh(),
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

    /// Structural replacement of `target` by `repl` anywhere in `expr`.
    fn replace_atom<'a>(
        ctx: &'a AtomArena<'a>,
        expr: Atom<'a>,
        target: Atom<'a>,
        repl: Atom<'a>,
    ) -> Atom<'a> {
        if expr == target {
            return repl;
        }
        match expr.node() {
            AtomNode::Num(_) | AtomNode::Var(_) => expr,
            AtomNode::Add(args) => {
                let out: Vec<Atom<'a>> = args
                    .iter()
                    .map(|a| replace_atom(ctx, *a, target, repl))
                    .collect();
                ctx.add(&out)
            }
            AtomNode::Mul(args) => {
                let out: Vec<Atom<'a>> = args
                    .iter()
                    .map(|a| replace_atom(ctx, *a, target, repl))
                    .collect();
                ctx.mul(&out)
            }
            AtomNode::Pow(b, e) => {
                let nb = replace_atom(ctx, *b, target, repl);
                let ne = replace_atom(ctx, *e, target, repl);
                ctx.pow(nb, ne)
            }
            AtomNode::Fun(name, args) => {
                let out: Vec<Atom<'a>> = args
                    .iter()
                    .map(|a| replace_atom(ctx, *a, target, repl))
                    .collect();
                ctx.fun(name.as_str(), &out)
            }
        }
    }

    /// Verify an antiderivative of an integrand carrying a branch-shifted
    /// composition kernel `K = F(G(u))` (whose value is `u` plus a branch
    /// constant). Replacing `K` by the affine proxy `u` in *both* the
    /// integrand and the answer leaves a consistent pair — the branch constant
    /// cancels out of every derivative and `K' = u' = σ` still holds — so the
    /// numeric check exercises exactly the kernel algebra.
    fn assert_antiderivative_branch_shifted<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        kernel: Atom<'a>,
        proxy: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_inverse_trig(ctx, integrand, var).expect("mechanism declined");
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        // The pass's fresh proxy symbol must never leak into the answer.
        for leaked in ["w", "k", "r", "z"] {
            assert!(
                !contains_symbol(result, Symbol::new(leaked)),
                "proxy symbol {leaked} leaked: {result}"
            );
        }
        let lhs_expr = replace_atom(ctx, result, kernel, proxy);
        let rhs_expr = replace_atom(ctx, integrand, kernel, proxy);
        let d = crate::diff(ctx, lhs_expr, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(rhs_expr, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
    }

    // ------------------------- M4: inverse-function substitution -------------------------

    /// M4 solve assertion: solved with no `Integral` residue, no `F(G(·))`
    /// composite, and `diff(result) == integrand` numerically.
    fn m4_assert_solved<'a>(
        ctx: &'a AtomArena<'a>,
        input: &str,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) -> Atom<'a> {
        let x = Symbol::new("x");
        let expr = parse_norm(ctx, input);
        let r = integrate_inverse_trig(ctx, expr, x).expect("M4 declined");
        assert!(!r.to_string().contains("Integral"), "{input}: residue: {r}");
        assert!(
            !contains_inv_composite(r),
            "{input}: F(G(·)) composite left: {r}"
        );
        let d = crate::diff(ctx, r, x);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((x, xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(expr, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{input} at x={xv}: diff={lhs} integrand={rhs} (result: {r})"
            );
        }
        r
    }

    #[test]
    fn m4_symbolic_slope_bare_power() {
        // rubi-00554: x²·acos(a·x)² — the symbolic slope makes x² map to
        // `(cos(u)/a)²`, a power of a product; without the flatten_powers
        // step the trig scanner cannot see the factors (regression guard).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 1.4)];
        m4_assert_solved(&ctx, "x^2*acos(a*x)^2", &env, &[0.2, 0.4, 0.6]);
    }

    #[test]
    fn m4_z1_zero_intercept_asin() {
        // §5 Z1: x²·asin(x) → t-form u²·sin²u·cos u (multiple angles).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x^2*asin(x)", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_z2_zero_intercept_asinh() {
        // §5 Z2: x²·asinh(x) — M2 keeps first claim (its caps cover k=1,
        // deg 2); the assertion guards the solved-ness and correctness.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x^2*asinh(x)", &[], &[0.3, 0.8, 1.4]);
    }

    #[test]
    fn m4_z3_zero_intercept_acos() {
        // §5 Z3: x·acos(x) → t-form −u·cos u·sin u.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x*acos(x)", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_k1_kernel_power_trig() {
        // §5 K1: (a + b·asin(x))² → t-form (a + b·u)²·cos u.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 0.7), (Symbol::new("b"), 1.1)];
        m4_assert_solved(&ctx, "(a + b*asin(x))^2", &env, &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_k2_kernel_power_times_poly() {
        // §5 K2 (rubi-01205 numeric shape): x³·(1 + 2·asin(x))².
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x^3*(1 + 2*asin(x))^2", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_k3_kernel_power_hyp() {
        // §5 K3: x·(a + b·asinh(c·x))² — t-form (a+bu)²·sinh u·cosh u/c².
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [
            (Symbol::new("a"), 0.9),
            (Symbol::new("b"), -0.8),
            (Symbol::new("c"), 1.4),
        ];
        m4_assert_solved(&ctx, "x*(a + b*asinh(c*x))^2", &env, &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn m4_t3_bare_power_atan() {
        // §5 T3: x·atan(x)² → t-form u²·tan u·sec²u, finished by the
        // ratio tabular (`I(2,1) = u²tan²u/2 − J(1,2)`).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x*atan(x)^2", &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn m4_b2_bare_power_acos() {
        // §5 B2 (rubi-00554 numeric shape): x²·acos(x)².
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x^2*acos(x)^2", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_h2_m1_first_claim_regression() {
        // §5 H2: atanh(x)²/(1−x²) is M1's kernel shape — M1 claims it
        // before M4 and the answer is atanh(x)³/3.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let r = m4_assert_solved(&ctx, "atanh(x)^2/(1 - x^2)", &[], &[0.2, 0.5, 0.8]);
        assert!(
            r.to_string().contains("(atanh(x))^3"),
            "expected M1's atanh^3/3 shape: {r}"
        );
    }

    #[test]
    fn m4_a1_m2_first_claim_regression() {
        // §5 A1: (a + b·asinh(c + d·x))² is M2's shape — first claim and
        // answer shape stay with M2.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [
            (Symbol::new("a"), 1.1),
            (Symbol::new("b"), 0.7),
            (Symbol::new("c"), 0.4),
            (Symbol::new("d"), 1.2),
        ];
        m4_assert_solved(&ctx, "(a + b*asinh(c + d*x))^2", &env, &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn m4_atah_tabular_weight() {
        // atanh family: x·atanh(x) → t-form u·tanh u·sech²u, and
        // x²·atanh(x) → u·tanh²u·sech²u — both elementary via the
        // I/J tabular (chain emits exp-forms here, so M4 integrates
        // in-module).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x*atanh(x)", &[], &[0.2, 0.5, 0.8]);
        m4_assert_solved(&ctx, "x^2*atanh(x)", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_affine_argument_and_hyp_k5() {
        // Beyond M2's caps: k = 5 kernel power (asinh) and an affine
        // argument with a polynomial weight.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("c"), 1.3)];
        m4_assert_solved(&ctx, "(1 + asinh(c*x))^5", &env, &[0.2, 0.6, 1.0]);
        let env2 = [(Symbol::new("c"), 0.5), (Symbol::new("d"), 0.8)];
        m4_assert_solved(&ctx, "x*asin(c + d*x)", &env2, &[0.1, 0.3, 0.45]);
    }

    #[test]
    fn m4_declines_spec_r_group() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // R1: F in the denominator — t-form cos u/(1+u), non-elementary.
        assert_declined(&ctx, "1/(1 + asin(x))");
        // R2: negative bare power — non-elementary.
        assert_declined(&ctx, "1/acos(x)^3");
        // R3: non-integer kernel power (rubi-01643 family, Fresnel-shaped).
        assert_declined(&ctx, "(1 + asin(x))^(-5/2)");
        // R4 (rubi-00531): polylog emerges in the t-form.
        assert_declined(&ctx, "atanh(1 + x)/(2 + 2*x)");
        // R5 (rubi-01225 family): non-companion rational factor —
        // t-form (a+bu)²·tanh u, polylog. (The radical `(1+x²)^(−1)` IS
        // the asinh companion with m = −2 since W5b; the decline emerges
        // from the t-form's negative cosh power, not the shape gate.)
        assert_declined(&ctx, "x*(a + b*asinh(x))^2/(1 + x^2)");
        // W5b: companion radicals are accepted — CR1/CR2/CR3 below. The
        // non-companion and multi-base radical declines live in
        // `m4_w5b_declines_non_companion_radicals`.
        // u^j·tanh u with j ≥ 1 is polylog: honest decline.
        assert_declined(&ctx, "x^2*(a + b*atanh(c*x))^2/(1 - c^2*x^2)^(-1)");
    }

    #[test]
    fn m4_idempotent_reentry() {
        // §5 W2: the answers re-enter the chain without loops. For the
        // Z1/K1 group every answer term is again elementary and the
        // re-entry must leave no residue. T3's answer provably contains
        // non-elementary answer terms (`atan(x)²`, `x²·atan(x)²` — their
        // t-forms `u²·sec²u` and `u²·tan²u·sec²u` reduce to `∫u·tan u du`,
        // the polylog/Ti₂ territory of §2.5), so for T3 the assertion is
        // termination plus: whatever M4 itself re-claims from the answer
        // terms is complete and correct (no partial guesses).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        for input in ["x^2*asin(x)", "(a + b*asin(x))^2"] {
            let first = m4_assert_solved(
                &ctx,
                input,
                &[(Symbol::new("a"), 0.7), (Symbol::new("b"), 1.1)],
                &[0.3],
            );
            let second = crate::integral::integrate(&ctx, first, x);
            assert!(
                !contains_integral(second),
                "{input}: answer re-entry left a residue: {second}"
            );
        }
        let t3 = m4_assert_solved(&ctx, "x*atan(x)^2", &[], &[0.3]);
        let terms: Vec<Atom> = match t3.node() {
            AtomNode::Add(args) => args.to_vec(),
            _ => vec![t3],
        };
        for term in terms {
            // Termination is the assertion for the full chain; M4's own
            // re-claims must be residue-free and differentiate back.
            if let Some(r) = integrate_inverse_trig(&ctx, term, x) {
                assert!(!contains_integral(r), "partial guess on {term}: {r}");
            }
        }
    }

    #[test]
    fn m4_certificate_gate() {
        // §5: the answers carry a symbolic certificate where the
        // certificate's layers can see the shape.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        for input in ["x^2*asin(x)", "x*atan(x)^2", "x*acos(x)", "x^2*acos(x)^2"] {
            let expr = parse_norm(&ctx, input);
            let r = integrate_inverse_trig(&ctx, expr, x).expect("M4 declined");
            crate::integral::certify::certify(&ctx, expr, r, x)
                .unwrap_or_else(|e| panic!("{input}: certificate declined ({e:?}): {r}"));
        }
        let expr = parse_norm(&ctx, "(a + b*asin(x))^2");
        let r = integrate_inverse_trig(&ctx, expr, x).expect("M4 declined");
        crate::integral::certify::certify(&ctx, expr, r, x)
            .unwrap_or_else(|e| panic!("K1: certificate declined ({e:?}): {r}"));
        // W5b companion-radical answers certify through the radical layer's
        // `r² = E(w)` relation (§4.2/§4.3).
        for input in [
            "x*asin(x)/sqrt(1 - x^2)",
            "x^3*asin(x)/sqrt(1 - x^2)",
            "x*(1 - x^2)^(3/2)*atanh(x)",
        ] {
            let expr = parse_norm(&ctx, input);
            let r = integrate_inverse_trig(&ctx, expr, x).expect("M4 declined");
            crate::integral::certify::certify(&ctx, expr, r, x)
                .unwrap_or_else(|e| panic!("{input}: certificate declined ({e:?}): {r}"));
        }
    }

    // ------------------------- M4/W5b: companion radicals -------------------------

    #[test]
    fn m4_cr1_companion_sqrt_inverse() {
        // §5 CR1: x·asin(x)/√(1−x²) → t-form u·sin u (the dx `cos u` and
        // the radical's `cos⁻¹ u` cancel) → x − √(1−x²)·asin x.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let r = m4_assert_solved(&ctx, "x*asin(x)/sqrt(1 - x^2)", &[], &[0.2, 0.5, 0.8]);
        let rs = r.to_string();
        assert!(
            rs.contains("sqrt") && rs.contains("asin(x)"),
            "expected the √E literal × asin(x) in the answer: {r}"
        );
    }

    #[test]
    fn m4_cr2_companion_with_cubic_weight() {
        // §5 CR2 (rubi-01568 family): x³·asin(x)/√(1−x²) → t-form u·sin³u.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x^3*asin(x)/sqrt(1 - x^2)", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_cr3_atanh_companion_half3() {
        // §5 CR3 (rubi-01733 numeric shape; SymPy declines): the t-form is
        // u·tanh u·sech⁵u — odd sech power, integrated by the `d(S^B)`
        // parts tabular with the `∫sech u du = asin(tanh u)` base, which
        // back-substitutes to `asin(x)`.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let r = m4_assert_solved(&ctx, "x*(1 - x^2)^(3/2)*atanh(x)", &[], &[0.2, 0.5, 0.8]);
        assert!(
            r.to_string().contains("asin(x)"),
            "expected the cross-family asin(w) term: {r}"
        );
    }

    #[test]
    fn m4_w5b_symbolic_companion_scale() {
        // rubi-01568: symbolic normalization s = d (assumed positive) under
        // an odd half power; rubi-01823: k = 2 kernel power × x³ × m = −1.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [
            (Symbol::new("a"), 0.9),
            (Symbol::new("b"), 1.1),
            (Symbol::new("c"), 0.8),
            (Symbol::new("d"), 1.7),
        ];
        m4_assert_solved(
            &ctx,
            "x^3*(d - c^2*d*x^2)^(1/2)*(a + b*asin(c*x))",
            &env,
            &[0.2, 0.5, 0.8],
        );
        m4_assert_solved(
            &ctx,
            "x^3*(a + b*asin(c*x))^2/(d - c^2*d*x^2)^(1/2)",
            &env,
            &[0.2, 0.5, 0.8],
        );
    }

    #[test]
    fn m4_w5b_companion_all_families() {
        // One companion-radical solve per remaining family: acos (sheet
        // and dx signs), asinh/acosh (hyperbolic, m = −1 and m = +1), and
        // atan (odd sec power through the `d(S^B)` tabular).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        m4_assert_solved(&ctx, "x*acos(x)/sqrt(1 - x^2)", &[], &[0.2, 0.5, 0.8]);
        m4_assert_solved(&ctx, "x*asinh(x)/sqrt(1 + x^2)", &[], &[0.3, 0.8, 1.4]);
        m4_assert_solved(&ctx, "x*sqrt(1 + x^2)*asinh(x)", &[], &[0.3, 0.8, 1.4]);
        m4_assert_solved(&ctx, "x*acosh(x)/sqrt(x^2 - 1)", &[], &[1.4, 1.9, 2.4]);
        m4_assert_solved(&ctx, "x*(1 + x^2)^(3/2)*atan(x)", &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn m4_w5b_negative_s_power_conversion() {
        // rubi-00689: atanh companion with m = −5 — the t-form's
        // `sech^(−3)` is the exact reciprocal `cosh³`, so the term
        // `u·tanh·cosh²` converts to the hyperbolic machinery.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 0.9), (Symbol::new("c"), 1.7)];
        m4_assert_solved(
            &ctx,
            "atanh(a*x)/(c - a^2*c*x^2)^(5/2)",
            &env,
            &[0.2, 0.5, 0.8],
        );
        // m = −3 with a weight: `u·cosh u` after the cancellation.
        m4_assert_solved(&ctx, "atanh(x)/(1 - x^2)^(3/2)", &[], &[0.2, 0.5, 0.8]);
    }

    #[test]
    fn m4_w5b_declines_non_companion_radicals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Companion exponent but an unrelated base (`d + e·x²` is not
        // `s·(1 − x²)` for symbolic e, d).
        assert_declined(&ctx, "x*asin(x)*sqrt(d + e*x^2)");
        // Two radical bases: one companion, one not — declined.
        assert_declined(&ctx, "x*asin(x)*sqrt(1 - x^2)*sqrt(1 + x^2)");
        // Imaginary-carrying radicand (rubi-00408 family): no `I`.
        assert_declined(&ctx, "x*asinh(x)*sqrt(i + i*x^2)");
        // Companion with m = −2: the combined head power is negative
        // (`u·sin u·sec u` — polylog territory), honest decline.
        assert_declined(&ctx, "x*asin(x)/(1 - x^2)");
        // Odd half power with a negative numeric normalization: no real
        // square root of s.
        assert_declined(&ctx, "x*asinh(x)*sqrt(-4 - 4*x^2)");
        // Symbolic exponent on the radicand: undecidable, decline (§2.5-5).
        assert_declined(&ctx, "x*asin(x)*(1 - x^2)^(n/2)");
    }

    #[test]
    fn m4_w5b_idempotent_reentry() {
        // §5 W2 for the W5b answers: re-entering them terminates with no
        // residue introduced by M4's own re-claims.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        for input in ["x*asin(x)/sqrt(1 - x^2)", "x*(1 - x^2)^(3/2)*atanh(x)"] {
            let first = m4_assert_solved(&ctx, input, &[], &[0.3]);
            let terms: Vec<Atom> = match first.node() {
                AtomNode::Add(args) => args.to_vec(),
                _ => vec![first],
            };
            for term in terms {
                // Termination is the assertion for the full chain; M4's own
                // re-claims must be complete and correct (no partial guesses).
                if let Some(r) = integrate_inverse_trig(&ctx, term, x) {
                    assert!(
                        !contains_integral(r),
                        "{input}: partial guess on {term}: {r}"
                    );
                }
            }
        }
    }

    // ------------------------- M1 -------------------------

    /// Simpson's rule on `[a, b]` for the tabular self-check below.
    fn simpson<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
        let h = (b - a) / n as f64;
        let mut acc = f(a) + f(b);
        for i in 1..n {
            let x = a + h * i as f64;
            acc += if i % 2 == 1 { 4.0 * f(x) } else { 2.0 * f(x) };
        }
        acc * h / 3.0
    }

    /// Verify the parameterised ratio tabular of all four families against
    /// numerical quadrature of `u^j·R^a·S^b` on a pole-free window. This is
    /// the guard that makes the sign algebra (`ε`, `δ`, `η`, `γ`) of the two
    /// Wave 5c mirror families trustworthy: every `(j, a, b)` combination the
    /// tabular accepts must reproduce the quadrature to `1e-6`.
    ///
    /// (It is also what caught the pre-existing `ratio_s_pow(j ≥ 1, b ≥ 3)`
    /// reduction emitting a non-antiderivative; see the helper's docs.)
    #[test]
    fn m4_w5c_ratio_tabular_numeric() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let t_sym = Symbol::new("t");
        let t = ctx.var("t");
        let fams = [
            ("tan/sec", RatioFam::TanSec),
            ("tanh/sech", RatioFam::TanhSech),
            ("cot/csc", RatioFam::CotCsc),
            ("coth/csch", RatioFam::CothCsch),
        ];
        let (x0, x1) = (0.7f64, 1.3f64);
        let mut accepted = 0usize;
        for (name, fam) in fams {
            let (r_name, s_name) = fam.heads();
            for j in 0..=3i64 {
                for a in 0..=4i64 {
                    for b in 0..=5i64 {
                        let Some(cand) = ratio_rs(&ctx, j, a, b, t, fam) else {
                            continue;
                        };
                        accepted += 1;
                        let integrand = |u: f64| {
                            let e = vec![(t_sym, u)];
                            let rv = eval_f64(ctx.fun(r_name, &[t]), &e).unwrap();
                            let sv = eval_f64(ctx.fun(s_name, &[t]), &e).unwrap();
                            u.powi(j as i32) * rv.powi(a as i32) * sv.powi(b as i32)
                        };
                        let numeric = simpson(integrand, x0, x1, 2000);
                        let f0 = eval_f64(cand, &[(t_sym, x0)]).unwrap();
                        let f1 = eval_f64(cand, &[(t_sym, x1)]).unwrap();
                        assert!(
                            f0.is_finite() && f1.is_finite(),
                            "{name}: non-finite {cand}"
                        );
                        let got = f1 - f0;
                        let tol = 1e-6 * numeric.abs().max(1.0);
                        assert!(
                            (got - numeric).abs() < tol,
                            "{name}: ∫u^{j}·{r_name}^{a}·{s_name}^{b} on [{x0},{x1}]: \
                             tabular {got} vs numeric {numeric} (cand {cand})"
                        );
                    }
                }
            }
        }
        // Sanity: the grid must actually exercise the tabular (a silently
        // declining tabular would pass every assertion above).
        assert!(accepted > 150, "only {accepted} tabular cases accepted");
    }

    // ------------------------- M4 / W5c: reciprocal heads -------------------------

    /// `(a, b, c, x-samples)`: the constants and the sample window of one
    /// reciprocal-head domain profile.
    type W5cProfile = (f64, f64, f64, [f64; 3]);
    /// One probe row: `(id, input, profile, solved?)`.
    type W5cProbe = (&'static str, &'static str, W5cProfile, bool);

    /// The three domain profiles the W5c solve assertions sample: the
    /// reciprocal heads have disjoint real domains, so one shared sample
    /// window cannot exercise them (`|w| > 1` for acsc/asec/acoth, `0 < w < 1`
    /// for asech, `w > 0` for acsch).
    const W5C_GT1: W5cProfile = (1.2, 0.8, 1.6, [0.8, 1.2, 1.6]);
    const W5C_LT1: W5cProfile = (0.2, 0.4, 0.4, [0.5, 0.9, 1.4]);
    const W5C_POS: W5cProfile = (0.5, 0.4, 1.2, [0.4, 0.8, 1.2]);

    /// `(a, b, c)` constants of a profile plus the `x` samples.
    fn w5c_env(profile: W5cProfile) -> ([(Symbol, f64); 3], [f64; 3]) {
        (
            [
                (Symbol::new("a"), profile.0),
                (Symbol::new("b"), profile.1),
                (Symbol::new("c"), profile.2),
            ],
            profile.3,
        )
    }

    /// Numeric-checked W5c solve assertion: solved with no `Integral`
    /// residue, no `F(G(·))` composite and no leaked substitution symbol, and
    /// `diff(result) == integrand` at every sample. Unlike the W5a/W5b helper
    /// this one *rejects* non-finite samples — a NaN comparison would compare
    /// equal to nothing and silently pass a sample outside the head's domain.
    fn m4_w5c_solved<'a>(
        ctx: &'a AtomArena<'a>,
        input: &str,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) -> Atom<'a> {
        let x = Symbol::new("x");
        let expr = parse_norm(ctx, input);
        let r = integrate_inverse_trig(ctx, expr, x).expect("W5c declined");
        assert!(!r.to_string().contains("Integral"), "{input}: residue: {r}");
        assert!(
            !contains_inv_composite(r),
            "{input}: F(G(·)) composite left: {r}"
        );
        for leaked in ["t", "u"] {
            assert!(
                !contains_symbol(r, Symbol::new(leaked)),
                "{input}: substitution symbol {leaked} leaked: {r}"
            );
        }
        let d = crate::diff(ctx, r, x);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((x, xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(expr, &env).expect("eval integrand");
            assert!(
                lhs.is_finite() && rhs.is_finite(),
                "{input} at x={xv}: sample outside the head's real domain: \
                 diff={lhs} integrand={rhs} (result: {r})"
            );
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{input} at x={xv}: diff={lhs} integrand={rhs} (result: {r})"
            );
        }
        r
    }

    /// Symbolic certificate gate for one W5c input.
    fn m4_w5c_certifies<'a>(ctx: &'a AtomArena<'a>, input: &str) {
        let x = Symbol::new("x");
        let expr = parse_norm(ctx, input);
        let r = integrate_inverse_trig(ctx, expr, x).expect("W5c declined");
        crate::integral::certify::certify(ctx, expr, r, x)
            .unwrap_or_else(|e| panic!("{input}: certificate declined ({e:?}): {r}"));
    }

    #[test]
    fn m4_w5c_bare_and_affine_heads() {
        // Every reciprocal head solves its own bare and affine-argument
        // shapes: the substitution `u = F(ρ + σx)` turns the integrand into
        // `u·G′(u)/σ`, and the new ratio families integrate it back into the
        // textbook `w·F(w) ± log(w + √E(w))` / `w·F(w) ± acos(w)` forms
        // (spec §2.1 table rows acsc/asec/acoth/asech/acsch).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Bare forms `F(x)`: the domain is the head's own.
        let bare: &[(&str, &[f64])] = &[
            ("acsc(x)", &[1.2, 1.8, 2.4]),
            ("asec(x)", &[1.2, 1.8, 2.4]),
            ("acoth(x)", &[1.2, 1.8, 2.4]),
            ("asech(x)", &[0.2, 0.5, 0.8]),
            ("acsch(x)", &[0.2, 0.6, 1.2]),
        ];
        for (input, samples) in bare {
            let r = m4_w5c_solved(&ctx, input, &[], samples);
            let head = &input[..input.find('(').unwrap()];
            assert!(
                r.to_string().contains(&format!("{head}(")),
                "{input}: answer lost its own head: {r}"
            );
        }
        // Affine forms `F(a + b·x)`, each with the domain profile its head
        // needs.
        let affine: &[(&str, W5cProfile)] = &[
            ("acsc(a + b*x)", W5C_GT1),
            ("asec(a + b*x)", W5C_GT1),
            ("acoth(a + b*x)", W5C_GT1),
            ("asech(a + b*x)", W5C_LT1),
            ("acsch(a + b*x)", W5C_POS),
        ];
        for (input, profile) in affine {
            let (env, samples) = w5c_env(*profile);
            let r = m4_w5c_solved(&ctx, input, &env, &samples);
            let head = &input[..input.find('(').unwrap()];
            assert!(
                r.to_string().contains(&format!("{head}(")),
                "{input}: answer lost its own head: {r}"
            );
        }
    }

    #[test]
    fn m4_w5c_weighted_heads_and_bucket_rows() {
        // Polynomial weights on every reciprocal head (k = 1 kernel power),
        // plus the two bucket rows the 0.30.0 attribution puts in W5c that
        // are inside the caps: rubi-00511 `x^3·acsc(a + b·x)` and
        // rubi-01708 `x^4·asec(a + b·x)` (spec §5 row A4).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let heads: &[(&str, W5cProfile)] = &[
            ("x^2*(a + b*acsc(c*x))", W5C_GT1),
            ("x^2*(a + b*asec(c*x))", W5C_GT1),
            ("x^2*(a + b*acoth(c*x))", W5C_GT1),
            ("x*(a + b*acoth(c*x))", W5C_GT1),
            ("x^2*(a + b*asech(c*x))", W5C_LT1),
            ("x*(a + b*asech(c*x))", W5C_LT1),
            ("x^2*(a + b*acsch(c*x))", W5C_POS),
            ("x*(a + b*acsch(c*x))", W5C_POS),
        ];
        for (input, profile) in heads {
            let (env, samples) = w5c_env(*profile);
            m4_w5c_solved(&ctx, input, &env, &samples);
        }
        // rubi-00511 / rubi-01708: affine argument with an intercept, so the
        // substitution carries the `(G(u) − ρ)/σ` numerator.
        let env = [(Symbol::new("a"), 0.9), (Symbol::new("b"), 1.1)];
        m4_w5c_solved(&ctx, "x^3*acsc(a + b*x)", &env, &[0.2, 0.5, 0.8]);
        m4_w5c_solved(&ctx, "x^4*asec(a + b*x)", &env, &[0.2, 0.5, 0.8]);
        // k = 2 with a linear weight also closes for asec/asec/acsc/acsch
        // (the `u^2` term rides on `R·S²`-shaped terms, which the tabular
        // finishes); the bare `(a + b·F)^2` forms decline (see the declines
        // test) because their `u^2` terms land on the polylog base.
        for (input, profile) in [
            ("x*(a + b*acsc(c*x))^2", W5C_GT1),
            ("x*(a + b*asec(c*x))^2", W5C_GT1),
            ("x*(a + b*asech(c*x))^2", W5C_LT1),
            ("x*(a + b*acsch(c*x))^2", W5C_POS),
        ] {
            let (env, samples) = w5c_env(profile);
            m4_w5c_solved(&ctx, input, &env, &samples);
        }
    }

    #[test]
    fn m4_w5c_companion_radicals() {
        // W5b's companion absorption reaches the reciprocal heads through
        // `family_kernel_head`: `(s·E(w))^(m/2) = s^(m/2)·kernel_head^m`, and
        // for these five the kernel head is one of the two `G′` factors
        // (cot/tan/csch/tanh/coth), so the radical merges into that factor's
        // power. Radicands: `σ²x²−1` for acsc/asec/acoth, `1−σ²x²` for asech,
        // `1+σ²x²` for acsch.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let cases: &[(&str, W5cProfile)] = &[
            ("x*acsc(c*x)*sqrt(c^2*x^2 - 1)", W5C_GT1),
            ("x*asec(c*x)*sqrt(c^2*x^2 - 1)", W5C_GT1),
            ("x*acoth(c*x)*sqrt(c^2*x^2 - 1)", W5C_GT1),
            ("x*asech(c*x)*sqrt(1 - c^2*x^2)", W5C_LT1),
            ("x*acsch(c*x)*sqrt(1 + c^2*x^2)", W5C_POS),
        ];
        for (input, profile) in cases {
            let (env, samples) = w5c_env(*profile);
            m4_w5c_solved(&ctx, input, &env, &samples);
        }
    }

    #[test]
    fn m4_w5c_declines_reciprocal_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // The 13 reciprocal-head rows of the 0.30.0 attribution plus the
        // per-head narrowing guards. Every decline here is a *shape* gate or
        // the polylog base case, never a wrong closed form.
        //
        // (a) Symbolic exponent on the weight or on the whole product
        // (`(f·x)^m`, `x^n`, `(d·x)^m`): undecidable, decline (§2.5-5).
        assert_declined(&ctx, "(f*x)^m*(a + b*asec(c*x))"); // 00441
        assert_declined(&ctx, "(d + e*x)^m*(a + b*asech(c*x))"); // 01317
        assert_declined(&ctx, "(d*x)^m*(a + b*asec(c*x))"); // 01786
        assert_declined(&ctx, "x^n*(a + b*acsc(c*x))");
        // (b) Negative powers of `x` / `F` in the denominator (F-in-denom
        // and reciprocal-power families, §2.5-2).
        assert_declined(&ctx, "acoth(a + b*x)/x^3"); // 00026
        assert_declined(&ctx, "(d + e*x^2)*(a + b*acsch(c*x))/x^6"); // 01376
        assert_declined(&ctx, "x^4*(a + b*asech(c*x))/(d + e*x^2)^3"); // 00635
        assert_declined(&ctx, "(a + b*acsch(c + d*x))^2/(e + f*x)^3"); // 00781
        // (c) Non-companion radicals: `d + e·x²` / `d + e·x` is not the
        // family kernel `σ²x² − 1` for symbolic coefficients (§2.5-4).
        assert_declined(&ctx, "x*(d + e*x^2)^(3/2)*(a + b*acsc(c*x))"); // 00091
        assert_declined(&ctx, "(d + e*x)^(3/2)*(a + b*acsc(c*x))"); // 00821
        assert_declined(&ctx, "x*(a + b*acsch(c*x))/(d + e*x)^(3/2)"); // 00544
        assert_declined(&ctx, "(a + b*acsch(c*x))/(x^3*(d + e*x^2)^(5/2))"); // 01304
        assert_declined(&ctx, "x^2*(a + b*acsc(c*x))/(d + e*x^2)^(5/2)"); // 01351
        assert_declined(&ctx, "(f*x)^m*(a + b*asec(c*x))/(d + e*x^2)^(3/2)"); // 01793
        // (d) Polynomial weight above `MAX_IS_POLY_DEG = 6`: total degree 7.
        assert_declined(&ctx, "x^5*(d + e*x^2)*(a + b*acsc(c*x))"); // 00167
        // (e) Non-affine argument (`f^(c + d·x)`): `linear_form` declines.
        assert_declined(&ctx, "acoth(a + b*f^((c + d*x)))"); // 01681
        // (f) Polylog base case: `∫u^j·S(u) du` / `∫u^j·S^odd·…` with
        // `j ≥ 1`. These are the reciprocal-head kernel powers whose `u^j`
        // terms land on the non-elementary base, so the tabular declines
        // (the residue-check analogue of §2.3-3).
        assert_declined(&ctx, "(a + b*acoth(c*x))^2");
        assert_declined(&ctx, "(a + b*acoth(c + d*x))^3"); // 00372
        assert_declined(&ctx, "x^2*asech(a*x)^3"); // 01347
        assert_declined(&ctx, "x^3*(a + b*acsc(c*x))^3"); // 01132
        assert_declined(&ctx, "x*(a + b*asec(c*x))^3"); // 01218
        assert_declined(&ctx, "x^2*(a + b*acsch(c*x))^2"); // 01467
        assert_declined(&ctx, "(a + b*asec(c*x))^2/(c^2*x^2 - 1)");
        // (g) The new heads must stay out of M1/P3 (their `F′` carries the
        // reciprocal forward head, so the kernel-derivative power rule does
        // not apply): a companion-shaped kernel factor in the denominator.
        assert_declined(&ctx, "x*(a + b*asec(c*x))^2/(c^2*x^2 - 1)");
        assert_declined(&ctx, "x*(a + b*acsch(c*x))^2/(1 + c^2*x^2)");
    }

    #[test]
    fn m4_w5c_suppressed_inside_parts_scope() {
        // Spec §3.1: inside an integration-by-parts scope M4 is all or
        // nothing and returns `None` outright — the substitution must never
        // upgrade a parts residue into a complete transcendental answer.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let cases: &[(&str, &[f64])] = &[
            ("acsc(x)", &[1.2, 1.8, 2.4]),
            ("asech(x)", &[0.2, 0.5, 0.8]),
            ("x^2*(a + b*asec(c*x))", &[0.8, 1.2, 1.6]),
        ];
        for (input, samples) in cases {
            let expr = parse_norm(&ctx, input);
            let inside =
                super::super::chain::with_parts_scope(|| integrate_inverse_trig(&ctx, expr, x));
            assert!(
                inside.is_none(),
                "{input}: emitted inside a parts scope: {inside:?}"
            );
            // Sanity: outside the scope the same input solves.
            m4_w5c_solved(&ctx, input, &w5c_env(W5C_GT1).0, samples);
        }
    }

    #[test]
    fn m4_w5c_certificate_gate() {
        // Every W5c answer that the certificate's layers can see is emitted
        // through them (spec §1.3-4): the bare/affine/weighted shapes certify
        // in the elementary field or through the radical layer's `r² = E(w)`
        // relation.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for input in [
            "acsc(x)",
            "asec(x)",
            "acoth(x)",
            "asech(x)",
            "acsch(x)",
            "acsc(a + b*x)",
            "asec(a + b*x)",
            "acoth(a + b*x)",
            "asech(a + b*x)",
            "acsch(a + b*x)",
            "x^2*(a + b*acsc(c*x))",
            "x^2*(a + b*asec(c*x))",
            "x*(a + b*acoth(c*x))",
            "x^2*(a + b*asech(c*x))",
            "x^2*(a + b*acsch(c*x))",
            "x^3*acsc(a + b*x)",
            "x^4*asec(a + b*x)",
        ] {
            m4_w5c_certifies(&ctx, input);
        }
        // The companion-radical answers certify once the slope is numeric
        // (with a *symbolic* slope the certificate declines on these shapes
        // already for the W5b families — `x*asin(c*x)/sqrt(1 - c^2*x^2)` is
        // the same NonZero decline — so this is a pre-existing certificate
        // gap of the radical layer, not a W5c regression).
        for input in [
            "x*asec(x)*sqrt(x^2 - 1)",
            "x*acsc(x)*sqrt(x^2 - 1)",
            "x*acsch(x)*sqrt(1 + x^2)",
            "x*asech(x)*sqrt(1 - x^2)",
            "x*acoth(x)*sqrt(x^2 - 1)",
        ] {
            m4_w5c_certifies(&ctx, input);
        }
    }

    #[test]
    fn m4_w5c_idempotent_reentry() {
        // Spec §1.3-5 / §5 W2. The answers re-enter the chain without loops
        // and M4 never re-claims a shape it cannot finish. The full chain may
        // still leave a residue on the answers' *own* elementary pieces
        // (`∫log(w + √E(w)) dw`, `∫√((a+bx)² − 1) dx` — outside this
        // module's shape space and already declined before W5c), so, exactly
        // like the W5a T3 case, the assertion is: re-entry terminates, and
        // whatever M4 itself re-claims from an answer term is complete and
        // correct (no partial guesses, no residue).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let cases: &[(&str, W5cProfile)] = &[
            ("acsc(x)", (1.2, 0.8, 1.0, [1.2, 1.8, 2.4])),
            ("asec(x)", (1.2, 0.8, 1.0, [1.2, 1.8, 2.4])),
            ("acoth(x)", (1.2, 0.8, 1.0, [1.2, 1.8, 2.4])),
            ("asech(x)", (0.2, 0.4, 0.4, [0.2, 0.5, 0.8])),
            ("acsch(x)", W5C_POS),
            ("x^2*(a + b*acsc(c*x))", W5C_GT1),
            ("x^2*(a + b*asec(c*x))", W5C_GT1),
            ("x*(a + b*acoth(c*x))", W5C_GT1),
            ("x^2*(a + b*asech(c*x))", W5C_LT1),
            ("x^2*(a + b*acsch(c*x))", W5C_POS),
        ];
        for (input, profile) in cases {
            let (env, samples) = w5c_env(*profile);
            let first = m4_w5c_solved(&ctx, input, &env, &samples);
            // Termination + completeness of the mechanism's own re-claims.
            let terms: Vec<Atom> = match first.node() {
                AtomNode::Add(args) => args.to_vec(),
                _ => vec![first],
            };
            for term in terms {
                if let Some(r) = integrate_inverse_trig(&ctx, term, x) {
                    assert!(
                        !contains_integral(r),
                        "{input}: partial guess on {term}: {r}"
                    );
                    assert!(
                        !contains_inv_composite(r),
                        "{input}: composite on {term}: {r}"
                    );
                    // The re-claim must still be an antiderivative.
                    let d = crate::diff(&ctx, r, x);
                    for &xv in &samples {
                        let mut e = env.to_vec();
                        e.push((x, xv));
                        let (Some(lhs), Some(rhs)) = (eval_f64(d, &e), eval_f64(term, &e)) else {
                            continue;
                        };
                        if !lhs.is_finite() || !rhs.is_finite() {
                            continue;
                        }
                        assert!(
                            (lhs - rhs).abs() < 1e-6 * rhs.abs().max(1.0),
                            "{input}: re-claim of {term} is not an antiderivative: {r}"
                        );
                    }
                }
            }
            // The full chain terminates on the answer (no hang, no loop).
            let second = crate::integral::integrate(&ctx, first, x);
            assert!(
                !contains_symbol(second, Symbol::new("t")),
                "{input}: substitution symbol leaked through re-entry"
            );
        }
    }

    /// W5c probe over the reciprocal-head rows of
    /// `ocas-tests/data/invtrig_attribution_030.csv` plus the per-head
    /// synthetic shapes (run with `--nocapture` and `OCAS_PROBE=1` for the
    /// table). Solved rows are verified by numeric differentiation at their
    /// domain samples; every row asserts the refusal path's cost stays well
    /// under the 2 s debug budget (spec §8 risk 7: refusing fast matters as
    /// much as solving right).
    #[test]
    fn m4_w5c_bucket_probe_lock() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // (rubi id, input, domain profile, solved?)
        let csv = |c: f64| (0.9, 1.1, c, [0.2, 0.5, 0.8]);
        let probes: &[W5cProbe] = &[
            ("00026", "acoth(a + b*x)/x^3", csv(0.8), false),
            (
                "00091",
                "x*(d + e*x^2)^(3/2)*(a + b*acsc(c*x))",
                csv(0.8),
                false,
            ),
            (
                "00167",
                "x^5*(d + e*x^2)*(a + b*acsc(c*x))",
                csv(0.8),
                false,
            ),
            ("00372", "(a + b*acoth(c + d*x))^3", csv(0.8), false),
            (
                "00441",
                "(f*x)^m*(a + b*asec(c*x))/(d + e*x^2)^2",
                csv(0.8),
                false,
            ),
            ("00511", "x^3*acsc(a + b*x)", csv(0.8), true),
            (
                "00544",
                "x*(a + b*acsch(c*x))/(d + e*x)^(3/2)",
                csv(0.8),
                false,
            ),
            (
                "00635",
                "x^4*(a + b*asech(c*x))/(d + e*x^2)^3",
                csv(0.8),
                false,
            ),
            (
                "00781",
                "(a + b*acsch(c + d*x))^2/(e + f*x)^3",
                csv(0.8),
                false,
            ),
            (
                "00821",
                "(d + e*x)^(3/2)*(a + b*acsc(c*x))",
                csv(0.8),
                false,
            ),
            ("01132", "x^3*(a + b*acsc(c*x))^3", csv(0.8), false),
            ("01218", "x*(a + b*asec(c*x))^3", csv(0.8), false),
            (
                "01304",
                "(a + b*acsch(c*x))/(x^3*(d + e*x^2)^(5/2))",
                csv(0.8),
                false,
            ),
            ("01317", "(d + e*x)^m*(a + b*asech(c*x))", csv(0.8), false),
            ("01347", "x^2*asech(a*x)^3", W5C_LT1, false),
            (
                "01351",
                "x^2*(a + b*acsc(c*x))/(d + e*x^2)^(5/2)",
                csv(0.8),
                false,
            ),
            (
                "01376",
                "(d + e*x^2)*(a + b*acsch(c*x))/x^6",
                csv(0.8),
                false,
            ),
            ("01467", "x^2*(a + b*acsch(c*x))^2", W5C_POS, false),
            ("01681", "acoth(a + b*f^((c + d*x)))", csv(0.8), false),
            ("01708", "x^4*asec(a + b*x)", csv(0.8), true),
            ("01786", "(d*x)^m*(a + b*asec(c*x))", csv(0.8), false),
            (
                "01793",
                "(f*x)^m*(a + b*asec(c*x))/(d + e*x^2)^(3/2)",
                csv(0.8),
                false,
            ),
            // Per-head synthetic shapes inside the caps (solved).
            ("w5c-acsc", "x^2*(a + b*acsc(c*x))", W5C_GT1, true),
            ("w5c-asec", "x^2*(a + b*asec(c*x))", W5C_GT1, true),
            ("w5c-acoth", "x*(a + b*acoth(c*x))", W5C_GT1, true),
            ("w5c-acoth2", "x^2*(a + b*acoth(c*x))", W5C_GT1, true),
            ("w5c-asech", "x^2*(a + b*asech(c*x))", W5C_LT1, true),
            ("w5c-acsch", "x^2*(a + b*acsch(c*x))", W5C_POS, true),
            ("w5c-acsc-k2", "x*(a + b*acsc(c*x))^2", W5C_GT1, true),
            ("w5c-asec-k2", "x*(a + b*asec(c*x))^2", W5C_GT1, true),
            ("w5c-asech-k2", "x*(a + b*asech(c*x))^2", W5C_LT1, true),
            ("w5c-acsch-k2", "x*(a + b*acsch(c*x))^2", W5C_POS, true),
            // Narrowing guards beside the rows.
            (
                "w5c-asec-den",
                "(a + b*asec(c*x))^2/(c^2*x^2 - 1)",
                W5C_GT1,
                false,
            ),
            ("w5c-acsc-sympow", "x^n*(a + b*acsc(c*x))", W5C_GT1, false),
            ("w5c-asech-bare-k3", "x^2*asech(a*x)^3", W5C_LT1, false),
        ];
        for (id, input, profile, expect_solved) in probes {
            let expr = parse_norm(&ctx, input);
            let t0 = std::time::Instant::now();
            let r = integrate_inverse_trig(&ctx, expr, x);
            let dt = t0.elapsed();
            assert!(
                dt.as_secs_f64() < 2.0,
                "{id}: M4 took {dt:?} (refusal budget is 2 s in debug)"
            );
            let solved = r.is_some_and(|r| !contains_integral(r));
            if std::env::var_os("OCAS_PROBE").is_some() {
                let shape = match &r {
                    Some(a) => a.to_string(),
                    None => "<declined>".to_string(),
                };
                eprintln!("[probe] rubi-{id}: solved={solved} {dt:?} {input} => {shape}");
            }
            assert_eq!(solved, *expect_solved, "rubi-{id}: {input}");
            if solved {
                // The pipeline itself (not just this module) must return the
                // answer: nothing downstream may eat M4's result.
                let e2e = crate::integral::integrate(&ctx, expr, x);
                assert!(
                    !contains_integral(e2e),
                    "rubi-{id}: full chain left a residue: {e2e}"
                );
            }
            if let Some(ans) = r {
                assert!(!contains_inv_composite(ans), "rubi-{id}: composite: {ans}");
                let (env, samples) = w5c_env(*profile);
                let d = crate::diff(&ctx, ans, x);
                for &xv in &samples {
                    let mut e = env.to_vec();
                    e.push((x, xv));
                    let lhs = eval_f64(d, &e).expect("eval diff");
                    let rhs = eval_f64(expr, &e).expect("eval integrand");
                    assert!(
                        lhs.is_finite() && rhs.is_finite(),
                        "rubi-{id} at x={xv}: non-finite sample"
                    );
                    let tol = 1e-6 * rhs.abs().max(1.0);
                    assert!(
                        (lhs - rhs).abs() < tol,
                        "rubi-{id} at x={xv}: diff={lhs} integrand={rhs} (result: {ans})"
                    );
                }
            }
        }
    }

    /// W5b probe over the companion-radical rows of
    /// `ocas-tests/data/invtrig_attribution_030.csv` (run with
    /// `--nocapture` and `OCAS_PROBE=1` to see the table). Solved rows are
    /// verified by numeric differentiation at the samples; the assertions
    /// lock the probed outcomes as regression guards.
    #[test]
    fn m4_w5b_bucket_probe_lock() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let env = [
            (Symbol::new("a"), 0.9),
            (Symbol::new("b"), 1.1),
            (Symbol::new("c"), 0.8),
            (Symbol::new("d"), 1.7),
            (Symbol::new("f"), 1.3),
            (Symbol::new("g"), 0.7),
        ];
        let samples = [0.2, 0.5, 0.8];
        // (rubi id, input, solved?)
        let probes: &[(&str, &str, bool)] = &[
            ("01568", "x^3*(d - c^2*d*x^2)^(1/2)*(a + b*asin(c*x))", true),
            (
                "01823",
                "x^3*(a + b*asin(c*x))^2/(d - c^2*d*x^2)^(1/2)",
                true,
            ),
            ("01733", "x*(1 - a^2*x^2)^(3/2)*atanh(a*x)", true),
            ("00689", "atanh(a*x)/(c - a^2*c*x^2)^(5/2)", true),
            (
                "00320",
                "(f + g*x)^2*(d - c^2*d*x^2)^(3/2)*(a + b*asin(c*x))^2",
                true,
            ),
            (
                "00459",
                "x^5*(a + b*asin(c*x))/(d - c^2*d*x^2)^(3/2)",
                false,
            ),
            (
                "00228",
                "x^3*(a + b*asin(c*x))/(d - c^2*d*x^2)^(5/2)",
                false,
            ),
            ("00811", "asin(a*x)^2/(c - a^2*c*x^2)^(3/2)", false),
            (
                "01666",
                "x^2*(d - c^2*d*x^2)^(3/2)*(a + b*asin(c*x))^n",
                false,
            ),
            ("00610", "(1 - a^2*x^2)^(1/2)/atanh(a*x)^2", false),
            ("01571", "(c - a^2*c*x^2)^(5/2)/acosh(a*x)^(3/2)", false),
            (
                "00408",
                "(f - i*c*f*x)^(5/2)*(a + b*asinh(c*x))/(d + i*c*d*x)^(5/2)",
                false,
            ),
        ];
        for (id, input, expect_solved) in probes {
            let expr = parse_norm(&ctx, input);
            let r = integrate_inverse_trig(&ctx, expr, x);
            let solved = r.is_some_and(|r| !contains_integral(r));
            if std::env::var_os("OCAS_PROBE").is_some() {
                let shape = match &r {
                    Some(a) => a.to_string(),
                    None => "<declined>".to_string(),
                };
                eprintln!("[probe] rubi-{id}: solved={solved} {input} => {shape}");
            }
            assert_eq!(solved, *expect_solved, "rubi-{id}: {input}");
            if let Some(ans) = r {
                // Numeric diff verification (the same discipline as the
                // m4_assert_solved tests).
                assert!(!contains_inv_composite(ans), "rubi-{id}: composite: {ans}");
                let d = crate::diff(&ctx, ans, x);
                for &xv in &samples {
                    let mut e = env.to_vec();
                    e.push((x, xv));
                    let lhs = eval_f64(d, &e).expect("eval diff");
                    let rhs = eval_f64(expr, &e).expect("eval integrand");
                    let tol = 1e-6 * rhs.abs().max(1.0);
                    assert!(
                        (lhs - rhs).abs() < tol,
                        "rubi-{id} at x={xv}: diff={lhs} integrand={rhs} (result: {ans})"
                    );
                }
            }
        }
    }

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
        // k = 6 exceeds the M2 cap — but 0.30.0 W5a's M4 accepts it
        // (k ≤ MAX_IS_K = 6), so the shape is now solved.
        {
            let expr = parse_norm(&ctx, "(d + e*x)*(a + b*acosh(c*x))^6");
            let env = [
                (Symbol::new("a"), 0.6),
                (Symbol::new("b"), 1.1),
                (Symbol::new("c"), 0.8),
                (Symbol::new("d"), 1.3),
                (Symbol::new("e"), 0.7),
            ];
            // acosh domain: c·x > 1.
            assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[1.6, 2.2, 3.0]);
        }
        // Numeric negative normalization under an odd square-root power.
        assert_declined(&ctx, "(-4 - 4*c^2*x^2)^(3/2)*(a + b*asinh(c*x))");
        // asinh needs E = 1 + u²; the 1 − c²x² radicand is the wrong kernel.
        assert_declined(&ctx, "(a + b*asinh(c*x))^2/sqrt(1 - c^2*x^2)");
        // atanh is outside M2's asinh/acosh substitution scope.
        assert_declined(&ctx, "(1 + c*x)^3*(a + b*atanh(c*x))^3");
    }

    // ------------------------- P1: composition cancellation -------------------------

    #[test]
    fn p1_atanh_tanh_numeric_and_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // Numeric coefficient.
        let expr = parse_norm(&ctx, "atanh(tanh(b*x))^2");
        let env = [(Symbol::new("b"), 1.2)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.1, 0.4, 0.9]);
        // Symbolic coefficients (corpus rubi-01885).
        let expr = parse_norm(&ctx, "atanh(tanh(a + b*x))^2");
        let env = [(Symbol::new("a"), 0.4), (Symbol::new("b"), 1.1)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.1, 0.3, 0.5]);
    }

    #[test]
    fn p1_corpus_atanh_tanh_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let env = [(Symbol::new("a"), 0.4), (Symbol::new("b"), 0.9)];
        // rubi-01578, rubi-01885.
        for input in ["atanh(tanh(a + b*x))^3", "atanh(tanh(a + b*x))^2"] {
            let expr = parse_norm(&ctx, input);
            assert_antiderivative_num(&ctx, expr, x, &env, &[0.1, 0.3, 0.5]);
        }
        // rubi-00090.
        let expr = parse_norm(&ctx, "atanh(tanh(a + b*x))^4/x^4");
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.3, 0.6, 1.1]);
        // rubi-00177.
        let expr = parse_norm(&ctx, "atanh(tanh(a + b*x))^2*sqrt(x)");
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.2, 0.5, 0.9]);
    }

    #[test]
    fn p1_acoth_coth_numeric_and_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // Numeric coefficient; samples keep |u| bounded and u ≠ 0.
        let expr = parse_norm(&ctx, "acoth(coth(b*x))^2");
        let env = [(Symbol::new("b"), 1.1)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.3, 0.7, 1.1]);
        // Symbolic coefficients.
        let expr = parse_norm(&ctx, "acoth(coth(a + b*x))^3");
        let env = [(Symbol::new("a"), 0.5), (Symbol::new("b"), 1.1)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.3, 0.6, 0.9]);
    }

    #[test]
    fn p1_asinh_sinh_numeric_and_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let expr = parse_norm(&ctx, "asinh(sinh(b*x))^3");
        let env = [(Symbol::new("b"), 1.3)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.2, 0.6, 1.0]);
        let expr = parse_norm(&ctx, "asinh(sinh(a + b*x))^2/x^2");
        let env = [(Symbol::new("a"), -0.3), (Symbol::new("b"), 0.8)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.4, 0.9, 1.4]);
    }

    #[test]
    fn p1_log_exp_numeric_and_symbolic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let expr = parse_norm(&ctx, "log(exp(b*x))^2");
        let env = [(Symbol::new("b"), 1.4)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.2, 0.5, 0.9]);
        let expr = parse_norm(&ctx, "log(exp(a + b*x))^3/x");
        let env = [(Symbol::new("a"), 0.6), (Symbol::new("b"), -1.1)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn p1_band_guards_certify_only_inside_the_band() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // Certified constant arguments inside the principal band cancel.
        for input in ["atan(tan(1/2))", "acosh(cosh(1/2))"] {
            let e = parse_norm(&ctx, input);
            let mut budget = MAX_CANCEL_BUDGET;
            let r = cancel_compositions(&ctx, e, x, &mut budget).expect("certified band");
            assert_eq!(r, parse_norm(&ctx, "1/2"), "{input}");
        }
        // A constant outside the band, and any symbolic argument, decline.
        for input in ["atan(tan(2))", "atan(tan(a + b*x))", "acosh(cosh(a + b*x))"] {
            let e = parse_norm(&ctx, input);
            let mut budget = MAX_CANCEL_BUDGET;
            assert!(
                cancel_compositions(&ctx, e, x, &mut budget).is_none(),
                "expected no cancellation for {input}"
            );
        }
    }

    #[test]
    fn p1_nested_pairs_reach_the_fixpoint_in_one_pass() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let e = parse_norm(&ctx, "atanh(tanh(asinh(sinh(x))))");
        let mut budget = MAX_CANCEL_BUDGET;
        let r = cancel_compositions(&ctx, e, x, &mut budget).expect("nested pair");
        assert_eq!(r, ctx.var("x"));
        // Idempotent: the rewritten form holds no pair any more.
        let mut budget = MAX_CANCEL_BUDGET;
        assert!(cancel_compositions(&ctx, r, x, &mut budget).is_none());
    }

    // ------------------------- P2: branch-shifted kernels -------------------------

    #[test]
    fn p2_acoth_tanh_corpus_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let kernel = parse_norm(&ctx, "acoth(tanh(a + b*x))");
        let proxy = parse_norm(&ctx, "a + b*x");
        let env = [(Symbol::new("a"), 0.4), (Symbol::new("b"), 1.1)];
        // rubi-01580 and rubi-01502.
        for input in ["acoth(tanh(a + b*x))^3/x^3", "x*acoth(tanh(a + b*x))^3"] {
            let expr = parse_norm(&ctx, input);
            assert_antiderivative_branch_shifted(
                &ctx,
                expr,
                kernel,
                proxy,
                x,
                &env,
                &[0.4, 0.8, 1.3],
            );
        }
        // rubi-01361 (x^(−6) cofactor).
        let expr = parse_norm(&ctx, "acoth(tanh(a + b*x))^3/x^6");
        assert_antiderivative_branch_shifted(&ctx, expr, kernel, proxy, x, &env, &[0.4, 0.8, 1.3]);
        // Numeric-coefficient variant.
        let expr = parse_norm(&ctx, "acoth(tanh(2*x))^2/x^2");
        let kernel = parse_norm(&ctx, "acoth(tanh(2*x))");
        let proxy = parse_norm(&ctx, "2*x");
        assert_antiderivative_branch_shifted(&ctx, expr, kernel, proxy, x, &[], &[0.5, 1.0, 1.6]);
    }

    #[test]
    fn p2_atan_tan_is_branch_independent() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let kernel = parse_norm(&ctx, "atan(tan(a + b*x))");
        let proxy = parse_norm(&ctx, "a + b*x");
        let env = [(Symbol::new("a"), 0.7), (Symbol::new("b"), 1.3)];
        // Samples deliberately cross the ±π/2 band edges: the kernel algebra
        // is valid on every band, unlike a naive `atan(tan(u)) → u` rewrite.
        for input in ["x^2*atan(tan(a + b*x))^2", "atan(tan(a + b*x))"] {
            let expr = parse_norm(&ctx, input);
            assert_antiderivative_branch_shifted(
                &ctx,
                expr,
                kernel,
                proxy,
                x,
                &env,
                &[1.0, 2.0, 3.0],
            );
        }
    }

    // ------------------------- P3: companion radicands -------------------------

    #[test]
    fn p3_companion_radical_with_intercept() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // asin's derivative exponent IS −1/2 and its kernel is 1 − u², so
        // ∫ asin(u)/sqrt(1 − u²) dx = asin(u)²/(2b) for u = a + b·x.
        let expr = parse_norm(&ctx, "asin(a + b*x)/sqrt(1 - (a + b*x)^2)");
        let env = [(Symbol::new("a"), 0.2), (Symbol::new("b"), 0.7)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.1, 0.4, 0.8]);
        // Symbolic-coefficient variant (zero intercept: M1 territory).
        let expr = parse_norm(&ctx, "asin(c*x)/sqrt(1 - c^2*x^2)");
        let env = [(Symbol::new("c"), 0.8)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.1, 0.5, 1.0]);
    }

    #[test]
    fn p3_rational_inverse_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        // rubi-00338: sqrt(asin(a·x))/(c − a²c·x²)^(1/2).
        let expr = parse_norm(&ctx, "sqrt(asin(a*x))/(c - a^2*c*x^2)^(1/2)");
        let env = [(Symbol::new("a"), 0.6), (Symbol::new("c"), 1.4)];
        assert_antiderivative_num(&ctx, expr, x, &env, &[0.2, 0.5, 0.9]);
        // Same family, numeric: sqrt(acosh(x))/sqrt(x² − 1) with x > 1.
        let expr = parse_norm(&ctx, "sqrt(acosh(x))/(x^2 - 1)^(1/2)");
        assert_antiderivative_num(&ctx, expr, x, &[], &[1.4, 1.9, 2.4]);
    }

    #[test]
    fn p3_declines_non_companion_radicals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // rubi-00198: atan's derivative exponent is −1, not −1/2, so the
        // radical is NOT the derivative-shaped companion of the inverse
        // function and the antiderivative is not elementary (it needs the
        // inverse tangent integral Ti₂). Declining is the correct answer.
        assert_declined(&ctx, "atan(a + b*x)/sqrt(1 + a^2 + 2*a*b*x + b^2*x^2)");
        // atanh belongs to the same (Recip) family: same mismatch.
        assert_declined(&ctx, "atanh(a + b*x)/sqrt(1 - (a + b*x)^2)");
        // Right exponent, wrong kernel (asinh needs 1 + u², not 1 − u²).
        assert_declined(&ctx, "asinh(a + b*x)/sqrt(1 - (a + b*x)^2)");
    }

    // ------------------------- new declines -------------------------

    #[test]
    fn declines_out_of_family_compositions() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Not inverse-composition pairs.
        assert_declined(&ctx, "atan(x^2)");
        assert_declined(&ctx, "atanh(sin(x))");
        assert_declined(&ctx, "sqrt(asin(x))");
        assert_declined(&ctx, "asin(x)*exp(x)");
        // acosh∘cosh is not a constant-derivative kernel (its derivative flips
        // sign at u = 0) and its cancellation band is uncertifiable.
        assert_declined(&ctx, "acosh(cosh(a + b*x))");
        // rubi-00717: atanh of an exponential-affine argument — no
        // composition pair applies, and the antiderivative is not elementary
        // (it needs polylogarithms), so declining is the only correct answer.
        assert_declined(&ctx, "x^2*atanh(a + b*f^((c + d*x)))");
        assert_declined(&ctx, "acoth(a + b*f^((c + d*x)))");
    }

    #[test]
    fn declines_hang_attributed_shapes_quickly() {
        // The 0.27.1 timeout attribution blamed `inverse_trig` for these two;
        // a traced run shows `inverse_trig` declining and the hang landing in
        // `symbolic_rational` / the heuristic chain instead. The module must
        // still reject both shapes structurally (no unbounded searching).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_declined(&ctx, "(c + d*x)^(1/2)/(a + b*x)^2");
        assert_declined(&ctx, "(e + f*x)^2*cos(c + d*x)/(a + b*sin(c + d*x))^3");
    }

    // ------------------------- stress -------------------------

    #[test]
    fn stress_repeat_is_stable_and_budget_bounded() {
        // 300 iterations over a mixed shape set: solved-ness, result strings
        // and result sizes must be byte-identical in every round.
        //
        // The calls go straight to the module entry. Every shape here resolves
        // without touching the pipeline's global chain-entry budget (a
        // per-top-level-call counter that public entry points reset):
        // `integrate_expanded_terms` integrates the expanded P1/P2 rewrites
        // with the power rule itself, so no term ever reaches the guarded
        // fallback, and the declines never re-enter at all. That makes the
        // direct loop a real statement — a module that leaked budget or grew
        // state would trip the 256-entry cap and break the stability
        // comparison long before round 300 (which is exactly what an earlier
        // version of this mechanism did, restoring the leak as a regression
        // guard).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let shapes = [
            "atanh(tanh(a + b*x))^2",
            "atanh(tanh(a + b*x))^4/x^4",
            "atanh(tanh(a + b*x))^2*sqrt(x)",
            "x*acoth(tanh(a + b*x))^3",
            "acoth(tanh(a + b*x))^3/x^3",
            "atan(x^2)",
            "atanh(sin(x))",
            "(c + d*x)^(1/2)/(a + b*x)^2",
        ];
        let solved = [true, true, true, true, true, false, false, false];
        let mut first: Vec<Option<(String, usize)>> = Vec::with_capacity(shapes.len());
        for round in 0..300 {
            for (i, input) in shapes.iter().enumerate() {
                let expr = parse_norm(&ctx, input);
                let observed =
                    integrate_inverse_trig(&ctx, expr, x).map(|r| (r.to_string(), node_count(r)));
                if round == 0 {
                    first.push(observed);
                } else {
                    assert_eq!(observed, first[i], "{input} diverged at round {round}");
                }
            }
        }
        for (i, input) in shapes.iter().enumerate() {
            if solved[i] {
                let r = first[i].as_ref().expect("expected a closed form");
                assert!(!r.0.contains("Integral"), "{input} left a residue: {}", r.0);
            } else {
                assert!(first[i].is_none(), "{input}: {:?}", first[i]);
            }
        }
    }
}
