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
//! - **M4 inverse-function substitution** (0.30.0 Wave 5,
//!   `docs/planning/INV_SUBST_SPEC_CN.md`): `(a + b·F(w))^k` times a
//!   polynomial weight, `w = ρ + σ·x` affine, integer `1 ≤ k ≤ 6`, all six
//!   main heads, plus (W5b) at most one companion radical
//!   `(s·E(w))^(m/2)`, integer `|m| ≤ 5`, absorbed as
//!   `s^(m/2)·G′(u)^m·σ` with the sheet factor σ = +1 on the principal
//!   band. `u = F(w)` maps the integrand to
//!   `P(u)·T(u)`; trig products reduce to multiple angles and re-enter the
//!   chain, hyperbolic families use the M2-style tabular closed form
//!   directly, and atanh uses a dedicated parts/tanh-reduction tabular.
//!   The back-substitution is a whitelist rewrite (`u → F(w)`,
//!   `G(u) → w`, `G′(u) → +√E(w)` on the principal band, multiple angles
//!   via Chebyshev); any residue, leftover `u`, or `F(G(·))` composite
//!   declines the whole shape. See the M4 section docs for the real-branch
//!   semantics and the deviations from the spec.
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
fn family_kernel_of_u<'a>(ctx: &'a AtomArena<'a>, f: InvFun, u: Atom<'a>) -> Atom<'a> {
    let u2 = int_pow(ctx, u, 2);
    match f {
        InvFun::Asin | InvFun::Acos | InvFun::Atanh => {
            normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), u2])]))
        }
        InvFun::Atan | InvFun::Asinh => normalize(ctx, ctx.add(&[ctx.num(1), u2])),
        InvFun::Acosh => normalize(ctx, ctx.add(&[ctx.num(-1), u2])),
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
    if ke != family_kernel_exp(f) {
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
    if family_sign(f) < 0 {
        out.push(ctx.num(-1));
    }
    out.push(s_factor);
    out.push(ctx.pow(base, rat_atom(ctx, p + q, q)));
    let denom = normalize(ctx, ctx.mul(&[b, sigma, rat_atom(ctx, p + q, q)]));
    out.push(inv(ctx, denom));
    Some(normalize(ctx, ctx.mul(&out)))
}

// =========================================================================
// M4: inverse-function substitution (0.30.0 Wave 5a)
// =========================================================================
//
// Implementation of `docs/planning/INV_SUBST_SPEC_CN.md` §2.1/2.3/2.4/2.6:
// `∫ C·(a + b·F(w))^k · P(x) dx` with `F` one of the six main inverse
// heads, `w = ρ + σ·x` affine (`σ ≢ 0`), integer `1 ≤ k ≤ MAX_IS_K` and
// `P` polynomial of total degree ≤ `MAX_IS_POLY_DEG`. The substitution
// `u = F(w)` maps the integrand to `P(u)·T(u)` with `T` a trig/hyperbolic
// monomial (`x = (G(u) − ρ)/σ`, `dx = G′(u)/σ·du`, per the spec §2.1
// table); elementarity is **not** prejudged — a re-integration residue
// declines the whole shape (§2.5).
//
// W5a/W5b scope limits (deviations from the full spec are all narrowings or
// measured reorderings):
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
//   - ratio families (atan/atanh): `ratio_i`/`ratio_j`, parts on the
//     exact differential `dv = R^a·S² du` plus the `R^c → R^{c−2}` peel
//     (`tan² = sec² − 1` / `tanh² = 1 − sech²`); the `J(j≥1, 1)` base
//     case is the polylog territory of §2.5 and declines.
// - **tan/sec even-power pre-reduction** (§2.3 second clause) is
//   subsumed by the ratio tabular, which handles `tan^A·sec²` directly.
// - Reciprocal heads (acsc/asec/acoth/asech/acsch) are W5c, not handled.
//
// Real-branch semantics (spec §4): every answer is exact on the principal
// band of `F` (asin: `u ∈ [−π/2, π/2]`; acos: `[0, π]`; atan:
// `(−π/2, π/2)`; asinh/atanh: ℝ; acosh: `[0, ∞)`). On the band `G′(u)`
// has a definite non-negative sign for all six families — the
// back-substitution (`cos u → +√E(w)`) already assumes it. The radical
// absorption extends the same convention: for odd `m`,
// `(s·E)^(m/2) = s^(m/2)·G′^m·σ` folds the sheet factor σ to `+1` with
// the principal-band sign evidence (§4.2), and the `√E(w)` answers
// re-enter the certificate's radical layer through the `r² = E(w)`
// relation.
// The in-module tabulars never produce `F(G(u))` composites; the terminal
// structural check (`contains_inv_composite`) is the belt-and-braces
// enforcement of §2.4's ban. No `I` is introduced anywhere.

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

/// `G` head of the family (`x = G(u)` at zero intercept), per spec §2.1.
fn family_g_head(f: InvFun) -> &'static str {
    match f {
        InvFun::Asin => "sin",
        InvFun::Acos => "cos",
        InvFun::Atan => "tan",
        InvFun::Asinh => "sinh",
        InvFun::Acosh => "cosh",
        InvFun::Atanh => "tanh",
    }
}

/// `G′(u)` per family as `(head, sign, power)`: asin `cos`, acos `−sin`,
/// atan `sec²`, asinh `cosh`, acosh `sinh`, atanh `sech²`.
fn family_gp_head(f: InvFun) -> (&'static str, i64, i64) {
    match f {
        InvFun::Asin => ("cos", 1, 1),
        InvFun::Acos => ("sin", -1, 1),
        InvFun::Atan => ("sec", 1, 2),
        InvFun::Asinh => ("cosh", 1, 1),
        InvFun::Acosh => ("sinh", 1, 1),
        InvFun::Atanh => ("sech", 1, 2),
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
/// `+1` on the band), and the remaining heads of the class by the
/// algebraic identities. `None` for heads of the other class.
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
            _ => return None,
        },
        InvFun::Acos => match name {
            "cos" => w,
            "sin" => r,
            "tan" => r_over_w(ctx),
            "sec" => inv(ctx, w),
            _ => return None,
        },
        InvFun::Atan => match name {
            "tan" => w,
            "sec" => r,
            "sin" => w_over_r(ctx),
            "cos" => inv(ctx, r),
            _ => return None,
        },
        InvFun::Asinh => match name {
            "sinh" => w,
            "cosh" => r,
            "tanh" => w_over_r(ctx),
            "sech" => inv(ctx, r),
            _ => return None,
        },
        InvFun::Acosh => match name {
            "cosh" => w,
            "sinh" => r,
            "tanh" => r_over_w(ctx),
            "sech" => inv(ctx, w),
            _ => return None,
        },
        InvFun::Atanh => match name {
            "tanh" => w,
            "sech" => r,
            "sinh" => w_over_r(ctx),
            "cosh" => inv(ctx, r),
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
/// the same six pairs [`contains_inv_composite`] enforces on the answer.
fn is_inv_pair(f: &str, g: &str) -> bool {
    matches!(
        (f, g),
        ("asin", "sin")
            | ("acos", "cos")
            | ("atan", "tan")
            | ("asinh", "sinh")
            | ("acosh", "cosh")
            | ("atanh", "tanh")
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
                    "sin" | "cos" | "tan" | "sec" | "sinh" | "cosh" | "tanh" | "sech"
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
/// an `F(G(·))` inverse-of-forward composite for any of the six pairs.
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
                && matches!(
                    (name.as_str(), inner.as_str()),
                    ("asin", "sin")
                        | ("acos", "cos")
                        | ("atan", "tan")
                        | ("asinh", "sinh")
                        | ("acosh", "cosh")
                        | ("atanh", "tanh")
                )
            {
                return true;
            }
            args.iter().any(|a| contains_inv_composite(*a))
        }
    }
}

/// The ratio-family heads: `(R, S)` with `R′ = ±S²` — `(tan, sec)` for
/// atan, `(tanh, sech)` for atanh.
fn ratio_heads(hyper: bool) -> (&'static str, &'static str) {
    if hyper {
        ("tanh", "sech")
    } else {
        ("tan", "sec")
    }
}

/// Decompose an expanded ratio-family t-term into `(constants, j, A, B)`
/// for `C·u^j·R(u)^A·S(u)^B`. Without a companion radical the `dx` factor
/// contributes exactly `S²`; W5b's absorption `(s·E)^(m/2) → s^(m/2)·S^m·σ`
/// shifts it to `S^(2+m)`, so `B` is any integer here (negative `B` is an
/// exact reciprocal power of `cos`/`cosh` — see [`integrate_ratio_term`]).
fn scan_ratio_term<'a>(
    term: Atom<'a>,
    t_sym: Symbol,
    hyper: bool,
) -> Option<(Vec<Atom<'a>>, i64, i64, i64)> {
    let (tan, sec) = ratio_heads(hyper);
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
/// `dv = R^a·S² du` (`v = R(u)`, since `R′ = ±S²`), reducing to
/// [`ratio_j`]. `None` when the recursion reaches a polylog term (§2.5).
fn ratio_i<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    a: i64,
    t: Atom<'a>,
    hyper: bool,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(0..=2 * MAX_IS_MU).contains(&a) {
        return None;
    }
    let (tan, _) = ratio_heads(hyper);
    let r_t = ctx.fun(tan, &[t]);
    let boundary = normalize(
        ctx,
        ctx.mul(&[
            int_pow(ctx, t, j),
            int_pow(ctx, r_t, a + 1),
            inv(ctx, ctx.num(a + 1)),
        ]),
    );
    if j == 0 {
        return Some(boundary);
    }
    let inner = ratio_j(ctx, j - 1, a + 1, t, hyper)?;
    let coeff = rat_atom(ctx, -j, a + 1);
    Some(normalize(
        ctx,
        ctx.add(&[boundary, ctx.mul(&[coeff, inner])]),
    ))
}

/// `J(j, c) = ∫ u^j·R(u)^c du`. tan: `tan^c = tan^{c−2}·(sec²−1)` gives
/// `J = I(c−2) − J(c−2)`; tanh: `tanh^c = tanh^{c−2}·(1−sech²)` gives
/// `J = J(c−2) − I(c−2)`. Bases: `J(j, 0) = u^{j+1}/(j+1)` and
/// `J(0, 1) = −log(cos u)` (tan) resp. `log(cosh u)` (tanh);
/// `J(j ≥ 1, 1)` is polylog territory and declines (§2.5 — elementarity
/// emerges from the reduction base case, the residue-check analogue).
fn ratio_j<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    c: i64,
    t: Atom<'a>,
    hyper: bool,
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
            let co = ctx.fun(if hyper { "cosh" } else { "cos" }, &[t]);
            let log = ctx.fun("log", &[co]);
            Some(if hyper {
                log
            } else {
                normalize(ctx, ctx.mul(&[ctx.num(-1), log]))
            })
        }
        _ => {
            let jm = ratio_j(ctx, j, c - 2, t, hyper)?;
            let im = ratio_i(ctx, j, c - 2, t, hyper)?;
            let (pos, neg) = if hyper { (jm, im) } else { (im, jm) };
            Some(normalize(
                ctx,
                ctx.add(&[pos, ctx.mul(&[ctx.num(-1), neg])]),
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
    hyper: bool,
) -> Option<Atom<'a>> {
    let (coeff, j, a, b) = scan_ratio_term(term, t_sym, hyper)?;
    let t = ctx.var(t_sym.as_str());
    if b < 0 {
        // `sec = 1/cos` / `sech = 1/cosh` hold exactly (no branch), so
        // `R^A·S^B = s(u)^A·c(u)^(−A−B)`: a pure sin/cos (sinh/cosh) term,
        // integrable by the trig/hyperbolic tabulars when `−A−B ≥ 0`.
        // Otherwise the reciprocal cosine power survives — decline (§2.5).
        let c_pow = -a - b;
        if c_pow < 0 {
            return None;
        }
        let (s_name, c_name) = if hyper {
            ("sinh", "cosh")
        } else {
            ("sin", "cos")
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
        return integrate_sincos_term(ctx, conv, t_sym, hyper);
    }
    let core = ratio_rs(ctx, j, a, b, t, hyper)?;
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

/// `∫ u^j·R(u)^A·S(u)^B du` with `B ≥ 0` (the `R′ = S²` family: `(R, S)` is
/// `(tan, sec)` resp. `(tanh, sech)`; `S′ = ±R·S` with `−` hyperbolic).
///
/// - even `B ≥ 2`: `S^B = S²·(1 ± R²)^((B−2)/2)` expands to `R^a·S²` terms
///   for the [`ratio_i`] tabular;
/// - odd `B`, odd `A ≥ 1`: parts on `d(S^B) = ±B·R·S^B du` with
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
    hyper: bool,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(0..=2 * MAX_IS_MU).contains(&a) {
        return None;
    }
    if b == 0 {
        return ratio_j(ctx, j, a, t, hyper);
    }
    if !(1..=2 * MAX_IS_MU).contains(&b) {
        return None;
    }
    let (_, sec) = ratio_heads(hyper);
    let s_t = ctx.fun(sec, &[t]);
    if b % 2 == 0 {
        // S² = 1 ∓ R²: hyper `(1 − R²)^k`, trig `(1 + R²)^k`, k = (B−2)/2.
        let k = (b - 2) / 2;
        let mut terms = Vec::new();
        for i in 0..=k {
            let c = binom(k, i) * if hyper && i % 2 == 1 { -1 } else { 1 };
            let core = ratio_i(ctx, j, a + 2 * i, t, hyper)?;
            terms.push(normalize(ctx, ctx.mul(&[ctx.num(c), core])));
        }
        return Some(normalize(ctx, ctx.add(&terms)));
    }
    if a == 0 {
        return ratio_s_pow(ctx, j, b, t, hyper);
    }
    // Sign of `S′ = sgn·R·S` and of `(R²)^k` expanded in `S²`:
    // hyperbolic `R² = 1 − S²` gives `(−1)^i`; trig `R² = S² − 1` gives
    // `(−1)^(k−i)`.
    let sq_sign = |k: i64, i: i64| {
        if if hyper { i % 2 == 1 } else { (k - i) % 2 == 1 } {
            -1i64
        } else {
            1
        }
    };
    if a % 2 == 1 {
        // dv = R^A·S^B du = R^(A−1)·d(S^B)/(sgn·B): with
        // R^(A−1) = Σ_i c_i·S^(2i) the primitive is
        // v = Σ_i c_i·S^(B+2i)/(sgn·(B+2i)).
        let k = (a - 1) / 2;
        let sgn: i64 = if hyper { -1 } else { 1 };
        let mut v_terms = Vec::new();
        for i in 0..=k {
            let c = binom(k, i) * sq_sign(k, i);
            let p = b + 2 * i;
            v_terms.push((c, p));
        }
        let mut v_atoms = Vec::with_capacity(v_terms.len());
        let mut inner = Vec::with_capacity(v_terms.len());
        for (c, p) in v_terms {
            let coeff = rat_atom(ctx, c, sgn * p);
            v_atoms.push(normalize(ctx, ctx.mul(&[coeff, int_pow(ctx, s_t, p)])));
            if j > 0 {
                let red = ratio_s_pow(ctx, j - 1, p, t, hyper)?;
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
    for i in 0..=k {
        let c = binom(k, i) * sq_sign(k, i);
        let red = ratio_s_pow(ctx, j, b + 2 * i, t, hyper)?;
        terms.push(normalize(ctx, ctx.mul(&[ctx.num(c), red])));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// `∫ u^j·S(u)^b du`, odd `b ≥ 1`. The `j = 0` reduction
/// `∫S^b = (R·S^(b−2) + (b−2)·∫S^(b−2))/(b−1)` is family-independent
/// (`d(R·S^(b−2)) = (b−1)·S^b − (b−2)·S^(b−2)` in both families); the bases
/// are `∫sec u du = log(sec u + tan u)` and `∫sech u du = asin(tanh u)`
/// (both exact on the principal band; the latter back-substitutes to
/// `asin(w)` through the whitelist's cross-family head map). `j ≥ 1`
/// integrates by parts on `v = ∫S^b du`; the surviving `u^j·S(u)` (`j ≥ 1`)
/// base is polylog territory and declines (§2.5).
fn ratio_s_pow<'a>(
    ctx: &'a AtomArena<'a>,
    j: i64,
    b: i64,
    t: Atom<'a>,
    hyper: bool,
) -> Option<Atom<'a>> {
    if !(0..=MAX_IS_K).contains(&j) || !(1..=2 * MAX_IS_MU).contains(&b) {
        return None;
    }
    let (tan, sec) = ratio_heads(hyper);
    let r_t = ctx.fun(tan, &[t]);
    let s_t = ctx.fun(sec, &[t]);
    if b == 1 {
        if j >= 1 {
            return None;
        }
        return Some(if hyper {
            // d/du asin(tanh u) = sech²u/√(1 − tanh²u) = sech u (sech > 0).
            ctx.fun("asin", &[r_t])
        } else {
            ctx.fun("log", &[normalize(ctx, ctx.add(&[s_t, r_t]))])
        });
    }
    if j == 0 {
        let rec = ratio_s_pow(ctx, 0, b - 2, t, hyper)?;
        return Some(normalize(
            ctx,
            ctx.mul(&[
                inv(ctx, ctx.num(b - 1)),
                ctx.add(&[
                    ctx.mul(&[r_t, int_pow(ctx, s_t, b - 2)]),
                    ctx.mul(&[ctx.num(b - 2), rec]),
                ]),
            ]),
        ));
    }
    // Parts with `v = ∫S^b du = (R·S^(b−2) + (b−2)·H_{b−2})/(b−1)`:
    // `∫ u^j·S^b = u^j·H_b − j·∫ u^(j−1)·H_b du`.
    let h_b = ratio_s_pow(ctx, 0, b, t, hyper)?;
    let boundary = normalize(ctx, ctx.mul(&[int_pow(ctx, t, j), h_b]));
    let t1 = ratio_rs(ctx, j - 1, 1, b - 2, t, hyper)?; // ∫ u^(j−1)·R·S^(b−2)
    let t2 = ratio_s_pow(ctx, j - 1, b - 2, t, hyper)?; // ∫ u^(j−1)·S^(b−2)
    let inner = normalize(
        ctx,
        ctx.mul(&[
            inv(ctx, ctx.num(b - 1)),
            ctx.add(&[t1, ctx.mul(&[ctx.num(b - 2), t2])]),
        ]),
    );
    Some(normalize(
        ctx,
        ctx.add(&[boundary, ctx.mul(&[ctx.num(-j), inner])]),
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
    // dx = gp_sign·G′(u)^gp_pow/σ·du. A companion radical contributes
    // `s^(rm/2)·G′(u)^rm·σ` with the *same* unsigned head (E(G(u)) is
    // `G′²` for the trig/hyperbolic families and `G′` itself for the ratio
    // families, so `E^(rm/2) = head^rm·σ` in both), and the two merge into
    // a single head power (§2.2). The sheet factor σ (odd `rm`) is `+1`
    // on the principal band of every family — see the section docs for the
    // §4.2 emission record.
    let (gp_name, gp_sign, gp_pow) = family_gp_head(m.f);
    let gp_t = ctx.fun(gp_name, &[t]);
    if gp_sign < 0 {
        t_factors.push(ctx.num(-1));
    }
    let mut gp_total = gp_pow;
    if let Some((rm, s_atom)) = &rad_data {
        gp_total += rm;
        t_factors.push(ctx.pow(*s_atom, rat_atom(ctx, *rm, 2)));
        if std::env::var_os("OCAS_M4_DEBUG").is_some() {
            eprintln!(
                "[m4] companion radical (s·E)^({rm}/2) absorbed; sheet factor σ folded to +1 (principal band)"
            );
        }
    }
    t_factors.push(normalize(
        ctx,
        ctx.mul(&[int_pow(ctx, gp_t, gp_total), inv(ctx, sigma)]),
    ));
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
        // Ratio families (atan/atanh): the t-terms are
        // `C·u^j·R(u)^A·S(u)²`, integrated by the parts/R-reduction
        // tabular (a polylog base case declines).
        InvFun::Atan => {
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                out.push(integrate_ratio_term(ctx, term, t_sym, false)?);
            }
            normalize(ctx, ctx.add(&out))
        }
        InvFun::Atanh => {
            let mut out = Vec::with_capacity(terms.len());
            for term in terms {
                if is_work_exceeded(&mut work, node_count(term)) {
                    return None;
                }
                out.push(integrate_ratio_term(ctx, term, t_sym, true)?);
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
