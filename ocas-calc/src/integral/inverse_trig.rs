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
///    radicands (P3), branch-shifted composition kernels (P2), and finally
///    the inv-hyp substitution (M2).
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
