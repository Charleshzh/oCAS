//! exp/log-kernel substitution mechanisms (0.27.1 Phase 1C).
//!
//! Three mechanisms for integrands built from a single exponential or
//! logarithm kernel:
//!
//! - **E1 exp-kernel rationalization**: integrands rational in `e^(a·x)`
//!   are mapped by `t = e^(a·x)` to a rational function of `t` (numeric
//!   coefficients are reduced to a rational-gcd base, so `exp(x)` and
//!   `exp(2·x)` share `t = e^x`). Hyperbolic-rational forms
//!   `R(sinh u, cosh u, …)` with `u = a·x + b` linear (symbolic `a`
//!   allowed) take `t = e^u` and rewrite `sinh u = (t²−1)/(2t)` etc.,
//!   which is again rational in `t`. The whole integrand must reduce to a
//!   rational function of the kernel: a bare `x` factor (rule B2's
//!   `x^n·e^(a·x)`) or a `sin`/`cos` factor (B3/B4 — routing those here
//!   once caused parts loops, the 0.27 lesson) declines structurally.
//! - **E2 log-kernel substitution**: `f(log x)/x` maps by `t = log x` to
//!   `f(t)`. Every occurrence of `x` must be inside `log(x)` kernels plus
//!   exactly one `x^−1` factor. Extension: `log(c·(d+e·x)^n)` and
//!   `log(x^m)` with small integer exponents expand formally to
//!   `n·log(d+e·x) + log c` / `m·log x` and reintegrate the expansion.
//! - **E3** fills the `m = −1` gap of rule B6: `log(x)^k/x` has the closed
//!   form `log(x)^(k+1)/(k+1)` for integer `k ≠ −1` (`k = −1` is B7's
//!   `1/(x·log x)`).
//! - **E4 exp-of-inverse algebraization**: `exp(k·F(u))` with `F` an inverse
//!   trigonometric/hyperbolic function is *algebraic* in `u`, not
//!   transcendental, so the site is rewritten by the exact identities
//!   `exp(n·atanh u) = ((1+u)/(1−u))^(n/2)`,
//!   `exp(n·acoth u) = ((1+u)/(u−1))^(n/2)`,
//!   `exp(i·r·atan u) = ((1+i·u)/(1−i·u))^(r/2)`,
//!   `exp(k·asinh u) = (u+√(u²+1))^k` and
//!   `exp(k·acosh u) = (u+√(u²−1))^k`, and the resulting
//!   algebraic/rational integrand is handed back to the chain. Only a linear
//!   argument `u = σ·x + τ` is in scope, and only exponents that keep the
//!   rewrite exact are accepted (`k` rational for `asinh`/`acosh`, `k = i·r`
//!   with `r` rational for `atan`); anything else declines.
//!
//! All chain re-entries go through `integrate_raw` and are declined on any
//! `Integral` residue. The substituted t-forms are rational in `t` (E1) or
//! kernel-free (E2), so the module cannot re-match its own output. E4's
//! rewritten form contains no `exp` of an inverse function either, so it
//! cannot re-match E4.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    contains_integral, int_pow, integrate_raw, inv, is_constant, lcm_i64, linear_form, node_count,
    pick_subst_symbol, rat_atom, replace_symbol,
};

/// Cap on the substituted t-form's numerator/denominator degree in `t`, on
/// the integer exponent of any single kernel power, and on log powers (E3).
const MAX_KERNEL_DEG: i64 = 8;
/// Degree cap for the hyperbolic-rational t-form of E1, whose integrand is a
/// product of kernel powers over a common denominator before substitution
/// (see `try_exp_kernel`). Only the degree check uses it; the kernel-exponent
/// cap stays `MAX_KERNEL_DEG`.
const MAX_HYP_DEG: i64 = 12;
/// Node budget for the input integrand and for the substituted t-form.
const MAX_SUBST_NODES: usize = 200;
/// Cap on the number of exp/hyperbolic/log kernel sites in one integrand.
const MAX_KERNELS: usize = 16;
/// `|n|` cap for the `log(c·g^n)` formal power expansion.
const MAX_LOG_EXPAND: i64 = 4;
/// Numerator cap `|p| ≤ 8` for a rational E4 kernel exponent `k = p/q`.
const MAX_INV_EXP_NUM: i64 = 8;
/// Denominator cap `q ≤ 4` for a rational E4 kernel exponent `k = p/q`.
const MAX_INV_EXP_DEN: i64 = 4;
/// The imaginary-unit spelling used by the corpus (`exp(4*i*atan(a*x))`).
const IMAG_UNIT: &str = "i";

/// Integrate `expr` via the exp/log-kernel mechanisms (see module docs).
/// Returns `None` when no mechanism applies or a budget is exceeded.
pub(crate) fn integrate_exp_log<'a>(
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
    if node_count(expr) > MAX_SUBST_NODES {
        return None;
    }
    try_log_power_over_x(ctx, expr, var)
        .or_else(|| try_log_kernel_subst(ctx, expr, var))
        .or_else(|| try_log_power_expand(ctx, expr, var))
        .or_else(|| try_exp_kernel(ctx, expr, var))
        .or_else(|| try_exp_inverse(ctx, expr, var))
}

// =========================================================================
// E3: log(x)^k / x closed form
// =========================================================================

/// ∫ `C·log(x)^k/x` dx = `C·log(x)^(k+1)/(k+1)` for integer `k ∉ {−1, 0}`,
/// `|k| ≤ 8`. Covers both normalized shapes of the corpus forms
/// `log(x)^n/x` and `1/(x·log(x)^n)`.
fn try_log_power_over_x<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (consts, k) = match_log_over_x(ctx, expr, var)?;
    if k == -1 || k == 0 || k.abs() > MAX_KERNEL_DEG {
        return None;
    }
    let log_x = ctx.fun("log", &[ctx.var(var.as_str())]);
    let k1 = k.checked_add(1)?;
    let mut factors = consts;
    factors.push(int_pow(ctx, log_x, k1));
    factors.push(inv(ctx, ctx.num(k1)));
    Some(normalize(ctx, ctx.mul(&factors)))
}

/// Decompose `expr` as `C · x^−1 · log(x)^k`, returning the constant
/// factors and `k`. Handles flat products carrying an `x^−1` factor and
/// `(C·x·log(x)^n)^−1` wrapping a product with a bare `x` factor.
fn match_log_over_x<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Vec<Atom<'a>>, i64)> {
    let (factors, sign): (&[Atom<'a>], i64) = match expr.node() {
        AtomNode::Mul(args) => (args, 1),
        AtomNode::Pow(b, e) if matches!(e.node(), AtomNode::Num(-1)) => match b.node() {
            AtomNode::Mul(args) => (args, -1),
            _ => return None,
        },
        _ => return None,
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut x_count = 0u32;
    let mut log_pow: Option<i64> = None;
    for f in factors {
        classify_log_factor(ctx, *f, var, sign, &mut consts, &mut x_count, &mut log_pow)?;
    }
    if x_count != 1 {
        return None;
    }
    Some((consts, log_pow?))
}

/// Classify one factor of the `log(x)^k/x` product. `sign` is `1` for a
/// flat product (the `x^−1` site is a `Pow(x, −1)` factor) and `−1` for a
/// product wrapped in `(…)^−1` (the site is a bare `x` factor, and log
/// powers contribute negated exponents).
#[allow(clippy::too_many_arguments)]
fn classify_log_factor<'a>(
    ctx: &'a AtomArena<'a>,
    f: Atom<'a>,
    var: Symbol,
    sign: i64,
    consts: &mut Vec<Atom<'a>>,
    x_count: &mut u32,
    log_pow: &mut Option<i64>,
) -> Option<()> {
    // Bare x in a flat product is `x^1·log^k`, rule B6 territory: decline.
    if matches!(f.node(), AtomNode::Var(v) if *v == var) {
        if sign == 1 {
            return None;
        }
        *x_count += 1;
        return Some(());
    }
    // The x^−1 site of a flat product.
    if sign == 1
        && let AtomNode::Pow(b, e) = f.node()
        && matches!(b.node(), AtomNode::Var(v) if *v == var)
        && matches!(e.node(), AtomNode::Num(-1))
    {
        *x_count += 1;
        if *x_count > 1 {
            return None;
        }
        return Some(());
    }
    if let Some(k) = log_kernel_power(f, var) {
        if log_pow.is_some() {
            return None;
        }
        *log_pow = Some(k.checked_mul(sign)?);
        return Some(());
    }
    if is_constant(f, var) {
        consts.push(if sign == 1 { f } else { inv(ctx, f) });
        return Some(());
    }
    None
}

/// The integer exponent of a `log(x)` power: `log(x)` → 1,
/// `log(x)^m` → m, `(log(x)^n)^m` → n·m.
fn log_kernel_power(f: Atom<'_>, var: Symbol) -> Option<i64> {
    match f.node() {
        AtomNode::Fun(name, args)
            if name.as_str() == "log"
                && args.len() == 1
                && matches!(args[0].node(), AtomNode::Var(v) if *v == var) =>
        {
            Some(1)
        }
        AtomNode::Pow(b, e) => {
            let inner = log_kernel_power(*b, var);
            match (inner, e.node()) {
                (Some(k), AtomNode::Num(m)) => k.checked_mul(*m),
                _ => None,
            }
        }
        _ => None,
    }
}

// =========================================================================
// E2: f(log x)/x substitution
// =========================================================================

/// ∫ `f(log x)/x` dx with `t = log x`: every occurrence of `x` must be
/// inside `log(x)` kernels plus exactly one `x^−1` factor.
fn try_log_kernel_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let mut x_inv = 0u32;
    let mut logs = 0u32;
    let f_t = log_subst(ctx, expr, var, t, &mut x_inv, &mut logs)?;
    if x_inv != 1 || logs == 0 {
        return None;
    }
    let integrand_t = normalize(ctx, f_t);
    if node_count(integrand_t) > MAX_SUBST_NODES {
        return None;
    }
    let result_t = integrate_raw(ctx, integrand_t, t_sym, 0, true, 0, 0);
    if contains_integral(result_t) {
        return None;
    }
    let log_x = ctx.fun("log", &[ctx.var(var.as_str())]);
    Some(normalize(ctx, replace_symbol(ctx, result_t, t_sym, log_x)))
}

/// Substitute `log(x) → t` throughout `expr` and erase the single `x^−1`
/// factor (`dt = dx/x`). Returns `None` on any other occurrence of `var`.
fn log_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    t: Atom<'a>,
    x_inv: &mut u32,
    logs: &mut u32,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) {
        return Some(expr);
    }
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        // A bare occurrence of the variable is only allowed as the `x^−1`
        // factor handled in the `Pow` branch below.
        AtomNode::Var(_) => None,
        AtomNode::Pow(b, e) => {
            // The `x^−1` site: flat `x^−1`, or `(x·g)^−1` carrying exactly
            // one bare `x` factor (the reciprocal-of-product parse shape).
            if matches!(e.node(), AtomNode::Num(-1)) {
                if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                    *x_inv += 1;
                    if *x_inv > 1 {
                        return None;
                    }
                    return Some(ctx.num(1));
                }
                if let AtomNode::Mul(margs) = b.node() {
                    let mut rebuilt: Vec<Atom<'a>> = Vec::with_capacity(margs.len());
                    let mut found_x = false;
                    for a in margs.iter() {
                        if matches!(a.node(), AtomNode::Var(v) if *v == var) {
                            if found_x {
                                return None;
                            }
                            found_x = true;
                            continue;
                        }
                        rebuilt.push(log_subst(ctx, *a, var, t, x_inv, logs)?);
                    }
                    if found_x {
                        if rebuilt.is_empty() {
                            return None;
                        }
                        *x_inv += 1;
                        if *x_inv > 1 {
                            return None;
                        }
                        let inner = if rebuilt.len() == 1 {
                            rebuilt[0]
                        } else {
                            ctx.mul(&rebuilt)
                        };
                        return Some(ctx.pow(inner, ctx.num(-1)));
                    }
                }
            }
            // log(x)^k with integer k.
            if let AtomNode::Fun(name, fargs) = b.node()
                && name.as_str() == "log"
                && fargs.len() == 1
                && matches!(fargs[0].node(), AtomNode::Var(v) if *v == var)
            {
                let AtomNode::Num(k) = e.node() else {
                    return None;
                };
                if k.abs() > MAX_KERNEL_DEG || *logs as usize >= MAX_KERNELS {
                    return None;
                }
                *logs += 1;
                return Some(int_pow(ctx, t, *k));
            }
            if !is_constant(*e, var) {
                return None;
            }
            let nb = log_subst(ctx, *b, var, t, x_inv, logs)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Fun(name, args) => {
            if name.as_str() == "log"
                && args.len() == 1
                && matches!(args[0].node(), AtomNode::Var(v) if *v == var)
            {
                if *logs as usize >= MAX_KERNELS {
                    return None;
                }
                *logs += 1;
                return Some(t);
            }
            // Other functions of the variable (sin(x), exp(x), log(2·x))
            // fail here through their arguments.
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(log_subst(ctx, *a, var, t, x_inv, logs)?);
            }
            Some(ctx.fun(name.as_str(), &rebuilt))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(log_subst(ctx, *a, var, t, x_inv, logs)?);
            }
            Some(if matches!(expr.node(), AtomNode::Add(_)) {
                ctx.add(&rebuilt)
            } else {
                ctx.mul(&rebuilt)
            })
        }
    }
}

// =========================================================================
// E2 extension: log of a small integer power
// =========================================================================

/// `C · log(c·g^n)` with `g = d + e·x` linear in `var` and `n` a small
/// integer uses the formal expansion `log(c·g^n) = n·log g + log c` and
/// the standard linear-log antiderivative `∫ log(g) dx = (g·log g − g)/e`.
/// (Verified: the chain does not solve a bare `log(d+e·x)` — the parts
/// heuristic requires a genuine product — so the closed form is built
/// here.) `log(x^m)` is the `g = x`, `c = 1` case.
fn try_log_power_expand<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut log_arg: Option<Atom<'a>> = None;
    for f in factors {
        if is_constant(f, var) {
            consts.push(f);
            continue;
        }
        if log_arg.is_none()
            && let AtomNode::Fun(name, args) = f.node()
            && name.as_str() == "log"
            && args.len() == 1
        {
            log_arg = Some(args[0]);
            continue;
        }
        return None;
    }
    let (base, n, inner) = match_log_pow_arg(log_arg?, var)?;
    let (e, _d) = linear_form(ctx, base, var)?;
    if matches!(e.node(), AtomNode::Num(0)) {
        return None;
    }
    // ∫ n·log(g) dx = n·(g·log g − g)/e.
    let log_g = ctx.fun("log", &[base]);
    let g_term = ctx.mul(&[
        ctx.num(n),
        inv(ctx, e),
        ctx.add(&[ctx.mul(&[base, log_g]), ctx.mul(&[ctx.num(-1), base])]),
    ]);
    let mut body = vec![g_term];
    if !inner.is_empty() {
        // ∫ log(c) dx = log(c)·x.
        let c = normalize(ctx, ctx.mul(&inner));
        body.push(ctx.mul(&[ctx.fun("log", &[c]), ctx.var(var.as_str())]));
    }
    let mut all = consts;
    all.push(if body.len() == 1 {
        body[0]
    } else {
        ctx.add(&body)
    });
    Some(normalize(ctx, ctx.mul(&all)))
}

/// Match the argument of `log(arg)` as `g^n` or `c·g^n` with
/// `2 ≤ |n| ≤ MAX_LOG_EXPAND`; returns `(g, n, c's factors)`.
fn match_log_pow_arg<'a>(arg: Atom<'a>, var: Symbol) -> Option<(Atom<'a>, i64, Vec<Atom<'a>>)> {
    match arg.node() {
        AtomNode::Pow(b, e) => {
            let AtomNode::Num(n) = e.node() else {
                return None;
            };
            if !(2..=MAX_LOG_EXPAND).contains(&n.abs()) {
                return None;
            }
            Some((*b, *n, Vec::new()))
        }
        AtomNode::Mul(args) => {
            let mut consts: Vec<Atom<'a>> = Vec::new();
            let mut pow: Option<(Atom<'a>, i64)> = None;
            for f in args.iter() {
                if is_constant(*f, var) {
                    consts.push(*f);
                    continue;
                }
                if pow.is_none()
                    && let AtomNode::Pow(b, e) = f.node()
                    && let AtomNode::Num(n) = e.node()
                    && (2..=MAX_LOG_EXPAND).contains(&n.abs())
                {
                    pow = Some((*b, *n));
                    continue;
                }
                return None;
            }
            if consts.is_empty() {
                return None;
            }
            let (b, n) = pow?;
            Some((b, n, consts))
        }
        _ => None,
    }
}

// =========================================================================
// E1: exp-kernel rationalization (incl. hyperbolic-rational forms)
// =========================================================================

/// The substitution base: `t = exp(a·x + b)` with `a` nonzero.
struct ExpBase<'a> {
    a: Atom<'a>,
    b: Atom<'a>,
    /// The exponent expression `a·x + b`, for back-substitution.
    arg: Atom<'a>,
}

fn try_exp_kernel<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let mut exps: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    let mut hyps: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    collect_kernels(ctx, expr, var, &mut exps, &mut hyps)?;
    if exps.is_empty() && hyps.is_empty() {
        return None;
    }
    let base = decide_base(ctx, &exps, &hyps, var)?;
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let substituted = exp_subst(ctx, expr, var, &base, t)?;
    // t = e^(a·x+b) → dt = a·t dx → dx = dt/(a·t).
    let integrand_t = normalize(
        ctx,
        ctx.mul(&[substituted, inv(ctx, base.a), ctx.pow(t, ctx.num(-1))]),
    );
    if node_count(integrand_t) > MAX_SUBST_NODES {
        return None;
    }
    // Rationality check and degree budget in one structural pass.
    let (num_deg, den_deg) = t_degree_bounds(integrand_t, t_sym)?;
    // Hyperbolic-rational integrands combine several kernel powers over a
    // common denominator, so their legitimate t-degree is larger than a
    // single exp kernel's: `sinh(x)^4/(1+tanh(x))` needs numerator degree 10.
    let cap = if hyps.is_empty() {
        MAX_KERNEL_DEG
    } else {
        MAX_HYP_DEG
    };
    if num_deg > cap || den_deg > cap {
        return None;
    }
    let result_t = integrate_raw(ctx, integrand_t, t_sym, 0, true, 0, 0);
    if contains_integral(result_t) {
        return None;
    }
    let back = ctx.fun("exp", &[base.arg]);
    Some(normalize(ctx, replace_symbol(ctx, result_t, t_sym, back)))
}

/// Collect the normalized linear exponent arguments `(a, b)` of every
/// `exp` and hyperbolic kernel. Any occurrence of `var` outside such a
/// kernel — a bare `x`, `sin(x)`, `log(x)`, a nonlinear exponent — makes
/// the integrand non-rational in the kernel, so the scan returns `None`.
fn collect_kernels<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    exps: &mut Vec<(Atom<'a>, Atom<'a>)>,
    hyps: &mut Vec<(Atom<'a>, Atom<'a>)>,
) -> Option<()> {
    if is_constant(expr, var) {
        return Some(());
    }
    if exps.len() + hyps.len() >= MAX_KERNELS {
        return None;
    }
    match expr.node() {
        AtomNode::Num(_) => Some(()),
        AtomNode::Var(_) => None,
        AtomNode::Fun(name, args) => {
            let n = name.as_str();
            if args.len() == 1 && (n == "exp" || is_hyperbolic(n)) {
                let (a, b) = linear_form(ctx, args[0], var)?;
                let a = normalize(ctx, a);
                if matches!(a.node(), AtomNode::Num(0)) {
                    return None;
                }
                let b = normalize(ctx, b);
                if n == "exp" {
                    exps.push((a, b));
                } else {
                    hyps.push((a, b));
                }
                return Some(());
            }
            for a in args.iter() {
                collect_kernels(ctx, *a, var, exps, hyps)?;
            }
            Some(())
        }
        AtomNode::Pow(b, e) => {
            // A kernel power records the kernel; the exponent is validated
            // at substitution time (and must not be scanned for `var`).
            if let AtomNode::Fun(name, fargs) = b.node()
                && fargs.len() == 1
                && (name.as_str() == "exp" || is_hyperbolic(name.as_str()))
                && !is_constant(fargs[0], var)
            {
                return collect_kernels(ctx, *b, var, exps, hyps);
            }
            collect_kernels(ctx, *b, var, exps, hyps)?;
            collect_kernels(ctx, *e, var, exps, hyps)
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_kernels(ctx, *a, var, exps, hyps)?;
            }
            Some(())
        }
    }
}

/// Pick the substitution base. Hyperbolic kernels pin the base to their
/// shared linear argument (symbolic coefficients allowed). Pure-exp
/// integrands use the rational gcd of numeric coefficients (so `exp(x)`
/// and `exp(2·x)` share `t = e^x`), or the common symbolic coefficient.
fn decide_base<'a>(
    ctx: &'a AtomArena<'a>,
    exps: &[(Atom<'a>, Atom<'a>)],
    hyps: &[(Atom<'a>, Atom<'a>)],
    var: Symbol,
) -> Option<ExpBase<'a>> {
    let (a, b) = if let Some(&(a0, b0)) = hyps.first() {
        if !hyps.iter().all(|&(a, b)| a == a0 && b == b0) {
            return None;
        }
        (a0, b0)
    } else {
        let a0 = exps.first()?.0;
        let rats: Option<Vec<(i64, i64)>> = exps.iter().map(|&(a, _)| rat_of(a)).collect();
        match rats {
            Some(rs) => {
                let g = rs.iter().map(|&(p, _)| p.unsigned_abs()).reduce(gcd_u64)?;
                let g = i64::try_from(g).ok()?;
                let l = rs.iter().map(|&(_, q)| q).try_fold(1i64, lcm_i64)?;
                let sign = rs.first()?.0.signum();
                (rat_atom(ctx, g.checked_mul(sign)?, l), ctx.num(0))
            }
            None => {
                if !exps.iter().all(|&(a, _)| a == a0) {
                    return None;
                }
                (a0, ctx.num(0))
            }
        }
    };
    let x = ctx.var(var.as_str());
    let arg = normalize(ctx, ctx.add(&[ctx.mul(&[a, x]), b]));
    Some(ExpBase { a, b, arg })
}

/// `a_i / base` as a nonzero integer `|k| ≤ MAX_KERNEL_DEG`. Atom equality
/// (the symbolic-coefficient case) gives `k = 1`; rational coefficients
/// divide exactly by the gcd base construction.
fn kernel_ratio(a_i: Atom<'_>, base: Atom<'_>) -> Option<i64> {
    if a_i == base {
        return Some(1);
    }
    let (p_i, q_i) = rat_of(a_i)?;
    let (p, q) = rat_of(base)?;
    let num = p_i.checked_mul(q)?;
    let den = q_i.checked_mul(p)?;
    if den == 0 || num % den != 0 {
        return None;
    }
    let k = num / den;
    if k == 0 || k.abs() > MAX_KERNEL_DEG {
        return None;
    }
    Some(k)
}

/// Replace every kernel in `expr` by its rational `t`-form. Any leftover
/// occurrence of `var` outside a kernel returns `None`.
fn exp_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    base: &ExpBase<'a>,
    t: Atom<'a>,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) {
        return Some(expr);
    }
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        AtomNode::Var(_) => None,
        AtomNode::Fun(name, args) => {
            let n = name.as_str();
            if args.len() == 1 && n == "exp" {
                let (a_i, b_i) = linear_form(ctx, args[0], var)?;
                let a_i = normalize(ctx, a_i);
                let b_i = normalize(ctx, b_i);
                let k = kernel_ratio(a_i, base.a)?;
                // exp(a_i·x + b_i) = e^c · t^k with c = b_i − k·b.
                let c = normalize(ctx, ctx.add(&[b_i, ctx.mul(&[ctx.num(-k), base.b])]));
                let mut factors = Vec::with_capacity(2);
                if !matches!(c.node(), AtomNode::Num(0)) {
                    factors.push(ctx.fun("exp", &[c]));
                }
                factors.push(int_pow(ctx, t, k));
                return Some(ctx.mul(&factors));
            }
            if args.len() == 1 && is_hyperbolic(n) {
                let (a_i, b_i) = linear_form(ctx, args[0], var)?;
                if normalize(ctx, a_i) != base.a || normalize(ctx, b_i) != base.b {
                    return None;
                }
                return hyp_to_t(ctx, n, t);
            }
            // Any other function of the variable is not kernel-rational.
            None
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, fargs) = b.node()
                && fargs.len() == 1
                && (name.as_str() == "exp" || is_hyperbolic(name.as_str()))
                && !is_constant(fargs[0], var)
            {
                let AtomNode::Num(k) = e.node() else {
                    return None;
                };
                if k.abs() > MAX_KERNEL_DEG {
                    return None;
                }
                let inner = exp_subst(ctx, *b, var, base, t)?;
                return Some(int_pow(ctx, inner, *k));
            }
            if !is_constant(*e, var) {
                return None;
            }
            let nb = exp_subst(ctx, *b, var, base, t)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(exp_subst(ctx, *a, var, base, t)?);
            }
            Some(if matches!(expr.node(), AtomNode::Add(_)) {
                ctx.add(&rebuilt)
            } else {
                ctx.mul(&rebuilt)
            })
        }
    }
}

/// Hyperbolic functions in `t = e^u`, kept as single fractions:
/// `sinh u = (t²−1)/(2t)`, `tanh u = (t²−1)/(t²+1)`, etc.
fn hyp_to_t<'a>(ctx: &'a AtomArena<'a>, name: &str, t: Atom<'a>) -> Option<Atom<'a>> {
    let two_t = ctx.mul(&[ctx.num(2), t]);
    let t2 = ctx.pow(t, ctx.num(2));
    let t2_minus = ctx.add(&[t2, ctx.num(-1)]);
    let t2_plus = ctx.add(&[t2, ctx.num(1)]);
    let r = match name {
        "sinh" => ctx.mul(&[t2_minus, ctx.pow(two_t, ctx.num(-1))]),
        "cosh" => ctx.mul(&[t2_plus, ctx.pow(two_t, ctx.num(-1))]),
        "tanh" => ctx.mul(&[t2_minus, ctx.pow(t2_plus, ctx.num(-1))]),
        "coth" => ctx.mul(&[t2_plus, ctx.pow(t2_minus, ctx.num(-1))]),
        "sech" => ctx.mul(&[two_t, ctx.pow(t2_plus, ctx.num(-1))]),
        "csch" => ctx.mul(&[two_t, ctx.pow(t2_minus, ctx.num(-1))]),
        _ => return None,
    };
    Some(r)
}

fn is_hyperbolic(name: &str) -> bool {
    matches!(name, "sinh" | "cosh" | "tanh" | "coth" | "sech" | "csch")
}

/// Upper bounds `(num_deg, den_deg)` on the `t`-degree of the numerator
/// and denominator of a rational t-form. Returns `None` when `expr` is not
/// rational in `t` (a function applied to a `t`-expression, or a
/// non-integer power of one) — doubling as the post-substitution
/// rationality check.
fn t_degree_bounds(expr: Atom<'_>, t: Symbol) -> Option<(i64, i64)> {
    match expr.node() {
        AtomNode::Num(_) => Some((0, 0)),
        AtomNode::Var(v) => Some(if *v == t { (1, 0) } else { (0, 0) }),
        AtomNode::Add(args) => {
            let mut acc = (0i64, 0i64);
            for a in args.iter() {
                let (n, d) = t_degree_bounds(*a, t)?;
                acc.0 = acc.0.max(n);
                acc.1 = acc.1.max(d);
            }
            Some(acc)
        }
        AtomNode::Mul(args) => {
            let mut acc = (0i64, 0i64);
            for a in args.iter() {
                let (n, d) = t_degree_bounds(*a, t)?;
                acc.0 = acc.0.checked_add(n)?;
                acc.1 = acc.1.checked_add(d)?;
            }
            Some(acc)
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Num(k) = e.node() {
                let (n, d) = t_degree_bounds(*b, t)?;
                return if *k >= 0 {
                    Some((k.checked_mul(n)?, k.checked_mul(d)?))
                } else {
                    let m = k.checked_neg()?;
                    Some((m.checked_mul(d)?, m.checked_mul(n)?))
                };
            }
            let (bn, bd) = t_degree_bounds(*b, t)?;
            let (en, ed) = t_degree_bounds(*e, t)?;
            if bn == 0 && bd == 0 && en == 0 && ed == 0 {
                Some((0, 0))
            } else {
                None
            }
        }
        AtomNode::Fun(_, args) => {
            for a in args.iter() {
                let (n, d) = t_degree_bounds(*a, t)?;
                if n != 0 || d != 0 {
                    return None;
                }
            }
            Some((0, 0))
        }
    }
}

// =========================================================================
// E4: exp of an inverse-function kernel → algebraic rewrite
// =========================================================================

/// Rewrite every `exp(k·F(u))` site (`F` an inverse kernel, `u` linear) into
/// its algebraic equivalent and re-integrate the rewritten integrand.
///
/// Declines when no site matches, when the rewritten integrand exceeds the
/// node budget, when the re-entered chain leaves an `Integral` residue, or
/// when the candidate fails the numeric self-check below.
fn try_exp_inverse<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let mut sites = 0u32;
    let rewritten = rewrite_inverse_exps(ctx, expr, var, &mut sites)?;
    if sites == 0 {
        return None;
    }
    let rewritten = normalize(ctx, rewritten);
    if node_count(rewritten) > MAX_SUBST_NODES {
        return None;
    }
    let result = integrate_raw(ctx, rewritten, var, 0, true, 0, 0);
    if contains_integral(result) {
        return None;
    }
    if !numerically_verified(ctx, rewritten, result, var) {
        return None;
    }
    Some(result)
}

/// Abscissae of the E4 numeric self-check.
const INV_SAMPLES: [f64; 4] = [-1.7, -0.37, 0.83, 2.41];
/// Synthetic real values for the free symbols, applied in sorted-name order.
const INV_PARAMS: [f64; 6] = [0.7, 1.5, 0.5, 1.3, 2.0, 1.1];

/// Deterministic numeric self-check of an E4 candidate: differentiate it and
/// compare with the rewritten integrand at fixed sample points with fixed
/// synthetic parameter values.
///
/// The chain's symbolic-rational backend returns wrong answers for some
/// denominators carrying non-numeric coefficients (`symbolic_rational::
/// rational_square_root`, under separate repair) and the E4 re-entry reaches
/// that path, so the engine's answer is not trusted on its own. A sample is
/// *usable* only when both sides evaluate to finite real values; a candidate
/// with no usable sample — the `i·atan` family, or an integrand that leaves
/// the real domain everywhere — is accepted without a verdict, which is the
/// documented limit of this guard. A single usable disagreement declines.
fn numerically_verified<'a>(
    ctx: &'a AtomArena<'a>,
    rewritten: Atom<'a>,
    candidate: Atom<'a>,
    var: Symbol,
) -> bool {
    let derivative = crate::diff(ctx, candidate, var);
    let mut params: Vec<Symbol> = Vec::new();
    collect_free_symbols(rewritten, var, &mut params);
    collect_free_symbols(derivative, var, &mut params);
    params.sort_unstable_by(|a, b| a.as_str().cmp(b.as_str()));
    params.dedup();
    if params.len() > INV_PARAMS.len() {
        return true;
    }
    let env: Vec<(Symbol, f64)> = params
        .iter()
        .zip(INV_PARAMS.iter())
        .map(|(s, v)| (*s, *v))
        .collect();
    for &xv in INV_SAMPLES.iter() {
        let mut e = env.clone();
        e.push((var, xv));
        // `eval_real` returns `None` for non-finite or out-of-domain values,
        // so a singularity or a complex branch simply drops the sample.
        let (Some(lhs), Some(rhs)) = (eval_real(derivative, &e), eval_real(rewritten, &e)) else {
            continue;
        };
        if (lhs - rhs).abs() > 1e-6 * rhs.abs().max(1.0) {
            return false;
        }
    }
    true
}

/// Real `f64` evaluator for the E4 self-check: `log` maps through the absolute
/// value (antiderivatives carry `log|·|`) and every non-finite or out-of-domain
/// value becomes `None`, so the caller drops the sample instead of comparing
/// meaningless numbers. Any head outside this table — in particular the
/// unresolved imaginary unit — makes the sample unusable.
fn eval_real(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
    let out = match expr.node() {
        AtomNode::Num(n) => *n as f64,
        AtomNode::Var(v) => *env.iter().find(|(s, _)| s == v).map(|(_, x)| x)?,
        AtomNode::Add(args) => {
            let mut acc = 0.0;
            for a in args.iter() {
                acc += eval_real(*a, env)?;
            }
            acc
        }
        AtomNode::Mul(args) => {
            let mut acc = 1.0;
            for a in args.iter() {
                acc *= eval_real(*a, env)?;
            }
            acc
        }
        AtomNode::Pow(b, e) => {
            let base = eval_real(*b, env)?;
            let exp = eval_real(*e, env)?;
            if base < 0.0 && exp.fract() != 0.0 {
                return None;
            }
            base.powf(exp)
        }
        AtomNode::Fun(name, args) => {
            let v = eval_real(*args.first()?, env)?;
            match name.as_str() {
                "sin" => v.sin(),
                "cos" => v.cos(),
                "tan" => v.tan(),
                "cot" => 1.0 / v.tan(),
                "sec" => 1.0 / v.cos(),
                "csc" => 1.0 / v.sin(),
                "exp" => v.exp(),
                "log" => v.abs().ln(),
                "sqrt" => {
                    if v < 0.0 {
                        return None;
                    }
                    v.sqrt()
                }
                "asin" => {
                    if !(-1.0..=1.0).contains(&v) {
                        return None;
                    }
                    v.asin()
                }
                "acos" => {
                    if !(-1.0..=1.0).contains(&v) {
                        return None;
                    }
                    v.acos()
                }
                "atan" => v.atan(),
                "sinh" => v.sinh(),
                "cosh" => v.cosh(),
                "tanh" => v.tanh(),
                "coth" => 1.0 / v.tanh(),
                "asinh" => v.asinh(),
                "acosh" => {
                    if v < 1.0 {
                        return None;
                    }
                    v.acosh()
                }
                "atanh" => {
                    if v.abs() >= 1.0 {
                        return None;
                    }
                    v.atanh()
                }
                "acoth" => {
                    if v.abs() <= 1.0 {
                        return None;
                    }
                    0.5 * ((v + 1.0) / (v - 1.0)).abs().ln()
                }
                _ => return None,
            }
        }
    };
    if out.is_finite() { Some(out) } else { None }
}

/// Free symbols of `expr` other than `var` and the imaginary unit, appended to
/// `out`. The imaginary unit is skipped so that an `i`-carrying integrand gets
/// no synthetic real value (its samples are then unusable, and the check
/// abstains instead of comparing two unrelated real functions).
fn collect_free_symbols(expr: Atom<'_>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Num(_) => {}
        AtomNode::Var(v) => {
            if *v != var && v.as_str() != IMAG_UNIT {
                out.push(*v);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_free_symbols(*a, var, out);
            }
        }
        AtomNode::Pow(b, e) => {
            collect_free_symbols(*b, var, out);
            collect_free_symbols(*e, var, out);
        }
    }
}

/// Structural pass of E4: replace every `exp(k·F(u))` site by its algebraic
/// form, leaving everything else intact. Any occurrence of `var` outside such
/// a site is fine (the rewrite stays an algebraic function of `var`), so this
/// does not need E1's kernel-rationality scan.
fn rewrite_inverse_exps<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sites: &mut u32,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(expr),
        AtomNode::Fun(name, args) => {
            if name.as_str() == "exp"
                && args.len() == 1
                && let Some(algebraic) = exp_inverse_algebraic(ctx, args[0], var)
            {
                *sites += 1;
                if *sites as usize > MAX_KERNELS {
                    return None;
                }
                return Some(algebraic);
            }
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(rewrite_inverse_exps(ctx, *a, var, sites)?);
            }
            Some(ctx.fun(name.as_str(), &rebuilt))
        }
        AtomNode::Add(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(rewrite_inverse_exps(ctx, *a, var, sites)?);
            }
            Some(ctx.add(&rebuilt))
        }
        AtomNode::Mul(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(rewrite_inverse_exps(ctx, *a, var, sites)?);
            }
            Some(ctx.mul(&rebuilt))
        }
        AtomNode::Pow(b, e) => {
            let nb = rewrite_inverse_exps(ctx, *b, var, sites)?;
            let ne = rewrite_inverse_exps(ctx, *e, var, sites)?;
            Some(ctx.pow(nb, ne))
        }
    }
}

/// The algebraic equivalent of `exp(arg)` when `arg = k·F(u)` with `F` one of
/// the inverse kernels, `k` constant in `var` and `u = σ·x + τ` linear with
/// `σ ≠ 0`. Returns `None` for every other argument, including the
/// unsupported exponents — the identity is applied only where it is exact.
fn exp_inverse_algebraic<'a>(
    ctx: &'a AtomArena<'a>,
    arg: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let factors: &[Atom<'a>] = match arg.node() {
        AtomNode::Mul(args) => args,
        _ => std::slice::from_ref(&arg),
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut kernel: Option<(&str, Atom<'a>)> = None;
    for f in factors {
        if let AtomNode::Fun(name, fargs) = f.node()
            && fargs.len() == 1
            && is_inverse_kernel(name.as_str())
            && !is_constant(fargs[0], var)
        {
            // A second (or a mixed) kernel makes the exponent a sum of
            // inverse functions, which is not algebraic: decline.
            if kernel.is_some() {
                return None;
            }
            kernel = Some((name.as_str(), fargs[0]));
            continue;
        }
        if is_constant(*f, var) {
            consts.push(*f);
            continue;
        }
        return None;
    }
    let (name, u) = kernel?;
    // Only a linear kernel argument is in scope: `exp(n·atanh(x²))` declines.
    let (sigma, tau) = linear_form(ctx, u, var)?;
    if matches!(normalize(ctx, sigma).node(), AtomNode::Num(0)) {
        return None;
    }
    let k = if consts.is_empty() {
        ctx.num(1)
    } else {
        normalize(ctx, ctx.mul(&consts))
    };
    let algebraic = match name {
        "atanh" | "acoth" => reciprocal_algebraic(ctx, name, u, k)?,
        "atan" => atan_algebraic(ctx, u, k)?,
        "asinh" => root_algebraic(ctx, u, (sigma, tau), var, k, 1)?,
        "acosh" => root_algebraic(ctx, u, (sigma, tau), var, k, -1)?,
        _ => return None,
    };
    Some(normalize(ctx, algebraic))
}

/// `exp(k·atanh u) = (1+u)^(k/2)·(1−u)^(−k/2)` and
/// `exp(k·acoth u) = (1+u)^(k/2)·(u−1)^(−k/2)`.
///
/// The identities are exponent-linear and exact for every real `k`, so a
/// symbolic `k` (the corpus `exp((n·atanh(a·x)))` shape) is accepted; a `k`
/// carrying the imaginary unit is not (it would put `i` in the exponent and
/// the "algebraic" rewrite would no longer be a function of `var` alone).
/// The two powers are emitted *split* rather than as one `((1+u)/(1−u))^(k/2)`
/// site so that `normalize` can merge them with structurally identical
/// `(1±u)`-factors elsewhere in the integrand — that merging is what turns the
/// corpus shapes into rational functions.
fn reciprocal_algebraic<'a>(
    ctx: &'a AtomArena<'a>,
    name: &str,
    u: Atom<'a>,
    k: Atom<'a>,
) -> Option<Atom<'a>> {
    if mentions_imaginary(k) {
        return None;
    }
    let one_plus = normalize(ctx, ctx.add(&[ctx.num(1), u]));
    let other = if name == "acoth" {
        normalize(ctx, ctx.add(&[u, ctx.num(-1)]))
    } else {
        normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), u])]))
    };
    let half_k = half_of(ctx, k)?;
    let neg_half_k = normalize(ctx, ctx.mul(&[ctx.num(-1), half_k]));
    Some(ctx.mul(&[ctx.pow(one_plus, half_k), ctx.pow(other, neg_half_k)]))
}

/// `k/2`, folded to a numeric literal when `k` is rational. The fold matters:
/// the chain's integer-exponent engines do not recognise `Mul(2, Pow(2,−1))`
/// as `1`, so an unfolded `k/2` makes `exp(2·atanh u) = (1+u)/(1−u)` look
/// like a symbolic power and every engine declines it. A symbolic `k` keeps
/// the formal `k/2` (and is declined downstream).
fn half_of<'a>(ctx: &'a AtomArena<'a>, k: Atom<'a>) -> Option<Atom<'a>> {
    match rat_of(k) {
        Some((p, q)) => {
            if p.abs() > 2 * MAX_INV_EXP_NUM || q > 2 * MAX_INV_EXP_DEN {
                return None;
            }
            Some(rat_atom(ctx, p, q.checked_mul(2)?))
        }
        None => Some(normalize(ctx, ctx.mul(&[k, inv(ctx, ctx.num(2))]))),
    }
}

/// `exp(i·r·atan u) = (1+i·u)^(r/2)·(1−i·u)^(−r/2)` with `r = k/i` rational.
///
/// Only the exact power is produced; a `k` that is not a rational multiple of
/// the imaginary unit (`exp(2·atan u)`, `exp(n·atan u)`) declines, as does an
/// out-of-budget `r`.
fn atan_algebraic<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, k: Atom<'a>) -> Option<Atom<'a>> {
    let (p, q) = imaginary_multiple(ctx, k)?;
    if p == 0 || p.abs() > MAX_INV_EXP_NUM || q > MAX_INV_EXP_DEN {
        return None;
    }
    let half_r = rat_atom(ctx, p, q.checked_mul(2)?);
    if matches!(half_r.node(), AtomNode::Num(0)) {
        return None;
    }
    let i = ctx.var(IMAG_UNIT);
    let iu = normalize(ctx, ctx.mul(&[i, u]));
    let one_plus = normalize(ctx, ctx.add(&[ctx.num(1), iu]));
    let one_minus = normalize(ctx, ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), iu])]));
    let neg_half_r = normalize(ctx, ctx.mul(&[ctx.num(-1), half_r]));
    Some(ctx.mul(&[ctx.pow(one_plus, half_r), ctx.pow(one_minus, neg_half_r)]))
}

/// `exp(k·asinh u) = (u + √(u²+1))^k` and
/// `exp(k·acosh u) = (u + √(u²−1))^k` for rational `k` (`e0 = ±1` picks the
/// radicand constant). A symbolic or `i`-carrying `k` declines: `(…)^k` is
/// then not an algebraic function.
fn root_algebraic<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    lin: (Atom<'a>, Atom<'a>),
    var: Symbol,
    k: Atom<'a>,
    e0: i64,
) -> Option<Atom<'a>> {
    let (p, q) = rat_of(k)?;
    if p == 0 || p.abs() > MAX_INV_EXP_NUM || q > MAX_INV_EXP_DEN {
        return None;
    }
    let radicand = expanded_square_plus(ctx, lin, var, e0);
    let root = ctx.fun("sqrt", &[radicand]);
    let base = ctx.add(&[u, root]);
    Some(ctx.pow(base, rat_atom(ctx, p, q)))
}

/// `(σ·x + τ)² + e0`, expanded to `σ²·x² + 2στ·x + τ² + e0`.
///
/// The radical engines match a *polynomial* radicand, so `Pow(Mul(σ,x),2)` is
/// not recognised; spelling `u²` out as a polynomial is what makes
/// `exp(k·asinh(a·x))` reachable.
fn expanded_square_plus<'a>(
    ctx: &'a AtomArena<'a>,
    lin: (Atom<'a>, Atom<'a>),
    var: Symbol,
    e0: i64,
) -> Atom<'a> {
    let (sigma, tau) = lin;
    let x = ctx.var(var.as_str());
    let mut terms = vec![
        ctx.num(e0),
        ctx.mul(&[ctx.pow(sigma, ctx.num(2)), ctx.pow(x, ctx.num(2))]),
    ];
    if !matches!(normalize(ctx, tau).node(), AtomNode::Num(0)) {
        terms.push(ctx.mul(&[ctx.num(2), sigma, tau, x]));
        terms.push(ctx.pow(tau, ctx.num(2)));
    }
    ctx.add(&terms)
}

/// The rational `r` when the constant `k` is exactly `i·r` with `i` the
/// imaginary-unit symbol; `None` for every other constant (a bare `r`, an
/// `i²`, a symbolic multiple, or a product of several `i` factors).
fn imaginary_multiple<'a>(ctx: &'a AtomArena<'a>, k: Atom<'a>) -> Option<(i64, i64)> {
    let factors: &[Atom<'a>] = match k.node() {
        AtomNode::Mul(args) => args,
        _ => std::slice::from_ref(&k),
    };
    let mut i_count = 0u32;
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        if matches!(f.node(), AtomNode::Var(v) if v.as_str() == IMAG_UNIT) {
            i_count += 1;
            continue;
        }
        rest.push(*f);
    }
    if i_count != 1 {
        return None;
    }
    if rest.is_empty() {
        return Some((1, 1));
    }
    rat_of(normalize(ctx, ctx.mul(&rest)))
}

/// True when the `i` symbol occurs anywhere in `expr`.
fn mentions_imaginary(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => v.as_str() == IMAG_UNIT,
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|a| mentions_imaginary(*a))
        }
        AtomNode::Pow(b, e) => mentions_imaginary(*b) || mentions_imaginary(*e),
    }
}

/// The inverse functions whose exponential is algebraic in their argument.
fn is_inverse_kernel(name: &str) -> bool {
    matches!(name, "atanh" | "acoth" | "atan" | "asinh" | "acosh")
}

// =========================================================================
// Small local utilities
// =========================================================================

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification. Antiderivatives of
    /// real integrands carry `log|·|`; evaluating `log` through the
    /// absolute value makes the check robust to the sign convention of
    /// the partial-fraction form.
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
                    "log" => v.abs().ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "tanh" => v.tanh(),
                    _ => return None,
                })
            }
        }
    }

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, s).unwrap()
    }

    /// Run the mechanism directly, require `Some` without residue, and
    /// check `diff(result) == integrand` numerically at the sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_exp_log(ctx, integrand, var).expect("mechanism declined");
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

    /// Minimal complex `f64` for the `exp(i·k·atan u)` family: `i` is a
    /// *symbol* in the corpus expressions, so the derivative check has to run
    /// in ℂ where `i` evaluates to the imaginary unit.
    #[derive(Clone, Copy, Debug)]
    struct C64 {
        re: f64,
        im: f64,
    }

    impl C64 {
        const fn num(re: f64) -> Self {
            Self { re, im: 0.0 }
        }

        const fn imag_unit() -> Self {
            Self { re: 0.0, im: 1.0 }
        }

        fn add(self, o: Self) -> Self {
            Self {
                re: self.re + o.re,
                im: self.im + o.im,
            }
        }

        fn sub(self, o: Self) -> Self {
            Self {
                re: self.re - o.re,
                im: self.im - o.im,
            }
        }

        fn mul(self, o: Self) -> Self {
            Self {
                re: self.re * o.re - self.im * o.im,
                im: self.re * o.im + self.im * o.re,
            }
        }

        fn div(self, o: Self) -> Self {
            let d = o.re * o.re + o.im * o.im;
            Self {
                re: (self.re * o.re + self.im * o.im) / d,
                im: (self.im * o.re - self.re * o.im) / d,
            }
        }

        fn abs(self) -> f64 {
            self.re.hypot(self.im)
        }

        /// Integer powers by repeated multiplication: exact, and free of the
        /// branch ambiguity of `exp(e·log z)`.
        fn ipow(self, n: i64) -> Self {
            let base = if n < 0 {
                Self::num(1.0).div(self)
            } else {
                self
            };
            let mut acc = Self::num(1.0);
            for _ in 0..n.unsigned_abs() {
                acc = acc.mul(base);
            }
            acc
        }

        fn exp(self) -> Self {
            let r = self.re.exp();
            Self {
                re: r * self.im.cos(),
                im: r * self.im.sin(),
            }
        }

        fn log(self) -> Self {
            if self.im == 0.0 && self.re > 0.0 {
                return Self::num(self.re.ln());
            }
            Self {
                re: self.abs().ln(),
                im: self.im.atan2(self.re),
            }
        }

        fn sqrt(self) -> Self {
            if self.im == 0.0 && self.re >= 0.0 {
                return Self::num(self.re.sqrt());
            }
            let r = self.abs().sqrt();
            let theta = self.im.atan2(self.re) * 0.5;
            Self {
                re: r * theta.cos(),
                im: r * theta.sin(),
            }
        }

        fn sin(self) -> Self {
            Self {
                re: self.re.sin() * self.im.cosh(),
                im: self.re.cos() * self.im.sinh(),
            }
        }

        fn cos(self) -> Self {
            Self {
                re: self.re.cos() * self.im.cosh(),
                im: -self.re.sin() * self.im.sinh(),
            }
        }

        fn sinh(self) -> Self {
            Self {
                re: self.re.sinh() * self.im.cos(),
                im: self.re.cosh() * self.im.sin(),
            }
        }

        fn cosh(self) -> Self {
            Self {
                re: self.re.cosh() * self.im.cos(),
                im: self.re.sinh() * self.im.sin(),
            }
        }

        fn pow(self, e: Self) -> Self {
            if e.im == 0.0 && e.re.fract() == 0.0 && e.re.abs() <= 64.0 {
                return self.ipow(e.re as i64);
            }
            if e.im == 0.0 && self.im == 0.0 && self.re > 0.0 {
                return Self::num(self.re.powf(e.re));
            }
            self.log().mul(e).exp()
        }
    }

    /// `eval_f64` in ℂ, so the imaginary unit can take its actual value and
    /// antiderivatives carrying `i` (or a complex branch) can be checked.
    /// Returns `None` for functions outside the evaluator, which makes the
    /// caller treat the case as inconclusive rather than failed.
    fn eval_c64(expr: Atom<'_>, env: &[(Symbol, C64)]) -> Option<C64> {
        match expr.node() {
            AtomNode::Num(n) => Some(C64::num(*n as f64)),
            AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, x)| *x),
            AtomNode::Add(args) => args
                .iter()
                .try_fold(C64::num(0.0), |acc, a| Some(acc.add(eval_c64(*a, env)?))),
            AtomNode::Mul(args) => args
                .iter()
                .try_fold(C64::num(1.0), |acc, a| Some(acc.mul(eval_c64(*a, env)?))),
            AtomNode::Pow(b, e) => {
                let base = eval_c64(*b, env)?;
                match e.node() {
                    AtomNode::Num(n) => Some(base.ipow(*n)),
                    _ => Some(base.pow(eval_c64(*e, env)?)),
                }
            }
            AtomNode::Fun(name, args) => {
                let v = eval_c64(*args.first()?, env)?;
                let one = C64::num(1.0);
                let half = C64::num(0.5);
                let i = C64::imag_unit();
                Some(match name.as_str() {
                    "exp" => v.exp(),
                    "log" => v.log(),
                    "sqrt" => v.sqrt(),
                    "sin" => v.sin(),
                    "cos" => v.cos(),
                    "tan" => v.sin().div(v.cos()),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "tanh" => v.sinh().div(v.cosh()),
                    // atan z = (i/2)·(log(1 − i·z) − log(1 + i·z)).
                    "atan" => i
                        .div(C64::num(2.0))
                        .mul(one.sub(i.mul(v)).log().sub(one.add(i.mul(v)).log())),
                    "asin" => i
                        .mul(one.sub(v.mul(v)).sqrt().add(i.mul(v)))
                        .log()
                        .mul(C64::num(-1.0)),
                    "acos" => C64::num(std::f64::consts::FRAC_PI_2).sub(
                        i.mul(one.sub(v.mul(v)).sqrt().add(i.mul(v)))
                            .log()
                            .mul(C64::num(-1.0)),
                    ),
                    "atanh" => half.mul(one.add(v).log().sub(one.sub(v).log())),
                    "acoth" => half.mul(v.add(one).log().sub(v.sub(one).log())),
                    "asinh" => v.add(v.mul(v).add(one).sqrt()).log(),
                    "acosh" => v.add(v.mul(v).sub(one).sqrt()).log(),
                    _ => return None,
                })
            }
        }
    }

    /// Parameter environment shared by the corpus-shape checks.
    fn corpus_env() -> Vec<(Symbol, C64)> {
        vec![
            (Symbol::new("a"), C64::num(0.7)),
            (Symbol::new("b"), C64::num(0.5)),
            (Symbol::new("c"), C64::num(1.5)),
            (Symbol::new("n"), C64::num(1.7)),
            (Symbol::new("m"), C64::num(2.0)),
            (Symbol::new("p"), C64::num(1.3)),
            (Symbol::new("i"), C64::imag_unit()),
        ]
    }

    /// Sample points per kernel domain: `atanh` needs `|a·x| < 1`, `acoth`
    /// `|a·x| > 1`, `acosh(a + b·x)` `a + b·x > 1`, the rest are unrestricted.
    fn samples_for(domain: &str) -> &'static [f64] {
        match domain {
            "acoth" => &[-2.6, 2.1, 3.0],
            "acosh" => &[1.2, 2.0, 3.0],
            // `acosh(a·x)` with `a` as small as 0.4 still needs `a·x ≥ 1`.
            "acosh0" => &[2.6, 3.5, 4.4],
            _ => &[-0.8, 0.35, 1.1],
        }
    }

    /// Run the full pipeline on a corpus shape and, when it produces an
    /// answer, verify `diff(answer) == integrand` numerically (complex
    /// arithmetic, so the `i·atan` family and off-domain branches are covered
    /// too). Returns whether an answer was produced.
    ///
    /// The public `integrate` entry is used rather than a direct
    /// `integrate_exp_log` call because it resets the thread-local chain
    /// budget: a long loop of direct calls would exhaust that budget and
    /// degrade later shapes to fallback for reasons unrelated to the shape.
    fn check_shape<'a>(ctx: &'a AtomArena<'a>, integrand: &str, domain: &str) -> bool {
        let var = Symbol::new("x");
        let expr = parse(ctx, integrand);
        let result = crate::integrate(ctx, expr, var);
        if contains_integral(result) {
            return false;
        }
        let d = crate::diff(ctx, result, var);
        let env = corpus_env();
        for &xv in samples_for(domain) {
            let mut e = env.clone();
            e.push((var, C64::num(xv)));
            let Some(lhs) = eval_c64(d, &e) else {
                continue;
            };
            let Some(rhs) = eval_c64(expr, &e) else {
                continue;
            };
            let tol = 1e-5 * rhs.abs().max(1.0);
            assert!(
                lhs.sub(rhs).abs() < tol,
                "{integrand} at x={xv}: diff={lhs:?} integrand={rhs:?} (result: {result})"
            );
        }
        true
    }

    /// Run `f` on a thread with a 32 MiB stack.
    ///
    /// The integrator recurses deeply on some corpus shapes; libtest's default
    /// test-thread stack (and the 1 MiB Windows main-thread stack used by
    /// `--test-threads=1`) overflows, which aborts the whole test binary
    /// instead of failing a single test.
    fn with_big_stack<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
        std::thread::Builder::new()
            .stack_size(32 << 20)
            .spawn(f)
            .expect("spawn test thread")
            .join()
            .expect("test thread panicked")
    }

    /// `(id, integrand, domain)` for the corpus shapes the mechanism is
    /// measured on. The full 65-shape sweep is run with the corpus harness
    /// (`OCAS_1892_CASES`); this is the bounded unit-test subset: every shape
    /// the mechanism unlocks plus the representative declines (the
    /// quadratic-denominator shape the numeric guard rejects, and the
    /// symbolic-exponent / mixed-function / `i`-family declinations).
    const TARGET_SHAPES: &[(&str, &str, &str)] = &[
        // Solved by E4 (verified numerically below).
        (
            "rubi-00196",
            "exp((2*atanh(a*x)))*x*(c - a^2*c*x^2)^2",
            "atanh",
        ),
        ("rubi-00359", "exp((2*atanh(a*x)))*(c - a*c*x)^2", "atanh"),
        (
            "rubi-00614",
            "1/(exp((2*atanh(a*x)))*(c - a*c*x)^(3/2))",
            "atanh",
        ),
        (
            "rubi-00662",
            "exp((2*acoth(a*x)))*(c - c/(a^2*x^2))",
            "acoth",
        ),
        (
            "rubi-00684",
            "exp((2*atanh(a*x)))*(c - a^2*c*x^2)^3/x",
            "atanh",
        ),
        ("rubi-00789", "exp((2*atanh(a*x)))*(c - a*c*x)^5", "atanh"),
        ("rubi-00980", "1/(exp((2*i*atan(a*x)))*x^2)", "iatan"),
        (
            "rubi-01137",
            "(c - c/(a^2*x^2))/exp((2*atanh(a*x)))",
            "atanh",
        ),
        ("rubi-01539", "exp((2*acoth(a*x)))/(c - c/(a*x))^2", "acoth"),
        ("rubi-01551", "exp((2*acoth(a*x)))/(c - c/(a*x))^4", "acoth"),
        ("rubi-01760", "exp((2*acoth(a*x)))*(c - a*c*x)^5", "acoth"),
        ("rubi-01890", "exp((4*i*atan(a*x)))", "iatan"),
        // Declined: the wrong-answer shape the guard must keep out.
        (
            "rubi-01096",
            "exp((2*atanh(a*x)))*x^2/(c - a^2*c*x^2)",
            "atanh",
        ),
        // Declined: symbolic exponents, quadratic radicands, symbolic `m`.
        ("rubi-01223", "exp((n*atanh(a*x)))/x^4", "atanh"),
        (
            "rubi-00471",
            "exp((3*atanh(a*x)))*(c - a^2*c*x^2)^p",
            "atanh",
        ),
        (
            "rubi-00692",
            "exp((3*acoth(a*x)))/(c - a^2*c*x^2)^4",
            "acoth",
        ),
        ("rubi-00780", "exp((3*atanh(a*x)))/(c - c/(a*x))^3", "atanh"),
        (
            "rubi-00899",
            "sqrt(c - a^2*c*x^2)/(exp((2*acoth(a*x)))*x^4)",
            "acoth",
        ),
        ("rubi-01018", "x^m/exp((4*i*atan(a*x)))", "iatan"),
        ("rubi-00727", "exp(asinh(a + b*x))*x^2", "asinh"),
        ("rubi-00856", "exp(acosh(a + b*x))*x", "acosh"),
    ];

    #[test]
    fn exp_rational_corpus() {
        // Corpus fallback: ∫ e^x/(1 − e^{2x}) dx → t = e^x gives ∫ dt/(1 − t²).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "exp(x)/(1-exp(2*x))");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn exp_rational_symbolic_coeff() {
        // ∫ e^{a·x}/(1 + e^{a·x}) dx = log(1 + e^{a·x})/a, symbolic a.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "exp(a*x)/(1+exp(a*x))");
        let env = [(Symbol::new("a"), 1.7)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.2, 0.6, 1.1]);
    }

    #[test]
    fn hyperbolic_sinh_over_cosh() {
        // Corpus shape: ∫ sinh(x)/(2 + 3·cosh(x)) dx = log(2 + 3·cosh(x))/3.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sinh(x)/(2+3*cosh(x))");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.8, 1.5]);
    }

    #[test]
    fn hyperbolic_tanh_symbolic() {
        // Corpus shape: ∫ dx/(a + b·tanh(x)), rationalized by t = e^x.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "1/(a+b*tanh(x))");
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 3.0)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.9, 1.7]);
    }

    #[test]
    fn log_kernel_subst_atan() {
        // Corpus fallback: ∫ dx/(x·(1 + log(x)²)) = atan(log(x)). Checked
        // on the raw parse and on the normalized shape.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for expr in [
            parse(&ctx, "1/(x*(1+log(x)^2))"),
            normalize(&ctx, parse(&ctx, "1/(x*(1+log(x)^2))")),
        ] {
            assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.5, 1.5, 3.0]);
        }
    }

    #[test]
    fn log_kernel_subst_function_of_kernel() {
        // ∫ sin(log(x))/x dx = −cos(log(x)): functions of the kernel pass.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sin(log(x))/x");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.4, 1.2, 2.5]);
    }

    #[test]
    fn log_pow_over_x() {
        // Corpus shape: ∫ log(x)³/x dx = log(x)⁴/4 (the B6 m = −1 gap).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "log(x)^3/x");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.5, 1.5, 3.0]);
    }

    #[test]
    fn inv_x_log_pow() {
        // ∫ dx/(x·log(x)²) = −log(x)⁻¹.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "1/(x*log(x)^2)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.5, 2.0, 3.0]);
    }

    #[test]
    fn log_power_expand_linear_base() {
        // ∫ log(2·(3 + 4·x)²) dx → formal expansion 2·log(3+4x) + log(2).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "log(2*(3+4*x)^2)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.2, 0.7, 1.3]);
    }

    #[test]
    fn log_power_expand_monomial() {
        // ∫ log(x³) dx → 3·log(x) → 3·x·(log(x) − 1).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "log(x^3)");
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.5, 1.0, 2.5]);
    }

    #[test]
    fn decline_exp_times_sin() {
        // B3's shape; the structural kernel scan rejects the sin factor
        // (routing exp·sin here once caused parts loops, the 0.27 lesson).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "exp(x)*sin(x)");
        assert!(integrate_exp_log(&ctx, expr, Symbol::new("x")).is_none());
        // Mixed hyperbolic/trig likewise declines.
        let expr2 = parse(&ctx, "sinh(x)*cos(x)");
        assert!(integrate_exp_log(&ctx, expr2, Symbol::new("x")).is_none());
    }

    #[test]
    fn decline_b2_shape_and_nonlinear_exponent() {
        // B2's shape: a bare `x` factor is not rational in the kernel.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "x*exp(x)");
        assert!(integrate_exp_log(&ctx, expr, Symbol::new("x")).is_none());
        // Nonlinear exponent: no linear kernel base exists.
        let expr2 = parse(&ctx, "exp(x^2)/(1+exp(x^2))");
        assert!(integrate_exp_log(&ctx, expr2, Symbol::new("x")).is_none());
    }

    #[test]
    fn decline_log_power_budget() {
        // log(x)^9/x exceeds the log-power budget (8).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "log(x)^9/x");
        assert!(integrate_exp_log(&ctx, expr, Symbol::new("x")).is_none());
    }

    #[test]
    fn exp_inverse_rational_collapse() {
        // `exp(2*atanh(a*x)) = (1+u)/(1-u)` is rational, so the shapes whose
        // remaining factors carry matching `(1+-u)`-powers collapse to a
        // rational function and integrate. This is the corpus set the
        // mechanism actually unlocks (measured on the 1892 benchmark).
        // Big stack: these are the deeply-recursing symbolic-rational shapes.
        with_big_stack(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            for (s, domain) in [
                ("exp((2*atanh(a*x)))*(c - a*c*x)^2", "atanh"), // 00359
                ("exp((2*atanh(a*x)))*(c - a*c*x)^5", "atanh"), // 00789
                ("exp((2*atanh(a*x)))*x*(c - a^2*c*x^2)^2", "atanh"), // 00196
                ("exp((2*atanh(a*x)))*(c - a^2*c*x^2)^3/x", "atanh"), // 00684
                ("(c - c/(a^2*x^2))/exp((2*atanh(a*x)))", "atanh"), // 01137
                ("1/(exp((2*atanh(a*x)))*(c - a*c*x)^(3/2))", "atanh"), // 00614
                ("exp((2*acoth(a*x)))*(c - c/(a^2*x^2))", "acoth"), // 00662
                ("exp((2*acoth(a*x)))/(c - c/(a*x))^2", "acoth"), // 01539
                ("exp((2*acoth(a*x)))/(c - c/(a*x))^4", "acoth"), // 01551
                ("exp((2*acoth(a*x)))*(c - a*c*x)^5", "acoth"), // 01760
            ] {
                assert!(check_shape(&ctx, s, domain), "expected E4 to solve {s}");
            }
        });
    }

    #[test]
    fn exp_inverse_mechanism_fires() {
        // The stage itself (not a later engine) produces the answer: three
        // direct calls keep the chain budget comfortably in range.
        with_big_stack(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x = Symbol::new("x");
            for s in [
                "exp((2*atanh(a*x)))*(c - a*c*x)^2",
                "exp((2*acoth(a*x)))/(c - c/(a*x))^2",
                "exp((4*i*atan(a*x)))",
            ] {
                let expr = parse(&ctx, s);
                let r = integrate_exp_log(&ctx, expr, x);
                assert!(r.is_some(), "E4 declined {s}");
            }
        });
    }

    #[test]
    fn exp_inverse_numeric_guard() {
        // The guard is what keeps a wrong engine answer out of the corpus
        // (the `symbolic_rational` quadratic-denominator bug is reachable from
        // the E4 re-entry): a correct candidate passes, a wrong one is
        // rejected.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        let f = parse(&ctx, "x^2");
        assert!(numerically_verified(&ctx, f, parse(&ctx, "(3^-1)*x^3"), x));
        assert!(!numerically_verified(&ctx, f, parse(&ctx, "x^3"), x));
        assert!(!numerically_verified(
            &ctx,
            f,
            parse(&ctx, "(3^-1)*x^3 + x"),
            x
        ));
    }

    #[test]
    fn exp_inverse_atanh_asinh_acosh_identities() {
        // Real identities at several parameter values, checked numerically
        // through the mechanism itself. Big stack: 36 full-chain integrations
        // of the deeply-recursing symbolic-rational shapes.
        with_big_stack(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let var = Symbol::new("x");
            for a in [0.4, 0.7, 1.3] {
                for c in [1.0, 1.5, 2.25] {
                    let env = [
                        (Symbol::new("a"), C64::num(a)),
                        (Symbol::new("c"), C64::num(c)),
                    ];
                    for (s, domain) in [
                        ("exp((2*atanh(a*x)))*(c - a*c*x)^2", "atanh"),
                        ("exp((2*acoth(a*x)))/(c - c/(a*x))^2", "acoth"),
                        ("exp(asinh(a*x))*x", "asinh"),
                        ("exp(acosh(a*x))*x", "acosh0"),
                    ] {
                        let expr = parse(&ctx, s);
                        // Public entry: resets the chain budget between rounds.
                        let result = crate::integrate(&ctx, expr, var);
                        assert!(
                            !contains_integral(result),
                            "declined {s} (a={a}, c={c}): {result}"
                        );
                        let d = crate::diff(&ctx, result, var);
                        for &xv in samples_for(domain) {
                            let mut e = env.to_vec();
                            e.push((var, C64::num(xv)));
                            let lhs = eval_c64(d, &e).expect("eval diff");
                            let rhs = eval_c64(expr, &e).expect("eval integrand");
                            assert!(
                                lhs.sub(rhs).abs() < 1e-6 * rhs.abs().max(1.0),
                                "{s} (a={a}, c={c}) at x={xv}: diff={lhs:?} integrand={rhs:?}"
                            );
                        }
                    }
                }
            }
        });
    }

    #[test]
    fn exp_inverse_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        for s in [
            // Non-linear kernel argument: out of scope by construction.
            "exp(n*atanh(x^2))",
            "exp(2*acoth(x^3))",
            "exp(asinh(x^2))*x",
            // Not an inverse kernel at all.
            "exp(sin(x))",
            "exp(x*atanh(x))",
            // Mixed kernels: the exponent is a sum, not a single kernel.
            "exp(atanh(x) + atan(x))",
            "exp(atanh(x))*exp(acoth(x))",
            // Unsupported exponents: no `i` factor (atan), symbolic
            // asinh/acosh coefficients.
            "exp(2*atan(x))",
            "exp(n*atan(x))*x",
            "exp(n*asinh(x))",
            "exp(n*acosh(1 + x))",
            // Imaginary multiple whose half is not an exact power.
            "exp(i*atan(x))",
            "exp((3*i*atan(x)))*x",
        ] {
            let expr = parse(&ctx, s);
            assert!(
                integrate_exp_log(&ctx, expr, x).is_none(),
                "expected a decline for {s}"
            );
        }
    }

    #[test]
    fn exp_inverse_iatan_family() {
        // `exp(2*i*atan u) = (1+i*u)/(1-i*u)`: the rewrite is a rational
        // function of x with `i` as a symbolic constant, checked in C.
        with_big_stack(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            for (s, domain) in [
                ("exp((4*i*atan(a*x)))", "iatan"),
                ("1/(exp((2*i*atan(a*x)))*x^2)", "iatan"),
                ("exp((-2*i*atan(a*x)))*x^3", "iatan"),
            ] {
                assert!(check_shape(&ctx, s, domain), "expected E4 to solve {s}");
            }
            // The symbolic-`m` corpus shape stays out of reach: `x^m` is not a
            // rational function and no engine accepts it.
            let expr = parse(&ctx, "x^m/exp((4*i*atan(a*x)))");
            assert!(integrate_exp_log(&ctx, expr, Symbol::new("x")).is_none());
        });
    }

    #[test]
    fn exp_inverse_repeat_stability() {
        // 300 repeats on both the declined and the solved path: the outcome
        // must be identical every round (no leakage through the stage, the
        // chain budget or the arena) and the answer must not grow.
        //
        // The decline loop calls the stage directly: those shapes produce no
        // rewrite, so no chain re-entry happens and nothing accumulates. The
        // solved loop needs the public entry, which resets the thread-local
        // chain budget per call.
        {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x = Symbol::new("x");
            for s in [
                "exp(sin(x))",
                "exp(n*atanh(x^2))",
                "exp(atanh(x) + atan(x))",
                "exp(2*atan(x))*x",
            ] {
                let expr = parse(&ctx, s);
                for i in 0..300 {
                    assert!(
                        integrate_exp_log(&ctx, expr, x).is_none(),
                        "{s} round {i} unexpectedly succeeded"
                    );
                }
            }
        }
        // The solved rounds run on a big stack (see `with_big_stack`).
        with_big_stack(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x = Symbol::new("x");
            let expr = parse(&ctx, "1/(exp((2*atanh(a*x)))*(c - a*c*x)^(3/2))");
            let mut first_len = None;
            for i in 0..300 {
                let r = crate::integrate(&ctx, expr, x);
                assert!(!contains_integral(r), "round {i}: {r}");
                let len = r.to_string().len();
                match first_len {
                    None => first_len = Some(len),
                    Some(f) => assert_eq!(f, len, "round {i}: answer grew"),
                }
                assert!(node_count(r) < 4096, "round {i}: answer blow-up");
            }
        });
    }

    /// Every corpus shape the mechanism targets: each answer it does produce
    /// must pass the numerical derivative check (the guard against the
    /// symbolic-rational quadratic-denominator bug), and declinations are
    /// allowed for the shapes no engine can finish.
    #[test]
    fn exp_inverse_targets_verified_or_declined() {
        // Runs on a big stack: the integrator recurses deeply on several of
        // these shapes (see `with_big_stack`).
        with_big_stack(|| {
            let mut solved = 0usize;
            let mut declined = Vec::new();
            // A fresh arena per shape, mirroring the corpus harness (one child
            // process per case).
            for (id, integrand, domain) in TARGET_SHAPES {
                let arena = Arena::new();
                let ctx = AtomArena::new(&arena);
                if check_shape(&ctx, integrand, domain) {
                    solved += 1;
                } else {
                    declined.push(*id);
                }
            }
            println!(
                "E4 solved {solved}/{} shapes; declined: {declined:?}",
                TARGET_SHAPES.len()
            );
        });
    }
}
