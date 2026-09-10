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
//!
//! All chain re-entries go through `integrate_raw` and are declined on any
//! `Integral` residue. The substituted t-forms are rational in `t` (E1) or
//! kernel-free (E2), so the module cannot re-match its own output.

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
/// Node budget for the input integrand and for the substituted t-form.
const MAX_SUBST_NODES: usize = 200;
/// Cap on the number of exp/hyperbolic/log kernel sites in one integrand.
const MAX_KERNELS: usize = 16;
/// `|n|` cap for the `log(c·g^n)` formal power expansion.
const MAX_LOG_EXPAND: i64 = 4;

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
    if num_deg > MAX_KERNEL_DEG || den_deg > MAX_KERNEL_DEG {
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
}
