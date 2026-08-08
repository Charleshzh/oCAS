//! Heuristic integration techniques.
//!
//! This module provides integration strategies that go beyond the lookup table,
//! Risch algorithm, and special-function forms. The techniques are:
//!
//! - **Integration by parts** (LIATE heuristic)
//! - **Trigonometric substitution** (`√(a²−x²)`, `√(a²+x²)`, `√(x²−a²)`)
//! - **Weierstrass substitution** (`t = tan(u/2)` for trig-rational integrands)
//! - **Euler substitution** (`√(ax²+bx+c)`)

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use crate::derivative::diff;
use crate::integral::{integrate_raw, is_constant, is_fallback};

/// Maximum recursion depth for integration by parts.
const PARTS_MAX_DEPTH: usize = 2;

// =========================================================================
// Public API
// =========================================================================

/// Try heuristic integration techniques on `expr`.
///
/// `parts_depth` counts integration-by-parts expansions performed so far
/// (threaded through [`integrate_raw`]); it bounds the parts loop even when
/// the integrand alternates between two forms (e.g. `exp(ax)sin(bx)` ↔
/// `exp(ax)cos(bx)`), which previously reset the parts depth on every
/// re-entry and looped forever. Legitimate multi-level parts chains
/// (`x²sin(x)` → `2x·cos(x)` → `−2sin(x)`) stay within the budget.
///
/// Returns `Some(result)` if any technique succeeds, `None` otherwise.
pub(crate) fn heuristic_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    parts_depth: usize,
) -> Option<Atom<'a>> {
    // 1. Integration by parts
    if let Some(r) = try_parts(ctx, expr, var, parts_depth) {
        if !is_fallback(&r) {
            return Some(r);
        }
    }
    // 2. Trigonometric substitution
    let ts = try_trig_substitution(ctx, expr, var);
    if let Some(r) = ts {
        if !is_fallback(&r) {
            return Some(r);
        }
    }
    // 3. Weierstrass substitution
    if let Some(r) = try_weierstrass(ctx, expr, var) {
        if !is_fallback(&r) {
            return Some(r);
        }
    }
    // 4. Euler substitution
    if let Some(r) = try_euler_substitution(ctx, expr, var) {
        if !is_fallback(&r) {
            return Some(r);
        }
    }
    None
}

// =========================================================================
// LIATE scoring
// =========================================================================

/// LIATE priority score for integration by parts.
///
/// Lower score → higher priority to be chosen as `u`.
/// - Log = 0, Inverse trig = 1, Algebraic = 2, Trig = 3, Exponential = 4, Other = 5
fn liate_score(expr: &Atom<'_>) -> u32 {
    match expr.node() {
        AtomNode::Fun(name, _) => match name.as_str() {
            "log" | "ln" => 0,
            "asin" | "acos" | "atan" | "acot" | "asec" | "acsc" | "asinh" | "acosh" | "atanh" => 1,
            "sin" | "cos" | "tan" | "sec" | "csc" | "cot" | "sinh" | "cosh" | "tanh" => 3,
            "exp" => 4,
            _ => 5,
        },
        AtomNode::Pow(base, exp) => {
            // x^n is algebraic (score 2)
            if matches!(base.node(), AtomNode::Var(_))
                && matches!(exp.node(), AtomNode::Num(n) if *n > 0)
            {
                2
            } else {
                // Recurse into base for scoring
                liate_score(base).min(5)
            }
        }
        AtomNode::Var(_) => 2, // plain variable is algebraic
        AtomNode::Mul(args) => args.iter().map(liate_score).min().unwrap_or(5),
        _ => 5,
    }
}

// =========================================================================
// Integration by parts (1a)
// =========================================================================

/// Try integration by parts using LIATE heuristic.
///
/// For `∫ u·v' dx`-type products, picks `u` by LIATE priority,
/// computes `v = ∫ v' dx`, then applies `∫ u·v' = u·V - ∫ u'·V`.
fn try_parts<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    parts_depth: usize,
) -> Option<Atom<'a>> {
    if parts_depth >= PARTS_MAX_DEPTH {
        return None;
    }

    // Collect multiplicative factors
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => return None,
    };

    if factors.len() < 2 {
        return None;
    }

    // Split into constant and non-constant factors
    let mut constants = Vec::new();
    let mut non_constants = Vec::new();
    for f in &factors {
        if is_constant(*f, var) {
            constants.push(*f);
        } else {
            non_constants.push(*f);
        }
    }

    if non_constants.len() < 2 {
        return None;
    }

    // Find the factor with the lowest LIATE score (best candidate for u)
    let (u_idx, _) = non_constants
        .iter()
        .enumerate()
        .min_by_key(|(_, e)| liate_score(e))
        .unwrap();

    let u = non_constants[u_idx];
    let v_prime_factors: Vec<Atom<'a>> = non_constants
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != u_idx)
        .map(|(_, e)| *e)
        .collect();
    let v_prime = if v_prime_factors.len() == 1 {
        v_prime_factors[0]
    } else {
        ctx.mul(&v_prime_factors)
    };

    // Compute V = ∫ v' dx (recursive, at higher depth)
    let v = integrate_raw(ctx, v_prime, var, parts_depth + 2, true, 0, parts_depth + 1);
    if is_fallback(&v) {
        return None;
    }

    // Compute u' = du/dx
    let u_prime = diff(ctx, u, var);
    if is_fallback(&u_prime) {
        return None;
    }

    // Compute u' * V
    let u_prime_v = ctx.mul(&[u_prime, v]);

    // Compute ∫ u' * V dx (recursive)
    let integral_u_prime_v = integrate_raw(
        ctx,
        u_prime_v,
        var,
        parts_depth + 2,
        true,
        0,
        parts_depth + 1,
    );

    // Build result: u * V - ∫ u' * V
    let u_times_v = ctx.mul(&[u, v]);

    // Combine with constant factors
    let core_result = if is_fallback(&integral_u_prime_v) {
        // If inner integral failed, try anyway with the parts formula
        // but only if the result is simpler
        return None;
    } else {
        ctx.add(&[u_times_v, ctx.mul(&[ctx.num(-1), integral_u_prime_v])])
    };

    if constants.is_empty() {
        Some(core_result)
    } else {
        let mut result_factors = constants;
        result_factors.push(core_result);
        Some(ctx.mul(&result_factors))
    }
}

// =========================================================================
// Trigonometric substitution (1b)
// =========================================================================

/// Check if exponent is 1/2 (square root marker).
#[allow(dead_code)]
fn is_half_exponent(e: &AtomNode) -> bool {
    match e {
        // Direct 2^(-1) form (arena normalizes 1 * x → x)
        AtomNode::Pow(b, exp) => {
            matches!(b.node(), AtomNode::Num(2)) && matches!(exp.node(), AtomNode::Num(-1))
        }
        // Mul(1, 2^(-1)) form
        AtomNode::Mul(args) if args.len() == 2 => {
            let has_half = args.iter().any(|a| {
                matches!(a.node(), AtomNode::Pow(b, e)
                if matches!(b.node(), AtomNode::Num(2)) && matches!(e.node(), AtomNode::Num(-1)))
            });
            let has_one = args.iter().any(|a| matches!(a.node(), AtomNode::Num(1)));
            has_half && has_one
        }
        _ => false,
    }
}

/// Check if exponent is -1/2 (1/sqrt marker).
fn is_neg_half_exponent(e: &AtomNode) -> bool {
    match e {
        // Mul(-1, 2^(-1)) form
        AtomNode::Mul(args) if args.len() == 2 => {
            let has_neg = args.iter().any(|a| matches!(a.node(), AtomNode::Num(-1)));
            let has_half = args.iter().any(|a| {
                matches!(a.node(), AtomNode::Pow(b, e)
                if matches!(b.node(), AtomNode::Num(2)) && matches!(e.node(), AtomNode::Num(-1)))
            });
            has_neg && has_half
        }
        _ => false,
    }
}

/// Match `positive_term - x²` pattern, returning `a` from `a²`.
fn match_subtracted_squares<'a>(
    ctx: &'a AtomArena<'a>,
    positive: Atom<'a>,
    negative: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // negative should be -x² or -1*x²
    let x_squared = match negative.node() {
        // -1 * x^2
        AtomNode::Mul(args) if args.len() == 2 => {
            if matches!(args[0].node(), AtomNode::Num(-1)) {
                args[1]
            } else if matches!(args[1].node(), AtomNode::Num(-1)) {
                args[0]
            } else {
                return None;
            }
        }
        _ => return None,
    };

    // Check x_squared is x^2
    if !matches!(x_squared.node(), AtomNode::Pow(b, e)
        if matches!(b.node(), AtomNode::Var(v) if *v == var) && matches!(e.node(), AtomNode::Num(2)))
    {
        return None;
    }

    // positive should be a² — extract a
    extract_square_base(ctx, positive, var)
}

/// Match `a² + x²` pattern, returning `a`.
fn match_sum_squares<'a>(
    ctx: &'a AtomArena<'a>,
    t1: Atom<'a>,
    t2: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // One term is x², the other is a² (constant w.r.t. var)
    let is_x2 = |a: &Atom<'a>| {
        matches!(a.node(), AtomNode::Pow(b, e)
            if matches!(b.node(), AtomNode::Var(v) if *v == var) && matches!(e.node(), AtomNode::Num(2)))
    };

    if is_x2(&t1) && is_constant(t2, var) {
        extract_square_base(ctx, t2, var)
    } else if is_x2(&t2) && is_constant(t1, var) {
        extract_square_base(ctx, t1, var)
    } else {
        None
    }
}

/// Match `x² - a²` pattern, returning `a`.
fn match_x_minus_a_squared<'a>(
    ctx: &'a AtomArena<'a>,
    positive: Atom<'a>,
    negative: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // positive should be x²
    if !matches!(positive.node(), AtomNode::Pow(b, e)
        if matches!(b.node(), AtomNode::Var(v) if *v == var) && matches!(e.node(), AtomNode::Num(2)))
    {
        return None;
    }

    // negative should be -a²
    let a_squared = match negative.node() {
        AtomNode::Mul(args) if args.len() == 2 => {
            if matches!(args[0].node(), AtomNode::Num(-1)) {
                args[1]
            } else if matches!(args[1].node(), AtomNode::Num(-1)) {
                args[0]
            } else {
                return None;
            }
        }
        _ => return None,
    };

    if is_constant(a_squared, var) {
        extract_square_base(ctx, a_squared, var)
    } else {
        None
    }
}

/// Extract `a` from `a²`.
fn extract_square_base<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if matches!(e.node(), AtomNode::Num(2)) && is_constant(*b, var) {
                Some(*b)
            } else {
                None
            }
        }
        AtomNode::Num(n) => {
            if *n > 0 {
                let sqrt = (*n as f64).sqrt();
                if sqrt == sqrt.floor() {
                    Some(ctx.num(sqrt as i64))
                } else {
                    None
                }
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Try trigonometric substitution on `expr`.
///
/// Instead of performing a full symbolic substitution + integration (which
/// requires trig identity simplification we don't have), we match the most
/// common patterns directly and return their known antiderivatives.
fn try_trig_substitution<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // Pattern: 1/√(a² - x²) → asin(x/a)
    if let Some(antideriv) = match_inv_sqrt_a2_minus_x2(ctx, expr, var) {
        return Some(antideriv);
    }
    // Pattern: 1/√(a² + x²) → asinh(x/a)  [= log(x + √(a²+x²))]
    if let Some(antideriv) = match_inv_sqrt_a2_plus_x2(ctx, expr, var) {
        return Some(antideriv);
    }
    // Pattern: 1/√(x² - a²) → acosh(x/a)  [= log(x + √(x²-a²))]
    if let Some(antideriv) = match_inv_sqrt_x2_minus_a2(ctx, expr, var) {
        return Some(antideriv);
    }
    // Pattern: √(a² - x²) → (x·√(a²-x²) + a²·asin(x/a)) / 2
    if let Some(antideriv) = match_sqrt_a2_minus_x2(ctx, expr, var) {
        return Some(antideriv);
    }
    // Pattern: √(a² + x²) → (x·√(a²+x²) + a²·asinh(x/a)) / 2
    if let Some(antideriv) = match_sqrt_a2_plus_x2(ctx, expr, var) {
        return Some(antideriv);
    }
    // Pattern: √(x² - a²) → (x·√(x²-a²) - a²·acosh(x/a)) / 2
    if let Some(antideriv) = match_sqrt_x2_minus_a2(ctx, expr, var) {
        return Some(antideriv);
    }

    None
}

/// Match `1/√(a² - x²)` → `asin(x/a)`.
fn match_inv_sqrt_a2_minus_x2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // expr = Pow(inner, -1/2) where inner = a² - x²
    let base = get_sqrt_base(expr, var)?;
    if let Some(_a) = match_subtracted_squares_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        if is_one_atom(a) {
            return Some(ctx.fun("asin", &[x]));
        }
        return Some(ctx.fun("asin", &[ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])]));
    }
    None
}

/// Match `1/√(a² + x²)` → `asinh(x/a)`.
fn match_inv_sqrt_a2_plus_x2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let base = get_sqrt_base(expr, var)?;
    if let Some(_a) = match_sum_squares_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        if is_one_atom(a) {
            return Some(ctx.fun("asinh", &[x]));
        }
        return Some(ctx.fun("asinh", &[ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])]));
    }
    None
}

/// Match `1/√(x² - a²)` → `acosh(x/a)`.
fn match_inv_sqrt_x2_minus_a2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let base = get_sqrt_base(expr, var)?;
    if let Some(_a) = match_x_minus_a_sq_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        if is_one_atom(a) {
            return Some(ctx.fun("acosh", &[x]));
        }
        return Some(ctx.fun("acosh", &[ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])]));
    }
    None
}

/// Match `√(a² - x²)` → `(x·√(a²-x²) + a²·asin(x/a)) / 2`.
fn match_sqrt_a2_minus_x2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // expr = Pow(inner, 1/2) where inner = a² - x²
    let base = get_sqrt_base_positive(expr, var)?;
    if let Some(_a) = match_subtracted_squares_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        let sqrt_part = ctx.mul(&[x, expr]); // x · √(a²-x²)
        let a_sq = ctx.pow(a, ctx.num(2));
        let asin_arg = if is_one_atom(a) {
            x
        } else {
            ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])
        };
        let asin_part = ctx.mul(&[a_sq, ctx.fun("asin", &[asin_arg])]);
        return Some(ctx.mul(&[
            ctx.add(&[sqrt_part, asin_part]),
            ctx.pow(ctx.num(2), ctx.num(-1)),
        ]));
    }
    None
}

/// Match `√(a² + x²)` → `(x·√(a²+x²) + a²·asinh(x/a)) / 2`.
fn match_sqrt_a2_plus_x2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let base = get_sqrt_base_positive(expr, var)?;
    if let Some(_a) = match_sum_squares_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        let sqrt_part = ctx.mul(&[x, expr]);
        let a_sq = ctx.pow(a, ctx.num(2));
        let asinh_arg = if is_one_atom(a) {
            x
        } else {
            ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])
        };
        let asinh_part = ctx.mul(&[a_sq, ctx.fun("asinh", &[asinh_arg])]);
        return Some(ctx.mul(&[
            ctx.add(&[sqrt_part, asinh_part]),
            ctx.pow(ctx.num(2), ctx.num(-1)),
        ]));
    }
    None
}

/// Match `√(x² - a²)` → `(x·√(x²-a²) - a²·acosh(x/a)) / 2`.
fn match_sqrt_x2_minus_a2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let base = get_sqrt_base_positive(expr, var)?;
    if let Some(_a) = match_x_minus_a_sq_pattern(ctx, base, var) {
        let x = ctx.var(var.as_str());
        let a = _a;
        let sqrt_part = ctx.mul(&[x, expr]);
        let a_sq = ctx.pow(a, ctx.num(2));
        let acosh_arg = if is_one_atom(a) {
            x
        } else {
            ctx.mul(&[x, ctx.pow(a, ctx.num(-1))])
        };
        let acosh_part = ctx.mul(&[a_sq, ctx.fun("acosh", &[acosh_arg])]);
        return Some(ctx.mul(&[
            ctx.add(&[sqrt_part, ctx.mul(&[ctx.num(-1), acosh_part])]),
            ctx.pow(ctx.num(2), ctx.num(-1)),
        ]));
    }
    None
}

// --- Helper: extract the base of a sqrt or 1/sqrt ---

/// If expr is `(something)^(-1/2)` or `(something)^(-1·2^(-1))`, return `something`.
/// Also handles `(something^(1/2))^(-1)` which is `1/sqrt(something)`.
fn get_sqrt_base<'a>(expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if is_neg_half_exponent(e.node()) {
                return Some(*b);
            }
            // Handle (sqrt(inner))^(-1) = 1/sqrt(inner)
            if matches!(e.node(), AtomNode::Num(-1)) {
                if let Some(inner) = get_sqrt_base_positive(*b, var) {
                    return Some(inner);
                }
            }
            None
        }
        _ => None,
    }
}

/// If expr is `(something)^(1/2)` or `(something)^(2^(-1))`, return `something`.
fn get_sqrt_base_positive<'a>(expr: Atom<'a>, _var: Symbol) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if is_half_exponent(e.node()) {
                Some(*b)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Check if expr is the number 1.
fn is_one_atom(expr: Atom<'_>) -> bool {
    matches!(expr.node(), AtomNode::Num(1))
}

/// Match inner as `positive - x²` pattern (returns `a` from `a²`).
fn match_subtracted_squares_pattern<'a>(
    ctx: &'a AtomArena<'a>,
    inner: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let args = match inner.node() {
        AtomNode::Add(a) if a.len() == 2 => a,
        _ => return None,
    };
    let (t1, t2) = (args[0], args[1]);
    match_subtracted_squares(ctx, t1, t2, var)
        .or_else(|| match_subtracted_squares(ctx, t2, t1, var))
}

/// Match inner as `a² + x²` pattern (returns `a`).
fn match_sum_squares_pattern<'a>(
    ctx: &'a AtomArena<'a>,
    inner: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let args = match inner.node() {
        AtomNode::Add(a) if a.len() == 2 => a,
        _ => return None,
    };
    let (t1, t2) = (args[0], args[1]);
    match_sum_squares(ctx, t1, t2, var)
}

/// Match inner as `x² - a²` pattern (returns `a`).
fn match_x_minus_a_sq_pattern<'a>(
    ctx: &'a AtomArena<'a>,
    inner: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let args = match inner.node() {
        AtomNode::Add(a) if a.len() == 2 => a,
        _ => return None,
    };
    let (t1, t2) = (args[0], args[1]);
    match_x_minus_a_squared(ctx, t1, t2, var).or_else(|| match_x_minus_a_squared(ctx, t2, t1, var))
}

/// Substitute every occurrence of `var` in `expr` with `replacement`.
fn substitute_atom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Var(v) => {
            if *v == var {
                replacement
            } else {
                expr
            }
        }
        AtomNode::Num(_) => expr,
        AtomNode::Add(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, var, replacement))
                .collect();
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, var, replacement))
                .collect();
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = substitute_atom(ctx, *base, var, replacement);
            let new_exp = substitute_atom(ctx, *exp, var, replacement);
            ctx.pow(new_base, new_exp)
        }
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, var, replacement))
                .collect();
            ctx.fun(name.as_str(), &new_args)
        }
    }
}

// =========================================================================
// Weierstrass substitution (1c)
// =========================================================================

/// Check if `expr` is a rational function of `sin(u)` and `cos(u)` where `u` is linear in `var`.
fn is_trig_rational<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        // A bare occurrence of the integration variable (outside sin/cos) is
        // NOT trig-rational: Weierstrass would treat it as a constant of the
        // t-integral and produce a wrong answer.
        AtomNode::Var(_) => is_constant(expr, var),
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|a| is_trig_rational(*a, var)),
        AtomNode::Pow(base, exp) => {
            // Allow integer powers
            matches!(exp.node(), AtomNode::Num(_)) && is_trig_rational(*base, var)
        }
        AtomNode::Fun(name, args) => {
            let n = name.as_str();
            if n == "sin" || n == "cos" {
                // Argument must be linear in var (or just var itself)
                args.len() == 1 && is_linear_in(args[0], var)
            } else {
                // Other functions: only if entirely constant
                is_constant(expr, var)
            }
        }
    }
}

/// Check if `expr` is linear in `var` (i.e., of the form `a·var + b`).
fn is_linear_in<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Var(v) => *v == var,
        AtomNode::Num(_) => true,
        AtomNode::Add(args) => args
            .iter()
            .all(|a| is_linear_in(*a, var) || is_constant(*a, var)),
        AtomNode::Mul(args) => {
            let mut has_var = false;
            for a in args.iter() {
                if let AtomNode::Var(v) = a.node() {
                    if *v == var {
                        has_var = true;
                    }
                } else if !is_constant(*a, var) {
                    return false;
                }
            }
            has_var
        }
        _ => false,
    }
}

/// Count sin/cos function applications anywhere in `expr`.
///
/// Weierstrass must only apply to expressions that actually contain trig
/// functions: a "trig-rational" with none (e.g. the t-integrand of a
/// previous Weierstrass pass, which is vacuously rational in `_t`) would be
/// substituted again — multiplying by `dx/dt` once more — and loop forever.
fn trig_count(expr: Atom<'_>) -> u32 {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            let n = name.as_str();
            let base = if n == "sin" || n == "cos" { 1 } else { 0 };
            base + args.iter().map(|a| trig_count(*a)).sum::<u32>()
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().map(|a| trig_count(*a)).sum(),
        AtomNode::Pow(base, exp) => trig_count(*base) + trig_count(*exp),
        AtomNode::Num(_) | AtomNode::Var(_) => 0,
    }
}

/// True when every occurrence of `var` in `expr` sits inside a `sin`/`cos`
/// call whose argument is exactly `var` (the only shapes
/// [`substitute_trig`] replaces).
///
/// `is_trig_rational` accepts linear arguments (`cos(c+d*x)`), but the
/// substitution does not replace them — the integration variable would then
/// survive inside the t-integrand as a "constant", and the result would be
/// a bogus `x·f(x)`-style antiderivative (observed as
/// `2·atan(tan(x/2))·f(x)`). Rejecting such expressions keeps Weierstrass
/// sound.
/// The common linear argument `u = a·x + b` of every `sin`/`cos` in
/// `expr`, with every occurrence of `var` inside such a call (the only
/// shape Weierstrass can substitute). Returns `None` when `var` appears
/// outside a substituted position (a bare factor, a non-trig function, or
/// sin/cos with different arguments) — substituting then would treat `var`
/// as a constant of the t-integral and produce a wrong `x·f(x)`-style
/// antiderivative.
fn trig_linear_arg(expr: Atom<'_>, var: Symbol) -> Option<Atom<'_>> {
    fn visit<'a>(expr: Atom<'a>, var: Symbol, found: &mut Option<Atom<'a>>) -> Option<()> {
        match expr.node() {
            AtomNode::Fun(name, args) if args.len() == 1 => {
                let n = name.as_str();
                if n == "sin" || n == "cos" {
                    let u = args[0];
                    if !contains_var(u, var) {
                        return Some(());
                    }
                    match found {
                        Some(prev) if *prev != u => return None,
                        _ => *found = Some(u),
                    }
                    Some(())
                } else {
                    // Other functions must not contain var.
                    if contains_var(expr, var) {
                        return None;
                    }
                    Some(())
                }
            }
            AtomNode::Fun(_, _) => {
                if contains_var(expr, var) {
                    return None;
                }
                Some(())
            }
            AtomNode::Add(args) | AtomNode::Mul(args) => {
                for a in args.iter() {
                    visit(*a, var, found)?;
                }
                Some(())
            }
            AtomNode::Pow(base, exp) => {
                visit(*base, var, found)?;
                visit(*exp, var, found)?;
                Some(())
            }
            AtomNode::Num(_) => Some(()),
            AtomNode::Var(v) => {
                if *v == var {
                    return None;
                }
                Some(())
            }
        }
    }
    let mut found: Option<Atom<'_>> = None;
    visit(expr, var, &mut found)?;
    found
}

/// True if `expr` contains `var` anywhere.
fn contains_var(expr: Atom<'_>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Var(v) => *v == var,
        AtomNode::Num(_) => false,
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_var(*a, var)),
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_var(*a, var)),
        AtomNode::Pow(base, exp) => contains_var(*base, var) || contains_var(*exp, var),
    }
}

/// Replace `sin(u)`/`cos(u)` in `expr` with the given t-forms; any other
/// occurrence of `var` (or a different trig argument) aborts with `None`.
fn substitute_trig_arg<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    u: Atom<'a>,
    var: Symbol,
    sin_u: Atom<'a>,
    cos_u: Atom<'a>,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Fun(name, args) if args.len() == 1 => {
            let n = name.as_str();
            if n == "sin" || n == "cos" {
                if args[0] == u {
                    return Some(if n == "sin" { sin_u } else { cos_u });
                }
                return None;
            }
            if contains_var(expr, var) {
                return None;
            }
            Some(expr)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => {
            if contains_var(expr, var) {
                return None;
            }
            Some(expr)
        }
        AtomNode::Add(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(substitute_trig_arg(ctx, *a, u, var, sin_u, cos_u)?);
            }
            Some(ctx.add(&rebuilt))
        }
        AtomNode::Mul(args) => {
            let mut rebuilt = Vec::with_capacity(args.len());
            for a in args.iter() {
                rebuilt.push(substitute_trig_arg(ctx, *a, u, var, sin_u, cos_u)?);
            }
            Some(ctx.mul(&rebuilt))
        }
        AtomNode::Pow(base, exp) => {
            let b = substitute_trig_arg(ctx, *base, u, var, sin_u, cos_u)?;
            let e = substitute_trig_arg(ctx, *exp, u, var, sin_u, cos_u)?;
            Some(ctx.pow(b, e))
        }
        AtomNode::Fun(_, _) => {
            if contains_var(expr, var) {
                return None;
            }
            Some(expr)
        }
    }
}

/// Try Weierstrass substitution: `t = tan(u/2)` for the linear argument
/// `u = a·x + b`.
fn try_weierstrass<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if !is_trig_rational(expr, var) {
        return None;
    }
    // No sin/cos → nothing to substitute; applying the substitution would
    // loop on the result.
    if trig_count(expr) == 0 {
        return None;
    }
    let u = trig_linear_arg(expr, var)?;
    let (a, _b) = crate::integral::linear_form(ctx, u, var)?;
    if matches!(a.node(), AtomNode::Num(0)) {
        return None;
    }
    // Symbolic (non-numeric) `a` is allowed: the t-integral then has
    // symbolic coefficients, handled by the symbolic-constant rational
    // backend.
    let t = ctx.var("_t");
    let one_plus_t2 = ctx.add(&[ctx.num(1), ctx.pow(t, ctx.num(2))]);
    let sin_u = ctx.mul(&[ctx.num(2), t, ctx.pow(one_plus_t2, ctx.num(-1))]);
    let cos_u = ctx.mul(&[
        ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(t, ctx.num(2))])]),
        ctx.pow(one_plus_t2, ctx.num(-1)),
    ]);
    // u = a·x + b ⇒ dx = du/a, and du = 2/(1+t²) dt.
    let dx_dt = ctx.mul(&[
        ctx.num(2),
        ctx.pow(a, ctx.num(-1)),
        ctx.pow(one_plus_t2, ctx.num(-1)),
    ]);

    let substituted = substitute_trig_arg(ctx, expr, u, var, sin_u, cos_u)?;
    let integrand = ctx.mul(&[substituted, dx_dt]);

    // Feasibility gate: if the t-integrand is beyond the symbolic rational
    // backend's reach, skip Weierstrass entirely — feeding it to the rest
    // of the chain grinds (rational/parts on the t-form) without solving.
    let t_sym = Symbol::new("_t");
    if !crate::integral::symbolic_rational::rational_complexity_ok(ctx, integrand, t_sym) {
        return None;
    }
    let result_t = integrate_raw(ctx, integrand, t_sym, 2, true, 0, 0);
    if is_fallback(&result_t) {
        return None;
    }

    // Back-substitute t = tan(u/2).
    let back = ctx.fun("tan", &[ctx.mul(&[u, ctx.pow(ctx.num(2), ctx.num(-1))])]);
    Some(substitute_atom(ctx, result_t, t_sym, back))
}

// =========================================================================
// Euler substitution (1d)
// =========================================================================

/// Substitute `t` for the square-root term in `expr`: after substituting
/// `x → x(t)`, the square root `√(a·x² + b·x + c)` becomes a rational
/// expression in `t`. Every `Fun(sqrt, [q])` / `q^(±1/2)` node whose
/// argument equals the substituted quadratic is replaced by that rational
/// expression; other square roots are left untouched (and make the t-form
/// non-rational, so [`crate::integral::rational::integrate_rational`] then
/// declines and the substitution is abandoned).
fn replace_sqrt<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    quad_t: Atom<'a>,
    sqrt_t: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            let q = ocas_atom::normalize::normalize(ctx, args[0]);
            if q == quad_t { sqrt_t } else { expr }
        }
        AtomNode::Pow(base, exp) => {
            // sqrt(q)^-1 and q^(±1/2) forms.
            if let (AtomNode::Fun(name, args), AtomNode::Num(-1)) = (base.node(), exp.node())
                && name.as_str() == "sqrt"
                && args.len() == 1
                && ocas_atom::normalize::normalize(ctx, args[0]) == quad_t
            {
                return ctx.pow(sqrt_t, ctx.num(-1));
            }
            let b = replace_sqrt(ctx, *base, quad_t, sqrt_t);
            let e = replace_sqrt(ctx, *exp, quad_t, sqrt_t);
            ctx.pow(b, e)
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt(ctx, *a, quad_t, sqrt_t))
                .collect();
            ctx.add(&rebuilt)
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt(ctx, *a, quad_t, sqrt_t))
                .collect();
            ctx.mul(&rebuilt)
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_sqrt(ctx, *a, quad_t, sqrt_t))
                .collect();
            ctx.fun(name.as_str(), &rebuilt)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Try Euler substitution for `√(a·x² + b·x + c)` with rational
/// coefficients.
///
/// - **Euler I** (a is a positive rational square): `t = √(ax²+bx+c) − √a·x`
///   makes the integrand a rational function of `t`.
/// - **Euler II** (c is a positive rational square, a not): `t = (√(...) − √c)/x`.
///
/// Non-rational coefficients or no eligible square root → `None` (the
/// heuristic chain continues). `integrate_rational` never re-enters the
/// heuristic, so there is no recursion risk.
fn try_euler_substitution<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (a, b, c) = crate::integral::rules::quadratic_coeffs(ctx, expr, var)?;
    let (pa, qa) = crate::integral::rules::rat_of(a)?;
    let (pc, qc) = crate::integral::rules::rat_of(c)?;

    let t = ctx.var("_t");
    let x = ctx.var(var.as_str());

    // Rational value helper: build Atom(p)/Atom(q) as p*q^-1.
    let rat = |ctx: &'a AtomArena<'a>, p: i64, q: i64| -> Atom<'a> {
        ctx.mul(&[ctx.num(p), ctx.pow(ctx.num(q), ctx.num(-1))])
    };

    // Euler I: a = (r/s)² > 0.
    let (x_t, sqrt_t, dx_dt): (Atom<'a>, Atom<'a>, Atom<'a>) =
        if let Some((r, s)) = crate::integral::rules::rational_sqrt(pa, qa) {
            let sq_a = rat(ctx, r, s);
            // x = (t² − c)/(b − 2√a·t)
            let denom = ctx.add(&[b, ctx.mul(&[ctx.num(-2), sq_a, t])]);
            let x_t = ctx.mul(&[
                ctx.add(&[ctx.pow(t, ctx.num(2)), ctx.mul(&[ctx.num(-1), c])]),
                ctx.pow(denom, ctx.num(-1)),
            ]);
            // √(ax²+bx+c) = √a·x + t (rational in t after x → x(t))
            let sqrt_t = ctx.add(&[ctx.mul(&[sq_a, x_t]), t]);
            // dx/dt = [2t(b−2√a·t) + 2√a·(t²−c)] / (b−2√a·t)²
            let num = ctx.add(&[
                ctx.mul(&[ctx.num(2), t, denom]),
                ctx.mul(&[
                    ctx.mul(&[ctx.num(2), sq_a]),
                    ctx.add(&[ctx.pow(t, ctx.num(2)), ctx.mul(&[ctx.num(-1), c])]),
                ]),
            ]);
            let dx_dt = ctx.mul(&[num, ctx.pow(ctx.pow(denom, ctx.num(2)), ctx.num(-1))]);
            (x_t, sqrt_t, dx_dt)
        } else if let Some((r, s)) = crate::integral::rules::rational_sqrt(pc, qc) {
            // Euler II: c = (r/s)² > 0, a not a square.
            let sq_c = rat(ctx, r, s);
            // x = (2√c·t − b)/(a − t²)
            let denom = ctx.add(&[a, ctx.mul(&[ctx.num(-1), ctx.pow(t, ctx.num(2))])]);
            let x_t = ctx.mul(&[
                ctx.add(&[ctx.mul(&[ctx.num(2), sq_c, t]), ctx.mul(&[ctx.num(-1), b])]),
                ctx.pow(denom, ctx.num(-1)),
            ]);
            // √(ax²+bx+c) = √c + t·x
            let sqrt_t = ctx.add(&[sq_c, ctx.mul(&[t, x_t])]);
            // dx/dt = [2√c(a−t²) + (2√c t − b)·2t] / (a−t²)²
            let num = ctx.add(&[
                ctx.mul(&[ctx.num(2), sq_c, denom]),
                ctx.mul(&[
                    ctx.add(&[ctx.mul(&[ctx.num(2), sq_c, t]), ctx.mul(&[ctx.num(-1), b])]),
                    ctx.mul(&[ctx.num(2), t]),
                ]),
            ]);
            let dx_dt = ctx.mul(&[num, ctx.pow(ctx.pow(denom, ctx.num(2)), ctx.num(-1))]);
            (x_t, sqrt_t, dx_dt)
        } else {
            return None;
        };

    // Substitute x → x(t) throughout, then replace the sqrt node.
    let substituted = substitute_atom(ctx, expr, var, x_t);
    // The substituted quadratic, normalized (structural match target).
    let quad_t = ocas_atom::normalize::normalize(
        ctx,
        ctx.add(&[
            ctx.mul(&[a, ctx.pow(x_t, ctx.num(2))]),
            ctx.mul(&[b, x_t]),
            c,
        ]),
    );
    let substituted = replace_sqrt(ctx, substituted, quad_t, sqrt_t);
    let integrand = ctx.mul(&[substituted, dx_dt]);

    let t_sym = Symbol::new("_t");
    let result_t = crate::integral::rational::integrate_rational(ctx, integrand, t_sym)?;

    // Back-substitute t = √(a·x²+b·x+c) − √a·x (Euler I) or
    // t = (√(a·x²+b·x+c) − √c)/x (Euler II): rebuild from the original
    // quadratic and the x-expression.
    let x_back = x;
    let sqrt_back = if crate::integral::rules::rational_sqrt(pa, qa).is_some() {
        let (r, s) = crate::integral::rules::rational_sqrt(pa, qa).unwrap();
        let sq_a = rat(ctx, r, s);
        let quad = ctx.add(&[
            ctx.mul(&[a, ctx.pow(x_back, ctx.num(2))]),
            ctx.mul(&[b, x_back]),
            c,
        ]);
        let sqrt = ctx.fun("sqrt", &[quad]);
        ctx.add(&[sqrt, ctx.mul(&[ctx.num(-1), sq_a, x_back])])
    } else {
        let (r, s) = crate::integral::rules::rational_sqrt(pc, qc).unwrap();
        let sq_c = rat(ctx, r, s);
        let quad = ctx.add(&[
            ctx.mul(&[a, ctx.pow(x_back, ctx.num(2))]),
            ctx.mul(&[b, x_back]),
            c,
        ]);
        let sqrt = ctx.fun("sqrt", &[quad]);
        // t = (√(..) − √c)/x
        ctx.mul(&[
            ctx.add(&[sqrt, ctx.mul(&[ctx.num(-1), sq_c])]),
            ctx.pow(x_back, ctx.num(-1)),
        ])
    };
    let back = substitute_atom(ctx, result_t, t_sym, sqrt_back);
    Some(ocas_atom::normalize::normalize(ctx, back))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrate;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn parts_x_exp() {
        // ∫ x·exp(x) dx — Risch handles this, but parts should also work
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.fun("exp", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn parts_x_sin() {
        // ∫ x·sin(x) dx = sin(x) - x·cos(x)
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.fun("sin", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn parts_x2_sin() {
        // ∫ x²·sin(x) dx
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.pow(x, ctx.num(2)), ctx.fun("sin", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn parts_log() {
        // ∫ log(x) dx = x·log(x) - x (parts with u=log, v'=1)
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("log", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn parts_x_log() {
        // ∫ x·log(x) dx = x²·log(x)/2 - x²/4
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.fun("log", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn trig_sub_asin_direct() {
        // ∫ 1/√(1-x²) dx = asin(x)
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let inner = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let half = ctx.pow(ctx.num(2), ctx.num(-1));
        let sqrt_expr = ctx.pow(inner, half);
        let expr = ctx.pow(sqrt_expr, ctx.num(-1));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "asin(x)");
    }

    #[test]
    fn trig_sub_sqrt_1_minus_x2() {
        // ∫ √(1-x²) dx
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let inner = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let sqrt_expr = ctx.pow(
            inner,
            ctx.mul(&[ctx.num(1), ctx.pow(ctx.num(2), ctx.num(-1))]),
        );
        let result = integrate(&ctx, sqrt_expr, Symbol::new("x"));
        // May or may not succeed depending on simplification
        let _ = result;
    }

    #[test]
    fn weierstrass_1_over_sin_plus_1() {
        // ∫ 1/(sin(x)+1) dx
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(ctx.add(&[ctx.fun("sin", &[x]), ctx.num(1)]), ctx.num(-1));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        // This is a challenging integral; may return unevaluated
        let _ = result;
    }

    #[test]
    fn weierstrass_1_over_2_plus_cos() {
        // ∫ 1/(2+cos(x)) dx
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(ctx.add(&[ctx.num(2), ctx.fun("cos", &[x])]), ctx.num(-1));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let _ = result;
    }

    #[test]
    fn liate_scoring() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");

        assert_eq!(liate_score(&ctx.fun("log", &[x])), 0);
        assert_eq!(liate_score(&ctx.fun("asin", &[x])), 1);
        assert_eq!(liate_score(&x), 2);
        assert_eq!(liate_score(&ctx.fun("sin", &[x])), 3);
        assert_eq!(liate_score(&ctx.fun("exp", &[x])), 4);
    }

    #[test]
    fn heuristic_none_for_unknown() {
        // All heuristics should return None for a completely unknown function
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("unknown_func", &[x]);
        let result = heuristic_integrate(&ctx, expr, Symbol::new("x"), 0);
        assert!(result.is_none());
    }

    /// Euler I: ∫ dx/√(x²+2x+3) — no Integral residue (plan acceptance).
    #[test]
    fn euler_i_inv_sqrt_quadratic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "1/sqrt(x^2+2*x+3)").unwrap();
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let s = result.to_string();
        assert!(!s.contains("Integral("), "got {s}");
        assert!(s.contains("asinh") || s.contains("log"), "got {s}");
    }

    /// ∫ √(x²+1) dx — solved (Euler or rule G8).
    #[test]
    fn euler_i_sqrt_x2_plus_1() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "sqrt(x^2+1)").unwrap();
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral("), "got {result}");
    }

    /// ∫ dx/(x·√(x²+x+1)) — Euler I with the sqrt inside a product.
    #[test]
    fn euler_i_inv_x_sqrt_quadratic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "1/(x*sqrt(x^2+x+1))").unwrap();
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral("), "got {result}");
    }

    /// Euler II: ∫ dx/√(2x²+3x+1) (a = 2 not a square, c = 1 a square).
    #[test]
    fn euler_ii_inv_sqrt_quadratic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "1/sqrt(2*x^2+3*x+1)").unwrap();
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral("), "got {result}");
    }
}
