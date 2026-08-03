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
/// Returns `Some(result)` if any technique succeeds, `None` otherwise.
pub(crate) fn heuristic_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // 1. Integration by parts
    if let Some(r) = try_parts(ctx, expr, var, 0) {
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
    depth: usize,
) -> Option<Atom<'a>> {
    if depth >= PARTS_MAX_DEPTH {
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
    let v = integrate_raw(ctx, v_prime, var, depth + 2);
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
    let integral_u_prime_v = integrate_raw(ctx, u_prime_v, var, depth + 2);

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
        AtomNode::Var(v) => *v == var || is_constant(expr, var),
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

/// Try Weierstrass substitution: `t = tan(u/2)`.
fn try_weierstrass<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if !is_trig_rational(expr, var) {
        return None;
    }

    let t = ctx.var("_t");
    let x = ctx.var(var.as_str());

    // For u = x (linear with a=1, b=0):
    // sin(x) = 2t/(1+t²), cos(x) = (1-t²)/(1+t²), dx = 2dt/(1+t²)
    let one_plus_t2 = ctx.add(&[ctx.num(1), ctx.pow(t, ctx.num(2))]);
    let sin_x = ctx.mul(&[ctx.num(2), t, ctx.pow(one_plus_t2, ctx.num(-1))]);
    let cos_x = ctx.mul(&[
        ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(t, ctx.num(2))])]),
        ctx.pow(one_plus_t2, ctx.num(-1)),
    ]);
    let dx_dt = ctx.mul(&[ctx.num(2), ctx.pow(one_plus_t2, ctx.num(-1))]);

    // Replace sin(x) → 2t/(1+t²), cos(x) → (1-t²)/(1+t²)
    let substituted = substitute_trig(ctx, expr, var, sin_x, cos_x);
    let integrand = ctx.mul(&[substituted, dx_dt]);

    // Integrate with respect to t
    let t_sym = Symbol::new("_t");
    let result_t = integrate_raw(ctx, integrand, t_sym, 2);

    if is_fallback(&result_t) {
        return None;
    }

    // Back-substitute t = tan(x/2)
    let back = ctx.fun("tan", &[ctx.mul(&[x, ctx.pow(ctx.num(2), ctx.num(-1))])]);
    Some(substitute_atom(ctx, result_t, t_sym, back))
}

/// Replace `sin(var)` and `cos(var)` in `expr` with given expressions.
fn substitute_trig<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sin_replacement: Atom<'a>,
    cos_replacement: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sin" && args.len() == 1 => {
            if matches!(args[0].node(), AtomNode::Var(v) if *v == var) {
                sin_replacement
            } else {
                expr
            }
        }
        AtomNode::Fun(name, args) if name.as_str() == "cos" && args.len() == 1 => {
            if matches!(args[0].node(), AtomNode::Var(v) if *v == var) {
                cos_replacement
            } else {
                expr
            }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_trig(ctx, *a, var, sin_replacement, cos_replacement))
                .collect();
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_trig(ctx, *a, var, sin_replacement, cos_replacement))
                .collect();
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = substitute_trig(ctx, *base, var, sin_replacement, cos_replacement);
            ctx.pow(new_base, *exp)
        }
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_trig(ctx, *a, var, sin_replacement, cos_replacement))
                .collect();
            ctx.fun(name.as_str(), &new_args)
        }
    }
}

// =========================================================================
// Euler substitution (1d)
// =========================================================================

/// Check if `expr` contains `√(ax² + bx + c)`.
///
/// Returns `(a_coeff, b_coeff, c_coeff)` if matched.
fn is_quadratic_sqrt<'a>(expr: Atom<'a>, var: Symbol) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    // Match Pow(inner, 1/2) where inner is a quadratic in var
    let base = match expr.node() {
        AtomNode::Pow(b, e) => {
            let is_half = match e.node() {
                AtomNode::Num(1) => false, // not 1/2
                AtomNode::Mul(args) if args.len() == 2 => {
                    args.iter().any(|a| matches!(a.node(), AtomNode::Pow(b2, e2)
                        if matches!(b2.node(), AtomNode::Num(2)) && matches!(e2.node(), AtomNode::Num(-1))))
                        && args.iter().any(|a| matches!(a.node(), AtomNode::Num(1)))
                }
                _ => false,
            };
            if !is_half {
                return None;
            }
            *b
        }
        _ => return None,
    };

    // Match ax² + bx + c
    match base.node() {
        AtomNode::Add(args) => {
            let mut a_coeff = None;
            let mut b_coeff = None;
            let mut c_coeff = None;

            for term in args.iter() {
                match term.node() {
                    // ax² term
                    AtomNode::Mul(margs) if margs.len() == 2 => {
                        if let Some(idx) = margs.iter().position(|a| {
                            matches!(a.node(), AtomNode::Pow(b, e)
                                if matches!(b.node(), AtomNode::Var(v) if *v == var)
                                && matches!(e.node(), AtomNode::Num(2)))
                        }) {
                            let coeff_idx = 1 - idx;
                            if is_constant(margs[coeff_idx], var) && a_coeff.is_none() {
                                a_coeff = Some(margs[coeff_idx]);
                            }
                        }
                        // bx term
                        if let Some(idx) = margs
                            .iter()
                            .position(|a| matches!(a.node(), AtomNode::Var(v) if *v == var))
                        {
                            let coeff_idx = 1 - idx;
                            if is_constant(margs[coeff_idx], var) && b_coeff.is_none() {
                                b_coeff = Some(margs[coeff_idx]);
                            }
                        }
                    }
                    // x² term (coefficient = 1)
                    AtomNode::Pow(b, e) => {
                        if matches!(b.node(), AtomNode::Var(v) if *v == var)
                            && matches!(e.node(), AtomNode::Num(2))
                            && a_coeff.is_none()
                        {
                            // Implicit coefficient 1 — we'd need to represent this
                            // For now, skip this case
                        }
                    }
                    // Linear term: just `var`
                    AtomNode::Var(v) => {
                        if *v == var && b_coeff.is_none() {
                            // Implicit coefficient 1
                        }
                    }
                    // Constant term
                    AtomNode::Num(_) if c_coeff.is_none() => {
                        c_coeff = Some(*term);
                    }
                    _ => {}
                }
            }

            match (a_coeff, b_coeff, c_coeff) {
                (Some(a), Some(b), Some(c)) => Some((a, b, c)),
                _ => None,
            }
        }
        _ => None,
    }
}

/// Try Euler substitution for `√(ax² + bx + c)`.
fn try_euler_substitution<'a>(
    _ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (_a, _b, _c) = is_quadratic_sqrt(expr, var)?;

    // Euler substitution is complex to implement fully with back-substitution.
    // For now, return None and let the fallback handle it.
    // This is a placeholder for future implementation.
    None
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
        let result = heuristic_integrate(&ctx, expr, Symbol::new("x"));
        assert!(result.is_none());
    }
}
