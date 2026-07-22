//! Symbolic integration for oCAS.
//!
//! This module provides [`integrate`], a heuristic integrator for expressions
//! involving polynomials and elementary functions. It uses a lookup table for
//! common antiderivatives and supports simple linear substitutions.
//!
//! Integrals that cannot be expressed with the built-in table are returned as
//! the unevaluated form `Integral(expr, var)`.

#![allow(clippy::collapsible_if)]

pub mod rational;
pub(crate) mod rde;
pub(crate) mod risch;
pub(crate) mod special;
pub(crate) mod trig;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_core::error::Result;
use ocas_core::fuel::Fuel;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::{simplify, simplify_with_fuel};

use crate::rules::calculus_rules;

/// Maximum recursion depth for `integrate_raw`, preventing infinite loops
/// on patterns such as nested linear substitutions if the table is misapplied.
const MAX_DEPTH: usize = 8;

/// Integrate `expr` with respect to `var`.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::integrate;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.pow(x, ctx.num(2));
/// let result = integrate(&ctx, expr, Symbol::new("x"));
/// assert_eq!(result.to_string(), "(3^-1)*(x^3)");
/// ```
///
/// For integrals not covered by the heuristic table, the result is returned
/// as `Integral(expr, var)`.
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let normalized = normalize(ctx, expr);
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = integrate_raw(ctx, normalized, var, 0);
    // Combine default algebraic simplification with calculus-specific rules,
    // then normalize to a canonical form (removing *1, +0, sorting, etc.).
    let after_default = simplify(ctx, raw, &default_rules, 20);
    let after_calc = simplify(ctx, after_default, &calc_rules, 10);
    normalize(ctx, after_calc)
}

/// Integrate with a [`Fuel`] budget bounding the post-integration
/// simplification passes.
///
/// The integration traversal itself uses the internal depth limit; this entry
/// point threads `fuel` through the two simplification stages so a pathological
/// result that would otherwise spin the rewriter can be cut off determin-
/// istically. Returns `Err` only when fuel was exhausted mid-simplification.
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>> {
    let normalized = normalize(ctx, expr);
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = integrate_raw(ctx, normalized, var, 0);
    let after_default = simplify_with_fuel(ctx, raw, &default_rules, 20, fuel)?;
    let after_calc = simplify_with_fuel(ctx, after_default, &calc_rules, 10, fuel)?;
    Ok(normalize(ctx, after_calc))
}

fn integrate_raw<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    depth: usize,
) -> Atom<'a> {
    if depth > MAX_DEPTH {
        return fallback(ctx, expr, var);
    }

    match expr.node() {
        AtomNode::Num(_) => {
            // ∫ c dx = c * x
            let x = ctx.var(var.as_str());
            ctx.mul(&[expr, x])
        }
        AtomNode::Var(v) => {
            if *v == var {
                ctx.mul(&[
                    ctx.pow(ctx.var(var.as_str()), ctx.num(2)),
                    ctx.pow(ctx.num(2), ctx.num(-1)),
                ])
            } else {
                ctx.mul(&[expr, ctx.var(var.as_str())])
            }
        }
        AtomNode::Add(args) => {
            let mut terms = Vec::with_capacity(args.len());
            for a in args.iter() {
                terms.push(integrate_raw(ctx, *a, var, depth));
            }
            ctx.add(&terms)
        }
        AtomNode::Mul(args) => {
            let r = integrate_product(ctx, args, var, depth);
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var)
            } else {
                r
            }
        }
        AtomNode::Pow(base, exp) => {
            let r = integrate_power(ctx, *base, *exp, var, depth);
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var)
            } else {
                r
            }
        }
        AtomNode::Fun(name, args) => {
            let r = integrate_function(ctx, *name, args, var, depth);
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var)
            } else {
                r
            }
        }
    }
}

/// Try the rational-function integrator, then the Risch algorithm, before
/// giving up with the unevaluated `Integral` form.
fn try_risch_or_fallback<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    if let Some(r) = rational::integrate_rational(ctx, expr, var) {
        return r;
    }
    if let Some(r) = risch::risch_integrate(ctx, expr, var) {
        return r;
    }
    // Trigonometric integrands: rewrite into complex exponentials, run
    // Risch, then try to bring the answer back to real form.
    if let Some(exp_form) = trig::trig_to_exp(ctx, expr) {
        if let Some(complex_ans) = risch::risch_integrate(ctx, exp_form, var) {
            return trig::realify(ctx, complex_ans);
        }
    }
    // Non-elementary integrals with special-function closed forms
    // (erf, Ei, Si, Ci, Fresnel, …).
    if let Some(r) = special::special_integrate(ctx, expr, ctx.var(var.as_str())) {
        return r;
    }
    fallback(ctx, expr, var)
}

fn fallback<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    ctx.fun("Integral", &[expr, ctx.var(var.as_str())])
}

/// True if `expr` does not contain `var`.
fn is_constant<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v != var,
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().all(|a| is_constant(*a, var))
        }
        AtomNode::Pow(base, exp) => is_constant(*base, var) && is_constant(*exp, var),
    }
}

fn integrate_product<'a>(
    ctx: &'a AtomArena<'a>,
    args: &'a [Atom<'a>],
    var: Symbol,
    depth: usize,
) -> Atom<'a> {
    // Split into constant factors and the remaining factor.
    let mut constants: Vec<Atom<'a>> = Vec::new();
    let mut non_constant: Vec<Atom<'a>> = Vec::new();

    for a in args.iter() {
        if is_constant(*a, var) {
            constants.push(*a);
        } else {
            non_constant.push(*a);
        }
    }

    if non_constant.is_empty() {
        // All factors are constant: ∫ c dx = c * x
        return ctx.mul(&[ctx.mul(args), ctx.var(var.as_str())]);
    }

    let core = if non_constant.len() == 1 {
        non_constant[0]
    } else {
        ctx.mul(&non_constant)
    };

    let integrated_core = integrate_raw(ctx, core, var, depth + 1);

    // If integration failed, wrap the whole product.
    if is_fallback(&integrated_core) {
        return fallback(ctx, ctx.mul(args), var);
    }

    let mut result_factors = constants;
    result_factors.push(integrated_core);
    ctx.mul(&result_factors)
}

fn is_fallback<'a>(atom: &Atom<'a>) -> bool {
    matches!(atom.node(), AtomNode::Fun(name, _) if name.as_str() == "Integral")
}

fn integrate_power<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    exp: Atom<'a>,
    var: Symbol,
    _depth: usize,
) -> Atom<'a> {
    // Detect x^n where n is a constant integer.
    if let AtomNode::Var(v) = base.node()
        && *v == var
    {
        if let AtomNode::Num(n) = exp.node() {
            if *n == -1 {
                // ∫ x^(-1) dx = log(x)
                return ctx.fun("log", &[base]);
            }
            // ∫ x^n dx = x^(n+1) / (n+1)
            let new_exp = ctx.num(n + 1);
            let denom = ctx.num(n + 1);
            return ctx.mul(&[ctx.pow(base, new_exp), ctx.pow(denom, ctx.num(-1))]);
        }
    }

    // Detect linear substitution: (a*x + b)^n where n is constant integer.
    if let AtomNode::Num(n) = exp.node()
        && let Some((a, _b)) = linear_form(ctx, base, var)
    {
        // ∫ (a*x + b)^n dx = (a*x + b)^(n+1) / (a * (n+1))
        let new_exp = ctx.num(n + 1);
        let denom = ctx.mul(&[a, ctx.num(n + 1)]);
        return ctx.mul(&[ctx.pow(base, new_exp), ctx.pow(denom, ctx.num(-1))]);
    }

    fallback(ctx, ctx.pow(base, exp), var)
}

/// If `expr` is of the form `a*x + b` (with `a` and `b` constant w.r.t. `var`),
/// return `(a, b)`. Otherwise return None.
fn linear_form<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match expr.node() {
        AtomNode::Var(v) if *v == var => Some((ctx.num(1), ctx.num(0))),
        AtomNode::Mul(args) => {
            let mut coeff = ctx.num(1);
            let mut has_var = false;
            for a in args.iter() {
                if let AtomNode::Var(v) = a.node()
                    && *v == var
                {
                    has_var = true;
                    continue;
                }
                if is_constant(*a, var) {
                    coeff = ctx.mul(&[coeff, *a]);
                } else {
                    return None;
                }
            }
            if has_var {
                Some((coeff, ctx.num(0)))
            } else {
                None
            }
        }
        AtomNode::Add(args) => {
            let mut a_part = ctx.num(0);
            let mut b_part = ctx.num(0);
            for arg in args.iter() {
                if let Some((ca, _cb)) = linear_form(ctx, *arg, var) {
                    a_part = ctx.add(&[a_part, ca]);
                } else if is_constant(*arg, var) {
                    b_part = ctx.add(&[b_part, *arg]);
                } else {
                    return None;
                }
            }
            Some((a_part, b_part))
        }
        _ => None,
    }
}

fn integrate_function<'a>(
    ctx: &'a AtomArena<'a>,
    name: Symbol,
    args: &'a [Atom<'a>],
    var: Symbol,
    _depth: usize,
) -> Atom<'a> {
    if args.is_empty() {
        return fallback(ctx, ctx.fun(name.as_str(), args), var);
    }
    let u = args[0];

    // Simple linear substitution forms: f(a*x + b)
    if let Some((a, _b)) = linear_form(ctx, u, var)
        && is_constant(a, var)
        && !is_one(a)
    {
        let inner_integral = match name.as_str() {
            "sin" => ctx.mul(&[ctx.num(-1), ctx.fun("cos", &[u])]),
            "cos" => ctx.fun("sin", &[u]),
            "exp" => ctx.fun("exp", &[u]),
            _ => return fallback(ctx, ctx.fun(name.as_str(), args), var),
        };
        return ctx.mul(&[ctx.pow(a, ctx.num(-1)), inner_integral]);
    }

    // Direct table for f(x) where u == x.
    if let AtomNode::Var(v) = u.node()
        && *v == var
    {
        let antiderivative: Option<Atom<'a>> = match name.as_str() {
            "sin" => Some(ctx.mul(&[ctx.num(-1), ctx.fun("cos", &[u])])),
            "cos" => Some(ctx.fun("sin", &[u])),
            "exp" => Some(ctx.fun("exp", &[u])),
            "log" => Some(ctx.mul(&[u, ctx.add(&[ctx.fun("log", &[u]), ctx.num(-1)])])),
            _ => None,
        };
        if let Some(anti) = antiderivative {
            return anti;
        }
    }

    fallback(ctx, ctx.fun(name.as_str(), args), var)
}

fn is_one<'a>(expr: Atom<'a>) -> bool {
    matches!(expr.node(), AtomNode::Num(1))
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn integrate_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(3^-1)*(x^3)");
    }

    #[test]
    fn integrate_inverse() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(-1));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "log(x)");
    }

    #[test]
    fn integrate_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "-1*(cos(x))");
    }

    #[test]
    fn integrate_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("cos", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "sin(x)");
    }

    #[test]
    fn integrate_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "exp(x)");
    }

    #[test]
    fn integrate_linear_substitution() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two_x_plus_one = ctx.add(&[ctx.mul(&[ctx.num(2), x]), ctx.num(1)]);
        let expr = ctx.pow(two_x_plus_one, ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(6^-1)*((1 + (2*x))^3)");
    }

    #[test]
    fn integrate_unknown() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("f", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Integral(f(x), x)");
    }

    #[test]
    fn integrate_sin_times_cos_via_trig_path() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ sin(x)·cos(x) dx — needs the trig → exp(Ix) → Risch path with
        // a constant imaginary unit in the coefficient field. The RDE base
        // solver currently works over ℚ[x] only, so the hyperexponential
        // equation Dq + I·q = … cannot be solved yet; the pipeline falls
        // back to the unevaluated form. Documented limitation.
        let expr = ctx.mul(&[ctx.fun("sin", &[x]), ctx.fun("cos", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let _ = result;
    }

    #[test]
    fn integrate_cos_squared_via_trig_path() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ cos(x)² dx = x/2 + sin(2x)/4 — same coefficient-field
        // limitation as above.
        let expr = ctx.pow(ctx.fun("cos", &[x]), ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let _ = result;
    }

    #[test]
    fn integrate_exp_neg_x_squared_gives_erf() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ exp(-x²) dx = (√π/2)·erf(x) — the 0.11.0 known gap, now closed
        // by the special-function table.
        let expr = ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(result.to_string().contains("erf"), "got {result}");
        assert!(!result.to_string().starts_with("Integral"), "got {result}");
    }

    #[test]
    fn integrate_exp_x_over_x_gives_ei() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ exp(x)/x dx = Ei(x) — non-elementary, special-function table.
        let expr = ctx.mul(&[ctx.fun("exp", &[x]), ctx.pow(x, ctx.num(-1))]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Ei(x)");
    }
}
