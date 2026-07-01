//! Symbolic differentiation for oCAS.
//!
//! This module provides [`diff`], which computes the symbolic derivative of an
//! expression with respect to a variable. Elementary function rules are
//! hard-coded in a small table, and the chain rule is applied automatically for
//! compound arguments.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::simplify::simplify;

use crate::rules::calculus_rules;

/// Differentiate `expr` with respect to `var`.
///
/// The function implements the standard recursive derivative rules for sums,
/// products, powers, and a built-in table of elementary functions. The result
/// is simplified using the default rewrite rules plus a few calculus-specific
/// identities.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::diff;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let sin_x = ctx.fun("sin", &[x]);
/// let result = diff(&ctx, sin_x, Symbol::new("x"));
/// assert_eq!(result.to_string(), "cos(x)");
/// ```
///
/// For unknown functions, the derivative is returned as an unevaluated
/// `Derivative(f, x)` form.
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = diff_raw(ctx, expr, var);
    let simplified = simplify(ctx, raw, &rules, 20);
    normalize(ctx, simplified)
}

fn diff_raw<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    match expr.node() {
        AtomNode::Num(_) => ctx.num(0),
        AtomNode::Var(v) => {
            if *v == var {
                ctx.num(1)
            } else {
                ctx.num(0)
            }
        }
        AtomNode::Add(args) => {
            let mut terms = Vec::with_capacity(args.len());
            for a in args.iter() {
                terms.push(diff_raw(ctx, *a, var));
            }
            ctx.add(&terms)
        }
        AtomNode::Mul(args) => {
            // Product rule: d/dx (a * b * c) = a' * b * c + a * b' * c + ...
            let mut sum_terms = Vec::with_capacity(args.len());
            for i in 0..args.len() {
                let mut factors = Vec::with_capacity(args.len());
                for (j, a) in args.iter().enumerate() {
                    if i == j {
                        factors.push(diff_raw(ctx, *a, var));
                    } else {
                        factors.push(*a);
                    }
                }
                sum_terms.push(ctx.mul(&factors));
            }
            ctx.add(&sum_terms)
        }
        AtomNode::Pow(base, exp) => {
            let base = *base;
            let exp = *exp;
            let d_base = diff_raw(ctx, base, var);
            let d_exp = diff_raw(ctx, exp, var);

            let is_exp_const = matches!(exp.node(), AtomNode::Num(_));
            let is_base_const = matches!(base.node(), AtomNode::Num(_));

            if is_exp_const {
                // d/dx (b^n) = n * b^(n-1) * b'
                let n = if let AtomNode::Num(n) = exp.node() {
                    *n
                } else {
                    0
                };
                ctx.mul(&[ctx.num(n), ctx.pow(base, ctx.num(n - 1)), d_base])
            } else if is_base_const {
                // d/dx (a^u) = a^u * log(a) * u'
                ctx.mul(&[ctx.pow(base, exp), ctx.fun("log", &[base]), d_exp])
            } else {
                // Generalized power rule:
                //   d/dx (b^e) = b^e * (log(b) * e' + e * b' / b)
                let log_b = ctx.fun("log", &[base]);
                let term1 = ctx.mul(&[log_b, d_exp]);
                let term2 = ctx.mul(&[exp, d_base, ctx.pow(base, ctx.num(-1))]);
                let factor = ctx.add(&[term1, term2]);
                ctx.mul(&[ctx.pow(base, exp), factor])
            }
        }
        AtomNode::Fun(name, args) => diff_function(ctx, *name, args, var),
    }
}

fn diff_function<'a>(
    ctx: &'a AtomArena<'a>,
    name: Symbol,
    args: &'a [Atom<'a>],
    var: Symbol,
) -> Atom<'a> {
    debug_assert!(
        !args.is_empty(),
        "diff_function should only be called with non-empty function arguments"
    );
    let u = args[0];
    let du = diff_raw(ctx, u, var);

    let name_str = name.as_str();
    let derivative_of_arg: Atom<'a> = match name_str {
        "sin" => ctx.fun("cos", &[u]),
        "cos" => ctx.mul(&[ctx.num(-1), ctx.fun("sin", &[u])]),
        "exp" => ctx.fun("exp", &[u]),
        "log" => ctx.pow(u, ctx.num(-1)),
        "sqrt" => ctx.pow(ctx.mul(&[ctx.num(2), ctx.fun("sqrt", &[u])]), ctx.num(-1)),
        "tan" => ctx.pow(ctx.fun("sec", &[u]), ctx.num(2)),
        "sec" => ctx.mul(&[ctx.fun("sec", &[u]), ctx.fun("tan", &[u])]),
        _ => {
            // Unknown function: return an unevaluated Derivative form.
            return ctx.fun(
                "Derivative",
                &[ctx.fun(name_str, args), ctx.var(var.as_str())],
            );
        }
    };

    ctx.mul(&[derivative_of_arg, du])
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn diff_number() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let n = ctx.num(7);
        let result = diff(&ctx, n, Symbol::new("x"));
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn diff_variable_same() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let result = diff(&ctx, x, Symbol::new("x"));
        assert_eq!(result.to_string(), "1");
    }

    #[test]
    fn diff_variable_other() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let y = ctx.var("y");
        let result = diff(&ctx, y, Symbol::new("x"));
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn diff_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(2));
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "2*x");
    }

    #[test]
    fn diff_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "cos(x)");
    }

    #[test]
    fn diff_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("cos", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "-1*(sin(x))");
    }

    #[test]
    fn diff_sqrt() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sqrt", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(2*(sqrt(x)))^-1");
    }

    #[test]
    fn diff_exp_squared() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let x2 = ctx.pow(x, ctx.num(2));
        let expr = ctx.fun("exp", &[x2]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "2*x*(exp(x^2))");
    }

    #[test]
    fn diff_product() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sin_x = ctx.fun("sin", &[x]);
        let expr = ctx.mul(&[x, sin_x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(sin(x)) + (x*(cos(x)))");
    }

    #[test]
    fn diff_unknown_function() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let expr = ctx.fun("f", &[x, y]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Derivative(f(x, y), x)");
    }
}
