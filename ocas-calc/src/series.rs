//! Symbolic Taylor series expansion for oCAS.
//!
//! This module provides [`taylor`], which expands an expression around a point
//! using repeated symbolic differentiation and evaluation at the expansion
//! point.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::simplify::simplify;
use ocas_rewrite::transformer::transform;

use crate::derivative::diff;
use crate::rules::calculus_rules;

/// Compute the Taylor expansion of `expr` around `point` with respect to `var`,
/// up to the given `order` (inclusive).
///
/// The result is the truncated polynomial
/// `sum_{n=0}^{order} f^(n)(point) * (var - point)^n / n!`.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::taylor;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.fun("exp", &[x]);
/// let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
/// assert_eq!(result.to_string(), "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))");
/// ```
pub fn taylor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    point: Atom<'a>,
    order: usize,
) -> Atom<'a> {
    let rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);

    let mut current = expr;
    let mut sum: Option<Atom<'a>> = None;
    let x_minus_p = ctx.add(&[ctx.var(var.as_str()), ctx.mul(&[ctx.num(-1), point])]);

    for n in 0..=order {
        let value_at_point = substitute(ctx, current, var, point);
        let coeff = mul_by_factorial_inverse(ctx, value_at_point, n);
        let term = if n == 0 {
            coeff
        } else {
            ctx.mul(&[coeff, ctx.pow(x_minus_p, ctx.num(n as i64))])
        };

        sum = Some(match sum {
            Some(prev) => ctx.add(&[prev, term]),
            None => term,
        });

        if n < order {
            current = diff(ctx, current, var);
        }
    }

    let raw = sum.expect("order >= 0 guarantees at least one term");
    let simplified = simplify(ctx, raw, &rules, 20);
    normalize(ctx, simplified)
}

/// Replace every occurrence of `var` inside `expr` with `replacement`.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::substitute;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), x]);
/// let result = substitute(&ctx, expr, Symbol::new("x"), y);
/// assert_eq!(result.to_string(), "(y^2) + y");
/// ```
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    transform(ctx, expr, |a| match a.node() {
        AtomNode::Var(v) if *v == var => Some(replacement),
        _ => None,
    })
}

fn mul_by_factorial_inverse<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, n: usize) -> Atom<'a> {
    if n == 0 {
        return expr;
    }
    let mut fact: i64 = 1;
    for i in 2..=n {
        fact = fact.checked_mul(i as i64).expect("factorial fits in i64");
    }
    ctx.mul(&[expr, ctx.pow(ctx.num(fact), ctx.num(-1))])
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn taylor_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[x]);
        let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
        assert_eq!(
            result.to_string(),
            "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))"
        );
    }

    #[test]
    fn taylor_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[x]);
        let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 5);
        assert_eq!(
            result.to_string(),
            "x + (-1*(6^-1)*(x^3)) + ((120^-1)*(x^5))"
        );
    }

    #[test]
    fn substitute_variable() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let expr = ctx.add(&[x, ctx.fun("sin", &[x])]);
        let result = substitute(&ctx, expr, Symbol::new("x"), y);
        assert_eq!(result.to_string(), "y + (sin(y))");
    }
}
