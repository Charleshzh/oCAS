//! Calculus primitives for oCAS.
//!
//! This crate provides symbolic differentiation, integration, and series
//! expansion for [`Atom`](ocas_atom::Atom) expression trees. Results are
//! returned as new atoms and are automatically simplified using the
//! rewrite engine.
//!
//! Unresolved or partially-resolved derivatives and integrals are represented
//! as the reserved function forms `Derivative(expr, var)` and
//! `Integral(expr, var)`.

#![warn(missing_docs)]

pub mod derivative;
pub mod integral;
pub mod series;

mod pattern_alloc;
mod rules;

pub use derivative::diff;
pub use integral::integrate;
pub use series::{substitute, taylor};

#[cfg(test)]
mod proptests {
    use ocas_atom::AtomArena;
    use ocas_atom::Symbol;
    use ocas_core::arena::Arena;
    use ocas_rewrite::rules::default_rules;
    use ocas_rewrite::simplify::simplify;
    use proptest::prelude::*;

    use super::*;

    fn simplify_atom<'a>(ctx: &'a AtomArena<'a>, expr: ocas_atom::Atom<'a>) -> ocas_atom::Atom<'a> {
        let rules = default_rules(ctx, &());
        simplify(ctx, expr, &rules, 20)
    }

    /// Generate a simple polynomial in `x` with small integer coefficients.
    fn poly_expr() -> impl Strategy<Value = Vec<(i64, usize)>> {
        prop::collection::vec((-5..5i64, 0..6usize), 1..5)
    }

    proptest! {
        #[test]
        fn diff_of_power(terms in poly_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x_sym = Symbol::new("x");
            let x = ctx.var("x");

            // Restrict to a single power term with positive coefficient and exponent >= 1.
            let (_, power) = terms[0];
            prop_assume!(power >= 1);
            let n = power as i64;
            let expr = ctx.pow(x, ctx.num(n));

            let d = diff(&ctx, expr, x_sym);
            let expected = ctx.mul(&[ctx.num(n), ctx.pow(x, ctx.num(n - 1))]);
            assert_eq!(
                simplify_atom(&ctx, d).to_string(),
                simplify_atom(&ctx, expected).to_string()
            );
        }

        #[test]
        fn integrate_constant(terms in poly_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x_sym = Symbol::new("x");
            let x = ctx.var("x");

            let (c, _) = terms[0];
            let expr = ctx.num(c);
            let integrated = integrate(&ctx, expr, x_sym);
            let expected = ctx.mul(&[ctx.num(c), x]);
            assert_eq!(
                simplify_atom(&ctx, integrated).to_string(),
                simplify_atom(&ctx, expected).to_string()
            );
        }
    }
}
