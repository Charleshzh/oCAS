//! Pattern matching and rewriting for oCAS.
//!
//! This crate builds on [`ocas_atom::Atom`] and provides:
//!
//! - A [`Pattern`] type with wildcards of three levels.
//! - A [`match_pattern`] engine that binds
//!   wildcards to sub-expressions.
//! - A [`transform`] function for bottom-up
//!   tree rewriting.
//! - A [`Rule`] abstraction and a [`simplify()`](crate::simplify::simplify)
//!   entry point.

pub mod matcher;
pub mod pattern;
pub mod rules;
pub mod simplify;
pub mod transformer;

pub use matcher::{Bindings, MatchError, match_pattern};
pub use pattern::{Pattern, WildcardLevel};
pub use rules::Rule;
pub use simplify::simplify;
pub use transformer::transform;

#[cfg(test)]
mod proptests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use proptest::prelude::*;
    use rules::default_rules;
    use simplify::simplify;

    #[derive(Debug, Clone)]
    enum PropExpr {
        Num(i64),
        Var(&'static str),
        Add(Vec<PropExpr>),
        Mul(Vec<PropExpr>),
        Pow(Box<PropExpr>, Box<PropExpr>),
    }

    fn build<'a>(ctx: &'a AtomArena<'a>, expr: &'a PropExpr) -> ocas_atom::Atom<'a> {
        match expr {
            PropExpr::Num(n) => ctx.num(*n),
            PropExpr::Var(name) => ctx.var(name),
            PropExpr::Add(args) => {
                let atoms: Vec<_> = args.iter().map(|a| build(ctx, a)).collect();
                ctx.add(&atoms)
            }
            PropExpr::Mul(args) => {
                let atoms: Vec<_> = args.iter().map(|a| build(ctx, a)).collect();
                ctx.mul(&atoms)
            }
            PropExpr::Pow(base, exp) => ctx.pow(build(ctx, base), build(ctx, exp)),
        }
    }

    fn prop_expr() -> impl Strategy<Value = PropExpr> {
        let leaf = prop_oneof! {
            (-100..100i64).prop_map(PropExpr::Num),
            Just(PropExpr::Var("x")),
            Just(PropExpr::Var("y")),
        };
        leaf.prop_recursive(4, 64, 4, |inner| {
            prop_oneof! {
                prop::collection::vec(inner.clone(), 1..4).prop_map(PropExpr::Add),
                prop::collection::vec(inner.clone(), 1..4).prop_map(PropExpr::Mul),
                (inner.clone(), inner.clone())
                    .prop_map(|(b, e)| PropExpr::Pow(Box::new(b), Box::new(e))),
            }
        })
    }

    fn simplify_expr<'a>(ctx: &'a AtomArena<'a>, expr: ocas_atom::Atom<'a>) -> ocas_atom::Atom<'a> {
        let rules = default_rules(ctx, &());
        simplify(ctx, expr, &rules, 20)
    }

    proptest! {
        #[test]
        fn simplify_is_idempotent(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build(&ctx, &expr);
            let once = simplify_expr(&ctx, atom);
            let twice = simplify_expr(&ctx, once);
            assert_eq!(once.to_string(), twice.to_string());
        }

        #[test]
        fn simplify_adds_zero(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build(&ctx, &expr);
            let with_zero = ctx.add(&[atom, ctx.num(0)]);
            let simplified = simplify_expr(&ctx, with_zero);
            let expected = simplify_expr(&ctx, atom);
            assert_eq!(simplified.to_string(), expected.to_string());
        }

        #[test]
        fn simplify_multiplies_zero(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build(&ctx, &expr);
            let with_zero = ctx.mul(&[atom, ctx.num(0)]);
            let simplified = simplify_expr(&ctx, with_zero);
            assert_eq!(simplified.to_string(), "0");
        }
    }
}
