//! Rewrite rules used by the calculus modules.
//!
//! These rules complement the default algebraic simplifier with a few
//! identities that are useful when simplifying derivatives, integrals, and
//! Taylor coefficients, such as `exp(0) -> 1` and `sin(0) -> 0`.

use ocas_atom::{AtomArena, Symbol};
use ocas_rewrite::matcher::MatchValue;
use ocas_rewrite::pattern::{Pattern, PatternAlloc, WildcardLevel};
use ocas_rewrite::rules::Rule;

/// Build the calculus-specific rule set used to clean up results from
/// [`crate::derivative::diff`], [`crate::integral::integrate`], and
/// [`crate::series::taylor`].
pub fn calculus_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>> {
    let mut rules = ocas_rewrite::rules::default_rules(ctx, alloc);

    // exp(0) -> 1
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("exp"),
            alloc.alloc_slice(&[Pattern::Literal(ctx.num(0))]).to_vec(),
        ),
        |_bindings, ctx| ctx.num(1),
    ));

    // log(1) -> 0
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("log"),
            alloc.alloc_slice(&[Pattern::Literal(ctx.num(1))]).to_vec(),
        ),
        |_bindings, ctx| ctx.num(0),
    ));

    // sin(0) -> 0
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("sin"),
            alloc.alloc_slice(&[Pattern::Literal(ctx.num(0))]).to_vec(),
        ),
        |_bindings, ctx| ctx.num(0),
    ));

    // cos(0) -> 1
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("cos"),
            alloc.alloc_slice(&[Pattern::Literal(ctx.num(0))]).to_vec(),
        ),
        |_bindings, ctx| ctx.num(1),
    ));

    // tan(0) -> 0
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("tan"),
            alloc.alloc_slice(&[Pattern::Literal(ctx.num(0))]).to_vec(),
        ),
        |_bindings, ctx| ctx.num(0),
    ));

    // polylog(1, z) -> -log(1 - z): Li₁ is the elementary logarithm, so the
    // order-1 head never needs to survive simplification (0.31.0).
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("polylog"),
            alloc
                .alloc_slice(&[
                    Pattern::Literal(ctx.num(1)),
                    Pattern::Wildcard(Symbol::new("z"), WildcardLevel::Single),
                ])
                .to_vec(),
        ),
        |bindings, ctx| {
            let bound = bindings
                .get(Symbol::new("z"))
                .expect("polylog(1, z) binds z");
            let MatchValue::Single(z) = bound else {
                panic!("the z wildcard is a single-atom binding");
            };
            let one_minus_z = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), *z])]);
            ctx.mul(&[ctx.num(-1), ctx.fun("log", &[one_minus_z])])
        },
    ));

    // polylog(2, 0) -> 0 and polylog(2, 1) -> pi^2/6: the two constant values
    // the reduction and the tests rely on (0.31.0). Numeric evaluation for
    // general arguments stays in the integration oracle, not in the rules.
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("polylog"),
            alloc
                .alloc_slice(&[Pattern::Literal(ctx.num(2)), Pattern::Literal(ctx.num(0))])
                .to_vec(),
        ),
        |_bindings, ctx| ctx.num(0),
    ));
    rules.push(Rule::new(
        Pattern::Fun(
            Symbol::new("polylog"),
            alloc
                .alloc_slice(&[Pattern::Literal(ctx.num(2)), Pattern::Literal(ctx.num(1))])
                .to_vec(),
        ),
        |_bindings, ctx| {
            ctx.mul(&[
                ctx.pow(ctx.var("pi"), ctx.num(2)),
                ctx.pow(ctx.num(6), ctx.num(-1)),
            ])
        },
    ));

    // (-1) * (-1) -> 1
    rules.push(Rule::new(
        Pattern::Mul(
            alloc
                .alloc_slice(&[Pattern::Literal(ctx.num(-1)), Pattern::Literal(ctx.num(-1))])
                .to_vec(),
        ),
        |_bindings, ctx| ctx.num(1),
    ));

    // 1 ^ _ -> 1
    rules.push(Rule::new(
        Pattern::Pow(Box::new((
            Pattern::Literal(ctx.num(1)),
            Pattern::Wildcard(Symbol::new("n"), WildcardLevel::Single),
        ))),
        |_bindings, ctx| ctx.num(1),
    ));

    rules
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_rewrite::simplify::simplify;

    use super::*;
    use crate::pattern_alloc::VecAlloc;

    #[test]
    fn exp_zero_simplifies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = calculus_rules(&ctx, &VecAlloc);
        let expr = ctx.fun("exp", &[ctx.num(0)]);
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "1");
    }

    #[test]
    fn sin_zero_simplifies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = calculus_rules(&ctx, &VecAlloc);
        let expr = ctx.fun("sin", &[ctx.num(0)]);
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn polylog_simplifications() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = calculus_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        // Li₁(z) = −log(1−z): the order-1 head is elementary.
        let expr = ctx.fun("polylog", &[ctx.num(1), x]);
        let got = simplify(&ctx, expr, &rules, 10).to_string();
        assert!(
            got.contains("log(1 + (-1*x))") && !got.contains("polylog"),
            "polylog(1, x) simplified to {got:?}"
        );
        // The two constant values.
        let expr = ctx.fun("polylog", &[ctx.num(2), ctx.num(0)]);
        assert_eq!(simplify(&ctx, expr, &rules, 10).to_string(), "0");
        let expr = ctx.fun("polylog", &[ctx.num(2), ctx.num(1)]);
        let got = simplify(&ctx, expr, &rules, 10).to_string();
        assert!(
            got.contains("pi") && got.contains("6"),
            "polylog(2, 1) simplified to {got:?}"
        );
        // A generic order-2 head is left untouched (no over-rewriting).
        let expr = ctx.fun("polylog", &[ctx.num(2), x]);
        assert_eq!(
            simplify(&ctx, expr, &rules, 10).to_string(),
            "polylog(2, x)"
        );
    }
}
