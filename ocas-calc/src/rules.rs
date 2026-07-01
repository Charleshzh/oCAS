//! Rewrite rules used by the calculus modules.
//!
//! These rules complement the default algebraic simplifier with a few
//! identities that are useful when simplifying derivatives, integrals, and
//! Taylor coefficients, such as `exp(0) -> 1` and `sin(0) -> 0`.

use ocas_atom::{AtomArena, Symbol};
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
}
