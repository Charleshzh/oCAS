//! Rule-based simplification for oCAS.
//!
//! The [`simplify`] function applies a list of rewrite rules to an atom tree in
//! a bottom-up traversal, repeatedly, until no rule fires or a limit is
//! reached. It composes with the [`Transformer`](crate::transformer) and
//! [`Rule`] machinery.

use ocas_atom::{Atom, AtomArena};
use ocas_core::error::Result;
use ocas_core::fuel::Fuel;

use crate::rules::Rule;
use crate::transformer::transform;

/// Simplify an atom using the supplied rewrite rules.
///
/// Rules are applied bottom-up and repeatedly until a fixpoint is reached or
/// `iter_limit` iterations have been performed. For the default rule set, a
/// small limit (e.g., 20) is usually sufficient.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::{default_rules, Rule};
/// use ocas_rewrite::simplify::simplify;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let rules = default_rules(&ctx, &VecAlloc);
///
/// let x = ctx.var("x");
/// let expr = ctx.mul(&[x, ctx.num(0)]);
/// let result = simplify(&ctx, expr, &rules, 10);
///
/// assert_eq!(result.to_string(), "0");
/// # struct VecAlloc;
/// # impl<'a> ocas_rewrite::pattern::PatternAlloc<'a> for VecAlloc {
/// #     fn alloc_slice(&self, _: &[ocas_rewrite::pattern::Pattern<'a>]) -> &'a [ocas_rewrite::pattern::Pattern<'a>] { unreachable!() }
/// # }
/// ```
pub fn simplify<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    rules: &[Rule<'a>],
    iter_limit: usize,
) -> Atom<'a> {
    let mut current = transform(ctx, atom, |a| apply_rules(ctx, a, rules));
    for _ in 1..iter_limit {
        let next = transform(ctx, current, |a| apply_rules(ctx, a, rules));
        if next == current {
            return next;
        }
        current = next;
    }
    current
}

fn apply_rules<'a>(ctx: &'a AtomArena<'a>, atom: Atom<'a>, rules: &[Rule<'a>]) -> Option<Atom<'a>> {
    for rule in rules {
        if let Some(replacement) = rule.apply(ctx, atom) {
            return Some(replacement);
        }
    }
    None
}

/// Simplify with a [`Fuel`] budget.
///
/// Like [`simplify`], but consumes one fuel unit per bottom-up traversal pass
/// and stops early (returning the current expression) when the budget is
/// exhausted. Use this to bound runaway simplification of pathological inputs.
/// The returned [`Result`] is `Err` only if the fuel ran out before a fixpoint;
/// on success it is `Ok` with the simplified atom (which may still differ from
/// the input if some rules fired before exhaustion).
///
/// Pass [`Fuel::default`](ocas_core::fuel::Fuel::default) for an effectively
/// unlimited budget that still participates in nested accounting.
pub fn simplify_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    rules: &[Rule<'a>],
    iter_limit: usize,
    fuel: &Fuel,
) -> Result<Atom<'a>> {
    fuel.check()?;
    let mut current = transform(ctx, atom, |a| apply_rules(ctx, a, rules));
    for _ in 1..iter_limit {
        fuel.consume(1);
        fuel.check()?;
        let next = transform(ctx, current, |a| apply_rules(ctx, a, rules));
        if next == current {
            return Ok(next);
        }
        current = next;
    }
    Ok(current)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use crate::pattern::PatternAlloc;
    use crate::rules::default_rules;

    struct VecAlloc;

    impl<'a> PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(
            &self,
            _items: &[crate::pattern::Pattern<'a>],
        ) -> &'a [crate::pattern::Pattern<'a>] {
            unreachable!()
        }
    }

    #[test]
    fn simplify_mul_zero() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.num(0)]);
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn simplify_with_fuel_completes_when_budget_generous() {
        use ocas_core::fuel::Fuel;
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.num(0)]);
        let fuel = Fuel::new(50);
        let result = simplify_with_fuel(&ctx, expr, &rules, 10, &fuel).unwrap();
        assert_eq!(result.to_string(), "0");
        // Most of the budget should remain since simplification converges quickly.
        assert!(fuel.remaining() > 40);
    }

    #[test]
    fn simplify_with_fuel_errors_when_exhausted_upfront() {
        use ocas_core::error::OcasError;
        use ocas_core::fuel::Fuel;
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.num(0)]);
        let fuel = Fuel::new(0);
        let err = simplify_with_fuel(&ctx, expr, &rules, 10, &fuel).unwrap_err();
        assert!(matches!(err, OcasError::OutOfFuel), "got {err:?}");
    }

    #[test]
    fn simplify_add_same() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let expr = ctx.add(&[x, x]);
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "2*x");
    }

    #[test]
    fn simplify_nested() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let inner = ctx.mul(&[x, ctx.num(0)]);
        let expr = ctx.add(&[inner, x]);
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "x");
    }

    #[test]
    fn simplify_pow_zero() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let rules = default_rules(&ctx, &VecAlloc);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(0));
        let result = simplify(&ctx, expr, &rules, 10);
        assert_eq!(result.to_string(), "1");
    }
}
