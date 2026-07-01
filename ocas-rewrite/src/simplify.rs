//! Rule-based simplification for oCAS.
//!
//! The [`simplify`] function applies a list of rewrite rules to an atom tree in
//! a bottom-up traversal, repeatedly, until no rule fires or a limit is
//! reached. It composes with the [`Transformer`](crate::transformer) and
//! [`Rule`](crate::rules::Rule) machinery.

use ocas_atom::{Atom, AtomArena};

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
    for _ in 0..iter_limit {
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
