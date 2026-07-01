//! Rewrite rules for oCAS.
//!
//! A [`Rule`] pairs a pattern with a replacement builder. When the pattern
//! matches a sub-expression, the builder receives the wildcard bindings and
//! produces a replacement atom. Rules are applied by the [`simplify`]
//! function in a bottom-up traversal until no more changes occur.

use ocas_atom::{Atom, AtomArena};

use crate::matcher::{Bindings, MatchError, match_pattern};
use crate::pattern::{Pattern, PatternAlloc};

/// A rewrite rule.
///
/// Rules are typically created from parsed patterns via [`Rule::new`], or by
/// using the convenience constructors in the [`rules`](crate::rules) module.
#[derive(Clone)]
pub struct Rule<'a> {
    pattern: Pattern<'a>,
    replacement: fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
    condition: Option<fn(&Bindings<'a>) -> bool>,
}

impl<'a> Rule<'a> {
    /// Create a rule from a pattern and a replacement builder.
    ///
    /// The `replacement` function receives the bindings produced by a successful
    /// match and the arena context so it can construct new atoms.
    pub fn new(
        pattern: Pattern<'a>,
        replacement: fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
    ) -> Self {
        Self {
            pattern,
            replacement,
            condition: None,
        }
    }

    /// Add a condition that must hold for the rule to fire.
    ///
    /// The condition is evaluated after a successful match but before the
    /// replacement is built.
    pub fn with_condition(mut self, condition: fn(&Bindings<'a>) -> bool) -> Self {
        self.condition = Some(condition);
        self
    }

    /// Try to apply this rule to `atom`. Returns `Some` if the rule matched and
    /// the condition (if any) was satisfied.
    pub fn apply(&self, ctx: &AtomArena<'a>, atom: Atom<'a>) -> Option<Atom<'a>> {
        match match_pattern(self.pattern.clone(), atom) {
            Ok(bindings) => {
                if let Some(cond) = self.condition
                    && !cond(&bindings)
                {
                    return None;
                }
                Some((self.replacement)(&bindings, ctx))
            }
            Err(MatchError::NoMatch) | Err(MatchError::InconsistentBinding) => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Built-in algebraic rules
// ---------------------------------------------------------------------------

fn pattern_from_str<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
    s: &'a str,
) -> Pattern<'a> {
    let atom = ocas_parse::parse(ctx, s).expect("built-in rule pattern is valid");
    Pattern::from_atom(alloc, atom)
}

macro_rules! binding_single {
    ($bindings:expr, $name:expr) => {
        match $bindings.get(ocas_atom::Symbol::new($name)) {
            Some(crate::matcher::MatchValue::Single(atom)) => *atom,
            _ => panic!(concat!("expected single binding for '", $name, "'")),
        }
    };
}

/// `x + 0 -> x`
pub fn add_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ + 0"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `0 + x -> x`
pub fn add_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "0 + x_"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `x * 0 -> 0`
pub fn mul_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ * 0"), |_bindings, ctx| {
        ctx.num(0)
    })
}

/// `0 * x -> 0`
pub fn mul_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "0 * x_"), |_bindings, ctx| {
        ctx.num(0)
    })
}

/// `x * 1 -> x`
pub fn mul_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ * 1"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `1 * x -> x`
pub fn mul_one_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "1 * x_"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `x + x -> 2*x`
pub fn add_same<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ + x_"), |bindings, ctx| {
        let x = binding_single!(bindings, "x");
        ctx.mul(&[ctx.num(2), x])
    })
}

/// `x ^ 0 -> 1`
pub fn pow_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ ^ 0"), |_bindings, ctx| {
        ctx.num(1)
    })
}

/// `x ^ 1 -> x`
pub fn pow_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ ^ 1"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// Return the default set of algebraic rewrite rules.
pub fn default_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>> {
    vec![
        add_zero(ctx, alloc),
        add_zero_left(ctx, alloc),
        mul_zero(ctx, alloc),
        mul_zero_left(ctx, alloc),
        mul_one(ctx, alloc),
        mul_one_left(ctx, alloc),
        add_same(ctx, alloc),
        pow_zero(ctx, alloc),
        pow_one(ctx, alloc),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    struct VecAlloc;

    impl<'a> PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(&self, _items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
            // Not used by the matcher with the current Vec-based Pattern design.
            unreachable!()
        }
    }

    fn pat_atom<'a>(ctx: &'a AtomArena<'a>, alloc: &'a VecAlloc, s: &'a str) -> Pattern<'a> {
        let atom = ocas_parse::parse(ctx, s).unwrap();
        Pattern::from_atom(alloc, atom)
    }

    #[test]
    fn add_zero_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = add_zero(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.add(&[x, ctx.num(0)]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result, x);
    }

    #[test]
    fn mul_zero_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = mul_zero(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.mul(&[x, ctx.num(0)]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result, ctx.num(0));
    }

    #[test]
    fn add_same_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = add_same(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.add(&[x, x]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result.to_string(), "2*x");
    }

    #[test]
    fn rule_with_condition_can_block() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_");
        let rule = Rule::new(pat, |bindings, _ctx| binding_single!(bindings, "x"))
        .with_condition(|bindings| {
            matches!(
                bindings.get(ocas_atom::Symbol::new("x")),
                Some(crate::matcher::MatchValue::Single(a)) if matches!(a.node(), ocas_atom::AtomNode::Num(_))
            )
        });

        let x = ctx.var("x");
        assert!(rule.apply(&ctx, x).is_none());
        assert!(rule.apply(&ctx, ctx.num(5)).is_some());
    }
}
