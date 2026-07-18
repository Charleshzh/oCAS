//! The imaginary unit `I` and its rewrite rules.
//!
//! Complex intermediate forms arise when trigonometric integrands are
//! rewritten into exponentials before Risch integration. The imaginary
//! unit is represented as the constant `Var("I")`; the rules in
//! [`complex_rules`] reduce powers of `I` via `I^2 = -1`.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::matcher::MatchValue;
use ocas_rewrite::pattern::{Pattern, PatternAlloc, WildcardLevel};
use ocas_rewrite::rules::Rule;

/// Return the imaginary unit atom `I`.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_calc::complex::i;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// assert_eq!(i(&ctx).to_string(), "I");
/// ```
pub fn i<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    ctx.var("I")
}

/// Return whether `atom` is the imaginary unit `I`.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_calc::complex::{i, is_i};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// assert!(is_i(i(&ctx)));
/// assert!(!is_i(ctx.var("x")));
/// ```
pub fn is_i(atom: Atom) -> bool {
    matches!(atom.node(), AtomNode::Var(s) if s.as_str() == "I")
}

/// Build the rule set reducing powers and products of the imaginary unit.
///
/// The rules implement `I^2 = -1`:
///
/// - `I * I -> -1`
/// - `I^n -> 1`, `I`, `-1`, or `-I` according to `n mod 4`, for any
///   integer exponent `n`.
///
/// These rules are intended to be appended to the default simplifier rules
/// whenever complex intermediate forms may occur.
pub fn complex_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>> {
    let i_atom = i(ctx);

    // I * I -> -1
    let mul_rule = Rule::new(
        Pattern::Mul(
            alloc
                .alloc_slice(&[Pattern::Literal(i_atom), Pattern::Literal(i_atom)])
                .to_vec(),
        ),
        |_bindings, ctx| ctx.num(-1),
    );

    // I^n -> canonical value for integer n (n mod 4).
    let pow_rule = Rule::new(
        Pattern::Pow(Box::new((
            Pattern::Literal(i_atom),
            Pattern::Wildcard(Symbol::new("n"), WildcardLevel::Single),
        ))),
        |bindings, ctx| {
            let n = match bindings.get(Symbol::new("n")) {
                Some(MatchValue::Single(v)) => match v.node() {
                    AtomNode::Num(n) => *n,
                    _ => 1,
                },
                _ => 1,
            };
            match n.rem_euclid(4) {
                0 => ctx.num(1),
                1 => ctx.var("I"),
                2 => ctx.num(-1),
                _ => ctx.mul(&[ctx.num(-1), ctx.var("I")]),
            }
        },
    )
    .with_condition(|bindings| match bindings.get(Symbol::new("n")) {
        Some(MatchValue::Single(v)) => matches!(v.node(), AtomNode::Num(_)),
        _ => false,
    });

    vec![mul_rule, pow_rule]
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_rewrite::simplify::simplify;

    use super::*;

    fn run_rules<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
        let rules = complex_rules(ctx, &());
        simplify(ctx, expr, &rules, 20)
    }

    #[test]
    fn imaginary_unit_basics() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let i = i(&ctx);
        assert!(is_i(i));
        assert!(!is_i(ctx.var("x")));
        assert!(!is_i(ctx.num(2)));
    }

    #[test]
    fn i_times_i_is_minus_one() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.mul(&[i(&ctx), i(&ctx)]);
        assert_eq!(run_rules(&ctx, expr).to_string(), "-1");
    }

    #[test]
    fn i_squared_is_minus_one() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.pow(i(&ctx), ctx.num(2));
        assert_eq!(run_rules(&ctx, expr).to_string(), "-1");
    }

    #[test]
    fn i_powers_cycle_mod_four() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_eq!(
            run_rules(&ctx, ctx.pow(i(&ctx), ctx.num(3))).to_string(),
            "-1*I"
        );
        assert_eq!(
            run_rules(&ctx, ctx.pow(i(&ctx), ctx.num(4))).to_string(),
            "1"
        );
        assert_eq!(
            run_rules(&ctx, ctx.pow(i(&ctx), ctx.num(5))).to_string(),
            "I"
        );
    }

    #[test]
    fn i_negative_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // I^-1 = -I
        assert_eq!(
            run_rules(&ctx, ctx.pow(i(&ctx), ctx.num(-1))).to_string(),
            "-1*I"
        );
    }
}
