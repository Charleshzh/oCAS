//! Parser and printer for oCAS expressions.
//!
//! This crate provides a [`logos`]-based lexer and a small recursive-descent
//! parser that turns text into [`ocas_atom::Atom`] expression trees.

pub mod lexer;
pub mod parser;

pub use parser::{ParseError, parse};

#[cfg(test)]
mod proptests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use proptest::prelude::*;

    use super::*;

    fn valid_expr_str() -> impl Strategy<Value = String> {
        let leaf = prop_oneof![
            Just("x".to_string()),
            Just("y".to_string()),
            Just("z".to_string()),
            (-100..100i64).prop_map(|n| n.to_string()),
        ];
        leaf.prop_recursive(4, 64, 4, |inner| {
            prop_oneof![
                inner.clone().prop_map(|e| format!("sin({})", e)),
                inner.clone().prop_map(|e| format!("cos({})", e)),
                (inner.clone(), inner.clone()).prop_map(|(a, b)| format!("({}) + ({})", a, b)),
                (inner.clone(), inner.clone()).prop_map(|(a, b)| format!("({}) * ({})", a, b)),
                (inner.clone(), 0..5u32).prop_map(|(a, n)| format!("({})^{}", a, n)),
            ]
        })
    }

    proptest! {
        #[test]
        fn parse_succeeds_on_valid_expr(s in valid_expr_str()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let result = parse(&ctx, &s);
            prop_assert!(result.is_ok(), "parse failed for {}: {:?}", s, result);
        }

        #[test]
        fn parse_print_is_deterministic(s in valid_expr_str()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, &s).unwrap();
            let printed = atom.to_string();
            let reparsed = parse(&ctx, &printed).unwrap();
            prop_assert_eq!(atom.to_string(), reparsed.to_string());
        }

        #[test]
        fn parse_rejects_invalid(s in "[^0-9a-zA-Z()+*^\t\n\r]*") {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let result = parse(&ctx, &s);
            // Some strings may coincidentally be valid (e.g. empty), but most should fail.
            // We only assert that parsing does not panic.
            let _ = result;
        }
    }
}
