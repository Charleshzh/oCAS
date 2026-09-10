#[test]
fn parse_simple() {
    let result = crate::parse_to_string("x + 1");
    assert_eq!(result, "x + 1");
}

#[test]
fn parse_medium() {
    let result = crate::parse_to_string("x^2 + 2*x + 1");
    assert_eq!(result, "((x^2) + (2*x)) + 1");
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn parse_complex() {
    let result = crate::parse_to_string("(x + y + z)^3 - 3*(x + y)*(y + z)*(z + x)");
    // The parser represents n-ary sums as left-nested binary additions.
    assert!(result.contains("((x + y) + z)"));
    assert!(result.contains("(x + y)"));
    assert!(result.contains("(y + z)"));
    assert!(result.contains("(z + x)"));
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn parse_very_complex() {
    let input = "sin(cos(tan(log(exp(sqrt(x))))))";
    let result = crate::parse_to_string(input);
    assert_eq!(result, "sin(cos(tan(log(exp(sqrt(x))))))");
}

#[test]
fn parse_roundtrip_simple() {
    let input = "x^2 + 2*x + 1";
    let first = crate::parse_to_string(input);
    let second = crate::parse_to_string(&first);
    assert_eq!(first, second);
}

#[test]
fn parse_roundtrip_medium() {
    let input = "sin(x) + cos(y) * (x^2 + 1)";
    let first = crate::parse_to_string(input);
    let second = crate::parse_to_string(&first);
    assert_eq!(first, second);
}

// ---------------------------------------------------------------------------
// Unary minus and bare binary subtraction.
//
// The lexer used to carry the sign in the integer token (`-?[0-9]+`), so `2-1`
// lexed as the two adjacent operands `Integer(2), Integer(-1)` and failed to
// parse, while `a-b` worked only because `b` is not a digit. The integer token
// is now unsigned and the parser owns the sign; `crate::parse_to_string` and
// `crate::normalize_to_string` unwrap, so every call below also asserts that
// parsing succeeds.
// ---------------------------------------------------------------------------

#[test]
fn parse_subtraction_without_spaces() {
    // normalize sorts `Add` terms, so the negative constant comes first.
    assert_eq!(crate::normalize_to_string("2-1"), "1");
    assert_eq!(crate::normalize_to_string("x-1"), "-1 + x");
    assert_eq!(crate::normalize_to_string("x^2-1"), "-1 + (x^2)");
    assert_eq!(crate::normalize_to_string("a1-b2"), "a1 + (-1*b2)");
    assert_eq!(
        crate::normalize_to_string("(x+1)-(x-1)"),
        "1 + x + (-1*(-1 + x))"
    );
    // Spellings that already worked must keep their printed form.
    assert_eq!(crate::normalize_to_string("x - 1"), "-1 + x");
    assert_eq!(crate::normalize_to_string("x- 1"), "-1 + x");
    assert_eq!(crate::normalize_to_string("a-b"), "a + (-1*b)");
}

#[test]
fn parse_unary_minus_forms() {
    assert_eq!(crate::parse_to_string("-x"), "-1*x");
    assert_eq!(crate::parse_to_string("-7"), "-7");
    assert_eq!(crate::parse_to_string("- 2"), "-2");
    assert_eq!(crate::parse_to_string("-(x+1)"), "-1*(x + 1)");
    assert_eq!(crate::parse_to_string("sin(-x)"), "sin(-1*x)");
    assert_eq!(crate::parse_to_string("f(-1, x)"), "f(-1, x)");
}

#[test]
fn parse_negative_exponent_and_factor() {
    assert_eq!(crate::parse_to_string("x^-1"), "x^-1");
    assert_eq!(crate::parse_to_string("x^(-1)"), "x^-1");
    assert_eq!(crate::parse_to_string("2^-1"), "2^-1");
    assert_eq!(crate::parse_to_string("x*-2"), "x*-2");
    assert_eq!(crate::parse_to_string("2*-1"), "2*-1");
    assert_eq!(crate::parse_to_string("x--1"), "x + (-1*-1)");
}

#[test]
fn parse_unary_minus_precedence() {
    // `^` binds tighter than negation, so `-x^2` is `-(x^2)`.
    assert_eq!(crate::parse_to_string("-x^2"), "-1*(x^2)");
    assert_eq!(crate::normalize_to_string("-x^2"), "-1*(x^2)");
    // `^` stays right-associative.
    assert_eq!(crate::parse_to_string("2^3^2"), "2^(3^2)");
    // A sign written directly on an integer literal stays part of the literal:
    // the printer emits a negative numeric base without parentheses, so
    // `-2^2` has to read back as `(-2)^2` for printing to stay lossless.
    assert_eq!(crate::parse_to_string("-2^2"), "-2^2");
    assert_eq!(crate::parse_to_string("(-2)^2"), "-2^2");
}

#[test]
fn parse_signed_forms_roundtrip() {
    let inputs = [
        "-x",
        "-7",
        "- 2",
        "--7",
        "-(x+1)",
        "sin(-x)",
        "f(-1, x)",
        "x^-1",
        "x^(-1)",
        "2^-1",
        "x*-2",
        "2*-1",
        "x--1",
        "2-1",
        "x-1",
        "x^2-1",
        "a1-b2",
        "(x+1)-(x-1)",
        "-x^2",
        "-2^2",
        "1/(x^2-1)^2",
    ];
    for input in inputs {
        let first = crate::parse_to_string(input);
        let second = crate::parse_to_string(&first);
        assert_eq!(first, second, "printing {input} did not round-trip");
    }
}

#[test]
fn parse_negative_denominator() {
    // End-to-end shape that used to be `PARSE_ERR` in the integration harness.
    assert_eq!(crate::normalize_to_string("1/(x^2-1)^2"), "(-1 + (x^2))^-2");
}
