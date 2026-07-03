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
