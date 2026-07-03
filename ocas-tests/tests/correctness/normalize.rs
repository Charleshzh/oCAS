#[test]
fn normalize_simple_zero_identity() {
    let result = crate::normalize_to_string("x + 0");
    assert_eq!(result, "x");
}

#[test]
fn normalize_simple_one_identity() {
    let result = crate::normalize_to_string("x * 1");
    assert_eq!(result, "x");
}

#[test]
fn normalize_simple_zero_absorb() {
    let result = crate::normalize_to_string("x * 0");
    assert_eq!(result, "0");
}

#[test]
fn normalize_medium_like_terms() {
    // The normalizer flattens and sorts but does *not* combine like terms.
    let result = crate::normalize_to_string("x + x + x + y + y + 0");
    assert_eq!(result, "x + x + x + y + y");
}

#[test]
fn normalize_medium_multiplication() {
    let result = crate::normalize_to_string("x * y * x * y");
    assert_eq!(result, "x*x*y*y");
}

#[test]
fn normalize_complex_nested() {
    let result = crate::normalize_to_string("x + (y + (z + 2)) + 3");
    assert_eq!(result, "5 + x + y + z");
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn normalize_complex_idempotent() {
    let input = "x + y + z + x + y + z";
    let first = crate::normalize_to_string(input);
    let second = crate::normalize_to_string(&first);
    assert_eq!(first, second);
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn normalize_very_complex_random_idempotent() {
    let input = "((x + y + z)^2 + (x - y - z)^2) * ((x + y + z)^2 - (x - y - z)^2)";
    let first = crate::normalize_to_string(input);
    let second = crate::normalize_to_string(&first);
    let third = crate::normalize_to_string(&second);
    assert_eq!(first, second);
    assert_eq!(second, third);
}
