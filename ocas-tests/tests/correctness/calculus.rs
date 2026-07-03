#[test]
fn diff_simple_power() {
    let result = crate::diff_to_string("x^2", "x");
    assert_eq!(result, "2*x");
    crate::assert_eq_sympy(&result, "diff", "x^2");
}

#[test]
fn diff_simple_sin() {
    let result = crate::diff_to_string("sin(x)", "x");
    assert_eq!(result, "cos(x)");
    crate::assert_eq_sympy(&result, "diff", "sin(x)");
}

#[test]
fn diff_medium_product() {
    let result = crate::diff_to_string("x*sin(x)", "x");
    crate::assert_eq_sympy(&result, "diff", "x*sin(x)");
}

#[test]
fn diff_medium_chain() {
    let result = crate::diff_to_string("sin(x^2)", "x");
    crate::assert_eq_sympy(&result, "diff", "sin(x^2)");
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn diff_complex_nested() {
    let result = crate::diff_to_string("sin(exp(x^2))*cos(log(x))", "x");
    crate::assert_eq_sympy(&result, "diff", "sin(exp(x^2))*cos(log(x))");
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn diff_very_complex_high_order() {
    let result = crate::diff_to_string("x^10*exp(x)", "x");
    crate::assert_eq_sympy(&result, "diff", "x^10*exp(x)");
}

#[test]
fn integrate_simple_power() {
    let result = crate::integrate_to_string("x^2", "x");
    assert_eq!(result, "(3^-1)*(x^3)");
    crate::assert_eq_sympy(&result, "integrate", "x^2");
}

#[test]
fn integrate_simple_sin() {
    let result = crate::integrate_to_string("sin(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "sin(x)");
}

#[test]
fn integrate_medium_by_parts() {
    let result = crate::integrate_to_string("x*sin(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "x*sin(x)");
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn integrate_complex_rational() {
    let result = crate::integrate_to_string("1/(x^2+1)", "x");
    crate::assert_eq_sympy(&result, "integrate", "1/(x^2+1)");
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn integrate_very_complex_known_gap() {
    // This is expected to fail or return an Integral form; it documents a gap.
    let result = crate::integrate_to_string("exp(-x^2)", "x");
    assert!(result.starts_with("Integral"));
}

#[test]
fn taylor_simple_exp() {
    let result = crate::taylor_to_string("exp(x)", "x", 0, 3);
    assert_eq!(result, "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))");
    crate::assert_eq_sympy(&result, "series", "exp(x):3");
}

#[test]
fn taylor_simple_sin() {
    let result = crate::taylor_to_string("sin(x)", "x", 0, 5);
    crate::assert_eq_sympy(&result, "series", "sin(x):5");
}

#[test]
fn taylor_medium_product() {
    let result = crate::taylor_to_string("cos(x)*exp(x)", "x", 0, 5);
    crate::assert_eq_sympy(&result, "series", "cos(x)*exp(x):5");
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn taylor_complex_tan() {
    let result = crate::taylor_to_string("tan(x)", "x", 0, 7);
    crate::assert_eq_sympy(&result, "series", "tan(x):7");
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn taylor_very_complex_nested() {
    let result = crate::taylor_to_string("exp(sin(x))", "x", 0, 10);
    crate::assert_eq_sympy(&result, "series", "exp(sin(x)):10");
}
