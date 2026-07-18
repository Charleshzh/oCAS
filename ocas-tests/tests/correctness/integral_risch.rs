//! Risch-algorithm and special-function integration correctness tests.
//!
//! Every test verifies the antiderivative against SymPy's `integrate`
//! (skipped gracefully when `uv` is unavailable), so the suite doubles as
//! a compatibility check with the reference CAS.

// ------------------------------------------------------------------
//  Rational functions (Hermite + logarithmic part)
// ------------------------------------------------------------------

#[test]
fn risch_rational_inverse() {
    let result = crate::integrate_to_string("1/x", "x");
    crate::assert_eq_sympy(&result, "integrate", "1/x");
}

#[test]
fn risch_rational_atan() {
    let result = crate::integrate_to_string("1/(x^2 + 1)", "x");
    crate::assert_eq_sympy(&result, "integrate", "1/(x^2 + 1)");
}

#[test]
fn risch_rational_hermite_repeated() {
    let result = crate::integrate_to_string("1/(x + 1)^2", "x");
    crate::assert_eq_sympy(&result, "integrate", "1/(x + 1)^2");
}

#[test]
fn risch_rational_log_derivative() {
    let result = crate::integrate_to_string("(2*x + 3)/(x^2 + 3*x + 5)", "x");
    crate::assert_eq_sympy(&result, "integrate", "(2*x + 3)/(x^2 + 3*x + 5)");
}

// ------------------------------------------------------------------
//  Elementary transcendental (log/exp towers)
// ------------------------------------------------------------------

#[test]
fn risch_log_x() {
    let result = crate::integrate_to_string("log(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "log(x)");
}

#[test]
fn risch_log_over_x() {
    let result = crate::integrate_to_string("log(x)/x", "x");
    crate::assert_eq_sympy(&result, "integrate", "log(x)/x");
}

#[test]
fn risch_x_log_x() {
    let result = crate::integrate_to_string("x*log(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "x*log(x)");
}

#[test]
fn risch_exp_x() {
    let result = crate::integrate_to_string("exp(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "exp(x)");
}

#[test]
fn risch_x_exp_x() {
    let result = crate::integrate_to_string("x*exp(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "x*exp(x)");
}

#[test]
fn risch_x2_exp_x() {
    let result = crate::integrate_to_string("x^2*exp(x)", "x");
    crate::assert_eq_sympy(&result, "integrate", "x^2*exp(x)");
}

#[test]
fn risch_x_exp_x2() {
    let result = crate::integrate_to_string("x*exp(x^2)", "x");
    crate::assert_eq_sympy(&result, "integrate", "x*exp(x^2)");
}

// ------------------------------------------------------------------
//  Special functions (Meijer-G endpoints)
// ------------------------------------------------------------------

#[test]
fn risch_exp_neg_x_squared_erf() {
    let result = crate::integrate_to_string("exp(-x^2)", "x");
    crate::assert_eq_sympy(&result, "integrate", "exp(-x^2)");
}

#[test]
fn risch_exp_x_over_x_ei() {
    let result = crate::integrate_to_string("exp(x)/x", "x");
    crate::assert_eq_sympy(&result, "integrate", "exp(x)/x");
}

#[test]
fn risch_sin_over_x_si() {
    let result = crate::integrate_to_string("sin(x)/x", "x");
    crate::assert_eq_sympy(&result, "integrate", "sin(x)/x");
}

#[test]
fn risch_cos_over_x_ci() {
    let result = crate::integrate_to_string("cos(x)/x", "x");
    crate::assert_eq_sympy(&result, "integrate", "cos(x)/x");
}
