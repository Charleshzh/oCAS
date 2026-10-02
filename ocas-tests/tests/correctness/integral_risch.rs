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
//  Hyperexponential levels: Laurent (negative-power) layers (0.29.0)
// ------------------------------------------------------------------

/// The 0.28.0 latent bug's minimal repro: the exp-level rational part
/// mis-scaled negative powers of `t = exp(x)` (2× on `t⁻²`, 4/3× on
/// `t⁻⁴`). 0.28.0 converted it to an honest fallback via the tower check;
/// 0.29.0 solves it through the Laurent split and the answer is verified
/// against SymPy.
#[test]
fn risch_exp_negative_power_layers() {
    let result = crate::integrate_to_string("(exp(x)^2 + 1)^3/(8*exp(x)^4)", "x");
    assert!(
        !result.contains("Integral("),
        "the latent-bug repro must now be solved, got: {result}"
    );
    crate::assert_eq_sympy(&result, "integrate", "(exp(x)^2 + 1)^3/(8*exp(x)^4)");
}

#[test]
fn risch_exp_reciprocal_layer() {
    // ∫ (1 + exp(x))/exp(x) dx = x − exp(−x): a single negative layer.
    let result = crate::integrate_to_string("(1 + exp(x))/exp(x)", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "(1 + exp(x))/exp(x)");
}

/// rubi-01524's shape: the dependent generators `exp(x³)` / `exp(4x³)`
/// merge (0.29.0), leaving a polynomial in `t` at the hyperexp level.
#[test]
fn risch_exp_integer_multiple_layers() {
    let result = crate::integrate_to_string("exp(x^3)*(1 - exp(4*x^3))^2*x^2", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "exp(x^3)*(1 - exp(4*x^3))^2*x^2");
}

// ------------------------------------------------------------------
//  Rothstein–Trager logarithmic part (0.29.0)
// ------------------------------------------------------------------

/// Single rational root at an exp level, with a `k`-valued remainder:
/// ∫ dx/(1 + exp(x)) = x − log(1 + exp(x)).
#[test]
fn risch_rt_single_root_exp() {
    let result = crate::integrate_to_string("1/(1 + exp(x))", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "1/(1 + exp(x))");
}

/// Double root counted once: ∫ dx/(exp(2x) − 1) = log(exp(2x) − 1)/2 − x.
#[test]
fn risch_rt_double_root_exp() {
    let result = crate::integrate_to_string("1/(exp(2*x) - 1)", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "1/(exp(2*x) - 1)");
}

/// rubi-00184: dependent exponential pair merged, then RT roots ±1/2.
#[test]
fn risch_rt_rubi_00184() {
    let result = crate::integrate_to_string("exp(x)/(1 - exp(2*x))", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "exp(x)/(1 - exp(2*x))");
}

/// Two rational roots at a log level, with the resultant coefficients in
/// the lower field: ∫ dx/(x·log(x)·(log(x)+1)) = log(log x) − log(log x + 1).
#[test]
fn risch_rt_two_roots_log() {
    let result = crate::integrate_to_string("1/(x*log(x)*(log(x) + 1))", "x");
    assert!(!result.contains("Integral("), "got: {result}");
    crate::assert_eq_sympy(&result, "integrate", "1/(x*log(x)*(log(x) + 1))");
}

/// A non-constant resultant root means non-elementary: the engine must
/// decline honestly (SymPy leaves it unevaluated too).
#[test]
fn risch_rt_non_constant_root_declines() {
    let result = crate::integrate_to_string("1/(exp(x) + x)", "x");
    assert!(
        result.contains("Integral("),
        "non-elementary case must fall back, got: {result}"
    );
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
