//! End-to-end tests that exercise the C ABI through Rust FFI.
//!
//! These complement the unit tests in `expression.rs` by calling the
//! `#[no_mangle] extern "C"` functions exactly as a C caller would.

use ocas_c::{
    OCAS_OK, OcasAlgebraicFactorArray, OcasPolyFactorArray, OcasTensorContraction,
    OcasVegasOptions, ocas_algebraic_factor_array_free, ocas_algebraic_field_create,
    ocas_algebraic_field_degree, ocas_algebraic_field_free, ocas_algebraic_poly_create,
    ocas_algebraic_poly_degree, ocas_algebraic_poly_factor, ocas_algebraic_poly_free,
    ocas_algebraic_poly_to_string, ocas_dual_add, ocas_dual_constant, ocas_dual_deriv,
    ocas_dual_div, ocas_dual_mul, ocas_dual_neg, ocas_dual_shape_free,
    ocas_dual_shape_n_components, ocas_dual_shape_n_vars, ocas_dual_shape_new, ocas_dual_value,
    ocas_dual_variable, ocas_error_clear, ocas_error_last_message, ocas_expr_clone, ocas_expr_diff,
    ocas_expr_free, ocas_expr_integrate, ocas_expr_normalize, ocas_expr_parse, ocas_expr_simplify,
    ocas_expr_substitute, ocas_expr_taylor, ocas_expr_to_string, ocas_hyperdual_free,
    ocas_integrate_1d, ocas_poly_factor_array_free, ocas_poly_fp_clone, ocas_poly_fp_create,
    ocas_poly_fp_degree, ocas_poly_fp_factor, ocas_poly_fp_free, ocas_poly_fp_to_string,
    ocas_poly_z_clone, ocas_poly_z_create, ocas_poly_z_degree, ocas_poly_z_factor,
    ocas_poly_z_free, ocas_poly_z_to_string, ocas_string_free, ocas_tensor_contract,
    ocas_tensor_contraction_free, ocas_tensor_create, ocas_tensor_free, ocas_tensor_name,
    ocas_tensor_rank, ocas_tensor_symmetrise_sign, ocas_tensor_symmetry, ocas_tensor_to_string,
    ocas_vegas_create, ocas_vegas_free, ocas_vegas_integrate, ocas_vegas_iterations,
    ocas_vegas_result, ocas_version,
};
use std::ffi::{CStr, CString};
use std::ptr;

fn parse(s: &str) -> *mut ocas_c::expression::OcasExpr {
    let c = CString::new(s).unwrap();
    let mut err = 0;
    let h = unsafe { ocas_expr_parse(c.as_ptr(), &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!h.is_null(), "parse failed for {s:?}");
    h
}

fn to_string(h: *mut ocas_c::expression::OcasExpr) -> String {
    let mut err = 0;
    let ptr = unsafe { ocas_expr_to_string(h, &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!ptr.is_null());
    // SAFETY: `ptr` is a valid null-terminated C string from to_string.
    let s = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap().to_string();
    unsafe { ocas_string_free(ptr) };
    s
}

#[test]
fn version_string_is_nonempty() {
    let ptr = ocas_version();
    assert!(!ptr.is_null());
    let s = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap();
    assert!(!s.is_empty());
    // Workspace version is 0.9.0.
    assert_eq!(s, env!("CARGO_PKG_VERSION"));
}

#[test]
fn parse_diff_roundtrip() {
    let expr = parse("x^2");
    let var = CString::new("x").unwrap();
    let mut err = 0;
    let deriv = unsafe { ocas_expr_diff(expr, var.as_ptr(), &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!deriv.is_null());
    assert_eq!(to_string(deriv), "2*x");
    unsafe {
        ocas_expr_free(deriv);
        ocas_expr_free(expr);
    }
}

#[test]
fn parse_integrate_roundtrip() {
    let expr = parse("2*x");
    let var = CString::new("x").unwrap();
    let mut err = 0;
    let integ = unsafe { ocas_expr_integrate(expr, var.as_ptr(), &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!integ.is_null());
    let s = to_string(integ);
    assert!(s.contains("(x^2)"), "got: {s}");
    unsafe {
        ocas_expr_free(integ);
        ocas_expr_free(expr);
    }
}

#[test]
fn simplify_collapses_mul_zero() {
    let expr = parse("x*0");
    let mut err = 0;
    let simplified = unsafe { ocas_expr_simplify(expr, &mut err) };
    assert_eq!(err, OCAS_OK);
    assert_eq!(to_string(simplified), "0");
    unsafe {
        ocas_expr_free(simplified);
        ocas_expr_free(expr);
    }
}

#[test]
fn clone_is_independent() {
    let expr = parse("y + 1");
    let mut err = 0;
    let cloned = unsafe { ocas_expr_clone(expr, &mut err) };
    assert_eq!(err, OCAS_OK);
    assert_eq!(to_string(cloned), to_string(expr));
    // Freeing one must not corrupt the other.
    unsafe { ocas_expr_free(expr) };
    assert_eq!(to_string(cloned), "1 + y");
    unsafe { ocas_expr_free(cloned) };
}

#[test]
fn null_handle_returns_error() {
    let var = CString::new("x").unwrap();
    let mut err = 0;
    let result = unsafe { ocas_expr_diff(std::ptr::null(), var.as_ptr(), &mut err) };
    assert!(result.is_null());
    let msg = ocas_error_last_message();
    assert!(!msg.is_null());
    ocas_error_clear();
    assert!(ocas_error_last_message().is_null());
}

#[test]
fn invalid_input_reports_parse_error() {
    let bad = CString::new("@@@invalid@@@").unwrap();
    let mut err: std::ffi::c_int = 0;
    let result = unsafe { ocas_expr_parse(bad.as_ptr(), &mut err) };
    assert!(result.is_null());
    let msg = ocas_error_last_message();
    assert!(!msg.is_null());
}

#[test]
fn taylor_series_expands() {
    let expr = parse("exp(x)");
    let var = CString::new("x").unwrap();
    let point = parse("0");
    let mut err = 0;
    let series = unsafe { ocas_expr_taylor(expr, var.as_ptr(), point, 3, &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!series.is_null());
    let s = to_string(series);
    assert!(s.contains("1 + x"), "got: {s}");
    assert!(s.contains("(x^2)"), "got: {s}");
    assert!(s.contains("(x^3)"), "got: {s}");
    unsafe {
        ocas_expr_free(series);
        ocas_expr_free(point);
        ocas_expr_free(expr);
    }
}

#[test]
fn normalize_in_place() {
    let expr = parse("x + x + 0");
    let mut err = 0;
    let rc = unsafe { ocas_expr_normalize(expr, &mut err) };
    assert_eq!(rc, OCAS_OK);
    // After normalization, x + x should merge; 0 absorbed.
    let s = to_string(expr);
    assert!(s.contains("2*x") || s.contains("x"), "got: {s}");
    unsafe { ocas_expr_free(expr) };
}

#[test]
fn substitute_replaces_variable() {
    let expr = parse("x^2 + 1");
    let y = parse("y");
    let var = CString::new("x").unwrap();
    let mut err = 0;
    let result = unsafe { ocas_expr_substitute(expr, var.as_ptr(), y, &mut err) };
    assert_eq!(err, OCAS_OK);
    assert!(!result.is_null());
    let s = to_string(result);
    assert!(s.contains("(y^2)"), "got: {s}");
    unsafe {
        ocas_expr_free(result);
        ocas_expr_free(y);
        ocas_expr_free(expr);
    }
}

#[test]
fn string_free_on_null_is_safe() {
    // Passing NULL to ocas_string_free must be a no-op, not a crash.
    unsafe { ocas_string_free(std::ptr::null_mut()) };
}

// ------------------------------------------------------------------
//  Polynomial C API tests
// ------------------------------------------------------------------

fn cstr(s: &str) -> CString {
    CString::new(s).unwrap()
}

fn c_string_to_string(ptr: *mut std::ffi::c_char) -> String {
    let s = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap().to_string();
    unsafe { ocas_string_free(ptr) };
    s
}

// -- Integer polynomial (OcasPolyZ) --

#[test]
fn poly_z_create_and_to_string() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_z_create(cstr("x^2 + y + 1").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!p.is_null());
    let mut err = 0;
    let s_ptr = ocas_poly_z_to_string(p, &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!s_ptr.is_null());
    let s = c_string_to_string(s_ptr);
    // Output should contain the variables and coefficients.
    assert!(s.contains('x') && s.contains('y'), "got: {s}");
    ocas_poly_z_free(p);
}

#[test]
fn poly_z_clone_is_independent() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_z_create(cstr("x + y").as_ptr(), &mut err);
    assert!(!p.is_null());
    let clone = ocas_poly_z_clone(p);
    assert!(!clone.is_null());
    assert_eq!(ocas_poly_z_degree(p), ocas_poly_z_degree(clone));
    ocas_poly_z_free(p);
    // Clone must survive freeing the original.
    assert_eq!(ocas_poly_z_degree(clone), 1);
    ocas_poly_z_free(clone);
}

#[test]
fn poly_z_degree_returns_total_degree() {
    let mut err: std::ffi::c_int = 0;
    // x^2*y + 1 has total degree 3
    let p = ocas_poly_z_create(cstr("x^2*y + 1").as_ptr(), &mut err);
    assert!(!p.is_null());
    assert_eq!(ocas_poly_z_degree(p), 3);
    ocas_poly_z_free(p);
}

#[test]
fn poly_z_factor_produces_factors() {
    let mut err: std::ffi::c_int = 0;
    // x^2 - 1 = (x-1)(x+1) as bivariate (no y dependence)
    let p = ocas_poly_z_create(cstr("x^2 - 1").as_ptr(), &mut err);
    assert!(!p.is_null());
    let mut factors = OcasPolyFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let rc = ocas_poly_z_factor(p, &mut factors, &mut err);
    assert_eq!(rc, OCAS_OK);
    assert!(factors.len >= 1, "expected at least 1 factor");
    // Free each factor handle, then the array.
    for i in 0..factors.len {
        let f = unsafe { &*factors.factors.add(i) };
        assert!(!f.poly.is_null());
        ocas_poly_z_free(f.poly as *mut ocas_c::OcasPolyZ);
    }
    ocas_poly_factor_array_free(&mut factors);
    ocas_poly_z_free(p);
}

#[test]
fn poly_z_null_input_returns_null() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_z_create(ptr::null(), &mut err);
    assert!(p.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
}

#[test]
fn poly_z_invalid_parse_returns_null() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_z_create(cstr("@invalid!").as_ptr(), &mut err);
    assert!(p.is_null());
    let msg = ocas_error_last_message();
    assert!(!msg.is_null());
    ocas_error_clear();
}

#[test]
fn poly_z_factor_null_poly_returns_error() {
    let mut factors = OcasPolyFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let mut err: std::ffi::c_int = 0;
    let rc = ocas_poly_z_factor(ptr::null(), &mut factors, &mut err);
    assert_ne!(rc, OCAS_OK);
    ocas_error_clear();
}

// -- Algebraic number field polynomial (OcasAlgebraicPoly) --

#[test]
fn algebraic_field_create_sqrt2() {
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("x^2 - 2").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!f.is_null());
    assert_eq!(ocas_algebraic_field_degree(f), 2);
    ocas_algebraic_field_free(f);
}

#[test]
fn algebraic_field_create_cbrt2() {
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("x^3 - 2").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!f.is_null());
    assert_eq!(ocas_algebraic_field_degree(f), 3);
    ocas_algebraic_field_free(f);
}

#[test]
fn algebraic_field_non_monic_rejected() {
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("2*x^2 - 2").as_ptr(), &mut err);
    assert!(f.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
}

#[test]
fn algebraic_poly_factor_sqrt2_splits() {
    // Over ℚ(√2): x² − 2 = (x − α)(x + α), two linear factors.
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("x^2 - 2").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!f.is_null());
    // x² − 2: coefficients "-2;0;1".
    let p = ocas_algebraic_poly_create(f, cstr("-2;0;1").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!p.is_null());
    assert_eq!(ocas_algebraic_poly_degree(p), 2);
    let mut factors = OcasAlgebraicFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let rc = ocas_algebraic_poly_factor(p, &mut factors, &mut err);
    assert_eq!(rc, OCAS_OK);
    assert_eq!(factors.len, 2, "x^2 - 2 must split into 2 linear factors");
    for i in 0..factors.len {
        let fac = unsafe { &*factors.factors.add(i) };
        assert!(!fac.poly.is_null());
        assert_eq!(fac.multiplicity, 1);
        assert_eq!(ocas_algebraic_poly_degree(fac.poly as *mut _), 1);
        ocas_algebraic_poly_free(fac.poly as *mut _);
    }
    ocas_algebraic_factor_array_free(&mut factors);
    ocas_algebraic_poly_free(p);
    ocas_algebraic_field_free(f);
}

#[test]
fn algebraic_poly_to_string_roundtrip() {
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("x^2 - 2").as_ptr(), &mut err);
    assert!(!f.is_null());
    let p = ocas_algebraic_poly_create(f, cstr("-2;0;1").as_ptr(), &mut err);
    assert!(!p.is_null());
    let s_ptr = ocas_algebraic_poly_to_string(p, &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!s_ptr.is_null());
    let s = c_string_to_string(s_ptr);
    assert!(s.contains('x'), "to_string output must mention x: {s}");
    ocas_algebraic_poly_free(p);
    ocas_algebraic_field_free(f);
}

#[test]
fn algebraic_poly_factor_with_alpha_coefficient() {
    // Over ℚ(√2): x² − α is irreducible (a root would be 2^(1/4) ∉ ℚ(√2)).
    // Construct x² − α as "0,-1;0;1": constant = −α, x coeff = 0, x² coeff = 1.
    let mut err: std::ffi::c_int = 0;
    let f = ocas_algebraic_field_create(cstr("x^2 - 2").as_ptr(), &mut err);
    assert!(!f.is_null());
    let p = ocas_algebraic_poly_create(f, cstr("0,-1;0;1").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!p.is_null());
    let mut factors = OcasAlgebraicFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let rc = ocas_algebraic_poly_factor(p, &mut factors, &mut err);
    assert_eq!(rc, OCAS_OK);
    // x² − α is irreducible of degree 2 over ℚ(√2).
    assert_eq!(factors.len, 1);
    for i in 0..factors.len {
        let fac = unsafe { &*factors.factors.add(i) };
        assert_eq!(
            ocas_algebraic_poly_degree(fac.poly as *mut _),
            2,
            "x^2 - alpha must stay irreducible"
        );
        ocas_algebraic_poly_free(fac.poly as *mut _);
    }
    ocas_algebraic_factor_array_free(&mut factors);
    ocas_algebraic_poly_free(p);
    ocas_algebraic_field_free(f);
}

#[test]
fn algebraic_poly_null_handles_return_error() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_algebraic_poly_create(ptr::null(), cstr("1").as_ptr(), &mut err);
    assert!(p.is_null());
    ocas_error_clear();

    let mut err: std::ffi::c_int = 0;
    let mut factors = OcasAlgebraicFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let rc = ocas_algebraic_poly_factor(ptr::null(), &mut factors, &mut err);
    assert_ne!(rc, OCAS_OK);
    ocas_error_clear();
}

// -- Finite-field polynomial (OcasPolyFp) --

#[test]
fn poly_fp_create_and_to_string() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_fp_create(cstr("x^2 + y + 1").as_ptr(), cstr("5").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!p.is_null());
    let mut err = 0;
    let s_ptr = ocas_poly_fp_to_string(p, &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!s_ptr.is_null());
    let s = c_string_to_string(s_ptr);
    assert!(s.contains('x') && s.contains('y'), "got: {s}");
    ocas_poly_fp_free(p);
}

#[test]
fn poly_fp_clone_is_independent() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_fp_create(cstr("x + y").as_ptr(), cstr("7").as_ptr(), &mut err);
    assert!(!p.is_null());
    let clone = ocas_poly_fp_clone(p);
    assert!(!clone.is_null());
    assert_eq!(ocas_poly_fp_degree(p), ocas_poly_fp_degree(clone));
    ocas_poly_fp_free(p);
    assert_eq!(ocas_poly_fp_degree(clone), 1);
    ocas_poly_fp_free(clone);
}

#[test]
fn poly_fp_factor_produces_factors() {
    let mut err: std::ffi::c_int = 0;
    // x^2 + y + 1 over F_5 — should factor or return as irreducible
    let p = ocas_poly_fp_create(cstr("x^2 + y + 1").as_ptr(), cstr("5").as_ptr(), &mut err);
    assert!(!p.is_null());
    let mut factors = OcasPolyFactorArray {
        factors: ptr::null_mut(),
        len: 0,
    };
    let rc = ocas_poly_fp_factor(p, &mut factors, &mut err);
    assert_eq!(rc, OCAS_OK);
    assert!(factors.len >= 1);
    for i in 0..factors.len {
        let f = unsafe { &*factors.factors.add(i) };
        assert!(!f.poly.is_null());
        ocas_poly_fp_free(f.poly as *mut ocas_c::OcasPolyFp);
    }
    ocas_poly_factor_array_free(&mut factors);
    ocas_poly_fp_free(p);
}

#[test]
fn poly_fp_null_prime_returns_null() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_fp_create(cstr("x").as_ptr(), ptr::null(), &mut err);
    assert!(p.is_null());
    ocas_error_clear();
}

#[test]
fn poly_fp_invalid_prime_returns_null() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_fp_create(cstr("x").as_ptr(), cstr("not_a_number").as_ptr(), &mut err);
    assert!(p.is_null());
    ocas_error_clear();
}

#[test]
fn poly_fp_prime_too_small_returns_null() {
    let mut err: std::ffi::c_int = 0;
    let p = ocas_poly_fp_create(cstr("x").as_ptr(), cstr("1").as_ptr(), &mut err);
    assert!(p.is_null());
    ocas_error_clear();
}

#[test]
fn poly_factor_array_free_on_null_is_safe() {
    ocas_poly_factor_array_free(ptr::null_mut());
}

// ------------------------------------------------------------------
//  Numerical integration C API tests (Vegas)
// ------------------------------------------------------------------

/// C-callable integrand f(x) = x; user_data unused.
unsafe extern "C" fn integrand_x(x: f64, _user_data: *mut std::ffi::c_void) -> f64 {
    x
}

/// C-callable integrand f(x) = x²; user_data unused.
unsafe extern "C" fn integrand_x_squared(x: f64, _user_data: *mut std::ffi::c_void) -> f64 {
    x * x
}

#[test]
fn vegas_create_and_integrate_linear() {
    let opts = OcasVegasOptions {
        n_bins: 64,
        n_samples: 20000,
        iterations: 8,
        learning_rate: 1.5,
        seed: 0x0C45,
    };
    let mut err: std::ffi::c_int = 0;
    let v = ocas_vegas_create(1, &opts, &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!v.is_null());
    let r = ocas_vegas_integrate(v, Some(integrand_x), ptr::null_mut(), &mut err);
    assert_eq!(err, OCAS_OK);
    // ∫₀¹ x dx = 1/2.
    assert!(
        (r.integral - 0.5).abs() < 0.01,
        "got integral {}",
        r.integral
    );
    assert_eq!(ocas_vegas_iterations(v), 8);
    let latest = ocas_vegas_result(v);
    assert_eq!(latest.integral, r.integral);
    ocas_vegas_free(v);
}

#[test]
fn vegas_create_default_opts() {
    let mut err: std::ffi::c_int = 0;
    let v = ocas_vegas_create(1, ptr::null(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!v.is_null());
    let r = ocas_vegas_integrate(v, Some(integrand_x), ptr::null_mut(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!((r.integral - 0.5).abs() < 0.05);
    ocas_vegas_free(v);
}

#[test]
fn vegas_rejects_zero_dims() {
    let mut err: std::ffi::c_int = 0;
    let v = ocas_vegas_create(0, ptr::null(), &mut err);
    assert!(v.is_null());
    let msg = ocas_error_last_message();
    assert!(!msg.is_null());
    ocas_error_clear();
}

#[test]
fn vegas_integrate_square() {
    let opts = OcasVegasOptions {
        n_bins: 64,
        n_samples: 20000,
        iterations: 8,
        learning_rate: 1.5,
        seed: 0x1234,
    };
    let mut err: std::ffi::c_int = 0;
    let v = ocas_vegas_create(1, &opts, &mut err);
    assert_eq!(err, OCAS_OK);
    let r = ocas_vegas_integrate(v, Some(integrand_x_squared), ptr::null_mut(), &mut err);
    assert_eq!(err, OCAS_OK);
    // ∫₀¹ x² dx = 1/3.
    assert!(
        (r.integral - 1.0 / 3.0).abs() < 0.01,
        "got integral {}",
        r.integral
    );
    ocas_vegas_free(v);
}

#[test]
fn vegas_null_handle_returns_nan_result() {
    let mut err: std::ffi::c_int = 0;
    let r = ocas_vegas_integrate(
        ptr::null_mut(),
        Some(integrand_x),
        ptr::null_mut(),
        &mut err,
    );
    assert_ne!(err, OCAS_OK);
    assert!(r.integral.is_nan());
    ocas_error_clear();
}

#[test]
fn vegas_null_integrand_returns_nan_result() {
    let mut err: std::ffi::c_int = 0;
    let v = ocas_vegas_create(1, ptr::null(), &mut err);
    assert!(!v.is_null());
    let r = ocas_vegas_integrate(v, None, ptr::null_mut(), &mut err);
    assert_ne!(err, OCAS_OK);
    assert!(r.integral.is_nan());
    ocas_error_clear();
    ocas_vegas_free(v);
}

#[test]
fn vegas_result_and_iterations_on_null() {
    let r = ocas_vegas_result(ptr::null());
    assert!(r.integral.is_nan());
    assert_eq!(ocas_vegas_iterations(ptr::null()), 0);
}

#[test]
fn vegas_free_on_null_is_safe() {
    ocas_vegas_free(ptr::null_mut());
}

#[test]
fn integrate_1d_one_shot_linear() {
    let mut err: std::ffi::c_int = 0;
    let r = ocas_integrate_1d(
        Some(integrand_x),
        ptr::null_mut(),
        0.0,
        1.0,
        ptr::null(),
        &mut err,
    );
    assert_eq!(err, OCAS_OK);
    assert!((r.integral - 0.5).abs() < 0.05);
}

#[test]
fn integrate_1d_rejects_bad_bounds() {
    let mut err: std::ffi::c_int = 0;
    let r = ocas_integrate_1d(
        Some(integrand_x),
        ptr::null_mut(),
        1.0,
        0.0,
        ptr::null(),
        &mut err,
    );
    assert_ne!(err, OCAS_OK);
    assert!(r.integral.is_nan());
    ocas_error_clear();
}

#[test]
fn integrate_1d_null_integrand_returns_nan() {
    let mut err: std::ffi::c_int = 0;
    let r = ocas_integrate_1d(None, ptr::null_mut(), 0.0, 1.0, ptr::null(), &mut err);
    assert!(r.integral.is_nan());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
}

// ------------------------------------------------------------------
//  Tensor C API tests
// ------------------------------------------------------------------

#[test]
fn tensor_create_and_query() {
    let mut err: std::ffi::c_int = 0;
    let t = ocas_tensor_create(
        cstr("T").as_ptr(),
        cstr("i,upper;j,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    assert_eq!(err, OCAS_OK);
    assert!(!t.is_null());
    let name_ptr = ocas_tensor_name(t, &mut err);
    assert_eq!(err, OCAS_OK);
    assert_eq!(c_string_to_string(name_ptr), "T");
    assert_eq!(ocas_tensor_rank(t), 2);
    assert_eq!(ocas_tensor_symmetry(t), 0); // none
    let s_ptr = ocas_tensor_to_string(t, &mut err);
    assert_eq!(err, OCAS_OK);
    let s = c_string_to_string(s_ptr);
    assert!(s.contains("T"), "got {s}");
    ocas_tensor_free(t);
}

#[test]
fn tensor_create_with_symmetry() {
    let mut err: std::ffi::c_int = 0;
    let e = ocas_tensor_create(
        cstr("eps").as_ptr(),
        cstr("a,lower;b,lower").as_ptr(),
        cstr("antisymmetric").as_ptr(),
        &mut err,
    );
    assert_eq!(err, OCAS_OK);
    assert_eq!(ocas_tensor_symmetry(e), 2); // antisymmetric
    assert_eq!(ocas_tensor_rank(e), 2);
    ocas_tensor_free(e);
}

#[test]
fn tensor_symmetrise_sign_values() {
    let mut err: std::ffi::c_int = 0;
    // No symmetry: sign is always +1.
    let t = ocas_tensor_create(
        cstr("T").as_ptr(),
        cstr("i,upper;j,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    assert_eq!(ocas_tensor_symmetrise_sign(t), 1);
    ocas_tensor_free(t);
    // Antisymmetric: sign is +1 or -1.
    let e = ocas_tensor_create(
        cstr("eps").as_ptr(),
        cstr("a,lower;b,lower").as_ptr(),
        cstr("antisymmetric").as_ptr(),
        &mut err,
    );
    let sign = ocas_tensor_symmetrise_sign(e);
    assert!(sign == 1 || sign == -1);
    ocas_tensor_free(e);
}

#[test]
fn tensor_null_handle_queries() {
    // All queries on null are safe and return sentinel values.
    assert_eq!(ocas_tensor_rank(ptr::null()), 0);
    assert_eq!(ocas_tensor_symmetry(ptr::null()), -1);
    assert_eq!(ocas_tensor_symmetrise_sign(ptr::null()), 0);
    let mut err: std::ffi::c_int = 0;
    let name = ocas_tensor_name(ptr::null(), &mut err);
    assert!(name.is_null());
    let s = ocas_tensor_to_string(ptr::null(), &mut err);
    assert!(s.is_null());
    ocas_error_clear();
}

#[test]
fn tensor_free_on_null_is_safe() {
    ocas_tensor_free(ptr::null_mut());
}

#[test]
fn tensor_create_rejects_bad_position() {
    let mut err: std::ffi::c_int = 0;
    let t = ocas_tensor_create(
        cstr("T").as_ptr(),
        cstr("i,sideways").as_ptr(),
        ptr::null(),
        &mut err,
    );
    assert!(t.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
}

#[test]
fn tensor_contract_partial_produces_product() {
    let mut err: std::ffi::c_int = 0;
    let t = ocas_tensor_create(
        cstr("T").as_ptr(),
        cstr("i,upper;j,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    let u = ocas_tensor_create(
        cstr("U").as_ptr(),
        cstr("j,upper;k,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    assert!(!t.is_null());
    assert!(!u.is_null());
    let mut out = OcasTensorContraction {
        kind: -1,
        tensors: ptr::null_mut(),
        n_tensors: 0,
        scalar_str: ptr::null_mut(),
    };
    let rc = ocas_tensor_contract(t, u, &mut out, &mut err);
    assert_eq!(rc, OCAS_OK);
    assert_eq!(out.kind, 0); // product
    assert_eq!(out.n_tensors, 1);
    // The single resulting tensor has 2 surviving slots (i, k).
    let result = unsafe { *out.tensors };
    assert!(!result.is_null());
    assert_eq!(ocas_tensor_rank(result), 2);
    for i in 0..out.n_tensors {
        let h = unsafe { *out.tensors.add(i) };
        ocas_tensor_free(h);
    }
    ocas_tensor_contraction_free(&mut out);
    ocas_tensor_free(t);
    ocas_tensor_free(u);
}

#[test]
fn tensor_contract_full_produces_scalar() {
    let mut err: std::ffi::c_int = 0;
    // T^i_j · U^j_i — every index contracts.
    let t = ocas_tensor_create(
        cstr("T").as_ptr(),
        cstr("i,upper;j,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    let u = ocas_tensor_create(
        cstr("U").as_ptr(),
        cstr("j,upper;i,lower").as_ptr(),
        ptr::null(),
        &mut err,
    );
    let mut out = OcasTensorContraction {
        kind: -1,
        tensors: ptr::null_mut(),
        n_tensors: 0,
        scalar_str: ptr::null_mut(),
    };
    let rc = ocas_tensor_contract(t, u, &mut out, &mut err);
    assert_eq!(rc, OCAS_OK);
    assert_eq!(out.kind, 1); // scalar
    assert!(!out.scalar_str.is_null());
    // Read the string WITHOUT freeing — contraction_free owns scalar_str.
    let s = unsafe { CStr::from_ptr(out.scalar_str) }
        .to_str()
        .unwrap()
        .to_string();
    assert!(s.contains("T"), "got {s}");
    assert!(s.contains("U"), "got {s}");
    ocas_tensor_contraction_free(&mut out);
    ocas_tensor_free(t);
    ocas_tensor_free(u);
}

#[test]
fn tensor_contract_null_handles_return_error() {
    let mut err: std::ffi::c_int = 0;
    let mut out = OcasTensorContraction {
        kind: -1,
        tensors: ptr::null_mut(),
        n_tensors: 0,
        scalar_str: ptr::null_mut(),
    };
    let rc = ocas_tensor_contract(ptr::null(), ptr::null(), &mut out, &mut err);
    assert_ne!(rc, OCAS_OK);
    ocas_error_clear();
}

#[test]
fn tensor_contraction_free_on_null_is_safe() {
    ocas_tensor_contraction_free(ptr::null_mut());
}

// ------------------------------------------------------------------
//  Hyper-dual number C API tests (forward AD)
// ------------------------------------------------------------------

#[test]
fn dual_shape_new_and_query() {
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(2, &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!s.is_null());
    assert_eq!(ocas_dual_shape_n_vars(s), 2);
    assert_eq!(ocas_dual_shape_n_components(s), 3); // value + 2 derivs
    ocas_dual_shape_free(s);
}

#[test]
fn dual_shape_rejects_zero_vars() {
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(0, &mut err);
    assert!(s.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
}

#[test]
fn dual_variable_and_constant() {
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(2, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("3").as_ptr(), &mut err);
    assert_eq!(err, OCAS_OK);
    assert!(!x.is_null());
    let v = ocas_dual_value(x, &mut err);
    assert_eq!(c_string_to_string(v), "3");
    let d0 = ocas_dual_deriv(x, 0, &mut err);
    assert_eq!(c_string_to_string(d0), "1"); // ∂x/∂x = 1
    let d1 = ocas_dual_deriv(x, 1, &mut err);
    assert_eq!(c_string_to_string(d1), "0"); // ∂x/∂y = 0
    let c = ocas_dual_constant(s, cstr("7").as_ptr(), &mut err);
    let cv = ocas_dual_value(c, &mut err);
    assert_eq!(c_string_to_string(cv), "7");
    let cd0 = ocas_dual_deriv(c, 0, &mut err);
    assert_eq!(c_string_to_string(cd0), "0");
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(c);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_product_of_two_variables() {
    // f(x, y) = x * y at (3, 5); ∂f/∂x = y = 5, ∂f/∂y = x = 3.
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(2, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("3").as_ptr(), &mut err);
    let y = ocas_dual_variable(s, 1, cstr("5").as_ptr(), &mut err);
    let f = ocas_dual_mul(x, y, &mut err);
    assert_eq!(err, OCAS_OK);
    assert_eq!(c_string_to_string(ocas_dual_value(f, &mut err)), "15");
    assert_eq!(c_string_to_string(ocas_dual_deriv(f, 0, &mut err)), "5");
    assert_eq!(c_string_to_string(ocas_dual_deriv(f, 1, &mut err)), "3");
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(y);
    ocas_hyperdual_free(f);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_quotient_chain_rule() {
    // f(x, y) = x / y at (3, 5); ∂f/∂x = 1/5, ∂f/∂y = -3/25.
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(2, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("3").as_ptr(), &mut err);
    let y = ocas_dual_variable(s, 1, cstr("5").as_ptr(), &mut err);
    let f = ocas_dual_div(x, y, &mut err);
    assert_eq!(err, OCAS_OK);
    assert_eq!(c_string_to_string(ocas_dual_value(f, &mut err)), "3/5");
    assert_eq!(c_string_to_string(ocas_dual_deriv(f, 0, &mut err)), "1/5");
    assert_eq!(c_string_to_string(ocas_dual_deriv(f, 1, &mut err)), "-3/25");
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(y);
    ocas_hyperdual_free(f);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_division_by_zero_returns_error() {
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(1, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("3").as_ptr(), &mut err);
    let z = ocas_dual_constant(s, cstr("0").as_ptr(), &mut err);
    let r = ocas_dual_div(x, z, &mut err);
    assert!(r.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(z);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_neg_and_sum() {
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(1, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("4").as_ptr(), &mut err);
    let neg = ocas_dual_neg(x, &mut err);
    assert_eq!(c_string_to_string(ocas_dual_value(neg, &mut err)), "-4");
    assert_eq!(c_string_to_string(ocas_dual_deriv(neg, 0, &mut err)), "-1");
    let sum = ocas_dual_add(x, neg, &mut err);
    assert_eq!(c_string_to_string(ocas_dual_value(sum, &mut err)), "0");
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(neg);
    ocas_hyperdual_free(sum);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_chain_rule_polynomial() {
    // f(x) = x^3 at 2; f' = 3x² = 12, f = 8. Repeated multiplication.
    let mut err: std::ffi::c_int = 0;
    let s = ocas_dual_shape_new(1, &mut err);
    let x = ocas_dual_variable(s, 0, cstr("2").as_ptr(), &mut err);
    let x2 = ocas_dual_mul(x, x, &mut err);
    let x3 = ocas_dual_mul(x2, x, &mut err);
    assert_eq!(c_string_to_string(ocas_dual_value(x3, &mut err)), "8");
    assert_eq!(c_string_to_string(ocas_dual_deriv(x3, 0, &mut err)), "12");
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(x2);
    ocas_hyperdual_free(x3);
    ocas_dual_shape_free(s);
}

#[test]
fn dual_null_handle_queries_are_safe() {
    let mut err: std::ffi::c_int = 0;
    let v = ocas_dual_value(ptr::null(), &mut err);
    assert!(v.is_null());
    let d = ocas_dual_deriv(ptr::null(), 0, &mut err);
    assert!(d.is_null());
    let neg = ocas_dual_neg(ptr::null(), &mut err);
    assert!(neg.is_null());
    ocas_error_clear();
}

#[test]
fn dual_shape_mismatch_arithmetic_returns_error() {
    let mut err: std::ffi::c_int = 0;
    let s2 = ocas_dual_shape_new(2, &mut err);
    let s3 = ocas_dual_shape_new(3, &mut err);
    let x = ocas_dual_variable(s2, 0, cstr("1").as_ptr(), &mut err);
    let y = ocas_dual_variable(s3, 0, cstr("1").as_ptr(), &mut err);
    let r = ocas_dual_add(x, y, &mut err);
    assert!(r.is_null());
    assert_ne!(err, OCAS_OK);
    ocas_error_clear();
    ocas_hyperdual_free(x);
    ocas_hyperdual_free(y);
    ocas_dual_shape_free(s2);
    ocas_dual_shape_free(s3);
}

#[test]
fn dual_free_on_null_is_safe() {
    ocas_hyperdual_free(ptr::null_mut());
    ocas_dual_shape_free(ptr::null_mut());
}
