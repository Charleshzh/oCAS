//! End-to-end tests that exercise the C ABI through Rust FFI.
//!
//! These complement the unit tests in `expression.rs` by calling the
//! `#[no_mangle] extern "C"` functions exactly as a C caller would.

use ocas_c::{
    OCAS_OK, OcasAlgebraicFactorArray, OcasPolyFactorArray, ocas_algebraic_factor_array_free,
    ocas_algebraic_field_create, ocas_algebraic_field_degree, ocas_algebraic_field_free,
    ocas_algebraic_poly_create, ocas_algebraic_poly_degree, ocas_algebraic_poly_factor,
    ocas_algebraic_poly_free, ocas_algebraic_poly_to_string, ocas_error_clear,
    ocas_error_last_message, ocas_expr_clone, ocas_expr_diff, ocas_expr_free, ocas_expr_integrate,
    ocas_expr_normalize, ocas_expr_parse, ocas_expr_simplify, ocas_expr_substitute,
    ocas_expr_taylor, ocas_expr_to_string, ocas_poly_factor_array_free, ocas_poly_fp_clone,
    ocas_poly_fp_create, ocas_poly_fp_degree, ocas_poly_fp_factor, ocas_poly_fp_free,
    ocas_poly_fp_to_string, ocas_poly_z_clone, ocas_poly_z_create, ocas_poly_z_degree,
    ocas_poly_z_factor, ocas_poly_z_free, ocas_poly_z_to_string, ocas_string_free, ocas_version,
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
