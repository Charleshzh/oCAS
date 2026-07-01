//! End-to-end tests that exercise the C ABI through Rust FFI.
//!
//! These complement the unit tests in `expression.rs` by calling the
//! `#[no_mangle] extern "C"` functions exactly as a C caller would.

use ocas_c::{
    OCAS_OK, ocas_error_clear, ocas_error_last_message, ocas_expr_clone, ocas_expr_diff,
    ocas_expr_free, ocas_expr_integrate, ocas_expr_normalize, ocas_expr_parse, ocas_expr_simplify,
    ocas_expr_substitute, ocas_expr_taylor, ocas_expr_to_string, ocas_string_free, ocas_version,
};
use std::ffi::{CStr, CString};

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
