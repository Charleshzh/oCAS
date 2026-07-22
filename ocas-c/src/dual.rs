//! C/C++ bindings for hyper-dual numbers (forward automatic differentiation).
//!
//! Exposes [`ocas_domain::dual`] restricted to [`Rational`](ocas_domain::Rational)
//! coefficients. Only polynomial/rational arithmetic is supported.
//!
//! # Coefficient string format
//!
//! Rational coefficients are passed as strings `"num"` (denominator 1) or
//! `"num/den"`. Returned strings use the same format (e.g. `"5"`, `"1/2"`).

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{CStr, CString, c_char, c_int};
use std::ptr;
use std::sync::Arc;

use ocas_domain::Rational;
use ocas_domain::dual::{DualShape, HyperDual, new_first_order};

use crate::error::{
    OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER, OCAS_ERROR_RUNTIME, OCAS_OK, set,
};

/// Owned storage for a shape handle. Wrapping the `Arc` in a named struct
/// avoids clippy's `Box<Arc<_>>` lint while keeping a stable opaque layout.
struct DualShapeStore {
    shape: Arc<DualShape>,
}

/// Opaque handle for a dual-number shape.
#[repr(C)]
pub struct OcasDualShape {
    _private: [u8; 0],
}

/// Opaque handle for a hyper-dual number.
#[repr(C)]
pub struct OcasHyperDual {
    _private: [u8; 0],
}

/// Storage behind an [`OcasHyperDual`]: the dual value plus the shape it
/// shares (for shape-equality checks during arithmetic).
struct DualStore {
    value: HyperDual<Rational>,
    shape: Arc<DualShape>,
}

// ------------------------------------------------------------------
//  Opaque-handle helpers
// ------------------------------------------------------------------

fn shape_ptr(s: Box<DualShapeStore>) -> *mut OcasDualShape {
    Box::into_raw(s) as *mut OcasDualShape
}

fn shape_ref<'a>(s: *const OcasDualShape) -> Option<&'a Arc<DualShape>> {
    if s.is_null() {
        return None;
    }
    // SAFETY: valid while the handle is alive.
    Some(unsafe { &(*(s as *const DualShapeStore)).shape })
}

fn dual_ptr(d: Box<DualStore>) -> *mut OcasHyperDual {
    Box::into_raw(d) as *mut OcasHyperDual
}

fn dual_ref<'a>(d: *const OcasHyperDual) -> Option<&'a DualStore> {
    if d.is_null() {
        return None;
    }
    Some(unsafe { &*(d as *const DualStore) })
}

/// Parse a coefficient string `"num"` or `"num/den"` into a [`Rational`].
fn parse_rational_str(s: &str) -> Result<Rational, String> {
    let s = s.trim();
    if let Some((n_str, d_str)) = s.split_once('/') {
        let n: i64 = n_str
            .trim()
            .parse()
            .map_err(|_| format!("invalid numerator {n_str:?}"))?;
        let d: i64 = d_str
            .trim()
            .parse()
            .map_err(|_| format!("invalid denominator {d_str:?}"))?;
        if d == 0 {
            return Err("denominator cannot be zero".to_string());
        }
        Ok(Rational::new(n, d))
    } else {
        let n: i64 = s.parse().map_err(|_| format!("invalid integer {s:?}"))?;
        Ok(Rational::new(n, 1))
    }
}

/// Format a [`Rational`] as `"num/den"` (or `"num"` when the denominator is 1).
fn rational_to_string(r: &Rational) -> String {
    let n = r.numer().to_i64().unwrap_or(0);
    let d = r.denom().to_i64().unwrap_or(0);
    if d == 1 {
        n.to_string()
    } else {
        format!("{n}/{d}")
    }
}

fn build_dual(value: HyperDual<Rational>) -> Box<DualStore> {
    let shape = value.shape().clone();
    Box::new(DualStore { value, shape })
}

// ------------------------------------------------------------------
//  DualShape C API
// ------------------------------------------------------------------

/// Build a first-order shape tracking one derivative per variable for
/// `n_vars` variables. Returns an opaque handle or `NULL` on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_shape_new(n_vars: usize, err: *mut c_int) -> *mut OcasDualShape {
    crate::error::clear();
    if n_vars == 0 {
        set(OCAS_ERROR_INVALID_ARGUMENT, "n_vars must be >= 1");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let shape = new_first_order::<Rational>(n_vars);
    crate::error::write_last_code(err);
    shape_ptr(Box::new(DualShapeStore { shape }))
}

/// Free a dual-number shape handle. Safe to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_shape_free(s: *mut OcasDualShape) {
    if s.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(s as *mut DualShapeStore));
    }
}

/// Return the number of differentiation variables, or `0` on a null handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_shape_n_vars(s: *const OcasDualShape) -> usize {
    match shape_ref(s) {
        Some(sh) => sh.n_vars(),
        None => 0,
    }
}

/// Return the total number of components, or `0` on a null handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_shape_n_components(s: *const OcasDualShape) -> usize {
    match shape_ref(s) {
        Some(sh) => sh.n_components(),
        None => 0,
    }
}

// ------------------------------------------------------------------
//  HyperDual C API
// ------------------------------------------------------------------

/// Create an independent variable `x_i = value` (coefficient string
/// `"num"` or `"num/den"`). Derivative is 1 w.r.t. variable `i`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_variable(
    shape: *const OcasDualShape,
    i: usize,
    coeff: *const c_char,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    crate::error::clear();
    let sh = match shape_ref(shape) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "shape handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    if i >= sh.n_vars() {
        set(OCAS_ERROR_INVALID_ARGUMENT, "variable index out of range");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    if coeff.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "coefficient string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let s = unsafe { CStr::from_ptr(coeff) };
    let s = match s.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(
                OCAS_ERROR_INVALID_ARGUMENT,
                "coefficient string is not valid UTF-8",
            );
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match parse_rational_str(s) {
        Ok(v) => {
            let hd = HyperDual::variable(sh, i, v);
            crate::error::write_last_code(err);
            dual_ptr(build_dual(hd))
        }
        Err(msg) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, &msg);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Create a constant dual number (all derivatives zero).
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_constant(
    shape: *const OcasDualShape,
    coeff: *const c_char,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    crate::error::clear();
    let sh = match shape_ref(shape) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "shape handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    if coeff.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "coefficient string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let s = unsafe { CStr::from_ptr(coeff) };
    let s = match s.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(
                OCAS_ERROR_INVALID_ARGUMENT,
                "coefficient string is not valid UTF-8",
            );
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match parse_rational_str(s) {
        Ok(v) => {
            let hd = HyperDual::constant(sh, v);
            crate::error::write_last_code(err);
            dual_ptr(build_dual(hd))
        }
        Err(msg) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, &msg);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Free a hyper-dual handle. Safe to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_hyperdual_free(d: *mut OcasHyperDual) {
    if d.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(d as *mut DualStore));
    }
}

/// Return the scalar value component as a heap-allocated string
/// (`"num"` or `"num/den"`). The caller must release it with
/// [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_value(d: *const OcasHyperDual, err: *mut c_int) -> *mut c_char {
    crate::error::clear();
    match dual_ref(d) {
        Some(store) => {
            let s = rational_to_string(store.value.value());
            match CString::new(s) {
                Ok(cs) => {
                    crate::error::write_last_code(err);
                    cs.into_raw()
                }
                Err(_) => {
                    set(OCAS_ERROR_RUNTIME, "value string contains a NUL byte");
                    crate::error::write_last_code(err);
                    ptr::null_mut()
                }
            }
        }
        None => {
            set(OCAS_ERROR_NULL_POINTER, "hyper-dual handle is null");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Return the derivative w.r.t. variable `i` as a heap-allocated string, or
/// `NULL` if the shape has no first-order component for `i`. On a null
/// handle or other error the error code is set and `NULL` is returned.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_deriv(
    d: *const OcasHyperDual,
    i: usize,
    err: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let store = match dual_ref(d) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "hyper-dual handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match store.value.deriv(i) {
        Some(r) => {
            let s = rational_to_string(r);
            match CString::new(s) {
                Ok(cs) => {
                    crate::error::write_last_code(err);
                    cs.into_raw()
                }
                Err(_) => {
                    set(OCAS_ERROR_RUNTIME, "deriv string contains a NUL byte");
                    crate::error::write_last_code(err);
                    ptr::null_mut()
                }
            }
        }
        None => {
            // No derivative component: not an error, just no result.
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Compute `a + b` and return a new handle. Both operands must share a shape.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_add(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    binary_op(a, b, err, |x, y| x + y)
}

/// Compute `a - b` and return a new handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_sub(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    binary_op(a, b, err, |x, y| x - y)
}

/// Compute `a * b` and return a new handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_mul(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    binary_op(a, b, err, |x, y| x * y)
}

/// Compute `a / b` and return a new handle. Sets an error if the divisor's
/// value component is zero (which would otherwise panic).
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_div(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
) -> *mut OcasHyperDual {
    crate::error::clear();
    let (av, bv) = match fetch_pair(a, b, err) {
        Some(pair) => pair,
        None => return ptr::null_mut(),
    };
    if bv.value() == &Rational::new(0, 1) {
        set(
            OCAS_ERROR_INVALID_ARGUMENT,
            "division by zero (value component is zero)",
        );
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let result = av / bv;
    crate::error::write_last_code(err);
    dual_ptr(build_dual(result))
}

/// Compute `-a` (negation) and return a new handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_dual_neg(a: *const OcasHyperDual, err: *mut c_int) -> *mut OcasHyperDual {
    crate::error::clear();
    let store = match dual_ref(a) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "hyper-dual handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let result = -store.value.clone();
    crate::error::write_last_code(err);
    dual_ptr(build_dual(result))
}

// ------------------------------------------------------------------
//  Internal binary-op helper
// ------------------------------------------------------------------

/// Fetch the two operand values, validating both handles and shape equality.
fn fetch_pair(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
) -> Option<(HyperDual<Rational>, HyperDual<Rational>)> {
    crate::error::clear();
    let sa = match dual_ref(a) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "first hyper-dual handle is null");
            crate::error::write_last_code(err);
            return None;
        }
    };
    let sb = match dual_ref(b) {
        Some(s) => s,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "second hyper-dual handle is null");
            crate::error::write_last_code(err);
            return None;
        }
    };
    if !Arc::ptr_eq(&sa.shape, &sb.shape) {
        set(
            OCAS_ERROR_INVALID_ARGUMENT,
            "operands have different shapes",
        );
        crate::error::write_last_code(err);
        return None;
    }
    Some((sa.value.clone(), sb.value.clone()))
}

fn binary_op<F>(
    a: *const OcasHyperDual,
    b: *const OcasHyperDual,
    err: *mut c_int,
    op: F,
) -> *mut OcasHyperDual
where
    F: FnOnce(HyperDual<Rational>, HyperDual<Rational>) -> HyperDual<Rational>,
{
    match fetch_pair(a, b, err) {
        Some((av, bv)) => {
            let result = op(av, bv);
            crate::error::write_last_code(err);
            dual_ptr(build_dual(result))
        }
        None => ptr::null_mut(),
    }
}

#[allow(unused_imports)]
use OCAS_OK as _OcasOkImport;
