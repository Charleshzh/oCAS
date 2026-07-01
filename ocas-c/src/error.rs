//! Error reporting for the C API.
//!
//! Errors are reported through a thread-local slot accessed via
//! [`last_message`](crate::error::last_message) and
//! [`clear`](crate::error::clear). Each C-callable function in this crate
//! clears the slot on entry and sets it on failure before returning a
//! non-zero error code.

use std::cell::RefCell;
use std::ffi::{CString, c_char, c_int};
use std::ptr;

/// Success error code.
pub const OCAS_OK: c_int = 0;
/// A null pointer was passed where a non-null pointer was required.
pub const OCAS_ERROR_NULL_POINTER: c_int = 1;
/// An operation failed inside the oCAS runtime.
pub const OCAS_ERROR_RUNTIME: c_int = 2;
/// A string could not be parsed as a valid expression.
pub const OCAS_ERROR_PARSE: c_int = 3;
/// An argument had an unexpected type or value.
pub const OCAS_ERROR_INVALID_ARGUMENT: c_int = 4;
/// A division by zero or similar undefined operation was attempted.
pub const OCAS_ERROR_DIVISION_BY_ZERO: c_int = 5;
/// Memory allocation failed.
pub const OCAS_ERROR_OUT_OF_MEMORY: c_int = 6;

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
    static LAST_ERROR_CODE: RefCell<c_int> = const { RefCell::new(OCAS_OK) };
}

/// Record an error on the calling thread.
pub fn set(code: c_int, message: &str) {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = CString::new(message).ok();
    });
    LAST_ERROR_CODE.with(|c| *c.borrow_mut() = code);
}

/// Record a runtime error with the message from an error display.
pub fn set_runtime<E: std::fmt::Display>(err: &E) -> c_int {
    set(OCAS_ERROR_RUNTIME, &err.to_string());
    OCAS_ERROR_RUNTIME
}

/// Clear the last error on the calling thread.
pub fn clear() {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = None;
    });
    LAST_ERROR_CODE.with(|c| *c.borrow_mut() = OCAS_OK);
}

/// Return the message pointer for the last error, or null if none.
pub fn last_message_ptr() -> *const c_char {
    LAST_ERROR.with(|e| {
        e.borrow()
            .as_ref()
            .map(|s| s.as_ptr())
            .unwrap_or(ptr::null())
    })
}

/// Return the last error code on the calling thread.
pub fn last_code() -> std::ffi::c_int {
    LAST_ERROR_CODE.with(|c| *c.borrow())
}

/// Write the last error code to `err_out` if non-null.
/// Safe to call with a null pointer (no-op).
#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub fn write_last_code(err_out: *mut std::ffi::c_int) {
    if !err_out.is_null() {
        // SAFETY: caller-provided pointer; we only write if non-null.
        // This matches the C convention where err_out is an out-parameter.
        unsafe { std::ptr::write(err_out, last_code()) };
    }
}
