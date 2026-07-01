//! C/C++ bindings for oCAS.
//!
//! This crate exposes a stable C ABI for the oCAS computer algebra system.
//! It is organized into two modules:
//!
//! - [`error`]: thread-local error reporting with integer error codes.
//! - [`expression`]: expression lifecycle (parse / free / clone / to_string)
//!   and operations (differentiation, integration, series, simplification,
//!   substitution).
//!
//! ## FFI conventions
//!
//! - Opaque pointers (`OcasExpr`, `OcasArena`) for Rust-owned objects.
//! - Thread-local error state for recoverable failures; query with
//!   [`ocas_error_last_message`] and clear with [`ocas_error_clear`].
//! - Strings returned by [`ocas_expr_to_string`] are heap-allocated and
//!   owned by the caller; release them with [`ocas_string_free`].
//!
//! ## Example
//!
//! ```c
//! #include <ocas.h>
//! #include <stdio.h>
//!
//! int main(void) {
//!     int err = 0;
//!     OcasExpr *e = ocas_expr_parse("x^2", &err);
//!     if (e == NULL) {
//!         fprintf(stderr, "parse failed: %s\n", ocas_error_last_message());
//!         return 1;
//!     }
//!     OcasExpr *d = ocas_expr_diff(e, "x", &err);
//!     char *s = ocas_expr_to_string(d, &err);
//!     printf("d/dx(x^2) = %s\n", s);
//!     ocas_string_free(s);
//!     ocas_expr_free(d);
//!     ocas_expr_free(e);
//!     return 0;
//! }
//! ```

#![warn(missing_docs)]
#![allow(clippy::missing_const_for_thread_local)]

pub mod error;
pub mod expression;

use std::ffi::c_char;

// Re-export error codes and the expression C API at the crate root so both
// Rust callers (integration tests) and cbindgen can find them easily.
pub use error::{
    OCAS_ERROR_DIVISION_BY_ZERO, OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER,
    OCAS_ERROR_OUT_OF_MEMORY, OCAS_ERROR_PARSE, OCAS_ERROR_RUNTIME, OCAS_OK,
};
pub use expression::{
    OcasExpr, ocas_expr_clone, ocas_expr_diff, ocas_expr_free, ocas_expr_integrate,
    ocas_expr_normalize, ocas_expr_parse, ocas_expr_simplify, ocas_expr_substitute,
    ocas_expr_taylor, ocas_expr_to_string, ocas_string_free,
};

/// Return the oCAS version string.
///
/// The returned pointer is valid for the lifetime of the program and must
/// not be freed by the caller.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0")
        .as_ptr()
        .cast::<c_char>()
}

/// Return the message for the last error on the calling thread, or `NULL`
/// if no error has occurred.
///
/// The returned string is owned by the library and must not be freed or
/// modified by the caller. It remains valid until the next call that sets
/// an error on the same thread or until [`ocas_error_clear`] is called.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_error_last_message() -> *const c_char {
    error::last_message_ptr()
}

/// Clear the last error on the calling thread.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_error_clear() {
    error::clear();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::{CStr, CString};

    #[test]
    fn version_returns_workspace_version() {
        let ptr = ocas_version();
        assert!(!ptr.is_null());
        let cstr = unsafe { CStr::from_ptr(ptr) };
        assert_eq!(cstr.to_str().unwrap(), env!("CARGO_PKG_VERSION"));
    }

    #[test]
    fn error_clear_works() {
        // Trigger an error via the expression API, then clear it.
        let bad = CString::new("@@@").unwrap();
        let mut err: std::ffi::c_int = 0;
        let result = unsafe { ocas_expr_parse(bad.as_ptr(), &mut err) };
        assert!(result.is_null());
        assert!(!ocas_error_last_message().is_null());
        ocas_error_clear();
        assert!(ocas_error_last_message().is_null());
    }
}
