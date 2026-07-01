//! C/C++ bindings for oCAS.
//!
//! This crate exposes a minimal, stable C ABI for the 0.1.0 release. It
//! demonstrates the FFI glue conventions described in ARCHITECTURE.md:
//!
//! - Opaque pointers for Rust-owned objects.
//! - Thread-local error state for recoverable failures.
//! - C strings returned as `const char*` that remain valid until the next
//!   call on the same thread or until `ocas_error_clear` is called.
//!
//! Currently only runtime objects (version, arena) are exposed. Expression
//! operations will be added in later releases.

#![allow(clippy::missing_const_for_thread_local)]

use std::cell::RefCell;
use std::ffi::{CString, c_char, c_int};
use std::ptr;

use ocas_core::arena::Arena;

/// Success error code.
pub const OCAS_OK: c_int = 0;
/// A null pointer was passed where a non-null pointer was required.
pub const OCAS_ERROR_NULL_POINTER: c_int = 1;
/// An operation failed inside the oCAS runtime.
pub const OCAS_ERROR_RUNTIME: c_int = 2;

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
}

fn set_error(message: &str) {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = CString::new(message).ok();
    });
}

fn clear_error() {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = None;
    });
}

fn error_message_ptr() -> *const c_char {
    LAST_ERROR.with(|e| {
        e.borrow()
            .as_ref()
            .map(|s| s.as_ptr())
            .unwrap_or(ptr::null())
    })
}

/// Opaque arena handle.
#[repr(C)]
pub struct OcasArena {
    _private: [u8; 0],
}

/// Return the oCAS version string.
///
/// The returned pointer is valid for the lifetime of the program and must not
/// be freed by the caller.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_version() -> *const c_char {
    concat!("0.1.0", "\0").as_ptr().cast::<c_char>()
}

/// Return the message for the last error on the calling thread, or `NULL` if
/// no error has occurred.
///
/// The returned string is owned by the library and must not be freed or
/// modified by the caller. It remains valid until the next call that sets an
/// error on the same thread or until `ocas_error_clear` is called.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_error_last_message() -> *const c_char {
    error_message_ptr()
}

/// Clear the last error on the calling thread.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_error_clear() {
    clear_error();
}

/// Create a new arena and return an opaque pointer to it.
///
/// Returns `NULL` if allocation fails. Use `ocas_error_last_message` to
/// retrieve the error message.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_arena_new() -> *mut OcasArena {
    clear_error();
    let arena: Box<Arena> = match std::panic::catch_unwind(|| Box::new(Arena::new())).ok() {
        Some(a) => a,
        None => {
            set_error("failed to allocate arena");
            return ptr::null_mut();
        }
    };
    let raw = Box::into_raw(arena);
    raw.cast::<OcasArena>()
}

/// Free an arena previously created with `ocas_arena_new`.
///
/// Passing `NULL` is a no-op.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_arena_free(arena: *mut OcasArena) {
    if arena.is_null() {
        return;
    }
    clear_error();
    // SAFETY: `arena` was created by `ocas_arena_new` and is not null.
    unsafe {
        let _ = Box::from_raw(arena.cast::<Arena>());
    }
}

/// Allocate a single `i64` in the arena and return its value.
///
/// This is a trivial demonstration of the arena lifetime model. Returns
/// `OCAS_ERROR_NULL_POINTER` if `arena` is null; otherwise returns `OCAS_OK`.
///
/// # Safety
///
/// `arena` must be a non-null pointer returned by `ocas_arena_new`. `out`
/// must be a valid, non-null pointer to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_arena_alloc_i64(
    arena: *mut OcasArena,
    value: i64,
    out: *mut i64,
) -> c_int {
    clear_error();
    if arena.is_null() {
        set_error("arena is null");
        return OCAS_ERROR_NULL_POINTER;
    }
    if out.is_null() {
        set_error("output pointer is null");
        return OCAS_ERROR_NULL_POINTER;
    }
    // SAFETY: `arena` is non-null and was created by `ocas_arena_new`.
    let arena_ref = unsafe { &*arena.cast::<Arena>() };
    let allocated = arena_ref.allocate_with(|| value);
    // SAFETY: `out` is non-null and points to writable memory.
    unsafe {
        ptr::write(out, *allocated);
    }
    OCAS_OK
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::ffi::CStr;

    #[test]
    fn version_returns_expected_string() {
        let ptr = ocas_version();
        assert!(!ptr.is_null());
        let cstr = unsafe { CStr::from_ptr(ptr) };
        assert_eq!(cstr.to_str().unwrap(), "0.1.0");
    }

    #[test]
    fn arena_lifecycle() {
        let arena = ocas_arena_new();
        assert!(!arena.is_null());
        ocas_arena_free(arena);
    }

    #[test]
    fn arena_null_pointer_errors() {
        let mut out = 0i64;
        let rc = unsafe { ocas_arena_alloc_i64(ptr::null_mut(), 42, &mut out) };
        assert_eq!(rc, OCAS_ERROR_NULL_POINTER);
        assert!(!ocas_error_last_message().is_null());
    }

    #[test]
    fn arena_alloc_i64_roundtrip() {
        let arena = ocas_arena_new();
        assert!(!arena.is_null());
        let mut out = 0i64;
        let rc = unsafe { ocas_arena_alloc_i64(arena, 123, &mut out) };
        assert_eq!(rc, OCAS_OK);
        assert_eq!(out, 123);
        ocas_arena_free(arena);
    }

    #[test]
    fn error_clear_works() {
        unsafe { ocas_arena_alloc_i64(ptr::null_mut(), 0, ptr::null_mut()) };
        assert!(!ocas_error_last_message().is_null());
        ocas_error_clear();
        assert!(ocas_error_last_message().is_null());
    }
}
