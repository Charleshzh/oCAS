//! Expression lifecycle and operation C API.
//!
//! Each [`OcasExpr`] handle owns a private arena, an [`AtomArena`] borrowing
//! that arena, and a normalized root [`Atom`]. The arena is heap-allocated
//! and the handle stores a raw pointer plus references into that allocation,
//! so the handle is movable; the heap address never changes for the
//! lifetime of the handle.
//!
//! Strings returned from `ocas_expr_to_string` are heap-allocated,
//! null-terminated, and must be released with [`ocas_string_free`].

use std::ffi::{CStr, CString, c_char, c_int};
use std::ptr;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, Symbol};
use ocas_calc::{diff, integrate, substitute, taylor};
use ocas_core::arena::Arena;
use ocas_parse::parse;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use crate::error::{OCAS_ERROR_NULL_POINTER, OCAS_ERROR_PARSE, OCAS_ERROR_RUNTIME, set};

/// Opaque expression handle. Each handle owns its arena; use
/// [`ocas_expr_free`] to release it.
#[repr(C)]
pub struct OcasExpr {
    _private: [u8; 0],
}

/// Internal representation behind an [`OcasExpr`] handle. Owns a leaked
/// `Arena` and a leaked `AtomArena<'static>` borrowing that arena; both are
/// recovered on `Drop`.
struct ExprBox {
    /// Raw pointer to the heap-allocated arena. Recovered in `Drop`.
    arena_ptr: *mut Arena,
    /// Raw pointer to the heap-allocated `AtomArena<'static>`. Recovered in
    /// `Drop`. Must be dropped before `arena_ptr` because it borrows it.
    ctx_ptr: *mut AtomArena<'static>,
    /// Normalized root atom, borrowing `&*ctx_ptr`.
    atom: Atom<'static>,
}

// SAFETY: `ExprBox` owns two heap allocations (`arena_ptr`, `ctx_ptr`).
// References inside `ctx` and `atom` point into those allocations, which
// remain valid until `Drop`. The heap allocations are not tied to any
// thread.
unsafe impl Send for ExprBox {}

impl Drop for ExprBox {
    fn drop(&mut self) {
        // Drop order matters: `atom` (a field) is dropped first by the
        // compiler, then we free `ctx_ptr` (which borrows `arena_ptr`),
        // then `arena_ptr`.
        // SAFETY: both pointers were created by `Box::into_raw` and are
        // valid. `ctx_ptr` borrows `arena_ptr`, so it must drop first.
        unsafe {
            let _ = Box::from_raw(self.ctx_ptr);
            let _ = Box::from_raw(self.arena_ptr);
        }
    }
}

/// Extend a string slice's lifetime to `'static`.
///
/// # Safety
///
/// Safe in practice for our use because oCAS atoms never retain borrows of
/// the input string — the parser copies all characters into arena-owned
/// nodes. We only need the longer lifetime to satisfy `parse`'s signature.
unsafe fn extend_str_lifetime(s: &str) -> &'static str {
    unsafe { std::mem::transmute::<&str, &'static str>(s) }
}

/// Leak a fresh arena and a fresh `AtomArena<'static>` borrowing it. Returns
/// raw pointers suitable for storage in [`ExprBox`].
fn leak_arena_and_ctx() -> (*mut Arena, *mut AtomArena<'static>) {
    let arena_box: Box<Arena> = Box::new(Arena::new());
    let arena_ptr = Box::into_raw(arena_box);
    // SAFETY: `arena_ptr` is valid for as long as the caller keeps the box
    // alive. We extend the borrowed lifetime to `'static` to match the
    // intended ownership model of `ExprBox`.
    let arena_ref: &'static Arena = unsafe { &*arena_ptr };
    let ctx = AtomArena::new(arena_ref);
    let ctx_box: Box<AtomArena<'static>> = Box::new(ctx);
    let ctx_ptr = Box::into_raw(ctx_box);
    (arena_ptr, ctx_ptr)
}

/// Borrow both leaked pointers as `&'static AtomArena<'static>`.
///
/// # Safety
///
/// `ctx_ptr` must come from [`leak_arena_and_ctx`] and be valid.
unsafe fn static_ctx(ctx_ptr: *mut AtomArena<'static>) -> &'static AtomArena<'static> {
    unsafe { &*ctx_ptr }
}

/// Free a leaked arena and its companion `AtomArena`.
///
/// # Safety
///
/// Both pointers must come from [`leak_arena_and_ctx`] and not have been
/// freed. `ctx_ptr` is dropped first because it borrows `arena_ptr`.
unsafe fn free_leaked(arena_ptr: *mut Arena, ctx_ptr: *mut AtomArena<'static>) {
    unsafe {
        let _ = Box::from_raw(ctx_ptr);
        let _ = Box::from_raw(arena_ptr);
    }
}

/// RAII guard that frees the leaked arena pair unless explicitly disarmed.
///
/// Used inside [`ExprBox::build`] to guarantee cleanup on panic: if the
/// builder closure or `normalize` panics, the guard's `Drop` recovers the
/// allocations. On the success path, `disarm()` is called before
/// transferring ownership to `ExprBox`.
struct ArenaGuard {
    arena_ptr: *mut Arena,
    ctx_ptr: *mut AtomArena<'static>,
    armed: bool,
}

impl ArenaGuard {
    fn new(arena_ptr: *mut Arena, ctx_ptr: *mut AtomArena<'static>) -> Self {
        ArenaGuard {
            arena_ptr,
            ctx_ptr,
            armed: true,
        }
    }

    /// Mark the guard as disarmed; `Drop` will no longer free the arenas.
    /// Ownership has been transferred elsewhere.
    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for ArenaGuard {
    fn drop(&mut self) {
        if self.armed {
            // SAFETY: same contract as `free_leaked`.
            unsafe {
                free_leaked(self.arena_ptr, self.ctx_ptr);
            }
        }
    }
}

impl ExprBox {
    /// Build an `ExprBox` from a builder closure that receives a
    /// `&'static AtomArena<'static>` and returns an `Atom<'static>`. On
    /// failure (including panic) the leaked allocations are freed.
    fn build<F>(f: F) -> Result<Box<Self>, String>
    where
        F: FnOnce(&'static AtomArena<'static>) -> Result<Atom<'static>, String>,
    {
        let (arena_ptr, ctx_ptr) = leak_arena_and_ctx();
        let mut guard = ArenaGuard::new(arena_ptr, ctx_ptr);
        let ctx = unsafe { static_ctx(ctx_ptr) };
        // If `f` or `normalize` panics, `guard` is dropped and frees the
        // arenas. On success we disarm and transfer ownership to ExprBox.
        let atom = f(ctx)?;
        let normalized = normalize(ctx, atom);
        guard.disarm();
        Ok(Box::new(ExprBox {
            arena_ptr,
            ctx_ptr,
            atom: normalized,
        }))
    }

    /// Parse `input` into a new `ExprBox`.
    fn from_parse(input: &str) -> Result<Box<Self>, String> {
        let static_input = unsafe { extend_str_lifetime(input) };
        Self::build(|ctx| match parse(ctx, static_input) {
            Ok(a) => Ok(a),
            Err(e) => Err(e.to_string()),
        })
    }

    /// Rebuild from the string form of `src` into a fresh arena.
    fn from_expr_string(src: String) -> Result<Box<Self>, String> {
        let static_src = unsafe { extend_str_lifetime(&src) };
        Self::build(|ctx| match parse(ctx, static_src) {
            Ok(a) => Ok(a),
            Err(e) => Err(e.to_string()),
        })
    }

    /// Apply a calculus operation `op` taking `(ctx, atom, var)`.
    fn apply_calculus(
        &self,
        var: &str,
        op: fn(&'static AtomArena<'static>, Atom<'static>, Symbol) -> Atom<'static>,
    ) -> Result<Box<Self>, String> {
        let var_sym = Symbol::new(var);
        let src = self.atom.to_string();
        let static_src = unsafe { extend_str_lifetime(&src) };
        Self::build(|ctx| match parse(ctx, static_src) {
            Ok(a) => Ok(op(ctx, a, var_sym)),
            Err(e) => Err(e.to_string()),
        })
    }

    /// Apply a substitution `self[var -> replacement]`.
    fn apply_substitute(&self, var: &str, replacement: &ExprBox) -> Result<Box<Self>, String> {
        let var_sym = Symbol::new(var);
        let expr_src = self.atom.to_string();
        let repl_src = replacement.atom.to_string();
        let static_expr = unsafe { extend_str_lifetime(&expr_src) };
        let static_repl = unsafe { extend_str_lifetime(&repl_src) };
        Self::build(|ctx| {
            let e = parse(ctx, static_expr).map_err(|e| e.to_string())?;
            let r = parse(ctx, static_repl).map_err(|e| e.to_string())?;
            Ok(substitute(ctx, e, var_sym, r))
        })
    }

    /// Compute the Taylor series of `self` around `point` up to `order`.
    fn apply_taylor(&self, var: &str, point: &ExprBox, order: usize) -> Result<Box<Self>, String> {
        let var_sym = Symbol::new(var);
        let expr_src = self.atom.to_string();
        let point_src = point.atom.to_string();
        let static_expr = unsafe { extend_str_lifetime(&expr_src) };
        let static_point = unsafe { extend_str_lifetime(&point_src) };
        Self::build(|ctx| {
            let e = parse(ctx, static_expr).map_err(|e| e.to_string())?;
            let p = parse(ctx, static_point).map_err(|e| e.to_string())?;
            Ok(taylor(ctx, e, var_sym, p, order))
        })
    }

    /// Simplify using the default rule set.
    fn simplify_default(&self) -> Result<Box<Self>, String> {
        let src = self.atom.to_string();
        let static_src = unsafe { extend_str_lifetime(&src) };
        Self::build(|ctx| {
            let a = parse(ctx, static_src).map_err(|e| e.to_string())?;
            let rules = default_rules(ctx, &());
            Ok(simplify(ctx, a, &rules, 20))
        })
    }
}

/// Convert a C string pointer to a Rust `&str`, setting an error on failure.
fn cstr_to_str<'a>(s: *const c_char, what: &str) -> Option<&'a str> {
    if s.is_null() {
        set(OCAS_ERROR_NULL_POINTER, &format!("{what} is null"));
        return None;
    }
    // SAFETY: caller guarantees `s` is a valid null-terminated C string.
    let cstr = unsafe { CStr::from_ptr(s) };
    match cstr.to_str() {
        Ok(s) => Some(s),
        Err(_) => {
            set(
                crate::error::OCAS_ERROR_INVALID_ARGUMENT,
                &format!("{what} is not valid UTF-8"),
            );
            None
        }
    }
}

/// Unwrap an [`OcasExpr`] pointer, setting a null error if needed.
fn as_expr<'a>(handle: *const OcasExpr) -> Option<&'a ExprBox> {
    if handle.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "expression handle is null");
        return None;
    }
    // SAFETY: caller guarantees `handle` is valid and not freed.
    Some(unsafe { &*handle.cast::<ExprBox>() })
}

/// Convert a `catch_unwind` result into an opaque handle (or null).
fn finish_op(
    result: std::thread::Result<Result<Box<ExprBox>, String>>,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    let boxed = match result {
        Ok(Ok(b)) => b,
        Ok(Err(msg)) => {
            set(OCAS_ERROR_RUNTIME, &msg);
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            return ptr::null_mut();
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during operation");
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            return ptr::null_mut();
        }
    };
    if !err_out.is_null() {
        unsafe { ptr::write(err_out, crate::error::OCAS_OK) };
    }
    Box::into_raw(boxed).cast::<OcasExpr>()
}

// ---------------------------------------------------------------------------
// C-callable functions
// ---------------------------------------------------------------------------

/// Parse `input` into a new expression handle.
///
/// Returns `NULL` on failure; use [`ocas_error_last_message`](crate::ocas_error_last_message)
/// to retrieve the error message. On success, `*err_out` (if non-null) is
/// set to [`OCAS_OK`](crate::OCAS_OK).
///
/// # Safety
///
/// `input` must be a valid null-terminated C string. `err_out` if non-null
/// must point to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_parse(
    input: *const c_char,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(s) = cstr_to_str(input, "input") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let result = std::panic::catch_unwind(|| ExprBox::from_parse(s));
    match result {
        Ok(Ok(b)) => {
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, crate::error::OCAS_OK) };
            }
            Box::into_raw(b).cast::<OcasExpr>()
        }
        Ok(Err(msg)) => {
            set(OCAS_ERROR_PARSE, &msg);
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_PARSE) };
            }
            ptr::null_mut()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during parse");
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            ptr::null_mut()
        }
    }
}

/// Free an expression handle. Passing `NULL` is a no-op.
///
/// # Safety
///
/// `handle` must be either null or a pointer returned by a function in this
/// module, and must not have been freed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_free(handle: *mut OcasExpr) {
    if handle.is_null() {
        return;
    }
    crate::error::clear();
    // SAFETY: caller guarantees `handle` is valid and not double-freed.
    unsafe {
        let _ = Box::from_raw(handle.cast::<ExprBox>());
    }
}

/// Clone an expression into a new arena.
///
/// # Safety
///
/// `handle` must be a valid non-null expression handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_clone(
    handle: *const OcasExpr,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let src = expr.atom.to_string();
    finish_op(
        std::panic::catch_unwind(|| ExprBox::from_expr_string(src)),
        err_out,
    )
}

/// Render `handle` to a null-terminated C string.
///
/// The returned string is heap-allocated and owned by the caller. Release it
/// with [`ocas_string_free`]. Returns `NULL` on failure.
///
/// # Safety
///
/// `handle` must be a valid non-null expression handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_to_string(
    handle: *const OcasExpr,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let s = expr.atom.to_string();
    match CString::new(s) {
        Ok(c) => {
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, crate::error::OCAS_OK) };
            }
            c.into_raw()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "expression contained a NUL byte");
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            ptr::null_mut()
        }
    }
}

/// Free a string returned by [`ocas_expr_to_string`]. Passing `NULL` is a
/// no-op.
///
/// # Safety
///
/// `s` must be either null or a pointer returned by [`ocas_expr_to_string`],
/// and must not have been freed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_string_free(s: *mut c_char) {
    if s.is_null() {
        return;
    }
    // SAFETY: caller guarantees `s` came from `CString::into_raw` here.
    unsafe {
        let _ = CString::from_raw(s);
    }
}

/// Re-normalize `handle` in place.
///
/// # Safety
///
/// `handle` must be a valid non-null expression handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_normalize(handle: *mut OcasExpr, err_out: *mut c_int) -> c_int {
    crate::error::clear();
    if handle.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "expression handle is null");
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return OCAS_ERROR_NULL_POINTER;
    }
    // SAFETY: caller guarantees `handle` is valid.
    let expr = unsafe { &mut *handle.cast::<ExprBox>() };
    let ctx = unsafe { static_ctx(expr.ctx_ptr) };
    expr.atom = normalize(ctx, expr.atom);
    if !err_out.is_null() {
        unsafe { ptr::write(err_out, crate::error::OCAS_OK) };
    }
    crate::error::OCAS_OK
}

/// Differentiate `handle` with respect to `var`. Returns a new expression
/// handle (caller owns it) or `NULL` on failure.
///
/// # Safety
///
/// `handle` must be a valid non-null expression handle. `var` must be a
/// valid null-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_diff(
    handle: *const OcasExpr,
    var: *const c_char,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    unary_op(handle, var, err_out, |ctx, a, sym| diff(ctx, a, sym))
}

/// Integrate `handle` with respect to `var`. Returns a new expression handle
/// or `NULL` on failure. If the integral cannot be solved analytically, the
/// result is the unevaluated form `Integral(expr, var)`.
///
/// # Safety
///
/// See [`ocas_expr_diff`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_integrate(
    handle: *const OcasExpr,
    var: *const c_char,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    unary_op(handle, var, err_out, |ctx, a, sym| integrate(ctx, a, sym))
}

/// Compute the Taylor series of `handle` around `point` up to `order`.
///
/// # Safety
///
/// `handle` and `point` must be valid non-null expression handles. `var`
/// must be a valid null-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_taylor(
    handle: *const OcasExpr,
    var: *const c_char,
    point: *const OcasExpr,
    order: u32,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let Some(point_expr) = as_expr(point) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let Some(var_str) = cstr_to_str(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    finish_op(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            expr.apply_taylor(var_str, point_expr, order as usize)
        })),
        err_out,
    )
}

/// Simplify `handle` using the default rule set. Returns a new expression
/// handle or `NULL` on failure.
///
/// # Safety
///
/// `handle` must be a valid non-null expression handle.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_simplify(
    handle: *const OcasExpr,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    finish_op(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| expr.simplify_default())),
        err_out,
    )
}

/// Substitute every occurrence of `var` in `handle` with `replacement`.
///
/// # Safety
///
/// `handle` and `replacement` must be valid non-null expression handles.
/// `var` must be a valid null-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_expr_substitute(
    handle: *const OcasExpr,
    var: *const c_char,
    replacement: *const OcasExpr,
    err_out: *mut c_int,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let Some(repl) = as_expr(replacement) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let Some(var_str) = cstr_to_str(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    finish_op(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            expr.apply_substitute(var_str, repl)
        })),
        err_out,
    )
}

/// Shared body for [`ocas_expr_diff`] and [`ocas_expr_integrate`].
fn unary_op(
    handle: *const OcasExpr,
    var: *const c_char,
    err_out: *mut c_int,
    op: fn(&'static AtomArena<'static>, Atom<'static>, Symbol) -> Atom<'static>,
) -> *mut OcasExpr {
    crate::error::clear();
    let Some(expr) = as_expr(handle) else {
        if !err_out.is_null() {
            unsafe { ptr::write(err_out, OCAS_ERROR_NULL_POINTER) };
        }
        return ptr::null_mut();
    };
    let Some(var_str) = cstr_to_str(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    // SAFETY: we are at an FFI boundary and panics are converted to error
    // codes below; the `&ExprBox` reference is valid for the duration of
    // the call.
    let expr_ref = std::panic::AssertUnwindSafe(expr);
    let var_owned = var_str.to_string();
    finish_op(
        std::panic::catch_unwind(move || {
            let expr: &ExprBox = *expr_ref;
            expr.apply_calculus(&var_owned, op)
        }),
        err_out,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_for_test(s: &str) -> Box<ExprBox> {
        ExprBox::from_parse(s).expect("parse should succeed")
    }

    #[test]
    fn parse_and_roundtrip() {
        let expr = parse_for_test("x^2 + 2*x + 1");
        assert!(expr.atom.to_string().contains("x"));
    }

    #[test]
    fn diff_basic() {
        // d/dx(x^2) = 2*x
        let expr = parse_for_test("x^2");
        let result = expr.apply_calculus("x", diff).unwrap();
        assert_eq!(result.atom.to_string(), "2*x");
    }

    #[test]
    fn integrate_basic() {
        // ∫ 2*x dx. The integrator produces `2*(2^-1)*(x^2)` (coefficient
        // times the power-rule factor `1/(n+1)`); mathematically equal to
        // `x^2`. We assert the structure rather than the unsimplified form.
        let expr = parse_for_test("2*x");
        let result = expr.apply_calculus("x", integrate).unwrap();
        let s = result.atom.to_string();
        assert!(s.contains("(x^2)"), "got: {s}");
        assert!(s.contains("2*(2^-1)"), "got: {s}");
    }

    #[test]
    fn simplify_mul_zero() {
        let expr = parse_for_test("x*0");
        let simplified = expr.simplify_default().unwrap();
        assert_eq!(simplified.atom.to_string(), "0");
    }

    #[test]
    fn clone_produces_equivalent_string() {
        let expr = parse_for_test("sin(x) + cos(x)");
        let src = expr.atom.to_string();
        let cloned = ExprBox::from_expr_string(src).unwrap();
        assert_eq!(expr.atom.to_string(), cloned.atom.to_string());
    }

    #[test]
    fn substitute_replaces_variable() {
        // x^2 + 1 with x -> y. Normalized form puts the constant first:
        // "1 + (y^2)".
        let expr = parse_for_test("x^2 + 1");
        let repl = parse_for_test("y");
        let result = expr.apply_substitute("x", &repl).unwrap();
        assert_eq!(result.atom.to_string(), "1 + (y^2)");
    }
}
