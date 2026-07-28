//! C/C++ bindings for ordinary differential equation solvers.
//!
//! All functions take the ODE as a string expression equal to zero (e.g.
//! `"Derivative(y(x), x) - y(x)"`), the unknown function name (e.g. `"y"`),
//! and the independent variable (e.g. `"x"`). Solutions are returned as
//! heap-allocated strings which must be released with
//! [`ocas_string_free`](crate::expression::ocas_string_free).

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{CString, c_char, c_int};
use std::ptr;

use ocas_atom::Symbol;
use ocas_calc::ode::{
    ODE, ODESolution, ODEType, classify_ode as rs_classify, dsolve as rs_dsolve,
    dsolve_ivp as rs_dsolve_ivp,
};

use crate::error::{OCAS_ERROR_RUNTIME, set};
use crate::expression::{cstr_to_str_pub, expr_ctx_atom, extend_str_lifetime_pub};

/// Classify an ODE and return the applicable method names as a
/// comma-separated string (e.g. `"LinearFirst,Separable,PowerSeries"`).
///
/// Returns `NULL` on failure; the returned string (on success) must be
/// freed with `ocas_string_free`.
///
/// # Safety
///
/// `equation`, `func`, `var` must be valid null-terminated C strings.
/// `err_out` if non-null must point to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ode_classify(
    equation: *const c_char,
    func: *const c_char,
    var: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let Some(eq_s) = cstr_to_str_pub(equation, "equation") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(func_s) = cstr_to_str_pub(func, "func") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(var_s) = cstr_to_str_pub(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let result = std::panic::catch_unwind(|| {
        with_ode_str(eq_s, func_s, var_s, |ctx, ode| {
            rs_classify(ctx, ode)
                .iter()
                .map(ode_type_name)
                .collect::<Vec<_>>()
                .join(",")
        })
    });
    finish_string(result, err_out)
}

/// Solve an ODE symbolically. `hint` may be null (auto classification) or
/// one of the method names returned by [`ocas_ode_classify`].
///
/// Returns a heap-allocated solution string such as `"y = C1*exp(x)"`,
/// or `"unsolved"`. The caller must free it with `ocas_string_free`.
///
/// # Safety
///
/// `equation`, `func`, `var` (and `hint` when non-null) must be valid
/// null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ode_dsolve(
    equation: *const c_char,
    func: *const c_char,
    var: *const c_char,
    hint: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let Some(eq_s) = cstr_to_str_pub(equation, "equation") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(func_s) = cstr_to_str_pub(func, "func") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(var_s) = cstr_to_str_pub(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let hint_type = if hint.is_null() {
        None
    } else {
        match cstr_to_str_pub(hint, "hint").map(parse_ode_type) {
            Some(Some(t)) => Some(t),
            _ => {
                set(OCAS_ERROR_RUNTIME, "unknown ODE type hint");
                crate::error::write_last_code(err_out);
                return ptr::null_mut();
            }
        }
    };
    let result = std::panic::catch_unwind(|| {
        with_ode_str(eq_s, func_s, var_s, |ctx, ode| {
            let sol = rs_dsolve(ctx, ode, hint_type);
            format_solution(&sol)
        })
    });
    finish_string(result, err_out)
}

/// Solve a first- or second-order linear constant-coefficient IVP via the
/// Laplace transform. `y0` is `y(0)` as an expression string (e.g. `"1"`);
/// `y1` is `y'(0)` (may be null for first-order problems).
///
/// Returns a heap-allocated explicit solution string with no free
/// constants; free with `ocas_string_free`.
///
/// # Safety
///
/// All string pointers must be valid null-terminated C strings (`y1` may
/// be null).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ode_dsolve_ivp(
    equation: *const c_char,
    func: *const c_char,
    var: *const c_char,
    y0: *const c_char,
    y1: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let Some(eq_s) = cstr_to_str_pub(equation, "equation") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(func_s) = cstr_to_str_pub(func, "func") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(var_s) = cstr_to_str_pub(var, "var") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let Some(y0_s) = cstr_to_str_pub(y0, "y0") else {
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    let y1_s = if y1.is_null() {
        None
    } else {
        cstr_to_str_pub(y1, "y1")
    };

    let result = std::panic::catch_unwind(|| {
        with_ode_str(eq_s, func_s, var_s, |ctx, ode| {
            let y0_atom = match ocas_parse::parse(ctx, y0_s) {
                Ok(a) => a,
                Err(e) => return format!("error: {e}"),
            };
            let y1_atom = y1_s.and_then(|s| ocas_parse::parse(ctx, s).ok());
            let sol = rs_dsolve_ivp(ctx, ode, y0_atom, y1_atom);
            format_solution(&sol)
        })
    });
    finish_string(result, err_out)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Run `f` with an `ODE` built inside a fresh arena from string inputs.
/// The result must be a plain `String` (no arena borrows).
fn with_ode_str(
    equation: &str,
    func: &str,
    var: &str,
    f: impl for<'a> FnOnce(&'a ocas_atom::AtomArena<'a>, ODE<'a>) -> String,
) -> String {
    let equation = unsafe { extend_str_lifetime_pub(equation) };
    expr_ctx_atom(|ctx| {
        let x = ctx.var(var);
        let func_atom = ctx.fun(func, &[x]);
        let eq_atom = ocas_parse::parse(ctx, equation).unwrap_or_else(|_| ctx.num(0));
        let ode = ODE {
            equation: eq_atom,
            func: func_atom,
            var: Symbol::new(var),
        };
        f(ctx, ode)
    })
}

/// Convert a panic-catching result into a heap-allocated C string.
fn finish_string(result: std::thread::Result<String>, err_out: *mut c_int) -> *mut c_char {
    let s = match result {
        Ok(s) => s,
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during ODE operation");
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            return ptr::null_mut();
        }
    };
    match CString::new(s) {
        Ok(cs) => {
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, crate::error::OCAS_OK) };
            }
            cs.into_raw()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "solution string contains NUL");
            if !err_out.is_null() {
                unsafe { ptr::write(err_out, OCAS_ERROR_RUNTIME) };
            }
            ptr::null_mut()
        }
    }
}

fn ode_type_name(t: &ODEType) -> &'static str {
    match t {
        ODEType::Separable => "Separable",
        ODEType::LinearFirst => "LinearFirst",
        ODEType::Bernoulli => "Bernoulli",
        ODEType::Exact => "Exact",
        ODEType::Homogeneous => "Homogeneous",
        ODEType::LinearConstantCoeff => "LinearConstantCoeff",
        ODEType::CauchyEuler => "CauchyEuler",
        ODEType::ReductionOfOrder => "ReductionOfOrder",
        ODEType::PowerSeries => "PowerSeries",
    }
}

fn parse_ode_type(name: &str) -> Option<ODEType> {
    match name {
        "Separable" => Some(ODEType::Separable),
        "LinearFirst" => Some(ODEType::LinearFirst),
        "Bernoulli" => Some(ODEType::Bernoulli),
        "Exact" => Some(ODEType::Exact),
        "Homogeneous" => Some(ODEType::Homogeneous),
        "LinearConstantCoeff" => Some(ODEType::LinearConstantCoeff),
        "CauchyEuler" => Some(ODEType::CauchyEuler),
        "ReductionOfOrder" => Some(ODEType::ReductionOfOrder),
        "PowerSeries" => Some(ODEType::PowerSeries),
        _ => None,
    }
}

fn format_solution(sol: &ODESolution<'_>) -> String {
    match sol {
        ODESolution::Explicit(e) => format!("y = {e}"),
        ODESolution::Implicit(e) => format!("{e} = C"),
        ODESolution::Parametric(a, b) => format!("x = {a}, y = {b}"),
        ODESolution::Series(e, n) => format!("series({n} terms): y = {e}"),
        ODESolution::System(comps) => {
            let joined: Vec<String> = comps
                .iter()
                .enumerate()
                .map(|(i, c)| format!("y{} = {c}", i + 1))
                .collect();
            joined.join(", ")
        }
        ODESolution::Unsolved(_) => "unsolved".to_string(),
    }
}
