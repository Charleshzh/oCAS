//! C/C++ bindings for numerical integration (Vegas adaptive Monte Carlo).
//!
//! This module exposes the Vegas adaptive Monte Carlo integrator as an
//! opaque [`OcasVegas`] handle plus a one-shot [`ocas_integrate_1d`]
//! convenience function.
//!
//! # Integrand convention
//!
//! C callers pass an integrand as a function pointer plus a `user_data`
//! pointer which is forwarded untouched to the integrand:
//!
//! ```c
//! typedef double (*ocas_integrand_t)(double x, void *user_data);
//! ```
//!
//! For the multi-dimensional [`ocas_vegas_integrate`], `x` is the first
//! coordinate of the sampled point only; this matches the 1-D focus of the
//! current API. (A future revision may expose a `const double*`/`size_t`
//! signature for full n-D Vegas.)
//!
//! # Example
//!
//! ```c
//! static double f(double x, void *ud) { return x; }
//!
//! int err = 0;
//! OcasVegas *v = ocas_vegas_create(1, &(OcasVegasOptions){.n_samples=20000,
//!                                                          .iterations=8},
//!                                   &err);
//! OcasIntegrateResult r = ocas_vegas_integrate(v, f, NULL, &err);
//! printf("integral = %g, error = %g\n", r.integral, r.error);
//! ocas_vegas_free(v);
//! ```

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::c_int;
use std::ptr;

use ocas_eval::numeric::{IntegrateResult, Integrator, Vegas, VegasOptions};

use crate::error::{OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER, set};

/// Opaque handle for a Vegas integrator.
#[repr(C)]
pub struct OcasVegas {
    _private: [u8; 0],
}

/// Result of a numerical integration: the estimate and its standard error.
#[repr(C)]
pub struct OcasIntegrateResult {
    /// Best estimate of the integral.
    pub integral: f64,
    /// Estimated standard error on `integral`.
    pub error: f64,
}

/// Tuning knobs for [`ocas_vegas_create`]. Pass zeros via the pointer to use
/// the library defaults (same as [`ocas_vegas_create_default`]).
#[repr(C)]
pub struct OcasVegasOptions {
    /// Number of bins per dimension (default 64).
    pub n_bins: usize,
    /// Number of samples per iteration (default 10000).
    pub n_samples: usize,
    /// Number of adaptive iterations (default 10).
    pub iterations: usize,
    /// Grid smoothing / learning rate (default 1.5).
    pub learning_rate: f64,
    /// RNG seed (default 0x0C45).
    pub seed: u64,
}

/// Integrand function pointer for one-dimensional integration.
///
/// `user_data` is passed through untouched from the caller of
/// [`ocas_vegas_integrate`] or [`ocas_integrate_1d`].
#[allow(non_camel_case_types)]
pub type ocas_integrand_t =
    Option<unsafe extern "C" fn(x: f64, user_data: *mut std::ffi::c_void) -> f64>;

// ------------------------------------------------------------------
//  Opaque-handle helpers
// ------------------------------------------------------------------

fn vegas_ptr(v: Box<Vegas>) -> *mut OcasVegas {
    Box::into_raw(v) as *mut OcasVegas
}

fn vegas_ref<'a>(v: *const OcasVegas) -> Option<&'a Vegas> {
    if v.is_null() {
        return None;
    }
    Some(unsafe { &*(v as *const Vegas) })
}

fn vegas_mut_ref<'a>(v: *mut OcasVegas) -> Option<&'a mut Vegas> {
    if v.is_null() {
        return None;
    }
    Some(unsafe { &mut *(v as *mut Vegas) })
}

fn opts_from_ptr(opts: *const OcasVegasOptions) -> VegasOptions {
    if opts.is_null() {
        return VegasOptions::default();
    }
    let raw = unsafe { &*opts };
    let mut o = VegasOptions::default();
    if raw.n_bins != 0 {
        o.n_bins = raw.n_bins;
    }
    if raw.n_samples != 0 {
        o.n_samples = raw.n_samples;
    }
    if raw.iterations != 0 {
        o.iterations = raw.iterations;
    }
    if raw.learning_rate > 0.0 {
        o.learning_rate = raw.learning_rate;
    }
    if raw.seed != 0 {
        o.seed = raw.seed;
    }
    o
}

impl From<IntegrateResult> for OcasIntegrateResult {
    fn from(r: IntegrateResult) -> Self {
        Self {
            integral: r.integral,
            error: r.error,
        }
    }
}

// ------------------------------------------------------------------
//  Vegas C API
// ------------------------------------------------------------------

/// Create a Vegas integrator for `n_dims` dimensions.
///
/// `opts` may be `NULL` to use library defaults. Returns an opaque handle,
/// or `NULL` on failure (with the error code written to `*err`).
#[unsafe(no_mangle)]
pub extern "C" fn ocas_vegas_create(
    n_dims: usize,
    opts: *const OcasVegasOptions,
    err: *mut c_int,
) -> *mut OcasVegas {
    crate::error::clear();
    if n_dims == 0 {
        set(OCAS_ERROR_INVALID_ARGUMENT, "n_dims must be >= 1");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let options = opts_from_ptr(opts);
    if options.n_bins == 0 {
        set(OCAS_ERROR_INVALID_ARGUMENT, "n_bins must be >= 1");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let vegas = Vegas::new(n_dims, options);
    crate::error::write_last_code(err);
    vegas_ptr(Box::new(vegas))
}

/// Free a Vegas integrator handle. Safe to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_vegas_free(v: *mut OcasVegas) {
    if v.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(v as *mut Vegas));
    }
}

/// Integrate `f` over the unit hypercube, invoking `f(x, user_data)` for each
/// sample. Returns the result; writes the error code to `*err` on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_vegas_integrate(
    v: *mut OcasVegas,
    f: ocas_integrand_t,
    user_data: *mut std::ffi::c_void,
    err: *mut c_int,
) -> OcasIntegrateResult {
    crate::error::clear();
    let integrand = match f {
        Some(g) => g,
        None => {
            set(
                OCAS_ERROR_NULL_POINTER,
                "integrand function pointer is NULL",
            );
            crate::error::write_last_code(err);
            return OcasIntegrateResult {
                integral: f64::NAN,
                error: f64::NAN,
            };
        }
    };
    let vegas = match vegas_mut_ref(v) {
        Some(v) => v,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "Vegas handle is NULL");
            crate::error::write_last_code(err);
            return OcasIntegrateResult {
                integral: f64::NAN,
                error: f64::NAN,
            };
        }
    };
    let wrapped = |x: &[f64]| -> f64 {
        // SAFETY: caller asserts `integrand` is a valid function and is
        // responsible for `user_data`'s validity for the duration of the call.
        unsafe { integrand(x[0], user_data) }
    };
    let r = vegas.integrate(&wrapped);
    crate::error::write_last_code(err);
    r.into()
}

/// Latest accumulated estimate and error after [`ocas_vegas_integrate`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_vegas_result(v: *const OcasVegas) -> OcasIntegrateResult {
    match vegas_ref(v) {
        Some(vegas) => vegas.result().into(),
        None => OcasIntegrateResult {
            integral: f64::NAN,
            error: f64::NAN,
        },
    }
}

/// Number of completed iterations, or `0` on a null handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_vegas_iterations(v: *const OcasVegas) -> usize {
    match vegas_ref(v) {
        Some(vegas) => vegas.iterations(),
        None => 0,
    }
}

// ------------------------------------------------------------------
//  One-shot integrate_1d
// ------------------------------------------------------------------

/// Numerically integrate `f` over `[a, b]` using Vegas in one shot.
///
/// `opts` may be `NULL` to use library defaults. Returns the result; writes
/// the error code to `*err` on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_integrate_1d(
    f: ocas_integrand_t,
    user_data: *mut std::ffi::c_void,
    a: f64,
    b: f64,
    opts: *const OcasVegasOptions,
    err: *mut c_int,
) -> OcasIntegrateResult {
    crate::error::clear();
    let integrand = match f {
        Some(g) => g,
        None => {
            set(
                OCAS_ERROR_NULL_POINTER,
                "integrand function pointer is NULL",
            );
            crate::error::write_last_code(err);
            return OcasIntegrateResult {
                integral: f64::NAN,
                error: f64::NAN,
            };
        }
    };
    if a.partial_cmp(&b) != Some(std::cmp::Ordering::Less) {
        set(
            OCAS_ERROR_INVALID_ARGUMENT,
            "integration upper bound b must be > a",
        );
        crate::error::write_last_code(err);
        return OcasIntegrateResult {
            integral: f64::NAN,
            error: f64::NAN,
        };
    }
    let options = opts_from_ptr(opts);
    let wrapped = |x: f64| -> f64 {
        // SAFETY: see [`ocas_vegas_integrate`].
        unsafe { integrand(x, user_data) }
    };
    let r = ocas_eval::numeric::integrate_1d(wrapped, a, b, options);
    crate::error::write_last_code(err);
    r.into()
}
