//! C bindings for Gröbner basis computation and ideal operations.
//!
//! Provides FFI functions for computing Gröbner bases and performing
//! ideal arithmetic over the rationals.
//!
//! ## Design
//!
//! Polynomials are passed as coefficient arrays: each polynomial is specified
//! by its number of variables, number of terms, a flat exponent matrix, and
//! rational coefficient arrays (numerator/denominator pairs).

use ocas_domain::{Rational, RationalDomain};
use ocas_poly::ideal;
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, GroebnerBasis, SparseMultivariatePolynomial, groebner_basis};

use crate::error;

/// Opaque handle to a Gröbner basis.
pub struct OcasGroebnerBasis {
    pub(crate) inner: GroebnerBasis<RationalDomain, Lex>,
}

/// Opaque handle to a polynomial system solution.
pub struct OcasSystemSolution {
    pub(crate) inner: ideal::PolynomialSystemSolution,
}

/// Parse polynomial data from arrays into Rust types.
#[allow(clippy::needless_range_loop)]
fn build_polys(
    n_polys: usize,
    n_vars_array: *const usize,
    n_terms_array: *const usize,
    exponents: *const usize,
    coeff_nums: *const i64,
    coeff_dens: *const i64,
) -> Result<Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>, &'static str> {
    let mut result = Vec::with_capacity(n_polys);
    let mut exp_offset = 0usize;
    let mut coeff_offset = 0usize;

    for p in 0..n_polys {
        let n_vars = unsafe { *n_vars_array.add(p) };
        let n_terms = unsafe { *n_terms_array.add(p) };
        let mut terms = Vec::with_capacity(n_terms);

        for _ in 0..n_terms {
            let mut exp = vec![0usize; n_vars];
            for v in 0..n_vars {
                exp[v] = unsafe { *exponents.add(exp_offset + v) };
            }
            exp_offset += n_vars;

            let num = unsafe { *coeff_nums.add(coeff_offset) };
            let den = unsafe { *coeff_dens.add(coeff_offset) };
            coeff_offset += 1;

            if den == 0 {
                return Err("zero denominator");
            }
            terms.push((exp, Rational::new(num, den)));
        }

        result.push(SparseMultivariatePolynomial::from_terms(
            RationalDomain,
            n_vars,
            terms,
        ));
    }
    Ok(result)
}

/// Compute a Gröbner basis from polynomial data arrays.
///
/// # Safety
/// All array pointers must be valid for the specified lengths.
/// The returned handle must be freed with `ocas_groebner_basis_free`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_groebner_basis(
    n_polys: usize,
    n_vars_array: *const usize,
    n_terms_array: *const usize,
    exponents: *const usize,
    coeff_nums: *const i64,
    coeff_dens: *const i64,
    algorithm: i32,
    err: *mut i32,
) -> *mut OcasGroebnerBasis {
    error::clear();
    if err.is_null() {
        return std::ptr::null_mut();
    }

    let algo = match algorithm {
        0 => Algorithm::Auto,
        1 => Algorithm::F4,
        2 => Algorithm::F5,
        3 => Algorithm::Buchberger,
        4 => Algorithm::MultiModular,
        _ => {
            error::set(-1, "invalid algorithm: expected 0-4");
            unsafe {
                *err = -1;
            }
            return std::ptr::null_mut();
        }
    };

    let gens = match build_polys(
        n_polys,
        n_vars_array,
        n_terms_array,
        exponents,
        coeff_nums,
        coeff_dens,
    ) {
        Ok(g) => g,
        Err(e) => {
            error::set(-1, e);
            unsafe {
                *err = -1;
            }
            return std::ptr::null_mut();
        }
    };

    let gb = groebner_basis(&gens, algo);
    unsafe {
        *err = 0;
    }
    Box::into_raw(Box::new(OcasGroebnerBasis { inner: gb }))
}

/// Free a Gröbner basis handle.
///
/// # Safety
/// `gb` must have been returned by `ocas_groebner_basis` and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_groebner_basis_free(gb: *mut OcasGroebnerBasis) {
    if !gb.is_null() {
        drop(unsafe { Box::from_raw(gb) });
    }
}

/// Get the number of elements in a Gröbner basis.
///
/// # Safety
/// `gb` must be a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_groebner_basis_len(gb: *const OcasGroebnerBasis) -> usize {
    if gb.is_null() {
        0
    } else {
        unsafe { &*gb }.inner.basis.len()
    }
}

/// Check if an ideal is zero-dimensional.
///
/// # Safety
/// `gb` must be a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_is_zero_dimensional(gb: *const OcasGroebnerBasis) -> bool {
    if gb.is_null() {
        false
    } else {
        ideal::is_zero_dimensional(unsafe { &(*gb).inner })
    }
}

/// Solve a polynomial system from data arrays.
///
/// Returns a handle to the solution, or NULL on error.
/// Free with `ocas_system_solution_free`.
///
/// # Safety
/// All array pointers must be valid.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_solve_polynomial_system(
    n_polys: usize,
    n_vars_array: *const usize,
    n_terms_array: *const usize,
    exponents: *const usize,
    coeff_nums: *const i64,
    coeff_dens: *const i64,
    algorithm: i32,
    err: *mut i32,
) -> *mut OcasSystemSolution {
    error::clear();
    if err.is_null() {
        return std::ptr::null_mut();
    }

    let algo = match algorithm {
        0 => Algorithm::Auto,
        1 => Algorithm::F4,
        2 => Algorithm::F5,
        3 => Algorithm::Buchberger,
        4 => Algorithm::MultiModular,
        _ => {
            error::set(-1, "invalid algorithm: expected 0-4");
            unsafe {
                *err = -1;
            }
            return std::ptr::null_mut();
        }
    };

    let gens = match build_polys(
        n_polys,
        n_vars_array,
        n_terms_array,
        exponents,
        coeff_nums,
        coeff_dens,
    ) {
        Ok(g) => g,
        Err(e) => {
            error::set(-1, e);
            unsafe {
                *err = -1;
            }
            return std::ptr::null_mut();
        }
    };

    let sol = ideal::solve_polynomial_system(&gens, algo);
    unsafe {
        *err = 0;
    }
    Box::into_raw(Box::new(OcasSystemSolution { inner: sol }))
}

/// Get the number of solutions.
///
/// Returns 0 for non-zero-dimensional or empty solutions.
///
/// # Safety
/// `sol` must be a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_system_solution_count(sol: *const OcasSystemSolution) -> usize {
    if sol.is_null() {
        return 0;
    }
    match &unsafe { &*sol }.inner {
        ideal::PolynomialSystemSolution::ZeroDimensional(z) => z.solutions.len(),
        _ => 0,
    }
}

/// Get a specific solution value.
///
/// Returns the value as f64, or 0.0 on error.
///
/// # Safety
/// `sol` must be valid; indices must be in bounds.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_system_solution_value(
    sol: *const OcasSystemSolution,
    sol_idx: usize,
    var_idx: usize,
) -> f64 {
    if sol.is_null() {
        return 0.0;
    }
    match &unsafe { &*sol }.inner {
        ideal::PolynomialSystemSolution::ZeroDimensional(z) => z
            .solutions
            .get(sol_idx)
            .and_then(|s| s.values.get(var_idx))
            .copied()
            .unwrap_or(0.0),
        _ => 0.0,
    }
}

/// Free a system solution handle.
///
/// # Safety
/// `sol` must have been returned by `ocas_solve_polynomial_system`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_system_solution_free(sol: *mut OcasSystemSolution) {
    if !sol.is_null() {
        drop(unsafe { Box::from_raw(sol) });
    }
}
