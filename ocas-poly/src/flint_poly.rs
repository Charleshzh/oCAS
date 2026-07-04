//! Experimental FLINT 3 backend for dense univariate polynomials.
//!
//! This module is only compiled when the `flint` feature is enabled. It
//! provides a thin wrapper around FLINT's `fmpz_poly_t` (integer polynomial)
//! type and delegates basic arithmetic to FLINT.
//!
//! Build requirements:
//! - Linux: `libflint-dev` (Debian/Ubuntu) or equivalent.
//! - Windows: MSYS2 MINGW64 with `mingw-w64-x86_64-flint`, or vcpkg.
//!
//! On platforms where FLINT is unavailable, omit the `flint` feature to use
//! the pure-Rust fallback in [`DenseUnivariatePolynomial`].

#![cfg(feature = "flint")]

use std::ffi::{CStr, CString};
use std::mem::MaybeUninit;

use flint3_sys::{
    c_char, c_void, fmpz, fmpz_poly_add, fmpz_poly_clear, fmpz_poly_init, fmpz_poly_mul,
    fmpz_poly_scalar_mul_fmpz, fmpz_poly_set_coeff_fmpz, fmpz_poly_sub, fmpz_poly_t,
};

use crate::dense::DenseUnivariatePolynomial;
#[cfg(not(feature = "gmp"))]
use num_bigint::BigInt;
use ocas_domain::{Integer, IntegerDomain};

/// A dense univariate polynomial backed by FLINT's `fmpz_poly_t`.
///
/// This type is intentionally separate from [`DenseUnivariatePolynomial`] so
/// that the pure-Rust API remains available when the `flint` feature is off.
pub struct FlintUnivariatePolynomial {
    raw: fmpz_poly_t,
    domain: IntegerDomain,
}

impl std::fmt::Debug for FlintUnivariatePolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FlintUnivariatePolynomial")
            .field("domain", &self.domain)
            .finish_non_exhaustive()
    }
}

impl FlintUnivariatePolynomial {
    /// Create the zero polynomial.
    pub fn new() -> Self {
        // SAFETY: fmpz_poly_t is a C array type; zeroed memory is a valid
        // initial state for fmpz_poly_init.
        let mut raw: fmpz_poly_t = unsafe { std::mem::zeroed() };
        unsafe {
            fmpz_poly_init(raw.as_mut_ptr());
            Self {
                raw,
                domain: IntegerDomain,
            }
        }
    }

    /// Create a FLINT polynomial from the pure-Rust dense representation.
    pub fn from_dense(poly: &DenseUnivariatePolynomial<IntegerDomain>) -> Self {
        let mut result = Self::new();
        for (i, coeff) in poly.coeffs().iter().enumerate() {
            let mut c = fmpz_from_integer(coeff);
            unsafe {
                fmpz_poly_set_coeff_fmpz(result.raw.as_mut_ptr(), i as i64, &c);
                flint3_sys::fmpz_clear(&mut c);
            }
        }
        result
    }

    /// Convert back to the pure-Rust dense representation.
    pub fn to_dense(&self) -> DenseUnivariatePolynomial<IntegerDomain> {
        let degree = unsafe { flint3_sys::fmpz_poly_degree(self.raw.as_ptr()) };
        let mut coeffs = Vec::with_capacity(if degree < 0 { 0 } else { degree as usize + 1 });
        for i in 0..=degree {
            let mut c = MaybeUninit::<fmpz>::uninit();
            unsafe {
                flint3_sys::fmpz_poly_get_coeff_fmpz(c.as_mut_ptr(), self.raw.as_ptr(), i);
                let mut c = c.assume_init();
                coeffs.push(integer_from_fmpz(&c));
                flint3_sys::fmpz_clear(&mut c);
            }
        }
        DenseUnivariatePolynomial::from_coeffs(IntegerDomain, coeffs)
    }

    /// Add two FLINT polynomials, returning a new polynomial.
    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::new();
        unsafe {
            fmpz_poly_add(
                result.raw.as_mut_ptr(),
                self.raw.as_ptr(),
                other.raw.as_ptr(),
            );
        }
        result
    }

    /// Subtract two FLINT polynomials, returning a new polynomial.
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = Self::new();
        unsafe {
            fmpz_poly_sub(
                result.raw.as_mut_ptr(),
                self.raw.as_ptr(),
                other.raw.as_ptr(),
            );
        }
        result
    }

    /// Multiply two FLINT polynomials, returning a new polynomial.
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = Self::new();
        unsafe {
            fmpz_poly_mul(
                result.raw.as_mut_ptr(),
                self.raw.as_ptr(),
                other.raw.as_ptr(),
            );
        }
        result
    }

    /// Multiply by a scalar integer.
    pub fn mul_scalar(&self, scalar: &Integer) -> Self {
        let mut result = Self::new();
        let mut c = fmpz_from_integer(scalar);
        unsafe {
            fmpz_poly_scalar_mul_fmpz(result.raw.as_mut_ptr(), self.raw.as_ptr(), &c);
            flint3_sys::fmpz_clear(&mut c);
        }
        result
    }
}

impl Drop for FlintUnivariatePolynomial {
    fn drop(&mut self) {
        unsafe {
            fmpz_poly_clear(self.raw.as_mut_ptr());
        }
    }
}

fn fmpz_from_integer(i: &Integer) -> fmpz {
    let s = i.to_string();
    let c_str = CString::new(s).expect("BigInt decimal string contains null byte");
    unsafe {
        let mut z = MaybeUninit::<fmpz>::uninit();
        flint3_sys::fmpz_init(z.as_mut_ptr());
        // Base 10. fmpz_set_str returns 0 on success; the string is always valid decimal.
        let _ = flint3_sys::fmpz_set_str(z.as_mut_ptr(), c_str.as_ptr() as *const c_char, 10);
        z.assume_init()
    }
}

fn integer_from_fmpz(z: &fmpz) -> Integer {
    unsafe {
        // fmpz_get_str with a null buffer allocates a string using flint_malloc.
        let c_str = flint3_sys::fmpz_get_str(std::ptr::null_mut(), 10, z);
        let bytes = CStr::from_ptr(c_str).to_bytes();
        let s = std::str::from_utf8(bytes).expect("FLINT produced a valid decimal integer string");
        #[cfg(not(feature = "gmp"))]
        let value =
            BigInt::parse_bytes(bytes, 10).expect("FLINT produced a valid decimal integer string");
        #[cfg(feature = "gmp")]
        let value = rug::Integer::from_str_radix(s, 10)
            .expect("FLINT produced a valid decimal integer string");
        flint3_sys::flint_free(c_str as *mut c_void);
        Integer::new(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use ocas_domain::Integer;

    #[test]
    fn flint_integer_round_trip() {
        fn check(value: i64) {
            let i = Integer::from(value);
            let mut z = fmpz_from_integer(&i);
            let round = integer_from_fmpz(&z);
            unsafe {
                flint3_sys::fmpz_clear(&mut z);
            }
            assert_eq!(round, i);
        }
        check(0);
        check(1);
        check(-1);
        check(42);
        check(-42);
        check(i64::MAX);
        check(i64::MIN);
    }

    #[test]
    fn flint_large_coefficient_round_trip() {
        let big = Integer::from(BigInt::from(10).pow(50));
        let mut z = fmpz_from_integer(&big);
        let round = integer_from_fmpz(&z);
        unsafe {
            flint3_sys::fmpz_clear(&mut z);
        }
        assert_eq!(round, big);
    }

    #[test]
    fn flint_add_matches_dense() {
        let a = int_poly(&[1, 2]);
        let b = int_poly(&[3, 0, 4]);
        let fa = FlintUnivariatePolynomial::from_dense(&a);
        let fb = FlintUnivariatePolynomial::from_dense(&b);
        assert_eq!(fa.add(&fb).to_dense(), a.add(&b));
    }

    #[test]
    fn flint_sub_matches_dense() {
        let a = int_poly(&[1, 2]);
        let b = int_poly(&[3, 0, 4]);
        let fa = FlintUnivariatePolynomial::from_dense(&a);
        let fb = FlintUnivariatePolynomial::from_dense(&b);
        assert_eq!(fa.sub(&fb).to_dense(), a.sub(&b));
    }

    #[test]
    fn flint_mul_matches_dense() {
        let a = int_poly(&[1, 2]);
        let b = int_poly(&[3, 0, 4]);
        let fa = FlintUnivariatePolynomial::from_dense(&a);
        let fb = FlintUnivariatePolynomial::from_dense(&b);
        assert_eq!(fa.mul(&fb).to_dense(), a.mul(&b));
    }

    #[test]
    fn flint_scalar_mul_matches_dense() {
        let a = int_poly(&[1, 2, 3]);
        let s = Integer::from(7);
        let fa = FlintUnivariatePolynomial::from_dense(&a);
        assert_eq!(fa.mul_scalar(&s).to_dense(), a.mul_scalar(&s));
    }

    fn int_poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<IntegerDomain> {
        DenseUnivariatePolynomial::from_coeffs(
            IntegerDomain,
            coeffs.iter().copied().map(Integer::from).collect(),
        )
    }
}
