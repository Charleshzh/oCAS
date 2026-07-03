#![allow(clippy::not_unsafe_ptr_arg_deref)]

//! C/C++ bindings for multivariate polynomials.
//!
//! This module exposes sparse multivariate polynomial objects as opaque C
//! handles. Currently supported:
//!
//! - bivariate integer polynomials (`OcasPolyZ`),
//! - bivariate polynomials over a prime finite field (`OcasPolyFp`).
//!
//! Polynomials are created from ASCII strings such as `"x^2 + y + 1"` and
//! can be factored, copied, printed, and freed.

use std::ffi::{CStr, CString, c_char, c_int, c_void};
use std::ptr;

use num_bigint::BigInt;
use num_traits::One;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_core::arena::Arena;
use ocas_domain::{Domain, FiniteField, Integer, IntegerDomain};
use ocas_parse::parse;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::factor::multivariate::{bivariate_factor_fp, bivariate_factor_z};
use ocas_poly::sparse::Lex;

use crate::error::{
    OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER, OCAS_ERROR_PARSE, OCAS_ERROR_RUNTIME,
    OCAS_OK, set,
};

type ZMPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;
type FpMPoly = SparseMultivariatePolynomial<FiniteField, Lex>;

/// Opaque handle for a bivariate integer polynomial.
#[repr(C)]
pub struct OcasPolyZ {
    _private: [u8; 0],
}

/// Opaque handle for a bivariate polynomial over a prime finite field.
#[repr(C)]
pub struct OcasPolyFp {
    _private: [u8; 0],
}

/// A factor returned by [`ocas_poly_z_factor`] or [`ocas_poly_fp_factor`].
///
/// The `poly` pointer is a generic `void*`; the caller must cast it to the
/// appropriate concrete type depending on which factorization function was
/// used (`OcasPolyZ*` for Z, `OcasPolyFp*` for Fp). Each returned polynomial
/// handle must be freed with the corresponding type-specific free function.
#[repr(C)]
pub struct OcasPolyFactor {
    /// The polynomial factor. The actual type is determined by the caller.
    pub poly: *mut c_void,
    /// Multiplicity of the factor.
    pub multiplicity: usize,
}

/// An array of factors returned by [`ocas_poly_z_factor`] or
/// [`ocas_poly_fp_factor`]. Free with [`ocas_poly_factor_array_free`].
#[repr(C)]
pub struct OcasPolyFactorArray {
    /// Pointer to the first factor. May be `NULL` if `len == 0`.
    pub factors: *mut OcasPolyFactor,
    /// Number of factors in the array.
    pub len: usize,
}

// ------------------------------------------------------------------
//  Opaque-handle helpers
// ------------------------------------------------------------------

fn z_ptr(p: Box<ZMPoly>) -> *mut OcasPolyZ {
    Box::into_raw(p) as *mut OcasPolyZ
}

fn z_ref<'a>(p: *const OcasPolyZ) -> Option<&'a ZMPoly> {
    if p.is_null() {
        return None;
    }
    Some(unsafe { &*(p as *const ZMPoly) })
}

fn fp_ptr(p: Box<FpMPoly>) -> *mut OcasPolyFp {
    Box::into_raw(p) as *mut OcasPolyFp
}

fn fp_ref<'a>(p: *const OcasPolyFp) -> Option<&'a FpMPoly> {
    if p.is_null() {
        return None;
    }
    Some(unsafe { &*(p as *const FpMPoly) })
}

// ------------------------------------------------------------------
//  Atom to polynomial conversion
// ------------------------------------------------------------------

fn parse_to_zpoly(input: &str) -> Result<ZMPoly, String> {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).map_err(|e| e.to_string())?;
    atom_to_zpoly(&atom)
}

fn parse_to_fpoly(input: &str, prime: BigInt) -> Result<FpMPoly, String> {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).map_err(|e| e.to_string())?;
    let domain = FiniteField::new(prime);
    atom_to_fpoly(&atom, domain)
}

fn atom_to_zpoly(atom: &Atom<'_>) -> Result<ZMPoly, String> {
    let domain = IntegerDomain;
    let mut result = ZMPoly::new(domain, 2);
    match atom.node() {
        AtomNode::Num(n) => {
            result.set_term_external(vec![0, 0], Integer::from(*n));
        }
        AtomNode::Var(s) => {
            let exp = var_exponent(s)?;
            result.set_term_external(exp, Integer::from(1));
        }
        AtomNode::Add(children) => {
            for child in children.iter() {
                result = result.add(&atom_to_zpoly(child)?);
            }
        }
        AtomNode::Mul(children) => {
            result = ZMPoly::from_terms(domain, 2, vec![(vec![0, 0], Integer::from(1))]);
            for child in children.iter() {
                result = result.mul(&atom_to_zpoly(child)?);
            }
        }
        AtomNode::Pow(base, exp) => {
            let base_poly = atom_to_zpoly(base)?;
            let exp_atom = *exp;
            match exp_atom.node() {
                AtomNode::Num(e) if *e >= 0 => {
                    let power = *e as u64;
                    result = ZMPoly::from_terms(domain, 2, vec![(vec![0, 0], Integer::from(1))]);
                    for _ in 0..power {
                        result = result.mul(&base_poly);
                    }
                }
                _ => return Err("non-negative integer exponent required".to_string()),
            }
        }
        AtomNode::Fun(s, _) => {
            return Err(format!(
                "unsupported function in polynomial: {}",
                s.as_str()
            ));
        }
    }
    Ok(result)
}

fn atom_to_fpoly(atom: &Atom<'_>, domain: FiniteField) -> Result<FpMPoly, String> {
    let mut result = FpMPoly::new(domain.clone(), 2);
    match atom.node() {
        AtomNode::Num(n) => {
            result.set_term_external(vec![0, 0], domain.element(*n));
        }
        AtomNode::Var(s) => {
            let exp = var_exponent(s)?;
            result.set_term_external(exp, domain.one());
        }
        AtomNode::Add(children) => {
            for child in children.iter() {
                result = result.add(&atom_to_fpoly(child, domain.clone())?);
            }
        }
        AtomNode::Mul(children) => {
            result = FpMPoly::from_terms(domain.clone(), 2, vec![(vec![0, 0], domain.one())]);
            for child in children.iter() {
                result = result.mul(&atom_to_fpoly(child, domain.clone())?);
            }
        }
        AtomNode::Pow(base, exp) => {
            let base_poly = atom_to_fpoly(base, domain.clone())?;
            let exp_atom = *exp;
            match exp_atom.node() {
                AtomNode::Num(e) if *e >= 0 => {
                    let power = *e as u64;
                    result =
                        FpMPoly::from_terms(domain.clone(), 2, vec![(vec![0, 0], domain.one())]);
                    for _ in 0..power {
                        result = result.mul(&base_poly);
                    }
                }
                _ => return Err("non-negative integer exponent required".to_string()),
            }
        }
        AtomNode::Fun(s, _) => {
            return Err(format!(
                "unsupported function in polynomial: {}",
                s.as_str()
            ));
        }
    }
    Ok(result)
}

fn var_exponent(s: &Symbol) -> Result<Vec<usize>, String> {
    match s.as_str() {
        "x" => Ok(vec![1, 0]),
        "y" => Ok(vec![0, 1]),
        name => Err(format!("unsupported polynomial variable: {}", name)),
    }
}

fn poly_to_string_z(poly: &ZMPoly) -> String {
    let sorted = poly.sorted_terms();
    if sorted.is_empty() {
        return "0".to_string();
    }
    let mut parts = Vec::new();
    for (exp, coeff) in sorted.iter().rev() {
        let mut s = String::new();
        if coeff.inner().is_one() && exp.iter().any(|&e| e > 0) {
            // omit coefficient 1 for non-constant terms
        } else if coeff.inner() == &BigInt::from(-1) && exp.iter().any(|&e| e > 0) {
            s.push('-');
        } else {
            s.push_str(&coeff.to_string());
        }
        for (var, &e) in exp.iter().enumerate() {
            if e == 0 {
                continue;
            }
            let name = if var == 0 { "x" } else { "y" };
            if e == 1 {
                s.push_str(name);
            } else {
                s.push_str(&format!("{}^{}", name, e));
            }
        }
        if s.is_empty() {
            s.push('1');
        }
        parts.push(s);
    }
    parts.join(" + ").replace(" + -", " - ")
}

fn poly_to_string_fp(poly: &FpMPoly) -> String {
    let sorted = poly.sorted_terms();
    if sorted.is_empty() {
        return "0".to_string();
    }
    let mut parts = Vec::new();
    for (exp, coeff) in sorted.iter().rev() {
        let mut s = String::new();
        let c_val = coeff.value().to_string();
        if c_val == "1" && exp.iter().any(|&e| e > 0) {
            // omit coefficient 1 for non-constant terms
        } else {
            s.push_str(&c_val);
        }
        for (var, &e) in exp.iter().enumerate() {
            if e == 0 {
                continue;
            }
            let name = if var == 0 { "x" } else { "y" };
            if e == 1 {
                s.push_str(name);
            } else {
                s.push_str(&format!("{}^{}", name, e));
            }
        }
        if s.is_empty() {
            s.push('1');
        }
        parts.push(s);
    }
    parts.join(" + ")
}

// ------------------------------------------------------------------
//  Integer polynomial API
// ------------------------------------------------------------------

/// Create a bivariate integer polynomial from a string expression.
///
/// The input may contain the variables `x` and `y`, integer coefficients,
/// addition, multiplication, and non-negative integer powers.
///
/// Returns an opaque handle, or `NULL` on parse failure. On failure the
/// error message can be read with [`ocas_error_last_message`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_create(input: *const c_char, err: *mut c_int) -> *mut OcasPolyZ {
    crate::error::clear();
    if input.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "input string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let input = unsafe { CStr::from_ptr(input) };
    let input = match input.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, "input is not valid UTF-8");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match parse_to_zpoly(input) {
        Ok(poly) => {
            crate::error::write_last_code(err);
            z_ptr(Box::new(poly))
        }
        Err(e) => {
            set(OCAS_ERROR_PARSE, &e);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Free a polynomial handle created with [`ocas_poly_z_create`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_free(poly: *mut OcasPolyZ) {
    if !poly.is_null() {
        unsafe {
            let _ = Box::from_raw(poly as *mut ZMPoly);
        }
    }
}

/// Clone a polynomial handle. Returns a new handle that must be freed.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_clone(poly: *const OcasPolyZ) -> *mut OcasPolyZ {
    match z_ref(poly) {
        Some(p) => z_ptr(Box::new(p.clone())),
        None => ptr::null_mut(),
    }
}

/// Return the total degree of the polynomial, or `0` for the zero polynomial.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_degree(poly: *const OcasPolyZ) -> usize {
    match z_ref(poly) {
        Some(p) => p.total_degree().unwrap_or(0),
        None => 0,
    }
}

/// Return a heap-allocated string representation of the polynomial.
/// The caller must release the returned string with [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_to_string(poly: *const OcasPolyZ, err: *mut c_int) -> *mut c_char {
    crate::error::clear();
    match z_ref(poly) {
        Some(p) => {
            let s = poly_to_string_z(p);
            match CString::new(s) {
                Ok(cs) => {
                    crate::error::write_last_code(err);
                    cs.into_raw()
                }
                Err(_) => {
                    set(OCAS_ERROR_RUNTIME, "failed to create string");
                    crate::error::write_last_code(err);
                    ptr::null_mut()
                }
            }
        }
        None => {
            set(OCAS_ERROR_NULL_POINTER, "polynomial is null");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Factor a bivariate integer polynomial.
///
/// On success, `out` is filled with a heap-allocated array of factors and
/// multiplicities. The caller must release it with
/// [`ocas_poly_factor_array_free`]. On failure `out` is unchanged and the
/// error message can be queried with [`ocas_error_last_message`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_z_factor(
    poly: *const OcasPolyZ,
    out: *mut OcasPolyFactorArray,
    err: *mut c_int,
) -> c_int {
    crate::error::clear();
    if out.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "output array pointer is null");
        crate::error::write_last_code(err);
        return OCAS_ERROR_NULL_POINTER;
    }
    let p = match z_ref(poly) {
        Some(p) => p,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "polynomial is null");
            crate::error::write_last_code(err);
            return OCAS_ERROR_NULL_POINTER;
        }
    };
    let factors = bivariate_factor_z(p, 0, 1);
    let mut array: Vec<OcasPolyFactor> = Vec::with_capacity(factors.len());
    for (poly, multiplicity) in factors {
        array.push(OcasPolyFactor {
            poly: z_ptr(Box::new(poly)) as *mut c_void,
            multiplicity,
        });
    }
    let len = array.len();
    let factors_ptr = array.as_mut_ptr();
    std::mem::forget(array);
    unsafe {
        ptr::write(
            out,
            OcasPolyFactorArray {
                factors: factors_ptr,
                len,
            },
        );
    }
    crate::error::write_last_code(err);
    OCAS_OK
}

// ------------------------------------------------------------------
//  Finite-field polynomial API
// ------------------------------------------------------------------

/// Create a bivariate polynomial over the prime field `F_p` from a string.
///
/// `prime` is the field modulus (must be prime). Coefficients in the string
/// are reduced modulo `p`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_fp_create(
    input: *const c_char,
    prime: *const c_char,
    err: *mut c_int,
) -> *mut OcasPolyFp {
    crate::error::clear();
    if input.is_null() || prime.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "input or prime is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let input = unsafe { CStr::from_ptr(input) };
    let input = match input.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, "input is not valid UTF-8");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let prime = unsafe { CStr::from_ptr(prime) };
    let prime = match prime.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, "prime is not valid UTF-8");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let prime = match prime.parse::<BigInt>() {
        Ok(p) => p,
        Err(_) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, "prime is not a valid integer");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    if prime < BigInt::from(2) {
        set(OCAS_ERROR_INVALID_ARGUMENT, "prime must be at least 2");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    match parse_to_fpoly(input, prime) {
        Ok(poly) => {
            crate::error::write_last_code(err);
            fp_ptr(Box::new(poly))
        }
        Err(e) => {
            set(OCAS_ERROR_PARSE, &e);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Free a finite-field polynomial handle.
#[unsafe(no_mangle)]
#[allow(clippy::not_unsafe_ptr_arg_deref)]
pub extern "C" fn ocas_poly_fp_free(poly: *mut OcasPolyFp) {
    if !poly.is_null() {
        unsafe {
            let _ = Box::from_raw(poly as *mut FpMPoly);
        }
    }
}

/// Clone a finite-field polynomial handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_fp_clone(poly: *const OcasPolyFp) -> *mut OcasPolyFp {
    match fp_ref(poly) {
        Some(p) => fp_ptr(Box::new(p.clone())),
        None => ptr::null_mut(),
    }
}

/// Return the total degree of the finite-field polynomial.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_fp_degree(poly: *const OcasPolyFp) -> usize {
    match fp_ref(poly) {
        Some(p) => p.total_degree().unwrap_or(0),
        None => 0,
    }
}

/// Return a heap-allocated string representation of the finite-field
/// polynomial. The caller must release it with [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_fp_to_string(poly: *const OcasPolyFp, err: *mut c_int) -> *mut c_char {
    crate::error::clear();
    match fp_ref(poly) {
        Some(p) => {
            let s = poly_to_string_fp(p);
            match CString::new(s) {
                Ok(cs) => {
                    crate::error::write_last_code(err);
                    cs.into_raw()
                }
                Err(_) => {
                    set(OCAS_ERROR_RUNTIME, "failed to create string");
                    crate::error::write_last_code(err);
                    ptr::null_mut()
                }
            }
        }
        None => {
            set(OCAS_ERROR_NULL_POINTER, "polynomial is null");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Factor a bivariate polynomial over a prime finite field.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_fp_factor(
    poly: *const OcasPolyFp,
    out: *mut OcasPolyFactorArray,
    err: *mut c_int,
) -> c_int {
    crate::error::clear();
    if out.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "output array pointer is null");
        crate::error::write_last_code(err);
        return OCAS_ERROR_NULL_POINTER;
    }
    let p = match fp_ref(poly) {
        Some(p) => p,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "polynomial is null");
            crate::error::write_last_code(err);
            return OCAS_ERROR_NULL_POINTER;
        }
    };
    let factors = bivariate_factor_fp(p, 0, 1);
    let mut array: Vec<OcasPolyFactor> = Vec::with_capacity(factors.len());
    for (poly, multiplicity) in factors {
        array.push(OcasPolyFactor {
            poly: fp_ptr(Box::new(poly)) as *mut c_void,
            multiplicity,
        });
    }
    let len = array.len();
    let factors_ptr = array.as_mut_ptr();
    std::mem::forget(array);
    unsafe {
        ptr::write(
            out,
            OcasPolyFactorArray {
                factors: factors_ptr,
                len,
            },
        );
    }
    crate::error::write_last_code(err);
    OCAS_OK
}

// ------------------------------------------------------------------
//  Factor array utilities
// ------------------------------------------------------------------

/// Free a factor array returned by [`ocas_poly_z_factor`] or
/// [`ocas_poly_fp_factor`].
///
/// This frees the array structure and the `OcasPolyFactor` objects themselves,
/// but *not* the individual polynomial handles. The caller must free each
/// returned polynomial factor with the appropriate type-specific free function
/// (e.g. [`ocas_poly_z_free`] or [`ocas_poly_fp_free`]).
#[unsafe(no_mangle)]
pub extern "C" fn ocas_poly_factor_array_free(arr: *mut OcasPolyFactorArray) {
    if arr.is_null() {
        return;
    }
    unsafe {
        let len = (*arr).len;
        if !(*arr).factors.is_null() && len > 0 {
            let _ = Vec::from_raw_parts((*arr).factors, len, len);
        }
        (*arr).factors = ptr::null_mut();
        (*arr).len = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn z_create_and_factor() {
        let input = CString::new("x^2 + y + 1").unwrap();
        let mut err: c_int = 0;
        let p = ocas_poly_z_create(input.as_ptr(), &mut err);
        assert!(!p.is_null());
        assert_eq!(err, OCAS_OK);
        let mut factors = OcasPolyFactorArray {
            factors: ptr::null_mut(),
            len: 0,
        };
        let rc = ocas_poly_z_factor(p, &mut factors, &mut err);
        assert_eq!(rc, OCAS_OK);
        ocas_poly_factor_array_free(&mut factors);
        ocas_poly_z_free(p);
    }

    #[test]
    fn fp_create_and_factor() {
        let input = CString::new("x^2 + y + 1").unwrap();
        let prime = CString::new("5").unwrap();
        let mut err: c_int = 0;
        let p = ocas_poly_fp_create(input.as_ptr(), prime.as_ptr(), &mut err);
        assert!(!p.is_null());
        assert_eq!(err, OCAS_OK);
        let mut factors = OcasPolyFactorArray {
            factors: ptr::null_mut(),
            len: 0,
        };
        let rc = ocas_poly_fp_factor(p, &mut factors, &mut err);
        assert_eq!(rc, OCAS_OK);
        ocas_poly_factor_array_free(&mut factors);
        ocas_poly_fp_free(p);
    }
}
