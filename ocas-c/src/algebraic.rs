//! C/C++ bindings for univariate polynomials over algebraic number fields.
//!
//! This module exposes Trager factorization over $\mathbb{Q}(\alpha)$ as
//! opaque C handles. An [`OcasAlgebraicField`] is created from a minimal
//! polynomial string (e.g. `"x^2 - 2"` for $\mathbb{Q}(\sqrt{2})$); an
//! [`OcasAlgebraicPoly`] is created from a coefficient-list string.
//!
//! # Coefficient string format
//!
//! Polynomial coefficients are separated by `;` (constant term first). Each
//! coefficient is a comma-separated list of rationals `n/d` giving the
//! $\alpha$-polynomial in ascending degree order. For example, over
//! $\mathbb{Q}(\sqrt{2})$:
//!
//! - `"-2;0;1"` — $x^2 - 2$ (all base-domain constants).
//! - `"0;0;0,1"` — $x^2 - \alpha$ (the $x^2$ coefficient is $0 + 1\cdot\alpha$).
//!
//! A singleton rational may omit the comma.

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{CStr, CString, c_char, c_int, c_void};
use std::ptr;

use ocas_atom::{Atom, AtomArena, AtomNode};
use ocas_core::arena::Arena;
use ocas_domain::{AlgebraicElement, AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_parse::parse;
use ocas_poly::DenseUnivariatePolynomial;

use crate::error::{
    OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER, OCAS_ERROR_PARSE, OCAS_ERROR_RUNTIME,
    OCAS_OK, set,
};

type AnfPoly = DenseUnivariatePolynomial<AlgebraicNumberField>;

// ------------------------------------------------------------------
//  Opaque handles
// ------------------------------------------------------------------

/// Opaque handle for an algebraic number field $\mathbb{Q}(\alpha)$.
#[repr(C)]
pub struct OcasAlgebraicField {
    _private: [u8; 0],
}

/// Opaque handle for a univariate polynomial over an algebraic number field.
#[repr(C)]
pub struct OcasAlgebraicPoly {
    _private: [u8; 0],
}

/// A factor returned by [`ocas_algebraic_poly_factor`]. The `poly` pointer
/// must be cast to `OcasAlgebraicPoly*` and freed with
/// [`ocas_algebraic_poly_free`].
#[repr(C)]
pub struct OcasAlgebraicFactor {
    /// The polynomial factor.
    pub poly: *mut c_void,
    /// Multiplicity of the factor.
    pub multiplicity: usize,
}

/// An array of factors returned by [`ocas_algebraic_poly_factor`]. Free with
/// [`ocas_algebraic_factor_array_free`].
#[repr(C)]
pub struct OcasAlgebraicFactorArray {
    /// Pointer to the first factor. May be `NULL` if `len == 0`.
    pub factors: *mut OcasAlgebraicFactor,
    /// Number of factors in the array.
    pub len: usize,
}

// ------------------------------------------------------------------
//  Opaque-handle helpers
// ------------------------------------------------------------------

fn field_ptr(f: Box<AlgebraicNumberField>) -> *mut OcasAlgebraicField {
    Box::into_raw(f) as *mut OcasAlgebraicField
}

fn field_ref<'a>(f: *const OcasAlgebraicField) -> Option<&'a AlgebraicNumberField> {
    if f.is_null() {
        return None;
    }
    Some(unsafe { &*(f as *const AlgebraicNumberField) })
}

fn poly_ptr(p: Box<AnfPoly>) -> *mut OcasAlgebraicPoly {
    Box::into_raw(p) as *mut OcasAlgebraicPoly
}

fn poly_ref<'a>(p: *const OcasAlgebraicPoly) -> Option<&'a AnfPoly> {
    if p.is_null() {
        return None;
    }
    Some(unsafe { &*(p as *const AnfPoly) })
}

// ------------------------------------------------------------------
//  Minimal-polynomial parsing (Atom -> ascending rational coefficients)
// ------------------------------------------------------------------

/// Add the contribution of `atom` (scaled by `mult`) into `coeffs`.
fn accumulate(
    dom: &RationalDomain,
    coeffs: &mut Vec<Rational>,
    atom: &Atom<'_>,
    mult: Rational,
) -> Result<(), String> {
    match atom.node() {
        AtomNode::Num(n) => {
            let c = dom.mul(&mult, &Rational::new(*n, 1));
            set_coeff(dom, coeffs, 0, &c);
        }
        AtomNode::Var(s) if s.as_str() == "x" => {
            let c = dom.mul(&mult, &dom.one());
            set_coeff(dom, coeffs, 1, &c);
        }
        AtomNode::Pow(base, exp) => {
            if let (AtomNode::Var(s), AtomNode::Num(e)) = (base.node(), exp.node())
                && s.as_str() == "x"
                && *e >= 0
            {
                let c = dom.mul(&mult, &dom.one());
                set_coeff(dom, coeffs, *e as usize, &c);
                return Ok(());
            }
            return Err("unsupported power term in minimal polynomial".to_string());
        }
        AtomNode::Add(children) => {
            for child in children.iter() {
                accumulate(dom, coeffs, child, mult.clone())?;
            }
        }
        AtomNode::Mul(children) => {
            let mut coef = mult;
            let mut xpow: usize = 0;
            for child in children.iter() {
                match child.node() {
                    AtomNode::Num(n) => {
                        coef = dom.mul(&coef, &Rational::new(*n, 1));
                    }
                    AtomNode::Var(s) if s.as_str() == "x" => {
                        xpow += 1;
                    }
                    AtomNode::Pow(base, exp)
                        if matches!(
                            (base.node(), exp.node()),
                            (AtomNode::Var(s), AtomNode::Num(e)) if s.as_str() == "x" && *e >= 0
                        ) =>
                    {
                        if let AtomNode::Num(e) = exp.node() {
                            xpow += *e as usize;
                        }
                    }
                    _ => return Err("unsupported factor in minimal polynomial".to_string()),
                }
            }
            set_coeff(dom, coeffs, xpow, &coef);
        }
        _ => return Err("unsupported term in minimal polynomial".to_string()),
    }
    Ok(())
}

fn set_coeff(dom: &RationalDomain, coeffs: &mut Vec<Rational>, i: usize, c: &Rational) {
    while coeffs.len() <= i {
        coeffs.push(dom.zero());
    }
    coeffs[i] = dom.add(&coeffs[i], c);
}

/// Parse a univariate (in `x`) polynomial string into ascending rational
/// coefficients, trimming trailing zeros.
fn parse_min_poly_str(s: &str) -> Result<Vec<Rational>, String> {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, s).map_err(|e| e.to_string())?;
    let dom = RationalDomain;
    let mut coeffs: Vec<Rational> = Vec::new();
    accumulate(&dom, &mut coeffs, &atom, dom.one())?;
    while let Some(last) = coeffs.last() {
        if dom.is_zero(last) {
            coeffs.pop();
        } else {
            break;
        }
    }
    Ok(coeffs)
}

/// Parse a single rational of the form `n` or `n/d`.
fn parse_single_rational(s: &str) -> Result<Rational, String> {
    let s = s.trim();
    if let Some((n, d)) = s.split_once('/') {
        let num: i64 = n
            .trim()
            .parse()
            .map_err(|_| format!("invalid rational numerator: {n}"))?;
        let den: i64 = d
            .trim()
            .parse()
            .map_err(|_| format!("invalid rational denominator: {d}"))?;
        if den == 0 {
            return Err("rational denominator cannot be zero".to_string());
        }
        Ok(Rational::new(num, den))
    } else {
        let n: i64 = s.parse().map_err(|_| format!("invalid rational: {s}"))?;
        Ok(Rational::new(n, 1))
    }
}

/// Parse a coefficient-list string into algebraic-field elements.
///
/// Format: polynomial coefficients separated by `;` (constant first); each
/// coefficient is a comma-separated list of `n/d` rationals (ascending
/// $\alpha$-polynomial).
fn parse_anf_coeffs(
    field: &AlgebraicNumberField,
    s: &str,
) -> Result<Vec<AlgebraicElement<Rational>>, String> {
    let mut out = Vec::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut alpha_coeffs: Vec<Rational> = Vec::new();
        for rc in part.split(',') {
            alpha_coeffs.push(parse_single_rational(rc)?);
        }
        out.push(field.element(alpha_coeffs));
    }
    Ok(out)
}

/// Render an algebraic-field polynomial as `[c0] + [c1]*x + [c2]*x^2 + ...`,
/// where each `[ci]` is a comma-separated list of $\alpha$-polynomial rationals.
fn anf_poly_to_string(p: &AnfPoly) -> String {
    let coeffs = p.coeffs();
    if coeffs.is_empty() {
        return "0".to_string();
    }
    let mut parts = Vec::new();
    for (i, c) in coeffs.iter().enumerate() {
        if c.coeffs().is_empty() {
            continue;
        }
        let joined: Vec<String> = c.coeffs().iter().map(|r| r.to_string()).collect();
        let cstr = joined.join(",");
        match i {
            0 => parts.push(format!("[{}]", cstr)),
            1 => parts.push(format!("[{}]*x", cstr)),
            _ => parts.push(format!("[{}]*x^{}", cstr, i)),
        }
    }
    if parts.is_empty() {
        "0".to_string()
    } else {
        parts.join(" + ")
    }
}

// ------------------------------------------------------------------
//  Algebraic number field API
// ------------------------------------------------------------------

/// Create an algebraic number field from its monic minimal polynomial.
///
/// `min_poly` is a string such as `"x^2 - 2"` (variable must be `x`). Returns
/// an opaque handle, or `NULL` on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_field_create(
    min_poly: *const c_char,
    err: *mut c_int,
) -> *mut OcasAlgebraicField {
    crate::error::clear();
    if min_poly.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "minimal polynomial string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let s = unsafe { CStr::from_ptr(min_poly) };
    let s = match s.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(
                OCAS_ERROR_INVALID_ARGUMENT,
                "minimal polynomial is not valid UTF-8",
            );
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match parse_min_poly_str(s) {
        Ok(coeffs) => {
            if coeffs.len() < 2 {
                set(
                    OCAS_ERROR_INVALID_ARGUMENT,
                    "minimal polynomial must have degree at least 1",
                );
                crate::error::write_last_code(err);
                return ptr::null_mut();
            }
            if coeffs.last() != Some(&Rational::new(1, 1)) {
                set(
                    OCAS_ERROR_INVALID_ARGUMENT,
                    "minimal polynomial must be monic",
                );
                crate::error::write_last_code(err);
                return ptr::null_mut();
            }
            let field = AlgebraicNumberField::new(RationalDomain, coeffs);
            crate::error::write_last_code(err);
            field_ptr(Box::new(field))
        }
        Err(msg) => {
            set(OCAS_ERROR_PARSE, &msg);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Free an algebraic number field handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_field_free(field: *mut OcasAlgebraicField) {
    if field.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(field as *mut AlgebraicNumberField));
    }
}

/// Return the extension degree $\deg(m)$, or `0` on a null handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_field_degree(field: *const OcasAlgebraicField) -> usize {
    match field_ref(field) {
        Some(f) => f.extension_degree(),
        None => 0,
    }
}

// ------------------------------------------------------------------
//  Algebraic polynomial API
// ------------------------------------------------------------------

/// Create a polynomial over an algebraic number field from a coefficient-list
/// string (see the [module docs](self) for the format).
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_poly_create(
    field: *const OcasAlgebraicField,
    coeffs: *const c_char,
    err: *mut c_int,
) -> *mut OcasAlgebraicPoly {
    crate::error::clear();
    let f = match field_ref(field) {
        Some(f) => f,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "algebraic field handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    if coeffs.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "coefficient string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let s = unsafe { CStr::from_ptr(coeffs) };
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
    match parse_anf_coeffs(f, s) {
        Ok(cs) => {
            let p = DenseUnivariatePolynomial::from_coeffs(f.clone(), cs);
            crate::error::write_last_code(err);
            poly_ptr(Box::new(p))
        }
        Err(msg) => {
            set(OCAS_ERROR_PARSE, &msg);
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Free an algebraic-field polynomial handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_poly_free(poly: *mut OcasAlgebraicPoly) {
    if poly.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(poly as *mut AnfPoly));
    }
}

/// Return the degree of the polynomial, or `0` for the zero polynomial.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_poly_degree(poly: *const OcasAlgebraicPoly) -> usize {
    match poly_ref(poly) {
        Some(p) => p.degree().unwrap_or(0),
        None => 0,
    }
}

/// Return a heap-allocated string representation of the polynomial.
/// The caller must release it with [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_poly_to_string(
    poly: *const OcasAlgebraicPoly,
    err: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    match poly_ref(poly) {
        Some(p) => {
            let s = anf_poly_to_string(p);
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
            set(OCAS_ERROR_NULL_POINTER, "polynomial handle is null");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Factor a polynomial over an algebraic number field (Trager's algorithm).
///
/// On success `out` is filled with a heap-allocated array of factors. The
/// caller must free each factor's `poly` via [`ocas_algebraic_poly_free`]
/// and then the array via [`ocas_algebraic_factor_array_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_poly_factor(
    poly: *const OcasAlgebraicPoly,
    out: *mut OcasAlgebraicFactorArray,
    err: *mut c_int,
) -> c_int {
    crate::error::clear();
    if out.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "output array pointer is null");
        crate::error::write_last_code(err);
        return OCAS_ERROR_NULL_POINTER;
    }
    let p = match poly_ref(poly) {
        Some(p) => p,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "polynomial handle is null");
            crate::error::write_last_code(err);
            return OCAS_ERROR_NULL_POINTER;
        }
    };
    let factors = p.factor();
    let mut array: Vec<OcasAlgebraicFactor> = Vec::with_capacity(factors.len());
    for (poly, multiplicity) in factors {
        array.push(OcasAlgebraicFactor {
            poly: poly_ptr(Box::new(poly)) as *mut c_void,
            multiplicity,
        });
    }
    let len = array.len();
    let factors_ptr = array.as_mut_ptr();
    std::mem::forget(array);
    unsafe {
        ptr::write(
            out,
            OcasAlgebraicFactorArray {
                factors: factors_ptr,
                len,
            },
        );
    }
    crate::error::write_last_code(err);
    OCAS_OK
}

/// Free a factor array returned by [`ocas_algebraic_poly_factor`].
///
/// This releases only the array storage; each factor's `poly` handle must be
/// freed separately via [`ocas_algebraic_poly_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_algebraic_factor_array_free(arr: *mut OcasAlgebraicFactorArray) {
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
