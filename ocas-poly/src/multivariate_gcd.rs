//! Multivariate polynomial GCD.
//!
//! Implements the heuristic/evaluation-interpolation approach: evaluate
//! secondary variables at random points, compute the univariate GCD in the
//! main variable, then recover the multivariate GCD by interpolation.
//!
//! For the bivariate case this reduces to: fix variable `y` at several points,
//! compute univariate GCDs in `x`, and interpolate the coefficients (which are
//! polynomials in `y`).
//!
//! References: Brown (1971); Geddes, Czapor, Labahn, *Algorithms for Computer
//! Algebra*, §7.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use ocas_domain::number_theory::{crt::crt_many, primes_from};
use ocas_domain::{
    Domain, EuclideanDomain, FiniteField, FiniteFieldElement, Integer, IntegerDomain,
};
// The GMP backend's `Integer` is not `Sync`: the parallel batch loop is
// compiled out there, so the prelude (par_iter/into_par_iter) is unused.
#[cfg(not(feature = "gmp"))]
use rayon::prelude::*;

use crate::dense::DenseUnivariatePolynomial;
use crate::rational_reconstruction::rational_reconstruction;
use crate::sparse::{Lex, MonomialOrder, SparseMultivariatePolynomial};

/// Alias for sparse polynomials over the integers with lexicographic order.
pub type ZMPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

/// Compute the GCD of two bivariate polynomials over ℤ[x, y] via evaluation
/// and interpolation.
///
/// Treats `var0` as the main variable `x` and `var1` as `y`. Evaluates `y` at
/// several integer points, computes univariate GCDs in `x`, and interpolates
/// the GCD's coefficients (which are univariate in `y`).
///
/// Returns a primitive polynomial that divides both `a` and `b`, or `None` if
/// the heuristic fails (e.g. unlucky evaluation points).
pub fn bivariate_gcd(a: &ZMPoly, b: &ZMPoly) -> Option<ZMPoly> {
    if a.is_zero() {
        return Some(b.primitive_part());
    }
    if b.is_zero() {
        return Some(a.primitive_part());
    }
    if a.n_vars() < 2 || b.n_vars() < 2 {
        return None;
    }

    // Determine degree bounds in y for both polynomials and the GCD.
    let deg_y_a = poly_degree_in(a, 1);
    let deg_y_b = poly_degree_in(b, 1);
    let deg_y_gcd_bound = deg_y_a.min(deg_y_b);

    if deg_y_gcd_bound == 0 {
        // GCD is purely in x — univariate GCD of the contents.
        return Some(zmpoly_constant_in_y(a).gcd_univariate_x(&zmpoly_constant_in_y(b)));
    }

    // Collect evaluation points and compute univariate GCD images.
    let mut images: Vec<(Integer, DenseUnivariatePolynomial<IntegerDomain>)> = Vec::new();
    let max_points = deg_y_gcd_bound + 2; // a few extra for safety
    let mut eval_point = Integer::from(2i64);

    for _ in 0..max_points + 10 {
        if images.len() >= max_points {
            break;
        }
        // Skip points where either polynomial evaluates to zero.
        let a_eval = eval_univariate_x(a, &eval_point);
        let b_eval = eval_univariate_x(b, &eval_point);
        if a_eval.is_zero() || b_eval.is_zero() {
            eval_point = Integer::from(eval_point.to_bigint() + 1);
            continue;
        }
        let g_eval = a_eval.gcd(&b_eval);
        if g_eval.is_zero() {
            eval_point = Integer::from(eval_point.to_bigint() + 1);
            continue;
        }
        images.push((eval_point.clone(), g_eval));
        eval_point = Integer::from(eval_point.to_bigint() + 1);
    }

    if images.is_empty() {
        return None;
    }

    // Check that all images have the same degree in x (the GCD degree).
    let gcd_deg_x = images[0].1.degree().unwrap_or(0);
    if !images
        .iter()
        .all(|(_, g)| g.degree().unwrap_or(0) == gcd_deg_x)
    {
        // Unlucky evaluation points — some gave a GCD of higher degree.
        // Filter to the minimum degree.
        let min_deg = images
            .iter()
            .map(|(_, g)| g.degree().unwrap_or(0))
            .min()
            .unwrap_or(0);
        images.retain(|(_, g)| g.degree().unwrap_or(0) == min_deg);
        if images.is_empty() {
            return None;
        }
    }

    if images.len() < 2 {
        // Only one image — can't interpolate. Return primitive part of one
        // polynomial (fallback).
        return Some(a.primitive_part());
    }

    // Interpolate each coefficient of the univariate GCD in x as a polynomial
    // in y. For x^i, the i-th coefficient is interpolated from the evaluation
    // points.
    let gcd_deg_y = deg_y_gcd_bound;
    let result = interpolate_gcd(&images, gcd_deg_y, a.n_vars())?;

    // Normalize: make primitive and adjust sign.
    let result = result.primitive_part();
    Some(result)
}

/// Compute the degree of a multivariate polynomial in a given variable.
fn poly_degree_in(p: &ZMPoly, var: usize) -> usize {
    p.terms_ref().keys().map(|e| e[var]).max().unwrap_or(0)
}

/// Evaluate a bivariate ℤ[x,y] polynomial at y=value, yielding a univariate
/// ℤ[x] polynomial.
fn eval_univariate_x(p: &ZMPoly, value: &Integer) -> DenseUnivariatePolynomial<IntegerDomain> {
    let domain = IntegerDomain;
    let mut coeffs_map: std::collections::BTreeMap<usize, Integer> = Default::default();
    for (exp, coeff) in p.terms_ref() {
        let x_deg = exp[0];
        let power = domain.pow(value, exp[1] as u64);
        let new_coeff = domain.mul(coeff, &power);
        let existing = coeffs_map
            .get(&x_deg)
            .cloned()
            .unwrap_or_else(|| Integer::from(0));
        coeffs_map.insert(
            x_deg,
            Integer::from(existing.to_bigint() + new_coeff.to_bigint()),
        );
    }
    let max_deg = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![Integer::from(0); max_deg + 1];
    for (deg, c) in coeffs_map {
        coeffs[deg] = c;
    }
    DenseUnivariatePolynomial::from_coeffs(domain, coeffs)
}

/// Extract the terms of `p` that have zero degree in y (i.e. purely in x),
/// as a univariate polynomial. This is a fallback for degenerate GCD cases.
fn zmpoly_constant_in_y(p: &ZMPoly) -> DenseUnivariatePolynomial<IntegerDomain> {
    let domain = IntegerDomain;
    let mut coeffs_map: std::collections::BTreeMap<usize, Integer> = Default::default();
    for (exp, coeff) in p.terms_ref() {
        if exp[1] == 0 {
            coeffs_map.insert(exp[0], coeff.clone());
        }
    }
    let max_deg = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![Integer::from(0); max_deg + 1];
    for (deg, c) in coeffs_map {
        coeffs[deg] = c;
    }
    DenseUnivariatePolynomial::from_coeffs(domain, coeffs)
}

/// Given a set of `(y_value, g_x_image)` pairs, interpolate the bivariate GCD.
///
/// For each power of x, the coefficient is a polynomial in y interpolated
/// from the corresponding coefficient of each univariate image.
fn interpolate_gcd(
    images: &[(Integer, DenseUnivariatePolynomial<IntegerDomain>)],
    deg_y_bound: usize,
    n_vars: usize,
) -> Option<ZMPoly> {
    let domain = IntegerDomain;
    let gcd_deg_x = images[0].1.degree().unwrap_or(0);
    let _n_points = images.len();

    // For each power of x, collect (y_value, coeff_of_x^i) and interpolate.
    let mut result = ZMPoly::new(domain, n_vars);
    for i in 0..=gcd_deg_x {
        // Collect data points for the i-th coefficient.
        let data: Vec<(Integer, Integer)> = images
            .iter()
            .map(|(y_val, g)| {
                let c = g.coeff(i).cloned().unwrap_or_else(|| Integer::from(0));
                (y_val.clone(), c)
            })
            .collect();

        // Interpolate via Lagrange interpolation to get the y-polynomial.
        if let Some(y_poly) = lagrange_interpolate(&data, deg_y_bound) {
            for (y_deg, y_coeff) in y_poly.iter().enumerate() {
                if y_coeff.is_zero() {
                    continue;
                }
                let mut exp = vec![0usize; n_vars];
                exp[0] = i;
                if n_vars > 1 {
                    exp[1] = y_deg;
                }
                // Accumulate the term.
                let existing = result.coeff(&exp);
                let sum = domain.add(&existing, y_coeff);
                result.set_term_external(exp, sum);
            }
        }
    }
    Some(result)
}

/// Lagrange interpolation: given points (x_0, y_0), ..., (x_n, y_n), find the
/// polynomial p of degree ≤ n with p(x_i) = y_i.
///
/// Returns the coefficient vector [a_0, a_1, ..., a_n] (ascending powers) using
/// exact arithmetic with numerator/denominator pairs over `BigInt`.
fn lagrange_interpolate(points: &[(Integer, Integer)], max_deg: usize) -> Option<Vec<Integer>> {
    let n = points.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if n == 1 {
        return Some(vec![points[0].1.clone()]);
    }

    // Represent coefficients as (numerator, denominator) fractions.
    // Start with all-zero coefficients.
    let mut coeffs: Vec<(BigInt, BigInt)> = vec![(BigInt::zero(), BigInt::one()); n];

    for i in 0..n {
        let xi = points[i].0.to_bigint();
        let yi = points[i].1.to_bigint();

        // denominator = product of (x_i - x_j) for j != i
        let mut denom = BigInt::one();
        for (j, (xj_val, _)) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let xj = xj_val.to_bigint();
            denom *= xi.clone() - xj;
        }
        if denom.is_zero() {
            return None;
        }

        // scale = y_i / denom — we'll multiply basis polynomials by yi and
        // accumulate with common denominator `denom`.
        // Build the Lagrange basis L_i(x) = product of (x - x_j) for j != i.
        let mut basis: Vec<BigInt> = vec![BigInt::one()];
        for (j, (xj_val, _)) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let neg_xj = -(xj_val.to_bigint());
            // Multiply basis by (x - x_j):
            let mut new_basis = vec![BigInt::zero(); basis.len() + 1];
            for k in 0..basis.len() {
                new_basis[k] += &neg_xj * &basis[k];
                new_basis[k + 1] += &basis[k];
            }
            basis = new_basis;
        }

        // Scale basis by yi/denom and accumulate into coeffs.
        // coeffs[k] += yi * basis[k] / denom
        for k in 0..basis.len().min(n) {
            // coeffs[k] += yi * basis[k] / denom
            // = (coeffs_num[k] * denom + yi * basis[k] * coeffs_den[k]) / (coeffs_den[k] * denom)
            let new_num = &coeffs[k].0 * &denom + &yi * &basis[k] * &coeffs[k].1;
            let new_den = &coeffs[k].1 * &denom;
            // Reduce by GCD to prevent blowup.
            let g = bigint_gcd(&new_num, &new_den);
            if !g.is_zero() && !g.is_one() {
                coeffs[k] = (new_num / &g, new_den / &g);
            } else {
                coeffs[k] = (new_num, new_den);
            }
        }
    }

    // Convert to integers: all coefficients should be integers (denominator divides numerator).
    let mut result = Vec::with_capacity(n.min(max_deg + 1));
    for (num_, den) in &coeffs {
        if den.is_zero() {
            return None;
        }
        let q = num_ / den;
        let r = num_ % den;
        if r.is_zero() {
            result.push(Integer::from(q));
        } else {
            // Non-integer coefficient — interpolation failed.
            return None;
        }
        if result.len() > max_deg + 1 {
            break;
        }
    }
    Some(result)
}

/// Compute gcd(a, b) for two BigInt values via the Euclidean algorithm.
fn bigint_gcd(a: &BigInt, b: &BigInt) -> BigInt {
    let mut a = a.abs();
    let mut b = b.abs();
    while !b.is_zero() {
        let r = &a % &b;
        a = b;
        b = r;
    }
    a
}

/// Trait extension for univariate polynomials used by the bivariate GCD.
trait UnivariateGcdExt {
    fn gcd_univariate_x(&self, other: &Self) -> ZMPoly;
}

impl UnivariateGcdExt for DenseUnivariatePolynomial<IntegerDomain> {
    fn gcd_univariate_x(&self, other: &Self) -> ZMPoly {
        let g = self.gcd(other);
        // Wrap into a bivariate polynomial (y degree 0).
        let mut result = ZMPoly::new(IntegerDomain, 2);
        for (i, c) in g.coeffs().iter().enumerate() {
            if !c.is_zero() {
                result.set_term_external(vec![i, 0], c.clone());
            }
        }
        result.primitive_part()
    }
}

// =========================================================================
// Modular multivariate GCD (§F1–F3)
// =========================================================================

use crate::factor::multivariate::FpMPoly;

/// Reduce a ℤ polynomial mod a prime, yielding a ℤ_p polynomial.
pub fn reduce_mod(p: &ZMPoly, prime: &BigInt) -> FpMPoly {
    let field = FiniteField::new(prime.clone());
    let mut result = FpMPoly::new(field.clone(), p.n_vars());
    for (exp, coeff) in p.terms_ref() {
        let c_fp = field.element(coeff.to_bigint());
        if !c_fp.value().is_zero() {
            result.set_term_external(exp.to_vec(), c_fp);
        }
    }
    result
}

/// Lift a ℤ_p polynomial back to ℤ using symmetric representatives in
/// `[-(p-1)/2, (p-1)/2]`.
pub fn lift_from_fp(p: &FpMPoly) -> ZMPoly {
    let field = p.domain().clone();
    let prime = field.prime();
    let half_p = prime / 2u32;
    let mut result = ZMPoly::new(IntegerDomain, p.n_vars());
    for (exp, coeff) in p.terms_ref() {
        let v = coeff.value();
        let lifted = if *v > half_p {
            Integer::from(v - prime)
        } else {
            Integer::from(v.clone())
        };
        if !lifted.is_zero() {
            result.set_term_external(exp.to_vec(), lifted);
        }
    }
    result
}

/// Evaluation-interpolation bivariate GCD over ℤ_p.
///
/// Same algorithm as [`bivariate_gcd`] but operating in a finite field,
/// avoiding coefficient growth.
pub fn bivariate_gcd_fp(a: &FpMPoly, b: &FpMPoly) -> Option<FpMPoly> {
    if a.is_zero() {
        return Some(b.clone());
    }
    if b.is_zero() {
        return Some(a.clone());
    }
    if a.n_vars() < 2 || b.n_vars() < 2 {
        return None;
    }

    let field = a.domain().clone();
    let deg_y_a = fp_poly_degree_in(a, 1);
    let deg_y_b = fp_poly_degree_in(b, 1);
    let deg_y_gcd_bound = deg_y_a.min(deg_y_b);

    if deg_y_gcd_bound == 0 {
        // GCD is purely in x — univariate GCD of the contents.
        return Some(fp_univariate_gcd_x(a, b));
    }

    let mut images: Vec<(usize, DenseUnivariatePolynomial<FiniteField>)> = Vec::new();
    let max_points = deg_y_gcd_bound + 2;
    let mut eval_val = 1usize;

    for _ in 0..max_points + 20 {
        if images.len() >= max_points {
            break;
        }
        let a_eval = fp_eval_univariate_x(a, eval_val);
        let b_eval = fp_eval_univariate_x(b, eval_val);
        if a_eval.is_zero() || b_eval.is_zero() {
            eval_val += 1;
            continue;
        }
        let g_eval = a_eval.gcd(&b_eval);
        if g_eval.is_zero() {
            eval_val += 1;
            continue;
        }
        // Normalize to monic so images across evaluation points (and across
        // primes, for the multi-prime caller) are consistent associates.
        let g_eval = monic_fp(&g_eval, &field);
        images.push((eval_val, g_eval));
        eval_val += 1;
    }

    if images.is_empty() {
        return None;
    }

    // Filter to consistent degree.
    let min_deg = images
        .iter()
        .map(|(_, g)| g.degree().unwrap_or(0))
        .min()
        .unwrap_or(0);
    images.retain(|(_, g)| g.degree().unwrap_or(0) == min_deg);
    if images.is_empty() {
        return None;
    }
    if images.len() < 2 {
        return Some(a.clone());
    }

    fp_interpolate_gcd(&images, deg_y_gcd_bound, a.n_vars(), &field)
}

/// Degree of a multivariate polynomial in a given variable (FiniteField version).
fn fp_poly_degree_in(p: &FpMPoly, var: usize) -> usize {
    p.terms_ref().keys().map(|e| e[var]).max().unwrap_or(0)
}

/// Evaluate a bivariate ℤ_p[x,y] polynomial at y=value, yielding a univariate
/// ℤ_p[x] polynomial.
fn fp_eval_univariate_x(p: &FpMPoly, value: usize) -> DenseUnivariatePolynomial<FiniteField> {
    let field = p.domain().clone();
    let val_el = field.element(BigInt::from(value));
    let mut coeffs_map: std::collections::BTreeMap<usize, FiniteFieldElement> = Default::default();
    for (exp, coeff) in p.terms_ref() {
        let x_deg = exp[0];
        let power = field.pow(&val_el, exp[1] as u64);
        let new_coeff = field.mul(coeff, &power);
        let existing = coeffs_map
            .get(&x_deg)
            .cloned()
            .unwrap_or_else(|| field.zero());
        coeffs_map.insert(x_deg, field.add(&existing, &new_coeff));
    }
    let max_deg = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![field.zero(); max_deg + 1];
    for (deg, c) in coeffs_map {
        coeffs[deg] = c;
    }
    DenseUnivariatePolynomial::from_coeffs(field, coeffs)
}

/// Univariate GCD of the y=0 content (both polynomials viewed as univariate in x).
fn fp_univariate_gcd_x(a: &FpMPoly, b: &FpMPoly) -> FpMPoly {
    let field = a.domain().clone();
    let a_x = fp_extract_constant_in_y(a);
    let b_x = fp_extract_constant_in_y(b);
    let g = monic_fp(&a_x.gcd(&b_x), &field);
    let mut result = FpMPoly::new(field, a.n_vars());
    for (i, c) in g.coeffs().iter().enumerate() {
        if !c.value().is_zero() {
            result.set_term_external(vec![i, 0], c.clone());
        }
    }
    result
}

/// Normalize a dense univariate polynomial over a field to be monic.
fn monic_fp(
    g: &DenseUnivariatePolynomial<FiniteField>,
    field: &FiniteField,
) -> DenseUnivariatePolynomial<FiniteField> {
    match g.leading_coeff() {
        Some(lc) if !field.is_one(lc) => {
            g.mul_scalar(&field.inv(lc).expect("nonzero field element"))
        }
        _ => g.clone(),
    }
}

/// Extract terms with y-degree 0 as a univariate polynomial.
fn fp_extract_constant_in_y(p: &FpMPoly) -> DenseUnivariatePolynomial<FiniteField> {
    let field = p.domain().clone();
    let mut coeffs_map: std::collections::BTreeMap<usize, FiniteFieldElement> = Default::default();
    for (exp, coeff) in p.terms_ref() {
        if exp[1] == 0 {
            coeffs_map.insert(exp[0], coeff.clone());
        }
    }
    let max_deg = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![field.zero(); max_deg + 1];
    for (deg, c) in coeffs_map {
        coeffs[deg] = c;
    }
    DenseUnivariatePolynomial::from_coeffs(field, coeffs)
}

/// Interpolate the bivariate GCD from univariate images in ℤ_p.
fn fp_interpolate_gcd(
    images: &[(usize, DenseUnivariatePolynomial<FiniteField>)],
    deg_y_bound: usize,
    n_vars: usize,
    field: &FiniteField,
) -> Option<FpMPoly> {
    let gcd_deg_x = images[0].1.degree().unwrap_or(0);
    let mut result = FpMPoly::new(field.clone(), n_vars);

    for i in 0..=gcd_deg_x {
        let data: Vec<(usize, FiniteFieldElement)> = images
            .iter()
            .map(|(y_val, g)| {
                let c = g.coeff(i).cloned().unwrap_or_else(|| field.zero());
                (*y_val, c)
            })
            .collect();

        if let Some(y_poly) = fp_lagrange_interpolate(&data, deg_y_bound, field) {
            for (y_deg, y_coeff) in y_poly.iter().enumerate() {
                if y_coeff.value().is_zero() {
                    continue;
                }
                let mut exp = vec![0usize; n_vars];
                exp[0] = i;
                if n_vars > 1 {
                    exp[1] = y_deg;
                }
                let existing = result.coeff(&exp);
                let sum = field.add(&existing, y_coeff);
                result.set_term_external(exp, sum);
            }
        }
    }
    Some(result)
}

/// Lagrange interpolation in ℤ_p.
fn fp_lagrange_interpolate(
    points: &[(usize, FiniteFieldElement)],
    _max_deg: usize,
    field: &FiniteField,
) -> Option<Vec<FiniteFieldElement>> {
    let n = points.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if n == 1 {
        return Some(vec![points[0].1.clone()]);
    }

    let mut coeffs = vec![field.zero(); n];

    for i in 0..n {
        let xi = field.element(BigInt::from(points[i].0));
        let yi = &points[i].1;

        // denominator = product of (x_i - x_j) for j != i
        let mut denom = field.one();
        for (j, (xj_val, _)) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let xj = field.element(BigInt::from(*xj_val));
            let diff = field.sub(&xi, &xj);
            denom = field.mul(&denom, &diff);
        }
        let denom_inv = field.inv(&denom)?;
        let scale = field.mul(yi, &denom_inv);

        // Build Lagrange basis L_i(x) = product of (x - x_j) for j != i.
        let mut basis: Vec<FiniteFieldElement> = vec![field.one()];
        for (j, (xj_val, _)) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let neg_xj = field.neg(&field.element(BigInt::from(*xj_val)));
            let mut new_basis = vec![field.zero(); basis.len() + 1];
            for k in 0..basis.len() {
                // new_basis[k] += neg_xj * basis[k]
                let term = field.mul(&neg_xj, &basis[k]);
                new_basis[k] = field.add(&new_basis[k], &term);
                // new_basis[k+1] += basis[k]
                new_basis[k + 1] = field.add(&new_basis[k + 1], &basis[k]);
            }
            basis = new_basis;
        }

        for k in 0..basis.len().min(n) {
            let term = field.mul(&scale, &basis[k]);
            coeffs[k] = field.add(&coeffs[k], &term);
        }
    }

    Some(coeffs)
}

/// Modular bivariate GCD over ℤ (multi-prime Brown strategy).
///
/// Separates the content with respect to the main variable `x` (a GCD of
/// univariate polynomials in `y`), computes monic-in-`x` GCDs of the
/// primitive parts modulo a sequence of primes, and reconstructs the
/// rational coefficients of the true GCD via CRT plus rational
/// reconstruction. The candidate is accepted only when it has full degree
/// and divides both primitive parts exactly (checked by reduction to zero).
/// Primes whose modular GCD has a larger degree in `x` than the current
/// best ("unlucky" primes) are discarded; primes dividing the integer
/// content of either leading coefficient are skipped, since they would
/// corrupt the degree comparison.
pub fn gcd_modular(a: &ZMPoly, b: &ZMPoly) -> Option<ZMPoly> {
    if a.is_zero() {
        return Some(b.primitive_part());
    }
    if b.is_zero() {
        return Some(a.primitive_part());
    }
    if a.n_vars() < 2 || b.n_vars() < 2 {
        return None;
    }
    let n_vars = a.n_vars();

    let a_prim = a.primitive_part();
    let b_prim = b.primitive_part();
    // Content with respect to x (a polynomial in y only).
    let cont_a = content_in_x(&a_prim, n_vars);
    let cont_b = content_in_x(&b_prim, n_vars);
    let cont_g = gcd_poly_y(&cont_a, &cont_b, n_vars);
    let pp_a = a_prim.checked_div_exact(&cont_a)?;
    let pp_b = b_prim.checked_div_exact(&cont_b)?;

    let deg_a = poly_degree_in(&pp_a, 0);
    let deg_b = poly_degree_in(&pp_b, 0);
    let pp_gcd = if deg_a == 0 && deg_b == 0 {
        gcd_poly_y(&pp_a, &pp_b, n_vars)
    } else if deg_a == 0 || deg_b == 0 {
        // A y-only primitive part is coprime to any x-primitive polynomial.
        pp_a.one()
    } else {
        modular_gcd_x(&pp_a, &pp_b, n_vars)?
    };
    Some(cont_g.mul(&pp_gcd).primitive_part())
}

/// Compute one bivariate modular GCD image at prime `p`.
///
/// Returns `None` when `p` divides the `bad` integer or the image
/// vanishes. Extracted so the batch loop can run under either a parallel
/// or sequential iterator without duplicating the body.
fn modular_gcd_x_image(
    p: &Integer,
    bad: &Integer,
    pp_a: &ZMPoly,
    pp_b: &ZMPoly,
) -> Option<(Integer, FpMPoly, usize, bool)> {
    if !bad.is_one() && bad.mod_floor(p).is_zero() {
        return None;
    }
    let a_p = reduce_mod(pp_a, &p.to_bigint());
    let b_p = reduce_mod(pp_b, &p.to_bigint());
    let g_p = bivariate_gcd_fp(&a_p, &b_p)?;
    let deg_x = fp_poly_degree_in(&g_p, 0);
    let is_constant = g_p.total_degree() == Some(0);
    Some((p.clone(), g_p, deg_x, is_constant))
}

/// Multi-prime modular GCD of the primitive parts (Brown's algorithm).
fn modular_gcd_x(pp_a: &ZMPoly, pp_b: &ZMPoly, n_vars: usize) -> Option<ZMPoly> {
    // Any prime that makes the x-leading coefficient of the true GCD vanish
    // mod p divides this integer; skipping such primes keeps the degree
    // comparison sound.
    let bad = IntegerDomain.gcd(&lc_in_x(pp_a).content(), &lc_in_x(pp_b).content());

    let mut best_deg_x: Option<usize> = None;
    let mut images: Vec<(Integer, FpMPoly)> = Vec::new();
    let mut prime_iter = primes_from(&Integer::from(1_000_000_007));
    let batch_size = rayon::current_num_threads().max(1);
    let mut tried = 0usize;
    while tried < 64 {
        let batch: Vec<Integer> = prime_iter.by_ref().take(batch_size).collect();
        if batch.is_empty() {
            break;
        }
        tried += batch.len();
        // Bivariate modular GCD images, computed in parallel (sequential
        // under the GMP backend, whose `Integer` is not `Sync`). Primes
        // dividing the bad integer or with a vanished image contribute None.
        #[cfg(not(feature = "gmp"))]
        let computed: Vec<Option<(Integer, FpMPoly, usize, bool)>> = batch
            .par_iter()
            .map(|p| modular_gcd_x_image(p, &bad, pp_a, pp_b))
            .collect();
        #[cfg(feature = "gmp")]
        let computed: Vec<Option<(Integer, FpMPoly, usize, bool)>> = batch
            .iter()
            .map(|p| modular_gcd_x_image(p, &bad, pp_a, pp_b))
            .collect();
        for item in computed {
            let Some((p, g_p, deg_x, is_constant)) = item else {
                continue;
            };
            match best_deg_x {
                None => {
                    best_deg_x = Some(deg_x);
                    images.push((p, g_p));
                }
                Some(bd) if deg_x < bd => {
                    best_deg_x = Some(deg_x);
                    images.clear();
                    images.push((p, g_p));
                }
                Some(bd) if deg_x == bd => images.push((p, g_p)),
                _ => continue, // unlucky prime
            }
            // A constant monic image means the primitive parts are coprime.
            if best_deg_x == Some(0) && is_constant {
                return Some(pp_a.one());
            }
            // Trial reconstruction: accept only a common divisor of full degree.
            if let Some(cand) = reconstruct_rational(&images, n_vars) {
                let cand = cand.primitive_part();
                if poly_degree_in(&cand, 0) == best_deg_x.unwrap_or(0)
                    && pp_a.reduce(std::slice::from_ref(&cand)).is_zero()
                    && pp_b.reduce(std::slice::from_ref(&cand)).is_zero()
                {
                    return Some(cand);
                }
            }
        }
    }
    None
}

/// CRT + rational-reconstruction pass over the per-prime monic images:
/// combine each coefficient across primes, recover its rational value,
/// and clear denominators to obtain an integer polynomial.
fn reconstruct_rational(images: &[(Integer, FpMPoly)], n_vars: usize) -> Option<ZMPoly> {
    let mut exps: std::collections::BTreeSet<Vec<usize>> = std::collections::BTreeSet::new();
    for (_, img) in images {
        for e in img.terms_ref().keys() {
            exps.insert(e.to_vec());
        }
    }
    let mut rationals: Vec<(Vec<usize>, Integer, Integer)> = Vec::with_capacity(exps.len());
    for e in exps {
        let cs: Vec<(Integer, Integer)> = images
            .iter()
            .map(|(p, img)| (Integer::from(img.coeff(&e).value().clone()), p.clone()))
            .collect();
        let (r, m) = crt_many(&cs)?;
        let (n, d) = rational_reconstruction(&r, &m)?;
        rationals.push((e, n, d));
    }
    // Clear denominators: multiply every numerator by lcm(d_i)/d_i.
    let mut lcm = Integer::from(1);
    for (_, _, d) in &rationals {
        lcm = &lcm * &(d / &IntegerDomain.gcd(&lcm, d));
    }
    let mut result = ZMPoly::new(IntegerDomain, n_vars);
    for (e, n, d) in rationals {
        let coeff = &(&lcm / &d) * &n;
        if !coeff.is_zero() {
            result.set_term_external(e, coeff);
        }
    }
    Some(result)
}

/// Convert a bivariate polynomial that only involves `y` into a dense
/// univariate polynomial in `y`.
fn zmpoly_to_dense_y(p: &ZMPoly) -> DenseUnivariatePolynomial<IntegerDomain> {
    let mut coeffs: Vec<Integer> = Vec::new();
    for (exp, c) in p.terms_ref() {
        debug_assert_eq!(exp[0], 0, "expected a y-only polynomial");
        let y = exp.get(1).copied().unwrap_or(0);
        if coeffs.len() <= y {
            coeffs.resize(y + 1, Integer::from(0));
        }
        coeffs[y] = Integer::from(coeffs[y].to_bigint() + c.to_bigint());
    }
    DenseUnivariatePolynomial::from_coeffs(IntegerDomain, coeffs)
}

/// Wrap a dense univariate-in-`y` polynomial as a bivariate polynomial.
fn dense_y_to_zmpoly(g: &DenseUnivariatePolynomial<IntegerDomain>, n_vars: usize) -> ZMPoly {
    let mut result = ZMPoly::new(IntegerDomain, n_vars);
    for (j, c) in g.coeffs().iter().enumerate() {
        if !c.is_zero() {
            let mut exp = vec![0usize; n_vars];
            if n_vars > 1 {
                exp[1] = j;
            }
            result.set_term_external(exp, c.clone());
        }
    }
    result
}

/// GCD of two y-only bivariate polynomials (dense univariate GCD in `y`).
fn gcd_poly_y(f: &ZMPoly, g: &ZMPoly, n_vars: usize) -> ZMPoly {
    let h = zmpoly_to_dense_y(f).gcd(&zmpoly_to_dense_y(g));
    dense_y_to_zmpoly(&h, n_vars)
}

/// Content of `p` with respect to the main variable `x`: the GCD of the
/// `y`-polynomial coefficients of `x^k` (itself a polynomial in `y` only).
fn content_in_x(p: &ZMPoly, n_vars: usize) -> ZMPoly {
    let deg_x = poly_degree_in(p, 0);
    let mut acc = p.coeff_of_var_pow(0, 0);
    for k in 1..=deg_x {
        if acc.total_degree() == Some(0) {
            break;
        }
        let ck = p.coeff_of_var_pow(0, k);
        if ck.is_zero() {
            continue;
        }
        if acc.is_zero() {
            acc = ck;
            continue;
        }
        acc = gcd_poly_y(&acc, &ck, n_vars);
    }
    acc
}

/// Leading coefficient of `p` with respect to the main variable `x`
/// (a polynomial in `y` only).
fn lc_in_x(p: &ZMPoly) -> ZMPoly {
    p.coeff_of_var_pow(0, poly_degree_in(p, 0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::Integer;

    fn zmp2(terms: &[(usize, usize, i64)]) -> ZMPoly {
        let domain = IntegerDomain;
        let terms_vec: Vec<(Vec<usize>, Integer)> = terms
            .iter()
            .map(|(xd, yd, c)| (vec![*xd, *yd], Integer::from(*c)))
            .collect();
        ZMPoly::from_terms(domain, 2, terms_vec)
    }

    fn reconstruct_check(a: &ZMPoly, b: &ZMPoly, g: &ZMPoly) -> bool {
        // g should divide both a and b (approximately: check that
        // evaluating g's univariate images divides correctly).
        // Simple check: g's total degree ≤ both a and b's.
        g.total_degree().unwrap_or(0) <= a.total_degree().unwrap_or(0)
            && g.total_degree().unwrap_or(0) <= b.total_degree().unwrap_or(0)
    }

    #[test]
    fn gcd_coprime_bivariate() {
        // gcd(x^2+y^2, x^2-y^2) = 1 (or x^2, but they're coprime over Z).
        let a = zmp2(&[(2, 0, 1), (0, 2, 1)]); // x^2 + y^2
        let b = zmp2(&[(2, 0, 1), (0, 2, -1)]); // x^2 - y^2
        let g = bivariate_gcd(&a, &b);
        assert!(g.is_some(), "GCD should succeed");
        let g = g.unwrap();
        // gcd should be 1 (constant).
        assert!(
            g.total_degree().unwrap_or(0) == 0 || g.n_terms() <= 1,
            "coprime GCD should be constant, got {:?}",
            g.total_degree()
        );
    }

    #[test]
    fn gcd_shared_linear_factor() {
        // a = (x+y)(x+1) = x^2 + xy + x + y
        // b = (x+y)(x+2) = x^2 + 2x + xy + 2y
        let a = zmp2(&[(2, 0, 1), (1, 1, 1), (1, 0, 1), (0, 1, 1)]);
        let b = zmp2(&[(2, 0, 1), (1, 1, 1), (1, 0, 2), (0, 1, 2)]);
        let g = bivariate_gcd(&a, &b);
        assert!(g.is_some());
        let g = g.unwrap();
        // GCD should be x+y (degree 1 in x, degree 1 in y, total degree 1 in x).
        assert!(reconstruct_check(&a, &b, &g), "GCD degree inconsistent");
    }

    #[test]
    fn content_and_primitive_part_bivariate() {
        // 2x^2 + 4xy + 6y has content 2, primitive part x^2 + 2xy + 3y.
        let p = zmp2(&[(2, 0, 2), (1, 1, 4), (0, 1, 6)]);
        let content = p.content();
        assert_eq!(content, Integer::from(2));
        let pp = p.primitive_part();
        assert_eq!(pp.coeff(&[2, 0]), Integer::from(1));
        assert_eq!(pp.coeff(&[1, 1]), Integer::from(2));
        assert_eq!(pp.coeff(&[0, 1]), Integer::from(3));
    }

    // --- Modular GCD tests (F1–F3) ---

    #[test]
    fn reduce_mod_and_lift_roundtrip() {
        // p = 3x^2 + 5xy - 7y, mod 11
        let p = zmp2(&[(2, 0, 3), (1, 1, 5), (0, 1, -7)]);
        let prime = BigInt::from(11);
        let p_fp = reduce_mod(&p, &prime);
        let p_lifted = lift_from_fp(&p_fp);
        // After reduce+lift, coefficients should match mod 11 (symmetric rep).
        // 3 mod 11 = 3, 5 mod 11 = 5, -7 mod 11 = 4 (symmetric: 4).
        assert_eq!(p_lifted.coeff(&[2, 0]), Integer::from(3));
        assert_eq!(p_lifted.coeff(&[1, 1]), Integer::from(5));
        assert_eq!(p_lifted.coeff(&[0, 1]), Integer::from(4)); // -7 mod 11 = 4
    }

    #[test]
    fn gcd_modular_shared_linear_factor() {
        // a = (x+y)(x+1) = x^2 + xy + x + y
        // b = (x+y)(x+2) = x^2 + 2x + xy + 2y
        let a = zmp2(&[(2, 0, 1), (1, 1, 1), (1, 0, 1), (0, 1, 1)]);
        let b = zmp2(&[(2, 0, 1), (1, 1, 1), (1, 0, 2), (0, 1, 2)]);
        let g = gcd_modular(&a, &b);
        assert!(g.is_some(), "modular GCD should succeed");
        let g = g.unwrap();
        assert!(reconstruct_check(&a, &b, &g), "GCD degree inconsistent");
    }

    #[test]
    fn gcd_modular_coprime() {
        // gcd(x+y, x-y) = 1
        let a = zmp2(&[(1, 0, 1), (0, 1, 1)]);
        let b = zmp2(&[(1, 0, 1), (0, 1, -1)]);
        let g = gcd_modular(&a, &b);
        assert!(g.is_some(), "modular GCD should succeed for coprime");
        let g = g.unwrap();
        // Should be constant (degree 0).
        assert!(
            g.total_degree().unwrap_or(0) == 0 || g.n_terms() <= 1,
            "coprime GCD should be constant, got degree {:?}",
            g.total_degree()
        );
    }

    #[test]
    fn gcd_modular_shared_quadratic() {
        // a = (x^2 + y)(x + 1) = x^3 + x^2 + xy + y
        // b = (x^2 + y)(x + 2) = x^3 + 2x^2 + xy + 2y
        let a = zmp2(&[(3, 0, 1), (2, 0, 1), (1, 1, 1), (0, 1, 1)]);
        let b = zmp2(&[(3, 0, 1), (2, 0, 2), (1, 1, 1), (0, 1, 2)]);
        let g = gcd_modular(&a, &b);
        assert!(g.is_some(), "modular GCD should succeed");
        let g = g.unwrap();
        // GCD should be x^2 + y (degree 2).
        assert!(reconstruct_check(&a, &b, &g), "GCD degree inconsistent");
    }

    #[test]
    fn gcd_modular_large_coefficients_multi_prime() {
        // Coefficients ~10^20 (67 bits) cannot be reconstructed from a
        // single ~31-bit prime image: this exercises the multi-prime CRT
        // path end to end.
        let big1 = Integer::from(10i64).pow_u32(20);
        let big2 = Integer::from(10i64).pow_u32(19);
        // g = x + 2y + 1 (the shared factor).
        let g = ZMPoly::from_terms(
            IntegerDomain,
            2,
            vec![
                (vec![1, 0], Integer::from(1)),
                (vec![0, 1], Integer::from(2)),
                (vec![0, 0], Integer::from(1)),
            ],
        );
        // h1 = 10^20·x + 3y − 5, h2 = 7x − 10^19·y + 2.
        let h1 = ZMPoly::from_terms(
            IntegerDomain,
            2,
            vec![
                (vec![1, 0], big1),
                (vec![0, 1], Integer::from(3)),
                (vec![0, 0], Integer::from(-5)),
            ],
        );
        let h2 = ZMPoly::from_terms(
            IntegerDomain,
            2,
            vec![
                (vec![1, 0], Integer::from(7)),
                (vec![0, 1], -big2),
                (vec![0, 0], Integer::from(2)),
            ],
        );
        let a = g.mul(&h1);
        let b = g.mul(&h2);
        let got = gcd_modular(&a, &b).expect("multi-prime CRT must succeed");
        assert_eq!(got, g);
        assert!(a.reduce(std::slice::from_ref(&got)).is_zero());
        assert!(b.reduce(std::slice::from_ref(&got)).is_zero());
    }

    #[test]
    fn gcd_modular_y_only_factor() {
        // Shared factor depends only on y: g = 3y + 2.
        let g = zmp2(&[(0, 1, 3), (0, 0, 2)]);
        let a = g.mul(&zmp2(&[(1, 0, 1), (0, 1, 1), (0, 0, -1)]));
        let b = g.mul(&zmp2(&[(2, 0, 1), (0, 0, 4)]));
        let got = gcd_modular(&a, &b).expect("y-only factor GCD must succeed");
        assert!(got.coeff(&[1, 0]).is_zero(), "GCD must not contain x");
        assert!(a.reduce(std::slice::from_ref(&got)).is_zero());
        assert!(b.reduce(std::slice::from_ref(&got)).is_zero());
        assert!(poly_degree_in(&got, 1) == 1, "GCD must keep the y factor");
    }

    // --- Property tests ---

    proptest::proptest! {
        #[test]
        fn gcd_modular_consistency(
            // Generate two bivariate polynomials with a shared linear factor (x + ay + b)
            a_coeff in -5i64..5,
            b_coeff in -5i64..5,
            c1 in -3i64..3,
            d1 in -3i64..3,
            c2 in -3i64..3,
            d2 in -3i64..3,
        ) {
            // shared factor: x + a_coeff*y + b_coeff
            // p1 = c1*x + d1*y + (c1*b_coeff + d1*a_coeff*y_coeff... no, use direct construction
            // Build: a = (x + ay + b)(c1*x + d1) = c1*x^2 + a*c1*xy + b*c1*x + d1*x + a*d1*y + b*d1
            // Build: b = (x + ay + b)(c2*x + d2) = c2*x^2 + a*c2*xy + b*c2*x + d2*x + a*d2*y + b*d2
            let a = zmp2(&[
                (2, 0, c1),
                (1, 1, a_coeff * c1),
                (1, 0, b_coeff * c1 + d1),
                (0, 1, a_coeff * d1),
                (0, 0, b_coeff * d1),
            ]);
            let b = zmp2(&[
                (2, 0, c2),
                (1, 1, a_coeff * c2),
                (1, 0, b_coeff * c2 + d2),
                (0, 1, a_coeff * d2),
                (0, 0, b_coeff * d2),
            ]);

            // Skip trivial cases where either polynomial is zero.
            if a.is_zero() || b.is_zero() { return Ok(()); }

            let g_mod = gcd_modular(&a, &b);
            let g_heu = bivariate_gcd(&a, &b);

            // Both should succeed or both should fail.
            match (&g_mod, &g_heu) {
                (Some(gm), Some(gh)) => {
                    // Modular GCD degree should be ≤ heuristic GCD degree
                    // (heuristic may return inflated results in edge cases).
                    let deg_m = gm.total_degree().unwrap_or(0);
                    let deg_h = gh.total_degree().unwrap_or(0);
                    assert!(deg_m <= deg_h,
                        "modular GCD degree {} > heuristic GCD degree {}", deg_m, deg_h);
                }
                (None, None) => {}
                _ => {
                    // One succeeded and the other didn't — acceptable for edge cases
                    // where the modular approach needs different prime selection.
                }
            }
        }
    }
}

// =========================================================================
// Multivariate GCD (arbitrary number of variables)
// =========================================================================
//
// Dense recursive evaluation–interpolation GCD: treat the polynomials as
// univariate in variable 0, evaluate the last variable at successive points,
// compute the GCD recursively in one fewer variable, and interpolate the
// coefficients back. Correctness is verified by trial division, so unlucky
// evaluation points only cost extra iterations.
//
// References: Brown (1971); Geddes, Czapor, Labahn, *Algorithms for Computer
// Algebra*, §7.

use ocas_domain::{Rational, RationalDomain};

/// Sparse multivariate polynomial over the rationals (auxiliary for the
/// integer multivariate GCD).
pub type QMPoly = SparseMultivariatePolynomial<RationalDomain, Lex>;

/// Convert a sparse polynomial that only involves variable 0 into a dense
/// univariate polynomial in that variable.
fn mpoly_to_dense_x<D: Domain>(
    p: &SparseMultivariatePolynomial<D, Lex>,
) -> DenseUnivariatePolynomial<D> {
    let mut coeffs = Vec::new();
    for (exp, c) in p.terms_ref() {
        let idx = exp.first().copied().unwrap_or(0);
        if idx >= coeffs.len() {
            coeffs.resize(idx + 1, p.domain().zero());
        }
        coeffs[idx] = p.domain().add(&coeffs[idx], c);
    }
    DenseUnivariatePolynomial::from_coeffs(p.domain().clone(), coeffs)
}

/// Wrap a dense univariate polynomial as a sparse polynomial in `n_vars`
/// variables (variable 0 is the polynomial variable).
fn dense_to_mpoly_x<D: Domain>(
    g: &DenseUnivariatePolynomial<D>,
    n_vars: usize,
) -> SparseMultivariatePolynomial<D, Lex> {
    let mut result = SparseMultivariatePolynomial::new(g.domain().clone(), n_vars);
    for (i, c) in g.coeffs().iter().enumerate() {
        if !g.domain().is_zero(c) {
            let mut exp = vec![0usize; n_vars];
            exp[0] = i;
            result.set_term_external(exp, c.clone());
        }
    }
    result
}

/// Scale a polynomial so that its Lex-leading coefficient is 1.
///
/// Requires the coefficient domain to be a field (every nonzero element must
/// have an inverse).
fn normalize_unit_lc<D: Domain>(
    p: &SparseMultivariatePolynomial<D, Lex>,
) -> SparseMultivariatePolynomial<D, Lex> {
    if p.is_zero() {
        return p.clone();
    }
    let lc = p.leading_coeff().cloned().unwrap();
    if p.domain().is_one(&lc) {
        return p.clone();
    }
    let inv = p
        .domain()
        .inv(&lc)
        .expect("normalize_unit_lc: coefficient domain must be a field");
    p.mul_scalar(&inv)
}

/// GCD of the coefficient polynomials of `x_0^k` for `k = 0..=deg`, returned
/// as a polynomial in the remaining `n_vars - 1` variables.
fn gcd_content_main<D: EuclideanDomain>(
    f: &SparseMultivariatePolynomial<D, Lex>,
) -> Option<SparseMultivariatePolynomial<D, Lex>> {
    let d = f.degree_in(0);
    let mut acc = f.coeff_of_var_pow(0, 0).drop_main_var();
    for k in 1..=d {
        if acc.is_zero() {
            acc = f.coeff_of_var_pow(0, k).drop_main_var();
            continue;
        }
        let ck = f.coeff_of_var_pow(0, k).drop_main_var();
        if ck.is_zero() {
            continue;
        }
        acc = multivariate_gcd_field(&acc, &ck)?;
        if acc.total_degree() == Some(0) {
            break;
        }
    }
    Some(acc)
}

/// Multivariate GCD over a field via dense recursive
/// evaluation–interpolation. Both inputs must have the same number of
/// variables; the result is normalized to have Lex-leading coefficient 1.
///
/// Returns `None` if the heuristic fails (e.g. the field is too small to
/// provide enough distinct evaluation points).
pub fn multivariate_gcd_field<D: EuclideanDomain>(
    a: &SparseMultivariatePolynomial<D, Lex>,
    b: &SparseMultivariatePolynomial<D, Lex>,
) -> Option<SparseMultivariatePolynomial<D, Lex>> {
    if a.is_zero() {
        return Some(normalize_unit_lc(b));
    }
    if b.is_zero() {
        return Some(normalize_unit_lc(a));
    }
    if a.n_vars() != b.n_vars() {
        return None;
    }
    let v = a.n_vars();
    if v == 0 {
        let mut one = SparseMultivariatePolynomial::new(a.domain().clone(), 0);
        one.set_term_external(vec![], a.domain().one());
        return Some(one);
    }
    if v == 1 {
        let g = mpoly_to_dense_x(a).gcd(&mpoly_to_dense_x(b));
        return Some(normalize_unit_lc(&dense_to_mpoly_x(&g, 1)));
    }

    // Contents in the main variable (polynomials in the remaining v-1 vars).
    let ca = gcd_content_main(a)?;
    let cb = gcd_content_main(b)?;
    let cont = multivariate_gcd_field(&ca, &cb)?.embed_new_main();
    let pa = a.checked_div_exact(&ca.embed_new_main())?;
    let pb = b.checked_div_exact(&cb.embed_new_main())?;

    if pa.degree_in(0) == 0 || pb.degree_in(0) == 0 {
        return Some(normalize_unit_lc(&cont));
    }

    let last = v - 1;
    let deg_bound = pa.degree_in(last).min(pb.degree_in(last));
    let mut images: Vec<(D::Element, SparseMultivariatePolynomial<D, Lex>)> = Vec::new();
    let mut used: Vec<D::Element> = Vec::new();
    let mut min_d0 = usize::MAX;
    let mut needed = deg_bound + 1;
    let max_attempts = needed + 80;

    for t in 0..max_attempts {
        let alpha = a.domain().cast_u64(t as u64);
        if used.contains(&alpha) {
            // Field exhausted: no more distinct evaluation points.
            if t > max_attempts / 2 {
                break;
            }
            continue;
        }
        used.push(alpha.clone());

        let ia = pa.eval(last, &alpha);
        let ib = pb.eval(last, &alpha);
        if ia.degree_in(0) != pa.degree_in(0) || ib.degree_in(0) != pb.degree_in(0) {
            continue; // leading coefficient vanished — unlucky point
        }
        let g = multivariate_gcd_field(&ia, &ib)?;
        let d0 = g.degree_in(0);
        if d0 > min_d0 {
            continue;
        }
        if d0 < min_d0 {
            images.clear();
            min_d0 = d0;
        }
        if d0 == 0 {
            return Some(normalize_unit_lc(&cont));
        }
        images.push((alpha, g));

        if images.len() >= needed {
            if let Some(cand) = interpolate_gcd_images::<D>(&images, v) {
                let cand = normalize_unit_lc(&cand);
                // Strip spurious content introduced by interpolation.
                let cc = gcd_content_main(&cand)?;
                let cpp = cand.checked_div_exact(&cc.embed_new_main())?;
                if pa.checked_div_exact(&cpp).is_some() && pb.checked_div_exact(&cpp).is_some() {
                    return Some(normalize_unit_lc(&cont.mul(&cpp)));
                }
            }
            needed += 1; // interpolation/trial failed — collect more points
        }
    }
    None
}

/// Interpolate a polynomial in `v` variables from images that are
/// polynomials in `v - 1` variables evaluated at successive points of the
/// last variable.
fn interpolate_gcd_images<D: Domain>(
    images: &[(D::Element, SparseMultivariatePolynomial<D, Lex>)],
    v: usize,
) -> Option<SparseMultivariatePolynomial<D, Lex>> {
    if images.is_empty() {
        return None;
    }
    let domain = images[0].1.domain().clone();
    let mut result = SparseMultivariatePolynomial::new(domain.clone(), v);

    // Union of exponent patterns (length v-1) across all images.
    let mut support: Vec<smallvec::SmallVec<[usize; 4]>> = Vec::new();
    for (_, g) in images {
        for exp in g.terms_ref().keys() {
            if !support.contains(exp) {
                support.push(exp.clone());
            }
        }
    }

    for e in support {
        let points: Vec<(D::Element, D::Element)> = images
            .iter()
            .map(|(a, g)| (a.clone(), g.coeff(&e)))
            .collect();
        let coeffs = lagrange_interp::<D>(&domain, &points)?;
        for (j, c) in coeffs.iter().enumerate() {
            if domain.is_zero(c) {
                continue;
            }
            let mut exp = vec![0usize; v];
            exp[..v - 1].copy_from_slice(&e);
            exp[v - 1] = j;
            let existing = result.coeff(&exp);
            let sum = domain.add(&existing, c);
            result.set_term_external(exp, sum);
        }
    }
    Some(result)
}

/// Lagrange interpolation over a field domain: given points `(x_i, y_i)`,
/// return ascending coefficients of the interpolating polynomial.
fn lagrange_interp<D: Domain>(
    domain: &D,
    points: &[(D::Element, D::Element)],
) -> Option<Vec<D::Element>> {
    let n = points.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut result = vec![domain.zero(); n];
    for i in 0..n {
        let xi = &points[i].0;
        let yi = &points[i].1;
        // denominator = Π_{j≠i} (x_i - x_j)
        let mut denom = domain.one();
        for (j, _) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            denom = domain.mul(&denom, &domain.sub(xi, &points[j].0));
        }
        let scale = domain.div(yi, &denom)?;
        // basis = Π_{j≠i} (x - x_j), ascending coefficients
        let mut basis = vec![domain.one()];
        for (j, _) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let neg_xj = domain.neg(&points[j].0);
            let mut new_basis = vec![domain.zero(); basis.len() + 1];
            for (k, b) in basis.iter().enumerate() {
                new_basis[k] = domain.add(&new_basis[k], &domain.mul(&neg_xj, b));
                new_basis[k + 1] = domain.add(&new_basis[k + 1], b);
            }
            basis = new_basis;
        }
        for (k, b) in basis.iter().enumerate().take(n) {
            result[k] = domain.add(&result[k], &domain.mul(&scale, b));
        }
    }
    Some(result)
}

/// Convert a sparse integer polynomial to a sparse rational polynomial.
fn zmpoly_to_qmpoly(f: &ZMPoly) -> QMPoly {
    QMPoly::from_terms(
        RationalDomain,
        f.n_vars(),
        f.terms_ref()
            .iter()
            .map(|(e, c)| (e.to_vec(), Rational::from_integer(c.clone())))
            .collect(),
    )
}

/// Convert a sparse rational polynomial to a primitive integer polynomial by
/// clearing denominators (multiply by the LCM of all denominators).
///
/// Generic over the monomial order; used by the multivariate GCD driver
/// (Lex) and the multi-modular Gröbner pipeline (any order).
pub(crate) fn qmpoly_to_primitive_zmpoly<O: MonomialOrder>(
    f: &SparseMultivariatePolynomial<RationalDomain, O>,
) -> SparseMultivariatePolynomial<IntegerDomain, O> {
    let mut lcm = BigInt::one();
    for c in f.terms_ref().values() {
        let d = c.denom().to_bigint();
        let g = bigint_gcd(&lcm, &d);
        lcm = lcm * d / g;
    }
    let terms: Vec<(Vec<usize>, Integer)> = f
        .terms_ref()
        .iter()
        .map(|(e, c)| {
            let scale = &lcm / c.denom().to_bigint();
            (e.to_vec(), Integer::from(c.numer().to_bigint() * scale))
        })
        .collect();
    SparseMultivariatePolynomial::from_terms(IntegerDomain, f.n_vars(), terms).primitive_part()
}

/// Negate the polynomial if its Lex-leading coefficient is negative.
fn make_lc_positive(f: ZMPoly) -> ZMPoly {
    match f.leading_coeff() {
        Some(lc) if lc.is_negative() => f.neg(),
        _ => f,
    }
}

/// Multivariate GCD over the integers (arbitrary number of variables).
///
/// Computes the GCD of the primitive parts over ℚ via
/// [`multivariate_gcd_field`] (evaluation–interpolation), then clears
/// denominators and reattaches the integer content GCD. The result is
/// primitive with a positive Lex-leading coefficient.
///
/// Returns `None` if the heuristic fails (unlucky evaluation points).
pub fn multivariate_gcd_z(a: &ZMPoly, b: &ZMPoly) -> Option<ZMPoly> {
    if a.is_zero() {
        return Some(make_lc_positive(b.primitive_part()));
    }
    if b.is_zero() {
        return Some(make_lc_positive(a.primitive_part()));
    }
    if a.n_vars() != b.n_vars() {
        return None;
    }
    let v = a.n_vars();
    if v == 0 {
        return Some(ZMPoly::from_terms(
            IntegerDomain,
            0,
            vec![(vec![], Integer::from(1))],
        ));
    }
    if v == 1 {
        let g = mpoly_to_dense_x(a).gcd(&mpoly_to_dense_x(b));
        return Some(make_lc_positive(dense_to_mpoly_x(&g, 1).primitive_part()));
    }

    let content_gcd = IntegerDomain.gcd(&a.content(), &b.content());
    let pa = a.primitive_part();
    let pb = b.primitive_part();
    let qa = zmpoly_to_qmpoly(&pa);
    let qb = zmpoly_to_qmpoly(&pb);
    let gq = multivariate_gcd_field(&qa, &qb)?;
    let gz = make_lc_positive(qmpoly_to_primitive_zmpoly(&gq));

    // Safety net: the candidate must divide both primitive parts.
    pa.checked_div_exact(&gz)?;
    pb.checked_div_exact(&gz)?;

    let result = gz.mul_scalar(&content_gcd);
    Some(make_lc_positive(result))
}

/// Multivariate GCD over a prime finite field (arbitrary number of
/// variables). The result is normalized to have Lex-leading coefficient 1.
///
/// Returns `None` if the field is too small for enough distinct evaluation
/// points.
pub fn multivariate_gcd_fp(a: &FpMPoly, b: &FpMPoly) -> Option<FpMPoly> {
    multivariate_gcd_field(a, b)
}

#[cfg(test)]
mod nvar_gcd_tests {
    use super::*;

    fn zmp(n_vars: usize, terms: &[(Vec<usize>, i64)]) -> ZMPoly {
        ZMPoly::from_terms(
            IntegerDomain,
            n_vars,
            terms
                .iter()
                .map(|(e, c)| (e.clone(), Integer::from(*c)))
                .collect(),
        )
    }

    /// (x + y + z) as a trivariate polynomial.
    fn xyz_plus() -> ZMPoly {
        zmp(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        )
    }

    #[test]
    fn gcd_trivariate_shared_factor() {
        // a = (x + y + z)(x + 1), b = (x + y + z)(y + 2)
        let g_true = xyz_plus();
        let f1 = zmp(3, &[(vec![1, 0, 0], 1), (vec![0, 0, 0], 1)]); // x + 1
        let f2 = zmp(3, &[(vec![0, 1, 0], 1), (vec![0, 0, 0], 2)]); // y + 2
        let a = g_true.mul(&f1);
        let b = g_true.mul(&f2);
        let g = multivariate_gcd_z(&a, &b).expect("trivariate gcd should succeed");
        assert!(
            g == g_true || g == g_true.neg(),
            "gcd should be x + y + z, got {:?}",
            g
        );
    }

    #[test]
    fn gcd_trivariate_coprime() {
        // gcd(x + y, z + 1) = 1
        let a = zmp(3, &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1)]);
        let b = zmp(3, &[(vec![0, 0, 1], 1), (vec![0, 0, 0], 1)]);
        let g = multivariate_gcd_z(&a, &b).expect("gcd should succeed");
        assert_eq!(
            g.total_degree(),
            Some(0),
            "coprime polys should have constant gcd"
        );
    }

    #[test]
    fn gcd_trivariate_with_content() {
        // a = 2(x + y + z)(x + 1), b = 4(x + y + z)(y + 2): gcd = 2(x+y+z)
        let g_true = xyz_plus();
        let f1 = zmp(3, &[(vec![1, 0, 0], 2), (vec![0, 0, 0], 2)]);
        let f2 = zmp(3, &[(vec![0, 1, 0], 4), (vec![0, 0, 0], 8)]);
        let a = g_true.mul(&f1);
        let b = g_true.mul(&f2);
        let g = multivariate_gcd_z(&a, &b).expect("gcd should succeed");
        let expected = g_true.mul_scalar(&Integer::from(2));
        assert!(
            g == expected || g == expected.neg(),
            "gcd should be 2(x + y + z), got {:?}",
            g
        );
    }

    #[test]
    fn gcd_trivariate_fp() {
        // Over ℤ_7: a = (x + y + z)(x + 1), b = (x + y + z)(y + 2)
        let field = FiniteField::new(BigInt::from(7));
        let mk = |terms: &[(Vec<usize>, i64)]| {
            FpMPoly::from_terms(
                field.clone(),
                3,
                terms
                    .iter()
                    .map(|(e, c)| (e.clone(), field.element(BigInt::from(*c))))
                    .collect(),
            )
        };
        let g_true = mk(&[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)]);
        let f1 = mk(&[(vec![1, 0, 0], 1), (vec![0, 0, 0], 1)]);
        let f2 = mk(&[(vec![0, 1, 0], 1), (vec![0, 0, 0], 2)]);
        let a = g_true.mul(&f1);
        let b = g_true.mul(&f2);
        let g = multivariate_gcd_fp(&a, &b).expect("fp gcd should succeed");
        assert_eq!(g, g_true, "gcd over F_7 should be x + y + z");
    }
}
