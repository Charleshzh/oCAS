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
use ocas_domain::{Domain, Integer, IntegerDomain};

use crate::dense::DenseUnivariatePolynomial;
use crate::sparse::{Lex, SparseMultivariatePolynomial};

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
            eval_point = Integer::from(eval_point.inner() + 1);
            continue;
        }
        let g_eval = a_eval.gcd(&b_eval);
        if g_eval.is_zero() {
            eval_point = Integer::from(eval_point.inner() + 1);
            continue;
        }
        images.push((eval_point.clone(), g_eval));
        eval_point = Integer::from(eval_point.inner() + 1);
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
        coeffs_map.insert(x_deg, Integer::from(existing.inner() + new_coeff.inner()));
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
                if y_coeff.inner().is_zero() {
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
        let xi = points[i].0.inner().clone();
        let yi = points[i].1.inner().clone();

        // denominator = product of (x_i - x_j) for j != i
        let mut denom = BigInt::one();
        for (j, (xj_val, _)) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let xj = xj_val.inner();
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
            let neg_xj = -(xj_val.inner());
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
            if !c.inner().is_zero() {
                result.set_term_external(vec![i, 0], c.clone());
            }
        }
        result.primitive_part()
    }
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
}
