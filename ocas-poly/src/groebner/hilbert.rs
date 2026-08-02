//! Hilbert-series bounds for monomial ideals.
//!
//! For a monomial ideal `⟨m₁, …, mₛ⟩` the Hilbert series of `R/I` is
//! `H(t) = (Σₖ (-1)ᵏ Σ_{|S|=k} t^{deg lcm(S)}) / (1-t)ⁿ`. The degree of
//! the numerator (the *regularity* of the staircase) bounds the highest
//! degree F4 must reach before all remaining S-polynomials reduce to
//! zero — a sound early-termination hint (Bayer–Stillman).
//!
//! This module computes the staircase Hilbert function incrementally.
//! It is used experimentally by F4; the bound is advisory and never
//! changes the computed basis.

use crate::sparse::monomial_lcm;

/// The Hilbert numerator of a monomial ideal: coefficients of
/// `Σₖ (-1)ᵏ Σ_{|S|=k} t^{deg lcm(S)}` as a sparse map (degree → coeff).
///
/// Computed by the inclusion-exclusion principle over the generators
/// (practical for up to ~20 generators).
pub fn hilbert_numerator(generators: &[Vec<usize>]) -> Vec<(usize, i64)> {
    use std::collections::BTreeMap;
    let mut coeffs: BTreeMap<usize, i64> = BTreeMap::new();
    let s = generators.len();
    // Inclusion-exclusion over non-empty subsets.
    for mask in 1..(1u64 << s) {
        let mut lcm: Option<Vec<usize>> = None;
        let mut bits = 0;
        for (i, g) in generators.iter().enumerate() {
            if mask & (1 << i) != 0 {
                bits += 1;
                lcm = Some(match lcm {
                    None => g.clone(),
                    Some(prev) => monomial_lcm(&prev, g).to_vec(),
                });
            }
        }
        let deg: usize = lcm.map(|l| l.iter().sum()).unwrap_or(0);
        let sign: i64 = if bits % 2 == 1 { 1 } else { -1 };
        *coeffs.entry(deg).or_insert(0) += sign;
    }
    coeffs.into_iter().filter(|&(_, c)| c != 0).collect()
}

/// The regularity bound of the staircase: the highest degree `d` for
/// which the Hilbert numerator has a non-zero coefficient. F4 may stop
/// selecting pairs above this degree when the ideal is zero-dimensional.
pub fn regularity_bound(generators: &[Vec<usize>]) -> usize {
    hilbert_numerator(generators)
        .iter()
        .map(|&(d, _)| d)
        .max()
        .unwrap_or(0)
}

/// The dimension of the staircase (vector-space dimension of `R/I` for
/// zero-dimensional ideals), from the Hilbert numerator evaluated at 1.
/// Returns `None` when the ideal is positive-dimensional (numerator sums
/// to 0).
pub fn staircase_dimension(generators: &[Vec<usize>]) -> Option<usize> {
    let sum: i64 = hilbert_numerator(generators).iter().map(|&(_, c)| c).sum();
    if sum == 0 {
        None
    } else {
        Some(sum.unsigned_abs() as usize)
    }
}

// ------------------------------------------------------------------
//  Complete Hilbert series for arbitrary ideals
// ------------------------------------------------------------------

/// The Hilbert series of a quotient ring $R/I$, represented as a rational
/// function $H(t) = N(t) / (1-t)^n$.
///
/// The numerator $N(t)$ is stored as a vector of coefficients:
/// `numerator[i]` is the coefficient of $t^i$.
/// The denominator is $(1-t)^n$ where $n$ is the number of variables.
#[derive(Debug, Clone)]
pub struct HilbertSeries {
    /// Numerator coefficients (from constant term upward).
    pub numerator: Vec<i64>,
    /// The power of $(1-t)$ in the denominator (= number of variables).
    pub denominator_power: usize,
}

impl HilbertSeries {
    /// Evaluate the Hilbert function at degree $d$: $\dim_k (R/I)_d$.
    ///
    /// Uses the formula $H(d) = [t^d] N(t) / (1-t)^n$.
    /// The coefficient of $t^k$ in $(1-t)^{-n}$ is $\binom{n+k-1}{k}$.
    pub fn hilbert_function(&self, degree: usize) -> i64 {
        let n = self.denominator_power as i64;
        let mut result = 0i64;
        for (i, &coeff) in self.numerator.iter().enumerate() {
            if i > degree {
                break;
            }
            let k = degree - i;
            let binom = binomial_general(n, k);
            result += coeff * binom;
        }
        result
    }

    /// Compute the Krull dimension of $R/I$.
    ///
    /// This is the degree of the Hilbert polynomial, computed by checking
    /// how many factors of $(1-t)$ divide the numerator.
    pub fn dimension(&self) -> usize {
        // If N(1) != 0, dimension = n.
        // Otherwise, we need to find the order of vanishing at t=1.
        let n = self.denominator_power;
        // Evaluate derivatives at t=1 to find multiplicity of root.
        let mut poly = self.numerator.clone();
        for dim in 0..n {
            let sum: i64 = poly.iter().sum();
            if sum != 0 {
                return n - dim;
            }
            // Differentiate: if poly = Σ aᵢ tⁱ, then poly' = Σ i·aᵢ t^{i-1}
            // Equivalently, shift and scale.
            let mut new_poly = Vec::with_capacity(poly.len().saturating_sub(1));
            for (i, &c) in poly.iter().enumerate().skip(1) {
                new_poly.push(c * i as i64);
            }
            poly = new_poly;
        }
        0
    }

    /// The degree of the projective variety.
    ///
    /// For a well-formed Hilbert series, this equals the value of the
    /// numerator at $t=1$ after dividing out the dimensional factors.
    pub fn degree(&self) -> i64 {
        let mut poly = self.numerator.clone();
        let dim = self.dimension();
        // Differentiate `n - dim` times and evaluate at t=1.
        for _ in 0..(self.denominator_power - dim) {
            let mut new_poly = Vec::with_capacity(poly.len().saturating_sub(1));
            for (i, &c) in poly.iter().enumerate().skip(1) {
                new_poly.push(c * i as i64);
            }
            poly = new_poly;
        }
        let sum: i64 = poly.iter().sum();
        // Divide by (n-dim)! to get the degree.
        let factorial: i64 = (1..=(self.denominator_power - dim) as i64).product();
        if factorial == 0 { sum } else { sum / factorial }
    }

    /// Compute the coefficients of the Hilbert polynomial $P(d)$ such that
    /// $H(d) = P(d)$ for $d \gg 0$, where $H(d) = \dim_k (R/I)_d$.
    ///
    /// Returns coefficients in ascending-degree order: `result[i]` is the
    /// coefficient of $d^i$. The polynomial has degree `self.dimension()`.
    ///
    /// Uses Lagrange interpolation on $d = 0, 1, \ldots, \dim$ since the
    /// Hilbert function agrees with the polynomial for $d \geq 0$ when
    /// the ideal is homogeneous.
    pub fn hilbert_polynomial(&self) -> Vec<f64> {
        let dim = self.dimension();
        let n = dim + 1; // number of interpolation points

        // Evaluate the Hilbert function at d = 0, 1, ..., dim.
        let xs: Vec<f64> = (0..n).map(|d| d as f64).collect();
        let ys: Vec<f64> = (0..n).map(|d| self.hilbert_function(d) as f64).collect();

        // Lagrange interpolation: compute coefficients in the monomial basis.
        // Start with the zero polynomial and add each Lagrange basis polynomial
        // scaled by its y-value.
        let mut coeffs = vec![0.0f64; n]; // coeffs[i] = coefficient of d^i

        for j in 0..n {
            // Build the Lagrange basis polynomial L_j(d) =
            // Π_{m ≠ j} (d - x_m) / (x_j - x_m).
            //
            // Start with the numerator polynomial Π_{m ≠ j} (d - x_m)
            // represented as coefficients [c_0, c_1, ...].
            let mut basis = vec![1.0f64]; // constant 1
            let mut denom = 1.0f64;
            for m in 0..n {
                if m == j {
                    continue;
                }
                // Multiply basis by (d - x_m).
                let xm = xs[m];
                let mut new_basis = vec![0.0f64; basis.len() + 1];
                for (k, &c) in basis.iter().enumerate() {
                    new_basis[k] -= c * xm;
                    new_basis[k + 1] += c;
                }
                basis = new_basis;
                denom *= xs[j] - xm;
            }
            // Scale by y_j / denom and accumulate.
            let scale = ys[j] / denom;
            for (k, &c) in basis.iter().enumerate() {
                coeffs[k] += c * scale;
            }
        }

        // Trim trailing near-zero coefficients.
        while coeffs.len() > 1 && coeffs.last().unwrap().abs() < 1e-10 {
            coeffs.pop();
        }

        coeffs
    }
}

/// Compute the Hilbert series of $R/I$ from its Gröbner basis.
///
/// Uses Macaulay's theorem: the Hilbert series of $R/I$ equals that of
/// $R/\text{LM}(I)$, the quotient by the leading-term ideal.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::Lex;
/// use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
/// use ocas_poly::groebner::hilbert::hilbert_series;
///
/// let d = RationalDomain;
/// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![2, 0], Rational::new(1, 1)),
/// ]);
/// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 1], Rational::new(1, 1)),
/// ]);
/// let gb = groebner_basis(&[f1, f2], Algorithm::F4);
/// let hs = hilbert_series(&gb);
/// assert!(hs.dimension() <= 2);
/// ```
pub fn hilbert_series(
    gb: &crate::groebner::GroebnerBasis<ocas_domain::RationalDomain, crate::sparse::Lex>,
) -> HilbertSeries {
    let n_vars = gb.basis.first().map(|p| p.n_vars()).unwrap_or(0);

    // Extract leading monomials.
    let lms: Vec<Vec<usize>> = gb
        .basis
        .iter()
        .filter_map(|p| p.leading_monomial().map(|m| m.to_vec()))
        .collect();

    if lms.is_empty() {
        return HilbertSeries {
            numerator: vec![1],
            denominator_power: n_vars,
        };
    }

    let num = hilbert_numerator(&lms);
    let max_deg = num.iter().map(|&(d, _)| d).max().unwrap_or(0);
    let mut numerator = vec![0i64; max_deg + 1];
    // Leading term: +1 at degree 0.
    numerator[0] = 1;
    // hilbert_numerator uses (-1)^{k+1} signs; negate for the correct
    // R/I Hilbert numerator which uses (-1)^k.
    for (deg, coeff) in num {
        numerator[deg] -= coeff;
    }

    HilbertSeries {
        numerator,
        denominator_power: n_vars,
    }
}

/// Compute $\binom{n+k-1}{k}$ for integer $n$ and non-negative integer $k$.
/// This is the coefficient of $t^k$ in $(1-t)^{-n}$.
fn binomial_general(n: i64, k: usize) -> i64 {
    if k == 0 {
        return 1;
    }
    let mut result = 1i64;
    for i in 0..k {
        result = result * (n + i as i64) / (i as i64 + 1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hilbert_numerator_single_generator() {
        // <x²> in 1 variable: numerator 1 - t².
        let coeffs = hilbert_numerator(&[vec![2]]);
        assert_eq!(coeffs, vec![(2, 1)]);
        // dim of staircase = 2 (monomials 1, x).
    }

    #[test]
    fn hilbert_numerator_two_generators() {
        // <x², y²> in 2 variables: numerator 1 - t² - t² + t⁴ = 1 - 2t² + t⁴.
        let coeffs = hilbert_numerator(&[vec![2, 0], vec![0, 2]]);
        assert_eq!(coeffs, vec![(2, 2), (4, -1)]);
        assert_eq!(regularity_bound(&[vec![2, 0], vec![0, 2]]), 4);
        assert_eq!(staircase_dimension(&[vec![2, 0], vec![0, 2]]), Some(1));
    }

    #[test]
    fn hilbert_numerator_linear() {
        // <x, y> in 2 variables: numerator 1 - t - t + t² = (1-t)².
        let coeffs = hilbert_numerator(&[vec![1, 0], vec![0, 1]]);
        assert_eq!(coeffs, vec![(1, 2), (2, -1)]);
        assert_eq!(regularity_bound(&[vec![1, 0], vec![0, 1]]), 2);
        assert_eq!(staircase_dimension(&[vec![1, 0], vec![0, 1]]), Some(1));
    }

    #[test]
    fn hilbert_series_xy_xz() {
        // Ideal (x^2, xy) in k[x,y].
        // LM = {x^2, xy}. Hilbert numerator: 1 - t^2.
        // (LMs: x^2 deg 2, xy deg 2. lcm(x^2,xy)=x^2y deg 3.)
        // N(t) = 1 - 2t^2 + t^3. H(t) = (1 - 2t^2 + t^3)/(1-t)^2.
        use crate::groebner::hilbert::hilbert_series;
        use crate::sparse::Lex;
        use crate::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
        use ocas_domain::{Rational, RationalDomain};

        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![2, 0], Rational::new(1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 1], Rational::new(1, 1))],
        );
        let gb = groebner_basis(&[f1, f2], Algorithm::F4);
        let hs = hilbert_series(&gb);
        // H(0) = 1 (the constant monomial).
        assert_eq!(hs.hilbert_function(0), 1);
        // H(1) = 2 (monomials x, y in degree 1).
        assert_eq!(hs.hilbert_function(1), 2);
    }

    #[test]
    fn hilbert_series_linear_ideal() {
        // Ideal (x, y) in k[x,y] → R/I has dim 1 (just constants).
        use crate::groebner::hilbert::hilbert_series;
        use crate::sparse::Lex;
        use crate::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
        use ocas_domain::{Rational, RationalDomain};

        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 0], Rational::new(1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![0, 1], Rational::new(1, 1))],
        );
        let gb = groebner_basis(&[f1, f2], Algorithm::F4);
        let hs = hilbert_series(&gb);
        // H(0) = 1, H(d) = 0 for d > 0.
        assert_eq!(hs.hilbert_function(0), 1);
        assert_eq!(hs.hilbert_function(1), 0);
    }

    #[test]
    fn hilbert_polynomial_linear_ideal() {
        // Ideal (x, y) in k[x,y]: dim=0, degree=1.
        // Hilbert polynomial = 1 (constant).
        use crate::groebner::hilbert::hilbert_series;
        use crate::sparse::Lex;
        use crate::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
        use ocas_domain::{Rational, RationalDomain};

        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 0], Rational::new(1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![0, 1], Rational::new(1, 1))],
        );
        let gb = groebner_basis(&[f1, f2], Algorithm::F4);
        let hs = hilbert_series(&gb);
        let hp = hs.hilbert_polynomial();
        // P(d) = 1 for all d.
        assert_eq!(hp.len(), 1);
        assert!((hp[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn hilbert_polynomial_empty_ideal() {
        // Empty ideal in k[x,y]: R/I = k[x,y], H(d) = d+1.
        // Hilbert polynomial = 1 + d → coeffs [1, 1].
        let hs = HilbertSeries {
            numerator: vec![1],
            denominator_power: 2,
        };
        let hp = hs.hilbert_polynomial();
        assert_eq!(hp.len(), 2);
        assert!((hp[0] - 1.0).abs() < 1e-10);
        assert!((hp[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn hilbert_polynomial_x_squared() {
        // Ideal (x²) in k[x]: H(0)=1, H(1)=1, H(d)=1 for d≥0.
        // dim=0, polynomial = 1.
        let hs = HilbertSeries {
            numerator: vec![1, 0, -1], // 1 - t²
            denominator_power: 1,
        };
        let hp = hs.hilbert_polynomial();
        assert_eq!(hp.len(), 1);
        assert!((hp[0] - 1.0).abs() < 1e-10);
    }
}
