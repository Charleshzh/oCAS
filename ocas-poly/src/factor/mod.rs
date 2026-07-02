//! Polynomial factorization algorithms.
//!
//! This module groups square-free factorization with complete factorization
//! over finite fields ([`finite_field`]) and, for lifting back to the integers,
//! Hensel lifting ([`hensel`]).
//!
//! The top-level entry point for factoring a univariate polynomial over
//! $\mathbb{Z}$ is [`DenseUnivariatePolynomial::factor`](crate::DenseUnivariatePolynomial::factor),
//! and over a finite field
//! [`factor_over_finite_field`](finite_field::factor_over_finite_field).

use ocas_domain::EuclideanDomain;
use ocas_domain::{FiniteField, IntegerDomain};

use crate::dense::DenseUnivariatePolynomial;

pub mod finite_field;
pub mod hensel;

/// Result of a square-free factorization: list of (factor, multiplicity) pairs.
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

/// Result of a complete factorization: list of (factor, multiplicity) pairs
/// where each factor is irreducible (or, over the integers, primitive and
/// irreducible over $\mathbb{Q}$).
pub type Factors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

impl<D: EuclideanDomain> DenseUnivariatePolynomial<D> {
    /// Compute the square-free factorization of this polynomial.
    ///
    /// Returns a list of (factor, multiplicity) pairs.
    /// For example, `(x+1)^2 * (x-1)` yields `[(x+1, 2), (x-1, 1)]`.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let d = IntegerDomain;
    /// // (x+1)^2*(x-1) = x^3 + x^2 - x - 1
    /// let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(-1), Integer::from(-1), Integer::from(1), Integer::from(1),
    /// ]);
    /// let factors = p.square_free_factorization();
    /// assert_eq!(factors.len(), 2);
    /// ```
    pub fn square_free_factorization(&self) -> SquareFreeFactors<D> {
        let mut factors = SquareFreeFactors::new();
        if self.is_zero() {
            return factors;
        }

        // Step 1: make polynomial primitive.
        let f = self.primitive_part();
        let f_deriv = f.derivative();

        // g = gcd(f, f')
        let mut g = f.gcd(&f_deriv);
        if g.is_zero() {
            return factors;
        }

        // w = f / g contains each distinct irreducible factor exactly once.
        let mut w = match f.div_rem(&g) {
            Some((q, _)) => q,
            None => return factors,
        };

        let mut k = 1;
        while !w.is_one() {
            // h = gcd(w, g)
            let h = w.gcd(&g);
            // z = w / h is the factor with multiplicity k.
            if let Some((z, _)) = w.div_rem(&h)
                && !z.is_one()
                && !z.is_zero()
            {
                factors.push((z, k));
            }

            // Prepare for next iteration.
            if let Some((q, _)) = g.div_rem(&h) {
                g = q;
            } else {
                break;
            }
            w = h;
            k += 1;
        }

        factors
    }

    /// Check whether this polynomial is square-free.
    ///
    /// A polynomial is square-free if gcd(p, p') = 1.
    pub fn is_square_free(&self) -> bool {
        if self.degree().unwrap_or(0) <= 1 {
            return true;
        }
        let deriv = self.derivative();
        let g = self.gcd(&deriv);
        g.degree() == Some(0)
    }
}

// ── factor() for integer polynomials ──────────────────────────────

impl DenseUnivariatePolynomial<IntegerDomain> {
    /// Completely factor this primitive integer polynomial into monic
    /// irreducible factors with multiplicities.
    ///
    /// The input must be primitive (coefficient content = 1). Use
    /// [`primitive_part`](crate::DenseUnivariatePolynomial::primitive_part)
    /// to prepare an arbitrary integer polynomial before factoring.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let d = IntegerDomain;
    /// // x^2 - 1 = (x-1)(x+1)
    /// let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(-1), Integer::from(0), Integer::from(1),
    /// ]);
    /// let factors = p.factor();
    /// assert_eq!(factors.len(), 2);
    /// ```
    pub fn factor(&self) -> Factors<IntegerDomain> {
        hensel::factor_primitive(self)
    }
}

// ── factor() for finite-field polynomials ─────────────────────────

impl DenseUnivariatePolynomial<FiniteField> {
    /// Completely factor this univariate polynomial over $\mathbb{F}_p$
    /// into monic irreducible factors with multiplicities.
    ///
    /// # Example
    ///
    /// ```
    /// use num_bigint::BigInt;
    /// use ocas_domain::{Domain, FiniteField};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let f = FiniteField::new(BigInt::from(5));
    /// // x^2 - 1 over F_5
    /// let p = DenseUnivariatePolynomial::from_coeffs(
    ///     f.clone(), vec![f.element(4), f.element(0), f.element(1)]);
    /// let factors = p.factor();
    /// assert!(!factors.is_empty());
    /// ```
    pub fn factor(&self) -> Factors<FiniteField> {
        finite_field::factor_over_finite_field(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain};

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn square_free_x_plus_1_cubed() {
        let d = IntegerDomain;
        // (x+1)^3 = x^3 + 3x^2 + 3x + 1
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(3), i(3), i(1)]);
        let factors = p.square_free_factorization();
        assert!(!factors.is_empty());
        for (factor, mult) in &factors {
            if factor.degree() == Some(1) {
                assert_eq!(*mult, 3);
            }
        }
    }

    #[test]
    fn is_square_free_linear() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(1)]); // x+1
        assert!(p.is_square_free());
    }

    #[test]
    fn is_not_square_free_perfect_square() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]); // (x+1)^2
        assert!(!p.is_square_free());
    }

    #[test]
    fn square_free_x2_minus_1() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]); // x^2-1
        let factors = p.square_free_factorization();
        assert!(p.is_square_free());
        assert!(!factors.is_empty());
    }
}
