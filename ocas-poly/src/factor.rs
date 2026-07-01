//! Polynomial factorization algorithms.
//!
//! Provides square-free factorization and a skeleton for full factorization
//! over Z[x] and finite fields. Full factorization (Hensel lifting) is
//! planned for a future release.

use ocas_domain::EuclideanDomain;

use crate::dense::DenseUnivariatePolynomial;

/// Result of a square-free factorization: list of (factor, multiplicity) pairs.
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

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
