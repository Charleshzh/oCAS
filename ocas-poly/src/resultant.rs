//! Polynomial resultant computation.
//!
//! Implements Brown's Polynomial Remainder Sequence (PRS) algorithm for
//! computing the resultant of two univariate polynomials over any
//! [`EuclideanDomain`].
//!
//! The resultant of two polynomials $a$ and $b$ is zero if and only if
//! they share a common root (or equivalently, a non-trivial GCD).

use ocas_domain::EuclideanDomain;

use crate::dense::DenseUnivariatePolynomial;

impl<D: EuclideanDomain> DenseUnivariatePolynomial<D> {
    /// Compute the resultant of `self` and `other` using Brown's PRS algorithm.
    ///
    /// The resultant $\operatorname{Res}(a, b)$ is a scalar in the coefficient
    /// domain. It is zero if and only if $\gcd(a, b)$ is non-constant.
    ///
    /// Ported from Symbolica's `resultant_prs` (`src/poly/resultant.rs`):
    /// subresultant PRS with exact division by `beta` at every step (the
    /// division is exact in any UFD by the subresultant theorem).
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let d = IntegerDomain;
    /// // Res(x - 1, x - 2) = 1 - 2 = -1
    /// let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(-1), Integer::from(1),
    /// ]);
    /// let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(-2), Integer::from(1),
    /// ]);
    /// assert_eq!(a.resultant(&b), Integer::from(-1));
    /// ```
    pub fn resultant(&self, other: &Self) -> D::Element {
        let d = self.domain();

        // Ensure deg(self) >= deg(other); swap with the sign
        // Res(a, b) = (-1)^(deg a · deg b) Res(b, a).
        match (self.degree(), other.degree()) {
            (None, _) | (_, None) => return d.zero(),
            (Some(ds), Some(do_)) if ds < do_ => {
                let r = other.resultant(self);
                if ds % 2 == 1 && do_ % 2 == 1 {
                    return d.neg(&r);
                }
                return r;
            }
            _ => {}
        }

        let deg_a = self.degree().expect("nonzero polynomial");
        let deg_b = other.degree().expect("nonzero polynomial");

        // If the smaller polynomial is constant, the resultant is
        // `constant^(deg of the larger)`.
        if deg_b == 0 {
            return d.pow(&other.constant(), deg_a as u64);
        }

        let mut a = self.clone();
        let mut a_new = other.clone();

        let mut deg = (a.degree().expect("nonzero") - a_new.degree().expect("nonzero")) as u64;
        let mut neg_lc = d.one(); // set before use
        let mut init = false;
        let mut beta = d.pow(&d.neg(&d.one()), deg + 1);
        let mut psi = d.neg(&d.one());

        // Collect (leading_coeff, degree) at each step.
        let mut lcs: Vec<(D::Element, u64)> =
            vec![(a.lcoeff(), a.degree().expect("nonzero") as u64)];

        while a_new.degree().unwrap_or(0) > 0 {
            if init {
                // Update psi and beta.
                psi = if deg == 0 {
                    // Can only happen on the first iteration.
                    psi
                } else if deg == 1 {
                    neg_lc.clone()
                } else {
                    let num = d.pow(&neg_lc, deg);
                    let den = d.pow(&psi, deg - 1);
                    let (q, r) = d
                        .div_rem(&num, &den)
                        .expect("subresultant psi division is exact");
                    debug_assert!(d.is_zero(&r));
                    q
                };
                deg = (a.degree().expect("nonzero") - a_new.degree().expect("nonzero")) as u64;
                beta = d.mul(&neg_lc, &d.pow(&psi, deg));
            } else {
                init = true;
            }

            neg_lc = d.neg(a_new.leading_coeff().expect("nonzero"));

            // Pseudo-remainder: a · (−lc(b))^(deg+1) mod b, with sign.
            let factor = d.pow(&neg_lc, deg + 1);
            let (_, mut r) = a
                .mul_scalar(&factor)
                .div_rem(&a_new)
                .expect("pseudo-division succeeds after scaling");
            if (deg + 1) % 2 == 1 {
                r = r.neg();
            }

            lcs.push((a_new.lcoeff(), a_new.degree().expect("nonzero") as u64));

            // Exact scalar division by beta (subresultant theorem).
            let r_reduced = Self::from_coeffs(
                d.clone(),
                r.coeffs()
                    .iter()
                    .map(|c| {
                        d.div(c, &beta)
                            .expect("subresultant beta division is exact")
                    })
                    .collect(),
            );
            a = a_new;
            a_new = r_reduced;
        }

        // A zero remainder before reaching a constant means a common factor.
        if a_new.is_zero() {
            return d.zero();
        }
        lcs.push((a_new.lcoeff(), 0));

        // Compute the resultant from the PRS using the fundamental theorem.
        let mut rho = d.one();
        let mut den = d.one();

        for k in 1..lcs.len() {
            let mut exponent: i64 = lcs[k - 1].1 as i64 - lcs[k].1 as i64;
            // Multiply by (deg differences from remaining steps).
            for l in k..lcs.len() - 1 {
                let dl = lcs[l].1 as i64;
                let dl1 = lcs[l + 1].1 as i64;
                exponent *= 1 - (dl - dl1);
            }

            if exponent > 0 {
                let pow_val = d.pow(&lcs[k].0, exponent as u64);
                rho = d.mul(&rho, &pow_val);
            } else if exponent < 0 {
                let pow_val = d.pow(&lcs[k].0, (-exponent) as u64);
                den = d.mul(&den, &pow_val);
            }
        }

        d.div_rem(&rho, &den)
            .expect("resultant reconstruction is exact")
            .0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain};

    fn int(i: i64) -> Integer {
        Integer::from(i)
    }

    fn poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<IntegerDomain> {
        DenseUnivariatePolynomial::from_coeffs(
            IntegerDomain,
            coeffs.iter().map(|&c| int(c)).collect(),
        )
    }

    #[test]
    fn resultant_linear_different_roots() {
        // Res(x - 1, x - 2) = 1 - 2 = -1 (product of (α_i - β_j))
        let a = poly(&[-1, 1]); // x - 1
        let b = poly(&[-2, 1]); // x - 2
        assert_eq!(a.resultant(&b), int(-1));
    }

    #[test]
    fn resultant_common_root() {
        // Res(x^2 - 1, x - 1) = 0 (share root x=1)
        let a = poly(&[-1, 0, 1]); // x^2 - 1
        let b = poly(&[-1, 1]); // x - 1
        assert_eq!(a.resultant(&b), int(0));
    }

    #[test]
    fn resultant_no_common_root() {
        // Res(x^2 + 1, (x+1)^2) = 4
        let a = poly(&[1, 0, 1]); // x^2 + 1
        let b = poly(&[1, 2, 1]); // x^2 + 2x + 1
        assert_eq!(a.resultant(&b), int(4));
    }

    #[test]
    fn resultant_shared_factor() {
        // Res((x-1)(x-2), (x-1)(x-3)) = 0
        let a = poly(&[2, -3, 1]); // x^2 - 3x + 2
        let b = poly(&[3, -4, 1]); // x^2 - 4x + 3
        assert_eq!(a.resultant(&b), int(0));
    }

    #[test]
    fn resultant_constant_poly() {
        // Res(x^2 + 1, 3) = 3^2 = 9
        let a = poly(&[1, 0, 1]);
        let b = poly(&[3]);
        assert_eq!(a.resultant(&b), int(9));
    }

    #[test]
    fn resultant_constant_constant() {
        // Res(2, 3): deg_a=0, deg_b=0, b^deg_a = 3^0 = 1
        let a = poly(&[2]);
        let b = poly(&[3]);
        assert_eq!(a.resultant(&b), int(1));
    }

    #[test]
    fn resultant_symmetric_up_to_sign() {
        // Res(a, b) = (-1)^(deg_a * deg_b) * Res(b, a)
        let a = poly(&[-1, 0, 1]); // x^2 - 1, deg=2
        let b = poly(&[-2, 1]); // x - 2, deg=1
        // deg_a * deg_b = 2, so Res(a,b) = Res(b,a)
        let r1 = a.resultant(&b);
        let r2 = b.resultant(&a);
        assert_eq!(r1, r2);
    }

    #[test]
    fn resultant_zero_poly() {
        let a = poly(&[0]); // zero polynomial
        let b = poly(&[1, 1]); // x + 1
        assert_eq!(a.resultant(&b), int(0));
    }

    #[test]
    fn resultant_quartic_cubic() {
        // SymPy: resultant(x^4 - 3, 3x^3 - x^2 + 2x + 1, x) == -2243.
        // Regression: the previous implementation skipped the beta division
        // unless beta was a unit, which is not a valid resultant algorithm
        // beyond trivial degrees.
        let a = poly(&[-3, 0, 0, 0, 1]);
        let b = poly(&[1, 2, -1, 3]);
        assert_eq!(a.resultant(&b), int(-2243));
        // And swapped (both degrees odd? 4 and 3 — no sign flip).
        assert_eq!(b.resultant(&a), int(-2243));
    }
}
