//! Polynomial GCD (greatest common divisor) algorithms.
//!
//! Implements the Euclidean algorithm for dense univariate polynomials
//! over any [`EuclideanDomain`]. For non-field domains (e.g. Z[x]),
//! pseudo-remainders are used to avoid fractional coefficients.
//!
//! For large integer coefficients, prefer the modular Brown algorithm in
//! [`modular`], which avoids the coefficient explosion of pseudo-remainders.

use ocas_domain::EuclideanDomain;

use crate::dense::DenseUnivariatePolynomial;

pub mod modular;

impl<D: EuclideanDomain> DenseUnivariatePolynomial<D> {
    /// Compute the pseudo-remainder of `self` divided by `other`.
    ///
    /// For polynomials over a non-field ring, standard division may fail
    /// because leading coefficients do not divide. Pseudo-division
    /// multiplies the dividend by `lc(divisor)^(deg(dividend) - deg(divisor) + 1)`
    /// before dividing, guaranteeing exact coefficient division.
    ///
    /// Returns `None` if `other` is zero or if the degree of `self` is
    /// less than the degree of `other`.
    pub(crate) fn pseudo_remainder(&self, divisor: &Self) -> Option<Self> {
        let self_deg = self.degree()?;
        let div_deg = divisor.degree()?;
        if self_deg < div_deg {
            return Some(self.clone());
        }

        let d = self.domain();
        let div_lc = divisor.leading_coeff()?;
        let mut remainder = self.clone();

        let exponent = self_deg - div_deg + 1;

        // Multiply by lc(divisor)^exponent.
        let factor = d.pow(div_lc, exponent as u64);
        remainder = remainder.mul_scalar(&factor);

        // Now perform standard polynomial division.
        let mut quot_coeffs = vec![d.zero(); self_deg - div_deg + 1];

        while let Some(deg) = remainder.degree() {
            if deg < div_deg {
                break;
            }
            let lc = remainder.leading_coeff().unwrap().clone();
            let (q, _) = d.div_rem(&lc, div_lc)?;
            let term_degree = deg - div_deg;
            quot_coeffs[term_degree] = d.add(&quot_coeffs[term_degree], &q);

            let mut sub_coeffs = vec![d.zero(); term_degree];
            sub_coeffs.extend(divisor.coeffs().iter().map(|c| d.mul(c, &q)));
            let sub = Self::from_coeffs(d.clone(), sub_coeffs);
            remainder = remainder.sub(&sub);

            if let Some(rem_deg) = remainder.degree() {
                if rem_deg >= deg {
                    break;
                }
            } else {
                break;
            }
        }

        Some(remainder)
    }

    /// Compute the greatest common divisor of `self` and `other`.
    ///
    /// Uses the subresultant PRS (Brown–Traub): the same recurrence as
    /// [`Self::resultant`], with exact division by `beta` at every step.
    /// The naive pseudo-remainder sequence without content control explodes
    /// coefficient sizes at moderate degrees (≥ 16 over ℤ, and over ℚ where
    /// `primitive_part` is only a unit scaling); the subresultant scaling
    /// keeps intermediate coefficients at the theoretical subresultant
    /// bound. The result is always primitive (content-free).
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let d = IntegerDomain;
    /// let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(-1), Integer::from(0), Integer::from(1),
    /// ]); // x^2 - 1 = (x-1)(x+1)
    /// let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    ///     Integer::from(1), Integer::from(2), Integer::from(1),
    /// ]); // x^2 + 2x + 1 = (x+1)^2
    /// let g = a.gcd(&b);
    /// assert_eq!(g.coeffs(), &[Integer::from(1), Integer::from(1)]); // x + 1
    /// ```
    pub fn gcd(&self, other: &Self) -> Self {
        if other.is_zero() {
            return self.primitive_part();
        }
        if self.is_zero() {
            return other.primitive_part();
        }

        let d = self.domain();
        let mut a = self.primitive_part();
        let mut a_new = other.primitive_part();
        if a.degree() < a_new.degree() {
            std::mem::swap(&mut a, &mut a_new);
        }
        // A constant divisor divides every polynomial up to content.
        if a_new.degree() == Some(0) {
            return self.one();
        }

        let mut deg = (a.degree().expect("nonzero") - a_new.degree().expect("nonzero")) as u64;
        let mut neg_lc = d.one(); // set before use
        let mut init = false;
        let mut beta = d.pow(&d.neg(&d.one()), deg + 1);
        let mut psi = d.neg(&d.one());

        loop {
            if init {
                // Update psi and beta (same recurrence as `resultant`).
                psi = if deg == 0 {
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

            if r_reduced.is_zero() {
                return a_new.primitive_part();
            }
            a = a_new;
            a_new = r_reduced;
            if a_new.degree() == Some(0) {
                // Nonzero constant remainder: coprime.
                return self.one();
            }
        }
    }

    /// Compute the content of this polynomial: the GCD of all its coefficients.
    ///
    /// For the zero polynomial the content is zero.
    pub fn content(&self) -> D::Element {
        if self.is_zero() {
            return self.domain().zero();
        }
        let coeffs = self.coeffs();
        let mut g = coeffs[0].clone();
        for c in &coeffs[1..] {
            g = self.domain().gcd(&g, c);
            if self.domain().is_one(&g) {
                break;
            }
        }
        g
    }

    /// Return the primitive part of this polynomial (polynomial / content).
    pub fn primitive_part(&self) -> Self {
        if self.is_zero() {
            return self.zero();
        }
        let content = self.content();
        let coeffs: Vec<D::Element> = self
            .coeffs()
            .iter()
            .map(|c| self.domain().div(c, &content).unwrap_or_else(|| c.clone()))
            .collect();
        Self::from_coeffs(self.domain().clone(), coeffs)
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
    fn gcd_x2_minus_1_and_x_plus_1() {
        let d = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
        let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(1)]);
        let g = a.gcd(&b);
        assert_eq!(g.coeffs(), &[i(1), i(1)]);
    }

    #[test]
    fn gcd_x2_minus_1_and_x2_plus_2x_plus_1() {
        let d = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
        let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);
        let g = a.gcd(&b);
        assert_eq!(g.coeffs(), &[i(1), i(1)]);
    }

    #[test]
    fn gcd_coprime() {
        let d = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(1)]);
        let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(2), i(1)]);
        let g = a.gcd(&b);
        assert_eq!(g.degree(), Some(0));
        assert!(!g.is_zero());
    }

    #[test]
    fn gcd_with_zero() {
        let d = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(2), i(4), i(2)]);
        let g = a.gcd(&a.zero());
        assert_eq!(g.coeffs(), &[i(1), i(2), i(1)]);
    }

    #[test]
    fn primitive_part_of_scaled_polynomial() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(2), i(4), i(6)]);
        let prim = p.primitive_part();
        assert_eq!(prim.coeffs(), &[i(1), i(2), i(3)]);
    }
}
