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

        // Ensure deg(a) >= deg(b).
        let (a, b, swapped) = match (self.degree(), other.degree()) {
            (None, _) | (_, None) => return d.zero(),
            (Some(da), Some(db)) => {
                if da >= db {
                    (self.clone(), other.clone(), false)
                } else {
                    (other.clone(), self.clone(), true)
                }
            }
        };

        let deg_a = a.degree().unwrap();
        let deg_b = b.degree().unwrap();

        // If b is constant, return b^(deg a).
        if deg_b == 0 {
            let val = b.constant();
            let res = d.pow(&val, deg_a as u64);
            // Sign correction: (-1)^(deg_a * deg_b) when swapped.
            if swapped && (deg_a * deg_b) % 2 == 1 {
                return d.neg(&res);
            }
            return res;
        }

        // Run Brown's PRS.
        let mut a_cur = a;
        let mut b_cur = b;

        let deg_diff = a_cur.degree().unwrap() - b_cur.degree().unwrap();
        let neg_lc = d.neg(b_cur.leading_coeff().unwrap());
        let mut beta = d.pow(&d.neg(&d.one()), (deg_diff + 1) as u64);
        let mut psi = d.neg(&d.one());

        // Collect (leading_coeff, degree) at each step.
        let mut lcs: Vec<(D::Element, usize)> = Vec::new();
        lcs.push((a_cur.lcoeff(), a_cur.degree().unwrap()));

        let mut first = true;

        while !b_cur.is_zero() {
            let b_deg = b_cur.degree().unwrap();

            if !first {
                // Update psi and beta.
                let cur_deg_diff = a_cur.degree().unwrap() - b_deg;
                psi = if cur_deg_diff == 0 {
                    psi
                } else if cur_deg_diff == 1 {
                    neg_lc.clone()
                } else {
                    let a_part = d.pow(&neg_lc, cur_deg_diff as u64);
                    let psi_old = d.pow(&psi, (cur_deg_diff - 1) as u64);
                    d.div_rem(&a_part, &psi_old).unwrap().0
                };

                let new_deg_diff = a_cur.degree().unwrap() - b_deg;
                beta = d.mul(&neg_lc, &d.pow(&psi, new_deg_diff as u64));
            }
            first = false;

            let neg_lc_new = d.neg(b_cur.leading_coeff().unwrap());
            let deg_diff_now = a_cur.degree().unwrap() - b_deg;

            // Compute pseudo-remainder: a * (-lc)^(deg_diff+1) mod b.
            let factor = d.pow(&neg_lc_new, (deg_diff_now + 1) as u64);
            let scaled = a_cur.mul_scalar(&factor);
            let (_, mut r) = scaled.div_rem(&b_cur).unwrap();

            // Sign correction: (-1)^(deg_diff + 1).
            if (deg_diff_now + 1) % 2 == 1 {
                r = r.neg();
            }

            // Normalize by beta.
            if !d.is_zero(&beta) {
                let inv_beta = d.div_rem(&d.one(), &beta);
                if let Some((q, rem)) = inv_beta
                    && d.is_zero(&rem)
                {
                    r = r.mul_scalar(&q);
                }
            }

            lcs.push((b_cur.lcoeff(), b_deg));

            a_cur = b_cur;
            b_cur = r;
        }

        // If the last non-zero polynomial is not constant, the GCD is
        // non-trivial and the resultant is zero.
        if let Some(last_deg) = b_cur.degree()
            && last_deg > 0
        {
            return d.zero();
        }
        // b_cur is now zero; check if a_cur is a non-constant GCD.
        if a_cur.degree().unwrap_or(0) > 0 {
            return d.zero();
        }

        // Compute resultant from PRS using the fundamental theorem.
        lcs.push((a_cur.lcoeff(), 0));

        let mut rho = d.one();
        let mut den = d.one();

        for k in 1..lcs.len() {
            let deg_k_prev = lcs[k - 1].1 as i64;
            let deg_k = lcs[k].1 as i64;
            #[allow(unused_variables)]
            let deg_k_next = if k + 1 < lcs.len() {
                lcs[k + 1].1 as i64
            } else {
                0
            };

            let mut exponent: i64 = deg_k_prev - deg_k;
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

        let result = d.div_rem(&rho, &den).unwrap().0;

        // Sign correction for swapping.
        if swapped && (deg_a * deg_b) % 2 == 1 {
            d.neg(&result)
        } else {
            result
        }
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
}
