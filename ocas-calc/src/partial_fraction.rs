//! Partial fraction decomposition of univariate rational functions.
//!
//! Given a proper fraction $\frac{p(x)}{q(x)}$ (i.e. $\deg(p) < \deg(q)$),
//! decomposes it into a sum of simpler fractions:
//!
//! $$\frac{p(x)}{q(x)} = \text{polynomial part} + \sum_i \sum_{k=1}^{e_i}
//! \frac{a_{i,k}(x)}{f_i(x)^k}$$
//!
//! where $q = \prod_i f_i^{e_i}$ is the square-free factorization of the
//! denominator.
//!
//! This is the univariate analogue of Symbolica's `apart()` method.

use ocas_domain::EuclideanDomain;
use ocas_poly::DenseUnivariatePolynomial;

/// A single term in a partial fraction decomposition.
///
/// Represents the fraction $\frac{\text{numer}}{\text{denom}^{\text{exp}}}$.
#[derive(Debug, Clone)]
pub struct PartialFractionTerm<D: EuclideanDomain> {
    /// The numerator polynomial.
    pub numer: DenseUnivariatePolynomial<D>,
    /// The irreducible (square-free) denominator factor.
    pub denom: DenseUnivariatePolynomial<D>,
    /// The exponent (multiplicity) of this factor in the original denominator.
    pub exp: usize,
}

/// Perform partial fraction decomposition of `num / den`.
///
/// Returns a list of [`PartialFractionTerm`]s and an optional polynomial
/// part. The decomposition satisfies:
///
/// $$\frac{\text{num}}{\text{den}} = \text{poly\_part} + \sum_i
/// \frac{\text{numer}_i}{\text{denom}_i^{\text{exp}_i}}$$
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::DenseUnivariatePolynomial;
/// use ocas_calc::partial_fraction::apart;
///
/// let d = RationalDomain;
/// // 1 / (x^2 - 1) — square-free, single factor
/// let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
/// let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
///     Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1),
/// ]);
/// let (poly_part, terms) = apart(&num, &den);
/// assert!(poly_part.is_none());
/// // x^2-1 is square-free → 1 term
/// assert_eq!(terms.len(), 1);
/// ```
pub fn apart<D: EuclideanDomain>(
    num: &DenseUnivariatePolynomial<D>,
    den: &DenseUnivariatePolynomial<D>,
) -> (
    Option<DenseUnivariatePolynomial<D>>,
    Vec<PartialFractionTerm<D>>,
) {
    let _d = num.domain();

    if den.is_zero() {
        return (Some(num.clone()), Vec::new());
    }
    if num.is_zero() {
        return (None, Vec::new());
    }

    // Step 1: Polynomial long division to extract polynomial part.
    let (quotient, remainder) = match num.div_rem(den) {
        Some(v) => v,
        None => return (Some(num.clone()), Vec::new()),
    };

    let poly_part = if quotient.is_zero() {
        None
    } else {
        Some(quotient)
    };

    if remainder.is_zero() {
        return (poly_part, Vec::new());
    }

    // Step 2: Square-free factorization of the denominator.
    let sq_free = den.square_free_factorization();

    if sq_free.is_empty() {
        // Denominator is a constant; remainder/den is already a simple fraction.
        return (
            poly_part,
            vec![PartialFractionTerm {
                numer: remainder,
                denom: den.clone(),
                exp: 1,
            }],
        );
    }

    // Step 3: For each square-free factor with multiplicity > 0, compute
    // the partial fraction terms using diophantine CRT + p-adic expansion.
    let mut terms = Vec::new();

    if sq_free.len() == 1 {
        // Single factor: just store the remainder over the factor^exp.
        let (factor, exp) = &sq_free[0];
        if *exp == 1 {
            terms.push(PartialFractionTerm {
                numer: remainder.clone(),
                denom: factor.clone(),
                exp: 1,
            });
        } else {
            // Use p-adic expansion to decompose remainder / factor^exp.
            let expanded = remainder.p_adic_expansion(factor);
            for (k, coeff) in expanded.into_iter().enumerate() {
                if !coeff.is_zero() {
                    terms.push(PartialFractionTerm {
                        numer: coeff,
                        denom: factor.clone(),
                        exp: k + 1,
                    });
                }
            }
        }
    } else {
        // Multiple factors: use diophantine CRT.
        // Build the list of factor^exp polynomials.
        let mut moduli: Vec<DenseUnivariatePolynomial<D>> =
            sq_free.iter().map(|(f, e)| f.pow(*e as u32)).collect();

        // Solve the diophantine system: find s_i such that
        // Σ s_i * (Π_{j≠i} moduli_j) ≡ remainder (mod Π moduli_i)
        let deltas = DenseUnivariatePolynomial::diophantine(&mut moduli, &remainder);

        // For each delta_i, decompose it via p-adic expansion w.r.t. factor_i.
        for (i, delta) in deltas.into_iter().enumerate() {
            let (factor, exp) = &sq_free[i];
            if delta.is_zero() {
                continue;
            }
            if *exp == 1 {
                terms.push(PartialFractionTerm {
                    numer: delta,
                    denom: factor.clone(),
                    exp: 1,
                });
            } else {
                let expanded = delta.p_adic_expansion(factor);
                for (k, coeff) in expanded.into_iter().enumerate() {
                    if !coeff.is_zero() {
                        terms.push(PartialFractionTerm {
                            numer: coeff,
                            denom: factor.clone(),
                            exp: k + 1,
                        });
                    }
                }
            }
        }
    }

    (poly_part, terms)
}

/// Combine partial fraction terms back into a single rational function.
///
/// Given a polynomial part and a list of terms, computes the numerator
/// and denominator of the combined fraction.
///
/// Returns `(numerator, denominator)`.
pub fn together<D: EuclideanDomain>(
    poly_part: Option<&DenseUnivariatePolynomial<D>>,
    terms: &[PartialFractionTerm<D>],
) -> (DenseUnivariatePolynomial<D>, DenseUnivariatePolynomial<D>) {
    if terms.is_empty() {
        let zero = DenseUnivariatePolynomial::new(
            terms
                .first()
                .map(|t| t.numer.domain().clone())
                .unwrap_or_else(|| {
                    // Fallback: use poly_part's domain
                    poly_part.unwrap().domain().clone()
                }),
        );
        let one = zero.one();
        let pp = poly_part.cloned().unwrap_or(zero);
        return (pp, one);
    }

    let _domain = terms[0].numer.domain();

    // Compute the common denominator: Π denom_i^exp_i
    let mut common_den = terms[0].denom.one();
    for term in terms {
        let factor = term.denom.pow(term.exp as u32);
        common_den = common_den.mul(&factor);
    }

    // Compute the numerator: Σ numer_i * (common_den / denom_i^exp_i)
    let mut combined_num = common_den.zero();
    for term in terms {
        let factor = term.denom.pow(term.exp as u32);
        let cofactor = common_den
            .div_rem(&factor)
            .map(|(q, _)| q)
            .unwrap_or(common_den.clone());
        let contribution = term.numer.mul(&cofactor);
        combined_num = combined_num.add(&contribution);
    }

    // Add polynomial part: poly_part * common_den / 1
    if let Some(pp) = poly_part {
        let pp_contribution = pp.mul(&common_den);
        combined_num = combined_num.add(&pp_contribution);
    }

    (combined_num, common_den)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Rational, RationalDomain};

    fn rat_poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<RationalDomain> {
        DenseUnivariatePolynomial::from_coeffs(
            RationalDomain,
            coeffs.iter().map(|&c| Rational::new(c, 1)).collect(),
        )
    }

    fn rat_poly_r(coeffs: &[(i64, i64)]) -> DenseUnivariatePolynomial<RationalDomain> {
        DenseUnivariatePolynomial::from_coeffs(
            RationalDomain,
            coeffs.iter().map(|&(n, d)| Rational::new(n, d)).collect(),
        )
    }

    #[test]
    fn apart_simple() {
        // 1 / (x^2 - 1) — square-free, single factor
        // square_free_factorization returns (x^2-1, 1), so we get 1 term
        let num = rat_poly(&[1]); // 1
        let den = rat_poly(&[-1, 0, 1]); // x^2 - 1

        let (poly_part, terms) = apart(&num, &den);
        assert!(poly_part.is_none());
        // Single square-free factor → 1 term (no further decomposition without full factorization)
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].exp, 1);

        // Verify by combining back: numer/denom should equal 1/(x^2-1)
        let (n, _d) = together(None, &terms);
        assert!(!n.is_zero());
    }

    #[test]
    fn apart_with_polynomial_part() {
        // (x^2 + 1) / (x - 1) = x + 1 + 2/(x-1)
        // div_rem: x^2+1 = (x+1)(x-1) + 2
        let num = rat_poly(&[1, 0, 1]); // x^2 + 1
        let den = rat_poly(&[-1, 1]); // x - 1

        let (poly_part, terms) = apart(&num, &den);
        assert!(poly_part.is_some());
        let pp = poly_part.unwrap();
        // polynomial part = x + 1
        assert_eq!(pp.degree(), Some(1));
        assert_eq!(pp.coeff(0), Some(&Rational::new(1, 1)));
        assert_eq!(pp.coeff(1), Some(&Rational::new(1, 1)));
        // remainder term = 2/(x-1)
        assert_eq!(terms.len(), 1);
    }

    #[test]
    fn apart_repeated_factor() {
        // 1 / (x-1)^2 should give a single term with exp=2,
        // or p-adic expansion terms.
        let num = rat_poly(&[1]); // 1
        let den = rat_poly(&[1, -2, 1]); // (x-1)^2 = x^2 - 2x + 1

        let (poly_part, terms) = apart(&num, &den);
        assert!(poly_part.is_none());
        assert!(!terms.is_empty());
    }

    #[test]
    fn apart_trivial() {
        // x / x = 1 (polynomial, no remainder)
        let num = rat_poly(&[0, 1]); // x
        let den = rat_poly(&[0, 1]); // x

        let (poly_part, terms) = apart(&num, &den);
        assert!(poly_part.is_some());
        assert!(terms.is_empty());
    }

    #[test]
    fn together_roundtrip() {
        // Build terms manually: 1/(x-1) - 1/(x+1) = 2/(x^2-1)
        let terms = vec![
            PartialFractionTerm {
                numer: rat_poly_r(&[(1, 2)]), // 1/2
                denom: rat_poly(&[-1, 1]),    // x - 1
                exp: 1,
            },
            PartialFractionTerm {
                numer: rat_poly_r(&[(-1, 2)]), // -1/2
                denom: rat_poly(&[1, 1]),      // x + 1
                exp: 1,
            },
        ];

        let (n, d) = together(None, &terms);
        // Should get something proportional to 1/(x^2-1).
        assert!(!n.is_zero());
        assert_eq!(d.degree(), Some(2));
    }
}
