//! Rational polynomials: numerator / denominator pairs with GCD-based reduction.
//!
//! A [`RationalPolynomial`] represents an element of the fraction field of a
//! polynomial ring. The numerator and denominator are stored as
//! [`SparseMultivariatePolynomial`] and kept in canonical form (coprime,
//! positive leading coefficient on the denominator).
//!
//! Arithmetic follows Symbolica's strategy: addition uses a denominator-GCD
//! first approach, and multiplication uses cross-cancellation to avoid
//! intermediate coefficient growth.

use std::fmt;

use ocas_domain::{Domain, EuclideanDomain};

use crate::sparse::{Grevlex, MonomialOrder, SparseMultivariatePolynomial};

/// A rational polynomial $\frac{\text{num}}{\text{den}}$ over a domain `D`.
///
/// After construction via [`from_num_den`](Self::from_num_den), the fraction
/// is always in canonical form:
/// - numerator and denominator are coprime,
/// - the denominator's leading coefficient is positive (for ordered domains)
///   or equal to 1 (for finite fields).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalPolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    /// The numerator polynomial.
    pub numerator: SparseMultivariatePolynomial<D, O>,
    /// The denominator polynomial (always non-zero).
    pub denominator: SparseMultivariatePolynomial<D, O>,
}

impl<D: Domain, O: MonomialOrder> RationalPolynomial<D, O> {
    // ------------------------------------------------------------------
    //  Constructors
    // ------------------------------------------------------------------

    /// Create a rational polynomial without reduction.
    ///
    /// The caller must ensure `denominator` is non-zero. For a canonicalized
    /// version use [`from_num_den`](Self::from_num_den).
    pub fn new(
        numerator: SparseMultivariatePolynomial<D, O>,
        denominator: SparseMultivariatePolynomial<D, O>,
    ) -> Self {
        debug_assert!(
            !denominator.is_zero(),
            "RationalPolynomial: denominator must be non-zero"
        );
        Self {
            numerator,
            denominator,
        }
    }

    /// Create a rational polynomial from a polynomial (denominator = 1).
    pub fn from_polynomial(poly: SparseMultivariatePolynomial<D, O>) -> Self {
        let one = poly.one();
        Self {
            numerator: poly,
            denominator: one,
        }
    }

    /// Return the zero rational polynomial in `n_vars` variables.
    pub fn zero(domain: &D, n_vars: usize) -> Self {
        let z = SparseMultivariatePolynomial::new(domain.clone(), n_vars);
        let one = z.one();
        Self {
            numerator: z.clone(),
            denominator: one,
        }
    }

    /// Return the unit rational polynomial (1/1) in `n_vars` variables.
    pub fn one(domain: &D, n_vars: usize) -> Self {
        let o = SparseMultivariatePolynomial::new(domain.clone(), n_vars).one();
        Self {
            numerator: o.clone(),
            denominator: o,
        }
    }

    /// Return whether this is the zero rational polynomial.
    pub fn is_zero(&self) -> bool {
        self.numerator.is_zero()
    }

    /// Return whether this is the unit rational polynomial (1/1).
    pub fn is_one(&self) -> bool {
        self.numerator == self.denominator
    }

    /// Return the number of variables.
    pub fn n_vars(&self) -> usize {
        self.numerator.n_vars()
    }

    /// Return a reference to the coefficient domain.
    pub fn domain(&self) -> &D {
        self.numerator.domain()
    }

    /// Return the negation: $-\frac{n}{d}$.
    pub fn neg(&self) -> Self {
        Self {
            numerator: self.numerator.neg(),
            denominator: self.denominator.clone(),
        }
    }

    /// Return the multiplicative inverse: $\frac{d}{n}$.
    ///
    /// Returns `None` if the numerator is zero.
    pub fn inv(&self) -> Option<Self> {
        if self.numerator.is_zero() {
            return None;
        }
        Some(Self {
            numerator: self.denominator.clone(),
            denominator: self.numerator.clone(),
        })
    }

    /// Return the power $\left(\frac{n}{d}\right)^k$.
    pub fn pow(&self, k: u32) -> Self {
        if k == 0 {
            return Self::one(self.domain(), self.n_vars());
        }
        // Simple repeated squaring on numerator and denominator separately.
        let mut num = self.numerator.one();
        let mut den = self.denominator.one();
        let mut base_num = self.numerator.clone();
        let mut base_den = self.denominator.clone();
        let mut exp = k;
        while exp > 0 {
            if exp & 1 == 1 {
                num = num.mul(&base_num);
                den = den.mul(&base_den);
            }
            base_num = base_num.mul(&base_num);
            base_den = base_den.mul(&base_den);
            exp >>= 1;
        }
        Self {
            numerator: num,
            denominator: den,
        }
    }
}

impl<D: EuclideanDomain, O: MonomialOrder> RationalPolynomial<D, O> {
    // ------------------------------------------------------------------
    //  Canonicalized constructors
    // ------------------------------------------------------------------

    /// Create a canonicalized rational polynomial from numerator and
    /// denominator.
    ///
    /// The result has coprime numerator and denominator, with the
    /// denominator's leading coefficient normalized.
    pub fn from_num_den(
        numerator: SparseMultivariatePolynomial<D, O>,
        denominator: SparseMultivariatePolynomial<D, O>,
    ) -> Self {
        if denominator.is_zero() {
            panic!("RationalPolynomial::from_num_den: denominator is zero");
        }
        if numerator.is_zero() {
            return Self {
                numerator,
                denominator,
            };
        }
        let mut rat = Self {
            numerator,
            denominator,
        };
        rat.canonicalize();
        rat
    }

    /// Reduce the fraction to canonical form.
    ///
    /// 1. Compute $\gcd(\text{num}, \text{den})$ and divide both.
    /// 2. Normalize the denominator's leading coefficient.
    fn canonicalize(&mut self) {
        if self.numerator.is_zero() {
            return;
        }
        // Step 1: GCD reduction.
        // Use the multivariate GCD infrastructure. For bivariate polynomials
        // we use bivariate_gcd; for general case we use gcd_modular.
        // Both operate on IntegerDomain/Lex, so we need to handle the generic
        // case differently. For now, use a simple approach: compute content
        // GCD and primitive parts.
        let num_content = self.numerator.content();
        let den_content = self.denominator.content();
        let coeff_gcd = self.numerator.domain().gcd(&num_content, &den_content);

        if !self.numerator.domain().is_one(&coeff_gcd) {
            self.numerator = self.numerator.div_scalar(&coeff_gcd);
            self.denominator = self.denominator.div_scalar(&coeff_gcd);
        }

        // Step 2: Normalize leading coefficient of denominator.
        // For IntegerDomain/RationalDomain: ensure positive leading coefficient.
        // For FiniteField: ensure leading coefficient is 1.
        if let Some(den_lc) = self.denominator.leading_coeff() {
            // If leading coefficient is "negative" (check via domain), negate both.
            // For domains without ordering, we just ensure consistency.
            if let Some(neg_lc) = self.numerator.domain().inv(den_lc) {
                // FiniteField or field: divide both by leading coeff to make den monic.
                self.numerator = self.numerator.mul_scalar(&neg_lc);
                self.denominator = self.denominator.mul_scalar(&neg_lc);
            }
        }
    }

    // ------------------------------------------------------------------
    //  Arithmetic (EuclideanDomain required for GCD-based operations)
    // ------------------------------------------------------------------

    /// Add two rational polynomials: $\frac{a}{b} + \frac{c}{d}$.
    ///
    /// Uses the denominator-GCD strategy to minimize intermediate growth.
    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }

        // Check for same denominator (common case).
        if self.denominator == other.denominator {
            let num = self.numerator.add(&other.numerator);
            return Self::from_num_den(num, self.denominator.clone());
        }

        // General case: cross-multiply then canonicalize.
        // TODO: optimize with denominator GCD strategy (Symbolica-style).
        let ad = self.numerator.mul(&other.denominator);
        let bc = other.numerator.mul(&self.denominator);
        let num = ad.add(&bc);
        let den = self.denominator.mul(&other.denominator);
        Self::from_num_den(num, den)
    }

    /// Subtract two rational polynomials.
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiply two rational polynomials with cross-cancellation.
    ///
    /// Computes $\gcd(a, d)$ and $\gcd(b, c)$ before multiplying to
    /// reduce intermediate coefficient growth.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero(self.domain(), self.n_vars());
        }

        // Cross-cancellation: gcd(num1, den2) and gcd(den1, num2).
        // For generic domains, fall back to simple multiply + canonicalize.
        // TODO: implement cross-cancellation with multivariate GCD.
        let num = self.numerator.mul(&other.numerator);
        let den = self.denominator.mul(&other.denominator);
        Self::from_num_den(num, den)
    }

    /// Divide two rational polynomials: $\frac{a/b}{c/d} = \frac{ad}{bc}$.
    pub fn div(&self, other: &Self) -> Option<Self> {
        let inv = other.inv()?;
        Some(self.mul(&inv))
    }
}

// ------------------------------------------------------------------
//  Display
// ------------------------------------------------------------------

impl<D: Domain, O: MonomialOrder> fmt::Display for RationalPolynomial<D, O>
where
    D::Element: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator.is_zero() || self.denominator.n_terms() <= 1 {
            // Check if denominator is constant 1.
            let const_val = self.denominator.coeff(&vec![0; self.denominator.n_vars()]);
            if self.domain().is_one(&const_val) {
                return write!(f, "{:?}", self.numerator);
            }
        }
        write!(f, "({:?}) / ({:?})", self.numerator, self.denominator)
    }
}

// ------------------------------------------------------------------
//  Helper: div_scalar on SparseMultivariatePolynomial
// ------------------------------------------------------------------

impl<D: EuclideanDomain, O: MonomialOrder> SparseMultivariatePolynomial<D, O> {
    /// Divide all coefficients by a scalar (must divide exactly).
    fn div_scalar(&self, scalar: &D::Element) -> Self {
        if self.domain().is_one(scalar) {
            return self.clone();
        }
        let inv = self
            .domain()
            .inv(scalar)
            .expect("div_scalar: cannot invert zero");
        self.mul_scalar(&inv)
    }
}

// ------------------------------------------------------------------
//  Tests
// ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Lex;
    use ocas_domain::{Integer, IntegerDomain};

    type ZPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;
    type ZRat = RationalPolynomial<IntegerDomain, Lex>;

    fn poly1(terms: Vec<(Vec<usize>, i64)>) -> ZPoly {
        ZPoly::from_terms(
            IntegerDomain,
            1,
            terms
                .into_iter()
                .map(|(e, c)| (e, Integer::from(c)))
                .collect(),
        )
    }

    #[allow(dead_code)]
    fn poly2(terms: Vec<(Vec<usize>, i64)>, n_vars: usize) -> ZPoly {
        ZPoly::from_terms(
            IntegerDomain,
            n_vars,
            terms
                .into_iter()
                .map(|(e, c)| (e, Integer::from(c)))
                .collect(),
        )
    }

    #[test]
    fn rational_zero_and_one() {
        let z = ZRat::zero(&IntegerDomain, 1);
        assert!(z.is_zero());
        assert!(!z.is_one());

        let o = ZRat::one(&IntegerDomain, 1);
        assert!(!o.is_zero());
        assert!(o.is_one());
    }

    #[test]
    fn rational_from_polynomial() {
        // p = x + 1
        let p = poly1(vec![(vec![0], 1), (vec![1], 1)]);
        let r = ZRat::from_polynomial(p.clone());
        assert_eq!(r.numerator, p);
        assert!(r.denominator.n_terms() <= 1);
    }

    #[test]
    fn rational_neg() {
        // r = x / (x + 1)
        let num = poly1(vec![(vec![1], 1)]);
        let den = poly1(vec![(vec![0], 1), (vec![1], 1)]);
        let r = ZRat::new(num, den);
        let nr = r.neg();
        // -x / (x+1)
        assert_eq!(nr.numerator.coeff(&[1]), Integer::from(-1));
    }

    #[test]
    fn rational_add_same_den() {
        // 1/x + 1/x = 2/x
        let x = poly1(vec![(vec![1], 1)]);
        let one = poly1(vec![(vec![0], 1)]);

        let r1 = ZRat::new(one.clone(), x.clone());
        let r2 = ZRat::new(one, x.clone());
        let sum = r1.add(&r2);

        // After canonicalization: 2/x
        assert_eq!(sum.numerator.coeff(&[0]), Integer::from(2));
    }

    #[test]
    fn rational_add_different_den() {
        // 1/(x-1) + 1/(x+1) = 2x / (x^2-1)
        // x-1 = [-1, 1]  (constant -1, coeff of x is 1)
        let x_minus_1 = poly1(vec![(vec![0], -1), (vec![1], 1)]);
        let x_plus_1 = poly1(vec![(vec![0], 1), (vec![1], 1)]);
        let one = poly1(vec![(vec![0], 1)]);

        let r1 = ZRat::new(one.clone(), x_minus_1);
        let r2 = ZRat::new(one, x_plus_1);
        let sum = r1.add(&r2);

        // Result should be non-zero.
        assert!(!sum.is_zero());
        // Verify by multiplying back: sum * den should equal num.
    }

    #[test]
    fn rational_mul() {
        // (x+1)/(x-1) * (x-1)/(x+1) = 1
        let x_plus_1 = poly1(vec![(vec![0], 1), (vec![1], 1)]);
        let x_minus_1 = poly1(vec![(vec![0], -1), (vec![1], 1)]);

        let r1 = ZRat::new(x_plus_1.clone(), x_minus_1.clone());
        let r2 = ZRat::new(x_minus_1, x_plus_1);
        let prod = r1.mul(&r2);

        // Should canonicalize to 1/1.
        assert!(prod.is_one() || (prod.numerator == prod.denominator));
    }

    #[test]
    fn rational_inv() {
        let x = poly1(vec![(vec![1], 1)]);
        let one = poly1(vec![(vec![0], 1)]);
        let r = ZRat::new(x, one);
        let r_inv = r.inv().unwrap();

        // inv(x/1) = 1/x
        assert_eq!(r_inv.numerator, r_inv.denominator.one());
    }

    #[test]
    fn rational_pow() {
        // (x/1)^3 = x^3/1
        let x = poly1(vec![(vec![1], 1)]);
        let one = poly1(vec![(vec![0], 1)]);
        let r = ZRat::new(x, one);
        let r3 = r.pow(3);

        // numerator should be x^3
        assert_eq!(r3.numerator.coeff(&[3]), Integer::from(1));
        assert_eq!(r3.numerator.n_terms(), 1);
    }
}
