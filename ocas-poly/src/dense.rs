//! Dense univariate polynomial implementation.
//!
//! A [`DenseUnivariatePolynomial`] stores all coefficients from the constant
//! term up to the leading coefficient in a contiguous vector. This is well
//! suited for univariate arithmetic with moderate degree.

use ocas_domain::Domain;

/// A dense univariate polynomial with coefficients in a domain `D`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseUnivariatePolynomial<D: Domain> {
    /// Coefficients from constant term upward. Trailing zeros are removed so
    /// the zero polynomial is represented by an empty vector.
    coeffs: Vec<D::Element>,
    /// The coefficient domain. Stored in the polynomial so all operations can
    /// access it without passing it explicitly.
    domain: D,
}

impl<D: Domain> DenseUnivariatePolynomial<D> {
    /// Create the zero polynomial over `domain`.
    pub fn new(domain: D) -> Self {
        Self {
            coeffs: Vec::new(),
            domain,
        }
    }

    /// Create a polynomial from a vector of coefficients `[a0, a1, ..., an]`.
    ///
    /// Trailing zero coefficients are stripped automatically.
    pub fn from_coeffs(domain: D, coeffs: Vec<D::Element>) -> Self {
        let mut poly = Self { coeffs, domain };
        poly.trim_trailing_zeros();
        poly
    }

    /// Return a reference to the coefficient domain.
    pub fn domain(&self) -> &D {
        &self.domain
    }

    /// Return the coefficients from constant term upward.
    pub fn coeffs(&self) -> &[D::Element] {
        &self.coeffs
    }

    /// Return whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Return the degree of the polynomial, or `None` for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// Return the coefficient of `x^n`, or `None` if the term is absent.
    pub fn coeff(&self, n: usize) -> Option<&D::Element> {
        self.coeffs.get(n)
    }

    /// Return the leading coefficient, or `None` for the zero polynomial.
    pub fn leading_coeff(&self) -> Option<&D::Element> {
        self.coeffs.last()
    }

    /// Return the zero polynomial with the same domain.
    pub fn zero(&self) -> Self {
        Self::new(self.domain.clone())
    }

    /// Return the constant polynomial `1` over the same domain.
    pub fn one(&self) -> Self {
        Self::from_coeffs(self.domain.clone(), vec![self.domain.one()])
    }

    /// Return the negation of this polynomial.
    pub fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter().map(|c| self.domain.neg(c)).collect();
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Add another polynomial.
    pub fn add(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(len);
        let zero = self.domain.zero();
        for i in 0..len {
            let a = self.coeffs.get(i).unwrap_or(&zero);
            let b = other.coeffs.get(i).unwrap_or(&zero);
            coeffs.push(self.domain.add(a, b));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Subtract another polynomial.
    pub fn sub(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(len);
        let zero = self.domain.zero();
        for i in 0..len {
            let a = self.coeffs.get(i).unwrap_or(&zero);
            let b = other.coeffs.get(i).unwrap_or(&zero);
            coeffs.push(self.domain.sub(a, b));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Multiply by a scalar coefficient.
    pub fn mul_scalar(&self, scalar: &D::Element) -> Self {
        if self.domain.is_zero(scalar) {
            return self.zero();
        }
        let coeffs = self
            .coeffs
            .iter()
            .map(|c| self.domain.mul(c, scalar))
            .collect();
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Multiply two polynomials.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return self.zero();
        }
        let mut coeffs = vec![self.domain.zero(); self.coeffs.len() + other.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in other.coeffs.iter().enumerate() {
                let prod = self.domain.mul(a, b);
                coeffs[i + j] = self.domain.add(&coeffs[i + j], &prod);
            }
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    fn trim_trailing_zeros(&mut self) {
        while let Some(last) = self.coeffs.last() {
            if self.domain.is_zero(last) {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{FiniteField, IntegerDomain, Rational, RationalDomain};

    #[test]
    fn zero_polynomial_has_no_degree() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::new(domain);
        assert!(p.is_zero());
        assert_eq!(p.degree(), None);
    }

    #[test]
    fn degree_and_coeffs() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![3.into(), 0.into(), 2.into(), 0.into()],
        );
        assert_eq!(p.degree(), Some(2));
        assert_eq!(p.coeff(0).cloned(), Some(3.into()));
        assert_eq!(p.coeff(2).cloned(), Some(2.into()));
        assert_eq!(p.coeff(3), None);
        assert_eq!(p.leading_coeff().cloned(), Some(2.into()));
    }

    #[test]
    fn trailing_zeros_are_trimmed() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![1.into(), 2.into(), 0.into(), 0.into()],
        );
        assert_eq!(p.degree(), Some(1));
        assert_eq!(p.coeffs().len(), 2);
    }

    #[test]
    fn add_polynomials() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let sum = a.add(&b);
        assert_eq!(sum.coeffs().to_vec(), vec![4.into(), 2.into(), 4.into()]);
    }

    #[test]
    fn sub_polynomials() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let diff = a.sub(&b);
        assert_eq!(
            diff.coeffs().to_vec(),
            vec![(-2).into(), 2.into(), (-4).into()]
        );
    }

    #[test]
    fn mul_polynomials() {
        let domain = IntegerDomain;
        // (1 + 2x) * (3 + 4x^2) = 3 + 6x + 4x^2 + 8x^3
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let prod = a.mul(&b);
        assert_eq!(
            prod.coeffs().to_vec(),
            vec![3.into(), 6.into(), 4.into(), 8.into()]
        );
    }

    #[test]
    fn mul_by_zero_yields_zero() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let zero = DenseUnivariatePolynomial::new(domain);
        let prod = a.mul(&zero);
        assert!(prod.is_zero());
    }

    #[test]
    fn rational_polynomial_multiplication() {
        let domain = RationalDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(1, 2), Rational::new(1, 1)],
        );
        let b = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(2, 1), Rational::new(1, 1)],
        );
        let prod = a.mul(&b);
        // (1/2 + x) * (2 + x) = 1 + (5/2)x + x^2
        assert_eq!(prod.coeff(0).cloned(), Some(Rational::new(1, 1)));
        assert_eq!(prod.coeff(1).cloned(), Some(Rational::new(5, 2)));
        assert_eq!(prod.coeff(2).cloned(), Some(Rational::new(1, 1)));
    }

    #[test]
    fn finite_field_polynomial_arithmetic() {
        let domain = FiniteField::new(num_bigint::BigInt::from(7));
        let a = DenseUnivariatePolynomial::from_coeffs(
            domain.clone(),
            vec![domain.element(3), domain.element(1)],
        );
        let b = DenseUnivariatePolynomial::from_coeffs(
            domain.clone(),
            vec![domain.element(2), domain.element(0), domain.element(1)],
        );
        let prod = a.mul(&b);
        // (3 + x) * (2 + x^2) = 6 + 2x + 3x^2 + x^3  (mod 7)
        assert_eq!(prod.coeff(0).cloned(), Some(domain.element(6)));
        assert_eq!(prod.coeff(1).cloned(), Some(domain.element(2)));
        assert_eq!(prod.coeff(2).cloned(), Some(domain.element(3)));
        assert_eq!(prod.coeff(3).cloned(), Some(domain.element(1)));
    }
}
