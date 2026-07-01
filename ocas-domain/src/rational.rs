//! Rational domain implementation.
//!
//! Supports arbitrary-precision rational numbers. The default build uses
//! [`num_rational::BigRational`]. When the `gmp` feature is enabled, the
//! implementation moves to `gmp_backend`.

#[cfg(not(feature = "gmp"))]
use num_rational::BigRational;
#[cfg(not(feature = "gmp"))]
use num_traits::{One, Zero};

#[cfg(not(feature = "gmp"))]
use crate::domain::{Domain, EuclideanDomain};

/// The rational number domain.
///
/// # Example
///
/// ```
/// use ocas_domain::{Domain, Rational, RationalDomain};
///
/// let domain = RationalDomain;
/// let a = Rational::new(1, 2);
/// let b = Rational::new(1, 3);
/// let sum = domain.add(&a, &b);
/// assert_eq!(sum, Rational::new(5, 6));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RationalDomain;

#[cfg(not(feature = "gmp"))]
/// Arbitrary-precision rational number backed by `num-rational`.
///
/// # Example
///
/// ```
/// use ocas_domain::Rational;
///
/// let a = Rational::new(1, 2);
/// let b = Rational::new(3, 4);
/// let sum = a.inner() + b.inner();
/// assert_eq!(sum.to_string(), "5/4");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Rational(BigRational);

#[cfg(not(feature = "gmp"))]
impl Rational {
    /// Create a rational number from a numerator and denominator.
    pub fn new(numer: i64, denom: i64) -> Self {
        Self(BigRational::from_integer(numer.into()) / BigRational::from_integer(denom.into()))
    }

    /// Access the underlying `BigRational`.
    pub fn inner(&self) -> &BigRational {
        &self.0
    }
}

#[cfg(not(feature = "gmp"))]
impl Domain for RationalDomain {
    type Element = Rational;

    fn zero(&self) -> Self::Element {
        Rational(BigRational::zero())
    }

    fn one(&self) -> Self::Element {
        Rational(BigRational::one())
    }

    fn pow(&self, a: &Self::Element, n: u64) -> Self::Element {
        Rational(a.0.clone().pow(n as i32))
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Rational(a.0.clone() + b.0.clone())
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Rational(a.0.clone() - b.0.clone())
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        Rational(-a.0.clone())
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Rational(a.0.clone() * b.0.clone())
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        if b.0.is_zero() {
            return None;
        }
        Some(Rational(a.0.clone() / b.0.clone()))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.0.is_zero() {
            return None;
        }
        Some(Rational(a.0.clone().recip()))
    }
}

#[cfg(not(feature = "gmp"))]
impl EuclideanDomain for RationalDomain {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        let q = self.div(a, b)?;
        Some((q, self.zero()))
    }
}

#[cfg(not(feature = "gmp"))]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_addition() {
        let domain = RationalDomain;
        let a = Rational::new(1, 2);
        let b = Rational::new(1, 3);
        let sum = domain.add(&a, &b);
        assert_eq!(sum, Rational::new(5, 6));
    }

    #[test]
    fn rational_division() {
        let domain = RationalDomain;
        let a = Rational::new(2, 3);
        let b = Rational::new(4, 5);
        let q = domain.div(&a, &b).unwrap();
        assert_eq!(q, Rational::new(5, 6));
    }

    #[test]
    fn rational_division_by_zero() {
        let domain = RationalDomain;
        let a = Rational::new(1, 2);
        let b = domain.zero();
        assert!(domain.div(&a, &b).is_none());
    }
}
