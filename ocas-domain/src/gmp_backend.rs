//! Integer and rational types backed by [`rug`] when the `gmp` feature is enabled.
//!
//! When `gmp` is enabled, this module defines the public [`Integer`](crate::Integer) and
//! [`Rational`](crate::Rational) types wrapping [`rug::Integer`] and [`rug::Rational`].
//! This module is not compiled otherwise.

#![cfg(feature = "gmp")]

use rug::Integer as RugInteger;
use rug::Rational as RugRational;

use crate::domain::{Domain, EuclideanDomain};
use crate::integer::IntegerDomain;
use crate::rational::RationalDomain;

/// Arbitrary-precision integer backed by `rug::Integer`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Integer(RugInteger);

/// Arbitrary-precision rational number backed by `rug::Rational`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Rational(RugRational);

impl std::fmt::Display for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Domain for IntegerDomain {
    type Element = Integer;

    fn zero(&self) -> Self::Element {
        Integer(RugInteger::from(0))
    }

    fn one(&self) -> Self::Element {
        Integer(RugInteger::from(1))
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Integer(a.0.clone() + b.0.clone())
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Integer(a.0.clone() - b.0.clone())
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        Integer(-a.0.clone())
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Integer(a.0.clone() * b.0.clone())
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        if b.0 == 0 {
            return None;
        }
        let (q, r) = a.0.clone().div_rem(b.0.clone());
        if r == 0 { Some(Integer(q)) } else { None }
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.0 == 1 {
            Some(self.one())
        } else if a.0 == -1 {
            Some(Integer(RugInteger::from(-1)))
        } else {
            None
        }
    }
}

impl EuclideanDomain for IntegerDomain {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        if b.0 == 0 {
            return None;
        }
        let (q, r) = a.0.clone().div_rem(b.0.clone());
        Some((Integer(q), Integer(r)))
    }
}

impl Domain for RationalDomain {
    type Element = Rational;

    fn zero(&self) -> Self::Element {
        Rational(RugRational::from(0))
    }

    fn one(&self) -> Self::Element {
        Rational(RugRational::from(1))
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
        if b.0 == 0 {
            return None;
        }
        Some(Rational(a.0.clone() / b.0.clone()))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.0 == 0 {
            return None;
        }
        Some(Rational(a.0.clone().recip()))
    }
}

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

impl From<RugInteger> for Integer {
    fn from(value: RugInteger) -> Self {
        Self(value)
    }
}

impl From<RugRational> for Rational {
    fn from(value: RugRational) -> Self {
        Self(value)
    }
}

impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Self::new(value)
    }
}

impl Integer {
    /// Create an integer from a machine integer or another `Into<RugInteger>`.
    pub fn new<T: Into<RugInteger>>(value: T) -> Self {
        Self(value.into())
    }

    /// Access the underlying [`rug::Integer`].
    pub fn inner(&self) -> &RugInteger {
        &self.0
    }
}

impl Rational {
    /// Create a rational number from a numerator and denominator.
    pub fn new(numer: i64, denom: i64) -> Self {
        Self(RugRational::from((numer, denom)))
    }

    /// Access the underlying [`rug::Rational`].
    pub fn inner(&self) -> &RugRational {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_addition() {
        let domain = IntegerDomain;
        let a = Integer::from(3);
        let b = Integer::from(5);
        assert_eq!(domain.add(&a, &b), Integer::from(8));
    }

    #[test]
    fn integer_div_exact() {
        let domain = IntegerDomain;
        let a = Integer::from(10);
        let b = Integer::from(3);
        assert!(domain.div(&a, &b).is_none());
        let c = Integer::from(2);
        assert_eq!(domain.div(&a, &c), Some(Integer::from(5)));
    }

    #[test]
    fn integer_div_rem() {
        let domain = IntegerDomain;
        let a = Integer::from(17);
        let b = Integer::from(5);
        let (q, r) = domain.div_rem(&a, &b).unwrap();
        assert_eq!(q, Integer::from(3));
        assert_eq!(r, Integer::from(2));
    }

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
