//! Integer and rational types backed by [`rug`] when the `gmp` feature is enabled.
//!
//! When `gmp` is enabled, this module defines the public [`Integer`](crate::Integer) and
//! [`Rational`](crate::Rational) types wrapping [`rug::Integer`] and [`rug::Rational`].
//! This module is not compiled otherwise.

#![cfg(feature = "gmp")]

use rug::Integer as RugInteger;
use rug::Rational as RugRational;
use std::ops::{Add, Div, Mul, Rem, Sub};

use crate::domain::{Domain, EuclideanDomain};
use crate::integer::IntegerDomain;
use crate::rational::RationalDomain;

/// Arbitrary-precision integer backed by `rug::Integer`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

impl From<num_bigint::BigInt> for Integer {
    fn from(value: num_bigint::BigInt) -> Self {
        use std::str::FromStr;
        Self(
            rug::Integer::from_str(&value.to_string())
                .expect("BigInt to rug::Integer conversion should never fail"),
        )
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

    /// Convert to a `BigInt` regardless of the backend.
    pub fn to_bigint(&self) -> num_bigint::BigInt {
        use std::str::FromStr;
        num_bigint::BigInt::from_str(&self.0.to_string())
            .expect("rug::Integer to BigInt conversion should never fail")
    }

    /// Raise to a `u32` power.
    pub fn pow_u32(&self, exp: u32) -> Self {
        use rug::ops::Pow;
        Integer(self.0.clone().pow(exp))
    }

    /// Modular exponentiation: `self^exp mod modulus`.
    pub fn modpow(&self, exp: &Integer, modulus: &Integer) -> Integer {
        Integer(
            self.0
                .clone()
                .pow_mod(&exp.0, &modulus.0)
                .expect("modpow: modulus must be positive"),
        )
    }

    /// Floor modulo: result `r` satisfies `0 ≤ r < |modulus|` for positive modulus.
    pub fn mod_floor(&self, modulus: &Integer) -> Integer {
        let r = self.clone() % modulus.clone();
        if r.is_negative() { r + modulus } else { r }
    }

    /// Division with remainder: `(quotient, remainder)`.
    pub fn div_rem(&self, other: &Integer) -> (Integer, Integer) {
        let (q, r) = self.0.clone().div_rem(other.0.clone());
        (Integer(q), Integer(r))
    }

    /// Returns `true` if the value is even.
    pub fn is_even(&self) -> bool {
        self.0.is_even()
    }

    /// Returns `true` if the value is negative.
    pub fn is_negative(&self) -> bool {
        self.0.is_negative()
    }

    /// Returns `true` if the value is zero.
    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Returns `true` if the value is one.
    pub fn is_one(&self) -> bool {
        self.0 == 1
    }

    /// Absolute value.
    pub fn abs(&self) -> Integer {
        Integer(self.0.clone().abs())
    }

    /// Integer square root (floor).
    pub fn sqrt(&self) -> Integer {
        Integer(self.0.clone().sqrt())
    }
}

// ---------------------------------------------------------------------------
// Arithmetic operators — forward to the inner rug::Integer.
// Note: rug's &Integer op &Integer returns incomplete types, so we clone
// to always use owned arithmetic.
// ---------------------------------------------------------------------------

macro_rules! impl_gmp_int_op_owned {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait for Integer {
            type Output = Integer;
            fn $method(self, rhs: Integer) -> Integer {
                Integer(self.0 $op rhs.0)
            }
        }
        impl $trait<&Integer> for Integer {
            type Output = Integer;
            fn $method(self, rhs: &Integer) -> Integer {
                Integer(self.0 $op rhs.0.clone())
            }
        }
        impl $trait<Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: Integer) -> Integer {
                Integer(self.0.clone() $op rhs.0)
            }
        }
        impl $trait<&Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: &Integer) -> Integer {
                Integer(self.0.clone() $op rhs.0.clone())
            }
        }
    };
}

impl_gmp_int_op_owned!(Add, add, +);
impl_gmp_int_op_owned!(Sub, sub, -);
impl_gmp_int_op_owned!(Mul, mul, *);
impl_gmp_int_op_owned!(Div, div, /);
impl_gmp_int_op_owned!(Rem, rem, %);

impl std::ops::Neg for Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        Integer(-self.0)
    }
}
impl std::ops::Neg for &Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        Integer(-self.0.clone())
    }
}

impl std::ops::ShrAssign<u32> for Integer {
    fn shr_assign(&mut self, shift: u32) {
        self.0 >>= shift;
    }
}
impl std::ops::Shr<u32> for Integer {
    type Output = Integer;
    fn shr(self, shift: u32) -> Integer {
        Integer(self.0 >> shift)
    }
}
impl std::ops::Shr<u32> for &Integer {
    type Output = Integer;
    fn shr(self, shift: u32) -> Integer {
        Integer(self.0.clone() >> shift)
    }
}

impl std::ops::AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, rhs: &Integer) {
        self.0 += &rhs.0;
    }
}
impl std::ops::SubAssign<&Integer> for Integer {
    fn sub_assign(&mut self, rhs: &Integer) {
        self.0 -= &rhs.0;
    }
}
impl std::ops::MulAssign<&Integer> for Integer {
    fn mul_assign(&mut self, rhs: &Integer) {
        self.0 *= &rhs.0;
    }
}
impl std::ops::DivAssign<&Integer> for Integer {
    fn div_assign(&mut self, rhs: &Integer) {
        self.0 /= &rhs.0;
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

    /// Create a rational number from an integer (denominator = 1).
    pub fn from_integer(n: Integer) -> Self {
        Self(RugRational::from(n.0))
    }

    /// Numerator as an [`Integer`].
    pub fn numer(&self) -> Integer {
        Integer(self.0.numer().clone())
    }

    /// Denominator as an [`Integer`].
    pub fn denom(&self) -> Integer {
        Integer(self.0.denom().clone())
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
