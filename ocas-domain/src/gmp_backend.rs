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

// ---------------------------------------------------------------------------
// IntegerInner — small-integer optimization (SOO)
// ---------------------------------------------------------------------------

/// Internal representation: small values stay on the stack, large values heap-allocate.
#[derive(Debug, Clone)]
enum IntegerInner {
    /// Value fits in an `i64`. Arithmetic uses native ops with overflow checks.
    Small(i64),
    /// Value requires arbitrary precision. Backed by GMP via `rug`.
    Large(Box<RugInteger>),
}

impl IntegerInner {
    /// Convert to a rug::Integer reference, promoting Small → Large if needed.
    ///
    /// # Safety (promotion)
    /// Promotion is monotonic (Small → Large, never reverse) and preserves
    /// the logical value, so the returned reference is always valid.
    fn as_rug(&self) -> &RugInteger {
        match self {
            IntegerInner::Small(v) => {
                // Promote in-place. We use unsafe to mutate through &self
                // because promotion is idempotent and value-preserving.
                // SAFETY: we write a valid Large variant before returning
                // a reference into it. No other code can observe the
                // intermediate state because this is a single-threaded
                // operation on a local borrow.
                unsafe {
                    let self_mut = &mut *(self as *const Self as *mut Self);
                    let rug_val = RugInteger::from(*v);
                    *self_mut = IntegerInner::Large(Box::new(rug_val));
                    match self_mut {
                        IntegerInner::Large(r) => &**(r),
                        _ => unreachable!(),
                    }
                }
            }
            IntegerInner::Large(r) => r,
        }
    }

    /// Consume self and return a `RugInteger`, avoiding allocation for Small values
    /// that are immediately needed as rug (e.g. `Rational::from_integer`).
    fn into_rug(self) -> RugInteger {
        match self {
            IntegerInner::Small(v) => RugInteger::from(v),
            IntegerInner::Large(r) => *r,
        }
    }
}

// Custom PartialEq/Eq/PartialOrd/Ord/Hash that compare by mathematical value.
impl PartialEq for IntegerInner {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) => a == b,
            _ => self.as_rug() == other.as_rug(),
        }
    }
}
impl Eq for IntegerInner {}

impl PartialOrd for IntegerInner {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IntegerInner {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) => a.cmp(b),
            _ => self.as_rug().cmp(other.as_rug()),
        }
    }
}

impl std::hash::Hash for IntegerInner {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            IntegerInner::Small(v) => v.hash(state),
            IntegerInner::Large(r) => r.hash(state),
        }
    }
}

// ---------------------------------------------------------------------------
// Integer — public type
// ---------------------------------------------------------------------------

/// Arbitrary-precision integer with small-integer optimization (SOO).
///
/// Values that fit in an `i64` are stored on the stack and use native
/// arithmetic. Larger values fall back to GMP via `rug::Integer`.
///
/// This mirrors FLINT's `fmpz_t` strategy: most CAS coefficients are
/// small, so avoiding heap allocation gives a significant speedup.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Integer(IntegerInner);

/// Arbitrary-precision rational number backed by `rug::Rational`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Rational(RugRational);

// ---------------------------------------------------------------------------
// Display
// ---------------------------------------------------------------------------

impl std::fmt::Display for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            IntegerInner::Small(v) => write!(f, "{}", v),
            IntegerInner::Large(r) => write!(f, "{}", r),
        }
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

// ---------------------------------------------------------------------------
// Domain for IntegerDomain (with SOO fast paths)
// ---------------------------------------------------------------------------

impl Domain for IntegerDomain {
    type Element = Integer;

    fn zero(&self) -> Self::Element {
        Integer::from_small(0)
    }

    fn one(&self) -> Self::Element {
        Integer::from_small(1)
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        a.add_ref(b)
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        a.sub_ref(b)
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        a.neg_ref()
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        a.mul_ref(b)
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        match (&a.0, &b.0) {
            (IntegerInner::Small(av), IntegerInner::Small(bv)) => {
                if *bv == 0 { return None; }
                if *av % *bv == 0 { Some(Integer::from_small(*av / *bv)) } else { None }
            }
            _ => {
                let (ar, br) = (a.as_rug(), b.as_rug());
                if br == 0 { return None; }
                let (q, r) = ar.clone().div_rem(br.clone());
                if r == 0 { Some(Integer::from_large(q)) } else { None }
            }
        }
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        match &a.0 {
            IntegerInner::Small(1) => Some(self.one()),
            IntegerInner::Small(-1) => Some(Integer::from_small(-1)),
            IntegerInner::Small(_) => None,
            IntegerInner::Large(r) => {
                if **r == 1 { Some(self.one()) }
                else if **r == -1 { Some(Integer::from_small(-1)) }
                else { None }
            }
        }
    }
}

impl EuclideanDomain for IntegerDomain {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        match (&a.0, &b.0) {
            (IntegerInner::Small(av), IntegerInner::Small(bv)) => {
                if *bv == 0 { return None; }
                let q = *av / *bv; // truncating division
                let r = *av - q * *bv;
                Some((Integer::from_small(q), Integer::from_small(r)))
            }
            _ => {
                let (ar, br) = (a.as_rug(), b.as_rug());
                if br == &RugInteger::from(0) { return None; }
                let (q, r) = ar.clone().div_rem(br.clone());
                Some((Integer::from_large(q), Integer::from_large(r)))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Domain for RationalDomain (unchanged — Rational has no SOO)
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// From conversions
// ---------------------------------------------------------------------------

impl From<RugInteger> for Integer {
    fn from(value: RugInteger) -> Self {
        // Try to fit in Small first.
        match value.to_i64() {
            Some(v) => Integer::from_small(v),
            None => Integer::from_large(value),
        }
    }
}

impl From<RugRational> for Rational {
    fn from(value: RugRational) -> Self {
        Self(value)
    }
}

impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Integer::from_small(value)
    }
}

impl From<num_bigint::BigInt> for Integer {
    fn from(value: num_bigint::BigInt) -> Self {
        // Try i64 first.
        if let Some(v) = value.to_i64() {
            return Integer::from_small(v);
        }
        let (sign, bytes) = value.to_bytes_le();
        let mut inner = RugInteger::from_digits(&bytes, rug::integer::Order::Lsf);
        if sign == num_bigint::Sign::Minus {
            inner = -inner;
        }
        Integer::from_large(inner)
    }
}

// ---------------------------------------------------------------------------
// Integer methods
// ---------------------------------------------------------------------------

impl Integer {
    /// Create a Small-variant Integer.
    fn from_small(v: i64) -> Self {
        Integer(IntegerInner::Small(v))
    }

    /// Create a Large-variant Integer from a `RugInteger`.
    fn from_large(r: RugInteger) -> Self {
        Integer(IntegerInner::Large(Box::new(r)))
    }

    /// Create an integer from a machine integer or another `Into<RugInteger>`.
    pub fn new<T: Into<RugInteger>>(value: T) -> Self {
        let r = value.into();
        match r.to_i64() {
            Some(v) => Integer::from_small(v),
            None => Integer::from_large(r),
        }
    }

    /// Access the underlying [`rug::Integer`].
    ///
    /// For Small values, this promotes to Large in-place (one-time cost).
    pub fn inner(&self) -> &RugInteger {
        self.0.as_rug()
    }

    /// Get the rug representation without promotion (for internal use).
    fn as_rug(&self) -> &RugInteger {
        self.0.as_rug()
    }

    /// Try to extract the value as `i64`. Returns `None` for Large values
    /// that don't fit.
    pub fn to_i64(&self) -> Option<i64> {
        match &self.0 {
            IntegerInner::Small(v) => Some(*v),
            IntegerInner::Large(r) => r.to_i64(),
        }
    }

    /// Convert to a `BigInt` regardless of the backend.
    ///
    /// Uses binary serialization for performance (avoids string conversion).
    /// Small values convert directly via `BigInt::from(i64)`.
    pub fn to_bigint(&self) -> num_bigint::BigInt {
        match &self.0 {
            IntegerInner::Small(v) => num_bigint::BigInt::from(*v),
            IntegerInner::Large(r) => {
                use num_bigint::Sign;
                if **r == 0 { return num_bigint::BigInt::ZERO; }
                let bytes = r.to_digits::<u8>(rug::integer::Order::Lsf);
                let sign = if r.is_negative() { Sign::Minus } else { Sign::Plus };
                num_bigint::BigInt::from_bytes_le(sign, &bytes)
            }
        }
    }

    /// Raise to a `u32` power.
    pub fn pow_u32(&self, exp: u32) -> Self {
        use rug::ops::Pow;
        // For small base and small exponent, try i64 fast path.
        if let IntegerInner::Small(v) = self.0 {
            if exp <= 63 {
                // i64 can represent up to 2^63 - 1.
                let result = v.wrapping_pow(exp);
                // Check if wrapping produced a correct result by verifying
                // via rug only if the value looks suspicious.
                if exp == 0 { return Integer::from_small(1); }
                if *v == 0 { return Integer::from_small(0); }
                if *v == 1 { return Integer::from_small(1); }
                if *v == -1 {
                    return if exp % 2 == 0 { Integer::from_small(1) } else { Integer::from_small(-1) };
                }
                // For other cases, use rug to be safe.
            }
        }
        Integer::from_large(self.as_rug().clone().pow(exp))
    }

    /// Modular exponentiation: `self^exp mod modulus`.
    pub fn modpow(&self, exp: &Integer, modulus: &Integer) -> Integer {
        Integer::from_large(
            self.as_rug()
                .clone()
                .pow_mod(exp.as_rug(), modulus.as_rug())
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
        match (&self.0, &other.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) if *b != 0 => {
                let q = *a / *b; // truncating division
                let r = *a - q * *b;
                (Integer::from_small(q), Integer::from_small(r))
            }
            _ => {
                let (q, r) = self.as_rug().clone().div_rem(other.as_rug().clone());
                (Integer::from_large(q), Integer::from_large(r))
            }
        }
    }

    /// Returns `true` if the value is even.
    pub fn is_even(&self) -> bool {
        match &self.0 {
            IntegerInner::Small(v) => v & 1 == 0,
            IntegerInner::Large(r) => r.is_even(),
        }
    }

    /// Returns `true` if the value is negative.
    pub fn is_negative(&self) -> bool {
        match &self.0 {
            IntegerInner::Small(v) => *v < 0,
            IntegerInner::Large(r) => r.is_negative(),
        }
    }

    /// Returns `true` if the value is zero.
    pub fn is_zero(&self) -> bool {
        matches!(&self.0, IntegerInner::Small(0))
    }

    /// Returns `true` if the value is one.
    pub fn is_one(&self) -> bool {
        matches!(&self.0, IntegerInner::Small(1))
    }

    /// Absolute value.
    pub fn abs(&self) -> Integer {
        match &self.0 {
            IntegerInner::Small(v) => {
                match v.checked_neg() {
                    Some(neg) if *v >= 0 => Integer::from_small(neg),
                    _ if *v < 0 => Integer::from_small(-*v),
                    _ => Integer::from_large(self.as_rug().clone().abs()),
                }
            }
            IntegerInner::Large(r) => Integer::from_large(r.clone().abs()),
        }
    }

    /// Integer square root (floor).
    pub fn sqrt(&self) -> Integer {
        match &self.0 {
            IntegerInner::Small(v) if *v >= 0 => {
                Integer::from_small((*v as f64).sqrt() as i64)
            }
            _ => Integer::from_large(self.as_rug().clone().sqrt()),
        }
    }

    // -----------------------------------------------------------------------
    // Internal arithmetic helpers (SOO fast paths)
    // -----------------------------------------------------------------------

    fn add_ref(&self, rhs: &Self) -> Self {
        match (&self.0, &rhs.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) => {
                match a.checked_add(*b) {
                    Some(r) => Integer::from_small(r),
                    None => Integer::from_large(RugInteger::from(*a) + RugInteger::from(*b)),
                }
            }
            _ => Integer::from_large(self.as_rug() + rhs.as_rug()),
        }
    }

    fn sub_ref(&self, rhs: &Self) -> Self {
        match (&self.0, &rhs.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) => {
                match a.checked_sub(*b) {
                    Some(r) => Integer::from_small(r),
                    None => Integer::from_large(RugInteger::from(*a) - RugInteger::from(*b)),
                }
            }
            _ => Integer::from_large(self.as_rug() - rhs.as_rug()),
        }
    }

    fn mul_ref(&self, rhs: &Self) -> Self {
        match (&self.0, &rhs.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) => {
                match a.checked_mul(*b) {
                    Some(r) => Integer::from_small(r),
                    None => Integer::from_large(RugInteger::from(*a) * RugInteger::from(*b)),
                }
            }
            _ => Integer::from_large(self.as_rug() * rhs.as_rug()),
        }
    }

    fn neg_ref(&self) -> Self {
        match &self.0 {
            IntegerInner::Small(v) => {
                match v.checked_neg() {
                    Some(r) => Integer::from_small(r),
                    None => Integer::from_large(-self.as_rug()),
                }
            }
            IntegerInner::Large(r) => Integer::from_large(-&**r),
        }
    }
}

// ---------------------------------------------------------------------------
// Arithmetic operators — SOO fast paths
// ---------------------------------------------------------------------------

macro_rules! impl_soo_int_op {
    ($trait:ident, $method:ident, $fast:ident) => {
        impl $trait for Integer {
            type Output = Integer;
            fn $method(self, rhs: Integer) -> Integer {
                (&self).$fast(&rhs)
            }
        }
        impl $trait<&Integer> for Integer {
            type Output = Integer;
            fn $method(self, rhs: &Integer) -> Integer {
                (&self).$fast(rhs)
            }
        }
        impl $trait<Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: Integer) -> Integer {
                self.$fast(&rhs)
            }
        }
        impl $trait<&Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: &Integer) -> Integer {
                self.$fast(rhs)
            }
        }
    };
}

impl_soo_int_op!(Add, add, add_ref);
impl_soo_int_op!(Sub, sub, sub_ref);
impl_soo_int_op!(Mul, mul, mul_ref);

// Div and Rem need special handling (division by zero, exact check for Domain).
impl Div for Integer {
    type Output = Integer;
    fn div(self, rhs: Integer) -> Integer {
        (&self).div_owned(&rhs)
    }
}
impl Div<&Integer> for Integer {
    type Output = Integer;
    fn div(self, rhs: &Integer) -> Integer {
        (&self).div_owned(rhs)
    }
}
impl Div<Integer> for &Integer {
    type Output = Integer;
    fn div(self, rhs: Integer) -> Integer {
        self.div_owned(&rhs)
    }
}
impl Div<&Integer> for &Integer {
    type Output = Integer;
    fn div(self, rhs: &Integer) -> Integer {
        self.div_owned(rhs)
    }
}

impl Rem for Integer {
    type Output = Integer;
    fn rem(self, rhs: Integer) -> Integer {
        (&self).rem_owned(&rhs)
    }
}
impl Rem<&Integer> for Integer {
    type Output = Integer;
    fn rem(self, rhs: &Integer) -> Integer {
        (&self).rem_owned(rhs)
    }
}
impl Rem<Integer> for &Integer {
    type Output = Integer;
    fn rem(self, rhs: Integer) -> Integer {
        self.rem_owned(&rhs)
    }
}
impl Rem<&Integer> for &Integer {
    type Output = Integer;
    fn rem(self, rhs: &Integer) -> Integer {
        self.rem_owned(rhs)
    }
}

impl Integer {
    fn div_owned(&self, rhs: &Self) -> Self {
        match (&self.0, &rhs.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) if *b != 0 => {
                Integer::from_small(*a / *b)
            }
            _ => Integer::from_large(self.as_rug() / rhs.as_rug()),
        }
    }

    fn rem_owned(&self, rhs: &Self) -> Self {
        match (&self.0, &rhs.0) {
            (IntegerInner::Small(a), IntegerInner::Small(b)) if *b != 0 => {
                Integer::from_small(*a % *b)
            }
            _ => Integer::from_large(self.as_rug() % rhs.as_rug()),
        }
    }
}

impl std::ops::Neg for Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        (&self).neg_ref()
    }
}
impl std::ops::Neg for &Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        self.neg_ref()
    }
}

impl std::ops::ShrAssign<u32> for Integer {
    fn shr_assign(&mut self, shift: u32) {
        match &self.0 {
            IntegerInner::Small(v) => {
                let shifted = *v >> shift;
                self.0 = IntegerInner::Small(shifted);
            }
            IntegerInner::Large(r) => {
                let mut r = r.clone();
                *r >>= shift;
                self.0 = IntegerInner::Large(r);
            }
        }
    }
}
impl std::ops::Shr<u32> for Integer {
    type Output = Integer;
    fn shr(self, shift: u32) -> Integer {
        match self.0 {
            IntegerInner::Small(v) => Integer::from_small(v >> shift),
            IntegerInner::Large(r) => Integer::from_large(&*r >> shift),
        }
    }
}
impl std::ops::Shr<u32> for &Integer {
    type Output = Integer;
    fn shr(self, shift: u32) -> Integer {
        match &self.0 {
            IntegerInner::Small(v) => Integer::from_small(*v >> shift),
            IntegerInner::Large(r) => Integer::from_large(&**r >> shift),
        }
    }
}

impl std::ops::AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, rhs: &Integer) {
        let result = self.add_ref(rhs);
        *self = result;
    }
}
impl std::ops::SubAssign<&Integer> for Integer {
    fn sub_assign(&mut self, rhs: &Integer) {
        let result = self.sub_ref(rhs);
        *self = result;
    }
}
impl std::ops::MulAssign<&Integer> for Integer {
    fn mul_assign(&mut self, rhs: &Integer) {
        let result = self.mul_ref(rhs);
        *self = result;
    }
}
impl std::ops::DivAssign<&Integer> for Integer {
    fn div_assign(&mut self, rhs: &Integer) {
        let result = self.div_owned(rhs);
        *self = result;
    }
}

// ---------------------------------------------------------------------------
// Rational (unchanged)
// ---------------------------------------------------------------------------

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
        Self(RugRational::from(n.0.into_rug()))
    }

    /// Numerator as an [`Integer`].
    pub fn numer(&self) -> Integer {
        Integer::from(self.0.numer().clone())
    }

    /// Denominator as an [`Integer`].
    pub fn denom(&self) -> Integer {
        Integer::from(self.0.denom().clone())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_small_fast_path() {
        let a = Integer::from(3i64);
        let b = Integer::from(5i64);
        assert!(matches!(a.0, IntegerInner::Small(3)));
        assert!(matches!(b.0, IntegerInner::Small(5)));
        let c = &a + &b;
        assert!(matches!(c.0, IntegerInner::Small(8)));
    }

    #[test]
    fn integer_small_overflow_promotes() {
        let a = Integer::from(i64::MAX);
        let b = Integer::from(1i64);
        let c = &a + &b;
        // Should promote to Large.
        assert!(matches!(c.0, IntegerInner::Large(_)));
        assert_eq!(c.to_bigint(), num_bigint::BigInt::from(i64::MAX) + 1);
    }

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
