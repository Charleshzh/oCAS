//! Integer domain implementation.
//!
//! Supports arbitrary-precision integers. The default build uses
//! [`num_bigint::BigInt`]. When the `gmp` feature is enabled, the
//! implementation moves to `gmp_backend`.

#[cfg(not(feature = "gmp"))]
use num_bigint::BigInt;
#[cfg(not(feature = "gmp"))]
use num_integer::Integer as _;
#[cfg(not(feature = "gmp"))]
use num_traits::{One, Zero};
#[cfg(not(feature = "gmp"))]
use std::ops::{Add, Div, Mul, Rem, Sub};

#[cfg(not(feature = "gmp"))]
use crate::domain::{Domain, EuclideanDomain};

/// The integer domain.
///
/// # Example
///
/// ```
/// use ocas_domain::{Domain, Integer, IntegerDomain};
///
/// let domain = IntegerDomain;
/// let a = Integer::from(3);
/// let b = Integer::from(5);
/// assert_eq!(domain.add(&a, &b), Integer::from(8));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntegerDomain;

#[cfg(not(feature = "gmp"))]
/// Arbitrary-precision integer backed by `num-bigint`.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
///
/// let a = Integer::from(42);
/// let b = Integer::new(100);
/// assert_eq!(a.inner().to_string(), "42");
/// assert_eq!(b.inner().to_string(), "100");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Integer(BigInt);

#[cfg(not(feature = "gmp"))]
impl Integer {
    /// Create an integer from a machine integer.
    pub fn new<T: Into<BigInt>>(value: T) -> Self {
        Self(value.into())
    }

    /// Access the underlying `BigInt`.
    pub fn inner(&self) -> &BigInt {
        &self.0
    }

    /// Convert to a `BigInt` regardless of the backend.
    pub fn to_bigint(&self) -> BigInt {
        self.0.clone()
    }

    /// Raise to a `u32` power.
    pub fn pow_u32(&self, exp: u32) -> Self {
        use num_traits::Pow;
        Integer(self.0.clone().pow(exp))
    }
}

#[cfg(not(feature = "gmp"))]
impl std::fmt::Display for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(not(feature = "gmp"))]
impl Domain for IntegerDomain {
    type Element = Integer;

    fn zero(&self) -> Self::Element {
        Integer(BigInt::zero())
    }

    fn one(&self) -> Self::Element {
        Integer(BigInt::one())
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
        if b.0.is_zero() {
            return None;
        }
        let (q, r) = a.0.clone().div_rem(&b.0);
        if r.is_zero() { Some(Integer(q)) } else { None }
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.0.is_one() {
            Some(self.one())
        } else if a.0 == -BigInt::one() {
            Some(Integer(-BigInt::one()))
        } else {
            None
        }
    }
}

#[cfg(not(feature = "gmp"))]
impl EuclideanDomain for IntegerDomain {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        if b.0.is_zero() {
            return None;
        }
        let (q, r) = a.0.clone().div_rem(&b.0);
        Some((Integer(q), Integer(r)))
    }
}

#[cfg(not(feature = "gmp"))]
impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Self::new(value)
    }
}

#[cfg(not(feature = "gmp"))]
impl From<BigInt> for Integer {
    fn from(value: BigInt) -> Self {
        Self(value)
    }
}

// ---------------------------------------------------------------------------
// Arithmetic operators — forward to the inner BigInt.
// These allow `number_theory.rs` and factorization code to use `Integer`
// directly instead of reaching through `.inner()`.
// ---------------------------------------------------------------------------
macro_rules! impl_int_op_owned {
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
                Integer(self.0 $op &rhs.0)
            }
        }
        impl $trait<Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: Integer) -> Integer {
                Integer(&self.0 $op rhs.0)
            }
        }
        impl $trait<&Integer> for &Integer {
            type Output = Integer;
            fn $method(self, rhs: &Integer) -> Integer {
                Integer(&self.0 $op &rhs.0)
            }
        }
    };
}

impl_int_op_owned!(Add, add, +);
impl_int_op_owned!(Sub, sub, -);
impl_int_op_owned!(Mul, mul, *);
impl_int_op_owned!(Div, div, /);
impl_int_op_owned!(Rem, rem, %);

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
        Integer(&self.0 >> shift)
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

// ---------------------------------------------------------------------------
// Number-theory helper methods on Integer.
// These replace the old `BigInt`-based APIs and enable `number_theory.rs`
// to work with `Integer` directly.
// ---------------------------------------------------------------------------

#[cfg(not(feature = "gmp"))]
impl Integer {
    /// Modular exponentiation: `self^exp mod modulus`.
    pub fn modpow(&self, exp: &Integer, modulus: &Integer) -> Integer {
        Integer(self.0.modpow(&exp.0, &modulus.0))
    }

    /// Floor modulo: result `r` satisfies `0 ≤ r < |modulus|`.
    pub fn mod_floor(&self, modulus: &Integer) -> Integer {
        use num_integer::Integer as _;
        Integer(self.0.mod_floor(&modulus.0))
    }

    /// Division with remainder: `(quotient, remainder)`.
    pub fn div_rem(&self, other: &Integer) -> (Integer, Integer) {
        use num_integer::Integer as _;
        let (q, r) = self.0.div_rem(&other.0);
        (Integer(q), Integer(r))
    }

    /// Returns `true` if the value is even.
    pub fn is_even(&self) -> bool {
        use num_integer::Integer as _;
        self.0.is_even()
    }

    /// Returns `true` if the value is negative.
    pub fn is_negative(&self) -> bool {
        use num_traits::Signed;
        self.0.is_negative()
    }

    /// Returns `true` if the value is zero.
    pub fn is_zero(&self) -> bool {
        use num_traits::Zero;
        self.0.is_zero()
    }

    /// Returns `true` if the value is one.
    pub fn is_one(&self) -> bool {
        use num_traits::One;
        self.0.is_one()
    }

    /// Absolute value.
    pub fn abs(&self) -> Integer {
        use num_traits::Signed;
        Integer(self.0.abs())
    }

    /// Integer square root (floor).
    pub fn sqrt(&self) -> Integer {
        Integer(self.0.sqrt())
    }
}

#[cfg(not(feature = "gmp"))]
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
}

#[cfg(not(feature = "gmp"))]
#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn any_integer() -> impl Strategy<Value = Integer> {
        any::<i64>().prop_map(Integer::from)
    }

    proptest! {
        #[test]
        fn addition_is_commutative(a in any_integer(), b in any_integer()) {
            let domain = IntegerDomain;
            assert_eq!(domain.add(&a, &b), domain.add(&b, &a));
        }

        #[test]
        fn multiplication_is_commutative(a in any_integer(), b in any_integer()) {
            let domain = IntegerDomain;
            assert_eq!(domain.mul(&a, &b), domain.mul(&b, &a));
        }

        #[test]
        fn zero_is_additive_identity(a in any_integer()) {
            let domain = IntegerDomain;
            let zero = domain.zero();
            assert_eq!(domain.add(&a, &zero), a);
        }

        #[test]
        fn subtraction_cancels_addition(a in any_integer(), b in any_integer()) {
            let domain = IntegerDomain;
            let sum = domain.add(&a, &b);
            assert_eq!(domain.sub(&sum, &b), a);
        }
    }
}
