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
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
