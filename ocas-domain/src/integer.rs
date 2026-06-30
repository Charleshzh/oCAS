//! Integer domain implementation.

use num_bigint::BigInt;
use num_integer::Integer as _;
use num_traits::{One, Zero};

use crate::domain::{Domain, EuclideanDomain};

/// Arbitrary-precision integer.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Integer(BigInt);

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

/// The integer domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntegerDomain;

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

impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Self::new(value)
    }
}

impl From<BigInt> for Integer {
    fn from(value: BigInt) -> Self {
        Self(value)
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
}
