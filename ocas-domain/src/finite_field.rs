//! Finite field domain implementation.
//!
//! Supports prime fields $\mathbb{Z}/p\mathbb{Z}$ where $p$ is prime.
//! Arithmetic is performed with arbitrary-precision integers so that large
//! primes can be used, though the modulus must fit in memory.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::domain::Domain;

/// An element of a prime finite field.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use ocas_domain::{Domain, FiniteField};
///
/// let f = FiniteField::new(BigInt::from(7));
/// let a = f.element(10);
/// assert_eq!(a.value().to_string(), "3");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FiniteFieldElement {
    value: BigInt,
}

impl FiniteFieldElement {
    /// Access the canonical representative in `[0, p-1]`.
    pub fn value(&self) -> &BigInt {
        &self.value
    }
}

/// A prime finite field $\mathbb{Z}/p\mathbb{Z}$.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use ocas_domain::{Domain, FiniteField};
///
/// let f = FiniteField::new(BigInt::from(7));
/// let a = f.element(3);
/// let b = f.element(5);
/// assert_eq!(f.add(&a, &b), f.element(1));
/// assert_eq!(f.mul(&a, &b), f.element(1));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteField {
    prime: BigInt,
}

impl FiniteField {
    /// Create a finite field with the given prime modulus.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `prime` is less than 2.
    pub fn new(prime: BigInt) -> Self {
        debug_assert!(prime > BigInt::one(), "modulus must be at least 2");
        Self { prime }
    }

    /// Create a field element from an arbitrary integer.
    ///
    /// The value is reduced into the canonical range `[0, p-1]`.
    pub fn element(&self, value: impl Into<BigInt>) -> FiniteFieldElement {
        let value = value.into();
        let value = value.mod_floor(&self.prime);
        FiniteFieldElement { value }
    }

    /// Return the field modulus.
    pub fn prime(&self) -> &BigInt {
        &self.prime
    }

    fn normalize(&self, value: BigInt) -> FiniteFieldElement {
        FiniteFieldElement {
            value: value.mod_floor(&self.prime),
        }
    }
}

impl Domain for FiniteField {
    type Element = FiniteFieldElement;

    fn zero(&self) -> Self::Element {
        self.element(BigInt::zero())
    }

    fn one(&self) -> Self::Element {
        self.element(BigInt::one())
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        self.normalize(a.value.clone() + b.value.clone())
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        self.normalize(a.value.clone() - b.value.clone())
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        self.normalize(-a.value.clone())
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        self.normalize(a.value.clone() * b.value.clone())
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        self.inv(b).map(|inv| self.mul(a, &inv))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.value.is_zero() {
            return None;
        }
        // Fermat's little theorem: a^(p-2) ≡ a^(-1) (mod p) for prime p.
        let exp = &self.prime - BigInt::one() - BigInt::one();
        Some(self.normalize(a.value.modpow(&exp, &self.prime)))
    }

    fn is_zero(&self, a: &Self::Element) -> bool {
        a.value.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_field_addition() {
        let f = FiniteField::new(BigInt::from(7));
        let a = f.element(3);
        let b = f.element(5);
        assert_eq!(f.add(&a, &b), f.element(1));
    }

    #[test]
    fn finite_field_subtraction_wraps() {
        let f = FiniteField::new(BigInt::from(7));
        let a = f.element(3);
        let b = f.element(5);
        assert_eq!(f.sub(&a, &b), f.element(5));
    }

    #[test]
    fn finite_field_multiplication() {
        let f = FiniteField::new(BigInt::from(7));
        let a = f.element(3);
        let b = f.element(5);
        assert_eq!(f.mul(&a, &b), f.element(1));
    }

    #[test]
    fn finite_field_inverse() {
        let f = FiniteField::new(BigInt::from(7));
        let a = f.element(3);
        let inv = f.inv(&a).unwrap();
        assert_eq!(f.mul(&a, &inv), f.one());
    }

    #[test]
    fn finite_field_inverse_of_zero_is_none() {
        let f = FiniteField::new(BigInt::from(7));
        assert!(f.inv(&f.zero()).is_none());
    }

    #[test]
    fn finite_field_division() {
        let f = FiniteField::new(BigInt::from(7));
        let a = f.element(5);
        let b = f.element(3);
        let q = f.div(&a, &b).unwrap();
        assert_eq!(f.mul(&q, &b), a);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn small_field() -> FiniteField {
        FiniteField::new(BigInt::from(7))
    }

    proptest! {
        #[test]
        fn multiplication_is_commutative(a in 0u64..7, b in 0u64..7) {
            let f = small_field();
            let x = f.element(a);
            let y = f.element(b);
            assert_eq!(f.mul(&x, &y), f.mul(&y, &x));
        }

        #[test]
        fn add_then_sub_is_identity(a in 0u64..7, b in 0u64..7) {
            let f = small_field();
            let x = f.element(a);
            let y = f.element(b);
            let sum = f.add(&x, &y);
            assert_eq!(f.sub(&sum, &y), x);
        }

        #[test]
        fn non_zero_inverse_exists(a in 1u64..7) {
            let f = small_field();
            let x = f.element(a);
            let inv = f.inv(&x).unwrap();
            assert_eq!(f.mul(&x, &inv), f.one());
        }
    }
}
