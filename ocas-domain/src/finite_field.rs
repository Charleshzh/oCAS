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

impl std::fmt::Display for FiniteFieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
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
    /// Cached `prime - 2` for fast modular inversion via Fermat's little theorem.
    prime_minus_two: BigInt,
    /// GMP-backed prime for fast modular arithmetic (only with `gmp` feature).
    #[cfg(feature = "gmp")]
    prime_gmp: rug::Integer,
    /// GMP-backed `prime - 2` for fast modular inversion (only with `gmp` feature).
    #[cfg(feature = "gmp")]
    prime_minus_two_gmp: rug::Integer,
}

impl FiniteField {
    /// Create a finite field with the given prime modulus.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `prime` is less than 2.
    pub fn new(prime: BigInt) -> Self {
        debug_assert!(prime > BigInt::one(), "modulus must be at least 2");
        let prime_minus_two = &prime - BigInt::one() - BigInt::one();
        #[cfg(feature = "gmp")]
        {
            let (sign, bytes) = prime.to_bytes_le();
            let mut prime_gmp = rug::Integer::from_digits(&bytes, rug::integer::Order::Lsf);
            if sign == num_bigint::Sign::Minus {
                prime_gmp = -prime_gmp;
            }
            let (sign2, bytes2) = prime_minus_two.to_bytes_le();
            let mut pmt_gmp = rug::Integer::from_digits(&bytes2, rug::integer::Order::Lsf);
            if sign2 == num_bigint::Sign::Minus {
                pmt_gmp = -pmt_gmp;
            }
            Self {
                prime,
                prime_minus_two,
                prime_gmp,
                prime_minus_two_gmp: pmt_gmp,
            }
        }
        #[cfg(not(feature = "gmp"))]
        Self {
            prime,
            prime_minus_two,
        }
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

    /// Convert a `BigInt` to a `rug::Integer` using binary serialization.
    #[cfg(feature = "gmp")]
    fn bigint_to_rug(value: &BigInt) -> rug::Integer {
        let (sign, bytes) = value.to_bytes_le();
        let mut inner = rug::Integer::from_digits(&bytes, rug::integer::Order::Lsf);
        if sign == num_bigint::Sign::Minus {
            inner = -inner;
        }
        inner
    }

    /// Convert a `rug::Integer` back to `BigInt` using binary serialization.
    #[cfg(feature = "gmp")]
    fn rug_to_bigint(value: &rug::Integer) -> BigInt {
        use num_bigint::Sign;
        if *value == 0 {
            return BigInt::ZERO;
        }
        let num_bytes = value.significant_digits::<u8>();
        let mut bytes = vec![0u8; num_bytes];
        value.write_digits(&mut bytes, rug::integer::Order::Lsf);
        let sign = if value.is_negative() {
            Sign::Minus
        } else {
            Sign::Plus
        };
        BigInt::from_bytes_le(sign, &bytes)
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
        #[cfg(feature = "gmp")]
        {
            let base = Self::bigint_to_rug(&a.value);
            let result = base
                .pow_mod(&self.prime_minus_two_gmp, &self.prime_gmp)
                .expect("modpow: modulus must be positive");
            Some(FiniteFieldElement {
                value: Self::rug_to_bigint(&result),
            })
        }
        #[cfg(not(feature = "gmp"))]
        {
            Some(self.normalize(a.value.modpow(&self.prime_minus_two, &self.prime)))
        }
    }

    fn is_zero(&self, a: &Self::Element) -> bool {
        a.value.is_zero()
    }

    /// Modular exponentiation using `modpow`, much faster than the
    /// default binary exponentiation for large exponents.
    fn pow(&self, a: &Self::Element, n: u64) -> Self::Element {
        if n == 0 {
            return self.one();
        }
        if self.is_zero(a) {
            return self.zero();
        }
        #[cfg(feature = "gmp")]
        {
            let base = Self::bigint_to_rug(&a.value);
            let exp = rug::Integer::from(n);
            let modulus = &self.prime_gmp;
            let result = base
                .pow_mod(&exp, modulus)
                .expect("modpow: modulus must be positive");
            FiniteFieldElement {
                value: Self::rug_to_bigint(&result),
            }
        }
        #[cfg(not(feature = "gmp"))]
        {
            let exp = BigInt::from(n);
            self.normalize(a.value.modpow(&exp, &self.prime))
        }
    }
}

impl crate::domain::EuclideanDomain for FiniteField {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        // Every nonzero element of a field is a unit, so division is exact
        // and the remainder is always zero.
        self.div(a, b).map(|q| (q, self.zero()))
    }

    fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        // In a field the GCD is degenerate: 0 if both are zero, else 1.
        if self.is_zero(a) && self.is_zero(b) {
            self.zero()
        } else {
            self.one()
        }
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
