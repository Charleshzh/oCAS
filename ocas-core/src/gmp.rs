//! Optional GMP backend integration for oCAS.
//!
//! This module is only available when the `gmp` feature is enabled. It provides
//! a thin wrapper around the `rug` crate's arbitrary-precision integer type,
//! exposing a small, stable surface that higher-level oCAS crates can use as a
//! backend.
//!
//! In the default configuration `ocas-core` does not link against GMP, keeping
//! the workspace buildable on platforms where GMP is unavailable (for example,
//! MSVC on Windows).

#[cfg(feature = "gmp")]
use std::fmt;

#[cfg(feature = "gmp")]
use gmp::Integer as RugInteger;

#[cfg(feature = "gmp")]
use crate::error::Result;

/// An arbitrary-precision integer backed by GMP.
///
/// Only available when the `gmp` feature is enabled. The wrapper intentionally
/// hides most of `rug`'s API, exposing only the operations needed by oCAS's
/// numerical backends. This keeps the dependency on `rug` replaceable in the
/// future if necessary.
///
/// # Example
///
/// ```
/// use ocas_core::gmp::GmpInteger;
///
/// let a = GmpInteger::from_i64(21);
/// let b = GmpInteger::from_i64(21);
/// let sum = a.add(&b);
/// assert_eq!(sum.to_decimal_string(), "42");
/// ```
#[cfg(feature = "gmp")]
#[derive(Debug, Clone)]
pub struct GmpInteger {
    inner: RugInteger,
}

#[cfg(feature = "gmp")]
impl GmpInteger {
    /// Construct a `GmpInteger` from an `i64`.
    pub fn from_i64(value: i64) -> Self {
        Self {
            inner: RugInteger::from(value),
        }
    }

    /// Construct a `GmpInteger` from a `num-bigint::BigInt`.
    ///
    /// # Errors
    ///
    /// Returns `OcasError::BackendError` if the conversion fails.
    pub fn from_bigint(value: &num_bigint::BigInt) -> Result<Self> {
        let (sign, bytes) = value.to_bytes_le();
        let mut inner = RugInteger::from_digits(&bytes, gmp::integer::Order::Lsf);
        if sign == num_bigint::Sign::Minus {
            inner = -inner;
        }
        Ok(Self { inner })
    }

    /// Return the underlying `rug::Integer` reference.
    pub fn inner(&self) -> &RugInteger {
        &self.inner
    }

    /// Add two `GmpInteger`s, returning a new value.
    pub fn add(&self, other: &Self) -> Self {
        Self {
            inner: (&self.inner + &other.inner).into(),
        }
    }

    /// Subtract `other` from `self`, returning a new value.
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            inner: (&self.inner - &other.inner).into(),
        }
    }

    /// Multiply two `GmpInteger`s, returning a new value.
    pub fn mul(&self, other: &Self) -> Self {
        Self {
            inner: (&self.inner * &other.inner).into(),
        }
    }

    /// Compare two `GmpInteger`s.
    pub fn compare(&self, other: &Self) -> std::cmp::Ordering {
        self.inner.cmp(&other.inner)
    }

    /// Return the decimal string representation.
    pub fn to_decimal_string(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(feature = "gmp")]
impl fmt::Display for GmpInteger {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.inner)
    }
}

#[cfg(feature = "gmp")]
impl PartialEq for GmpInteger {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

#[cfg(feature = "gmp")]
impl Eq for GmpInteger {}

#[cfg(feature = "gmp")]
impl PartialOrd for GmpInteger {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(std::cmp::Ord::cmp(self, other))
    }
}

#[cfg(feature = "gmp")]
impl Ord for GmpInteger {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.compare(other)
    }
}

#[cfg(all(test, feature = "gmp"))]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use proptest::prelude::*;

    #[test]
    fn gmp_construction_from_i64() {
        let zero = GmpInteger::from_i64(0);
        let max = GmpInteger::from_i64(i64::MAX);
        let min = GmpInteger::from_i64(i64::MIN);

        assert_eq!(zero.to_decimal_string(), "0");
        assert_eq!(max.to_decimal_string(), i64::MAX.to_string());
        assert_eq!(min.to_decimal_string(), i64::MIN.to_string());
    }

    #[test]
    fn gmp_arithmetic_basic() {
        let a = GmpInteger::from_i64(21);
        let b = GmpInteger::from_i64(21);
        let sum = a.add(&b);
        let diff = sum.sub(&a);
        let prod = a.mul(&b);

        assert_eq!(sum.to_decimal_string(), "42");
        assert_eq!(diff.to_decimal_string(), "21");
        assert_eq!(prod.to_decimal_string(), "441");
    }

    #[test]
    fn gmp_comparison() {
        let a = GmpInteger::from_i64(-5);
        let b = GmpInteger::from_i64(5);
        assert!(a.compare(&b) == std::cmp::Ordering::Less);
        assert!(b.compare(&a) == std::cmp::Ordering::Greater);
        assert_eq!(
            a.compare(&GmpInteger::from_i64(-5)),
            std::cmp::Ordering::Equal
        );
    }

    #[test]
    fn gmp_from_bigint_roundtrip() {
        let big = BigInt::from(1234567890123456789i64) * BigInt::from(9876543210987654321u64);
        let gmp = GmpInteger::from_bigint(&big).expect("conversion should succeed");
        assert_eq!(gmp.to_decimal_string(), big.to_string());
    }

    #[test]
    fn gmp_from_bigint_negative() {
        let big = BigInt::from(-1_000_000_000_000i64);
        let gmp = GmpInteger::from_bigint(&big).expect("conversion should succeed");
        assert_eq!(gmp.to_decimal_string(), "-1000000000000");
    }

    proptest! {
        #[test]
        fn gmp_arithmetic_matches_num_bigint(a in any::<i64>(), b in any::<i64>()) {
            let ga = GmpInteger::from_i64(a);
            let gb = GmpInteger::from_i64(b);

            let expected_add = BigInt::from(a) + BigInt::from(b);
            let expected_sub = BigInt::from(a) - BigInt::from(b);
            let expected_mul = BigInt::from(a) * BigInt::from(b);

            let actual_add = ga.add(&gb).to_decimal_string().parse::<BigInt>().unwrap();
            let actual_sub = ga.sub(&gb).to_decimal_string().parse::<BigInt>().unwrap();
            let actual_mul = ga.mul(&gb).to_decimal_string().parse::<BigInt>().unwrap();

            prop_assert_eq!(actual_add, expected_add);
            prop_assert_eq!(actual_sub, expected_sub);
            prop_assert_eq!(actual_mul, expected_mul);
        }
    }
}
