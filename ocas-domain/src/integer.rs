//! Integer domain implementation.

use num_bigint::BigInt;

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
