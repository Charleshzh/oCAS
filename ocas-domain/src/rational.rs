//! Rational domain implementation.

use num_rational::BigRational;

/// Arbitrary-precision rational number.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Rational(BigRational);

impl Rational {
    /// Create a rational number from a numerator and denominator.
    pub fn new(numer: i64, denom: i64) -> Self {
        Self(BigRational::from_integer(numer.into()) / BigRational::from_integer(denom.into()))
    }

    /// Access the underlying `BigRational`.
    pub fn inner(&self) -> &BigRational {
        &self.0
    }
}
