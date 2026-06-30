//! Core domain trait for generic algebraic algorithms.
//!
//! A [`Domain`] describes a set of values together with the basic arithmetic
//! operations needed by polynomial and matrix algorithms. Implementations are
//! provided for integers, rationals, and finite fields.

/// A coefficient domain for generic computer-algebra routines.
///
/// Domains describe operations on their elements. The domain object itself may
/// carry parameters (such as a finite-field modulus), so every operation
/// takes `&self`. This mirrors the conventional "domain object" pattern used
/// by Flint, SymPy's `Domain`, and other CAS libraries.
pub trait Domain: Clone + PartialEq + Eq + std::fmt::Debug + Sized {
    /// The type of elements in the domain.
    type Element: Clone + PartialEq + Eq + std::fmt::Debug;

    /// The additive identity.
    fn zero(&self) -> Self::Element;

    /// The multiplicative identity.
    fn one(&self) -> Self::Element;

    /// Add two elements.
    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;

    /// Subtract `b` from `a`.
    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;

    /// Negate an element.
    fn neg(&self, a: &Self::Element) -> Self::Element;

    /// Multiply two elements.
    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;

    /// Divide `a` by `b`.
    ///
    /// Returns `None` if division is not exact or `b` is zero.
    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element>;

    /// Return the multiplicative inverse of `a`.
    ///
    /// Returns `None` if `a` is zero.
    fn inv(&self, a: &Self::Element) -> Option<Self::Element>;

    /// Test whether an element is the additive identity.
    fn is_zero(&self, a: &Self::Element) -> bool {
        *a == self.zero()
    }

    /// Test whether an element is the multiplicative identity.
    fn is_one(&self, a: &Self::Element) -> bool {
        *a == self.one()
    }

    /// Return `a` raised to the non-negative integer power `n`.
    ///
    /// The default implementation uses binary exponentiation. Domains that
    /// can do better (e.g. modular exponentiation) may override it.
    fn pow(&self, a: &Self::Element, n: u64) -> Self::Element {
        let mut base = a.clone();
        let mut result = self.one();
        let mut exp = n;
        while exp > 0 {
            if exp & 1 == 1 {
                result = self.mul(&result, &base);
            }
            base = self.mul(&base, &base);
            exp >>= 1;
        }
        result
    }

    /// Convert a `u64` into an element of the domain.
    ///
    /// This is used by generic algorithms (e.g. polynomial differentiation)
    /// that need small positive integer coefficients. Domains that cannot
    /// represent every `u64` may wrap or truncate as appropriate for their
    /// semantics.
    fn cast_u64(&self, n: u64) -> Self::Element {
        let mut result = self.zero();
        let one = self.one();
        for _ in 0..n {
            result = self.add(&result, &one);
        }
        result
    }
}

/// Marker trait for domains that support exact division with remainder.
///
/// Euclidean domains provide `div_rem`, which returns the quotient and
/// remainder. The remainder must satisfy `rem == 0` or `deg(rem) < deg(b)`.
pub trait EuclideanDomain: Domain {
    /// Divide `a` by `b` returning `(quotient, remainder)`.
    ///
    /// Returns `None` if `b` is zero.
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)>;
}
