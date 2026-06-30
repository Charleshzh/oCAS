//! Complex number domain implementation.
//!
//! This module provides a generic complex domain built on top of any base
//! [`Domain`]. The default build uses [`num_complex::Complex`] as the storage
//! representation; a future MPC backend may be added behind a feature flag.

use std::marker::PhantomData;

use num_complex::Complex as NumComplex;

use crate::domain::Domain;

/// A complex number whose real and imaginary parts live in a base domain.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Complex<D: Domain> {
    inner: NumComplex<D::Element>,
}

impl<D: Domain> Complex<D> {
    /// Create a complex number from real and imaginary parts.
    pub fn new(real: D::Element, imag: D::Element) -> Self {
        Self {
            inner: NumComplex::new(real, imag),
        }
    }

    /// Access the real part.
    pub fn re(&self) -> &D::Element {
        &self.inner.re
    }

    /// Access the imaginary part.
    pub fn im(&self) -> &D::Element {
        &self.inner.im
    }

    /// Access the underlying `num_complex::Complex`.
    pub fn inner(&self) -> &NumComplex<D::Element> {
        &self.inner
    }
}

/// The complex domain over a base domain `D`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComplexDomain<D: Domain> {
    base: D,
    _marker: PhantomData<D>,
}

impl<D: Domain> ComplexDomain<D> {
    /// Create the complex domain over `base`.
    pub fn new(base: D) -> Self {
        Self {
            base,
            _marker: PhantomData,
        }
    }

    /// Return a reference to the base domain.
    pub fn base(&self) -> &D {
        &self.base
    }

    /// Create a purely real element in this domain.
    pub fn real_element(&self, re: D::Element) -> Complex<D> {
        Complex::new(re, self.base.zero())
    }

    /// Create a purely imaginary element in this domain.
    pub fn imag_element(&self, im: D::Element) -> Complex<D> {
        Complex::new(self.base.zero(), im)
    }
}

impl<D: Domain> Domain for ComplexDomain<D> {
    type Element = Complex<D>;

    fn zero(&self) -> Self::Element {
        Complex::new(self.base.zero(), self.base.zero())
    }

    fn one(&self) -> Self::Element {
        Complex::new(self.base.one(), self.base.zero())
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Complex::new(self.base.add(a.re(), b.re()), self.base.add(a.im(), b.im()))
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        Complex::new(self.base.sub(a.re(), b.re()), self.base.sub(a.im(), b.im()))
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        Complex::new(self.base.neg(a.re()), self.base.neg(a.im()))
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        let ac = self.base.mul(a.re(), b.re());
        let bd = self.base.mul(a.im(), b.im());
        let ad = self.base.mul(a.re(), b.im());
        let bc = self.base.mul(a.im(), b.re());
        Complex::new(self.base.sub(&ac, &bd), self.base.add(&ad, &bc))
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        // (a + bi)/(c + di) = ((ac + bd) + (bc - ad)i) / (c^2 + d^2)
        let c2 = self.base.mul(b.re(), b.re());
        let d2 = self.base.mul(b.im(), b.im());
        let denom = self.base.add(&c2, &d2);
        let ac = self.base.mul(a.re(), b.re());
        let bd = self.base.mul(a.im(), b.im());
        let bc = self.base.mul(a.im(), b.re());
        let ad = self.base.mul(a.re(), b.im());
        let real_num = self.base.add(&ac, &bd);
        let imag_num = self.base.sub(&bc, &ad);
        let real = self.base.div(&real_num, &denom)?;
        let imag = self.base.div(&imag_num, &denom)?;
        Some(Complex::new(real, imag))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        self.div(&self.one(), a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{IntegerDomain, Rational, RationalDomain};

    #[test]
    fn complex_addition() {
        let base = IntegerDomain;
        let domain = ComplexDomain::new(base);
        let a = Complex::new(1.into(), 2.into());
        let b = Complex::new(3.into(), 4.into());
        let sum = domain.add(&a, &b);
        assert_eq!(sum.re(), &4.into());
        assert_eq!(sum.im(), &6.into());
    }

    #[test]
    fn complex_multiplication() {
        let base = IntegerDomain;
        let domain = ComplexDomain::new(base);
        let a = Complex::new(1.into(), 2.into());
        let b = Complex::new(3.into(), 4.into());
        let prod = domain.mul(&a, &b);
        // (1 + 2i)(3 + 4i) = -5 + 10i
        assert_eq!(prod.re(), &(-5).into());
        assert_eq!(prod.im(), &10.into());
    }

    #[test]
    fn complex_i_squared_is_minus_one() {
        let base = IntegerDomain;
        let domain = ComplexDomain::new(base);
        let i = domain.imag_element(1.into());
        let i2 = domain.mul(&i, &i);
        assert_eq!(i2.re(), &(-1).into());
        assert_eq!(i2.im(), &0.into());
    }

    #[test]
    fn complex_inverse_over_rationals() {
        let base = RationalDomain;
        let domain = ComplexDomain::new(base);
        let z = Complex::new(Rational::new(1, 1), Rational::new(1, 1));
        let inv = domain.inv(&z).expect("non-zero complex is invertible");
        let prod = domain.mul(&z, &inv);
        assert_eq!(prod, domain.one());
    }
}
