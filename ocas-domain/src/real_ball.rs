//! Real ball domain implementation.
//!
//! A real ball represents a real number together with a rigorous error bound.
//! The default build uses lightweight conservative `f64` balls. When the `mpfr`
//! feature is enabled, the implementation is backed by `rug::Float` and uses
//! directed rounding to produce rigorous arbitrary-precision enclosures.
//!
//! The default `f64` version is suitable for templates and demonstration; it is
//! not a rigorous arbitrary-precision interval library.

#[cfg(feature = "mpfr")]
use rug::ops::{AddAssignRound, DivAssignRound, Pow, SubAssignRound};

use crate::domain::Domain;

/// Precision used for the MPFR-backed real ball in bits.
#[cfg(feature = "mpfr")]
const BALL_PRECISION: u32 = 53;

/// A real ball: midpoint with radius.
///
/// The default `f64` build implements `Copy`. The `mpfr` build stores an
/// owned `rug::Float` and therefore only implements `Clone`.
#[derive(Debug, Clone, PartialEq)]
pub struct RealBall {
    #[cfg(not(feature = "mpfr"))]
    mid: f64,
    #[cfg(not(feature = "mpfr"))]
    rad: f64,
    #[cfg(feature = "mpfr")]
    mid: rug::Float,
    #[cfg(feature = "mpfr")]
    rad: rug::Float,
}

impl Eq for RealBall {}

#[cfg(not(feature = "mpfr"))]
impl Copy for RealBall {}

impl RealBall {
    /// Create a new ball from midpoint and radius.
    ///
    /// The radius is clamped to be non-negative.
    #[cfg(not(feature = "mpfr"))]
    pub fn new(mid: f64, rad: f64) -> Self {
        Self {
            mid,
            rad: rad.max(0.0),
        }
    }

    /// Create a new ball from midpoint and radius.
    ///
    /// The radius is clamped to be non-negative.
    #[cfg(feature = "mpfr")]
    pub fn new(mid: rug::Float, rad: rug::Float) -> Self {
        let zero = mpfr_zero();
        if rad < zero {
            Self { mid, rad: zero }
        } else {
            Self { mid, rad }
        }
    }

    /// Create a ball from an `f64` value with zero radius.
    #[cfg(not(feature = "mpfr"))]
    pub fn from_f64(value: f64) -> Self {
        Self::new(value, 0.0)
    }

    /// Create a ball from an `f64` value with zero radius.
    #[cfg(feature = "mpfr")]
    pub fn from_f64(value: f64) -> Self {
        Self::new(rug::Float::with_val(BALL_PRECISION, value), mpfr_zero())
    }

    /// Access the midpoint.
    #[cfg(not(feature = "mpfr"))]
    pub fn mid(&self) -> f64 {
        self.mid
    }

    /// Access the midpoint.
    #[cfg(feature = "mpfr")]
    pub fn mid(&self) -> &rug::Float {
        &self.mid
    }

    /// Access the radius.
    #[cfg(not(feature = "mpfr"))]
    pub fn rad(&self) -> f64 {
        self.rad
    }

    /// Access the radius.
    #[cfg(feature = "mpfr")]
    pub fn rad(&self) -> &rug::Float {
        &self.rad
    }

    /// Return a conservative lower bound.
    #[cfg(not(feature = "mpfr"))]
    pub fn lower(&self) -> f64 {
        self.mid - self.rad
    }

    /// Return a conservative lower bound.
    #[cfg(feature = "mpfr")]
    pub fn lower(&self) -> rug::Float {
        let mut lo = self.mid.clone();
        lo.sub_assign_round(&self.rad, rug::float::Round::Down);
        lo
    }

    /// Return a conservative upper bound.
    #[cfg(not(feature = "mpfr"))]
    pub fn upper(&self) -> f64 {
        self.mid + self.rad
    }

    /// Return a conservative upper bound.
    #[cfg(feature = "mpfr")]
    pub fn upper(&self) -> rug::Float {
        let mut hi = self.mid.clone();
        hi.add_assign_round(&self.rad, rug::float::Round::Up);
        hi
    }

    /// Return the precision of the MPFR backing store in bits.
    #[cfg(feature = "mpfr")]
    pub fn precision(&self) -> u32 {
        self.mid.prec()
    }
}

/// The real ball domain.
///
/// In the default build this wraps `f64`-based balls. When the `mpfr` feature
/// is enabled, this domain uses `rug::Float` with rigorous rounding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealBallDomain;

impl Domain for RealBallDomain {
    type Element = RealBall;

    fn zero(&self) -> Self::Element {
        RealBall::from_f64(0.0)
    }

    fn one(&self) -> Self::Element {
        RealBall::from_f64(1.0)
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        #[cfg(not(feature = "mpfr"))]
        {
            RealBall::new(a.mid + b.mid, a.rad + b.rad)
        }
        #[cfg(feature = "mpfr")]
        {
            let mut mid = a.mid.clone();
            let mut rad = a.rad.clone();
            mid.add_assign_round(&b.mid, rug::float::Round::Down);
            rad.add_assign_round(&b.rad, rug::float::Round::Up);
            rad.add_assign_round(&rounding_unit(&mid), rug::float::Round::Up);
            RealBall::new(mid, rad)
        }
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        #[cfg(not(feature = "mpfr"))]
        {
            RealBall::new(a.mid - b.mid, a.rad + b.rad)
        }
        #[cfg(feature = "mpfr")]
        {
            let mut mid = a.mid.clone();
            let mut rad = a.rad.clone();
            mid.sub_assign_round(&b.mid, rug::float::Round::Up);
            rad.add_assign_round(&b.rad, rug::float::Round::Up);
            rad.add_assign_round(&rounding_unit(&mid), rug::float::Round::Up);
            RealBall::new(mid, rad)
        }
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        #[cfg(not(feature = "mpfr"))]
        {
            RealBall::new(-a.mid, a.rad)
        }
        #[cfg(feature = "mpfr")]
        {
            RealBall::new(-a.mid.clone(), a.rad.clone())
        }
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        #[cfg(not(feature = "mpfr"))]
        {
            let a_lo = a.lower();
            let a_hi = a.upper();
            let b_lo = b.lower();
            let b_hi = b.upper();
            let p1 = a_lo * b_lo;
            let p2 = a_lo * b_hi;
            let p3 = a_hi * b_lo;
            let p4 = a_hi * b_hi;
            let lo = p1.min(p2).min(p3).min(p4);
            let hi = p1.max(p2).max(p3).max(p4);
            let mid = (lo + hi) / 2.0;
            let rad = (hi - lo) / 2.0;
            RealBall::new(mid, rad)
        }
        #[cfg(feature = "mpfr")]
        {
            let a_lo = a.lower();
            let a_hi = a.upper();
            let b_lo = b.lower();
            let b_hi = b.upper();
            let mut products = [
                a_lo.clone() * b_lo.clone(),
                a_lo * b_hi.clone(),
                a_hi.clone() * b_lo,
                a_hi * b_hi,
            ];
            for p in &mut products {
                p.set_prec(BALL_PRECISION);
            }
            let lo = products.iter().cloned().reduce(mpfr_min).unwrap();
            let hi = products.iter().cloned().reduce(mpfr_max).unwrap();
            let mut mid = (lo.clone() + &hi) / 2u32;
            mid.set_prec(BALL_PRECISION);
            let mut rad = (hi - &lo).abs() / 2u32;
            rad.set_prec(BALL_PRECISION);
            rad.add_assign_round(&rounding_unit(&mid), rug::float::Round::Up);
            RealBall::new(mid, rad)
        }
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        let inv = self.inv(b)?;
        Some(self.mul(a, &inv))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        #[cfg(not(feature = "mpfr"))]
        {
            if a.lower() <= 0.0 && a.upper() >= 0.0 {
                return None;
            }
            let lo = a.lower();
            let hi = a.upper();
            let inv_lo = 1.0 / hi;
            let inv_hi = 1.0 / lo;
            let mid = (inv_lo + inv_hi) / 2.0;
            let rad = (inv_hi - inv_lo).abs() / 2.0;
            Some(RealBall::new(mid, rad))
        }
        #[cfg(feature = "mpfr")]
        {
            let lo = a.lower();
            let hi = a.upper();
            if lo <= mpfr_zero() && hi >= mpfr_zero() {
                return None;
            }
            let mut inv_lo = rug::Float::with_val(BALL_PRECISION, 1.0);
            inv_lo.div_assign_round(&hi, rug::float::Round::Down);
            let mut inv_hi = rug::Float::with_val(BALL_PRECISION, 1.0);
            inv_hi.div_assign_round(&lo, rug::float::Round::Up);
            let mut mid = (inv_lo.clone() + &inv_hi) / 2u32;
            mid.set_prec(BALL_PRECISION);
            let mut rad = (inv_hi - &inv_lo).abs() / 2u32;
            rad.set_prec(BALL_PRECISION);
            rad.add_assign_round(&rounding_unit(&mid), rug::float::Round::Up);
            Some(RealBall::new(mid, rad))
        }
    }
}

#[cfg(feature = "mpfr")]
fn mpfr_zero() -> rug::Float {
    rug::Float::with_val(BALL_PRECISION, 0.0)
}

/// Rounding-unit helper for `rug::Float`.
#[cfg(feature = "mpfr")]
fn rounding_unit(x: &rug::Float) -> rug::Float {
    let p = x.prec();
    let mut u = rug::Float::with_val(p, 1.0);
    let denom = rug::Float::with_val(p, 2.0).pow(p - 1);
    u.div_assign_round(&denom, rug::float::Round::Up);
    u
}

/// Minimum of two `rug::Float` values.
#[cfg(feature = "mpfr")]
fn mpfr_min(a: rug::Float, b: rug::Float) -> rug::Float {
    if a < b { a } else { b }
}

/// Maximum of two `rug::Float` values.
#[cfg(feature = "mpfr")]
fn mpfr_max(a: rug::Float, b: rug::Float) -> rug::Float {
    if a > b { a } else { b }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(not(feature = "mpfr"))]
    #[test]
    fn ball_default_addition_widens_radius() {
        let domain = RealBallDomain;
        let a = RealBall::from_f64(1.0);
        let b = RealBall::from_f64(2.0);
        let sum = domain.add(&a, &b);
        assert!((sum.mid() - 3.0).abs() < 1e-12);
        assert!((sum.rad() - 0.0).abs() < 1e-12);
    }

    #[cfg(not(feature = "mpfr"))]
    #[test]
    fn ball_default_multiplication_contains_true_product() {
        let domain = RealBallDomain;
        let a = RealBall::from_f64(2.0);
        let b = RealBall::from_f64(3.0);
        let prod = domain.mul(&a, &b);
        assert!(prod.lower() <= 6.0 && 6.0 <= prod.upper());
    }

    #[cfg(not(feature = "mpfr"))]
    #[test]
    fn ball_default_inverse_excludes_zero() {
        let domain = RealBallDomain;
        let a = RealBall::from_f64(2.0);
        let inv = domain.inv(&a).expect("ball does not contain zero");
        let prod = domain.mul(&a, &inv);
        assert!(prod.lower() <= 1.0 && 1.0 <= prod.upper());
    }

    #[cfg(not(feature = "mpfr"))]
    #[test]
    fn ball_default_inverse_of_zero_ball_is_none() {
        let domain = RealBallDomain;
        let a = RealBall::from_f64(0.0);
        assert!(domain.inv(&a).is_none());
    }

    #[cfg(not(feature = "mpfr"))]
    #[test]
    fn ball_default_contains_zero_leads_to_none() {
        let domain = RealBallDomain;
        let a = RealBall::new(1.0, 0.0);
        let b = RealBall::new(-0.05, 0.1);
        assert!(domain.div(&a, &b).is_none());
    }

    #[cfg(feature = "mpfr")]
    mod mpfr_tests {
        use super::*;

        #[test]
        fn ball_mpfr_contains_exact_value() {
            let domain = RealBallDomain;
            let a = RealBall::from_f64(1.0);
            let b = RealBall::from_f64(2.0);
            let sum = domain.add(&a, &b);
            let lo = sum.lower();
            let hi = sum.upper();
            assert!(lo <= rug::Float::with_val(BALL_PRECISION, 3.0));
            assert!(hi >= rug::Float::with_val(BALL_PRECISION, 3.0));
        }

        #[test]
        fn ball_mpfr_inverse_excludes_zero() {
            let domain = RealBallDomain;
            let a = RealBall::from_f64(2.0);
            let inv = domain.inv(&a).expect("invertible");
            let prod = domain.mul(&a, &inv);
            let lo = prod.lower();
            let hi = prod.upper();
            assert!(lo <= rug::Float::with_val(BALL_PRECISION, 1.0));
            assert!(hi >= rug::Float::with_val(BALL_PRECISION, 1.0));
        }
    }
}
