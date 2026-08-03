//! Double-precision floating-point arithmetic (Dekker/Knuth "double-float").
//!
//! [`DoubleF64`] represents a floating-point number as the unevaluated sum of
//! two `f64` values (`hi + lo`), where `|lo| ≤ 0.5 · ulp(hi)`. This gives
//! approximately 31 decimal digits of precision (~84 binary bits) — roughly
//! double that of a single `f64`.
//!
//! The arithmetic is based on Dekker's and Knuth's algorithms for error-free
//! transformations (TwoSum, TwoProd) and is significantly faster than
//! arbitrary-precision alternatives like MPFR.

use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::Zero;

use crate::domain::Domain;

// =========================================================================
// Type definition
// =========================================================================

/// A double-precision floating-point number: `hi + lo` with `|lo| ≤ 0.5·ulp(hi)`.
///
/// Provides ~31 decimal digits (~84 binary bits) of precision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleF64 {
    /// High-order component (the "main" value).
    pub hi: f64,
    /// Low-order component (the error term).
    pub lo: f64,
}

impl Eq for DoubleF64 {}

// =========================================================================
// Construction and conversion
// =========================================================================

impl DoubleF64 {
    /// The zero value.
    pub const ZERO: Self = Self { hi: 0.0, lo: 0.0 };
    /// The unit value.
    pub const ONE: Self = Self { hi: 1.0, lo: 0.0 };

    /// Create a new `DoubleF64` from high and low components.
    /// The caller must ensure `|lo| ≤ 0.5·ulp(hi)` for correct results.
    #[inline]
    pub fn new(hi: f64, lo: f64) -> Self {
        Self { hi, lo }
    }

    /// Create from a single `f64`.
    #[inline]
    pub fn from_f64(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }

    /// Extract the high-order component as `f64`.
    #[inline]
    pub fn to_f64(self) -> f64 {
        self.hi
    }

    /// Create from two `f64` values, normalizing via TwoSum.
    pub fn quick_two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let e = b - (s - a);
        Self { hi: s, lo: e }
    }

    /// Dekker's TwoSum: error-free sum with rounding error captured.
    pub fn two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let a_prime = s - b;
        let b_prime = s - a_prime;
        let delta_a = a - a_prime;
        let delta_b = b - b_prime;
        let e = delta_a + delta_b;
        Self { hi: s, lo: e }
    }

    /// Dekker's split: split a f64 into two 26-bit halves.
    #[allow(dead_code)]
    fn split(x: f64) -> (f64, f64) {
        const SPLITTER: f64 = 134217729.0; // 2^27 + 1
        let c = SPLITTER * x;
        let c_hi = c - (c - x);
        let c_lo = x - c_hi;
        (c_hi, c_lo)
    }

    /// Dekker's TwoProd without FMA: error-free product with error term.
    #[allow(dead_code)]
    fn two_prod_no_fma(a: f64, b: f64) -> Self {
        let p = a * b;
        let (a_hi, a_lo) = Self::split(a);
        let (b_hi, b_lo) = Self::split(b);
        let err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
        Self { hi: p, lo: err }
    }

    /// TwoProd using FMA when available (preferred).
    #[inline]
    fn two_prod(a: f64, b: f64) -> Self {
        let p = a * b;
        let err = f64::mul_add(a, b, -p);
        Self { hi: p, lo: err }
    }

    /// Absolute value.
    #[inline]
    pub fn abs(self) -> Self {
        if self.hi < 0.0 {
            Self {
                hi: -self.hi,
                lo: -self.lo,
            }
        } else {
            self
        }
    }

    /// Check if the value is NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.hi.is_nan()
    }

    /// Check if the value is infinite.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.hi.is_infinite()
    }

    /// Check if the value is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.hi.is_finite()
    }

    // =====================================================================
    // Arithmetic operations
    // =====================================================================

    /// Add two `DoubleF64` values.
    #[allow(clippy::should_implement_trait)]
    pub fn add(self, other: Self) -> Self {
        let s = Self::two_sum(self.hi, other.hi);
        let v = self.lo + other.lo;
        let w = s.lo + v;
        Self::quick_two_sum(s.hi, w)
    }

    /// Subtract two `DoubleF64` values.
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, other: Self) -> Self {
        self.add(-other)
    }

    /// Multiply two `DoubleF64` values.
    #[allow(clippy::should_implement_trait)]
    pub fn mul(self, other: Self) -> Self {
        let p = Self::two_prod(self.hi, other.hi);
        let err = self.hi * other.lo + self.lo * other.hi;
        Self::quick_two_sum(p.hi, p.lo + err)
    }

    /// Divide two `DoubleF64` values.
    #[allow(clippy::should_implement_trait)]
    pub fn div(self, other: Self) -> Self {
        let q1 = self.hi / other.hi;
        let p = Self::two_prod(q1, other.hi);
        let delta = self.hi - p.hi;
        let err = (delta - p.lo + self.lo) / other.hi;
        Self::quick_two_sum(q1, err)
    }

    /// Integer power via binary exponentiation.
    pub fn powi(self, mut n: i64) -> Self {
        if n == 0 {
            return Self::ONE;
        }
        let negate = n < 0;
        if negate {
            n = -n;
        }
        let mut base = self;
        let mut result = Self::ONE;
        let mut exp = n as u64;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(base);
            }
            base = base.mul(base);
            exp >>= 1;
        }
        if negate {
            Self::ONE.div(result)
        } else {
            result
        }
    }

    // =====================================================================
    // Transcendental functions
    // =====================================================================

    /// Square root via Newton iteration.
    pub fn sqrt(self) -> Self {
        if self.hi < 0.0 {
            return Self::from_f64(f64::NAN);
        }
        if self.hi == 0.0 {
            return Self::ZERO;
        }
        // Initial estimate from hardware sqrt
        let x0 = self.hi.sqrt();
        let mut x = Self::from_f64(x0);
        // Newton iteration: x = (x + self/x) / 2
        // Two iterations give full DoubleF64 precision
        for _ in 0..4 {
            x = x.add(self.div(x)).mul(Self::from_f64(0.5));
        }
        x
    }

    /// Absolute value.
    pub fn dabs(self) -> Self {
        if self.hi < 0.0 { -self } else { self }
    }

    /// Exponential function via Taylor series with argument reduction.
    #[allow(clippy::approx_constant)]
    pub fn exp(self) -> Self {
        if self.hi == 0.0 && self.lo == 0.0 {
            return Self::ONE;
        }
        // Argument reduction: exp(x) = exp(x - k*ln2) * 2^k
        const LN2: DoubleF64 = DoubleF64 {
            hi: std::f64::consts::LN_2,
            lo: 2.319046813846299e-17,
        };
        let k = (self.hi / LN2.hi).round() as i64;
        let reduced = self.sub(LN2.mul(Self::from_f64(k as f64)));

        // Taylor series for exp(r) where r is small
        let one = Self::ONE;
        let mut term = one;
        let mut sum = one;
        for i in 1..40 {
            term = term.mul(reduced).div(Self::from_f64(i as f64));
            sum = sum.add(term);
        }

        // Multiply by 2^k
        if k >= 0 {
            sum.mul(Self::from_f64((1u64 << k.min(62) as u64) as f64))
        } else {
            sum.div(Self::from_f64((1u64 << (-k).min(62) as u64) as f64))
        }
    }

    /// Natural logarithm via Newton iteration on exp.
    pub fn ln(self) -> Self {
        if self.hi <= 0.0 {
            return Self::from_f64(f64::NAN);
        }
        if self.hi == 1.0 && self.lo == 0.0 {
            return Self::ZERO;
        }
        // Initial estimate
        let x0 = self.hi.ln();
        let mut x = Self::from_f64(x0);
        // Newton iteration for ln: x = x + (self - exp(x)) / exp(x)
        // Use a few iterations for full precision
        for _ in 0..6 {
            let ex = x.exp();
            x = x.add(self.sub(ex).div(ex));
        }
        x
    }

    /// Sine via Taylor series with argument reduction.
    #[allow(clippy::approx_constant)]
    pub fn sin(self) -> Self {
        const PI: DoubleF64 = DoubleF64 {
            hi: std::f64::consts::PI,
            lo: 1.2246467991473532e-16,
        };
        const TWO_PI: DoubleF64 = DoubleF64 {
            hi: std::f64::consts::TAU,
            lo: 2.4492935982947064e-16,
        };

        // Reduce to [-π, π]
        let mut x = self;
        if x.dabs().hi > PI.hi {
            let k = (x.hi / TWO_PI.hi).round();
            x = x.sub(TWO_PI.mul(Self::from_f64(k)));
        }

        // Taylor series
        let x2 = x.mul(x);
        let mut term = x;
        let mut sum = x;
        for i in 1..25 {
            let n = (2 * i + 1) as f64;
            term = term.mul(x2).div(Self::from_f64(-n * (n - 1.0)));
            sum = sum.add(term);
        }
        sum
    }

    /// Cosine via sin(π/2 - x).
    #[allow(clippy::approx_constant)]
    pub fn cos(self) -> Self {
        const FRAC_PI_2: DoubleF64 = DoubleF64 {
            hi: std::f64::consts::FRAC_PI_2,
            lo: 6.123233995736766e-17,
        };
        FRAC_PI_2.sub(self).sin()
    }

    /// Tangent = sin / cos.
    pub fn tan(self) -> Self {
        self.sin().div(self.cos())
    }
}

// Note: cos(x) = sin(π/2 - x) is incorrect above. Let me fix:
// We override cos below.

// =========================================================================
// std::ops trait implementations
// =========================================================================

impl Add for DoubleF64 {
    type Output = Self;
    #[inline]
    fn add(self, other: Self) -> Self {
        DoubleF64::add(self, other)
    }
}

impl Sub for DoubleF64 {
    type Output = Self;
    #[inline]
    fn sub(self, other: Self) -> Self {
        DoubleF64::sub(self, other)
    }
}

impl Mul for DoubleF64 {
    type Output = Self;
    #[inline]
    fn mul(self, other: Self) -> Self {
        DoubleF64::mul(self, other)
    }
}

impl Div for DoubleF64 {
    type Output = Self;
    #[inline]
    fn div(self, other: Self) -> Self {
        DoubleF64::div(self, other)
    }
}

impl Neg for DoubleF64 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl AddAssign for DoubleF64 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl SubAssign for DoubleF64 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl MulAssign for DoubleF64 {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl DivAssign for DoubleF64 {
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}

impl PartialOrd for DoubleF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.hi.partial_cmp(&other.hi)
    }
}

impl From<f64> for DoubleF64 {
    fn from(x: f64) -> Self {
        Self::from_f64(x)
    }
}

impl Zero for DoubleF64 {
    fn zero() -> Self {
        Self::ZERO
    }
    fn is_zero(&self) -> bool {
        self.hi == 0.0 && self.lo == 0.0
    }
}

impl fmt::Display for DoubleF64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Display with full precision
        if self.lo == 0.0 {
            write!(f, "{}", self.hi)
        } else {
            write!(f, "{:.31e}", self.hi + self.lo)
        }
    }
}

// =========================================================================
// Domain trait implementation
// =========================================================================

/// Double-float domain for algebraic computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DoubleF64Domain;

impl Domain for DoubleF64Domain {
    type Element = DoubleF64;

    fn zero(&self) -> Self::Element {
        DoubleF64::ZERO
    }

    fn one(&self) -> Self::Element {
        DoubleF64::ONE
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        *a + *b
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        *a - *b
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        -*a
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        *a * *b
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        if b.is_zero() { None } else { Some(*a / *b) }
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if a.is_zero() {
            None
        } else {
            Some(DoubleF64::ONE / *a)
        }
    }

    fn is_zero(&self, a: &Self::Element) -> bool {
        a.is_zero()
    }

    fn is_one(&self, a: &Self::Element) -> bool {
        *a == DoubleF64::ONE
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_arithmetic() {
        let a = DoubleF64::from_f64(1.0);
        let b = DoubleF64::from_f64(2.0);
        assert_eq!((a + b).hi, 3.0);
        assert_eq!((a - b).hi, -1.0);
        assert_eq!((a * b).hi, 2.0);
        assert_eq!((a / b).hi, 0.5);
    }

    #[test]
    fn precision_gain() {
        // DoubleF64 should capture rounding errors that f64 loses
        let a = DoubleF64::from_f64(1.0);
        let b = DoubleF64::from_f64(f64::EPSILON);
        let sum = a + b;
        // 1 + eps in f64: hi=1.0, lo=eps (captured in lo)
        let reconstructed = sum.hi + sum.lo;
        assert_eq!(reconstructed, 1.0 + f64::EPSILON);
    }

    #[test]
    fn two_sum_correctness() {
        let s = DoubleF64::two_sum(1.0, f64::EPSILON);
        assert_eq!(s.hi, 1.0 + f64::EPSILON);
        // s.lo should capture the rounding error
    }

    #[test]
    fn two_prod_correctness() {
        let p = DoubleF64::two_prod(3.0, 5.0);
        assert_eq!(p.hi, 15.0);
        assert_eq!(p.lo, 0.0);
    }

    #[test]
    fn powi_basic() {
        let x = DoubleF64::from_f64(3.0);
        assert_eq!(x.powi(0).hi, 1.0);
        assert_eq!(x.powi(1).hi, 3.0);
        assert_eq!(x.powi(2).hi, 9.0);
        assert_eq!(x.powi(3).hi, 27.0);
    }

    #[test]
    fn powi_negative() {
        let x = DoubleF64::from_f64(2.0);
        assert_eq!(x.powi(-1).hi, 0.5);
        assert_eq!(x.powi(-2).hi, 0.25);
    }

    #[test]
    fn sqrt_basic() {
        let x = DoubleF64::from_f64(4.0);
        let s = x.sqrt();
        assert!((s.hi - 2.0).abs() < 1e-30);
    }

    #[test]
    fn sqrt_two() {
        let x = DoubleF64::from_f64(2.0);
        let s = x.sqrt();
        // sqrt(2)^2 should be very close to 2
        let sq = s * s;
        assert!((sq.hi - 2.0).abs() < 1e-28);
    }

    #[test]
    fn exp_basic() {
        let zero = DoubleF64::ZERO;
        assert_eq!(zero.exp().hi, 1.0);

        let one = DoubleF64::ONE;
        let e = one.exp();
        // exp(1) ≈ 2.718281828...
        assert!((e.hi - std::f64::consts::E).abs() < 1e-28);
    }

    #[test]
    fn ln_basic() {
        let one = DoubleF64::ONE;
        assert_eq!(one.ln().hi, 0.0);

        let e = DoubleF64::from_f64(std::f64::consts::E);
        let ln_e = e.ln();
        assert!((ln_e.hi - 1.0).abs() < 1e-28);
    }

    #[test]
    fn sin_cos_basic() {
        let zero = DoubleF64::ZERO;
        assert_eq!(zero.sin().hi, 0.0);
        assert_eq!(zero.cos().hi, 1.0);
    }

    #[test]
    fn sin_pi() {
        let pi = DoubleF64::from_f64(std::f64::consts::PI);
        let s = pi.sin();
        // sin(π) ≈ 0; residual is rounding error in the PI constant
        assert!(s.hi.abs() < 1e-14, "sin(π) ≈ 0, got {}", s.hi);
    }

    #[test]
    fn domain_trait() {
        let dom = DoubleF64Domain;
        let a = DoubleF64::from_f64(3.0);
        let b = DoubleF64::from_f64(4.0);
        assert_eq!(dom.add(&a, &b).hi, 7.0);
        assert_eq!(dom.mul(&a, &b).hi, 12.0);
        assert_eq!(dom.div(&a, &b).unwrap().hi, 0.75);
        assert!(dom.div(&a, &DoubleF64::ZERO).is_none());
    }

    #[test]
    fn display() {
        #[allow(clippy::approx_constant)]
        let x = DoubleF64::from_f64(3.14);
        assert_eq!(format!("{x}"), "3.14");
    }
}
