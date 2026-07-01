//! Numeric evaluation domain trait.
//!
//! The [`EvaluationDomain`] trait abstracts over numeric types that can be
//! used for expression evaluation. It provides arithmetic operations and
//! a table of built-in mathematical functions (sin, cos, exp, etc.).
//!
//! # Implementations
//!
//! - `f64` — standard double-precision floating point (always available)
//! - `ocas_domain::Integer` — arbitrary-precision integers
//! - `ocas_domain::Rational` — arbitrary-precision rationals
//! - `ocas_domain::RealBall` — rigorous real interval arithmetic

use crate::error::{EvaluationError, Result};

/// Trait for types that can serve as evaluation domains.
///
/// Unlike the algebraic [`Domain`](ocas_domain::Domain) trait which is
/// object-safe and uses `&self` receivers, `EvaluationDomain` uses static
/// methods and `&self`/`&other` for value operations. This makes it
/// compatible with `f64` and other `Copy` types.
///
/// # Example
///
/// ```
/// use ocas_eval::EvaluationDomain;
///
/// let x = f64::from_f64(3.0);
/// let y = f64::from_f64(2.0);
/// assert_eq!(x.add_ref(&y), 5.0);
/// assert_eq!(f64::resolve_builtin("sin", &std::f64::consts::FRAC_PI_2).unwrap(), 1.0);
/// ```
pub trait EvaluationDomain: Sized + Clone + 'static {
    /// Create a value from an `f64`.
    fn from_f64(value: f64) -> Self;

    /// The additive identity (0).
    fn zero() -> Self;

    /// The multiplicative identity (1).
    fn one() -> Self;

    /// `self + other`
    fn add_ref(&self, other: &Self) -> Self;

    /// `self - other`
    fn sub_ref(&self, other: &Self) -> Self;

    /// `self * other`
    fn mul_ref(&self, other: &Self) -> Self;

    /// `self / other`, or [`EvaluationError::DivisionByZero`] if `other` is zero.
    fn div_ref(&self, other: &Self) -> Result<Self>;

    /// `-self`
    fn neg_ref(&self) -> Self;

    /// `self^exp` for integer exponents. Returns 1 when exp == 0.
    fn powi_ref(&self, exp: i64) -> Self;

    /// Resolve a built-in mathematical function.
    ///
    /// Accepts both lowercase (`sin`) and capitalized (`Sin`) names.
    /// The following functions are supported:
    ///
    /// | Function | Description |
    /// |---|---|
    /// | `sin` / `Sin` | Sine |
    /// | `cos` / `Cos` | Cosine |
    /// | `tan` / `Tan` | Tangent |
    /// | `sec` / `Sec` | Secant |
    /// | `csc` / `Csc` | Cosecant |
    /// | `cot` / `Cot` | Cotangent |
    /// | `exp` / `Exp` | Exponential (eˣ) |
    /// | `log` / `Log` | Natural logarithm |
    /// | `sqrt` / `Sqrt` | Square root |
    /// | `abs` / `Abs` | Absolute value |
    fn resolve_builtin(name: &str, arg: &Self) -> Result<Self>;
}

// ---------------------------------------------------------------------------
// f64 implementation
// ---------------------------------------------------------------------------

impl EvaluationDomain for f64 {
    #[inline]
    fn from_f64(value: f64) -> Self {
        value
    }

    #[inline]
    fn zero() -> Self {
        0.0
    }

    #[inline]
    fn one() -> Self {
        1.0
    }

    #[inline]
    fn add_ref(&self, other: &Self) -> Self {
        self + other
    }

    #[inline]
    fn sub_ref(&self, other: &Self) -> Self {
        self - other
    }

    #[inline]
    fn mul_ref(&self, other: &Self) -> Self {
        self * other
    }

    #[inline]
    fn div_ref(&self, other: &Self) -> Result<Self> {
        if *other == 0.0 {
            Err(EvaluationError::DivisionByZero)
        } else {
            Ok(self / other)
        }
    }

    #[inline]
    fn neg_ref(&self) -> Self {
        -self
    }

    #[inline]
    fn powi_ref(&self, exp: i64) -> Self {
        self.powi(exp as i32)
    }

    fn resolve_builtin(name: &str, arg: &Self) -> Result<Self> {
        match name.to_lowercase().as_str() {
            "sin" => Ok(arg.sin()),
            "cos" => Ok(arg.cos()),
            "tan" => Ok(arg.tan()),
            "sec" => Ok(1.0 / arg.cos()),
            "csc" => Ok(1.0 / arg.sin()),
            "cot" => Ok(1.0 / arg.tan()),
            "exp" => Ok(arg.exp()),
            "log" => {
                if *arg <= 0.0 {
                    Err(EvaluationError::UnsupportedOperation {
                        message: "log of non-positive number".into(),
                    })
                } else {
                    Ok(arg.ln())
                }
            }
            "sqrt" => {
                if *arg < 0.0 {
                    Err(EvaluationError::UnsupportedOperation {
                        message: "sqrt of negative number".into(),
                    })
                } else {
                    Ok(arg.sqrt())
                }
            }
            "abs" => Ok(arg.abs()),
            _ => Err(EvaluationError::FunctionNotFound {
                name: name.to_string(),
            }),
        }
    }
}

// ---------------------------------------------------------------------------
// PowfExtension
// ---------------------------------------------------------------------------

/// Extension to [`EvaluationDomain`] for floating-point exponentiation.
///
/// This is split from the main trait because integer domains cannot
/// meaningfully compute `a^b` for non-integer `b`.
pub trait PowfExtension: EvaluationDomain {
    /// `self^exp` for floating-point exponents.
    fn powf_ref(&self, exp: &Self) -> Result<Self>;
}

impl PowfExtension for f64 {
    fn powf_ref(&self, exp: &Self) -> Result<Self> {
        Ok(self.powf(*exp))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_arithmetic() {
        assert_eq!(f64::zero(), 0.0);
        assert_eq!(f64::one(), 1.0);
        assert_eq!(3.0f64.add_ref(&2.0), 5.0);
        assert_eq!(3.0f64.sub_ref(&2.0), 1.0);
        assert_eq!(3.0f64.mul_ref(&2.0), 6.0);
        assert_eq!(6.0f64.div_ref(&2.0).unwrap(), 3.0);
        assert!(6.0f64.div_ref(&0.0).is_err());
        assert_eq!(3.0f64.neg_ref(), -3.0);
        assert_eq!(2.0f64.powi_ref(3), 8.0);
        assert_eq!(2.0f64.powi_ref(0), 1.0);
    }

    #[test]
    fn f64_builtin_sin_lowercase() {
        let result = f64::resolve_builtin("sin", &std::f64::consts::FRAC_PI_2).unwrap();
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_sin_capitalized() {
        let result = f64::resolve_builtin("Sin", &std::f64::consts::FRAC_PI_2).unwrap();
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_cos() {
        let result = f64::resolve_builtin("cos", &std::f64::consts::PI).unwrap();
        assert!((result + 1.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_exp() {
        let result = f64::resolve_builtin("exp", &1.0).unwrap();
        assert!((result - std::f64::consts::E).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_log() {
        let result = f64::resolve_builtin("log", &std::f64::consts::E).unwrap();
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_log_negative() {
        assert!(f64::resolve_builtin("log", &(-1.0)).is_err());
    }

    #[test]
    fn f64_builtin_sqrt() {
        let result = f64::resolve_builtin("sqrt", &4.0).unwrap();
        assert!((result - 2.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_sqrt_negative() {
        assert!(f64::resolve_builtin("sqrt", &(-1.0)).is_err());
    }

    #[test]
    fn f64_builtin_abs() {
        assert_eq!(f64::resolve_builtin("abs", &(-3.0)).unwrap(), 3.0);
        assert_eq!(f64::resolve_builtin("abs", &3.0).unwrap(), 3.0);
    }

    #[test]
    fn f64_builtin_tan() {
        let result = f64::resolve_builtin("tan", &0.0).unwrap();
        assert!((result - 0.0).abs() < 1e-10);
    }

    #[test]
    fn f64_builtin_unknown() {
        assert!(f64::resolve_builtin("unknown_fn", &0.0).is_err());
    }

    #[test]
    fn f64_from_f64() {
        assert_eq!(f64::from_f64(42.0), 42.0);
    }
}
