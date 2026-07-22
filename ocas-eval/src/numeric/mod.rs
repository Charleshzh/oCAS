//! Numerical integration (adaptive Monte Carlo / Vegas) and deterministic
//! quadrature bridges.
//!
//! This module provides a [`Vegas`] integrator implementing the classic
//! Leppler/Kleisser/R. (Vegas) adaptive-grid Monte Carlo algorithm, together
//! with a convenience [`integrate_1d`] entry point. Both share the
//! [`Integrator`] trait so callers can swap methods uniformly.
//!
//! The integrands are plain `Fn(&[f64]) -> f64` closures; combine with the
//! crate's [`ExpressionEvaluator`](crate::ExpressionEvaluator) to integrate
//! symbolic expressions by wrapping `evaluate` in a closure.

pub mod statistics;
pub mod vegas;

pub use statistics::StatisticsAccumulator;
pub use vegas::{IntegrateResult, Integrator, Vegas, VegasOptions};

/// Numerically integrate a one-dimensional function `f` over `[a, b]` using
/// Vegas with default options. Returns the estimate and standard error.
///
/// The integrand receives `x` directly (not a unit-hypercube coordinate): the
/// linear change of variables is applied internally, so `jacobian = (b − a)`
/// is folded into the result.
///
/// ```
/// use ocas_eval::numeric::integrate_1d;
///
/// let r = integrate_1d(|x| x, 0.0, 1.0, Default::default());
/// assert!((r.integral - 0.5).abs() < 0.01);
/// ```
pub fn integrate_1d<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: VegasOptions,
) -> IntegrateResult {
    let width = b - a;
    let wrapped = move |u: &[f64]| f(a + u[0] * width) * width;
    let mut vegas = Vegas::new(1, opts);
    vegas.integrate(&wrapped)
}
