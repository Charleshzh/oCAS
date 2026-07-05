//! Fast polynomial evaluation using Estrin's scheme.
//!
//! This module wraps the [`fast_polynomial`] crate to provide
//! instruction-level parallelism (ILP) for evaluating polynomials
//! at a single `f64` point. Estrin's scheme reorganizes the
//! computation so that independent sub-expressions can execute
//! in parallel on modern CPUs.
//!
//! Enabled with the `fast-poly` feature flag.

use fast_polynomial;

/// Evaluate a polynomial at `x` using Estrin's scheme.
///
/// `coeffs` are given from constant term upward: `a[0] + a[1]*x + a[2]*x^2 + ...`
///
/// This is typically 1.5–2× faster than Horner's method for degree ≥ 8
/// on CPUs with FMA support, due to instruction-level parallelism.
///
/// # Example
///
/// ```
/// use ocas_eval::poly_eval::eval_estrin;
///
/// // 1 + 2x + 3x^2 at x = 2.0 → 1 + 4 + 12 = 17
/// let result = eval_estrin(&[1.0, 2.0, 3.0], 2.0);
/// assert!((result - 17.0).abs() < 1e-10);
/// ```
pub fn eval_estrin(coeffs: &[f64], x: f64) -> f64 {
    match coeffs.len() {
        0 => 0.0,
        1 => coeffs[0],
        2 => coeffs[0] + coeffs[1] * x,
        _ => fast_polynomial::poly(x, coeffs),
    }
}

/// Evaluate a polynomial at multiple `x` values using Estrin's scheme.
///
/// Returns a vector of results, one per input value.
pub fn eval_estrin_batch(coeffs: &[f64], xs: &[f64]) -> Vec<f64> {
    xs.iter().map(|&x| eval_estrin(coeffs, x)).collect()
}

/// Horner's method for comparison/baseline.
pub fn eval_horner(coeffs: &[f64], x: f64) -> f64 {
    let mut result = 0.0f64;
    for &c in coeffs.iter().rev() {
        result = result * x + c;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estrin_constant() {
        assert!((eval_estrin(&[42.0], 7.0) - 42.0).abs() < 1e-15);
    }

    #[test]
    fn estrin_linear() {
        // 1 + 2x at x = 3 → 7
        assert!((eval_estrin(&[1.0, 2.0], 3.0) - 7.0).abs() < 1e-15);
    }

    #[test]
    fn estrin_quadratic() {
        // 1 + 2x + 3x^2 at x = 2 → 1 + 4 + 12 = 17
        assert!((eval_estrin(&[1.0, 2.0, 3.0], 2.0) - 17.0).abs() < 1e-15);
    }

    #[test]
    fn estrin_vs_horner_random() {
        // Cross-check Estrin against Horner for a degree-15 polynomial
        let coeffs: Vec<f64> = (0..16).map(|i| (i as f64) * 0.1 + 1.0).collect();
        let x = 1.5;
        let estrin_result = eval_estrin(&coeffs, x);
        let horner_result = eval_horner(&coeffs, x);
        assert!(
            (estrin_result - horner_result).abs() < 1e-10,
            "estrin={estrin_result}, horner={horner_result}"
        );
    }

    #[test]
    fn estrin_batch() {
        let coeffs = vec![1.0, 2.0, 3.0]; // 1 + 2x + 3x^2
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let results = eval_estrin_batch(&coeffs, &xs);
        assert_eq!(results.len(), 4);
        assert!((results[0] - 1.0).abs() < 1e-15);
        assert!((results[1] - 6.0).abs() < 1e-15);
        assert!((results[2] - 17.0).abs() < 1e-15);
        assert!((results[3] - 34.0).abs() < 1e-15);
    }

    #[test]
    fn estrin_empty() {
        assert_eq!(eval_estrin(&[], 5.0), 0.0);
    }
}
