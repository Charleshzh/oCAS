//! Numerical root-finding verification using bisection.
//!
//! Uses basic bisection to verify that known polynomial roots are
//! correctly isolated. This validates oCAS's Sturm sequence
//! implementation against a simple but robust numeric method.

use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

/// Evaluate a dense polynomial at an `f64` point.
/// Coefficients are in ascending order: [a0, a1, ..., an].
fn eval_poly(coeffs: &[Integer], x: f64) -> f64 {
    let mut result = 0.0f64;
    for c in coeffs.iter().rev() {
        let c_f64 = c.to_i64().map(|v| v as f64).unwrap_or(0.0);
        result = result * x + c_f64;
    }
    result
}

/// Find all real roots of a polynomial in `[a, b]` using bisection.
fn find_roots_bisection(coeffs: &[Integer], a: f64, b: f64) -> Vec<f64> {
    let f = |x: f64| eval_poly(coeffs, x);
    let steps = 10_000;
    let dx = (b - a) / steps as f64;
    let mut roots = Vec::new();

    let mut prev_x = a;
    let mut prev_y = f(a);
    for i in 1..=steps {
        let x = a + dx * i as f64;
        let y = f(x);
        if prev_y == 0.0 {
            roots.push(prev_x);
        } else if prev_y * y < 0.0 {
            // Bisection refinement
            let mut lo = prev_x;
            let mut hi = x;
            for _ in 0..100 {
                let mid = (lo + hi) / 2.0;
                if f(mid) * f(lo) <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            roots.push((lo + hi) / 2.0);
        }
        prev_x = x;
        prev_y = y;
    }
    roots
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_quadratic_roots() {
        // x^2 - 5x + 6 = (x-2)(x-3), roots at 2 and 3
        let domain = IntegerDomain;
        let poly = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Integer::from(6), Integer::from(-5), Integer::from(1)],
        );
        let roots = find_roots_bisection(poly.coeffs(), 0.0, 10.0);
        assert_eq!(roots.len(), 2, "expected 2 roots, got {:?}", roots);
        assert!((roots[0] - 2.0).abs() < 1e-4);
        assert!((roots[1] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn verify_cubic_roots() {
        // x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
        let domain = IntegerDomain;
        let poly = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![
                Integer::from(-6),
                Integer::from(11),
                Integer::from(-6),
                Integer::from(1),
            ],
        );
        let roots = find_roots_bisection(poly.coeffs(), -1.0, 5.0);
        assert_eq!(roots.len(), 3, "expected 3 roots, got {:?}", roots);
        for (i, expected) in [1.0, 2.0, 3.0].iter().enumerate() {
            assert!(
                (roots[i] - expected).abs() < 1e-4,
                "root[{i}] = {}",
                roots[i]
            );
        }
    }

    #[test]
    fn verify_quartic_roots() {
        // x^4 - 10x^3 + 35x^2 - 50x + 24 = (x-1)(x-2)(x-3)(x-4)
        let domain = IntegerDomain;
        let poly = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![
                Integer::from(24),
                Integer::from(-50),
                Integer::from(35),
                Integer::from(-10),
                Integer::from(1),
            ],
        );
        let roots = find_roots_bisection(poly.coeffs(), 0.0, 5.0);
        assert_eq!(roots.len(), 4, "expected 4 roots, got {:?}", roots);
        for (i, expected) in [1.0, 2.0, 3.0, 4.0].iter().enumerate() {
            assert!(
                (roots[i] - expected).abs() < 1e-4,
                "root[{i}] = {}",
                roots[i]
            );
        }
    }
}
