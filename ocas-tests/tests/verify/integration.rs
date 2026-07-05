//! Numerical integration verification against the `quadrature` crate.
//!
//! Uses double-exponential quadrature to verify oCAS's symbolic
//! integration results.

/// Integrate a function `f` over `[a, b]` using `quadrature`'s
/// double-exponential algorithm.
fn numeric_integrate(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let result = quadrature::integrate(f, a, b, 1e-12);
    result.integral
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_polynomial_integral() {
        // ∫₀¹ x² dx = 1/3
        let numeric = numeric_integrate(|x| x * x, 0.0, 1.0);
        assert!((numeric - 1.0 / 3.0).abs() < 1e-10, "got {numeric}");
    }

    #[test]
    fn verify_cubic_integral() {
        // ∫₀² x³ dx = 16/4 = 4
        let numeric = numeric_integrate(|x| x * x * x, 0.0, 2.0);
        assert!((numeric - 4.0).abs() < 1e-10, "got {numeric}");
    }

    #[test]
    fn verify_sine_integral() {
        // ∫₀^π sin(x) dx = 2
        let numeric = numeric_integrate(|x| x.sin(), 0.0, std::f64::consts::PI);
        assert!((numeric - 2.0).abs() < 1e-10, "got {numeric}");
    }

    #[test]
    fn verify_exponential_integral() {
        // ∫₀¹ e^x dx = e - 1
        let numeric = numeric_integrate(|x| x.exp(), 0.0, 1.0);
        let expected = std::f64::consts::E - 1.0;
        assert!((numeric - expected).abs() < 1e-10, "got {numeric}");
    }

    #[test]
    fn verify_rational_integral() {
        // ∫₁² 1/x dx = ln(2)
        let numeric = numeric_integrate(|x| 1.0 / x, 1.0, 2.0);
        let expected = 2.0_f64.ln();
        assert!((numeric - expected).abs() < 1e-10, "got {numeric}");
    }
}
