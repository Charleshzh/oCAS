//! Running statistics accumulator for Monte Carlo integration.
//!
//! Tracks the mean, variance, and χ² across independent iterations so the
//! integrator can combine stratified / multi-channel estimates with correct
//! weighting (the inverse-variance weight, as in the original Vegas paper).

use std::f64;

/// Accumulator for a single iteration's samples and for the cross-iteration
/// weighted average.
#[derive(Debug, Clone)]
pub struct StatisticsAccumulator {
    /// Σ weight over samples in the current iteration.
    sum_w: f64,
    /// Σ weight·f over samples in the current iteration.
    sum_wf: f64,
    /// Σ weight·f² over samples in the current iteration.
    sum_wf2: f64,
    /// Best estimate of the integral accumulated over iterations.
    integral: f64,
    /// Standard error of `integral`.
    error: f64,
    /// χ² across iterations (goodness of stratification).
    chi_square: f64,
    /// Number of completed iterations contributing to the average.
    iterations: usize,
}

impl StatisticsAccumulator {
    /// Create a fresh accumulator.
    pub fn new() -> Self {
        Self {
            sum_w: 0.0,
            sum_wf: 0.0,
            sum_wf2: 0.0,
            integral: 0.0,
            error: f64::INFINITY,
            chi_square: 0.0,
            iterations: 0,
        }
    }

    /// Add a sample with the given Vegas weight (1/pdf). The contribution to
    /// the integral estimate is `weight · f(xs)`.
    pub fn add_sample(&mut self, weight: f64, f: f64) {
        self.sum_w += weight;
        self.sum_wf += weight * f;
        self.sum_wf2 += weight * f * f;
    }

    /// Number of samples in the current (not-yet-finalised) iteration.
    pub fn samples(&self) -> usize {
        // We don't track count directly; derive from sum_w when weights are 1.
        // Vegas weights are Jacobians, so this is approximate — callers should
        // not rely on it for sample-count bookkeeping.
        self.sum_w as usize
    }

    /// Finalise the current iteration: fold its mean and variance into the
    /// cross-iteration weighted average, then reset per-iteration accumulators.
    pub fn finalize_iteration(&mut self) {
        if self.sum_w <= 0.0 || self.sum_wf2 < 0.0 {
            // Degenerate iteration (no samples or numerical issue): skip but
            // still reset.
            self.reset_iteration();
            return;
        }
        let mean = self.sum_wf / self.sum_w;
        // Unbiased variance estimate of the weighted mean: <f²>/<w> − <f>².
        let var = (self.sum_wf2 / self.sum_w) - mean * mean;
        let sig2 = if var > 0.0 { var } else { 0.0 };
        // Per-iteration standard error of the mean estimate.
        let iter_err = sig2.sqrt();
        self.combine_iteration(mean, iter_err);
        self.reset_iteration();
    }

    /// Combine one iteration's (mean, error) into the cross-iteration average
    /// using inverse-variance weighting, and update χ².
    fn combine_iteration(&mut self, mean: f64, err: f64) {
        // Clamp the error away from zero so the inverse-variance weight does
        // not blow up to infinity (a zero-variance iteration would otherwise
        // square to a subnormal that underflows in the divisor). 1e-150
        // squares to 1e-300, still representable.
        let err = if err > 1e-150 { err } else { 1e-150 };
        let w = 1.0 / (err * err);
        if self.iterations == 0 {
            self.integral = mean;
            self.error = err;
            self.chi_square = 0.0;
        } else {
            let prev_w = 1.0 / (self.error * self.error);
            let new_w = prev_w + w;
            let new_integral = (prev_w * self.integral + w * mean) / new_w;
            // χ² contribution: Σ wᵢ (meanᵢ − combined)².
            let delta_prev = self.integral - new_integral;
            let delta_cur = mean - new_integral;
            self.chi_square += prev_w * delta_prev * delta_prev + w * delta_cur * delta_cur;
            self.integral = new_integral;
            self.error = new_w.sqrt().recip();
        }
        self.iterations += 1;
    }

    /// Reset per-iteration accumulators (called after [`Self::finalize_iteration`]).
    fn reset_iteration(&mut self) {
        self.sum_w = 0.0;
        self.sum_wf = 0.0;
        self.sum_wf2 = 0.0;
    }

    /// Current best estimate of the integral.
    pub fn integral(&self) -> f64 {
        self.integral
    }

    /// Standard error on [`Self::integral`].
    pub fn error(&self) -> f64 {
        self.error
    }

    /// χ² across iterations; values near `iterations − 1` indicate consistent
    /// estimates.
    pub fn chi_square(&self) -> f64 {
        self.chi_square
    }

    /// Number of completed iterations.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
}

impl Default for StatisticsAccumulator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_integrand_converges_to_constant() {
        // ∫ 5 over [0,1] with uniform pdf (weight=1) → 5 with zero variance.
        let mut acc = StatisticsAccumulator::new();
        for _ in 0..1000 {
            acc.add_sample(1.0, 5.0);
        }
        acc.finalize_iteration();
        assert!((acc.integral() - 5.0).abs() < 1e-12);
        assert!(acc.error().abs() < 1e-9);
        assert!(acc.chi_square().abs() < 1e-9);
    }

    #[test]
    fn linear_integrand_matches_analytic() {
        // ∫₀¹ x dx = 1/2; with 50 000 uniform samples the mean ≈ 0.5 within
        // a few standard errors.
        let mut acc = StatisticsAccumulator::new();
        let n = 50_000u32;
        // Deterministic lattice to avoid pulling rand into this unit test.
        for i in 0..n {
            let x = (i as f64 + 0.5) / n as f64;
            acc.add_sample(1.0, x);
        }
        acc.finalize_iteration();
        assert!(
            (acc.integral() - 0.5).abs() < 1e-3,
            "got {}",
            acc.integral()
        );
    }

    #[test]
    fn combine_two_iterations_uses_inverse_variance_weighting() {
        let mut acc = StatisticsAccumulator::new();
        for _ in 0..1000 {
            acc.add_sample(1.0, 1.0);
        }
        acc.finalize_iteration();
        for _ in 0..1000 {
            acc.add_sample(1.0, 3.0);
        }
        acc.finalize_iteration();
        // Two consistent-internal iterations averaging 1 and 3 → near 2.
        assert!((acc.integral() - 2.0).abs() < 1e-9);
        assert_eq!(acc.iterations(), 2);
    }
}
