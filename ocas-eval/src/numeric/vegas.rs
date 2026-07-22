//! Adaptive Monte Carlo integration (Vegas).
//!
//! Implements the classic Vegas algorithm of Lepage: a product grid over the
//! unit hypercube whose bin boundaries are iteratively refined so that each
//! bin captures an equal share of the integrand's variance. The estimate and
//! its error are combined across iterations using inverse-variance weighting
//! via [`StatisticsAccumulator`](super::statistics::StatisticsAccumulator).
//!
//! Integrands are closures `Fn(&[f64]) -> f64` taking a point in the unit
//! hypercube; map physical bounds with a linear change of variables in the
//! closure.

use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::numeric::statistics::StatisticsAccumulator;

/// Result of a numerical integration: the estimate and its standard error.
#[derive(Debug, Clone, Copy)]
pub struct IntegrateResult {
    /// Best estimate of the integral.
    pub integral: f64,
    /// Estimated standard error on `integral`.
    pub error: f64,
}

/// A numerical integrator produces an [`IntegrateResult`] from a closure.
pub trait Integrator {
    /// Integrate `f` (which receives a point in the unit hypercube of the
    /// integrator's native domain) and return the estimate and error.
    fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult;
}

/// Tuning knobs for [`Vegas`].
#[derive(Debug, Clone, Copy)]
pub struct VegasOptions {
    /// Number of bins per dimension (the grid is a product of 1-D grids).
    pub n_bins: usize,
    /// Number of samples per iteration.
    pub n_samples: usize,
    /// Number of adaptive iterations.
    pub iterations: usize,
    /// Grid smoothing / learning rate (typical 1.0–2.0).
    pub learning_rate: f64,
    /// RNG seed (deterministic across runs).
    pub seed: u64,
}

impl Default for VegasOptions {
    fn default() -> Self {
        Self {
            n_bins: 64,
            n_samples: 10_000,
            iterations: 10,
            learning_rate: 1.5,
            seed: 0x0C45,
        }
    }
}

/// One-dimensional Vegas grid: bin boundaries on [0,1] plus a per-bin
/// accumulator for the current iteration's importance estimate.
#[derive(Debug, Clone)]
struct GridAxis {
    /// Bin boundary positions, `n_bins + 1` of them, in `[0,1]`. Starts uniform.
    boundaries: Vec<f64>,
    /// Per-bin accumulator of `f²·w` (the importance training signal).
    bin_accum: Vec<f64>,
}

impl GridAxis {
    fn new(n_bins: usize) -> Self {
        let boundaries = (0..=n_bins).map(|i| i as f64 / n_bins as f64).collect();
        Self {
            boundaries,
            bin_accum: vec![0.0; n_bins],
        }
    }

    /// Sample a coordinate: pick a bin uniformly and map a uniform deviate
    /// inside it. Returns `(x, jacobian)` where `jacobian = n_bins · bin_width`
    /// is the inverse-pdf contribution from this axis.
    fn sample<R: Rng>(&self, rng: &mut R) -> (f64, f64) {
        let n = self.bin_accum.len();
        let b = rng.random_range(0..n);
        let lo = self.boundaries[b];
        let hi = self.boundaries[b + 1];
        let u = rng.random::<f64>();
        let x = lo + (hi - lo) * u;
        (x, (hi - lo) * n as f64)
    }

    /// Find the bin containing `x` and add `weight · f²` to that bin.
    fn add_training(&mut self, x: f64, weight: f64, f2: f64) {
        // Binary-search the bin boundaries for x.
        let b = match self
            .boundaries
            .binary_search_by(|v| v.partial_cmp(&x).unwrap_or(std::cmp::Ordering::Equal))
        {
            Ok(i) => i.min(self.bin_accum.len().saturating_sub(1)),
            Err(i) => i
                .saturating_sub(1)
                .min(self.bin_accum.len().saturating_sub(1)),
        };
        if b < self.bin_accum.len() {
            self.bin_accum[b] += weight * f2;
        }
    }

    /// Refine bin boundaries from the accumulated importance, using the Vegas
    /// cumulative-arc-length update: bin boundaries are redistributed so that
    /// each bin carries an equal share of the smoothed importance. The
    /// `learning_rate` damps the update via `d^(1/lr)` (1.0 = full step).
    fn update(&mut self, learning_rate: f64) {
        let n = self.bin_accum.len();
        if n == 0 {
            return;
        }
        let total: f64 = self.bin_accum.iter().sum();
        if total <= 0.0 {
            return;
        }
        // Smoothed, average-normalised importance per bin (1.0 = average).
        let avg = total / n as f64;
        let mut d = vec![0.0; n];
        for (i, d_slot) in d.iter_mut().enumerate() {
            let prev = if i > 0 { self.bin_accum[i - 1] } else { 0.0 };
            let next = if i + 1 < n {
                self.bin_accum[i + 1]
            } else {
                0.0
            };
            let smooth = (prev + self.bin_accum[i] + next) / 3.0;
            *d_slot = smooth / avg;
        }
        // Damp the update (learning_rate > 1 softens the grid change).
        if (learning_rate - 1.0).abs() > 1e-12 {
            for v in d.iter_mut() {
                *v = v.max(1e-30).powf(1.0 / learning_rate);
            }
        }
        // Cumulative arc length; redistribute boundaries at equal spacing.
        let mut cum = vec![0.0; n + 1];
        for i in 0..n {
            cum[i + 1] = cum[i] + d[i];
        }
        let final_cum = cum[n];
        if final_cum <= 0.0 {
            return;
        }
        let mut new_boundaries = vec![0.0; n + 1];
        new_boundaries[0] = 0.0;
        new_boundaries[n] = 1.0;
        let mut j = 0;
        for (i, boundary) in new_boundaries.iter_mut().enumerate().take(n).skip(1) {
            let target = i as f64 / n as f64 * final_cum;
            while j < n && cum[j + 1] < target {
                j += 1;
            }
            let lo = cum[j];
            let hi = cum[j + 1];
            let frac = if hi > lo {
                (target - lo) / (hi - lo)
            } else {
                0.0
            };
            *boundary = (j as f64 + frac) / n as f64;
        }
        // Enforce monotone non-decreasing; clamp tiny regressions.
        for i in 1..=n {
            if new_boundaries[i] < new_boundaries[i - 1] {
                new_boundaries[i] = new_boundaries[i - 1];
            }
        }
        new_boundaries[n] = 1.0;
        self.boundaries = new_boundaries;
        self.bin_accum.fill(0.0);
    }
}

/// Adaptive Monte Carlo integrator (Vegas) over the unit hypercube.
pub struct Vegas {
    opts: VegasOptions,
    axes: Vec<GridAxis>,
    accumulator: StatisticsAccumulator,
}

impl Vegas {
    /// Create a Vegas integrator for `n_dims` dimensions with the given options.
    pub fn new(n_dims: usize, opts: VegasOptions) -> Self {
        let axes = (0..n_dims).map(|_| GridAxis::new(opts.n_bins)).collect();
        Self {
            opts,
            axes,
            accumulator: StatisticsAccumulator::new(),
        }
    }

    /// Latest accumulated estimate and error after [`Integrator::integrate`].
    pub fn result(&self) -> IntegrateResult {
        IntegrateResult {
            integral: self.accumulator.integral(),
            error: self.accumulator.error(),
        }
    }

    /// Number of completed iterations.
    pub fn iterations(&self) -> usize {
        self.accumulator.iterations()
    }
}

impl Integrator for Vegas {
    fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult {
        let n_dims = self.axes.len();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(self.opts.seed);
        for _ in 0..self.opts.iterations {
            for _ in 0..self.opts.n_samples {
                // Sample each axis; collect x and the total jacobian.
                let mut x = Vec::with_capacity(n_dims);
                let mut jac = 1.0;
                for axis in self.axes.iter_mut() {
                    let (xi, wi) = axis.sample(&mut rng);
                    x.push(xi);
                    jac *= wi;
                }
                let fx = f(&x);
                self.accumulator.add_sample(jac, fx);
                let f2 = fx * fx;
                for (i, xi) in x.iter().enumerate() {
                    self.axes[i].add_training(*xi, jac, f2);
                }
            }
            self.accumulator.finalize_iteration();
            for axis in self.axes.iter_mut() {
                axis.update(self.opts.learning_rate);
            }
        }
        self.result()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_constant_exactly() {
        // ∫₀¹ 7 dx = 7 with zero variance.
        let mut v = Vegas::new(1, VegasOptions::default());
        let r = v.integrate(&|_x: &[f64]| 7.0);
        assert!((r.integral - 7.0).abs() < 1e-9, "got {}", r.integral);
        assert!(r.error < 1e-6, "error {}", r.error);
    }

    #[test]
    fn integrates_linear_to_one_percent() {
        // ∫₀¹ x dx = 1/2.
        let opts = VegasOptions {
            n_samples: 20_000,
            iterations: 8,
            ..VegasOptions::default()
        };
        let mut v = Vegas::new(1, opts);
        let r = v.integrate(&|x: &[f64]| x[0]);
        assert!((r.integral - 0.5).abs() < 0.01, "got {}", r.integral);
    }

    #[test]
    fn integrates_gaussian_peak() {
        // ∫ exp(-50 (x-0.5)²) dx over [0,1] ≈ sqrt(π/50) ≈ 0.2507.
        // Vegas' adaptive grid should resolve the narrow peak.
        let opts = VegasOptions {
            n_bins: 128,
            n_samples: 20_000,
            iterations: 12,
            ..VegasOptions::default()
        };
        let mut v = Vegas::new(1, opts);
        let r = v.integrate(&|x: &[f64]| (-50.0 * (x[0] - 0.5).powi(2)).exp());
        let analytic = (std::f64::consts::PI / 50.0).sqrt();
        assert!(
            (r.integral - analytic).abs() < 0.02 * analytic,
            "got {}, expected {}",
            r.integral,
            analytic
        );
    }

    #[test]
    fn integrates_two_dimensional_product() {
        // ∫₀¹∫₀¹ x·y dx dy = 1/4.
        let opts = VegasOptions {
            n_samples: 20_000,
            iterations: 8,
            ..VegasOptions::default()
        };
        let mut v = Vegas::new(2, opts);
        let r = v.integrate(&|x: &[f64]| x[0] * x[1]);
        assert!((r.integral - 0.25).abs() < 0.01, "got {}", r.integral);
    }

    #[test]
    fn deterministic_across_runs_with_same_seed() {
        let opts = VegasOptions {
            n_samples: 5000,
            iterations: 4,
            seed: 42,
            ..VegasOptions::default()
        };
        let mut a = Vegas::new(1, opts);
        let ra = a.integrate(&|x: &[f64]| x[0] * x[0]);
        let mut b = Vegas::new(1, opts);
        let rb = b.integrate(&|x: &[f64]| x[0] * x[0]);
        assert_eq!(ra.integral, rb.integral);
        assert_eq!(ra.error, rb.error);
    }

    #[test]
    fn integrate_1d_over_physical_bounds() {
        use super::super::integrate_1d;
        // ∫₀² x dx = 2.
        let r = integrate_1d(|x| x, 0.0, 2.0, Default::default());
        assert!((r.integral - 2.0).abs() < 0.02, "got {}", r.integral);
        // ∫₁² x² dx = 7/3.
        let r2 = integrate_1d(|x| x * x, 1.0, 2.0, VegasOptions::default());
        assert!(
            (r2.integral - 7.0 / 3.0).abs() < 0.03,
            "got {}",
            r2.integral
        );
    }
}
