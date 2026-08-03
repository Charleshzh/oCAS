# Numeric Integration

oCAS ships an adaptive Monte Carlo integrator (Vegas) for numerical
evaluation of definite integrals. This is useful when a closed-form
symbolic antiderivative does not exist or is too expensive to compute.

---

## The Vegas Algorithm

Vegas (Lepage 1978) is an importance-sampling Monte Carlo method that
iteratively refines a piecewise-constant approximation of the integrand
to concentrate samples where the function varies most. oCAS implements
the standard adaptive-grid variant with configurable iteration count,
sample budget, and learning rate.

| Entry point | Description |
|---|---|
| `integrate_1d(f, a, b, opts)` | Convenience wrapper for 1-D integrals |
| `Vegas::new(n_dims, opts)` | Multi-dimensional integrator |
| `Integrator::integrate(&mut self, f)` | Run integration on a closure |

Results are returned as `IntegrateResult { integral, error }`.

---

## Quick Start: 1-D Integration

```rust
use ocas_eval::numeric::{integrate_1d, VegasOptions};

// Integrate sin(x) from 0 to π
let opts = VegasOptions {
    n_bins: 50,
    n_samples: 10_000,
    iterations: 5,
    learning_rate: 1.5,
    seed: 42,
};
let result = integrate_1d(|x| x.sin(), 0.0, std::f64::consts::PI, opts);
println!("integral = {:.6} ± {:.6}", result.integral, result.error);
// ≈ 2.000000
```

---

## Multi-Dimensional Integration

For integrals over multiple variables, use `Vegas` directly:

```rust
use ocas_eval::numeric::{Vegas, VegasOptions, Integrator};

let opts = VegasOptions {
    n_bins: 50,
    n_samples: 10_000,
    iterations: 10,
    learning_rate: 1.5,
    seed: 0,
};
let mut vegas = Vegas::new(2, opts);

// Integrate x*y over [0,1]×[0,1]  → exact result = 0.25
let result = vegas.integrate(&|coords| coords[0] * coords[1]);
println!("integral = {:.6} ± {:.6}", result.integral, result.error);
```

The closure receives a `&[f64]` slice with one entry per dimension, each
already mapped to the unit hypercube `[0, 1]ⁿ`.

---

## Tuning Parameters

| Field | Default | Effect |
|---|---|---|
| `n_bins` | 64 | Number of grid bins per dimension |
| `n_samples` | 10,000 | Samples per iteration |
| `iterations` | 10 | Number of adaptive iterations |
| `learning_rate` | 1.5 | Grid adaptation speed (1.0–2.0) |
| `seed` | `0x0C45` (`u64`) | RNG seed for reproducibility |

More iterations and samples reduce the variance estimate but increase
runtime. The learning rate controls how aggressively the grid adapts to
the integrand structure.

---

## Statistics Accumulator

`StatisticsAccumulator` is the internal inverse-variance weighted
accumulator used by Vegas. It tracks per-iteration integral estimates,
chi-square diagnostics, and the final combined result.

```rust
use ocas_eval::numeric::StatisticsAccumulator;

let mut acc = StatisticsAccumulator::new();
acc.add_sample(1.0, 1.5);
acc.add_sample(1.0, 1.8);
acc.add_sample(1.0, 2.1);
acc.finalize_iteration();

println!("integral = {:.6}", acc.integral());
println!("error    = {:.6}", acc.error());
println!("chi²     = {:.2}", acc.chi_square());
```

---

## Python & C Usage

### Python

```python
import ocas

# 1-D convenience
result = ocas.integrate_1d(lambda x: x**2, 0, 1, n_samples=10000, iterations=10)
print(result.integral, result.error)

# Multi-dimensional
vegas = ocas.Vegas(n_dims=2, n_samples=10000, iterations=10)
result = vegas.integrate(lambda coords: coords[0] * coords[1])
print(result.integral)
```

### C

```c
#include <ocas.h>

/* 1-D convenience (pass NULL for opts to use the library defaults) */
int err = 0;
struct ocas_OcasIntegrateResult result =
    ocas_integrate_1d(my_fn, NULL, 0.0, 1.0, NULL, &err);
printf("integral = %f ± %f\n", result.integral, result.error);
```

See the [Python API](./api/python.md) and [C/C++ API](./api/c.md)
chapters for full binding documentation.

---

## Limitations

- Vegas uses Monte Carlo sampling; results are statistical estimates with
  error bars, not exact symbolic values.
- Convergence is slow for high-dimensional integrals with strong
  correlations between variables.
- The integrand must be a plain `f64 → f64` function; symbolic expressions
  must be compiled to a numeric evaluator first (see [Evaluation & JIT](./evaluation.md)).
