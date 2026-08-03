# Rust API Reference: Numerical Integration

This chapter covers the numerical integration functionality of the `ocas-eval` crate, centered on the **Vegas adaptive Monte Carlo integrator** (Lepage 1978). All interfaces are exported through the `ocas_eval::numeric` module.

**Module overview**

```rust
use ocas_eval::numeric::{
    Vegas,                // adaptive Monte Carlo integrator
    VegasOptions,         // configuration parameters
    IntegrateResult,      // integration result
    Integrator,           // integrator trait
    integrate_1d,         // 1-D convenience function
    StatisticsAccumulator, // inverse-variance weighted statistics accumulator
};
```

---

## IntegrateResult

**Signature**:

```rust
#[derive(Debug, Clone, Copy)]
pub struct IntegrateResult {
    pub integral: f64,
    pub error: f64,
}
```

**Description**: The return result of numerical integration, containing the integral estimate and its standard error.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `integral` | `f64` | The best estimate of the integral |
| `error` | `f64` | The standard error estimate of `integral` |

**Example**:

```rust
use ocas_eval::numeric::integrate_1d;

let r = integrate_1d(|x| x * x, 0.0, 1.0, Default::default());
println!("∫₀¹ x² dx ≈ {:.6} ± {:.6}", r.integral, r.error);
// Output: ∫₀¹ x² dx ≈ 0.3334 ± 0.0003 (exact value 1/3)
```

---

## Integrator

**Signature**:

```rust
pub trait Integrator {
    fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult;
}
```

**Description**: The unified trait for numerical integrators. The integrand `f` receives a point in the unit hypercube $[0,1]^d$ and returns a scalar value. Integrating over a physical domain requires performing the linear change of variables manually inside the closure.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&F where F: Fn(&[f64]) -> f64` | The integrand; the input is a coordinate point in $[0,1]^d$ |

**Returns**: `IntegrateResult` — the integral estimate and its standard error.

**Note**: `integrate` takes `&F` rather than `F`, so the same closure can be called multiple times without reallocation.

---

## VegasOptions

**Signature**:

```rust
#[derive(Debug, Clone, Copy)]
pub struct VegasOptions {
    pub n_bins: usize,
    pub n_samples: usize,
    pub iterations: usize,
    pub learning_rate: f64,
    pub seed: u64,
}
```

**Description**: Tuning parameters for the Vegas integrator.

**Fields**:

| Field | Type | Default | Description |
|---|---|---|---|
| `n_bins` | `usize` | `64` | The number of bins per dimension (the overall grid is the direct product of per-dimension 1-D grids) |
| `n_samples` | `usize` | `10_000` | The number of samples per iteration |
| `iterations` | `usize` | `10` | The number of adaptive iterations |
| `learning_rate` | `f64` | `1.5` | Grid smoothing/learning rate (typical range 1.0–2.0; larger values make grid updates gentler) |
| `seed` | `u64` | `0x0C45` | Random seed (deterministic across runs) |

**Implementation details**:

- The `Default` implementation provides the defaults above.
- Role of the learning rate: during grid updates, each bin's importance value $d_i$ is transformed into $d_i^{1/\text{learning\_rate}}$. `learning_rate = 1.0` is the full step; larger values make grid changes smoother and more stable.
- Uses the `Xoshiro256PlusPlus` pseudo-random number generator (initialized via `seed`) for reproducible results.

**Example**:

```rust
use ocas_eval::numeric::VegasOptions;

// High-precision configuration
let opts = VegasOptions {
    n_bins: 100,
    n_samples: 100_000,
    iterations: 20,
    learning_rate: 1.5,
    seed: 42,
};

// Use the defaults
let default_opts = VegasOptions::default();
```

**See also**: [`Vegas::new`](#vegasnew)

---

## Vegas

**Signature**:

```rust
pub struct Vegas { /* private fields */ }
```

**Description**: An adaptive Monte Carlo integrator based on the Lepage Vegas algorithm, operating on the unit hypercube $[0,1]^d$.

**Design rationale**: Vegas maintains a product grid (independent 1-D bins per dimension) and refines it iteratively so that each bin captures an approximately equal share of the integrand's variance. Estimates from multiple iterations are combined with inverse-variance weighting (see [`StatisticsAccumulator`](#statisticsaccumulator)).

### Vegas::new

**Signature**:

```rust
pub fn new(n_dims: usize, opts: VegasOptions) -> Self
```

**Description**: Creates a Vegas integrator in `n_dims` dimensions.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n_dims` | `usize` | The number of integration dimensions (≥ 1) |
| `opts` | `VegasOptions` | Tuning parameters |

**Returns**: An initialized `Vegas` instance with uniform bins in every dimension.

**Example**:

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};

// Two-dimensional integrator with default options
let mut vegas = Vegas::new(2, VegasOptions::default());
```

---

### Vegas::integrate

**Signature**:

```rust
pub fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult
```

**Description**: Runs the integration. For each iteration: sample `n_samples` points on the unit hypercube, evaluate the integrand, accumulate statistics, then update the grid boundaries to improve sampling efficiency.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&F where F: Fn(&[f64]) -> f64` | The integrand; the input is a coordinate in $[0,1]^d$ |

**Returns**: `IntegrateResult` — the inverse-variance weighted average estimate over all iterations and its standard error.

**Algorithm steps** (per iteration):

1. **Sampling**: for each dimension, pick a bin uniformly and sample uniformly within the bin; the Jacobian is `n_bins × bin_width` (i.e., the inverse-pdf contribution).
2. **Evaluation**: compute $f(\mathbf{x})$ and accumulate it into the [`StatisticsAccumulator`](#statisticsaccumulator) weighted by the total Jacobian.
3. **Training signal**: accumulate $f^2 \cdot w$ into the training slot of the corresponding bin.
4. **Statistics finalization**: call `finalize_iteration()` to fold this iteration's estimate into the cross-iteration inverse-variance weighted average.
5. **Grid update**: for each dimension, perform a smooth cumulative-arc-length redistribution — 3-bin smoothing, mean normalization, damping $d^{1/\text{lr}}$, redistributing boundaries at equal-arc-length quantiles, and monotonicity clamping.

**Example**:

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};

// Compute ∫₀¹∫₀¹ sin(x·y) dx dy
let opts = VegasOptions {
    n_samples: 50_000,
    iterations: 15,
    ..Default::default()
};
let mut vegas = Vegas::new(2, opts);
let r = vegas.integrate(&|x: &[f64]| (x[0] * x[1]).sin());
println!("Result: {:.6} ± {:.6}", r.integral, r.error);
// Result: ≈ 0.2398 ± 0.0002
```

**Note**: The integrand receives coordinates in $[0,1]^d$. To integrate over a physical domain $[a,b]$, perform the linear change $x_{\text{phys}} = a + u \cdot (b-a)$ inside the closure and multiply by the Jacobian $(b-a)$. In one dimension, [`integrate_1d`](#integrate_1d) can be used directly.

---

### Vegas::result

**Signature**:

```rust
pub fn result(&self) -> IntegrateResult
```

**Description**: Returns the accumulated estimate after the most recent `integrate` call.

**Returns**: `IntegrateResult` — the current best integral estimate and its standard error.

---

### Vegas::iterations

**Signature**:

```rust
pub fn iterations(&self) -> usize
```

**Description**: Returns the number of completed iterations.

---

## integrate_1d

**Signature**:

```rust
pub fn integrate_1d<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: VegasOptions,
) -> IntegrateResult
```

**Description**: A convenience function for one-dimensional numerical integration. Integrates $f(x)$ over the physical interval $[a, b]$, automatically handling the change of variables $u \mapsto a + u(b-a)$ and folding the Jacobian $(b-a)$ into the result.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `F where F: Fn(f64) -> f64` | The integrand; receives the physical coordinate $x \in [a, b]$ directly |
| `a` | `f64` | The lower limit of integration |
| `b` | `f64` | The upper limit of integration (may be less than `a`, in which case the result is negative) |
| `opts` | `VegasOptions` | Vegas parameters (use `Default::default()` for the defaults) |

**Returns**: `IntegrateResult` — the integral estimate and its standard error.

**Implementation details**: creates a one-dimensional `Vegas` instance, wraps `f` as `|u| f(a + u[0] * width) * width`, and calls `integrate`.

**Example**:

```rust
use ocas_eval::numeric::{integrate_1d, VegasOptions};

// ∫₀¹ x dx = 1/2
let r = integrate_1d(|x| x, 0.0, 1.0, Default::default());
assert!((r.integral - 0.5).abs() < 0.01);

// ∫₀² x² dx = 8/3 ≈ 2.667
let r = integrate_1d(|x| x * x, 0.0, 2.0, Default::default());
assert!((r.integral - 8.0 / 3.0).abs() < 0.05);

// High-precision configuration
let opts = VegasOptions {
    n_samples: 100_000,
    iterations: 20,
    ..Default::default()
};
let r = integrate_1d(|x| (-x * x).exp(), -5.0, 5.0, opts);
// ∫₋₅⁵ e^{-x²} dx ≈ √π ≈ 1.7725
assert!((r.integral - std::f64::consts::PI.sqrt()).abs() < 0.01);
```

**See also**: [`Vegas::new`](#vegasnew), [`Vegas::integrate`](#vegasintegrate)

---

## StatisticsAccumulator

**Signature**:

```rust
#[derive(Debug, Clone)]
pub struct StatisticsAccumulator { /* private fields */ }
```

**Description**: An inverse-variance weighted statistics accumulator for combining Monte Carlo integration estimates across iterations. Implements the weighting scheme from Lepage's original Vegas paper.

**Design rationale**: Each iteration produces an integral estimate $I_i$ and a standard error $\sigma_i$. Iterations are combined with inverse-variance weights $w_i = 1/\sigma_i^2$:

$$
I_{\text{combined}} = \frac{\sum_i w_i I_i}{\sum_i w_i}, \quad \sigma_{\text{combined}} = \frac{1}{\sqrt{\sum_i w_i}}
$$

A $\chi^2$ statistic $\sum_i w_i (I_i - I_{\text{combined}})^2$ is also maintained to diagnose the quality of the grid stratification — $\chi^2 \approx \text{iterations} - 1$ indicates that the per-iteration estimates are consistent.

**Internal fields** (private):

| Field | Type | Description |
|---|---|---|
| `sum_w` | `f64` | Current iteration: $\sum \text{weight}$ |
| `sum_wf` | `f64` | Current iteration: $\sum \text{weight} \cdot f$ |
| `sum_wf2` | `f64` | Current iteration: $\sum \text{weight} \cdot f^2$ |
| `integral` | `f64` | The cross-iteration weighted average integral estimate |
| `error` | `f64` | The standard error of the integral estimate |
| `chi_square` | `f64` | The cross-iteration $\chi^2$ statistic |
| `iterations` | `usize` | The number of completed iterations |

### StatisticsAccumulator::new

**Signature**:

```rust
pub fn new() -> Self
```

**Description**: Creates an empty accumulator. Initial state: `integral = 0.0`, `error = INFINITY`, `iterations = 0`.

---

### StatisticsAccumulator::add_sample

**Signature**:

```rust
pub fn add_sample(&mut self, weight: f64, f: f64)
```

**Description**: Adds a sample point. In Vegas, `weight` is the Jacobian (inverse pdf) at the point and `f` is the integrand value. The sample's contribution to this iteration's integral estimate is `weight * f`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `weight` | `f64` | The sample weight (in Vegas, the Jacobian = the product of `n_bins × bin_width` over all dimensions) |
| `f` | `f64` | The value of the integrand at the point |

**Internal accumulation**: `sum_w += weight`, `sum_wf += weight * f`, `sum_wf2 += weight * f²`.

---

### StatisticsAccumulator::finalize_iteration

**Signature**:

```rust
pub fn finalize_iteration(&mut self)
```

**Description**: Finalizes the current iteration — computes this iteration's mean and variance, folds them into the cross-iteration average with inverse-variance weighting, then resets this iteration's accumulators.

**This iteration's mean**: $\bar{f} = \sum(w \cdot f) / \sum w$

**This iteration's variance**: $\text{Var} = \sum(w \cdot f^2) / \sum w - \bar{f}^2$

**Degenerate handling**: if `sum_w ≤ 0` or the variance is negative (numerical issues), the iteration is skipped but still reset.

**Error clamping**: the standard error is clamped to be no smaller than $10^{-150}$ (to avoid overflow in the inverse-variance weights).

---

### StatisticsAccumulator::samples

**Signature**:

```rust
pub fn samples(&self) -> usize
```

**Description**: Returns the approximate number of samples in the current (unfinalized) iteration.

**Note**: derived internally from `sum_w` (`sum_w as usize`); exact when the weights are 1. Vegas's Jacobian weights make this value approximate.

---

### StatisticsAccumulator::integral

**Signature**:

```rust
pub fn integral(&self) -> f64
```

**Description**: The current best integral estimate (the cross-iteration inverse-variance weighted average).

---

### StatisticsAccumulator::error

**Signature**:

```rust
pub fn error(&self) -> f64
```

**Description**: The standard error of the integral estimate, i.e., $\sigma_{\text{combined}} = 1/\sqrt{\sum w_i}$.

---

### StatisticsAccumulator::chi_square

**Signature**:

```rust
pub fn chi_square(&self) -> f64
```

**Description**: The cross-iteration $\chi^2$ statistic. Values close to `iterations - 1` indicate that the per-iteration estimates are consistent (good grid stratification); values significantly larger indicate structure in the integrand not captured by the grid.

---

### StatisticsAccumulator::iterations

**Signature**:

```rust
pub fn iterations(&self) -> usize
```

**Description**: The number of completed iterations.

---

**Full example**:

```rust
use ocas_eval::numeric::StatisticsAccumulator;

let mut acc = StatisticsAccumulator::new();

// Simulate two iterations with three samples each
// First iteration
acc.add_sample(1.0, 1.5);
acc.add_sample(1.0, 1.8);
acc.add_sample(1.0, 2.1);
acc.finalize_iteration();
assert_eq!(acc.iterations(), 1);

// Second iteration
acc.add_sample(1.0, 1.6);
acc.add_sample(1.0, 1.9);
acc.add_sample(1.0, 2.0);
acc.finalize_iteration();
assert_eq!(acc.iterations(), 2);

// Cross-iteration weighted average
println!("Integral: {:.4} ± {:.4}", acc.integral(), acc.error());
println!("χ² = {:.2}", acc.chi_square());
```

**See also**: [`Vegas::integrate`](#vegasintegrate) (uses this accumulator internally)

---

## Working with ExpressionEvaluator

The Vegas integrator can be combined with [`ExpressionEvaluator`](./rust-evaluation.md) to numerically integrate symbolic expressions:

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sin", &[ctx.var("x")]);

let evaluator: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut vegas = Vegas::new(1, VegasOptions::default());
let r = vegas.integrate(&|x: &[f64]| {
    evaluator.evaluate(&[x[0]]).unwrap()[0]
});
```

For one-dimensional cases, the `integrate_1d` wrapper can be used directly:

```rust
use ocas_eval::numeric::integrate_1d;
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sin", &[ctx.var("x")]);
let evaluator: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

let r = integrate_1d(
    |x| evaluator.evaluate(&[x]).unwrap()[0],
    0.0,
    std::f64::consts::PI,
    Default::default(),
);
// ∫₀^π sin(x) dx ≈ 2.0
```

**See also**: [Evaluation and JIT](./rust-evaluation.md), [Monte Carlo integration foundations](../math/monte-carlo-integration.md)
