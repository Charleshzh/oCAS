//! Python bindings for numerical integration (Vegas adaptive Monte Carlo).
//!
//! Wraps [`ocas_eval::numeric::Vegas`] and the [`integrate_1d`](ocas_eval::numeric::integrate_1d)
//! convenience entry point. Python callables are wrapped into Rust closures
//! so users can integrate arbitrary Python functions.
//!
//! ```python
//! import ocas
//!
//! # One-shot helper: integrate x over [0, 1].
//! r = ocas.integrate_1d(lambda x: x, 0.0, 1.0)
//! print(r.integral, r.error)  # ~0.5, small error
//!
//! # Multi-dimensional Vegas with explicit control.
//! v = ocas.Vegas(2, n_samples=20000, iterations=8, seed=1)
//! r = v.integrate(lambda xs: xs[0] * xs[1])
//! ```

use ocas_eval::numeric::{IntegrateResult, Integrator, Vegas, VegasOptions};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyFloat, PyList};

/// Result of a numerical integration: the estimate and its standard error.
///
/// Instances are returned by [`Vegas.integrate`][PyVegas.integrate] and
/// [`ocas.integrate_1d`]. The fields `integral` and `error` are also
/// accessible by index (`result[0]`, `result[1]`) and by unpacking
/// (`integral, error = result`).
#[pyclass(name = "IntegrateResult")]
pub struct PyIntegrateResult {
    /// Best estimate of the integral.
    #[pyo3(get)]
    pub integral: f64,
    /// Estimated standard error on `integral`.
    #[pyo3(get)]
    pub error: f64,
}

#[pymethods]
impl PyIntegrateResult {
    #[new]
    fn new(integral: f64, error: f64) -> Self {
        Self { integral, error }
    }

    fn __getitem__(&self, idx: usize) -> PyResult<f64> {
        match idx {
            0 => Ok(self.integral),
            1 => Ok(self.error),
            _ => Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "IntegrateResult index {idx} out of range (only 0, 1 valid)"
            ))),
        }
    }

    fn __len__(&self) -> usize {
        2
    }

    fn __repr__(&self) -> String {
        format!(
            "IntegrateResult(integral={:?}, error={:?})",
            self.integral, self.error
        )
    }
}

impl From<IntegrateResult> for PyIntegrateResult {
    fn from(r: IntegrateResult) -> Self {
        Self {
            integral: r.integral,
            error: r.error,
        }
    }
}

/// Parse a non-negative integer option from a Python kwarg, validating it.
fn parse_usize_opt(value: &Bound<'_, PyAny>, name: &str) -> PyResult<usize> {
    let n: usize = value
        .extract()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be a non-negative integer")))?;
    Ok(n)
}

fn parse_f64_opt(value: &Bound<'_, PyAny>, name: &str) -> PyResult<f64> {
    value
        .extract()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be a float")))
}

/// Build [`VegasOptions`] from Python kwargs. All keys are optional.
fn kwargs_to_opts(
    n_bins: Option<&Bound<'_, PyAny>>,
    n_samples: Option<&Bound<'_, PyAny>>,
    iterations: Option<&Bound<'_, PyAny>>,
    learning_rate: Option<&Bound<'_, PyAny>>,
    seed: Option<&Bound<'_, PyAny>>,
) -> PyResult<VegasOptions> {
    let mut opts = VegasOptions::default();
    if let Some(v) = n_bins {
        opts.n_bins = parse_usize_opt(v, "n_bins")?;
        if opts.n_bins == 0 {
            return Err(PyValueError::new_err("n_bins must be >= 1"));
        }
    }
    if let Some(v) = n_samples {
        opts.n_samples = parse_usize_opt(v, "n_samples")?;
    }
    if let Some(v) = iterations {
        opts.iterations = parse_usize_opt(v, "iterations")?;
    }
    if let Some(v) = learning_rate {
        opts.learning_rate = parse_f64_opt(v, "learning_rate")?;
        if opts.learning_rate.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return Err(PyValueError::new_err("learning_rate must be positive"));
        }
    }
    if let Some(v) = seed {
        opts.seed = v
            .extract()
            .map_err(|_| PyTypeError::new_err("seed must be an integer"))?;
    }
    Ok(opts)
}

/// Adaptive Monte Carlo integrator (Vegas) over the unit hypercube.
///
/// Construct with the number of dimensions and optional tuning knobs, then
/// call [`integrate`][PyVegas.integrate] with a Python callable taking a
/// list of `n_dims` floats in `[0, 1]` and returning a float.
#[pyclass(name = "Vegas")]
pub struct PyVegas {
    inner: Vegas,
}

#[pymethods]
impl PyVegas {
    /// Create a Vegas integrator for `n_dims` dimensions.
    ///
    /// All options are keyword-only and optional:
    ///
    /// - `n_bins` (default 64): bins per dimension.
    /// - `n_samples` (default 10000): samples per iteration.
    /// - `iterations` (default 10): adaptive iterations.
    /// - `learning_rate` (default 1.5): grid smoothing rate (1.0–2.0).
    /// - `seed` (default 0x0C45): RNG seed for reproducibility.
    #[new]
    #[pyo3(signature = (n_dims, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None))]
    fn new(
        n_dims: usize,
        n_bins: Option<&Bound<'_, PyAny>>,
        n_samples: Option<&Bound<'_, PyAny>>,
        iterations: Option<&Bound<'_, PyAny>>,
        learning_rate: Option<&Bound<'_, PyAny>>,
        seed: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        if n_dims == 0 {
            return Err(PyValueError::new_err("n_dims must be >= 1"));
        }
        let opts = kwargs_to_opts(n_bins, n_samples, iterations, learning_rate, seed)?;
        Ok(Self {
            inner: Vegas::new(n_dims, opts),
        })
    }

    /// Integrate a Python callable. The callable receives a list of `n_dims`
    /// floats in `[0, 1]` and must return a float.
    fn integrate(&mut self, f: &Bound<'_, PyAny>) -> PyResult<PyIntegrateResult> {
        let r = Python::attach(|py| -> PyResult<IntegrateResult> {
            let cb = f.clone();
            let wrapped = |x: &[f64]| -> f64 {
                match PyList::new(py, x.iter().copied()) {
                    Ok(list) => {
                        let arg = list.into_any();
                        match cb.call1((arg,)) {
                            Ok(value) => match value.extract::<f64>() {
                                Ok(v) => v,
                                Err(e) => {
                                    e.restore(py);
                                    f64::NAN
                                }
                            },
                            Err(e) => {
                                e.restore(py);
                                f64::NAN
                            }
                        }
                    }
                    Err(e) => {
                        e.restore(py);
                        f64::NAN
                    }
                }
            };
            let result = self.inner.integrate(&wrapped);
            if PyErr::take(py).is_some() {
                return Err(PyValueError::new_err(
                    "integrand raised an exception (or returned non-float)",
                ));
            }
            Ok(result)
        })?;
        Ok(r.into())
    }

    /// Latest accumulated estimate and error after `integrate`.
    #[getter]
    fn result(&self) -> PyIntegrateResult {
        self.inner.result().into()
    }

    /// Number of completed iterations.
    #[getter]
    fn iterations(&self) -> usize {
        self.inner.iterations()
    }

    fn __repr__(&self) -> String {
        format!(
            "Vegas(iterations={}, integral={:?})",
            self.inner.iterations(),
            self.inner.result().integral
        )
    }
}

/// Numerically integrate a one-dimensional Python callable `f` over `[a, b]`.
///
/// All options are keyword-only and optional (see
/// [`Vegas`][PyVegas] for their meaning). Returns an
/// [`IntegrateResult`][PyIntegrateResult].
#[pyfunction]
#[pyo3(signature = (f, a, b, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None))]
#[allow(clippy::too_many_arguments)]
pub fn integrate_1d(
    f: &Bound<'_, PyAny>,
    a: f64,
    b: f64,
    n_bins: Option<&Bound<'_, PyAny>>,
    n_samples: Option<&Bound<'_, PyAny>>,
    iterations: Option<&Bound<'_, PyAny>>,
    learning_rate: Option<&Bound<'_, PyAny>>,
    seed: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyIntegrateResult> {
    if a.partial_cmp(&b) != Some(std::cmp::Ordering::Less) {
        return Err(PyValueError::new_err(
            "integration upper bound b must be > a",
        ));
    }
    let opts = kwargs_to_opts(n_bins, n_samples, iterations, learning_rate, seed)?;
    let r = Python::attach(|py| -> PyResult<IntegrateResult> {
        let cb = f.clone();
        let wrapped = |x: f64| -> f64 {
            let arg = PyFloat::new(py, x);
            match cb.call1((&arg,)) {
                Ok(value) => match value.extract::<f64>() {
                    Ok(v) => v,
                    Err(e) => {
                        e.restore(py);
                        f64::NAN
                    }
                },
                Err(e) => {
                    e.restore(py);
                    f64::NAN
                }
            }
        };
        let result = ocas_eval::numeric::integrate_1d(wrapped, a, b, opts);
        if PyErr::take(py).is_some() {
            return Err(PyValueError::new_err(
                "integrand raised an exception (or returned non-float)",
            ));
        }
        Ok(result)
    })?;
    Ok(r.into())
}
