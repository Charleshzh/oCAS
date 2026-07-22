//! Python bindings for hyper-dual numbers (forward automatic differentiation).
//!
//! Wraps [`ocas_domain::dual`] restricted to [`Rational`](ocas_domain::Rational)
//! coefficients. Exposes [`DualShape`] (a first-order derivative layout) and
//! [`HyperDual`] (a value plus its partial derivatives). Arithmetic operators
//! propagate derivatives exactly via rational arithmetic.
//!
//! Only polynomial/rational arithmetic is supported (`+`, `-`, `*`, `/`,
//! unary `-`). Transcendental functions (sin/exp/log) are out of scope; use
//! repeated multiplication for integer powers.
//!
//! ```python
//! from ocas import DualShape, HyperDual
//!
//! shape = DualShape.first_order(2)   # track ∂/∂x₀ and ∂/∂x₁
//! # f(x, y) = x * y at point (3, 5).
//! x = HyperDual.variable(shape, 0, 3)
//! y = HyperDual.variable(shape, 1, 5)
//! f = x * y
//! print(f.value())         # 15
//! print(f.deriv(0))        # 5  (∂f/∂x = y)
//! print(f.deriv(1))        # 3  (∂f/∂y = x)
//! ```

use std::sync::Arc;

use ocas_domain::Rational;
use ocas_domain::dual::{DualShape, HyperDual, new_first_order};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

/// Parse a Python `int` or `(num, denom)` tuple into a [`Rational`].
fn py_to_rational(obj: &Bound<'_, PyAny>) -> PyResult<Rational> {
    if let Ok(n) = obj.extract::<i64>() {
        Ok(Rational::new(n, 1))
    } else if let Ok((num, den)) = obj.extract::<(i64, i64)>() {
        if den == 0 {
            Err(PyValueError::new_err("rational denominator cannot be zero"))
        } else {
            Ok(Rational::new(num, den))
        }
    } else {
        Err(PyTypeError::new_err("expected int or (num, denom) tuple"))
    }
}

/// Format a [`Rational`] as a Python-friendly string `"num/den"` (or `"num"`
/// when the denominator is 1).
fn rational_to_string(r: &Rational) -> String {
    let n = r.numer().to_i64().unwrap_or(0);
    let d = r.denom().to_i64().unwrap_or(0);
    if d == 1 {
        n.to_string()
    } else {
        format!("{n}/{d}")
    }
}

/// A first-order dual-number shape: declares how many differentiation
/// variables are tracked. Cheap to share; build once via
/// [`first_order`][PyDualShape::first_order] and reuse.
#[pyclass(name = "DualShape")]
pub struct PyDualShape {
    shape: Arc<DualShape>,
}

#[pymethods]
impl PyDualShape {
    /// Build a first-order shape tracking one derivative per variable for
    /// `n_vars` variables.
    #[staticmethod]
    fn first_order(n_vars: usize) -> PyResult<Self> {
        if n_vars == 0 {
            return Err(PyValueError::new_err("n_vars must be >= 1"));
        }
        Ok(PyDualShape {
            shape: new_first_order::<Rational>(n_vars),
        })
    }

    /// Number of differentiation variables.
    #[getter]
    fn n_vars(&self) -> usize {
        self.shape.n_vars()
    }

    /// Total number of components (value + derivative slots).
    #[getter]
    fn n_components(&self) -> usize {
        self.shape.n_components()
    }

    fn __repr__(&self) -> String {
        format!(
            "DualShape(n_vars={}, n_components={})",
            self.n_vars(),
            self.n_components()
        )
    }
}

/// A hyper-dual number over the rationals: a value plus its partial
/// derivatives with respect to the variables of a [`DualShape`].
///
/// Construct via [`variable`][PyHyperDual::variable] or
/// [`constant`][PyHyperDual::constant]; combine with arithmetic operators.
#[pyclass(name = "HyperDual")]
pub struct PyHyperDual {
    inner: HyperDual<Rational>,
    shape: Arc<DualShape>,
}

impl PyHyperDual {
    fn new(inner: HyperDual<Rational>) -> Self {
        let shape = inner.shape().clone();
        PyHyperDual { inner, shape }
    }
}

#[pymethods]
impl PyHyperDual {
    /// Create an independent variable `x_i = value` (derivative 1 w.r.t. `i`).
    #[staticmethod]
    fn variable(shape: &PyDualShape, i: usize, value: &Bound<'_, PyAny>) -> PyResult<Self> {
        if i >= shape.shape.n_vars() {
            return Err(PyValueError::new_err(format!(
                "variable index {i} out of range (n_vars = {})",
                shape.shape.n_vars()
            )));
        }
        let v = py_to_rational(value)?;
        Ok(PyHyperDual::new(HyperDual::variable(&shape.shape, i, v)))
    }

    /// Create a constant `value` (all derivatives zero).
    #[staticmethod]
    fn constant(shape: &PyDualShape, value: &Bound<'_, PyAny>) -> PyResult<Self> {
        let v = py_to_rational(value)?;
        Ok(PyHyperDual::new(HyperDual::constant(&shape.shape, v)))
    }

    /// The scalar value component, as a string (`"n"` or `"n/d"`).
    fn value(&self) -> String {
        rational_to_string(self.inner.value())
    }

    /// The derivative w.r.t. variable `i` as a string, or `None` if the
    /// shape has no first-order component for `i`.
    fn deriv(&self, i: usize) -> Option<String> {
        self.inner.deriv(i).map(rational_to_string)
    }

    /// The number of differentiation variables.
    #[getter]
    fn n_vars(&self) -> usize {
        self.shape.n_vars()
    }

    fn __repr__(&self) -> String {
        format!(
            "HyperDual(value={}, n_vars={})",
            self.value(),
            self.n_vars()
        )
    }

    fn __add__(&self, other: &PyHyperDual) -> PyResult<PyHyperDual> {
        if !Arc::ptr_eq(&self.shape, &other.shape) {
            return Err(PyValueError::new_err(
                "cannot add HyperDuals with different shapes",
            ));
        }
        Ok(PyHyperDual::new(self.inner.clone() + other.inner.clone()))
    }

    fn __sub__(&self, other: &PyHyperDual) -> PyResult<PyHyperDual> {
        if !Arc::ptr_eq(&self.shape, &other.shape) {
            return Err(PyValueError::new_err(
                "cannot subtract HyperDuals with different shapes",
            ));
        }
        Ok(PyHyperDual::new(self.inner.clone() - other.inner.clone()))
    }

    fn __mul__(&self, other: &PyHyperDual) -> PyResult<PyHyperDual> {
        if !Arc::ptr_eq(&self.shape, &other.shape) {
            return Err(PyValueError::new_err(
                "cannot multiply HyperDuals with different shapes",
            ));
        }
        Ok(PyHyperDual::new(self.inner.clone() * other.inner.clone()))
    }

    fn __truediv__(&self, other: &PyHyperDual) -> PyResult<PyHyperDual> {
        if !Arc::ptr_eq(&self.shape, &other.shape) {
            return Err(PyValueError::new_err(
                "cannot divide HyperDuals with different shapes",
            ));
        }
        // Div panics when the divisor's value component is zero; guard first.
        if other.inner.value() == &Rational::new(0, 1) {
            return Err(PyValueError::new_err(
                "division by zero (value component is zero)",
            ));
        }
        Ok(PyHyperDual::new(self.inner.clone() / other.inner.clone()))
    }

    fn __neg__(&self) -> PyHyperDual {
        PyHyperDual::new(-self.inner.clone())
    }
}
