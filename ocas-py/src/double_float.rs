//! Python bindings for [`DoubleF64`].

use ocas_domain::DoubleF64;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Python wrapper for double-precision floating-point arithmetic.
///
/// Provides ~31 decimal digits of precision (~84 binary bits) using
/// Dekker/Knuth "double-float" arithmetic.
///
/// ```python
/// from ocas import DoubleF64
///
/// a = DoubleF64(1.0)
/// b = DoubleF64(2.0)
/// print(a + b)  # 3.0
/// ```
#[pyclass(name = "DoubleF64", skip_from_py_object)]
#[derive(Debug, Clone, Copy)]
pub struct PyDoubleF64 {
    pub(crate) inner: DoubleF64,
}

#[pymethods]
impl PyDoubleF64 {
    /// Create a new `DoubleF64` from a high and optional low component.
    ///
    /// ```python
    /// x = DoubleF64(3.14)       # from single float
    /// y = DoubleF64(1.0, 1e-20) # from hi, lo pair
    /// ```
    #[new]
    #[pyo3(signature = (hi, lo = 0.0))]
    fn new(hi: f64, lo: f64) -> Self {
        Self {
            inner: DoubleF64::new(hi, lo),
        }
    }

    /// Return the high-order component as a Python `float`.
    #[allow(clippy::wrong_self_convention)]
    fn to_f64(&self) -> f64 {
        self.inner.to_f64()
    }

    /// Return `(hi, lo)` as a Python tuple.
    fn components(&self) -> (f64, f64) {
        (self.inner.hi, self.inner.lo)
    }

    fn __repr__(&self) -> String {
        format!("DoubleF64({}, {})", self.inner.hi, self.inner.lo)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __add__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner + other.inner,
        }
    }

    fn __sub__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner - other.inner,
        }
    }

    fn __mul__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner * other.inner,
        }
    }

    fn __truediv__(&self, other: &Self) -> PyResult<Self> {
        if other.inner.hi == 0.0 && other.inner.lo == 0.0 {
            Err(PyValueError::new_err("division by zero"))
        } else {
            Ok(Self {
                inner: self.inner / other.inner,
            })
        }
    }

    fn __neg__(&self) -> Self {
        Self { inner: -self.inner }
    }

    fn __abs__(&self) -> Self {
        Self {
            inner: self.inner.dabs(),
        }
    }

    fn __pow__(&self, exp: i64, _mod: Option<i64>) -> Self {
        Self {
            inner: self.inner.powi(exp),
        }
    }

    fn __richcmp__(&self, other: &Self, op: pyo3::basic::CompareOp) -> bool {
        match op {
            pyo3::basic::CompareOp::Lt => self.inner < other.inner,
            pyo3::basic::CompareOp::Le => self.inner <= other.inner,
            pyo3::basic::CompareOp::Eq => self.inner == other.inner,
            pyo3::basic::CompareOp::Ne => self.inner != other.inner,
            pyo3::basic::CompareOp::Gt => self.inner > other.inner,
            pyo3::basic::CompareOp::Ge => self.inner >= other.inner,
        }
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hi.to_bits().hash(&mut hasher);
        self.inner.lo.to_bits().hash(&mut hasher);
        hasher.finish()
    }

    // Transcendental functions

    /// Sine of this value.
    fn sin(&self) -> Self {
        Self {
            inner: self.inner.sin(),
        }
    }

    /// Cosine of this value.
    fn cos(&self) -> Self {
        Self {
            inner: self.inner.cos(),
        }
    }

    /// Tangent of this value.
    fn tan(&self) -> Self {
        Self {
            inner: self.inner.tan(),
        }
    }

    /// Natural exponential (e^x).
    fn exp(&self) -> Self {
        Self {
            inner: self.inner.exp(),
        }
    }

    /// Natural logarithm.
    fn ln(&self) -> PyResult<Self> {
        if self.inner.hi <= 0.0 {
            Err(PyValueError::new_err("log of non-positive number"))
        } else {
            Ok(Self {
                inner: self.inner.ln(),
            })
        }
    }

    /// Square root.
    fn sqrt(&self) -> PyResult<Self> {
        if self.inner.hi < 0.0 {
            Err(PyValueError::new_err("sqrt of negative number"))
        } else {
            Ok(Self {
                inner: self.inner.sqrt(),
            })
        }
    }
}
