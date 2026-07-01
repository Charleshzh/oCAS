//! Python `Domain` classes — coefficient-domain selectors.
//!
//! Provides [`PyIntegerDomain`], [`PyRationalDomain`], and [`PyFiniteField`]
//! wrapper classes. These mirror oCAS's Rust-side domain objects and are used
//! by the Python `Polynomial` and `Matrix` classes to select the coefficient
//! ring.

use num_bigint::BigInt;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Selector enum bridging the three supported Python domains to their Rust
/// counterparts. Used internally by `Polynomial` and `Matrix`.
pub(crate) enum DomainKind {
    Integer,
    Rational,
    FiniteField(BigInt),
}

impl DomainKind {
    /// Parse a Python object into a [`DomainKind`].
    ///
    /// Accepts either a domain-name string (`"integer"`, `"int"`, `"Z"`,
    /// `"rational"`, `"rat"`, `"Q"`) or a `FiniteField` instance.
    pub(crate) fn from_py(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(s) = obj.extract::<String>() {
            match s.as_str() {
                "integer" | "int" | "Z" => Ok(DomainKind::Integer),
                "rational" | "rat" | "Q" => Ok(DomainKind::Rational),
                other => Err(PyValueError::new_err(format!(
                    "unknown domain string: {other:?} (expected one of \
                     'integer'/'int'/'Z', 'rational'/'rat'/'Q', or a FiniteField)"
                ))),
            }
        } else if let Ok(fq) = obj.extract::<PyRef<'_, PyFiniteField>>() {
            Ok(DomainKind::FiniteField(fq.modulus.clone()))
        } else {
            Err(PyValueError::new_err(
                "domain must be a string ('integer'/'rational') or a FiniteField instance",
            ))
        }
    }
}

/// The integer domain ℤ.
///
/// ```python
/// from ocas import IntegerDomain
/// d = IntegerDomain()
/// print(repr(d))  # IntegerDomain()
/// ```
#[pyclass(name = "IntegerDomain", skip_from_py_object)]
#[derive(Clone)]
pub struct PyIntegerDomain;

#[pymethods]
impl PyIntegerDomain {
    #[new]
    fn new() -> Self {
        PyIntegerDomain
    }

    fn __repr__(&self) -> String {
        "IntegerDomain()".to_string()
    }
}

/// The rational number domain ℚ.
///
/// ```python
/// from ocas import RationalDomain
/// d = RationalDomain()
/// ```
#[pyclass(name = "RationalDomain", skip_from_py_object)]
#[derive(Clone)]
pub struct PyRationalDomain;

#[pymethods]
impl PyRationalDomain {
    #[new]
    fn new() -> Self {
        PyRationalDomain
    }

    fn __repr__(&self) -> String {
        "RationalDomain()".to_string()
    }
}

/// A finite field GF(p) for a prime modulus `p`.
///
/// ```python
/// from ocas import FiniteField
/// gf5 = FiniteField(5)
/// print(repr(gf5))  # FiniteField(5)
/// ```
#[pyclass(name = "FiniteField", from_py_object)]
#[derive(Clone)]
pub struct PyFiniteField {
    pub(crate) modulus: BigInt,
}

#[pymethods]
impl PyFiniteField {
    /// Create `GF(p)` with prime modulus `p` (an int ≥ 2).
    #[new]
    fn new(modulus: i64) -> PyResult<Self> {
        if modulus < 2 {
            return Err(PyValueError::new_err(format!(
                "finite-field modulus must be a prime ≥ 2, got {modulus}"
            )));
        }
        Ok(PyFiniteField {
            modulus: BigInt::from(modulus),
        })
    }

    /// Return the prime modulus as a decimal string.
    #[getter]
    fn modulus(&self) -> String {
        self.modulus.to_string()
    }

    fn __repr__(&self) -> String {
        format!("FiniteField({})", self.modulus)
    }
}
