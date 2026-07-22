//! Python bindings for algebraic number fields and polynomials over them.
//!
//! Wraps [`ocas_domain::AlgebraicNumberField`] (an extension
//! $\mathbb{Q}(\alpha)$ defined by a monic minimal polynomial) and
//! [`ocas_poly::DenseUnivariatePolynomial`] over it, exposing construction,
//! arithmetic on elements, and Trager factorization.
//!
//! ```python
//! from ocas import AlgebraicExtension, AlgebraicPolynomial
//!
//! # ℚ(√2): minimal polynomial α² − 2 (ascending coefficients).
//! field = AlgebraicExtension([-2, 0, 1])
//! print(field.extension_degree())   # 2
//!
//! # x² − 2 splits over ℚ(√2) as (x − α)(x + α).
//! p = AlgebraicPolynomial(field, [-2, 0, 1])
//! for fac in p.factor():
//!     print(fac.factor.to_string(), fac.multiplicity)
//! ```

use ocas_domain::{AlgebraicElement, AlgebraicNumberField, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

/// An algebraic number field $K = \mathbb{Q}(\alpha)$ defined by a monic
/// minimal polynomial with rational coefficients.
///
/// The minimal polynomial is given as a list of coefficients in ascending
/// degree order; the leading (last) coefficient must be `1`. For example,
/// `AlgebraicExtension([-2, 0, 1])` defines $\alpha^2 - 2$ (i.e. $\mathbb{Q}(\sqrt{2})$).
#[pyclass(name = "AlgebraicExtension")]
pub struct PyAlgebraicExtension {
    pub(crate) field: AlgebraicNumberField,
}

/// An element of an [`AlgebraicExtension`], stored as a polynomial in $\alpha$
/// with rational coefficients (ascending degree).
#[pyclass(name = "AlgebraicElement")]
pub struct PyAlgebraicElement {
    pub(crate) elem: AlgebraicElement<Rational>,
}

/// A dense univariate polynomial over an [`AlgebraicExtension`].
#[pyclass(name = "AlgebraicPolynomial", skip_from_py_object)]
#[derive(Clone)]
pub struct PyAlgebraicPolynomial {
    pub(crate) inner: DenseUnivariatePolynomial<AlgebraicNumberField>,
}

/// A single (polynomial, multiplicity) factor returned by
/// [`AlgebraicPolynomial.factor`][PyAlgebraicPolynomial::factor].
#[pyclass(name = "AlgebraicFactor", skip_from_py_object)]
pub struct PyAlgebraicFactor {
    #[pyo3(get)]
    pub factor: PyAlgebraicPolynomial,
    #[pyo3(get)]
    pub multiplicity: usize,
}

// ------------------------------------------------------------------
//  Parsing helpers
// ------------------------------------------------------------------

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

/// Parse minimal-polynomial coefficients: a list of ints or `(num, denom)` tuples.
fn parse_min_poly(coeffs: &Bound<'_, PyAny>) -> PyResult<Vec<Rational>> {
    let iter = coeffs.try_iter().map_err(|_| {
        PyTypeError::new_err(
            "minimal polynomial coefficients must be a list of ints or (num, denom) tuples",
        )
    })?;
    iter.map(|c| py_to_rational(&c?)).collect()
}

/// Build an algebraic-field element from a Python value.
///
/// Accepts:
/// - an `int` or `(num, denom)` tuple (embedded as a base-domain constant),
/// - an [`AlgebraicElement`] instance, or
/// - a list of ints / `(num, denom)` tuples (ascending $\alpha$-polynomial).
fn py_to_anf_element(
    field: &AlgebraicNumberField,
    obj: &Bound<'_, PyAny>,
) -> PyResult<AlgebraicElement<Rational>> {
    if let Ok(r) = py_to_rational(obj) {
        return Ok(field.from_base(r));
    }
    if let Ok(elem) = obj.extract::<PyRef<'_, PyAlgebraicElement>>() {
        return Ok(elem.elem.clone());
    }
    if obj.try_iter().is_ok() {
        let cs: PyResult<Vec<Rational>> = obj
            .try_iter()
            .unwrap()
            .map(|c| py_to_rational(&c?))
            .collect();
        return Ok(field.element(cs?));
    }
    Err(PyTypeError::new_err(
        "coefficient must be int, (num, denom), list, or AlgebraicElement",
    ))
}

// ------------------------------------------------------------------
//  AlgebraicExtension
// ------------------------------------------------------------------

#[pymethods]
impl PyAlgebraicExtension {
    /// Create an algebraic number field from its monic minimal polynomial.
    ///
    /// `min_poly` is the list of rational coefficients in ascending degree
    /// order; the leading coefficient must be `1`.
    #[new]
    fn new(min_poly: &Bound<'_, PyAny>) -> PyResult<Self> {
        let coeffs = parse_min_poly(min_poly)?;
        if coeffs.len() < 2 {
            return Err(PyValueError::new_err(
                "minimal polynomial must have degree at least 1",
            ));
        }
        if coeffs.last() != Some(&Rational::new(1, 1)) {
            return Err(PyValueError::new_err("minimal polynomial must be monic"));
        }
        Ok(Self {
            field: AlgebraicNumberField::new(RationalDomain, coeffs),
        })
    }

    /// Return the extension degree $\deg(m)$.
    fn extension_degree(&self) -> usize {
        self.field.extension_degree()
    }

    /// Return the generator $\alpha$ of the extension.
    fn alpha(&self) -> PyAlgebraicElement {
        PyAlgebraicElement {
            elem: self.field.alpha(),
        }
    }

    /// Embed a rational constant (int or `(num, denom)`) into the field.
    #[allow(clippy::wrong_self_convention)]
    fn from_base(&self, c: &Bound<'_, PyAny>) -> PyResult<PyAlgebraicElement> {
        let r = py_to_rational(c)?;
        Ok(PyAlgebraicElement {
            elem: self.field.from_base(r),
        })
    }

    /// Create an element from $\alpha$-polynomial coefficients (ascending).
    fn element(&self, coeffs: &Bound<'_, PyAny>) -> PyResult<PyAlgebraicElement> {
        let iter = coeffs.try_iter().map_err(|_| {
            PyTypeError::new_err(
                "element coefficients must be a list of ints or (num, denom) tuples",
            )
        })?;
        let cs: PyResult<Vec<Rational>> = iter.map(|c| py_to_rational(&c?)).collect();
        Ok(PyAlgebraicElement {
            elem: self.field.element(cs?),
        })
    }

    fn __repr__(&self) -> String {
        format!("AlgebraicExtension(deg={})", self.field.extension_degree())
    }
}

// ------------------------------------------------------------------
//  AlgebraicElement
// ------------------------------------------------------------------

#[pymethods]
impl PyAlgebraicElement {
    /// Return the $\alpha$-polynomial coefficients (ascending) as decimal
    /// strings (rationals render as `n/d`).
    fn coeffs(&self) -> Vec<String> {
        self.elem.coeffs().iter().map(|c| c.to_string()).collect()
    }

    fn __str__(&self) -> String {
        format!("{}", self.elem)
    }

    fn __repr__(&self) -> String {
        format!("AlgebraicElement({})", self.elem)
    }
}

// ------------------------------------------------------------------
//  AlgebraicPolynomial
// ------------------------------------------------------------------

#[pymethods]
impl PyAlgebraicPolynomial {
    /// Create a polynomial over an algebraic number field.
    ///
    /// `coeffs` is a list with the constant term first. Each item is one of:
    /// - an `int` or `(num, denom)` tuple (a rational constant),
    /// - a list of ints / `(num, denom)` tuples (ascending $\alpha$-polynomial), or
    /// - an `AlgebraicElement`.
    #[new]
    fn new(field: PyRef<'_, PyAlgebraicExtension>, coeffs: &Bound<'_, PyAny>) -> PyResult<Self> {
        let f = &field.field;
        let iter = coeffs
            .try_iter()
            .map_err(|_| PyTypeError::new_err("polynomial coefficients must be a list"))?;
        let mut out = Vec::new();
        for c in iter {
            out.push(py_to_anf_element(f, &c?)?);
        }
        Ok(Self {
            inner: DenseUnivariatePolynomial::from_coeffs(f.clone(), out),
        })
    }

    /// Return the degree, or `None` for the zero polynomial.
    fn degree(&self) -> Option<usize> {
        self.inner.degree()
    }

    /// Return the number of stored coefficients.
    fn len(&self) -> usize {
        self.inner.coeffs().len()
    }

    /// Return `True` if this is the zero polynomial.
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    /// Return the coefficients (constant term first). Each coefficient is a
    /// list of decimal strings giving the $\alpha$-polynomial in ascending
    /// degree order.
    fn coeffs(&self) -> Vec<Vec<String>> {
        self.inner
            .coeffs()
            .iter()
            .map(|c| c.coeffs().iter().map(|r| r.to_string()).collect())
            .collect()
    }

    fn __str__(&self) -> String {
        format!("{}", PolyDisplay(&self.inner))
    }

    /// Return the list of irreducible factors with multiplicities.
    fn factor(&self) -> Vec<PyAlgebraicFactor> {
        self.inner
            .factor()
            .into_iter()
            .map(|(f, m)| PyAlgebraicFactor {
                factor: PyAlgebraicPolynomial { inner: f },
                multiplicity: m,
            })
            .collect()
    }
}

// ------------------------------------------------------------------
//  Display helper
// ------------------------------------------------------------------

/// Wrapper to render an algebraic-field polynomial as `c0 + c1*x + ...`.
struct PolyDisplay<'a>(&'a DenseUnivariatePolynomial<AlgebraicNumberField>);

impl std::fmt::Display for PolyDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coeffs = self.0.coeffs();
        if coeffs.is_empty() {
            return write!(f, "0");
        }
        let mut first = true;
        for (i, c) in coeffs.iter().enumerate() {
            // Skip zero coefficients for compactness.
            if c.coeffs().is_empty() {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            match i {
                0 => write!(f, "({})", c)?,
                1 => write!(f, "({})*x", c)?,
                _ => write!(f, "({})*x^{}", c, i)?,
            }
        }
        if first {
            write!(f, "0")?;
        }
        Ok(())
    }
}
