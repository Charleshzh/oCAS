//! Python `Matrix` class — dense matrices over ℤ, ℚ, or GF(p).
//!
//! Wraps [`ocas_poly::matrix::Matrix`] with an enum-erasure strategy,
//! mirroring the approach used by [`crate::polynomial`].

use crate::domain::DomainKind;
use ocas_domain::{Domain, FiniteField, Integer, IntegerDomain, Rational, RationalDomain};
use ocas_poly::matrix::{Matrix, MatrixError};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

/// Type-erased matrix over one of the three supported domains.
pub(crate) enum MatrixErased {
    Int(Matrix<IntegerDomain>),
    Rat(Matrix<RationalDomain>),
    Fq(Matrix<FiniteField>),
}

/// A dense matrix.
///
/// `rows` is a list of rows, each a list of coefficients. The `domain`
/// argument selects the coefficient ring exactly as for `Polynomial`.
///
/// ```python
/// from ocas import Matrix
///
/// a = Matrix([[1, 2], [3, 4]])
/// print(a.determinant())   # -2
/// print((a @ a).rows())    # [[7, 10], [15, 22]]
/// ```
#[pyclass(name = "Matrix", skip_from_py_object)]
pub struct PyMatrix {
    pub(crate) inner: MatrixErased,
}

fn map_matrix_err(e: MatrixError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Extract a 2-D integer coefficient matrix from a Python list-of-lists.
fn extract_int_rows(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Vec<Integer>>> {
    let data: Vec<Vec<i64>> = obj
        .extract()
        .map_err(|_| PyTypeError::new_err("integer matrix entries must be ints"))?;
    Ok(data
        .into_iter()
        .map(|r| r.into_iter().map(Integer::from).collect())
        .collect())
}

/// Extract a 2-D rational coefficient matrix. Each entry is an int or a
/// `(num, denom)` tuple; the whole matrix must be uniformly one form.
fn extract_rat_rows(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Vec<Rational>>> {
    if let Ok(int_rows) = obj.extract::<Vec<Vec<i64>>>() {
        Ok(int_rows
            .into_iter()
            .map(|r| r.into_iter().map(|n| Rational::new(n, 1)).collect())
            .collect())
    } else {
        let pair_rows: Vec<Vec<(i64, i64)>> = obj.extract().map_err(|_| {
            PyTypeError::new_err("rational entries must be ints or (num, denom) tuples")
        })?;
        pair_rows
            .into_iter()
            .map(|r| {
                r.into_iter()
                    .map(|(num, den)| {
                        if den == 0 {
                            Err(PyValueError::new_err("rational denominator cannot be zero"))
                        } else {
                            Ok(Rational::new(num, den))
                        }
                    })
                    .collect()
            })
            .collect()
    }
}

/// Build a `PyMatrix` from a Python list-of-lists and a domain kind.
pub(crate) fn build_matrix(rows: &Bound<'_, PyAny>, domain: &DomainKind) -> PyResult<PyMatrix> {
    let inner = match domain {
        DomainKind::Integer => {
            let r = extract_int_rows(rows)?;
            MatrixErased::Int(Matrix::from_rows(r, IntegerDomain))
        }
        DomainKind::Rational => {
            let r = extract_rat_rows(rows)?;
            MatrixErased::Rat(Matrix::from_rows(r, RationalDomain))
        }
        DomainKind::FiniteField(p) => {
            let field = FiniteField::new(p.clone());
            let data: Vec<Vec<i64>> = rows
                .extract()
                .map_err(|_| PyTypeError::new_err("finite-field matrix entries must be ints"))?;
            let rows: Vec<Vec<_>> = data
                .into_iter()
                .map(|r| r.into_iter().map(|v| field.element(v)).collect())
                .collect();
            MatrixErased::Fq(Matrix::from_rows(rows, field))
        }
    };
    Ok(PyMatrix { inner })
}

#[pymethods]
impl PyMatrix {
    /// Create a matrix from a list of rows.
    ///
    /// `domain` selects the coefficient ring (`"integer"` is the default).
    #[new]
    #[pyo3(signature = (rows, domain=None))]
    fn new(rows: &Bound<'_, PyAny>, domain: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let kind = match domain {
            Some(d) => DomainKind::from_py(d)?,
            None => DomainKind::Integer,
        };
        build_matrix(rows, &kind)
    }

    /// Number of rows.
    #[getter]
    fn nrows(&self) -> usize {
        match &self.inner {
            MatrixErased::Int(m) => m.nrows(),
            MatrixErased::Rat(m) => m.nrows(),
            MatrixErased::Fq(m) => m.nrows(),
        }
    }

    /// Number of columns.
    #[getter]
    fn ncols(&self) -> usize {
        match &self.inner {
            MatrixErased::Int(m) => m.ncols(),
            MatrixErased::Rat(m) => m.ncols(),
            MatrixErased::Fq(m) => m.ncols(),
        }
    }

    /// Return the shape as `(nrows, ncols)`.
    fn shape(&self) -> (usize, usize) {
        (self.nrows(), self.ncols())
    }

    /// Return `self[(i, j)]`.
    fn __getitem__(&self, idx: (usize, usize)) -> PyResult<String> {
        let (i, j) = idx;
        let out = match &self.inner {
            MatrixErased::Int(m) => {
                if i >= m.nrows() || j >= m.ncols() {
                    return Err(PyValueError::new_err("index out of bounds"));
                }
                m[(i, j)].to_string()
            }
            MatrixErased::Rat(m) => {
                if i >= m.nrows() || j >= m.ncols() {
                    return Err(PyValueError::new_err("index out of bounds"));
                }
                m[(i, j)].to_string()
            }
            MatrixErased::Fq(m) => {
                if i >= m.nrows() || j >= m.ncols() {
                    return Err(PyValueError::new_err("index out of bounds"));
                }
                m[(i, j)].value().to_string()
            }
        };
        Ok(out)
    }

    /// Return all rows as a list of lists of decimal strings (rational
    /// entries are rendered as `n/d`). Wrap each entry in `int(...)` to
    /// obtain Python integers.
    fn rows(&self) -> Vec<Vec<String>> {
        match &self.inner {
            MatrixErased::Int(m) => (0..m.nrows())
                .map(|i| (0..m.ncols()).map(|j| m[(i, j)].to_string()).collect())
                .collect(),
            MatrixErased::Rat(m) => (0..m.nrows())
                .map(|i| (0..m.ncols()).map(|j| m[(i, j)].to_string()).collect())
                .collect(),
            MatrixErased::Fq(m) => (0..m.nrows())
                .map(|i| {
                    (0..m.ncols())
                        .map(|j| m[(i, j)].value().to_string())
                        .collect()
                })
                .collect(),
        }
    }

    /// Return the transpose.
    fn transpose(&self) -> PyMatrix {
        match &self.inner {
            MatrixErased::Int(m) => PyMatrix {
                inner: MatrixErased::Int(m.transpose()),
            },
            MatrixErased::Rat(m) => PyMatrix {
                inner: MatrixErased::Rat(m.transpose()),
            },
            MatrixErased::Fq(m) => PyMatrix {
                inner: MatrixErased::Fq(m.transpose()),
            },
        }
    }

    /// Return the trace (sum of the diagonal) of a square matrix.
    fn trace(&self) -> PyResult<String> {
        match &self.inner {
            MatrixErased::Int(m) => Ok(m.trace().map_err(map_matrix_err)?.to_string()),
            MatrixErased::Rat(m) => Ok(m.trace().map_err(map_matrix_err)?.to_string()),
            MatrixErased::Fq(m) => Ok(m.trace().map_err(map_matrix_err)?.value().to_string()),
        }
    }

    /// Return the rank.
    fn rank(&self) -> usize {
        match &self.inner {
            MatrixErased::Int(m) => m.rank(),
            MatrixErased::Rat(m) => m.rank(),
            MatrixErased::Fq(m) => m.rank(),
        }
    }

    /// Return the determinant of a square matrix.
    fn determinant(&self) -> PyResult<String> {
        match &self.inner {
            MatrixErased::Int(m) => Ok(m.determinant().map_err(map_matrix_err)?.to_string()),
            MatrixErased::Rat(m) => Ok(m.determinant().map_err(map_matrix_err)?.to_string()),
            MatrixErased::Fq(m) => Ok(m.determinant().map_err(map_matrix_err)?.value().to_string()),
        }
    }

    /// Return the inverse, or raise `ValueError` if singular/non-square.
    fn inverse(&self) -> PyResult<PyMatrix> {
        match &self.inner {
            MatrixErased::Int(m) => Ok(PyMatrix {
                inner: MatrixErased::Int(m.inverse().map_err(map_matrix_err)?),
            }),
            MatrixErased::Rat(m) => Ok(PyMatrix {
                inner: MatrixErased::Rat(m.inverse().map_err(map_matrix_err)?),
            }),
            MatrixErased::Fq(m) => Ok(PyMatrix {
                inner: MatrixErased::Fq(m.inverse().map_err(map_matrix_err)?),
            }),
        }
    }

    /// Solve `self * x = rhs` for the vector `rhs`.
    ///
    /// `rhs` is a list of ints (integer/fq) or `(num, denom)` tuples
    /// (rational). Returns the solution as a list of decimal strings.
    fn solve(&self, rhs: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
        match &self.inner {
            MatrixErased::Int(m) => {
                let b: Vec<Integer> = extract_int_vector(rhs)?;
                let sol = m.solve(&b).map_err(map_matrix_err)?;
                Ok(sol.into_iter().map(|c| c.to_string()).collect())
            }
            MatrixErased::Rat(m) => {
                let b: Vec<Rational> = extract_rat_vector(rhs)?;
                let sol = m.solve(&b).map_err(map_matrix_err)?;
                Ok(sol.into_iter().map(|c| c.to_string()).collect())
            }
            MatrixErased::Fq(m) => {
                let field = m.domain().clone();
                let ints: Vec<i64> = rhs
                    .extract()
                    .map_err(|_| PyTypeError::new_err("finite-field rhs entries must be ints"))?;
                let b: Vec<_> = ints.into_iter().map(|v| field.element(v)).collect();
                let sol = m.solve(&b).map_err(map_matrix_err)?;
                Ok(sol.into_iter().map(|c| c.value().to_string()).collect())
            }
        }
    }

    /// Matrix product `self @ other`.
    fn __matmul__(&self, other: &PyMatrix) -> PyResult<PyMatrix> {
        match (&self.inner, &other.inner) {
            (MatrixErased::Int(a), MatrixErased::Int(b)) => Ok(PyMatrix {
                inner: MatrixErased::Int(a.matmul(b).map_err(map_matrix_err)?),
            }),
            (MatrixErased::Rat(a), MatrixErased::Rat(b)) => Ok(PyMatrix {
                inner: MatrixErased::Rat(a.matmul(b).map_err(map_matrix_err)?),
            }),
            (MatrixErased::Fq(a), MatrixErased::Fq(b)) => Ok(PyMatrix {
                inner: MatrixErased::Fq(a.matmul(b).map_err(map_matrix_err)?),
            }),
            _ => Err(PyTypeError::new_err(
                "@ requires both matrices to share the same coefficient domain",
            )),
        }
    }

    /// Element-wise addition `self + other`.
    fn __add__(&self, other: &PyMatrix) -> PyResult<PyMatrix> {
        match (&self.inner, &other.inner) {
            (MatrixErased::Int(a), MatrixErased::Int(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = *a.domain();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.add(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Int(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            (MatrixErased::Rat(a), MatrixErased::Rat(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = *a.domain();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.add(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Rat(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            (MatrixErased::Fq(a), MatrixErased::Fq(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = a.domain().clone();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.add(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Fq(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            _ => Err(PyTypeError::new_err(
                "+ requires both matrices to share the same coefficient domain",
            )),
        }
    }

    /// Element-wise subtraction `self - other`.
    fn __sub__(&self, other: &PyMatrix) -> PyResult<PyMatrix> {
        match (&self.inner, &other.inner) {
            (MatrixErased::Int(a), MatrixErased::Int(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = *a.domain();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.sub(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Int(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            (MatrixErased::Rat(a), MatrixErased::Rat(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = *a.domain();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.sub(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Rat(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            (MatrixErased::Fq(a), MatrixErased::Fq(b)) => {
                check_shape(a.nrows(), a.ncols(), b.nrows(), b.ncols())?;
                let d = a.domain().clone();
                let mut data = Vec::with_capacity(a.nrows() * a.ncols());
                for i in 0..a.nrows() {
                    for j in 0..a.ncols() {
                        data.push(d.sub(&a[(i, j)], &b[(i, j)]));
                    }
                }
                Ok(PyMatrix {
                    inner: MatrixErased::Fq(Matrix::new(a.nrows(), a.ncols(), data, d)),
                })
            }
            _ => Err(PyTypeError::new_err(
                "- requires both matrices to share the same coefficient domain",
            )),
        }
    }

    fn __eq__(&self, other: &PyMatrix) -> bool {
        match (&self.inner, &other.inner) {
            (MatrixErased::Int(a), MatrixErased::Int(b)) => a == b,
            (MatrixErased::Rat(a), MatrixErased::Rat(b)) => a == b,
            (MatrixErased::Fq(a), MatrixErased::Fq(b)) => a == b,
            _ => false,
        }
    }

    fn __repr__(&self) -> String {
        let dom = match &self.inner {
            MatrixErased::Int(_) => "integer",
            MatrixErased::Rat(_) => "rational",
            MatrixErased::Fq(_) => "finite-field",
        };
        format!(
            "Matrix({}x{}, domain='{}')",
            self.nrows(),
            self.ncols(),
            dom
        )
    }
}

fn check_shape(r1: usize, c1: usize, r2: usize, c2: usize) -> PyResult<()> {
    if r1 != r2 || c1 != c2 {
        Err(PyValueError::new_err(format!(
            "shape mismatch: {r1}x{c1} vs {r2}x{c2}"
        )))
    } else {
        Ok(())
    }
}

/// Extract an integer vector from a Python iterable of ints.
fn extract_int_vector(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Integer>> {
    let ints: Vec<i64> = obj
        .extract()
        .map_err(|_| PyTypeError::new_err("integer rhs entries must be ints"))?;
    Ok(ints.into_iter().map(Integer::from).collect())
}

/// Extract a rational vector from a Python iterable. Entries may be ints or
/// `(num, denom)` tuples, uniformly.
fn extract_rat_vector(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Rational>> {
    if let Ok(ints) = obj.extract::<Vec<i64>>() {
        Ok(ints.into_iter().map(|n| Rational::new(n, 1)).collect())
    } else {
        let pairs: Vec<(i64, i64)> = obj.extract().map_err(|_| {
            PyTypeError::new_err("rational rhs entries must be ints or (num, denom) tuples")
        })?;
        pairs
            .into_iter()
            .map(|(num, den)| {
                if den == 0 {
                    Err(PyValueError::new_err("rational denominator cannot be zero"))
                } else {
                    Ok(Rational::new(num, den))
                }
            })
            .collect()
    }
}
