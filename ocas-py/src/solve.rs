//! Python `solve` module — linear and Diophantine equation solvers.

use ocas_calc::solve::{
    DiophantineSolution, SolveError, solve_diophantine as rs_diophantine,
    solve_linear_integer as rs_linear_integer, solve_linear_rational as rs_linear_rational,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn map_solve_err(e: SolveError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Solve the linear system `a * x = b` over the rationals.
///
/// Each element of `a` is a row (list of ints). Returns a list of
/// `(numerator, denominator)` tuples, or raises `ValueError` if no
/// unique solution exists.
#[pyfunction]
#[pyo3(name = "solve_linear_rational")]
pub fn py_solve_linear_rational(a: Vec<Vec<i64>>, b: Vec<i64>) -> PyResult<Vec<(i64, i64)>> {
    rs_linear_rational(&a, &b).map_err(map_solve_err)
}

/// Solve the linear system `a * x = b` over the integers.
///
/// Returns a list of integer solutions, or raises `ValueError` if no
/// integer solution exists.
#[pyfunction]
#[pyo3(name = "solve_linear_integer")]
pub fn py_solve_linear_integer(a: Vec<Vec<i64>>, b: Vec<i64>) -> PyResult<Vec<i64>> {
    rs_linear_integer(&a, &b).map_err(map_solve_err)
}

/// A solution to the Diophantine equation `a*x + b*y = c`.
#[pyclass(name = "DiophantineSolution")]
pub struct PyDiophantineSolution {
    /// Particular solution `(x0, y0)`.
    #[pyo3(get)]
    pub particular: (i64, i64),
    /// Homogeneous direction `(tx, ty)`; general solution is
    /// `(x0 + k*tx, y0 + k*ty)` for any integer `k`.
    #[pyo3(get)]
    pub general: (i64, i64),
}

impl From<DiophantineSolution> for PyDiophantineSolution {
    fn from(s: DiophantineSolution) -> Self {
        PyDiophantineSolution {
            particular: s.particular,
            general: s.general,
        }
    }
}

/// Solve the linear Diophantine equation `a*x + b*y = c`.
///
/// Returns `None` if no integer solution exists.
#[pyfunction]
#[pyo3(name = "solve_diophantine")]
pub fn py_solve_diophantine(a: i64, b: i64, c: i64) -> Option<PyDiophantineSolution> {
    rs_diophantine(a, b, c).map(Into::into)
}
