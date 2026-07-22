//! Python bindings for oCAS.
//!
//! The top-level Python module is named `ocas` (see `pyproject.toml`).
//!
//! # Quick start
//!
//! ```python
//! import ocas
//!
//! e = ocas.Expression("x^2 + 2*x + 1")
//! print(e.diff("x"))                 # derivative
//! print(e.simplify())
//!
//! print(ocas.solve_diophantine(3, 5, 1))
//! ```

use pyo3::prelude::*;

pub mod algebraic;
pub mod domain;
pub mod eval;
pub mod expression;
pub mod matrix;
pub mod polynomial;
pub mod solve;

pub use algebraic::{
    PyAlgebraicElement, PyAlgebraicExtension, PyAlgebraicFactor, PyAlgebraicPolynomial,
};
pub use domain::{PyFiniteField, PyIntegerDomain, PyRationalDomain};
pub use eval::PyExpressionEvaluator;
pub use expression::Expression;
pub use matrix::PyMatrix;
pub use polynomial::{PyPolynomial, PyPolynomialFactor};
pub use solve::{
    PyDiophantineSolution, py_solve_diophantine, py_solve_linear_integer, py_solve_linear_rational,
};

/// The oCAS Python module entry point.
///
/// The function name `ocas` determines the exported symbol `PyInit_ocas`,
/// which Python looks for when importing the module named `ocas`.
#[pymodule]
fn ocas(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_class::<Expression>()?;
    m.add_class::<PyExpressionEvaluator>()?;
    m.add_class::<PyDiophantineSolution>()?;
    m.add_class::<PyPolynomial>()?;
    m.add_class::<PyPolynomialFactor>()?;
    m.add_class::<PyAlgebraicExtension>()?;
    m.add_class::<PyAlgebraicElement>()?;
    m.add_class::<PyAlgebraicPolynomial>()?;
    m.add_class::<PyAlgebraicFactor>()?;
    m.add_class::<PyMatrix>()?;
    m.add_class::<PyIntegerDomain>()?;
    m.add_class::<PyRationalDomain>()?;
    m.add_class::<PyFiniteField>()?;

    m.add_function(wrap_pyfunction!(py_solve_linear_rational, m)?)?;
    m.add_function(wrap_pyfunction!(py_solve_linear_integer, m)?)?;
    m.add_function(wrap_pyfunction!(py_solve_diophantine, m)?)?;

    Ok(())
}
