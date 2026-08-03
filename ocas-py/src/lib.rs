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
pub mod double_float;
pub mod dual;
pub mod eval;
pub mod expression;
pub mod groebner;
pub mod matrix;
pub mod ntheory;
pub mod numeric;
pub mod ode;
pub mod polynomial;
pub mod solve;
pub mod tensor;

pub use algebraic::{
    PyAlgebraicElement, PyAlgebraicExtension, PyAlgebraicFactor, PyAlgebraicPolynomial,
};
pub use domain::{PyFiniteField, PyIntegerDomain, PyRationalDomain};
pub use double_float::PyDoubleF64;
pub use dual::{PyDualShape, PyHyperDual};
pub use eval::PyExpressionEvaluator;
pub use expression::Expression;
pub use groebner::{
    PyGroebnerBasis, PyHilbertSeries, PyMultivariatePolynomial, PyPolynomialSystemSolution,
    PyPrimaryComponent, PyRealSolution, py_eliminate, py_groebner_basis, py_hilbert_series,
    py_ideal_contains, py_ideal_radical, py_is_zero_dimensional, py_primary_decomposition,
    py_solve_polynomial_system,
};
pub use matrix::PyMatrix;
pub use ntheory::{
    py_crt, py_discrete_log, py_divisor_count, py_divisor_sigma, py_factorint, py_isprime,
    py_isprime_u64, py_jacobi_symbol, py_liouville_lambda, py_mobius, py_nextprime, py_totient,
};
pub use numeric::{PyIntegrateResult, PyVegas, integrate_1d};
pub use ode::{py_classify_ode, py_dsolve, py_dsolve_ivp};
pub use polynomial::{PyPolynomial, PyPolynomialFactor};
pub use solve::{
    PyDiophantineSolution, py_solve_diophantine, py_solve_linear_integer, py_solve_linear_rational,
};
pub use tensor::{
    PyTensor, canonicalize_tensors, contract_tensors, refresh_dummies, tensor_symmetrise_sign,
    young_project,
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
    m.add_class::<PyVegas>()?;
    m.add_class::<PyIntegrateResult>()?;
    m.add_class::<PyTensor>()?;
    m.add_class::<PyDualShape>()?;
    m.add_class::<PyHyperDual>()?;
    m.add_class::<PyDoubleF64>()?;
    m.add_class::<PyGroebnerBasis>()?;
    m.add_class::<PyRealSolution>()?;
    m.add_class::<PyPolynomialSystemSolution>()?;
    m.add_class::<PyHilbertSeries>()?;
    m.add_class::<PyPrimaryComponent>()?;
    m.add_class::<PyMultivariatePolynomial>()?;

    m.add_function(wrap_pyfunction!(py_solve_linear_rational, m)?)?;
    m.add_function(wrap_pyfunction!(py_solve_linear_integer, m)?)?;
    m.add_function(wrap_pyfunction!(py_solve_diophantine, m)?)?;
    m.add_function(wrap_pyfunction!(integrate_1d, m)?)?;
    m.add_function(wrap_pyfunction!(contract_tensors, m)?)?;
    m.add_function(wrap_pyfunction!(tensor_symmetrise_sign, m)?)?;
    m.add_function(wrap_pyfunction!(canonicalize_tensors, m)?)?;
    m.add_function(wrap_pyfunction!(young_project, m)?)?;
    m.add_function(wrap_pyfunction!(refresh_dummies, m)?)?;
    m.add_function(wrap_pyfunction!(py_classify_ode, m)?)?;
    m.add_function(wrap_pyfunction!(py_dsolve, m)?)?;
    m.add_function(wrap_pyfunction!(py_dsolve_ivp, m)?)?;
    m.add_function(wrap_pyfunction!(py_factorint, m)?)?;
    m.add_function(wrap_pyfunction!(py_isprime, m)?)?;
    m.add_function(wrap_pyfunction!(py_isprime_u64, m)?)?;
    m.add_function(wrap_pyfunction!(py_nextprime, m)?)?;
    m.add_function(wrap_pyfunction!(py_discrete_log, m)?)?;
    m.add_function(wrap_pyfunction!(py_crt, m)?)?;
    m.add_function(wrap_pyfunction!(py_jacobi_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(py_totient, m)?)?;
    m.add_function(wrap_pyfunction!(py_mobius, m)?)?;
    m.add_function(wrap_pyfunction!(py_divisor_count, m)?)?;
    m.add_function(wrap_pyfunction!(py_divisor_sigma, m)?)?;
    m.add_function(wrap_pyfunction!(py_liouville_lambda, m)?)?;
    m.add_function(wrap_pyfunction!(py_groebner_basis, m)?)?;
    m.add_function(wrap_pyfunction!(py_ideal_contains, m)?)?;
    m.add_function(wrap_pyfunction!(py_solve_polynomial_system, m)?)?;
    m.add_function(wrap_pyfunction!(py_hilbert_series, m)?)?;
    m.add_function(wrap_pyfunction!(py_ideal_radical, m)?)?;
    m.add_function(wrap_pyfunction!(py_primary_decomposition, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_zero_dimensional, m)?)?;
    m.add_function(wrap_pyfunction!(py_eliminate, m)?)?;

    Ok(())
}
