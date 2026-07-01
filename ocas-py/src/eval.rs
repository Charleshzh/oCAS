//! Python evaluation module — numeric evaluation of expressions.
//!
//! Provides [`PyExpressionEvaluator`] wrapping
//! [`ocas_eval::ExpressionEvaluator<f64>`]. The evaluator compiles an
//! expression to an instruction sequence once, then evaluates it many
//! times with different parameter values.

use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_eval::{ExpressionEvaluator, compile_atom};
use ocas_parse::parse;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A compiled expression evaluator for fast numeric evaluation.
///
/// Compile once, evaluate many times:
///
/// ```python
/// from ocas import ExpressionEvaluator
///
/// evaluator = ExpressionEvaluator("x^2 + y", ["x", "y"])
/// print(evaluator.evaluate([3.0, 1.0]))  # [10.0]
/// print(evaluator.evaluate([2.0, 0.0]))  # [4.0]
/// ```
#[pyclass(name = "ExpressionEvaluator")]
pub struct PyExpressionEvaluator {
    evaluator: ExpressionEvaluator<f64>,
    param_names: Vec<String>,
}

#[pymethods]
impl PyExpressionEvaluator {
    /// Compile `input` with the given parameter names (in order).
    #[new]
    fn new(input: &str, param_names: Vec<String>) -> PyResult<Self> {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom =
            parse(&ctx, input).map_err(|e| PyValueError::new_err(format!("parse error: {e}")))?;
        let evaluator = compile_atom::<f64>(atom)
            .map_err(|e| PyValueError::new_err(format!("compile error: {e}")))?;
        Ok(PyExpressionEvaluator {
            evaluator,
            param_names,
        })
    }

    /// Evaluate the compiled expression with the given parameter values.
    ///
    /// `values` must be a list of floats with the same length as the
    /// parameter names passed to the constructor. Returns a list of result
    /// floats.
    fn evaluate(&self, values: Vec<f64>) -> PyResult<Vec<f64>> {
        if values.len() != self.param_names.len() {
            return Err(PyValueError::new_err(format!(
                "expected {} values, got {}",
                self.param_names.len(),
                values.len()
            )));
        }
        self.evaluator
            .evaluate(&values)
            .map_err(|e| PyValueError::new_err(format!("evaluation error: {e}")))
    }

    /// Number of parameters.
    #[getter]
    fn n_params(&self) -> usize {
        self.param_names.len()
    }

    fn __repr__(&self) -> String {
        format!("ExpressionEvaluator(params={:?})", self.param_names)
    }
}
