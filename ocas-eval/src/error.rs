//! Evaluation-specific error types.
//!
//! These errors are produced during numeric evaluation, compilation,
//! and JIT code generation.

use thiserror::Error;

/// Errors that can occur during expression evaluation.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum EvaluationError {
    /// A variable referenced in the expression was not provided.
    #[error("undefined variable '{name}'")]
    UndefinedVariable {
        /// Name of the undefined variable.
        name: String,
    },

    /// A type mismatch occurred (e.g. integer where float was expected).
    #[error("type mismatch: expected {expected}, found {found}")]
    TypeMismatch {
        /// Expected type description.
        expected: String,
        /// Actual type description.
        found: String,
    },

    /// Division by zero.
    #[error("division by zero")]
    DivisionByZero,

    /// A user-defined function was not found in the registry.
    #[error("function '{name}' not found")]
    FunctionNotFound {
        /// Name of the missing function.
        name: String,
    },

    /// A function was called with the wrong number of arguments.
    #[error("wrong arity for '{name}': expected {expected}, got {got}")]
    WrongArity {
        /// Name of the function.
        name: String,
        /// Expected number of arguments.
        expected: usize,
        /// Actual number of arguments.
        got: usize,
    },

    /// JIT compilation failed.
    #[error("JIT compilation error: {message}")]
    JitCompilationError {
        /// Description of the compilation failure.
        message: String,
    },

    /// The requested operation is not supported for this domain.
    #[error("unsupported operation: {message}")]
    UnsupportedOperation {
        /// Description of the unsupported operation.
        message: String,
    },
}

/// A convenient result type for evaluation operations.
pub type Result<T> = std::result::Result<T, EvaluationError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_undefined_variable() {
        let err = EvaluationError::UndefinedVariable { name: "x".into() };
        assert_eq!(err.to_string(), "undefined variable 'x'");
    }

    #[test]
    fn display_division_by_zero() {
        let err = EvaluationError::DivisionByZero;
        assert_eq!(err.to_string(), "division by zero");
    }

    #[test]
    fn display_wrong_arity() {
        let err = EvaluationError::WrongArity {
            name: "f".into(),
            expected: 2,
            got: 1,
        };
        assert_eq!(err.to_string(), "wrong arity for 'f': expected 2, got 1");
    }
}
