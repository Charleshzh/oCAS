//! Unified error types for oCAS.

use thiserror::Error;

/// The primary error type returned by oCAS operations.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum OcasError {
    /// A parsing error with an optional source span.
    #[error("parse error{}: {message}", match span {
        Some((start, end)) => format!(" at bytes {start}..{end}"),
        None => String::new(),
    })]
    ParseError {
        /// Human-readable error message.
        message: String,
        /// Optional byte offset into the source string.
        span: Option<(usize, usize)>,
    },

    /// An operation was requested on an incompatible domain.
    #[error("domain error: expected {expected}, found {found}")]
    DomainError {
        /// Expected domain or type.
        expected: String,
        /// Actual value or type encountered.
        found: String,
    },

    /// A numeric overflow or underflow occurred.
    #[error("numeric overflow")]
    NumericOverflow,

    /// The requested operation is not yet implemented or supported.
    #[error("unsupported operation: {message}")]
    UnsupportedOperation {
        /// Description of the unsupported operation.
        message: String,
    },

    /// A backend library returned an error.
    #[error("backend error ({backend}): {message}")]
    BackendError {
        /// Name of the backend.
        backend: String,
        /// Backend-specific error message.
        message: String,
    },

    /// An invalid argument was supplied.
    #[error("invalid argument '{name}': {reason}")]
    InvalidArgument {
        /// Name of the argument.
        name: String,
        /// Reason the argument is invalid.
        reason: String,
    },
}

/// A convenient result type for oCAS operations.
pub type Result<T> = std::result::Result<T, OcasError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_parse_error_with_span() {
        let err = OcasError::ParseError {
            message: "unexpected token".into(),
            span: Some((0, 3)),
        };
        assert_eq!(
            err.to_string(),
            "parse error at bytes 0..3: unexpected token"
        );
    }

    #[test]
    fn display_parse_error_without_span() {
        let err = OcasError::ParseError {
            message: "unexpected end of input".into(),
            span: None,
        };
        assert_eq!(err.to_string(), "parse error: unexpected end of input");
    }

    #[test]
    fn display_domain_error() {
        let err = OcasError::DomainError {
            expected: "integer".into(),
            found: "rational".into(),
        };
        assert_eq!(
            err.to_string(),
            "domain error: expected integer, found rational"
        );
    }

    #[test]
    fn display_numeric_overflow() {
        let err = OcasError::NumericOverflow;
        assert_eq!(err.to_string(), "numeric overflow");
    }

    #[test]
    fn display_unsupported_operation() {
        let err = OcasError::UnsupportedOperation {
            message: "symbolic integration not implemented".into(),
        };
        assert_eq!(
            err.to_string(),
            "unsupported operation: symbolic integration not implemented"
        );
    }

    #[test]
    fn display_backend_error() {
        let err = OcasError::BackendError {
            backend: "gmp".into(),
            message: "division by zero".into(),
        };
        assert_eq!(err.to_string(), "backend error (gmp): division by zero");
    }

    #[test]
    fn display_invalid_argument() {
        let err = OcasError::InvalidArgument {
            name: "threads".into(),
            reason: "must be greater than zero".into(),
        };
        assert_eq!(
            err.to_string(),
            "invalid argument 'threads': must be greater than zero"
        );
    }

    #[test]
    fn error_implements_std_error() {
        let err = OcasError::NumericOverflow;
        let dyn_err: &dyn std::error::Error = &err;
        assert_eq!(dyn_err.to_string(), "numeric overflow");
    }

    #[test]
    fn error_clone_and_equality() {
        let err = OcasError::ParseError {
            message: "unexpected token".into(),
            span: Some((0, 3)),
        };
        let cloned = err.clone();
        assert_eq!(err, cloned);
    }

    #[test]
    fn result_type_alias_compiles() {
        fn returns_result() -> Result<i32> {
            Ok(42)
        }
        assert_eq!(returns_result().unwrap(), 42);
    }
}
