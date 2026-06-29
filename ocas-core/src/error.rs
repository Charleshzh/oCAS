//! Unified error types for oCAS.

use std::fmt;

/// The primary error type returned by oCAS operations.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum OcasError {
    /// A parsing error with an optional source span.
    ParseError {
        /// Human-readable error message.
        message: String,
        /// Optional byte offset into the source string.
        span: Option<(usize, usize)>,
    },

    /// An operation was requested on an incompatible domain.
    DomainError {
        /// Expected domain or type.
        expected: String,
        /// Actual value or type encountered.
        found: String,
    },

    /// A numeric overflow or underflow occurred.
    NumericOverflow,

    /// The requested operation is not yet implemented or supported.
    UnsupportedOperation {
        /// Description of the unsupported operation.
        message: String,
    },

    /// A backend library returned an error.
    BackendError {
        /// Name of the backend.
        backend: String,
        /// Backend-specific error message.
        message: String,
    },

    /// An invalid argument was supplied.
    InvalidArgument {
        /// Name of the argument.
        name: String,
        /// Reason the argument is invalid.
        reason: String,
    },
}

impl fmt::Display for OcasError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OcasError::ParseError { message, span } => match span {
                Some((start, end)) => write!(
                    f,
                    "parse error at bytes {}..{}: {}",
                    start, end, message
                ),
                None => write!(f, "parse error: {}", message),
            },
            OcasError::DomainError { expected, found } => {
                write!(f, "domain error: expected {}, found {}", expected, found)
            }
            OcasError::NumericOverflow => write!(f, "numeric overflow"),
            OcasError::UnsupportedOperation { message } => {
                write!(f, "unsupported operation: {}", message)
            }
            OcasError::BackendError { backend, message } => {
                write!(f, "backend error ({}): {}", backend, message)
            }
            OcasError::InvalidArgument { name, reason } => {
                write!(f, "invalid argument '{}': {}", name, reason)
            }
        }
    }
}

impl std::error::Error for OcasError {}

/// A convenient result type for oCAS operations.
pub type Result<T> = std::result::Result<T, OcasError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_parse_error() {
        let err = OcasError::ParseError {
            message: "unexpected token".into(),
            span: Some((0, 3)),
        };
        assert_eq!(
            err.to_string(),
            "parse error at bytes 0..3: unexpected token"
        );
    }
}
