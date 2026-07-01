//! Unified error types for oCAS.

use thiserror::Error;

/// The primary error type returned by oCAS operations.
///
/// # Example
///
/// ```
/// use ocas_core::error::OcasError;
///
/// let err = OcasError::ParseError {
///     message: "unexpected token".into(),
///     span: Some((0, 3)),
/// };
/// assert_eq!(err.to_string(), "parse error at bytes 0..3: unexpected token");
/// ```
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
    use proptest::prelude::*;

    mod simple {
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
        fn display_numeric_overflow() {
            let err = OcasError::NumericOverflow;
            assert_eq!(err.to_string(), "numeric overflow");
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
    }

    mod medium {
        use super::*;

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
    }

    mod complex {
        use super::*;

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

        #[test]
        fn all_variants_round_trip_through_display() {
            let errors: Vec<OcasError> = vec![
                OcasError::ParseError {
                    message: "m".into(),
                    span: Some((1, 2)),
                },
                OcasError::ParseError {
                    message: "m".into(),
                    span: None,
                },
                OcasError::DomainError {
                    expected: "e".into(),
                    found: "f".into(),
                },
                OcasError::NumericOverflow,
                OcasError::UnsupportedOperation {
                    message: "u".into(),
                },
                OcasError::BackendError {
                    backend: "b".into(),
                    message: "m".into(),
                },
                OcasError::InvalidArgument {
                    name: "n".into(),
                    reason: "r".into(),
                },
            ];
            for err in errors {
                assert!(!err.to_string().is_empty());
                assert_eq!(err.clone(), err);
            }
        }
    }

    mod extreme {
        use super::*;

        proptest! {
            #[test]
            fn unsupported_operation_display_contains_message(message in "[a-zA-Z0-9_ ]{1,64}") {
                let err = OcasError::UnsupportedOperation { message: message.clone() };
                let text = err.to_string();
                prop_assert!(!text.is_empty());
                prop_assert!(text.contains(&message), "{text} should contain {message}");
            }

            #[test]
            fn parse_error_span_is_in_display((start, end) in (0usize..1000, 0usize..1000)) {
                let err = OcasError::ParseError {
                    message: "test".into(),
                    span: Some((start, end)),
                };
                let text = err.to_string();
                prop_assert!(text.contains("parse error"));
            }
        }
    }
}
