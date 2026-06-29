//! Core runtime, memory management, and backend glue for oCAS.
//!
//! This crate provides the foundation used by all other oCAS crates:
//! - arena allocation for expression nodes
//! - unified error types
//! - thread pool utilities
//! - FFI glue conventions
//! - thin wrappers around numerical backends (GMP, MPFR, FLINT, etc.)

#![warn(missing_docs)]

pub mod arena;
pub mod error;
pub mod thread_pool;

#[cfg(feature = "gmp")]
pub mod gmp;
