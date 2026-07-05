//! Numerical verification infrastructure for symbolic computation results.
//!
//! Cross-validates oCAS symbolic results against numeric libraries:
//! - `roots` crate for root-finding verification (feature `verify-roots`)
//! - `quadrature` crate for integration verification (feature `verify-quadrature`)

#[cfg(feature = "verify-roots")]
mod root_finding;

#[cfg(feature = "verify-quadrature")]
mod integration;
