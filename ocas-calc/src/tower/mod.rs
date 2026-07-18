//! Differential field towers for Risch integration.
//!
//! A tower represents the differential field `ℚ(x, t₁, …, tₙ)` in which
//! the integrand lives. Elements are stored as multivariate rational
//! functions over `ℚ` whose variables are the tower generators (see
//! [`convert`]).

pub(crate) mod build;
pub mod convert;
pub(crate) mod elem;
