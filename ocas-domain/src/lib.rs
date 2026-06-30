//! Algebraic domains and number types for oCAS.

#![warn(missing_docs)]

pub mod domain;
pub mod integer;
pub mod rational;

pub use domain::{Domain, EuclideanDomain};
