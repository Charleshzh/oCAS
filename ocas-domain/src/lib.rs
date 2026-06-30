//! Algebraic domains and number types for oCAS.

#![warn(missing_docs)]

pub mod domain;
pub mod finite_field;
pub mod integer;
pub mod rational;

pub use domain::{Domain, EuclideanDomain};
pub use finite_field::{FiniteField, FiniteFieldElement};
pub use integer::{Integer, IntegerDomain};
pub use rational::{Rational, RationalDomain};
