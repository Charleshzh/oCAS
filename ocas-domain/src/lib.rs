//! Algebraic domains and number types for oCAS.

#![warn(missing_docs)]

pub mod algebraic;
pub mod assumptions;
pub mod complex;
pub mod domain;
pub mod dual;
pub mod finite_field;
pub mod integer;
pub mod number_theory;
pub mod rational;
pub mod real_ball;

#[cfg(feature = "gmp")]
pub mod gmp_backend;

pub use algebraic::{AlgebraicElement, AlgebraicExtension, AlgebraicNumberField};
pub use assumptions::{Assumption, Assumptions, SymbolAssumptions};
pub use complex::{Complex, ComplexDomain};
pub use domain::{Domain, EuclideanDomain};
pub use finite_field::{FiniteField, FiniteFieldElement};
pub use integer::IntegerDomain;
pub use rational::RationalDomain;
pub use real_ball::{RealBall, RealBallDomain};

#[cfg(not(feature = "gmp"))]
pub use integer::Integer;
#[cfg(not(feature = "gmp"))]
pub use rational::Rational;

#[cfg(feature = "gmp")]
pub use gmp_backend::{Integer, Rational};
