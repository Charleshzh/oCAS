//! Calculus primitives for oCAS.
//!
//! This crate provides symbolic differentiation, integration, and series
//! expansion for [`Atom`](ocas_atom::Atom) expression trees. Results are
//! returned as new atoms and are automatically simplified using the
//! rewrite engine.
//!
//! Unresolved or partially-resolved derivatives and integrals are represented
//! as the reserved function forms `Derivative(expr, var)` and
//! `Integral(expr, var)`.

#![warn(missing_docs)]

pub mod derivative;
pub mod integral;
pub mod series;

mod pattern_alloc;
mod rules;

pub use derivative::diff;
pub use integral::integrate;
pub use series::taylor;
