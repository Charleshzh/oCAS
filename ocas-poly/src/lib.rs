//! Polynomial data structures and algorithms for oCAS.

#![warn(missing_docs)]

pub mod dense;
pub mod sparse;

#[cfg(feature = "flint")]
pub mod flint_poly;

pub use dense::DenseUnivariatePolynomial;
pub use sparse::{Grevlex, Lex, MonomialOrder, SparseMultivariatePolynomial};
