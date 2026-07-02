//! Polynomial data structures and algorithms for oCAS.

#![warn(missing_docs)]

pub mod dense;
pub mod factor;
pub mod gcd;
pub mod groebner;
pub mod matrix;
pub mod multivariate_gcd;
pub mod roots;
pub mod sparse;

#[cfg(feature = "flint")]
pub mod flint_poly;

pub use dense::DenseUnivariatePolynomial;
pub use groebner::{GroebnerBasis, buchberger};
pub use matrix::{Matrix, MatrixError};
pub use roots::RootInterval;
pub use sparse::{
    Grevlex, Lex, MonomialOrder, SparseMultivariatePolynomial, monomial_are_coprime,
    monomial_divides, monomial_lcm,
};
