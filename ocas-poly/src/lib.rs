//! Polynomial data structures and algorithms for oCAS.

#![warn(missing_docs)]

pub mod dense;
pub mod factor;
pub mod gcd;
pub mod groebner;
pub mod matrix;
pub mod multivariate_gcd;
pub mod rational;
pub mod rational_reconstruction;
pub mod resultant;
pub mod roots;
pub mod sparse;

#[cfg(feature = "flint")]
pub mod flint_poly;

#[cfg(feature = "ntt")]
pub mod ntt;

#[cfg(feature = "sprs")]
pub mod sprs_backend;

pub use dense::DenseUnivariatePolynomial;
pub use groebner::{Algorithm, GroebnerBasis, buchberger, f4, f5, groebner_basis};
pub use matrix::{Matrix, MatrixError};
pub use multivariate_gcd::{bivariate_gcd, gcd_modular, lift_from_fp, reduce_mod};
pub use rational::RationalPolynomial;
pub use roots::RootInterval;
pub use sparse::{
    Grevlex, Grlex, Lex, MonomialOrder, SparseMultivariatePolynomial, monomial_are_coprime,
    monomial_divides, monomial_lcm,
};
