//! oCAS — open Computer Algebra System.
//!
//! This crate is the top-level Rust API for oCAS. It re-exports the public
//! items from the lower-level crates so that most users only need one import:
//!
//! ```
//! use ocas::prelude::*;
//! ```
//!
//! Backends (GMP, MPFR, FLINT, LLVM, Python, GPL) are opt-in via feature flags
//! to keep the default build portable, including on Windows MSVC where GMP is
//! not readily available.

#![warn(missing_docs)]

/// Convenience prelude that brings commonly used oCAS types into scope.
///
/// # Examples
///
/// ```
/// use ocas::prelude::*;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// assert_eq!(x.to_string(), "x");
/// ```
pub mod prelude {
    pub use ocas_atom::{Atom, AtomArena, AtomNode, Symbol, normalize};
    pub use ocas_core::arena::Arena;
    pub use ocas_core::error::{OcasError, Result};
    pub use ocas_domain::{
        Complex, ComplexDomain, Domain, EuclideanDomain, FiniteField, FiniteFieldElement, Integer,
        IntegerDomain, Rational, RationalDomain, RealBall, RealBallDomain,
    };
    pub use ocas_parse::{ParseError, parse};
    pub use ocas_poly::{
        DenseUnivariatePolynomial, Grevlex, Lex, MonomialOrder, SparseMultivariatePolynomial,
    };
    pub use ocas_rewrite::{
        Bindings, MatchError, Pattern, Rule, WildcardLevel, match_pattern, simplify, transform,
    };
}

// Re-export crates for users who prefer fully qualified names.
pub use ocas_atom;
pub use ocas_core;
pub use ocas_domain;
pub use ocas_parse;
pub use ocas_poly;
pub use ocas_rewrite;

// Re-export the most common types at the crate root as well.
pub use prelude::{
    Arena, Atom, AtomArena, AtomNode, Bindings, Complex, ComplexDomain, DenseUnivariatePolynomial,
    Domain, EuclideanDomain, FiniteField, Integer, IntegerDomain, Lex, MatchError, MonomialOrder,
    OcasError, ParseError, Pattern, Rational, RationalDomain, RealBall, RealBallDomain, Result,
    Rule, SparseMultivariatePolynomial, Symbol, WildcardLevel, match_pattern, simplify, transform,
};
