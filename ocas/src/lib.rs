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

// Use mimalloc as the global allocator when the feature is enabled.
// This can significantly speed up workloads with many small allocations
// (e.g. polynomial arithmetic, expression tree construction).
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Convenience prelude that brings commonly used oCAS types and functions into scope.
///
/// This is the recommended way to use oCAS in application code:
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
///
/// The prelude contains the following groups of items:
///
/// - **Expression trees**: [`Atom`], [`AtomArena`], [`AtomNode`], [`Symbol`], [`normalize`]
/// - **Calculus**: [`diff`], [`integrate`], [`taylor`], [`substitute`], [`apart`]
/// - **Solving**: [`solve_linear_rational`], [`solve_linear_integer`], [`solve_diophantine`]
/// - **Parsing**: [`parse`], [`ParseError`]
/// - **Polynomials**: [`DenseUnivariatePolynomial`], [`SparseMultivariatePolynomial`],
///   [`RationalPolynomial`], [`MonomialOrder`], [`Lex`], [`Grevlex`], `Grlex`
/// - **Domains**: [`Integer`], [`Rational`], [`RealBall`], [`Complex`], [`FiniteField`],
///   [`Domain`], [`EuclideanDomain`]
/// - **Rewriting**: [`Pattern`], [`Rule`], [`match_pattern`], [`simplify()`], [`transform`],
///   [`Bindings`], [`MatchError`], [`WildcardLevel`]
/// - **Evaluation**: [`ExpressionEvaluator`], [`FunctionMap`], [`EvaluationDomain`],
///   [`EvaluationError`], [`EvalTree`], [`Instr`], [`Instruction`], [`Slot`]
/// - **Runtime**: [`Arena`], [`OcasError`], [`Result`]
///
/// Optional backends (GMP, MPFR, FLINT, LLVM, etc.) are enabled via the
/// corresponding feature flags on this crate.
pub mod prelude {
    pub use ocas_atom::{Atom, AtomArena, AtomNode, Symbol, normalize};
    pub use ocas_calc::solve::{
        self, DiophantineSolution, SolveError, solve_diophantine, solve_linear_integer,
        solve_linear_rational, solve_polynomial_system,
    };
    pub use ocas_calc::{apart, diff, integrate, substitute, taylor};
    pub use ocas_core::arena::Arena;
    pub use ocas_core::error::{OcasError, Result};
    pub use ocas_domain::{
        AlgebraicElement, AlgebraicExtension, AlgebraicNumberField, Assumption, Assumptions,
        Complex, ComplexDomain, Domain, EuclideanDomain, FiniteField, FiniteFieldElement, Integer,
        IntegerDomain, Rational, RationalDomain, RealBall, RealBallDomain, SymbolAssumptions,
    };
    #[cfg(feature = "simd")]
    pub use ocas_eval::VectorEvaluator;
    #[cfg(feature = "jit")]
    pub use ocas_eval::jit::{JitCompiledFunction, JitEngine};
    pub use ocas_eval::{
        EvalTree, EvaluationDomain, EvaluationError, ExpressionEvaluator, FunctionMap, Instr,
        Instruction, PowfExtension, Slot,
    };
    pub use ocas_parse::{ParseError, parse};
    pub use ocas_poly::{
        DenseUnivariatePolynomial, Grevlex, Grlex, GroebnerBasis, Lex, Matrix, MatrixError,
        MonomialOrder, RationalPolynomial, RootInterval, SparseMultivariatePolynomial, buchberger,
        f4, monomial_are_coprime, monomial_divides, monomial_lcm,
    };
    pub use ocas_rewrite::{
        Bindings, MatchError, Pattern, Rule, WildcardLevel, match_pattern, simplify, transform,
    };
}

// Re-export crates for users who prefer fully qualified names. These modules are
// kept public for advanced use cases but are hidden from the generated docs; the
// stable entry points are `ocas::prelude::*` and `ocas::*` crate-root items.
#[doc(hidden)]
pub use ocas_atom;
#[doc(hidden)]
pub use ocas_calc;
#[doc(hidden)]
pub use ocas_core;
#[doc(hidden)]
pub use ocas_domain;
#[doc(hidden)]
pub use ocas_eval;
#[doc(hidden)]
pub use ocas_parse;
#[doc(hidden)]
pub use ocas_poly;
#[doc(hidden)]
pub use ocas_rewrite;

// Re-export the most common types and functions at the crate root as well.
pub use prelude::{
    Arena, Assumption, Assumptions, Atom, AtomArena, AtomNode, Bindings, Complex, ComplexDomain,
    DenseUnivariatePolynomial, DiophantineSolution, Domain, EuclideanDomain, EvalTree,
    EvaluationDomain, EvaluationError, ExpressionEvaluator, FiniteField, FiniteFieldElement,
    FunctionMap, Grevlex, GroebnerBasis, Instr, Instruction, Integer, IntegerDomain, Lex,
    MatchError, Matrix, MatrixError, MonomialOrder, OcasError, ParseError, Pattern, PowfExtension,
    Rational, RationalDomain, RationalPolynomial, RealBall, RealBallDomain, Result, RootInterval,
    Rule, Slot, SolveError, SparseMultivariatePolynomial, Symbol, SymbolAssumptions, WildcardLevel,
    apart, buchberger, diff, integrate, match_pattern, monomial_are_coprime, monomial_divides,
    monomial_lcm, normalize, parse, simplify, solve_diophantine, solve_linear_integer,
    solve_linear_rational, solve_polynomial_system, substitute, taylor, transform,
};
