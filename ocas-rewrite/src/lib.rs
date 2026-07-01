//! Pattern matching and rewriting for oCAS.
//!
//! This crate builds on [`ocas_atom::Atom`] and provides:
//!
//! - A [`Pattern`] type with wildcards of three levels.
//! - A [`match_pattern`](crate::matcher::match_pattern) engine that binds
//!   wildcards to sub-expressions.
//! - A [`transform`](crate::transformer::transform) function for bottom-up
//!   tree rewriting.
//! - A [`Rule`](crate::rules::Rule) abstraction and a [`simplify`](crate::simplify::simplify)
//!   entry point.

pub mod matcher;
pub mod pattern;
pub mod rules;
pub mod simplify;
pub mod transformer;

pub use matcher::{Bindings, MatchError, match_pattern};
pub use pattern::{Pattern, WildcardLevel};
pub use rules::Rule;
pub use simplify::simplify;
pub use transformer::transform;
