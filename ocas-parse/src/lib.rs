//! Parser and printer for oCAS expressions.
//!
//! This crate provides a [`logos`]-based lexer and a small recursive-descent
//! parser that turns text into [`ocas_atom::Atom`] expression trees.

pub mod lexer;
