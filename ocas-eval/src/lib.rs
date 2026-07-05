//! Evaluation and JIT code generation for oCAS.
//!
//! This crate provides a stack-based virtual machine for numeric evaluation
//! of symbolic expressions, an AST-to-instruction compiler, a Cranelift JIT
//! backend, and SIMD vectorized evaluation.
//!
//! # Architecture
//!
//! The evaluation pipeline:
//!
//! ```text
//! Atom (arena-backed)
//!   → EvalTree (owned intermediate)
//!   → Instr sequence (compile + optimize)
//!   → ExpressionEvaluator (stack VM)
//!   → or JitCompiledFunction (Cranelift, feature = "jit")
//!   → or VectorEvaluator (SIMD batch)
//! ```
//!
//! # Example
//!
//! ```ignore
//! use ocas_eval::{ExpressionEvaluator, EvaluationDomain};
//! use ocas_atom::AtomArena;
//! use ocas_core::arena::Arena;
//!
//! let arena = Arena::new();
//! let ctx = AtomArena::new(&arena);
//! let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);
//!
//! let eval: ExpressionEvaluator<f64> =
//!     ExpressionEvaluator::compile(&expr).unwrap();
//! let result = eval.evaluate(&[2.0]).unwrap();
//! assert_eq!(result[0], 3.0);
//! ```

pub mod domain;
pub mod error;
pub mod evaluator;
pub mod function_map;
pub mod instruction;
pub mod tree;

pub mod compile;
mod optimize;

#[cfg(feature = "jit")]
pub mod jit;

#[cfg(feature = "simd")]
pub mod simd;

#[cfg(feature = "fast-poly")]
pub mod poly_eval;

pub use compile::{compile_atom, compile_atom_with, compile_tree, compile_tree_with};
pub use domain::{EvaluationDomain, PowfExtension};
pub use error::EvaluationError;
pub use evaluator::ExpressionEvaluator;
pub use function_map::FunctionMap;
pub use instruction::{Instr, Instruction, Slot};
pub use tree::EvalTree;

#[cfg(feature = "simd")]
pub use simd::VectorEvaluator;
