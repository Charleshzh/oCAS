//! Instruction set for the stack-based evaluation VM.
//!
//! The instruction set has two representations:
//!
//! - [`Instr`]: Internal, index-based instructions used by the evaluator.
//!   All slot references are absolute stack indices.
//! - [`Instruction`]: Public, [`Slot`]-based instructions for inspection
//!   and serialization.
//!
//! # Stack layout
//!
//! ```text
//! [params (param_count)]
//! [constants (const_count)]
//! [temporaries (temp_count)]
//! [outputs (result_count)]
//! ```
//!
//! The evaluator pre-fills params from user input and constants from the
//! compiled expression, then executes instructions to compute temporaries
//! and results.

use ocas_atom::Symbol;

/// A named slot in the evaluator stack.
///
/// Used by the public [`Instruction`] type to refer to stack positions
/// semantically rather than by raw index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Slot {
    /// A parameter slot. Index is 0-based within the params area.
    Param(usize),
    /// A constant value slot. Index is 0-based within the constants area.
    Const(usize),
    /// A temporary value slot. Index is 0-based within the temporaries area.
    Temp(usize),
}

/// A public, self-documenting instruction.
///
/// Unlike [`Instr`] (which uses raw stack indices), `Instruction` uses
/// [`Slot`] to make the instruction stream readable and serializable.
#[derive(Debug, Clone, PartialEq)]
pub enum Instruction {
    /// `dst = sum of sources`
    Add(Slot, Vec<Slot>),
    /// `dst = product of sources`
    Mul(Slot, Vec<Slot>),
    /// `dst = base^exp` where exp is an integer exponent
    Pow(Slot, Slot, i64),
    /// `dst = base^exp` where exp is another slot
    Powf(Slot, Slot, Slot),
    /// `dst = builtin_function(src)`
    Fun(Slot, Symbol, Slot),
    /// `dst = external_function(srcs...)`
    ExternalFun(Slot, usize, Vec<Slot>),
    /// `dst = src` (copy)
    Assign(Slot, Slot),
}

/// An internal, index-based instruction executed by [`super::ExpressionEvaluator`].
///
/// All indices are absolute positions in the evaluator's flat stack:
/// indices `0..param_count` are parameters, `param_count..param_count+const_count`
/// are constants, and the remainder are temporaries and outputs.
#[derive(Debug, Clone)]
pub enum Instr {
    /// `stack[dst] = sum(stack[srcs[0]], stack[srcs[1]], ...)`
    Add { dst: usize, srcs: Vec<usize> },
    /// `stack[dst] = product(stack[srcs[0]], stack[srcs[1]], ...)`
    Mul { dst: usize, srcs: Vec<usize> },
    /// `stack[dst] = stack[base]^exp` where exp is an integer
    Pow { dst: usize, base: usize, exp: i64 },
    /// `stack[dst] = stack[base]^stack[exp]` (floating-point exponent)
    Powf { dst: usize, base: usize, exp: usize },
    /// `stack[dst] = builtin(stack[src])`
    BuiltinFun {
        dst: usize,
        name: Symbol,
        src: usize,
    },
    /// `stack[dst] = fns[fn_idx](&stack[srcs[0]..])`
    ExternalFun {
        dst: usize,
        fn_idx: usize,
        srcs: Vec<usize>,
    },
    /// `stack[dst] = stack[src]`
    Copy { dst: usize, src: usize },
}

impl Instr {
    /// Return the destination stack index of this instruction.
    pub fn dst(&self) -> usize {
        match self {
            Instr::Add { dst, .. }
            | Instr::Mul { dst, .. }
            | Instr::Pow { dst, .. }
            | Instr::Powf { dst, .. }
            | Instr::BuiltinFun { dst, .. }
            | Instr::ExternalFun { dst, .. }
            | Instr::Copy { dst, .. } => *dst,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_equality() {
        assert_eq!(Slot::Param(0), Slot::Param(0));
        assert_ne!(Slot::Param(0), Slot::Const(0));
        assert_eq!(Slot::Temp(3), Slot::Temp(3));
    }

    #[test]
    fn instr_dst() {
        let add = Instr::Add {
            dst: 5,
            srcs: vec![1, 2],
        };
        assert_eq!(add.dst(), 5);

        let fun = Instr::BuiltinFun {
            dst: 3,
            name: Symbol::new("sin"),
            src: 2,
        };
        assert_eq!(fun.dst(), 3);
    }
}
