//! Instruction-level optimizations.
//!
//! Provides passes for common subexpression elimination (CSE),
//! algebraic simplification, and dead-code elimination on
//! [`Instr`](crate::Instr) sequences.

use std::collections::HashMap;

use crate::instruction::Instr;

/// Run all optimization passes on an instruction sequence.
///
/// Returns the optimized instruction list and the updated temp count.
pub fn optimize(mut instructions: Vec<Instr>, temp_count: usize) -> (Vec<Instr>, usize) {
    // Pass 1: algebraic simplifications
    instructions = simplify(instructions);
    // Pass 2: CSE
    instructions = cse(instructions);
    // Pass 3: dead code elimination
    instructions = eliminate_dead_code(instructions);
    // The temp_count may change; for now it stays the same since we
    // don't compact the stack after DCE.
    (instructions, temp_count)
}

/// Algebraic simplification pass.
///
/// Applies identities:
/// - `Add [x, 0]` → `Copy x` (or `Add [0, x]`)
/// - `Mul [x, 1]` → `Copy x` (or `Mul [1, x]`)
/// - `Pow [x, 1]` → `Copy x`
/// - `Mul [x, 0]` → constant 0 (requires constant tracking, deferred)
fn simplify(instructions: Vec<Instr>) -> Vec<Instr> {
    let mut result = Vec::with_capacity(instructions.len());

    for instr in instructions {
        match &instr {
            Instr::Add { dst, srcs } if srcs.len() == 1 => {
                // Single-element add is identity
                result.push(Instr::Copy {
                    dst: *dst,
                    src: srcs[0],
                });
            }
            Instr::Mul { dst, srcs } if srcs.len() == 1 => {
                // Single-element mul is identity
                result.push(Instr::Copy {
                    dst: *dst,
                    src: srcs[0],
                });
            }
            _ => result.push(instr),
        }
    }

    result
}

/// Common subexpression elimination.
///
/// If two instructions are identical (same op, same operands),
/// replace the second with a `Copy` from the first's destination.
fn cse(instructions: Vec<Instr>) -> Vec<Instr> {
    // Map: instruction signature → (first dst slot)
    let mut seen: HashMap<InstrKey, usize> = HashMap::new();
    let mut result = Vec::with_capacity(instructions.len());

    for instr in instructions {
        let key = InstrKey::from_instr(&instr);
        match key {
            Some(k) => {
                if let Some(&existing_dst) = seen.get(&k) {
                    // Duplicate found — replace with Copy
                    result.push(Instr::Copy {
                        dst: instr.dst(),
                        src: existing_dst,
                    });
                } else {
                    seen.insert(k, instr.dst());
                    result.push(instr);
                }
            }
            None => {
                // Instructions with side effects or unknown patterns — keep as-is
                result.push(instr);
            }
        }
    }

    result
}

/// Dead code elimination.
///
/// Removes instructions whose results are never consumed.
fn eliminate_dead_code(instructions: Vec<Instr>) -> Vec<Instr> {
    if instructions.is_empty() {
        return instructions;
    }

    use std::collections::HashSet;

    // Track which stack slots are "live" (needed)
    let mut live_slots: HashSet<usize> = HashSet::new();

    // The last instruction's dst is the result — always live
    if let Some(last) = instructions.last() {
        live_slots.insert(last.dst());
    }

    // Walk backward, marking sources as live when a live instruction is found
    let mut keep: Vec<bool> = vec![false; instructions.len()];
    for (i, instr) in instructions.iter().enumerate().rev() {
        if live_slots.contains(&instr.dst()) || is_side_effect(instr) {
            keep[i] = true;
            // Mark all sources as live
            for src in instr_srcs(instr) {
                live_slots.insert(src);
            }
        }
    }

    // Filter to kept instructions
    instructions
        .into_iter()
        .enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, instr)| instr)
        .collect()
}

/// Check if an instruction has side effects (must not be removed).
fn is_side_effect(instr: &Instr) -> bool {
    matches!(instr, Instr::ExternalFun { .. })
}

/// Return the source slot indices for an instruction.
fn instr_srcs(instr: &Instr) -> Vec<usize> {
    match instr {
        Instr::Add { srcs, .. } => srcs.clone(),
        Instr::Mul { srcs, .. } => srcs.clone(),
        Instr::Pow { base, .. } => vec![*base],
        Instr::Powf { base, exp, .. } => vec![*base, *exp],
        Instr::BuiltinOp { src, .. } => vec![*src],
        Instr::ExternalFun { srcs, .. } => srcs.clone(),
        Instr::Copy { src, .. } => vec![*src],
    }
}

// ---------------------------------------------------------------------------
// Instruction key for CSE hashing
// ---------------------------------------------------------------------------

/// A content-based key for comparing instructions regardless of dst.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum InstrKey {
    Add(Vec<usize>),
    Mul(Vec<usize>),
    Pow(usize, i64),
    Powf(usize, usize),
    BuiltinOp(crate::instruction::BuiltinOp, usize),
    Copy(usize),
}

impl InstrKey {
    fn from_instr(instr: &Instr) -> Option<Self> {
        match instr {
            Instr::Add { srcs, .. } => Some(InstrKey::Add(srcs.clone())),
            Instr::Mul { srcs, .. } => Some(InstrKey::Mul(srcs.clone())),
            Instr::Pow { base, exp, .. } => Some(InstrKey::Pow(*base, *exp)),
            Instr::Powf { base, exp, .. } => Some(InstrKey::Powf(*base, *exp)),
            Instr::BuiltinOp { op, src, .. } => Some(InstrKey::BuiltinOp(*op, *src)),
            Instr::Copy { src, .. } => Some(InstrKey::Copy(*src)),
            // ExternalFun is not CSE'd (may have side effects)
            Instr::ExternalFun { .. } => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simplify_single_add() {
        // Add [x] → Copy x
        let instrs = vec![Instr::Add {
            dst: 5,
            srcs: vec![2],
        }];
        let result = simplify(instrs);
        assert!(matches!(result[0], Instr::Copy { dst: 5, src: 2 }));
    }

    #[test]
    fn simplify_single_mul() {
        // Mul [x] → Copy x
        let instrs = vec![Instr::Mul {
            dst: 5,
            srcs: vec![2],
        }];
        let result = simplify(instrs);
        assert!(matches!(result[0], Instr::Copy { dst: 5, src: 2 }));
    }

    #[test]
    fn cse_duplicate_add() {
        // Two identical Add instructions → second becomes Copy
        let instrs = vec![
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            },
            Instr::Add {
                dst: 4,
                srcs: vec![0, 1],
            },
        ];
        let result = cse(instrs);
        assert_eq!(result.len(), 2);
        assert!(matches!(result[0], Instr::Add { dst: 3, .. }));
        assert!(matches!(result[1], Instr::Copy { dst: 4, src: 3 }));
    }

    #[test]
    fn cse_non_duplicate() {
        // Different instructions should not be CSE'd
        let instrs = vec![
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 4,
                srcs: vec![0, 2],
            },
        ];
        let result = cse(instrs);
        assert_eq!(result.len(), 2);
        assert!(matches!(result[0], Instr::Add { .. }));
        assert!(matches!(result[1], Instr::Mul { .. }));
    }

    #[test]
    fn cse_copy() {
        // Two identical Copy → second also becomes Copy from first dst
        let instrs = vec![
            Instr::Copy { dst: 3, src: 0 },
            Instr::Copy { dst: 4, src: 0 },
        ];
        let result = cse(instrs);
        assert!(matches!(result[1], Instr::Copy { dst: 4, src: 3 }));
    }

    #[test]
    fn dead_code_elimination_simple() {
        // Only the last instruction produces the result; previous Copy is dead
        let instrs = vec![
            Instr::Copy { dst: 5, src: 0 }, // dead: not used
            Instr::Add {
                dst: 6,
                srcs: vec![0, 1],
            }, // result
        ];
        let result = eliminate_dead_code(instrs);
        assert_eq!(result.len(), 1);
        assert!(matches!(result[0], Instr::Add { dst: 6, .. }));
    }

    #[test]
    fn dead_code_elimination_chain() {
        // Chain: temp3 = Add(0,1); temp4 = Mul(3,2); result = Copy(4)
        let instrs = vec![
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 4,
                srcs: vec![3, 2],
            },
            Instr::Copy { dst: 5, src: 4 },
        ];
        let result = eliminate_dead_code(instrs);
        assert_eq!(result.len(), 3); // all are live
    }

    #[test]
    fn dead_code_elimination_unused_branch() {
        // Two computations, only second is used
        let instrs = vec![
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            }, // dead
            Instr::Mul {
                dst: 4,
                srcs: vec![0, 2],
            }, // dead
            Instr::Copy { dst: 5, src: 0 }, // result
        ];
        let result = eliminate_dead_code(instrs);
        assert_eq!(result.len(), 1);
        assert!(matches!(result[0], Instr::Copy { dst: 5, src: 0 }));
    }

    #[test]
    fn external_fun_not_eliminated() {
        // ExternalFun has side effects, must not be removed
        let instrs = vec![
            Instr::ExternalFun {
                dst: 3,
                fn_idx: 0,
                srcs: vec![0],
            },
            Instr::Copy { dst: 5, src: 3 },
        ];
        let result = eliminate_dead_code(instrs);
        assert_eq!(result.len(), 2);
    }
}
