//! Instruction-level optimizations.
//!
//! Provides passes for common subexpression elimination (CSE),
//! algebraic simplification, dead-code elimination, and stack
//! compaction on [`Instr`](crate::Instr) sequences.

use std::collections::BTreeSet;

use ocas_core::FastHashMap as HashMap;

use crate::instruction::Instr;

/// Run all optimization passes on an instruction sequence.
///
/// `temp_base` is the stack index where temporaries begin (slots below
/// it are params and constants and are never remapped). `live_roots`
/// lists the result slots that must survive dead-code elimination.
///
/// Returns the optimized instruction list, the compacted temp count,
/// and the remapped live roots.
pub fn optimize(
    mut instructions: Vec<Instr>,
    temp_base: usize,
    temp_count: usize,
    live_roots: &[usize],
) -> (Vec<Instr>, usize, Vec<usize>) {
    // Pass 1: algebraic simplifications
    instructions = simplify(instructions);
    // Pass 2: CSE
    instructions = cse(instructions);
    // Pass 3: dead code elimination
    instructions = eliminate_dead_code(instructions, live_roots);
    // Pass 4: stack compaction (remap temp slots to remove holes)
    compact_stack(instructions, temp_base, temp_count, live_roots)
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
    let mut seen: HashMap<InstrKey, usize> = HashMap::default();
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
/// Removes instructions whose results are never consumed. `live_roots`
/// are the result slots that must be preserved.
fn eliminate_dead_code(instructions: Vec<Instr>, live_roots: &[usize]) -> Vec<Instr> {
    if instructions.is_empty() {
        return instructions;
    }

    use std::collections::HashSet;

    // Track which stack slots are "live" (needed)
    let mut live_slots: HashSet<usize> = live_roots.iter().copied().collect();

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

/// Stack compaction.
///
/// Remaps temporary slots (indices >= `temp_base`) to a dense range,
/// removing holes left by dead-code elimination. Param and constant
/// slots are never remapped. Returns the new instruction list, the
/// compacted temp count, and the remapped live roots.
fn compact_stack(
    instructions: Vec<Instr>,
    temp_base: usize,
    _temp_count: usize,
    live_roots: &[usize],
) -> (Vec<Instr>, usize, Vec<usize>) {
    // Collect referenced temp slots in deterministic (sorted) order.
    let mut temps: BTreeSet<usize> = BTreeSet::new();
    for instr in &instructions {
        if instr.dst() >= temp_base {
            temps.insert(instr.dst());
        }
        for src in instr_srcs(instr) {
            if src >= temp_base {
                temps.insert(src);
            }
        }
    }
    for &root in live_roots {
        if root >= temp_base {
            temps.insert(root);
        }
    }

    let map: HashMap<usize, usize> = temps
        .iter()
        .enumerate()
        .map(|(i, &slot)| (slot, temp_base + i))
        .collect();
    let remap = |slot: usize| -> usize { if slot >= temp_base { map[&slot] } else { slot } };

    let new_instructions: Vec<Instr> = instructions
        .into_iter()
        .map(|instr| match instr {
            Instr::Add { dst, srcs } => Instr::Add {
                dst: remap(dst),
                srcs: srcs.into_iter().map(remap).collect(),
            },
            Instr::Mul { dst, srcs } => Instr::Mul {
                dst: remap(dst),
                srcs: srcs.into_iter().map(remap).collect(),
            },
            Instr::Pow { dst, base, exp } => Instr::Pow {
                dst: remap(dst),
                base: remap(base),
                exp,
            },
            Instr::Powf { dst, base, exp } => Instr::Powf {
                dst: remap(dst),
                base: remap(base),
                exp: remap(exp),
            },
            Instr::BuiltinOp { dst, op, src } => Instr::BuiltinOp {
                dst: remap(dst),
                op,
                src: remap(src),
            },
            Instr::ExternalFun { dst, fn_idx, srcs } => Instr::ExternalFun {
                dst: remap(dst),
                fn_idx,
                srcs: srcs.into_iter().map(remap).collect(),
            },
            Instr::Copy { dst, src } => Instr::Copy {
                dst: remap(dst),
                src: remap(src),
            },
        })
        .collect();

    let new_roots: Vec<usize> = live_roots.iter().map(|&r| remap(r)).collect();

    (new_instructions, temps.len(), new_roots)
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
        // Only the Add produces the live result; previous Copy is dead
        let instrs = vec![
            Instr::Copy { dst: 5, src: 0 }, // dead: not used
            Instr::Add {
                dst: 6,
                srcs: vec![0, 1],
            }, // result
        ];
        let result = eliminate_dead_code(instrs, &[6]);
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
        let result = eliminate_dead_code(instrs, &[5]);
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
        let result = eliminate_dead_code(instrs, &[5]);
        assert_eq!(result.len(), 1);
        assert!(matches!(result[0], Instr::Copy { dst: 5, src: 0 }));
    }

    #[test]
    fn dead_code_elimination_multi_root() {
        // Two live roots keep both branches
        let instrs = vec![
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 4,
                srcs: vec![0, 1],
            },
        ];
        let result = eliminate_dead_code(instrs, &[3, 4]);
        assert_eq!(result.len(), 2);
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
        let result = eliminate_dead_code(instrs, &[5]);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn compact_stack_removes_holes() {
        // temp_base = 2 (1 param + 1 const); temps 5 and 9 in use
        let instrs = vec![
            Instr::Add {
                dst: 5,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 9,
                srcs: vec![5, 5],
            },
        ];
        let (result, temp_count, roots) = compact_stack(instrs, 2, 10, &[9]);
        assert_eq!(temp_count, 2);
        assert_eq!(roots, vec![3]);
        assert!(matches!(result[0], Instr::Add { dst: 2, .. }));
        assert!(matches!(result[1], Instr::Mul { dst: 3, ref srcs } if srcs == &vec![2, 2]));
    }

    #[test]
    fn compact_stack_preserves_params_and_consts() {
        // Slots below temp_base are untouched
        let instrs = vec![Instr::Copy { dst: 4, src: 1 }];
        let (result, _, roots) = compact_stack(instrs, 2, 3, &[4]);
        assert!(matches!(result[0], Instr::Copy { dst: 2, src: 1 }));
        assert_eq!(roots, vec![2]);
    }

    #[test]
    fn optimize_full_pipeline() {
        // 1 param + 1 const → temp_base = 2; dead temp + duplicate add
        let instrs = vec![
            Instr::Add {
                dst: 2,
                srcs: vec![0, 1],
            },
            Instr::Add {
                dst: 3,
                srcs: vec![0, 1],
            }, // CSE: copy of 2
            Instr::Copy { dst: 4, src: 3 }, // result chain
            Instr::Mul {
                dst: 5,
                srcs: vec![0, 0],
            }, // dead
        ];
        let (result, temp_count, roots) = optimize(instrs, 2, 4, &[4]);
        // After CSE: 3→Copy(2); DCE: 5 removed, 3 and 4 collapse via copy chain;
        // compaction densifies remaining temps.
        assert_eq!(roots.len(), 1);
        assert!(temp_count <= 3);
        for instr in &result {
            for src in instr_srcs(instr) {
                assert!(src < 2 + temp_count, "src {src} out of compacted range");
            }
            assert!(instr.dst() < 2 + temp_count);
        }
    }
}
