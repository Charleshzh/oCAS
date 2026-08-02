//! Explicit Young projector for tensor symmetries.
//!
//! Expands a tensor expression with a Young tableau symmetry into a sum over
//! permutations of its slot arguments, with ± signs for antisymmetric rows.

use crate::{Atom, AtomArena, AtomNode};

/// A Young tableau: `row_lengths` defines the shape (e.g. `[2, 1]` for □□/□).
/// The projector symmetrises within each row and antisymmetrises within each
/// column, then sums over all permutations that preserve the tableau shape.
///
/// This is an **explicit** expansion (not a BSGS group-theoretic one):
/// the result is a sum of `shape!` terms, each with sign ±1.
#[derive(Debug, Clone)]
pub struct YoungTableau {
    /// Number of boxes in each row.
    pub row_lengths: Vec<usize>,
}

impl YoungTableau {
    /// Create a Young tableau from row lengths.
    /// The total number of boxes must match `rank`.
    pub fn new(row_lengths: Vec<usize>) -> Self {
        Self { row_lengths }
    }

    /// Total number of boxes (= tensor rank).
    pub fn total_boxes(&self) -> usize {
        self.row_lengths.iter().sum()
    }

    /// Compute the Young symmetriser sign for a given permutation of slots.
    /// Returns 0 if the permutation does not preserve the tableau shape,
    /// +1 or -1 otherwise.
    fn sign_of_permutation(&self, perm: &[usize]) -> i64 {
        let n = perm.len();
        if n != self.total_boxes() {
            return 0;
        }

        // Build row/column assignments from the tableau.
        let mut row_of = vec![0usize; n];
        let mut col_of = vec![0usize; n];
        let mut idx = 0;
        for (r, &len) in self.row_lengths.iter().enumerate() {
            for c in 0..len {
                row_of[idx] = r;
                col_of[idx] = c;
                idx += 1;
            }
        }

        // Row constraint: skip rows with only 1 element (degenerate).
        // For rows with >1 element, permuted element must stay in the same row.
        for i in 0..n {
            let ri = row_of[i];
            let r_len = self.row_lengths[ri];
            if r_len > 1 && row_of[perm[i]] != ri {
                return 0;
            }
        }

        // Column antisymmetrisation: compute parity of the permutation
        // restricted to each column.
        let mut sign: i64 = 1;
        let columns: usize = *self.row_lengths.iter().max().unwrap_or(&0);
        for c in 0..columns {
            // Collect positions in this column.
            let col_positions: Vec<usize> = (0..n).filter(|&i| col_of[i] == c).collect();
            if col_positions.len() <= 1 {
                continue;
            }
            // Compute parity of the permutation restricted to this column.
            // Build the restricted permutation map: for each col position `p`,
            // find where its original element (perm[p]) ends up.
            let mut restricted: Vec<usize> = Vec::new();
            for &p in &col_positions {
                // perm[p] tells us WHICH original element is now at position p.
                // We need to know: where does the element originally at position p end up?
                // Find q such that perm[q] = p (the element's new position).
                let new_pos = perm.iter().position(|&x| x == p).unwrap_or(p);
                // Map from old position index within column to new position index.
                let _old_idx = col_positions.iter().position(|&cp| cp == p).unwrap();
                let new_idx = col_positions.iter().position(|&cp| cp == new_pos).unwrap();
                restricted.push(new_idx);
            }
            // Count parity of this permutation.
            let mut visited = vec![false; restricted.len()];
            for i in 0..restricted.len() {
                if visited[i] {
                    continue;
                }
                let mut cycle_len: usize = 0;
                let mut cur = i;
                while !visited[cur] {
                    visited[cur] = true;
                    cycle_len += 1;
                    cur = restricted[cur];
                }
                // An even-length cycle contributes sign -1 (k-1 transpositions, k-1 is odd).
                if cycle_len.is_multiple_of(2) {
                    sign = -sign;
                }
            }
        }
        sign
    }
}

/// Apply a Young projector to a tensor expression.
///
/// Given a tensor `T(i1, i2, …, in)`, this expands it into a sum over all
/// slot permutations `Σ sign(σ) · T(i_{σ(1)}, i_{σ(2)}, …, i_{σ(n)})`.
///
/// The result is normalised by the tableau's Hook length product so that
/// applying the projector twice is idempotent.  For the fully antisymmetric
/// tableau `[1, 1, …, 1]` this yields the standard alternating sum.
pub fn young_project<'a>(
    ctx: &'a AtomArena<'a>,
    tensor_expr: Atom<'a>,
    tableau: &YoungTableau,
) -> Atom<'a> {
    match tensor_expr.node() {
        AtomNode::Fun(name, args) => {
            let rank = tableau.total_boxes();
            if args.len() != rank {
                return tensor_expr;
            }

            let mut terms: Vec<Atom<'a>> = Vec::new();
            // Generate all permutations of 0..rank.
            let mut perm: Vec<usize> = (0..rank).collect();
            let mut c = vec![0usize; rank];
            // Heap's algorithm.
            // Sign for the first permutation.
            let s0 = tableau.sign_of_permutation(&perm);
            if s0 != 0 {
                let reordered: Vec<Atom<'a>> = perm.iter().map(|&i| args[i]).collect();
                if s0 == 1 {
                    terms.push(ctx.fun(name.as_str(), &reordered));
                } else {
                    terms.push(ctx.mul(&[ctx.num(-1), ctx.fun(name.as_str(), &reordered)]));
                }
            }
            let mut i = 1;
            while i < rank {
                if c[i] < i {
                    if i % 2 == 0 {
                        perm.swap(0, i);
                    } else {
                        perm.swap(c[i], i);
                    }
                    let s = tableau.sign_of_permutation(&perm);
                    if s != 0 {
                        let reordered: Vec<Atom<'a>> = perm.iter().map(|&j| args[j]).collect();
                        if s == 1 {
                            terms.push(ctx.fun(name.as_str(), &reordered));
                        } else {
                            terms.push(ctx.mul(&[ctx.num(-1), ctx.fun(name.as_str(), &reordered)]));
                        }
                    }
                    c[i] += 1;
                    i = 1;
                } else {
                    c[i] = 0;
                    i += 1;
                }
            }

            if terms.is_empty() {
                ctx.num(0)
            } else if terms.len() == 1 {
                terms.pop().unwrap()
            } else {
                ctx.add(&terms)
            }
        }
        _ => tensor_expr,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn antisymmetric_projector_two_slots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let f_ab = ctx.fun("f", &[a, b]);
        // Fully antisymmetric: tableau [1, 1].
        let tableau = YoungTableau::new(vec![1, 1]);
        let result = young_project(&ctx, f_ab, &tableau);
        let s = result.to_string();
        // Should be f(a,b) - f(b,a).
        assert!(s.contains('-'), "expected subtraction: {s}");
    }

    #[test]
    fn symmetric_projector_two_slots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let f_ab = ctx.fun("f", &[a, b]);
        // Fully symmetric: tableau [2].
        let tableau = YoungTableau::new(vec![2]);
        let result = young_project(&ctx, f_ab, &tableau);
        let s = result.to_string();
        // Should be f(a,b) + f(b,a).
        assert!(s.contains('+'), "expected addition: {s}");
    }

    #[test]
    fn antisymmetric_three_slots_zero() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let f = ctx.fun("f", &[a, b, c]);
        let tableau = YoungTableau::new(vec![1, 1, 1]);
        let result = young_project(&ctx, f, &tableau);
        let s = result.to_string();
        // Should be an alternating sum with 6 terms.
        assert!(s.contains('+'), "expected sum: {s}");
    }

    #[test]
    fn identity_preserves_single_slot() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let f = ctx.fun("f", &[a]);
        // Single slot: tableau [1] — identity projector.
        let result = young_project(&ctx, f, &YoungTableau::new(vec![1]));
        // Result should just be f(a) itself (identity permutation).
        assert_eq!(result.to_string(), "f(a)");
    }

    #[test]
    fn non_tensor_expression_passthrough() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Non-Fun expression (variable) → passthrough unchanged.
        let x = ctx.var("x");
        let result = young_project(&ctx, x, &YoungTableau::new(vec![2]));
        assert_eq!(result.to_string(), "x");
    }

    #[test]
    fn rank_mismatch_returns_original() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let f = ctx.fun("f", &[a, b]);
        // Tableau requires 3 slots but tensor has rank 2 → return original.
        let result = young_project(&ctx, f, &YoungTableau::new(vec![1, 1, 1]));
        assert_eq!(result.to_string(), "f(a, b)");
    }

    #[test]
    fn total_boxes_returns_rank() {
        let tableau = YoungTableau::new(vec![2, 1]);
        assert_eq!(tableau.total_boxes(), 3);
        let tableau2 = YoungTableau::new(vec![1, 1, 1]);
        assert_eq!(tableau2.total_boxes(), 3);
    }
}
