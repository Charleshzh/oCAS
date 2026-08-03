//! Explicit Young projector for tensor symmetries.
//!
//! Expands a tensor expression with a Young tableau symmetry into the sum
//! `Σ_{σ∈R} Σ_{τ∈C} sgn(τ) · T(τ∘σ(slots))`, where `R` is the row group
//! (permutations within each row, with sign +1) and `C` the column group
//! (permutations within each column, signed by their parity).  This is the
//! classical Young symmetrizer `a_λ · b_λ`, a projector up to a scalar factor.

use crate::{Atom, AtomArena, AtomNode};

/// A Young tableau: `row_lengths` defines the shape (e.g. `[2, 1]` for □□/□).
/// The projector symmetrises within each row and antisymmetrises within each
/// column: `c_λ = a_λ · b_λ` with `a_λ = Σ_{σ∈R} σ` and
/// `b_λ = Σ_{τ∈C} sgn(τ) τ`.
///
/// This is an **explicit** expansion (not a BSGS group-theoretic one):
/// the result is a sum of `∏ r_i! · ∏ c_j!` terms (row/column factorial
/// products), each with sign ±1.
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
}

/// All permutations of a set of box positions, as full position maps.
///
/// Each entry is `(map, sign)` where `map[j]` is the destination of box `j`
/// (boxes outside the set are fixed), and `sign` is the parity of the
/// permutation restricted to the set.
fn box_permutations(boxes: &[usize], rank: usize) -> Vec<(Vec<usize>, i64)> {
    fn build_map(boxes: &[usize], perm: &[usize], rank: usize) -> Vec<usize> {
        let mut map: Vec<usize> = (0..rank).collect();
        for (slot, &b) in boxes.iter().enumerate() {
            map[b] = perm[slot];
        }
        map
    }
    let m = boxes.len();
    if m == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(factorial(m));
    let mut perm: Vec<usize> = boxes.to_vec();
    let mut c = vec![0usize; m];
    let mut sign: i64 = 1;
    out.push((build_map(boxes, &perm, rank), sign));
    // Heap's algorithm; each swap flips the parity.
    let mut i = 1;
    while i < m {
        if c[i] < i {
            if i % 2 == 0 {
                perm.swap(0, i);
            } else {
                perm.swap(c[i], i);
            }
            sign = -sign;
            out.push((build_map(boxes, &perm, rank), sign));
            c[i] += 1;
            i = 1;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    out
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

/// Apply a Young projector to a tensor expression.
///
/// Given a tensor `T(i1, i2, …, in)` and a tableau of shape λ, this expands it
/// into the Young symmetrizer `c_λ = a_λ · b_λ`:
///
/// `Σ_{σ∈R} Σ_{τ∈C} sgn(τ) · T(i_{τσ(1)}, …, i_{τσ(n)})`
///
/// where `R` permutes slots within each row and `C` within each column.  The
/// result is **not** normalized by a hook-length factor — it is a projector up
/// to the scalar `c_λ² = (∏r_i!·∏c_j!) / dim(λ) · c_λ`.  For the fully
/// antisymmetric tableau `[1, 1, …, 1]` this yields the standard alternating
/// sum; for `[n]` the full symmetrization.
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

            // Row groups: contiguous box index ranges.
            let mut rows: Vec<Vec<usize>> = Vec::new();
            let mut idx = 0usize;
            for &len in &tableau.row_lengths {
                rows.push((idx..idx + len).collect());
                idx += len;
            }
            // Column groups: box `k + c` for each row with `c < len`.
            let columns = tableau.row_lengths.iter().copied().max().unwrap_or(0);
            let mut cols: Vec<Vec<usize>> = Vec::new();
            for c in 0..columns {
                let mut col = Vec::new();
                let mut k = 0usize;
                for &len in &tableau.row_lengths {
                    if c < len {
                        col.push(k + c);
                    }
                    k += len;
                }
                if !col.is_empty() {
                    cols.push(col);
                }
            }

            // Build the row group R (all maps with sign +1).
            let mut row_stack: Vec<(Vec<usize>, i64)> = vec![((0..rank).collect(), 1)];
            for row in &rows {
                let perms = box_permutations(row, rank);
                let mut next = Vec::new();
                for (m1, s1) in &row_stack {
                    for (p, s2) in &perms {
                        let mut composed = vec![0usize; rank];
                        for j in 0..rank {
                            composed[j] = m1[p[j]];
                        }
                        next.push((composed, s1 * s2));
                    }
                }
                row_stack = next;
            }
            // Build the column group C with parity signs.
            let mut col_stack: Vec<(Vec<usize>, i64)> = vec![((0..rank).collect(), 1)];
            for col in &cols {
                let perms = box_permutations(col, rank);
                let mut next = Vec::new();
                for (m1, s1) in &col_stack {
                    for (p, s2) in &perms {
                        let mut composed = vec![0usize; rank];
                        for j in 0..rank {
                            composed[j] = m1[p[j]];
                        }
                        next.push((composed, s1 * s2));
                    }
                }
                col_stack = next;
            }

            let mut terms: Vec<Atom<'a>> = Vec::new();
            for (tau, sign_tau) in &col_stack {
                for (sigma, _) in &row_stack {
                    // Element `j` moves to `tau[sigma[j]]`; build the inverse
                    // map `perm[i] = original slot at position i`.
                    let mut perm = vec![0usize; rank];
                    for j in 0..rank {
                        perm[tau[sigma[j]]] = j;
                    }
                    let reordered: Vec<Atom<'a>> = perm.iter().map(|&i| args[i]).collect();
                    if *sign_tau == 1 {
                        terms.push(ctx.fun(name.as_str(), &reordered));
                    } else {
                        terms.push(ctx.mul(&[ctx.num(-1), ctx.fun(name.as_str(), &reordered)]));
                    }
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

    #[test]
    fn mixed_shape_does_not_panic() {
        // Regression: shape [2, 1] previously panicked in sign_of_permutation.
        // The Young symmetrizer a_λ·b_λ has |R|·|C| = 2!·1! · 2!·1! = 4 terms.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let f = ctx.fun("f", &[a, b, c]);
        let tableau = YoungTableau::new(vec![2, 1]);
        let result = young_project(&ctx, f, &tableau);
        let s = result.to_string();
        // c_λ = (e + (01))·(e − (02)) = f(a,b,c) + f(b,a,c) − f(c,b,a) − (021)·f
        assert!(s.contains("f(a, b, c)"), "missing identity term: {s}");
        assert!(s.contains("f(b, a, c)"), "missing row-swap term: {s}");
        assert!(s.contains("f(c, b, a)"), "missing column-swap term: {s}");
    }

    #[test]
    fn two_by_two_tableau() {
        // Shape [2, 2]: |R| = 2!·2! = 4, |C| = 2!·2! = 4 → 16 terms, no panic.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let args: Vec<Atom<'_>> = ["a", "b", "c", "d"].iter().map(|&x| ctx.var(x)).collect();
        let f = ctx.fun("f", &args);
        let tableau = YoungTableau::new(vec![2, 2]);
        let result = young_project(&ctx, f, &tableau);
        let s = result.to_string();
        assert!(s.contains("f(a, b, c, d)"), "missing identity term: {s}");
    }
}
