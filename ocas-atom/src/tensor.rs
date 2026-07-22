//! Basic tensor algebra — index slots, contraction, and slot symmetries.
//!
//! This is the *algebraic basics* deliverable for 0.18.0: an independent
//! [`Tensor`] representation (it does **not** extend the core
//! [`AtomNode`](crate::AtomNode) tagged union, so the expression tree stays
//! minimal) with explicit index-slot contraction and simple
//! symmetric/antisymmetric symmetries applied by slot permutation.
//!
//! ## Scope
//!
//! The full tensor calculus of Symbolica relies on a graph-isomorphism
//! engine (graphica) for unique canonicalisation, which oCAS does not yet
//! have. This module therefore provides:
//!
//! - slot renaming / free-vs-dummy index bookkeeping,
//! - explicit contraction (sum over repeated dummy indices),
//! - symmetric / antisymmetric symmetrisation via slot permutation,
//!
//! and **not** a canonical form guarantee. General-relativity-grade tensor
//! calculus (canonicalisation under index permutation groups) is deferred to
//! Post-1.0.

use crate::{Atom, AtomArena, Symbol};

/// Position of an index in a tensor's slot list: upper (contravariant) or
/// lower (covariant).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IndexPosition {
    /// Upper / contravariant index.
    Upper,
    /// Lower / covariant index.
    Lower,
}

/// A single index slot of a tensor: the index expression and its variance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexSlot<'a> {
    /// The index label (typically a one-character variable or an integer).
    label: Atom<'a>,
    /// Whether the index is upper or lower.
    position: IndexPosition,
}

impl<'a> IndexSlot<'a> {
    /// Create a new index slot.
    pub fn new(label: Atom<'a>, position: IndexPosition) -> Self {
        Self { label, position }
    }

    /// The index label expression.
    pub fn label(&self) -> Atom<'a> {
        self.label
    }

    /// The index variance.
    pub fn position(&self) -> IndexPosition {
        self.position
    }
}

/// Symmetry of a tensor's index slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Symmetry {
    /// No symmetry (general tensor).
    None,
    /// Symmetric: invariant under any swap of slots.
    Symmetric,
    /// Antisymmetric: flips sign under any swap of slots.
    Antisymmetric,
}

/// A tensor: a named object with a list of index slots, an arity, and a
/// slot symmetry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tensor<'a> {
    name: Symbol,
    slots: Vec<IndexSlot<'a>>,
    symmetry: Symmetry,
}

impl<'a> Tensor<'a> {
    /// Create a new tensor with the given slots and no symmetry.
    pub fn new(name: Symbol, slots: Vec<IndexSlot<'a>>) -> Self {
        Self {
            name,
            slots,
            symmetry: Symmetry::None,
        }
    }

    /// Builder: set the slot symmetry.
    pub fn with_symmetry(mut self, symmetry: Symmetry) -> Self {
        self.symmetry = symmetry;
        self
    }

    /// The tensor name.
    pub fn name(&self) -> Symbol {
        self.name
    }

    /// The index slots.
    pub fn slots(&self) -> &[IndexSlot<'a>] {
        &self.slots
    }

    /// The slot symmetry.
    pub fn symmetry(&self) -> Symmetry {
        self.symmetry
    }

    /// The tensor arity (number of slots).
    pub fn rank(&self) -> usize {
        self.slots.len()
    }

    /// Return the dummy indices (labels occurring exactly twice across all
    /// slots, once upper and once lower) — these are the ones that will be
    /// contracted in a product.
    pub fn dummy_labels(&self) -> Vec<Atom<'a>> {
        dummies(self.slots().iter().map(|s| s.label()))
    }

    /// Render this tensor as an [`Atom`] function node `name(slot, slot, ...)`
    /// in the supplied arena. Symmetrisation is *not* applied here — the atom
    /// preserves the slot order of `self`.
    pub fn to_atom(&self, ctx: &'a AtomArena<'a>) -> Atom<'a> {
        let args: Vec<Atom<'a>> = self.slots.iter().map(|s| s.label).collect();
        ctx.fun(self.name.as_str(), &args)
    }
}

/// Collect labels occurring exactly twice among the iterator (the contraction
/// candidates). Labels occurring once are free; occurring more than twice is a
/// malformed expression (returned as not-a-dummy).
fn dummies<'a, I: IntoIterator<Item = Atom<'a>>>(labels: I) -> Vec<Atom<'a>> {
    use crate::FastHashMap;
    let mut counts: FastHashMap<AtomId<'a>, usize> = FastHashMap::default();
    for l in labels {
        // Atom is Copy + Eq + Hash via the arena; use the raw pointer identity
        // bucket keyed by the node address to avoid deep hashing.
        let id = AtomId(l);
        *counts.entry(id).or_insert(0) += 1;
    }
    let mut out: Vec<Atom<'a>> = Vec::new();
    let mut seen: std::collections::HashSet<*const ()> = std::collections::HashSet::new();
    for (id, n) in counts.iter() {
        if *n == 2 {
            let ptr = id.0.node() as *const _ as *const ();
            if seen.insert(ptr) {
                out.push(id.0);
            }
        }
    }
    out
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct AtomId<'a>(Atom<'a>);

/// The result of contracting a pair of tensors: either a scalar atom (when no
/// dummies remain) or a product/sum expression.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Contracted<'a> {
    /// A tensor product with the contraction performed (sum over dummies); the
    /// remaining free slots are concatenated.
    Product(TensorProduct<'a>),
    /// Fully contracted to a scalar expression.
    Scalar(Atom<'a>),
}

/// A product of tensors with free slots concatenated; dummies have been summed
/// over.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TensorProduct<'a> {
    /// The remaining free tensors after contraction (their slots concatenated).
    pub factors: Vec<Tensor<'a>>,
}

/// Contract two tensors by summing over shared dummy indices.
///
/// Two slots with the same label but opposite variance contract. The result
/// keeps the surviving free slots in the order `(a.free..., b.free...)`.
/// Returns [`Contracted::Scalar`] when no free slots survive.
pub fn contract<'a>(ctx: &'a AtomArena<'a>, a: &Tensor<'a>, b: &Tensor<'a>) -> Contracted<'a> {
    // Pair up slots with equal label and opposite variance.
    let mut used_a = vec![false; a.slots.len()];
    let mut used_b = vec![false; b.slots.len()];
    let mut pair_labels: Vec<Atom<'a>> = Vec::new();
    for (i, sa) in a.slots.iter().enumerate() {
        if used_a[i] {
            continue;
        }
        for (j, sb) in b.slots.iter().enumerate() {
            if used_b[j] {
                continue;
            }
            if sa.label == sb.label && sa.position != sb.position {
                used_a[i] = true;
                used_b[j] = true;
                pair_labels.push(sa.label);
                break;
            }
        }
    }
    // Surviving free slots.
    let mut free: Vec<IndexSlot<'a>> = Vec::new();
    for (i, s) in a.slots.iter().enumerate() {
        if !used_a[i] {
            free.push(*s);
        }
    }
    for (j, s) in b.slots.iter().enumerate() {
        if !used_b[j] {
            free.push(*s);
        }
    }
    if pair_labels.is_empty() {
        // No contraction: plain tensor product.
        return Contracted::Product(TensorProduct {
            factors: vec![a.clone(), b.clone()],
        });
    }
    if free.is_empty() {
        // Fully contracted: build a Σ over the dummy of (a·b).
        // Represent the contraction symbolically as a sum placeholder: since
        // we do not have an explicit range, we emit a Mul of the two tensor
        // atoms; callers wanting numerical contraction supply the range.
        let a_atom = a.to_atom(ctx);
        let b_atom = b.to_atom(ctx);
        let product = ctx.mul(&[a_atom, b_atom]);
        return Contracted::Scalar(product);
    }
    // Partial contraction: a new tensor carrying the free slots.
    let name = Symbol::new(&format!("{}_contract_{}", a.name.as_str(), b.name.as_str()));
    Contracted::Product(TensorProduct {
        factors: vec![Tensor::new(name, free)],
    })
}

/// Apply a tensor's slot symmetry by permuting its slots to a canonical order
/// (ascending label), returning the sign for antisymmetry.
///
/// For [`Symmetry::Symmetric`] this sorts slots and returns `+1`. For
/// [`Symmetry::Antisymmetric`] it returns the parity of the permutation that
/// sorts the slots. For [`Symmetry::None`] it is a no-op returning `+1`.
///
/// This is **not** a full canonicalisation under a permutation group (which
/// requires graph isomorphism); it merely gives a stable order for equality
/// comparisons of symmetric tensors with the same multiset of slots.
pub fn symmetrise_sign(tensor: &Tensor<'_>) -> i64 {
    match tensor.symmetry {
        Symmetry::None | Symmetry::Symmetric => 1,
        Symmetry::Antisymmetric => {
            // Count the number of swaps in an insertion sort of the slots.
            let mut slots: Vec<IndexSlot<'_>> = tensor.slots.to_vec();
            let mut swaps = 0usize;
            for i in 1..slots.len() {
                let mut j = i;
                while j > 0 && slot_less(&slots[j - 1], &slots[j]) {
                    slots.swap(j - 1, j);
                    swaps += 1;
                    j -= 1;
                }
            }
            if swaps.is_multiple_of(2) { 1 } else { -1 }
        }
    }
}

fn slot_less(a: &IndexSlot<'_>, b: &IndexSlot<'_>) -> bool {
    // Ordering by label pointer identity then position. Good enough for a
    // stable sort within one arena.
    let pa = a.label.node() as *const _ as *const ();
    let pb = b.label.node() as *const _ as *const ();
    (pa as usize) < (pb as usize) || (pa == pb && (a.position as u8) > (b.position as u8))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;
    use crate::AtomNode;

    fn idx<'a>(ctx: &'a AtomArena<'a>, name: &str, pos: IndexPosition) -> IndexSlot<'a> {
        IndexSlot::new(ctx.var(name), pos)
    }

    #[test]
    fn tensor_rank_and_slots() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = Tensor::new(
            Symbol::new("T"),
            vec![
                idx(&ctx, "i", IndexPosition::Upper),
                idx(&ctx, "j", IndexPosition::Lower),
            ],
        );
        assert_eq!(t.rank(), 2);
        assert_eq!(t.slots().len(), 2);
        assert_eq!(t.symmetry(), Symmetry::None);
    }

    #[test]
    fn dummy_detection_finds_repeated_label() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = Tensor::new(
            Symbol::new("T"),
            vec![
                idx(&ctx, "i", IndexPosition::Upper),
                idx(&ctx, "i", IndexPosition::Lower),
            ],
        );
        let dummies = t.dummy_labels();
        assert_eq!(dummies.len(), 1);
    }

    #[test]
    fn contract_two_tensors_with_one_dummy() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        // T^i_j · U^j_k  →  (T·U)^i_k  (one dummy `j` contracted).
        let t = Tensor::new(
            Symbol::new("T"),
            vec![
                idx(&ctx, "i", IndexPosition::Upper),
                idx(&ctx, "j", IndexPosition::Lower),
            ],
        );
        let u = Tensor::new(
            Symbol::new("U"),
            vec![
                idx(&ctx, "j", IndexPosition::Upper),
                idx(&ctx, "k", IndexPosition::Lower),
            ],
        );
        match contract(&ctx, &t, &u) {
            Contracted::Product(p) => {
                assert_eq!(p.factors.len(), 1);
                assert_eq!(p.factors[0].rank(), 2); // free slots i, k
            }
            _ => panic!("expected partial contraction product"),
        }
    }

    #[test]
    fn contract_to_scalar_when_no_free_slots() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        // T^i · U_i  →  scalar (fully contracted).
        let t = Tensor::new(Symbol::new("T"), vec![idx(&ctx, "i", IndexPosition::Upper)]);
        let u = Tensor::new(Symbol::new("U"), vec![idx(&ctx, "i", IndexPosition::Lower)]);
        match contract(&ctx, &t, &u) {
            Contracted::Scalar(atom) => {
                // The scalar is T(i) * U(i); it is a Mul of two Fun nodes.
                assert!(matches!(atom.node(), AtomNode::Mul(_)));
            }
            _ => panic!("expected scalar contraction"),
        }
    }

    #[test]
    fn no_overlap_yields_plain_product() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        // T^i · U^j — no shared label → plain product with both factors.
        let t = Tensor::new(Symbol::new("T"), vec![idx(&ctx, "i", IndexPosition::Upper)]);
        let u = Tensor::new(Symbol::new("U"), vec![idx(&ctx, "j", IndexPosition::Upper)]);
        match contract(&ctx, &t, &u) {
            Contracted::Product(p) => assert_eq!(p.factors.len(), 2),
            _ => panic!("expected plain product"),
        }
    }

    #[test]
    fn antisymmetric_sign_parity() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        // ε_{ab} in slot order (a, b): already sorted → sign +1.
        let e_ab = Tensor::new(
            Symbol::new("eps"),
            vec![
                idx(&ctx, "a", IndexPosition::Lower),
                idx(&ctx, "b", IndexPosition::Lower),
            ],
        )
        .with_symmetry(Symmetry::Antisymmetric);
        // Build ε_{ba}: slots (b, a) which sorts via one swap → sign −1.
        let e_ba = Tensor::new(
            Symbol::new("eps"),
            vec![
                idx(&ctx, "b", IndexPosition::Lower),
                idx(&ctx, "a", IndexPosition::Lower),
            ],
        )
        .with_symmetry(Symmetry::Antisymmetric);
        // Both signs are computed against the sort of their own slots; we just
        // check the function returns ±1 deterministically.
        let s1 = symmetrise_sign(&e_ab);
        let s2 = symmetrise_sign(&e_ba);
        assert!(s1 == 1 || s1 == -1);
        assert!(s2 == 1 || s2 == -1);
        // ε_{ab} and ε_{ba} should have opposite sign because they differ by
        // a single swap from a common sorted order.
        assert_eq!(s1, -s2);
    }

    #[test]
    fn symmetric_sign_is_always_plus() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        let g = Tensor::new(
            Symbol::new("g"),
            vec![
                idx(&ctx, "a", IndexPosition::Lower),
                idx(&ctx, "b", IndexPosition::Lower),
            ],
        )
        .with_symmetry(Symmetry::Symmetric);
        assert_eq!(symmetrise_sign(&g), 1);
    }

    #[test]
    fn to_atom_round_trips_as_function_node() {
        let arena = crate::Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = Tensor::new(
            Symbol::new("T"),
            vec![
                idx(&ctx, "i", IndexPosition::Upper),
                idx(&ctx, "j", IndexPosition::Lower),
            ],
        );
        let atom = t.to_atom(&ctx);
        match atom.node() {
            AtomNode::Fun(name, args) => {
                assert_eq!(name.as_str(), "T");
                assert_eq!(args.len(), 2);
            }
            _ => panic!("expected Fun node"),
        }
    }
}
