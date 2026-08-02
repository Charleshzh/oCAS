//! Tensor symmetry specifications for canonicalisation.
//!
//! A [`TensorSpec`] declares which function heads represent tensors and how
//! their slots behave under permutation (symmetric, antisymmetric subsets,
//! cyclic).  Used by [`super::canon`] to encode tensor expressions into the
//! graph-isomorphism engine.

use std::collections::HashMap;

use crate::Symbol;

/// Slot symmetry for a tensor: a set of slot-index subsets that are
/// symmetric, antisymmetric, or form a cycle.
#[derive(Debug, Clone, Default)]
pub struct SymmetrySpec {
    /// Slot subsets whose members are interchangeable (symmetric).
    pub symmetric_subsets: Vec<Vec<usize>>,
    /// Slot subsets whose members are antisymmetric (swap flips sign).
    pub antisymmetric_subsets: Vec<Vec<usize>>,
    /// Cyclic permutation on a subset of slots.
    pub cyclic: Option<Vec<usize>>,
}

impl SymmetrySpec {
    /// No symmetry at all — every slot is independent.
    pub fn none() -> Self {
        Self {
            symmetric_subsets: Vec::new(),
            antisymmetric_subsets: Vec::new(),
            cyclic: None,
        }
    }

    /// All slots are fully symmetric.
    pub fn fully_symmetric(rank: usize) -> Self {
        Self {
            symmetric_subsets: vec![(0..rank).collect()],
            antisymmetric_subsets: Vec::new(),
            cyclic: None,
        }
    }

    /// All slots are fully antisymmetric.
    pub fn fully_antisymmetric(rank: usize) -> Self {
        Self {
            symmetric_subsets: Vec::new(),
            antisymmetric_subsets: vec![(0..rank).collect()],
            cyclic: None,
        }
    }

    /// Check whether this spec implies the given slot position should
    /// have its location "hidden" in the graph encoding (i.e. not
    /// participate in canonicalisation comparison).
    pub fn is_slot_hidden(&self, pos: usize) -> bool {
        for subset in &self.symmetric_subsets {
            if subset.contains(&pos) {
                return true;
            }
        }
        self.cyclic.as_ref().is_some_and(|c| c.contains(&pos))
    }
}

/// Complete tensor specification: which function heads are tensors and
/// what symmetries their slots have.
#[derive(Debug, Clone, Default)]
pub struct TensorRegistry {
    specs: HashMap<Symbol, SymmetrySpec>,
    /// Index group assignment: symbol → group identifier.
    /// Different groups prevent dummy index renaming across dimensions.
    index_groups: HashMap<Symbol, u64>,
}

impl TensorRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a tensor with its slot symmetry spec.
    pub fn register(&mut self, name: Symbol, spec: SymmetrySpec) {
        self.specs.insert(name, spec);
    }

    /// Set the index group for a label (e.g. "mu" → 1 for spacetime,
    /// "i" → 2 for internal).
    pub fn set_index_group(&mut self, label: Symbol, group: u64) {
        self.index_groups.insert(label, group);
    }

    /// Look up a tensor's symmetry spec.
    pub fn spec(&self, name: Symbol) -> Option<&SymmetrySpec> {
        self.specs.get(&name)
    }

    /// Look up an index label's group (0 = ungrouped/default).
    pub fn index_group(&self, label: Symbol) -> u64 {
        self.index_groups.get(&label).copied().unwrap_or(0)
    }
}
