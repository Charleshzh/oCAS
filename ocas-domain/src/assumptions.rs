//! Assumptions and domain restrictions for symbols in oCAS.
//!
//! This module provides a lightweight predicate system for declaring
//! properties of symbolic variables (e.g. "x is real", "n is a positive
//! integer"). Solvers and simplifiers consult assumptions to choose
//! algorithms and validate solutions.

use std::fmt;
use std::ops::BitOr;

/// A single predicate that can be asserted about a symbolic variable.
///
/// Assumptions are independent — a variable can carry multiple
/// assumptions simultaneously (e.g. `Positive | Integer`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Assumption {
    /// The variable is a real number.
    Real,
    /// The variable is a complex number.
    Complex,
    /// The variable is an integer.
    Integer,
    /// The variable is a rational number.
    Rational,
    /// The variable is strictly positive (> 0).
    Positive,
    /// The variable is strictly negative (< 0).
    Negative,
    /// The variable is non-negative (≥ 0).
    NonNegative,
    /// The variable is non-positive (≤ 0).
    NonPositive,
    /// The variable is non-zero.
    NonZero,
    /// The variable is finite (not ±∞).
    Finite,
    /// The variable is even.
    Even,
    /// The variable is odd.
    Odd,
    /// The variable is prime.
    Prime,
}

impl fmt::Display for Assumption {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Assumption::Real => f.write_str("real"),
            Assumption::Complex => f.write_str("complex"),
            Assumption::Integer => f.write_str("integer"),
            Assumption::Rational => f.write_str("rational"),
            Assumption::Positive => f.write_str("positive"),
            Assumption::Negative => f.write_str("negative"),
            Assumption::NonNegative => f.write_str("non-negative"),
            Assumption::NonPositive => f.write_str("non-positive"),
            Assumption::NonZero => f.write_str("non-zero"),
            Assumption::Finite => f.write_str("finite"),
            Assumption::Even => f.write_str("even"),
            Assumption::Odd => f.write_str("odd"),
            Assumption::Prime => f.write_str("prime"),
        }
    }
}

impl Assumption {
    /// Return the assumptions that are logically implied by this one.
    ///
    /// For example, `Positive` implies `NonNegative`, `NonZero`, and `Real`.
    pub fn implied(&self) -> &'static [Assumption] {
        match self {
            Assumption::Real => &[],
            Assumption::Complex => &[Assumption::Real],
            Assumption::Integer => &[Assumption::Rational, Assumption::Real],
            Assumption::Rational => &[Assumption::Real],
            Assumption::Positive => &[
                Assumption::NonNegative,
                Assumption::NonZero,
                Assumption::Real,
            ],
            Assumption::Negative => &[
                Assumption::NonPositive,
                Assumption::NonZero,
                Assumption::Real,
            ],
            Assumption::NonNegative => &[Assumption::Real],
            Assumption::NonPositive => &[Assumption::Real],
            Assumption::NonZero => &[],
            Assumption::Finite => &[],
            Assumption::Even => &[Assumption::Integer],
            Assumption::Odd => &[Assumption::Integer],
            Assumption::Prime => &[Assumption::Integer, Assumption::Positive],
        }
    }

    /// Return the assumptions that conflict with this one.
    ///
    /// Any set containing both an assumption and one of its conflicts
    /// is inconsistent.
    pub fn conflicts(&self) -> &'static [Assumption] {
        match self {
            Assumption::Real => &[],
            Assumption::Complex => &[],
            Assumption::Integer => &[],
            Assumption::Rational => &[],
            Assumption::Positive => &[Assumption::Negative, Assumption::NonPositive],
            Assumption::Negative => &[Assumption::Positive, Assumption::NonNegative],
            Assumption::NonNegative => &[Assumption::Negative],
            Assumption::NonPositive => &[Assumption::Positive],
            Assumption::NonZero => &[],
            Assumption::Finite => &[],
            Assumption::Even => &[Assumption::Odd],
            Assumption::Odd => &[Assumption::Even],
            Assumption::Prime => &[],
        }
    }
}

/// A set of assumptions about a symbolic variable.
///
/// Internally represented as a sorted, deduplicated vector for small-N
/// efficiency. Operations close over logical implication: inserting
/// `Positive` also makes `NonNegative` and `Real` available.
///
/// # Example
///
/// ```
/// use ocas_domain::assumptions::{Assumption, Assumptions};
///
/// let mut a = Assumptions::new();
/// a.insert(Assumption::Positive);
/// a.insert(Assumption::Integer);
/// assert!(a.contains(Assumption::Real));
/// assert!(a.implies(Assumption::NonZero));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Assumptions {
    inner: Vec<Assumption>,
}

impl Assumptions {
    /// Create an empty set of assumptions.
    pub fn new() -> Self {
        Self { inner: Vec::new() }
    }

    /// Create a set containing a single assumption (and its implications).
    pub fn single(a: Assumption) -> Self {
        let mut s = Self::new();
        s.insert(a);
        s
    }

    /// Return the number of explicitly stored assumptions.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Return true if no assumptions are stored.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Check whether an assumption is entailed by this set.
    ///
    /// This checks both direct membership and logical implication.
    pub fn contains(&self, a: Assumption) -> bool {
        self.inner.contains(&a)
    }

    /// Check whether this set logically implies `other`.
    ///
    /// Returns true if every assumption in `other` is entailed by this set.
    pub fn implies(&self, other: Assumption) -> bool {
        if self.contains(other) {
            return true;
        }
        // Check if any stored assumption implies `other`.
        for &a in &self.inner {
            if a.implied().contains(&other) {
                return true;
            }
        }
        false
    }

    /// Insert an assumption and all its logical implications.
    ///
    /// Returns `false` if the insertion created a contradiction (i.e., the
    /// set is now inconsistent). The inconsistent assumptions are still
    /// stored; callers should check [`is_consistent`](Self::is_consistent).
    pub fn insert(&mut self, a: Assumption) -> bool {
        if self.inner.contains(&a) {
            return !self.conflicts_with(a);
        }
        self.inner.push(a);
        // Transitive closure: insert everything implied by `a`.
        let implied: Vec<Assumption> = a.implied().to_vec();
        for imp in implied {
            if !self.inner.contains(&imp) {
                self.inner.push(imp);
            }
        }
        self.inner.sort_unstable_by_key(|a| *a as u8);
        self.inner.dedup();
        !self.conflicts_with(a)
    }

    /// Remove an assumption from the set.
    ///
    /// Note: this does not remove assumptions that were implied by the
    /// removed one, as they may also be implied by other stored assumptions.
    pub fn remove(&mut self, a: Assumption) {
        self.inner.retain(|&x| x != a);
    }

    /// Check whether the set is consistent (no contradictory assumptions).
    pub fn is_consistent(&self) -> bool {
        for &a in &self.inner {
            if self.conflicts_with(a) {
                return false;
            }
        }
        true
    }

    /// Iterate over stored assumptions.
    pub fn iter(&self) -> impl Iterator<Item = Assumption> + '_ {
        self.inner.iter().copied()
    }

    /// Check whether `a` conflicts with any assumption already in the set.
    fn conflicts_with(&self, a: Assumption) -> bool {
        for &conflict in a.conflicts() {
            if self.inner.contains(&conflict) {
                return true;
            }
        }
        false
    }
}

impl BitOr<Assumption> for Assumption {
    type Output = Assumptions;

    fn bitor(self, rhs: Assumption) -> Assumptions {
        let mut s = Assumptions::single(self);
        s.insert(rhs);
        s
    }
}

impl BitOr<Assumption> for Assumptions {
    type Output = Assumptions;

    fn bitor(mut self, rhs: Assumption) -> Assumptions {
        self.insert(rhs);
        self
    }
}

impl BitOr<Assumptions> for Assumptions {
    type Output = Assumptions;

    fn bitor(mut self, rhs: Assumptions) -> Assumptions {
        for a in rhs.inner {
            self.insert(a);
        }
        self
    }
}

impl FromIterator<Assumption> for Assumptions {
    fn from_iter<I: IntoIterator<Item = Assumption>>(iter: I) -> Self {
        let mut s = Self::new();
        for a in iter {
            s.insert(a);
        }
        s
    }
}

impl fmt::Display for Assumptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return f.write_str("(none)");
        }
        let parts: Vec<String> = self.inner.iter().map(|a| a.to_string()).collect();
        f.write_str(&parts.join(", "))
    }
}

/// A mapping from symbol names to their domain and assumptions.
///
/// This is used by solvers and simplifiers to determine valid
/// transformations. For example, `sqrt(x^2) → x` is only valid when
/// `x` is assumed `NonNegative`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SymbolAssumptions {
    // Small-N map: linear scan over a sorted vec of (name, assumptions).
    entries: Vec<(String, Assumptions)>,
}

impl SymbolAssumptions {
    /// Create an empty symbol-assumption map.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Set the assumptions for a symbol, replacing any previous entry.
    pub fn set(&mut self, symbol: &str, assumptions: Assumptions) {
        match self
            .entries
            .binary_search_by(|(s, _)| s.as_str().cmp(symbol))
        {
            Ok(idx) => {
                self.entries[idx].1 = assumptions;
            }
            Err(idx) => {
                self.entries.insert(idx, (symbol.to_owned(), assumptions));
            }
        }
    }

    /// Get the assumptions for a symbol, if any.
    pub fn get(&self, symbol: &str) -> Option<&Assumptions> {
        match self
            .entries
            .binary_search_by(|(s, _)| s.as_str().cmp(symbol))
        {
            Ok(idx) => Some(&self.entries[idx].1),
            Err(_) => None,
        }
    }

    /// Remove assumptions for a symbol.
    pub fn remove(&mut self, symbol: &str) {
        if let Ok(idx) = self
            .entries
            .binary_search_by(|(s, _)| s.as_str().cmp(symbol))
        {
            self.entries.remove(idx);
        }
    }

    /// Check whether a symbol satisfies a specific assumption.
    pub fn check(&self, symbol: &str, assumption: Assumption) -> bool {
        self.get(symbol)
            .map(|a| a.implies(assumption))
            .unwrap_or(false)
    }

    /// Return the number of symbols with assumptions.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return true if no symbols have assumptions.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate over (symbol, assumptions) pairs.
    pub fn iter(&self) -> impl Iterator<Item = &(String, Assumptions)> {
        self.entries.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_assumptions() {
        let a = Assumptions::new();
        assert!(a.is_empty());
        assert!(a.is_consistent());
        assert!(!a.implies(Assumption::Real));
    }

    #[test]
    fn positive_implies_real_and_nonzero() {
        let a = Assumptions::single(Assumption::Positive);
        assert!(a.implies(Assumption::Real));
        assert!(a.implies(Assumption::NonNegative));
        assert!(a.implies(Assumption::NonZero));
        assert!(!a.implies(Assumption::Integer));
    }

    #[test]
    fn integer_implies_rational_and_real() {
        let a = Assumptions::single(Assumption::Integer);
        assert!(a.implies(Assumption::Rational));
        assert!(a.implies(Assumption::Real));
    }

    #[test]
    fn complex_implies_real() {
        let a = Assumptions::single(Assumption::Complex);
        assert!(a.implies(Assumption::Real));
    }

    #[test]
    fn positive_and_integer() {
        let mut a = Assumptions::new();
        a.insert(Assumption::Positive);
        a.insert(Assumption::Integer);
        assert!(a.implies(Assumption::Real));
        assert!(a.implies(Assumption::NonZero));
        assert!(a.implies(Assumption::Rational));
    }

    #[test]
    fn conflict_positive_negative() {
        let mut a = Assumptions::new();
        a.insert(Assumption::Positive);
        a.insert(Assumption::Negative);
        assert!(!a.is_consistent());
    }

    #[test]
    fn conflict_even_odd() {
        let mut a = Assumptions::new();
        a.insert(Assumption::Even);
        a.insert(Assumption::Odd);
        assert!(!a.is_consistent());
    }

    #[test]
    fn no_conflict_real_integer() {
        let mut a = Assumptions::new();
        a.insert(Assumption::Real);
        a.insert(Assumption::Integer);
        assert!(a.is_consistent());
    }

    #[test]
    fn bitor_operator() {
        let a = Assumption::Positive | Assumption::Integer;
        assert!(a.implies(Assumption::Real));
        assert!(a.implies(Assumption::Rational));
        assert!(a.implies(Assumption::NonZero));
    }

    #[test]
    fn symbol_assumptions_basics() {
        let mut sa = SymbolAssumptions::new();
        sa.set("x", Assumptions::single(Assumption::Positive));
        assert!(sa.check("x", Assumption::Real));
        assert!(sa.check("x", Assumption::NonNegative));
        assert!(!sa.check("x", Assumption::Integer));
        assert!(!sa.check("y", Assumption::Real));
    }

    #[test]
    fn symbol_assumptions_override() {
        let mut sa = SymbolAssumptions::new();
        sa.set("x", Assumptions::single(Assumption::Positive));
        sa.set("x", Assumptions::single(Assumption::Integer));
        assert!(sa.check("x", Assumption::Integer));
        assert!(!sa.check("x", Assumption::Positive));
    }

    #[test]
    fn prime_implies_integer_and_positive() {
        let a = Assumptions::single(Assumption::Prime);
        assert!(a.implies(Assumption::Integer));
        assert!(a.implies(Assumption::Positive));
        assert!(a.implies(Assumption::Real));
    }
}
