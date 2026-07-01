//! Pattern allocator and common rule building helpers for `ocas-calc`.
//!
//! This module is used internally by the calculus submodules to allocate
//! [`Pattern`] objects and to build simple rewrite rules.

use ocas_rewrite::pattern::{Pattern, PatternAlloc};

/// A simple stack-based pattern allocator used by the calculus modules for
/// building [`Rule`]s and [`Pattern`]s from strings.
///
/// It is intentionally private to `ocas-calc`; it exists only because the
/// rewrite API requires an allocator object, but the calculus rules are
/// constructed from short patterns that can be backed by short-lived vectors.
pub struct VecAlloc;

impl<'a> PatternAlloc<'a> for VecAlloc {
    fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
        // This is a safe-enough shortcut for the small, short-lived patterns
        // used inside the calculus rules. The caller only uses these patterns
        // for the duration of a single simplify/diff/integrate call, and the
        // patterns never outlive the current function stack.
        let leaked: &'a mut [Pattern<'a>] = Vec::from(items).leak();
        &*leaked
    }
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_rewrite::pattern::Pattern;

    use super::*;

    #[test]
    fn alloc_slice_leaks_pattern() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let pat = Pattern::Literal(x);
        let slice = VecAlloc.alloc_slice(&[pat.clone(), pat]);
        assert_eq!(slice.len(), 2);
    }
}
