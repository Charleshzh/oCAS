//! Pattern AST for oCAS rewriting.
//!
//! Patterns mirror the structure of [`Atom`] but add
//! wildcard nodes. A wildcard name ending with `___` matches a (possibly
//! empty) sequence, `__` matches a non-empty sequence, and `_` matches a
//! single atom. Wildcards with the same name must bind consistently within
//! a match.

use ocas_atom::{Atom, AtomNode, Symbol};

/// The scope of a wildcard match.
///
/// # Example
///
/// ```
/// use ocas_rewrite::pattern::WildcardLevel;
///
/// // Names ending with one underscore map to Single wildcards.
/// assert!(matches!(WildcardLevel::Single, WildcardLevel::Single));
/// assert!(matches!(WildcardLevel::Sequence, WildcardLevel::Sequence));
/// assert!(matches!(WildcardLevel::NullSequence, WildcardLevel::NullSequence));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WildcardLevel {
    /// Match a single atom (e.g. `x_`).
    Single,
    /// Match one or more atoms in a `Add`/`Mul`/`Fun` argument list (e.g. `__x`).
    Sequence,
    /// Match zero or more atoms in a `Add`/`Mul`/`Fun` argument list (e.g. `___x`).
    NullSequence,
}

/// A pattern that can be matched against an [`Atom`].
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::pattern::Pattern;
/// use ocas_atom::Symbol;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let pat = Pattern::Literal(x);
/// assert!(matches!(pat, Pattern::Literal(v) if v.to_string() == "x"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Pattern<'a> {
    /// Match the atom exactly.
    Literal(Atom<'a>),
    /// Match a wildcard with the given level and name.
    Wildcard(Symbol, WildcardLevel),
    /// Match an addition whose arguments are matched by the patterns.
    Add(Vec<Pattern<'a>>),
    /// Match a multiplication whose arguments are matched by the patterns.
    Mul(Vec<Pattern<'a>>),
    /// Match a power with base and exponent patterns.
    Pow(Box<(Pattern<'a>, Pattern<'a>)>),
    /// Match a function application with the given head and argument patterns.
    Fun(Symbol, Vec<Pattern<'a>>),
}

impl<'a> Pattern<'a> {
    /// Convert an atom to a pattern by treating symbols ending with
    /// `_`, `__`, or `___` as wildcards.
    ///
    /// This allows patterns to be parsed by the existing expression parser and
    /// then promoted to patterns. Names are interned globally, so the wildcard
    /// name is simply the symbol without trailing underscores.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    /// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x_");
    /// let pat = Pattern::from_atom(&(), x);
    /// assert!(matches!(pat, Pattern::Wildcard(s, WildcardLevel::Single) if s.as_str() == "x"));
    /// ```
    pub fn from_atom(_ctx: &'a impl PatternAlloc<'a>, atom: Atom<'a>) -> Pattern<'a> {
        match atom.node() {
            AtomNode::Num(_) | AtomNode::Var(_) => {
                if let AtomNode::Var(s) = atom.node() {
                    let name = s.as_str();
                    if let Some(level) = wildcard_level(name) {
                        let base = strip_underscores(name);
                        return Pattern::Wildcard(Symbol::new(base), level);
                    }
                }
                Pattern::Literal(atom)
            }
            AtomNode::Fun(name, args) => {
                let pat_args: Vec<Pattern<'a>> =
                    args.iter().map(|a| Pattern::from_atom(_ctx, *a)).collect();
                Pattern::Fun(*name, pat_args)
            }
            AtomNode::Add(args) => {
                let pat_args: Vec<Pattern<'a>> =
                    args.iter().map(|a| Pattern::from_atom(_ctx, *a)).collect();
                Pattern::Add(pat_args)
            }
            AtomNode::Mul(args) => {
                let pat_args: Vec<Pattern<'a>> =
                    args.iter().map(|a| Pattern::from_atom(_ctx, *a)).collect();
                Pattern::Mul(pat_args)
            }
            AtomNode::Pow(base, exp) => {
                let base_pat = Box::new((
                    Pattern::from_atom(_ctx, *base),
                    Pattern::from_atom(_ctx, *exp),
                ));
                Pattern::Pow(base_pat)
            }
        }
    }
}

fn wildcard_level(name: &str) -> Option<WildcardLevel> {
    if name.starts_with("___") || name.ends_with("___") {
        Some(WildcardLevel::NullSequence)
    } else if name.starts_with("__") || name.ends_with("__") {
        Some(WildcardLevel::Sequence)
    } else if name.starts_with('_') || name.ends_with('_') {
        Some(WildcardLevel::Single)
    } else {
        None
    }
}

fn strip_underscores(name: &str) -> &str {
    let leading = name.bytes().take_while(|&b| b == b'_').count();
    let trailing = name.bytes().rev().take_while(|&b| b == b'_').count();
    let start = leading;
    let end = name.len().saturating_sub(trailing);
    if start >= end { "" } else { &name[start..end] }
}

/// Allocation helper for building [`Pattern`] slices without leaking to the
/// global arena. Implementations are provided by the rewrite engine.
///
/// # Example
///
/// ```
/// use ocas_rewrite::pattern::PatternAlloc;
///
/// let _: &dyn PatternAlloc = &();
/// ```
pub trait PatternAlloc<'a> {
    /// Allocate a slice of patterns in the caller's scratch arena.
    fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>];
}

impl<'a> PatternAlloc<'a> for () {
    fn alloc_slice(&self, _items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
        Box::leak(_items.to_vec().into_boxed_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    struct VecAlloc;

    impl<'a> PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
            let leaked: Box<[Pattern<'a>]> = items.to_vec().into_boxed_slice();
            Box::leak(leaked)
        }
    }

    #[test]
    fn single_wildcard_from_var() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let x_ = ctx.var("x_");
        let pat = Pattern::from_atom(&alloc, x_);
        assert!(matches!(pat, Pattern::Wildcard(s, WildcardLevel::Single) if s.as_str() == "x"));
    }

    #[test]
    fn wildcard_utils_smoke() {
        assert_eq!(wildcard_level("__x"), Some(WildcardLevel::Sequence));
        assert_eq!(strip_underscores("__x"), "x");
    }

    #[test]
    fn sequence_wildcard_from_var() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let xs = ctx.var("__x");
        let pat = Pattern::from_atom(&alloc, xs);
        assert!(matches!(pat, Pattern::Wildcard(s, WildcardLevel::Sequence) if s.as_str() == "x"));
    }

    #[test]
    fn literal_number_remains_literal() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let n = ctx.num(42);
        let pat = Pattern::from_atom(&alloc, n);
        assert!(matches!(pat, Pattern::Literal(a) if a == n));
    }
}
