//! Symbolic expression tree and rewriting for oCAS.
//!
//! This crate provides the [`Atom`] type: an arena-backed, tagged-union
//! representation of symbolic expressions. Atoms are immutable, copyable
//! references into an [`Arena`] and form the core data structure used by the
//! parser, printer, and rewrite engine.

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use ocas_core::arena::Arena;

pub mod normalize;

/// An interned symbolic name (variable, function, or constant).
///
/// Symbols are deduplicated globally and live for the remainder of the
/// process. This keeps [`Atom`] small and comparable by identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Symbol(&'static str);

impl Symbol {
    /// Create a symbol from a name, interning it globally.
    pub fn new(name: &str) -> Self {
        Self(intern(name))
    }

    /// Return the symbol's string representation.
    pub fn as_str(&self) -> &str {
        self.0
    }
}

fn intern(name: &str) -> &'static str {
    static TABLE: OnceLock<Mutex<HashMap<String, &'static str>>> = OnceLock::new();
    let table = TABLE.get_or_init(|| Mutex::new(HashMap::new()));
    let mut table = table.lock().expect("symbol interner lock poisoned");
    table.entry(name.to_owned()).or_insert_with(|| {
        let boxed = name.to_owned().into_boxed_str();
        Box::leak(boxed)
    })
}

/// A reference to an expression node allocated in an arena.
///
/// `Atom` is a small copyable handle. The actual node data lives in the
/// [`Arena`] and is freed when the arena is dropped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Atom<'a>(&'a AtomNode<'a>);

impl<'a> Atom<'a> {
    /// Access the underlying node data.
    pub fn node(&self) -> &'a AtomNode<'a> {
        self.0
    }
}

/// The concrete data stored for each expression node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AtomNode<'a> {
    /// A 64-bit signed integer literal.
    Num(i64),
    /// A named variable or constant.
    Var(Symbol),
    /// A sum of sub-expressions.
    Add(&'a [Atom<'a>]),
    /// A product of sub-expressions.
    Mul(&'a [Atom<'a>]),
    /// A power with base and exponent.
    Pow(Atom<'a>, Atom<'a>),
}

/// A context that allocates [`Atom`]s in an [`Arena`].
///
/// All construction methods are immutable from the caller's perspective;
/// mutation happens through the arena's interior mutability. Identical
/// sub-expressions are hash-consed so that structural equality implies
/// pointer equality.
pub struct AtomArena<'a> {
    arena: &'a Arena,
    cons_table: RefCell<HashMap<AtomNode<'a>, Atom<'a>>>,
}

impl<'a> AtomArena<'a> {
    /// Create an `AtomArena` backed by the given arena.
    pub fn new(arena: &'a Arena) -> Self {
        Self {
            arena,
            cons_table: RefCell::new(HashMap::new()),
        }
    }

    fn intern(&self, candidate: AtomNode<'a>) -> Atom<'a> {
        let mut table = self.cons_table.borrow_mut();
        *table
            .entry(candidate)
            .or_insert_with(|| Atom(self.arena.allocate_with(|| candidate)))
    }

    /// Create an integer literal atom.
    pub fn num(&self, value: i64) -> Atom<'a> {
        self.intern(AtomNode::Num(value))
    }

    /// Create a variable atom from a name.
    pub fn var(&self, name: &str) -> Atom<'a> {
        self.intern(AtomNode::Var(Symbol::new(name)))
    }

    /// Create an addition atom.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `args` is empty.
    pub fn add(&self, args: &[Atom<'a>]) -> Atom<'a> {
        debug_assert!(!args.is_empty(), "Add node requires at least one argument");
        let slice = self.arena.allocate_slice(args);
        self.intern(AtomNode::Add(slice))
    }

    /// Create a multiplication atom.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `args` is empty.
    pub fn mul(&self, args: &[Atom<'a>]) -> Atom<'a> {
        debug_assert!(!args.is_empty(), "Mul node requires at least one argument");
        let slice = self.arena.allocate_slice(args);
        self.intern(AtomNode::Mul(slice))
    }

    /// Create a power atom.
    pub fn pow(&self, base: Atom<'a>, exp: Atom<'a>) -> Atom<'a> {
        self.intern(AtomNode::Pow(base, exp))
    }
}

impl std::fmt::Display for Atom<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.node() {
            AtomNode::Num(n) => write!(f, "{n}"),
            AtomNode::Var(s) => write!(f, "{}", s.as_str()),
            AtomNode::Add(args) => {
                for (i, arg) in args.iter().enumerate() {
                    if i > 0 {
                        write!(f, " + ")?;
                    }
                    write_parenthesized(arg, f)?;
                }
                Ok(())
            }
            AtomNode::Mul(args) => {
                for (i, arg) in args.iter().enumerate() {
                    if i > 0 {
                        write!(f, "*")?;
                    }
                    write_parenthesized(arg, f)?;
                }
                Ok(())
            }
            AtomNode::Pow(base, exp) => {
                write_parenthesized(base, f)?;
                write!(f, "^")?;
                write_parenthesized(exp, f)
            }
        }
    }
}

fn write_parenthesized(atom: &Atom<'_>, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => write!(f, "{atom}"),
        _ => write!(f, "({atom})"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_num() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let n = ctx.num(42);
        assert_eq!(n.to_string(), "42");
        assert!(matches!(n.node(), AtomNode::Num(42)));
    }

    #[test]
    fn construct_var() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert_eq!(x.to_string(), "x");
        assert!(matches!(x.node(), AtomNode::Var(s) if s.as_str() == "x"));
    }

    #[test]
    fn construct_add() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        assert_eq!(sum.to_string(), "x + y");
    }

    #[test]
    fn construct_mul() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two = ctx.num(2);
        let prod = ctx.mul(&[x, two]);
        assert_eq!(prod.to_string(), "x*2");
    }

    #[test]
    fn construct_pow() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two = ctx.num(2);
        let pow = ctx.pow(x, two);
        assert_eq!(pow.to_string(), "x^2");
    }

    #[test]
    fn nested_expression_prints_with_parentheses() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let two = ctx.num(2);
        let squared = ctx.pow(sum, two);
        assert_eq!(squared.to_string(), "(x + y)^2");
    }

    #[test]
    fn atom_equality_uses_structure() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let a = ctx.add(&[x, y]);
        let b = ctx.add(&[x, y]);
        let c = ctx.add(&[y, x]);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn symbol_identity_is_preserved() {
        let a = Symbol::new("x");
        let b = Symbol::new("x");
        let c = Symbol::new("y");
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_eq!(a.as_str(), "x");
    }

    #[test]
    fn atom_is_copyable() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let copied = x;
        assert_eq!(x, copied);
    }

    #[test]
    fn hash_consing_reuses_identical_nodes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let a = ctx.add(&[x, y]);
        let b = ctx.add(&[x, y]);
        // Hash-consing should return the same arena pointer.
        assert!(std::ptr::eq(a.node(), b.node()));
    }

    #[test]
    fn hash_consing_distinguishes_different_nodes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let a = ctx.add(&[x, y]);
        let b = ctx.add(&[y, x]);
        assert!(!std::ptr::eq(a.node(), b.node()));
    }
}
