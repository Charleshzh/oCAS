//! Symbolic expression tree and rewriting for oCAS.
//!
//! This crate provides the [`Atom`] type: an arena-backed, tagged-union
//! representation of symbolic expressions. Atoms are immutable, copyable
//! references into an [`Arena`] and form the core data structure used by the
//! parser, printer, and rewrite engine.

use std::cell::RefCell;
use std::sync::{Mutex, OnceLock};

use ocas_core::FastHashMap;
use ocas_core::arena::Arena;

pub mod normalize;
pub mod tensor;
pub mod walk;
pub mod workspace;

/// An interned symbolic name (variable, function, or constant).
///
/// Symbols are deduplicated globally and live for the remainder of the
/// process. This keeps [`Atom`] small and comparable by identity.
///
/// # Example
///
/// ```
/// use ocas_atom::Symbol;
///
/// let x = Symbol::new("x");
/// let also_x = Symbol::new("x");
/// assert_eq!(x, also_x);
/// assert_eq!(x.as_str(), "x");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Symbol(&'static str);

impl Symbol {
    /// Create a symbol from a name, interning it globally.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::Symbol;
    ///
    /// let sym = Symbol::new("my_var");
    /// assert_eq!(sym.as_str(), "my_var");
    /// ```
    pub fn new(name: &str) -> Self {
        Self(intern(name))
    }

    /// Return the symbol's string representation.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::Symbol;
    ///
    /// let sym = Symbol::new("y");
    /// assert_eq!(sym.as_str(), "y");
    /// ```
    pub fn as_str(&self) -> &str {
        self.0
    }
}

fn intern(name: &str) -> &'static str {
    static TABLE: OnceLock<Mutex<FastHashMap<String, &'static str>>> = OnceLock::new();
    let table = TABLE.get_or_init(|| Mutex::new(FastHashMap::default()));
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
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let two = ctx.num(2);
/// let expr = ctx.pow(x, two);
/// assert_eq!(expr.to_string(), "x^2");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Atom<'a>(&'a AtomNode<'a>);

impl<'a> Atom<'a> {
    /// Access the underlying node data.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::{AtomArena, AtomNode};
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// assert!(matches!(x.node(), AtomNode::Var(_)));
    /// ```
    pub fn node(&self) -> &'a AtomNode<'a> {
        self.0
    }

    /// Returns the direct children of this atom, in left-to-right order.
    ///
    /// `Num`, `Var`, and `Pow` report no children through this API; use
    /// [`Self::binary_children`] for the two operands of `Pow`.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let y = ctx.var("y");
    /// let sum = ctx.add(&[x, y, ctx.num(1)]);
    /// assert_eq!(sum.children().len(), 3);
    /// ```
    pub fn children(&self) -> &'a [Atom<'a>] {
        match self.node() {
            AtomNode::Num(_) | AtomNode::Var(_) => &[],
            AtomNode::Fun(_, args) | AtomNode::Add(args) | AtomNode::Mul(args) => args,
            AtomNode::Pow(base, exp) => {
                // This function cannot return a dynamically-allocated slice,
                // so callers that need the two-element slice should use
                // [`Self::binary_children`]. For now, `Pow` reports no children
                // through this API to keep the return type a plain slice.
                let _ = (base, exp);
                &[]
            }
        }
    }

    /// If this atom is a binary operator, return its two operands.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let y = ctx.var("y");
    /// let power = ctx.pow(x, y);
    /// let (base, exp) = power.binary_children().unwrap();
    /// assert_eq!(base.to_string(), "x");
    /// assert_eq!(exp.to_string(), "y");
    /// ```
    pub fn binary_children(&self) -> Option<(Atom<'a>, Atom<'a>)> {
        match self.node() {
            AtomNode::Pow(base, exp) => Some((*base, *exp)),
            _ => None,
        }
    }
}

/// The concrete data stored for each expression node.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, AtomNode};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// match x.node() {
///     AtomNode::Var(s) => assert_eq!(s.as_str(), "x"),
///     _ => panic!("expected variable"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AtomNode<'a> {
    /// A 64-bit signed integer literal.
    Num(i64),
    /// A named variable or constant.
    Var(Symbol),
    /// A named function applied to a list of arguments.
    Fun(Symbol, &'a [Atom<'a>]),
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
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let sum = ctx.add(&[x, y]);
/// assert_eq!(sum.to_string(), "x + y");
/// ```
pub struct AtomArena<'a> {
    arena: &'a Arena,
    cons_table: RefCell<FastHashMap<AtomNode<'a>, Atom<'a>>>,
}

impl<'a> AtomArena<'a> {
    /// Create an `AtomArena` backed by the given arena.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let n = ctx.num(42);
    /// assert_eq!(n.to_string(), "42");
    /// ```
    pub fn new(arena: &'a Arena) -> Self {
        Self {
            arena,
            cons_table: RefCell::new(FastHashMap::default()),
        }
    }

    fn intern(&self, candidate: AtomNode<'a>) -> Atom<'a> {
        let mut table = self.cons_table.borrow_mut();
        *table
            .entry(candidate)
            .or_insert_with(|| Atom(self.arena.allocate_with(|| candidate)))
    }

    /// Create an integer literal atom.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let n = ctx.num(7);
    /// assert_eq!(n.to_string(), "7");
    /// ```
    pub fn num(&self, value: i64) -> Atom<'a> {
        self.intern(AtomNode::Num(value))
    }

    /// Create a variable atom from a name.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// assert_eq!(x.to_string(), "x");
    /// ```
    pub fn var(&self, name: &str) -> Atom<'a> {
        self.intern(AtomNode::Var(Symbol::new(name)))
    }

    /// Create a function application atom.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `args` is empty.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let f = ctx.fun("sin", &[x]);
    /// assert_eq!(f.to_string(), "sin(x)");
    /// ```
    pub fn fun(&self, name: &str, args: &[Atom<'a>]) -> Atom<'a> {
        debug_assert!(!args.is_empty(), "Fun node requires at least one argument");
        let slice = self.arena.allocate_slice(args);
        self.intern(AtomNode::Fun(Symbol::new(name), slice))
    }

    /// Create an addition atom.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if `args` is empty.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let y = ctx.var("y");
    /// let sum = ctx.add(&[x, y]);
    /// assert_eq!(sum.to_string(), "x + y");
    /// ```
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
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let y = ctx.var("y");
    /// let product = ctx.mul(&[x, y]);
    /// assert_eq!(product.to_string(), "x*y");
    /// ```
    pub fn mul(&self, args: &[Atom<'a>]) -> Atom<'a> {
        debug_assert!(!args.is_empty(), "Mul node requires at least one argument");
        let slice = self.arena.allocate_slice(args);
        self.intern(AtomNode::Mul(slice))
    }

    /// Create a power atom.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let p = ctx.pow(x, ctx.num(3));
    /// assert_eq!(p.to_string(), "x^3");
    /// ```
    pub fn pow(&self, base: Atom<'a>, exp: Atom<'a>) -> Atom<'a> {
        self.intern(AtomNode::Pow(base, exp))
    }
}

impl std::fmt::Display for Atom<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.node() {
            AtomNode::Num(n) => write!(f, "{n}"),
            AtomNode::Var(s) => write!(f, "{}", s.as_str()),
            AtomNode::Fun(name, args) => {
                write!(f, "{}(", name.as_str())?;
                for (i, arg) in args.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{arg}")?;
                }
                write!(f, ")")
            }
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
    fn construct_fun() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sin = ctx.fun("sin", &[x]);
        assert_eq!(sin.to_string(), "sin(x)");
    }

    #[test]
    fn fun_with_multiple_args() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let f = ctx.fun("f", &[x, y]);
        assert_eq!(f.to_string(), "f(x, y)");
    }

    #[test]
    fn children_returns_direct_subexpressions() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        assert_eq!(sum.children(), &[x, y]);
        assert_eq!(x.children(), &[]);
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

#[cfg(test)]
mod proptests {
    use super::*;
    use ocas_core::arena::Arena;
    use proptest::prelude::*;

    /// Owned expression tree used for property-test generation.
    #[derive(Debug, Clone)]
    enum PropExpr {
        Num(i64),
        Var(&'static str),
        Fun(&'static str, Vec<PropExpr>),
        Add(Vec<PropExpr>),
        Mul(Vec<PropExpr>),
        Pow(Box<PropExpr>, Box<PropExpr>),
    }

    fn build_atom<'a>(ctx: &AtomArena<'a>, expr: &PropExpr) -> Atom<'a> {
        match expr {
            PropExpr::Num(n) => ctx.num(*n),
            PropExpr::Var(name) => ctx.var(name),
            PropExpr::Fun(name, args) => {
                let atoms: Vec<Atom<'a>> = args.iter().map(|a| build_atom(ctx, a)).collect();
                ctx.fun(name, &atoms)
            }
            PropExpr::Add(args) => {
                let atoms: Vec<Atom<'a>> = args.iter().map(|a| build_atom(ctx, a)).collect();
                ctx.add(&atoms)
            }
            PropExpr::Mul(args) => {
                let atoms: Vec<Atom<'a>> = args.iter().map(|a| build_atom(ctx, a)).collect();
                ctx.mul(&atoms)
            }
            PropExpr::Pow(base, exp) => ctx.pow(build_atom(ctx, base), build_atom(ctx, exp)),
        }
    }

    fn prop_expr() -> impl Strategy<Value = PropExpr> {
        let leaf = prop_oneof![
            (-100..100i64).prop_map(PropExpr::Num),
            Just(PropExpr::Var("x")),
            Just(PropExpr::Var("y")),
            Just(PropExpr::Var("z")),
        ];
        leaf.prop_recursive(4, 64, 4, |inner| {
            prop_oneof![
                inner.clone().prop_map(|e| PropExpr::Fun("sin", vec![e])),
                inner.clone().prop_map(|e| PropExpr::Fun("cos", vec![e])),
                prop::collection::vec(inner.clone(), 1..4).prop_map(PropExpr::Add),
                prop::collection::vec(inner.clone(), 1..4).prop_map(PropExpr::Mul),
                (inner.clone(), inner.clone())
                    .prop_map(|(b, e)| PropExpr::Pow(Box::new(b), Box::new(e))),
            ]
        })
    }

    proptest! {
        #[test]
        fn normalize_is_idempotent(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build_atom(&ctx, &expr);
            let once = normalize::normalize(&ctx, atom);
            let twice = normalize::normalize(&ctx, once);
            assert_eq!(once.to_string(), twice.to_string());
        }

        #[test]
        fn add_identity(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build_atom(&ctx, &expr);
            let zero = ctx.num(0);
            let with_zero = ctx.add(&[atom, zero]);
            let normalized = normalize::normalize(&ctx, with_zero);
            assert_eq!(normalized.to_string(), normalize::normalize(&ctx, atom).to_string());
        }

        #[test]
        fn mul_identity(expr in prop_expr()) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = build_atom(&ctx, &expr);
            let one = ctx.num(1);
            let with_one = ctx.mul(&[atom, one]);
            let normalized = normalize::normalize(&ctx, with_one);
            assert_eq!(normalized.to_string(), normalize::normalize(&ctx, atom).to_string());
        }
    }
}
