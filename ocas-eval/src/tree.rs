//! Owned intermediate representation for expression tree compilation.
//!
//! [`EvalTree`] is an arena-free, owned expression tree that can be
//! constructed from an [`Atom`](ocas_atom::Atom) and then optimized
//! before instruction generation. It decouples the compilation pipeline
//! from the arena lifetime.

use ocas_atom::{Atom, AtomNode};

/// An owned intermediate representation of a symbolic expression.
///
/// Unlike [`Atom`](ocas_atom::Atom), `EvalTree` owns all its data and
/// does not depend on an arena. This makes it suitable for multi-pass
/// compilation and optimization.
#[derive(Debug, Clone, PartialEq)]
pub enum EvalTree {
    /// A numeric constant.
    Num(f64),
    /// A variable reference.
    Var(String),
    /// A named function applied to arguments.
    Fun(String, Vec<EvalTree>),
    /// A sum of terms.
    Add(Vec<EvalTree>),
    /// A product of factors.
    Mul(Vec<EvalTree>),
    /// A power expression: base^exponent.
    Pow(Box<EvalTree>, Box<EvalTree>),
}

impl EvalTree {
    /// Convert an [`Atom`] into an owned `EvalTree`.
    ///
    /// This dissociates the expression from the arena, allowing the
    /// compilation pipeline to work without lifetime constraints.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use ocas_atom::AtomArena;
    /// use ocas_eval::EvalTree;
    ///
    /// let arena = ocas_core::arena::Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let atom = ctx.add(&[ctx.var("x"), ctx.num(2)]);
    /// let tree = EvalTree::from_atom(atom);
    /// ```
    pub fn from_atom(atom: Atom<'_>) -> Self {
        match atom.node() {
            AtomNode::Num(n) => EvalTree::Num(*n as f64),
            AtomNode::Var(s) => EvalTree::Var(s.as_str().to_string()),
            AtomNode::Fun(s, args) => {
                let converted: Vec<EvalTree> =
                    args.iter().map(|a| EvalTree::from_atom(*a)).collect();
                EvalTree::Fun(s.as_str().to_string(), converted)
            }
            AtomNode::Add(terms) => {
                let converted: Vec<EvalTree> =
                    terms.iter().map(|a| EvalTree::from_atom(*a)).collect();
                EvalTree::Add(converted)
            }
            AtomNode::Mul(factors) => {
                let converted: Vec<EvalTree> =
                    factors.iter().map(|a| EvalTree::from_atom(*a)).collect();
                EvalTree::Mul(converted)
            }
            AtomNode::Pow(base, exp) => EvalTree::Pow(
                Box::new(EvalTree::from_atom(*base)),
                Box::new(EvalTree::from_atom(*exp)),
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn from_atom_num() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = ctx.num(42);
        let tree = EvalTree::from_atom(atom);
        assert_eq!(tree, EvalTree::Num(42.0));
    }

    #[test]
    fn from_atom_var() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = ctx.var("x");
        let tree = EvalTree::from_atom(atom);
        assert_eq!(tree, EvalTree::Var("x".into()));
    }

    #[test]
    fn from_atom_add() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = ctx.add(&[ctx.var("x"), ctx.num(2)]);
        let tree = EvalTree::from_atom(atom);
        assert_eq!(
            tree,
            EvalTree::Add(vec![EvalTree::Var("x".into()), EvalTree::Num(2.0)])
        );
    }

    #[test]
    fn from_atom_pow() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let atom = ctx.pow(ctx.var("x"), ctx.num(2));
        let tree = EvalTree::from_atom(atom);
        assert_eq!(
            tree,
            EvalTree::Pow(
                Box::new(EvalTree::Var("x".into())),
                Box::new(EvalTree::Num(2.0))
            )
        );
    }
}
