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

    /// Fold constant subtrees and apply algebraic identities.
    ///
    /// Rules applied:
    /// - `Add`: drop `0` terms, sum all-constant terms, collapse single terms
    /// - `Mul`: absorb `0`, drop `1` factors, multiply all-constant factors,
    ///   collapse single factors
    /// - `Pow`: `x^1 → x`, `x^0 → 1`, `Num^Num` evaluated
    /// - `Fun`: builtin functions of all-constant arguments evaluated
    ///   (external functions are never folded — they may have side effects)
    pub fn fold_constants(&self) -> EvalTree {
        match self {
            EvalTree::Num(_) | EvalTree::Var(_) => self.clone(),
            EvalTree::Add(terms) => {
                let mut folded = Vec::with_capacity(terms.len());
                let mut const_sum = 0.0f64;
                let mut has_const = false;
                for t in terms {
                    match t.fold_constants() {
                        EvalTree::Num(n) => {
                            const_sum += n;
                            has_const = true;
                        }
                        other => folded.push(other),
                    }
                }
                if folded.is_empty() {
                    return EvalTree::Num(const_sum);
                }
                if has_const && const_sum != 0.0 {
                    folded.push(EvalTree::Num(const_sum));
                }
                if folded.len() == 1 {
                    folded.pop().expect("len checked")
                } else {
                    EvalTree::Add(folded)
                }
            }
            EvalTree::Mul(factors) => {
                let mut folded = Vec::with_capacity(factors.len());
                let mut const_prod = 1.0f64;
                let mut has_const = false;
                for f in factors {
                    match f.fold_constants() {
                        EvalTree::Num(n) => {
                            const_prod *= n;
                            has_const = true;
                        }
                        other => folded.push(other),
                    }
                }
                if has_const && const_prod == 0.0 {
                    return EvalTree::Num(0.0);
                }
                if folded.is_empty() {
                    return EvalTree::Num(const_prod);
                }
                if has_const && const_prod != 1.0 {
                    folded.push(EvalTree::Num(const_prod));
                }
                if folded.len() == 1 {
                    folded.pop().expect("len checked")
                } else {
                    EvalTree::Mul(folded)
                }
            }
            EvalTree::Pow(base, exp) => {
                let base = base.fold_constants();
                let exp = exp.fold_constants();
                match (&base, &exp) {
                    (EvalTree::Num(b), EvalTree::Num(e)) => EvalTree::Num(b.powf(*e)),
                    (_, EvalTree::Num(e)) if *e == 0.0 => EvalTree::Num(1.0),
                    (_, EvalTree::Num(e)) if *e == 1.0 => base,
                    _ => EvalTree::Pow(Box::new(base), Box::new(exp)),
                }
            }
            EvalTree::Fun(name, args) => {
                let folded: Vec<EvalTree> = args.iter().map(|a| a.fold_constants()).collect();
                // Fold builtins with all-constant arguments; external
                // functions may have side effects and are never folded.
                if folded.len() == 1
                    && let (Some(op), EvalTree::Num(x)) =
                        (crate::instruction::BuiltinOp::from_name(name), &folded[0])
                {
                    return EvalTree::Num(apply_builtin_f64(op, *x));
                }
                EvalTree::Fun(name.clone(), folded)
            }
        }
    }
}

/// Apply a builtin operation to an f64 value (used by constant folding).
fn apply_builtin_f64(op: crate::instruction::BuiltinOp, x: f64) -> f64 {
    use crate::instruction::BuiltinOp;
    match op {
        BuiltinOp::Sin => x.sin(),
        BuiltinOp::Cos => x.cos(),
        BuiltinOp::Tan => x.tan(),
        BuiltinOp::Sec => 1.0 / x.cos(),
        BuiltinOp::Csc => 1.0 / x.sin(),
        BuiltinOp::Cot => 1.0 / x.tan(),
        BuiltinOp::Exp => x.exp(),
        BuiltinOp::Log => x.ln(),
        BuiltinOp::Sqrt => x.sqrt(),
        BuiltinOp::Abs => x.abs(),
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

    #[test]
    fn fold_add_all_constants() {
        let tree = EvalTree::Add(vec![EvalTree::Num(2.0), EvalTree::Num(3.0)]);
        assert_eq!(tree.fold_constants(), EvalTree::Num(5.0));
    }

    #[test]
    fn fold_add_drops_zero() {
        let tree = EvalTree::Add(vec![EvalTree::Var("x".into()), EvalTree::Num(0.0)]);
        assert_eq!(tree.fold_constants(), EvalTree::Var("x".into()));
    }

    #[test]
    fn fold_add_merges_constants() {
        let tree = EvalTree::Add(vec![
            EvalTree::Var("x".into()),
            EvalTree::Num(2.0),
            EvalTree::Num(3.0),
        ]);
        assert_eq!(
            tree.fold_constants(),
            EvalTree::Add(vec![EvalTree::Var("x".into()), EvalTree::Num(5.0)])
        );
    }

    #[test]
    fn fold_mul_absorbs_zero() {
        let tree = EvalTree::Mul(vec![
            EvalTree::Var("x".into()),
            EvalTree::Num(0.0),
            EvalTree::Var("y".into()),
        ]);
        assert_eq!(tree.fold_constants(), EvalTree::Num(0.0));
    }

    #[test]
    fn fold_mul_drops_one() {
        let tree = EvalTree::Mul(vec![EvalTree::Var("x".into()), EvalTree::Num(1.0)]);
        assert_eq!(tree.fold_constants(), EvalTree::Var("x".into()));
    }

    #[test]
    fn fold_pow_one_and_zero() {
        let x = EvalTree::Var("x".into());
        assert_eq!(
            EvalTree::Pow(Box::new(x.clone()), Box::new(EvalTree::Num(1.0))).fold_constants(),
            x
        );
        assert_eq!(
            EvalTree::Pow(Box::new(x), Box::new(EvalTree::Num(0.0))).fold_constants(),
            EvalTree::Num(1.0)
        );
    }

    #[test]
    fn fold_pow_constants() {
        let tree = EvalTree::Pow(Box::new(EvalTree::Num(2.0)), Box::new(EvalTree::Num(10.0)));
        assert_eq!(tree.fold_constants(), EvalTree::Num(1024.0));
    }

    #[test]
    fn fold_builtin_constant() {
        let tree = EvalTree::Fun("sin".into(), vec![EvalTree::Num(0.0)]);
        assert_eq!(tree.fold_constants(), EvalTree::Num(0.0));
    }

    #[test]
    fn fold_keeps_external_fun() {
        // Unknown (external) functions are not folded even with constant args
        let tree = EvalTree::Fun("my_callback".into(), vec![EvalTree::Num(1.0)]);
        assert_eq!(tree.fold_constants(), tree);
    }

    #[test]
    fn fold_nested() {
        // (2 + 3) * x + 0 → x * 5
        let tree = EvalTree::Add(vec![
            EvalTree::Mul(vec![
                EvalTree::Add(vec![EvalTree::Num(2.0), EvalTree::Num(3.0)]),
                EvalTree::Var("x".into()),
            ]),
            EvalTree::Num(0.0),
        ]);
        assert_eq!(
            tree.fold_constants(),
            EvalTree::Mul(vec![EvalTree::Var("x".into()), EvalTree::Num(5.0)])
        );
    }
}
