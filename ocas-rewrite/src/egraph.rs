//! Optional equality saturation integration using `egg`.
//!
//! This module is only available when the `egg` feature is enabled. It provides
//! an [`AtomLanguage`] implementation of [`egg::Language`] and a helper to
//! simplify expressions via equality saturation.

use egg::{define_language, rewrite as rw, AstSize, Extractor, Id, RecExpr, Runner};

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

define_language! {
    /// An `egg` language node for oCAS atoms.
    pub enum AtomLanguage {
        "num" = Num(i64),
        "var" = Var(Symbol),
        "fun" = Fun([Id]),
        "add" = Add([Id]),
        "mul" = Mul([Id]),
        "pow" = Pow([Id; 2]),
    }
}

impl AtomLanguage {
    /// Convert an oCAS [`Atom`] into an `egg` [`RecExpr`].
    ///
    /// `cache` maps already-visited `Atom` nodes to their `Id` so that shared
    /// sub-expressions are represented once in the `RecExpr`.
    pub fn to_recexpr(
        atom: Atom,
        egraph: &mut egg::EGraph<Self, ()>,
        cache: &mut Vec<(Atom, Id)>,
    ) -> Id {
        for (a, id) in cache.iter() {
            if *a == atom {
                return *id;
            }
        }

        let node = match atom.node() {
            AtomNode::Num(n) => AtomLanguage::Num(*n),
            AtomNode::Var(s) => AtomLanguage::Var(*s),
            AtomNode::Fun(name, args) => {
                let mut ids = vec![egraph.add(AtomLanguage::Var(*name))];
                ids.extend(args.iter().map(|a| Self::to_recexpr(*a, egraph, cache)));
                AtomLanguage::Fun(ids.into_boxed_slice())
            }
            AtomNode::Add(args) => {
                let ids: Vec<Id> = args
                    .iter()
                    .map(|a| Self::to_recexpr(*a, egraph, cache))
                    .collect();
                AtomLanguage::Add(ids.into_boxed_slice())
            }
            AtomNode::Mul(args) => {
                let ids: Vec<Id> = args
                    .iter()
                    .map(|a| Self::to_recexpr(*a, egraph, cache))
                    .collect();
                AtomLanguage::Mul(ids.into_boxed_slice())
            }
            AtomNode::Pow(base, exp) => {
                let base_id = Self::to_recexpr(*base, egraph, cache);
                let exp_id = Self::to_recexpr(*exp, egraph, cache);
                AtomLanguage::Pow([base_id, exp_id])
            }
        };

        let id = egraph.add(node);
        cache.push((atom, id));
        id
    }

    /// Convert an `egg` [`RecExpr`] back to an oCAS [`Atom`].
    ///
    /// The `ocas` arena is the one that will own the new atoms. The `egg` arena
    /// is only needed for the mapping from `Id` to `Language` node.
    pub fn from_recexpr<'a>(
        expr: &RecExpr<Self>,
        id: Id,
        ocas_arena: &'a AtomArena<'a>,
    ) -> Atom<'a> {
        let node = expr[id].clone();
        match node {
            AtomLanguage::Num(n) => ocas_arena.num(n),
            AtomLanguage::Var(s) => ocas_arena.var(s.as_str()),
            AtomLanguage::Add(ids) => {
                let args: Vec<Atom> = ids
                    .iter()
                    .map(|i| Self::from_recexpr(expr, *i, ocas_arena))
                    .collect();
                ocas_arena.add(&args)
            }
            AtomLanguage::Mul(ids) => {
                let args: Vec<Atom> = ids
                    .iter()
                    .map(|i| Self::from_recexpr(expr, *i, ocas_arena))
                    .collect();
                ocas_arena.mul(&args)
            }
            AtomLanguage::Pow([base, exp]) => {
                let base_atom = Self::from_recexpr(expr, base, ocas_arena);
                let exp_atom = Self::from_recexpr(expr, exp, ocas_arena);
                ocas_arena.pow(base_atom, exp_atom)
            }
            AtomLanguage::Fun(ids) => {
                let mut iter = ids.iter();
                let head_id = iter.next().expect("fun must have at least head");
                let head = match expr[*head_id] {
                    AtomLanguage::Var(s) => s,
                    _ => panic!("fun head must be a variable"),
                };
                let args: Vec<Atom> = iter
                    .map(|i| Self::from_recexpr(expr, *i, ocas_arena))
                    .collect();
                ocas_arena.fun(head.as_str(), &args)
            }
        }
    }
}

fn rules() -> Vec<egg::Rewrite<AtomLanguage, ()>> {
    vec![
        rw!("add-zero"; "(add 0 ?a)" => "?a"),
        rw!("mul-zero"; "(mul ?a 0)" => "0"),
        rw!("mul-one"; "(mul 1 ?a)" => "?a"),
        rw!("pow-zero"; "(pow ?a 0)" => "1"),
        rw!("pow-one"; "(pow ?a 1)" => "?a"),
        rw!("pythagorean"; "(add (pow (sin ?x) 2) (pow (cos ?x) 2))" => "1"),
    ]
    .unwrap()
}

/// Simplify an expression using equality saturation.
///
/// This is a minimal integration: it converts the expression into an egg
/// e-graph, runs a small set of built-in rules, and extracts the best
/// expression using AST size as the cost function.
pub fn simplify_with_egraph<'a>(
    atom: Atom<'a>,
    ocas_arena: &'a AtomArena<'a>,
    iter_limit: usize,
) -> Atom<'a> {
    let mut egraph = egg::EGraph::<AtomLanguage, ()>::default();
    let mut cache = Vec::new();
    let root = AtomLanguage::to_recexpr(atom, &mut egraph, &mut cache);

    let runner = Runner::default()
        .with_iter_limit(iter_limit)
        .with_egraph(egraph)
        .run(&rules());

    let extractor = Extractor::new(&runner.egraph, AstSize);
    let (_, best_expr) = extractor.find_best(root);
    AtomLanguage::from_recexpr(&best_expr, runner.roots[0], ocas_arena)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn pythagorean_identity() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);

        let x = ctx.var("x");
        let two = ctx.num(2);
        let sin_x = ctx.fun("sin", &[x]);
        let cos_x = ctx.fun("cos", &[x]);
        let sum = ctx.add(&[ctx.pow(sin_x, two), ctx.pow(cos_x, two)]);

        let result = simplify_with_egraph(sum, &ctx, 5);
        assert_eq!(result.to_string(), "1");
    }
}
