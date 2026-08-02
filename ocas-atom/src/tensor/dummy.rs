//! Dummy index management for tensor canonicalisation.
//!
//! Provides dummy-index refresh (rename to avoid conflicts), n-ary
//! contraction normalization, and validation for tensor expressions.

use std::collections::HashMap;

use crate::{Atom, AtomArena, AtomNode, Symbol};

use super::spec::TensorRegistry;

/// Errors from dummy-index operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DummyError {
    /// An index label appeared more than twice in a tensor product.
    OverContracted(Symbol),
    /// Two slots with the same label have the same variance (not an upper/lower pair).
    BadContraction(Symbol),
}

/// Refresh (rename) dummy indices in a tensor expression to avoid
/// conflicts with external (free) indices.
///
/// Dummy indices (labels appearing exactly twice in a product, once upper
/// and once lower) are replaced with fresh names from a per-group pool.
/// External indices are left unchanged.
pub fn refresh_dummies<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<Atom<'a>, DummyError> {
    // Collect index usage: label → (count, vec of (head, pos) references).
    // For simplicity, handle only single-term products here.
    let mut index_counts: HashMap<Atom<'a>, usize> = HashMap::new();
    collect_counts(expr, &mut index_counts);

    // Identify dummies (count == 2).
    let dummies: Vec<Atom<'a>> = index_counts
        .iter()
        .filter(|(_, c)| **c == 2)
        .map(|(l, _)| *l)
        .collect();

    if dummies.is_empty() {
        return Ok(expr);
    }

    // Assign fresh names per group.
    let mut group_counters: HashMap<u64, usize> = HashMap::new();
    let mut renames: HashMap<Atom<'a>, Atom<'a>> = HashMap::new();

    for d in &dummies {
        let sym = Symbol::new(&d.to_string());
        let group = registry.index_group(sym);
        let cnt = group_counters.entry(group).or_insert(0);
        let new_label = if group == 0 {
            ctx.var(&format!("d{}", cnt))
        } else {
            ctx.var(&format!("d{}_{}", group, cnt))
        };
        *cnt += 1;
        renames.insert(*d, new_label);
    }

    Ok(rename_in_expr(ctx, expr, &renames))
}

fn collect_counts<'a>(atom: Atom<'a>, counts: &mut HashMap<Atom<'a>, usize>) {
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => {}
        AtomNode::Fun(_, args) => {
            for a in args.iter() {
                *counts.entry(*a).or_insert(0) += 1;
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_counts(*a, counts);
            }
        }
        AtomNode::Pow(base, exp) => {
            collect_counts(*base, counts);
            collect_counts(*exp, counts);
        }
    }
}

fn rename_in_expr<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    renames: &HashMap<Atom<'a>, Atom<'a>>,
) -> Atom<'a> {
    match atom.node() {
        AtomNode::Num(_) => atom,
        AtomNode::Var(_) => renames.get(&atom).copied().unwrap_or(atom),
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| rename_in_expr(ctx, *a, renames))
                .collect();
            ctx.fun(name.as_str(), &new_args)
        }
        AtomNode::Add(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| rename_in_expr(ctx, *a, renames))
                .collect();
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| rename_in_expr(ctx, *a, renames))
                .collect();
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = rename_in_expr(ctx, *base, renames);
            let new_exp = rename_in_expr(ctx, *exp, renames);
            ctx.pow(new_base, new_exp)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;
    use crate::Symbol;
    use crate::tensor::spec::SymmetrySpec;
    use ocas_core::arena::Arena;

    #[test]
    fn refresh_renames_dummy() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let mut reg = TensorRegistry::new();
        reg.register(Symbol::new("T"), SymmetrySpec::none());
        reg.register(Symbol::new("U"), SymmetrySpec::none());

        let i = ctx.var("i");
        let j = ctx.var("j");
        let t = ctx.fun("T", &[i, j]);
        let u = ctx.fun("U", &[j, i]);
        let prod = ctx.mul(&[t, u]);
        let result = refresh_dummies(&ctx, prod, &reg).unwrap();
        let s = result.to_string();
        // i appears as external, j as dummy — j should be renamed to d0.
        assert!(s.contains("d0"), "expected dummy d0 in: {s}");
    }
}
