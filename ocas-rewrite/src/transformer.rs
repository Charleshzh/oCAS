use ocas_atom::{Atom, AtomArena, AtomNode};

#[cfg(test)]
use ocas_core::arena::Arena;

use crate::combinatorics;

/// Partition an `arg(a₁, a₂, …, aₙ)` expression into named bins and return a
/// sum of products `Σ coeff · f₁(…)·f₂(…)·…`.
///
/// Parameters mirror Symbolica's `Transformer::Partition`:
/// * `bins` — list of `(function_name, capacity)`.
/// * `fill_last` — surplus elements absorbed into the last bin.
/// * `repeat` — repeat the bin pattern until all elements consumed.
///
/// Returns `ctx.num(0)` when no valid partition exists.
pub fn partition_expr<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    bins: &[(ocas_atom::Symbol, usize)],
    fill_last: bool,
    repeat: bool,
) -> Atom<'a> {
    // Extract `arg(...)` args.
    let args: &[Atom<'a>] = match expr.node() {
        AtomNode::Fun(name, a) if name.as_str() == "arg" => a,
        _ => return expr,
    };

    if args.is_empty() || bins.is_empty() {
        return ctx.num(0);
    }

    let elements: Vec<i64> = args
        .iter()
        .filter_map(|a| match a.node() {
            AtomNode::Num(n) => Some(*n),
            _ => None,
        })
        .collect();

    // If any arg is not a number, bail.
    if elements.len() != args.len() {
        return expr;
    }

    let bin_specs: Vec<(ocas_atom::Symbol, usize)> = bins.to_vec();
    let sols = combinatorics::partitions(&elements, &bin_specs, fill_last, repeat);

    if sols.is_empty() {
        return ctx.num(0);
    }

    let mut terms: Vec<Atom<'a>> = Vec::new();
    for sol in &sols {
        let coeff_atom = ctx.num(sol.coefficient as i64);
        let mut factors: Vec<Atom<'a>> = vec![coeff_atom];
        for (name, content) in &sol.bins {
            let content_atoms: Vec<Atom<'a>> = content.iter().map(|&n| ctx.num(n)).collect();
            factors.push(ctx.fun(name.as_str(), &content_atoms));
        }
        terms.push(ctx.mul(&factors));
    }

    if terms.is_empty() {
        ctx.num(0)
    } else if terms.len() == 1 {
        terms.pop().unwrap()
    } else {
        ctx.add(&terms)
    }
}

/// Transform an atom tree bottom-up.
///
/// The supplied function `f` is called on each node **after** its children
/// have already been transformed. If `f` returns `Some(atom)`, that atom is
/// used in place of the original; if it returns `None`, the original (with
/// transformed children) is kept.
///
/// This is the standard "rewriting traversal" used by the oCAS rule engine
/// and simplifier.
///
/// # Example
///
/// ```
/// use ocas_atom::{Atom, AtomArena, AtomNode};
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::transformer::transform;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let sum = ctx.add(&[x, y]);
///
/// let result = transform(&ctx, sum, |a| {
///     if let AtomNode::Add(args) = a.node() {
///         if args.len() == 2 && args[0] == x && args[1] == y {
///             return Some(ctx.add(&[y, x]));
///         }
///     }
///     None
/// });
///
/// assert_eq!(result.to_string(), "y + x");
/// ```
fn recurse<'a, F>(ctx: &'a AtomArena<'a>, atom: Atom<'a>, f: &mut F) -> Atom<'a>
where
    F: FnMut(Atom<'a>) -> Option<Atom<'a>>,
{
    let rebuilt = match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => atom,
        AtomNode::Add(args) => {
            let new_args: Vec<Atom<'a>> = args.iter().map(|a| recurse(ctx, *a, f)).collect();
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args: Vec<Atom<'a>> = args.iter().map(|a| recurse(ctx, *a, f)).collect();
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = recurse(ctx, *base, f);
            let new_exp = recurse(ctx, *exp, f);
            ctx.pow(new_base, new_exp)
        }
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom<'a>> = args.iter().map(|a| recurse(ctx, *a, f)).collect();
            ctx.fun(name.as_str(), &new_args)
        }
    };
    f(rebuilt).unwrap_or(rebuilt)
}

pub fn transform<'a, F>(ctx: &'a AtomArena<'a>, atom: Atom<'a>, mut f: F) -> Atom<'a>
where
    F: FnMut(Atom<'a>) -> Option<Atom<'a>>,
{
    recurse(ctx, atom, &mut f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transform_add_children() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let sum = ctx.add(&[x, y, z]);

        let result = transform(&ctx, sum, |a| match a.node() {
            AtomNode::Var(s) if s.as_str() == "x" => Some(ctx.var("a")),
            _ => None,
        });

        assert_eq!(result.to_string(), "a + y + z");
    }

    #[test]
    fn transform_mul_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two = ctx.num(2);
        let three = ctx.num(3);
        let pow = ctx.pow(x, two);
        let prod = ctx.mul(&[pow, three]);

        let result = transform(&ctx, prod, |a| match a.node() {
            AtomNode::Num(2) => Some(ctx.num(7)),
            _ => None,
        });

        assert_eq!(result.to_string(), "(x^7)*3");
    }
}
