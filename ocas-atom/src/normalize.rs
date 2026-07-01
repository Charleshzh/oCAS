//! Normalization for [`Atom`] expression trees.
//!
//! The normalizer puts expressions into a deterministic canonical form:
//! nested additions and multiplications are flattened, arguments are sorted,
//! and numeric coefficients are merged.

use crate::{Atom, AtomArena, AtomNode};

/// Normalize an atom into canonical form.
///
/// The result is allocated in the same arena as the input via `ctx`.
///
/// # Example
///
/// ```
/// use ocas_atom::normalize::normalize;
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let z = ctx.var("z");
/// let inner = ctx.add(&[x, y]);
/// let outer = ctx.add(&[inner, z, ctx.num(2), ctx.num(3)]);
/// let result = normalize(&ctx, outer);
/// assert_eq!(result.to_string(), "5 + x + y + z");
/// ```
pub fn normalize<'a>(ctx: &AtomArena<'a>, atom: Atom<'a>) -> Atom<'a> {
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => atom,
        AtomNode::Fun(name, args) => {
            let mut normalized: Vec<Atom<'a>> = args.iter().map(|a| normalize(ctx, *a)).collect();
            // Preserve argument order for calculus forms where order is semantic.
            if !matches!(name.as_str(), "Derivative" | "Integral") {
                normalized.sort();
            }
            ctx.fun(name.as_str(), &normalized)
        }
        AtomNode::Add(args) => {
            let mut flat = Vec::new();
            collect_add(args, &mut flat);
            let mut normalized: Vec<Atom<'a>> =
                flat.into_iter().map(|a| normalize(ctx, a)).collect();
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(0)));
            normalized.sort();
            merge_numbers(ctx, &mut normalized, true);
            if normalized.is_empty() {
                ctx.num(0)
            } else if normalized.len() == 1 {
                normalized[0]
            } else {
                ctx.add(&normalized)
            }
        }
        AtomNode::Mul(args) => {
            let mut flat = Vec::new();
            collect_mul(args, &mut flat);
            let mut normalized: Vec<Atom<'a>> =
                flat.into_iter().map(|a| normalize(ctx, a)).collect();
            if normalized
                .iter()
                .any(|a| matches!(a.node(), AtomNode::Num(0)))
            {
                return ctx.num(0);
            }
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(1)));
            normalized.sort();
            merge_numbers(ctx, &mut normalized, false);
            if normalized.is_empty() {
                ctx.num(1)
            } else if normalized.len() == 1 {
                normalized[0]
            } else {
                ctx.mul(&normalized)
            }
        }
        AtomNode::Pow(base, exp) => {
            let base = normalize(ctx, *base);
            let exp = normalize(ctx, *exp);
            ctx.pow(base, exp)
        }
    }
}

fn collect_add<'a>(args: &[Atom<'a>], out: &mut Vec<Atom<'a>>) {
    for &arg in args {
        match arg.node() {
            AtomNode::Add(inner) => collect_add(inner, out),
            _ => out.push(arg),
        }
    }
}

fn collect_mul<'a>(args: &[Atom<'a>], out: &mut Vec<Atom<'a>>) {
    for &arg in args {
        match arg.node() {
            AtomNode::Mul(inner) => collect_mul(inner, out),
            _ => out.push(arg),
        }
    }
}

fn merge_numbers<'a>(ctx: &AtomArena<'a>, args: &mut Vec<Atom<'a>>, is_add: bool) {
    let count = args
        .iter()
        .take_while(|a| matches!(a.node(), AtomNode::Num(_)))
        .count();

    if count >= 2 {
        let merged = if is_add {
            args[0..count]
                .iter()
                .map(|a| match a.node() {
                    AtomNode::Num(n) => *n,
                    _ => unreachable!(),
                })
                .sum()
        } else {
            args[0..count]
                .iter()
                .map(|a| match a.node() {
                    AtomNode::Num(n) => *n,
                    _ => unreachable!(),
                })
                .product()
        };
        args.drain(0..count);
        args.insert(0, ctx.num(merged));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn normalize_leaves_atom_unchanged() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert_eq!(normalize(&ctx, x).to_string(), "x");
    }

    #[test]
    fn normalize_flattens_nested_add() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let inner = ctx.add(&[x, y]);
        let outer = ctx.add(&[inner, z]);
        assert_eq!(normalize(&ctx, outer).to_string(), "x + y + z");
    }

    #[test]
    fn normalize_sorts_arguments() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let expr = ctx.add(&[z, x, y]);
        assert_eq!(normalize(&ctx, expr).to_string(), "x + y + z");
    }

    #[test]
    fn normalize_merges_numeric_literals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let one = ctx.num(1);
        let two = ctx.num(2);
        let x = ctx.var("x");
        let expr = ctx.add(&[one, x, two]);
        assert_eq!(normalize(&ctx, expr).to_string(), "3 + x");
    }

    #[test]
    fn normalize_pow() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two = ctx.num(2);
        let pow = ctx.pow(x, two);
        assert_eq!(normalize(&ctx, pow).to_string(), "x^2");
    }

    #[test]
    fn normalize_sorts_fun_arguments() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let f = ctx.fun("f", &[y, x]);
        assert_eq!(normalize(&ctx, f).to_string(), "f(x, y)");
    }
}
