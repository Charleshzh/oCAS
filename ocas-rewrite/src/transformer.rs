use ocas_atom::{Atom, AtomArena, AtomNode};

#[cfg(test)]
use ocas_core::arena::Arena;

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
