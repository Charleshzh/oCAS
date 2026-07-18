//! Read-only traversal utilities over expression trees.
//!
//! These helpers walk an [`Atom`] tree in post-order and collect subtrees
//! such as function applications or variables. They are used, for example,
//! by the calculus tower construction to discover which elementary
//! functions occur in an integrand.

use crate::{Atom, AtomNode, Symbol};

/// Collect all function applications occurring in `atom`.
///
/// Functions are returned in post-order (innermost applications first) with
/// duplicates removed — since atoms are hash-consed, structurally equal
/// applications appear only once.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, walk::collect_funs};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let sin_x = ctx.fun("sin", &[x]);
/// let expr = ctx.fun("cos", &[sin_x]);
/// let funs = collect_funs(expr);
/// assert_eq!(funs.len(), 2);
/// assert_eq!(funs[0].0.as_str(), "sin");
/// assert_eq!(funs[1].0.as_str(), "cos");
/// ```
pub fn collect_funs<'a>(atom: Atom<'a>) -> Vec<(Symbol, Atom<'a>)> {
    let mut out = Vec::new();
    collect_funs_into(atom, &mut out);
    out
}

fn collect_funs_into<'a>(atom: Atom<'a>, out: &mut Vec<(Symbol, Atom<'a>)>) {
    for child in atom.children() {
        collect_funs_into(*child, out);
    }
    if let Some((base, exp)) = atom.binary_children() {
        collect_funs_into(base, out);
        collect_funs_into(exp, out);
    }
    if let AtomNode::Fun(name, _) = atom.node()
        && !out.iter().any(|&(_, existing)| existing == atom)
    {
        out.push((*name, atom));
    }
}

/// Collect the distinct variable names occurring in `atom`.
///
/// Variables are returned in order of first appearance (depth-first,
/// left-to-right).
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, walk::collect_vars};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let expr = ctx.add(&[ctx.mul(&[x, y]), x]);
/// let vars = collect_vars(expr);
/// assert_eq!(vars.len(), 2);
/// assert_eq!(vars[0].as_str(), "x");
/// assert_eq!(vars[1].as_str(), "y");
/// ```
pub fn collect_vars(atom: Atom) -> Vec<Symbol> {
    let mut out = Vec::new();
    collect_vars_into(atom, &mut out);
    out
}

fn collect_vars_into(atom: Atom, out: &mut Vec<Symbol>) {
    for child in atom.children() {
        collect_vars_into(*child, out);
    }
    if let Some((base, exp)) = atom.binary_children() {
        collect_vars_into(base, out);
        collect_vars_into(exp, out);
    }
    if let AtomNode::Var(name) = atom.node()
        && !out.contains(name)
    {
        out.push(*name);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn collect_funs_post_order_and_dedup() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sin_x = ctx.fun("sin", &[x]);
        // sin(x) + cos(sin(x)) + sin(x)
        let expr = ctx.add(&[sin_x, ctx.fun("cos", &[sin_x]), sin_x]);
        let funs = collect_funs(expr);
        assert_eq!(funs.len(), 2);
        assert_eq!(funs[0].0.as_str(), "sin");
        assert_eq!(funs[0].1, sin_x);
        assert_eq!(funs[1].0.as_str(), "cos");
    }

    #[test]
    fn collect_funs_reaches_pow_operands() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // x^sin(x): the exponent is hidden behind Pow's binary children.
        let expr = ctx.pow(x, ctx.fun("sin", &[x]));
        let funs = collect_funs(expr);
        assert_eq!(funs.len(), 1);
        assert_eq!(funs[0].0.as_str(), "sin");
    }

    #[test]
    fn collect_funs_plain_polynomial_is_empty() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), x, ctx.num(1)]);
        assert!(collect_funs(expr).is_empty());
    }

    #[test]
    fn collect_vars_first_appearance_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let expr = ctx.add(&[ctx.mul(&[y, x]), ctx.pow(y, ctx.num(2))]);
        let vars = collect_vars(expr);
        assert_eq!(vars.len(), 2);
        assert_eq!(vars[0].as_str(), "y");
        assert_eq!(vars[1].as_str(), "x");
    }

    #[test]
    fn collect_vars_ignores_function_names() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("log", &[x]);
        let vars = collect_vars(expr);
        assert_eq!(vars.len(), 1);
        assert_eq!(vars[0].as_str(), "x");
    }
}
