//! Multi-pattern replacement for oCAS.
//!
//! Provides `replace_once`, `replace_all`, and `replace_all_multiple` with
//! configurable traversal direction and condition guards.

use ocas_atom::{Atom, AtomArena};

use crate::matcher::{self, Bindings};

/// Settings controlling replacement traversal.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReplaceSettings {
    /// If true, replace exactly once and stop.
    pub once: bool,
    /// If true, traverse bottom-up (children first); otherwise top-down.
    pub bottom_up: bool,
    /// If true, replace inside already-replaced sub-expressions.
    pub nested: bool,
}

/// A condition on match bindings.
#[derive(Clone)]
pub enum Condition<'a> {
    /// Evaluates to true if the closure returns true.
    Predicate(std::sync::Arc<dyn Fn(&Bindings<'a>) -> bool + 'a>),
}

impl<'a> Condition<'a> {
    /// Create a condition from a predicate closure.
    pub fn new<F: Fn(&Bindings<'a>) -> bool + 'a>(f: F) -> Self {
        Condition::Predicate(std::sync::Arc::new(f))
    }

    fn eval(&self, bindings: &Bindings<'a>) -> bool {
        match self {
            Condition::Predicate(f) => f(bindings),
        }
    }
}

/// A single replacement: pattern → replacement with optional condition.
pub struct Replacement<'a, F>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    /// The pattern to match.
    pub pattern: crate::pattern::Pattern<'a>,
    /// Function to produce the replacement atom from bindings and arena.
    pub replacement: F,
    /// Optional condition guard.
    pub condition: Option<Condition<'a>>,
}

/// Replace the first match of `pattern` → `replacement` found in a top-down
/// traversal.  Returns the (possibly unchanged) atom.
pub fn replace_once<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    pattern: crate::pattern::Pattern<'a>,
    replacement: F,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    replace_top_down(ctx, atom, &pattern, &replacement, &None, true)
}

/// Replace all matches of `pattern` → `replacement` found in the atom tree.
/// By default traverses top-down without nested replacement.
pub fn replace_all<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    pattern: crate::pattern::Pattern<'a>,
    replacement: F,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    replace_top_down(ctx, atom, &pattern, &replacement, &None, false)
}

/// Replace all matches using multiple replacements tried in order.
/// First replacement that matches is applied at each node.
pub fn replace_all_multiple<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    replacements: &[Replacement<'a, F>],
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    replace_multiple_top_down(ctx, atom, replacements, false)
}

fn replace_top_down<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    pattern: &crate::pattern::Pattern<'a>,
    replacement: &F,
    condition: &Option<Condition<'a>>,
    once: bool,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    // Try matching this node.
    if let Ok(bindings) = matcher::match_pattern(pattern.clone(), atom)
        && condition.as_ref().is_none_or(|c| c.eval(&bindings))
    {
        return replacement(&bindings, ctx);
    }

    // Recurse into children.
    use ocas_atom::AtomNode;
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => atom,
        AtomNode::Add(args) => {
            let new_args = try_replace_child(args, ctx, pattern, replacement, condition, once);
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args = try_replace_child(args, ctx, pattern, replacement, condition, once);
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = replace_top_down(ctx, *base, pattern, replacement, condition, once);
            let new_exp = replace_top_down(ctx, *exp, pattern, replacement, condition, once);
            ctx.pow(new_base, new_exp)
        }
        AtomNode::Fun(name, args) => {
            let new_args = try_replace_child(args, ctx, pattern, replacement, condition, once);
            ctx.fun(name.as_str(), &new_args)
        }
    }
}

/// Try to replace the first matching child. Returns the (possibly unchanged)
/// child list.  If `once` is true and a child is replaced, subsequent children
/// are kept unchanged.
fn try_replace_child<'a, F>(
    args: &'a [Atom<'a>],
    ctx: &'a AtomArena<'a>,
    pattern: &crate::pattern::Pattern<'a>,
    replacement: &F,
    condition: &Option<Condition<'a>>,
    once: bool,
) -> Vec<Atom<'a>>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    if !once {
        return args
            .iter()
            .map(|a| replace_top_down(ctx, *a, pattern, replacement, condition, false))
            .collect();
    }
    // once mode: scan for the first match.
    let mut found = false;
    args.iter()
        .map(|a| {
            if found {
                *a
            } else {
                let result = replace_top_down(ctx, *a, pattern, replacement, condition, true);
                if result != *a {
                    found = true;
                }
                result
            }
        })
        .collect()
}

fn replace_multiple_top_down<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    replacements: &[Replacement<'a, F>],
    once: bool,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    // Try each replacement in order.
    for repl in replacements {
        if let Ok(bindings) = matcher::match_pattern(repl.pattern.clone(), atom)
            && repl.condition.as_ref().is_none_or(|c| c.eval(&bindings))
        {
            return (repl.replacement)(&bindings, ctx);
        }
    }
    if once {
        return atom;
    }
    use ocas_atom::AtomNode;
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => atom,
        AtomNode::Add(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_multiple_top_down(ctx, *a, replacements, once))
                .collect();
            ctx.add(&new_args)
        }
        AtomNode::Mul(args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_multiple_top_down(ctx, *a, replacements, once))
                .collect();
            ctx.mul(&new_args)
        }
        AtomNode::Pow(base, exp) => {
            let new_base = replace_multiple_top_down(ctx, *base, replacements, once);
            let new_exp = replace_multiple_top_down(ctx, *exp, replacements, once);
            ctx.pow(new_base, new_exp)
        }
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_multiple_top_down(ctx, *a, replacements, once))
                .collect();
            ctx.fun(name.as_str(), &new_args)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pattern::Pattern;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn replace_once_simple() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let pat = Pattern::Literal(x);
        let result = replace_once(&ctx, sum, pat, |_, ctx| ctx.num(42));
        assert_eq!(result.to_string(), "42 + y");
    }

    #[test]
    fn replace_all_replaces_all() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, ctx.add(&[y, x])]);
        let pat = Pattern::Literal(x);
        let result = replace_all(&ctx, sum, pat, |_, ctx| ctx.num(1));
        // Both x's replaced by 1; nested add structure preserved.
        let s = result.to_string();
        assert!(s.contains("1") && !s.contains('x'), "result: {s}");
    }
}
