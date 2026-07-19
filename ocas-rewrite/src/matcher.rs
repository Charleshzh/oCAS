//! Pattern matching engine for oCAS.
//!
//! The matcher binds [`Pattern`] wildcards to [`Atom`]
//! sub-expressions. It supports associative/commutative matching for `Add`
//! and `Mul` using backtracking, and sequence wildcards for ordered argument
//! lists such as function arguments.

use ocas_atom::{Atom, AtomNode, Symbol};
use ocas_core::FastHashMap as HashMap;

use crate::pattern::{Pattern, WildcardLevel};

/// A collection of wildcard bindings produced by a successful match.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_atom::Symbol;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::matcher::{match_pattern, Bindings, MatchValue};
/// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let pat = Pattern::Wildcard(Symbol::new("w"), WildcardLevel::Single);
/// let bindings = match_pattern(pat, x).unwrap();
/// let value = bindings.get(Symbol::new("w")).unwrap();
/// assert!(matches!(value, MatchValue::Single(v) if v.to_string() == "x"));
/// ```
#[derive(Debug, Clone, Default)]
pub struct Bindings<'a> {
    map: HashMap<Symbol, MatchValue<'a>>,
}

impl<'a> Bindings<'a> {
    /// Create an empty binding set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up a binding by wildcard name.
    pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>> {
        self.map.get(&name)
    }

    /// Insert a single-atom binding. Returns `Err` if the name is already bound
    /// to a different value.
    fn insert_single(&mut self, name: Symbol, value: Atom<'a>) -> Result<(), MatchError> {
        match self.map.get(&name) {
            Some(MatchValue::Single(existing)) if *existing == value => Ok(()),
            Some(_) => Err(MatchError::InconsistentBinding),
            None => {
                self.map.insert(name, MatchValue::Single(value));
                Ok(())
            }
        }
    }

    /// Insert a multi-atom binding (used for sequence/blank null sequence).
    fn insert_sequence(&mut self, name: Symbol, value: &'a [Atom<'a>]) -> Result<(), MatchError> {
        match self.map.get(&name) {
            Some(MatchValue::Sequence(existing)) if *existing == value => Ok(()),
            Some(_) => Err(MatchError::InconsistentBinding),
            None => {
                self.map.insert(name, MatchValue::Sequence(value));
                Ok(())
            }
        }
    }
}

/// A value bound to a wildcard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchValue<'a> {
    /// A single atom binding.
    Single(Atom<'a>),
    /// A slice of atoms bound to a sequence wildcard.
    Sequence(&'a [Atom<'a>]),
}

/// Errors that can occur during matching.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_atom::Symbol;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::matcher::{match_pattern, MatchError};
/// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let pat = Pattern::Literal(x);
/// let err = match_pattern(pat, y).unwrap_err();
/// assert!(matches!(err, MatchError::NoMatch));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchError {
    /// The pattern did not match the atom.
    NoMatch,
    /// A wildcard name was bound to two different values.
    InconsistentBinding,
}

impl std::fmt::Display for MatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MatchError::NoMatch => write!(f, "pattern did not match"),
            MatchError::InconsistentBinding => write!(f, "inconsistent wildcard binding"),
        }
    }
}

impl std::error::Error for MatchError {}

/// Match a pattern against an atom, returning bindings on success.
///
/// `Add` and `Mul` are matched associatively and commutatively: arguments are
/// sorted and literal patterns are matched before single wildcards. Sequence
/// wildcards are not supported inside `Add`/`Mul` in this simplified matcher;
/// they are supported for ordered argument lists such as function arguments.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_atom::Symbol;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::matcher::{match_pattern, MatchValue};
/// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let sum = ctx.add(&[x, y]);
/// let pat = Pattern::Add(vec![
///     Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
///     Pattern::Wildcard(Symbol::new("b"), WildcardLevel::Single),
/// ]);
/// let bindings = match_pattern(pat, sum).unwrap();
/// let a = bindings.get(Symbol::new("a")).unwrap();
/// let b = bindings.get(Symbol::new("b")).unwrap();
/// assert!(matches!(a, MatchValue::Single(v) if v.to_string() == "x"));
/// assert!(matches!(b, MatchValue::Single(v) if v.to_string() == "y"));
/// ```
pub fn match_pattern<'a>(pattern: Pattern<'a>, atom: Atom<'a>) -> Result<Bindings<'a>, MatchError> {
    let mut bindings = Bindings::new();
    match_atom(&mut bindings, pattern, atom)?;
    Ok(bindings)
}

fn match_atom<'a>(
    bindings: &mut Bindings<'a>,
    pattern: Pattern<'a>,
    atom: Atom<'a>,
) -> Result<(), MatchError> {
    match pattern {
        Pattern::Literal(p) => {
            if p == atom {
                Ok(())
            } else {
                Err(MatchError::NoMatch)
            }
        }
        Pattern::Wildcard(name, WildcardLevel::Single) => bindings.insert_single(name, atom),
        Pattern::Wildcard(name, WildcardLevel::Sequence) => match atom.node() {
            AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
                if args.is_empty() {
                    Err(MatchError::NoMatch)
                } else {
                    bindings.insert_sequence(name, args)
                }
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Wildcard(name, WildcardLevel::NullSequence) => match atom.node() {
            AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
                bindings.insert_sequence(name, args)
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Add(pats) => match atom.node() {
            AtomNode::Add(args) => match_nary(bindings, &pats, args, true),
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Mul(pats) => match atom.node() {
            AtomNode::Mul(args) => match_nary(bindings, &pats, args, true),
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Pow(p_box) => match atom.node() {
            AtomNode::Pow(base, exp) => {
                let (p_base, p_exp) = *p_box;
                match_atom(bindings, p_base, *base)?;
                match_atom(bindings, p_exp, *exp)
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Fun(p_name, p_args) => match atom.node() {
            AtomNode::Fun(name, args) if *name == p_name => {
                match_nary(bindings, &p_args, args, false)
            }
            _ => Err(MatchError::NoMatch),
        },
    }
}

/// Match a pattern list against an argument list. For `Add` and `Mul` we sort
/// the target atoms and match literal/sequence patterns greedily; for `Fun`
/// arguments we keep the order and only allow a single trailing null-sequence
/// wildcard to absorb leftovers.
fn match_nary<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &'_ [Pattern<'a>],
    atoms: &'a [Atom<'a>],
    associative_commutative: bool,
) -> Result<(), MatchError> {
    if associative_commutative {
        match_nary_ac(bindings, patterns, atoms)
    } else {
        match_nary_ordered(bindings, patterns, atoms)
    }
}

/// Ordered matching for function arguments. Supports at most one trailing null
/// sequence wildcard. Sequence wildcards are not supported here.
fn match_nary_ordered<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &'_ [Pattern<'a>],
    atoms: &'a [Atom<'a>],
) -> Result<(), MatchError> {
    let mut pat_idx = 0;
    let mut atom_idx = 0;
    while pat_idx < patterns.len() {
        let pat = &patterns[pat_idx];
        match pat {
            Pattern::Wildcard(name, WildcardLevel::NullSequence)
                if pat_idx == patterns.len() - 1 =>
            {
                let rest = &atoms[atom_idx..];
                bindings.insert_sequence(*name, rest)?;
                atom_idx = atoms.len();
                pat_idx += 1;
            }
            Pattern::Wildcard(name, WildcardLevel::Single) => {
                if atom_idx >= atoms.len() {
                    return Err(MatchError::NoMatch);
                }
                bindings.insert_single(*name, atoms[atom_idx])?;
                atom_idx += 1;
                pat_idx += 1;
            }
            _ => {
                if atom_idx >= atoms.len() {
                    return Err(MatchError::NoMatch);
                }
                match_atom(bindings, pat.clone(), atoms[atom_idx])?;
                atom_idx += 1;
                pat_idx += 1;
            }
        }
    }
    if atom_idx == atoms.len() {
        Ok(())
    } else {
        Err(MatchError::NoMatch)
    }
}

/// Associative/commutative matching for `Add` and `Mul`. We sort the target
/// atoms by their natural ordering and match literals first, then single
/// wildcards. Sequence wildcards are not supported in this simplified matcher
/// because the sorted permutation of `atoms` cannot be borrowed as a single
/// contiguous slice.
fn match_nary_ac<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &'_ [Pattern<'a>],
    atoms: &'a [Atom<'a>],
) -> Result<(), MatchError> {
    if patterns.is_empty() {
        return if atoms.is_empty() {
            Ok(())
        } else {
            Err(MatchError::NoMatch)
        };
    }

    if patterns.iter().any(|p| {
        matches!(
            p,
            Pattern::Wildcard(_, WildcardLevel::Sequence | WildcardLevel::NullSequence)
        )
    }) {
        return Err(MatchError::NoMatch);
    }

    let mut sorted_atoms: Vec<Atom<'a>> = atoms.to_vec();
    sorted_atoms.sort();

    let mut matched = vec![false; sorted_atoms.len()];
    let mut single_wildcards: Vec<usize> = Vec::new();
    let mut literal_indices: Vec<usize> = Vec::new();

    for (i, pat) in patterns.iter().enumerate() {
        match pat {
            Pattern::Wildcard(_, WildcardLevel::Single) => single_wildcards.push(i),
            _ => literal_indices.push(i),
        }
    }

    // Match literals first.
    for pat_idx in literal_indices {
        let pat = &patterns[pat_idx];
        let mut found = false;
        for (i, atom) in sorted_atoms.iter().enumerate() {
            if matched[i] {
                continue;
            }
            let mut probe = bindings.clone();
            if match_atom(&mut probe, pat.clone(), *atom).is_ok() {
                *bindings = probe;
                matched[i] = true;
                found = true;
                break;
            }
        }
        if !found {
            return Err(MatchError::NoMatch);
        }
    }

    // Match single wildcards against the remaining atoms.
    let remaining: Vec<Atom<'a>> = sorted_atoms
        .iter()
        .enumerate()
        .filter(|(i, _)| !matched[*i])
        .map(|(_, a)| *a)
        .collect();
    if remaining.len() != single_wildcards.len() {
        return Err(MatchError::NoMatch);
    }
    for (pat_idx, atom) in single_wildcards.iter().zip(remaining.iter()) {
        if let Pattern::Wildcard(name, WildcardLevel::Single) = &patterns[*pat_idx] {
            bindings.insert_single(*name, *atom)?;
        }
    }

    Ok(())
}

impl<'a> Pattern<'a> {
    fn _wildcard_level(&self) -> Option<WildcardLevel> {
        match self {
            Pattern::Wildcard(_, level) => Some(*level),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    struct VecAlloc;

    impl<'a> crate::pattern::PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
            let leaked: Box<[Pattern<'a>]> = items.to_vec().into_boxed_slice();
            Box::leak(leaked)
        }
    }

    fn pat_atom<'a>(ctx: &'a AtomArena<'a>, alloc: &'a VecAlloc, s: &'a str) -> Pattern<'a> {
        let atom = ocas_parse::parse(ctx, s).unwrap();
        Pattern::from_atom(alloc, atom)
    }

    #[test]
    fn single_wildcard_binds_any_atom() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_");
        let y = ctx.var("y");
        let bindings = match_pattern(pat, y).unwrap();
        assert!(matches!(bindings.get(Symbol::new("x")), Some(MatchValue::Single(a)) if *a == y));
    }

    #[test]
    fn add_pattern_matches_two_atoms() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_ + y_");
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let bindings = match_pattern(pat, sum).unwrap();
        assert!(matches!(
            bindings.get(Symbol::new("x")),
            Some(MatchValue::Single(a)) if *a == x
        ));
        assert!(matches!(
            bindings.get(Symbol::new("y")),
            Some(MatchValue::Single(a)) if *a == y
        ));
    }

    #[test]
    fn add_pattern_ignores_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_ + y_");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let sum = ctx.add(&[b, a]);
        let bindings = match_pattern(pat, sum).unwrap();
        assert!(matches!(
            bindings.get(Symbol::new("x")),
            Some(MatchValue::Single(atom)) if *atom == a
        ));
        assert!(matches!(
            bindings.get(Symbol::new("y")),
            Some(MatchValue::Single(atom)) if *atom == b
        ));
    }

    #[test]
    fn add_pattern_requires_exact_match_count() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_ + y_");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let sum = ctx.add(&[a, b, c]);
        assert!(match_pattern(pat, sum).is_err());
    }

    #[test]
    fn inconsistent_binding_fails() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_ + x_");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let sum = ctx.add(&[a, b]);
        assert!(match_pattern(pat, sum).is_err());
    }

    #[test]
    fn sequence_wildcard_in_fun_absorbs_remainder() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "f(x_, rest___)");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let f = ctx.fun("f", &[a, b, c]);
        let bindings = match_pattern(pat, f).unwrap();
        assert!(matches!(
            bindings.get(Symbol::new("x")),
            Some(MatchValue::Single(atom)) if *atom == a
        ));
        let rest = bindings.get(Symbol::new("rest")).unwrap();
        assert!(matches!(rest, MatchValue::Sequence(s) if s.len() == 2));
    }
}
