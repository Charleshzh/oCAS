//! Pattern matching engine for oCAS.
//!
//! The matcher binds [`Pattern`] wildcards to [`Atom`] sub-expressions.
//! Associative/commutative matching for `Add`/`Mul` uses full backtracking
//! search with a budget to prevent pathological explosion.  Sequence wildcards
//! are supported in all argument-list contexts (`Add`, `Mul`, `Fun`).

use ocas_atom::{Atom, AtomNode, Symbol};
use ocas_core::FastHashMap as HashMap;

use crate::pattern::{Pattern, WildcardLevel};

/// Default maximum number of backtrack attempts per AC match.
pub const DEFAULT_MAX_BACKTRACKS: usize = 10_000;

/// A collection of wildcard bindings produced by a successful match.
#[derive(Debug, Clone, Default)]
pub struct Bindings<'a> {
    map: HashMap<Symbol, MatchValue<'a>>,
}

impl<'a> Bindings<'a> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>> {
        self.map.get(&name)
    }
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
    fn remove(&mut self, name: Symbol) {
        self.map.remove(&name);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchValue<'a> {
    Single(Atom<'a>),
    Sequence(&'a [Atom<'a>]),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchError {
    NoMatch,
    InconsistentBinding,
    BudgetExhausted,
}

impl std::fmt::Display for MatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MatchError::NoMatch => write!(f, "pattern did not match"),
            MatchError::InconsistentBinding => write!(f, "inconsistent wildcard binding"),
            MatchError::BudgetExhausted => write!(f, "backtrack budget exhausted"),
        }
    }
}
impl std::error::Error for MatchError {}

pub fn match_pattern<'a>(pattern: Pattern<'a>, atom: Atom<'a>) -> Result<Bindings<'a>, MatchError> {
    match_pattern_with_budget(pattern, atom, DEFAULT_MAX_BACKTRACKS)
}

pub fn match_pattern_with_budget<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
    max_backtracks: usize,
) -> Result<Bindings<'a>, MatchError> {
    let mut bindings = Bindings::new();
    let mut backtrack_count = 0usize;
    match_atom(
        &mut bindings,
        pattern,
        atom,
        &mut backtrack_count,
        max_backtracks,
    )?;
    Ok(bindings)
}

fn match_atom<'a>(
    bindings: &mut Bindings<'a>,
    pattern: Pattern<'a>,
    atom: Atom<'a>,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> Result<(), MatchError> {
    if *backtrack_count >= max_backtracks {
        return Err(MatchError::BudgetExhausted);
    }
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
            AtomNode::Add(args) => {
                match_nary(bindings, &pats, args, true, backtrack_count, max_backtracks)
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Mul(pats) => match atom.node() {
            AtomNode::Mul(args) => {
                match_nary(bindings, &pats, args, true, backtrack_count, max_backtracks)
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Pow(p_box) => match atom.node() {
            AtomNode::Pow(base, exp) => {
                let (p_base, p_exp) = *p_box;
                match_atom(bindings, p_base, *base, backtrack_count, max_backtracks)?;
                match_atom(bindings, p_exp, *exp, backtrack_count, max_backtracks)
            }
            _ => Err(MatchError::NoMatch),
        },
        Pattern::Fun(p_name, p_args) => match atom.node() {
            AtomNode::Fun(name, args) if *name == p_name => match_nary(
                bindings,
                &p_args,
                args,
                false,
                backtrack_count,
                max_backtracks,
            ),
            _ => Err(MatchError::NoMatch),
        },
    }
}

fn match_nary<'a>(
    bindings: &'_ mut Bindings<'a>,
    patterns: &[Pattern<'a>],
    atoms: &'a [Atom<'a>],
    associative_commutative: bool,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> Result<(), MatchError> {
    if *backtrack_count >= max_backtracks {
        return Err(MatchError::BudgetExhausted);
    }
    if patterns.is_empty() {
        return if atoms.is_empty() {
            Ok(())
        } else {
            Err(MatchError::NoMatch)
        };
    }
    if associative_commutative {
        let mut sorted: Vec<Atom<'a>> = atoms.to_vec();
        sorted.sort();
        let mut used = vec![false; sorted.len()];
        match_nary_ac(
            bindings,
            patterns,
            &sorted,
            &mut used,
            0,
            backtrack_count,
            max_backtracks,
        )
    } else {
        match_nary_ordered(
            bindings,
            patterns,
            atoms,
            0,
            0,
            backtrack_count,
            max_backtracks,
        )
    }
}

// ---- ordered matching (Fun args) with sequence-wildcard support at any position ----

fn match_nary_ordered<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &[Pattern<'a>],
    atoms: &'a [Atom<'a>],
    pat_idx: usize,
    atom_idx: usize,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> Result<(), MatchError> {
    if *backtrack_count >= max_backtracks {
        return Err(MatchError::BudgetExhausted);
    }
    if pat_idx >= patterns.len() {
        return if atom_idx >= atoms.len() {
            Ok(())
        } else {
            Err(MatchError::NoMatch)
        };
    }
    let pat = &patterns[pat_idx];
    match pat {
        Pattern::Wildcard(name, WildcardLevel::NullSequence) => {
            let remaining = atoms.len().saturating_sub(atom_idx);
            for len in 0..=remaining {
                let slice = &atoms[atom_idx..atom_idx + len];
                let mut probe = bindings.clone();
                if probe.insert_sequence(*name, slice).is_err() {
                    continue;
                }
                *backtrack_count += 1;
                match match_nary_ordered(
                    &mut probe,
                    patterns,
                    atoms,
                    pat_idx + 1,
                    atom_idx + len,
                    backtrack_count,
                    max_backtracks,
                ) {
                    Ok(()) => {
                        *bindings = probe;
                        return Ok(());
                    }
                    Err(MatchError::BudgetExhausted) => return Err(MatchError::BudgetExhausted),
                    Err(_) => {}
                }
            }
            Err(MatchError::NoMatch)
        }
        Pattern::Wildcard(name, WildcardLevel::Sequence) => {
            let remaining = atoms.len().saturating_sub(atom_idx);
            if remaining == 0 {
                return Err(MatchError::NoMatch);
            }
            for len in 1..=remaining {
                let slice = &atoms[atom_idx..atom_idx + len];
                let mut probe = bindings.clone();
                if probe.insert_sequence(*name, slice).is_err() {
                    continue;
                }
                *backtrack_count += 1;
                match match_nary_ordered(
                    &mut probe,
                    patterns,
                    atoms,
                    pat_idx + 1,
                    atom_idx + len,
                    backtrack_count,
                    max_backtracks,
                ) {
                    Ok(()) => {
                        *bindings = probe;
                        return Ok(());
                    }
                    Err(MatchError::BudgetExhausted) => return Err(MatchError::BudgetExhausted),
                    Err(_) => {}
                }
            }
            Err(MatchError::NoMatch)
        }
        Pattern::Wildcard(name, WildcardLevel::Single) => {
            if atom_idx >= atoms.len() {
                return Err(MatchError::NoMatch);
            }
            let atom = atoms[atom_idx];
            let mut probe = bindings.clone();
            if probe.insert_single(*name, atom).is_err() {
                return Err(MatchError::NoMatch);
            }
            *backtrack_count += 1;
            match_nary_ordered(
                &mut probe,
                patterns,
                atoms,
                pat_idx + 1,
                atom_idx + 1,
                backtrack_count,
                max_backtracks,
            )
            .map(|()| {
                *bindings = probe;
            })
        }
        _ => {
            if atom_idx >= atoms.len() {
                return Err(MatchError::NoMatch);
            }
            let mut probe = bindings.clone();
            match_atom(
                &mut probe,
                pat.clone(),
                atoms[atom_idx],
                backtrack_count,
                max_backtracks,
            )?;
            *backtrack_count += 1;
            match_nary_ordered(
                &mut probe,
                patterns,
                atoms,
                pat_idx + 1,
                atom_idx + 1,
                backtrack_count,
                max_backtracks,
            )
            .map(|()| {
                *bindings = probe;
            })
        }
    }
}

// ---- AC matching (Add/Mul) — backtracking over sorted atoms ----

fn match_nary_ac<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &[Pattern<'a>],
    sorted_atoms: &[Atom<'a>],
    used: &mut [bool],
    pat_idx: usize,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> Result<(), MatchError> {
    if *backtrack_count >= max_backtracks {
        return Err(MatchError::BudgetExhausted);
    }
    if pat_idx >= patterns.len() {
        return if used.iter().all(|&u| u) {
            Ok(())
        } else {
            Err(MatchError::NoMatch)
        };
    }

    let pat = &patterns[pat_idx];
    match pat {
        Pattern::Literal(v) => {
            for i in 0..sorted_atoms.len() {
                if used[i] || sorted_atoms[i] != *v {
                    continue;
                }
                used[i] = true;
                *backtrack_count += 1;
                match match_nary_ac(
                    bindings,
                    patterns,
                    sorted_atoms,
                    used,
                    pat_idx + 1,
                    backtrack_count,
                    max_backtracks,
                ) {
                    Ok(()) => return Ok(()),
                    Err(MatchError::BudgetExhausted) => return Err(MatchError::BudgetExhausted),
                    Err(_) => {
                        used[i] = false;
                    }
                }
            }
            Err(MatchError::NoMatch)
        }
        Pattern::Wildcard(name, WildcardLevel::Single) => {
            for i in 0..sorted_atoms.len() {
                if used[i] {
                    continue;
                }
                let atom = sorted_atoms[i];
                match bindings.get(*name).copied() {
                    Some(MatchValue::Single(existing)) if existing == atom => {
                        used[i] = true;
                        match match_nary_ac(
                            bindings,
                            patterns,
                            sorted_atoms,
                            used,
                            pat_idx + 1,
                            backtrack_count,
                            max_backtracks,
                        ) {
                            Ok(()) => return Ok(()),
                            Err(MatchError::BudgetExhausted) => {
                                used[i] = false;
                                return Err(MatchError::BudgetExhausted);
                            }
                            Err(_) => {
                                used[i] = false;
                            }
                        }
                    }
                    Some(_) => continue,
                    None => {
                        used[i] = true;
                        let _ = bindings.insert_single(*name, atom);
                        *backtrack_count += 1;
                        match match_nary_ac(
                            bindings,
                            patterns,
                            sorted_atoms,
                            used,
                            pat_idx + 1,
                            backtrack_count,
                            max_backtracks,
                        ) {
                            Ok(()) => return Ok(()),
                            Err(MatchError::BudgetExhausted) => {
                                bindings.remove(*name);
                                used[i] = false;
                                return Err(MatchError::BudgetExhausted);
                            }
                            Err(_) => {
                                bindings.remove(*name);
                                used[i] = false;
                            }
                        }
                    }
                }
            }
            Err(MatchError::NoMatch)
        }
        Pattern::Wildcard(name, WildcardLevel::NullSequence) => {
            let free: Vec<usize> = (0..sorted_atoms.len()).filter(|&i| !used[i]).collect();
            enumerate_subsets(
                bindings,
                patterns,
                sorted_atoms,
                used,
                pat_idx,
                *name,
                &free,
                0,
                free.len(),
                0,
                true,
                backtrack_count,
                max_backtracks,
            )
        }
        Pattern::Wildcard(name, WildcardLevel::Sequence) => {
            let free: Vec<usize> = (0..sorted_atoms.len()).filter(|&i| !used[i]).collect();
            if free.is_empty() {
                return Err(MatchError::NoMatch);
            }
            enumerate_subsets(
                bindings,
                patterns,
                sorted_atoms,
                used,
                pat_idx,
                *name,
                &free,
                0,
                free.len(),
                1,
                false,
                backtrack_count,
                max_backtracks,
            )
        }
        _ => {
            for i in 0..sorted_atoms.len() {
                if used[i] {
                    continue;
                }
                let mut probe = bindings.clone();
                if let Ok(()) = match_atom(
                    &mut probe,
                    pat.clone(),
                    sorted_atoms[i],
                    backtrack_count,
                    max_backtracks,
                ) {
                    used[i] = true;
                    *backtrack_count += 1;
                    match match_nary_ac(
                        &mut probe,
                        patterns,
                        sorted_atoms,
                        used,
                        pat_idx + 1,
                        backtrack_count,
                        max_backtracks,
                    ) {
                        Ok(()) => {
                            *bindings = probe;
                            return Ok(());
                        }
                        Err(MatchError::BudgetExhausted) => {
                            used[i] = false;
                            return Err(MatchError::BudgetExhausted);
                        }
                        Err(_) => {
                            used[i] = false;
                        }
                    }
                }
            }
            Err(MatchError::NoMatch)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn enumerate_subsets<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &[Pattern<'a>],
    sorted_atoms: &[Atom<'a>],
    used: &mut [bool],
    pat_idx: usize,
    name: Symbol,
    free: &[usize],
    start: usize,
    total_free: usize,
    min_size: usize,
    allow_empty: bool,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> Result<(), MatchError> {
    if *backtrack_count >= max_backtracks {
        return Err(MatchError::BudgetExhausted);
    }
    let rem_pats = patterns.len().saturating_sub(pat_idx + 1);
    let max_size = total_free
        .saturating_sub(start)
        .saturating_sub(rem_pats)
        .min(total_free.saturating_sub(start));
    let min = if allow_empty { 0 } else { min_size };
    for size in min..=max_size {
        let mut chosen: Vec<usize> = Vec::with_capacity(size);
        let result = enumerate_combinations_step(
            bindings,
            patterns,
            sorted_atoms,
            used,
            pat_idx,
            name,
            free,
            start,
            size,
            0,
            &mut chosen,
            backtrack_count,
            max_backtracks,
        );
        match result {
            SubsetResult::Found => return Ok(()),
            SubsetResult::BudgetExhausted => return Err(MatchError::BudgetExhausted),
            SubsetResult::NoMatch => {}
        }
    }
    Err(MatchError::NoMatch)
}

enum SubsetResult {
    Found,
    NoMatch,
    BudgetExhausted,
}

#[allow(clippy::too_many_arguments)]
fn enumerate_combinations_step<'a>(
    bindings: &mut Bindings<'a>,
    patterns: &[Pattern<'a>],
    sorted_atoms: &[Atom<'a>],
    used: &mut [bool],
    pat_idx: usize,
    name: Symbol,
    free: &[usize],
    start: usize,
    remaining: usize,
    _depth: usize,
    chosen: &mut Vec<usize>,
    backtrack_count: &mut usize,
    max_backtracks: usize,
) -> SubsetResult {
    if *backtrack_count >= max_backtracks {
        return SubsetResult::BudgetExhausted;
    }
    if remaining == 0 {
        for &idx in chosen.iter() {
            used[idx] = true;
        }
        let slice: Vec<Atom<'a>> = chosen.iter().map(|&i| sorted_atoms[i]).collect();
        let leaked: &'a [Atom<'a>] = Vec::leak(slice);
        let mut probe = bindings.clone();
        if let Ok(()) = probe.insert_sequence(name, leaked) {
            *backtrack_count += 1;
            match match_nary_ac(
                &mut probe,
                patterns,
                sorted_atoms,
                used,
                pat_idx + 1,
                backtrack_count,
                max_backtracks,
            ) {
                Ok(()) => {
                    *bindings = probe;
                    return SubsetResult::Found;
                }
                Err(MatchError::BudgetExhausted) => {
                    for &idx in chosen.iter() {
                        used[idx] = false;
                    }
                    return SubsetResult::BudgetExhausted;
                }
                Err(_) => {}
            }
        }
        for &idx in chosen.iter() {
            used[idx] = false;
        }
        return SubsetResult::NoMatch;
    }
    let needed = remaining;
    let available = free.len().saturating_sub(start);
    if available < needed {
        return SubsetResult::NoMatch;
    }
    for i in start..=free.len().saturating_sub(needed) {
        chosen.push(free[i]);
        match enumerate_combinations_step(
            bindings,
            patterns,
            sorted_atoms,
            used,
            pat_idx,
            name,
            free,
            i + 1,
            remaining - 1,
            _depth + 1,
            chosen,
            backtrack_count,
            max_backtracks,
        ) {
            SubsetResult::Found => return SubsetResult::Found,
            SubsetResult::BudgetExhausted => {
                chosen.pop();
                return SubsetResult::BudgetExhausted;
            }
            SubsetResult::NoMatch => {}
        }
        chosen.pop();
    }
    SubsetResult::NoMatch
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    struct VecAlloc;
    impl<'a> crate::pattern::PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
            Box::leak(items.to_vec().into_boxed_slice())
        }
    }

    fn pat_expr<'a>(ctx: &'a AtomArena<'a>, _alloc: &'a VecAlloc, s: &'a str) -> Pattern<'a> {
        use ocas_parse;
        let atom = ocas_parse::parse(ctx, s).expect("parse");
        Pattern::from_atom(&(), atom)
    }

    #[test]
    fn match_single_wildcard() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let pat = Pattern::Wildcard(Symbol::new("w"), WildcardLevel::Single);
        let bindings = match_pattern(pat, x).unwrap();
        assert!(matches!(bindings.get(Symbol::new("w")), Some(MatchValue::Single(v)) if *v == x));
    }

    #[test]
    fn match_add_two_singles() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("b"), WildcardLevel::Single),
        ]);
        let bindings = match_pattern(pat, sum).unwrap();
        assert!(matches!(bindings.get(Symbol::new("a")), Some(MatchValue::Single(v)) if *v == x));
        assert!(matches!(bindings.get(Symbol::new("b")), Some(MatchValue::Single(v)) if *v == y));
    }

    #[test]
    fn ac_match_with_literal() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
            Pattern::Literal(y),
        ]);
        let bindings = match_pattern(pat, sum).unwrap();
        assert!(matches!(bindings.get(Symbol::new("a")), Some(MatchValue::Single(v)) if *v == x));
    }

    #[test]
    fn ac_mismatched_count_fails() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let sum = ctx.add(&[x, y, z]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("b"), WildcardLevel::Single),
        ]);
        assert!(matches!(match_pattern(pat, sum), Err(MatchError::NoMatch)));
    }

    #[test]
    fn inconsistent_binding_fails() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let sum = ctx.add(&[x, y]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("w"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("w"), WildcardLevel::Single),
        ]);
        assert!(matches!(
            match_pattern(pat, sum),
            Err(MatchError::InconsistentBinding | MatchError::NoMatch)
        ));
    }

    #[test]
    fn fun_trailing_null_sequence() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let f = ctx.fun("f", &[x, y, z]);
        let alloc = VecAlloc;
        let pattern_str = "f(x_, ___rest)";
        let pat = pat_expr(&ctx, &alloc, pattern_str);
        let bindings = match_pattern(pat, f).unwrap();
        let rest = bindings.get(Symbol::new("rest")).unwrap();
        assert!(matches!(rest, MatchValue::Sequence(s) if s.len() == 2));
    }

    #[test]
    fn ac_three_wildcards_four_terms() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let w = ctx.var("w");
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let sum = ctx.add(&[w, x, y, z]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("b"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("c"), WildcardLevel::Single),
        ]);
        assert!(matches!(match_pattern(pat, sum), Err(MatchError::NoMatch)));
    }

    #[test]
    fn ac_sequence_wildcard_in_add() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let sum = ctx.add(&[a, b, c]);
        let pat = Pattern::Add(vec![
            Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single),
            Pattern::Wildcard(Symbol::new("rest"), WildcardLevel::Sequence),
        ]);
        let bindings = match_pattern(pat, sum).unwrap();
        let rest = bindings.get(Symbol::new("rest")).unwrap();
        assert!(matches!(rest, MatchValue::Sequence(s) if s.len() == 2));
    }

    #[test]
    fn ac_null_sequence_consumes_all() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let sum = ctx.add(&[a]);
        let pat = Pattern::Add(vec![Pattern::Wildcard(
            Symbol::new("rest"),
            WildcardLevel::NullSequence,
        )]);
        let bindings = match_pattern(pat, sum).unwrap();
        let rest = bindings.get(Symbol::new("rest")).unwrap();
        assert!(matches!(rest, MatchValue::Sequence(s) if s.len() == 1));
    }

    #[test]
    fn ordered_sequence_mid_function() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let a = ctx.var("a");
        let b = ctx.var("b");
        let c = ctx.var("c");
        let f = ctx.fun("f", &[a, b, c]);
        let alloc = VecAlloc;
        let pattern_str = "f(x_, __mid, z_)";
        let pat = pat_expr(&ctx, &alloc, pattern_str);
        let bindings = match_pattern(pat, f).unwrap();
        let mid = bindings.get(Symbol::new("mid")).unwrap();
        assert!(matches!(mid, MatchValue::Sequence(s) if s.len() == 1));
    }
}
