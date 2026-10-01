//! Merging of algebraically dependent tower generators (0.28.0).
//!
//! `tower/build.rs` used to *reject* an integrand whose generators were
//! algebraically dependent (`log(x)` with `log(2x)`, `exp(x)` with
//! `exp(x+1)`). That rejection had a large measured cost: rewriting
//! `sinh(u) = (eᵘ − e⁻ᵘ)/2` — which every hyperbolic integrand does on its
//! way into Risch — necessarily produces the dependent pair `exp(u)` /
//! `exp(−u)`, so the whole hyperbolic family was locked out of the general
//! engine (0.28.0 baseline: 246 unsolved corpus integrands mention a
//! hyperbolic head).
//!
//! This module expresses a dependent candidate generator in terms of the
//! generators already in the tower, using exact field identities:
//!
//! | candidate | relation | rewrite |
//! |---|---|---|
//! | `exp(u)` vs `exp(v)` | `u + v` constant | `exp(u) = exp(u+v)/exp(v)` |
//! | `exp(u)` vs `exp(v)` | `u − v` constant | `exp(u) = exp(v)·exp(u−v)` |
//! | `log(u)` vs `log(v)` | `u = v^k` | `log(u) = k·log(v)` |
//! | `log(u)` vs `log(v)` | `v = u^k` | `log(u) = log(v)/k` |
//! | `log(u)` vs `log(v)` | `u/v` constant | `log(u) = log(v) + log(u/v)` |
//! | `log(u)` vs `log(v)` | `u·v` constant | `log(u) = log(u·v) − log(v)` |
//! | `log(exp(u))` | — | `u` |
//! | `exp(log(u))` | — | `u` |
//!
//! Each rewrite may need a *constant* atom such as `log(2)` or `exp(1)`;
//! those are registered as `GenKind::Constant` generators (`D t = 0`),
//! which is exactly what the constant generator slot already exists for.
//!
//! The rewrites are exact in the differential field up to an additive
//! constant of integration, which is harmless: a constant generator has
//! zero derivative, so every `dt` relation the tower stores stays true.
//!
//! A candidate that matches no rule is treated as independent, exactly as
//! before — the old structural detector was already incomplete, so no
//! dependent pair that used to be accepted can slip through, and pairs that
//! used to be rejected now have a chance to merge.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::build::{GenInfo, GenKind};
use crate::integral::is_constant;
use crate::ode::util::collect_terms;

/// The rewrite of a dependent candidate generator.
pub(crate) struct Merge<'a> {
    /// Expression the candidate's application atom is replaced by.
    pub replacement: Atom<'a>,
    /// Constant atom that must be registered as a `GenKind::Constant`
    /// generator (e.g. `log(2)`), when the rewrite introduced one.
    pub constant: Option<Atom<'a>>,
}

/// Try to express `app` (with kind `kind` and argument `arg`) in terms of
/// the generators already collected in `gens`.
///
/// Returns `None` when no relation is recognised — the caller then treats
/// the candidate as a new, independent generator.
pub(crate) fn merge_candidate<'a>(
    ctx: &'a AtomArena<'a>,
    kind: GenKind,
    arg: Atom<'a>,
    app: Atom<'a>,
    gens: &[GenInfo<'a>],
    var: Symbol,
) -> Option<Merge<'a>> {
    for g in gens {
        match (kind, g.kind) {
            (GenKind::Exp, GenKind::Exp) => {
                // exp(u)·exp(v) = exp(u+v).
                let sum = collect_terms(ctx, ctx.add(&[arg, g.arg]));
                if is_constant(sum, var) {
                    if is_zero(sum) {
                        // exp(u) = 1/exp(v)
                        return Some(Merge {
                            replacement: ctx.pow(g.atom, ctx.num(-1)),
                            constant: None,
                        });
                    }
                    let constant = ctx.fun("exp", &[sum]);
                    return Some(Merge {
                        replacement: ctx.mul(&[constant, ctx.pow(g.atom, ctx.num(-1))]),
                        constant: Some(constant),
                    });
                }
                // exp(u) = exp(v)·exp(u−v).
                let minus_one = ctx.num(-1);
                let diff = collect_terms(ctx, ctx.add(&[arg, ctx.mul(&[minus_one, g.arg])]));
                if is_constant(diff, var) {
                    if is_zero(diff) {
                        return Some(Merge {
                            replacement: g.atom,
                            constant: None,
                        });
                    }
                    let constant = ctx.fun("exp", &[diff]);
                    return Some(Merge {
                        replacement: ctx.mul(&[g.atom, constant]),
                        constant: Some(constant),
                    });
                }
            }
            (GenKind::Log, GenKind::Log) => {
                if let Some(k) = integer_power_of(arg, g.arg) {
                    // log(v^k) = k·log(v)
                    return Some(Merge {
                        replacement: ctx.mul(&[ctx.num(k), g.atom]),
                        constant: None,
                    });
                }
                if let Some(k) = integer_power_of(g.arg, arg) {
                    // log(u) = log(v)/k when v = u^k.
                    return Some(Merge {
                        replacement: ctx.mul(&[g.atom, ctx.pow(ctx.num(k), ctx.num(-1))]),
                        constant: None,
                    });
                }
                let ratio = collect_terms(ctx, ctx.mul(&[arg, ctx.pow(g.arg, ctx.num(-1))]));
                if is_constant(ratio, var) {
                    if is_one(ratio) {
                        return Some(Merge {
                            replacement: g.atom,
                            constant: None,
                        });
                    }
                    let constant = ctx.fun("log", &[ratio]);
                    return Some(Merge {
                        replacement: ctx.add(&[g.atom, constant]),
                        constant: Some(constant),
                    });
                }
                let product = collect_terms(ctx, ctx.mul(&[arg, g.arg]));
                if is_constant(product, var) {
                    let minus_g = ctx.mul(&[ctx.num(-1), g.atom]);
                    if is_one(product) {
                        return Some(Merge {
                            replacement: minus_g,
                            constant: None,
                        });
                    }
                    let constant = ctx.fun("log", &[product]);
                    return Some(Merge {
                        replacement: ctx.add(&[constant, minus_g]),
                        constant: Some(constant),
                    });
                }
            }
            // log(exp(u)) = u and exp(log(u)) = u.
            (GenKind::Log, GenKind::Exp) | (GenKind::Exp, GenKind::Log) if arg == g.atom => {
                return Some(Merge {
                    replacement: g.arg,
                    constant: None,
                });
            }
            _ => {}
        }
        let _ = app;
    }
    None
}

/// `Some(k)` when `base^k == powered` structurally for an integer `k`.
fn integer_power_of(powered: Atom<'_>, base: Atom<'_>) -> Option<i64> {
    match powered.node() {
        AtomNode::Pow(b, e) if *b == base => match e.node() {
            AtomNode::Num(k) if *k != 0 => Some(*k),
            _ => None,
        },
        _ => None,
    }
}

fn is_zero(atom: Atom<'_>) -> bool {
    matches!(atom.node(), AtomNode::Num(0))
}

fn is_one(atom: Atom<'_>) -> bool {
    matches!(atom.node(), AtomNode::Num(1))
}

/// Replace every occurrence of `old` in `expr` by `new`.
pub(crate) fn substitute_atom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    old: Atom<'a>,
    new: Atom<'a>,
) -> Atom<'a> {
    if expr == old {
        return new;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, old, new))
                .collect();
            let rebuilt = ctx.fun(name.as_str(), &rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, old, new))
                .collect();
            let rebuilt = ctx.add(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| substitute_atom(ctx, *a, old, new))
                .collect();
            let rebuilt = ctx.mul(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Pow(base, exp) => {
            let b = substitute_atom(ctx, *base, old, new);
            let e = substitute_atom(ctx, *exp, old, new);
            let rebuilt = ctx.pow(b, e);
            if rebuilt == expr { expr } else { rebuilt }
        }
    }
}

#[cfg(test)]
mod tests {
    use ocas_core::arena::Arena;
    use ocas_parse::parse;

    use super::*;

    fn rebuild<'a>(ctx: &'a AtomArena<'a>, src: &str) -> Atom<'a> {
        ocas_atom::normalize::normalize(ctx, parse(ctx, src).expect("parse"))
    }

    #[test]
    fn substitute_replaces_nested_occurrences() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let old = ctx.fun("log", &[ctx.var("x")]);
        let new = ctx.var("u");
        let expr = rebuild(&ctx, "log(x)^2 + exp(log(x)) + 1");
        let out = substitute_atom(&ctx, expr, old, new);
        assert_eq!(out.to_string(), "1 + (exp(u)) + (u^2)");
    }

    #[test]
    fn exp_of_negation_merges_to_an_inverse() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let exp_x = ctx.fun("exp", &[x]);
        let gens = vec![GenInfo {
            kind: GenKind::Exp,
            atom: exp_x,
            arg: x,
            dt: crate::tower::elem::KElem::zero(2),
        }];
        let neg_x = rebuild(&ctx, "-x");
        let merge = merge_candidate(
            &ctx,
            GenKind::Exp,
            neg_x,
            ctx.fun("exp", &[neg_x]),
            &gens,
            Symbol::new("x"),
        )
        .expect("merge");
        assert_eq!(merge.replacement.to_string(), "(exp(x))^-1");
        assert!(merge.constant.is_none());
    }

    #[test]
    fn exp_shift_needs_a_constant_generator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let exp_x = ctx.fun("exp", &[x]);
        let gens = vec![GenInfo {
            kind: GenKind::Exp,
            atom: exp_x,
            arg: x,
            dt: crate::tower::elem::KElem::zero(2),
        }];
        let shifted = rebuild(&ctx, "x + 1");
        let merge = merge_candidate(
            &ctx,
            GenKind::Exp,
            shifted,
            ctx.fun("exp", &[shifted]),
            &gens,
            Symbol::new("x"),
        )
        .expect("merge");
        assert!(merge.constant.is_some());
        assert_eq!(merge.constant.unwrap().to_string(), "exp(1)");
    }

    #[test]
    fn dependent_logs_merge_via_the_ratio() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let log_x = ctx.fun("log", &[x]);
        let gens = vec![GenInfo {
            kind: GenKind::Log,
            atom: log_x,
            arg: x,
            dt: crate::tower::elem::KElem::zero(2),
        }];
        let two_x = rebuild(&ctx, "2*x");
        let merge = merge_candidate(
            &ctx,
            GenKind::Log,
            two_x,
            ctx.fun("log", &[two_x]),
            &gens,
            Symbol::new("x"),
        )
        .expect("merge");
        assert_eq!(merge.constant.unwrap().to_string(), "log(2)");
    }

    #[test]
    fn independent_candidates_do_not_merge() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let exp_x = ctx.fun("exp", &[x]);
        let gens = vec![GenInfo {
            kind: GenKind::Exp,
            atom: exp_x,
            arg: x,
            dt: crate::tower::elem::KElem::zero(2),
        }];
        let x_squared = rebuild(&ctx, "x^2");
        assert!(
            merge_candidate(
                &ctx,
                GenKind::Exp,
                x_squared,
                ctx.fun("exp", &[x_squared]),
                &gens,
                Symbol::new("x"),
            )
            .is_none()
        );
        let log_x = ctx.fun("log", &[x]);
        let log_gens = vec![GenInfo {
            kind: GenKind::Log,
            atom: log_x,
            arg: x,
            dt: crate::tower::elem::KElem::zero(2),
        }];
        let x_plus_one = rebuild(&ctx, "x + 1");
        assert!(
            merge_candidate(
                &ctx,
                GenKind::Log,
                x_plus_one,
                ctx.fun("log", &[x_plus_one]),
                &log_gens,
                Symbol::new("x"),
            )
            .is_none()
        );
    }
}
