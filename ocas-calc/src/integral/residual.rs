//! Budgeted residue resolution for the integration chain (0.28.0).
//!
//! Several stages return a **partial** result: a valid decomposition whose
//! remaining `Integral(g, v)` term is itself an elementary integral the
//! pipeline already knows how to do. `1/(x²(x²+5))` is the canonical
//! example — the partial fraction is computed, and the leftover
//! `Integral(-1/5·(x²+5)⁻¹, x)` is left in the answer even though
//! `1/(x²+5)` integrates fine.
//!
//! 0.27.3 measured a prototype that resolved such residues *inside*
//! `rational`/`symbolic_rational`: net +4 solved, but one regression and
//! three new timeouts. This module therefore exposes a **budgeted** resolver
//! that the *top-level entry points* drive, after the chain is done:
//!
//! - [`super::chain::note_partial`] records that the pipeline finished with
//!   a partial result (0.27.3's flow is preserved untouched);
//! - the top-level entry resolves that result's residues once.
//!
//! Resolving *inside* the chain was measured to change the stage budgets the
//! outer recursion depends on (`rubi-01646` regressed from solved to
//! fallback). Resolving at the top preserves every solve.
//!
//! Four invariants keep resolution from doing harm:
//!
//! 1. **Strict progress** — a resolution is accepted only if the number of
//!    `Integral(·, v)` nodes strictly decreases; otherwise the original
//!    expression is returned unchanged.
//! 2. **Shape gate** — only residues that are *rational functions* over at
//!    most [`MAX_RESIDUAL_VARS`] symbols are re-integrated. That is exactly
//!    the measured defect; residues carrying function heads cost seconds
//!    each and are left alone.
//! 3. **Size gate** — residues above [`MAX_RESIDUAL_NODES`] are left alone.
//! 4. **Step and depth budgets** — a per-top-level-call step budget plus a
//!    nesting depth cap, so nested resolutions cannot recurse without bound.
//!
//! The plain fallback `Integral(f, x)` (the pipeline's "nothing worked"
//! answer) is never resolved: it would re-run the whole chain — and, with
//! `rules = false`, would let the rule table solve a case the caller
//! explicitly disabled it for.

use std::cell::Cell;

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::node_count;

/// Resolution steps charged by one top-level call.
pub(crate) const MAX_RESIDUAL_STEPS: u32 = 8;

/// Nesting cap for resolutions (a resolution may itself produce residues).
pub(crate) const MAX_RESIDUAL_DEPTH: u32 = 2;

/// Node ceiling for a residue this resolver will re-integrate.
pub(crate) const MAX_RESIDUAL_NODES: usize = 64;

thread_local! {
    /// `(steps_left, depth)` for the current top-level call.
    static STATE: Cell<(u32, u32)> = const { Cell::new((MAX_RESIDUAL_STEPS, 0)) };
}

/// Reset the residue budget; called from the public entry points together
/// with [`super::chain::reset`].
pub(crate) fn reset() {
    STATE.with(|state| state.set((MAX_RESIDUAL_STEPS, 0)));
}

/// Whether `expr` contains an `Integral(g, v)` residue anywhere.
///
/// **Any** integration variable counts, not just the top-level one: the
/// mechanism substitutions (`t = √x`, `t = e^u`, Weierstrass `_t`, …) leave
/// residues in their own variable, and the corpus harness counts any
/// `Integral(` in the printed result as unresolved.
pub(crate) fn contains_residue(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "Integral" && args.len() == 2 => {
            if matches!(args[1].node(), AtomNode::Var(_)) {
                return true;
            }
            args.iter().any(|a| contains_residue(*a))
        }
        AtomNode::Fun(_, args) | AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_residue(*a))
        }
        AtomNode::Pow(base, exp) => contains_residue(*base) || contains_residue(*exp),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Whether `expr` is a rational function of its variables — i.e. contains no
/// function applications and no non-integer powers.
///
/// Residues of this shape are the 0.28.0 defect: the partial fraction is
/// already computed and the leftover piece is a plain rational integral.
/// Residues carrying function heads (`cosh`, `log`, …) are the shapes that
/// made resolution expensive (measured: `1/(a + b·cosh + c·cosh²)` went from
/// 5.6 s to 48 s when every residue was re-integrated), so they are left
/// alone.
pub(crate) fn is_resolvable_shape(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Fun(_, _) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|a| is_resolvable_shape(*a)),
        AtomNode::Pow(base, exp) => {
            matches!(exp.node(), AtomNode::Num(_))
                && is_resolvable_shape(*base)
                && is_resolvable_shape(*exp)
        }
    }
}

/// Maximum number of distinct free symbols a residue may mention.
///
/// The defect's residues are univariate (or carry one symbolic coefficient).
/// Measured: Weierstrass/hyperbolic t-form residues with four symbols
/// (`t, a, b, c`) are rational too, but re-integrating them costs ~0.5 s
/// each — 44 of them took 24 s on `1/(a + b·cosh + c·cosh²)` — so they are
/// left alone.
const MAX_RESIDUAL_VARS: usize = 2;

/// Whether `expr` mentions at most [`MAX_RESIDUAL_VARS`] distinct symbols.
fn few_symbols(expr: Atom<'_>) -> bool {
    let mut seen: Vec<Symbol> = Vec::with_capacity(MAX_RESIDUAL_VARS);
    fn walk(expr: Atom<'_>, seen: &mut Vec<Symbol>) -> bool {
        match expr.node() {
            AtomNode::Num(_) => true,
            AtomNode::Var(v) => {
                if seen.contains(v) {
                    return true;
                }
                if seen.len() >= MAX_RESIDUAL_VARS {
                    return false;
                }
                seen.push(*v);
                true
            }
            AtomNode::Fun(_, args) | AtomNode::Add(args) | AtomNode::Mul(args) => {
                args.iter().all(|a| walk(*a, seen))
            }
            AtomNode::Pow(base, exp) => walk(*base, seen) && walk(*exp, seen),
        }
    }
    walk(expr, &mut seen)
}

/// Number of `Integral(g, v)` residues in `expr`.
pub(crate) fn count_residues(expr: Atom<'_>) -> usize {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "Integral" && args.len() == 2 => {
            let own = usize::from(matches!(args[1].node(), AtomNode::Var(_)));
            own + args.iter().map(|a| count_residues(*a)).sum::<usize>()
        }
        AtomNode::Fun(_, args) | AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().map(|a| count_residues(*a)).sum()
        }
        AtomNode::Pow(base, exp) => count_residues(*base) + count_residues(*exp),
        AtomNode::Num(_) | AtomNode::Var(_) => 0,
    }
}

/// Re-integrate the `Integral(g, v)` residues of `expr`, each with respect
/// to its own variable `v`.
///
/// `max_steps` caps how many residues this call may charge, and
/// `parts_depth` is the caller's integration-by-parts recursion budget,
/// forwarded unchanged: a resolution that reset it to zero would let the
/// parts↔Weierstrass ping-pong start over from the top on every residue
/// (measured: `rubi-01646` burned the whole chain backstop before the trace
/// showed the loop).
///
/// Returns `expr` unchanged when there is no budget left, when the depth cap
/// is reached, or when resolution did not strictly reduce the residue count.
pub(crate) fn resolve<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
    max_steps: u32,
    node_ceiling: usize,
) -> Atom<'a> {
    if !contains_residue(expr) {
        return expr;
    }
    let (steps_left, depth) = STATE.with(|state| state.get());
    if steps_left == 0 || depth >= MAX_RESIDUAL_DEPTH {
        return expr;
    }
    let allowance = steps_left.min(max_steps);
    STATE.with(|state| state.set((steps_left, depth + 1)));
    if super::trace_enabled() {
        eprintln!(
            "[trace] resolve-entry residues={} steps={} depth={}",
            count_residues(expr),
            steps_left,
            depth
        );
    }

    let before = count_residues(expr);
    let mut walker = ResolveWalk {
        rules_enabled,
        rule_depth,
        parts_depth,
        node_ceiling,
        allowance,
        used: 0,
    };
    // Charge the nested chain entries to the resolution budget, so a
    // resolution attempt cannot starve the primary chain's stage retries.
    let out = super::chain::with_resolve_scope(|| walker.walk(ctx, expr));
    let used = walker.used;

    STATE.with(|state| state.set((steps_left.saturating_sub(used), depth)));

    if count_residues(out) < before {
        out
    } else {
        expr
    }
}

/// Shared state for one resolution pass.
struct ResolveWalk {
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
    node_ceiling: usize,
    allowance: u32,
    used: u32,
}

impl ResolveWalk {
    fn walk<'a>(&mut self, ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
        match expr.node() {
            AtomNode::Fun(name, args) if name.as_str() == "Integral" && args.len() == 2 => {
                let integrand = args[0];
                let AtomNode::Var(var) = args[1].node() else {
                    return expr;
                };
                if self.used >= self.allowance || node_count(integrand) > self.node_ceiling {
                    if super::trace_enabled() {
                        eprintln!(
                            "[trace] resolve-skip size nodes={} used={} allowance={}",
                            node_count(integrand),
                            self.used,
                            self.allowance
                        );
                    }
                    // Out of budget, or a residue big enough to be one of the
                    // measured timeout shapes: leave it in place.
                    return expr;
                }
                if !is_resolvable_shape(integrand) || !few_symbols(integrand) {
                    if super::trace_enabled() {
                        eprintln!(
                            "[trace] resolve-skip shape rational={} symbols={} :: {integrand}",
                            is_resolvable_shape(integrand),
                            few_symbols(integrand)
                        );
                    }
                    // A non-rational residue (function heads, radicals) is
                    // not the defect this resolver exists for, and resolving
                    // it is what made the measured cases expensive.
                    return expr;
                }
                if super::trace_enabled() {
                    eprintln!(
                        "[trace] resolve-attempt {} :: {integrand}",
                        node_count(integrand)
                    );
                }
                self.used += 1;
                let resolved = crate::integral::integrate_raw(
                    ctx,
                    integrand,
                    *var,
                    0,
                    self.rules_enabled,
                    self.rule_depth,
                    self.parts_depth,
                );
                if resolved == expr {
                    // No progress: keep the original residue so the strict
                    // decrease gate sees it.
                    return expr;
                }
                resolved
            }
            AtomNode::Fun(name, args) => {
                let rebuilt: Vec<Atom<'a>> = args.iter().map(|a| self.walk(ctx, *a)).collect();
                let rebuilt = ctx.fun(name.as_str(), &rebuilt);
                if rebuilt == expr { expr } else { rebuilt }
            }
            AtomNode::Add(args) => {
                let rebuilt: Vec<Atom<'a>> = args.iter().map(|a| self.walk(ctx, *a)).collect();
                let rebuilt = ctx.add(&rebuilt);
                if rebuilt == expr { expr } else { rebuilt }
            }
            AtomNode::Mul(args) => {
                let rebuilt: Vec<Atom<'a>> = args.iter().map(|a| self.walk(ctx, *a)).collect();
                let rebuilt = ctx.mul(&rebuilt);
                if rebuilt == expr { expr } else { rebuilt }
            }
            AtomNode::Pow(base, exp) => {
                let b = self.walk(ctx, *base);
                let e = self.walk(ctx, *exp);
                let rebuilt = ctx.pow(b, e);
                if rebuilt == expr { expr } else { rebuilt }
            }
            AtomNode::Num(_) | AtomNode::Var(_) => expr,
        }
    }
}

#[cfg(test)]
mod tests {
    use ocas_core::arena::Arena;
    use ocas_parse::parse;

    use super::*;

    fn prep<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_atom::normalize::normalize(ctx, parse(ctx, s).expect("parse"))
    }

    #[test]
    fn counts_residues_in_any_variable() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x_residue = prep(&ctx, "Integral(x^2, x)");
        let y_residue = prep(&ctx, "Integral(x^2, y)");
        assert_eq!(count_residues(x_residue), 1);
        assert_eq!(count_residues(y_residue), 1);
        assert!(contains_residue(x_residue));
        assert!(contains_residue(y_residue));
        // Nesting inside a sum and a power is counted too.
        let nested = prep(&ctx, "Integral(x^2, x) + (Integral(x, y) + 1)^2");
        assert_eq!(count_residues(nested), 2);
        assert!(!contains_residue(prep(&ctx, "x + 1")));
    }

    #[test]
    fn resolves_a_solvable_residue() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        let expr = prep(&ctx, "Integral(x^2, x)");
        let out = resolve(
            &ctx,
            expr,
            true,
            0,
            0,
            MAX_RESIDUAL_STEPS,
            MAX_RESIDUAL_NODES,
        );
        assert!(!contains_residue(out));
    }

    #[test]
    fn resolves_a_substitution_variable_residue() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        // The mechanism substitutions leave residues in their own variable.
        let expr = prep(&ctx, "atan(t) + Integral(1/(1 + t^2), t)");
        let out = resolve(
            &ctx,
            expr,
            true,
            0,
            0,
            MAX_RESIDUAL_STEPS,
            MAX_RESIDUAL_NODES,
        );
        assert!(!contains_residue(out), "unresolved: {out}");
    }

    #[test]
    fn keeps_an_unsolvable_residue_unchanged() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        // An unknown function head has no mechanism and no special-function
        // table entry, so resolution cannot make progress.
        let expr = prep(&ctx, "Integral(foo(x), x)");
        let out = resolve(
            &ctx,
            expr,
            true,
            0,
            0,
            MAX_RESIDUAL_STEPS,
            MAX_RESIDUAL_NODES,
        );
        assert_eq!(out, expr);
    }

    #[test]
    fn non_rational_residues_are_left_alone() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        let expr = prep(&ctx, "Integral(log(x)/(1 + x^2), x)");
        let out = resolve(
            &ctx,
            expr,
            true,
            0,
            0,
            MAX_RESIDUAL_STEPS,
            MAX_RESIDUAL_NODES,
        );
        assert_eq!(out, expr);
    }

    #[test]
    fn budget_exhaustion_leaves_the_residue() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        STATE.with(|state| state.set((0, 0)));
        let expr = prep(&ctx, "Integral(x^2, x)");
        let out = resolve(
            &ctx,
            expr,
            true,
            0,
            0,
            MAX_RESIDUAL_STEPS,
            MAX_RESIDUAL_NODES,
        );
        assert_eq!(out, expr);
        reset();
    }

    #[test]
    fn oversized_residues_are_left_alone() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        let expr = prep(&ctx, "Integral(x^2, x)");
        let out = resolve(&ctx, expr, true, 0, 0, MAX_RESIDUAL_STEPS, 1);
        assert_eq!(out, expr);
        reset();
    }
}
