//! Chain-entry bookkeeping for the integration pipeline (0.28.0).
//!
//! The pipeline re-enters [`crate::integral::integrate_raw`] from many
//! stages (Weierstrass/parts interactions, substitution back-substitution,
//! bounded expansion retries, rule residuals…). Two failure modes have to
//! be caught deterministically:
//!
//! - **True cycles**: the same expression re-enters the chain while it is
//!   already an ancestor on the call stack. 0.27.3 measured one corpus case
//!   (`rubi-00008`) re-entering the chain 306 times over 18 rounds before the
//!   *global* entry cap cut it off — a cap that also killed legitimate
//!   re-entrant work.
//! - **Runaway total work**: an absolute backstop, kept only as defence
//!   against loops that are mathematically the same shape but structurally
//!   different.
//!
//! The cycle test is done at the **expression level**: `Atom` is
//! hash-consed, so within one arena structural equality is pointer
//! equality, and an ancestor stack of node addresses decides identity in
//! O(1) per entry. Only ancestors count — repeated sibling sub-integrals
//! (e.g. two equal terms of an expanded sum) are legitimate and must pass.
//!
//! Because a `thread_local!` cannot hold the arena's lifetime, the stack
//! stores the addresses as `usize`. That is sound here: every node lives in
//! the arena for the whole top-level call, and the stack is only compared
//! while that call is active.

use std::cell::RefCell;

use ocas_atom::Atom;

/// Absolute backstop on chain entries per top-level `integrate` call.
///
/// **Kept at 0.27.3's 256.** 0.28.0 initially raised it to 1024 on the theory
/// that the expression-level cycle test would handle the ping-pongs; it
/// does, but grinders that are *not* ancestor cycles still burn every entry
/// they are given — measured: corpus wall clock 507 s → 806 s and timeouts
/// 12 → 18, with `rubi-00260` alone going from ~7 s to 30 s at 34 082 trace
/// lines. The backstop bounds work; it is not a correctness device.
pub(crate) const MAX_CHAIN_ENTRIES: u32 = 256;

/// Separate entry budget for nested residue resolution.
///
/// Resolution re-enters the chain, and charging those entries to the same
/// 256-entry pool starves the *primary* chain of the stage retries that used
/// to solve the case (measured: `rubi-01646` is solved with the pool to
/// itself and a fallback once resolution has spent part of it).
pub(crate) const MAX_RESIDUAL_ENTRIES: u32 = 128;

/// Chain entry identity: the expression's node address plus the recursion
/// budgets that were active at entry.
///
/// **Why the budgets are part of the key**: re-entering the *same* shape with
/// a smaller `parts_depth`/`rule_depth` is how the pipeline terminates its
/// legitimate recursions (integration by parts re-enters the chain on its own
/// sub-integrals with a reduced budget and eventually declines). Treating
/// shape repetition alone as a cycle kills that progress — measured:
/// `rubi-01646` is solved by exactly such a budget-decreasing re-entry and
/// regressed to a fallback with a shape-only key. A repeat with the *same*
/// budgets is a true cycle.
#[derive(Clone, Copy, PartialEq, Eq)]
struct ChainKey {
    node: usize,
    rule_depth: usize,
    parts_depth: usize,
}

thread_local! {
    static STATE: RefCell<ChainState> = const { RefCell::new(ChainState::new()) };
    /// Non-zero while a residue resolution is in progress.
    static RESOLVE_DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Non-zero while an integration-by-parts sub-integral is being
    /// computed (0.29.0 E2).
    static PARTS_DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Whether a stage returned a **partial** (residue-carrying) result
    /// during the current top-level call.
    ///
    /// The pipeline finishes with that partial unchanged — 0.27.3's flow —
    /// and the top-level entry point resolves its residues once, after the
    /// chain is done. Resolving *inside* the chain changes the stage budgets
    /// the outer recursion depends on and regressed `rubi-01646` from solved
    /// to fallback (measured), while resolving at the top keeps every solve.
    static PARTIAL_SEEN: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Record that a stage returned a partial result.
pub(crate) fn note_partial() {
    PARTIAL_SEEN.with(|p| p.set(true));
}

/// Whether a stage returned a partial result in this top-level call.
pub(crate) fn partial_seen() -> bool {
    PARTIAL_SEEN.with(|p| p.get())
}

struct ChainState {
    /// Keys of the `try_risch_or_fallback` entries currently on the call
    /// stack, outermost first.
    stack: Vec<ChainKey>,
    /// Total entries charged to the current top-level call.
    entries: u32,
    /// Entries charged to the current top-level call while resolving
    /// residues.
    resolve_entries: u32,
}

impl ChainState {
    const fn new() -> Self {
        Self {
            stack: Vec::new(),
            entries: 0,
            resolve_entries: 0,
        }
    }
}

/// Run `f` in a residue-resolution scope, so the nested chain entries are
/// charged to [`MAX_RESIDUAL_ENTRIES`] instead of [`MAX_CHAIN_ENTRIES`].
pub(crate) fn with_resolve_scope<R>(f: impl FnOnce() -> R) -> R {
    struct Guard;
    impl Drop for Guard {
        fn drop(&mut self) {
            RESOLVE_DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
        }
    }
    RESOLVE_DEPTH.with(|d| d.set(d.get() + 1));
    let _guard = Guard;
    f()
}

/// Run `f` in an integration-by-parts sub-integral scope (0.29.0 E2).
///
/// In-chain residue resolution is suppressed while a parts sub-integral is
/// being computed: resolution changes a partial answer's *shape* (a residue
/// becomes a complete answer with `atan`/`log` terms), and the parts
/// continuation `∫ v·du` then has to integrate those transcendental pieces,
/// which burns the chain budget (measured: `rubi-01646` regresses from
/// solved to fallback). The flag deliberately survives substitution
/// boundaries (which reset the numeric `parts_depth` budget).
pub(crate) fn with_parts_scope<R>(f: impl FnOnce() -> R) -> R {
    struct Guard;
    impl Drop for Guard {
        fn drop(&mut self) {
            PARTS_DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
        }
    }
    PARTS_DEPTH.with(|d| d.set(d.get() + 1));
    let _guard = Guard;
    f()
}

/// Whether a parts sub-integral is anywhere on the current call stack.
pub(crate) fn parts_active() -> bool {
    PARTS_DEPTH.with(|d| d.get() > 0)
}

/// Reset the chain bookkeeping. Called by every public integration entry
/// point (and by the residue resolver, which budgets separately).
pub(crate) fn reset() {
    STATE.with(|state| {
        let mut state = state.borrow_mut();
        state.stack.clear();
        state.entries = 0;
        state.resolve_entries = 0;
    });
    PARTIAL_SEEN.with(|p| p.set(false));
}

/// Result of trying to enter the chain for one expression.
pub(crate) enum Enter {
    /// Entry granted; the guard pops the stack on drop.
    Ok(ChainGuard),
    /// `expr` is already an ancestor: a true cycle.
    Cycle,
    /// The absolute backstop was exceeded.
    Exhausted,
}

/// RAII guard for one chain entry.
pub(crate) struct ChainGuard(());

impl Drop for ChainGuard {
    fn drop(&mut self) {
        STATE.with(|state| {
            state.borrow_mut().stack.pop();
        });
    }
}

/// Try to enter the chain for `expr` with the given recursion budgets.
pub(crate) fn enter(expr: Atom<'_>, rule_depth: usize, parts_depth: usize) -> Enter {
    let key = ChainKey {
        node: expr.node() as *const _ as usize,
        rule_depth,
        parts_depth,
    };
    let resolving = RESOLVE_DEPTH.with(|d| d.get() > 0);
    STATE.with(|state| {
        let mut state = state.borrow_mut();
        if state.stack.contains(&key) {
            return Enter::Cycle;
        }
        if resolving {
            if state.resolve_entries >= MAX_RESIDUAL_ENTRIES {
                return Enter::Exhausted;
            }
            state.resolve_entries = state.resolve_entries.saturating_add(1);
        } else {
            if state.entries >= MAX_CHAIN_ENTRIES {
                return Enter::Exhausted;
            }
            state.entries = state.entries.saturating_add(1);
        }
        state.stack.push(key);
        Enter::Ok(ChainGuard(()))
    })
}

/// Entries charged so far by the current top-level call (diagnostics only).
pub(crate) fn entries() -> u32 {
    STATE.with(|state| state.borrow().entries)
}

/// Diagnostic helper: the depth of the current chain stack.
#[cfg(test)]
pub(crate) fn depth() -> usize {
    STATE.with(|state| state.borrow().stack.len())
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn ancestor_reentry_is_a_cycle() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        reset();
        let _outer = match enter(x, 0, 0) {
            Enter::Ok(guard) => guard,
            _ => panic!("first entry must be granted"),
        };
        assert!(matches!(enter(x, 0, 0), Enter::Cycle));
        drop(_outer);
        // Popped: the same expression can be entered again.
        assert!(matches!(enter(x, 0, 0), Enter::Ok(_)));
        reset();
    }

    #[test]
    fn budget_decreasing_reentry_is_allowed() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        reset();
        // The pipeline terminates its legitimate recursions by shrinking the
        // parts/rule budgets; re-entering the same shape with a smaller
        // budget is progress, not a cycle.
        let _outer = match enter(x, 0, 4) {
            Enter::Ok(guard) => guard,
            _ => panic!("first entry must be granted"),
        };
        assert!(matches!(enter(x, 0, 3), Enter::Ok(_)));
        assert!(matches!(enter(x, 1, 0), Enter::Ok(_)));
        reset();
    }

    #[test]
    fn sibling_repeats_are_allowed() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sum = ctx.add(&[x, x]);
        reset();
        // Two equal children are entered one after the other (not nested):
        // each must be granted.
        assert!(matches!(enter(x, 0, 0), Enter::Ok(_)));
        assert!(matches!(enter(sum, 0, 0), Enter::Ok(_)));
        reset();
    }

    #[test]
    fn backstop_is_per_top_level_call() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        reset();
        // Enter/leave the same expression more times than the backstop.
        let mut granted = 0;
        for _ in 0..(MAX_CHAIN_ENTRIES + 8) {
            if let Enter::Ok(guard) = enter(ctx.var("x"), 0, 0) {
                granted += 1;
                drop(guard);
            }
        }
        assert_eq!(granted, MAX_CHAIN_ENTRIES as usize);
        // The backstop is checked before the counter advances, so the
        // counter stops exactly at the cap.
        assert_eq!(entries(), MAX_CHAIN_ENTRIES);
        reset();
        assert_eq!(entries(), 0);
        assert_eq!(depth(), 0);
    }
}
