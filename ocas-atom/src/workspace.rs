//! Thread-local arena pool for high-churn computations.
//!
//! Normalization, rewriting, and other temporary-heavy computations
//! repeatedly build and discard expression arenas. [`WorkspaceArena`]
//! recycles [`Arena`] allocations through a thread-local pool (the
//! RecycledAtom/Workspace pattern from Symbolica's `state.rs`), so
//! steady-state workloads allocate no new memory between generations.
//!
//! # Usage
//!
//! ```
//! use ocas_atom::workspace::WorkspaceArena;
//! use ocas_atom::AtomArena;
//!
//! let ws = WorkspaceArena::acquire();
//! let ctx = AtomArena::new(ws.arena());
//! let x = ctx.var("x");
//! let _ = ctx.add(&[x, ctx.num(1)]);
//! // `ctx` must be dropped before `ws`; on drop, `ws` returns the arena
//! // to the thread-local pool for the next computation.
//! ```

use std::cell::RefCell;

use ocas_core::arena::Arena;

/// Maximum number of arenas retained in each thread's pool. Excess
/// arenas are dropped to bound peak memory.
const MAX_POOLED: usize = 4;

thread_local! {
    #[allow(clippy::missing_const_for_thread_local)]
    static ARENA_POOL: RefCell<Vec<Arena>> = RefCell::new(Vec::new());
}

/// A handle to an [`Arena`] borrowed from the thread-local pool.
///
/// On drop, the arena is reset (all allocations invalidated) and
/// returned to the pool. Because the pool is thread-local, this type
/// is `!Send`.
pub struct WorkspaceArena {
    arena: Option<Arena>,
    // Make the handle !Send/!Sync: the pool is thread-local, so moving
    // the handle across threads would return the arena to the wrong pool.
    _not_send: std::marker::PhantomData<*const ()>,
}

impl WorkspaceArena {
    /// Acquire an arena from the thread-local pool, or create a fresh
    /// one if the pool is empty.
    pub fn acquire() -> Self {
        let arena = ARENA_POOL
            .with(|pool| pool.borrow_mut().pop())
            .unwrap_or_default();
        Self {
            arena: Some(arena),
            _not_send: std::marker::PhantomData,
        }
    }

    /// Access the pooled arena.
    pub fn arena(&self) -> &Arena {
        self.arena.as_ref().expect("arena present until drop")
    }

    /// Reset the arena in place, invalidating all values allocated so
    /// far. References handed out before the reset must not be used
    /// afterwards.
    pub fn reset(&mut self) {
        self.arena
            .as_ref()
            .expect("arena present until drop")
            .reset();
    }
}

impl Drop for WorkspaceArena {
    fn drop(&mut self) {
        if let Some(arena) = self.arena.take() {
            // Return the arena to the pool in a clean state for the next
            // computation.
            arena.reset();
            ARENA_POOL.with(|pool| {
                let mut pool = pool.borrow_mut();
                if pool.len() < MAX_POOLED {
                    pool.push(arena);
                }
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;

    #[test]
    fn acquire_and_use() {
        let ws = WorkspaceArena::acquire();
        let ctx = AtomArena::new(ws.arena());
        let x = ctx.var("x");
        let expr = ctx.add(&[x, ctx.num(1)]);
        assert_eq!(expr.to_string(), "x + 1");
    }

    #[test]
    fn pool_reuses_arena_memory() {
        // First generation: allocate and drop, returning the arena to the pool.
        {
            let ws = WorkspaceArena::acquire();
            let ctx = AtomArena::new(ws.arena());
            let x = ctx.var("x");
            let _ = ctx.add(&[x, ctx.num(1)]);
            drop(ctx);
            // ws drops here, returning the arena to the pool.
        }
        // Second generation: an arena comes back from the pool (not a fresh
        // allocation). We verify reuse semantically: the pool is non-empty
        // after the first drop, and the second acquire pops from it.
        ARENA_POOL.with(|pool| assert!(!pool.borrow().is_empty()));
        let ws2 = WorkspaceArena::acquire();
        let ctx2 = AtomArena::new(ws2.arena());
        let y = ctx2.var("y");
        let expr = ctx2.add(&[y, ctx2.num(2)]);
        assert_eq!(expr.to_string(), "y + 2");
    }

    #[test]
    fn reset_generations_are_independent() {
        let mut ws = WorkspaceArena::acquire();
        for round in 0..100 {
            {
                let ctx = AtomArena::new(ws.arena());
                let x = ctx.var("x");
                let expr = ctx.add(&[x, ctx.num(round)]);
                assert_eq!(expr.to_string(), format!("x + {round}"));
            }
            ws.reset();
        }
    }

    #[test]
    fn pool_is_bounded() {
        // Dropping more handles than MAX_POOLED must not grow the pool.
        let handles: Vec<WorkspaceArena> = (0..16).map(|_| WorkspaceArena::acquire()).collect();
        drop(handles);
        ARENA_POOL.with(|pool| assert!(pool.borrow().len() <= MAX_POOLED));
    }

    #[test]
    fn stress_many_generations() {
        // 10k acquire/build/drop cycles; pool keeps memory bounded.
        for i in 0..10_000u64 {
            let ws = WorkspaceArena::acquire();
            {
                let ctx = AtomArena::new(ws.arena());
                let x = ctx.var("x");
                let mut acc = ctx.num(0);
                for j in 0..(i % 50) {
                    acc = ctx.add(&[acc, ctx.mul(&[x, ctx.num(j as i64)])]);
                }
                let _ = acc.to_string();
            }
        }
        ARENA_POOL.with(|pool| assert!(pool.borrow().len() <= MAX_POOLED));
    }
}
