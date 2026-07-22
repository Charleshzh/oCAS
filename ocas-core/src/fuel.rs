//! Resource budget for long-running algorithms.
//!
//! [`Fuel`] is a shared, decrementing budget that algorithms consume as they
//! work. When the budget is exhausted, callers observe an
//! [`OcasError::OutOfFuel`](crate::error::OcasError::OutOfFuel) and can stop
//! the computation deterministically.
//!
//! # Design
//!
//! `Fuel` wraps an [`AtomicUsize`] counter so that a single budget can be
//! shared (cheaply cloned as a handle) across nested calls — for example a
//! top-level `simplify` that drives many pattern matches. Consumption is
//! infallible: once the counter hits zero every subsequent [`Fuel::consume`]
//! or [`Fuel::check`] reports exhaustion. The budget is monotone
//! non-increasing, never refilled except by constructing a new `Fuel`.
//!
//! Algorithms that already keep per-call limits (e.g. `simplify`'s
//! `iter_limit`, `integrate`'s `MAX_DEPTH`) can opt into `Fuel` without
//! breaking their existing API by adding a `*_with_fuel` entry point.
//!
//! # Example
//!
//! ```
//! use ocas_core::fuel::Fuel;
//!
//! let fuel = Fuel::new(100);
//! for _ in 0..50 {
//!     fuel.consume(1);
//! }
//! assert!(fuel.check().is_ok());
//! fuel.consume(100);
//! assert!(fuel.check().is_err());
//! ```

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::error::{OcasError, Result};

/// A decrementing resource budget shared across algorithm invocations.
///
/// Clone produces a new handle to the **same** underlying counter; consumption
/// on any clone is observed by all clones. Use this to thread a single budget
/// through nested calls without passing `&Fuel` everywhere.
#[derive(Debug, Clone)]
pub struct Fuel {
    remaining: Arc<AtomicUsize>,
}

impl Fuel {
    /// Create a new budget with `budget` units.
    pub fn new(budget: usize) -> Self {
        Self {
            remaining: Arc::new(AtomicUsize::new(budget)),
        }
    }

    /// Consume `n` units. If the remaining budget would drop to zero or below,
    /// the counter saturates at zero (subsequent [`Self::check`] fails). This
    /// call is infallible — query the outcome with [`Self::check`] or read the
    /// residual via [`Self::remaining`].
    pub fn consume(&self, n: usize) {
        // Saturating subtraction on the atomic. Relaxed ordering is fine: fuel
        // is a best-effort cut-off, not a synchronisation primitive. We only
        // need eventual visibility across threads, which Relaxed provides.
        loop {
            let cur = self.remaining.load(Ordering::Relaxed);
            let next = cur.saturating_sub(n);
            // Try to commit; if another thread beat us, retry with their value.
            if self
                .remaining
                .compare_exchange_weak(cur, next, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                break;
            }
        }
    }

    /// Return `Ok(())` if budget remains, `Err(OutOfFuel)` otherwise. Cheap
    /// hot-path probe; use between small units of work to avoid an explicit
    /// [`Self::consume`] on every iteration.
    pub fn check(&self) -> Result<()> {
        if self.remaining.load(Ordering::Relaxed) == 0 {
            Err(OcasError::OutOfFuel)
        } else {
            Ok(())
        }
    }

    /// The number of units still available. `0` means exhausted.
    pub fn remaining(&self) -> usize {
        self.remaining.load(Ordering::Relaxed)
    }
}

impl Default for Fuel {
    fn default() -> Self {
        // A generous default; algorithms that opt in may override.
        Self::new(1_000_000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consume_then_check() {
        let f = Fuel::new(3);
        assert!(f.check().is_ok());
        f.consume(2);
        assert!(f.check().is_ok());
        f.consume(1);
        assert!(f.check().is_err());
        // Stays exhausted.
        f.consume(10);
        assert!(f.check().is_err());
    }

    #[test]
    fn remaining_is_monotone() {
        let f = Fuel::new(10);
        assert_eq!(f.remaining(), 10);
        f.consume(3);
        assert_eq!(f.remaining(), 7);
        f.consume(100); // saturates
        assert_eq!(f.remaining(), 0);
    }

    #[test]
    fn clone_shares_counter_decrements() {
        let f = Fuel::new(5);
        let g = f.clone();
        f.consume(2);
        // g observes the decrement made on f.
        assert_eq!(g.remaining(), 3);
        g.consume(3);
        assert!(f.check().is_err());
    }

    #[test]
    fn zero_budget_immediately_exhausted() {
        let f = Fuel::new(0);
        assert!(f.check().is_err());
        assert_eq!(f.remaining(), 0);
    }

    #[test]
    fn consume_zero_is_noop() {
        let f = Fuel::new(5);
        f.consume(0);
        assert_eq!(f.remaining(), 5);
        assert!(f.check().is_ok());
    }

    #[test]
    fn default_has_generous_budget() {
        let f = Fuel::default();
        assert!(f.remaining() > 0);
        assert!(f.check().is_ok());
    }
}
