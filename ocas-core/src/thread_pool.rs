//! Thread pool utilities.
//!
//! oCAS uses `rayon` for data parallelism. This module provides convenience
//! helpers and configuration defaults.

pub use rayon::prelude::*;

use crate::error::{OcasError, Result};

/// A scoped thread pool wrapper around `rayon::ThreadPool`.
///
/// Unlike the global pool, a `ThreadPool` can be created and destroyed
/// independently, which is useful for short-lived parallel tasks or tests.
pub struct ThreadPool {
    inner: rayon::ThreadPool,
}

impl ThreadPool {
    /// Create a new thread pool with the given number of worker threads.
    ///
    /// # Errors
    ///
    /// Returns `OcasError::BackendError` if rayon fails to build the pool
    /// (for example, because the global pool was already initialized with a
    /// different configuration).
    pub fn new(threads: usize) -> Result<Self> {
        let inner = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| OcasError::BackendError {
                backend: "rayon".into(),
                message: e.to_string(),
            })?;
        Ok(Self { inner })
    }

    /// Execute a closure in this thread pool.
    pub fn install<F, R>(&self, op: F) -> R
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        self.inner.install(op)
    }

    /// Return the number of worker threads in this pool.
    pub fn current_num_threads(&self) -> usize {
        self.inner.current_num_threads()
    }
}

/// Initialize the global thread pool with the default configuration.
///
/// This is normally called automatically by `rayon`, but explicit
/// initialization allows setting the number of threads.
///
/// # Errors
///
/// Returns `OcasError::BackendError` if the global pool was already
/// initialized with a different configuration.
pub fn init(threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|e| OcasError::BackendError {
            backend: "rayon".into(),
            message: e.to_string(),
        })
}

/// Execute a closure in the global thread pool.
pub fn install<F, R>(op: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    rayon::join(|| {}, op).1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_pool_install() {
        let result = install(|| 21 + 21);
        assert_eq!(result, 42);
    }

    #[test]
    fn custom_thread_pool() {
        let pool = ThreadPool::new(2).expect("should create a 2-thread pool");
        let result = pool.install(|| {
            let (a, b) = rayon::join(|| 1 + 2, || 3 + 4);
            a + b
        });
        assert_eq!(result, 10);
        assert_eq!(pool.current_num_threads(), 2);
    }

    #[test]
    fn thread_pool_maps_error_to_ocas_error() {
        // Initialize the global pool first so that building a new pool cannot
        // conflict. Then attempt to re-initialize with a different size, which
        // must fail.
        init(2).ok();
        let err = init(4).expect_err("reinitializing global pool should fail");
        assert_eq!(
            err,
            OcasError::BackendError {
                backend: "rayon".into(),
                message: "The global thread pool has already been initialized.".into(),
            }
        );
    }

    #[test]
    fn thread_pool_parallel_sum() {
        let pool = ThreadPool::new(2).expect("should create a 2-thread pool");
        let data: Vec<i64> = (0..1_000).collect();
        let sum = pool.install(|| data.into_par_iter().sum::<i64>());
        assert_eq!(sum, 499_500);
    }
}
