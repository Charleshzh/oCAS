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
///
/// # Example
///
/// ```
/// use ocas_core::thread_pool::ThreadPool;
///
/// let pool = ThreadPool::new(2).expect("failed to create pool");
/// let result = pool.install(|| {
///     let (a, b) = rayon::join(|| 1 + 2, || 3 + 4);
///     a + b
/// });
/// assert_eq!(result, 10);
/// ```
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
///
/// # Example
///
/// ```
/// use ocas_core::thread_pool::install;
///
/// let result = install(|| 21 + 21);
/// assert_eq!(result, 42);
/// ```
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

    mod simple {
        use super::*;

        #[test]
        #[cfg(not(miri))]
        fn global_pool_install() {
            let result = install(|| 21 + 21);
            assert_eq!(result, 42);
        }

        #[test]
        #[cfg(not(miri))]
        fn thread_pool_creation_and_size() {
            let pool = ThreadPool::new(2).expect("should create a 2-thread pool");
            assert_eq!(pool.current_num_threads(), 2);
        }
    }

    mod medium {
        use super::*;

        #[test]
        #[cfg(not(miri))]
        fn custom_thread_pool_join() {
            let pool = ThreadPool::new(2).expect("should create a 2-thread pool");
            let result = pool.install(|| {
                let (a, b) = rayon::join(|| 1 + 2, || 3 + 4);
                a + b
            });
            assert_eq!(result, 10);
        }

        #[test]
        #[cfg(not(miri))]
        fn thread_pool_parallel_sum() {
            let pool = ThreadPool::new(2).expect("should create a 2-thread pool");
            let data: Vec<i64> = (0..1_000).collect();
            let sum = pool.install(|| data.into_par_iter().sum::<i64>());
            assert_eq!(sum, 499_500);
        }
    }

    mod complex {
        use super::*;

        #[test]
        #[cfg(not(miri))]
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
        #[cfg(not(miri))]
        fn nested_pools_run_on_different_threads() {
            let outer = ThreadPool::new(2).expect("outer pool");
            let result = outer.install(|| {
                let inner = ThreadPool::new(2).expect("inner pool");
                inner.install(|| {
                    let (a, b) = rayon::join(|| 1 + 1, || 2 + 2);
                    a + b
                })
            });
            assert_eq!(result, 6);
        }
    }

    mod extreme {
        use super::*;

        #[test]
        #[cfg(not(miri))]
        fn stress_many_parallel_tasks() {
            let pool = ThreadPool::new(4).expect("should create a 4-thread pool");
            let result = pool.install(|| (0..10_000).into_par_iter().map(|x| x * x).sum::<i64>());
            // Sum of squares 0^2 + ... + (n-1)^2 = (n-1)n(2n-1)/6
            let n = 10_000i64;
            assert_eq!(result, (n - 1) * n * (2 * n - 1) / 6);
        }
    }
}
