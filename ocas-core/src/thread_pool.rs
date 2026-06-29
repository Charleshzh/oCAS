//! Thread pool utilities.
//!
//! oCAS uses `rayon` for data parallelism. This module provides convenience
//! helpers and configuration defaults.

pub use rayon::prelude::*;

/// Initialize the global thread pool with the default configuration.
///
/// This is normally called automatically by `rayon`, but explicit
/// initialization allows setting the number of threads.
pub fn init(threads: usize) -> Result<(), rayon::ThreadPoolBuildError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
}

/// Execute a closure in the global thread pool.
pub fn install<F, R>(op: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    rayon::join(|| {}, op).1
}
