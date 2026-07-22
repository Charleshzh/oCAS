//! Core runtime, memory management, and backend glue for oCAS.
//!
//! This crate provides the foundation used by all other oCAS crates:
//! - arena allocation for expression nodes
//! - unified error types
//! - thread pool utilities
//! - FFI glue conventions
//! - thin wrappers around numerical backends (GMP, MPFR, FLINT, etc.)

#![warn(missing_docs)]

pub mod arena;
pub mod error;
pub mod fuel;
pub mod thread_pool;

/// Hash map using the DoS-resistant, high-performance `ahash` hasher.
///
/// Prefer this over [`std::collections::HashMap`] on hot paths
/// (polynomial term tables, hash-consing, CSE maps); std's SipHash is
/// needlessly slow for trusted internal keys.
pub type FastHashMap<K, V> = std::collections::HashMap<K, V, ahash::RandomState>;

/// Hash set using the `ahash` hasher. See [`FastHashMap`].
pub type FastHashSet<K> = std::collections::HashSet<K, ahash::RandomState>;

#[cfg(feature = "gmp")]
pub mod gmp;
