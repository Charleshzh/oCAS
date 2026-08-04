//! Shared test helpers and polynomial-system generators for oCAS.
//!
//! Integration tests and benchmarks under `tests/` and `benches/` use
//! these generators so that identical systems are exercised everywhere
//! (cyclic-n and Katsura-n over ℚ and ℤ_p).

pub mod systems;
