# oCAS Roadmap

This document outlines the development roadmap of oCAS from pre-alpha
experiments to a stable 1.0 release, with each 0.x version carrying concrete
deliverables. For the Chinese edition, see [ROADMAP_CN.md](ROADMAP_CN.md).
Companion documents: [EVOLUTION_PLAN_EN.md](EVOLUTION_PLAN_EN.md) (fine-grained
per-version plan) and [GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md) (gap snapshot).

---

## Legend

| Tag | Meaning |
|---|---|
| `API` | Public API surface |
| `ALG` | Algebraic algorithms |
| `NUM` | Numerical backends |
| `PERF` | Performance and optimization |
| `BIND` | Language bindings |
| `DOC` | Documentation and examples |
| `TEST` | Testing and quality |

---

## Phase 1: Pre-Alpha — Foundation

> **Goal**: Establish the workspace, runtime, and basic expression core. Prove
> that the layered architecture compiles and runs.

### 0.1.0 — Workspace & Runtime

**Target**: Month 1

**Deliverables**:

- [x] Workspace structure with all 12 crates
- [x] CI pipeline: `cargo test`, `cargo clippy`, `cargo-deny`, formatting, Miri
- [x] Unified error type `OcasError`
- [x] Arena / bump allocator with Miri-safe API
- [x] Thread pool wrapper around `rayon`
- [x] FFI glue conventions (minimal C ABI example)
- [x] GMP bindings (via `rug`) behind `gmp` feature
- [x] Initial benchmark harness

**Success Criteria**:

- `cargo build --workspace` succeeds on Linux/macOS/Windows (no-default-features on MSVC).
- Arena passes Miri and valgrind/ASan checks.
- GMP integer arithmetic is callable from Rust on supported platforms.

### 0.2.0 — Expression Tree Core

**Target**: Month 2

**Deliverables**:

- [x] `ocas-atom` crate
- [x] `Atom` tagged-union design
- [x] Arena-backed AST with safe public API
- [x] Hash consing for common subexpressions
- [x] Lexer using `logos`
- [x] Recursive-descent / Pratt parser
- [x] Printer: ASCII and compact forms
- [x] Normalizer: flatten `Add`/`Mul`, sort terms, merge coefficients

**Success Criteria**:

- `parse("x^2 + 2*x + 1")` produces the expected AST.
- `to_string(parse(s)) == s` for a broad set of expressions.
- Normalization is deterministic and property-tested.

---

## Phase 2: Alpha — Symbolic Engine

> **Goal**: A usable Rust API for parsing, simplification, differentiation, and
> basic polynomial operations.

### 0.3.0 — Domains & Polynomials

**Target**: Month 4

**Deliverables**:

- [x] `ocas-domain` crate
- [x] Domains: `Integer`, `Rational`, `FiniteField`
- [x] Domain trait for generic algorithms
- [x] `ocas-poly` crate
- [x] Dense univariate polynomial
- [x] Domains: `RealBall`, `Complex`
- [x] Sparse multivariate polynomial
- [x] Division with remainder
- [x] FLINT 3 integration behind `flint` feature
- [x] Optional GMP backend for `Integer`/`Rational` via `rug`
- [x] Optional MPFR backend for `RealBall` via `rug`

  > **Note**: The `flint` feature is experimental. It builds and runs on
  > Linux/WSL where system FLINT is available, but it is not yet supported on
  > Windows because `flint3-sys` depends on POSIX-only types such as
  > `pthread_mutex_t`. The default recommended Windows backend for
  > arbitrary-precision integers, rationals, and rigorous real balls is
  > `gmp`/`mpfr` via `rug` with system GMP/MPFR installed through MSYS2.

**Success Criteria**:

- Polynomial operations match SymPy outputs on regression suite.
- FLINT path produces identical results to pure-Rust fallback for supported operations.

### 0.4.0 — Pattern Matching & Rewriting

**Target**: Month 5

**Deliverables**:

- [x] Pattern matching engine with wildcards and conditions
- [x] `Transformer` visitor API
- [x] Basic built-in rewrite rules
- [x] `egg` integration for equality saturation
- [x] Rule-based simplifier

**Success Criteria**:

- Common identities (e.g., `x + x -> 2*x`, `x * 0 -> 0`) are applied automatically.
- E-graph can simplify `sin(x)^2 + cos(x)^2` to `1` under assumptions.

### 0.5.0 — Calculus Basics

**Target**: Month 6

**Deliverables**:

- [x] Symbolic differentiation
- [x] Derivative table for elementary functions
- [x] Taylor series expansion
- [x] Partial integration with heuristic table
- [x] `ocas-calc` crate initial release

**Success Criteria**:

- Differentiation passes a comprehensive test suite.
- Integration succeeds on standard calculus problems.

### 0.6.0 — First Rust API Release Candidate

**Target**: Month 7

**Deliverables**:

- [x] Stable `ocas` prelude
- [x] Rustdoc examples for all public APIs
- [x] Property tests with `proptest`
- [x] Initial benchmark suite
- [x] SymPy comparison harness via `uv`
- [x] crates.io publish preparation (internal workspace deps versioned)

**Success Criteria**:

- `cargo test --workspace --exclude ocas-py` passes.
- Benchmarks demonstrate parity with SymPy on basic polynomial, calculus, and rewriting operations.
- `cargo publish --dry-run -p ocas-core` succeeds; top-level `ocas` is ready to publish once the internal crates are uploaded.

---

## Phase 3: Beta — Solvers, JIT, Bindings

> **Goal**: Multi-language availability and performance. Core algebra is
> feature-complete for a CAS beta.

### 0.7.0 — Equation Solvers

**Target**: Month 9

**Deliverables**:

- [x] Linear system solver (`faer` / `LinBox`)
- [x] Polynomial system solver (Gröbner + root isolation)
- [x] Single-variable root finding via Arb
- [x] Diophantine solver basics
- [x] Assumptions / domain system

**Success Criteria**:

- Linear and polynomial solvers produce correct results verified against SageMath.

### 0.8.0 — Evaluation & JIT

**Target**: Month 11

**Deliverables**:

- [x] Tree interpreter for scalar and vector evaluation
- [x] AST-to-instruction compiler
- [x] Function registry for user-defined functions
- [x] Cranelift JIT backend
- [x] SIMD vectorized evaluation

**Success Criteria**:

- JIT evaluates repeated expressions at least 10x faster than interpreter.
- SIMD path works for dense polynomial evaluation.

### 0.9.0 — Python & C/C++ Bindings

**Target**: Month 13

**Deliverables**:

- [x] `ocas-py` crate with PyO3
- [~] Python classes: `Expression` (done), `Polynomial`/`Matrix`/`Domain` (deferred to 0.10.0)
- [x] Maturin wheel build for Linux/macOS/Windows
- [x] `ocas-c` crate with cbindgen
- [x] Stable C API for expression lifecycle
- [x] C++ RAII wrapper

**Success Criteria**:

- `pip install ocas` works on supported platforms.
- C example compiles and runs against the shared library.
- No memory leaks in binding tests (tracemalloc + RAII-guarded arenas).

### 0.10.0 — Beta Release

**Target**: Month 14

**Deliverables**:

- [x] Python classes deferred from 0.9.0: `Polynomial`, `Matrix`, `Domain`
- [x] Feature freeze for 1.0
- [x] Comprehensive documentation site
- [x] Performance comparison with Symbolica and SageMath
- [x] Community feedback integration
- [x] Bug-fix only period

**Success Criteria**:

- All public APIs documented.
- CI green on all supported platforms.

---

## Phase 4: Stable 1.0

> **Goal**: A production-ready CAS library with stable APIs and broad backend
> support.

### 1.0.0 — Stable Release

**Target**: Month 16

**Deliverables**:

- [ ] Stable semantic versioning guarantee
- [ ] Full Rust, Python, and C/C++ API coverage
- [ ] Comprehensive test suite (>80% line coverage)
- [ ] Published benchmarks
- [ ] Migration guide from Symbolica/SymPy
- [ ] Signed release artifacts

**Success Criteria**:

- No breaking API changes planned for 1.x.
- Performance parity or better with Symbolica on core benchmarks.

> The fine-grained per-version plan from Beta to 1.0 (0.11 factorization →
> 0.12 rational functions → 0.13 Gröbner F4 → 0.14 Risch integration → 0.15
> multi-output JIT) is detailed in [EVOLUTION_PLAN_EN.md](EVOLUTION_PLAN_EN.md).

---

## Post-1.0

After 1.0, development will focus on:

- Advanced symbolic integration (Risch algorithm)
- Differential equation solvers
- Optional GPL backends (`ocas-gpl`)
- GPU acceleration (CUDA / HIP / Vulkan compute)
- LLVM/Inkwell JIT backend
- Domain-specific toolkits (physics, robotics, machine learning)

---

## Milestones

| Version | Phase | Target | Key Deliverable |
|---|---|---|---|
| 0.1.0 | Pre-Alpha | Month 1 | Workspace + runtime |
| 0.2.0 | Pre-Alpha | Month 2 | Expression core |
| 0.3.0 | Alpha | Month 4 | Domains & polynomials |
| 0.4.0 | Alpha | Month 5 | Pattern matching & rewriting |
| 0.5.0 | Alpha | Month 6 | Calculus basics |
| 0.6.0 | Alpha | Month 7 | Rust API RC |
| 0.7.0 | Beta | Month 9 | Equation solvers |
| 0.8.0 | Beta | Month 11 | JIT & evaluation |
| 0.9.0 | Beta | Month 13 | Python & C/C++ bindings |
| 0.10.0 | Beta | Month 14 | Feature freeze |
| 0.11.0 | Beta | Month 15 | Polynomial factorization (univariate) |
| 0.11.1 | Beta | Month 15 | Polynomial factorization (bivariate + bindings + docs) |
| 0.11.2 | Beta | Month 16 | Computation acceleration (SOO Integer, mimalloc, modular GCD) |
| 0.12.0 | Beta | Month 17 | Rational polynomials + resultant + partial fractions + Karatsuba + rational reconstruction |
| 0.13.0 | Beta | Month 19 | Gröbner F4 matrix algorithm |
| 0.14.0 | 1.0 Candidate | Month 22 | Risch symbolic integration |
| 0.15.0 | 1.0 Candidate | Month 24 | Performance + JIT + memory optimization |
| 1.0.0 | Stable | Month 26 | Stable release |

---

## How to Read This Roadmap

- Each version represents a **potentially publishable** increment.
- Dates are approximate and depend on contributor availability.
- Features may shift between versions based on user feedback and technical discoveries.

---

## Contributing to the Roadmap

If you want to work on a specific version or feature, please open a GitHub issue
and we will assign a tracking issue to you.
