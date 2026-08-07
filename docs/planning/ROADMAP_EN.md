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

## Phase 4: Closing the Competitive Gap (0.24–0.26) — COMPLETE

> **Goal**: close the key gaps found by the competitive survey
> (GAP_ANALYSIS_EN.md §5) — Symbolica 2.2 Rubi integration, msolve Gröbner
> performance, SymPy 1.14 DomainMatrix — before freezing 1.0.0.
>
> Background: Phase B++ "Competitive Alignment" (0.19–0.23) completed on
> 2026-08-02. Competitors evolved significantly meanwhile — Symbolica 2.2
> ported 7000+ Rubi integration rules, SymPy 1.14 DomainMatrix became
> 10000× faster, msolve set a cyclic-6 Gröbner benchmark of ~0.04 s — so the
> original "1.0 is freeze-and-polish only" plan was no longer sufficient.
>
> **(This phase is COMPLETE: 0.24 heuristic integration + DoubleF64, 0.25
> MultiModular Gröbner + parallel modular GCD, 0.26 packed monomial F5 fast
> channel; cyclic-6 ℤ₁₃ grevlex measured 55.04 ms. 0.26.0 shipped a different
> scope than originally planned here — the matrix engine / Smith normal form
> were deferred to 0.30.0 in Phase 5.)**

### 0.24.0 — Symbolic Integration Breadth + DoubleFloat

**Goal**: narrow the integration-coverage gap vs Symbolica Rubi (P0); introduce
the DoubleFloat evaluation path (P2).

**Deliverables**:

- [x] Heuristic integration pool behind Risch fallback (`heuristic_integrate`)
  - Integration by parts (LIATE/ILATE heuristic)
  - Trigonometric substitution ($\sqrt{a^2 - x^2}$, $\sqrt{a^2 + x^2}$, $\sqrt{x^2 - a^2}$)
  - Rational parameter substitution (Weierstrass $t = \tan(x/2)$)
  - Euler substitution (rationalising quadratic radicals, placeholder)
  - Reference: SymPy `manualintegrate` heuristic pool
- [x] DoubleFloat evaluation path (`DoubleF64`: ~31 digits, >3× faster than
  arbitrary precision)
  - Reference: Symbolica 2.0 `double-float` implementation
  - New `DoubleFloat` type in `ocas-domain`
  - JIT/SIMD evaluator DoubleFloat pipeline
- [x] Python/C bindings: `integrate_heuristic`, `DoubleFloat` type
- [ ] Rubi 1892-problem subset benchmark vs symbolica-integrate (deferred to 0.27.0)

**Success Criteria**:

- Rubi 1892-problem subset coverage improved ≥30 percentage points over the
  Risch-only baseline (Risch + heuristics)
- DoubleFloat evaluation ≥3× faster than arbitrary precision
- `cargo test --workspace` passes

### 0.25.0 — Gröbner Performance at Scale (Multi-Modular)

**Goal**: align Gröbner performance with msolve (P1); cyclic-6 ℤ₁₃ from 2.63 s
to < 0.5 s.

**Deliverables**:

- [x] Multi-modular strategy
  - Parallel Gröbner basis computation over several primes
  - CRT reconstruction of integer-coefficient bases
  - Rational reconstruction to recover ℚ coefficients
  - Reference: msolve F4 + multi-modular + Hensel + BM
- [x] Hensel lifting of Gröbner bases
  - Lift from $\mathbb{F}_p$ basis to $\mathbb{Z}$ basis
  - Fewer primes needed for CRT reconstruction
- [x] Large-coefficient polynomial GCD acceleration
  - Brown modular GCD further accelerated by multi-modular arithmetic
- [ ] Benchmarks: cyclic-6/7, katsura-6/7 vs msolve (katsura deferred to 0.28.0)

**Success Criteria**:

- cyclic-6 ℤ₁₃ < 0.5 s (was 2.63 s; msolve 0.04 s)
- cyclic-7 ℤ₁₃ tractable (previously untested)
- Benchmark results within one order of magnitude of msolve (< 10× gap)

### 0.26.0 — Packed F5 Fast Channel + grevlex Benchmarks (as shipped)

**Goal**: push the F5 main loop into a u128 SWAR fast channel, closing in on
msolve performance; add grevlex benchmark variants. (The originally planned
domain-aware matrix engine + Smith/Hermite normal forms were not shipped in
0.26.0 — deferred to 0.30.0.)

**Deliverables**:

- [x] Packed-monomial F5 fast channel (u128 SWAR)
  - Auto-routed when n_vars ≤ 8 and exponents < 2¹⁵; falls back to the generic
    path out of bounds
- [x] Echelon i32 / clone-free two-phase rework
- [x] grevlex benchmark variants (measurement baseline beyond Lex)
- [x] Fixed pre-existing Graded-order degree-direction inversion bug
- [ ] Domain-aware matrix engine (`DomainMatrix` analogue) → 0.30.0
- [ ] Smith/Hermite normal forms → 0.30.0
- [ ] Matrix performance benchmarks → 0.30.0
- [ ] Pre-1.0 freeze preparation (API audit / migration guide / cross-platform CI) → 0.30.0

**Success Criteria** (measured 2026-08-06):

- cyclic-6 ℤ₁₃ grevlex 52.07 ms (criterion median), Lex 936 ms
- cyclic-7 ℤ₁₃ grevlex single round 5.755 s (209 basis elements)
- Packed fast channel and generic path produce identical results (random cross-checks)

---

## Phase 5: Competitive Gap Closure (0.27–0.30)

> **Goal**: close the remaining P0–P3 gaps before freezing 1.0.0, per the
> 2026-08-06 priority re-ranking (GAP_ANALYSIS_EN.md §5): P0 symbolic
> integration breadth, P1 Gröbner performance at scale (katsura + cyclic-7),
> P1 LLVM JIT code generation, P2 matrix/linear algebra (DomainMatrix analogue
> + Smith/Hermite normal forms), P2 Windows FLINT, P3 quadratic sieve and
> tensor handling inside nested functions. Phase B+++ (0.24–0.26) delivered
> heuristic integration / DoubleF64, MultiModular Gröbner, and the packed F5
> fast channel (cyclic-6 grevlex 55.04 ms); this phase closes the rest, then
> 1.0.0 freezes.

### 0.27.0 — Symbolic Integration Breadth (Rubi-Grade Rule Set)

**Goal**: close the largest functional gap (P0) vs `symbolica-integrate`
(Rubi 7000+ rules, 72,944-problem corpus); lift the 1892-problem subset
coverage substantially.

**Deliverables**:

- [ ] Rule-table-driven integration engine (match → template substitution)
  - Power/polynomial/exponential/logarithm rule families
  - Trigonometric/hyperbolic/inverse-trigonometric/inverse-hyperbolic rule families
  - Radical and quadratic-form substitutions (extending the 0.24
    trig-substitution/Weierstrass/Euler framework, completing the Euler placeholder)
  - Special-function rule families (erf/Ei/Si/Ci/Fresnel, bridging the 0.14 table)
- [ ] Strategy dispatch chain: Risch (0.14) → heuristic four techniques (0.24)
  → rule library → `Integral(...)` fallback
- [ ] Rule provenance strategy (per GAP_ANALYSIS_EN.md §7.3 licence risk):
  - Preferred: self-developed rule set (Option C hybrid: Risch + heuristics +
    rule structure informed by Rubi's classification)
  - Evaluate integrating `symbolica-integrate` (MIT) as an optional feature
- [ ] 1892-problem coverage benchmark harness: coverage report + failure taxonomy
- [ ] Python/C bindings: `integrate` rule-path toggle

**Success Criteria**:

- 1892-problem subset coverage ≥30 percentage points above the current level
- Rule path agrees with SymPy `manualintegrate`/`integrate` on sampled cases
- `cargo test --workspace` passes

### 0.28.0 — Gröbner Performance at Scale (katsura + cyclic-7)

**Goal**: align with measured msolve 0.10.1 (katsura 3–7 ms, cyclic-7 55 ms)
(P1): katsura-6 < 1 s, cyclic-7 grevlex within one order of magnitude.

**Deliverables**:

- [ ] Extend the u128 packed F5 fast channel to katsura and cyclic-7
  (exponent-range / sparsity adaptation)
- [ ] Scale the MultiModular ℚ pipeline (0.25) to large instances
  - Parallel lucky-prime scheduling + CRT + rational reconstruction + traceless
    p-adic Hensel lifting
- [ ] Sparsity-aware echelon optimisation (successor of the 0.15.2 sparse
  echelon: row/column pruning)
- [ ] katsura-6/7 and cyclic-7 grevlex/Lex benchmarks vs measured msolve (WSL2)

**Success Criteria**:

- katsura-6 ℤ₁₃ < 1 s (currently not completed); katsura-7 tractable
- cyclic-7 grevlex within 10× of msolve (currently ~70×)
- Multi-modular path agrees with the single-prime path on 100 random cases;
  `is_groebner_basis` verified

### 0.29.0 — Code Generation Extension (LLVM/inkwell JIT)

**Goal**: land a second JIT backend — LLVM (via `inkwell`, already a workspace
dependency) — narrowing the code-generation gap vs Symbolica SymJIT (P1).

**Deliverables**:

- [ ] `ocas-eval::jit_llvm`: AST → LLVM IR + function registry + multi-output
- [ ] Evaluation pipeline coverage: f64/f32 mixed precision + DoubleF64 + SIMD
  vectorisation
- [ ] Runtime backend selection: Cranelift (default, fast compile) / LLVM
  (optimised code)
- [ ] Performance benchmarks: LLVM vs Cranelift vs interpreter (hold the
  multi-output 97×/21× baseline)
- [ ] Python/C bindings exposing the backend-selection parameter

**Success Criteria**:

- LLVM JIT on par with Cranelift or better; ≥10× vs interpreter maintained
- LLVM builds green on Linux/macOS/Windows CI
- Output identical to the Cranelift path (1000 random expressions)

### 0.30.0 — Matrix Engine + Platform Close-Out + 1.0 Freeze Preparation

**Goal**: close the P2/P3 gaps and finish pre-1.0 freeze preparation:
domain-aware matrix engine (DomainMatrix analogue) + Smith/Hermite normal
forms, Windows FLINT, quadratic sieve, tensor handling inside nested functions.

**Deliverables**:

- [ ] Domain-aware matrix engine (`DomainMatrix` analogue, deferred from 0.26.0)
  - `Matrix<D>` generic over `IntegerDomain`/`FiniteField`/`RationalDomain`
  - Domain-specialised paths for dense matrices (avoid generic `Domain` trait
    overhead)
  - Reference: SymPy DomainMatrix + FLINT backend
- [ ] Smith normal form (integer matrices; for module-structure analysis and
  homological algebra)
- [ ] Hermite normal form (integer matrices; for linear Diophantine equations)
- [ ] Matrix performance benchmarks: 20×20/30×30 integer rref/inv/det vs SymPy
  DomainMatrix
- [ ] Windows FLINT support (flint3-sys Windows build assessment + CI)
- [ ] Quadratic-sieve integer factorisation (vs SymPy `qs_factor`; next level
  above ECM)
- [ ] Tensor handling inside nested functions (vs Symbolica Graphica; 0.22
  delivered basic canonicalisation)
- [ ] Pre-1.0 freeze preparation
  - API audit: documentation completeness for all public types/functions
  - Migration guide finalised (Symbolica/SymPy → oCAS)
  - Cross-platform CI verification (Linux/macOS/Windows)
  - Published benchmarks (per BENCHMARK_SUITE_EN.md)

**Success Criteria**:

- Smith/Hermite normal forms agree with SymPy on 100 random cases
- 20×20 integer-matrix rref within one order of magnitude of SymPy DomainMatrix
- Quadratic-sieve benchmark recorded vs SymPy `qs_factor`
- Windows FLINT available on three platforms (or a documented hard blocker)
- Pre-1.0 freeze checklist ≥80% complete

---

## Phase 6: Stable 1.0

> **Goal**: A production-ready CAS library with stable APIs and broad backend
> support.

### 1.0.0 — Stable Release

**Target**: after 0.30.0

**Deliverables**:

- [ ] Stable semantic versioning guarantee
- [ ] Full Rust, Python, and C/C++ API coverage
- [ ] Comprehensive test suite (>80% line coverage)
- [ ] Published benchmarks
- [ ] Migration guide from Symbolica/SymPy
- [ ] Signed release artifacts
- [ ] Competitive comparison report (per the final COMPETITIVE_MATRIX_EN.md)

**Success Criteria**:

- No breaking API changes planned for 1.x.
- P0 gap (symbolic integration breadth) significantly narrowed (1892-problem
  subset coverage target met, GAP_ANALYSIS_EN.md §5).
- P1 gap (Gröbner performance) aligned with msolve within one order of
  magnitude (katsura-6 < 1 s, cyclic-7 grevlex < 10× msolve).
- P1 gap (code generation) closed with an LLVM/inkwell JIT backend.
- P2 gap (matrix/linear algebra) closed with a domain-aware matrix engine +
  Smith/Hermite normal forms.
- Performance ahead of SymPy across core benchmarks.

> The fine-grained per-version plan from Beta to 1.0 (0.11 factorization →
> 0.12 rational functions → 0.13 Gröbner F4 → 0.14 Risch integration → 0.15
> multi-output JIT → 0.15.2 Gröbner performance at scale → 0.16 arbitrary
> multivariate factorization → 0.16.1 non-constant leading-coefficient
> imposition → 0.17 algebraic-number-field factorization →
> 0.18 numerical integration / duals / tensors / fuel) is detailed in
> [EVOLUTION_PLAN_EN.md](EVOLUTION_PLAN_EN.md). Versions 0.15.2–0.18.0 form
> Phase B+ "Closing the Symbolica Gap" (complete); 0.19–0.23 form Phase B++
> "Competitive Alignment" (F5 Gröbner → ODE solvers → number theory → tensor
> canonicalisation → algebraic geometry); 0.24–0.26 form Phase B+++
> "Competitive Gap Bridging" (complete); 0.27–0.30 form Phase B++++
> "Competitive Gap Closure" (P0–P3). After Phase B++++, 1.0.0 freezes.

---

## Post-1.0

After 1.0, development will focus on:

- Partial differential equation (PDE) solvers (Poisson, heat, wave)
- Differential Galois theory (research prelude)
- Optional GPL backends (`ocas-gpl`)
- GPU acceleration + code export (CUDA / HIP / Vulkan compute, CUDA/WASM)
- Domain-specific toolkits (physics, robotics, machine learning)

> The LLVM/Inkwell JIT backend moved to 0.29.0, and the quadratic-sieve integer
> factorisation moved to 0.30.0 (both pre-1.0).

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
| 0.13.1 | Beta | Month 19 | docs.rs build fix |
| 0.13.2 | Beta | Month 19 | PyPI release (`pip install ocas`) + dependency upgrades + CI hardening |
| 0.14.0 | 1.0 Candidate | Month 22 | Risch symbolic integration + rational-function integration + special-function table + FGLM/F5/Hilbert + trigonometric integration |
| 0.15.0 | 1.0 Candidate | Month 24 | Multi-output JIT + f32 mixed precision + streaming evaluation + Arena/workspace pool + ahash + native i64 F4 |
| 0.15.1 | 1.0 Candidate | Month 24 | F4 real linear algebra fix (cyclic-5 ~85,000× faster, cyclic-6 tractable) |
| 0.15.2 | 1.0 Candidate | Month 25 | Gröbner performance at scale (LM index + sparse echelon, cyclic-6 ℤ₁₃ 9970 s → 3670 s; <5 s needs F5) |
| 0.16.0 | 1.0 Candidate | Month 26 | Arbitrary multivariate factorization (Wang EEZ, ≥3 variables, ℤ and ℤ_p) ✅ |
| 0.16.1 | 1.0 Candidate | Month 26 | Non-constant leading-coefficient imposition (mod-p Hensel) + multivariate sparsity + sparse Diophantine ✅ |
| 0.16.2 | 1.0 Candidate | Month 26 | $\mathbb{F}_p$-path non-constant LC preprocessing (field Wang) + sampling performance |
| 0.17.0 | 1.0 Candidate | Month 27 | Algebraic number field & extension-field factorization (Trager) ✅ (univariate path; multivariate extension deferred) |
| 0.18.0 | 1.0 Candidate | Month 28 | Numerical integration (Vegas) + dual-number AD + tensor basics + fuel resource control |
| 0.18.1 | 1.0 Candidate | Month 28 | Python/C bindings backfill for the three 0.18.0 capabilities (numeric integration + tensor + dual) + prelude completeness ✅ |
| 0.19.0 | 1.0 Candidate | Month 30 | F5 Gröbner basis signature reduction (cyclic-6 ℤ₁₃ <5 s target) ✅ (2.63 s, ~1400×; multi-order deferred to 0.19.1) |
| 0.20.0 | 1.0 Candidate | Month 33 | Ordinary differential equation solvers (5 first-order + 2 second-order + power-series framework + classifier) ✅ (core complete; Laplace/systems/bindings deferred) |
| 0.20.1 | 1.0 Candidate | Month 33 | ODE backfill: integrating factors + VOP + reduction of order + series recursion + Frobenius + Laplace IVP + 2×2 systems + Python/C bindings + 31 substitution-verified tests ✅ |
| 0.21.0 | 1.0 Candidate | Month 36 | Number theory & computational algebra (modular GCD + integer factorization + primality + discrete log + CRT + number-theoretic functions) ✅ (incl. Python/C bindings; ECM factors 30-digit semiprimes in 1.1 s) |
| 0.22.0 | 1.0 Candidate | Month 39 | Tensor canonicalisation (graph-isomorphism engine) + advanced pattern matching (`Transformer::Partition`) ✅ |
| 0.23.0 | 1.0 Candidate | Month 42 | Advanced Gröbner & algebraic-geometry tooling (ideal ops + RUR + primary decomposition + Hilbert series) ✅ |
| 0.24.0 | Beta | Month 45 | Symbolic integration breadth (heuristic expansion) + DoubleFloat evaluation path (P0 integration + P2 DoubleFloat) ✅ |
| 0.25.0 | Beta | Month 47 | Gröbner performance at scale (multi-modular vs msolve, cyclic-6 < 0.5 s) (P1) ✅ |
| 0.26.0 | Beta | Month 49 | Packed-monomial F5 fast channel + grevlex benchmarks (cyclic-6 grevlex 55.04 ms measured) ✅ |
| 0.27.0 | Beta | Month 51 | Symbolic integration breadth (Rubi-grade rule set + 1892-problem coverage benchmark) (P0) |
| 0.28.0 | Beta | Month 53 | Gröbner performance at scale (katsura-6 < 1 s, cyclic-7 within 10× of msolve) (P1) |
| 0.29.0 | Beta | Month 55 | Code generation extension (LLVM/inkwell JIT backend) (P1) |
| 0.30.0 | Beta | Month 57 | Matrix engine (DomainMatrix analogue + Smith/Hermite) + Windows FLINT + quadratic sieve + 1.0 freeze preparation (P2/P3) |
| 1.0.0 | Stable | Month 59 | Stable release (frozen after Phase B++++ competitive gap closure: P0 integration breadth met + P1 Gröbner aligned with msolve + LLVM JIT landed + performance ahead of SymPy) |

---

## How to Read This Roadmap

- Each version represents a **potentially publishable** increment.
- Dates are approximate and depend on contributor availability.
- Features may shift between versions based on user feedback and technical discoveries.

---

## Contributing to the Roadmap

If you want to work on a specific version or feature, please open a GitHub issue
and we will assign a tracking issue to you.
