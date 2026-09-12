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
> were deferred to 0.35.0 in Phase 5.)**

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
- [ ] Benchmarks: cyclic-6/7, katsura-6/7 vs msolve (katsura deferred to 0.33.0)

**Success Criteria**:

- cyclic-6 ℤ₁₃ < 0.5 s (was 2.63 s; msolve 0.04 s)
- cyclic-7 ℤ₁₃ tractable (previously untested)
- Benchmark results within one order of magnitude of msolve (< 10× gap)

### 0.26.0 — Packed F5 Fast Channel + grevlex Benchmarks (as shipped)

**Goal**: push the F5 main loop into a u128 SWAR fast channel, closing in on
msolve performance; add grevlex benchmark variants. (The originally planned
domain-aware matrix engine + Smith/Hermite normal forms were not shipped in
0.26.0 — deferred to 0.35.0.)

**Deliverables**:

- [x] Packed-monomial F5 fast channel (u128 SWAR)
  - Auto-routed when n_vars ≤ 8 and exponents < 2¹⁵; falls back to the generic
    path out of bounds
- [x] Echelon i32 / clone-free two-phase rework
- [x] grevlex benchmark variants (measurement baseline beyond Lex)
- [x] Fixed pre-existing Graded-order degree-direction inversion bug
- [ ] Domain-aware matrix engine (`DomainMatrix` analogue) → 0.35.0
- [ ] Smith/Hermite normal forms → 0.35.0
- [ ] Matrix performance benchmarks → 0.35.0
- [ ] Pre-1.0 freeze preparation (API audit / migration guide / cross-platform CI) → 0.35.0

**Success Criteria** (measured 2026-08-06):

- cyclic-6 ℤ₁₃ grevlex 52.07 ms (criterion median), Lex 936 ms
- cyclic-7 ℤ₁₃ grevlex single round 5.755 s (209 basis elements)
- Packed fast channel and generic path produce identical results (random cross-checks)

---

## Phase 5: Competitive Gap Closure and Mechanism Push (0.27–0.35)

> **Goal**: close the remaining P0–P3 gaps before freezing 1.0.0, per the
> 2026-08-06 priority re-ranking (GAP_ANALYSIS_EN.md §5): P0 symbolic
> integration breadth, P1 Gröbner performance at scale (katsura + cyclic-7),
> P1 LLVM JIT code generation, P2 matrix/linear algebra (DomainMatrix analogue
> + Smith/Hermite normal forms), P2 Windows FLINT, P3 quadratic sieve and
> tensor handling inside nested functions. Phase B+++ (0.24–0.26) delivered
> heuristic integration / DoubleF64, MultiModular Gröbner, and the packed F5
> fast channel (cyclic-6 grevlex 55.04 ms); this phase closes the rest, then
> 1.0.0 freezes.
>
> **2026-09-12 re-plan (after the 0.27.3 close-out)**: the 0.27.x failure
> attribution (BENCHMARK_RESULTS_CN.md §"0.27.3 follow-up") and the
> general-mechanism feasibility assessment
> (GENERAL_MECHANISM_FEASIBILITY_EN.md) **redefined the nature of P0**.
> 0.27.0's original acceptance line ("Rubi-grade rule set + 30 pp on the 1892")
> was a **rule-enumeration** metric, whereas the marginal return of six waves of
> mechanism work fell from +6.82 pp to +1.11 pp, and only 0.7% of the remaining
> 1522 unsolved cases share a shape skeleton with a solved one. Therefore:
>
> - **0.28.0–0.32.0 become the mechanism push**: first fix the two measured
>   architectural defects and establish the "symbolic certificate + three-valued
>   outcome" correctness discipline, then complete the transcendental Risch, the
>   algebraic extensions, the non-elementary layer and complexity, following
>   Bronstein's ladder;
> - **the former 0.28.0 / 0.29.0 / 0.30.0 slide in order to
>   0.33.0 / 0.34.0 / 0.35.0** (Gröbner performance at scale, LLVM JIT, matrix
>   engine + platform close-out + 1.0 freeze preparation);
> - 1.0.0 slides to Month 69.

### 0.27.0 — Symbolic Integration Breadth (Rubi-Grade Rule Set)

**Goal**: close the largest functional gap (P0) vs `symbolica-integrate`
(Rubi 7000+ rules, 72,944-problem corpus); lift the 1892-problem subset
coverage substantially.

**Deliverables**:

- [x] Rule-table-driven integration engine (match → template substitution)
  - Power/polynomial/exponential/logarithm rule families
  - Trigonometric/hyperbolic/inverse-trigonometric/inverse-hyperbolic rule families
  - Radical and quadratic-form substitutions (extending the 0.24
    trig-substitution/Weierstrass/Euler framework, completing the Euler placeholder)
  - Special-function rule families (erf/Ei/Si/Ci/Fresnel, bridging the 0.14 table)
- [x] Strategy dispatch chain: Risch (0.14) → heuristic four techniques (0.24)
  → rule library → `Integral(...)` fallback
- [x] Rule provenance strategy (per GAP_ANALYSIS_EN.md §7.3 licence risk):
  - Preferred: self-developed rule set (Option C hybrid: Risch + heuristics +
    rule structure informed by Rubi's classification)
  - Evaluate integrating `symbolica-integrate` (MIT) as an optional feature
- [x] 1892-problem coverage benchmark harness: coverage report + failure taxonomy
- [x] Python/C bindings: `integrate` rule-path toggle
- [x] (phase-(c) additions) symbolic-constant rational integrator, Weierstrass
  linear arguments, bounded distributive expansion retry, trig
  product-to-sum/power reduction, subresultant dense GCD, global chain-entry
  budget for the integration pipeline, exact real-root isolation
  (Wilkinson 10/10)

**Success Criteria** (honest record):

- 1892-problem subset coverage ≥30 percentage points above the current level
  — **NOT met**: baseline 5.87% → 7.66% after phase (b) → see
  BENCHMARK_RESULTS_CN.md for phase (c); root causes quantified
  (Rubi-grade rule volume / symbolic-coefficient rational backend / nested
  radicals all need much larger investment); alternative paths evaluated in
  GAP_ANALYSIS §7.3
- Rule path agrees with SymPy `manualintegrate`/`integrate` on sampled cases — met
- `cargo test --workspace` passes — met

### 0.27.1 — Integration-Breadth Mechanism Push (0.27.0 acceptance-line sequel)

**Goal**: continue raising 1892-problem coverage along the 0.27.0 +30pp
acceptance line, via mechanism-level upgrades (not rule volume).

**Deliverables** (all landed):

- [x] Chebyshev binomial differentials + fractional-power rationalization
  (`binomial.rs`)
- [x] Trig-denominator power reductions / linear-numerator decomposition /
  polynomial×trig closed forms (`trig_reduction.rs`)
- [x] exp/log kernel substitutions (exp-kernel rationalization, hyperbolic
  t=e^u, f(log x)/x) (`exp_log.rs`)
- [x] General sqrt-quadratic engine + Euler III (`sqrt_quadratic.rs`)
- [x] Inverse-trig/hyperbolic kernel-derivative powers and substitutions
  (`inverse_trig.rs`)
- [x] Single-trig-kernel rational normalization + tan/sec-family reductions
  (`trig_kernel.rs`)
- [x] Denominator-power recurrences + two-linear-factor partial fractions
  (`quad_power.rs`)
- [x] Wrong-answer fixes: C14/D7b linear-argument reduction residual
  coefficients; rational.rs √(p/q) dropping 1/q
- [x] Stability: symbolic-rational coefficient budgets + many-symbol entry
  gate (timeouts 49→33, wall clock −19%), deep-residue checks in the
  heuristic stage, expansion moved before it to stop budget starvation

**Success Criteria** (honest record):

- 1892-problem subset +30pp — **NOT met**: 9.62% → 16.44% (+6.82pp,
  311/1892, 129 newly solved, zero regressions, zero crashes; the
  quantified gap — elliptic families / high-symbol-count quotients /
  composite shells — is recorded in BENCHMARK_RESULTS_CN.md §0.27.1)
- Every new mechanism verified by eval_f64 numeric differentiation sampling
  plus the SymPy cross-check suite — met
- Quality gates (fmt / clippy -D warnings / workspace test / deny) — met

### 0.27.2 — Hang Elimination, Verified Coverage and the Elliptic Foundation

**Goal**: turn every remaining per-case hang into a deterministic decline,
add an independent numerical-verification criterion alongside the string
coverage metric, close the elementary mechanism families the failure dump
points at, and lay the elliptic-integral foundation.

**Deliverables** (all landed):

- [x] Diagnostics: `OCAS_INTEGRATE_TRACE` stage tracing; all 33 baseline
  timeouts attributed to pipeline stages (symbolic_rational 18, heuristic 4,
  trig_kernel 3, inverse_trig 2, rational 2, trig_reduction 2,
  sqrt_quadratic 1, untraced 1)
- [x] Bounded-expansion pre-pass ahead of the symbolic-rational backend
  (`expand_prepass` in `integral/mod.rs`) for unexpanded products whose
  expansion is a small Laurent polynomial (restricted to purely algebraic
  expansions so it cannot steal trig/hyperbolic cases from their mechanisms)
- [x] Deterministic work budgets in the hanging stages (`symbolic_rational`,
  `trig_kernel`, `heuristic`, `rational`, `sqrt_quadratic`)
- [x] Wrong-answer fix: `symbolic_rational::rational_square_root` treated a
  polynomial **sum** as a monomial square, splitting quadratic denominators at
  bogus roots and emitting wrong logarithms (0.27.1 shipped wrong answers for
  e.g. `1/(b*x^2+2*a*x-b)` and, through the `t = e^x` hyperbolic path,
  `1/(a+b*sinh(x))`); five further wrong-answer classes were fixed
  (`complete_square` monic assumption, Chebyshev case-3 branch sign, Risch
  out-of-field residuals, rule-A4 sequence wildcard, resonant product-to-sum
  zero denominators)
- [x] Numerical-verification oracle (`ocas-tests/src/integral_eval.rs`) and
  `verified_solved` / `unverified_solved` / `verify_indeterminate` reporting
  in the 1892 harness, plus `OCAS_1892_VERIFY` / `OCAS_1892_BUCKET` switches
- [x] Hyperbolic closed-form family (`hyperbolic_reduction.rs`)
- [x] Rational-derivative kernel substitution (`kernel_subst.rs`)
- [x] Trig phase-shift normalization (`trig_reduction.rs` extension; the
  phase/ratio path is implemented but gated behind
  `PHASE_RATIO_ENABLED = false` because its composed atom does not yet survive
  `crate::diff`)
- [x] Inverse-composition cancellation (`inverse_trig.rs` extension)
- [x] `exp(inverse function)` algebraization (`exp_log.rs` extension)
- [x] Half-power front-end (`halfpower.rs`) and elliptic reduction
  (`elliptic.rs`): Legendre reduction to `EllipticF`/`EllipticE`
  (`EllipticPi` is registered but has no producer yet), SymPy's `m = k²`
  convention, argument order preserved by
  `ocas_atom::normalize::preserves_argument_order`
- [x] Wrong-answer regression guard in the correctness suite
  (`tests/correctness/integral_verify.rs`)
- [x] Extra: `ocas-parse` unary minus (`2-1`, `x-1`, `x^2-1` no longer
  `PARSE_ERR`; breaking `lex` API change)

**Success Criteria** (honest record):

- Both metrics over the 1892-problem subset: **solved 349 (18.45%), verified
  325/349 (93.1%)**; per-case diff: 44 newly solved, 6 deliberate declines
  (all were 0.27.1 wrong answers), net +38
- Zero wrong answers in the solved set: **`verify_mismatches` = 0 — met**
- Timeouts 13 (target ≤ 5 — **not met**), crashes 0, wall clock 445.7 s
  (−22% against the 572 s baseline)
- Verified ratio 93.1% (target ≥ 95% — **not met**: of 24 inconclusive cases,
  11 are domain-restricted, 10 carry the imaginary unit, 3 fail the step-size
  self-consistency gate; all are "undecidable", none is a wrong answer)
- Quality gates: fmt / both clippy tiers / deny green; workspace tests green
  (the `ocas-c` rules-toggle probe now uses `csc(x)^5`, since kernel
  substitution took over the old `tan(x)^4` probe)
- Next wave (0.27.3): the 13 timeout families, composite-shell decomposition,
  the special-function family, elliptic-family breadth

### 0.27.3 — Composite Shells, Special-Function Breadth and Elliptic Family Coverage

**Goal**: continue the mechanism push on the largest clusters the 0.27.2 failure dump
leaves behind.

**Deliverables** (each line records what was actually measured):

- [x] Prerequisite: special-function derivative table (`derivative.rs`) and oracle heads
  (`Ei`/`Ei(n,z)`/`Si`/`Ci`/`Shi`/`Chi`/`fresnels`/`fresnelc`). Every algorithm was
  cross-checked against `mpmath` at 40 digits before being written down; an unimplemented
  partial (the order slot of `Ei(n,z)`) declines to an unevaluated `Derivative` instead of
  silently dropping a term
- [x] Special-function families (`integral/special.rs`): four reduction families —
  polynomial × `F(a+bx)`, polynomial × `Ei(n,a+bx)`, `F(bx)/x^m` descent, polynomial ×
  `F(a+bx)²`; step budgets `MAX_SPECIAL_STEPS=16`, `MAX_SPECIAL_DEG=8`; every residual is
  closed-form, so no family ever emits `Integral(...)`. **special bucket: 15 of 19 solved
  (all numerically verified), 4 honestly declined** (two are `Unintegrable` in Rubi itself,
  two are Fresnel denominator chains deliberately not implemented)
- [x] Half-power affine arguments (`halfpower.rs`): accepts `cos(c + d·x)`, emits the sheet
  factor on both branches (`cos(u/2)·(1−sin²(u/2))^(−1/2)`), and folds same-base
  `S^a·(√S)^b` into a single half power first. **4 newly solved**
  (`rubi-00377/00445/01598/00291`; all four read `indeterminate` at the corpus oracle only
  because their parameter is pinned at `m=2`, outside its `|m| ≤ 0.9` window — the module's
  own independent numeric oracle verifies them at 1e-9)
- [x] Exact linear-square fold (`integral/mod.rs`): `p²+2pq+q² → (p+q)²`, accepted only when
  the base is **affine in the integration variable** and re-expansion reproduces the input
  sum; never applied under a non-integer power (`((p+q)²)^{1/2} = |p+q| ≠ p+q`). **2 newly
  solved**, and `rubi-00854` goes from a 10 s per-case timeout to a 0.01 s solve
- [ ] Composite-shell decomposition (mixed-other cluster): **partly delivered, partly
  withdrawn**. The concrete implementation of factor-level `F(g(x))·w(x)` splitting — a
  log-of-a-linear-fraction reduction (`C·(A+B·log Q)·(F+Gx)^p` by parts) — exposed a real
  wrong answer once integrated (its residual was returned unwrapped), and its residual was
  not reliably integrable at an acceptable wall-clock cost; per the "keep only what adds no
  regressions and no wall-clock cost" rule it was **withdrawn wholesale**
  (`log_fraction.rs` deleted). The attempt did locate and harden the A4 regression guard
  (from "must leave a residue" to "if solved, check the derivative numerically" — neither
  outcome lets a wrong answer through). The `trig_linear_arg` direction of this item is
  delivered by the affine-argument work above
- [ ] Elliptic-family breadth: full trig-quadratic routing, cubic radicands, complex
  `EllipticPi` — **not achieved, recorded honestly**. 0.27.3 delivered only the affine
  argument layer; 136 of the 140-case affine-cos half-power cluster still decline, with a
  quantified root cause: they need **products/quotients of two different radical bases**, a
  rational prefactor in `cos u`, `|p| > 5`, or positive powers ≥ 3 — i.e. a multi-radical
  rational-prefactor reduction engine (a new engine, not this wave). Cubic radicands and
  complex `EllipticPi` were not touched

**Success criteria** (honest record):

- 1892-problem dual metric: **solved 370 (19.56%), verified 343/370 (92.7%),
  mismatches 0**; per-case diff **21 newly solved, 0 regressed** (bucket deltas:
  special +15, trig +3, power-binomial +2, mixed-other +1)
- Zero wrong answers in the solved set: **`verify_mismatches` = 0 — met** (this is the
  wave's hard red line, and the direct reason `log_fraction` was withdrawn)
- Timeouts 13 → **12** (target ≤ 5 — **not met**), crashes 0 (**met**), wall clock 461.2 s
  (target < 445.7 s — **not met**, +3.5%: the new reduction families' scan in the `special`
  stage and the exact fold's entry traversal are the main sources)
- Verified ratio 92.7% (target ≥ 95% — **not met**: all 27 inconclusive cases are
  "undecidable" rather than wrong, and 4 of the new ones come from the corpus oracle's
  `|m| ≤ 0.9` window)
- Quality gates: fmt / both clippy tiers / deny / workspace test all green
  (`ocas-calc` 416 passed)
- **0.27-line freeze determination**: by the literal `EVOLUTION_PLAN` rule ("two consecutive
  waves adding < 60 solved problems freeze the 0.27 line"), 0.27.2 (+38) and 0.27.3 (+21)
  are **two consecutive waves each below 60** → **the 0.27 line is frozen; 0.28.0 is the
  next active line.** Under that determination 0.27.3 is the last release of the 0.27 series
- Next wave (0.28.0): the mechanism push (correctness foundation → transcendental Risch →
  algebraic extensions → non-elementary layer → complexity), absorbing the 12 timeout cases,
  the 136-case multi-radical half-power cluster, elliptic-family breadth and the two measured
  defects (residue resolution, cycle detection) that this wave quantified; the former 0.28.0
  (Gröbner) slides to 0.33.0 (see the Phase-5 note and
  GENERAL_MECHANISM_FEASIBILITY_EN.md)

### 0.28.0 — Integration-Mechanism Correctness Foundation (defects + certificates + regular towers)

**Goal**: fix the two architectural defects the 0.27.3 investigation located, upgrade
"correct" from sampled evidence to **symbolic certificates**, and relax the general Risch
engine's most expensive entry condition (dependent generators being rejected).
(Basis: GENERAL_MECHANISM_FEASIBILITY_EN.md §3, §5, §6 P0.)

**Deliverables**:

- [ ] **Residue resolution** (defect): the `Integral(...)` residues produced by
  `rational::integral_fallback` and by `symbolic_rational`'s high-degree branch are resolved
  at the **chain tail with a deterministic budget** (prototype measured net +4: +5 solved /
  −1 regressed / +3 timeouts — see BENCHMARK_RESULTS_CN.md §"0.27.3 follow-up" §3.1)
- [ ] **Expression-level cycle detection** (defect) replacing/narrowing the global
  `MAX_CHAIN_ENTRIES` total cap (`rubi-00008` measured 306 stage entries / 18 cycles)
- [ ] **Symbolic certificates**: every output must pass `normalize(D(F) − f) == 0` inside the
  differential field; the numerical oracle is demoted to a probe/regression tool
- [ ] **Three-valued `Outcome`**: `Found { value, certificate }` /
  `ProvedNonElementary { witness }` / `Unknown`; `Unknown` honestly returns `Integral(...)`
- [ ] **Regular towers**: merge algebraically dependent generators (`log x`/`log 2x`,
  `exp x`/`exp(x+1)`, the `exp(±u)` pairs produced by hyperbolic rewriting); return `Unknown`
  when the dependency is undecidable
- [ ] Metric: the harness gains `certified_rate = certificates passed / solved`, gated at
  exactly 1.0 in CI

**Success Criteria**:

- Symbolic certificates are exactly 0 for **100%** of the 1892 solved set; `certified_rate = 1.0`
- New solves in the hyperbolic family (baseline: 111 unsolved cases contain hyperbolic heads)
- `verify_mismatches` stays 0; per-case diff **0 regressions**; the timeout count does not rise

### 0.29.0 — Completing the Transcendental Risch

**Goal**: widen `integral/rde.rs` from "polynomial solutions only" to the complete fragment,
and add the structure theorem for the logarithmic part.
(Basis: GENERAL_MECHANISM_FEASIBILITY_EN.md §3.4, §6 P1.)

**Deliverables**:

- [ ] RDE **rational solutions** (denominator bounds + `D`-rational solutions, the full
  Bronstein ch. 6)
- [ ] **Coupled differential systems** (`D y + A y = b`)
- [ ] The general **structure theorem for the logarithmic part** (no longer only the
  logarithmic-derivative identity)
- [ ] Generality testing: proptest generates random elements **from the tower** → integrate →
  check the symbolic certificate

**Success Criteria**:

- Clear improvement in the `exp-log` bucket (baseline 79/83 unsolved); the certificate gate holds
- 100% certificate pass rate on the random tower-element family; 0 regressions on the 1892

### 0.30.0 — Algebraic Extensions and Inverse-Function Substitution

**Goal**: wire in the algebraic integration chain (reusing the existing Trager assets) and
cover the inverse-function families.
(Basis: GENERAL_MECHANISM_FEASIBILITY_EN.md §3.2, §3.5, §6 P2.)

**Deliverables**:

- [ ] Algebraic generators entering the tower (`√x` and general radicals; `tower/build.rs`
  currently rejects them outright)
- [ ] **Integral basis** (Trager) + **algebraic Hermite reduction** + algebraic residues
  (Bronstein ch. 7–8)
- [ ] Multi-radical-base reduction (an estimated 222 cases with ≥2 distinct radical bases)
- [ ] Inverse-function substitution engine: polynomial weight × `(a+b·f(ax+b))^k`
  (the `inverse-trig-hyper` bucket, baseline 138/146 unsolved)
- [ ] Reuse `ocas-poly`'s Trager factorisation + resultants + algebraic-field GCD

**Success Criteria**:

- The `1/(1+x⁴)`, `1/(1−3x²+x⁴)` and irreducible sextic/octic denominator families solve
- Clear improvement in the radical bucket (baseline 356/407 unsolved); the certificate gate holds

### 0.31.0 — The Non-Elementary Layer

**Goal**: cover what is **completely unreachable** today: antiderivatives that need
dilogarithms or Meijer G (baseline: 122 cases whose reference answer contains `polylog`,
8.0% of the unsolved set).
(Basis: GENERAL_MECHANISM_FEASIBILITY_EN.md §3.6, §6 P3.)

**Deliverables**:

- [ ] `polylog`/`Li₂` function heads + numerical oracle evaluators (add the head and its
  numeric verification **before** the reduction — the pattern 0.27.3 already proved out for
  `Ei`/`Si`/`Ci`/`Shi`/`Chi`/Fresnel)
- [ ] Li₂ reduction (logarithmic-rational integration)
- [ ] Meijer G / Slater expansion (or the holonomic/D-finite route) as the general
  non-elementary mechanism
- [ ] `ProvedNonElementary` becomes available for the first time (an engineered subset of
  Singer's structure theorems)

**Success Criteria**:

- The reachable part of the 122 `polylog` cases unlocks; the rest of the `exp-log` bucket improves
- The new function heads verify numerically in the oracle; the certificate gate holds

### 0.32.0 — Complexity and Performance (treating complexity as an algorithmic problem)

**Goal**: solve coefficient blow-up under symbolic coefficients with **algorithms** instead of
more wall-clock budgets (baseline: 9 of the 12 timeouts sit in `symbolic_rational`).
(Basis: GENERAL_MECHANISM_FEASIBILITY_EN.md §3.7, §6 P4.)

**Deliverables**:

- [ ] Modular arithmetic under symbolic coefficients (reusing the 0.25 MultiModular ℚ pipeline)
- [ ] Parallel Risch
- [ ] Lazy series / degree-bound pruning replacing "add budget" mitigation
- [ ] Stage-level treatment of the 12 remaining timeouts (`symbolic_rational` 9, `binomial` 1,
  `rational` 1, `trig_kernel` 1 — see `timeout_attribution_0273.csv`)

**Success Criteria**:

- Wall clock no worse than the 0.28.0 baseline (0.27.3: 461.2 s); the timeout count falls
- High-level/many-symbol instances no longer fall back through coefficient blow-up; the
  certificate gate holds

### 0.33.0 — Gröbner Performance at Scale (katsura + cyclic-7) (formerly 0.28.0)

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

### 0.34.0 — Code Generation Extension (LLVM/inkwell JIT) (formerly 0.29.0)

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

### 0.35.0 — Matrix Engine + Platform Close-Out + 1.0 Freeze Preparation (formerly 0.30.0)

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

**Target**: after 0.35.0

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

> The LLVM/Inkwell JIT backend moved to 0.34.0, and the quadratic-sieve integer
> factorisation moved to 0.35.0 (both pre-1.0).

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
| 0.28.0 | Beta | Month 53 | Integration-mechanism correctness foundation (residue resolution + expression-level cycle detection + symbolic certificates + three-valued outcome + regular towers) |
| 0.29.0 | Beta | Month 55 | Completing the transcendental Risch (rational RDE solutions + coupled systems + log-part structure theorem) |
| 0.30.0 | Beta | Month 57 | Algebraic extensions and inverse-function substitution (integral basis + algebraic Hermite + residues + multi-radical bases + inverse-function engine) |
| 0.31.0 | Beta | Month 59 | Non-elementary layer (`polylog`/`Li₂` heads + Li₂/Meijer G reductions + `ProvedNonElementary`) |
| 0.32.0 | Beta | Month 61 | Complexity and performance (modular algorithms + parallel Risch + lazy series) |
| 0.33.0 | Beta | Month 63 | Gröbner performance at scale (katsura-6 < 1 s, cyclic-7 within 10× of msolve) (P1; formerly 0.28.0) |
| 0.34.0 | Beta | Month 65 | Code generation extension (LLVM/inkwell JIT backend) (P1; formerly 0.29.0) |
| 0.35.0 | Beta | Month 67 | Matrix engine (DomainMatrix analogue + Smith/Hermite) + Windows FLINT + quadratic sieve + 1.0 freeze preparation (P2/P3; formerly 0.30.0) |
| 1.0.0 | Stable | Month 69 | Stable release (frozen after Phase B++++ competitive gap closure: P0 integration-mechanism correctness and generality met + P1 Gröbner aligned with msolve + LLVM JIT landed + performance ahead of SymPy) |

---

## How to Read This Roadmap

- Each version represents a **potentially publishable** increment.
- Dates are approximate and depend on contributor availability.
- Features may shift between versions based on user feedback and technical discoveries.

---

## Contributing to the Roadmap

If you want to work on a specific version or feature, please open a GitHub issue
and we will assign a tracking issue to you.
