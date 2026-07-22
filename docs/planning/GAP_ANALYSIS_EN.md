# Gap Analysis: oCAS vs Reference Systems

This document tracks the implementation completeness of oCAS milestone by
milestone (0.1 → 1.0+) and the gap against the three reference systems:
**Symbolica** (Rust), **SageMath** (Python ecosystem), and **SymPy** (pure
Python). It is a living document and must be refreshed at every version bump.
For the Chinese edition, see [GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md).

> Last evaluated: **0.18.1 @ 2026-07-23** (Python/C bindings backfill for the three 0.18.0 capabilities + prelude completeness)

---

## Legend

| Mark | Meaning |
|---|---|
| ✅ | Complete |
| 🟡 | Basic / partial |
| 🔴 | Missing or major gap |
| ⚠️ | Complete with caveats |

---

## 1. Version Completion Status (0.1–0.15.2)

| Version | Phase | Roadmap | Verified Status |
|---|---|---|---|
| 0.1.0 | Pre-Alpha | ✅ | ✅ 12-crate workspace, CI, `OcasError`, arena (Miri-aware), rayon pool, FFI glue, `gmp` feature via `rug` |
| 0.2.0 | Pre-Alpha | ✅ | ✅ `ocas-atom`, `Atom` tagged union, arena AST, hash consing, logos lexer, Pratt parser, normalizer |
| 0.3.0 | Alpha | ✅ | ⚠️ `Integer/Rational/FiniteField/RealBall/Complex`; dense/sparse poly, div-rem; `flint` Linux/WSL only, GMP/MPFR via `rug` |
| 0.4.0 | Alpha | ✅ | ✅ matcher, pattern, rules, simplify, transformer, `egraph.rs` (egg integration) |
| 0.5.0 | Alpha | ✅ | ⚠️ derivative, integral (heuristic), Taylor series; integration is table-based, no Risch |
| 0.6.0 | Alpha | ✅ | ✅ stable `ocas` prelude, rustdoc examples, proptest, criterion, SymPy harness, crates.io prep |
| 0.7.0 | Beta | ✅ | ⚠️ linear (rational/integer) + Diophantine + polynomial system (Gröbner); Sturm root isolation; assumptions — algorithms are basic |
| 0.8.0 | Beta | ✅ | ✅ tree interpreter, AST→instruction compiler, function registry, Cranelift JIT, SIMD vectorized eval |
| 0.9.0 | Beta | ✅ | ⚠️ PyO3 `Expression`/`Evaluator`/`solve_*`; cbindgen + C++ RAII wrapper — some classes deferred to 0.10 |
| 0.10.0 | Beta | ✅ | ✅ Python `Polynomial/Matrix/Domain`, Matrix linear algebra (Bareiss), mdBook docs site, 3-platform wheels CI, version frozen at 0.10.0 |
| 0.11.0 | Beta | ✅ | ✅ Complete polynomial factorization over ℤ and ℤ_p (Yun SFF → CZ → Hensel → Zassenhaus), multivariate GCD, 500 proptest round-trip cases, version bumped to 0.11.0 |
| 0.11.1 | Beta | ✅ | ✅ Bivariate factorization over ℤ and ℤ_p (monic-in-x Wang Hensel), sparse multivariate `factor()` entry points, C polynomial bindings, mdBook factorization chapter, version bumped to 0.11.1 |
| 0.12.0 | Beta | ✅ | ✅ Rational polynomial `RationalPolynomial<D,O>`, Brown PRS resultant, Karatsuba multiplication, extended GCD, polynomial CRT/Diophantine, p-adic expansion, partial fraction decomposition, rational reconstruction, version bumped to 0.12.0 |
| 0.12.1 | Beta | ✅ | ✅ Self-implemented NTT over ℤ_p, `pulp` SIMD dispatch, Estrin polynomial evaluation, sparse matrix backend for F4, numerical verification features, version bumped to 0.12.1 |
| 0.13.0 | Beta | ✅ | ✅ F4 Gröbner basis algorithm with Gebauer-Moeller pair filtering and simplification cache, `Grlex` ordering, `Domain` trait extensions, `FiniteField` ℤ_p fast-path utilities, version bumped to 0.13.0 |
| 0.14.0 | 1.0 Candidate | ✅ | ✅ Risch symbolic integration (Hermite reduction, log-derivative identity, primitive undetermined coefficients, hyperexponential RDE), rational-function integration (Hermite + Rothstein–Trager), special-function table (erf/Ei/Si/Ci/Fresnel), trigonometric integration (exp(I·x) + realify), FGLM/F5/Hilbert, `reorder`, two mdBook chapters |
| 0.15.0 | 1.0 Candidate | ✅ | ✅ Multi-output JIT (97×/21×), f32 mixed precision (JIT + SIMD 16 lanes), streaming evaluation (constant memory over 1M rows), const-folding + stack compaction, Arena reset + workspace pool, ahash hot-path replacement, native i64 F4 pipeline; cyclic-6 <5s deferred to 0.15.1 (needs RREF/F5) |
| 0.15.1 | 1.0 Candidate | ✅ | ✅ Real F4 linear algebra fix: descending matrix column order + echelon write-back condition + Symbolica GM criteria port + classic extraction (separate multiples + input-heads, zero reduction). cyclic-5 ℤ₁₃ 2609 s → 31 ms (~85,000×) with first-ever `is_groebner_basis` pass; cyclic-6 tractable (9970 s); <5s deferred to 0.15.2 (LM index + sparse echelon) |
| 0.15.2 | 1.0 Candidate | ✅ | ✅ Reducer LM hash index (support-mask buckets + submask enumeration) + sparse-row echelon (two-pointer merge cancellation, O(nnz)/op) + hashed extraction dedup + worklist preprocessing + row-template cache. cyclic-6 ℤ₁₃ 9970 s → 3670 s (2.7×, basis=20 correct); phase profile shifted to elimination-dominated (echelon ≈89%); <5s not reached (264k rows is F4's intrinsic size, needs F5 signature reduction) |

All 0.1–0.15.2 deliverables landed. The workspace is pinned at 0.15.2. Quality
gates are green: `cargo fmt`, `clippy -D warnings`, workspace tests,
`cargo deny`, pytest cases, `mdbook build`.

---

## 2. Code Scale

Snapshot of `src/` Rust lines (excluding tests and benches).

| Crate | Files | Lines |
|---|---|---|
| ocas-poly | 22 | ~10,560 |
| ocas-calc | 18 | ~5,649 |
| ocas-eval | 13 | ~3,855 |
| ocas-domain | 10 | ~3,337 |
| ocas-rewrite | 7 | ~1,593 |
| ocas-py | 7 | ~1,461 |
| ocas-c | 4 | ~1,454 |
| ocas-core | 5 | ~1,115 |
| ocas-atom | 4 | ~1,111 |
| ocas-parse | 3 | ~495 |
| ocas (prelude) | 1 | ~115 |
| ocas-gpl | 1 | 0 (placeholder) |
| **Total src** | **95** | **~30.7k** |

Up ~70% from the 0.10 snapshot (66 files / ~18k lines); growth comes mainly
from Risch and rational-function integration (ocas-calc), F4/FGLM/F5 and
factorization (ocas-poly), and multi-output JIT / streaming evaluation
(ocas-eval).

`ocas-gpl` is a placeholder; GPL-exclusive backends are Post-1.0 work, in line
with the roadmap.

---

## 3. Algorithm Depth Audit

This section is the single most decisive factor in CAS maturity and the main
source of the gap.

| Algorithm Area | oCAS Status | Maturity |
|---|---|---|
| Polynomial factorization | `factor()` on `DenseUnivariatePolynomial` over ℤ and ℤ_p, arbitrary multivariate `factor()` on `SparseMultivariatePolynomial` over ℤ and ℤ_p (0.16.x Wang EEZ + non-constant LC imposition), plus univariate `factor()` over `AlgebraicNumberField` (0.17.0 Trager: shifted norm + modular GCD) | 🟢 Univariate/bivariate/multivariate/ANF (univariate) |
| Gröbner basis | F4 with real linear algebra (0.15.1: descending column order + Symbolica GM criteria + classic extraction) + FGLM + experimental F5 + native i64 ℤ_p pipeline; cyclic-5 ℤ₁₃ 23 ms (re-measured 2026-07-21) | 🟢 F4 complete |
| Symbolic integration | Risch (elementary transcendental towers + RDE polynomial fragment) + rational-function Hermite + trig exp(I·x) + special-function table (erf/Ei/Si/Ci/Fresnel); falls back to `Integral(...)` | 🟢 Risch done |
| Real root isolation | Sturm sequence + interval isolation + refine (univariate); known gap: only 8/10 roots isolated on expanded Wilkinson n=10 | 🟡 Fairly complete |
| Polynomial GCD | GCD + primitive part + extended GCD (0.12); no modular GCD for large integer coefficients | 🟡 Usable |
| Linear solving | Rational/integer linear systems + bivariate Diophantine (`ax+by=c`) | 🟡 Usable, limited scale |
| JIT evaluation | Cranelift backend; ≥10x speedup target met (per roadmap criterion) | 🟢 Complete |

---

## 4. Gap Analysis vs Reference Systems

### 4.1 vs Symbolica (Rust, AGPL)

Symbolica's `examples/` directory reveals the maturity gap. After 0.11–0.15,
oCAS covers most of Symbolica's core feature surface; the gap has narrowed to
breadth and large-scale performance.

| Capability | oCAS | Symbolica |
|---|---|---|
| Polynomial factorization | � univariate ℤ and ℤ_p (CZ + Hensel + Zassenhaus) + arbitrary multivariate (0.16.0 Wang EEZ, constant-LC preprocessing); algebraic number fields missing; non-constant LC imposition 0.16.1 | ✅ full (arbitrary multivariate + algebraic number fields, `factorization.rs`) |
| Rational polynomials | ✅ `RationalPolynomial<D,O>` with GCD canonicalization | ✅ `rational_polynomial.rs` |
| Partial fractions | ✅ `apart()` / `together()` over Euclidean domains | ✅ `partial_fraction.rs` |
| Rational reconstruction | ✅ `rational_reconstruction(a, m)` via extended Euclidean | ✅ `rational_reconstruction.rs` |
| Numerical integration | 🔴 none | ✅ `numerical_integration.rs` |
| Streaming API | ✅ `streaming.rs` (`StreamingEvaluator`: chunked input + reused stack, constant memory over 1M rows) | ✅ `streaming.rs` |
| Tensors / dual numbers | 🔴 none | ✅ `tensors.rs` / `dual.rs` |
| Optimization / codegen | ✅ multi-output JIT (`compile_multi` + CSE + const folding + stack compaction) + f32 mixed precision | ✅ `optimize.rs` / multi-output |
| Gröbner bases | 🟡 F4 complete + large-scale perf (0.15.2: LM index + sparse echelon + row-template cache, cyclic-6 ℤ₁₃ 9970 s → 3670 s); cyclic-6 <5s not reached, needs F5 signature reduction | ✅ industrial-grade |
| Resource control (fuel) | 🔴 none | ✅ `fuel_backend.rs` |

Symbolica 2.1.0's core strengths — industrial factorization, rational function
arithmetic, multi-output optimization, streaming — have largely been closed by
oCAS during 0.11–0.15. The remaining gaps are: arbitrary multivariate (≥3
variables) and algebraic-number-field factorization, numerical integration,
tensors / dual numbers, fuel-based resource control, and Gröbner performance
at scale (cyclic-6 class, where Symbolica is still far ahead).

### 4.2 vs SageMath (Python ecosystem)

SageMath is a "Swiss-army-knife" scientific environment. The gap is
**breadth-level**.

| Domain | oCAS | SageMath |
|---|---|---|
| Algebraic geometry | 🟡 basic Gröbner | ✅ Singular integration |
| Number theory | 🟡 basic Diophantine | ✅ PARI/FLINT full stack |
| Differential equations | 🔴 none | ✅ full ODE/PDE solvers |
| Group / representation theory | 🔴 none | ✅ GAP integration |
| Combinatorics | 🔴 none | ✅ complete |
| Plotting / visualization | 🔴 none | ✅ matplotlib integration |
| Database interfaces | 🔴 none | ✅ OEIS / LMFDB |

SageMath achieves breadth by wrapping 80+ specialized libraries; oCAS is a
self-contained kernel. The two have different positioning — oCAS targets a
high-performance **library**, SageMath a full **environment**. Comparison is
meaningful mainly on core algebra performance, not feature breadth.

### 4.3 vs SymPy (pure Python)

SymPy is oCAS's most realistic target for both feature parity and performance
leadership.

| Domain | oCAS vs SymPy | Note |
|---|---|---|
| Parsing / simplification | 🟢 parity | both complete |
| Differentiation | 🟢 parity | chain/product/power rules |
| Integration | � rough parity | both have Risch (oCAS since 0.14); SymPy's heuristic/manual fallbacks are broader, oCAS returns `Integral(...)` when uncovered |
| Factorization | 🟡 oCAS slightly weaker | univariate ℤ and ℤ_p parity (CZ + Hensel + Zassenhaus); oCAS supports bivariate, SymPy arbitrary multivariate |
| Gröbner | 🟢 oCAS advantage | oCAS F4 matrix algorithm with real linear algebra (cyclic-5 ℤ₁₃ 23 ms) outperforms SymPy's Buchberger implementation |
| Matrix / linear algebra | 🟢 parity | oCAS has Bareiss determinant/inverse |
| **Performance** | 🟢 **oCAS advantage** | Rust + Cranelift JIT + arena vs pure Python; measured x³⁰−1 square-free factorization 39 µs vs SymPy full factor ~0.9 ms (~24×, 2026-07-21) |
| Python ergonomics | 🟢 parity | oCAS has `ocas-py` bindings |

The 0.6.0 success criterion — "parity with SymPy on basic polynomial,
calculus, and rewriting" — is met and exceeded on the **performance** axis,
and **integration** was closed by Risch in 0.14; the remaining feature gaps
against SymPy are **arbitrary multivariate factorization** and the **breadth
of integration heuristic fallbacks**.

---

## 5. Key Gaps & Priorities

Ranked by impact × implementation cost. All hard-algorithm gaps planned before
1.0 are closed; the remaining Symbolica gaps are scheduled into Phase B+
(0.15.2–0.18.0, see EVOLUTION_PLAN) with the goal of closing them all before
1.0.

| # | Gap | Priority |
|---|---|---|
| 1 | ~~Full polynomial factorization~~ (completed 0.11.0–0.11.1) | ✅ done — univariate and bivariate (monic-in-x) closed; ≥3 variables see #7 |
| 2 | ~~Risch symbolic integration~~ (completed 0.14) | ✅ done — elementary transcendental towers + RDE fragment + rational Hermite + special-function table |
| 3 | ~~Gröbner F4/F5~~ (completed 0.13 / 0.14 / 0.15.1) | ✅ F4 with real linear algebra + FGLM + experimental F5; scale performance see #6 |
| 4 | ~~Rational polynomials / partial fractions~~ (completed 0.12) | ✅ done — `RationalPolynomial` type + partial fractions + resultant + Karatsuba multiplication |
| 5 | ~~Multi-output optimization / codegen~~ (done in 0.15) | ✅ done — multi-output JIT (97×/21×) + f32 mixed precision + CSE/const-folding/stack-compaction |
| 6 | Gröbner performance at scale (cyclic-6 ℤ_p < 5 s) | � 0.15.2 done (9970 s → 3670 s, 2.7×); <5 s needs F5 signature reduction (post-1.0) |
| 7 | ~~Arbitrary multivariate (≥3 variables) factorization~~ (completed 0.16) | ✅ done — Wang EEZ lifting + LC preprocessing (constant LC) + Zassenhaus; non-constant LC imposition see #7a |
| 7a | ~~Non-constant leading-coefficient imposition + multivariate sparsity~~ (completed 0.16.1/0.16.2) | ✅ done — mod-p Hensel imposition + sparse Diophantine + field Wang preprocessing on the Fp path |
| 8 | ~~Algebraic-number-field factorization~~ (completed 0.17) | ✅ done — Trager algorithm (shifted norm + ℚ factorization + GF(p^d) modular GCD), univariate path; multivariate extension deferred |
| 9 | ~~Numerical integration / dual numbers / tensor basics / fuel~~ (done in 0.18) | ✅ Done — Vegas + HyperDual + index contraction + fuel; 0.18.1 backfilled the Python/C bindings |
| 10 | ODE/PDE solvers (Post-1.0) | 🟢 high user demand |

---

## 6. Overall Assessment

Execution quality of 0.1 → 0.15.1 is high: every roadmap deliverable shipped,
the layered architecture is clean (no cycles), the 12-crate workspace is
strictly layered, quality gates are strict (`-D warnings` + deny + Miri
awareness), and docs/bindings/CI are well-engineered. The three hard
algorithms planned before 1.0 — polynomial factorization (0.11), Gröbner F4
(0.13, real linear algebra fixed in 0.15.1), and Risch symbolic integration
(0.14) — are all closed and continuously regressed via the SymPy/Symbolica
cross-verification framework.

Realistic positioning: oCAS today is "a high-performance SymPy core, with
Risch symbolic integration, univariate/bivariate factorization and rational
functions, Gröbner F4 with real linear algebra, and multi-output JIT /
streaming evaluation". Re-measured 0.15.1 performance: F4 cyclic-5 ℤ₁₃
23 ms; x³⁰−1 square-free factorization 39 µs (SymPy full factor ~0.9 ms,
~24×); JIT 97× single-output, 21× three-output.

The remaining pre-1.0 work is Phase B+ "Closing the Symbolica Gap" (0.15.2
Gröbner performance at scale → 0.16 arbitrary multivariate factorization ✅ →
0.16.1 non-constant leading-coefficient imposition ✅ → 0.17
algebraic-number-field factorization ✅ → 0.18 numerical integration /
duals / tensors / fuel, see EVOLUTION_PLAN), after which 1.0.0 is
stabilization and release engineering only (API freeze, coverage, migration
guides, signed artifacts). ODE/PDE and full tensor calculus remain Post-1.0
topics.

---

## Update Log

Record every refresh here (version, date, evaluator, deltas).

| Version | Date | Deltas |
|---|---|---|
| 0.10.0 | 2026-07-02 | Initial assessment. All 0.1–0.10 deliverables verified complete; gaps against Symbolica / SageMath / SymPy documented; factorization + Risch integration identified as top priorities. |
| 0.11.0 | 2026-07-03 | Polynomial factorization completed (univariate ℤ and ℤ_p); multivariate GCD added; SymPy comparison updated to parity for univariate factorization; highest-priority gap shifted to rational functions / partial fractions (0.12). |
| 0.11.1 | 2026-07-04 | Bivariate factorization over ℤ and ℤ_p (monic-in-x Wang Hensel) added; sparse multivariate `factor()` entry points and C polynomial bindings landed; mdBook factorization chapter added; highest-priority gap remains rational functions / partial fractions (0.12). |
| 0.12.0 | 2026-07-04 | Rational function stack completed (`RationalPolynomial` + partial fractions + Brown PRS resultant + Karatsuba multiplication + rational reconstruction); parity with Symbolica for rational functions; highest-priority gap shifted to Gröbner F4 (0.13) and Risch integration (0.14). |
| 0.13.0 | 2026-07-06 | Gröbner F4 matrix algorithm completed (Faugère 1999); Gebauer-Moeller pair filtering + simplification cache + ℤ_p fast path; `minimize()` bug fix; Gröbner upgraded from 🟡 to 🟢; highest-priority gap shifted to Risch integration (0.14). |
| 0.13.1 | 2026-07-17 | Patch release: docs.rs builds now use portable features only (no gmp/mpfr/flint/python/gpl), restoring hosted documentation; no algorithm changes, gap conclusions unchanged from 0.13.0. |
| 0.13.2 | 2026-07-18 | Engineering & distribution milestone: `pip install ocas` live on PyPI (5 platform wheels + sdist, incl. both macOS archs); OIDC trusted publishing pipeline established; crossbeam-epoch RUSTSEC-2026-0204 fixed; cranelift/chumsky/logos/cbindgen/criterion/hashbrown/flint3-sys/egg upgraded; no algorithm changes, gap conclusions unchanged. |
| 0.14.0 | 2026-07-18 | Risch symbolic integration completed (elementary transcendental towers + RDE polynomial fragment); rational-function integration (Hermite + logarithmic part); special-function table (erf/Ei/Si/Ci/Fresnel) closing the 0.11.0 known gap `exp(-x²)→erf`; trigonometric exp(I·x) + realify; Gröbner wrap-up (FGLM zero-dimensional conversion + experimental F5 + Hilbert bounds + reorder); parser `-x^2` precedence fix; symbolic integration upgraded from 🟡 to 🟢; highest-priority gap shifted to 0.15 performance / multi-output JIT. |
| 0.15.0 | 2026-07-20 | Multi-output JIT (97×/21×) + f32 mixed precision + streaming evaluation (constant memory over 1M rows) + const-folding/stack-compaction + Arena reset/workspace pool + ahash + native i64 F4 pipeline; JIT calling-convention Windows fix; F4 bottleneck localized via section timing (extraction = 99.98%); cyclic-6 <5s deferred to 0.15.1 (needs RREF/F5); highest-priority gap shifted to 1.0 stable release. |
| 0.15.1 | 2026-07-20 | Real F4 linear algebra fix: descending matrix column order (was ascending — echelon was decorative, F4 was effectively Buchberger) + echelon write-back condition + Symbolica GM criteria port + classic extraction (separate multiples + input-heads, zero reduction at extraction); cyclic-5 ℤ₁₃ 2609 s → 31 ms (~85,000×) with first-ever `is_groebner_basis` pass; cyclic-6 tractable (9970 s, basis=20); <5s deferred to 0.15.2 (LM index + sparse echelon). |
| 0.16.0–0.16.2 | 2026-07-21 | Arbitrary multivariate factorization stack (Wang EEZ + Hensel + non-constant leading-coefficient imposition + sparse Diophantine small-prime escalation), covering both ℤ and 𝔽ₚ multivariate paths; multivariate factorization upgraded from 🔴 to 🟢. |
| 0.17.0 | 2026-07-22 | Algebraic number field factorization (Trager) completed: `AlgebraicNumberField` + modular GCD over number fields (GF(p^d) + CRT + rational reconstruction) + shifted norm; fixed general-degree Brown PRS resultant bug; algebraic-number-field factorization upgraded from 🔴 to 🟢 (univariate path). |
| 0.17.1 | 2026-07-22 | Patch: algebraic-number-field Python/C bindings completion (`AlgebraicExtension`/`AlgebraicElement`/`AlgebraicPolynomial` Python classes + `OcasAlgebraicField`/`OcasAlgebraicPoly` opaque handles and `ocas_algebraic_*` C ABI + `RootOf` parse confirmation); no algorithm changes, gap conclusions unchanged. |
| 0.18.0 | 2026-07-23 | Numerical integration (Vegas adaptive Monte Carlo + `integrate_1d`), forward automatic differentiation (`HyperDual<T>` runtime shape), fuel resource control (`Fuel` + `simplify_with_fuel`/`integrate_with_fuel`), and tensor basics (independent `Tensor` type + explicit contraction + symmetrisation sign) landed; added `rand`/`rand_xoshiro` dependencies; full tensor canonicalisation and the deterministic quadrature bridge deferred. |
| 0.18.1 | 2026-07-23 | Patch: backfilled Python/C bindings for the three 0.18.0 capabilities (numerical integration / dual AD / tensor basics) — `ocas-py::{numeric,tensor,dual}` modules + `ocas-c::{numeric,tensor,dual}` opaque handles and C ABI + `include/ocas.h` synced + prelude re-exports for tensor / dual / `StatisticsAccumulator`; 41 Python tests + 31 C API tests; no algorithmic change, gap conclusions unchanged. |
| 0.15.1 | 2026-07-21 | Re-evaluation: code-scale snapshot updated to 95 files / ~30.7k lines (+~70% vs 0.10's ~18k); F4 cyclic-5 ℤ₁₃ re-measured at 23 ms; new measurement x³⁰−1 square-free factorization 39 µs vs SymPy full factor ~0.9 ms (~24×); stale post-0.14/0.15 statements fixed (§3 GCD/root-isolation, §4.1 "largely absent" paragraph, §4.3 integration/factorization/Gröbner, §5 Risch priority, mojibake characters); gaps re-ranked — all pre-1.0 hard algorithms closed, remaining items moved to Post-1.0: arbitrary multivariate (≥3 variables) + algebraic-number-field factorization, numerical integration, tensors / dual numbers, ODE/PDE; cyclic-6 <5s scoped to 0.15.2. |
| 0.15.2 | 2026-07-21 | Gröbner performance at scale: reducer LM hash index (support-mask buckets + submask enumeration, removing the O(monomials × basis) linear scan) + sparse-row echelon (two-pointer merge cancellation, O(nnz)/op, replacing the dense buffer) + hashed extraction dedup + worklist preprocessing + row-template cache; cyclic-6 ℤ₁₃ 9970 s → 3670 s (2.7×, basis=20 correct), phase profile shifted to elimination-dominated (echelon ≈89%); <5s not reached — the cyclic-6 F4 matrix hits 264k rows × 284k cols at round 22 (intrinsic to F4), a further order-of-magnitude win needs F5 signature reduction (eliminating zero-reducing rows), moved to post-1.0; version bumped to 0.15.2. |
| 0.16.0 | 2026-07-21 | Arbitrary multivariate factorization (Wang EEZ) done: landed `factor::eez` (generic multivariate Diophantine + per-variable EEZ Hensel lifting + $n$-variate GCD + characteristic-$p$ $p$-th powers + Wang LC preprocessing [constant LC] + Zassenhaus recombination); `factor()` generalized to any arity; three pre-existing bugs fixed (`div_rem_sparse` divisibility order, Diophantine loop bound, non-monic univariate factorization); factorization upgraded 🟡 → 🟢 (univariate/bivariate/arbitrary multivariate); 0.16.1 added (non-constant LC imposition + sparsity); version bumped to 0.16.0. |
| 0.17.0 | 2026-07-22 | Algebraic-number-field factorization (Trager) done: new `ocas-domain::algebraic` (`AlgebraicExtension<D>` — one implementation for ℚ(α) and GF(p^d), EEA inversion) + `ocas-poly::factor::algebraic` (shifted norm via evaluation–interpolation resultants + modular number-field GCD [GF(p^d) + CRT + rational reconstruction + trial division] + rational fast path); fixed the Brown PRS resultant bug for general degrees (beta division was applied only for unit betas — not a valid resultant algorithm; re-ported from Symbolica's `resultant_prs`); 0.16.2 sparse Diophantine small-prime escalation completed; factorization now covers univariate/bivariate/multivariate/ANF (univariate); performance target met (degree ≤ 12 at 8–32 ms < 100 ms); version bumped to 0.17.0. |
