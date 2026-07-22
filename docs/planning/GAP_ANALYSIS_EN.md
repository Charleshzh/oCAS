# Gap Analysis: oCAS vs Reference Systems

This document tracks the implementation completeness of oCAS milestone by
milestone (0.1 → 1.0+) and the gap against the three reference systems:
**Symbolica** (Rust), **SageMath** (Python ecosystem), and **SymPy** (pure
Python). It is a living document and must be refreshed at every version bump.
For the Chinese edition, see [GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md).

> Last evaluated: **0.18.1 @ 2026-07-23** (full re-evaluation: code-scale snapshot refreshed to 0.18.1, Symbolica gap table re-checked after 0.16–0.18 landed)

---

## Legend

| Mark | Meaning |
|---|---|
| ✅ | Complete |
| 🟡 | Basic / partial |
| 🔴 | Missing or major gap |
| ⚠️ | Complete with caveats |

---

## 1. Version Completion Status (0.1–0.18.1)

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
| 0.16.0 | 1.0 Candidate | ✅ | ✅ Arbitrary multivariate factorization (Wang EEZ): generic multivariate Diophantine + per-variable EEZ Hensel lifting + n-variate GCD + characteristic-p p-th powers + Wang LC preprocessing (constant LC) + Zassenhaus recombination; `factor()` generalized to any arity; three pre-existing bugs fixed (`div_rem_sparse` divisibility order, Diophantine loop bound, non-monic univariate factorization) |
| 0.16.1 | 1.0 Candidate | ✅ | ✅ Non-constant leading-coefficient imposition (mod-p Hensel) + multivariate sparsity improvements on the ℤ path |
| 0.16.2 | 1.0 Candidate | ✅ | ✅ 𝔽_p non-constant LC factorization (Fp Wang LC reconstruction + small-prime escalation for sparse Diophantine) on the Fp path |
| 0.17.0 | 1.0 Candidate | ✅ | ✅ Algebraic-number-field factorization (Trager): `AlgebraicExtension<D>` (one implementation for ℚ(α) and GF(p^d)) + shifted norm via evaluation–interpolation resultants + modular number-field GCD (GF(p^d) + CRT + rational reconstruction); Brown PRS resultant general-degree bug re-ported from Symbolica; degree ≤ 12 ANF factorization at 8–32 ms |
| 0.17.1 | 1.0 Candidate | ✅ | ✅ Algebraic-number Python/C bindings: `AlgebraicExtension`/`AlgebraicElement`/`AlgebraicPolynomial` Python classes + `OcasAlgebraicField`/`OcasAlgebraicPoly` opaque handles and `ocas_algebraic_*` C ABI; `RootOf(poly, idx)` parse confirmation |
| 0.18.0 | 1.0 Candidate | ✅ | ✅ Numerical integration (Vegas adaptive Monte Carlo + `integrate_1d` + `StatisticsAccumulator` + `Integrator` trait), forward automatic differentiation (`HyperDual<T>` runtime shape + truncated product table + geometric-series inverse + `DualCoeff` trait, Rational dual-path), fuel resource control (`Fuel = Arc<AtomicUsize>` + `OutOfFuel` + `simplify_with_fuel`/`integrate_with_fuel`), tensor basics (independent `Tensor` type + index slots + explicit contraction + `symmetrise_sign`); added `rand`/`rand_xoshiro` |
| 0.18.1 | 1.0 Candidate | ✅ | ✅ Patch: Python/C bindings backfill for the three 0.18.0 capabilities (`ocas-py::{numeric,tensor,dual}` + `ocas-c::{numeric,tensor,dual}` opaque handles and C ABI + `include/ocas.h` synced) + prelude re-exports for tensor / dual / `StatisticsAccumulator`; 41 Python tests + 31 C API tests added; `normalize` idempotency bug fixed (drop Num(0)/Num(1) after `merge_numbers` in Add/Mul) |

All 0.1–0.18.1 deliverables landed. The workspace is pinned at 0.18.1. Quality
gates are green: `cargo fmt`, `clippy -D warnings`, workspace tests,
`cargo deny`, pytest cases, `mdbook build`.

---

## 2. Code Scale

Snapshot of `src/` Rust lines (non-blank, excluding tests and benches).

| Crate | Files | Lines |
|---|---|---|
| ocas-poly | 24 | ~15,587 |
| ocas-calc | 18 | ~5,672 |
| ocas-domain | 12 | ~4,475 |
| ocas-eval | 16 | ~4,379 |
| ocas-c | 8 | ~3,195 |
| ocas-py | 11 | ~2,570 |
| ocas-rewrite | 7 | ~1,653 |
| ocas-atom | 5 | ~1,558 |
| ocas-core | 6 | ~1,269 |
| ocas-parse | 3 | ~495 |
| ocas (prelude) | 1 | ~125 |
| ocas-gpl | 1 | 0 (placeholder) |
| **Total src** | **112** | **~40.9k** |

Up ~33% from the 0.15.1 snapshot (95 files / ~30.7k lines) and ~127% from the
0.10 snapshot (66 files / ~18k lines). Growth 0.15.1 → 0.18.1 comes mainly
from arbitrary multivariate + algebraic-number-field factorization
(ocas-poly, +~5.0k), the domain layer (`algebraic` + `dual`, ocas-domain
+~1.1k), numerical integration / streaming (ocas-eval, +~0.5k), and the
Python/C binding expansion for the three 0.18.0 capabilities (ocas-c +~1.7k,
ocas-py +~1.1k).

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
| Polynomial GCD | GCD + primitive part + extended GCD (0.12) + arbitrary-arity multivariate GCD via EEZ (0.16) + modular number-field GCD over GF(p^d) with CRT + rational reconstruction (0.17); no modular-GCD fast path for very large integer coefficients | 🟢 Usable, no HEVMGCD |
| Linear solving | Rational/integer linear systems + bivariate Diophantine (`ax+by=c`) | 🟡 Usable, limited scale |
| JIT evaluation | Cranelift backend; ≥10x speedup target met (per roadmap criterion) | 🟢 Complete |

---

## 4. Gap Analysis vs Reference Systems

### 4.1 vs Symbolica (Rust, AGPL)

Symbolica's `examples/` directory (30 examples) is the maturity benchmark.
After 0.11–0.18, oCAS covers **all** of Symbolica's core example surface;
the gap has narrowed to **scale performance** (cyclic-6 Gröbner) and a few
**specialized pattern transformers** (e.g. `Transformer::Partition` for
argument-sequence partitioning).

| Capability | oCAS | Symbolica |
|---|---|---|
| Polynomial factorization | ✅ univariate ℤ/ℤ_p (CZ + Hensel + Zassenhaus) + arbitrary multivariate (0.16 Wang EEZ + non-constant LC imposition 0.16.1/0.16.2) + algebraic-number-field (0.17 Trager, univariate) | ✅ full (arbitrary multivariate + algebraic number fields, `factorization.rs`) |
| Rational polynomials | ✅ `RationalPolynomial<D,O>` with GCD canonicalization | ✅ `rational_polynomial.rs` |
| Partial fractions | ✅ `apart()` / `together()` over Euclidean domains | ✅ `partial_fraction.rs` |
| Rational reconstruction | ✅ `rational_reconstruction(a, m)` via extended Euclidean | ✅ `rational_reconstruction.rs` |
| Numerical integration | ✅ Vegas adaptive Monte Carlo + `integrate_1d` + `StatisticsAccumulator` (0.18) | ✅ `numerical_integration.rs` |
| Streaming API | ✅ `StreamingEvaluator`: chunked input + reused stack, constant memory over 1M rows | ✅ `streaming.rs` |
| Tensors / dual numbers | ✅ independent `Tensor` type + index contraction + `symmetrise_sign` (0.18 basics; full canonicalisation Post-1.0); `HyperDual<T>` forward AD (0.18) | ✅ `tensors.rs` / `dual.rs` (full graphica-based canonicalisation) |
| Optimization / codegen | ✅ multi-output JIT (`compile_multi` + CSE + const folding + stack compaction) + f32 mixed precision | ✅ `optimize.rs` / multi-output |
| Gröbner bases | 🟡 F4 with real linear algebra complete (0.15.1) + LM index + sparse echelon (0.15.2); cyclic-6 ℤ₁₃ 9970 s → 3670 s (2.7×); <5s not reached, needs F5 signature reduction | ✅ industrial-grade |
| Resource control (fuel) | ✅ `Fuel = Arc<AtomicUsize>` + `simplify_with_fuel`/`integrate_with_fuel` (0.18) | ✅ `fuel_backend.rs` |
| Pattern transformers | 🟡 matcher/replace/transformer complete; specialized sequence `Transformer::Partition` not implemented | ✅ full transformer set |

Symbolica 2.1.0's core strengths — industrial factorization (including
algebraic number fields), rational-function arithmetic, multi-output
optimization, streaming, numerical integration, dual numbers, tensors, and
fuel-based resource control — have **all been closed** by oCAS during
0.11–0.18. The remaining gaps are: **Gröbner-basis performance at scale**
(cyclic-6 class, where Symbolica's F5/signature machinery is still ahead),
**full tensor canonicalisation** (oCAS ships basics only; Symbolica uses the
graphica graph-isomorphism engine), and a few **specialized pattern
transformers** (e.g. `Transformer::Partition`).

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
| Integration | 🟢 rough parity | both have Risch (oCAS since 0.14); SymPy's heuristic/manual fallbacks are broader, oCAS returns `Integral(...)` when uncovered |
| Factorization | 🟢 parity | univariate ℤ/ℤ_p + arbitrary multivariate (0.16 Wang EEZ) + algebraic-number-field (0.17 Trager); SymPy has broader ANF coverage |
| Gröbner | 🟢 oCAS advantage | oCAS F4 matrix algorithm with real linear algebra (cyclic-5 ℤ₁₃ 23 ms) outperforms SymPy's Buchberger implementation |
| Matrix / linear algebra | 🟢 parity | oCAS has Bareiss determinant/inverse |
| **Performance** | 🟢 **oCAS advantage** | Rust + Cranelift JIT + arena vs pure Python; measured x³⁰−1 square-free factorization 39 µs vs SymPy full factor ~0.9 ms (~24×, 2026-07-21) |
| Python ergonomics | 🟢 parity | oCAS has `ocas-py` bindings |

The 0.6.0 success criterion — "parity with SymPy on basic polynomial,
calculus, and rewriting" — is met and exceeded on the **performance** axis,
and **integration** was closed by Risch in 0.14 while **factorization**
reached arbitrary-multivariate parity in 0.16 (plus algebraic-number-field
via Trager in 0.17). The remaining feature gap against SymPy is the **breadth
of integration heuristic fallbacks** (SymPy's `manualintegrate` /
heuristic pool is wider than oCAS's Risch + table path).

---

## 5. Key Gaps & Priorities

Ranked by impact × implementation cost. All hard-algorithm gaps planned before
1.0 are **closed**; Phase B+ "Closing the Symbolica Gap" (0.15.2–0.18.0) is
complete — see EVOLUTION_PLAN. The remaining items are scale performance,
breadth, and Post-1.0 topics.

| # | Gap | Priority |
|---|---|---|
| 1 | ~~Full polynomial factorization~~ (completed 0.11.0–0.11.1) | ✅ done — univariate and bivariate (monic-in-x) closed; ≥3 variables see #7 |
| 2 | ~~Risch symbolic integration~~ (completed 0.14) | ✅ done — elementary transcendental towers + RDE polynomial fragment + rational Hermite + special-function table |
| 3 | ~~Gröbner F4/F5~~ (completed 0.13 / 0.14 / 0.15.1) | ✅ F4 with real linear algebra + FGLM + experimental F5; scale performance see #6 |
| 4 | ~~Rational polynomials / partial fractions~~ (completed 0.12) | ✅ done — `RationalPolynomial` type + partial fractions + resultant + Karatsuba multiplication |
| 5 | ~~Multi-output optimization / codegen~~ (done in 0.15) | ✅ done — multi-output JIT (97×/21×) + f32 mixed precision + CSE/const-folding/stack-compaction |
| 6 | Gröbner performance at scale (cyclic-6 ℤ_p < 5 s) | 🟡 0.15.2 done (9970 s → 3670 s, 2.7×); <5 s scoped to 0.19 F5 signature reduction (Phase B++) |
| 7 | ~~Arbitrary multivariate (≥3 variables) factorization~~ (completed 0.16) | ✅ done — Wang EEZ lifting + LC preprocessing (constant LC) + Zassenhaus; non-constant LC imposition see #7a |
| 7a | ~~Non-constant leading-coefficient imposition + multivariate sparsity~~ (completed 0.16.1/0.16.2) | ✅ done — mod-p Hensel imposition + sparse Diophantine + field Wang preprocessing on the Fp path |
| 8 | ~~Algebraic-number-field factorization~~ (completed 0.17) | ✅ done — Trager algorithm (shifted norm + ℚ factorization + GF(p^d) modular GCD), univariate path; multivariate extension deferred |
| 9 | ~~Numerical integration / dual numbers / tensor basics / fuel~~ (done in 0.18) | ✅ Done — Vegas + HyperDual + index contraction + fuel; 0.18.1 backfilled the Python/C bindings |
| 10 | ODE solvers (Phase B++ 0.20) | 🟢 SageMath/SymPy parity; first/second-order + systems + series + Laplace |
| 11 | Number theory stack (Phase B++ 0.21) | 🟢 SageMath/PARI parity; modular GCD + integer factorization + primality + discrete log + CRT |
| 12 | Full tensor canonicalisation + specialized pattern transformers (Phase B++ 0.22) | 🟡 Symbolica's last bastion; needs graph-isomorphism engine |
| 13 | Algebraic-geometry tooling (Phase B++ 0.23) | 🟢 SageMath/Singular parity; ideal ops + RUR + primary decomposition + Hilbert series |
| 14 | PDE solvers (Post-1.0) | 🟢 high user demand; Poisson/heat/wave |

---

## 6. Overall Assessment

Execution quality of 0.1 → 0.18.1 is high: every roadmap deliverable shipped,
the layered architecture is clean (no cycles), the 13-crate workspace is
strictly layered, quality gates are strict (`-D warnings` + deny + Miri
awareness), and docs/bindings/CI are well-engineered. The three hard
algorithms planned before 1.0 — polynomial factorization (0.11), Gröbner F4
(0.13, real linear algebra fixed in 0.15.1), and Risch symbolic integration
(0.14) — are all closed and continuously regressed via the SymPy/Symbolica
cross-verification framework.

Realistic positioning: oCAS today is "a high-performance, self-contained
algebra kernel with feature parity against SymPy and near-complete coverage
of Symbolica's example surface". Concretely it ships Risch symbolic
integration, univariate/bivariate/arbitrary-multivariate factorization (plus
algebraic-number-field via Trager), rational functions, Gröbner F4 with real
linear algebra, multi-output JIT / streaming evaluation, Vegas numerical
integration, hyper-dual forward AD, tensor basics, and fuel-based resource
control. Re-measured 0.15.1 performance (still representative): F4 cyclic-5
ℤ₁₃ 23 ms; x³⁰−1 square-free factorization 39 µs (SymPy full factor ~0.9 ms,
~24×); JIT 97× single-output, 21× three-output.

Phase B+ "Closing the Symbolica Gap" (0.15.2 → 0.18.0) is **complete**: every
Symbolica example-domain gap that was open at 0.15.1 — arbitrary multivariate
factorization, algebraic-number-field factorization, numerical integration,
dual numbers, tensors, and fuel — is now closed. Phase B++ "Competitive
Alignment" (0.19.0 → 0.23.0, see EVOLUTION_PLAN) then targets the remaining
gaps before the 1.0.0 freeze: Gröbner performance at the cyclic-6 scale (F5
signature reduction, 0.19), ODE solvers (SageMath/SymPy parity, 0.20), number
theory (SageMath/PARI parity, 0.21), full tensor canonicalisation + advanced
pattern matching (Symbolica's last bastion, 0.22), and algebraic-geometry
tooling (SageMath/Singular parity, 0.23). After Phase B++, 1.0.0 is strictly
**stabilization and release engineering only** (API freeze, coverage,
migration guides, signed artifacts).

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
| 0.18.1 | 2026-07-23 | **Full re-evaluation** after 0.16–0.18 landed. Code-scale snapshot refreshed to 112 files / ~40.9k lines (+33% vs 0.15.1's 95 files / ~30.7k; +127% vs 0.10's ~18k). §1 version table extended through 0.18.1 (0.16.0–0.18.1 rows added). §3 polynomial GCD upgraded 🟡→🟢 (arbitrary-arity multivariate GCD via EEZ [0.16] + modular number-field GCD [0.17]). §4.1 Symbolica gap table rewritten: numerical integration / tensors / duals / fuel all upgraded 🔴→✅ (closed in 0.18); factorization row notes ANF done (0.17); pattern-transformer row added (🟡, `Transformer::Partition` missing); closing paragraph rewritten — all Symbolica example-domain gaps closed except scale Gröbner + full tensor canonicalisation. §4.3 SymPy factorization upgraded 🟡→🟢 (arbitrary-multivariate parity, 0.16). §5 added #11 (tensor canonicalisation + specialized pattern transformers, Post-1.0); header rewritten — Phase B+ declared complete. §6 overall assessment rewritten — 1.0 is stabilization/release-engineering only. Multiple mojibake characters fixed throughout. |
