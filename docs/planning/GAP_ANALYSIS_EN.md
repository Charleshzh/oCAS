# Gap Analysis: oCAS vs Reference Systems

This document tracks the implementation completeness of oCAS milestone by
milestone (0.1 → 1.0+) and the gap against the three reference systems:
**Symbolica** (Rust), **SageMath** (Python ecosystem), and **SymPy** (pure
Python). It is a living document and must be refreshed at every version bump.
For the Chinese edition, see [GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md).

> Last evaluated: **0.11.0 @ 2026-07-03**

---

## Legend

| Mark | Meaning |
|---|---|
| ✅ | Complete |
| 🟡 | Basic / partial |
| 🔴 | Missing or major gap |
| ⚠️ | Complete with caveats |

---

## 1. Version Completion Status (0.1–0.10)

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

All 0.1–0.11.1 deliverables landed. The workspace is pinned at 0.11.1. Quality
gates are green: `cargo fmt`, `clippy -D warnings`, workspace tests,
`cargo deny`, 77 pytest cases, `mdbook build`.

---

## 2. Code Scale

Snapshot of `src/` Rust lines (excluding tests and benches).

| Crate | Files | Lines |
|---|---|---|
| ocas-poly | 9 | ~3,200 |
| ocas-eval | 11 | ~2,525 |
| ocas-domain | 9 | ~2,115 |
| ocas-rewrite | 7 | ~1,719 |
| ocas-py | 7 | ~1,546 |
| ocas-calc | 7 | ~1,393 |
| ocas-c | 3 | ~894 |
| ocas-core | 5 | ~1,150 |
| ocas-atom | 2 | ~864 |
| ocas-parse | 3 | ~565 |
| ocas (prelude) | 1 | ~113 |
| ocas-gpl | 1 | 1 (placeholder) |
| **Total src** | **64** | **~16k** |

`ocas-gpl` is a placeholder; GPL-exclusive backends are Post-1.0 work, in line
with the roadmap.

---

## 3. Algorithm Depth Audit

This section is the single most decisive factor in CAS maturity and the main
source of the gap.

| Algorithm Area | oCAS Status | Maturity |
|---|---|---|
| Polynomial factorization | `factor()` on `DenseUnivariatePolynomial` over ℤ and ℤ_p, plus bivariate `factor()` on `SparseMultivariatePolynomial` over ℤ and ℤ_p (monic-in-x Wang Hensel) | 🟢 Fairly complete |
| Gröbner basis | Classic Buchberger + minimize/auto-reduce; **no** F4/F5, no heuristics | 🟡 Basic |
| Symbolic integration | Heuristic table (power/inverse/sin/cos/exp/linear subst); falls back to `Integral(...)`; **no** Risch | 🟡 Basic |
| Real root isolation | Sturm sequence + interval isolation + refine (univariate) | 🟢 Fairly complete |
| Polynomial GCD | GCD + primitive part; no modular GCD / EEA optimization | 🟡 Usable |
| Linear solving | Rational/integer linear systems + bivariate Diophantine (`ax+by=c`) | 🟡 Usable, limited scale |
| JIT evaluation | Cranelift backend; ≥10x speedup target met (per roadmap criterion) | 🟢 Complete |

---

## 4. Gap Analysis vs Reference Systems

### 4.1 vs Symbolica (Rust, AGPL)

Symbolica's `examples/` directory reveals the maturity gap. oCAS is roughly an
early functional subset of Symbolica.

| Capability | oCAS | Symbolica |
|---|---|---|
| Polynomial factorization | ✅ `factor()` over ℤ and ℤ_p (CZ + Hensel + Zassenhaus); bivariate factorization over ℤ and ℤ_p (Wang Hensel, monic-in-x) | ✅ full (`factorization.rs`) |
| Rational polynomials | 🟡 partial | ✅ `rational_polynomial.rs` |
| Partial fractions | 🔴 none | ✅ `partial_fraction.rs` |
| Rational reconstruction | 🔴 none | ✅ `rational_reconstruction.rs` |
| Numerical integration | 🔴 none | ✅ `numerical_integration.rs` |
| Streaming API | 🔴 none | ✅ `streaming.rs` |
| Tensors / dual numbers | 🔴 none | ✅ `tensors.rs` / `dual.rs` |
| Optimization / codegen | 🟡 JIT, f64 only | ✅ `optimize.rs` / multi-output |
| Gröbner basis | 🟡 Buchberger | ✅ industrial grade |

Symbolica's core strengths — industrial factorization, rational function
arithmetic, multi-output optimization, streaming — are largely absent in oCAS.
Symbolica has been refined over years; oCAS must close the hard-algorithm gap
in the ALG layer.

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
| Integration | 🟡 oCAS weaker | SymPy has Risch + heuristics; oCAS heuristic only |
| Factorization | � parity | univariate ℤ and ℤ_p via CZ + Hensel + Zassenhaus; multivariate still in 0.11.1 |
| Gröbner | 🟡 oCAS slightly weaker | both mid-tier, SymPy richer |
| Matrix / linear algebra | 🟢 parity | oCAS has Bareiss determinant/inverse |
| **Performance** | 🟢 **oCAS advantage** | Rust + Cranelift JIT + arena vs pure Python |
| Python ergonomics | 🟢 parity | oCAS has `ocas-py` bindings |

The 0.6.0 success criterion — "parity with SymPy on basic polynomial,
calculus, and rewriting" — is met and exceeded on the **performance** axis,
and factorization is now closed on the univariate side; only **integration**
remains a notable hard-algorithm trail on the SymPy comparison.

---

## 5. Key Gaps & Priorities

Ranked by impact × implementation cost, the hard problems on the road to 1.0.

| # | Gap | Priority |
|---|---|---|
| 1 | ~~Full polynomial factorization~~ (completed 0.11.0–0.11.1) | ✅ done — univariate and bivariate (monic-in-x) closed; unblocks rational functions, partial fractions, solvers |
| 2 | Risch symbolic integration (roadmap: 0.14) | 🔴 hallmark of "can it integrate" |
| 3 | Gröbner F4/F5 (roadmap: 0.13) | 🟡 current Buchberger too slow on large cyclic-n |
| 4 | Rational polynomials / partial fractions (roadmap: 0.12) | 🟡 Symbolica core feature, depends on factorization |
| 5 | Multi-output optimization / codegen (roadmap: 0.15) | 🟡 JIT is f64 single-output; extend to multi-output/multi-precision |

---

## 6. Overall Assessment

Execution quality of 0.1 → 0.10 is high: every roadmap deliverable shipped,
the layered architecture is clean (no cycles), the 12-crate workspace is
strictly layered, quality gates are strict (`-D warnings` + deny + Miri
awareness), and docs/bindings/CI are well-engineered. As a self-developed CAS
of ~16k lines over ~14 months, the foundation is solid.

However, the 1.0 goal of "performance parity or better with Symbolica" still
has a core hard-algorithm shortfall: Risch integration and F4/F5 Gröbner are
the remaining "rites of passage" of a CAS. Factorization has been closed through
0.11.1 (univariate and bivariate monic-in-x), so the next high-value leap is
rational function arithmetic and partial fractions (0.12). On performance (arena + JIT + SIMD) oCAS has a structural
edge, but algorithm depth decides "what can be computed", not "how fast".

Realistic positioning: oCAS today is closer to "a high-performance subset of
SymPy's core, with evaluation performance exceeding SymPy and univariate
factorization parity", rather than a direct Symbolica replacement. Closing
rational function arithmetic before 1.0 is the highest-value leap.

---

## Update Log

Record every refresh here (version, date, evaluator, deltas).

| Version | Date | Deltas |
|---|---|---|
| 0.10.0 | 2026-07-02 | Initial assessment. All 0.1–0.10 deliverables verified complete; gaps against Symbolica / SageMath / SymPy documented; factorization + Risch integration identified as top priorities. |
| 0.11.0 | 2026-07-03 | Polynomial factorization completed (univariate ℤ and ℤ_p); multivariate GCD added; SymPy comparison updated to parity for univariate factorization; highest-priority gap shifted to rational functions / partial fractions (0.12). |
| 0.11.1 | 2026-07-04 | Bivariate factorization over ℤ and ℤ_p (monic-in-x Wang Hensel) added; sparse multivariate `factor()` entry points and C polynomial bindings landed; mdBook factorization chapter added; highest-priority gap remains rational functions / partial fractions (0.12). |
