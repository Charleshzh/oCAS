# Gap Analysis: oCAS vs Reference Systems

This document tracks the implementation completeness of oCAS milestone by
milestone (0.1 → 1.0+) and the gap against the three reference systems:
**Symbolica** (Rust), **SageMath** (Python ecosystem), and **SymPy** (pure
Python). It is a living document and must be refreshed at every version bump.
For the Chinese edition, see [GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md).

> Last evaluated: **0.26.0 @ 2026-08-06** (competitor versions re-verified from
> authoritative sources: Symbolica 2.2.0 (unchanged), SymPy 1.14.0 (unchanged),
> SageMath 10.9 (2026-05-05), FLINT 3.6.0 (2026-06-29; Kinoshita-Li series
> composition, padic_radix, subresultant resultants), msolve 0.10.1 (2026-07-08;
> Gebauer-Möller improvements, QQ lifting fixes), GiNaC 1.8.10 (unchanged),
> mathcore 0.3.1 (previous 0.5.0 record corrected — it does not exist on
> crates.io or GitHub), Numerica (no tagged releases; active development);
> full local re-benchmark of oCAS/Symbolica/SymPy + measured msolve 0.10.1 in
> WSL2: cyclic-6 grevlex 55.04 ms meets the <0.5 s milestone, msolve measured
> 4 ms; DoubleFloat gap closed by 0.24.0 DoubleF64; §5 reprioritized with
> DoubleFloat and cyclic-6 moved to completed items)

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
| 0.19.0/0.19.1 | 1.0 Candidate | ✅ | ✅ F5 Gröbner signature reduction (cyclic-6 ℤ₁₃ 3670 s → 2.63 s, ~1400×); `MonomialOrder` trait refactor + `WeightOrder`/`BlockOrder` |
| 0.20.0/0.20.1 | 1.0 Candidate | ✅ | ✅ Full ODE solver: five first-order methods + integrating factors; second-order constant-coefficient/Cauchy-Euler + VOP + reduction of order + extended undetermined coefficients; series recursion + Frobenius; Laplace IVP (`dsolve_ivp`); 2×2 systems (`dsolve_system`); Python/C bindings; 31 substitution-verified correctness tests |
| 0.21.0 | 1.0 Candidate | ✅ | ✅ Number theory & computational algebra stack: multi-modulus CRT accumulator, BPSW primality + deterministic MR below 2⁶⁴, integer factorization (trial / Brent rho / Pollard p−1 / Williams p+1 / ECM Suyama-Montgomery), BSGS + Pohlig-Hellman discrete logarithms, φ/μ/τ/σ_k/λ functions; univariate Brown modular GCD (`gcd::modular::gcd_modular_z`) + full Brown rewrite of the bivariate `gcd_modular` (content separation + monic interpolation images + multi-prime CRT + rational reconstruction); Python/C bindings (`ocas::ntheory` / `ocas_ntheory_*`); latent integer-sqrt performance bomb in `rational_reconstruction` fixed; ECM factors a 30-digit semiprime in 1.1 s (<10 s) |
| 0.22.0 | 1.0 Candidate | ✅ | ✅ Tensor canonicalisation + advanced pattern matching: graph-iso canonical labelling, tensor canonical form, Young projector, Partition transformer, multi-pattern replace, backtracking matcher |
| 0.23.0 | 1.0 Candidate | ✅ | ✅ Algebraic-geometry tooling: ideal ops (contains/sum/product/quotient/saturate/intersection), MatrixOrder elimination + eliminate(), zero-dimensional solving (Sturm), primary decomposition, radical, Hilbert series/dimension/degree, rational-root theorem; Python/C bindings |
| 0.24.0 | 1.0 Candidate | ✅ | ✅ Heuristic integration module (four techniques: integration by parts LIATE, trig substitution, Weierstrass, Euler [placeholder]) wired into `try_risch_or_fallback` + `integrate_heuristic` API; **DoubleF64** (Dekker/Knuth double-float, ~31 decimal digits, transcendentals, `EvaluationDomain`); Python/C bindings |
| 0.25.0 | 1.0 Candidate | ✅ | ✅ **Multi-modular Gröbner bases** (ℚ ideals: parallel F5 lucky-prime images + CRT + rational reconstruction + exact ℚ verification + trace-free p-adic Hensel lift + fallback); `Algorithm::Auto` routing; F5 speedups (DivisorIndex, bucketed syzygies, parallel row construction, two-phase echelon; cyclic-6 ℤ₁₃ 2.63 s → 1.415 s); **parallel modular GCD**; katsura-6/7 pre-existing gap recorded |
| 0.26.0 | 1.0 Candidate | ✅ | ✅ **SWAR-packed monomial F5 fast path** (u128 packed monomials, n_vars ≤ 8, exponent < 2¹⁵, auto fallback) + echelon rework (i32 coefficients, clone-free two-phase) + **grevlex benchmark variants**; fixed inverted degree/weight direction in graded monomial orders (pre-existing); cyclic-6 ℤ₁₃ grevlex **52.07 ms** (criterion median), Lex 936 ms; cyclic-7 grevlex single round 5.755 s (209 basis elements) |

All 0.1–0.26.0 deliverables landed. The workspace is pinned at 0.26.0. Quality
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
| Gröbner basis | F4 with real linear algebra (0.15.1) + F5 signature reduction (0.19.0: `Signature`/`SyzygySet` + native ℤ_p fast path `f5_fp`) + FGLM + unified `groebner_basis()` dispatch + native i64 ℤ_p pipeline + multi-modular ℚ pipeline (0.25) + u128 SWAR-packed F5 fast path (0.26); cyclic-6 ℤ₁₃ **55.04 ms** grevlex (2026-08-06 measured, criterion median; baseline 2.63 s at 0.19.0, ~48×); cyclic-5 ℤ₁₃ 8.97 ms grevlex | 🟢 F4 + F5 + MultiModular complete |
| Symbolic integration | Risch (elementary transcendental towers + RDE polynomial fragment) + rational-function Hermite + trig exp(I·x) + special-function table (erf/Ei/Si/Ci/Fresnel) + 0.24 heuristic module (parts/trig-sub/Weierstrass/Euler [placeholder]) in `try_risch_or_fallback`; falls back to `Integral(...)`; **gap**: Symbolica 2.2 Rubi port (7000+ rules, 72,944-problem corpus) far wider | 🟢 Risch + heuristics done, Rubi breadth gap |
| Algebraic-geometry tooling | ideal ops (contains/sum/product/quotient/saturate/intersection), elimination orders + eliminate(), zero-dimensional solving, primary decomposition, radical, Hilbert series/dimension/degree, rational-root theorem (0.23) | 🟢 Complete (0.23) |
| Advanced pattern matching | backtracking matcher + multi-pattern replace + Partition transformer + graph-iso canonical labelling (0.22); **gap**: Symbolica `opt`/`alt`/attribute filters more mature (Rubi-grade) | 🟡 Basic-usable, below Symbolica Rubi grade |
| Real root isolation | Sturm sequence + interval isolation + refine (univariate); known gap: only 8/10 roots isolated on expanded Wilkinson n=10 | 🟡 Fairly complete |
| Polynomial GCD | GCD + primitive part + extended GCD (0.12) + arbitrary-arity multivariate GCD via EEZ (0.16) + modular number-field GCD over GF(p^d) with CRT + rational reconstruction (0.17) + univariate Brown modular GCD and bivariate multi-prime modular GCD (0.21, no coefficient explosion on large inputs) | 🟢 Complete (incl. modular fast path, no HEVMGCD) |
| Linear solving | Rational/integer linear systems + bivariate Diophantine (`ax+by=c`) | 🟡 Usable, limited scale |
| JIT evaluation | Cranelift backend; ≥10x speedup target met (per roadmap criterion) | 🟢 Complete |
| ODE | `ocas-calc::ode`: `dsolve()` entry + `classify_ode()` classifier; 5 first-order + integrating factors; 2 second-order + VOP + reduction of order; series recursion + Frobenius; Laplace IVP (`dsolve_ivp`); 2×2 systems (`dsolve_system`); Python/C bindings | 🟢 Complete (0.20.1) |
| Number theory | Multi-modulus CRT, BPSW primality + deterministic MR below 2⁶⁴, integer factorization (rho/p−1/p+1/ECM, 30-digit semiprime in 1.1 s), BSGS + Pohlig-Hellman discrete logarithms, φ/μ/τ/σ_k/λ, quadratic-residue symbols and modular square roots; Python/C bindings (0.21) | 🟢 Core stack complete (0.21) |

---

## 4. Gap Analysis vs Reference Systems

### 4.1 vs Symbolica (Rust, source-available commercial)

Symbolica 2.2.0 (2026-07-24) is oCAS's closest competitor. It moved from
AGPL-3.0 to a source-available commercial license in early 2026 (free for
single-core non-commercial use) and split out MIT crates (Numerica, Graphica,
symbolica-integrate). 2.2 ships the Rubi rule port (7000+ rules, MIT crate)
and SymJIT/CUDA/WASM/C++/ASM code generation.

| Capability | oCAS | Symbolica |
|---|---|---|
| Polynomial factorization | ✅ univariate ℤ/ℤ_p (CZ + Hensel + Zassenhaus) + arbitrary multivariate (0.16 Wang EEZ + non-constant LC imposition 0.16.1/0.16.2) + algebraic-number-field (0.17 Trager, univariate) | ✅ full (arbitrary multivariate + algebraic number fields, `factorization.rs`) |
| Rational polynomials | ✅ `RationalPolynomial<D,O>` with GCD canonicalization | ✅ `rational_polynomial.rs` |
| Partial fractions | ✅ `apart()` / `together()` over Euclidean domains | ✅ `partial_fraction.rs` |
| Rational reconstruction | ✅ `rational_reconstruction(a, m)` via extended Euclidean | ✅ `rational_reconstruction.rs` |
| **Symbolic integration** | 🟢 **Risch + 0.24 heuristics, narrower coverage** | ✅ **Rubi port (7000+ rules, 72,944-problem corpus)** |
| Numerical integration | ✅ Vegas adaptive Monte Carlo + `integrate_1d` + `StatisticsAccumulator` (0.18) | ✅ `numerical_integration.rs` |
| Streaming API | ✅ `StreamingEvaluator`: chunked input + reused stack, constant memory over 1M rows | ✅ `streaming.rs` |
| Tensors / dual numbers | ✅ graph-iso canonical labelling + tensor canonical form + Young projector (0.22); `HyperDual<T>` forward AD (0.18) | ✅ full graphica-based canonicalisation |
| Optimization / codegen | ✅ multi-output JIT (`compile_multi` + CSE + const folding + stack compaction) + f32 mixed precision | ✅ SymJIT + CUDA/WASM/C++/ASM export |
| **DoubleFloat** | ✅ **DoubleF64 (0.24, ~31 decimal digits)** | ✅ **~31 digits, >3× faster than arbitrary precision** |
| Gröbner bases | ✅ F4 + F5 + MultiModular (0.25) + packed F5 (0.26); cyclic-6 ℤ₁₃ grevlex **55.04 ms** (2026-08-06 measured) | ✅ industrial-grade, ~1 s |
| ODE solvers | ✅ complete (0.20.1) | 🔴 none |
| Number theory | ✅ core stack complete (0.21) | 🔴 none |
| Algebraic-geometry tooling | ✅ ideal ops + primary decomposition + Hilbert series (0.23) | 🔴 none |
| Equality saturation | ✅ egg integration | 🔴 none |
| Pattern matching | 🟡 backtracking matcher (0.22); `opt`/`alt`/attribute filters missing | ✅ Rubi-grade matcher |
| Resource control (fuel) | ✅ `Fuel = Arc<AtomicUsize>` + `simplify_with_fuel`/`integrate_with_fuel` (0.18) | ✅ `fuel_backend.rs` |
| **License** | ✅ **LGPL-3.0+, embeddable in commercial products** | ⚠️ **source-available commercial** |

Symbolica 2.2's core strengths — Rubi integration, SymJIT/CUDA/WASM code
generation, DoubleFloat — have been partially closed (DoubleFloat → DoubleF64
in 0.24). Remaining gaps: **symbolic-integration breadth** (Rubi), **code
generation targets** (CUDA/WASM/C++/ASM), and **Gröbner performance at scale**
(see §4.4 vs msolve). oCAS leads in ODE solvers, number theory, algebraic
geometry, equality saturation, and license flexibility (LGPL vs
source-available commercial).

### 4.2 vs SageMath (Python ecosystem)

SageMath is a "Swiss-army-knife" scientific environment. The gap is
**breadth-level**.

| Domain | oCAS | SageMath |
|---|---|---|
| Algebraic geometry | 🟡 basic Gröbner | ✅ Singular integration |
| Number theory | � core stack complete (0.21: CRT + factorization + primality + discrete log + number-theoretic functions) | ✅ PARI/FLINT full stack |
| Differential equations | 🟢 first/second-order/systems/series/Laplace/bindings complete (0.20.1) | ✅ full ODE/PDE solvers |
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
| **Performance** | 🟢 **oCAS advantage** | Rust + Cranelift JIT + arena vs pure Python; measured 2026-08-06: parse 100×, simplify 124×, series 2,550×, integrate 39–76×, eval 183× (SymPy 1.14.0, same inputs); **exception**: `factor(x^30-1)` — SymPy ~50× faster (cyclotomic fast path) |
| Python ergonomics | 🟢 parity | oCAS has `ocas-py` bindings |

The 0.6.0 success criterion — "parity with SymPy on basic polynomial,
calculus, and rewriting" — is met and exceeded on the **performance** axis,
and **integration** was closed by Risch in 0.14 while **factorization**
reached arbitrary-multivariate parity in 0.16 (plus algebraic-number-field
via Trager in 0.17). The remaining feature gap against SymPy is the **breadth
of integration heuristic fallbacks** (SymPy's `manualintegrate` /
heuristic pool is wider than oCAS's Risch + table path).

### 4.4 vs msolve (Gröbner performance benchmark)

msolve is the open-source Gröbner performance benchmark (F4/F5 + modular
arithmetic + Hensel lifting + Berlekamp-Massey). **Measured 2026-08-06 on this
machine** (WSL2, msolve 0.10.1 built from source, `-g 2` GB-only mode, DRL
order, single-threaded; basis sizes match oCAS exactly):

| Benchmark | msolve 0.10.1 (measured) | oCAS 0.26 (measured) | Ratio |
|---|---|---|---|
| cyclic-5 ℤ₁₃ | 3 ms | 8.97 ms (grevlex) | 3.0× |
| cyclic-6 ℤ₁₃ | **4 ms** | **55.04 ms** (grevlex) | 13.8× |
| cyclic-7 ℤ₁₃ | **55 ms** | 3.829 s (grevlex, single round) | ~70× |
| katsura-6 | 3 ms | not measured (pre-existing gap) | — |
| katsura-7 | 7 ms | not measured (pre-existing gap) | — |

The 0.26.0 packed-monomial F5 brought cyclic-6 to 55 ms (~48× over the 0.19.0
baseline of 2.63 s and meeting the <0.5 s milestone); cyclic-7 (~70×) and the
katsura family (oCAS unfinished vs 3–7 ms) remain the largest single gaps.
The previously cited msolve values (cyclic-6 0.04 s) have been replaced by
these measurements (4 ms).

**Windows availability**: msolve builds on Windows only via MSYS2 (GMP/MPFR/FLINT);
oCAS's native Windows support remains a differentiation point.

---

## 5. Key Gaps & Priorities

Ranked by impact × implementation cost. All hard-algorithm gaps planned before
1.0 are **closed**; Phase B+ "Closing the Symbolica Gap" (0.15.2–0.18.0) and
Phase B++ "Competitive Alignment" (0.19.0–0.23.0) plus the 0.24–0.26 round
(heuristics/DoubleF64, MultiModular, packed F5) are complete. Remaining open
items are re-ranked below (2026-08-06).

| # | Gap | Priority |
|---|---|---|
| 1 | ~~Full polynomial factorization~~ (completed 0.11.0–0.11.1) | ✅ done — univariate and bivariate (monic-in-x) closed; ≥3 variables see #7 |
| 2 | ~~Risch symbolic integration~~ (completed 0.14) | ✅ done — elementary transcendental towers + RDE polynomial fragment + rational Hermite + special-function table |
| 3 | ~~Gröbner F4/F5~~ (completed 0.13 / 0.14 / 0.15.1) | ✅ F4 with real linear algebra + FGLM + experimental F5; scale performance see #6 |
| 4 | ~~Rational polynomials / partial fractions~~ (completed 0.12) | ✅ done — `RationalPolynomial` type + partial fractions + resultant + Karatsuba multiplication |
| 5 | ~~Multi-output optimization / codegen~~ (done in 0.15) | ✅ done — multi-output JIT (97×/21×) + f32 mixed precision + CSE/const-folding/stack-compaction |
| 6 | ~~Gröbner performance at scale (cyclic-6 ℤ_p < 5 s)~~ (completed 0.19) | ✅ done — F5 signature reduction (0.19.0): cyclic-6 ℤ₁₃ 3670 s → **2.63 s** (~1400×); F4/F5 unified dispatch; generic-domain + native ℤ_p fast path both verified |
| 7 | ~~Arbitrary multivariate (≥3 variables) factorization~~ (completed 0.16) | ✅ done — Wang EEZ lifting + LC preprocessing (constant LC) + Zassenhaus; non-constant LC imposition see #7a |
| 7a | ~~Non-constant leading-coefficient imposition + multivariate sparsity~~ (completed 0.16.1/0.16.2) | ✅ done — mod-p Hensel imposition + sparse Diophantine + field Wang preprocessing on the Fp path |
| 8 | ~~Algebraic-number-field factorization~~ (completed 0.17) | ✅ done — Trager algorithm (shifted norm + ℚ factorization + GF(p^d) modular GCD), univariate path; multivariate extension deferred |
| 9 | ~~Numerical integration / dual numbers / tensor basics / fuel~~ (done in 0.18) | ✅ Done — Vegas + HyperDual + index contraction + fuel; 0.18.1 backfilled the Python/C bindings |
| 10 | ~~ODE solvers~~ (0.20.1 complete) | ✅ Done — 5 first-order + integrating factors; 2 second-order + VOP + reduction of order; series recursion + Frobenius; Laplace IVP; 2×2 systems; Python/C bindings |
| 11 | ~~Number theory stack~~ (done 0.21) | ✅ done — modular GCD (univariate Brown + bivariate multi-prime) + integer factorization (ECM: 30-digit semiprime in 1.1 s) + BPSW primality + discrete log + CRT + number-theoretic functions + Python/C bindings |
| 12 | Full tensor canonicalisation + specialized pattern transformers (Phase B++ 0.22) | � delivered in 0.22.0: graph-iso canonical labelling + tensor canonical form + Young projector + Partition + multi-pattern replace |
| 13 | Algebraic-geometry tooling (Phase B++ 0.23) | 🟢 SageMath/Singular parity; ideal ops + RUR + primary decomposition + Hilbert series |
| 14 | PDE solvers (Post-1.0) | 🟢 high user demand; Poisson/heat/wave |
| 15 | ~~DoubleFloat~~ (completed 0.24) | ✅ done — **DoubleF64** (Dekker/Knuth double-float, ~31 decimal digits, transcendentals, `EvaluationDomain`, Python/C bindings); former P2 gap closed |
| 16 | ~~Gröbner cyclic-6 < 0.5 s~~ (completed 0.26) | ✅ done — packed-monomial F5 + echelon rework; cyclic-6 ℤ₁₃ grevlex **55.04 ms** (2026-08-06 criterion median); ratio vs measured msolve 4 ms = 13.8× tracked separately |
| 17 | **Gröbner scale: katsura-6/7 + cyclic-7** (open, P1) | ⚠️ katsura-6/7 unfinished (single round >30 min, pre-existing); cyclic-7 Lex >2 h; cyclic-7 grevlex 3.829 s vs measured msolve 55 ms (~70×); target katsura-6 < 1 s, cyclic-7 tractable — extend packed pipeline + multi-modular strategy |
| 18 | **Symbolic-integration breadth** (open, P0) | ⚠️ Symbolica Rubi port (7000+ rules, 72,944-problem corpus) vs oCAS Risch + 0.24 heuristics |
| 19 | **Code generation targets** (open, P1) | ⚠️ LLVM JIT + CUDA/WASM export vs Symbolica SymJIT/CUDA/WASM/C++/ASM |
| 20 | **Matrix / linear algebra** (open, P2) | ⚠️ SymPy 1.14 DomainMatrix 10000× rref speedup + Smith normal form vs oCAS Bareiss |
| 21 | **Windows FLINT support** (open, P2) | ⚠️ `flint` feature Linux/WSL only; target 3-platform |
| 22 | **Tensor nested-function handling** (open, P3) | 🟡 nested-function tensors vs Symbolica Graphica engine maturity |
| 23 | **Quadratic sieve factorization** (open, P3) | 🟡 ECM 30-digit in 1.1 s vs SymPy `qs_factor` on large composites |

---

## 6. Overall Assessment

Execution quality of 0.1 → 0.26.0 is high: every roadmap deliverable shipped,
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
integration (plus 0.24 heuristics), univariate/bivariate/arbitrary-
multivariate factorization (plus algebraic-number-field via Trager), rational
functions, Gröbner F4/F5/MultiModular with packed fast paths, multi-output JIT
/ streaming evaluation, Vegas numerical integration, hyper-dual forward AD,
tensor canonicalisation, ODE solvers, a number-theory stack, algebraic-geometry
tooling, and fuel-based resource control. Re-measured 2026-08-06: cyclic-6
ℤ₁₃ grevlex 55.04 ms (msolve 4 ms, ratio 13.8×); JIT 73.8× single-output,
22.6× three-output; parse/simplify/series vs SymPy 100×/124×/2,550×.

Phase B+ "Closing the Symbolica Gap" (0.15.2 → 0.18.0) and Phase B++
"Competitive Alignment" (0.19.0 → 0.23.0) are **complete**, and the 0.24–0.26
round closed DoubleFloat (→ DoubleF64) and the cyclic-6 <0.5 s milestone
(55 ms). Remaining gaps before the 1.0.0 freeze: symbolic-integration breadth
(Rubi), Gröbner scale (katsura-6/7 unfinished; cyclic-7 ~70× vs msolve),
code generation targets (CUDA/WASM/C++), matrix/linear algebra (DomainMatrix),
and Windows FLINT support.

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
| 0.19.0 | 2026-07-23 | **F5 Gröbner basis released — cyclic-6 scale gap closed.** §3 Gröbner row upgraded 🟡→🟢 (F5 signature reduction). §4.1 Gröbner bases competitor row upgraded 🟡→🟢. §5 #6 (Gröbner performance at scale) marked ✅ done — cyclic-6 ℤ₁₃ 3670 s → **2.63 s** (~1400×) via `f5_fp` native ℤ_p fast path; cyclic-5 0.05 s; generic-domain + ℤ_p paths both verified. Unified `groebner_basis()` dispatch (`Algorithm::{Auto,F4,F5,Buchberger}`). Multi-order (`WeightOrder`/`BlockOrder`) deferred to 0.19.1 (trait refactor). |
| 0.19.1 | 2026-07-23 | **MonomialOrder trait refactor + WeightOrder/BlockOrder released.** `Copy` + static dispatch → `Clone + Default` + method dispatch (`&self`); `PhantomData<O>` → `order: O` field; new `WeightOrder` (weighted) and `BlockOrder` (block) orderings with `SubOrder` enum; all 11 `O::cmp` call sites updated; `Signature::cmp_pot` updated with `order: &O` param. Multi-order support upgraded `[~]`→`[x]`. |
| 0.20.0 | 2026-07-27 | **ODE solver released.** §3 added ODE row (🟡). §4.2 SageMath ODE row upgraded 🔴→🟡. §5 #10 (ODE solver) marked 🟡 partially complete — 5 first-order (separable/linear/Bernoulli/exact/homogeneous) + 2 second-order (const-coeff/Cauchy-Euler) + power-series framework + `classify_ode()` classifier + `dsolve()` entry; Laplace transform, ODE systems, Python/C bindings deferred. Version bumped to 0.20.0. |
| 0.20.1 | 2026-07-27 | **ODE solver backfill complete.** Integrating factors (μ(x)/μ(y)); variation of parameters (VOP, fixing silently-dropped Cauchy-Euler forcing); reduction of order; power-series coefficient recursion + Frobenius (real rational indicial roots); Laplace IVP (`dsolve_ivp`); 2×2 constant-coefficient systems (`dsolve_system`); Python/C bindings (`classify_ode`/`dsolve`/`dsolve_ivp`). Fixed `real_roots` isqrt/formula bugs, hardcoded `is_exact`, unnormalized Cauchy-Euler coefficients, integrator (ax+b)⁻¹ and fractional-power gaps, `substitute_solution` missing bare y(x), and series-coefficient `diff` pollution. Added the `collect_terms` like-term collector + `expand`. 31 substitution-verified correctness tests (3 known gaps ignored). ODE gap upgraded 🟡→🟢. Version bumped to 0.20.1. |
| 0.26.0 re-eval | 2026-08-06 | **Competitor versions re-verified + full local re-benchmark.** Header updated to 0.26.0. Competitor re-verification from authoritative sources: FLINT 3.5.0→3.6.0 (Kinoshita-Li series composition, padic_radix, subresultant resultants), msolve 0.7.x→0.10.1 (Gebauer-Möller, QQ lifting fixes), mathcore corrected to 0.3.1 (0.5.0 does not exist), Numerica has no tagged releases, SageMath date corrected to 10.9@2026-05-05; Symbolica/SymPy/GiNaC unchanged. §1 version table extended through 0.26.0 (0.22.0–0.26.0 rows added). §3 Gröbner/integration rows updated (packed F5, MultiModular, 0.24 heuristics); algebraic-geometry + pattern-matching rows added. §4.1 rewritten (source-available license, Rubi, DoubleFloat closed by DoubleF64 in 0.24, measured cyclic-6 55.04 ms). §4.3 performance row refreshed with 2026-08-06 measurements (parse 100×, series 2,550×, factor(x^30-1) exception: SymPy ~50× faster). New §4.4 vs msolve with measured WSL2 values (cyclic-6 4 ms, cyclic-7 55 ms, katsura 3–7 ms; basis sizes match oCAS). §5 items 15–16 marked complete (DoubleFloat, cyclic-6 <0.5 s), items 17–22 re-ranked (katsura/cyclic-7 P1, integration breadth P0, codegen P1, matrix P2, Windows FLINT P2, tensor+quad-sieve P3). §6 overall assessment refreshed with 2026-08-06 measurements. |
