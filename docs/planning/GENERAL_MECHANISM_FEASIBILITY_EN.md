# General and Correct Integration Mechanisms: Feasibility and Method (2026-09-12)

> This document assesses the "theorem-driven mechanisms instead of rule enumeration"
> route: its feasibility, its method, and its cost. The conclusions rest on three
> kinds of evidence: facts from this repository's code and module documentation,
> the quantitative data from the 0.27.3 investigation (see
> [BENCHMARK_RESULTS_CN.md](BENCHMARK_RESULTS_CN.md) §"0.27.3 follow-up"), and the
> maturity of the published algorithms. Estimates are marked as such.

---

## 1. Criteria first

| Concept | Operational definition |
|---|---|
| **General** | Output comes from a **theorem-covered parameterised family**, not from case-by-case patterns; code size is **decoupled** from the size of the covered input family (Rubi's 6700 rules are the exact opposite) |
| **Correct** | Every output carries a **machine-checkable proof** (a symbolic certificate `D(F) − f ≡ 0` inside the differential field), not sampled numerical evidence; wrong answers become structurally impossible |
| **Complete (relative)** | For the implemented extension fragment the algorithm is a **decision procedure**: either an antiderivative or a **proof** of non-elementarity |
| **Honest** | A three-valued outcome: `Found` / `ProvedNonElementary` / `Unknown`; oCAS currently conflates the last two into "fallback", which loses information |

"General and correct" does not mean "covers everything": it means **whatever fragment you
implement is a decision procedure and its results are trustworthy**. Coverage is a separate
engineering problem and must not be bought by sacrificing correctness (the six 0.27.2
wrong-answer classes are the lesson).

## 2. Theoretical backbone: the algorithms exist and are decision procedures

- **Liouville's theorem + the Risch algorithm**: for an elementary function over an
  elementary transcendental extension, decide whether an elementary antiderivative exists,
  and construct it when it does. **Generality and correctness follow from the theorem**;
  no per-case human verification is needed.
- The complete textbook ladder (Bronstein, *Symbolic Integration I*, with pseudocode):

| Ch. | Content | oCAS today |
|---|---|---|
| 2 | Rational functions: Hermite reduction + Rothstein–Trager | ✅ `integral/rational.rs` |
| 3–5 | Elementary transcendental extensions (primitive `log` / hyperexponential `exp`) | 🟡 partial (see §3) |
| 6 | The Risch differential equation `Dq + fq = g` and coupled systems | 🟡 **polynomial solutions only** (`integral/rde.rs` says so itself) |
| 7–8 | **Algebraic extensions** and regular towers (Trager's integral basis) | ❌ not wired in (factorisation exists, see §3.5) |
| 9 | Parallel Risch / engineering the structure theorems | ❌ |
| (non-elementary) | Singer's structure theorems; the most general decidable engineering route is **Meijer G / Slater expansion** (the SymPy `meijerg` route) or **holonomic** (Gosper–Zeilberger style); the logarithmic case is **Li₂ reduction** (Lewin's formulas plus the algorithmic work of Raab/Baddoura) | ❌ |

**Key judgement**: this route's theoretical maturity is "textbook plus published algorithm",
not research frontier. The real difficulty is **engineering scale** and **complexity
control**, not "does an algorithm exist".

## 3. What oCAS already has, and what is missing (from code/docs, not speculation)

### 3.1 What exists — the main reason feasibility is high

- **Differential tower** `ℚ(x, t₁, …, tₙ)`: `tower/build.rs` has `GenKind::{Constant, Log, Exp}`
  with each generator's derivative, and `tower/convert.rs` converts Atom ↔ tower elements.
- **Tower element algebra**: `tower/elem.rs` provides `KElem` (field elements), `KPoly`
  (polynomials in the top generator with lower-field coefficients) and `KRat`, with
  `derivative_dt`, `div_rem`, `gcd`, `eea`, `square_free`, `partial_deriv` — every field
  operation Risch needs.
- **The Risch recursion**: `integral/risch.rs` descends level by level (Hermite reduction per
  level, logarithmic-derivative identity, undetermined coefficients at primitive levels,
  `rde.rs` at hyperexponential levels, the base level delegated to `rational.rs`).
- **Algebraic number fields**: `ocas-poly/src/factor/algebraic.rs` already implements
  **Trager's algorithm** (norm via evaluation–interpolation of a scalar resultant, Trager's
  shift for square-freeness, factor over ℚ then GCD over K) — the core tool the algebraic case needs.
- **Modular algorithms**: the 0.25 `MultiModular ℚ` pipeline (parallel lucky primes + CRT +
  rational reconstruction + traceless Hensel) — the answer to **coefficient blow-up with
  symbolic coefficients** (9 of the current 12 timeouts sit in `symbolic_rational`).
- **Verification infrastructure**: the 0.27.2/0.27.3 numerical oracle
  (`ocas-tests/src/integral_eval.rs`), the 1892 harness, SymPy cross-checks, proptest, and the
  `preserves_argument_order` semantic guardrail.

### 3.2 Gap 1: algebraic generators are rejected outright

`tower/build.rs`: "algebraic functions such as `√x` (non-integer exponents) are rejected".
Consequence: the whole **radical bucket (356 unsolved — one of the largest)** and the algebraic
family never reach the general engine.

### 3.3 Gap 2: dependent generators are rejected instead of merged (impact measured)

`tower/build.rs`: "algebraically dependent generators (e.g. `log(x)` and `log(2x)`, or
`exp(x)` and `exp(x+1)`) are rejected rather than merged".

The real impact of this is badly under-appreciated: trigonometric/hyperbolic integrands are
rewritten to exponentials before entering Risch, so `sinh(u) = (eᵘ − e⁻ᵘ)/2` necessarily
produces `exp(u)` and `exp(−u)`, and `eᵘ·e⁻ᵘ = 1` is **algebraically dependent** → the whole
family is rejected. **Measured: 111 unsolved cases (7.3%) contain hyperbolic heads**
(mixed-other 76 + hyperbolic 35). To a significant degree `hyperbolic_reduction.rs`,
`kernel_subst.rs` and `exp_log.rs` in 0.27.x are **patches around this entry-condition limit** —
patching mechanisms instead of fixing the general engine's front door.

### 3.4 Gap 3: the RDE seeks polynomial solutions only

`integral/rde.rs` states that only polynomial solutions are sought and that rational-function
solutions (which need denominator bounds) are outside the fragment. Consequence: a substantial
part of the `exp`-level integrals is abandoned (formally part of Bronstein ch. 6's full version).

### 3.5 Gap 4: the algebraic integration chain is not wired in

Trager **factorisation** exists, but the integral basis, algebraic Hermite reduction, algebraic
residues and regular-tower construction do not. This is the subject matter of ch. 7–8, and it is
the key to whole families: `1/(1+x⁴)`, `1/(1−3x²+x⁴)`, irreducible sextic denominators (the 0.27.3
probes confirmed all of these decline, and SymPy itself returns `RootSum` for several).

### 3.6 Gap 5: no non-elementary decision, and the function heads are missing

- There is no Meijer-G / Li₂ reduction, and not even a `polylog`/`Li₂` function head.
- **Measured: 122 unsolved cases (8.0%) have reference answers containing `polylog`** → unreachable
  today, independently of mechanism quality.
- With no non-elementarity decision the engine cannot distinguish "proved non-elementary" from
  "I don't know"; the user sees the same `Integral(...)` for both.

### 3.7 Gap 6: the bottleneck is complexity, not budget

9 of 12 timeouts sit in `symbolic_rational`. That is **coefficient blow-up in the field Euclidean
algorithm with symbolic coefficients** — an algorithmic-complexity problem. Repeatedly "adding
budgets/gates" in 0.27.x only converts it into a fast, honest refusal. The right direction is
modular arithmetic (MultiModular already exists) + parallel Risch + lazy series, not more wall clock.

## 4. Component-by-component feasibility

| Component | Theory maturity | oCAS foundation | Effort | Risk | Expected gain (measurable) |
|---|---|---|---|---|---|
| **Certificates + three-valued Outcome** | high (evaluation in a formal differential field) | tower normalisation + numeric oracle exist | **small** | low | hard invariant: wrong answers can no longer escape; provides the gate for every later step |
| **Dependent-generator merging (regular tower)** | high (Trager regularisation) | tower + Trager factorisation exist | medium | medium (must return `Unknown` when dependency is undecidable) | the hyperbolic family, ~111 cases; unblocks the real problem the patch mechanisms were routed around |
| **Complete the transcendental Risch (rational RDE solutions / coupled systems / log-part structure theorem)** | high (textbook pseudocode) | `rde.rs` + logarithmic-derivative identity exist | medium–large | medium | the exp/log families; the elementary part of the `exp-log` bucket (79/83 unsolved) |
| **Algebraic extension integration (integral basis + algebraic Hermite + residues)** | high (ch. 7–8 + Trager) | Trager factorisation + resultants + algebraic-field GCD | **large** | medium–high | the radical bucket (356 unsolved) + the irreducible-denominator families |
| **Non-elementary layer (Li₂/Meijer G + new heads + oracle)** | high (Slater/Lewin classics; SymPy as a reference) | the 0.27.3 "add the head + verify numerically first" pattern is proven | **large** | medium | the 122 `polylog` cases + the rest of `exp-log` + composite `special` shapes |
| **Performance (modular + parallel Risch)** | high | MultiModular exists | medium | low | timeouts and wall clock (currently 12 cases / 461 s) |
| **Inverse-function substitution engine** | medium (rule-ified substitution classification, Rubi 4.x–5.x's idea) | new | large | medium | the `inverse-trig-hyper` bucket (138/146 unsolved) |
| **Generality testing (proptest + certificates)** | high | proptest + oracle exist | small | low | a new metric `certified_rate`, targeted to stay exactly 1.0 |

## 5. Method: make "correct" a discipline and "general" a recursion

1. **One uniform interface, structural recursion**: `trait RischIntegrate` with one implementation
   per extension type (primitive / hyperexponential / algebraic); recursion descends the tower.
   Pattern matching must not replace the structure theorems — that is how the project degenerates
   into Rubi.
2. **Certificates first (the most important point in this assessment)**: the acceptance criterion
   for every output is the **symbolic** certificate `normalize_in_field(D(F) − f) == 0`. The
   numerical oracle is demoted to a probe/regression tool. Only this makes "correct" scale beyond
   human review of each rule.
3. **Three-valued Outcome threaded through the chain**: `Found { value, certificate }` /
   `ProvedNonElementary { witness }` / `Unknown`; `Unknown` must honestly return `Integral(...)`
   and never masquerade.
4. **Regular towers**: merge algebraically dependent generators rather than rejecting them; when
   the dependency is **undecidable**, return `Unknown`.
5. **Reuse the existing assets for the algebraic case**: defining polynomial → integral basis
   (Trager) → algebraic Hermite → residues; `ocas-poly`'s Trager factorisation and resultants are
   directly usable.
6. **Build the non-elementary layer in two steps**: add the function heads plus oracle evaluators
   first (a pattern 0.27.3 already proved out: `Ei`/`Si`/`Ci`/`Shi`/`Chi`/Fresnel all landed that
   way), then add the reductions (Li₂ → Meijer G).
7. **Performance by algorithm, not by budget**: modular arithmetic (MultiModular) + parallel Risch
   + lazy series; wall-clock budgets are for **defence** only, never to hide complexity.
8. **Generality testing is a new metric** (running the Rubi corpus alone is not enough):
   - proptest generates random elements **from the tower** (not random strings) → integrate →
     check the **symbolic certificate**;
   - adversarial families: dependent generators, branch/domain degeneracies, parameter
     degeneracies (e.g. `m = −1` invalidating a formula);
   - a new metric **`certified_rate` = certificates passed / solved**, targeted to stay exactly
     `1.0` and enforced as a CI gate.
9. **Real-valued semantics are annotated separately**: certificates are formal (inside the
   differential field); branch/domain handling lives outside the symbolic layer (`halfpower`'s
   sheet factor is the precedent) and must be recorded explicitly rather than assuming
   "formally correct = real-valued correct".

## 6. Phased route and acceptance gates

**Version mapping (written into ROADMAP / EVOLUTION_PLAN on 2026-09-12)**: P0→**0.28.0**,
P1→**0.29.0**, P2→**0.30.0**, P3→**0.31.0**, P4→**0.32.0**. The former
0.28.0/0.29.0/0.30.0 (Gröbner performance at scale, LLVM JIT, matrix engine + platform
close-out) slide in order to **0.33.0/0.34.0/0.35.0**, and 1.0.0 slides to Month 69.
The two **defects** the 0.27.3 investigation located (residue resolution, expression-level
cycle detection) ship together with 0.28.0.

| Phase | Version | Deliverable | Gate (existing harness + the new metric) |
|---|---|---|---|
| **P0** | **0.28.0** | Defect fixes (residue resolution + expression-level cycle detection) + symbolic certificates + three-valued `Outcome` + dependent-generator merging | certificates exactly 0 for 100% of the 1892 solved set; `certified_rate = 1.0`; new solves in the hyperbolic family (measured baseline: 111 unsolved); the timeout count does not rise |
| **P1** | **0.29.0** | Complete the transcendental Risch (rational RDE solutions + coupled systems + log-part structure theorem) | clear improvement in the `exp-log` bucket; **0 regressions**; certificate gate holds |
| **P2** | **0.30.0** | Algebraic extensions (regular tower + integral basis + algebraic Hermite + residues) + multi-radical bases + inverse-function substitution engine | the `1/(1+x⁴)` / `1/(1−3x²+x⁴)` families solve; the radical and `inverse-trig-hyper` buckets improve; certificate gate holds |
| **P3** | **0.31.0** | Non-elementary layer (`polylog`/`Li₂` heads + oracle + Li₂/Meijer reductions + `ProvedNonElementary`) | the reachable part of the 122 `polylog` cases; the rest of `exp-log`; certificate gate holds |
| **P4** | **0.32.0** | Performance (modular + parallel Risch + lazy series) | wall clock no worse than the 0.28.0 baseline; timeout count falls |

Ordering rationale: P0 is a **zero-mathematical-risk prerequisite** that immediately yields the
hard "correct" invariant; P1 reuses the most existing code; P2 reuses the Trager assets; P3 carries
the largest risk and effort but the most unique payoff (8% of the corpus is completely unreachable
today); P4 comes after capability so that performance budgets cannot mask missing capability.

## 7. Honest limits and risks

1. **Scale**: Bronstein's book is research-length; a complete Risch implementation is measured in
   **person-years**. Delivery must therefore be fragment-by-fragment, while always preserving the
   hard invariant of **never emitting an uncertified result** — capability can be incremental,
   correctness cannot.
2. **Undecidability**: constant identification in the tower and algebraic-dependency testing are
   undecidable in general → `Unknown`, never a guess. This matches a phase that needs zero wrong
   answers rather than higher coverage.
3. **Real-valued branches**: a differential-field certificate does not guarantee that the real
   antiderivative is correct on a given interval (logarithm branches, `|p+q|` cases). Extra
   domain/sheet handling and documentation are required.
4. **Performance**: full Risch is slow; Symbolica's Rubi port is about 4.1× faster than oCAS on the
   1892 subset. The architecture should therefore be **general engine as the correctness backbone +
   rules/heuristics as the fast path**, with **both certified**.
5. **Relationship to the rule route**: this route does not exclude a rule layer, but the rule layer
   must run after the mechanisms (see [GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md) §7.3), otherwise the
   project degrades back into enumeration.
6. **Relationship to the 0.27.3 defects**: items 1–2 of §7.3's priority list (residue resolution,
   expression-level cycle detection) are **defects**, not capabilities, and should be fixed first
   (item 1 has a measured net +4); they do not conflict with P0 here.

## 8. Conclusion

- **Feasibility: high.** The theory is a textbook-level decision procedure, and oCAS already owns
  the tower, the field operations, the RDE, rational integration, Trager factorisation, modular
  arithmetic, the oracle and the harness — **more than half the foundation**. What is missing is
  wiring in the algebraic case, relaxing the conservative refusals (dependent generators, rational
  RDE solutions), and adding a non-elementary layer.
- **Method: theorem-driven structural recursion + a certified-correctness discipline + an
  independent generality metric.** All three are required: recursion alone degenerates, certificates
  alone have no capability, and the corpus alone overfits.
- **Cost**: multi-year, fragment-by-fragment engineering, but every fragment delivers immediately
  measurable coverage plus a permanent correctness gain, and nothing conflicts with the LGPL
  position (all self-implemented; algorithms themselves are not copyrightable).
- **What not to do**: do not replace structure theorems with rule enumeration; do not hide
  complexity behind wall-clock budgets; do not widen the emission surface without certificates.
