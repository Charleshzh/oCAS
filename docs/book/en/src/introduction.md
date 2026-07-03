# Introduction

oCAS (open Computer Algebra System) is a modern, high-performance computer
algebra system written in Rust. It aims to match or exceed the core
performance of Symbolica and SageMath while remaining free and open under the
**LGPL-3.0-or-later** license.

---

## Why oCAS?

| Feature | oCAS | Symbolica | SymPy | SageMath |
|---|---|---|---|---|
| Language | Rust | Rust | Python | Python/Cython |
| License | LGPL-3.0+ | Proprietary / source-available | BSD | GPL |
| Native speed | ✅ | ✅ | ❌ | ⚠️ |
| Rust API | ✅ | ✅ | ❌ | ❌ |
| Python API | ✅ | ✅ | ✅ | ✅ |
| C/C++ API | ✅ | ❌ | ❌ | ❌ |
| No GPL contamination | ✅ | ❌ | ✅ | ❌ |

---

## Key Features

- **Layered Rust architecture** — 12 crates from the arena runtime up to
  language bindings, with strict downward dependencies.
- **Multiple coefficient domains** — arbitrary-precision integers, rationals,
  finite fields, real balls, and complex numbers.
- **Polynomial algebra** — dense/sparse multivariate polynomials, GCD,
  multivariate GCD, factorization (Hensel lifting), square-free factorization,
  Gröbner bases (Buchberger), root isolation.
- **Symbolic calculus** — differentiation, Taylor series, heuristic
  integration, expression substitution.
- **Linear algebra** — matrices with Bareiss determinant, rank, inverse,
  transpose, trace, and linear system solving.
- **Equation solvers** — linear systems (ℚ, ℤ), Diophantine equations,
  polynomial systems via Gröbner bases.
- **JIT evaluation** — Cranelift backend and SIMD-vectorized batch evaluation.
- **Rewrite & simplification** — pattern matching with wildcards, rule-based
  fixed-point simplification, optional egg e-graph equality saturation.
- **Tri-language bindings** — Rust, Python (PyO3), and C/C++ (cbindgen).
- **Correctness framework** — 82 automated cross-validation tests against
  SymPy, SageMath, and Symbolica across 16 modules.
- **Optional numerical backends** — GMP/MPFR/FLINT behind feature flags,
  isolated GPL backends in `ocas-gpl`.

---

## Project Status

oCAS is currently at version **0.11.0 (Beta)**. The core symbolic engine,
polynomial algebra, solvers, JIT evaluation, tri-language bindings, and
correctness comparison framework are feature-complete for a beta release. See
the [roadmap](https://github.com/charleshzh/ocas/blob/main/docs/planning/ROADMAP_EN.md)
for the path to stable 1.0.
