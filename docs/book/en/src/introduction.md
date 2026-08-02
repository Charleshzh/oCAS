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
| License | LGPL-3.0+ | source-available / commercial | BSD | GPL |
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
  multivariate GCD, factorization (Hensel lifting), algebraic number fields
  (Trager), Gröbner bases (Buchberger, F4, F5), root isolation.
- **Symbolic calculus** — differentiation, Taylor series, Risch algorithm,
  heuristic integration, expression substitution.
- **ODE solvers** — first-order (separable, linear, Bernoulli, exact,
  integrating factors), second-order (constant coefficient, Cauchy–Euler,
  reduction of order, undetermined coefficients, variation of parameters),
  series solutions (power series, Frobenius), Laplace IVP, 2×2 systems.
- **Linear algebra** — matrices with Bareiss determinant, rank, inverse,
  transpose, trace, and linear system solving.
- **Numerical integration** — adaptive Monte Carlo (Vegas) for
  high-dimensional definite integrals.
- **Automatic differentiation** — forward-mode via hyper-dual numbers
  (`HyperDual`), first-order and higher-order derivatives.
- **Tensor algebra** — index slots, explicit contraction, symmetrization.
- **Equation solvers** — linear systems (ℚ, ℤ), Diophantine equations,
  polynomial systems via Gröbner bases.
- **JIT evaluation** — Cranelift backend and SIMD-vectorized batch evaluation.
- **Rewrite & simplification** — pattern matching with wildcards, rule-based
  fixed-point simplification, fuel-bounded simplification, optional egg
  e-graph equality saturation.
- **Tri-language bindings** — Rust, Python (PyO3), and C/C++ (cbindgen).
- **Correctness framework** — cross-validation tests against SymPy,
  SageMath, and Symbolica across multiple mathematical modules.
- **Optional numerical backends** — GMP/MPFR/FLINT behind feature flags,
  isolated GPL backends in `ocas-gpl`.

---

## Project Status

oCAS is currently at version **0.20.1 (Beta)**. The core symbolic engine,
polynomial algebra, ODE solvers, JIT evaluation, tri-language bindings,
and correctness comparison framework are feature-complete for a beta release. See
the [roadmap](https://github.com/charleshzh/ocas/blob/main/docs/planning/ROADMAP_EN.md)
for the path to stable 1.0.
