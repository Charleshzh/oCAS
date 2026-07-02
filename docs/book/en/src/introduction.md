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
  square-free factorization, Gröbner bases (Buchberger).
- **Symbolic calculus** — differentiation, Taylor series, heuristic
  integration.
- **JIT evaluation** — Cranelift backend and SIMD-vectorized evaluation.
- **Tri-language bindings** — Rust, Python (PyO3), and C/C++ (cbindgen).
- **Optional numerical backends** — GMP/MPFR/FLINT behind feature flags,
  isolated GPL backends in `ocas-gpl`.

---

## Project Status

oCAS is currently at version **0.10.0 (Beta)**. The core symbolic engine,
polynomial algebra, solvers, JIT evaluation, and tri-language bindings are
feature-complete for a beta release. See the
[roadmap](https://github.com/charleshzh/ocas/blob/main/ROADMAP_EN.md) for the
path to stable 1.0.
