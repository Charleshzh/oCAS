# Introduction / 简介

**English**

oCAS (open Computer Algebra System) is a modern, high-performance computer
algebra system written in Rust. It aims to match or exceed the core
performance of Symbolica and SageMath while remaining free and open under the
**LGPL-3.0-or-later** license.

oCAS（开源计算机代数系统）是一个使用 Rust 编写的现代化高性能计算机代数系统，目标是在核心性能上达到或超越 Symbolica 与 SageMath，同时在 LGPL-3.0-or-later 许可证下保持自由开源。

---

## Why oCAS? / 为什么选择 oCAS？

**English**

| Feature | oCAS | Symbolica | SymPy | SageMath |
|---|---|---|---|---|
| Language | Rust | Rust | Python | Python/Cython |
| License | LGPL-3.0+ | Proprietary / source-available | BSD | GPL |
| Native speed | ✅ | ✅ | ❌ | ⚠️ |
| Rust API | ✅ | ✅ | ❌ | ❌ |
| Python API | ✅ | ✅ | ✅ | ✅ |
| C/C++ API | ✅ | ❌ | ❌ | ❌ |
| No GPL contamination | ✅ | ❌ | ✅ | ❌ |

**中文**

| 特性 | oCAS | Symbolica | SymPy | SageMath |
|---|---|---|---|---|
| 语言 | Rust | Rust | Python | Python/Cython |
| 许可证 | LGPL-3.0+ | 专有/源码可用 | BSD | GPL |
| 原生性能 | ✅ | ✅ | ❌ | ⚠️ |
| Rust API | ✅ | ✅ | ❌ | ❌ |
| Python API | ✅ | ✅ | ✅ | ✅ |
| C/C++ API | ✅ | ❌ | ❌ | ❌ |
| 无 GPL 污染 | ✅ | ❌ | ✅ | ❌ |

---

## Key Features / 关键特性

**English**

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

**中文**

- **分层 Rust 架构** —— 从 arena 运行时到语言绑定的 12 个 crate，依赖严格向下。
- **多种系数域** —— 任意精度整数、有理数、有限域、实数球、复数。
- **多项式代数** —— 稠密/稀疏多元多项式、GCD、无平方因子分解、Gröbner 基（Buchberger）。
- **符号微积分** —— 微分、Taylor 级数、启发式积分。
- **JIT 求值** —— Cranelift 后端与 SIMD 向量化求值。
- **三语言绑定** —— Rust、Python（PyO3）、C/C++（cbindgen）。
- **可选数值后端** —— GMP/MPFR/FLINT 隐藏在 feature flag 后，GPL 后端隔离在 `ocas-gpl`。

---

## Project Status / 项目状态

**English**

oCAS is currently at version **0.10.0 (Beta)**. The core symbolic engine,
polynomial algebra, solvers, JIT evaluation, and tri-language bindings are
feature-complete for a beta release. See the
[roadmap](https://github.com/charleshzh/ocas/blob/main/ROADMAP.md) for the
path to stable 1.0.

**中文**

oCAS 当前版本为 **0.10.0（Beta）**。核心符号引擎、多项式代数、求解器、JIT 求值与三语言绑定已达到 Beta 功能完备。通往稳定 1.0 的路线见
[路线图](https://github.com/charleshzh/ocas/blob/main/ROADMAP.md)。
