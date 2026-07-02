# 简介

oCAS（开源计算机代数系统）是一个使用 Rust 编写的现代化高性能计算机代数系统，目标是在核心性能上达到或超越 Symbolica 与 SageMath，同时在 LGPL-3.0-or-later 许可证下保持自由开源。

---

## 为什么选择 oCAS？

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

## 关键特性

- **分层 Rust 架构** —— 从 arena 运行时到语言绑定的 12 个 crate，依赖严格向下。
- **多种系数域** —— 任意精度整数、有理数、有限域、实数球、复数。
- **多项式代数** —— 稠密/稀疏多元多项式、GCD、无平方因子分解、Gröbner 基（Buchberger）。
- **符号微积分** —— 微分、Taylor 级数、启发式积分。
- **JIT 求值** —— Cranelift 后端与 SIMD 向量化求值。
- **三语言绑定** —— Rust、Python（PyO3）、C/C++（cbindgen）。
- **可选数值后端** —— GMP/MPFR/FLINT 隐藏在 feature flag 后，GPL 后端隔离在 `ocas-gpl`。

---

## 项目状态

oCAS 当前版本为 **0.10.0（Beta）**。核心符号引擎、多项式代数、求解器、JIT 求值与三语言绑定已达到 Beta 功能完备。通往稳定 1.0 的路线见
[路线图](https://github.com/charleshzh/ocas/blob/main/docs/planning/ROADMAP_CN.md)。
