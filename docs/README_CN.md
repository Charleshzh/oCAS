# oCAS — 开源计算机代数系统

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.97%2B-orange.svg)](https://www.rust-lang.org)
[![CI](https://github.com/charleshzh/ocas/actions/workflows/ci.yml/badge.svg)](https://github.com/charleshzh/ocas/actions/workflows/ci.yml)
[![docs.rs](https://docs.rs/ocas/badge.svg)](https://docs.rs/ocas)
[![PyPI](https://img.shields.io/pypi/v/ocas.svg)](https://pypi.org/project/ocas/)

> **语言**: [English](../README.md) · [中文](README_CN.md)

**oCAS** 是一个使用 Rust 编写的现代化高性能计算机代数系统。它的目标是在核心
性能上达到或超越 Symbolica 与 SageMath，同时在 **LGPL-3.0-or-later** 许可证下
保持自由与开放。

> **状态**：Beta (0.13.2)。完整的符号引擎、多项式代数、方程求解器、JIT/SIMD
> 求值与三语言绑定（Rust、Python、C/C++）已功能完备。`pip install ocas`
> 已在 PyPI 上线。详见
> [文档](https://charleshzh.github.io/ocas/latest/zh/)
>（[English](https://charleshzh.github.io/ocas/latest/en/)）、
> [Rust API 文档](https://docs.rs/ocas)与
> [路线图](planning/ROADMAP_CN.md)。
>
> **docs.rs 说明**：托管的 Rust API 文档仅使用可移植特性构建（不包含系统
> GMP/MPFR/FLINT 后端）。如需查看包含后端特性的完整 API，请本地构建：
> `cargo doc -p ocas --features gmp,mpfr,flint --no-deps`。

---

## 为什么选择 oCAS？

| 能力 | oCAS | Symbolica | SageMath |
|---|---|---|---|
| 核心语言 | Rust | Rust | Python/Cython |
| 许可证 | **LGPL-3.0+** | 专有 / 源码可用 | GPL-3.0 |
| Python API | ✅ PyO3 | ✅ | ✅ 原生 |
| C/C++ API | ✅ cbindgen | ✅ | ⚠️ 有限 |
| Rust API | ✅ 原生 | ✅ | ❌ |
| 高性能后端 | GMP / FLINT 3 / Arb / LinBox | 自托管 + 可选 | FLINT / Singular / PARI |
| JIT 代码生成 | Cranelift / LLVM | 自研 | 有限 |
| E-图化简 | ✅ `egg`（可选） | ❌ | 部分 |

oCAS 面向需要快速、可嵌入且许可证清晰的符号计算引擎的研究人员、工程师与
开发者。

---

## 主要特性

- **符号表达式**：变量、函数、有理数、任意精度整数、浮点数与复数。
- **快速解析与打印**：类 Mathematica 与 Python 语法，支持 LaTeX 输出。
- **化简与重写**：规范化、hash consing，以及基于 `egg` 的等式饱和。
- **多项式运算**：稠密与稀疏多元多项式、最大公因式、因式分解、Gröbner 基、
  结式。
- **微积分**：符号微分、部分积分、Taylor/Laurent 级数。
- **方程求解**：线性方程组、多项式方程组、代数数。
- **数值求值**：任意精度浮点、区间算术、数值积分。
- **JIT 编译**：将重复表达式编译为 Cranelift 或 LLVM 原生代码。
- **多语言绑定**：基于同一代码库提供 Rust、Python 与 C/C++ 接口。
- **优化整数算术**：SOO（小对象优化）将 fit i64 的值栈分配，大值回退 GMP。
- **模方法多变量 GCD**：通过模素数约化与重构计算 $\mathbb{Z}[x,y]$ 上的 GCD。

---

## 架构

```text
┌─────────────────────────────────────────────────────────────┐
│ 第 5 层：多语言绑定                                          │
│  Python (PyO3) │ C/C++ (cbindgen) │ Rust API                │
├─────────────────────────────────────────────────────────────┤
│ 第 4 层：应用层                                              │
│  Parser │ Printer │ Solver │ Calculus │ Series              │
├─────────────────────────────────────────────────────────────┤
│ 第 3 层：符号引擎                                            │
│  Expression (Atom) │ Transformer │ Pattern Matcher │ Rewriting │ E-Graph │
├─────────────────────────────────────────────────────────────┤
│ 第 2 层：代数核                                              │
│  Polynomial │ Matrix │ Domain │ Number Theory │ Factorize  │
├─────────────────────────────────────────────────────────────┤
│ 第 1 层：数值后端                                            │
│  GMP/MPFR/MPC │ FLINT 3 / Arb │ NTL │ LinBox │ faer        │
├─────────────────────────────────────────────────────────────┤
│ 第 0 层：运行时                                              │
│  Arena Allocator │ Error Handling │ Thread Pool │ FFI Glue  │
└─────────────────────────────────────────────────────────────┘
```

完整设计请参阅 [ARCHITECTURE_CN.md](planning/ARCHITECTURE_CN.md)。

---

## 快速开始

### Rust

```toml
[dependencies]
ocas = "0.13"
```

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

fn main() -> Result<(), ParseError> {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    let expr = parse(&ctx, "x^2 + 2*x + 1")?;
    println!("{}", expr); // 1 + 2*x + x^2

    Ok(())
}
```

基于规则的化简、模式匹配与可选的 `egg` 等式饱和可通过 `ocas::prelude` 使用
（参见 `ocas-rewrite` 示例）。

---

## 从源码构建

### 环境要求

- Rust 1.97 或更高版本
- C 编译器（MSVC、GCC 或 Clang）
- GMP、MPFR 与 FLINT 3 开发库

在 Ubuntu/Debian 上：

```bash
sudo apt-get install libgmp-dev libmpfr-dev libflint-dev
```

在 macOS 上：

```bash
brew install gmp mpfr flint
```

Windows 用户请参考 [build-windows.md](build-windows.md)。

### 构建

```bash
git clone https://github.com/charleshzh/ocas.git
cd ocas
cargo build --release
```

### 运行测试

```bash
cargo test --workspace --exclude ocas-py
```

---

## 性能

oCAS 致力于在核心工作负载上与 Symbolica 和 SageMath 竞争性能：

- 整数与有理数运算通过 **GMP**。
- 多项式乘法、最大公因式与因式分解通过 **FLINT 3**。
- 精确线性代数默认使用 **faer**，可选 **LinBox**。
- JIT 求值通过 **Cranelift**。

基准测试正在积极开发中，参见 [ocas-tests/benches](../ocas-tests/benches)。与
Symbolica、SageMath、SymPy 的差距分析见 [GAP_ANALYSIS_CN.md](planning/GAP_ANALYSIS_CN.md)。

---

## 路线图

版本化计划见 [ROADMAP_CN.md](planning/ROADMAP_CN.md)，细粒度逐版本演进计划（Beta → 1.0）
见 [EVOLUTION_PLAN_CN.md](planning/EVOLUTION_PLAN_CN.md)。

---

## 许可证

本项目采用 **GNU 宽通用公共许可证第 3 版或更新版本**（LGPL-3.0-or-later）。
完整文本参见 [LICENSE](../LICENSE)。

### 可选 GPL 后端

某些可选后端（如 NTL、Singular、PARI/GP、SageMath 接口）采用 GPL 兼容许可证。
它们被隔离在 `ocas-gpl` crate 中，默认禁用。启用后将产生受 GPL 条款约束的合并
作品。

---

## 贡献

我们欢迎贡献。在提交 Pull Request 之前，请先阅读
[CONTRIBUTING.md](book/zh/src/contributing.md)。

---

## 致谢

oCAS 建立在众多优秀开源项目的基础之上：

- [GMP](https://gmplib.org/)、[MPFR](https://www.mpfr.org/)、[MPC](https://www.multiprecision.org/mpc/)
- [FLINT](https://flintlib.org/) 与 [Arb](https://arblib.org/)（现属 FLINT）
- [LinBox](https://linalg.org/) 与 [Givaro](https://github.com/linbox-team/givaro)
- [SymEngine](https://github.com/symengine/symengine)
- [faer](https://codeberg.org/sarah-quinones/faer)（主开发已从 [GitHub 镜像](https://github.com/sarah-quinones/faer-rs)迁移）
- [egg](https://github.com/egraphs-good/egg)
- [Cranelift](https://github.com/bytecodealliance/wasmtime/tree/main/cranelift)
- [PyO3](https://github.com/PyO3/pyo3)
