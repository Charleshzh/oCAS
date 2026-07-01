# oCAS — open Computer Algebra System / 开源计算机代数系统

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.89%2B-orange.svg)](https://www.rust-lang.org)
[![CI](https://github.com/charleshzh/ocas/actions/workflows/ci.yml/badge.svg)](https://github.com/charleshzh/ocas/actions/workflows/ci.yml)
[![docs.rs](https://docs.rs/ocas/badge.svg)](https://docs.rs/ocas)

**English**

**oCAS** is a modern, high-performance computer algebra system written in Rust. It is designed to match or exceed the core performance of Symbolica and SageMath while remaining free and open under the **LGPL-3.0-or-later** license.

> **Status**: Beta (0.10.0). The full symbolic engine, polynomial algebra, equation solvers, JIT/SIMD evaluation, and tri-language bindings (Rust, Python, C/C++) are feature-complete. See the [documentation](https://github.com/charleshzh/ocas/tree/main/ocas/docs/book) and [roadmap](ROADMAP.md).

**中文**

**oCAS** 是一个使用 Rust 编写的现代化高性能计算机代数系统。它的目标是在核心性能上达到或超越 Symbolica 与 SageMath，同时在 **LGPL-3.0-or-later** 许可证下保持自由与开放。

> **状态**：Beta (0.10.0)。完整的符号引擎、多项式代数、方程求解器、JIT/SIMD 求值与三语言绑定（Rust、Python、C/C++）已功能完备。详见[文档](https://github.com/charleshzh/ocas/tree/main/ocas/docs/book)与[路线图](ROADMAP.md)。

---

## Why oCAS? / 为什么选择 oCAS？

| Capability / 能力 | oCAS | Symbolica | SageMath |
|---|---|---|---|
| Core language / 核心语言 | Rust | Rust | Python/Cython |
| License / 许可证 | **LGPL-3.0+** | Proprietary / source-available | GPL-3.0 |
| Python API | ✅ PyO3 | ✅ | ✅ Native |
| C/C++ API | ✅ cbindgen | ✅ | ⚠️ Limited |
| Rust API | ✅ Native | ✅ | ❌ |
| High-performance backend / 高性能后端 | GMP / FLINT 3 / Arb / LinBox | Self-hosted + optional | FLINT / Singular / PARI |
| JIT code generation / JIT 代码生成 | Cranelift / LLVM | Custom | Limited |
| E-graph simplification / E-图简化 | ✅ `egg` (optional) | ❌ | Partial |

**English**

oCAS targets researchers, engineers, and developers who need a fast, embeddable, and license-clean symbolic computation engine.

**中文**

oCAS 面向需要快速、可嵌入且许可证清晰的符号计算引擎的研究人员、工程师和开发者。

---

## Key Features / 主要特性

**English**

- **Symbolic expressions**: variables, functions, rational numbers, arbitrary-precision integers, floats, and complex numbers.
- **Fast parsing & printing**: Mathematica-like and Python-like syntax, LaTeX output.
- **Simplification & rewriting**: normalization, hash consing, and equality saturation via `egg`.
- **Polynomial arithmetic**: dense and sparse multivariate polynomials, GCD, factorization, Gröbner bases, resultants.
- **Calculus**: symbolic differentiation, partial integration, Taylor/Laurent series.
- **Equation solving**: linear systems, polynomial systems, algebraic numbers.
- **Numerical evaluation**: arbitrary-precision floating point, interval arithmetic, numerical integration.
- **JIT compilation**: compile repeated expressions to native code with Cranelift or LLVM.
- **Multi-language bindings**: Rust, Python, and C/C++ from a single codebase.

**中文**

- **符号表达式**：变量、函数、有理数、任意精度整数、浮点数与复数。
- **快速解析与打印**：类 Mathematica 与 Python 语法，支持 LaTeX 输出。
- **化简与重写**：规范化、hash consing，以及基于 `egg` 的等式饱和。
- **多项式运算**：稠密与稀疏多元多项式、最大公因式、因式分解、Gröbner 基、结式。
- **微积分**：符号微分、部分积分、Taylor/Laurent 级数。
- **方程求解**：线性方程组、多项式方程组、代数数。
- **数值求值**：任意精度浮点、区间算术、数值积分。
- **JIT 编译**：将重复表达式编译为 Cranelift 或 LLVM 原生代码。
- **多语言绑定**：基于同一代码库提供 Rust、Python 与 C/C++ 接口。

---

## Architecture / 架构

```text
┌─────────────────────────────────────────────────────────────┐
│ Layer 5: Language Bindings / 多语言绑定                     │
│  Python (PyO3) │ C/C++ (cbindgen) │ Rust API                │
├─────────────────────────────────────────────────────────────┤
│ Layer 4: Applications / 应用层                              │
│  Parser │ Printer │ Solver │ Calculus │ Series              │
├─────────────────────────────────────────────────────────────┤
│ Layer 3: Symbolic Engine / 符号引擎                         │
│  Expression (Atom) │ Transformer │ Pattern Matcher │ Rewriting │ E-Graph │
├─────────────────────────────────────────────────────────────┤
│ Layer 2: Algebraic Kernel / 代数核                          │
│  Polynomial │ Matrix │ Domain │ Number Theory │ Factorize  │
├─────────────────────────────────────────────────────────────┤
│ Layer 1: Numeric Backend / 数值后端                         │
│  GMP/MPFR/MPC │ FLINT 3 / Arb │ NTL │ LinBox │ faer        │
├─────────────────────────────────────────────────────────────┤
│ Layer 0: Runtime / 运行时                                   │
│  Arena Allocator │ Error Handling │ Thread Pool │ FFI Glue  │
└─────────────────────────────────────────────────────────────┘
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for the full design.

完整设计请参阅 [ARCHITECTURE.md](ARCHITECTURE.md)。

---

## Quick Start / 快速开始

### Rust

```toml
[dependencies]
ocas = "0.4"
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

Rule-based simplification, pattern matching, and optional `egg` equality saturation are available through `ocas::prelude` (see `ocas-rewrite` examples).

---

## Building from Source / 从源码构建

### Requirements / 环境要求

**English**

- Rust 1.89 or later
- C compiler (MSVC, GCC, or Clang)
- GMP, MPFR, and FLINT 3 development libraries

**中文**

- Rust 1.89 或更高版本
- C 编译器（MSVC、GCC 或 Clang）
- GMP、MPFR 与 FLINT 3 开发库

On Ubuntu/Debian:

```bash
sudo apt-get install libgmp-dev libmpfr-dev libflint-dev
```

On macOS:

```bash
brew install gmp mpfr flint
```

On Windows, see [docs/build-windows.md](docs/build-windows.md).

Windows 用户请参考 [docs/build-windows.md](docs/build-windows.md)。

### Build / 构建

```bash
git clone https://github.com/charleshzh/ocas.git
cd ocas
cargo build --release
```

### Run Tests / 运行测试

```bash
cargo test --workspace --exclude ocas-py
```

---

## Performance / 性能

**English**

oCAS aims for competitive performance with Symbolica and SageMath on core workloads:

- Integer and rational arithmetic via **GMP**.
- Polynomial multiplication, GCD, and factorization via **FLINT 3**.
- Exact linear algebra via **faer** (default) and **LinBox** (optional).
- JIT evaluation via **Cranelift**.

Benchmarks are under active development. See [ocas-tests/benches](ocas-tests/benches).

**中文**

oCAS 致力于在核心工作负载上与 Symbolica 和 SageMath 竞争性能：

- 整数与有理数运算通过 **GMP**。
- 多项式乘法、最大公因式与因式分解通过 **FLINT 3**。
- 精确线性代数默认使用 **faer**，可选 **LinBox**。
- JIT 求值通过 **Cranelift**。

基准测试正在积极开发中，参见 [ocas-tests/benches](ocas-tests/benches)。

---

## Roadmap / 路线图

See [ROADMAP.md](ROADMAP.md) for the detailed versioned plan.

详细版本化计划请参阅 [ROADMAP.md](ROADMAP.md)。

---

## License / 许可证

**English**

This project is licensed under the **GNU Lesser General Public License v3.0 or later** (LGPL-3.0-or-later).

See [LICENSE](LICENSE) for the full text.

### Optional GPL Backends

Some optional backends (e.g., NTL, Singular, PARI/GP, SageMath interfaces) are licensed under GPL-compatible licenses. These are isolated in the `ocas-gpl` crate and disabled by default. Enabling them will produce a combined work under GPL terms.

**中文**

本项目采用 **GNU 宽通用公共许可证第 3 版或更新版本**（LGPL-3.0-or-later）。

完整文本参见 [LICENSE](LICENSE)。

### 可选 GPL 后端

某些可选后端（如 NTL、Singular、PARI/GP、SageMath 接口）采用 GPL 兼容许可证。它们被隔离在 `ocas-gpl` crate 中，默认禁用。启用后将产生受 GPL 条款约束的合并作品。

---

## Contributing / 贡献

**English**

We welcome contributions. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

**中文**

我们欢迎贡献。在提交 Pull Request 之前，请先阅读 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## Acknowledgments / 致谢

oCAS builds on the shoulders of many excellent open-source projects:

oCAS 建立在众多优秀开源项目的基础之上：

- [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/), [MPC](https://www.multiprecision.org/mpc/)
- [FLINT](https://flintlib.org/) and [Arb](https://arblib.org/) (now part of FLINT)
- [LinBox](https://linalg.org/) and [Givaro](https://github.com/linbox-team/givaro)
- [SymEngine](https://github.com/symengine/symengine)
- [faer](https://codeberg.org/sarah-quinones/faer) (primary development moved from [GitHub mirror](https://github.com/sarah-quinones/faer-rs))
- [egg](https://github.com/egraphs-good/egg)
- [Cranelift](https://github.com/bytecodealliance/wasmtime/tree/main/cranelift)
- [PyO3](https://github.com/PyO3/pyo3)
