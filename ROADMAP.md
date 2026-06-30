# oCAS Roadmap / oCAS 路线图

**English**

This document outlines the development roadmap of oCAS from pre-alpha experiments to a stable 1.0 release, with each 0.x version carrying concrete deliverables.

**中文**

本文档概述 oCAS 从 pre-alpha 实验到稳定 1.0 版本的开发路线图，每个 0.x 版本都包含具体的交付物。

---

## Legend / 图例

| Tag / 标签 | Meaning / 含义 |
|---|---|
| `API` | Public API surface / 公共 API |
| `ALG` | Algebraic algorithms / 代数算法 |
| `NUM` | Numerical backends / 数值后端 |
| `PERF` | Performance and optimization / 性能与优化 |
| `BIND` | Language bindings / 语言绑定 |
| `DOC` | Documentation and examples / 文档与示例 |
| `TEST` | Testing and quality / 测试与质量 |

---

## Phase 1: Pre-Alpha — Foundation / Pre-Alpha 阶段 — 基础

> **Goal / 目标**: Establish the workspace, runtime, and basic expression core. Prove that the layered architecture compiles and runs.
> 建立 workspace、运行时与基础表达式核心，证明分层架构可以编译并运行。

### 0.1.0 — Workspace & Runtime

**Target / 目标日期**: Month 1

**Deliverables / 交付物**:

- [x] Workspace structure with all 12 crates / 包含全部 12 个 crate 的 workspace 结构
- [x] CI pipeline: `cargo test`, `cargo clippy`, `cargo-deny`, formatting, Miri / CI 流水线
- [x] Unified error type `OcasError` / 统一错误类型
- [x] Arena / bump allocator with Miri-safe API / 通过 Miri 安全验证的 arena
- [x] Thread pool wrapper around `rayon` / 基于 `rayon` 的线程池包装
- [x] FFI glue conventions (minimal C ABI example) / FFI 胶水约定（最小 C ABI 示例）
- [x] GMP bindings (via `rug`) behind `gmp` feature / `gmp` feature 后的 GMP 绑定（基于 `rug`）
- [x] Initial benchmark harness / 初始基准测试框架

**Success Criteria / 成功标准**:

- `cargo build --workspace` succeeds on Linux/macOS/Windows (no-default-features on MSVC).
- Arena passes Miri and valgrind/ASan checks.
- GMP integer arithmetic is callable from Rust on supported platforms.

### 0.2.0 — Expression Tree Core

**Target / 目标日期**: Month 2

**Deliverables / 交付物**:

- [x] `ocas-atom` crate / `ocas-atom` crate
- [x] `Atom` tagged-union design / `Atom` 标签联合设计
- [x] Arena-backed AST with safe public API / 带安全公共 API 的 arena 后端 AST
- [x] Hash consing for common subexpressions / 公共子表达式 hash consing
- [x] Lexer using `logos` / 基于 `logos` 的词法分析器
- [x] Recursive-descent / Pratt parser / 递归下降 / Pratt 语法分析器
- [x] Printer: ASCII and compact forms / ASCII 与紧凑形式打印器
- [x] Normalizer: flatten `Add`/`Mul`, sort terms, merge coefficients / 规范化器

**Success Criteria / 成功标准**:

- `parse("x^2 + 2*x + 1")` produces the expected AST.
- `to_string(parse(s)) == s` for a broad set of expressions.
- Normalization is deterministic and property-tested.

---

## Phase 2: Alpha — Symbolic Engine / Alpha 阶段 — 符号引擎

> **Goal / 目标**: A usable Rust API for parsing, simplification, differentiation, and basic polynomial operations.
> 提供可用的 Rust API，支持解析、化简、微分与基础多项式运算。

### 0.3.0 — Domains & Polynomials

**Target / 目标日期**: Month 4

**Deliverables / 交付物**:

- [ ] `ocas-domain` crate / `ocas-domain` crate
- [ ] Domains: `Integer`, `Rational`, `FiniteField`, `RealBall`, `Complex` / 域实现
- [ ] Domain trait for generic algorithms / 泛型算法的 Domain trait
- [ ] `ocas-poly` crate / `ocas-poly` crate
- [ ] Dense univariate polynomial / 稠密单变量多项式
- [ ] Sparse multivariate polynomial / 稀疏多元多项式
- [ ] Addition, multiplication, division with remainder / 加减乘除与带余除法
- [ ] FLINT 3 integration behind `flint` feature / `flint` feature 后的 FLINT 3 集成

**Success Criteria / 成功标准**:

- Polynomial operations match SymPy outputs on regression suite.
- FLINT path produces identical results to pure-Rust fallback for supported operations.

### 0.4.0 — Pattern Matching & Rewriting

**Target / 目标日期**: Month 5

**Deliverables / 交付物**:

- [ ] Pattern matching engine with wildcards and conditions / 模式匹配引擎
- [ ] `Transformer` visitor API / `Transformer` 访问者 API
- [ ] Basic built-in rewrite rules / 基础内置重写规则
- [ ] `egg` integration for equality saturation / `egg` 等式饱和集成
- [ ] Rule-based simplifier / 基于规则的化简器

**Success Criteria / 成功标准**:

- Common identities (e.g., `x + x -> 2*x`, `x * 0 -> 0`) are applied automatically.
- E-graph can simplify `sin(x)^2 + cos(x)^2` to `1` under assumptions.

### 0.5.0 — Calculus Basics

**Target / 目标日期**: Month 6

**Deliverables / 交付物**:

- [ ] Symbolic differentiation / 符号微分
- [ ] Derivative table for elementary functions / 初等函数导数表
- [ ] Taylor series expansion / Taylor 级数展开
- [ ] Partial integration with heuristic table / 基于启发式表的部分积分
- [ ] `ocas-calc` crate initial release / `ocas-calc` 初始版本

**Success Criteria / 成功标准**:

- Differentiation passes a comprehensive test suite.
- Integration succeeds on standard calculus problems.

### 0.6.0 — First Rust API Release Candidate

**Target / 目标日期**: Month 7

**Deliverables / 交付物**:

- [ ] Stable `ocas` prelude / 稳定的 `ocas` prelude
- [ ] Rustdoc examples for all public APIs / 所有公共 API 的 rustdoc 示例
- [ ] Property tests with `proptest` / `proptest` 属性测试
- [ ] Initial benchmark suite / 初始基准测试套件
- [ ] crates.io publish (alpha) / crates.io 发布（alpha）

**Success Criteria / 成功标准**:

- `cargo test --workspace` passes.
- Benchmarks demonstrate parity with SymPy on basic polynomial operations.

---

## Phase 3: Beta — Solvers, JIT, Bindings / Beta 阶段 — 求解器、JIT、绑定

> **Goal / 目标**: Multi-language availability and performance. Core algebra is feature-complete for a CAS beta.
> 实现多语言可用性与性能，核心代数功能达到 CAS beta 标准。

### 0.7.0 — Equation Solvers

**Target / 目标日期**: Month 9

**Deliverables / 交付物**:

- [ ] Linear system solver (`faer` / `LinBox`) / 线性方程组求解器
- [ ] Polynomial system solver (Gröbner + root isolation) / 多项式方程组求解器
- [ ] Single-variable root finding via Arb / 基于 Arb 的单变量求根
- [ ] Diophantine solver basics / 丢番图方程基础
- [ ] Assumptions / domain system / 假设/域系统

**Success Criteria / 成功标准**:

- Linear and polynomial solvers produce correct results verified against SageMath.

### 0.8.0 — Evaluation & JIT

**Target / 目标日期**: Month 11

**Deliverables / 交付物**:

- [ ] Tree interpreter for scalar and vector evaluation / 标量与向量求值解释器
- [ ] AST-to-instruction compiler / AST 到指令序列编译器
- [ ] Function registry for user-defined functions / 用户自定义函数注册表
- [ ] Cranelift JIT backend / Cranelift JIT 后端
- [ ] Optional LLVM/Inkwell backend / 可选 LLVM/Inkwell 后端
- [ ] SIMD vectorized evaluation / SIMD 向量化求值

**Success Criteria / 成功标准**:

- JIT evaluates repeated expressions at least 10x faster than interpreter.
- SIMD path works for dense polynomial evaluation.

### 0.9.0 — Python & C/C++ Bindings

**Target / 目标日期**: Month 13

**Deliverables / 交付物**:

- [ ] `ocas-py` crate with PyO3 / 基于 PyO3 的 `ocas-py`
- [ ] Python classes: `Expression`, `Polynomial`, `Matrix`, `Domain` / Python 类
- [ ] Maturin wheel build for Linux/macOS/Windows / Maturin 轮子构建
- [ ] `ocas-c` crate with cbindgen / 基于 cbindgen 的 `ocas-c`
- [ ] Stable C API for expression lifecycle / 稳定的表达式生命周期 C API
- [ ] C++ RAII wrapper / C++ RAII 包装

**Success Criteria / 成功标准**:

- `pip install ocas` works on supported platforms.
- C example compiles and runs against the shared library.
- No memory leaks in binding tests.

### 0.10.0 — Beta Release

**Target / 目标日期**: Month 14

**Deliverables / 交付物**:

- [ ] Feature freeze for 1.0 / 1.0 功能冻结
- [ ] Comprehensive documentation site / 综合文档站点
- [ ] Performance comparison with Symbolica and SageMath / 与 Symbolica/SageMath 的性能对比
- [ ] Community feedback integration / 社区反馈整合
- [ ] Bug-fix only period / 仅修复 bug 阶段

**Success Criteria / 成功标准**:

- All public APIs documented.
- CI green on all supported platforms.

---

## Phase 4: Stable 1.0 / 稳定 1.0

> **Goal / 目标**: A production-ready CAS library with stable APIs and broad backend support.
> 发布 API 稳定、后端支持广泛的成熟 CAS 库。

### 1.0.0 — Stable Release

**Target / 目标日期**: Month 16

**Deliverables / 交付物**:

- [ ] Stable semantic versioning guarantee / 稳定语义化版本保证
- [ ] Full Rust, Python, and C/C++ API coverage / 完整的 Rust、Python 与 C/C++ API
- [ ] Comprehensive test suite (>80% line coverage) / 综合测试套件（行覆盖率 >80%）
- [ ] Published benchmarks / 已发布基准测试
- [ ] Migration guide from Symbolica/SymPy / 从 Symbolica/SymPy 迁移指南
- [ ] Signed release artifacts / 签名发布产物

**Success Criteria / 成功标准**:

- No breaking API changes planned for 1.x.
- Performance parity or better with Symbolica on core benchmarks.

---

## Post-1.0 / 1.0 之后

**English**

After 1.0, development will focus on:

- Advanced symbolic integration (Risch algorithm)
- Differential equation solvers
- Optional GPL backends (`ocas-gpl`)
- GPU acceleration (CUDA / HIP / Vulkan compute)
- Domain-specific toolkits (physics, robotics, machine learning)

**中文**

1.0 之后，开发重点将转向：

- 高级符号积分（Risch 算法）
- 微分方程求解器
- 可选 GPL 后端（`ocas-gpl`）
- GPU 加速（CUDA / HIP / Vulkan compute）
- 领域专用工具包（物理、机器人、机器学习）

---

## Milestones / 里程碑

| Version / 版本 | Phase / 阶段 | Target / 目标日期 | Key Deliverable / 关键交付物 |
|---|---|---|---|
| 0.1.0 | Pre-Alpha | Month 1 | Workspace + runtime / 工作空间 + 运行时 |
| 0.2.0 | Pre-Alpha | Month 2 | Expression core / 表达式核心 |
| 0.3.0 | Alpha | Month 4 | Domains & polynomials / 域与多项式 |
| 0.4.0 | Alpha | Month 5 | Pattern matching & rewriting / 模式匹配与重写 |
| 0.5.0 | Alpha | Month 6 | Calculus basics / 微积分基础 |
| 0.6.0 | Alpha | Month 7 | Rust API RC / Rust API 候选版 |
| 0.7.0 | Beta | Month 9 | Equation solvers / 方程求解器 |
| 0.8.0 | Beta | Month 11 | JIT & evaluation / JIT 与求值 |
| 0.9.0 | Beta | Month 13 | Python & C/C++ bindings / Python 与 C/C++ 绑定 |
| 0.10.0 | Beta | Month 14 | Feature freeze / 功能冻结 |
| 1.0.0 | Stable | Month 16 | Stable release / 稳定版发布 |

---

## How to Read This Roadmap / 如何阅读本路线图

**English**

- Each version represents a **potentially publishable** increment.
- Dates are approximate and depend on contributor availability.
- Features may shift between versions based on user feedback and technical discoveries.

**中文**

- 每个版本代表一个**可发布**的增量。
- 日期为预估值，取决于贡献者可用时间。
- 功能可能根据用户反馈与技术发现在不同版本间调整。

---

## Contributing to the Roadmap / 参与路线图

**English**

If you want to work on a specific version or feature, please open a GitHub issue and we will assign a tracking issue to you.

**中文**

如果你想参与某个特定版本或功能，请创建 GitHub issue，我们会为你分配跟踪 issue。
