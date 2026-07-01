# Changelog / 变更日志

All notable changes to the oCAS project will be documented in this file.

oCAS 项目的所有重大变更都将记录在此文件中。

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---
## [0.5.0] - 2026-07-01

### Added / 新增

- `ocas-calc` crate initial release / `ocas-calc` crate 初始版本
- Symbolic differentiation (`diff`) / 符号微分 (`diff`)
- Derivative table for elementary functions (`sin`, `cos`, `exp`, `log`, `sqrt`, `tan`, `sec`) / 初等函数导数表
- Chain rule, product rule, and generalized power rule / 链式法则、乘积法则与广义幂法则
- Taylor series expansion (`taylor`) around a point / 在某点处的 Taylor 级数展开 (`taylor`)
- Heuristic integration (`integrate`) with table lookup and linear substitution / 基于查表与线性替换的启发式积分 (`integrate`)
- Linear substitution support for integrals / 积分中的线性替换支持
- `Derivative(expr, var)` and `Integral(expr, var)` unevaluated forms / 未求出的 `Derivative(expr, var)` 与 `Integral(expr, var)` 形式
- Re-export `diff`, `integrate`, and `taylor` from the top-level `ocas` crate / 在顶层 `ocas` crate 中重新导出 `diff`、`integrate`、`taylor`
- End-to-end calculus integration tests / 端到端微积分集成测试

### Changed / 变更

- `normalize` now removes `+0` and `*1` identity terms, absorbs `*0` into `0`, and preserves argument order for `Derivative` / `Integral` / `normalize` 现在会移除 `+0`、`*1` 单位元，将 `*0` 吸收为 `0`，并对 `Derivative` / `Integral` 保持参数顺序
- Workspace version bumped to `0.5.0` / 工作区版本提升至 `0.5.0`

---

## [0.4.0] - 2026-07-01

### Added / 新增

- `ocas-rewrite` crate / `ocas-rewrite` crate
- Pattern AST with wildcards and conditional matching / 带通配符与条件匹配的模式 AST
- AC (associative-commutative) matching for `Add`/`Mul` / `Add`/`Mul` 的 AC 匹配
- `Transformer` visitor API for expression traversal / 用于表达式遍历的 `Transformer` 访问者 API
- `Rule` type with closure-based replacers and conditions / 支持闭包替换器/条件的 `Rule` 类型
- Built-in algebraic rewrite rules (`x + x -> 2*x`, `x * 0 -> 0`, etc.) / 内置代数重写规则
- Rule-based `simplify` engine / 基于规则的 `simplify` 化简引擎
- Optional `egg` feature for equality saturation / 可选的 `egg` 等式饱和特性
- `sin(x)^2 + cos(x)^2 -> 1` e-graph simplification test / `sin(x)^2 + cos(x)^2 -> 1` e-graph 化简测试
- Re-export `Pattern`, `Rule`, `simplify`, and `transform` from the top-level `ocas` crate / 在顶层 `ocas` crate 中重新导出 `Pattern`、`Rule`、`simplify` 与 `transform`

### Notes / 说明

- The `egg` feature is optional and disabled by default. It enables equality
  saturation as an additional simplification backend on supported platforms.
- `egg` 特性为可选，默认未启用。它在支持的平台上作为额外的化简后端提供等式饱和。

---

## [0.3.0] - 2026-06-30

### Added / 新增

- `ocas-domain` crate / `ocas-domain` crate
- `Domain` and `EuclideanDomain` traits / `Domain` 与 `EuclideanDomain` trait
- Domains: `Integer`, `Rational`, `FiniteField` / 域实现
- Domains: `RealBall` (with optional `mpfr` backend) and `Complex` / 域实现：可选 `mpfr` 后端的 `RealBall` 与 `Complex`
- Optional `gmp` feature using `rug` for GMP-backed `Integer` and `Rational` / 可选 `gmp` 特性：基于 `rug` 的 GMP 后端 `Integer` 与 `Rational`
- `ocas-poly` crate / `ocas-poly` crate
- Dense univariate polynomial with `div_rem` / 带 `div_rem` 的稠密单变量多项式
- Sparse multivariate polynomial with `Lex` and `Grevlex` orderings / 支持 `Lex` 与 `Grevlex` 序的稀疏多元多项式
- Experimental FLINT 3 backend for integer polynomials behind `flint` feature / `flint` feature 后用于整数多项式的实验性 FLINT 3 后端
- Re-export `RealBall`, `Complex`, and `SparseMultivariatePolynomial` in `ocas::prelude` / 在 `ocas::prelude` 中重新导出 `RealBall`、`Complex` 与 `SparseMultivariatePolynomial`

### Notes / 说明

- The `flint` feature is experimental and requires system FLINT. It is not
  yet supported on Windows because `flint3-sys` depends on POSIX-only types.
  Use Linux, macOS, or WSL for FLINT-backed tests.

- `flint` 特性为实验性，需要系统 FLINT。由于 `flint3-sys` 依赖仅 POSIX 的
  类型，目前尚不支持 Windows。请在 Linux、macOS 或 WSL 下运行 FLINT 后端测试。

---

## [0.2.0] - 2026-06-30

### Added / 新增

- `ocas-atom` crate / `ocas-atom` crate
- `Atom` tagged-union design / `Atom` 标签联合设计
- Arena-backed AST with safe public API / 带安全公共 API 的 arena 后端 AST
- Hash consing for common subexpressions / 公共子表达式 hash consing
- Lexer using `logos` / 基于 `logos` 的词法分析器
- Recursive-descent / Pratt parser / 递归下降 / Pratt 语法分析器
- Printer: ASCII and compact forms / ASCII 与紧凑形式打印器
- Normalizer: flatten `Add`/`Mul`, sort terms, merge coefficients / 规范化器

### Notes / 说明

- This release is **Pre-Alpha** and focuses on the expression tree core.
  Algebraic-domain crates (`ocas-domain`, `ocas-poly`, etc.) are still
  placeholders.

- 本版本为 **Pre-Alpha**，仅聚焦表达式树核心。代数域相关 crate
  （`ocas-domain`、`ocas-poly` 等）仍为占位。

---

## [0.1.0] - 2026-06-30

### Added / 新增

- Workspace with all 12 crates / 包含全部 12 个 crate 的工作空间
- Cross-platform CI pipeline (fmt, clippy, test, backend test, bindings, Miri) /
  跨平台 CI 流水线（格式化、Clippy、测试、后端测试、绑定构建、Miri）
- `OcasError` unified error type based on `thiserror` /
  基于 `thiserror` 的统一错误类型 `OcasError`
- Bump allocator `Arena` with `allocate_with` API and Miri-validated safety /
  带 `allocate_with` API 的 bump allocator `Arena`，通过 Miri 安全验证
- `OwnedExpr<T>` for self-contained arena-backed expressions /
  自包含的 arena 后端表达式类型 `OwnedExpr<T>`
- `ThreadPool` wrapper around `rayon` with `OcasError` propagation /
  带 `OcasError` 传播的 `rayon` 线程池包装 `ThreadPool`
- Optional `GmpInteger` backend behind the `gmp` feature /
  `gmp` feature 后的可选 `GmpInteger` 后端
- Minimal C ABI example (`ocas_version`, `ocas_arena_new/free`, error handling) /
  最小 C ABI 示例（版本、arena 生命周期、错误处理）
- C example program in `ocas-c/examples/basic.c` /
  `ocas-c/examples/basic.c` 中的 C 示例程序
- Real micro-benchmarks for arena allocation and GMP integer arithmetic /
  arena 分配与 GMP 整数运算的真实微基准测试

### Notes / 说明

- This release is **Pre-Alpha** and focuses on runtime foundations only.
  Symbolic computation crates (`ocas-atom`, `ocas-poly`, etc.) are still
  placeholders.

- 本版本为 **Pre-Alpha**，仅聚焦运行时基础。符号计算 crate
  （`ocas-atom`、`ocas-poly` 等）仍为占位。

- The `gmp` feature is not supported on Windows MSVC because `rug` cannot
  build GMP in that environment. Use MSYS2/MINGW64 or Linux/macOS instead.

- `gmp` feature 在 Windows MSVC 上不受支持，因为 `rug` 无法在该环境下构建
  GMP。请改用 MSYS2/MINGW64 或 Linux/macOS。
