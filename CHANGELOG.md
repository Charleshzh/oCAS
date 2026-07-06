# Changelog / 变更日志

All notable changes to the oCAS project will be documented in this file.

oCAS 项目的所有重大变更都将记录在此文件中。

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.13.0] - 2026-07-06

### Added / 新增

- **F4 Gröbner basis algorithm** (`ocas-poly`): Matrix-based F4 algorithm from
  Faugère (1999) replacing sequential Buchberger for batched S-polynomial
  reductions. Entry point: `ocas_poly::groebner::f4::f4()`. Supports any
  `Domain` with specialized ℤ_p fast path using `i64` arithmetic /
  **F4 Gröbner 基算法**：基于 Faugère 1999 论文的矩阵化 F4 算法，批量处理
  S-多项式约化。入口：`ocas_poly::groebner::f4::f4()`，支持任意域，ℤ_p 专用
  快速路径
- **Gebauer-Moeller pair filtering** (`ocas-poly`): `CriticalPair` struct with
  precomputed `lcm`/`degree`. `update_pairs()` implements first criterion
  (coprime skip), second criterion (lcm minimality), and redundant pair cleanup /
  **Gebauer-Moeller 临界对筛选**：预计算 lcm/degree，第一判据（互素跳过）、
  第二判据（lcm 最小性）、冗余对清理
- **Simplification cache** (`ocas-poly`): `SimpCache` per-basis-element cache
  for `basis[i].mul_monomial(diff)` products, avoiding redundant computation
  in symbolic preprocessing / **简化缓存**：每个基元素的乘积缓存，避免符号
  预处理中的重复计算
- **`Grlex` monomial ordering** (`ocas-poly`): Graded lexicographic order
  alongside existing `Lex` and `Grevlex` / **`Grlex` 单项式序**：分次字典序
- **`Domain` trait extensions** (`ocas-domain`): `mul_assign()` and
  `sub_mul_assign()` default methods for in-place arithmetic /
  **`Domain` trait 扩展**：原地乘法和减乘融合操作
- **`FiniteField` utilities** (`ocas-domain`): `prime_u64()`, `to_i64()`,
  `from_i64()` for ℤ_p fast path conversion / **`FiniteField` 工具方法**

### Performance / 性能

- F4 cyclic-3 ℚ: **147 µs** (was 308 µs, **-52%**), now faster than Buchberger
- F4 cyclic-4 ℚ: **2.13 ms** (was 3.99 ms, **-47%**)
- F4 cyclic-3 ℤ₁₃: **276 µs** (was 582 µs, **-53%**)
- F4 cyclic-4 ℤ₁₃: **2.80 ms** (was 6.19 ms, **-55%**)
- F4 cyclic-3 ℤ₁₀₁: **270 µs** (was 517 µs, **-48%**)
- F4 cyclic-4 ℤ₁₀₁: **2.89 ms** (was 4.87 ms, **-41%**)

### Fixed / 修复

- **`minimize()` bug** (`ocas-poly`): `monomial_divides` arguments were
  swapped (`&lms[j], &lms[i]` instead of `&lms[i], &lms[j]`), preventing
  redundant LM removal. This was the root cause of incorrect Gröbner basis
  output / **`minimize()` bug**：`monomial_divides` 参数顺序错误，导致冗余
  首项未被移除
- **`auto_reduce()` direction** (`ocas-poly`): Now processes elements in
  ascending LM order and reduces only by smaller-LM elements, matching the
  standard reduced Gröbner basis algorithm / **`auto_reduce()` 方向**：按
  升序处理，只用更小首项约化
- **`reduce()` iteration limit** (`ocas-poly`): Increased from 200 to 10000
  for complex ideals / **`reduce()` 迭代上限**：从 200 提升到 10000
- **cyclic-4 test** (`ocas-tests`): Fixed incomplete generator set (was missing
  `abc+bcd+cda+dab` and `abcd-1`) / **cyclic-4 测试**：修复不完整的生成元集合

---

## [0.12.1] - 2026-07-06

### Added / 新增

- **Self-implemented NTT** (`ocas-poly`): Number Theoretic Transform for fast
  polynomial multiplication over ℤ_p in $O(n \log n)$. Radix-2 Cooley-Tukey
  with bit-reversal permutation. `DenseUnivariatePolynomial<FiniteField>::mul_ntt()`
  automatically selects NTT for degree ≥ 256 on NTT-friendly primes / **自研
  NTT**：ℤ_p 上快速多项式乘法，radix-2 Cooley-Tukey，degree ≥ 256 时自动启用
- **`pulp` SIMD dispatch** (`ocas-eval`): replaced `wide` with `pulp` for
  portable SIMD. Runtime CPU feature detection (SSE2/AVX2/AVX-512) with
  automatic lane width selection / **`pulp` SIMD 分派**：替换 `wide`，运行时
  CPU 特性检测，自动选择 SIMD 宽度
- **Estrin polynomial evaluation** (`ocas-eval`): `eval_estrin()` and
  `eval_estrin_batch()` via `fast_polynomial` crate for ILP-accelerated
  polynomial evaluation. Feature `fast-poly` / **Estrin 多项式求值**：通过
  `fast_polynomial` 利用指令级并行加速
- **Sparse matrix backend** (`ocas-poly`): `SprsMacaulayMatrix` adapter using
  `sprs` crate for F4 Macaulay matrix storage. Feature `sprs` / **稀疏矩阵后端**：
  使用 `sprs` 的 F4 Macaulay 矩阵适配器
- **Numerical verification** (`ocas-tests`): integration verification via
  `quadrature` crate (feature `verify-quadrature`); root-finding verification
  via bisection (feature `verify-roots`) / **数值验证**：`quadrature` 积分验证
  + 二分法求根验证
- **Feature matrix**: `ntt`, `sprs`, `pulp` (replaces `simd`'s `wide`),
  `fast-poly`, `verify-roots`, `verify-quadrature` / **Feature 矩阵**

### Changed / 变更

- **`simd` feature** now uses `pulp` instead of `wide`. The `wide` dependency
  has been removed from the workspace / **`simd` feature** 改用 `pulp` 替代 `wide`
- **`BuiltinOp` enum** (`ocas-eval`): `Instr::BuiltinFun { name: Symbol }` replaced
  by `Instr::BuiltinOp { op: BuiltinOp }`. Built-in functions are resolved at compile
  time, eliminating `to_lowercase()` + string matching on the SIMD hot path.
  SIMD trig throughput improved ~68% (1.9× → 3.2× on batch-4k) / **`BuiltinOp` 枚举**：
  内置函数编译时预分派，消除 SIMD 热路径上的字符串匹配，trig 吞吐提升 ~68%
- **SIMD stack buffer pre-allocation** (`ocas-eval`): `eval_simd_chunks` reuses a
  pre-allocated `Vec<[f64; 8]>` across chunks instead of allocating per chunk.
  SIMD poly throughput improved ~52% (6.6× → 10.0× on batch-4k) / **SIMD 栈缓冲区
  预分配**：chunk 间复用预分配缓冲区，poly 吞吐提升 ~52%
- **Montgomery modular multiplication** (`ocas-poly`): NTT hot path replaces
  `u128 % p` with Montgomery reduction (multiply + shift). NTT degree-1024
  throughput improved ~34% (999µs → 663µs, 90× vs Karatsuba) / **Montgomery 模乘**：
  NTT 热路径用 Montgomery 约减替代 128-bit 除法，degree-1024 提速 ~34%
- **NTT twiddle factor precomputation** (`ocas-poly`): `ntt_butterfly_mont`
  precomputes all stage roots once to avoid repeated `modpow` / **NTT 旋转因子
  预计算**：预计算所有层旋转因子，避免重复 modpow

---

## [0.12.0] - 2026-07-04

### Added / 新增

- **`RationalPolynomial<D, O>` type** (`ocas-poly`): rational function type
  with numerator/denominator as `SparseMultivariatePolynomial`, GCD-based
  canonicalization, and arithmetic (`add`, `sub`, `mul`, `div`, `neg`,
  `inv`, `pow`) / **有理多项式类型**：分子/分母表示 + GCD 规范化 + 四则运算
- **Brown PRS resultant** (`ocas-poly`): `DenseUnivariatePolynomial::resultant()`
  computes the resultant of two polynomials using Brown's Polynomial Remainder
  Sequence algorithm / **Brown PRS 结式**：多项式结式计算
- **Karatsuba multiplication** (`ocas-poly`): `DenseUnivariatePolynomial::mul_into`
  now uses Karatsuba fast multiplication (threshold=32) for large polynomials,
  replacing pure schoolbook O(n·m) / **Karatsuba 快乘法**：大次数多项式乘法加速
- **Polynomial extended GCD** (`ocas-poly`): `DenseUnivariatePolynomial::extended_gcd_poly()`
  returns `(g, s, t)` such that `s·self + t·other = g` / **多项式扩展 GCD**
- **Diophantine CRT** (`ocas-poly`): `DenseUnivariatePolynomial::diophantine()`
  solves the polynomial Chinese Remainder Theorem / **多项式 CRT 求解器**
- **p-adic expansion** (`ocas-poly`): `DenseUnivariatePolynomial::p_adic_expansion()`
  decomposes a polynomial with respect to another / **p-adic 展开**
- **Polynomial `pow()`** (`ocas-poly`): `DenseUnivariatePolynomial::pow(n)` by
  repeated squaring / **多项式幂运算**
- **Partial fraction decomposition** (`ocas-calc`): `apart()` decomposes a
  rational function into simpler fractions; `together()` combines them back /
  **部分分式分解**：`apart()` 分解 + `together()` 合并
- **Rational reconstruction** (`ocas-poly`): `rational_reconstruction(a, m)`
  recovers `(n, d)` from `a ≡ n/d (mod m)` using the extended Euclidean
  algorithm / **有理重构**：从模表示恢复有理数
- **Sparse polynomial helpers** (`ocas-poly`): `div_exact()`, `degree_in()`,
  `div_rem_sparse()` on `SparseMultivariatePolynomial` / **稀疏多项式辅助方法**
- **Dense polynomial helpers** (`ocas-poly`): `lcoeff()`, `constant()`,
  `mul_coeff()`, `div_coeff()`, public `content()` / **稠密多项式辅助方法**
- **Prelude expansion** (`ocas`): `RationalPolynomial` and `apart` now
  available via `use ocas::prelude::*` / **Prelude 扩展**：新增有理多项式和部分分式

### Changed / 变更

- **Dense polynomial multiplication** now routes through Karatsuba for
  polynomials with ≥32 coefficients / **稠密多项式乘法**改为 Karatsuba 路由

---

## [0.11.2] - 2026-07-04

### Added / 新增

- **Small-integer optimization (SOO) for `Integer`** (`ocas-domain`): values
  fitting in `i64` are stored on the stack in an `enum { Small(i64),
  Large(Box<rug::Integer>) }` with `UnsafeCell`-based lazy promotion,
  avoiding heap allocation for the common case. Arithmetic uses
  `i64::checked_add/sub/mul` with overflow fallback to GMP /
  **小整数 SOO 优化**：fit i64 的值栈分配，算术走 i64 快速路径，溢出回退 GMP
- **`to_i64()` method** on `Integer` (both GMP and non-GMP backends): extracts
  the value as `Option<i64>`, replacing the `.inner().try_into()` pattern /
  **`to_i64()` 方法**：替代 `.inner().try_into()` 模式
- **GMP backend binary serialization** (`ocas-domain`): `to_bigint()` and
  `From<BigInt>` use `write_digits`/`from_digits` instead of string conversion /
  **GMP 后端二进制序列化**：`to_bigint()` 和 `From<BigInt>` 使用二进制而非字符串
- **FiniteField optimization** (`ocas-domain`): cached `prime_minus_two` for
  fast `inv()` via Fermat's little theorem; `pow()` overrides default with
  `modpow`; GMP path caches `rug::Integer` versions for `pow_mod` /
  **有限域优化**：缓存 `prime-2` 加速求逆；`pow()` 使用 `modpow`；GMP 路径
  缓存 `rug::Integer` 版本
- **mimalloc global allocator** (`ocas` crate): optional `mimalloc` feature
  flag configures `mimalloc::MiMalloc` as `#[global_allocator]` /
  **mimalloc 全局分配器**：可选 `mimalloc` feature 配置全局分配器
- **Dense polynomial `mul_into()` buffer reuse** (`ocas-poly`):
  `DenseUnivariatePolynomial::mul_into(&self, other, buf)` writes into a
  caller-provided buffer, avoiding repeated allocation in hot loops /
  **稠密多项式 `mul_into()` 缓冲区复用**：热路径可复用缓冲区避免重复分配
- **Modular multivariate GCD** (`ocas-poly`): `gcd_modular(a, b)` reduces
  polynomials mod a prime, computes GCD in $\mathbb{F}_p$ via
  evaluation-interpolation, and lifts back to $\mathbb{Z}$ /
  **模方法多变量 GCD**：`gcd_modular` 在 $\mathbb{F}_p$ 中计算 GCD 后提升回 $\mathbb{Z}$
- **SOO benchmarks** (`ocas-tests/benches/gmp_integer.rs`): small/large
  add/mul/to_bigint/is_zero benchmarks for Integer /
  **SOO 基准**：小整数/大整数 add/mul/to_bigint/is_zero 基准
- **Modular GCD benchmark** (`ocas-tests/benches/poly_gcd.rs`): heuristic
  vs modular bivariate GCD comparison /
  **模方法 GCD 基准**：启发式 vs 模方法二元 GCD 对比

### Changed / 变更

- **Migrated all `.inner()` calls** across oCAS codebase to direct methods
  (`to_i64()`, `numer()`, `denom()`, `to_string()`, `Display`) /
  **迁移所有 `.inner()` 调用**为直接方法调用
- **rug 1.30 API compatibility**: `to_digits` → `write_digits`,
  `from_digits` via `impl Into<RugInteger>`, `ShrAssign` uses `unsafe`
  `UnsafeCell::get()` for sound mutation /
  **rug 1.30 API 适配**：`to_digits` → `write_digits` 等

---
## [0.11.1] - 2026-07-04

### Added / 新增

- **Bivariate integer polynomial factorization** (`ocas-poly`): monic-in-x
  bivariate factorization over $\mathbb{Z}[x,y]$ using Wang's Hensel lifting,
  with rational Bézout coefficients and integral correction reconstruction /
  **二元整数多项式因式分解**：基于 Wang Hensel 提升，使用有理 Bézout 系数
  与整系数修正重建
- **Bivariate finite-field polynomial factorization** (`ocas-poly`):
  monic-in-x bivariate factorization over $\mathbb{F}_p[x,y]$ using Hensel
  lifting / **二元有限域多项式因式分解**：基于有限域上的 Hensel 提升
- **Sparse multivariate `factor()` entry points** (`ocas-poly::sparse`):
  `SparseMultivariatePolynomial<IntegerDomain, Lex>::factor` and
  `SparseMultivariatePolynomial<FiniteField, Lex>::factor` /
  **稀疏多元 `factor()` 入口**：支持整数域和有限域
- **C/C++ polynomial bindings** (`ocas-c`): opaque `OcasPolyZ` and `OcasPolyFp`
  handles, string-based creation, factorization, degree, string output, and
  lifecycle functions / **C/C++ 多项式绑定**：不透明句柄、字符串创建、因式分解、
  次数、字符串输出与生命周期函数
- **mdBook factorization chapter** (`docs/book/{en,zh}/src/algorithms/factorization.md`):
  documentation covering univariate and bivariate factorization over $\mathbb{Z}$
  and $\mathbb{F}_p$, plus the C API / **mdBook 因式分解章节**：涵盖整数与有限域
  上的一元、二元因式分解及 C API 文档
- **Berlekamp algorithm enabled** for univariate factorization over prime finite
  fields / **Berlekamp 算法启用**：用于素有限域上一元因式分解

### Fixed / 修复

- **Unlucky evaluation points in Wang Hensel lifting**: when a chosen
  $y = \alpha$ produces a univariate factorization that is inconsistent with the
  bivariate factorization, the implementation now tries additional candidates
  and falls back to irreducible instead of panicking / **Wang Hensel 提升中的
  不幸赋值点**：现在会尝试额外候选点，并回退为不可约而非 panic

---
## [0.11.0] - 2026-07-03

### Added / 新增

- **Correctness comparison framework** (`ocas-tests/tests/correctness/`):
  16 modules and 82 tests covering parse, normalize, rewrite, calculus,
  evaluation, polynomial arithmetic, GCD, factorization, Gröbner bases,
  finite fields, matrices, root isolation, and linear solvers / **正确性对比框架**：
  16 个模块 82 个测试，覆盖解析、化简、重写、微积分、求值、多项式运算等
- **SymPy reference harness** (`scripts/compare_sympy.py`): supports
  `check`, `verify` (JSON via stdin), and `time` modes to compare oCAS
  output against SymPy 1.14 / **SymPy 参考工具**：支持 check/verify/time 三种模式
- **Symbolica comparison harness** (`scripts/compare_symbolica.py` and
  `symbolica_runner/`): isolated Rust crate running Symbolica via subprocess,
  keeping the AGPL dependency out of the main build / **Symbolica 对比工具**：
  独立子进程调用，AGPL 不链入主构建
- **Audit report generator** (`scripts/generate_audit_report.py`): runs
  simple/medium and complex/very_complex tests separately and writes a
  Markdown summary / **审计报告生成器**：分别运行不同难度测试并生成 Markdown
- **Difficulty tier annotations**: `#[ignore]` marks complex and very complex
  tests for manual/audit runs while simple+medium tests run in CI /
  **难度分级**：complex/very_complex 测试标 `#[ignore]` 供手动审计运行
- **Complete polynomial factorization over Z**: `factor()` on
  `DenseUnivariatePolynomial<IntegerDomain>` (Yun → CZ mod p → Hensel → Zassenhaus).
  Handles `x^100 - 1` into 9 cyclotomic factors. / **Z 上完整因式分解**
- **Complete factorization over F_p**: `factor()` on
  `DenseUnivariatePolynomial<FiniteField>` (Cantor–Zassenhaus DDF+EDF, char-p
  Bernardin square-free). / **F_p 上完整因式分解**
- **Number-theory primitives**: `ocas_domain::number_theory` — Miller–Rabin,
  Chinese remainder, Legendre/Jacobi, Tonelli–Shanks, modular inverse. /
  **数论原语**
- **Multivariate GCD**: `heuristic_gcd` on `SparseMultivariatePolynomial`,
  plus `content`/`primitive_part`/`eval`. / **多元 GCD**
- **Python `Polynomial.factor()`** / **Proptest 500-case roundtrip** /
  **Criterion benchmarks** (`poly_factor_z`, `poly_factor_fp`)

### Changed / 变更

- `bench_sympy.py` now delegates to `compare_sympy.py time` mode so timing
  and correctness share the same SymPy task definitions /
  `bench_sympy.py` 改为委托 `compare_sympy.py` 的 `time` 模式

### Notes / 说明

- All 82 correctness tests pass (`--include-ignored`) / 82 个正确性测试全部通过。
- Two known gaps are documented in tests: `sin(x)^2 + cos(x)^2 -> 1` is not
  handled by the default simplifier, and the real-root isolator finds 8 of
  10 roots for the expanded Wilkinson n=10 polynomial /
  已记录两个已知差距：默认化简器不处理 `sin²+cos²->1`；Wilkinson n=10
  展开多项式实根隔离只找到 8 个根。

---
## [0.10.0] - 2026-07-02

### Added / 新增

- **Python `Polynomial` class**: dense univariate polynomials over ℤ, ℚ, or
  GF(p) via an enum-erasure strategy. Exposes construction, `coeffs`,
  `degree`, `eval`, arithmetic operators, `derivative`, `integral`, `gcd`,
  `div_rem`, `square_free_factorization`, `primitive_part`, and finite-field
  evaluation / **Python `Polynomial` 类**：基于枚举擦除的稠密一元多项式，
  支持整数、有理数、有限域 GF(p) 三种系数域
- **Python `Matrix` class**: dense matrices over ℤ, ℚ, GF(p) with
  `transpose`, `trace`, `determinant` (Bareiss), `rank`, `inverse`, `matmul`
  (`@`), `solve`, and arithmetic operators / **Python `Matrix` 类**：支持
  整数、有理数、有限域的稠密矩阵
- **Python `IntegerDomain` / `RationalDomain` / `FiniteField` classes**:
  coefficient-domain selectors consumed by `Polynomial` and `Matrix` /
  **Python 系数域类**：供 `Polynomial` 与 `Matrix` 使用的域选择器
- **`Matrix` Rust API**: new methods `transpose`, `trace`, `matmul`, `rank`,
  `determinant` (Bareiss fraction-free with partial pivoting), `inverse`,
  `row`, `column` / **`Matrix` Rust API**：新增矩阵方法
- **`FiniteField` now implements `EuclideanDomain`**: enables
  `Matrix<FiniteField>` and polynomial operations over finite fields /
  **`FiniteField` 实现 `EuclideanDomain`**：支持有限域上的矩阵与多项式运算
- **`Display` for `Rational` and `FiniteFieldElement`**: consistent string
  rendering across gmp/non-gmp builds / **`Rational` 与
  `FiniteFieldElement` 的 `Display` 实现**
- **New benchmarks**: `poly_gcd` (univariate GCD), `poly_factor`
  (square-free factorization), `groebner` (cyclic-n Gröbner bases) /
  **新基准**：多项式 GCD、因式分解、Gröbner 基
- **Extended SymPy comparison**: `factor`, `gcd`, `series`, and large
  expansion tasks / **扩展 SymPy 对比**：新增 factor/gcd/series 任务
- **SageMath comparison harness** (`scripts/bench_sage.py`): local manual
  comparison mirroring `bench_sympy.py` / **SageMath 对比工具**
- **mdBook documentation site**: bilingual introduction, getting-started,
  architecture, Python/C bindings, performance, backends, and contributing
  chapters, deployed via `docs.yml` GitHub Actions / **mdBook 文档站点**
- **`docs.rs` metadata**: `all-features = true` so the online API docs show
  every backend / **`docs.rs` 元数据**：全 feature 在线文档
- **README badges**: CI and docs.rs status badges / **README 徽章**

### Changed / 变更

- Workspace version bumped from `0.9.0` to `0.10.0` / 工作区版本提升
- README status updated from "Alpha (0.4.0)" to "Beta (0.10.0)" / README 状态
  更新
- 0.10.0 ROADMAP deliverables marked complete / ROADMAP 0.10.0 交付物标记完成

### Notes / 说明

- Symbolica comparison is documented as a local manual workflow (running
  Symbolica's example binaries) rather than linked into the build, owing to
  Symbolica's AGPL license and separate workspace / 因 Symbolica 的 AGPL
  许可证与独立 workspace，Symbolica 对比以本地手动工作流文档化，不链接进构建
- `Polynomial` and `Matrix` Python objects define `__eq__` but are not
  hashable (a pyo3 0.29 limitation); they behave like Python `list`/`dict`
  in this respect / `Polynomial` 与 `Matrix` Python 对象定义了 `__eq__` 但
  不可哈希（pyo3 0.29 限制），行为与 Python `list`/`dict` 一致

---
## [0.9.0] - 2026-07-02

### Added / 新增

- **Python bindings** (`ocas-py`): PyO3-based Python module exposing
  `Expression`, `ExpressionEvaluator`, and `DiophantineSolution` classes,
  plus `solve_linear_rational`, `solve_linear_integer`, and
  `solve_diophantine` functions. Each `Expression` owns a private arena
  for self-contained lifetime management. `Polynomial`, `Matrix`, and
  `Domain` classes are deferred to 0.10.0 / **Python 绑定**：基于 PyO3
  的 Python 模块，暴露 `Expression`、`ExpressionEvaluator`、
  `DiophantineSolution` 类与 `solve_*` 函数
- **`Expression` Python class**: Parse, `diff`, `integrate`, `taylor`,
  `simplify`, `substitute`, `normalize`, `clone`, plus operator overloads
  (`__add__`, `__sub__`, `__mul__`, `__pow__`, `__neg__`, `__eq__`) /
  **`Expression` Python 类**
- **`ExpressionEvaluator` Python class**: Compile-once / evaluate-many
  numeric evaluator over `f64` with parameter binding /
  **`ExpressionEvaluator` Python 类**
- **`pyproject.toml`**: Maturin build backend; module name `ocas`;
  `pip install ocas` ready on Linux/macOS/Windows (pure Rust wheel) /
  **`pyproject.toml`**：Maturin 构建后端
- **C expression lifecycle API** (`ocas-c`): `ocas_expr_parse`,
  `ocas_expr_free`, `ocas_expr_clone`, `ocas_expr_to_string`,
  `ocas_expr_normalize`, `ocas_string_free` / **C 表达式生命周期 API**
- **C calculus API**: `ocas_expr_diff`, `ocas_expr_integrate`,
  `ocas_expr_taylor`, `ocas_expr_simplify`, `ocas_expr_substitute` /
  **C 微积分 API**
- **Extended error model** (`ocas-c`): New error codes `OCAS_ERROR_PARSE`,
  `OCAS_ERROR_INVALID_ARGUMENT`, `OCAS_ERROR_DIVISION_BY_ZERO`,
  `OCAS_ERROR_OUT_OF_MEMORY` / **扩展错误模型**
- **C++ RAII wrapper** (`ocas-c/include/ocas.hpp`): `ocas::Expression`
  with automatic resource management, move/copy semantics, and exception
  translation / **C++ RAII 包装**
- **New C examples**: `examples/expression.c` and
  `examples/cpp_example.cpp` demonstrating the full expression lifecycle /
  **新 C/C++ 示例**
- **C API integration tests** (`ocas-c/tests/c_api.rs`): End-to-end tests
  exercising the `#[no_mangle] extern "C"` functions through Rust FFI /
  **C API 集成测试**
- **C example compilation test** (`ocas-c/tests/examples_compile.rs`):
  Compiles and runs `examples/expression.c` against the built static
  library, verifying the C example compiles and runs (ROADMAP success
  criterion) / **C 示例编译验证测试**
- **Python pytest suite** (`ocas-tests/tests/python/`): 33 tests covering
  parsing, calculus, simplification, substitution, operators, solvers,
  numeric evaluation, hash/eq semantics, and memory pressure cycles /
  **Python pytest 套件**：33 个测试
- **Wheel build CI** (`.github/workflows/wheels.yml`): Maturin-based
  matrix build for Linux/macOS/Windows with PyPI publishing on tag /
  **Wheel 构建 CI**

### Changed / 变更

- Workspace version bumped to `0.9.0` / 工作区版本提升至 `0.9.0`
- `ocas-c` refactored into `error.rs` + `expression.rs` modules;
  `crate-type` now includes `rlib` for integration testing / `ocas-c`
  重构为模块化结构
- `ocas_version()` now uses `env!("CARGO_PKG_VERSION")` instead of a
  hardcoded string / `ocas_version()` 改用编译期版本
- `cbindgen.toml`: Added `include_guard = "OCAS_H"`, `usize_is_size_t`,
  `style = "tag"` / `cbindgen.toml` 补全
- `ocas-eval`: `compile` module made public; `compile_atom`,
  `compile_atom_with`, `compile_tree`, `compile_tree_with` re-exported /
  `ocas-eval` 编译模块公开
- Top-level `ocas` prelude: Added `DiophantineSolution`, `PowfExtension`;
  crate-root flat exports now include `solve_*`, `Assumption*`, `Matrix`,
  `GroebnerBasis`, `buchberger`, `monomial_*`, `RootInterval` / 顶层
  prelude 一致性修复
- `ocas-py` and `ocas-c` now depend on `ocas-rewrite` (previously
  missing, blocking `simplify`/`transform`) / `ocas-py` 与 `ocas-c`
  补齐 `ocas-rewrite` 依赖

### Fixed / 修复

- **Panic-safe arena recovery**: `build()` in both `ocas-c` and `ocas-py`
  now uses an `ArenaGuard` RAII wrapper, ensuring leaked arenas are
  recovered even if `normalize` or the builder closure panics /
  **Panic 安全的 arena 回收**
- **C++ namespace mismatch**: Removed `namespace = "ocas"` from
  `cbindgen.toml` so generated types match the global-scope references in
  `ocas.hpp` (the C++ RAII wrapper now compiles) / **C++ 命名空间不匹配**
- **C error code consistency**: `cstr_to_str` failures now propagate the
  correct error code via `error::write_last_code()` instead of
  hardcoding `OCAS_ERROR_NULL_POINTER` for UTF-8 errors /
  **C 错误码一致性**
- Removed unused dependencies from `ocas-py` (`ocas-domain`,
  `ocas-poly`) and `ocas-c` (`ocas-domain`, `ocas-poly`, `ocas-eval`) /
  删除未使用的依赖

---
## [0.8.0] - 2026-07-02

### Added / 新增

- **`ocas-eval` crate**: Expression evaluation engine with stack-based VM
  interpreter, AST compiler, function registry, optimizer, SIMD vectorization,
  and Cranelift JIT backend / **表达式求值引擎**：栈式 VM 解释器、AST 编译器、
  函数注册表、优化器、SIMD 向量化、Cranelift JIT 后端
- **`ExpressionEvaluator<T>`**: Generic stack VM that compiles `Atom` trees
  into linear instruction sequences and evaluates them with user-provided
  parameter values / **泛型栈式 VM**：将 `Atom` 树编译为线性指令序列并按
  用户提供的参数值求值
- **`VectorEvaluator`** (`simd` feature): Batch evaluation of expressions
  using `wide::f64x4` SIMD primitives, processing 4 lanes in parallel with
  scalar fallback for remainders / **SIMD 向量化求值器**：使用 `wide::f64x4`
  进行 4 路并行批量求值，余数回退标量计算
- **`EvalTree`**: Owned intermediate representation that decouples compilation
  from the arena lifetime / **自有中间表示**：将编译与 arena 生命周期解耦
- **`FunctionMap<T>`**: User-defined function registry with name resolution,
  case-insensitive lookup, aliases, and index-based calling /
  **用户自定义函数注册表**：支持名称解析、大小写不敏感查找、别名、索引调用
- **Instruction optimizer**: Common subexpression elimination (CSE), dead code
  elimination, and algebraic simplification (single-element Add/Mul → Copy) /
  **指令优化器**：公共子表达式消除、死代码消除、代数简化
- **Cranelift JIT backend** (`jit` feature): Compiles instruction sequences
  to native machine code via Cranelift 0.117 (experimental; runtime tuning
  in progress) / **Cranelift JIT 后端**：通过 Cranelift 将指令序列编译为
  原生机器码（实验性，运行时调优进行中）
- **`EvaluationDomain` trait**: Numeric evaluation trait with built-in
  function table (sin/cos/exp/log/sqrt/abs/tan/sec/csc/cot), case-insensitive
  function names, and `f64` implementation / **数值求值 trait**：含内置函数
  表、大小写兼容函数名、`f64` 实现
- Top-level prelude exports: `ExpressionEvaluator`, `VectorEvaluator`,
  `FunctionMap`, `EvalTree`, `EvaluationDomain`, `EvaluationError`,
  `Instr`, `Instruction`, `Slot` / **顶层 prelude 导出**

### Changed / 变更

- Workspace version bumped to `0.8.0` / 工作区版本提升至 `0.8.0`
- Removed `llvm` feature from `ocas` and `ocas-eval`; LLVM backend deferred
  to Post-1.0 / 从 `ocas` 和 `ocas-eval` 中移除 `llvm` feature；LLVM 后端
  推迟到 1.0 之后
- `ROADMAP.md`: 0.8.0 deliverables marked complete; LLVM moved to Post-1.0

---
## [0.7.0] - 2026-07-01

### Added / 新增

- **Assumptions system** (`ocas-domain`): `Assumption` enum with 14 predicates
  (`Real`, `Positive`, `Integer`, `Even`, `Prime`, …), `Assumptions` set with
  logical implication and conflict detection, `SymbolAssumptions` map /
  **假设系统**：包含 14 种谓词的 `Assumption` 枚举、带逻辑蕴含与冲突检测的
  `Assumptions` 集合、`SymbolAssumptions` 映射
- **Matrix types** (`ocas-poly`): `Matrix<D>` with fraction-free Gaussian
  elimination, back-substitution, and `solve()` for exact linear systems /
  **矩阵类型**：带分数无关高斯消元、回代与 `solve()` 的 `Matrix<D>`
- **Polynomial GCD** (`ocas-poly`): `gcd()` with pseudo-remainder algorithm,
  `primitive_part`, `content` for `DenseUnivariatePolynomial` /
  **多项式 GCD**：基于伪余式算法的 `gcd()`、`primitive_part`、`content`
- **Square-free factorization** (`ocas-poly`): Yun's algorithm for
  `DenseUnivariatePolynomial` / **无平方因子分解**：Yun 算法
- **Real root isolation** (`ocas-poly`): Sturm sequences, `count_real_roots`,
  `isolate_real_roots`, bisection refinement / **实根隔离**：Sturm 序列、
  `count_real_roots`、`isolate_real_roots`、二分精化
- **Gröbner bases** (`ocas-poly`): `GroebnerBasis` with Buchberger's algorithm,
  coprime criterion, minimization, auto-reduction, `is_groebner_basis` /
  **Gröbner 基**：`GroebnerBasis` 含 Buchberger 算法、互质准则、最小化、自归约
- **Linear system solver** (`ocas-calc`): `solve_linear_rational` and
  `solve_linear_integer` / **线性方程组求解器**
- **Diophantine equation solver** (`ocas-calc`): `solve_diophantine` for
  `ax + by = c` / **丢番图方程求解器**
- **Polynomial system solver** (`ocas-calc`): `solve_polynomial_system` using
  lexicographic Gröbner bases / **多项式系统求解器**：基于 Lex Gröbner 基
- Sparse polynomial Gröbner helpers: `leading_term`, `reduce`, `spoly`,
  `mul_monomial`, `monomial_divides`, `monomial_lcm`, `monomial_are_coprime` /
  稀疏多项式 Gröbner 辅助方法

### Changed / 变更

- `EuclideanDomain` trait: added `gcd()` and `extended_gcd()` default methods /
  `EuclideanDomain` trait 新增 `gcd()` 和 `extended_gcd()` 默认方法
- `Integer` and `Rational` types: added `Display` implementation
- `DenseUnivariatePolynomial`: added `is_one()`
- `SparseMultivariatePolynomial`: added Gröbner-basis support methods
- `ocas` prelude: exposed all new types and functions / 暴露所有新类型与函数
- Workspace version bumped to `0.7.0` / 工作区版本提升至 `0.7.0`

---
## [0.6.0] - 2026-07-08

### Added / 新增

- Stable top-level `ocas` prelude API / 稳定的顶层 `ocas` prelude API
- Rustdoc examples for all public items in `ocas`, `ocas-atom`, `ocas-calc`, `ocas-parse`, `ocas-rewrite`, and `ocas-core` / 为 `ocas`、`ocas-atom`、`ocas-calc`、`ocas-parse`、`ocas-rewrite` 和 `ocas-core` 的所有公共项添加 rustdoc 示例
- Property-based tests via `proptest` in `ocas-atom`, `ocas-rewrite`, `ocas-calc`, `ocas-parse`, and `ocas-domain` / 在 `ocas-atom`、`ocas-rewrite`、`ocas-calc`、`ocas-parse` 和 `ocas-domain` 中通过 `proptest` 添加基于属性的测试
- Criterion benchmarks for parsing, normalization, dense/sparse polynomials, calculus, rewriting, and a SymPy comparison harness via `uv` / 用于解析、规范化、稠密/稀疏多项式、微积分、重写的 Criterion 基准，以及通过 `uv` 运行的 SymPy 对比基准
- `substitute` exported from `ocas-calc` and the top-level prelude / 从 `ocas-calc` 和顶层 prelude 导出 `substitute`
- `#[doc(hidden)]` subcrate re-exports so the top-level `ocas` API is the documented public surface / 子 crate 重新导出标记为 `#[doc(hidden)]`，使顶层 `ocas` API 成为文档中的公共接口

### Changed / 变更

- Workspace version bumped to `0.6.0` / 工作区版本提升至 `0.6.0`
- Internal workspace crates now reference each other through `[workspace.dependencies]` so they are ready for publication / 内部工作区 crate 现在通过 `[workspace.dependencies]` 互相引用，以便发布

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

## [0.1.0] - 2026-06-29

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