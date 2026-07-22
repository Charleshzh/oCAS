# Changelog / 变更日志

All notable changes to the oCAS project will be documented in this file.

oCAS 项目的所有重大变更都将记录在此文件中。

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.18.0] - 2026-07-23

### Added / 新增

- **fuel 资源控制**：新增 `ocas-core::fuel` 模块，`Fuel` 为 `Arc<AtomicUsize>`
  共享递减预算（clone 共享计数器），`consume`/`check`/`remaining` + `OcasError::OutOfFuel`；
  `simplify_with_fuel` 与 `integrate_with_fuel` 接入（保留旧 API 冻结）/ **fuel
  resource control**: new `ocas-core::fuel` module with `Fuel` as a shared
  decrementing budget, `OutOfFuel` error, and `simplify_with_fuel` /
  `integrate_with_fuel` entry points (legacy APIs preserved).
- **超对偶数自动微分**：新增 `ocas-domain::dual` 模块，`HyperDual<T>` 运行时形状
  + 截断乘法表 + 几何级数求逆；`DualCoeff` trait + `new_first_order` 构造器；
  为 `Rational` 补 `std::ops` impl（gmp 与非 gmp 双路径）/ **Hyper-dual forward
  AD**: `ocas-domain::dual` with runtime-shaped `HyperDual<T>`, truncated
  multiplication table, geometric-series inverse, and `std::ops` impls for
  `Rational` on both backends.
- **Vegas 自适应蒙特卡洛积分**：新增 `ocas-eval::numeric` 模块，`Vegas` 积分器
  （累积弧长网格更新，数值稳定版）+ `StatisticsAccumulator`（inverse-variance
  跨迭代合并）+ `Integrator` trait + `integrate_1d` 顶层入口 / **Vegas adaptive
  Monte Carlo integration**: `ocas-eval::numeric` with the `Vegas` integrator
  (numerically stable cumulative-arc-length grid update),
  `StatisticsAccumulator`, `Integrator` trait, and `integrate_1d`.
- **张量基础**：新增 `ocas-atom::tensor` 模块，独立 `Tensor` 类型（指标槽 +
  `IndexPosition` + `Symmetry`）、`contract` 显式收缩、`symmetrise_sign` 对称化
  符号（完整图规范化推迟 Post-1.0）/ **Tensor basics**: `ocas-atom::tensor`
  with an independent `Tensor` type (index slots + variance + symmetry),
  explicit `contract`, and `symmetrise_sign` (full graph canonicalisation
  deferred to Post-1.0).
- 测试：fuel 8 项、dual 11 项单元 + 3 项跨域 proptest（`dual_vs_diff.rs`，
  与 `ocas_calc::diff` 对照）、numeric 9 项（含高斯峰自适应）、tensor 8 项 /
  Tests: 8 fuel, 11 dual unit + 3 cross-domain proptests vs `diff`, 9 numeric
  (incl. Gaussian-peak adaptation), 8 tensor.

### Changed / 变更

- workspace 版本 0.17.1 → 0.18.0（13 crate）；新增 `rand`/`rand_xoshiro` 依赖 /
  workspace 0.17.1 → 0.18.0; added `rand`/`rand_xoshiro` dependencies.

## [0.17.1] - 2026-07-22

### Added / 新增

- **Python 绑定 — 代数数域**：新增 `ocas-py::algebraic` 模块，暴露
  `AlgebraicExtension`（由首一极小多项式升序系数构造，如
  `AlgebraicExtension([-2, 0, 1])` 即 $\mathbb{Q}(\sqrt{2})$）、
  `AlgebraicElement`、`AlgebraicPolynomial`（系数项可为整数、
  `(num, denom)` 元组、$\alpha$-多项式系数列表或 `AlgebraicElement`）、
  `AlgebraicFactor`，并支持 Trager 因式分解（`factor()`）/ **Python
  bindings for algebraic number fields**: new `ocas-py::algebraic` module
  exposing `AlgebraicExtension`, `AlgebraicElement`, `AlgebraicPolynomial`,
  and `AlgebraicFactor` with Trager factorization.
- **C/C++ 绑定 — 代数数域**：新增 `ocas-c::algebraic` 模块，不透明句柄
  `OcasAlgebraicField` / `OcasAlgebraicPoly` / `OcasAlgebraicFactorArray`，
  以及 `ocas_algebraic_field_create`（极小多项式字符串如 `"x^2 - 2"`）、
  `ocas_algebraic_poly_create`（系数列表字符串，多项式系数以 `;` 分隔、
  每项 $\alpha$-多项式以 `,` 分隔）、`ocas_algebraic_poly_factor` 等
  C ABI；`include/ocas.h` 已重新生成 / **C/C++ bindings for algebraic
  number fields**: opaque handles and `ocas_algebraic_*` C ABI for
  creating fields and polynomials and factoring via Trager; `include/ocas.h`
  regenerated.
- `RootOf(poly, root_index)` 经现有递归下降解析器自动产生
  `Fun("RootOf", ...)` 节点（无需词法/语法改动）/ `RootOf(poly, idx)`
  parses directly to an `Atom` function node with no parser changes.
- 测试：13 项 Python 测试（`test_algebraic.py`）与 7 项 C API 测试
  （含 $\mathbb{Q}(\sqrt{2})$、$\mathbb{Q}(\sqrt[3]{2})$ 因式分解、
  $x^2-\alpha$ 不可约性验证）/ Tests: 13 Python and 7 C API cases
  covering ANF factorization and irreducibility over $\mathbb{Q}(\sqrt{2})$.

## [0.17.0] - 2026-07-22

### Added / 新增

- **代数数域因式分解（Trager 算法）**：新增 `ocas-domain::algebraic` 模块，
  `AlgebraicExtension<D>` 表示 $D[\alpha]/(m)$（$m$ 首一；$D=\mathbb{Q}$
  即代数数域 `AlgebraicNumberField`，$D=\mathbb{F}_p$ 即
  $\mathrm{GF}(p^d)$），求逆经自包含的扩展 Euclid 实现 / **Algebraic
  number field factorization (Trager's algorithm)**: new
  `ocas-domain::algebraic` module with `AlgebraicExtension<D>` for
  $D[\alpha]/(m)$ (`AlgebraicNumberField` over $\mathbb{Q}$,
  $\mathrm{GF}(p^d)$ over $\mathbb{F}_p$), inversion via a self-contained
  extended Euclidean algorithm.
- `DenseUnivariatePolynomial<AlgebraicNumberField>::factor()`：Yun 无平方
  分解（模 GCD）→ 平移范数（求值–插值结式，mod-$p$ 无平方检验）→
  $\mathbb{Q}$ 上分解 → 数域模 GCD（$\mathrm{GF}(p^d)$ + CRT + 有理重构
  + 试除）→ 回代首一化；有理系数输入走快速通道（先在 $\mathbb{Q}$ 上
  分解）/ `factor()` over $\mathbb{Q}(\alpha)$: Yun square-free stage
  with modular GCD, shifted norm via evaluation–interpolation resultants,
  rational factorization, modular GCD over $\mathrm{GF}(p^d)$ with CRT and
  rational reconstruction, and a rational fast path.
- 稀疏 Diophantine 小素数升级启发式：骨架组大小超过 $p-1$ 时
  `padic_lift_factors` 第一轮素数扫描自动升级到更大素数，第二轮才允许
  回退稠密 Diophantine / Sparse Diophantine small-prime escalation: the
  first prime pass escalates past primes too small for skeleton
  interpolation instead of degrading to the dense solver.
- 顶层 `ocas` crate prelude 导出 `AlgebraicElement`、`AlgebraicExtension`、
  `AlgebraicNumberField` / Prelude re-exports for the new ANF types.
- 测试与基准：13 项单元测试（含 Symbolica `algebraic_extension`、
  `gcd_number_field` 镜像）、21 项 correctness 用例（SymPy
  `factor(extension=...)` 交叉验证）、proptest 往返（ignore，手动运行）、
  criterion 组 `poly_factor_anf`（$\mathbb{Q}(\sqrt2)$ deg 12 ≈ 8 ms、
  $\mathbb{Q}(\sqrt[3]{2})$ deg 9 ≈ 24 ms、$\mathbb{Q}(\zeta_5)$ deg 9
  ≈ 32 ms，均远低于 100 ms 指标）/ Tests, correctness cases, proptest
  roundtrip, and `poly_factor_anf` benchmarks (all well under the 100 ms
  target).

### Fixed / 修复

- **结式（Brown PRS）在一般次数下结果错误**：原实现仅在 $\beta$ 为单位
  时才执行 $\beta$ 除法，不是合法的结式算法（小次数凑巧正确）。已按
  Symbolica `resultant_prs` 移植为逐步精确除以 $\beta$ 的
  subresultant PRS；新增回归测试
  $\operatorname{Res}(x^4-3,\,3x^3-x^2+2x+1)=-2243$（SymPy 验证） /
  **Resultant (Brown PRS) wrong beyond trivial degrees**: the beta division
  was applied only when beta was a unit; ported Symbolica's `resultant_prs`
  with exact per-step division. Regression test against SymPy added.

---

## [0.16.2] - 2026-07-22

### Added / 新增

- **$\mathbb{F}_p$ 非常数首项系数因式分解**：$\mathbb{F}_p$ 多元因式分解现支持
  非常数首项系数，通过 Wang 首项系数重建（`wang_reconstruct_lcoeffs_fp`）
  分配 LC 因子后以 `eez_lift_imposed` 提升 / **$\mathbb{F}_p$ non-constant
  leading-coefficient factorization** via Wang LC reconstruction (`wang_reconstruct_lcoeffs_fp`)
  distributing LC factors among polynomial factors, then `eez_lift_imposed`.
- `find_sample_fp` 新增 `lc_filter` 参数：过滤 LC 因子在采样点为零的样本 /
  `find_sample_fp` gains `lc_filter` parameter to skip samples where any
  LC factor evaluates to zero.
- `find_sample_z` 新增分解上限（200 次有效单变量分解）+ 自适应值域中间档
  （7→15→25）/ `find_sample_z` decomposition cap (200) and intermediate
  value-bound level (7→15→25).
- 测试：3 项 $\mathbb{F}_p$ 非常数 LC 用例（二元、三元、可约 LC）+ 1 项 4
  变量（ignored）/ 3 Fp non-constant LC tests + 1 ignored 4-variable test.

### Changed / 变更

- `multivariate_factor_fp` 文档更新：移除 "LC must be constant" 限制声明 /
  `multivariate_factor_fp` doc updated: removed LC-must-be-constant restriction.
- `factor_square_free_fp` 移除非常数 LC 放弃分支，改走 Wang LC + `eez_lift_imposed`
  路径 / `factor_square_free_fp` non-constant LC now routes through Wang LC +
  `eez_lift_imposed` instead of returning unfactored.

---

## [0.16.1] - 2026-07-22

### Added / 新增

- **非常数首项系数强加**：$\mathbb{Z}$ 多元因式分解现可通过 p-adic
  系数 Hensel 提升实现非平凡（非常数）首项系数的完整强加，移植自
  Symbolica `sparse_coefficient_hensel_lift_mod_prime`（密 Diophantine
  变体） / **Non-constant leading-coefficient imposition** for $\mathbb{Z}$
  multivariate factorization via a p-adic coefficient Hensel lift porting
  Symbolica's `sparse_coefficient_hensel_lift_mod_prime`.
- **稀疏多元 Diophantine 求解器**：骨架插值（Vandermonde + 移位
  解）+ 两/多因式 EEA 序列求解，移植自 Symbolica
  `sparse_multivariate_diophantine_*`（含 ≥512 项上限保护） / **Sparse
  multivariate Diophantine solver** via skeleton interpolation + Vandermonde,
  porting Symbolica's `sparse_multivariate_diophantine_*`.
- 自适应采样：ℤ 路径候选排序引入 `(因式数, content)` 键（content=1
  优先）；候选去重（`HashSet`）；非首一 LC 因子镜像过滤；值域自适应
  递增（7 → 25）；Fp 路径两轮递增范围（8 → 32） / Adaptive sampling:
  dedup, LC-factor filter, content-aware candidate ranking, value-bound
  escalation (ℤ 7→25, $\mathbb{F}_p$ 8→32).
- 二元非常数 LC 改走 EEZ 路径（`sparse.rs::factor()` 自动分派） / Bivariate
  non-constant-LC dispatches to the EEZ path via `factor()`.
- 审计报告扩展：`generate_audit_report.py` 新增 Symbolica 因式分解计时
  对比表；`symbolica_runner` 新增 `factor_time` 任务 / Audit report
  extended with Symbolica factorization timing comparison.
- 测试/基准：新增 `z_nonconstant_lcoeff_shared_monomial`、
  `z_nonconstant_lcoeff_reducible_lc`、`z_coefficient_lift_uses_sparse_diophantine`
  回归测试；correctness 框架新增 4 个非常数 LC 用例（二元、三元、可约
  LC、4 变量 ≥50 项）；criterion 组新增 `trivariate_nonconstant_lcoeff`
  与 `sparse_4var_nonconstant_lcoeff` 基准；proptest `factor_product_of_two_nonconstant_lc` / Tests, correctness cases, and benchmarks for non-constant LC paths.

### Fixed / 修复

- **EEZ 多变量 Diophantine 契约违反**：第 $k$ 变量提升中修正项引入
  $x_k$ 后 `diophantine` 基例 `mpoly_to_dense` 在 $x_k=1$ 处折叠求和。
  修复：求解前对因子做 `eval_keep(k, a_k)`（对齐 Symbolica
  `factors_mod = replace(last_var, 0)`）。0.16.0 的全部线性因子测试
  未触发此路径 / **EEZ multivariate Diophantine contract violation**:
  correction terms introducing $x_k$ were incorrectly summed by
  `mpoly_to_dense` at $x_k = 1$; fixed by evaluating the factors at the
  sample before solving.
- **稀疏采样样本项系数在幂次提升时被错误平方**：per-term 求值含系数
  时，采样 $s≥1$ 会将系数一并提升。修复：求值仅存单项式值，系数在
  构建像时乘入 / **Term coefficient squaring in sparse sampling**: fixed
  by separating monomial evaluation from coefficient multiplication.
- `find_sample_z` 早退（2 候选即停）导致非常数 LC 的候选全部无效；
  修复：候选去重 + 移除早退 + LC 因子镜像过滤（$|g_j(s)|>1$）+
  content 排序（content=1 优先） / `find_sample_z` early exit with
  only 2 candidates could yield all invalid candidates for non-constant LC;
  fixed with dedup, content-aware ranking, and LC-factor filtering.
- `div_rem_sparse` 整除判断 `monomial_divides` 参数方向错误（无调用方
  未暴露） / `monomial_divides` argument order bug (latent).

---

## [0.16.0] - 2026-07-21

### Added / 新增

- **任意多元因式分解**（`ocas-poly::factor::eez`）：$\mathbb{Z}$ 与
  $\mathbb{F}_p$ 上任意变量数的多项式因式分解，核心是 Wang EEZ
  （求值–提升）算法——无平方分解（Yun + 特征 $p$ 的 $p$ 次幂处理）、
  采样点搜索、Wang 首项系数预处理、逐变量 Hensel 提升（多元
  Diophantine 方程）、Zassenhaus 重组 / **Arbitrary multivariate
  factorization** over $\mathbb{Z}$ and $\mathbb{F}_p$ via Wang's EEZ
  (evaluation–lifting) algorithm: square-free factorization (Yun with
  characteristic-$p$ $p$-th power handling), sample-point search, Wang
  leading-coefficient preprocessing, variable-by-variable Hensel lifting via
  multivariate Diophantine equations, and Zassenhaus recombination.
- **$n$ 元 GCD**（`ocas-poly::multivariate_gcd`）：稠密递归求值–插值的
  `multivariate_gcd_field`/`_z`/`_fp`，作为无平方分解与多元分解的前置 /
  **$n$-variate GCD** via dense recursive evaluation–interpolation.
- `SparseMultivariatePolynomial` 多元分解辅助 API：`coeff_of_var_pow`、
  `leading_coeff_in`、`derivative(var)`、`taylor_coefficients`、
  `drop_main_var`、`embed_new_main`、`permute_variables`、
  `checked_div_exact`、`eval_keep` / multivariate-factorization helper
  APIs on `SparseMultivariatePolynomial`.
- `factor()` 入口泛化：≥3 变量走任意多元路径，2 变量保留原二元路径 /
  `factor()` dispatch: ≥3 variables use the multivariate path.
- 一元非首一分解：`factor_square_free` 首项系数变换，修复
  `factor_square_free_monic` 不能分解非首一多项式的问题 /
  **Non-monic univariate factorization** via a leading-coefficient
  transformation (`factor_square_free`).
- 测试/基准：多元 correctness 往返用例、criterion 多元基准组、
  proptest 往返（标记 ignore，手动运行）/ multivariate correctness
  roundtrip cases, criterion benchmarks, and an (ignored) proptest.

### Fixed / 修复

- `div_rem_sparse` 的整除判断方向（潜在 bug，此前无调用方）/
  `div_rem_sparse` monomial-divisibility argument order (latent bug).
- 多元 Diophantine 求解的循环上界取原始误差次数（对齐 Symbolica），
  消除伪无限循环 / multivariate Diophantine loop bound uses the original
  error degree (matching Symbolica), eliminating a pseudo-infinite loop.
- `zassenhaus_combine` 重写为按 Symbolica 算法（候选乘余因子首项系数后
  取本原部分）/ `zassenhaus_combine` rewritten per Symbolica.

### Known Limitations / 已知限制

- 多元路径对**首项系数非常数且无法施加**（需模 $p$ Hensel 提升施加真
  首项系数）的输入会保守地按不可约返回；该增强计划在 0.16.1 /
  Multivariate inputs whose non-constant leading coefficient cannot be
  imposed are reported irreducible; the mod-$p$ imposition lift is planned
  for 0.16.1.

---

## [0.15.2] - 2026-07-21

### Performance / 性能

- **reducer 首项单项式哈希索引**（`ocas-poly::groebner::f4`）：基首项按
  support-mask 分桶 + 子掩码枚举的 `DivisorIndex`/`find_reducer`，消除
  符号预处理中 O(单项式 × 基) 的线性扫描；保留"n_terms 最小 + 低下标
  平局"选择语义 / **Reducer LM divisor index**: basis leading monomials
  are bucketed by exact support mask; reducer queries enumerate the
  submasks of the query support instead of scanning the whole basis —
  eliminating the O(monomials × basis) linear scan in symbolic
  preprocessing while preserving the "fewest terms, lowest index"
  selection semantics
- **稀疏行 echelon**（`echelonize_fp`/`echelonize_generic`）：消元由稠密
  buffer + 全列扫描改为稀疏双指针归并相消（`sub_scaled_fp`/
  `sub_scaled_generic`），每次相消 O(nnz) 而非 O(ncols)；scratch 行复用、
  头部相消跳过、去零 / **Sparse-row echelon**: elimination now merges
  rows with a two-pointer sparse AXPY instead of a dense buffer + full
  column scan — O(nnz) per cancellation; reused scratch row, head
  cancellation skipped, zeros dropped
- **单项式表行模板缓存**（fp 路径）：基倍式内容经 `row_cache`/
  `row_store` 缓存，避免重跑 `get_simplified` 简化缓存 / **Row-template
  cache** (fp path): basis-multiple content is cached by `(basis_idx,
  diff)`, avoiding re-running the simplification cache
- **提取阶段查重哈希化**：基首项查重由 O(基) 线性扫描改为 `HashSet`
  / **Extraction dedup hashed**: basis-LM duplicate check uses a `HashSet`
  instead of a linear scan
- **符号预处理 worklist 驱动**：新注册单项式即入队、LIFO 处理，替代
  每轮对全单项式表的 `!present` 扫描 / **Worklist-driven symbolic
  preprocessing**: newly-registered monomials are queued and processed
  LIFO instead of rescanning the whole monomial table each pass
- cyclic-6 ℤ₁₃：9970 s → **3670 s**（2.7×，basis=20 正确且通过
  `is_groebner_basis`）；阶段占比转为消除主导（echelon ≈89%）。**<5 s
  未达成**——cyclic-6 的 F4 矩阵第 22 轮达 264k 行 × 284k 列，为 F4 对
  该理想的固有规模；进一步数量级提升需 F5 签名约简（消除零约化行），
  列入 post-1.0 / cyclic-6 over ℤ₁₃: 9970 s → **3670 s** (2.7×,
  basis=20, `is_groebner_basis` pass); the phase profile shifted to
  elimination-dominated (echelon ≈89%). **<5 s not reached** — the F4
  matrix reaches 264k rows × 284k cols at round 22, intrinsic to F4 for
  this ideal; a further order-of-magnitude win needs F5 signature
  reduction (post-1.0)

---

## [0.15.1] - 2026-07-20

### Fixed / 修复

- **F4 矩阵列序**（`ocas-poly::groebner::f4`）：列改按单项式**降序**排列
  ——此前为升序，导致高斯消元在尾项上建主元，echelon 形同虚设，全部
  约化工作退化为多项式除法（F4 实为 Buchberger）/ **F4 matrix column
  order**: columns are now sorted in *descending* monomial order — they
  were ascending, so elimination pivoted on trailing terms and the echelon
  step was decorative (F4 was effectively Buchberger)
- **echelon 回写条件**：被消元行未写回、新主元行被清空的反向条件已修正 /
  **echelon write-back condition**: eliminated rows were never written back
  and new pivot rows were emptied — the inverted condition is fixed
- **Gebauer–Moeller 判据**（`update_pairs`）：移植 Symbolica 的正确实现
  ——旧实现对 lcm(i,new) == lcm(i,j) 的配对误删，cyclic-5 上基不完整 /
  **Gebauer–Moeller criteria**: ported Symbolica's correct `update` — the
  old version dropped pairs whose lcm is reproduced by (i, new), producing
  incomplete bases on cyclic-5
- **F4 提取判据**：S 部分改以两个独立倍式入行（经典 F4 形式），提取仅当
  行首项不在输入行首项集合时加入，**不再做基约化** / **F4 extraction**:
  S-parts enter as two separate multiples (classic F4); a row joins the
  basis only when its LM differs from every input row head — no basis
  reduction during extraction

### Performance / 性能

- **cyclic-5 ℤ₁₃：2609 s → 31 ms（≈85 000×）**，且首次通过
  `is_groebner_basis` 验证 / **cyclic-5 over ℤ₁₃: 2609 s → 31 ms
  (≈85,000×)**, passing `is_groebner_basis` for the first time
- cyclic-6/7 现为可解（见 `ocas-tests/tests/groebner_timing.rs`）/
  cyclic-6/7 are now tractable

---

## [0.15.0] - 2026-07-20

### Added / 新增

- **多输出 JIT**（`ocas-eval::jit`）：`ExpressionEvaluator::compile_jit`
  支持任意多输出；`JitCompiledFunction::call_into` 零分配调用；整数幂
  经平方幂链内联（|exp| ≤ 16），常量在编译期内联 / **Multi-output
  JIT**: `compile_jit` supports any number of outputs; `call_into` for
  allocation-free calls; integer powers inlined via exponentiation by
  squaring (|exp| ≤ 16); constants embedded at compile time
- **f32 混合精度**（`ocas-eval`）：`JitCompiledF32`/`compile_jit_f32`（libm
  `*f` 符号）与 `VectorEvaluatorF32`/`compile_vector_evaluator_f32`
  （16 lane）/ **f32 mixed precision**: `JitCompiledF32`/`compile_jit_f32`
  (libm `*f` symbols) and `VectorEvaluatorF32`/`compile_vector_evaluator_f32`
  (16 lanes)
- **流式求值 API**（`ocas-eval::streaming`）：`StreamingEvaluator` 分块
  输入 + 复用栈缓冲，百万行内存恒定 / **Streaming evaluation**:
  `StreamingEvaluator` with chunked input and reused stack, constant
  memory over million-row streams
- **多输出编译**（`ocas-eval::compile`）：`compile_atoms_multi` /
  `compile_trees_multi` / `ExpressionEvaluator::compile_multi`，跨输出
  共享 CSE 与常量表 / **Multi-output compilation** sharing CSE and
  constant tables across outputs
- **常量折叠与栈压缩**（`ocas-eval`）：`EvalTree::fold_constants`（树级
  代数恒等式 + 常量求值）、DCE 后栈压缩消除空洞 / **Constant folding
  & stack compaction**: tree-level identities + DCE stack compaction
- **栈复用**（`ocas-eval::evaluator`）：`evaluate_with_stack` 复用调用方
  缓冲，消除每次求值的堆分配 / **Stack reuse**: `evaluate_with_stack`
  reuses caller buffers, eliminating per-evaluation heap allocation
- **Arena reset + workspace 池**（`ocas-core`/`ocas-atom`）：
  `Arena::reset()` 保留首块复用；`ocas-atom::workspace` 线程本地
  Arena 池（RecycledAtom 模式）/ **`Arena::reset()`** (keep first chunk)
  and a thread-local workspace pool (RecycledAtom pattern)
- **ahash 哈希**（`ocas-core`）：`FastHashMap`/`FastHashSet` 别名替换
  热点 std HashMap（sparse terms、F4 列映射、hash-consing、CSE）/
  **ahash**: `FastHashMap`/`FastHashSet` replacing std HashMap on hot paths
- **原生 i64 F4 管线**（`ocas-poly::groebner::f4`）：ℤ_p 上全程 i64
  残基运算（`FpPoly`），消除 BigInt 往返；`OCAS_F4_STATS` 分段插装 /
  **Native i64 F4 pipeline**: full `i64` residue arithmetic over ℤ_p
  (`FpPoly`), eliminating BigInt round-trips; `OCAS_F4_STATS` section timing

### Fixed / 修复

- **JIT 调用约定**（`ocas-eval::jit`）：Windows 上改用 target 默认调用
  约定（`default_call_conv()`），取消 `#[ignore]` 测试 / **JIT calling
  convention**: use target-default call conv on Windows, un-ignore tests
- **JIT 多输出写回**：`result_indices` 确定性映射替代 HashMap 顺序
  遍历 / **Multi-output writeback**: deterministic `result_indices`
  mapping instead of HashMap iteration order
- **Arena 对齐**：`Chunk::try_alloc` 对齐检查修复（超对齐分配独占
  chunk）/ **Arena alignment**: `Chunk::try_alloc` alignment check
  (over-aligned allocations get their own chunk)

### Performance / 性能

- JIT 求值：多项式 1000 次调用 221 µs → 2.27 µs（**97×**），三输出
  479 µs → 22.4 µs（**21×**）/ **JIT**: 97×/21× over interpreter
- 流式：100k 行 23.4 ms → 16.8 ms（**28%↑**）/ **Streaming**: 28% faster
- F4 瓶颈定位：extract 占 2608 s/2609 s（99.98%），echelon 仅 12.87 ms
  ——cyclic-6 < 5 s 目标需 RREF/F5，推迟到 0.15.1 / **F4 bottleneck**:
  extraction dominates (99.98%); cyclic-6 < 5 s deferred to 0.15.1

### Changed / 变更

- `VectorEvaluator::evaluate` 返回 `Vec<Vec<f64>>`（多输出，API 破坏
  变更）/ `VectorEvaluator::evaluate` now returns `Vec<Vec<f64>>`
  (breaking change for multi-output support)
- workspace 版本 0.14.0 → 0.15.0

---

## [0.14.0] - 2026-07-18

### Added / 新增

- **Risch 符号积分**（`ocas-calc`）：初等超越塔上的递归积分器
  （Bronstein 第 5 章）。Hermite 约化、对数导数恒等式、primitive 待定
  系数、hyperexponential Risch 微分方程（第 6 章多项式片段）。入口
  `integrate()` 管线第三层 / **Risch symbolic integration**: recursive
  integrator over elementary transcendental towers (Bronstein ch. 5).
  Hermite reduction, log-derivative identity, primitive undetermined
  coefficients, hyperexponential RDE (ch. 6 polynomial fragment)
- **有理函数积分**（`ocas-calc::integral::rational`）：Hermite 约化 +
  对数部分（`c·f'/f`、二次配方 log/atan、Rothstein–Trager 插值 +
  ℤ 因式分解）/ **Rational-function integration**: Hermite reduction +
  logarithmic part (log-derivative identity, completing the square,
  Rothstein–Trager)
- **特殊函数积分表**（`ocas-calc::integral::special`）：非初等积分的
  闭式端点——`exp(-x²)→erf`、`exp(c·x²)→erf/erfi`、`exp(x)/x→Ei`、
  `sin(x)/x→Si`、`cos(x)/x→Ci`、`sinh/cosh→Shi/Chi`、
  `sin(x²)/cos(x²)→fresnels/fresnelc`（Meijer-G 端点，替代完整管线）/
  **Special-function antiderivatives**: closed-form endpoints for
  non-elementary integrals (erf/erfi/Ei/Si/Ci/Shi/Chi/Fresnel)
- **三角积分**（`ocas-calc::integral::trig`）：sin/cos/tan/cot/sec/csc
  重写为 `exp(I·x)` 后 Risch 积分，共轭对数对合并回实 log/atan /
  **Trigonometric integration**: rewrite to exp(I·x), Risch, realify
- **FGLM 换序**（`ocas-poly::groebner::fglm`）：零维理想 Gröbner 基
  换序（BFS 阶梯 + 高斯消元），远快于重跑 F4 / **FGLM order
  conversion** for zero-dimensional ideals
- **F5 签名判据**（`ocas-poly::groebner::f5`，实验性）/ **F5 signature
  criterion** (experimental)
- **Hilbert 界**（`ocas-poly::groebner::hilbert`）：单项式理想 Hilbert
  分子 + 正则性界 / **Hilbert bounds** for monomial ideals
- **`GroebnerBasis::reorder`**（`ocas-poly`）：换单项式序（简单路径，
  重跑 F4）/ **monomial-order conversion** (re-run F4)
- **mdBook 章节**：`algorithms/groebner.md` 与 `algorithms/integration.md`
  （中英双语）/ **mdBook chapters**: groebner.md + integration.md (en/zh)

### Fixed / 修复

- **解析器优先级**（`ocas-parse`）：`-x^2` 现在正确解析为 `-(x^2)`
  （幂优先于负号），而非 `(-x)^2` / **parser precedence**: `-x^2` now
  parses as `-(x^2)`, not `(-x)^2`
- **`RationalPolynomial::canonicalize`**（`ocas-poly`）：单变量情形现在
  执行真实多项式 GCD 约分（此前仅系数 content 约分），修复非规范形 /
  **canonicalize**: exact polynomial GCD reduction for univariate
- **0.11.0 已知差距**：`exp(-x²)→erf` 经特殊函数表闭合 /
  **0.11.0 known gap**: `exp(-x²)→erf` closed

### Changed / 变更

- **`KElem` 自动约分**（`ocas-calc::tower::elem`）：比例/单项式约分
  消除垃圾分母（`t/t`），Risch 塔算术规范化 / **KElem auto-reduction**:
  ratio/monomial cancellation
- **积分管线**：`integrate()` = 启发式 → 有理 → Risch → 三角 → 特殊
  函数 → `Integral(...)` / **integration pipeline**: heuristic → rational
  → Risch → trig → special → unevaluated
- **测试**：correctness 套件新增 `integral_risch.rs`（15 项，与 SymPy
  对比）；积分基准 `benches/integrate.rs` / **tests**: integral_risch
  correctness suite (15 SymPy-checked) + integrate benchmark

---

## [0.13.2] - 2026-07-18

### Added / 新增

- **PyPI 发布**（`ocas-py`）：`pip install ocas` 正式上线。通过 OIDC trusted
  publishing 发布，覆盖 5 个平台 wheel（macOS x86_64/arm64、Linux
  x86_64/aarch64、Windows x86_64）+ sdist / **PyPI 发布**：`pip install ocas`
  正式上线，OIDC trusted publishing，5 平台 wheel + sdist
- **macOS wheels**：补齐三平台 wheel 覆盖（0.9.0 路线图承诺），macOS Intel
  wheel 在 Apple Silicon runner 上交叉编译产出 / **macOS wheels**：补齐三平台
  覆盖，Intel wheel 经交叉编译产出

### Fixed / 修复

- **CI lint**：清理 `clippy --all-targets`（含 bench/test 目标）下的冗余闭包与
  未使用 import；纯 Rust CI job 移除 `system-libs` 特性（需系统 MPC ≥ 1.4.1，
  Ubuntu 24.04 仅有 1.3.1）/ **CI lint**：清理 `--all-targets` lint；纯 Rust job
  移除 `system-libs`
- **Security audit**：固定 stable 工具链（cargo-audit 依赖需 Rust ≥ 1.96）/
  **安全审计**：固定 stable 工具链
- **ocas-c 发布**：build.rs 改为写入 `OUT_DIR` 而非源目录，修复 `cargo publish`
  "Source directory was modified by build.rs" 错误 / **ocas-c 发布**：build.rs
  写入 `OUT_DIR`
- **Trusted publishing**：移除 `release.yml` 对 `wheels.yml` 的 reusable
  workflow 调用（PyPI trusted publishing 不支持 reusable workflow，OIDC
  `workflow_ref` 会错位导致 `invalid-publisher`）/ **Trusted publishing**：
  移除 reusable workflow 调用
- **Release tag 校验**：预发布 tag 允许 base 版本 ≥ workspace 版本 /
  **Release tag 校验**：放宽预发布 tag

### Changed / 变更

- **依赖升级**：cranelift 0.117→0.127（锁步组）、chumsky 0.10→0.13、logos
  0.15→0.16、cbindgen 0.28→0.29、criterion 0.5→0.8、hashbrown 0.15→0.17、
  flint3-sys 3.5→3.6、egg 0.10→0.11 / **依赖升级**：cranelift、chumsky、logos、
  cbindgen、criterion、hashbrown、flint3-sys、egg
- **安全修复**：crossbeam-epoch 0.9.18→0.9.20（RUSTSEC-2026-0204，`fmt::Pointer`
  无效指针解引用）/ **安全修复**：crossbeam-epoch（RUSTSEC-2026-0204）
- **dependabot 降噪**：cranelift 锁步组 + patch/actions 批量组；忽略
  dtolnay/rust-toolchain 误升级 / **dependabot 降噪**：分组 + 忽略误升级

## [0.13.1] - 2026-07-06

### Fixed / 修复

- **docs.rs build**: avoid system C library dependencies in documentation
  builds. The hosted docs are now built with portable features only (no
  `gmp`, `mpfr`, `flint`, `python`, `gpl`). Users who need the full backend API
  docs can build locally with `cargo doc -p ocas --features gmp,mpfr,flint
  --no-deps` / **docs.rs 文档构建**：避免在文档构建中依赖系统 C 库。托管文档
  现在仅使用可移植特性构建（不含 `gmp`、`mpfr`、`flint`、`python`、
  `gpl`）。需要完整后端 API 文档的用户可本地运行 `cargo doc -p ocas
  --features gmp,mpfr,flint --no-deps`。

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
  **mimalloc 全局分配器**：可选 `mimalloc` 特性配置全局分配器
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

- `gmp` 特性在 Windows MSVC 上不受支持，因为 `rug` 无法在该环境下构建
  GMP。请改用 MSYS2/MINGW64 或 Linux/macOS。