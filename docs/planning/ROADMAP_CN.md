# oCAS 路线图

本文档概述 oCAS 从 pre-alpha 实验到稳定 1.0 版本的开发路线图，每个 0.x 版本
都包含具体的交付物。英文版见 [ROADMAP_EN.md](ROADMAP_EN.md)。配套文档：
[EVOLUTION_PLAN_CN.md](EVOLUTION_PLAN_CN.md)（细粒度逐版本计划）与
[GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md)（差距快照）。

---

## 图例

| 标签 | 含义 |
|---|---|
| `API` | 公共 API 表面 |
| `ALG` | 代数算法 |
| `NUM` | 数值后端 |
| `PERF` | 性能与优化 |
| `BIND` | 语言绑定 |
| `DOC` | 文档与示例 |
| `TEST` | 测试与质量 |

---

## 阶段 1：Pre-Alpha — 基础

> **目标**：建立 workspace、运行时与基础表达式核心，证明分层架构可以编译并运行。

### 0.1.0 — Workspace 与运行时

**目标日期**：第 1 个月

**交付物**：

- [x] 包含全部 12 个 crate 的 workspace 结构
- [x] CI 流水线：`cargo test`、`cargo clippy`、`cargo-deny`、格式化、Miri
- [x] 统一错误类型 `OcasError`
- [x] 通过 Miri 安全验证的 arena / bump 分配器
- [x] 基于 `rayon` 的线程池包装
- [x] FFI 胶水约定（最小 C ABI 示例）
- [x] `gmp` feature 后的 GMP 绑定（基于 `rug`）
- [x] 初始基准测试框架

**成功标准**：

- `cargo build --workspace` 在 Linux/macOS/Windows 上成功（MSVC 上无默认特性）。
- Arena 通过 Miri 与 valgrind/ASan 检查。
- 在受支持平台上可从 Rust 调用 GMP 整数运算。

### 0.2.0 — 表达式树核心

**目标日期**：第 2 个月

**交付物**：

- [x] `ocas-atom` crate
- [x] `Atom` 标签联合设计
- [x] 带安全公共 API 的 arena 后端 AST
- [x] 公共子表达式 hash consing
- [x] 基于 `logos` 的词法分析器
- [x] 递归下降 / Pratt 语法分析器
- [x] 打印器：ASCII 与紧凑形式
- [x] 规范化器：展平 `Add`/`Mul`、排序项、合并系数

**成功标准**：

- `parse("x^2 + 2*x + 1")` 产生预期的 AST。
- 对大量表达式满足 `to_string(parse(s)) == s`。
- 规范化具有确定性且通过属性测试。

---

## 阶段 2：Alpha — 符号引擎

> **目标**：提供可用的 Rust API，支持解析、化简、微分与基础多项式运算。

### 0.3.0 — 域与多项式

**目标日期**：第 4 个月

**交付物**：

- [x] `ocas-domain` crate
- [x] 域实现：`Integer`、`Rational`、`FiniteField`
- [x] 泛型算法的 Domain trait
- [x] `ocas-poly` crate
- [x] 稠密单变量多项式
- [x] 域实现：`RealBall`、`Complex`
- [x] 稀疏多元多项式
- [x] 带余除法
- [x] `flint` feature 后的 FLINT 3 集成
- [x] 通过 `rug` 提供的可选 GMP `Integer`/`Rational` 后端
- [x] 通过 `rug` 提供的可选 MPFR `RealBall` 后端

  > **说明**：该特性为实验性。在提供系统 FLINT 的 Linux/WSL 上可构建运行，
  > 但目前尚不支持 Windows，因为 `flint3-sys` 依赖 `pthread_mutex_t` 等仅
  > POSIX 的类型。Windows 上默认推荐的大整数、有理数与严格实数后端为通过
  > MSYS2 安装系统 GMP/MPFR 后使用 `rug` 的 `gmp`/`mpfr` 特性。

**成功标准**：

- 多项式运算在回归套件上与 SymPy 输出一致。
- FLINT 路径在受支持运算上与纯 Rust 回退产生相同结果。

### 0.4.0 — 模式匹配与重写

**目标日期**：第 5 个月

**交付物**：

- [x] 带通配符与条件的模式匹配引擎
- [x] `Transformer` 访问者 API
- [x] 基础内置重写规则
- [x] `egg` 等式饱和集成
- [x] 基于规则的化简器

**成功标准**：

- 常见恒等式（如 `x + x -> 2*x`、`x * 0 -> 0`）自动应用。
- E-graph 在假设下可将 `sin(x)^2 + cos(x)^2` 化简为 `1`。

### 0.5.0 — 微积分基础

**目标日期**：第 6 个月

**交付物**：

- [x] 符号微分
- [x] 初等函数导数表
- [x] Taylor 级数展开
- [x] 基于启发式表的部分积分
- [x] `ocas-calc` crate 初始版本

**成功标准**：

- 微分通过综合测试套件。
- 积分在标准微积分问题上成功。

### 0.6.0 — 首个 Rust API 候选版

**目标日期**：第 7 个月

**交付物**：

- [x] 稳定的 `ocas` prelude
- [x] 所有公共 API 的 rustdoc 示例
- [x] `proptest` 属性测试
- [x] 初始基准测试套件
- [x] 通过 `uv` 的 SymPy 对比基准
- [x] crates.io 发布准备（内部工作区依赖已版本化）

**成功标准**：

- `cargo test --workspace --exclude ocas-py` 通过。
- 基准在基础多项式、微积分与重写运算上展示与 SymPy 持平。
- `cargo publish --dry-run -p ocas-core` 成功；内部 crate 上传后顶层 `ocas` 即可发布。

---

## 阶段 3：Beta — 求解器、JIT、绑定

> **目标**：实现多语言可用性与性能，核心代数功能达到 CAS beta 标准。

### 0.7.0 — 方程求解器

**目标日期**：第 9 个月

**交付物**：

- [x] 线性方程组求解器（`faer` / `LinBox`）
- [x] 多项式方程组求解器（Gröbner + 根隔离）
- [x] 基于 Arb 的单变量求根
- [x] 丢番图方程基础
- [x] 假设/域系统

**成功标准**：

- 线性与多项式求解器产生经 SageMath 验证的正确结果。

### 0.8.0 — 求值与 JIT

**目标日期**：第 11 个月

**交付物**：

- [x] 标量与向量求值的树解释器
- [x] AST 到指令序列编译器
- [x] 用户自定义函数注册表
- [x] Cranelift JIT 后端
- [x] SIMD 向量化求值

**成功标准**：

- JIT 求值重复表达式比解释器至少快 10 倍。
- SIMD 路径对稠密多项式求值有效。

### 0.9.0 — Python 与 C/C++ 绑定

**目标日期**：第 13 个月

**交付物**：

- [x] 基于 PyO3 的 `ocas-py` crate
- [~] Python 类：`Expression`（完成），`Polynomial`/`Matrix`/`Domain`（推迟到 0.10.0）
- [x] Linux/macOS/Windows 的 Maturin 轮子构建
- [x] 基于 cbindgen 的 `ocas-c` crate
- [x] 稳定的表达式生命周期 C API
- [x] C++ RAII 包装

**成功标准**：

- `pip install ocas` 在受支持平台上可用。
- C 示例针对共享库编译并运行。
- 绑定测试无内存泄漏（tracemalloc + RAII 守护的 arena）。

### 0.10.0 — Beta 发布

**目标日期**：第 14 个月

**交付物**：

- [x] 从 0.9.0 推迟的 Python 类：`Polynomial`、`Matrix`、`Domain`
- [x] 1.0 功能冻结
- [x] 综合文档站点
- [x] 与 Symbolica/SageMath 的性能对比
- [x] 社区反馈整合
- [x] 仅修复 bug 阶段

**成功标准**：

- 所有公共 API 均有文档。
- CI 在所有受支持平台上绿灯。

---

## 阶段 4：竞品差距弥合（0.24–0.26）

> **目标**：弥合本次竞品调研发现的关键差距（GAP_ANALYSIS_CN.md §5），
> 在 1.0.0 冻结前将 P0–P2 缺口降至可接受水平。
> 
> 背景：阶段 B++ “竞品全面对齐”（0.19–0.23）已于 2026-08-02 完成。
> 但竞品在此期间有重大演进——Symbolica 2.2 移植 Rubi 7000+ 积分规则、
> SymPy 1.14 DomainMatrix 10000× 加速、msolve Gröbner 性能标杆
> cyclic-6 仅 0.04 s——导致原有“1.0 仅做冻结”的计划不再充分。
> 新增阶段 B+++ 三个版本弥合差距，然后进入 1.0.0 冻结。
>
> **（本阶段已完成：0.24 启发式积分四技术 + DoubleF64、0.25 MultiModular
> Gröbner + 并行模 GCD、0.26 打包单项式 F5 快通道；cyclic-6 ℤ₁₃ grevlex
> 实测 55.04 ms。0.26.0 实际交付与本节原计划不同——矩阵引擎/Smith 标准形
> 顺延至阶段 5 的 0.35.0。）**

### 0.24.0 — 符号积分广度 + DoubleFloat

**目标**：缩小与 Symbolica Rubi 的积分覆盖面差距；引入 DoubleFloat 求值路径。

**交付物**：

- [x] 积分启发式扩展：Risch 回退后的 `heuristic_integrate` 池
  - 分部积分（LIATE/ILATE 启发式）
  - 三角替换（$\sqrt{a^2 - x^2}$、$\sqrt{a^2 + x^2}$、$\sqrt{x^2 - a^2}$）
  - 有理参数替换（Weierstrass $t = \tan(x/2)$）
  - Euler 代换（二次根式下的有理化，占位）
  - 参考：SymPy `manualintegrate` 启发式池
- [x] DoubleFloat 求值路径（`DoubleF64`：~31 位，>3× 快于任意精度）
  - 参考：Symbolica 2.0 `double-float` 实现
  - 在 `ocas-domain` 新增 `DoubleFloat` 类型
  - JIT/SIMD 求值器支持 DoubleFloat 管线
- [x] Python/C 绑定：`integrate_heuristic`、`DoubleFloat` 类型
- [ ] 基准：Rubi 1892 题子集对标 symbolica-integrate（推迟到 0.27.0）

**成功标准**：
- Rubi 1892 题子集覆盖率从当前水平提升 ≥30%（从 Risch-only 到 Risch + 启发式）
- DoubleFloat 求值比任意精度快 ≥3×
- `cargo test --workspace` 通过

### 0.25.0 — Gröbner 大规模性能（Multi-Modular）

**目标**：对齐 msolve 的 Gröbner 性能，cyclic-6 ℤ₁₃ 从 2.63 s 降至 < 0.5 s。

**交付物**：

- [x] 多模算术（multi-modular）策略
  - 多个素数并行计算 Gröbner 基
  - 中国剩余定理（CRT）重建整数系数基
  - 有理重构（rational reconstruction）恢复 $\mathbb{Q}$ 系数
  - 参考：msolve F4 + multi-modular + Hensel + BM
- [x] Hensel 提升 Gröbner 基
  - 从 $\mathbb{F}_p$ 基提升到 $\mathbb{Z}$ 基
  - 减少 CRT 重建的素数数量
- [x] 大系数多项式 GCD 加速
  - Brown 模 GCD 利用 multi-modular 进一步加速
- [ ] 基准：cyclic-6/7、katsura-6/7 对标 msolve（katsura 系推迟到 0.33.0）

**成功标准**：
- cyclic-6 ℤ₁₃ < 0.5 s（当前 2.63 s，msolve 0.04 s）
- cyclic-7 ℤ₁₃ 可解（当前未测）
- 基准结果与 msolve 在同一数量级（< 10× 差距）

### 0.26.0 — 打包 F5 快通道 + grevlex 性能（实际交付）

**目标**：将 F5 主循环压入 u128 SWAR 快通道，进一步对齐 msolve 性能；
补充 grevlex 基准变体。（原计划的域感知矩阵引擎 + Smith/Hermite 标准形
未在 0.26.0 交付，顺延至 0.35.0。）

**交付物**：

- [x] 打包单项式 F5 快通道（u128 SWAR）
  - n_vars ≤ 8 且指数 < 2¹⁵ 自动路由；超界回落通用路径
- [x] echelon i32 / 免克隆两阶段改造
- [x] grevlex 基准变体（补 Lex 之外的实测基线）
- [x] 修复 Graded 序度方向反置的预存在 bug
- [ ] 域感知矩阵引擎（`DomainMatrix` 类似物）→ 0.35.0
- [ ] Smith/Hermite 标准形 → 0.35.0
- [ ] 矩阵性能基准 → 0.35.0
- [ ] 1.0 冻结前准备（API 审计/迁移指南/跨平台 CI）→ 0.35.0

**成功标准**（实测，2026-08-06）：
- cyclic-6 ℤ₁₃ grevlex 52.07 ms（criterion 中位数）、Lex 936 ms
- cyclic-7 ℤ₁₃ grevlex 单轮 5.755 s（209 基元素）
- 打包快通道与通用路径结果一致（随机基准交叉验证）

---

## 阶段 5：竞品差距收尾与机制攻坚（0.27–0.35）

> **目标**：依据 2026-08-06 竞品复测的优先级重排（GAP_ANALYSIS_CN.md §5），
> 在 1.0.0 冻结前收尾剩余 P0–P3 缺口：P0 符号积分广度、P1 Gröbner
> 大规模性能（katsura 系 + cyclic-7）、P1 LLVM JIT 代码生成、P2 矩阵/
> 线性代数（DomainMatrix 类似物 + Smith/Hermite 标准形）、P2 Windows
> FLINT、P3 二次筛与张量嵌套函数内处理。阶段 B+++（0.24–0.26）已交付
> 启发式积分/DoubleF64、MultiModular Gröbner、打包 F5 快通道
> （cyclic-6 grevlex 55.04 ms）；本阶段收尾其余差距后进入 1.0.0 冻结。
>
> **2026-09-12 重排（0.27.3 收尾后）**：0.27.x 的失败归因调研
> （BENCHMARK_RESULTS_CN.md §「0.27.3 后续调研」）与通用机制可行性评估
> （GENERAL_MECHANISM_FEASIBILITY_CN.md）把 P0 的**性质**重新定义了：
> 0.27.0 的原始验收线（「Rubi 级规则集 + 1892 题 +30pp」）是**规则枚举**
> 指标，而 0.27.x 六波机制攻坚的边际收益已从 +6.82pp 衰减到 +1.11pp，
> 且剩余 1522 例中仅 0.7% 与已解题目共享形状骨架。因此：
>
> - **0.28.0–0.32.0 改为机制攻坚**：先修两个已实测的架构性缺陷并建立
>   「符号证书 + 三值输出」的正确性纪律，再按 Bronstein 阶梯补超越 Risch、
>   代数扩张、非初等层与复杂度；
> - **原 0.28.0 / 0.29.0 / 0.30.0 依序顺延为 0.33.0 / 0.34.0 / 0.35.0**
>   （Gröbner 大规模性能、LLVM JIT、矩阵引擎 + 平台收尾 + 1.0 冻结准备）；
> - 1.0.0 顺延至第 69 月。

### 0.27.0 — 符号积分广度（Rubi 级规则集）

**目标**：弥合与 `symbolica-integrate`（Rubi 7000+ 规则、72,944 题库）的
最大功能缺口（P0），1892 题子集覆盖率显著提升。

**交付物**：

- [x] 规则表驱动的积分规则引擎（match → 模板替换）
  - 幂/多项式/指数/对数规则族
  - 三角/双曲/反三角/反双曲规则族
  - 根式与二次型代换（在 0.24 三角换元/Weierstrass/Euler 框架上扩展，
    补齐 Euler 占位）
  - 特殊函数规则族（erf/Ei/Si/Ci/Fresnel，衔接 0.14 函数表）
- [x] 策略调度链：Risch（0.14）→ 启发式四技术（0.24）→ 规则库 →
  `Integral(...)` 回退
- [x] 规则来源策略（参考 GAP_ANALYSIS_CN.md §7.3 许可证风险分析）：
  - 首选自研规则集（方案 C 混合：Risch + 启发式 + 规则结构参考 Rubi 分类）
  - 评估 `symbolica-integrate`（MIT）作为可选 feature 的集成可行性
- [x] 1892 题覆盖率基准 harness：覆盖率报告 + 失败分类分析
- [x] Python/C 绑定：`integrate` 规则路径开关
- [x] （(c) 阶段增补）符号常数有理积分器、Weierstrass 线性变元、
  有界分配展开重试、三角积化和差/降幂、subresultant 稠密 GCD、
  积分链全局条目预算、Wilkinson 精确实根隔离（10/10）

**成功标准**（诚实记录）：

- 1892 题子集覆盖率从当前水平提升 ≥30 个百分点 —— **未达成**：
  基线 5.87% → (b) 阶段 7.66% → (c) 阶段见 BENCHMARK_RESULTS_CN.md；
  根因量化（Rubi 级规则量级/符号系数有理后端/嵌套根式均需更大投入），
  备选路径评估见 GAP_ANALYSIS §7.3
- 规则路径与 SymPy `manualintegrate`/`integrate` 抽样交叉验证一致 —— 达成
- `cargo test --workspace` 通过 —— 达成

### 0.27.1 — 积分广度机制攻坚（0.27.0 验收线续作）

**目标**：沿 0.27.0 的 +30pp 验收线继续提升 1892 题覆盖率，以机制级
升级（非规则堆量）解锁整族题型。

**交付物**（全部落地）：

- [x] Chebyshev 二项微分 + 分数幂有理化（`binomial.rs`）
- [x] 三角分母幂递推 / 线性分子分解 / 多项式×三角闭式（`trig_reduction.rs`）
- [x] exp/log 核代换（exp 核有理化、双曲 t=e^u、f(log x)/x）（`exp_log.rs`）
- [x] 一般二次根式引擎 + Euler III（`sqrt_quadratic.rs`）
- [x] 反三角/反双曲核导数幂与代换（`inverse_trig.rs`）
- [x] 单三角核有理式归一化 + tan/sec 族剥项（`trig_kernel.rs`）
- [x] 分母幂递推 + 双线性因子部分分式（`quad_power.rs`）
- [x] 错案修复：C14/D7b 线性变元递推残项系数、rational.rs √(p/q) 丢 1/q
- [x] 稳定性：symbolic_rational 系数预算与多符号闸门（超时 49→33、
  墙钟 −19%）、heuristic 深度残项检查、expand 前移防预算饥饿

**成功标准**（诚实记录）：

- 1892 题子集 +30pp —— **未达成**：9.62% → 16.44%（+6.82pp，311/1892，
  129 新解、0 回归、0 崩溃；量化缺口：椭圆族/高符号商式/复合壳层，
  见 BENCHMARK_RESULTS_CN.md 0.27.1 段）
- 全部新增机制经 eval_f64 数值求导抽样核验 + SymPy 对拍全绿 —— 达成
- 质量门（fmt/clippy -D warnings/workspace test/deny）全绿 —— 达成

### 0.27.2 — 悬挂消除、已验证覆盖率与椭圆积分基础

**目标**：把每一例残留超时变成确定性拒绝；在字符串覆盖率之外增加独立的
数值验证口径；补齐失败转储指向的初等机制族；落地椭圆积分基础。

**交付物**（全部落地）：

- [x] 诊断：`OCAS_INTEGRATE_TRACE` 阶段打点；33 例基线超时全部归因到阶段
  （symbolic_rational 18、heuristic 4、trig_kernel 3、inverse_trig 2、
  rational 2、trig_reduction 2、sqrt_quadratic 1、未打点 1）
- [x] 符号有理后端前的有界展开预通道（`integral/mod.rs` 的 `expand_prepass`）：
  展开为小型 Laurent 多项式的未展开发散积不再进入该后端（仅限纯代数展开，
  避免抢走三角/双曲机制应有的题目）
- [x] 悬挂阶段（`symbolic_rational`、`trig_kernel`、`heuristic`、`rational`、
  `sqrt_quadratic`）的确定性工作量预算
- [x] 错案修复：`symbolic_rational::rational_square_root` 把多项式**和**当作
  单项式平方，用伪根拆分二次分母并输出错误对数（0.27.1 因此对
  `1/(b*x^2+2*a*x-b)` 以及经 `t = e^x` 双曲路径的 `1/(a+b*sinh(x))`
  给出错误答案）；另修 `complete_square` 首一假定、Chebyshev 回代丢符号、
  Risch 非域项、规则 A4 序列通配符、共振积化和差零分母共六组错案
- [x] 数值验证 oracle（`ocas-tests/src/integral_eval.rs`）与 1892 harness 的
  `verified_solved` / `unverified_solved` / `verify_indeterminate` 口径，
  以及 `OCAS_1892_VERIFY` / `OCAS_1892_BUCKET` 开关
- [x] 双曲闭式族（`hyperbolic_reduction.rs`）
- [x] 有理导数核代换（`kernel_subst.rs`）
- [x] 三角相位归一（`trig_reduction.rs` 扩展；相位/比值路径已实现但以
  `PHASE_RATIO_ENABLED = false` 门控，因组合原子尚未通过 `crate::diff`）
- [x] 逆函数复合消去（`inverse_trig.rs` 扩展）
- [x] `exp(逆函数)` 代数化（`exp_log.rs` 扩展）
- [x] 半幂前端（`halfpower.rs`）与椭圆约化（`elliptic.rs`）：Legendre 约化到
  `EllipticF`/`EllipticE`（`EllipticPi` 头已注册但尚无生产者），采用 SymPy 的
  `m = k²` 约定，实参顺序由 `ocas_atom::normalize::preserves_argument_order` 保持
- [x] 错案回归护栏（正确性套件 `tests/correctness/integral_verify.rs`）
- [x] 额外落地：`ocas-parse` 一元负号（`2-1`/`x-1`/`x^2-1` 不再 `PARSE_ERR`，
  破坏性 `lex` API 变更）

**成功标准**（诚实记录）：

- 1892 题双口径：**solved 349（18.45%）、verified 325/349（93.1%）**；
  逐题 diff 新解 44、主动回退 6（全部是 0.27.1 的错案）、净 +38
- 已解题目零错案：**`verify_mismatches` = 0 —— 达成**
- 超时 13（目标 ≤5 **未达**）、崩溃 0、墙钟 445.7 s（较 572 s 基线 −22%）
- 已验证比例 93.1%（目标 ≥95% **未达**：24 例 inconclusive 中 11 例定义域受限、
  10 例含虚数单位、3 例步长自洽闸门判为不可判；均为「不可判」而非错案）
- 质量门：fmt/clippy 两档/deny 全绿；workspace test 全绿（`ocas-c` 的
  规则开关探针已按 kernel_subst 接手 `tan(x)^4` 的事实更新为 `csc(x)^5`）
- 下一波（0.27.3）承接：13 例超时族、复合壳层拆分、特殊函数族、椭圆族广度

### 0.27.3 — 复合壳层、特殊函数广度与椭圆族覆盖

**目标**：沿 0.27.2 失败转储留存的最大簇继续机制攻坚。

**交付物**（实测结果逐条记录）：

- [x] 前置：特殊函数导数表（`derivative.rs`）与数值 oracle 扩头
  （`Ei`/`Ei(n,z)`/`Si`/`Ci`/`Shi`/`Chi`/`fresnels`/`fresnelc`）。所有算法先与
  `mpmath` 40 位对拍再落码；`Ei(n,z)` 的未实现偏导（阶数槽）一律回退为未求值
  `Derivative`，绝不静默丢项
- [x] 特殊函数族扩展（`integral/special.rs`）：四个归约族——多项式 × `F(a+bx)`、
  多项式 × `Ei(n,a+bx)`、`F(bx)/x^m` 递降、多项式 × `F(a+bx)²`；步数预算
  `MAX_SPECIAL_STEPS=16`、`MAX_SPECIAL_DEG=8`，残项全部闭式，不留 `Integral(...)`。
  **19 例 special 桶：15 解（全部经数值验证）、4 诚实拒绝**（两例 Rubi 自身
  `Unintegrable`，两例 Fresnel 分母链刻意不实现）
- [x] 半幂前端仿射变元（`halfpower.rs`）：接受 `cos(c + d·x)`，两支都发射片层因子
  （`cos(u/2)·(1−sin²(u/2))^(−1/2)`），并把同基的 `S^a·(√S)^b` 先折叠为单一半幂。
  **新解 4 例**（`rubi-00377/00445/01598/00291`；其中 4 例在赛题 oracle 判为
  inconclusive，因为其参数固定为 `m=2`，超出 oracle 的 `|m| ≤ 0.9` 窗口，模块自带
  独立数值 oracle 在 1e-9 上核验通过）
- [x] 精确线性平方折叠（`integral/mod.rs`）：`p²+2pq+q² → (p+q)²`，**仅当底对积分变量
  仿射**、且重新展开能精确复现原和时接受；绝不在非整数幂之下生效
  （`((p+q)²)^{1/2} = |p+q| ≠ p+q`）。**新解 2 例**，并把 `rubi-00854`
  从 10 s 超时变为 0.01 s 求解
- [ ] 复合壳层拆分（mixed-other 簇）：**部分达成 / 部分撤回**。因子级
  `F(g(x))·w(x)` 拆分的具体实现「线性分式对数归约」（`C·(A+B·log Q)·(F+Gx)^p`
  分部）在集成后暴露一个真错案（残项未包 `Integral(...)` 直接返回）且残项积分不可靠、
  墙钟代价高，按「只保留无回归且无墙钟代价的代码」原则**整体撤回**；`log_fraction.rs`
  已删除。撤回过程反而定位并加固了 A4 回归护栏（由「必须留残项」改为「若求解则数值求导
  核验」，两种结局都不放过错案）。仿射变元泛化（`trig_linear_arg` 的方向）由上一项交付
- [ ] 椭圆族广度：三角二次根式全量路由、三次根式、`EllipticPi` 复特征值 ——
  **未达成，如实记录**。0.27.3 只交付「仿射变元」这一层；140 例仿射 `cos` 半幂簇中
  仍有 136 例拒绝，根因已量化：它们需要**两个不同根式基的积/商**、`cos u` 的有理前因子、
  `|p| > 5` 或正幂 ≥ 3，即一个多根式有理前因子归约引擎（新引擎，非本波范围）。
  三次根式与 `EllipticPi` 复特征值未动

**成功标准**（诚实记录）：

- 1892 题双口径：**solved 370（19.56%）、verified 343/370（92.7%）、mismatches 0**；
  逐题 diff **新解 21、回归 0**（桶 delta：special +15、trig +3、power-binomial +2、
  mixed-other +1）
- 已解题目零错案：**`verify_mismatches` = 0 —— 达成**（这是本波硬性红线，也是
  0.27.3 撤回 `log_fraction` 的直接原因）
- 超时 13 → **12**（目标 ≤5 **未达**）、崩溃 0（**达成**）、
  墙钟 461.2 s（目标 < 445.7 s **未达**，+3.5%：新增归约族在 `special` 阶段的扫描与
  精确折叠的前置遍历是主要来源）
- 已验证比例 92.7%（目标 ≥95% **未达**：27 例 inconclusive 全部是「不可判」而非错案，
  其中新增 4 例是赛题 oracle 的 `|m| ≤ 0.9` 窗口所致）
- 质量门：fmt / clippy 两档 / deny / workspace test 全绿（`ocas-calc` 416 passed）
- **0.27 线冻结判定**：按 `EVOLUTION_PLAN` 的「连续两波净增 < 60 题则冻结」字面口径，
  0.27.2（+38）与 0.27.3（+21）**连续两波各自净增均 < 60** → **冻结 0.27 线，
  下一活动线为 0.28.0**。0.27.3 是该判定下的 0.27 系列最后一版
- 下一波（0.28.0）承接：机制攻坚（正确性地基 → 超越 Risch → 代数扩张 →
  非初等层 → 复杂度），并把本波量化但未解决的 12 例超时族、136 例多根式半幂簇、
  椭圆族广度、残项解析与循环检测两个缺陷纳入其中；原 0.28.0（Gröbner）
  顺延为 0.33.0（见 §5 阶段说明与 GENERAL_MECHANISM_FEASIBILITY_CN.md）

### 0.28.0 — 积分机制正确性地基（缺陷修复 + 证书化 + 正则塔）

**目标**：把 0.27.3 调研定位的两个**架构性缺陷**修掉，把「正确」从抽样验证
升级为**符号证书**，并放开通用 Risch 引擎最贵的入口限制（依赖生成元被拒）。
（依据：GENERAL_MECHANISM_FEASIBILITY_CN.md §3、§5、§6 P0）

**交付物**：

- [ ] **残项解析**（缺陷）：`rational::integral_fallback` 与 `symbolic_rational`
  高阶因子分支产生的 `Integral(...)` 残项，改为**链尾 + 确定性预算**的解析
  （原型实测净 +4：+5 新解 / −1 回归 / +3 新超时，
  见 BENCHMARK_RESULTS_CN.md §「0.27.3 后续调研」§3.1）
- [ ] **循环检测改为表达式级**（缺陷），替换/收窄 `MAX_CHAIN_ENTRIES`
  全局总量上限（`rubi-00008` 实测 306 次阶段进入 / 18 轮循环）
- [ ] **符号证书**：每个输出必须通过微分域内的 `normalize(D(F) − f) == 0`；
  数值 oracle 降级为探测/回归工具
- [ ] **三值 Outcome**：`Found { value, certificate }` /
  `ProvedNonElementary { witness }` / `Unknown`；`Unknown` 诚实返回 `Integral(...)`
- [ ] **正则塔**：合并代数相关生成元（`log x`/`log 2x`、`exp x`/`exp(x+1)`、
  双曲重写产生的 `exp(±u)`），相关但不可判定时返回 `Unknown`
- [ ] 指标：harness 新增 `certified_rate = 证书通过 / 已解`，CI 门禁要求恒为 1.0

**成功标准**：

- 1892 已解集合的符号证书 **100% 为 0**；`certified_rate = 1.0`
- 双曲族新增解（基线：111 例未解含双曲函数）
- `verify_mismatches` 保持 0；逐题 diff **0 回归**；超时数不上升

### 0.29.0 — 超越 Risch 补全

**目标**：把 `integral/rde.rs` 从「仅多项式解」扩到完整片段，并补对数部分的
结构定理。（依据：GENERAL_MECHANISM_FEASIBILITY_CN.md §3.4、§6 P1）

**交付物**：

- [ ] RDE **有理解**（分母界 + `D`-有理解，Bronstein ch. 6 完整版）
- [ ] **耦合微分系统**（`D y + A y = b`）
- [ ] 对数部分的**一般结构定理**（不再只依赖对数导数恒等式）
- [ ] 通用性测试：proptest **从塔**随机生成元素 → 积分 → 符号证书校验

**成功标准**：

- `exp-log` 桶（基线 79/83 未解）显著改善；证书门保持
- 随机塔元素族的证书通过率 100%；1892 逐题 diff 0 回归

### 0.30.0 — 代数扩张与反函数代换

**目标**：接上代数积分链（复用已有 Trager 资产），并覆盖反函数族。
（依据：GENERAL_MECHANISM_FEASIBILITY_CN.md §3.2、§3.5、§6 P2）

**交付物**：

- [ ] 代数生成元进入塔（`√x` 与一般根式；当前 `tower/build.rs` 直接拒绝）
- [ ] **积分基**（Trager）+ **代数 Hermite 约化** + 代数留数（Bronstein ch. 7–8）
- [ ] 多根式基归约（估算 222 例含 ≥2 个不同根式基）
- [ ] 反函数代换引擎：多项式权 × `(a+b·f(ax+b))^k`
  （`inverse-trig-hyper` 桶基线 138/146 未解）
- [ ] 复用 `ocas-poly` 的 Trager 因式分解 + 结式 + 代数数域 GCD

**成功标准**：

- `1/(1+x⁴)`、`1/(1−3x²+x⁴)`、不可约六次/八次分母族可解
- radical 桶（基线 356/407 未解）显著改善；证书门保持

### 0.31.0 — 非初等层

**目标**：补上今天**完全不可达**的那一类：需要二重对数/Meijer G 的原函数
（基线：122 例参考答案含 `polylog`，占未解 8.0%）。
（依据：GENERAL_MECHANISM_FEASIBILITY_CN.md §3.6、§6 P3）

**交付物**：

- [ ] `polylog`/`Li₂` 函数头 + 数值 oracle 核验
  （先补头与数值核验、再补归约——沿用 0.27.3 对 `Ei`/`Si`/`Ci`/`Shi`/`Chi`/Fresnel 的落地模式）
- [ ] Li₂ 归约（对数-有理积分）
- [ ] Meijer G / Slater 展开（或 holonomic/D-finite 路线）作为一般非初等机制
- [ ] `ProvedNonElementary` 首次可用（Singer 结构定理的工程化子集）

**成功标准**：

- 122 例 `polylog` 题中可及部分解锁；`exp-log` 桶剩余改善
- 新增函数头在 oracle 中数值核验通过；证书门保持

### 0.32.0 — 复杂度与性能（把复杂度当算法问题）

**目标**：符号系数下的系数爆炸用**算法**解决，而不是继续加墙钟预算
（基线：12 例超时中 9 例挂在 `symbolic_rational`）。
（依据：GENERAL_MECHANISM_FEASIBILITY_CN.md §3.7、§6 P4）

**交付物**：

- [ ] 符号系数下的模算法（复用 0.25 的 MultiModular ℚ 管线）
- [ ] 并行 Risch
- [ ] 惰性级数 / 度界剪枝，替换「加预算」式缓解
- [ ] 12 例残存超时的阶段级处理（`symbolic_rational` 9、`binomial` 1、
  `rational` 1、`trig_kernel` 1，见 `timeout_attribution_0273.csv`）

**成功标准**：

- 墙钟不高于 0.28.0 基线（0.27.3 为 461.2 s）；超时数下降
- 高层数/多符号实例不再出现系数爆炸式回退；证书门保持

### 0.33.0 — Gröbner 大规模性能（katsura 系 + cyclic-7）（原 0.28.0 顺延）

**目标**：对齐 msolve 0.10.1 实测（katsura 3–7 ms、cyclic-7 55 ms）（P1），
katsura-6 < 1 s、cyclic-7 grevlex 进入同数量级。

**交付物**：

- [ ] u128 打包 F5 快通道扩展到 katsura 系与 cyclic-7（指数域/稀疏度适配）
- [ ] MultiModular ℚ 管线（0.25）扩展到大规模实例
  - 并行幸运素数调度 + CRT + 有理重构 + 无迹 p-adic Hensel 提升
- [ ] echelon 稀疏度感知优化（0.15.2 稀疏 echelon 的后续：行/列剪枝）
- [ ] katsura-6/7、cyclic-7 grevlex/Lex 基准对标 msolve 实测（WSL2）

**成功标准**：

- katsura-6 ℤ₁₃ < 1 s（当前未完成）；katsura-7 可完成
- cyclic-7 grevlex 与 msolve 差距 < 10×（当前 ~70×）
- 多模路径与单素数路径随机 100 例一致；`is_groebner_basis` 验证

### 0.34.0 — 代码生成扩展（LLVM/inkwell JIT）（原 0.29.0 顺延）

**目标**：落地第二个 JIT 后端——LLVM（经 `inkwell`，已在 workspace 依赖），
缩小与 Symbolica SymJIT 的代码生成差距（P1）。

**交付物**：

- [ ] `ocas-eval::jit_llvm`：AST → LLVM IR + 函数注册表 + 多输出
- [ ] 求值管线覆盖：f64/f32 混合精度 + DoubleF64 + SIMD 向量化
- [ ] 运行时后端选择：Cranelift（默认，编译快）/ LLVM（优化代码）
- [ ] 性能基准：LLVM vs Cranelift vs 解释器（保持多输出 97×/21× 基线）
- [ ] Python/C 绑定暴露后端选择参数

**成功标准**：

- LLVM JIT 与 Cranelift 持平或更优；相对解释器 ≥10× 保持
- 三平台（Linux/macOS/Windows）LLVM 构建 CI 全绿
- 与 Cranelift 路径输出一致（随机 1000 表达式）

### 0.35.0 — 矩阵引擎 + 平台收尾 + 1.0 冻结准备（原 0.30.0 顺延）

**目标**：收尾 P2/P3 差距并完成 1.0 冻结前准备：域感知矩阵引擎
（DomainMatrix 类似物）+ Smith/Hermite 标准形、Windows FLINT、二次筛、
张量嵌套函数内处理。

**交付物**：

- [ ] 域感知矩阵引擎（`DomainMatrix` 类似物，从 0.26.0 顺延）
  - `Matrix<D>` 泛型化：支持 `IntegerDomain`、`FiniteField`、`RationalDomain`
  - Dense 矩阵的域特化路径（避免通用 `Domain` trait 开销）
  - 参考：SymPy DomainMatrix + FLINT 后端
- [ ] Smith 标准形（整数矩阵，用于模结构分析与同调代数）
- [ ] Hermite 标准形（整数矩阵，用于线性丢番图方程）
- [ ] 矩阵性能基准：20×20/30×30 整数矩阵 rref/inv/det 对标 SymPy DomainMatrix
- [ ] Windows FLINT 支持（flint3-sys Windows 构建评估 + CI）
- [ ] 二次筛大整数分解（对标 SymPy `qs_factor`，ECM 之上的下一级）
- [ ] 张量嵌套函数内处理（对标 Symbolica Graphica；0.22 已交付基础规范化）
- [ ] 1.0 冻结前准备
  - API 审计：所有公共类型/函数的文档完整性
  - 迁移指南定稿（从 Symbolica/SymPy 迁移到 oCAS）
  - 跨平台 CI 验证（Linux/macOS/Windows）
  - 已发布基准（基于 BENCHMARK_SUITE_CN.md）

**成功标准**：

- Smith/Hermite 标准形与 SymPy 交叉验证一致（随机 100 例）
- 20×20 整数矩阵 rref 性能与 SymPy DomainMatrix 在同一数量级
- 二次筛基准与 SymPy `qs_factor` 对比记录
- Windows FLINT 三平台可用（或记录明确的技术阻塞）
- 1.0 冻结前准备清单完成 ≥80%

---

## 阶段 6：稳定 1.0

> **目标**：发布 API 稳定、后端支持广泛的成熟 CAS 库。

### 1.0.0 — 稳定发布

**目标日期**：0.35.0 之后

**交付物**：

- [ ] 稳定语义化版本保证
- [ ] 完整的 Rust、Python 与 C/C++ API 覆盖
- [ ] 综合测试套件（行覆盖率 >80%）
- [ ] 已发布基准测试（基于 BENCHMARK_SUITE_CN.md）
- [ ] 从 Symbolica/SymPy 迁移指南
- [ ] 签名发布产物
- [ ] 竞品对标报告（基于 COMPETITIVE_MATRIX_CN.md 最终版）

**成功标准**：

- 1.x 期间无计划中的破坏性 API 变更。
- P0 差距（符号积分广度）已显著缩小（1892 题子集覆盖率目标达成，
  GAP_ANALYSIS_CN.md §5）。
- P1 差距（Gröbner 性能）已对齐 msolve 至同一数量级（katsura-6 < 1 s、
  cyclic-7 grevlex < 10× msolve）。
- P1 差距（代码生成）已落地 LLVM/inkwell JIT 后端。
- P2 差距（矩阵/线性代数）已落地域感知矩阵引擎 + Smith/Hermite 标准形。
- 在核心基准上性能全面领先 SymPy。

> 细粒度逐版本计划详见 [EVOLUTION_PLAN_CN.md](EVOLUTION_PLAN_CN.md)。
> 阶段 A（Beta 硬代数 0.11–0.13）、阶段 B+（Symbolica 差距清零 0.15.2–0.18.1）、
> 阶段 B++（竞品全面对齐 0.19–0.23）、阶段 B+++（竞品差距弥合 0.24–0.26）
> 均已完成。
> 阶段 B++++（竞品差距收尾 0.27–0.30）收尾本次调研剩余的 P0–P3 差距。

---

## 1.0 之后

1.0 之后，开发重点将转向：

- 偏微分方程（PDE）求解器（Poisson、热传导、波动）
- 微分 Galois 理论（研究序章）
- 可选 GPL 后端（`ocas-gpl`）
- CUDA/WASM 代码导出（对标 Symbolica SymJIT CUDA/WASM）
- 领域专用工具包（物理、机器人、机器学习）

> LLVM/Inkwell JIT 已前置至 0.34.0，二次筛整数分解已前置至 0.35.0。

---

## 里程碑

| 版本 | 阶段 | 目标日期 | 关键交付物 |
|---|---|---|---|
| 0.1.0 | Pre-Alpha | 第 1 个月 | 工作空间 + 运行时 |
| 0.2.0 | Pre-Alpha | 第 2 个月 | 表达式核心 |
| 0.3.0 | Alpha | 第 4 个月 | 域与多项式 |
| 0.4.0 | Alpha | 第 5 个月 | 模式匹配与重写 |
| 0.5.0 | Alpha | 第 6 个月 | 微积分基础 |
| 0.6.0 | Alpha | 第 7 个月 | Rust API 候选版 |
| 0.7.0 | Beta | 第 9 个月 | 方程求解器 |
| 0.8.0 | Beta | 第 11 个月 | JIT 与求值 |
| 0.9.0 | Beta | 第 13 个月 | Python 与 C/C++ 绑定 |
| 0.10.0 | Beta | 第 14 个月 | 功能冻结 |
| 0.11.0 | Beta | 第 15 月 | 多项式因式分解（一元） |
| 0.11.1 | Beta | 第 15 月 | 多项式因式分解（二元 + 绑定 + 文档） |
| 0.11.2 | Beta | 第 16 月 | 计算加速基础设施（SOO Integer、mimalloc、模方法 GCD） |
| 0.12.0 | Beta | 第 17 月 | 有理多项式 + 结式 + 部分分式 + Karatsuba 乘法 + 有理重构 |
| 0.13.0 | Beta | 第 19 月 | Gröbner F4 矩阵化算法 |
| 0.13.1 | Beta | 第 19 月 | docs.rs 构建修复 |
| 0.13.2 | Beta | 第 19 月 | PyPI 发布（`pip install ocas`）+ 依赖升级 + CI 加固 |
| 0.14.0 | 1.0 候选 | 第 22 月 | Risch 符号积分 + 有理函数积分 + 特殊函数表 + FGLM/F5/Hilbert + 三角积分 |
| 0.15.0 | 1.0 候选 | 第 24 月 | 多输出 JIT + f32 混合精度 + 流式求值 + Arena/workspace 池 + ahash + 原生 i64 F4 |
| 0.15.1 | 1.0 候选 | 第 24 月 | F4 真实线性代数修复（cyclic-5 提速 ≈85 000×，cyclic-6 可解） |
| 0.15.2 | 1.0 候选 | 第 25 月 | Gröbner 大规模性能（LM 索引 + 稀疏 echelon，cyclic-6 ℤ₁₃ 9970 s → 3670 s；<5 s 需 F5） |
| 0.16.0 | 1.0 候选 | 第 26 月 | 任意多元因式分解（Wang EEZ，≥3 变量，ℤ 与 ℤ_p）✅ |
| 0.16.1 | 1.0 候选 | 第 26 月 | 非常数首项系数强加（模 p Hensel）+ 多元稀疏化 + 稀疏 Diophantine | ✅ |
| 0.16.2 | 1.0 候选 | 第 26 月 | $\mathbb{F}_p$ 路径非常数 LC 预处理（域版 Wang）+ 采样性能优化 |
| 0.17.0 | 1.0 候选 | 第 27 月 | 代数数域与扩域因式分解（Trager 算法）✅（一元路径；多元扩域留待后续） |
| 0.18.0 | 1.0 候选 | 第 28 月 | 数值积分（Vegas）+ 双数自动微分 + 张量基础 + fuel 资源控制 |
| 0.18.1 | 1.0 候选 | 第 28 月 | 0.18.0 三项能力的 Python/C 绑定补齐（数值积分 + 张量 + 双数）+ prelude 补齐 ✅ |
| 0.19.0 | 1.0 候选 | 第 30 月 | F5 Gröbner 基签名约简（cyclic-6 ℤ₁₃ <5 s 目标）✅（2.63 s，≈1400×；多序推迟到 0.19.1） |
| 0.20.0 | 1.0 候选 | 第 33 月 | 常微分方程求解器（一阶 5 种 + 二阶 2 种 + 幂级数框架 + 分类引擎）✅（核心完成；Laplace/系统/绑定推迟） |
| 0.20.1 | 1.0 候选 | 第 33 月 | ODE 补齐：积分因子 + 常数变易法 + 降阶法 + 级数递推 + Frobenius + Laplace IVP + 2×2 系统 + Python/C 绑定 + 31 项代入验证测试 ✅ |
| 0.21.0 | 1.0 候选 | 第 36 月 | 数论与计算代数（模 GCD + 整数分解 + 素性 + 离散对数 + CRT + 数论函数）✅（另含 Python/C 绑定；ECM 30 位半素数 1.1 s） |
| 0.22.0 | 1.0 候选 | 第 39 月 | 张量规范化（图同构引擎）+ 高级模式匹配（`Transformer::Partition`）✅ |
| 0.23.0 | 1.0 候选 | 第 42 月 | 高级 Gröbner 与代数几何工具（理想运算 + 准素分解 + Hilbert 级数）✅ |
| 0.24.0 | Beta | 第 45 月 | 符号积分广度（启发式扩展）+ DoubleFloat 求值路径（P0 积分 + P2 DoubleFloat）✅ |
| 0.25.0 | Beta | 第 47 月 | Gröbner 大规模性能（multi-modular 对标 msolve，cyclic-6 < 0.5 s）（P1）✅ |
| 0.26.0 | Beta | 第 49 月 | 打包单项式 F5 快通道 + grevlex 基准（cyclic-6 grevlex 55.04 ms 实测）✅ |
| 0.27.0 | Beta | 第 51 月 | 符号积分广度（Rubi 级规则集 + 1892 题覆盖率基准）（P0） |
| 0.28.0 | Beta | 第 53 月 | 积分机制正确性地基（残项解析 + 表达式级循环检测 + 符号证书 + 三值输出 + 正则塔） |
| 0.29.0 | Beta | 第 55 月 | 超越 Risch 补全（RDE 有理解 + 耦合系统 + 对数部分结构定理） |
| 0.30.0 | Beta | 第 57 月 | 代数扩张与反函数代换（积分基 + 代数 Hermite + 留数 + 多根式基 + 反函数代换引擎） |
| 0.31.0 | Beta | 第 59 月 | 非初等层（`polylog`/`Li₂` 头 + Li₂/Meijer G 归约 + `ProvedNonElementary`） |
| 0.32.0 | Beta | 第 61 月 | 复杂度与性能（模算法 + 并行 Risch + 惰性级数） |
| 0.33.0 | Beta | 第 63 月 | Gröbner 大规模性能（katsura-6 < 1 s、cyclic-7 同数量级 msolve）（P1；原 0.28.0 顺延） |
| 0.34.0 | Beta | 第 65 月 | 代码生成扩展（LLVM/inkwell JIT 后端）（P1；原 0.29.0 顺延） |
| 0.35.0 | Beta | 第 67 月 | 矩阵引擎（DomainMatrix 类似物 + Smith/Hermite）+ Windows FLINT + 二次筛 + 1.0 冻结准备（P2/P3；原 0.30.0 顺延） |
| 1.0.0 | Stable | 第 69 月 | 稳定版发布（阶段 B++++ 竞品差距收尾完成后冻结：P0 积分机制正确性与通用性达成 + P1 Gröbner 对齐 msolve + LLVM JIT 落地 + 性能全面领先 SymPy） |

---

## 如何阅读本路线图

- 每个版本代表一个**可发布**的增量。
- 日期为预估值，取决于贡献者可用时间。
- 功能可能根据用户反馈与技术发现在不同版本间调整。

---

## 参与路线图

如果你想参与某个特定版本或功能，请创建 GitHub issue，我们会为你分配跟踪 issue。
