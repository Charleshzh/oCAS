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
> 顺延至阶段 5 的 0.30.0。）**

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
- [ ] 基准：cyclic-6/7、katsura-6/7 对标 msolve（katsura 系推迟到 0.28.0）

**成功标准**：
- cyclic-6 ℤ₁₃ < 0.5 s（当前 2.63 s，msolve 0.04 s）
- cyclic-7 ℤ₁₃ 可解（当前未测）
- 基准结果与 msolve 在同一数量级（< 10× 差距）

### 0.26.0 — 打包 F5 快通道 + grevlex 性能（实际交付）

**目标**：将 F5 主循环压入 u128 SWAR 快通道，进一步对齐 msolve 性能；
补充 grevlex 基准变体。（原计划的域感知矩阵引擎 + Smith/Hermite 标准形
未在 0.26.0 交付，顺延至 0.30.0。）

**交付物**：

- [x] 打包单项式 F5 快通道（u128 SWAR）
  - n_vars ≤ 8 且指数 < 2¹⁵ 自动路由；超界回落通用路径
- [x] echelon i32 / 免克隆两阶段改造
- [x] grevlex 基准变体（补 Lex 之外的实测基线）
- [x] 修复 Graded 序度方向反置的预存在 bug
- [ ] 域感知矩阵引擎（`DomainMatrix` 类似物）→ 0.30.0
- [ ] Smith/Hermite 标准形 → 0.30.0
- [ ] 矩阵性能基准 → 0.30.0
- [ ] 1.0 冻结前准备（API 审计/迁移指南/跨平台 CI）→ 0.30.0

**成功标准**（实测，2026-08-06）：
- cyclic-6 ℤ₁₃ grevlex 52.07 ms（criterion 中位数）、Lex 936 ms
- cyclic-7 ℤ₁₃ grevlex 单轮 5.755 s（209 基元素）
- 打包快通道与通用路径结果一致（随机基准交叉验证）

---

## 阶段 5：竞品差距收尾（0.27–0.30）

> **目标**：依据 2026-08-06 竞品复测的优先级重排（GAP_ANALYSIS_CN.md §5），
> 在 1.0.0 冻结前收尾剩余 P0–P3 缺口：P0 符号积分广度、P1 Gröbner
> 大规模性能（katsura 系 + cyclic-7）、P1 LLVM JIT 代码生成、P2 矩阵/
> 线性代数（DomainMatrix 类似物 + Smith/Hermite 标准形）、P2 Windows
> FLINT、P3 二次筛与张量嵌套函数内处理。阶段 B+++（0.24–0.26）已交付
> 启发式积分/DoubleF64、MultiModular Gröbner、打包 F5 快通道
> （cyclic-6 grevlex 55.04 ms）；本阶段收尾其余差距后进入 1.0.0 冻结。

### 0.27.0 — 符号积分广度（Rubi 级规则集）

**目标**：弥合与 `symbolica-integrate`（Rubi 7000+ 规则、72,944 题库）的
最大功能缺口（P0），1892 题子集覆盖率显著提升。

**交付物**：

- [ ] 规则表驱动的积分规则引擎（match → 模板替换）
  - 幂/多项式/指数/对数规则族
  - 三角/双曲/反三角/反双曲规则族
  - 根式与二次型代换（在 0.24 三角换元/Weierstrass/Euler 框架上扩展，
    补齐 Euler 占位）
  - 特殊函数规则族（erf/Ei/Si/Ci/Fresnel，衔接 0.14 函数表）
- [ ] 策略调度链：Risch（0.14）→ 启发式四技术（0.24）→ 规则库 →
  `Integral(...)` 回退
- [ ] 规则来源策略（参考 GAP_ANALYSIS_CN.md §7.3 许可证风险分析）：
  - 首选自研规则集（方案 C 混合：Risch + 启发式 + 规则结构参考 Rubi 分类）
  - 评估 `symbolica-integrate`（MIT）作为可选 feature 的集成可行性
- [ ] 1892 题覆盖率基准 harness：覆盖率报告 + 失败分类分析
- [ ] Python/C 绑定：`integrate` 规则路径开关

**成功标准**：

- 1892 题子集覆盖率从当前水平提升 ≥30 个百分点
- 规则路径与 SymPy `manualintegrate`/`integrate` 抽样交叉验证一致
- `cargo test --workspace` 通过

### 0.28.0 — Gröbner 大规模性能（katsura 系 + cyclic-7）

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

### 0.29.0 — 代码生成扩展（LLVM/inkwell JIT）

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

### 0.30.0 — 矩阵引擎 + 平台收尾 + 1.0 冻结准备

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

**目标日期**：0.30.0 之后

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

> LLVM/Inkwell JIT 已前置至 0.29.0，二次筛整数分解已前置至 0.30.0。

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
| 0.28.0 | Beta | 第 53 月 | Gröbner 大规模性能（katsura-6 < 1 s、cyclic-7 同数量级 msolve）（P1） |
| 0.29.0 | Beta | 第 55 月 | 代码生成扩展（LLVM/inkwell JIT 后端）（P1） |
| 0.30.0 | Beta | 第 57 月 | 矩阵引擎（DomainMatrix 类似物 + Smith/Hermite）+ Windows FLINT + 二次筛 + 1.0 冻结准备（P2/P3） |
| 1.0.0 | Stable | 第 59 月 | 稳定版发布（阶段 B++++ 竞品差距收尾完成后冻结：P0 积分广度达成 + P1 Gröbner 对齐 msolve + LLVM JIT 落地 + 性能全面领先 SymPy） |

---

## 如何阅读本路线图

- 每个版本代表一个**可发布**的增量。
- 日期为预估值，取决于贡献者可用时间。
- 功能可能根据用户反馈与技术发现在不同版本间调整。

---

## 参与路线图

如果你想参与某个特定版本或功能，请创建 GitHub issue，我们会为你分配跟踪 issue。
