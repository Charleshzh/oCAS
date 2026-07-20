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

## 阶段 4：稳定 1.0

> **目标**：发布 API 稳定、后端支持广泛的成熟 CAS 库。

### 1.0.0 — 稳定发布

**目标日期**：第 16 个月

**交付物**：

- [ ] 稳定语义化版本保证
- [ ] 完整的 Rust、Python 与 C/C++ API 覆盖
- [ ] 综合测试套件（行覆盖率 >80%）
- [ ] 已发布基准测试
- [ ] 从 Symbolica/SymPy 迁移指南
- [ ] 签名发布产物

**成功标准**：

- 1.x 期间无计划中的破坏性 API 变更。
- 在核心基准上与 Symbolica 持平或更优。

> 从 Beta 到 1.0 的细粒度逐版本计划（0.11 因式分解 → 0.12 有理函数 →
> 0.13 Gröbner F4 → 0.14 Risch 积分 → 0.15 多输出 JIT）详见
> [EVOLUTION_PLAN_CN.md](EVOLUTION_PLAN_CN.md)。

---

## 1.0 之后

1.0 之后，开发重点将转向：

- 高级符号积分（Risch 算法）
- 微分方程求解器
- 可选 GPL 后端（`ocas-gpl`）
- GPU 加速（CUDA / HIP / Vulkan compute）
- LLVM/Inkwell JIT 后端
- 领域专用工具包（物理、机器人、机器学习）

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
| 1.0.0 | Stable | 第 26 月 | 稳定版发布 |

---

## 如何阅读本路线图

- 每个版本代表一个**可发布**的增量。
- 日期为预估值，取决于贡献者可用时间。
- 功能可能根据用户反馈与技术发现在不同版本间调整。

---

## 参与路线图

如果你想参与某个特定版本或功能，请创建 GitHub issue，我们会为你分配跟踪 issue。
