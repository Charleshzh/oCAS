# 差距分析：oCAS 与参考系统

本文档逐里程碑（0.1 → 1.0+）跟踪 oCAS 的实现完成度，并对照三大参考系统
评估差距：**Symbolica**（Rust）、**SageMath**（Python 生态）、**SymPy**
（纯 Python）。本文档为活文档，每次版本发布时必须更新。英文版见
[GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md)。

> 最后评估：**0.13.1 @ 2026-07-17**

---

## 图例

| 标记 | 含义 |
|---|---|
| ✅ | 完成 |
| 🟡 | 基础可用或部分完成 |
| 🔴 | 缺失或重大缺口 |
| ⚠️ | 完成但有保留 |

---

## 1. 版本完成状态（0.1–0.10）

| 版本 | 阶段 | 路线图 | 核验状态 |
|---|---|---|---|
| 0.1.0 | Pre-Alpha | ✅ | ✅ 12-crate workspace、CI、`OcasError`、arena（Miri 感知）、rayon 线程池、FFI 胶水、`gmp` feature（via `rug`） |
| 0.2.0 | Pre-Alpha | ✅ | ✅ `ocas-atom`、`Atom` 标签联合、arena AST、hash consing、logos 词法、Pratt 解析、规范化器 |
| 0.3.0 | Alpha | ✅ | ⚠️ `Integer/Rational/FiniteField/RealBall/Complex`；稠密/稀疏多项式、带余除法；`flint` 仅 Linux/WSL，GMP/MPFR 经 `rug` |
| 0.4.0 | Alpha | ✅ | ✅ matcher、pattern、rules、simplify、transformer、`egraph.rs`（egg 集成） |
| 0.5.0 | Alpha | ✅ | ⚠️ 微分、积分（启发式）、Taylor 级数；积分为查表式，无 Risch |
| 0.6.0 | Alpha | ✅ | ✅ 稳定的 `ocas` prelude、rustdoc 示例、proptest、criterion、SymPy harness、crates.io 准备 |
| 0.7.0 | Beta | ✅ | ⚠️ 线性（有理/整数）+ 丢番图 + 多项式组（Gröbner）；Sturm 根隔离；assumptions——算法偏基础 |
| 0.8.0 | Beta | ✅ | ✅ 树解释器、AST→指令编译器、函数注册表、Cranelift JIT、SIMD 向量化求值 |
| 0.9.0 | Beta | ✅ | ⚠️ PyO3 `Expression`/`Evaluator`/`solve_*`；cbindgen + C++ RAII 包装——部分类推迟到 0.10 |
| 0.10.0 | Beta | ✅ | ✅ Python `Polynomial/Matrix/Domain`、Matrix 线性代数（Bareiss）、mdBook 文档站、三平台 wheels CI、版本锁定 0.10.0 |
| 0.11.0 | Beta | ✅ | ✅ 完整多项式因式分解（ℤ 与 ℤ_p：Yun SFF → CZ → Hensel → Zassenhaus）、多元 GCD、500 例 proptest 往返测试、版本提升至 0.11.0 |
| 0.11.1 | Beta | ✅ | ✅ 二元因式分解（ℤ 与 ℤ_p：关于 x 首一的 Wang Hensel）、稀疏多元 `factor()` 入口、C 多项式绑定、mdBook 因式分解章节、版本提升至 0.11.1 |
| 0.12.0 | Beta | ✅ | ✅ 有理多项式 `RationalPolynomial<D,O>`、Brown PRS 结式、Karatsuba 快乘、扩展 GCD、多项式 CRT/丢番图、p-adic 展开、部分分式分解、有理重构、版本提升至 0.12.0 |
| 0.12.1 | Beta | ✅ | ✅ 自研 ℤ_p 上 NTT、`pulp` SIMD 分派、Estrin 多项式求值、F4 稀疏矩阵后端、数值验证特性、版本提升至 0.12.1 |
| 0.13.0 | Beta | ✅ | ✅ F4 Gröbner 基算法（含 Gebauer-Moeller 临界对筛选与简化缓存）、`Grlex` 单项式序、`Domain` trait 扩展、`FiniteField` ℤ_p 快速路径工具、版本提升至 0.13.0 |

0.1–0.13.0 交付物全部落地，workspace 版本锁定 0.13.0。质量门全绿：
`cargo fmt`、`clippy -D warnings`、workspace 测试、`cargo deny`、pytest、
`mdbook build`。

---

## 2. 代码规模

`src/` 下 Rust 行数快照（不含测试与基准）。

| Crate | 文件数 | 行数 |
|---|---|---|
| ocas-poly | 10 | ~4,250 |
| ocas-eval | 11 | ~2,525 |
| ocas-domain | 9 | ~2,115 |
| ocas-rewrite | 7 | ~1,719 |
| ocas-py | 7 | ~1,546 |
| ocas-calc | 7 | ~1,393 |
| ocas-c | 4 | ~1,550 |
| ocas-core | 5 | ~1,150 |
| ocas-atom | 2 | ~864 |
| ocas-parse | 3 | ~565 |
| ocas (prelude) | 1 | ~113 |
| ocas-gpl | 1 | 1（占位） |
| **src 合计** | **66** | **~18k** |

`ocas-gpl` 为占位；GPL 专属后端属 Post-1.0 工作，符合路线图。

---

## 3. 算法深度核验

本节是决定 CAS 成熟度最关键的因素，也是差距的主要来源。

| 算法领域 | oCAS 现状 | 成熟度 |
|---|---|---|
| 多项式因式分解 | `DenseUnivariatePolynomial` 上 ℤ 与 ℤ_p 的 `factor()`，以及 `SparseMultivariatePolynomial` 上二元 ℤ 与 ℤ_p 的 `factor()`（关于 x 首一的 Wang Hensel） | 🟢 较完整 |
| Gröbner 基 | F4 矩阵化算法（Faugère 1999）+ Gebauer-Moeller + 简化缓存 + ℤ_p 快速路径 | 🟢 F4 完成 |
| 符号积分 | 启发式查表（幂/逆/sin/cos/exp/线性替换）；回退为 `Integral(...)`；**无** Risch | 🟡 基础 |
| 实根隔离 | Sturm 序列 + 区间隔离 + refine（单变量） | 🟢 较完整 |
| 多项式 GCD | GCD + 本原部分；无模 GCD / EEA 优化 | 🟡 可用 |
| 线性求解 | 有理/整数线性方程组 + 二元丢番图（`ax+by=c`） | 🟡 可用，规模有限 |
| JIT 求值 | Cranelift 后端；≥10x 加速目标达成（按路线图标准） | 🟢 完整 |

---

## 4. 与参考系统的差距

### 4.1 对照 Symbolica（Rust，AGPL）

Symbolica 的 `examples/` 目录揭示了成熟度差距。oCAS 大致相当于 Symbolica
早期的功能子集。

| 能力 | oCAS | Symbolica |
|---|---|---|
| 多项式因式分解 | ✅ ℤ 与 ℤ_p 上 `factor()`（CZ + Hensel + Zassenhaus）；二元 ℤ 与 ℤ_p 上因式分解（关于 x 首一的 Wang Hensel） | ✅ 完整（`factorization.rs`） |
| 有理多项式 | ✅ 含 GCD 规范化的 `RationalPolynomial<D,O>` | ✅ `rational_polynomial.rs` |
| 部分分式 | ✅ 任意 `EuclideanDomain` 上的 `apart()` / `together()` | ✅ `partial_fraction.rs` |
| 有理重构 | ✅ 基于扩展欧几里得的 `rational_reconstruction(a, m)` | ✅ `rational_reconstruction.rs` |
| 数值积分 | 🔴 无 | ✅ `numerical_integration.rs` |
| 流式 API | 🔴 无 | ✅ `streaming.rs` |
| 张量 / 双数 | 🔴 无 | ✅ `tensors.rs` / `dual.rs` |
| 优化 / 代码生成 | 🟡 JIT，仅 f64 | ✅ `optimize.rs` / 多输出 |
| Gröbner 基 | � F4 完成 | ✅ 工业级 |

Symbolica 的核心竞争力——工业级因式分解、有理函数运算、多输出优化、流式
API——oCAS 基本缺失。Symbolica 经多年打磨，oCAS 需在 ALG 层补齐硬算法。

### 4.2 对照 SageMath（Python 生态）

SageMath 是"瑞士军刀"式科学计算环境，差距是**广度级**的。

| 领域 | oCAS | SageMath |
|---|---|---|
| 代数几何 | 🟡 基础 Gröbner | ✅ Singular 集成 |
| 数论 | 🟡 基础丢番图 | ✅ PARI/FLINT 全栈 |
| 微分方程 | 🔴 无 | ✅ 完整 ODE/PDE 求解器 |
| 群论/表示论 | 🔴 无 | ✅ GAP 集成 |
| 组合数学 | 🔴 无 | ✅ 完整 |
| 绘图/可视化 | 🔴 无 | ✅ matplotlib 集成 |
| 数据库接口 | 🔴 无 | ✅ OEIS / LMFDB |

SageMath 通过包装 80+ 专用库实现广度；oCAS 是自包含内核。两者定位不同——
oCAS 面向高性能**库**，SageMath 面向完整**环境**。可比性主要集中在核心
代数性能，而非功能广度。

### 4.3 对照 SymPy（纯 Python）

SymPy 是 oCAS 最现实的"功能对标 + 性能超越"目标。

| 领域 | oCAS vs SymPy | 说明 |
|---|---|---|
| 解析/化简 | 🟢 持平 | 双方都完备 |
| 微分 | 🟢 持平 | 链式/乘积/幂法则 |
| 积分 | 🟡 oCAS 较弱 | SymPy 有 Risch + 启发式；oCAS 仅启发式 |
| 因式分解 | � 持平 | 单变量 ℤ 与 ℤ_p 已通过 CZ + Hensel + Zassenhaus；多元推迟到 0.11.1 |
| Gröbner | 🟡 oCAS 略弱 | 双方都非顶级，SymPy 略丰富 |
| 矩阵/线性代数 | 🟢 持平 | oCAS 有 Bareiss 行列式/逆 |
| **性能** | 🟢 **oCAS 优势** | Rust + Cranelift JIT + arena 对纯 Python |
| Python 易用性 | 🟢 持平 | oCAS 有 `ocas-py` 绑定 |

0.6.0 成功标准"基础多项式/微积分/重写与 SymPy 持平"——在**性能**维度已
达成并领先，单变量**因式分解**已实现持平；**积分**仍是与 SymPy 硬算法差距
中的主要短板。

---

## 5. 关键缺口与优先级

按"影响面 × 实现成本"排序，通往 1.0 的硬骨头。

| # | 缺口 | 优先级 |
|---|---|---|
| 1 | ~~完整多项式因式分解~~（0.11.0–0.11.1 完成） | ✅ 已完成——一元与二元（关于 x 首一）闭合，解阻塞有理函数、部分分式、求解器 |
| 2 | Risch 符号积分（路线图：0.14） | 🔴 "能否积分"的标志 |
| 3 | Gröbner F4/F5（路线图：0.13） | � F4 核心完成（0.13.0），F5 推迟 |
| 4 | ~~有理多项式/部分分式~~（0.12 完成） | ✅ 已完成——`RationalPolynomial` 类型 + 部分分式 + 结式 + Karatsuba 乘法 |
| 5 | 多输出优化/代码生成 | 🟡 JIT 为 f64 单输出；扩展到多输出/多精度 |
| 6 | ODE/PDE 求解器（Post-1.0） | 🟢 用户期望高 |

---

## 6. 总评

0.1 → 0.12 执行质量很高：每个路线图交付物均兑现，分层架构干净（无环依赖），
12 crate workspace 严格分层，质量门严格（`-D warnings` + deny + Miri 意识），
文档/绑定/CI 工程化完备。

0.12 完成了有理函数运算栈（`RationalPolynomial` 类型 + 四则运算 + 部分分式
+ 结式 + Karatsuba 乘法 + 有理重构），弥补了 GAP_ANALYSIS 中标记的三大
🔴 缺口。至此 oCAS 在有理函数能力上与 Symbolica 持平（单变量层面）。

剩余硬算法：Risch 积分（0.14）与 Gröbner F4（0.13）是通往 1.0 的最后
两个"成人礼"。

务实定位：当前 oCAS 更接近"高性能 SymPy 核心子集 + 优于 SymPy 的求值性能
+ 因式分解与有理函数持平 + Karatsuba 加速"。1.0 发布前补齐 F4 与 Risch，
是性价比最高的跃迁路径。

---

## 更新日志

每次更新在此记录（版本、日期、评估人、变更点）。

| 版本 | 日期 | 变更 |
|---|---|---|
| 0.10.0 | 2026-07-02 | 初始评估。0.1–0.10 交付物核验完成；记录与 Symbolica / SageMath / SymPy 的差距；因式分解 + Risch 积分列为最高优先级。 |
| 0.11.0 | 2026-07-03 | 多项式因式分解完成（单变量 ℤ 与 ℤ_p）；多元 GCD 加入；与 SymPy 的因式分解对比更新为持平；最高优先级缺口转为 0.12 有理函数/部分分式。 |
| 0.11.1 | 2026-07-04 | 新增二元 ℤ 与 ℤ_p 因式分解（关于 x 首一的 Wang Hensel）；稀疏多元 `factor()` 入口与 C 多项式绑定落地；新增 mdBook 因式分解章节；最高优先级缺口仍为 0.12 有理函数/部分分式。 |
| 0.12.0 | 2026-07-04 | 有理函数运算栈完成（`RationalPolynomial` + 部分分式 + Brown PRS 结式 + Karatsuba 乘法 + 有理重构）；与 Symbolica 有理函数能力持平；最高优先级缺口转为 0.13 Gröbner F4 与 0.14 Risch 积分。 |
| 0.13.0 | 2026-07-06 | Gröbner F4 矩阵化算法完成（Faugère 1999）；Gebauer-Moeller 临界对筛选 + 简化缓存 + ℤ_p 快速路径；`minimize()` bug 修复；Gröbner 从 🟡 升级为 🟢；最高优先级缺口转为 0.14 Risch 积分。 |
| 0.13.1 | 2026-07-17 | 补丁发布：docs.rs 构建改为纯 Rust 特性（不含 gmp/mpfr/flint/python/gpl），托管文档恢复构建；功能与算法层面与 0.13.0 一致，差距结论不变。 |
