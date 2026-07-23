# 差距分析：oCAS 与参考系统

本文档逐里程碑（0.1 → 1.0+）跟踪 oCAS 的实现完成度，并对照三大参考系统
评估差距：**Symbolica**（Rust）、**SageMath**（Python 生态）、**SymPy**
（纯 Python）。本文档为活文档，每次版本发布时必须更新。英文版见
[GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md)。

> 最后评估：**0.19.1 @ 2026-07-23**（0.19.1 MonomialOrder trait 重构 + WeightOrder/BlockOrder 发布；多序支持从 `[~]` 升级为 `[x]`）

---

## 图例

| 标记 | 含义 |
|---|---|
| ✅ | 完成 |
| 🟡 | 基础可用或部分完成 |
| 🔴 | 缺失或重大缺口 |
| ⚠️ | 完成但有保留 |

---

## 1. 版本完成状态（0.1–0.18.1）

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
| 0.14.0 | 1.0 候选 | ✅ | ✅ Risch 符号积分（Hermite 约化、对数导数恒等式、primitive 待定系数、hyperexponential RDE）、有理函数积分（Hermite + Rothstein–Trager）、特殊函数表（erf/Ei/Si/Ci/Fresnel）、三角积分（exp(I·x) + realify）、FGLM/F5/Hilbert、`reorder`、mdBook 双章节 |
| 0.15.0 | 1.0 候选 | ✅ | ✅ 多输出 JIT（97×/21×）、f32 混合精度（JIT + SIMD 16 lane）、流式求值（百万行恒定内存）、常量折叠 + 栈压缩、Arena reset + workspace 池、ahash 热点替换、原生 i64 F4 管线；cyclic-6 <5s 推迟到 0.15.1（需 RREF/F5） |
| 0.15.1 | 1.0 候选 | ✅ | ✅ F4 真实线性代数修复：矩阵列序降序 + echelon 回写条件 + Symbolica GM 判据移植 + 经典提取（独立倍式 + input_heads、零约化）。cyclic-5 ℤ₁₃ 2609 s → 31 ms（≈85 000×）且首次通过 `is_groebner_basis`；cyclic-6 可解（9970 s）；<5s 推迟到 0.15.2（LM 索引 + 稀疏 echelon） |
| 0.15.2 | 1.0 候选 | ✅ | ✅ reducer LM 哈希索引（support-mask 桶 + 子掩码枚举）+ 稀疏行 echelon（双指针归并相消，O(nnz)/次）+ 提取查重哈希化 + worklist 预处理 + 行模板缓存。cyclic-6 ℤ₁₃ 9970 s → 3670 s（2.7×，basis=20 正确）；阶段占比转为消除主导（echelon ≈89%）；<5s 未达（cyclic-6 F4 矩阵第 22 轮达 264k 行 × 284k 列，为 F4 固有规模，进一步数量级提升需 F5 签名约简，列入 post-1.0） |
| 0.16.0 | 1.0 候选 | ✅ | ✅ 任意多元因式分解（Wang EEZ）：泛型多元 Diophantine + 逐变量 EEZ Hensel 提升 + $n$ 元 GCD + 特征 $p$ $p$ 次幂 + Wang 首项系数预处理[常数 LC] + Zassenhaus 重组；`factor()` 泛化到任意变量数；顺手修复 3 个既有 bug（`div_rem_sparse` 整除方向、Diophantine 循环上界、单变量非首一分解）；因式分解从 🟡 升级为 🟢（一元/二元/任意多元）；新增 0.16.1（非常数 LC 强加 + 稀疏化）；版本提升 0.16.0。 |
| 0.16.1 | 1.0 候选 | ✅ | ✅ 非常数首项系数强加（模 p Hensel）+ 多元稀疏化改进（ℤ 路径） |
| 0.16.2 | 1.0 候选 | ✅ | ✅ 𝔽_p 非常数 LC 因式分解（Fp Wang LC 重建 + 稀疏 Diophantine 小素数升级，Fp 路径） |
| 0.17.0 | 1.0 候选 | ✅ | ✅ 代数数域因式分解（Trager）：`AlgebraicExtension<D>`（ℚ(α) 与 GF(p^d) 同一实现）+ 平移范数（求值–插值结式）+ 数域模 GCD（GF(p^d) + CRT + 有理重构）；Brown PRS 结式一般次数 bug 按 Symbolica 重移植；deg ≤ 12 ANF 分解 8–32 ms |
| 0.17.1 | 1.0 候选 | ✅ | ✅ 代数数域 Python/C 绑定：`AlgebraicExtension`/`AlgebraicElement`/`AlgebraicPolynomial` Python 类 + `OcasAlgebraicField`/`OcasAlgebraicPoly` 不透明句柄与 `ocas_algebraic_*` C ABI；`RootOf(poly, idx)` 解析确认 |
| 0.18.0 | 1.0 候选 | ✅ | ✅ 数值积分（Vegas 自适应蒙特卡洛 + `integrate_1d` + `StatisticsAccumulator` + `Integrator` trait）、前向自动微分（`HyperDual<T>` 运行时形状 + 截断乘法表 + 几何级数求逆 + `DualCoeff` trait，Rational 双路径）、fuel 资源控制（`Fuel = Arc<AtomicUsize>` + `OutOfFuel` + `simplify_with_fuel`/`integrate_with_fuel`）、张量基础（独立 `Tensor` 类型 + 指标槽 + 显式收缩 + `symmetrise_sign`）；新增 `rand`/`rand_xoshiro` |
| 0.18.1 | 1.0 候选 | ✅ | ✅ 补丁：0.18.0 三项能力（数值积分/双数 AD/张量基础）的 Python/C 绑定补齐（`ocas-py::{numeric,tensor,dual}` + `ocas-c::{numeric,tensor,dual}` 不透明句柄与 C ABI + `include/ocas.h` 同步）+ prelude 补齐张量/双数/`StatisticsAccumulator` 导出；新增 41 Python 测试 + 31 C API 测试；无算法变更，差距结论不变。 |

0.1–0.18.1 交付物全部落地，workspace 版本锁定 0.18.1。质量门全绿：
`cargo fmt`、`clippy -D warnings`、workspace 测试、`cargo deny`、pytest、
`mdbook build`。

---

## 2. 代码规模

`src/` 下 Rust 行数快照（非空行，不含测试与基准）。

| Crate | 文件数 | 行数 |
|---|---|---|
| ocas-poly | 24 | ~15,587 |
| ocas-calc | 18 | ~5,672 |
| ocas-domain | 12 | ~4,475 |
| ocas-eval | 16 | ~4,379 |
| ocas-c | 8 | ~3,195 |
| ocas-py | 11 | ~2,570 |
| ocas-rewrite | 7 | ~1,653 |
| ocas-atom | 5 | ~1,558 |
| ocas-core | 6 | ~1,269 |
| ocas-parse | 3 | ~495 |
| ocas (prelude) | 1 | ~125 |
| ocas-gpl | 1 | 0（占位） |
| **src 合计** | **112** | **~40.9k** |

较 0.15.1 快照（95 文件 / ~30.7k 行）增长约 33%，较 0.10 快照
（66 文件 / ~18k 行）增长约 127%。0.15.1 → 0.18.1 的增量主要来自任意
多元 + 代数数域因式分解（ocas-poly，+~5.0k）、域层（`algebraic` + `dual`，
ocas-domain +~1.1k）、数值积分/流式求值（ocas-eval，+~0.5k），以及三项
0.18.0 能力的 Python/C 绑定扩展（ocas-c +~1.7k，ocas-py +~1.1k）。

`ocas-gpl` 为占位；GPL 专属后端属 Post-1.0 工作，符合路线图。

---

## 3. 算法深度核验

本节是决定 CAS 成熟度最关键的因素，也是差距的主要来源。

| 算法领域 | oCAS 现状 | 成熟度 |
|---|---|---|
| 多项式因式分解 | `DenseUnivariatePolynomial` 上 ℤ 与 ℤ_p 的 `factor()`，`SparseMultivariatePolynomial` 上任意多元 ℤ 与 ℤ_p 的 `factor()`（0.16.x Wang EEZ + 非常数 LC 强加），以及 `AlgebraicNumberField` 上的一元 `factor()`（0.17.0 Trager：平移范数 + 模 GCD） | 🟢 一元/二元/任意多元/代数数域（一元） |
| Gröbner 基 | F4 真实线性代数（0.15.1）+ F5 签名约简（0.19.0：`Signature`/`SyzygySet` + ℤ_p 原生快速路径 `f5_fp`）+ FGLM + 统一 `groebner_basis()` 分派 + ℤ_p 原生 i64 管线；cyclic-6 ℤ₁₃ **2.63 s**（基线 3670 s，≈1400×）；cyclic-5 ℤ₁₃ 0.05 s | 🟢 F4 + F5 完成 |
| 符号积分 | Risch（初等超越塔 + RDE 多项式片段）+ 有理函数 Hermite + 三角 exp(I·x) + 特殊函数表（erf/Ei/Si/Ci/Fresnel）；回退 `Integral(...)` | 🟢 Risch 完成 |
| 实根隔离 | Sturm 序列 + 区间隔离 + refine（单变量）；已知缺口：Wilkinson n=10 展开多项式仅隔离 8/10 根 | 🟡 较完整 |
| 多项式 GCD | GCD + 本原部分 + 扩展 GCD（0.12）+ 经 EEZ 的任意元数多元 GCD（0.16）+ GF(p^d) 上模数域 GCD（CRT + 有理重构，0.17）；大整数系数尚无模 GCD 快速路径 | 🟢 可用，无 HEVMGCD |
| 线性求解 | 有理/整数线性方程组 + 二元丢番图（`ax+by=c`） | 🟡 可用，规模有限 |
| JIT 求值 | Cranelift 后端；≥10x 加速目标达成（按路线图标准） | 🟢 完整 |

---

## 4. 与参考系统的差距

### 4.1 对照 Symbolica（Rust，AGPL）

Symbolica 的 `examples/` 目录（30 个示例）是成熟度基准。0.11–0.18 之后
oCAS 已覆盖 Symbolica 核心示例面的**全部**；差距收敛到**大规模性能**
（cyclic-6 Gröbner）与少数**专用模式变换器**（如参数序列分拆用的
`Transformer::Partition`）。

| 能力 | oCAS | Symbolica |
|---|---|---|
| 多项式因式分解 | ✅ 一元 ℤ/ℤ_p（CZ + Hensel + Zassenhaus）+ 任意多元（0.16 Wang EEZ + 非常数 LC 强加 0.16.1/0.16.2）+ 代数数域（0.17 Trager，一元） | ✅ 完整（任意多元 + 代数数域，`factorization.rs`） |
| 有理多项式 | ✅ 含 GCD 规范化的 `RationalPolynomial<D,O>` | ✅ `rational_polynomial.rs` |
| 部分分式 | ✅ 任意 `EuclideanDomain` 上的 `apart()` / `together()` | ✅ `partial_fraction.rs` |
| 有理重构 | ✅ 基于扩展欧几里得的 `rational_reconstruction(a, m)` | ✅ `rational_reconstruction.rs` |
| 数值积分 | ✅ Vegas 自适应蒙特卡洛 + `integrate_1d` + `StatisticsAccumulator`（0.18） | ✅ `numerical_integration.rs` |
| 流式 API | ✅ `StreamingEvaluator`：分块输入 + 复用栈，百万行恒定内存 | ✅ `streaming.rs` |
| 张量 / 双数 | ✅ 独立 `Tensor` 类型 + 指标收缩 + `symmetrise_sign`（0.18 基础版，完整规范化 Post-1.0）；`HyperDual<T>` 前向 AD（0.18） | ✅ `tensors.rs` / `dual.rs`（基于 graphica 的完整规范化） |
| 优化 / 代码生成 | ✅ 多输出 JIT（`compile_multi` + CSE + 常量折叠 + 栈压缩）+ f32 混合精度 | ✅ `optimize.rs` / 多输出 |
| Gröbner 基 | � F4（0.15.1 真实线性代数）+ F5 签名约简（0.19.0）；cyclic-6 ℤ₁₃ **2.63 s**（较 0.15.2 ≈1400×）；F4/F5 统一分派；通用域 + ℤ_p 快速路径 | ✅ 工业级 |
| 资源控制（fuel） | ✅ `Fuel = Arc<AtomicUsize>` + `simplify_with_fuel`/`integrate_with_fuel`（0.18） | ✅ `fuel_backend.rs` |
| 模式变换器 | 🟡 matcher/replace/transformer 完备；专用序列 `Transformer::Partition` 未实现 | ✅ 完整变换器集 |

Symbolica 2.1.0 的核心竞争力——工业级因式分解（含代数数域）、有理函数
运算、多输出优化、流式 API、数值积分、双数、张量、fuel 资源控制——
oCAS 已在 0.11–0.18 期间**全部补齐**。剩余差距为：**Gröbner 大规模性能**
（cyclic-6 量级，Symbolica 的 F5/签名机制仍领先）、**张量完整规范化**
（oCAS 仅基础版，Symbolica 用 graphica 图同构引擎）、以及少数**专用模式
变换器**（如 `Transformer::Partition`）。

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
| 积分 | 🟢 基本持平 | 双方均有 Risch（oCAS 0.14 补齐）；SymPy 的启发式/manual 回退覆盖更广，oCAS 未覆盖情形返回 `Integral(...)` |
| 因式分解 | 🟢 持平 | 一元 ℤ/ℤ_p + 任意多元（0.16 Wang EEZ）+ 代数数域（0.17 Trager）；SymPy 的 ANF 覆盖更广 |
| Gröbner | 🟢 oCAS 优势 | oCAS F4 矩阵化 + 真实线性代数（cyclic-5 ℤ₁₃ 23 ms），优于 SymPy 的 Buchberger 实现 |
| 矩阵/线性代数 | 🟢 持平 | oCAS 有 Bareiss 行列式/逆 |
| **性能** | 🟢 **oCAS 优势** | Rust + Cranelift JIT + arena 对纯 Python；实测 x³⁰−1 无平方分解 39 µs vs SymPy 完全分解 ~0.9 ms（~24×，2026-07-21） |
| Python 易用性 | 🟢 持平 | oCAS 有 `ocas-py` 绑定 |

0.6.0 成功标准"基础多项式/微积分/重写与 SymPy 持平"——在**性能**维度已
达成并领先，**积分**经 0.14 Risch 补齐，**因式分解**经 0.16 达任意多元
持平（0.17 另补 Trager 代数数域）。与 SymPy 的剩余功能差距集中在**积分
启发式回退的广度**（SymPy 的 `manualintegrate`/启发式池比 oCAS 的
Risch + 查表路径更宽）。

---

## 5. 关键缺口与优先级

按"影响面 × 实现成本"排序。1.0 前规划的硬算法缺口已**全部闭合**；阶段
B+ "Symbolica 差距清零"（0.15.2–0.18.0）已完成——详见 EVOLUTION_PLAN。
剩余项为大规模性能、广度与 Post-1.0 议题。

| # | 缺口 | 优先级 |
|---|---|---|
| 1 | ~~完整多项式因式分解~~（0.11.0–0.11.1 完成） | ✅ 已完成——一元与二元（关于 x 首一）闭合；≥3 变量见 #7 |
| 2 | ~~Risch 符号积分~~（0.14 完成） | ✅ 已完成——初等超越塔 + RDE 片段 + 有理函数 Hermite + 特殊函数表 |
| 3 | ~~Gröbner F4/F5~~（0.13 / 0.14 / 0.15.1 完成） | ✅ F4 真实线性代数 + FGLM + 实验性 F5；大规模性能见 #6 |
| 4 | ~~有理多项式/部分分式~~（0.12 完成） | ✅ 已完成——`RationalPolynomial` 类型 + 部分分式 + 结式 + Karatsuba 乘法 |
| 5 | ~~多输出优化/代码生成~~（0.15 完成） | ✅ 已完成——多输出 JIT（97×/21×）+ f32 混合精度 + CSE/常量折叠/栈压缩 |
| 6 | ~~Gröbner 大规模性能（cyclic-6 ℤ_p < 5 s）~~（0.19 完成） | ✅ 完成 — F5 签名约简（0.19.0）：cyclic-6 ℤ₁₃ 3670 s → **2.63 s**（≈1400×）；F4/F5 统一分派；通用域 + ℤ_p 原生快速路径均验证 |
| 7 | ~~任意多元（≥3 变量）因式分解~~（0.16 完成） | ✅ 已完成——Wang EEZ 提升 + 首项系数预处理（常数 LC）+ Zassenhaus 重组；非常数 LC 强加见 #7a |
| 7a | ~~非常数首项系数强加 + 多元稀疏化~~（0.16.1/0.16.2 完成） | ✅ 已完成——模 p Hensel 强加 + 稀疏 Diophantine + Fp 路径域版 Wang 预处理 |
| 8 | ~~代数数域因式分解~~（0.17 完成） | ✅ 已完成——Trager 算法（平移范数 + ℚ 分解 + GF(p^d) 模 GCD），一元路径；多元扩域留待后续 |
| 9 | ~~数值积分 / 双数 / 张量基础 / fuel 资源控制~~（0.18 完成） | ✅ 已完成——Vegas + HyperDual + 指标收缩 + fuel；0.18.1 补齐 Python/C 绑定 |
| 10 | ODE 求解器（阶段 B++ 0.20） | 🟢 SageMath/SymPy 对齐；一阶/二阶 + 系统 + 级数 + Laplace |
| 11 | 数论栈（阶段 B++ 0.21） | 🟢 SageMath/PARI 对齐；模 GCD + 整数分解 + 素性 + 离散对数 + CRT |
| 12 | 张量完整规范化 + 专用模式变换器（阶段 B++ 0.22） | 🟡 Symbolica 最后阵地；需图同构引擎 |
| 13 | 代数几何工具（阶段 B++ 0.23） | 🟢 SageMath/Singular 对齐；理想运算 + RUR + 准素分解 + Hilbert 级数 |
| 14 | PDE 求解器（Post-1.0） | 🟢 用户期望高；Poisson/热传导/波动 |

---

## 6. 总评

0.1 → 0.18.1 执行质量很高：每个路线图交付物均兑现，分层架构干净（无环
依赖），13 crate workspace 严格分层，质量门严格（`-D warnings` + deny +
Miri 意识），文档/绑定/CI 工程化完备。1.0 前规划的三大硬算法——多项式
因式分解（0.11）、Gröbner F4（0.13，真实线性代数于 0.15.1 修复）、
Risch 符号积分（0.14）——全部闭合，并经 SymPy/Symbolica 交叉验证框架
持续回归。

务实定位：当前 oCAS 是"高性能、自包含的代数内核，功能对标 SymPy，并
近全覆盖 Symbolica 的示例面"。具体交付：Risch 符号积分、一元/二元/任意
多元因式分解（另含 Trager 代数数域）、有理函数、真实线性代数的 Gröbner
F4、多输出 JIT/流式求值、Vegas 数值积分、hyper-dual 前向 AD、张量基础、
fuel 资源控制。0.15.1 复测性能（仍具代表性）：F4 cyclic-5 ℤ₁₃ 23 ms；
x³⁰−1 无平方分解 39 µs（SymPy 完全分解 ~0.9 ms，~24×）；JIT 单输出 97×、
三输出 21×。

阶段 B+ "Symbolica 差距清零"（0.15.2 → 0.18.0）**已完成**：0.15.1 时仍
开着的每个 Symbolica 示例域缺口——任意多元因式分解、代数数域因式分解、
数值积分、双数、张量、fuel——现已全部闭合。随后阶段 B++ "竞品全面对齐"
（0.19.0 → 0.23.0，详见 EVOLUTION_PLAN）在 1.0.0 冻结前瞄准剩余缺口：
cyclic-6 量级 Gröbner 性能（F5 签名约简，0.19）、ODE 求解器
（SageMath/SymPy 对齐，0.20）、数论（SageMath/PARI 对齐，0.21）、张量
完整规范化 + 高级模式匹配（Symbolica 最后阵地，0.22）、代数几何工具
（SageMath/Singular 对齐，0.23）。阶段 B++ 之后，1.0.0 严格**仅做稳定性
与发布工程**（API 冻结、覆盖率、迁移指南、签名产物）。

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
| 0.13.2 | 2026-07-18 | 工程与发布里程碑：`pip install ocas` 上线 PyPI（5 平台 wheel + sdist，含 macOS 双架构）；打通 OIDC trusted publishing；修复 crossbeam-epoch RUSTSEC-2026-0204；cranelift/chumsky/logos/cbindgen/criterion/hashbrown/flint3-sys/egg 依赖升级；无算法变更，差距结论不变。 |
| 0.14.0 | 2026-07-18 | Risch 符号积分完成（初等超越塔 + RDE 多项式片段）；有理函数积分（Hermite + 对数部分）；特殊函数表（erf/Ei/Si/Ci/Fresnel）闭合 0.11.0 已知差距 `exp(-x²)→erf`；三角 exp(I·x) + realify；Gröbner 收尾（FGLM 零维换序 + F5 实验性 + Hilbert 界 + reorder）；解析器 `-x^2` 优先级修复；符号积分从 🟡 升级为 🟢；最高优先级缺口转为 0.15 性能/多输出 JIT。 |
| 0.15.0 | 2026-07-20 | 多输出 JIT（97×/21×）+ f32 混合精度 + 流式求值（百万行恒定内存）+ 常量折叠/栈压缩 + Arena reset/workspace 池 + ahash + 原生 i64 F4 管线；JIT 调用约定 Windows 修复；分段插装定位 F4 瓶颈（extract 99.98%）；cyclic-6 <5s 推迟到 0.15.1（需 RREF/F5）；最高优先级缺口转为 1.0 稳定版。 |
| 0.15.1 | 2026-07-20 | F4 真实线性代数修复：矩阵列序降序（此前升序致 echelon 形同虚设，F4 实为 Buchberger）+ echelon 回写条件 + Symbolica GM 判据移植 + 经典提取（独立倍式 + input_heads、提取零约化）；cyclic-5 ℤ₁₃ 2609 s → 31 ms（≈85 000×）且首次通过 `is_groebner_basis`；cyclic-6 可解（9970 s，basis=20）；<5s 推迟到 0.15.2（LM 索引 + 稀疏 echelon）。 |
| 0.16.0–0.16.2 | 2026-07-21 | 任意多元因式分解栈（Wang EEZ + Hensel + 非常数首项系数强加 + 稀疏 Diophantine 小素数升级），覆盖 ℤ 与 𝔽ₚ 多元路径；多元因式分解从 🔴 升级为 🟢。 |
| 0.17.0 | 2026-07-22 | 代数数域因式分解（Trager）完成：`AlgebraicNumberField` + 数域模 GCD（GF(p^d) + CRT + 有理重构）+ 平移范数；修复结式 Brown PRS 一般次数 bug；代数数域因式分解从 🔴 升级为 🟢（单变量路径）。 |
| 0.17.1 | 2026-07-22 | 补丁：代数数域 Python/C 绑定收尾（`AlgebraicExtension`/`AlgebraicElement`/`AlgebraicPolynomial` Python 类 + `OcasAlgebraicField`/`OcasAlgebraicPoly` 不透明句柄与 `ocas_algebraic_*` C ABI + `RootOf` 解析确认）；无算法变更，差距结论不变。 |
| 0.18.0 | 2026-07-23 | 数值积分（Vegas 自适应蒙特卡洛 + `integrate_1d`）、前向自动微分（`HyperDual<T>` 运行时形状）、fuel 资源控制（`Fuel` + `simplify_with_fuel`/`integrate_with_fuel`）、张量基础（独立 `Tensor` 类型 + 显式收缩 + 对称化符号）落地；新增 `rand`/`rand_xoshiro` 依赖；张量完整规范化与确定性 quadrature 桥接推迟。 |
| 0.18.1 | 2026-07-23 | 补丁：0.18.0 三项能力（数值积分/双数 AD/张量基础）的 Python/C 绑定补齐——`ocas-py::{numeric,tensor,dual}` 模块 + `ocas-c::{numeric,tensor,dual}` 不透明句柄与 C ABI + `include/ocas.h` 同步 + prelude 补齐张量/双数/`StatisticsAccumulator` 导出；41 Python 测试 + 31 C API 测试；无算法变更，差距结论不变。 |
| 0.15.1 | 2026-07-21 | 重新评估：代码规模快照更新至 95 文件 / ~30.7k 行（较 0.10 的 ~18k 增长 ~70%）；F4 cyclic-5 ℤ₁₃ 复测 23 ms；新增实测 x³⁰−1 无平方分解 39 µs vs SymPy 完全分解 ~0.9 ms（~24×）；修正 0.14/0.15 后的过时表述（§3 GCD/实根隔离、§4.1 "基本缺失"段落、§4.3 积分/因式分解/Gröbner、§5 Risch 优先级、乱码字符）；缺口重排——1.0 前硬算法全部闭合，剩余项转为 Post-1.0：任意多元（≥3 变量）与代数数域因式分解、数值积分、张量/双数、ODE/PDE，cyclic-6 <5s 定界 0.15.2。 |
| 0.15.2 | 2026-07-21 | Gröbner 大规模性能：reducer LM 哈希索引（support-mask 桶 + 子掩码枚举，消除 O(单项式×基) 线性扫描）+ 稀疏行 echelon（双指针归并相消 O(nnz)/次，替代稠密 buffer）+ 提取查重哈希化 + worklist 预处理 + 行模板缓存；cyclic-6 ℤ₁₃ 9970 s → 3670 s（2.7×，basis=20 正确），阶段占比转为消除主导（echelon ≈89%）；<5s 未达——cyclic-6 F4 矩阵第 22 轮达 264k 行 × 284k 列，为 F4 固有规模，进一步数量级提升需 F5 签名约简（消除零约化行），列入 post-1.0；版本提升 0.15.2。 |
| 0.16.0 | 2026-07-21 | 任意多元因式分解（Wang EEZ）完成：落地 `factor::eez`（泛型多元 Diophantine + 逐变量 EEZ Hensel 提升 + $n$ 元 GCD + 特征 $p$ $p$ 次幂 + Wang 首项系数预处理[常数 LC] + Zassenhaus 重组）；`factor()` 泛化到任意变量数；顺手修复 3 个既有 bug（`div_rem_sparse` 整除方向、Diophantine 循环上界、单变量非首一分解）；因式分解从 🟡 升级为 🟢（一元/二元/任意多元）；新增 0.16.1（非常数 LC 强加 + 稀疏化）；版本提升 0.16.0。 |
| 0.17.0 | 2026-07-22 | 代数数域因式分解（Trager）完成：新增 `ocas-domain::algebraic`（`AlgebraicExtension<D>`：ℚ(α) 与 GF(p^d) 同一实现，EEA 求逆）+ `ocas-poly::factor::algebraic`（平移范数[求值–插值结式] + 数域模 GCD[GF(p^d) + CRT + 有理重构 + 试除] + 有理快速通道）；修复结式 Brown PRS 一般次数 bug（β 除法仅在单位时执行的非法实现，按 Symbolica `resultant_prs` 重移植）；0.16.2 稀疏 Diophantine 小素数升级启发式补齐；因式分解能力覆盖 一元/二元/任意多元/代数数域（一元）；性能指标达成（deg≤12 实测 8–32 ms < 100 ms）；版本提升 0.17.0。 |
| 0.18.1 | 2026-07-23 | **全面重新评估**（0.16–0.18 落地后）。代码规模快照刷新至 112 文件 / ~40.9k 行（较 0.15.1 的 95 文件 / ~30.7k 增长 33%；较 0.10 的 ~18k 增长 127%）。§1 版本表扩展至 0.18.1（新增 0.16.0–0.18.1 行）。§3 多项式 GCD 从 🟡 升级为 🟢（经 EEZ 的任意元数多元 GCD [0.16] + 模数域 GCD [0.17]）。§4.1 Symbolica 差距表重写：数值积分/张量/双数/fuel 全部从 🔴 升级为 ✅（0.18 闭合）；因式分解行补注 ANF 已完成（0.17）；新增模式变换器行（🟡，缺 `Transformer::Partition`）；收尾段落重写——除大规模 Gröbner + 张量完整规范化外，Symbolica 示例域缺口全部闭合。§4.3 SymPy 因式分解从 🟡 升级为 🟢（任意多元持平，0.16）。§5 新增 #11（张量规范化 + 专用模式变换器，Post-1.0）；表头重写——阶段 B+ 宣告完成。§6 总评重写——1.0 仅剩稳定性/发布工程。全文多处乱码字符修复。 |
| 0.19.0 | 2026-07-23 | **F5 Gröbner 基发布——cyclic-6 规模缺口闭合。** §3 Gröbner 行从 🟡 升级为 🟢（F5 签名约简）。§4.1 Gröbner 基竞品行从 🟡 升级为 🟢。§5 #6（Gröbner 大规模性能）标记 ✅ 完成——cyclic-6 ℤ₁₃ 3670 s → **2.63 s**（≈1400×），经 `f5_fp` ℤ_p 原生快速路径；cyclic-5 0.05 s；通用域 + ℤ_p 路径均验证。统一 `groebner_basis()` 分派（`Algorithm::{Auto,F4,F5,Buchberger}`）。多序（`WeightOrder`/`BlockOrder`）推迟到 0.19.1（trait 重构）。 |
| 0.19.1 | 2026-07-23 | **MonomialOrder trait 重构 + WeightOrder/BlockOrder 发布。** `Copy` + 静态分派 → `Clone + Default` + 方法分派（`&self`）；`PhantomData<O>` → `order: O` 字段；新增 `WeightOrder`（加权序）与 `BlockOrder`（分块序）+ `SubOrder` 枚举；11 处 `O::cmp` 调用点全部更新；`Signature::cmp_pot` 签名新增 `order: &O` 参数。多序支持标记从 `[~]` 升级为 `[x]`。 |