# 差距分析：oCAS 与参考系统

本文档逐里程碑（0.1 → 1.0+）跟踪 oCAS 的实现完成度，并对照三大参考系统
评估差距：**Symbolica**（Rust）、**SageMath**（Python 生态）、**SymPy**
（纯 Python）。本文档为活文档，每次版本发布时必须更新。英文版见
[GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md)。

> 最后评估：**0.18.1 @ 2026-07-23**（0.18.0 三项能力的 Python/C 绑定补齐 + prelude 补齐）

---

## 图例

| 标记 | 含义 |
|---|---|
| ✅ | 完成 |
| 🟡 | 基础可用或部分完成 |
| 🔴 | 缺失或重大缺口 |
| ⚠️ | 完成但有保留 |

---

## 1. 版本完成状态（0.1–0.15.2）

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
| 0.15.2 | 1.0 候选 | ✅ | ✅ reducer LM 哈希索引（support-mask 桶 + 子掩码枚举）+ 稀疏行 echelon（双指针归并相消，O(nnz)/次）+ 提取查重哈希化 + worklist 预处理 + 行模板缓存。cyclic-6 ℤ₁₃ 9970 s → 3670 s（2.7×，basis=20 正确）；阶段占比转为消除主导（echelon ≈89%）；<5s 未达（264k 行是 F4 固有规模，需 F5 签名约简） |

0.1–0.15.2 交付物全部落地，workspace 版本锁定 0.15.2。质量门全绿：
`cargo fmt`、`clippy -D warnings`、workspace 测试、`cargo deny`、pytest、
`mdbook build`。

---

## 2. 代码规模

`src/` 下 Rust 行数快照（不含测试与基准）。

| Crate | 文件数 | 行数 |
|---|---|---|
| ocas-poly | 22 | ~10,560 |
| ocas-calc | 18 | ~5,649 |
| ocas-eval | 13 | ~3,855 |
| ocas-domain | 10 | ~3,337 |
| ocas-rewrite | 7 | ~1,593 |
| ocas-py | 7 | ~1,461 |
| ocas-c | 4 | ~1,454 |
| ocas-core | 5 | ~1,115 |
| ocas-atom | 4 | ~1,111 |
| ocas-parse | 3 | ~495 |
| ocas (prelude) | 1 | ~115 |
| ocas-gpl | 1 | 0（占位） |
| **src 合计** | **95** | **~30.7k** |

较 0.10 快照（66 文件 / ~18k 行）增长约 70%，增量主要来自 Risch 积分与
有理函数积分（ocas-calc）、F4/FGLM/F5 与因式分解（ocas-poly）、多输出
JIT / 流式求值（ocas-eval）。

`ocas-gpl` 为占位；GPL 专属后端属 Post-1.0 工作，符合路线图。

---

## 3. 算法深度核验

本节是决定 CAS 成熟度最关键的因素，也是差距的主要来源。

| 算法领域 | oCAS 现状 | 成熟度 |
|---|---|---|
| 多项式因式分解 | `DenseUnivariatePolynomial` 上 ℤ 与 ℤ_p 的 `factor()`，`SparseMultivariatePolynomial` 上任意多元 ℤ 与 ℤ_p 的 `factor()`（0.16.x Wang EEZ + 非常数 LC 强加），以及 `AlgebraicNumberField` 上的一元 `factor()`（0.17.0 Trager：平移范数 + 模 GCD） | 🟢 一元/二元/任意多元/代数数域（一元） |
| Gröbner 基 | F4 真实线性代数（0.15.1：降序列序 + Symbolica GM 判据 + 经典提取）+ FGLM + 实验性 F5 + ℤ_p 原生 i64 管线；cyclic-5 ℤ₁₃ 23 ms（2026-07-21 复测） | 🟢 F4 完成 |
| 符号积分 | Risch（初等超越塔 + RDE 多项式片段）+ 有理函数 Hermite + 三角 exp(I·x) + 特殊函数表（erf/Ei/Si/Ci/Fresnel）；回退 `Integral(...)` | 🟢 Risch 完成 |
| 实根隔离 | Sturm 序列 + 区间隔离 + refine（单变量）；已知缺口：Wilkinson n=10 展开多项式仅隔离 8/10 根 | 🟡 较完整 |
| 多项式 GCD | GCD + 本原部分 + 扩展 GCD（0.12）；无模方法 GCD（大整数系数）优化 | 🟡 可用 |
| 线性求解 | 有理/整数线性方程组 + 二元丢番图（`ax+by=c`） | 🟡 可用，规模有限 |
| JIT 求值 | Cranelift 后端；≥10x 加速目标达成（按路线图标准） | 🟢 完整 |

---

## 4. 与参考系统的差距

### 4.1 对照 Symbolica（Rust，AGPL）

Symbolica 的 `examples/` 目录揭示了成熟度差距。0.11–0.15 之后 oCAS 已
覆盖 Symbolica 核心功能面的大部，差距收敛到广度与大规模性能维度。

| 能力 | oCAS | Symbolica |
|---|---|---|
| 多项式因式分解 | � 一元 ℤ 与 ℤ_p（CZ + Hensel + Zassenhaus）+ 任意多元（0.16.0 Wang EEZ，含常数 LC 预处理）；代数数域缺失；非常数 LC 强加 0.16.1 | ✅ 完整（任意多元 + 代数数域，`factorization.rs`） |
| 有理多项式 | ✅ 含 GCD 规范化的 `RationalPolynomial<D,O>` | ✅ `rational_polynomial.rs` |
| 部分分式 | ✅ 任意 `EuclideanDomain` 上的 `apart()` / `together()` | ✅ `partial_fraction.rs` |
| 有理重构 | ✅ 基于扩展欧几里得的 `rational_reconstruction(a, m)` | ✅ `rational_reconstruction.rs` |
| 数值积分 | 🔴 无 | ✅ `numerical_integration.rs` |
| 流式 API | ✅ `streaming.rs`（`StreamingEvaluator`：分块输入 + 复用栈，百万行恒定内存） | ✅ `streaming.rs` |
| 张量 / 双数 | 🔴 无 | ✅ `tensors.rs` / `dual.rs` |
| 优化 / 代码生成 | ✅ 多输出 JIT（`compile_multi` + CSE + 常量折叠 + 栈压缩）+ f32 混合精度 | ✅ `optimize.rs` / 多输出 |
| Gröbner 基 | 🟡 F4 完成 + 大规模性能优化（0.15.2：LM 索引 + 稀疏 echelon + 行模板缓存，cyclic-6 ℤ₁₃ 9970 s → 3670 s）；cyclic-6 <5s 未达，需 F5 签名约简 | ✅ 工业级 |
| 资源控制（fuel） | 🔴 无 | ✅ `fuel_backend.rs` |

Symbolica 2.1.0 的核心竞争力——工业级因式分解、有理函数运算、多输出
优化、流式 API——oCAS 已在 0.11–0.15 期间大部补齐。剩余差距集中在：
任意多元（≥3 变量）与代数数域因式分解、数值积分、张量/双数、fuel
资源控制，以及 Gröbner 的大规模性能（cyclic-6 量级，Symbolica 仍显著
领先）。

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
| 积分 | � 基本持平 | 双方均有 Risch（oCAS 0.14 补齐）；SymPy 的启发式/manual 回退覆盖更广，oCAS 未覆盖情形返回 `Integral(...)` |
| 因式分解 | 🟡 oCAS 略弱 | 一元 ℤ 与 ℤ_p 持平（CZ + Hensel + Zassenhaus）；oCAS 支持二元，SymPy 支持任意多元 |
| Gröbner | 🟢 oCAS 优势 | oCAS F4 矩阵化 + 真实线性代数（cyclic-5 ℤ₁₃ 23 ms），优于 SymPy 的 Buchberger 实现 |
| 矩阵/线性代数 | 🟢 持平 | oCAS 有 Bareiss 行列式/逆 |
| **性能** | 🟢 **oCAS 优势** | Rust + Cranelift JIT + arena 对纯 Python；实测 x³⁰−1 无平方分解 39 µs vs SymPy 完全分解 ~0.9 ms（~24×，2026-07-21） |
| Python 易用性 | 🟢 持平 | oCAS 有 `ocas-py` 绑定 |

0.6.0 成功标准"基础多项式/微积分/重写与 SymPy 持平"——在**性能**维度已
达成并领先，**积分**经 0.14 Risch 补齐；与 SymPy 的剩余功能差距集中在
**任意多元因式分解**与**积分启发式回退的广度**。

---

## 5. 关键缺口与优先级

按"影响面 × 实现成本"排序。1.0 前规划的硬算法缺口已全部闭合；与
Symbolica 的剩余差距已排入阶段 B+（0.15.2–0.18.0，详见
EVOLUTION_PLAN），目标 1.0 前清零。

| # | 缺口 | 优先级 |
|---|---|---|
| 1 | ~~完整多项式因式分解~~（0.11.0–0.11.1 完成） | ✅ 已完成——一元与二元（关于 x 首一）闭合；≥3 变量见 #7 |
| 2 | ~~Risch 符号积分~~（0.14 完成） | ✅ 已完成——初等超越塔 + RDE 片段 + 有理函数 Hermite + 特殊函数表 |
| 3 | ~~Gröbner F4/F5~~（0.13 / 0.14 / 0.15.1 完成） | ✅ F4 真实线性代数 + FGLM + 实验性 F5；大规模性能见 #6 |
| 4 | ~~有理多项式/部分分式~~（0.12 完成） | ✅ 已完成——`RationalPolynomial` 类型 + 部分分式 + 结式 + Karatsuba 乘法 |
| 5 | ~~多输出优化/代码生成~~（0.15 完成） | ✅ 已完成——多输出 JIT（97×/21×）+ f32 混合精度 + CSE/常量折叠/栈压缩 |
| 6 | Gröbner 大规模性能（cyclic-6 ℤ_p < 5 s） | 🔴 0.15.2——LM 哈希索引 + 稀疏 echelon |
| 7 | ~~任意多元（≥3 变量）因式分解~~（0.16 完成） | ✅ 已完成——Wang EEZ 提升 + 首项系数预处理（常数 LC）+ Zassenhaus 重组；非常数 LC 强加见 #7a |
| 7a | ~~非常数首项系数强加 + 多元稀疏化~~（0.16.1/0.16.2 完成） | ✅ 已完成——模 p Hensel 强加 + 稀疏 Diophantine + Fp 路径域版 Wang 预处理 |
| 8 | ~~代数数域因式分解~~（0.17 完成） | ✅ 已完成——Trager 算法（平移范数 + ℚ 分解 + GF(p^d) 模 GCD），一元路径；多元扩域留待后续 |
| 9 | ~~数值积分 / 双数 / 张量基础 / fuel 资源控制~~（0.18 完成） | ✅ 已完成——Vegas + HyperDual + 指标收缩 + fuel；0.18.1 补齐 Python/C 绑定 |
| 10 | ODE/PDE 求解器（Post-1.0） | 🟢 用户期望高 |

---

## 6. 总评

0.1 → 0.15.1 执行质量很高：每个路线图交付物均兑现，分层架构干净（无环
依赖），12 crate workspace 严格分层，质量门严格（`-D warnings` + deny +
Miri 意识），文档/绑定/CI 工程化完备。1.0 前规划的三大硬算法——多项式
因式分解（0.11）、Gröbner F4（0.13，真实线性代数于 0.15.1 修复）、
Risch 符号积分（0.14）——全部闭合，并经 SymPy/Symbolica 交叉验证框架
持续回归。

务实定位：当前 oCAS 是"高性能 SymPy 核心 + Risch 符号积分 + 一元/二元
因式分解与有理函数 + 真实线性代数的 Gröbner F4 + 多输出 JIT / 流式
求值"。0.15.1 复测性能：F4 cyclic-5 ℤ₁₃ 23 ms；x³⁰−1 无平方分解
39 µs（SymPy 完全分解 ~0.9 ms，~24×）；JIT 单输出 97×、三输出 21×。

1.0 前剩余工作：阶段 B+ "Symbolica 差距清零"（0.15.2 Gröbner 大规模性能
→ 0.16 任意多元因式分解 ✅ → 0.16.1 非常数首项系数强加 ✅ → 0.17 代数数域
因式分解 ✅ → 0.18 数值积分/双数/张量/fuel，详见 EVOLUTION_PLAN），随后
1.0.0 仅做稳定性与发布工程（API 冻结、覆盖率、迁移指南、签名产物）。
ODE/PDE 与完整张量微积分仍为 Post-1.0 议题。

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
| 0.13.1 | 2026-07-17 | 补丁发布：docs.rs 构建改为纯 Rust 特性（不含 gmp/mpfr/flint/python/gpl），托管文档恢复构建；功能与算法层面与 0.13.0 一致，差距结论不变。 || 0.13.2 | 2026-07-18 | 工程与发布里程碑：`pip install ocas` 上线 PyPI（5 平台 wheel + sdist，含 macOS 双架构）；打通 OIDC trusted publishing；修复 crossbeam-epoch RUSTSEC-2026-0204；cranelift/chumsky/logos/cbindgen/criterion/hashbrown/flint3-sys/egg 依赖升级；无算法变更，差距结论不变。 |
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