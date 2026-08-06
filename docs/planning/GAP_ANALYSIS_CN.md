# 差距分析：oCAS 与参考系统

本文档逐里程碑（0.1 → 1.0+）跟踪 oCAS 的实现完成度，并对照参考系统
评估差距：**Symbolica**（Rust，source-available 商业）、**SageMath**
（Python 生态）、**SymPy**（纯 Python）、**msolve**（Gröbner 性能标杆）、
**FLINT**（多项式/数论性能库）。本文档为活文档，每次版本发布时必须更新。
英文版见 [GAP_ANALYSIS_EN.md](GAP_ANALYSIS_EN.md)。
配套文档：[COMPETITIVE_MATRIX_CN.md](COMPETITIVE_MATRIX_CN.md)（竞品能力矩阵）、
[BENCHMARK_SUITE_CN.md](BENCHMARK_SUITE_CN.md)（基准测试套件）。

> 最后评估：**0.26.0 @ 2026-08-06**（竞品版本重新核实：Symbolica 2.2.0（无更新）、
> SymPy 1.14.0（无更新）、SageMath 10.9（2026-05-05）、FLINT 3.6.0（2026-06-29，
> Kinoshita-Li 级数复合/padic_radix/subresultant 结式）、msolve 0.10.1（2026-07-08，
> Gebauer-Möller 改进 + QQ 提升修复）、GiNaC 1.8.10（无更新）、mathcore 0.3.1
> （更正上一版误记的 0.5.0）、Numerica（无 tag 发布，开发活跃）；本机全量复测
> oCAS/Symbolica/SymPy + WSL2 实测 msolve 0.10.1：cyclic-6 grevlex 55 ms 达成
> <0.5 s 里程碑，msolve 实测 4 ms；DoubleFloat 已由 0.24.0 DoubleF64 兑现；
> §5 优先级重排，DoubleFloat 与 cyclic-6 移入已完成项）

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
| 0.19.0/0.19.1 | 1.0 候选 | ✅ | ✅ F5 Gröbner 签名约简（cyclic-6 ℤ₁₃ 3670 s → 2.63 s，≈1400×）；`MonomialOrder` trait 重构 + `WeightOrder`/`BlockOrder` |
| 0.20.0/0.20.1 | 1.0 候选 | ✅ | ✅ ODE 求解器全量：一阶 5 种 + 积分因子；二阶常系数/Cauchy-Euler + VOP + 降阶法 + 待定系数扩展；级数递推 + Frobenius；Laplace IVP（`dsolve_ivp`）；2×2 系统（`dsolve_system`）；Python/C 绑定；31 项代入验证正确性测试 |
| 0.21.0 | 1.0 候选 | ✅ | ✅ 数论与计算代数栈：CRT 多模累加、BPSW 素性 + 2⁶⁴ 确定性 MR、整数分解（试除/Brent rho/Pollard p−1/Williams p+1/ECM Suyama-Montgomery）、BSGS + Pohlig-Hellman 离散对数、φ/μ/τ/σ_k/λ 数论函数；单变量 Brown 模 GCD（`gcd::modular::gcd_modular_z`）+ 二元 `gcd_modular` 完整 Brown 重写（内容分离 + monic 插值像 + 多素数 CRT + 有理重构）；Python/C 绑定（`ocas::ntheory` / `ocas_ntheory_*`）；修复 `rational_reconstruction` 整数平方根性能炸弹；ECM 30 位半素数 1.1 s（<10 s） |
| 0.22.0 | 1.0 候选 | ✅ | ✅ 张量完整规范化 + 高级模式匹配：图同构正则标号、张量规范形、Young 投影子、Partition 变换器、多模式替换、回溯匹配器优化 |
| 0.23.0 | 1.0 候选 | ✅ | ✅ 高级 Gröbner 与代数几何工具：`ocas-poly::ideal`（ideal_contains/sum/product/quotient/saturate/intersection）、MatrixOrder 消元序 + eliminate()、零维求解（Sturm 根隔离）、准素分解（Lex GB 因式分解 + 饱和分离）、根式（无平方 + Jacobian 饱和）、Hilbert 级数/维数/次数/多项式、有理根定理、Python MultivariatePolynomial + C FFI 绑定；212 测试通过 |

0.1–0.23.0 交付物全部落地，workspace 版本锁定 0.21.0。质量门全绿：
`cargo fmt`、`clippy -D warnings`、workspace 测试、`cargo deny`、pytest、
`mdbook build`。

| 版本 | 阶段 | 路线图 | 核验状态 |
|---|---|---|---|
| 0.24.0 | 1.0 候选 | ✅ | ✅ 启发式积分四技术（分部 LIATE/三角换元/Weierstrass/Euler 占位）接入 `try_risch_or_fallback` 调度链 + `integrate_heuristic` 公共 API；**DoubleF64**（Dekker/Knuth 双精度 ~31 位十进制，含超越函数 + EvaluationDomain）；Python/C 绑定 |
| 0.25.0 | 1.0 候选 | ✅ | ✅ **MultiModular Gröbner**（ℚ 理想：并行 F5 幸运素数像 + CRT + 有理重构 + ℚ 精确验证 + 无迹 p-adic Hensel 提升 + 回落）；`Algorithm::Auto` 自动路由；F5 提速（DivisorIndex/分桶 syzygy/并行行构造/两阶段 echelon，cyclic-6 ℤ₁₃ 2.63 s → 1.415 s）；**并行模 GCD**；katsura-6/7 预存在差距记录 |
| 0.26.0 | 1.0 候选 | ✅ | ✅ **打包单项式 F5 快通道**（u128 SWAR，n_vars≤8 且指数<2¹⁵ 自动路由，超界回落）+ echelon i32/免克隆两阶段改造 + **grevlex 基准变体**；修复 Graded 序度方向反置预存在 bug；cyclic-6 ℤ₁₃ grevlex 52.07 ms（criterion 中位数）、Lex 936 ms；cyclic-7 grevlex 单轮 5.755 s（209 基元素） |

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
| Gröbner 基 | F4 真实线性代数（0.15.1）+ F5 签名约简（0.19.0：`Signature`/`SyzygySet` + ℤ_p 原生快速路径 `f5_fp`）+ FGLM + 统一 `groebner_basis()` 分派 + ℤ_p 原生 i64 管线 + MultiModular ℚ 管线（0.25）+ u128 打包 F5 快通道（0.26）；cyclic-6 ℤ₁₃ **55.04 ms** grevlex（2026-08-06 实测 criterion 中位数；0.19.0 基线 2.63 s，约 48×）；cyclic-5 ℤ₁₃ 8.97 ms grevlex | 🟢 F4 + F5 + MultiModular 完成 |
| 符号积分 | Risch（初等超越塔 + RDE 多项式片段）+ 有理函数 Hermite + 三角 exp(I·x) + 特殊函数表（erf/Ei/Si/Ci/Fresnel）+ 0.24 启发式模块（分部/三角换元/Weierstrass/Euler 占位）接入 `try_risch_or_fallback`；回退 `Integral(...)`；**关键差距**：Symbolica 2.2 移植 Rubi 4.17（7000+ 规则、72,944 题库、MIT crate），覆盖面远超 Risch | 🟢 Risch + 启发式完成，Rubi 广度差距大 |
| 实根隔离 | Sturm 序列 + 区间隔离 + refine（单变量）；已知缺口：Wilkinson n=10 展开多项式仅隔离 8/10 根 | 🟡 较完整 |
| 多项式 GCD | GCD + 本原部分 + 扩展 GCD（0.12）+ 经 EEZ 的任意元数多元 GCD（0.16）+ GF(p^d) 上模数域 GCD（CRT + 有理重构，0.17）+ 单变量 Brown 模 GCD 与二元多素数模 GCD（0.21，大整数系数无爆炸） | 🟢 完整（含模快速路径，无 HEVMGCD） |
| 线性求解 | 有理/整数线性方程组 + 二元丢番图（`ax+by=c`） | 🟡 可用，规模有限 |
| JIT 求值 | Cranelift 后端；≥10x 加速目标达成（按路线图标准） | 🟢 完整 |
| 常微分方程 | `ocas-calc::ode`：`dsolve()` 入口 + `classify_ode()` 分类引擎；一阶 5 种（可分离/线性/Bernoulli/恰当/齐次）+ 积分因子检测；二阶 2 种（常系数/Cauchy-Euler）+ VOP + 降阶法；幂级数系数递推 + Frobenius（实有理指标根）；Laplace IVP（`dsolve_ivp`）；2×2 常系数系统（`dsolve_system`）；Python/C 绑定 | 🟢 完整（0.20.1） |
| 数论 | CRT 多模累加、BPSW 素性 + 2⁶⁴ 确定性 MR、整数分解（rho/p−1/p+1/ECM，30 位半素数 1.1 s）、BSGS + Pohlig-Hellman 离散对数、φ/μ/τ/σ_k/λ、二次剩余符号与模平方根；Python/C 绑定（0.21） | 🟢 核心栈完整（0.21） |
| 代数几何工具 | ideal_contains/sum/product/quotient/saturate/intersection、MatrixOrder 消元序 + eliminate()、零维求解、准素分解、根式、Hilbert 级数/维数/次数/多项式、有理根定理（0.23） | 🟢 完整（0.23） |
| 高级模式匹配 | 回溯匹配器 + 多模式替换 + Partition 变换器（0.22）；**差距**：Symbolica `opt` 可选通配符、`alt` 备选、属性过滤更成熟 | 🟡 基本可用，但不及 Symbolica Rubi 级别 |

---

## 4. 与参考系统的差距

### 4.1 对照 Symbolica（Rust，source-available 商业）

> **重大变更**：Symbolica 于 2026 年初从 AGPL-3.0 变更为 source-available 商业许可。
> 免费提供单核非商业用途；商业用途需付费许可。同时拆分出 Numerica（MIT）和
> Graphica（MIT）两个开源 crate。

Symbolica 2.2.0（2026-07-24）是 oCAS 最直接的竞品。2.0 引入可编程符号（normalization/
printing/derivative/series/evaluation hooks）、SymJIT（默认 Python JIT 后端）、
DoubleFloat（~31 位，>3× 快于任意精度）、CUDA/WASM/C++/ASM 代码生成。2.2 移植
Rubi 4.17（7000+ 规则，MIT crate `symbolica-integrate`），是其杀手特性。

| 能力 | oCAS | Symbolica |
|---|---|---|
| 多项式因式分解 | ✅ 一元/多元/代数数域 | ✅ 完整 |
| 有理多项式 | ✅ `RationalPolynomial<D,O>` | ✅ |
| 部分分式/有理重构 | ✅ | ✅ |
| **符号积分** | 🟢 **Risch 完成，覆盖面窄** | ✅ **Rubi 4.17（7000+ 规则、72,944 题、MIT crate）** |
| 数值积分 | ✅ Vegas（0.18） | ✅ |
| 流式 API | ✅ 百万行恒定内存 | ✅ |
| 张量 | ✅ 完整规范化 + Young 投影子（0.22） | ✅ Graphica 图同构引擎 |
| **代码生成** | 🟢 **Cranelift JIT（97×/21×）** | ✅ **SymJIT + CUDA/WASM/C++/ASM** |
| **DoubleFloat** | ✅ **DoubleF64（0.24，~31 位十进制）** | ✅ **~31 位，>3× 快于任意精度** |
| Gröbner 基 | ✅ F5 + u128 打包快通道，cyclic-6 ℤ₁₃ grevlex **55 ms**（2026-08-06 实测） | ✅ 工业级，~1 s |
| ODE 求解器 | ✅ 完整（0.20.1） | 🔴 无 |
| 数论 | ✅ 核心栈完整（0.21） | 🔴 无 |
| 代数几何工具 | ✅ 理想运算+准素分解+Hilbert（0.23） | 🔴 无 |
| 等式饱和 | ✅ egg 集成 | 🔴 无 |
| 模式匹配 | 🟡 回溯匹配器（0.22） | ✅ `opt`/`alt`/属性过滤 |
| 资源控制 | ✅ fuel（0.18） | ✅ |
| **许可证** | ✅ **LGPL-3.0+，可嵌入商业** | ⚠️ **source-available 商业** |

Symbolica 2.2 的核心竞争力——Rubi 积分（7000+ 规则）、SymJIT/CUDA/WASM 代码生成、
DoubleFloat、可编程符号——构成了与 oCAS 的主要功能差距。但 oCAS 在**ODE 求解器**、
**数论**、**代数几何工具**、**等式饱和**上领先（Symbolica 不涉足这些领域），且
**LGPL-3.0 许可证**在商业嵌入场景中是关键差异化优势。

**许可证变更风险**：Symbolica 从 AGPL 变更为 source-available 商业，意味着
竞品格局发生根本性变化——oCAS 现在是唯一同时支持商业嵌入和开源的高性能 Rust CAS。
这一定位差异在 1.0 前后应被明确强调。

### 4.2 对照 SageMath（Python 生态）

SageMath 是"瑞士军刀"式科学计算环境，差距是**广度级**的。

| 领域 | oCAS | SageMath |
|---|---|---|
| 代数几何 | 🟡 基础 Gröbner | ✅ Singular 集成 |
| 数论 | � 核心栈完整（0.21：CRT + 分解 + 素性 + 离散对数 + 数论函数） | ✅ PARI/FLINT 全栈 |
| 微分方程 | 🟢 一阶/二阶/系统/级数/Laplace/绑定 完整（0.20.1） | ✅ 完整 ODE/PDE 求解器 |
| 群论/表示论 | 🔴 无 | ✅ GAP 集成 |
| 组合数学 | 🔴 无 | ✅ 完整 |
| 绘图/可视化 | 🔴 无 | ✅ matplotlib 集成 |
| 数据库接口 | 🔴 无 | ✅ OEIS / LMFDB |

SageMath 通过包装 80+ 专用库实现广度；oCAS 是自包含内核。两者定位不同——
oCAS 面向高性能**库**，SageMath 面向完整**环境**。可比性主要集中在核心
代数性能，而非功能广度。

### 4.3 对照 SymPy 1.14（纯 Python）

SymPy 1.14.0（2025-04-27）引入 DomainMatrix（FLINT 后端）大幅提升矩阵性能。

| 领域 | oCAS vs SymPy | 说明 |
|---|---|---|
| 解析/化简 | 🟢 持平 | 双方都完备 |
| 微分 | 🟢 持平 | 链式/乘积/幂法则 |
| 积分 | 🟢 基本持平 | 双方均有 Risch；SymPy 的 `manualintegrate` 启发式覆盖更广 |
| 因式分解 | 🟢 持平 | oCAS 任意多元+代数数域已补齐 |
| Gröbner | 🟢 oCAS 优势 | oCAS F5 打包快通道 cyclic-6 grevlex 55 ms（实测），SymPy 仅 Buchberger |
| **矩阵/线性代数** | 🟡 **SymPy 追赶** | SymPy 1.14 DomainMatrix + FLINT 后端 rref **10000× 加速** + Smith 标准形；oCAS 仍为 Bareiss |
| 数论 | 🟢 持平 | oCAS ECM，SymPy 二次筛（大整数可能更快） |
| **性能** | 🟢 **oCAS 优势** | Rust + JIT 对纯 Python；2026-08-06 实测：parse 100×、simplify 124×、series 2,550×、integrate 39–76×、eval 183×；**例外**：`factor(x^30-1)` SymPy 快约 50×（分圆幂的快路径） |
| Python 易用性 | 🟢 持平 | oCAS 有 `ocas-py` 绑定 |

0.6.0 成功标准“基础多项式/微积分/重写与 SymPy 持平”——在**性能**维度已
达成并领先，**积分**经 0.14 Risch 补齐，**因式分解**经 0.16 达任意多元
持平（0.17 另补 Trager 代数数域）。与 SymPy 的剩余功能差距集中在**积分
启发式回退的广度**（SymPy 的 `manualintegrate`/启发式池比 oCAS 的
Risch + 查表路径更宽）和**矩阵/线性代数**（DomainMatrix 10000× 加速后
差距扩大，Smith 标准形 oCAS 缺失）。

### 4.4 对照 msolve（Gröbner 性能标杆）

msolve 是 Gröbner 基计算的开源性能标杆，使用 F4/F5 + 多模算术 + Hensel 提升 +
Berlekamp-Massey 优化。**2026-08-06 本机实测**（WSL2，msolve 0.10.1 源码构建，
`-g 2` 仅 GB 模式、DRL 序、单线程；基元素数与 oCAS 逐项一致）：

| 基准 | msolve 0.10.1（实测） | oCAS 0.26（实测） | 比值 |
|---|---|---|---|
| cyclic-5 ℤ₁₃ | 3 ms | 8.97 ms（grevlex） | 3.0× |
| cyclic-6 ℤ₁₃ | **4 ms** | **55.04 ms**（grevlex） | 13.8× |
| cyclic-7 ℤ₁₃ | **55 ms** | 3.829 s（grevlex 单轮） | ~70× |
| katsura-6 | 3 ms | 未测（预存在差距） | — |
| katsura-7 | 7 ms | 未测（预存在差距） | — |

**差距根源**：msolve 的多模算术 + Hensel 提升 + BM 策略在大规模 Gröbner 上仍领先
oCAS 1–2 个数量级。0.26.0 的打包单项式 F5 已将 cyclic-6 拉到 55 ms（0.19.0 基线
2.63 s，约 48×），但 cyclic-7（~70×）与 katsura 系（oCAS 未完成 vs msolve 3–7 ms）
仍是最大单项差距。上一版文档引用的 msolve 数值（cyclic-6 0.04 s）已被本轮实测
（4 ms）替换。

**Windows 可用性**：msolve 可通过 MSYS2 在 Windows 上编译（需 GMP/MPFR/FLINT），
但非原生体验。oCAS 在 Windows 上原生支持是差异化优势。

### 4.5 新兴竞品

| 竞品 | 语言 | 许可证 | 功能 | 对 oCAS 的威胁 |
|---|---|---|---|---|
| **Numerica** | Rust | MIT | 高性能数类型、误差跟踪浮点、有限域、矩阵、Vegas 积分、双数（2026-08-06 核实：无 release/tag，开发活跃） | 直接竞争 `ocas-domain`/`ocas-eval` 的数值层 |
| **Graphica** | Rust | MIT | 图同构（McKay 算法）、Feynman 图生成 | 竞争 `ocas-rewrite` 的张量规范化 |
| **mathcore** | Rust | MIT | 泛型代数结构（环/域/多项式/矩阵）；2026-08-06 核实：最新 **0.3.1**（2025-08-30，crates.io），上一版误记 0.5.0 | 底层代数基础设施竞品，功能远不及 oCAS |
| **cas-rs** | Rust | — | Rust CAS 工作进行中 | 功能远不及 oCAS，但验证 Rust CAS 市场需求 |

**评估**：Numerica 和 Graphica 是 Symbolica 的 MIT 拆分，功能专精但不构成完整
CAS 威胁。mathcore/cas-rs 是早期项目。oCAS 的自包含内核 + LGPL 许可证 + 三语言
绑定形成了明确的差异化定位。

---

## 5. 关键缺口与优先级（2026-08-06 重排）

基于本轮（2026-08-06）竞品版本重新核实 + 本机全量复测（msolve 0.10.1 实测、
oCAS/Symbolica/SymPy 基准、cyclic-6 grevlex 55 ms 达成 <0.5 s 里程碑），
重排 1.0.0 前的优先级。阶段 B++ “竞品全面对齐”（0.19–0.23）+ 0.24–0.26
（启发式积分/DoubleF64、MultiModular、打包 F5）已完成。

| 优先级 | 缺口 | 现状 | 目标 | 理由 |
|---|---|---|---|---|
| **P0** | 符号积分广度（Rubi 规则集成或等效） | Risch + 0.24 启发式四技术，覆盖面仍窄 | 对标 symbolica-integrate 1892 题 | Symbolica 2.2 Rubi 7000+ 规则（72,944 题库）仍是最大功能缺口 |
| **P1** | Gröbner 大规模性能（katsura 系 + cyclic-7） | katsura-6/7 未完成（单轮 >30 min）；cyclic-7 Lex >2 h 未完成；cyclic-7 grevlex 3.829 s vs msolve 55 ms（~70×） | katsura-6 < 1 s；cyclic-7 可完成 | msolve 0.10.1 实测 katsura 3–7 ms、cyclic-7 55 ms；打包管线 + 多模策略向 katsura/cyclic-7 扩展 |
| **P1** | 代码生成扩展（LLVM JIT + CUDA/WASM 导出） | 仅 Cranelift JIT | 至少 LLVM JIT | Symbolica SymJIT/CUDA/WASM 形成代差 |
| **P2** | 矩阵/线性代数增强 | Bareiss 行列式/逆 | DomainMatrix 类似引擎 + Smith 标准形 | SymPy 1.14 DomainMatrix 10000× 加速后差距扩大 |
| **P2** | Windows FLINT 支持 | `flint` feature 仅 Linux/WSL | 三平台可用 | 平台覆盖完整性 |
| **P3** | 张量嵌套函数内处理 | 基础版 | 对标 Symbolica Graphica | 0.22 已补齐基础规范化 |
| **P3** | 大整数分解（二次筛） | ECM 30 位 1.1 s | 对标 SymPy qs_factor | SymPy 二次筛在大整数上可能更快 |

### 已完成项（历史记录）

| # | 缺口 | 完成版本 |
|---|---|---|
| 1 | 完整多项式因式分解 | 0.11.0–0.11.1 |
| 2 | Risch 符号积分 | 0.14.0 |
| 3 | Gröbner F4/F5 | 0.13.0/0.14.0/0.15.1/0.19.0 |
| 4 | 有理多项式/部分分式 | 0.12.0 |
| 5 | 多输出优化/代码生成 | 0.15.0 |
| 6 | Gröbner 大规模性能 cyclic-6 < 5 s | 0.19.0 (2.63 s) |
| 7 | 任意多元因式分解 | 0.16.0–0.16.2 |
| 8 | 代数数域因式分解 | 0.17.0–0.17.1 |
| 9 | 数值积分/双数/张量基础/fuel | 0.18.0–0.18.1 |
| 10 | ODE 求解器 | 0.20.0–0.20.1 |
| 11 | 数论栈 | 0.21.0 |
| 12 | 张量规范化+高级模式匹配 | 0.22.0 |
| 13 | 代数几何工具 | 0.23.0 |
| 14 | PDE 求解器 | Post-1.0 |
| 15 | **DoubleFloat（DoubleF64，~31 位扩展精度）** | **0.24.0**（原 P2，本轮移入） |
| 16 | **Gröbner cyclic-6 < 0.5 s**（grevlex 实测 55.04 ms） | **0.26.0**（原 P1，本轮移入；msolve 比值 13.8× 另行跟踪） |

---

## 6. 总评

0.1 → 0.26.0 执行质量极高：每个路线图交付物均兑现，分层架构干净（无环依赖），
13 crate workspace 严格分层，质量门严格（`-D warnings` + deny + Miri 意识），
文档/绑定/CI 工程化完备。阶段 A（Beta 硬代数）、阶段 B+（Symbolica 差距清零）、
阶段 B++（竞品全面对齐）、0.24–0.26（启发式积分/DoubleF64、MultiModular、
打包 F5）全部完成。2026-08-06 全量复测确认：cyclic-6 ℤ₁₃ grevlex 55.04 ms
（0.19.0 基线的 ~48×），达成「< 0.5 s」里程碑。

**务实定位**：当前 oCAS 是“高性能、自包含的代数内核，功能对标 SymPy，许可证
优于 Symbolica”。核心优势：

1. **唯一 LGPL-3.0 高性能 Rust CAS**——Symbolica 已转为 source-available 商业
2. **功能广度领先 Symbolica**——ODE、数论、代数几何、等式饱和（Symbolica 不涉足）
3. **性能全面领先 SymPy**——Rust + JIT 对纯 Python（parse 100×、series 2,550× 等；
   例外：`factor(x^n−1)` 类分圆输入 SymPy 快约 50×）
4. **三语言绑定**——Rust + Python + C/C++，比 Symbolica 更广

**关键差距**（需 1.0.0 前或紧随其后解决）：

1. 符号积分广度——Symbolica Rubi 7000+ 规则 vs oCAS Risch + 0.24 启发式
2. Gröbner 大规模性能——katsura-6/7 未完成、cyclic-7 grevlex ~70× msolve
   （msolve 0.10.1 实测 4–55 ms；cyclic-6 已收敛至 55 ms / 13.8×）
3. 代码生成——Symbolica CUDA/WASM/C++ vs oCAS 仅 Cranelift
4. 矩阵/线性代数——SymPy DomainMatrix 10000× 加速

---

## 7. 许可证与生态位分析

### 7.1 许可证格局变更

| 系统 | 原许可证 | 现许可证 | 变更日期 | 影响 |
|---|---|---|---|---|
| Symbolica | AGPL-3.0 | source-available 商业 | 2026 年初 | 商业用途需付费；单核非商业免费 |
| Numerica | （从 Symbolica 拆分） | MIT | 2025-11 | 可自由嵌入商业项目 |
| Graphica | （从 Symbolica 拆分） | MIT | 2025-11 | 可自由嵌入商业项目 |
| symbolica-integrate | （从 Symbolica 拆分） | MIT | 2026-07 | Rubi 规则的 MIT Rust 移植 |

### 7.2 oCAS 的许可证优势

**LGPL-3.0+** 允许：
- ✅ 嵌入商业项目（无需开源商业代码）
- ✅ 静态/动态链接
- ✅ 修改 oCAS 本身需回馈（copyleft 仅限库本身）
- ✅ 与 MIT/Apache/BSD 项目兼容

**对比**：
- Symbolica：source-available 商业——商业用途需付费许可
- SymPy：BSD-3-Clause——最宽松，但纯 Python 性能差
- SageMath：GPL-2.0+——GPL 传染，不可嵌入商业项目
- msolve：GPL-3.0——GPL 传染
- GiNaC：GPL-2.0+——GPL 传染
- FLINT：LGPL-3.0+——与 oCAS 相同，但仅是多项式/数论库

**结论**：oCAS 是唯一同时满足“高性能 Rust CAS”和“可嵌入商业项目”的系统。
这是 1.0 前后应明确强调的核心差异化定位。

### 7.3 Rubi 规则集许可证风险

symbolica-integrate 是 MIT 许可证的 Rust 移植，但原始 Rubi 规则来自 Wolfram
Language（Mathematica）。Rubi 本身是开源的（CC BY-NC-SA 3.0），但其许可证
与商业使用存在潜在争议。建议：

1. **方案 A**：集成 symbolica-integrate（MIT crate）——最直接，但需评估 Rubi
   原始规则的许可证兼容性
2. **方案 B**：自研规则集——工作量大但完全可控
3. **方案 C**：混合——Risch 算法 + 启发式扩展 + 少量 Rubi 规则作为参考

---

## 8. 战略建议

### 8.1 1.0.0 前（冻结前）

| 建议 | 版本 | 交付物 | 理由 |
|---|---|---|---|
| 积分广度扩展 | 1.0-rc | Risch + 启发式扩展 或 Rubi 规则集成 | Symbolica 2.2 杀手特性；覆盖面差距是最大功能缺口 |
| Gröbner katsura/cyclic-7 扩展 | 1.0-rc | katsura-6 < 1 s；cyclic-7 可完成 | msolve 实测 3–55 ms；打包 F5 已验证 cyclic-6 收敛 |
| 矩阵增强 | 1.0-rc | DomainMatrix 类似引擎 | SymPy 1.14 差距扩大 |

> 已兑现（0.24–0.26）：DoubleFloat（→DoubleF64）、Gröbner cyclic-6 < 0.5 s
> （grevlex 55 ms）、MultiModular ℚ 管线、启发式积分四技术。

### 8.2 Post-1.0

| 建议 | 优先级 | 理由 |
|---|---|---|
| LLVM/inkwell JIT 后端 | P1 | Symbolica SymJIT 代差；Cranelift 已到性能天花板 |
| CUDA/WASM 代码导出 | P1 | Symbolica 已支持；GPU/浏览器场景需求 |
| PDE 求解器 | P2 | 用户期望高；Poisson/热传导/波动 |
| 二次筛分解 | P2 | SymPy qs_factor 可能更快 |
| Windows FLINT | P2 | 平台覆盖完整性 |

### 8.3 定位建议

**一句话定位**：
> oCAS 是唯一同时支持商业嵌入的高性能 Rust CAS，功能广度领先 Symbolica，
> 性能全面领先 SymPy。

**差异化强调**：
1. 许可证优势（LGPL-3.0 vs Symbolica 商业化）
2. 功能广度（ODE + 数论 + 代数几何 + 等式饱和，Symbolica 无）
3. 三语言绑定（Rust + Python + C/C++）
4. Windows 原生支持（msolve/SageMath 需 WSL）

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
| 0.20.0 | 2026-07-27 | **常微分方程求解器发布。** §3 新增 ODE 行（🟡）。§4.1 新增 ODE 竞品行（🟡）。§5 #10（ODE 求解器）标记 🟡 部分完成——一阶 5 种（可分离/线性/Bernoulli/恰当/齐次）+ 二阶 2 种（常系数/Cauchy-Euler）+ 幂级数框架 + ODE 分类引擎 `classify_ode()` + `dsolve()` 入口；Laplace 变换、ODE 系统、Python/C 绑定推迟。版本提升 0.20.0。 |
| 0.20.1 | 2026-07-27 | **ODE 求解器全量收尾。** 积分因子检测（μ(x)/μ(y)）；常数变易法（VOP，修复 Cauchy-Euler forcing 静默丢弃）；降阶法；幂级数系数递推 + Frobenius（实有理指标根）；Laplace IVP（`dsolve_ivp`）；2×2 常系数系统（`dsolve_system`）；Python/C 绑定（`classify_ode`/`dsolve`/`dsolve_ivp`）。修复 real_roots isqrt/公式 bug、is_exact 硬编码、Cauchy-Euler 系数归一化、积分器 (ax+b)^-1 与分数次幂、substitute_solution 裸 y(x)、级数系数 diff 污染。新增同类项收集器 collect_terms + expand。31 项代入验证正确性测试（3 项已知限制 ignore）。ODE 缺口从 🟡 升级为 🟢。版本提升 0.20.1。 |
| 0.23.0 竞品调研 | 2026-08-03 | **全面竞品差距调研与评估。** 头部更新至 0.23.0。新增配套文档 COMPETITIVE_MATRIX_CN.md（竞品能力矩阵）+ BENCHMARK_SUITE_CN.md（基准测试套件设计）。§1 版本表扩展至 0.23.0。§3 算法深度新增代数几何工具、高级模式匹配行，符号积分行标注 Rubi 广度差距。§4.1 Symbolica 对照重写（2.2.0 source-available + Rubi + SymJIT/CUDA/WASM + DoubleFloat）；新增 §4.4 msolve 对照（cyclic-6 0.04 s 标杆）；新增 §4.5 新兴竞品（Numerica/Graphica/mathcore/cas-rs）。§5 优先级重排（P0 积分广度、P1 Gröbner+代码生成、P2 矩阵+DoubleFloat+FLINT、P3 张量+二次筛）。新增 §7 许可证与生态位分析（Symbolica 许可证变更 + Rubi 许可证风险 + oCAS LGPL 优势）。新增 §8 战略建议（1.0 前 + Post-1.0 + 定位建议）。§6 总评重写。 |
| 0.26.0 复测 | 2026-08-06 | **竞品版本重新核实 + 本机全量复测。** 头部更新至 0.26.0。竞品核实：FLINT 3.5.0→3.6.0（Kinoshita-Li 级数复合、padic_radix、subresultant 结式）、msolve 0.7.x→0.10.1（GM 改进、QQ 提升修复）、mathcore 更正为 0.3.1（0.5.0 不存在）、Numerica 无 tag、SageMath 日期更正 10.9@2026-05-05；Symbolica/SymPy/GiNaC 无更新。§1 版本表追加 0.24.0/0.25.0/0.26.0 行。§4.1 DoubleFloat 行 🔴→✅（0.24 DoubleF64）；§4.4 msolve 表以 WSL2 实测值替换引用值（cyclic-6 4 ms、cyclic-7 55 ms、katsura 3–7 ms）。§5 优先级重排：DoubleFloat 与 cyclic-6<0.5 s（grevlex 55.04 ms 实测）移入已完成项，Gröbner 剩余差距重定为 katsura+cyclic-7。§6/§8 同步改写（含 factor(x^n−1) 类 SymPy 快 ~50× 的诚实记录）。 |