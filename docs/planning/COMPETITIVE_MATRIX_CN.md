# oCAS 0.23.0 竞品能力矩阵

> **编制日期**：2026-08-03
> **数据来源**：Symbolica 官方博客（symbolica.io）、SymPy 1.14 文档、SageMath GitLab 标签、
> FLINT GitHub/changelog、msolve GitHub/Groebner.jl 论文、GiNaC 官网/Codeberg、
> Numerica GitHub、mathcore crates.io
> **评估基准**：oCAS 0.23.0（2026-08-02 发布）

---

## 图例

| 标记 | 含义 |
|---|---|
| ✅ | 完成 / 工业级 |
| 🟡 | 基础可用 / 部分完成 |
| 🔴 | 缺失 / 重大缺口 |
| ⚠️ | 完成但有保留（平台限制等） |

---

## 1. 竞品版本快照

| 系统 | 最新版本 | 发布日期 | 语言 | 许可证 | 定位 |
|---|---|---|---|---|---|
| **oCAS** | 0.23.0 | 2026-08-02 | Rust | LGPL-3.0+ | 高性能自包含 CAS 内核 |
| **Symbolica** | 2.2.0 | 2026-07-24 | Rust | source-available 商业（单核免费非商业） | 高性能符号框架 |
| **SymPy** | 1.14.0 | 2025-04-27 | Python | BSD-3-Clause | 纯 Python CAS |
| **SageMath** | 10.9 | 2026-05-04 | Python/Cython | GPL-2.0+ | 科学计算环境（80+ 库集成） |
| **FLINT** | 3.5.0 | 2026-04-24 | C | LGPL-3.0+ | 多项式/数论性能库 |
| **msolve** | 0.7.x | 2026 | C | GPL-3.0 | Gröbner 基性能标杆 |
| **Numerica** | — | 2025-11 | Rust | MIT | Symbolica 数值拆分 |
| **GiNaC** | 1.8.10 | 2026-02-11 | C++ | GPL-2.0+ | C++ 嵌入式符号库 |
| **mathcore** | 0.5.0 | 2026-03 | Rust | MIT | Rust CAS 代数基础 |

---

## 2. 逐维度能力矩阵

### 2.1 基础代数（多项式运算、GCD、因式分解、结式）

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | FLINT 3.5 | msolve | GiNaC 1.8.10 |
|---|---|---|---|---|---|---|---|
| 稠密单变量多项式 | ✅ | ✅ | ✅ | ✅ | ✅ 性能标杆 | — | ✅ |
| 稀疏多元多项式 | ✅ | ✅ | ✅ | ✅ (Singular) | ✅ | — | 🟡 有限 |
| 多项式 GCD（ℤ） | ✅ 模方法 | ✅ | ✅ | ✅ FLINT | ✅ 性能标杆 | — | ✅ |
| 多项式 GCD（ℤ_p） | ✅ 模方法 | ✅ | ✅ | ✅ | ✅ | — | ✅ |
| 多元 GCD | ✅ Brown 模 GCD | ✅ | ✅ | ✅ | ✅ | — | 🟡 |
| 因式分解（一元 ℤ/ℤ_p） | ✅ CZ+Hensel | ✅ | ✅ | ✅ FLINT | ✅ | — | ✅ |
| 因式分解（多元） | ✅ Wang EEZ | ✅ | ✅ | ✅ Singular | ✅ | — | 🔴 |
| 因式分解（代数数域） | ✅ Trager | ✅ | 🟡 | ✅ | — | — | 🔴 |
| 结式 | ✅ Brown PRS | ✅ | ✅ | ✅ | ✅ | — | ✅ |
| 无平方分解 | ✅ Yun SFF | ✅ | ✅ | ✅ | ✅ | — | ✅ |
| 有理函数 | ✅ RationalPolynomial | ✅ | ✅ | ✅ | — | — | 🟡 |
| 部分分式 | ✅ | ✅ | ✅ | ✅ | — | — | ✅ |

**评估**：oCAS 0.23 在多项式代数上已达到与 Symbolica、SymPy 持平的水平。FLINT 在纯
数值性能上仍是天花板，但 oCAS 已有 `flint` feature 集成。SageMath 通过 Singular/FLINT
组合覆盖更广。

### 2.2 Gröbner 基

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | msolve | Singular (via SageMath) |
|---|---|---|---|---|---|---|
| F4 算法 | ✅ | ✅ | 🔴 | ✅ (Singular) | ✅ 多模+Hensel | ✅ |
| F5 签名约简 | ✅ f5_fp | ✅ | 🔴 | 🟡 | ✅ F4/F5 | 🟡 |
| Buchberger 经典 | ✅ | ✅ | ✅ | ✅ | — | ✅ |
| 多单项式序 | ✅ Grlex/Weight/Block | ✅ | ✅ | ✅ | ✅ | ✅ |
| FGLM 换序 | ✅ | ✅ | — | ✅ | ✅ | ✅ |
| Hilbert 级数 | ✅ | ✅ | — | ✅ | — | ✅ |
| cyclic-6 ℤ₁₃ | **2.63 s** | ~1 s | 🔴 | ~3.4 s (Oscar) | **0.04 s** | ~0.5 s |
| cyclic-7 ℤ₁₃ | 🔴 未测 | ✅ | 🔴 | ✅ | **~1 s** | ✅ |
| 理想运算 | ✅ sum/product/quotient/saturate | ✅ | 🟡 | ✅ | — | ✅ |
| 准素分解 | ✅ 零维 | ✅ | 🟡 | ✅ | — | ✅ |

**评估**：oCAS F5 在 cyclic-6 上 2.63 s，与 Symbolica ~1 s 在同一量级，但 msolve 的
0.04 s（多模+Hensel+BM 优化）领先约 65×。cyclic-7 需要 multi-modular 策略进一步
加速。oCAS 理想运算和准素分解已补齐。

### 2.3 符号积分

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | Mathematica (参考) |
|---|---|---|---|---|---|
| Risch 算法 | ✅ 初等超越塔+RDE | 🟡 部分 | ✅ | ✅ (Maxima) | ✅ |
| Rubi 规则集 | 🔴 | ✅ 7000+ 规则 (MIT crate) | 🔴 | 🔴 | ✅ 原生 |
| 启发式查表 | ✅ 特殊函数表 | ✅ | ✅ manualintegrate | ✅ | ✅ |
| 有理函数积分 | ✅ Hermite+Rothstein-Trager | ✅ Trager | ✅ | ✅ | ✅ |
| 三角积分 | ✅ exp(I·x)+realify | ✅ | ✅ | ✅ | ✅ |
| 特殊函数 | ✅ erf/Ei/Si/Ci/Fresnel | ✅ gamma/polylog/Bessel/zeta | ✅ | ✅ | ✅ |
| 步骤推导 | 🔴 | ✅ rule-by-rule trace | 🔴 | 🔴 | ✅ |
| Rubi 语料库覆盖 | — | ✅ 72,944 题 99.9% | — | — | ✅ |

**评估**：这是 oCAS 与 Symbolica 2.2 最大的功能差距。Symbolica 的 Rubi 移植
（7000+ 规则、72,944 题库、MIT crate `symbolica-integrate`）是其杀手特性。
oCAS 的 Risch 算法在理论上更完备（decidable），但实际覆盖面远窄于 Rubi。
SymPy 的 `manualintegrate` 启发式覆盖面也比 oCAS 广。

### 2.4 常微分方程

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 |
|---|---|---|---|---|
| 一阶（可分离/线性/Bernoulli/恰当） | ✅ 5 种 | 🟡 | ✅ | ✅ |
| 二阶（常系数/Cauchy-Euler） | ✅ | 🔴 | ✅ | ✅ |
| 降阶法/待定系数/VOP | ✅ | 🔴 | ✅ | ✅ |
| 级数解（幂级数/Frobenius） | ✅ | 🔴 | ✅ | ✅ |
| Laplace 变换 IVP | ✅ | 🔴 | ✅ | ✅ |
| ODE 系统 | ✅ 2×2 | 🔴 | ✅ | ✅ |
| PDE 求解 | 🔴 Post-1.0 | 🔴 | 🟡 | ✅ |

**评估**：oCAS 在 ODE 上已超越 Symbolica（Symbolica 未投入 ODE）并与 SymPy 持平。
SageMath 通过 Maxima/desolve 覆盖更广。PDE 是 Post-1.0 议题。

### 2.5 数论

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | FLINT 3.5 |
|---|---|---|---|---|---|
| 素性检测 | ✅ BPSW+确定性 MR | 🔴 | ✅ BPSW | ✅ PARI | ✅ |
| 整数分解 | ✅ rho/p±1/ECM | 🔴 | ✅ rho/p-1/二次筛 | ✅ PARI/ECM/QS | ✅ |
| 离散对数 | ✅ BSGS+Pohlig-Hellman | 🔴 | ✅ BSGS+index calc | ✅ PARI | ✅ |
| 数论函数 | ✅ φ/μ/τ/σ_k/λ | 🔴 | ✅ | ✅ | ✅ |
| CRT 多模 | ✅ | 🔴 | ✅ | ✅ | ✅ |
| 模 GCD | ✅ Brown 单变量+二元 | ✅ | 🟡 单变量 | ✅ FLINT | ✅ 性能标杆 |
| 椭圆曲线 | 🔴 | 🔴 | 🟡 | ✅ PARI | 🟡 |
| 二次筛 | 🔴 | 🔴 | ✅ qs_factor | ✅ | 🔴 |

**评估**：oCAS 数论核心栈完整（0.21），但 SymPy 的 `qs_factor()`（二次筛）在
大整数分解上可能优于 oCAS 的 ECM。SageMath 通过 PARI 2.17 提供完整的椭圆曲线/
类域论。Symbolica 不涉足数论。

### 2.6 求值与 JIT

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SymJIT (Symbolica) |
|---|---|---|---|---|
| 解释器求值 | ✅ 树解释器 | ✅ | ✅ Python | ✅ |
| JIT 后端 | ✅ Cranelift | ✅ SymJIT (默认 Python) | 🔴 | ✅ |
| SIMD 向量化 | ✅ AVX2 via pulp | ✅ | 🔴 | ✅ |
| 多输出批处理 | ✅ compile_multi | ✅ | 🔴 | ✅ |
| f32 混合精度 | ✅ | ✅ | 🔴 | ✅ |
| DoubleFloat (~31 位) | 🔴 | ✅ >3× 快于任意精度 | 🔴 | ✅ |
| 流式求值 | ✅ 百万行恒定内存 | ✅ | 🔴 | ✅ |
| 常量折叠+栈压缩 | ✅ | ✅ | 🔴 | ✅ |
| CUDA 代码生成 | 🔴 Post-1.0 | ✅ | 🔴 | ✅ |
| WASM 代码生成 | 🔴 Post-1.0 | ✅ | 🔴 | ✅ |
| C++/ASM 导出 | 🔴 Post-1.0 | ✅ | 🔴 | ✅ |
| LLVM 后端 | 🔴 Post-1.0 (inkwell) | 🔴 | 🔴 | 🔴 |

**评估**：Symbolica 2.0 的 SymJIT（默认 Python 后端）+ CUDA/WASM/C++/ASM 导出形成
显著代差。oCAS 的 Cranelift JIT 在单/多输出上表现良好（97×/21×），但缺少
CUDA/WASM 导出和 DoubleFloat。DoubleFloat 是低垂果实（~31 位，>3× 快于任意精度）。

### 2.7 张量

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 |
|---|---|---|---|---|
| 张量类型+指标槽 | ✅ | ✅ | ✅ | ✅ |
| 显式收缩 | ✅ | ✅ | ✅ | ✅ |
| 对称化/反对称化 | ✅ symmetrise_sign | ✅ | ✅ | ✅ |
| Young 投影子 | ✅ 0.22 | ✅ | 🟡 | ✅ |
| 图同构正则标号 | ✅ (自研) | ✅ (Graphica) | 🔴 | 🔴 |
| 完整规范形 | 🟡 基础版 | ✅ (Graphica 图同构) | 🟡 | ✅ |
| 嵌套函数内张量 | 🟡 | ✅ | 🟡 | ✅ |

**评估**：oCAS 0.22 已补齐 Young 投影子和图同构正则标号，但 Symbolica 通过
Graphica（MIT）提供更完整的图同构引擎。差距在嵌套函数内张量处理和完整规范形。

### 2.8 线性代数

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | FLINT 3.5 |
|---|---|---|---|---|---|
| 稠密矩阵 | ✅ | ✅ | ✅ DomainMatrix | ✅ | ✅ |
| 行列式 | ✅ Bareiss | ✅ | ✅ | ✅ | ✅ |
| RREF | ✅ | ✅ | ✅ 10000× 加速 | ✅ | ✅ |
| 矩阵求逆 | ✅ | ✅ | ✅ | ✅ | ✅ |
| 特征值 | 🔴 | 🔴 | 🟡 数值 | ✅ | ✅ |
| Smith 标准形 | 🔴 | 🔴 | ✅ | ✅ | ✅ |
| 稀疏矩阵 | 🔴 | 🔴 | ✅ DomainMatrix | ✅ | ✅ |
| 线性系统求解 | ✅ ℚ/ℤ | ✅ | ✅ | ✅ | ✅ |

**评估**：SymPy 1.14 的 DomainMatrix 引入 FLINT 后端后 rref 性能提升 10000×，
Smith 标准形也已支持。oCAS 在线性代数上规模有限，缺少特征值和 Smith 标准形。

### 2.9 模式匹配与重写

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 |
|---|---|---|---|
| 通配符匹配 | ✅ | ✅ | ✅ |
| 条件守卫 | ✅ | ✅ | ✅ |
| 交换/结合匹配 | ✅ | ✅ | 🟡 |
| 回溯匹配器 | ✅ 0.22 | ✅ | 🔴 |
| 多模式替换 | ✅ 0.22 | ✅ | 🟡 |
| 可选通配符 (opt) | 🔴 | ✅ | 🔴 |
| 模式备选 (alt) | 🔴 | ✅ | 🔴 |
| 属性过滤 | 🟡 | ✅ | 🟡 |
| Partition 变换器 | ✅ 0.22 | ✅ | 🔴 |
| egg 等式饱和 | ✅ | 🔴 | 🔴 |

**评估**：Symbolica 2.0 的模式匹配更成熟（`opt` 可选通配符、`alt` 备选、属性过滤），
支撑了 Rubi 7000+ 规则的复杂匹配需求。oCAS 的回溯匹配器是 0.22 新实现，
功能覆盖度仍有差距。但 oCAS 独有 egg 等式饱和。

### 2.10 代码生成与绑定

| 能力 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | GiNaC 1.8.10 |
|---|---|---|---|---|---|
| Rust API | ✅ | ✅ | 🔴 | 🔴 | 🔴 |
| Python API | ✅ PyO3 | ✅ PyO3 | ✅ 原生 | ✅ 原生 | 🔴 |
| C API | ✅ cbindgen | 🔴 | 🔴 | 🔴 | ✅ |
| C++ API | ✅ RAII 包装 | 🔴 | 🔴 | 🔴 | ✅ 原生 |
| C 代码生成 | 🔴 | ✅ | 🔴 | 🔴 | 🔴 |
| C++ 代码生成 | 🔴 | ✅ | 🔴 | 🔴 | 🔴 |
| CUDA 代码生成 | 🔴 | ✅ | 🔴 | 🔴 | 🔴 |
| WASM 代码生成 | 🔴 | ✅ | 🔴 | 🔴 | 🔴 |
| ASM 代码生成 | 🔴 | ✅ | 🔴 | 🔴 | 🔴 |
| PyPI 分发 | ✅ maturin wheels | ✅ | ✅ | ✅ | 🔴 |
| crates.io | ✅ | ✅ source-available | 🔴 | 🔴 | 🔴 |

**评估**：oCAS 的三语言绑定（Rust/Python/C/C++）与 Symbolica 持平甚至更广（C/C++）。
但 Symbolica 的代码生成目标远超 oCAS（C++/ASM/SIMD/CUDA/WASM）。oCAS 的
LGPL-3.0 许可证在嵌入商业项目方面有优势。

### 2.11 许可证与生态

| 维度 | oCAS 0.23 | Symbolica 2.2 | SymPy 1.14 | SageMath 10.9 | FLINT 3.5 | msolve | GiNaC 1.8.10 | Numerica |
|---|---|---|---|---|---|---|---|---|
| 许可证 | LGPL-3.0+ | source-available 商业 | BSD-3-Clause | GPL-2.0+ | LGPL-3.0+ | GPL-3.0 | GPL-2.0+ | MIT |
| 商业嵌入 | ✅ | ⚠️ 需付费许可 | ✅ | 🔴 GPL 传染 | ✅ | 🔴 GPL 传染 | 🔴 GPL 传染 | ✅ |
| 免费使用 | ✅ | ⚠️ 单核非商业 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Rust 生态 | ✅ crates.io | ✅ crates.io | 🔴 | 🔴 | 🟡 flint3-sys | 🔴 | 🔴 | ✅ crates.io |

**评估**：oCAS 的 LGPL-3.0 在商业嵌入场景中是关键差异化优势——Symbolica 已转为
source-available 商业模式。Numerica（MIT）和 mathcore（MIT）是新兴 Rust 竞品，
但功能远不及 oCAS。

---

## 3. 总结：oCAS 竞争位置

### 优势领域
1. **许可证**：LGPL-3.0 是唯一同时支持商业嵌入和开源的高性能 Rust CAS
2. **ODE 求解器**：超越 Symbolica，与 SymPy 持平
3. **三语言绑定**：Rust + Python + C/C++，比 Symbolica 更广
4. **Gröbner F5**：cyclic-6 2.63 s，与 Symbolica ~1 s 同一量级
5. **数论核心栈**：完整的 CRT + 分解 + 素性 + 离散对数
6. **代数几何工具**：理想运算 + 准素分解 + Hilbert 级数（0.23）

### 关键差距（按优先级）
1. **P0 — 符号积分广度**：Symbolica Rubi 7000+ 规则 vs oCAS Risch（覆盖面窄）
2. **P1 — Gröbner 大规模性能**：msolve 0.04 s vs oCAS 2.63 s（~65×）
3. **P1 — 代码生成扩展**：Symbolica CUDA/WASM/C++/ASM vs oCAS 仅 Cranelift JIT
4. **P2 — 矩阵/线性代数**：SymPy DomainMatrix 10000× 加速 + Smith 标准形
5. **P2 — DoubleFloat**：Symbolica ~31 位 >3× 快于任意精度
6. **P3 — 张量完整规范形**：Symbolica Graphica 图同构引擎更成熟
7. **P3 — 大整数分解**：SymPy 二次筛 vs oCAS ECM

### 无差距 / oCAS 领先
- ODE 求解器（Symbolica 无）
- 数论（Symbolica 无）
- 代数几何工具（Symbolica 无）
- egg 等式饱和（Symbolica 无）
- 许可证灵活性（Symbolica 商业化）
