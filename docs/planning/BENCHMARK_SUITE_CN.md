# oCAS 跨 CAS 基准测试套件设计

> **编制日期**：2026-08-03
> **关联文档**：[COMPETITIVE_MATRIX_CN.md](COMPETITIVE_MATRIX_CN.md)、[GAP_ANALYSIS_CN.md](GAP_ANALYSIS_CN.md)

---

## 1. 目标

设计一套**可复现、跨 CAS** 的基准测试套件，量化 oCAS 与竞品在核心算法维度的
性能差距，为 1.0.0 前的优先级调整提供数据支撑。

**原则**：
- 每个基准有明确的输入 → 输出 → 计时方法
- oCAS 侧使用 criterion，外部竞品使用 Python `timeit` 或 CLI 计时
- 记录完整环境信息（OS、CPU、Rust 版本、Python 版本、竞品版本）
- 结果以"输入 → oCAS 耗时 → 竞品耗时 → 比值"格式记录

---

## 2. 基准测试用例

### 2.1 多项式因式分解

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| F1 | x^100-1（ℤ 完全分解） | `factor(x^100 - 1)` | < 50 ms | Symbolica、FLINT |
| F2 | 多元 3 变量（Wang EEZ） | `factor(x^4*y^3 + x^3*y^4 + x^2*y^2)` | < 100 ms | Symbolica |
| F3 | 代数数域一元 | `factor(x^4 - 2, QQ(sqrt(2)))` | < 50 ms | Symbolica |
| F4 | x^200-1（ℤ，大规模） | `factor(x^200 - 1)` | 标记基准 | Symbolica、FLINT |

### 2.2 Gröbner 基

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| G1 | cyclic-5 ℤ₁₃ | 标准 cyclic-5 5 变量 | < 0.1 s | msolve（实测 3 ms @ 0.10.1，2026-08-06）、Symbolica |
| G2 | cyclic-6 ℤ₁₃ | 标准 cyclic-6 6 变量 | < 5 s | msolve（实测 4 ms @ 0.10.1，2026-08-06）、Symbolica (~1 s) |
| G3 | cyclic-7 ℤ₁₃ | 标准 cyclic-7 7 变量 | 标记基准 | msolve（实测 55 ms @ 0.10.1，2026-08-06）、Symbolica |
| G4 | katsura-6 ℤ₁₃ | 标准 katsura-6 | < 1 s | msolve（实测 3 ms @ 0.10.1，2026-08-06）、Symbolica |
| G5 | katsura-7 ℤ₁₃ | 标准 katsura-7 | 标记基准 | msolve（实测 7 ms @ 0.10.1，2026-08-06）、Symbolica |

### 2.3 符号积分

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| I1 | 有理函数 | `∫(x^3+2x+1)/(x^2+1) dx` | < 1 ms | Symbolica Rubi |
| I2 | 三角积分 | `∫sin(x)^3*cos(x)^2 dx` | < 1 ms | Symbolica Rubi |
| I3 | 特殊函数 | `∫exp(-x^2) dx` → erf | < 1 ms | Symbolica Rubi |
| I4 | 对数-代数混合 | `∫(log(x)/(x+1)^2) dx` | < 5 ms | Symbolica Rubi |
| I5 | Rubi 1892 题子集 | 从 Rubi 语料库采样 1892 题 | 总时间标记 | symbolica-integrate (111 s) |

### 2.4 多项式 GCD

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| C1 | gcd(x^50-1, x^30-1)（ℤ 大系数） | `gcd(x^50-1, x^30-1)` | < 1 s | Symbolica (~4 s)、FLINT |
| C2 | 多元 GCD | 三元多项式对 | < 100 ms | Symbolica |
| C3 | 模 GCD 性能 | 大系数单变量 | < 500 ms | FLINT |

### 2.5 JIT 求值

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| J1 | 单输出 JIT | `x^4 + 2x^3 + x^2 + 1`，batch 4096 | > 10× 解释器 | Symbolica SymJIT |
| J2 | 多输出 JIT | 3 输出表达式，batch 4096 | > 15× 解释器 | Symbolica SymJIT |
| J3 | SIMD 向量化 | 多项式求值，batch 4096 | > 8× 标量 | Symbolica SymJIT |

### 2.6 ODE 求解

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| O1 | 一阶可分离 | `y' = x*y` | < 10 ms | SymPy dsolve |
| O2 | 二阶常系数 | `y'' + 2y' + y = sin(x)` | < 50 ms | SymPy dsolve |
| O3 | Frobenius 级数 | 正则奇点 ODE | < 100 ms | SymPy dsolve |
| O4 | Laplace IVP | `y'' + y = 0, y(0)=1, y'(0)=0` | < 20 ms | SymPy dsolve |

### 2.7 数论

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| N1 | 30 位半素数分解 | `factor(9999999967 * 9999999943)` | < 10 s | SymPy factorint |
| N2 | BPSW 素性 | 2^61-1 (Mersenne 素数) | < 1 ms | SymPy isprime |
| N3 | 离散对数 | `dlog(3, 5, 101)` | < 10 ms | SymPy |
| N4 | 大系数模 GCD | 1000 度多项式 ℤ_p | < 500 ms | FLINT |

### 2.8 矩阵运算

| ID | 测试用例 | 输入 | 目标 | 对标竞品 |
|---|---|---|---|---|
| M1 | 20×20 整数矩阵 rref | 随机整数矩阵 | < 100 ms | SymPy DomainMatrix |
| M2 | 30×30 整数矩阵行列式 | 随机整数矩阵 | < 500 ms | SymPy DomainMatrix |
| M3 | 线性系统求解 | 20×20 系统 | < 100 ms | SymPy DomainMatrix |

---

## 3. 环境记录格式

每个基准结果必须附带：

```
环境：
- OS: Windows 11 Pro (26200)
- CPU: Intel Core Ultra 7 255H
- Rust: 1.97.1 (toolchain)
- Python: 3.12
- oCAS: 0.26.0
- Symbolica: 2.2.0 (pip show symbolica)
- SymPy: 1.14.0 (pip show sympy)
- FLINT: 3.6.0 (if applicable)
- msolve: 0.10.1 (if applicable)
```

---

## 4. 实现计划

### 4.1 oCAS 侧

- 扩展 `ocas-tests/benches/` 下的 criterion 基准
- 每个 ID 对应一个 criterion benchmark group
- 结果输出为 JSON（criterion `--output-format`）

### 4.2 竞品侧

- 扩展 `ocas-tests/scripts/compare_sympy.py` 框架
- 新增 `compare_multi.py`：支持 Symbolica、SymPy、SageMath（可选）
- 每个 ID 对应一个 Python 函数，使用 `timeit` 计时
- 结果输出为 JSON

### 4.3 结果汇总

- 生成 Markdown 表格：输入 → oCAS → 竞品 → 比值
- 归档至 `docs/planning/BENCHMARK_RESULTS_CN.md`
- 每次版本发布时更新

---

## 5. 优先执行顺序

1. **G2（cyclic-6 ℤ₁₃）**：直接对标 msolve，量化多模策略收益
2. **F1（x^100-1）**：因式分解经典基准
3. **I5（Rubi 子集）**：量化积分广度差距
4. **C1（大系数 GCD）**：模 GCD 性能
5. **J1/J2（JIT 求值）**：JIT 性能对标
6. 其余按需执行
