# 正确性

oCAS 包含自动化正确性交叉验证框架，将结果与三个参考系统进行对比：SymPy、SageMath 和 Symbolica。
本章介绍该框架、如何运行及其当前已知限制。

---

## 框架概述

正确性套件位于 `ocas-tests/tests/correctness/`，包含 **201 项 `#[test]` 测试，覆盖 19 个数学模块**（截至 0.24.0；其中 35 项标记 `#[ignore]` 记录已知差距）。每项测试：

1. 生成输入（表达式、多项式、方程组等）
2. 使用 oCAS 计算结果
3. 使用参考系统计算等价结果
4. 断言两者语义相等

模块覆盖 oCAS 的全部功能：

| 模块 | 测试数 | 覆盖 |
|---|---|---|
| `calculus` | 16 | 微分、Taylor、积分（SymPy 对比） |
| `evaluation` | 6 | 数值求值 |
| `finite_field` | 5 | 有限域运算 |
| `groebner` | 18 | Gröbner 基计算 |
| `integral_risch` | 15 | Risch 符号积分 |
| `linear_solve` | 5 | 线性求解器 |
| `matrix` | 5 | 矩阵运算 |
| `normalize` | 8 | 表达式规范化 |
| `ntheory` | 11 | 数论（与 SymPy `ntheory` 交叉验证） |
| `ode` | 34 | ODE 求解 |
| `parse` | 6 | 表达式解析与输出 |
| `partial_fraction` | 8 | 部分分式分解 |
| `poly_arithmetic` | 6 | 稠密/稀疏多项式运算 |
| `poly_factor` | 12 | 无平方与完全因式分解 |
| `poly_factor_anf` | 21 | 代数数域因式分解 |
| `poly_gcd` | 5 | 多项式 GCD |
| `resultant` | 8 | 结式 |
| `rewrite` | 8 | 重写与化简 |
| `root_isolation` | 4 | 实根隔离 |

---

## 忽略的测试（已知差距）

部分测试标记 `#[ignore]`（当前 36 项），用于记录已知差距——它们预期失败，
仅在需要复现/推进时手动运行（`cargo test -p ocas-tests --test correctness
-- --ignored`）。例如 Wilkinson 多项式求根：10 个实根中仅找到 8 个。

---

## 运行测试

```bash
# 运行全部正确性测试
cargo test -p ocas-tests --test correctness

# 运行特定模块（如 ODE 求解）
cargo test -p ocas-tests --test correctness ode

# 运行被忽略的测试（已知差距）
cargo test -p ocas-tests --test correctness -- --ignored

# 详细输出以检查失败
cargo test -p ocas-tests --test correctness -- --nocapture
```

测试无需外部依赖 —— 所有参考计算通过 `uv` 管理的 Python 子进程使用 SymPy，子进程自动自举。

---

## 对比工具

独立脚本提供针对 SageMath 和 Symbolica 的手动交叉检查，用于深入调查：

```bash
# SageMath（需本地安装 `sage`）
cd ocas-tests
sage scripts/bench_sage.py factor "x^30 - 1" 100

# Symbolica（需 Symbolica 源码检出）
cd ../symbolica
cargo run --release --example factorization
```

这些工具在开发期间用于验证自动化测试套件，并维护用于手动回归测试。

---

## 审计报告

`ocas-tests/scripts/generate_audit_report.py` 会运行全部测试（含 `--ignored`
部分）并生成审计报告：

```bash
cd ocas-tests
python scripts/generate_audit_report.py
```

报告写入 `docs/planning/correctness/audit-<日期>.md`，包含：
- 普通测试与被忽略测试的通过/失败/忽略计数摘要
- 失败测试列表
- 与 Symbolica 的因式分解计时对比

---

## 已知限制

| 问题 | 模块 | 状态 |
|---|---|---|
| Wilkinson n=10：10 个实根中仅找到 8 个 | `root_isolation` | 已知差距（`#[ignore]`） |
| `sin(x)^2 + cos(x)^2 → 1` 需要 `egg` feature | `rewrite` | 启用 `egg` feature 后可化简 |
| Bernoulli forcing y^n 使线性系数提取混淆 | `ode` | 已知限制（`#[ignore]`） |
| 积分器缺 tan/sec 表项 | `ode` | 计划中（`#[ignore]`） |

---

## 参见

- [基准与性能对比](./performance.md) — 基准套件详情
- [Rust API](./api/rust.md) — 测试中使用的核心类型
- [贡献](./contributing.md) — 如何添加新的正确性测试
