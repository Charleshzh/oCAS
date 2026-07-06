# 正确性

oCAS 包含自动化正确性交叉验证框架，将结果与三个参考系统进行对比：SymPy、SageMath 和 Symbolica。
本章介绍该框架、如何运行及其当前已知限制。

---

## 框架概述

正确性套件位于 `ocas-tests/tests/correctness/`，包含 **82 项测试，覆盖 16 个数学模块**。每项测试：

1. 生成输入（表达式、多项式、方程组等）
2. 使用 oCAS 计算结果
3. 使用参考系统计算等价结果
4. 断言两者语义相等

模块覆盖 oCAS 的全部功能：

| 模块 | 测试数 | 覆盖 |
|---|---|---|
| `algebra` | 8 | 化简、恒等律 |
| `calculus_diff` | 7 | 符号微分 |
| `calculus_int` | 5 | 启发式积分 |
| `calculus_series` | 4 | Taylor 级数展开 |
| `domain_integer` | 5 | 整数运算、GCD |
| `domain_rational` | 4 | 有理数运算 |
| `domain_finite_field` | 3 | 有限域运算 |
| `evaluation` | 6 | 数值求值 |
| `groebner` | 4 | Gröbner 基计算 |
| `linear_algebra` | 6 | 矩阵运算、线性求解 |
| `parsing` | 5 | 表达式解析与输出 |
| `poly_dense` | 5 | 稠密多项式运算 |
| `poly_factor` | 4 | 无平方因子与完全因式分解 |
| `poly_gcd` | 5 | 多项式 GCD |
| `poly_sparse` | 5 | 稀疏多元多项式运算 |
| `solvers` | 6 | 线性与丢番图求解器 |

---

## 难度分级

测试按难度分类，便于定位调试目标：

| 级别 | 描述 | 数量 |
|---|---|---|
| Trivial | 基本健全性检查（如 `x + 0 = x`） | ~20 |
| Easy | 单步操作（如 `d/dx x^3`） | ~30 |
| Medium | 多步或中等复杂度 | ~20 |
| Hard | 大型表达式、边界情况 | ~8 |
| Extreme | 已知会触及限制 | ~4 |

Extreme 级测试**预期失败**，用于记录已知差距（如 Wilkinson 多项式求根：10 个实根中仅找到 8 个）。

---

## 运行测试

```bash
# 运行全部正确性测试
cargo test -p ocas-tests --test correctness

# 运行特定模块
cargo test -p ocas-tests --test correctness algebra

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

运行完整套件可生成 `correctness_report.md`，包含每个模块的通过/失败统计及已知限制注释。

```bash
cargo test -p ocas-tests --test correctness -- --generate-report
```

报告包含：
- 每个模块的通过/失败计数
- 失败测试列表及预期与实际结果对比
- 难度级别分解
- 已知差距注释及跟踪 issue 链接

---

## 已知限制

| 问题 | 模块 | 状态 |
|---|---|---|
| Wilkinson n=10：10 个实根中仅找到 8 个 | `poly_sparse` / roots | 调查中 |
| `sin(x)^2 + cos(x)^2 → 1` 需要 `egg` feature | `algebra` | 启用 `egg` feature 后正常 |
| 部分启发式积分尚未实现 | `calculus_int` | 0.14（Risch）中扩展 |
| 多项式因式分解限于 Z[x] | `poly_factor` | Q[x] 尚未支持 |

---

## 参见

- [基准与性能对比](./performance.md) — 基准套件详情
- [Rust API](./rust-api.md) — 测试中使用的核心类型
- [贡献](./contributing.md) — 如何添加新的正确性测试
