# 基准与性能对比

oCAS 在 `ocas-tests/benches/` 中提供基于 criterion 的基准套件，以及面向 SymPy、SageMath、Symbolica 的跨语言对比工具。本章说明如何运行各项基准及其测量内容。

---

## 运行 oCAS 基准

```bash
# 所有基准
cargo bench --workspace

# 单项基准
cargo bench --bench poly_gcd
cargo bench --bench poly_factor
cargo bench --bench groebner

# 更快但精度较低的运行
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

| 基准 | 覆盖 |
|---|---|
| `arena` | Arena 分配吞吐 |
| `parse` / `normalize` | 表达式解析与规范化 |
| `poly_dense` / `poly_sparse` | 多项式运算 |
| `poly_gcd` | 一元多项式 GCD |
| `poly_factor` | 无平方因子分解 |
| `groebner` | Gröbner 基（cyclic-n 理想） |
| `calculus` / `rewrite` | 微分、Taylor、规则化简 |
| `eval_interpreter` / `eval_jit` / `eval_simd` | 数值求值路径 |
| `sympy_comparison` | 与 SymPy 直接对比 |

---

## SymPy 对比（自动化）

`sympy_comparison` 基准通过 `uv` 管理的 Python 子进程（`scripts/bench_sympy.py`）驱动 SymPy，用 `iter_custom` 将耗时（纳秒）接入 criterion，使 oCAS 与 SymPy 在同一报告中并列展示。

```bash
# 需要 PATH 中有 `uv`；Python 环境会自动准备
cargo bench --bench sympy_comparison
```

支持的任务：`parse`、`diff`、`expand`、`factor`、`gcd`、`series`。

---

## SageMath 对比（本地手动）

SageMath 安装体积过大，不适合进 CI，因此对比在本地通过 `sage` 解释器运行。该工具与 `bench_sympy.py` 的输出契约一致。

```bash
# 在 ocas-tests 目录下运行（需已安装 SageMath）
sage scripts/bench_sage.py factor "x^30 - 1" 100
```

任务：`parse`、`diff`、`expand`、`factor`、`gcd`。注意 SageMath 用 `^`
表示乘方（与 oCAS 相同），无需语法转换。

---

## Symbolica 对比（本地手动）

[Symbolica](https://github.com/symbolica-dev/symbolica) 是 oCAS 的主要性能参考。由于 Symbolica 采用 AGPL 类许可证且作为独立 Cargo workspace 发布，**不会**链接进 oCAS 构建。请从源码检出运行 Symbolica 自带的示例二进制，手动对比耗时。

推荐对比矩阵（Symbolica 示例 → oCAS 基准）：

| Symbolica 示例 | oCAS 基准 | 工作负载 |
|---|---|---|
| `polynomial_gcd` | `poly_gcd` | 整数与有理多项式 GCD |
| `factorization` | `poly_factor` | `x^n - 1` 无平方/完整因式分解 |
| `groebner_basis` | `groebner` | cyclic-4 理想 |
| `derivative` | `calculus` | 符号微分 |
| `series` | `calculus` | Taylor 展开 |

```bash
# 运行 Symbolica 示例（在 symbolica 源码根目录）
cd ../symbolica
cargo run --release --example polynomial_gcd

# 再运行对应的 oCAS 基准
cd ../ocas
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

---

## 报告结果

criterion 将 HTML 报告写入 `target/criterion/`。在浏览器打开
`target/criterion/poly_gcd/index.html` 可查看分布、回归与并列对比。
