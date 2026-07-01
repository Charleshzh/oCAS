# Benchmarks & Comparison / 基准与性能对比

**English**

oCAS ships a [criterion](https://bheisler.github.io/criterion.rs/)-based
benchmark suite in `ocas-tests/benches/`, plus cross-language comparison
harnesses for SymPy, SageMath, and Symbolica. This chapter explains how to run
each and what they measure.

**中文**

oCAS 在 `ocas-tests/benches/` 中提供基于 criterion 的基准套件，以及面向 SymPy、SageMath、Symbolica 的跨语言对比工具。本章说明如何运行各项基准及其测量内容。

---

## Running oCAS benchmarks / 运行 oCAS 基准

```bash
# All benchmarks
cargo bench --workspace

# A specific benchmark
cargo bench --bench poly_gcd
cargo bench --bench poly_factor
cargo bench --bench groebner

# Faster, less precise runs
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

**English**

| Benchmark | Covers |
|---|---|
| `arena` | Arena allocation throughput |
| `parse` / `normalize` | Expression parsing and normalization |
| `poly_dense` / `poly_sparse` | Polynomial arithmetic |
| `poly_gcd` | Univariate polynomial GCD |
| `poly_factor` | Square-free factorization |
| `groebner` | Gröbner bases (cyclic-n ideals) |
| `calculus` / `rewrite` | Differentiation, Taylor, rule-based simplification |
| `eval_interpreter` / `eval_jit` / `eval_simd` | Numeric evaluation paths |
| `sympy_comparison` | Head-to-head vs SymPy |

**中文**

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

## SymPy comparison (automated) / SymPy 对比（自动化）

**English**

The `sympy_comparison` benchmark drives SymPy through a `uv`-managed Python
subprocess (`scripts/bench_sympy.py`) and feeds the elapsed nanoseconds into
criterion via `iter_custom`, so oCAS and SymPy appear side-by-side in the same
report.

**中文**

`sympy_comparison` 基准通过 `uv` 管理的 Python 子进程（`scripts/bench_sympy.py`）驱动 SymPy，用 `iter_custom` 将耗时（纳秒）接入 criterion，使 oCAS 与 SymPy 在同一报告中并列展示。

```bash
# Requires `uv` on PATH; the Python env is provisioned automatically.
cargo bench --bench sympy_comparison
```

Supported tasks: `parse`, `diff`, `expand`, `factor`, `gcd`, `series`.

---

## SageMath comparison (local, manual) / SageMath 对比（本地手动）

**English**

SageMath is too heavy to install in CI, so the comparison runs locally via the
`sage` interpreter. The harness mirrors `bench_sympy.py`'s output contract.

**中文**

SageMath 安装体积过大，不适合进 CI，因此对比在本地通过 `sage` 解释器运行。该工具与 `bench_sympy.py` 的输出契约一致。

```bash
# From the ocas-tests directory (requires SageMath installed)
sage scripts/bench_sage.py factor "x^30 - 1" 100
```

Tasks: `parse`, `diff`, `expand`, `factor`, `gcd`. Note that SageMath uses `^`
for exponentiation (same as oCAS), so no syntax translation is needed.

---

## Symbolica comparison (local, manual) / Symbolica 对比（本地手动）

**English**

[Symbolica](https://github.com/symbolica-dev/symbolica) is the primary
performance reference for oCAS. Because Symbolica uses an AGPL-style license
and ships as a separate Cargo workspace, it is **not** linked into the oCAS
build. Instead, run Symbolica's own example binaries from the source checkout
and compare timings manually.

**中文**

[Symbolica](https://github.com/symbolica-dev/symbolica) 是 oCAS 的主要性能参考。由于 Symbolica 采用 AGPL 类许可证且作为独立 Cargo workspace 发布，**不会**链接进 oCAS 构建。请从源码检出运行 Symbolica 自带的示例二进制，手动对比耗时。

**English**

Recommended comparison matrix (Symbolica examples → oCAS benchmarks):

| Symbolica example | oCAS benchmark | Workload |
|---|---|---|
| `polynomial_gcd` | `poly_gcd` | Integer & rational polynomial GCD |
| `factorization` | `poly_factor` | `x^n - 1` square-free / full factorization |
| `groebner_basis` | `groebner` | cyclic-4 ideal |
| `derivative` | `calculus` | Symbolic differentiation |
| `series` | `calculus` | Taylor expansion |

**中文**

推荐对比矩阵（Symbolica 示例 → oCAS 基准）：

| Symbolica 示例 | oCAS 基准 | 工作负载 |
|---|---|---|
| `polynomial_gcd` | `poly_gcd` | 整数与有理多项式 GCD |
| `factorization` | `poly_factor` | `x^n - 1` 无平方/完整因式分解 |
| `groebner_basis` | `groebner` | cyclic-4 理想 |
| `derivative` | `calculus` | 符号微分 |
| `series` | `calculus` | Taylor 展开 |

```bash
# Run a Symbolica example (from the symbolica source root)
cd ../symbolica
cargo run --release --example polynomial_gcd

# Then run the corresponding oCAS benchmark
cd ../ocas
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

---

## Reporting results / 报告结果

**English**

criterion writes HTML reports to `target/criterion/`. Open
`target/criterion/poly_gcd/index.html` in a browser to inspect distributions,
regressions, and side-by-side comparisons.

**中文**

criterion 将 HTML 报告写入 `target/criterion/`。在浏览器打开
`target/criterion/poly_gcd/index.html` 可查看分布、回归与并列对比。
