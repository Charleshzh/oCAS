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
cargo bench --bench poly_multivariate_gcd
cargo bench --bench roots

# 更快但精度较低的运行
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

| 基准 | 覆盖 |
|---|---|
| `arena` | Arena 分配吞吐 |
| `parse` / `normalize` | 表达式解析与规范化 |
| `poly_dense` / `poly_sparse` | 多项式运算 |
| `poly_gcd` | 一元多项式 GCD |
| `poly_multivariate_gcd` | 多元多项式 GCD |
| `poly_factor` | 无平方因子分解 |
| `hensel_factor` | Hensel 提升完全因式分解 |
| `roots` | 实根隔离 |
| `groebner` | Gröbner 基（cyclic-n 理想） |
| `calculus` / `rewrite` | 微分、Taylor、规则化简 |
| `eval_interpreter` / `eval_jit` / `eval_simd` | 数值求值路径 |
| `sympy_comparison` | 与 SymPy 直接对比 |

---

## SymPy 对比（自动化）

`sympy_comparison` 基准通过 `uv` 管理的 Python 子进程（`scripts/compare_sympy.py`）驱动 SymPy，用 `iter_custom` 将耗时（纳秒）接入 criterion，使 oCAS 与 SymPy 在同一报告中并列展示。

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

## JIT 与求值

`eval_jit` 基准对比 Cranelift JIT 与栈式解释器在单输出与多输出工作负载上的表现（各 1000 次调用，criterion `iter_custom` 计时）：

| 工作负载 | 解释器 | JIT | 加速比 |
|---|---|---|---|
| 多项式（单输出） | 221 µs | 2.27 µs | **97×** |
| 三角函数 3 输出 | 479 µs | 22.4 µs | **21×** |

```bash
cargo bench --bench eval_jit --features jit
```

多输出用例通过 `compile_multi` 将三个表达式（`sin(x)`、`cos(x)`、`sin(x)/cos(x)`）编译为一个求值器，跨输出共享 `sin(x)` 子表达式；`call_into` 将结果写入栈分配缓冲，每次调用零堆分配。

### 流式求值

`StreamingEvaluator` 跨行复用内部缓冲，处理百万行数据集时内存恒定：

| 工作负载 | 逐行 `evaluate` | `StreamingEvaluator` | 加速比 |
|---|---|---|---|
| 100k 行，多项式 | 23.4 ms | 16.8 ms | **28%** |

```bash
cargo bench --bench eval_streaming
```

### f32 混合精度

`compile_jit_f32` / `compile_vector_evaluator_f32` 生成单精度代码。同一硬件上 SIMD 求值器通道数翻倍（AVX-512 上 f64 为 8、f32 为 16）。当 f32 精度足够时使用。

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

## 正确性对比

除性能外，oCAS 还在 `ocas-tests/tests/correctness/` 中提供了正确性交叉验证框架。
它运行 82 项自动化测试，覆盖 16 个数学模块，将 oCAS 结果与 SymPy、SageMath、Symbolica
进行对比。详见[正确性](./correctness.md)章节。

---

## 报告结果

criterion 将 HTML 报告写入 `target/criterion/`。在浏览器打开
`target/criterion/poly_gcd/index.html` 可查看分布、回归与并列对比。
