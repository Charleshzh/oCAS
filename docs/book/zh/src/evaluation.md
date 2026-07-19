# 求值与 JIT

oCAS 为符号表达式的数值求值提供三条路径：解释器（栈 VM）、Cranelift JIT 编译器和 SIMD 向量化批量求值器。
本章解释每条路径及何时使用。

---

## 栈 VM 解释器

默认求值路径。`ExpressionEvaluator` 将 `Atom` 表达式树编译为栈机指令序列，然后在扁平操作数栈上执行。

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// 解析并编译
let e = parse(&ctx, "x^2 + y").unwrap();
let ev = ExpressionEvaluator::<f64>::compile(e).unwrap();

// 用参数值求值
let result = ev.evaluate(&[3.0, 1.0]).unwrap();  // [10.0]
let result = ev.evaluate(&[2.0, 0.0]).unwrap();  // [4.0]
```

编译器自动检测自由变量，按排序顺序分配参数槽，并优化指令序列（常量折叠、拷贝链消除）。

Python：

```python
ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
print(ev.evaluate([3.0, 1.0]))  # [10.0]
print(ev.evaluate([2.0, 0.0]))  # [4.0]
```

---

## Cranelift JIT

启用 `jit` feature 后，oCAS 通过 Cranelift 将表达式编译为本机机器码。
这适用于使用不同输入重复求值同一表达式 —— 将编译开销分摊到多次调用。

```bash
cargo build -p ocas --features jit
```

```rust
#[cfg(feature = "jit")]
{
    use ocas::prelude::*;

    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let e = parse(&ctx, "sin(x) * cos(y)").unwrap();

    let ev: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(e).unwrap();
    let compiled = ev.compile_jit().unwrap();
    let result = compiled.call(&[0.5, 1.0]);  // ~0.2590
}
```

JIT 路径将解释器使用的同一 IR 翻译为本机 x86-64 或 aarch64 代码。对于需要数千次求值的表达式，
最高可比解释器快 **97 倍**（见[基准](./performance.md#jit--evaluation)）。

### 多输出 JIT

`compile_multi` 将多个表达式编译为一个求值器，跨输出共享公共子表达式；`call_into` 将结果写入
调用方提供的缓冲（每次调用零分配）。

```rust
#[cfg(feature = "jit")]
{
    use ocas::prelude::*;

    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let e1 = parse(&ctx, "sin(x) + 1").unwrap();
    let e2 = parse(&ctx, "sin(x) * 2").unwrap();

    let ev: ExpressionEvaluator<f64> =
        ExpressionEvaluator::compile_multi(&[e1, e2]).unwrap();
    let compiled = ev.compile_jit().unwrap();
    let mut out = [0.0f64; 2];
    compiled.call_into(&[1.0], &mut out);
}
```

### f32 混合精度

`compile_jit_f32` 生成单精度代码（libm `*f` 符号）；`compile_vector_evaluator_f32` 在同一硬件上将
SIMD 通道数翻倍。当 f32 精度足够时使用。

---

## SIMD 批量求值

`simd` feature 使用 `pulp` 启用向量化求值，通过运行时检测的 SIMD 宽度（SSE2/AVX2/AVX-512）同时
计算多个输入。

```bash
cargo build -p ocas --features simd
```

```rust
#[cfg(feature = "simd")]
{
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let e = parse(&ctx, "x^2 + y").unwrap();

    let ev = VectorEvaluator::<f64>::compile(e).unwrap();
    let results = ev.evaluate_batch(&[
        vec![1.0, 0.0],   // 1^2 + 0 = 1
        vec![2.0, 1.0],   // 2^2 + 1 = 5
        vec![3.0, 2.0],   // 3^2 + 2 = 11
        vec![0.0, 5.0],   // 0^2 + 5 = 5
    ]).unwrap();
    // results = [1.0, 5.0, 11.0, 5.0]
}
```

SIMD 求值并行处理四组参数，适用于参数扫描、绘图和 Monte Carlo 工作负载。

---

## 对比

| 路径 | Feature | 延迟 | 吞吐 | 最适合 |
|---|---|---|---|---|
| 解释器 | （默认） | 低（无编译开销） | 中 | 单次求值、交互使用 |
| JIT | `jit` | 高（编译）+ 低（运行） | 高 | 同一表达式的重复求值 |
| SIMD | `simd` | 低–中 | 极高（4×） | 批量求值、参数扫描 |

### 使用建议

- **解释器**：交互式 REPL 会话、一次性求值、Python API。
- **JIT**：重复数值求值（求根迭代、优化循环）。
- **SIMD**：大量输入的批量处理（绘图、灵敏度分析、ML 管道）。

---

## 用户自定义函数

`FunctionMap` 允许注册自定义 Rust 函数，使其可在求值表达式中调用。

```rust
let mut fns = FunctionMap::<f64>::new();
fns.insert("f", |args: &[f64]| args[0].sin() * args[1].cos());

let e = parse(&ctx, "f(x, y)").unwrap();
let ev = ExpressionEvaluator::<f64>::compile_with(e, fns).unwrap();
let result = ev.evaluate(&[0.5, 1.0]).unwrap();
```

---

## 参见

- [Rust API](./rust-api.md) — 构建用于求值的表达式
- [基准与性能对比](./performance.md) — 三条路径的基准结果
