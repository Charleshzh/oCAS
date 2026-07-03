# Evaluation & JIT

oCAS provides three paths for numeric evaluation of symbolic expressions:
an interpreter (stack VM), a Cranelift JIT compiler, and a SIMD-vectorized
batch evaluator. This chapter explains each path and when to use it.

---

## Stack VM interpreter

The default evaluation path. `ExpressionEvaluator` compiles an `Atom`
expression tree into a sequence of stack-machine instructions, then
executes them on a flat operand stack.

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// Parse and compile
let e = parse(&ctx, "x^2 + y").unwrap();
let ev = ExpressionEvaluator::<f64>::compile(e).unwrap();

// Evaluate with parameter values
let result = ev.evaluate(&[3.0, 1.0]).unwrap();  // [10.0]
let result = ev.evaluate(&[2.0, 0.0]).unwrap();  // [4.0]
```

The compiler automatically detects free variables, assigns them parameter
slots in sorted order, and optimizes the instruction sequence (constant
folding, copy-chain removal).

Python:

```python
ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
print(ev.evaluate([3.0, 1.0]))  # [10.0]
print(ev.evaluate([2.0, 0.0]))  # [4.0]
```

---

## Cranelift JIT

With the `jit` feature, oCAS compiles expressions to native machine code
via Cranelift. This is ideal for repeatedly evaluating the same expression
with different inputs — amortizing the compilation cost over many calls.

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

    let mut engine = JitEngine::new();
    let compiled = engine.compile::<f64, _>(e).unwrap();
    let result = compiled.call(&[0.5, 1.0]);  // ~0.2590
}
```

The JIT path translates the same IR used by the interpreter into native
x86-64 or aarch64 code. For expressions evaluated thousands of times, this
can yield a 10–50× speedup over the interpreter.

---

## SIMD batch evaluation

The `simd` feature enables vectorized evaluation using `wide::f64x4`,
computing four inputs simultaneously with SIMD instructions.

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

SIMD evaluation processes batches of four parameter sets in parallel,
ideal for parameter sweeps, plotting, and Monte Carlo workloads.

---

## Comparison

| Path | Feature | Latency | Throughput | Best for |
|---|---|---|---|---|
| Interpreter | (default) | Low (no compile overhead) | Medium | Single-shot, interactive use |
| JIT | `jit` | High (compile) + Low (run) | High | Repeated evaluation of same expression |
| SIMD | `simd` | Low–Medium | Very High (4×) | Batch evaluation, parameter sweeps |

### When to use each

- **Interpreter**: Interactive REPL sessions, one-off evaluations, Python API.
- **JIT**: Repetitive numeric evaluation (root-finding iterations, optimization loops).
- **SIMD**: Batching many inputs (plotting, sensitivity analysis, ML pipelines).

---

## User-defined functions

`FunctionMap` allows registering custom Rust functions callable from
evaluated expressions.

```rust
let mut fns = FunctionMap::<f64>::new();
fns.insert("f", |args: &[f64]| args[0].sin() * args[1].cos());

let e = parse(&ctx, "f(x, y)").unwrap();
let ev = ExpressionEvaluator::<f64>::compile_with(e, fns).unwrap();
let result = ev.evaluate(&[0.5, 1.0]).unwrap();
```

---

## See also

- [Rust API](./rust-api.md) — building expressions for evaluation
- [Performance](./performance.md) — benchmark results across all three paths
