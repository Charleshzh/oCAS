# Benchmarks & Comparison

oCAS ships a [criterion](https://bheisler.github.io/criterion.rs/)-based
benchmark suite in `ocas-tests/benches/`, plus cross-language comparison
harnesses for SymPy, SageMath, and Symbolica. This chapter explains how to run
each and what they measure.

---

## Running oCAS benchmarks

```bash
# All benchmarks
cargo bench --workspace

# A specific benchmark
cargo bench --bench poly_gcd
cargo bench --bench poly_factor
cargo bench --bench groebner
cargo bench --bench poly_multivariate_gcd
cargo bench --bench roots

# Faster, less precise runs
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

| Benchmark | Covers |
|---|---|
| `arena` | Arena allocation throughput |
| `parse` / `normalize` | Expression parsing and normalization |
| `poly_dense` / `poly_sparse` | Polynomial arithmetic |
| `poly_gcd` | Univariate polynomial GCD |
| `poly_multivariate_gcd` | Multivariate polynomial GCD |
| `poly_factor` | Square-free factorization |
| `hensel_factor` | Hensel-lifting full factorization |
| `roots` | Real root isolation |
| `groebner` | Gröbner bases (cyclic-n ideals) |
| `calculus` / `rewrite` | Differentiation, Taylor, rule-based simplification |
| `eval_interpreter` / `eval_jit` / `eval_simd` | Numeric evaluation paths |
| `sympy_comparison` | Head-to-head vs SymPy |

---

## SymPy comparison (automated)

The `sympy_comparison` benchmark drives SymPy through a `uv`-managed Python
subprocess (`scripts/compare_sympy.py`) and feeds the elapsed nanoseconds into
criterion via `iter_custom`, so oCAS and SymPy appear side-by-side in the same
report.

```bash
# Requires `uv` on PATH; the Python env is provisioned automatically.
cargo bench --bench sympy_comparison
```

Supported tasks: `parse`, `diff`, `expand`, `factor`, `gcd`, `series`.

---

## SageMath comparison (local, manual)

SageMath is too heavy to install in CI, so the comparison runs locally via the
`sage` interpreter. The harness mirrors `bench_sympy.py`'s output contract.

```bash
# From the ocas-tests directory (requires SageMath installed)
sage scripts/bench_sage.py factor "x^30 - 1" 100
```

Tasks: `parse`, `diff`, `expand`, `factor`, `gcd`. Note that SageMath uses `^`
for exponentiation (same as oCAS), so no syntax translation is needed.

---

## JIT & evaluation

The `eval_jit` benchmark compares the Cranelift JIT against the stack-based
interpreter on single- and multi-output workloads (1000 calls each, criterion
`iter_custom` timing):

| Workload | Interpreter | JIT | Speedup |
|---|---|---|---|
| Polynomial (single output) | 221 µs | 2.27 µs | **97×** |
| Trig 3-output | 479 µs | 22.4 µs | **21×** |

```bash
cargo bench --bench eval_jit --features jit
```

The multi-output case compiles three expressions (`sin(x)`, `cos(x)`,
`sin(x)/cos(x)`) into one evaluator via `compile_multi`, sharing the `sin(x)`
subexpression across outputs; `call_into` writes results into a stack-allocated
buffer so each call performs zero heap allocation.

### Streaming evaluation

`StreamingEvaluator` reuses internal buffers across rows, so processing a
million-row dataset uses constant memory regardless of stream length:

| Workload | Per-row `evaluate` | `StreamingEvaluator` | Speedup |
|---|---|---|---|
| 100k rows, polynomial | 23.4 ms | 16.8 ms | **28%** |

```bash
cargo bench --bench eval_streaming
```

### f32 mixed precision

`compile_jit_f32` / `compile_vector_evaluator_f32` generate single-precision
code. On the same hardware the SIMD evaluator doubles its lane count (16 lanes
vs 8 for f64 on AVX-512). Use when f32 accuracy is sufficient.

---

## Symbolica comparison (local, manual)

[Symbolica](https://github.com/symbolica-dev/symbolica) is the primary
performance reference for oCAS. Because Symbolica uses an AGPL-style license
and ships as a separate Cargo workspace, it is **not** linked into the oCAS
build. Instead, run Symbolica's own example binaries from the source checkout
and compare timings manually.

Recommended comparison matrix (Symbolica examples → oCAS benchmarks):

| Symbolica example | oCAS benchmark | Workload |
|---|---|---|
| `polynomial_gcd` | `poly_gcd` | Integer & rational polynomial GCD |
| `factorization` | `poly_factor` | `x^n - 1` square-free / full factorization |
| `groebner_basis` | `groebner` | cyclic-4 ideal |
| `derivative` | `calculus` | Symbolic differentiation |
| `series` | `calculus` | Taylor expansion |

```bash
# Run a Symbolica example (from the symbolica source root)
cd ../symbolica
cargo run --release --example polynomial_gcd

# Then run the corresponding oCAS benchmark
cd ../ocas
cargo bench --bench poly_gcd -- --warm-up-time 0.5 --measurement-time 1
```

---

## Correctness comparison

Beyond performance, oCAS ships a correctness cross-validation framework in
`ocas-tests/tests/correctness/`. It runs 82 automated tests across 16
mathematical modules, comparing oCAS results against SymPy, SageMath, and
Symbolica. See the [Correctness](./correctness.md) chapter for details.

---

## Reporting results

criterion writes HTML reports to `target/criterion/`. Open
`target/criterion/poly_gcd/index.html` in a browser to inspect distributions,
regressions, and side-by-side comparisons.
