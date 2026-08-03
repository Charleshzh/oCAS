# Rust API Reference: Evaluation and JIT

oCAS's evaluation system compiles symbolic expressions into stack-machine instruction sequences and then executes them efficiently on numeric types. Three execution backends are supported: an interpreter (general purpose), JIT compilation (Cranelift generates native machine code), and SIMD batch evaluation (vectorized processing of many parameter sets).

**Module path**: `ocas_eval`

**Compilation pipeline**:

```text
Atom (arena-backed)
  → EvalTree (owned intermediate, optional constant-folding optimization)
  → Instr sequence (compilation + optimization passes)
  → ExpressionEvaluator (stack VM interpreter)
  → or JitCompiledFunction (Cranelift, feature = "jit")
  → or VectorEvaluator (SIMD batch, feature = "simd")
```

**Feature flags**:

| Flag | Purpose |
|---|---|
| `jit` | enables the Cranelift JIT compilation backend |
| `simd` | enables SIMD batch evaluation (`pulp`, runtime detection of SSE2/AVX2/AVX-512) |
| `fast-poly` | enables fast polynomial evaluation optimizations |

**Imports**:

```rust
use ocas_eval::{
    ExpressionEvaluator, EvaluationDomain, PowfExtension,
    FunctionMap, EvalTree,
    Instr, Instruction, Slot,
    EvaluationError,
    StreamingEvaluator,
    compile_atom, compile_atom_with,
    compile_atoms_multi, compile_atoms_multi_with,
};
use ocas_eval::instruction::BuiltinOp;

// JIT (feature = "jit")
use ocas_eval::jit::{JitCompiledFunction, JitCompiledF32, JitEngine, FloatWidth};

// SIMD (feature = "simd")
use ocas_eval::simd::{VectorEvaluator, VectorEvaluatorF32};
```

---

## EvalTree

**Signature**:

```rust
pub enum EvalTree {
    Num(f64),
    Var(String),
    Fun(String, Vec<EvalTree>),
    Add(Vec<EvalTree>),
    Mul(Vec<EvalTree>),
    Pow(Box<EvalTree>, Box<EvalTree>),
}
```

**Description**: an arena-independent, owned intermediate representation. Converted from `Atom`, it is decoupled from the arena lifetime, so the compilation pipeline can run multiple optimization passes freely.

**Variants**:

| Variant | Description |
|---|---|
| `Num(f64)` | numeric constant |
| `Var(String)` | variable reference |
| `Fun(String, Vec<EvalTree>)` | named function application |
| `Add(Vec<EvalTree>)` | summation |
| `Mul(Vec<EvalTree>)` | product |
| `Pow(Box<EvalTree>, Box<EvalTree>)` | power base^exp |

#### EvalTree::from_atom

**Signature**: `pub fn from_atom(atom: Atom<'_>) -> Self`

**Description**: converts an `Atom` into an owned `EvalTree`. Recursively traverses the expression tree, turning arena handles into independent data.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom<'_>` | the source expression |

**Returns**: an owned `EvalTree`.

**Example**:

```rust
use ocas_atom::AtomArena;
use ocas_eval::EvalTree;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let atom = ctx.add(&[ctx.var("x"), ctx.num(2)]);
let tree = EvalTree::from_atom(atom);
// tree = Add([Var("x"), Num(2.0)])
```

#### EvalTree::fold_constants

**Signature**: `pub fn fold_constants(&self) -> EvalTree`

**Description**: folds constant subtrees and applies algebraic identities. Runs before compilation to reduce the instruction count.

**Optimization rules**:

| Pattern | Rule |
|---|---|
| `Add` | remove zero terms; sum fully constant terms; fold single terms |
| `Mul` | absorb zero factors; remove unit factors; multiply fully constant factors; fold single factors |
| `Pow` | `x^1 → x`, `x^0 → 1`, evaluate `Num^Num` directly |
| `Fun` | evaluate built-in functions with fully constant arguments directly (external functions are not folded; they may have side effects) |

**Returns**: an optimized new `EvalTree`.

**Example**:

```rust
use ocas_eval::EvalTree;

let tree = EvalTree::Add(vec![
    EvalTree::Num(3.0),
    EvalTree::Mul(vec![
        EvalTree::Num(0.0),
        EvalTree::Var("x".into()),
    ]),
]);
let folded = tree.fold_constants();
assert_eq!(folded, EvalTree::Num(3.0));
// Output: 3 + 0*x → 3
```

---

## ExpressionEvaluator\<T\>

**Signature**:

```rust
pub struct ExpressionEvaluator<T: EvaluationDomain> { /* private fields */ }
```

**Description**: a compiled expression that can be evaluated numerically many times. It holds the instruction sequence, a preallocated stack, the constant table, and an optional user-defined function registry. The type parameter `T` implements `EvaluationDomain`, usually `f64`.

### Compilation methods

#### ExpressionEvaluator::compile

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile(atom: Atom<'_>) -> Result<Self>
}
```

**Description**: compiles an `Atom` into an executable evaluator. Internally it first converts to an `EvalTree`, then folds constants and generates the instruction sequence.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom<'_>` | the source expression |

**Returns**: `Ok(ExpressionEvaluator<T>)` or a compilation error.

**Errors**: see [`EvaluationError`](#evaluationerror).

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.mul(&[ctx.var("x"), ctx.var("x")]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
assert_eq!(eval.param_count(), 1);
assert_eq!(eval.result_count(), 1);
```

#### ExpressionEvaluator::compile_with

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_with(atom: Atom<'_>, map: FunctionMap<T>) -> Result<Self>
}
```

**Description**: compiles an expression with a user-defined function registry. `ExternalFun` nodes appearing in the expression are resolved through `map`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom<'_>` | the source expression |
| `map` | `FunctionMap<T>` | the user-defined function registry (pass an empty map when no external functions are needed) |

**Returns**: `Ok(ExpressionEvaluator<T>)` or a compilation/resolution error.

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, FunctionMap};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("myfunc", &[ctx.var("x")]);

let mut fmap = FunctionMap::<f64>::new();
fmap.register("myfunc", 1, Box::new(|args| args[0] * args[0]));

let eval = ExpressionEvaluator::compile_with(expr, fmap).unwrap();
let result = eval.evaluate(&[3.0]).unwrap();
assert_eq!(result, vec![9.0]);
// Output: myfunc(3) = 3² = 9
```

#### ExpressionEvaluator::compile_multi

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_multi(atoms: &[Atom<'_>]) -> Result<Self>
}
```

**Description**: compiles multiple `Atom`s into a single multi-output evaluator. Constants, common subexpressions (CSE), and stack slots are shared across all outputs to avoid redundant computation.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atoms` | `&[Atom<'_>]` | the multiple expressions |

**Returns**: `Ok(ExpressionEvaluator<T>)` with `result_count() == atoms.len()`.

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let expr1 = ctx.mul(&[x, x]);       // x²
let expr2 = ctx.add(&[x, ctx.num(1)]); // x+1

let eval: ExpressionEvaluator<f64> =
    ExpressionEvaluator::compile_multi(&[expr1, expr2]).unwrap();
let results = eval.evaluate(&[3.0]).unwrap();
assert_eq!(results, vec![9.0, 4.0]);
// Output: x²=9, x+1=4
```

#### ExpressionEvaluator::compile_multi_with

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_multi_with(atoms: &[Atom<'_>], map: FunctionMap<T>) -> Result<Self>
}
```

**Description**: the `compile_multi` variant with a function map (the map is passed by value, not as an `Option`).

### Execution methods

#### ExpressionEvaluator::evaluate

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn evaluate(&self, params: &[T]) -> Result<Vec<T>>
}
```

**Description**: executes the compiled expression with the given parameter values. Allocates a fresh stack and result buffer on every call.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `params` | `&[T]` | the parameter values; length must equal `param_count()` |

**Returns**: `Ok(Vec<T>)` with length equal to `result_count()`.

**Errors**:

| Error | Condition |
|---|---|
| `WrongArity` | `params.len() != param_count()` |
| `DivisionByZero` | division by zero encountered during evaluation |
| `FunctionNotFound` | an external function is not registered |

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let result = eval.evaluate(&[3.0]).unwrap();
assert_eq!(result[0], 10.0);
// Output: 3² + 1 = 10
```

#### ExpressionEvaluator::evaluate_with_stack

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn evaluate_with_stack(
        &self,
        params: &[T],
        stack: &mut Vec<T>,
        results: &mut Vec<T>,
    ) -> Result<()>
}
```

**Description**: evaluates using caller-provided buffers, avoiding per-call heap allocation. Prefer this method in streaming and batch scenarios.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `params` | `&[T]` | the parameter values |
| `stack` | `&mut Vec<T>` | the stack buffer (automatically resized to `stack_size()`) |
| `results` | `&mut Vec<T>` | the result buffer (filled automatically) |

**Returns**: `Ok(())` or an evaluation error.

### Accessors

#### ExpressionEvaluator::param_count

**Signature**: `pub fn param_count(&self) -> usize`

**Description**: returns the expected number of parameters.

#### ExpressionEvaluator::result_count

**Signature**: `pub fn result_count(&self) -> usize`

**Description**: returns the number of outputs (1 for a single-output expression, possibly greater for one created by `compile_multi`).

#### ExpressionEvaluator::stack_size

**Signature**: `pub fn stack_size(&self) -> usize`

**Description**: returns the stack size needed for evaluation (parameters + constants + temporaries + outputs). The `stack` argument of `evaluate_with_stack` should be preallocated with this capacity.

---

## JIT compilation (feature = "jit")

The JIT backend uses Cranelift to compile the instruction sequence into native machine code. Constants are embedded as immediates, and libm is used for transcendental functions. Performance is typically more than 10× faster than the interpreter.

### ExpressionEvaluator::compile_jit

**Signature**:

```rust
#[cfg(feature = "jit")]
impl ExpressionEvaluator<f64> {
    pub fn compile_jit(&self) -> Result<JitCompiledFunction>
}
```

**Description**: compiles the current evaluator's instruction sequence into double-precision (f64) native code.

**Returns**: `Ok(JitCompiledFunction)` or a compilation error.

**Errors**:

| Error | Condition |
|---|---|
| `UnsupportedOperation` | contains instructions the JIT cannot lower (e.g., external functions, `sec`/`csc`/`cot`) |
| `JitCompilationError` | the Cranelift backend failed to compile |

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let jit = eval.compile_jit().unwrap();
let result = jit.call(&[3.0]);
assert_eq!(result, vec![10.0]);
// Output: 3² + 1 = 10 (executed as native machine code)
```

### ExpressionEvaluator::compile_jit_f32

**Signature**:

```rust
#[cfg(feature = "jit")]
impl ExpressionEvaluator<f64> {
    pub fn compile_jit_f32(&self) -> Result<JitCompiledF32>
}
```

**Description**: compiles the instruction sequence into single-precision (f32) native code. Constants are narrowed from f64 to f32; results have f32 precision.

**Errors**: same as `compile_jit`.

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.mul(&[ctx.var("x"), ctx.num(2)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let jit_f32 = eval.compile_jit_f32().unwrap();
let result = jit_f32.call(&[3.0_f32]);
assert_eq!(result, vec![6.0_f32]);
// Output: 3 × 2 = 6 (f32 precision)
```

### JitCompiledFunction

**Signature**:

```rust
pub struct JitCompiledFunction { /* private fields */ }
```

**Description**: a JIT-compiled double-precision function that can be called repeatedly.

#### JitCompiledFunction::call

**Signature**: `pub fn call(&self, params: &[f64]) -> Vec<f64>`

**Description**: calls the JIT-compiled function with the given parameters. Allocates a new result vector on every call.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `params` | `&[f64]` | the parameter values |

**Returns**: a result vector of length equal to `result_count()`.

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let jit = eval.compile_jit().unwrap();
let result = jit.call(&[3.0]);
assert_eq!(result, vec![10.0]);
// Output: value of x²+1 at x=3 → 10
```

#### JitCompiledFunction::call_into

**Signature**: `pub fn call_into(&self, params: &[f64], results: &mut [f64])`

**Description**: calls the JIT function, writing results into a caller-provided buffer to avoid per-call heap allocation. `results.len()` must be at least `result_count()`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `params` | `&[f64]` | the parameter values |
| `results` | `&mut [f64]` | the result buffer |

**Example**:

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let jit = eval.compile_jit().unwrap();
let mut buf = vec![0.0_f64; jit.result_count()];
jit.call_into(&[42.0], &mut buf);
assert_eq!(buf, vec![1765.0]);
// Results are in buf (42² + 1 = 1765), no heap allocation
```

#### JitCompiledFunction::param_count

**Signature**: `pub fn param_count(&self) -> usize`

**Description**: returns the expected number of parameters.

#### JitCompiledFunction::result_count

**Signature**: `pub fn result_count(&self) -> usize`

**Description**: returns the number of results.

### JitCompiledF32

**Signature**:

```rust
pub struct JitCompiledF32 { /* private fields */ }
```

**Description**: a JIT-compiled single-precision function. The API is fully symmetric with `JitCompiledFunction`; parameters and return values are `f32`.

#### JitCompiledF32::call

**Signature**: `pub fn call(&self, params: &[f32]) -> Vec<f32>`

#### JitCompiledF32::call_into

**Signature**: `pub fn call_into(&self, params: &[f32], results: &mut [f32])`

#### JitCompiledF32::param_count

**Signature**: `pub fn param_count(&self) -> usize`

#### JitCompiledF32::result_count

**Signature**: `pub fn result_count(&self) -> usize`

### JitEngine

**Signature**:

```rust
pub struct JitEngine;
```

**Description**: the entry point of the JIT compilation engine. Provides the low-level `compile` method that compiles directly from an instruction sequence.

#### JitEngine::compile

**Signature**:

```rust
impl JitEngine {
    pub fn compile(
        instructions: &[Instr],
        param_count: usize,
        constants: &[f64],
        result_indices: &[usize],
    ) -> Result<JitCompiledFunction>

    pub fn compile_f32(
        instructions: &[Instr],
        param_count: usize,
        constants: &[f32],
        result_indices: &[usize],
    ) -> Result<JitCompiledF32>
}
```

**Description**: compiles an instruction sequence into a callable JIT function. `compile` produces a double-precision `JitCompiledFunction`; `compile_f32` produces a single-precision `JitCompiledF32` (libm calls use the `*f` symbol variants). This is the low-level implementation behind `ExpressionEvaluator::compile_jit` / `compile_jit_f32`. Internally, `compile_module` returns the JIT module (kept alive by the caller) and the entry function pointer, but only the compiled function objects are exposed publicly.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `instructions` | `&[Instr]` | the instruction sequence |
| `param_count` | `usize` | the number of parameters |
| `constants` | `&[f64]` (`compile`) / `&[f32]` (`compile_f32`) | the constant values (embedded as immediates) |
| `result_indices` | `&[usize]` | the result slot indices |

### FloatWidth

**Signature**:

```rust
pub enum FloatWidth {
    F32,
    F64,
}
```

**Description**: the floating-point width selection for JIT code generation. `F32` uses the `*f` variants of libm functions (e.g., `sinf`); `F64` uses standard libm.

---

## SIMD batch evaluation (feature = "simd")

`VectorEvaluator` uses the `pulp` crate to provide runtime-detected SIMD-vectorized evaluation. On x86_64 it automatically selects SSE2 (2-lane f64), AVX2 (4-lane f64), or AVX-512 (8-lane f64). The f32 variant has twice the lane count.

### VectorEvaluator

**Signature**:

```rust
pub struct VectorEvaluator { /* private fields */ }
```

**Description**: a double-precision SIMD batch evaluator. Splits multiple sets of input parameters into SIMD-width chunks processed in parallel; the remainder falls back to scalar.

#### ExpressionEvaluator::compile_vector_evaluator

**Signature**:

```rust
#[cfg(feature = "simd")]
impl ExpressionEvaluator<f64> {
    pub fn compile_vector_evaluator(&self) -> Result<VectorEvaluator>
}
```

**Description**: creates a SIMD batch evaluator from a compiled `ExpressionEvaluator`. Detects the SIMD width at construction and preallocates the stack and constant buffers. Returns `UnsupportedOperation` if the expression contains external functions (not supported in SIMD mode).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `self` | `&ExpressionEvaluator<f64>` | the source evaluator |

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, simd::VectorEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let simd_eval = eval.compile_vector_evaluator().unwrap();
```

#### VectorEvaluator::evaluate

**Signature**:

```rust
impl VectorEvaluator {
    pub fn evaluate(&self, params: &[Vec<f64>]) -> Result<Vec<Vec<f64>>>
}
```

**Description**: batch evaluation. `params[i]` holds all sample points of the `i`-th parameter. All `params[i]` must have the same length. Returns `Vec<Vec<f64>>` where `result[j][i]` is the `j`-th output at the `i`-th sample point.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `params` | `&[Vec<f64>]` | the parameter vectors, each `Vec` of the same length |

**Returns**: `Ok(Vec<Vec<f64>>)` or an evaluation error.

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, simd::VectorEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.mul(&[ctx.var("x"), ctx.var("x")]), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let simd_eval = eval.compile_vector_evaluator().unwrap();
let x_values: Vec<f64> = (0..1000).map(|i| i as f64 * 0.01).collect();
let results = simd_eval.evaluate(&[x_values]).unwrap();
assert_eq!(results.len(), 1);
assert_eq!(results[0].len(), 1000);
assert_eq!(results[0][0], 1.0);
// results[0] holds the 1000 values of x²+1, computed in parallel via SIMD
```

### VectorEvaluatorF32

**Signature**:

```rust
pub struct VectorEvaluatorF32 { /* private fields */ }
```

**Description**: a single-precision SIMD batch evaluator. Its lane count is twice that of `VectorEvaluator` (e.g., 8-lane f32 vs 4-lane f64 on AVX2).

#### ExpressionEvaluator::compile_vector_evaluator_f32

**Signature**:

```rust
#[cfg(feature = "simd")]
impl ExpressionEvaluator<f64> {
    pub fn compile_vector_evaluator_f32(&self) -> Result<VectorEvaluatorF32>
}
```

**Description**: creates an f32 SIMD evaluator from an f64 evaluator, narrowing the constants to f32. External functions are likewise unsupported (returns `UnsupportedOperation`).

#### VectorEvaluatorF32::evaluate

**Signature**: `pub fn evaluate(&self, params: &[Vec<f32>]) -> Result<Vec<Vec<f32>>>`

**Description**: single-precision batch evaluation; the interface is symmetric with `VectorEvaluator::evaluate`.

---

## StreamingEvaluator\<'a, T\>

**Signature**:

```rust
pub struct StreamingEvaluator<'a, T: EvaluationDomain> {
    evaluator: &'a ExpressionEvaluator<T>,
    params: Vec<T>,
    stack: Vec<T>,
    results: Vec<T>,
}
```

**Description**: a streaming evaluator that reuses internal buffers, keeping memory usage constant for input streams of arbitrary length. All scratch memory (parameter staging, evaluation stack, result buffers) is allocated once at construction and reused for every row.

### StreamingEvaluator::new

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn new(evaluator: &ExpressionEvaluator<T>) -> StreamingEvaluator<'_, T>
}
```

**Description**: creates a streaming evaluator, preallocating all buffers.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `evaluator` | `&ExpressionEvaluator<T>` | the source evaluator (borrowed) |

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, StreamingEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut stream = StreamingEvaluator::new(&eval);
```

### StreamingEvaluator::for_each

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn for_each<I, S, F>(&mut self, rows: I, mut sink: F) -> Result<usize>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[T]>,
        F: FnMut(&[T]),
}
```

**Description**: processes every row of the input stream, calling the `sink` callback with each row's results. Memory usage is constant and does not grow with the number of rows.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rows` | `I: IntoIterator<Item = S>` | the input row iterator; each row is `AsRef<[T]>` of length `param_count()` |
| `sink` | `FnMut(&[T])` | the result callback; its `&[T]` has length `result_count()` and is only valid inside the callback |

**Returns**: `Ok(number of rows processed)` or an evaluation error.

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, StreamingEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut stream = StreamingEvaluator::new(&eval);
let rows = (0..1_000_000).map(|i| [i as f64]);
let n = stream.for_each(rows, |results| {
    // results[0] is x+1 for the current row, only valid inside the callback
}).unwrap();
assert_eq!(n, 1_000_000);
```

**See also**: [StreamingEvaluator::evaluate_chunk](#streamingevaluatorevaluate_chunk)

### StreamingEvaluator::evaluate_chunk

**Signature**:

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn evaluate_chunk<S: AsRef<[T]>>(&mut self, rows: &[S]) -> Result<Vec<Vec<T>>>
}
```

**Description**: processes a batch of rows and collects all results. Suitable for bounded batches; for unbounded streams use `for_each`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rows` | `&[S]` | the input row slice |

**Returns**: `Ok(Vec<Vec<T>>)`, each `Vec<T>` holding one row's results.

**Example**:

```rust
use ocas_eval::{ExpressionEvaluator, StreamingEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut stream = StreamingEvaluator::new(&eval);
let rows: Vec<[f64; 1]> = (0..100).map(|i| [i as f64]).collect();
let all_results = stream.evaluate_chunk(&rows).unwrap();
assert_eq!(all_results.len(), 100);
assert_eq!(all_results[0][0], 1.0); // x=0 → 0+1=1
assert_eq!(all_results[99][0], 100.0); // x=99 → 99+1=100
```

---

## FunctionMap\<T\>

**Signature**:

```rust
pub struct FunctionMap<T: EvaluationDomain> {
    entries: Vec<(String, FunctionEntry<T>)>,
    name_to_idx: HashMap<String, usize>,
    aliases: HashMap<String, String>,
}
```

**Description**: a user-defined function registry. Allows calling registered external functions during expression evaluation. Name lookup is case-insensitive and supports aliases.

### FunctionMap::new

**Signature**:

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn new() -> Self
}
```

**Description**: creates an empty function map.

### FunctionMap::register

**Signature**:

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn register(&mut self, name: &str, arity: usize, func: Box<dyn Fn(&[T]) -> T + Send + Sync>)
}
```

**Description**: registers an external function.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `&str` | the function name (case-insensitive) |
| `arity` | `usize` | the expected number of arguments |
| `func` | `Box<dyn Fn(&[T]) -> T + Send + Sync>` | the function implementation |

**Example**:

```rust
use ocas_eval::FunctionMap;

let mut fmap = FunctionMap::<f64>::new();
fmap.register("square", 1, Box::new(|args| args[0] * args[0]));
fmap.register("add3", 3, Box::new(|args| args[0] + args[1] + args[2]));
```

### FunctionMap::register_alias

**Signature**:

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn register_alias(&mut self, alias: &str, canonical: &str)
}
```

**Description**: adds an alias for a registered function name. Looking up `alias` is redirected to `canonical`.

### FunctionMap::resolve

**Signature**:

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn resolve(&self, name: &str) -> Option<&FunctionEntry<T>>
}
```

**Description**: looks up a function by name (resolving aliases and case).

### FunctionMap::index_of

**Signature**: `pub fn index_of(&self, name: &str) -> Option<usize>`

**Description**: returns the function's index in the map (used for `fn_idx` in instructions).

### FunctionMap::call_by_index

**Signature**: `pub fn call_by_index(&self, idx: usize, args: &[T]) -> Option<T>`

**Description**: calls a function by index.

### FunctionMap::len / is_empty

**Signature**: `pub fn len(&self) -> usize` / `pub fn is_empty(&self) -> bool`

**Description**: returns the number of registered functions / whether the map is empty.

### FunctionEntry\<T\>

**Signature**:

```rust
pub struct FunctionEntry<T: EvaluationDomain> {
    pub arity: usize,
    // private func field
}
```

**Description**: a registered external function entry. `arity` is the expected number of arguments.

---

## Instruction set

### Slot

**Signature**:

```rust
pub enum Slot {
    Param(usize),
    Const(usize),
    Temp(usize),
}
```

**Description**: a named slot in the evaluator stack. Used by the public `Instruction` type to reference stack positions semantically rather than by raw index.

| Variant | Description |
|---|---|
| `Param(usize)` | parameter slot, 0-based index |
| `Const(usize)` | constant slot, 0-based index |
| `Temp(usize)` | temporary value slot, 0-based index |

### Instruction

**Signature**:

```rust
pub enum Instruction {
    Add(Slot, Vec<Slot>),
    Mul(Slot, Vec<Slot>),
    Pow(Slot, Slot, i64),
    Powf(Slot, Slot, Slot),
    Fun(Slot, Symbol, Slot),
    ExternalFun(Slot, usize, Vec<Slot>),
    Assign(Slot, Slot),
}
```

**Description**: a self-documenting public instruction. Unlike the internal `Instr`, it uses `Slot` instead of raw indices, making it convenient for inspection and serialization.

| Variant | Format | Description |
|---|---|---|
| `Add` | `dst = sum(sources)` | summation |
| `Mul` | `dst = product(sources)` | product |
| `Pow` | `dst = base^exp` | integer power |
| `Powf` | `dst = base^exp` | floating-point power |
| `Fun` | `dst = builtin(src)` | built-in function |
| `ExternalFun` | `dst = fns[idx](srcs...)` | external function |
| `Assign` | `dst = src` | copy |

### Instr

**Signature**:

```rust
pub enum Instr {
    Add { dst: usize, srcs: Vec<usize> },
    Mul { dst: usize, srcs: Vec<usize> },
    Pow { dst: usize, base: usize, exp: i64 },
    Powf { dst: usize, base: usize, exp: usize },
    BuiltinOp { dst: usize, op: BuiltinOp, src: usize },
    ExternalFun { dst: usize, fn_idx: usize, srcs: Vec<usize> },
    Copy { dst: usize, src: usize },
}
```

**Description**: an internal index-based instruction, executed by `ExpressionEvaluator`. All indices are absolute positions in the stack.

**Stack layout**:

```text
[params (param_count)] [constants (const_count)] [temporaries (temp_count)] [outputs (result_count)]
```

| Variant | Description |
|---|---|
| `Add { dst, srcs }` | `stack[dst] = stack[srcs[0]] + stack[srcs[1]] + ...` |
| `Mul { dst, srcs }` | `stack[dst] = stack[srcs[0]] * stack[srcs[1]] * ...` |
| `Pow { dst, base, exp }` | `stack[dst] = stack[base]^exp`, with `exp` an `i64` |
| `Powf { dst, base, exp }` | `stack[dst] = stack[base]^stack[exp]`, floating-point exponent |
| `BuiltinOp { dst, op, src }` | `stack[dst] = op(stack[src])` |
| `ExternalFun { dst, fn_idx, srcs }` | `stack[dst] = fns[fn_idx](&stack[srcs...])` |
| `Copy { dst, src }` | `stack[dst] = stack[src]` |

---

## BuiltinOp

**Signature**:

```rust
pub enum BuiltinOp {
    Sin, Cos, Tan,
    Sec, Csc, Cot,
    Exp, Log,
    Sqrt, Abs,
}
```

**Description**: pre-resolved built-in function opcodes. The compiler converts `Symbol` names into `BuiltinOp` variants at compile time, avoiding string matching on the hot path.

| Variant | Math | Description |
|---|---|---|
| `Sin` | $\sin(x)$ | sine |
| `Cos` | $\cos(x)$ | cosine |
| `Tan` | $\tan(x)$ | tangent |
| `Sec` | $\sec(x) = 1/\cos(x)$ | secant |
| `Csc` | $\csc(x) = 1/\sin(x)$ | cosecant |
| `Cot` | $\cot(x) = 1/\tan(x)$ | cotangent |
| `Exp` | $e^x$ | exponential function |
| `Log` | $\ln(x)$ | natural logarithm |
| `Sqrt` | $\sqrt{x}$ | square root |
| `Abs` | $|x|$ | absolute value |

#### BuiltinOp::from_name

**Signature**: `pub fn from_name(name: &str) -> Option<Self>`

**Description**: tries to resolve a function name into a built-in opcode. Case-insensitive. Both `"log"` and `"ln"` map to `Log`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `&str` | the function name |

**Returns**: `Some(BuiltinOp)` or `None` (not a built-in function).

**Example**:

```rust
use ocas_eval::instruction::BuiltinOp;

assert_eq!(BuiltinOp::from_name("sin"), Some(BuiltinOp::Sin));
assert_eq!(BuiltinOp::from_name("Log"), Some(BuiltinOp::Log));
assert_eq!(BuiltinOp::from_name("ln"), Some(BuiltinOp::Log));
assert_eq!(BuiltinOp::from_name("myfunc"), None);
```

---

## EvaluationDomain trait

**Signature**:

```rust
pub trait EvaluationDomain: Sized + Clone + 'static {
    fn from_f64(value: f64) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn add_ref(&self, other: &Self) -> Self;
    fn sub_ref(&self, other: &Self) -> Self;
    fn mul_ref(&self, other: &Self) -> Self;
    fn div_ref(&self, other: &Self) -> Result<Self>;
    fn neg_ref(&self) -> Self;
    fn powi_ref(&self, exp: i64) -> Self;
    fn resolve_builtin(name: &str, arg: &Self) -> Result<Self>;
}
```

**Description**: the trait for types that can serve as evaluation domains. Unlike the algebraic `Domain` trait (`Domain` is object-safe and uses `&self` receivers), `EvaluationDomain` uses static methods and value operations, making it compatible with `Copy` types such as `f64`.

**Methods**:

| Method | Signature | Description |
|---|---|---|
| `from_f64` | `fn from_f64(f64) -> Self` | create a value from an `f64` |
| `zero` | `fn zero() -> Self` | additive identity |
| `one` | `fn one() -> Self` | multiplicative identity |
| `add_ref` | `fn add_ref(&self, &Self) -> Self` | addition |
| `sub_ref` | `fn sub_ref(&self, &Self) -> Self` | subtraction |
| `mul_ref` | `fn mul_ref(&self, &Self) -> Self` | multiplication |
| `div_ref` | `fn div_ref(&self, &Self) -> Result<Self>` | division; division by zero returns `DivisionByZero` |
| `neg_ref` | `fn neg_ref(&self) -> Self` | negation |
| `powi_ref` | `fn powi_ref(&self, i64) -> Self` | integer power; returns 1 when `exp == 0` |
| `resolve_builtin` | `fn resolve_builtin(name: &str, &Self) -> Result<Self>` | execute a built-in function |

**Implementations**:

| Type | Description |
|---|---|
| `f64` | standard double precision (always available) |
| `DoubleF64` | oCAS double-precision float (~31 significant digits) |

**Example**:

```rust
use ocas_eval::EvaluationDomain;

let x = f64::from_f64(3.0);
let y = f64::from_f64(2.0);
assert_eq!(x.add_ref(&y), 5.0);
assert_eq!(f64::resolve_builtin("sin", &std::f64::consts::FRAC_PI_2).unwrap(), 1.0);
```

---

## PowfExtension trait

**Signature**:

```rust
pub trait PowfExtension: EvaluationDomain {
    fn powf_ref(&self, exp: &Self) -> Result<Self>;
}
```

**Description**: an extension trait for powers with floating-point exponents. It is separated from `EvaluationDomain` because integer domains cannot meaningfully compute $a^b$ (with non-integer $b$).

**Implementations**: `f64`, `DoubleF64`.

**Example**:

```rust
use ocas_eval::PowfExtension;

let base = 2.0_f64;
let exp = 0.5_f64;
let result = base.powf_ref(&exp).unwrap();
assert!((result - std::f64::consts::SQRT_2).abs() < 1e-15);
// Output: 2^0.5 ≈ 1.4142135623730951
```

---

## Compilation functions

The following free functions are the low-level entry points behind `ExpressionEvaluator`'s compilation methods.

### compile_atom

**Signature**:

```rust
pub fn compile_atom<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
) -> Result<ExpressionEvaluator<T>>
```

**Description**: compiles an `Atom` into an evaluator. Equivalent to `ExpressionEvaluator::compile`.

### compile_atom_with

**Signature**:

```rust
pub fn compile_atom_with<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**Description**: compiles with an optional function map.

### compile_tree / compile_tree_with

**Signature**:

```rust
pub fn compile_tree<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
) -> Result<ExpressionEvaluator<T>>

pub fn compile_tree_with<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**Description**: compiles from an `EvalTree`. Suitable when you need to apply custom optimizations to the `EvalTree` before compiling.

### compile_atoms_multi / compile_atoms_multi_with

**Signature**:

```rust
pub fn compile_atoms_multi<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
) -> Result<ExpressionEvaluator<T>>

pub fn compile_atoms_multi_with<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**Description**: compiles multiple `Atom`s into a single multi-output evaluator, sharing CSEs and constants across outputs.

### compile_trees_multi

**Signature**:

```rust
pub fn compile_trees_multi<T: EvaluationDomain + PowfExtension>(
    trees: &[&EvalTree],
) -> Result<ExpressionEvaluator<T>>
```

**Description**: compiles multiple `EvalTree`s into a single multi-output evaluator.

---

## EvaluationError

**Signature**:

```rust
#[non_exhaustive]
pub enum EvaluationError {
    UndefinedVariable { name: String },
    TypeMismatch { expected: String, found: String },
    DivisionByZero,
    FunctionNotFound { name: String },
    WrongArity { name: String, expected: usize, got: usize },
    JitCompilationError { message: String },
    UnsupportedOperation { message: String },
}
```

**Description**: the errors that can occur during evaluation.

| Variant | Description |
|---|---|
| `UndefinedVariable { name }` | the expression references a variable that was not supplied |
| `TypeMismatch { expected, found }` | type mismatch (e.g., expected a float but got an integer) |
| `DivisionByZero` | division by zero |
| `FunctionNotFound { name }` | a user-defined function was not found in the registry |
| `WrongArity { name, expected, got }` | wrong number of arguments in a function call |
| `JitCompilationError { message }` | JIT compilation failed (Cranelift backend error) |
| `UnsupportedOperation { message }` | the requested operation is not supported in the current domain |

**Display format**:

```text
undefined variable 'x'
type mismatch: expected float, found integer
division by zero
function 'f' not found
wrong arity for 'f': expected 2, got 1
JIT compilation error: <message>
unsupported operation: <message>
```

**Example**:

```rust
use ocas_eval::EvaluationError;

let err = EvaluationError::UndefinedVariable { name: "x".into() };
assert_eq!(err.to_string(), "undefined variable 'x'");

let err = EvaluationError::WrongArity {
    name: "f".into(),
    expected: 2,
    got: 1,
};
assert_eq!(err.to_string(), "wrong arity for 'f': expected 2, got 1");
```

---

## Usage patterns

### Basic evaluation

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// Build x² + 2x + 1
let x = ctx.var("x");
let expr = ctx.add(&[
    ctx.mul(&[x, x]),
    ctx.mul(&[ctx.num(2), x]),
    ctx.num(1),
]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

// Evaluate
let result = eval.evaluate(&[3.0]).unwrap();
assert_eq!(result[0], 16.0); // 9 + 6 + 1 = 16
```

### JIT hot loop

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sin", &[ctx.var("x")]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let jit = eval.compile_jit().unwrap();

let mut results = vec![0.0; jit.result_count()];
for i in 0..1_000_000 {
    jit.call_into(&[i as f64 * 0.001], &mut results);
    // use results[0]...
}
// Zero-allocation hot loop, native machine code
```

### SIMD batch processing

```rust
use ocas_eval::{ExpressionEvaluator, simd::VectorEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("exp", &[ctx.var("x")]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let simd = eval.compile_vector_evaluator().unwrap();

let x: Vec<f64> = (0..10000).map(|i| i as f64 * 0.001).collect();
let results = simd.evaluate(&[x]).unwrap();
// results[0] holds the 10000 values of e^x
```

### Streaming large datasets

```rust
use ocas_eval::{ExpressionEvaluator, StreamingEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut stream = StreamingEvaluator::new(&eval);

// Process 1 million rows with constant memory usage
let rows = (0..1_000_000).map(|i| [i as f64]);
let n = stream.for_each(rows, |_results| {
    // process results...
}).unwrap();
assert_eq!(n, 1_000_000);
```

### Custom functions

```rust
use ocas_eval::{ExpressionEvaluator, FunctionMap};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sigmoid", &[ctx.var("x")]);

let mut fmap = FunctionMap::<f64>::new();
fmap.register("sigmoid", 1, Box::new(|args| {
    1.0 / (1.0 + (-args[0]).exp())
}));

let eval = ExpressionEvaluator::compile_with(expr, fmap).unwrap();
let result = eval.evaluate(&[0.0]).unwrap();
assert!((result[0] - 0.5).abs() < 1e-15);
// Output: sigmoid(0) = 0.5
```

**See also**: [Expression system](./rust-expressions.md) (`Atom`, `AtomArena`), [Numerical integration](./rust-numeric-integration.md) (`Vegas` uses the evaluator), [Automatic differentiation](./rust-autodiff.md) (evaluation over the `HyperDual` domain)
