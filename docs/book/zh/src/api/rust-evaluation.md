# Rust API 参考：求值与 JIT

oCAS 的求值系统将符号表达式编译为栈式虚拟机指令序列，然后在数值类型上高效执行。支持三种执行后端：解释器（通用）、JIT 编译（Cranelift 生成原生机器码）、SIMD 批量求值（向量化处理多组参数）。

**模块路径**：`ocas_eval`

**编译管线**：

```text
Atom (arena-backed)
  → EvalTree (owned intermediate, 可选常量折叠优化)
  → Instr sequence (编译 + 优化 pass)
  → ExpressionEvaluator (栈式 VM 解释器)
  → 或 JitCompiledFunction (Cranelift, feature = "jit")
  → 或 VectorEvaluator (SIMD 批量, feature = "simd")
```

**Feature flags**：

| Flag | 功能 |
|---|---|
| `jit` | 启用 Cranelift JIT 编译后端 |
| `simd` | 启用 SIMD 批量求值（`pulp`，运行时检测 SSE2/AVX2/AVX-512） |
| `fast-poly` | 启用快速多项式求值优化 |

**导入方式**：

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

**签名**：

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

**功能**：arena 无关的、拥有所有权的中间表示。从 `Atom` 转换而来，断开与 arena 生命周期的绑定，使编译管线可自由进行多遍优化。

**变体**：

| 变体 | 说明 |
|---|---|
| `Num(f64)` | 数值常量 |
| `Var(String)` | 变量引用 |
| `Fun(String, Vec<EvalTree>)` | 命名函数应用 |
| `Add(Vec<EvalTree>)` | 求和 |
| `Mul(Vec<EvalTree>)` | 乘积 |
| `Pow(Box<EvalTree>, Box<EvalTree>)` | 幂运算 base^exp |

#### EvalTree::from_atom

**签名**：`pub fn from_atom(atom: Atom<'_>) -> Self`

**功能**：从 `Atom` 转换为拥有的 `EvalTree`。递归遍历表达式树，将 arena 句柄转为独立数据。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom<'_>` | 源表达式 |

**返回值**：拥有的 `EvalTree`。

**示例**：

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

**签名**：`pub fn fold_constants(&self) -> EvalTree`

**功能**：折叠常量子树并应用代数恒等式。在编译前运行，减少指令数。

**优化规则**：

| 模式 | 规则 |
|---|---|
| `Add` | 删除零项；全常量项求和；单项折叠 |
| `Mul` | 吸收零因子；删除单位因子；全常量因子求积；单因子折叠 |
| `Pow` | `x^1 → x`，`x^0 → 1`，`Num^Num` 直接求值 |
| `Fun` | 内建函数的全常量参数直接求值（外部函数不折叠，可能有副作用） |

**返回值**：优化后的新 `EvalTree`。

**示例**：

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
// 输出：3 + 0*x → 3
```

---

## ExpressionEvaluator\<T\>

**签名**：

```rust
pub struct ExpressionEvaluator<T: EvaluationDomain> { /* private fields */ }
```

**功能**：编译好的表达式，可反复执行数值求值。持有指令序列、预分配栈、常量表和可选的用户自定义函数注册表。类型参数 `T` 实现 `EvaluationDomain`，通常为 `f64`。

### 编译方法

#### ExpressionEvaluator::compile

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile(atom: Atom<'_>) -> Result<Self>
}
```

**功能**：将 `Atom` 编译为可执行求值器。内部先转为 `EvalTree`，再折叠常量、生成指令序列。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom<'_>` | 源表达式 |

**返回值**：`Ok(ExpressionEvaluator<T>)` 或编译错误。

**错误**：见 [`EvaluationError`](#evaluationerror)。

**示例**：

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

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_with(atom: Atom<'_>, map: FunctionMap<T>) -> Result<Self>
}
```

**功能**：编译表达式，并传入用户自定义函数注册表。表达式中出现的 `ExternalFun` 节点通过 `map` 解析。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom<'_>` | 源表达式 |
| `map` | `FunctionMap<T>` | 自定义函数注册表（无外部函数时传入空映射） |

**返回值**：`Ok(ExpressionEvaluator<T>)` 或编译/解析错误。

**示例**：

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
// 输出：myfunc(3) = 3² = 9
```

#### ExpressionEvaluator::compile_multi

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_multi(atoms: &[Atom<'_>]) -> Result<Self>
}
```

**功能**：将多个 `Atom` 编译为单个多输出求值器。跨所有输出共享常量、公共子表达式（CSE）和栈槽位，避免重复计算。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atoms` | `&[Atom<'_>]` | 多个表达式 |

**返回值**：`Ok(ExpressionEvaluator<T>)`，`result_count() == atoms.len()`。

**示例**：

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
// 输出：x²=9, x+1=4
```

#### ExpressionEvaluator::compile_multi_with

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn compile_multi_with(atoms: &[Atom<'_>], map: FunctionMap<T>) -> Result<Self>
}
```

**功能**：`compile_multi` 的带函数映射版本（`map` 直接传入，而非 `Option`）。

### 执行方法

#### ExpressionEvaluator::evaluate

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    pub fn evaluate(&self, params: &[T]) -> Result<Vec<T>>
}
```

**功能**：用给定参数值执行编译好的表达式。每次调用分配新的栈和结果缓冲区。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `params` | `&[T]` | 参数值，长度必须等于 `param_count()` |

**返回值**：`Ok(Vec<T>)`，长度等于 `result_count()`。

**错误**：

| 错误 | 条件 |
|---|---|
| `WrongArity` | `params.len() != param_count()` |
| `DivisionByZero` | 求值过程中遇到除零 |
| `FunctionNotFound` | 外部函数未注册 |

**示例**：

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
// 输出：3² + 1 = 10
```

#### ExpressionEvaluator::evaluate_with_stack

**签名**：

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

**功能**：使用调用者提供的缓冲区执行求值，避免每次调用的堆分配。流式和批处理场景应优先使用此方法。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `params` | `&[T]` | 参数值 |
| `stack` | `&mut Vec<T>` | 栈缓冲区（自动 resize 到 `stack_size()`） |
| `results` | `&mut Vec<T>` | 结果缓冲区（自动填充） |

**返回值**：`Ok(())` 或求值错误。

### 访问器

#### ExpressionEvaluator::param_count

**签名**：`pub fn param_count(&self) -> usize`

**功能**：返回期望的参数个数。

#### ExpressionEvaluator::result_count

**签名**：`pub fn result_count(&self) -> usize`

**功能**：返回输出结果个数（单输出表达式为 1，`compile_multi` 创建的可能大于 1）。

#### ExpressionEvaluator::stack_size

**签名**：`pub fn stack_size(&self) -> usize`

**功能**：返回求值所需栈大小（参数 + 常量 + 临时值 + 输出）。`evaluate_with_stack` 的 `stack` 参数应预分配此容量。

---

## JIT 编译（feature = "jit"）

JIT 后端使用 Cranelift 将指令序列编译为原生机器码，常量嵌入为立即数，调用 libm 进行超越函数计算。性能通常比解释器快 10 倍以上。

### ExpressionEvaluator::compile_jit

**签名**：

```rust
#[cfg(feature = "jit")]
impl ExpressionEvaluator<f64> {
    pub fn compile_jit(&self) -> Result<JitCompiledFunction>
}
```

**功能**：将当前求值器的指令序列编译为双精度（f64）原生代码。

**返回值**：`Ok(JitCompiledFunction)` 或编译错误。

**错误**：

| 错误 | 条件 |
|---|---|
| `UnsupportedOperation` | 包含 JIT 无法降低的指令（如外部函数、`sec`/`csc`/`cot`） |
| `JitCompilationError` | Cranelift 后端编译失败 |

**示例**：

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
// 输出：3² + 1 = 10（原生机器码执行）
```

### ExpressionEvaluator::compile_jit_f32

**签名**：

```rust
#[cfg(feature = "jit")]
impl ExpressionEvaluator<f64> {
    pub fn compile_jit_f32(&self) -> Result<JitCompiledF32>
}
```

**功能**：将指令序列编译为单精度（f32）原生代码。常量从 f64 窄化为 f32；结果精度为 f32。

**错误**：与 `compile_jit` 相同。

**示例**：

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
// 输出：3 × 2 = 6（f32 精度）
```

### JitCompiledFunction

**签名**：

```rust
pub struct JitCompiledFunction { /* private fields */ }
```

**功能**：JIT 编译好的双精度函数，可反复调用。

#### JitCompiledFunction::call

**签名**：`pub fn call(&self, params: &[f64]) -> Vec<f64>`

**功能**：用给定参数调用 JIT 编译的函数。每次分配新的结果向量。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `params` | `&[f64]` | 参数值 |

**返回值**：结果向量，长度等于 `result_count()`。

**示例**：

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
// 输出：x²+1 在 x=3 时的值 10
```

#### JitCompiledFunction::call_into

**签名**：`pub fn call_into(&self, params: &[f64], results: &mut [f64])`

**功能**：调用 JIT 函数，结果写入调用者提供的缓冲区，避免每次调用的堆分配。`results.len()` 必须至少为 `result_count()`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `params` | `&[f64]` | 参数值 |
| `results` | `&mut [f64]` | 结果缓冲区 |

**示例**：

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
// 结果在 buf 中（42² + 1 = 1765），无堆分配
```

#### JitCompiledFunction::param_count

**签名**：`pub fn param_count(&self) -> usize`

**功能**：返回期望的参数个数。

#### JitCompiledFunction::result_count

**签名**：`pub fn result_count(&self) -> usize`

**功能**：返回结果个数。

### JitCompiledF32

**签名**：

```rust
pub struct JitCompiledF32 { /* private fields */ }
```

**功能**：JIT 编译好的单精度函数。API 与 `JitCompiledFunction` 完全对称，参数和返回值为 `f32`。

#### JitCompiledF32::call

**签名**：`pub fn call(&self, params: &[f32]) -> Vec<f32>`

#### JitCompiledF32::call_into

**签名**：`pub fn call_into(&self, params: &[f32], results: &mut [f32])`

#### JitCompiledF32::param_count

**签名**：`pub fn param_count(&self) -> usize`

#### JitCompiledF32::result_count

**签名**：`pub fn result_count(&self) -> usize`

### JitEngine

**签名**：

```rust
pub struct JitEngine;
```

**功能**：JIT 编译引擎的入口点。提供底层 `compile` 方法，可直接从指令序列编译。

#### JitEngine::compile

**签名**：

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

**功能**：将指令序列编译为可调用的 JIT 函数。`compile` 生成双精度 `JitCompiledFunction`，`compile_f32` 生成单精度 `JitCompiledF32`（libm 调用使用 `*f` 符号变体）。这是 `ExpressionEvaluator::compile_jit` / `compile_jit_f32` 的底层实现。内部 `compile_module` 返回 JIT 模块（由调用者保持存活）与入口函数指针，但对外只暴露编译好的函数对象。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `instructions` | `&[Instr]` | 指令序列 |
| `param_count` | `usize` | 参数个数 |
| `constants` | `&[f64]`（`compile`）/ `&[f32]`（`compile_f32`） | 常量值（嵌入为立即数） |
| `result_indices` | `&[usize]` | 结果槽索引 |

### FloatWidth

**签名**：

```rust
pub enum FloatWidth {
    F32,
    F64,
}
```

**功能**：JIT 代码生成的浮点宽度选择。`F32` 使用 `*f` 变体的 libm 函数（如 `sinf`），`F64` 使用标准 libm。

---

## SIMD 批量求值（feature = "simd"）

`VectorEvaluator` 使用 `pulp` 库实现运行时检测的 SIMD 向量化求值。在 x86_64 上自动选择 SSE2（2 路 f64）、AVX2（4 路 f64）或 AVX-512（8 路 f64）。f32 版本的通道数翻倍。

### VectorEvaluator

**签名**：

```rust
pub struct VectorEvaluator { /* private fields */ }
```

**功能**：双精度 SIMD 批量求值器。将多组输入参数分为 SIMD 宽度的 chunk 并行处理，余量回退到标量。

#### ExpressionEvaluator::compile_vector_evaluator

**签名**：

```rust
#[cfg(feature = "simd")]
impl ExpressionEvaluator<f64> {
    pub fn compile_vector_evaluator(&self) -> Result<VectorEvaluator>
}
```

**功能**：从编译好的 `ExpressionEvaluator` 创建 SIMD 批量求值器。构造时检测 SIMD 宽度，预分配栈和常量缓冲区。若表达式包含外部函数（SIMD 模式不支持），返回 `UnsupportedOperation`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `self` | `&ExpressionEvaluator<f64>` | 源求值器 |

**示例**：

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

**签名**：

```rust
impl VectorEvaluator {
    pub fn evaluate(&self, params: &[Vec<f64>]) -> Result<Vec<Vec<f64>>>
}
```

**功能**：批量求值。`params[i]` 是第 `i` 个参数的所有采样点。所有 `params[i]` 长度必须相同。返回 `Vec<Vec<f64>>`，其中 `result[j][i]` 是第 `i` 个采样点的第 `j` 个输出。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `params` | `&[Vec<f64>]` | 参数向量，每个 `Vec` 长度相同 |

**返回值**：`Ok(Vec<Vec<f64>>)` 或求值错误。

**示例**：

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
// results[0] 包含 1000 个 x²+1 的值，SIMD 并行计算
```

### VectorEvaluatorF32

**签名**：

```rust
pub struct VectorEvaluatorF32 { /* private fields */ }
```

**功能**：单精度 SIMD 批量求值器。通道数是 `VectorEvaluator` 的两倍（如 AVX2 上为 8 路 f32 vs 4 路 f64）。

#### ExpressionEvaluator::compile_vector_evaluator_f32

**签名**：

```rust
#[cfg(feature = "simd")]
impl ExpressionEvaluator<f64> {
    pub fn compile_vector_evaluator_f32(&self) -> Result<VectorEvaluatorF32>
}
```

**功能**：从 f64 求值器创建 f32 SIMD 求值器，常量窄化为 f32。外部函数同样不支持（返回 `UnsupportedOperation`）。

#### VectorEvaluatorF32::evaluate

**签名**：`pub fn evaluate(&self, params: &[Vec<f32>]) -> Result<Vec<Vec<f32>>>`

**功能**：单精度批量求值，接口与 `VectorEvaluator::evaluate` 对称。

---

## StreamingEvaluator\<'a, T\>

**签名**：

```rust
pub struct StreamingEvaluator<'a, T: EvaluationDomain> {
    evaluator: &'a ExpressionEvaluator<T>,
    params: Vec<T>,
    stack: Vec<T>,
    results: Vec<T>,
}
```

**功能**：流式求值器，复用内部缓冲区，处理任意长的输入流时内存使用恒定。所有 scratch 内存（参数暂存、求值栈、结果缓冲区）在创建时一次性分配，每行复用。

### StreamingEvaluator::new

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn new(evaluator: &ExpressionEvaluator<T>) -> StreamingEvaluator<'_, T>
}
```

**功能**：创建流式求值器，预分配所有缓冲区。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `evaluator` | `&ExpressionEvaluator<T>` | 源求值器（借用） |

**示例**：

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

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn for_each<I, S, F>(&mut self, rows: I, mut sink: F) -> Result<usize>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[T]>,
        F: FnMut(&[T]),
}
```

**功能**：处理输入流的每一行，对每行结果调用 `sink` 回调。内存使用恒定，不随行数增长。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rows` | `I: IntoIterator<Item = S>` | 输入行迭代器，每行 `AsRef<[T]>` 长度等于 `param_count()` |
| `sink` | `FnMut(&[T])` | 结果回调，`&[T]` 长度等于 `result_count()`，仅在回调内有效 |

**返回值**：`Ok(处理的行数)` 或求值错误。

**示例**：

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
    // results[0] 为当前行的 x+1，仅在回调内有效
}).unwrap();
assert_eq!(n, 1_000_000);
```

**参见**：[StreamingEvaluator::evaluate_chunk](#streamingevaluatorevaluate_chunk)

### StreamingEvaluator::evaluate_chunk

**签名**：

```rust
impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    pub fn evaluate_chunk<S: AsRef<[T]>>(&mut self, rows: &[S]) -> Result<Vec<Vec<T>>>
}
```

**功能**：处理一批行并收集所有结果。适用于有界批次；无界流应使用 `for_each`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rows` | `&[S]` | 输入行切片 |

**返回值**：`Ok(Vec<Vec<T>>)`，每个 `Vec<T>` 对应一行的结果。

**示例**：

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

**签名**：

```rust
pub struct FunctionMap<T: EvaluationDomain> {
    entries: Vec<(String, FunctionEntry<T>)>,
    name_to_idx: HashMap<String, usize>,
    aliases: HashMap<String, String>,
}
```

**功能**：用户自定义函数注册表。允许在表达式求值期间调用注册的外部函数。名称查找不区分大小写，支持别名。

### FunctionMap::new

**签名**：

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn new() -> Self
}
```

**功能**：创建空的函数映射。

### FunctionMap::register

**签名**：

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn register(&mut self, name: &str, arity: usize, func: Box<dyn Fn(&[T]) -> T + Send + Sync>)
}
```

**功能**：注册一个外部函数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `&str` | 函数名（不区分大小写） |
| `arity` | `usize` | 期望的参数个数 |
| `func` | `Box<dyn Fn(&[T]) -> T + Send + Sync>` | 函数实现 |

**示例**：

```rust
use ocas_eval::FunctionMap;

let mut fmap = FunctionMap::<f64>::new();
fmap.register("square", 1, Box::new(|args| args[0] * args[0]));
fmap.register("add3", 3, Box::new(|args| args[0] + args[1] + args[2]));
```

### FunctionMap::register_alias

**签名**：

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn register_alias(&mut self, alias: &str, canonical: &str)
}
```

**功能**：为已注册的函数名添加别名。查找 `alias` 时会重定向到 `canonical`。

### FunctionMap::resolve

**签名**：

```rust
impl<T: EvaluationDomain> FunctionMap<T> {
    pub fn resolve(&self, name: &str) -> Option<&FunctionEntry<T>>
}
```

**功能**：按名称查找函数（解析别名和大小写）。

### FunctionMap::index_of

**签名**：`pub fn index_of(&self, name: &str) -> Option<usize>`

**功能**：返回函数在映射中的索引（用于指令中的 `fn_idx`）。

### FunctionMap::call_by_index

**签名**：`pub fn call_by_index(&self, idx: usize, args: &[T]) -> Option<T>`

**功能**：按索引调用函数。

### FunctionMap::len / is_empty

**签名**：`pub fn len(&self) -> usize` / `pub fn is_empty(&self) -> bool`

**功能**：返回已注册函数数量 / 判断是否为空。

### FunctionEntry\<T\>

**签名**：

```rust
pub struct FunctionEntry<T: EvaluationDomain> {
    pub arity: usize,
    // private func field
}
```

**功能**：已注册的外部函数条目。`arity` 为期望的参数个数。

---

## 指令集

### Slot

**签名**：

```rust
pub enum Slot {
    Param(usize),
    Const(usize),
    Temp(usize),
}
```

**功能**：求值器栈中的命名槽位。用于公开的 `Instruction` 类型，以语义而非原始索引引用栈位置。

| 变体 | 说明 |
|---|---|
| `Param(usize)` | 参数槽，0 基索引 |
| `Const(usize)` | 常量槽，0 基索引 |
| `Temp(usize)` | 临时值槽，0 基索引 |

### Instruction

**签名**：

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

**功能**：自文档化的公开指令。与内部 `Instr` 不同，使用 `Slot` 而非原始索引，便于检查和序列化。

| 变体 | 格式 | 说明 |
|---|---|---|
| `Add` | `dst = sum(sources)` | 求和 |
| `Mul` | `dst = product(sources)` | 求积 |
| `Pow` | `dst = base^exp` | 整数次幂 |
| `Powf` | `dst = base^exp` | 浮点次幂 |
| `Fun` | `dst = builtin(src)` | 内建函数 |
| `ExternalFun` | `dst = fns[idx](srcs...)` | 外部函数 |
| `Assign` | `dst = src` | 复制 |

### Instr

**签名**：

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

**功能**：内部索引式指令，由 `ExpressionEvaluator` 执行。所有索引是栈中的绝对位置。

**栈布局**：

```text
[params (param_count)] [constants (const_count)] [temporaries (temp_count)] [outputs (result_count)]
```

| 变体 | 说明 |
|---|---|
| `Add { dst, srcs }` | `stack[dst] = stack[srcs[0]] + stack[srcs[1]] + ...` |
| `Mul { dst, srcs }` | `stack[dst] = stack[srcs[0]] * stack[srcs[1]] * ...` |
| `Pow { dst, base, exp }` | `stack[dst] = stack[base]^exp`，`exp` 为 `i64` |
| `Powf { dst, base, exp }` | `stack[dst] = stack[base]^stack[exp]`，浮点指数 |
| `BuiltinOp { dst, op, src }` | `stack[dst] = op(stack[src])` |
| `ExternalFun { dst, fn_idx, srcs }` | `stack[dst] = fns[fn_idx](&stack[srcs...])` |
| `Copy { dst, src }` | `stack[dst] = stack[src]` |

---

## BuiltinOp

**签名**：

```rust
pub enum BuiltinOp {
    Sin, Cos, Tan,
    Sec, Csc, Cot,
    Exp, Log,
    Sqrt, Abs,
}
```

**功能**：预解析的内建函数操作码。编译器在编译时将 `Symbol` 名称转换为 `BuiltinOp` 变体，避免在热路径上进行字符串匹配。

| 变体 | 数学 | 说明 |
|---|---|---|
| `Sin` | $\sin(x)$ | 正弦 |
| `Cos` | $\cos(x)$ | 余弦 |
| `Tan` | $\tan(x)$ | 正切 |
| `Sec` | $\sec(x) = 1/\cos(x)$ | 正割 |
| `Csc` | $\csc(x) = 1/\sin(x)$ | 余割 |
| `Cot` | $\cot(x) = 1/\tan(x)$ | 余切 |
| `Exp` | $e^x$ | 指数函数 |
| `Log` | $\ln(x)$ | 自然对数 |
| `Sqrt` | $\sqrt{x}$ | 平方根 |
| `Abs` | $|x|$ | 绝对值 |

#### BuiltinOp::from_name

**签名**：`pub fn from_name(name: &str) -> Option<Self>`

**功能**：尝试从函数名解析为内建操作码。不区分大小写。`"log"` 和 `"ln"` 均映射到 `Log`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `&str` | 函数名 |

**返回值**：`Some(BuiltinOp)` 或 `None`（非内建函数）。

**示例**：

```rust
use ocas_eval::instruction::BuiltinOp;

assert_eq!(BuiltinOp::from_name("sin"), Some(BuiltinOp::Sin));
assert_eq!(BuiltinOp::from_name("Log"), Some(BuiltinOp::Log));
assert_eq!(BuiltinOp::from_name("ln"), Some(BuiltinOp::Log));
assert_eq!(BuiltinOp::from_name("myfunc"), None);
```

---

## EvaluationDomain trait

**签名**：

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

**功能**：可作为求值域的类型的 trait。与代数 `Domain` trait 不同（`Domain` 是对象安全的、使用 `&self` 接收器），`EvaluationDomain` 使用静态方法和值操作，兼容 `f64` 等 `Copy` 类型。

**方法**：

| 方法 | 签名 | 说明 |
|---|---|---|
| `from_f64` | `fn from_f64(f64) -> Self` | 从 `f64` 创建值 |
| `zero` | `fn zero() -> Self` | 加法单位元 |
| `one` | `fn one() -> Self` | 乘法单位元 |
| `add_ref` | `fn add_ref(&self, &Self) -> Self` | 加法 |
| `sub_ref` | `fn sub_ref(&self, &Self) -> Self` | 减法 |
| `mul_ref` | `fn mul_ref(&self, &Self) -> Self` | 乘法 |
| `div_ref` | `fn div_ref(&self, &Self) -> Result<Self>` | 除法，除零返回 `DivisionByZero` |
| `neg_ref` | `fn neg_ref(&self) -> Self` | 取反 |
| `powi_ref` | `fn powi_ref(&self, i64) -> Self` | 整数次幂，`exp == 0` 时返回 1 |
| `resolve_builtin` | `fn resolve_builtin(name: &str, &Self) -> Result<Self>` | 执行内建函数 |

**实现**：

| 类型 | 说明 |
|---|---|
| `f64` | 标准双精度浮点（始终可用） |
| `DoubleF64` | oCAS 双精度浮点（~31 位有效数字） |

**示例**：

```rust
use ocas_eval::EvaluationDomain;

let x = f64::from_f64(3.0);
let y = f64::from_f64(2.0);
assert_eq!(x.add_ref(&y), 5.0);
assert_eq!(f64::resolve_builtin("sin", &std::f64::consts::FRAC_PI_2).unwrap(), 1.0);
```

---

## PowfExtension trait

**签名**：

```rust
pub trait PowfExtension: EvaluationDomain {
    fn powf_ref(&self, exp: &Self) -> Result<Self>;
}
```

**功能**：浮点指数幂的扩展 trait。从 `EvaluationDomain` 分离出来，因为整数域无法有意义地计算 $a^b$（$b$ 非整数）。

**实现**：`f64`、`DoubleF64`。

**示例**：

```rust
use ocas_eval::PowfExtension;

let base = 2.0_f64;
let exp = 0.5_f64;
let result = base.powf_ref(&exp).unwrap();
assert!((result - std::f64::consts::SQRT_2).abs() < 1e-15);
// 输出：2^0.5 ≈ 1.4142135623730951
```

---

## 编译函数

以下自由函数是 `ExpressionEvaluator` 编译方法的底层入口。

### compile_atom

**签名**：

```rust
pub fn compile_atom<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
) -> Result<ExpressionEvaluator<T>>
```

**功能**：将 `Atom` 编译为求值器。等价于 `ExpressionEvaluator::compile`。

### compile_atom_with

**签名**：

```rust
pub fn compile_atom_with<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**功能**：带可选函数映射的编译。

### compile_tree / compile_tree_with

**签名**：

```rust
pub fn compile_tree<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
) -> Result<ExpressionEvaluator<T>>

pub fn compile_tree_with<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**功能**：从 `EvalTree` 编译。适用于需要先对 `EvalTree` 进行自定义优化再编译的场景。

### compile_atoms_multi / compile_atoms_multi_with

**签名**：

```rust
pub fn compile_atoms_multi<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
) -> Result<ExpressionEvaluator<T>>

pub fn compile_atoms_multi_with<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>>
```

**功能**：将多个 `Atom` 编译为单个多输出求值器。跨输出共享 CSE 和常量。

### compile_trees_multi

**签名**：

```rust
pub fn compile_trees_multi<T: EvaluationDomain + PowfExtension>(
    trees: &[&EvalTree],
) -> Result<ExpressionEvaluator<T>>
```

**功能**：将多个 `EvalTree` 编译为单个多输出求值器。

---

## EvaluationError

**签名**：

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

**功能**：求值过程中可能产生的错误。

| 变体 | 说明 |
|---|---|
| `UndefinedVariable { name }` | 表达式引用了未提供的变量 |
| `TypeMismatch { expected, found }` | 类型不匹配（如期望浮点却得到整数） |
| `DivisionByZero` | 除零 |
| `FunctionNotFound { name }` | 用户自定义函数未在注册表中找到 |
| `WrongArity { name, expected, got }` | 函数调用参数个数错误 |
| `JitCompilationError { message }` | JIT 编译失败（Cranelift 后端错误） |
| `UnsupportedOperation { message }` | 请求的操作在当前域不支持 |

**Display 格式**：

```text
undefined variable 'x'
type mismatch: expected float, found integer
division by zero
function 'f' not found
wrong arity for 'f': expected 2, got 1
JIT compilation error: <message>
unsupported operation: <message>
```

**示例**：

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

## 使用模式

### 基本求值

```rust
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// 构建 x² + 2x + 1
let x = ctx.var("x");
let expr = ctx.add(&[
    ctx.mul(&[x, x]),
    ctx.mul(&[ctx.num(2), x]),
    ctx.num(1),
]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

// 求值
let result = eval.evaluate(&[3.0]).unwrap();
assert_eq!(result[0], 16.0); // 9 + 6 + 1 = 16
```

### JIT 热循环

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
    // 使用 results[0]...
}
// 零分配热循环，原生机器码执行
```

### SIMD 批量处理

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
// results[0] 包含 10000 个 e^x 的值
```

### 流式处理大数据集

```rust
use ocas_eval::{ExpressionEvaluator, StreamingEvaluator};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.add(&[ctx.var("x"), ctx.num(1)]);

let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut stream = StreamingEvaluator::new(&eval);

// 处理 100 万行，内存使用恒定
let rows = (0..1_000_000).map(|i| [i as f64]);
let n = stream.for_each(rows, |_results| {
    // 处理 results...
}).unwrap();
assert_eq!(n, 1_000_000);
```

### 自定义函数

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
// 输出：sigmoid(0) = 0.5
```

**参见**：[表达式系统](./rust-expressions.md)（`Atom`、`AtomArena`）、[数值积分](./rust-numeric-integration.md)（`Vegas` 使用求值器）、[自动微分](./rust-autodiff.md)（`HyperDual` 域上的求值）
