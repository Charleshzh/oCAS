# Rust API 参考：张量

oCAS 提供带显式指标管理的基础张量代数。`Tensor` 携带命名指标槽（上标/下标）、可选的对称性元数据，并可转换为 `Atom` 表达式进行符号处理。自 0.22.0 起，额外提供基于图同构的规范化引擎和 Young 投影子。

**模块路径**：`ocas_atom::tensor`（`mod.rs`）、`ocas_atom::tensor::spec`、`ocas_atom::tensor::canon`、`ocas_atom::tensor::dummy`、`ocas_atom::tensor::young`

**导入方式**：

```rust
use ocas_atom::tensor::{
    IndexPosition, IndexSlot, Symmetry, Tensor,
    contract, Contracted, TensorProduct,
    symmetrise_sign,
};
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::{canonicalize_tensors, CanonicalTensor, TensorCanonError};
use ocas_atom::tensor::dummy::{refresh_dummies, DummyError};
use ocas_atom::tensor::young::{YoungTableau, young_project};
```

---

## IndexPosition

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IndexPosition {
    Upper,
    Lower,
}
```

**功能**：指标在张量中的位置（方差）。`Upper` 表示逆变指标（上标），`Lower` 表示协变指标（下标）。缩并要求两个相同标签的指标具有相反的位置。

**变体**：

| 变体 | 含义 |
|---|---|
| `Upper` | 逆变指标（上标），物理中对应反变向量 |
| `Lower` | 协变指标（下标），物理中对应协变向量 |

**示例**：

```rust
use ocas_atom::tensor::IndexPosition;

assert_eq!(IndexPosition::Upper, IndexPosition::Upper);
assert_ne!(IndexPosition::Upper, IndexPosition::Lower);
// 输出：所有断言通过
```

---

## IndexSlot

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexSlot<'a> { /* private fields */ }
```

**功能**：张量的单个指标槽，由指标标签（`Atom`）和位置（`IndexPosition`）组成。`Copy` 类型，可廉价复制。

### IndexSlot::new

**签名**：`pub fn new(label: Atom<'a>, position: IndexPosition) -> Self`

**功能**：创建一个新的指标槽。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `label` | `Atom<'a>` | 指标标签表达式（通常为变量名如 `mu`、`i`） |
| `position` | `IndexPosition` | 指标位置：`Upper` 或 `Lower` |

**返回值**：`IndexSlot<'a>`

**示例**：

```rust
use ocas_atom::tensor::{IndexSlot, IndexPosition};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mu = IndexSlot::new(ctx.var("mu"), IndexPosition::Upper);
let nu = IndexSlot::new(ctx.var("nu"), IndexPosition::Lower);
assert_eq!(mu.position(), IndexPosition::Upper);
assert_eq!(nu.position(), IndexPosition::Lower);
// 输出：所有断言通过
```

### IndexSlot::label

**签名**：`pub fn label(&self) -> Atom<'a>`

**功能**：返回指标标签表达式。

**返回值**：`Atom<'a>` — 指标标签的 `Atom` 句柄。

### IndexSlot::position

**签名**：`pub fn position(&self) -> IndexPosition`

**功能**：返回指标的位置（`Upper` 或 `Lower`）。

**返回值**：`IndexPosition`

**参见**：[IndexPosition](#indexposition)

---

## Symmetry

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Symmetry {
    None,
    Symmetric,
    Antisymmetric,
}
```

**功能**：张量指标槽的对称性。这是**建议性**元数据——`contract` 等运算不会自动对称化；对称性由下游消费者（如规范化引擎或 Young 投影子）使用。

**变体**：

| 变体 | 含义 |
|---|---|
| `None` | 无对称性（一般张量） |
| `Symmetric` | 对称：任意交换指标槽不改变张量 |
| `Antisymmetric` | 反对称：交换任意两个指标槽改变符号 |

**示例**：

```rust
use ocas_atom::tensor::Symmetry;

assert_eq!(Symmetry::None, Symmetry::None);
assert_ne!(Symmetry::Symmetric, Symmetry::Antisymmetric);
// 输出：所有断言通过
```

---

## Tensor

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tensor<'a> { /* private fields */ }
```

**功能**：命名张量对象，由名称（`Symbol`）、指标槽列表（`Vec<IndexSlot>`）和对称性（`Symmetry`）组成。可通过 `to_atom` 降级为标准 `Atom` 表达式节点。

### Tensor::new

**签名**：`pub fn new(name: Symbol, slots: Vec<IndexSlot<'a>>) -> Self`

**功能**：创建一个新张量，默认对称性为 `Symmetry::None`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `Symbol` | 张量名称（interned 字符串，如 `"g"`、`"Riemann"`） |
| `slots` | `Vec<IndexSlot<'a>>` | 指标槽列表，顺序有意义 |

**返回值**：`Tensor<'a>`

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
];
let t = Tensor::new(Symbol::new("g"), slots);
assert_eq!(t.rank(), 2);
assert_eq!(t.name().as_str(), "g");
// 输出：所有断言通过
```

### Tensor::with_symmetry

**签名**：`pub fn with_symmetry(mut self, symmetry: Symmetry) -> Self`

**功能**：Builder 方法，设置张量的对称性。消费并返回 `self`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `symmetry` | `Symmetry` | 目标对称性 |

**返回值**：`Self` — 设置了对称性的张量。

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
];
let sym = Tensor::new(Symbol::new("g"), slots.clone())
    .with_symmetry(Symmetry::Symmetric);
let anti = Tensor::new(Symbol::new("epsilon"), slots)
    .with_symmetry(Symmetry::Antisymmetric);

assert_eq!(sym.symmetry(), Symmetry::Symmetric);
assert_eq!(anti.symmetry(), Symmetry::Antisymmetric);
// 输出：所有断言通过
```

### Tensor::name

**签名**：`pub fn name(&self) -> Symbol`

**功能**：返回张量名称。

**返回值**：`Symbol` — interned 字符串句柄。

### Tensor::slots

**签名**：`pub fn slots(&self) -> &[IndexSlot<'a>]`

**功能**：返回指标槽的切片引用。

**返回值**：`&[IndexSlot<'a>]`

### Tensor::symmetry

**签名**：`pub fn symmetry(&self) -> Symmetry`

**功能**：返回张量的对称性。

**返回值**：`Symmetry`

### Tensor::rank

**签名**：`pub fn rank(&self) -> usize`

**功能**：返回张量的阶（指标槽数量）。

**返回值**：`usize`

### Tensor::dummy_labels

**签名**：`pub fn dummy_labels(&self) -> Vec<Atom<'a>>`

**功能**：返回哑指标标签列表。哑指标定义为在所有槽中恰好出现两次的标签（一次 `Upper`、一次 `Lower`），即缩并候选。

**返回值**：`Vec<Atom<'a>>` — 哑指标标签列表。

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ_μ：mu 出现两次（Upper + Lower）→ 哑指标
let t = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
]);
assert_eq!(t.dummy_labels().len(), 1);
// 输出：1 个哑指标
```

### Tensor::to_atom

**签名**：`pub fn to_atom(&self, ctx: &'a AtomArena<'a>) -> Atom<'a>`

**功能**：将张量降级为标准 `Atom` 函数节点 `name(slot₁, slot₂, …)`。不应用对称化——原子保留 `self` 的槽顺序。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |

**返回值**：`Atom<'a>` — 函数节点表达式。

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let t = Tensor::new(Symbol::new("g"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
]);
let atom = t.to_atom(&ctx);
println!("{}", atom);  // g(mu, nu)
// 输出：g(mu, nu)
```

**参见**：[Atom](./rust-expressions.md#atom)

---

## contract

**签名**：

```rust
pub fn contract<'a>(
    ctx: &'a AtomArena<'a>,
    a: &Tensor<'a>,
    b: &Tensor<'a>,
) -> Contracted<'a>
```

**功能**：通过求和共享哑指标来缩并两个张量。两个具有相同标签但相反位置的槽会缩并。结果保留存活的自由指标，顺序为 `(a 的自由槽, b 的自由槽)`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `a` | `&Tensor<'a>` | 第一个张量 |
| `b` | `&Tensor<'a>` | 第二个张量 |

**返回值**：`Contracted<'a>`
- 无缩并对时返回 `Contracted::Product`，包含原始两个张量
- 部分缩并时返回 `Contracted::Product`，包含一个携带自由槽的新张量
- 完全缩并（无自由槽）时返回 `Contracted::Scalar`

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, contract, Contracted};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ B_μ → 标量
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
]);
match contract(&ctx, &a, &b) {
    Contracted::Scalar(expr) => println!("scalar: {}", expr),
    Contracted::Product(tp) => {
        for f in &tp.factors {
            println!("free: {} (rank {})", f.name().as_str(), f.rank());
        }
    }
}
// 输出：scalar: A(mu)*B(mu)
```

**参见**：[Contracted](#contracted)、[TensorProduct](#tensorproduct)

---

## Contracted

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Contracted<'a> {
    Product(TensorProduct<'a>),
    Scalar(Atom<'a>),
}
```

**功能**：缩并两个张量的结果。当无自由指标存活时为标量表达式；否则为含自由槽的张量积。

**变体**：

| 变体 | 类型 | 说明 |
|---|---|---|
| `Product` | `TensorProduct<'a>` | 部分缩并或无缩并，含自由指标的张量积 |
| `Scalar` | `Atom<'a>` | 完全缩并为标量表达式 |

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, contract, Contracted};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ_ν B^ν_ρ → 部分缩并，剩余 A^μ_ρ
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("rho"), IndexPosition::Lower),
]);
match contract(&ctx, &a, &b) {
    Contracted::Product(tp) => {
        assert_eq!(tp.factors.len(), 1);
        assert_eq!(tp.factors[0].rank(), 2); // mu, rho
    }
    Contracted::Scalar(_) => panic!("expected partial contraction"),
}
// 输出：部分缩并，剩余 2 个自由指标
```

---

## TensorProduct

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TensorProduct<'a> {
    pub factors: Vec<Tensor<'a>>,
}
```

**功能**：缩并后保留自由槽的张量积。`factors` 包含缩并后的存活张量，其槽按 `(a 的自由槽, b 的自由槽)` 顺序拼接。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `factors` | `Vec<Tensor<'a>>` | 缩并后的张量因子列表 |

**参见**：[contract](#contract)、[Contracted](#contracted)

---

## symmetrise_sign

**签名**：`pub fn symmetrise_sign(tensor: &Tensor<'_>) -> i64`

**功能**：对张量的指标槽按标签升序排列，并返回反对称化的符号。

- `Symmetry::None` 或 `Symmetry::Symmetric`：返回 `+1`
- `Symmetry::Antisymmetric`：返回排列奇偶性（偶排列 `+1`，奇排列 `-1`）

**注意**：这不是置换群下的完全规范化（那需要图同构）；仅为对称张量提供稳定的等价比较顺序。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `tensor` | `&Tensor<'_>` | 输入张量 |

**返回值**：`i64` — `+1` 或 `-1`。

**示例**：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry, symmetrise_sign};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("a"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("b"), IndexPosition::Upper),
];
let anti = Tensor::new(Symbol::new("epsilon"), slots)
    .with_symmetry(Symmetry::Antisymmetric);

let sign = symmetrise_sign(&anti);
// 若 a < b 则排列不变，sign = +1
// 输出：+1
```

---

## SymmetrySpec

**签名**：

```rust
#[derive(Debug, Clone, Default)]
pub struct SymmetrySpec {
    pub symmetric_subsets: Vec<Vec<usize>>,
    pub antisymmetric_subsets: Vec<Vec<usize>>,
    pub cyclic: Option<Vec<usize>>,
}
```

**功能**：张量的细粒度对称性规格，用于规范化引擎。比 `Symmetry` 枚举更灵活——可以为不同的槽位子集指定不同的对称行为。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `symmetric_subsets` | `Vec<Vec<usize>>` | 对称的槽位子集列表。同一子集内的槽位可互换 |
| `antisymmetric_subsets` | `Vec<Vec<usize>>` | 反对称的槽位子集列表。交换子集内两个槽位翻转符号 |
| `cyclic` | `Option<Vec<usize>>` | 循环置换的槽位子集 |

### SymmetrySpec::none

**签名**：`pub fn none() -> Self`

**功能**：无对称性——所有槽位独立。

**返回值**：`SymmetrySpec`

**示例**：

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::none();
assert!(spec.symmetric_subsets.is_empty());
assert!(spec.antisymmetric_subsets.is_empty());
assert!(spec.cyclic.is_none());
// 输出：所有断言通过
```

### SymmetrySpec::fully_symmetric

**签名**：`pub fn fully_symmetric(rank: usize) -> Self`

**功能**：所有槽位完全对称。创建一个包含全部 `0..rank` 槽位的对称子集。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rank` | `usize` | 张量阶数（槽位总数） |

**返回值**：`SymmetrySpec`

**示例**：

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_symmetric(3);
assert_eq!(spec.symmetric_subsets, vec![vec![0, 1, 2]]);
// 输出：对称子集 = [[0, 1, 2]]
```

### SymmetrySpec::fully_antisymmetric

**签名**：`pub fn fully_antisymmetric(rank: usize) -> Self`

**功能**：所有槽位完全反对称。创建一个包含全部 `0..rank` 槽位的反对称子集。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rank` | `usize` | 张量阶数（槽位总数） |

**返回值**：`SymmetrySpec`

**示例**：

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_antisymmetric(4);
assert_eq!(spec.antisymmetric_subsets, vec![vec![0, 1, 2, 3]]);
// 输出：反对称子集 = [[0, 1, 2, 3]]
```

### SymmetrySpec::is_slot_hidden

**签名**：`pub fn is_slot_hidden(&self, pos: usize) -> bool`

**功能**：检查给定槽位位置是否应在图编码中"隐藏"（即不参与规范化比较）。对称槽位和循环槽位被隐藏。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pos` | `usize` | 槽位位置索引 |

**返回值**：`bool` — `true` 表示该槽位隐藏。

**示例**：

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_symmetric(2);
assert!(spec.is_slot_hidden(0));   // 在对称子集中
assert!(spec.is_slot_hidden(1));   // 在对称子集中

let none = SymmetrySpec::none();
assert!(!none.is_slot_hidden(0));  // 不在任何子集中
// 输出：所有断言通过
```

---

## TensorRegistry

**签名**：

```rust
#[derive(Debug, Clone, Default)]
pub struct TensorRegistry { /* private fields */ }
```

**功能**：张量规范化的完整注册表。记录哪些函数头是张量以及它们的对称性规格，同时管理指标维度组分配（防止跨维度的哑指标重命名冲突）。

### TensorRegistry::new

**签名**：`pub fn new() -> Self`

**功能**：创建空注册表。

**返回值**：`TensorRegistry`

### TensorRegistry::register

**签名**：`pub fn register(&mut self, name: Symbol, spec: SymmetrySpec)`

**功能**：注册一个张量名称及其对称性规格。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `Symbol` | 张量函数头名称 |
| `spec` | `SymmetrySpec` | 对称性规格 |

**示例**：

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::Symbol;

let mut reg = TensorRegistry::new();
reg.register(Symbol::new("g"), SymmetrySpec::fully_symmetric(2));
reg.register(Symbol::new("Riemann"), SymmetrySpec::none());
reg.register(Symbol::new("epsilon"), SymmetrySpec::fully_antisymmetric(4));
// 输出：注册了 3 个张量
```

### TensorRegistry::set_index_group

**签名**：`pub fn set_index_group(&mut self, label: Symbol, group: u64)`

**功能**：设置指标标签的维度组标识符。不同组的哑指标不会被重命名为相同的规范名，从而避免跨维度冲突。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `label` | `Symbol` | 指标标签（如 `"mu"`、`"i"`） |
| `group` | `u64` | 组标识符。`0` 表示默认/未分组 |

**示例**：

```rust
use ocas_atom::tensor::spec::TensorRegistry;
use ocas_atom::Symbol;

let mut reg = TensorRegistry::new();
reg.set_index_group(Symbol::new("mu"), 1);  // 时空指标
reg.set_index_group(Symbol::new("i"), 2);   // 内部指标
assert_eq!(reg.index_group(Symbol::new("mu")), 1);
assert_eq!(reg.index_group(Symbol::new("i")), 2);
assert_eq!(reg.index_group(Symbol::new("other")), 0); // 默认组
// 输出：所有断言通过
```

### TensorRegistry::spec

**签名**：`pub fn spec(&self, name: Symbol) -> Option<&SymmetrySpec>`

**功能**：查找已注册张量的对称性规格。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `Symbol` | 张量名称 |

**返回值**：`Option<&SymmetrySpec>` — 未注册时返回 `None`。

### TensorRegistry::index_group

**签名**：`pub fn index_group(&self, label: Symbol) -> u64`

**功能**：查找指标标签的维度组。未设置的标签返回 `0`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `label` | `Symbol` | 指标标签 |

**返回值**：`u64` — 组标识符。

**参见**：[SymmetrySpec](#symmetryspec)

---

## canonicalize_tensors

**签名**：

```rust
pub fn canonicalize_tensors<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<CanonicalTensor<'a>, TensorCanonError>
```

**功能**：将张量表达式规范化为与指标命名无关的唯一规范形。内部将表达式编码为图（张量头→头顶点、指标槽→槽顶点、缩并→边），运行 McKay 细化-个体化算法计算规范标号，然后重建表达式并重命名哑指标。

对于和表达式（`Add`），逐项规范化并验证所有项具有相同的自由指标集。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 待规范化的张量表达式（`Mul` of `Fun` 节点） |
| `registry` | `&TensorRegistry` | 张量注册表（对称性规格 + 指标组） |

**返回值**：`Result<CanonicalTensor<'a>, TensorCanonError>`

**错误**：

| 错误变体 | 说明 |
|---|---|
| `ContractedMoreThanOnce(Symbol)` | 某个指标标签被缩并超过两次 |
| `BadContraction(Symbol)` | 两个相同标签的槽具有相同位置（不是上/下对） |
| `NotATensor(Symbol)` | 表达式中出现未在注册表中注册的函数头 |
| `InconsistentOpenIndices` | 和表达式中各项的自由指标不一致 |
| `UnsupportedPower` | 表达式中出现不支持的幂次节点 |

**示例**：

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());

// T(i, j) 和 T(j, i) 是不同的（无对称性）
let t_ij = ctx.fun("T", &[ctx.var("i"), ctx.var("j")]);
let t_ji = ctx.fun("T", &[ctx.var("j"), ctx.var("i")]);

let ct1 = canonicalize_tensors(&ctx, t_ij, &reg).unwrap();
let ct2 = canonicalize_tensors(&ctx, t_ji, &reg).unwrap();
// 无对称性时，T(i,j) ≠ T(j,i)
// 注意：无对称性时 T(i,j) 与 T(j,i) 的规范形不同
```

**参见**：[CanonicalTensor](#canonicaltensor)、[TensorCanonError](#tensorcanonerror)、[TensorRegistry](#tensorregistry)

---

## CanonicalTensor

**签名**：

```rust
#[derive(Debug, Clone)]
pub struct CanonicalTensor<'a> {
    pub canonical_form: Atom<'a>,
    pub external_indices: Vec<Atom<'a>>,
    pub dummy_indices: Vec<Atom<'a>>,
}
```

**功能**：张量规范化的结果。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `canonical_form` | `Atom<'a>` | 规范化后的表达式（哑指标已重命名） |
| `external_indices` | `Vec<Atom<'a>>` | 自由（外部）指标列表 |
| `dummy_indices` | `Vec<Atom<'a>>` | 哑指标列表（已重命名为 `d0`、`d1`、…） |

**示例**：

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());
reg.register(Symbol::new("U"), SymmetrySpec::none());

let prod = ctx.mul(&[ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
                      ctx.fun("U", &[ctx.var("j"), ctx.var("k")])]);
let ct = canonicalize_tensors(&ctx, prod, &reg).unwrap();

println!("规范形: {}", ct.canonical_form);
println!("外部指标: {:?}", ct.external_indices);
println!("哑指标: {:?}", ct.dummy_indices);
// 输出：规范形含重命名后的哑指标，外部指标为 i, k
```

---

## TensorCanonError

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TensorCanonError {
    ContractedMoreThanOnce(Symbol),
    BadContraction(Symbol),
    NotATensor(Symbol),
    InconsistentOpenIndices,
    UnsupportedPower,
}
```

**功能**：张量规范化过程中可能产生的错误。

**变体**：

| 变体 | 说明 |
|---|---|
| `ContractedMoreThanOnce(Symbol)` | 指标标签被缩并超过两次（出现 3 次以上） |
| `BadContraction(Symbol)` | 两个相同标签的槽具有相同方差（不是有效的上/下对） |
| `NotATensor(Symbol)` | 出现未在 `TensorRegistry` 中注册的函数头 |
| `InconsistentOpenIndices` | 和表达式中各项的自由指标集不一致 |
| `UnsupportedPower` | 表达式中出现不支持的幂次节点 |

---

## refresh_dummies

**签名**：

```rust
pub fn refresh_dummies<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<Atom<'a>, DummyError>
```

**功能**：重命名张量表达式中的哑指标以避免与外部（自由）指标冲突。恰好出现两次的标签（一次上标、一次下标）被替换为按维度组分配的新名称（`d0`、`d1`、…）。外部指标保持不变。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 待处理的张量表达式 |
| `registry` | `&TensorRegistry` | 张量注册表（用于维度组分配） |

**返回值**：`Result<Atom<'a>, DummyError>`

**错误**：

| 错误变体 | 说明 |
|---|---|
| `OverContracted(Symbol)` | 某个指标标签出现超过两次 |
| `BadContraction(Symbol)` | 两个相同标签的槽具有相同方差 |

**示例**：

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::dummy::refresh_dummies;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());
reg.register(Symbol::new("U"), SymmetrySpec::none());

// T(i,j) * U(j,i)：i 和 j 都恰好出现两次 → 哑指标
let expr = ctx.mul(&[ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
                      ctx.fun("U", &[ctx.var("j"), ctx.var("i")])]);
let refreshed = refresh_dummies(&ctx, expr, &reg).unwrap();
println!("{}", refreshed);
// 输出：哑指标被重命名为 d0, d1（按维度组分配）
```

**参见**：[DummyError](#dummyerror)、[TensorRegistry](#tensorregistry)

---

## DummyError

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DummyError {
    OverContracted(Symbol),
    BadContraction(Symbol),
}
```

**功能**：哑指标操作的错误类型。

**变体**：

| 变体 | 说明 |
|---|---|
| `OverContracted(Symbol)` | 指标标签出现超过两次（不是合法的缩并候选） |
| `BadContraction(Symbol)` | 两个相同标签的槽具有相同方差（不是上/下对） |

---

## YoungTableau

**签名**：

```rust
#[derive(Debug, Clone)]
pub struct YoungTableau {
    pub row_lengths: Vec<usize>,
}
```

**功能**：Young 盘，由行长度列表定义形状（如 `[2, 1]` 表示 □□/□）。Young 投影子实现经典对称子 $c_\lambda = a_\lambda \cdot b_\lambda$：对行内置换 $\sigma \in R$ 与列内置换 $\tau \in C$ 的所有组合求和，每项带列置换奇偶符号 $\operatorname{sgn}(\tau)$。

这是**显式**展开（非 BSGS 群论实现）：结果是 $\prod r_i! \cdot \prod c_j!$ 个项的和（行、列阶乘之积），每项符号 ±1，不做 Hook 长度归一化。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `row_lengths` | `Vec<usize>` | 每行的格子数 |

### YoungTableau::new

**签名**：`pub fn new(row_lengths: Vec<usize>) -> Self`

**功能**：从行长度创建 Young 盘。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `row_lengths` | `Vec<usize>` | 每行的格子数列表。总格子数应等于张量阶数 |

**返回值**：`YoungTableau`

**示例**：

```rust
use ocas_atom::tensor::young::YoungTableau;

// [2, 1] = □□/□（3 格）
let tab = YoungTableau::new(vec![2, 1]);
assert_eq!(tab.total_boxes(), 3);

// [1, 1, 1] = 全反对称
let anti = YoungTableau::new(vec![1, 1, 1]);
assert_eq!(anti.total_boxes(), 3);

// [3] = 全对称
let sym = YoungTableau::new(vec![3]);
assert_eq!(sym.total_boxes(), 3);
// 输出：所有断言通过
```

### YoungTableau::total_boxes

**签名**：`pub fn total_boxes(&self) -> usize`

**功能**：返回盘中总格子数（等于张量阶数）。

**返回值**：`usize`

---

## young_project

**签名**：

```rust
pub fn young_project<'a>(
    ctx: &'a AtomArena<'a>,
    tensor_expr: Atom<'a>,
    tableau: &YoungTableau,
) -> Atom<'a>
```

**功能**：对张量表达式应用 Young 投影子。将 $T(i_1, i_2, \dots, i_n)$ 展开为置换和：

$$\sum_\sigma \text{sign}(\sigma) \cdot T(i_{\sigma(1)}, i_{\sigma(2)}, \dots, i_{\sigma(n)})$$

当前实现**不做** Hook 长度积归一化——每项系数仅为 $\pm 1$。对于全反对称盘 `[1, 1, …, 1]`，得到标准交替和。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `tensor_expr` | `Atom<'a>` | 张量表达式（`Fun` 节点） |
| `tableau` | `&YoungTableau` | Young 盘形状 |

**返回值**：`Atom<'a>` — 投影后的表达式（`Add` 节点，各项带符号）。

**示例**：

```rust
use ocas_atom::tensor::young::{YoungTableau, young_project};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// 全反对称 [1, 1]：f(a, b) → f(a, b) - f(b, a)
let f_ab = ctx.fun("f", &[ctx.var("a"), ctx.var("b")]);
let result = young_project(&ctx, f_ab, &YoungTableau::new(vec![1, 1]));
println!("{}", result);
// 输出：f(a,b) - f(b,a)

// 全对称 [2]：g(a, b) → g(a, b) + g(b, a)
let g_ab = ctx.fun("g", &[ctx.var("a"), ctx.var("b")]);
let result = young_project(&ctx, g_ab, &YoungTableau::new(vec![2]));
println!("{}", result);
// 输出：g(a,b) + g(b,a)

// 全反对称 [1, 1, 1]：h(a, b, c) → 6 项交错和（系数 ±1，未归一化）
let h = ctx.fun("h", &[ctx.var("a"), ctx.var("b"), ctx.var("c")]);
let result = young_project(&ctx, h, &YoungTableau::new(vec![1, 1, 1]));
println!("{}", result);
// 输出：h(a,b,c) - h(b,a,c) - ...（6 项，各项系数 ±1）
```

**参见**：[YoungTableau](#youngtableau)、[SymmetrySpec](#symmetryspec)

---

## 设计不变量

### Atom 是 Copy 的 arena 句柄

`Atom<'a>` 是指向 arena 中节点的 Copy 句柄。`IndexSlot<'a>` 也是 Copy。这意味着张量可以廉价复制和比较——结构相等等价于指针相等（hash-consing）。

### 显式指标匹配

oCAS 的张量代数使用**显式指标匹配**而非爱因斯坦求和约定。缩并必须手动调用 `contract`，指标按标签和位置（上/下）配对。

### 对称性是建议性的

`Symmetry` 枚举和 `SymmetrySpec` 是元数据，不会被 `contract` 等基础运算自动使用。完全的对称性处理需要：
1. 规范化（`canonicalize_tensors`）—— 图同构级别的对称性强制
2. Young 投影（`young_project`）—— 显式置换展开

### 图同构引擎

0.22.0 的规范化基于 McKay 细化-个体化图同构算法（实现在 `tensor::graph` 模块）。该引擎是一个独立的 nauty 实现，提供：
- 1-WL 颜色细化（邻居签名 + 边数据 + 方向）
- 公平划分 + 个体化-细化 DFS
- 路径不变量剪枝 + 自同构轨道
- 规范形 = 字典序最大证书

---

## 完整示例

```rust
use ocas_atom::tensor::{
    IndexPosition, IndexSlot, Symmetry, Tensor,
    contract, Contracted, symmetrise_sign,
};
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::tensor::dummy::refresh_dummies;
use ocas_atom::tensor::young::{YoungTableau, young_project};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    // 1. 创建带对称性的度规张量 g_μν
    let g = Tensor::new(Symbol::new("g"), vec![
        IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
        IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
    ]).with_symmetry(Symmetry::Symmetric);
    println!("g rank = {}, symmetry = {:?}", g.rank(), g.symmetry());
    // 输出：g rank = 2, symmetry = Symmetric

    // 2. 缩并 g^μν A_ν → 部分缩并
    let g_inv = Tensor::new(Symbol::new("g"), vec![
        IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
        IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
    ]);
    let a = Tensor::new(Symbol::new("A"), vec![
        IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
    ]);
    match contract(&ctx, &g_inv, &a) {
        Contracted::Product(tp) => {
            println!("缩并后自由指标数: {}", tp.factors[0].rank());
            // 输出：缩并后自由指标数: 1 (mu)
        }
        _ => {}
    }

    // 3. 规范化：T(i,j)*U(j,k) 中 j 被重命名
    let mut reg = TensorRegistry::new();
    reg.register(Symbol::new("T"), SymmetrySpec::none());
    reg.register(Symbol::new("U"), SymmetrySpec::none());

    let prod = ctx.mul(&[
        ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
        ctx.fun("U", &[ctx.var("j"), ctx.var("k")]),
    ]);
    let ct = canonicalize_tensors(&ctx, prod, &reg).unwrap();
    println!("规范形: {}", ct.canonical_form);
    // 输出：规范形含 d0, d1 等哑指标

    // 4. Young 投影：全反对称化 f(a, b)
    let f = ctx.fun("f", &[ctx.var("a"), ctx.var("b")]);
    let anti = young_project(&ctx, f, &YoungTableau::new(vec![1, 1]));
    println!("反对称化: {}", anti);
    // 输出：f(a,b) - f(b,a)
}
```

---

## 参见

- [表达式系统](./rust-expressions.md) — `Atom`、`AtomArena`、`Symbol` 基础类型
- [重写与化简](./rust-rewrite.md) — 基于模式匹配的表达式变换
- [线性代数](../math/linear-algebra.md) — 矩阵运算（张量的低阶特例）
- [张量代数与规范化](../math/tensor-canonicalization.md) — 数学理论基础
