# Rust API 参考：表达式系统

本章记录 oCAS 表达式系统的核心类型与函数，涵盖 Arena 分配器、符号、表达式节点、构造上下文、解析、规范化以及模式匹配。

---

## 核心类型概览

表达式系统的核心类型层次如下：

```
Arena (bump 分配器)
 └── AtomArena (hash-consing 构造器)
      └── Atom (Copy 句柄)
           └── AtomNode (enum: Num/Var/Fun/Add/Mul/Pow)
                └── Symbol (interned 字符串)
```

**设计不变量**：`Atom` 是一个 `Copy` 的 Arena 句柄，指向 `Arena` 中分配的 `AtomNode`。由于 `AtomArena` 使用 hash-consing（相同结构的子表达式在同一个 `AtomArena` 中返回同一个指针），**结构相等等价于指针相等**——`==` 的结果与指针比较一致。注意实现细节：`Atom` 的 `PartialEq` 是派生的结构比较（递归比较 `AtomNode`），hash-consing 保证其真值与指针相等一致，而非直接比较指针。减法和除法在解析时脱糖：$x - y$ 变为 `Add([x, Mul([Num(-1), y])])`，$x / y$ 变为 `Mul([x, Pow(y, Num(-1))])`。

---

## Arena

```rust
use ocas_core::arena::Arena;
```

`Arena` 是 bump 分配器，为表达式节点提供批量内存管理。所有表达式节点分配在 Arena 中，Arena 销毁时一次性释放全部内存。当前版本不运行析构函数，因此仅安全存储 `Copy` 类型。

### Arena::new

```rust
pub fn new() -> Self
```

**功能**：创建使用默认块大小（64 KiB）的 Arena。

**参数**：无。

**返回值**：`Arena`。

**示例**：

```rust
use ocas_core::arena::Arena;

let arena = Arena::new();
let value = arena.allocate_with(|| 42);
assert_eq!(*value, 42);
```

### Arena::with_capacity

```rust
pub fn with_capacity(block_size: usize) -> Self
```

**功能**：创建使用指定块大小的 Arena。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `block_size` | `usize` | 每个内存块的字节数。 |

**返回值**：`Arena`。

### Arena::allocate_with

```rust
pub fn allocate_with<T>(&self, init: impl FnOnce() -> T) -> &mut T
```

**功能**：在 Arena 中分配一个值，通过闭包构造，返回与 Arena 生命周期绑定的可变引用。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `init` | `impl FnOnce() -> T` | 构造值的闭包。 |

**返回值**：`&mut T`。

**错误**：若 `T` 的大小为零则 panic。

**示例**：

```rust
use ocas_core::arena::Arena;

let arena = Arena::new();
let value = arena.allocate_with(|| "hello");
assert_eq!(*value, "hello");
```

### Arena::allocate_slice

```rust
pub fn allocate_slice<T: Copy>(&self, values: &[T]) -> &[T]
```

**功能**：将一个切片复制到 Arena 中，返回与 Arena 生命周期绑定的切片引用。空切片返回 `&[]`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `values` | `&[T]` | 待分配的切片，`T: Copy`。 |

**返回值**：`&[T]`。

**示例**：

```rust
use ocas_core::arena::Arena;

let arena = Arena::new();
let data = [1, 2, 3];
let slice = arena.allocate_slice(&data);
assert_eq!(slice, &[1, 2, 3]);
```

### Arena::reset

```rust
pub fn reset(&self)
```

**功能**：重置 Arena——保留第一个内存块并重置偏移量，释放其余块。重置后，之前分配的任何引用**不得继续使用**。

**参数**：无。

**返回值**：无。

**示例**：

```rust
use ocas_core::arena::Arena;

let arena = Arena::new();
let _ = arena.allocate_with(|| 1);
arena.reset();
let value = arena.allocate_with(|| 2);
assert_eq!(*value, 2);
```

### Arena::chunk_count

```rust
pub fn chunk_count(&self) -> usize
```

**功能**：返回 Arena 当前持有的内存块数量。

**参数**：无。

**返回值**：`usize`。

---

## Symbol

```rust
use ocas_atom::Symbol;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Symbol(&'static str);
```

`Symbol` 是全局 interned 的符号名称，用于变量名、函数名和常量名。相同内容的 `Symbol` 在进程中指向同一块静态内存，比较为 $O(1)$ 指针比较。`Symbol` 实现了 `Copy`。

### Symbol::new

```rust
pub fn new(name: &str) -> Self
```

**功能**：创建或查找一个已存在的 `Symbol`。首次调用时 intern 字符串，后续调用返回已有的 `Symbol`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `&str` | 符号名称。 |

**返回值**：`Symbol`。

**示例**：

```rust
use ocas_atom::Symbol;

let x = Symbol::new("x");
let also_x = Symbol::new("x");
assert_eq!(x, also_x);
assert_eq!(x.as_str(), "x");
```

### Symbol::as_str

```rust
pub fn as_str(&self) -> &str
```

**功能**：返回符号的字符串切片。

**参数**：无。

**返回值**：`&str`。

**参见**：[Atom](#atom)、[AtomArena](#atomarena)

---

## Atom

```rust
use ocas_atom::Atom;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Atom<'a>(&'a AtomNode<'a>);
```

`Atom` 是指向 Arena 中 `AtomNode` 的轻量 `Copy` 句柄。复制 `Atom` 仅复制一个指针，不复制底层数据。由于 hash-consing，`a == b` 的真值与指针比较一致（同一 `AtomArena` 内）；其实现是派生的结构比较。

`Atom` 实现了 `Display`，输出格式为：
- `Num(42)` → `42`
- `Var("x")` → `x`
- `Fun("sin", [x])` → `sin(x)`
- `Add([x, y])` → `x + y`（子表达式若非叶子节点则加括号）
- `Mul([x, y])` → `x*y`
- `Pow(x, n)` → `x^n`（底数/指数若非叶子则加括号）

### Atom::node

```rust
pub fn node(&self) -> &'a AtomNode<'a>
```

**功能**：获取底层节点数据的引用。用于模式匹配表达式结构。

**参数**：无。

**返回值**：`&'a AtomNode<'a>`。

**示例**：

```rust
use ocas_atom::{AtomArena, AtomNode};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
assert!(matches!(x.node(), AtomNode::Var(_)));
```

### Atom::children

```rust
pub fn children(&self) -> &'a [Atom<'a>]
```

**功能**：返回直接子表达式的切片（从左到右）。`Num` 和 `Var` 返回空切片；`Fun`、`Add`、`Mul` 返回参数切片；`Pow` 返回空切片（需使用 `binary_children` 获取两个操作数）。

**参数**：无。

**返回值**：`&'a [Atom<'a>]`。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let sum = ctx.add(&[x, y, ctx.num(1)]);
assert_eq!(sum.children().len(), 3);
assert_eq!(x.children().len(), 0);
```

### Atom::binary_children

```rust
pub fn binary_children(&self) -> Option<(Atom<'a>, Atom<'a>)>
```

**功能**：若此 Atom 是 `Pow` 节点，返回 `(base, exp)`；否则返回 `None`。

**参数**：无。

**返回值**：`Option<(Atom<'a>, Atom<'a>)>`。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let power = ctx.pow(x, y);
let (base, exp) = power.binary_children().unwrap();
assert_eq!(base.to_string(), "x");
assert_eq!(exp.to_string(), "y");
```

**参见**：[Atom::children](#atomchildren)、[AtomNode](#atomnode)

---

## AtomNode

```rust
use ocas_atom::AtomNode;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AtomNode<'a> {
    Num(i64),
    Var(Symbol),
    Fun(Symbol, &'a [Atom<'a>]),
    Add(&'a [Atom<'a>]),
    Mul(&'a [Atom<'a>]),
    Pow(Atom<'a>, Atom<'a>),
}
```

表达式树中每个节点的具化数据。`Atom` 句柄通过 `node()` 方法获取 `AtomNode` 引用。

### 变体说明

| 变体 | 说明 |
|---|---|
| `Num(i64)` | 64 位有符号整数字面量。 |
| `Var(Symbol)` | 命名变量或常量。 |
| `Fun(Symbol, &'a [Atom<'a>])` | 函数应用。第一个参数是函数名，第二个是参数列表（至少一个元素）。 |
| `Add(&'a [Atom<'a>])` | 加法。参数列表至少一个元素。 |
| `Mul(&'a [Atom<'a>])` | 乘法。参数列表至少一个元素。 |
| `Pow(Atom<'a>, Atom<'a>)` | 幂运算。第一个参数是底数，第二个是指数。 |

**设计说明**：
- 减法/除法在**解析时**脱糖：$x - y$ 变为 `Add([x, Mul([Num(-1), y])])`，$x / y$ 变为 `Mul([x, Pow(y, Num(-1))])`。AST 中不存在独立的减法/除法节点。
- `Add` 和 `Mul` 的参数列表通过 `normalize` 后是排序的，确保结构相等的表达式产生相同的 AST。

**示例**：

```rust
use ocas_atom::{AtomArena, AtomNode};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
match x.node() {
    AtomNode::Var(s) => assert_eq!(s.as_str(), "x"),
    _ => panic!("expected variable"),
}
```

**参见**：[Atom](#atom)、[Symbol](#symbol)

---

## AtomArena

```rust
use ocas_atom::AtomArena;
```

```rust
pub struct AtomArena<'a> {
    arena: &'a Arena,
    cons_table: RefCell<FastHashMap<AtomNode<'a>, Atom<'a>>>,
}
```

`AtomArena` 是构造 `Atom` 的唯一入口。它封装了一个 `Arena` 引用和一个 hash-consing 表（内部可变性通过 `RefCell` 实现）。所有构造方法从调用者角度看是不可变的。相同结构的子表达式在同一个 `AtomArena` 中始终返回同一个 `Atom` 句柄——这使得结构相等等价于指针相等。

### AtomArena::new

```rust
pub fn new(arena: &'a Arena) -> Self
```

**功能**：创建以给定 Arena 为后端的 `AtomArena`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `arena` | `&'a Arena` | bump 分配器引用。 |

**返回值**：`AtomArena<'a>`。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let n = ctx.num(42);
assert_eq!(n.to_string(), "42");
```

### AtomArena::num

```rust
pub fn num(&self, value: i64) -> Atom<'a>
```

**功能**：创建整数字面量 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `value` | `i64` | 整数值。 |

**返回值**：`Atom<'a>`——`Num(value)` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let n = ctx.num(7);
assert_eq!(n.to_string(), "7");
```

### AtomArena::var

```rust
pub fn var(&self, name: &str) -> Atom<'a>
```

**功能**：创建变量 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `&str` | 变量名称。 |

**返回值**：`Atom<'a>`——`Var(Symbol)` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
assert_eq!(x.to_string(), "x");
```

### AtomArena::fun

```rust
pub fn fun(&self, name: &str, args: &[Atom<'a>]) -> Atom<'a>
```

**功能**：创建函数应用 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `&str` | 函数名。 |
| `args` | `&[Atom<'a>]` | 参数列表。debug 模式下为空会 panic。 |

**返回值**：`Atom<'a>`——`Fun(Symbol, &[Atom])` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let f = ctx.fun("sin", &[x]);
assert_eq!(f.to_string(), "sin(x)");
```

### AtomArena::add

```rust
pub fn add(&self, args: &[Atom<'a>]) -> Atom<'a>
```

**功能**：创建加法 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `args` | `&[Atom<'a>]` | 操作数列表。debug 模式下为空会 panic。 |

**返回值**：`Atom<'a>`——`Add(&[Atom])` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let sum = ctx.add(&[x, y]);
assert_eq!(sum.to_string(), "x + y");
```

### AtomArena::mul

```rust
pub fn mul(&self, args: &[Atom<'a>]) -> Atom<'a>
```

**功能**：创建乘法 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `args` | `&[Atom<'a>]` | 操作数列表。debug 模式下为空会 panic。 |

**返回值**：`Atom<'a>`——`Mul(&[Atom])` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let product = ctx.mul(&[x, y]);
assert_eq!(product.to_string(), "x*y");
```

### AtomArena::pow

```rust
pub fn pow(&self, base: Atom<'a>, exp: Atom<'a>) -> Atom<'a>
```

**功能**：创建幂运算 `Atom`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `base` | `Atom<'a>` | 底数。 |
| `exp` | `Atom<'a>` | 指数。 |

**返回值**：`Atom<'a>`——`Pow(base, exp)` 节点。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let p = ctx.pow(x, ctx.num(3));
assert_eq!(p.to_string(), "x^3");
```

### AtomArena::slice

```rust
pub fn slice(&self, atoms: &[Atom<'a>]) -> &'a [Atom<'a>]
```

**功能**：将 Atom 切片分配到 Arena 中并返回引用。用于需要共享 Arena 生命周期的多结果场景（如 ODE 系统的分量解）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atoms` | `&[Atom<'a>]` | 待分配的 Atom 切片。 |

**返回值**：`&'a [Atom<'a>]`。

**参见**：[Arena](#arena)、[Atom](#atom)

---

## 函数

### parse

```rust
use ocas_parse::parse;
```

```rust
pub fn parse<'a>(ctx: &'a AtomArena<'a>, input: &str) -> Result<Atom<'a>, ParseError>
```

**功能**：将数学表达式字符串解析为 `Atom` 树。支持整数、变量、函数调用（`f(x)`）、加减乘除幂运算和括号。减法脱糖为加法+负系数，除法脱糖为乘法+负指数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 用于分配节点的 arena 上下文。 |
| `input` | `&str` | 待解析的表达式字符串。 |

**返回值**：`Result<Atom<'a>, ParseError>`。

**错误**：

| 变体 | 说明 |
|---|---|
| `ParseError::Lex(LexError)` | 输入包含非法字符。 |
| `ParseError::UnexpectedEof` | 输入意外结束（如 `"x +"`）。 |
| `ParseError::UnexpectedToken` | 遇到意外的 token（如 `"*x"`）。 |

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_parse::parse;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = parse(&ctx, "x^2 + 2*x + 1").unwrap();
assert_eq!(expr.to_string(), "((x^2) + (2*x)) + 1");
```

**参见**：[normalize](#normalize)

### normalize

```rust
use ocas_atom::normalize::normalize;
```

```rust
pub fn normalize<'a>(ctx: &AtomArena<'a>, atom: Atom<'a>) -> Atom<'a>
```

**功能**：将表达式规范化为确定性标准形式。具体操作：
- 展平嵌套的 `Add` 和 `Mul`（如 `Add([Add([x, y]), z])` → `Add([x, y, z])`）
- 对操作数排序
- 合并数值系数（如 `Add([x, Num(2), Num(3)])` → `Add([Num(5), x])`）

结果分配在与输入相同的 arena 中。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 用于分配结果的 arena 上下文。 |
| `atom` | `Atom<'a>` | 待规范化的表达式。 |

**返回值**：`Atom<'a>`——规范化后的表达式。

**示例**：

```rust
use ocas_atom::normalize::normalize;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let z = ctx.var("z");
let inner = ctx.add(&[x, y]);
let outer = ctx.add(&[inner, z, ctx.num(2), ctx.num(3)]);
let result = normalize(&ctx, outer);
assert_eq!(result.to_string(), "5 + x + y + z");
```

**参见**：[parse](#parse)、[transform](#transform)

### substitute

```rust
use ocas_calc::series::substitute;
```

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a>
```

**功能**：将表达式中所有出现的变量 `var` 替换为 `replacement`。基于 `transform` 实现的便捷函数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | arena 上下文。 |
| `expr` | `Atom<'a>` | 待替换的表达式。 |
| `var` | `Symbol` | 被替换的变量名。 |
| `replacement` | `Atom<'a>` | 替换表达式。 |

**返回值**：`Atom<'a>`——替换后的表达式。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::series::substitute;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(1)]);
let result = substitute(&ctx, expr, Symbol::new("x"), y);
assert_eq!(result.to_string(), "(y^2) + 1");
```

**参见**：[transform](#transform)

### transform

```rust
use ocas_rewrite::transform;
```

```rust
pub fn transform<'a, F>(ctx: &'a AtomArena<'a>, atom: Atom<'a>, mut f: F) -> Atom<'a>
where
    F: FnMut(Atom<'a>) -> Option<Atom<'a>>,
```

**功能**：自底向上遍历表达式树并对每个节点应用变换函数 `f`。`f` 在子节点已被变换之后调用。返回 `Some(atom)` 替换原节点，返回 `None` 保留原节点（已变换子节点的版本）。这是 oCAS 规则引擎和化简器使用的标准重写遍历模式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | arena 上下文。 |
| `atom` | `Atom<'a>` | 待变换的表达式。 |
| `f` | `FnMut(Atom<'a>) -> Option<Atom<'a>>` | 变换函数。 |

**返回值**：`Atom<'a>`——变换后的表达式。

**示例**：

```rust
use ocas_atom::{Atom, AtomArena, AtomNode};
use ocas_core::arena::Arena;
use ocas_rewrite::transform;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let sum = ctx.add(&[x, y]);

let result = transform(&ctx, sum, |a| {
    if let AtomNode::Add(args) = a.node() {
        if args.len() == 2 && args[0] == x && args[1] == y {
            return Some(ctx.add(&[y, x]));
        }
    }
    None
});

assert_eq!(result.to_string(), "y + x");
```

**参见**：[substitute](#substitute)、[normalize](#normalize)

### collect_funs

```rust
use ocas_atom::walk::collect_funs;
```

```rust
pub fn collect_funs<'a>(atom: Atom<'a>) -> Vec<(Symbol, Atom<'a>)>
```

**功能**：收集表达式中所有函数应用，按后序（最内层优先）返回，去重（hash-consing 保证结构相同的应用只出现一次）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom<'a>` | 待遍历的表达式。 |

**返回值**：`Vec<(Symbol, Atom<'a>)>`——`(函数名, 函数应用节点)` 列表。

**示例**：

```rust
use ocas_atom::{AtomArena, walk::collect_funs};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let sin_x = ctx.fun("sin", &[x]);
let expr = ctx.fun("cos", &[sin_x]);
let funs = collect_funs(expr);
assert_eq!(funs.len(), 2);
assert_eq!(funs[0].0.as_str(), "sin");
assert_eq!(funs[1].0.as_str(), "cos");
```

**参见**：[collect_vars](#collect_vars)

### collect_vars

```rust
use ocas_atom::walk::collect_vars;
```

```rust
pub fn collect_vars(atom: Atom) -> Vec<Symbol>
```

**功能**：收集表达式中所有不同的变量名，按首次出现顺序（深度优先、从左到右）返回。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom` | 待遍历的表达式。 |

**返回值**：`Vec<Symbol>`——变量名列表。

**示例**：

```rust
use ocas_atom::{AtomArena, walk::collect_vars};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let expr = ctx.add(&[ctx.mul(&[x, y]), x]);
let vars = collect_vars(expr);
assert_eq!(vars.len(), 2);
assert_eq!(vars[0].as_str(), "x");
assert_eq!(vars[1].as_str(), "y");
```

**参见**：[collect_funs](#collect_funs)

---

## 模式匹配

模式匹配系统允许用带通配符的模式（`Pattern`）匹配 `Atom` 表达式树，支持结合律/交换律（AC）匹配和回溯预算控制。

### Pattern

```rust
use ocas_rewrite::pattern::Pattern;
```

```rust
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Pattern<'a> {
    Literal(Atom<'a>),
    Wildcard(Symbol, WildcardLevel),
    Add(Vec<Pattern<'a>>),
    Mul(Vec<Pattern<'a>>),
    Pow(Box<(Pattern<'a>, Pattern<'a>)>),
    Fun(Symbol, Vec<Pattern<'a>>),
}
```

模式镜像了 `AtomNode` 的结构，但增加通配符节点。`Add` 和 `Mul` 的匹配支持结合律-交换律（AC）：模式中的子模式可以匹配参数列表的任意子集。

| 变体 | 说明 |
|---|---|
| `Literal(Atom)` | 精确匹配给定的 Atom。 |
| `Wildcard(Symbol, WildcardLevel)` | 通配符匹配，名称和级别决定绑定行为。 |
| `Add(Vec<Pattern>)` | 匹配 `Add` 节点，参数列表 AC 匹配。 |
| `Mul(Vec<Pattern>)` | 匹配 `Mul` 节点，参数列表 AC 匹配。 |
| `Pow(Box<(Pattern, Pattern)>)` | 匹配 `Pow` 节点，底数和指数分别匹配。 |
| `Fun(Symbol, Vec<Pattern>)` | 匹配 `Fun` 节点，函数名和参数列表分别匹配。 |

#### Pattern::from_atom

```rust
pub fn from_atom(_ctx: &'a impl PatternAlloc<'a>, atom: Atom<'a>) -> Pattern<'a>
```

**功能**：将 `Atom` 转换为 `Pattern`。名称以 `_` 结尾（或以 `_` 开头）的变量被视为通配符：
- `x_` → `Wildcard(Symbol("x"), Single)`
- `x__` → `Wildcard(Symbol("x"), Sequence)`
- `x___` → `Wildcard(Symbol("x"), NullSequence)`

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `_ctx` | `&'a impl PatternAlloc<'a>` | 模式分配器（通常传 `&()`）。 |
| `atom` | `Atom<'a>` | 待转换的表达式。 |

**返回值**：`Pattern<'a>`。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x_");
let pat = Pattern::from_atom(&(), x);
assert!(matches!(pat, Pattern::Wildcard(s, WildcardLevel::Single) if s.as_str() == "x"));
```

**参见**：[WildcardLevel](#wildcardlevel)、[match_pattern](#match_pattern)

### WildcardLevel

```rust
use ocas_rewrite::pattern::WildcardLevel;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WildcardLevel {
    Single,
    Sequence,
    NullSequence,
}
```

通配符的匹配范围。

| 变体 | 命名约定 | 说明 |
|---|---|---|
| `Single` | 名称以单个 `_` 结尾（如 `x_`；以 `_` 开头亦可，如 `_x`） | 匹配恰好一个 Atom。 |
| `Sequence` | 名称以 `__` 结尾（如 `x__`；以 `__` 开头亦可，如 `__x`） | 匹配 `Add`/`Mul`/`Fun` 参数列表中的一个或多个 Atom。 |
| `NullSequence` | 名称以 `___` 结尾（如 `x___`；以 `___` 开头亦可，如 `___x`） | 匹配零个或多个 Atom。 |

**示例**：

```rust
use ocas_rewrite::pattern::WildcardLevel;

assert!(matches!(WildcardLevel::Single, WildcardLevel::Single));
assert!(matches!(WildcardLevel::Sequence, WildcardLevel::Sequence));
assert!(matches!(WildcardLevel::NullSequence, WildcardLevel::NullSequence));
```

### match_pattern

```rust
use ocas_rewrite::matcher::match_pattern;
```

```rust
pub fn match_pattern<'a>(pattern: Pattern<'a>, atom: Atom<'a>) -> Result<Bindings<'a>, MatchError>
```

**功能**：尝试将 `pattern` 与 `atom` 匹配，成功时返回绑定集合。使用默认回溯预算（`DEFAULT_MAX_BACKTRACKS = 10_000`）。`Add`/`Mul` 使用 AC 匹配（全回溯搜索），`Fun` 使用有序匹配（支持序列通配符）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 匹配模式。 |
| `atom` | `Atom<'a>` | 被匹配的表达式。 |

**返回值**：`Result<Bindings<'a>, MatchError>`。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::matcher::match_pattern;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let sum = ctx.add(&[x, ctx.num(1)]);

let pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
    Pattern::Literal(ctx.num(1)),
]);
let bindings = match_pattern(pat, sum).unwrap();
```

### match_pattern_with_budget

```rust
use ocas_rewrite::matcher::match_pattern_with_budget;
```

```rust
pub fn match_pattern_with_budget<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
    max_backtracks: usize,
) -> Result<Bindings<'a>, MatchError>
```

**功能**：与 `match_pattern` 相同，但允许自定义回溯预算，防止病态输入导致指数时间爆炸。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 匹配模式。 |
| `atom` | `Atom<'a>` | 被匹配的表达式。 |
| `max_backtracks` | `usize` | 最大回溯次数。 |

**返回值**：`Result<Bindings<'a>, MatchError>`。

**错误**：若回溯次数超限，返回 `MatchError::BudgetExhausted`。

**参见**：[match_pattern](#match_pattern)

### Bindings

```rust
use ocas_rewrite::matcher::Bindings;
```

```rust
#[derive(Debug, Clone, Default)]
pub struct Bindings<'a> { /* 私有字段 */ }
```

成功匹配后返回的通配符绑定集合。

#### Bindings::new

```rust
pub fn new() -> Self
```

**功能**：创建空的绑定集合。

#### Bindings::get

```rust
pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>>
```

**功能**：按通配符名称查找绑定值。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `Symbol` | 通配符名称（不含尾部下划线）。 |

**返回值**：`Option<&MatchValue<'a>>`。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::matcher::match_pattern;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");
let sum = ctx.add(&[x, y]);

let pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("a"), WildcardLevel::Single),
    Pattern::Wildcard(Symbol::new("b"), WildcardLevel::Single),
]);
let bindings = match_pattern(pat, sum).unwrap();
let val = bindings.get(Symbol::new("a")).unwrap();
```

### MatchValue

```rust
use ocas_rewrite::matcher::MatchValue;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchValue<'a> {
    Single(Atom<'a>),
    Sequence(&'a [Atom<'a>]),
}
```

通配符绑定的值类型。

| 变体 | 说明 |
|---|---|
| `Single(Atom<'a>)` | `Single` 级别通配符绑定到单个 Atom。 |
| `Sequence(&'a [Atom<'a>])` | `Sequence` 或 `NullSequence` 级别通配符绑定到一个 Atom 切片。 |

### MatchError

```rust
use ocas_rewrite::matcher::MatchError;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchError {
    NoMatch,
    InconsistentBinding,
    BudgetExhausted,
}
```

模式匹配错误。

| 变体 | Display | 说明 |
|---|---|---|
| `NoMatch` | `"pattern did not match"` | 模式不匹配目标表达式。 |
| `InconsistentBinding` | `"inconsistent wildcard binding"` | 同一通配符名称在不同位置绑定到不同值。 |
| `BudgetExhausted` | `"backtrack budget exhausted"` | AC 匹配的回溯预算耗尽。表示输入可能导致指数时间匹配，应增大预算或简化模式。 |

---

## 参见

- [重写与化简](./rust-rewrite.md)——`Rule`、`simplify`、`replace_all` 等完整重写引擎
- [微积分](./rust-calculus.md)——`diff`、`integrate`、`taylor` 等基于表达式系统的运算
- [系数域](./rust-domains.md)——`Domain` trait 和各种域类型
