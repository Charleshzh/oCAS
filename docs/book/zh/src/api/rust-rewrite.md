# Rust API 参考：重写与化简

本章记录 `ocas-rewrite` crate 的全部公共 API。该 crate 构建于 `ocas-atom` 之上，提供：

- **模式匹配**：支持通配符（`_`、`__`、`___`）与 AC（交换-结合）回溯匹配
- **重写规则**：`Rule` 抽象将模式与替换闭包配对
- **化简引擎**：`simplify` / `simplify_with_fuel` 迭代应用规则直到不动点
- **替换**：`replace_once` / `replace_all` / `replace_all_multiple` 受控替换
- **变换**：`transform` 自底向上遍历，`partition_expr` 组合分区
- **E-图**：基于 `egg` crate 的等式饱和化简（需 `egg` feature）

模块层级：

```
ocas_rewrite
├── pattern      — Pattern, WildcardLevel, PatternAlloc
├── matcher      — match_pattern, Bindings, MatchValue, MatchError
├── rules        — Rule, default_rules, 9 条内置规则
├── simplify     — simplify, simplify_with_fuel
├── replace      — replace_once, replace_all, replace_all_multiple, Replacement, Condition
├── transformer  — transform, partition_expr
├── combinatorics — partitions, PartitionSolution
└── egraph       — AtomLanguage, simplify_with_egraph (feature = "egg")
```

---

## 模式系统

### WildcardLevel

```rust
pub enum WildcardLevel {
    Single,       // x_  — 匹配单个原子
    Sequence,     // __x — 匹配 1 个或多个原子（在 Add/Mul/Fun 参数列表中）
    NullSequence, // ___x — 匹配 0 个或多个原子
}
```

**功能**：定义通配符的匹配范围。

| 变体 | 说明 | 命名约定 |
|---|---|---|
| `Single` | 恰好匹配一个原子 | 变量名以 `_` 结尾，如 `x_` |
| `Sequence` | 匹配一个或多个原子（非空序列） | 变量名以 `__` 开头，如 `__x` |
| `NullSequence` | 匹配零个或多个原子（可空序列） | 变量名以 `___` 结尾，如 `x___` |

**参见**：`Pattern`、`match_pattern`

---

### Pattern

```rust
pub enum Pattern<'a> {
    Literal(Atom<'a>),
    Wildcard(Symbol, WildcardLevel),
    Add(Vec<Pattern<'a>>),
    Mul(Vec<Pattern<'a>>),
    Pow(Box<(Pattern<'a>, Pattern<'a>)>),
    Fun(Symbol, Vec<Pattern<'a>>),
}
```

**功能**：模式 AST，镜像 `Atom` 结构但增加通配符节点。

| 变体 | 说明 |
|---|---|
| `Literal(Atom)` | 精确匹配给定原子 |
| `Wildcard(name, level)` | 通配符，按 `WildcardLevel` 匹配并绑定到 `name` |
| `Add(pats)` | 匹配加法节点，参数由 `pats` 匹配（AC 匹配） |
| `Mul(pats)` | 匹配乘法节点，参数由 `pats` 匹配（AC 匹配） |
| `Pow(Box<(base, exp)>)` | 匹配幂节点 |
| `Fun(head, pats)` | 匹配函数应用，`head` 为函数名，参数按序匹配 |

#### Pattern::from_atom

```rust
pub fn from_atom(_ctx: &'a impl PatternAlloc<'a>, atom: Atom<'a>) -> Pattern<'a>
```

**功能**：将 `Atom` 转换为 `Pattern`。变量名以 `_`、`__`、`___` 结尾或开头的自动识别为对应级别的通配符（如 `x_`、`__x`、`x___`）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&impl PatternAlloc<'a>` | 模式分配器（传 `&()` 即可） |
| `atom` | `Atom<'a>` | 待转换的表达式 |

**返回值**：`Pattern<'a>` — 对应的模式。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let atom = ctx.var("x_");
let pat = Pattern::from_atom(&(), atom);
// pat == Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single)
```

**参见**：`WildcardLevel`、`match_pattern`

---

### PatternAlloc

```rust
pub trait PatternAlloc<'a> {
    fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>];
}
```

**功能**：模式切片分配辅助 trait，用于在调用方持有的暂存 arena 中分配模式切片，避免泄漏到全局 arena。`()` 实现直接通过 `Box::leak` 返回切片（适用于一次性示例）。

**参见**：`Rule`、`default_rules`

---

## 模式匹配

### MatchValue

```rust
pub enum MatchValue<'a> {
    Single(Atom<'a>),
    Sequence(&'a [Atom<'a>]),
}
```

**功能**：通配符绑定值。

| 变体 | 说明 |
|---|---|
| `Single(Atom)` | 单个通配符（`_`）绑定到一个原子 |
| `Sequence(&[Atom])` | 序列通配符（`__` 或 `___`）绑定到原子切片 |

**参见**：`Bindings`

---

### Bindings

```rust
pub struct Bindings<'a> { /* 私有字段 */ }
```

**功能**：模式匹配成功后产生的通配符绑定集合。

#### Bindings::new

```rust
pub fn new() -> Self
```

**功能**：创建空绑定集合。

#### Bindings::get

```rust
pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>>
```

**功能**：按名称查询绑定值。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `Symbol` | 通配符名称（不含尾部下划线） |

**返回值**：`Option<&MatchValue<'a>>` — 存在则返回绑定值引用。

**示例**：

```rust
use ocas_rewrite::matcher::Bindings;

let bindings = Bindings::new();
// 通常由 match_pattern 产生
```

**参见**：`match_pattern`、`MatchValue`

---

### MatchError

```rust
pub enum MatchError {
    NoMatch,
    InconsistentBinding,
    BudgetExhausted,
}
```

**功能**：模式匹配失败原因。

| 变体 | 说明 |
|---|---|
| `NoMatch` | 模式结构不匹配 |
| `InconsistentBinding` | 同名通配符绑定到不同值 |
| `BudgetExhausted` | AC 回溯次数超出预算（默认 10,000 次） |

实现 `Display` 和 `Error` trait。

**参见**：`match_pattern_with_budget`

---

### match_pattern

```rust
pub fn match_pattern<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
) -> Result<Bindings<'a>, MatchError>
```

**功能**：将模式与原子进行匹配，返回绑定集合。使用默认回溯预算（10,000 次）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 待匹配的模式 |
| `atom` | `Atom<'a>` | 目标表达式 |

**返回值**：
- `Ok(Bindings)` — 匹配成功，包含所有通配符绑定
- `Err(MatchError)` — 匹配失败

**匹配语义**：
- `Add`/`Mul` 节点使用 AC（交换-结合）匹配：参数顺序无关，自动处理结合性
- `Fun` 节点按序匹配参数，支持序列通配符
- 同名通配符必须绑定到相同值

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_atom::Symbol;
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::matcher::{match_pattern, MatchValue};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

let pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single),
    Pattern::Literal(ctx.num(0)),
]);

let expr = ctx.add(&[ctx.var("y"), ctx.num(0)]);
let bindings = match_pattern(pat, expr).unwrap();
let x = bindings.get(Symbol::new("x")).unwrap();
// MatchValue::Single(Atom("y"))
```

**参见**：`match_pattern_with_budget`、`Pattern`、`Bindings`

---

### match_pattern_with_budget

```rust
pub fn match_pattern_with_budget<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
    max_backtracks: usize,
) -> Result<Bindings<'a>, MatchError>
```

**功能**：带自定义回溯预算的模式匹配。适用于模式复杂度较高、需要更大搜索空间的场景。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 待匹配的模式 |
| `atom` | `Atom<'a>` | 目标表达式 |
| `max_backtracks` | `usize` | 最大回溯次数 |

**返回值**：同 `match_pattern`。

**错误**：预算耗尽时返回 `MatchError::BudgetExhausted`。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_atom::Symbol;
use ocas_core::arena::Arena;
use ocas_rewrite::matcher::match_pattern_with_budget;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
let atom = ctx.var("y");

// 使用更大预算处理复杂模式
let bindings = match_pattern_with_budget(pat, atom, 100_000).unwrap();
assert!(bindings.get(Symbol::new("x")).is_some());
```

**参见**：`match_pattern`、`DEFAULT_MAX_BACKTRACKS`

---

### DEFAULT_MAX_BACKTRACKS

```rust
pub const DEFAULT_MAX_BACKTRACKS: usize = 10_000;
```

**功能**：`match_pattern` 的默认回溯预算常量。

---

## 重写规则

### Rule

```rust
pub struct Rule<'a> { /* 私有字段 */ }
```

**功能**：重写规则，将模式与替换闭包配对。当模式匹配子表达式时，替换闭包接收绑定并产生新原子。

#### Rule::new

```rust
pub fn new<F>(pattern: Pattern<'a>, replacement: F) -> Self
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a> + 'a,
```

**功能**：从模式和替换闭包创建规则。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 匹配模式 |
| `replacement` | `Fn(&Bindings, &AtomArena) -> Atom` | 替换闭包，接收绑定和 arena 上下文 |

**返回值**：`Rule<'a>` — 无条件规则。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_atom::Symbol;
use ocas_core::arena::Arena;
use ocas_rewrite::matcher::{Bindings, MatchValue};
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::rules::Rule;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
let rule = Rule::new(pat, |bindings: &Bindings, _ctx: &AtomArena| {
    let x = bindings.get(Symbol::new("x")).unwrap();
    let MatchValue::Single(v) = x else { panic!("expected single"); };
    _ctx.mul(&[_ctx.num(2), *v])
});

let y = ctx.var("y");
let result = rule.apply(&ctx, y).unwrap();
assert_eq!(result.to_string(), "2*y");
```

#### Rule::with_condition

```rust
pub fn with_condition<F>(self, condition: F) -> Self
where
    F: Fn(&Bindings<'a>) -> bool + 'a,
```

**功能**：为规则添加前置条件。条件在匹配成功后、替换之前求值。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `condition` | `Fn(&Bindings) -> bool` | 条件闭包，返回 `true` 时才执行替换 |

**返回值**：`Rule<'a>` — 带条件的规则。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_atom::Symbol;
use ocas_core::arena::Arena;
use ocas_rewrite::matcher::{Bindings, MatchValue};
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::rules::Rule;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
let rule = Rule::new(pat, |_bindings: &Bindings, ctx: &AtomArena| {
    ctx.num(99)
}).with_condition(|bindings: &Bindings| {
    let x = bindings.get(Symbol::new("x")).unwrap();
    let MatchValue::Single(v) = x else { return false; };
    v.to_string() == "y"
});

let y = ctx.var("y");
let z = ctx.var("z");
assert_eq!(rule.apply(&ctx, y).unwrap().to_string(), "99");
assert!(rule.apply(&ctx, z).is_none());
```

#### Rule::apply

```rust
pub fn apply(&self, ctx: &AtomArena<'a>, atom: Atom<'a>) -> Option<Atom<'a>>
```

**功能**：尝试将规则应用于原子。匹配成功且条件满足时返回 `Some(替换结果)`，否则返回 `None`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 目标表达式 |

**返回值**：`Option<Atom<'a>>` — `Some` 表示规则已应用，`None` 表示未匹配或条件不满足。

**参见**：`match_pattern`、`simplify`

---

### default_rules

```rust
pub fn default_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>>
```

**功能**：返回内置代数化简规则集（共 9 条）。

**内置规则列表**：

| 规则 | 化简 | 说明 |
|---|---|---|
| `add_zero` | `x + 0 → x` | 加法零元 |
| `add_zero_left` | `0 + x → x` | 加法零元（左侧） |
| `mul_zero` | `x × 0 → 0` | 乘法零元 |
| `mul_zero_left` | `0 × x → 0` | 乘法零元（左侧） |
| `mul_one` | `x × 1 → x` | 乘法单位元 |
| `mul_one_left` | `1 × x → x` | 乘法单位元（左侧） |
| `add_same` | `x + x → 2x` | 合并同类项 |
| `pow_zero` | `x⁰ → 1` | 零次幂 |
| `pow_one` | `x¹ → x` | 一次幂 |

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `alloc` | `&impl PatternAlloc<'a>` | 模式分配器（传 `&()` 即可） |

**返回值**：`Vec<Rule<'a>>` — 规则列表。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::default_rules;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let rules = default_rules(&ctx, &());
// rules.len() == 9
```

**参见**：`simplify`、`Rule`

---

### 内置规则函数

以下函数各自返回一条内置规则，可单独使用或组合自定义规则集：

```rust
pub fn add_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn add_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn mul_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn mul_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn mul_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn mul_one_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn add_same<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn pow_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
pub fn pow_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a>
```

每条规则的签名和用法相同：传入 `ctx` 和 `alloc`，返回 `Rule<'a>`。

**示例**（单独使用 `add_same`）：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::add_same;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let expr = ctx.add(&[x, x]);
let rule = add_same(&ctx, &());
let result = rule.apply(&ctx, expr).unwrap();
assert_eq!(result.to_string(), "2*x");
```

---

## 化简引擎

### simplify

```rust
pub fn simplify<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    rules: &[Rule<'a>],
    iter_limit: usize,
) -> Atom<'a>
```

**功能**：对表达式反复应用重写规则，自底向上遍历，直到达到不动点或迭代上限。这是 oCAS 化简引擎的核心入口。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 待化简的表达式 |
| `rules` | `&[Rule<'a>]` | 重写规则列表 |
| `iter_limit` | `usize` | 最大迭代次数（推荐 10–20） |

**返回值**：`Atom<'a>` — 化简后的表达式（已注册到 arena）。

**行为**：
1. 每轮迭代自底向上遍历表达式树
2. 对每个节点依次尝试所有规则
3. 任何规则触发则该轮有变化，继续下一轮
4. 无规则触发或达到 `iter_limit` 则停止

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let rules = default_rules(&ctx, &());

let x = ctx.var("x");
let expr = ctx.mul(&[x, ctx.num(0)]);
let result = simplify(&ctx, expr, &rules, 10);
assert_eq!(result.to_string(), "0");
```

**参见**：`simplify_with_fuel`、`default_rules`、`Rule`

---

### simplify_with_fuel

```rust
pub fn simplify_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    rules: &[Rule<'a>],
    iter_limit: usize,
    fuel: &Fuel,
) -> Result<Atom<'a>, OcasError>
```

**功能**：带燃料预算的化简。每轮自底向上遍历消耗一个燃料单位，燃料耗尽时提前终止。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 待化简的表达式 |
| `rules` | `&[Rule<'a>]` | 重写规则列表 |
| `iter_limit` | `usize` | 最大迭代次数 |
| `fuel` | `&Fuel` | 燃料预算（`Fuel::default()` 为有效无限） |

**返回值**：
- `Ok(Atom)` — 化简完成（不动点达成）
- `Err(OcasError::OutOfFuel)` — 燃料耗尽（此时不返回表达式，调用方需自行保存中间结果）

**用途**：防止病态输入导致化简引擎无限循环。在嵌套调用链中参与燃料记账。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_core::fuel::Fuel;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify_with_fuel;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.var("x");
let rules = default_rules(&ctx, &());
let fuel = Fuel::new(100);
let result = simplify_with_fuel(&ctx, expr, &rules, 20, &fuel);
match result {
    Ok(e) => println!("simplified: {}", e),
    Err(_) => println!("fuel exhausted"),
}
```

**参见**：`simplify`、`Fuel`（`ocas-core` crate）

---

## 替换

### replace_once

```rust
pub fn replace_once<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    pattern: Pattern<'a>,
    replacement: F,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
```

**功能**：自顶向下遍历，替换第一个匹配的子表达式后停止。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 目标表达式 |
| `pattern` | `Pattern<'a>` | 匹配模式 |
| `replacement` | `Fn(&Bindings, &AtomArena) -> Atom` | 替换闭包 |

**返回值**：`Atom<'a>` — 替换后的表达式（若无匹配则原样返回）。

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_atom::Symbol;
use ocas_core::arena::Arena;
use ocas_rewrite::pattern::{Pattern, WildcardLevel};
use ocas_rewrite::replace::replace_once;
use ocas_rewrite::matcher::MatchValue;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
let expr = ctx.add(&[ctx.var("a"), ctx.var("b")]);

let result = replace_once(&ctx, expr, pat, |bindings, ctx| {
    let x = bindings.get(Symbol::new("x")).unwrap();
    let MatchValue::Single(v) = x else { panic!() };
    ctx.mul(&[ctx.num(2), *v])
});
// 第一个匹配的子表达式被替换为 2*a
```

**参见**：`replace_all`、`replace_all_multiple`

---

### replace_all

```rust
pub fn replace_all<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    pattern: Pattern<'a>,
    replacement: F,
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
```

**功能**：自顶向下遍历，替换所有匹配的子表达式。不进行嵌套替换（即替换结果不再参与后续匹配）。

**参数**：同 `replace_once`。

**返回值**：`Atom<'a>` — 所有匹配替换后的表达式。

**参见**：`replace_once`、`replace_all_multiple`

---

### replace_all_multiple

```rust
pub fn replace_all_multiple<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    replacements: &[Replacement<'a, F>],
) -> Atom<'a>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
```

**功能**：使用多条替换规则，按顺序尝试，第一个匹配的规则生效。自顶向下遍历所有节点。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 目标表达式 |
| `replacements` | `&[Replacement<'a, F>]` | 替换规则列表（按优先级排序） |

**返回值**：`Atom<'a>` — 替换后的表达式。

**参见**：`Replacement`、`replace_all`

---

### Replacement

```rust
pub struct Replacement<'a, F>
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a>,
{
    pub pattern: Pattern<'a>,
    pub replacement: F,
    pub condition: Option<Condition<'a>>,
}
```

**功能**：单条替换规则，包含模式、替换闭包和可选条件。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `pattern` | `Pattern<'a>` | 匹配模式 |
| `replacement` | `F` | 替换闭包 |
| `condition` | `Option<Condition<'a>>` | 可选前置条件 |

**参见**：`replace_all_multiple`、`Condition`

---

### Condition

```rust
pub enum Condition<'a> {
    Predicate(Arc<dyn Fn(&Bindings<'a>) -> bool + 'a>),
}
```

**功能**：替换条件，当前仅支持谓词形式。

#### Condition::new

```rust
pub fn new<F: Fn(&Bindings<'a>) -> bool + 'a>(f: F) -> Self
```

**功能**：从闭包创建条件。

**示例**：

```rust
use ocas_rewrite::replace::Condition;

let cond = Condition::new(|bindings| {
    // 仅当绑定变量 x 为正数时替换
    true
});
```

**参见**：`Replacement`、`Rule::with_condition`

---

### ReplaceSettings

```rust
pub struct ReplaceSettings {
    pub once: bool,
    pub bottom_up: bool,
    pub nested: bool,
}
```

**功能**：替换遍历设置。

| 字段 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `once` | `bool` | `false` | 是否只替换第一个匹配 |
| `bottom_up` | `bool` | `false` | 是否自底向上遍历（默认自顶向下） |
| `nested` | `bool` | `false` | 是否在已替换的子表达式中继续替换 |

**参见**：`replace_all`

---

## 变换

### transform

```rust
pub fn transform<'a, F>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    f: F,
) -> Atom<'a>
where
    F: FnMut(Atom<'a>) -> Option<Atom<'a>>,
```

**功能**：自底向上遍历表达式树，对每个节点调用闭包 `f`。闭包返回 `Some(atom)` 时替换该节点，返回 `None` 时保留原节点（子表达式已变换）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `atom` | `Atom<'a>` | 目标表达式 |
| `f` | `FnMut(Atom) -> Option<Atom>` | 变换闭包 |

**返回值**：`Atom<'a>` — 变换后的表达式。

**行为**：
1. 递归遍历子表达式（叶子节点直接返回）
2. 用变换后的子表达式重建当前节点
3. 对重建后的节点调用 `f`
4. `f` 返回 `Some` 则使用新节点，`None` 则使用重建节点

**示例**：

```rust
use ocas_atom::{Atom, AtomArena, AtomNode};
use ocas_core::arena::Arena;
use ocas_rewrite::transformer::transform;

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

**参见**：`simplify`、`partition_expr`

---

### partition_expr

```rust
pub fn partition_expr<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    bins: &[(Symbol, usize)],
    fill_last: bool,
    repeat: bool,
) -> Atom<'a>
```

**功能**：将 `arg(a₁, a₂, …, aₙ)` 表达式分入命名桶，返回乘积之和 $\sum \text{coeff} \cdot f_1(\cdots) \cdot f_2(\cdots) \cdots$。参数设计镜像 Symbolica 的 `Transformer::Partition`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena 上下文 |
| `expr` | `Atom<'a>` | 输入表达式（必须是 `arg(...)` 形式） |
| `bins` | `&[(Symbol, usize)]` | 桶列表，每个元素为 `(函数名, 容量)` |
| `fill_last` | `bool` | 多余元素是否吸收到最后一个桶 |
| `repeat` | `bool` | 桶模式是否重复直到所有元素消耗完毕 |

**返回值**：`Atom<'a>` — 分区结果之和；无有效分区时返回 `ctx.num(0)`；输入不是 `arg(...)` 函数应用或参数含非数值时原样返回 `expr`。

**行为**：
- 输入必须是 `arg(...)` 函数应用，参数为数值
- 枚举所有合法分配方式
- 每种分配的系数为多项式系数 $\binom{n}{k_1, k_2, \dots}$

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_rewrite::transformer::partition_expr;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let args = ctx.fun("arg", &[ctx.num(1), ctx.num(2), ctx.num(3)]);
let bins = &[(Symbol::new("f"), 2), (Symbol::new("g"), 1)];
let result = partition_expr(&ctx, args, bins, false, false);
// 返回 Σ coeff · f(…, …) · g(…)
```

**参见**：`transform`、`combinatorics::partitions`

---

## 组合数学辅助

### partitions

```rust
pub fn partitions<T, B>(
    elements: &[T],
    bins: &[(B, usize)],
    fill_last: bool,
    repeat: bool,
) -> Vec<PartitionSolution<T, B>>
where
    T: Clone + Ord + Hash,
    B: Clone + Ord + Hash,
```

**功能**：枚举将 `elements` 分入 `bins` 的所有方式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `elements` | `&[T]` | 待分配元素 |
| `bins` | `&[(B, usize)]` | 桶列表 `(名称, 容量)` |
| `fill_last` | `bool` | 多余元素吸收到最后一个桶 |
| `repeat` | `bool` | 桶模式重复直到元素耗尽 |

**返回值**：`Vec<PartitionSolution<T, B>>` — 所有合法分区。

**参见**：`partition_expr`、`PartitionSolution`

---

### PartitionSolution

```rust
pub struct PartitionSolution<T, B> {
    pub coefficient: usize,
    pub bins: Vec<(B, Vec<T>)>,
}
```

**功能**：单个分区解。

| 字段 | 类型 | 说明 |
|---|---|---|
| `coefficient` | `usize` | 多项式系数 |
| `bins` | `Vec<(B, Vec<T>)>` | 每个桶的名称和内容 |

---

## E-图化简

> **Feature Gate**：需要启用 `egg` feature。在 `Cargo.toml` 中添加：
> ```toml
> [dependencies]
> ocas-rewrite = { version = "0.23", features = ["egg"] }
> ```

### simplify_with_egraph

```rust
pub fn simplify_with_egraph<'a>(
    atom: Atom<'a>,
    ocas_arena: &'a AtomArena<'a>,
    iter_limit: usize,
) -> Atom<'a>
```

**功能**：使用 `egg` 等式饱和引擎化简表达式。将 oCAS 表达式转换为 `egg` E-图，运行内置规则集，使用 AST 大小作为代价函数提取最优表达式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `atom` | `Atom<'a>` | 待化简的表达式 |
| `ocas_arena` | `&AtomArena<'a>` | oCAS 表达式 arena |
| `iter_limit` | `usize` | 等式饱和迭代上限 |

**返回值**：`Atom<'a>` — 化简后的表达式。

**内置 E-图规则**：

| 规则名 | 化简 |
|---|---|
| `add-zero` | `(add 0 ?a) → ?a` |
| `mul-zero` | `(mul ?a 0) → 0` |
| `mul-one` | `(mul 1 ?a) → ?a` |
| `pow-zero` | `(pow ?a 0) → 1` |
| `pow-one` | `(pow ?a 1) → ?a` |
| `pythagorean` | $\sin^2 x + \cos^2 x \to 1$ |

**示例**：

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::egraph::simplify_with_egraph;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

let x = ctx.var("x");
let two = ctx.num(2);
let sin_x = ctx.fun("sin", &[x]);
let cos_x = ctx.fun("cos", &[x]);
let sum = ctx.add(&[ctx.pow(sin_x, two), ctx.pow(cos_x, two)]);

let result = simplify_with_egraph(sum, &ctx, 5);
assert_eq!(result.to_string(), "1");
```

**参见**：`AtomLanguage`（`ocas-rewrite::egraph`）

---

### AtomLanguage

```rust
pub enum AtomLanguage {
    Num(i64),
    Var(Symbol),
    Fun(Vec<Id>),
    Add(Vec<Id>),
    Mul(Vec<Id>),
    Pow([Id; 2]),
}
```

**功能**：`egg::Language` trait 的 oCAS 实现，用于 E-图节点表示。

**说明**：手动实现 `Language` trait，因为 `egg` 的 `define_language!` 宏不支持 `i64` 和 `Symbol` 作为叶子类型。

#### AtomLanguage::to_recexpr

```rust
pub fn to_recexpr<'b>(
    atom: Atom<'b>,
    egraph: &mut egg::EGraph<Self, ()>,
    cache: &mut Vec<(Atom<'b>, Id)>,
) -> Id
```

**功能**：将 oCAS `Atom` 转换为 `egg` `RecExpr`。共享子表达式通过 `cache` 去重。

#### AtomLanguage::from_recexpr

```rust
pub fn from_recexpr<'a>(
    expr: &RecExpr<Self>,
    id: Id,
    ocas_arena: &'a AtomArena<'a>,
) -> Atom<'a>
```

**功能**：将 `egg` `RecExpr` 转换回 oCAS `Atom`。

**参见**：`simplify_with_egraph`、`egg::Language`

---

## 模块依赖关系

```mermaid
graph TD
    A[ocas-atom] --> B[ocas-rewrite]
    C[ocas-core] --> B
    B --> D[pattern]
    B --> E[matcher]
    B --> F[rules]
    B --> G[simplify]
    B --> H[replace]
    B --> I[transformer]
    B --> J[combinatorics]
    B --> K[egraph]
    D --> E
    D --> F
    E --> F
    E --> G
    E --> H
    F --> G
    I --> G
    J --> I
    K -.->|feature = "egg"| B
```
