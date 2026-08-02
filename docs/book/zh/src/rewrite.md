# 重写与化简

oCAS 提供模式匹配引擎和基于规则的表达简化器。本章介绍核心概念与 API。

---

## 模式

`Pattern` 用通配符描述表达式结构。通配符名称后缀决定匹配级别：

| 级别 | 名称后缀 | 匹配 |
|---|---|---|
| `WildcardLevel::Single` | `_`（如 `x_`） | 任意单个子表达式 |
| `WildcardLevel::Sequence` | `__`（如 `__x`） | 有序列表中的一个或多个操作数 |
| `WildcardLevel::NullSequence` | `___`（如 `___x`） | 有序列表中的零个或多个操作数 |

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// 模式：x + y_  其中 y_ 匹配任意单个子表达式
let x = ctx.var("x");
let pat = Pattern::Add(vec![
    Pattern::Literal(x),
    Pattern::Wildcard(Symbol::new("y_"), WildcardLevel::Single),
]);

// 匹配：x + 5  →  将 y_ 绑定到 5
let e = parse(&ctx, "x + 5").unwrap();
let bindings = match_pattern(pat, e).unwrap();
match bindings.get(Symbol::new("y_")).unwrap() {
    MatchValue::Single(v) => assert_eq!(v.to_string(), "5"),
    _ => {}
}
```

在 `Add` 和 `Mul` 内部，匹配是**结合可交换**的：参数被排序并按规范顺序匹配。

---

## 匹配与绑定

`match_pattern` 返回 `Result<Bindings, MatchError>`。`Bindings` 将通配符 `Symbol` 名称映射到匹配值。

```rust
// 模式：a_ + b_ + ___rest — 捕获两项及剩余部分
let pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("a_"), WildcardLevel::Single),
    Pattern::Wildcard(Symbol::new("b_"), WildcardLevel::Single),
    Pattern::Wildcard(Symbol::new("___rest"), WildcardLevel::NullSequence),
]);

let e = parse(&ctx, "x + y + z + 5").unwrap();
let bindings = match_pattern(pat, e).unwrap();

// 绑定使用 MatchValue::Single(atom) 或 MatchValue::Sequence(slice)
use ocas_rewrite::MatchValue;
match bindings.get(Symbol::new("a_")).unwrap() {
    MatchValue::Single(a) => println!("a = {}", a),  // 例如 "x" 或 "5"
    _ => {}
}
```

`MatchError` 变体：
- `NoMatch` — 模式不匹配
- `InconsistentBinding` — 同一通配符名称绑定到不同值

---

## 规则

`Rule` 将模式与替换闭包配对，闭包接收匹配绑定和 arena 上下文。

```rust
use ocas_rewrite::rules::default_rules;

// 内置规则集
let rules = default_rules(&ctx, &());

// 自定义规则：x_ + 0 → x_
let custom_pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("x_"), WildcardLevel::Single),
    Pattern::Literal(ctx.num(0)),
]);
let custom_rule = Rule::new(custom_pat, |bindings, _ctx| {
    match bindings.get(Symbol::new("x_")).unwrap() {
        MatchValue::Single(x) => *x,
        _ => unreachable!(),
    }
});
```

---

## 化简

`simplify()` 反复应用规则集，直到不动点（或达到最大迭代次数）。

```rust
let e = parse(&ctx, "x + 0 + y*0 + z*1").unwrap();
let rules = default_rules(&ctx, &());
let simplified = simplify(&ctx, e, &rules, 20);
println!("{}", simplified);  // x + z
```

默认规则集处理：
- **恒等元移除**：`x + 0 → x`、`x * 1 → x`、`x * 0 → 0`
- **常量折叠**：`2 + 3 → 5`、`2 * 3 → 6`
- **幂化简**：`x^0 → 1`、`x^1 → x`、`0^x → 0`、`1^x → 1`
- **数字运算**：`2*3 + 4*5 → 26`

---

## 自底向上变换

`transform()` 自底向上遍历表达式树，对每个节点应用函数。适用于不适合模式匹配模型的自定义遍历。

```rust
// 将每个变量 "x" 替换为 "t"
let replacer = |_ctx: &AtomArena, atom: Atom| {
    if let AtomNode::Var(sym) = _ctx.get(atom) {
        if sym.as_str() == "x" {
            return _ctx.var("t");
        }
    }
    atom
};

let e = parse(&ctx, "x^2 + x + 1").unwrap();
let result = transform(&ctx, e, &replacer);
println!("{}", result);  // t^2 + t + 1
```

---

## E-graph 化简（egg feature）

启用 `egg` feature 后，oCAS 可使用等式饱和进行更强大的化简，这是纯规则重写无法实现的。

```bash
cargo build -p ocas --features egg
```

```rust
// 需要 `egg` feature
#[cfg(feature = "egg")]
{
    use ocas_rewrite::egraph::egg_simplify;
    let e = parse(&ctx, "sin(x)^2 + cos(x)^2").unwrap();
    let result = egg_simplify(&ctx, e).unwrap();
    println!("{}", result);  // 1
}
```

E-graph 方法同时探索多种等价形式，通过同余闭包组合重写。可处理需要特定多步
重写顺序的三角恒等式和代数恒等式。

---

## Fuel 受限化简

`simplify_with_fuel` 是 `simplify` 的变体，每次自底向上遍历消耗一个
fuel 单元，预算耗尽时提前停止。适用于处理不可信或病态表达式，防止重写
器无限循环。

```rust
use ocas_core::fuel::Fuel;
use ocas_rewrite::simplify::simplify_with_fuel;

let fuel = Fuel::new(100);
let result = simplify_with_fuel(&ctx, e, &rules, 20, &fuel);
match result {
    Ok(expr) => println!("simplified: {}", expr),
    Err(_) => println!("fuel exhausted before fixpoint"),
}
```

`simplify_with_fuel` 仅在 fuel 耗尽未达不动点时返回 `Err`。旧的
`simplify` API 保持可用，语义不变（无 fuel 限制）。

---

## 局限性

默认基于规则的化简器有意保持可预测性：
- **不**处理 `sin(x)^2 + cos(x)^2 → 1`（需 `egg` feature）
- 默认**不**执行完全多项式展开
- **不**应用三角恒等式或对数恒等式

如需高级化简，请启用 `egg` feature 或编写自定义规则。

---

## 参见

- [Rust API](./rust-api.md) — 从 Rust 构建表达式与模式
- [求值与 JIT](./evaluation.md) — 化简后的数值求值

---

## 回溯 AC 匹配（0.22.0）

匹配器现在对 `Add`/`Mul` 上下文使用**全回溯搜索**，取代之前的贪心
算法。这带来：

- **Add/Mul 内序列通配符**：`x_ + __rest` 匹配任意和，至少捕获一项。
- **AC 内多通配符**：`x_ + y_ + z_` 对 `a + b + c + d` 正确绑定。
- **回溯预算**：防止病态模式的组合爆炸（默认 10,000 次尝试）。

```rust
use ocas_rewrite::matcher::match_pattern_with_budget;
let bindings = match_pattern_with_budget(pattern, atom, 50_000)?;
```

---

## 多模式替换（0.22.0）

`ocas_rewrite::replace` 提供带条件守卫和遍历设置的受控替换。

### `replace_once` / `replace_all`

```rust
use ocas_rewrite::replace::{replace_once, replace_all};

// 替换树中第一个 x 出现：
let result = replace_once(&ctx, expr, Pattern::Literal(x), |_, ctx| ctx.num(42));
// 替换所有 x：
let result = replace_all(&ctx, expr, Pattern::Literal(x), |_, ctx| ctx.num(42));
```

### `replace_all_multiple` — 首匹配获胜

```rust
use ocas_rewrite::replace::{replace_all_multiple, Replacement};
let replacements = vec![
    Replacement { pattern: pat1, replacement: rhs1, condition: None },
    Replacement { pattern: pat2, replacement: rhs2, condition: Some(cond) },
];
let result = replace_all_multiple(&ctx, expr, &replacements);
```

### 条件

`Condition` 是对 `Bindings` 的谓词：

```rust
use ocas_rewrite::replace::Condition;
let cond = Condition::new(|bindings| {
    match bindings.get(Symbol::new("x_")) {
        Some(MatchValue::Single(a)) => matches!(a.node(), AtomNode::Num(_)),
        _ => false,
    }
});
```

---

## Transformer::Partition（0.22.0）

将 `arg(a, b, c, …)` 表达式按指定容量的命名桶分拆，返回所有合法
分配方式的求和。

```rust
use ocas_rewrite::transformer::partition_expr;
let result = partition_expr(&ctx, expr, &[(Symbol::new("f"), 2), (Symbol::new("g"), 2), (Symbol::new("f"), 1)], false, false);
```

三种模式：
- **exact**：桶总容量必须等于元素数量。
- **fill_last**：剩余元素吸收到最后一个桶。
- **repeat**：桶模式重复直到所有元素消耗完毕。
