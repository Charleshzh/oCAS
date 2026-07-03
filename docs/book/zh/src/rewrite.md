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
