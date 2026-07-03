# Rewrite & Simplification

oCAS provides a pattern-matching engine and a rule-based simplifier for
symbolic expressions. This chapter covers the core concepts and APIs.

---

## Patterns

A `Pattern` describes an expression structure with wildcards. Three wildcard
levels are available, triggered by the wildcard name suffix:

| Level | Name suffix | Matches |
|---|---|---|
| `WildcardLevel::Single` | `_` (e.g. `x_`) | Any single sub-expression |
| `WildcardLevel::Sequence` | `__` (e.g. `__x`) | One or more operands in an ordered list |
| `WildcardLevel::NullSequence` | `___` (e.g. `___x`) | Zero or more operands in an ordered list |

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// Pattern: x + y_  where y_ matches any single sub-expression
let x = ctx.var("x");
let pat = Pattern::Add(vec![
    Pattern::Literal(x),
    Pattern::Wildcard(Symbol::new("y_"), WildcardLevel::Single),
]);

// Match: x + 5  →  binds y_ to 5
let e = parse(&ctx, "x + 5").unwrap();
let bindings = match_pattern(pat, e).unwrap();
match bindings.get(Symbol::new("y_")).unwrap() {
    MatchValue::Single(v) => assert_eq!(v.to_string(), "5"),
    _ => {}
}
```

Within `Add` and `Mul`, matching is **associative and commutative**:
arguments are sorted and matched in a canonical order.

---

## Match & Bindings

`match_pattern` returns `Result<Bindings, MatchError>`. `Bindings` maps
wildcard `Symbol` names to matched values.

```rust
// Pattern: a_ + b_ + ___rest  — capture two terms and the rest
let pat = Pattern::Add(vec![
    Pattern::Wildcard(Symbol::new("a_"), WildcardLevel::Single),
    Pattern::Wildcard(Symbol::new("b_"), WildcardLevel::Single),
    Pattern::Wildcard(Symbol::new("___rest"), WildcardLevel::NullSequence),
]);

let e = parse(&ctx, "x + y + z + 5").unwrap();
let bindings = match_pattern(pat, e).unwrap();

// Bindings use MatchValue::Single(atom) or MatchValue::Sequence(slice)
use ocas_rewrite::MatchValue;
match bindings.get(Symbol::new("a_")).unwrap() {
    MatchValue::Single(a) => println!("a = {}", a),  // e.g. "x" or "5"
    _ => {}
}
```

`MatchError` variants:
- `NoMatch` — pattern does not match
- `InconsistentBinding` — same wildcard name bound to different values

---

## Rules

A `Rule` pairs a pattern with a replacement closure that receives the
match bindings and the arena context.

```rust
use ocas_rewrite::rules::default_rules;

// Built-in rule set
let rules = default_rules(&ctx, &());

// Custom rule: x_ + 0 → x_
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

## Simplification

`simplify()` applies a rule set repeatedly until a fixed point (or max
iterations) is reached.

```rust
let e = parse(&ctx, "x + 0 + y*0 + z*1").unwrap();
let rules = default_rules(&ctx, &());
let simplified = simplify(&ctx, e, &rules, 20);
println!("{}", simplified);  // x + z
```

The default rule set handles:
- **Identity removal**: `x + 0 → x`, `x * 1 → x`, `x * 0 → 0`
- **Constant folding**: `2 + 3 → 5`, `2 * 3 → 6`
- **Power simplifications**: `x^0 → 1`, `x^1 → x`, `0^x → 0`, `1^x → 1`
- **Arithmetic on numbers**: `2*3 + 4*5 → 26`

---

## Bottom-up transformation

`transform()` walks the expression tree bottom-up, applying a function to
each node. This is useful for custom traversals that do not fit the
pattern-matching model.

```rust
// Replace every variable "x" with "t"
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

## E-graph simplification (egg feature)

With the `egg` feature enabled, oCAS can use equality saturation for more
powerful simplifications that rule-based rewriting alone cannot achieve.

```bash
cargo build -p ocas --features egg
```

```rust
// Requires `egg` feature
#[cfg(feature = "egg")]
{
    use ocas_rewrite::egraph::egg_simplify;
    let e = parse(&ctx, "sin(x)^2 + cos(x)^2").unwrap();
    let result = egg_simplify(&ctx, e).unwrap();
    println!("{}", result);  // 1
}
```

The e-graph approach explores multiple equivalent forms simultaneously,
combining rewrites via congruence closure. This handles trigonometric
identities and algebraic equalities that require multiple rewrite steps
in a specific order.

---

## Limitations

The default rule-based simplifier intentionally keeps things predictable:
- Does **not** handle `sin(x)^2 + cos(x)^2 → 1` (requires the `egg` feature)
- Does **not** perform full polynomial expansion by default
- Does **not** apply trigonometric or logarithmic identities

For advanced simplification, enable the `egg` feature or write custom rules.

---

## See also

- [Rust API](./rust-api.md) — building expressions and patterns from Rust
- [Evaluation & JIT](./evaluation.md) — numeric evaluation after simplification
