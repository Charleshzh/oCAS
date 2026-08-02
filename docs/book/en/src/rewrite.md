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

## Fuel-bounded simplification

`simplify_with_fuel` is a variant of `simplify` that consumes one fuel unit
per bottom-up traversal pass and stops early when the budget is exhausted.
This is useful when processing untrusted or pathological expressions that
could cause the rewriter to spin indefinitely.

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

`simplify_with_fuel` returns `Err` only when fuel ran out before a fixpoint
was reached. The old `simplify` API remains available with identical semantics
(no fuel limit).

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

---

## Backtracking AC Matching (0.22.0)

The matcher now uses **full backtracking search** for `Add`/`Mul` contexts,
replacing the previous greedy algorithm.  This enables:

- **Sequence wildcards inside Add/Mul**: `x_ + __rest` matches any sum
  with at least one captured term.
- **Multiple wildcards in AC**: `x_ + y_ + z_` against `a + b + c + d`
  correctly binds each wildcard to a distinct term.
- **Backtrack budget**: prevents pathological blow-up on adversarial
  patterns (default 10,000 backtrack attempts).

The default `match_pattern` uses the budget; to adjust:

```rust
use ocas_rewrite::matcher::match_pattern_with_budget;
let bindings = match_pattern_with_budget(pattern, atom, 50_000)?;
```

---

## Multi-pattern Replacement (0.22.0)

`ocas_rewrite::replace` provides controlled replacements with condition
guards and traversal settings.

### `replace_once` / `replace_all`

```rust
use ocas_rewrite::replace::{replace_once, replace_all};

// Replace the first occurrence of x anywhere in the tree:
let result = replace_once(&ctx, expr, Pattern::Literal(x), |_, ctx| ctx.num(42));

// Replace all occurrences of x:
let result = replace_all(&ctx, expr, Pattern::Literal(x), |_, ctx| ctx.num(42));
```

### `replace_all_multiple` — first-match-wins

```rust
use ocas_rewrite::replace::{replace_all_multiple, Replacement};

let replacements = vec![
    Replacement { pattern: pat1, replacement: rhs1, condition: None },
    Replacement { pattern: pat2, replacement: rhs2, condition: Some(cond) },
];
let result = replace_all_multiple(&ctx, expr, &replacements);
```

At each node, replacements are tried in order; the first match wins.

### Conditions

A `Condition` is a predicate on `Bindings`:

```rust
use ocas_rewrite::replace::Condition;

let cond = Condition::new(|bindings| {
    // Only apply when bound atom is a number
    match bindings.get(Symbol::new("x_")) {
        Some(MatchValue::Single(a)) => matches!(a.node(), AtomNode::Num(_)),
        _ => false,
    }
});
```

---

## Transformer::Partition (0.22.0)

Partitions a `arg(a, b, c, …)` expression into named bins of given
capacities, returning a sum over all valid ways to distribute the elements.

```rust
use ocas_rewrite::transformer::partition_expr;
use ocas_atom::Symbol;

// arg(1, 3, 2, 3, 1) with bins [(f,2), (g,2), (f,1)]
let result = partition_expr(&ctx, expr, &[(Symbol::new("f"), 2), (Symbol::new("g"), 2), (Symbol::new("f"), 1)], false, false);
// Returns Σ coeff · f(…)*g(…)*f(…) for all valid partitions.
```

The three modes:
- **exact** (both flags false): total bin capacity must equal element count.
- **fill_last**: surplus elements absorbed into the last bin.
- **repeat**: bin pattern repeated until all elements consumed.
