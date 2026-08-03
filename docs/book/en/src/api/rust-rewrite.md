# Rust API Reference: Rewrite & Simplification

This chapter documents the full public API of the `ocas-rewrite` crate. Built on top of `ocas-atom`, it provides:

- **Pattern matching**: wildcards (`_`, `__`, `___`) with AC (associative-commutative) backtracking matching
- **Rewrite rules**: the `Rule` abstraction pairing a pattern with a replacement closure
- **Simplification engine**: `simplify` / `simplify_with_fuel` apply rules iteratively until a fixed point
- **Replacement**: `replace_once` / `replace_all` / `replace_all_multiple` controlled replacement
- **Transformation**: `transform` bottom-up traversal, `partition_expr` combinatorial partitioning
- **E-graphs**: equality saturation simplification based on the `egg` crate (requires the `egg` feature)

Module layout:

```
ocas_rewrite
├── pattern      — Pattern, WildcardLevel, PatternAlloc
├── matcher      — match_pattern, Bindings, MatchValue, MatchError
├── rules        — Rule, default_rules, 9 built-in rules
├── simplify     — simplify, simplify_with_fuel
├── replace      — replace_once, replace_all, replace_all_multiple, Replacement, Condition
├── transformer  — transform, partition_expr
├── combinatorics — partitions, PartitionSolution
└── egraph       — AtomLanguage, simplify_with_egraph (feature = "egg")
```

---

## Pattern System

### WildcardLevel

```rust
pub enum WildcardLevel {
    Single,       // x_  — matches a single atom
    Sequence,     // __x — matches 1 or more atoms (in an Add/Mul/Fun argument list)
    NullSequence, // ___x — matches 0 or more atoms
}
```

**Function**: Defines the matching scope of a wildcard.

| Variant | Description | Naming convention |
|---|---|---|
| `Single` | Matches exactly one atom | Variable name ends with `_`, e.g. `x_` |
| `Sequence` | Matches one or more atoms (non-empty sequence) | Variable name starts with `__`, e.g. `__x` |
| `NullSequence` | Matches zero or more atoms (nullable sequence) | Variable name ends with `___`, e.g. `x___` |

**See also**: `Pattern`, `match_pattern`

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

**Function**: A pattern AST that mirrors the `Atom` structure but adds wildcard nodes.

| Variant | Description |
|---|---|
| `Literal(Atom)` | Matches the given atom exactly |
| `Wildcard(name, level)` | Wildcard, matched per `WildcardLevel` and bound to `name` |
| `Add(pats)` | Matches an addition node; arguments matched by `pats` (AC matching) |
| `Mul(pats)` | Matches a multiplication node; arguments matched by `pats` (AC matching) |
| `Pow(Box<(base, exp)>)` | Matches an exponentiation node |
| `Fun(head, pats)` | Matches a function application; `head` is the function name, arguments matched in order |

#### Pattern::from_atom

```rust
pub fn from_atom(_ctx: &'a impl PatternAlloc<'a>, atom: Atom<'a>) -> Pattern<'a>
```

**Function**: Converts an `Atom` into a `Pattern`. Variable names ending or starting with `_`, `__`, or `___` are automatically recognized as wildcards of the corresponding level (e.g. `x_`, `__x`, `x___`).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&impl PatternAlloc<'a>` | The pattern allocator (pass `&()`) |
| `atom` | `Atom<'a>` | The expression to convert |

**Return value**: `Pattern<'a>` — the corresponding pattern.

**Example**:

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

**See also**: `WildcardLevel`, `match_pattern`

---

### PatternAlloc

```rust
pub trait PatternAlloc<'a> {
    fn alloc_slice(&self, items: &[Pattern<'a>]) -> &'a [Pattern<'a>];
}
```

**Function**: Helper trait for allocating pattern slices in a caller-owned scratch arena, avoiding leaks into the global arena. The `()` implementation simply returns the slice via `Box::leak` (fine for one-off examples).

**See also**: `Rule`, `default_rules`

---

## Pattern Matching

### MatchValue

```rust
pub enum MatchValue<'a> {
    Single(Atom<'a>),
    Sequence(&'a [Atom<'a>]),
}
```

**Function**: The value of a wildcard binding.

| Variant | Description |
|---|---|
| `Single(Atom)` | A single wildcard (`_`) bound to one atom |
| `Sequence(&[Atom])` | A sequence wildcard (`__` or `___`) bound to a slice of atoms |

**See also**: `Bindings`

---

### Bindings

```rust
pub struct Bindings<'a> { /* private fields */ }
```

**Function**: The set of wildcard bindings produced by a successful pattern match.

#### Bindings::new

```rust
pub fn new() -> Self
```

**Function**: Creates an empty binding set.

#### Bindings::get

```rust
pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>>
```

**Function**: Looks up a bound value by name.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `Symbol` | The wildcard name (without trailing underscores) |

**Return value**: `Option<&MatchValue<'a>>` — a reference to the bound value if present.

**Example**:

```rust
use ocas_rewrite::matcher::Bindings;

let bindings = Bindings::new();
// Usually produced by match_pattern
```

**See also**: `match_pattern`, `MatchValue`

---

### MatchError

```rust
pub enum MatchError {
    NoMatch,
    InconsistentBinding,
    BudgetExhausted,
}
```

**Function**: The reason a pattern match failed.

| Variant | Description |
|---|---|
| `NoMatch` | The pattern structure does not match |
| `InconsistentBinding` | The same-named wildcard is bound to different values |
| `BudgetExhausted` | AC backtracking exceeded the budget (10,000 by default) |

Implements the `Display` and `Error` traits.

**See also**: `match_pattern_with_budget`

---

### match_pattern

```rust
pub fn match_pattern<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
) -> Result<Bindings<'a>, MatchError>
```

**Function**: Matches a pattern against an atom and returns the binding set. Uses the default backtracking budget (10,000 backtracks).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The pattern to match |
| `atom` | `Atom<'a>` | The target expression |

**Return value**:
- `Ok(Bindings)` — match succeeded, containing all wildcard bindings
- `Err(MatchError)` — match failed

**Matching semantics**:
- `Add`/`Mul` nodes use AC (associative-commutative) matching: argument order does not matter, associativity is handled automatically
- `Fun` nodes match arguments in order, with sequence wildcard support
- Same-named wildcards must bind to the same value

**Example**:

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

**See also**: `match_pattern_with_budget`, `Pattern`, `Bindings`

---

### match_pattern_with_budget

```rust
pub fn match_pattern_with_budget<'a>(
    pattern: Pattern<'a>,
    atom: Atom<'a>,
    max_backtracks: usize,
) -> Result<Bindings<'a>, MatchError>
```

**Function**: Pattern matching with a custom backtracking budget. Suitable for patterns with higher complexity that need a larger search space.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The pattern to match |
| `atom` | `Atom<'a>` | The target expression |
| `max_backtracks` | `usize` | The maximum number of backtracks |

**Return value**: Same as `match_pattern`.

**Errors**: Returns `MatchError::BudgetExhausted` when the budget is exhausted.

**Example**:

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

// Use a larger budget for complex patterns
let bindings = match_pattern_with_budget(pat, atom, 100_000).unwrap();
assert!(bindings.get(Symbol::new("x")).is_some());
```

**See also**: `match_pattern`, `DEFAULT_MAX_BACKTRACKS`

---

### DEFAULT_MAX_BACKTRACKS

```rust
pub const DEFAULT_MAX_BACKTRACKS: usize = 10_000;
```

**Function**: The default backtracking budget constant used by `match_pattern`.

---

## Rewrite Rules

### Rule

```rust
pub struct Rule<'a> { /* private fields */ }
```

**Function**: A rewrite rule pairing a pattern with a replacement closure. When the pattern matches a subexpression, the replacement closure receives the bindings and produces a new atom.

#### Rule::new

```rust
pub fn new<F>(pattern: Pattern<'a>, replacement: F) -> Self
where
    F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a> + 'a,
```

**Function**: Creates a rule from a pattern and a replacement closure.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The matching pattern |
| `replacement` | `Fn(&Bindings, &AtomArena) -> Atom` | The replacement closure, receiving bindings and the arena context |

**Return value**: `Rule<'a>` — an unconditional rule.

**Example**:

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

**Function**: Adds a precondition to the rule. The condition is evaluated after a successful match, before the replacement.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `condition` | `Fn(&Bindings) -> bool` | The condition closure; the replacement runs only when it returns `true` |

**Return value**: `Rule<'a>` — a rule with a condition.

**Example**:

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

**Function**: Attempts to apply the rule to an atom. Returns `Some(replacement result)` when the pattern matches and the condition holds, otherwise `None`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The target expression |

**Return value**: `Option<Atom<'a>>` — `Some` means the rule was applied, `None` means no match or the condition failed.

**See also**: `match_pattern`, `simplify`

---

### default_rules

```rust
pub fn default_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>>
```

**Function**: Returns the built-in set of algebraic simplification rules (9 rules in total).

**Built-in rule list**:

| Rule | Rewrite | Description |
|---|---|---|
| `add_zero` | `x + 0 → x` | Additive identity |
| `add_zero_left` | `0 + x → x` | Additive identity (left) |
| `mul_zero` | `x × 0 → 0` | Multiplicative zero |
| `mul_zero_left` | `0 × x → 0` | Multiplicative zero (left) |
| `mul_one` | `x × 1 → x` | Multiplicative identity |
| `mul_one_left` | `1 × x → x` | Multiplicative identity (left) |
| `add_same` | `x + x → 2x` | Combine like terms |
| `pow_zero` | `x⁰ → 1` | Zero-th power |
| `pow_one` | `x¹ → x` | First power |

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `alloc` | `&impl PatternAlloc<'a>` | The pattern allocator (pass `&()`) |

**Return value**: `Vec<Rule<'a>>` — the rule list.

**Example**:

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::default_rules;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let rules = default_rules(&ctx, &());
// rules.len() == 9
```

**See also**: `simplify`, `Rule`

---

### Built-in rule functions

Each of the following functions returns a single built-in rule and can be used individually or combined into custom rule sets:

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

Each rule has the same signature and usage: pass `ctx` and `alloc`, get back a `Rule<'a>`.

**Example** (using `add_same` alone):

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

## Simplification Engine

### simplify

```rust
pub fn simplify<'a>(
    ctx: &'a AtomArena<'a>,
    atom: Atom<'a>,
    rules: &[Rule<'a>],
    iter_limit: usize,
) -> Atom<'a>
```

**Function**: Applies the rewrite rules to the expression repeatedly, traversing bottom-up, until a fixed point or the iteration limit is reached. This is the core entry point of the oCAS simplification engine.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The expression to simplify |
| `rules` | `&[Rule<'a>]` | The list of rewrite rules |
| `iter_limit` | `usize` | The maximum number of iterations (10–20 recommended) |

**Return value**: `Atom<'a>` — the simplified expression (registered in the arena).

**Behavior**:
1. Each iteration traverses the expression tree bottom-up
2. All rules are tried in turn on each node
3. If any rule fires, the round produced a change and the next round runs
4. Stops when no rule fires or when `iter_limit` is reached

**Example**:

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

**See also**: `simplify_with_fuel`, `default_rules`, `Rule`

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

**Function**: Simplification with a fuel budget. Each bottom-up traversal round consumes one fuel unit; when the fuel runs out, the operation stops early.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The expression to simplify |
| `rules` | `&[Rule<'a>]` | The list of rewrite rules |
| `iter_limit` | `usize` | The maximum number of iterations |
| `fuel` | `&Fuel` | The fuel budget (`Fuel::default()` is effectively unlimited) |

**Return value**:
- `Ok(Atom)` — simplification complete (fixed point reached)
- `Err(OcasError::OutOfFuel)` — fuel exhausted (no expression is returned in this case; callers must save intermediate results themselves)

**Purpose**: Prevents pathological inputs from making the simplification engine loop forever. Participates in fuel accounting within nested call chains.

**Example**:

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

**See also**: `simplify`, `Fuel` (from the `ocas-core` crate)

---

## Replacement

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

**Function**: Traverses top-down and replaces the first matching subexpression, then stops.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The target expression |
| `pattern` | `Pattern<'a>` | The matching pattern |
| `replacement` | `Fn(&Bindings, &AtomArena) -> Atom` | The replacement closure |

**Return value**: `Atom<'a>` — the replaced expression (returned unchanged if nothing matches).

**Example**:

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
// The first matching subexpression is replaced by 2*a
```

**See also**: `replace_all`, `replace_all_multiple`

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

**Function**: Traverses top-down and replaces all matching subexpressions. Does not perform nested replacement (i.e., replacement results do not participate in further matching).

**Parameters**: Same as `replace_once`.

**Return value**: `Atom<'a>` — the expression with all matches replaced.

**See also**: `replace_once`, `replace_all_multiple`

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

**Function**: Uses multiple replacement rules, tried in order; the first matching rule takes effect. Traverses all nodes top-down.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The target expression |
| `replacements` | `&[Replacement<'a, F>]` | The list of replacement rules (ordered by priority) |

**Return value**: `Atom<'a>` — the replaced expression.

**See also**: `Replacement`, `replace_all`

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

**Function**: A single replacement rule, containing a pattern, a replacement closure, and an optional condition.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The matching pattern |
| `replacement` | `F` | The replacement closure |
| `condition` | `Option<Condition<'a>>` | The optional precondition |

**See also**: `replace_all_multiple`, `Condition`

---

### Condition

```rust
pub enum Condition<'a> {
    Predicate(Arc<dyn Fn(&Bindings<'a>) -> bool + 'a>),
}
```

**Function**: A replacement condition; currently only a predicate form is supported.

#### Condition::new

```rust
pub fn new<F: Fn(&Bindings<'a>) -> bool + 'a>(f: F) -> Self
```

**Function**: Creates a condition from a closure.

**Example**:

```rust
use ocas_rewrite::replace::Condition;

let cond = Condition::new(|bindings| {
    // Replace only when the bound variable x is positive
    true
});
```

**See also**: `Replacement`, `Rule::with_condition`

---

### ReplaceSettings

```rust
pub struct ReplaceSettings {
    pub once: bool,
    pub bottom_up: bool,
    pub nested: bool,
}
```

**Function**: Replacement traversal settings.

| Field | Type | Default | Description |
|---|---|---|---|
| `once` | `bool` | `false` | Whether to replace only the first match |
| `bottom_up` | `bool` | `false` | Whether to traverse bottom-up (default is top-down) |
| `nested` | `bool` | `false` | Whether to continue replacing inside already-replaced subexpressions |

**See also**: `replace_all`

---

## Transformation

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

**Function**: Traverses the expression tree bottom-up and calls the closure `f` on each node. When the closure returns `Some(atom)`, the node is replaced; when it returns `None`, the original node is kept (with transformed subexpressions).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `atom` | `Atom<'a>` | The target expression |
| `f` | `FnMut(Atom) -> Option<Atom>` | The transformation closure |

**Return value**: `Atom<'a>` — the transformed expression.

**Behavior**:
1. Recursively traverses the subexpressions (leaf nodes are returned directly)
2. Rebuilds the current node with the transformed subexpressions
3. Calls `f` on the rebuilt node
4. If `f` returns `Some`, the new node is used; if `None`, the rebuilt node is used

**Example**:

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

**See also**: `simplify`, `partition_expr`

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

**Function**: Partitions an `arg(a₁, a₂, …, aₙ)` expression into named bins, returning the sum of products $\sum \text{coeff} \cdot f_1(\cdots) \cdot f_2(\cdots) \cdots$. The parameter design mirrors Symbolica's `Transformer::Partition`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The expression arena context |
| `expr` | `Atom<'a>` | The input expression (must be of the form `arg(...)`) |
| `bins` | `&[(Symbol, usize)]` | The bin list, each element `(function name, capacity)` |
| `fill_last` | `bool` | Whether surplus elements are absorbed into the last bin |
| `repeat` | `bool` | Whether the bin pattern repeats until all elements are consumed |

**Return value**: `Atom<'a>` — the sum of partition results; `ctx.num(0)` when no valid partition exists; the input `expr` is returned unchanged when it is not an `arg(...)` application or an argument is not numeric.

**Behavior**:
- The input must be an `arg(...)` function application with numeric arguments
- Enumerates all legal assignment ways
- The coefficient of each assignment is the multinomial coefficient $\binom{n}{k_1, k_2, \dots}$

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_rewrite::transformer::partition_expr;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let args = ctx.fun("arg", &[ctx.num(1), ctx.num(2), ctx.num(3)]);
let bins = &[(Symbol::new("f"), 2), (Symbol::new("g"), 1)];
let result = partition_expr(&ctx, args, bins, false, false);
// Returns Σ coeff · f(…, …) · g(…)
```

**See also**: `transform`, `combinatorics::partitions`

---

## Combinatorics Helpers

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

**Function**: Enumerates all ways to partition `elements` into `bins`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `elements` | `&[T]` | The elements to distribute |
| `bins` | `&[(B, usize)]` | The bin list `(name, capacity)` |
| `fill_last` | `bool` | Whether surplus elements are absorbed into the last bin |
| `repeat` | `bool` | Whether the bin pattern repeats until the elements are exhausted |

**Return value**: `Vec<PartitionSolution<T, B>>` — all legal partitions.

**See also**: `partition_expr`, `PartitionSolution`

---

### PartitionSolution

```rust
pub struct PartitionSolution<T, B> {
    pub coefficient: usize,
    pub bins: Vec<(B, Vec<T>)>,
}
```

**Function**: A single partition solution.

| Field | Type | Description |
|---|---|---|
| `coefficient` | `usize` | The multinomial coefficient |
| `bins` | `Vec<(B, Vec<T>)>` | The name and contents of each bin |

---

## E-Graph Simplification

> **Feature Gate**: requires the `egg` feature. Add to `Cargo.toml`:
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

**Function**: Simplifies an expression using the `egg` equality saturation engine. Converts the oCAS expression into an `egg` E-graph, runs the built-in rule set, and extracts the optimal expression using the AST size as the cost function.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom<'a>` | The expression to simplify |
| `ocas_arena` | `&AtomArena<'a>` | The oCAS expression arena |
| `iter_limit` | `usize` | The equality saturation iteration limit |

**Return value**: `Atom<'a>` — the simplified expression.

**Built-in E-graph rules**:

| Rule name | Rewrite |
|---|---|
| `add-zero` | `(add 0 ?a) → ?a` |
| `mul-zero` | `(mul ?a 0) → 0` |
| `mul-one` | `(mul 1 ?a) → ?a` |
| `pow-zero` | `(pow ?a 0) → 1` |
| `pow-one` | `(pow ?a 1) → ?a` |
| `pythagorean` | $\sin^2 x + \cos^2 x \to 1$ |

**Example**:

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

**See also**: `AtomLanguage` (`ocas-rewrite::egraph`)

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

**Function**: The oCAS implementation of the `egg::Language` trait, used for E-graph node representation.

**Note**: The `Language` trait is implemented manually because `egg`'s `define_language!` macro does not support `i64` and `Symbol` as leaf types.

#### AtomLanguage::to_recexpr

```rust
pub fn to_recexpr<'b>(
    atom: Atom<'b>,
    egraph: &mut egg::EGraph<Self, ()>,
    cache: &mut Vec<(Atom<'b>, Id)>,
) -> Id
```

**Function**: Converts an oCAS `Atom` into an `egg` `RecExpr`. Shared subexpressions are deduplicated via `cache`.

#### AtomLanguage::from_recexpr

```rust
pub fn from_recexpr<'a>(
    expr: &RecExpr<Self>,
    id: Id,
    ocas_arena: &'a AtomArena<'a>,
) -> Atom<'a>
```

**Function**: Converts an `egg` `RecExpr` back into an oCAS `Atom`.

**See also**: `simplify_with_egraph`, `egg::Language`

---

## Module Dependencies

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
