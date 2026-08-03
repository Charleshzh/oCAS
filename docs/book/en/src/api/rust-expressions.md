# Rust API Reference: Expression System

This chapter documents the core types and functions of the oCAS expression system: the Arena allocator, symbols, expression nodes, construction contexts, parsing, normalization, and pattern matching.

---

## Core Type Overview

The core type hierarchy of the expression system is:

```
Arena (bump allocator)
 └── AtomArena (hash-consing constructor)
      └── Atom (Copy handle)
           └── AtomNode (enum: Num/Var/Fun/Add/Mul/Pow)
                └── Symbol (interned string)
```

**Design invariants**: `Atom` is a `Copy` arena handle pointing to an `AtomNode` allocated in the `Arena`. Because `AtomArena` uses hash-consing (structurally identical subexpressions return the same pointer within the same `AtomArena`), **structural equality is equivalent to pointer equality** — the result of `==` agrees with pointer comparison. Note the implementation detail: `Atom`'s `PartialEq` is a derived structural comparison (recursively comparing `AtomNode`s); hash-consing guarantees its truth value agrees with pointer equality rather than comparing pointers directly. Subtraction and division are desugared at parse time: $x - y$ becomes `Add([x, Mul([Num(-1), y])])`, and $x / y$ becomes `Mul([x, Pow(y, Num(-1))])`.

---

## Arena

```rust
use ocas_core::arena::Arena;
```

`Arena` is a bump allocator that provides bulk memory management for expression nodes. All expression nodes are allocated in an `Arena`, and the entire memory is freed at once when the `Arena` is dropped. The current version does not run destructors, so it is only safe to store `Copy` types.

### Arena::new

```rust
pub fn new() -> Self
```

**Function**: Creates an `Arena` using the default block size (64 KiB).

**Parameters**: None.

**Return value**: `Arena`.

**Example**:

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

**Function**: Creates an `Arena` using the given block size.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `block_size` | `usize` | Size of each memory block, in bytes. |

**Return value**: `Arena`.

### Arena::allocate_with

```rust
pub fn allocate_with<T>(&self, init: impl FnOnce() -> T) -> &mut T
```

**Function**: Allocates a value in the `Arena`, constructed via the closure, and returns a mutable reference bound to the `Arena`'s lifetime.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `init` | `impl FnOnce() -> T` | Closure that constructs the value. |

**Return value**: `&mut T`.

**Errors**: Panics if `T` has zero size.

**Example**:

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

**Function**: Copies a slice into the `Arena` and returns a slice reference bound to the `Arena`'s lifetime. An empty slice returns `&[]`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `values` | `&[T]` | Slice to allocate, `T: Copy`. |

**Return value**: `&[T]`.

**Example**:

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

**Function**: Resets the `Arena` — keeps the first memory block and resets the offset, freeing the remaining blocks. After a reset, any references allocated previously **must not be used**.

**Parameters**: None.

**Return value**: None.

**Example**:

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

**Function**: Returns the number of memory blocks currently held by the `Arena`.

**Parameters**: None.

**Return value**: `usize`.

---

## Symbol

```rust
use ocas_atom::Symbol;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Symbol(&'static str);
```

`Symbol` is a globally interned symbol name used for variable names, function names, and constant names. `Symbol`s with identical contents point to the same static memory within the process, making comparisons $O(1)$ pointer comparisons. `Symbol` implements `Copy`.

### Symbol::new

```rust
pub fn new(name: &str) -> Self
```

**Function**: Creates or looks up an existing `Symbol`. The string is interned on first call; subsequent calls return the existing `Symbol`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `&str` | The symbol name. |

**Return value**: `Symbol`.

**Example**:

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

**Function**: Returns the string slice of the symbol.

**Parameters**: None.

**Return value**: `&str`.

**See also**: [Atom](#atom), [AtomArena](#atomarena)

---

## Atom

```rust
use ocas_atom::Atom;
```

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Atom<'a>(&'a AtomNode<'a>);
```

`Atom` is a lightweight `Copy` handle to an `AtomNode` in the `Arena`. Copying an `Atom` only copies a pointer, not the underlying data. Thanks to hash-consing, the truth value of `a == b` agrees with pointer comparison (within the same `AtomArena`); the implementation is a derived structural comparison.

`Atom` implements `Display`, with the following output formats:
- `Num(42)` → `42`
- `Var("x")` → `x`
- `Fun("sin", [x])` → `sin(x)`
- `Add([x, y])` → `x + y` (subexpressions are parenthesized unless they are leaf nodes)
- `Mul([x, y])` → `x*y`
- `Pow(x, n)` → `x^n` (base/exponent parenthesized unless they are leaves)

### Atom::node

```rust
pub fn node(&self) -> &'a AtomNode<'a>
```

**Function**: Returns a reference to the underlying node data. Used for pattern-matching the expression structure.

**Parameters**: None.

**Return value**: `&'a AtomNode<'a>`.

**Example**:

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

**Function**: Returns the slice of direct child expressions (left to right). `Num` and `Var` return an empty slice; `Fun`, `Add`, and `Mul` return their argument slices; `Pow` returns an empty slice (use `binary_children` to get the two operands).

**Parameters**: None.

**Return value**: `&'a [Atom<'a>]`.

**Example**:

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

**Function**: If this `Atom` is a `Pow` node, returns `(base, exp)`; otherwise returns `None`.

**Parameters**: None.

**Return value**: `Option<(Atom<'a>, Atom<'a>)>`.

**Example**:

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

**See also**: [Atom::children](#atomchildren), [AtomNode](#atomnode)

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

The concrete data of each node in the expression tree. An `Atom` handle obtains an `AtomNode` reference via the `node()` method.

### Variant descriptions

| Variant | Description |
|---|---|
| `Num(i64)` | 64-bit signed integer literal. |
| `Var(Symbol)` | Named variable or constant. |
| `Fun(Symbol, &'a [Atom<'a>])` | Function application. The first argument is the function name; the second is the argument list (at least one element). |
| `Add(&'a [Atom<'a>])` | Addition. The argument list has at least one element. |
| `Mul(&'a [Atom<'a>])` | Multiplication. The argument list has at least one element. |
| `Pow(Atom<'a>, Atom<'a>)` | Exponentiation. The first argument is the base; the second is the exponent. |

**Design notes**:
- Subtraction/division are desugared **at parse time**: $x - y$ becomes `Add([x, Mul([Num(-1), y])])`, and $x / y$ becomes `Mul([x, Pow(y, Num(-1))])`. There are no standalone subtraction or division nodes in the AST.
- The argument lists of `Add` and `Mul` are sorted after `normalize`, ensuring that structurally equal expressions produce the same AST.

**Example**:

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

**See also**: [Atom](#atom), [Symbol](#symbol)

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

`AtomArena` is the only entry point for constructing `Atom`s. It wraps an `Arena` reference and a hash-consing table (interior mutability via `RefCell`). All construction methods are immutable from the caller's perspective. Structurally identical subexpressions always return the same `Atom` handle within the same `AtomArena` — this makes structural equality equivalent to pointer equality.

### AtomArena::new

```rust
pub fn new(arena: &'a Arena) -> Self
```

**Function**: Creates an `AtomArena` backed by the given `Arena`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `arena` | `&'a Arena` | Reference to the bump allocator. |

**Return value**: `AtomArena<'a>`.

**Example**:

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

**Function**: Creates an integer literal `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `value` | `i64` | The integer value. |

**Return value**: `Atom<'a>` — a `Num(value)` node.

**Example**:

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

**Function**: Creates a variable `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `&str` | The variable name. |

**Return value**: `Atom<'a>` — a `Var(Symbol)` node.

**Example**:

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

**Function**: Creates a function-application `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `&str` | The function name. |
| `args` | `&[Atom<'a>]` | Argument list. Panics on an empty list in debug mode. |

**Return value**: `Atom<'a>` — a `Fun(Symbol, &[Atom])` node.

**Example**:

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

**Function**: Creates an addition `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `args` | `&[Atom<'a>]` | Operand list. Panics on an empty list in debug mode. |

**Return value**: `Atom<'a>` — an `Add(&[Atom])` node.

**Example**:

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

**Function**: Creates a multiplication `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `args` | `&[Atom<'a>]` | Operand list. Panics on an empty list in debug mode. |

**Return value**: `Atom<'a>` — a `Mul(&[Atom])` node.

**Example**:

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

**Function**: Creates an exponentiation `Atom`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `base` | `Atom<'a>` | The base. |
| `exp` | `Atom<'a>` | The exponent. |

**Return value**: `Atom<'a>` — a `Pow(base, exp)` node.

**Example**:

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

**Function**: Allocates a slice of `Atom`s into the `Arena` and returns the reference. Used for multi-result scenarios that share the `Arena` lifetime (e.g., component solutions of ODE systems).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atoms` | `&[Atom<'a>]` | Slice of `Atom`s to allocate. |

**Return value**: `&'a [Atom<'a>]`.

**See also**: [Arena](#arena), [Atom](#atom)

---

## Functions

### parse

```rust
use ocas_parse::parse;
```

```rust
pub fn parse<'a>(ctx: &'a AtomArena<'a>, input: &str) -> Result<Atom<'a>, ParseError>
```

**Function**: Parses a mathematical expression string into an `Atom` tree. Supports integers, variables, function calls (`f(x)`), addition, subtraction, multiplication, division, exponentiation, and parentheses. Subtraction desugars to addition with a negative coefficient; division desugars to multiplication with a negative exponent.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The arena context used to allocate nodes. |
| `input` | `&str` | The expression string to parse. |

**Return value**: `Result<Atom<'a>, ParseError>`.

**Errors**:

| Variant | Description |
|---|---|
| `ParseError::Lex(LexError)` | The input contains an illegal character. |
| `ParseError::UnexpectedEof` | The input ends unexpectedly (e.g., `"x +"`). |
| `ParseError::UnexpectedToken` | An unexpected token is encountered (e.g., `"*x"`). |

**Example**:

```rust
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_parse::parse;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = parse(&ctx, "x^2 + 2*x + 1").unwrap();
assert_eq!(expr.to_string(), "((x^2) + (2*x)) + 1");
```

**See also**: [normalize](#normalize)

### normalize

```rust
use ocas_atom::normalize::normalize;
```

```rust
pub fn normalize<'a>(ctx: &AtomArena<'a>, atom: Atom<'a>) -> Atom<'a>
```

**Function**: Normalizes an expression into a deterministic canonical form. Specifically:
- Flattens nested `Add` and `Mul` (e.g., `Add([Add([x, y]), z])` → `Add([x, y, z])`)
- Sorts the operands
- Combines numeric coefficients (e.g., `Add([x, Num(2), Num(3)])` → `Add([Num(5), x])`)

The result is allocated in the same arena as the input.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | The arena context used to allocate the result. |
| `atom` | `Atom<'a>` | The expression to normalize. |

**Return value**: `Atom<'a>` — the normalized expression.

**Example**:

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

**See also**: [parse](#parse), [transform](#transform)

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

**Function**: Replaces every occurrence of the variable `var` in the expression with `replacement`. A convenience function implemented on top of `transform`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The arena context. |
| `expr` | `Atom<'a>` | The expression to substitute into. |
| `var` | `Symbol` | The variable name to replace. |
| `replacement` | `Atom<'a>` | The replacement expression. |

**Return value**: `Atom<'a>` — the substituted expression.

**Example**:

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

**See also**: [transform](#transform)

### transform

```rust
use ocas_rewrite::transform;
```

```rust
pub fn transform<'a, F>(ctx: &'a AtomArena<'a>, atom: Atom<'a>, mut f: F) -> Atom<'a>
where
    F: FnMut(Atom<'a>) -> Option<Atom<'a>>,
```

**Function**: Traverses the expression tree bottom-up and applies the transform function `f` to each node. `f` is called after the node's children have been transformed. Returning `Some(atom)` replaces the node; returning `None` keeps the node (with its transformed children). This is the standard rewrite traversal pattern used by the oCAS rule engine and simplifier.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The arena context. |
| `atom` | `Atom<'a>` | The expression to transform. |
| `f` | `FnMut(Atom<'a>) -> Option<Atom<'a>>` | The transform function. |

**Return value**: `Atom<'a>` — the transformed expression.

**Example**:

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

**See also**: [substitute](#substitute), [normalize](#normalize)

### collect_funs

```rust
use ocas_atom::walk::collect_funs;
```

```rust
pub fn collect_funs<'a>(atom: Atom<'a>) -> Vec<(Symbol, Atom<'a>)>
```

**Function**: Collects all function applications in the expression in post-order (innermost first), deduplicated (hash-consing guarantees that structurally identical applications appear only once).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom<'a>` | The expression to traverse. |

**Return value**: `Vec<(Symbol, Atom<'a>)>` — a list of `(function name, function application node)`.

**Example**:

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

**See also**: [collect_vars](#collect_vars)

### collect_vars

```rust
use ocas_atom::walk::collect_vars;
```

```rust
pub fn collect_vars(atom: Atom) -> Vec<Symbol>
```

**Function**: Collects all distinct variable names in the expression, in order of first appearance (depth-first, left to right).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `atom` | `Atom` | The expression to traverse. |

**Return value**: `Vec<Symbol>` — the list of variable names.

**Example**:

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

**See also**: [collect_funs](#collect_funs)

---

## Pattern Matching

The pattern-matching system matches `Atom` expression trees against patterns with wildcards (`Pattern`), with support for associative/commutative (AC) matching and backtracking-budget control.

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

Patterns mirror the structure of `AtomNode` but add wildcard nodes. Matching of `Add` and `Mul` supports associative-commutative (AC) matching: subpatterns can match any subset of the argument list.

| Variant | Description |
|---|---|
| `Literal(Atom)` | Matches the given `Atom` exactly. |
| `Wildcard(Symbol, WildcardLevel)` | Wildcard match; the name and level determine the binding behavior. |
| `Add(Vec<Pattern>)` | Matches an `Add` node; the argument list is matched AC-wise. |
| `Mul(Vec<Pattern>)` | Matches a `Mul` node; the argument list is matched AC-wise. |
| `Pow(Box<(Pattern, Pattern)>)` | Matches a `Pow` node; base and exponent are matched separately. |
| `Fun(Symbol, Vec<Pattern>)` | Matches a `Fun` node; function name and argument list are matched separately. |

#### Pattern::from_atom

```rust
pub fn from_atom(_ctx: &'a impl PatternAlloc<'a>, atom: Atom<'a>) -> Pattern<'a>
```

**Function**: Converts an `Atom` into a `Pattern`. Variables whose names end in `_` (or begin with `_`) are treated as wildcards:
- `x_` → `Wildcard(Symbol("x"), Single)`
- `x__` → `Wildcard(Symbol("x"), Sequence)`
- `x___` → `Wildcard(Symbol("x"), NullSequence)`

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `_ctx` | `&'a impl PatternAlloc<'a>` | The pattern allocator (usually `&()`). |
| `atom` | `Atom<'a>` | The expression to convert. |

**Return value**: `Pattern<'a>`.

**Example**:

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

**See also**: [WildcardLevel](#wildcardlevel), [match_pattern](#match_pattern)

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

The matching scope of a wildcard.

| Variant | Naming convention | Description |
|---|---|---|
| `Single` | One trailing `_` in the name (e.g., `x_`; a leading `_` also works, e.g. `_x`) | Matches exactly one `Atom`. |
| `Sequence` | Two trailing `_` in the name (e.g., `x__`; two leading `_` also work, e.g. `__x`) | Matches one or more `Atom`s in an `Add`/`Mul`/`Fun` argument list. |
| `NullSequence` | Three trailing `_` in the name (e.g., `x___`; three leading `_` also work, e.g. `___x`) | Matches zero or more `Atom`s. |

**Example**:

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

**Function**: Attempts to match `pattern` against `atom`, returning the binding set on success. Uses the default backtracking budget (`DEFAULT_MAX_BACKTRACKS = 10_000`). `Add`/`Mul` use AC matching (full backtracking search); `Fun` uses ordered matching (with sequence wildcard support).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The matching pattern. |
| `atom` | `Atom<'a>` | The expression to match against. |

**Return value**: `Result<Bindings<'a>, MatchError>`.

**Example**:

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

**Function**: Same as `match_pattern`, but allows a custom backtracking budget to prevent exponential-time blowups on pathological inputs.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pattern` | `Pattern<'a>` | The matching pattern. |
| `atom` | `Atom<'a>` | The expression to match against. |
| `max_backtracks` | `usize` | The maximum number of backtracks. |

**Return value**: `Result<Bindings<'a>, MatchError>`.

**Errors**: Returns `MatchError::BudgetExhausted` if the backtracking limit is exceeded.

**See also**: [match_pattern](#match_pattern)

### Bindings

```rust
use ocas_rewrite::matcher::Bindings;
```

```rust
#[derive(Debug, Clone, Default)]
pub struct Bindings<'a> { /* private fields */ }
```

The set of wildcard bindings produced by a successful match.

#### Bindings::new

```rust
pub fn new() -> Self
```

**Function**: Creates an empty binding set.

#### Bindings::get

```rust
pub fn get(&self, name: Symbol) -> Option<&MatchValue<'a>>
```

**Function**: Looks up the bound value by wildcard name.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `Symbol` | The wildcard name (without trailing underscores). |

**Return value**: `Option<&MatchValue<'a>>`.

**Example**:

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

The value type of a wildcard binding.

| Variant | Description |
|---|---|
| `Single(Atom<'a>)` | A `Single`-level wildcard bound to a single `Atom`. |
| `Sequence(&'a [Atom<'a>])` | A `Sequence`- or `NullSequence`-level wildcard bound to a slice of `Atom`s. |

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

Pattern-matching errors.

| Variant | Display | Description |
|---|---|---|
| `NoMatch` | `"pattern did not match"` | The pattern does not match the target expression. |
| `InconsistentBinding` | `"inconsistent wildcard binding"` | The same wildcard name is bound to different values at different positions. |
| `BudgetExhausted` | `"backtrack budget exhausted"` | The backtracking budget of AC matching is exhausted. Indicates the input may cause exponential-time matching; increase the budget or simplify the pattern. |

---

## See also

- [Rewrite & Simplification](./rust-rewrite.md) — the full rewrite engine: `Rule`, `simplify`, `replace_all`, and more
- [Calculus](./rust-calculus.md) — `diff`, `integrate`, `taylor`, and other operations built on the expression system
- [Coefficient Domains](./rust-domains.md) — the `Domain` trait and the various domain types
