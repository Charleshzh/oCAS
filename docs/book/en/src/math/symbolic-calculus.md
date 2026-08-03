# Advanced: Symbolic Calculus

## Prerequisites

- Calculus basics: definition of the derivative, basic differentiation rules, the concept of Taylor series
- Expression trees and recursive data structures (see [Polynomial Algebra](./polynomial-algebra.md))
- The oCAS expression system (`Atom`, `AtomNode`, `AtomArena`)

Suggested reading: [Polynomial Algebra](./polynomial-algebra.md), [Linear Algebra](./linear-algebra.md).

---

## Basic Concepts

### Symbolic Calculus on Expression Trees

In symbolic computation, a function $f(x)$ is represented as an **expression tree**. The leaf nodes are constants or variables; the internal nodes are operators (`Add`, `Mul`, `Pow`) or function calls (`sin`, `exp`, `log`, etc.).

For example, the expression tree of $f(x) = x^2 + \sin(x)$ is:

```
Add
├── Pow
│   ├── Var("x")
│   └── Num(2)
└── Fun("sin", [Var("x")])
```

The core idea of symbolic differentiation is **pattern matching**: for each node of the expression tree, apply the differentiation rule corresponding to its type, recursively building the derivative expression tree.

### Basic Differentiation Rules

Given differentiable functions $f(x)$ and $g(x)$:

| Rule | Formula |
|---|---|
| Constant | $\frac{d}{dx}[c] = 0$ |
| Identity | $\frac{d}{dx}[x] = 1$ |
| Sum | $\frac{d}{dx}[f + g] = f' + g'$ |
| Product | $\frac{d}{dx}[f \cdot g] = f' \cdot g + f \cdot g'$ |
| Chain | $\frac{d}{dx}[f(g(x))] = f'(g(x)) \cdot g'(x)$ |
| Power (constant exponent $n$) | $\frac{d}{dx}[f^n] = n \cdot f^{n-1} \cdot f'$ |
| Exponential (constant base $a$) | $\frac{d}{dx}[a^g] = a^g \cdot \ln a \cdot g'$ |
| General power | $\frac{d}{dx}[f^g] = f^g \cdot (\ln f \cdot g' + g \cdot f'/f)$ |

### Derivatives of Elementary Functions

| $f(u)$ | $f'(u)$ |
|---|---|
| $\sin u$ | $\cos u$ |
| $\cos u$ | $-\sin u$ |
| $\exp u$ | $\exp u$ |
| $\ln u$ | $u^{-1}$ |
| $\sqrt{u}$ | $(2\sqrt{u})^{-1}$ |
| $\tan u$ | $\sec^2 u$ |
| $\sec u$ | $\sec u \cdot \tan u$ |
| $\arctan u$ | $(1 + u^2)^{-1}$ |

Note: the table gives $f'(u)$ (the derivative with respect to $u$); in fact $\frac{d}{dx}[f(u(x))] = f'(u) \cdot u'(x)$, i.e. one still has to multiply by the chain-rule factor $u'(x)$.

### Taylor Series

Suppose $f$ is infinitely differentiable in a neighbourhood of the point $a$. The **Taylor series** of $f$ at $a$ is:

$$
f(x) = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x - a)^n
$$

where $f^{(n)}(a)$ is the $n$-th derivative of $f$ at $a$. Truncating at order $N$ yields the **Taylor polynomial**:

$$
T_N(x) = \sum_{n=0}^{N} \frac{f^{(n)}(a)}{n!} (x - a)^n
$$

Common expansions (at $a = 0$, i.e. Maclaurin series):

$$
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots
$$

$$
\sin x = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \cdots
$$

$$
\cos x = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \cdots
$$

$$
\ln(1+x) = x - \frac{x^2}{2} + \frac{x^3}{3} - \cdots \quad (|x| < 1)
$$

---

## Core Theory

### Recursive Differentiation on Expression Trees

The core algorithm of symbolic differentiation is **structural recursion** over the expression tree. For each node type, the corresponding differentiation rule is applied and the subexpressions are differentiated recursively.

**Algorithm** `diff(expr, var)`:

1. **Constant** `Num(c)`: return `Num(0)`.
2. **Variable** `Var(v)`: if $v = \text{var}$, return `Num(1)`; otherwise return `Num(0)`.
3. **Sum** `Add(a₁, a₂, …, aₙ)`: differentiate term by term, returning `Add(diff(a₁), diff(a₂), …, diff(aₙ))`.
4. **Product** `Mul(a₁, a₂, …, aₙ)`: apply the generalised product rule. For the product $a_1 a_2 \cdots a_n$ of $n$ factors:

$$
\frac{d}{dx}\left[\prod_{i=1}^n a_i\right] = \sum_{i=1}^n \left(\frac{d a_i}{dx} \cdot \prod_{j \neq i} a_j\right)
$$

That is, for each position $i$, replace the $i$-th factor by its derivative, keep the other factors unchanged, and sum the $n$ results.

5. **Power** `Pow(base, exp)`: three cases —
   - Constant exponent $n$: $\frac{d}{dx}[b^n] = n \cdot b^{n-1} \cdot b'$
   - Constant base $a$: $\frac{d}{dx}[a^u] = a^u \cdot \ln a \cdot u'$
   - General case (generalised power rule): $\frac{d}{dx}[b^e] = b^e \cdot (\ln b \cdot e' + e \cdot b'/b)$

6. **Function call** `Fun(name, [u])`: look up $f'(u)$ in the table, then multiply by $u'$ (chain rule).
7. **Unknown function**: return the unevaluated form `Derivative(f(x), x)`.

**Key property**: because oCAS expressions use hash-consing (structural sharing), the intermediate expressions produced by recursive differentiation automatically share identical subtrees. After differentiation, the rewrite engine simplifies the result (e.g. $0 \cdot f \to 0$, $1 \cdot f \to f$), eliminating redundant terms.

### The Recurrence Algorithm for Taylor Expansion

A direct implementation of Taylor expansion can exploit **successive differentiation**:

**Algorithm** `taylor(expr, var, point, order)`:

1. Let `current = expr`, `sum = 0`.
2. For $n = 0, 1, \dots, \text{order}$:
   - Compute $f^{(n)}(\text{point})$: replace `var` in `current` by `point` (i.e. evaluate).
   - Compute the coefficient $c_n = f^{(n)}(\text{point}) / n!$.
   - If $n = 0$ the term is $c_0$; otherwise the term is $c_n \cdot (x - \text{point})^n$.
   - Accumulate: $\text{sum} \mathrel{+}= \text{term}$.
   - If $n < \text{order}$, update `current = diff(current, var)`.
3. Simplify and return `sum`.

**Handling factorials**: $1/n!$ is represented as a power $n!^{-1}$ to avoid introducing floating-point numbers. The `mul_by_factorial_inverse` function computes $n!$ (using `i64` multiplication, safe for $n \leq 20$) and then constructs $c \cdot (n!)^{-1}$.

**Complexity**: in the loop $n = 0, \dots, N$, each step performs one substitution (evaluation), but a symbolic differentiation is performed only when $n < N$. Symbolic differentiation itself is $O(|\text{expr}|)$ ($|\text{expr}|$ being the size of the expression tree), but simplification may add nodes. Overall, a Taylor expansion of order $N$ requires $N$ differentiations and $N+1$ substitutions.

### Semantics of Substitution

`substitute(expr, var, replacement)` replaces every free occurrence of the variable `var` in the expression by `replacement`.

**Implementation strategy**: oCAS uses a **bottom-up transform** — the `transform` function walks the expression tree, and for each leaf node checks whether it is the variable to be replaced. If so it returns `replacement`; otherwise the node is left unchanged.

**Deep copy vs. references**:

- oCAS's `Atom` is a **Copy handle** (a pointer) into an arena, not an owning value. `transform` allocates new nodes on the arena to build the result, but unmodified subtrees reuse the original handles (zero copy).
- This makes `substitute` **purely functional**: the original expression is unchanged, and the result is a newly allocated expression tree in which the replaced paths point to new nodes while unaffected subtrees are shared with the original tree.
- Unlike "destructive modification", this guarantees referential transparency: substituting into the same expression repeatedly does not interfere.

**Limitations of substitution**: `substitute` is a **syntactic** substitution, not a semantic one. It replaces the literal variable name, not its "value". For example, replacing $t \to y$ inside $\int_0^x f(t)\, dt$ would wrongly change the integration variable — such semantic replacements require finer-grained binding analysis.

---

## Implementation in oCAS

### `diff`: Symbolic Differentiation

**File**: `ocas-calc/src/derivative.rs`

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**Structure of the implementation**:

`diff` is the entry point: it calls the internal `diff_raw` to perform the actual recursion, then simplifies the result with `simplify`, and finally normalizes it with `normalize`.

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = diff_raw(ctx, expr, var);
    let simplified = simplify(ctx, raw, &rules, 20);
    normalize(ctx, simplified)
}
```

`diff_raw` pattern-matches on every variant of `AtomNode`:

- **`Num(_)`**: constant, derivative 0.
- **`Var(v)`**: derivative 1 if `v == var`, otherwise 0.
- **`Add(args)`**: call `diff_raw` on each argument and collect the results into a new `Add` node.
- **`Mul(args)`**: the generalised product rule — for $n$ factors, generate $n$ product terms (each with one factor replaced by its derivative), then add them.
- **`Pow(base, exp)`**: three cases (constant exponent / constant base / general), as above.
- **`Fun(name, args)`**: delegate to `diff_function`.

`diff_function` maintains a **hard-coded derivative table** covering the eight functions `sin`, `cos`, `exp`, `log`, `sqrt`, `tan`, `sec`, `atan`. For each function $f(u)$ the table stores $f'(u)$ (the derivative with respect to $u$), which is then multiplied by the chain-rule factor $\frac{du}{dx}$.

For functions not in the table, an **unevaluated form** is returned:

```rust
ctx.fun("Derivative", &[ctx.fun(name_str, args), ctx.var(var.as_str())])
```

namely `Derivative(f(x), x)`, representing a derivative that could not be computed automatically.

**Simplification pipeline**: the raw derivative produced by `diff` may contain redundancy (such as `0 * sin(x)`, `1 * cos(x)`). `simplify` runs `calculus_rules` (calculus-specific rules plus the default simplification rules) for at most 20 rounds of fixpoint simplification, and then `normalize` brings the result into canonical form.

### `taylor`: Taylor Expansion

**File**: `ocas-calc/src/series.rs`

```rust
pub fn taylor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    point: Atom<'a>,
    order: usize,
) -> Atom<'a>
```

**Implementation details**:

1. Precompute `(x - point)` as the common subexpression `x_minus_p`.
2. Loop $n = 0, \dots, \text{order}$:
   - `substitute(ctx, current, var, point)` evaluates $f^{(n)}(a)$.
   - `mul_by_factorial_inverse` computes $f^{(n)}(a) / n!$.
   - Build the term $c_n \cdot (x - a)^n$ (the power factor is omitted for $n = 0$).
   - Accumulate into `sum`.
   - If $n < \text{order}$, `current = diff(ctx, current, var)` computes the next derivative.
3. Simplify and normalize.

**`mul_by_factorial_inverse`**:

```rust
fn mul_by_factorial_inverse<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, n: usize) -> Atom<'a> {
    if n == 0 { return expr; }
    let mut fact: i64 = 1;
    for i in 2..=n {
        fact = fact.checked_mul(i as i64).expect("factorial fits in i64");
    }
    ctx.mul(&[expr, ctx.pow(ctx.num(fact), ctx.num(-1))])
}
```

$1/n!$ is represented as $(n!)^{-1}$, i.e. `Pow(Num(n!), Num(-1))`, preserving an exact rational representation. `checked_mul` guarantees $n \leq 20$ ($20! = 2\,432\,902\,008\,176\,640\,000 < 2^{63}$).

**Example**:

```rust
// expand e^x at x=0 to order 3
let x = ctx.var("x");
let expr = ctx.fun("exp", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
// result = "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))"
// i.e. 1 + x + x²/2 + x³/6

// expand sin(x) at x=0 to order 5
let sin_x = ctx.fun("sin", &[x]);
let result = taylor(&ctx, sin_x, Symbol::new("x"), ctx.num(0), 5);
// result = "x + (-1*(6^-1)*(x^3)) + ((120^-1)*(x^5))"
// i.e. x - x³/6 + x⁵/120 (the even-order coefficients vanish automatically)
```

### `substitute`: Variable Substitution

**File**: `ocas-calc/src/series.rs`

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a>
```

**Implementation**:

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    transform(ctx, expr, |a| match a.node() {
        AtomNode::Var(v) if *v == var => Some(replacement),
        _ => None,
    })
}
```

The `transform` function (bottom-up transform) walks the expression tree. The closure checks each node: if it is `Var(v)` with `v == var`, it returns `Some(replacement)` (replace); otherwise `None` (leave unchanged).

**Key properties**:

- **Purely functional**: the original expression is not modified; a new tree is returned. Unmodified subtrees are shared through `Atom`'s Copy semantics.
- **Deep substitution**: all free occurrences are replaced, including those nested in function arguments (such as the `x` in `sin(x)`).
- **No capture check**: there is no check whether the substitution causes variable capture — this is a syntactic-level substitution.

### `integrate`: Symbolic Integration

**File**: `ocas-calc/src/integral/mod.rs`

```rust
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

oCAS's integrator uses a **layered pipeline** architecture (see [Symbolic Integration Algorithms](../algorithms/integration.md) and [The Risch Algorithm](./risch-algorithm.md)):

1. **Direct lookup layer** (`integrate_raw`): direct integration rules for constants, variables, sums, and simple products/powers/functions.
2. **Rational function layer** (`rational.rs`): partial fraction decomposition + Hermite reduction + the logarithmic part (logarithmic-derivative identity, completing the square, Rothstein–Trager).
3. **Risch algorithm layer** (`risch.rs`): structural theorems over the differential field tower $\mathbb{Q}(x, t_1, \dots, t_n)$.
4. **Trigonometric layer** (`trig.rs`): rewriting via Euler's formula into complex exponentials ($t = e^{ix}$) → retry Risch on the rewritten form → `realify` back to real form.
5. **Special functions layer** (`special.rs`): tables of special functions such as erf, erfi, Ei, Si, Ci, Shi, Chi, Fresnel.
6. **Heuristic layer** (`heuristic.rs`): integration by parts (LIATE ordering), trigonometric substitution ($\sqrt{a^2 - x^2} \to a\sin\theta$), Weierstrass rationalization $t = \tan(x/2)$, Euler substitution.
7. **Fallback**: return the unevaluated form `Integral(expr, var)`.

The `integrate_with_fuel` variant accepts a `Fuel` budget to prevent infinite loops in the simplification phase. The integration traversal itself is bounded by the internal recursion depth limit `MAX_DEPTH = 8` (the Risch layer has its own `MAX_RISCH_DEPTH = 16`):

```rust
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>>
```

---

## Advanced Topics

### The Interaction of Differentiation and Simplification

The correctness of `diff` depends not only on the differentiation rules themselves but also on the quality of simplification. For example:

- $\frac{d}{dx}[x^2 + 3x] = 2x + 3$ — this requires simplifying $1 \cdot 3$ to $3$.
- $\frac{d}{dx}[e^{\ln x}] = e^{\ln x} \cdot \frac{1}{x} = 1$ — this requires simplifying $e^{\ln x} \to x$.

oCAS handles such cases through multi-round simplification with `calculus_rules` and the default rules. However, some identities (such as $\sin^2 x + \cos^2 x = 1$) may require dedicated trigonometric simplification rules.

### Unevaluated Derivatives and Integrals

When `diff` encounters an unknown function (a function not in the derivative table), it returns `Derivative(f, x)`. Similarly, `integrate` returns `Integral(f, x)` for integrals it cannot handle.

These unevaluated forms are **valid expression nodes** — they can participate in subsequent symbolic operations. For example, $\frac{d}{dx}[\text{Derivative}(f, x)]$ returns `Derivative(Derivative(f, x), x)`, representing the second derivative.

### Bridging to Numerical Evaluation

Expressions produced by symbolic differentiation can be passed to the `ExpressionEvaluator` (see [Evaluation and JIT](../api/rust-evaluation.md)) for efficient numerical evaluation. A typical workflow:

1. Symbolic differentiation: `diff(&ctx, f, x)` → the derivative expression $f'(x)$.
2. Compile the evaluator: `ExpressionEvaluator::compile(f')`.
3. Numerical evaluation: `evaluator.evaluate(&[x_val])` → $f'(x_{\text{val}})$.

This is far more accurate than numerical differentiation (finite differences) and avoids the problem of step-size selection.

---

## References

1. **Geddes, K. O., Czapor, S. R. & Labahn, G.** *Algorithms for Computer Algebra.* Kluwer Academic Publishers, 1992. — Chapter 11, "Symbolic Differentiation", covers recursive differentiation on expression trees, simplification strategies, and derivative-table design.
2. **Bronstein, M.** *Symbolic Integration I: Transcendental Functions.* Springer, 2005. — The authoritative reference on the Risch algorithm and the design basis of oCAS's integrator.
3. **Cohen, H.** *A Course in Computational Algebraic Number Theory.* Springer, 1993. — Chapter 1 contains general methodology for complexity analysis of algorithms.
