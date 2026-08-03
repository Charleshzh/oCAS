# Calculus

> Source: `ocas-calc/src/`

oCAS's calculus module provides symbolic differentiation, integration, series expansion, and partial fraction decomposition. All operations execute on an `AtomArena` and produce new `Atom` handles. Derivatives or integrals that cannot be simplified are returned in unevaluated form as `Derivative(expr, var)` or `Integral(expr, var)`.

**Module structure**:

| Submodule | Purpose |
|---|---|
| `derivative` | symbolic differentiation `diff` |
| `integral` | layered integration pipeline `integrate`, `integrate_heuristic`, `integrate_with_fuel` |
| `series` | Taylor expansion `taylor`, substitution `substitute` |
| `partial_fraction` | partial fraction decomposition `apart`, recomposition `together` |
| `tower` | differential field tower construction for the Risch algorithm |
| `ode` | ordinary differential equation solving (see [Solvers](./rust-solvers.md)) |

---

## `diff`

**Signature**:

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**Description**: computes the symbolic derivative of `expr` with respect to `var`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena used to construct the result |
| `expr` | `Atom<'a>` | the expression to differentiate |
| `var` | `Symbol` | the differentiation variable (interned string) |

**Returns**: `Atom<'a>` — the differentiated and simplified expression.

**Algorithm**:

`diff` recursively traverses the expression tree and applies the standard differentiation rules to each node type:

| Node type | Differentiation rule |
|---|---|
| `Num(_)` | $\frac{d}{dx} c = 0$ |
| `Var(s)` | $1$ if $s = \text{var}$, otherwise $0$ |
| `Add([a₁, …, aₙ])` | derivative of a sum = sum of derivatives: $\sum \frac{d a_i}{dx}$ |
| `Mul([a₁, …, aₙ])` | product rule: $\frac{d}{dx}\prod a_i = \sum_i \left(\frac{d a_i}{dx} \prod_{j \neq i} a_j\right)$ |
| `Pow([base, exp])` | power rule (with chain rule), supports constant and non-constant exponents |
| `Fun(name, args)` | built-in function table + chain rule |

**Built-in derivative table**:

| Function $f(u)$ | Derivative $f'(u) \cdot u'$ |
|---|---|
| $\sin(u)$ | $\cos(u) \cdot u'$ |
| $\cos(u)$ | $-\sin(u) \cdot u'$ |
| $\tan(u)$ | $\sec^2(u) \cdot u'$ |
| $\sec(u)$ | $\sec(u)\tan(u) \cdot u'$ |
| $\exp(u)$ | $\exp(u) \cdot u'$ |
| $\log(u)$ | $u^{-1} \cdot u'$ |
| $\sqrt{u}$ | $(2\sqrt{u})^{-1} \cdot u'$ |
| $\text{atan}(u)$ | $(1 + u^2)^{-1} \cdot u'$ |

The table is hard-coded in `diff_function` in `derivative.rs`. Functions **not** in the table (such as $\text{asin}$, $\sinh$, $\cosh$, $\tanh$, etc.) also return the unevaluated form `Derivative(f, var)`.

**Unevaluated form**: when the function being differentiated is not in the built-in table, the result is `Derivative(expr, var)`, which can be detected by later pattern matching (match a `Fun` node whose name is `"Derivative"`; note that `is_fallback` only detects `"Integral"`, not `Derivative`).

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::diff;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// d/dx [sin(x)] = cos(x)
let sin_x = ctx.fun("sin", &[x]);
let result = diff(&ctx, sin_x, Symbol::new("x"));
assert_eq!(result.to_string(), "cos(x)");

// d/dx [x^3] = 3*x^2
let x_cubed = ctx.pow(x, ctx.num(3));
let result = diff(&ctx, x_cubed, Symbol::new("x"));
assert_eq!(result.to_string(), "3*(x^2)");

// d/dx [exp(x^2)] = 2*x*exp(x^2) (chain rule)
let inner = ctx.pow(x, ctx.num(2));
let expr = ctx.fun("exp", &[inner]);
let result = diff(&ctx, expr, Symbol::new("x"));
// Output contains a simplified form of 2*x*exp(x^2)

// Unknown functions return the Derivative form
let f = ctx.fun("my_func", &[x]);
let result = diff(&ctx, f, Symbol::new("x"));
// result contains Derivative(my_func(x), x)
```

**See also**: [taylor](#taylor), [integrate](#integrate)

---

## `integrate`

**Signature**:

```rust
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**Description**: computes the symbolic integral of `expr` with respect to `var`. Uses a layered integration pipeline that tries different algorithms layer by layer.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena |
| `expr` | `Atom<'a>` | the integrand |
| `var` | `Symbol` | the integration variable |

**Returns**: `Atom<'a>` — the integration result. Returns `Integral(expr, var)` when unsolved.

**Integration pipeline**:

`integrate` internally calls `integrate_raw` (recursive, depth limit `MAX_DEPTH = 8`). The pipeline has two phases:

```
┌───────────────────────────────────────────┐
│ Phase 1: integrate_raw dispatches on node │
│   type                                    │
├───────────────────────────────────────────┤
│ Num / Var          → constant integration:│
│                      c·x; ∫x dx = (1/2)·x² │
│ Add                → integrate term-wise  │
│ Mul (integrate_product)  → extract the    │
│                      constant factor      │
│ Pow (integrate_power)    → power rule     │
│                      ∫xⁿ, ∫(a·x+b)ⁿ,      │
│                      fractional (a·x+b)^(p/q)│
│ Fun (integrate_function) → table + linear │
│                      substitution         │
│                      ∫sin(u), ∫cos(u),    │
│                      ∫exp(u), ∫log(u)     │
│                      if u = a·x+b, apply  │
│                      the chain rule       │
├───────────────────────────────────────────┤
│ Phase 2: when the dispatch above fails    │
│   (returns the fallback form)             │
│   → try_risch_or_fallback, in order:      │
├───────────────────────────────────────────┤
│ 1. rational function  (integrate_rational)│
│ 2. Risch algorithm    (risch_integrate)   │
│ 3. trig→exponential   (trig_to_exp +      │
│    + Risch + realify)                     │
│ 4. special functions  (special_integrate) │
│ 5. heuristics         (heuristic_integrate)│
│    parts (LIATE), trig substitution,      │
│    Weierstrass, Euler substitution        │
│ 6. unevaluated form   (fallback)          │
│    Integral(expr, var)                    │
└───────────────────────────────────────────┘
```

**Phase 1 — node-type dispatch**:

- `Num(_)`: the integral of a constant is the linear function `c * x`
- `Var(var)`: $\int x\,dx = \frac{1}{2} x^2$; other variables are treated as constants
- `Add`: integrate each term recursively and sum
- `Mul` (`integrate_product`): recognizes the form `c · f(x)` (a constant times a non-constant factor), extracts the constant factor, and integrates recursively
- `Pow` (`integrate_power`): power-function integration, supporting $\int x^n\,dx$, linear forms $(a \cdot x + b)^n$ (including $n = -1$: $\int \frac{dx}{ax+b} = \frac{\log(ax+b)}{a}$), and fractional exponents $(a \cdot x + b)^{p/q}$
- `Fun` (`integrate_function`): table lookup + linear substitution

**Table lookup and linear substitution**:

The built-in integral table covers power functions (`integrate_power`) plus `sin`, `cos`, `exp`, `log`. When the argument of the integrand (or the base of a power) has the linear form $a \cdot x + b$, the chain rule $\int f(ax+b)\,dx = \frac{1}{a} F(ax+b)$ is applied automatically.

**Phase 2 — fallback chain** (`try_risch_or_fallback`):

1. **Rational function integration** (`integrate_rational`): when the expression is a rational function of $\var$ (with no other variables or function applications), perform:
   1. Integrate the polynomial part term by term
   2. Hermite reduction to separate the rational part
   3. Logarithmic part: completing the square (quadratic denominators → $\log$ or $\text{atan}$) + Rothstein–Trager resultant

   See [rational function integration](#rational-function-integration).

2. **Risch algorithm** (`risch_integrate`): constructs the differential field tower $\mathbb{Q}(x, t_1, \dots, t_n)$ and integrates recursively layer by layer. See [Risch algorithm](#risch-algorithm).

3. **Trigonometric rewriting**: when Risch fails on the original expression and the integrand contains trigonometric functions, they are rewritten into complex exponential form via `trig_to_exp`, Risch is tried again, and `realify` converts the answer back to real form.

4. **Special functions** (`special_integrate`): when both Risch and the trigonometric rewriting fail, tries to match known special-function antiderivatives. See [special function integration](#special-function-integration).

5. **Heuristics** (`heuristic_integrate`): integration by parts (LIATE), trigonometric substitution, Weierstrass substitution, Euler substitution. See [integrate_heuristic](#integrate_heuristic).

6. **Unevaluated form** (`fallback`): when everything above fails, returns `Integral(expr, var)`. Detectable via `is_fallback`:

```rust
pub(crate) fn is_fallback<'a>(atom: &Atom<'a>) -> bool {
    matches!(atom.node(), AtomNode::Fun(name, _) if name.as_str() == "Integral")
}
```

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integrate;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ x² dx = (1/3) x³
let expr = ctx.pow(x, ctx.num(2));
let result = integrate(&ctx, expr, Symbol::new("x"));
assert_eq!(result.to_string(), "(3^-1)*(x^3)");

// ∫ sin(x) dx = -cos(x)
let expr = ctx.fun("sin", &[x]);
let result = integrate(&ctx, expr, Symbol::new("x"));
// Output: a simplified form of -cos(x)

// ∫ 1/x dx = log(x) (via the rational function integrator or Risch)
let expr = ctx.pow(x, ctx.num(-1));
let result = integrate(&ctx, expr, Symbol::new("x"));
// Output: log(x)

// Non-integrable expressions return the Integral form
let expr = ctx.fun("my_func", &[x]);
let result = integrate(&ctx, expr, Symbol::new("x"));
// result is Integral(my_func(x), x)
```

**See also**: [integrate_heuristic](#integrate_heuristic), [integrate_with_fuel](#integrate_with_fuel), [rational function integration](#rational-function-integration), [Risch algorithm](#risch-algorithm)

---

## `integrate_heuristic`

**Signature**:

```rust
pub fn integrate_heuristic<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**Description**: attempts heuristic integration techniques. Uses only integration by parts, trigonometric substitution, Weierstrass substitution, and Euler substitution; does not invoke the Risch algorithm or the rational function integrator.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena |
| `expr` | `Atom<'a>` | the integrand |
| `var` | `Symbol` | the integration variable |

**Returns**: `Atom<'a>` — the integration result if a heuristic succeeds, otherwise `Integral(expr, var)`.

**Heuristic methods**:

### 1. Integration by parts (LIATE heuristic)

For products $\int u \cdot v' \, dx$, choose $u$ by **LIATE** priority (highest to lowest):

| Priority | Type | LIATE score | Example |
|---|---|---|---|
| 0 (highest) | Logarithmic (L) | 0 | $\log(x)$ |
| 1 | Inverse trigonometric (I) | 1 | $\text{asin}(x)$ |
| 2 | Algebraic (A) | 2 | $x^2$ |
| 3 | Trigonometric (T) | 3 | $\sin(x)$ |
| 4 | Exponential (E) | 4 | $\exp(x)$ |
| 5 (lowest) | Other | 5 | — |

Choose the factor with the lowest score as $u$ and the rest as $v'$. The recursion depth limit is `PARTS_MAX_DEPTH = 2`.

### 2. Trigonometric substitution

Matches the following patterns and returns the known antiderivative:

| Integrand | Substitution | Result |
|---|---|---|
| $\frac{1}{\sqrt{a^2 - x^2}}$ | $x = a\sin\theta$ | $\text{asin}(x/a)$ |
| $\frac{1}{\sqrt{a^2 + x^2}}$ | $x = a\sinh t$ | $\text{asinh}(x/a)$ |
| $\frac{1}{\sqrt{x^2 - a^2}}$ | $x = a\cosh t$ | $\text{acosh}(x/a)$ |
| $\sqrt{a^2 - x^2}$ | $x = a\sin\theta$ | $\frac{x\sqrt{a^2-x^2} + a^2\,\text{asin}(x/a)}{2}$ |
| $\sqrt{a^2 + x^2}$ | $x = a\sinh t$ | $\frac{x\sqrt{a^2+x^2} + a^2\,\text{asinh}(x/a)}{2}$ |
| $\sqrt{x^2 - a^2}$ | $x = a\cosh t$ | $\frac{x\sqrt{x^2-a^2} - a^2\,\text{acosh}(x/a)}{2}$ |

### 3. Weierstrass substitution

When the integrand is a rational function of $\sin(u)$ and $\cos(u)$ ($u$ linear in $\text{var}$), apply the substitution $t = \tan(u/2)$:

$$\sin(u) = \frac{2t}{1+t^2}, \quad \cos(u) = \frac{1-t^2}{1+t^2}$$

Rationalize the trigonometric rational function into a rational function of $t$, then integrate recursively.

### 4. Euler substitution

When the integrand contains $\sqrt{ax^2 + bx + c}$, try the Euler substitution to eliminate the radical.

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integrate_heuristic;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ x·exp(x) dx — integration by parts (LIATE: x=algebraic, exp=exponential)
let expr = ctx.mul(&[x, ctx.fun("exp", &[x])]);
let result = integrate_heuristic(&ctx, expr, Symbol::new("x"));
// Result contains a simplified form of (x-1)*exp(x)

// ∫ 1/√(1 - x²) dx = asin(x)
// Build (1 - x²)^(-1/2): the exponent is -1·(2^-1)
let one = ctx.num(1);
let sqrt_arg = ctx.add(&[one, ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
let neg_half = ctx.mul(&[ctx.num(-1), ctx.pow(ctx.num(2), ctx.num(-1))]);
let expr = ctx.pow(sqrt_arg, neg_half);
let result = integrate_heuristic(&ctx, expr, Symbol::new("x"));
// Output: a simplified form of asin(x)
```

**See also**: [integrate](#integrate)

---

## `integrate_with_fuel`

**Signature**:

```rust
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>, FuelError>
```

**Description**: integration with a fuel budget. The integration traversal itself uses an internal depth limit (`MAX_DEPTH = 8`); this entry point threads `fuel` into the two simplification phases after integration, preventing pathological results from making the rewrite engine loop indefinitely.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena |
| `expr` | `Atom<'a>` | the integrand |
| `var` | `Symbol` | the integration variable |
| `fuel` | `&Fuel` | the fuel budget (rewrite step limit) |

**Returns**:
- `Ok(Atom<'a>)` — integration simplification succeeded
- `Err(FuelError)` — fuel exhausted during simplification

**Difference from `integrate`**:

`integrate` imposes no limit on the simplification phases; `integrate_with_fuel` applies a fuel constraint to them. The depth limit (`MAX_DEPTH`) of the integration traversal itself is the same for both.

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_core::fuel::Fuel;
use ocas_calc::integrate_with_fuel;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, x, Symbol::new("x"), &fuel);
match result {
    Ok(expr) => assert_eq!(expr.to_string(), "(2^-1)*(x^2)"),
    Err(_) => panic!("fuel exhausted"),
}
```

**See also**: [integrate](#integrate)

---

## `taylor`

**Signature**:

```rust
pub fn taylor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    point: Atom<'a>,
    order: usize,
) -> Atom<'a>
```

**Description**: computes the Taylor expansion of `expr` in `var` around `point`, up to and including `order`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena |
| `expr` | `Atom<'a>` | the expression to expand |
| `var` | `Symbol` | the expansion variable |
| `point` | `Atom<'a>` | the expansion point (any expression, usually `0` or a constant) |
| `order` | `usize` | the expansion order (inclusive) |

**Returns**: `Atom<'a>` — the truncated polynomial:

$$\sum_{n=0}^{\text{order}} \frac{f^{(n)}(\text{point})}{n!} \cdot (\text{var} - \text{point})^n$$

**Algorithm**:

Computes each coefficient by repeated symbolic differentiation and evaluation at the expansion point:

1. For $n = 0, 1, \dots, \text{order}$:
   - Compute $f^{(n)}(x)$ (via `diff`)
   - Evaluate at $x = \text{point}$ with `substitute` to obtain $f^{(n)}(\text{point})$
   - Multiply by $\frac{1}{n!}$ (via `mul_by_factorial_inverse`)
   - Multiply by $(x - \text{point})^n$
2. Sum all terms and apply the calculus simplification rules

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::taylor;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// exp(x) expanded at x=0 to order 3
let expr = ctx.fun("exp", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
assert_eq!(
    result.to_string(),
    "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))"
);

// sin(x) expanded at x=0 to order 5
let expr = ctx.fun("sin", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 5);
// Output: x + (-1*(6^-1)*(x^3)) + ((120^-1)*(x^5))
```

**See also**: [diff](#diff), [substitute](#substitute)

---

## `substitute`

**Signature**:

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a>
```

**Description**: replaces every occurrence of the variable `var` in `expr` with `replacement`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | the expression arena |
| `expr` | `Atom<'a>` | the expression to rewrite |
| `var` | `Symbol` | the variable name to replace |
| `replacement` | `Atom<'a>` | the replacement expression |

**Returns**: `Atom<'a>` — the new expression after substitution. The original expression is unaffected (immutable data structure on the arena).

**Semantics**:

- Traverses the expression tree depth-first
- Replaces only leaf nodes `Var(Symbol)` that exactly match `var`
- No automatic simplification after substitution (call `simplify` separately if needed)

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::substitute;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");

// (x² + x)[x → y] = y² + y
let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), x]);
let result = substitute(&ctx, expr, Symbol::new("x"), y);
assert_eq!(result.to_string(), "(y^2) + y");

// sin(x)[x → 2*y] = sin(2*y)
let expr = ctx.fun("sin", &[x]);
let two_y = ctx.mul(&[ctx.num(2), y]);
let result = substitute(&ctx, expr, Symbol::new("x"), two_y);
assert_eq!(result.to_string(), "sin(2*y)");
```

**See also**: [taylor](#taylor), [diff](#diff)

---

## `apart`

**Signature**:

```rust
pub fn apart<D: EuclideanDomain>(
    num: &DenseUnivariatePolynomial<D>,
    den: &DenseUnivariatePolynomial<D>,
) -> (
    Option<DenseUnivariatePolynomial<D>>,
    Vec<PartialFractionTerm<D>>,
)
```

**Description**: performs partial fraction decomposition of the rational function $\frac{\text{num}(x)}{\text{den}(x)}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `num` | `&DenseUnivariatePolynomial<D>` | the numerator polynomial |
| `den` | `&DenseUnivariatePolynomial<D>` | the denominator polynomial |

**Returns**: `(Option<poly>, Vec<PartialFractionTerm<D>>)` — the polynomial part (when the numerator degree ≥ the denominator degree) and the list of partial fraction terms.

**Mathematical background**:

Given a proper fraction $\frac{p(x)}{q(x)}$ (i.e., $\deg(p) < \deg(q)$), factor the denominator square-free as $q = \prod_i f_i^{e_i}$, then decompose:

$$\frac{p(x)}{q(x)} = \text{poly\_part} + \sum_i \sum_{k=1}^{e_i} \frac{a_{i,k}(x)}{f_i(x)^k}$$

If $\deg(p) \geq \deg(q)$, perform polynomial division first; the quotient becomes `poly_part` and the remainder is decomposed.

**The `PartialFractionTerm<D>` struct**:

```rust
pub struct PartialFractionTerm<D: EuclideanDomain> {
    pub numer: DenseUnivariatePolynomial<D>,  // numerator a_{i,k}(x)
    pub denom: DenseUnivariatePolynomial<D>,  // irreducible factor f_i(x)
    pub exp: usize,                            // exponent k
}
```

Represents the fraction $\frac{\text{numer}(x)}{\text{denom}(x)^{\text{exp}}}$.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;

// 1 / (x² - 1) — the denominator is square-free (a single square-free factor x²-1)
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1),  // constant term
    Rational::new(0, 1),   // x coefficient
    Rational::new(1, 1),   // x² coefficient
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_none()); // deg(num) < deg(den), no polynomial part
// x²-1 is square-free → a single term with denominator x²-1 and exponent 1
// (apart works on the square-free factorization, not the irreducible
// factorization, so it does not split into (x-1)(x+1))

// (x³ + 1) / (x² - 1) — requires polynomial division
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(0, 1),
    Rational::new(0, 1), Rational::new(1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_some()); // the quotient polynomial x
```

**See also**: [together](#together), [rational function integration](#rational-function-integration)

---

## `together`

**Signature**:

```rust
pub fn together<D: EuclideanDomain>(
    poly_part: Option<&DenseUnivariatePolynomial<D>>,
    terms: &[PartialFractionTerm<D>],
) -> (DenseUnivariatePolynomial<D>, DenseUnivariatePolynomial<D>)
```

**Description**: recomposes partial fraction terms back into a single rational function. The inverse of `apart`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `poly_part` | `Option<&DenseUnivariatePolynomial<D>>` | the polynomial part (the first element returned by `apart`) |
| `terms` | `&[PartialFractionTerm<D>]` | the list of partial fraction terms |

**Returns**: `(numerator, denominator)` — the merged numerator and denominator polynomials.

**Algorithm**:

Combines the polynomial part and all fraction terms over a common denominator: compute the common denominator $\text{lcm}$, expand each term to the common denominator, sum, and return the resulting numerator and denominator.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::{apart, together};

let d = RationalDomain;
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1),
]);

// apart → round-trip → together should restore the original fraction
let (poly_part, terms) = apart(&num, &den);
let (result_num, result_den) = together(poly_part.as_ref(), &terms);
// result_num/result_den is equivalent to num/den (possibly up to a constant factor)
```

**See also**: [apart](#apart)

---

## Rational function integration

**Source**: `ocas-calc/src/integral/rational.rs`

The rational function integrator `integrate_rational` handles rational functions of $\var$ (with no other variables or function applications). It follows the standard three-step method from Chapter 2 of Bronstein's *Symbolic Integration I*:

### 1. Polynomial part

Integrate $\int \sum c_k x^k \, dx = \sum \frac{c_k}{k+1} x^{k+1}$ term by term.

### 2. Hermite reduction

For a proper fraction $\frac{a(x)}{d(x)}$ ($d$ monic), decompose:

$$\frac{a}{d} = \frac{d}{dx}\left(\frac{g_{\text{num}}}{g_{\text{den}}}\right) + \frac{a_1}{d_1}$$

where $d_1$ is square-free. The derivative of the rational part $\frac{g_{\text{num}}}{g_{\text{den}}}$ can be computed directly.

### 3. Logarithmic part

For the remaining part $\frac{a_1}{d_1}$ with square-free denominator, choose a strategy by denominator degree:

| Denominator degree | Strategy |
|---|---|
| 0 | constant integration |
| 1 | $c \cdot \log(ax + b)$ |
| 2 | complete the square → $\log$ or $\text{atan}$ |
| ≥ 3 | Rothstein–Trager resultant method |

**The Rothstein–Trager method**:

Compute $R(t) = \text{Res}_x(d, a - t \cdot d')$; the rational roots $c_i$ of $R(t)$ give:

$$\int \frac{a}{d} = \sum_i c_i \cdot \log(\gcd(d, a - c_i \cdot d'))$$

When $R(t)$ cannot be completely factored over $\mathbb{Q}$, the corresponding terms are returned in unevaluated form `Integral(term, var)`.

**Returns**:
- `Some(Atom)` — the expression is a rational function of $\var$; returns the integration result
- `None` — the expression is not a pure rational function of $\var$ (contains other variables or functions)

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integral::rational::integrate_rational;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ 1/x dx = log(x)
let result = integrate_rational(&ctx, ctx.pow(x, ctx.num(-1)), Symbol::new("x"));
assert_eq!(result.unwrap().to_string(), "log(x)");

// ∫ 1/(x²+1) dx = atan(x) (via completing the square)
// The expression must be built as a rational function of x
```

**See also**: [apart](#apart), [integrate](#integrate)

---

## Risch algorithm

**Source**: `ocas-calc/src/integral/risch.rs`, `ocas-calc/src/tower/`

The Risch algorithm handles integrals of elementary transcendental functions, based on Chapters 5–6 of Bronstein's *Symbolic Integration I*.

### Differential field tower

The algorithm first constructs the differential field tower $\mathbb{Q}(x, t_1, \dots, t_n)$, where each generator $t_i$ is a logarithm or an exponential:

```rust
pub(crate) enum GenKind {
    Constant,  // a constant symbol (e.g. the imaginary unit I), D t = 0
    Log,       // t_i = log(u)
    Exp,       // t_i = exp(u)
}
```

Tower construction (`build_tower`) walks the function applications in the expression, creating a new layer for every `log` or `exp`. Constraints:

- Only `log` and `exp` are accepted (trigonometric functions must first be rewritten via `trig_to_exp`)
- Algebraic functions are rejected (non-integer exponents, e.g. $\sqrt{x}$)
- Algebraically dependent generators are rejected (e.g., $\log(x)$ and $\log(2x)$)

### Layer-by-layer integration

At each tower layer $k(t_i)$, the algorithm performs:

1. **Hermite reduction** (`hermite_tower`): decompose $\frac{a}{d}$ into $D(g) + \frac{a_1}{d_1}$, where $d_1$ is square-free
2. **Polynomial part**:
   - **Primitive layer** ($t = \log(u)$): method of undetermined coefficients; the top-level constant is fixed by the logarithmic constraint
   - **Hyperexponential layer** ($t = \exp(u)$): solve the Risch differential equation $Dq + fq = g$ at each layer
3. **Logarithmic part**: match the logarithmic derivative identity $c \cdot \frac{D d}{d} \to c \cdot \log(d)$

### Handling trigonometric functions

`integrate` first tries Risch on the original expression; when that fails and the integrand contains trigonometric functions, they are rewritten into complex exponential form via `trig_to_exp` and Risch is attempted again:

$$\sin(u) \to \frac{e^{iu} - e^{-iu}}{2i}, \quad \cos(u) \to \frac{e^{iu} + e^{-iu}}{2}$$

The integration result is then converted back to real form where possible via `realify`:

| Complex pattern | Real form |
|---|---|
| $c \cdot \log(u + iv) + c \cdot \log(u - iv)$ | $c \cdot \log(u^2 + v^2)$ |
| $c \cdot \log(u + iv) - c \cdot \log(u - iv)$ | $2c \cdot \text{atan}(v/u)$ |
| $\exp(Iu) \cdot \exp(Iv)$ | $\exp(I(u+v))$ (combines the exponents so cancellation is visible, e.g. $e^{Ix}e^{-Ix} \to e^0$) |

If `realify` cannot match a pattern, the complex form is returned (still mathematically correct; verifiable by differentiation).

### Depth limit

`MAX_RISCH_DEPTH = 16`, enforced via the thread-local counter `RISCH_DEPTH` and the guard `RischDepthGuard` (RAII, decrements automatically).

**Returns**:
- `Some(Atom)` — integration succeeded
- `None` — the expression is beyond the implemented fragment (the caller falls back to other integrators)

**Example**:

```rust
// The Risch algorithm is invoked automatically through integrate; usually not called directly
// ∫ exp(x²) · 2x dx = exp(x²) — a typical example the Risch algorithm handles
// ∫ sin(x)/x dx — no elementary antiderivative; Risch returns None and the special function layer handles it
```

**See also**: [integrate](#integrate), [special function integration](#special-function-integration)

---

## Special function integration

**Source**: `ocas-calc/src/integral/special.rs`

When the Risch algorithm proves that the integral has no elementary antiderivative, the special function integrator tries to match known non-elementary antiderivatives.

### Supported special functions

| Integrand | Antiderivative | Function name |
|---|---|---|
| $e^{-x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erf}(x)$ | error function |
| $e^{x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erfi}(x)$ | imaginary error function |
| $\frac{e^x}{x}$ | $\text{Ei}(x)$ | exponential integral |
| $\frac{\sin x}{x}$ | $\text{Si}(x)$ | sine integral |
| $\frac{\cos x}{x}$ | $\text{Ci}(x)$ | cosine integral |
| $\frac{\sinh x}{x}$ | $\text{Shi}(x)$ | hyperbolic sine integral |
| $\frac{\cosh x}{x}$ | $\text{Chi}(x)$ | hyperbolic cosine integral |
| $\sin(x^2)$ | $\sqrt{\pi/2}\,\text{fresnels}(\sqrt{2/\pi}\,x)$ | Fresnel S integral (function name `fresnels`) |
| $\cos(x^2)$ | $\sqrt{\pi/2}\,\text{fresnelc}(\sqrt{2/\pi}\,x)$ | Fresnel C integral (function name `fresnelc`) |

**Matching strategy**:

`special_integrate` examines the factor structure of the integrand:
- `erf_family`: matches forms $e^{c \cdot x^2}$
- `ei_family`: matches forms $e^{cx} / x$
- `trig_integral_family`: matches $\sin(x)/x$, $\cos(x)/x$, $\sinh(x)/x$, $\cosh(x)/x$ (the argument must be exactly $x$)
- `fresnel_family`: matches $\sin(cx^2)$, $\cos(cx^2)$

**Returns**:
- `Some(Atom)` — matched successfully; returns the antiderivative containing special functions
- `None` — no known pattern matches; the caller returns `Integral(expr, var)`

**Design notes**:

The special function definitions are consistent with SymPy, and the results can be cross-validated with `sympy.integrate`.

**See also**: [Risch algorithm](#risch-algorithm), [integrate](#integrate)

---

## Helper functions

### `is_fallback`

```rust
pub(crate) fn is_fallback<'a>(atom: &Atom<'a>) -> bool
```

**Description**: checks whether `atom` is an unevaluated integral form (function name `"Integral"`). Note that unevaluated derivatives (the `Derivative` form) are constructed directly by `diff` and are not detected by this function.

**See also**: [integrate](#integrate), [diff](#diff)

### `is_constant`

```rust
pub(crate) fn is_constant<'a>(expr: Atom<'a>, var: Symbol) -> bool
```

**Description**: checks whether `expr` does not contain the variable `var` (i.e., is constant with respect to `var`).

### `fraction_exponent`

```rust
fn fraction_exponent<'a>(exp: Atom<'a>) -> Option<(i64, i64)>
```

**Description**: parses an exponent atom into a fraction $p/q$ (small integers). Accepts the form $p \cdot q^{-1}$; used for power-function integration.

### `linear_form`

```rust
fn linear_form<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol)
    -> Option<(Atom<'a>, Atom<'a>)>
```

**Description**: if `expr` is of the form $a \cdot \text{var} + b$ ($a, b$ constant with respect to `var`), returns `(a, b)`; otherwise returns `None`.

---

## Internal architecture

### Simplification pipeline

- The result of `integrate` passes through two simplification stages: first the default rewrite rules `default_rules` (algebraic identities such as $0 + x = x$, $1 \cdot x = x$, $x^1 = x$) for 20 steps, then `calculus_rules` for 10 steps, and finally `normalize` to canonical form
- `diff` and `taylor` simplify with `calculus_rules` directly, then `normalize`

`calculus_rules` (`rules.rs`) builds on `default_rules` and adds calculus-specific identities: `exp(0) → 1`, `log(1) → 0`, `sin(0) → 0`, `cos(0) → 1`, `tan(0) → 0`, $(-1) \cdot (-1) \to 1$, $1^x \to 1$, etc.

### Depth limits

| Limit | Value | Purpose |
|---|---|---|
| `MAX_DEPTH` | 8 | `integrate_raw` recursion depth |
| `MAX_RISCH_DEPTH` | 16 | Risch algorithm recursion depth |
| `PARTS_MAX_DEPTH` | 2 | integration-by-parts recursion depth |

### Module dependency graph

```
ocas-calc
├── derivative::diff
├── integral
│   ├── mod.rs          — integrate, integrate_with_fuel, integrate_heuristic
│   ├── heuristic.rs    — integration by parts, trig/Weierstrass/Euler substitutions
│   ├── rational.rs     — Hermite reduction, Rothstein–Trager
│   ├── risch.rs        — Risch algorithm core
│   ├── rde.rs          — Risch differential equation solver
│   ├── trig.rs         — trig→exponential rewriting, realify
│   ├── special.rs      — erf/Ei/Si/Ci/Fresnel
│   └── tower/
│       ├── build.rs    — differential field tower construction
│       ├── elem.rs     — tower element types (KElem, KPoly, KRat)
│       └── convert.rs  — Atom ↔ tower element conversion
├── series.rs           — taylor, substitute
├── partial_fraction.rs — apart, together
└── rules.rs            — calculus simplification rules
```

---

## See also

- [Math: Symbolic calculus](../math/symbolic-calculus.md) — the mathematics of differentiation and Taylor expansion
- [Math: The Risch algorithm](../math/risch-algorithm.md) — full mathematical derivation of the differential field tower and the Risch algorithm
- [Math: Polynomial GCD and factorization](../math/poly-gcd-factoring.md) — the mathematical foundations of Hermite reduction and square-free factorization
- [Solvers API](./rust-solvers.md) — ODE solving (`dsolve`, `classify_ode`)
- [Expression system](./rust-expressions.md) — basic usage of `Atom`, `AtomArena`, `Symbol`
