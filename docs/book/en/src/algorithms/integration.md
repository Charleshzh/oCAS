# Symbolic Integration

oCAS integrates expressions by a layered pipeline: a fast heuristic table,
then a rational-function integrator, then the Risch algorithm over
elementary towers, then a trigonometric-to-exponential rewrite, and
finally a special-function table. The first layer that produces an answer
wins. This chapter explains each layer and when the unevaluated
`Integral(expr, var)` form is returned.

---

## The Pipeline

`integrate(expr, var)` tries, in order:

1. **Fast table** — inline structural rules: power rules (`x^n`,
   `(ax+b)^n`, including fractional exponents), `sin`/`cos`/`exp` of
   linear arguments, and direct table entries such as `log(x)`. Fast and
   always attempted first.
2. **Rational-function integrator** — Hermite reduction plus the
   logarithmic part (logarithmic-derivative identity, completing the
   square, Rothstein–Trager). Handles every rational function of `x`.
3. **Risch algorithm** — elementary transcendental towers built from
   `log` and `exp` (the tower recursion is capped by
   `MAX_RISCH_DEPTH = 16`).
4. **Trigonometric rewrite** — `sin`/`cos`/`tan`/… rewritten into
   `exp(I·x)` and re-integrated by Risch, then converted back to real
   form on a best-effort basis.
5. **Special-function table** — non-elementary integrals with closed
   forms in terms of `erf`, `Ei`, `Si`, `Ci`, Fresnel `S`/`C`, …
6. **Heuristic techniques** — integration by parts, trigonometric
   substitution, Weierstrass $t = \tan(x/2)$, Euler substitutions. Tried
   after the special-function table, just before the unevaluated form.
7. **Unevaluated form** — `Integral(expr, var)`.

---

## Rational Functions

Every rational function of the integration variable is integrated
exactly. The polynomial part is integrated termwise; the proper fraction
is split by Hermite reduction into a rational part plus a remainder with
a squarefree denominator; the remainder yields logarithms (via the
identity `c·f'/f → c·log(f)`), arctangents (degree-2 denominators by
completing the square), or Rothstein–Trager logarithms.

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = parse(&ctx, "(2*x + 3)/(x^2 + 3*x + 5)").unwrap();
let result = integrate(&ctx, expr, Symbol::new("x"));
// log(x^2 + 3*x + 5)
```

---

## The Risch Algorithm

Elementary transcendental integrands are handled by building a
*differential field tower* `ℚ(x, t₁, …, tₙ)` where each `tᵢ` is a
`log` or `exp` over the field below, and integrating recursively
(Bronstein, *Symbolic Integration I*, ch. 5):

- at each level, the rational part is split off by Hermite reduction;
- the logarithmic part uses the logarithmic-derivative identity;
- the polynomial part is integrated by undetermined coefficients at
  `log` levels and by the Risch differential equation `Dq + f·q = g`
  at `exp` levels;
- the base `ℚ(x)` delegates to the rational-function integrator.

The tower recursion is capped by `MAX_RISCH_DEPTH = 16`: beyond that
depth the Risch layer gives up and hands the integrand to the next
layer (this keeps pathological integrands from recursing forever).

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
// ∫ x·exp(x) dx = (x - 1)·exp(x)
let result = integrate(&ctx, parse(&ctx, "x*exp(x)").unwrap(), Symbol::new("x"));
```

### Scope limits

The current fragment seeks only **polynomial** solutions of the Risch
differential equation and uses only the logarithmic-derivative identity
for logarithmic parts. Consequences:

- `∫ exp(x)/x dx` has no elementary antiderivative — it is answered by
  the special-function table as `Ei(x)`.
- Some `log`-tower cases needing a free-constant choice that makes lower
  layers integrable (e.g. `log(x+1)`) are not yet decided and fall back.

When no layer succeeds, the result is the unevaluated form
`Integral(expr, var)` — a deliberate answer, not an error.

---

## Trigonometric Integrands

`sin`, `cos`, `tan`, `cot`, `sec`, `csc` are rewritten into complex
exponentials via `t = exp(I·x)` and integrated by Risch. The imaginary
unit is carried as a constant tower generator (`D I = 0`). Results are
converted back to real form where possible: conjugate logarithm pairs
merge into real `log`/`atan` terms.

The Risch differential-equation solver currently works over `ℚ[x]`, so
hyperexponential equations whose coefficients contain `I` (e.g. the ones
produced by `sin(x)·cos(x)` or `cos(x)²`) cannot be solved yet; those
integrands return the unevaluated form. Simple `sin`/`cos` of linear
arguments are covered by the heuristic table.

---

## Special Functions

Integrals with no elementary antiderivative but a standard closed form
are answered directly (definitions match SymPy):

| Integrand | Result |
|---|---|
| `exp(-x²)` | `(√π/2)·erf(x)` |
| `exp(x²)` | `(√π/2)·erfi(x)` |
| `exp(c·x²)`, `c < 0` | `√π/(2√(-c))·erf(√(-c)·x)` |
| `exp(x)/x` | `Ei(x)` |
| `sin(x)/x` | `Si(x)` |
| `cos(x)/x` | `Ci(x)` |
| `sinh(x)/x` | `Shi(x)` |
| `cosh(x)/x` | `Chi(x)` |
| `sin(x²)` | `√(π/2)·fresnels(√(2/π)·x)` |
| `cos(x²)` | `√(π/2)·fresnelc(√(2/π)·x)` |

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
// ∫ exp(-x^2) dx = (√π/2)·erf(x)
let result = integrate(&ctx, parse(&ctx, "exp(-x^2)").unwrap(), Symbol::new("x"));
```

---

## Fuel-bounded integration

`integrate_with_fuel` wraps the same pipeline but threads a [`Fuel`] budget
through the two post-integration simplification stages. A pathological result
that would otherwise spin the rewriter can be cut off deterministically.

```rust
use ocas_core::arena::Arena;
use ocas_core::fuel::Fuel;
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integral::integrate_with_fuel;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.var("x");
let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, expr, Symbol::new("x"), &fuel);
match result {
    Ok(r) => println!("{}", r),
    Err(_) => println!("fuel exhausted during simplification"),
}
```

Returns `Err` only when fuel was exhausted mid-simplification. The
integration traversal itself uses the internal depth limit; fuel bounds
only the post-simplification passes.

---

## Bindings

The same pipeline backs the Python and C APIs:

- Python: `Expression.integrate(var)`
- C: `ocas_expr_integrate(...)`

Both return the unevaluated `Integral(...)` form when no closed form is
found, exactly like the Rust API.
