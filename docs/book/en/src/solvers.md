# Solvers

oCAS provides solvers for linear systems, Diophantine equations, and
polynomial systems. This chapter covers the available solvers and their usage.

---

## Linear systems over ℚ

`solve_linear_rational` solves an $n \times n$ system $Ax = b$ over the
rational numbers. Input coefficients are `i64` values; the solution is
returned as `(numerator, denominator)` pairs.

```rust
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
// x = [(2, 1), (1, 1)]  → 2, 1
```

Errors: `EmptySystem`, `NonSquare`, `Inconsistent`, `Underdetermined { rank }`.

Python:

```python
print(ocas.solve_linear_rational([[2, 1], [1, -1]], [5, 1]))
# [(2, 1), (1, 1)]
```

---

## Linear systems over ℤ

`solve_linear_integer` finds integer solutions to $Ax = b$. It returns an
error if no integer solution exists.

```rust
// 2x + y = 3
let a = vec![vec![2, 1]];
let b = vec![3];
let x = solve_linear_integer(&a, &b).unwrap();
// x = [1, 1]  (2·1 + 1·1 = 3)
```

Errors include `ResultNotInDomain` when the solution involves fractions.

---

## Diophantine equations

`solve_diophantine` solves the linear Diophantine equation
$a \cdot x + b \cdot y = c$ for integer $x, y$.

```rust
let sol = solve_diophantine(3, 5, 1).unwrap();
// sol = DiophantineSolution { x0: 2, y0: -1, x_step: 5, y_step: -3 }
```

The result gives a particular solution $(x_0, y_0)$ and step values.
The general solution is:

$$
\begin{aligned}
x &= x_0 + x_{step} \cdot t \\
y &= y_0 + y_{step} \cdot t
\end{aligned}
$$

for any integer $t$.

---

## Polynomial systems (via Gröbner bases)

`solve_polynomial_system` solves systems of polynomial equations: it computes
a Gröbner basis (default F4; the `Algorithm` enum selects Buchberger/F4/F5),
then performs Sturm real-root isolation and back-substitution for
zero-dimensional ideals.

```rust
use ocas_domain::{Rational, RationalDomain};
use ocas_poly::groebner::Algorithm;
use ocas_poly::ideal::solve_polynomial_system;
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};

let d = RationalDomain;
// x + y = 0, x*y - 1 = 0
let eq1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let eq2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);

let sol = solve_polynomial_system(&[eq1, eq2], Algorithm::Auto);
```

The result is a `PolynomialSystemSolution`: a list of real solutions
(`ZeroDimSolutions`) for zero-dimensional ideals, the Gröbner basis
(`PositiveDimensional`) for positive-dimensional ones, or `Empty` when no
solution exists.

### Ideal operations (since 0.23.0)

Ideal arithmetic is provided by the `ocas_poly::ideal` and
`ocas_poly::groebner` modules:

| Operation | Function | Description |
|---|---|---|
| Membership | `ideal_contains(gens, f, algo)` | Test if $f \in I$ |
| Sum | `ideal_sum(I, J)` | $I + J$ |
| Product | `ideal_product(I, J)` | $I \cdot J$ |
| Quotient | `ideal_quotient(I, J)` | $I : J$ |
| Saturation | `ideal_saturate(I, J)` | $I : J^\infty$ |
| Intersection | `ideal_intersection(I, J)` | $I \cap J$ |
| Elimination | `eliminate(gens, elim_vars, algo)` (`ocas_poly::groebner`) | Eliminate variables |
| Radical | `ideal_radical(gens)` | $\sqrt{I}$ |
| Primary decomposition | `primary_decomposition(gens)` | Decompose into primary components |
| Hilbert series | `hilbert_series(&gb)` (`ocas_poly::groebner::hilbert`) | Compute Hilbert series |

```rust
use ocas_domain::{Rational, RationalDomain};
use ocas_poly::groebner::Algorithm;
use ocas_poly::ideal::*;
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};

let d = RationalDomain;
let circle = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 2], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);
let line = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);

// Solve circle ∩ line
let sol = solve_polynomial_system(&[circle, line], Algorithm::Auto);
```

---

## Ordinary differential equations

`dsolve` solves ordinary differential equations analytically. The ODE is
given as an expression equal to zero, the unknown function (e.g. `y(x)`),
and the independent variable (e.g. `x`).

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);
// y' - y = 0
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
let sol = dsolve(&ctx, ode, None);
// ODESolution::Explicit(C1*exp(x))
```

`classify_ode` returns the applicable solving methods without solving:

```rust
let types = classify_ode(&ctx, ode);
// [LinearFirst, Separable, Homogeneous, PowerSeries]
```

### Supported ODE types

| Type | Form | Method |
|---|---|---|
| Separable | $f(x)dx = g(y)dy$ | Direct integration |
| Linear first-order | $y' + p(x)y = q(x)$ | Integrating factor $\mu = e^{\int p}$ |
| Bernoulli | $y' + p(x)y = q(x)y^n$ | Substitution $v = y^{1-n}$ |
| Exact | $M dx + N dy = 0$, $\partial M/\partial y = \partial N/\partial x$ | Potential function |
| Homogeneous | $y' = f(y/x)$ | Substitution $v = y/x$ |
| Integrating factors | non-exact first-order | $\mu(x)$ or $\mu(y)$ detection |
| Constant coefficients | $ay'' + by' + cy = f(x)$ | Characteristic equation |
| Cauchy-Euler | $ax^2y'' + bxy' + cy = f(x)$ | Indicial equation |
| Reduction of order | second-order linear | $y_2 = y_1\int e^{-\int p}/y_1^2$ |
| Variation of parameters | second-order non-homogeneous | Wronskian formula |
| Undetermined coefficients | polynomial/exponential/trigonometric forcing | Coefficient matching + resonance |
| Power series | ordinary point | Coefficient recursion |
| Frobenius | regular singular point | Indicial equation + recursion |
| Laplace IVP | first/second-order linear IVP | `dsolve_ivp` |
| 2×2 systems | $\mathbf{Y}' = A\mathbf{Y}$ | `dsolve_system` (eigen-decomposition) |

ODEs that cannot be solved analytically are returned as unevaluated
`ODE(equation, func)` forms.

### Initial value problems

`dsolve_ivp` solves linear constant-coefficient IVPs via the Laplace
transform:

```rust
// y' - y = 0, y(0) = 1  =>  y = exp(x)
let sol = dsolve_ivp(&ctx, ode, ctx.num(1), None);
```

### Systems

`dsolve_system` solves 2×2 constant-coefficient systems $\mathbf{Y}' = A\mathbf{Y}$:

```rust
// y1' = y2, y2' = -y1 (harmonic oscillator)
let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
// ODESolution::System([C1*sin(x) + C2*cos(x), C1*cos(x) - C2*sin(x)])
```

### Python

```python
import ocas

e = ocas.Expression("Derivative(y(x), x) - y(x)")
print(ocas.classify_ode(e, "y", "x"))     # ['LinearFirst', 'Separable', ...]
print(ocas.dsolve(e, "y", "x"))            # y = C1*exp(x)
print(ocas.dsolve_ivp(e, "y", "x", "1"))   # y = exp(x)
```

### C

```c
int err = 0;
char *types = ocas_ode_classify("Derivative(y(x), x) - y(x)", "y", "x", &err);
char *sol = ocas_ode_dsolve("Derivative(y(x), x) - y(x)", "y", "x", NULL, &err);
char *ivp = ocas_ode_dsolve_ivp("Derivative(y(x), x) - y(x)", "y", "x", "1", NULL, &err);
ocas_string_free(types);
ocas_string_free(sol);
ocas_string_free(ivp);
```

---

## Errors

The linear solvers return `Result<T, SolveError>`; `solve_diophantine`
returns `Option<DiophantineSolution>` (`None` when no solution exists); and
`solve_polynomial_system` returns `PolynomialSystemSolution`. Common error
variants:

| Error | Meaning |
|---|---|
| `EmptySystem` | No equations provided |
| `NonLinear` | System is not linear in the requested variables |
| `NonSquare` | Number of equations ≠ number of unknowns |
| `Inconsistent` | No solution exists |
| `Underdetermined { rank }` | Infinitely many solutions |
| `ResultNotInDomain` | Solution contains fractions when integers are required |

---

## See also

- [Rust API](./api/rust.md) — domain types and polynomial operations
- [Rewrite & Simplification](./rewrite.md) — simplifying solved expressions
- [Performance](./performance.md) — Gröbner basis benchmark results
