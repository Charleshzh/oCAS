# Rust API Reference: Solvers

This chapter covers the equation solving functionality in oCAS, organized into three modules: linear systems and Diophantine equations (`ocas-calc::solve`), polynomial system solving (`ocas-poly::ideal`), and ordinary differential equation solving (`ocas-calc::ode`).

---

## Linear systems

### solve_linear_rational

**Signature**: `pub fn solve_linear_rational(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<(i64, i64)>, SolveError>`

**Description**: Solves the linear system $A\mathbf{x} = \mathbf{b}$ over the field of rational numbers $\mathbb{Q}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `&[Vec<i64>]` | The $n \times n$ coefficient matrix stored row-major; each `Vec<i64>` is one row |
| `b` | `&[i64]` | The right-hand side vector of length $n$ |

**Returns**: `Result<Vec<(i64, i64)>, SolveError>`
- `Ok(vec)`: the solution vector; each element `(num, den)` represents the rational number $\frac{\text{num}}{\text{den}}$, with $\text{den} > 0$ and $\gcd(|\text{num}|, \text{den}) = 1$ guaranteed
- `Err(SolveError)`: see the error list below

**Errors**:
- `SolveError::EmptySystem` — $A$ or $\mathbf{b}$ is empty
- `SolveError::NonSquare` — $A$ is not square (or the length of $\mathbf{b}$ does not match the number of rows)
- `SolveError::Inconsistent` — the system has no solution (rank of the augmented matrix > rank of the coefficient matrix)
- `SolveError::Underdetermined { rank }` — the system has infinitely many solutions (rank deficient)
- `SolveError::ResultNotInDomain` — the solution is not in the target domain
- `SolveError::Matrix(_)` — an internal matrix operation failed (converted via `From<MatrixError>`)
- `SolveError::Other(_)` — any other error

**Example**:

```rust
use ocas_calc::solve::solve_linear_rational;

// Solve:
//   2x + y = 5
//   x - y = 1
// → x = 2, y = 1
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
assert_eq!(x, vec![(2, 1), (1, 1)]);
```

**See also**: [`solve_linear_integer`](#solve_linear_integer), [`Matrix::solve`](./rust-matrix.md#solve)

---

### solve_linear_integer

**Signature**: `pub fn solve_linear_integer(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<i64>, SolveError>`

**Description**: Solves the linear system $A\mathbf{x} = \mathbf{b}$ over the ring of integers $\mathbb{Z}$. Requires the solution to be integral; otherwise an error is returned.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `&[Vec<i64>]` | The $n \times n$ coefficient matrix stored row-major |
| `b` | `&[i64]` | The right-hand side vector of length $n$ |

**Returns**: `Result<Vec<i64>, SolveError>`
- `Ok(vec)`: the integral solution vector
- `Err(SolveError)`: the error types are the same as for [`solve_linear_rational`](#solve_linear_rational), with a possible additional `ResultNotInDomain` (the solution contains fractions)

**Errors**: same as [`solve_linear_rational`](#solve_linear_rational).

**Example**:

```rust
use ocas_calc::solve::solve_linear_integer;

// Solve:
//   x + y = 3
//   x - y = 1
// → x = 2, y = 1
let a = vec![vec![1, 1], vec![1, -1]];
let b = vec![3, 1];
let x = solve_linear_integer(&a, &b).unwrap();
assert_eq!(x, vec![2, 1]);
```

**See also**: [`solve_linear_rational`](#solve_linear_rational)

---

## Diophantine equations

### DiophantineSolution

**Signature**:
```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiophantineSolution {
    /// Particular solution (x0, y0)
    pub particular: (i64, i64),
    /// Homogeneous solution step (x_step, y_step)
    pub general: (i64, i64),
}
```

**Description**: Representation of the solutions of the linear Diophantine equation $ax + by = c$.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `particular` | `(i64, i64)` | A particular solution $(x_0, y_0)$ |
| `general` | `(i64, i64)` | A basis solution $(x_{\text{step}}, y_{\text{step}})$ of the homogeneous equation $ax + by = 0$; the general solution is $(x_0 + k \cdot x_{\text{step}},\; y_0 + k \cdot y_{\text{step}})$ with $k \in \mathbb{Z}$ |

---

### solve_diophantine

**Signature**: `pub fn solve_diophantine(a: i64, b: i64, c: i64) -> Option<DiophantineSolution>`

**Description**: Solves the linear Diophantine equation $ax + by = c$ using the extended Euclidean algorithm. A solution exists if and only if $\gcd(a, b) \mid c$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `i64` | The coefficient of $x$ |
| `b` | `i64` | The coefficient of $y$ |
| `c` | `i64` | The constant on the right-hand side |

**Returns**: `Option<DiophantineSolution>`
- `Some(sol)`: the particular/general solution structure when a solution exists
- `None`: $\gcd(a, b) \nmid c$; the equation has no integer solution

**Errors**: does not return a `Result`; returns `None` when unsolvable.

**Example**:

```rust
use ocas_calc::solve::solve_diophantine;

// Solve 3x + 5y = 1
// Particular solution: x = 2, y = -1 (since 3×2 + 5×(-1) = 1)
// General solution: x = 2 + 5k, y = -1 - 3k
let sol = solve_diophantine(3, 5, 1).unwrap();
assert_eq!(sol.particular, (2, -1));
assert_eq!(sol.general, (5, -3));

// No solution: 2x + 4y = 3 (gcd(2,4)=2 does not divide 3)
assert!(solve_diophantine(2, 4, 3).is_none());
```

**See also**: [`solve_linear_integer`](#solve_linear_integer)

---

## Polynomial system solving

### RealSolution

**Signature**:
```rust
#[derive(Debug, Clone)]
pub struct RealSolution {
    /// The value of each variable (in variable order)
    pub values: Vec<f64>,
    /// Algebraic multiplicity
    pub multiplicity: usize,
}
```

**Description**: A real solution (numerical approximation) of a polynomial system of equations.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `values` | `Vec<f64>` | The variable values $[x_1, x_2, \dots, x_n]$, in the order of the variables in the original system |
| `multiplicity` | `usize` | The algebraic multiplicity of the solution (the current implementation always sets this to 1; the actual multiplicity is not computed yet) |

---

### ZeroDimSolutions

**Signature**:
```rust
#[derive(Debug, Clone)]
pub struct ZeroDimSolutions {
    /// All real solutions
    pub solutions: Vec<RealSolution>,
    /// Vector space dimension of the quotient ring k[x₁,...,xₙ]/I
    pub vector_space_dimension: usize,
}
```

**Description**: The solution result of a zero-dimensional polynomial system.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `solutions` | `Vec<RealSolution>` | All real solutions (may contain approximations) |
| `vector_space_dimension` | `usize` | The vector space dimension of the quotient ring $k[x_1, \dots, x_n]/I$, equal to the product of the degrees of the univariate polynomials in the Lex GB |

---

### PolynomialSystemSolution

**Signature**:
```rust
#[derive(Debug, Clone)]
pub enum PolynomialSystemSolution {
    /// Zero-dimensional system: finitely many solutions
    ZeroDimensional(ZeroDimSolutions),
    /// Positive-dimensional system: the solution set has a positive-dimensional component; returns a Gröbner basis
    PositiveDimensional(GroebnerBasis<RationalDomain, Lex>),
    /// Inconsistent system (the ideal is ⟨1⟩), no solutions
    Empty,
}
```

**Description**: The enum type of the result of solving a polynomial system. The system automatically determines the dimension and selects the appropriate solving strategy.

**Variants**:

| Variant | Description |
|---|---|
| `ZeroDimensional(z)` | The system has finitely many solutions; `z` contains all real solutions and the vector space dimension |
| `PositiveDimensional(gb)` | The solution set of the system has a positive-dimensional component (infinitely many solutions); returns a Lex-order Gröbner basis |
| `Empty` | The system is inconsistent (the Gröbner basis is {1}); no solutions. Note: an empty input equation list returns `PositiveDimensional` (with an empty basis), not `Empty` |

---

### solve_polynomial_system

**Signature**:
```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

**Description**: Solves the polynomial system $f_1 = f_2 = \cdots = f_m = 0$. Automatically computes a Lex-order Gröbner basis, determines the dimension of the system, and for zero-dimensional systems finds all real solutions by triangular decomposition and back substitution.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equations` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The polynomial system; each polynomial represents one equation $f_i = 0$ |
| `algo` | `Algorithm` | Gröbner basis algorithm selection: `Algorithm::Auto` (automatic), `Algorithm::Buchberger`, `Algorithm::F4`, `Algorithm::F5` |

**Returns**: [`PolynomialSystemSolution`](#polynomialsystemsolution)

**Errors**: does not return a `Result`; all possible outcomes are expressed through the enum variants (zero-dimensional solutions, positive-dimensional basis, empty system).

**Algorithm description**:
1. Compute a Gröbner basis of the ideal generated by the equations (using the specified algorithm). Since the input polynomials are already typed in Lex order (`SparseMultivariatePolynomial<RationalDomain, Lex>`), the GB is computed directly in Lex order; no order change is needed
2. If the GB is {1}, return `Empty` (the system is inconsistent, no solutions)
3. Determine whether it is zero-dimensional: whether each variable $x_i$ has a pure-power leading monomial $x_i^{N_i}$ in the GB
4. **Zero-dimensional**: extract the triangular structure of the Lex GB and back-substitute variable by variable starting from the last one (Sturm root isolation + bisection refinement), solving the real roots of the univariate polynomial for each variable; the quotient-ring dimension is estimated as the product of the univariate polynomial degrees
5. **Positive-dimensional**: return the Lex GB for further analysis

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};
use ocas_poly::ideal::{solve_polynomial_system, PolynomialSystemSolution};

let d = RationalDomain;

// System:
//   x² + y² - 1 = 0  (unit circle)
//   x - y = 0        (line y = x)
// Solutions: (±1/√2, ±1/√2)
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),   // x²
    (vec![0, 2], Rational::new(1, 1)),   // y²
    (vec![0, 0], Rational::new(-1, 1)),  // -1
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),   // x
    (vec![0, 1], Rational::new(-1, 1)),  // -y
]);

let sol = solve_polynomial_system(&[f1, f2], Algorithm::Auto);
match sol {
    PolynomialSystemSolution::ZeroDimensional(z) => {
        assert_eq!(z.solutions.len(), 2);
        // Two solutions: approximately (0.707, 0.707) and (-0.707, -0.707)
    }
    _ => panic!("expected zero-dimensional system"),
}
```

**See also**: [`GroebnerBasis`](./rust-groebner.md), [`Algorithm`](./rust-groebner.md#algorithm), [`is_zero_dimensional`](./rust-groebner.md#is_zero_dimensional)

---

## Ordinary differential equation (ODE) solving

### ODE

**Signature**:
```rust
#[derive(Debug, Clone, Copy)]
pub struct ODE<'a> {
    /// Equation in standard form lhs - rhs = 0
    pub equation: Atom<'a>,
    /// Unknown function, e.g. y(x)
    pub func: Atom<'a>,
    /// Independent variable, e.g. x
    pub var: Symbol,
}
```

**Description**: Describes an ordinary differential equation. The equation is stored in the standard form `equation = 0`, where `func` is the unknown function (e.g., $y(x)$) and `var` is the independent variable (e.g., $x$).

**Fields**:

| Field | Type | Description |
|---|---|---|
| `equation` | `Atom<'a>` | The normalized form of the equation (`lhs - rhs`) |
| `func` | `Atom<'a>` | The unknown function, e.g., `y(x)` |
| `var` | `Symbol` | The independent variable symbol, e.g., `"x"` |

**Construction example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::ODE;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);

// Represents y' - y = 0
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
```

---

### ODESolution

**Signature**:
```rust
#[derive(Debug, Clone, Copy)]
pub enum ODESolution<'a> {
    /// Explicit solution y = expr
    Explicit(Atom<'a>),
    /// Implicit solution F(x, y) = 0
    Implicit(Atom<'a>),
    /// Parametric solution (x(t), y(t))
    Parametric(Atom<'a>, Atom<'a>),
    /// Series solution (truncated expression + number of terms)
    Series(Atom<'a>, usize),
    /// System solution (components)
    System(&'a [Atom<'a>]),
    /// Failed to solve; returns the original ODE
    Unsolved(ODE<'a>),
}
```

**Description**: The enum type of the result of solving an ODE.

**Variants**:

| Variant | Description |
|---|---|
| `Explicit(expr)` | An explicit solution $y = \text{expr}$, containing free constants C1, C2, ... |
| `Implicit(expr)` | An implicit solution $F(x, y) = 0$ (stores the expression $F$) |
| `Parametric(x(t), y(t))` | A parametric solution: two expressions giving $x(t)$ and $y(t)$ |
| `Series(expr, n)` | A series solution; `expr` is the truncated series expression and `n` the number of terms |
| `System(components)` | A system solution; `components[i]` is the solution for the $i$-th unknown function |
| `Unsolved(ode)` | No analytic solution found; returns the original ODE |

---

### ODEType

**Signature**:
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ODEType {
    Separable,
    LinearFirst,
    Bernoulli,
    Exact,
    Homogeneous,
    LinearConstantCoeff,
    CauchyEuler,
    ReductionOfOrder,
    PowerSeries,
}
```

**Description**: The ODE type enum, which determines which solution method is used.

**Variants**:

| Variant | Equation form | Description |
|---|---|---|
| `Separable` | $g(y)\,y' = f(x)$ | Separable: $\int g(y)\,dy = \int f(x)\,dx + C$ |
| `LinearFirst` | $y' + p(x)\,y = q(x)$ | First-order linear: integrating factor $\mu = e^{\int p\,dx}$ |
| `Bernoulli` | $y' + p(x)\,y = q(x)\,y^n$ | Bernoulli equation: substitution $v = y^{1-n}$ reduces to linear |
| `Exact` | $M\,dx + N\,dy = 0$, $\frac{\partial M}{\partial y} = \frac{\partial N}{\partial x}$ | Exact equation: potential function $F(x,y) = C$ |
| `Homogeneous` | $y' = f(y/x)$ | Homogeneous equation: substitution $v = y/x$ reduces to separable |
| `LinearConstantCoeff` | $a\,y'' + b\,y' + c\,y = f(x)$ | Linear with constant coefficients: characteristic equation + undetermined coefficients (polynomial/exponential/trigonometric forcing) |
| `CauchyEuler` | $a\,x^2\,y'' + b\,x\,y' + c\,y = f(x)$ | Cauchy-Euler: substitution $x = e^t$ reduces to constant coefficients |
| `ReductionOfOrder` | $a(x)\,y'' + b(x)\,y' + c(x)\,y = f(x)$ | Reduction of order: tries simple candidate solutions $y_1$; $y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$; nonzero forcing uses variation of parameters for the particular solution |
| `PowerSeries` | linear ODE, ordinary point $x_0 = 0$ | Power series: $y = \sum_{n=0}^{N-1} a_n x^n$; falls back to the Frobenius method automatically on failure |

> Note: variation of parameters and the Frobenius method are used internally as solving routines (the former for particular solutions in reduction of order, the latter as the fallback when the power series fails), but they are **not** separate `ODEType` variants.

---

### classify_ode

**Signature**: `pub fn classify_ode<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>) -> Vec<ODEType>`

**Description**: Analyzes the ODE and returns all applicable solution methods, ordered by priority. The order actually returned (implemented in `classify.rs`): for first-order equations, `LinearFirst` → `Bernoulli` → `Separable` → `Exact` → `Homogeneous` are checked in that order; for second-order (and higher) linear equations, `LinearConstantCoeff` → `CauchyEuler` → `ReductionOfOrder` (second order only); finally, `PowerSeries` is appended as a fallback for any linear equation of order ≥ 1.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `ode` | `ODE<'a>` | The ODE to classify |

**Returns**: `Vec<ODEType>` — the list of all applicable solution methods (possibly empty).

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_calc::ode::{classify_ode, ODE, ODEType};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);

// y' - y = 0 is first-order linear
let dy = ctx.fun("Derivative", &[y, x]);
let eq = ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]);
let ode = ODE { equation: eq, func: y, var: Symbol::new("x") };
let methods = classify_ode(&ctx, ode);
assert!(methods.contains(&ODEType::LinearFirst));
```

**See also**: [`dsolve`](#dsolve), [`ODEType`](#odetype)

---

### dsolve

**Signature**: `pub fn dsolve<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>, hint: Option<ODEType>) -> ODESolution<'a>`

**Description**: Solves an ordinary differential equation. Automatically classifies the ODE type and tries the methods in priority order; a specific method can also be requested via `hint`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `ode` | `ODE<'a>` | The ODE to solve |
| `hint` | `Option<ODEType>` | Optional: a specific solution method. When `None`, the ODE is classified automatically and the methods are tried in priority order |

**Returns**: [`ODESolution<'a>`](#odesolution)

**Solution strategy** (when `hint = None`):
1. Call [`classify_ode`](#classify_ode) to obtain the list of candidate methods
2. Try each method in priority order
3. The solution of the first successful method is returned
4. `ODESolution::Unsolved` is returned when all methods fail

**First-order solvers**:
- **Separable**: $\int g(y)\,dy = \int f(x)\,dx + C$
- **First-order linear**: integrating factor $\mu(x) = e^{\int p(x)\,dx}$, general solution $y = \frac{1}{\mu}\left(\int \mu\,q\,dx + C\right)$
- **Bernoulli**: substitution $v = y^{1-n}$ reduces to first-order linear
- **Exact equations**: find a potential function $F$ such that $\frac{\partial F}{\partial x} = M$ and $\frac{\partial F}{\partial y} = N$; when not exact, an integrating factor $\mu(x)$ or $\mu(y)$ is attempted
- **Homogeneous**: substitution $v = y/x$ reduces to separable

**Second-order solvers**:
- **Constant coefficients**: characteristic equation $ar^2 + br + c = 0$; build the fundamental solution set according to the discriminant $\Delta = b^2 - 4ac$; find a particular solution by undetermined coefficients
- **Cauchy-Euler**: substitution $x = e^t$ converts to a constant-coefficient equation; indicial equation $ar^2 + (b-a)r + c = 0$
- **Reduction of order**: tries simple candidate solutions ($1, x, x^2, e^x, e^{-x}, e^{2x}$); once one is found, constructs a second solution via $y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$
- **Variation of parameters**: given fundamental solutions $y_1, y_2$ of the homogeneous equation, finds a particular solution using the Wronskian $W = y_1 y_2' - y_1' y_2$

**Series solvers** (dispatched by `PowerSeries`, fixed expansion point $x_0 = 0$, 8 terms):
- **Power series**: expand $y = \sum_{n=0}^{N-1} a_n x^n$ at the ordinary point $x_0 = 0$; substitute into the ODE and recurse to determine the coefficients
- **Frobenius** (internal fallback): when the power series is not applicable at the ordinary point, falls back to the regular singular point method, expanding $y = x^r \sum_{n=0}^{N-1} a_n x^n$; the indicial equation determines $r$

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::{dsolve, ODE, ODESolution};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);

// Solve y' - y = 0 → C1*exp(x)
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
let sol = dsolve(&ctx, ode, None);
assert!(matches!(sol, ODESolution::Explicit(_)));

// Specify a method
let sol_hint = dsolve(&ctx, ode, Some(ODEType::LinearFirst));
assert!(matches!(sol_hint, ODESolution::Explicit(_)));
```

**See also**: [`classify_ode`](#classify_ode), [`dsolve_ivp`](#dsolve_ivp), [`dsolve_system`](#dsolve_system)

---

### dsolve_ivp

**Signature**:
```rust
pub fn dsolve_ivp<'a>(
    ctx: &'a AtomArena<'a>,
    ode: ODE<'a>,
    y0: Atom<'a>,
    y1: Option<Atom<'a>>,
) -> ODESolution<'a>
```

**Description**: Solves the initial value problem (IVP) for first- or second-order linear constant-coefficient ODEs using the Laplace transform. The result is an explicit solution without free constants.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `ode` | `ODE<'a>` | The ODE (must be linear with constant coefficients) |
| `y0` | `Atom<'a>` | The initial condition $y(0)$ |
| `y1` | `Option<Atom<'a>>` | $y'(0)$ for second-order problems; ignored for first-order problems |

**Returns**: [`ODESolution<'a>`](#odesolution)
- `Explicit(expr)`: an explicit solution without free constants
- `Unsolved`: the Laplace transform method is not applicable or the inverse cannot be computed

**Method description**:
1. Take the Laplace transform of both sides of the ODE, using $\mathcal{L}\{y'\} = sY - y(0)$ and $\mathcal{L}\{y''\} = s^2 Y - sy(0) - y'(0)$
2. Substitute the initial conditions to obtain an algebraic equation for $Y(s)$
3. Solve for $Y(s)$
4. Perform the inverse Laplace transform via partial fraction decomposition + lookup of standard pairs

**Supported forcing terms**: polynomials, exponentials $e^{kx}$, sines/cosines $\sin(\omega x)$/$\cos(\omega x)$, and linear combinations of $e^{kx}\sin(\omega x)$/$e^{kx}\cos(\omega x)$.

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::{dsolve_ivp, ODE, ODESolution};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);

// y' - y = 0, y(0) = 1  →  y = exp(x)
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
let sol = dsolve_ivp(&ctx, ode, ctx.num(1), None);
assert!(matches!(sol, ODESolution::Explicit(_)));
```

**See also**: [`dsolve`](#dsolve)

---

### dsolve_system

**Signature**:
```rust
pub fn dsolve_system<'a>(
    ctx: &'a AtomArena<'a>,
    equations: &[Atom<'a>],
    funcs: &[Atom<'a>],
    var: Symbol,
) -> ODESolution<'a>
```

**Description**: Solves a $2 \times 2$ constant-coefficient linear ODE system $\mathbf{Y}' = A\mathbf{Y}$. The general solution is obtained via eigenvalue decomposition.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `equations` | `&[Atom<'a>]` | The list of equations, each of the form `Derivative(y_i, x) - (a_i1*y1 + a_i2*y2) = 0` |
| `funcs` | `&[Atom<'a>]` | The list of unknown functions, e.g., `[y1(x), y2(x)]` |
| `var` | `Symbol` | The independent variable symbol |

**Returns**: [`ODESolution<'a>`](#odesolution)
- `System(&[Atom])`: the general solution of each component (containing free constants C1, C2)
- `Unsolved`: the system is not supported (not 2×2, or the eigenvalues have no closed form)

**Supported cases**:
- Distinct real eigenvalues $\lambda_1 \neq \lambda_2$: $\mathbf{Y} = C_1 \mathbf{v}_1 e^{\lambda_1 x} + C_2 \mathbf{v}_2 e^{\lambda_2 x}$
- Repeated real eigenvalues (with a generalized eigenvector): $\mathbf{Y} = (C_1 \mathbf{v} + C_2(\mathbf{w} + x\mathbf{v}))e^{\lambda x}$
- Conjugate complex eigenvalues $\alpha \pm \beta i$: real-valued fundamental solutions such as $e^{\alpha x}(\mathbf{p}\cos\beta x - \mathbf{q}\sin\beta x)$

**Example**:

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::{dsolve_system, ODESolution};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y1 = ctx.fun("y1", &[x]);
let y2 = ctx.fun("y2", &[x]);
let dy1 = ctx.fun("Derivative", &[y1, x]);
let dy2 = ctx.fun("Derivative", &[y2, x]);

// Harmonic oscillator: y1' = y2, y2' = -y1
let eq1 = ctx.add(&[dy1, ctx.mul(&[ctx.num(-1), y2])]);
let eq2 = ctx.add(&[dy2, y1]);
let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
assert!(matches!(sol, ODESolution::System(_)));
```

**See also**: [`dsolve`](#dsolve)

---

## SolveError

**Signature**:
```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolveError {
    /// The system of equations is empty
    EmptySystem,
    /// The system is nonlinear
    NonLinear,
    /// The coefficient matrix is not square (or the equation count does not match the unknown count)
    NonSquare,
    /// The system is inconsistent (no solution)
    Inconsistent,
    /// The system is underdetermined (infinitely many solutions)
    Underdetermined { rank: usize },
    /// The result is not in the target domain
    ResultNotInDomain,
    /// An internal matrix operation failed
    Matrix(MatrixError),
    /// Other error with a description
    Other(String),
}
```

**Description**: Errors that can occur during equation solving. Implements the `Display` and `Error` traits. Implements `From<MatrixError>`: `ShapeMismatch` → `NonSquare`, `Inconsistent` → `Inconsistent`, `Underdetermined { rank }` → `Underdetermined { rank }`, `ResultNotInDomain` → `ResultNotInDomain`, `RightHandSideIsNotVector` → `Other("right-hand side is not a vector")`.

**Variants**:

| Variant | Description |
|---|---|
| `EmptySystem` | The input system of equations is empty |
| `NonLinear` | The system is nonlinear in the target variables |
| `NonSquare` | The coefficient matrix is not square (equation count ≠ unknown count) |
| `Inconsistent` | The rank of the augmented matrix is greater than the rank of the coefficient matrix |
| `Underdetermined { rank }` | The rank of the coefficient matrix is less than the number of variables; `rank` is the actual rank |
| `ResultNotInDomain` | The solution is not in the target number domain (e.g., an integral solution was requested but the solution contains fractions) |
| `Matrix(MatrixError)` | An internal matrix operation failed (converted automatically via `From<MatrixError>`) |
| `Other(String)` | Other error with a description |

**Example**:

```rust
use ocas_calc::solve::{solve_linear_rational, SolveError};

// Inconsistent system: x + y = 1, x + y = 2
let a = vec![vec![1, 1], vec![1, 1]];
let b = vec![1, 2];
let result = solve_linear_rational(&a, &b);
assert_eq!(result, Err(SolveError::Inconsistent));
```

---

## Helper functions

### substitute_solution_collected

**Signature**:
```rust
pub fn substitute_solution_collected<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a>
```

**Description**: Substitutes a candidate solution into the ODE equation and returns the residual after collecting like terms. A zero residual means the candidate solution satisfies the ODE. Mainly used for verification.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `equation` | `Atom<'a>` | The ODE equation (in `lhs - rhs` form) |
| `func` | `Atom<'a>` | The unknown function (e.g., `y(x)`) |
| `sol` | `Atom<'a>` | The candidate solution expression |
| `var` | `Symbol` | The independent variable symbol |

**Returns**: `Atom<'a>` — the residual expression after substitution and simplification. If zero, the candidate solution is correct.

---

## Source modules

| Module | Path | Contents |
|---|---|---|
| Linear solving | `ocas-calc/src/solve.rs` | `solve_linear_rational`, `solve_linear_integer`, `solve_diophantine`, `SolveError` |
| Polynomial systems | `ocas-poly/src/ideal.rs` | `solve_polynomial_system`, `PolynomialSystemSolution`, `RealSolution`, `ZeroDimSolutions` |
| ODE classification | `ocas-calc/src/ode/classify.rs` | `classify_ode`, `ODEType` |
| First-order ODEs | `ocas-calc/src/ode/first_order.rs` | separable, linear, Bernoulli, exact, and homogeneous solvers |
| Second-order ODEs | `ocas-calc/src/ode/second_order.rs` | constant-coefficient, Cauchy-Euler, reduction-of-order, and variation-of-parameters solvers |
| Series solutions for ODEs | `ocas-calc/src/ode/series.rs` | power series and Frobenius solvers |
| Laplace transforms for ODEs | `ocas-calc/src/ode/laplace.rs` | the Laplace transform implementation behind `dsolve_ivp` |
| ODE systems | `ocas-calc/src/ode/systems.rs` | the eigenvalue decomposition implementation behind `dsolve_system` |
