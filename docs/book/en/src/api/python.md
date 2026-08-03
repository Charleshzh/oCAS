# Python API Reference

This chapter documents the complete API of the oCAS Python bindings (the `ocas` module). Every class and function includes its signature, parameter descriptions, return value, exceptions, and complete examples.

> **Version**: This document corresponds to oCAS 0.24.x.

## Import

```python
import ocas
```

The module exports all classes and functions at the top level; no submodule imports are required.

---

## Expressions

### Expression

Symbolic expression class. Constructed by parsing a string; supports arithmetic operations, calculus, simplification, and more.

**Signature**:

```python
ocas.Expression(input: str)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `input` | `str` | Symbolic expression string; `^` denotes exponentiation and `*` denotes multiplication |

**Returns**: an `Expression` instance.

**Exceptions**: `ValueError` — the expression string failed to parse.

**Supported operations**:

| Operation | Method | Description |
|---|---|---|
| Addition | `__add__` | `e1 + e2` |
| Subtraction | `__sub__` | `e1 - e2` |
| Multiplication | `__mul__` | `e1 * e2` |
| Power | `__pow__` | `e1 ** e2` |
| Negation | `__neg__` | `-e` |
| Equality | `__eq__` | `e1 == e2` (comparison based on the normalized form) |
| Hash | `__hash__` | Usable as a dictionary key or set element |
| String | `__str__` | Outputs the internal representation |
| Representation | `__repr__` | `Expression("...")` format |

**Methods**:

---

#### Expression.clone

```python
Expression.clone() -> Expression
```

Returns a deep copy of the expression.

**Example**:

```python
>>> e = ocas.Expression("x + 1")
>>> f = e.clone()
>>> e == f
True
```

---

#### Expression.simplify

```python
Expression.simplify() -> Expression
```

Simplifies the expression using the built-in rule set (fixed-point iteration).

**Returns**: the simplified `Expression`.

**Exceptions**: `ValueError` — an internal error occurred during simplification.

**Example**:

```python
>>> e = ocas.Expression("x + x")
>>> print(e.simplify())
2*x
```

---

#### Expression.diff

```python
Expression.diff(var: str) -> Expression
```

Differentiates the expression with respect to the variable `var`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var` | `str` | The differentiation variable |

**Returns**: the derivative expression.

**Example**:

```python
>>> e = ocas.Expression("x^3 + 2*x")
>>> print(e.diff("x"))
2 + (3*(x^2))
```

---

#### Expression.integrate

```python
Expression.integrate(var: str) -> Expression
```

Symbolically integrates the expression with respect to the variable `var`. Uses a layered integration pipeline: fast lookup table → rational functions → Risch → trigonometric rewrite + realify → special functions → heuristic → unevaluated form.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var` | `str` | The integration variable |

**Returns**: the resulting expression. If no closed form can be found, returns the unevaluated `Integral(expr, var)` form.

**Example**:

```python
>>> e = ocas.Expression("x^2")
>>> print(e.integrate("x"))
(3^-1)*(x^3)
```

---

#### Expression.integrate_heuristic

```python
Expression.integrate_heuristic(var: str) -> Expression
```

Integrates using heuristic methods only (integration by parts, trigonometric substitutions, etc.), without invoking the full Risch algorithm. Returns a closed form on success and an unevaluated form on failure.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var` | `str` | The integration variable |

**Returns**: the integral result or the unevaluated `Integral(expr, var)`.

**Example**:

```python
>>> e = ocas.Expression("sin(x)")
>>> print(e.integrate_heuristic("x"))
-1*(cos(x))
```

---

#### Expression.taylor

```python
Expression.taylor(var: str, point: Expression, order: int) -> Expression
```

Computes the Taylor expansion of the expression in `var` about `point` up to order `order`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var` | `str` | The expansion variable |
| `point` | `Expression` | The expansion center |
| `order` | `int` | Expansion order (a non-negative integer) |

**Returns**: the Taylor polynomial.

**Example**:

```python
>>> e = ocas.Expression("exp(x)")
>>> p = ocas.Expression("0")
>>> print(e.taylor("x", p, 4))
1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3)) + ((24^-1)*(x^4))
```

---

#### Expression.substitute

```python
Expression.substitute(var: str, replacement: Expression) -> Expression
```

Replaces every occurrence of `var` in the expression with `replacement`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var` | `str` | The variable to be replaced |
| `replacement` | `Expression` | The replacement expression |

**Returns**: the substituted expression.

**Example**:

```python
>>> e = ocas.Expression("x^2 + y")
>>> r = ocas.Expression("3")
>>> print(e.substitute("x", r))
y + (3^2)
```

---

## Polynomials

### Polynomial

Dense univariate polynomial supporting three coefficient domains: integers ℤ, rationals ℚ, and finite fields GF(p).

**Signature**:

```python
ocas.Polynomial(coeffs, domain=None)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `coeffs` | `list` | Coefficient list, ordered from lowest to highest degree |
| `domain` | `str`, `FiniteField`, or `None` | Coefficient domain: `"integer"` (default), `"rational"`, or `FiniteField(p)` |

**Coefficient formats**:
- Integer domain: a list of integers, e.g. `[1, 2, 1]` represents $1 + 2x + x^2$
- Rational domain: a list of integers or a list of `(numerator, denominator)` tuples
- Finite field: a list of integers; arithmetic is automatically reduced modulo p

**Exceptions**: `TypeError` — coefficient types do not match the selected domain.

**Example**:

```python
>>> from ocas import Polynomial
>>> p = Polynomial([1, 2, 1])
>>> print(p.degree())
2
```

---

#### Polynomial.coeffs

```python
Polynomial.coeffs() -> list[str]
```

Returns the coefficient list (as strings), ordered from lowest to highest degree.

**Example**:

```python
>>> p = Polynomial([1, 2, 1])
>>> p.coeffs()
['1', '2', '1']
```

---

#### Polynomial.degree

```python
Polynomial.degree() -> int | None
```

Returns the degree of the polynomial. The zero polynomial returns `None`.

**Example**:

```python
>>> Polynomial([1, 2, 1]).degree()
2
>>> Polynomial([0]).degree() is None
True
```

---

#### Polynomial.len

```python
Polynomial.len() -> int
```

Returns the number of stored coefficients.

---

#### Polynomial.is_zero

```python
Polynomial.is_zero() -> bool
```

Whether the polynomial is the zero polynomial.

---

#### Polynomial.eval

```python
Polynomial.eval(x) -> str
```

Evaluates the polynomial at the point `x`. Over the integer domain `x` is an integer; over the rational domain `x` is an integer or a `(numerator, denominator)` tuple.

**Returns**: the string representation of the evaluation result.

**Example**:

```python
>>> p = Polynomial([1, 2, 1])
>>> p.eval(2)
'9'
```

---

#### Polynomial.derivative

```python
Polynomial.derivative() -> Polynomial
```

Returns the formal derivative.

**Example**:

```python
>>> p = Polynomial([1, 2, 1])
>>> q = p.derivative()
>>> q.coeffs()
['2', '2']
```

---

#### Polynomial.integral

```python
Polynomial.integral() -> Polynomial
```

Returns the formal integral (with zero constant term). Available on all three coefficient domains: exact over the rationals; over the integers and finite fields, coefficients whose division by $i+1$ does not divide exactly are set to zero (no exception is raised).

---

#### Polynomial.primitive_part

```python
Polynomial.primitive_part() -> Polynomial
```

Returns the primitive part (content removed). Supported on the integer domain only.

**Exceptions**: `ValueError` — the polynomial is not over the integer domain.

---

#### Polynomial.factor

```python
Polynomial.factor() -> list[PolynomialFactor]
```

Complete factorization. Returns a list of `(factor, multiplicity)` pairs: over the integers each factor is primitive; over a finite field each factor is monic.

**Factorization strategy**:
- ℤ[x]: square-free factorization + Berlekamp–Zassenhaus + Hensel lifting
- GF(p)[x]: Berlekamp's algorithm

**Example**:

```python
>>> p = Polynomial([1, 2, 1])
>>> for f in p.factor():
...     print(f.factor.coeffs(), f.multiplicity)
['1', '1'] 2
```

---

#### Polynomial.square_free_factorization

```python
Polynomial.square_free_factorization() -> list[PolynomialFactor]
```

Square-free factorization. Returns the pairwise-distinct square-free factors together with their multiplicities.

**Example**:

```python
>>> p = Polynomial([1, 2, 1])
>>> for f in p.square_free_factorization():
...     print(f.factor.coeffs(), f.multiplicity)
['1', '1'] 2
```

---

#### Polynomial.is_square_free

```python
Polynomial.is_square_free() -> bool
```

Whether the polynomial is square-free (has no repeated irreducible factors).

---

#### Polynomial.gcd

```python
Polynomial.gcd(other: Polynomial) -> Polynomial
```

Computes the greatest common divisor with `other`. Both polynomials must be over the same coefficient domain.

**Exceptions**: `TypeError` — the coefficient domains do not match.

**Example**:

```python
>>> a = Polynomial([0, 1])         # x
>>> b = Polynomial([0, 0, 1])     # x^2
>>> a.gcd(b).coeffs()
['0', '1']
```

---

#### Polynomial.div_rem

```python
Polynomial.div_rem(other: Polynomial) -> tuple[Polynomial, Polynomial] | None
```

Division with remainder, returning `(quotient, remainder)`. Returns `None` when `other` is zero.

**Exceptions**: `TypeError` — the coefficient domains do not match.

**Example**:

```python
>>> a = Polynomial([1, 2, 1])
>>> b = Polynomial([1, 1])
>>> q, r = a.div_rem(b)
>>> q.coeffs()
['1', '1']
>>> r.is_zero()
True
```

---

#### Polynomial arithmetic

| Operation | Syntax | Description |
|---|---|---|
| Addition | `p + q` | Polynomial addition |
| Subtraction | `p - q` | Polynomial subtraction |
| Multiplication | `p * q` | Polynomial multiplication |
| Negation | `-p` | Negates the coefficients |
| Equality | `p == q` | Comparison of normalized coefficients |

---

### PolynomialFactor

A single factor in a factorization result.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `factor` | `Polynomial` | The factor polynomial (primitive over ℤ, monic over GF(p)) |
| `multiplicity` | `int` | The multiplicity |

---

### MultivariatePolynomial

Sparse multivariate polynomial ($\mathbb{Q}[x_1, \dots, x_n]$), using the lexicographic order (Lex).

**Signature**:

```python
ocas.MultivariatePolynomial(terms: dict, n_vars: int)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `terms` | `dict[tuple[int,...], int]` | Keys are exponent tuples (of length `n_vars`), values are coefficients |
| `n_vars` | `int` | The number of variables |

**Example**:

```python
>>> from ocas import MultivariatePolynomial
>>> # x² + y² - 1
>>> p = MultivariatePolynomial({(2, 0): 1, (0, 2): 1, (0, 0): -1}, n_vars=2)
>>> p.n_vars()
2
```

---

#### MultivariatePolynomial.n_vars

```python
MultivariatePolynomial.n_vars() -> int
```

Returns the number of variables.

---

## Matrices

### Matrix

Dense matrix over a field, supporting ℤ, ℚ, and GF(p).

**Signature**:

```python
ocas.Matrix(rows, domain=None)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rows` | `list[list]` | Row-major two-dimensional list |
| `domain` | `str`, `FiniteField`, or `None` | Coefficient domain: `"integer"` (default), `"rational"`, or `FiniteField(p)` |

**Exceptions**: `TypeError` — element types do not match the selected domain.

**Example**:

```python
>>> from ocas import Matrix
>>> a = Matrix([[1, 2], [3, 4]])
>>> a.determinant()
'-2'
```

---

#### Matrix.nrows

```python
Matrix.nrows -> int  # property
```

The number of rows.

---

#### Matrix.ncols

```python
Matrix.ncols -> int  # property
```

The number of columns.

---

#### Matrix.shape

```python
Matrix.shape() -> tuple[int, int]
```

Returns `(nrows, ncols)`.

---

#### Matrix.__getitem__

```python
Matrix[i, j] -> str
```

Returns the element at position `(i, j)` (as a string).

**Exceptions**: `ValueError` — index out of bounds.

**Example**:

```python
>>> a = Matrix([[1, 2], [3, 4]])
>>> a[1, 0]
'3'
```

---

#### Matrix.rows

```python
Matrix.rows() -> list[list[str]]
```

Returns the two-dimensional list of all rows (as strings).

---

#### Matrix.transpose

```python
Matrix.transpose() -> Matrix
```

Returns the transposed matrix.

---

#### Matrix.trace

```python
Matrix.trace() -> str
```

Returns the trace of the square matrix (sum of diagonal elements).

**Exceptions**: `ValueError` — the matrix is not square.

---

#### Matrix.rank

```python
Matrix.rank() -> int
```

Returns the rank of the matrix.

---

#### Matrix.determinant

```python
Matrix.determinant() -> str
```

Returns the determinant of the square matrix (Bareiss fraction-free algorithm).

**Exceptions**: `ValueError` — the matrix is not square.

**Example**:

```python
>>> a = Matrix([[1, 2], [3, 4]])
>>> a.determinant()
'-2'
```

---

#### Matrix.inverse

```python
Matrix.inverse() -> Matrix
```

Returns the inverse matrix.

**Exceptions**: `ValueError` — the matrix is singular or not square.

---

#### Matrix.solve

```python
Matrix.solve(rhs) -> list[str]
```

Solves the linear system `self · x = rhs`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rhs` | `list` | The right-hand side vector |

**Returns**: the solution vector (as a list of strings).

**Exceptions**: `ValueError` — no solution or the solution is not unique.

**Example**:

```python
>>> a = Matrix([[2, 0], [0, 3]])
>>> a.solve([4, 9])
['2', '3']
```

---

#### Matrix arithmetic

| Operation | Syntax | Description |
|---|---|---|
| Matrix multiplication | `a @ b` | `__matmul__` |
| Addition | `a + b` | Element-wise addition |
| Subtraction | `a - b` | Element-wise subtraction |
| Equality | `a == b` | Element-wise comparison |

---

## Coefficient Domains

### IntegerDomain

Selector class for the integer domain ℤ.

**Signature**:

```python
ocas.IntegerDomain()
```

> **Note**: the `domain` parameter of `Polynomial` and `Matrix` currently accepts strings (`"integer"` / `"rational"`) or a `FiniteField` instance; `IntegerDomain()` / `RationalDomain()` instances cannot yet be passed directly as that parameter.

**Example**:

```python
>>> from ocas import IntegerDomain, Polynomial
>>> d = IntegerDomain()
>>> repr(d)
'IntegerDomain()'
>>> p = Polynomial([1, 2, 1], domain="integer")
>>> p.coeffs()
['1', '2', '1']
```

---

### RationalDomain

Selector for the rational number domain ℚ.

**Signature**:

```python
ocas.RationalDomain()
```

**Example**:

```python
>>> from ocas import RationalDomain, Polynomial
>>> d = RationalDomain()
>>> repr(d)
'RationalDomain()'
>>> p = Polynomial([(1, 2), (3, 4)], domain="rational")
>>> p.coeffs()
['1/2', '3/4']
```

---

### FiniteField

The finite field GF(p), where $p$ is a prime.

**Signature**:

```python
ocas.FiniteField(modulus: int)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `modulus` | `int` | Prime modulus, $p \geq 2$ |

**Returns**: a `FiniteField` instance.

**Exceptions**: `ValueError` — `modulus < 2`.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `modulus` | `str` | Decimal string of the prime modulus |

**Example**:

```python
>>> from ocas import FiniteField, Polynomial
>>> gf5 = FiniteField(5)
>>> p = Polynomial([1, 2, 1], domain=gf5)
>>> print(repr(gf5))
FiniteField(5)
```

---

## Solvers

### solve_linear_rational

```python
ocas.solve_linear_rational(a: list[list[int]], b: list[int]) -> list[tuple[int, int]]
```

Solves the linear system $A \mathbf{x} = \mathbf{b}$ over ℚ.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `list[list[int]]` | Coefficient matrix (row-major, integers) |
| `b` | `list[int]` | Right-hand side vector (integers) |

**Returns**: the solution vector, each component being a `(numerator, denominator)` tuple.

**Exceptions**: `ValueError` — the solution is not unique.

**Example**:

```python
>>> ocas.solve_linear_rational([[2, 0], [0, 3]], [4, 9])
[(2, 1), (3, 1)]
```

---

### solve_linear_integer

```python
ocas.solve_linear_integer(a: list[list[int]], b: list[int]) -> list[int]
```

Solves the linear system $A \mathbf{x} = \mathbf{b}$ over ℤ.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `list[list[int]]` | Coefficient matrix (row-major, integers) |
| `b` | `list[int]` | Right-hand side vector (integers) |

**Returns**: the integer solution vector.

**Exceptions**: `ValueError` — no integer solution exists.

**Example**:

```python
>>> ocas.solve_linear_integer([[2, 0], [0, 3]], [4, 9])
[2, 3]
```

---

### solve_diophantine

```python
ocas.solve_diophantine(a: int, b: int, c: int) -> DiophantineSolution | None
```

Solves the linear Diophantine equation $ax + by = c$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `int` | Coefficient of $x$ |
| `b` | `int` | Coefficient of $y$ |
| `c` | `int` | Right-hand side constant |

**Returns**: a `DiophantineSolution` (when a solution exists) or `None` (when no integer solution exists).

**Example**:

```python
>>> sol = ocas.solve_diophantine(3, 5, 1)
>>> sol.particular
(2, -1)
>>> sol.general
(5, -3)
```

---

### DiophantineSolution

The solution of a Diophantine equation.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `particular` | `tuple[int, int]` | A particular solution $(x_0, y_0)$ |
| `general` | `tuple[int, int]` | The general-solution direction $(t_x, t_y)$; the general solution is $(x_0 + k \cdot t_x, \; y_0 + k \cdot t_y)$, $k \in \mathbb{Z}$ |

---

## ODE

### classify_ode

```python
ocas.classify_ode(equation: Expression, func: str, var: str) -> list[str]
```

Classifies an ordinary differential equation and returns the list of available solver method names.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equation` | `Expression` | The ODE expression (in the form equal to zero), e.g. `Derivative(y(x), x) - y(x)` |
| `func` | `str` | Name of the unknown function, e.g. `"y"` |
| `var` | `str` | Name of the independent variable, e.g. `"x"` |

**Returns**: a list of method names, e.g. `["LinearFirst", "PowerSeries"]`.

**Supported ODE types**: `Separable`, `LinearFirst`, `Bernoulli`, `Exact`, `Homogeneous`, `LinearConstantCoeff`, `CauchyEuler`, `ReductionOfOrder`, `PowerSeries`.

**Example**:

```python
>>> eq = ocas.Expression("Derivative(y(x), x) - y(x)")
>>> ocas.classify_ode(eq, "y", "x")
['LinearFirst', 'PowerSeries']
```

---

### dsolve

```python
ocas.dsolve(equation: Expression, func: str, var: str, hint: str | None = None) -> str
```

Symbolically solves an ODE.

**Signature** (with defaults): `dsolve(equation, func, var, hint=None)`

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equation` | `Expression` | The ODE expression (in the form equal to zero) |
| `func` | `str` | Name of the unknown function |
| `var` | `str` | Name of the independent variable |
| `hint` | `str` or `None` | A specific solver method (one of the names returned by `classify_ode`) |

**Returns**: the string representation of the solution (explicit, implicit, series, or unevaluated form).

**Example**:

```python
>>> eq = ocas.Expression("Derivative(y(x), x) - y(x)")
>>> ocas.dsolve(eq, "y", "x")
'y = C1*exp(x)'
```

---

### dsolve_ivp

```python
ocas.dsolve_ivp(equation: Expression, func: str, var: str, y0: str, y1: str | None = None) -> str
```

Solves an initial value problem for a constant-coefficient linear ODE via the Laplace transform.

**Signature** (with defaults): `dsolve_ivp(equation, func, var, y0, y1=None)`

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equation` | `Expression` | The ODE expression (in the form equal to zero) |
| `func` | `str` | Name of the unknown function |
| `var` | `str` | Name of the independent variable |
| `y0` | `str` | The value of $y(0)$ (string expression, e.g. `"1"`) |
| `y1` | `str` or `None` | The value of $y'(0)$ (required for second-order problems) |

**Returns**: the explicit solution string without arbitrary constants.

**Example**:

```python
>>> eq = ocas.Expression("Derivative(y(x), x, 2) + y(x)")
>>> ocas.dsolve_ivp(eq, "y", "x", "0", "1")
'y = sin(x)'
```

---

## Number Theory

### factorint

```python
ocas.factorint(n: int) -> list[tuple[str, int]]
```

Factors $|n|$ into prime factors. Negative inputs start with `("-1", 1)`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | The integer to be factored (arbitrary precision) |

**Returns**: a list of `(prime, exponent)` tuples, sorted by prime in ascending order.

**Example**:

```python
>>> ocas.factorint(60)
[('2', 2), ('3', 1), ('5', 1)]
>>> ocas.factorint(-12)
[('-1', 1), ('2', 2), ('3', 1)]
```

---

### isprime

```python
ocas.isprime(n: int) -> bool
```

BPSW probabilistic primality test. Deterministic for $n < 2^{64}$; no composite numbers of any size are currently known to pass.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | The integer to test (arbitrary precision) |

**Returns**: `True` if (probably) prime.

**Example**:

```python
>>> ocas.isprime(97)
True
>>> ocas.isprime(100)
False
```

---

### isprime_u64

```python
ocas.isprime_u64(n: int) -> bool
```

Deterministic primality test for $n < 2^{64}$ (u64 fast path).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | The integer to test (must fit in the u64 range) |

**Returns**: `True` if prime.

**Example**:

```python
>>> ocas.isprime_u64(2**61 - 1)
True
```

---

### nextprime

```python
ocas.nextprime(n: int) -> int
```

Returns the smallest prime strictly greater than $n$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | The starting integer |

**Returns**: the next prime (arbitrary precision).

**Example**:

```python
>>> ocas.nextprime(10)
11
>>> ocas.nextprime(11)
13
```

---

### discrete_log

```python
ocas.discrete_log(p: int, base: int, target: int) -> int
```

Solves the discrete logarithm $\text{base}^x \equiv \text{target} \pmod{p}$. Uses the Pohlig–Hellman algorithm for prime $p$, otherwise BSGS.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `p` | `int` | The modulus |
| `base` | `int` | The base |
| `target` | `int` | The target value |

**Returns**: the discrete logarithm $x$.

**Exceptions**: `ValueError` — the logarithm does not exist.

**Example**:

```python
>>> ocas.discrete_log(13, 2, 9)
8
>>> pow(2, 8, 13)
9
```

---

### crt

```python
ocas.crt(moduli: list[int], residues: list[int]) -> tuple[int, int]
```

Chinese remainder theorem: solves the system of congruences $x \equiv r_i \pmod{m_i}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `moduli` | `list[int]` | List of moduli (need not be pairwise coprime) |
| `residues` | `list[int]` | List of residues |

**Returns**: `(r, m)` such that $x \equiv r \pmod{m}$.

**Exceptions**: `ValueError` — the system is inconsistent.

**Example**:

```python
>>> ocas.crt([3, 5, 7], [2, 3, 2])
(23, 105)
```

---

### jacobi_symbol

```python
ocas.jacobi_symbol(a: int, n: int) -> int
```

Computes the Jacobi symbol $\left(\frac{a}{n}\right)$, where $n$ is a positive odd integer.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `int` | Numerator |
| `n` | `int` | Denominator (positive odd) |

**Returns**: $-1$, $0$, or $1$.

**Exceptions**: `ValueError` — `n` is not a positive odd integer.

**Example**:

```python
>>> ocas.jacobi_symbol(2, 7)
1
>>> ocas.jacobi_symbol(3, 7)
-1
```

---

### totient

```python
ocas.totient(n: int) -> int
```

Euler totient function $\varphi(n)$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | A positive integer |

**Returns**: $\varphi(n)$ (arbitrary precision).

**Example**:

```python
>>> ocas.totient(12)
4
```

---

### mobius

```python
ocas.mobius(n: int) -> int
```

Möbius function $\mu(n)$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | A positive integer |

**Returns**: $-1$, $0$, or $1$.

**Example**:

```python
>>> ocas.mobius(6)
1
>>> ocas.mobius(4)
0
```

---

### divisor_count

```python
ocas.divisor_count(n: int) -> int
```

Number of positive divisors $\tau(n)$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | A positive integer |

**Returns**: the number of positive divisors.

**Example**:

```python
>>> ocas.divisor_count(12)
6
```

---

### divisor_sigma

```python
ocas.divisor_sigma(n: int, k: int = 1) -> int
```

Sum of the $k$-th powers of the divisors $\sigma_k(n) = \sum_{d \mid n} d^k$.

**Signature** (with defaults): `divisor_sigma(n, k=1)`

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | A positive integer |
| `k` | `int` | The exponent (default 1, i.e. the sum of divisors $\sigma_1$) |

**Returns**: $\sigma_k(n)$.

**Example**:

```python
>>> ocas.divisor_sigma(12)      # 1+2+3+4+6+12
28
>>> ocas.divisor_sigma(12, 2)   # 1+4+9+16+36+144
210
```

---

### liouville_lambda

```python
ocas.liouville_lambda(n: int) -> int
```

Liouville function $\lambda(n) = (-1)^{\Omega(n)}$, where $\Omega(n)$ is the number of prime factors of $n$ counted with multiplicity.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | A positive integer |

**Returns**: $-1$ or $1$.

**Example**:

```python
>>> ocas.liouville_lambda(12)  # Ω(12)=3
-1
>>> ocas.liouville_lambda(9)   # Ω(9)=2
1
```

---

## Gröbner Bases & Ideals

> **Note**: The following functions keep their `py_` prefix (the Rust side does not use `#[pyo3(name=...)]` renaming). The actual call names in Python are `ocas.py_groebner_basis`, `ocas.py_ideal_contains`, etc.

### py_groebner_basis

```python
ocas.py_groebner_basis(generators, n_vars: int = 1, algorithm: str = "auto") -> GroebnerBasis
```

Computes a Gröbner basis of the polynomial ideal.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `list` | List of generators (`MultivariatePolynomial` or `Polynomial` objects) |
| `n_vars` | `int` | Number of variables (default 1; auto-detected if the generators are `MultivariatePolynomial`) |
| `algorithm` | `str` | Algorithm selection: `"auto"` (default), `"f4"`, `"f5"`, `"buchberger"` |

**Returns**: a `GroebnerBasis` object.

**Exceptions**: `ValueError` — invalid algorithm name.

**Example**:

```python
>>> from ocas import MultivariatePolynomial, py_groebner_basis
>>> f = MultivariatePolynomial({(2, 0): 1, (0, 0): -1}, n_vars=2)   # x² - 1
>>> g = MultivariatePolynomial({(0, 2): 1, (0, 0): -1}, n_vars=2)   # y² - 1
>>> gb = py_groebner_basis([f, g], n_vars=2)
>>> len(gb)
2
```

---

### py_ideal_contains

```python
ocas.py_ideal_contains(generators, f, n_vars: int = 1, algorithm: str = "auto") -> bool
```

Determines whether the polynomial $f$ belongs to the ideal generated by `generators`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `list` | List of generators of the ideal |
| `f` | `MultivariatePolynomial` or `Polynomial` | The polynomial to test |
| `n_vars` | `int` | Number of variables (default 1) |
| `algorithm` | `str` | Algorithm selection (default `"auto"`) |

**Returns**: `True` if $f \in I$.

**Example**:

```python
>>> from ocas import MultivariatePolynomial, py_ideal_contains
>>> f = MultivariatePolynomial({(2, 0): 1, (0, 0): -1}, n_vars=2)
>>> g = MultivariatePolynomial({(0, 2): 1, (0, 0): -1}, n_vars=2)
>>> x = MultivariatePolynomial({(1, 0): 1}, n_vars=2)
>>> py_ideal_contains([f, g], f, n_vars=2)
True
>>> py_ideal_contains([f, g], x, n_vars=2)
False
```

---

### py_solve_polynomial_system

```python
ocas.py_solve_polynomial_system(equations, n_vars: int = 1, algorithm: str = "auto") -> PolynomialSystemSolution
```

Solves a system of polynomial equations.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equations` | `list` | List of equations (each a polynomial in the form equal to zero) |
| `n_vars` | `int` | Number of variables (default 1) |
| `algorithm` | `str` | Algorithm selection (default `"auto"`) |

**Returns**: a `PolynomialSystemSolution` object.

**Example**:

```python
>>> from ocas import MultivariatePolynomial, py_solve_polynomial_system
>>> f = MultivariatePolynomial({(2, 0): 1, (0, 2): 1, (0, 0): -1}, n_vars=2)
>>> g = MultivariatePolynomial({(1, 0): 1, (0, 0): -1}, n_vars=2)
>>> sol = py_solve_polynomial_system([f, g], n_vars=2)
>>> sol.kind
'zero_dimensional'
>>> for s in sol.solutions():
...     print(s.values, s.multiplicity)
[1.0, 0.0] 1
```

---

### py_hilbert_series

```python
ocas.py_hilbert_series(gb: GroebnerBasis) -> HilbertSeries
```

Computes the Hilbert series of the Gröbner basis.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `gb` | `GroebnerBasis` | A Gröbner basis |

**Returns**: a `HilbertSeries` object.

**Example**:

```python
>>> from ocas import MultivariatePolynomial, py_groebner_basis, py_hilbert_series
>>> f = MultivariatePolynomial({(2, 0): 1, (0, 0): -1}, n_vars=2)
>>> gb = py_groebner_basis([f], n_vars=2)
>>> hs = py_hilbert_series(gb)
>>> hs.dimension
1
```

---

### py_ideal_radical

```python
ocas.py_ideal_radical(generators, n_vars: int = 1) -> GroebnerBasis
```

Computes the radical $\sqrt{I}$ of the ideal.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `list` | List of generators of the ideal |
| `n_vars` | `int` | Number of variables (default 1) |

**Returns**: a Gröbner basis of $\sqrt{I}$.

**Algorithm**: square-free factorization for zero-dimensional ideals; Jacobian saturation $\sqrt{I} = I : h^\infty$ for positive-dimensional ideals.

---

### py_primary_decomposition

```python
ocas.py_primary_decomposition(generators, n_vars: int = 1) -> list[PrimaryComponent]
```

Computes the primary decomposition of the ideal.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `list` | List of generators of the ideal |
| `n_vars` | `int` | Number of variables (default 1) |

**Returns**: a list of `PrimaryComponent`.

**Example**:

```python
>>> from ocas import MultivariatePolynomial, py_primary_decomposition
>>> f = MultivariatePolynomial({(2, 0): 1, (0, 0): -1}, n_vars=1)
>>> comps = py_primary_decomposition([f], n_vars=1)
>>> len(comps)
2
```

---

### py_is_zero_dimensional

```python
ocas.py_is_zero_dimensional(gb: GroebnerBasis) -> bool
```

Determines whether the ideal of the Gröbner basis is zero-dimensional.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `gb` | `GroebnerBasis` | A Gröbner basis |

**Returns**: `True` if zero-dimensional.

---

### py_eliminate

```python
ocas.py_eliminate(generators, elim_vars: int, n_vars: int = 1, algorithm: str = "auto") -> GroebnerBasis
```

Elimination: computes a Gröbner basis of the elimination ideal $I \cap k[x_{e+1}, \dots, x_n]$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `list` | List of generators of the ideal |
| `elim_vars` | `int` | Number of variables to eliminate (from the front) |
| `n_vars` | `int` | Total number of variables (default 1) |
| `algorithm` | `str` | Algorithm selection (default `"auto"`) |

**Returns**: a Gröbner basis of the elimination ideal.

---

### MultivariatePolynomial (Gröbner helper class)

Already introduced in the Polynomials section. The Gröbner functions accept both `MultivariatePolynomial` and `Polynomial` as generators; the latter are automatically mapped to univariate multivariate polynomials.

---

### GroebnerBasis

The result of a Gröbner basis computation.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `n_vars` | `int` | Number of variables |

**Methods**:

| Method | Returns | Description |
|---|---|---|
| `__len__()` | `int` | Number of polynomials in the basis |
| `is_groebner_basis()` | `bool` | Verifies whether it is indeed a Gröbner basis |

---

### RealSolution

A real solution of a polynomial system.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `values` | `list[float]` | Coordinates of the solution |
| `multiplicity` | `int` | The multiplicity |

---

### PolynomialSystemSolution

The solving result of a polynomial system.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `kind` | `str` | `"zero_dimensional"`, `"positive_dimensional"`, or `"empty"` |
| `vector_space_dimension` | `int` or `None` | Dimension of the vector space (only set for the zero-dimensional case) |

**Methods**:

| Method | Returns | Description |
|---|---|---|
| `solutions()` | `list[RealSolution]` | List of real solutions (non-empty only in the zero-dimensional case) |

---

### HilbertSeries

The Hilbert series.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `dimension` | `int` | Krull dimension |
| `degree` | `int` | Multiplicity (degree) |
| `numerator` | `list[int]` | Coefficients of the numerator polynomial |

**Methods**:

| Method | Parameters | Returns | Description |
|---|---|---|---|
| `hilbert_function(d)` | `d: int` | `int` | Value of the Hilbert function at degree $d$ |

---

### PrimaryComponent

A component of a primary decomposition.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `n_vars` | `int` | Number of variables |

---

## Algebraic Numbers

### AlgebraicExtension

The algebraic number field $K = \mathbb{Q}(\alpha)$, defined by a monic minimal polynomial.

**Signature**:

```python
ocas.AlgebraicExtension(min_poly: list)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `min_poly` | `list[int]` or `list[tuple[int,int]]` | Coefficients of the minimal polynomial (ascending powers); the last entry (leading coefficient) must be 1 |

**Exceptions**: `ValueError` — the coefficient list is too short, the leading coefficient is not 1, or irreducibility verification fails.

**Example**:

```python
>>> from ocas import AlgebraicExtension
>>> # α² - 2 (i.e. ℚ(√2))
>>> field = AlgebraicExtension([-2, 0, 1])
>>> field.extension_degree()
2
```

---

#### AlgebraicExtension.extension_degree

```python
AlgebraicExtension.extension_degree() -> int
```

Returns the degree of the extension $[K:\mathbb{Q}] = \deg(m)$.

---

#### AlgebraicExtension.alpha

```python
AlgebraicExtension.alpha() -> AlgebraicElement
```

Returns the generator $\alpha$ of the extension.

**Example**:

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> a = field.alpha()
>>> a.coeffs()
['0', '1']
```

---

#### AlgebraicExtension.from_base

```python
AlgebraicExtension.from_base(c) -> AlgebraicElement
```

Embeds the rational number $c$ as a constant element of the field.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `c` | `int` or `tuple[int, int]` | A rational number (an integer or a `(numerator, denominator)` tuple) |

---

#### AlgebraicExtension.element

```python
AlgebraicExtension.element(coeffs: list) -> AlgebraicElement
```

Creates an element from the $\alpha$-polynomial coefficients (ascending powers).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `coeffs` | `list[int]` or `list[tuple[int,int]]` | List of rational coefficients (ascending powers) |

**Example**:

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> # 3 + 2α
>>> e = field.element([3, 2])
>>> e.coeffs()
['3', '2']
```

---

### AlgebraicElement

An element of the algebraic number field $\mathbb{Q}(\alpha)$, stored as a polynomial in $\alpha$.

#### AlgebraicElement.coeffs

```python
AlgebraicElement.coeffs() -> list[str]
```

Returns the $\alpha$-polynomial coefficients (list of strings, ascending powers).

**Example**:

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> e = field.element([1, 1])  # 1 + α = 1 + √2
>>> e.coeffs()
['1', '1']
```

---

### AlgebraicPolynomial

Dense univariate polynomial over the algebraic number field $\mathbb{Q}(\alpha)$.

**Signature**:

```python
ocas.AlgebraicPolynomial(field: AlgebraicExtension, coeffs: list)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `field` | `AlgebraicExtension` | The algebraic number field containing the coefficients |
| `coeffs` | `list` | Coefficient list (ascending powers), each entry being an `int`, an `(int, int)` tuple, or an `AlgebraicElement` |

**Example**:

```python
>>> from ocas import AlgebraicExtension, AlgebraicPolynomial
>>> field = AlgebraicExtension([-2, 0, 1])
>>> # x² - 2 over ℚ(√2)
>>> p = AlgebraicPolynomial(field, [-2, 0, 1])
>>> p.degree()
2
```

---

#### AlgebraicPolynomial.degree

```python
AlgebraicPolynomial.degree() -> int | None
```

Returns the degree. The zero polynomial returns `None`.

---

#### AlgebraicPolynomial.len

```python
AlgebraicPolynomial.len() -> int
```

Returns the number of stored coefficients.

---

#### AlgebraicPolynomial.is_zero

```python
AlgebraicPolynomial.is_zero() -> bool
```

Whether the polynomial is the zero polynomial.

---

#### AlgebraicPolynomial.coeffs

```python
AlgebraicPolynomial.coeffs() -> list[list[str]]
```

Returns the coefficient list. Each coefficient is a list of rationals (as strings) representing the $\alpha$-polynomial coefficients.

---

#### AlgebraicPolynomial.factor

```python
AlgebraicPolynomial.factor() -> list[AlgebraicFactor]
```

Factors the polynomial over the algebraic number field (Trager's algorithm).

**Returns**: a list of `AlgebraicFactor`, each containing an irreducible factor and its multiplicity.

**Example**:

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> p = AlgebraicPolynomial(field, [-2, 0, 1])  # x² - 2 = (x-α)(x+α)
>>> facs = p.factor()
>>> len(facs)
2
>>> [f.multiplicity for f in facs]
[1, 1]
```

---

### AlgebraicFactor

A factor in the factorization of an algebraic-number-field polynomial.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `factor` | `AlgebraicPolynomial` | The irreducible factor |
| `multiplicity` | `int` | The multiplicity |

---

## Automatic Differentiation

### DualShape

First-order dual-number layout description, declaring the number of tracked differentiation variables. Construct once, share everywhere.

**Signature**:

```python
# Not constructed directly; use the static method
shape = DualShape.first_order(n_vars)
```

---

#### DualShape.first_order (static method)

```python
DualShape.first_order(n_vars: int) -> DualShape
```

Creates a layout with `n_vars` differentiation variables.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n_vars` | `int` | Number of differentiation variables ($\geq 1$) |

**Exceptions**: `ValueError` — `n_vars == 0`.

**Example**:

```python
>>> from ocas import DualShape
>>> shape = DualShape.first_order(3)
>>> shape.n_vars
3
>>> shape.n_components
4
```

---

#### DualShape.n_vars

```python
DualShape.n_vars -> int  # property
```

Number of differentiation variables.

---

#### DualShape.n_components

```python
DualShape.n_components -> int  # property
```

Total number of components (value + derivative slots).

---

### HyperDual

Hyper-dual number: a scalar value together with partial derivatives with respect to the variables of a `DualShape`. Computes exactly with rational numbers.

**Limitations**: supports only polynomial/rational arithmetic (`+`, `-`, `*`, `/`, negation). Transcendental functions (sin/exp/log) are not implemented; use repeated multiplication for integer powers.

---

#### HyperDual.variable (static method)

```python
HyperDual.variable(shape: DualShape, i: int, value) -> HyperDual
```

Creates an independent variable $x_i = \text{value}$ (derivative 1 with respect to the $i$-th variable, 0 with respect to the others).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `shape` | `DualShape` | The dual-number layout |
| `i` | `int` | Variable index |
| `value` | `int` or `tuple[int, int]` | The variable value (an integer or a rational tuple) |

**Exceptions**: `ValueError` — index out of bounds.

---

#### HyperDual.constant (static method)

```python
HyperDual.constant(shape: DualShape, value) -> HyperDual
```

Creates a constant (all derivatives zero).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `shape` | `DualShape` | The dual-number layout |
| `value` | `int` or `tuple[int, int]` | The constant value |

---

#### HyperDual.value

```python
HyperDual.value() -> str
```

Returns the scalar value (as a string, e.g. `"5"` or `"3/7"`).

---

#### HyperDual.deriv

```python
HyperDual.deriv(i: int) -> str | None
```

Returns the partial derivative with respect to the $i$-th variable. Returns `None` if `i` is out of bounds.

---

#### HyperDual.n_vars

```python
HyperDual.n_vars -> int  # property
```

Number of differentiation variables.

---

#### HyperDual arithmetic

| Operation | Syntax | Description |
|---|---|---|
| Addition | `a + b` | Component-wise addition |
| Subtraction | `a - b` | Component-wise subtraction |
| Multiplication | `a * b` | Product rule |
| Division | `a / b` | Quotient rule |
| Negation | `-a` | Component-wise negation |

All operations require both operands to share the same `DualShape`. Different shapes raise `ValueError`.

**Complete example**:

```python
>>> from ocas import DualShape, HyperDual
>>> shape = DualShape.first_order(2)
>>> x = HyperDual.variable(shape, 0, 3)
>>> y = HyperDual.variable(shape, 1, 5)
>>> f = x * y
>>> f.value()
'15'
>>> f.deriv(0)   # ∂f/∂x = y
'5'
>>> f.deriv(1)   # ∂f/∂y = x
'3'
```

---

## Tensors

### Tensor

Named tensor object with index slots and optional symmetry.

**Signature**:

```python
ocas.Tensor(name: str, slots: list[tuple[str, str]], symmetry: str = "none")
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `str` | Tensor name |
| `slots` | `list[tuple[str, str]]` | List of index slots, each being `(label, position)`. Position is `"upper"` (contravariant) or `"lower"` (covariant) |
| `symmetry` | `str` | Symmetry: `"none"` (default), `"symmetric"`, or `"antisymmetric"` |

**Exceptions**: `ValueError` — invalid position or symmetry string.

**Example**:

```python
>>> from ocas import Tensor
>>> t = Tensor("T", [("i", "upper"), ("j", "lower")])
>>> t.name
'T'
>>> t.rank
2
```

---

#### Tensor.name

```python
Tensor.name -> str  # property
```

The tensor name.

---

#### Tensor.rank

```python
Tensor.rank -> int  # property
```

The tensor rank (number of slots).

---

#### Tensor.symmetry

```python
Tensor.symmetry -> str  # property
```

The symmetry string (`"none"`, `"symmetric"`, `"antisymmetric"`).

---

#### Tensor.slots

```python
Tensor.slots() -> list[tuple[str, str]]
```

Returns the list of slots, each being `(label, position)`.

---

#### Tensor.dummy_labels

```python
Tensor.dummy_labels() -> list[str]
```

Returns the labels of the current dummy indices (those appearing twice).

---

#### Tensor.to_string_atom

```python
Tensor.to_string_atom() -> str
```

Renders the tensor as the string form of an expression function node `name(slot, slot, ...)`.

---

### contract_tensors

```python
ocas.contract_tensors(a: Tensor, b: Tensor) -> tuple[str, list | str]
```

Contracts two tensors: sums over shared dummy indices (same label, opposite variance).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `Tensor` | The first tensor |
| `b` | `Tensor` | The second tensor |

**Returns**: a `(kind, payload)` tuple:
- when `kind = "product"`, `payload` is the list of tensors after contraction (free indices concatenated)
- when `kind = "scalar"`, `payload` is the string of the contracted scalar expression

**Example**:

```python
>>> from ocas import Tensor, contract_tensors
>>> t = Tensor("T", [("i", "upper"), ("j", "lower")])
>>> u = Tensor("U", [("j", "upper"), ("k", "lower")])
>>> kind, payload = contract_tensors(t, u)
>>> kind
'product'
```

---

### tensor_symmetrise_sign

```python
ocas.tensor_symmetrise_sign(tensor: Tensor) -> int
```

Returns the symmetrization sign of the tensor ($+1$ or $-1$). `"none"` and `"symmetric"` always return $+1$; `"antisymmetric"` returns the parity of the slot permutation.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `tensor` | `Tensor` | The input tensor |

**Returns**: $+1$ or $-1$.

**Example**:

```python
>>> from ocas import Tensor, tensor_symmetrise_sign
>>> eps = Tensor("eps", [("a", "lower"), ("b", "lower")], symmetry="antisymmetric")
>>> tensor_symmetrise_sign(eps) in (1, -1)
True
```

---

### canonicalize_tensors

```python
ocas.canonicalize_tensors(expr: str, specs: dict[str, str], index_groups: dict[str, int] | None = None) -> str
```

Canonicalizes a tensor expression using the graph isomorphism engine.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `expr` | `str` | The tensor expression string |
| `specs` | `dict[str, str]` | Tensor name → symmetry specification (`"none"`, `"symmetric"`, `"antisymmetric"`) |
| `index_groups` | `dict[str, int]` or `None` | Index dimension groups (optional) |

**Returns**: the canonicalized expression string.

**Exceptions**: `ValueError` — parsing or canonicalization failed.

---

### young_project

```python
ocas.young_project(expr: str, tableau: list[int]) -> str
```

Applies a Young projection to a tensor expression.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `expr` | `str` | The tensor expression string |
| `tableau` | `list[int]` | Row lengths of the Young diagram, e.g. `[2, 1]` denotes □□/□ |

**Returns**: the projected expression string.

**Example**:

```python
>>> proj = ocas.young_project("f(a,b,c)", [1, 1, 1])  # □/□/□ fully antisymmetric projection
>>> "f(a,b,c)" in proj
True
```

---

### refresh_dummies

```python
ocas.refresh_dummies(expr: str, specs: dict[str, str]) -> str
```

Renames (refreshes) the dummy indices in a tensor expression.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `expr` | `str` | The tensor expression string |
| `specs` | `dict[str, str]` | Tensor name → symmetry specification |

**Returns**: the expression string with the dummy indices renamed.

---

## Evaluation

### ExpressionEvaluator

Compiles an expression for fast numerical evaluation. Compile once, evaluate many times.

**Signature**:

```python
ocas.ExpressionEvaluator(input: str, param_names: list[str])
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `input` | `str` | The expression string |
| `param_names` | `list[str]` | List of parameter names (in evaluation order) |

**Exceptions**: `ValueError` — parsing or compilation error.

**Example**:

```python
>>> from ocas import ExpressionEvaluator
>>> ev = ExpressionEvaluator("x^2 + y", ["x", "y"])
>>> ev.evaluate([3.0, 1.0])
[10.0]
>>> ev.evaluate([2.0, 0.0])
[4.0]
```

---

#### ExpressionEvaluator.evaluate

```python
ExpressionEvaluator.evaluate(values: list[float]) -> list[float]
```

Evaluates with the given parameter values.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `values` | `list[float]` | List of floating-point values; its length must match the number of parameters |

**Returns**: the list of result floats.

**Exceptions**: `ValueError` — parameter count mismatch or evaluation error.

---

#### ExpressionEvaluator.n_params

```python
ExpressionEvaluator.n_params -> int  # property
```

The number of parameters.

---

## Numerical Integration

### Vegas

Adaptive Monte Carlo integrator (the Vegas algorithm), supporting multi-dimensional integrals.

**Signature**:

```python
ocas.Vegas(n_dims: int, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n_dims` | `int` | The integration dimension |
| `n_bins` | `int` or `None` | Number of bins per dimension (default 64) |
| `n_samples` | `int` or `None` | Number of samples per iteration (default 10000) |
| `iterations` | `int` or `None` | Number of iterations (default 10) |
| `learning_rate` | `float` or `None` | Grid adaptation learning rate (default 1.5) |
| `seed` | `int` or `None` | RNG seed (default `0x0C45`) |

**Example**:

```python
>>> from ocas import Vegas
>>> v = Vegas(2, n_samples=20000, iterations=8, seed=1)
>>> r = v.integrate(lambda xs: xs[0] * xs[1])
>>> abs(r.integral - 0.25) < 0.01
True
```

---

#### Vegas.integrate

```python
Vegas.integrate(f: Callable[[list[float]], float]) -> IntegrateResult
```

Integrates the callable $f$ over the unit hypercube $[0, 1]^n$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `Callable` | Accepts a list of `n_dims` floats in `[0, 1]` and returns a single float |

**Returns**: an `IntegrateResult`.

---

#### Vegas.result

```python
Vegas.result -> IntegrateResult  # property
```

The latest accumulated estimate and error.

---

#### Vegas.iterations

```python
Vegas.iterations -> int  # property
```

The number of completed iterations.

---

### integrate_1d

```python
ocas.integrate_1d(f: Callable[[float], float], a: float, b: float, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None) -> IntegrateResult
```

Convenience function for one-dimensional numerical integration.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `Callable` | A unary floating-point function |
| `a` | `float` | Lower integration limit |
| `b` | `float` | Upper integration limit |
| `n_bins` etc. | Same as `Vegas` | Optional tuning parameters |

**Returns**: an `IntegrateResult`.

**Example**:

```python
>>> import ocas
>>> r = ocas.integrate_1d(lambda x: x, 0.0, 1.0)
>>> abs(r.integral - 0.5) < 0.01
True
```

---

### IntegrateResult

The result of a numerical integration.

**Attributes**:

| Attribute | Type | Description |
|---|---|---|
| `integral` | `float` | The integral estimate |
| `error` | `float` | The estimated standard error |

Supports indexing: `result[0]` = `integral`, `result[1]` = `error`; supports unpacking: `integral, error = result`.

**Example**:

```python
>>> import ocas
>>> r = ocas.integrate_1d(lambda x: x**2, 0.0, 1.0, seed=42)
>>> r[0]       # integral
0.333...
>>> r[1]       # error
0.00...
>>> i, e = r   # unpacking
```

---

## Double-Precision Floats

### DoubleF64

Double-precision floating-point arithmetic (~31 significant decimal digits, ~84 binary bits), using the Dekker/Knuth "double-float" algorithm.

**Signature**:

```python
ocas.DoubleF64(hi: float, lo: float = 0.0)
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `hi` | `float` | The primary value |
| `lo` | `float` | The error term (default 0.0) |

**Example**:

```python
>>> from ocas import DoubleF64
>>> a = DoubleF64(1.0)
>>> b = DoubleF64(2.0)
>>> str(a + b)
'3'
```

---

#### DoubleF64.to_f64

```python
DoubleF64.to_f64() -> float
```

Converts to a standard `float` (loses the precision of the error term).

---

#### DoubleF64.components

```python
DoubleF64.components() -> tuple[float, float]
```

Returns the `(hi, lo)` tuple.

---

#### DoubleF64 arithmetic

| Operation | Syntax | Description |
|---|---|---|
| Addition | `a + b` | Double-precision addition |
| Subtraction | `a - b` | Double-precision subtraction |
| Multiplication | `a * b` | Double-precision multiplication |
| Division | `a / b` | Double-precision division |
| Negation | `-a` | Component-wise negation |
| Absolute value | `abs(a)` | Double-precision absolute value |
| Power | `a ** n` | Integer power (`n` is an `int`) |
| Comparison | `==`, `<`, `<=`, `>`, `>=` | Total order comparison |

**Exceptions**: `ValueError` — division by zero (`__truediv__`).

---

#### DoubleF64 transcendental functions

| Method | Description | Exceptions |
|---|---|---|
| `sin()` | Sine | — |
| `cos()` | Cosine | — |
| `tan()` | Tangent | — |
| `exp()` | Natural exponential $e^x$ | — |
| `ln()` | Natural logarithm | `ValueError`: non-positive input |
| `sqrt()` | Square root | `ValueError`: negative input |

**Complete example**:

```python
>>> from ocas import DoubleF64
>>> DoubleF64(0.0).sin().to_f64()
0.0
>>> DoubleF64(4.0).sqrt().components()
(2.0, 0.0)
```

---

## Function Quick Reference

| Function | Signature | Description |
|---|---|---|
| `solve_linear_rational` | `(a, b) -> list[tuple]` | Solves Ax=b over ℚ |
| `solve_linear_integer` | `(a, b) -> list[int]` | Solves Ax=b over ℤ |
| `solve_diophantine` | `(a, b, c) -> DiophantineSolution \| None` | ax+by=c |
| `classify_ode` | `(equation, func, var) -> list[str]` | ODE classification |
| `dsolve` | `(equation, func, var, hint=None) -> str` | Symbolically solves ODEs |
| `dsolve_ivp` | `(equation, func, var, y0, y1=None) -> str` | Laplace IVP |
| `factorint` | `(n) -> list[tuple]` | Integer prime factorization |
| `isprime` | `(n) -> bool` | BPSW primality test |
| `isprime_u64` | `(n) -> bool` | Deterministic primality (u64) |
| `nextprime` | `(n) -> int` | Next prime |
| `discrete_log` | `(p, base, target) -> int` | Discrete logarithm |
| `crt` | `(moduli, residues) -> tuple` | Chinese remainder theorem |
| `jacobi_symbol` | `(a, n) -> int` | Jacobi symbol |
| `totient` | `(n) -> int` | Euler totient |
| `mobius` | `(n) -> int` | Möbius function |
| `divisor_count` | `(n) -> int` | Number of divisors |
| `divisor_sigma` | `(n, k=1) -> int` | Sum of divisor powers |
| `liouville_lambda` | `(n) -> int` | Liouville function |
| `py_groebner_basis` | `(generators, n_vars=1, algorithm="auto")` | Gröbner basis |
| `py_ideal_contains` | `(generators, f, n_vars=1, algorithm="auto")` | Ideal membership test |
| `py_solve_polynomial_system` | `(equations, n_vars=1, algorithm="auto")` | Solves polynomial systems |
| `py_hilbert_series` | `(gb) -> HilbertSeries` | Hilbert series |
| `py_ideal_radical` | `(generators, n_vars=1)` | Ideal radical |
| `py_primary_decomposition` | `(generators, n_vars=1)` | Primary decomposition |
| `py_is_zero_dimensional` | `(gb) -> bool` | Zero-dimensionality test |
| `py_eliminate` | `(generators, elim_vars, n_vars=1, algorithm="auto")` | Elimination |
| `contract_tensors` | `(a, b) -> tuple` | Tensor contraction |
| `tensor_symmetrise_sign` | `(tensor) -> int` | Symmetrization sign |
| `canonicalize_tensors` | `(expr, specs, index_groups=None)` | Tensor canonicalization |
| `young_project` | `(expr, tableau) -> str` | Young projection |
| `refresh_dummies` | `(expr, specs) -> str` | Dummy index refresh |
| `integrate_1d` | `(f, a, b, **opts) -> IntegrateResult` | One-dimensional numerical integration |

---

## Class Quick Reference

| Class | Constructor | Description |
|---|---|---|
| `Expression` | `Expression(input)` | Symbolic expression |
| `Polynomial` | `Polynomial(coeffs, domain=None)` | Univariate polynomial |
| `MultivariatePolynomial` | `MultivariatePolynomial(terms, n_vars)` | Multivariate polynomial |
| `Matrix` | `Matrix(rows, domain=None)` | Matrix over a domain |
| `IntegerDomain` | `IntegerDomain()` | ℤ selector |
| `RationalDomain` | `RationalDomain()` | ℚ selector |
| `FiniteField` | `FiniteField(modulus)` | GF(p) |
| `AlgebraicExtension` | `AlgebraicExtension(min_poly)` | Algebraic number field |
| `AlgebraicElement` | (created by `AlgebraicExtension`) | Field element |
| `AlgebraicPolynomial` | `AlgebraicPolynomial(field, coeffs)` | Polynomial over an algebraic field |
| `DualShape` | `DualShape.first_order(n_vars)` | Dual-number layout |
| `HyperDual` | `HyperDual.variable(...)` / `.constant(...)` | Hyper-dual number |
| `Tensor` | `Tensor(name, slots, symmetry="none")` | Named tensor |
| `ExpressionEvaluator` | `ExpressionEvaluator(input, param_names)` | Numerical evaluator |
| `Vegas` | `Vegas(n_dims, **opts)` | Monte Carlo integrator |
| `IntegrateResult` | `IntegrateResult(integral, error)` | Integration result |
| `DoubleF64` | `DoubleF64(hi, lo=0.0)` | Double-precision float |
| `GroebnerBasis` | (returned by functions) | Gröbner basis result |
| `RealSolution` | (returned by functions) | Real solution |
| `PolynomialSystemSolution` | (returned by functions) | System solving result |
| `HilbertSeries` | (returned by functions) | Hilbert series |
| `PrimaryComponent` | (returned by functions) | Primary decomposition component |
| `DiophantineSolution` | (returned by functions) | Diophantine equation solution |
| `PolynomialFactor` | (returned by functions) | Factorization factor |
| `AlgebraicFactor` | (returned by functions) | Algebraic-field factorization factor |
