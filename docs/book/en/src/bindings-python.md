# Python API

The `ocas` Python package (built with PyO3) exposes symbolic expressions,
polynomials, matrices, coefficient domains, and numeric evaluation. Install
with `pip install ocas`.

---

## Expression

```python
import ocas

e = ocas.Expression("sin(x)^2 + cos(x)^2")
print(e.simplify())          # 1
print(e.diff("x"))           # 2*cos(x)*sin(x) - 2*sin(x)*cos(x)  (pre-simplify)
print(e.taylor("x", 0, 4))   # Taylor expansion around 0

# Operator overloads
f = ocas.Expression("x^2") + ocas.Expression("2*x") + ocas.Expression("1")
print(f == ocas.Expression("x^2 + 2*x + 1"))  # True
```

## Polynomial

`Polynomial` wraps a dense univariate polynomial. The coefficient domain is
selected by the `domain` argument: `"integer"` (default), `"rational"`, or a
`FiniteField` instance.

```python
# Over the integers (default)
p = ocas.Polynomial([1, 2, 1])     # 1 + 2x + x^2
print(p.coeffs())                   # ['1', '2', '1']
print(p.degree())                   # 2
print(p.eval(2))                    # '9'

q = ocas.Polynomial([1, 1])         # 1 + x
print((p * q).coeffs())             # ['1', '3', '3', '1']

# GCD and factorization
a = ocas.Polynomial([-1, 0, 1])     # x^2 - 1
b = ocas.Polynomial([1, 1])         # x + 1
print(a.gcd(b).coeffs())            # ['1', '1']

for fac, mult in a.square_free_factorization():
    print(fac.coeffs(), mult)

# Over the rationals: pass ints or (num, denom) tuples
r = ocas.Polynomial([(1, 2), 3], domain="rational")  # 1/2 + 3x
print(r.coeffs())                   # ['1/2', '3']

# Over a finite field
gf5 = ocas.FiniteField(5)
fq = ocas.Polynomial([1, 2, 1], domain=gf5)
print(fq.eval(3))                   # '4'  (1 + 6 + 9 = 16 ≡ 4 mod 5)
```

Coefficient values are returned as decimal **strings** to preserve
arbitrary precision across the gmp/non-gmp builds; wrap them in `int(...)`
to obtain Python integers. Rational values are rendered as `n/d`.

## Matrix

```python
m = ocas.Matrix([[1, 2], [3, 4]])
print(m.nrows, m.ncols)             # 2 2
print(m.shape())                    # (2, 2)
print(m[0, 1])                      # '2'
print(m.determinant())              # '-2'
print(m.rank())                     # 2
print(m.transpose().rows())         # [['1', '3'], ['2', '4']]
print(m.trace())                    # '5'

# Arithmetic
a = ocas.Matrix([[1, 2], [3, 5]])   # det = -1, integer inverse exists
inv = a.inverse()
print((a @ inv).rows())             # [['1','0'],['0','1']]

# Solve Ax = b
A = ocas.Matrix([[2, 1], [1, 1]])
print(A.solve([4, 3]))              # ['1', '2']
```

## Domains

```python
ocas.IntegerDomain()        # ℤ
ocas.RationalDomain()       # ℚ
ocas.FiniteField(7)         # GF(7); modulus must be a prime ≥ 2
```

## Numeric evaluation

```python
ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
print(ev.evaluate([3.0, 1.0]))      # [10.0]
print(ev.evaluate([2.0, 0.0]))      # [4.0]
```

## Solvers

```python
# Linear systems over ℚ and ℤ
print(ocas.solve_linear_rational([[1, 1], [1, -1]], [3, 1]))
print(ocas.solve_linear_integer([[2, 1]], [3]))

# Diophantine: a*x + b*y = c
print(ocas.solve_diophantine(3, 5, 1))

# ODE: dsolve with optional hint, Laplace IVP
e = ocas.Expression("Derivative(y(x), x) - y(x)")
print(ocas.classify_ode(e, "y", "x"))     # ['LinearFirst', 'Separable', ...]
print(ocas.dsolve(e, "y", "x"))            # y = C1*exp(x)
print(ocas.dsolve_ivp(e, "y", "x", "1"))   # y = exp(x)
```
## Numeric Integration (Vegas)

```python
import ocas

# 1-D convenience
result = ocas.integrate_1d(lambda x: x**2, 0, 1, n_samples=10000, iterations=10)
print(result.integral, result.error)

# Multi-dimensional
vegas = ocas.Vegas(n_dims=2, n_samples=10000, iterations=10)
result = vegas.integrate(lambda coords: coords[0] * coords[1])
print(result.integral, result.error)
```

`integrate_1d` is a convenience wrapper; `Vegas` supports arbitrary
dimensions. See [Numeric Integration](./numeric-integration.md) for details.

## Automatic Differentiation

```python
from ocas import DualShape, HyperDual, Rational

shape = DualShape.first_order(2)
x = HyperDual.variable(shape, 0, Rational(1))
y = HyperDual.variable(shape, 1, Rational(2))

f = x * y
print(f.value())    # 2
print(f.deriv(0))   # 2  (∂f/∂x)
print(f.deriv(1))   # 1  (∂f/∂y)

# Arithmetic: __add__, __sub__, __mul__, __truediv__, __neg__
g = x + y * x
print(g.deriv(0))   # 3
```

Only `Rational` coefficients are supported. See
[Automatic Differentiation](./autodiff.md).

## Tensors

```python
from ocas import Tensor, contract_tensors, tensor_symmetrise_sign

A = Tensor("A", [("mu", "upper")])
B = Tensor("B", [("mu", "lower")])

# Contract matching indices
result = contract_tensors(A, B)
print(result)  # scalar expression

# Rank-2 tensor with symmetry
g = Tensor("g", [("mu", "upper"), ("nu", "lower")])
print(g.rank())    # 2
print(g.symmetry())  # None

sign = tensor_symmetrise_sign(g)
```

Indices are matched by label and opposite position. See
[Tensors](./tensors.md).

## Algebraic Numbers

```python
from ocas import AlgebraicExtension, AlgebraicPolynomial, FiniteField

# Q(√2)
ext = AlgebraicExtension("x^2 - 2")  # minimal polynomial
alpha = ext.alpha()
print(ext.extension_degree())          # 2

# Factor over Q(√2)
p = AlgebraicPolynomial(ext, [1, 0, -1])  # x^2 - 1
for fac, mult in p.factor():
    print(fac, mult)
```

`AlgebraicExtension` supports extensions over ℚ and GF(p^d). See the
[Factorization](./algorithms/factorization.md) chapter for algorithmic details.

## Gröbner Bases and Ideal Operations

```python
from ocas import MultivariatePolynomial, groebner_basis, ideal_contains
from ocas import solve_polynomial_system, ideal_radical, primary_decomposition
from ocas import hilbert_series, is_zero_dimensional, eliminate

# Create multivariate polynomials via dict format
# x² + y² - 1 in k[x,y]
circle = MultivariatePolynomial({(2, 0): 1, (0, 2): 1, (0, 0): -1}, n_vars=2)

# x - y
line = MultivariatePolynomial({(1, 0): 1, (0, 1): -1}, n_vars=2)

# Compute Gröbner basis
gb = groebner_basis([circle, line], n_vars=2)
print(len(gb))           # 2

# Membership testing
print(ideal_contains([line], circle, n_vars=2))  # True

# Solve polynomial system
sol = solve_polynomial_system([circle, line], n_vars=2)
print(sol.kind)          # "ZeroDimensional"
for s in sol.real_solutions:
    print(s.values)      # [0.707107, 0.707107], [-0.707107, -0.707107]

# Hilbert series
hs = hilbert_series(gb)
print(hs.dimension())    # 0
print(hs.degree())       # 2

# Radical
x2 = MultivariatePolynomial({(2, 0): 1}, n_vars=2)
xy = MultivariatePolynomial({(1, 1): 1}, n_vars=2)
rad = ideal_radical([x2, xy], n_vars=2)  # √(x², xy) = (x)

# Primary decomposition
decomp = primary_decomposition([x2, xy], n_vars=2)
print(len(decomp))       # 2

# Elimination
# Eliminate x from {x² - y, x - z} → {y - z²}
f1 = MultivariatePolynomial({(2, 0): 1, (0, 1): -1}, n_vars=2)  # x² - y
f2 = MultivariatePolynomial({(1, 0): 1, (0, 1): -1}, n_vars=2)  # x - y
result = eliminate([f1, f2], elim_vars=1, n_vars=2)
```

Coefficients can be `int`, `float`, or `(num, denom)` tuples for rationals.
See [Gröbner Bases](./algorithms/groebner.md) for algorithmic details.