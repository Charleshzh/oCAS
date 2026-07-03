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
```
