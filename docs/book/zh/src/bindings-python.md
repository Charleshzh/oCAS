# Python API

`ocas` Python 包（基于 PyO3）提供符号表达式、多项式、矩阵、系数域与数值求值。用 `pip install ocas` 安装。

---

## 表达式

```python
import ocas

e = ocas.Expression("sin(x)^2 + cos(x)^2")
print(e.simplify())          # 1
print(e.diff("x"))           # 2*cos(x)*sin(x) - 2*sin(x)*cos(x)  (化简前)
print(e.taylor("x", 0, 4))   # 在 0 处的 Taylor 展开

# 运算符重载
f = ocas.Expression("x^2") + ocas.Expression("2*x") + ocas.Expression("1")
print(f == ocas.Expression("x^2 + 2*x + 1"))  # True
```

## 多项式

`Polynomial` 封装稠密一元多项式。系数域由 `domain` 参数选择：`"integer"`（默认）、`"rational"` 或 `FiniteField` 实例。

```python
# 整数系数（默认）
p = ocas.Polynomial([1, 2, 1])     # 1 + 2x + x^2
print(p.coeffs())                   # ['1', '2', '1']
print(p.degree())                   # 2
print(p.eval(2))                    # '9'

q = ocas.Polynomial([1, 1])         # 1 + x
print((p * q).coeffs())             # ['1', '3', '3', '1']

# GCD 与因式分解
a = ocas.Polynomial([-1, 0, 1])     # x^2 - 1
b = ocas.Polynomial([1, 1])         # x + 1
print(a.gcd(b).coeffs())            # ['1', '1']

for fac, mult in a.square_free_factorization():
    print(fac.coeffs(), mult)

# 有理系数：传入整数或 (分子, 分母) 元组
r = ocas.Polynomial([(1, 2), 3], domain="rational")  # 1/2 + 3x
print(r.coeffs())                   # ['1/2', '3']

# 有限域
gf5 = ocas.FiniteField(5)
fq = ocas.Polynomial([1, 2, 1], domain=gf5)
print(fq.eval(3))                   # '4'  (1 + 6 + 9 = 16 ≡ 4 mod 5)
```

系数以十进制**字符串**返回，以在 gmp/非 gmp 构建间保持任意精度；用 `int(...)` 转换为 Python 整数。有理数值以 `n/d` 形式表示。

## 矩阵

```python
m = ocas.Matrix([[1, 2], [3, 4]])
print(m.nrows, m.ncols)             # 2 2
print(m.shape())                    # (2, 2)
print(m[0, 1])                      # '2'
print(m.determinant())              # '-2'
print(m.rank())                     # 2
print(m.transpose().rows())         # [['1', '3'], ['2', '4']]
print(m.trace())                    # '5'

# 运算
a = ocas.Matrix([[1, 2], [3, 5]])   # det = -1，整数逆存在
inv = a.inverse()
print((a @ inv).rows())             # [['1','0'],['0','1']]

# 解 Ax = b
A = ocas.Matrix([[2, 1], [1, 1]])
print(A.solve([4, 3]))              # ['1', '2']
```

## 系数域

```python
ocas.IntegerDomain()        # ℤ
ocas.RationalDomain()       # ℚ
ocas.FiniteField(7)         # GF(7)；模数必须是 ≥ 2 的素数
```

## 数值求值

```python
ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
print(ev.evaluate([3.0, 1.0]))      # [10.0]
print(ev.evaluate([2.0, 0.0]))      # [4.0]
```

## 求解器

```python
# ℚ 与 ℤ 上的线性方程组
print(ocas.solve_linear_rational([[1, 1], [1, -1]], [3, 1]))
print(ocas.solve_linear_integer([[2, 1]], [3]))

# 丢番图方程：a*x + b*y = c
print(ocas.solve_diophantine(3, 5, 1))

# 常微分方程：dsolve（可选 hint）、Laplace 初值问题
e = ocas.Expression("Derivative(y(x), x) - y(x)")
print(ocas.classify_ode(e, "y", "x"))     # ['LinearFirst', 'Separable', ...]
print(ocas.dsolve(e, "y", "x"))            # y = C1*exp(x)
print(ocas.dsolve_ivp(e, "y", "x", "1"))   # y = exp(x)
```

## 数值积分（Vegas）

```python
import ocas

# 一维便捷函数
result = ocas.integrate_1d(lambda x: x**2, 0, 1, n_samples=10000, iterations=10)
print(result.integral, result.error)

# 多维积分
vegas = ocas.Vegas(n_dims=2, n_samples=10000, iterations=10)
result = vegas.integrate(lambda coords: coords[0] * coords[1])
print(result.integral, result.error)
```

`integrate_1d` 是便捷函数；`Vegas` 支持任意维度。详见
[数值积分](./numeric-integration.md)。

## 自动微分

```python
from ocas import DualShape, HyperDual, Rational

shape = DualShape.first_order(2)
x = HyperDual.variable(shape, 0, Rational(1))
y = HyperDual.variable(shape, 1, Rational(2))

f = x * y
print(f.value())    # 2
print(f.deriv(0))   # 2  (∂f/∂x)
print(f.deriv(1))   # 1  (∂f/∂y)

# 算术：__add__、__sub__、__mul__、__truediv__、__neg__
g = x + y * x
print(g.deriv(0))   # 3
```

仅支持 `Rational` 系数。详见[自动微分](./autodiff.md)。

## 张量

```python
from ocas import Tensor, contract_tensors, tensor_symmetrise_sign

A = Tensor("A", [("mu", "upper")])
B = Tensor("B", [("mu", "lower")])

# 缩并匹配指标
result = contract_tensors(A, B)
print(result)  # 标量表达式

# 2 阶张量与对称性
g = Tensor("g", [("mu", "upper"), ("nu", "lower")])
print(g.rank())    # 2
print(g.symmetry())  # None

sign = tensor_symmetrise_sign(g)
```

指标按标签和相反位置匹配。详见[张量](./tensors.md)。

## 代数数

```python
from ocas import AlgebraicExtension, AlgebraicPolynomial, FiniteField

# Q(√2)
ext = AlgebraicExtension("x^2 - 2")  # 极小多项式
alpha = ext.alpha()
print(ext.extension_degree())          # 2

# 在 Q(√2) 上因式分解
p = AlgebraicPolynomial(ext, [1, 0, -1])  # x^2 - 1
for fac, mult in p.factor():
    print(fac, mult)
```

`AlgebraicExtension` 支持 ℚ 和 GF(p^d) 上的扩域。算法细节见
[因式分解](./algorithms/factorization.md) 章节。
