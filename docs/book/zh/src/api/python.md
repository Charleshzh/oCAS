# Python API 参考

本章记录 oCAS Python 绑定 (`ocas` 模块) 的完整 API。每个类和函数均包含签名、参数说明、返回值、异常和完整示例。

> **版本**：本文档对应 oCAS 0.24.x。

## 导入

```python
import ocas
```

模块顶层导出所有类和函数，无需子模块导入。

---

## 表达式

### Expression

符号表达式类。从字符串解析构造，支持算术运算、微积分、化简等操作。

**签名**：

```python
ocas.Expression(input: str)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `input` | `str` | 符号表达式字符串，使用 `^` 表示幂次，`*` 表示乘法 |

**返回值**：`Expression` 实例。

**异常**：`ValueError` — 表达式字符串解析失败。

**支持的运算**：

| 运算 | 方法 | 说明 |
|---|---|---|
| 加法 | `__add__` | `e1 + e2` |
| 减法 | `__sub__` | `e1 - e2` |
| 乘法 | `__mul__` | `e1 * e2` |
| 幂 | `__pow__` | `e1 ** e2` |
| 取负 | `__neg__` | `-e` |
| 相等比较 | `__eq__` | `e1 == e2`（基于规范化形式比较） |
| 哈希 | `__hash__` | 可用作字典键或集合元素 |
| 字符串化 | `__str__` | 输出内部表示 |
| 表示 | `__repr__` | `Expression("...")` 格式 |

**方法**：

---

#### Expression.clone

```python
Expression.clone() -> Expression
```

返回表达式的深拷贝。

**示例**：

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

使用内置规则集对表达式进行化简（不动点迭代）。

**返回值**：化简后的 `Expression`。

**异常**：`ValueError` — 化简过程中发生内部错误。

**示例**：

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

对表达式关于变量 `var` 求导。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var` | `str` | 微分变量名 |

**返回值**：导数表达式。

**示例**：

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

对表达式关于变量 `var` 进行符号积分。采用分层积分管线：快速查表 → 有理函数 → Risch → 三角重写+realify → 特殊函数 → 启发式 → 未求值形式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var` | `str` | 积分变量名 |

**返回值**：积分结果表达式。若无法求出闭式解，返回未求值的 `Integral(expr, var)` 形式。

**示例**：

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

仅使用启发式方法积分（分部积分、三角替换等），不调用完整 Risch 算法。成功时返回闭式解，失败时返回未求值形式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var` | `str` | 积分变量名 |

**返回值**：积分结果或未求值 `Integral(expr, var)`。

**示例**：

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

计算表达式在 `point` 处关于 `var` 的 Taylor 展开至 `order` 阶。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var` | `str` | 展开变量名 |
| `point` | `Expression` | 展开中心点 |
| `order` | `int` | 展开阶数（非负整数） |

**返回值**：Taylor 多项式。

**示例**：

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

将表达式中所有 `var` 的出现替换为 `replacement`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var` | `str` | 被替换的变量名 |
| `replacement` | `Expression` | 替换表达式 |

**返回值**：替换后的表达式。

**示例**：

```python
>>> e = ocas.Expression("x^2 + y")
>>> r = ocas.Expression("3")
>>> print(e.substitute("x", r))
y + (3^2)
```

---

## 多项式

### Polynomial

稠密一元多项式，支持整数 ℤ、有理数 ℚ 和有限域 GF(p) 三种系数域。

**签名**：

```python
ocas.Polynomial(coeffs, domain=None)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `coeffs` | `list` | 系数列表，从低次到高次排列 |
| `domain` | `str` 或 `FiniteField` 或 `None` | 系数域：`"integer"`（默认）、`"rational"`、`FiniteField(p)` |

**系数格式**：
- 整数域：整数列表，如 `[1, 2, 1]` 表示 $1 + 2x + x^2$
- 有理数域：整数列表或 `(分子, 分母)` 元组列表
- 有限域：整数列表，运算自动取模

**异常**：`TypeError` — 系数类型不匹配所选域。

**示例**：

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

返回系数列表（字符串形式），从低次到高次。

**示例**：

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

返回多项式的次数。零多项式返回 `None`。

**示例**：

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

返回存储的系数个数。

---

#### Polynomial.is_zero

```python
Polynomial.is_zero() -> bool
```

是否为零多项式。

---

#### Polynomial.eval

```python
Polynomial.eval(x) -> str
```

在点 `x` 处求值。整数域下 `x` 为整数；有理数域下 `x` 为整数或 `(分子, 分母)` 元组。

**返回值**：求值结果的字符串表示。

**示例**：

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

返回形式导数。

**示例**：

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

返回形式积分（常数项为零）。三个系数域均可调用：在有理数域上精确；在整数域和有限域上，当系数除以 $i+1$ 不能整除时该系数被置零（不抛异常）。

---

#### Polynomial.primitive_part

```python
Polynomial.primitive_part() -> Polynomial
```

返回本原部分（去除内容）。仅支持整数域。

**异常**：`ValueError` — 非整数域多项式。

---

#### Polynomial.factor

```python
Polynomial.factor() -> list[PolynomialFactor]
```

完整因式分解。返回 `(因子, 重数)` 对的列表：整数域上每个因子为本原多项式，有限域上为首一多项式。

**因式分解策略**：
- ℤ[x]：无平方分解 + Berlekamp–Zassenhaus + Hensel 提升
- GF(p)[x]：Berlekamp 算法

**示例**：

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

无平方分解。返回互不相同的无平方因子及其重数。

**示例**：

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

是否无平方（无重复不可约因子）。

---

#### Polynomial.gcd

```python
Polynomial.gcd(other: Polynomial) -> Polynomial
```

计算与 `other` 的最大公因式。两个多项式必须在同一系数域上。

**异常**：`TypeError` — 系数域不匹配。

**示例**：

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

带余除法，返回 `(商, 余式)`。`other` 为零时返回 `None`。

**异常**：`TypeError` — 系数域不匹配。

**示例**：

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

#### Polynomial 算术运算

| 运算 | 语法 | 说明 |
|---|---|---|
| 加法 | `p + q` | 多项式加法 |
| 减法 | `p - q` | 多项式减法 |
| 乘法 | `p * q` | 多项式乘法 |
| 取负 | `-p` | 系数取负 |
| 相等 | `p == q` | 归一化系数比较 |

---

### PolynomialFactor

因式分解结果中的单个因子。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `factor` | `Polynomial` | 因子多项式（整数域上本原，有限域上首一） |
| `multiplicity` | `int` | 重数 |

---

### MultivariatePolynomial

稀疏多元多项式（$\mathbb{Q}[x_1, \dots, x_n]$），使用字典序 (Lex)。

**签名**：

```python
ocas.MultivariatePolynomial(terms: dict, n_vars: int)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `terms` | `dict[tuple[int,...], int]` | 键为指数元组（长度为 `n_vars`），值为系数 |
| `n_vars` | `int` | 变量个数 |

**示例**：

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

返回变量个数。

---

## 矩阵

### Matrix

域上的稠密矩阵，支持 ℤ、ℚ 和 GF(p)。

**签名**：

```python
ocas.Matrix(rows, domain=None)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rows` | `list[list]` | 行优先的二维列表 |
| `domain` | `str` 或 `FiniteField` 或 `None` | 系数域：`"integer"`（默认）、`"rational"`、`FiniteField(p)` |

**异常**：`TypeError` — 元素类型不匹配所选域。

**示例**：

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

行数。

---

#### Matrix.ncols

```python
Matrix.ncols -> int  # property
```

列数。

---

#### Matrix.shape

```python
Matrix.shape() -> tuple[int, int]
```

返回 `(nrows, ncols)`。

---

#### Matrix.__getitem__

```python
Matrix[i, j] -> str
```

返回 `(i, j)` 位置的元素（字符串形式）。

**异常**：`ValueError` — 索引越界。

**示例**：

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

返回所有行的二维列表（字符串形式）。

---

#### Matrix.transpose

```python
Matrix.transpose() -> Matrix
```

返回转置矩阵。

---

#### Matrix.trace

```python
Matrix.trace() -> str
```

返回方阵的迹（对角线元素之和）。

**异常**：`ValueError` — 非方阵。

---

#### Matrix.rank

```python
Matrix.rank() -> int
```

返回矩阵的秩。

---

#### Matrix.determinant

```python
Matrix.determinant() -> str
```

返回方阵的行列式（Bareiss 无分数算法）。

**异常**：`ValueError` — 非方阵。

**示例**：

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

返回逆矩阵。

**异常**：`ValueError` — 奇异矩阵或非方阵。

---

#### Matrix.solve

```python
Matrix.solve(rhs) -> list[str]
```

求解线性方程组 `self · x = rhs`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rhs` | `list` | 右端向量 |

**返回值**：解向量（字符串列表）。

**异常**：`ValueError` — 无解或不唯一。

**示例**：

```python
>>> a = Matrix([[2, 0], [0, 3]])
>>> a.solve([4, 9])
['2', '3']
```

---

#### Matrix 算术运算

| 运算 | 语法 | 说明 |
|---|---|---|
| 矩阵乘法 | `a @ b` | `__matmul__` |
| 加法 | `a + b` | 逐元素加法 |
| 减法 | `a - b` | 逐元素减法 |
| 相等 | `a == b` | 逐元素比较 |

---

## 系数域

### IntegerDomain

整数域 ℤ 的选择器类。

**签名**：

```python
ocas.IntegerDomain()
```

> **注意**：`Polynomial` 和 `Matrix` 的 `domain` 参数目前接受字符串（`"integer"` / `"rational"`）或 `FiniteField` 实例；`IntegerDomain()` / `RationalDomain()` 实例尚不能直接作为该参数传递。

**示例**：

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

有理数域 ℚ 的选择器。

**签名**：

```python
ocas.RationalDomain()
```

**示例**：

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

有限域 GF(p)，$p$ 为素数。

**签名**：

```python
ocas.FiniteField(modulus: int)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `modulus` | `int` | 素数模数，$p \geq 2$ |

**返回值**：`FiniteField` 实例。

**异常**：`ValueError` — `modulus < 2`。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `modulus` | `str` | 素数模数的十进制字符串 |

**示例**：

```python
>>> from ocas import FiniteField, Polynomial
>>> gf5 = FiniteField(5)
>>> p = Polynomial([1, 2, 1], domain=gf5)
>>> print(repr(gf5))
FiniteField(5)
```

---

## 求解器

### solve_linear_rational

```python
ocas.solve_linear_rational(a: list[list[int]], b: list[int]) -> list[tuple[int, int]]
```

在 ℚ 上求解线性方程组 $A \mathbf{x} = \mathbf{b}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `list[list[int]]` | 系数矩阵（行优先，整数） |
| `b` | `list[int]` | 右端向量（整数） |

**返回值**：解向量，每个分量为 `(分子, 分母)` 元组。

**异常**：`ValueError` — 无唯一解。

**示例**：

```python
>>> ocas.solve_linear_rational([[2, 0], [0, 3]], [4, 9])
[(2, 1), (3, 1)]
```

---

### solve_linear_integer

```python
ocas.solve_linear_integer(a: list[list[int]], b: list[int]) -> list[int]
```

在 ℤ 上求解线性方程组 $A \mathbf{x} = \mathbf{b}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `list[list[int]]` | 系数矩阵（行优先，整数） |
| `b` | `list[int]` | 右端向量（整数） |

**返回值**：整数解向量。

**异常**：`ValueError` — 无整数解。

**示例**：

```python
>>> ocas.solve_linear_integer([[2, 0], [0, 3]], [4, 9])
[2, 3]
```

---

### solve_diophantine

```python
ocas.solve_diophantine(a: int, b: int, c: int) -> DiophantineSolution | None
```

求解线性丢番图方程 $ax + by = c$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `int` | $x$ 的系数 |
| `b` | `int` | $y$ 的系数 |
| `c` | `int` | 右端常数 |

**返回值**：`DiophantineSolution`（有解时）或 `None`（无整数解时）。

**示例**：

```python
>>> sol = ocas.solve_diophantine(3, 5, 1)
>>> sol.particular
(2, -1)
>>> sol.general
(5, -3)
```

---

### DiophantineSolution

丢番图方程的解。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `particular` | `tuple[int, int]` | 特解 $(x_0, y_0)$ |
| `general` | `tuple[int, int]` | 通解方向 $(t_x, t_y)$；通解为 $(x_0 + k \cdot t_x, \; y_0 + k \cdot t_y)$，$k \in \mathbb{Z}$ |

---

## ODE

### classify_ode

```python
ocas.classify_ode(equation: Expression, func: str, var: str) -> list[str]
```

分类常微分方程，返回可用求解方法名列表。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `Expression` | ODE 表达式（等于零的形式），如 `Derivative(y(x), x) - y(x)` |
| `func` | `str` | 未知函数名，如 `"y"` |
| `var` | `str` | 自变量名，如 `"x"` |

**返回值**：方法名列表，如 `["LinearFirst", "PowerSeries"]`。

**支持的 ODE 类型**：`Separable`、`LinearFirst`、`Bernoulli`、`Exact`、`Homogeneous`、`LinearConstantCoeff`、`CauchyEuler`、`ReductionOfOrder`、`PowerSeries`。

**示例**：

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

符号求解 ODE。

**签名**（含默认值）：`dsolve(equation, func, var, hint=None)`

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `Expression` | ODE 表达式（等于零的形式） |
| `func` | `str` | 未知函数名 |
| `var` | `str` | 自变量名 |
| `hint` | `str` 或 `None` | 指定求解方法（`classify_ode` 返回的名称之一） |

**返回值**：解的字符串表示（显式、隐式、级数或未求值形式）。

**示例**：

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

通过 Laplace 变换求解常系数线性 ODE 的初值问题。

**签名**（含默认值）：`dsolve_ivp(equation, func, var, y0, y1=None)`

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `Expression` | ODE 表达式（等于零的形式） |
| `func` | `str` | 未知函数名 |
| `var` | `str` | 自变量名 |
| `y0` | `str` | $y(0)$ 的值（字符串表达式，如 `"1"`） |
| `y1` | `str` 或 `None` | $y'(0)$ 的值（二阶问题需要） |

**返回值**：不含自由常数的显式解字符串。

**示例**：

```python
>>> eq = ocas.Expression("Derivative(y(x), x, 2) + y(x)")
>>> ocas.dsolve_ivp(eq, "y", "x", "0", "1")
'y = sin(x)'
```

---

## 数论

### factorint

```python
ocas.factorint(n: int) -> list[tuple[str, int]]
```

将 $|n|$ 分解为素因子。负输入的首项为 `("-1", 1)`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 待分解的整数（任意精度） |

**返回值**：`(素数, 指数)` 元组列表，按素数升序排列。

**示例**：

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

BPSW 概率素性测试。对 $n < 2^{64}$ 为确定性判定；目前已知对任意大小无合数通过。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 待测整数（任意精度） |

**返回值**：`True` 表示（可能是）素数。

**示例**：

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

对 $n < 2^{64}$ 的确定性素性判定（u64 快速路径）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 待测整数（必须在 u64 范围内） |

**返回值**：`True` 表示素数。

**示例**：

```python
>>> ocas.isprime_u64(2**61 - 1)
True
```

---

### nextprime

```python
ocas.nextprime(n: int) -> int
```

返回严格大于 $n$ 的最小素数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 起始整数 |

**返回值**：下一个素数（任意精度）。

**示例**：

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

求解离散对数 $\text{base}^x \equiv \text{target} \pmod{p}$。对素数 $p$ 使用 Pohlig–Hellman 算法，否则使用 BSGS。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `p` | `int` | 模数 |
| `base` | `int` | 底数 |
| `target` | `int` | 目标值 |

**返回值**：离散对数 $x$。

**异常**：`ValueError` — 对数不存在。

**示例**：

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

中国剩余定理：求解同余方程组 $x \equiv r_i \pmod{m_i}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `moduli` | `list[int]` | 模数列表（不要求互素） |
| `residues` | `list[int]` | 余数列表 |

**返回值**：`(r, m)` 表示 $x \equiv r \pmod{m}$。

**异常**：`ValueError` — 方程组不相容。

**示例**：

```python
>>> ocas.crt([3, 5, 7], [2, 3, 2])
(23, 105)
```

---

### jacobi_symbol

```python
ocas.jacobi_symbol(a: int, n: int) -> int
```

计算 Jacobi 符号 $\left(\frac{a}{n}\right)$，其中 $n$ 为正奇数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `int` | 分子 |
| `n` | `int` | 分母（正奇数） |

**返回值**：$-1$、$0$ 或 $1$。

**异常**：`ValueError` — $n$ 非正奇数。

**示例**：

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

Euler totient 函数 $\varphi(n)$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 正整数 |

**返回值**：$\varphi(n)$（任意精度）。

**示例**：

```python
>>> ocas.totient(12)
4
```

---

### mobius

```python
ocas.mobius(n: int) -> int
```

Möbius 函数 $\mu(n)$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 正整数 |

**返回值**：$-1$、$0$ 或 $1$。

**示例**：

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

正因数个数 $\tau(n)$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 正整数 |

**返回值**：正因数个数。

**示例**：

```python
>>> ocas.divisor_count(12)
6
```

---

### divisor_sigma

```python
ocas.divisor_sigma(n: int, k: int = 1) -> int
```

因数的 $k$ 次幂之和 $\sigma_k(n) = \sum_{d \mid n} d^k$。

**签名**（含默认值）：`divisor_sigma(n, k=1)`

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 正整数 |
| `k` | `int` | 幂次（默认 1，即因数之和 $\sigma_1$） |

**返回值**：$\sigma_k(n)$。

**示例**：

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

Liouville 函数 $\lambda(n) = (-1)^{\Omega(n)}$，其中 $\Omega(n)$ 是 $n$ 的素因子个数（计重数）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `int` | 正整数 |

**返回值**：$-1$ 或 $1$。

**示例**：

```python
>>> ocas.liouville_lambda(12)  # Ω(12)=3
-1
>>> ocas.liouville_lambda(9)   # Ω(9)=2
1
```

---

## Gröbner 基与理想

> **注意**：以下函数保留 `py_` 前缀（Rust 端未使用 `#[pyo3(name=...)]` 重命名）。在 Python 中的实际调用名为 `ocas.py_groebner_basis`、`ocas.py_ideal_contains` 等。

### py_groebner_basis

```python
ocas.py_groebner_basis(generators, n_vars: int = 1, algorithm: str = "auto") -> GroebnerBasis
```

计算多项式理想的 Gröbner 基。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `list` | 生成元列表（`MultivariatePolynomial` 或 `Polynomial` 对象） |
| `n_vars` | `int` | 变量个数（默认 1；若生成元为 `MultivariatePolynomial` 则自动检测） |
| `algorithm` | `str` | 算法选择：`"auto"`（默认）、`"f4"`、`"f5"`、`"buchberger"` |

**返回值**：`GroebnerBasis` 对象。

**异常**：`ValueError` — 算法名称无效。

**示例**：

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

判断多项式 $f$ 是否属于由 `generators` 生成的理想。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `list` | 理想的生成元列表 |
| `f` | `MultivariatePolynomial` 或 `Polynomial` | 待检测的多项式 |
| `n_vars` | `int` | 变量个数（默认 1） |
| `algorithm` | `str` | 算法选择（默认 `"auto"`） |

**返回值**：`True` 表示 $f \in I$。

**示例**：

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

求解多项式方程组。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equations` | `list` | 方程列表（每个方程为等于零的多项式） |
| `n_vars` | `int` | 变量个数（默认 1） |
| `algorithm` | `str` | 算法选择（默认 `"auto"`） |

**返回值**：`PolynomialSystemSolution` 对象。

**示例**：

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

计算 Gröbner 基对应的 Hilbert 级数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `gb` | `GroebnerBasis` | Gröbner 基 |

**返回值**：`HilbertSeries` 对象。

**示例**：

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

计算理想的根式 $\sqrt{I}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `list` | 理想的生成元列表 |
| `n_vars` | `int` | 变量个数（默认 1） |

**返回值**：$\sqrt{I}$ 的 Gröbner 基。

**算法**：零维使用无平方分解；正维使用 Jacobian 饱和 $\sqrt{I} = I : h^\infty$。

---

### py_primary_decomposition

```python
ocas.py_primary_decomposition(generators, n_vars: int = 1) -> list[PrimaryComponent]
```

计算理想的准素分解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `list` | 理想的生成元列表 |
| `n_vars` | `int` | 变量个数（默认 1） |

**返回值**：`PrimaryComponent` 列表。

**示例**：

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

判断 Gröbner 基对应的理想是否为零维的。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `gb` | `GroebnerBasis` | Gröbner 基 |

**返回值**：`True` 表示零维。

---

### py_eliminate

```python
ocas.py_eliminate(generators, elim_vars: int, n_vars: int = 1, algorithm: str = "auto") -> GroebnerBasis
```

消元：计算消元理想 $I \cap k[x_{e+1}, \dots, x_n]$ 的 Gröbner 基。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `list` | 理想的生成元列表 |
| `elim_vars` | `int` | 消除的变量数（从前向后） |
| `n_vars` | `int` | 总变量个数（默认 1） |
| `algorithm` | `str` | 算法选择（默认 `"auto"`） |

**返回值**：消元理想的 Gröbner 基。

---

### MultivariatePolynomial（Gröbner 辅助类）

已在多项式一节介绍。Gröbner 函数同时接受 `MultivariatePolynomial` 和 `Polynomial` 作为生成元；后者自动映射为单变量多元多项式。

---

### GroebnerBasis

Gröbner 基计算结果。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `n_vars` | `int` | 变量个数 |

**方法**：

| 方法 | 返回值 | 说明 |
|---|---|---|
| `__len__()` | `int` | 基中多项式个数 |
| `is_groebner_basis()` | `bool` | 验证是否确实为 Gröbner 基 |

---

### RealSolution

多项式方程组的一个实数解。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `values` | `list[float]` | 解的坐标值 |
| `multiplicity` | `int` | 重数 |

---

### PolynomialSystemSolution

多项式方程组的求解结果。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `kind` | `str` | `"zero_dimensional"`、`"positive_dimensional"` 或 `"empty"` |
| `vector_space_dimension` | `int` 或 `None` | 向量空间维数（仅零维时有值） |

**方法**：

| 方法 | 返回值 | 说明 |
|---|---|---|
| `solutions()` | `list[RealSolution]` | 实数解列表（仅零维时非空） |

---

### HilbertSeries

Hilbert 级数。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `dimension` | `int` | Krull 维数 |
| `degree` | `int` | 重数（degree） |
| `numerator` | `list[int]` | 分子多项式系数 |

**方法**：

| 方法 | 参数 | 返回值 | 说明 |
|---|---|---|---|
| `hilbert_function(d)` | `d: int` | `int` | 在次数 $d$ 处的 Hilbert 函数值 |

---

### PrimaryComponent

准素分解的一个分量。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `n_vars` | `int` | 变量个数 |

---

## 代数数

### AlgebraicExtension

代数数域 $K = \mathbb{Q}(\alpha)$，由首一极小多项式定义。

**签名**：

```python
ocas.AlgebraicExtension(min_poly: list)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `min_poly` | `list[int]` 或 `list[tuple[int,int]]` | 极小多项式系数（升幂排列），末项（首项系数）必须为 1 |

**异常**：`ValueError` — 系数列表过短、首项系数不为 1 或多项式不可约性验证失败。

**示例**：

```python
>>> from ocas import AlgebraicExtension
>>> # α² - 2（即 ℚ(√2)）
>>> field = AlgebraicExtension([-2, 0, 1])
>>> field.extension_degree()
2
```

---

#### AlgebraicExtension.extension_degree

```python
AlgebraicExtension.extension_degree() -> int
```

返回扩张次数 $[K:\mathbb{Q}] = \deg(m)$。

---

#### AlgebraicExtension.alpha

```python
AlgebraicExtension.alpha() -> AlgebraicElement
```

返回扩张的生成元 $\alpha$。

**示例**：

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

将有理数 $c$ 嵌入为域中的常量元素。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `c` | `int` 或 `tuple[int, int]` | 有理数（整数或 `(分子, 分母)` 元组） |

---

#### AlgebraicExtension.element

```python
AlgebraicExtension.element(coeffs: list) -> AlgebraicElement
```

从 $\alpha$-多项式系数（升幂）创建元素。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `coeffs` | `list[int]` 或 `list[tuple[int,int]]` | 有理系数列表（升幂） |

**示例**：

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> # 3 + 2α
>>> e = field.element([3, 2])
>>> e.coeffs()
['3', '2']
```

---

### AlgebraicElement

代数数域 $\mathbb{Q}(\alpha)$ 中的一个元素，存储为 $\alpha$ 的多项式。

#### AlgebraicElement.coeffs

```python
AlgebraicElement.coeffs() -> list[str]
```

返回 $\alpha$-多项式系数（字符串列表，升幂）。

**示例**：

```python
>>> field = AlgebraicExtension([-2, 0, 1])
>>> e = field.element([1, 1])  # 1 + α = 1 + √2
>>> e.coeffs()
['1', '1']
```

---

### AlgebraicPolynomial

代数数域 $\mathbb{Q}(\alpha)$ 上的稠密一元多项式。

**签名**：

```python
ocas.AlgebraicPolynomial(field: AlgebraicExtension, coeffs: list)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `field` | `AlgebraicExtension` | 系数所在的代数数域 |
| `coeffs` | `list` | 系数列表（升幂），每个元素为 `int`、`(int, int)` 元组或 `AlgebraicElement` |

**示例**：

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

返回次数。零多项式返回 `None`。

---

#### AlgebraicPolynomial.len

```python
AlgebraicPolynomial.len() -> int
```

返回存储的系数个数。

---

#### AlgebraicPolynomial.is_zero

```python
AlgebraicPolynomial.is_zero() -> bool
```

是否为零多项式。

---

#### AlgebraicPolynomial.coeffs

```python
AlgebraicPolynomial.coeffs() -> list[list[str]]
```

返回系数列表。每个系数为一个有理数列表（字符串形式），表示 $\alpha$-多项式系数。

---

#### AlgebraicPolynomial.factor

```python
AlgebraicPolynomial.factor() -> list[AlgebraicFactor]
```

在代数数域上进行因式分解（Trager 算法）。

**返回值**：`AlgebraicFactor` 列表，每个包含不可约因子和重数。

**示例**：

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

代数数域多项式因式分解中的一个因子。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `factor` | `AlgebraicPolynomial` | 不可约因子 |
| `multiplicity` | `int` | 重数 |

---

## 自动微分

### DualShape

一阶对偶数布局描述，声明追踪的微分变量个数。构建一次，多处共享。

**签名**：

```python
# 不直接构造，使用静态方法
shape = DualShape.first_order(n_vars)
```

---

#### DualShape.first_order（静态方法）

```python
DualShape.first_order(n_vars: int) -> DualShape
```

创建包含 `n_vars` 个微分变量的布局。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n_vars` | `int` | 微分变量个数（$\geq 1$） |

**异常**：`ValueError` — `n_vars == 0`。

**示例**：

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

微分变量个数。

---

#### DualShape.n_components

```python
DualShape.n_components -> int  # property
```

总分量数（值 + 导数槽位）。

---

### HyperDual

超对偶数：一个标量值加上对 `DualShape` 中各变量的偏导数。使用有理数精确计算。

**限制**：仅支持多项式/有理算术（`+`、`-`、`*`、`/`、取负）。超越函数（sin/exp/log）未实现；整数幂请用重复乘法。

---

#### HyperDual.variable（静态方法）

```python
HyperDual.variable(shape: DualShape, i: int, value) -> HyperDual
```

创建独立变量 $x_i = \text{value}$（对第 $i$ 个变量的导数为 1，其余为 0）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `DualShape` | 对偶数布局 |
| `i` | `int` | 变量索引 |
| `value` | `int` 或 `tuple[int, int]` | 变量值（整数或有理数元组） |

**异常**：`ValueError` — 索引越界。

---

#### HyperDual.constant（静态方法）

```python
HyperDual.constant(shape: DualShape, value) -> HyperDual
```

创建常量（所有导数为零）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `DualShape` | 对偶数布局 |
| `value` | `int` 或 `tuple[int, int]` | 常量值 |

---

#### HyperDual.value

```python
HyperDual.value() -> str
```

返回标量值（字符串形式，如 `"5"` 或 `"3/7"`）。

---

#### HyperDual.deriv

```python
HyperDual.deriv(i: int) -> str | None
```

返回对第 $i$ 个变量的偏导数。若 `i` 越界返回 `None`。

---

#### HyperDual.n_vars

```python
HyperDual.n_vars -> int  # property
```

微分变量个数。

---

#### HyperDual 算术运算

| 运算 | 语法 | 说明 |
|---|---|---|
| 加法 | `a + b` | 逐分量加法 |
| 减法 | `a - b` | 逐分量减法 |
| 乘法 | `a * b` | 乘积法则 |
| 除法 | `a / b` | 商法则 |
| 取负 | `-a` | 逐分量取负 |

所有运算要求两个操作数共享同一 `DualShape`。不同形状会抛出 `ValueError`。

**完整示例**：

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

## 张量

### Tensor

命名张量对象，带指标槽位和可选对称性。

**签名**：

```python
ocas.Tensor(name: str, slots: list[tuple[str, str]], symmetry: str = "none")
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `str` | 张量名称 |
| `slots` | `list[tuple[str, str]]` | 指标槽位列表，每项为 `(标签, 位置)`。位置取 `"upper"`（逆变）或 `"lower"`（协变） |
| `symmetry` | `str` | 对称性：`"none"`（默认）、`"symmetric"`、`"antisymmetric"` |

**异常**：`ValueError` — 位置或对称性字符串无效。

**示例**：

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

张量名称。

---

#### Tensor.rank

```python
Tensor.rank -> int  # property
```

张量秩（槽位数）。

---

#### Tensor.symmetry

```python
Tensor.symmetry -> str  # property
```

对称性字符串（`"none"`、`"symmetric"`、`"antisymmetric"`）。

---

#### Tensor.slots

```python
Tensor.slots() -> list[tuple[str, str]]
```

返回槽位列表，每项为 `(标签, 位置)`。

---

#### Tensor.dummy_labels

```python
Tensor.dummy_labels() -> list[str]
```

返回当前哑指标（出现两次的指标）标签列表。

---

#### Tensor.to_string_atom

```python
Tensor.to_string_atom() -> str
```

将张量渲染为表达式函数节点 `name(slot, slot, ...)` 的字符串形式。

---

### contract_tensors

```python
ocas.contract_tensors(a: Tensor, b: Tensor) -> tuple[str, list | str]
```

对两个张量进行缩并：对共享哑指标（相同标签、相反变异性）求和。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `Tensor` | 第一个张量 |
| `b` | `Tensor` | 第二个张量 |

**返回值**：`(kind, payload)` 元组：
- `kind = "product"` 时，`payload` 为缩并后张量列表（自由指标拼接）
- `kind = "scalar"` 时，`payload` 为缩并后标量表达式的字符串

**示例**：

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

返回张量的对称化符号（$+1$ 或 $-1$）。`"none"` 和 `"symmetric"` 永远返回 $+1$；`"antisymmetric"` 返回槽位排列的奇偶性。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `tensor` | `Tensor` | 输入张量 |

**返回值**：$+1$ 或 $-1$。

**示例**：

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

使用图同构引擎对张量表达式进行规范化。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `expr` | `str` | 张量表达式字符串 |
| `specs` | `dict[str, str]` | 张量名 → 对称性规格（`"none"`、`"symmetric"`、`"antisymmetric"`） |
| `index_groups` | `dict[str, int]` 或 `None` | 指标维度分组（可选） |

**返回值**：规范化后的表达式字符串。

**异常**：`ValueError` — 解析或规范化失败。

---

### young_project

```python
ocas.young_project(expr: str, tableau: list[int]) -> str
```

对张量表达式应用 Young 投影。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `expr` | `str` | 张量表达式字符串 |
| `tableau` | `list[int]` | Young 图各行长度，如 `[2, 1]` 表示 □□/□ |

**返回值**：投影后的表达式字符串。

**示例**：

```python
>>> proj = ocas.young_project("f(a,b,c)", [1, 1, 1])  # □/□/□ 全反对称投影
>>> "f(a,b,c)" in proj
True
```

---

### refresh_dummies

```python
ocas.refresh_dummies(expr: str, specs: dict[str, str]) -> str
```

重命名（刷新）张量表达式中的哑指标。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `expr` | `str` | 张量表达式字符串 |
| `specs` | `dict[str, str]` | 张量名 → 对称性规格 |

**返回值**：哑指标重命名后的表达式字符串。

---

## 求值

### ExpressionEvaluator

编译表达式以进行快速数值求值。编译一次，多次求值。

**签名**：

```python
ocas.ExpressionEvaluator(input: str, param_names: list[str])
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `input` | `str` | 表达式字符串 |
| `param_names` | `list[str]` | 参数名列表（求值时的顺序） |

**异常**：`ValueError` — 解析或编译错误。

**示例**：

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

用给定参数值求值。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `values` | `list[float]` | 浮点数值列表，长度必须与参数个数匹配 |

**返回值**：结果浮点数列表。

**异常**：`ValueError` — 参数个数不匹配或求值错误。

---

#### ExpressionEvaluator.n_params

```python
ExpressionEvaluator.n_params -> int  # property
```

参数个数。

---

## 数值积分

### Vegas

自适应蒙特卡洛积分器（Vegas 算法），支持多维积分。

**签名**：

```python
ocas.Vegas(n_dims: int, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n_dims` | `int` | 积分维度 |
| `n_bins` | `int` 或 `None` | 每维的 bin 数（默认 64） |
| `n_samples` | `int` 或 `None` | 每次迭代的采样数（默认 10000） |
| `iterations` | `int` 或 `None` | 迭代次数（默认 10） |
| `learning_rate` | `float` 或 `None` | 网格自适应学习率（默认 1.5） |
| `seed` | `int` 或 `None` | RNG 种子（默认 `0x0C45`） |

**示例**：

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

对可调用对象 $f$ 在单位超立方体 $[0, 1]^n$ 上积分。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `Callable` | 接受 `n_dims` 个 `[0, 1]` 内浮点数的列表，返回一个浮点数 |

**返回值**：`IntegrateResult`。

---

#### Vegas.result

```python
Vegas.result -> IntegrateResult  # property
```

最新的累计估计值和误差。

---

#### Vegas.iterations

```python
Vegas.iterations -> int  # property
```

已完成的迭代次数。

---

### integrate_1d

```python
ocas.integrate_1d(f: Callable[[float], float], a: float, b: float, *, n_bins=None, n_samples=None, iterations=None, learning_rate=None, seed=None) -> IntegrateResult
```

一维数值积分便捷函数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `Callable` | 一元浮点函数 |
| `a` | `float` | 积分下限 |
| `b` | `float` | 积分上限 |
| `n_bins` 等 | 同 `Vegas` | 可选调参 |

**返回值**：`IntegrateResult`。

**示例**：

```python
>>> import ocas
>>> r = ocas.integrate_1d(lambda x: x, 0.0, 1.0)
>>> abs(r.integral - 0.5) < 0.01
True
```

---

### IntegrateResult

数值积分结果。

**属性**：

| 属性 | 类型 | 说明 |
|---|---|---|
| `integral` | `float` | 积分估计值 |
| `error` | `float` | 估计标准误差 |

支持索引访问：`result[0]` = `integral`，`result[1]` = `error`；支持解包：`integral, error = result`。

**示例**：

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

## 双精度浮点

### DoubleF64

双精度浮点数运算（~31 位有效数字，~84 二进制位），使用 Dekker/Knuth "double-float" 算法。

**签名**：

```python
ocas.DoubleF64(hi: float, lo: float = 0.0)
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `hi` | `float` | 主值 |
| `lo` | `float` | 误差项（默认 0.0） |

**示例**：

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

转换为标准 `float`（丢失误差项精度）。

---

#### DoubleF64.components

```python
DoubleF64.components() -> tuple[float, float]
```

返回 `(hi, lo)` 元组。

---

#### DoubleF64 算术运算

| 运算 | 语法 | 说明 |
|---|---|---|
| 加法 | `a + b` | 双精度加法 |
| 减法 | `a - b` | 双精度减法 |
| 乘法 | `a * b` | 双精度乘法 |
| 除法 | `a / b` | 双精度除法 |
| 取负 | `-a` | 逐分量取负 |
| 绝对值 | `abs(a)` | 双精度绝对值 |
| 幂 | `a ** n` | 整数幂（`n` 为 `int`） |
| 比较 | `==`, `<`, `<=`, `>`, `>=` | 全序比较 |

**异常**：`ValueError` — 除以零（`__truediv__`）。

---

#### DoubleF64 超越函数

| 方法 | 说明 | 异常 |
|---|---|---|
| `sin()` | 正弦 | — |
| `cos()` | 余弦 | — |
| `tan()` | 正切 | — |
| `exp()` | 自然指数 $e^x$ | — |
| `ln()` | 自然对数 | `ValueError`：非正数 |
| `sqrt()` | 平方根 | `ValueError`：负数 |

**完整示例**：

```python
>>> from ocas import DoubleF64
>>> DoubleF64(0.0).sin().to_f64()
0.0
>>> DoubleF64(4.0).sqrt().components()
(2.0, 0.0)
```

---

## 函数速查表

| 函数 | 签名 | 说明 |
|---|---|---|
| `solve_linear_rational` | `(a, b) -> list[tuple]` | ℚ 上求解 Ax=b |
| `solve_linear_integer` | `(a, b) -> list[int]` | ℤ 上求解 Ax=b |
| `solve_diophantine` | `(a, b, c) -> DiophantineSolution \| None` | ax+by=c |
| `classify_ode` | `(equation, func, var) -> list[str]` | ODE 分类 |
| `dsolve` | `(equation, func, var, hint=None) -> str` | 符号求解 ODE |
| `dsolve_ivp` | `(equation, func, var, y0, y1=None) -> str` | Laplace IVP |
| `factorint` | `(n) -> list[tuple]` | 整数素因子分解 |
| `isprime` | `(n) -> bool` | BPSW 素性测试 |
| `isprime_u64` | `(n) -> bool` | 确定性素性（u64） |
| `nextprime` | `(n) -> int` | 下一个素数 |
| `discrete_log` | `(p, base, target) -> int` | 离散对数 |
| `crt` | `(moduli, residues) -> tuple` | 中国剩余定理 |
| `jacobi_symbol` | `(a, n) -> int` | Jacobi 符号 |
| `totient` | `(n) -> int` | Euler totient |
| `mobius` | `(n) -> int` | Möbius 函数 |
| `divisor_count` | `(n) -> int` | 因数个数 |
| `divisor_sigma` | `(n, k=1) -> int` | 因数幂和 |
| `liouville_lambda` | `(n) -> int` | Liouville 函数 |
| `py_groebner_basis` | `(generators, n_vars=1, algorithm="auto")` | Gröbner 基 |
| `py_ideal_contains` | `(generators, f, n_vars=1, algorithm="auto")` | 理想成员判定 |
| `py_solve_polynomial_system` | `(equations, n_vars=1, algorithm="auto")` | 多项式方程组求解 |
| `py_hilbert_series` | `(gb) -> HilbertSeries` | Hilbert 级数 |
| `py_ideal_radical` | `(generators, n_vars=1)` | 理想根式 |
| `py_primary_decomposition` | `(generators, n_vars=1)` | 准素分解 |
| `py_is_zero_dimensional` | `(gb) -> bool` | 零维判定 |
| `py_eliminate` | `(generators, elim_vars, n_vars=1, algorithm="auto")` | 消元 |
| `contract_tensors` | `(a, b) -> tuple` | 张量缩并 |
| `tensor_symmetrise_sign` | `(tensor) -> int` | 对称化符号 |
| `canonicalize_tensors` | `(expr, specs, index_groups=None)` | 张量规范化 |
| `young_project` | `(expr, tableau) -> str` | Young 投影 |
| `refresh_dummies` | `(expr, specs) -> str` | 哑指标刷新 |
| `integrate_1d` | `(f, a, b, **opts) -> IntegrateResult` | 一维数值积分 |

---

## 类速查表

| 类 | 构造函数 | 说明 |
|---|---|---|
| `Expression` | `Expression(input)` | 符号表达式 |
| `Polynomial` | `Polynomial(coeffs, domain=None)` | 一元多项式 |
| `MultivariatePolynomial` | `MultivariatePolynomial(terms, n_vars)` | 多元多项式 |
| `Matrix` | `Matrix(rows, domain=None)` | 域上矩阵 |
| `IntegerDomain` | `IntegerDomain()` | ℤ 选择器 |
| `RationalDomain` | `RationalDomain()` | ℚ 选择器 |
| `FiniteField` | `FiniteField(modulus)` | GF(p) |
| `AlgebraicExtension` | `AlgebraicExtension(min_poly)` | 代数数域 |
| `AlgebraicElement` | （由 `AlgebraicExtension` 创建） | 域元素 |
| `AlgebraicPolynomial` | `AlgebraicPolynomial(field, coeffs)` | 代数域上多项式 |
| `DualShape` | `DualShape.first_order(n_vars)` | 对偶数布局 |
| `HyperDual` | `HyperDual.variable(...)` / `.constant(...)` | 超对偶数 |
| `Tensor` | `Tensor(name, slots, symmetry="none")` | 命名张量 |
| `ExpressionEvaluator` | `ExpressionEvaluator(input, param_names)` | 数值求值器 |
| `Vegas` | `Vegas(n_dims, **opts)` | Monte Carlo 积分器 |
| `IntegrateResult` | `IntegrateResult(integral, error)` | 积分结果 |
| `DoubleF64` | `DoubleF64(hi, lo=0.0)` | 双精度浮点 |
| `GroebnerBasis` | （由函数返回） | Gröbner 基结果 |
| `RealSolution` | （由函数返回） | 实数解 |
| `PolynomialSystemSolution` | （由函数返回） | 方程组求解结果 |
| `HilbertSeries` | （由函数返回） | Hilbert 级数 |
| `PrimaryComponent` | （由函数返回） | 准素分解分量 |
| `DiophantineSolution` | （由函数返回） | 丢番图方程解 |
| `PolynomialFactor` | （由函数返回） | 因式分解因子 |
| `AlgebraicFactor` | （由函数返回） | 代数域因式分解因子 |
