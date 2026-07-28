# 求解器

oCAS 提供线性方程组、丢番图方程与多项式系统的求解器。本章介绍所有可用求解器及其用法。

---

## ℚ 上的线性方程组

`solve_linear_rational` 在有理数域上求解 $n \times n$ 系统 $Ax = b$。
输入系数为 `i64` 值；解以 `(分子, 分母)` 对返回。

```rust
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
// x = [(2, 1), (1, 1)]  → 2, 1
```

错误：`EmptySystem`、`NonSquare`、`Inconsistent`、`Underdetermined { rank }`。

Python：

```python
print(ocas.solve_linear_rational([[2, 1], [1, -1]], [5, 1]))
# [(2, 1), (1, 1)]
```

---

## ℤ 上的线性方程组

`solve_linear_integer` 求 $Ax = b$ 的整数解。若无整数解则返回错误。

```rust
// 2x + y = 3
let a = vec![vec![2, 1]];
let b = vec![3];
let x = solve_linear_integer(&a, &b).unwrap();
// x = [1, 1]  (2·1 + 1·1 = 3)
```

当解含分数时返回 `ResultNotInDomain` 错误。

---

## 丢番图方程

`solve_diophantine` 求解线性丢番图方程 $a \cdot x + b \cdot y = c$ 的整数解 $x, y$。

```rust
let sol = solve_diophantine(3, 5, 1).unwrap();
// sol = DiophantineSolution { x0: 2, y0: -1, x_step: 5, y_step: -3 }
```

结果给出特解 $(x_0, y_0)$ 和步长值。通解为：

$$
\begin{aligned}
x &= x_0 + x_{step} \cdot t \\
y &= y_0 + y_{step} \cdot t
\end{aligned}
$$

其中 $t$ 为任意整数。

---

## 多项式系统（基于 Gröbner 基）

`solve_polynomial_system` 首先计算 Gröbner 基，然后进行回代，求解多项式方程组。
它使用 Buchberger 算法，支持可配置的单项式序。

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// x + y = 0, x*y - 1 = 0  →  x + y = 0, y^2 + 1 = 0
let eq1 = parse(&ctx, "x + y").unwrap();
let eq2 = parse(&ctx, "x*y - 1").unwrap();
let sol = solve_polynomial_system(&ctx, &[eq1, eq2], &[Symbol::new("x"), Symbol::new("y")]);
```

结果为三角形多项式系统，可通过回代求解。

---

## 常微分方程

`dsolve` 解析求解常微分方程。方程以等于零的表达式、未知函数
（如 `y(x)`）与自变量（如 `x`）给出。

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

`classify_ode` 只返回可用的求解方法而不求解：

```rust
let types = classify_ode(&ctx, ode);
// [LinearFirst, Separable, PowerSeries]
```

### 支持的 ODE 类型

| 类型 | 形式 | 方法 |
|---|---|---|
| 可分离 | $f(x)dx = g(y)dy$ | 直接积分 |
| 一阶线性 | $y' + p(x)y = q(x)$ | 积分因子 $\mu = e^{\int p}$ |
| Bernoulli | $y' + p(x)y = q(x)y^n$ | 替换 $v = y^{1-n}$ |
| 恰当 | $M dx + N dy = 0$，$\partial M/\partial y = \partial N/\partial x$ | 势函数 |
| 齐次 | $y' = f(y/x)$ | 替换 $v = y/x$ |
| 积分因子 | 非恰当一阶 | $\mu(x)$ 或 $\mu(y)$ 检测 |
| 常系数 | $ay'' + by' + cy = f(x)$ | 特征方程 |
| Cauchy-Euler | $ax^2y'' + bxy' + cy = f(x)$ | 指标方程 |
| 降阶法 | 二阶线性 | $y_2 = y_1\int e^{-\int p}/y_1^2$ |
| 常数变易法 | 二阶非齐次 | Wronskian 公式 |
| 待定系数法 | 多项式/指数/三角 forcing | 系数匹配 + 共振 |
| 幂级数 | 常点 | 系数递推 |
| Frobenius | 正则奇点 | 指标方程 + 递推 |
| Laplace IVP | 一阶/二阶线性初值问题 | `dsolve_ivp` |
| 2×2 系统 | $\mathbf{Y}' = A\mathbf{Y}$ | `dsolve_system`（特征分解） |

无法解析求解的 ODE 返回未求值的 `ODE(equation, func)` 形式。

### 初值问题

`dsolve_ivp` 通过 Laplace 变换求解线性常系数初值问题：

```rust
// y' - y = 0, y(0) = 1  =>  y = exp(x)
let sol = dsolve_ivp(&ctx, ode, ctx.num(1), None);
```

### 系统

`dsolve_system` 求解 2×2 常系数系统 $\mathbf{Y}' = A\mathbf{Y}$：

```rust
// y1' = y2, y2' = -y1（谐振子）
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

## 错误

所有求解器返回 `Result<T, SolveError>`。常见错误变体：

| 错误 | 含义 |
|---|---|
| `EmptySystem` | 未提供方程 |
| `NonLinear` | 系统对请求变量非线程 |
| `NonSquare` | 方程数与未知数个数不匹配 |
| `Inconsistent` | 无解 |
| `Underdetermined { rank }` | 无穷多解 |
| `ResultNotInDomain` | 解含分数但要求整数 |

---

## 参见

- [Rust API](./rust-api.md) — 系数域类型与多项式操作
- [重写与化简](./rewrite.md) — 化简求解结果
- [基准与性能对比](./performance.md) — Gröbner 基基准结果
