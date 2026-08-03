# Rust API 参考：求解器

本章涵盖 oCAS 中的方程求解功能，分为三大模块：线性方程组与丢番图方程（`ocas-calc::solve`）、多项式系统求解（`ocas-poly::ideal`）、以及常微分方程求解（`ocas-calc::ode`）。

---

## 线性方程组

### solve_linear_rational

**签名**：`pub fn solve_linear_rational(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<(i64, i64)>, SolveError>`

**功能**：在有理数域 $\mathbb{Q}$ 上求解线性方程组 $A\mathbf{x} = \mathbf{b}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `&[Vec<i64>]` | $n \times n$ 系数矩阵，行优先存储，每个 `Vec<i64>` 为一行 |
| `b` | `&[i64]` | 长度为 $n$ 的右端向量 |

**返回值**：`Result<Vec<(i64, i64)>, SolveError>`
- `Ok(vec)`：解向量，每个元素 `(num, den)` 表示有理数 $\frac{\text{num}}{\text{den}}$，保证 $\text{den} > 0$ 且 $\gcd(|\text{num}|, \text{den}) = 1$
- `Err(SolveError)`：见下方错误列表

**错误**：
- `SolveError::EmptySystem` — $A$ 或 $\mathbf{b}$ 为空
- `SolveError::NonSquare` — $A$ 不是方阵（或 $\mathbf{b}$ 长度 ≠ 行数）
- `SolveError::Inconsistent` — 方程组无解（增广矩阵秩 > 系数矩阵秩）
- `SolveError::Underdetermined { rank }` — 方程组有无穷多解（秩不足）
- `SolveError::ResultNotInDomain` — 解不在目标数域中
- `SolveError::Matrix(_)` — 内部矩阵运算失败（经 `From<MatrixError>` 转换）
- `SolveError::Other(_)` — 其他错误

**示例**：

```rust
use ocas_calc::solve::solve_linear_rational;

// 求解：
//   2x + y = 5
//   x - y = 1
// → x = 2, y = 1
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
assert_eq!(x, vec![(2, 1), (1, 1)]);
```

**参见**：[`solve_linear_integer`](#solve_linear_integer)、[`Matrix::solve`](./rust-matrix.md#solve)

---

### solve_linear_integer

**签名**：`pub fn solve_linear_integer(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<i64>, SolveError>`

**功能**：在整数环 $\mathbb{Z}$ 上求解线性方程组 $A\mathbf{x} = \mathbf{b}$。要求解恰好为整数，否则返回错误。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `&[Vec<i64>]` | $n \times n$ 系数矩阵，行优先存储 |
| `b` | `&[i64]` | 长度为 $n$ 的右端向量 |

**返回值**：`Result<Vec<i64>, SolveError>`
- `Ok(vec)`：整数解向量
- `Err(SolveError)`：错误类型同 [`solve_linear_rational`](#solve_linear_rational)，额外可能返回 `ResultNotInDomain`（解含分数）

**错误**：同 [`solve_linear_rational`](#solve_linear_rational)。

**示例**：

```rust
use ocas_calc::solve::solve_linear_integer;

// 求解：
//   x + y = 3
//   x - y = 1
// → x = 2, y = 1
let a = vec![vec![1, 1], vec![1, -1]];
let b = vec![3, 1];
let x = solve_linear_integer(&a, &b).unwrap();
assert_eq!(x, vec![2, 1]);
```

**参见**：[`solve_linear_rational`](#solve_linear_rational)

---

## 丢番图方程

### DiophantineSolution

**签名**：
```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiophantineSolution {
    /// 特解 (x0, y0)
    pub particular: (i64, i64),
    /// 齐次解步长 (x_step, y_step)
    pub general: (i64, i64),
}
```

**功能**：线性丢番图方程 $ax + by = c$ 的解的表示。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `particular` | `(i64, i64)` | 一组特解 $(x_0, y_0)$ |
| `general` | `(i64, i64)` | 齐次方程 $ax + by = 0$ 的基解 $(x_{\text{step}}, y_{\text{step}})$，通解为 $(x_0 + k \cdot x_{\text{step}},\; y_0 + k \cdot y_{\text{step}})$，$k \in \mathbb{Z}$ |

---

### solve_diophantine

**签名**：`pub fn solve_diophantine(a: i64, b: i64, c: i64) -> Option<DiophantineSolution>`

**功能**：求解线性丢番图方程 $ax + by = c$。使用扩展 Euclid 算法，当且仅当 $\gcd(a, b) \mid c$ 时有解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `i64` | $x$ 的系数 |
| `b` | `i64` | $y$ 的系数 |
| `c` | `i64` | 右端常数 |

**返回值**：`Option<DiophantineSolution>`
- `Some(sol)`：解存在时返回特解和通解结构
- `None`：$\gcd(a, b) \nmid c$，方程无整数解

**错误**：无返回 `Result`；不可解时返回 `None`。

**示例**：

```rust
use ocas_calc::solve::solve_diophantine;

// 求解 3x + 5y = 1
// 特解：x = 2, y = -1（因为 3×2 + 5×(-1) = 1）
// 通解：x = 2 + 5k, y = -1 - 3k
let sol = solve_diophantine(3, 5, 1).unwrap();
assert_eq!(sol.particular, (2, -1));
assert_eq!(sol.general, (5, -3));

// 无解的情况：2x + 4y = 3（gcd(2,4)=2 不整除 3）
assert!(solve_diophantine(2, 4, 3).is_none());
```

**参见**：[`solve_linear_integer`](#solve_linear_integer)

---

## 多项式系统求解

### RealSolution

**签名**：
```rust
#[derive(Debug, Clone)]
pub struct RealSolution {
    /// 每个变量的值（按变量顺序）
    pub values: Vec<f64>,
    /// 代数重数
    pub multiplicity: usize,
}
```

**功能**：多项式方程组的一个实数解（数值近似）。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `values` | `Vec<f64>` | 变量值 $[x_1, x_2, \dots, x_n]$，按原方程组中的变量顺序排列 |
| `multiplicity` | `usize` | 解的代数重数（当前实现中固定为 1，尚未计算实际重数） |

---

### ZeroDimSolutions

**签名**：
```rust
#[derive(Debug, Clone)]
pub struct ZeroDimSolutions {
    /// 所有实数解
    pub solutions: Vec<RealSolution>,
    /// 商环 k[x₁,...,xₙ]/I 的向量空间维数
    pub vector_space_dimension: usize,
}
```

**功能**：零维多项式方程组的求解结果。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `solutions` | `Vec<RealSolution>` | 所有实数解（可能包含近似值） |
| `vector_space_dimension` | `usize` | 商环 $k[x_1, \dots, x_n]/I$ 的向量空间维数，等于 Lex GB 中各变量一元多项式次数之积 |

---

### PolynomialSystemSolution

**签名**：
```rust
#[derive(Debug, Clone)]
pub enum PolynomialSystemSolution {
    /// 零维系统：有限个解
    ZeroDimensional(ZeroDimSolutions),
    /// 正维系统：解集有正维数分量，返回 Gröbner 基
    PositiveDimensional(GroebnerBasis<RationalDomain, Lex>),
    /// 矛盾系统（理想为 ⟨1⟩），无解
    Empty,
}
```

**功能**：多项式方程组求解结果的枚举类型。系统自动判定维数并选择合适的求解策略。

**变体**：

| 变体 | 说明 |
|---|---|
| `ZeroDimensional(z)` | 系统有有限个解，`z` 包含所有实数解和向量空间维数 |
| `PositiveDimensional(gb)` | 系统的解集有正维分量（无穷多解），返回 Lex 序下的 Gröbner 基 |
| `Empty` | 方程组矛盾（Gröbner 基为 {1}），无解。注意：传入空方程组时返回的是 `PositiveDimensional`（空基），而非 `Empty` |

---

### solve_polynomial_system

**签名**：
```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

**功能**：求解多项式方程组 $f_1 = f_2 = \cdots = f_m = 0$。自动计算 Lex 序 Gröbner 基，判定系统维数，并对零维系统通过三角分解回代求出所有实数解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equations` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 多项式方程组，每个多项式表示一个方程 $f_i = 0$ |
| `algo` | `Algorithm` | Gröbner 基算法选择：`Algorithm::Auto`（自动）、`Algorithm::Buchberger`、`Algorithm::F4`、`Algorithm::F5` |

**返回值**：[`PolynomialSystemSolution`](#polynomialsystemsolution)

**错误**：不返回 `Result`；通过枚举变体表达所有可能结果（零维解、正维基、空系统）。

**算法说明**：
1. 计算方程组生成理想的 Gröbner 基（使用指定算法）。输入多项式本身已按 Lex 序类型化（`SparseMultivariatePolynomial<RationalDomain, Lex>`），因此 GB 直接在 Lex 序下计算，无需换序
2. 若 GB = {1}，返回 `Empty`（方程组矛盾、无解）
3. 判定是否零维：每个变量 $x_i$ 在 GB 中是否有纯幂首单项式 $x_i^{N_i}$
4. **零维**：提取 Lex GB 的三角结构，从最后一个变量开始逐个回代（Sturm 实根隔离 + 二分精化），对每个变量求解一元多项式的实根；商环维数按各变量一元多项式次数之积估计
5. **正维**：返回 Lex GB 供进一步分析

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};
use ocas_poly::ideal::{solve_polynomial_system, PolynomialSystemSolution};

let d = RationalDomain;

// 方程组：
//   x² + y² - 1 = 0  （单位圆）
//   x - y = 0         （直线 y = x）
// 解：(±1/√2, ±1/√2)
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
        // 两个解：约 (0.707, 0.707) 和 (-0.707, -0.707)
    }
    _ => panic!("期望零维系统"),
}
```

**参见**：[`GroebnerBasis`](./rust-groebner.md)、[`Algorithm`](./rust-groebner.md#algorithm)、[`is_zero_dimensional`](./rust-groebner.md#is_zero_dimensional)

---

## 常微分方程（ODE）求解

### ODE

**签名**：
```rust
#[derive(Debug, Clone, Copy)]
pub struct ODE<'a> {
    /// 方程，标准形式 lhs - rhs = 0
    pub equation: Atom<'a>,
    /// 未知函数，例如 y(x)
    pub func: Atom<'a>,
    /// 自变量，例如 x
    pub var: Symbol,
}
```

**功能**：描述一个常微分方程。方程以 `equation = 0` 的标准形式存储，其中 `func` 是未知函数（如 $y(x)$），`var` 是自变量（如 $x$）。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `equation` | `Atom<'a>` | 方程的规范化形式（`lhs - rhs`），已归一化 |
| `func` | `Atom<'a>` | 未知函数，例如 `y(x)` |
| `var` | `Symbol` | 自变量符号，例如 `"x"` |

**构造示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::ODE;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);

// 表示 y' - y = 0
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
```

---

### ODESolution

**签名**：
```rust
#[derive(Debug, Clone, Copy)]
pub enum ODESolution<'a> {
    /// 显式解 y = expr
    Explicit(Atom<'a>),
    /// 隐式解 F(x, y) = 0
    Implicit(Atom<'a>),
    /// 参数解 (x(t), y(t))
    Parametric(Atom<'a>, Atom<'a>),
    /// 级数解（截断表达式 + 项数）
    Series(Atom<'a>, usize),
    /// 系统解（各分量）
    System(&'a [Atom<'a>]),
    /// 未能求解，返回原 ODE
    Unsolved(ODE<'a>),
}
```

**功能**：ODE 求解结果的枚举类型。

**变体**：

| 变体 | 说明 |
|---|---|
| `Explicit(expr)` | 显式解 $y = \text{expr}$，含自由常数 C1, C2, ... |
| `Implicit(expr)` | 隐式解 $F(x, y) = 0$（存储表达式 $F$） |
| `Parametric(x(t), y(t))` | 参数解：两个表达式分别给出 $x(t)$ 和 $y(t)$ |
| `Series(expr, n)` | 级数解，`expr` 为截断后的级数表达式，`n` 为项数 |
| `System(components)` | 系统解，`components[i]` 为第 $i$ 个未知函数的解 |
| `Unsolved(ode)` | 未能找到解析解，返回原 ODE |

---

### ODEType

**签名**：
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

**功能**：ODE 类型枚举，决定使用哪种求解方法。

**变体**：

| 变体 | 方程形式 | 说明 |
|---|---|---|
| `Separable` | $g(y)\,y' = f(x)$ | 可分离变量：$\int g(y)\,dy = \int f(x)\,dx + C$ |
| `LinearFirst` | $y' + p(x)\,y = q(x)$ | 一阶线性：积分因子 $\mu = e^{\int p\,dx}$ |
| `Bernoulli` | $y' + p(x)\,y = q(x)\,y^n$ | Bernoulli 方程：替换 $v = y^{1-n}$ 化为线性 |
| `Exact` | $M\,dx + N\,dy = 0$，$\frac{\partial M}{\partial y} = \frac{\partial N}{\partial x}$ | 恰当方程：势函数 $F(x,y) = C$ |
| `Homogeneous` | $y' = f(y/x)$ | 齐次方程：替换 $v = y/x$ 化为可分离 |
| `LinearConstantCoeff` | $a\,y'' + b\,y' + c\,y = f(x)$ | 常系数线性：特征方程 + 待定系数法（多项式/指数/三角强迫项） |
| `CauchyEuler` | $a\,x^2\,y'' + b\,x\,y' + c\,y = f(x)$ | Cauchy-Euler：替换 $x = e^t$ 化为常系数 |
| `ReductionOfOrder` | $a(x)\,y'' + b(x)\,y' + c(x)\,y = f(x)$ | 降阶法：尝试简单候选解 $y_1$，$y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$，非零强迫项用参数变易求特解 |
| `PowerSeries` | 线性 ODE，常点 $x_0 = 0$ | 幂级数法：$y = \sum_{n=0}^{N-1} a_n x^n$；失败后自动回退到 Frobenius 方法 |

> 注：变易常数（variation of parameters）与 Frobenius 方法在内部作为求解例程使用（前者用于降阶法的特解，后者作为幂级数失败后的回退），但它们**不是** `ODEType` 的独立变体。

---

### classify_ode

**签名**：`pub fn classify_ode<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>) -> Vec<ODEType>`

**功能**：分析 ODE 并返回所有适用的求解方法，按优先级排列。实际返回顺序（实现于 `classify.rs`）：一阶方程依次检查 `LinearFirst` → `Bernoulli` → `Separable` → `Exact` → `Homogeneous`；二阶（及以上）线性方程依次检查 `LinearConstantCoeff` → `CauchyEuler` → `ReductionOfOrder`（仅二阶）；最后，任何线性方程（阶 ≥ 1）都会追加 `PowerSeries` 作为兜底。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `ode` | `ODE<'a>` | 待分类的 ODE |

**返回值**：`Vec<ODEType>` — 所有适用的求解方法列表（可能为空）。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_calc::ode::{classify_ode, ODE, ODEType};

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);

// y' - y = 0 是一阶线性
let dy = ctx.fun("Derivative", &[y, x]);
let eq = ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]);
let ode = ODE { equation: eq, func: y, var: Symbol::new("x") };
let methods = classify_ode(&ctx, ode);
assert!(methods.contains(&ODEType::LinearFirst));
```

**参见**：[`dsolve`](#dsolve)、[`ODEType`](#odetype)

---

### dsolve

**签名**：`pub fn dsolve<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>, hint: Option<ODEType>) -> ODESolution<'a>`

**功能**：求解常微分方程。自动分类 ODE 类型并按优先级尝试各种方法；也可通过 `hint` 指定特定方法。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `ode` | `ODE<'a>` | 待求解的 ODE |
| `hint` | `Option<ODEType>` | 可选：指定求解方法。`None` 时自动分类并按优先级尝试 |

**返回值**：[`ODESolution<'a>`](#odesolution)

**求解策略**（`hint = None` 时）：
1. 调用 [`classify_ode`](#classify_ode) 获取候选方法列表
2. 按优先级依次尝试每种方法
3. 第一个成功的方法的解被返回
4. 所有方法均失败时返回 `ODESolution::Unsolved`

**一阶求解器**：
- **可分离**：$\int g(y)\,dy = \int f(x)\,dx + C$
- **一阶线性**：积分因子 $\mu(x) = e^{\int p(x)\,dx}$，通解 $y = \frac{1}{\mu}\left(\int \mu\,q\,dx + C\right)$
- **Bernoulli**：替换 $v = y^{1-n}$ 化为一阶线性
- **恰当方程**：求势函数 $F$ 使 $\frac{\partial F}{\partial x} = M$，$\frac{\partial F}{\partial y} = N$；非恰当时尝试积分因子 $\mu(x)$ 或 $\mu(y)$
- **齐次**：替换 $v = y/x$ 化为可分离

**二阶求解器**：
- **常系数**：特征方程 $ar^2 + br + c = 0$，根据判别式 $\Delta = b^2 - 4ac$ 构造基本解集，待定系数法求特解
- **Cauchy-Euler**：替换 $x = e^t$ 转化为常系数方程，指标方程 $ar^2 + (b-a)r + c = 0$
- **降阶法**：尝试简单候选解（$1, x, x^2, e^x, e^{-x}, e^{2x}$），找到后用 $y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$ 构造第二解
- **参数变易**：已知齐次基本解 $y_1, y_2$ 时，用 Wronskian $W = y_1 y_2' - y_1' y_2$ 求特解

**级数求解器**（`PowerSeries` 分派，固定展开点 $x_0 = 0$，8 项）：
- **幂级数**：在常点 $x_0 = 0$ 展开 $y = \sum_{n=0}^{N-1} a_n x^n$，代入 ODE 递推求系数
- **Frobenius**（内部回退）：若幂级数在常点不适用，自动回退到正则奇点方法，展开 $y = x^r \sum_{n=0}^{N-1} a_n x^n$，指标方程确定 $r$

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::ode::{dsolve, ODE, ODESolution};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.fun("y", &[x]);
let dy = ctx.fun("Derivative", &[y, x]);

// 求解 y' - y = 0 → C1*exp(x)
let ode = ODE {
    equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
    func: y,
    var: Symbol::new("x"),
};
let sol = dsolve(&ctx, ode, None);
assert!(matches!(sol, ODESolution::Explicit(_)));

// 指定方法
let sol_hint = dsolve(&ctx, ode, Some(ODEType::LinearFirst));
assert!(matches!(sol_hint, ODESolution::Explicit(_)));
```

**参见**：[`classify_ode`](#classify_ode)、[`dsolve_ivp`](#dsolve_ivp)、[`dsolve_system`](#dsolve_system)

---

### dsolve_ivp

**签名**：
```rust
pub fn dsolve_ivp<'a>(
    ctx: &'a AtomArena<'a>,
    ode: ODE<'a>,
    y0: Atom<'a>,
    y1: Option<Atom<'a>>,
) -> ODESolution<'a>
```

**功能**：用 Laplace 变换求解一阶或二阶线性常系数 ODE 的初值问题（IVP）。结果为不含自由常数的显式解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `ode` | `ODE<'a>` | ODE（须为线性常系数） |
| `y0` | `Atom<'a>` | 初始条件 $y(0)$ |
| `y1` | `Option<Atom<'a>>` | 二阶问题时为 $y'(0)$；一阶时忽略 |

**返回值**：[`ODESolution<'a>`](#odesolution)
- `Explicit(expr)`：不含自由常数的显式解
- `Unsolved`：Laplace 变换法不适用或无法求逆

**方法说明**：
1. 对 ODE 两端取 Laplace 变换，利用 $\mathcal{L}\{y'\} = sY - y(0)$ 和 $\mathcal{L}\{y''\} = s^2 Y - sy(0) - y'(0)$
2. 代入初始条件，得到 $Y(s)$ 的代数方程
3. 解出 $Y(s)$
4. 通过部分分式分解 + 标准对查表进行逆 Laplace 变换

**支持的强迫项**：多项式、指数 $e^{kx}$、正弦/余弦 $\sin(\omega x)$/$\cos(\omega x)$，以及 $e^{kx}\sin(\omega x)$/$e^{kx}\cos(\omega x)$ 的线性组合。

**示例**：

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

**参见**：[`dsolve`](#dsolve)

---

### dsolve_system

**签名**：
```rust
pub fn dsolve_system<'a>(
    ctx: &'a AtomArena<'a>,
    equations: &[Atom<'a>],
    funcs: &[Atom<'a>],
    var: Symbol,
) -> ODESolution<'a>
```

**功能**：求解 $2 \times 2$ 常系数线性 ODE 系统 $\mathbf{Y}' = A\mathbf{Y}$。通过特征值分解得到通解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `equations` | `&[Atom<'a>]` | 方程列表，每个方程形式为 `Derivative(y_i, x) - (a_i1*y1 + a_i2*y2) = 0` |
| `funcs` | `&[Atom<'a>]` | 未知函数列表，如 `[y1(x), y2(x)]` |
| `var` | `Symbol` | 自变量符号 |

**返回值**：[`ODESolution<'a>`](#odesolution)
- `System(&[Atom])`：各分量的通解（含自由常数 C1, C2）
- `Unsolved`：系统不支持（非 2×2 或特征值无闭式解）

**支持情况**：
- 不同实特征值 $\lambda_1 \neq \lambda_2$：$\mathbf{Y} = C_1 \mathbf{v}_1 e^{\lambda_1 x} + C_2 \mathbf{v}_2 e^{\lambda_2 x}$
- 重复实特征值（含广义特征向量）：$\mathbf{Y} = (C_1 \mathbf{v} + C_2(\mathbf{w} + x\mathbf{v}))e^{\lambda x}$
- 共轭复特征值 $\alpha \pm \beta i$：实值基本解 $e^{\alpha x}(\mathbf{p}\cos\beta x - \mathbf{q}\sin\beta x)$ 等

**示例**：

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

// 谐振子：y1' = y2, y2' = -y1
let eq1 = ctx.add(&[dy1, ctx.mul(&[ctx.num(-1), y2])]);
let eq2 = ctx.add(&[dy2, y1]);
let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
assert!(matches!(sol, ODESolution::System(_)));
```

**参见**：[`dsolve`](#dsolve)

---

## SolveError

**签名**：
```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolveError {
    /// 方程组为空
    EmptySystem,
    /// 方程组非线性
    NonLinear,
    /// 系数矩阵不是方阵（或方程数 ≠ 未知数个数）
    NonSquare,
    /// 方程组矛盾（无解）
    Inconsistent,
    /// 方程组欠定（无穷多解）
    Underdetermined { rank: usize },
    /// 结果不在目标域中
    ResultNotInDomain,
    /// 内部矩阵运算失败
    Matrix(MatrixError),
    /// 带描述的其他错误
    Other(String),
}
```

**功能**：方程求解过程中可能出现的错误。实现 `Display` 和 `Error` trait。实现 `From<MatrixError>`：`ShapeMismatch` → `NonSquare`、`Inconsistent` → `Inconsistent`、`Underdetermined { rank }` → `Underdetermined { rank }`、`ResultNotInDomain` → `ResultNotInDomain`、`RightHandSideIsNotVector` → `Other("right-hand side is not a vector")`。

**变体**：

| 变体 | 说明 |
|---|---|
| `EmptySystem` | 输入方程组为空 |
| `NonLinear` | 方程组在目标变量中非线性 |
| `NonSquare` | 系数矩阵不是方阵（方程数 ≠ 未知数个数） |
| `Inconsistent` | 增广矩阵的秩大于系数矩阵的秩 |
| `Underdetermined { rank }` | 系数矩阵的秩小于变量数，`rank` 为实际秩 |
| `ResultNotInDomain` | 解不在目标数域中（如要求整数解但解为分数） |
| `Matrix(MatrixError)` | 内部矩阵运算失败（由 `From<MatrixError>` 自动转换） |
| `Other(String)` | 其他带描述的错误 |

**示例**：

```rust
use ocas_calc::solve::{solve_linear_rational, SolveError};

// 矛盾系统：x + y = 1, x + y = 2
let a = vec![vec![1, 1], vec![1, 1]];
let b = vec![1, 2];
let result = solve_linear_rational(&a, &b);
assert_eq!(result, Err(SolveError::Inconsistent));
```

---

## 辅助函数

### substitute_solution_collected

**签名**：
```rust
pub fn substitute_solution_collected<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a>
```

**功能**：将候选解代入 ODE 方程，合并同类项后返回残差。残差为零说明候选解满足 ODE。主要用于验证。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | 表达式 arena |
| `equation` | `Atom<'a>` | ODE 方程（`lhs - rhs` 形式） |
| `func` | `Atom<'a>` | 未知函数（如 `y(x)`） |
| `sol` | `Atom<'a>` | 候选解表达式 |
| `var` | `Symbol` | 自变量符号 |

**返回值**：`Atom<'a>` — 代入并化简后的残差表达式。为零则候选解正确。

---

## 来源模块

| 模块 | 路径 | 内容 |
|---|---|---|
| 线性求解 | `ocas-calc/src/solve.rs` | `solve_linear_rational`、`solve_linear_integer`、`solve_diophantine`、`SolveError` |
| 多项式系统 | `ocas-poly/src/ideal.rs` | `solve_polynomial_system`、`PolynomialSystemSolution`、`RealSolution`、`ZeroDimSolutions` |
| ODE 分类 | `ocas-calc/src/ode/classify.rs` | `classify_ode`、`ODEType` |
| ODE 一阶 | `ocas-calc/src/ode/first_order.rs` | 可分离、线性、Bernoulli、恰当、齐次求解器 |
| ODE 二阶 | `ocas-calc/src/ode/second_order.rs` | 常系数、Cauchy-Euler、降阶、参数变易求解器 |
| ODE 级数 | `ocas-calc/src/ode/series.rs` | 幂级数、Frobenius 求解器 |
| ODE Laplace | `ocas-calc/src/ode/laplace.rs` | `dsolve_ivp` 的 Laplace 变换实现 |
| ODE 系统 | `ocas-calc/src/ode/systems.rs` | `dsolve_system` 的特征值分解实现 |
