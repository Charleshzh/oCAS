# 微积分

> 来源：`ocas-calc/src/`

oCAS 的微积分模块提供符号求导、积分、级数展开和部分分式分解。所有运算在 `AtomArena` 上执行，结果为新的 `Atom` 句柄。无法化简的导数或积分以保留函数形式 `Derivative(expr, var)` 或 `Integral(expr, var)` 返回。

**模块结构**：

| 子模块 | 功能 |
|---|---|
| `derivative` | 符号微分 `diff` |
| `integral` | 分层积分管线 `integrate`、`integrate_heuristic`、`integrate_with_fuel` |
| `series` | Taylor 展开 `taylor`、替换 `substitute` |
| `partial_fraction` | 部分分式分解 `apart`、合并 `together` |
| `tower` | Risch 算法的微分域塔构造 |
| `ode` | 常微分方程求解（见 [求解器](./rust-solvers.md)） |

---

## `diff`

**签名**：

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**功能**：对 `expr` 关于 `var` 求符号导数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena，用于构造结果 |
| `expr` | `Atom<'a>` | 被求导的表达式 |
| `var` | `Symbol` | 求导变量（interned 字符串） |

**返回值**：`Atom<'a>` — 求导并化简后的表达式。

**算法说明**：

`diff` 递归遍历表达式树，对每种节点类型应用标准求导规则：

| 节点类型 | 求导规则 |
|---|---|
| `Num(_)` | $\frac{d}{dx} c = 0$ |
| `Var(s)` | 若 $s = \text{var}$ 则为 $1$，否则为 $0$ |
| `Add([a₁, …, aₙ])` | 和的导数 = 导数的和：$\sum \frac{d a_i}{dx}$ |
| `Mul([a₁, …, aₙ])` | 乘积法则：$\frac{d}{dx}\prod a_i = \sum_i \left(\frac{d a_i}{dx} \prod_{j \neq i} a_j\right)$ |
| `Pow([base, exp])` | 幂函数法则（含链式法则），支持常数和非常数指数 |
| `Fun(name, args)` | 内置函数查表 + 链式法则 |

**内置函数求导表**：

| 函数 $f(u)$ | 导数 $f'(u) \cdot u'$ |
|---|---|
| $\sin(u)$ | $\cos(u) \cdot u'$ |
| $\cos(u)$ | $-\sin(u) \cdot u'$ |
| $\tan(u)$ | $\sec^2(u) \cdot u'$ |
| $\sec(u)$ | $\sec(u)\tan(u) \cdot u'$ |
| $\exp(u)$ | $\exp(u) \cdot u'$ |
| $\log(u)$ | $u^{-1} \cdot u'$ |
| $\sqrt{u}$ | $(2\sqrt{u})^{-1} \cdot u'$ |
| $\text{atan}(u)$ | $(1 + u^2)^{-1} \cdot u'$ |

表在 `derivative.rs` 的 `diff_function` 中硬编码。**不在**表中的函数（如 $\text{asin}$、$\sinh$、$\cosh$、$\tanh$ 等）也返回未求值形式 `Derivative(f, var)`。

**未求值形式**：当求导函数不在内置表中时，结果为 `Derivative(expr, var)`，可通过后续模式匹配检测（匹配函数名 `"Derivative"` 的 `Fun` 节点；注意 `is_fallback` 只检测 `"Integral"`，不检测 `Derivative`）。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::diff;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// d/dx [sin(x)] = cos(x)
let sin_x = ctx.fun("sin", &[x]);
let result = diff(&ctx, sin_x, Symbol::new("x"));
assert_eq!(result.to_string(), "cos(x)");

// d/dx [x^3] = 3*x^2
let x_cubed = ctx.pow(x, ctx.num(3));
let result = diff(&ctx, x_cubed, Symbol::new("x"));
assert_eq!(result.to_string(), "3*(x^2)");

// d/dx [exp(x^2)] = 2*x*exp(x^2) （链式法则）
let inner = ctx.pow(x, ctx.num(2));
let expr = ctx.fun("exp", &[inner]);
let result = diff(&ctx, expr, Symbol::new("x"));
// 输出包含 2*x*exp(x^2) 的化简形式

// 未知函数返回 Derivative 形式
let f = ctx.fun("my_func", &[x]);
let result = diff(&ctx, f, Symbol::new("x"));
// result 包含 Derivative(my_func(x), x)
```

**参见**：[taylor](#taylor)、[integrate](#integrate)

---

## `integrate`

**签名**：

```rust
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**功能**：对 `expr` 关于 `var` 求符号积分。使用分层积分管线，逐层尝试不同算法。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 被积表达式 |
| `var` | `Symbol` | 积分变量 |

**返回值**：`Atom<'a>` — 积分结果。无法求解时返回 `Integral(expr, var)`。

**积分管线**：

`integrate` 内部调用 `integrate_raw`（递归，深度上限 `MAX_DEPTH = 8`）。管线分两个阶段：

```
┌───────────────────────────────────────────┐
│ 阶段一：integrate_raw 按节点类型分派        │
├───────────────────────────────────────────┤
│ Num / Var          → 常数积分：c·x；       │
│                      ∫x dx = (1/2)·x²      │
│ Add                → 逐项递归积分           │
│ Mul (integrate_product)  → 提取常数因子    │
│ Pow (integrate_power)    → 幂函数积分      │
│                      ∫xⁿ, ∫(a·x+b)ⁿ,       │
│                      分数指数 (a·x+b)^(p/q) │
│ Fun (integrate_function) → 查表 + 线性替换  │
│                      ∫sin(u), ∫cos(u),     │
│                      ∫exp(u), ∫log(u)      │
│                      若 u = a·x+b 自动应用  │
│                      链式法则              │
├───────────────────────────────────────────┤
│ 阶段二：以上分派失败（返回 fallback 形式）  │
│ → try_risch_or_fallback，按顺序尝试：      │
├───────────────────────────────────────────┤
│ 1. 有理函数积分   (integrate_rational)     │
│ 2. Risch 算法     (risch_integrate)        │
│ 3. 三角→指数重写  (trig_to_exp +           │
│    + Risch + realify)                     │
│ 4. 特殊函数       (special_integrate)      │
│ 5. 启发式方法     (heuristic_integrate)    │
│    分部积分 (LIATE)、三角替换、            │
│    Weierstrass、Euler 替换                 │
│ 6. 未求值形式     (fallback)               │
│    Integral(expr, var)                    │
└───────────────────────────────────────────┘
```

**阶段一 — 节点类型分派**：

- `Num(_)`：常数的积分为线性函数 `c * x`
- `Var(var)`：$\int x\,dx = \frac{1}{2} x^2$；其他变量视为常数
- `Add`：逐项递归积分再求和
- `Mul`（`integrate_product`）：识别 `c · f(x)` 形式（常数乘以非常数因子），提取常数因子后递归积分
- `Pow`（`integrate_power`）：幂函数积分，支持 $\int x^n\,dx$、线性形式 $(a \cdot x + b)^n$（含 $n = -1$ 时 $\int \frac{dx}{ax+b} = \frac{\log(ax+b)}{a}$）以及分数指数 $(a \cdot x + b)^{p/q}$
- `Fun`（`integrate_function`）：查表 + 线性替换

**查表与线性替换**：

内置积分表覆盖幂函数（`integrate_power`）以及 `sin`、`cos`、`exp`、`log`。当被积函数（或幂的底）的参数为线性形式 $a \cdot x + b$ 时，自动应用链式法则 $\int f(ax+b)\,dx = \frac{1}{a} F(ax+b)$。

**阶段二 — 回退链**（`try_risch_or_fallback`）：

1. **有理函数积分**（`integrate_rational`）：当表达式为 $\var$ 的有理函数（无其他变量或函数应用）时，执行：
   1. 多项式部分逐项积分
   2. Hermite 约化分离有理部分
   3. 对数部分：完全平方（二次分母 → $\log$ 或 $\text{atan}$）+ Rothstein–Trager 结式

   详见 [有理函数积分](#有理函数积分)。

2. **Risch 算法**（`risch_integrate`）：构造微分域塔 $\mathbb{Q}(x, t_1, \dots, t_n)$，逐层递归积分。详见 [Risch 算法](#risch-算法)。

3. **三角→指数重写**：若 Risch 对原表达式失败且表达式含三角函数，通过 `trig_to_exp` 重写为复指数形式，再次尝试 Risch，最后用 `realify` 转回实数形式。

4. **特殊函数**（`special_integrate`）：当 Risch 与三角重写均未成功时，尝试匹配已知的特殊函数反导数。详见 [特殊函数积分](#特殊函数积分)。

5. **启发式方法**（`heuristic_integrate`）：分部积分 (LIATE)、三角替换、Weierstrass 替换、Euler 替换。详见 [integrate_heuristic](#integrate_heuristic)。

6. **未求值形式**（`fallback`）：以上均失败时，返回 `Integral(expr, var)`。可通过 `is_fallback` 检测：

```rust
pub(crate) fn is_fallback<'a>(atom: &Atom<'a>) -> bool {
    matches!(atom.node(), AtomNode::Fun(name, _) if name.as_str() == "Integral")
}
```

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integrate;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ x² dx = (1/3) x³
let expr = ctx.pow(x, ctx.num(2));
let result = integrate(&ctx, expr, Symbol::new("x"));
assert_eq!(result.to_string(), "(3^-1)*(x^3)");

// ∫ sin(x) dx = -cos(x)
let expr = ctx.fun("sin", &[x]);
let result = integrate(&ctx, expr, Symbol::new("x"));
// 输出：-cos(x) 的化简形式

// ∫ 1/x dx = log(x) （通过有理函数积分器或 Risch）
let expr = ctx.pow(x, ctx.num(-1));
let result = integrate(&ctx, expr, Symbol::new("x"));
// 输出：log(x)

// 不可积表达式返回 Integral 形式
let expr = ctx.fun("my_func", &[x]);
let result = integrate(&ctx, expr, Symbol::new("x"));
// result 为 Integral(my_func(x), x)
```

**参见**：[integrate_heuristic](#integrate_heuristic)、[integrate_with_fuel](#integrate_with_fuel)、[有理函数积分](#有理函数积分)、[Risch 算法](#risch-算法)

---

## `integrate_heuristic`

**签名**：

```rust
pub fn integrate_heuristic<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**功能**：尝试启发式积分技巧。仅使用分部积分、三角替换、Weierstrass 替换和 Euler 替换，不调用 Risch 算法或有理函数积分器。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 被积表达式 |
| `var` | `Symbol` | 积分变量 |

**返回值**：`Atom<'a>` — 若启发式成功则返回积分结果，否则返回 `Integral(expr, var)`。

**启发式方法**：

### 1. 分部积分（LIATE 启发式）

对 $\int u \cdot v' \, dx$ 形式的乘积，按 **LIATE** 优先级选择 $u$（优先级从高到低）：

| 优先级 | 类型 | LIATE 评分 | 示例 |
|---|---|---|---|
| 0（最高） | 对数 (L) | 0 | $\log(x)$ |
| 1 | 反三角 (I) | 1 | $\text{asin}(x)$ |
| 2 | 代数 (A) | 2 | $x^2$ |
| 3 | 三角 (T) | 3 | $\sin(x)$ |
| 4 | 指数 (E) | 4 | $\exp(x)$ |
| 5（最低） | 其他 | 5 | — |

选择评分最低的因子作为 $u$，其余作为 $v'$。递归深度上限为 `PARTS_MAX_DEPTH = 2`。

### 2. 三角替换

匹配以下模式并返回已知反导数：

| 被积函数 | 替换 | 结果 |
|---|---|---|
| $\frac{1}{\sqrt{a^2 - x^2}}$ | $x = a\sin\theta$ | $\text{asin}(x/a)$ |
| $\frac{1}{\sqrt{a^2 + x^2}}$ | $x = a\sinh t$ | $\text{asinh}(x/a)$ |
| $\frac{1}{\sqrt{x^2 - a^2}}$ | $x = a\cosh t$ | $\text{acosh}(x/a)$ |
| $\sqrt{a^2 - x^2}$ | $x = a\sin\theta$ | $\frac{x\sqrt{a^2-x^2} + a^2\,\text{asin}(x/a)}{2}$ |
| $\sqrt{a^2 + x^2}$ | $x = a\sinh t$ | $\frac{x\sqrt{a^2+x^2} + a^2\,\text{asinh}(x/a)}{2}$ |
| $\sqrt{x^2 - a^2}$ | $x = a\cosh t$ | $\frac{x\sqrt{x^2-a^2} - a^2\,\text{acosh}(x/a)}{2}$ |

### 3. Weierstrass 替换

当被积函数是 $\sin(u)$ 和 $\cos(u)$ 的有理函数（$u$ 关于 $\text{var}$ 线性）时，应用 $t = \tan(u/2)$ 替换：

$$\sin(u) = \frac{2t}{1+t^2}, \quad \cos(u) = \frac{1-t^2}{1+t^2}$$

将三角有理函数化为 $t$ 的有理函数后递归积分。

### 4. Euler 替换

当被积函数包含 $\sqrt{ax^2 + bx + c}$ 时，尝试 Euler 替换消去根号。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integrate_heuristic;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ x·exp(x) dx — 分部积分 (LIATE: x=代数, exp=指数)
let expr = ctx.mul(&[x, ctx.fun("exp", &[x])]);
let result = integrate_heuristic(&ctx, expr, Symbol::new("x"));
// 结果包含 (x-1)*exp(x) 的化简形式

// ∫ 1/√(1 - x²) dx = asin(x)
// 构造 (1 - x²)^(-1/2)：指数为 -1·(2^-1)
let one = ctx.num(1);
let sqrt_arg = ctx.add(&[one, ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
let neg_half = ctx.mul(&[ctx.num(-1), ctx.pow(ctx.num(2), ctx.num(-1))]);
let expr = ctx.pow(sqrt_arg, neg_half);
let result = integrate_heuristic(&ctx, expr, Symbol::new("x"));
// 输出：asin(x) 的化简形式
```

**参见**：[integrate](#integrate)

---

## `integrate_with_fuel`

**签名**：

```rust
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>, FuelError>
```

**功能**：带燃料预算的积分。积分遍历本身使用内部深度限制（`MAX_DEPTH = 8`），此入口将 `fuel` 线程化到积分后的两个化简阶段，防止病态结果使重写引擎无限循环。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 被积表达式 |
| `var` | `Symbol` | 积分变量 |
| `fuel` | `&Fuel` | 燃料预算（重写步数上限） |

**返回值**：
- `Ok(Atom<'a>)` — 积分化简成功
- `Err(FuelError)` — 燃料在化简过程中耗尽

**与 `integrate` 的区别**：

`integrate` 在化简阶段不设上限；`integrate_with_fuel` 对化简阶段施加燃料约束。积分遍历本身的深度限制（`MAX_DEPTH`）两者相同。

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use ocas_core::fuel::Fuel;
use ocas_calc::integrate_with_fuel;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, x, Symbol::new("x"), &fuel);
match result {
    Ok(expr) => assert_eq!(expr.to_string(), "(2^-1)*(x^2)"),
    Err(_) => panic!("fuel exhausted"),
}
```

**参见**：[integrate](#integrate)

---

## `taylor`

**签名**：

```rust
pub fn taylor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    point: Atom<'a>,
    order: usize,
) -> Atom<'a>
```

**功能**：计算 `expr` 在 `point` 处关于 `var` 的 Taylor 展开，展开到 `order` 阶（含）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 被展开的表达式 |
| `var` | `Symbol` | 展开变量 |
| `point` | `Atom<'a>` | 展开点（可为任意表达式，通常为 `0` 或某常数） |
| `order` | `usize` | 展开阶数（含） |

**返回值**：`Atom<'a>` — 截断多项式：

$$\sum_{n=0}^{\text{order}} \frac{f^{(n)}(\text{point})}{n!} \cdot (\text{var} - \text{point})^n$$

**算法说明**：

通过重复符号微分 + 在展开点处求值计算各阶系数：

1. 对 $n = 0, 1, \dots, \text{order}$：
   - 计算 $f^{(n)}(x)$（通过 `diff`）
   - 用 `substitute` 在 $x = \text{point}$ 处求值得 $f^{(n)}(\text{point})$
   - 乘以 $\frac{1}{n!}$（通过 `mul_by_factorial_inverse`）
   - 乘以 $(x - \text{point})^n$
2. 将所有项求和，应用微积分化简规则

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::taylor;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// exp(x) 在 x=0 处展开到 3 阶
let expr = ctx.fun("exp", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
assert_eq!(
    result.to_string(),
    "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))"
);

// sin(x) 在 x=0 处展开到 5 阶
let expr = ctx.fun("sin", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 5);
// 输出：x + (-1*(6^-1)*(x^3)) + ((120^-1)*(x^5))
```

**参见**：[diff](#diff)、[substitute](#substitute)

---

## `substitute`

**签名**：

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a>
```

**功能**：将 `expr` 中所有出现的变量 `var` 替换为 `replacement`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ctx` | `&AtomArena<'a>` | 表达式 arena |
| `expr` | `Atom<'a>` | 被替换的表达式 |
| `var` | `Symbol` | 要替换的变量名 |
| `replacement` | `Atom<'a>` | 替换表达式 |

**返回值**：`Atom<'a>` — 替换后的新表达式。原表达式不受影响（Arena 上的不可变数据结构）。

**语义**：

- 深度遍历表达式树
- 仅替换精确匹配 `var` 的叶节点 `Var(Symbol)`
- 替换后不自动化简（如需化简，需额外调用 `simplify`）

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::substitute;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");
let y = ctx.var("y");

// (x² + x)[x → y] = y² + y
let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), x]);
let result = substitute(&ctx, expr, Symbol::new("x"), y);
assert_eq!(result.to_string(), "(y^2) + y");

// sin(x)[x → 2*y] = sin(2*y)
let expr = ctx.fun("sin", &[x]);
let two_y = ctx.mul(&[ctx.num(2), y]);
let result = substitute(&ctx, expr, Symbol::new("x"), two_y);
assert_eq!(result.to_string(), "sin(2*y)");
```

**参见**：[taylor](#taylor)、[diff](#diff)

---

## `apart`

**签名**：

```rust
pub fn apart<D: EuclideanDomain>(
    num: &DenseUnivariatePolynomial<D>,
    den: &DenseUnivariatePolynomial<D>,
) -> (
    Option<DenseUnivariatePolynomial<D>>,
    Vec<PartialFractionTerm<D>>,
)
```

**功能**：对有理函数 $\frac{\text{num}(x)}{\text{den}(x)}$ 执行部分分式分解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `num` | `&DenseUnivariatePolynomial<D>` | 分子多项式 |
| `den` | `&DenseUnivariatePolynomial<D>` | 分母多项式 |

**返回值**：`(Option<poly>, Vec<PartialFractionTerm<D>>)` — 多项式部分（若分子次数 ≥ 分母次数）和部分分式项列表。

**数学原理**：

给定适当分式 $\frac{p(x)}{q(x)}$（即 $\deg(p) < \deg(q)$），将分母作无平方分解 $q = \prod_i f_i^{e_i}$，然后分解为：

$$\frac{p(x)}{q(x)} = \text{poly\_part} + \sum_i \sum_{k=1}^{e_i} \frac{a_{i,k}(x)}{f_i(x)^k}$$

若 $\deg(p) \geq \deg(q)$，先执行多项式除法，商作为 `poly_part`，余式继续分解。

**`PartialFractionTerm<D>` 结构**：

```rust
pub struct PartialFractionTerm<D: EuclideanDomain> {
    pub numer: DenseUnivariatePolynomial<D>,  // 分子 a_{i,k}(x)
    pub denom: DenseUnivariatePolynomial<D>,  // 不可约因子 f_i(x)
    pub exp: usize,                            // 指数 k
}
```

表示分数 $\frac{\text{numer}(x)}{\text{denom}(x)^{\text{exp}}}$。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;

// 1 / (x² - 1) — 分母无平方（单一平方自由因子 x²-1）
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1),  // 常数项
    Rational::new(0, 1),   // x 系数
    Rational::new(1, 1),   // x² 系数
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_none()); // deg(num) < deg(den)，无多项式部分
// x²-1 无平方 → 返回单项：分母为 x²-1、指数为 1
// （apart 基于无平方分解而非不可约分解，不会拆成 (x-1)(x+1)）

// (x³ + 1) / (x² - 1) — 需要多项式除法
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(0, 1),
    Rational::new(0, 1), Rational::new(1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_some()); // 商多项式 x
```

**参见**：[together](#together)、[有理函数积分](#有理函数积分)

---

## `together`

**签名**：

```rust
pub fn together<D: EuclideanDomain>(
    poly_part: Option<&DenseUnivariatePolynomial<D>>,
    terms: &[PartialFractionTerm<D>],
) -> (DenseUnivariatePolynomial<D>, DenseUnivariatePolynomial<D>)
```

**功能**：将部分分式项合并回单一有理函数。是 `apart` 的逆操作。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `poly_part` | `Option<&DenseUnivariatePolynomial<D>>` | 多项式部分（`apart` 返回的第一项） |
| `terms` | `&[PartialFractionTerm<D>]` | 部分分式项列表 |

**返回值**：`(numerator, denominator)` — 合并后的分子和分母多项式。

**算法**：

将多项式部分和所有分式项通分：计算公分母 $\text{lcm}$，将每项扩展到公分母后求和，返回结果的分子和分母。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::{apart, together};

let d = RationalDomain;
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1),
]);

// apart → round-trip → together 应还原原始分数
let (poly_part, terms) = apart(&num, &den);
let (result_num, result_den) = together(poly_part.as_ref(), &terms);
// result_num/result_den 等价于 num/den（可能差一个常数因子）
```

**参见**：[apart](#apart)

---

## 有理函数积分

**来源**：`ocas-calc/src/integral/rational.rs`

有理函数积分器 `integrate_rational` 处理 $\var$ 的有理函数（不含其他变量或函数应用）。采用 Bronstein《Symbolic Integration I》第 2 章的标准三步法：

### 1. 多项式部分

对 $\int \sum c_k x^k \, dx = \sum \frac{c_k}{k+1} x^{k+1}$ 逐项积分。

### 2. Hermite 约化

对真分式 $\frac{a(x)}{d(x)}$（$d$ 为首一多项式），分解为：

$$\frac{a}{d} = \frac{d}{dx}\left(\frac{g_{\text{num}}}{g_{\text{den}}}\right) + \frac{a_1}{d_1}$$

其中 $d_1$ 是无平方的。有理部分 $\frac{g_{\text{num}}}{g_{\text{den}}}$ 的导数可直接计算。

### 3. 对数部分

对无平方分母的剩余部分 $\frac{a_1}{d_1}$，根据分母次数选择策略：

| 分母次数 | 策略 |
|---|---|
| 0 | 常数积分 |
| 1 | $c \cdot \log(ax + b)$ |
| 2 | 完全平方 → $\log$ 或 $\text{atan}$ |
| ≥ 3 | Rothstein–Trager 结式法 |

**Rothstein–Trager 方法**：

计算 $R(t) = \text{Res}_x(d, a - t \cdot d')$，$R(t)$ 的有理根 $c_i$ 给出：

$$\int \frac{a}{d} = \sum_i c_i \cdot \log(\gcd(d, a - c_i \cdot d'))$$

当 $R(t)$ 在 $\mathbb{Q}$ 上不能完全分解时，对应项以未求值形式 `Integral(term, var)` 返回。

**返回值**：
- `Some(Atom)` — 表达式是 $\var$ 的有理函数，返回积分结果
- `None` — 表达式不是 $\var$ 的纯有理函数（含其他变量或函数）

**示例**：

```rust
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integral::rational::integrate_rational;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let x = ctx.var("x");

// ∫ 1/x dx = log(x)
let result = integrate_rational(&ctx, ctx.pow(x, ctx.num(-1)), Symbol::new("x"));
assert_eq!(result.unwrap().to_string(), "log(x)");

// ∫ 1/(x²+1) dx = atan(x) （通过完全平方）
// 表达式需构建为 x 的有理函数
```

**参见**：[apart](#apart)、[integrate](#integrate)

---

## Risch 算法

**来源**：`ocas-calc/src/integral/risch.rs`、`ocas-calc/src/tower/`

Risch 算法处理初等超越函数的积分，基于 Bronstein《Symbolic Integration I》第 5–6 章。

### 微分域塔

算法首先构造微分域塔 $\mathbb{Q}(x, t_1, \dots, t_n)$，其中每个生成元 $t_i$ 是对数或指数：

```rust
pub(crate) enum GenKind {
    Constant,  // 常数符号（如虚数单位 I），D t = 0
    Log,       // t_i = log(u)
    Exp,       // t_i = exp(u)
}
```

塔的构造（`build_tower`）遍历表达式中的函数应用，为每个 `log` 或 `exp` 创建新层。限制条件：

- 仅接受 `log` 和 `exp` 函数（三角函数需先通过 `trig_to_exp` 重写）
- 拒绝代数函数（非整数指数，如 $\sqrt{x}$）
- 拒绝代数相关的生成元（如 $\log(x)$ 和 $\log(2x)$）

### 逐层积分

在塔的每一层 $k(t_i)$，算法执行：

1. **Hermite 约化**（`hermite_tower`）：将 $\frac{a}{d}$ 分解为 $D(g) + \frac{a_1}{d_1}$，其中 $d_1$ 无平方
2. **多项式部分**：
   - **原语层**（$t = \log(u)$）：待定系数法，顶层常数由对数约束确定
   - **超指数层**（$t = \exp(u)$）：每层解 Risch 微分方程 $Dq + fq = g$
3. **对数部分**：匹配对数导数恒等式 $c \cdot \frac{D d}{d} \to c \cdot \log(d)$

### 三角函数处理

`integrate` 先对原表达式直接尝试 Risch；失败且表达式含三角函数时，通过 `trig_to_exp` 将三角函数重写为复指数形式后再尝试：

$$\sin(u) \to \frac{e^{iu} - e^{-iu}}{2i}, \quad \cos(u) \to \frac{e^{iu} + e^{-iu}}{2}$$

积分结果通过 `realify` 尝试转回实数形式：

| 复数模式 | 实数形式 |
|---|---|
| $c \cdot \log(u + iv) + c \cdot \log(u - iv)$ | $c \cdot \log(u^2 + v^2)$ |
| $c \cdot \log(u + iv) - c \cdot \log(u - iv)$ | $2c \cdot \text{atan}(v/u)$ |
| $\exp(Iu) \cdot \exp(Iv)$ | $\exp(I(u+v))$（合并指数使抵消可见，如 $e^{Ix}e^{-Ix} \to e^0$） |

若 `realify` 无法匹配模式，返回复数形式（数学上仍正确，可通过微分验证）。

### 深度限制

`MAX_RISCH_DEPTH = 16`，通过线程本地计数器 `RISCH_DEPTH` 和守卫 `RischDepthGuard`（RAII 自动减计数）实现。

**返回值**：
- `Some(Atom)` — 积分成功
- `None` — 表达式超出已实现片段（调用者回退到其他积分器）

**示例**：

```rust
// Risch 算法通过 integrate 自动调用，通常不直接使用
// ∫ exp(x²) · 2x dx = exp(x²) — Risch 可处理的典型例子
// ∫ sin(x)/x dx — 无初等原函数，Risch 返回 None，特殊函数层处理
```

**参见**：[integrate](#integrate)、[特殊函数积分](#特殊函数积分)

---

## 特殊函数积分

**来源**：`ocas-calc/src/integral/special.rs`

当 Risch 算法证明积分无初等原函数时，特殊函数积分器尝试匹配已知的非初等反导数。

### 支持的特殊函数

| 被积函数 | 反导数 | 函数名 |
|---|---|---|
| $e^{-x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erf}(x)$ | 误差函数 |
| $e^{x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erfi}(x)$ | 虚误差函数 |
| $\frac{e^x}{x}$ | $\text{Ei}(x)$ | 指数积分 |
| $\frac{\sin x}{x}$ | $\text{Si}(x)$ | 正弦积分 |
| $\frac{\cos x}{x}$ | $\text{Ci}(x)$ | 余弦积分 |
| $\frac{\sinh x}{x}$ | $\text{Shi}(x)$ | 双曲正弦积分 |
| $\frac{\cosh x}{x}$ | $\text{Chi}(x)$ | 双曲余弦积分 |
| $\sin(x^2)$ | $\sqrt{\pi/2}\,\text{fresnels}(\sqrt{2/\pi}\,x)$ | Fresnel S 积分（函数名 `fresnels`） |
| $\cos(x^2)$ | $\sqrt{\pi/2}\,\text{fresnelc}(\sqrt{2/\pi}\,x)$ | Fresnel C 积分（函数名 `fresnelc`） |

**匹配策略**：

`special_integrate` 检查被积函数的因子结构：
- `erf_family`：匹配 $e^{c \cdot x^2}$ 形式
- `ei_family`：匹配 $e^{cx} / x$ 形式
- `trig_integral_family`：匹配 $\sin(x)/x$、$\cos(x)/x$、$\sinh(x)/x$、$\cosh(x)/x$（自变量必须恰为 $x$）
- `fresnel_family`：匹配 $\sin(cx^2)$、$\cos(cx^2)$

**返回值**：
- `Some(Atom)` — 匹配成功，返回含特殊函数的反导数
- `None` — 不匹配任何已知模式，调用者以 `Integral(expr, var)` 形式返回

**设计说明**：

特殊函数定义与 SymPy 一致，结果可通过 `sympy.integrate` 交叉验证。

**参见**：[Risch 算法](#risch-算法)、[integrate](#integrate)

---

## 辅助函数

### `is_fallback`

```rust
pub(crate) fn is_fallback<'a>(atom: &Atom<'a>) -> bool
```

**功能**：检查 `atom` 是否为未求值积分形式（函数名为 `"Integral"`）。注意：未求值的导数（`Derivative` 形式）由 `diff` 直接构造，不经由此函数检测。

**参见**：[integrate](#integrate)、[diff](#diff)

### `is_constant`

```rust
pub(crate) fn is_constant<'a>(expr: Atom<'a>, var: Symbol) -> bool
```

**功能**：检查 `expr` 是否不包含变量 `var`（即关于 `var` 为常数）。

### `fraction_exponent`

```rust
fn fraction_exponent<'a>(exp: Atom<'a>) -> Option<(i64, i64)>
```

**功能**：将指数原子解析为分数 $p/q$（小整数）。接受 $p \cdot q^{-1}$ 形式，用于幂函数积分。

### `linear_form`

```rust
fn linear_form<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol)
    -> Option<(Atom<'a>, Atom<'a>)>
```

**功能**：若 `expr` 形如 $a \cdot \text{var} + b$（$a, b$ 关于 `var` 为常数），返回 `(a, b)`；否则返回 `None`。

---

## 内部架构

### 化简管线

- `integrate` 的结果经过两级化简：先用默认重写规则 `default_rules`（$0 + x = x$、$1 \cdot x = x$、$x^1 = x$ 等代数恒等式）化简 20 步，再用 `calculus_rules` 化简 10 步，最后 `normalize` 到规范形式
- `diff` 与 `taylor` 直接使用 `calculus_rules` 化简后 `normalize`

`calculus_rules`（`rules.rs`）以 `default_rules` 为基础，追加微积分专用恒等式：`exp(0) → 1`、`log(1) → 0`、`sin(0) → 0`、`cos(0) → 1`、`tan(0) → 0`、$(-1) \cdot (-1) \to 1$、$1^x \to 1$ 等。

### 深度限制

| 限制 | 值 | 用途 |
|---|---|---|
| `MAX_DEPTH` | 8 | `integrate_raw` 递归深度 |
| `MAX_RISCH_DEPTH` | 16 | Risch 算法递归深度 |
| `PARTS_MAX_DEPTH` | 2 | 分部积分递归深度 |

### 模块依赖关系

```
ocas-calc
├── derivative::diff
├── integral
│   ├── mod.rs          — integrate, integrate_with_fuel, integrate_heuristic
│   ├── heuristic.rs    — 分部积分, 三角/Weierstrass/Euler 替换
│   ├── rational.rs     — Hermite 约化, Rothstein–Trager
│   ├── risch.rs        — Risch 算法核心
│   ├── rde.rs          — Risch 微分方程求解器
│   ├── trig.rs         — 三角→指数重写, realify
│   ├── special.rs      — erf/Ei/Si/Ci/Fresnel
│   └── tower/
│       ├── build.rs    — 微分域塔构造
│       ├── elem.rs     — 塔元素类型 (KElem, KPoly, KRat)
│       └── convert.rs  — Atom ↔ 塔元素转换
├── series.rs           — taylor, substitute
├── partial_fraction.rs — apart, together
└── rules.rs            — 微积分化简规则
```

---

## 参见

- [数学：符号微积分](../math/symbolic-calculus.md) — 微分与 Taylor 展开的数学原理
- [数学：Risch 积分算法](../math/risch-algorithm.md) — 微分域塔与 Risch 算法的完整数学推导
- [数学：多项式 GCD 与因式分解](../math/poly-gcd-factoring.md) — Hermite 约化与无平方分解的数学基础
- [求解器 API](./rust-solvers.md) — ODE 求解（`dsolve`、`classify_ode`）
- [表达式系统](./rust-expressions.md) — `Atom`、`AtomArena`、`Symbol` 的基本用法
