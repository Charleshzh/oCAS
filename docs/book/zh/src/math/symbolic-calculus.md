# 进阶：符号微积分

## 前提知识

- 微积分基础：导数定义、基本求导公式、Taylor 级数概念
- 表达式树与递归数据结构（见 [多项式代数](./polynomial-algebra.md)）
- oCAS 表达式系统（`Atom`、`AtomNode`、`AtomArena`）

建议先阅读：[多项式代数](./polynomial-algebra.md)、[线性代数](./linear-algebra.md)。

---

## 基础概念

### 表达式树上的符号微积分

在符号计算中，函数 $f(x)$ 被表示为一棵**表达式树**（expression tree）。树的叶节点是常数或变量，内部节点是运算符（`Add`、`Mul`、`Pow`）或函数调用（`sin`、`exp`、`log` 等）。

例如，$f(x) = x^2 + \sin(x)$ 的表达式树为：

```
Add
├── Pow
│   ├── Var("x")
│   └── Num(2)
└── Fun("sin", [Var("x")])
```

符号微分的核心思想是**模式匹配**（pattern matching）：对表达式树的每个节点，根据其类型应用对应的微分规则，递归地构造导数表达式树。

### 基本求导规则

给定可微函数 $f(x)$ 和 $g(x)$：

| 规则 | 公式 |
|---|---|
| 常数 | $\frac{d}{dx}[c] = 0$ |
| 恒等 | $\frac{d}{dx}[x] = 1$ |
| 和 | $\frac{d}{dx}[f + g] = f' + g'$ |
| 积 | $\frac{d}{dx}[f \cdot g] = f' \cdot g + f \cdot g'$ |
| 链式 | $\frac{d}{dx}[f(g(x))] = f'(g(x)) \cdot g'(x)$ |
| 幂（指数为常数 $n$） | $\frac{d}{dx}[f^n] = n \cdot f^{n-1} \cdot f'$ |
| 指数（底为常数 $a$） | $\frac{d}{dx}[a^g] = a^g \cdot \ln a \cdot g'$ |
| 一般幂 | $\frac{d}{dx}[f^g] = f^g \cdot (\ln f \cdot g' + g \cdot f'/f)$ |

### 初等函数导数表

| $f(u)$ | $f'(u)$ |
|---|---|
| $\sin u$ | $\cos u$ |
| $\cos u$ | $-\sin u$ |
| $\exp u$ | $\exp u$ |
| $\ln u$ | $u^{-1}$ |
| $\sqrt{u}$ | $(2\sqrt{u})^{-1}$ |
| $\tan u$ | $\sec^2 u$ |
| $\sec u$ | $\sec u \cdot \tan u$ |
| $\arctan u$ | $(1 + u^2)^{-1}$ |

注意：上表给出的是 $f'(u)$（关于 $u$ 的导数），实际 $\frac{d}{dx}[f(u(x))] = f'(u) \cdot u'(x)$，即还需乘以链式法则因子 $u'(x)$。

### Taylor 级数

设 $f$ 在点 $a$ 的某邻域内无穷可微，则 $f$ 在 $a$ 处的 **Taylor 级数**为：

$$
f(x) = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x - a)^n
$$

其中 $f^{(n)}(a)$ 是 $f$ 在 $a$ 处的 $n$ 阶导数。截断到 $N$ 阶得到 **Taylor 多项式**：

$$
T_N(x) = \sum_{n=0}^{N} \frac{f^{(n)}(a)}{n!} (x - a)^n
$$

常见展开（在 $a = 0$ 处，即 Maclaurin 级数）：

$$
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots
$$

$$
\sin x = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \cdots
$$

$$
\cos x = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \cdots
$$

$$
\ln(1+x) = x - \frac{x^2}{2} + \frac{x^3}{3} - \cdots \quad (|x| < 1)
$$

---

## 核心理论

### 表达式树上的递归微分

符号微分的核心算法是对表达式树进行**结构递归**（structural recursion）。对每个节点类型，应用对应的微分规则，递归地对子表达式求导。

**算法** `diff(expr, var)`：

1. **常数** `Num(c)`：返回 `Num(0)`。
2. **变量** `Var(v)`：若 $v = \text{var}$，返回 `Num(1)`；否则返回 `Num(0)`。
3. **和** `Add(a₁, a₂, …, aₙ)`：逐项求导，返回 `Add(diff(a₁), diff(a₂), …, diff(aₙ))`。
4. **积** `Mul(a₁, a₂, …, aₙ)`：应用广义乘积法则。对 $n$ 个因子的乘积 $a_1 a_2 \cdots a_n$：

$$
\frac{d}{dx}\left[\prod_{i=1}^n a_i\right] = \sum_{i=1}^n \left(\frac{d a_i}{dx} \cdot \prod_{j \neq i} a_j\right)
$$

即对每个位置 $i$，将第 $i$ 个因子替换为其导数，其余因子保持不变，然后将 $n$ 个结果相加。

5. **幂** `Pow(base, exp)`：分三种情况——
   - 指数为常数 $n$：$\frac{d}{dx}[b^n] = n \cdot b^{n-1} \cdot b'$
   - 底为常数 $a$：$\frac{d}{dx}[a^u] = a^u \cdot \ln a \cdot u'$
   - 一般情形（广义幂法则）：$\frac{d}{dx}[b^e] = b^e \cdot (\ln b \cdot e' + e \cdot b'/b)$

6. **函数调用** `Fun(name, [u])`：查表得到 $f'(u)$，再乘以 $u'$（链式法则）。
7. **未知函数**：返回未求值形式 `Derivative(f(x), x)`。

**关键性质**：由于 oCAS 的表达式使用 hash-consing（结构共享），递归微分产生的中间表达式自动共享相同子树。微分完成后，通过重写引擎化简（例如 $0 \cdot f \to 0$，$1 \cdot f \to f$），消除冗余项。

### Taylor 展开的递推算法

Taylor 展开的直接实现可以利用**逐次微分**：

**算法** `taylor(expr, var, point, order)`：

1. 令 `current = expr`，`sum = 0`。
2. 对 $n = 0, 1, \dots, \text{order}$：
   - 计算 $f^{(n)}(\text{point})$：将 `current` 中的 `var` 替换为 `point`（即求值）。
   - 计算系数 $c_n = f^{(n)}(\text{point}) / n!$。
   - 若 $n = 0$，项为 $c_0$；否则项为 $c_n \cdot (x - \text{point})^n$。
   - 累加：$\text{sum} \mathrel{+}= \text{term}$。
   - 若 $n < \text{order}$，更新 `current = diff(current, var)`。
3. 化简并返回 `sum`。

**阶乘的处理**：$1/n!$ 以幂的形式 $n!^{-1}$ 表示，避免引入浮点数。`mul_by_factorial_inverse` 函数计算 $n!$（使用 `i64` 乘法，对 $n \leq 20$ 安全），然后构造 $c \cdot (n!)^{-1}$。

**复杂度**：循环 $n = 0, \dots, N$ 中每一步执行一次替换（求值），但仅在 $n < N$ 时执行一次符号微分。符号微分本身是 $O(|\text{expr}|)$ 的（$|\text{expr}|$ 为表达式树大小），但化简可能增加节点数。总体而言，$N$ 阶 Taylor 展开需要 $N$ 次符号微分和 $N+1$ 次替换。

### 替换的语义

`substitute(expr, var, replacement)` 将表达式中所有自由出现的变量 `var` 替换为 `replacement`。

**实现策略**：oCAS 使用**自底向上变换**（bottom-up transform）——通过 `transform` 函数遍历表达式树，对每个叶节点检查是否为要替换的变量。若是，返回 `replacement`；否则保持不变。

**深拷贝 vs 引用**：

- oCAS 的 `Atom` 是 arena 上的 **Copy 句柄**（指针），而非拥有所有权的值。`transform` 在 arena 上分配新节点来构造结果，但未修改的子树直接复用原句柄（零拷贝）。
- 这意味着 `substitute` 的语义是**纯函数式**的：原表达式不变，结果是新分配的表达式树，其中被替换的路径指向新节点，未受影响的子树与原树共享。
- 与"破坏性修改"不同，这保证了引用透明性：对同一表达式多次替换不会互相干扰。

**替换的局限**：`substitute` 是**语法替换**而非语义替换。它替换的是字面上的变量名，而非"值"。例如，对 $\int_0^x f(t)\, dt$ 中的 $t$ 做替换 $t \to y$ 会错误地改变积分变量——这类语义替换需要更精细的绑定分析。

---

## 在 oCAS 中的实现

### `diff`：符号微分

**文件**：`ocas-calc/src/derivative.rs`

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

**实现结构**：

`diff` 是入口函数，它调用内部的 `diff_raw` 执行实际递归，然后用 `simplify` 化简结果，最后用 `normalize` 归一化。

```rust
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = diff_raw(ctx, expr, var);
    let simplified = simplify(ctx, raw, &rules, 20);
    normalize(ctx, simplified)
}
```

`diff_raw` 对 `AtomNode` 的每个变体进行模式匹配：

- **`Num(_)`**：常数，导数为 0。
- **`Var(v)`**：若 `v == var` 则导数为 1，否则为 0。
- **`Add(args)`**：逐项调用 `diff_raw`，结果收集到新的 `Add` 节点。
- **`Mul(args)`**：广义乘积法则——对 $n$ 个因子，生成 $n$ 个乘积项（每项将一个因子替换为其导数），然后相加。
- **`Pow(base, exp)`**：三种情况（指数常数 / 底常数 / 一般），见上文。
- **`Fun(name, args)`**：委托给 `diff_function`。

`diff_function` 维护一张**硬编码导数表**，覆盖 `sin`、`cos`、`exp`、`log`、`sqrt`、`tan`、`sec`、`atan` 八个函数。对每个函数 $f(u)$，表中存储 $f'(u)$（关于 $u$ 的导数），然后乘以链式法则因子 $\frac{du}{dx}$。

对于表中没有的函数，返回**未求值形式**：

```rust
ctx.fun("Derivative", &[ctx.fun(name_str, args), ctx.var(var.as_str())])
```

即 `Derivative(f(x), x)`，表示无法自动求导的导数。

**化简管线**：`diff` 产生的原始导数表达式可能包含冗余（如 `0 * sin(x)`、`1 * cos(x)`）。`simplify` 使用 `calculus_rules`（微积分专用规则 + 默认化简规则）进行最多 20 轮不动点化简，然后 `normalize` 将结果归一化为规范形式。

### `taylor`：Taylor 展开

**文件**：`ocas-calc/src/series.rs`

```rust
pub fn taylor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    point: Atom<'a>,
    order: usize,
) -> Atom<'a>
```

**实现细节**：

1. 预计算 `(x - point)` 作为公共子表达式 `x_minus_p`。
2. 循环 $n = 0, \dots, \text{order}$：
   - `substitute(ctx, current, var, point)` 求 $f^{(n)}(a)$。
   - `mul_by_factorial_inverse` 计算 $f^{(n)}(a) / n!$。
   - 构造项 $c_n \cdot (x - a)^n$（$n = 0$ 时省略幂因子）。
   - 累加到 `sum`。
   - 若 $n < \text{order}$，`current = diff(ctx, current, var)` 计算下一阶导数。
3. 化简 + 归一化。

**`mul_by_factorial_inverse`**：

```rust
fn mul_by_factorial_inverse<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, n: usize) -> Atom<'a> {
    if n == 0 { return expr; }
    let mut fact: i64 = 1;
    for i in 2..=n {
        fact = fact.checked_mul(i as i64).expect("factorial fits in i64");
    }
    ctx.mul(&[expr, ctx.pow(ctx.num(fact), ctx.num(-1))])
}
```

将 $1/n!$ 表示为 $(n!)^{-1}$，即 `Pow(Num(n!), Num(-1))`，保持精确有理数表示。`checked_mul` 确保 $n \leq 20$（$20! = 2\,432\,902\,008\,176\,640\,000 < 2^{63}$）。

**示例**：

```rust
// e^x 在 x=0 处展开到 3 阶
let x = ctx.var("x");
let expr = ctx.fun("exp", &[x]);
let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 3);
// result = "1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3))"
// 即 1 + x + x²/2 + x³/6

// sin(x) 在 x=0 处展开到 5 阶
let sin_x = ctx.fun("sin", &[x]);
let result = taylor(&ctx, sin_x, Symbol::new("x"), ctx.num(0), 5);
// result = "x + (-1*(6^-1)*(x^3)) + ((120^-1)*(x^5))"
// 即 x - x³/6 + x⁵/120（偶数阶系数为零，自动消失）
```

### `substitute`：变量替换

**文件**：`ocas-calc/src/series.rs`

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a>
```

**实现**：

```rust
pub fn substitute<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    transform(ctx, expr, |a| match a.node() {
        AtomNode::Var(v) if *v == var => Some(replacement),
        _ => None,
    })
}
```

使用 `transform`（自底向上变换）遍历表达式树。闭包对每个节点检查：若为 `Var(v)` 且 `v == var`，返回 `Some(replacement)`（替换）；否则返回 `None`（保持不变）。

**关键性质**：

- **纯函数式**：不修改原表达式，返回新树。未修改的子树通过 `Atom` 的 Copy 语义共享。
- **深替换**：替换所有自由出现，包括嵌套在函数参数中的（如 `sin(x)` 中的 `x`）。
- **无捕获检查**：不检查替换后是否产生变量捕获——这是语法层面的替换。

### `integrate`：符号积分

**文件**：`ocas-calc/src/integral/mod.rs`

```rust
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a>
```

oCAS 的积分器采用**分层管线**架构（详见 [符号积分算法](../algorithms/integration.md) 和 [Risch 积分算法](./risch-algorithm.md)）：

1. **直接查表层**（`integrate_raw`）：常数、变量、和式以及简单乘积/幂/函数的直接积分规则。
2. **有理函数层**（`rational.rs`）：部分分式分解 + Hermite 约化 + 对数部分（对数导数恒等式、配方、Rothstein–Trager）。
3. **Risch 算法层**（`risch.rs`）：微分域塔 $\mathbb{Q}(x, t_1, \dots, t_n)$ 上的结构定理。
4. **三角层**（`trig.rs`）：Euler 公式重写为复指数（$t = e^{ix}$）→ 对重写结果再次尝试 Risch → `realify` 回转实数形式。
5. **特殊函数层**（`special.rs`）：erf、erfi、Ei、Si、Ci、Shi、Chi、Fresnel 等特殊函数表。
6. **启发式层**（`heuristic.rs`）：分部积分（LIATE 排序）、三角替换（$\sqrt{a^2 - x^2} \to a\sin\theta$）、Weierstrass $t = \tan(x/2)$ 有理化、Euler 替换。
7. **回退**：返回未求值形式 `Integral(expr, var)`。

`integrate_with_fuel` 变体接受 `Fuel` 预算，防止化简阶段无限循环。积分遍历本身受 `MAX_DEPTH = 8` 的内部递归深度限制（Risch 层另有 `MAX_RISCH_DEPTH = 16`）：

```rust
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>>
```

---

## 进阶话题

### 微分与化简的交互

`diff` 的正确性不仅取决于微分规则本身，还取决于化简质量。例如：

- $\frac{d}{dx}[x^2 + 3x] = 2x + 3$——需要化简 $1 \cdot 3$ 为 $3$。
- $\frac{d}{dx}[e^{\ln x}] = e^{\ln x} \cdot \frac{1}{x} = 1$——需要化简 $e^{\ln x} \to x$。

oCAS 通过 `calculus_rules` 和默认规则的多轮化简处理这类情况。但某些恒等式（如 $\sin^2 x + \cos^2 x = 1$）可能需要专门的三角化简规则。

### 未求值导数与积分

当 `diff` 遇到未知函数（不在导数表中的函数）时，返回 `Derivative(f, x)`。类似地，`integrate` 对无法处理的积分返回 `Integral(f, x)`。

这些未求值形式是**合法的表达式节点**——它们可以参与后续的符号运算。例如，$\frac{d}{dx}[\text{Derivative}(f, x)]$ 会返回 `Derivative(Derivative(f, x), x)`，表示二阶导数。

### 与数值求值的衔接

符号微分产生的表达式可以传递给 `ExpressionEvaluator`（见 [求值与 JIT](../api/rust-evaluation.md)）进行高效数值求值。典型工作流：

1. 符号微分：`diff(&ctx, f, x)` → 导数表达式 $f'(x)$。
2. 编译求值器：`ExpressionEvaluator::compile(f')`。
3. 数值求值：`evaluator.evaluate(&[x_val])` → $f'(x_{\text{val}})$。

这比数值差分（有限差分）精确得多，且避免了步长选择问题。

---

## 参考文献

1. **Geddes, K. O., Czapor, S. R. & Labahn, G.** *Algorithms for Computer Algebra.* Kluwer Academic Publishers, 1992. — 第 11 章"Symbolic Differentiation"涵盖表达式树上的递归微分算法、化简策略和导数表设计。
2. **Bronstein, M.** *Symbolic Integration I: Transcendental Functions.* Springer, 2005. — Risch 算法的权威参考，oCAS 积分器的设计依据。
3. **Cohen, H.** *A Course in Computational Algebraic Number Theory.* Springer, 1993. — 第 1 章包含算法复杂度分析的一般方法论。
