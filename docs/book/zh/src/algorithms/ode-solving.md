# ODE 求解算法

oCAS 的 ODE 模块（`ocas-calc/src/ode/`）实现了常微分方程的符号求解引擎。本章详述
每种 ODE 类型的分类条件、求解算法的实现细节和设计选择。数学理论推导参见
[ODE 求解理论](../math/ode-theory.md)。

---

## 总体架构

ODE 求解的核心入口是 `dsolve`，它执行三步：

1. **规范化**（`normalize_ode`）：对输入方程做代数化简和标准化
2. **分类**（`classify_ode`）：按优先级返回所有适用的求解方法
3. **分派**：依次尝试各方法，返回第一个成功的结果

```
输入方程 ─→ normalize ─→ classify ─→ dispatch ─→ ODESolution
                              │
          ┌───────────────────┼───────────────────┐
          │                   │                   │
     一阶 ODE           二阶线性 ODE         级数解（回退）
    ┌────┴────┐          ┌────┴────┐
  可分离  线性  Bernoulli  常系数  Cauchy-Euler
  恰当    齐次   精确     降阶    参数变易
```

### 求解结果类型

`ODESolution` 枚举表示求解结果：

| 变体 | 含义 | 典型来源 |
|---|---|---|
| `Explicit(expr)` | 显式解 $y = f(x)$ | 线性、常系数、Laplace |
| `Implicit(expr)` | 隐式解 $F(x, y) = C$ | 可分离、恰当、齐次 |
| `Parametric(x(t), y(t))` | 参数化解 $(x(t), y(t))$ | 枚举预留（当前无求解器产生） |
| `Series(expr, n)` | 截断级数解 | 幂级数、Frobenius |
| `System(&[Atom])` | 系统解分量 $(y_1, y_2)$ | $2\times2$ 系统 |
| `Unsolved(ode)` | 无法解析求解 | 所有方法失败时 |

---

## 分类引擎

分类器（`classify.rs`）对 ODE 做结构分析，返回按优先级排列的候选方法列表。
不涉及求解——仅做模式识别。

### 检测流程

**一阶 ODE**（`order == 1`）：

1. 检查 `is_first_order_linear`：$y$ 和 $y'$ 仅以一次方出现 → `LinearFirst`
2. 检查 `is_bernoulli`：存在 $y^n$（$n \geq 2$）+ 线性 $y$ 项 → `Bernoulli`
3. 检查 `is_separable`：有不含 $y$ 的加法项 + 有含 $y$ 的加法项 → `Separable`
4. 检查 `is_exact`：$M + Ny' = 0$ 中 $\partial M/\partial y = \partial N/\partial x$ → `Exact`
5. 检查 `is_homogeneous`：所有项的 $(x, y)$ 总次数相同 → `Homogeneous`

**二阶及以上线性 ODE**：

6. `is_constant_coeff_linear`：$y, y', y''$ 的系数均不含 $x$ → `LinearConstantCoeff`
7. `is_cauchy_euler`：各项为 $c_k x^k y^{(k)}$ 形式 → `CauchyEuler`
8. 二阶线性 → `ReductionOfOrder`（尝试简单候选解）

**全局回退**：

9. 线性 ODE → `PowerSeries`（幂级数或 Frobenius）

### 精确方程的检测

精确方程判定 `is_exact` 实现如下：

1. 用 `split_mn` 将方程拆分为 $M(x, y) + N(x, y)\,y' = 0$
2. 将 $y(x)$ 替换为纯符号 $y$，使 $M, N$ 成为二元表达式
3. 计算 $\partial M/\partial y$ 和 $\partial N/\partial x$
4. 先比较规范化字符串（`normalize` 后的 `to_string()`）
5. 若不同，计算差值并用规则化简器检查是否为零

```rust
let dm_norm = normalize(ctx, dm_dy);
let dn_norm = normalize(ctx, dn_dx);
if dm_norm.to_string() == dn_norm.to_string() {
    return true;
}
// 回退：化简差值
let difference = simplify(ctx, ctx.add(&[dm_dy, ctx.mul(&[ctx.num(-1), dn_dx])]), &rules, 20);
matches!(normalize(ctx, difference).node(), AtomNode::Num(0))
```

这种双重检查策略弥补了化简器不合并同类项的局限。

### 齐次方程的检测

`is_homogeneous` 检查方程是否关于 $(x, y)$ 齐次（零次齐次性）。实现上，
`all_terms_homogeneous_degree` 遍历所有加法项，计算每项在 $x$ 和 `func` 上的总
次数（`term_degree`），断言所有项次数相同且为正。

对 $\text{Derivative}(y(x), x)$，按 $y$ 的次数 1 计算（不含 $x$ 的额外幂次）。

---

## 一阶 ODE 求解器

一阶求解器在 `first_order.rs` 中实现，共五种方法。

### 可分离方程

**判定**：方程可写为 $g(y)\,y' = f(x)$。

**算法**：

1. 用 `separate_by_func` 将方程拆为含 $y$ 的项（$y$-terms）和不含 $y$ 的项（$x$-terms）
2. 将 $y$-terms 中的 $x$ 替换为 $y$（`substitute_var`）
3. 对两侧分别积分：$\int g(y)\,dy$ 和 $\int f(x)\,dx$
4. 返回隐式解 $\int g(y)\,dy - \int f(x)\,dx = C$

**设计选择**：返回 `Implicit` 而非 `Explicit`，因为一般情形下需要反解 $y$，这
在符号计算中不可判定。积分调用 `crate::integral::integrate`，它启动完整的
积分管线（快速表→有理函数→Risch→三角重写→特殊函数→启发式技巧→未求值形式）。

### 一阶线性方程

**判定**：$y' + p(x)\,y = q(x)$，其中 $y$ 和 $y'$ 仅以一次方出现。

**算法**（积分因子法）：

1. 用 `extract_linear_coeffs` 提取 $p(x)$ 和 $q(x)$
2. 计算积分因子 $\mu(x) = e^{\int p\,dx}$
3. 通解：
$$y = e^{-\int p\,dx}\left(\int e^{\int p\,dx}\,q(x)\,dx + C_1\right)$$

```rust
let p_integral = integrate(ctx, p, var);
let mu = ctx.fun("exp", &[p_integral]);
let mu_q = ctx.mul(&[mu, q]);
let int_mu_q = integrate(ctx, mu_q, var);
let neg_p_integral = ctx.mul(&[ctx.num(-1), p_integral]);
let exp_neg_p = ctx.fun("exp", &[neg_p_integral]);
let particular = ctx.mul(&[exp_neg_p, int_mu_q]);
let c1 = ctx.var("C1");
let homogeneous = ctx.mul(&[c1, exp_neg_p]);
```

`extract_linear_coeffs` 将方程分解为加法项，分类为 $y'$-项、$y$-项和自由项。
对 $y$-项，通过过滤掉 `func` 因子提取 $p(x)$ 系数。

### Bernoulli 方程

**判定**：$y' + p(x)\,y = q(x)\,y^n$（$n \neq 0, 1$）。

**算法**：

1. 用 `find_bernoulli_power` 递归查找 $y^n$ 中的 $n$
2. 代换 $v = y^{1-n}$，方程化为 $v' + (1-n)\,p(x)\,v = (1-n)\,q(x)$
3. 用线性方程方法求解 $v$
4. 反代换 $y = v^{1/(1-n)}$

`find_power_inner` 递归遍历 AST，当遇到 `Pow(base, Num(n))` 且 `base` 含有
`func` 时返回 $n$。

**边界情况**：$n = 0$ 或 $n = 1$ 时退化为线性方程，返回 `None` 让 `LinearFirst`
处理。

### 恰当方程（Exact）

**判定**：$M(x,y) + N(x,y)\,y' = 0$ 且 $\partial M/\partial y = \partial N/\partial x$。

**算法**：

1. 用 `split_mn` 拆分 $M$ 和 $N$（$y'$ 必须线性出现，否则返回 `None`）
2. 检查恰当性（`partials_equal`）
3. 若不恰当，尝试积分因子（`find_integrating_factor`）
4. 对 $M$ 关于 $x$ 积分得 $F$ 的部分结果
5. 计算修正项 $g(y) = N - \partial F/\partial y$
6. 若 $g(y)$ 不含 $x$，对 $y$ 积分得完整解 $F + g = C$

**积分因子策略**：

- **候选 1**：$(M_y - N_x)/N$ 仅含 $x$ → $\mu(x) = \exp\!\int\!(M_y - N_x)/N\,dx$
- **候选 2**：$(N_x - M_y)/M$ 仅含 $y$ → $\mu(y) = \exp\!\int\!(N_x - M_y)/M\,dy$

`exp_simplify` 辅助化简 $\exp(k \cdot \log u) = u^k$，避免积分因子保持为
未化简的指数形式。

### 齐次方程

**判定**：$y' = f(y/x)$（所有加法项的 $(x,y)$ 总次数相同）。

**算法**：

1. 代换 $v = y/x$（即 $y = vx$），得 $y' = v + xv'$
2. 方程化为 $v + xv' = f(v)$，即 $xv' = f(v) - v$（关于 $v$ 和 $x$ 可分离）
3. 用 `separate_by_var` 拆分为含 $v$ 和含 $x$ 的项
4. 两侧积分，再反代换 $v = y/x$

返回 `Implicit` 解。

---

## 二阶线性 ODE 求解器

二阶求解器在 `second_order.rs` 中实现。

### 常系数方程

**形式**：$a\,y'' + b\,y' + c\,y = f(x)$，其中 $a, b, c$ 为常数。

**算法**：

1. 用 `extract_second_order_coeffs` 提取 $a, b, c$ 和 $f(x)$
2. 验证 $a, b, c$ 不含 $x$（`contains_x`）
3. 求特征方程 $ar^2 + br + c = 0$ 的判别式 $\Delta = b^2 - 4ac$
4. 构造齐次解 $y_c$
5. 构造特解 $y_p$（待定系数或参数变易）
6. 返回 $y = y_c + y_p$

#### 齐次解的三种情况

`constant_coeff_basis` 根据判别式构造基本解组：

| 判别式 | 基本解组 |
|---|---|
| $\Delta > 0$（两实根 $r_1, r_2$） | $e^{r_1 x},\; e^{r_2 x}$ |
| $\Delta = 0$（重根 $r$） | $e^{rx},\; x\,e^{rx}$ |
| $\Delta < 0$（复根 $\alpha \pm \beta i$） | $e^{\alpha x}\cos\beta x,\; e^{\alpha x}\sin\beta x$ |

实现中，判别式被折叠（`collect_terms`）为 `Num(d)` 时直接分支。对复根情况，
$\beta = \sqrt{-\Delta}/(2a)$：当 $-\Delta$ 是完全平方数时精确计算有理 $\beta$；
否则保留符号根式。

```rust
let sn = isqrt(num);
if sn * sn == num {
    // beta = sn/den — 精确有理数
    ctx.mul(&[ctx.num(sn), ctx.pow(ctx.num(den), ctx.num(-1))])
} else {
    // beta = sqrt(num)/den — 符号形式
    ctx.mul![
        ctx.pow(ctx.num(num), ctx.pow(ctx.num(2), ctx.num(-1))),
        ctx.pow(ctx.num(den), ctx.num(-1)),
    ]
}
```

#### 待定系数法

`particular_solution_undetermined` 按激励类型分派：

| 激励类型 | 检测方式 | 求解策略 |
|---|---|---|
| 多项式 $f(x)$ | `is_polynomial_in` | 设 $y_p = \sum A_k x^{k+s}$，回代求系数 |
| 指数 $F e^{kx}$ | `extract_exp_forcing` | 设 $y_p = A x^s e^{kx}$ |
| 三角 $f_c\cos\omega x + f_s\sin\omega x$ | `extract_trig_forcing` | 设 $y_p = x^s(A\cos\omega x + B\sin\omega x)$ |

$s$ 是共振移位：当 $k$（或 $0$）是特征方程的 $s$ 重根时，试探解乘以 $x^s$。

**多项式激励的回代算法**：

对 $y_p = \sum_{k=0}^{d} A_k x^{k+s}$，代入 $ay'' + by' + cy = f(x)$ 后，
$x$ 的同次幂系数给出递推关系。从最高次 $k = d$ 开始反向求解：

- $s = 0$（无共振）：$A_k = \frac{f_k - b(k+1)A_{k+1} - a(k+2)(k+1)A_{k+2}}{c}$
- $s = 1$（一重共振）：$A_k = \frac{f_k - a(k+2)(k+1)A_{k+1}}{b(k+1)}$
- $s = 2$（二重共振）：$A_k = \frac{f_k}{a(k+2)(k+1)}$

**叠加原理**：对多项式激励的加法项逐项求解再求和。

#### 参数变易法

当待定系数法不适用时，回退到参数变易（`variation_of_parameters`）：

给定齐次解 $y_1, y_2$ 和标准形方程 $y'' + py' + qy = g$：

1. 计算 Wronskian $W = y_1 y_2' - y_1' y_2$
2. $u_1' = -y_2 g / W$，$u_2' = y_1 g / W$
3. 积分得 $u_1, u_2$
4. $y_p = y_1 u_1 + y_2 u_2$

若积分返回未求值的 `Integral(...)` 形式，该方法失败（`is_integral_fallback`）。

### Cauchy–Euler 方程

**形式**：$a x^2 y'' + b x y' + c y = f(x)$。

**算法**：

1. 用 `extract_cauchy_euler_coeffs` 提取系数
2. 约化为常数系数（除以 $x$ 的幂次）
3. 求指标方程 $ar(r-1) + br + c = 0$，即 $ar^2 + (b-a)r + c = 0$
4. 根据判别式构造齐次解：

| 判别式 | 基本解组 |
|---|---|
| $\Delta > 0$（两实根 $r_1, r_2$） | $x^{r_1},\; x^{r_2}$ |
| $\Delta = 0$（重根 $r$） | $x^r,\; x^r \ln x$ |
| $\Delta < 0$（复根 $\alpha \pm \beta i$） | $x^\alpha\cos(\beta\ln x),\; x^\alpha\sin(\beta\ln x)$ |

5. 非齐次项用参数变易求特解（除以 $ax^2$ 化为标准形）

### 降阶法

**适用**：任何二阶线性 ODE（不要求常系数）。

**算法**：

1. 尝试候选解 $y_1 \in \{1, x, x^2, e^x, e^{-x}, e^{2x}\}$
2. 验证 $a y_1'' + b y_1' + c y_1 = 0$（`satisfies_extracted`）
3. 若 $y_1$ 满足齐次方程，第二解为：
$$y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$$
其中 $p = b/a$
4. 特解通过参数变易获得
5. 返回 $y = C_1 y_1 + C_2 y_2 + y_p$

**设计选择**：候选解集 $\{1, x, x^2, e^x, e^{-x}, e^{2x}\}$ 覆盖了教学和工程中最
常见的简单解形式。对更一般的方程，回退到 `PowerSeries`。

---

## Laplace 变换法

`laplace.rs` 实现初值问题的 Laplace 变换求解，通过 `dsolve_ivp` 入口调用。

### 适用范围

- 一阶：$a y' + b y = f(x)$，$y(0) = y_0$
- 二阶：$a y'' + b y' + c y = f(x)$，$y(0) = y_0$，$y'(0) = y_1$
- 系数 $a, b, c$ 必须为整数常数
- 激励 $f(x)$ 必须在变换表中

### 正变换

`laplace_kernel` 对 $x$-依赖的核函数做变换，支持的核：

| 核 | Laplace 变换 |
|---|---|
| $1$（常数） | $1/s$ |
| $x$ | $1/s^2$ |
| $x^n$（$n \leq 12$） | $n!/s^{n+1}$ |
| $e^{kx}$ | $1/(s-k)$ |
| $\sin\omega x$ | $\omega/(s^2 + \omega^2)$ |
| $\cos\omega x$ | $s/(s^2 + \omega^2)$ |
| $e^{kx}\sin\omega x$ | $\omega/((s-k)^2 + \omega^2)$ |
| $e^{kx}\cos\omega x$ | $(s-k)/((s-k)^2 + \omega^2)$ |

`split_const_factors` 将项拆分为 $x$-无关常数和 $x$-依赖核，仅变换核部分。
`linear_coeff_of` 从 $kx$ 中提取 $k$（$x$ 本身返回 $1$）。

### 代数方程构造

变换后的代数方程：

**一阶**：$a(sY - y_0) + bY = F(s)$
$$Y(s) = \frac{F(s) + a\,y_0}{as + b}$$

**二阶**：$a(s^2Y - sy_0 - y_1) + b(sY - y_0) + cY = F(s)$
$$Y(s) = \frac{F(s) + a(sy_0 + y_1) + by_0}{as^2 + bs + c}$$

### 逆变换

`inverse_laplace_term` 对有理分式 $Y(s) = \frac{n_1 s + n_0}{as^2 + bs + c}$ 做
部分分式分解和查表反变换。

先计算判别式 $\Delta = b^2 - 4ac$，按根的类型分三种情况：

**情况 1：两不同实根 $r_1 \neq r_2$**

$$\frac{n_1 s + n_0}{(s - r_1)(s - r_2)} = \frac{A}{s - r_1} + \frac{B}{s - r_2}$$

其中 $A = \frac{n_1 r_1 + n_0}{r_1 - r_2}$，$B = \frac{n_1 r_2 + n_0}{r_2 - r_1}$。

逆变换：$A\,e^{r_1 x} + B\,e^{r_2 x}$。

**情况 2：重根 $r$**（$\Delta = 0$）

$$\frac{n_1 s + n_0}{(s - r)^2} = \frac{n_1}{s - r} + \frac{n_1 r + n_0}{(s - r)^2}$$

逆变换：$n_1\,e^{rx} + (n_1 r + n_0)\,x\,e^{rx}$。

**情况 3：复根 $k \pm i\omega$**（$\Delta < 0$）

$$\frac{n_1 s + n_0}{(s - k)^2 + \omega^2} = e^{kx}\!\left[n_1\cos\omega x + \frac{n_0 + kn_1}{\omega}\sin\omega x\right]$$

$k = 0$ 时去掉 $e^{kx}$ 因子。要求 $(n_0 + kn_1) \bmod \omega = 0$（整数算术）。

### 整数算术约束

Laplace 模块全程使用 `i64` 整数运算（通过 `const_i64`、`quadratic_coeffs_i64`、
`linear_coeffs_i64`）。这限制了适用范围——系数含符号常数或无理数的方程会
回退到其他方法。设计理由：整数运算完全精确，不引入浮点误差。

---

## 幂级数与 Frobenius 方法

`series.rs` 实现两种级数解法，作为所有其他方法失败时的回退。

### 幂级数法

**适用**：线性 ODE 在常点 $x_0$ 附近的解。

**算法**（`solve_power_series`）：

1. 构造符号系数 $a_0, a_1, \ldots, a_{N-1}$
2. 手动构造级数 $S = \sum a_n (x - x_0)^n$ 及其一阶、二阶导数 $S', S''$
   （`build_series_triple` —— 不使用 `diff`，避免对 $a_n$ 微分）
3. 代入 ODE：$R(x) = a\,S'' + b\,S' + c\,S - f = 0$
4. 逐阶满足条件：对 $k = 0, 1, \ldots$，计算 $R^{(k)}(x_0) = 0$
5. 每个条件线性地确定一个新系数 $a_{k+\text{order}}$
6. 已求解的系数立即回代

**关键实现细节**：

- `substitute_series` 递归替换方程中所有 $y(x)$、$y'(x)$、$y''(x)$ 为级数形式
- `solve_linear_coeff` 从条件方程 $\text{coeff} \cdot a_n + \text{rest} = 0$ 解出 $a_n$
- 自由参数 $a_0, \ldots, a_{\text{order}-1}$ 保持为符号（初值条件）
- 若某条件不含待解系数但非零，说明 $x_0$ 不是常点 → 返回 `None`，交给 Frobenius

**限制**：默认 $N = 8$ 项。系数是初值的有理表达式（非数值）。

### Frobenius 方法

**适用**：二阶线性 ODE 在正则奇点 $x_0 = 0$ 附近的解。

**算法**（`solve_frobenius`）：

1. 提取系数 $a(x), b(x), c(x)$ 并验证 $f = 0$（仅齐次方程）
2. 设 $y = x^r \sum a_n x^n$，构造 $y, y', y''$
3. 代入 ODE，按 $x$ 的幂次分组残差
4. 最低幂次组给出**指标方程** $Ar^2 + Br + C = 0$
5. 求指标方程的有理根（`indicial_coeffs` → `isqrt` + 判别式）
6. 取较大根 $r_1$，代入后续组递推求 $a_n$

**$y$ 的构造**：

$$y = u \cdot S, \quad y' = r x^{-1} u S + u S', \quad y'' = r(r-1)x^{-2}uS + 2rx^{-1}uS' + uS''$$

其中 $u = x^r$ 是不透明占位符，$S = \sum a_n x^n$。`strip_x_and_u` 从残差项中
剥离 $x$ 的幂次和 $u$ 因子，用于分组。

**指标方程的提取**（`indicial_coeffs`）：

最低幂次组的项应为 $a_0$ 乘以 $r$ 的二次式。提取 $A, B, C$ 使得该组等于
$a_0(Ar^2 + Br + C)$。`expand_indicial` 处理含 $(r-1)$ 因子的展开。

**限制**：
- 仅支持 $x_0 = 0$
- 仅齐次方程
- 仅实有理根
- 仅取较大根的级数（较小根可能给出线性无关解）

---

## $2\times2$ 线性 ODE 系统

`systems.rs` 实现常系数线性系统的求解，通过 `dsolve_system` 入口调用。

### 输入格式

系统 $\mathbf{Y}' = A\mathbf{Y}$ 以两个方程表示：
$$\text{Derivative}(y_i, x) - (a_{i1}y_1 + a_{i2}y_2) = 0$$

`extract_coeff` 从每个方程提取 $2\times2$ 系数矩阵 $A$。

### 特征多项式

$$\lambda^2 - \text{tr}(A)\lambda + \det(A) = 0$$

其中 $\text{tr} = a_{11} + a_{22}$，$\det = a_{11}a_{22} - a_{12}a_{21}$。

判别式 $\Delta = \text{tr}^2 - 4\det$。

### 三种情况的求解

#### 不同实特征值

$\lambda_1 \neq \lambda_2$，均为整数。

1. `eigenvector` 求 $(A - \lambda I)\mathbf{v} = 0$：若 $a_{12} \neq 0$，取
   $\mathbf{v} = (a_{12},\; \lambda - a_{11})$
2. 通解：
$$\mathbf{Y} = C_1 \mathbf{v}_1 e^{\lambda_1 x} + C_2 \mathbf{v}_2 e^{\lambda_2 x}$$

#### 重特征值

$\lambda_1 = \lambda_2 = \lambda$。

**子情况 A**：$A = \lambda I$（完全矩阵）：
$$\mathbf{Y} = C_1 \mathbf{e}_1 e^{\lambda x} + C_2 \mathbf{e}_2 e^{\lambda x}$$

**子情况 B**：$A$ 亏损（defective）：
1. 求特征向量 $\mathbf{v}$（`eigenvector`）
2. 求广义特征向量 $\mathbf{w}$ 满足 $(A - \lambda I)\mathbf{w} = \mathbf{v}$
   （`generalized_eigenvector`）
3. 通解：
$$\mathbf{Y} = e^{\lambda x}\!\left[C_1 \mathbf{v} + C_2(x\mathbf{v} + \mathbf{w})\right]$$

`generalized_eigenvector` 通过解线性方程组 $(a_{11}-\lambda)w_1 + a_{12}w_2 = v_1$
寻找整数解，尝试 $w_1 = 0$ 和 $w_1 = 1$。

#### 复特征值

$\lambda = \alpha \pm \beta i$（要求 $\text{tr}$ 为偶数且 $-\Delta$ 为完全平方）。

1. $\alpha = \text{tr}/2$，$\beta = \sqrt{-\Delta}/2$
2. 特征向量 $(v + iw)$ 的实部 $p$ 和虚部 $q$：取 $v_1 = a_{12}$，
   $v_2 = (\alpha - a_{11}) + i\beta$，故 $p = (a_{12},\; \alpha - a_{11})$，
   $q = (0,\; \beta)$
3. 实基本解：
$$\mathbf{Y}_1 = e^{\alpha x}(p\cos\beta x - q\sin\beta x)$$
$$\mathbf{Y}_2 = e^{\alpha x}(p\sin\beta x + q\cos\beta x)$$

$\alpha = 0$ 时去掉 $e^{\alpha x}$ 因子，避免冗余的 `exp(0)`。

---

## 辅助工具

`util.rs` 提供 ODE 模块共用的基础设施。

### 阶数检测

`ode_order` 递归扫描 AST，查找 `Derivative(func, var, var, ...)` 节点，
返回最高微分阶数。

### 线性性检查

`is_linear_in` 检查表达式中 `func` 及其导数是否仅以一次方出现：
- `is_func_first_degree` 统计 `func`/`Derivative(func, ...)` 在乘法因子中的出现次数
- 若某因子含 `Pow(func, Num(n))` 且 $n \geq 2$，则非线性

### 同类项收集

`collect_terms` 解决化简器不合并同类项的问题：

1. 将加法项分解为有理系数 + 底数-指数因子对
2. 按底数-指数签名分组
3. 对每组求和有理系数
4. 系数为零的组丢弃

分解支持有理指数：$x \cdot x^{-1/2} = x^{1/2}$。

### 指数化简

`exp_simplify` 处理 $\exp(k \cdot \log u) = u^k$：
- 识别 `Mul(k, Fun("log", [u]))` 形式的指数
- 数值常数因子从 $\log$ 参数中丢弃（对 ODE 积分因子无影响）
- 回退到字面 `exp` 形式

### 解的替换与验证

`substitute_solution` 将候选解 $y = \text{sol}$ 代入方程：
- 替换 `func` → `sol`
- 替换 `Derivative(func, var)` → `diff(sol, var)`
- 替换 `Derivative(func, var, var)` → `diff(diff(sol, var), var)`

`verify_solution`（仅测试）代入后检查残差是否化简为零。

---

## 设计选择与权衡

### 化简策略

ODE 模块重度依赖 `ocas_rewrite::simplify`（规则化简，最多 20 轮迭代）
和 `ocas_atom::normalize`（结构规范化）。两者分工：

- `normalize`：快速的结构性等价判定（字符串比较）
- `simplify`：基于规则的代数化简（合并、展开、消除）

在恰当方程检测和解的构造中，两者交替使用以处理化简器无法自动发现的等价性。

### 积分依赖

ODE 求解器调用 `crate::integral::integrate` 来计算不定积分。这启动了完整的
积分管线（快速表→有理函数→Risch→三角重写→特殊函数→启发式技巧→未求值形式）。
当积分无法闭合求值时，返回 `Integral(...)` 形式——此时参数变易和降阶法标记为失败。

### 整数限制

Laplace 变换模块严格使用 `i64` 整数算术，拒绝符号系数。这是因为部分分式
分解需要精确的有理算术，而整数是 `i64` 上最安全的选择。对需要符号系数的
初值问题，应回退到通用 `dsolve` + 初始条件代入。

### 回退链

分类器返回的优先级顺序确保最具体的方法先尝试：

```
LinearFirst → Bernoulli → Separable → Exact → Homogeneous
    → LinearConstantCoeff → CauchyEuler → ReductionOfOrder
        → PowerSeries
```

`PowerSeries` 作为最后手段总是在列表末尾；其求解器内部先尝试幂级数法，
失败时再尝试 Frobenius 法。`dsolve` 依次尝试直到第一个成功。

---

## 源码位置

| 文件 | 职责 |
|---|---|
| `ocas-calc/src/ode/mod.rs` | 入口：`dsolve`、`dsolve_ivp`、`dsolve_system`、`ODE`、`ODESolution` |
| `ocas-calc/src/ode/classify.rs` | 分类器：`classify_ode`、`ODEType`、各项检测函数 |
| `ocas-calc/src/ode/first_order.rs` | 一阶求解器：可分离、线性、Bernoulli、恰当、齐次 |
| `ocas-calc/src/ode/second_order.rs` | 二阶求解器：常系数、Cauchy-Euler、降阶、参数变易 |
| `ocas-calc/src/ode/laplace.rs` | Laplace 变换：正变换表、逆变换（三种根类型）、IVP 求解 |
| `ocas-calc/src/ode/series.rs` | 幂级数和 Frobenius 方法 |
| `ocas-calc/src/ode/systems.rs` | $2\times2$ 系统：特征值分解、实/复/重根情况 |
| `ocas-calc/src/ode/util.rs` | 工具：阶数、线性性、同类项收集、指数化简、解替换 |

---

## 参见

- [ODE 求解理论](../math/ode-theory.md) — 各方法的数学推导和证明
- [符号微积分](../math/symbolic-calculus.md) — 微分和 Taylor 展开的理论基础
- [Rust API：求解器](../api/rust-solvers.md) — `dsolve`、`dsolve_ivp`、`dsolve_system` 的 API 参考
- [Rust API：微积分](../api/rust-calculus.md) — `diff`、`integrate` 的 API 参考
