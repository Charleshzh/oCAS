# 进阶：ODE 求解理论

## 前提知识

- [多项式代数](./polynomial-algebra.md) — 多项式运算与因式分解
- [线性代数](./linear-algebra.md) — 矩阵运算与特征值
- [符号微积分](./symbolic-calculus.md) — 求导与积分的基本规则

---

## 基础概念

### 常微分方程的定义

**常微分方程**（ordinary differential equation, ODE）是包含未知函数 $y(x)$ 及其导数的方程。一般形式为

$$F\!\left(x,\, y,\, y',\, y'',\, \dots,\, y^{(n)}\right) = 0$$

其中 $y^{(k)} = \dfrac{d^k y}{dx^k}$。

### 阶（order）

ODE 中出现的**最高阶导数**决定了方程的阶。例如：

| 方程 | 阶 |
|---|---|
| $y' + y = 0$ | 1 |
| $y'' + 2y' + y = \sin x$ | 2 |
| $y''' - y = 0$ | 3 |

### 线性与非线性

一个 $n$ 阶 ODE 是**线性**的，当且仅当它可以写成

$$a_n(x)\, y^{(n)} + a_{n-1}(x)\, y^{(n-1)} + \cdots + a_1(x)\, y' + a_0(x)\, y = g(x)$$

其中 $a_k(x)$ 和 $g(x)$ 仅依赖于自变量 $x$，不依赖于 $y$。若 $g(x) \equiv 0$，则称为**齐次**线性 ODE；否则称为**非齐次**。

含有 $y^2$、$\sin y$、$y\,y'$ 等项的方程是**非线性**的。

### 初值条件与边值条件

- **初值问题**（IVP）：给定 $y(x_0) = y_0$（一阶），或 $y(x_0) = y_0,\; y'(x_0) = y_1$（二阶），确定唯一解。
- **边值问题**（BVP）：在两个不同点指定条件，如 $y(0) = 0,\; y(\pi) = 0$。

### 存在唯一性定理（Picard–Lindelöf）

**定理**（一阶情形）. 设 $f(x, y)$ 在矩形区域 $R = \{(x,y) : |x - x_0| \leq a,\; |y - y_0| \leq b\}$ 上连续，且关于 $y$ 满足 Lipschitz 条件

$$|f(x, y_1) - f(x, y_2)| \leq L\, |y_1 - y_2|$$

则初值问题 $y' = f(x,y),\; y(x_0) = y_0$ 在 $|x - x_0| \leq h$（$h = \min(a,\, b/M)$，$M = \max_R |f|$）上存在**唯一**解。

**推广到 $n$ 阶**. 将 $n$ 阶方程化为一阶方程组 $\mathbf{Y}' = \mathbf{F}(x, \mathbf{Y})$，对 $\mathbf{F}$ 施加同样的 Lipschitz 条件即可。

这个定理保证了：当右端函数足够"好"时，ODE 的解不仅存在，而且唯一。这是所有符号求解器的理论基础——如果符号方法找到了一个解，Picard–Lindelöf 保证了它就是*唯一的*解。

---

## 核心理论

### 一阶 ODE

#### 可分离变量方程（Separable）

**形式**：

$$g(y)\, y' = f(x)$$

**求解步骤**：

1. 将方程写为 $g(y)\, dy = f(x)\, dx$
2. 两边积分：$\displaystyle\int g(y)\, dy = \int f(x)\, dx + C$

**例子**. 求解 $y' = xy$：

$$\frac{dy}{y} = x\, dx \quad\Longrightarrow\quad \ln|y| = \frac{x^2}{2} + C \quad\Longrightarrow\quad y = C_1\, e^{x^2/2}$$

**oCAS 中的实现**. 分类器 `is_separable` 通过检查方程是否能写成"所有含 $y$ 的项伴随 $y'$，其余项仅含 $x$"来判断可分离性。求解器 `solve_separable` 将方程拆为 $(y\text{-项},\, x\text{-项})$，分别积分后构造隐式解 $\int g\, dy - \int f\, dx = C$。

---

#### 一阶线性方程（Linear First-Order）

**形式**：

$$y' + p(x)\, y = q(x)$$

**积分因子法**. 乘以积分因子

$$\mu(x) = e^{\int p(x)\, dx}$$

方程化为

$$\frac{d}{dx}\!\bigl[\mu(x)\, y\bigr] = \mu(x)\, q(x)$$

积分得

$$y = \frac{1}{\mu(x)} \int \mu(x)\, q(x)\, dx + \frac{C}{\mu(x)}$$

**推导**. 观察到 $\mu' = p\,\mu$，因此

$$\mu\, y' + p\,\mu\, y = \mu\, y' + \mu'\, y = (\mu\, y)'$$

**例子**. 求解 $y' - \dfrac{2}{x}\, y = x^2$：

$$p(x) = -\frac{2}{x}, \quad \mu = e^{\int -2/x\, dx} = e^{-2\ln x} = x^{-2}$$

$$\frac{d}{dx}\!\bigl[x^{-2} y\bigr] = x^{-2} \cdot x^2 = 1$$

$$x^{-2}\, y = x + C \quad\Longrightarrow\quad y = x^3 + Cx^2$$

**oCAS 中的实现**. `solve_linear_first` 提取 $p(x)$ 和 $q(x)$，计算 $\mu = \exp\!\bigl(\int p\, dx\bigr)$，然后返回

$$y = \mu^{-1}\!\left(\int \mu\, q\, dx + C\right)$$

---

#### Bernoulli 方程

**形式**：

$$y' + p(x)\, y = q(x)\, y^n, \qquad n \neq 0,\, 1$$

**变换**. 令 $v = y^{1-n}$，则 $v' = (1-n)\, y^{-n}\, y'$。原方程两边乘以 $(1-n)\, y^{-n}$：

$$(1-n)\, y^{-n}\, y' + (1-n)\, p(x)\, y^{1-n} = (1-n)\, q(x)$$

即

$$v' + (1-n)\, p(x)\, v = (1-n)\, q(x)$$

这是关于 $v$ 的一阶线性方程，用积分因子法求解后回代 $v = y^{1-n}$。

**例子**. 求解 $y' + \dfrac{1}{x}\, y = x\, y^2$（$n = 2$）：

令 $v = y^{-1}$，$v' = -y^{-2}\, y'$，得

$$-v' + \frac{1}{x}\, v = x \quad\Longrightarrow\quad v' - \frac{1}{x}\, v = -x$$

积分因子 $\mu = x^{-1}$，解得 $v = -x^2 + Cx$，即 $y = \dfrac{1}{Cx - x^2}$。

**oCAS 中的实现**. `solve_bernoulli` 检测非线性幂 $n$，构造 $v = y^{1-n}$ 替换，调用 `solve_linear_first` 求解线性化后的方程，然后回代。

---

#### 恰当方程（Exact Equations）

**形式**：

$$M(x,y) + N(x,y)\, y' = 0$$

若存在势函数 $F(x,y)$ 使得 $\dfrac{\partial F}{\partial x} = M$，$\dfrac{\partial F}{\partial y} = N$，则方程等价于 $\dfrac{d}{dx} F(x, y(x)) = 0$，解为 $F(x,y) = C$。

**恰当性判据**. 方程恰当当且仅当

$$\frac{\partial M}{\partial y} = \frac{\partial N}{\partial x}$$

**求解步骤**：

1. 由 $\partial F/\partial x = M$ 对 $x$ 积分：$F = \int M\, dx + h(y)$
2. 由 $\partial F/\partial y = N$ 确定 $h(y)$：$h'(y) = N - \partial/\partial y \int M\, dx$
3. 代入得 $F(x,y) = C$

**非恰当方程的积分因子**. 若方程不恰当，可尝试寻找积分因子 $\mu$：

- 若 $\dfrac{M_y - N_x}{N}$ 仅依赖于 $x$，则 $\mu(x) = \exp\!\left(\int \dfrac{M_y - N_x}{N}\, dx\right)$
- 若 $\dfrac{N_x - M_y}{M}$ 仅依赖于 $y$，则 $\mu(y) = \exp\!\left(\int \dfrac{N_x - M_y}{M}\, dy\right)$

**例子**. 求解 $(2xy + 3) + (x^2 + 4y)\, y' = 0$：

$M_y = 2x = N_x$，方程恰当。

$$F = \int (2xy + 3)\, dx = x^2 y + 3x + h(y)$$

$$F_y = x^2 + h'(y) = x^2 + 4y \quad\Longrightarrow\quad h'(y) = 4y \quad\Longrightarrow\quad h = 2y^2$$

解：$x^2 y + 3x + 2y^2 = C$。

**oCAS 中的实现**. `solve_exact` 首先调用 `partials_equal` 验证恰当性。若不恰当，`find_integrating_factor` 尝试上述两种积分因子（依赖于 $x$ 或依赖于 $y$ 的情况）。找到势函数后返回隐式解。

---

#### 齐次方程（Homogeneous）

**形式**：

$$y' = f\!\left(\frac{y}{x}\right)$$

等价地，方程右端 $f(x,y)$ 是关于 $(x,y)$ 的零次齐次函数：$f(tx, ty) = f(x,y)$。

**变换**. 令 $v = y/x$，则 $y = vx$，$y' = v + xv'$。代入：

$$v + x\, v' = f(v) \quad\Longrightarrow\quad x\, v' = f(v) - v$$

这是关于 $v$ 和 $x$ 的可分离方程：

$$\frac{dv}{f(v) - v} = \frac{dx}{x}$$

积分后回代 $v = y/x$。

**例子**. 求解 $y' = \dfrac{y + x}{x} = \dfrac{y}{x} + 1$：

令 $v = y/x$，$v + xv' = v + 1$，$xv' = 1$，$v = \ln|x| + C$，$y = x\ln|x| + Cx$。

**oCAS 中的实现**. 分类器 `is_homogeneous` 要求方程是一阶**线性**的，且所有加法项关于 $(x, y)$ 具有相同总次数（零次齐次性），并存在不含 $y$ 的项。求解器 `solve_homogeneous` 执行 $v = y/x$ 替换，化为可分离方程后求解。

---

### 二阶线性 ODE

二阶线性 ODE 的一般形式为

$$a(x)\, y'' + b(x)\, y' + c(x)\, y = f(x)$$

**叠加原理**. 若 $y_1, y_2$ 是齐次方程的两个线性无关解，$y_p$ 是非齐次方程的一个特解，则通解为

$$y = C_1\, y_1 + C_2\, y_2 + y_p$$

#### 常系数方程（Constant Coefficients）

**形式**：

$$a\, y'' + b\, y' + c\, y = f(x), \qquad a, b, c \in \mathbb{R},\; a \neq 0$$

**特征方程**. 齐次方程 $ay'' + by' + cy = 0$ 的特征方程为

$$a\, r^2 + b\, r + c = 0$$

判别式 $\Delta = b^2 - 4ac$ 决定了解的形式：

| 判别式 | 特征根 | 齐次通解 |
|---|---|---|
| $\Delta > 0$ | $r_1, r_2 \in \mathbb{R},\; r_1 \neq r_2$ | $C_1\, e^{r_1 x} + C_2\, e^{r_2 x}$ |
| $\Delta = 0$ | $r_1 = r_2 = r$ | $(C_1 + C_2\, x)\, e^{rx}$ |
| $\Delta < 0$ | $\alpha \pm \beta i$ | $e^{\alpha x}(C_1 \cos\beta x + C_2 \sin\beta x)$ |

**待定系数法**（Undetermined Coefficients）. 当 $f(x)$ 是多项式、指数函数 $e^{kx}$、三角函数 $\sin\omega x$、$\cos\omega x$ 或其乘积/和时，可以猜测特解形式。

| $f(x)$ 的形式 | 猜测的 $y_p$ 形式 |
|---|---|
| $P_n(x)$（$n$ 次多项式） | $Q_n(x) = a_0 + a_1 x + \cdots + a_n x^n$ |
| $F\, e^{kx}$ | $A\, e^{kx}$ |
| $f_c \cos\omega x + f_s \sin\omega x$ | $A\cos\omega x + B\sin\omega x$ |
| $F\, e^{kx} \cos\omega x$ | $e^{kx}(A\cos\omega x + B\sin\omega x)$ |

**共振（resonance）**. 当猜测形式与齐次解线性相关时，需乘以 $x^s$（$s$ 为 $k$ 作为特征根的重数）：

- 若 $k$ 是单特征根：$y_p = Ax\, e^{kx}$
- 若 $k$ 是二重特征根：$y_p = Ax^2\, e^{kx}$

对于多项式强迫，若 $r = 0$ 是 $s$ 重特征根，则猜测形式需乘以 $x^s$。

**例子**. 求解 $y'' - 3y' + 2y = e^{3x}$：

特征方程 $r^2 - 3r + 2 = 0$，$r_1 = 1, r_2 = 2$。

齐次通解：$C_1 e^x + C_2 e^{2x}$。

$k = 3$ 不是特征根，猜测 $y_p = Ae^{3x}$。代入：$9A - 9A + 2A = 1$，$A = 1/2$。

通解：$y = C_1 e^x + C_2 e^{2x} + \dfrac{1}{2} e^{3x}$。

**oCAS 中的实现**. `solve_constant_coeff` 提取系数 $(a, b, c)$ 和强迫项 $f(x)$。`constant_coeff_basis` 根据判别式构造齐次基解。`particular_solution_undetermined` 处理多项式、指数和三角强迫，包含完整的共振检测。

---

#### Cauchy–Euler 方程

**形式**：

$$a\, x^2\, y'' + b\, x\, y' + c\, y = f(x)$$

系数中 $y^{(k)}$ 的系数恰好是 $x^k$ 的常数倍。

**变换**. 令 $x = e^t$（即 $t = \ln x$），利用

$$y' = \frac{1}{x}\, \frac{dy}{dt}, \qquad y'' = \frac{1}{x^2}\!\left(\frac{d^2 y}{dt^2} - \frac{dy}{dt}\right)$$

原方程化为常系数方程

$$a\, \ddot{y} + (b - a)\, \dot{y} + c\, y = f(e^t)$$

其中 $\dot{y} = dy/dt$。

**指标方程**. 对齐次方程，直接猜测 $y = x^r$，代入得指标方程

$$a\, r(r-1) + b\, r + c = 0 \quad\Longleftrightarrow\quad a\, r^2 + (b - a)\, r + c = 0$$

解的结构与常系数方程类似（根据判别式分三种情况），但基本解为 $x^r$ 而非 $e^{rx}$：

| 判别式 | 基本解 |
|---|---|
| $\Delta > 0$ | $x^{r_1},\; x^{r_2}$ |
| $\Delta = 0$ | $x^r,\; x^r \ln x$ |
| $\Delta < 0$（$\alpha \pm \beta i$） | $x^\alpha \cos(\beta\ln x),\; x^\alpha \sin(\beta\ln x)$ |

**例子**. 求解 $x^2 y'' - 2xy' + 2y = 0$：

指标方程：$r(r-1) - 2r + 2 = r^2 - 3r + 2 = 0$，$r_1 = 1, r_2 = 2$。

通解：$y = C_1 x + C_2 x^2$。

**oCAS 中的实现**. 分类器 `is_cauchy_euler` 检查方程是否符合 $c_k x^k y^{(k)}$ 的模式。`solve_cauchy_euler` 提取 Cauchy–Euler 系数，解指标方程，然后根据判别式构造通解。`cauchy_euler_basis` 处理三种判别式情况。

---

#### 降阶法（Reduction of Order）

当已知齐次方程的一个解 $y_1$ 时，第二个线性无关解可通过降阶获得。

**方法**. 设 $y_2 = v(x)\, y_1$，代入 $y'' + p\, y' + q\, y = 0$。利用 $y_1$ 满足齐次方程的条件，化简后得到关于 $v'$ 的一阶方程：

$$v'' y_1 + v'(2y_1' + p\, y_1) = 0$$

令 $w = v'$，分离变量：

$$\frac{w'}{w} = -\frac{2y_1' + p\, y_1}{y_1} = -2\, \frac{y_1'}{y_1} - p$$

积分得

$$w = \frac{e^{-\int p\, dx}}{y_1^2}$$

因此

$$y_2 = y_1 \int \frac{e^{-\int p\, dx}}{y_1^2}\, dx$$

**例子**. $x^2 y'' - 2xy' + 2y = 0$，已知 $y_1 = x$。

标准形：$y'' - \dfrac{2}{x}\, y' + \dfrac{2}{x^2}\, y = 0$，$p = -2/x$。

$$y_2 = x \int \frac{e^{\int 2/x\, dx}}{x^2}\, dx = x \int \frac{x^2}{x^2}\, dx = x \cdot x = x^2$$

**oCAS 中的实现**. `solve_reduction_of_order` 首先尝试简单候选解（$1, x, x^2, e^x, e^{-x}, e^{2x}$）代入齐次方程。找到 $y_1$ 后，用上述公式计算 $y_2$。对于非齐次方程，再用参数变易法求特解。

---

#### 参数变易法（Variation of Parameters）

当待定系数法不适用（如 $f(x)$ 不是标准形式）时，参数变易法是通用方法。

**方法**. 已知齐次方程的两个线性无关解 $y_1, y_2$，设非齐次方程的特解为

$$y_p = u_1(x)\, y_1 + u_2(x)\, y_2$$

其中 $u_1, u_2$ 满足

$$\begin{cases} u_1'\, y_1 + u_2'\, y_2 = 0 \\ u_1'\, y_1' + u_2'\, y_2' = g(x) \end{cases}$$

这里 $g(x) = f(x)/a$ 是标准形的右端。由 Cramer 法则：

$$u_1' = -\frac{y_2\, g}{W}, \qquad u_2' = \frac{y_1\, g}{W}$$

其中 $W = y_1\, y_2' - y_1'\, y_2$ 是 **Wronskian**。

因此

$$y_p = -y_1 \int \frac{y_2\, g}{W}\, dx + y_2 \int \frac{y_1\, g}{W}\, dx$$

**Wronskian 的意义**. $y_1, y_2$ 线性无关当且仅当 $W \neq 0$。对于二阶齐次线性 ODE，Abel 公式给出 $W(x) = W(x_0)\, \exp\!\bigl(-\int_{x_0}^x p\, dt\bigr)$，因此 $W$ 要么恒零要么恒非零。

**oCAS 中的实现**. `variation_of_parameters` 接受两个基解 $y_1, y_2$ 和标准形右端 $g$，计算 Wronskian，然后用 oCAS 的符号积分器（`integrate`）计算 $u_1, u_2$。若积分无法化为闭形（返回未求值的 `Integral(...)` 形式），则放弃此方法。

---

### 级数解

当 ODE 的系数不是常数且无法用初等方法求解时，级数解法提供了构造解的系统方法。

#### 幂级数解法（Ordinary Point）

**适用条件**. 考虑二阶线性 ODE

$$y'' + p(x)\, y' + q(x)\, y = 0$$

若 $p(x)$ 和 $q(x)$ 在 $x = x_0$ 处解析（即可以展开为收敛幂级数），则 $x_0$ 是**常点**（ordinary point），方程在 $x_0$ 附近有两个线性无关的解析解。

**方法**. 设

$$y = \sum_{n=0}^{\infty} a_n (x - x_0)^n$$

则

$$y' = \sum_{n=1}^{\infty} n\, a_n (x - x_0)^{n-1}, \qquad y'' = \sum_{n=2}^{\infty} n(n-1)\, a_n (x - x_0)^{n-2}$$

将 $y, y', y''$ 和 $p(x), q(x)$ 的幂级数展开代入方程，合并同次幂，令各次幂系数为零，得到关于 $a_n$ 的**递推关系**。

对于 $k$ 阶 ODE，初始系数 $a_0, a_1, \dots, a_{k-1}$ 是自由参数（对应通解中的 $k$ 个任意常数），后续系数由递推确定。

**例子**. 求解 $y'' - xy = 0$（Airy 方程）在 $x = 0$ 附近的级数解：

$$\sum_{n=2}^{\infty} n(n-1)a_n x^{n-2} - x \sum_{n=0}^{\infty} a_n x^n = 0$$

$$\sum_{n=0}^{\infty} (n+2)(n+1)a_{n+2} x^n - \sum_{n=0}^{\infty} a_n x^{n+1} = 0$$

$x^0$ 项：$2 \cdot 1 \cdot a_2 = 0$，$a_2 = 0$。

$x^n$（$n \geq 1$）项：$(n+2)(n+1)\, a_{n+2} - a_{n-1} = 0$，即

$$a_{n+2} = \frac{a_{n-1}}{(n+2)(n+1)}$$

取 $a_0 = 1, a_1 = 0$：$a_3 = a_0/(3 \cdot 2) = 1/6$，$a_6 = a_3/(6 \cdot 5) = 1/180$，……

取 $a_0 = 0, a_1 = 1$：$a_4 = a_1/(4 \cdot 3) = 1/12$，$a_7 = a_4/(7 \cdot 6) = 1/504$，……

**oCAS 中的实现**. `solve_power_series` 设置 $y = \sum a_n x^n$，构造 $y, y', y''$ 的级数三元组（`build_series_triple`），代入方程，按 $x$ 的幂次分组，解线性方程求各 $a_n$（`solve_linear_coeff`）。返回截断级数 $ODESolution::Series(expr, n\_terms)$。

---

#### Frobenius 方法（Regular Singular Point）

**适用条件**. 若 $x = x_0$ 不是常点，但 $(x - x_0)\, p(x)$ 和 $(x - x_0)^2\, q(x)$ 在 $x_0$ 处解析，则 $x_0$ 是**正则奇点**（regular singular point）。

**方法**. 设

$$y = x^r \sum_{n=0}^{\infty} a_n x^n = \sum_{n=0}^{\infty} a_n x^{n+r}, \qquad a_0 \neq 0$$

代入方程，最低次幂 $x^{r-2}$（对于 Cauchy–Euler 型的主部）的系数给出**指标方程**（indicial equation）：

$$A\, r(r-1) + B\, r + C = 0$$

其中 $A, B, C$ 由方程系数在 $x_0$ 处的值确定。后续幂次给出 $a_n$ 的递推关系。

**指标根的关系**. 设指标方程的两个根为 $r_1, r_2$（$r_1 \geq r_2$）：

| 根的情况 | 解的形式 |
|---|---|
| $r_1 - r_2 \notin \mathbb{Z}$ | 两个独立级数解 $y_1 = x^{r_1}\sum a_n x^n$，$y_2 = x^{r_2}\sum b_n x^n$ |
| $r_1 = r_2$ | $y_1 = x^{r_1}\sum a_n x^n$，$y_2 = y_1 \ln x + x^{r_1}\sum c_n x^n$ |
| $r_1 - r_2 \in \mathbb{Z}^+$ | $y_1$ 如上，$y_2$ 可能含 $\ln x$ 项（对数情形） |

**例子**. 求解 $x^2 y'' + xy' + (x^2 - 1)y = 0$（Bessel 方程，$\nu = 1$）在 $x = 0$ 附近：

$p(x) = 1/x$，$q(x) = 1 - 1/x^2$。$x = 0$ 是正则奇点。

设 $y = \sum a_n x^{n+r}$，代入后 $x^r$ 的系数为

$$a_0[r(r-1) + r - 1] = a_0(r^2 - 1) = 0$$

指标方程：$r^2 - 1 = 0$，$r_1 = 1, r_2 = -1$。

$r_1 - r_2 = 2 \in \mathbb{Z}$，需检查对数情形。

**oCAS 中的实现**. `solve_frobenius` 设置 $y = x^r \sum a_n x^n$（$r$ 为符号变量），代入方程，按 $x$ 的幂次分组。最低幂次组通过 `indicial_coeffs` 提取指标方程 $Ar^2 + Br + C = 0$ 的整数系数。解出 $r$ 后，后续组递推确定 $a_n$。目前仅处理指标方程有**有理实根**的情况，返回较大根对应的级数。

---

### Laplace 变换

Laplace 变换是求解线性常系数 ODE 初值问题的强有力工具，尤其适用于分段连续或含有脉冲函数的强迫项。

#### 定义

函数 $f(x)$（$x \geq 0$）的 Laplace 变换定义为

$$\mathcal{L}\{f\}(s) = F(s) = \int_0^{\infty} e^{-sx}\, f(x)\, dx$$

#### 导数的变换

利用分部积分：

$$\mathcal{L}\{y'\}(s) = s\, Y(s) - y(0)$$

$$\mathcal{L}\{y''\}(s) = s^2\, Y(s) - s\, y(0) - y'(0)$$

一般地

$$\mathcal{L}\{y^{(n)}\}(s) = s^n\, Y(s) - \sum_{k=0}^{n-1} s^{n-1-k}\, y^{(k)}(0)$$

#### 常用变换对

| $f(x)$ | $F(s) = \mathcal{L}\{f\}(s)$ |
|---|---|
| $1$ | $\dfrac{1}{s}$ |
| $x^n$ | $\dfrac{n!}{s^{n+1}}$ |
| $e^{kx}$ | $\dfrac{1}{s - k}$ |
| $\cos\omega x$ | $\dfrac{s}{s^2 + \omega^2}$ |
| $\sin\omega x$ | $\dfrac{\omega}{s^2 + \omega^2}$ |
| $e^{kx}\cos\omega x$ | $\dfrac{s - k}{(s-k)^2 + \omega^2}$ |
| $e^{kx}\sin\omega x$ | $\dfrac{\omega}{(s-k)^2 + \omega^2}$ |
| $x\, e^{kx}$ | $\dfrac{1}{(s-k)^2}$ |

#### IVP 求解步骤

以二阶方程 $ay'' + by' + cy = f(x)$，$y(0) = y_0$，$y'(0) = y_1$ 为例：

1. **变换**：对方程两边取 Laplace 变换

$$a[s^2 Y - sy_0 - y_1] + b[sY - y_0] + cY = F(s)$$

2. **解代数方程**：

$$Y(s) = \frac{F(s) + a(sy_0 + y_1) + by_0}{as^2 + bs + c}$$

3. **部分分式分解**：将 $Y(s)$ 分解为标准变换对的线性组合

4. **逆变换**：查表得到 $y(x) = \mathcal{L}^{-1}\{Y(s)\}(x)$

**逆变换的关键技术**. 对于分母为二次多项式的有理函数：

- **互异实根** $r_1 \neq r_2$：$\dfrac{n_1 s + n_0}{(s-r_1)(s-r_2)} = \dfrac{A}{s-r_1} + \dfrac{B}{s-r_2}$，$A = \dfrac{n_1 r_1 + n_0}{r_1 - r_2}$

- **重根** $r$：$\dfrac{n_1 s + n_0}{(s-r)^2} = \dfrac{n_1}{s-r} + \dfrac{n_1 r + n_0}{(s-r)^2}$

- **共轭复根** $k \pm i\omega$：$\dfrac{n_1 s + n_0}{(s-k)^2 + \omega^2} = e^{kx}\!\left[n_1 \cos\omega x + \dfrac{n_0 + kn_1}{\omega}\sin\omega x\right]$

**例子**. 求解 $y'' - 3y' + 2y = 4$，$y(0) = 2$，$y'(0) = 3$：

$$s^2 Y - 2s - 3 - 3(sY - 2) + 2Y = \frac{4}{s}$$

$$(s^2 - 3s + 2)Y = \frac{4}{s} + 2s - 3 = \frac{2s^2 - 3s + 4}{s}$$

$$Y = \frac{2s^2 - 3s + 4}{s(s-1)(s-2)}$$

部分分式：$Y = \dfrac{2}{s} - \dfrac{3}{s-1} + \dfrac{3}{s-2}$。

逆变换：$y(x) = 2 - 3e^x + 3e^{2x}$。

验证：$y(0) = 2$，$y'(0) = -3 + 6 = 3$，且 $y'' - 3y' + 2y = 4$ ✓。

**oCAS 中的实现**. `dsolve_ivp` 是 IVP 求解的入口，内部调用 `solve_laplace`。

`solve_laplace` 的工作流程：

1. 判断方程为一阶或二阶线性常系数方程
2. 提取系数 $(a, b, c)$ 和强迫项 $f(x)$
3. `laplace_transform` 计算 $\mathcal{L}\{f\}(s)$ — 支持多项式、指数、三角函数及其乘积
4. 构造 $Y(s) = \dfrac{F(s) + \text{初值项}}{as^2 + bs + c}$
5. `inverse_laplace` 通过 `split_fraction` 分解有理函数，对每项调用 `inverse_distinct_roots`、`inverse_repeated_root` 或 `inverse_complex_roots`
6. 返回 $ODESolution::Explicit(y(x))$（不含自由常数）

---

### 2×2 常系数线性方程组

考虑一阶常系数线性方程组

$$\mathbf{Y}' = A\, \mathbf{Y}, \qquad \mathbf{Y} = \begin{pmatrix} y_1 \\ y_2 \end{pmatrix}, \quad A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$$

#### 通解：矩阵指数

形式解为 $\mathbf{Y}(x) = e^{Ax}\, \mathbf{Y}(0)$。实际求解依赖于 $A$ 的特征值。

#### 特征值分解

$A$ 的特征方程为

$$\lambda^2 - (a_{11} + a_{22})\lambda + (a_{11}a_{22} - a_{12}a_{21}) = 0$$

设特征值为 $\lambda_1, \lambda_2$，对应特征向量为 $\mathbf{v}_1, \mathbf{v}_2$。

**情形 1：互异实特征值**（$\lambda_1 \neq \lambda_2 \in \mathbb{R}$）

$$\mathbf{Y} = C_1\, e^{\lambda_1 x}\, \mathbf{v}_1 + C_2\, e^{\lambda_2 x}\, \mathbf{v}_2$$

**情形 2：重特征值**（$\lambda_1 = \lambda_2 = \lambda$）

若 $A$ 不可对角化（只有一个线性无关特征向量 $\mathbf{v}$），需要广义特征向量 $\mathbf{w}$ 满足 $(A - \lambda I)\mathbf{w} = \mathbf{v}$：

$$\mathbf{Y} = C_1\, e^{\lambda x}\, \mathbf{v} + C_2\, e^{\lambda x}(x\, \mathbf{v} + \mathbf{w})$$

**情形 3：共轭复特征值**（$\lambda = \alpha \pm \beta i$，$\beta \neq 0$）

特征向量也有实部和虚部：$\mathbf{v} = \mathbf{p} + i\mathbf{q}$。实值通解为

$$\mathbf{Y} = e^{\alpha x}\bigl[C_1(\mathbf{p}\cos\beta x - \mathbf{q}\sin\beta x) + C_2(\mathbf{p}\sin\beta x + \mathbf{q}\cos\beta x)\bigr]$$

**例子**. 求解谐振子系统 $y_1' = y_2,\; y_2' = -y_1$：

$$A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$$

特征方程：$\lambda^2 + 1 = 0$，$\lambda = \pm i$（$\alpha = 0, \beta = 1$）。

$\lambda = i$ 的特征向量：$(A - iI)\mathbf{v} = 0$，$\mathbf{v} = (1, i)^T$，$\mathbf{p} = (1, 0)^T$，$\mathbf{q} = (0, 1)^T$。

$$\mathbf{Y} = C_1 \begin{pmatrix} \cos x \\ -\sin x \end{pmatrix} + C_2 \begin{pmatrix} \sin x \\ \cos x \end{pmatrix}$$

即 $y_1 = C_1\cos x + C_2\sin x$，$y_2 = -C_1\sin x + C_2\cos x$。

**oCAS 中的实现**. `solve_linear_system` 处理 $2 \times 2$ 常系数系统：

1. 从方程组中提取矩阵 $A$ 的各元素（`extract_coeff`）
2. 计算特征多项式，求特征值（整数或有理数范围内）
3. 对互异实特征值：`eigenvector` 求整数特征向量
4. 对重特征值：`generalized_eigenvector` 求广义特征向量 $(A - \lambda I)\mathbf{w} = \mathbf{v}$
5. 对复特征值：`solve_complex_2x2` 构造实值基本解
6. 返回 $ODESolution::System([y_1, y_2])$，其中含自由常数 $C_1, C_2$

---

## 在 oCAS 中的实现

### 分类器（classify.rs）

`classify_ode` 是 ODE 求解的入口决策器。它分析方程结构，返回按优先级排列的候选方法列表：

| 优先级 | `ODEType` | 检测条件 |
|---|---|---|
| 1 | `LinearFirst` | 一阶，$y$ 和 $y'$ 线性出现，无 $y^2, \sin y$ 等非线性项，且含 $y$ |
| 2 | `Bernoulli` | 一阶，含 $y^n$（$n \geq 2$）非线性项 |
| 3 | `Separable` | 一阶，$y$ 的项可与 $y'$ 分组，其余项仅含 $x$ |
| 4 | `Exact` | 一阶，$M + Ny' = 0$，$\partial M/\partial y = \partial N/\partial x$ |
| 5 | `Homogeneous` | 一阶**线性**，所有加法项在 $(x,y)$ 上具有相同总次数（零次齐次），且存在不含 $y$ 的项 |
| 6 | `LinearConstantCoeff` | $\geq 2$ 阶，线性，各阶导数的系数不含自变量 |
| 7 | `CauchyEuler` | $\geq 2$ 阶，$y^{(k)}$ 的系数为 $c_k x^k$ |
| 8 | `ReductionOfOrder` | 二阶线性，非上述类型 |
| 9 | `PowerSeries` | 线性 ODE 的回退选项（先尝试常点幂级数，再尝试 Frobenius） |

`dsolve` 按此顺序依次尝试各方法，第一个成功的即返回结果。用户也可通过 `hint` 参数指定具体方法。

### 求解器架构

```
dsolve(ctx, ode, hint?)
├── normalize_ode()          // 化简方程
├── classify_ode()           // 返回候选方法列表
└── 依次尝试：
    ├── first_order::solve_linear_first()
    ├── first_order::solve_bernoulli()
    ├── first_order::solve_separable()
    ├── first_order::solve_exact()
    ├── first_order::solve_homogeneous()
    ├── second_order::solve_constant_coeff()
    ├── second_order::solve_cauchy_euler()
    ├── second_order::solve_reduction_of_order()
    └── series::solve_power_series() 或 series::solve_frobenius()

dsolve_ivp(ctx, ode, y0, y1?)
└── laplace::solve_laplace()    // Laplace 变换法

dsolve_system(ctx, equations, funcs, var)
└── systems::solve_linear_system()  // 2×2 特征值分解
```

### 返回值类型

`ODESolution` 枚举涵盖所有可能的输出形式：

| 变体 | 含义 | 使用场景 |
|---|---|---|
| `Explicit(expr)` | $y = \text{expr}$ | 大多数成功求解 |
| `Implicit(expr)` | $F(x,y) = 0$ | 可分离方程（隐式解） |
| `Parametric(x_t, y_t)` | $(x(t), y(t))$ | 参数形式解 |
| `Series(expr, n_terms)` | 截断级数 | 幂级数 / Frobenius |
| `System([y1, y2])` | 各分量 | $2 \times 2$ 方程组 |
| `Unsolved(ode)` | 原方程 | 所有方法均失败 |

### 与符号积分的集成

ODE 求解器深度依赖 oCAS 的符号积分能力：

- **一阶线性**：积分因子 $\mu = \exp(\int p\, dx)$ 和最终积分 $\int \mu\, q\, dx$
- **Bernoulli**：线性化后的积分因子
- **恰当方程**：势函数的二重积分
- **齐次方程**：分离后的积分
- **参数变易**：$\int y_i\, g / W\, dx$
- **降阶法**：$\int e^{-\int p\, dx} / y_1^2\, dx$

这些积分调用 oCAS 的分层积分管线（有理函数 → Risch → 三角 → 特殊函数 → 启发式），详见[符号微积分](./symbolic-calculus.md)和 [Risch 积分算法](./risch-algorithm.md)。

### 局限性

当前实现的已知限制：

1. **二阶以上**：仅支持线性常系数和 Cauchy–Euler 型。一般的变系数高阶方程无专用求解器。
2. **非线性 ODE**：除 Bernoulli（可线性化）外，非线性方程通常返回 `Unsolved`。
3. **Laplace 变换**：仅处理整数系数的常系数方程，且强迫项必须有表格变换。
4. **级数解**：默认展开 $x_0 = 0$ 处 8 项。Frobenius 方法仅处理指标方程有有理实根的情况。
5. **2×2 系统**：仅支持常系数，且特征值需为整数或有理数。
6. **验证**：`verify_solution` 将解代入原方程检查残差是否为零（仅用于测试）。

---

## 参考文献

1. **Boyce, W. E. & DiPrima, R. C.** *Elementary Differential Equations and Boundary Value Problems*, 11th ed., Wiley, 2017. — Chapters 2–7 cover first-order equations (Ch.2), second-order linear (Ch.3–4), series solutions (Ch.5), Laplace transform (Ch.6), and systems (Ch.7).

2. **Coddington, E. A. & Levinson, N.** *Theory of Ordinary Differential Equations*, McGraw-Hill, 1955. — Rigorous treatment of existence, uniqueness, and stability.

3. **Ince, E. L.** *Ordinary Differential Equations*, Dover, 1956. — Classical reference on special functions and series methods.

4. **Zill, D. G.** *A First Course in Differential Equations with Modeling Applications*, 12th ed., Cengage, 2021. — Accessible introduction with applications.

5. **Tenenbaum, M. & Pollard, H.** *Ordinary Differential Equations*, Dover, 1985. — Comprehensive classical methods including exact equations and integrating factors.
