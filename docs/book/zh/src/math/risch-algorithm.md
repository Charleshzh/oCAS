# 高阶：Risch 积分算法

## 前提知识

- [多项式代数](./polynomial-algebra.md) — 多项式环、因式分解、结式
- [有限域与模算术](./finite-fields.md) — 模逆与扩展 Euclid 算法
- [线性代数](./linear-algebra.md) — 高斯消元与行列式
- [符号微积分](./symbolic-calculus.md) — 求导规则与链式法则
- [多项式 GCD 与因式分解](./poly-gcd-factoring.md) — 无平方分解与 Hensel 提升

---

## 基础概念

### 初等函数的精确定义

"初等函数"一词在微积分课本中往往凭直觉使用。Risch 算法需要严格定义：

**定义 1**（常数域）. 令 $k$ 为一个特征零的微分域，导子记为 $D$。最基本的常数域是 $\mathbb{Q}$，其上 $D = 0$。

**定义 2**（初等扩张）. 设 $k$ 为微分域，$K$ 为 $k$ 的扩域。称 $K/k$ 是一个**初等扩张**（elementary extension），如果 $K$ 可由以下三种类型的有限扩张逐层生成：

1. **代数扩张**：$K = k(\alpha)$，其中 $\alpha$ 在 $k$ 上代数（即 $\alpha$ 是某个 $k[x]$ 中多项式的根）。例如 $\mathbb{Q}(\sqrt{2})/\mathbb{Q}$。
2. **对数扩张**：$K = k(t)$，其中 $Dt = Du / u$ 对某个 $u \in k^*$，即 $t = \log(u)$。导子规则：$D(\log u) = Du/u$。
3. **指数扩张**：$K = k(t)$，其中 $Dt / t = Du$ 对某个 $u \in k$，即 $t = \exp(u)$。导子规则：$D(\exp u) = \exp(u) \cdot Du$。

**定义 3**（初等函数）. 域 $k$ 上的**初等函数**是某个初等扩张塔 $k = k_0 \subset k_1 \subset \cdots \subset k_n$ 中的元素，其中每个 $k_{i+1}/k_i$ 是上述三种类型之一。

**定理**（Liouville, 1835）. 设 $f$ 为 $k$ 上的初等函数。$f$ 在 $k$ 上有初等原函数，当且仅当它可以写成

$$f = Dv + \sum_{i=1}^{n} c_i \cdot \frac{Du_i}{u_i}$$

其中 $v, u_1, \dots, u_n \in k$，$c_1, \dots, c_n$ 为 $D$-常数（即 $Dc_i = 0$）。

这个定理的深刻意义在于：初等函数的原函数如果存在，其"新"的部分只可能是常数倍的对数。这是 Risch 算法的理论基础。

### 微分域

**定义 4**（微分域）. 一个**微分域**（differential field）是一个域 $k$ 配备一个导子（derivation）$D: k \to k$，满足：

- **加法性**：$D(a + b) = Da + Db$
- **Leibniz 法则**：$D(a \cdot b) = a \cdot Db + Da \cdot b$

导子的核 $\{c \in k : Dc = 0\}$ 称为**常数域**（field of constants），记为 $C_k$ 或 $k^D$。

在 Risch 算法的语境下：
- 基域 $k_0 = \mathbb{Q}(x)$，$D = d/dx$（关于积分变量 $x$ 的导数），常数域 $C_{k_0} = \mathbb{Q}$。
- 对数扩张 $k(t) = k(\log u)$，$Dt = Du/u$。
- 指数扩张 $k(t) = k(\exp u)$，$Dt = t \cdot Du$。

**关键事实**：在初等扩张塔中，常数域保持不变：$C_{k_n} = C_{k_0} = \mathbb{Q}$。这个性质在算法中被隐式使用——当我们说"有理常数"时，就是 $\mathbb{Q}$ 中的元素。

---

## 核心理论

### 微分域塔的构造

给定被积函数 $f$（一个初等函数表达式），Risch 算法的第一步是构造它所"生活"的微分域塔。

**算法** `build_tower(f, var)`：

1. **收集函数符号**：扫描 $f$ 的表达式树，收集所有函数应用（`log`、`exp`）。
2. **排序**：按依赖关系排序——若 $t_2 = \exp(\log(x) + 1)$ 则 $t_1 = \log(x)$ 必须在 $t_2$ 之前。
3. **代数依赖检测**：对每个候选生成元 $t_i$，检查它是否在已有塔上代数相关。保守策略：
   - $\log(c \cdot u)$ 与 $\log(u)$（$c$ 为有理常数）→ 代数相关（$\log(c \cdot u) = \log(u) + \log(c)$）
   - $\exp(u + c)$ 与 $\exp(u)$（$c$ 为有理常数）→ 代数相关（$\exp(u + c) = e^c \cdot \exp(u)$）
   - 否则拒绝（保守策略：宁可拒绝错误相关，也不错误合并）
4. **拒绝非整数幂**：若 $f$ 含有 $\sqrt{x}$ 等非整数幂（即代数函数），拒绝（返回 `None`）。
5. **计算导数**：对每个生成元 $t_i$，计算 $Dt_i$（关于 $D = d/dx$）：
   - $t_i = \log(u)$：$Dt_i = Du / u$（递归计算 $Du$）
   - $t_i = \exp(u)$：$Dt_i = t_i \cdot Du$

在 oCAS 中，`build_tower` 位于 `ocas-calc/src/tower/build.rs`。它返回一个 `Tower` 结构体，包含生成元列表 `gens: Vec<GenInfo>`（每个含 kind、原子引用、导数）和常数域信息。

**域元素的表示**：塔中的元素使用**稀疏多元多项式**（`SparseMultivariatePolynomial<RationalDomain, Lex>`）的分子/分母对表示：

- `KElem`（k-元素）：`{ num: Sparse, den: Sparse }`，表示 $k_\ell$ 中的一个元素 $p/q$
- `KPoly`（k-多项式）：关于顶层生成元 $t_\ell$ 的一元多项式，系数是 `KElem`
- `KRat`（k-有理函数）：`{ num: KPoly, den: KPoly }`，表示 $k_\ell(t_\ell)$ 中的一个有理函数

**例子**. 对 $f = x \cdot \exp(x) + \log(x)$：

| 层级 | 域 | 生成元 | 导子 |
|---|---|---|---|
| $k_0$ | $\mathbb{Q}(x)$ | $x$ | $Dx = 1$ |
| $k_1$ | $k_0(t_1)$，$t_1 = \log(x)$ | $t_1$ | $Dt_1 = 1/x$ |
| $k_2$ | $k_1(t_2)$，$t_2 = \exp(x)$ | $t_2$ | $Dt_2 = t_2$ |

### 逐层积分：总框架

Risch 算法的核心是**自顶向下逐层递归**。在塔的每一层 $\ell$，被积函数 $f$ 被视为 $k_{\ell-1}(t_\ell)$ 中的一个有理函数。

**算法** `integrate_level(tower, level, f)`：

设 $f = a/d \in k_\ell(t_\ell)$，其中 $a, d \in k_\ell[t_\ell]$，$\gcd(a, d) = 1$。

1. **Hermite 约化**（分离有理部分）：将 $a/d$ 分解为

$$\frac{a}{d} = Dg + \frac{a_1}{d_1}$$

其中 $g \in k_\ell(t_\ell)$ 有显式原函数，$d_1$ 无平方。

2. **积分多项式部分**：若 $a/d$ 的多项式部分为 $p(t_\ell) \in k_\ell[t_\ell]$，则按生成元类型积分：
   - **原始层**（$t_\ell = \log u$）：待定系数法
   - **超指数层**（$t_\ell = \exp u$）：Risch 微分方程

3. **积分对数部分**：对 $a_1/d_1$（$d_1$ 无平方），检查是否匹配**对数导数恒等式**。

4. **基域** $\mathbb{Q}(x)$：委托给有理函数积分器。

结果组装为 `LevelResult { elem, logs, extras }`，其中 `elem` 是域元素部分，`logs` 是 $\sum c_i \log(v_i)$ 形式的对数项，`extras` 是未能积分的部分。

### Hermite 约化

**目标**：将有理函数 $a/d$ 分解为"已积"部分加"更简单"的剩余（剩余部分的分母无平方）。

**算法**（Bronstein Ch. 2 的 Hermite 约化；oCAS 采用无平方分解递推形式）：

设 $a/d$ 是一个真分式（$\deg a < \deg d$），$d$ 为 monic 多项式。

1. 计算 $d$ 的无平方分解 $d = \prod_j v_j^{m_j}$。若最大重数 $m \leq 1$，则 $d$ 已无平方，直接返回 $(0, a, d)$。
2. 选取重数最大（$= m \geq 2$）的因子 $v$，写 $d = u \cdot v^m$。
3. 令 $B = u \cdot Dv$（塔导子下的导数；基域上即 $u \cdot v'$），用扩展 Euclid 算法求 $s, t$ 使得 $s \cdot B + t \cdot v = a$（$s$ 取对 $v$ 的模余式）。
4. 递推公式：
$$\frac{a}{u\, v^m} = D\!\left(\frac{-s}{(m-1)\, v^{m-1}}\right) + \frac{t + u \cdot D\bigl(s/(m-1)\bigr)}{u\, v^{m-1}}$$
5. 对剩余项 $\dfrac{t + u \cdot D(s/(m-1))}{u\, v^{m-1}}$ 递归执行，直到分母无平方。

其中 $D\!\left(\frac{-s}{(m-1) v^{m-1}}\right)$ 有显式原函数 $\frac{-s}{(m-1) v^{m-1}}$ 本身，最终剩余 $a_1/d_1$ 的分母无平方。每次迭代将最大重数降低一次，故至多迭代 $m$ 次。

在 oCAS 中，Hermite 约化在两个地方实现：
- `ocas-calc/src/integral/rational.rs`：基域 $\mathbb{Q}(x)$ 上的 Hermite 约化（`hermite_reduce`）
- `ocas-calc/src/integral/risch.rs`：任意塔层上的 Hermite 约化（`hermite_tower`）

**例子**. $\int \frac{1}{(x+1)^2}\, dx$：

$d = (x+1)^2$，其无平方分解为 $(x+1)^2$，最大重数 $m = 2$，取 $v = x+1$，$u = 1$。$B = u \cdot Dv = 1$，解 $s \cdot B + t \cdot v = 1$ 得 $s = 1$（对 $v$ 取模），$t = 0$。递推公式给出：

$$\frac{1}{(x+1)^2} = D\!\left(\frac{-1}{x+1}\right) + \frac{0}{x+1}$$

即 $\int = -1/(x+1)$。更一般地，$\int a/(x+1)^2$ 总能用 Hermite 约化处理。

### 对数导数恒等式

**核心思想**：形如 $c \cdot \frac{Du}{u}$ 的项有平凡原函数 $c \cdot \log(u)$。

**定理**（对数导数恒等式）. 设 $f \in k_\ell(t_\ell)$，$d_1$ 无平方。若

$$\frac{a_1}{d_1} = c \cdot \frac{D d_1}{d_1}$$

对某个 $D$-常数 $c$，则

$$\int \frac{a_1}{d_1} = c \cdot \log(d_1)$$

更一般地，对无平方分母 $d_1$，对数部分的所有贡献形如 $\sum_i c_i \log(v_i)$，其中 $v_i$ 由 $d_1$ 的不可约因子给出。

**判断方法**：在 $\mathbb{Q}(x)$ 基域上，检查 $a_1$ 是否等于某个有理常数 $c$ 乘以 $d_1'$。在高阶塔上，检查 $a_1$ 是否是系数域中某个常数乘以 $d_1'$。

**更精确的方法**（Rothstein–Trager 结式）：对一般无平方分母 $d(x)$，$\int a(x)/d(x)\,dx$ 的对数部分由 **Rothstein–Trager 结式**给出：

$$R(t) = \text{Res}_x\!\bigl(d(x),\; a(x) - t \cdot d'(x)\bigr)$$

设 $R(t)$ 在 $\mathbb{Q}$ 上的根为 $c_1, \dots, c_m$，则

$$\int \frac{a}{d}\,dx = \sum_{i=1}^{m} c_i \cdot \log\gcd(d,\; a - c_i \cdot d') + C$$

在 oCAS 的实现中（`rothstein_trager`），$R(t)$ 通过**插值**计算——在 $t = 0, 1, 2, \dots$ 处求结式值，然后 Lagrange 插值恢复多项式。当 $R(t)$ 在 $\mathbb{Q}$ 上不完全分裂时，对应项以未求值形式 `Integral(term, var)` 返回。

**二次分母的特殊情况**：当 $d(x) = x^2 + bx + c$（不可约）时，配方（completing the square）给出：

$$\int \frac{Ax + B}{x^2 + bx + c}\,dx = \frac{A}{2}\log(x^2 + bx + c) + \frac{2B - Ab}{\sqrt{4c - b^2}}\arctan\!\left(\frac{2x + b}{\sqrt{4c - b^2}}\right)$$

当 $4c - b^2 < 0$（实根情况），$\arctan$ 变为 $\text{artanh}$（反双曲正切），对应对数形式。

### Risch 微分方程

Risch 算法中最核心的子问题是求解 **Risch 微分方程**（Risch Differential Equation, RDE）。

**问题**. 给定微分域 $k_\ell$ 中的元素 $f, g$（两者都不含顶层变量 $t_\ell$），求 $q \in k_\ell[t_\ell]$ 满足：

$$Dq + f \cdot q = g$$

其中 $D$ 是塔导子（对 $t_\ell$ 的全导数）。注意 $f$ 和 $g$ 的系数在 $k_\ell$ 中，$q$ 是 $k_\ell[t_\ell]$ 中的多项式。

**为什么只求多项式解**：有理函数解 $q = p/d$ 需要额外的分母界（denominator bound）分析，这是当前 oCAS 实现未覆盖的片段。返回 `None` 时调用者回退到其他积分方法。

#### 基域 $\mathbb{Q}(x)$ 上的 RDE

在基域 $k_0 = \mathbb{Q}(x)$ 上，RDE 化为常微分方程 $q'(x) + f(x) \cdot q(x) = g(x)$。

**算法**（`base_rde`）：

1. **度界**：设 $\deg f = p$、$\deg g = r$。$f = 0$ 时直接积分，$\deg q = r + 1$；$f$ 为非零常数（$p = 0$）时 $\deg q = r$；$p \geq 1$ 时 $\deg(f \cdot q) = p + \deg q$ 的首项主导 $\deg(q') \leq \deg q - 1$，故 $\deg q = r - p$（唯一确定）。
2. **待定系数**：设 $q = \sum_{i=0}^{n} a_i x^i$，代入方程。
3. **从高次到低次消除**：比较 $x$ 的各次幂系数，逐个确定 $a_i$。每个 $a_i$ 要么被唯一确定，要么矛盾（原函数不存在）。

**例子**. 求解 $q' + q = x$（$f = 1, g = x$）：

设 $q = a_1 x + a_0$。$q' = a_1$。代入：$a_1 + a_1 x + a_0 = x$。

比较系数：$x^1: a_1 = 1$；$x^0: a_1 + a_0 = 0 \Rightarrow a_0 = -1$。

解：$q = x - 1$。验证：$(x-1)' + (x-1) = 1 + x - 1 = x$ ✓

这给出 $\int x \cdot e^x\,dx = (x - 1) e^x$（在超指数层中）。

#### 原始层（Primitive Level）的 RDE

在 $t_\ell = \log(u)$ 的层上，导子满足 $Dt_\ell = Du/u \in k_{\ell-1}$（不含 $t_\ell$）。

**性质**：$D(a_0 + a_1 t + \cdots + a_m t^m) = Da_0 + Da_1 \cdot t + \cdots + Da_m \cdot t^m + (a_1 + 2a_2 t + \cdots + m a_m t^{m-1}) \cdot Dt$

由于 $Dt$ 不含 $t$（原始扩张的关键性质），$D$ 作用在 $k[t]$ 上时**不改变关于 $t$ 的次数**（除了 $\partial/\partial t$ 部分降低一次）。

**算法**：从 $q$ 的最高次系数 $a_m$ 开始，**自顶向下消除**：

1. 由方程中 $t^m$ 的系数确定 $a_m$（一个关于 $a_m$ 的 RDE 在 $k_{\ell-1}$ 中）。
2. 代入后降低次数，对 $a_{m-1}$ 重复。
3. 最终 $a_0$ 需满足一个**对数约束**：$a_0$ 的常数部分必须使下层积分中出现的对数项一致。

**递归结构**：每步产生一个 $k_{\ell-1}$ 上的 RDE，递归调用 `rde_solve` 直到基域。

#### 超指数层（Hyperexponential Level）的 RDE

在 $t_\ell = \exp(u)$ 的层上，导子满足 $Dt_\ell = t_\ell \cdot Du$。

**关键性质**：$Dt = t \cdot Du$ 含 $t$，因此 $D$ 作用在 $k[t]$ 上时**混合次数**：

$$D(a_i t^i) = (Da_i) t^i + a_i \cdot i \cdot t^i \cdot Du = (Da_i + i \cdot a_i \cdot Du) t^i$$

即 $D(a_i t^i) = (Da_i + i \cdot a_i \cdot Du) \cdot t^i$——每个 $t^i$ 层**独立**！

**解耦性质**：将 $q = \sum a_i t^i$ 和 $g = \sum b_j t^j$ 代入 $Dq + fq = g$，比较 $t^k$ 的系数：

$$Da_k + k \cdot a_k \cdot Du + f \cdot a_k = b_k$$

即

$$Da_k + (f + k \cdot Du) \cdot a_k = b_k$$

这是 $k_{\ell-1}$ 上的一个 RDE（未知数为 $a_k$），各 $k$ 独立！

**算法**：对 $k = 0, 1, \dots, \deg g$，分别解 RDE：

$$Da_k + (f + k \cdot Du) \cdot a_k = b_k$$

若某个 $a_k$ 无解，则整体无解。$\deg q$ 的上界由 $\deg g$ 给出。

**例子**. $\int x \cdot e^x\,dx$：在 $t = e^x$（$Dt = t$）的塔上，被积函数为 $x \cdot t$。

多项式部分 $p(t) = x \cdot t$（$m = 1$）。比较 $t^1$：$Da_1 + a_1 \cdot 1 = x$，即 $a_1' + a_1 = x$。由基域 RDE 解得 $a_1 = x - 1$。$t^0$ 部分为 0。

结果：$\int = (x-1) \cdot t = (x-1) \cdot e^x$。

### 对数层的待定系数法

在原始（对数）层 $t = \log(u)$ 上积分多项式 $p(t) = \sum a_i t^i$ 时，使用**待定系数法**。

**算法** `integrate_kpoly_primitive`：

设待积多项式为 $p(t) = \sum_{i=0}^{m} a_i t^i$，猜测原函数 $q(t) = \sum_{i=0}^{m'} b_i t^i$。

1. 设 $m' = m + 1$（原函数的次数至多比被积函数高 1，因为 $D$ 的 $\partial/\partial t$ 部分会降次）。
2. 计算 $Dq = \sum (Db_i) t^i + \sum i \cdot b_i \cdot (Du/u) \cdot t^{i-1}$。
3. 比较 $p = Dq$ 的系数：
   - $t^m$：$a_m = Db_m + (m+1) b_{m+1} \cdot Du/u$
   - $\vdots$
   - $t^0$：$a_0 = Db_0 + b_1 \cdot Du/u$
4. 从 $t^{m+1}$ 开始向下求解。每个 $b_i$ 可能需要递归调用积分器（对 $Db_i$ 部分）。
5. $b_0$ 的确定受**对数约束**限制：在最终结果中，所有对数项 $\sum c_j \log(v_j)$ 必须一致。若 $b_0$ 含有来自下层的对数项且不满足 Liouville 条件，则无解。

### 三角到复指数的重写

Risch 算法原生只处理 `log` 和 `exp`。三角函数需要先重写为复指数形式。

**Euler 公式**：

$$\sin(u) = \frac{e^{iu} - e^{-iu}}{2i}, \qquad \cos(u) = \frac{e^{iu} + e^{-iu}}{2}$$

$$\tan(u) = \frac{\sin u}{\cos u} = \frac{e^{iu} - e^{-iu}}{i(e^{iu} + e^{-iu})}$$

在 oCAS 中，`trig_to_exp`（`ocas-calc/src/integral/trig.rs`）遍历表达式树，将 `sin`、`cos`、`tan`、`cot`、`sec`、`csc` 逐个替换为上述等价形式。虚数单位 $I$ 作为常数生成元加入塔中（$DI = 0$）。

**重写后的积分**：重写后的表达式在 $\mathbb{Q}(x, I, e^{Ix}, \dots)$ 上的微分域塔中，可由 Risch 算法处理。

**回转为实数形式**（`realify`）：Risch 积分的结果可能含有虚数单位。`realify` 尝试将结果转换回实数形式：

- **共轭对数合并**：若结果中含 $c \cdot \log(u + Iv) + c \cdot \log(u - Iv)$，合并为 $c \cdot \log(u^2 + v^2)$
- **共轭对数之差**：若含 $c \cdot \log(u + Iv) - c \cdot \log(u - Iv)$，合并为 $2c \cdot \arctan(v/u)$
- **指数合并**：$e^{Iu} \cdot e^{-Iu} = 1$，$e^{Iu} + e^{-Iu} = 2\cos(u)$ 等

这是"尽力而为"的过程。若无法匹配已知模式，保留复数形式（数学上仍然正确，因为求导验证仍然成立）。

**当前限制**：Risch 微分方程求解器在 $\mathbb{Q}[x]$ 上工作，因此当三角重写产生含 $I$ 的超指数方程（例如 $\sin(x)\cos(x)$ 或 $\cos^2(x)$ 积分时产生的方程）时，系数域扩展为 $\mathbb{Q}(I)$ 而非 $\mathbb{Q}$，求解器可能失败。这些被积函数返回未求值形式。

### 启发式层

在 Risch 算法之前，oCAS 尝试一组**启发式积分技术**（`ocas-calc/src/integral/heuristic.rs`）。这些技术快速且覆盖常见模式。

#### 分部积分（LIATE 排序）

对乘积 $f \cdot g$，尝试分部积分 $\int u\,dv = uv - \int v\,du$。

**LIATE 优先级**（选 $u$ 的启发式规则）：

| 优先级 | 类型 | 分数 | 例子 |
|---|---|---|---|
| 1 | **L**ogarithmic | 0 | $\log x$, $\arctan x$ |
| 2 | **I**nverse trig | 1 | $\arcsin x$, $\text{arccosh}\, x$ |
| 3 | **A**lgebraic | 2 | $x^2$, $\sqrt{x}$ |
| 4 | **T**rigonometric | 3 | $\sin x$, $\cos x$ |
| 5 | **E**xponential | 4 | $e^x$, $2^x$ |

**算法**：对乘积中的每个因子，按 LIATE 分数排序。选择分数最低（优先级最高）的作为 $u$，其余的（经积分得 $v$）作为 $dv$。

**深度限制**：`PARTS_MAX_DEPTH = 2`，防止无限递归。

#### 三角替换

匹配 $\sqrt{a^2 - x^2}$、$\sqrt{a^2 + x^2}$、$\sqrt{x^2 - a^2}$ 及其倒数，直接返回已知原函数：

| 被积函数 | 原函数 |
|---|---|
| $\frac{1}{\sqrt{a^2 - x^2}}$ | $\arcsin(x/a)$ |
| $\frac{1}{\sqrt{a^2 + x^2}}$ | $\text{arcsinh}(x/a)$ |
| $\frac{1}{\sqrt{x^2 - a^2}}$ | $\text{arccosh}(x/a)$ |
| $\sqrt{a^2 - x^2}$ | $\frac{x\sqrt{a^2-x^2} + a^2 \arcsin(x/a)}{2}$ |
| $\sqrt{a^2 + x^2}$ | $\frac{x\sqrt{a^2+x^2} + a^2 \text{arcsinh}(x/a)}{2}$ |
| $\sqrt{x^2 - a^2}$ | $\frac{x\sqrt{x^2-a^2} - a^2 \text{arccosh}(x/a)}{2}$ |

#### Weierstrass 替换

对 $\sin(u)$ 和 $\cos(u)$ 的有理函数（$u$ 是 $x$ 的线性函数），使用万能替换 $t = \tan(u/2)$：

$$\sin u = \frac{2t}{1+t^2}, \qquad \cos u = \frac{1-t^2}{1+t^2}, \qquad du = \frac{2\,dt}{1+t^2}$$

将三角有理积分化为 $t$ 的有理函数积分。

#### Euler 替换

对含 $\sqrt{ax^2 + bx + c}$ 的被积函数，根据 $a$ 的符号选择 Euler 替换：

- $a > 0$：$\sqrt{a}\,x + t = \sqrt{ax^2 + bx + c}$
- $a < 0$：$t(x - \alpha) = \sqrt{ax^2 + bx + c}$（$\alpha$ 为根之一）
- $c > 0$：$\sqrt{c} + xt = \sqrt{ax^2 + bx + c}$

这些替换将根号有理积分化为有理函数积分。

### 特殊函数表

当 Risch 算法证明某个积分**没有初等原函数**时，许多常见情形仍有特殊函数的闭式。oCAS 在 `ocas-calc/src/integral/special.rs` 中直接编码这些标准反导数（定义与 SymPy 一致）：

| 被积函数 | 原函数 | 特殊函数 |
|---|---|---|
| $e^{-x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erf}(x)$ | 误差函数 |
| $e^{x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erfi}(x)$ | 虚误差函数 |
| $e^{cx^2}$（$c < 0$） | $\frac{\sqrt{\pi}}{2\sqrt{-c}}\,\text{erf}(\sqrt{-c}\,x)$ | 误差函数 |
| $e^x / x$ | $\text{Ei}(x)$ | 指数积分 |
| $\sin(x)/x$ | $\text{Si}(x)$ | 正弦积分 |
| $\cos(x)/x$ | $\text{Ci}(x)$ | 余弦积分 |
| $\sinh(x)/x$ | $\text{Shi}(x)$ | 双曲正弦积分 |
| $\cosh(x)/x$ | $\text{Chi}(x)$ | 双曲余弦积分 |
| $\sin(x^2)$ | $\sqrt{\frac{\pi}{2}}\,S\!\left(\sqrt{\frac{2}{\pi}}\,x\right)$ | Fresnel $S$ |
| $\cos(x^2)$ | $\sqrt{\frac{\pi}{2}}\,C\!\left(\sqrt{\frac{2}{\pi}}\,x\right)$ | Fresnel $C$ |

**匹配逻辑**：`special_integrate` 将被积函数分解为乘积因子，尝试匹配以下模式族：

- **erf 族**：$e^{c \cdot x^2}$（$c < 0$），含 $e^{-x^2}$ 特例
- **Ei 族**：$e^x / x$，含 $e^{cx}/x$ 推广
- **Si/Ci/Shi/Chi 族**：三角或双曲函数除以 $x$
- **Fresnel 族**：$\sin(x^2)$、$\cos(x^2)$

注意：这些**不是** Risch 算法的一部分。Risch 证明"不存在初等原函数"后，特殊函数表提供"次优"答案。

---

## 在 oCAS 中的实现

### 积分管线架构

`integrate(expr, var)`（`ocas-calc/src/integral/mod.rs`）按顺序尝试以下层，第一个产生答案的层获胜：

```
┌─────────────────────────────────────────────────┐
│  integrate(expr, var)                           │
├─────────────────────────────────────────────────┤
│  1. 查找表（幂法则、线性参数函数）               │
│     ↓ 失败                                      │
│  2. 有理函数积分器（integrate_rational）         │
│     ↓ 失败                                      │
│  3. Risch 算法（risch_integrate）                │
│     ↓ 失败                                      │
│  4. 三角重写（trig_to_exp）→ 重试 Risch → realify │
│     ↓ 失败                                      │
│  5. 特殊函数表（special_integrate）              │
│     ↓ 失败                                      │
│  6. 启发式层（heuristic_integrate）              │
│     ↓ 失败                                      │
│  7. 返回未求值形式 Integral(expr, var)           │
└─────────────────────────────────────────────────┘
```

**深度限制**：`MAX_DEPTH = 8`（`integrate_raw` 的递归深度），`MAX_RISCH_DEPTH = 16`（`risch_integrate` 的递归深度，使用 `thread_local!` 计数器）。这些限制防止病态输入（如 `sec(x)` 触发 VOP 无限重试）导致栈溢出。

### 模块组织

```
ocas-calc/src/integral/
├── mod.rs         ← integrate 入口 + 查找表 + 管线调度
├── rational.rs    ← 基域 Q(x) 上的有理函数积分器
│                    Hermite 约化 + 对数部分 + Rothstein–Trager
├── risch.rs       ← Risch 算法主循环
│                    逐层 integrate_level + Hermite 层上约化
│                    多项式部分（原始/超指数）
├── rde.rs         ← Risch 微分方程求解器
│                    base_rde (Q[x]) + 递归塔
├── trig.rs        ← 三角↔复指数重写（trig_to_exp / realify）
├── heuristic.rs   ← 启发式：分部积分 / 三角替换 / Weierstrass / Euler
└── special.rs     ← 特殊函数表（erf, Ei, Si, Ci, Fresnel 等）
```

### Risch 算法的执行流程

以 $\int (x+1) e^x\,dx$ 为例，追踪完整执行：

1. **构建塔**：`build_tower` 收集 $t_1 = \exp(x)$，塔为 $\mathbb{Q}(x) \subset \mathbb{Q}(x)(t_1)$。
2. **`risch_integrate`**：将被积函数转换为 $k_1(t_1)$ 中的有理函数：$(x+1) \cdot t_1$。
3. **`integrate_level(tower, level=1, f)`**：
   - 分母为 1，无 Hermite 约化。
   - 多项式部分 $p(t_1) = (x+1) t_1$。
   - 调用 `integrate_kpoly_hyperexp`（因为 $t_1 = \exp(x)$）：
     - 比较 $t_1^1$ 系数：$Da_1 + a_1 \cdot 1 = x + 1$
     - 递归调用 `rde_solve`：$a_1' + a_1 = x + 1$
     - 基域 `base_rde`：设 $a_1 = ax + b$，$a_1' = a$，$a + ax + b = x + 1$
     - 比较：$a = 1$，$a + b = 1 \Rightarrow b = 0$
     - $a_1 = x$，验证：$x' + x = 1 + x = x + 1$ ✓
   - $t_1^0$ 系数为 0，$b_0 = 0$。
4. **结果组装**：$q = x \cdot t_1$，转为原子：$x \cdot e^x$。
5. **化简 + 归一化**：返回 $x \cdot e^x$。

**验证**：$\frac{d}{dx}[x \cdot e^x] = e^x + x \cdot e^x = (x+1)e^x$ ✓

### 域元素的内部表示

在塔中，域元素使用**平坦多元稀疏多项式**表示：

```rust
// KElem: k = Q(x, t_1, ..., t_{n-1}) 中的元素
struct KElem {
    num: SparseMultivariatePolynomial<RationalDomain, Lex>,
    den: SparseMultivariatePolynomial<RationalDomain, Lex>,
}
```

变量索引分配：`0` 对应 $x$，`1` 对应 $t_1$，…，`n-1` 对应 $t_{n-1}$。

`KPoly` 是关于顶层变量的稠密一元多项式，系数为 `KElem`：

```rust
// KPoly: k[t_n] 中的多项式，系数在 k 中
struct KPoly {
    coeffs: Vec<KElem>,  // 从低次到高次
    top: usize,          // 顶层变量的索引
    n: usize,            // 总变量数
}
```

`KRat` 是有理函数 $p/q$：

```rust
struct KRat {
    num: KPoly,
    den: KPoly,
}
```

**无 GCD 约简**：`KElem` 不执行多元 GCD 约简（因为通用多元 GCD 代价高昂）。零检测仅通过分子为零来判断，交叉相乘相等是唯一可靠的相等性判断。

### Fuel 受限积分

`integrate_with_fuel(ctx, expr, var, &fuel)` 将 `Fuel` 预算贯穿到两个积分后化简阶段。积分遍历本身使用 `MAX_DEPTH` / `MAX_RISCH_DEPTH` 限制；`Fuel` 仅约束化简阶段，防止病态结果（如产生极大表达式的积分）导致重写器无限循环。

```rust
use ocas_core::fuel::Fuel;

let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, expr, Symbol::new("x"), &fuel);
// Ok(result)  — 正常完成
// Err(_)      — 化简中途 fuel 耗尽
```

### 范围限制与回退策略

Risch 算法的当前实现有以下限制：

| 限制 | 原因 | 回退行为 |
|---|---|---|
| 只求 RDE 的多项式解 | 有理函数解需要分母界分析 | 返回 `None`，调用者尝试其他层 |
| 塔层上的对数部分只用对数导数恒等式 $a_1 = c \cdot Dd_1$ | 完整的对数部分需要塔层上的 Rothstein–Trager / 迹函数判定 | 返回未求值形式（基域 $\mathbb{Q}(x)$ 上仍使用完整处理） |
| 代数函数（$\sqrt{x}$ 等）不支持 | 需要代数函数域扩展 | 由启发式三角替换覆盖常见模式 |
| 代数相关生成元保守拒绝 | $\log(2x)$ 与 $\log(x)$ 的关系检测 | 返回 `None` |
| 含 $I$ 的超指数 RDE | RDE 求解器仅在 $\mathbb{Q}[x]$ 上工作 | 三角被积函数返回未求值形式 |

当所有层都失败时，返回 `Integral(expr, var)`——这是**有意的答案**，表示"该积分在当前实现中没有闭式解"，而非程序错误。

---

## 参考文献

1. **Bronstein, M.** *Symbolic Integration I: Transcendental Functions*, 2nd ed., Springer, 2005. — Risch 算法的权威参考。第 2 章：有理函数积分（Hermite 约化、Rothstein–Trager）；第 5 章：初等函数上的 Risch 算法；第 6 章：Risch 微分方程。oCAS 的实现严格遵循此书。
2. **Bronstein, M.** *Symbolic Integration II: Transcendental and Algebraic Functions*, Springer, 2004. — 扩展到代数函数和特殊函数。oCAS 尚未覆盖此书的内容。
3. **Geddes, K. O., Czapor, S. R., & Labahn, G.** *Algorithms for Computer Algebra*, Kluwer, 1992. — 第 12 章提供 Risch 算法的替代讲解。
4. **Liouville, J.** "Sur les transcendantes elliptiques de première et de seconde espèce considérées comme fonctions de leur amplitude." *Journal de l'École Polytechnique*, 1835. — 初等函数原函数的 Liouville 定理原始出处。
5. **Risch, R. H.** "The problem of integration in finite terms." *Transactions of the AMS*, 139:167–189, 1969. — Risch 算法的奠基论文。
6. **Rothstein, M.** "A new algorithm for the integration of exponential and logarithmic functions." *Proceedings of the 1977 MACSYMA Users Conference*, 1977. — Rothstein–Trager 结式方法。
7. **Trager, B. M.** "Algebraic factoring and rational function integration." *Proceedings of SYMSAC '76*, 1976. — Trager 算法（代数数域上的因式分解）。
8. **Lazard, D. & Rioboo, R.** "Integration of rational functions: Rational computation of the logarithmic part." *Journal of Symbolic Computation*, 9(2):113–129, 1990. — 对数部分的高效算法。
