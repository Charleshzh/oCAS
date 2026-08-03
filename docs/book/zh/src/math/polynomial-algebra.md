# 基础：多项式代数

多项式代数是符号计算的基石。本章从多项式环的基本定义出发，逐步引入单项式序、多元除法算法和 Hilbert 基定理，并说明 oCAS 如何在 Rust 中实现这些概念。

## 前提知识

阅读本章前，建议具备以下基础：

- **抽象代数**：群、环、域的基本概念（单位元、零因子、整环）
- **线性代数**：向量空间、线性无关、基
- **数学归纳法**：用于证明除法算法的终止性

如需复习，可参考 [数学基础总览](./overview.md) 中的学习路径推荐。

## 基础概念

### 一元多项式环

设 $\mathbb{F}$ 为一个域（如 $\mathbb{Q}$、$\mathbb{F}_p$）。**一元多项式环** $\mathbb{F}[x]$ 是所有形如

$$f(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n$$

的表达式构成的集合，其中 $a_i \in \mathbb{F}$，$n \geq 0$。

**定义（次数）**。对于非零多项式 $f = \sum_{i=0}^n a_i x^i$（$a_n \neq 0$），其**次数**（degree）为

$$\deg(f) = n$$

零多项式 $f = 0$ 的次数约定为 $-\infty$（或在实现中表示为 `None`）。

**定义（首项系数）**。$f$ 的**首项系数**（leading coefficient）为 $\text{lc}(f) = a_n$。当 $\text{lc}(f) = 1$ 时，称 $f$ 为**首一多项式**（monic）。

**性质**：

| 运算 | 次数关系 |
|---|---|
| $f + g$ | $\deg(f+g) \leq \max(\deg f, \deg g)$ |
| $f \cdot g$ | $\deg(f \cdot g) = \deg f + \deg g$ |
| $f / g$（带余除法） | $\deg(r) < \deg(g)$ |

### 多元多项式环

设 $\mathbb{F}$ 为域，$x_1, x_2, \ldots, x_n$ 为 $n$ 个不定元。**多元多项式环** $\mathbb{F}[x_1, \ldots, x_n]$ 是所有有限和

$$f = \sum_{\alpha \in \mathbb{N}^n} c_\alpha \, \mathbf{x}^\alpha$$

构成的集合，其中 $\alpha = (\alpha_1, \ldots, \alpha_n) \in \mathbb{N}^n$ 称为**指数向量**，$\mathbf{x}^\alpha = x_1^{\alpha_1} x_2^{\alpha_2} \cdots x_n^{\alpha_n}$ 称为**单项式**（monomial），$c_\alpha \in \mathbb{F}$ 称为**系数**。

**定义（全次数）**。单项式 $\mathbf{x}^\alpha$ 的**全次数**（total degree）为 $|\alpha| = \alpha_1 + \cdots + \alpha_n$。多项式 $f$ 的全次数为 $\deg(f) = \max\{|\alpha| : c_\alpha \neq 0\}$。

**定义（偏次数）**。多项式 $f$ 关于变量 $x_i$ 的**偏次数**为 $\deg_{x_i}(f) = \max\{\alpha_i : c_\alpha \neq 0\}$。

> **重要区别**：一元多项式环 $\mathbb{F}[x]$ 是主理想整环（PID），每个理想都由单个多项式生成。多元多项式环 $\mathbb{F}[x_1, \ldots, x_n]$（$n \geq 2$）不再是 PID，这是 Gröbner 基理论的出发点——参见 [Gröbner 基理论](./groebner-theory.md)。

### 多项式的内部表示

多项式有两种基本的计算机表示方式：

| 表示 | 存储方式 | 适用场景 |
|---|---|---|
| **稠密表示** | 系数向量 $[a_0, a_1, \ldots, a_n]$ | 一元、中等次数 |
| **稀疏表示** | 仅存储非零项的 $\{(\alpha, c_\alpha)\}$ 对 | 多元、大量项为零 |

在 oCAS 中，`DenseUnivariatePolynomial` 使用稠密表示，`SparseMultivariatePolynomial` 使用稀疏表示。详见 [§ 在 oCAS 中的实现](#在-ocas-中的实现)。

## 单项式序

在多元多项式环中，"谁是首项"不再一目了然——我们需要一种比较单项式的规则。

**定义（单项式序）**。$\mathbb{N}^n$ 上的**单项式序**（monomial ordering）是一个全序 $\succ$，满足：

1. **良序性**：$\mathbb{N}^n$ 的每个非空子集都有最小元（等价地：$0 \preceq \alpha$ 对所有 $\alpha$ 成立）
2. **相容性**：若 $\alpha \succ \beta$，则 $\alpha + \gamma \succ \beta + \gamma$ 对所有 $\gamma \in \mathbb{N}^n$ 成立

这两个条件保证了多项式除法的终止性。下面介绍三种最常用的单项式序。

### 字典序（Lex）

$x^\alpha \succ_{\text{lex}} x^\beta$ 当且仅当 $\alpha - \beta$ 的最左非零分量为正。

等价描述：从左到右逐分量比较指数向量，首次不等处决定大小。

**示例**：在 $\mathbb{F}[x, y, z]$ 中，按字典序排列：

$$x^2 y^3 z \succ_{\text{lex}} x^2 y^2 z^3 \succ_{\text{lex}} x y^5 \succ_{\text{lex}} y^{100} \succ_{\text{lex}} z^{100}$$

字典序的特点是"消除性"：若 $x_1 \succ x_2 \succ \cdots \succ x_n$，则按字典序的首项关于 $x_1$ 的次数最大。这使得字典序基天然具有消元性质。

### 分次字典序（Grlex）

$x^\alpha \succ_{\text{grlex}} x^\beta$ 当且仅当：

1. $|\alpha| > |\beta|$，或
2. $|\alpha| = |\beta|$ 且 $\alpha \succ_{\text{lex}} \beta$

即先按全次数排，同次数内按字典序排。

**示例**：在 $\mathbb{F}[x, y]$ 中，$\deg(x^2) = 2$、$\deg(xy) = 2$、$\deg(y^3) = 3$，故

$$y^3 \succ_{\text{grlex}} x^2 \succ_{\text{grlex}} xy$$

### 分次反字典序（Grevlex）

$x^\alpha \succ_{\text{grevlex}} x^\beta$ 当且仅当：

1. $|\alpha| > |\beta|$，或
2. $|\alpha| = |\beta|$ 且 $\alpha - \beta$ 的**最右**非零分量为**负**

即先按全次数排，同次数内从右往左比较，**小**的一方更大。

**示例**：在 $\mathbb{F}[x, y, z]$ 中，$x^2 y$ 和 $x y^2$ 的全次数都是 3。$\alpha - \beta = (1, -1, 0)$，最右非零分量是 $-1$（负），故

$$x^2 y \succ_{\text{grevlex}} x y^2$$

Grevlex 的关键优势在于计算：对相同全次数的单项式，Grevlex 倾向于让低次变量的指数更小，这在 Gröbner 基计算中产生更小的中间多项式。

### 三种序的比较

| 单项式 | Lex | Grlex | Grevlex |
|---|---|---|---|
| $x^2$ ($\alpha = (2,0)$) | 1 | 2 | 2 |
| $xy$ ($\alpha = (1,1)$) | 2 | 3 | 3 |
| $y^3$ ($\alpha = (0,3)$) | 3 | **1**（次数最高） | **1**（次数最高） |

> **注意**：上表中排名 1 表示在对应序下"最大"。Grlex 与 Grevlex 都先按
> 全次数排序，因此 $y^3$（全次数 3）排在 $x^2$、$xy$（全次数 2）之前。
> 在二元情形下，同次数内的平局决胜规则使 Grlex 与 Grevlex **重合**：
> 对 $x^2$ 与 $xy$，$(2,0) - (1,1) = (1, -1)$ 的最右非零分量为 $-1 < 0$，
> 故两种序下都有 $x^2 \succ xy$。两者的差异要到 $n \geq 3$ 个变量时才
> 显现，例如 $xz^2$ 与 $y^3$（全次数均为 3）：Grlex 下 $xz^2 \succ y^3$
> （字典序平局决胜：$(1,0,2) > (0,3,0)$），而 Grevlex 下
> $(0,3,0) - (1,0,2) = (-1,3,-2)$ 的最右非零分量为 $-2 < 0$，故
> $y^3 \succ_{\text{grevlex}} xz^2$。
>
> 再验证 Grevlex 的平局决胜：$x^2 y$ 与 $x y^2$：$(2,1) - (1,2) = (1, -1)$，
> 最右非零分量是 $-1 < 0$，所以 $x^2 y \succ_{\text{grevlex}} xy^2$。

### 扩展序

除了上述三种基本序，oCAS 还支持以下参数化序：

**加权序**（WeightOrder）：给每个变量赋予权重 $w_i$，按 $\sum w_i \alpha_i$ 降序排列。当权重向量为 $(1, 1, \ldots, 1)$ 时退化为全次数序。

**矩阵序**（MatrixOrder）：给定 $n \times n$ 整数矩阵 $M$，$\alpha \succ \beta$ 当且仅当 $M\alpha >_{\text{lex}} M\beta$。矩阵序是所有线性单项式序的完全一般化。特别地，消元序（elimination order）可通过 `MatrixOrder::elimination_order(elim_vars, n_vars)` 构造，用于从多项式系统中消除变量。

**分块序**（BlockOrder）：将变量划分为连续块，每块内使用独立的子序（Lex、Grlex 或 Grevlex）。消元序的标准形式是"前 $k$ 个变量用 Lex，其余用 Grevlex"。

## 核心理论

### 一元带余除法

**定理（带余除法）**。设 $\mathbb{F}$ 为域，$f, g \in \mathbb{F}[x]$，$g \neq 0$。存在唯一的 $q, r \in \mathbb{F}[x]$ 使得

$$f = q \cdot g + r, \quad \deg(r) < \deg(g)$$

**算法**（逐次消首项）：

```
输入: f, g (g ≠ 0)
输出: 商 q, 余数 r

q ← 0, r ← f
while r ≠ 0 and deg(r) ≥ deg(g):
    t ← lc(r) / lc(g) · x^(deg(r) - deg(g))
    q ← q + t
    r ← r - t · g
return (q, r)
```

**终止性**：每一步 $\deg(r)$ 严格递减（或 $r$ 变为零），故算法必然终止。

> **整环上的推广**：当系数环是欧几里得整环（如 $\mathbb{Z}$）时，带余除法仍然存在，但需要伪除法（pseudo-division）来避免分数。oCAS 中 `DenseUnivariatePolynomial::div_rem` 通过 `EuclideanDomain::div_rem` 对首项系数执行带余除法来实现。

### 多元带余除法

多元情况远比一元复杂——**商和余数不唯一**，依赖于单项式序和除法顺序。

**定理（多元除法算法）**。设 $\mathbb{F}$ 为域，$f, g_1, \ldots, g_s \in \mathbb{F}[x_1, \ldots, x_n]$，$g_i \neq 0$。固定一个单项式序。存在算法计算 $q_1, \ldots, q_s, r \in \mathbb{F}[x_1, \ldots, x_n]$ 使得

$$f = q_1 g_1 + q_2 g_2 + \cdots + q_s g_s + r$$

其中 $r$ 的每个非零项都不被任何 $\text{lt}(g_i)$ 整除。这样的 $r$ 称为 $f$ 对 $\{g_1, \ldots, g_s\}$ 的**余数**（或**正规形**，normal form）。

**算法**：

```
输入: f, G = [g1, ..., gs]
输出: q1, ..., qs, r

r ← 0, p ← f
while p ≠ 0:
    if ∃ i 使得 lt(gi) | lt(p):
        选择最小的这样的 i
        qi ← qi + lt(p) / lt(gi)
        p ← p - (lt(p) / lt(gi)) · gi
    else:
        r ← r + lt(p)
        p ← p - lt(p)
return (q1, ..., qs, r)
```

**关键性质**：

1. **终止性**：由单项式序的良序性保证
2. **不唯一性**：余数 $r$ 依赖于除法顺序 $g_1, \ldots, g_s$。不同的顺序可能产生不同的余数
3. **Gröbner 基的关键作用**：当 $\{g_1, \ldots, g_s\}$ 构成 Gröbner 基时，余数 $r$ 是**唯一**的——这正是 Gröbner 基的核心价值之一

> **示例**：取 $f = x^2 y + x y^2 + y^2$，$g_1 = xy - 1$，$g_2 = y^2 - 1$，使用 Grevlex 序。
>
> 第一步：$\text{lt}(f) = x^2 y$，$\text{lt}(g_1) = xy$ 整除 $x^2 y$。
>
> $$f \leftarrow f - x \cdot g_1 = x^2 y + x y^2 + y^2 - x(xy - 1) = x y^2 + y^2 + x$$
>
> 第二步：$\text{lt}(p) = x y^2$，$\text{lt}(g_1) = xy$ 整除 $x y^2$。
>
> $$p \leftarrow x y^2 + y^2 + x - y \cdot g_1 = x y^2 + y^2 + x - y(xy - 1) = y^2 + x + y$$
>
> 第三步：$\text{lt}(p) = y^2$，$\text{lt}(g_2) = y^2$ 整除 $y^2$。
>
> $$p \leftarrow y^2 + x + y - 1 \cdot g_2 = y^2 + x + y - (y^2 - 1) = x + y + 1$$
>
> 此时 $\text{lt}(p) = x$ 不被 $xy$ 或 $y^2$ 整除，移入余数：$r = x + y + 1$。

### S-多项式

S-多项式是 Gröbner 基算法的核心构件。

**定义**。设 $f, g$ 为非零多项式，$\gamma = \text{lcm}(\text{lm}(f), \text{lm}(g))$（指数向量逐分量取最大值）。$f$ 和 $g$ 的 **S-多项式**为

$$S(f, g) = \frac{\gamma}{\text{lt}(f)} \cdot f - \frac{\gamma}{\text{lt}(g)} \cdot g = \frac{x^\gamma}{\text{lm}(f)} \cdot f - \frac{x^\gamma}{\text{lm}(g)} \cdot g$$

S-多项式的直觉是"消除 $f$ 和 $g$ 的首项"。Buchberger 算法反复构造 S-多项式并化简，直到所有 S-多项式化简为零——此时的多项式集合即为 Gröbner 基。

### Dickson 引理

Dickson 引理是单项式理想的有限生成性定理，也是 Hilbert 基定理的基础。

**引理（Dickson）**。$\mathbb{N}^n$ 的每个子集 $S$ 都有有限的"最小元"集合 $M \subseteq S$，使得对每个 $\alpha \in S$，存在 $\beta \in M$ 满足 $\beta \leq \alpha$（逐分量）。

等价表述：$\mathbb{F}[x_1, \ldots, x_n]$ 中由单项式生成的每个理想都由有限个单项式生成。

**直觉**：在 $\mathbb{N}^n$ 中，"向下封闭"的集合只能由有限个"极小元"生成。这与 $\mathbb{N}$ 中的良序性类似，但推广到高维。

### Hilbert 基定理

**定理（Hilbert 基定理）**。$\mathbb{F}[x_1, \ldots, x_n]$ 中的每个理想 $I$ 都是有限生成的。即存在 $f_1, \ldots, f_s \in I$ 使得 $I = \langle f_1, \ldots, f_s \rangle$。

**证明思路**：设 $J = \langle f_1, \ldots, f_s \rangle$ 为待证对象。考虑首项理想 $L(I) = \langle \text{lm}(f) : f \in I \rangle$。由 Dickson 引理，$L(I) = \langle \text{lm}(f_1), \ldots, \text{lm}(f_s) \rangle$，故 $L(J) \supseteq L(I)$。若 $J \subsetneq I$，取 $f \in I \setminus J$ 使 $\text{lm}(f)$ 在 $\succ$ 下极小。由于 $\text{lm}(f) \in L(I) \subseteq L(J)$，单项式理想的性质保证存在 $f_i$ 使 $\text{lm}(f_i) \mid \text{lm}(f)$；用 $f_i$ 消去 $f$ 的首项得到 $f' = f - c\, \mathbf{x}^\gamma f_i \in I \setminus J$，且 $\text{lm}(f') \prec \text{lm}(f)$——与 $\text{lm}(f)$ 的极小性矛盾（单项式序的良序性禁止无限下降链）。故 $J = I$。

> **意义**：Hilbert 基定理保证了 Gröbner 基的存在性——每个理想的 Gröbner 基都是有限集。Buchberger 算法进一步提供了构造性的证明。

## 在 oCAS 中的实现

### 稠密一元多项式

`DenseUnivariatePolynomial<D>` 将系数从常数项到首项系数存储在一个连续向量中：

```rust
pub struct DenseUnivariatePolynomial<D: Domain> {
    coeffs: Vec<D::Element>,  // [a0, a1, ..., an]
    domain: D,
}
```

**设计决策**：

- **尾零裁剪**：构造时自动移除尾部零系数，零多项式用空向量表示
- **系数域内联**：`domain` 字段存储在多项式内部，避免每次运算时传入
- **多项式类型参数**：`D: Domain` 泛型允许在整数环 $\mathbb{Z}$、有理数域 $\mathbb{Q}$、有限域 $\mathbb{F}_p$ 等上统一实现

**乘法策略**：

| 条件 | 算法 | 复杂度 |
|---|---|---|
| $\min(\deg f, \deg g) + 1 < 32$（系数个数） | Schoolbook（逐项相乘） | $O(n \cdot m)$ |
| $\min(\deg f, \deg g) + 1 \geq 32$ | Karatsuba | $O(n^{1.585})$ |
| 有限域 + NTT 友好素数 + 系数个数 $\geq 256$（`mul_ntt`） | 数论变换 | $O(n \log n)$ |

> **注意**：通用 `mul` 走前两行（Karatsuba/Schoolbook，阈值按系数个数
> `KARATSUBA_THRESHOLD = 32` 判断，即次数 $\geq 31$）。NTT 通过
> `DenseUnivariatePolynomial<FiniteField>::mul_ntt` 启用：它内部检测
> `NTT_THRESHOLD = 256`（系数个数）与素数是否 NTT-friendly，否则回退到
> Karatsuba/Schoolbook。

Karatsuba 的核心思想是将多项式在中点 $m = n/2$ 处拆分：

$$a = a_0 + a_1 x^m, \quad b = b_0 + b_1 x^m$$

通过三次半规模乘法（而非四次）计算乘积：

$$z_0 = a_0 b_0, \quad z_2 = a_1 b_1, \quad z_1 = (a_0 + a_1)(b_0 + b_1) - z_0 - z_2$$

$$a \cdot b = z_0 + z_1 x^m + z_2 x^{2m}$$

**EuclideanDomain 方法**（当 $D: \text{EuclideanDomain}$ 时可用）：

- `div_rem(&self, divisor)` — 带余除法，返回 `Option<(quotient, remainder)>`
- `gcd(&self, other)` — 多项式 GCD（通过伪余式 PRS 链）
- `extended_gcd_poly(&self, other)` — 扩展 GCD，返回 $(g, s, t)$ 满足 $s \cdot f + t \cdot g = \text{gcd}(f, g)$
- `content()` — 内容（所有系数的 GCD）
- `primitive_part()` — 本原部分 $f / \text{content}(f)$
- `square_free_factorization()` — 无平方因式分解
- `resultant(&self, other)` — 结式（Brown PRS 算法）
- `pow(n)` — 快速幂（重复平方法）
- `p_adic_expansion(p)` — $p$-进展开

**一元多项式求值**使用 Horner 法则，$O(n)$ 时间：

$$f(x) = a_0 + x(a_1 + x(a_2 + \cdots + x \cdot a_n))$$

```rust
// Horner's method
let mut result = domain.zero();
for coeff in self.coeffs.iter().rev() {
    result = domain.mul(&result, x);
    result = domain.add(&result, coeff);
}
```

### 稀疏多元多项式

`SparseMultivariatePolynomial<D, O>` 仅存储非零项，用 HashMap 从指数向量映射到系数：

```rust
pub struct SparseMultivariatePolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    terms: HashMap<SmallVec<[usize; 4]>, D::Element>,
    domain: D,
    n_vars: usize,
    pub order: O,
}
```

**设计要点**：

- **指数向量**：`SmallVec<[usize; 4]>`——小变量数时内联存储（≤4 个变量不分配堆内存），大变量数时溢出到堆
- **默认序**：`O = Grevlex`——这是 Gröbner 基计算最常用的序，因为它产生最小的中间表达式
- **序作为类型参数**：单项式序在编译时确定，零开销抽象。`Lex`、`Grlex`、`Grevlex` 是零大小类型（ZST），不占运行时空间；`WeightOrder`、`MatrixOrder`、`BlockOrder` 携带运行时配置

**MonomialOrder trait**：

```rust
pub trait MonomialOrder: Clone + PartialEq + Eq + Debug + Default {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering;
}
```

单一方法 `cmp` 比较两个指数向量。各序的实现：

```rust
// 字典序：直接比较指数向量
impl MonomialOrder for Lex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> Ordering {
        lhs.cmp(rhs)
    }
}

// 分次反字典序：先比全次数（降序），再反字典序
impl MonomialOrder for Grevlex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> Ordering {
        let deg_lhs: usize = lhs.iter().sum();
        let deg_rhs: usize = rhs.iter().sum();
        deg_rhs.cmp(&deg_lhs)
            .then_with(|| rhs.iter().rev().cmp(lhs.iter().rev()))
    }
}

// 分次字典序：先比全次数（降序），再字典序
impl MonomialOrder for Grlex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> Ordering {
        let deg_lhs: usize = lhs.iter().sum();
        let deg_rhs: usize = rhs.iter().sum();
        deg_rhs.cmp(&deg_lhs)
            .then_with(|| lhs.cmp(rhs))
    }
}
```

> **注意实现细节**：`deg_rhs.cmp(&deg_lhs)` 而非 `deg_lhs.cmp(&deg_rhs)`——因为 Rust 的 `Ordering` 中 `Less` 表示"排在前面"，我们需要更高次数排在前面（即 `Less`），所以用 `deg_rhs` 对 `deg_lhs` 比较。类似地，Grevlex 中 `rhs.iter().rev().cmp(lhs.iter().rev())` 实现了"反"字典序。

**多元除法**通过 `reduce(&self, basis: &[Self])` 实现：

```rust
pub fn reduce(&self, basis: &[Self]) -> Self {
    let mut remainder = self.clone();
    let mut result = self.zero();

    // 缓存每个基元素的首项
    let basis_lts: Vec<_> = basis
        .iter()
        .filter_map(|g| g.leading_term().map(|(e, c)| (g, e.clone(), c.clone())))
        .collect();

    let max_iter = 10000;
    for _ in 0..max_iter {
        if remainder.is_zero() {
            break;
        }
        let (rm, rc) = match remainder.leading_term() {
            Some((e, c)) => (e, c),
            None => break,
        };

        let mut reduced = false;
        for (g, lm, lc) in &basis_lts {
            // lm | rm（逐分量 rm[i] ≥ lm[i]）
            if monomial_divides(&rm, lm) {
                // 计算商单项式
                let qm: SmallVec<[usize; 4]> =
                    rm.iter().zip(lm.iter()).map(|(a, b)| a - b).collect();
                let qc = match self.domain.div(&rc, lc) {
                    Some(q) => q,
                    None => break,
                };
                remainder = remainder.sub(&g.mul_monomial(&qm).mul_scalar(&qc));
                reduced = true;
                break;
            }
        }

        if !reduced {
            result.terms.insert(rm, rc);
            remainder.terms.remove(&rm);
        }
    }
    result
}
```

**辅助函数**（单项式操作）：

```rust
// 返回 x^b 是否整除 x^a（逐分量 a[i] ≥ b[i]）
pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool

// 最小公倍单项式（逐分量取 max）
pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]>

// 两个单项式是否互素（无共同变量）
pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool
```

### 稠密 vs 稀疏：选择指南

| 考虑因素 | 稠密 `DenseUnivariatePolynomial` | 稀疏 `SparseMultivariatePolynomial` |
|---|---|---|
| 变量数 | 一元 | 多元（$n$ 个变量） |
| 存储 | 连续 Vec，$O(\deg f + 1)$ | HashMap，$O(\text{非零项数})$ |
| 乘法 | Karatsuba / NTT | 逐项展开 |
| 除法 | 唯一商和余数 | 依赖单项式序，余数不唯一 |
| GCD | Euclid 算法 / PRS | 需要 Gröbner 基技术 |
| 适用系数域 | `EuclideanDomain` | `Domain`（除法需要域） |

**经验法则**：一元多项式用稠密表示；两个以上变量用稀疏表示。即使多元多项式只有少数非零项（如 $x^{100} + y^{100}$），稀疏表示也只需存储两项，而稠密表示需要 $100 \times 100$ 的系数矩阵。

## 进阶话题

### 单项式序的等价与统一

所有标准单项式序都可以用矩阵序统一表示。设 $n$ 个变量：

- **Lex** $=$ $M = I_n$（单位矩阵）
- **Grlex** $=$ 第一行全 1，其余行为字典序单位向量（$e_1, e_2, \ldots, e_n$，按序排列——先比全次数，再从左到右比较各分量）
- **Grevlex** $=$ 第一行全 1，其余行为反字典序的负单位向量（$-e_n, -e_{n-1}, \ldots, -e_1$——先比全次数，再从右到左、取较小者更大）

`MatrixOrder::elimination_order(k, n)` 构造的消元序等价于 `BlockOrder` 中前 $k$ 个变量用 Lex、其余用 Grevlex，但表达为矩阵形式更为简洁。

### 除法的局限性

多元带余除法的一个重要局限：余数依赖于除子的顺序。例如，对 $f$ 用 $(g_1, g_2)$ 除和用 $(g_2, g_1)$ 除可能得到不同的余数。Gröbner 基的引入正是为了解决这个问题——对 Gröbner 基的除法余数是唯一的。

### 整数环上的多项式

当系数环是 $\mathbb{Z}$（而非域）时，带余除法需要调整。`DenseUnivariatePolynomial<IntegerDomain>` 的 `div_rem` 通过在首项系数上执行整数带余除法来实现伪除法（pseudo-division）。具体地，对于 $f = q \cdot g + r$，$q$ 和 $r$ 的系数可能比 $f$ 的系数大，因为除法在 $\mathbb{Z}$ 中不总是精确的。GCD 计算则使用伪余式序列（Pseudo-Remainder Sequence, PRS）来避免分数膨胀；对很大的整数系数可改用模 GCD（`gcd/modular.rs` 的 Brown 算法）。结式（resultant）的计算才使用子结式 PRS（`resultant.rs`）。

## 参考文献

1. **Cox, D., Little, J., O'Shea, D.** *Ideals, Varieties, and Algorithms* (4th ed.). Springer, 2015.
   - 第 1 章：多项式环的基本性质
   - 第 2 章：Gröbner 基的引入与 Buchberger 算法

2. **Gathen, J. von zur, Gerhard, J.** *Modern Computer Algebra* (3rd ed.). Cambridge University Press, 2013.
   - 第 6 章：多项式乘法算法（Karatsuba、FFT、NTT）
   - 第 21 章：多元多项式 GCD

3. **Adams, W. W., Loustaunau, P.** *An Introduction to Gröbner Bases*. AMS, 1994.
   - 第 1 章：单项式序与多元除法

## 参见

- [Gröbner 基理论](./groebner-theory.md) — 除法不唯一性的解决方案
- [多项式 GCD 与因式分解](./poly-gcd-factoring.md) — 一元和多元 GCD 算法
- [有限域与模算术](./finite-fields.md) — 有限域上的多项式运算
- [线性代数](./linear-algebra.md) — Bareiss 无分数行列式与高斯消元
- [Rust API：多项式](../api/rust-polynomials.md) — `DenseUnivariatePolynomial` 和 `SparseMultivariatePolynomial` 的完整 API
- [Rust API：系数域](../api/rust-domains.md) — `Domain` 和 `EuclideanDomain` trait
