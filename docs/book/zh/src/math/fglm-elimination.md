# 高阶：FGLM 与消元理论

本章系统讲解零维理想的 Gröbner 基换序算法（FGLM）、消元定理，以及基于 Gröbner 基的理想运算（商、饱和、交集、准素分解、根式）。这些工具构成了 oCAS 多项式系统求解和代数几何计算的核心基础设施。

---

## 前提知识

阅读本章前，建议具备以下基础：

- **Gröbner 基理论**：Buchberger 算法、F4/F5 算法、S-多项式、约化 Gröbner 基——参见 [Gröbner 基理论](./groebner-theory.md)
- **多项式代数**：多元多项式环 $\mathbb{F}[x_1, \ldots, x_n]$、单项式序（Lex / Grlex / Grevlex）——参见 [多项式代数](./polynomial-algebra.md)
- **线性代数**：高斯消元、向量空间维数、线性相关性——参见 [线性代数](./linear-algebra.md)
- **多项式因式分解**：无平方分解、一元多项式因式分解——参见 [多项式 GCD 与因式分解](./poly-gcd-factoring.md)

如需系统学习路径，请参考 [数学基础总览](./overview.md)。

---

## 基础概念

### 零维理想

**定义**。设 $I \subseteq R = \mathbb{F}[x_1, \ldots, x_n]$ 是一个理想。$I$ 称为**零维的**（zero-dimensional），如果商环 $R/I$ 作为 $\mathbb{F}$-向量空间是**有限维**的。

等价条件（以下任一成立即为零维）：

1. **代数几何**：仿射簇 $V(I) \subseteq \overline{\mathbb{F}}^n$ 是有限集
2. **Gröbner 基刻画**：对每个变量 $x_i$（$1 \leq i \leq n$），Gröbner 基中存在某个多项式 $g$ 使得 $\text{lm}(g) = x_i^{N_i}$（纯幂）
3. **阶梯有限性**：标准单项式（不被任何首项单项式整除的单项式）集合有限

**例**。$I = \langle x^2 - 1, y - x \rangle \subseteq \mathbb{Q}[x, y]$。其 Lex Gröbner 基为 $\{x^2 - 1, y - x\}$，首项为 $\{x^2, y\}$。对 $x$ 有 $x^2$，对 $y$ 有 $y$——纯幂条件满足，故 $I$ 是零维的。实际上 $V(I) = \{(1, 1), (-1, -1)\}$，有限。

**反例**。$J = \langle x - y \rangle \subseteq \mathbb{Q}[x, y]$。首项为 $\{x\}$，对 $y$ 无纯幂首项。$V(J) = \{(t, t) : t \in \overline{\mathbb{Q}}\}$ 是无穷集，故 $J$ 是正维的。

### 阶梯与标准单项式

**定义**。设 $G = \{g_1, \ldots, g_t\}$ 是理想 $I$ 关于某个单项式序 $\succ$ 的 Gröbner 基。令 $\text{LT}(I) = \langle \text{lm}(g_1), \ldots, \text{lm}(g_t) \rangle$ 为首项理想。

**标准单项式**（standard monomials）集合定义为

$$\text{Std}(I) = \{x^\alpha \in \mathbb{F}[x_1, \ldots, x_n] : x^\alpha \notin \text{LT}(I)\}$$

即不被任何 $\text{lm}(g_i)$ 整除的单项式。标准单项式的集合也称为**阶梯**（staircase），因为其形状在二维时像一个阶梯。

**核心性质**。标准单项式构成 $R/I$ 的一组 $\mathbb{F}$-向量空间基。即每个陪集 $f + I \in R/I$ 都有**唯一**的表示

$$f + I = \sum_{x^\alpha \in \text{Std}(I)} c_\alpha \cdot x^\alpha + I$$

其中 $c_\alpha \in \mathbb{F}$。这等价于说 $f$ 对 Gröbner 基的正规形是标准单项式的唯一线性组合。

**例**。$I = \langle x^2, xy \rangle \subseteq \mathbb{F}[x, y]$（Grevlex 序）。首项为 $\{x^2, xy\}$。

- $1$：不被 $x^2$ 或 $xy$ 整除 → 标准单项式 ✓
- $y$：不被整除 → ✓
- $y^2$：不被整除 → ✓
- $y^k$（$k \geq 0$）：不被整除 → ✓
- $x$：不被整除 → ✓
- $x^2$：被 $x^2$ 整除 → ✗
- $xy$：被 $xy$ 整除 → ✗

阶梯为 $\{1, x, y, y^2, y^3, \ldots\}$——无穷集！这反映了 $I$ 是正维理想（$V(I) = \{(0, t) : t \in \overline{\mathbb{F}}\}$）。

**维数**。对零维理想，阶梯的大小

$$D = \dim_{\mathbb{F}}(R/I) = |\text{Std}(I)|$$

是一个有限正整数，称为理想 $I$ 的**向量空间维数**或 **Hilbert 函数的稳定值**。这个维数 $D$ 在 FGLM 算法的复杂度分析中扮演核心角色。

### 正规形

**定义**。多项式 $f$ 关于 Gröbner 基 $G$ 的**正规形**（normal form）$\text{NF}(f, G)$ 是 $f$ 对 $G$ 做带余除法后的余数。它满足：

1. $\text{NF}(f, G) \in f + I$（与 $f$ 同余模 $I$）
2. $\text{NF}(f, G)$ 的每一项都是标准单项式
3. 表示唯一（对约化 Gröbner 基）

对零维理想，正规形可以表示为长度为 $D$ 的系数向量——每个分量对应一个标准单项式。这正是 FGLM 算法将正规形计算转化为线性代数操作的关键。

---

## 核心理论

### 消元定理

**消元理想**。设 $I \subseteq \mathbb{F}[x_1, \ldots, x_n]$，固定 $1 \leq \ell \leq n$。第 $\ell$ 个**消元理想**为

$$I_\ell = I \cap \mathbb{F}[x_{\ell+1}, \ldots, x_n]$$

即 $I$ 中不涉及前 $\ell$ 个变量的多项式构成的理想。

**消元定理**（Elimination Theorem）。设 $G$ 是 $I$ 关于 **Lex 序**的 Gröbner 基（$x_1 \succ x_2 \succ \cdots \succ x_n$）。则

$$G_\ell = G \cap \mathbb{F}[x_{\ell+1}, \ldots, x_n]$$

是 $I_\ell$ 的 Gröbner 基。

换言之，Lex Gröbner 基中不涉及 $x_1, \ldots, x_\ell$ 的那些多项式恰好生成了消元理想 $I_\ell$。

**证明思路**。Lex 序的关键性质是：若 $x^\alpha \succ x^\beta$ 且 $\alpha_1 = \cdots = \alpha_\ell = 0$，则 $\beta_1 = \cdots = \beta_\ell = 0$。因此首项理想 $\langle \text{lt}(g) : g \in G \rangle$ 与 $\mathbb{F}[x_{\ell+1}, \ldots, x_n]$ 的交集恰好由 $G_\ell$ 的首项生成。$\square$

**应用**。消元定理是以下所有构造的理论基础：

- **隐函数消元**：从方程组中消去变量，得到仅含剩余变量的关系
- **理想商的 Rabinowitsch 技巧**：引入辅助变量后消元
- **理想交集**：引入辅助变量 $t$ 后消元
- **多项式系统求解**：从 Lex GB 中逐变量提取约束

**注意**。Grevlex 序下的 Gröbner 基**不**具有消元性质。因此所有需要消元的操作都必须在 Lex 序下进行。这正是 FGLM 换序算法的核心价值：先用 Grevlex 快速计算基，再用 FGLM 转为 Lex 用于消元。

### FGLM 算法

FGLM 算法（Faugère–Gianni–Lazard–Mora, 1993）对**零维理想**提供高效的 Gröbner 基换序。其核心思想是：利用源序下的正规形计算，按目标序遍历单项式，通过线性代数方法直接构造目标序下的 Gröbner 基——完全避免在目标序下重新运行 F4。

#### 动机

假设我们已经用 Grevlex 序（最高效的计算序）得到了理想 $I$ 的 Gröbner 基 $G_{\text{grevlex}}$。为了做消元或求解，需要 Lex 序下的基 $G_{\text{lex}}$。

两种方法的对比：

| 方法 | 复杂度 | 说明 |
|------|--------|------|
| 重跑 F4（Lex） | 指数级（最坏） | 完全重新计算，中间多项式可能膨胀 |
| FGLM | $O(n \cdot D^3)$ | $D = \dim(R/I)$，纯线性代数 |

对零维理想，FGLM 几乎总是更快——尤其是当 $D$ 适中但 Lex 基的中间多项式很庞大时。

#### 算法描述

**输入**：零维理想 $I$ 在源序 $O_1$ 下的约化 Gröbner 基 $G_{\text{src}}$，目标序 $O_2$。

**输出**：$I$ 在目标序 $O_2$ 下的约化 Gröbner 基 $G_{\text{tgt}}$。

**初始化**：

1. 从 $G_{\text{src}}$ 提取所有首项单项式 $\text{lm}(g_1), \ldots, \text{lm}(g_t)$
2. 计算阶梯 $\text{Std}(I)$（BFS：从 $1$ 开始，逐变量扩展，跳过被首项整除者，直到不再有新标准单项式——零维保证终止）
3. 令 $D = |\text{Std}(I)|$，将标准单项式编号为 $s_1, s_2, \ldots, s_D$
4. 对每个变量 $x_i$，构造**乘法矩阵** $M_i$：对每个标准单项式 $s_j$，计算 $x_i \cdot s_j$ 的正规形（在 $G_{\text{src}}$ 下），得到 $D \times D$ 矩阵

**主循环**：按目标序 $O_2$ 递增遍历单项式：

```
seen_nfs ← []     // 已见正规形（D 维向量）
seen_mons ← []    // 对应的单项式
boundary ← {1}    // 边界集合
visited ← ∅       // 已访问单项式
GB_tgt ← ∅

while boundary ≠ ∅:
    m ← boundary 中 O2-最小的单项式
    visited ← visited ∪ {m}

    // 计算 m 在源序基下的正规形（D 维坐标向量）
    nf ← normal_form(m, G_src, staircase)

    if nf ∈ span(seen_nfs):      // 线性相关？
        // 找到关系 nf = Σ c_i · seen_nfs[i]
        coeffs ← Gaussian_solve(seen_nfs, nf)
        // 构造新多项式：m - Σ c_i · seen_mons[i]
        new_poly ← (m, 1) + Σ (-c_i, seen_mons[i])
        GB_tgt ← GB_tgt ∪ {new_poly}
        // m 的所有倍都在理想中，标记为已访问
        mark_multiples(visited, m)
    else:
        // 线性无关：m 是目标序下的新标准单项式
        seen_nfs.append(nf)
        seen_mons.append(m)
        // 将 m 的邻居加入边界
        for i = 1 to n:
            if m · x_i ∉ visited:
                boundary ← boundary ∪ {m · x_i}

return minimize(GB_tgt).auto_reduce()
```

#### 关键子程序

**正规形计算**。单项式 $m$ 在 $G_{\text{src}}$ 下的正规形是一个 $D$ 维向量。实现上，将 $m$ 视为系数为 $1$ 的单项式多项式，对 $G_{\text{src}}$ 做带余除法，然后将余数的每个项映射到阶梯的位置索引上。

**线性关系检测**。维护一个增广矩阵 $[\text{seen\_nfs} \mid \text{nf}]$（$D$ 行，$k+1$ 列，其中 $k$ 是已见正规形数）。每次加入新 nf 时做高斯消元：

- 若消元后增广列全零 → 线性相关，提取系数向量
- 若消元后增广列有非零元 → 线性无关，加入 seen_nfs

这可以用增量式高斯消元（每次只处理一列）高效实现，避免重复消元。

**倍标记**。当发现线性关系 $\text{nf}(m) = \sum c_i \cdot \text{nf}(m_i)$ 时，构造的多项式 $m - \sum c_i m_i \in I$ 的首项为 $m$（按目标序）。$m$ 的所有倍 $m \cdot x^\alpha$ 也必然在 $I$ 中，因此不需要再考虑——将它们从边界中移除。

#### 复杂度分析

- **阶梯计算**：$O(D)$ 个单项式，每个需要对基中每个首项做一次整除检查（每次检查 $O(n)$ 次比较）
- **乘法矩阵**：$n$ 个 $D \times D$ 矩阵，每个元素需一次正规形计算 → $O(n \cdot D^2 \cdot T_{\text{NF}})$
- **主循环**：至多发现 $D$ 个线性无关的标准单项式，每发现一个线性相关就输出一个目标基元素；总迭代数 $\le D + |G_{\text{tgt}}|$，实现中以 $D(D+1) + 4n$ 为安全上限
- **每次迭代**：一次正规形计算 $O(T_{\text{NF}})$ + 一次高斯消元 $O(D^2)$
- **总复杂度**：$O(n \cdot D^3)$ 次域运算

其中 $T_{\text{NF}}$ 是一次正规形计算的开销，对稀疏多项式通常与 $D$ 线性相关。

#### 正确性

**定理**（Faugère–Gianni–Lazard–Mora 1993）。设 $I$ 是零维理想，$G_{\text{src}}$ 是 $I$ 在序 $O_1$ 下的约化 Gröbner 基。FGLM 算法输出 $G_{\text{tgt}}$ 是 $I$ 在序 $O_2$ 下的约化 Gröbner 基。

**证明思路**：

1. 正规形计算与序无关——$G_{\text{src}}$ 的正规形映射与 $G_{\text{tgt}}$ 的正规形映射是同一个线性映射在不同基下的表示
2. 主循环按目标序递增遍历单项式，恰好覆盖所有标准单项式和产生首项的单项式
3. 线性关系检测保证：输出的多项式恰好生成首项理想
4. minimize 和 auto_reduce 保证约化性

### 理想运算

Gröbner 基使得理想的基本算术运算可以算法化。本节介绍 oCAS 实现的核心理想运算，每种都依赖消元定理。

#### 理想商（Ideal Quotient）

**定义**。两个理想 $I, J \subseteq R$ 的**理想商**为

$$I : J = \{f \in R : f \cdot g \in I, \, \forall g \in J\}$$

对单生成元 $J = \langle g \rangle$，简记 $I : g$。

**Rabinowitsch 技巧**。计算 $I : g$ 的关键技巧是引入一个新变量 $w$。但要注意：辅助变量消元构造实际给出的是**饱和理想**：

$$I : g^\infty = \left( I \cdot R[w] + \langle 1 - wg \rangle \right) \cap R$$

即在扩展环 $R[w]$ 中，将 $I$ 的生成元提升、添加 $1 - wg$、计算 Gröbner 基，然后取不涉及 $w$ 的多项式（消元）。当 $I$ 关于 $g$ 已饱和（即 $I : g = I : g^\infty$，例如无平方、无嵌入分支的情形）时，它恰等于 $I : g$。oCAS 的 `ideal_quotient` 实现的正是这个消元构造（数学上即 $I : g^\infty$）。

**正确性证明**。局部化观点：$R[w]/\langle 1 - wg \rangle \cong R[1/g]$（$w \mapsto 1/g$）。于是 $f \in (I \cdot R[w] + \langle 1-wg\rangle) \cap R$ 当且仅当 $f \in I \cdot R[1/g]$，即存在 $N$ 使得 $g^N f \in I$，亦即 $f \in I : g^N \subseteq I : g^\infty$。反之，若 $g^N f \in I$，则

$$f = f\bigl(1-(wg)^N\bigr) + w^N \cdot (g^N f)$$

第一项在 $\langle 1-wg\rangle$ 中（$1-(wg)^N$ 被 $1-wg$ 整除），第二项在 $I \cdot R[w]$ 中。$\square$

**多个生成元**。对 $J = \langle g_1, \ldots, g_m \rangle$，有

$$I : J^\infty = \bigcap_{j=1}^{m} (I : g_j^\infty)$$

oCAS 通过逐个计算每个 $I : g_j^\infty$ 再求交集实现（这正是 `ideal_quotient` 的实际行为）。

**例**。$I = \langle x^2, xy \rangle$，$J = \langle x \rangle$。

$$I : x = \langle x^2, xy \rangle : \langle x \rangle$$

在 $\mathbb{Q}[x, y, w]$ 中，计算 $\text{GB}(x^2, xy, 1 - wx)$ 并消去 $w$：

- $1 - wx \implies w = 1/x$
- $x^2 \cdot w = x$，$xy \cdot w = y$（由 $x^2 w - x = -x(1-wx)$ 与 $xy w - y = -y(1-wx)$）

消元后得到 $\langle x, y \rangle$，故 $\langle x^2, xy \rangle : \langle x \rangle = \langle x, y \rangle$。

#### 理想交集（Ideal Intersection）

**定义**。$I \cap J$ 是同时属于 $I$ 和 $J$ 的多项式构成的理想。

**辅助变量法**。引入新变量 $t$：

$$I \cap J = \langle t \cdot f_i, \, (1-t) \cdot g_j \rangle_{i,j} \cap R$$

其中 $\{f_i\}$ 和 $\{g_j\}$ 分别是 $I$ 和 $J$ 的生成元。

**直觉**。$t = 1$ 时约束退化为 $f_i = 0$（$I$ 的条件），$t = 0$ 时退化为 $g_j = 0$（$J$ 的条件）。消去 $t$ 后，得到同时满足两组条件的多项式。

**正确性**。$h \in I \cap J$ ⟺ $h = t \cdot h + (1-t) \cdot h$，其中 $t \cdot h \in t \cdot I$ 且 $(1-t) \cdot h \in (1-t) \cdot J$。反之，消元理想中的多项式在 $t = 0$ 和 $t = 1$ 处都为零，故同时属于 $I$ 和 $J$。$\square$

**例**。$\langle x \rangle \cap \langle y \rangle = \langle xy \rangle$。

在 $\mathbb{Q}[x, y, t]$ 中，计算 $\text{GB}(tx, (1-t)y)$ 并消去 $t$：

- $tx = 0 \implies t = 0$ 或 $x = 0$
- $(1-t)y = 0 \implies t = 1$ 或 $y = 0$

消元后得到 $\langle xy \rangle$。$\square$

#### 理想饱和（Ideal Saturation）

**定义**。$I$ 关于 $J$ 的**饱和**为

$$I : J^\infty = \bigcup_{k=1}^{\infty} (I : J^k)$$

等价地，$I : J^\infty$ 是 $I$ 在 $V(J)$ 之外的"限制"——在 $J$ 的零点集之外与 $I$ 一致。

**迭代计算**。由升链条件，存在 $k_0$ 使得 $I : J^{k_0} = I : J^{k_0+1} = \cdots = I : J^\infty$。算法反复计算 $I : J$，直到稳定：

```
I_old ← I
loop:
    I_new ← I_old : J
    if I_new == I_old:  // 比较 Gröbner 基
        return I_new
    I_old ← I_new
```

**应用**。饱和在代数几何中极为重要：

- **根式计算**：$\sqrt{I} = I : h^\infty$（正维情形，$h$ 与 Jacobian 相关）
- **准素分解**：分离不同维度的分支
- **闭包计算**：射影闭包、完备交

**例**。$\langle x^2 y, xy^2 \rangle : \langle x \rangle^\infty$。

按数学定义逐轮计算商：

第一轮：$\langle x^2 y, xy^2 \rangle : \langle x \rangle = \langle xy, y^2 \rangle$

第二轮：$\langle xy, y^2 \rangle : \langle x \rangle = \langle y \rangle$

第三轮：$\langle y \rangle : \langle x \rangle = \langle y \rangle$（稳定）

故 $\langle x^2 y, xy^2 \rangle : \langle x \rangle^\infty = \langle y \rangle$。$\square$

（注：oCAS 的 `ideal_quotient` 用消元构造单次计算 $I : g^\infty$——本例中一次调用即得到 $\langle y \rangle$，`ideal_saturate` 的迭代仅用于确认稳定。）

#### 理想和、积与成员判定

这些是最基本的理想运算，实现相对简单：

| 运算 | 定义 | 实现 |
|------|------|------|
| $I + J$ | $\langle f_1, \ldots, f_m, g_1, \ldots, g_k \rangle$ | 合并生成元，计算 GB |
| $I \cdot J$ | $\langle f_i g_j \rangle_{i,j}$ | 所有生成元两两乘积，计算 GB |
| $f \in I$ ? | — | 计算 $\text{NF}(f, G_I)$，是否为零 |

### 准素分解

**定义**。理想 $Q \subseteq R$ 称为**准素理想**（primary ideal），如果 $fg \in Q$ 蕴含 $f \in Q$ 或 $g^n \in Q$（对某个 $n \geq 1$）。

**准素分解定理**（Lasker–Noether）。每个理想 $I \subseteq \mathbb{F}[x_1, \ldots, x_n]$（$\mathbb{F}$ 特征零）都有准素分解

$$I = Q_1 \cap Q_2 \cap \cdots \cap Q_r$$

其中每个 $Q_i$ 是准素理想。对应的**关联素理想**为 $\mathfrak{p}_i = \sqrt{Q_i}$。

#### 零维准素分解

对零维理想，准素分解可以通过 Lex Gröbner 基中一元多项式的因式分解来实现。

**核心观察**。零维理想的 Lex GB 包含一个仅含 $x_1$ 的多项式 $p_1(x_1)$、一个含 $x_1, x_2$ 的多项式 $p_2(x_1, x_2)$，以此类推。$p_1(x_1)$ 的因式分解对应了理想在 $x_1$ 方向上的不同分支。

**算法**（对应 oCAS 实现 `primary_decomp_zero_dim`）：

1. 计算 $I$ 的 Lex Gröbner 基 $G$
2. 提取仅含 $x_1$ 的多项式 $p_1(x_1)$，计算其**无平方部分** $\tilde{p}_1 = p_1 / \gcd(p_1, p_1')$
3. 将 $\tilde{p}_1$ 在 $\mathbb{Q}$ 上分解为互异不可约因子 $\tilde{p}_1 = q_1 \cdot q_2 \cdots q_s$
4. 对每个因子 $q_i$，用**饱和**分离分支：$I_i = I : \left(\prod_{j \neq i} q_j\right)^\infty$（源码逐个因子依次饱和）
5. 每个分支 $I_i$ 的关联素理想取为 $\mathfrak{p}_i = \mathrm{GB}(I + \langle q_i \rangle)$

**注意**。若 $\tilde{p}_1$ 只有一个不可约因子（或 $G$ 中不存在仅含 $x_1$ 的多项式），算法返回单一分支 $\{\text{primary} = I,\; \text{prime} = \sqrt{I}\}$。当前实现只处理第一个变量 $x_1$ 的因子分解，**不递归**处理后续变量。

**例**。$I = \langle x^2 - 1, y - x \rangle \subseteq \mathbb{Q}[x, y]$。

Lex GB：$\{x^2 - 1, y - x\}$。$p_1(x) = x^2 - 1$ 已无平方，分解为 $q_1 = x - 1$、$q_2 = x + 1$。

- 分支 1：$I_1 = I : \langle x + 1 \rangle^\infty = \langle x - 1, y - 1 \rangle$（关联素理想 $\langle x - 1, y - 1 \rangle$）
- 分支 2：$I_2 = I : \langle x - 1 \rangle^\infty = \langle x + 1, y + 1 \rangle$（关联素理想 $\langle x + 1, y + 1 \rangle$）

检查：$\langle x^2 - 1, y - x \rangle = \langle x - 1, y - 1 \rangle \cap \langle x + 1, y + 1 \rangle$，两个分支都是素理想。

作为对照，$I = \langle x^2, xy \rangle$ 的 $p_1(x) = x^2$ 无平方部分只有 $x$ 一个因子，故算法返回单一分支 $\{\text{primary} = \langle x^2, xy \rangle,\; \text{prime} = \sqrt{I} = \langle x \rangle\}$。不过该理想仍有更细的分解 $\langle x^2, xy \rangle = \langle x \rangle \cap \langle x^2, y \rangle$，其中 $\langle x \rangle$ 是素理想，$\langle x^2, y \rangle$ 是关联素理想为 $\langle x, y \rangle$ 的准素理想（嵌入分支）。

#### 正维准素分解

正维理想的准素分解更加复杂，oCAS 目前标记为待实现（`TODO`）。完整实现需要：

- Gianni–Trager–Zacharias 算法
- 或 Eisenbud–Huneke–Vasconcelos 的特征方法

### 根式（Radical）

**定义**。理想 $I$ 的**根式**为

$$\sqrt{I} = \{f \in R : f^n \in I \text{ 对某个 } n \geq 1\}$$

等价地，$\sqrt{I}$ 是包含 $I$ 的最大理想使得 $V(\sqrt{I}) = V(I)$（相同的零点集）。

#### 零维根式

对零维理想，根式计算利用 Lex GB 中一元多项式的无平方分解。

**算法**：

1. 计算 Lex GB $G$
2. 对 $G$ 中每个仅含 $x_i$ 的多项式 $p_i(x_i)$，计算其无平方部分 $\tilde{p}_i = p_i / \gcd(p_i, p_i')$
3. 用 $\tilde{p}_i$ 替换 $p_i$，得到 $\sqrt{I}$ 的生成元

**正确性**。零维理想的根式对应于将每个关联素理想的幂次降为 $1$。无平方分解恰好实现了这一点——$p_i^{e_i}$ 的无平方部分是 $p_i$。

**例**。$\sqrt{\langle x^2, xy \rangle}$。

Lex GB：$\{x^2, xy\}$。$x^2$ 的无平方部分：$x^2 / \gcd(x^2, 2x) = x^2 / x = x$。

根式：$\langle x, xy \rangle = \langle x \rangle$（因为 $x$ 已经整除 $xy$）。$\square$

#### 正维根式：Jacobian 饱和

对正维理想，oCAS 使用基于 Jacobian 的饱和方法（简化版 Kemper 算法）：

$$\sqrt{I} = I : h^\infty$$

其中 $h$ 与 Jacobian 行列式相关。

**算法**（特征零，对应实现 `radical_via_jacobian`）：

1. 对 $I$ 的每个生成元 $f_i$ 和每个变量 $x_j$，计算偏导数 $\partial f_i / \partial x_j$（跳过常值/零偏导数）
2. 启发式选取 $h$：源码**并非**计算真正的 $\gcd$，而是对非平凡偏导数做 `reduce` 折叠、保留**总次数最小**的那个。这是对经典 Jacobian 饱和公式（$\sqrt{I} = I : h^\infty$，其中 $h$ 需取合适的 Jacobian 相关元素）的**启发式近似**——$h$ 只是某个偏导数（Jacobian 矩阵的一个元素），既不保证整除 Jacobian 行列式，也不保证饱和结果恰好等于 $\sqrt{I}$
3. 计算 $I : h^\infty$（迭代饱和）

**理论基础**。在正则点处，Jacobian 矩阵的秩等于簇的维数（Jacobian 判据）；奇异点正是秩不足（小于 $n - \dim$）的点——在这些点处幂次结构非平凡。饱和 $I : h^\infty$ 的作用是去除这些幂次。

**限制**。当所有偏导数都是常数/零，或 $h$ 是常值多项式（总次数为 0）时，算法回退到返回原始 GB。由于 $h$ 是启发式选取而非真正的 Jacobian 相关元素，饱和结果只是 $\sqrt{I}$ 的近似（可能偏大，也可能偏小、漏掉分量）——对正维理想，完整根式计算需要更精细的算法（如 Gianni–Trager–Zacharias 或 Eisenbud–Huneke–Vasconcelos）。

### 素性与准素性判定

**零维素性判定**。零维理想 $I$ 是素理想的充要条件是：Lex GB 中每个仅含 $x_i$ 的多项式都是不可约的。

oCAS 实现（`is_prime_zero_dim`）：

1. 计算 Lex GB
2. 提取每个一元多项式 $p_i(x_i)$
3. 对次数 $\le 3$ 的多项式用**有理根定理**检查可约性（$\mathbb{Q}$ 上次 $d \le 3$ 的多项式不可约 ⟺ 无有理根）；次数 $> 3$ 的多项式不做完全分解检验

**注意**。正维理想目前保守返回 `false`——素的正维理想会被误判为非素（**假阴性**）。而零维情形下，次数 $> 3$ 的一元多项式若可约但无有理根（如 $(x^2-2)(x^2-3)$），会被误判为不可约，从而可能产生**假阳性**（把非素理想报为素）。

**准素性判定**。理想 $I$ 是准素理想当且仅当它恰好有一个关联素理想。实现上通过准素分解检查分量数。

---

## 在 oCAS 中的实现

### FGLM 实现

FGLM 算法实现在 `ocas-poly/src/groebner/fglm.rs` 中。

#### 入口函数

```rust
pub fn fglm<D: Domain, O2: MonomialOrder>(
    gb: &GroebnerBasis<D, impl MonomialOrder>,
) -> Option<GroebnerBasis<D, O2>>
```

泛型参数：
- `D`：系数域（实现 `Domain` trait）
- `O2`：目标单项式序
- 输入基的源序由 `gb` 的类型参数 `impl MonomialOrder` 确定

返回 `None` 当理想不是零维的（阶梯计算发现无穷多个标准单项式）。

#### 阶梯计算

```rust
fn compute_staircase(lms: &[Vec<usize>], n_vars: usize) -> Option<Vec<Vec<usize>>>
```

BFS 算法：

- 从全零指数向量 $1 = x_1^0 \cdots x_n^0$ 开始
- 对队列中每个单项式 $m$，检查是否被某个 $\text{lm}(g_i)$ 整除
- 若不被整除 → 加入阶梯，将其 $n$ 个邻居（$m$ 的每个变量指数 $+1$）加入队列
- 安全阈值 100,000：若 BFS 访问超过此数，判定为正维

整除检查使用逐分量比较：

```rust
fn monomial_divides_big(lm: &[usize], big: &[usize]) -> bool {
    lm.iter().zip(big.iter()).all(|(a, b)| a <= b)
}
```

#### 正规形计算

```rust
fn normal_form_monomial<D: Domain>(
    m: &[usize],
    gb: &GroebnerBasis<D, impl MonomialOrder>,
    staircase: &[Vec<usize>],
    domain: &D,
) -> Vec<D::Element>
```

将单项式 $m$ 构造为单项式多项式 $1 \cdot x^m$，对 GB 做带余除法，然后将余数的每个项映射到阶梯的位置索引上。返回长度为 $D$ 的坐标向量。

#### 线性关系检测

```rust
fn find_relation<D: Domain>(
    seen: &[Vec<D::Element>],
    nf: &[D::Element],
    domain: &D,
) -> Option<Vec<D::Element>>
```

构造增广矩阵 $[\text{seen}^T \mid \text{nf}]$（$D$ 行，$k+1$ 列），执行前向消元 + 回代：

1. 逐列找主元（非零元素），交换行
2. 归一化主元行，消去其余行的该列
3. 检查一致性：若某行左侧全零但右侧非零 → 无解（线性无关）
4. 回代提取系数向量

返回 `Some(coeffs)` 当线性相关，`None` 当线性无关。

#### 倍标记

```rust
fn mark_multiples(
    visited: &mut HashMap<Vec<usize>, bool>,
    m: &[usize],
    n_vars: usize,
    max_deg: usize,
)
```

BFS 标记 $m$ 的所有倍（总次数不超过 `max_deg`）为已访问。`max_deg` 设为 $2D$——足够覆盖所有可能产生新首项的单项式。

#### 主循环细节

主循环维护：

- `boundary`：候选单项式的集合（初始为 $\{1\}$，即全零指数向量）
- `seen_nfs`：已见正规形向量的列表
- `seen_mons`：对应的单项式列表
- `new_basis`：输出的 Gröbner 基

每步取 boundary 中目标序最小的单项式，计算正规形，检测线性关系。若相关 → 输出新多项式并标记倍；若无关 → 加入 seen 并扩展邻居。

循环终止条件：boundary 为空或步数超过 $D(D+1) + 4n$（安全上限）。

最终对 `new_basis` 执行 `minimize()`（移除冗余）和 `auto_reduce()`（互相约化），得到约化 Gröbner 基。

### 理想运算实现

理想运算实现在 `ocas-poly/src/ideal.rs` 中，统一使用 Lex 序以保证消元性质。

#### 成员判定

```rust
pub fn ideal_contains<D: Domain + 'static>(
    generators: &[SparseMultivariatePolynomial<D, Lex>],
    f: &SparseMultivariatePolynomial<D, Lex>,
    algo: Algorithm,
) -> bool
```

计算 GB 后对 $f$ 求正规形，检查是否为零。

#### 理想商

```rust
pub fn ideal_quotient<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

内部流程：

1. 对 $J$ 的每个生成元 $g_j$，调用 `quotient_single_generator`（Rabinowitsch 技巧）
2. `quotient_single_generator` 实现：
   - 将 $I$ 的生成元提升到 $R[w]$
   - 添加 $1 - wg$
   - 计算 Lex GB
   - 消去 $w$（取不涉及 $w$ 的多项式）
3. 对所有 $I : g_j^\infty$ 的结果取交集（`intersect_generators`）

#### 理想交集

```rust
pub fn ideal_intersection<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

内部调用 `intersect_generators`：

1. 引入辅助变量 $t$（index 0）
2. 构造生成元 $\{t \cdot f_i\} \cup \{(1-t) \cdot g_j\}$
3. 在扩展环 $R[t]$ 中计算 Lex GB
4. 消去 $t$，得到 $I \cap J$ 的生成元

#### 理想饱和

```rust
pub fn ideal_saturate<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

迭代调用 `ideal_quotient` 直到 GB 稳定（双向包含判定：两边基元素互相约化为零），上限 20 轮；超限后返回当前结果。

#### 准素分解

```rust
pub fn primary_decomposition(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> Vec<PrimaryComponent>
```

目前仅实现零维情况（`primary_decomp_zero_dim`），且只处理第一个变量：

1. 计算 Lex GB
2. 提取仅含 $x_1$ 的一元多项式 $p_1(x_1)$，计算其无平方部分
3. 对无平方部分做 $\mathbb{Q}$ 上因式分解
4. 用饱和分离各分支（对每个其他因子依次 `ideal_saturate`）
5. 分支的关联素理想取 $\mathrm{GB}(I + \langle q_i \rangle)$；只有一个因子时返回单一分支（$\text{prime} = \sqrt{I}$）

`PrimaryComponent` 结构体包含：

```rust
pub struct PrimaryComponent {
    pub primary: GroebnerBasis<RationalDomain, Lex>,  // 准素理想
    pub prime: GroebnerBasis<RationalDomain, Lex>,     // 关联素理想
}
```

#### 根式

```rust
pub fn ideal_radical(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> GroebnerBasis<RationalDomain, Lex>
```

分两种情况：

- **零维**（`radical_zero_dim`）：对 Lex GB 中的每个一元多项式做无平方分解，替换后重新生成 GB
- **正维**（`radical_via_jacobian`）：计算所有生成元对所有变量的非平凡偏导数，取 $h$ 为其中**总次数最小**者（启发式，非真正的 $\gcd$），计算 $I : h^\infty$

#### 其他判定

```rust
pub fn is_zero_dimensional(gb: &GroebnerBasis<RationalDomain, Lex>) -> bool
pub fn is_prime_ideal(generators: &[...]) -> bool
pub fn is_primary_ideal(generators: &[...]) -> bool
```

- `is_zero_dimensional`：检查每个变量是否有纯幂首项
- `is_prime_ideal`：零维检查不可约性，正维保守返回 `false`
- `is_primary_ideal`：检查是否恰好有一个关联素理想

### 零维系统求解

```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

算法流程：

1. 计算 Lex GB（若输入不是 GB）
2. 检查零维性
3. 三角分解求解（`solve_triangular`）——从 Lex 序中**最小的变量**（$x_n$）开始向前回代：
   - 提取仅含 $x_n$ 的多项式 → 实根隔离 + 数值求根
   - 对每个根，代入涉及 $x_{n-1}$ 的多项式 → 求 $x_{n-1}$ 的根
   - 递归直到所有变量确定（结果按 $x_1, \ldots, x_n$ 顺序返回）

返回类型：

```rust
pub enum PolynomialSystemSolution {
    ZeroDimensional(ZeroDimSolutions),  // 有限个解
    PositiveDimensional(GroebnerBasis),  // 正维（无穷多解）
    Empty,                               // 无解
}
```

---

## 进阶话题

### FGLM 与 F4 的选择

| 场景 | 推荐方法 | 原因 |
|------|----------|------|
| 零维 + 需要 Lex GB | Grevlex F4 → FGLM | $O(n \cdot D^3)$ 远快于 Lex F4 |
| 零维 + 需要 Grevlex GB | 直接 F4 | 无需换序 |
| 正维 + 需要 Lex GB | Lex F4 或 `reorder` | FGLM 不适用 |
| 消元 + 零维 | FGLM → 消元 | 消元定理保证正确性 |
| 消元 + 正维 | Lex F4 | FGLM 不适用 |

### 倍域上的 FGLM

oCAS 的 FGLM 实现对任意域 $D$ 泛型工作，正规形计算完全经由 `Domain` trait 的运算（有限域元素以 $[0, p)$ 的规范代表存储）。注意：F4/F5 中针对 $\mathbb{F}_p$ 的 `FpPoly` `i64` 快速路径并不用于 FGLM。

### 局限性

1. **仅限零维**：正维理想的阶梯无穷，FGLM 不适用
2. **$D^3$ 瓶颈**：当向量空间维数 $D$ 很大时（如 $D > 10^4$），高斯消元成为瓶颈
3. **内存**：需要存储 $D \times D$ 乘法矩阵和增广矩阵

### 与文献的关系

oCAS 的实现直接基于原始论文：

- FGLM：Faugère, Gianni, Lazard, Mora (1993), *Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering*, JSC
- 理想运算：Cox, Little, O'Shea, *Ideals, Varieties, and Algorithms*, Chapters 4, 8
- 准素分解：Gianni, Trager, Zacharias (1988), *Gröbner Bases and Primary Decomposition of Polynomial Ideals*
- 根式计算：Kemper (2002), *A Course in Commutative Algebra*

---

## 参考文献

1. J.-C. Faugère, P. Gianni, D. Lazard, T. Mora. "Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering." *Journal of Symbolic Computation*, 16(4):329–344, 1993.
2. D. Cox, J. Little, D. O'Shea. *Ideals, Varieties, and Algorithms*. 4th ed., Springer, 2015. Chapters 3 (消元), 4 (理想商与饱和), 8 (准素分解).
3. W. W. Adams, P. Loustaunau. *An Introduction to Gröbner Bases*. AMS, 1994.
4. P. Gianni, B. Trager, G. Zacharias. "Gröbner Bases and Primary Decomposition of Polynomial Ideals." *Journal of Symbolic Computation*, 6(2–3):149–167, 1988.
5. G. Kemper. *A Course in Commutative Algebra*. Springer, 2011.
6. T. Becker, V. Weispfenning. *Gröbner Bases: A Computational Approach to Commutative Algebra*. Springer, 1993.

---

**参见**：

- [Rust API：Gröbner 基与理想](../api/rust-groebner.md) — 函数签名与完整示例
- [Gröbner 基理论](./groebner-theory.md) — Buchberger/F4/F5 算法的理论基础
- [多项式代数](./polynomial-algebra.md) — 单项式序与多元带余除法
- [多项式 GCD 与因式分解](./poly-gcd-factoring.md) — 无平方分解与因式分解算法
- [求解器](../solvers.md) — 多项式系统求解与 ODE
