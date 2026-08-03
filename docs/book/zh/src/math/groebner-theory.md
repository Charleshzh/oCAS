# 进阶：Gröbner 基理论

Gröbner 基是多元多项式理想论的核心工具。本章从单项式理想和 S-多项式的
基本定义出发，系统讲解 Buchberger 算法、F4 矩阵方法和 F5 签名方法的理论
基础，并说明 oCAS 如何在 Rust 中实现这些算法。

---

## 前提知识

阅读本章前，建议具备以下基础：

- **多项式代数**：多元多项式环 $\mathbb{F}[x_1, \ldots, x_n]$、单项式序（Lex / Grlex / Grevlex）、多元带余除法——参见[多项式代数](./polynomial-algebra.md)
- **线性代数**：矩阵运算、高斯消元、行阶梯形——参见[线性代数](./linear-algebra.md)
- **有限域**：$\mathbb{F}_p$ 的构造与基本运算——参见[有限域与模算术](./finite-fields.md)
- **理想与商环**：理想 $I \subseteq \mathbb{F}[x_1, \ldots, x_n]$、商环 $R/I$ 的基本概念

如需系统学习路径，请参考[数学基础总览](./overview.md)中的路径 A。

---

## 基础概念

### 单项式理想

**定义**。设 $R = \mathbb{F}[x_1, \ldots, x_n]$。由单项式生成的理想

$$M = \langle x^{\alpha^{(1)}}, x^{\alpha^{(2)}}, \ldots, x^{\alpha^{(s)}} \rangle$$

称为**单项式理想**（monomial ideal）。$f \in M$ 当且仅当 $f$ 的每个非零项的单项式都被某个 $x^{\alpha^{(i)}}$ 整除。

**例**。$M = \langle x^2, xy \rangle \subseteq \mathbb{F}[x, y]$。多项式 $f = x^2 y + x y^2$ 的各项：$x^2 y$ 被 $x^2$ 整除，$xy^2$ 被 $xy$ 整除，故 $f \in M$。而 $g = y^3$ 的单项式 $y^3$ 不被 $x^2$ 或 $xy$ 整除，故 $g \notin M$。

### Dickson 引理

**引理（Dickson）**。$\mathbb{N}^n$ 的每个子集 $S$ 都有有限的**极小元**集合 $M \subseteq S$，使得对每个 $\alpha \in S$，存在 $\beta \in M$ 满足 $\beta \leq \alpha$（逐分量比较）。

等价表述：$\mathbb{F}[x_1, \ldots, x_n]$ 中的每个单项式理想都由有限个单项式生成。

**证明思路**（对 $n$ 归纳）。

- $n = 1$：$\mathbb{N}$ 的非空子集有最小元，取该最小元即得 $M$。
- $n \to n+1$：对 $S \subseteq \mathbb{N}^{n+1}$，对每个 $a \in \mathbb{N}$ 考虑截面 $S_a = \{(\alpha_2, \ldots, \alpha_{n+1}) : (a, \alpha_2, \ldots, \alpha_{n+1}) \in S\}$。由归纳假设，每个 $S_a$ 有有限极小元。进一步分析"被更小的第一分量支配"的情况，可构造 $S$ 的有限极小元集。

> **意义**：Dickson 引理保证了单项式理想的有限生成性，这直接推出 Hilbert 基定理——每个多项式理想都是有限生成的。Gröbner 基理论正是在此基础上建立了"找到合适生成元"的系统方法。

### S-多项式

**定义**。设 $f, g \in R$ 为非零多项式。令 $\gamma = \text{lcm}(\text{lm}(f), \text{lm}(g))$——指数向量逐分量取最大值。$f$ 和 $g$ 的 **S-多项式**（S-polynomial）为

$$S(f, g) = \frac{x^\gamma}{\text{lt}(f)} \cdot f - \frac{x^\gamma}{\text{lt}(g)} \cdot g = \frac{x^\gamma}{\text{lm}(f)} \cdot f - \frac{x^\gamma}{\text{lm}(g)} \cdot g$$

S-多项式的关键性质：它的首项被**消去**了。直觉上，$S(f,g)$ 度量了 $f$ 和 $g$ 在首项之外的"不一致性"。当 $f$ 和 $g$ 的首项互素（即 $\text{lcm} = \text{lm}(f) \cdot \text{lm}(g)$）时，$S(f,g)$ 通常可以直接化简为零——这是 Buchberger 的第一条优化判据。

**例**。在 $\mathbb{F}[x, y]$ 中，取 $f = x^2 - y$，$g = xy - 1$（Grevlex 序）。

$$\gamma = \text{lcm}(x^2, xy) = x^2 y$$

$$S(f, g) = \frac{x^2 y}{x^2} \cdot (x^2 - y) - \frac{x^2 y}{xy} \cdot (xy - 1) = y(x^2 - y) - x(xy - 1) = -y^2 + x$$

### Gröbner 基的定义

**定义**。固定单项式序 $\succ$。理想 $I$ 的生成集 $\{g_1, \ldots, g_t\}$ 称为 $I$ 关于 $\succ$ 的 **Gröbner 基**（Gröbner basis），如果

$$\langle \text{lt}(g_1), \ldots, \text{lt}(g_t) \rangle = \langle \text{lt}(f) : f \in I \rangle$$

即生成元的首项生成的理想恰好等于 $I$ 中所有多项式首项生成的理想（$I$ 的**首项理想** $L(I)$）。

**等价刻画**（Buchberger 判据）。$\{g_1, \ldots, g_t\}$ 是 $I$ 的 Gröbner 基当且仅当对每对 $i \neq j$，$S(g_i, g_j)$ 对 $\{g_1, \ldots, g_t\}$ 的余数为零。

**Gröbner 基的核心价值**：

1. **成员判定**：$f \in I$ 当且仅当 $f$ 对 Gröbner 基的余数为零（正规形唯一性）
2. **理想相等判定**：两个理想 $I = J$ 当且仅当它们的（约化）Gröbner 基相同
3. **消元**：Lex 序下的 Gröbner 基天然具有消元性质（消元理想的生成元出现在基中）
4. **商环结构**：$R/I$ 的 $k$-向量空间基由不被首项理想整除的单项式给出（阶梯）

**约化 Gröbner 基**。Gröbner 基 $\{g_1, \ldots, g_t\}$ 是**约化的**（reduced），如果：

1. 每个 $g_i$ 是首一的（$\text{lc}(g_i) = 1$）
2. 对 $i \neq j$，$g_i$ 的任何非零项都不被 $\text{lm}(g_j)$ 整除

约化 Gröbner 基是唯一的（固定序和域下）。oCAS 的所有算法都输出约化 Gröbner 基。

---

## 核心理论

### Buchberger 算法

Buchberger 算法（1965）是 Gröbner 基计算的经典方法。其思想朴素而直接：
反复构造 S-多项式、化简，直到所有 S-多项式化简为零。

```
输入: F = {f1, ..., fs}（理想生成元）
输出: G（Gröbner 基）

G ← F
P ← {(fi, fj) | i < j}  // 初始临界对集合
while P ≠ ∅:
    从 P 中取出一对 (f, g)
    h ← NF(S(f, g), G)   // S-多项式对 G 的余数
    if h ≠ 0:
        P ← P ∪ {(h, g') | g' ∈ G}  // 新增临界对
        G ← G ∪ {h}
return G
```

**终止性**：由 Buchberger 判据保证——每增加一个新元素，首项理想严格扩大，而 Hilbert 基定理保证链终止。实际实现中还需 `minimize`（移除冗余元素）和 `auto_reduce`（互相约化），得到约化基。

**效率问题**：朴素 Buchberger 算法的主要瓶颈在于——

1. **临界对数量**：$O(|G|^2)$ 个临界对，大部分是冗余的
2. **逐个处理**：每次只化简一个 S-多项式，无法利用批处理
3. **反复扫描**：寻找约化器需要线性扫描整个基

以下三节分别介绍解决这些问题的方法：Gebauer–Moeller 筛选（Buchberger 加速）、F4（批处理）、F5（避免零约化）。

### Buchberger 判据与优化

Buchberger 的两个判据可以大幅减少需要处理的临界对数量。

**第一判据（首项互素判据）**。若 $\text{lm}(f_i)$ 和 $\text{lm}(f_j)$ 互素（即 $\gcd = 1$），则

$$S(f_i, f_j) \xrightarrow{G} 0$$

即该 S-多项式直接化简为零，可以从临界对集合中删除。

**证明**：设 $\text{lm}(f_i) = m_i$，$\text{lm}(f_j) = m_j$，$\gcd(m_i, m_j) = 1$。则 $\text{lcm}(m_i, m_j) = m_i m_j$，且 $S(f_i, f_j) = m_j f_i - m_i f_j$。对其逐项化简时，每个非零项 $c \cdot x^\alpha$ 满足 $x^\alpha$ 被 $m_i$ 或 $m_j$ 之一整除，因此可以被对应的 $f_i$ 或 $f_j$ 化简，最终余数为零。$\square$

**第二判据（链判据 / Buchberger 判据的推广）**。设 $f_k \in G$。若 $\text{lm}(f_k)$ 整除 $\text{lcm}(\text{lm}(f_i), \text{lm}(f_j))$，且 $\text{lcm}(\text{lm}(f_i), \text{lm}(f_k))$ 和 $\text{lcm}(\text{lm}(f_j), \text{lm}(f_k))$ 都**严格**整除 $\text{lcm}(\text{lm}(f_i), \text{lm}(f_j))$，则

$$S(f_i, f_j) \xrightarrow{G} 0$$

可以不处理临界对 $(f_i, f_j)$。

> **冗余对**：更一般地，若临界对 $(f_i, f_j)$ 的 $\text{lcm}$ 被另一个已有临界对的 $\text{lcm}$ 整除，则 $(f_i, f_j)$ 是冗余的——化简另一个已隐含此对的处理。

### Gebauer–Moeller 临界对管理

Gebauer 和 Moeller（1988）将 Buchberger 判据系统化，提供了高效的临界对管理框架。在新多项式 $f_{\text{new}}$ 加入基 $G = \{f_1, \ldots, f_m\}$ 时：

1. **生成候选对**：$P_{\text{new}} = \{(f_{\text{new}}, f_i) \mid 1 \leq i \leq m\}$
2. **第一判据过滤**：删除 $\text{lm}(f_{\text{new}})$ 与 $\text{lm}(f_i)$ 互素的对
3. **第二判据过滤**：删除 $\text{lcm}$ 被其他（仍存活）对的 $\text{lcm}$ 整除的对（整除关系含相等，重复对会被丢弃）
4. **冗余对清理**：将被新对的 $\text{lcm}$ 整除的旧对标记为冗余
5. **合并**：将筛选后的 $P_{\text{new}}$ 合并到主临界对集合

关键优化：按 $\text{lcm}$ 的大小（用单项式序比较）对临界对排序，优先处理
小的——这使得更大的对更可能被第二判据过滤。

### F4 算法：矩阵行阶梯批处理

F4 算法（Faugère 1999）的核心思想是：**将多个 S-多项式的化简批处理为一次稀疏矩阵的行阶梯运算**。

#### 符号预处理

F4 不直接对 S-多项式做除法，而是先进行**符号预处理**（symbolic preprocessing）：对选定的临界对集合，找出所有需要的约化器（基元素的单项式倍），将它们的非零项展开并映射到列编号。

具体步骤：

1. **选择临界对**：从队列中选取一批临界对（按 $\text{lcm}$ 排序）
2. **收集单项式**：对每个选出的临界对，计算 S-多项式涉及的单项式；对每个需要约化的单项式，查找基中的约化器并将其乘以相应的单项式倍展开
3. **建立列映射**：将所有出现的单项式按序排列，编号为 $0, 1, \ldots, N-1$
4. **构造矩阵**：每行对应一个多项式（S-多项式或约化器的倍），每列对应一个单项式，矩阵元素为系数

这产生一个**稀疏矩阵** $M$，其中每行的非零元素对应多项式的非零项。

#### 行阶梯运算

对矩阵 $M$ 执行高斯消元（行阶梯化）：

1. 按行的首非零列排序
2. 逐列消元：找到首非零元素为当前列的行作为主元，将其余行的该列消为零
3. 消元后的行中，非零行对应的多项式即为新基元素

**ℓ_p 原生快速路径**。当系数域为有限域 $\mathbb{F}_p$（$p < 2^{31}$）时，所有矩阵运算使用 `i64` 模算术——无 BigInt 开销。这包括：

- 稀疏行减法 `sub_scaled_fp`：归并两列升序行，跳过首列（已消去），逐项 $(row_j - c \cdot pivot_j) \bmod p$，丢弃零项
- 模逆查找 `mod_inv`：扩展 Euclid
- 模归约 `norm_mod`：确保系数在 $[0, p)$

#### 提取新基元素

消元后，矩阵的每个非零行 $r$ 对应一个多项式 $h_r$。对 $h_r$：

1. 若 $h_r$ 的首项是"新"的（不在之前已收集的首项集合中），则加入基
2. 用尾部约化（tail reduction）进一步简化：用已有基元素化简 $h_r$ 的尾项

#### F4 的复杂度

F4 的主要开销在于行阶梯运算——对 $m$ 行 $N$ 列的矩阵，高斯消元需要 $O(m \cdot N \cdot \text{nnz\_per\_row})$ 次域运算。实际中矩阵是稀疏的（$\text{nnz\_per\_row} \ll N$），加上 Gebauer–Moeller 筛选减少了临界对数量，F4 比朴素 Buchberger 算法快一到两个数量级。

### F5 算法：基于签名的 Gröbner 基

F5 算法（Faugère 2002）的革命性贡献在于：**在矩阵构造之前就识别并拒绝会产生零约化的 S-多项式**。这通过给每个多项式附加一个**签名**（signature）实现。

#### 签名与模表示

将理想生成问题转化为**模元素**问题。设 $I = \langle f_1, \ldots, f_s \rangle$。考虑自由模 $R^s$ 的标准基 $e_1, \ldots, e_s$。每个多项式 $h = \sum q_i f_i$ 对应一个模元素

$$\vec{h} = (q_1, \ldots, q_s) = \sum q_i e_i$$

$h$ 的**签名**定义为模元素中"最高"的单项式倍：

$$\text{sig}(h) = \max_\succ \{x^{\alpha_i} e_i : q_i \neq 0\}$$

其中比较使用 **pot**（position-over-term）序：先比较模块位置（$e_i$ 的索引 $i$，索引小者签名更小、更先被处理），再比较单项式。即 $x^\alpha e_i \succ x^\beta e_j$ 当且仅当 $i > j$，或 $i = j$ 且 $\alpha \succ \beta$。

**直觉**：签名记录了多项式 $h$ 的"来源"——它是由哪个生成元（$f_i$）的哪个倍（$x^\alpha$）开始推导的。多项式的签名是其推导历史的最简表示。

#### F5 准则

**签名判据**（Signature Criterion）。若一个候选 S-多项式的签名是某个**已知 syzygy** 的单项式倍，则该 S-多项式必化简为零——可直接拒绝，无需进入矩阵。

**syzygy 的含义**：当矩阵中某行化简为零时，该行对应多项式 $h = 0$，但其模表示 $(q_1, \ldots, q_s) \neq 0$。这给出一个 syzygy（关系）$\sum q_i f_i = 0$。该行的签名 $\text{sig}(h)$ 记录了这个 syzygy。

**SyzygySet 数据结构**：对每个模块位置 $k$，维护已知 syzygy 的首单项式集合。判断签名 $(k, t)$ 是否为 syzygy 时，检查是否有已存储的 $k$-位置 syzygy 的首单项式整除 $t$。

#### 正则序列与零约化的完全避免

**定理（Faugère 2002）**。若输入生成元 $f_1, \ldots, f_s$ 构成一个**正则序列**（regular sequence），则 F5 的签名判据足以拒绝**全部**零约化——即矩阵中不会出现化简为零的行。

正则序列的定义：$f_1, \ldots, f_s$ 是正则序列，如果 $f_1$ 不是零因子，且对 $i > 1$，$f_i$ 在 $R/\langle f_1, \ldots, f_{i-1} \rangle$ 中不是零因子。

实践中，即使生成元不是正则序列，F5 的签名判据也能过滤掉大量零约化。文献中在困难的 benchmark 理想上，F5 相比 F4 可获得数量级的加速（oCAS 的 F5 实现仍在持续优化中）。

#### F5 的矩阵构造

F5 的矩阵构造与 F4 类似，但每行额外携带签名信息：

1. **增量处理**：逐个添加生成元 $f_1, f_2, \ldots, f_s$
2. **带标签的 S-多项式**：每个 S-多项式都带着签名进入矩阵
3. **符号预处理**：与 F4 相同，但约化器的签名也需传播——当用基元素 $g$ 的倍 $x^\delta g$ 约化时，新行的签名为 $x^\delta \cdot \text{sig}(g)$
4. **syzygy 更新**：化简为零的行的签名加入 SyzygySet

**标签行（LabeledRow）**：矩阵中每行存储为 `(Signature, Vec<(coefficient, column)>)`——签名与稀疏行数据的配对。行阶梯运算时，签名不变（因为行运算是行之间的线性组合，签名取诸参与行中的最大者）。

#### F5 的复杂度优势

理论上，F5 的签名判据保证了：

- 对正则序列：矩阵中无零行，矩阵大小最小化
- 对一般理想：大幅减少零行，节省行阶梯运算开销

实践中的主要瓶颈是签名比较和 SyzygySet 查询。oCAS 的实现通过 SmallVec 和 HashMap 优化了这些操作。

### 单项式序的等价与换序

不同的计算任务需要不同的单项式序：

| 任务 | 需要的序 |
|------|---------|
| 消元 | Lex |
| 求解多项式系统 | Lex（结合 Sturm 根隔离） |
| 理想成员判定 | 任意序均可 |
| 高效计算 | Grevlex（中间多项式最小） |

标准策略是"先快后换"：

1. 用 **Grevlex** 计算基（最快——中间多项式次数最低）
2. 转换为 **Lex**（消元、求解需要）

两种转换方法：

**方法一：重新计算**（`reorder`）。将 Grevlex 基的多项式重新解释为 Lex 序，然后在 Lex 序下重跑 F4。对一般理想，这是默认方法——因为基的元素在新序下可能需要额外约化。

**方法二：FGLM**（零维理想专用）。FGLM 算法（下一小节详述）对零维理想提供 $O(n \cdot D^3)$ 的换序，远快于重新计算。

### FGLM 算法

FGLM（Faugère–Gianni–Lazard–Mora 1993）对**零维理想**提供高效的换序算法。

**前提条件**：理想 $I$ 是零维的，即 $R/I$ 作为 $k$-向量空间是有限维的。

#### 阶梯与正规形

**阶梯**（staircase）：不被任何首项单项式 $\text{lm}(g_i)$ 整除的单项式集合。对零维理想，阶梯是有限集，大小 $D = \dim_k(R/I)$。

阶梯中的单项式构成 $R/I$ 的一组 $k$-向量空间基。任何多项式 $f$ 的**正规形**（normal form）$f \bmod I$ 可以唯一表示为阶梯单项式的 $k$-线性组合。

#### 算法思路

FGLM 按目标序递增遍历单项式，计算其正规形，并检测线性关系：

```
输入: GB_src（源序下的约化 Gröbner 基）
输出: GB_tgt（目标序下的约化 Gröbner 基）

S ← {1}                    // 已处理的单项式（阶梯中）
B ← {x1, x2, ..., xn}     // 边界：已处理单项式的"邻居"
GB_tgt ← ∅
seen_nf ← {}               // 已见正规形

for m in 目标序递增:
    nf ← normal_form(m, GB_src)  // 在源序基下求正规形
    if nf ∈ span(seen_nf):        // 线性相关？
        coeffs ← find_relation(seen_nf, nf)
        new_poly ← m - Σ coeffs_i · m_i
        GB_tgt ← GB_tgt ∪ {new_poly}
    else:
        S ← S ∪ {m}
        seen_nf ← seen_nf ∪ {nf}
        B ← B ∪ {m · xi : i = 1..n}  // 扩展边界
    if |S| = D:                     // 阶梯满了
        break
return minimize(GB_tgt)
```

**关键性质**：

- 每个线性相关产生目标序基中的一个多项式
- 每个线性无关扩展阶梯，直到阶梯大小达到 $D$
- 总操作数 $O(n \cdot D^3)$——$D$ 次正规形计算（每次 $O(D^2)$），加上 $D$ 次高斯消元（每次 $O(D^2)$）

**非零维情况**：当阶梯无限时，FGLM 无法终止——`fglm` 返回 `None`，此时必须使用 `reorder`。

### Hilbert 函数与正则性

Hilbert 函数是描述分次环结构的强大工具，在 Gröbner 基计算中提供提前终止界。

#### Hilbert 函数与 Hilbert 级数

**定义**。设 $I \subseteq R = \mathbb{F}[x_1, \ldots, x_n]$ 是齐次理想。**Hilbert 函数**

$$H_{R/I}(d) = \dim_\mathbb{F} (R/I)_d$$

度量商环在次数 $d$ 处的齐次分量的向量空间维数。

**Hilbert 级数**是 Hilbert 函数的生成函数：

$$\text{HS}_{R/I}(t) = \sum_{d=0}^{\infty} H_{R/I}(d) \cdot t^d$$

**定理（Macaulay）**。$R/I$ 的 Hilbert 级数只依赖于首项理想 $\langle \text{lm}(g) : g \in \text{GB}(I) \rangle$，即

$$\text{HS}_{R/I}(t) = \text{HS}_{R/\langle\text{LM}(I)\rangle}(t)$$

这意味着我们可以仅从 Gröbner 基的首项计算 Hilbert 级数，无需知道理想的其他信息。

#### 单项式理想的 Hilbert 级数

对由单项式 $m_1, \ldots, m_s$ 生成的理想 $M$，容斥原理给出：

$$\text{HS}_{R/M}(t) = \frac{N(t)}{(1-t)^n}$$

其中分子 $N(t)$ 为

$$N(t) = \sum_{k=0}^{s} (-1)^k \sum_{\substack{S \subseteq \{1,\ldots,s\} \\ |S| = k}} t^{\deg \text{lcm}(S)}$$

这里 $\text{lcm}(\emptyset) = 1$（次数 0），$\text{lcm}(\{m_i\}) = m_i$。

**展开**：

$$N(t) = 1 - \sum_{i} t^{\deg m_i} + \sum_{i < j} t^{\deg \text{lcm}(m_i, m_j)} - \sum_{i < j < k} t^{\deg \text{lcm}(m_i, m_j, m_k)} + \cdots$$

**示例**。$M = \langle x^2, xy \rangle \subseteq \mathbb{F}[x, y]$（$n = 2$）。

- $|S| = 0$: $t^0 = 1$
- $|S| = 1$: $-t^2 - t^2 = -2t^2$（$\deg(x^2) = \deg(xy) = 2$）
- $|S| = 2$: $+t^{\deg \text{lcm}(x^2, xy)} = +t^{\deg(x^2 y)} = +t^3$

$$N(t) = 1 - 2t^2 + t^3$$

$$\text{HS}_{R/M}(t) = \frac{1 - 2t^2 + t^3}{(1-t)^2}$$

#### 正则性与提前终止

**正则性**（regularity）定义为 Hilbert 分子 $N(t)$ 的最高次非零项的次数：

$$\text{reg}(M) = \max\{d : [t^d] N(t) \neq 0\}$$

**定理（Bayer–Stillman）**。对零维理想，Gröbner 基中所有多项式的次数不超过 $\text{reg}(M)$。

这意味着 F4 可以安全地**忽略**次数超过正则性界的临界对——它们必然化简为零。oCAS 的 `hilbert` 模块通过 `regularity_bound` 计算此界，作为 F4 的提前终止提示（advisory bound，不影响正确性）。

#### 完整 Hilbert 级数

对一般（非齐次）理想的 Gröbner 基，Hilbert 级数也可计算：

$$\text{HS}_{R/I}(t) = \frac{N(t)}{(1-t)^n}$$

其中 $N(t)$ 由首项理想的容斥公式给出。oCAS 提供 `HilbertSeries` 结构体，支持：

- `hilbert_function(d)`：$\dim_k (R/I)_d$（通过展开 $\frac{N(t)}{(1-t)^n}$ 的 $t^d$ 系数）
- `dimension()`：$N(1) \neq 0$ 时维数为 $n - \text{ord}_{t=1}(N(t))$，即 $(1-t)$ 的消去次数
- `degree()`：射影簇的次数 $= h(1)$，其中 $h(t) = N(t)/(1-t)^{n-\dim}$ 是约去 $(1-t)^{n-\dim}$ 后的分子；对零维理想（$\dim = 0$）这等于 $\dim_k(R/I)$（点的个数）
- `hilbert_polynomial()`：对充分大的 $d$，$H_{R/I}(d)$ 是 $d$ 的多项式，通过 Lagrange 插值计算

---

## 在 oCAS 中的实现

oCAS 在 `ocas-poly/src/groebner/` 中实现三种 Gröbner 基算法和 Hilbert 级数计算。统一入口 `groebner_basis(ideal, Algorithm::Auto)` 目前默认路由到 F4。

### f4.rs：F4 矩阵行阶梯算法

#### Gebauer–Moeller 临界对管理

`update_pairs` 函数在每次新多项式加入基时执行 Gebauer–Moeller 筛选：

1. **生成候选**：新多项式与基中每个现有元素配对
2. **第一判据**（`monomial_are_coprime`）：若首项互素，直接跳过（乘积判据）
3. **第二判据**（链判据，顺序消去）：若另一个（仍存活）候选对的 $\text{lcm}$ 整除当前对的 $\text{lcm}$（整除关系含相等，重复对被丢弃），当前对冗余、不保留
4. **更新判据**（针对旧对）：仅当新首项整除旧对的 $\text{lcm}$，且 $\text{lcm}(i, \text{new})$、$\text{lcm}(j, \text{new})$ 都**严格小于** $\text{lcm}(i,j)$ 时，旧对才被移除；等号情形必须保留（否则会丢掉完整性所需的 S-多项式）

`CriticalPair` 结构体预计算了 $\text{lcm}$ 和单项式差（S-多项式的构造信息），避免在后续步骤中重复计算。

#### DivisorIndex：O(1)-ish 约化器查找

符号预处理的核心操作是"查找基中首项整除给定单项式的元素"。朴素实现需要 $O(\text{monomials} \times |\text{basis}|)$ 的线性扫描。oCAS 使用**支持位掩码**（support bitmask）索引加速此操作。

**support_mask**：对指数向量 $\alpha = (\alpha_1, \ldots, \alpha_n)$，定义位掩码

$$\text{mask}(\alpha) = \sum_{v : \alpha_v > 0} 2^v$$

即第 $v$ 位被设置当且仅当变量 $x_v$ 出现在单项式中（变量索引不超过 63）。

**关键观察**：若 $\text{lm}(g)$ 整除 $x^\alpha$，则 $\text{support}(\text{lm}(g)) \subseteq \text{support}(x^\alpha)$。因此查询时只需枚举 $\text{mask}(x^\alpha)$ 的**子掩码**，检查对应桶中的基元素。

**DivisorIndex** 结构：

```rust
struct DivisorIndex {
    // support_mask → 拥有该 mask 的基元素索引列表
    buckets: HashMap<u64, Vec<usize>>,
}
```

查询 `find_reducer(index, basis, exp)`：

1. 计算 $m = \text{mask}(\text{exp})$
2. 枚举 $m$ 的所有子掩码 $s \subseteq m$
3. 对 `buckets[s]` 中的每个基元素 $g$，精确检查 $\text{lm}(g) \mid x^\text{exp}$
4. 在所有满足条件的约化器中，选择**项数最少**的（加速约化）

子掩码枚举的复杂度 $O(3^v)$（$v$ 为变量数），但实际中 $v$ 通常很小（$\leq 10$），且桶内元素极少。

#### SimpCache：约化器倍的缓存

符号预处理中，同一基元素可能被多个单项式乘以不同的差来约化。`SimpCache<P>` 缓存已计算的多项式倍：

```rust
type SimpCache<P> = Vec<(SmallVec<[usize; 4]>, P)>;
```

`get_simplified(cache, diff, basis_poly)` 先查缓存，命中则直接返回；否则计算 $x^{\text{diff}} \cdot \text{basis\_poly}$ 并存入缓存。缓存跨轮次持久化。

#### ℤ_p 原生快速路径

当系数域为 `FiniteField(p)`（$p < 2^{31}$）时，F4 使用 `FpPoly` 替代泛型多项式：

```rust
struct FpPoly {
    terms: Vec<(i64, SmallVec<[usize; 4]>)>,  // (coeff, exponent)
    // coeff ∈ [0, p)，按降序排列
}
```

所有矩阵运算在 `i64` 上进行——无 BigInt 开销。`sub_scaled_fp` 通过**双指针归并**实现 O(nnz) 的稀疏行减法：

```rust
fn sub_scaled_fp(
    row: &mut Vec<(i64, usize)>,     // 列升序
    pivot: &[(i64, usize)],           // 列升序
    c: i64,                           // 乘数
    p: i64,                           // 模数
    scratch: &mut Vec<(i64, usize)>,  // 临时缓冲
) {
    // 两个指针归并：row -= c * pivot
    // 首列（主元列）自动消去（c = row[head_col] / pivot[head_col]）
    // 所有结果取 mod p，零系数丢弃
}
```

`register_row_fp` 进一步缓存 S-多项式和约化器倍的 `FpPoly` 表示，跨轮次复用。

#### 矩阵消元流程

完整 F4 流程：

1. **选择临界对**：从队列按 $\text{lcm}$ 大小排序选取
2. **符号预处理**：展开每个对的 S-多项式和所有约化器，建立列映射
3. **行阶梯化**：`echelonize_fp`（ℤ_p 快速路径）或 `echelonize_generic`（泛型）
4. **提取新元素**：非零行转回多项式，检查首项是否为"新的"
5. **更新基和对**：新元素加入基，`update_pairs` 更新临界对队列

### f5.rs：F5 签名算法

#### 签名结构

```rust
struct Signature {
    module_pos: usize,              // 生成元索引 k
    monomial: SmallVec<[usize; 4]>, // 单项式指数
}
```

签名按 **pot**（position-over-term）序比较：先比较 `module_pos`（小者优先），再按单项式序 `O` 比较 `monomial`。

#### LabeledPoly：带签名的基元素

```rust
struct LabeledPoly<D: Domain, O: MonomialOrder> {
    poly: SparseMultivariatePolynomial<D, O>,
    sig: Signature,
}
```

每个基元素携带签名。当用 $g$ 的倍 $x^\delta g$ 约化时，新多项式的签名为 $x^\delta \cdot \text{sig}(g)$。

#### SyzygySet：零约化追踪

```rust
struct SyzygySet {
    // module_pos → 该位置已知 syzygy 的首单项式列表
    lms: HashMap<usize, Vec<SmallVec<[usize; 4]>>>,
}
```

`insert(sig)`：当矩阵行化简为零时，其签名 $(k, t)$ 被记录——在 `lms[k]` 中存入 $t$。

`contains(sig)`：检查签名 $(k, t)$——若 `lms[k]` 中存在 $s$ 使 $s \mid t$，则为 syzygy。使用 `monomial_divides(t, s)` 做逐分量比较（$t_i \geq s_i$ 对所有 $i$）。

#### 矩阵构造与消元

F5 的 `build_and_reduce` 流程：

1. **选择对**：与 F4 类似的临界对选择（使用 `update_pairs` 共享实现）
2. **syzygy 过滤**：在将对加入矩阵前，检查其签名——若为 syzygy 则**跳过**
3. **带标签的行**：每行携带 `(Signature, Vec<(coefficient, column)>)`
4. **符号预处理**：与 F4 相同的列映射构建，但约化器的签名需传播
5. **行阶梯化**：按签名升序排列行，执行高斯消元
6. **新元素提取**：非零行提取为 `LabeledPoly`，签名保留
7. **syzygy 更新**：零行的签名加入 `SyzygySet`

#### ℤ_p 快速路径

F5 复用 F4 的 `FpPoly`、`CriticalPair`、`update_pairs` 等基础设施，通过 `LabeledFpPoly` 和 `LabeledFpRow` 添加签名标签。`f5_fp` 的结构与泛型 `f5` 完全对称——仅多项式表示从泛型域切换为 `i64` 模算术。

### hilbert.rs：Hilbert 级数计算

#### 容斥原理实现

`hilbert_numerator(generators)` 计算单项式理想 $\langle m_1, \ldots, m_s \rangle$ 的 Hilbert 分子：

```rust
pub fn hilbert_numerator(generators: &[Vec<usize>]) -> Vec<(usize, i64)> {
    // 对所有子集 S ⊆ {0, ..., s-1} 遍历
    // |S| = k 贡献 (-1)^k · t^{deg lcm(S)}
    // 返回 (degree, coefficient) 对的稀疏表示
}
```

枚举所有 $2^s$ 个子集（实用上限约 $s \leq 20$），对每个子集计算 $\text{lcm}$ 的次数（逐分量取最大值后求和），按符号累加到对应次数。

#### 正则性界

`regularity_bound(generators)` 返回 Hilbert 分子的最高次非零项次数：

$$\text{reg} = \max\{d : [t^d] N(t) \neq 0\}$$

F4 将此界用作**提前终止提示**——选择临界对时，若 $\deg(\text{lcm}) > \text{reg}$，可安全跳过。这不影响正确性（只影响效率），且仅对零维理想有效。

#### 完整 Hilbert 级数

`hilbert_series(gb)` 从 Gröbner 基计算 $R/I$ 的完整 Hilbert 级数 $\frac{N(t)}{(1-t)^n}$：

1. **提取首项**：从基中取出每个多项式的首项指数向量
2. **容斥计算**：调用 `hilbert_numerator` 得到 $N(t)$
3. **展开级数**：$\frac{N(t)}{(1-t)^n} = N(t) \cdot \sum_{k \geq 0} \binom{n+k-1}{k} t^k$

`HilbertSeries` 结构体的方法：

| 方法 | 计算方式 |
|------|---------|
| `hilbert_function(d)` | 展开 $\frac{N(t)}{(1-t)^n}$ 的 $t^d$ 系数 |
| `dimension()` | $N(1) \neq 0$ 时为 $n$，否则递归检查 $(1-t)$ 的消去 |
| `degree()` | 对 $N(t)$ 求 $(n-\dim)$ 阶导数并在 $t=1$ 求值，再除以 $(n-\dim)!$；零维（$\dim=0$）时等于 $\dim_k(R/I) = \text{HS}(1)$（点的个数） |
| `hilbert_polynomial()` | Lagrange 插值：在足够多的 $d$ 处计算 $H(d)$，拟合多项式 |

`binomial_general(n, k)` 计算广义二项式系数 $\binom{n+k-1}{k}$（$n$ 可为负数或零），用于展开 $(1-t)^{-n}$。

### fglm.rs：零维换序算法

#### 阶梯计算

`compute_staircase(lms, n_vars)` 枚举不被任何首项单项式整除的单项式：

1. 从 $1 = x_1^0 \cdots x_n^0$ 开始 BFS
2. 对每个单项式 $m$，检查是否被某个 $\text{lm}(g_i)$ 整除
3. 不被整除的加入阶梯，其所有"邻居"（$m \cdot x_i$）加入队列
4. 若阶梯大小超过阈值仍未停止，返回 `None`（正维理想）

#### 正规形计算

`normal_form_monomial(m, gb, staircase, domain)` 在源序基下计算单项式 $m$ 的正规形，结果表示为阶梯坐标的 $k$-向量：

$$m \bmod I = \sum_{i=1}^{D} c_i \cdot s_i$$

其中 $s_1, \ldots, s_D$ 是阶梯单项式。计算方式是标准的多元除法算法。

#### 线性相关检测

`find_relation(seen, nf, domain)` 检测新正规形 `nf` 是否在已见正规形的张成空间中。使用**增量高斯消元**：维护一个阶梯形矩阵，每次新向量尝试化简——化简为零则找到关系。

#### 边界管理

`mark_multiples(visited, m, n_vars, max_deg)` 标记单项式 $m$ 的所有倍（至全次数 `max_deg`）为已访问，避免重复插入边界。使用递归枚举。

### 统一入口与算法选择

```rust
pub fn groebner_basis<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    algo: Algorithm,
) -> GroebnerBasis<D, O>
```

`Algorithm` 枚举：

| 变体 | 行为 |
|------|------|
| `Auto` | 当前路由到 F4；未来将基于 cyclic-n 基准校准 F5 交叉点 |
| `Buchberger` | 经典 S-多项式迭代（教学用） |
| `F4` | 矩阵行阶梯批处理 |
| `F5` | 签名判据 + 矩阵构造 |

所有路径最终输出**约化 Gröbner 基**：先 `minimize`（移除被其他元素首项整除的冗余元素），再 `auto_reduce`（互相约化尾项）。

---

## 进阶话题

### F4 vs Buchberger 的实际性能

oCAS 的基准（`ocas-tests/benches/groebner.rs`）在 cyclic-n 理想上比较 Buchberger 与 F4：

| 基准 | Buchberger | F4 | 说明 |
|------|-----|-----|------|
| cyclic-3（$\mathbb{Q}$） | 可行 | 快 | Buchberger 在此规模已很慢 |
| cyclic-4（$\mathbb{Q}$） | 过慢 | 可行 | 只有 F4 可行 |
| cyclic-3/4/5（$\mathbb{Z}_{13}$、$\mathbb{Z}_{101}$） | — | 快 | i64 快速路径可延伸到 cyclic-5 |

F5 的签名判据在理论上能避免零约化，对困难理想有数量级加速的潜力，但 oCAS 的 F5 实现较新，尚无正式基准数据。

**实践建议**：对一般理想使用 `Algorithm::Auto`（当前 = F4）；对已知为正则序列的系统，显式指定 `Algorithm::F5` 可能更快。

### F5 的理论保证

**定理**。F5 算法在以下意义下是"最优"的：

1. **无冗余约化**：对正则序列，矩阵中无化简为零的行
2. **完备性**：产生与 F4 完全相同的约化 Gröbner 基
3. **签名唯一性**：每个多项式的签名唯一确定其来源

这些保证使 F5 成为处理困难理想的理论首选，尽管实际中 F4 的简洁实现常在中小规模问题上更快。

### Hilbert 级数的应用

Hilbert 级数不仅是理论工具，在 oCAS 中有以下实际应用：

1. **提前终止**：正则性界 $\text{reg}(M)$ 为 F4 提供可靠的次数上界
2. **维数判定**：$R/I$ 的 Krull 维数 = $\text{ord}_{t=1} \text{HS}(t)$ 的阶数
3. **次数计算**：射影簇的次数 = $N(1)$（零维时为阶梯大小）
4. **零维判定**：$N(1) \neq 0$ 当且仅当理想为零维

### 计算代数几何中的 Gröbner 基

Gröbner 基是计算代数几何的通用工具。在 oCAS 中，它们支撑着以下应用（详见 [FGLM 与消元理论](./fglm-elimination.md) 和 [Gröbner 基实现](../algorithms/groebner.md)）：

- **多项式方程组求解**：Gröbner 基 → 消元 → 一元方程 → 根隔离
- **理想运算**：成员判定、和、积、商、饱和、交集、消元、根式、准素分解
- **代数簇的几何性质**：维数、次数、不可约分量

---

## 参考文献

1. **Cox, D., Little, J. & O'Shea, D.** *Ideals, Varieties, and Algorithms.* 4th ed., Springer, 2015.
   - 第 2 章：Gröbner 基的定义、Buchberger 算法、Buchberger 判据
   - 第 3 章：消元定理与理想运算
   - 第 8 章：FGLM 与换序

2. **Adams, W.W. & Loustaunau, P.** *An Introduction to Gröbner Bases.* AMS, 1994.
   - 第 4–5 章：Buchberger 算法的详细分析、Gebauer–Moeller 优化

3. **Faugère, J.-C.** "A New Efficient Algorithm for Computing Gröbner Bases (F4)." *Journal of Pure and Applied Algebra*, 139(1–3):61–88, 1999.
   - F4 算法的原始论文：矩阵行阶梯批处理

4. **Faugère, J.-C.** "A New Efficient Algorithm for Computing Gröbner Bases without Reduction to Zero (F5)." *ISSAC 2002*, pp. 75–83.
   - F5 算法的原始论文：签名判据与正则序列

5. **Faugère, J.-C., Gianni, P., Lazard, D. & Mora, T.** "Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering." *Journal of Symbolic Computation*, 16(4):329–344, 1993.
   - FGLM 算法的原始论文

6. **Gebauer, R. & Möller, H.M.** "On an Installation of Buchberger's Algorithm." *Journal of Symbolic Computation*, 6(2–3):275–286, 1988.
   - Gebauer–Moeller 临界对筛选准则

7. **Eder, C. & Perry, J.** "Signature-based Algorithms to Compute Gröbner Bases." *ISSAC 2009*, pp. 139–146.
   - Eder–Perry 的 F5 形式化

8. **Bayer, D. & Stillman, M.** "A Criterion for Detecting m-Regularity." *Inventiones Mathematicae*, 87(1):1–11, 1987.
   - 正则性界与 Hilbert 函数的关系

9. **Gathen, J. von zur & Gerhard, J.** *Modern Computer Algebra.* 3rd ed., Cambridge University Press, 2013.
   - 第 21 章：Gröbner 基的全面综述

**参见**：

- [多项式代数](./polynomial-algebra.md) — 单项式序、多元除法算法的基础
- [线性代数](./linear-algebra.md) — Bareiss 算法与高斯消元
- [FGLM 与消元理论](./fglm-elimination.md) — FGLM 详细分析与理想运算
- [Gröbner 基实现](../algorithms/groebner.md) — 使用指南、基准与理想运算 API
- [Rust API：Gröbner 基与理想](../api/rust-groebner.md) — 函数签名与完整示例
