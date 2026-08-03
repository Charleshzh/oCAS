# 基础：线性代数

本章介绍 oCAS 中矩阵运算的数学基础。oCAS 的 `Matrix` 类型工作在
`EuclideanDomain` 上——不只是浮点数或有理数，而是任意满足欧几里得除法
的代数结构。这意味着我们需要一种**不产生分数**的消元算法。

---

## 前提知识

在阅读本章之前，你需要了解：

- 环与域的基本概念（参见[有限域与模算术](./finite-fields.md)）
- 欧几里得域：支持 `add`、`mul`、`sub`、`neg`、`div`（带余除法）、`gcd`、`rem` 的代数结构
- 向量空间的基本直觉：向量加法与标量乘法

不需要线性代数的系统训练——本章从零开始构建。

---

## 基础概念

### 矩阵

一个 $m \times n$ **矩阵**是一个按 $m$ 行 $n$ 列排列的元素矩形阵列：

$$
A = \begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

元素 $a_{ij}$ 位于第 $i$ 行第 $j$ 列。在 oCAS 中，元素存储在行优先的
一维数组中：`data[i * ncols + j]` 对应 $a_{i+1,\,j+1}$（零索引）。

### 矩阵运算

**加法**：两个同形矩阵逐元素相加：

$$
(A + B)_{ij} = a_{ij} + b_{ij}
$$

**标量乘法**：

$$
(cA)_{ij} = c \cdot a_{ij}
$$

**矩阵乘法**：$A$（$m \times n$）与 $B$（$n \times p$）的积 $C = AB$
是 $m \times p$ 矩阵：

$$
c_{ij} = \sum_{k=1}^{n} a_{ik} \cdot b_{kj}
$$

乘法要求 $A$ 的列数等于 $B$ 的行数。注意矩阵乘法一般**不满足交换律**
（$AB \neq BA$），但满足**结合律**（$(AB)C = A(BC)$）。

**转置**：交换行与列：

$$
(A^\top)_{ij} = a_{ji}
$$

性质：$(AB)^\top = B^\top A^\top$，$(A^\top)^\top = A$。

**迹**（trace）：方阵对角元素之和：

$$
\operatorname{tr}(A) = \sum_{i=1}^{n} a_{ii}
$$

性质：$\operatorname{tr}(AB) = \operatorname{tr}(BA)$（当乘积有定义时）。

### 行列式

**行列式**是方阵的一个标量值，是线性代数中最核心的概念之一。对
$n \times n$ 矩阵 $A$，行列式 $\det(A)$ 有以下等价定义：

**递归定义（Laplace 展开）**：沿第 $i$ 行展开：

$$
\det(A) = \sum_{j=1}^{n} (-1)^{i+j} \, a_{ij} \, M_{ij}
$$

其中 $M_{ij}$ 是删去第 $i$ 行第 $j$ 列后的 $(n-1)\times(n-1)$ 子矩阵
的行列式（**余子式**）。对 $1\times 1$ 矩阵，$\det(a_{11}) = a_{11}$。

**几何定义**：$|\det(A)|$ 是矩阵列向量所张成的平行多面体的"有向体积"。
$\det(A) > 0$ 表示保持定向，$\det(A) < 0$ 表示翻转定向。

#### 行列式的关键性质

1. **乘法性**：$\det(AB) = \det(A) \cdot \det(B)$
2. **转置不变性**：$\det(A^\top) = \det(A)$
3. **行交换变号**：交换两行，行列式变号
4. **行倍乘**：某行乘以 $c$，行列式乘以 $c$
5. **行加法不变**：将一行的倍数加到另一行，行列式不变
6. **奇异判定**：$\det(A) = 0$ 当且仅当 $A$ 不可逆（奇异）
7. **上/下三角矩阵**：$\det(A) = \prod_{i=1}^{n} a_{ii}$（对角元素之积）

这些性质直接引出了用**消元法**计算行列式的思路：通过行变换将矩阵化为
三角矩阵，行列式就是对角元素之积（乘以行交换的符号因子）。

---

## 核心理论

### 问题：域上的分数膨胀

经典的高斯消元在有理数域上工作良好——我们可以自由地进行除法。但 oCAS 的
`Matrix<D: EuclideanDomain>` 要求在**任意**欧几里得域上工作，包括：

- **整数环 $\mathbb{Z}$**：$3 \div 2$ 无意义（不存在整数商）
- **多项式环 $\mathbb{F}[x]$**：$x^2 \div (x+1)$ 不是多项式
- **有限域 $\mathbb{F}_p$**：可以除，但每次除法产生模逆运算，且中间结果可能溢出

即使在有理数上，经典消元也会产生越来越大的分母。例如对整数矩阵：

$$
A = \begin{pmatrix} 2 & 3 \\ 5 & 7 \end{pmatrix}
$$

经典消元需要计算 $7 - \frac{5 \cdot 3}{2} = 7 - \frac{15}{2} = -\frac{1}{2}$——
出现了分数。对 $n \times n$ 矩阵，分数的分子分母可能增长到 $O(n!)$ 量级。

### Bareiss 无分数行列式算法

**Bareiss 算法**（1968）的核心洞察：在整数（或更一般的欧几里得域）上，
消元过程中的所有中间值都可以通过除以前一个主元来保持为整数。

#### 算法描述

对 $n \times n$ 矩阵 $A$，Bareiss 算法产生上三角矩阵，其最后一个对角
元素 $a_{nn}^{(n-1)}$ 就是 $\det(A)$（可能需要乘以 $-1$ 来补偿行交换）。

初始化：令 $a_{ij}^{(0)} = a_{ij}$，$p_0 = 1$。

对 $k = 0, 1, \ldots, n-2$：

1. **选主元**：若 $a_{kk}^{(k)} = 0$，寻找 $i > k$ 使 $a_{ik}^{(k)} \neq 0$，
   交换第 $k$ 行与第 $i$ 行（记录符号翻转）。若找不到，$\det(A) = 0$。

2. **消元**：对 $i = k+1, \ldots, n-1$ 和 $j = k+1, \ldots, n-1$：

$$
a_{ij}^{(k+1)} = \frac{a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}}{p_k}
$$

3. **更新除数**：$p_{k+1} = a_{kk}^{(k)}$

最终：

$$
\det(A) = \begin{cases}
a_{nn}^{(n-1)} & \text{若无行交换或偶数次交换} \\
-a_{nn}^{(n-1)} & \text{若有奇数次行交换}
\end{cases}
$$

#### 为什么除法是精确的

**定理**（Bareiss 1968）：上述公式中，$p_k$ 精确整除分子
$a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}$。

直觉上，这来自 Sylvester 行列式恒等式。消元矩阵的 $(k+1)$ 步更新
等价于计算一个 $(k+2) \times (k+2)$ 子矩阵的行列式，而 Schur 补公式
保证除以前一个主元的结果仍是整数。

#### 复杂度与数值行为

- **时间复杂度**：$O(n^3)$ 次域运算——与经典高斯消元相同
- **中间值大小**：对整数矩阵，中间值 $a_{ij}^{(k)}$ 的位数增长为
  $O(n \log n)$——远优于经典消元的 $O(n!)$ 量级增长
- **无舍入误差**：在整数和多项式上，所有运算精确，无浮点近似

### 高斯消元与行阶梯形

高斯消元将矩阵化为**行阶梯形**（Row Echelon Form, REF）：

$$
\begin{pmatrix}
\boxed{a_{11}} & a_{12} & a_{13} & a_{14} \\
0 & \boxed{a_{22}} & a_{23} & a_{24} \\
0 & 0 & 0 & \boxed{a_{34}} \\
0 & 0 & 0 & 0
\end{pmatrix}
$$

带框的元素是**主元**（pivots）。行阶梯形满足：
1. 所有非零行在所有零行之上
2. 每个主元所在列在其下方行中全为零
3. 每个主元严格在其上方行的主元右侧

继续消元可得**简化行阶梯形**（Reduced Row Echelon Form, RREF）——
每个主元为 1，且主元所在列的其余元素全为 0。

#### 主元选取

欧几里得域上未必有"绝对值"概念（例如 $\mathbb{F}[x]$），因此 oCAS 的消元不使用数值线性代数中的**部分主元选取**（partial pivoting）。`row_echelon` 与 `determinant` 都选择当前列中**第一个非零元素**作为主元；若当前列全为零则跳过该列（不交换行）。中间值的膨胀由**内容剥离**（content stripping，见下文）控制——消元后立即提取主元行的 GCD 并整行除以它。

### 分数无关回代

消元后需要**回代**（back-substitution）来求解线性系统。oCAS 的回代也是
分数无关的：从最后一行开始向上，用 GCD 缩放避免引入分数。

对行阶梯形中的第 $i$ 行（主元在第 $j$ 列）和第 $k < i$ 行（第 $j$ 列
有非零元素 $a_{kj}$）：

$$
g = \gcd(a_{ij}, a_{kj}), \quad s_p = \frac{a_{kj}}{g}, \quad s_r = \frac{a_{ij}}{g}
$$

更新第 $k$ 行：对每个 $l > j$：

$$
a_{kl} \leftarrow s_r \cdot a_{kl} - s_p \cdot a_{il}
$$

$$
a_{kj} \leftarrow 0
$$

这与 Bareiss 消元的无分数性质完全一致。

### 内容剥离

消元过程中，行的系数可能包含公共因子。oCAS 在每一步主元处理后执行
**内容剥离**（content stripping）：计算主元行从主元列到最后一列的
GCD，然后整行除以该 GCD。

$$
g = \gcd(a_{ij}, a_{i,j+1}, \ldots, a_{i,n-1})
$$

$$
a_{il} \leftarrow \frac{a_{il}}{g}, \quad l = j, j+1, \ldots, n-1
$$

这个优化对于多项式环尤其关键——否则中间多项式的次数会指数增长。

### 逆矩阵

方阵 $A$ 的**逆矩阵** $A^{-1}$ 满足：

$$
AA^{-1} = A^{-1}A = I
$$

其中 $I$ 是单位矩阵。逆存在的充要条件是 $\det(A) \neq 0$（$A$ 非奇异）。

计算方法：将 $[A \mid I]$ 做行变换化为 $[I \mid A^{-1}]$。在 oCAS 中，
`inverse()` 通过对单位矩阵的每一列 $e_j$ 求解 $Ax_j = e_j$ 来得到
$A^{-1}$ 的第 $j$ 列。

---

## 在 oCAS 中的实现

### `Matrix<D: EuclideanDomain>` 结构

oCAS 的矩阵是**泛型**的——它接受任何 `EuclideanDomain` 作为系数域：

```rust
use ocas_domain::{EuclideanDomain, IntegerDomain, Integer};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::new(2, 2, vec![
    Integer::from(2), Integer::from(3),
    Integer::from(5), Integer::from(7),
], d);
```

内部存储为行优先的 `Vec<D::Element>`，加上行数、列数和系数域实例
（`domain: D` 按值存储，而非引用）。这使得所有算术操作通过域的 trait
方法（`add`、`mul`、`sub`、`div`、`gcd`）进行，完全不依赖具体类型。

### 为什么选择 Bareiss 而非 Laplace 展开

| 方法 | 时间复杂度 | 空间复杂度 | 分数膨胀 | 域兼容性 |
|---|---|---|---|---|
| Laplace 展开 | $O(n!)$ | $O(n^2)$（递归深度 $O(n)$） | 无 | 任意交换环（实际不可行） |
| 经典高斯消元 | $O(n^3)$ | $O(n^2)$ | 严重 | 仅域 |
| **Bareiss** | $O(n^3)$ | $O(n^2)$ | **无** | **欧几里得域** |

- **Laplace 展开**的 $O(n!)$ 复杂度使其对 $n > 10$ 完全不可用
- **经典消元**需要分数运算，不适用于 $\mathbb{Z}$ 或 $\mathbb{F}[x]$
- **Bareiss 算法**兼得两者之长：$O(n^3)$ 效率 + 无分数 + 支持任意欧几里得域

这是 oCAS 选择 Bareiss 的根本原因——它是一个**通用符号计算系统**，
矩阵元素可能是整数、多项式、有限域元素，甚至代数数域元素。

### 行阶梯化：`row_echelon`

```rust
// 分数无关高斯消元，返回秩
pub fn row_echelon(&mut self, max_col: usize) -> usize
```

`row_echelon` 对前 `max_col` 列执行分数无关消元，返回矩阵的秩。算法步骤：

1. 对每列 $j$，选择当前列中第一个非零元素作为主元（若全列为零则跳过该列）
2. **内容剥离**：提取主元行的 GCD 并整行除以它
3. **消元**：对主元下方的每行，用 GCD 缩放后的减法消除第 $j$ 列（$a_{ij}$ 为主元，$a_{kj}$ 为被消元元素，$s_r$、$s_p$ 须在清零前计算）：

$$
g = \gcd(a_{ij}, a_{kj}), \quad s_r = \frac{a_{ij}}{g}, \quad s_p = \frac{a_{kj}}{g}
$$

$$
a_{kj} \leftarrow 0, \qquad a_{kl} \leftarrow s_r \cdot a_{kl} - s_p \cdot a_{il} \quad (l > j)
$$

这保证所有中间值保持在域内，无分数产生。

### 行列式：`determinant`

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(2)],
    vec![Integer::from(3), Integer::from(4)],
], d);
assert_eq!(a.determinant().unwrap(), Integer::from(-2));
```

`determinant` 内部直接实现 Bareiss 算法：

1. 初始化：复制数据到工作数组，$p_0 = 1$，符号标记 $\epsilon = 1$
2. 对 $k = 0, \ldots, n-2$：
   - 若主元为零，寻找下方非零行并交换（翻转 $\epsilon$）
   - 更新：$a_{ij}^{(k+1)} = (a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}) / p_k$
   - 更新除数：$p_{k+1} = a_{kk}^{(k)}$
3. 返回 $\epsilon \cdot a_{n-1,n-1}^{(n-1)}$

特殊情况：$n = 0$ 返回 $1$（空矩阵行列式），$n = 1$ 直接返回唯一元素。

### 线性系统求解：`solve`

```rust
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(1)],
    vec![Integer::from(1), Integer::from(-1)],
], d);
let b = vec![Integer::from(3), Integer::from(-1)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2)]);
```

`solve` 的完整流程：

1. **增广**：将 $b$ 拼接为 $[A \mid b]$
2. **消元**：`row_echelon(nvars)` 化为行阶梯形
3. **一致性检查**：若主元行之后有非零的 $b$ 分量，系统矛盾，返回
   `MatrixError::Inconsistent`
4. **回代**：`back_substitution(nvars)` 分数无关回代
5. **秩检查**：若秩 $< n$，系统欠定，返回 `MatrixError::Underdetermined`
6. **最终除法**：每个解分量除以其对应主元

#### 错误类型

| 错误 | 含义 |
|---|---|
| `ShapeMismatch` | $A$ 不是方阵，或 $b$ 的长度不等于行数 |
| `RightHandSideIsNotVector` | $b$ 不是列向量 |
| `Inconsistent` | 系统无解（增广矩阵秩 > 系数矩阵秩） |
| `Underdetermined { rank }` | 系统有无穷多解（秩 < 未知数个数） |
| `ResultNotInDomain` | 解不在域中（整数矩阵可能产生有理数解） |

### 秩与逆矩阵

```rust
// 秩：行阶梯形的主元个数
let rank = a.rank();

// 逆矩阵：逐列求解 Ax_j = e_j
let inv = a.inverse()?;
```

`rank` 将矩阵克隆后调用 `row_echelon`，返回主元个数。

`inverse` 对单位矩阵的每一列求解线性系统，得到逆矩阵的各列。这要求
矩阵非奇异且解在域内——例如整数矩阵的逆一般不在 $\mathbb{Z}$ 中，
除非 $\det(A) = \pm 1$（幺模矩阵）：

```rust
// 幺模矩阵：det = 1，整数逆存在
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(2)],
    vec![Integer::from(0), Integer::from(1)],
], d);
let inv = a.inverse().unwrap();
assert_eq!(inv[(0, 0)], Integer::from(1));
assert_eq!(inv[(0, 1)], Integer::from(-2));
assert_eq!(inv[(1, 0)], Integer::from(0));
assert_eq!(inv[(1, 1)], Integer::from(1));
```

### 矩阵乘法与转置

```rust
// 乘法：C = A · B
let c = a.matmul(&b)?;

// 转置
let at = a.transpose();

// 迹（方阵对角元素之和）
let tr = a.trace()?;
```

`matmul` 使用标准的三重循环 $O(n^3)$ 算法，所有运算通过域的 trait 方法。
`transpose` 通过重新排列数据实现。`trace` 对方阵求对角元素之和。

### 完整示例：求解整数线性系统

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 方程组：
//   2x + 3y = 8
//   5x + 7y = 19
let a = Matrix::from_rows(vec![
    vec![Integer::from(2), Integer::from(3)],
    vec![Integer::from(5), Integer::from(7)],
], d);
let b = vec![Integer::from(8), Integer::from(19)];

let x = a.solve(&b).unwrap();
// x = [1, 2]
assert_eq!(x[0], Integer::from(1));
assert_eq!(x[1], Integer::from(2));

// 验证：Ax = b
let x_matrix = Matrix::new(2, 1, x, d);
let result = a.matmul(&x_matrix).unwrap();
assert_eq!(result[(0, 0)], Integer::from(8));
assert_eq!(result[(1, 0)], Integer::from(19));

// 行列式：2·7 - 3·5 = -1
assert_eq!(a.determinant().unwrap(), Integer::from(-1));
```

消元过程中的中间值全部保持为整数——Bareiss 算法保证了这一点。

---

## 进阶话题

### Bareiss 算法的推广

Bareiss 算法可以推广到计算**子矩阵的行列式**（所有 $k \times k$ 阶子式），
这在计算 Smith 标准形和 Hermite 标准形时很有用。

对多项式矩阵，Bareiss 通过 GCD 缩放有效地控制了系数膨胀，使得
多项式矩阵的行列式计算在实践中可行。

### 与 LU 分解的关系

Bareiss 算法可以看作是**无分数 LU 分解**的变体。经典 LU 分解将矩阵
分解为 $A = LU$，其中 $L$ 是下三角、$U$ 是上三角。Bareiss 产生
的是类似的分解，但所有中间值保持在域内。

$\det(A) = \det(L) \cdot \det(U) = \prod_{i} l_{ii} \cdot \prod_{i} u_{ii}$。

### 符号计算中的矩阵

在符号计算中，矩阵元素可能是**多项式**或**有理函数**。此时 Bareiss
算法的 GCD 缩放尤为关键——它将多项式系数的次数增长控制在 $O(n)$ 而非
指数增长。

oCAS 的 Gröbner 基 F4 算法内部大量使用矩阵行阶梯化——将多元多项式的
消元问题转化为系数域上的矩阵问题，再用 Bareiss 或其 `ℤ_p` 快速路径
（惰性模算术）来高效处理。

---

## 参考文献

- **Axler, S.** *Linear Algebra Done Right*, 3rd ed. Springer, 2015.
  Chapters 4–6: 行列式、特征值与内积空间。

- **Bareiss, E. H.** "Sylvester's Identity and Multistep
  Integer-Preserving Gaussian Elimination." *Mathematics of Computation*,
  22(103):565–578, 1968.
  Bareiss 算法的原始论文。

- **Cox, D., Little, J., O'Shea, D.** *Ideals, Varieties, and Algorithms*,
  4th ed. Springer, 2015.
  多项式矩阵与 Gröbner 基中的线性代数。

- **Geddes, K. O., Czapor, S. R., Labahn, G.** *Algorithms for Computer
  Algebra*. Kluwer, 1992.
  第 5 章：矩阵上的符号消元算法。
