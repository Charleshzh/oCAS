# Rust API 参考：多项式

本章记录 oCAS 的多项式系统，涵盖三类核心数据结构：

- **`DenseUnivariatePolynomial<D>`** — 稠密一元多项式
- **`SparseMultivariatePolynomial<D, O>`** — 稀疏多元多项式
- **`RationalPolynomial<D, O>`** — 多项式分式域元素

以及单项式序 trait `MonomialOrder` 和多种序的实现。

**模块路径**：`ocas_poly`

---

## 目录

- [单项式序](#单项式序)
  - [MonomialOrder trait](#monomialorder-trait)
  - [Lex](#lex)
  - [Grlex](#grlex)
  - [Grevlex](#grevlex)
  - [WeightOrder](#weightorder)
  - [BlockOrder 与 SubOrder](#blockorder-与-suborder)
  - [MatrixOrder](#matrixorder)
- [DenseUnivariatePolynomial](#denseunivariatepolynomial)
  - [构造与属性](#dense-构造与属性)
  - [基本算术](#dense-基本算术)
  - [EuclideanDomain 操作](#euclideandomain-操作)
  - [因式分解与结式](#因式分解与结式)
- [SparseMultivariatePolynomial](#sparsemultivariatepolynomial)
  - [构造与属性](#sparse-构造与属性)
  - [基本算术](#sparse-基本算术)
  - [Gröbner 基支持](#gröbner-基支持)
  - [多元因式分解](#多元因式分解)
- [RationalPolynomial](#rationalpolynomial)
- [辅助函数](#辅助函数)

---

## 单项式序

### MonomialOrder trait

**签名**：

```rust
pub trait MonomialOrder: Clone + PartialEq + Eq + std::fmt::Debug + Default {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering;
}
```

**功能**：定义单项式之间的全序关系。多项式的首项、排序和 Gröbner 基计算都依赖于此。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `self` | `&Self` | 序的实例（零大小类型或带参数的类型） |
| `lhs` | `&[usize]` | 左侧单项式的指数向量 $[\alpha_1, \alpha_2, \dots]$ |
| `rhs` | `&[usize]` | 右侧单项式的指数向量 |

**返回值**：`std::cmp::Ordering` — `Less` 表示 `lhs` 在序中排在 `rhs` **之前**（即 `lhs` 更大），`Greater` 表示 `lhs` 更小。

**设计说明**：简单序（Lex、Grevlex、Grlex）是零大小类型，无运行时开销。参数化序（WeightOrder、BlockOrder、MatrixOrder）在构造时存储配置。

**参见**：[Lex](#lex)、[Grevlex](#grevlex)、[Grlex](#grlex)、[WeightOrder](#weightorder)、[BlockOrder](#blockorder-与-suborder)、[MatrixOrder](#matrixorder)

---

### Lex

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Lex;
```

**功能**：字典序。从左到右逐分量比较指数向量。

**比较规则**：`lhs > rhs` 当且仅当存在 $i$ 使得 $\alpha_j = \beta_j$ 对所有 $j < i$ 且 $\alpha_i > \beta_i$。

**示例**：

```rust
use ocas_poly::sparse::{Lex, MonomialOrder};

let a = [2, 1]; // x^2 y
let b = [1, 1]; // x y
assert_eq!(Lex.cmp(&a, &b), std::cmp::Ordering::Greater);
// a 在 Lex 序中更大（首分量 2 > 1）
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

### Grlex

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grlex;
```

**功能**：分次字典序。先按总次数降序排列，同次再按字典序。

**比较规则**：
1. 总次数 $\sum \alpha_i$ 更大者排在前面
2. 总次数相同时，按字典序比较

**示例**：

```rust
use ocas_poly::sparse::{Grlex, MonomialOrder};

let a = [2, 0]; // x^2, 次数 2
let b = [1, 1]; // xy, 次数 2
let c = [0, 3]; // y^3, 次数 3
// c 次数最高，排在最前
assert_eq!(Grlex.cmp(&c, &a), std::cmp::Ordering::Less);
// a 和 b 次数相同，按字典序 a > b
assert_eq!(Grlex.cmp(&a, &b), std::cmp::Ordering::Greater);
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

### Grevlex

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grevlex;
```

**功能**：分次反字典序。先按总次数降序排列，同次再按**反向**字典序比较（从最后一个分量开始，方向反转）。

**比较规则**：
1. 总次数更大者排在前面
2. 总次数相同时，从最后一个分量开始比较，**更小**者排在前面

**特点**：Grevlex 序在 Gröbner 基计算中通常产生最小的中间矩阵，是默认序。

**示例**：

```rust
use ocas_poly::sparse::{Grevlex, Lex, MonomialOrder};

let a = [2, 1];
let b = [1, 1];
assert_eq!(Lex.cmp(&a, &b), std::cmp::Ordering::Greater);
assert_eq!(Grevlex.cmp(&a, &b), std::cmp::Ordering::Less);
// Grevlex 下 a < b（反向字典序反转了方向）
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

### WeightOrder

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightOrder {
    weights: SmallVec<[i64; 4]>,
}
```

**功能**：加权序。按 $\sum_i w_i \cdot \alpha_i$ 降序排列。适用于无法用零大小类型表达的消元序。

**构造方法**：

```rust
pub fn new(weights: SmallVec<[i64; 4]>) -> Self
pub fn from_slice(weights: &[i64]) -> Self
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `weights` | `SmallVec<[i64; 4]>` 或 `&[i64]` | 每个变量的权重，长度应等于变量数 |

**Default**：全 1 权重（即总次数序）。

**示例**：

```rust
use ocas_poly::sparse::{MonomialOrder, WeightOrder};
use smallvec::smallvec;

let ord = WeightOrder::new(smallvec![2, 1]);
// [1,0] → 权重 2, [0,1] → 权重 1 → [1,0] 更大
assert_eq!(ord.cmp(&[1, 0], &[0, 1]), std::cmp::Ordering::Less);
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

### BlockOrder 与 SubOrder

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockOrder {
    boundaries: SmallVec<[usize; 4]>,
    orders: SmallVec<[SubOrder; 4]>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubOrder {
    Lex,
    Grevlex,
    Grlex,
}
```

**功能**：分块消元序。将变量分为连续的块，每个块内使用独立的子序。

**构造方法**：

```rust
pub fn new(boundaries: SmallVec<[usize; 4]>, orders: SmallVec<[SubOrder; 4]>) -> Self
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `boundaries` | `SmallVec<[usize; 4]>` | 排序的分割点列表（排他上界，不含 `n_vars`） |
| `orders` | `SmallVec<[SubOrder; 4]>` | 每个块的子序，长度必须等于 `boundaries.len() + 1` |

**说明**：`boundaries = [2]` 加 `orders = [Lex, Grevlex]` 在 4 变量多项式上表示：先比较变量 0–1（Lex），相等再比较变量 2–3（Grevlex）。

**Default**：单块 Grevlex。

**示例**：

```rust
use ocas_poly::sparse::{BlockOrder, MonomialOrder, SubOrder};
use smallvec::smallvec;

let ord = BlockOrder::new(
    smallvec![2],
    smallvec![SubOrder::Lex, SubOrder::Grevlex],
);
let a = [1, 0, 0, 0]; // x₀
let b = [0, 1, 0, 0]; // x₁
// 块内 Lex：[1,0] > [0,1]
assert_eq!(ord.cmp(&a, &b), std::cmp::Ordering::Greater);
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

### MatrixOrder

**签名**：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixOrder {
    matrix: Vec<Vec<i64>>,
    n_vars: usize,
}
```

**功能**：矩阵序。将指数向量乘以整数矩阵后按字典序比较。给定 $n \times n$ 矩阵 $M$，单项式 $\alpha > \beta$ 当且仅当 $M\alpha >_{\text{lex}} M\beta$。泛化所有标准序，特别适合构造消元序。

**构造方法**：

```rust
pub fn new(matrix: Vec<Vec<i64>>) -> Self
pub fn elimination_order(elim_vars: usize, n_vars: usize) -> Self
```

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `matrix` | `Vec<Vec<i64>>` | $n \times n$ 权重矩阵（行优先） |
| `elim_vars` | `usize` | 要消元的前几个变量数 |
| `n_vars` | `usize` | 总变量数 |

**`elimination_order`**：构造等价于 `BlockOrder([elim_vars in Lex, rest in Grevlex])` 的矩阵序。

**Default**：$1 \times 1$ 单位矩阵。

**示例**：

```rust
use ocas_poly::sparse::{MatrixOrder, MonomialOrder};

// 2×2 单位矩阵 = Lex 序
let ord = MatrixOrder::new(vec![vec![1, 0], vec![0, 1]]);
assert_eq!(ord.cmp(&[1, 0], &[0, 1]), std::cmp::Ordering::Greater);
```

**参见**：[MonomialOrder trait](#monomialorder-trait)

---

## DenseUnivariatePolynomial

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseUnivariatePolynomial<D: Domain> {
    coeffs: Vec<D::Element>,  // 从常数项开始，已去除尾零
    domain: D,
}
```

稠密一元多项式。系数存储在连续向量中，从常数项 $a_0$ 到最高次项 $a_n$。零多项式用空向量表示。乘法自动在 Karatsuba（两个多项式系数个数均 $\geq 32$ 时）和学校乘法之间选择。

---

### Dense 构造与属性

#### DenseUnivariatePolynomial::new

**签名**：`pub fn new(domain: D) -> Self`

**功能**：创建零多项式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `domain` | `D` | 系数域 |

**返回值**：零多项式（空系数向量）。

**参见**：[from_coeffs](#denseunivariatepolynomialfrom_coeffs)

---

#### DenseUnivariatePolynomial::from_coeffs

**签名**：`pub fn from_coeffs(domain: D, coeffs: Vec<D::Element>) -> Self`

**功能**：从系数向量构造多项式。自动去除尾零系数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `domain` | `D` | 系数域 |
| `coeffs` | `Vec<D::Element>` | 系数 $[a_0, a_1, \dots, a_n]$，常数项在前 |

**返回值**：多项式 $a_0 + a_1 x + \cdots + a_n x^n$。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let domain = IntegerDomain;
let p = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(0), Integer::from(2)],
);
assert_eq!(p.degree(), Some(2));
assert_eq!(p.coeff(2), Some(&Integer::from(2)));
// p(x) = 1 + 2x^2，中间的零系数被保留（非尾零）
```

**参见**：[new](#denseunivariatepolynomialnew)

---

#### DenseUnivariatePolynomial::domain

**签名**：`pub fn domain(&self) -> &D`

**功能**：返回系数域的引用。

---

#### DenseUnivariatePolynomial::coeffs

**签名**：`pub fn coeffs(&self) -> &[D::Element]`

**功能**：返回系数切片，从常数项开始。

---

#### DenseUnivariatePolynomial::is_zero

**签名**：`pub fn is_zero(&self) -> bool`

**功能**：判断是否为零多项式。

---

#### DenseUnivariatePolynomial::degree

**签名**：`pub fn degree(&self) -> Option<usize>`

**功能**：返回多项式的次数。零多项式返回 `None`。

---

#### DenseUnivariatePolynomial::coeff

**签名**：`pub fn coeff(&self, n: usize) -> Option<&D::Element>`

**功能**：返回 $x^n$ 的系数，不存在则返回 `None`。

---

#### DenseUnivariatePolynomial::leading_coeff

**签名**：`pub fn leading_coeff(&self) -> Option<&D::Element>`

**功能**：返回首项系数。零多项式返回 `None`。

---

#### DenseUnivariatePolynomial::lcoeff

**签名**：`pub fn lcoeff(&self) -> D::Element`

**功能**：返回首项系数的便捷别名。零多项式返回域的零元。

**返回值**：首项系数，或域的零元。

**参见**：[leading_coeff](#denseunivariatepolynomialleading_coeff)

---

#### DenseUnivariatePolynomial::constant

**签名**：`pub fn constant(&self) -> D::Element`

**功能**：返回常数项（$x^0$ 的系数）。零多项式返回域的零元。

---

#### DenseUnivariatePolynomial::zero

**签名**：`pub fn zero(&self) -> Self`

**功能**：返回同域的零多项式。

---

#### DenseUnivariatePolynomial::one

**签名**：`pub fn one(&self) -> Self`

**功能**：返回常数多项式 $1$。

---

#### DenseUnivariatePolynomial::is_one

**签名**：`pub fn is_one(&self) -> bool`

**功能**：判断是否为常数多项式 $1$。

---

### Dense 基本算术

以下方法在 `D: Domain` 上可用：

#### DenseUnivariatePolynomial::neg

**签名**：`pub fn neg(&self) -> Self`

**功能**：返回 $-p(x)$。

---

#### DenseUnivariatePolynomial::add

**签名**：`pub fn add(&self, other: &Self) -> Self`

**功能**：多项式加法。

---

#### DenseUnivariatePolynomial::sub

**签名**：`pub fn sub(&self, other: &Self) -> Self`

**功能**：多项式减法。

---

#### DenseUnivariatePolynomial::mul_scalar

**签名**：`pub fn mul_scalar(&self, scalar: &D::Element) -> Self`

**功能**：标量乘法。所有系数乘以 `scalar`。

---

#### DenseUnivariatePolynomial::mul

**签名**：`pub fn mul(&self, other: &Self) -> Self`

**功能**：多项式乘法。根据次数自动选择学校乘法或 Karatsuba。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let domain = IntegerDomain;
let a = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(1)],
);
let b = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(-1)],
);
let c = a.mul(&b);
assert_eq!(c.coeffs(), &[
    Integer::from(1),
    Integer::from(0),
    Integer::from(-1),
]);
// (1 + x)(1 - x) = 1 - x^2
```

**参见**：[mul_into](#denseunivariatepolynomialmul_into)

---

#### DenseUnivariatePolynomial::mul_into

**签名**：`pub fn mul_into(&self, other: &Self, buf: &mut Vec<D::Element>)`

**功能**：乘法，将结果写入缓冲区（避免热循环中的重复堆分配）。缓冲区被清空后重用。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Self` | 另一个多项式 |
| `buf` | `&mut Vec<D::Element>` | 输出缓冲区，调用后包含乘积系数 |

---

#### DenseUnivariatePolynomial::eval

**签名**：`pub fn eval(&self, x: &D::Element) -> D::Element`

**功能**：用 Horner 法在 $x$ 处求值。零多项式返回域的零元。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let domain = IntegerDomain;
let p = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(2), Integer::from(3)],
);
let value = p.eval(&Integer::from(2));
assert_eq!(value, Integer::from(17));
// p(2) = 1 + 2·2 + 3·4 = 17
```

---

#### DenseUnivariatePolynomial::derivative

**签名**：`pub fn derivative(&self) -> Self`

**功能**：返回形式导数 $p'(x)$。

---

#### DenseUnivariatePolynomial::integral

**签名**：`pub fn integral(&self) -> Self`

**功能**：返回形式积分 $\int p(x)\,dx$，常数项为零。要求域支持除法（如 `RationalDomain`）。

---

### EuclideanDomain 操作

以下方法要求 `D: EuclideanDomain`：

#### DenseUnivariatePolynomial::mul_coeff

**签名**：`pub fn mul_coeff(&self, c: &D::Element) -> Self`

**功能**：所有系数乘以常数 $c$。等价于 `mul_scalar`，但限于 `EuclideanDomain`。

**参见**：[mul_scalar](#denseunivariatepolynomialmul_scalar)

---

#### DenseUnivariatePolynomial::div_coeff

**签名**：`pub fn div_coeff(&self, c: &D::Element) -> Self`

**功能**：所有系数除以常数 $c$（必须整除）。

**错误**：若 $c$ 在域中不可逆（例如 $\mathbb{Z}$ 上 $c$ 不是 $\pm 1$），则无条件 panic（即使所有系数都能被 $c$ 整除）。

---

#### DenseUnivariatePolynomial::div_rem

**签名**：`pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>`

**功能**：带余除法。返回 $(q, r)$ 使得 $p = q \cdot \text{divisor} + r$，其中 $\deg(r) < \deg(\text{divisor})$。

**返回值**：`Some((quotient, remainder))`，除数为零时返回 `None`。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let domain = IntegerDomain;
let p = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(0), Integer::from(-1)],
);
let q = DenseUnivariatePolynomial::from_coeffs(
    domain,
    vec![Integer::from(1), Integer::from(1)],
);
let (quot, rem) = p.div_rem(&q).unwrap();
assert_eq!(quot.coeffs(), &[Integer::from(1), Integer::from(-1)]);
assert!(rem.is_zero());
// (x^2 - 1) / (x + 1) = x - 1，余数为 0
```

**参见**：[gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::gcd

**签名**：`pub fn gcd(&self, other: &Self) -> Self`

**功能**：计算两个多项式的最大公因式。使用 Euclid 算法（非域时用伪余式）。结果总是本原的（content-free）。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]); // x^2 - 1
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]); // (x+1)^2
let g = a.gcd(&b);
assert_eq!(g.coeffs(), &[Integer::from(1), Integer::from(1)]);
// gcd = x + 1
```

**参见**：[extended_gcd_poly](#denseunivariatepolynomialextended_gcd_poly)

---

#### DenseUnivariatePolynomial::content

**签名**：`pub fn content(&self) -> D::Element`

**功能**：返回所有系数的最大公因式。零多项式返回域的零元。

**返回值**：系数的 GCD。

**参见**：[primitive_part](#denseunivariatepolynomialprimitive_part)

---

#### DenseUnivariatePolynomial::primitive_part

**签名**：`pub fn primitive_part(&self) -> Self`

**功能**：返回本原部分 $p / \text{content}(p)$。结果的 content 为 1（或为零多项式）。

**参见**：[content](#denseunivariatepolynomialcontent)

---

#### DenseUnivariatePolynomial::extended_gcd_poly

**签名**：`pub fn extended_gcd_poly(&self, other: &Self) -> (Self, Self, Self)`

**功能**：扩展 Euclid 算法。返回 $(g, s, t)$ 使得 $s \cdot p + t \cdot q = g$，其中 $g = \gcd(p, q)$，$g$ 为首一多项式。

**返回值**：`(gcd, bezout_s, bezout_t)`

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;

let d = RationalDomain;
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(0, 1), Rational::new(1, 1),
]); // x^2 + 1
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(1, 1),
]); // x + 1
let (g, s, t) = a.extended_gcd_poly(&b);
// s·a + t·b = g
```

**参见**：[gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::pow

**签名**：`pub fn pow(&self, n: u32) -> Self`

**功能**：快速幂（反复平方法）。$p^0 = 1$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `u32` | 非负整数指数 |

**返回值**：$p(x)^n$。

---

#### DenseUnivariatePolynomial::p_adic_expansion

**签名**：`pub fn p_adic_expansion(&self, p: &Self) -> Vec<Self>`

**功能**：$p$-adic 展开。返回 $[a_0, a_1, a_2, \dots]$ 使得 $\text{self} = a_0 + a_1 \cdot p + a_2 \cdot p^2 + \cdots$，每个 $a_k$ 的次数小于 $\deg(p)$。

**实现**：通过反复带余除法。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// 对 f(x) = x^3 关于 p(x) = x + 1 做 p-adic 展开
let f = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(0), Integer::from(0), Integer::from(0), Integer::from(1),
]);
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(1),
]);
let expansion = f.p_adic_expansion(&p);
// expansion = [a0, a1, a2, ...] 使得 f = a0 + a1*p + a2*p^2 + ...
```

**参见**：[div_rem](#denseunivariatepolynomialdiv_rem)

---

#### DenseUnivariatePolynomial::diophantine

**签名**：`pub fn diophantine(polys: &mut [Self], b: &Self) -> Vec<Self>`

**功能**：多项式 CRT（丢番图求解器）。给定两两互素的多项式列表 `polys` 和目标 $b$，返回 $[s_0, \dots, s_n]$ 使得：

$$\sum_i s_i \cdot \prod_{j \neq i} p_j \equiv b \pmod{\prod_i p_i}$$

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `polys` | `&mut [Self]` | 两两互素的多项式列表 |
| `b` | `&Self` | 目标多项式 |

**返回值**：余式列表 $[s_0, \dots, s_n]$。

**错误**：若多项式非两两互素则 panic。

---

### 因式分解与结式

#### DenseUnivariatePolynomial::square_free_factorization

**签名**：`pub fn square_free_factorization(&self) -> SquareFreeFactors<D>`

**功能**：无平方分解。使用 Yun 算法：$g = \gcd(f, f')$，$w = f/g$，迭代 $h = \gcd(w, g)$，$z = w/h$ 收集重数 $k$ 因子。

**返回值**：`Vec<(DenseUnivariatePolynomial<D>, usize)>` — `(因子, 重数)` 对的列表。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// (x+1)^2*(x-1) = x^3 + x^2 - x - 1
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(-1), Integer::from(1), Integer::from(1),
]);
let factors = p.square_free_factorization();
assert_eq!(factors.len(), 2);
```

**参见**：[factor](#denseunivariatepolynomialfactor)、[gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::is_square_free

**签名**：`pub fn is_square_free(&self) -> bool`

**功能**：判断多项式是否无平方（$\gcd(p, p') = 1$）。

---

#### DenseUnivariatePolynomial::factor（IntegerDomain）

**签名**：`impl DenseUnivariatePolynomial<IntegerDomain> { pub fn factor(&self) -> Factors<IntegerDomain> }`

**功能**：将本原整数多项式完全分解为不可约因子。使用无平方分解 + Berlekamp–Zassenhaus + Hensel 提升算法。

**前置条件**：输入必须是本原的（content = 1）。对任意多项式先调用 `primitive_part()`。

**返回值**：`Vec<(DenseUnivariatePolynomial<IntegerDomain>, usize)>` — `(不可约因子, 重数)` 对。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x^2 - 1 = (x-1)(x+1)
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
let factors = p.factor();
assert_eq!(factors.len(), 2);
```

**参见**：[square_free_factorization](#denseunivariatepolynomialsquare_free_factorization)

---

#### DenseUnivariatePolynomial::factor（FiniteField）

**签名**：`impl DenseUnivariatePolynomial<FiniteField> { pub fn factor(&self) -> Factors<FiniteField> }`

**功能**：将 $\mathbb{F}_p$ 上的多项式完全分解。使用 Berlekamp 算法（或 Cantor–Zassenhaus）。

**返回值**：`Vec<(DenseUnivariatePolynomial<FiniteField>, usize)>` — `(不可约因子, 重数)` 对。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::DenseUnivariatePolynomial;

let f = FiniteField::new(BigInt::from(5));
// x^2 - 1 over F_5
let p = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(4), f.element(0), f.element(1)]);
let factors = p.factor();
assert!(!factors.is_empty());
```

**参见**：[factor（IntegerDomain）](#denseunivariatepolynomialfactorintegerdomain)

---

#### DenseUnivariatePolynomial::resultant

**签名**：`pub fn resultant(&self, other: &Self) -> D::Element`

**功能**：用 Brown PRS 算法计算结式 $\operatorname{Res}(a, b)$。结式为零当且仅当 $a$ 和 $b$ 有非常数公因子。

**实现细节**：子结式 PRS，每步用 $\beta$ 精确除法（子结式定理保证整除）。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// Res(x - 1, x - 2) = 1 - 2 = -1
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(1),
]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-2), Integer::from(1),
]);
assert_eq!(a.resultant(&b), Integer::from(-1));
```

**性质**：$\operatorname{Res}(a, b) = (-1)^{\deg a \cdot \deg b} \operatorname{Res}(b, a)$。

**参见**：[gcd](#denseunivariatepolynomialgcd)

---

## SparseMultivariatePolynomial

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultivariatePolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    terms: HashMap<SmallVec<[usize; 4]>, D::Element>,
    domain: D,
    n_vars: usize,
    pub order: O,
}
```

稀疏多元多项式。只存储非零项，以指数向量 $\vec{e} = [e_1, e_2, \dots]$ 为键、系数为值的 HashMap。指数向量 $\vec{e}$ 表示单项式 $x_1^{e_1} x_2^{e_2} \cdots$。单项式序由类型参数 `O` 控制，默认为 `Grevlex`。

---

### Sparse 构造与属性

#### SparseMultivariatePolynomial::new

**签名**：`pub fn new(domain: D, n_vars: usize) -> Self`

**功能**：创建零多项式，使用默认单项式序。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `domain` | `D` | 系数域 |
| `n_vars` | `usize` | 变量数 |

**参见**：[new_with_order](#sparsemultivariatepolynomialnew_with_order)

---

#### SparseMultivariatePolynomial::new_with_order

**签名**：`pub fn new_with_order(domain: D, n_vars: usize, order: O) -> Self`

**功能**：创建零多项式，指定单项式序。

**示例**：

```rust
use ocas_domain::IntegerDomain;
use ocas_poly::sparse::{SparseMultivariatePolynomial, WeightOrder};

let order = WeightOrder::from_slice(&[2, 1]);
let p = SparseMultivariatePolynomial::<_, WeightOrder>::new_with_order(
    IntegerDomain, 2, order,
);
assert_eq!(p.n_vars(), 2);
```

---

#### SparseMultivariatePolynomial::from_terms

**签名**：`pub fn from_terms(domain: D, n_vars: usize, terms: Vec<(Vec<usize>, D::Element)>) -> Self`

**功能**：从（指数向量, 系数）对的列表构造多项式。零系数项自动丢弃。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `domain` | `D` | 系数域 |
| `n_vars` | `usize` | 变量数 |
| `terms` | `Vec<(Vec<usize>, D::Element)>` | (指数向量, 系数) 对 |

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::sparse::Grevlex;
use ocas_poly::SparseMultivariatePolynomial;

let domain = IntegerDomain;
let p = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
    domain,
    2,
    vec![(vec![1, 0], Integer::from(2)), (vec![0, 1], Integer::from(3))],
);
assert_eq!(p.n_terms(), 2);
assert_eq!(p.coeff(&[1, 0]), Integer::from(2));
// p = 2x + 3y
```

**参见**：[new](#sparsemultivariatepolynomialnew)

---

#### SparseMultivariatePolynomial::domain

**签名**：`pub fn domain(&self) -> &D`

**功能**：返回系数域的引用。

---

#### SparseMultivariatePolynomial::n_vars

**签名**：`pub fn n_vars(&self) -> usize`

**功能**：返回变量数。

---

#### SparseMultivariatePolynomial::n_terms

**签名**：`pub fn n_terms(&self) -> usize`

**功能**：返回非零项数。

---

#### SparseMultivariatePolynomial::is_zero

**签名**：`pub fn is_zero(&self) -> bool`

**功能**：判断是否为零多项式。

---

#### SparseMultivariatePolynomial::terms_ref

**签名**：`pub fn terms_ref(&self) -> &HashMap<SmallVec<[usize; 4]>, D::Element>`

**功能**：返回内部项映射的引用（指数 → 系数）。

---

#### SparseMultivariatePolynomial::set_term_external

**签名**：`pub fn set_term_external(&mut self, exp: Vec<usize>, coeff: D::Element)`

**功能**：设置单项式的系数。零系数会移除该项。

---

#### SparseMultivariatePolynomial::total_degree

**签名**：`pub fn total_degree(&self) -> Option<usize>`

**功能**：返回总次数（所有单项式次数的最大值）。零多项式返回 `None`。

---

#### SparseMultivariatePolynomial::coeff

**签名**：`pub fn coeff(&self, exp: &[usize]) -> D::Element`

**功能**：返回给定单项式的系数，不存在则返回域的零元。

---

#### SparseMultivariatePolynomial::degree_in

**签名**：`pub fn degree_in(&self, var_index: usize) -> usize`

**功能**：返回关于变量 `var_index` 的次数。零多项式返回 0。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var_index` | `usize` | 变量索引（0-based） |

**返回值**：该变量在所有项中的最大指数。

---

#### SparseMultivariatePolynomial::zero / one

**签名**：

```rust
pub fn zero(&self) -> Self
pub fn one(&self) -> Self
```

**功能**：返回同形状的零多项式或常数 $1$。

---

### Sparse 基本算术

#### SparseMultivariatePolynomial::neg

**签名**：`pub fn neg(&self) -> Self`

**功能**：返回 $-p$。

---

#### SparseMultivariatePolynomial::add

**签名**：`pub fn add(&self, other: &Self) -> Self`

**功能**：多项式加法。

**错误**：若变量数不同则 panic。

---

#### SparseMultivariatePolynomial::sub

**签名**：`pub fn sub(&self, other: &Self) -> Self`

**功能**：多项式减法。

**错误**：若变量数不同则 panic。

---

#### SparseMultivariatePolynomial::mul_scalar

**签名**：`pub fn mul_scalar(&self, scalar: &D::Element) -> Self`

**功能**：标量乘法。

---

#### SparseMultivariatePolynomial::mul

**签名**：`pub fn mul(&self, other: &Self) -> Self`

**功能**：多项式乘法。对每对项分别相乘后合并同类项。

**错误**：若变量数不同则 panic。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::sparse::Grevlex;
use ocas_poly::SparseMultivariatePolynomial;

let domain = IntegerDomain;
let p = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
    domain,
    2,
    vec![(vec![1, 0], Integer::from(2)), (vec![0, 1], Integer::from(3))],
);
let q = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
    domain,
    2,
    vec![(vec![1, 0], Integer::from(1)), (vec![0, 0], Integer::from(1))],
);
let r = p.mul(&q);
assert_eq!(r.coeff(&[1, 0]), Integer::from(2));
assert_eq!(r.coeff(&[0, 1]), Integer::from(3));
assert_eq!(r.coeff(&[2, 0]), Integer::from(2));
// (2x + 3y)(x + 1) = 2x^2 + 2x + 3xy + 3y
```

---

#### SparseMultivariatePolynomial::mul_monomial

**签名**：`pub fn mul_monomial(&self, exp: &[usize]) -> Self`

**功能**：每个项的指数向量逐分量加上 `exp`。用于 Gröbner 基归约。

---

#### SparseMultivariatePolynomial::sorted_terms

**签名**：`pub fn sorted_terms(&self) -> Vec<(&SmallVec<[usize; 4]>, &D::Element)>`

**功能**：按单项式序返回排序后的（指数, 系数）对。

---

#### SparseMultivariatePolynomial::eval

**签名**：`pub fn eval(&self, var_index: usize, value: &D::Element) -> Self`

**功能**：将变量 `var_index` 替换为值 `value`，返回少一个变量的多项式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `var_index` | `usize` | 要替换的变量索引 |
| `value` | `&D::Element` | 替换值 |

**返回值**：少一个变量的多项式（剩余变量保持相对顺序）。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 2,
    vec![
        (vec![1, 1], Integer::from(1)), // xy
        (vec![0, 1], Integer::from(2)), // 2y
    ],
);
// 代入 x=3：结果 = 3y + 2y = 5y
let r = p.eval(0, &Integer::from(3));
assert_eq!(r.coeff(&[1]), Integer::from(5));
```

---

#### SparseMultivariatePolynomial::eval_keep

**签名**：`pub fn eval_keep(&self, var_index: usize, value: &D::Element) -> Self`

**功能**：将变量 `var_index` 替换为值，但保持变量总数不变（被替换变量的指数置零）。用于 Hensel 提升中变量位置需保持固定的场景。

---

#### SparseMultivariatePolynomial::leading_term

**签名**：`pub fn leading_term(&self) -> Option<(&SmallVec<[usize; 4]>, &D::Element)>`

**功能**：返回首项（指数向量, 系数）。零多项式返回 `None`。$O(n)$ 扫描 HashMap。

---

#### SparseMultivariatePolynomial::leading_monomial

**签名**：`pub fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>>`

**功能**：返回首单项式的指数向量。

---

#### SparseMultivariatePolynomial::leading_coeff

**签名**：`pub fn leading_coeff(&self) -> Option<&D::Element>`

**功能**：返回首项系数。

---

#### SparseMultivariatePolynomial::content

**签名**：

```rust
pub fn content(&self) -> D::Element
where D: EuclideanDomain
```

**功能**：返回所有系数的 GCD。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 1,
    vec![(vec![2], Integer::from(6)), (vec![1], Integer::from(9)), (vec![0], Integer::from(3))],
);
assert_eq!(p.content(), Integer::from(3));
```

**参见**：[primitive_part](#sparsemultivariatepolynomialprimitive_part)

---

#### SparseMultivariatePolynomial::primitive_part

**签名**：

```rust
pub fn primitive_part(&self) -> Self
where D: EuclideanDomain
```

**功能**：返回本原部分（多项式 / content）。

**参见**：[content](#sparsemultivariatepolynomialcontent)

---

#### SparseMultivariatePolynomial::div_exact

**签名**：`pub fn div_exact(&self, divisor: &Self) -> Self`

**功能**：精确除法（假设无余式）。用于有理函数规范化中已知 GCD 整除的场景。

**错误**：debug 模式下若除法不精确则 panic。

**参见**：[checked_div_exact](#sparsemultivariatepolynomialchecked_div_exact)

---

#### SparseMultivariatePolynomial::checked_div_exact

**签名**：`pub fn checked_div_exact(&self, divisor: &Self) -> Option<Self>`

**功能**：精确除法，不精确时返回 `None`。

**返回值**：`Some(quotient)` 或 `None`。

---

#### SparseMultivariatePolynomial::derivative

**签名**：`pub fn derivative(&self, var_index: usize) -> Self`

**功能**：对变量 `var_index` 求形式偏导。

---

#### SparseMultivariatePolynomial::taylor_coefficients

**签名**：`pub fn taylor_coefficients(&self, var_index: usize, a: &D::Element) -> Vec<Self>`

**功能**：计算关于变量 `var_index` 在点 $a$ 处的 Taylor 系数。返回 $[t_0, t_1, \dots, t_d]$ 使得 $f = \sum_j t_j (x_{\text{var}} - a)^j$。

---

### Gröbner 基支持

#### SparseMultivariatePolynomial::reduce

**签名**：`pub fn reduce(&self, basis: &[Self]) -> Self`

**功能**：多元多项式除法。反复寻找能除尽当前首项的基元素，减去适当倍数；否则将首项移入余式。要求域为域（div 总是成功）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `basis` | `&[Self]` | 除数列表 |

**返回值**：归约后的余式。

---

#### SparseMultivariatePolynomial::spoly

**签名**：`pub fn spoly(&self, other: &Self) -> Self`

**功能**：计算 S-多项式。$S(f, g) = \frac{\text{lcm}}{\text{lt}(f)} \cdot f - \frac{\text{lcm}}{\text{lt}(g)} \cdot g$。

---

#### SparseMultivariatePolynomial::make_monic_inplace

**签名**：`pub fn make_monic_inplace(&mut self) -> bool`

**功能**：原地将多项式化为首一（除以首项系数）。首项系数不可逆时返回 `false`。

---

#### SparseMultivariatePolynomial::exponents_iter

**签名**：`pub fn exponents_iter(&self) -> impl Iterator<Item = &SmallVec<[usize; 4]>>`

**功能**：按单项式序（降序）迭代所有指数向量。用于 F4 算法的符号预处理。

---

### 多元因式分解

#### SparseMultivariatePolynomial::factor（IntegerDomain, Lex）

**签名**：

```rust
impl SparseMultivariatePolynomial<IntegerDomain, Lex> {
    pub fn factor(&self) -> Vec<(Self, usize)>
}
```

**功能**：将多元整数多项式分解为不可约因子。二元且主变量首项系数为常数时使用双变量 Hensel 提升（Wang）；三元及以上、或主变量首项系数非常数的二元情形使用 EEZ 算法。

**返回值**：`(不可约因子, 重数)` 对的列表。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

// (x^2 + y + 1)(x + y + 2)
let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 2,
    vec![
        (vec![3, 0], Integer::from(1)),
        (vec![2, 1], Integer::from(1)),
        (vec![2, 0], Integer::from(2)),
        (vec![1, 1], Integer::from(1)),
        (vec![1, 0], Integer::from(1)),
        (vec![0, 2], Integer::from(1)),
        (vec![0, 1], Integer::from(3)),
        (vec![0, 0], Integer::from(2)),
    ],
);
let factors = f.factor();
assert!(factors.len() >= 2);
```

**参见**：[factor（FiniteField）](#sparsemultivariatepolynomialfactorfinitefield-lex)

---

#### SparseMultivariatePolynomial::factor（FiniteField, Lex）

**签名**：

```rust
impl SparseMultivariatePolynomial<FiniteField, Lex> {
    pub fn factor(&self) -> Vec<(Self, usize)>
}
```

**功能**：将多元 $\mathbb{F}_p$ 多项式分解。二元用求值–Hensel 路径，三元及以上用 EEZ。

**参见**：[factor（IntegerDomain）](#sparsemultivariatepolynomialfactorintegerdomain-lex)

---

## RationalPolynomial

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalPolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    pub numerator: SparseMultivariatePolynomial<D, O>,
    pub denominator: SparseMultivariatePolynomial<D, O>,
}
```

多项式分式域元素 $\frac{\text{num}}{\text{den}}$。通过 `from_num_den` 构造时自动规范化为最简形式（分子分母互素，分母首项系数归一化）。

---

#### RationalPolynomial::new

**签名**：

```rust
pub fn new(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

**功能**：不化简地构造有理多项式。调用者须确保分母非零。

**参见**：[from_num_den](#rationalpolynomialfrom_num_den)

---

#### RationalPolynomial::from_num_den

**签名**：

```rust
impl<D: EuclideanDomain, O: MonomialOrder> RationalPolynomial<D, O> {
    pub fn from_num_den(
        numerator: SparseMultivariatePolynomial<D, O>,
        denominator: SparseMultivariatePolynomial<D, O>,
    ) -> Self
}
```

**功能**：从分子分母构造并规范化。结果为最简形式：分子分母互素，分母首项系数归一化（有限域下为 1，整数域下为正）。

**前置条件**：分母非零。

**规范化步骤**：
1. 计算 $\gcd(\text{num}, \text{den})$ 并约去
2. 一元时用稠密 Euclid 算法精确 GCD 约简
3. 归一化分母首项系数

**错误**：分母为零时 panic。

**参见**：[new](#rationalpolynomialnew)

---

#### RationalPolynomial::from_polynomial

**签名**：`pub fn from_polynomial(poly: SparseMultivariatePolynomial<D, O>) -> Self`

**功能**：从多项式构造（分母 = 1）。

---

#### RationalPolynomial::zero / one

**签名**：

```rust
pub fn zero(domain: &D, n_vars: usize) -> Self
pub fn one(domain: &D, n_vars: usize) -> Self
```

**功能**：返回零或单位有理多项式。

---

#### RationalPolynomial::is_zero / is_one

**签名**：

```rust
pub fn is_zero(&self) -> bool
pub fn is_one(&self) -> bool
```

**功能**：判断是否为零或 $1/1$。

---

#### RationalPolynomial::n_vars

**签名**：`pub fn n_vars(&self) -> usize`

**功能**：返回变量数。

---

#### RationalPolynomial::domain

**签名**：`pub fn domain(&self) -> &D`

**功能**：返回系数域引用。

---

#### RationalPolynomial::neg

**签名**：`pub fn neg(&self) -> Self`

**功能**：返回 $-\frac{n}{d}$。

---

#### RationalPolynomial::inv

**签名**：`pub fn inv(&self) -> Option<Self>`

**功能**：返回乘法逆 $\frac{d}{n}$。分子为零时返回 `None`。

**返回值**：`Some(逆)` 或 `None`（分子为零）。

---

#### RationalPolynomial::pow

**签名**：`pub fn pow(&self, k: u32) -> Self`

**功能**：快速幂。$\left(\frac{n}{d}\right)^k$，分子分母分别用反复平方法。

---

#### RationalPolynomial::add

**签名**：

```rust
impl<D: EuclideanDomain, O: MonomialOrder> RationalPolynomial<D, O> {
    pub fn add(&self, other: &Self) -> Self
}
```

**功能**：有理多项式加法。同分母时直接加分子；异分母时交叉相乘后规范化。

---

#### RationalPolynomial::sub

**签名**：`pub fn sub(&self, other: &Self) -> Self`

**功能**：有理多项式减法。等价于 `self.add(&other.neg())`。

---

#### RationalPolynomial::mul

**签名**：`pub fn mul(&self, other: &Self) -> Self`

**功能**：有理多项式乘法。交叉相乘后规范化。

---

#### RationalPolynomial::div

**签名**：`pub fn div(&self, other: &Self) -> Option<Self>`

**功能**：有理多项式除法。$\frac{a/b}{c/d} = \frac{ad}{bc}$。

**返回值**：`Some(商)` 或 `None`（除数分子为零）。

---

## 辅助函数

以下函数位于 `ocas_poly::sparse` 模块中：

### monomial_divides

**签名**：`pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool`

**功能**：判断单项式 $b$ 是否整除 $a$（即 $a$ 是 $b$ 的倍式）。返回 `true` 当且仅当 $a_i \geq b_i$ 对所有 $i$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `&[usize]` | 被除数的指数向量（倍式） |
| `b` | `&[usize]` | 除数的指数向量 |

**返回值**：`true` 若 $b$ 整除 $a$。

---

### monomial_lcm

**签名**：`pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]>`

**功能**：计算两个单项式的最小公倍式：逐分量取最大值。

$$\text{lcm}(x^a y^b, x^c y^d) = x^{\max(a,c)} y^{\max(b,d)}$$

---

### monomial_are_coprime

**签名**：`pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool`

**功能**：判断两个单项式是否互素（没有变量同时出现在两者中）。

---

## 类型别名

```rust
/// 无平方分解结果：(因子, 重数) 对的列表。
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

/// 完全分解结果：(不可约因子, 重数) 对的列表。
pub type Factors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;
```
