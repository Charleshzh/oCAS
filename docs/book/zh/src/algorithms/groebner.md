# Gröbner 基

oCAS 可在任意域上计算多元多项式理想的 Gröbner 基。提供三种算法以及
换序工具，本章对比它们并说明何时使用哪一种。

---

## 范围

| 算法 | 入口 | 适用场景 |
|---|---|---|
| Buchberger | `GroebnerBasis::buchberger` | 小理想、教学 |
| **F4** | `f4::f4` | 生产环境 —— 默认 |
| **F5** | `f5::f5` | 基于签名；正则序列 |

换序：

| 工具 | 入口 | 适用场景 |
|---|---|---|
| 重跑 F4 | `GroebnerBasis::reorder` | 一般理想 |
| **FGLM** | `fglm::fglm` | 零维理想（快得多） |

---

## Buchberger 与 F4 对比

Buchberger 算法逐个处理 S-多项式。F4（Faugère 1999）把大量 S-多项式
约化批处理成一次稀疏矩阵行阶梯计算，对中等和大理想显著更快——因为
线性代数占主导，并且可以优化（缓存友好的稀疏行、ℤ_p 原生 `i64`
快速路径）。

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{SparseMultivariatePolynomial, f4};

let d = RationalDomain;
// cyclic-3 方程组
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
    (vec![1, 0, 0], Rational::new(1, 1)),
    (vec![0, 1, 0], Rational::new(1, 1)),
    (vec![0, 0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
    (vec![1, 1, 0], Rational::new(1, 1)),
    (vec![0, 1, 1], Rational::new(1, 1)),
    (vec![1, 0, 1], Rational::new(1, 1)),
]);
let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
    (vec![1, 1, 1], Rational::new(1, 1)),
    (vec![0, 0, 0], Rational::new(-1, 1)),
]);
let gb = f4::f4(&[f1, f2, f3]);
assert!(gb.is_groebner_basis());
```

F4 使用 Gebauer–Moeller 临界对筛选（第一、第二判据加冗余对清理）
和逐基元素简化缓存，因此它构造的矩阵接近最小。

---

## 单项式序与 `reorder`

支持 `Lex`、`Grlex` 与 `Grevlex` 序。消元理论需要 `Lex` 基，但它
通常计算代价最高。标准策略是：

1. 先计算 `Grevlex` 基（最快），
2. 再转换为 `Lex`。

对一般理想，`reorder` 在新序下重新解释基并重跑 F4：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::{Grevlex, Lex};
use ocas_poly::{SparseMultivariatePolynomial, f4};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb_grevlex = f4::f4(&[f1, f2]);
let gb_lex = gb_grevlex.reorder::<Lex>();
assert!(gb_lex.is_groebner_basis());
```

---

## FGLM：零维理想的快速换序

零维理想（有限多个公共根）具有有限的*阶梯*——不被任何首项单项式
整除的单项式集合。FGLM 算法（Faugère–Gianni–Lazard–Mora 1993）按
目标序递增遍历单项式，计算它们对现有基的正规形，并检测线性相关。
每次相关产生新基的一个多项式。代价为 `O(n·D³)` 次域运算，其中 `D`
是阶梯维数，与产生原基的 F4 代价无关。

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::{Grevlex, Lex};
use ocas_poly::{SparseMultivariatePolynomial, f4};
use ocas_poly::groebner::fglm::fglm;

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb_lex = f4::f4(&[f1, f2]);
let gb_grevlex = fglm::<_, Grevlex>(&gb_lex).expect("zero-dimensional");
assert!(gb_grevlex.is_groebner_basis());
```

当理想为正维（无限阶梯）时 `fglm` 返回 `None`，此时请改用
`reorder`。

---

## F5：基于签名的 Gröbner 基

`f5::f5` 实现 Faugère 的签名判据（2002）：签名已出现的 S-对会被
跳过，可证明对正则序列避免全部零约化。自 0.19.0 版本起，F5 已是
生产级实现，在大理想上性能与 F4 竞争。

`hilbert` 模块通过容斥原理计算单项式理想的 Hilbert 分子，给出阶梯
的正则性——一个可靠的次数界，F4 可将其用作提前终止提示。

---

## 基准

Criterion 计时（cyclic 方程组，ℚ 与 ℤ₁₃ 上，本机）：

| 方程组 | Buchberger | F4 | 加速比 |
|---|---|---|---|
| cyclic-3 ℚ | 308 µs | 147 µs | 2.1× |
| cyclic-4 ℚ | 3.99 ms | 2.13 ms | 1.9× |
| cyclic-3 ℤ₁₃ | 582 µs | 276 µs | 2.1× |
| cyclic-4 ℤ₁₃ | 6.19 ms | 2.80 ms | 2.2× |

ℤ_p 原生 `i64` 快速路径（行阶梯步骤中的惰性模算术）使有限域计时
接近有理数计时，尽管系数更小。

---

## 单项式序

0.19.1 版本扩展了单项式序系统，在内置 `Lex`、`Grlex` 和 `Grevlex`
之外支持运行时可配置的序：

| 序 | 入口 | 用途 |
|---|---|---|
| `WeightOrder` | `WeightOrder::from_slice(&[2, 1])` | 通过变量权重实现消元序 |
| `BlockOrder` | `BlockOrder::new(boundaries, orders)` | 分块消元：将变量分组，各组独立排序 |

`WeightOrder` 按加权和 $\sum w_i e_i$ 降序比较单项式。`BlockOrder`
将变量划分为连续块，各块按独立子序（Lex 或 Grevlex）比较。

```rust
use ocas_poly::sparse::WeightOrder;

let order = WeightOrder::from_slice(&[2, 1]);
let p = SparseMultivariatePolynomial::new_with_order(d, 2, order);
```

### MatrixOrder 与消元

0.23.0 版本新增 `MatrixOrder`，基于权重矩阵的通用单项式序，支持消元：

| 序 | 入口 | 用途 |
|---|---|---|
| `MatrixOrder` | `MatrixOrder::new(matrix)` | 通用权重矩阵序 |
| 消元序 | `MatrixOrder::elimination_order(elim_vars, n_vars)` | 消去前 `elim_vars` 个变量 |

```rust
use ocas_poly::sparse::{MatrixOrder, MonomialOrder};

// 消元序：先消去 x_0，剩余变量按 Grevlex 比较。
let ord = MatrixOrder::elimination_order(1, 3);
```

---

## 理想运算

0.23.0 版本在 `ocas_poly::ideal` 中引入完整理想运算库。所有运算使用
`Lex` 序以保持一致性。

| 运算 | 入口 | 描述 |
|---|---|---|
| 成员判定 | `ideal_contains(gens, f, algo)` | 测试 $f \in I$ |
| 和 | `ideal_sum(I, J)` | $I + J$ |
| 积 | `ideal_product(I, J)` | $I \cdot J$ |
| 商 | `ideal_quotient(I, J)` | $I : J$（Rabinowitsch 技巧）|
| 饱和 | `ideal_saturate(I, J)` | $I : J^\infty$ |
| 交集 | `ideal_intersection(I, J)` | $I \cap J$（辅助变量法）|
| 消元 | `eliminate(gens, elim_vars, algo)` | 通过 Lex GB 消元 |

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::{ideal_contains, ideal_saturate};
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 2], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]); // x² + y² - 1
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]); // x - y

// 成员判定：x² + y² - 1 ∈ ⟨x - y, x² + y² - 1⟩
assert!(ideal_contains(&[f2.clone()], &f1, Algorithm::Auto));
```

---

## 零维求解

`solve_polynomial_system` 对方程组分类并求实数解：

| 解类型 | 描述 |
|---|---|
| `ZeroDimensional` | 有限个实数解，通过 Sturm 根隔离 |
| `PositiveDimensional` | 无限解集；返回 Lex GB |
| `Empty` | 无解（理想为 $\langle 1 \rangle$）|

```rust
use ocas_poly::ideal::solve_polynomial_system;

// 圆 ∩ 直线：x² + y² = 1, x = y
let sol = solve_polynomial_system(&[f1, f2], Algorithm::Auto);
// 返回两个实数解 (±1/√2, ±1/√2)
```

使用 `is_zero_dimensional(&gb)` 检查维数而不求解。

---

## 准素分解与根式

对零维理想，准素分解通过 Lex GB 中一元多项式的因式分解和饱和分离分量：

```rust
use ocas_poly::ideal::{primary_decomposition, ideal_radical};

// (x², xy) = (x) ∩ (x², y)
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let decomp = primary_decomposition(&[f1, f2]);
assert_eq!(decomp.len(), 2);
```

`ideal_radical` 计算 $\sqrt{I}$：
- **零维**：一元多项式无平方分解
- **正维**：Jacobian 饱和法（$\sqrt{I} = I : h^\infty$）

```rust
let rad = ideal_radical(&[f1, f2]);
// √(x², xy) = (x)
```

---

## Hilbert 级数

`hilbert` 模块提供从 Gröbner 基计算完整 Hilbert 级数的功能：

| 方法 | 描述 |
|---|---|
| `hilbert_function(d)` | $\dim_k (R/I)_d$ 在次数 $d$ 处的值 |
| `dimension()` | $R/I$ 的 Krull 维数 |
| `degree()` | 射影簇的次数 |
| `hilbert_polynomial()` | 完整多项式系数（Lagrange 插值）|

```rust
use ocas_poly::groebner::hilbert::hilbert_series;

let hs = hilbert_series(&gb);
println!("H(5) = {}", hs.hilbert_function(5));
println!("dim = {}", hs.dimension());
println!("polynomial = {:?}", hs.hilbert_polynomial());
```
