# Gröbner 基

oCAS 可在任意域上计算多元多项式理想的 Gröbner 基。提供三种算法以及
换序工具，本章对比它们并说明何时使用哪一种。

---

## 范围

| 算法 | 入口 | 适用场景 |
|---|---|---|
| Buchberger | `GroebnerBasis::buchberger` | 小理想、教学 |
| **F4** | `f4::f4` | 生产环境 —— 默认 |
| F5（实验性） | `f5::f5` | 研究、签名剪枝 |

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

## F5 与 Hilbert 界（实验性）

`f5::f5` 实现 Faugère 的签名判据（2002）：签名已出现的 S-对会被
跳过，可证明对正则序列避免全部零约化。该实现保留用于研究；在测试
套件的理想上签名规则很少触发，因此 F4 仍是推荐默认。

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
