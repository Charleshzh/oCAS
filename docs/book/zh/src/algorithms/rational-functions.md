# 有理函数、结式与部分分式

oCAS 提供完整的有理函数运算、结式计算和部分分式分解栈，支持任何欧几里得域。
这些功能在 0.12.0 版本中添加，弥补了与 Symbolica 的 `rational_polynomial.rs`、
`resultant.rs` 和 `partial_fraction.rs` 的差距。

---

## 有理多项式

`RationalPolynomial<D, O>` 表示多项式环分式域中的元素——即 $\frac{p}{q}$，
其中 $p$ 和 $q$ 是域 `D` 上的多元多项式。

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::{RationalPolynomial, SparseMultivariatePolynomial, Lex};

// 创建多项式：x + 1 和 x - 1
let x_plus_1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 1,
    vec![(vec![0], Integer::from(1)), (vec![1], Integer::from(1))],
);
let x_minus_1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 1,
    vec![(vec![0], Integer::from(-1)), (vec![1], Integer::from(1))],
);

// (x+1) / (x-1)
let rat = RationalPolynomial::from_num_den(x_plus_1, x_minus_1);
```

### 规范化

通过 `from_num_den` 构造时，分数自动约简为规范形式：

1. 分子分母的 GCD 被约去。
2. 分母的首项系数被归一化（有序域为正，有限域为 1）。

### 算术运算

支持所有标准运算：

| 运算 | 方法 | 策略 |
|---|---|---|
| 加法 | `a.add(&b)` | 交叉相乘，然后规范化 |
| 减法 | `a.sub(&b)` | 通过取反 + 加法 |
| 乘法 | `a.mul(&b)` | 交叉消去 GCD，然后相乘 |
| 除法 | `a.div(&b)` | 通过逆元 + 乘法 |
| 取反 | `a.neg()` | 取反分子 |
| 逆元 | `a.inv()` | 交换分子分母 |
| 幂 | `a.pow(n)` | 快速幂 |

---

## 结式

两个多项式 $a$ 和 $b$ 的结式是一个标量，当且仅当 $a$ 和 $b$ 有公共根
（或等价地，有非平凡 GCD）时为零。

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// Res(x^2 + 1, (x+1)^2) = 4
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(0), Integer::from(1),
]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
assert_eq!(a.resultant(&b), Integer::from(4));
```

### 算法：Brown PRS

oCAS 使用 **Brown 多项式余式序列**算法，避免构造完整的 Sylvester 矩阵。
该算法跟踪首项系数和次数差，用基本定理公式从 PRS 计算结式。

对于两个 15 次多项式，结式在 20ms 内完成。

---

## 部分分式

给定真分式 $\frac{p(x)}{q(x)}$（其中 $\deg(p) < \deg(q)$），
部分分式分解将其表示为更简单分式的和。

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(0, 1), Rational::new(-1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
// poly_part 为 None（真分式）
// terms 包含分解后的分式
```

### 算法

分解步骤：

1. **多项式除法**：若 $\deg(p) \geq \deg(q)$，通过 `div_rem` 提取整式部分。
2. **无平方分解**：将分母分解为无平方因子 $f_i^{e_i}$。
3. **Diophantine CRT**：对多个因子，求解多项式中国剩余定理以拆分分式。
4. **p-adic 展开**：对重因子（$e_i > 1$），对分子做 p-adic 展开得到各项。

`together` 函数执行逆操作，将各项合并回单个有理函数。

---

## 辅助方法

为支持这些算法，`DenseUnivariatePolynomial` 新增了以下方法：

| 方法 | 说明 |
|---|---|
| `extended_gcd_poly(&self, other)` | 扩展 GCD：返回 `(g, s, t)` 使得 $s \cdot a + t \cdot b = g$ |
| `diophantine(polys, b)` | 多项式 CRT 求解器 |
| `p_adic_expansion(&self, p)` | 重复除法展开 |
| `pow(n)` | 多项式幂（快速幂） |
| `lcoeff()` | 首项系数（便捷） |
| `constant()` | 常数项（便捷） |
| `mul_coeff(c)` / `div_coeff(c)` | 标量乘/除所有系数 |
| `content()` | 系数的 GCD（现为 public） |

`SparseMultivariatePolynomial` 新增：

| 方法 | 说明 |
|---|---|
| `div_exact(&self, divisor)` | 精确多项式除法（无余数） |
| `degree_in(var_index)` | 指定变量的次数 |

---

## 有理重构

`rational_reconstruction` 函数从模表示恢复有理数 $\frac{n}{d}$：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::rational_reconstruction::rational_reconstruction;

// 3/7 mod 101: 7^{-1} mod 101 = 29, 所以 a = 3 * 29 mod 101 = 87
let m = Integer::from(101);
let a = Integer::from(87);
let result = rational_reconstruction(&a, &m);
assert!(result.is_some());
```

使用扩展 Euclidean 算法（连分数方法），适用于需要将结果提升回 $\mathbb{Q}$
的模方法算法。

---

## Karatsuba 乘法

稠密多项式乘法现在对 32 个或更多系数的多项式使用 **Karatsuba 算法**。
这将乘法复杂度从 $O(n^2)$ 降低到 $O(n^{1.585})$。

阈值（32）通过经验选择——低于此大小时，Karatsuba 额外加减法的开销超过
减少乘法次数的收益。

对于两个 500 次 $\mathbb{Z}$ 上的多项式，Karatsuba 比之前的 schoolbook
实现有显著加速。

---

## SymPy 迁移

| SymPy | oCAS |
|---|---|
| `sp.apart(expr, x)` | `apart(&num, &den)` |
| `sp.together(expr)` | `together(poly_part, &terms)` |
| `sp.resultant(a, b, x)` | `a.resultant(&b)` |
| `sp.Rational(n, d)` | `Rational::new(n, d)` |

---

## 限制

- 部分分式分解目前使用无平方分解，适用于任何 `EuclideanDomain`。
  完全分解为不可约因子需要 `IntegerDomain` 或 `FiniteField`。
- 多元部分分式（`apart_multivariate`）推迟到 0.13+，需要 Gröbner F4。
- FFT/NTT 乘法推迟到未来版本。
