# 因式分解

本章记录 oCAS 中多项式因式分解的完整 API，覆盖一元 $\mathbb{Z}[x]$、一元 $\mathbb{F}_p[x]$、多元 $\mathbb{Z}[x_1,\dots,x_n]$、多元 $\mathbb{F}_p[x_1,\dots,x_n]$、代数数域 $\mathbb{Q}(\alpha)[x]$ 上的因式分解，以及有理函数算术与结式。

相关模块：

| 模块路径 | 功能 |
|---|---|
| `ocas_poly::factor` | 因式分解顶层入口 |
| `ocas_poly::factor::hensel` | Hensel 提升与 Zassenhaus 重组（$\mathbb{Z}[x]$） |
| `ocas_poly::factor::finite_field` | $\mathbb{F}_p[x]$ 因式分解（Cantor–Zassenhaus + Berlekamp） |
| `ocas_poly::factor::multivariate` | 双变量 Hensel 提升（Wang 算法） |
| `ocas_poly::factor::eez` | 多元 EEZ Hensel 提升（Wang + 首项系数重建） |
| `ocas_poly::factor::algebraic` | 代数数域因式分解（Trager 算法） |
| `ocas_poly::gcd` | GCD（伪余式 + Euclid 算法） |
| `ocas_poly::gcd::modular` | 模 GCD（Brown 1971） |
| `ocas_poly::resultant` | 结式（Brown PRS） |
| `ocas_poly::rational` | 有理函数分式域 |
| `ocas_calc::partial_fraction` | 部分分式分解 |

---

## 类型别名

```rust
/// 无平方分解结果：(因子, 重数) 对列表。
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

/// 完全因式分解结果：(不可约因子, 重数) 对列表。
pub type Factors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;
```

---

## 一、一元整数多项式因式分解 $\mathbb{Z}[x]$

### DenseUnivariatePolynomial::factor

**签名**：`pub fn factor(&self) -> Factors<IntegerDomain>`

**功能**：将一元整数多项式分解为不可约因子（在 $\mathbb{Q}$ 上不可约的本原多项式）及其重数。

**参数**：无（`self` 为待分解多项式）。

**返回值**：`Vec<(DenseUnivariatePolynomial<IntegerDomain>, usize)>` —— 不可约因子与重数对列表。各因子为本原多项式且首项系数为正。所有因子的乘积（取各自重数幂）等于 `self` 的本原部分 `primitive_part(self)`；输入必须为本原多项式（content = 1），否则内容会被丢弃。

**算法流程**：
1. 输入须为本原多项式（content = 1，可用 `primitive_part()` 预处理）。
2. 无平方分解（Yun 算法）：$\gcd(f, f')$ 分离重数。
3. 对每个无平方分量调用 `factor_square_free`（内部对非首一情形使用首项系数变换，再调用 `factor_square_free_monic`）：
   - 选取素数 $p$（$p \nmid \mathrm{lc}(f)$，$f \bmod p$ 无平方）。
   - 在 $\mathbb{F}_p[x]$ 中分解为首一不可约因子（Cantor–Zassenhaus）。
   - 计算 Mignotte 界 $B = 2^n \|f\|_2$。
   - 线性 Hensel 提升 $p \to p^k$（$p^k > 2B$）。
   - Zassenhaus 子集重组：枚举提升因子的子集，试除验证。
4. 非首一多项式通过首项系数变换 $a^{d-1} f(x/a)$ 处理。

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
// 每个因子次数为 1，重数为 1
for (g, m) in &factors {
    assert_eq!(g.degree(), Some(1));
    assert_eq!(*m, 1);
}
// 输出：[(x - 1, 1), (x + 1, 1)]
```

**参见**：[`square_free_factorization`](#square_free_factorization)、[`factor_over_finite_field`](#factor_over_finite_field)、[多项式 GCD 与因式分解](../math/poly-gcd-factoring.md)

---

### DenseUnivariatePolynomial::square_free_factorization

**签名**：`pub fn square_free_factorization(&self) -> SquareFreeFactors<D>`

**功能**：计算多项式的无平方分解。返回互不相同的无平方因子及其在原多项式中的重数。

**参数**：无（`self` 为待分解多项式）。

**返回值**：`Vec<(DenseUnivariatePolynomial<D>, usize)>` —— 无平方因子与重数对列表。

**算法**：Yun 算法（特征 0）。令 $g = \gcd(f, f')$，$w = f/g$，迭代 $h = \gcd(w, g)$，$z = w/h$ 收集重数 $k$ 的因子。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// (x+1)^2 * (x-1) = x^3 + x^2 - x - 1
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(-1), Integer::from(1), Integer::from(1),
]);
let factors = p.square_free_factorization();
assert_eq!(factors.len(), 2);
// 输出：[(x - 1, 1), (x + 1, 2)]（Yun 算法先收集重数 1 的因子）
```

**参见**：[`is_square_free`](#is_square_free)、[`factor`](#denseunivariatepolynomialfactor)

---

### DenseUnivariatePolynomial::is_square_free

**签名**：`pub fn is_square_free(&self) -> bool`

**功能**：判断多项式是否无平方（即 $\gcd(f, f') = 1$）。

**参数**：无。

**返回值**：`bool` —— 无平方返回 `true`。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x^2 - 1 无平方
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
assert!(p.is_square_free());

// (x+1)^2 不是无平方的
let q = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
assert!(!q.is_square_free());
```

---

## 二、一元有限域多项式因式分解 $\mathbb{F}_p[x]$

### factor_over_finite_field

**签名**：`pub fn factor_over_finite_field(f: &FpPoly) -> Vec<(FpPoly, usize)>`

**功能**：将 $\mathbb{F}_p[x]$ 中的多项式分解为首一不可约因子及重数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&DenseUnivariatePolynomial<FiniteField>` | 待分解多项式 |

**返回值**：`Vec<(FpPoly, usize)>` —— 首一不可约因子与重数对列表。首项系数作为常数因子（重数 1）单独列出。

**算法流程**：
1. 无平方分解（Musser/Bernardin 算法，处理特征 $p$ 的 $p$ 次根情况）。
2. 对每个无平方分量：
   - **DDF**（分次分解）：利用 Frobenius 映射 $x \mapsto x^{p^d} \bmod f$，将因子按次数分组。
   - **EDF**（等次分解）：对奇特征 $p$，随机选取 $a$，计算 $\gcd(f, a^{(p^d-1)/2} - 1)$；对特征 2，使用迹映射。
3. 小素数时也可使用 Berlekamp 算法（Frobenius 矩阵 $Q^T - I$ 的核空间）。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::factor::finite_field::factor_over_finite_field;

let f = FiniteField::new(BigInt::from(5));
// x^2 - 1 = (x-1)(x+1) over F_5
let p = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(4), f.element(0), f.element(1)]);
let factors = factor_over_finite_field(&p);
let linear_count = factors.iter()
    .filter(|(g, _)| g.degree() == Some(1)).count();
assert_eq!(linear_count, 2);
// 输出：两个一次首一因子
```

**参见**：[Cantor–Zassenhaus 算法详解](../algorithms/factorization.md)

---

### DenseUnivariatePolynomial::factor (FiniteField)

**签名**：`pub fn factor(&self) -> Factors<FiniteField>`

**功能**：`DenseUnivariatePolynomial<FiniteField>` 上的 `factor()` 方法，等价于调用 `factor_over_finite_field`。

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

---

### berlekamp

**签名**：`pub fn berlekamp(f: &FpPoly) -> Vec<FpPoly>`

**功能**：Berlekamp 因式分解算法，适用于小素数域。构造 Frobenius 矩阵 $Q$（$Q[i][j]$ 为 $x^{ip} \bmod f$ 的 $x^j$ 系数），求 $Q^T - I$ 的核空间；每个非零核向量 $v$ 满足 $v^p \equiv v \pmod{f}$，通过 $\gcd(f, v - a)$（$a \in \mathbb{F}_p$）分裂因子。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&FpPoly` | 首一无平方多项式 |

**返回值**：`Vec<FpPoly>` —— 首一不可约因子列表。

**参见**：Berlekamp (1970)

---

### cantor_zassenhaus

**签名**：`pub fn cantor_zassenhaus(f: &FpPoly) -> Vec<FpPoly>`

**功能**：Cantor–Zassenhaus 因式分解算法。先 DDF（分次分解），再 EDF（等次分解）。返回首一无平方多项式的首一不可约因子列表。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&FpPoly` | 首一无平方多项式 |

**返回值**：`Vec<FpPoly>` —— 首一不可约因子列表。

**参见**：Cantor & Zassenhaus (1981)

---

### poly_pow_mod

**签名**：`pub fn poly_pow_mod(base: &FpPoly, exp: &BigInt, modulus: &FpPoly) -> FpPoly`

**功能**：计算 $\text{base}^{\text{exp}} \bmod \text{modulus}$，使用快速幂（重复平方），每次乘法后取模。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `base` | `&FpPoly` | 底数多项式 |
| `exp` | `&BigInt` | 指数（非负） |
| `modulus` | `&FpPoly` | 模多项式（非零） |

**返回值**：`FpPoly` —— 模幂结果。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::factor::finite_field::poly_pow_mod;

let f = FiniteField::new(BigInt::from(7));
let m = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(1), f.element(0), f.element(1)]); // x^2 + 1
let base = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(0), f.element(1)]); // x
// x^2 mod (x^2+1) = -1 = 6 in F_7
let r = poly_pow_mod(&base, &BigInt::from(2), &m);
assert_eq!(r.coeff(0).cloned(), Some(f.element(6)));
// 输出：6（即 -1 mod 7）
```

---

## 三、多元整数多项式因式分解 $\mathbb{Z}[x_1,\dots,x_n]$

### SparseMultivariatePolynomial::factor (IntegerDomain, Lex)

**签名**：`pub fn factor(&self) -> Vec<(Self, usize)>`

**功能**：将稀疏多元整数多项式分解为不可约因子及重数。

**参数**：无（`self` 为待分解的 `SparseMultivariatePolynomial<IntegerDomain, Lex>`）。

**返回值**：`Vec<(SparseMultivariatePolynomial<IntegerDomain, Lex>, usize)>` —— 不可约因子与重数对列表。

**算法选择**：
- **$n \geq 3$ 变量**：使用 EEZ Hensel 提升（Wang 首项系数重建 + p-adic 系数 Hensel 提升 + Zassenhaus 重组）。
- **双变量 + 常数首项系数**：使用双变量 Hensel 提升（评估–提升路径）。
- **双变量 + 非常数首项系数**：回退到 EEZ 路径。

**算法流程**（EEZ 路径，$n \geq 3$）：
1. 提取主变量内容，无平方分解。
2. 选取样本点 $(a_1, \dots, a_n)$，使单变量像 $f(x_0, a_1, \dots)$ 次数保持且无平方。
3. 在 $\mathbb{Z}$ 上分解单变量像。
4. Wang 首项系数重建：将总首项系数的不可约因子分配给各单变量因子。
5. EEZ 逐变量 Hensel 提升：依次恢复 $x_1, x_2, \dots$，每步解多元 Diophantine 方程。
6. p-adic 系数 Hensel 提升：从模 $p$ 提升到足够大的 $p^k$。
7. Zassenhaus 子集重组：枚举模因子子集，试除验证。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

// (x^2 + y + 1)(x + y + 2)  在 Z[x,y] 中
let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 2,
    vec![
        (vec![3, 0], Integer::from(1)),   // x^3
        (vec![2, 1], Integer::from(1)),   // x^2 y
        (vec![2, 0], Integer::from(2)),   // 2x^2
        (vec![1, 1], Integer::from(1)),   // xy
        (vec![1, 0], Integer::from(1)),   // x
        (vec![0, 2], Integer::from(1)),   // y^2
        (vec![0, 1], Integer::from(3)),   // 3y
        (vec![0, 0], Integer::from(2)),   // 2
    ],
);
let factors = f.factor();
assert!(factors.len() >= 2);
// 输出：[(x^2 + y + 1, 1), (x + y + 2, 1)]
```

**参见**：[Wang EEZ 算法](../math/poly-gcd-factoring.md)、[`bivariate_factor_z`](#bivariate_factor_z)

---

### bivariate_factor_z

**签名**：`pub fn bivariate_factor_z(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<(ZMPoly, usize)>`

**功能**：将双变量整数多项式分解为不可约因子及重数。要求 $x$ 的首项系数为非零整数常数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<IntegerDomain, Lex>` | 待分解多项式 |
| `x_var` | `usize` | 主变量索引 |
| `y_var` | `usize` | 次变量索引 |

**返回值**：`Vec<(ZMPoly, usize)>` —— 不可约因子与重数对。

**算法**：
1. 双变量无平方分解（启发式双变量 GCD）。
2. 选取 $y = \alpha$ 使单变量像 $f(x, \alpha)$ 无平方且因子数最少。
3. Hensel 提升：从 $f(x, \alpha)$ 的因子提升回双变量，使用 Taylor 展开逐次修正。

**参见**：Wang (1978)

---

### multivariate_factor_z

**签名**：`pub fn multivariate_factor_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)>`

**功能**：EEZ Hensel 提升多元因式分解入口。支持非恒定首项系数（通过 Wang 首项系数重建）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<IntegerDomain, Lex>` | 待分解多项式 |

**返回值**：`Vec<(ZmPoly, usize)>` —— 不可约因子与重数对。

**参见**：[Wang EEZ 算法详解](../math/poly-gcd-factoring.md)

---

## 四、多元有限域多项式因式分解 $\mathbb{F}_p[x_1,\dots,x_n]$

### SparseMultivariatePolynomial::factor (FiniteField, Lex)

**签名**：`pub fn factor(&self) -> Vec<(Self, usize)>`

**功能**：将稀疏多元有限域多项式分解为不可约因子及重数。

**参数**：无。

**返回值**：`Vec<(SparseMultivariatePolynomial<FiniteField, Lex>, usize)>` —— 不可约因子与重数对列表。

**算法选择**：
- **$n \geq 3$ 变量**：EEZ Hensel 提升。
- **双变量**：评估–Hensel 路径。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

let fp = FiniteField::new(BigInt::from(7));
// x*y + 1 over F_7[x,y] — 已不可约
let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    fp.clone(), 2,
    vec![
        (vec![1, 1], fp.element(1)),  // xy
        (vec![0, 0], fp.element(1)),  // 1
    ],
);
let factors = f.factor();
assert_eq!(factors.len(), 1);
assert_eq!(factors[0].1, 1); // 重数 1
```

---

### bivariate_factor_fp

**签名**：`pub fn bivariate_factor_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<(FpMPoly, usize)>`

**功能**：将双变量 $\mathbb{F}_p$ 多项式分解为不可约因子及重数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<FiniteField, Lex>` | 待分解多项式 |
| `x_var` | `usize` | 主变量索引 |
| `y_var` | `usize` | 次变量索引 |

**返回值**：`Vec<(FpMPoly, usize)>` —— 不可约因子与重数对。

---

### multivariate_factor_fp

**签名**：`pub fn multivariate_factor_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)>`

**功能**：多元 $\mathbb{F}_p$ 因式分解入口，支持非恒定首项系数（Wang 首项系数重建 + EEZ 提升）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<FiniteField, Lex>` | 待分解多项式 |

**返回值**：`Vec<(FpMPoly, usize)>` —— 不可约因子与重数对。

---

## 五、代数数域因式分解 $\mathbb{Q}(\alpha)[x]$

### DenseUnivariatePolynomial::factor (AlgebraicNumberField)

**签名**：`pub fn factor(&self) -> Factors<AlgebraicNumberField>`

**功能**：将代数数域 $\mathbb{Q}(\alpha)$ 上的多项式分解为首一不可约因子及重数（Trager 算法）。

**参数**：无（`self` 为 `DenseUnivariatePolynomial<AlgebraicNumberField>`）。

**返回值**：`Vec<(DenseUnivariatePolynomial<AlgebraicNumberField>, usize)>` —— 首一不可约因子与重数对。各因子的乘积（取重数幂）等于 `self` 除以首项系数（$K$ 的单位）。

**算法流程**（Trager）：
1. 无平方分解（Yun 算法，使用模 GCD `gcd_anf` 避免系数爆炸）。
2. 对每个无平方分量调用 `factor_square_free_anf`：
   - **Trager 平移**：找 $s \geq 0$ 使 $f(x - s\alpha)$ 的范数无平方。
   - **范数计算**：$\operatorname{Res}_\alpha(m(\alpha), f(x, \alpha))$，通过求值–插值。
   - **$\mathbb{Q}$ 上分解范数**：使用 Hensel 路径。
   - **恢复 $K$ 上因子**：对每个范数因子 $g_i$ 计算 $\gcd_K(f, g_i(\alpha))$。
3. GCD 使用模方法（Encarnación）：映射到 $\mathrm{GF}(p^d)$，CRT 合并，有理重构，试除验证。

**示例**：

```rust
use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

// 构造 Q(√2)，极小多项式 x^2 - 2
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
// x^2 - 2 在 Q(√2) 上可分解为 (x - √2)(x + √2)
let f = DenseUnivariatePolynomial::from_coeffs(
    field.clone(),
    vec![
        field.from_base(Rational::new(-2, 1)),
        field.zero(),
        field.one(),
    ],
);
let factors = f.factor();
assert_eq!(factors.len(), 2);
assert!(factors.iter().all(|(g, m)| *m == 1 && g.degree() == Some(1)));
// 输出：[(x - √2, 1), (x + √2, 1)]
```

**参见**：[代数数域与 Galois 理论](../math/algebraic-number-fields.md)、[Trager 算法](../math/poly-gcd-factoring.md)

---

### norm_with_shift

**签名**：`pub(crate) fn norm_with_shift(field: &AlgebraicNumberField, f: &UP<AlgebraicNumberField>) -> Option<(u64, UP<AlgebraicNumberField>, UP<RationalDomain>)>`

**功能**：Trager 平移：找 $s \geq 0$ 使 $f(x - s\alpha)$ 的范数在 $\mathbb{Q}$ 上无平方。返回 `(s, g, norm)`。最多尝试 `MAX_TRAGER_SHIFTS`（100）次。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `field` | `&AlgebraicNumberField` | 代数数域 $\mathbb{Q}(\alpha)$ |
| `f` | `&DenseUnivariatePolynomial<AlgebraicNumberField>` | 待分解多项式 |

**返回值**：`Option<(u64, UP<AlgebraicNumberField>, UP<RationalDomain>)>` —— `(s, f(x - s\alpha), \operatorname{Norm}(f(x - s\alpha)))`。

---

### gcd_anf

**签名**：`pub(crate) fn gcd_anf(field: &AlgebraicNumberField, a: &UP<AlgebraicNumberField>, b: &UP<AlgebraicNumberField>) -> UP<AlgebraicNumberField>`

**功能**：代数数域上两个一元多项式的 GCD。使用模方法（Encarnación）：映射到 $\mathrm{GF}(p^d)$，组合首一模 GCD（CRT），有理重构系数，试除验证。最多使用 `MAX_ANF_GCD_PRIMES`（1000）个素数，超出后回退到稠密 Euclid GCD。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `field` | `&AlgebraicNumberField` | 代数数域 |
| `a` | `&UP<AlgebraicNumberField>` | 第一个多项式 |
| `b` | `&UP<AlgebraicNumberField>` | 第二个多项式 |

**返回值**：首一 GCD 多项式。

---

## 六、多项式 GCD

### DenseUnivariatePolynomial::gcd

**签名**：`pub fn gcd(&self, other: &Self) -> Self`

**功能**：计算两个一元多项式的最大公因子。对非域系数使用伪余式 Euclid 算法。结果为本原多项式。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Self` | 另一个多项式 |

**返回值**：本原 GCD 多项式。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x^2 - 1 = (x-1)(x+1)
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
// x^2 + 2x + 1 = (x+1)^2
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
let g = a.gcd(&b);
assert_eq!(g.coeffs(), &[Integer::from(1), Integer::from(1)]);
// 输出：x + 1
```

---

### gcd_modular_z

**签名**：`pub fn gcd_modular_z(a: &ZPoly, b: &ZPoly) -> ZPoly`

**功能**：Brown 模 GCD 算法。在多个素数 $p$ 上计算 $\mathbb{F}_p[x]$ 中的首一 GCD，用 CRT 合并为 $\mathbb{Z}[x]$ 中的本原 GCD，试除验证。比伪余式 GCD 在大系数时高效得多。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `&DenseUnivariatePolynomial<IntegerDomain>` | 第一个多项式 |
| `b` | `&DenseUnivariatePolynomial<IntegerDomain>` | 第二个多项式 |

**返回值**：本原 GCD 多项式。

**算法细节**：
1. 选取素数 $p > 2^{30}$（避免小素数问题）。
2. 在 $\mathbb{F}_p[x]$ 中计算首一 GCD，缩放 $\gamma = \gcd(\mathrm{lc}(a), \mathrm{lc}(b))$。
3. CRT 合并（对称代表），丢弃"不幸素数"（GCD 次数高于真值）。
4. 精确试除验证。最多尝试 `MAX_PRIMES`（10000）个素数，超出后回退到伪余式。

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::gcd::modular::gcd_modular_z;

let d = IntegerDomain;
let i = |v: i64| Integer::from(v);
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);
let g = gcd_modular_z(&a, &b);
assert_eq!(g.coeffs(), &[i(1), i(1)]); // x + 1
// 输出：x + 1
```

**参见**：[模 GCD 算法（Brown 1971）](../math/poly-gcd-factoring.md)

---

### DenseUnivariatePolynomial::content

**签名**：`pub fn content(&self) -> D::Element`

**功能**：计算多项式的内容（所有系数的 GCD）。零多项式的内容为零。

---

### DenseUnivariatePolynomial::primitive_part

**签名**：`pub fn primitive_part(&self) -> Self`

**功能**：返回本原部分（多项式除以内容）。

---

## 七、结式

### DenseUnivariatePolynomial::resultant

**签名**：`pub fn resultant(&self, other: &Self) -> D::Element`

**功能**：使用 Brown PRS（多项式余式序列）算法计算两个一元多项式的结式 $\operatorname{Res}(a, b)$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Self` | 另一个多项式 |

**返回值**：`D::Element` —— 结式标量。当且仅当 $\gcd(a, b)$ 非常数时为零。

**算法**：子结式 PRS（Brown），每步精确除以 $\beta$（在 UFD 中子结式定理保证整除）。

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
// 输出：-1
```

**示例（有公共根）**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x - 1 和 x^2 - 1 共享根 x = 1
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(1),
]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
assert_eq!(a.resultant(&b), Integer::from(0));
// 输出：0（有公共因子 x - 1）
```

**参见**：[结式与子结式 PRS](../math/poly-gcd-factoring.md)

---

## 八、有理函数分式域

### RationalPolynomial

**签名**：

```rust
pub struct RationalPolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    pub numerator: SparseMultivariatePolynomial<D, O>,
    pub denominator: SparseMultivariatePolynomial<D, O>,
}
```

**功能**：多项式分式域元素 $\frac{n}{d}$，始终维持规范形式：分子与分母互素，分母首项系数为正（有序域）或 1（有限域）。

#### 构造方法

##### RationalPolynomial::new

```rust
pub fn new(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

创建有理多项式。**不**自动约简——需要规范形式时使用 `from_num_den`。

##### RationalPolynomial::from_num_den

```rust
pub fn from_num_den(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

从分子分母创建并自动约简（GCD 化简 + 分母首项系数归一化）。

##### RationalPolynomial::from_polynomial

```rust
pub fn from_polynomial(poly: SparseMultivariatePolynomial<D, O>) -> Self
```

从多项式创建（分母 = 1）。

#### 查询方法

| 方法 | 签名 | 说明 |
|---|---|---|
| `is_zero` | `&self -> bool` | 分子是否为零 |
| `is_one` | `&self -> bool` | 是否为 1/1 |
| `n_vars` | `&self -> usize` | 变量数 |
| `domain` | `&self -> &D` | 系数域引用 |

#### 算术运算

| 方法 | 签名 | 说明 |
|---|---|---|
| `add` | `(&self, &Self) -> Self` | 加法（分母 GCD 策略，减少中间膨胀） |
| `sub` | `(&self, &Self) -> Self` | 减法 |
| `mul` | `(&self, &Self) -> Self` | 乘法（交叉消去） |
| `div` | `(&self, &Self) -> Option<Self>` | 除法（除数分子为零时返回 `None`） |
| `neg` | `&self -> Self` | 取反 $-\frac{n}{d}$ |
| `inv` | `&self -> Option<Self>` | 取逆（分子为零时返回 `None`） |
| `pow` | `(&self, k: u32) -> Self` | 幂 $\left(\frac{n}{d}\right)^k$ |

**示例**：

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::rational::RationalPolynomial;
use ocas_poly::Grevlex;

let d = IntegerDomain;
let n_vars = 2;

// x / y
let num: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![1, 0], Integer::from(1))]);
let den: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![0, 1], Integer::from(1))]);
let f = RationalPolynomial::from_num_den(num, den);

// y / x
let num2: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![0, 1], Integer::from(1))]);
let den2: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![1, 0], Integer::from(1))]);
let g = RationalPolynomial::from_num_den(num2, den2);

// (x/y) * (y/x) = 1
let h = f.mul(&g);
assert!(h.is_one());
// 输出：1
```

**参见**：[RationalPolynomial 定义](../api/rust-polynomials.md)

---

## 九、部分分式分解

### apart

**签名**：

```rust
pub fn apart<D: EuclideanDomain>(
    num: &DenseUnivariatePolynomial<D>,
    den: &DenseUnivariatePolynomial<D>,
) -> (
    Option<DenseUnivariatePolynomial<D>>,
    Vec<PartialFractionTerm<D>>,
)
```

**功能**：对 $\frac{\text{num}}{\text{den}}$ 进行部分分式分解。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `num` | `&DenseUnivariatePolynomial<D>` | 分子多项式 |
| `den` | `&DenseUnivariatePolynomial<D>` | 分母多项式（非零） |

**返回值**：`(Option<poly_part>, Vec<PartialFractionTerm<D>>)` —— 可选的多项式部分与部分分式项列表，满足：

$$\frac{\text{num}}{\text{den}} = \text{poly\_part} + \sum_i \frac{\text{numer}_i}{\text{denom}_i^{\text{exp}_i}}$$

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;
// 1 / (x^2 - 1)
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_none()); // 真分式，无多项式部分
// x^2 - 1 = (x-1)(x+1) 无平方 → 输出按因子分组
// 输出：poly_part = None, terms = [...]
```

---

### PartialFractionTerm

**签名**：

```rust
pub struct PartialFractionTerm<D: EuclideanDomain> {
    pub numer: DenseUnivariatePolynomial<D>,  // 分子多项式
    pub denom: DenseUnivariatePolynomial<D>,  // 不可约（无平方）分母因子
    pub exp: usize,                            // 该因子在原分母中的重数
}
```

**功能**：部分分式分解中的单项，表示 $\frac{\text{numer}}{\text{denom}^{\text{exp}}}$。

---

### together

**签名**：

```rust
pub fn together<D: EuclideanDomain>(
    poly_part: Option<&DenseUnivariatePolynomial<D>>,
    terms: &[PartialFractionTerm<D>],
) -> (DenseUnivariatePolynomial<D>, DenseUnivariatePolynomial<D>)
```

**功能**：将多项式部分和部分分式项合并为单个有理函数 $\frac{n}{d}$。是 `apart` 的逆操作。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `poly_part` | `Option<&DenseUnivariatePolynomial<D>>` | 多项式部分（可选） |
| `terms` | `&[PartialFractionTerm<D>]` | 部分分式项列表 |

**返回值**：`(numerator, denominator)` —— 合并后的分子分母。

---

## 十、单变量多项式辅助方法

### DenseUnivariatePolynomial::extended_gcd_poly

**签名**：`pub fn extended_gcd_poly(&self, other: &Self) -> (Self, Self, Self)`

**功能**：扩展 Euclid 算法。返回 $(g, s, t)$ 满足 $s \cdot a + t \cdot b = g = \gcd(a, b)$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Self` | 另一个多项式 |

**返回值**：`(gcd, s, t)`。

---

### DenseUnivariatePolynomial::div_rem

**签名**：`pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>`

**功能**：带余除法。返回 `(商, 余式)`，使 `self = quotient * divisor + remainder`。除数为零时返回 `None`。

---

### DenseUnivariatePolynomial::pow

**签名**：`pub fn pow(&self, n: u32) -> Self`

**功能**：重复平方计算 $\text{self}^n$。

---

### DenseUnivariatePolynomial::p_adic_expansion

**签名**：`pub fn p_adic_expansion(&self, p: &Self) -> Vec<Self>`

**功能**：以多项式 $p$ 为基展开 $p$-adic 表示。返回系数列表 $[c_0, c_1, \dots]$ 使 $\text{self} = \sum_i c_i \cdot p^i$，其中 $\deg(c_i) < \deg(p)$。

---

## 十一、单项式工具函数

### monomial_divides

**签名**：`pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool`

**功能**：判断单项式 $b$ 是否整除 $a$（即 $a$ 是 $b$ 的倍式：$a_i \geq b_i$ 对所有 $i$）。

---

### monomial_lcm

**签名**：`pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]>`

**功能**：计算两个单项式的最小公倍数（逐分量取最大值）。

---

### monomial_are_coprime

**签名**：`pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool`

**功能**：判断两个单项式是否互素（无变量同时出现在两者中）。

---

## 十二、内部算法概览

### Hensel 提升与 Zassenhaus 重组

文件：`ocas-poly/src/factor/hensel.rs`

**mignotte_bound**

```rust
pub(crate) fn mignotte_bound(f: &ZPoly) -> Integer
```

Landau–Mignotte 界：对 $n$ 次多项式 $f$，其任何因子 $g$ 满足 $\|g\|_\infty \leq 2^n \|f\|_2$。

**factor_square_free_monic**

```rust
pub fn factor_square_free_monic(f: &ZPoly) -> Vec<ZPoly>
```

将首一无平方 $\mathbb{Z}[x]$ 多项式分解为首一不可约因子。

**factor_square_free**

```rust
pub fn factor_square_free(f: &ZPoly) -> Vec<ZPoly>
```

将无平方本原 $\mathbb{Z}[x]$ 多项式分解为不可约因子。非首一输入通过首项系数变换 $a^{d-1} f(x/a)$ 处理。

**factor_primitive**

```rust
pub fn factor_primitive(f: &ZPoly) -> Vec<(ZPoly, usize)>
```

将本原 $\mathbb{Z}[x]$ 多项式分解为不可约因子及重数（Yun 无平方分解 + `factor_square_free_monic`）。这是 `DenseUnivariatePolynomial::factor()` 的内部实现。

---

### EEZ Hensel 提升

文件：`ocas-poly/src/factor/eez.rs`

关键内部函数：

| 函数 | 说明 |
|---|---|
| `eez_lift` | 泛型 EEZ 提升（域上，首一） |
| `eez_lift_imposed` | 非首一 EEZ 提升（Wang 首项系数施加） |
| `eez_lift_z` | 整数 EEZ 提升（$\mathbb{Q}$ 上求解 + 积分性检查） |
| `coefficient_hensel_lift_z` | p-adic 系数 Hensel 提升 |
| `diophantine` | 递归多元 Diophantine 方程求解器 |
| `sparse_diophantine_fp` | 骨架插值稀疏 Diophantine 求解器 |
| `wang_reconstruct_lcoeffs` | Wang 首项系数重建 |
| `zassenhaus_multivariate` | 多元 Zassenhaus 子集重组 |

---

### 有限域因式分解

文件：`ocas-poly/src/factor/finite_field.rs`

| 函数 | 说明 |
|---|---|
| `factor_over_finite_field` | 顶层入口：无平方 + DDF + EDF |
| `distinct_degree_factorization` | DDF：按不可约因子次数分组 |
| `equal_degree_factorization` | EDF：随机分裂等次因子积 |
| `berlekamp` | Berlekamp 算法（小素数） |
| `cantor_zassenhaus` | DDF + EDF 组合 |
| `square_free_factorization_ff` | 特征 $p$ 无平方分解（Musser/Bernardin） |
| `pth_root_prime` | $\mathbb{F}_p$ 上的 $p$ 次根 |

---

### 代数数域因式分解

文件：`ocas-poly/src/factor/algebraic.rs`

| 函数 | 说明 |
|---|---|
| `factor_anf` | 顶层入口：无平方 + Trager |
| `factor_square_free_anf` | Trager：范数 → $\mathbb{Q}$ 分解 → GCD 恢复 |
| `norm_with_shift` | Trager 平移找无平方范数 |
| `norm_eval_interp` | 求值–插值计算范数 |
| `gcd_anf` | 模 GCD over $\mathbb{Q}(\alpha)$ |
| `square_free_anf` | Yun 无平方分解（模 GCD） |
| `factor_square_free_rationals` | $\mathbb{Q}$ 上无平方因子分解（清除分母 → $\mathbb{Z}$ Hensel） |

---

## 参见

- [多项式 GCD 与因式分解](../math/poly-gcd-factoring.md) — 数学理论与算法详解
- [代数数域与 Galois 理论](../math/algebraic-number-fields.md) — Trager 算法的数学背景
- [多项式 API](./rust-polynomials.md) — `DenseUnivariatePolynomial` 与 `SparseMultivariatePolynomial` 完整 API
- [系数域 API](./rust-domains.md) — `Integer`、`FiniteField`、`AlgebraicNumberField` 定义
