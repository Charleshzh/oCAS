# Gröbner 基与理想

本章记录 oCAS 中 Gröbner 基计算和理想运算的完整 Rust API。所有算法通过统一入口 [`groebner_basis`] 访问，理想运算在 `ocas_poly::ideal` 模块中。

## 目录

- [`Algorithm`](#algorithm) — 算法选择枚举
- [`GroebnerBasis`](#groebnerbasis) — Gröbner 基结构体
- [`groebner_basis`](#groebner_basis) — 统一入口函数
- [`buchberger`](#buchberger) — Buchberger 算法
- [`f4::f4`](#f4f4) — F4 矩阵算法
- [`f5::f5`](#f5f5) — F5 签名算法
- [`fglm`](#fglm) — FGLM 换序算法
- [`HilbertSeries`](#hilbertseries) — Hilbert 级数
- [`hilbert_series`](#hilbert_series) — 计算 Hilbert 级数
- [`eliminate`](#eliminate) — 消元
- [理想运算](#理想运算) — `ideal_contains`、`ideal_sum`、`ideal_product` 等
- [`PrimaryComponent`](#primarycomponent) — 准素分量
- [`PolynomialSystemSolution`](#polynomialsystemsolution) — 多项式系统求解结果

---

## `Algorithm`

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Algorithm {
    Auto,
    F4,
    F5,
    Buchberger,
}
```

**功能**：Gröbner 基计算的算法选择器，传递给 [`groebner_basis`] 统一入口。

**变体**：

| 变体 | 说明 |
|---|---|
| `Auto`（默认） | 根据理想大小和结构自动选择算法。当前路由到 F4，后续将根据 cyclic-n 基准测试调整 F5 的切换阈值。 |
| `F4` | 强制使用 F4 矩阵算法（Faugère 1999）。批量处理 S-多项式归约为稀疏矩阵行操作，对较大理想显著快于 Buchberger。 |
| `F5` | 强制使用 F5 签名算法（Faugère 2002）。通过 syzygy 判据在矩阵构造前拒绝零归约器，对困难理想（如 cyclic-n）可实现数量级加速。 |
| `Buchberger` | 强制使用经典 Buchberger S-多项式迭代算法，适合小型理想。 |

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, groebner_basis, SparseMultivariatePolynomial};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::F4);
assert!(gb.is_groebner_basis());
```

**参见**：[`groebner_basis`](#groebner_basis)

---

## `GroebnerBasis`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerBasis<D: Domain, O: MonomialOrder> {
    pub basis: Vec<SparseMultivariatePolynomial<D, O>>,
}
```

**功能**：多项式理想的 Gröbner 基。`basis` 字段存储基的多项式列表。通过 `buchberger`、`f4::f4`、`f5::f5` 或 `groebner_basis` 构造，所有入口函数均返回**既约 Gröbner 基**（minimized + auto-reduced）。

**类型参数**：

| 参数 | 约束 | 说明 |
|---|---|---|
| `D` | `Domain` | 系数域（如 `RationalDomain`、`FiniteField`） |
| `O` | `MonomialOrder` | 单项式序（如 `Lex`、`Grevlex`） |

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `basis` | `Vec<SparseMultivariatePolynomial<D, O>>` | 基的多项式列表 |

### `GroebnerBasis::buchberger`

```rust
pub fn buchberger(ideal: &[SparseMultivariatePolynomial<D, O>]) -> Self
```

**功能**：使用 Buchberger 算法从生成元计算 Gröbner 基。内部过滤零多项式，应用 Buchberger 第一判据（首项互素时跳过 S-多项式），最大迭代 10000 次。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | 理想的生成元列表 |

**返回值**：`GroebnerBasis<D, O>` — 未最小化/未自归约的原始基。

**注意**：要求系数域支持精确除法（即为域）。除法失败时会 panic。便捷函数 [`buchberger`](#buchberger)（自由函数）会额外调用 `minimize().auto_reduce()`。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = GroebnerBasis::buchberger(&[f1, f2]);
assert!(gb.basis.len() >= 2);
```

**参见**：[`buchberger`](#buchberger)（便捷自由函数）、[`f4::f4`](#f4f4)、[`f5::f5`](#f5f5)

### `minimize`

```rust
pub fn minimize(mut self) -> Self
```

**功能**：最小化基——移除首项单项式可被其他元素首项整除的多项式。

**返回值**：`Self` — 最小化后的基（消费 self）。

**参见**：[`auto_reduce`](#auto_reduce)

### `auto_reduce`

```rust
pub fn auto_reduce(mut self) -> Self
```

**功能**：自归约基——将每个元素对其他元素归约并使其为首一多项式。按首项单项式升序处理，确保标准既约 Gröbner 基性质：任何基元素的单项式不被其他基元素的首项整除。

**返回值**：`Self` — 既约基（消费 self）。

**参见**：[`minimize`](#minimize)

### `is_groebner_basis`

```rust
pub fn is_groebner_basis(&self) -> bool
```

**功能**：验证此基确实是 Gröbner 基——检查所有 S-多项式是否归约为零。

**返回值**：`bool` — `true` 表示是合法的 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::Auto);
assert!(gb.is_groebner_basis());
```

**参见**：[`groebner_basis`](#groebner_basis)

### `reorder`

```rust
pub fn reorder<O2: MonomialOrder>(&self) -> GroebnerBasis<D, O2>
where
    D: 'static,
```

**功能**：改变 Gröbner 基的单项式序。将多项式重新解释为目标序并重新运行 F4。这是简单换序路径；对零维理想，使用 [`fglm`](#fglm) 获得更快的 $O(n \cdot D^3)$ 转换。

**类型参数**：

| 参数 | 说明 |
|---|---|
| `O2` | 目标单项式序 |

**返回值**：`GroebnerBasis<D, O2>` — 目标序下的 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::{Grevlex, Lex};
use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial, f4};

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
let gb_grevlex = gb_lex.reorder::<Grevlex>();
assert!(gb_grevlex.is_groebner_basis());
```

**参见**：[`fglm`](#fglm)（零维理想更快的换序方法）

---

## `groebner_basis`

```rust
pub fn groebner_basis<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    algo: Algorithm,
) -> GroebnerBasis<D, O>
```

**功能**：Gröbner 基计算的统一入口。根据 `algo` 参数选择算法（`Auto`/`F4`/`F5`/`Buchberger`）。零多项式在内部由各后端过滤。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | 理想的生成元列表 |
| `algo` | `Algorithm` | 算法选择 |

**返回值**：`GroebnerBasis<D, O>` — 既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, groebner_basis, SparseMultivariatePolynomial};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::Auto);
assert!(gb.is_groebner_basis());
```

**参见**：[`Algorithm`](#algorithm)、[`GroebnerBasis`](#groebnerbasis)

---

## `buchberger`

```rust
pub fn buchberger<D: Domain, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**功能**：便捷函数——计算 Gröbner 基并执行最小化和自归约（等价于 `GroebnerBasis::buchberger(ideal).minimize().auto_reduce()`）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | 理想的生成元列表 |

**返回值**：`GroebnerBasis<D, O>` — 既约 Gröbner 基。

**参见**：[`GroebnerBasis::buchberger`](#groebnerbasisbuchberger)、[`groebner_basis`](#groebner_basis)

---

## `f4::f4`

```rust
pub fn f4<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**功能**：使用 F4 矩阵算法计算 Gröbner 基（Faugère 1999）。将顺序 S-多项式归替换为批量稀疏矩阵行操作（Gaussian 消元）。对 `FiniteField` 系数域自动使用原生 ℤ_p 快速路径（`FpPoly` i64 模算术），避免 BigInt 开销。

**内部优化**：

- **Gebauer–Moeller 对筛选**：chain criterion + update criterion，严格整除守卫
- **`SimpCache`**：缓存多项式乘单项式的结果，避免重复计算
- **`DivisorIndex`**：基于支持位掩码的 O(1)-ish 归约器查找（变量索引 ≤ 63 时位掩码精确，>63 时退化为正确性保守的过滤）
- **双指针稀疏行减法** `sub_scaled_fp`：O(nnz) 复杂度
- **原生 ℤ_p 快速路径**：`FpPoly` 使用 i64 模算术，仅在读取输入和输出结果时进行 BigInt 转换

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | 理想的生成元列表 |

**返回值**：`GroebnerBasis<D, O>` — 既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::groebner::f4::f4;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = f4(&[f1, f2]);
assert!(!gb.basis.is_empty());
```

**参见**：[`Algorithm::F4`](#algorithm)、[`f5::f5`](#f5f5)

---

## `f5::f5`

```rust
pub fn f5<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**功能**：使用 F5 签名算法计算 Gröbner 基（Faugère 2002）。为每个多项式附加签名（`module_pos`, `monomial`），通过 syzygy 判据在矩阵构造前拒绝零归约器，对困难理想（如 cyclic-n）可实现数量级加速。

**算法要点**：

- **签名**：`(module_pos, monomial)` 记录多项式的历史——`module_pos` 是输入生成元的索引，`monomial` 是应用的单项式倍数。签名按 **pot**（position-over-term）序比较：先比较 module position（小者优先），再按单项式序 `O` 比较。
- **Syzygy 跟踪**：矩阵行归约为零时，其签名为 syzygy。未来签名是已知 syzygy 的单项式倍数的行将被立即跳过（F5 syzygy 判据）。
- **增量处理**：逐个处理生成元，每个新生成元触发一轮矩阵构造和归约。

当前实现提供通用域 F5 核心。对 `FiniteField` 系数域自动使用原生 ℤ_p 快速路径。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | 理想的生成元列表 |

**返回值**：`GroebnerBasis<D, O>` — 既约 Gröbner 基（与 F4 输出完全相同）。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::groebner::f5::f5;

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let gb = f5(&[f1, f2]);
assert!(gb.is_groebner_basis());
```

**参见**：[`Algorithm::F5`](#algorithm)、[`f4::f4`](#f4f4)

---

## `fglm`

```rust
pub fn fglm<D: Domain, O2: MonomialOrder>(
    gb: &GroebnerBasis<D, impl MonomialOrder>,
) -> Option<GroebnerBasis<D, O2>>
```

**功能**：FGLM 换序算法（Faugère–Gianni–Lazard–Mora 1993）。将**零维**理想的 Gröbner 基从一种单项式序转换为另一种，复杂度 $O(n \cdot D^3)$（$D$ 为 $R/I$ 的向量空间维数），对大型零维理想远快于重新运行 F4。

**算法步骤**：

1. 计算阶梯（staircase）——不被任何首项单项式整除的单项式集合
2. 按目标序遍历阶梯单项式
3. 对每个单项式计算当前基下的正规形
4. 若正规形与已见向量线性相关，从系数关系构造新基元素
5. 否则将其添加到已见向量集合

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `gb` | `&GroebnerBasis<D, impl MonomialOrder>` | 输入 Gröbner 基（必须是既约的） |

**返回值**：`Option<GroebnerBasis<D, O2>>` — `Some(gb)` 为转换后的基，`None` 表示理想非零维（阶梯无限）。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::{Grevlex, Lex};
use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial, f4};
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

**参见**：[`GroebnerBasis::reorder`](#reorder)（简单换序，适用于任意维数）、[`is_zero_dimensional`](#is_zero_dimensional)

---

## `HilbertSeries`

```rust
#[derive(Debug, Clone)]
pub struct HilbertSeries {
    pub numerator: Vec<i64>,
    pub denominator_power: usize,
}
```

**功能**：商环 $R/I$ 的 Hilbert 级数，表示为有理函数 $H(t) = N(t) / (1-t)^n$。分子 $N(t)$ 存储为系数向量（`numerator[i]` 为 $t^i$ 的系数），分母为 $(1-t)^n$（$n$ 为变量数）。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `numerator` | `Vec<i64>` | 分子 $N(t)$ 的系数（从常数项起） |
| `denominator_power` | `usize` | $(1-t)$ 的幂次（= 变量数 $n$） |

### `hilbert_function`

```rust
pub fn hilbert_function(&self, degree: usize) -> i64
```

**功能**：计算度数 $d$ 处的 Hilbert 函数值 $\dim_k (R/I)_d$。使用公式 $H(d) = [t^d] N(t) / (1-t)^n$，其中 $(1-t)^{-n}$ 中 $t^k$ 的系数为 $\binom{n+k-1}{k}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `degree` | `usize` | 度数 $d$ |

**返回值**：`i64` — $\dim_k (R/I)_d$ 的值。

### `dimension`

```rust
pub fn dimension(&self) -> usize
```

**功能**：计算 $R/I$ 的 Krull 维数。通过检查分子在 $t=1$ 处的零点阶数（逐次求导直到在 $t=1$ 处非零）得到。

**返回值**：`usize` — Krull 维数。零维理想返回 0。

### `degree`

```rust
pub fn degree(&self) -> i64
```

**功能**：计算射影簇的次数。对良构 Hilbert 级数，等于分子在除去维数因子后在 $t=1$ 处的值。

**返回值**：`i64` — 射影簇的次数。

### `hilbert_polynomial`

```rust
pub fn hilbert_polynomial(&self) -> Vec<f64>
```

**功能**：计算 Hilbert 多项式 $P(d)$ 的系数（$H(d) = P(d)$ 当 $d \gg 0$ 时），使用 Lagrange 插值。返回升幂排列的系数（`result[i]` 为 $d^i$ 的系数），多项式次数为 `self.dimension()`。

**返回值**：`Vec<f64>` — Hilbert 多项式系数（升幂排列）。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::groebner::hilbert::hilbert_series;

let d = RationalDomain;
// 理想 (x² - 1, y² - 1)：LM 理想 (x², y²)，分子 N(t) = 1 - 2t² + t⁴
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 2], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::F4);
let hs = hilbert_series(&gb);
assert_eq!(hs.dimension(), 0); // 零维：4 个解
assert_eq!(hs.degree(), 4); // 次数 = 解的个数
println!("维数: {}, 次数: {}", hs.dimension(), hs.degree());
for d in 0..10 {
    println!("H({}) = {}", d, hs.hilbert_function(d));
}
let hp = hs.hilbert_polynomial();
println!("Hilbert 多项式系数: {:?}", hp);
```

**参见**：[`hilbert_series`](#hilbert_series)

---

## `hilbert_series`

```rust
pub fn hilbert_series(
    gb: &GroebnerBasis<RationalDomain, Lex>,
) -> HilbertSeries
```

**功能**：从 Gröbner 基计算 $R/I$ 的 Hilbert 级数。使用 Macaulay 定理：$R/I$ 的 Hilbert 级数等于 $R/\text{LM}(I)$（首项理想）的 Hilbert 级数。通过容斥原理计算单项式理想 $\langle m_1, \dots, m_s \rangle$ 的 Hilbert 分子：$N(t) = \sum_k (-1)^k \sum_{|S|=k} t^{\deg \text{lcm}(S)}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `gb` | `&GroebnerBasis<RationalDomain, Lex>` | 有理数域上 Lex 序的 Gröbner 基 |

**返回值**：`HilbertSeries` — 商环的 Hilbert 级数。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::groebner::hilbert::hilbert_series;

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::F4);
let hs = hilbert_series(&gb);
assert!(hs.dimension() <= 2);
```

**参见**：[`HilbertSeries`](#hilbertseries)、[`is_zero_dimensional`](#is_zero_dimensional)

---

## `eliminate`

```rust
pub fn eliminate<D: Domain + 'static>(
    ideal: &[SparseMultivariatePolynomial<D, Lex>],
    elim_vars: usize,
    algo: Algorithm,
) -> GroebnerBasis<D, Lex>
```

**功能**：从理想中消元。返回 $I \cap k[x_{\text{elim\_vars}}, \dots, x_{n-1}]$ 的 Gröbner 基，即不包含前 `elim_vars` 个变量的多项式。使用 Lex 序——在 Lex 序下，理想的既约 Gröbner 基自动包含消元理想的生成元。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想的生成元（必须为 Lex 序） |
| `elim_vars` | `usize` | 要消去的变量数（消去 $x_0, \dots, x_{\text{elim\_vars}-1}$） |
| `algo` | `Algorithm` | Gröbner 基算法选择 |

**返回值**：`GroebnerBasis<D, Lex>` — 消元理想的 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{SparseMultivariatePolynomial, eliminate, Algorithm};

let d = RationalDomain;
// 理想: x + y + z, x*y + x*z 于 k[x,y,z]; 消去 x
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
    (vec![1, 0, 0], Rational::new(1, 1)),
    (vec![0, 1, 0], Rational::new(1, 1)),
    (vec![0, 0, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
    (vec![1, 1, 0], Rational::new(1, 1)),
    (vec![1, 0, 1], Rational::new(1, 1)),
]);
let elim = eliminate(&[f1, f2], 1, Algorithm::Auto);
// 结果应在 k[y,z] 中
for p in &elim.basis {
    assert!(p.degree_in(0) == 0, "消去的变量 x 不应出现");
}
```

**参见**：[`groebner_basis`](#groebner_basis)、[`ideal_quotient`](#ideal_quotient)

---

## 理想运算

所有理想运算位于 `ocas_poly::ideal` 模块，使用 `Lex` 序以保持与消元计算的一致性。

### `ideal_contains`

```rust
pub fn ideal_contains<D: Domain + 'static>(
    generators: &[SparseMultivariatePolynomial<D, Lex>],
    f: &SparseMultivariatePolynomial<D, Lex>,
    algo: Algorithm,
) -> bool
```

**功能**：测试 $f$ 是否属于生成元张成的理想。计算理想的 Gröbner 基并对 $f$ 归约，$f \in I$ 当且仅当余式为零。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想的生成元 |
| `f` | `&SparseMultivariatePolynomial<D, Lex>` | 待测试的多项式 |
| `algo` | `Algorithm` | Gröbner 基算法选择 |

**返回值**：`bool` — `true` 表示 $f \in I$。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_contains;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};

let d = RationalDomain;
let x = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let y = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 1], Rational::new(1, 1)),
]);
assert!(ideal_contains(&[x.clone(), y.clone()], &x, Algorithm::Auto));
assert!(!ideal_contains(&[y], &x, Algorithm::Auto));  // x ∉ ⟨y⟩
```

**参见**：[`groebner_basis`](#groebner_basis)

---

### `ideal_sum`

```rust
pub fn ideal_sum<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**功能**：两个理想的和 $I + J = \langle f_1, \dots, f_m, g_1, \dots, g_n \rangle$。合并生成元并计算 Gröbner 基。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $I$ 的生成元 |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $J$ 的生成元 |

**返回值**：`GroebnerBasis<D, Lex>` — $I + J$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_sum;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};

let d = RationalDomain;
let x = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let y = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 1], Rational::new(1, 1)),
]);
let gb = ideal_sum(&[x], &[y]);
// ⟨x⟩ + ⟨y⟩ = ⟨x, y⟩
assert!(gb.basis.len() >= 2);
```

**参见**：[`ideal_intersection`](#ideal_intersection)

---

### `ideal_product`

```rust
pub fn ideal_product<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**功能**：两个理想的积 $I \cdot J = \langle f_i \cdot g_j \rangle$。计算所有 $f_i g_j$ 的乘积并求 Gröbner 基。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $I$ 的生成元 |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $J$ 的生成元 |

**返回值**：`GroebnerBasis<D, Lex>` — $I \cdot J$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_product;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
let x = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let y = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 1], Rational::new(1, 1)),
]);
let gb = ideal_product(&[x], &[y]);
// ⟨x⟩ · ⟨y⟩ = ⟨xy⟩
assert_eq!(gb.basis.len(), 1);
```

**参见**：[`ideal_quotient`](#ideal_quotient)

---

### `ideal_quotient`

```rust
pub fn ideal_quotient<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**功能**：理想的商 $I : J = \{f : f \cdot g \in I, \forall g \in J\}$。对 $J$ 的每个生成元 $g$ 使用 Rabinowitsch 技巧：在扩展环 $k[x_1, \dots, x_n, w]$ 中计算 $\text{GB}(I \cup \{1 - w \cdot g\})$ 并消去 $w$，然后对结果求交集。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators_i` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $I$ 的生成元 |
| `generators_j` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $J$ 的生成元 |

**返回值**：`GroebnerBasis<D, Lex>` — $I : J$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_quotient;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
// ⟨x², xy⟩ : ⟨x⟩ = ⟨x⟩
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let g = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let gb = ideal_quotient(&[f1, f2], &[g]);
assert!(!gb.basis.is_empty());
```

**参见**：[`ideal_saturate`](#ideal_saturate)、[`ideal_contains`](#ideal_contains)

---

### `ideal_saturate`

```rust
pub fn ideal_saturate<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**功能**：理想的饱和 $I : J^\infty = \bigcup_k (I : J^k)$。迭代计算 $I : J$、$(I : J) : J$ 等，直到稳定。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators_i` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $I$ 的生成元 |
| `generators_j` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $J$ 的生成元 |

**返回值**：`GroebnerBasis<D, Lex>` — $I : J^\infty$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_saturate;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
// ⟨x²y, xy²⟩ : ⟨x⟩^∞ = ⟨y⟩
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 1], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 2], Rational::new(1, 1)),
]);
let g = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let gb = ideal_saturate(&[f1, f2], &[g]);
assert!(!gb.basis.is_empty());
```

**参见**：[`ideal_quotient`](#ideal_quotient)

---

### `ideal_intersection`

```rust
pub fn ideal_intersection<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**功能**：两个理想的交 $I \cap J$。使用辅助变量 $t$：$I \cap J = \langle t \cdot f_i, (1-t) \cdot g_j \rangle \cap k[x_1, \dots, x_n]$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $I$ 的生成元 |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | 理想 $J$ 的生成元 |

**返回值**：`GroebnerBasis<D, Lex>` — $I \cap J$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_intersection;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
let x = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
]);
let y = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 1], Rational::new(1, 1)),
]);
let gb = ideal_intersection(&[x], &[y]);
// ⟨x⟩ ∩ ⟨y⟩ = ⟨xy⟩
assert_eq!(gb.basis.len(), 1);
```

**参见**：[`ideal_sum`](#ideal_sum)、[`eliminate`](#eliminate)

---

### `ideal_radical`

```rust
pub fn ideal_radical(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> GroebnerBasis<RationalDomain, Lex>
```

**功能**：计算理想的根式 $\sqrt{I}$。

- **零维理想**：通过 Lex GB 中一元多项式的无平方分解计算（$\sqrt{I}$ 的生成元由各变量的无平方一元多项式和非单变量基元素组成）。
- **正维理想**：使用 Jacobian 饱和方法（简化 Kemper 算法）：$\sqrt{I} = I : h^\infty$。当前实现以启发式选取 $h$——取总次数最小的非平凡偏导数（并非精确 GCD）。若所有偏导数均为常数或零（Jacobian 平凡），保守返回原 GB。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 理想的生成元 |

**返回值**：`GroebnerBasis<RationalDomain, Lex>` — $\sqrt{I}$ 的既约 Gröbner 基。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::ideal_radical;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
// √(x², xy) = (x)
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let rad = ideal_radical(&[f1, f2]);
assert!(!rad.basis.is_empty());
```

**参见**：[`primary_decomposition`](#primary_decomposition)

---

### `is_zero_dimensional`

```rust
pub fn is_zero_dimensional(gb: &GroebnerBasis<RationalDomain, Lex>) -> bool
```

**功能**：检查理想是否为零维。理想零维当且仅当对每个变量 $x_i$，GB 中某个首项单项式是纯幂 $x_i^N$（等价于阶梯/标准单项式集合有限）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `gb` | `&GroebnerBasis<RationalDomain, Lex>` | 有理数域上 Lex 序的 Gröbner 基 |

**返回值**：`bool` — `true` 表示零维。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, GroebnerBasis, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::ideal::is_zero_dimensional;

let d = RationalDomain;
// x² - 1, y - x → 零维
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![0, 1], Rational::new(1, 1)),
    (vec![1, 0], Rational::new(-1, 1)),
]);
let gb = groebner_basis(&[f1, f2], Algorithm::F4);
assert!(is_zero_dimensional(&gb));
```

**参见**：[`solve_polynomial_system`](#solve_polynomial_system)、[`fglm`](#fglm)

---

### `solve_polynomial_system`

```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

**功能**：求解零维多项式系统。将 GB 转换为 Lex 序，提取各变量的一元多项式，通过回代求解。返回 `PolynomialSystemSolution` 枚举，区分零维（有限解）、正维（无穷解集）和空集（$\langle 1 \rangle$）。

**求解过程**：

1. 计算 Gröbner 基
2. 检查是否为 $\langle 1 \rangle$（空集）
3. 检查是否零维
4. 转换为 Lex 序的三角形分解
5. 从最后一个变量开始回代：对每个变量求一元多项式的实根（Sturm 定理隔离 + 精化至 $10^{-14}$），递归代入已知值

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `equations` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 方程组（多项式 = 0） |
| `algo` | `Algorithm` | Gröbner 基算法选择 |

**返回值**：`PolynomialSystemSolution` — 求解结果。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::ideal::{solve_polynomial_system, PolynomialSystemSolution};

let d = RationalDomain;
// x² + y² - 1, x - y → 解在 (±1/√2, ±1/√2)
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 2], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]);
let sol = solve_polynomial_system(&[f1, f2], Algorithm::Auto);
match sol {
    PolynomialSystemSolution::ZeroDimensional(z) => {
        assert_eq!(z.solutions.len(), 2);
    }
    _ => panic!("expected zero-dimensional"),
}
```

**参见**：[`PolynomialSystemSolution`](#polynomialsystemsolution)、[`is_zero_dimensional`](#is_zero_dimensional)

---

### `primary_decomposition`

```rust
pub fn primary_decomposition(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> Vec<PrimaryComponent>
```

**功能**：计算理想的准素分解。

- **零维理想**：对 Lex GB 中第一个变量的一元多项式进行因式分解，通过饱和分离准素分量（对每个因子 $f_i$，计算 $I : (\prod_{j \neq i} f_j)^\infty$）。
- **正维理想**：保守返回单一准素分量（原 GB 自身作为 primary 和 prime）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 理想的生成元 |

**返回值**：`Vec<PrimaryComponent>` — 准素分量列表。

**示例**：

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::primary_decomposition;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
// (x², xy) 为正维理想 → 保守返回单一分量
// （理想论分解为 (x², xy) = (x) ∩ (x², y)）
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let decomp = primary_decomposition(&[f1, f2]);
assert!(decomp.len() >= 1);
```

**参见**：[`PrimaryComponent`](#primarycomponent)、[`ideal_radical`](#ideal_radical)

---

### `is_prime_ideal`

```rust
pub fn is_prime_ideal(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> bool
```

**功能**：测试理想是否为素理想。

- **零维理想**：检查 Lex GB 中的一元多项式是否不可约（对次数 ≤ 3 的多项式，使用有理根定理检测）。
- **正维理想**：保守返回 `false`（完整实现需要检查簇的不可约性，尚未实现）。

**注意**：这是一个保守近似——永远不会返回假阳性（非素理想报告为素），只可能返回假阴性（素理想报告为非素）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 理想的生成元 |

**返回值**：`bool` — `true` 表示（已确认的）素理想。

**参见**：[`is_primary_ideal`](#is_primary_ideal)、[`primary_decomposition`](#primary_decomposition)

---

### `is_primary_ideal`

```rust
pub fn is_primary_ideal(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> bool
```

**功能**：测试理想是否为准素理想。准素理想恰有一个伴随素理想，即 `primary_decomposition` 返回的分量数 ≤ 1。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | 理想的生成元 |

**返回值**：`bool` — `true` 表示准素理想。

**参见**：[`primary_decomposition`](#primary_decomposition)、[`is_prime_ideal`](#is_prime_ideal)

---

## `PrimaryComponent`

```rust
#[derive(Debug, Clone)]
pub struct PrimaryComponent {
    pub primary: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
    pub prime: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
}
```

**功能**：理想的准素分量——一个准素理想及其伴随素理想（即准素理想的根式）。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `primary` | `Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>` | 准素理想的生成元 |
| `prime` | `Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>` | 伴随素理想（根式）的生成元 |

**参见**：[`primary_decomposition`](#primary_decomposition)

---

## `PolynomialSystemSolution`

```rust
#[derive(Debug, Clone)]
pub enum PolynomialSystemSolution {
    ZeroDimensional(ZeroDimSolutions),
    PositiveDimensional(GroebnerBasis<RationalDomain, Lex>),
    Empty,
}
```

**功能**：多项式系统求解的结果，按解集类型区分。

**变体**：

| 变体 | 说明 |
|---|---|
| `ZeroDimensional(ZeroDimSolutions)` | 有限个实解（零维理想） |
| `PositiveDimensional(GroebnerBasis<RationalDomain, Lex>)` | 无穷解集（正维理想），返回 Lex 序的 Gröbner 基 |
| `Empty` | 无解（理想为 $\langle 1 \rangle$） |

**参见**：[`solve_polynomial_system`](#solve_polynomial_system)

### `ZeroDimSolutions`

```rust
#[derive(Debug, Clone)]
pub struct ZeroDimSolutions {
    pub solutions: Vec<RealSolution>,
    pub vector_space_dimension: usize,
}
```

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `solutions` | `Vec<RealSolution>` | 找到的实解列表 |
| `vector_space_dimension` | `usize` | 商环 $k[x_1,\dots,x_n]/I$ 的向量空间维数（$\mathbb{C}$ 上解的总数，计重数） |

### `RealSolution`

```rust
#[derive(Debug, Clone)]
pub struct RealSolution {
    pub values: Vec<f64>,
    pub multiplicity: usize,
}
```

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `values` | `Vec<f64>` | 各变量的值（每个变量一个） |
| `multiplicity` | `usize` | 解的代数重数 |

---

## 模块路径速查

| 函数/类型 | 完整路径 |
|---|---|
| `Algorithm` | `ocas_poly::Algorithm` |
| `GroebnerBasis` | `ocas_poly::GroebnerBasis` |
| `groebner_basis` | `ocas_poly::groebner_basis` |
| `buchberger` | `ocas_poly::buchberger` |
| `eliminate` | `ocas_poly::eliminate` |
| `f4::f4` | `ocas_poly::groebner::f4::f4` |
| `f5::f5` | `ocas_poly::groebner::f5::f5` |
| `fglm` | `ocas_poly::groebner::fglm::fglm` |
| `HilbertSeries` | `ocas_poly::groebner::hilbert::HilbertSeries` |
| `hilbert_series` | `ocas_poly::groebner::hilbert::hilbert_series` |
| `ideal_contains` | `ocas_poly::ideal::ideal_contains` |
| `ideal_sum` | `ocas_poly::ideal::ideal_sum` |
| `ideal_product` | `ocas_poly::ideal::ideal_product` |
| `ideal_quotient` | `ocas_poly::ideal::ideal_quotient` |
| `ideal_saturate` | `ocas_poly::ideal::ideal_saturate` |
| `ideal_intersection` | `ocas_poly::ideal::ideal_intersection` |
| `ideal_radical` | `ocas_poly::ideal::ideal_radical` |
| `primary_decomposition` | `ocas_poly::ideal::primary_decomposition` |
| `is_zero_dimensional` | `ocas_poly::ideal::is_zero_dimensional` |
| `solve_polynomial_system` | `ocas_poly::ideal::solve_polynomial_system` |
| `is_prime_ideal` | `ocas_poly::ideal::is_prime_ideal` |
| `is_primary_ideal` | `ocas_poly::ideal::is_primary_ideal` |
| `PrimaryComponent` | `ocas_poly::ideal::PrimaryComponent` |
| `PolynomialSystemSolution` | `ocas_poly::ideal::PolynomialSystemSolution` |
| `ZeroDimSolutions` | `ocas_poly::ideal::ZeroDimSolutions` |
| `RealSolution` | `ocas_poly::ideal::RealSolution` |
