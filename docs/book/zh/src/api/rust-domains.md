# Rust API 参考：系数域

系数域是 oCAS 所有多项式、矩阵和求解器算法的基础抽象。每个域描述一组元素及其算术运算，通过 `Domain` 和 `EuclideanDomain` trait 统一接口。

**模块路径**：`ocas_domain`

**导入方式**：

```rust
use ocas_domain::{
    Domain, EuclideanDomain,
    Integer, IntegerDomain,
    Rational, RationalDomain,
    FiniteField, FiniteFieldElement,
    RealBall, RealBallDomain,
    Complex, ComplexDomain,
    DoubleF64, DoubleF64Domain,
    AlgebraicExtension, AlgebraicNumberField, AlgebraicElement,
};
use ocas_domain::assumptions::{Assumption, Assumptions, SymbolAssumptions};
```

---

## Domain trait

**签名**：

```rust
pub trait Domain: Clone + PartialEq + Eq + std::fmt::Debug + Sized {
    type Element: Clone + PartialEq + Eq + std::fmt::Debug + 'static;

    fn zero(&self) -> Self::Element;
    fn one(&self) -> Self::Element;
    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn neg(&self, a: &Self::Element) -> Self::Element;
    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element>;
    fn inv(&self, a: &Self::Element) -> Option<Self::Element>;
    fn is_zero(&self, a: &Self::Element) -> bool;
    fn is_one(&self, a: &Self::Element) -> bool;
    fn mul_assign(&self, a: &mut Self::Element, b: &Self::Element);
    fn sub_mul_assign(&self, a: &mut Self::Element, b: &Self::Element, c: &Self::Element);
    fn pow(&self, a: &Self::Element, n: u64) -> Self::Element;
    fn cast_u64(&self, n: u64) -> Self::Element;
}
```

**功能**：系数域的核心 trait。域对象本身可以携带参数（如有限域的模数），因此所有运算通过 `&self` 进行。这与 Flint、SymPy `Domain` 等 CAS 库的"域对象"模式一致。

**关联类型**：

| 关联类型 | 约束 | 说明 |
|---|---|---|
| `Element` | `Clone + PartialEq + Eq + Debug + 'static` | 域中元素的类型 |

**方法**：

| 方法 | 签名 | 说明 |
|---|---|---|
| `zero` | `fn zero(&self) -> Self::Element` | 加法单位元 |
| `one` | `fn one(&self) -> Self::Element` | 乘法单位元 |
| `add` | `fn add(&self, a, b) -> Self::Element` | 加法 |
| `sub` | `fn sub(&self, a, b) -> Self::Element` | 减法 |
| `neg` | `fn neg(&self, a) -> Self::Element` | 取反 |
| `mul` | `fn mul(&self, a, b) -> Self::Element` | 乘法 |
| `div` | `fn div(&self, a, b) -> Option<Self::Element>` | 除法。`b` 为零或除法不精确时返回 `None` |
| `inv` | `fn inv(&self, a) -> Option<Self::Element>` | 乘法逆元。`a` 为零时返回 `None` |
| `is_zero` | `fn is_zero(&self, a) -> bool` | 判断是否为加法单位元（默认 `*a == self.zero()`） |
| `is_one` | `fn is_one(&self, a) -> bool` | 判断是否为乘法单位元（默认 `*a == self.one()`） |
| `mul_assign` | `fn mul_assign(&self, a: &mut E, b: &E)` | 原地乘法 `*a *= b`。默认创建新元素；GMP 等高性能域可覆盖 |
| `sub_mul_assign` | `fn sub_mul_assign(&self, a: &mut E, b: &E, c: &E)` | 融合减乘 `*a -= b * c`。F4 行阶梯化中广泛使用 |
| `pow` | `fn pow(&self, a, n: u64) -> Self::Element` | 非负整数次幂。默认使用二进制快速幂；有限域覆盖为 `modpow` |
| `cast_u64` | `fn cast_u64(&self, n: u64) -> Self::Element` | 将 `u64` 转换为域元素。默认逐次加一 |

**示例**：

```rust
use ocas_domain::{Domain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(3);
let b = Integer::from(5);
assert_eq!(domain.add(&a, &b), Integer::from(8));
assert_eq!(domain.mul(&a, &b), Integer::from(15));
assert_eq!(domain.pow(&a, 3), Integer::from(27));
// 输出：所有断言通过
```

**参见**：[EuclideanDomain](#euclideandomain-trait)

---

## EuclideanDomain trait

**签名**：

```rust
pub trait EuclideanDomain: Domain {
    fn div_rem(&self, a: &Self::Element, b: &Self::Element)
        -> Option<(Self::Element, Self::Element)>;
    fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn extended_gcd(&self, a: &Self::Element, b: &Self::Element)
        -> (Self::Element, Self::Element, Self::Element);
}
```

**功能**：支持带余除法的欧几里得域。`Domain` 的扩展，额外提供 `div_rem`（带余除法）、`gcd`（最大公因子）和 `extended_gcd`（扩展欧几里得算法）。

**方法**：

### div_rem

**签名**：`fn div_rem(&self, a: &Self::Element, b: &Self::Element) -> Option<(Self::Element, Self::Element)>`

**功能**：带余除法，返回 `(商, 余数)`，满足 `a = 商 * b + 余数`，且 `余数 == 0` 或 `deg(余数) < deg(b)`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `&Self::Element` | 被除数 |
| `b` | `&Self::Element` | 除数 |

**返回值**：`Some((quotient, remainder))`，`b` 为零时返回 `None`。

**示例**：

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(17);
let b = Integer::from(5);
let (q, r) = domain.div_rem(&a, &b).unwrap();
assert_eq!(q, Integer::from(3));
assert_eq!(r, Integer::from(2));
// 输出：17 = 3 × 5 + 2
```

### gcd

**签名**：`fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element`

**功能**：计算最大公因子。默认实现使用欧几里得算法。

**返回值**：`gcd(a, b)`。对于域（如 `FiniteField`），两个非零元素的 GCD 退化为 1。

**示例**：

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(12);
let b = Integer::from(8);
let g = domain.gcd(&a, &b);
assert_eq!(g, Integer::from(4));
// 输出：gcd(12, 8) = 4
```

### extended_gcd

**签名**：`fn extended_gcd(&self, a: &Self::Element, b: &Self::Element) -> (Self::Element, Self::Element, Self::Element)`

**功能**：扩展欧几里得算法，返回 `(g, x, y)` 满足 `g = gcd(a, b) = a * x + b * y`。

**返回值**：三元组 `(g, x, y)`。

**示例**：

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(35);
let b = Integer::from(15);
let (g, x, y) = domain.extended_gcd(&a, &b);
// g = 5, x = 1, y = -2，满足 5 = 35×1 + 15×(-2)
assert_eq!(g, Integer::from(5));
```

**参见**：[Domain](#domain-trait)

---

## Integer / IntegerDomain

### IntegerDomain

**签名**：`pub struct IntegerDomain;`

**功能**：整数域 $\mathbb{Z}$，元素类型为 `Integer`。实现 `Domain` 和 `EuclideanDomain`。

**特征**：

- `div` 要求精确除法（余数为零），否则返回 `None`
- `inv` 仅对 $\pm 1$ 返回结果，其余返回 `None`
- `pow` 使用二进制快速幂

**示例**：

```rust
use ocas_domain::{Domain, EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(10);
let b = Integer::from(3);

// 精确除法：10/3 不精确
assert!(domain.div(&a, &b).is_none());

// 带余除法
let (q, r) = domain.div_rem(&a, &b).unwrap();
assert_eq!(q, Integer::from(3));
assert_eq!(r, Integer::from(1));

// 逆元：只有 ±1 有逆
assert!(domain.inv(&Integer::from(1)).is_some());
assert!(domain.inv(&Integer::from(5)).is_none());
// 输出：所有断言通过
```

### Integer

**签名**：`pub struct Integer(BigInt);`

**功能**：任意精度整数。默认构建使用 `num-bigint` 的 `BigInt`；启用 `gmp` feature 时使用 GMP 后端。

**Derive**：`Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash`

#### 构造与转换

##### Integer::new

**签名**：`pub fn new<T: Into<BigInt>>(value: T) -> Self`

**功能**：从机器整数或 `BigInt` 创建任意精度整数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `value` | `T: Into<BigInt>` | 整数值（`i32`、`i64`、`u64`、`BigInt` 等） |

**示例**：

```rust
use ocas_domain::Integer;

let a = Integer::new(42);
let b = Integer::new(100_i64);
assert_eq!(a.to_string(), "42");
// 输出：42
```

##### Integer::from

**签名**：`impl From<i64> for Integer` / `impl From<BigInt> for Integer`

**功能**：从 `i64` 或 `BigInt` 转换。

**示例**：

```rust
use ocas_domain::Integer;

let a = Integer::from(42);
assert_eq!(a.to_string(), "42");
// 输出：42
```

#### 访问器

##### Integer::inner

**签名**：`pub fn inner(&self) -> &BigInt`

**功能**：访问底层 `BigInt` 引用。

##### Integer::to_bigint

**签名**：`pub fn to_bigint(&self) -> BigInt`

**功能**：克隆为 `BigInt`（不论后端）。

##### Integer::to_i64

**签名**：`pub fn to_i64(&self) -> Option<i64>`

**功能**：尝试转换为 `i64`。溢出时返回 `None`。

#### 运算方法

##### Integer::pow_u32

**签名**：`pub fn pow_u32(&self, exp: u32) -> Self`

**功能**：计算 $n^{exp}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `exp` | `u32` | 非负指数 |

**示例**：

```rust
use ocas_domain::Integer;

let a = Integer::from(2);
assert_eq!(a.pow_u32(10).to_string(), "1024");
// 输出：1024
```

##### Integer::modpow

**签名**：`pub fn modpow(&self, exp: &Integer, modulus: &Integer) -> Integer`

**功能**：模幂运算 $self^{exp} \bmod modulus$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `exp` | `&Integer` | 指数（非负） |
| `modulus` | `&Integer` | 模数（正） |

**示例**：

```rust
use ocas_domain::Integer;

let base = Integer::from(3);
let exp = Integer::from(100);
let modulus = Integer::from(7);
let result = base.modpow(&exp, &modulus);
assert_eq!(result.to_string(), "4");
// 输出：3^100 mod 7 = 4
```

##### Integer::mod_floor

**签名**：`pub fn mod_floor(&self, modulus: &Integer) -> Integer`

**功能**：地板模运算，结果 $r$ 满足 $0 \leq r < |modulus|$。

**示例**：

```rust
use ocas_domain::Integer;

let a = Integer::from(-7);
let m = Integer::from(3);
assert_eq!(a.mod_floor(&m).to_string(), "2");
// 输出：-7 mod_floor 3 = 2
```

##### Integer::div_rem

**签名**：`pub fn div_rem(&self, other: &Integer) -> (Integer, Integer)`

**功能**：带余除法 `(商, 余数)`。注意：与 `EuclideanDomain::div_rem` 不同，此方法不要求 `other` 非零（行为由底层 `BigInt` 决定）。

##### Integer::is_even

**签名**：`pub fn is_even(&self) -> bool`

**功能**：判断是否为偶数。

##### Integer::is_negative

**签名**：`pub fn is_negative(&self) -> bool`

**功能**：判断是否为负数。

##### Integer::is_zero

**签名**：`pub fn is_zero(&self) -> bool`

**功能**：判断是否为零。

##### Integer::is_one

**签名**：`pub fn is_one(&self) -> bool`

**功能**：判断是否为一。

##### Integer::abs

**签名**：`pub fn abs(&self) -> Integer`

**功能**：绝对值。

##### Integer::sqrt

**签名**：`pub fn sqrt(&self) -> Integer`

**功能**：整数平方根（向下取整）。

**示例**：

```rust
use ocas_domain::Integer;

let a = Integer::from(10);
assert_eq!(a.sqrt().to_string(), "3");
// 输出：sqrt(10) = 3
```

#### 运算符

`Integer` 实现了以下标准运算符 trait（支持所有引用组合：owned × owned、owned × &、& × owned、& × &）：

| Trait | 运算 |
|---|---|
| `Add` | `+` |
| `Sub` | `-` |
| `Mul` | `*` |
| `Div` | `/`（整数除法，截断向零） |
| `Rem` | `%` |
| `Neg` | 一元 `-` |
| `Shr<u32>` / `ShrAssign<u32>` | 右移 |
| `AddAssign<&Integer>` / `SubAssign<&Integer>` / `MulAssign<&Integer>` / `DivAssign<&Integer>` | 复合赋值 |

**参见**：[Domain](#domain-trait)、[Rational](#rational--rationaldomain)、[数论函数](rust-ntheory.md)

---

## Rational / RationalDomain

### RationalDomain

**签名**：`pub struct RationalDomain;`

**功能**：有理数域 $\mathbb{Q}$，元素类型为 `Rational`。实现 `Domain` 和 `EuclideanDomain`。

**特征**：

- `div` 对所有非零除数精确返回结果
- `inv` 对所有非零元素返回结果
- `EuclideanDomain` 的 `div_rem` 对有理数域退化：余数恒为零

**示例**：

```rust
use ocas_domain::{Domain, Rational, RationalDomain};

let domain = RationalDomain;
let a = Rational::new(1, 2);
let b = Rational::new(1, 3);
let sum = domain.add(&a, &b);
assert_eq!(sum, Rational::new(5, 6));
// 输出：1/2 + 1/3 = 5/6
```

### Rational

**签名**：`pub struct Rational(BigRational);`

**功能**：任意精度有理数。默认构建使用 `num-rational` 的 `BigRational`；启用 `gmp` feature 时使用 GMP 后端。

**Derive**：`Debug, Clone, PartialEq, Eq, Hash`

#### 构造

##### Rational::new

**签名**：`pub fn new(numer: i64, denom: i64) -> Self`

**功能**：从分子和分母（`i64`）创建有理数。自动约简至最简分数。分母为零时行为由底层库决定。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `numer` | `i64` | 分子 |
| `denom` | `i64` | 分母（非零） |

**示例**：

```rust
use ocas_domain::Rational;

let a = Rational::new(3, 6);
assert_eq!(a.to_string(), "1/2");
// 输出：3/6 自动约简为 1/2
```

##### Rational::from_bigints

**签名**：`pub fn from_bigints(numer: BigInt, denom: BigInt) -> Self`

**功能**：从任意精度整数分子和分母创建有理数。

##### Rational::from_integer

**签名**：`pub fn from_integer(n: Integer) -> Self`

**功能**：从整数创建有理数（分母 = 1）。

#### 访问器

##### Rational::inner

**签名**：`pub fn inner(&self) -> &BigRational`

**功能**：访问底层 `BigRational` 引用。

##### Rational::numer

**签名**：`pub fn numer(&self) -> Integer`

**功能**：返回分子（作为 `Integer`）。

##### Rational::denom

**签名**：`pub fn denom(&self) -> Integer`

**功能**：返回分母（作为 `Integer`，始终为正）。

**示例**：

```rust
use ocas_domain::Rational;

let r = Rational::new(3, 4);
assert_eq!(r.numer().to_string(), "3");
assert_eq!(r.denom().to_string(), "4");
// 输出：分子 = 3，分母 = 4
```

#### 运算符

`Rational` 实现了以下标准运算符 trait：

| Trait | 运算 |
|---|---|
| `Add` | `+` |
| `Sub` | `-` |
| `Mul` | `*` |
| `Div` | `/` |
| `Neg` | 一元 `-` |
| `AddAssign` / `SubAssign` / `MulAssign` / `DivAssign` | 复合赋值 |

**参见**：[Domain](#domain-trait)、[Integer](#integer--integerdomain)

---

## FiniteField / FiniteFieldElement

### FiniteField

**签名**：`pub struct FiniteField { /* prime: BigInt, ... */ }`

**功能**：素数有限域 $\mathbb{Z}/p\mathbb{Z}$。算术使用任意精度整数，支持大素数。实现 `Domain` 和 `EuclideanDomain`。

#### 构造

##### FiniteField::new

**签名**：`pub fn new(prime: BigInt) -> Self`

**功能**：以给定素数为模创建有限域。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `prime` | `BigInt` | 模数（要求 $\geq 2$，不验证素性） |

**错误**：debug 模式下若 `prime < 2` 会 panic。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};

let f = FiniteField::new(BigInt::from(7));
let a = f.element(3);
let b = f.element(5);
assert_eq!(f.add(&a, &b), f.element(1));   // 3 + 5 = 1 (mod 7)
assert_eq!(f.mul(&a, &b), f.element(1));   // 3 × 5 = 1 (mod 7)
// 输出：所有断言通过
```

#### 元素构造

##### FiniteField::element

**签名**：`pub fn element(&self, value: impl Into<BigInt>) -> FiniteFieldElement`

**功能**：从任意整数创建域元素。值自动规约到 $[0, p-1]$ 范围。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `value` | `impl Into<BigInt>` | 整数值（可超出 $[0, p-1]$，自动取模） |

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};

let f = FiniteField::new(BigInt::from(7));
let a = f.element(10);
assert_eq!(a.value().to_string(), "3");
// 输出：10 mod 7 = 3
```

##### FiniteField::from_i64

**签名**：`pub fn from_i64(&self, val: i64) -> FiniteFieldElement`

**功能**：从 `i64` 创建域元素（规约 mod p）。

#### 访问器

##### FiniteField::prime

**签名**：`pub fn prime(&self) -> &BigInt`

**功能**：返回域的模数。

##### FiniteField::prime_u64

**签名**：`pub fn prime_u64(&self) -> u64`

**功能**：将模数作为 `u64` 返回。

**错误**：若素数不适用 `u64` 表示则 panic。

##### FiniteField::to_i64

**签名**：`pub fn to_i64(&self, a: &FiniteFieldElement) -> i64`

**功能**：将域元素转换为 `i64`（范围 $[0, p)$）。

**错误**：若素数不适用 `u64` 表示则 panic。

#### Domain 实现

- `div` 通过 Fermat 小定理求逆：$a^{-1} \equiv a^{p-2} \pmod{p}$
- `inv` 同上，零元素返回 `None`
- `pow` 使用 `modpow` 优化（远快于默认二进制幂）

#### EuclideanDomain 实现

- `div_rem`：域中除法精确，余数恒为零
- `gcd`：域中 GCD 退化——两者均为零时返回 0，否则返回 1

**参见**：[Domain](#domain-trait)、[有限域数学基础](../math/finite-fields.md)

---

### FiniteFieldElement

**签名**：`pub struct FiniteFieldElement { value: BigInt }`

**功能**：素数有限域的元素。值始终在 $[0, p-1]$ 范围内。

**Derive**：`Debug, Clone, PartialEq, Eq, Hash`

#### FiniteFieldElement::value

**签名**：`pub fn value(&self) -> &BigInt`

**功能**：返回 $[0, p-1]$ 范围内的标准代表元。

**示例**：

```rust
use num_bigint::BigInt;
use ocas_domain::FiniteField;

let f = FiniteField::new(BigInt::from(7));
let a = f.element(-3);
assert_eq!(a.value().to_string(), "4");
// 输出：-3 mod 7 = 4
```

**参见**：[FiniteField](#finitefield--finitefieldelement)

---

## RealBall / RealBallDomain

### RealBallDomain

**签名**：`pub struct RealBallDomain;`

**功能**：实数球（区间）域。元素类型为 `RealBall`。仅实现 `Domain`（不实现 `EuclideanDomain`，因为实数球不支持精确带余除法）。

**注意**：默认构建使用轻量 `f64` 球，适合模板和演示。启用 `mpfr` feature 后使用 `rug::Float` 和有向舍入产生严格区间。

### RealBall

**签名**：`pub struct RealBall { mid: f64, rad: f64 }`（默认构建）

**功能**：实数球：中点 ± 半径。真实值保证包含在 $[mid - rad, mid + rad]$ 中。

**Derive**：`Debug, Clone, PartialEq, Eq`；默认构建额外实现 `Copy`。

#### 构造

##### RealBall::new

**签名**（默认）：`pub fn new(mid: f64, rad: f64) -> Self`

**功能**：从中点和半径创建球。半径被钳制为非负。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `mid` | `f64`（或 `rug::Float`） | 中点 |
| `rad` | `f64`（或 `rug::Float`） | 半径（$\geq 0$） |

##### RealBall::from_f64

**签名**：`pub fn from_f64(value: f64) -> Self`

**功能**：从 `f64` 值创建零半径球（精确值）。

**示例**：

```rust
use ocas_domain::RealBall;

let ball = RealBall::from_f64(3.14);
assert_eq!(ball.mid(), 3.14);
assert_eq!(ball.rad(), 0.0);
// 输出：精确球
```

#### 访问器

##### RealBall::mid

**签名**：`pub fn mid(&self) -> f64`（默认） / `pub fn mid(&self) -> &rug::Float`（mpfr）

**功能**：返回中点。

##### RealBall::rad

**签名**：`pub fn rad(&self) -> f64`（默认） / `pub fn rad(&self) -> &rug::Float`（mpfr）

**功能**：返回半径。

##### RealBall::lower

**签名**：`pub fn lower(&self) -> f64`（默认） / `pub fn lower(&self) -> rug::Float`（mpfr）

**功能**：返回保守下界 $mid - rad$（mpfr 版本使用向下舍入）。

##### RealBall::upper

**签名**：`pub fn upper(&self) -> f64`（默认） / `pub fn upper(&self) -> rug::Float`（mpfr）

**功能**：返回保守上界 $mid + rad$（mpfr 版本使用向上舍入）。

##### RealBall::precision

**签名**：`pub fn precision(&self) -> u32`（仅 `mpfr` feature）

**功能**：返回 MPFR 后端的精度（位数）。

#### Domain 实现

- `add`：$(a \pm r_a) + (b \pm r_b) = (a+b) \pm (r_a + r_b)$（mpfr 版本额外加舍入误差）
- `sub`：类似加法，半径相加
- `mul`：四角法——计算四个极端乘积 $a_lo \cdot b_lo$、$a_lo \cdot b_hi$ 等，取最小/最大值得到新区间
- `div`：$a / b = a \cdot b^{-1}$，当球包含零时返回 `None`
- `inv`：$1/(mid \pm rad)$，当球包含零时返回 `None`

**示例**：

```rust
use ocas_domain::{Domain, RealBall, RealBallDomain};

let domain = RealBallDomain;
let a = RealBall::from_f64(2.0);
let b = RealBall::from_f64(3.0);
let prod = domain.mul(&a, &b);
assert!(prod.lower() <= 6.0 && 6.0 <= prod.upper());
// 输出：乘积球包含真值 6.0
```

**参见**：[Domain](#domain-trait)

---

## Complex / ComplexDomain

### ComplexDomain

**签名**：`pub struct ComplexDomain<D: Domain> { base: D, ... }`

**功能**：在任意基域 $D$ 上构造的复数域。元素类型为 `Complex<D>`。仅实现 `Domain`。

#### 构造

##### ComplexDomain::new

**签名**：`pub fn new(base: D) -> Self`

**功能**：在基域 `base` 上创建复数域。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `base` | `D: Domain` | 基域（如 `IntegerDomain`、`RationalDomain`） |

##### ComplexDomain::base

**签名**：`pub fn base(&self) -> &D`

**功能**：返回基域引用。

##### ComplexDomain::real_element

**签名**：`pub fn real_element(&self, re: D::Element) -> Complex<D>`

**功能**：创建纯实数元素。

##### ComplexDomain::imag_element

**签名**：`pub fn imag_element(&self, im: D::Element) -> Complex<D>`

**功能**：创建纯虚数元素。

#### Domain 实现

- `mul`：$(a+bi)(c+di) = (ac-bd) + (ad+bc)i$
- `div`：$(a+bi)/(c+di) = \frac{(ac+bd)+(bc-ad)i}{c^2+d^2}$，分母为零时返回 `None`
- `inv`：通过 `div(one, a)` 实现

**示例**：

```rust
use ocas_domain::{Complex, ComplexDomain, Domain, Integer, IntegerDomain};

let domain = ComplexDomain::new(IntegerDomain);
let a = Complex::new(Integer::from(1), Integer::from(2));
let b = Complex::new(Integer::from(3), Integer::from(4));
let sum = domain.add(&a, &b);
assert_eq!(*sum.re(), Integer::from(4));
assert_eq!(*sum.im(), Integer::from(6));
// 输出：(1+2i) + (3+4i) = 4+6i
```

**参见**：[Domain](#domain-trait)

---

### Complex

**签名**：`pub struct Complex<D: Domain> { inner: NumComplex<D::Element> }`

**功能**：复数，实部和虚部属于基域 $D$。

**Derive**：`Debug, Clone, PartialEq, Eq, Hash`

#### 构造与访问

##### Complex::new

**签名**：`pub fn new(real: D::Element, imag: D::Element) -> Self`

**功能**：从实部和虚部创建复数。

##### Complex::re

**签名**：`pub fn re(&self) -> &D::Element`

**功能**：返回实部引用。

##### Complex::im

**签名**：`pub fn im(&self) -> &D::Element`

**功能**：返回虚部引用。

##### Complex::inner

**签名**：`pub fn inner(&self) -> &NumComplex<D::Element>`

**功能**：返回底层 `num_complex::Complex` 引用。

**参见**：[ComplexDomain](#complexdomain)

---

## DoubleF64 / DoubleF64Domain

### DoubleF64Domain

**签名**：`pub struct DoubleF64Domain;`

**功能**：双精度浮点域。元素类型为 `DoubleF64`。仅实现 `Domain`。

### DoubleF64

**签名**：

```rust
pub struct DoubleF64 {
    pub hi: f64,  // 高阶分量（主值）
    pub lo: f64,  // 低阶分量（误差项）
}
```

**功能**：双精度浮点数，表示为 $hi + lo$，满足 $|lo| \leq 0.5 \cdot \text{ulp}(hi)$。提供约 31 位十进制有效数字（~84 位二进制），约为单个 `f64` 精度的两倍。

**Derive**：`Debug, Clone, Copy, PartialEq`；实现 `Eq`。

算术基于 Dekker 和 Knuth 的无误差变换算法（TwoSum、TwoProd），显著快于 MPFR 等任意精度替代方案。

#### 常量

| 常量 | 值 | 说明 |
|---|---|---|
| `DoubleF64::ZERO` | `{ hi: 0.0, lo: 0.0 }` | 零 |
| `DoubleF64::ONE` | `{ hi: 1.0, lo: 0.0 }` | 一 |

#### 构造与转换

##### DoubleF64::new

**签名**：`pub fn new(hi: f64, lo: f64) -> Self`

**功能**：从高阶和低阶分量创建。调用者需确保 $|lo| \leq 0.5 \cdot \text{ulp}(hi)$ 以保证正确性。

##### DoubleF64::from_f64

**签名**：`pub fn from_f64(x: f64) -> Self`

**功能**：从单个 `f64` 创建（$lo = 0$）。

##### DoubleF64::to_f64

**签名**：`pub fn to_f64(self) -> f64`

**功能**：提取高阶分量。

##### From\<f64\>

**签名**：`impl From<f64> for DoubleF64`

**功能**：等价于 `from_f64`。

#### 无误差变换

##### DoubleF64::quick_two_sum

**签名**：`pub fn quick_two_sum(a: f64, b: f64) -> Self`

**功能**：当 $|a| \geq |b|$ 时的无误差求和。比 `two_sum` 快但有前提条件。

##### DoubleF64::two_sum

**签名**：`pub fn two_sum(a: f64, b: f64) -> Self`

**功能**：Dekker TwoSum——无误差求和，舍入误差精确捕获到 `lo` 分量。

**示例**：

```rust
use ocas_domain::DoubleF64;

let s = DoubleF64::two_sum(1.0, f64::EPSILON);
assert_eq!(s.hi, 1.0 + f64::EPSILON);
// s.lo 捕获了舍入误差
```

#### 查询

| 方法 | 签名 | 说明 |
|---|---|---|
| `abs` | `pub fn abs(self) -> Self` | 绝对值 |
| `is_nan` | `pub fn is_nan(self) -> bool` | 是否为 NaN |
| `is_infinite` | `pub fn is_infinite(self) -> bool` | 是否为无穷 |
| `is_finite` | `pub fn is_finite(self) -> bool` | 是否有限 |

#### 算术运算

##### DoubleF64::add / sub / mul / div

**签名**：

```rust
pub fn add(self, other: Self) -> Self
pub fn sub(self, other: Self) -> Self
pub fn mul(self, other: Self) -> Self
pub fn div(self, other: Self) -> Self
```

**功能**：双精度加减乘除。内部使用 TwoSum/TwoProd 捕获舍入误差。

##### DoubleF64::powi

**签名**：`pub fn powi(self, n: i64) -> Self`

**功能**：整数次幂，使用二进制快速幂。支持负指数。

**示例**：

```rust
use ocas_domain::DoubleF64;

let x = DoubleF64::from_f64(3.0);
assert_eq!(x.powi(3).hi, 27.0);
assert_eq!(x.powi(-1).hi, 1.0 / 3.0);
// 输出：3^3 = 27, 3^(-1) = 0.333...
```

#### 超越函数

| 方法 | 签名 | 说明 |
|---|---|---|
| `sqrt` | `pub fn sqrt(self) -> Self` | Newton 迭代平方根。负数返回 NaN |
| `exp` | `pub fn exp(self) -> Self` | 指数函数，Taylor 级数 + 参数规约（$\ln 2$ 缩减） |
| `ln` | `pub fn ln(self) -> Self` | 自然对数，Newton 迭代。非正数返回 NaN |
| `sin` | `pub fn sin(self) -> Self` | 正弦，Taylor 级数 + $[-\pi, \pi]$ 规约 |
| `cos` | `pub fn cos(self) -> Self` | 余弦，$\cos(x) = \sin(\pi/2 - x)$ |
| `tan` | `pub fn tan(self) -> Self` | 正切，$\sin(x)/\cos(x)$ |

**示例**：

```rust
use ocas_domain::DoubleF64;

let x = DoubleF64::from_f64(1.0);
let e = x.exp();
// e ≈ 2.718281828...
assert!((e.hi - std::f64::consts::E).abs() < 1e-15);

let pi = DoubleF64::from_f64(std::f64::consts::PI);
let s = pi.sin();
assert!(s.hi.abs() < 1e-15); // sin(π) ≈ 0
// 输出：exp(1) ≈ e, sin(π) ≈ 0
```

#### 运算符

`DoubleF64` 实现了 `Add`、`Sub`、`Mul`、`Div`、`Neg`、`AddAssign`、`SubAssign`、`MulAssign`、`DivAssign`、`PartialOrd`、`Zero`（`num_traits`）。

#### Display

当 `lo == 0` 时显示 `hi` 值；否则以 31 位科学计数法显示 $hi + lo$ 的和。

**参见**：[Domain](#domain-trait)、[RealBall](#realball--realballdomain)

---

## AlgebraicExtension / AlgebraicNumberField / AlgebraicElement

### AlgebraicExtension

**签名**：`pub struct AlgebraicExtension<D: Domain> { base: D, min_poly: Vec<D::Element> }`

**功能**：代数扩张 $D[\alpha]/(m(\alpha))$，其中 $m$ 是首一多项式。当基域 $D$ 是域且 $m$ 不可约时，商环是域：

- `AlgebraicExtension<RationalDomain>` = 代数数域 $\mathbb{Q}(\alpha)$
- `AlgebraicExtension<FiniteField>` = Galois 域 $\mathrm{GF}(p^d)$

元素是剩余类，由次数小于 $\deg(m)$ 的唯一多项式代表表示。逆元使用基域上的扩展欧几里得算法。

**⚠️ 注意**：不检查极小多项式的不可约性。在可约模下，环有零因子，`Domain::inv` 对非单位返回 `None`。

**实现**：`Domain`、`EuclideanDomain`

#### 构造

##### AlgebraicExtension::new

**签名**：`pub fn new(base: D, min_poly: Vec<D::Element>) -> Self`

**功能**：从基域和首一极小多项式创建代数扩张。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `base` | `D: Domain` | 基域 |
| `min_poly` | `Vec<D::Element>` | 极小多项式，升序排列，首一（最高次系数为 1），次数 $\geq 1$ |

**错误**：debug 模式下检查：(1) 次数 $\geq 1$（即 `min_poly.len() >= 2`）；(2) 首一性。

**示例**：

```rust
use ocas_domain::{AlgebraicExtension, Domain, Rational, RationalDomain};

// ℚ(√2)：极小多项式 α² − 2
let two = Rational::new(2, 1);
let neg_two = RationalDomain.neg(&two);
let field = AlgebraicExtension::new(
    RationalDomain,
    vec![neg_two, Rational::new(0, 1), Rational::new(1, 1)],
);
let sqrt2 = field.alpha();
// √2 · √2 = 2
assert_eq!(field.mul(&sqrt2, &sqrt2), field.from_base(two));
// 输出：α² = 2
```

#### 访问器

##### AlgebraicExtension::base_domain

**签名**：`pub fn base_domain(&self) -> &D`

**功能**：返回基域引用。

##### AlgebraicExtension::min_poly

**签名**：`pub fn min_poly(&self) -> &[D::Element]`

**功能**：返回极小多项式系数（升序，首一）。

##### AlgebraicExtension::extension_degree

**签名**：`pub fn extension_degree(&self) -> usize`

**功能**：返回扩张次数 $\deg(m)$。

#### 元素构造

##### AlgebraicExtension::from_base

**签名**：`pub fn from_base(&self, c: D::Element) -> AlgebraicElement<D::Element>`

**功能**：将基域常量嵌入扩张域。

##### AlgebraicExtension::alpha

**签名**：`pub fn alpha(&self) -> AlgebraicElement<D::Element>`

**功能**：返回扩张的生成元 $\alpha$。

##### AlgebraicExtension::element

**签名**：`pub fn element(&self, coeffs: Vec<D::Element>) -> AlgebraicElement<D::Element>`

**功能**：从系数（升序）创建元素，自动规约模极小多项式。

#### Domain 实现

- `mul`：先多项式乘法，再规约模 $m(\alpha)$
- `inv`：使用扩展欧几里得算法——若 $\gcd(a, m) = 1$（常数），则 $a^{-1} = s \bmod m$
- `div`：$a/b = a \cdot b^{-1}$
- `is_zero`：系数向量为空

#### EuclideanDomain 实现

- `div_rem`：域上除法精确，余数恒为零
- `gcd`：域中退化——两者均零返回 0，否则返回 1

**示例**（Galois 域 $\mathrm{GF}(3^2)$）：

```rust
use ocas_domain::{AlgebraicExtension, Domain};
use ocas_domain::FiniteField;
use num_bigint::BigInt;

let base = FiniteField::new(BigInt::from(3));
let field = AlgebraicExtension::new(
    base.clone(),
    vec![base.element(1), base.element(0), base.element(1)], // α² + 1
);
let alpha = field.alpha();
// α² = −1 = 2 (mod 3)
assert_eq!(field.mul(&alpha, &alpha), field.from_base(base.element(2)));
// 乘法群阶为 8：(1+α)⁸ = 1
let a = field.add(&field.one(), &alpha);
assert_eq!(field.pow(&a, 8), field.one());
// 输出：α² = 2, (1+α)⁸ = 1
```

**参见**：[Domain](#domain-trait)、[代数数域数学基础](../math/algebraic-number-fields.md)

---

### AlgebraicNumberField

**签名**：`pub type AlgebraicNumberField = AlgebraicExtension<RationalDomain>;`

**功能**：$\mathbb{Q}(\alpha)$ 的类型别名。

---

### AlgebraicElement

**签名**：`pub struct AlgebraicElement<E> { coeffs: Vec<E> }`

**功能**：代数扩张的元素——$\alpha$ 的多项式剩余类，次数小于扩张次数。系数升序存储，尾随零被裁剪。零元素的系数向量为空。

**Derive**：`Debug, Clone, PartialEq, Eq, Hash`

#### AlgebraicElement::coeffs

**签名**：`pub fn coeffs(&self) -> &[E]`

**功能**：返回系数切片（升序，尾随零已裁剪）。

#### Display

显示格式：`0`（零元素）、`c`（纯常数）、`(c1)·α + c0`、`(c2)·α^2 + (c1)·α + c0` 等（常数项不加括号，各次项之间以 ` + ` 连接）。

**示例**：

```rust
use ocas_domain::{AlgebraicExtension, Domain, Rational, RationalDomain};

// ℚ(i)：极小多项式 α² + 1
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(1, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
let i = field.alpha();
let one_plus_i = field.add(&field.one(), &i);
// (1+i)⁻¹ = (1-i)/2
let inv = field.inv(&one_plus_i).unwrap();
assert_eq!(inv.coeffs(), &[Rational::new(1, 2), Rational::new(-1, 2)]);
// 输出：1/(1+i) = (1/2) - (1/2)·α
```

**参见**：[AlgebraicExtension](#algebraicextension--algebraicnumberfield--algebraicelement)

---

## Assumptions 系统

假设系统用于声明符号变量的性质（如"x 是实数"、"n 是正整数"）。求解器和化简器通过假设来选择算法和验证解。

**模块路径**：`ocas_domain::assumptions`

---

### Assumption

**签名**：

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Assumption {
    Real, Complex, Integer, Rational,
    Positive, Negative, NonNegative, NonPositive, NonZero, Finite,
    Even, Odd, Prime,
}
```

**功能**：可对符号变量声明的单一谓词。假设相互独立——一个变量可同时携带多个假设（如 `Positive | Integer`）。

#### 变体说明

| 变体 | 含义 | Display |
|---|---|---|
| `Real` | 实数 | `"real"` |
| `Complex` | 复数 | `"complex"` |
| `Integer` | 整数 | `"integer"` |
| `Rational` | 有理数 | `"rational"` |
| `Positive` | 严格正（$> 0$） | `"positive"` |
| `Negative` | 严格负（$< 0$） | `"negative"` |
| `NonNegative` | 非负（$\geq 0$） | `"non-negative"` |
| `NonPositive` | 非正（$\leq 0$） | `"non-positive"` |
| `NonZero` | 非零 | `"non-zero"` |
| `Finite` | 有限（非 $\pm\infty$） | `"finite"` |
| `Even` | 偶数 | `"even"` |
| `Odd` | 奇数 | `"odd"` |
| `Prime` | 素数 | `"prime"` |

#### Assumption::implied

**签名**：`pub fn implied(&self) -> &'static [Assumption]`

**功能**：返回此假设逻辑蕴含的其他假设列表。插入假设时，其蕴含会自动传播。

**蕴含关系表**：

| 假设 | 蕴含 |
|---|---|
| `Positive` | `NonNegative`, `NonZero`, `Real` |
| `Negative` | `NonPositive`, `NonZero`, `Real` |
| `NonNegative` | `Real` |
| `NonPositive` | `Real` |
| `Integer` | `Rational`, `Real` |
| `Rational` | `Real` |
| `Complex` | `Real` |
| `Even` | `Integer` |
| `Odd` | `Integer` |
| `Prime` | `Integer`, `Positive` |
| `Real`, `NonZero`, `Finite` | （无额外蕴含） |

#### Assumption::conflicts

**签名**：`pub fn conflicts(&self) -> &'static [Assumption]`

**功能**：返回与之矛盾的假设列表。同时包含一个假设及其冲突项的集合是不一致的。

**冲突关系表**：

| 假设 | 与之冲突 |
|---|---|
| `Positive` | `Negative`, `NonPositive` |
| `Negative` | `Positive`, `NonNegative` |
| `NonNegative` | `Negative` |
| `NonPositive` | `Positive` |
| `Even` | `Odd` |
| `Odd` | `Even` |

其他假设无冲突项。

#### BitOr 支持

`Assumption | Assumption` 返回 `Assumptions`：

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let a = Assumption::Positive | Assumption::Integer;
assert!(a.implies(Assumption::Real));
// 输出：Positive | Integer 蕴含 Real
```

**参见**：[Assumptions](#assumptions)、[SymbolAssumptions](#symbolassumptions)

---

### Assumptions

**签名**：`pub struct Assumptions { inner: Vec<Assumption> }`

**功能**：关于一个符号变量的假设集合。内部以排序、去重向量存储以优化小规模效率。操作在逻辑蕴含下封闭——插入 `Positive` 会使 `NonNegative` 和 `Real` 可用。

**Derive**：`Debug, Clone, PartialEq, Eq, Default`

#### 构造

##### Assumptions::new

**签名**：`pub fn new() -> Self`

**功能**：创建空假设集。

##### Assumptions::single

**签名**：`pub fn single(a: Assumption) -> Self`

**功能**：创建包含单个假设（及其蕴含）的集合。

#### 查询

##### Assumptions::len

**签名**：`pub fn len(&self) -> usize`

**功能**：返回显式存储的假设数量。

##### Assumptions::is_empty

**签名**：`pub fn is_empty(&self) -> bool`

**功能**：是否为空集。

##### Assumptions::contains

**签名**：`pub fn contains(&self, a: Assumption) -> bool`

**功能**：检查假设是否被此集合蕴含（直接成员或逻辑蕴含）。

##### Assumptions::implies

**签名**：`pub fn implies(&self, other: Assumption) -> bool`

**功能**：检查此集合是否逻辑蕴含 `other`。

**示例**：

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let mut a = Assumptions::new();
a.insert(Assumption::Positive);
a.insert(Assumption::Integer);
assert!(a.contains(Assumption::Real));      // Positive 蕴含 Real
assert!(a.implies(Assumption::NonZero));    // Positive 蕴含 NonZero
assert!(!a.implies(Assumption::Even));      // 不蕴含 Even
// 输出：所有断言通过
```

#### 修改

##### Assumptions::insert

**签名**：`pub fn insert(&mut self, a: Assumption) -> bool`

**功能**：插入假设及其所有逻辑蕴含。返回 `false` 若插入导致矛盾（不一致假设仍被存储，应检查 `is_consistent`）。

##### Assumptions::remove

**签名**：`pub fn remove(&mut self, a: Assumption)`

**功能**：移除假设。不移除由被移除假设所蕴含的其他假设（它们可能仍由其他存储的假设蕴含）。

##### Assumptions::is_consistent

**签名**：`pub fn is_consistent(&self) -> bool`

**功能**：检查集合是否一致（无矛盾假设）。

##### Assumptions::iter

**签名**：`pub fn iter(&self) -> impl Iterator<Item = Assumption> + '_`

**功能**：遍历存储的假设。

#### 运算符

- `Assumption | Assumption` → `Assumptions`
- `Assumptions | Assumption` → `Assumptions`
- `Assumptions | Assumptions` → `Assumptions`（并集）
- `FromIterator<Assumption>` → `Assumptions`

**示例**：

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let a: Assumptions = [Assumption::Positive, Assumption::Integer].into_iter().collect();
assert!(a.is_consistent());

let b = Assumption::Positive | Assumption::Negative;
assert!(!b.is_consistent());
// 输出：Positive + Integer 一致；Positive + Negative 不一致
```

**参见**：[Assumption](#assumption)、[SymbolAssumptions](#symbolassumptions)

---

### SymbolAssumptions

**签名**：`pub struct SymbolAssumptions { entries: Vec<(String, Assumptions)> }`

**功能**：符号名到假设的映射。求解器和化简器用此确定合法变换。例如 $\sqrt{x^2} \to x$ 仅在 $x$ 被假设为 `NonNegative` 时有效。

**Derive**：`Debug, Clone, PartialEq, Eq, Default`

#### 构造

##### SymbolAssumptions::new

**签名**：`pub fn new() -> Self`

**功能**：创建空映射。

#### 查询与修改

##### SymbolAssumptions::set

**签名**：`pub fn set(&mut self, symbol: &str, assumptions: Assumptions)`

**功能**：为符号设置假设（替换已有条目）。

##### SymbolAssumptions::get

**签名**：`pub fn get(&self, symbol: &str) -> Option<&Assumptions>`

**功能**：获取符号的假设。

##### SymbolAssumptions::remove

**签名**：`pub fn remove(&mut self, symbol: &str)`

**功能**：移除符号的假设。

##### SymbolAssumptions::check

**签名**：`pub fn check(&self, symbol: &str, assumption: Assumption) -> bool`

**功能**：检查符号是否满足特定假设。

##### SymbolAssumptions::len

**签名**：`pub fn len(&self) -> usize`

**功能**：返回有假设的符号数量。

##### SymbolAssumptions::is_empty

**签名**：`pub fn is_empty(&self) -> bool`

**功能**：是否无符号有假设。

##### SymbolAssumptions::iter

**签名**：`pub fn iter(&self) -> impl Iterator<Item = &(String, Assumptions)>`

**功能**：遍历所有 `(符号名, 假设集)` 对。

**示例**：

```rust
use ocas_domain::assumptions::{Assumption, Assumptions, SymbolAssumptions};

let mut sa = SymbolAssumptions::new();
sa.set("x", Assumptions::single(Assumption::Positive));
sa.set("n", Assumption::Integer | Assumption::Positive);

assert!(sa.check("x", Assumption::Real));      // Positive 蕴含 Real
assert!(sa.check("n", Assumption::NonZero));   // Positive 蕴含 NonZero
assert!(!sa.check("x", Assumption::Integer));  // 不蕴含 Integer
// 输出：所有断言通过
```

**参见**：[Assumptions](#assumptions)、[Assumption](#assumption)
