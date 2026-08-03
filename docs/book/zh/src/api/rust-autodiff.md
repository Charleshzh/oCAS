# 自动微分

oCAS 通过超对偶数（hyper-dual numbers）实现前向模式自动微分。一个 `HyperDual<T>` 携带一个标量值和一组导数分量，分量布局由共享的 `DualShape`（通过 `Arc` 共享）预先计算，并在算术运算中使用预构建的乘法表。

> **源文件**：`ocas-domain/src/dual.rs`

## DualCoeff trait

```rust
pub trait DualCoeff:
    Clone
    + PartialEq
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
    + std::ops::Neg<Output = Self>
    + std::ops::AddAssign
    + std::ops::MulAssign
{
    fn zero() -> Self;
    fn one() -> Self;
}
```

**功能**：定义超对偶数所需的系数类型约束。要求完整的四则运算（含乘法逆元）以及加法/乘法单位元。

**参数**：无（trait bound）。

**内置实现**：

| 类型 | 说明 |
|---|---|
| `Rational` | 任意精度有理数，`zero()` = `0/1`，`one()` = `1/1` |

**限制**：当前仅 `Rational` 实现了 `DualCoeff`。不支持 `f64` 或其他浮点类型作为系数。

## DualShape

```rust
#[derive(Debug, Clone)]
pub struct DualShape { /* private fields */ }
```

**功能**：描述超对偶数的分量布局。每个分量由一个多指标（multi-index）标识——长度为变量数的非负整数向量，记录各变量的微分阶数。布局必须满足**祖先封闭性**：若多指标 $\mathbf{m}$ 存在，则所有分量逐项不大于 $\mathbf{m}$ 的多指标也必须存在。

分量 0 始终是全零多指标（标量值分量）。

### DualShape::new

```rust
pub fn new(components: Vec<Vec<usize>>) -> Option<Self>
```

**功能**：从祖先封闭的多指标列表构建布局。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `components` | `Vec<Vec<usize>>` | 多指标列表，每个 `Vec<usize>` 长度相同（不足部分自动补零） |

**返回值**：`Some(DualShape)` 若布局合法，否则 `None`。

**错误条件**：
- 列表为空 → `None`
- 不含全零多指标 → `None`
- 不满足祖先封闭性 → `None`

**示例**：

```rust
use ocas_domain::dual::DualShape;

// 一阶二变量：值 + ∂/∂x₀ + ∂/∂x₁
let shape = DualShape::new(vec![
    vec![0, 0],
    vec![1, 0],
    vec![0, 1],
]).unwrap();
assert_eq!(shape.n_components(), 3);
assert_eq!(shape.n_vars(), 2);

// 非法：缺少 [1]（不祖先封闭）
assert!(DualShape::new(vec![vec![0], vec![2]]).is_none());
```

### DualShape::n_vars

```rust
pub fn n_vars(&self) -> usize
```

**功能**：返回微分变量数（多指标的长度）。

**返回值**：`usize`

### DualShape::n_components

```rust
pub fn n_components(&self) -> usize
```

**功能**：返回分量总数（含标量值分量）。

**返回值**：`usize`

### DualShape::components

```rust
pub fn components(&self) -> &[Vec<usize>]
```

**功能**：返回分量多指标切片。索引 0 为标量值分量。

**返回值**：`&[Vec<usize>]`

### DualShape::index_of

```rust
pub fn index_of(&self, multi_index: &[usize]) -> Option<usize>
```

**功能**：查找给定多指标的分量索引。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `multi_index` | `&[usize]` | 待查找的多指标 |

**返回值**：`Some(index)` 若存在，否则 `None`。

### DualShape::mult_table

```rust
pub fn mult_table(&self) -> &[(usize, usize, usize)]
```

**功能**：返回乘法表（用于调试/测试）。每个三元组 `(a, b, c)` 表示分量 `a` 与分量 `b` 的乘积贡献到分量 `c`。不包含标量分量（索引 0）参与的对。

**返回值**：`&[(usize, usize, usize)]`

## new_first_order

```rust
pub fn new_first_order<T: DualCoeff>(nvars: usize) -> Arc<DualShape>
```

**功能**：构建一阶形状，追踪每个变量的偏导数 $\frac{\partial}{\partial x_i}$。不追踪高阶或混合偏导。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `nvars` | `usize` | 微分变量数 |

**返回值**：`Arc<DualShape>` — 包含 $n+1$ 个分量：值分量 $[0,\dots,0]$ 加每个变量的 $[0,\dots,\underset{i}{1},\dots,0]$。

**示例**：

```rust
use ocas_domain::dual::new_first_order;
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(3);
assert_eq!(shape.n_components(), 4); // 值 + 3 个偏导
assert_eq!(shape.n_vars(), 3);
```

## HyperDual\<T\>

```rust
#[derive(Debug, Clone)]
pub struct HyperDual<T: DualCoeff> { /* private fields */ }
```

**功能**：超对偶数类型，携带标量值和按 `DualShape` 布局的导数分量。所有算术运算自动传播导数。

**类型约束**：`T: DualCoeff`

### HyperDual::variable

```rust
pub fn variable(shape: &Arc<DualShape>, i: usize, c: T) -> Self
```

**功能**：构造第 $i$ 个独立变量，取值为 $c$。导数分量 $[0,\dots,\underset{i}{1},\dots,0]$ 设为 1，其余为 0。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `&Arc<DualShape>` | 分量布局 |
| `i` | `usize` | 变量编号（0-indexed），须 < `n_vars()`（越界时静默忽略，不 panic） |
| `c` | `T` | 该变量在求值点的值 |

**返回值**：`HyperDual<T>`

**示例**：

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(2);
let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));
assert_eq!(x.value(), &Rational::new(3, 1));
assert_eq!(x.deriv(0), Some(&Rational::new(1, 1)));
assert_eq!(x.deriv(1), Some(&Rational::new(0, 1)));
```

### HyperDual::constant

```rust
pub fn constant(shape: &Arc<DualShape>, c: T) -> Self
```

**功能**：构造常量——值为 $c$，所有导数分量为零。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `&Arc<DualShape>` | 分量布局 |
| `c` | `T` | 常量值 |

**返回值**：`HyperDual<T>`

**示例**：

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(1);
let c = HyperDual::constant(&shape, Rational::new(7, 1));
assert_eq!(c.value(), &Rational::new(7, 1));
assert_eq!(c.deriv(0), Some(&Rational::new(0, 1)));
```

### HyperDual::value

```rust
pub fn value(&self) -> &T
```

**功能**：返回标量值分量（分量 0）的引用。

**返回值**：`&T`

### HyperDual::deriv

```rust
pub fn deriv(&self, i: usize) -> Option<&T>
```

**功能**：返回对变量 $i$ 的一阶偏导数 $\frac{\partial f}{\partial x_i}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `i` | `usize` | 变量编号（0-indexed） |

**返回值**：`Some(&T)` 若一阶分量 $[0,\dots,\underset{i}{1},\dots,0]$ 存在于形状中（对 `new_first_order` 形状等价于 $i < n\_vars$），否则 `None`。

### HyperDual::values

```rust
pub fn values(&self) -> &[T]
```

**功能**：返回所有分量的切片，按形状顺序排列。

**返回值**：`&[T]`

### HyperDual::shape

```rust
pub fn shape(&self) -> &Arc<DualShape>
```

**功能**：返回共享的形状引用。

**返回值**：`&Arc<DualShape>`

### HyperDual::from_values

```rust
pub fn from_values(shape: Arc<DualShape>, values: Vec<T>) -> Option<Self>
```

**功能**：从完整的分量向量构造。长度必须匹配形状的分量数。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `Arc<DualShape>` | 分量布局 |
| `values` | `Vec<T>` | 分量值，长度须为 `shape.n_components()` |

**返回值**：`Some(HyperDual<T>)` 若长度匹配，否则 `None`。

### HyperDual::zero / HyperDual::one

```rust
pub fn zero(shape: &Arc<DualShape>) -> Self
pub fn one(shape: &Arc<DualShape>) -> Self
```

**功能**：构造加法单位元（所有分量为零）和乘法单位元（值为 1，导数分量为零）。

### HyperDual::inv

```rust
pub fn inv(&self) -> Option<Self>
```

**功能**：乘法逆元 $\frac{1}{f}$，通过几何级数 $\frac{1}{v+\varepsilon} = \frac{1}{v}\sum_{k \geq 0}\left(-\frac{\varepsilon}{v}\right)^k$ 截断计算。

**返回值**：`Some(1/self)` 若值分量非零，否则 `None`（除以零）。

## 算术运算

`HyperDual<T>` 实现了以下标准 trait，所有运算自动传播导数：

| Trait | 运算 | 导数规则 |
|---|---|---|
| `Add` | `a + b` | $\frac{\partial}{\partial x_i}(a+b) = a'_i + b'_i$ |
| `Sub` | `a - b` | $\frac{\partial}{\partial x_i}(a-b) = a'_i - b'_i$ |
| `Neg` | `-a` | $\frac{\partial}{\partial x_i}(-a) = -a'_i$ |
| `Mul` | `a * b` | $\frac{\partial}{\partial x_i}(ab) = a'b + ab'$（通过乘法表处理高阶项） |
| `Div` | `a / b` | 等价于 `a * b.inv()`，要求 `b` 的值分量非零 |

**约束**：`Add`、`Sub`、`Mul`、`Div` 要求两个操作数共享同一 `DualShape`。debug 模式下校验分量数是否一致，不一致时 panic。

**除零**：`Div` 调用 `inv()`，若除数的值分量为零则 panic。

## 完整示例

### 示例 1：二元函数偏导数

计算 $f(x, y) = x \cdot y$ 在 $(3, 5)$ 处的值和偏导：

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(2);
let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));
let y = HyperDual::variable(&shape, 1, Rational::new(5, 1));

let f = x * y;

assert_eq!(f.value(), &Rational::new(15, 1));     // f(3,5) = 15
assert_eq!(f.deriv(0), Some(&Rational::new(5, 1))); // ∂f/∂x = y = 5
assert_eq!(f.deriv(1), Some(&Rational::new(3, 1))); // ∂f/∂y = x = 3
```

### 示例 2：幂函数导数

计算 $f(x) = x^3$ 在 $x=3$ 处的导数 $\frac{d}{dx}x^3 = 3x^2 = 27$：

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(1);
let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));

// x^3 = x * x * x
let x2 = x.clone() * x.clone();
let x3 = x2 * x;

assert_eq!(x3.value(), &Rational::new(27, 1));
assert_eq!(x3.deriv(0), Some(&Rational::new(27, 1)));
```

### 示例 3：商的导数

计算 $f(x,y) = \frac{x}{y}$ 在 $(6, 3)$ 处：$\frac{\partial f}{\partial x} = \frac{1}{y} = \frac{1}{3}$，$\frac{\partial f}{\partial y} = -\frac{x}{y^2} = -\frac{2}{3}$。

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(2);
let x = HyperDual::variable(&shape, 0, Rational::new(6, 1));
let y = HyperDual::variable(&shape, 1, Rational::new(3, 1));

let f = x / y;

assert_eq!(f.value(), &Rational::new(2, 1));
assert_eq!(f.deriv(0), Some(&Rational::new(1, 3)));
assert_eq!(f.deriv(1), Some(&Rational::new(-2, 3)));
```

### 示例 4：三变量乘积

计算 $f(x,y,z) = xyz$ 在 $(2, 3, 5)$ 处：$\frac{\partial f}{\partial x} = yz = 15$，$\frac{\partial f}{\partial y} = xz = 10$，$\frac{\partial f}{\partial z} = xy = 6$。

```rust
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(3);
let x = HyperDual::variable(&shape, 0, Rational::new(2, 1));
let y = HyperDual::variable(&shape, 1, Rational::new(3, 1));
let z = HyperDual::variable(&shape, 2, Rational::new(5, 1));

let f = x * y * z;

assert_eq!(f.value(), &Rational::new(30, 1));
assert_eq!(f.deriv(0), Some(&Rational::new(15, 1)));
assert_eq!(f.deriv(1), Some(&Rational::new(10, 1)));
assert_eq!(f.deriv(2), Some(&Rational::new(6, 1)));
```

### 示例 5：二阶导数

使用自定义形状追踪二阶导数。对于 $f(x) = x^2$，有 $f''(x) = 2$。注意：多指标为 $[2]$ 的分量存储的是 $\varepsilon^2$ 的系数，即 $f''(x)/2!$，而非 $f''(x)$ 本身。

```rust
use ocas_domain::dual::{DualShape, HyperDual};
use ocas_domain::Rational;

// 形状：[0]（值）、[1]（一阶）、[2]（二阶）
let shape = std::sync::Arc::new(
    DualShape::new(vec![vec![0], vec![1], vec![2]]).unwrap()
);
let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));

// f = x * x
let f = x.clone() * x;

assert_eq!(f.value(), &Rational::new(9, 1));     // 3^2 = 9
assert_eq!(f.deriv(0), Some(&Rational::new(6, 1))); // 2x = 6
// 分量 [2] 存储 ε² 的系数 = f''(3)/2! = 2/2 = 1
assert_eq!(f.values()[2], Rational::new(1, 1));  // f''(3)/2! = 1
```

## 实现细节

### 乘法表

`DualShape` 预计算一个乘法表：对于每对非标量分量 $(a, b)$，若其多指标之和也存在于布局中，记录三元组 $(a, b, c)$。乘法运算时，分量 $k$ 的结果为：

$$\text{result}[k] = a[k] \cdot b[0] + a[0] \cdot b[k] + \sum_{(i,j,k) \in \text{table}} a[i] \cdot b[j]$$

超出布局的高阶项被截断（truncation）。

### 逆元的几何级数

`inv()` 通过几何级数 $\frac{1}{v+\varepsilon} = \frac{1}{v}\sum_{p=0}^{\infty}\left(-\frac{\varepsilon}{v}\right)^p$ 计算，迭代直到高阶项为零（自动截断至形状精度）。

### 共享形状

所有参与同一表达式计算的 `HyperDual` 必须共享同一 `Arc<DualShape>`。形状通过 `Arc` 克隆成本极低。

## 限制

| 限制 | 说明 |
|---|---|
| **仅支持有理系数** | 当前 `DualCoeff` 仅对 `Rational` 实现，不支持 `f64` 或其他浮点类型 |
| **无超越函数** | `sin`、`cos`、`exp`、`ln` 等超越函数未实现——它们需要实系数 trait |
| **除零 panic** | `Div` 调用 `inv()`，值分量为零时 panic |
| **形状不匹配 panic** | `Add`/`Sub`/`Mul`/`Div` 要求操作数共享同一形状（debug 模式下校验分量数），不一致时 panic |
| **高阶导数需手动布局** | 默认的 `new_first_order` 仅追踪一阶偏导；二阶及以上需通过 `DualShape::new` 手动构建布局 |

## 参见

- [系数域](./rust-domains.md) — `Rational` 及其他域类型
- [求值与 JIT](./rust-evaluation.md) — 数值求值（浮点）
- [微积分](./rust-calculus.md) — 符号微分 `diff`、积分 `integrate`
