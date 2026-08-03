# Automatic Differentiation

oCAS implements forward-mode automatic differentiation via hyper-dual numbers. A `HyperDual<T>` carries a scalar value and a set of derivative components. The component layout is precomputed by a shared `DualShape` (shared through an `Arc`), and arithmetic operations use a prebuilt multiplication table.

> **Source file**: `ocas-domain/src/dual.rs`

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

**Description**: defines the coefficient type constraints required by hyper-dual numbers. Requires the full four arithmetic operations (including multiplicative inverses) plus additive and multiplicative identities.

**Parameters**: none (trait bound).

**Built-in implementations**:

| Type | Description |
|---|---|
| `Rational` | arbitrary-precision rational; `zero()` = `0/1`, `one()` = `1/1` |

**Limitation**: currently only `Rational` implements `DualCoeff`. `f64` and other floating-point types are not supported as coefficients.

## DualShape

```rust
#[derive(Debug, Clone)]
pub struct DualShape { /* private fields */ }
```

**Description**: describes the component layout of a hyper-dual number. Each component is identified by a multi-index — a vector of non-negative integers of length equal to the number of variables, recording the order of differentiation with respect to each variable. The layout must satisfy **ancestor closure**: if the multi-index $\mathbf{m}$ is present, then every multi-index that is component-wise ≤ $\mathbf{m}$ must also be present.

Component 0 is always the all-zero multi-index (the scalar value component).

### DualShape::new

```rust
pub fn new(components: Vec<Vec<usize>>) -> Option<Self>
```

**Description**: builds a layout from an ancestor-closed list of multi-indices.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `components` | `Vec<Vec<usize>>` | the list of multi-indices; each `Vec<usize>` must have the same length (missing entries are zero-padded automatically) |

**Returns**: `Some(DualShape)` if the layout is valid, otherwise `None`.

**Error conditions**:
- empty list → `None`
- no all-zero multi-index → `None`
- ancestor closure violated → `None`

**Example**:

```rust
use ocas_domain::dual::DualShape;

// First order, two variables: value + ∂/∂x₀ + ∂/∂x₁
let shape = DualShape::new(vec![
    vec![0, 0],
    vec![1, 0],
    vec![0, 1],
]).unwrap();
assert_eq!(shape.n_components(), 3);
assert_eq!(shape.n_vars(), 2);

// Invalid: [1] missing (not ancestor-closed)
assert!(DualShape::new(vec![vec![0], vec![2]]).is_none());
```

### DualShape::n_vars

```rust
pub fn n_vars(&self) -> usize
```

**Description**: returns the number of differentiation variables (the length of the multi-indices).

**Returns**: `usize`

### DualShape::n_components

```rust
pub fn n_components(&self) -> usize
```

**Description**: returns the total number of components (including the scalar value component).

**Returns**: `usize`

### DualShape::components

```rust
pub fn components(&self) -> &[Vec<usize>]
```

**Description**: returns the slice of component multi-indices. Index 0 is the scalar value component.

**Returns**: `&[Vec<usize>]`

### DualShape::index_of

```rust
pub fn index_of(&self, multi_index: &[usize]) -> Option<usize>
```

**Description**: looks up the component index for a given multi-index.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `multi_index` | `&[usize]` | the multi-index to look up |

**Returns**: `Some(index)` if present, otherwise `None`.

### DualShape::mult_table

```rust
pub fn mult_table(&self) -> &[(usize, usize, usize)]
```

**Description**: returns the multiplication table (for debugging/testing). Each triple `(a, b, c)` means that the product of component `a` and component `b` contributes to component `c`. Pairs involving the scalar component (index 0) are not included.

**Returns**: `&[(usize, usize, usize)]`

## new_first_order

```rust
pub fn new_first_order<T: DualCoeff>(nvars: usize) -> Arc<DualShape>
```

**Description**: builds a first-order shape tracking the partial derivative $\frac{\partial}{\partial x_i}$ of each variable. Does not track higher-order or mixed partials.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `nvars` | `usize` | the number of differentiation variables |

**Returns**: `Arc<DualShape>` — contains $n+1$ components: the value component $[0,\dots,0]$ plus one component $[0,\dots,\underset{i}{1},\dots,0]$ per variable.

**Example**:

```rust
use ocas_domain::dual::new_first_order;
use ocas_domain::Rational;

let shape = new_first_order::<Rational>(3);
assert_eq!(shape.n_components(), 4); // value + 3 partials
assert_eq!(shape.n_vars(), 3);
```

## HyperDual\<T\>

```rust
#[derive(Debug, Clone)]
pub struct HyperDual<T: DualCoeff> { /* private fields */ }
```

**Description**: the hyper-dual number type, carrying a scalar value and derivative components laid out according to a `DualShape`. All arithmetic operations propagate derivatives automatically.

**Type constraint**: `T: DualCoeff`

### HyperDual::variable

```rust
pub fn variable(shape: &Arc<DualShape>, i: usize, c: T) -> Self
```

**Description**: constructs the $i$-th independent variable with value $c$. The derivative component $[0,\dots,\underset{i}{1},\dots,0]$ is set to 1; all others are 0.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `shape` | `&Arc<DualShape>` | the component layout |
| `i` | `usize` | the variable index (0-based), must be < `n_vars()` (out-of-range indices are silently ignored, no panic) |
| `c` | `T` | the variable's value at the evaluation point |

**Returns**: `HyperDual<T>`

**Example**:

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

**Description**: constructs a constant — value $c$, all derivative components zero.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `shape` | `&Arc<DualShape>` | the component layout |
| `c` | `T` | the constant value |

**Returns**: `HyperDual<T>`

**Example**:

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

**Description**: returns a reference to the scalar value component (component 0).

**Returns**: `&T`

### HyperDual::deriv

```rust
pub fn deriv(&self, i: usize) -> Option<&T>
```

**Description**: returns the first-order partial derivative $\frac{\partial f}{\partial x_i}$ with respect to variable $i$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `i` | `usize` | the variable index (0-based) |

**Returns**: `Some(&T)` if the first-order component $[0,\dots,\underset{i}{1},\dots,0]$ exists in the shape (for `new_first_order` shapes this is equivalent to $i < n\_vars$), otherwise `None`.

### HyperDual::values

```rust
pub fn values(&self) -> &[T]
```

**Description**: returns a slice of all components, in shape order.

**Returns**: `&[T]`

### HyperDual::shape

```rust
pub fn shape(&self) -> &Arc<DualShape>
```

**Description**: returns the shared shape reference.

**Returns**: `&Arc<DualShape>`

### HyperDual::from_values

```rust
pub fn from_values(shape: Arc<DualShape>, values: Vec<T>) -> Option<Self>
```

**Description**: constructs from a full component vector. The length must match the shape's component count.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `shape` | `Arc<DualShape>` | the component layout |
| `values` | `Vec<T>` | the component values; length must be `shape.n_components()` |

**Returns**: `Some(HyperDual<T>)` if the length matches, otherwise `None`.

### HyperDual::zero / HyperDual::one

```rust
pub fn zero(shape: &Arc<DualShape>) -> Self
pub fn one(shape: &Arc<DualShape>) -> Self
```

**Description**: construct the additive identity (all components zero) and the multiplicative identity (value 1, derivative components zero).

### HyperDual::inv

```rust
pub fn inv(&self) -> Option<Self>
```

**Description**: the multiplicative inverse $\frac{1}{f}$, computed by truncating the geometric series $\frac{1}{v+\varepsilon} = \frac{1}{v}\sum_{k \geq 0}\left(-\frac{\varepsilon}{v}\right)^k$.

**Returns**: `Some(1/self)` if the value component is nonzero, otherwise `None` (division by zero).

## Arithmetic operations

`HyperDual<T>` implements the following standard traits; all operations propagate derivatives automatically:

| Trait | Operation | Derivative rule |
|---|---|---|
| `Add` | `a + b` | $\frac{\partial}{\partial x_i}(a+b) = a'_i + b'_i$ |
| `Sub` | `a - b` | $\frac{\partial}{\partial x_i}(a-b) = a'_i - b'_i$ |
| `Neg` | `-a` | $\frac{\partial}{\partial x_i}(-a) = -a'_i$ |
| `Mul` | `a * b` | $\frac{\partial}{\partial x_i}(ab) = a'b + ab'$ (higher-order terms handled via the multiplication table) |
| `Div` | `a / b` | equivalent to `a * b.inv()`; requires the value component of `b` to be nonzero |

**Constraint**: `Add`, `Sub`, `Mul`, `Div` require both operands to share the same `DualShape`. In debug mode the component counts are checked, and mismatches panic.

**Division by zero**: `Div` calls `inv()`, which panics if the divisor's value component is zero.

## Complete examples

### Example 1: partial derivatives of a two-variable function

Compute the value and partial derivatives of $f(x, y) = x \cdot y$ at $(3, 5)$:

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

### Example 2: derivative of a power function

Compute the derivative of $f(x) = x^3$ at $x=3$: $\frac{d}{dx}x^3 = 3x^2 = 27$:

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

### Example 3: derivative of a quotient

Compute $f(x,y) = \frac{x}{y}$ at $(6, 3)$: $\frac{\partial f}{\partial x} = \frac{1}{y} = \frac{1}{3}$, $\frac{\partial f}{\partial y} = -\frac{x}{y^2} = -\frac{2}{3}$.

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

### Example 4: three-variable product

Compute $f(x,y,z) = xyz$ at $(2, 3, 5)$: $\frac{\partial f}{\partial x} = yz = 15$, $\frac{\partial f}{\partial y} = xz = 10$, $\frac{\partial f}{\partial z} = xy = 6$.

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

### Example 5: second derivative

Use a custom shape to track second derivatives. For $f(x) = x^2$, $f''(x) = 2$. Note: the component with multi-index $[2]$ stores the coefficient of $\varepsilon^2$, i.e. $f''(x)/2!$, not $f''(x)$ itself.

```rust
use ocas_domain::dual::{DualShape, HyperDual};
use ocas_domain::Rational;

// Shape: [0] (value), [1] (first order), [2] (second order)
let shape = std::sync::Arc::new(
    DualShape::new(vec![vec![0], vec![1], vec![2]]).unwrap()
);
let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));

// f = x * x
let f = x.clone() * x;

assert_eq!(f.value(), &Rational::new(9, 1));     // 3^2 = 9
assert_eq!(f.deriv(0), Some(&Rational::new(6, 1))); // 2x = 6
// Component [2] stores the coefficient of ε² = f''(3)/2! = 2/2 = 1
assert_eq!(f.values()[2], Rational::new(1, 1));  // f''(3)/2! = 1
```

## Implementation details

### Multiplication table

`DualShape` precomputes a multiplication table: for every pair of non-scalar components $(a, b)$, if the sum of their multi-indices is also present in the layout, record the triple $(a, b, c)$. During multiplication, component $k$ is computed as:

$$\text{result}[k] = a[k] \cdot b[0] + a[0] \cdot b[k] + \sum_{(i,j,k) \in \text{table}} a[i] \cdot b[j]$$

Higher-order terms beyond the layout are truncated.

### Geometric series for the inverse

`inv()` computes via the geometric series $\frac{1}{v+\varepsilon} = \frac{1}{v}\sum_{p=0}^{\infty}\left(-\frac{\varepsilon}{v}\right)^p$, iterating until the higher-order terms vanish (automatically truncated to the shape's precision).

### Shared shape

All `HyperDual`s participating in the same expression must share the same `Arc<DualShape>`. Cloning the shape via `Arc` is very cheap.

## Limitations

| Limitation | Description |
|---|---|
| **Rational coefficients only** | `DualCoeff` is currently implemented only for `Rational`; `f64` and other floating-point types are not supported |
| **No transcendental functions** | transcendental functions such as `sin`, `cos`, `exp`, `ln` are not implemented — they require real-coefficient traits |
| **Division-by-zero panic** | `Div` calls `inv()`, which panics when the value component is zero |
| **Shape-mismatch panic** | `Add`/`Sub`/`Mul`/`Div` require the operands to share the same shape (component counts are checked in debug mode); mismatches panic |
| **Higher-order derivatives need a manual layout** | the default `new_first_order` tracks only first-order partials; second and higher orders require building the layout manually via `DualShape::new` |

## See also

- [Coefficient domains](./rust-domains.md) — `Rational` and other domain types
- [Evaluation and JIT](./rust-evaluation.md) — numerical evaluation (floating point)
- [Calculus](./rust-calculus.md) — symbolic differentiation `diff`, integration `integrate`
