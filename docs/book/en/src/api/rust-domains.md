# Rust API Reference: Coefficient Domains

Coefficient domains are the foundational abstraction behind all of oCAS's polynomial, matrix, and solver algorithms. Each domain describes a set of elements and their arithmetic operations, unified through the `Domain` and `EuclideanDomain` traits.

**Module path**: `ocas_domain`

**Imports**:

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

**Signature**:

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

**Description**: The core trait for coefficient domains. The domain object itself may carry parameters (such as the modulus of a finite field), so all operations go through `&self`. This follows the "domain object" pattern used by CAS libraries such as Flint and SymPy's `Domain`.

**Associated types**:

| Associated type | Bounds | Description |
|---|---|---|
| `Element` | `Clone + PartialEq + Eq + Debug + 'static` | The type of elements in the domain |

**Methods**:

| Method | Signature | Description |
|---|---|---|
| `zero` | `fn zero(&self) -> Self::Element` | Additive identity |
| `one` | `fn one(&self) -> Self::Element` | Multiplicative identity |
| `add` | `fn add(&self, a, b) -> Self::Element` | Addition |
| `sub` | `fn sub(&self, a, b) -> Self::Element` | Subtraction |
| `neg` | `fn neg(&self, a) -> Self::Element` | Negation |
| `mul` | `fn mul(&self, a, b) -> Self::Element` | Multiplication |
| `div` | `fn div(&self, a, b) -> Option<Self::Element>` | Division. Returns `None` if `b` is zero or the division is not exact |
| `inv` | `fn inv(&self, a) -> Option<Self::Element>` | Multiplicative inverse. Returns `None` if `a` is zero |
| `is_zero` | `fn is_zero(&self, a) -> bool` | Whether `a` is the additive identity (default: `*a == self.zero()`) |
| `is_one` | `fn is_one(&self, a) -> bool` | Whether `a` is the multiplicative identity (default: `*a == self.one()`) |
| `mul_assign` | `fn mul_assign(&self, a: &mut E, b: &E)` | In-place multiplication `*a *= b`. Defaults to creating a new element; high-performance domains such as GMP may override |
| `sub_mul_assign` | `fn sub_mul_assign(&self, a: &mut E, b: &E, c: &E)` | Fused subtract-multiply `*a -= b * c`. Used heavily in F4 row echelonization |
| `pow` | `fn pow(&self, a, n: u64) -> Self::Element` | Non-negative integer power. Defaults to binary exponentiation; finite fields override with `modpow` |
| `cast_u64` | `fn cast_u64(&self, n: u64) -> Self::Element` | Converts a `u64` to a domain element. Defaults to adding one repeatedly |

**Example**:

```rust
use ocas_domain::{Domain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(3);
let b = Integer::from(5);
assert_eq!(domain.add(&a, &b), Integer::from(8));
assert_eq!(domain.mul(&a, &b), Integer::from(15));
assert_eq!(domain.pow(&a, 3), Integer::from(27));
// Output: all assertions pass
```

**See also**: [EuclideanDomain](#euclideandomain-trait)

---

## EuclideanDomain trait

**Signature**:

```rust
pub trait EuclideanDomain: Domain {
    fn div_rem(&self, a: &Self::Element, b: &Self::Element)
        -> Option<(Self::Element, Self::Element)>;
    fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn extended_gcd(&self, a: &Self::Element, b: &Self::Element)
        -> (Self::Element, Self::Element, Self::Element);
}
```

**Description**: A Euclidean domain supporting division with remainder. An extension of `Domain` that additionally provides `div_rem` (division with remainder), `gcd` (greatest common divisor), and `extended_gcd` (extended Euclidean algorithm).

**Methods**:

### div_rem

**Signature**: `fn div_rem(&self, a: &Self::Element, b: &Self::Element) -> Option<(Self::Element, Self::Element)>`

**Description**: Division with remainder, returning `(quotient, remainder)` such that `a = quotient * b + remainder` and `remainder == 0` or `deg(remainder) < deg(b)`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `&Self::Element` | Dividend |
| `b` | `&Self::Element` | Divisor |

**Return value**: `Some((quotient, remainder))`, or `None` when `b` is zero.

**Example**:

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(17);
let b = Integer::from(5);
let (q, r) = domain.div_rem(&a, &b).unwrap();
assert_eq!(q, Integer::from(3));
assert_eq!(r, Integer::from(2));
// Output: 17 = 3 × 5 + 2
```

### gcd

**Signature**: `fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element`

**Description**: Computes the greatest common divisor. The default implementation uses the Euclidean algorithm.

**Return value**: `gcd(a, b)`. For a field (such as `FiniteField`), the GCD of two non-zero elements degenerates to 1.

**Example**:

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(12);
let b = Integer::from(8);
let g = domain.gcd(&a, &b);
assert_eq!(g, Integer::from(4));
// Output: gcd(12, 8) = 4
```

### extended_gcd

**Signature**: `fn extended_gcd(&self, a: &Self::Element, b: &Self::Element) -> (Self::Element, Self::Element, Self::Element)`

**Description**: Extended Euclidean algorithm, returning `(g, x, y)` such that `g = gcd(a, b) = a * x + b * y`.

**Return value**: The triple `(g, x, y)`.

**Example**:

```rust
use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(35);
let b = Integer::from(15);
let (g, x, y) = domain.extended_gcd(&a, &b);
// g = 5, x = 1, y = -2, so that 5 = 35×1 + 15×(-2)
assert_eq!(g, Integer::from(5));
```

**See also**: [Domain](#domain-trait)

---

## Integer / IntegerDomain

### IntegerDomain

**Signature**: `pub struct IntegerDomain;`

**Description**: The integer domain $\mathbb{Z}$, with element type `Integer`. Implements `Domain` and `EuclideanDomain`.

**Characteristics**:

- `div` requires exact division (zero remainder), otherwise returns `None`
- `inv` returns a result only for $\pm 1$, otherwise `None`
- `pow` uses binary exponentiation

**Example**:

```rust
use ocas_domain::{Domain, EuclideanDomain, Integer, IntegerDomain};

let domain = IntegerDomain;
let a = Integer::from(10);
let b = Integer::from(3);

// Exact division: 10/3 is not exact
assert!(domain.div(&a, &b).is_none());

// Division with remainder
let (q, r) = domain.div_rem(&a, &b).unwrap();
assert_eq!(q, Integer::from(3));
assert_eq!(r, Integer::from(1));

// Inverse: only ±1 have inverses
assert!(domain.inv(&Integer::from(1)).is_some());
assert!(domain.inv(&Integer::from(5)).is_none());
// Output: all assertions pass
```

### Integer

**Signature**: `pub struct Integer(BigInt);`

**Description**: An arbitrary-precision integer. The default build uses `num-bigint`'s `BigInt`; with the `gmp` feature enabled, a GMP backend is used.

**Derives**: `Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash`

#### Construction and conversion

##### Integer::new

**Signature**: `pub fn new<T: Into<BigInt>>(value: T) -> Self`

**Description**: Creates an arbitrary-precision integer from a machine integer or a `BigInt`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `value` | `T: Into<BigInt>` | The integer value (`i32`, `i64`, `u64`, `BigInt`, etc.) |

**Example**:

```rust
use ocas_domain::Integer;

let a = Integer::new(42);
let b = Integer::new(100_i64);
assert_eq!(a.to_string(), "42");
// Output: 42
```

##### Integer::from

**Signature**: `impl From<i64> for Integer` / `impl From<BigInt> for Integer`

**Description**: Conversion from `i64` or `BigInt`.

**Example**:

```rust
use ocas_domain::Integer;

let a = Integer::from(42);
assert_eq!(a.to_string(), "42");
// Output: 42
```

#### Accessors

##### Integer::inner

**Signature**: `pub fn inner(&self) -> &BigInt`

**Description**: Accesses the underlying `BigInt` reference.

##### Integer::to_bigint

**Signature**: `pub fn to_bigint(&self) -> BigInt`

**Description**: Clones into a `BigInt` (regardless of backend).

##### Integer::to_i64

**Signature**: `pub fn to_i64(&self) -> Option<i64>`

**Description**: Attempts conversion to `i64`. Returns `None` on overflow.

#### Arithmetic methods

##### Integer::pow_u32

**Signature**: `pub fn pow_u32(&self, exp: u32) -> Self`

**Description**: Computes $n^{exp}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `exp` | `u32` | Non-negative exponent |

**Example**:

```rust
use ocas_domain::Integer;

let a = Integer::from(2);
assert_eq!(a.pow_u32(10).to_string(), "1024");
// Output: 1024
```

##### Integer::modpow

**Signature**: `pub fn modpow(&self, exp: &Integer, modulus: &Integer) -> Integer`

**Description**: Modular exponentiation $self^{exp} \bmod modulus$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `exp` | `&Integer` | Exponent (non-negative) |
| `modulus` | `&Integer` | Modulus (positive) |

**Example**:

```rust
use ocas_domain::Integer;

let base = Integer::from(3);
let exp = Integer::from(100);
let modulus = Integer::from(7);
let result = base.modpow(&exp, &modulus);
assert_eq!(result.to_string(), "4");
// Output: 3^100 mod 7 = 4
```

##### Integer::mod_floor

**Signature**: `pub fn mod_floor(&self, modulus: &Integer) -> Integer`

**Description**: Floor modulus; the result $r$ satisfies $0 \leq r < |modulus|$.

**Example**:

```rust
use ocas_domain::Integer;

let a = Integer::from(-7);
let m = Integer::from(3);
assert_eq!(a.mod_floor(&m).to_string(), "2");
// Output: -7 mod_floor 3 = 2
```

##### Integer::div_rem

**Signature**: `pub fn div_rem(&self, other: &Integer) -> (Integer, Integer)`

**Description**: Division with remainder, `(quotient, remainder)`. Note: unlike `EuclideanDomain::div_rem`, this method does not require `other` to be non-zero (behavior is determined by the underlying `BigInt`).

##### Integer::is_even

**Signature**: `pub fn is_even(&self) -> bool`

**Description**: Whether the integer is even.

##### Integer::is_negative

**Signature**: `pub fn is_negative(&self) -> bool`

**Description**: Whether the integer is negative.

##### Integer::is_zero

**Signature**: `pub fn is_zero(&self) -> bool`

**Description**: Whether the integer is zero.

##### Integer::is_one

**Signature**: `pub fn is_one(&self) -> bool`

**Description**: Whether the integer is one.

##### Integer::abs

**Signature**: `pub fn abs(&self) -> Integer`

**Description**: Absolute value.

##### Integer::sqrt

**Signature**: `pub fn sqrt(&self) -> Integer`

**Description**: Integer square root (floored).

**Example**:

```rust
use ocas_domain::Integer;

let a = Integer::from(10);
assert_eq!(a.sqrt().to_string(), "3");
// Output: sqrt(10) = 3
```

#### Operators

`Integer` implements the following standard operator traits (all reference combinations are supported: owned × owned, owned × &, & × owned, & × &):

| Trait | Operation |
|---|---|
| `Add` | `+` |
| `Sub` | `-` |
| `Mul` | `*` |
| `Div` | `/` (integer division, truncated toward zero) |
| `Rem` | `%` |
| `Neg` | unary `-` |
| `Shr<u32>` / `ShrAssign<u32>` | right shift |
| `AddAssign<&Integer>` / `SubAssign<&Integer>` / `MulAssign<&Integer>` / `DivAssign<&Integer>` | compound assignment |

**See also**: [Domain](#domain-trait), [Rational](#rational--rationaldomain), [Number-theoretic functions](rust-ntheory.md)

---

## Rational / RationalDomain

### RationalDomain

**Signature**: `pub struct RationalDomain;`

**Description**: The rational number field $\mathbb{Q}$, with element type `Rational`. Implements `Domain` and `EuclideanDomain`.

**Characteristics**:

- `div` returns an exact result for every non-zero divisor
- `inv` returns a result for every non-zero element
- `EuclideanDomain`'s `div_rem` degenerates over the rationals: the remainder is always zero

**Example**:

```rust
use ocas_domain::{Domain, Rational, RationalDomain};

let domain = RationalDomain;
let a = Rational::new(1, 2);
let b = Rational::new(1, 3);
let sum = domain.add(&a, &b);
assert_eq!(sum, Rational::new(5, 6));
// Output: 1/2 + 1/3 = 5/6
```

### Rational

**Signature**: `pub struct Rational(BigRational);`

**Description**: An arbitrary-precision rational number. The default build uses `num-rational`'s `BigRational`; with the `gmp` feature enabled, a GMP backend is used.

**Derives**: `Debug, Clone, PartialEq, Eq, Hash`

#### Construction

##### Rational::new

**Signature**: `pub fn new(numer: i64, denom: i64) -> Self`

**Description**: Creates a rational from a numerator and denominator (`i64`). Automatically reduced to lowest terms. Behavior for a zero denominator is determined by the underlying library.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `numer` | `i64` | Numerator |
| `denom` | `i64` | Denominator (non-zero) |

**Example**:

```rust
use ocas_domain::Rational;

let a = Rational::new(3, 6);
assert_eq!(a.to_string(), "1/2");
// Output: 3/6 automatically reduces to 1/2
```

##### Rational::from_bigints

**Signature**: `pub fn from_bigints(numer: BigInt, denom: BigInt) -> Self`

**Description**: Creates a rational from an arbitrary-precision integer numerator and denominator.

##### Rational::from_integer

**Signature**: `pub fn from_integer(n: Integer) -> Self`

**Description**: Creates a rational from an integer (denominator = 1).

#### Accessors

##### Rational::inner

**Signature**: `pub fn inner(&self) -> &BigRational`

**Description**: Accesses the underlying `BigRational` reference.

##### Rational::numer

**Signature**: `pub fn numer(&self) -> Integer`

**Description**: Returns the numerator (as an `Integer`).

##### Rational::denom

**Signature**: `pub fn denom(&self) -> Integer`

**Description**: Returns the denominator (as an `Integer`, always positive).

**Example**:

```rust
use ocas_domain::Rational;

let r = Rational::new(3, 4);
assert_eq!(r.numer().to_string(), "3");
assert_eq!(r.denom().to_string(), "4");
// Output: numerator = 3, denominator = 4
```

#### Operators

`Rational` implements the following standard operator traits:

| Trait | Operation |
|---|---|
| `Add` | `+` |
| `Sub` | `-` |
| `Mul` | `*` |
| `Div` | `/` |
| `Neg` | unary `-` |
| `AddAssign` / `SubAssign` / `MulAssign` / `DivAssign` | compound assignment |

**See also**: [Domain](#domain-trait), [Integer](#integer--integerdomain)

---

## FiniteField / FiniteFieldElement

### FiniteField

**Signature**: `pub struct FiniteField { /* prime: BigInt, ... */ }`

**Description**: The prime finite field $\mathbb{Z}/p\mathbb{Z}$. Arithmetic uses arbitrary-precision integers, supporting large primes. Implements `Domain` and `EuclideanDomain`.

#### Construction

##### FiniteField::new

**Signature**: `pub fn new(prime: BigInt) -> Self`

**Description**: Creates a finite field with the given prime modulus.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `prime` | `BigInt` | The modulus (must be $\geq 2$; primality is not verified) |

**Errors**: panics in debug mode if `prime < 2`.

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};

let f = FiniteField::new(BigInt::from(7));
let a = f.element(3);
let b = f.element(5);
assert_eq!(f.add(&a, &b), f.element(1));   // 3 + 5 = 1 (mod 7)
assert_eq!(f.mul(&a, &b), f.element(1));   // 3 × 5 = 1 (mod 7)
// Output: all assertions pass
```

#### Element construction

##### FiniteField::element

**Signature**: `pub fn element(&self, value: impl Into<BigInt>) -> FiniteFieldElement`

**Description**: Creates a field element from an arbitrary integer. The value is automatically reduced to $[0, p-1]$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `value` | `impl Into<BigInt>` | The integer value (may lie outside $[0, p-1]$; reduced automatically) |

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};

let f = FiniteField::new(BigInt::from(7));
let a = f.element(10);
assert_eq!(a.value().to_string(), "3");
// Output: 10 mod 7 = 3
```

##### FiniteField::from_i64

**Signature**: `pub fn from_i64(&self, val: i64) -> FiniteFieldElement`

**Description**: Creates a field element from an `i64` (reduced mod p).

#### Accessors

##### FiniteField::prime

**Signature**: `pub fn prime(&self) -> &BigInt`

**Description**: Returns the field's modulus.

##### FiniteField::prime_u64

**Signature**: `pub fn prime_u64(&self) -> u64`

**Description**: Returns the modulus as a `u64`.

**Errors**: panics if the prime does not fit in a `u64`.

##### FiniteField::to_i64

**Signature**: `pub fn to_i64(&self, a: &FiniteFieldElement) -> i64`

**Description**: Converts a field element to an `i64` (in the range $[0, p)$).

**Errors**: panics if the prime does not fit in a `u64`.

#### Domain implementation

- `div` inverts via Fermat's little theorem: $a^{-1} \equiv a^{p-2} \pmod{p}$
- `inv` as above; returns `None` for the zero element
- `pow` uses the `modpow` optimization (much faster than the default binary exponentiation)

#### EuclideanDomain implementation

- `div_rem`: division is exact in a field; the remainder is always zero
- `gcd`: GCD degenerates in a field — returns 0 when both are zero, otherwise 1

**See also**: [Domain](#domain-trait), [Finite field mathematics](../math/finite-fields.md)

---

### FiniteFieldElement

**Signature**: `pub struct FiniteFieldElement { value: BigInt }`

**Description**: An element of a prime finite field. The value is always in $[0, p-1]$.

**Derives**: `Debug, Clone, PartialEq, Eq, Hash`

#### FiniteFieldElement::value

**Signature**: `pub fn value(&self) -> &BigInt`

**Description**: Returns the canonical representative in $[0, p-1]$.

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::FiniteField;

let f = FiniteField::new(BigInt::from(7));
let a = f.element(-3);
assert_eq!(a.value().to_string(), "4");
// Output: -3 mod 7 = 4
```

**See also**: [FiniteField](#finitefield--finitefieldelement)

---

## RealBall / RealBallDomain

### RealBallDomain

**Signature**: `pub struct RealBallDomain;`

**Description**: The real ball (interval) domain. Element type `RealBall`. Implements only `Domain` (not `EuclideanDomain`, because real balls do not support exact division with remainder).

**Note**: The default build uses lightweight `f64` balls, suitable for templates and demos. With the `mpfr` feature enabled, `rug::Float` with directed rounding produces rigorous intervals.

### RealBall

**Signature**: `pub struct RealBall { mid: f64, rad: f64 }` (default build)

**Description**: A real ball: midpoint ± radius. The true value is guaranteed to lie in $[mid - rad, mid + rad]$.

**Derives**: `Debug, Clone, PartialEq, Eq`; the default build additionally implements `Copy`.

#### Construction

##### RealBall::new

**Signature** (default): `pub fn new(mid: f64, rad: f64) -> Self`

**Description**: Creates a ball from a midpoint and radius. The radius is clamped to be non-negative.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `mid` | `f64` (or `rug::Float`) | Midpoint |
| `rad` | `f64` (or `rug::Float`) | Radius ($\geq 0$) |

##### RealBall::from_f64

**Signature**: `pub fn from_f64(value: f64) -> Self`

**Description**: Creates a zero-radius (exact) ball from an `f64` value.

**Example**:

```rust
use ocas_domain::RealBall;

let ball = RealBall::from_f64(3.14);
assert_eq!(ball.mid(), 3.14);
assert_eq!(ball.rad(), 0.0);
// Output: exact ball
```

#### Accessors

##### RealBall::mid

**Signature**: `pub fn mid(&self) -> f64` (default) / `pub fn mid(&self) -> &rug::Float` (mpfr)

**Description**: Returns the midpoint.

##### RealBall::rad

**Signature**: `pub fn rad(&self) -> f64` (default) / `pub fn rad(&self) -> &rug::Float` (mpfr)

**Description**: Returns the radius.

##### RealBall::lower

**Signature**: `pub fn lower(&self) -> f64` (default) / `pub fn lower(&self) -> rug::Float` (mpfr)

**Description**: Returns the conservative lower bound $mid - rad$ (the mpfr version rounds downward).

##### RealBall::upper

**Signature**: `pub fn upper(&self) -> f64` (default) / `pub fn upper(&self) -> rug::Float` (mpfr)

**Description**: Returns the conservative upper bound $mid + rad$ (the mpfr version rounds upward).

##### RealBall::precision

**Signature**: `pub fn precision(&self) -> u32` (`mpfr` feature only)

**Description**: Returns the precision (in bits) of the MPFR backend.

#### Domain implementation

- `add`: $(a \pm r_a) + (b \pm r_b) = (a+b) \pm (r_a + r_b)$ (the mpfr version additionally adds rounding error)
- `sub`: like addition; the radii add
- `mul`: four-corner method — compute the four extreme products $a_lo \cdot b_lo$, $a_lo \cdot b_hi$, etc., and take the min/max for the new interval
- `div`: $a / b = a \cdot b^{-1}$, returns `None` when the ball contains zero
- `inv`: $1/(mid \pm rad)$, returns `None` when the ball contains zero

**Example**:

```rust
use ocas_domain::{Domain, RealBall, RealBallDomain};

let domain = RealBallDomain;
let a = RealBall::from_f64(2.0);
let b = RealBall::from_f64(3.0);
let prod = domain.mul(&a, &b);
assert!(prod.lower() <= 6.0 && 6.0 <= prod.upper());
// Output: the product ball contains the true value 6.0
```

**See also**: [Domain](#domain-trait)

---

## Complex / ComplexDomain

### ComplexDomain

**Signature**: `pub struct ComplexDomain<D: Domain> { base: D, ... }`

**Description**: The complex field built over an arbitrary base domain $D$. Element type `Complex<D>`. Implements only `Domain`.

#### Construction

##### ComplexDomain::new

**Signature**: `pub fn new(base: D) -> Self`

**Description**: Creates the complex field over the base domain `base`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `base` | `D: Domain` | The base domain (e.g. `IntegerDomain`, `RationalDomain`) |

##### ComplexDomain::base

**Signature**: `pub fn base(&self) -> &D`

**Description**: Returns a reference to the base domain.

##### ComplexDomain::real_element

**Signature**: `pub fn real_element(&self, re: D::Element) -> Complex<D>`

**Description**: Creates a purely real element.

##### ComplexDomain::imag_element

**Signature**: `pub fn imag_element(&self, im: D::Element) -> Complex<D>`

**Description**: Creates a purely imaginary element.

#### Domain implementation

- `mul`: $(a+bi)(c+di) = (ac-bd) + (ad+bc)i$
- `div`: $(a+bi)/(c+di) = \frac{(ac+bd)+(bc-ad)i}{c^2+d^2}$, returns `None` when the denominator is zero
- `inv`: implemented via `div(one, a)`

**Example**:

```rust
use ocas_domain::{Complex, ComplexDomain, Domain, Integer, IntegerDomain};

let domain = ComplexDomain::new(IntegerDomain);
let a = Complex::new(Integer::from(1), Integer::from(2));
let b = Complex::new(Integer::from(3), Integer::from(4));
let sum = domain.add(&a, &b);
assert_eq!(*sum.re(), Integer::from(4));
assert_eq!(*sum.im(), Integer::from(6));
// Output: (1+2i) + (3+4i) = 4+6i
```

**See also**: [Domain](#domain-trait)

---

### Complex

**Signature**: `pub struct Complex<D: Domain> { inner: NumComplex<D::Element> }`

**Description**: A complex number whose real and imaginary parts belong to the base domain $D$.

**Derives**: `Debug, Clone, PartialEq, Eq, Hash`

#### Construction and access

##### Complex::new

**Signature**: `pub fn new(real: D::Element, imag: D::Element) -> Self`

**Description**: Creates a complex number from its real and imaginary parts.

##### Complex::re

**Signature**: `pub fn re(&self) -> &D::Element`

**Description**: Returns a reference to the real part.

##### Complex::im

**Signature**: `pub fn im(&self) -> &D::Element`

**Description**: Returns a reference to the imaginary part.

##### Complex::inner

**Signature**: `pub fn inner(&self) -> &NumComplex<D::Element>`

**Description**: Returns a reference to the underlying `num_complex::Complex`.

**See also**: [ComplexDomain](#complexdomain)

---

## DoubleF64 / DoubleF64Domain

### DoubleF64Domain

**Signature**: `pub struct DoubleF64Domain;`

**Description**: The double-double floating-point domain. Element type `DoubleF64`. Implements only `Domain`.

### DoubleF64

**Signature**:

```rust
pub struct DoubleF64 {
    pub hi: f64,  // high component (principal value)
    pub lo: f64,  // low component (error term)
}
```

**Description**: A double-double floating-point number represented as $hi + lo$ with $|lo| \leq 0.5 \cdot \text{ulp}(hi)$. Provides about 31 decimal digits of precision (~84 bits), roughly twice the precision of a single `f64`.

**Derives**: `Debug, Clone, Copy, PartialEq`; implements `Eq`.

Arithmetic is based on the error-free transformation algorithms of Dekker and Knuth (TwoSum, TwoProd), significantly faster than arbitrary-precision alternatives such as MPFR.

#### Constants

| Constant | Value | Description |
|---|---|---|
| `DoubleF64::ZERO` | `{ hi: 0.0, lo: 0.0 }` | Zero |
| `DoubleF64::ONE` | `{ hi: 1.0, lo: 0.0 }` | One |

#### Construction and conversion

##### DoubleF64::new

**Signature**: `pub fn new(hi: f64, lo: f64) -> Self`

**Description**: Creates from high and low components. The caller must ensure $|lo| \leq 0.5 \cdot \text{ulp}(hi)$ for correctness.

##### DoubleF64::from_f64

**Signature**: `pub fn from_f64(x: f64) -> Self`

**Description**: Creates from a single `f64` ($lo = 0$).

##### DoubleF64::to_f64

**Signature**: `pub fn to_f64(self) -> f64`

**Description**: Extracts the high component.

##### From\<f64\>

**Signature**: `impl From<f64> for DoubleF64`

**Description**: Equivalent to `from_f64`.

#### Error-free transformations

##### DoubleF64::quick_two_sum

**Signature**: `pub fn quick_two_sum(a: f64, b: f64) -> Self`

**Description**: Error-free summation for $|a| \geq |b|$. Faster than `two_sum` but requires the precondition.

##### DoubleF64::two_sum

**Signature**: `pub fn two_sum(a: f64, b: f64) -> Self`

**Description**: Dekker TwoSum — error-free summation in which the rounding error is captured exactly in the `lo` component.

**Example**:

```rust
use ocas_domain::DoubleF64;

let s = DoubleF64::two_sum(1.0, f64::EPSILON);
assert_eq!(s.hi, 1.0 + f64::EPSILON);
// s.lo captures the rounding error
```

#### Queries

| Method | Signature | Description |
|---|---|---|
| `abs` | `pub fn abs(self) -> Self` | Absolute value |
| `is_nan` | `pub fn is_nan(self) -> bool` | Whether the value is NaN |
| `is_infinite` | `pub fn is_infinite(self) -> bool` | Whether the value is infinite |
| `is_finite` | `pub fn is_finite(self) -> bool` | Whether the value is finite |

#### Arithmetic operations

##### DoubleF64::add / sub / mul / div

**Signature**:

```rust
pub fn add(self, other: Self) -> Self
pub fn sub(self, other: Self) -> Self
pub fn mul(self, other: Self) -> Self
pub fn div(self, other: Self) -> Self
```

**Description**: Double-double addition, subtraction, multiplication, and division. Internally uses TwoSum/TwoProd to capture rounding errors.

##### DoubleF64::powi

**Signature**: `pub fn powi(self, n: i64) -> Self`

**Description**: Integer power using binary exponentiation. Supports negative exponents.

**Example**:

```rust
use ocas_domain::DoubleF64;

let x = DoubleF64::from_f64(3.0);
assert_eq!(x.powi(3).hi, 27.0);
assert_eq!(x.powi(-1).hi, 1.0 / 3.0);
// Output: 3^3 = 27, 3^(-1) = 0.333...
```

#### Transcendental functions

| Method | Signature | Description |
|---|---|---|
| `sqrt` | `pub fn sqrt(self) -> Self` | Square root by Newton iteration. Returns NaN for negative input |
| `exp` | `pub fn exp(self) -> Self` | Exponential function, Taylor series + argument reduction (by $\ln 2$) |
| `ln` | `pub fn ln(self) -> Self` | Natural logarithm, Newton iteration. Returns NaN for non-positive input |
| `sin` | `pub fn sin(self) -> Self` | Sine, Taylor series + reduction to $[-\pi, \pi]$ |
| `cos` | `pub fn cos(self) -> Self` | Cosine, $\cos(x) = \sin(\pi/2 - x)$ |
| `tan` | `pub fn tan(self) -> Self` | Tangent, $\sin(x)/\cos(x)$ |

**Example**:

```rust
use ocas_domain::DoubleF64;

let x = DoubleF64::from_f64(1.0);
let e = x.exp();
// e ≈ 2.718281828...
assert!((e.hi - std::f64::consts::E).abs() < 1e-15);

let pi = DoubleF64::from_f64(std::f64::consts::PI);
let s = pi.sin();
assert!(s.hi.abs() < 1e-15); // sin(π) ≈ 0
// Output: exp(1) ≈ e, sin(π) ≈ 0
```

#### Operators

`DoubleF64` implements `Add`, `Sub`, `Mul`, `Div`, `Neg`, `AddAssign`, `SubAssign`, `MulAssign`, `DivAssign`, `PartialOrd`, and `Zero` (from `num_traits`).

#### Display

Displays the `hi` value when `lo == 0`; otherwise displays the sum $hi + lo$ in 31-digit scientific notation.

**See also**: [Domain](#domain-trait), [RealBall](#realball--realballdomain)

---

## AlgebraicExtension / AlgebraicNumberField / AlgebraicElement

### AlgebraicExtension

**Signature**: `pub struct AlgebraicExtension<D: Domain> { base: D, min_poly: Vec<D::Element> }`

**Description**: The algebraic extension $D[\alpha]/(m(\alpha))$, where $m$ is a monic polynomial. When the base domain $D$ is a field and $m$ is irreducible, the quotient ring is a field:

- `AlgebraicExtension<RationalDomain>` = algebraic number field $\mathbb{Q}(\alpha)$
- `AlgebraicExtension<FiniteField>` = Galois field $\mathrm{GF}(p^d)$

Elements are residue classes represented by the unique polynomial representative of degree less than $\deg(m)$. Inverses use the extended Euclidean algorithm over the base domain.

**⚠️ Note**: irreducibility of the minimal polynomial is not checked. Under a reducible modulus the ring has zero divisors, and `Domain::inv` returns `None` for non-units.

**Implements**: `Domain`, `EuclideanDomain`

#### Construction

##### AlgebraicExtension::new

**Signature**: `pub fn new(base: D, min_poly: Vec<D::Element>) -> Self`

**Description**: Creates an algebraic extension from a base domain and a monic minimal polynomial.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `base` | `D: Domain` | The base domain |
| `min_poly` | `Vec<D::Element>` | The minimal polynomial, in ascending order, monic (leading coefficient 1), of degree $\geq 1$ |

**Errors**: debug-mode checks: (1) degree $\geq 1$ (i.e. `min_poly.len() >= 2`); (2) monic.

**Example**:

```rust
use ocas_domain::{AlgebraicExtension, Domain, Rational, RationalDomain};

// ℚ(√2): minimal polynomial α² − 2
let two = Rational::new(2, 1);
let neg_two = RationalDomain.neg(&two);
let field = AlgebraicExtension::new(
    RationalDomain,
    vec![neg_two, Rational::new(0, 1), Rational::new(1, 1)],
);
let sqrt2 = field.alpha();
// √2 · √2 = 2
assert_eq!(field.mul(&sqrt2, &sqrt2), field.from_base(two));
// Output: α² = 2
```

#### Accessors

##### AlgebraicExtension::base_domain

**Signature**: `pub fn base_domain(&self) -> &D`

**Description**: Returns a reference to the base domain.

##### AlgebraicExtension::min_poly

**Signature**: `pub fn min_poly(&self) -> &[D::Element]`

**Description**: Returns the minimal polynomial coefficients (ascending order, monic).

##### AlgebraicExtension::extension_degree

**Signature**: `pub fn extension_degree(&self) -> usize`

**Description**: Returns the extension degree $\deg(m)$.

#### Element construction

##### AlgebraicExtension::from_base

**Signature**: `pub fn from_base(&self, c: D::Element) -> AlgebraicElement<D::Element>`

**Description**: Embeds a base-domain constant into the extension.

##### AlgebraicExtension::alpha

**Signature**: `pub fn alpha(&self) -> AlgebraicElement<D::Element>`

**Description**: Returns the generator $\alpha$ of the extension.

##### AlgebraicExtension::element

**Signature**: `pub fn element(&self, coeffs: Vec<D::Element>) -> AlgebraicElement<D::Element>`

**Description**: Creates an element from coefficients (ascending order), reducing modulo the minimal polynomial automatically.

#### Domain implementation

- `mul`: multiply first, then reduce modulo $m(\alpha)$
- `inv`: extended Euclidean algorithm — if $\gcd(a, m) = 1$ (a constant), then $a^{-1} = s \bmod m$
- `div`: $a/b = a \cdot b^{-1}$
- `is_zero`: the coefficient vector is empty

#### EuclideanDomain implementation

- `div_rem`: division is exact in a field; the remainder is always zero
- `gcd`: degenerates in a field — 0 when both are zero, otherwise 1

**Example** (Galois field $\mathrm{GF}(3^2)$):

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
// The multiplicative group has order 8: (1+α)⁸ = 1
let a = field.add(&field.one(), &alpha);
assert_eq!(field.pow(&a, 8), field.one());
// Output: α² = 2, (1+α)⁸ = 1
```

**See also**: [Domain](#domain-trait), [Algebraic number fields](../math/algebraic-number-fields.md)

---

### AlgebraicNumberField

**Signature**: `pub type AlgebraicNumberField = AlgebraicExtension<RationalDomain>;`

**Description**: Type alias for $\mathbb{Q}(\alpha)$.

---

### AlgebraicElement

**Signature**: `pub struct AlgebraicElement<E> { coeffs: Vec<E> }`

**Description**: An element of an algebraic extension — a polynomial residue class in $\alpha$ of degree less than the extension degree. Coefficients are stored in ascending order with trailing zeros trimmed. The zero element has an empty coefficient vector.

**Derives**: `Debug, Clone, PartialEq, Eq, Hash`

#### AlgebraicElement::coeffs

**Signature**: `pub fn coeffs(&self) -> &[E]`

**Description**: Returns the coefficient slice (ascending order, trailing zeros trimmed).

#### Display

Display format: `0` (zero element), `c` (pure constant), `(c1)·α + c0`, `(c2)·α^2 + (c1)·α + c0`, etc. (the constant term is not parenthesized; terms are joined with ` + `).

**Example**:

```rust
use ocas_domain::{AlgebraicExtension, Domain, Rational, RationalDomain};

// ℚ(i): minimal polynomial α² + 1
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(1, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
let i = field.alpha();
let one_plus_i = field.add(&field.one(), &i);
// (1+i)⁻¹ = (1-i)/2
let inv = field.inv(&one_plus_i).unwrap();
assert_eq!(inv.coeffs(), &[Rational::new(1, 2), Rational::new(-1, 2)]);
// Output: 1/(1+i) = (1/2) - (1/2)·α
```

**See also**: [AlgebraicExtension](#algebraicextension--algebraicnumberfield--algebraicelement)

---

## Assumptions system

The assumptions system declares properties of symbolic variables (such as "x is real" or "n is a positive integer"). Solvers and simplifiers use assumptions to choose algorithms and validate solutions.

**Module path**: `ocas_domain::assumptions`

---

### Assumption

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Assumption {
    Real, Complex, Integer, Rational,
    Positive, Negative, NonNegative, NonPositive, NonZero, Finite,
    Even, Odd, Prime,
}
```

**Description**: A single predicate that can be declared about a symbolic variable. Assumptions are independent — a variable may carry several assumptions at once (e.g. `Positive | Integer`).

#### Variant descriptions

| Variant | Meaning | Display |
|---|---|---|
| `Real` | real number | `"real"` |
| `Complex` | complex number | `"complex"` |
| `Integer` | integer | `"integer"` |
| `Rational` | rational number | `"rational"` |
| `Positive` | strictly positive ($> 0$) | `"positive"` |
| `Negative` | strictly negative ($< 0$) | `"negative"` |
| `NonNegative` | non-negative ($\geq 0$) | `"non-negative"` |
| `NonPositive` | non-positive ($\leq 0$) | `"non-positive"` |
| `NonZero` | non-zero | `"non-zero"` |
| `Finite` | finite (not $\pm\infty$) | `"finite"` |
| `Even` | even | `"even"` |
| `Odd` | odd | `"odd"` |
| `Prime` | prime | `"prime"` |

#### Assumption::implied

**Signature**: `pub fn implied(&self) -> &'static [Assumption]`

**Description**: Returns the list of assumptions logically implied by this assumption. When an assumption is inserted, its implications propagate automatically.

**Implication table**:

| Assumption | Implies |
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
| `Real`, `NonZero`, `Finite` | (no further implications) |

#### Assumption::conflicts

**Signature**: `pub fn conflicts(&self) -> &'static [Assumption]`

**Description**: Returns the list of assumptions this assumption conflicts with. A set containing both an assumption and one of its conflicts is inconsistent.

**Conflict table**:

| Assumption | Conflicts with |
|---|---|
| `Positive` | `Negative`, `NonPositive` |
| `Negative` | `Positive`, `NonNegative` |
| `NonNegative` | `Negative` |
| `NonPositive` | `Positive` |
| `Even` | `Odd` |
| `Odd` | `Even` |

All other assumptions have no conflicts.

#### BitOr support

`Assumption | Assumption` returns `Assumptions`:

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let a = Assumption::Positive | Assumption::Integer;
assert!(a.implies(Assumption::Real));
// Output: Positive | Integer implies Real
```

**See also**: [Assumptions](#assumptions), [SymbolAssumptions](#symbolassumptions)

---

### Assumptions

**Signature**: `pub struct Assumptions { inner: Vec<Assumption> }`

**Description**: A set of assumptions about one symbolic variable. Stored internally as a sorted, de-duplicated vector for small-scale efficiency. Operations are closed under logical implication — inserting `Positive` makes `NonNegative` and `Real` available.

**Derives**: `Debug, Clone, PartialEq, Eq, Default`

#### Construction

##### Assumptions::new

**Signature**: `pub fn new() -> Self`

**Description**: Creates an empty assumption set.

##### Assumptions::single

**Signature**: `pub fn single(a: Assumption) -> Self`

**Description**: Creates a set containing a single assumption (and its implications).

#### Queries

##### Assumptions::len

**Signature**: `pub fn len(&self) -> usize`

**Description**: Returns the number of explicitly stored assumptions.

##### Assumptions::is_empty

**Signature**: `pub fn is_empty(&self) -> bool`

**Description**: Whether the set is empty.

##### Assumptions::contains

**Signature**: `pub fn contains(&self, a: Assumption) -> bool`

**Description**: Checks whether this set implies the assumption (as a direct member or via logical implication).

##### Assumptions::implies

**Signature**: `pub fn implies(&self, other: Assumption) -> bool`

**Description**: Checks whether this set logically implies `other`.

**Example**:

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let mut a = Assumptions::new();
a.insert(Assumption::Positive);
a.insert(Assumption::Integer);
assert!(a.contains(Assumption::Real));      // Positive implies Real
assert!(a.implies(Assumption::NonZero));    // Positive implies NonZero
assert!(!a.implies(Assumption::Even));      // does not imply Even
// Output: all assertions pass
```

#### Modification

##### Assumptions::insert

**Signature**: `pub fn insert(&mut self, a: Assumption) -> bool`

**Description**: Inserts the assumption and all of its logical implications. Returns `false` if the insertion creates a contradiction (the inconsistent assumption is still stored; check `is_consistent`).

##### Assumptions::remove

**Signature**: `pub fn remove(&mut self, a: Assumption)`

**Description**: Removes the assumption. Does not remove assumptions implied by it (they may still be implied by other stored assumptions).

##### Assumptions::is_consistent

**Signature**: `pub fn is_consistent(&self) -> bool`

**Description**: Checks whether the set is consistent (no conflicting assumptions).

##### Assumptions::iter

**Signature**: `pub fn iter(&self) -> impl Iterator<Item = Assumption> + '_`

**Description**: Iterates over the stored assumptions.

#### Operators

- `Assumption | Assumption` → `Assumptions`
- `Assumptions | Assumption` → `Assumptions`
- `Assumptions | Assumptions` → `Assumptions` (union)
- `FromIterator<Assumption>` → `Assumptions`

**Example**:

```rust
use ocas_domain::assumptions::{Assumption, Assumptions};

let a: Assumptions = [Assumption::Positive, Assumption::Integer].into_iter().collect();
assert!(a.is_consistent());

let b = Assumption::Positive | Assumption::Negative;
assert!(!b.is_consistent());
// Output: Positive + Integer is consistent; Positive + Negative is inconsistent
```

**See also**: [Assumption](#assumption), [SymbolAssumptions](#symbolassumptions)

---

### SymbolAssumptions

**Signature**: `pub struct SymbolAssumptions { entries: Vec<(String, Assumptions)> }`

**Description**: A mapping from symbol names to assumptions. Solvers and simplifiers use it to determine which transformations are legal. For example, $\sqrt{x^2} \to x$ is valid only when `x` is assumed `NonNegative`.

**Derives**: `Debug, Clone, PartialEq, Eq, Default`

#### Construction

##### SymbolAssumptions::new

**Signature**: `pub fn new() -> Self`

**Description**: Creates an empty mapping.

#### Queries and modification

##### SymbolAssumptions::set

**Signature**: `pub fn set(&mut self, symbol: &str, assumptions: Assumptions)`

**Description**: Sets the assumptions for a symbol (replacing any existing entry).

##### SymbolAssumptions::get

**Signature**: `pub fn get(&self, symbol: &str) -> Option<&Assumptions>`

**Description**: Gets the assumptions for a symbol.

##### SymbolAssumptions::remove

**Signature**: `pub fn remove(&mut self, symbol: &str)`

**Description**: Removes the assumptions for a symbol.

##### SymbolAssumptions::check

**Signature**: `pub fn check(&self, symbol: &str, assumption: Assumption) -> bool`

**Description**: Checks whether a symbol satisfies a particular assumption.

##### SymbolAssumptions::len

**Signature**: `pub fn len(&self) -> usize`

**Description**: Returns the number of symbols with assumptions.

##### SymbolAssumptions::is_empty

**Signature**: `pub fn is_empty(&self) -> bool`

**Description**: Whether no symbol has assumptions.

##### SymbolAssumptions::iter

**Signature**: `pub fn iter(&self) -> impl Iterator<Item = &(String, Assumptions)>`

**Description**: Iterates over all `(symbol name, assumption set)` pairs.

**Example**:

```rust
use ocas_domain::assumptions::{Assumption, Assumptions, SymbolAssumptions};

let mut sa = SymbolAssumptions::new();
sa.set("x", Assumptions::single(Assumption::Positive));
sa.set("n", Assumption::Integer | Assumption::Positive);

assert!(sa.check("x", Assumption::Real));      // Positive implies Real
assert!(sa.check("n", Assumption::NonZero));   // Positive implies NonZero
assert!(!sa.check("x", Assumption::Integer));  // does not imply Integer
// Output: all assertions pass
```

**See also**: [Assumptions](#assumptions), [Assumption](#assumption)
