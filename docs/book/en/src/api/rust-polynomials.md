# Rust API Reference: Polynomials

This chapter documents oCAS's polynomial system, covering three core data structures:

- **`DenseUnivariatePolynomial<D>`** — dense univariate polynomials
- **`SparseMultivariatePolynomial<D, O>`** — sparse multivariate polynomials
- **`RationalPolynomial<D, O>`** — elements of the polynomial fraction field

as well as the monomial order trait `MonomialOrder` and its implementations for various orders.

**Module path**: `ocas_poly`

---

## Table of Contents

- [Monomial Orders](#monomial-orders)
  - [MonomialOrder trait](#monomialorder-trait)
  - [Lex](#lex)
  - [Grlex](#grlex)
  - [Grevlex](#grevlex)
  - [WeightOrder](#weightorder)
  - [BlockOrder and SubOrder](#blockorder-and-suborder)
  - [MatrixOrder](#matrixorder)
- [DenseUnivariatePolynomial](#denseunivariatepolynomial)
  - [Construction and Properties](#dense-construction-and-properties)
  - [Basic Arithmetic](#dense-basic-arithmetic)
  - [EuclideanDomain Operations](#euclideandomain-operations)
  - [Factorization and Resultants](#factorization-and-resultants)
- [SparseMultivariatePolynomial](#sparsemultivariatepolynomial)
  - [Construction and Properties](#sparse-construction-and-properties)
  - [Basic Arithmetic](#sparse-basic-arithmetic)
  - [Gröbner Basis Support](#gröbner-basis-support)
  - [Multivariate Factorization](#multivariate-factorization)
- [RationalPolynomial](#rationalpolynomial)
- [Helper Functions](#helper-functions)

---

## Monomial Orders

### MonomialOrder trait

**Signature**:

```rust
pub trait MonomialOrder: Clone + PartialEq + Eq + std::fmt::Debug + Default {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering;
}
```

**Description**: Defines a total order on monomials. A polynomial's leading term, sorting, and Gröbner basis computations all depend on it.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `self` | `&Self` | The order instance (a zero-sized type or a type with parameters) |
| `lhs` | `&[usize]` | Exponent vector of the left monomial $[\alpha_1, \alpha_2, \dots]$ |
| `rhs` | `&[usize]` | Exponent vector of the right monomial |

**Returns**: `std::cmp::Ordering` — `Less` means `lhs` comes **before** `rhs` in the order (i.e., `lhs` is larger), `Greater` means `lhs` is smaller.

**Design notes**: The simple orders (Lex, Grevlex, Grlex) are zero-sized types with no runtime overhead. The parameterized orders (WeightOrder, BlockOrder, MatrixOrder) store their configuration at construction time.

**See also**: [Lex](#lex), [Grevlex](#grevlex), [Grlex](#grlex), [WeightOrder](#weightorder), [BlockOrder](#blockorder-and-suborder), [MatrixOrder](#matrixorder)

---

### Lex

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Lex;
```

**Description**: Lexicographic order. Compares exponent vectors component by component from left to right.

**Comparison rule**: `lhs > rhs` if and only if there exists an $i$ such that $\alpha_j = \beta_j$ for all $j < i$ and $\alpha_i > \beta_i$.

**Example**:

```rust
use ocas_poly::sparse::{Lex, MonomialOrder};

let a = [2, 1]; // x^2 y
let b = [1, 1]; // x y
assert_eq!(Lex.cmp(&a, &b), std::cmp::Ordering::Greater);
// a is larger in Lex order (first component 2 > 1)
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

### Grlex

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grlex;
```

**Description**: Graded lexicographic order. Orders by total degree descending first, then lexicographically among monomials of equal degree.

**Comparison rule**:
1. The monomial with the larger total degree $\sum \alpha_i$ comes first
2. Monomials of equal total degree are compared lexicographically

**Example**:

```rust
use ocas_poly::sparse::{Grlex, MonomialOrder};

let a = [2, 0]; // x^2, degree 2
let b = [1, 1]; // xy, degree 2
let c = [0, 3]; // y^3, degree 3
// c has the highest degree, so it comes first
assert_eq!(Grlex.cmp(&c, &a), std::cmp::Ordering::Less);
// a and b have the same degree; lexicographically a > b
assert_eq!(Grlex.cmp(&a, &b), std::cmp::Ordering::Greater);
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

### Grevlex

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grevlex;
```

**Description**: Graded reverse lexicographic order. Orders by total degree descending first, then compares in **reverse** lexicographic order (starting from the last component, with the direction reversed).

**Comparison rule**:
1. The monomial with the larger total degree comes first
2. For equal total degree, compare starting from the last component, where the **smaller** component comes first

**Properties**: Grevlex typically produces the smallest intermediate matrices in Gröbner basis computations and is the default order.

**Example**:

```rust
use ocas_poly::sparse::{Grevlex, Lex, MonomialOrder};

let a = [2, 1];
let b = [1, 1];
assert_eq!(Lex.cmp(&a, &b), std::cmp::Ordering::Greater);
assert_eq!(Grevlex.cmp(&a, &b), std::cmp::Ordering::Less);
// Under Grevlex, a < b (reverse lexicographic order flips the direction)
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

### WeightOrder

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightOrder {
    weights: SmallVec<[i64; 4]>,
}
```

**Description**: Weighted order. Orders by $\sum_i w_i \cdot \alpha_i$ descending. Suitable for elimination orders that cannot be expressed with a zero-sized type.

**Construction**:

```rust
pub fn new(weights: SmallVec<[i64; 4]>) -> Self
pub fn from_slice(weights: &[i64]) -> Self
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `weights` | `SmallVec<[i64; 4]>` or `&[i64]` | Weight of each variable; the length should equal the number of variables |

**Default**: all weights set to 1 (i.e., the total degree order).

**Example**:

```rust
use ocas_poly::sparse::{MonomialOrder, WeightOrder};
use smallvec::smallvec;

let ord = WeightOrder::new(smallvec![2, 1]);
// [1,0] → weight 2, [0,1] → weight 1 → [1,0] is larger
assert_eq!(ord.cmp(&[1, 0], &[0, 1]), std::cmp::Ordering::Less);
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

### BlockOrder and SubOrder

**Signature**:

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

**Description**: Block elimination order. Partitions the variables into consecutive blocks, each block using its own sub-order.

**Construction**:

```rust
pub fn new(boundaries: SmallVec<[usize; 4]>, orders: SmallVec<[SubOrder; 4]>) -> Self
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `boundaries` | `SmallVec<[usize; 4]>` | List of split points for the ordering (exclusive upper bounds, not including `n_vars`) |
| `orders` | `SmallVec<[SubOrder; 4]>` | Sub-order for each block; its length must equal `boundaries.len() + 1` |

**Notes**: `boundaries = [2]` with `orders = [Lex, Grevlex]` on a polynomial in 4 variables means: compare variables 0–1 first (Lex), and if equal, compare variables 2–3 (Grevlex).

**Default**: a single Grevlex block.

**Example**:

```rust
use ocas_poly::sparse::{BlockOrder, MonomialOrder, SubOrder};
use smallvec::smallvec;

let ord = BlockOrder::new(
    smallvec![2],
    smallvec![SubOrder::Lex, SubOrder::Grevlex],
);
let a = [1, 0, 0, 0]; // x₀
let b = [0, 1, 0, 0]; // x₁
// Lex within the block: [1,0] > [0,1]
assert_eq!(ord.cmp(&a, &b), std::cmp::Ordering::Greater);
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

### MatrixOrder

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixOrder {
    matrix: Vec<Vec<i64>>,
    n_vars: usize,
}
```

**Description**: Matrix order. Multiplies exponent vectors by an integer matrix and compares the results lexicographically. Given an $n \times n$ matrix $M$, a monomial $\alpha > \beta$ if and only if $M\alpha >_{\text{lex}} M\beta$. It generalizes all standard orders and is especially suitable for constructing elimination orders.

**Construction**:

```rust
pub fn new(matrix: Vec<Vec<i64>>) -> Self
pub fn elimination_order(elim_vars: usize, n_vars: usize) -> Self
```

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `matrix` | `Vec<Vec<i64>>` | $n \times n$ weight matrix (row-major) |
| `elim_vars` | `usize` | Number of leading variables to eliminate |
| `n_vars` | `usize` | Total number of variables |

**`elimination_order`**: Constructs a matrix order equivalent to `BlockOrder([elim_vars in Lex, rest in Grevlex])`.

**Default**: the $1 \times 1$ identity matrix.

**Example**:

```rust
use ocas_poly::sparse::{MatrixOrder, MonomialOrder};

// 2×2 identity matrix = Lex order
let ord = MatrixOrder::new(vec![vec![1, 0], vec![0, 1]]);
assert_eq!(ord.cmp(&[1, 0], &[0, 1]), std::cmp::Ordering::Greater);
```

**See also**: [MonomialOrder trait](#monomialorder-trait)

---

## DenseUnivariatePolynomial

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseUnivariatePolynomial<D: Domain> {
    coeffs: Vec<D::Element>,  // starts at the constant term, trailing zeros removed
    domain: D,
}
```

A dense univariate polynomial. Coefficients are stored in a contiguous vector, from the constant term $a_0$ to the highest-degree term $a_n$. The zero polynomial is represented by an empty vector. Multiplication automatically selects between Karatsuba (when both polynomials have at least 32 coefficients) and schoolbook multiplication.

---

### Dense Construction and Properties

#### DenseUnivariatePolynomial::new

**Signature**: `pub fn new(domain: D) -> Self`

**Description**: Creates the zero polynomial.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `domain` | `D` | The coefficient domain |

**Returns**: the zero polynomial (an empty coefficient vector).

**See also**: [from_coeffs](#denseunivariatepolynomialfrom_coeffs)

---

#### DenseUnivariatePolynomial::from_coeffs

**Signature**: `pub fn from_coeffs(domain: D, coeffs: Vec<D::Element>) -> Self`

**Description**: Constructs a polynomial from a coefficient vector. Trailing zero coefficients are removed automatically.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `domain` | `D` | The coefficient domain |
| `coeffs` | `Vec<D::Element>` | Coefficients $[a_0, a_1, \dots, a_n]$, constant term first |

**Returns**: the polynomial $a_0 + a_1 x + \cdots + a_n x^n$.

**Example**:

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
// p(x) = 1 + 2x^2, the middle zero coefficient is kept (not a trailing zero)
```

**See also**: [new](#denseunivariatepolynomialnew)

---

#### DenseUnivariatePolynomial::domain

**Signature**: `pub fn domain(&self) -> &D`

**Description**: Returns a reference to the coefficient domain.

---

#### DenseUnivariatePolynomial::coeffs

**Signature**: `pub fn coeffs(&self) -> &[D::Element]`

**Description**: Returns the coefficient slice, starting from the constant term.

---

#### DenseUnivariatePolynomial::is_zero

**Signature**: `pub fn is_zero(&self) -> bool`

**Description**: Returns whether the polynomial is zero.

---

#### DenseUnivariatePolynomial::degree

**Signature**: `pub fn degree(&self) -> Option<usize>`

**Description**: Returns the degree of the polynomial. Returns `None` for the zero polynomial.

---

#### DenseUnivariatePolynomial::coeff

**Signature**: `pub fn coeff(&self, n: usize) -> Option<&D::Element>`

**Description**: Returns the coefficient of $x^n$, or `None` if it does not exist.

---

#### DenseUnivariatePolynomial::leading_coeff

**Signature**: `pub fn leading_coeff(&self) -> Option<&D::Element>`

**Description**: Returns the leading coefficient. Returns `None` for the zero polynomial.

---

#### DenseUnivariatePolynomial::lcoeff

**Signature**: `pub fn lcoeff(&self) -> D::Element`

**Description**: Convenience alias for the leading coefficient. Returns the domain's zero for the zero polynomial.

**Returns**: the leading coefficient, or the domain's zero.

**See also**: [leading_coeff](#denseunivariatepolynomialleading_coeff)

---

#### DenseUnivariatePolynomial::constant

**Signature**: `pub fn constant(&self) -> D::Element`

**Description**: Returns the constant term (the coefficient of $x^0$). Returns the domain's zero for the zero polynomial.

---

#### DenseUnivariatePolynomial::zero

**Signature**: `pub fn zero(&self) -> Self`

**Description**: Returns the zero polynomial over the same domain.

---

#### DenseUnivariatePolynomial::one

**Signature**: `pub fn one(&self) -> Self`

**Description**: Returns the constant polynomial $1$.

---

#### DenseUnivariatePolynomial::is_one

**Signature**: `pub fn is_one(&self) -> bool`

**Description**: Returns whether the polynomial is the constant $1$.

---

### Dense Basic Arithmetic

The following methods are available for `D: Domain`:

#### DenseUnivariatePolynomial::neg

**Signature**: `pub fn neg(&self) -> Self`

**Description**: Returns $-p(x)$.

---

#### DenseUnivariatePolynomial::add

**Signature**: `pub fn add(&self, other: &Self) -> Self`

**Description**: Polynomial addition.

---

#### DenseUnivariatePolynomial::sub

**Signature**: `pub fn sub(&self, other: &Self) -> Self`

**Description**: Polynomial subtraction.

---

#### DenseUnivariatePolynomial::mul_scalar

**Signature**: `pub fn mul_scalar(&self, scalar: &D::Element) -> Self`

**Description**: Scalar multiplication. Multiplies every coefficient by `scalar`.

---

#### DenseUnivariatePolynomial::mul

**Signature**: `pub fn mul(&self, other: &Self) -> Self`

**Description**: Polynomial multiplication. Automatically selects between schoolbook multiplication and Karatsuba based on the degree.

**Example**:

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

**See also**: [mul_into](#denseunivariatepolynomialmul_into)

---

#### DenseUnivariatePolynomial::mul_into

**Signature**: `pub fn mul_into(&self, other: &Self, buf: &mut Vec<D::Element>)`

**Description**: Multiplication that writes the result into a buffer (avoiding repeated heap allocations in hot loops). The buffer is cleared and reused.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Self` | The other polynomial |
| `buf` | `&mut Vec<D::Element>` | Output buffer; after the call it contains the product coefficients |

---

#### DenseUnivariatePolynomial::eval

**Signature**: `pub fn eval(&self, x: &D::Element) -> D::Element`

**Description**: Evaluates the polynomial at $x$ using Horner's method. Returns the domain's zero for the zero polynomial.

**Example**:

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

**Signature**: `pub fn derivative(&self) -> Self`

**Description**: Returns the formal derivative $p'(x)$.

---

#### DenseUnivariatePolynomial::integral

**Signature**: `pub fn integral(&self) -> Self`

**Description**: Returns the formal integral $\int p(x)\,dx$ with zero constant term. Requires a domain that supports division (such as `RationalDomain`).

---

### EuclideanDomain Operations

The following methods require `D: EuclideanDomain`:

#### DenseUnivariatePolynomial::mul_coeff

**Signature**: `pub fn mul_coeff(&self, c: &D::Element) -> Self`

**Description**: Multiplies every coefficient by the constant $c$. Equivalent to `mul_scalar`, but restricted to `EuclideanDomain`.

**See also**: [mul_scalar](#denseunivariatepolynomialmul_scalar)

---

#### DenseUnivariatePolynomial::div_coeff

**Signature**: `pub fn div_coeff(&self, c: &D::Element) -> Self`

**Description**: Divides every coefficient by the constant $c$ (which must divide exactly).

**Errors**: panics unconditionally if $c$ is not invertible in the domain (e.g. over $\mathbb{Z}$ when $c$ is not $\pm 1$), even if every coefficient is divisible by $c$.

---

#### DenseUnivariatePolynomial::div_rem

**Signature**: `pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>`

**Description**: Division with remainder. Returns $(q, r)$ such that $p = q \cdot \text{divisor} + r$, where $\deg(r) < \deg(\text{divisor})$.

**Returns**: `Some((quotient, remainder))`, or `None` if the divisor is zero.

**Example**:

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
// (x^2 - 1) / (x + 1) = x - 1, remainder 0
```

**See also**: [gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::gcd

**Signature**: `pub fn gcd(&self, other: &Self) -> Self`

**Description**: Computes the greatest common divisor of two polynomials. Uses the Euclidean algorithm (with pseudo-remainders when the coefficient domain is not a field). The result is always primitive (content-free).

**Example**:

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

**See also**: [extended_gcd_poly](#denseunivariatepolynomialextended_gcd_poly)

---

#### DenseUnivariatePolynomial::content

**Signature**: `pub fn content(&self) -> D::Element`

**Description**: Returns the greatest common divisor of all coefficients. Returns the domain's zero for the zero polynomial.

**Returns**: the GCD of the coefficients.

**See also**: [primitive_part](#denseunivariatepolynomialprimitive_part)

---

#### DenseUnivariatePolynomial::primitive_part

**Signature**: `pub fn primitive_part(&self) -> Self`

**Description**: Returns the primitive part $p / \text{content}(p)$. The result has content 1 (or is the zero polynomial).

**See also**: [content](#denseunivariatepolynomialcontent)

---

#### DenseUnivariatePolynomial::extended_gcd_poly

**Signature**: `pub fn extended_gcd_poly(&self, other: &Self) -> (Self, Self, Self)`

**Description**: Extended Euclidean algorithm. Returns $(g, s, t)$ such that $s \cdot p + t \cdot q = g$, where $g = \gcd(p, q)$ and $g$ is monic.

**Returns**: `(gcd, bezout_s, bezout_t)`

**Example**:

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

**See also**: [gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::pow

**Signature**: `pub fn pow(&self, n: u32) -> Self`

**Description**: Fast exponentiation by repeated squaring. $p^0 = 1$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `u32` | Non-negative integer exponent |

**Returns**: $p(x)^n$.

---

#### DenseUnivariatePolynomial::p_adic_expansion

**Signature**: `pub fn p_adic_expansion(&self, p: &Self) -> Vec<Self>`

**Description**: $p$-adic expansion. Returns $[a_0, a_1, a_2, \dots]$ such that $\text{self} = a_0 + a_1 \cdot p + a_2 \cdot p^2 + \cdots$, where each $a_k$ has degree less than $\deg(p)$.

**Implementation**: by repeated division with remainder.

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// p-adic expansion of f(x) = x^3 with respect to p(x) = x + 1
let f = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(0), Integer::from(0), Integer::from(0), Integer::from(1),
]);
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(1),
]);
let expansion = f.p_adic_expansion(&p);
// expansion = [a0, a1, a2, ...] such that f = a0 + a1*p + a2*p^2 + ...
```

**See also**: [div_rem](#denseunivariatepolynomialdiv_rem)

---

#### DenseUnivariatePolynomial::diophantine

**Signature**: `pub fn diophantine(polys: &mut [Self], b: &Self) -> Vec<Self>`

**Description**: Polynomial CRT (Diophantine solver). Given a list of pairwise coprime polynomials `polys` and a target $b$, returns $[s_0, \dots, s_n]$ such that:

$$\sum_i s_i \cdot \prod_{j \neq i} p_j \equiv b \pmod{\prod_i p_i}$$

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `polys` | `&mut [Self]` | List of pairwise coprime polynomials |
| `b` | `&Self` | The target polynomial |

**Returns**: the list of remainders $[s_0, \dots, s_n]$.

**Errors**: panics if the polynomials are not pairwise coprime.

---

### Factorization and Resultants

#### DenseUnivariatePolynomial::square_free_factorization

**Signature**: `pub fn square_free_factorization(&self) -> SquareFreeFactors<D>`

**Description**: Square-free factorization. Uses Yun's algorithm: $g = \gcd(f, f')$, $w = f/g$, then iterate $h = \gcd(w, g)$, $z = w/h$ to collect the factors of multiplicity $k$.

**Returns**: `Vec<(DenseUnivariatePolynomial<D>, usize)>` — a list of `(factor, multiplicity)` pairs.

**Example**:

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

**See also**: [factor](#denseunivariatepolynomialfactor), [gcd](#denseunivariatepolynomialgcd)

---

#### DenseUnivariatePolynomial::is_square_free

**Signature**: `pub fn is_square_free(&self) -> bool`

**Description**: Returns whether the polynomial is square-free ($\gcd(p, p') = 1$).

---

#### DenseUnivariatePolynomial::factor (IntegerDomain)

**Signature**: `impl DenseUnivariatePolynomial<IntegerDomain> { pub fn factor(&self) -> Factors<IntegerDomain> }`

**Description**: Completely factors a primitive integer polynomial into irreducible factors. Uses the square-free factorization + Berlekamp–Zassenhaus + Hensel lifting algorithms.

**Preconditions**: the input must be primitive (content = 1). Call `primitive_part()` first on arbitrary polynomials.

**Returns**: `Vec<(DenseUnivariatePolynomial<IntegerDomain>, usize)>` — `(irreducible factor, multiplicity)` pairs.

**Example**:

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

**See also**: [square_free_factorization](#denseunivariatepolynomialsquare_free_factorization)

---

#### DenseUnivariatePolynomial::factor (FiniteField)

**Signature**: `impl DenseUnivariatePolynomial<FiniteField> { pub fn factor(&self) -> Factors<FiniteField> }`

**Description**: Completely factors a polynomial over $\mathbb{F}_p$. Uses Berlekamp's algorithm (or Cantor–Zassenhaus).

**Returns**: `Vec<(DenseUnivariatePolynomial<FiniteField>, usize)>` — `(irreducible factor, multiplicity)` pairs.

**Example**:

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

**See also**: [factor (IntegerDomain)](#denseunivariatepolynomialfactor-integerdomain)

---

#### DenseUnivariatePolynomial::resultant

**Signature**: `pub fn resultant(&self, other: &Self) -> D::Element`

**Description**: Computes the resultant $\operatorname{Res}(a, b)$ using the Brown PRS algorithm. The resultant is zero if and only if $a$ and $b$ have a non-constant common factor.

**Implementation details**: subresultant PRS, with an exact division by $\beta$ at each step (the subresultant theorem guarantees exactness).

**Example**:

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

**Properties**: $\operatorname{Res}(a, b) = (-1)^{\deg a \cdot \deg b} \operatorname{Res}(b, a)$.

**See also**: [gcd](#denseunivariatepolynomialgcd)

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

A sparse multivariate polynomial. Only nonzero terms are stored, in a HashMap keyed by the exponent vector $\vec{e} = [e_1, e_2, \dots]$ with the coefficient as value. An exponent vector $\vec{e}$ represents the monomial $x_1^{e_1} x_2^{e_2} \cdots$. The monomial order is controlled by the type parameter `O`, which defaults to `Grevlex`.

---

### Sparse Construction and Properties

#### SparseMultivariatePolynomial::new

**Signature**: `pub fn new(domain: D, n_vars: usize) -> Self`

**Description**: Creates the zero polynomial with the default monomial order.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `domain` | `D` | The coefficient domain |
| `n_vars` | `usize` | Number of variables |

**See also**: [new_with_order](#sparsemultivariatepolynomialnew_with_order)

---

#### SparseMultivariatePolynomial::new_with_order

**Signature**: `pub fn new_with_order(domain: D, n_vars: usize, order: O) -> Self`

**Description**: Creates the zero polynomial with the given monomial order.

**Example**:

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

**Signature**: `pub fn from_terms(domain: D, n_vars: usize, terms: Vec<(Vec<usize>, D::Element)>) -> Self`

**Description**: Constructs a polynomial from a list of (exponent vector, coefficient) pairs. Terms with zero coefficients are dropped automatically.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `domain` | `D` | The coefficient domain |
| `n_vars` | `usize` | Number of variables |
| `terms` | `Vec<(Vec<usize>, D::Element)>` | (exponent vector, coefficient) pairs |

**Example**:

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

**See also**: [new](#sparsemultivariatepolynomialnew)

---

#### SparseMultivariatePolynomial::domain

**Signature**: `pub fn domain(&self) -> &D`

**Description**: Returns a reference to the coefficient domain.

---

#### SparseMultivariatePolynomial::n_vars

**Signature**: `pub fn n_vars(&self) -> usize`

**Description**: Returns the number of variables.

---

#### SparseMultivariatePolynomial::n_terms

**Signature**: `pub fn n_terms(&self) -> usize`

**Description**: Returns the number of nonzero terms.

---

#### SparseMultivariatePolynomial::is_zero

**Signature**: `pub fn is_zero(&self) -> bool`

**Description**: Returns whether the polynomial is zero.

---

#### SparseMultivariatePolynomial::terms_ref

**Signature**: `pub fn terms_ref(&self) -> &HashMap<SmallVec<[usize; 4]>, D::Element>`

**Description**: Returns a reference to the internal term map (exponent → coefficient).

---

#### SparseMultivariatePolynomial::set_term_external

**Signature**: `pub fn set_term_external(&mut self, exp: Vec<usize>, coeff: D::Element)`

**Description**: Sets the coefficient of a monomial. A zero coefficient removes the term.

---

#### SparseMultivariatePolynomial::total_degree

**Signature**: `pub fn total_degree(&self) -> Option<usize>`

**Description**: Returns the total degree (the maximum degree over all monomials). Returns `None` for the zero polynomial.

---

#### SparseMultivariatePolynomial::coeff

**Signature**: `pub fn coeff(&self, exp: &[usize]) -> D::Element`

**Description**: Returns the coefficient of the given monomial, or the domain's zero if it is absent.

---

#### SparseMultivariatePolynomial::degree_in

**Signature**: `pub fn degree_in(&self, var_index: usize) -> usize`

**Description**: Returns the degree in the variable `var_index`. Returns 0 for the zero polynomial.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var_index` | `usize` | Variable index (0-based) |

**Returns**: the maximum exponent of that variable over all terms.

---

#### SparseMultivariatePolynomial::zero / one

**Signature**:

```rust
pub fn zero(&self) -> Self
pub fn one(&self) -> Self
```

**Description**: Returns the zero polynomial or the constant $1$ with the same shape.

---

### Sparse Basic Arithmetic

#### SparseMultivariatePolynomial::neg

**Signature**: `pub fn neg(&self) -> Self`

**Description**: Returns $-p$.

---

#### SparseMultivariatePolynomial::add

**Signature**: `pub fn add(&self, other: &Self) -> Self`

**Description**: Polynomial addition.

**Errors**: panics if the numbers of variables differ.

---

#### SparseMultivariatePolynomial::sub

**Signature**: `pub fn sub(&self, other: &Self) -> Self`

**Description**: Polynomial subtraction.

**Errors**: panics if the numbers of variables differ.

---

#### SparseMultivariatePolynomial::mul_scalar

**Signature**: `pub fn mul_scalar(&self, scalar: &D::Element) -> Self`

**Description**: Scalar multiplication.

---

#### SparseMultivariatePolynomial::mul

**Signature**: `pub fn mul(&self, other: &Self) -> Self`

**Description**: Polynomial multiplication. Multiplies each pair of terms and combines like terms.

**Errors**: panics if the numbers of variables differ.

**Example**:

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

**Signature**: `pub fn mul_monomial(&self, exp: &[usize]) -> Self`

**Description**: Adds `exp` component-wise to the exponent vector of every term. Used in Gröbner basis reduction.

---

#### SparseMultivariatePolynomial::sorted_terms

**Signature**: `pub fn sorted_terms(&self) -> Vec<(&SmallVec<[usize; 4]>, &D::Element)>`

**Description**: Returns the (exponent, coefficient) pairs sorted by the monomial order.

---

#### SparseMultivariatePolynomial::eval

**Signature**: `pub fn eval(&self, var_index: usize, value: &D::Element) -> Self`

**Description**: Substitutes the value `value` for the variable `var_index`, returning a polynomial in one fewer variable.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `var_index` | `usize` | Index of the variable to substitute |
| `value` | `&D::Element` | The substitution value |

**Returns**: a polynomial in one fewer variable (the remaining variables keep their relative order).

**Example**:

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
// Substitute x=3: the result is 3y + 2y = 5y
let r = p.eval(0, &Integer::from(3));
assert_eq!(r.coeff(&[1]), Integer::from(5));
```

---

#### SparseMultivariatePolynomial::eval_keep

**Signature**: `pub fn eval_keep(&self, var_index: usize, value: &D::Element) -> Self`

**Description**: Substitutes the variable `var_index` with a value but keeps the total number of variables unchanged (the substituted variable's exponents are zeroed). Used in scenarios where variable positions must stay fixed during Hensel lifting.

---

#### SparseMultivariatePolynomial::leading_term

**Signature**: `pub fn leading_term(&self) -> Option<(&SmallVec<[usize; 4]>, &D::Element)>`

**Description**: Returns the leading term (exponent vector, coefficient). Returns `None` for the zero polynomial. An $O(n)$ scan of the HashMap.

---

#### SparseMultivariatePolynomial::leading_monomial

**Signature**: `pub fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>>`

**Description**: Returns the exponent vector of the leading monomial.

---

#### SparseMultivariatePolynomial::leading_coeff

**Signature**: `pub fn leading_coeff(&self) -> Option<&D::Element>`

**Description**: Returns the leading coefficient.

---

#### SparseMultivariatePolynomial::content

**Signature**:

```rust
pub fn content(&self) -> D::Element
where D: EuclideanDomain
```

**Description**: Returns the GCD of all coefficients.

**Example**:

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

**See also**: [primitive_part](#sparsemultivariatepolynomialprimitive_part)

---

#### SparseMultivariatePolynomial::primitive_part

**Signature**:

```rust
pub fn primitive_part(&self) -> Self
where D: EuclideanDomain
```

**Description**: Returns the primitive part (polynomial / content).

**See also**: [content](#sparsemultivariatepolynomialcontent)

---

#### SparseMultivariatePolynomial::div_exact

**Signature**: `pub fn div_exact(&self, divisor: &Self) -> Self`

**Description**: Exact division (assumes no remainder). Used in rational function normalization where the GCD is known to divide.

**Errors**: panics in debug mode if the division is not exact.

**See also**: [checked_div_exact](#sparsemultivariatepolynomialchecked_div_exact)

---

#### SparseMultivariatePolynomial::checked_div_exact

**Signature**: `pub fn checked_div_exact(&self, divisor: &Self) -> Option<Self>`

**Description**: Exact division, returning `None` when the division is not exact.

**Returns**: `Some(quotient)` or `None`.

---

#### SparseMultivariatePolynomial::derivative

**Signature**: `pub fn derivative(&self, var_index: usize) -> Self`

**Description**: Computes the formal partial derivative with respect to the variable `var_index`.

---

#### SparseMultivariatePolynomial::taylor_coefficients

**Signature**: `pub fn taylor_coefficients(&self, var_index: usize, a: &D::Element) -> Vec<Self>`

**Description**: Computes the Taylor coefficients in the variable `var_index` at the point $a$. Returns $[t_0, t_1, \dots, t_d]$ such that $f = \sum_j t_j (x_{\text{var}} - a)^j$.

---

### Gröbner Basis Support

#### SparseMultivariatePolynomial::reduce

**Signature**: `pub fn reduce(&self, basis: &[Self]) -> Self`

**Description**: Multivariate polynomial division. Repeatedly finds a basis element whose leading term divides the current leading term and subtracts the appropriate multiple; otherwise the leading term is moved into the remainder. Requires the coefficient domain to be a field (division always succeeds).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `basis` | `&[Self]` | The list of divisors |

**Returns**: the reduced remainder.

---

#### SparseMultivariatePolynomial::spoly

**Signature**: `pub fn spoly(&self, other: &Self) -> Self`

**Description**: Computes the S-polynomial. $S(f, g) = \frac{\text{lcm}}{\text{lt}(f)} \cdot f - \frac{\text{lcm}}{\text{lt}(g)} \cdot g$.

---

#### SparseMultivariatePolynomial::make_monic_inplace

**Signature**: `pub fn make_monic_inplace(&mut self) -> bool`

**Description**: Normalizes the polynomial to be monic in place (divides by the leading coefficient). Returns `false` if the leading coefficient is not invertible.

---

#### SparseMultivariatePolynomial::exponents_iter

**Signature**: `pub fn exponents_iter(&self) -> impl Iterator<Item = &SmallVec<[usize; 4]>>`

**Description**: Iterates over all exponent vectors in monomial order (descending). Used for the symbolic preprocessing of the F4 algorithm.

---

### Multivariate Factorization

#### SparseMultivariatePolynomial::factor (IntegerDomain, Lex)

**Signature**:

```rust
impl SparseMultivariatePolynomial<IntegerDomain, Lex> {
    pub fn factor(&self) -> Vec<(Self, usize)>
}
```

**Description**: Factors a multivariate integer polynomial into irreducible factors. Bivariate polynomials with a constant leading coefficient in the main variable use bivariate Hensel lifting (Wang); polynomials in three or more variables, or bivariate ones with a non-constant leading coefficient, use the EEZ algorithm.

**Returns**: a list of `(irreducible factor, multiplicity)` pairs.

**Example**:

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

**See also**: [factor (FiniteField)](#sparsemultivariatepolynomialfactor-finitefield-lex)

---

#### SparseMultivariatePolynomial::factor (FiniteField, Lex)

**Signature**:

```rust
impl SparseMultivariatePolynomial<FiniteField, Lex> {
    pub fn factor(&self) -> Vec<(Self, usize)>
}
```

**Description**: Factors a multivariate $\mathbb{F}_p$ polynomial. Bivariate polynomials use the evaluation–Hensel path; three or more variables use EEZ.

**See also**: [factor (IntegerDomain)](#sparsemultivariatepolynomialfactor-integerdomain-lex)

---

## RationalPolynomial

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalPolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    pub numerator: SparseMultivariatePolynomial<D, O>,
    pub denominator: SparseMultivariatePolynomial<D, O>,
}
```

An element of the polynomial fraction field $\frac{\text{num}}{\text{den}}$. Constructing via `from_num_den` automatically normalizes to lowest terms (coprime numerator and denominator, with the denominator's leading coefficient normalized).

---

#### RationalPolynomial::new

**Signature**:

```rust
pub fn new(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

**Description**: Constructs a rational polynomial without normalization. The caller must ensure the denominator is nonzero.

**See also**: [from_num_den](#rationalpolynomialfrom_num_den)

---

#### RationalPolynomial::from_num_den

**Signature**:

```rust
impl<D: EuclideanDomain, O: MonomialOrder> RationalPolynomial<D, O> {
    pub fn from_num_den(
        numerator: SparseMultivariatePolynomial<D, O>,
        denominator: SparseMultivariatePolynomial<D, O>,
    ) -> Self
}
```

**Description**: Constructs from numerator and denominator and normalizes. The result is in lowest terms: coprime numerator and denominator, with the denominator's leading coefficient normalized (1 over a finite field, positive over the integers).

**Preconditions**: the denominator must be nonzero.

**Normalization steps**:
1. Compute $\gcd(\text{num}, \text{den})$ and cancel it
2. In the univariate case, reduce by an exact GCD via the dense Euclidean algorithm
3. Normalize the denominator's leading coefficient

**Errors**: panics when the denominator is zero.

**See also**: [new](#rationalpolynomialnew)

---

#### RationalPolynomial::from_polynomial

**Signature**: `pub fn from_polynomial(poly: SparseMultivariatePolynomial<D, O>) -> Self`

**Description**: Constructs from a polynomial (denominator = 1).

---

#### RationalPolynomial::zero / one

**Signature**:

```rust
pub fn zero(domain: &D, n_vars: usize) -> Self
pub fn one(domain: &D, n_vars: usize) -> Self
```

**Description**: Returns the zero or the identity rational polynomial.

---

#### RationalPolynomial::is_zero / is_one

**Signature**:

```rust
pub fn is_zero(&self) -> bool
pub fn is_one(&self) -> bool
```

**Description**: Returns whether the rational polynomial is zero or $1/1$.

---

#### RationalPolynomial::n_vars

**Signature**: `pub fn n_vars(&self) -> usize`

**Description**: Returns the number of variables.

---

#### RationalPolynomial::domain

**Signature**: `pub fn domain(&self) -> &D`

**Description**: Returns a reference to the coefficient domain.

---

#### RationalPolynomial::neg

**Signature**: `pub fn neg(&self) -> Self`

**Description**: Returns $-\frac{n}{d}$.

---

#### RationalPolynomial::inv

**Signature**: `pub fn inv(&self) -> Option<Self>`

**Description**: Returns the multiplicative inverse $\frac{d}{n}$. Returns `None` when the numerator is zero.

**Returns**: `Some(inverse)` or `None` (when the numerator is zero).

---

#### RationalPolynomial::pow

**Signature**: `pub fn pow(&self, k: u32) -> Self`

**Description**: Fast exponentiation. Computes $\left(\frac{n}{d}\right)^k$ by repeated squaring of the numerator and denominator separately.

---

#### RationalPolynomial::add

**Signature**:

```rust
impl<D: EuclideanDomain, O: MonomialOrder> RationalPolynomial<D, O> {
    pub fn add(&self, other: &Self) -> Self
}
```

**Description**: Addition of rational polynomials. With equal denominators the numerators are added directly; with different denominators it cross-multiplies and normalizes.

---

#### RationalPolynomial::sub

**Signature**: `pub fn sub(&self, other: &Self) -> Self`

**Description**: Subtraction of rational polynomials. Equivalent to `self.add(&other.neg())`.

---

#### RationalPolynomial::mul

**Signature**: `pub fn mul(&self, other: &Self) -> Self`

**Description**: Multiplication of rational polynomials. Cross-multiplies and normalizes.

---

#### RationalPolynomial::div

**Signature**: `pub fn div(&self, other: &Self) -> Option<Self>`

**Description**: Division of rational polynomials. $\frac{a/b}{c/d} = \frac{ad}{bc}$.

**Returns**: `Some(quotient)` or `None` (when the numerator of the divisor is zero).

---

## Helper Functions

The following functions live in the `ocas_poly::sparse` module:

### monomial_divides

**Signature**: `pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool`

**Description**: Checks whether the monomial $b$ divides $a$ (i.e., $a$ is a multiple of $b$). Returns `true` iff $a_i \geq b_i$ for all $i$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `&[usize]` | Exponent vector of the dividend (the multiple) |
| `b` | `&[usize]` | Exponent vector of the divisor |

**Returns**: `true` if $b$ divides $a$.

---

### monomial_lcm

**Signature**: `pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]>`

**Description**: Computes the least common multiple of two monomials: the component-wise maximum.

$$\text{lcm}(x^a y^b, x^c y^d) = x^{\max(a,c)} y^{\max(b,d)}$$

---

### monomial_are_coprime

**Signature**: `pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool`

**Description**: Checks whether two monomials are coprime (no variable appears in both).

---

## Type Aliases

```rust
/// Square-free factorization result: list of (factor, multiplicity) pairs.
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

/// Complete factorization result: list of (irreducible factor, multiplicity) pairs.
pub type Factors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;
```
