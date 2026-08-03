# Gröbner Bases and Ideals

This chapter documents the complete Rust API for Gröbner basis computation and ideal operations in oCAS. All algorithms are accessed through the unified entry point [`groebner_basis`]; ideal operations live in the `ocas_poly::ideal` module.

## Table of contents

- [`Algorithm`](#algorithm) — algorithm selection enum
- [`GroebnerBasis`](#groebnerbasis) — Gröbner basis struct
- [`groebner_basis`](#groebner_basis) — unified entry point
- [`buchberger`](#buchberger) — Buchberger algorithm
- [`f4::f4`](#f4f4) — F4 matrix algorithm
- [`f5::f5`](#f5f5) — F5 signature algorithm
- [`fglm`](#fglm) — FGLM order-change algorithm
- [`HilbertSeries`](#hilbertseries) — Hilbert series
- [`hilbert_series`](#hilbert_series) — computing the Hilbert series
- [`eliminate`](#eliminate) — elimination
- [Ideal operations](#ideal-operations) — `ideal_contains`, `ideal_sum`, `ideal_product`, etc.
- [`PrimaryComponent`](#primarycomponent) — primary component
- [`PolynomialSystemSolution`](#polynomialsystemsolution) — polynomial system solution result

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

**Description**: Algorithm selector for Gröbner basis computation; passed to the unified entry point [`groebner_basis`].

**Variants**:

| Variant | Description |
|---|---|
| `Auto` (default) | Automatically selects the algorithm based on the size and structure of the ideal. Currently routes to F4; the F5 switching threshold will be tuned against the cyclic-n benchmarks in the future. |
| `F4` | Forces the F4 matrix algorithm (Faugère 1999). Batches S-polynomial reduction into sparse matrix row operations, significantly faster than Buchberger for larger ideals. |
| `F5` | Forces the F5 signature algorithm (Faugère 2002). Rejects zero reductions before matrix construction via the syzygy criterion, giving order-of-magnitude speedups on hard ideals (e.g., cyclic-n). |
| `Buchberger` | Forces the classic Buchberger S-polynomial iteration algorithm; suitable for small ideals. |

**Example**:

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

**See also**: [`groebner_basis`](#groebner_basis)

---

## `GroebnerBasis`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerBasis<D: Domain, O: MonomialOrder> {
    pub basis: Vec<SparseMultivariatePolynomial<D, O>>,
}
```

**Description**: A Gröbner basis of a polynomial ideal. The `basis` field stores the list of polynomials of the basis. Constructed via `buchberger`, `f4::f4`, `f5::f5`, or `groebner_basis`; all entry points return a **reduced Gröbner basis** (minimized + auto-reduced).

**Type parameters**:

| Parameter | Constraint | Description |
|---|---|---|
| `D` | `Domain` | Coefficient domain (e.g., `RationalDomain`, `FiniteField`) |
| `O` | `MonomialOrder` | Monomial order (e.g., `Lex`, `Grevlex`) |

**Fields**:

| Field | Type | Description |
|---|---|---|
| `basis` | `Vec<SparseMultivariatePolynomial<D, O>>` | List of polynomials of the basis |

### `GroebnerBasis::buchberger`

```rust
pub fn buchberger(ideal: &[SparseMultivariatePolynomial<D, O>]) -> Self
```

**Description**: Computes a Gröbner basis from the generators using the Buchberger algorithm. Internally filters out zero polynomials, applies Buchberger's first criterion (skips S-polynomials whose leading terms are coprime), with a maximum of 10,000 iterations.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | List of generators of the ideal |

**Returns**: `GroebnerBasis<D, O>` — the raw basis, not minimized/auto-reduced.

**Note**: Requires the coefficient domain to support exact division (i.e., to be a field). Panics if division fails. The convenience function [`buchberger`](#buchberger) (a free function) additionally calls `minimize().auto_reduce()`.

**Example**:

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

**See also**: [`buchberger`](#buchberger) (convenience free function), [`f4::f4`](#f4f4), [`f5::f5`](#f5f5)

### `minimize`

```rust
pub fn minimize(mut self) -> Self
```

**Description**: Minimizes the basis — removes polynomials whose leading monomial is divisible by the leading monomial of another element.

**Returns**: `Self` — the minimized basis (consumes `self`).

**See also**: [`auto_reduce`](#auto_reduce)

### `auto_reduce`

```rust
pub fn auto_reduce(mut self) -> Self
```

**Description**: Auto-reduces the basis — reduces each element against the others and makes it monic. Processes elements in ascending order of leading monomial, guaranteeing the standard reduced Gröbner basis property: no monomial of any basis element is divisible by the leading monomial of another basis element.

**Returns**: `Self` — the reduced basis (consumes `self`).

**See also**: [`minimize`](#minimize)

### `is_groebner_basis`

```rust
pub fn is_groebner_basis(&self) -> bool
```

**Description**: Verifies that this basis is indeed a Gröbner basis — checks that all S-polynomials reduce to zero.

**Returns**: `bool` — `true` if it is a valid Gröbner basis.

**Example**:

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

**See also**: [`groebner_basis`](#groebner_basis)

### `reorder`

```rust
pub fn reorder<O2: MonomialOrder>(&self) -> GroebnerBasis<D, O2>
where
    D: 'static,
```

**Description**: Changes the monomial order of the Gröbner basis. Reinterprets the polynomials under the target order and re-runs F4. This is the simple order-change path; for zero-dimensional ideals, use [`fglm`](#fglm) for the faster $O(n \cdot D^3)$ conversion.

**Type parameters**:

| Parameter | Description |
|---|---|
| `O2` | Target monomial order |

**Returns**: `GroebnerBasis<D, O2>` — a Gröbner basis under the target order.

**Example**:

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

**See also**: [`fglm`](#fglm) (faster order change for zero-dimensional ideals)

---

## `groebner_basis`

```rust
pub fn groebner_basis<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    algo: Algorithm,
) -> GroebnerBasis<D, O>
```

**Description**: Unified entry point for Gröbner basis computation. Selects the algorithm according to the `algo` parameter (`Auto`/`F4`/`F5`/`Buchberger`). Zero polynomials are filtered out internally by the backends.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | List of generators of the ideal |
| `algo` | `Algorithm` | Algorithm selection |

**Returns**: `GroebnerBasis<D, O>` — a reduced Gröbner basis.

**Example**:

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

**See also**: [`Algorithm`](#algorithm), [`GroebnerBasis`](#groebnerbasis)

---

## `buchberger`

```rust
pub fn buchberger<D: Domain, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**Description**: Convenience function — computes a Gröbner basis and then minimizes and auto-reduces it (equivalent to `GroebnerBasis::buchberger(ideal).minimize().auto_reduce()`).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | List of generators of the ideal |

**Returns**: `GroebnerBasis<D, O>` — a reduced Gröbner basis.

**See also**: [`GroebnerBasis::buchberger`](#groebnerbasisbuchberger), [`groebner_basis`](#groebner_basis)

---

## `f4::f4`

```rust
pub fn f4<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**Description**: Computes a Gröbner basis using the F4 matrix algorithm (Faugère 1999). Replaces sequential S-polynomial reduction with batched sparse matrix row operations (Gaussian elimination). For `FiniteField` coefficient domains, the native ℤ_p fast path (`FpPoly` i64 modular arithmetic) is used automatically, avoiding BigInt overhead.

**Internal optimizations**:

- **Gebauer–Moeller pair selection**: chain criterion + update criterion with strict divisibility guards
- **`SimpCache`**: caches the results of multiplying polynomials by monomials, avoiding recomputation
- **`DivisorIndex`**: support-bitmask-based O(1)-ish reducer lookup (exact bitmasks for variable indices ≤ 63; degrades to a correctness-conservative filter above 63)
- **Two-pointer sparse row subtraction** `sub_scaled_fp`: O(nnz) complexity
- **Native ℤ_p fast path**: `FpPoly` uses i64 modular arithmetic, converting to/from BigInt only when reading inputs and writing outputs

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | List of generators of the ideal |

**Returns**: `GroebnerBasis<D, O>` — a reduced Gröbner basis.

**Example**:

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

**See also**: [`Algorithm::F4`](#algorithm), [`f5::f5`](#f5f5)

---

## `f5::f5`

```rust
pub fn f5<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O>
```

**Description**: Computes a Gröbner basis using the F5 signature algorithm (Faugère 2002). Attaches a signature (`module_pos`, `monomial`) to each polynomial and rejects zero reductions before matrix construction via the syzygy criterion, achieving order-of-magnitude speedups on hard ideals (e.g., cyclic-n).

**Algorithm highlights**:

- **Signatures**: `(module_pos, monomial)` records the history of a polynomial — `module_pos` is the index of the input generator and `monomial` is the monomial multiple applied. Signatures are compared under the **pot** (position-over-term) order: module position first (smaller preferred), then by the monomial order `O`.
- **Syzygy tracking**: when a matrix row reduces to zero, its signature is a syzygy. Rows whose future signatures are monomial multiples of known syzygies are skipped immediately (the F5 syzygy criterion).
- **Incremental processing**: generators are processed one at a time; each new generator triggers a round of matrix construction and reduction.

The current implementation provides a general-domain F5 core. For `FiniteField` coefficient domains, the native ℤ_p fast path is used automatically.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, O>]` | List of generators of the ideal |

**Returns**: `GroebnerBasis<D, O>` — a reduced Gröbner basis (identical to the F4 output).

**Example**:

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

**See also**: [`Algorithm::F5`](#algorithm), [`f4::f4`](#f4f4)

---

## `fglm`

```rust
pub fn fglm<D: Domain, O2: MonomialOrder>(
    gb: &GroebnerBasis<D, impl MonomialOrder>,
) -> Option<GroebnerBasis<D, O2>>
```

**Description**: The FGLM order-change algorithm (Faugère–Gianni–Lazard–Mora 1993). Converts a Gröbner basis of a **zero-dimensional** ideal from one monomial order to another, with complexity $O(n \cdot D^3)$ ($D$ is the vector space dimension of $R/I$), much faster than re-running F4 for large zero-dimensional ideals.

**Algorithm steps**:

1. Compute the staircase — the set of monomials not divisible by any leading monomial
2. Traverse the staircase monomials in the target order
3. For each monomial, compute its normal form under the current basis
4. If the normal form is linearly dependent on the vectors seen so far, construct a new basis element from the coefficient relation
5. Otherwise, add it to the set of seen vectors

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `gb` | `&GroebnerBasis<D, impl MonomialOrder>` | The input Gröbner basis (must be reduced) |

**Returns**: `Option<GroebnerBasis<D, O2>>` — `Some(gb)` is the converted basis; `None` means the ideal is not zero-dimensional (the staircase is infinite).

**Example**:

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

**See also**: [`GroebnerBasis::reorder`](#reorder) (simple order change, works for any dimension), [`is_zero_dimensional`](#is_zero_dimensional)

---

## `HilbertSeries`

```rust
#[derive(Debug, Clone)]
pub struct HilbertSeries {
    pub numerator: Vec<i64>,
    pub denominator_power: usize,
}
```

**Description**: The Hilbert series of the quotient ring $R/I$, represented as the rational function $H(t) = N(t) / (1-t)^n$. The numerator $N(t)$ is stored as a coefficient vector (`numerator[i]` is the coefficient of $t^i$), and the denominator is $(1-t)^n$ ($n$ is the number of variables).

**Fields**:

| Field | Type | Description |
|---|---|---|
| `numerator` | `Vec<i64>` | Coefficients of the numerator $N(t)$ (starting from the constant term) |
| `denominator_power` | `usize` | Power of $(1-t)$ (= the number of variables $n$) |

### `hilbert_function`

```rust
pub fn hilbert_function(&self, degree: usize) -> i64
```

**Description**: Computes the Hilbert function value $\dim_k (R/I)_d$ at degree $d$. Uses the formula $H(d) = [t^d] N(t) / (1-t)^n$, where the coefficient of $t^k$ in $(1-t)^{-n}$ is $\binom{n+k-1}{k}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `degree` | `usize` | The degree $d$ |

**Returns**: `i64` — the value of $\dim_k (R/I)_d$.

### `dimension`

```rust
pub fn dimension(&self) -> usize
```

**Description**: Computes the Krull dimension of $R/I$. Obtained by checking the order of vanishing of the numerator at $t=1$ (differentiating repeatedly until the value at $t=1$ is nonzero).

**Returns**: `usize` — the Krull dimension. Zero-dimensional ideals return 0.

### `degree`

```rust
pub fn degree(&self) -> i64
```

**Description**: Computes the degree of the projective variety. For a well-formed Hilbert series, this equals the value of the numerator at $t=1$ after removing the dimension factors.

**Returns**: `i64` — the degree of the projective variety.

### `hilbert_polynomial`

```rust
pub fn hilbert_polynomial(&self) -> Vec<f64>
```

**Description**: Computes the coefficients of the Hilbert polynomial $P(d)$ ($H(d) = P(d)$ for $d \gg 0$) using Lagrange interpolation. Returns the coefficients in ascending power order (`result[i]` is the coefficient of $d^i$); the degree of the polynomial is `self.dimension()`.

**Returns**: `Vec<f64>` — the Hilbert polynomial coefficients (ascending power order).

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::groebner::hilbert::hilbert_series;

let d = RationalDomain;
// Ideal (x² - 1, y² - 1): LM ideal (x², y²), numerator N(t) = 1 - 2t² + t⁴
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
assert_eq!(hs.dimension(), 0); // zero-dimensional: 4 solutions
assert_eq!(hs.degree(), 4); // degree = number of solutions
println!("dimension: {}, degree: {}", hs.dimension(), hs.degree());
for d in 0..10 {
    println!("H({}) = {}", d, hs.hilbert_function(d));
}
let hp = hs.hilbert_polynomial();
println!("Hilbert polynomial coefficients: {:?}", hp);
```

**See also**: [`hilbert_series`](#hilbert_series)

---

## `hilbert_series`

```rust
pub fn hilbert_series(
    gb: &GroebnerBasis<RationalDomain, Lex>,
) -> HilbertSeries
```

**Description**: Computes the Hilbert series of $R/I$ from a Gröbner basis. Uses Macaulay's theorem: the Hilbert series of $R/I$ equals the Hilbert series of $R/\text{LM}(I)$ (the leading-term ideal). The Hilbert numerator of the monomial ideal $\langle m_1, \dots, m_s \rangle$ is computed via inclusion–exclusion: $N(t) = \sum_k (-1)^k \sum_{|S|=k} t^{\deg \text{lcm}(S)}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `gb` | `&GroebnerBasis<RationalDomain, Lex>` | A Lex-order Gröbner basis over the rationals |

**Returns**: `HilbertSeries` — the Hilbert series of the quotient ring.

**Example**:

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

**See also**: [`HilbertSeries`](#hilbertseries), [`is_zero_dimensional`](#is_zero_dimensional)

---

## `eliminate`

```rust
pub fn eliminate<D: Domain + 'static>(
    ideal: &[SparseMultivariatePolynomial<D, Lex>],
    elim_vars: usize,
    algo: Algorithm,
) -> GroebnerBasis<D, Lex>
```

**Description**: Eliminates variables from an ideal. Returns a Gröbner basis of $I \cap k[x_{\text{elim\_vars}}, \dots, x_{n-1}]$, i.e., the polynomials that do not contain the first `elim_vars` variables. Uses the Lex order — under Lex, the reduced Gröbner basis of an ideal automatically contains generators of the elimination ideal.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ideal` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of the ideal (must be in Lex order) |
| `elim_vars` | `usize` | The number of variables to eliminate (eliminates $x_0, \dots, x_{\text{elim\_vars}-1}$) |
| `algo` | `Algorithm` | Gröbner basis algorithm selection |

**Returns**: `GroebnerBasis<D, Lex>` — a Gröbner basis of the elimination ideal.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{SparseMultivariatePolynomial, eliminate, Algorithm};

let d = RationalDomain;
// Ideal: x + y + z, x*y + x*z in k[x,y,z]; eliminate x
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
// The result should lie in k[y,z]
for p in &elim.basis {
    assert!(p.degree_in(0) == 0, "the eliminated variable x should not appear");
}
```

**See also**: [`groebner_basis`](#groebner_basis), [`ideal_quotient`](#ideal_quotient)

---

## Ideal operations

All ideal operations live in the `ocas_poly::ideal` module and use the `Lex` order for consistency with elimination computations.

### `ideal_contains`

```rust
pub fn ideal_contains<D: Domain + 'static>(
    generators: &[SparseMultivariatePolynomial<D, Lex>],
    f: &SparseMultivariatePolynomial<D, Lex>,
    algo: Algorithm,
) -> bool
```

**Description**: Tests whether $f$ belongs to the ideal generated by the generators. Computes a Gröbner basis of the ideal and reduces $f$; $f \in I$ if and only if the remainder is zero.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of the ideal |
| `f` | `&SparseMultivariatePolynomial<D, Lex>` | The polynomial to test |
| `algo` | `Algorithm` | Gröbner basis algorithm selection |

**Returns**: `bool` — `true` if $f \in I$.

**Example**:

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

**See also**: [`groebner_basis`](#groebner_basis)

---

### `ideal_sum`

```rust
pub fn ideal_sum<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**Description**: The sum of two ideals $I + J = \langle f_1, \dots, f_m, g_1, \dots, g_n \rangle$. Merges the generators and computes a Gröbner basis.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $I$ |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $J$ |

**Returns**: `GroebnerBasis<D, Lex>` — a reduced Gröbner basis of $I + J$.

**Example**:

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

**See also**: [`ideal_intersection`](#ideal_intersection)

---

### `ideal_product`

```rust
pub fn ideal_product<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**Description**: The product of two ideals $I \cdot J = \langle f_i \cdot g_j \rangle$. Computes all products $f_i g_j$ and then a Gröbner basis.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $I$ |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $J$ |

**Returns**: `GroebnerBasis<D, Lex>` — a reduced Gröbner basis of $I \cdot J$.

**Example**:

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

**See also**: [`ideal_quotient`](#ideal_quotient)

---

### `ideal_quotient`

```rust
pub fn ideal_quotient<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**Description**: The quotient of ideals $I : J = \{f : f \cdot g \in I, \forall g \in J\}$. Uses the Rabinowitsch trick for each generator $g$ of $J$: compute $\text{GB}(I \cup \{1 - w \cdot g\})$ in the extended ring $k[x_1, \dots, x_n, w]$, eliminate $w$, and intersect the results.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators_i` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $I$ |
| `generators_j` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $J$ |

**Returns**: `GroebnerBasis<D, Lex>` — a reduced Gröbner basis of $I : J$.

**Example**:

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

**See also**: [`ideal_saturate`](#ideal_saturate), [`ideal_contains`](#ideal_contains)

---

### `ideal_saturate`

```rust
pub fn ideal_saturate<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**Description**: The saturation of an ideal $I : J^\infty = \bigcup_k (I : J^k)$. Iteratively computes $I : J$, $(I : J) : J$, etc., until stabilization.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators_i` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $I$ |
| `generators_j` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $J$ |

**Returns**: `GroebnerBasis<D, Lex>` — a reduced Gröbner basis of $I : J^\infty$.

**Example**:

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

**See also**: [`ideal_quotient`](#ideal_quotient)

---

### `ideal_intersection`

```rust
pub fn ideal_intersection<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

**Description**: The intersection of two ideals $I \cap J$. Uses an auxiliary variable $t$: $I \cap J = \langle t \cdot f_i, (1-t) \cdot g_j \rangle \cap k[x_1, \dots, x_n]$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators_a` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $I$ |
| `generators_b` | `&[SparseMultivariatePolynomial<D, Lex>]` | The generators of ideal $J$ |

**Returns**: `GroebnerBasis<D, Lex>` — a reduced Gröbner basis of $I \cap J$.

**Example**:

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

**See also**: [`ideal_sum`](#ideal_sum), [`eliminate`](#eliminate)

---

### `ideal_radical`

```rust
pub fn ideal_radical(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> GroebnerBasis<RationalDomain, Lex>
```

**Description**: Computes the radical $\sqrt{I}$ of the ideal.

- **Zero-dimensional ideals**: computed via the square-free factorization of the univariate polynomials in the Lex GB (the generators of $\sqrt{I}$ consist of the square-free univariate polynomials for each variable together with the non-univariate basis elements).
- **Positive-dimensional ideals**: uses the Jacobian saturation method (a simplified Kemper algorithm): $\sqrt{I} = I : h^\infty$. The current implementation picks $h$ heuristically — the non-trivial partial derivative of smallest total degree (not the exact GCD). If all partial derivatives are constants or zero (trivial Jacobian), the original GB is returned conservatively.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The generators of the ideal |

**Returns**: `GroebnerBasis<RationalDomain, Lex>` — a reduced Gröbner basis of $\sqrt{I}$.

**Example**:

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

**See also**: [`primary_decomposition`](#primary_decomposition)

---

### `is_zero_dimensional`

```rust
pub fn is_zero_dimensional(gb: &GroebnerBasis<RationalDomain, Lex>) -> bool
```

**Description**: Checks whether the ideal is zero-dimensional. The ideal is zero-dimensional if and only if for every variable $x_i$ some leading monomial in the GB is a pure power $x_i^N$ (equivalently, the staircase/set of standard monomials is finite).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `gb` | `&GroebnerBasis<RationalDomain, Lex>` | A Lex-order Gröbner basis over the rationals |

**Returns**: `bool` — `true` if zero-dimensional.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, GroebnerBasis, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::ideal::is_zero_dimensional;

let d = RationalDomain;
// x² - 1, y - x → zero-dimensional
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

**See also**: [`solve_polynomial_system`](#solve_polynomial_system), [`fglm`](#fglm)

---

### `solve_polynomial_system`

```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

**Description**: Solves a zero-dimensional polynomial system. Converts the GB to Lex order, extracts the univariate polynomial for each variable, and solves by back substitution. Returns the `PolynomialSystemSolution` enum, distinguishing zero-dimensional (finitely many solutions), positive-dimensional (infinite solution set), and empty ($\langle 1 \rangle$).

**Solution procedure**:

1. Compute a Gröbner basis
2. Check whether it is $\langle 1 \rangle$ (empty set)
3. Check whether it is zero-dimensional
4. Convert to a triangular decomposition in Lex order
5. Back-substitute starting from the last variable: for each variable, find the real roots of the univariate polynomial (Sturm's theorem isolation + refinement to $10^{-14}$), recursively substituting the known values

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `equations` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The system of equations (polynomials = 0) |
| `algo` | `Algorithm` | Gröbner basis algorithm selection |

**Returns**: `PolynomialSystemSolution` — the solution result.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, SparseMultivariatePolynomial, groebner_basis};
use ocas_poly::ideal::{solve_polynomial_system, PolynomialSystemSolution};

let d = RationalDomain;
// x² + y² - 1, x - y → solutions at (±1/√2, ±1/√2)
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

**See also**: [`PolynomialSystemSolution`](#polynomialsystemsolution), [`is_zero_dimensional`](#is_zero_dimensional)

---

### `primary_decomposition`

```rust
pub fn primary_decomposition(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> Vec<PrimaryComponent>
```

**Description**: Computes the primary decomposition of the ideal.

- **Zero-dimensional ideals**: factor the univariate polynomial in the first variable of the Lex GB and separate the primary components by saturation (for each factor $f_i$, compute $I : (\prod_{j \neq i} f_j)^\infty$).
- **Positive-dimensional ideals**: conservatively return a single primary component (the original GB itself as both primary and prime).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The generators of the ideal |

**Returns**: `Vec<PrimaryComponent>` — the list of primary components.

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::primary_decomposition;
use ocas_poly::SparseMultivariatePolynomial;

let d = RationalDomain;
// (x², xy) is positive-dimensional → one conservative component
// (ideally (x², xy) = (x) ∩ (x², y))
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let decomp = primary_decomposition(&[f1, f2]);
assert!(decomp.len() >= 1);
```

**See also**: [`PrimaryComponent`](#primarycomponent), [`ideal_radical`](#ideal_radical)

---

### `is_prime_ideal`

```rust
pub fn is_prime_ideal(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> bool
```

**Description**: Tests whether the ideal is prime.

- **Zero-dimensional ideals**: checks whether the univariate polynomials in the Lex GB are irreducible (for polynomials of degree ≤ 3, the rational root theorem is used).
- **Positive-dimensional ideals**: conservatively returns `false` (a complete implementation would require checking the irreducibility of the variety, which is not yet implemented).

**Note**: This is a conservative approximation — it never returns a false positive (a non-prime ideal reported as prime); it can only return a false negative (a prime ideal reported as non-prime).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The generators of the ideal |

**Returns**: `bool` — `true` for a (confirmed) prime ideal.

**See also**: [`is_primary_ideal`](#is_primary_ideal), [`primary_decomposition`](#primary_decomposition)

---

### `is_primary_ideal`

```rust
pub fn is_primary_ideal(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> bool
```

**Description**: Tests whether the ideal is primary. A primary ideal has exactly one associated prime ideal, i.e., `primary_decomposition` returns at most one component.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `generators` | `&[SparseMultivariatePolynomial<RationalDomain, Lex>]` | The generators of the ideal |

**Returns**: `bool` — `true` if the ideal is primary.

**See also**: [`primary_decomposition`](#primary_decomposition), [`is_prime_ideal`](#is_prime_ideal)

---

## `PrimaryComponent`

```rust
#[derive(Debug, Clone)]
pub struct PrimaryComponent {
    pub primary: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
    pub prime: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
}
```

**Description**: A primary component of an ideal — a primary ideal together with its associated prime ideal (i.e., the radical of the primary ideal).

**Fields**:

| Field | Type | Description |
|---|---|---|
| `primary` | `Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>` | The generators of the primary ideal |
| `prime` | `Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>` | The generators of the associated prime ideal (the radical) |

**See also**: [`primary_decomposition`](#primary_decomposition)

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

**Description**: The result of solving a polynomial system, classified by the type of the solution set.

**Variants**:

| Variant | Description |
|---|---|
| `ZeroDimensional(ZeroDimSolutions)` | Finitely many real solutions (zero-dimensional ideal) |
| `PositiveDimensional(GroebnerBasis<RationalDomain, Lex>)` | Infinite solution set (positive-dimensional ideal); returns a Lex-order Gröbner basis |
| `Empty` | No solutions (the ideal is $\langle 1 \rangle$) |

**See also**: [`solve_polynomial_system`](#solve_polynomial_system)

### `ZeroDimSolutions`

```rust
#[derive(Debug, Clone)]
pub struct ZeroDimSolutions {
    pub solutions: Vec<RealSolution>,
    pub vector_space_dimension: usize,
}
```

**Fields**:

| Field | Type | Description |
|---|---|---|
| `solutions` | `Vec<RealSolution>` | The list of real solutions found |
| `vector_space_dimension` | `usize` | The vector space dimension of the quotient ring $k[x_1,\dots,x_n]/I$ (the total number of solutions over $\mathbb{C}$, counted with multiplicity) |

### `RealSolution`

```rust
#[derive(Debug, Clone)]
pub struct RealSolution {
    pub values: Vec<f64>,
    pub multiplicity: usize,
}
```

**Fields**:

| Field | Type | Description |
|---|---|---|
| `values` | `Vec<f64>` | The value of each variable (one per variable) |
| `multiplicity` | `usize` | The algebraic multiplicity of the solution |

---

## Module path quick reference

| Function/Type | Full path |
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
