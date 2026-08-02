# Gröbner Bases

oCAS computes Gröbner bases of multivariate polynomial ideals over any
field. Three algorithms are provided, plus order-conversion utilities.
This chapter compares them and explains when to use each.

---

## Scope

| Algorithm | Entry point | Best for |
|---|---|---|
| Buchberger | `GroebnerBasis::buchberger` | Small ideals, teaching |
| **F4** | `f4::f4` | Production use — default |
| **F5** | `f5::f5` | Signature-based; regular sequences |

Order conversion:

| Tool | Entry point | Best for |
|---|---|---|
| Re-run F4 | `GroebnerBasis::reorder` | General ideals |
| **FGLM** | `fglm::fglm` | Zero-dimensional ideals (much faster) |

---

## Buchberger vs F4

Buchberger's algorithm processes S-polynomials one at a time. F4 (Faugère
1999) batches many S-polynomial reductions into a single sparse-matrix
row echelon computation, which is dramatically faster for medium and
large ideals because the linear algebra dominates and can be optimized
(cache-friendly sparse rows, the ℤ_p native `i64` fast path).

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::{SparseMultivariatePolynomial, f4};

let d = RationalDomain;
// cyclic-3 system
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

F4 uses Gebauer–Moeller critical-pair filtering (first and second
criteria plus redundant-pair cleanup) and a per-basis-element
simplification cache, so the matrices it builds are close to minimal.

---

## Monomial Orders and `reorder`

`Lex`, `Grlex`, and `Grevlex` orders are supported. `Lex` bases are what
elimination theory needs, but they are usually the most expensive to
compute. The standard strategy is:

1. compute a `Grevlex` basis (fastest),
2. convert it to `Lex`.

For a general ideal, `reorder` re-interprets the basis under the new
order and re-runs F4:

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

## FGLM: Fast Conversion for Zero-Dimensional Ideals

A zero-dimensional ideal (finitely many common roots) has a finite
*staircase* — the monomials not divisible by any leading monomial. The
FGLM algorithm (Faugère–Gianni–Lazard–Mora 1993) walks monomials of the
target order, computes their normal forms against the existing basis,
and detects linear dependencies. Each dependency yields one polynomial
of the new basis. The cost is `O(n·D³)` field operations where `D` is
the staircase dimension, independent of the F4 cost that produced the
original basis.

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

`fglm` returns `None` when the ideal is positive-dimensional (infinite
staircase). Use `reorder` in that case.

---

## F5: Signature-Based Gröbner Bases

`f5::f5` implements Faugère's signature criterion (2002): S-pairs whose
signature is already present are skipped, which provably avoids all
reductions to zero for regular sequences. Since version 0.19.0, F5 is a
production-grade implementation with performance competitive with F4 on
large ideals.

The `hilbert` module computes the Hilbert numerator of a monomial ideal
by inclusion–exclusion, giving the regularity of the staircase — a sound
degree bound that F4 can use as an early-termination hint.

---

## Benchmarks

Criterion timings (cyclic systems over ℚ and ℤ₁₃, this machine):

| System | Buchberger | F4 | Speedup |
|---|---|---|---|
| cyclic-3 ℚ | 308 µs | 147 µs | 2.1× |
| cyclic-4 ℚ | 3.99 ms | 2.13 ms | 1.9× |
| cyclic-3 ℤ₁₃ | 582 µs | 276 µs | 2.1× |
| cyclic-4 ℤ₁₃ | 6.19 ms | 2.80 ms | 2.2× |

The ℤ_p native `i64` fast path (lazy modular arithmetic in the row
echelon step) is what keeps the finite-field timings close to the
rational ones despite the smaller coefficients.

---

## Monomial Orders

Version 0.19.1 extended the monomial order system with runtime-configurable
orderings beyond the built-in `Lex`, `Grlex`, and `Grevlex`:

| Order | Entry point | Use case |
|---|---|---|
| `WeightOrder` | `WeightOrder::from_slice(&[2, 1])` | Elimination-style orderings via variable weights |
| `BlockOrder` | `BlockOrder::new(boundaries, orders)` | Block elimination: partition variables into groups with independent orderings |

`WeightOrder` compares monomials by the weighted sum $\sum w_i e_i$ in
descending order. `BlockOrder` partitions variables into contiguous blocks,
each compared under its own sub-ordering (Lex or Grevlex).

```rust
use ocas_poly::sparse::WeightOrder;

let order = WeightOrder::from_slice(&[2, 1]);
let p = SparseMultivariatePolynomial::new_with_order(d, 2, order);
```

### MatrixOrder and Elimination

Version 0.23.0 added `MatrixOrder`, a general matrix-based monomial ordering
that supports elimination via weight matrices:

| Order | Entry point | Use case |
|---|---|---|
| `MatrixOrder` | `MatrixOrder::new(matrix)` | General weight-matrix ordering |
| Elimination | `MatrixOrder::elimination_order(elim_vars, n_vars)` | Eliminate first `elim_vars` variables |

```rust
use ocas_poly::sparse::{MatrixOrder, MonomialOrder};

// Elimination order: eliminate x_0, then compare remaining by Grevlex.
let ord = MatrixOrder::elimination_order(1, 3);
```

---

## Ideal Operations

Version 0.23.0 introduced a complete ideal arithmetic library in `ocas_poly::ideal`.
All operations work over `Lex` ordering for consistency with elimination.

| Operation | Entry point | Description |
|---|---|---|
| Membership | `ideal_contains(gens, f, algo)` | Test if $f \in I$ |
| Sum | `ideal_sum(I, J)` | $I + J$ |
| Product | `ideal_product(I, J)$ | $I \cdot J$ |
| Quotient | `ideal_quotient(I, J)` | $I : J$ (Rabinowitsch trick) |
| Saturation | `ideal_saturate(I, J)` | $I : J^\infty$ |
| Intersection | `ideal_intersection(I, J)` | $I \cap J$ (auxiliary variable) |
| Elimination | `eliminate(gens, elim_vars, algo)` | Eliminate variables via Lex GB |

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::sparse::Lex;
use ocas_poly::ideal::{ideal_contains, ideal_saturate};
use ocas_poly::{Algorithm, SparseMultivariatePolynomial};

let d = RationalDomain;
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
    (vec![0, 2], Rational::new(1, 1)),
    (vec![0, 0], Rational::new(-1, 1)),
]); // x² + y² - 1
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 0], Rational::new(1, 1)),
    (vec![0, 1], Rational::new(-1, 1)),
]); // x - y

// Test membership: x² + y² - 1 ∈ ⟨x - y, x² + y² - 1⟩
assert!(ideal_contains(&[f2.clone()], &f1, Algorithm::Auto));
```

---

## Zero-Dimensional Solving

`solve_polynomial_system` classifies systems and finds real solutions:

| Solution type | Description |
|---|---|
| `ZeroDimensional` | Finite number of real solutions via Sturm root isolation |
| `PositiveDimensional` | Infinite solution set; returns Lex GB |
| `Empty` | No solutions (ideal is $\langle 1 \rangle$) |

```rust
use ocas_poly::ideal::solve_polynomial_system;

// Circle ∩ line: x² + y² = 1, x = y
let sol = solve_polynomial_system(&[f1, f2], Algorithm::Auto);
// Returns two real solutions at (±1/√2, ±1/√2)
```

Use `is_zero_dimensional(&gb)` to check dimensionality without solving.

---

## Primary Decomposition and Radical

For zero-dimensional ideals, primary decomposition factors univariate polynomials
in the Lex GB and separates components via saturation:

```rust
use ocas_poly::ideal::{primary_decomposition, ideal_radical};

// (x², xy) = (x) ∩ (x², y)
let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![2, 0], Rational::new(1, 1)),
]);
let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    (vec![1, 1], Rational::new(1, 1)),
]);
let decomp = primary_decomposition(&[f1, f2]);
assert_eq!(decomp.len(), 2);
```

`ideal_radical` computes $\sqrt{I}$:
- **Zero-dimensional**: squarefree decomposition of univariate polynomials
- **Positive-dimensional**: Jacobian saturation ($\sqrt{I} = I : h^\infty$)

```rust
let rad = ideal_radical(&[f1, f2]);
// √(x², xy) = (x)
```

---

## Hilbert Series

The `hilbert` module provides complete Hilbert series computation from Gröbner bases:

| Method | Description |
|---|---|
| `hilbert_function(d)` | $\dim_k (R/I)_d$ at degree $d$ |
| `dimension()` | Krull dimension of $R/I$ |
| `degree()` | Degree of the projective variety |
| `hilbert_polynomial()` | Full polynomial coefficients (Lagrange interpolation) |

```rust
use ocas_poly::groebner::hilbert::hilbert_series;

let hs = hilbert_series(&gb);
println!("H(5) = {}", hs.hilbert_function(5));
println!("dim = {}", hs.dimension());
println!("polynomial = {:?}", hs.hilbert_polynomial());
```
