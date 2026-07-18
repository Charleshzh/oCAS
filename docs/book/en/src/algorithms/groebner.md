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
| F5 (experimental) | `f5::f5` | Research, signature pruning |

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

## F5 and Hilbert Bounds (Experimental)

`f5::f5` implements Faugère's signature criterion (2002): S-pairs whose
signature is already present are skipped, which provably avoids all
reductions to zero for regular sequences. The implementation is kept for
research; on the test-suite ideals the signature rule rarely fires, so
F4 remains the recommended default.

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
