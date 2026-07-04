# Rational Functions, Resultants & Partial Fractions

oCAS provides a complete stack for rational function arithmetic, resultant
computation, and partial fraction decomposition over any Euclidean domain.
These capabilities were added in version 0.12.0, closing the gap with
Symbolica's `rational_polynomial.rs`, `resultant.rs`, and
`partial_fraction.rs`.

---

## Rational Polynomials

A `RationalPolynomial<D, O>` represents an element of the fraction field of a
polynomial ring — that is, a ratio $\frac{p}{q}$ where $p$ and $q$ are
multivariate polynomials over a domain `D`.

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::{RationalPolynomial, SparseMultivariatePolynomial, Lex};

// Create polynomials: x + 1 and x - 1
let x_plus_1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 1,
    vec![(vec![0], Integer::from(1)), (vec![1], Integer::from(1))],
);
let x_minus_1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 1,
    vec![(vec![0], Integer::from(-1)), (vec![1], Integer::from(1))],
);

// (x+1) / (x-1)
let rat = RationalPolynomial::from_num_den(x_plus_1, x_minus_1);
```

### Canonicalization

When constructed via `from_num_den`, the fraction is automatically reduced to
canonical form:

1. The GCD of numerator and denominator is divided out.
2. The leading coefficient of the denominator is normalized (positive for
   ordered domains, 1 for finite fields).

### Arithmetic

All standard operations are supported:

| Operation | Method | Strategy |
|---|---|---|
| Addition | `a.add(&b)` | Cross-multiply, then canonicalize |
| Subtraction | `a.sub(&b)` | Via negation + addition |
| Multiplication | `a.mul(&b)` | Cross-cancel GCDs, then multiply |
| Division | `a.div(&b)` | Via inverse + multiplication |
| Negation | `a.neg()` | Negate numerator |
| Inverse | `a.inv()` | Swap numerator/denominator |
| Power | `a.pow(n)` | Repeated squaring |

---

## Resultants

The resultant of two polynomials $a$ and $b$ is a scalar that is zero if and
only if $a$ and $b$ share a common root (or equivalently, a non-trivial GCD).

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// Res(x^2 + 1, (x+1)^2) = 4
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(0), Integer::from(1),
]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
assert_eq!(a.resultant(&b), Integer::from(4));
```

### Algorithm: Brown's PRS

oCAS uses **Brown's Polynomial Remainder Sequence** algorithm, which avoids
constructing the full Sylvester matrix. The algorithm tracks leading
coefficients and degree differences to compute the resultant from the PRS
using the fundamental theorem formula.

For two polynomials of degree 15, the resultant completes in under 20 ms.

---

## Partial Fractions

Given a proper fraction $\frac{p(x)}{q(x)}$ (where $\deg(p) < \deg(q)$),
the partial fraction decomposition expresses it as a sum of simpler fractions.

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(1, 1), Rational::new(0, 1), Rational::new(-1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
// poly_part is None (proper fraction)
// terms contains the decomposed fractions
```

### Algorithm

The decomposition proceeds in steps:

1. **Polynomial division**: If $\deg(p) \geq \deg(q)$, extract the polynomial
   part via `div_rem`.
2. **Square-free factorization**: Decompose the denominator into square-free
   factors $f_i^{e_i}$.
3. **Diophantine CRT**: For multiple factors, solve the polynomial Chinese
   Remainder Theorem to split the fraction.
4. **p-adic expansion**: For repeated factors ($e_i > 1$), expand the
   numerator with respect to each factor to get individual terms.

The `together` function reverses the decomposition, combining terms back into
a single rational function.

---

## Helper Methods

Several new methods were added to `DenseUnivariatePolynomial` to support
these algorithms:

| Method | Description |
|---|---|
| `extended_gcd_poly(&self, other)` | Extended GCD: returns `(g, s, t)` with $s \cdot a + t \cdot b = g$ |
| `diophantine(polys, b)` | Polynomial CRT solver |
| `p_adic_expansion(&self, p)` | Repeated division expansion |
| `pow(n)` | Polynomial power by repeated squaring |
| `lcoeff()` | Leading coefficient (convenience) |
| `constant()` | Constant term (convenience) |
| `mul_coeff(c)` / `div_coeff(c)` | Scalar multiply/divide all coefficients |
| `content()` | GCD of all coefficients (now public) |

For `SparseMultivariatePolynomial`:

| Method | Description |
|---|---|
| `div_exact(&self, divisor)` | Exact polynomial division (no remainder) |
| `degree_in(var_index)` | Degree in a specific variable |

---

## Rational Reconstruction

The `rational_reconstruction` function recovers a rational number $\frac{n}{d}$
from its modular image $a \equiv \frac{n}{d} \pmod{m}$:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::rational_reconstruction::rational_reconstruction;

// 3/7 mod 101: 7^{-1} mod 101 = 29, so a = 3 * 29 mod 101 = 87
let m = Integer::from(101);
let a = Integer::from(87);
let result = rational_reconstruction(&a, &m);
assert!(result.is_some());
```

This uses the extended Euclidean algorithm (continued fraction approach) and
is useful for modular algorithms that need to lift results back to $\mathbb{Q}$.

---

## Karatsuba Multiplication

Dense polynomial multiplication now uses **Karatsuba's algorithm** for
polynomials with 32 or more coefficients. This reduces the multiplication
complexity from $O(n^2)$ to $O(n^{1.585})$.

The threshold (32) was chosen empirically — below this size, the overhead of
Karatsuba's extra additions and subtractions outweighs the savings from fewer
multiplications.

For two degree-500 polynomials over $\mathbb{Z}$, Karatsuba provides a
significant speedup over the previous schoolbook implementation.

---

## SymPy Migration

| SymPy | oCAS |
|---|---|
| `sp.apart(expr, x)` | `apart(&num, &den)` |
| `sp.together(expr)` | `together(poly_part, &terms)` |
| `sp.resultant(a, b, x)` | `a.resultant(&b)` |
| `sp.Rational(n, d)` | `Rational::new(n, d)` |

---

## Limitations

- Partial fraction decomposition currently works over any `EuclideanDomain`
  using square-free factorization. Full factorization into irreducibles
  requires `IntegerDomain` or `FiniteField`.
- Multivariate partial fractions (`apart_multivariate`) are deferred to 0.13+
  as they require Gröbner F4.
- FFT/NTT multiplication is deferred to a future version.
