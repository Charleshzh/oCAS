# Polynomial Factorization

oCAS implements polynomial factorization for both univariate and bivariate
polynomials over the integers and over prime finite fields. This chapter
describes the algorithms, their scope, and the public APIs used to access them.

---

## Scope

Current factorization support covers:

| Domain | Polynomials | Algorithms |
|---|---|---|
| $\mathbb{Z}[x]$ | Univariate | Square-free factorization, Berlekamp–Zassenhaus, Hensel lifting |
| $\mathbb{F}_p[x]$ | Univariate | Square-free factorization, Berlekamp |
| $\mathbb{Z}[x,y]$ | Bivariate (any LC) | Wang's Hensel lifting; constant-LC fast path, non-constant LC via EEZ |
| $\mathbb{F}_p[x,y]$ | Bivariate (monic in $x$) | Hensel lifting over $\mathbb{F}_p$ |
| $\mathbb{Z}[x_1,\dots,x_n]$ | Multivariate | Wang EEZ Hensel lifting + leading-coefficient preprocessing + p-adic coefficient Hensel lift (non-constant LC) + skeleton-interpolation Diophantine + Zassenhaus recombination |
| $\mathbb{F}_p[x_1,\dots,x_n]$ | Multivariate | EEZ Hensel lifting (with characteristic-$p$ $p$-th power handling; non-constant LC via field Wang preprocessing) |
| $\mathbb{Q}(\alpha)[x]$ | Univariate | Trager's algorithm: norm + factorization over $\mathbb{Q}$ + modular GCD over $\mathrm{GF}(p^d)$ |

Since 0.16.0, multivariate factorization with more than two variables is
supported. Since 0.16.1, non-constant leading coefficients in the main
variable are fully handled via a p-adic coefficient Hensel lift that imposes
Wang's reconstructed leading coefficients at every iteration; the bivariate
constant-LC fast path is retained for efficiency. Since 0.17.0, univariate
factorization over algebraic number fields is supported (Trager).

---

## Univariate Factorization over a Finite Field

For a prime field $\mathbb{F}_p$, oCAS uses **Berlekamp's algorithm**. The
polynomial is first made square-free, then the kernel of the Frobenius matrix
$Q - I$ is computed. Each basis vector of the kernel gives a non-trivial
factorization via gcds with elements of the kernel.

```rust
use ocas_domain::{FiniteField, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let p = FiniteField::new(Integer::from(101));
let mut f = DenseUnivariatePolynomial::<FiniteField>::from_coeffs(p.clone(), vec![
    p.element(1), // constant term
    p.element(0),
    p.element(1), // x^2
]);
let factors = f.factor();
```

The returned factorization is a list of `(factor, multiplicity)` pairs. Over
finite fields, multiplicities are always `1` after the square-free step.

---

## Univariate Factorization over the Integers

For $\mathbb{Z}[x]$, oCAS combines **square-free factorization** with
**Berlekamp–Zassenhaus-style Hensel lifting**. The high-level steps are:

1. Compute the content and reduce to a primitive polynomial.
2. Compute the square-free decomposition.
3. Choose a small prime $p$ such that the reduction stays square-free and has
the same degree as the input.
4. Factor modulo $p$ using Berlekamp.
5. Lift the modular factors to factors over $\mathbb{Z}[x]$ using Hensel
lifting.

```rust
use ocas_domain::IntegerDomain;
use ocas_poly::DenseUnivariatePolynomial;

let f = DenseUnivariatePolynomial::<IntegerDomain>::from_coeffs(
    IntegerDomain,
    vec![1.into(), 0.into(), 1.into()], // x^2 + 1
);
let factors = f.factor();
```

---

## Bivariate Factorization over the Integers

`ocas-poly` provides bivariate factorization over $\mathbb{Z}[x,y]$ using
**Wang's Hensel lifting**, assuming the polynomial is monic in $x$.

The algorithm:

1. Choose an evaluation point $y = \alpha$ so that the univariate image
$f(x, \alpha)$ is square-free and has the fewest irreducible factors.
2. Factor the univariate image over $\mathbb{Z}[x]$.
3. Lift the univariate factors back to bivariate factors by correcting the
Taylor coefficients of $f$ around $y = \alpha$.

The correction step uses rational Bézout coefficients for the univariate
factors, then reconstructs integral corrections with integer division. If the
reconstruction fails (non-integral remainder) or the lifted product does not
match the original polynomial, the implementation tries the next candidate
evaluation point and eventually falls back to returning the polynomial as
irreducible.

```rust
use ocas_domain::IntegerDomain;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::Lex;

type MPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

let domain = IntegerDomain;
let f = MPoly::from_terms(domain, 2, vec![
    (vec![3, 0], 1.into()),  // x^3
    (vec![2, 1], 1.into()),  // x^2*y
    (vec![2, 0], 2.into()),  // 2*x^2
    (vec![1, 1], 1.into()),  // x*y
    (vec![1, 0], 1.into()),  // x
    (vec![0, 2], 1.into()),  // y^2
    (vec![0, 1], 3.into()),  // 3*y
    (vec![0, 0], 2.into()),  // 2
]);

let factors = f.factor();
// factors contains (x^2 + y + 1, 1) and (x + y + 2, 1)
```

---

## Bivariate Factorization over a Finite Field

Over $\mathbb{F}_p[x,y]$, the same Hensel-lifting structure is used, but the
arithmetic is performed directly in the finite field. Bézout coefficients are
computed with finite-field gcds, and all corrections are guaranteed to remain in
the field, so no integral reconstruction is required.

```rust
use ocas_domain::FiniteField;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::Lex;

let domain = FiniteField::new(101.into());
type FpPoly = SparseMultivariatePolynomial<FiniteField, Lex>;

let f = FpPoly::from_terms(domain.clone(), 2, vec![
    (vec![2, 0], 1.into()), // x^2
    (vec![0, 1], 1.into()), // y
    (vec![0, 0], 1.into()), // 1
]);
let factors = f.factor();
```

---

## C/C++ Polynomial API

The C bindings expose opaque handles for bivariate integer and finite-field
polynomials. Polynomials are created from ASCII strings, factored, printed, and
freed through the C API.

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;
    OcasPolyZ *f = ocas_poly_z_create("x^2 + y + 1", &err);
    if (f == NULL) {
        fprintf(stderr, "parse error: %s\n", ocas_error_last_message());
        return 1;
    }

    OcasPolyFactorArray factors = {0};
    int rc = ocas_poly_z_factor(f, &factors, &err);
    if (rc != OCAS_OK) {
        fprintf(stderr, "factor error: %s\n", ocas_error_last_message());
        ocas_poly_z_free(f);
        return 1;
    }

    printf("factors: %zu\n", factors.len);
    for (size_t i = 0; i < factors.len; ++i) {
        OcasPolyZ *factor = (OcasPolyZ *)factors.factors[i].poly;
        char *s = ocas_poly_z_to_string(factor, &err);
        printf("  %s^%zu\n", s, factors.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_z_free(factor);
    }
    ocas_poly_factor_array_free(&factors);
    ocas_poly_z_free(f);
    return 0;
}
```

For finite fields, use `OcasPolyFp`, `ocas_poly_fp_create`, and
`ocas_poly_fp_factor`. The factor array stores `void*` polynomial handles, so
the caller must cast the pointer back to the correct concrete type before
printing or freeing it.

---

## Multivariate Factorization (Wang EEZ)

Since 0.16.0, oCAS factors polynomials in any number of variables over
$\mathbb{Z}$ and $\mathbb{F}_p$ using **Wang's EEZ (Evaluation and
EZ-lifting) algorithm**:

1. **Square-free factorization** (Yun): differentiate in the main variable
$x_1$ and peel repeated factors via $n$-variate GCD (dense recursive
evaluation–interpolation); characteristic-$p$ $p$-th powers are shrunk and
expanded.
2. **Sample-point search**: substitute the secondary variables
$x_2,\dots,x_n$ at sample points so the univariate image keeps its degree
and stays square-free.
3. **Wang leading-coefficient preprocessing**: factor the leading
coefficient $\ell(x_2,\dots)$ and distribute its factors among the
univariate factors using their pairwise-coprime integer images, reconstructing
the true leading coefficients $\ell_i$.
4. **Variable-by-variable Hensel lifting**: lift the univariate factors back
to multivariate ones through the ideals $(x_k - a_k)$, solving multivariate
Diophantine equations at each step.
5. **Zassenhaus recombination**: when the lifted factors are finer than the
true factors, enumerate subsets and trial-divide to obtain the irreducible
factors.

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};

// (x + y + z)(x - y + 2z)
let f = SparseMultivariatePolynomial::<IntegerDomain, Lex>::from_terms(
    IntegerDomain, 3,
    vec![
        (vec![2, 0, 0], Integer::from(1)),
        (vec![1, 0, 1], Integer::from(3)),
        (vec![0, 1, 1], Integer::from(1)),
        (vec![0, 2, 0], Integer::from(-1)),
        (vec![0, 0, 2], Integer::from(2)),
    ],
);
let factors = f.factor();
```

The result is a list of `(factor, multiplicity)` pairs, normalized to be
primitive with a positive leading coefficient.

---

## Factorization over Algebraic Number Fields (Trager)

Since 0.17.0, oCAS factors univariate polynomials over an algebraic number
field $K = \mathbb{Q}(\alpha)$, where $\alpha$ is given by a monic minimal
polynomial over $\mathbb{Q}$. The domain lives in `ocas-domain` as
`AlgebraicExtension<D>` (`AlgebraicNumberField` for $D = \mathbb{Q}$);
the same type with $D = \mathbb{F}_p$ represents $\mathrm{GF}(p^d)$.

The factorizer implements **Trager's algorithm**:

1. **Square-free factorization** (Yun) over $K$, using the modular GCD
   below instead of a pseudo-remainder sequence.
2. **Norm with shift**: for $s = 0, 1, 2, \dots$ compute the norm of
   $g(x) = f(x - s\alpha)$ down to $\mathbb{Q}[x]$ by
   evaluation–interpolation of the scalar resultant
   $\operatorname{Res}_\alpha(m, g)$, until the norm is square-free
   (checked modulo small primes; acceptance is exact).
3. **Rational factorization**: the square-free norm is factored over
   $\mathbb{Z}$ (Hensel path).
4. **Modular GCD over $K$**: each rational factor is mapped into $K[x]$
   and its GCD with $g$ is computed by the modular method — map to
   $\mathrm{GF}(p^d)$ for primes with $m$ irreducible mod $p$, combine
   monic modular GCDs by CRT, rational-reconstruct the coefficients, and
   verify by trial division.
5. **Shift back** with $x \mapsto x + s\alpha$ and normalize monic.

A **rational fast path** applies when all coefficients of $f$ are rational
constants: $f$ is first factored over $\mathbb{Q}$, and only the
$\mathbb{Q}$-irreducible factors go through Trager's norm machinery.

```rust
use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

// ℚ(√2): minimal polynomial α² − 2.
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
// x² − 2 = (x − √2)(x + √2) splits over ℚ(√2).
let f = DenseUnivariatePolynomial::from_coeffs(
    field.clone(),
    vec![field.from_base(Rational::new(-2, 1)), field.zero(), field.one()],
);
let factors = f.factor();
assert_eq!(factors.len(), 2);
```

Benchmarks (criterion group `poly_factor_anf`): degree ≤ 12 over
$\mathbb{Q}(\sqrt2)$, $\mathbb{Q}(\sqrt[3]{2})$, and the cyclotomic field
$\mathbb{Q}(\zeta_5)$ all factor in well under 100 ms.

---

## Limitations and Future Work

- Non-constant leading coefficients in the main variable are fully supported
  since 0.16.1 via a p-adic coefficient Hensel lift with skeleton-interpolation
  Diophantine.
- The bivariate Wang Hensel limitation that the leading coefficient be a
constant is subsumed by the 0.16.0 arbitrary-multivariate (EEZ) path; bivariate
non-constant-LC inputs now dispatch to the EEZ path automatically.
- Since 0.16.2 the $\mathbb{F}_p$ multivariate path also supports non-constant
  leading coefficients via field Wang LC preprocessing
  (`wang_reconstruct_lcoeffs_fp`) and `eez_lift_imposed`.
- Algebraic number field factorization (0.17.0) is currently **univariate
  only**; extension-field factorization of multivariate polynomials and a
  `RootOf`/radical syntax in the parser are future work. The minimal
  polynomial's irreducibility is the caller's responsibility (over a
  reducible modulus the ring has zero divisors and inversion fails).
- The adaptive evaluation-point search has been widened (7 → 25 for $\mathbb{Z}$,
8 → 32 for $\mathbb{F}_p$); very sparse or highly specialized polynomials may
still need an extended range or additional Diophantine interpolation passes.

These limitations are tracked in the project roadmap and will be lifted as the
algebra kernel matures.
