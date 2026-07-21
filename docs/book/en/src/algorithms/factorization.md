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
| $\mathbb{Z}[x,y]$ | Bivariate (monic in $x$) | Wang's Hensel lifting |
| $\mathbb{F}_p[x,y]$ | Bivariate (monic in $x$) | Hensel lifting over $\mathbb{F}_p$ |
| $\mathbb{Z}[x_1,\dots,x_n]$ | Multivariate | Wang EEZ Hensel lifting + leading-coefficient preprocessing + Zassenhaus recombination |
| $\mathbb{F}_p[x_1,\dots,x_n]$ | Multivariate | EEZ Hensel lifting (with characteristic-$p$ $p$-th power handling) |

Since 0.16.0, multivariate factorization with more than two variables is
supported. Non-monic univariate polynomials are factored via a leading-
coefficient transformation. The multivariate path conservatively reports
polynomials whose non-constant leading coefficient cannot be imposed as
irreducible for now (see "Limitations and Future Work").

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

## Limitations and Future Work

- The multivariate path conservatively reports inputs whose non-constant
leading coefficient cannot be imposed (requiring a mod-$p$ Hensel lift to
impose true leading coefficients) as irreducible; this enhancement is planned
for 0.16.1.
- The bivariate Wang Hensel limitation that the leading coefficient be a
constant is subsumed by the 0.16.0 arbitrary-multivariate (EEZ) path.
- The evaluation-point search is bounded; very sparse or highly specialized
polynomials may need an extended range or sparse Diophantine/interpolation in
the future.

These limitations are tracked in the project roadmap and will be lifted as the
algebra kernel matures.
