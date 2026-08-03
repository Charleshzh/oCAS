# Factorization

This chapter documents the complete polynomial factorization API in oCAS, covering factorization over the univariate rings $\mathbb{Z}[x]$ and $\mathbb{F}_p[x]$, the multivariate rings $\mathbb{Z}[x_1,\dots,x_n]$ and $\mathbb{F}_p[x_1,\dots,x_n]$, and algebraic number fields $\mathbb{Q}(\alpha)[x]$, as well as rational function arithmetic and resultants.

Related modules:

| Module path | Description |
|---|---|
| `ocas_poly::factor` | Top-level factorization entry point |
| `ocas_poly::factor::hensel` | Hensel lifting and Zassenhaus recombination ($\mathbb{Z}[x]$) |
| `ocas_poly::factor::finite_field` | $\mathbb{F}_p[x]$ factorization (Cantor–Zassenhaus + Berlekamp) |
| `ocas_poly::factor::multivariate` | Bivariate Hensel lifting (Wang's algorithm) |
| `ocas_poly::factor::eez` | Multivariate EEZ Hensel lifting (Wang + leading-coefficient reconstruction) |
| `ocas_poly::factor::algebraic` | Algebraic number field factorization (Trager's algorithm) |
| `ocas_poly::gcd` | GCD (pseudo-remainder + Euclidean algorithm) |
| `ocas_poly::gcd::modular` | Modular GCD (Brown 1971) |
| `ocas_poly::resultant` | Resultant (Brown PRS) |
| `ocas_poly::rational` | Rational function fraction field |
| `ocas_calc::partial_fraction` | Partial fraction decomposition |

---

## Type Aliases

```rust
/// Square-free factorization result: list of (factor, multiplicity) pairs.
pub type SquareFreeFactors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;

/// Complete factorization result: list of (irreducible factor, multiplicity) pairs.
pub type Factors<D> = Vec<(DenseUnivariatePolynomial<D>, usize)>;
```

---

## 1. Univariate Integer Polynomial Factorization ($\mathbb{Z}[x]$)

### DenseUnivariatePolynomial::factor

**Signature**: `pub fn factor(&self) -> Factors<IntegerDomain>`

**Description**: Factors a univariate integer polynomial into irreducible factors (primitive polynomials that are irreducible over $\mathbb{Q}$) together with their multiplicities.

**Parameters**: none (`self` is the polynomial to factor).

**Returns**: `Vec<(DenseUnivariatePolynomial<IntegerDomain>, usize)>` — a list of irreducible factor and multiplicity pairs. Each factor is a primitive polynomial with a positive leading coefficient. The product of all factors (each raised to its multiplicity) equals `primitive_part(self)`; the input must be primitive (content = 1), otherwise the content is discarded.

**Algorithm**:
1. The input must be primitive (content = 1; preprocess with `primitive_part()` if needed).
2. Square-free factorization (Yun's algorithm): use $\gcd(f, f')$ to separate multiplicities.
3. Call `factor_square_free` on each square-free component (internally it handles non-monic inputs via the leading-coefficient transformation, then calls `factor_square_free_monic`):
   - Choose a prime $p$ ($p \nmid \mathrm{lc}(f)$ and $f \bmod p$ square-free).
   - Factor into monic irreducible factors in $\mathbb{F}_p[x]$ (Cantor–Zassenhaus).
   - Compute the Mignotte bound $B = 2^n \|f\|_2$.
   - Linear Hensel lifting $p \to p^k$ ($p^k > 2B$).
   - Zassenhaus subset recombination: enumerate subsets of the lifted factors and verify by trial division.
4. Non-monic polynomials are handled via the leading-coefficient transformation $a^{d-1} f(x/a)$.

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
// Each factor has degree 1 and multiplicity 1
for (g, m) in &factors {
    assert_eq!(g.degree(), Some(1));
    assert_eq!(*m, 1);
}
// Output: [(x - 1, 1), (x + 1, 1)]
```

**See also**: [`square_free_factorization`](#square_free_factorization), [`factor_over_finite_field`](#factor_over_finite_field), [Polynomial GCD and Factorization](../math/poly-gcd-factoring.md)

---

### DenseUnivariatePolynomial::square_free_factorization

**Signature**: `pub fn square_free_factorization(&self) -> SquareFreeFactors<D>`

**Description**: Computes the square-free factorization of the polynomial. Returns the distinct square-free factors together with their multiplicities in the original polynomial.

**Parameters**: none (`self` is the polynomial to factor).

**Returns**: `Vec<(DenseUnivariatePolynomial<D>, usize)>` — a list of square-free factor and multiplicity pairs.

**Algorithm**: Yun's algorithm (characteristic 0). Let $g = \gcd(f, f')$, $w = f/g$, then iterate $h = \gcd(w, g)$, $z = w/h$ to collect the factors of multiplicity $k$.

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// (x+1)^2 * (x-1) = x^3 + x^2 - x - 1
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(-1), Integer::from(1), Integer::from(1),
]);
let factors = p.square_free_factorization();
assert_eq!(factors.len(), 2);
// Output: [(x - 1, 1), (x + 1, 2)] (Yun's algorithm collects the multiplicity-1 factor first)
```

**See also**: [`is_square_free`](#is_square_free), [`factor`](#denseunivariatepolynomialfactor)

---

### DenseUnivariatePolynomial::is_square_free

**Signature**: `pub fn is_square_free(&self) -> bool`

**Description**: Tests whether the polynomial is square-free (i.e., $\gcd(f, f') = 1$).

**Parameters**: none.

**Returns**: `bool` — `true` if the polynomial is square-free.

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x^2 - 1 is square-free
let p = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
assert!(p.is_square_free());

// (x+1)^2 is not square-free
let q = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
assert!(!q.is_square_free());
```

---

## 2. Univariate Finite Field Polynomial Factorization ($\mathbb{F}_p[x]$)

### factor_over_finite_field

**Signature**: `pub fn factor_over_finite_field(f: &FpPoly) -> Vec<(FpPoly, usize)>`

**Description**: Factors a polynomial in $\mathbb{F}_p[x]$ into monic irreducible factors with multiplicities.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&DenseUnivariatePolynomial<FiniteField>` | The polynomial to factor |

**Returns**: `Vec<(FpPoly, usize)>` — a list of monic irreducible factor and multiplicity pairs. The leading coefficient is listed separately as a constant factor (multiplicity 1).

**Algorithm**:
1. Square-free factorization (Musser/Bernardin algorithm, handling $p$-th roots in characteristic $p$).
2. For each square-free component:
   - **DDF** (distinct-degree factorization): groups factors by degree using the Frobenius map $x \mapsto x^{p^d} \bmod f$.
   - **EDF** (equal-degree factorization): for odd characteristic $p$, pick a random $a$ and compute $\gcd(f, a^{(p^d-1)/2} - 1)$; for characteristic 2, use the trace map.
3. For small primes, Berlekamp's algorithm (the kernel of the Frobenius matrix $Q^T - I$) can also be used.

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::factor::finite_field::factor_over_finite_field;

let f = FiniteField::new(BigInt::from(5));
// x^2 - 1 = (x-1)(x+1) over F_5
let p = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(4), f.element(0), f.element(1)]);
let factors = factor_over_finite_field(&p);
let linear_count = factors.iter()
    .filter(|(g, _)| g.degree() == Some(1)).count();
assert_eq!(linear_count, 2);
// Output: two linear monic factors
```

**See also**: [Cantor–Zassenhaus algorithm details](../algorithms/factorization.md)

---

### DenseUnivariatePolynomial::factor (FiniteField)

**Signature**: `pub fn factor(&self) -> Factors<FiniteField>`

**Description**: The `factor()` method on `DenseUnivariatePolynomial<FiniteField>`, equivalent to calling `factor_over_finite_field`.

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

---

### berlekamp

**Signature**: `pub fn berlekamp(f: &FpPoly) -> Vec<FpPoly>`

**Description**: Berlekamp's factorization algorithm, suitable for small prime fields. Builds the Frobenius matrix $Q$ (where $Q[i][j]$ is the coefficient of $x^j$ in $x^{ip} \bmod f$) and computes the kernel of $Q^T - I$; every nonzero kernel vector $v$ satisfies $v^p \equiv v \pmod{f}$ and splits factors via $\gcd(f, v - a)$ for $a \in \mathbb{F}_p$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&FpPoly` | A monic square-free polynomial |

**Returns**: `Vec<FpPoly>` — a list of monic irreducible factors.

**See also**: Berlekamp (1970)

---

### cantor_zassenhaus

**Signature**: `pub fn cantor_zassenhaus(f: &FpPoly) -> Vec<FpPoly>`

**Description**: Cantor–Zassenhaus factorization algorithm. Runs DDF (distinct-degree factorization) first, then EDF (equal-degree factorization). Returns the list of monic irreducible factors of a monic square-free polynomial.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&FpPoly` | A monic square-free polynomial |

**Returns**: `Vec<FpPoly>` — a list of monic irreducible factors.

**See also**: Cantor & Zassenhaus (1981)

---

### poly_pow_mod

**Signature**: `pub fn poly_pow_mod(base: &FpPoly, exp: &BigInt, modulus: &FpPoly) -> FpPoly`

**Description**: Computes $\text{base}^{\text{exp}} \bmod \text{modulus}$ using fast exponentiation (repeated squaring), reducing after each multiplication.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `base` | `&FpPoly` | The base polynomial |
| `exp` | `&BigInt` | The exponent (non-negative) |
| `modulus` | `&FpPoly` | The modulus polynomial (nonzero) |

**Returns**: `FpPoly` — the result of the modular exponentiation.

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::factor::finite_field::poly_pow_mod;

let f = FiniteField::new(BigInt::from(7));
let m = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(1), f.element(0), f.element(1)]); // x^2 + 1
let base = DenseUnivariatePolynomial::from_coeffs(
    f.clone(), vec![f.element(0), f.element(1)]); // x
// x^2 mod (x^2+1) = -1 = 6 in F_7
let r = poly_pow_mod(&base, &BigInt::from(2), &m);
assert_eq!(r.coeff(0).cloned(), Some(f.element(6)));
// Output: 6 (i.e., -1 mod 7)
```

---

## 3. Multivariate Integer Polynomial Factorization ($\mathbb{Z}[x_1,\dots,x_n]$)

### SparseMultivariatePolynomial::factor (IntegerDomain, Lex)

**Signature**: `pub fn factor(&self) -> Vec<(Self, usize)>`

**Description**: Factors a sparse multivariate integer polynomial into irreducible factors with multiplicities.

**Parameters**: none (`self` is the `SparseMultivariatePolynomial<IntegerDomain, Lex>` to factor).

**Returns**: `Vec<(SparseMultivariatePolynomial<IntegerDomain, Lex>, usize)>` — a list of irreducible factor and multiplicity pairs.

**Algorithm selection**:
- **$n \geq 3$ variables**: uses EEZ Hensel lifting (Wang leading-coefficient reconstruction + p-adic coefficient Hensel lifting + Zassenhaus recombination).
- **Bivariate with constant leading coefficient**: uses bivariate Hensel lifting (evaluation–lifting path).
- **Bivariate with non-constant leading coefficient**: falls back to the EEZ path.

**Algorithm** (EEZ path, $n \geq 3$):
1. Extract the content in the main variable and compute the square-free factorization.
2. Pick sample points $(a_1, \dots, a_n)$ such that the univariate image $f(x_0, a_1, \dots)$ keeps its degree and is square-free.
3. Factor the univariate image over $\mathbb{Z}$.
4. Wang leading-coefficient reconstruction: distribute the irreducible factors of the total leading coefficient among the univariate factors.
5. EEZ variable-by-variable Hensel lifting: recover $x_1, x_2, \dots$ in turn, solving a multivariate Diophantine equation at each step.
6. p-adic coefficient Hensel lifting: lift from modulus $p$ to a sufficiently large $p^k$.
7. Zassenhaus subset recombination: enumerate subsets of the modular factors and verify by trial division.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

// (x^2 + y + 1)(x + y + 2) in Z[x,y]
let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    IntegerDomain, 2,
    vec![
        (vec![3, 0], Integer::from(1)),   // x^3
        (vec![2, 1], Integer::from(1)),   // x^2 y
        (vec![2, 0], Integer::from(2)),   // 2x^2
        (vec![1, 1], Integer::from(1)),   // xy
        (vec![1, 0], Integer::from(1)),   // x
        (vec![0, 2], Integer::from(1)),   // y^2
        (vec![0, 1], Integer::from(3)),   // 3y
        (vec![0, 0], Integer::from(2)),   // 2
    ],
);
let factors = f.factor();
assert!(factors.len() >= 2);
// Output: [(x^2 + y + 1, 1), (x + y + 2, 1)]
```

**See also**: [Wang EEZ algorithm](../math/poly-gcd-factoring.md), [`bivariate_factor_z`](#bivariate_factor_z)

---

### bivariate_factor_z

**Signature**: `pub fn bivariate_factor_z(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<(ZMPoly, usize)>`

**Description**: Factors a bivariate integer polynomial into irreducible factors with multiplicities. Requires the leading coefficient in $x$ to be a nonzero integer constant.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<IntegerDomain, Lex>` | The polynomial to factor |
| `x_var` | `usize` | Index of the main variable |
| `y_var` | `usize` | Index of the secondary variable |

**Returns**: `Vec<(ZMPoly, usize)>` — irreducible factor and multiplicity pairs.

**Algorithm**:
1. Bivariate square-free factorization (heuristic bivariate GCD).
2. Choose $y = \alpha$ so that the univariate image $f(x, \alpha)$ is square-free with the fewest factors.
3. Hensel lifting: lift the factors of $f(x, \alpha)$ back to bivariate factors, correcting iteratively via Taylor expansion.

**See also**: Wang (1978)

---

### multivariate_factor_z

**Signature**: `pub fn multivariate_factor_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)>`

**Description**: Entry point for multivariate factorization via EEZ Hensel lifting. Supports non-constant leading coefficients (through Wang leading-coefficient reconstruction).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<IntegerDomain, Lex>` | The polynomial to factor |

**Returns**: `Vec<(ZmPoly, usize)>` — irreducible factor and multiplicity pairs.

**See also**: [Wang EEZ algorithm details](../math/poly-gcd-factoring.md)

---

## 4. Multivariate Finite Field Polynomial Factorization ($\mathbb{F}_p[x_1,\dots,x_n]$)

### SparseMultivariatePolynomial::factor (FiniteField, Lex)

**Signature**: `pub fn factor(&self) -> Vec<(Self, usize)>`

**Description**: Factors a sparse multivariate finite field polynomial into irreducible factors with multiplicities.

**Parameters**: none.

**Returns**: `Vec<(SparseMultivariatePolynomial<FiniteField, Lex>, usize)>` — a list of irreducible factor and multiplicity pairs.

**Algorithm selection**:
- **$n \geq 3$ variables**: EEZ Hensel lifting.
- **Bivariate**: evaluation–Hensel path.

**Example**:

```rust
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::Lex;

let fp = FiniteField::new(BigInt::from(7));
// x*y + 1 over F_7[x,y] — already irreducible
let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    fp.clone(), 2,
    vec![
        (vec![1, 1], fp.element(1)),  // xy
        (vec![0, 0], fp.element(1)),  // 1
    ],
);
let factors = f.factor();
assert_eq!(factors.len(), 1);
assert_eq!(factors[0].1, 1); // multiplicity 1
```

---

### bivariate_factor_fp

**Signature**: `pub fn bivariate_factor_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<(FpMPoly, usize)>`

**Description**: Factors a bivariate $\mathbb{F}_p$ polynomial into irreducible factors with multiplicities.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<FiniteField, Lex>` | The polynomial to factor |
| `x_var` | `usize` | Index of the main variable |
| `y_var` | `usize` | Index of the secondary variable |

**Returns**: `Vec<(FpMPoly, usize)>` — irreducible factor and multiplicity pairs.

---

### multivariate_factor_fp

**Signature**: `pub fn multivariate_factor_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)>`

**Description**: Entry point for multivariate $\mathbb{F}_p$ factorization, supporting non-constant leading coefficients (Wang leading-coefficient reconstruction + EEZ lifting).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `f` | `&SparseMultivariatePolynomial<FiniteField, Lex>` | The polynomial to factor |

**Returns**: `Vec<(FpMPoly, usize)>` — irreducible factor and multiplicity pairs.

---

## 5. Algebraic Number Field Factorization ($\mathbb{Q}(\alpha)[x]$)

### DenseUnivariatePolynomial::factor (AlgebraicNumberField)

**Signature**: `pub fn factor(&self) -> Factors<AlgebraicNumberField>`

**Description**: Factors a polynomial over the algebraic number field $\mathbb{Q}(\alpha)$ into monic irreducible factors with multiplicities (Trager's algorithm).

**Parameters**: none (`self` is a `DenseUnivariatePolynomial<AlgebraicNumberField>`).

**Returns**: `Vec<(DenseUnivariatePolynomial<AlgebraicNumberField>, usize)>` — monic irreducible factor and multiplicity pairs. The product of all factors (each raised to its multiplicity) equals `self` divided by its leading coefficient (a unit of $K$).

**Algorithm** (Trager):
1. Square-free factorization (Yun's algorithm, using the modular GCD `gcd_anf` to avoid coefficient blow-up).
2. Call `factor_square_free_anf` on each square-free component:
   - **Trager shift**: find $s \geq 0$ such that the norm of $f(x - s\alpha)$ is square-free.
   - **Norm computation**: $\operatorname{Res}_\alpha(m(\alpha), f(x, \alpha))$, via evaluation–interpolation.
   - **Factor the norm over $\mathbb{Q}$**: using the Hensel path.
   - **Recover the factors over $K$**: compute $\gcd_K(f, g_i(\alpha))$ for each norm factor $g_i$.
3. GCD uses the modular method (Encarnación): map to $\mathrm{GF}(p^d)$, combine via CRT, rational reconstruction, and trial-division verification.

**Example**:

```rust
use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

// Construct Q(√2) with minimal polynomial x^2 - 2
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
// x^2 - 2 factors over Q(√2) as (x - √2)(x + √2)
let f = DenseUnivariatePolynomial::from_coeffs(
    field.clone(),
    vec![
        field.from_base(Rational::new(-2, 1)),
        field.zero(),
        field.one(),
    ],
);
let factors = f.factor();
assert_eq!(factors.len(), 2);
assert!(factors.iter().all(|(g, m)| *m == 1 && g.degree() == Some(1)));
// Output: [(x - √2, 1), (x + √2, 1)]
```

**See also**: [Algebraic Number Fields and Galois Theory](../math/algebraic-number-fields.md), [Trager's algorithm](../math/poly-gcd-factoring.md)

---

### norm_with_shift

**Signature**: `pub(crate) fn norm_with_shift(field: &AlgebraicNumberField, f: &UP<AlgebraicNumberField>) -> Option<(u64, UP<AlgebraicNumberField>, UP<RationalDomain>)>`

**Description**: Trager shift: find $s \geq 0$ such that the norm of $f(x - s\alpha)$ is square-free over $\mathbb{Q}$. Returns `(s, g, norm)`. Tries at most `MAX_TRAGER_SHIFTS` (100) shifts.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `field` | `&AlgebraicNumberField` | The algebraic number field $\mathbb{Q}(\alpha)$ |
| `f` | `&DenseUnivariatePolynomial<AlgebraicNumberField>` | The polynomial to factor |

**Returns**: `Option<(u64, UP<AlgebraicNumberField>, UP<RationalDomain>)>` — `(s, f(x - s\alpha), \operatorname{Norm}(f(x - s\alpha)))`.

---

### gcd_anf

**Signature**: `pub(crate) fn gcd_anf(field: &AlgebraicNumberField, a: &UP<AlgebraicNumberField>, b: &UP<AlgebraicNumberField>) -> UP<AlgebraicNumberField>`

**Description**: GCD of two univariate polynomials over an algebraic number field. Uses the modular method (Encarnación): map to $\mathrm{GF}(p^d)$, combine the monic modular GCDs (CRT), rationally reconstruct the coefficients, and verify by trial division. Uses at most `MAX_ANF_GCD_PRIMES` (1000) primes, falling back to a dense Euclidean GCD beyond that.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `field` | `&AlgebraicNumberField` | The algebraic number field |
| `a` | `&UP<AlgebraicNumberField>` | The first polynomial |
| `b` | `&UP<AlgebraicNumberField>` | The second polynomial |

**Returns**: a monic GCD polynomial.

---

## 6. Polynomial GCD

### DenseUnivariatePolynomial::gcd

**Signature**: `pub fn gcd(&self, other: &Self) -> Self`

**Description**: Computes the greatest common divisor of two univariate polynomials. Uses the pseudo-remainder Euclidean algorithm for non-field coefficients. The result is a primitive polynomial.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Self` | The other polynomial |

**Returns**: a primitive GCD polynomial.

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x^2 - 1 = (x-1)(x+1)
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
// x^2 + 2x + 1 = (x+1)^2
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(1), Integer::from(2), Integer::from(1),
]);
let g = a.gcd(&b);
assert_eq!(g.coeffs(), &[Integer::from(1), Integer::from(1)]);
// Output: x + 1
```

---

### gcd_modular_z

**Signature**: `pub fn gcd_modular_z(a: &ZPoly, b: &ZPoly) -> ZPoly`

**Description**: Brown's modular GCD algorithm. Computes monic GCDs in $\mathbb{F}_p[x]$ for several primes $p$, combines them via CRT into a primitive GCD in $\mathbb{Z}[x]$, and verifies by trial division. It is much more efficient than pseudo-remainder GCD for large coefficients.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `a` | `&DenseUnivariatePolynomial<IntegerDomain>` | The first polynomial |
| `b` | `&DenseUnivariatePolynomial<IntegerDomain>` | The second polynomial |

**Returns**: a primitive GCD polynomial.

**Algorithm details**:
1. Choose primes $p > 2^{30}$ (to avoid small-prime issues).
2. Compute the monic GCD in $\mathbb{F}_p[x]$, scaled by $\gamma = \gcd(\mathrm{lc}(a), \mathrm{lc}(b))$.
3. Combine via CRT (symmetric representatives), discarding "unlucky primes" (where the GCD degree is higher than the true value).
4. Verify by exact trial division. Tries at most `MAX_PRIMES` (10000) primes, falling back to pseudo-remainder GCD beyond that.

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::gcd::modular::gcd_modular_z;

let d = IntegerDomain;
let i = |v: i64| Integer::from(v);
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);
let g = gcd_modular_z(&a, &b);
assert_eq!(g.coeffs(), &[i(1), i(1)]); // x + 1
// Output: x + 1
```

**See also**: [Modular GCD algorithm (Brown 1971)](../math/poly-gcd-factoring.md)

---

### DenseUnivariatePolynomial::content

**Signature**: `pub fn content(&self) -> D::Element`

**Description**: Computes the content of the polynomial (the GCD of all coefficients). The content of the zero polynomial is zero.

---

### DenseUnivariatePolynomial::primitive_part

**Signature**: `pub fn primitive_part(&self) -> Self`

**Description**: Returns the primitive part (the polynomial divided by its content).

---

## 7. Resultants

### DenseUnivariatePolynomial::resultant

**Signature**: `pub fn resultant(&self, other: &Self) -> D::Element`

**Description**: Computes the resultant $\operatorname{Res}(a, b)$ of two univariate polynomials using the Brown PRS (polynomial remainder sequence) algorithm.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Self` | The other polynomial |

**Returns**: `D::Element` — the resultant scalar. It is zero if and only if $\gcd(a, b)$ is non-constant.

**Algorithm**: subresultant PRS (Brown), with an exact division by $\beta$ at each step (the subresultant theorem guarantees exactness in a UFD).

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
// Output: -1
```

**Example (with a common root)**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let d = IntegerDomain;
// x - 1 and x^2 - 1 share the root x = 1
let a = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(1),
]);
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Integer::from(-1), Integer::from(0), Integer::from(1),
]);
assert_eq!(a.resultant(&b), Integer::from(0));
// Output: 0 (they share the factor x - 1)
```

**See also**: [Resultants and subresultant PRS](../math/poly-gcd-factoring.md)

---

## 8. Rational Function Fraction Field

### RationalPolynomial

**Signature**:

```rust
pub struct RationalPolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    pub numerator: SparseMultivariatePolynomial<D, O>,
    pub denominator: SparseMultivariatePolynomial<D, O>,
}
```

**Description**: An element of the polynomial fraction field $\frac{n}{d}$, always kept in canonical form: coprime numerator and denominator, with the denominator's leading coefficient positive (ordered domains) or 1 (finite fields).

#### Construction

##### RationalPolynomial::new

```rust
pub fn new(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

Creates a rational polynomial. Does **not** reduce automatically — use `from_num_den` when a canonical form is needed.

##### RationalPolynomial::from_num_den

```rust
pub fn from_num_den(
    numerator: SparseMultivariatePolynomial<D, O>,
    denominator: SparseMultivariatePolynomial<D, O>,
) -> Self
```

Creates from numerator and denominator and reduces automatically (GCD simplification + denominator leading-coefficient normalization).

##### RationalPolynomial::from_polynomial

```rust
pub fn from_polynomial(poly: SparseMultivariatePolynomial<D, O>) -> Self
```

Creates from a polynomial (denominator = 1).

#### Query Methods

| Method | Signature | Description |
|---|---|---|
| `is_zero` | `&self -> bool` | Whether the numerator is zero |
| `is_one` | `&self -> bool` | Whether it is 1/1 |
| `n_vars` | `&self -> usize` | Number of variables |
| `domain` | `&self -> &D` | Coefficient domain reference |

#### Arithmetic Operations

| Method | Signature | Description |
|---|---|---|
| `add` | `(&self, &Self) -> Self` | Addition (denominator-GCD strategy to reduce intermediate swell) |
| `sub` | `(&self, &Self) -> Self` | Subtraction |
| `mul` | `(&self, &Self) -> Self` | Multiplication (cross-cancellation) |
| `div` | `(&self, &Self) -> Option<Self>` | Division (returns `None` when the divisor's numerator is zero) |
| `neg` | `&self -> Self` | Negation $-\frac{n}{d}$ |
| `inv` | `&self -> Option<Self>` | Inverse (returns `None` when the numerator is zero) |
| `pow` | `(&self, k: u32) -> Self` | Power $\left(\frac{n}{d}\right)^k$ |

**Example**:

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::rational::RationalPolynomial;
use ocas_poly::Grevlex;

let d = IntegerDomain;
let n_vars = 2;

// x / y
let num: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![1, 0], Integer::from(1))]);
let den: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![0, 1], Integer::from(1))]);
let f = RationalPolynomial::from_num_den(num, den);

// y / x
let num2: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![0, 1], Integer::from(1))]);
let den2: SparseMultivariatePolynomial<_, Grevlex> = SparseMultivariatePolynomial::from_terms(d, n_vars,
    vec![(vec![1, 0], Integer::from(1))]);
let g = RationalPolynomial::from_num_den(num2, den2);

// (x/y) * (y/x) = 1
let h = f.mul(&g);
assert!(h.is_one());
// Output: 1
```

**See also**: [RationalPolynomial definition](../api/rust-polynomials.md)

---

## 9. Partial Fraction Decomposition

### apart

**Signature**:

```rust
pub fn apart<D: EuclideanDomain>(
    num: &DenseUnivariatePolynomial<D>,
    den: &DenseUnivariatePolynomial<D>,
) -> (
    Option<DenseUnivariatePolynomial<D>>,
    Vec<PartialFractionTerm<D>>,
)
```

**Description**: Computes the partial fraction decomposition of $\frac{\text{num}}{\text{den}}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `num` | `&DenseUnivariatePolynomial<D>` | The numerator polynomial |
| `den` | `&DenseUnivariatePolynomial<D>` | The denominator polynomial (nonzero) |

**Returns**: `(Option<poly_part>, Vec<PartialFractionTerm<D>>)` — an optional polynomial part and a list of partial fraction terms, satisfying:

$$\frac{\text{num}}{\text{den}} = \text{poly\_part} + \sum_i \frac{\text{numer}_i}{\text{denom}_i^{\text{exp}_i}}$$

**Example**:

```rust
use ocas_domain::{RationalDomain, Rational};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::apart;

let d = RationalDomain;
// 1 / (x^2 - 1)
let num = DenseUnivariatePolynomial::from_coeffs(d, vec![Rational::new(1, 1)]);
let den = DenseUnivariatePolynomial::from_coeffs(d, vec![
    Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1),
]);
let (poly_part, terms) = apart(&num, &den);
assert!(poly_part.is_none()); // proper fraction, no polynomial part
// x^2 - 1 = (x-1)(x+1) is square-free → the output groups by factor
// Output: poly_part = None, terms = [...]
```

---

### PartialFractionTerm

**Signature**:

```rust
pub struct PartialFractionTerm<D: EuclideanDomain> {
    pub numer: DenseUnivariatePolynomial<D>,  // numerator polynomial
    pub denom: DenseUnivariatePolynomial<D>,  // irreducible (square-free) denominator factor
    pub exp: usize,                            // multiplicity of this factor in the original denominator
}
```

**Description**: A single term of a partial fraction decomposition, representing $\frac{\text{numer}}{\text{denom}^{\text{exp}}}$.

---

### together

**Signature**:

```rust
pub fn together<D: EuclideanDomain>(
    poly_part: Option<&DenseUnivariatePolynomial<D>>,
    terms: &[PartialFractionTerm<D>],
) -> (DenseUnivariatePolynomial<D>, DenseUnivariatePolynomial<D>)
```

**Description**: Combines the polynomial part and partial fraction terms into a single rational function $\frac{n}{d}$. The inverse of `apart`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `poly_part` | `Option<&DenseUnivariatePolynomial<D>>` | The polynomial part (optional) |
| `terms` | `&[PartialFractionTerm<D>]` | The list of partial fraction terms |

**Returns**: `(numerator, denominator)` — the combined numerator and denominator.

---

## 10. Univariate Polynomial Helper Methods

### DenseUnivariatePolynomial::extended_gcd_poly

**Signature**: `pub fn extended_gcd_poly(&self, other: &Self) -> (Self, Self, Self)`

**Description**: Extended Euclidean algorithm. Returns $(g, s, t)$ such that $s \cdot a + t \cdot b = g = \gcd(a, b)$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Self` | The other polynomial |

**Returns**: `(gcd, s, t)`.

---

### DenseUnivariatePolynomial::div_rem

**Signature**: `pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)>`

**Description**: Division with remainder. Returns `(quotient, remainder)` such that `self = quotient * divisor + remainder`. Returns `None` if the divisor is zero.

---

### DenseUnivariatePolynomial::pow

**Signature**: `pub fn pow(&self, n: u32) -> Self`

**Description**: Computes $\text{self}^n$ by repeated squaring.

---

### DenseUnivariatePolynomial::p_adic_expansion

**Signature**: `pub fn p_adic_expansion(&self, p: &Self) -> Vec<Self>`

**Description**: Computes the $p$-adic expansion of the polynomial with respect to the base polynomial $p$. Returns the coefficient list $[c_0, c_1, \dots]$ such that $\text{self} = \sum_i c_i \cdot p^i$, where $\deg(c_i) < \deg(p)$.

---

## 11. Monomial Utility Functions

### monomial_divides

**Signature**: `pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool`

**Description**: Checks whether the monomial $b$ divides $a$ (i.e., $a$ is a multiple of $b$: $a_i \geq b_i$ for all $i$).

---

### monomial_lcm

**Signature**: `pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]>`

**Description**: Computes the least common multiple of two monomials (component-wise maximum).

---

### monomial_are_coprime

**Signature**: `pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool`

**Description**: Checks whether two monomials are coprime (no variable appears in both).

---

## 12. Internal Algorithm Overview

### Hensel Lifting and Zassenhaus Recombination

File: `ocas-poly/src/factor/hensel.rs`

**mignotte_bound**

```rust
pub(crate) fn mignotte_bound(f: &ZPoly) -> Integer
```

Landau–Mignotte bound: for a degree-$n$ polynomial $f$, any factor $g$ satisfies $\|g\|_\infty \leq 2^n \|f\|_2$.

**factor_square_free_monic**

```rust
pub fn factor_square_free_monic(f: &ZPoly) -> Vec<ZPoly>
```

Factors a monic square-free $\mathbb{Z}[x]$ polynomial into monic irreducible factors.

**factor_square_free**

```rust
pub fn factor_square_free(f: &ZPoly) -> Vec<ZPoly>
```

Factors a square-free primitive $\mathbb{Z}[x]$ polynomial into irreducible factors. Non-monic input is handled via the leading-coefficient transformation $a^{d-1} f(x/a)$.

**factor_primitive**

```rust
pub fn factor_primitive(f: &ZPoly) -> Vec<(ZPoly, usize)>
```

Factors a primitive $\mathbb{Z}[x]$ polynomial into irreducible factors with multiplicities (Yun square-free factorization + `factor_square_free_monic`). This is the internal implementation of `DenseUnivariatePolynomial::factor()`.

---

### EEZ Hensel Lifting

File: `ocas-poly/src/factor/eez.rs`

Key internal functions:

| Function | Description |
|---|---|
| `eez_lift` | Generic EEZ lifting (over a field, monic) |
| `eez_lift_imposed` | Non-monic EEZ lifting (Wang leading-coefficient imposition) |
| `eez_lift_z` | Integer EEZ lifting (solve over $\mathbb{Q}$ + integrality check) |
| `coefficient_hensel_lift_z` | p-adic coefficient Hensel lifting |
| `diophantine` | Recursive multivariate Diophantine equation solver |
| `sparse_diophantine_fp` | Sparse Diophantine solver with skeleton interpolation |
| `wang_reconstruct_lcoeffs` | Wang leading-coefficient reconstruction |
| `zassenhaus_multivariate` | Multivariate Zassenhaus subset recombination |

---

### Finite Field Factorization

File: `ocas-poly/src/factor/finite_field.rs`

| Function | Description |
|---|---|
| `factor_over_finite_field` | Top-level entry: square-free + DDF + EDF |
| `distinct_degree_factorization` | DDF: group by the degree of the irreducible factors |
| `equal_degree_factorization` | EDF: randomly split the product of equal-degree factors |
| `berlekamp` | Berlekamp's algorithm (small primes) |
| `cantor_zassenhaus` | DDF + EDF combination |
| `square_free_factorization_ff` | Characteristic-$p$ square-free factorization (Musser/Bernardin) |
| `pth_root_prime` | $p$-th roots over $\mathbb{F}_p$ |

---

### Algebraic Number Field Factorization

File: `ocas-poly/src/factor/algebraic.rs`

| Function | Description |
|---|---|
| `factor_anf` | Top-level entry: square-free + Trager |
| `factor_square_free_anf` | Trager: norm → $\mathbb{Q}$ factorization → GCD recovery |
| `norm_with_shift` | Trager shift to find a square-free norm |
| `norm_eval_interp` | Norm computation by evaluation–interpolation |
| `gcd_anf` | Modular GCD over $\mathbb{Q}(\alpha)$ |
| `square_free_anf` | Yun square-free factorization (modular GCD) |
| `factor_square_free_rationals` | Square-free factorization over $\mathbb{Q}$ (clear denominators → $\mathbb{Z}$ Hensel) |

---

## See also

- [Polynomial GCD and Factorization](../math/poly-gcd-factoring.md) — mathematical theory and algorithm details
- [Algebraic Number Fields and Galois Theory](../math/algebraic-number-fields.md) — mathematical background of Trager's algorithm
- [Polynomial API](./rust-polynomials.md) — complete API of `DenseUnivariatePolynomial` and `SparseMultivariatePolynomial`
- [Coefficient Domain API](./rust-domains.md) — definitions of `Integer`, `FiniteField`, `AlgebraicNumberField`
