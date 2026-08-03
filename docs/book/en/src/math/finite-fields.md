# Foundations: Finite Fields and Modular Arithmetic

## Prerequisites

- Set theory and basic algebraic structures (definitions of groups, rings, and fields)
- Integer division and remainders
- Basic polynomial operations

Suggested reading: [Polynomial Algebra](./polynomial-algebra.md), [Linear Algebra](./linear-algebra.md).

---

## Basic Concepts

### Modular Arithmetic

Given a positive integer $n \geq 2$, the integers $a$ and $b$ are **congruent modulo $n$** (written $a \equiv b \pmod{n}$) if and only if $n \mid (a - b)$. Modular arithmetic partitions the integers into $n$ equivalence classes:

$$
\mathbb{Z}/n\mathbb{Z} = \{\overline{0},\, \overline{1},\, \dots,\, \overline{n-1}\}
$$

where $\overline{a} = \{a + kn \mid k \in \mathbb{Z}\}$. Addition and multiplication are performed on representatives and the result is then reduced modulo $n$:

$$
\overline{a} + \overline{b} = \overline{a + b}, \quad \overline{a} \cdot \overline{b} = \overline{a \cdot b}
$$

This makes $\mathbb{Z}/n\mathbb{Z}$ a commutative ring.

### Modular Inverses and the Extended Euclid Algorithm

The element $\overline{a} \in \mathbb{Z}/n\mathbb{Z}$ has a multiplicative inverse if and only if $\gcd(a, n) = 1$. In that case there exist integers $x, y$ satisfying the Bézout identity:

$$
ax + ny = \gcd(a, n) = 1
$$

Reducing modulo $n$ gives $ax \equiv 1 \pmod{n}$, hence $\overline{x} = \overline{a}^{-1}$.

The **extended Euclid algorithm** solves this identity iteratively:

1. Initialize: $(r_0, r_1) = (a, n)$, $(s_0, s_1) = (1, 0)$, $(t_0, t_1) = (0, 1)$.
2. At each step compute $q = \lfloor r_{i-1}/r_i \rfloor$ and update:

$$
\begin{pmatrix} r_{i+1} \\ s_{i+1} \\ t_{i+1} \end{pmatrix} = \begin{pmatrix} r_{i-1} - q \cdot r_i \\ s_{i-1} - q \cdot s_i \\ t_{i-1} - q \cdot t_i \end{pmatrix}
$$

3. Terminate when $r_{i+1} = 0$; then $r_i = \gcd(a, n)$ and $s_i \cdot a + t_i \cdot n = r_i$.

### Prime Fields $\mathbb{Z}/p\mathbb{Z}$

When $n = p$ is **prime**, the structure of $\mathbb{Z}/p\mathbb{Z}$ changes fundamentally: every nonzero element is coprime to $p$ and therefore has an inverse. This makes $\mathbb{Z}/p\mathbb{Z}$ a **field**, denoted $\mathbb{F}_p$.

**Theorem (Fermat's Little Theorem)**: let $p$ be prime and $a \not\equiv 0 \pmod{p}$; then

$$
a^{p-1} \equiv 1 \pmod{p}
$$

It follows immediately that $a^{p-2} \equiv a^{-1} \pmod{p}$, which provides an implementation of finite-field inversion that **does not need the extended Euclid algorithm** — a single modular exponentiation suffices. With the `pow_mod` of the GMP backend this is much faster than repeated division.

**Basic properties of fields**:

- The additive group $(\mathbb{F}_p, +)$ is cyclic of order $p$, isomorphic to $(\mathbb{Z}/p\mathbb{Z}, +)$.
- The multiplicative group $(\mathbb{F}_p^*, \times)$ is cyclic of order $p-1$ — there exists a **primitive root** $g$ such that $\{g^0, g^1, \dots, g^{p-2}\} = \mathbb{F}_p^*$.
- The **characteristic** is $p$: $p \cdot 1 = \underbrace{1 + \cdots + 1}_{p} = 0$, and $p$ is the smallest positive integer with this property.

---

## Core Theory

### Construction of Finite Fields $\mathbb{F}_{p^d}$

The prime field $\mathbb{F}_p$ is not the only finite field. The **classification theorem for finite fields** states that for every prime power $q = p^d$ ($d \geq 1$) there exists a unique (up to isomorphism) finite field $\mathbb{F}_q$ with $q$ elements, and every finite field is of this form.

The standard way to construct $\mathbb{F}_{p^d}$ ($d > 1$) is by **irreducible polynomial extension**:

1. Pick a degree-$d$ **irreducible polynomial** $m(x)$ in $\mathbb{F}_p[x]$.
2. Define $\mathbb{F}_{p^d} = \mathbb{F}_p[x] / (m(x))$ — the quotient of the polynomial ring by the ideal $(m(x))$.
3. Elements are polynomials $a_0 + a_1 \alpha + \cdots + a_{d-1}\alpha^{d-1}$ of degree $< d$, where $\alpha = \overline{x}$ is a root of $m$.
4. Addition is coefficient-wise modulo $p$; multiplication is polynomial multiplication followed by reduction modulo $m(x)$.

**Existence of irreducible polynomials**: the number of **monic** degree-$d$ irreducible polynomials in $\mathbb{F}_p[x]$ is

$$
N_p(d) = \frac{1}{d}\sum_{k \mid d} \mu(k) \cdot p^{d/k}
$$

where $\mu$ is the Möbius function. This guarantees that irreducible polynomials exist for every $d \geq 1$.

**Example**: $\mathbb{F}_4 = \mathbb{F}_2[x]/(x^2 + x + 1)$. Let $\alpha$ be a root of $x^2 + x + 1$; then $\alpha^2 = \alpha + 1$ (in $\mathbb{F}_2$, $-1 = 1$). The four elements are $\{0, 1, \alpha, \alpha+1\}$.

### Cyclic Multiplicative Groups

**Theorem**: the multiplicative group $\mathbb{F}_{p^d}^*$ of a finite field is cyclic, of order $p^d - 1$.

This means there exists a **primitive element** $g \in \mathbb{F}_{p^d}^*$ such that

$$
\mathbb{F}_{p^d}^* = \{g^0, g^1, \dots, g^{p^d - 2}\}
$$

Cyclicity has far-reaching applications in cryptography and coding theory: the discrete logarithm problem (DLP), pseudorandom number generation, and the construction of error-correcting codes all rely on it.

### Characteristic and Prime Subfield

The **characteristic** of the finite field $\mathbb{F}_{p^d}$ is $p$. From the definition of the characteristic:

- $p \cdot a = 0$ for all $a \in \mathbb{F}_{p^d}$.
- $\mathbb{F}_p$ is the **prime subfield** of $\mathbb{F}_{p^d}$, i.e. the smallest subfield containing $1$.
- As an $\mathbb{F}_p$-vector space, $\mathbb{F}_{p^d}$ has dimension $d$.

The Frobenius endomorphism $\varphi: x \mapsto x^p$ is an important structural map on $\mathbb{F}_{p^d}$. Its $d$-fold iterate satisfies $\varphi^d = \mathrm{id}$, and its fixed points are exactly the elements of $\mathbb{F}_p$.

---

## Implementation in oCAS

oCAS's finite-field implementation lives in `ocas-domain/src/finite_field.rs` and supports prime fields $\mathbb{Z}/p\mathbb{Z}$. Extension fields $\mathbb{F}_{p^d}$ are built via `AlgebraicExtension<FiniteField>` (see `ocas-domain/src/algebraic.rs`).

### Core Structs

**`FiniteField`** represents a prime field $\mathbb{F}_p$:

```rust
pub struct FiniteField {
    prime: BigInt,
    prime_minus_two: BigInt,  // cached p-2, used for Fermat inversion
    #[cfg(feature = "gmp")]
    prime_gmp: rug::Integer,  // GMP backend to accelerate modular arithmetic
    #[cfg(feature = "gmp")]
    prime_minus_two_gmp: rug::Integer,
}
```

- The `prime` field stores the modulus $p$ as a `num-bigint` arbitrary-precision integer.
- `prime_minus_two` is precomputed at construction time so it does not have to be recomputed on every inversion.
- With the `gmp` feature enabled, the GMP representations are cached as well, exploiting `rug`'s native modular exponentiation.

**`FiniteFieldElement`** represents an element of the field:

```rust
pub struct FiniteFieldElement {
    value: BigInt,  // always kept in the range [0, p-1]
}
```

Elements are always stored as their **canonical representative** — a nonnegative integer in the range $[0, p-1]$. This guarantees that semantic equality coincides with structural equality (`PartialEq` compares the values directly) without reducing modulo $p$ on every comparison.

### Construction and Element Creation

```rust
// create F_7
let f = FiniteField::new(BigInt::from(7));

// create an element from an arbitrary integer (automatically reduced to [0, p-1])
let a = f.element(10);   // a = 10 mod 7 = 3
let b = f.element(-3);   // b = -3 mod 7 = 4
```

The `element()` method uses `mod_floor` (rather than plain `%`), ensuring that negative numbers are also mapped correctly into $[0, p-1]$.

### Arithmetic Operations

All operations go through the uniform interface of the `Domain` trait:

| Operation | Method | Implementation strategy |
|---|---|---|
| Addition | `f.add(&a, &b)` | $(a + b) \bmod p$ |
| Subtraction | `f.sub(&a, &b)` | $(a - b) \bmod p$ (`mod_floor` handles negatives) |
| Multiplication | `f.mul(&a, &b)` | $(a \cdot b) \bmod p$ |
| Negation | `f.neg(&a)` | $(-a) \bmod p$ |
| Inversion | `f.inv(&a)` | $a^{p-2} \bmod p$ (Fermat's little theorem); `None` when $a = 0$ |
| Division | `f.div(&a, &b)` | $a \cdot b^{-1}$; `None` when $b = 0$ |
| Exponentiation | `f.pow(&a, n)` | fast modular exponentiation (binary exponentiation) |

**Inversion implementation details**:

```rust
fn inv(&self, a: &FiniteFieldElement) -> Option<FiniteFieldElement> {
    if a.value.is_zero() {
        return None;  // zero has no inverse
    }
    // Fermat: a^(p-2) ≡ a^{-1} (mod p)
    Some(self.normalize(a.value.modpow(&self.prime_minus_two, &self.prime)))
}
```

Fermat's little theorem is used instead of the extended Euclid algorithm because binary exponentiation via `modpow` needs only $O(\log p)$ multiplications and is simple to implement. With the GMP backend enabled, `rug::Integer::pow_mod` uses GMP's highly optimized implementation underneath.

### The EuclideanDomain Implementation

Finite fields also implement the `EuclideanDomain` trait, but Euclidean division over a field is degenerate:

```rust
// division over a field is always exact; the remainder is always zero
fn div_rem(&self, a, b) -> Option<(Self::Element, Self::Element)> {
    self.div(a, b).map(|q| (q, self.zero()))
}

// GCD over a field: returns 0 when both are zero, otherwise 1
fn gcd(&self, a, b) -> Self::Element {
    if self.is_zero(a) && self.is_zero(b) {
        self.zero()
    } else {
        self.one()
    }
}
```

This is because every nonzero element of a field is a **unit**, so the notion of GCD degenerates to a trivial case.

### Extension Fields $\mathbb{F}_{p^d}$

Constructed via `AlgebraicExtension<FiniteField>`:

```rust
use ocas_domain::{AlgebraicExtension, Domain, FiniteField};
use num_bigint::BigInt;

// construct F_4 = F_2[x]/(x^2 + x + 1)
let f2 = FiniteField::new(BigInt::from(2));
let f4 = AlgebraicExtension::new(
    f2,
    vec![
        f2.element(1),  // constant term
        f2.element(1),  // x coefficient
        f2.element(1),  // x^2 coefficient (leading term)
    ],
);
let alpha = f4.alpha();  // α = x mod (x^2+x+1)
// α^2 = α + 1 (in F_2)
let alpha_sq = f4.mul(&alpha, &alpha);
```

Internally, `AlgebraicElement<E>` stores a vector of coefficients in ascending powers, with trailing zeros trimmed:

```rust
pub struct AlgebraicElement<E> {
    coeffs: Vec<E>,  // [a_0, a_1, ..., a_{d-1}] represents a_0 + a_1*α + ...
}
```

Inversion uses the **extended Euclid algorithm** (a self-contained implementation over dense polynomials of the base field): it computes $s(x) \cdot a(x) + t(x) \cdot m(x) = 1$, so that $s(x) \bmod m(x) = a(x)^{-1}$.

---

## Advanced Topics

### Number-Theoretic Transforms (NTT) and Fast Polynomial Multiplication

The NTT is the finite-field analogue of the FFT. Standard FFT uses the roots of unity $\omega_n = e^{2\pi i/n}$ over the complex numbers $\mathbb{C}$; the NTT instead uses **$n$-th roots of unity modulo $p$** in $\mathbb{F}_p$.

**Existence condition**: $\mathbb{F}_p$ contains a primitive $n$-th root of unity if and only if $n \mid (p - 1)$. This is because $\mathbb{F}_p^*$ is a cyclic group of order $p-1$, and an $n$-th root of unity exists $\iff$ $n$ divides the group order.

**Algorithm** (given $f, g \in \mathbb{F}_p[x]$, each of degree $< n$):

1. **Zero padding**: pad the coefficient vectors of $f, g$ to length $N = 2^{\lceil \log_2(2n) \rceil}$.
2. **Forward transform**: $\hat{f} = \mathrm{NTT}(f)$, $\hat{g} = \mathrm{NTT}(g)$.
3. **Pointwise multiplication**: $\hat{h}_i = \hat{f}_i \cdot \hat{g}_i$ (component-wise modulo $p$).
4. **Inverse transform**: $h = \mathrm{NTT}^{-1}(\hat{h})$, then multiply by $N^{-1} \bmod p$.

Total complexity is $O(N \log N)$, compared with $O(n^2)$ for naive convolution and $O(n^{1.585})$ for Karatsuba.

**The NTT implementation in oCAS** (`ocas-poly/src/ntt.rs`):

- **Algorithm**: Cooley–Tukey radix-2 DIT (decimation in time) with bit-reversal permutation.
- **Montgomery arithmetic**: all intermediate values are stored in Montgomery form $x \cdot R \bmod p$ ($R = 2^{64}$) to avoid expensive `u128 % p` operations. One conversion is performed at the entry and one at the exit.
- **Trigger condition**: activated automatically when the number of coefficients is $\geq 256$ (`NTT_THRESHOLD`, i.e. degree $\geq 255$) and $p$ is NTT-friendly ($p-1$ contains a sufficiently large power of 2); otherwise it falls back to Karatsuba.
- **Limitation**: currently only primes representable in a `u64` are supported ($p < 2^{64}$, checked in `try_ntt_mul_fp`); since Montgomery reduction uses 128-bit intermediates, $p \leq 2^{63}$ guarantees no overflow.

```rust
// NTT-friendly check
pub fn is_ntt_friendly(p: u64, n: usize) -> bool {
    // n must divide p-1
    if n == 0 { return true; }
    let pm1 = p - 1;
    pm1 % (n as u64) == 0
}
```

**Construction of NTT-friendly primes**: in practice one uses primes of the form $p = k \cdot 2^m + 1$ (e.g. the NTT prime $998244353 = 119 \cdot 2^{23} + 1$), which guarantees $2^{23} \mid (p-1)$ and hence supports transforms of length up to $2^{23}$.

**Montgomery multiplication** replaces modular reduction with multiplications and shifts:

$$
\mathrm{MontMul}(a, b) = a \cdot b \cdot R^{-1} \bmod p
$$

where $R = 2^{64}$, and $R^{-1} \bmod p$ together with $p' = -p^{-1} \bmod R$ are precomputed. The multiplication steps are:

1. $t = a \cdot b$ (128-bit intermediate value)
2. $m = (t \bmod R) \cdot p' \bmod R$
3. $u = (t + m \cdot p) / R$
4. if $u \geq p$ then $u = u - p$

Only additions, multiplications, and shifts are used throughout — no expensive division.

### NTT in Polynomial GCD

The NTT does not only accelerate multiplication; it also indirectly speeds up polynomial GCD and factorization. During Hensel lifting, multiplying the many small-coefficient polynomials is the bottleneck, and the NTT reduces it from $O(n^2)$ to $O(n \log n)$. oCAS's `DenseUnivariatePolynomial<FiniteField>::mul_ntt` automatically selects the NTT when the number of coefficients is $\geq 256$ and the prime is NTT-friendly (otherwise it falls back to Karatsuba/Schoolbook); the generic `mul` still uses the Karatsuba/Schoolbook path.

---

## References

1. **Shoup, V.** *A Computational Introduction to Number Theory and Algebra.* Cambridge University Press, 2nd edition, 2009. — Chapters 4–5 cover finite-field construction, cyclic groups, and irreducible polynomials.
2. **Lidl, R. & Niederreiter, H.** *Introduction to Finite Fields and Their Applications.* Cambridge University Press, 1994. — The classical textbook on finite-field theory.
3. **Gathen, J. von zur & Gerhard, J.** *Modern Computer Algebra.* Cambridge University Press, 3rd edition, 2013. — Chapter 8 (fast polynomial multiplication and the NTT) and Chapter 14 (finite-field arithmetic).
4. **Menezes, A., van Oorschot, P. & Vanstone, S.** *Handbook of Applied Cryptography.* CRC Press, 1996. — Chapter 2 covers modular arithmetic and the extended Euclid algorithm.
