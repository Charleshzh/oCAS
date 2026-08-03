# Advanced: Polynomial GCD and Factorization

This chapter systematically explains the core algorithms for computing polynomial greatest common divisors (GCDs) and factorization. These algorithms are the foundational building blocks of a computer algebra system — Gröbner basis computation, rational function simplification, and the partial fraction decomposition used in ODE solving all depend on them. We start from the classical Euclidean algorithm, gradually introduce modular methods, Hensel lifting, and multivariate extensions, and finally cover factorization over finite fields and algebraic number fields.

## Prerequisites

Before reading this chapter, we recommend studying:

- [Basic: Polynomial Algebra](./polynomial-algebra.md) — polynomial rings, degrees, monomial orders
- [Basic: Finite Fields and Modular Arithmetic](./finite-fields.md) — construction of $\mathbb{F}_p$, modular inverse, CRT
- [Basic: Linear Algebra](./linear-algebra.md) — matrix operations, determinants

### Review of the Euclidean Algorithm

Let $\mathbb{F}$ be a field and $f, g \in \mathbb{F}[x]$ with $g \ne 0$. The Euclidean algorithm computes the GCD by repeated division with remainder:

$$f = q_1 g + r_1, \quad g = q_2 r_1 + r_2, \quad r_1 = q_3 r_2 + r_3, \quad \ldots, \quad r_{k-2} = q_k r_{k-1} + r_k$$

The last nonzero remainder $r_k$ is $\gcd(f, g)$ (up to a constant factor).

In $\mathbb{Z}[x]$, the coefficients are not field elements, so division with remainder may not be possible. We need **pseudo-division**:

$$\text{lc}(g)^{d+1} \cdot f = q \cdot g + r, \quad d = \deg(f) - \deg(g)$$

where $\text{lc}(g)$ is the leading coefficient of $g$. This guarantees that the quotient and remainder lie in $\mathbb{Z}[x]$, at the cost of possible coefficient growth in the remainders.

### Pseudo-Remainder Sequences (PRS)

Applying the Euclidean algorithm directly to $\mathbb{Z}[x]$ with pseudo-division leads to **coefficient growth**: at degrees above about 16, the coefficients of intermediate results can reach hundreds of digits. To control this growth, the literature proposes several **pseudo-remainder sequences**:

- **Primitive PRS**: at each step, take the **primitive part** of the pseudo-remainder (divide by the content of the coefficients), keeping the intermediate remainders primitive and avoiding needless coefficient growth
- **Subresultant PRS**: carefully chosen scaling factors keep the intermediate coefficients at exactly the theoretical size of the subresultant determinants, avoiding the excessive growth of the primitive PRS

In oCAS, `DenseUnivariatePolynomial::gcd` uses the Euclidean algorithm with pseudo-remainders (`ocas-poly/src/gcd.rs`), returning a primitive GCD. Since pseudo-remainders explode at higher degrees (in practice roughly $\deg > 16$), `gcd_modular_z` (Brown's modular GCD, `ocas-poly/src/gcd/modular.rs`) should be used instead there; the two are chosen explicitly by the caller — there is no automatic degree-based dispatch.

### Resultant

The **resultant** $\text{Res}(f, g)$ of two polynomials $f, g \in \mathbb{F}[x]$ is a scalar satisfying:

$$\text{Res}(f, g) = 0 \iff \gcd(f, g) \ne 1$$

The resultant can be computed as the determinant of the Sylvester matrix, or more efficiently via a PRS (Brown's PRS algorithm). In oCAS, the `resultant()` method is implemented with the Brown PRS algorithm.

## Basic Concepts

### Definition of GCD

Let $D$ be a unique factorization domain (UFD) and $f, g \in D[x]$. The **greatest common divisor** $\gcd(f, g)$ of $f$ and $g$ is the monic polynomial $d$ satisfying:

1. $d \mid f$ and $d \mid g$ (common divisor)
2. for any common divisor $c$, we have $c \mid d$ (maximality)

In $\mathbb{F}[x]$ ($\mathbb{F}$ a field), the GCD is always monic (multiply by the inverse of the leading coefficient). In $\mathbb{Z}[x]$, the GCD is usually taken to be **primitive**, i.e. with coefficients whose GCD is 1.

### Content and Primitive Part

The **content** of a polynomial $f = \sum a_i x^i \in \mathbb{Z}[x]$ is defined as the GCD of its coefficients:

$$\text{cont}(f) = \gcd(a_0, a_1, \ldots, a_n)$$

The **primitive part** is defined as:

$$\text{pp}(f) = f / \text{cont}(f)$$

**Gauss's lemma**: the product of two primitive polynomials is again primitive. It follows that:

$$\gcd(f, g) = \gcd(\text{cont}(f),\, \text{cont}(g)) \cdot \gcd(\text{pp}(f),\, \text{pp}(g))$$

This allows the content (an integer GCD) and the primitive part (a polynomial GCD) to be computed separately.

In oCAS:

```rust
use ocas_poly::DenseUnivariatePolynomial;
use ocas_domain::IntegerDomain;

let f = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    vec![Integer::from(6), Integer::from(4), Integer::from(2)],
);
assert_eq!(f.content(), Integer::from(2));       // cont(f) = gcd(6,4,2) = 2
assert_eq!(f.primitive_part().coeffs(),           // pp(f) = 3 + 2x + x²
    &[Integer::from(3), Integer::from(2), Integer::from(1)]);
```

## Core Theory

### Square-free Factorization

#### Motivation

The first step of complete factorization is to decompose the polynomial into powers of square-free factors:

$$f = c \cdot g_1^1 \cdot g_2^2 \cdot g_3^3 \cdots$$

where the $g_i$ are pairwise coprime and square-free ($\gcd(g_i, g_i') = 1$). This greatly simplifies the subsequent irreducible factorization — each $g_i$ only needs to be factored once.

#### Yun's Algorithm

For fields of characteristic 0 (such as $\mathbb{Q}$ or $\mathbb{Z}$), Yun's algorithm is based on the observation that if $f = g_1 g_2^2 g_3^3 \cdots$, then

$$\gcd(f, f') = g_2 g_3^2 g_4^3 \cdots$$

**Algorithm steps**:

```
Input: f ∈ ℤ[x] (already made primitive)
Output: [(g₁, 1), (g₂, 2), …] such that f = ∏ gₖᵏ

1. f ← pp(f),  f' ← derivative(f)
2. g ← gcd(f, f')
3. w ← f / g
4. k ← 1
5. while w ≠ 1:
6.     h ← gcd(w, g)
7.     z ← w / h          // z = gₖ (square-free factor of multiplicity k)
8.     if z ≠ 1: output (z, k)
9.     g ← g / h
10.    w ← h
11.    k ← k + 1
```

**Key property**: the GCD computation at each iteration exactly removes the factors of higher multiplicity, so $z = w/h$ contains precisely the factors of multiplicity $k$.

**Complexity**: $O(n^2)$ coefficient operations ($n = \deg f$), with the GCD computations as the bottleneck.

#### The Special Case over Finite Fields

In $\mathbb{F}_p[x]$, when $p \mid \deg f$, the formal derivative $f'$ may vanish (for example $f(x) = x^p - 1$ has $f' = 0$ over $\mathbb{F}_p$). In that case $f$ is a $p$-th power of some polynomial: $f(x) = g(x^p) = g(x)^p$. The Musser/Bernardin algorithm detects this situation, takes the $p$-th root and recurses, multiplying the multiplicities by $p$.

In oCAS, `DenseUnivariatePolynomial::square_free_factorization()` implements Yun's algorithm (`ocas-poly/src/factor/mod.rs`), and `square_free_factorization_ff()` implements the finite-field version (`ocas-poly/src/factor/finite_field.rs`).

### Modular GCD (Brown 1971)

#### The Problem

Computing GCDs directly over $\mathbb{Z}[x]$ with pseudo-division leads to exponential coefficient growth. Brown's modular algorithm "projects" the problem onto several $\mathbb{F}_p[x]$, exploits the efficient GCD over fields, and reconstructs the answer with CRT.

#### Core Idea

Let $a, b \in \mathbb{Z}[x]$ and $g = \gcd(a, b)$ (primitive). The key observation:

$$\text{lc}(g) \mid \gcd(\text{lc}(a), \text{lc}(b)) \equiv \gamma$$

Therefore $\gamma^{-1} \cdot \text{lc}(g) \mid 1$, i.e. the monic GCD modulo $p$ multiplied by $\gamma$ gives a scaled version of the GCD over $\mathbb{Z}[x]$.

#### Algorithm Steps

```
Input: a, b ∈ ℤ[x] (nonzero)
Output: gcd(a, b) ∈ ℤ[x] (primitive)

1. a_p ← pp(a),  b_p ← pp(b)
2. γ ← gcd(lc(a_p), lc(b_p))    // a multiple of the leading coefficient of the GCD
3. images ← [],  best_deg ← ∞
4. for p in primes_from(2³⁰):    // start from primes greater than 2³⁰
5.     if p | γ: continue
6.     field ← FiniteField(p)
7.     g_p ← monic(gcd(a_p mod p, b_p mod p))
8.     g_scaled ← g_p · (γ mod p) · (lc(g_p)⁻¹ mod p)
9.     
10.    if deg(g_scaled) > best_deg: continue    // unlucky prime, discard
11.    if deg(g_scaled) < best_deg:             // a smaller GCD found
12.        best_deg ← deg(g_scaled)
13.        images ← [(p, g_scaled)]
14.    else: images.append((p, g_scaled))
15.    
16.    // CRT reconstruction + trial division verification
17.    candidate ← primitive_part(CRT_reconstruct(images))
18.    if candidate divides a_p AND candidate divides b_p:
19.        return candidate
20.    
21.    if the number of primes tried > MAX_PRIMES: break
22.
23. return fallback_pseudo_remainder_gcd(a, b)  // extreme fallback
```

#### Key Details

**Prime selection**: primes start above $2^{30}$. Each prime contributes roughly 30 bits of precision to the CRT reconstruction, so a polynomial whose coefficients have bit-length $L$ typically needs only about $L/30$ primes (in practice at most a few dozen).

**Landau-$\gamma$ scaling**: the monic modular GCD $g_p$ must be multiplied by $\gamma = \gcd(\text{lc}(a), \text{lc}(b))$ to align with the GCD over $\mathbb{Z}[x]$. Specifically, $\text{lc}(g_p)^{-1} \cdot \gamma \bmod p$ is the scaling factor.

**CRT symmetric representatives**: reconstruction uses the symmetric residue range $(-p/2, p/2]$ rather than $[0, p)$, keeping the absolute values of coefficients minimal.

**Unlucky prime detection**: if the degree of the modular GCD for some prime $p$ is higher than the current best, that prime is "unlucky" ($p$ divides the denominator of some coefficient of the true GCD) and is simply discarded. If the degree is lower, all previous primes were unlucky, so the images are cleared and computation restarts.

**Trial division verification**: the CRT-reconstructed candidate $c$ must satisfy $c \mid a_p$ and $c \mid b_p$ (exact division); otherwise more primes are collected.

**Safe fallback**: after at most 10000 primes (in practice rarely more than a few dozen are needed), the algorithm falls back to a pseudo-remainder GCD to guarantee termination.

In oCAS, `gcd_modular_z()` is implemented in `ocas-poly/src/gcd/modular.rs`:

```rust
use ocas_poly::gcd::modular::gcd_modular_z;

let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]); // x² - 1
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);  // x² + 2x + 1
let g = gcd_modular_z(&a, &b);
assert_eq!(g.coeffs(), &[i(1), i(1)]);  // x + 1
```

### Multivariate GCD (Brown Evaluation–Interpolation)

#### The Problem

Given $f, g \in \mathbb{Z}[x_1, \ldots, x_n]$, compute $\gcd(f, g)$. The challenges for multivariate GCD are: (1) multivariate polynomials have no canonical "division with remainder"; (2) the coefficient structure is more complex.

#### Brown's Evaluation–Interpolation Strategy (Bivariate Case)

Let $f, g \in \mathbb{Z}[x, y]$, viewed as polynomials in the main variable $x$ with coefficients in $\mathbb{Z}[y]$.

**Steps**:

1. **Evaluate**: choose $\alpha \in \mathbb{Z}$ and compute $f(x, \alpha)$ and $g(x, \alpha)$ (univariate polynomials)
2. **Univariate GCD**: $d_\alpha(x) = \gcd(f(x, \alpha), g(x, \alpha))$
3. **Verify**: check that the degree of $d_\alpha(x)$ equals the degree of the true GCD (take the minimal consistent degree over several $\alpha$)
4. **Lagrange interpolation**: for sufficiently many $\alpha_1, \ldots, \alpha_k$, compute the $d_{\alpha_i}(x)$ and treat each $x$-coefficient as a polynomial in $y$, recovered by Lagrange interpolation to obtain $d(x, y)$

**Key condition**: $\alpha$ must be "lucky" — i.e. the degree of $\gcd(f(x, \alpha), g(x, \alpha))$ equals the degree of the true GCD. If $\alpha$ is unlucky (for example it makes a leading coefficient vanish), it must be retried.

#### Recursion for the Multivariate Case

An $n$-variable GCD reduces recursively by one variable:

1. Choose an evaluation point $\alpha$ for $x_n$; compute $f|_{x_n = \alpha}$ and $g|_{x_n = \alpha}$
2. Recursively compute the GCD in $n-1$ variables
3. Interpolate to recover the dependence on $x_n$

#### Embedding $\mathbb{Z}$ into $\mathbb{Q}$

When computing GCDs over $\mathbb{Z}[x_1, \ldots, x_n]$, handling the content is crucial:

1. First extract the **content with respect to the main variable** $x_0$ (content in $x_0$), i.e. the GCD of all $x_0^k$ coefficients
2. Compute the multivariate GCD of the primitive parts
3. The final GCD is $\gcd(\text{cont}_x(f), \text{cont}_x(g)) \cdot \gcd(\text{pp}_x(f), \text{pp}_x(g))$

### Hensel Lifting

#### Core Idea

Hensel lifting is the key technique for "lifting" a factorization modulo $p$ to one modulo $p^k$ (or to $\mathbb{Z}$). Let $f \in \mathbb{Z}[x]$ be monic with

$$f \equiv g_0 \cdot h_0 \pmod{p}, \quad \gcd(g_0, h_0) = 1 \pmod{p}$$

The goal is to find $g, h \in \mathbb{Z}[x]$ such that $f = g \cdot h$, with $g \equiv g_0 \pmod{p}$ and $h \equiv h_0 \pmod{p}$.

#### Linear Hensel Lifting (Two-Factor Case)

Let $s, t$ be the Bézout coefficients over $\mathbb{F}_p[x]$: $s \cdot g_0 + t \cdot h_0 = 1 \pmod{p}$.

**Iteration steps** (from $m = p$ until $m > B$):

```
1. e ← f - g·h                    // the error
2. if e = 0: return (g, h)        // exact factorization
3. ē ← (e / m) mod p              // the "derivative" of the error
4. Δg ← (t · ē) mod g₀           // correction to g (deg < deg g₀)
5. Δh ← (ē - Δg · h₀) / g₀      // correction to h (exact division)
6. g ← g + Δg · m
7. h ← h + Δh · m
8. m ← m · p
```

Each iteration lifts the precision from $\bmod m$ to $\bmod m \cdot p$.

#### The Mignotte Bound

When should lifting stop? We need an upper bound on the size of the factor coefficients. The **Landau–Mignotte bound** states:

$$\|g\|_\infty \le 2^n \|f\|_2$$

where $n = \deg f$ and $\|f\|_2 = \sqrt{\sum a_i^2}$ is the 2-norm. It therefore suffices to lift until $p^k > 2 \cdot 2^n \|f\|_2$.

#### Multi-Factor Hensel Lifting

For $r > 2$ factors, a **factor-by-factor peeling** strategy is used:

```
Input: f, [g₁, g₂, …, gᵣ] mod p
Output: [G₁, G₂, …, Gᵣ] such that f = ∏ Gᵢ

1. for i = 1, …, r-1:
2.     h₀ ← g_{i+1} · g_{i+2} · … · gᵣ  (the product of the remaining factors)
3.     (Gᵢ, H) ← hensel_lift_pair(f_current, gᵢ, h₀, p, bound)
4.     f_current ← H
5. Gᵣ ← f_current
```

Each step is a two-factor Hensel lift.

#### Zassenhaus Factor Recombination

Hensel lifting produces $r$ factors $\tilde{g}_1, \ldots, \tilde{g}_r$ modulo $p^k$. The true factors over $\mathbb{Z}[x]$ are products of some subset of the $\tilde{g}_i$ (reduced modulo $p^k$ into the symmetric range $(-p^k/2, p^k/2]$).

**Recombination strategy**:

```
1. remaining ← [g̃₁, …, g̃ᵣ]
2. result ← []
3. for size = 1, 2, …, r:
4.     for each subset S ⊂ remaining of size |S|:
5.         candidate ← primitive_part(lc(rest) · ∏_{i∈S} g̃ᵢ)
6.         reduce candidate to symmetric range mod p^k
7.         if candidate divides f:
8.             result.append(candidate)
9.             remaining ← remaining \ S
10.            rest ← f / (∏ result)
11.            break  // restart from size=1
12. return result
```

**Key optimization**: multiplying by the leading coefficient of the current cofactor before taking the primitive part is Zassenhaus's classical trick — a true factor $h$ satisfies $\text{lc}(\text{rest}) \cdot \prod S = c \cdot h$.

**Complexity**: worst case $2^r$ subsets ($r$ = the number of mod-$p$ factors), but in practice the number of factors is usually small.

In oCAS, Hensel lifting and Zassenhaus recombination are implemented in `ocas-poly/src/factor/hensel.rs`:

```rust
// The complete ℤ[x] factorization pipeline:
// 1. Square-free factorization
// 2. For each square-free component:
//    a. Choose a prime p (not dividing the leading coefficient, and f mod p square-free)
//    b. Factor over 𝔽_p (Cantor-Zassenhaus)
//    c. Hensel lift to ℤ
//    d. Zassenhaus recombination
// Entry point: DenseUnivariatePolynomial::factor()
```

### The Berlekamp Algorithm

#### Applicability

Factoring square-free polynomials over the finite field $\mathbb{F}_p$ with small prime $p$. When $p$ is small ($p \le 1000$), the matrix method is more efficient than the randomized Cantor–Zassenhaus algorithm.

#### Theoretical Foundation

Let $f \in \mathbb{F}_p[x]$ be a square-free polynomial of degree $n$. Define the **Frobenius matrix** $Q \in \mathbb{F}_p^{n \times n}$:

$$Q_{ij} = [x^j] \text{ coefficient in } x^{ip} \bmod f$$

i.e. row $i$ of $Q$ is the coefficient vector of $x^{ip} \bmod f$.

**Key theorem**: $v \in \mathbb{F}_p^n$ satisfies $Q^T v = v$ (i.e. $v$ lies in the null space of $Q^T - I$) if and only if the corresponding polynomial $v(x) = \sum v_i x^i$ satisfies

$$v(x)^p \equiv v(x) \pmod{f}$$

This means $v(x)$ takes values in $\mathbb{F}_p$, and

$$f \mid \prod_{a \in \mathbb{F}_p} (v(x) - a)$$

Therefore $\gcd(f, v(x) - a)$ gives a nontrivial factor of $f$ for some $a \in \mathbb{F}_p$.

#### Algorithm Steps

```
Input: monic square-free f ∈ 𝔽_p[x], deg f = n
Output: the list of irreducible factors of f

1. Construct the Frobenius matrix Q
2. Compute a basis {v₁, v₂, …, vᵣ} of the null space of Q^T - I
3. factors ← [f]
4. for each vⱼ (r ≥ 1):
5.     new_factors ← []
6.     for each factor g in factors:
7.         if deg(g) ≤ 1:
8.             new_factors.append(g); continue
9.         for a = 0, 1, …, p-1:
10.            d ← gcd(g, vⱼ - a)
11.            if 0 < deg(d) < deg(g):
12.                new_factors.append(d)
13.                g ← g / d
14.        new_factors.append(g)
15.    factors ← new_factors
16. return factors
```

**Null space dimension**: if $f$ has $r$ irreducible factors, the null space has dimension $r$. Hence only $r - 1$ nontrivial null space vectors are needed for a complete split.

**Complexity**: $O(n^3)$ (matrix elimination) + $O(r \cdot n^2 \cdot p)$ (GCD splitting). When $p$ is large, the cost of the $a$-loop becomes prohibitive and Cantor–Zassenhaus should be used instead.

In oCAS, `berlekamp()` is implemented in `ocas-poly/src/factor/finite_field.rs`. It is used automatically when $p \le 1000$.

### The Cantor–Zassenhaus Algorithm

#### Applicability

Factoring square-free polynomials over the finite field $\mathbb{F}_p$ with large prime $p$. It proceeds in two stages: **distinct-degree factorization** (DDF) and **equal-degree factorization** (EDF).

#### Distinct-Degree Factorization (DDF)

**Goal**: decompose $f$ as $f = g_1 \cdot g_2 \cdots g_s$, where every irreducible factor of $g_d$ has degree exactly $d$.

**Theoretical foundation**: the elements of $\mathbb{F}_{p^d}$ are exactly the roots of $x^{p^d} - x = 0$. Hence the irreducible factors of $f$ of degree $d$ divide exactly

$$\gcd(f,\, x^{p^d} - x)$$

but do not divide the earlier $x^{p^{d'}} - x$ ($d' < d$).

**Algorithm**:

```
Input: monic square-free f ∈ 𝔽_p[x]
Output: [(g₁, 1), (g₂, 2), …] where the factors of gₖ have degree k

1. current ← f,  h ← x,  degree ← 1
2. while deg(current) ≥ 2·degree:
3.     h ← h^p mod current           // Frobenius iteration: h = x^(p^degree)
4.     g ← gcd(current, h - x)
5.     if deg(g) > 0:
6.          output (monic(g), degree)
7.         current ← current / g
8.         h ← h mod current
9.     degree ← degree + 1
10. if deg(current) > 0:
11.      output (monic(current), deg(current))
```

**Optimizing the Frobenius iteration**: $x^{p^d}$ need not be computed directly (the exponent is enormous); instead iterate $h \leftarrow h^p \bmod f$ step by step, using fast modular exponentiation each time.

**Complexity**: $O(n^2 \log p)$ per Frobenius iteration, $O(n^3 \log p)$ in total.

#### Equal-Degree Factorization (EDF)

**Goal**: completely split $f = g_1 g_2 \cdots g_r$ (each $g_i$ irreducible of degree $d$).

**Odd characteristic case** ($p > 2$):

Use the fact that exactly half of the elements of $\mathbb{F}_{p^d}$ are quadratic residues. For a random polynomial $a$, compute

$$b = a^{(p^d - 1)/2} \bmod f$$

Then $b$ takes the value $1$ on each irreducible factor $g_i$ where $a$ is a quadratic residue in that factor's field, and $-1$ where it is a non-residue. Therefore

$$\gcd(f,\, b - 1) = \prod_{i:\, a \text{ is a quadratic residue in } g_i} g_i$$

gives a nontrivial split of $f$ (unless all $g_i$ agree, which happens with probability $2^{1-r}$).

**Characteristic 2 case** ($p = 2$):

The $b - 1$ trick does not apply ($1 = -1$). Instead use the **trace map**:

$$T(a) = a + a^2 + a^{2^2} + \cdots + a^{2^{d-1}} \bmod f$$

$T(a)$ takes values in $\mathbb{F}_2$ on each $\mathbb{F}_{2^d}$ (a property of the trace), and for random $a$ the values on different factors are independent. Hence $\gcd(f, T(a))$ gives a nontrivial split.

```
Input: f (the product of the d-degree factors from DDF)
Output: the list of irreducible factors of f

1. factors ← [f]
2. while there exists a factor with deg > d:
3.     for each factor g with deg(g) > d:
4.         choose a random polynomial a (deg < deg(g))
5.         if p = 2:
6.             b ← T(a) = Σᵢ₌₀^{d-1} a^{2^i} mod g
7.         else:
8.             b ← a^{(p^d-1)/2} mod g
9.         d₁ ← gcd(g, b - 1)   // or gcd(g, b) for char 2
10.        if 0 < deg(d₁) < deg(g):
11.            replace g by d₁ and g/d₁ in factors
12. return factors
```

In oCAS, `cantor_zassenhaus()` is implemented in `ocas-poly/src/factor/finite_field.rs`; the top-level entry `factor_over_finite_field()` automatically chooses Berlekamp ($p \le 1000$) or Cantor–Zassenhaus according to the size of $p$.

### Wang's EEZ Multivariate Factorization

#### The Problem

Given $f \in \mathbb{Z}[x_1, \ldots, x_n]$ (or $\mathbb{F}_p[x_1, \ldots, x_n]$), factor it into irreducible factors.

#### Strategy Overview

Wang's (1978) EEZ (Evaluation, Exact division, Zassenhaus) algorithm generalizes univariate Hensel lifting to the multivariate case:

1. **Evaluate**: choose evaluation points for the auxiliary variables to obtain a univariate image $f(x_1, \alpha_2, \ldots, \alpha_n)$
2. **Univariate factorization**: factor the univariate image over $\mathbb{Z}$ (or $\mathbb{F}_p$)
3. **Variable-by-variable Hensel lifting**: lift back to the multivariate factorization one variable at a time
4. **Zassenhaus recombination**: combine the lifted factors

#### Wang's Leading-Coefficient Preprocessing

One difficulty of multivariate factorization is the nonconstant leading coefficient. Let $\ell(x_1, \ldots, x_n)$ be the leading coefficient of $f$ in the main variable $x_0$. If $\ell$ is not constant, the leading coefficients $\ell_i$ of the factors $f_i$ need not be constant either, and $\prod \ell_i = \ell$.

**Wang's greedy assignment**:

1. Factor $\ell$ into powers of irreducible factors: $\ell = \prod g_j^{e_j}$
2. At the evaluation point $\alpha = (\alpha_2, \ldots, \alpha_n)$, the image of $\ell$ splits into the images of the univariate factors
3. Assign the nontrivial factors $g_j$ of $\ell$ to the corresponding $f_i$ so that $g_j(\alpha) = \text{lc}(u_i)$ ($u_i$ the univariate image factor)
4. Ensure consistency of the assignment via pairwise coprimality conditions ($\alpha_j = |g_j(\alpha)| > 1$)

#### Variable-by-Variable Hensel Lifting

Let $f$ have $n$ variables and let the univariate image factors be $u_1, \ldots, u_r$.

**Lifting the variable $x_k$** ($k = 1, 2, \ldots, n-1$, one variable at a time in index order):

1. Evaluate the current factors $F_i^{(k-1)}$ (in the variables $x_1, \ldots, x_{k-1}$) at $x_k = \alpha_k$ to obtain the $u_i$
2. Compute Bézout coefficients $b_i$ such that $\sum b_i \prod_{j \ne i} u_j = 1$
3. For $t = 1, 2, \ldots$, successively solve the **multivariate Diophantine equation**:

$$\sum_i \sigma_i \cdot \prod_{j \ne i} F_j = e_t$$

where $e_t$ is the $t$-th term of the Taylor expansion of the current error $f - \prod F_i$

4. Correct $F_i \leftarrow F_i + \sigma_i \cdot (x_k - \alpha_k)^t$

#### The Multivariate Diophantine Equation Solver

When lifting the variable $x_k$, the equations to solve have the form:

$$\sum_{i=1}^{r} \sigma_i \cdot \prod_{j \ne i} u_j = e \pmod{(x_k - \alpha_k)^t}$$

with $\deg_{x_0}(\sigma_i) < \deg_{x_0}(u_i)$. This is a linear system (in the coefficients of the $\sigma_i$) that can be solved recursively:

- when $k = 1$ (univariate): the extended Euclidean algorithm
- when $k > 1$: recurse on $k-1$ variables — evaluate at $\alpha_k$, solve the smaller system, then interpolate

#### p-adic Coefficient Hensel Lifting

For the nonconstant leading coefficient case over $\mathbb{Z}[x_1, \ldots, x_n]$, Wang EEZ must be followed by **p-adic coefficient lifting**:

1. Reduce the coefficients of the multivariate factors modulo $p$ to obtain a skeleton over $\mathbb{F}_p$
2. Iteratively solve the mod-$p$ Diophantine equations, lifting the coefficients step by step from $\bmod p$ to $\bmod p^k$
3. Continue until the error vanishes or $p^k$ exceeds the coefficient bound (the Gelfond bound)

**The Gelfond coefficient bound**:

$$B = \left(\sqrt{\prod_v (d_v + 1) \cdot 2^{2 \sum d_v - n}} + 1\right) \cdot \|f\|_\infty \cdot |\text{lc}(f)|$$

where $d_v$ is the degree in the variable $v$ and $n$ is the number of variables.

#### Sparse Diophantine Solving

When the Diophantine correction terms have a sparse structure, dense recursive solving can be avoided. oCAS implements a **skeleton interpolation** strategy:

1. Extract the skeleton from the terms of the error $e$ — the possible exponent patterns of the correction terms
2. Evaluate at several random base points to obtain univariate Diophantine equations
3. Interpolate the coefficients through a Vandermonde system
4. Verify the interpolated result

This is especially effective when the factors have many variables but low degree in each.

In oCAS, the general EEZ algorithm is implemented in `ocas-poly/src/factor/eez.rs`, with a bivariate specialization in `ocas-poly/src/factor/multivariate.rs` (which currently requires the leading coefficient in the main variable to be constant):

```rust
// Multivariate factorization entry over ℤ (eez.rs):
pub fn multivariate_factor_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)>

// Multivariate factorization entry over 𝔽_p (eez.rs):
pub fn multivariate_factor_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)>

// Bivariate specialization entries (multivariate.rs):
pub fn bivariate_factor_z(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<(ZMPoly, usize)>
pub fn bivariate_factor_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<(FpMPoly, usize)>
```

### Trager Factorization over Algebraic Number Fields

#### The Problem

Given $f \in K[x]$ where $K = \mathbb{Q}(\alpha)$ is an algebraic number field ($\alpha$ a root of the minimal polynomial $m(\alpha) = 0$), factor $f$ into irreducible factors in $K[x]$.

#### Reduction via the Norm

The core idea of Trager's algorithm is to reduce the problem from $K[x]$ to $\mathbb{Q}[x]$ via the **norm**:

$$N(f) = \text{Res}_\alpha(m(\alpha),\, f(x, \alpha))$$

The norm $N(f) \in \mathbb{Q}[x]$ has degree $\deg_x(f) \cdot [K:\mathbb{Q}]$, and if $f$ is reducible in $K[x]$, then $N(f)$ is reducible in $\mathbb{Q}[x]$.

#### Computing the Norm by Evaluation–Interpolation

Constructing the Sylvester matrix directly to compute the resultant is expensive. oCAS uses **evaluation–interpolation**:

1. Compute the scalar resultants $\text{Res}_\alpha(m, f(x_j, \alpha))$ at $\deg_x(f) \cdot [K:\mathbb{Q}] + 1$ rational points $x_j$
2. Each scalar resultant is the norm of an element of $\mathbb{Q}(\alpha)$, obtained by evaluating $m$ at $[K:\mathbb{Q}]$ points and taking the product
3. Recover $N(f) \in \mathbb{Q}[x]$ by Newton divided-difference interpolation

#### The Trager Shift

If $N(f)$ has repeated factors, information about the factors of $f$ is lost. The **Trager shift** substitutes $x \mapsto x - s\alpha$ ($s \ge 0$) to make the norm square-free:

$$N(f(x - s\alpha)) \text{ is square-free}$$

Such an $s$ always exists (only finitely many $s$ are bad); in practice $s = 0$ or a small value suffices.

#### Factorization Steps

```
Input: f ∈ K[x] (square-free)
Output: the irreducible factors of f in K[x]

1. for s = 0, 1, 2, …:
2.     g ← f(x - s·α)
3.     N ← Res_α(m, g(x, α))  via evaluation-interpolation
4.     if N is square-free: break
5. 
6. N₁, …, Nₖ ← factor_over_Q(N)   // factor over ℚ
7. for each Nᵢ:
8.     N̂ᵢ ← embed_Q_to_K(Nᵢ)
9.     hᵢ ← gcd_K(g, N̂ᵢ)          // GCD over K
10.    output ← compose_linear(hᵢ, +s·α)   // inverse shift
11. return [monic(output₁), …]
```

#### Modular GCD over Algebraic Number Fields

The $K[x]$ GCD in step 9 is computed efficiently by modular methods (`gcd_anf`):

1. Choose a prime $p$ such that $m$ is irreducible over $\mathbb{F}_p$ (ensuring $\mathbb{F}_p[\alpha]/(m) \cong \text{GF}(p^d)$)
2. Map $a, b \in K[x]$ to $\text{GF}(p^d)[x]$
3. Compute the monic GCD over $\text{GF}(p^d)[x]$
4. CRT combination + rational reconstruction + trial division verification

Unlucky primes (where the modular GCD degree is too large) are discarded; after at most 1000 primes the algorithm falls back to a dense Euclidean GCD.

In oCAS, Trager's algorithm is implemented in `ocas-poly/src/factor/algebraic.rs`:

```rust
// Factorization entry over K = ℚ(α):
impl DenseUnivariatePolynomial<AlgebraicNumberField> {
    pub fn factor(&self) -> Factors<AlgebraicNumberField>
}

// Example: x² - 2 factors over ℚ(√2) as (x - √2)(x + √2)
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
let f = DenseUnivariatePolynomial::from_coeffs(field.clone(), vec![
    field.from_base(Rational::new(-2, 1)),
    field.zero(),
    field.one(),
]);
let factors = f.factor();
assert_eq!(factors.len(), 2);  // (x - √2) and (x + √2)
```

## Implementation in oCAS

### File Map

| Source file | Algorithm | Description |
|---|---|---|
| `ocas-poly/src/dense.rs` | Euclidean GCD, Karatsuba multiplication | Core data structure for univariate polynomials |
| `ocas-poly/src/gcd/modular.rs` | Brown modular GCD | Efficient GCD over $\mathbb{Z}[x]$ |
| `ocas-poly/src/factor/mod.rs` | Top-level entry, Yun square-free factorization | `factor()` and `square_free_factorization()` |
| `ocas-poly/src/factor/finite_field.rs` | Berlekamp, Cantor–Zassenhaus (DDF + EDF) | Factorization over $\mathbb{F}_p[x]$ |
| `ocas-poly/src/factor/hensel.rs` | Hensel lifting, Zassenhaus recombination, Mignotte bound | The $\mathbb{Z}[x]$ factorization pipeline |
| `ocas-poly/src/factor/multivariate.rs` | Bivariate factorization (Wang EEZ specialization; leading coefficient in the main variable must be constant) | $\mathbb{Z}[x,y]$ and $\mathbb{F}_p[x,y]$ |
| `ocas-poly/src/factor/eez.rs` | General multivariate EEZ Hensel lifting | $\mathbb{Z}[x_1,\ldots,x_n]$ and $\mathbb{F}_p[x_1,\ldots,x_n]$ |
| `ocas-poly/src/factor/algebraic.rs` | Trager's algorithm (factorization over algebraic number fields) | $\mathbb{Q}(\alpha)[x]$ |
| `ocas-poly/src/multivariate_gcd.rs` | Multivariate GCD (evaluation–interpolation) | Supports multivariate factorization |

### Algorithm Selection Strategy

oCAS automatically selects the best algorithm according to the input:

```
factor(f ∈ ℤ[x]):
  1. Square-free factorization (Yun's algorithm)
  2. For each square-free component g:
     a. if deg(g) ≤ 1: return it directly
     b. choose a prime p (not dividing lc(g), and g mod p square-free)
     c. factor over 𝔽_p (Cantor–Zassenhaus; the general factor_over_finite_field entry switches to Berlekamp automatically when p ≤ 1000)
     d. if the number of mod-p factors = 1: g is irreducible, return
     e. Hensel lifting (the Mignotte bound determines the precision)
     f. Zassenhaus factor recombination
     g. non-monic case: leading coefficient transformation a^{d-1}·f(x/a)

factor(f ∈ 𝔽_p[x]):
  1. Extract the leading coefficient
  2. Square-free factorization (Musser/Bernardin, handling characteristic p)
  3. For each square-free component: Berlekamp (small p) or Cantor-Zassenhaus

factor(f ∈ ℚ(α)[x]):
  1. Yun square-free factorization (using the modular GCD gcd_anf)
  2. For each component: Trager norm + ℚ factorization + K-GCD recovery

factor(f ∈ ℤ[x₁,…,xₙ]):
  1. Extract the content (in the main variable x₀)
  2. Square-free factorization
  3. Wang leading-coefficient preprocessing
  4. Choose evaluation points, factor the univariate image
  5. EEZ variable-by-variable Hensel lifting
  6. p-adic coefficient lifting (when the leading coefficient is nonconstant)
  7. Zassenhaus recombination
```

### Performance Characteristics

| Algorithm | Time complexity | Space complexity | Use case |
|---|---|---|---|
| Euclidean GCD (PRS) | $O(n^2)$ | $O(n)$ | Low degree ($\lesssim 16$) |
| Brown modular GCD | $O(n^2 \cdot k)$ | $O(n \cdot k)$ | High-degree $\mathbb{Z}[x]$, $k$ = number of primes |
| Hensel lifting | $O(n^2 \cdot \log B)$ | $O(n \cdot r)$ | $r$ factors, $B$ = Mignotte bound |
| Zassenhaus recombination | $O(2^r \cdot n)$ | $O(n \cdot r)$ | $r$ mod-$p$ factors |
| Berlekamp | $O(n^3 + r \cdot n^2 p)$ | $O(n^2)$ | Small $p$ |
| Cantor–Zassenhaus | $O(n^3 \log p)$ | $O(n)$ | Large $p$ |
| Wang EEZ | Exponential in the number of variables | Multivariate polynomials | Multivariate factorization |
| Trager | $O(n^3 [K:\mathbb{Q}]^2)$ | $O(n \cdot [K:\mathbb{Q}])$ | Algebraic number fields |

## References

- **[Gathen–Gerhard]** J. von zur Gathen and J. Gerhard, *Modern Computer Algebra*, 3rd ed., Cambridge University Press, 2013. Chapters 14–15 (factoring), 18 (GCD).
- **[Brown]** W. S. Brown, "On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors," *J. ACM*, 18(4):478–504, 1971.
- **[Zassenhaus]** H. Zassenhaus, "On Hensel Factorization I," *J. Number Theory*, 1(3):291–311, 1969.
- **[Wang]** P. S. Wang, "An Improved Multivariate Polynomial Factoring Algorithm," *Math. Comp.*, 32(144):1215–1231, 1978.
- **[Berlekamp]** E. R. Berlekamp, "Factoring Polynomials over Large Finite Fields," *Math. Comp.*, 24(111):713–735, 1970.
- **[Cantor–Zassenhaus]** D. G. Cantor and H. Zassenhaus, "A New Algorithm for Factoring Polynomials over Finite Fields," *Math. Comp.*, 36(154):587–592, 1981.
- **[Trager]** B. M. Trager, "Algebraic Factoring and Rational Function Integration," *Proc. SYMSAC '76*, pp. 219–226, 1976.
- **[Mignotte]** M. Mignotte, *Mathematics for Computer Algebra*, Springer, 1992.
- **[Geddes–Czapor–Labahn]** K. O. Geddes, S. R. Czapor, and G. Labahn, *Algorithms for Computer Algebra*, Kluwer, 1992. Chapters 6 (multivariate factoring), 7 (GCD).
- **[Knuth]** D. E. Knuth, *The Art of Computer Programming*, Vol. 2, §4.6.2, Addison-Wesley, 1997.
