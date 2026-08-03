# Advanced: Number Theory Algorithms

This chapter systematically explains the core number theory algorithms implemented in oCAS: primality testing, integer factorization, discrete logarithms, the Chinese remainder theorem, quadratic residues, and rational reconstruction. Each algorithm starts from its mathematical principle and goes into the concrete implementation details in oCAS.

## Prerequisites

Before reading this chapter, the reader should be familiar with the following concepts:

- **Primes and divisibility**: definition of primes, greatest common divisors, the Euclidean algorithm, the extended Euclidean algorithm
- **Modular arithmetic**: congruences, modular inverse (computing $a^{-1} \bmod m$ with the extended Euclidean algorithm), fast modular exponentiation (binary exponentiation)
- **Euler's φ function**: $\varphi(n) = n \prod_{p \mid n}(1 - 1/p)$, counting the integers in $[1, n]$ coprime to $n$
- **Multiplicative functions**: $f(mn) = f(m)f(n)$ when $\gcd(m, n) = 1$; completely multiplicative when this holds for all $m, n$
- **Fermat's little theorem**: $a^{p-1} \equiv 1 \pmod{p}$ for prime $p$ with $\gcd(a, p) = 1$
- **Basics of group theory**: cyclic groups, the order of an element, generators (primitive roots)

**Recommended reading**: Shoup, *A Computational Introduction to Number Theory and Algebra*, Ch.1–9.

## Basic Concepts

### The Primality Testing Problem

Given a positive integer $n$, decide whether $n$ is prime. This is a decision problem, but deterministic algorithms (such as AKS) are far slower in practice than probabilistic ones. Modern practice uses **probable prime tests**: the probability that a composite passes the test can be made extremely small, but it cannot be absolutely ruled out.

### The Integer Factorization Problem

Given a composite $n$, find a nontrivial factor $d$ ($1 < d < n$). Repeating recursively until all factors are prime yields $n = p_1^{e_1} \cdots p_k^{e_k}$. The difficulty of factorization is the foundational assumption behind the security of public-key cryptography such as RSA.

### The Discrete Logarithm Problem

In the group $\mathbb{Z}_p^*$, given $g$ and $h$, find $x$ such that $g^x \equiv h \pmod{p}$. When $p - 1$ is smooth (has only small prime factors), the Pohlig–Hellman algorithm solves this efficiently.

### Quadratic Residues

For an odd prime $p$ and an integer $a$, if there exists $x$ with $x^2 \equiv a \pmod{p}$, then $a$ is called a **quadratic residue** modulo $p$. The Legendre symbol $\left(\frac{a}{p}\right)$ equals $1$ for a residue, $-1$ for a non-residue, and $0$ when $p \mid a$.

### Multiplicative Number-Theoretic Functions

| Function | Definition | Multiplicative |
|---|---|---|
| $\varphi(n)$ | $\#\{k \in [1,n] : \gcd(k,n)=1\}$ | yes |
| $\mu(n)$ | $1$ ($n=1$), $(-1)^k$ ($n$ has exactly $k$ distinct prime factors), $0$ ($n$ has a square factor) | yes |
| $\tau(n)$ | number of positive divisors $= \prod(e_i + 1)$ | yes |
| $\sigma_k(n)$ | $\sum_{d \mid n} d^k$ | yes |
| $\lambda(n)$ | Liouville function $= (-1)^{\Omega(n)}$ | yes (completely multiplicative) |

## Core Theory

### BPSW Primality Testing

The **BPSW test** (Baillie–Pomerance–Selfridge–Wagstaff) combines two independent probable prime tests: a base-2 strong Miller–Rabin test plus a strong Lucas probable prime test. Since it was proposed in 1980, **no composite number is known to pass the BPSW test** (it has been exhaustively verified for $n < 2^{64}$).

#### The Miller–Rabin Test

**Principle**: write $n - 1 = d \cdot 2^r$ ($d$ odd). For a base $a$, compute $a^d \bmod n$. If the result is $1$, or $-1$ appears in the sequence $a^d, a^{2d}, a^{4d}, \ldots, a^{d \cdot 2^{r-1}} \pmod{n}$, then $n$ passes this round.

**Guarantee**: if $n$ is prime, it passes for every $1 < a < n$. If $n$ is composite, at least $3/4$ of the bases are witnesses, i.e. they expose $n$. After testing $k$ random bases, the probability that a composite passes all of them is at most $4^{-k}$.

**Determinism**: for $n < 3.317 \times 10^{24}$, the Miller–Rabin test with the fixed base set $\{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37\}$ is **deterministic**. oCAS uses these 12 bases as the `MR_WITNESSES` constant.

#### The Strong Lucas Probable Prime Test

**Lucas sequences**: for parameters $(P, Q)$, define

$$U_0 = 0,\quad U_1 = 1,\quad U_{m+1} = P \cdot U_m - Q \cdot U_{m-1}$$

$$V_0 = 2,\quad V_1 = P,\quad V_{m+1} = P \cdot V_m - Q \cdot V_{m-1}$$

**Binary ladder**: compute $(U_k, V_k, Q^k) \bmod n$ efficiently. Use the doubling formulas

$$U_{2m} = U_m \cdot V_m, \quad V_{2m} = V_m^2 - 2Q^m$$

and the increment formulas

$$2U_{m+1} = P \cdot U_m + V_m, \quad 2V_{m+1} = D \cdot U_m + P \cdot V_m$$

where $D = P^2 - 4Q$. Divisions are avoided through the half-mod operation modulo $n$ (for odd $n$, $\text{half\_mod}(x, n) = y$ satisfies $2y \equiv x \pmod{n}$).

**Selfridge parameter selection**: choose the first $D \in \{5, -7, 9, -11, 13, \ldots\}$ with Jacobi symbol $\left(\frac{D}{n}\right) = -1$, then set $P = 1$ and $Q = (1 - D)/4$.

**The strong Lucas test**: write $n + 1 = d \cdot 2^r$ ($d$ odd); $n$ passes if and only if $U_d \equiv 0 \pmod{n}$ or there exists $0 \leq i < r$ with $V_{d \cdot 2^i} \equiv 0 \pmod{n}$.

#### The BPSW Combination

```
is_prime_bpsw(n):
    1. if n < 2, return false; if n equals a prime in {2,3,5,…,37}, return true
    2. if n is divisible by some prime in 2..37, return false
    3. if n is a perfect square, return false
    4. strong Miller–Rabin test (base 2)
    5. if it passes, run the strong Lucas PRP test (Selfridge parameters)
    6. return true only if both pass
```

In oCAS, `is_prime()` uses the full 12-base Miller–Rabin, while `is_prime_bpsw()` uses the single-base Miller–Rabin plus the Lucas test. The two complement each other: `is_prime()` is deterministic for $n < 3.317 \times 10^{24}$; `is_prime_bpsw()` is theoretically stronger (it combines two tests based on different principles).

### The Integer Factorization Strategy

oCAS's `factor_integer` is a **layered driver** that combines several factorization algorithms into an escalating strategy.

#### Trial Division

First generate the list of primes $\leq 1000$ with the sieve of Eratosthenes and divide by each in turn. This step costs $O(\pi(1000)) = O(168)$, negligible, but it quickly strips off small factors.

$$n = \prod_{p \leq 1000} p^{e_p} \cdot C$$

The remaining cofactor $C$ has all prime factors $> 1000$.

#### The Pollard rho–Brent Variant

**Principle**: random collision detection based on the birthday paradox. Define a sequence $x_0, x_1, x_2, \ldots$ with $x_{i+1} = x_i^2 + c \pmod{n}$ ($c$ a random constant). Modulo $n$ the sequence eventually enters a cycle, while modulo a true factor $p$ a collision occurs earlier.

The core optimizations of the **Brent variant**:

1. **Exponential search**: detect collisions at $2^k$-step intervals (Floyd's tortoise-and-hare checks at every step; Brent only checks at powers of two)
2. **Batched GCD**: accumulate the product of several $(x_i - y_i)$ values and take a single $\gcd$, amortizing the expensive GCD
3. **Backtracking**: when a batched GCD returns $n$ itself, retry the individual differences within that batch

**Time complexity**: expected $O(n^{1/4})$ modular multiplications. In practice it is very efficient for factors of moderate size (< 20 digits).

#### Pollard's $p - 1$ Method (Stage 1)

**Principle**: if $n$ has a factor $p$ and $p - 1$ is $B_1$-smooth (i.e. all prime power factors of $p - 1$ are $\leq B_1$), then $p - 1 \mid M$ where $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$. By Fermat's little theorem, for any $\gcd(a, p) = 1$:

$$a^M \equiv 1 \pmod{p} \implies p \mid \gcd(a^M - 1, n)$$

**Implementation details**:

- Generate the primes $q \leq B_1$ with the sieve of Eratosthenes
- Modular exponentiation per prime: $a \leftarrow a^{q^e} \bmod n$ ($q^e \leq B_1 < q^{e+1}$)
- Check $\gcd(a - 1, n)$ after accumulating a certain number of steps
- If $\gcd = n$ (all factors are $p-1$-smooth), the method fails
- If $1 < \gcd < n$, a nontrivial factor is found

**Limitation**: requires $p - 1$ to be smooth. If $p - 1$ has a large prime factor, the method does not work.

#### Williams's $p + 1$ Method (Stage 1)

**Principle**: complementary to the $p - 1$ method. If $p + 1$ is $B_1$-smooth, use the Lucas $V$ sequence instead of modular exponentiation.

**The Lucas $V$ sequence**: take $Q = 1$ and choose $P$ randomly with $\left(\frac{P^2 - 4}{n}\right) = -1$ (Jacobi symbol). Compute $V_M \bmod n$, where $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$. If $p + 1$ is $B_1$-smooth, then $V_M \equiv V_0 = 2 \pmod{p}$, i.e. $p \mid \gcd(V_M - 2, n)$.

**Parameter selection**: choose $P$ at random and check $\text{jacobi}(P^2 - 4, n) = -1$. The oCAS implementation uses `lucas_uv_mod` (the Lucas chain code shared with BPSW).

**Limitation**: requires $p + 1$ to be smooth. It complements the $p - 1$ method (some $p$ satisfy one but not the other).

#### Lenstra's Elliptic Curve Method (ECM)

**Principle**: reduce integer factorization to a smoothness question about the order of an elliptic curve group. For a random elliptic curve $E/\mathbb{F}_p$, the group order $|E(\mathbb{F}_p)|$ is distributed roughly uniformly in the Hasse interval $[p + 1 - 2\sqrt{p}, p + 1 + 2\sqrt{p}]$ (Sato–Tate). ECM succeeds when the group order of some curve is $B_1$-smooth.

**Key advantage**: unlike rho and the $p \pm 1$ methods, ECM does not fix the group order range — each new curve gives a new group order. This makes ECM succeed with positive probability for **factors of any size**, with complexity depending on the smoothness of the smallest factor $p$ rather than on $n$ itself.

**Suyama parameterization**: construct a Montgomery curve from a random $\sigma \notin \{0, 1, 5\}$:

$$u = \sigma^2 - 5, \quad v = 4\sigma$$

$$A = \frac{(v - u)^3(3u + v)}{4u^3 v} - 2$$

The curve is $B y^2 = x^3 + A x^2 + x$ with base point $(u^3, v^3)$. The Suyama parameterization guarantees that the base point is always rational, and if during the construction a GCD of a denominator exposes a factor of $n$, it is returned immediately.

**Montgomery coordinates**: use projective coordinates $(X : Z)$ (omitting the $Y$ coordinate), with curve equation $By^2 = x^3 + Ax^2 + x$ and $x = X/Z$.

- **Doubling** $[2](X : Z)$: computed with the constant $a_{24} = (A + 2)/4$, needing only 4 modular multiplications and 4 modular additions/subtractions
- **Differential addition** $P + Q$ given $P - Q$: only 6 modular multiplications
- **Montgomery ladder**: scalar multiplication $[k]P$ via alternating doublings and additions, **always using the base point as the fixed difference**

**Stage 1 algorithm**: compute $[M]P$ where $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$. Take GCDs during the computation; return as soon as a nontrivial factor is found.

#### The `factor_integer` Driver

`factor_integer` combines the above algorithms into an escalating strategy:

```
factor_integer(n):
    1. Trial division: strip off all factors p ≤ 1000
    2. For the cofactor C (if C > 1):
       a. if is_prime_bpsw(C), record it as a prime factor
       b. otherwise call find_factor(C) to split it, and recurse
```

The escalation strategy of `find_factor(n)`:

1. One Pollard rho–Brent attempt first (cheapest for small factors; no smoothness requirement)
2. Loop (starting with $B_1 = 2000$):
   - Pollard $p - 1$ (stage 1, requires some factor with $p - 1$ $B_1$-smooth)
   - Williams $p + 1$ (stage 1, requires some factor with $p + 1$ $B_1$-smooth; complements $p - 1$)
   - ECM: number of curves $\approx B_1 / 550$ (clamped to 10–300)
   - if all fail, $B_1 \leftarrow B_1 \times 4$ and continue the loop

The curve budget of about $B_1 / 550$ is the standard empirical ratio for ECM smoothness search, ensuring enough random curves per $B_1$ level to cover the group orders in the Hasse interval.

### Discrete Logarithms

#### Baby-Step Giant-Step (BSGS)

**Problem**: given $g, h, p$, find $x$ such that $g^x \equiv h \pmod{p}$ ($0 \leq x < m$, where $m$ is an upper bound on the order of $g$).

**Algorithm**: let $m = \lfloor \sqrt{\text{bound}} \rfloor + 1$.

1. **Baby steps**: compute and store $\{(j, g^j \bmod p) : 0 \leq j < m\}$ in a HashMap
2. **Giant steps**: compute $g^{-m} \bmod p$; for $i = 0, 1, \ldots, m-1$, check whether $h \cdot (g^{-m})^i \bmod p$ is in the table
3. If a match $g^j = h \cdot g^{-im}$ is found, then $x = im + j$

**Complexity**: $O(\sqrt{m})$ modular multiplications in time, $O(\sqrt{m})$ in space.

oCAS's `dlog_bsgs` stores the baby steps in a `HashMap<Integer, Integer>` with $(j, g^j)$ key-value pairs, giving $O(1)$ lookups.

#### The Pohlig–Hellman Algorithm

**Scenario**: $p$ prime, $g$ of order $n = p - 1$, with $n$ smooth: $n = q_1^{e_1} \cdots q_k^{e_k}$ ($q_i$ small primes).

**Core idea**: decompose the problem of finding $x \bmod n$ into the subproblems of finding $x \bmod q_i^{e_i}$, then combine with CRT.

**Prime-power decomposition**: for each $q = q_i$, $e = e_i$, find $x \bmod q^e$:

Write $x = x_0 + x_1 q + x_2 q^2 + \cdots + x_{e-1} q^{e-1}$ ($0 \leq x_j < q$; recover one "digit" at a time).

1. Let $g_0 = g^{n/q} \bmod p$, $h_0 = h^{n/q} \bmod p$. Then $g_0^{x_0} = h_0 \pmod{p}$. Use BSGS to find $x_0 \bmod q$.
2. Update $h \leftarrow h \cdot g^{-x_0}$, and compute $h_1 = (h \cdot g^{-x_0})^{n/q^2}$. Find $x_1$.
3. Continue until $x_{e-1}$.
4. $x \bmod q^e = x_0 + x_1 q + \cdots + x_{e-1} q^{e-1}$.

**CRT combination**: combine the results for all $q_i^{e_i}$ with `crt_many` to obtain $x \bmod n$.

**Complexity**: BSGS per prime power costs $O(\sqrt{q_i})$, for a total of $O(\sum e_i \sqrt{q_i})$. When all $q_i$ are small, this is much faster than the direct BSGS cost of $O(\sqrt{n})$.

oCAS's `dlog_pohlig_hellman` requires $p$ to be prime. It uses `factor_integer` to factor $p - 1$ (reusing the integer factorization module), recovers the digits prime power by prime power with BSGS, and finally combines and verifies the result with `crt_many`.

### The Chinese Remainder Theorem

#### The Classical CRT

**Theorem**: let $m_1, \ldots, m_k$ be pairwise coprime. Given the system of congruences

$$x \equiv r_i \pmod{m_i}, \quad i = 1, \ldots, k$$

there exists a unique solution $x \bmod M$ where $M = m_1 \cdots m_k$.

**Construction**: $x = \sum_i r_i \cdot M_i \cdot M_i^{-1}$, where $M_i = M / m_i$ and $M_i^{-1} = M_i^{-1} \bmod m_i$.

#### Generalization: Non-Coprime Moduli

oCAS's `crt` (pairwise combination) and `crt_many` (combining many congruences) **do not require the moduli to be pairwise coprime**.

**Pairwise CRT**: for $x \equiv r_1 \pmod{m_1}$ and $x \equiv r_2 \pmod{m_2}$:

1. Compute $g = \gcd(m_1, m_2)$
2. Check $r_1 \equiv r_2 \pmod{g}$ (the consistency condition)
3. If inconsistent, return `None`
4. Otherwise solve with the extended Euclidean algorithm; the combined modulus is $\text{lcm}(m_1, m_2)$

**Combining many congruences** `crt_many`: a left fold — start from the first congruence and combine with the next one at a time:

$$x \equiv r_1 \pmod{m_1} \xrightarrow{\text{combine}} x \equiv R_2 \pmod{\text{lcm}(m_1, m_2)} \xrightarrow{\text{combine}} \cdots$$

The final result is $(R, M)$ where $M = \text{lcm}(m_1, \ldots, m_k)$ and $0 \leq R < M$. If any step detects an inconsistency, the whole computation returns `None`.

**Complexity**: $k$ pairwise combinations, each involving a GCD computation of cost $O(\log M)$.

### Quadratic Residues

#### The Legendre Symbol and the Jacobi Symbol

**Legendre symbol**: for an odd prime $p$,

$$\left(\frac{a}{p}\right) = \begin{cases} 0 & p \mid a \\ 1 & a \text{ is a quadratic residue mod } p \\ -1 & a \text{ is a quadratic non-residue mod } p \end{cases}$$

**Euler's criterion**: $\left(\frac{a}{p}\right) \equiv a^{(p-1)/2} \pmod{p}$.

**Jacobi symbol**: the generalization of the Legendre symbol to composite moduli. If $n = p_1^{e_1} \cdots p_k^{e_k}$, then

$$\left(\frac{a}{n}\right) = \prod_i \left(\frac{a}{p_i}\right)^{e_i}$$

**Key property**: the Jacobi symbol can be computed quickly via the **law of quadratic reciprocity**, without factoring $n$.

#### The Law of Quadratic Reciprocity

For odd primes $p, q$:

$$\left(\frac{p}{q}\right) \left(\frac{q}{p}\right) = (-1)^{\frac{p-1}{2} \cdot \frac{q-1}{2}}$$

Supplementary rules:

- $\left(\frac{2}{n}\right) = (-1)^{(n^2 - 1)/8}$, i.e. $1$ when $n \equiv \pm 1 \pmod 8$, and $-1$ when $n \equiv \pm 3 \pmod 8$
- $\left(\frac{-1}{n}\right) = (-1)^{(n-1)/2}$, i.e. $1$ when $n \equiv 1 \pmod 4$

**Computational algorithm** (oCAS's `jacobi` function):

```
jacobi(a, n):
    if n is even or non-positive, return 0
    repeat:
      1. 2-adic stripping: extract all factors of 2 from a, accumulating the sign with the mod-8 rule
      2. reciprocity flip: swap a ↔ n, adjusting the sign according to the mod-4 values of a and n
      3. a = a mod n
      4. if a == 0, return 0 (n > 1) or the accumulated sign (n == 1)
```

**Complexity**: $O(\log^2 n)$, the same order as the Euclidean algorithm.

#### The Tonelli–Shanks Algorithm

Find $x$ such that $x^2 \equiv a \pmod{p}$ ($p$ an odd prime, $\left(\frac{a}{p}\right) = 1$).

**Fast path**: when $p \equiv 3 \pmod{4}$,

$$x = a^{(p+1)/4} \bmod p$$

because $x^2 = a^{(p+1)/2} = a \cdot a^{(p-1)/2} = a \cdot \left(\frac{a}{p}\right) = a$.

**The general case** ($p \equiv 1 \pmod{4}$) — the Tonelli–Shanks algorithm:

1. **Decompose**: write $p - 1 = Q \cdot 2^S$ ($Q$ odd)
2. **Find a non-residue**: search for $z$ with $\left(\frac{z}{p}\right) = -1$
3. **Initialize**:
   - $M = S$
   - $c = z^Q \bmod p$ (an element of order $2^{S-1}$)
   - $t = a^Q \bmod p$
   - $R = a^{(Q+1)/2} \bmod p$
4. **Loop** (while $t \neq 1$):
   - find the smallest $i > 0$ such that $t^{2^i} \equiv 1 \pmod{p}$
   - update $c \leftarrow c^{2^{M-i-1}}$, $t \leftarrow t \cdot c^2$, $R \leftarrow R \cdot c$, $M \leftarrow i$
5. Return $R$

**Invariant**: $R^2 = a \cdot t \pmod{p}$ and $c$ has order $2^M$. When the loop terminates $t = 1$, hence $R^2 \equiv a$.

**Complexity**: $O(S^2 \log p)$, where $S = v_2(p - 1)$ (the power of 2 dividing $p - 1$). On average $S$ is small, so the algorithm is efficient.

### Rational Reconstruction

**Problem**: given $a \in \mathbb{Z}$ and a modulus $m$, find $n, d \in \mathbb{Z}$ such that

$$a \cdot d \equiv n \pmod{m}, \quad \gcd(n, d) = 1, \quad 2|n| \cdot |d| < m$$

**Application**: in the modular GCD algorithm, after combining the GCD coefficients from several $\mathbb{F}_p$'s with CRT, the result modulo $M$ must be "reconstructed" back to rational numbers over $\mathbb{Q}$.

**The Wang/extended Euclidean algorithm**:

Use the extended Euclidean algorithm tracking the sequence $(r_i, t_i)$:

1. Initialize: $(r_0, r_1) = (m, a)$, $(t_0, t_1) = (0, 1)$
2. Iterate: $q = \lfloor r_0 / r_1 \rfloor$, $(r_0, r_1) \leftarrow (r_1, r_0 - q \cdot r_1)$, $(t_0, t_1) \leftarrow (t_1, t_0 - q \cdot t_1)$
3. **Termination condition**: stop when $|r_1| \leq \sqrt{m/2}$ **and** $|t_1| \leq \sqrt{m/2}$
4. **Verification**: let $n = r_1$ and $d = t_1$ (take positive values); check $a \cdot d \equiv n \pmod{m}$ and $2|n| \cdot |d| < m$
5. If verification fails or $t_1 = 0$, return `None` (no rational reconstruction satisfies the conditions)

**Uniqueness theorem**: when $2|n| \cdot |d| < m$, the pair $(n, d)$ satisfying the conditions is unique. This condition guarantees that the reconstructed result is the "simplest" rational representation.

**Complexity**: $O(\log m)$ GCD steps, the same order as the Euclidean algorithm.

## Implementation in oCAS

oCAS's number theory algorithms live in the `number_theory` module of the `ocas-domain` crate and in the rational reconstruction module of `ocas-poly`.

### Module Structure

```
ocas-domain/src/
├── number_theory.rs          ← primality tests, modular inverse, CRT, Legendre/Jacobi, mod_sqrt
└── number_theory/
    ├── primes.rs             ← BPSW, Lucas sequences (lucas_uv_mod, strong_lucas_prp)
    ├── factor.rs             ← the integer factorization driver (factor_integer) and its methods
    ├── dlog.rs               ← BSGS and Pohlig–Hellman
    ├── crt.rs                ← multi-congruence CRT (crt_many)
    └── functions.rs          ← φ, μ, τ, σ, λ

ocas-poly/src/
└── rational_reconstruction.rs ← rational reconstruction (extended Euclidean method)
```

### Primality Testing

| Function | Location | Description |
|---|---|---|
| `is_prime(n)` | `number_theory.rs` | 12-base Miller–Rabin; deterministic for $n < 3.317 \times 10^{24}$ |
| `is_prime_bpsw(n)` | `primes.rs` | base-2 MR + strong Lucas PRP |
| `is_prime_u64(n)` | `primes.rs` | `u64`-specific; delegates to `is_prime` |
| `next_prime(n)` | `number_theory.rs` | test odd numbers one by one starting from $n + 1$ |
| `primes_from(n)` | `number_theory.rs` | a prime iterator |

Lucas sequence computation is implemented by `lucas_uv_mod` using the binary ladder, returning $(U_k, V_k, Q^k) \bmod n$. This function is shared by the BPSW test and Williams's $p + 1$ factorization method.

### Integer Factorization

| Function | Location | Description |
|---|---|---|
| `factor_integer(n)` | `factor.rs` | entry point: trial division + recursive splitting |
| `factor_integer_with_rng(n, rng)` | `factor.rs` | version with an explicit RNG |
| `factor_trial(n, limit)` | `factor.rs` | trial division, returns `(factor, cofactor)` |
| `pollard_rho_brent(n, rng)` | `factor.rs` | Brent variant, with retries |
| `pollard_pm1(n, b1, rng)` | `factor.rs` | Pollard $p - 1$, Stage 1 |
| `williams_pp1(n, b1, rng)` | `factor.rs` | Williams $p + 1$, Stage 1 |
| `ecm(n, b1, max_curves, rng)` | `factor.rs` | Lenstra ECM, Suyama parameterization |

Internal implementation details:

- `primes_up_to(limit)`: the sieve of Eratosthenes, generating prime tables for trial division and smoothness bounds
- `prime_power_le(q, bound)`: computes $q^e \leq \text{bound}$, used to construct the smooth exponent $M$
- `ProjPoint { x, z }`: a Montgomery projective-coordinate point
- `ecm_double`, `ecm_add`, `ecm_mul`: group operations on the Montgomery curve
- `suyama_curve(sigma, n)`: the Suyama parameterization, returning the `Suyama` enum (`Curve`/`Factor`/`Degenerate`)
- `find_factor(n, rng)`: the escalation driver (rho → p−1 → p+1 → ECM, with increasing smoothness bounds)

### Discrete Logarithms

| Function | Location | Description |
|---|---|---|
| `dlog_bsgs(base, target, modulus)` | `dlog.rs` | baby-step giant-step |
| `dlog_pohlig_hellman(base, target, p)` | `dlog.rs` | Pohlig–Hellman (requires $p$ prime) |

Internals:

- `bsgs_bounded(base, target, modulus, order_bound)`: the BSGS core, searching with a given order bound
- `dlog_bsgs` sets the bound to `modulus - 1` and delegates to `bsgs_bounded`
- `dlog_pohlig_hellman` factors $p - 1$ with `factor_integer`, recovers the digits prime power by prime power with BSGS, and combines them with `crt_many`

### The Chinese Remainder Theorem

| Function | Location | Description |
|---|---|---|
| `crt(r1, m1, r2, m2)` | `number_theory.rs` | pairwise combination; the moduli need not be coprime |
| `crt_many(congruences)` | `crt.rs` | left-fold combination of many congruences |

`crt_many` starts from the first congruence and combines with the next one at a time via `crt`. If any step is inconsistent, the whole computation returns `None`.

### Quadratic Residues

| Function | Location | Description |
|---|---|---|
| `legendre(a, p)` | `number_theory.rs` | the Legendre symbol, delegating to `jacobi` |
| `jacobi(a, n)` | `number_theory.rs` | the Jacobi symbol (quadratic reciprocity + 2-adic stripping) |
| `mod_sqrt(a, p)` | `number_theory.rs` | Tonelli–Shanks, including the $p \equiv 3 \pmod{4}$ fast path |

The `mod_sqrt` implementation:

1. First check $\left(\frac{a}{p}\right) = 1$ (otherwise return `None`)
2. If $p \equiv 3 \pmod{4}$, take the fast path $x = a^{(p+1)/4} \bmod p$
3. Otherwise run the full Tonelli–Shanks algorithm: find a non-residue $z$, decompose $p - 1 = Q \cdot 2^S$, and loop maintaining the $(c, t, r, m)$ invariants

### Auxiliary Number-Theoretic Functions

| Function | Location | Description |
|---|---|---|
| `euler_phi(n)` | `functions.rs` | $\varphi(n) = |n| \prod(1 - 1/p)$, based on `factor_integer` |
| `moebius_mu(n)` | `functions.rs` | $\mu(n)$, based on `factor_integer` |
| `divisor_tau(n)` | `functions.rs` | $\tau(n) = \prod(e_i + 1)$ |
| `divisor_sigma(n, k)` | `functions.rs` | $\sigma_k(n) = \prod \frac{p^{k(e+1)} - 1}{p^k - 1}$ |
| `liouville_lambda(n)` | `functions.rs` | $\lambda(n) = (-1)^{\Omega(n)}$ |

All multiplicative functions first obtain the prime factorization via `factor_integer` and then evaluate the multiplicative formula.

### Rational Reconstruction

| Function | Location | Description |
|---|---|---|
| `rational_reconstruction(a, m)` | `ocas-poly/src/rational_reconstruction.rs` | the Wang/extended Euclidean method |

Internally it uses `integer_sqrt` to compute $\lfloor\sqrt{m/2}\rfloor$ as the threshold for the termination bound.

### Basic Utility Functions

| Function | Location | Description |
|---|---|---|
| `mod_inv(a, m)` | `number_theory.rs` | modular inverse via the extended Euclidean algorithm |
| `extended_gcd(a, b)` | `number_theory.rs` | returns $(g, x, y)$ with $ax + by = g$ |
| `symmetric_mod(a, m)` | `number_theory.rs` | reduce into $(-m/2, m/2]$ |

## References

1. **Shoup, V.** *A Computational Introduction to Number Theory and Algebra.* Cambridge University Press, 2nd edition, 2009.
   - Ch.10: Primality testing (Miller–Rabin, randomized algorithms)
   - Ch.11: Finding discrete logarithms (BSGS, Pohlig–Hellman)
   - Ch.19: Factoring integers (Pollard rho, ECM)

2. **Crandall, R. & Pomerance, C.** *Prime Numbers: A Computational Perspective.* Springer, 2nd edition, 2005.
   - Ch.3: Recognizing primes and composites (BPSW, Lucas tests)
   - Ch.5: Exponential factoring algorithms (ECM, Pollard rho)
   - Ch.6: Subexponential factoring (Pollard $p-1$, Williams $p+1$)
   - Ch.7: Modern discrete logarithm algorithms

3. **Brent, R. P.** "An improved Monte Carlo factorization algorithm." *BIT Numerical Mathematics*, 20(2):176–184, 1980.
   - The Brent variant of the Pollard rho algorithm

4. **Williams, H. C.** "A p+1 method of factoring." *Mathematics of Computation*, 39(159):225–234, 1982.
   - Williams's $p + 1$ method

5. **Lenstra, H. W. Jr.** "Factoring integers with elliptic curves." *Annals of Mathematics*, 126(3):649–673, 1987.
   - The ECM method

6. **Montgomery, P. L.** "Speeding the Pollard and elliptic curve methods of factorization." *Mathematics of Computation*, 48(177):243–264, 1987.
   - Montgomery curve coordinates and ladder optimizations

7. **Shanks, D. & Tonelli, R.** See Cohen, H. *A Course in Computational Algebraic Number Theory*, Algorithm 1.5.1.
   - The Tonelli–Shanks modular square root algorithm

8. **Baillie, R. & Wagstaff, S. S.** "Lucas pseudoprimes." *Mathematics of Computation*, 35(152):1391–1417, 1980.
   - The Lucas part of the BPSW test

9. **Pohlig, S. & Hellman, M.** "An improved algorithm for computing logarithms over GF(p) and its cryptographic significance." *IEEE Transactions on Information Theory*, 24(1):106–110, 1978.
   - The Pohlig–Hellman discrete logarithm algorithm

10. **Wang, P. S.** "A p-adic algorithm for univariate partial fractions." *Proceedings of SYMSAC '81*, 212–217, 1981.
    - The extended Euclidean method for rational reconstruction
