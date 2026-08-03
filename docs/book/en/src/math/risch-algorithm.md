# The Risch Integration Algorithm

## Prerequisites

- [Polynomial Algebra](./polynomial-algebra.md) — polynomial rings, factorization, resultants
- [Finite Fields & Modular Arithmetic](./finite-fields.md) — modular inverses and the extended Euclidean algorithm
- [Linear Algebra](./linear-algebra.md) — Gaussian elimination and determinants
- [Symbolic Calculus](./symbolic-calculus.md) — differentiation rules and the chain rule
- [Polynomial GCD & Factorization](./poly-gcd-factoring.md) — square-free factorization and Hensel lifting

---

## Basic Concepts

### The Precise Definition of Elementary Functions

The term "elementary function" is usually used intuitively in calculus textbooks. The Risch algorithm requires a precise definition:

**Definition 1** (Field of constants). Let $k$ be a differential field of characteristic zero with derivation denoted by $D$. The most basic field of constants is $\mathbb{Q}$, on which $D = 0$.

**Definition 2** (Elementary extension). Let $k$ be a differential field and $K$ an extension field of $k$. $K/k$ is called an **elementary extension** if $K$ is generated from $k$ by a finite chain of extensions of the following three types:

1. **Algebraic extension**: $K = k(\alpha)$, where $\alpha$ is algebraic over $k$ (i.e. $\alpha$ is a root of some polynomial in $k[x]$). For example $\mathbb{Q}(\sqrt{2})/\mathbb{Q}$.
2. **Logarithmic extension**: $K = k(t)$, where $Dt = Du / u$ for some $u \in k^*$, i.e. $t = \log(u)$. Derivation rule: $D(\log u) = Du/u$.
3. **Exponential extension**: $K = k(t)$, where $Dt / t = Du$ for some $u \in k$, i.e. $t = \exp(u)$. Derivation rule: $D(\exp u) = \exp(u) \cdot Du$.

**Definition 3** (Elementary function). An **elementary function** over a field $k$ is an element of some elementary extension tower $k = k_0 \subset k_1 \subset \cdots \subset k_n$ in which each $k_{i+1}/k_i$ is one of the three types above.

**Theorem** (Liouville, 1835). Let $f$ be an elementary function over $k$. Then $f$ has an elementary antiderivative over $k$ if and only if it can be written as

$$f = Dv + \sum_{i=1}^{n} c_i \cdot \frac{Du_i}{u_i}$$

where $v, u_1, \dots, u_n \in k$ and $c_1, \dots, c_n$ are $D$-constants (i.e. $Dc_i = 0$).

The deep significance of this theorem is that if an elementary function has an antiderivative, the only "new" parts of it can be constant multiples of logarithms. This is the theoretical foundation of the Risch algorithm.

### Differential Fields

**Definition 4** (Differential field). A **differential field** is a field $k$ equipped with a derivation $D: k \to k$ satisfying:

- **Additivity**: $D(a + b) = Da + Db$
- **Leibniz rule**: $D(a \cdot b) = a \cdot Db + Da \cdot b$

The kernel of the derivation $\{c \in k : Dc = 0\}$ is called the **field of constants**, denoted $C_k$ or $k^D$.

In the context of the Risch algorithm:

- The base field is $k_0 = \mathbb{Q}(x)$ with $D = d/dx$ (differentiation with respect to the integration variable $x$); its field of constants is $C_{k_0} = \mathbb{Q}$.
- A logarithmic extension $k(t) = k(\log u)$ has $Dt = Du/u$.
- An exponential extension $k(t) = k(\exp u)$ has $Dt = t \cdot Du$.

**Key fact**: in an elementary extension tower the field of constants stays the same: $C_{k_n} = C_{k_0} = \mathbb{Q}$. This property is used implicitly by the algorithm — when we say "rational constant", we mean an element of $\mathbb{Q}$.

---

## Core Theory

### Constructing the Differential-Field Tower

Given an integrand $f$ (an elementary function expression), the first step of the Risch algorithm is to construct the differential-field tower in which it "lives".

**Algorithm** `build_tower(f, var)`:

1. **Collect function symbols**: scan the expression tree of $f$ and collect all function applications (`log`, `exp`).
2. **Order them**: sort by dependency — if $t_2 = \exp(\log(x) + 1)$, then $t_1 = \log(x)$ must come before $t_2$.
3. **Detect algebraic dependence**: for each candidate generator $t_i$, check whether it is algebraically dependent over the tower built so far. The conservative strategy is:
   - $\log(c \cdot u)$ and $\log(u)$ (with $c$ a rational constant) → algebraically dependent ($\log(c \cdot u) = \log(u) + \log(c)$)
   - $\exp(u + c)$ and $\exp(u)$ (with $c$ a rational constant) → algebraically dependent ($\exp(u + c) = e^c \cdot \exp(u)$)
   - otherwise reject (conservative strategy: rather reject a true dependence than merge incorrectly)
4. **Reject non-integral powers**: if $f$ contains non-integral powers such as $\sqrt{x}$ (i.e. algebraic functions), reject (return `None`).
5. **Compute derivatives**: for each generator $t_i$, compute $Dt_i$ (with respect to $D = d/dx$):
   - $t_i = \log(u)$: $Dt_i = Du / u$ (compute $Du$ recursively)
   - $t_i = \exp(u)$: $Dt_i = t_i \cdot Du$

In oCAS, `build_tower` lives in `ocas-calc/src/tower/build.rs`. It returns a `Tower` struct containing the list of generators `gens: Vec<GenInfo>` (each with kind, atom reference, and derivative) and the field-of-constants information.

**Representation of field elements**: elements of the tower are represented by numerator/denominator pairs of **sparse multivariate polynomials** (`SparseMultivariatePolynomial<RationalDomain, Lex>`):

- `KElem` (k-element): `{ num: Sparse, den: Sparse }`, representing an element $p/q$ of $k_\ell$
- `KPoly` (k-polynomial): a univariate polynomial in the top-level generator $t_\ell$ with coefficients in `KElem`
- `KRat` (k-rational function): `{ num: KPoly, den: KPoly }`, representing a rational function in $k_\ell(t_\ell)$

**Example.** For $f = x \cdot \exp(x) + \log(x)$:

| Level | Field | Generator | Derivation |
|---|---|---|---|
| $k_0$ | $\mathbb{Q}(x)$ | $x$ | $Dx = 1$ |
| $k_1$ | $k_0(t_1)$, $t_1 = \log(x)$ | $t_1$ | $Dt_1 = 1/x$ |
| $k_2$ | $k_1(t_2)$, $t_2 = \exp(x)$ | $t_2$ | $Dt_2 = t_2$ |

### Level-by-Level Integration: The Overall Framework

The core of the Risch algorithm is **top-down recursion through the levels**. At each level $\ell$ of the tower, the integrand $f$ is viewed as a rational function in $k_{\ell-1}(t_\ell)$.

**Algorithm** `integrate_level(tower, level, f)`:

Write $f = a/d \in k_\ell(t_\ell)$ with $a, d \in k_\ell[t_\ell]$ and $\gcd(a, d) = 1$.

1. **Hermite reduction** (split off the rational part): decompose $a/d$ as

$$\frac{a}{d} = Dg + \frac{a_1}{d_1}$$

where $g \in k_\ell(t_\ell)$ has an explicit antiderivative and $d_1$ is square-free.

2. **Integrate the polynomial part**: if the polynomial part of $a/d$ is $p(t_\ell) \in k_\ell[t_\ell]$, integrate it according to the type of the generator:
   - **Primitive level** ($t_\ell = \log u$): the method of undetermined coefficients
   - **Hyperexponential level** ($t_\ell = \exp u$): the Risch differential equation

3. **Integrate the logarithmic part**: for $a_1/d_1$ ($d_1$ square-free), check whether it matches the **logarithmic derivative identity**.

4. **Base field** $\mathbb{Q}(x)$: delegate to the rational-function integrator.

The results are assembled into `LevelResult { elem, logs, extras }`, where `elem` is the field-element part, `logs` is the logarithmic part of the form $\sum c_i \log(v_i)$, and `extras` is the part that could not be integrated.

### Hermite Reduction

**Goal**: decompose a rational function $a/d$ into an "already integrated" part plus a "simpler" remainder (the denominator of the remainder is square-free).

**Algorithm** (Hermite reduction as in Bronstein ch. 2; oCAS uses the square-free-factorization recursive form):

Let $a/d$ be a proper fraction ($\deg a < \deg d$) with $d$ monic.

1. Compute the square-free factorization $d = \prod_j v_j^{m_j}$. If the maximum multiplicity is $m \leq 1$, then $d$ is already square-free; return $(0, a, d)$.
2. Pick a factor $v$ of maximum multiplicity ($= m \geq 2$) and write $d = u \cdot v^m$.
3. Let $B = u \cdot Dv$ (the derivative under the tower derivation; over the base field this is $u \cdot v'$). Use the extended Euclidean algorithm to find $s, t$ with $s \cdot B + t \cdot v = a$ ($s$ reduced modulo $v$).
4. Recurrence formula:
$$\frac{a}{u\, v^m} = D\!\left(\frac{-s}{(m-1)\, v^{m-1}}\right) + \frac{t + u \cdot D\bigl(s/(m-1)\bigr)}{u\, v^{m-1}}$$
5. Recurse on the remainder $\dfrac{t + u \cdot D(s/(m-1))}{u\, v^{m-1}}$ until the denominator is square-free.

Here $D\!\left(\frac{-s}{(m-1) v^{m-1}}\right)$ has the explicit antiderivative $\frac{-s}{(m-1) v^{m-1}}$ itself, and the denominator $d_1$ of the final remainder $a_1/d_1$ is square-free. Each iteration lowers the maximum multiplicity by one, so at most $m$ iterations are needed.

In oCAS, Hermite reduction is implemented in two places:

- `ocas-calc/src/integral/rational.rs`: Hermite reduction over the base field $\mathbb{Q}(x)$ (`hermite_reduce`)
- `ocas-calc/src/integral/risch.rs`: Hermite reduction at an arbitrary tower level (`hermite_tower`)

**Example**. $\int \frac{1}{(x+1)^2}\, dx$:

$d = (x+1)^2$; its square-free factorization is $(x+1)^2$, so the maximum multiplicity is $m = 2$, $v = x+1$, and $u = 1$. Then $B = u \cdot Dv = 1$, and solving $s \cdot B + t \cdot v = 1$ gives $s = 1$ (mod $v$), $t = 0$. The recurrence yields:

$$\frac{1}{(x+1)^2} = D\!\left(\frac{-1}{x+1}\right) + \frac{0}{x+1}$$

i.e. $\int = -1/(x+1)$. More generally, $\int a/(x+1)^2$ can always be handled by Hermite reduction.

### The Logarithmic Derivative Identity

**Core idea**: a term of the form $c \cdot \frac{Du}{u}$ has the trivial antiderivative $c \cdot \log(u)$.

**Theorem** (logarithmic derivative identity). Let $f \in k_\ell(t_\ell)$ and let $d_1$ be square-free. If

$$\frac{a_1}{d_1} = c \cdot \frac{D d_1}{d_1}$$

for some $D$-constant $c$, then

$$\int \frac{a_1}{d_1} = c \cdot \log(d_1)$$

More generally, for a square-free denominator $d_1$, all contributions to the logarithmic part have the form $\sum_i c_i \log(v_i)$, where the $v_i$ are given by the irreducible factors of $d_1$.

**How to check it**: over the base field $\mathbb{Q}(x)$, check whether $a_1$ equals some rational constant $c$ times $d_1'$. At higher tower levels, check whether $a_1$ is a constant of the coefficient field times $d_1'$.

**A more precise method** (Rothstein–Trager resultant): for a general square-free denominator $d(x)$, the logarithmic part of $\int a(x)/d(x)\,dx$ is given by the **Rothstein–Trager resultant**:

$$R(t) = \text{Res}_x\!\bigl(d(x),\; a(x) - t \cdot d'(x)\bigr)$$

If the roots of $R(t)$ over $\mathbb{Q}$ are $c_1, \dots, c_m$, then

$$\int \frac{a}{d}\,dx = \sum_{i=1}^{m} c_i \cdot \log\gcd(d,\; a - c_i \cdot d') + C$$

In the oCAS implementation (`rothstein_trager`), $R(t)$ is computed by **interpolation** — evaluate the resultant at $t = 0, 1, 2, \dots$ and then recover the polynomial by Lagrange interpolation. When $R(t)$ does not split completely over $\mathbb{Q}$, the corresponding terms are returned in the unevaluated form `Integral(term, var)`.

**The special case of a quadratic denominator**: when $d(x) = x^2 + bx + c$ (irreducible), completing the square gives:

$$\int \frac{Ax + B}{x^2 + bx + c}\,dx = \frac{A}{2}\log(x^2 + bx + c) + \frac{2B - Ab}{\sqrt{4c - b^2}}\arctan\!\left(\frac{2x + b}{\sqrt{4c - b^2}}\right)$$

When $4c - b^2 < 0$ (real-root case), $\arctan$ becomes $\text{artanh}$ (inverse hyperbolic tangent), corresponding to the logarithmic form.

### The Risch Differential Equation

The central subproblem of the Risch algorithm is solving the **Risch differential equation** (RDE).

**Problem**. Given elements $f, g$ of the differential field $k_\ell$ (neither containing the top-level variable $t_\ell$), find $q \in k_\ell[t_\ell]$ satisfying:

$$Dq + f \cdot q = g$$

where $D$ is the tower derivation (total derivative with respect to $t_\ell$). Note that the coefficients of $f$ and $g$ lie in $k_\ell$ and $q$ is a polynomial in $k_\ell[t_\ell]$.

**Why only polynomial solutions?** Rational-function solutions $q = p/d$ require an additional denominator-bound analysis, a piece not covered by the current oCAS implementation. When `None` is returned, the caller falls back to other integration methods.

#### RDE over the Base Field $\mathbb{Q}(x)$

Over the base field $k_0 = \mathbb{Q}(x)$, the RDE becomes the ordinary differential equation $q'(x) + f(x) \cdot q(x) = g(x)$.

**Algorithm** (`base_rde`):

1. **Degree bound**: let $\deg f = p$ and $\deg g = r$. If $f = 0$, integrate directly, giving $\deg q = r + 1$; if $f$ is a nonzero constant ($p = 0$), then $\deg q = r$; if $p \geq 1$, the leading term of $f \cdot q$ (of degree $p + \deg q$) dominates $\deg(q') \leq \deg q - 1$, so $\deg q = r - p$ (uniquely determined).
2. **Undetermined coefficients**: write $q = \sum_{i=0}^{n} a_i x^i$ and substitute into the equation.
3. **Eliminate from high to low degree**: compare the coefficients of the powers of $x$ and determine the $a_i$ one by one. Each $a_i$ is either determined uniquely or contradictory (the antiderivative does not exist).

**Example**. Solve $q' + q = x$ ($f = 1, g = x$):

Write $q = a_1 x + a_0$. Then $q' = a_1$. Substituting: $a_1 + a_1 x + a_0 = x$.

Comparing coefficients: $x^1: a_1 = 1$; $x^0: a_1 + a_0 = 0 \Rightarrow a_0 = -1$.

Solution: $q = x - 1$. Verification: $(x-1)' + (x-1) = 1 + x - 1 = x$ ✓

This gives $\int x \cdot e^x\,dx = (x - 1) e^x$ (at the hyperexponential level).

#### RDE at a Primitive Level

At a level with $t_\ell = \log(u)$, the derivation satisfies $Dt_\ell = Du/u \in k_{\ell-1}$ (it does not contain $t_\ell$).

**Property**: $D(a_0 + a_1 t + \cdots + a_m t^m) = Da_0 + Da_1 \cdot t + \cdots + Da_m \cdot t^m + (a_1 + 2a_2 t + \cdots + m a_m t^{m-1}) \cdot Dt$

Since $Dt$ does not contain $t$ (the key property of a primitive extension), applying $D$ to $k[t]$ **does not change the degree in $t$** (except that the $\partial/\partial t$ part lowers the degree by one).

**Algorithm**: start from the top coefficient $a_m$ of $q$ and **eliminate top-down**:

1. Determine $a_m$ from the coefficient of $t^m$ in the equation (an RDE for $a_m$ in $k_{\ell-1}$).
2. Substitute and lower the degree; repeat for $a_{m-1}$.
3. Finally $a_0$ must satisfy a **logarithmic constraint**: the constant part of $a_0$ must make the logarithmic terms arising in the lower-level integration consistent.

**Recursive structure**: each step produces an RDE over $k_{\ell-1}$; call `rde_solve` recursively down to the base field.

#### RDE at a Hyperexponential Level

At a level with $t_\ell = \exp(u)$, the derivation satisfies $Dt_\ell = t_\ell \cdot Du$.

**Key property**: $Dt = t \cdot Du$ contains $t$, so applying $D$ to $k[t]$ **mixes degrees**:

$$D(a_i t^i) = (Da_i) t^i + a_i \cdot i \cdot t^i \cdot Du = (Da_i + i \cdot a_i \cdot Du) t^i$$

i.e. $D(a_i t^i) = (Da_i + i \cdot a_i \cdot Du) \cdot t^i$ — each power $t^i$ is **independent**!

**Decoupling property**: substituting $q = \sum a_i t^i$ and $g = \sum b_j t^j$ into $Dq + fq = g$ and comparing the coefficient of $t^k$:

$$Da_k + k \cdot a_k \cdot Du + f \cdot a_k = b_k$$

i.e.

$$Da_k + (f + k \cdot Du) \cdot a_k = b_k$$

This is an RDE over $k_{\ell-1}$ (unknown $a_k$), and the equations for different $k$ are **independent**!

**Algorithm**: for each $k = 0, 1, \dots, \deg g$, solve the RDE

$$Da_k + (f + k \cdot Du) \cdot a_k = b_k$$

independently. If some $a_k$ has no solution, the whole equation has no solution. An upper bound for $\deg q$ is given by $\deg g$.

**Example**. $\int x \cdot e^x\,dx$: on the tower with $t = e^x$ ($Dt = t$), the integrand is $x \cdot t$.

The polynomial part is $p(t) = x \cdot t$ ($m = 1$). Comparing the $t^1$ coefficient: $Da_1 + a_1 \cdot 1 = x$, i.e. $a_1' + a_1 = x$. The base-field RDE gives $a_1 = x - 1$. The $t^0$ part is 0.

Result: $\int = (x-1) \cdot t = (x-1) \cdot e^x$.

### The Method of Undetermined Coefficients at Logarithmic Levels

When integrating a polynomial $p(t) = \sum a_i t^i$ at a primitive (logarithmic) level $t = \log(u)$, the **method of undetermined coefficients** is used.

**Algorithm** `integrate_kpoly_primitive`:

Let the polynomial to be integrated be $p(t) = \sum_{i=0}^{m} a_i t^i$ and guess the antiderivative $q(t) = \sum_{i=0}^{m'} b_i t^i$.

1. Set $m' = m + 1$ (the degree of the antiderivative is at most one higher than that of the integrand, because the $\partial/\partial t$ part of $D$ lowers the degree).
2. Compute $Dq = \sum (Db_i) t^i + \sum i \cdot b_i \cdot (Du/u) \cdot t^{i-1}$.
3. Compare the coefficients of $p = Dq$:
   - $t^m$: $a_m = Db_m + (m+1) b_{m+1} \cdot Du/u$
   - $\vdots$
   - $t^0$: $a_0 = Db_0 + b_1 \cdot Du/u$
4. Solve downwards starting from $t^{m+1}$. Each $b_i$ may require a recursive call to the integrator (for the $Db_i$ part).
5. The determination of $b_0$ is constrained by the **logarithmic constraint**: in the final result, all logarithmic terms $\sum c_j \log(v_j)$ must be consistent. If $b_0$ contains logarithmic terms coming from lower levels that do not satisfy the Liouville condition, there is no solution.

### Rewriting Trigonometric Functions to Complex Exponentials

The Risch algorithm natively handles only `log` and `exp`. Trigonometric functions must first be rewritten in complex-exponential form.

**Euler's formula**:

$$\sin(u) = \frac{e^{iu} - e^{-iu}}{2i}, \qquad \cos(u) = \frac{e^{iu} + e^{-iu}}{2}$$

$$\tan(u) = \frac{\sin u}{\cos u} = \frac{e^{iu} - e^{-iu}}{i(e^{iu} + e^{-iu})}$$

In oCAS, `trig_to_exp` (`ocas-calc/src/integral/trig.rs`) walks the expression tree and replaces `sin`, `cos`, `tan`, `cot`, `sec`, `csc` one by one with the equivalent forms above. The imaginary unit $I$ is added to the tower as a constant generator ($DI = 0$).

**Integrating the rewritten form**: the rewritten expression lives in the differential-field tower over $\mathbb{Q}(x, I, e^{Ix}, \dots)$ and can be handled by the Risch algorithm.

**Conversion back to real form** (`realify`): the result of a Risch integration may contain the imaginary unit. `realify` attempts to convert the result back into real form:

- **Merging conjugate logarithms**: if the result contains $c \cdot \log(u + Iv) + c \cdot \log(u - Iv)$, merge it into $c \cdot \log(u^2 + v^2)$
- **Difference of conjugate logarithms**: if it contains $c \cdot \log(u + Iv) - c \cdot \log(u - Iv)$, merge it into $2c \cdot \arctan(v/u)$
- **Merging exponentials**: $e^{Iu} \cdot e^{-Iu} = 1$, $e^{Iu} + e^{-Iu} = 2\cos(u)$, etc.

This is a "best-effort" process. If no known pattern matches, the complex form is kept (it is still mathematically correct, since differentiation-verification still holds).

**Current limitation**: the Risch differential equation solver works over $\mathbb{Q}[x]$, so when trigonometric rewriting produces hyperexponential equations containing $I$ (e.g. the equations arising when integrating $\sin(x)\cos(x)$ or $\cos^2(x)$), the coefficient field extends to $\mathbb{Q}(I)$ rather than $\mathbb{Q}$, and the solver may fail. Such integrands are returned in unevaluated form.

### The Heuristic Layer

Before the Risch algorithm, oCAS tries a set of **heuristic integration techniques** (`ocas-calc/src/integral/heuristic.rs`). These techniques are fast and cover common patterns.

#### Integration by Parts (LIATE Ordering)

For a product $f \cdot g$, try integration by parts $\int u\,dv = uv - \int v\,du$.

**LIATE priority** (heuristic rule for choosing $u$):

| Priority | Type | Score | Example |
|---|---|---|---|
| 1 | **L**ogarithmic | 0 | $\log x$, $\arctan x$ |
| 2 | **I**nverse trig | 1 | $\arcsin x$, $\text{arccosh}\, x$ |
| 3 | **A**lgebraic | 2 | $x^2$, $\sqrt{x}$ |
| 4 | **T**rigonometric | 3 | $\sin x$, $\cos x$ |
| 5 | **E**xponential | 4 | $e^x$, $2^x$ |

**Algorithm**: for each factor of the product, sort by LIATE score. Choose the factor with the lowest score (highest priority) as $u$, and the rest (integrated to obtain $v$) as $dv$.

**Depth limit**: `PARTS_MAX_DEPTH = 2`, to prevent infinite recursion.

#### Trigonometric Substitution

Match $\sqrt{a^2 - x^2}$, $\sqrt{a^2 + x^2}$, $\sqrt{x^2 - a^2}$ and their reciprocals, and directly return the known antiderivatives:

| Integrand | Antiderivative |
|---|---|
| $\frac{1}{\sqrt{a^2 - x^2}}$ | $\arcsin(x/a)$ |
| $\frac{1}{\sqrt{a^2 + x^2}}$ | $\text{arcsinh}(x/a)$ |
| $\frac{1}{\sqrt{x^2 - a^2}}$ | $\text{arccosh}(x/a)$ |
| $\sqrt{a^2 - x^2}$ | $\frac{x\sqrt{a^2-x^2} + a^2 \arcsin(x/a)}{2}$ |
| $\sqrt{a^2 + x^2}$ | $\frac{x\sqrt{a^2+x^2} + a^2 \text{arcsinh}(x/a)}{2}$ |
| $\sqrt{x^2 - a^2}$ | $\frac{x\sqrt{x^2-a^2} - a^2 \text{arccosh}(x/a)}{2}$ |

#### Weierstrass Substitution

For rational functions of $\sin(u)$ and $\cos(u)$ (with $u$ a linear function of $x$), use the universal substitution $t = \tan(u/2)$:

$$\sin u = \frac{2t}{1+t^2}, \qquad \cos u = \frac{1-t^2}{1+t^2}, \qquad du = \frac{2\,dt}{1+t^2}$$

This reduces the integral of a rational function of trigonometric functions to the integral of a rational function in $t$.

#### Euler Substitution

For integrands containing $\sqrt{ax^2 + bx + c}$, choose the Euler substitution according to the sign of $a$:

- $a > 0$: $\sqrt{a}\,x + t = \sqrt{ax^2 + bx + c}$
- $a < 0$: $t(x - \alpha) = \sqrt{ax^2 + bx + c}$ ($\alpha$ one of the roots)
- $c > 0$: $\sqrt{c} + xt = \sqrt{ax^2 + bx + c}$

These substitutions reduce integrals of rational functions with square roots to integrals of rational functions.

### The Special Function Table

When the Risch algorithm proves that an integral has **no elementary antiderivative**, many common cases still have closed forms in terms of special functions. oCAS encodes these standard antiderivatives directly in `ocas-calc/src/integral/special.rs` (definitions consistent with SymPy):

| Integrand | Antiderivative | Special function |
|---|---|---|
| $e^{-x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erf}(x)$ | Error function |
| $e^{x^2}$ | $\frac{\sqrt{\pi}}{2}\,\text{erfi}(x)$ | Imaginary error function |
| $e^{cx^2}$ ($c < 0$) | $\frac{\sqrt{\pi}}{2\sqrt{-c}}\,\text{erf}(\sqrt{-c}\,x)$ | Error function |
| $e^x / x$ | $\text{Ei}(x)$ | Exponential integral |
| $\sin(x)/x$ | $\text{Si}(x)$ | Sine integral |
| $\cos(x)/x$ | $\text{Ci}(x)$ | Cosine integral |
| $\sinh(x)/x$ | $\text{Shi}(x)$ | Hyperbolic sine integral |
| $\cosh(x)/x$ | $\text{Chi}(x)$ | Hyperbolic cosine integral |
| $\sin(x^2)$ | $\sqrt{\frac{\pi}{2}}\,S\!\left(\sqrt{\frac{2}{\pi}}\,x\right)$ | Fresnel $S$ |
| $\cos(x^2)$ | $\sqrt{\frac{\pi}{2}}\,C\!\left(\sqrt{\frac{2}{\pi}}\,x\right)$ | Fresnel $C$ |

**Matching logic**: `special_integrate` decomposes the integrand into product factors and tries to match the following pattern families:

- **erf family**: $e^{c \cdot x^2}$ ($c < 0$), including the special case $e^{-x^2}$
- **Ei family**: $e^x / x$, including the generalization $e^{cx}/x$
- **Si/Ci/Shi/Chi family**: trigonometric or hyperbolic functions divided by $x$
- **Fresnel family**: $\sin(x^2)$, $\cos(x^2)$

Note: these are **not** part of the Risch algorithm. After Risch proves "no elementary antiderivative exists", the special-function table provides a "suboptimal" answer.

---

## Implementation in oCAS

### The Integration Pipeline Architecture

`integrate(expr, var)` (`ocas-calc/src/integral/mod.rs`) tries the following layers in order; the first layer that produces an answer wins:

```
┌────────────────────────────────────────────────────────┐
│  integrate(expr, var)                                  │
├────────────────────────────────────────────────────────┤
│ 1. Lookup table (power rule, linear-argument functions)│
│     ↓ failure                                          │
│  2. Rational-function integrator (integrate_rational)  │
│     ↓ failure                                          │
│  3. Risch algorithm (risch_integrate)                  │
│     ↓ failure                                          │
│  4. Trig rewrite (trig_to_exp) → retry Risch → realify   │
│     ↓ failure                                          │
│  5. Special-function table (special_integrate)         │
│     ↓ failure                                          │
│  6. Heuristic layer (heuristic_integrate)              │
│     ↓ failure                                          │
│  7. Return the unevaluated form Integral(expr, var)    │
└────────────────────────────────────────────────────────┘
```

**Depth limits**: `MAX_DEPTH = 8` (the recursion depth of `integrate_raw`) and `MAX_RISCH_DEPTH = 16` (the recursion depth of `risch_integrate`, using a `thread_local!` counter). These limits prevent stack overflow caused by pathological inputs (e.g. `sec(x)` triggering infinite VOP retries).

### Module Organization

```
ocas-calc/src/integral/
├── mod.rs         ← integrate entry + lookup tables + pipeline dispatch
├── rational.rs    ← rational-function integrator over the base field Q(x)
│                    Hermite reduction + logarithmic part + Rothstein–Trager
├── risch.rs       ← main loop of the Risch algorithm
│                    per-level integrate_level + Hermite reduction
│                    polynomial part (primitive / hyperexponential)
├── rde.rs         ← Risch differential equation solver
│                    base_rde (Q[x]) + recursive tower
├── trig.rs        ← trigonometric ↔ complex-exponential rewrite (trig_to_exp / realify)
├── heuristic.rs   ← heuristics: parts / trig substitution / Weierstrass / Euler
└── special.rs     ← special-function table (erf, Ei, Si, Ci, Fresnel, etc.)
```

### Execution Flow of the Risch Algorithm

Trace the full execution with $\int (x+1) e^x\,dx$ as an example:

1. **Build the tower**: `build_tower` collects $t_1 = \exp(x)$; the tower is $\mathbb{Q}(x) \subset \mathbb{Q}(x)(t_1)$.
2. **`risch_integrate`**: convert the integrand into a rational function in $k_1(t_1)$: $(x+1) \cdot t_1$.
3. **`integrate_level(tower, level=1, f)`**:
   - The denominator is 1; no Hermite reduction.
   - Polynomial part $p(t_1) = (x+1) t_1$.
   - Call `integrate_kpoly_hyperexp` (since $t_1 = \exp(x)$):
     - Compare the $t_1^1$ coefficient: $Da_1 + a_1 \cdot 1 = x + 1$
     - Recursively call `rde_solve`: $a_1' + a_1 = x + 1$
     - Base-field `base_rde`: write $a_1 = ax + b$; then $a_1' = a$ and $a + ax + b = x + 1$
     - Compare: $a = 1$, $a + b = 1 \Rightarrow b = 0$
     - $a_1 = x$; verification: $x' + x = 1 + x = x + 1$ ✓
   - The $t_1^0$ coefficient is 0, so $b_0 = 0$.
4. **Assemble the result**: $q = x \cdot t_1$, converted to an atom: $x \cdot e^x$.
5. **Simplify + normalize**: return $x \cdot e^x$.

**Verification**: $\frac{d}{dx}[x \cdot e^x] = e^x + x \cdot e^x = (x+1)e^x$ ✓

### Internal Representation of Field Elements

In the tower, field elements are represented by **flat multivariate sparse polynomials**:

```rust
// KElem: an element of k = Q(x, t_1, ..., t_{n-1})
struct KElem {
    num: SparseMultivariatePolynomial<RationalDomain, Lex>,
    den: SparseMultivariatePolynomial<RationalDomain, Lex>,
}
```

Variable index assignment: `0` corresponds to $x$, `1` to $t_1$, …, `n-1` to $t_{n-1}$.

`KPoly` is a dense univariate polynomial in the top-level variable with coefficients in `KElem`:

```rust
// KPoly: a polynomial in k[t_n] with coefficients in k
struct KPoly {
    coeffs: Vec<KElem>,  // from lowest to highest degree
    top: usize,          // index of the top-level variable
    n: usize,            // total number of variables
}
```

`KRat` is a rational function $p/q$:

```rust
struct KRat {
    num: KPoly,
    den: KPoly,
}
```

**No GCD reduction**: `KElem` does not perform multivariate GCD reduction (because a general multivariate GCD is expensive). Zero detection only checks whether the numerator is zero; cross-multiplied equality is the only reliable equality test.

### Fuel-Bounded Integration

`integrate_with_fuel(ctx, expr, var, &fuel)` threads a `Fuel` budget through the two post-integration simplification phases. The integration traversal itself uses the `MAX_DEPTH` / `MAX_RISCH_DEPTH` limits; `Fuel` constrains only the simplification phases, preventing pathological results (e.g. integrals producing enormous expressions) from sending the rewriter into an infinite loop.

```rust
use ocas_core::fuel::Fuel;

let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, expr, Symbol::new("x"), &fuel);
// Ok(result)  — completed normally
// Err(_)      — fuel exhausted during simplification
```

### Limitations and Fallback Strategy

The current implementation of the Risch algorithm has the following limitations:

| Limitation | Reason | Fallback behavior |
|---|---|---|
| Only polynomial solutions of the RDE | Rational solutions need denominator-bound analysis | Return `None`; the caller tries other layers |
| At tower levels the logarithmic part only uses the logarithmic-derivative identity $a_1 = c \cdot Dd_1$ | The full logarithmic part needs tower-level Rothstein–Trager / trace-function techniques | Return the unevaluated form (the base field $\mathbb{Q}(x)$ still gets the full treatment) |
| Algebraic functions ($\sqrt{x}$, etc.) not supported | Needs algebraic-function field extensions | Common patterns are covered by the heuristic trigonometric substitution |
| Conservative rejection of algebraically dependent generators | Detecting relations like $\log(2x)$ vs. $\log(x)$ | Return `None` |
| Hyperexponential RDEs containing $I$ | The RDE solver works only over $\mathbb{Q}[x]$ | Trigonometric integrands are returned in unevaluated form |

When all layers fail, `Integral(expr, var)` is returned — this is an **intentional answer**, meaning "this integral has no closed form in the current implementation", not a program error.

---

## References

1. **Bronstein, M.** *Symbolic Integration I: Transcendental Functions*, 2nd ed., Springer, 2005. — The authoritative reference on the Risch algorithm. Chapter 2: integration of rational functions (Hermite reduction, Rothstein–Trager); Chapter 5: the Risch algorithm for elementary functions; Chapter 6: the Risch differential equation. The oCAS implementation follows this book closely.
2. **Bronstein, M.** *Symbolic Integration II: Transcendental and Algebraic Functions*, Springer, 2004. — Extensions to algebraic and special functions. Not yet covered by oCAS.
3. **Geddes, K. O., Czapor, S. R., & Labahn, G.** *Algorithms for Computer Algebra*, Kluwer, 1992. — Chapter 12 offers an alternative presentation of the Risch algorithm.
4. **Liouville, J.** "Sur les transcendantes elliptiques de première et de seconde espèce considérées comme fonctions de leur amplitude." *Journal de l'École Polytechnique*, 1835. — The original source of Liouville's theorem on elementary antiderivatives.
5. **Risch, R. H.** "The problem of integration in finite terms." *Transactions of the AMS*, 139:167–189, 1969. — The founding paper of the Risch algorithm.
6. **Rothstein, M.** "A new algorithm for the integration of exponential and logarithmic functions." *Proceedings of the 1977 MACSYMA Users Conference*, 1977. — The Rothstein–Trager resultant method.
7. **Trager, B. M.** "Algebraic factoring and rational function integration." *Proceedings of SYMSAC '76*, 1976. — Trager's algorithm (factorization over algebraic number fields).
8. **Lazard, D. & Rioboo, R.** "Integration of rational functions: Rational computation of the logarithmic part." *Journal of Symbolic Computation*, 9(2):113–129, 1990. — An efficient algorithm for the logarithmic part.
