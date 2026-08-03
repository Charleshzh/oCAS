# FGLM & Elimination Theory

This chapter systematically presents the FGLM algorithm for changing the order of Gröbner bases of zero-dimensional ideals, the elimination theorem, and Gröbner-basis-based ideal operations (quotients, saturations, intersections, primary decomposition, radicals). Together these tools form the core infrastructure of oCAS's polynomial-system solving and algebraic-geometry computations.

---

## Prerequisites

Before reading this chapter, we recommend the following background:

- **Gröbner basis theory**: Buchberger's algorithm, the F4/F5 algorithms, S-polynomials, reduced Gröbner bases — see [Gröbner Basis Theory](./groebner-theory.md)
- **Polynomial algebra**: multivariate polynomial rings $\mathbb{F}[x_1, \ldots, x_n]$, monomial orders (Lex / Grlex / Grevlex) — see [Polynomial Algebra](./polynomial-algebra.md)
- **Linear algebra**: Gaussian elimination, vector-space dimension, linear dependence — see [Linear Algebra](./linear-algebra.md)
- **Polynomial factorization**: square-free factorization, univariate polynomial factorization — see [Polynomial GCD & Factorization](./poly-gcd-factoring.md)

For a systematic study path, see the [Mathematics Overview](./overview.md).

---

## Basic Concepts

### Zero-Dimensional Ideals

**Definition**. Let $I \subseteq R = \mathbb{F}[x_1, \ldots, x_n]$ be an ideal. $I$ is called **zero-dimensional** if the quotient ring $R/I$ is **finite-dimensional** as an $\mathbb{F}$-vector space.

Equivalent conditions (any of the following implies zero-dimensionality):

1. **Algebraic geometry**: the affine variety $V(I) \subseteq \overline{\mathbb{F}}^n$ is a finite set
2. **Gröbner-basis characterization**: for every variable $x_i$ ($1 \leq i \leq n$), the Gröbner basis contains a polynomial $g$ with $\text{lm}(g) = x_i^{N_i}$ (a pure power)
3. **Finiteness of the staircase**: the set of standard monomials (monomials not divisible by any leading monomial) is finite

**Example**. $I = \langle x^2 - 1, y - x \rangle \subseteq \mathbb{Q}[x, y]$. Its Lex Gröbner basis is $\{x^2 - 1, y - x\}$ with leading terms $\{x^2, y\}$. For $x$ we have $x^2$ and for $y$ we have $y$ — the pure-power condition is satisfied, so $I$ is zero-dimensional. Indeed $V(I) = \{(1, 1), (-1, -1)\}$ is finite.

**Counterexample**. $J = \langle x - y \rangle \subseteq \mathbb{Q}[x, y]$. Its leading terms are $\{x\}$; there is no pure-power leading term for $y$. $V(J) = \{(t, t) : t \in \overline{\mathbb{Q}}\}$ is infinite, so $J$ is positive-dimensional.

### The Staircase and Standard Monomials

**Definition**. Let $G = \{g_1, \ldots, g_t\}$ be a Gröbner basis of an ideal $I$ with respect to a monomial order $\succ$. Let $\text{LT}(I) = \langle \text{lm}(g_1), \ldots, \text{lm}(g_t) \rangle$ be the leading-term ideal.

The set of **standard monomials** is defined as

$$\text{Std}(I) = \{x^\alpha \in \mathbb{F}[x_1, \ldots, x_n] : x^\alpha \notin \text{LT}(I)\}$$

i.e. the monomials not divisible by any $\text{lm}(g_i)$. The set of standard monomials is also called the **staircase**, because in two dimensions its shape resembles a staircase.

**Core property**. The standard monomials form an $\mathbb{F}$-vector-space basis of $R/I$. That is, every coset $f + I \in R/I$ has a **unique** representation

$$f + I = \sum_{x^\alpha \in \text{Std}(I)} c_\alpha \cdot x^\alpha + I$$

with $c_\alpha \in \mathbb{F}$. This is equivalent to saying that the normal form of $f$ with respect to the Gröbner basis is a unique linear combination of standard monomials.

**Example**. $I = \langle x^2, xy \rangle \subseteq \mathbb{F}[x, y]$ (Grevlex order). The leading terms are $\{x^2, xy\}$.

- $1$: not divisible by $x^2$ or $xy$ → standard monomial ✓
- $y$: not divisible → ✓
- $y^2$: not divisible → ✓
- $y^k$ ($k \geq 0$): not divisible → ✓
- $x$: not divisible → ✓
- $x^2$: divisible by $x^2$ → ✗
- $xy$: divisible by $xy$ → ✗

The staircase is $\{1, x, y, y^2, y^3, \ldots\}$ — an infinite set! This reflects that $I$ is positive-dimensional ($V(I) = \{(0, t) : t \in \overline{\mathbb{F}}\}$).

**Dimension**. For a zero-dimensional ideal, the size of the staircase

$$D = \dim_{\mathbb{F}}(R/I) = |\text{Std}(I)|$$

is a finite positive integer, called the **vector-space dimension** of $I$ or the **stable value of the Hilbert function**. This dimension $D$ plays a central role in the complexity analysis of the FGLM algorithm.

### Normal Forms

**Definition**. The **normal form** $\text{NF}(f, G)$ of a polynomial $f$ with respect to a Gröbner basis $G$ is the remainder of dividing $f$ by $G$ with remainder. It satisfies:

1. $\text{NF}(f, G) \in f + I$ ($\text{NF}(f, G)$ is congruent to $f$ modulo $I$)
2. Every term of $\text{NF}(f, G)$ is a standard monomial
3. The representation is unique (for a reduced Gröbner basis)

For a zero-dimensional ideal, the normal form can be represented as a coefficient vector of length $D$ — one component per standard monomial. This is precisely the key that lets the FGLM algorithm turn normal-form computation into linear-algebra operations.

---

## Core Theory

### The Elimination Theorem

**Elimination ideal**. Let $I \subseteq \mathbb{F}[x_1, \ldots, x_n]$ and fix $1 \leq \ell \leq n$. The $\ell$-th **elimination ideal** is

$$I_\ell = I \cap \mathbb{F}[x_{\ell+1}, \ldots, x_n]$$

i.e. the ideal of the polynomials of $I$ that do not involve the first $\ell$ variables.

**Elimination Theorem**. Let $G$ be a Gröbner basis of $I$ with respect to the **Lex order** ($x_1 \succ x_2 \succ \cdots \succ x_n$). Then

$$G_\ell = G \cap \mathbb{F}[x_{\ell+1}, \ldots, x_n]$$

is a Gröbner basis of $I_\ell$.

In other words, the polynomials of a Lex Gröbner basis that do not involve $x_1, \ldots, x_\ell$ generate exactly the elimination ideal $I_\ell$.

**Sketch of proof**. The key property of the Lex order is: if $x^\alpha \succ x^\beta$ and $\alpha_1 = \cdots = \alpha_\ell = 0$, then $\beta_1 = \cdots = \beta_\ell = 0$. Hence the intersection of the leading-term ideal $\langle \text{lt}(g) : g \in G \rangle$ with $\mathbb{F}[x_{\ell+1}, \ldots, x_n]$ is generated exactly by the leading terms of $G_\ell$. $\square$

**Applications**. The elimination theorem is the theoretical basis of all of the following constructions:

- **Implicit-function elimination**: eliminate variables from a system of equations, obtaining relations in the remaining variables only
- **The Rabinowitsch trick for ideal quotients**: eliminate after introducing an auxiliary variable
- **Ideal intersection**: eliminate after introducing the auxiliary variable $t$
- **Polynomial-system solving**: extract constraints variable by variable from the Lex GB

**Note**. Gröbner bases in the Grevlex order do **not** have the elimination property. Hence every operation requiring elimination must be performed in the Lex order. This is exactly the core value of the FGLM order-change algorithm: first compute the basis quickly in Grevlex, then convert it to Lex with FGLM for elimination.

### The FGLM Algorithm

The FGLM algorithm (Faugère–Gianni–Lazard–Mora, 1993) provides an efficient change of order for Gröbner bases of **zero-dimensional ideals**. Its core idea: using normal-form computations in the source order, enumerate monomials in the target order and construct the target-order Gröbner basis directly by linear-algebra methods — completely avoiding a re-run of F4 in the target order.

#### Motivation

Suppose we already have a Gröbner basis $G_{\text{grevlex}}$ of $I$ in the Grevlex order (the most efficient computational order). For elimination or solving, we need a basis $G_{\text{lex}}$ in the Lex order.

Comparison of the two methods:

| Method | Complexity | Description |
|------|--------|------|
| Re-running F4 (Lex) | exponential (worst case) | recomputes everything; intermediate polynomials may blow up |
| FGLM | $O(n \cdot D^3)$ | $D = \dim(R/I)$; pure linear algebra |

For zero-dimensional ideals, FGLM is almost always faster — especially when $D$ is moderate but the intermediate polynomials of the Lex basis are enormous.

#### Algorithm Description

**Input**: a reduced Gröbner basis $G_{\text{src}}$ of a zero-dimensional ideal $I$ in the source order $O_1$, and a target order $O_2$.

**Output**: a reduced Gröbner basis $G_{\text{tgt}}$ of $I$ in the target order $O_2$.

**Initialization**:

1. Extract all leading monomials $\text{lm}(g_1), \ldots, \text{lm}(g_t)$ from $G_{\text{src}}$
2. Compute the staircase $\text{Std}(I)$ (BFS: start from $1$, extend variable by variable, skip monomials divisible by a leading term, until no new standard monomial appears — zero-dimensionality guarantees termination)
3. Let $D = |\text{Std}(I)|$ and number the standard monomials $s_1, s_2, \ldots, s_D$
4. For each variable $x_i$, construct the **multiplication matrix** $M_i$: for each standard monomial $s_j$, compute the normal form of $x_i \cdot s_j$ (under $G_{\text{src}}$), giving a $D \times D$ matrix

**Main loop**: enumerate the monomials in increasing target order $O_2$:

```
seen_nfs ← []     // seen normal forms (D-dimensional vectors)
seen_mons ← []    // corresponding monomials
boundary ← {1}    // boundary set
visited ← ∅       // visited monomials
GB_tgt ← ∅

while boundary ≠ ∅:
    m ← the O2-smallest monomial in boundary
    visited ← visited ∪ {m}

    // compute the normal form of m w.r.t. the source-order basis (D-dim coordinate vector)
    nf ← normal_form(m, G_src, staircase)

    if nf ∈ span(seen_nfs):      // linearly dependent?
        // found a relation nf = Σ c_i · seen_nfs[i]
        coeffs ← Gaussian_solve(seen_nfs, nf)
        // construct the new polynomial: m - Σ c_i · seen_mons[i]
        new_poly ← (m, 1) + Σ (-c_i, seen_mons[i])
        GB_tgt ← GB_tgt ∪ {new_poly}
        // all multiples of m lie in the ideal; mark them as visited
        mark_multiples(visited, m)
    else:
        // linearly independent: m is a new standard monomial in the target order
        seen_nfs.append(nf)
        seen_mons.append(m)
        // add the neighbours of m to the boundary
        for i = 1 to n:
            if m · x_i ∉ visited:
                boundary ← boundary ∪ {m · x_i}

return minimize(GB_tgt).auto_reduce()
```

#### Key Subroutines

**Normal-form computation**. The normal form of a monomial $m$ under $G_{\text{src}}$ is a $D$-dimensional vector. In the implementation, $m$ is turned into the monomial polynomial with coefficient $1$, divided by $G_{\text{src}}$ with remainder, and then each term of the remainder is mapped to its position index in the staircase.

**Linear-relation detection**. Maintain an augmented matrix $[\text{seen\_nfs} \mid \text{nf}]$ ($D$ rows, $k+1$ columns, where $k$ is the number of seen normal forms). Perform Gaussian elimination whenever a new nf is added:

- if the augmented column is all zero after elimination → linearly dependent; extract the coefficient vector
- if the augmented column has a non-zero entry after elimination → linearly independent; add to seen_nfs

This can be implemented efficiently with incremental Gaussian elimination (processing one column at a time), avoiding repeated elimination.

**Marking multiples**. When a linear relation $\text{nf}(m) = \sum c_i \cdot \text{nf}(m_i)$ is found, the constructed polynomial $m - \sum c_i m_i \in I$ has leading term $m$ (in the target order). All multiples $m \cdot x^\alpha$ of $m$ necessarily lie in $I$ as well, so they need not be considered further — they are removed from the boundary.

#### Complexity Analysis

- **Staircase computation**: $O(D)$ monomials, each requiring one divisibility check against every leading monomial of the basis (each check costs $O(n)$ comparisons)
- **Multiplication matrices**: $n$ matrices of size $D \times D$, one normal-form computation per entry → $O(n \cdot D^2 \cdot T_{\text{NF}})$
- **Main loop**: at most $D$ linearly independent standard monomials are discovered, and each linear-dependence discovery emits one target-basis element; the total number of iterations is $\le D + |G_{\text{tgt}}|$, with the implementation's safety cap at $D(D+1) + 4n$
- **Per iteration**: one normal-form computation $O(T_{\text{NF}})$ + one Gaussian elimination $O(D^2)$
- **Total complexity**: $O(n \cdot D^3)$ field operations

where $T_{\text{NF}}$ is the cost of one normal-form computation, typically linear in $D$ for sparse polynomials.

#### Correctness

**Theorem** (Faugère–Gianni–Lazard–Mora 1993). Let $I$ be a zero-dimensional ideal and $G_{\text{src}}$ a reduced Gröbner basis of $I$ in the order $O_1$. The FGLM algorithm outputs a $G_{\text{tgt}}$ that is a reduced Gröbner basis of $I$ in the order $O_2$.

**Sketch of proof**:

1. Normal-form computation is order-independent — the normal-form map w.r.t. $G_{\text{src}}$ and the one w.r.t. $G_{\text{tgt}}$ are representations of the same linear map in different bases
2. The main loop enumerates monomials in increasing target order and exactly covers all standard monomials and the monomials producing leading terms
3. The linear-relation detection guarantees that the output polynomials generate exactly the leading-term ideal
4. `minimize` and `auto_reduce` guarantee reducedness

### Ideal Operations

Gröbner bases turn the basic arithmetic operations on ideals into algorithms. This section presents the core ideal operations implemented in oCAS; each one relies on the elimination theorem.

#### Ideal Quotient

**Definition**. The **ideal quotient** of two ideals $I, J \subseteq R$ is

$$I : J = \{f \in R : f \cdot g \in I, \, \forall g \in J\}$$

For a principal ideal $J = \langle g \rangle$ it is written simply as $I : g$.

**Rabinowitsch trick**. The key trick for computing $I : g$ is to introduce a new variable $w$. Note, however, that the auxiliary-variable construction actually computes the **saturation**:

$$I : g^\infty = \left( I \cdot R[w] + \langle 1 - wg \rangle \right) \cap R$$

i.e. in the extended ring $R[w]$, lift the generators of $I$, add $1 - wg$, compute a Gröbner basis, and then take the polynomials not involving $w$ (elimination). When $I$ is already saturated with respect to $g$ (i.e. $I : g = I : g^\infty$, e.g. square-free ideals without embedded components), this equals $I : g$. oCAS's `ideal_quotient` implements exactly this elimination construction (mathematically, $I : g^\infty$).

**Proof of correctness**. By localization, $R[w]/\langle 1 - wg \rangle \cong R[1/g]$ ($w \mapsto 1/g$). Hence $f \in (I \cdot R[w] + \langle 1-wg\rangle) \cap R$ iff $f \in I \cdot R[1/g]$, i.e. iff $g^N f \in I$ for some $N$, i.e. $f \in I : g^N \subseteq I : g^\infty$. Conversely, if $g^N f \in I$, then

$$f = f\bigl(1-(wg)^N\bigr) + w^N \cdot (g^N f)$$

where the first term lies in $\langle 1-wg\rangle$ (since $1-(wg)^N$ is divisible by $1-wg$) and the second in $I \cdot R[w]$. $\square$

**Several generators**. For $J = \langle g_1, \ldots, g_m \rangle$, one has

$$I : J^\infty = \bigcap_{j=1}^{m} (I : g_j^\infty)$$

Implemented by computing each $I : g_j^\infty$ and then intersecting the results (which is precisely the behavior of `ideal_quotient`).

**Example**. $I = \langle x^2, xy \rangle$, $J = \langle x \rangle$.

$$I : x = \langle x^2, xy \rangle : \langle x \rangle$$

In $\mathbb{Q}[x, y, w]$, compute $\text{GB}(x^2, xy, 1 - wx)$ and eliminate $w$:

- $1 - wx \implies w = 1/x$
- $x^2 \cdot w = x$, $xy \cdot w = y$ (since $x^2 w - x = -x(1-wx)$ and $xy w - y = -y(1-wx)$)

Elimination yields $\langle x, y \rangle$, hence $\langle x^2, xy \rangle : \langle x \rangle = \langle x, y \rangle$.

#### Ideal Intersection

**Definition**. $I \cap J$ is the ideal of the polynomials that belong to both $I$ and $J$.

**The auxiliary-variable method**. Introduce a new variable $t$:

$$I \cap J = \langle t \cdot f_i, \, (1-t) \cdot g_j \rangle_{i,j} \cap R$$

where $\{f_i\}$ and $\{g_j\}$ are generators of $I$ and $J$ respectively.

**Intuition**. At $t = 1$ the constraints degenerate to $f_i = 0$ (the conditions of $I$); at $t = 0$ they degenerate to $g_j = 0$ (the conditions of $J$). After eliminating $t$, the polynomials satisfying both sets of conditions are obtained.

**Correctness**. $h \in I \cap J$ ⟺ $h = t \cdot h + (1-t) \cdot h$, where $t \cdot h \in t \cdot I$ and $(1-t) \cdot h \in (1-t) \cdot J$. Conversely, a polynomial in the elimination ideal vanishes at both $t = 0$ and $t = 1$, so it belongs to both $I$ and $J$. $\square$

**Example**. $\langle x \rangle \cap \langle y \rangle = \langle xy \rangle$.

In $\mathbb{Q}[x, y, t]$, compute $\text{GB}(tx, (1-t)y)$ and eliminate $t$:

- $tx = 0 \implies t = 0$ or $x = 0$
- $(1-t)y = 0 \implies t = 1$ or $y = 0$

Elimination yields $\langle xy \rangle$. $\square$

#### Ideal Saturation

**Definition**. The **saturation** of $I$ with respect to $J$ is

$$I : J^\infty = \bigcup_{k=1}^{\infty} (I : J^k)$$

Equivalently, $I : J^\infty$ is the "restriction" of $I$ away from $V(J)$ — it agrees with $I$ outside the zero set of $J$.

**Iterative computation**. By the ascending chain condition there exists $k_0$ with $I : J^{k_0} = I : J^{k_0+1} = \cdots = I : J^\infty$. The algorithm repeatedly computes $I : J$ until stabilization:

```
I_old ← I
loop:
    I_new ← I_old : J
    if I_new == I_old:  // compare Gröbner bases
        return I_new
    I_old ← I_new
```

**Applications**. Saturation is extremely important in algebraic geometry:

- **Radical computation**: $\sqrt{I} = I : h^\infty$ (in the positive-dimensional case, $h$ is related to the Jacobian)
- **Primary decomposition**: separating the components of different dimensions
- **Closure computation**: projective closures, complete intersections

**Example**. $\langle x^2 y, xy^2 \rangle : \langle x \rangle^\infty$.

Computing the quotients round by round according to the mathematical definition:

First round: $\langle x^2 y, xy^2 \rangle : \langle x \rangle = \langle xy, y^2 \rangle$

Second round: $\langle xy, y^2 \rangle : \langle x \rangle = \langle y \rangle$

Third round: $\langle y \rangle : \langle x \rangle = \langle y \rangle$ (stable)

Hence $\langle x^2 y, xy^2 \rangle : \langle x \rangle^\infty = \langle y \rangle$. $\square$

(Note: oCAS's `ideal_quotient` computes $I : g^\infty$ in a single elimination step — in this example one call already yields $\langle y \rangle$, and the iteration in `ideal_saturate` merely confirms stability.)

#### Sums, Products, and Membership

These are the most basic ideal operations; they are comparatively simple to implement:

| Operation | Definition | Implementation |
|------|------|------|
| $I + J$ | $\langle f_1, \ldots, f_m, g_1, \ldots, g_k \rangle$ | merge the generators and compute a GB |
| $I \cdot J$ | $\langle f_i g_j \rangle_{i,j}$ | all pairwise products of the generators, then compute a GB |
| $f \in I$ ? | — | compute $\text{NF}(f, G_I)$ and check whether it is zero |

### Primary Decomposition

**Definition**. An ideal $Q \subseteq R$ is called a **primary ideal** if $fg \in Q$ implies $f \in Q$ or $g^n \in Q$ (for some $n \geq 1$).

**Primary decomposition theorem** (Lasker–Noether). Every ideal $I \subseteq \mathbb{F}[x_1, \ldots, x_n]$ ($\mathbb{F}$ of characteristic zero) has a primary decomposition

$$I = Q_1 \cap Q_2 \cap \cdots \cap Q_r$$

where each $Q_i$ is a primary ideal. The corresponding **associated prime ideals** are $\mathfrak{p}_i = \sqrt{Q_i}$.

#### Zero-Dimensional Primary Decomposition

For zero-dimensional ideals, primary decomposition can be achieved by factoring the univariate polynomials in a Lex Gröbner basis.

**Key observation**. The Lex GB of a zero-dimensional ideal contains a polynomial $p_1(x_1)$ in $x_1$ alone, a polynomial $p_2(x_1, x_2)$ in $x_1, x_2$, and so on. The factorization of $p_1(x_1)$ corresponds to the different branches of the ideal in the $x_1$ direction.

**Algorithm** (matching the oCAS implementation `primary_decomp_zero_dim`):

1. Compute the Lex Gröbner basis $G$ of $I$
2. Extract the polynomial $p_1(x_1)$ in $x_1$ alone and compute its **square-free part** $\tilde{p}_1 = p_1 / \gcd(p_1, p_1')$
3. Factor $\tilde{p}_1$ over $\mathbb{Q}$ into pairwise distinct irreducible factors $\tilde{p}_1 = q_1 \cdot q_2 \cdots q_s$
4. For each factor $q_i$, separate the branch by saturation: $I_i = I : \left(\prod_{j \neq i} q_j\right)^\infty$ (the source saturates by each other factor sequentially)
5. The associated prime ideal of branch $I_i$ is taken as $\mathfrak{p}_i = \mathrm{GB}(I + \langle q_i \rangle)$

**Note**. If $\tilde{p}_1$ has only one irreducible factor (or $G$ contains no polynomial in $x_1$ alone), the algorithm returns the single branch $\{\text{primary} = I,\; \text{prime} = \sqrt{I}\}$. The current implementation handles only the factorization in the first variable $x_1$; it does **not** recurse into the remaining variables.

**Example**. $I = \langle x^2 - 1, y - x \rangle \subseteq \mathbb{Q}[x, y]$.

Lex GB: $\{x^2 - 1, y - x\}$. $p_1(x) = x^2 - 1$ is already square-free and factors as $q_1 = x - 1$, $q_2 = x + 1$.

- Branch 1: $I_1 = I : \langle x + 1 \rangle^\infty = \langle x - 1, y - 1 \rangle$ (associated prime $\langle x - 1, y - 1 \rangle$)
- Branch 2: $I_2 = I : \langle x - 1 \rangle^\infty = \langle x + 1, y + 1 \rangle$ (associated prime $\langle x + 1, y + 1 \rangle$)

Check: $\langle x^2 - 1, y - x \rangle = \langle x - 1, y - 1 \rangle \cap \langle x + 1, y + 1 \rangle$; both branches are prime.

For comparison, for $I = \langle x^2, xy \rangle$ the square-free part of $p_1(x) = x^2$ has only the single factor $x$, so the algorithm returns the single branch $\{\text{primary} = \langle x^2, xy \rangle,\; \text{prime} = \sqrt{I} = \langle x \rangle\}$. The ideal nevertheless admits the finer decomposition $\langle x^2, xy \rangle = \langle x \rangle \cap \langle x^2, y \rangle$, where $\langle x \rangle$ is prime and $\langle x^2, y \rangle$ is primary with associated prime $\langle x, y \rangle$ (an embedded component).

#### Positive-Dimensional Primary Decomposition

Primary decomposition of positive-dimensional ideals is much more complicated; oCAS currently marks it as not yet implemented (`TODO`). A complete implementation requires:

- the Gianni–Trager–Zacharias algorithm
- or the characteristic method of Eisenbud–Huneke–Vasconcelos

### Radicals

**Definition**. The **radical** of an ideal $I$ is

$$\sqrt{I} = \{f \in R : f^n \in I \text{ for some } n \geq 1\}$$

Equivalently, $\sqrt{I}$ is the largest ideal containing $I$ with $V(\sqrt{I}) = V(I)$ (the same zero set).

#### Zero-Dimensional Radicals

For zero-dimensional ideals, radical computation uses the square-free factorization of the univariate polynomials in the Lex GB.

**Algorithm**:

1. Compute the Lex GB $G$
2. For each polynomial $p_i(x_i)$ in $G$ involving only $x_i$, compute its square-free part $\tilde{p}_i = p_i / \gcd(p_i, p_i')$
3. Replace $p_i$ by $\tilde{p}_i$ to obtain generators of $\sqrt{I}$

**Correctness**. The radical of a zero-dimensional ideal corresponds to reducing the exponent of each associated prime ideal to $1$. Square-free factorization achieves exactly this — the square-free part of $p_i^{e_i}$ is $p_i$.

**Example**. $\sqrt{\langle x^2, xy \rangle}$.

Lex GB: $\{x^2, xy\}$. The square-free part of $x^2$: $x^2 / \gcd(x^2, 2x) = x^2 / x = x$.

Radical: $\langle x, xy \rangle = \langle x \rangle$ (since $x$ already divides $xy$). $\square$

#### Positive-Dimensional Radicals: Jacobian Saturation

For positive-dimensional ideals, oCAS uses a Jacobian-based saturation method (a simplified Kemper algorithm):

$$\sqrt{I} = I : h^\infty$$

where $h$ is related to the Jacobian determinant.

**Algorithm** (characteristic zero, matching `radical_via_jacobian`):

1. For each generator $f_i$ of $I$ and each variable $x_j$, compute the partial derivative $\partial f_i / \partial x_j$ (skipping constant/zero derivatives)
2. Choose $h$ heuristically: the source does **not** compute the true gcd, but folds the non-trivial derivatives with `reduce`, keeping the one of **smallest total degree**. This is a heuristic approximation of the classical Jacobian-saturation formula ($\sqrt{I} = I : h^\infty$, where $h$ must be a suitable Jacobian-related element) — $h$ is just a single partial derivative (one entry of the Jacobian matrix); it is not guaranteed to divide the Jacobian determinant, nor is the saturation result guaranteed to equal $\sqrt{I}$
3. Compute $I : h^\infty$ (iterative saturation)

**Theoretical basis**. At a regular point, the rank of the Jacobian matrix equals the dimension of the variety (Jacobian criterion); the singular points are exactly those where the rank is deficient (less than $n - \dim$) — where the exponent structure is non-trivial. The saturation $I : h^\infty$ serves to remove these powers.

**Limitation**. When all partial derivatives are constant/zero, or $h$ is a constant polynomial (total degree 0), the algorithm falls back to returning the original GB. Since $h$ is chosen heuristically rather than being a genuine Jacobian-related element, the saturation result is only an approximation of $\sqrt{I}$ (it may be too large, or too small and miss components) — for positive-dimensional ideals, a complete radical computation requires more refined algorithms (such as Gianni–Trager–Zacharias or Eisenbud–Huneke–Vasconcelos).

### Testing Primality and Primaryness

**Zero-dimensional primality test**. A zero-dimensional ideal $I$ is prime if and only if every polynomial in its Lex GB involving only $x_i$ is irreducible.

oCAS implementation (`is_prime_zero_dim`):

1. Compute the Lex GB
2. Extract each univariate polynomial $p_i(x_i)$
3. For polynomials of degree $\le 3$, check reducibility with the **rational-root theorem** (a polynomial of degree $d \le 3$ over $\mathbb{Q}$ is irreducible iff it has no rational root); polynomials of degree $> 3$ receive no full factorization check

**Note**. Positive-dimensional ideals currently conservatively return `false` — a prime positive-dimensional ideal is misreported as non-prime (**false negative**). In the zero-dimensional case, however, a reducible degree-$>3$ univariate polynomial without rational roots (such as $(x^2-2)(x^2-3)$) is misjudged as irreducible, so **false positives** (reporting a non-prime ideal as prime) are possible.

**Primaryness test**. An ideal $I$ is primary if and only if it has exactly one associated prime ideal. The implementation checks the number of components via primary decomposition.

---

## Implementation in oCAS

### The FGLM Implementation

The FGLM algorithm is implemented in `ocas-poly/src/groebner/fglm.rs`.

#### Entry Function

```rust
pub fn fglm<D: Domain, O2: MonomialOrder>(
    gb: &GroebnerBasis<D, impl MonomialOrder>,
) -> Option<GroebnerBasis<D, O2>>
```

Generic parameters:

- `D`: the coefficient domain (implementing the `Domain` trait)
- `O2`: the target monomial order
- the source order of the input basis is determined by the type parameter `impl MonomialOrder` of `gb`

Returns `None` when the ideal is not zero-dimensional (the staircase computation finds infinitely many standard monomials).

#### Staircase Computation

```rust
fn compute_staircase(lms: &[Vec<usize>], n_vars: usize) -> Option<Vec<Vec<usize>>>
```

BFS algorithm:

- start from the all-zero exponent vector $1 = x_1^0 \cdots x_n^0$
- for each monomial $m$ in the queue, check whether it is divisible by some $\text{lm}(g_i)$
- if not divisible → add it to the staircase and enqueue its $n$ neighbours (incrementing one variable exponent of $m$ at a time)
- safety threshold 100,000: if the BFS visits more than this many monomials, the ideal is judged positive-dimensional

The divisibility check compares component by component:

```rust
fn monomial_divides_big(lm: &[usize], big: &[usize]) -> bool {
    lm.iter().zip(big.iter()).all(|(a, b)| a <= b)
}
```

#### Normal-Form Computation

```rust
fn normal_form_monomial<D: Domain>(
    m: &[usize],
    gb: &GroebnerBasis<D, impl MonomialOrder>,
    staircase: &[Vec<usize>],
    domain: &D,
) -> Vec<D::Element>
```

Build the monomial $m$ as the monomial polynomial $1 \cdot x^m$, divide it by the GB with remainder, and then map each term of the remainder to its position index in the staircase. Returns a coordinate vector of length $D$.

#### Linear-Relation Detection

```rust
fn find_relation<D: Domain>(
    seen: &[Vec<D::Element>],
    nf: &[D::Element],
    domain: &D,
) -> Option<Vec<D::Element>>
```

Build the augmented matrix $[\text{seen}^T \mid \text{nf}]$ ($D$ rows, $k+1$ columns) and perform forward elimination + back substitution:

1. find a pivot (a non-zero element) column by column and swap rows
2. normalize the pivot row and eliminate that column in the remaining rows
3. check consistency: if some row has all zeros on the left but a non-zero right-hand side → no solution (linearly independent)
4. extract the coefficient vector by back substitution

Returns `Some(coeffs)` when linearly dependent and `None` when linearly independent.

#### Marking Multiples

```rust
fn mark_multiples(
    visited: &mut HashMap<Vec<usize>, bool>,
    m: &[usize],
    n_vars: usize,
    max_deg: usize,
)
```

BFS marks all multiples of $m$ (with total degree at most `max_deg`) as visited. `max_deg` is set to $2D$ — large enough to cover all monomials that could produce new leading terms.

#### Main-Loop Details

The main loop maintains:

- `boundary`: the set of candidate monomials (initially $\{1\}$, the all-zero exponent vector)
- `seen_nfs`: the list of seen normal-form vectors
- `seen_mons`: the corresponding list of monomials
- `new_basis`: the output Gröbner basis

Each step takes the target-order-smallest monomial in the boundary, computes its normal form, and tests for linear dependence. If dependent → output the new polynomial and mark multiples; if independent → add it to `seen` and expand its neighbours.

Loop termination: the boundary is empty or the number of steps exceeds $D(D+1) + 4n$ (a safety bound).

Finally, `minimize()` (removing redundant polynomials) and `auto_reduce()` (mutual reduction) are applied to `new_basis`, producing the reduced Gröbner basis.

### Implementing Ideal Operations

The ideal operations are implemented in `ocas-poly/src/ideal.rs`, uniformly using the Lex order to guarantee the elimination property.

#### Membership

```rust
pub fn ideal_contains<D: Domain + 'static>(
    generators: &[SparseMultivariatePolynomial<D, Lex>],
    f: &SparseMultivariatePolynomial<D, Lex>,
    algo: Algorithm,
) -> bool
```

Compute the GB, then compute the normal form of $f$ and check whether it is zero.

#### Ideal Quotient

```rust
pub fn ideal_quotient<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

Internal flow:

1. For each generator $g_j$ of $J$, call `quotient_single_generator` (the Rabinowitsch trick)
2. `quotient_single_generator` implements:
   - lift the generators of $I$ into $R[w]$
   - add $1 - wg$
   - compute the Lex GB
   - eliminate $w$ (take the polynomials not involving $w$)
3. Intersect the results of all $I : g_j^\infty$ (`intersect_generators`)

#### Ideal Intersection

```rust
pub fn ideal_intersection<D: Domain + 'static>(
    generators_a: &[SparseMultivariatePolynomial<D, Lex>],
    generators_b: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

Internally calls `intersect_generators`:

1. Introduce the auxiliary variable $t$ (index 0)
2. Construct the generators $\{t \cdot f_i\} \cup \{(1-t) \cdot g_j\}$
3. Compute the Lex GB in the extended ring $R[t]$
4. Eliminate $t$ to obtain generators of $I \cap J$

#### Ideal Saturation

```rust
pub fn ideal_saturate<D: Domain + 'static>(
    generators_i: &[SparseMultivariatePolynomial<D, Lex>],
    generators_j: &[SparseMultivariatePolynomial<D, Lex>],
) -> GroebnerBasis<D, Lex>
```

Iteratively calls `ideal_quotient` until the GB stabilizes (a two-way containment check: every basis element of each side reduces to zero against the other), with a cap of 20 rounds; the current result is returned when the cap is exceeded.

#### Primary Decomposition

```rust
pub fn primary_decomposition(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> Vec<PrimaryComponent>
```

Currently only the zero-dimensional case is implemented (`primary_decomp_zero_dim`), and only the first variable is handled:

1. Compute the Lex GB
2. Extract the univariate polynomial $p_1(x_1)$ in $x_1$ alone and compute its square-free part
3. Factor the square-free part over $\mathbb{Q}$
4. Separate the branches by saturation (sequentially saturating by each other factor via `ideal_saturate`)
5. The associated prime of a branch is taken as $\mathrm{GB}(I + \langle q_i \rangle)$; a single factor yields the single branch ($\text{prime} = \sqrt{I}$)

The `PrimaryComponent` struct contains:

```rust
pub struct PrimaryComponent {
    pub primary: GroebnerBasis<RationalDomain, Lex>,  // the primary ideal
    pub prime: GroebnerBasis<RationalDomain, Lex>,     // the associated prime ideal
}
```

#### Radicals

```rust
pub fn ideal_radical(
    generators: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
) -> GroebnerBasis<RationalDomain, Lex>
```

Two cases:

- **Zero-dimensional** (`radical_zero_dim`): apply square-free factorization to each univariate polynomial in the Lex GB, replace them, and regenerate the GB
- **Positive-dimensional** (`radical_via_jacobian`): compute the non-trivial partial derivatives of all generators with respect to all variables, take $h$ as the one of **smallest total degree** (a heuristic, not the true gcd), and compute $I : h^\infty$

#### Other Tests

```rust
pub fn is_zero_dimensional(gb: &GroebnerBasis<RationalDomain, Lex>) -> bool
pub fn is_prime_ideal(generators: &[...]) -> bool
pub fn is_primary_ideal(generators: &[...]) -> bool
```

- `is_zero_dimensional`: checks whether every variable has a pure-power leading term
- `is_prime_ideal`: zero-dimensional: check irreducibility; positive-dimensional: conservatively return `false`
- `is_primary_ideal`: checks whether there is exactly one associated prime ideal

### Solving Zero-Dimensional Systems

```rust
pub fn solve_polynomial_system(
    equations: &[SparseMultivariatePolynomial<RationalDomain, Lex>],
    algo: Algorithm,
) -> PolynomialSystemSolution
```

Algorithm flow:

1. Compute the Lex GB (if the input is not already a GB)
2. Check zero-dimensionality
3. Solve by triangular decomposition (`solve_triangular`) — back-substitution starting from the **smallest variable** in the Lex order ($x_n$):
   - extract the polynomial in $x_n$ alone → isolate real roots and find them numerically
   - for each root, substitute into the polynomial involving $x_{n-1}$ → find the roots of $x_{n-1}$
   - recurse until all variables are determined (the results are returned in the order $x_1, \ldots, x_n$)

Return type:

```rust
pub enum PolynomialSystemSolution {
    ZeroDimensional(ZeroDimSolutions),  // finitely many solutions
    PositiveDimensional(GroebnerBasis),  // positive-dimensional (infinitely many solutions)
    Empty,                               // no solutions
}
```

---

## Advanced Topics

### Choosing Between FGLM and F4

| Scenario | Recommended method | Reason |
|------|----------|------|
| Zero-dimensional + Lex GB needed | Grevlex F4 → FGLM | $O(n \cdot D^3)$ is far faster than Lex F4 |
| Zero-dimensional + Grevlex GB needed | F4 directly | no order change needed |
| Positive-dimensional + Lex GB needed | Lex F4 or `reorder` | FGLM does not apply |
| Elimination + zero-dimensional | FGLM → elimination | the elimination theorem guarantees correctness |
| Elimination + positive-dimensional | Lex F4 | FGLM does not apply |

### FGLM over Arbitrary Fields

oCAS's FGLM implementation works generically over an arbitrary field $D$; normal-form computation goes entirely through the operations of the `Domain` trait (finite-field elements are stored as canonical representatives in $[0, p)$). Note that the `FpPoly` `i64` fast path of F4/F5 for $\mathbb{F}_p$ is not used by FGLM.

### Limitations

1. **Zero-dimensional only**: the staircase of a positive-dimensional ideal is infinite, so FGLM does not apply
2. **The $D^3$ bottleneck**: when the vector-space dimension $D$ is very large (e.g. $D > 10^4$), Gaussian elimination becomes the bottleneck
3. **Memory**: the $D \times D$ multiplication matrices and the augmented matrix must be stored

### Relation to the Literature

The oCAS implementation is based directly on the original papers:

- FGLM: Faugère, Gianni, Lazard, Mora (1993), *Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering*, JSC
- Ideal operations: Cox, Little, O'Shea, *Ideals, Varieties, and Algorithms*, Chapters 4, 8
- Primary decomposition: Gianni, Trager, Zacharias (1988), *Gröbner Bases and Primary Decomposition of Polynomial Ideals*
- Radical computation: Kemper (2002), *A Course in Commutative Algebra*

---

## References

1. J.-C. Faugère, P. Gianni, D. Lazard, T. Mora. "Efficient Computation of Zero-dimensional Gröbner Bases by Change of Ordering." *Journal of Symbolic Computation*, 16(4):329–344, 1993.
2. D. Cox, J. Little, D. O'Shea. *Ideals, Varieties, and Algorithms*. 4th ed., Springer, 2015. Chapters 3 (elimination), 4 (ideal quotients and saturation), 8 (primary decomposition).
3. W. W. Adams, P. Loustaunau. *An Introduction to Gröbner Bases*. AMS, 1994.
4. P. Gianni, B. Trager, G. Zacharias. "Gröbner Bases and Primary Decomposition of Polynomial Ideals." *Journal of Symbolic Computation*, 6(2–3):149–167, 1988.
5. G. Kemper. *A Course in Commutative Algebra*. Springer, 2011.
6. T. Becker, V. Weispfenning. *Gröbner Bases: A Computational Approach to Commutative Algebra*. Springer, 1993.

---

**See also**:

- [Rust API: Gröbner Bases and Ideals](../api/rust-groebner.md) — function signatures and complete examples
- [Gröbner Basis Theory](./groebner-theory.md) — the theoretical basis of the Buchberger/F4/F5 algorithms
- [Polynomial Algebra](./polynomial-algebra.md) — monomial orders and multivariate division with remainder
- [Polynomial GCD & Factorization](./poly-gcd-factoring.md) — square-free factorization and factorization algorithms
- [Solvers](../solvers.md) — polynomial-system solving and ODEs
