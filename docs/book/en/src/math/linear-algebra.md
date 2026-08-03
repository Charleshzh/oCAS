# Foundations: Linear Algebra

This chapter presents the mathematical foundations of matrix operations in oCAS. The `Matrix` type in oCAS works over a `EuclideanDomain` — not just floats or rationals, but any algebraic structure supporting Euclidean division. This means we need an elimination algorithm that **does not introduce fractions**.

---

## Prerequisites

Before reading this chapter, you should be familiar with:

- Basic notions of rings and fields (see [Finite Fields and Modular Arithmetic](./finite-fields.md))
- Euclidean domains: algebraic structures supporting `add`, `mul`, `sub`, `neg`, `div` (division with remainder), `gcd`, and `rem`
- Basic intuition for vector spaces: vector addition and scalar multiplication

No systematic linear algebra training is required — this chapter builds everything from scratch.

---

## Basic Concepts

### Matrices

An $m \times n$ **matrix** is a rectangular array of elements arranged in $m$ rows and $n$ columns:

$$
A = \begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

The element $a_{ij}$ sits in row $i$ and column $j$. In oCAS, elements are stored in a row-major one-dimensional array: `data[i * ncols + j]` corresponds to $a_{i+1,\,j+1}$ (zero-indexed).

### Matrix Operations

**Addition**: two matrices of the same shape are added elementwise:

$$
(A + B)_{ij} = a_{ij} + b_{ij}
$$

**Scalar multiplication**:

$$
(cA)_{ij} = c \cdot a_{ij}
$$

**Matrix multiplication**: the product $C = AB$ of an $m \times n$ matrix $A$ and an $n \times p$ matrix $B$ is the $m \times p$ matrix:

$$
c_{ij} = \sum_{k=1}^{n} a_{ik} \cdot b_{kj}
$$

Multiplication requires the number of columns of $A$ to equal the number of rows of $B$. Note that matrix multiplication is generally **not commutative** ($AB \neq BA$), but it is **associative** ($(AB)C = A(BC)$).

**Transpose**: swapping rows and columns:

$$
(A^\top)_{ij} = a_{ji}
$$

Properties: $(AB)^\top = B^\top A^\top$, $(A^\top)^\top = A$.

**Trace**: the sum of the diagonal elements of a square matrix:

$$
\operatorname{tr}(A) = \sum_{i=1}^{n} a_{ii}
$$

Property: $\operatorname{tr}(AB) = \operatorname{tr}(BA)$ (whenever the product is defined).

### Determinants

The **determinant** is a scalar value attached to a square matrix, one of the most central concepts in linear algebra. For an $n \times n$ matrix $A$, the determinant $\det(A)$ has the following equivalent definitions:

**Recursive definition (Laplace expansion)**: expanding along row $i$:

$$
\det(A) = \sum_{j=1}^{n} (-1)^{i+j} \, a_{ij} \, M_{ij}
$$

where $M_{ij}$ is the determinant of the $(n-1)\times(n-1)$ submatrix obtained by deleting row $i$ and column $j$ (the **minor**). For a $1\times 1$ matrix, $\det(a_{11}) = a_{11}$.

**Geometric definition**: $|\det(A)|$ is the "oriented volume" of the parallelepiped spanned by the column vectors of the matrix. $\det(A) > 0$ means orientation is preserved; $\det(A) < 0$ means orientation is flipped.

#### Key Properties of Determinants

1. **Multiplicativity**: $\det(AB) = \det(A) \cdot \det(B)$
2. **Transpose invariance**: $\det(A^\top) = \det(A)$
3. **Row swap changes sign**: swapping two rows changes the sign of the determinant
4. **Row scaling**: multiplying a row by $c$ multiplies the determinant by $c$
5. **Row addition invariance**: adding a multiple of one row to another leaves the determinant unchanged
6. **Singularity test**: $\det(A) = 0$ if and only if $A$ is not invertible (singular)
7. **Upper/lower triangular matrices**: $\det(A) = \prod_{i=1}^{n} a_{ii}$ (the product of the diagonal elements)

These properties directly suggest computing determinants by **elimination**: transform the matrix into triangular form by row operations, and the determinant is the product of the diagonal elements (times a sign factor from row swaps).

---

## Core Theory

### The Problem: Fraction Explosion over Domains

Classical Gaussian elimination works well over the rationals — we can divide freely. But `Matrix<D: EuclideanDomain>` in oCAS is required to work over **arbitrary** Euclidean domains, including:

- **The integer ring $\mathbb{Z}$**: $3 \div 2$ is meaningless (no integer quotient)
- **Polynomial rings $\mathbb{F}[x]$**: $x^2 \div (x+1)$ is not a polynomial
- **Finite fields $\mathbb{F}_p$**: division is possible, but every division incurs a modular inverse, and intermediate results may overflow

Even over the rationals, classical elimination produces ever-growing denominators. For example, for the integer matrix:

$$
A = \begin{pmatrix} 2 & 3 \\ 5 & 7 \end{pmatrix}
$$

classical elimination needs to compute $7 - \frac{5 \cdot 3}{2} = 7 - \frac{15}{2} = -\frac{1}{2}$ — a fraction appears. For an $n \times n$ matrix, the numerators and denominators of the fractions can grow to magnitude $O(n!)$.

### Bareiss's Fraction-Free Determinant Algorithm

**Bareiss's algorithm** (1968) is built on the key insight that during elimination over the integers (or a more general Euclidean domain), all intermediate values can be kept integral by dividing by the previous pivot.

#### Algorithm Description

For an $n \times n$ matrix $A$, Bareiss's algorithm produces an upper triangular matrix whose last diagonal element $a_{nn}^{(n-1)}$ is $\det(A)$ (possibly multiplied by $-1$ to compensate for row swaps).

Initialization: set $a_{ij}^{(0)} = a_{ij}$ and $p_0 = 1$.

For $k = 0, 1, \ldots, n-2$:

1. **Pivot selection**: if $a_{kk}^{(k)} = 0$, find $i > k$ with $a_{ik}^{(k)} \neq 0$, and swap row $k$ with row $i$ (record the sign flip). If none exists, $\det(A) = 0$.

2. **Elimination**: for $i = k+1, \ldots, n-1$ and $j = k+1, \ldots, n-1$:

$$
a_{ij}^{(k+1)} = \frac{a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}}{p_k}
$$

3. **Update the divisor**: $p_{k+1} = a_{kk}^{(k)}$

Finally:

$$
\det(A) = \begin{cases}
a_{nn}^{(n-1)} & \text{if no row swap or an even number of swaps} \\
-a_{nn}^{(n-1)} & \text{if an odd number of row swaps}
\end{cases}
$$

#### Why the Division Is Exact

**Theorem** (Bareiss 1968): in the formula above, $p_k$ divides the numerator $a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}$ exactly.

Intuitively, this follows from Sylvester's determinant identity. The step-$(k+1)$ update of the elimination matrix is equivalent to computing the determinant of a $(k+2) \times (k+2)$ submatrix, and the Schur complement formula guarantees that dividing by the previous pivot still yields an integer.

#### Complexity and Numerical Behavior

- **Time complexity**: $O(n^3)$ domain operations — the same as classical Gaussian elimination
- **Size of intermediate values**: for integer matrices, the intermediate values $a_{ij}^{(k)}$ grow to $O(n \log n)$ digits — far better than the magnitude-$O(n!)$ growth of classical elimination
- **No rounding errors**: over the integers and polynomials, all operations are exact, with no floating-point approximation

### Gaussian Elimination and Row Echelon Form

Gaussian elimination transforms a matrix into **row echelon form** (REF):

$$
\begin{pmatrix}
\boxed{a_{11}} & a_{12} & a_{13} & a_{14} \\
0 & \boxed{a_{22}} & a_{23} & a_{24} \\
0 & 0 & 0 & \boxed{a_{34}} \\
0 & 0 & 0 & 0
\end{pmatrix}
$$

The boxed elements are the **pivots**. Row echelon form satisfies:
1. All nonzero rows are above all zero rows
2. The column of each pivot is all zeros below it
3. Each pivot lies strictly to the right of the pivots in the rows above it

Continuing the elimination yields the **reduced row echelon form** (RREF) — every pivot is $1$ and all remaining elements in pivot columns are $0$.

#### Pivot Selection

A Euclidean domain need not have a notion of "absolute value" (for instance $\mathbb{F}[x]$), so oCAS's elimination does not use the **partial pivoting** of numerical linear algebra. Both `row_echelon` and `determinant` select the **first nonzero element** in the current column as the pivot; if the column is entirely zero, it is skipped (no row swap). Intermediate value growth is controlled by **content stripping** (see below) — after elimination, the GCD of the pivot row is immediately factored out and the whole row is divided by it.

### Fraction-Free Back-Substitution

After elimination, **back-substitution** is needed to solve linear systems. oCAS's back-substitution is also fraction-free: starting from the last row and working upward, GCD scaling is used to avoid introducing fractions.

For row $i$ in row echelon form (pivot in column $j$) and row $k < i$ (with a nonzero element $a_{kj}$ in column $j$):

$$
g = \gcd(a_{ij}, a_{kj}), \quad s_p = \frac{a_{kj}}{g}, \quad s_r = \frac{a_{ij}}{g}
$$

Update row $k$: for every $l > j$:

$$
a_{kl} \leftarrow s_r \cdot a_{kl} - s_p \cdot a_{il}
$$

$$
a_{kj} \leftarrow 0
$$

This is fully consistent with the fraction-free nature of Bareiss elimination.

### Content Stripping

During elimination, the coefficients of a row may share a common factor. At each step after processing a pivot, oCAS performs **content stripping**: it computes the GCD of the pivot row from the pivot column to the last column, then divides the entire row by that GCD.

$$
g = \gcd(a_{ij}, a_{i,j+1}, \ldots, a_{i,n-1})
$$

$$
a_{il} \leftarrow \frac{a_{il}}{g}, \quad l = j, j+1, \ldots, n-1
$$

This optimization is especially critical for polynomial rings — otherwise the degrees of the intermediate polynomials would grow exponentially.

### The Inverse of a Matrix

The **inverse** $A^{-1}$ of a square matrix $A$ satisfies:

$$
AA^{-1} = A^{-1}A = I
$$

where $I$ is the identity matrix. The inverse exists if and only if $\det(A) \neq 0$ ($A$ is nonsingular).

Computation method: transform $[A \mid I]$ by row operations into $[I \mid A^{-1}]$. In oCAS, `inverse()` obtains the $j$-th column of $A^{-1}$ by solving $Ax_j = e_j$ for each column $e_j$ of the identity matrix.

---

## Implementation in oCAS

### The `Matrix<D: EuclideanDomain>` Structure

The oCAS matrix is **generic** — it accepts any `EuclideanDomain` as the coefficient domain:

```rust
use ocas_domain::{EuclideanDomain, IntegerDomain, Integer};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::new(2, 2, vec![
    Integer::from(2), Integer::from(3),
    Integer::from(5), Integer::from(7),
], d);
```

Internally it is stored as a row-major `Vec<D::Element>` together with the number of rows, the number of columns, and the coefficient domain value (`domain: D` is stored by value, not by reference). This way all arithmetic operations go through the domain's trait methods (`add`, `mul`, `sub`, `div`, `gcd`) and never depend on the concrete type.

### Why Bareiss and Not Laplace Expansion

| Method | Time complexity | Space complexity | Fraction explosion | Domain compatibility |
|---|---|---|---|---|
| Laplace expansion | $O(n!)$ | $O(n^2)$ (recursion depth $O(n)$) | None | Any commutative ring (infeasible in practice) |
| Classical Gaussian elimination | $O(n^3)$ | $O(n^2)$ | Severe | Fields only |
| **Bareiss** | $O(n^3)$ | $O(n^2)$ | **None** | **Euclidean domains** |

- The $O(n!)$ complexity of **Laplace expansion** makes it completely unusable for $n > 10$
- **Classical elimination** requires fraction arithmetic and does not work over $\mathbb{Z}$ or $\mathbb{F}[x]$
- **Bareiss's algorithm** combines the best of both: $O(n^3)$ efficiency + fraction-free + support for arbitrary Euclidean domains

This is the fundamental reason oCAS chooses Bareiss — it is a **general-purpose symbolic computation system**, whose matrix elements may be integers, polynomials, finite-field elements, or even algebraic number field elements.

### Row Echelon Form: `row_echelon`

```rust
// fraction-free Gaussian elimination, returns the rank
pub fn row_echelon(&mut self, max_col: usize) -> usize
```

`row_echelon` performs fraction-free elimination on the first `max_col` columns and returns the rank of the matrix. Algorithm steps:

1. For each column $j$, select the first nonzero element in the current column as the pivot (skip the column if it is entirely zero)
2. **Content stripping**: compute the GCD of the pivot row and divide the whole row by it
3. **Elimination**: for each row below the pivot, eliminate column $j$ by GCD-scaled subtraction ($a_{ij}$ is the pivot and $a_{kj}$ the entry being eliminated; $s_r$ and $s_p$ must be computed before zeroing):

$$
g = \gcd(a_{ij}, a_{kj}), \quad s_r = \frac{a_{ij}}{g}, \quad s_p = \frac{a_{kj}}{g}
$$

$$
a_{kj} \leftarrow 0, \qquad a_{kl} \leftarrow s_r \cdot a_{kl} - s_p \cdot a_{il} \quad (l > j)
$$

This guarantees that all intermediate values stay inside the domain, with no fractions produced.

### Determinant: `determinant`

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(2)],
    vec![Integer::from(3), Integer::from(4)],
], d);
assert_eq!(a.determinant().unwrap(), Integer::from(-2));
```

`determinant` implements Bareiss's algorithm directly:

1. Initialization: copy the data into a working array, $p_0 = 1$, sign flag $\epsilon = 1$
2. For $k = 0, \ldots, n-2$:
   - If the pivot is zero, find a nonzero row below it and swap (flip $\epsilon$)
   - Update: $a_{ij}^{(k+1)} = (a_{ij}^{(k)} \cdot a_{kk}^{(k)} - a_{ik}^{(k)} \cdot a_{kj}^{(k)}) / p_k$
   - Update the divisor: $p_{k+1} = a_{kk}^{(k)}$
3. Return $\epsilon \cdot a_{n-1,n-1}^{(n-1)}$

Special cases: $n = 0$ returns $1$ (the determinant of the empty matrix), and $n = 1$ returns the single element directly.

### Solving Linear Systems: `solve`

```rust
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(1)],
    vec![Integer::from(1), Integer::from(-1)],
], d);
let b = vec![Integer::from(3), Integer::from(-1)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2)]);
```

The complete flow of `solve`:

1. **Augmentation**: append $b$ to form $[A \mid b]$
2. **Elimination**: `row_echelon(nvars)` to reach row echelon form
3. **Consistency check**: if there is a nonzero $b$ component after the pivot rows, the system is inconsistent and `MatrixError::Inconsistent` is returned
4. **Back-substitution**: fraction-free `back_substitution(nvars)`
5. **Rank check**: if the rank $< n$, the system is underdetermined and `MatrixError::Underdetermined` is returned
6. **Final division**: divide each solution component by its corresponding pivot

#### Error Types

| Error | Meaning |
|---|---|
| `ShapeMismatch` | $A$ is not square, or the length of $b$ is not equal to the number of rows |
| `RightHandSideIsNotVector` | $b$ is not a column vector |
| `Inconsistent` | The system has no solution (rank of the augmented matrix > rank of the coefficient matrix) |
| `Underdetermined { rank }` | The system has infinitely many solutions (rank < number of unknowns) |
| `ResultNotInDomain` | The solution is not in the domain (an integer matrix may produce rational solutions) |

### Rank and Inverse

```rust
// rank: the number of pivots in row echelon form
let rank = a.rank();

// inverse: solve Ax_j = e_j column by column
let inv = a.inverse()?;
```

`rank` clones the matrix, calls `row_echelon`, and returns the number of pivots.

`inverse` solves a linear system for each column of the identity matrix to obtain the columns of the inverse. This requires the matrix to be nonsingular and the solutions to lie in the domain — for example, the inverse of an integer matrix is generally not in $\mathbb{Z}$ unless $\det(A) = \pm 1$ (a unimodular matrix):

```rust
// unimodular matrix: det = 1, the integer inverse exists
let a = Matrix::from_rows(vec![
    vec![Integer::from(1), Integer::from(2)],
    vec![Integer::from(0), Integer::from(1)],
], d);
let inv = a.inverse().unwrap();
assert_eq!(inv[(0, 0)], Integer::from(1));
assert_eq!(inv[(0, 1)], Integer::from(-2));
assert_eq!(inv[(1, 0)], Integer::from(0));
assert_eq!(inv[(1, 1)], Integer::from(1));
```

### Matrix Multiplication and Transpose

```rust
// multiplication: C = A · B
let c = a.matmul(&b)?;

// transpose
let at = a.transpose();

// trace (sum of the diagonal elements of a square matrix)
let tr = a.trace()?;
```

`matmul` uses the standard triple-loop $O(n^3)$ algorithm; all operations go through the domain's trait methods. `transpose` is implemented by rearranging the data. `trace` sums the diagonal elements of a square matrix.

### Complete Example: Solving an Integer Linear System

```rust
use ocas_domain::{IntegerDomain, Integer};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// System of equations:
//   2x + 3y = 8
//   5x + 7y = 19
let a = Matrix::from_rows(vec![
    vec![Integer::from(2), Integer::from(3)],
    vec![Integer::from(5), Integer::from(7)],
], d);
let b = vec![Integer::from(8), Integer::from(19)];

let x = a.solve(&b).unwrap();
// x = [1, 2]
assert_eq!(x[0], Integer::from(1));
assert_eq!(x[1], Integer::from(2));

// Verify: Ax = b
let x_matrix = Matrix::new(2, 1, x, d);
let result = a.matmul(&x_matrix).unwrap();
assert_eq!(result[(0, 0)], Integer::from(8));
assert_eq!(result[(1, 0)], Integer::from(19));

// determinant: 2·7 - 3·5 = -1
assert_eq!(a.determinant().unwrap(), Integer::from(-1));
```

All intermediate values during the elimination remain integers — Bareiss's algorithm guarantees this.

---

## Advanced Topics

### Generalizations of Bareiss's Algorithm

Bareiss's algorithm can be generalized to compute **determinants of submatrices** (all $k \times k$ minors), which is useful for computing Smith normal forms and Hermite normal forms.

For polynomial matrices, Bareiss effectively controls coefficient growth through GCD scaling, making determinant computation of polynomial matrices feasible in practice.

### Relation to LU Decomposition

Bareiss's algorithm can be viewed as a **fraction-free variant of LU decomposition**. Classical LU decomposition factors the matrix as $A = LU$, where $L$ is lower triangular and $U$ is upper triangular. Bareiss produces a similar decomposition, but all intermediate values stay within the domain.

$\det(A) = \det(L) \cdot \det(U) = \prod_{i} l_{ii} \cdot \prod_{i} u_{ii}$.

### Matrices in Symbolic Computation

In symbolic computation, matrix elements may be **polynomials** or **rational functions**. In this setting the GCD scaling of Bareiss's algorithm is especially critical — it keeps the degree growth of polynomial coefficients at $O(n)$ instead of exponential.

The F4 algorithm of oCAS's Gröbner basis machinery makes heavy internal use of matrix row echelonization — reducing the elimination problem for multivariate polynomials to a matrix problem over the coefficient domain, then processing it efficiently with Bareiss or its $\mathbb{Z}_p$ fast path (lazy modular arithmetic).

---

## References

- **Axler, S.** *Linear Algebra Done Right*, 3rd ed. Springer, 2015.
  Chapters 4–6: determinants, eigenvalues, and inner product spaces.

- **Bareiss, E. H.** "Sylvester's Identity and Multistep
  Integer-Preserving Gaussian Elimination." *Mathematics of Computation*,
  22(103):565–578, 1968.
  The original paper on Bareiss's algorithm.

- **Cox, D., Little, J., O'Shea, D.** *Ideals, Varieties, and Algorithms*,
  4th ed. Springer, 2015.
  Linear algebra in polynomial matrices and Gröbner bases.

- **Geddes, K. O., Czapor, S. R., Labahn, G.** *Algorithms for Computer
  Algebra*. Kluwer, 1992.
  Chapter 5: symbolic elimination algorithms on matrices.
