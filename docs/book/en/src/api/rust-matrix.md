# Matrices

> Source: `ocas-poly/src/matrix.rs`

`Matrix<D>` is the general dense matrix type in oCAS, with elements from an arbitrary `EuclideanDomain`. It uses row-major storage, and all elimination operations use **fraction-free** algorithms to avoid fraction blow-up of domain elements.

---

## `MatrixError`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    ShapeMismatch,
    RightHandSideIsNotVector,
    Inconsistent,
    Underdetermined { rank: usize },
    ResultNotInDomain,
}
```

Errors that may arise in matrix operations. Implements `Display` and `std::error::Error`.

### `ShapeMismatch`

**Meaning**: Matrix shapes are incompatible.

**Triggered when**: `trace` is used on a non-square matrix; `matmul` has `self.ncols != other.nrows`; `solve` has `self.nrows != b.len()`; `augment` has mismatched row counts.

### `RightHandSideIsNotVector`

**Meaning**: The right-hand side is not a column vector.

**Triggered when**: Not used directly by `solve` in the current version (`solve` accepts a `&[D::Element]` slice), but kept for extension.

### `Inconsistent`

**Meaning**: The linear system has no solution (inconsistent).

**Triggered when**: In `solve`, after elimination a free row of the augmented matrix has a nonzero entry in the last (constant) column.

### `Underdetermined { rank }`

**Meaning**: The linear system has infinitely many solutions (underdetermined).

**Triggered when**: In `solve`, the rank `rank` of the coefficient matrix is less than the number of variables `nvars`.

### `ResultNotInDomain`

**Meaning**: The solution does not lie in the target domain.

**Triggered when**: In the final division step of `solve`, `div` returns `None` (e.g., not divisible in an integer domain).

**Display output**:

| Variant | Output |
|---|---|
| `ShapeMismatch` | `"matrix shape mismatch"` |
| `RightHandSideIsNotVector` | `"right-hand side must be a column vector"` |
| `Inconsistent` | `"inconsistent linear system"` |
| `Underdetermined { rank }` | `"underdetermined system (rank N)"` |
| `ResultNotInDomain` | `"solution does not lie in the expected domain"` |

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;

// Inconsistent system: identical rows but different constant terms
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(1), Integer::from(2)];
assert_eq!(a.solve(&b), Err(MatrixError::Inconsistent));
```

**See also**: [solve](#solve)

---

## `Matrix<D>`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<D: EuclideanDomain> {
    data: Vec<D::Element>,   // row-major storage
    nrows: usize,
    ncols: usize,
    domain: D,
}
```

**Function**: A dense matrix over the domain `D`, with elements stored row-major in a `Vec<D::Element>`.

**Design invariants**:
- `data.len() == nrows * ncols` (guaranteed by constructors)
- Element `(i, j)` is at position `i * ncols + j` in `data`
- `D` is constrained to `EuclideanDomain`, supporting `div`, `gcd`, etc., which makes fraction-free elimination possible

**Trait implementations**:
- `Index<(usize, usize)>` / `IndexMut<(usize, usize)>` — element access via `matrix[(i, j)]`
- `Display` — one line per row, elements separated by spaces
- `Debug`, `Clone`, `PartialEq`, `Eq`

---

## Constructors

### `new`

**Signature**: `pub fn new(nrows: usize, ncols: usize, data: Vec<D::Element>, domain: D) -> Self`

**Function**: Creates a matrix from row-major one-dimensional data.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `nrows` | `usize` | Number of rows |
| `ncols` | `usize` | Number of columns |
| `data` | `Vec<D::Element>` | Elements in row-major order; length must equal `nrows * ncols` |
| `domain` | `D` | The coefficient domain |

**Return value**: `Matrix<D>`

**Errors**: Panics at runtime if `data.len() != nrows * ncols`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// Create a 2×3 matrix:
// [ 1  2  3 ]
// [ 4  5  6 ]
let m = Matrix::new(2, 3, vec![
    Integer::from(1), Integer::from(2), Integer::from(3),
    Integer::from(4), Integer::from(5), Integer::from(6),
], d);
assert_eq!(m.nrows(), 2);
assert_eq!(m.ncols(), 3);
assert_eq!(m[(0, 0)], Integer::from(1));
assert_eq!(m[(1, 2)], Integer::from(6));
```

**See also**: [from_rows](#from_rows), [zeros](#zeros), [identity](#identity)

---

### `zeros`

**Signature**: `pub fn zeros(nrows: usize, ncols: usize, domain: D) -> Self`

**Function**: Creates a zero matrix of the given shape.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `nrows` | `usize` | Number of rows |
| `ncols` | `usize` | Number of columns |
| `domain` | `D` | The coefficient domain |

**Return value**: A `Matrix<D>` with all elements equal to `domain.zero()`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let z = Matrix::zeros(2, 3, d);
assert_eq!(z[(0, 0)], Integer::from(0));
assert_eq!(z[(1, 2)], Integer::from(0));
```

**See also**: [identity](#identity), [new](#new)

---

### `identity`

**Signature**: `pub fn identity(n: usize, domain: D) -> Self`

**Function**: Creates the `n × n` identity matrix.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `n` | `usize` | Matrix dimension |
| `domain` | `D` | The coefficient domain |

**Return value**: A square matrix with `domain.one()` on the diagonal and `domain.zero()` elsewhere.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let id = Matrix::identity(3, d);
assert_eq!(id[(0, 0)], Integer::from(1));
assert_eq!(id[(0, 1)], Integer::from(0));
assert_eq!(id[(2, 2)], Integer::from(1));
```

**See also**: [zeros](#zeros)

---

### `from_rows`

**Signature**: `pub fn from_rows(rows: Vec<Vec<D::Element>>, domain: D) -> Self`

**Function**: Creates a matrix from nested row vectors.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rows` | `Vec<Vec<D::Element>>` | Each inner `Vec` is one row; all rows must have the same length |
| `domain` | `D` | The coefficient domain |

**Return value**: `Matrix<D>`

**Errors**: Panics at runtime if rows have inconsistent lengths.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(2), Integer::from(3)],
        vec![Integer::from(4), Integer::from(5), Integer::from(6)],
    ],
    d,
);
assert_eq!(m.nrows(), 2);
assert_eq!(m.ncols(), 3);
assert_eq!(m[(0, 2)], Integer::from(3));
```

**See also**: [new](#new), [into_rows](#into_rows)

---

## Accessor Methods

### `nrows`

**Signature**: `pub fn nrows(&self) -> usize`

**Function**: Returns the number of rows of the matrix.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 2, vec![Integer::from(0); 6], d);
assert_eq!(m.nrows(), 3);
```

---

### `ncols`

**Signature**: `pub fn ncols(&self) -> usize`

**Function**: Returns the number of columns of the matrix.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 2, vec![Integer::from(0); 6], d);
assert_eq!(m.ncols(), 2);
```

---

### `shape`

There is no separate `shape()` method. The shape is obtained by combining `nrows()` and `ncols()`:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 4, vec![Integer::from(0); 12], d);
let (rows, cols) = (m.nrows(), m.ncols());
assert_eq!((rows, cols), (3, 4));
```

---

### `domain`

**Signature**: `pub fn domain(&self) -> &D`

**Function**: Returns a reference to the matrix's coefficient domain.

**Example**:

```rust
use ocas_domain::{IntegerDomain, EuclideanDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::zeros(1, 1, d);
let dom = m.domain();
assert!(dom.is_zero(&dom.zero()));
```

---

### `data`

**Signature**: `pub fn data(&self) -> &[D::Element]`

**Function**: Returns a read-only slice of the underlying row-major storage.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let raw = m.data();
assert_eq!(raw.len(), 4);
assert_eq!(raw[0], Integer::from(1));  // (0,0)
assert_eq!(raw[3], Integer::from(4));  // (1,1)
```

---

### `row`

**Signature**: `pub fn row(&self, i: usize) -> Vec<D::Element>`

**Function**: Returns a copy of the `i`-th row.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `i` | `usize` | Row index, 0-based |

**Return value**: `Vec<D::Element>` of length `ncols`.

**Errors**: Panics if the index is out of bounds.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let r = m.row(1);
assert_eq!(r, vec![Integer::from(3), Integer::from(4)]);
```

**See also**: [column](#column), [into_rows](#into_rows)

---

### `column`

**Signature**: `pub fn column(&self, j: usize) -> Vec<D::Element>`

**Function**: Returns a copy of the `j`-th column.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `j` | `usize` | Column index, 0-based |

**Return value**: `Vec<D::Element>` of length `nrows`.

**Errors**: Panics if the index is out of bounds.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let c = m.column(0);
assert_eq!(c, vec![Integer::from(1), Integer::from(3)]);
```

**See also**: [row](#row)

---

### `into_rows`

**Signature**: `pub fn into_rows(self) -> Vec<Vec<D::Element>>`

**Function**: Consumes the matrix and converts it into nested row vectors.

**Return value**: Outer `Vec` of length `nrows`, each inner `Vec` of length `ncols`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let rows = m.into_rows();
assert_eq!(rows.len(), 2);
assert_eq!(rows[0], vec![Integer::from(1), Integer::from(2)]);
assert_eq!(rows[1], vec![Integer::from(3), Integer::from(4)]);
```

**See also**: [from_rows](#from_rows), [row](#row)

---

## Indexing & Display

### `Index<(usize, usize)>`

**Signature**: `impl Index<(usize, usize)> for Matrix<D>`

**Function**: Read-only access to element `(i, j)` via `matrix[(i, j)]`, equivalent to `data[i * ncols + j]`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(5), Integer::from(6)]], d,
);
assert_eq!(m[(0, 0)], Integer::from(5));
assert_eq!(m[(0, 1)], Integer::from(6));
```

---

### `IndexMut<(usize, usize)>`

**Signature**: `impl IndexMut<(usize, usize)> for Matrix<D>`

**Function**: Writable access to elements via `matrix[(i, j)] = val`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut m = Matrix::zeros(2, 2, d);
m[(0, 1)] = Integer::from(42);
m[(1, 0)] = Integer::from(-7);
assert_eq!(m[(0, 1)], Integer::from(42));
assert_eq!(m[(1, 0)], Integer::from(-7));
```

---

### `Display`

**Signature**: `impl Display for Matrix<D> where D::Element: Display`

**Function**: Formats the matrix with one row per line, elements separated by spaces.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let s = format!("{:?}", m);
// Output: Debug representation (IntegerDomain does not implement Display, and Matrix's Display requires it)
assert!(s.contains("1"));
```

---

## Matrix Operations

### `swap_rows`

**Signature**: `pub fn swap_rows(&mut self, i: usize, j: usize, start_col: usize)`

**Function**: Swaps row `i` and row `j`, starting from column `start_col`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `i` | `usize` | Index of the first row |
| `j` | `usize` | Index of the second row |
| `start_col` | `usize` | Start column (inclusive); elements before it are unchanged |

**Return value**: None (in-place modification).

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
m.swap_rows(0, 1, 0);
assert_eq!(m[(0, 0)], Integer::from(3));
assert_eq!(m[(0, 1)], Integer::from(4));
assert_eq!(m[(1, 0)], Integer::from(1));
assert_eq!(m[(1, 1)], Integer::from(2));
```

---

### `transpose`

**Signature**: `pub fn transpose(&self) -> Matrix<D>`

**Function**: Returns the transpose of the matrix. $(A^T)_{ij} = A_{ji}$.

**Return value**: A new `Matrix<D>` of shape `(ncols, nrows)`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2), Integer::from(3)],
         vec![Integer::from(4), Integer::from(5), Integer::from(6)]],
    d,
);
let t = a.transpose();
assert_eq!(t.nrows(), 3);
assert_eq!(t.ncols(), 2);
assert_eq!(t[(0, 0)], Integer::from(1));
assert_eq!(t[(0, 1)], Integer::from(4));
assert_eq!(t[(2, 1)], Integer::from(6));
```

**See also**: [trace](#trace)

---

### `trace`

**Signature**: `pub fn trace(&self) -> Result<D::Element, MatrixError>`

**Function**: Computes the trace of a square matrix, i.e. the sum of the diagonal elements: $\text{tr}(A) = \sum_{i=0}^{n-1} a_{ii}$.

**Return value**:
- `Ok(sum)` — the trace
- `Err(MatrixError::ShapeMismatch)` — the matrix is not square

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
// tr = 1 + 4 = 5
assert_eq!(a.trace().unwrap(), Integer::from(5));
```

**Non-square error**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let m = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
assert_eq!(m.trace(), Err(MatrixError::ShapeMismatch));
```

**See also**: [determinant](#determinant)

---

### `matmul`

**Signature**: `pub fn matmul(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError>`

**Function**: Computes the matrix product $C = AB$, where $C_{ij} = \sum_k A_{ik} \cdot B_{kj}$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Matrix<D>` | The matrix to multiply on the right |

**Return value**:
- `Ok(Matrix)` — the product matrix of shape `(self.nrows, other.ncols)`
- `Err(MatrixError::ShapeMismatch)` — `self.ncols != other.nrows`

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let b = Matrix::from_rows(
    vec![vec![Integer::from(5), Integer::from(6)],
         vec![Integer::from(7), Integer::from(8)]],
    d,
);
let c = a.matmul(&b).unwrap();
// C = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]]
//   = [[19, 22], [43, 50]]
assert_eq!(c[(0, 0)], Integer::from(19));
assert_eq!(c[(0, 1)], Integer::from(22));
assert_eq!(c[(1, 0)], Integer::from(43));
assert_eq!(c[(1, 1)], Integer::from(50));
```

**Shape mismatch**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
let b = Matrix::from_rows(vec![vec![Integer::from(3), Integer::from(4)]], d);
// a: 1×2, b: 1×2 → ncols(a)=2 ≠ nrows(b)=1
assert_eq!(a.matmul(&b), Err(MatrixError::ShapeMismatch));
```

**See also**: [solve](#solve)

---

### `determinant`

**Signature**: `pub fn determinant(&self) -> Result<D::Element, MatrixError>`

**Function**: Computes the determinant of a square matrix using the **Bareiss fraction-free algorithm** (with partial pivoting).

**Algorithm highlights**:
- Reduces the matrix to upper triangular form by elimination, dividing at each step by the previous pivot (Bareiss guarantees exact divisibility)
- Partial pivoting: when a zero pivot is encountered, searches downward for a nonzero row and swaps
- If an entire column becomes zero during elimination, the determinant is zero
- The sign is determined by the parity of the number of row swaps

**Return value**:
- `Ok(det)` — the determinant
- `Err(MatrixError::ShapeMismatch)` — the matrix is not square

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 2×2: det = 1*4 - 2*3 = -2
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
assert_eq!(a.determinant().unwrap(), Integer::from(-2));
```

**Singular matrix**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// The second row is 2× the first row
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(2), Integer::from(4)]],
    d,
);
assert_eq!(a.determinant().unwrap(), Integer::from(0));
```

**See also**: [inverse](#inverse), [rank](#rank)

---

### `rank`

**Signature**: `pub fn rank(&self) -> usize`

**Function**: Computes the rank of the matrix via fraction-free Gaussian elimination. Clones the matrix, runs `row_echelon`, and returns the number of nonzero rows.

**Return value**: `usize` — the rank of the matrix.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// Full-rank 2×2
let id = Matrix::identity(2, d);
assert_eq!(id.rank(), 2);

// Rank-deficient: second row = 2 × first row
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2), Integer::from(3)],
         vec![Integer::from(2), Integer::from(4), Integer::from(6)],
         vec![Integer::from(1), Integer::from(1), Integer::from(1)]],
    d,
);
assert_eq!(a.rank(), 2);
```

**See also**: [determinant](#determinant), [row_echelon](#row_echelon)

---

### `inverse`

**Signature**: `pub fn inverse(&self) -> Result<Matrix<D>, MatrixError>`

**Function**: Computes the inverse matrix $A^{-1}$ of a square matrix, by solving $A x_j = e_j$ column by column ($e_j$ is the $j$-th column of the identity matrix).

**Return value**:
- `Ok(Matrix)` — the inverse, satisfying $A \cdot A^{-1} = I$
- `Err(MatrixError::ShapeMismatch)` — the matrix is not square
- `Err(MatrixError::Inconsistent)` — the matrix is singular (no solution)
- `Err(MatrixError::Underdetermined { .. })` — singularity leads to an underdetermined system
- `Err(MatrixError::ResultNotInDomain)` — an entry of the inverse is not in the domain

**⚠️ Note**: Over the integer domain `IntegerDomain`, only unimodular matrices ($\det = \pm 1$) have integer inverses. Invert general matrices over the rational domain `RationalDomain`.

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// Unimodular matrix: det = 1*1 - 0*2 = 1
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(0), Integer::from(1)]],
    d,
);
let inv = a.inverse().unwrap();
assert_eq!(inv[(0, 0)], Integer::from(1));
assert_eq!(inv[(0, 1)], Integer::from(-2));
assert_eq!(inv[(1, 0)], Integer::from(0));
assert_eq!(inv[(1, 1)], Integer::from(1));

// Verify A * A^{-1} = I
let prod = a.matmul(&inv).unwrap();
assert_eq!(prod, Matrix::identity(2, IntegerDomain));
```

**Singular matrix error**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(2), Integer::from(4)]],
    d,
);
assert!(a.inverse().is_err());
```

**See also**: [solve](#solve), [determinant](#determinant)

---

### `augment`

**Signature**: `pub fn augment(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError>`

**Function**: Horizontally concatenates two matrices, producing the augmented matrix $[A \mid B]$.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `other` | `&Matrix<D>` | The matrix to concatenate on the right |

**Return value**:
- `Ok(Matrix)` — the augmented matrix of shape `(nrows, ncols_a + ncols_b)`
- `Err(MatrixError::ShapeMismatch)` — the two matrices have different row counts

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(3), Integer::from(2)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = Matrix::new(2, 1, vec![Integer::from(5), Integer::from(2)], d);
let aug = a.augment(&b).unwrap();
assert_eq!(aug.nrows(), 2);
assert_eq!(aug.ncols(), 3);
assert_eq!(aug[(0, 0)], Integer::from(3));
assert_eq!(aug[(0, 2)], Integer::from(5));
```

**See also**: [solve](#solve)

---

### `solve`

**Signature**: `pub fn solve(&self, b: &[D::Element]) -> Result<Vec<D::Element>, MatrixError>`

**Function**: Solves the linear system $Ax = b$, using fraction-free Gaussian elimination with back substitution.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `b` | `&[D::Element]` | The right-hand side vector; length must equal `self.nrows` |

**Return value**:
- `Ok(Vec<D::Element>)` — the solution vector of length `self.ncols`
- `Err(MatrixError::ShapeMismatch)` — `b.len() != self.nrows`
- `Err(MatrixError::Inconsistent)` — the system has no solution
- `Err(MatrixError::Underdetermined { rank })` — infinitely many solutions, `rank < ncols`
- `Err(MatrixError::ResultNotInDomain)` — the solution is not in the domain (e.g., integer division fails)

**Algorithm steps**:
1. Build the augmented matrix $[A \mid b]$
2. Run fraction-free row echelon reduction (`row_echelon`) on the first `ncols` columns
3. Check whether the constant terms of the rows below the rank are zero (consistency check)
4. Run fraction-free back substitution (`back_substitution`)
5. If the rank < number of variables, return `Underdetermined`
6. Divide by the pivots to obtain the solution

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 2x + y = 4
// x  + y = 3  → x = 1, y = 2
let a = Matrix::from_rows(
    vec![vec![Integer::from(2), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(4), Integer::from(3)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2)]);
```

**3×3 system**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// x + y + z = 6
// 2x - y + z = 3
// x + 2y - z = 2  → x=1, y=2, z=3
let a = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(1), Integer::from(1)],
        vec![Integer::from(2), Integer::from(-1), Integer::from(1)],
        vec![Integer::from(1), Integer::from(2), Integer::from(-1)],
    ],
    d,
);
let b = vec![Integer::from(6), Integer::from(3), Integer::from(2)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2), Integer::from(3)]);
```

**Inconsistent system**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(1), Integer::from(2)];
assert_eq!(a.solve(&b), Err(MatrixError::Inconsistent));
```

**Underdetermined system**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(3)];
assert!(matches!(a.solve(&b), Err(MatrixError::Underdetermined { .. })));
```

**See also**: [inverse](#inverse), [augment](#augment)

---

## Internal Elimination Methods

The following two methods are `pub`; they usually do not need to be called directly, but may be useful in advanced scenarios (e.g., custom elimination pipelines).

### `row_echelon`

**Signature**: `pub fn row_echelon(&mut self, max_col: usize) -> usize`

**Function**: Runs fraction-free Gaussian elimination on the first `max_col` columns, reducing the matrix to row echelon form. Modifies the matrix in place. Returns the rank (number of nonzero rows).

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `max_col` | `usize` | The maximum number of columns involved in the elimination (exclusive); automatically clamped to `ncols` |

**Return value**: `usize` — the rank of the matrix.

**Algorithm details**:
- Finds a nonzero pivot in each column, swapping rows if necessary
- Extracts the GCD content from the pivot row to prevent coefficient growth
- Elimination uses gcd-reduced fraction-free scaling: with $g = \gcd(\text{pivot}, a_{kj})$, the eliminated row is scaled by $\text{pivot} / g$ and the pivot row by $a_{kj} / g$ before subtracting, so no fractions appear

**Example**:

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut a = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(2), Integer::from(3)],
        vec![Integer::from(2), Integer::from(4), Integer::from(6)],
        vec![Integer::from(0), Integer::from(1), Integer::from(1)],
    ],
    d,
);
let rank = a.row_echelon(3);
assert_eq!(rank, 2);
```

**See also**: [back_substitution](#back_substitution), [rank](#rank)

---

### `back_substitution`

**Signature**: `pub fn back_substitution(&mut self, max_col: usize)`

**Function**: Runs fraction-free back substitution on a matrix already in row echelon form. Modifies the first `max_col` columns of the matrix in place.

**Algorithm details**:
- Processes rows from the bottom up
- For each pivot row, first extracts the GCD content
- For rows above the pivot, eliminates the pivot column with the same gcd-reduced scaling (the row above is scaled by $\text{pivot} / g$ and the pivot row by $a_{kj} / g$, where $g = \gcd(\text{pivot}, a_{kj})$) before subtracting

**See also**: [row_echelon](#row_echelon), [solve](#solve)

---

## Complete Example

The following comprehensive example demonstrates the basic matrix workflow: construction, operations, and solving.

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

fn main() {
    let d = IntegerDomain;

    // Construct matrix A
    let a = Matrix::from_rows(
        vec![
            vec![Integer::from(2), Integer::from(1), Integer::from(-1)],
            vec![Integer::from(-3), Integer::from(-1), Integer::from(2)],
            vec![Integer::from(-2), Integer::from(1), Integer::from(2)],
        ],
        d,
    );

    // Determinant
    let det = a.determinant().unwrap();
    assert_eq!(det, Integer::from(-1)); // det = -1 (unimodular)

    // Rank
    assert_eq!(a.rank(), 3);

    // Trace
    let tr = a.trace().unwrap();
    assert_eq!(tr, Integer::from(3)); // 2 + (-1) + 2

    // Transpose
    let at = a.transpose();
    assert_eq!(at[(0, 1)], Integer::from(-3));

    // Solve Ax = b
    let b = vec![Integer::from(8), Integer::from(-11), Integer::from(-3)];
    let x = a.solve(&b).unwrap();
    assert_eq!(x, vec![Integer::from(2), Integer::from(3), Integer::from(-1)]);

    // Inverse (det = -1, so an integer inverse exists)
    let inv = a.inverse().unwrap();
    // Verify A * A^{-1} = I
    let prod = a.matmul(&inv).unwrap();
    assert_eq!(prod, Matrix::identity(3, IntegerDomain));

    // Matrix multiplication
    let c = a.matmul(&a).unwrap();
    assert_eq!(c.nrows(), 3);
    assert_eq!(c.ncols(), 3);
}
```

---

## See also

- [Expression System](./rust-expressions.md) — expression tree construction and manipulation
- [Coefficient Domains](./rust-domains.md) — the `EuclideanDomain` trait and domain implementations
- [Polynomials](./rust-polynomials.md) — polynomial types using the same domain traits
- [Linear Algebra Foundations](../math/linear-algebra.md) — the mathematics of the Bareiss algorithm and Gaussian elimination
