//! Matrix types and linear algebra over algebraic domains.
//!
//! Provides a generic [`Matrix`] type parameterized by any [`EuclideanDomain`].
//! Supports Gaussian elimination, back-substitution, and solving linear
//! systems Ax = b with fraction-free arithmetic to avoid coefficient blow-up.

use std::fmt;
use std::ops::{Index, IndexMut};

use ocas_domain::EuclideanDomain;

/// Errors that can occur during matrix operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    /// Matrix shapes are incompatible for the requested operation.
    ShapeMismatch,
    /// The right-hand side is not a column vector.
    RightHandSideIsNotVector,
    /// The linear system is inconsistent (no solution).
    Inconsistent,
    /// The system is underdetermined (infinitely many solutions).
    Underdetermined {
        /// Rank of the coefficient matrix.
        rank: usize,
    },
    /// The solution does not lie in the expected domain.
    ResultNotInDomain,
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::ShapeMismatch => f.write_str("matrix shape mismatch"),
            MatrixError::RightHandSideIsNotVector => {
                f.write_str("right-hand side must be a column vector")
            }
            MatrixError::Inconsistent => f.write_str("inconsistent linear system"),
            MatrixError::Underdetermined { rank } => {
                write!(f, "underdetermined system (rank {})", rank)
            }
            MatrixError::ResultNotInDomain => {
                f.write_str("solution does not lie in the expected domain")
            }
        }
    }
}

impl std::error::Error for MatrixError {}

/// A dense matrix with elements from a [`EuclideanDomain`].
///
/// Elements are stored in row-major order.
///
/// # Example
///
/// ```
/// use ocas_domain::{EuclideanDomain, IntegerDomain, Integer};
/// use ocas_poly::matrix::Matrix;
///
/// let d = IntegerDomain;
/// let a = Matrix::from_rows(vec![
///     vec![Integer::from(1), Integer::from(1)],
///     vec![Integer::from(1), Integer::from(-1)],
/// ], d);
/// let b = vec![Integer::from(3), Integer::from(-1)];
/// let x = a.solve(&b).unwrap();
/// assert_eq!(x, vec![Integer::from(1), Integer::from(2)]);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<D: EuclideanDomain> {
    data: Vec<D::Element>,
    nrows: usize,
    ncols: usize,
    domain: D,
}

impl<D: EuclideanDomain> Matrix<D> {
    /// Create a matrix from row-major data.
    ///
    /// Panics if `data.len() != nrows * ncols`.
    pub fn new(nrows: usize, ncols: usize, data: Vec<D::Element>, domain: D) -> Self {
        assert_eq!(
            data.len(),
            nrows * ncols,
            "data length {} != {} * {}",
            data.len(),
            nrows,
            ncols
        );
        Self {
            data,
            nrows,
            ncols,
            domain,
        }
    }

    /// Create a zero matrix of the given shape.
    pub fn zeros(nrows: usize, ncols: usize, domain: D) -> Self {
        let data = vec![domain.zero(); nrows * ncols];
        Self {
            data,
            nrows,
            ncols,
            domain,
        }
    }

    /// Create an identity matrix of size `n`.
    pub fn identity(n: usize, domain: D) -> Self {
        let mut m = Self::zeros(n, n, domain);
        for i in 0..n {
            m[(i, i)] = m.domain.one();
        }
        m
    }

    /// Create a matrix from nested row vectors.
    ///
    /// Returns `None` if rows have inconsistent lengths.
    pub fn from_rows(rows: Vec<Vec<D::Element>>, domain: D) -> Self {
        let nrows = rows.len();
        let ncols = rows.first().map_or(0, |r| r.len());
        let mut data = Vec::with_capacity(nrows * ncols);
        for row in &rows {
            assert_eq!(
                row.len(),
                ncols,
                "inconsistent row lengths: expected {}, got {}",
                ncols,
                row.len()
            );
            data.extend(row.iter().cloned());
        }
        Self {
            data,
            nrows,
            ncols,
            domain,
        }
    }

    /// Return the number of rows.
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Return the number of columns.
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Return a reference to the domain.
    pub fn domain(&self) -> &D {
        &self.domain
    }

    /// Swap two rows, starting from column `start_col`.
    pub fn swap_rows(&mut self, i: usize, j: usize, start_col: usize) {
        if i == j {
            return;
        }
        for col in start_col..self.ncols {
            self.data.swap(i * self.ncols + col, j * self.ncols + col);
        }
    }

    /// Horizontally concatenate `self` with `other`.
    ///
    /// Both must have the same number of rows.
    pub fn augment(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError> {
        if self.nrows != other.nrows {
            return Err(MatrixError::ShapeMismatch);
        }
        let new_ncols = self.ncols + other.ncols;
        let mut data = Vec::with_capacity(self.nrows * new_ncols);
        for row in 0..self.nrows {
            let start = row * self.ncols;
            data.extend_from_slice(&self.data[start..start + self.ncols]);
            let ostart = row * other.ncols;
            data.extend_from_slice(&other.data[ostart..ostart + other.ncols]);
        }
        Ok(Matrix {
            data,
            nrows: self.nrows,
            ncols: new_ncols,
            domain: self.domain.clone(),
        })
    }

    /// Perform fraction-free Gaussian elimination on the first `max_col`
    /// columns, putting the matrix in row echelon form. Returns the rank.
    ///
    /// This mirrors Symbolica's `partial_row_reduce_fraction_free`.
    pub fn row_echelon(&mut self, max_col: usize) -> usize {
        let max_col = max_col.min(self.ncols);
        let mut i = 0;

        for j in 0..max_col {
            if i >= self.nrows {
                break;
            }

            if self.domain.is_zero(&self[(i, j)]) {
                // Select a non-zero pivot.
                let mut found = false;
                for k in i + 1..self.nrows {
                    if !self.domain.is_zero(&self[(k, j)]) {
                        self.swap_rows(i, k, j);
                        found = true;
                        break;
                    }
                }
                if !found {
                    continue; // zero column
                }
            }

            // Strip content from pivot row to prevent coefficient growth.
            let mut g = self[(i, j)].clone();
            for l in j + 1..self.ncols {
                if self.domain.is_one(&g) {
                    break;
                }
                g = self.domain.gcd(&g, &self[(i, l)]);
            }
            if !self.domain.is_one(&g) {
                for l in j..self.ncols {
                    self[(i, l)] = self.domain.div(&self[(i, l)], &g).unwrap();
                }
            }

            // Eliminate below.
            let pivot = self[(i, j)].clone();
            for k in i + 1..self.nrows {
                if !self.domain.is_zero(&self[(k, j)]) {
                    let g = self.domain.gcd(&pivot, &self[(k, j)]);
                    let scale_pivot = self.domain.div(&self[(k, j)], &g).unwrap();
                    let scale_row = self.domain.div(&pivot, &g).unwrap();

                    self[(k, j)] = self.domain.zero();
                    for l in j + 1..self.ncols {
                        let term1 = self.domain.mul(&self[(k, l)], &scale_row);
                        let term2 = self.domain.mul(&self[(i, l)], &scale_pivot);
                        self[(k, l)] = self.domain.sub(&term1, &term2);
                    }
                }
            }

            i += 1;
        }

        i
    }

    /// Perform fraction-free back substitution on a matrix already in
    /// row echelon form (mutating the first `max_col` columns).
    pub fn back_substitution(&mut self, max_col: usize) {
        let max_col = max_col.min(self.ncols);
        for i in (0..self.nrows).rev() {
            if let Some(j) = (0..max_col).find(|&j| !self.domain.is_zero(&self[(i, j)])) {
                // Strip content from pivot row.
                let mut g = self[(i, j)].clone();
                for l in j + 1..self.ncols {
                    if self.domain.is_one(&g) {
                        break;
                    }
                    g = self.domain.gcd(&g, &self[(i, l)]);
                }
                if !self.domain.is_one(&g) {
                    for l in j..self.ncols {
                        self[(i, l)] = self.domain.div(&self[(i, l)], &g).unwrap();
                    }
                }

                // Eliminate above.
                for k in 0..i {
                    if !self.domain.is_zero(&self[(k, j)]) {
                        let g = self.domain.gcd(&self[(i, j)], &self[(k, j)]);
                        let scale_pivot = self.domain.div(&self[(k, j)], &g).unwrap();
                        let scale_row = self.domain.div(&self[(i, j)], &g).unwrap();

                        if !self.domain.is_one(&scale_row) {
                            for l in 0..self.ncols {
                                if !self.domain.is_zero(&self[(k, l)]) {
                                    self[(k, l)] = self.domain.mul(&self[(k, l)], &scale_row);
                                }
                            }
                        }

                        self[(k, j)] = self.domain.zero();
                        for l in j + 1..self.ncols {
                            let term1 = self[(k, l)].clone();
                            let term2 = self.domain.mul(&self[(i, l)], &scale_pivot);
                            self[(k, l)] = self.domain.sub(&term1, &term2);
                        }
                    }
                }
            }
        }
    }

    /// Solve the linear system `Ax = b` where `A` is `self` and `b` is a
    /// column vector represented as a slice.
    ///
    /// Returns the solution vector `x` on success.
    pub fn solve(&self, b: &[D::Element]) -> Result<Vec<D::Element>, MatrixError> {
        if self.nrows != b.len() {
            return Err(MatrixError::ShapeMismatch);
        }

        let b_matrix = Matrix::new(b.len(), 1, b.to_vec(), self.domain.clone());

        let mut augmented = self.augment(&b_matrix)?;
        let nvars = self.ncols;

        let rank = augmented.row_echelon(nvars);

        // Check consistency.
        for k in rank..self.nrows {
            if !self.domain.is_zero(&augmented[(k, nvars)]) {
                return Err(MatrixError::Inconsistent);
            }
        }

        augmented.back_substitution(nvars);

        if rank < nvars {
            return Err(MatrixError::Underdetermined { rank });
        }

        // Divide by pivot to get the final solution.
        let mut solution = Vec::with_capacity(nvars);
        for i in 0..nvars {
            match self.domain.div(&augmented[(i, nvars)], &augmented[(i, i)]) {
                Some(val) => solution.push(val),
                None => return Err(MatrixError::ResultNotInDomain),
            }
        }

        Ok(solution)
    }

    /// Convert the matrix into a nested vector of rows.
    pub fn into_rows(self) -> Vec<Vec<D::Element>> {
        let mut rows = Vec::with_capacity(self.nrows);
        for i in 0..self.nrows {
            let start = i * self.ncols;
            rows.push(self.data[start..start + self.ncols].to_vec());
        }
        rows
    }

    /// Return a reference to the underlying row-major data.
    pub fn data(&self) -> &[D::Element] {
        &self.data
    }

    /// Return a copy of a single row as a vector.
    pub fn row(&self, i: usize) -> Vec<D::Element> {
        let start = i * self.ncols;
        self.data[start..start + self.ncols].to_vec()
    }

    /// Return a copy of a single column as a vector.
    pub fn column(&self, j: usize) -> Vec<D::Element> {
        (0..self.nrows).map(|i| self[(i, j)].clone()).collect()
    }

    /// Return the transpose of the matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// let a = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
    /// let t = a.transpose();
    /// assert_eq!(t.nrows(), 2);
    /// assert_eq!(t.ncols(), 1);
    /// assert_eq!(t[(0, 0)], Integer::from(1));
    /// assert_eq!(t[(1, 0)], Integer::from(2));
    /// ```
    pub fn transpose(&self) -> Matrix<D> {
        let mut data = Vec::with_capacity(self.nrows * self.ncols);
        for j in 0..self.ncols {
            for i in 0..self.nrows {
                data.push(self[(i, j)].clone());
            }
        }
        Matrix {
            data,
            nrows: self.ncols,
            ncols: self.nrows,
            domain: self.domain.clone(),
        }
    }

    /// Return the trace (sum of the diagonal) of a square matrix.
    ///
    /// Returns an error if the matrix is not square.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// let a = Matrix::from_rows(
    ///     vec![vec![Integer::from(1), Integer::from(2)], vec![Integer::from(3), Integer::from(4)]],
    ///     d,
    /// );
    /// assert_eq!(a.trace().unwrap(), Integer::from(5));
    /// ```
    pub fn trace(&self) -> Result<D::Element, MatrixError> {
        if self.nrows != self.ncols {
            return Err(MatrixError::ShapeMismatch);
        }
        let mut sum = self.domain.zero();
        for i in 0..self.nrows {
            sum = self.domain.add(&sum, &self[(i, i)]);
        }
        Ok(sum)
    }

    /// Compute the matrix product `self * other`.
    ///
    /// Returns an error if `self.ncols != other.nrows`.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// let a = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
    /// let b = Matrix::from_rows(vec![vec![Integer::from(3)], vec![Integer::from(4)]], d);
    /// let c = a.matmul(&b).unwrap();
    /// assert_eq!(c[(0, 0)], Integer::from(11));
    /// ```
    pub fn matmul(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError> {
        if self.ncols != other.nrows {
            return Err(MatrixError::ShapeMismatch);
        }
        let mut data = Vec::with_capacity(self.nrows * other.ncols);
        for i in 0..self.nrows {
            for j in 0..other.ncols {
                let mut acc = self.domain.zero();
                for k in 0..self.ncols {
                    let term = self.domain.mul(&self[(i, k)], &other[(k, j)]);
                    acc = self.domain.add(&acc, &term);
                }
                data.push(acc);
            }
        }
        Ok(Matrix {
            data,
            nrows: self.nrows,
            ncols: other.ncols,
            domain: self.domain.clone(),
        })
    }

    /// Compute the rank of the matrix via fraction-free Gaussian elimination.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// let a = Matrix::from_rows(
    ///     vec![vec![Integer::from(1), Integer::from(2)], vec![Integer::from(2), Integer::from(4)]],
    ///     d,
    /// );
    /// assert_eq!(a.rank(), 1);
    /// ```
    pub fn rank(&self) -> usize {
        let mut copy = self.clone();
        copy.row_echelon(self.ncols)
    }

    /// Compute the determinant of a square matrix using the Bareiss
    /// fraction-free algorithm with partial pivoting.
    ///
    /// Returns an error if the matrix is not square.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// let a = Matrix::from_rows(
    ///     vec![vec![Integer::from(1), Integer::from(2)], vec![Integer::from(3), Integer::from(4)]],
    ///     d,
    /// );
    /// assert_eq!(a.determinant().unwrap(), Integer::from(-2));
    /// ```
    pub fn determinant(&self) -> Result<D::Element, MatrixError> {
        if self.nrows != self.ncols {
            return Err(MatrixError::ShapeMismatch);
        }
        let n = self.nrows;
        if n == 0 {
            return Ok(self.domain.one());
        }
        if n == 1 {
            return Ok(self.data[0].clone());
        }
        // Bareiss fraction-free elimination with partial pivoting.
        let mut m = self.data.clone();
        let mut sign_pos = true;
        let mut prev = self.domain.one();
        for k in 0..n - 1 {
            let pivot = m[k * n + k].clone();
            if self.domain.is_zero(&pivot) {
                // Find a row below with a nonzero entry in column k.
                let mut swap_row = None;
                for i in k + 1..n {
                    if !self.domain.is_zero(&m[i * n + k]) {
                        swap_row = Some(i);
                        break;
                    }
                }
                match swap_row {
                    Some(i) => {
                        for j in 0..n {
                            m.swap(k * n + j, i * n + j);
                        }
                        sign_pos = !sign_pos;
                    }
                    None => return Ok(self.domain.zero()),
                }
            }
            let pivot = m[k * n + k].clone();
            for i in k + 1..n {
                for j in k + 1..n {
                    let term1 = self.domain.mul(&m[i * n + j], &pivot);
                    let term2 = self.domain.mul(&m[i * n + k], &m[k * n + j]);
                    let diff = self.domain.sub(&term1, &term2);
                    // Bareiss guarantees exact divisibility by prev.
                    m[i * n + j] = self.domain.div(&diff, &prev).unwrap_or(diff);
                }
            }
            prev = pivot;
        }
        let det = m[(n - 1) * n + (n - 1)].clone();
        if sign_pos {
            Ok(det)
        } else {
            Ok(self.domain.neg(&det))
        }
    }

    /// Compute the inverse of a square non-singular matrix.
    ///
    /// Returns an error if the matrix is not square or is singular over the
    /// coefficient domain.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::matrix::Matrix;
    ///
    /// let d = IntegerDomain;
    /// // Unimodular matrix: determinant 1, integer inverse exists.
    /// let a = Matrix::from_rows(
    ///     vec![vec![Integer::from(1), Integer::from(2)], vec![Integer::from(0), Integer::from(1)]],
    ///     d,
    /// );
    /// let inv = a.inverse().unwrap();
    /// assert_eq!(inv[(0, 0)], Integer::from(1));
    /// assert_eq!(inv[(0, 1)], Integer::from(-2));
    /// assert_eq!(inv[(1, 0)], Integer::from(0));
    /// assert_eq!(inv[(1, 1)], Integer::from(1));
    /// ```
    pub fn inverse(&self) -> Result<Matrix<D>, MatrixError> {
        if self.nrows != self.ncols {
            return Err(MatrixError::ShapeMismatch);
        }
        let n = self.nrows;
        let identity = Matrix::identity(n, self.domain.clone());
        // Solve A * x_j = e_j for each column j of the identity, giving
        // column j of A^{-1}. Stored in row-major order.
        let mut inv_data = vec![self.domain.zero(); n * n];
        for j in 0..n {
            let b = identity.column(j);
            let col = self.solve(&b)?;
            if col.len() != n {
                return Err(MatrixError::Underdetermined { rank: col.len() });
            }
            for i in 0..n {
                inv_data[i * n + j] = col[i].clone();
            }
        }
        Ok(Matrix {
            data: inv_data,
            nrows: n,
            ncols: n,
            domain: self.domain.clone(),
        })
    }
}

impl<D: EuclideanDomain> Index<(usize, usize)> for Matrix<D> {
    type Output = D::Element;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row * self.ncols + col]
    }
}

impl<D: EuclideanDomain> IndexMut<(usize, usize)> for Matrix<D> {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row * self.ncols + col]
    }
}

impl<D: EuclideanDomain + fmt::Display> fmt::Display for Matrix<D>
where
    D::Element: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nrows {
            if i > 0 {
                writeln!(f)?;
            }
            for j in 0..self.ncols {
                if j > 0 {
                    write!(f, " ")?;
                }
                write!(f, "{}", self[(i, j)])?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain};

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn identity_solve() {
        let d = IntegerDomain;
        let a = Matrix::identity(3, d);
        let b = vec![i(1), i(2), i(3)];
        let x = a.solve(&b).unwrap();
        assert_eq!(x, b);
    }

    #[test]
    fn solve_2x2() {
        let d = IntegerDomain;
        // 2x + y = 5
        // x + 3y = 6  => x=1, y=3... wait no: 2*1+3=5 ✓, 1+3*3=10 ≠ 6
        // Let's check: 2x+y=5, x+3y=6
        // From 1: y=5-2x. Into 2: x+3(5-2x)=6 → x+15-6x=6 → -5x=-9 → x=9/5 (not integer)
        // Let me use a solvable system: 2x+y=4, x+y=3 → x=1, y=2
        let a = Matrix::from_rows(vec![vec![i(2), i(1)], vec![i(1), i(1)]], d);
        let b = vec![i(4), i(3)];
        let x = a.solve(&b).unwrap();
        assert_eq!(x, vec![i(1), i(2)]);
    }

    #[test]
    fn solve_3x3() {
        let d = IntegerDomain;
        // x + y + z = 6
        // 2x - y + z = 3
        // x + 2y - z = 2  → x=1, y=2, z=3
        let a = Matrix::from_rows(
            vec![
                vec![i(1), i(1), i(1)],
                vec![i(2), i(-1), i(1)],
                vec![i(1), i(2), i(-1)],
            ],
            d,
        );
        let b = vec![i(6), i(3), i(2)];
        let x = a.solve(&b).unwrap();
        assert_eq!(x, vec![i(1), i(2), i(3)]);
    }

    #[test]
    fn inconsistent_system() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(1)], vec![i(1), i(1)]], d);
        let b = vec![i(1), i(2)];
        assert_eq!(a.solve(&b), Err(MatrixError::Inconsistent));
    }

    #[test]
    fn underdetermined_system() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(1), i(1)]], d);
        let b = vec![i(3)];
        assert!(matches!(
            a.solve(&b),
            Err(MatrixError::Underdetermined { .. })
        ));
    }

    #[test]
    fn row_echelon_rank() {
        let d = IntegerDomain;
        let mut a = Matrix::from_rows(
            vec![
                vec![i(1), i(2), i(3)],
                vec![i(2), i(4), i(6)],
                vec![i(0), i(1), i(1)],
            ],
            d,
        );
        let rank = a.row_echelon(3);
        assert_eq!(rank, 2);
    }

    #[test]
    fn augment_and_solve() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(3), i(2)], vec![i(1), i(1)]], d);
        let b = Matrix::new(2, 1, vec![i(5), i(2)], d);
        let aug = a.augment(&b).unwrap();
        assert_eq!(aug.nrows, 2);
        assert_eq!(aug.ncols, 3);
        assert_eq!(aug[(0, 0)], i(3));
        assert_eq!(aug[(0, 2)], i(5));
        assert_eq!(aug[(1, 0)], i(1));
        assert_eq!(aug[(1, 2)], i(2));
    }

    #[test]
    fn transpose_rect() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2), i(3)], vec![i(4), i(5), i(6)]], d);
        let t = a.transpose();
        assert_eq!(t.nrows(), 3);
        assert_eq!(t.ncols(), 2);
        assert_eq!(t[(0, 0)], i(1));
        assert_eq!(t[(0, 1)], i(4));
        assert_eq!(t[(2, 1)], i(6));
    }

    #[test]
    fn transpose_square() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(3), i(4)]], d);
        let t = a.transpose();
        assert_eq!(t[(0, 1)], i(3));
        assert_eq!(t[(1, 0)], i(2));
    }

    #[test]
    fn trace_square() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(
            vec![
                vec![i(1), i(2), i(3)],
                vec![i(4), i(5), i(6)],
                vec![i(7), i(8), i(9)],
            ],
            d,
        );
        assert_eq!(a.trace().unwrap(), i(15));
    }

    #[test]
    fn trace_nonsquare_errors() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)]], d);
        assert_eq!(a.trace(), Err(MatrixError::ShapeMismatch));
    }

    #[test]
    fn matmul_basic() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(3), i(4)]], d);
        let b = Matrix::from_rows(vec![vec![i(5), i(6)], vec![i(7), i(8)]], d);
        let c = a.matmul(&b).unwrap();
        assert_eq!(c[(0, 0)], i(19));
        assert_eq!(c[(0, 1)], i(22));
        assert_eq!(c[(1, 0)], i(43));
        assert_eq!(c[(1, 1)], i(50));
    }

    #[test]
    fn matmul_identity() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(3), i(4)]], d);
        let id = Matrix::identity(2, d);
        let c = a.matmul(&id).unwrap();
        assert_eq!(c, a);
    }

    #[test]
    fn matmul_shape_mismatch() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)]], d);
        let b = Matrix::from_rows(vec![vec![i(3), i(4)]], d);
        assert_eq!(a.matmul(&b), Err(MatrixError::ShapeMismatch));
    }

    #[test]
    fn rank_full() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(0)], vec![i(0), i(1)]], d);
        assert_eq!(a.rank(), 2);
    }

    #[test]
    fn rank_deficient() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(
            vec![
                vec![i(1), i(2), i(3)],
                vec![i(2), i(4), i(6)],
                vec![i(1), i(1), i(1)],
            ],
            d,
        );
        assert_eq!(a.rank(), 2);
    }

    #[test]
    fn determinant_2x2() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(3), i(4)]], d);
        assert_eq!(a.determinant().unwrap(), i(-2));
    }

    #[test]
    fn determinant_3x3() {
        let d = IntegerDomain;
        // det = 1*(0*7 - 6*6) - 2*(2*7 - 6*5) + 3*(2*6 - 0*5)
        //     = 1*(-36) - 2*(14-30) + 3*(12)
        //     = -36 - 2*(-16) + 36 = -36 + 32 + 36 = 32
        let a = Matrix::from_rows(
            vec![
                vec![i(1), i(2), i(3)],
                vec![i(4), i(5), i(6)],
                vec![i(5), i(6), i(7)],
            ],
            d,
        );
        // Recompute: 1*(5*7-6*6) - 2*(4*7-6*5) + 3*(4*6-5*5)
        //          = 1*(35-36) - 2*(28-30) + 3*(24-25)
        //          = -1 -2*(-2) + 3*(-1) = -1 +4 -3 = 0
        assert_eq!(a.determinant().unwrap(), i(0));
    }

    #[test]
    fn determinant_singular() {
        let d = IntegerDomain;
        // Singular: second row is 2x first.
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(2), i(4)]], d);
        assert_eq!(a.determinant().unwrap(), i(0));
    }

    #[test]
    fn determinant_nonsquare_errors() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)]], d);
        assert_eq!(a.determinant(), Err(MatrixError::ShapeMismatch));
    }

    #[test]
    fn inverse_unimodular() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(3), i(5)]], d);
        // det = 5 - 6 = -1, so integer inverse exists.
        let inv = a.inverse().unwrap();
        // A^{-1} = (1/det) * [[5,-2],[-3,1]] = (-1) * [[5,-2],[-3,1]] = [[-5,2],[3,-1]]
        assert_eq!(inv[(0, 0)], i(-5));
        assert_eq!(inv[(0, 1)], i(2));
        assert_eq!(inv[(1, 0)], i(3));
        assert_eq!(inv[(1, 1)], i(-1));
        // Verify A * A^{-1} = I.
        let prod = a.matmul(&inv).unwrap();
        assert_eq!(prod, Matrix::identity(2, IntegerDomain));
    }

    #[test]
    fn inverse_singular_errors() {
        let d = IntegerDomain;
        let a = Matrix::from_rows(vec![vec![i(1), i(2)], vec![i(2), i(4)]], d);
        assert!(a.inverse().is_err());
    }
}
