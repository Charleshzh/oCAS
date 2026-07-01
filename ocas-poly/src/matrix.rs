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
}
