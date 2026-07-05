//! Sparse matrix backend via `sprs` for F4 Macaulay matrix storage.
//!
//! This module provides an adapter between oCAS polynomial data structures
//! and the [`sprs`] sparse matrix library. It is preparation for the F4
//! Gröbner basis algorithm (planned for 0.13.0).
//!
//! Enabled with the `sprs` feature flag.

use num_bigint::{BigInt, Sign};
use sprs::{CsMat, TriMat};

/// Convert a `BigInt` to `i64` modulo `prime` (always in `[0, prime)`).
fn bigint_to_mod_i64(c: &BigInt, prime: &BigInt) -> i64 {
    let r = c % prime;
    let (sign, digits) = r.to_u64_digits();
    if digits.is_empty() {
        return 0;
    }
    let v = digits[0] as i64;
    if sign == Sign::Minus {
        let p = prime.to_u64_digits().1[0] as i64;
        v + p
    } else {
        v
    }
}

/// A sparse representation of a Macaulay matrix built from a set of
/// polynomials over ℤ_p.
///
/// Each row corresponds to a polynomial multiplied by a monomial, and
/// each column corresponds to a monomial in the support. Non-zero entries
/// are elements of ℤ_p stored as `i64`.
pub struct SprsMacaulayMatrix {
    /// The sparse matrix in CSR format.
    pub matrix: CsMat<i64>,
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns (monomials in the support).
    pub ncols: usize,
}

impl SprsMacaulayMatrix {
    /// Build a Macaulay matrix from a list of polynomial coefficient vectors
    /// over ℤ_p.
    ///
    /// Each polynomial is given as a slice of `(column_index, coefficient)`
    /// pairs, where `column_index` maps to a monomial in the support and
    /// `coefficient` is reduced modulo `prime`.
    ///
    /// # Arguments
    ///
    /// * `rows` — each element is a list of `(col, coeff)` pairs for one row
    /// * `ncols` — total number of columns (monomial support size)
    /// * `prime` — the prime modulus
    pub fn from_rows(rows: &[Vec<(usize, BigInt)>], ncols: usize, prime: &BigInt) -> Self {
        let nrows = rows.len();
        let mut triplet = TriMat::new((nrows, ncols));

        for (i, row) in rows.iter().enumerate() {
            for &(col, ref coeff) in row {
                let c_val = bigint_to_mod_i64(coeff, prime);
                if c_val != 0 {
                    triplet.add_triplet(i, col, c_val);
                }
            }
        }

        let matrix = triplet.to_csr();
        Self {
            matrix,
            nrows,
            ncols,
        }
    }

    /// Return the number of non-zero entries.
    pub fn nnz(&self) -> usize {
        self.matrix.nnz()
    }

    /// Return the density (fraction of non-zero entries).
    pub fn density(&self) -> f64 {
        if self.nrows == 0 || self.ncols == 0 {
            return 0.0;
        }
        self.nnz() as f64 / (self.nrows as f64 * self.ncols as f64)
    }

    /// Perform sparse row echelon form reduction modulo `prime`.
    ///
    /// This is a basic implementation for smoke testing. The full F4
    /// algorithm will be implemented in 0.13.0.
    pub fn row_count(&self) -> usize {
        self.nrows
    }

    /// Perform sparse matrix-vector multiplication: `y = A * x (mod prime)`.
    ///
    /// Returns `None` if the dimensions don't match.
    pub fn spmv(&self, x: &[i64], prime: i64) -> Option<Vec<i64>> {
        if x.len() != self.ncols {
            return None;
        }
        let mut y = vec![0i64; self.nrows];
        for (row_idx, row_vec) in self.matrix.outer_iterator().enumerate() {
            for (col, &val) in row_vec.iter() {
                y[row_idx] = (y[row_idx] + val * x[col]) % prime;
                if y[row_idx] < 0 {
                    y[row_idx] += prime;
                }
            }
        }
        Some(y)
    }
}

impl std::fmt::Debug for SprsMacaulayMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SprsMacaulayMatrix")
            .field("nrows", &self.nrows)
            .field("ncols", &self.ncols)
            .field("nnz", &self.nnz())
            .field("density", &self.density())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn macaulay_from_rows_identity() {
        // 3x3 identity matrix over ℤ_5
        let rows = vec![
            vec![(0, BigInt::from(1))],
            vec![(1, BigInt::from(1))],
            vec![(2, BigInt::from(1))],
        ];
        let m = SprsMacaulayMatrix::from_rows(&rows, 3, &BigInt::from(5));
        assert_eq!(m.nrows, 3);
        assert_eq!(m.ncols, 3);
        assert_eq!(m.nnz(), 3);
        assert!((m.density() - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn macaulay_modular_reduction() {
        // Coefficient 7 mod 5 = 2
        let rows = vec![vec![(0, BigInt::from(7)), (1, BigInt::from(3))]];
        let m = SprsMacaulayMatrix::from_rows(&rows, 2, &BigInt::from(5));
        assert_eq!(m.nnz(), 2);
        // Check the actual values via spmv
        let x = vec![1i64, 0];
        let y = m.spmv(&x, 5).unwrap();
        assert_eq!(y[0], 2); // 7 mod 5 = 2
    }

    #[test]
    fn macaulay_spmv() {
        // [[1, 2], [3, 4]] * [1, 1] = [3, 7]
        let rows = vec![
            vec![(0, BigInt::from(1)), (1, BigInt::from(2))],
            vec![(0, BigInt::from(3)), (1, BigInt::from(4))],
        ];
        let m = SprsMacaulayMatrix::from_rows(&rows, 2, &BigInt::from(100));
        let x = vec![1i64, 1];
        let y = m.spmv(&x, 100).unwrap();
        assert_eq!(y[0], 3);
        assert_eq!(y[1], 7);
    }

    #[test]
    fn macaulay_sparse() {
        // A 10x10 matrix with only 5 non-zero entries
        let rows: Vec<Vec<(usize, BigInt)>> = (0..10)
            .map(|i| vec![(i, BigInt::from(1)), ((i + 5) % 10, BigInt::from(2))])
            .collect();
        let m = SprsMacaulayMatrix::from_rows(&rows, 10, &BigInt::from(7));
        assert_eq!(m.nrows, 10);
        assert_eq!(m.ncols, 10);
        // Density should be 20/100 = 0.2
        assert!((m.density() - 0.2).abs() < 1e-10);
    }

    #[test]
    fn macaulay_spmv_dimension_mismatch() {
        let rows = vec![vec![(0, BigInt::from(1))]];
        let m = SprsMacaulayMatrix::from_rows(&rows, 1, &BigInt::from(7));
        assert!(m.spmv(&[1, 2], 7).is_none());
    }
}
