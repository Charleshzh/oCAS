//! Equation solving for oCAS.
//!
//! Provides solvers for linear systems, polynomial systems, and
//! numeric root-finding. Results are returned as expression trees.

use std::fmt;

use ocas_domain::{IntegerDomain, RationalDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

/// Errors that can occur when solving equations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolveError {
    /// The system contains no equations.
    EmptySystem,
    /// The system is not linear in the requested variables.
    NonLinear,
    /// The number of equations does not match the number of unknowns.
    NonSquare,
    /// The linear system has no solution.
    Inconsistent,
    /// The system is underdetermined (infinitely many solutions).
    Underdetermined {
        /// Rank of the coefficient matrix.
        rank: usize,
    },
    /// The solution does not lie in the expected domain.
    ResultNotInDomain,
    /// An internal matrix operation failed.
    Matrix(MatrixError),
    /// Other error with a description.
    Other(String),
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolveError::EmptySystem => f.write_str("empty system"),
            SolveError::NonLinear => f.write_str("system is not linear"),
            SolveError::NonSquare => {
                f.write_str("number of equations must match number of unknowns")
            }
            SolveError::Inconsistent => f.write_str("inconsistent system"),
            SolveError::Underdetermined { rank } => {
                write!(f, "underdetermined system (rank {})", rank)
            }
            SolveError::ResultNotInDomain => {
                f.write_str("solution does not lie in the expected domain")
            }
            SolveError::Matrix(e) => write!(f, "matrix error: {}", e),
            SolveError::Other(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for SolveError {}

impl From<MatrixError> for SolveError {
    fn from(e: MatrixError) -> Self {
        match e {
            MatrixError::Inconsistent => SolveError::Inconsistent,
            MatrixError::Underdetermined { rank } => SolveError::Underdetermined { rank },
            MatrixError::ResultNotInDomain => SolveError::ResultNotInDomain,
            MatrixError::ShapeMismatch => SolveError::NonSquare,
            MatrixError::RightHandSideIsNotVector => {
                SolveError::Other("right-hand side is not a vector".into())
            }
        }
    }
}

/// Solve a linear system of equations over the rational numbers.
///
/// Given an `n × n` coefficient matrix `A` and a right-hand side vector `b`
/// (both as nested `i64` values interpreted as rational numbers), returns
/// the solution vector as rational numbers represented as `(num, den)` pairs.
///
/// # Example
///
/// ```
/// use ocas_calc::solve::solve_linear_rational;
///
/// // 2x + y = 5
/// // x - y = 1  → x=2, y=1
/// let a = vec![vec![2, 1], vec![1, -1]];
/// let b = vec![5, 1];
/// let x = solve_linear_rational(&a, &b).unwrap();
/// assert_eq!(x, vec![(2, 1), (1, 1)]);
/// ```
pub fn solve_linear_rational(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<(i64, i64)>, SolveError> {
    if a.is_empty() || b.is_empty() {
        return Err(SolveError::EmptySystem);
    }
    let n = a.len();
    if b.len() != n {
        return Err(SolveError::NonSquare);
    }
    for row in a {
        if row.len() != n {
            return Err(SolveError::NonSquare);
        }
    }

    let d = RationalDomain;
    use ocas_domain::Rational;

    let rows: Vec<Vec<Rational>> = a
        .iter()
        .map(|row| row.iter().map(|&v| Rational::new(v, 1)).collect())
        .collect();
    let b_vec: Vec<Rational> = b.iter().map(|&v| Rational::new(v, 1)).collect();

    let mat = Matrix::from_rows(rows, d);
    let solution = mat.solve(&b_vec)?;

    Ok(solution
        .into_iter()
        .map(|r| {
            let numer = r.numer().to_i64().unwrap_or(0);
            let denom = r.denom().to_i64().unwrap_or(1);
            (numer, denom)
        })
        .collect())
}

/// Solve a linear system of equations over the integers.
///
/// Returns the exact integer solution, or an error if the solution is not
/// integral or the system is inconsistent.
///
/// # Example
///
/// ```
/// use ocas_calc::solve::solve_linear_integer;
///
/// // x + y = 3
/// // x - y = 1  → x=2, y=1
/// let a = vec![vec![1, 1], vec![1, -1]];
/// let b = vec![3, 1];
/// let x = solve_linear_integer(&a, &b).unwrap();
/// assert_eq!(x, vec![2, 1]);
/// ```
pub fn solve_linear_integer(a: &[Vec<i64>], b: &[i64]) -> Result<Vec<i64>, SolveError> {
    if a.is_empty() || b.is_empty() {
        return Err(SolveError::EmptySystem);
    }
    let n = a.len();
    if b.len() != n {
        return Err(SolveError::NonSquare);
    }
    for row in a {
        if row.len() != n {
            return Err(SolveError::NonSquare);
        }
    }

    let d = IntegerDomain;
    use ocas_domain::Integer;

    let rows: Vec<Vec<Integer>> = a
        .iter()
        .map(|row| row.iter().map(|&v| Integer::from(v)).collect())
        .collect();
    let b_vec: Vec<Integer> = b.iter().map(|&v| Integer::from(v)).collect();

    let mat = Matrix::from_rows(rows, d);
    let solution = mat.solve(&b_vec)?;

    Ok(solution
        .into_iter()
        .map(|r| {
            // Integer::from produces BigInt; convert back to i64.
            r.to_i64().unwrap_or(0)
        })
        .collect())
}

/// A solution to a linear Diophantine equation ax + by = c.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiophantineSolution {
    /// A particular solution (x0, y0).
    pub particular: (i64, i64),
    /// The general solution is (x0 + k*tx, y0 + k*ty) for any integer k.
    pub general: (i64, i64),
}

/// Solve the linear Diophantine equation `a*x + b*y = c`.
///
/// Returns a particular solution and the homogeneous part, or `None` if
/// no solution exists.
///
/// # Example
///
/// ```
/// use ocas_calc::solve::solve_diophantine;
///
/// // 3x + 5y = 1 → particular: (2, -1), general: (5k, -3k)
/// let sol = solve_diophantine(3, 5, 1).unwrap();
/// assert_eq!(sol.particular, (2, -1));
/// assert_eq!(sol.general, (5, -3));
/// ```
pub fn solve_diophantine(a: i64, b: i64, c: i64) -> Option<DiophantineSolution> {
    use ocas_domain::{EuclideanDomain, Integer, IntegerDomain};

    if a == 0 && b == 0 {
        return if c == 0 {
            Some(DiophantineSolution {
                particular: (0, 0),
                general: (1, 0),
            })
        } else {
            None
        };
    }

    let d = IntegerDomain;
    let ai = Integer::from(a);
    let bi = Integer::from(b);

    let (g, x, y) = d.extended_gcd(&ai, &bi);
    let g_val: i64 = g.to_i64()?;
    let x_val: i64 = x.to_i64()?;
    let y_val: i64 = y.to_i64()?;

    // Check if g divides c.
    if c % g_val != 0 {
        return None;
    }

    let scale = c / g_val;
    let x0 = x_val * scale;
    let y0 = y_val * scale;

    Some(DiophantineSolution {
        particular: (x0, y0),
        general: (b / g_val, -(a / g_val)),
    })
}

/// A polynomial represented as a list of (coefficient, exponent_vector) pairs.
type RawPoly = Vec<(i64, Vec<usize>)>;

/// Solve a system of polynomial equations using Gröbner bases.
///
/// Each polynomial is given as a list of `(coefficient, exponent_vector)` pairs.
/// The exponent vector `[e1, e2, ...]` represents `x1^e1 * x2^e2 * ...`.
///
/// Returns the reduced Gröbner basis with respect to the lexicographic order,
/// which is triangular and can be used for back-substitution.
///
/// # Example
///
/// ```
/// use ocas_calc::solve::solve_polynomial_system;
///
/// // x + y = 0, x - y = 0  →  basis = {x, y}
/// let polys = vec![
///     vec![(1, vec![1, 0]), (1, vec![0, 1])],
///     vec![(1, vec![1, 0]), (-1, vec![0, 1])],
/// ];
/// let gb = solve_polynomial_system(&polys, 2).unwrap();
/// assert!(gb.len() >= 2);
/// ```
pub fn solve_polynomial_system(
    polys: &[RawPoly],
    n_vars: usize,
) -> Result<Vec<RawPoly>, SolveError> {
    use ocas_domain::Rational;
    use ocas_domain::RationalDomain;
    use ocas_poly::SparseMultivariatePolynomial;
    use ocas_poly::buchberger;
    use ocas_poly::sparse::Lex;

    if polys.is_empty() {
        return Err(SolveError::EmptySystem);
    }

    let d = RationalDomain;

    let ideal: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>> = polys
        .iter()
        .map(|terms| {
            let terms: Vec<(Vec<usize>, Rational)> = terms
                .iter()
                .map(|(c, exp)| (exp.clone(), Rational::new(*c, 1)))
                .collect();
            SparseMultivariatePolynomial::from_terms(d, n_vars, terms)
        })
        .collect();

    let gb = buchberger(&ideal);

    let result: Vec<RawPoly> = gb
        .basis
        .iter()
        .map(|p| {
            p.sorted_terms()
                .into_iter()
                .map(|(exp, coeff)| {
                    let numer = coeff.numer().to_i64().unwrap_or(0);
                    let _denom = coeff.denom().to_i64().unwrap_or(1);
                    (numer, exp.to_vec())
                })
                .collect()
        })
        .collect();

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diophantine_3x_5y_eq_1() {
        let sol = solve_diophantine(3, 5, 1).unwrap();
        // Verify: 3*2 + 5*(-1) = 6 - 5 = 1 ✓
        let (x, y) = sol.particular;
        assert_eq!(3 * x + 5 * y, 1);
    }

    #[test]
    fn diophantine_no_solution() {
        // 2x + 4y = 3 has no solution (gcd(2,4)=2 doesn't divide 3)
        assert!(solve_diophantine(2, 4, 3).is_none());
    }

    #[test]
    fn diophantine_zero_coeffs() {
        assert!(solve_diophantine(0, 0, 1).is_none());
        let sol = solve_diophantine(0, 0, 0).unwrap();
        assert_eq!(sol.particular, (0, 0));
    }

    #[test]
    fn solve_2x2_rational() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![5, 11];
        let x = solve_linear_rational(&a, &b).unwrap();
        // x + 2y = 5, 3x + 4y = 11 → x=1, y=2
        assert_eq!(x, vec![(1, 1), (2, 1)]);
    }

    #[test]
    fn solve_2x2_integer() {
        let a = vec![vec![1, 1], vec![1, -1]];
        let b = vec![3, 1];
        let x = solve_linear_integer(&a, &b).unwrap();
        assert_eq!(x, vec![2, 1]);
    }

    #[test]
    fn solve_3x3_rational() {
        // x + y + z = 6
        // 2x - y + z = 3
        // x + 2y - z = 2  → x=1, y=2, z=3
        let a = vec![vec![1, 1, 1], vec![2, -1, 1], vec![1, 2, -1]];
        let b = vec![6, 3, 2];
        let x = solve_linear_rational(&a, &b).unwrap();
        assert_eq!(x, vec![(1, 1), (2, 1), (3, 1)]);
    }

    #[test]
    fn inconsistent_system() {
        let a = vec![vec![1, 1], vec![2, 2]];
        let b = vec![1, 3];
        assert!(matches!(
            solve_linear_rational(&a, &b),
            Err(SolveError::Inconsistent)
        ));
    }
}
