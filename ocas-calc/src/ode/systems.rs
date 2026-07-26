//! Linear ODE system solvers.
//!
//! Provides [`solve_linear_system`] for constant-coefficient systems of the
//! form $\mathbf{Y}' = A\mathbf{Y} + \mathbf{g}(x)$.

use ocas_atom::{Atom, AtomArena, Symbol};

use super::ODESolution;

/// Solve a constant-coefficient linear ODE system.
///
/// Given a system $\mathbf{Y}' = A\mathbf{Y} + \mathbf{g}(x)$ where $A$ is a
/// constant matrix and $\mathbf{g}(x)$ is a forcing vector, attempts to find
/// an analytical solution.
///
/// # Parameters
///
/// - `equations`: list of ODE equations in the form `y_i' - (A*Y + g_i) = 0`
/// - `funcs`: list of unknown functions `[y1(x), y2(x), ...]`
/// - `var`: the independent variable
///
/// # Limitations
///
/// - Only constant-coefficient systems are supported.
/// - The matrix $A$ must be diagonalizable (eigendecomposition approach).
/// - Non-diagonalizable systems fall back to `Unsolved`.
#[allow(dead_code)]
pub(crate) fn solve_linear_system<'a>(
    _ctx: &'a AtomArena<'a>,
    _equations: &[Atom<'a>],
    _funcs: &[Atom<'a>],
    _var: Symbol,
) -> Option<ODESolution<'a>> {
    // Full linear system solving requires:
    // 1. Extract the coefficient matrix A from the system
    // 2. Compute eigenvalues (polynomial root finding)
    // 3. Compute eigenvectors (null space computation)
    // 4. Build matrix exponential e^{Ax}
    // 5. Compute particular solution via variation of parameters
    //
    // This is a significant piece of infrastructure that requires
    // integration with ocas-poly (root finding) and ocas-atom::matrix.
    // For now, return None to indicate the system cannot be solved.
    //
    // Future implementation:
    // - 2x2 systems: direct characteristic polynomial approach
    // - NxN systems: eigendecomposition via companion matrix
    None
}
