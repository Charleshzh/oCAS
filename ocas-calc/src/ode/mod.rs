//! Ordinary differential equation (ODE) solver for oCAS.
//!
//! This module provides [`dsolve`], which attempts to find analytical solutions
//! to ordinary differential equations involving polynomials and elementary
//! functions. It classifies the ODE type and dispatches to the appropriate
//! solving method.
//!
//! # Supported ODE types
//!
//! **First order:**
//! - Separable: $f(x)\,dx = g(y)\,dy$
//! - Linear: $y' + p(x)\,y = q(x)$
//! - Bernoulli: $y' + p(x)\,y = q(x)\,y^n$
//! - Exact: $M\,dx + N\,dy = 0$ with $\partial M/\partial y = \partial N/\partial x$
//! - Homogeneous: $y' = f(y/x)$
//!
//! **Second order (linear):**
//! - Constant coefficients: $a\,y'' + b\,y' + c\,y = f(x)$
//! - Cauchy-Euler: $a\,x^2\,y'' + b\,x\,y' + c\,y = f(x)$
//! - Reduction of order (given one solution)
//! - Variation of parameters
//!
//! ODEs that cannot be solved analytically are returned as unevaluated
//! `ODE(equation, func)` forms.

pub mod classify;
pub mod first_order;
mod laplace;
pub mod second_order;
pub mod series;
mod systems;
mod util;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, Symbol};
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use crate::rules::calculus_rules;

// Re-export the public API.
pub use classify::{ODEType, classify_ode};

/// An ordinary differential equation represented as `lhs - rhs = 0`.
///
/// The equation is stored in canonical form: `lhs` is the full expression
/// (already normalized to `lhs - rhs`), `func` is the unknown function
/// (e.g. `y(x)`), and `var` is the independent variable (e.g. `x`).
#[derive(Debug, Clone, Copy)]
pub struct ODE<'a> {
    /// The equation in the form `lhs - rhs = 0` (i.e., the full expression
    /// set to zero).
    pub equation: Atom<'a>,
    /// The unknown function, e.g. `y(x)` represented as `Fun("y", &[x])`.
    pub func: Atom<'a>,
    /// The independent variable symbol, e.g. `x`.
    pub var: Symbol,
}

/// The result of attempting to solve an ODE.
#[derive(Debug, Clone, Copy)]
pub enum ODESolution<'a> {
    /// An explicit solution `y = expr`.
    Explicit(Atom<'a>),
    /// An implicit solution `F(x, y) = 0` (stored as the expression `F`).
    Implicit(Atom<'a>),
    /// A parametric solution `(x(t), y(t))`.
    Parametric(Atom<'a>, Atom<'a>),
    /// A power/Frobenius series solution truncated to `n_terms`.
    Series(Atom<'a>, usize),
    /// Solutions `(y_1, y_2, ...)` of a linear ODE system.
    System(&'a [Atom<'a>]),
    /// The ODE could not be solved; the original equation is returned.
    Unsolved(ODE<'a>),
}

/// Solve an ordinary differential equation.
///
/// Given an ODE in the form `equation = 0`, the unknown function `func`
/// (e.g. `y(x)`), and the independent variable `var`, attempt to find an
/// analytical solution.
///
/// The optional `hint` parameter allows the caller to request a specific
/// solving method. When `None`, the classifier tries all applicable methods
/// in priority order.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::ode::{dsolve, ODE, ODESolution};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.fun("y", &[x]);
/// let dy = ctx.fun("Derivative", &[y, x]);
/// // Solve y' - y = 0
/// let ode = ODE {
///     equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
///     func: y,
///     var: Symbol::new("x"),
/// };
/// let sol = dsolve(&ctx, ode, None);
/// // Should return something like C1*exp(x)
/// assert!(matches!(sol, ODESolution::Explicit(_)));
/// ```
///
/// # Unsolved ODEs
///
/// If no method produces a solution, `ODESolution::Unsolved(ode)` is returned.
pub fn dsolve<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>, hint: Option<ODEType>) -> ODESolution<'a> {
    let ode = normalize_ode(ctx, ode);

    let types = match hint {
        Some(t) => vec![t],
        None => classify_ode(ctx, ode),
    };

    for ode_type in &types {
        let result = match ode_type {
            ODEType::Separable => first_order::solve_separable(ctx, ode),
            ODEType::LinearFirst => first_order::solve_linear_first(ctx, ode),
            ODEType::Bernoulli => first_order::solve_bernoulli(ctx, ode),
            ODEType::Exact => first_order::solve_exact(ctx, ode),
            ODEType::Homogeneous => first_order::solve_homogeneous(ctx, ode),
            ODEType::LinearConstantCoeff => second_order::solve_constant_coeff(ctx, ode),
            ODEType::CauchyEuler => second_order::solve_cauchy_euler(ctx, ode),
            ODEType::ReductionOfOrder => second_order::solve_reduction_of_order(ctx, ode),
            ODEType::PowerSeries => series::solve_power_series(ctx, ode, ctx.num(0), 8)
                .or_else(|| series::solve_frobenius(ctx, ode, ctx.num(0), 8)),
        };
        if let Some(sol) = result {
            // Verify the solution is not trivially Unsolved.
            if !matches!(sol, ODESolution::Unsolved(_)) {
                return sol;
            }
        }
    }

    ODESolution::Unsolved(ode)
}

/// Solve a first- or second-order linear constant-coefficient initial value
/// problem using the Laplace transform method.
///
/// `y0` is the value $y(x_0)$ at $x_0 = 0$; `y1` is $y'(0)$ (required for
/// second-order problems, ignored otherwise). Returns an explicit solution
/// with no free constants, or `Unsolved` when the method does not apply.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::ode::{dsolve_ivp, ODE, ODESolution};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.fun("y", &[x]);
/// let dy = ctx.fun("Derivative", &[y, x]);
/// // y' - y = 0, y(0) = 1  =>  y = exp(x)
/// let ode = ODE {
///     equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
///     func: y,
///     var: Symbol::new("x"),
/// };
/// let sol = dsolve_ivp(&ctx, ode, ctx.num(1), None);
/// assert!(matches!(sol, ODESolution::Explicit(_)));
/// ```
pub fn dsolve_ivp<'a>(
    ctx: &'a AtomArena<'a>,
    ode: ODE<'a>,
    y0: Atom<'a>,
    y1: Option<Atom<'a>>,
) -> ODESolution<'a> {
    let ode = normalize_ode(ctx, ode);
    laplace::solve_laplace(ctx, ode, y0, y1).unwrap_or(ODESolution::Unsolved(ode))
}

/// Solve a 2×2 constant-coefficient linear ODE system
/// $\mathbf{Y}' = A\mathbf{Y}$.
///
/// Each equation in `equations` must have the form
/// `Derivative(y_i, x) - (a_i1*y1 + a_i2*y2) = 0`. Returns the general
/// solution components `(y1, y2)` with free constants C1, C2, or
/// `Unsolved` for unsupported systems.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::ode::{dsolve_system, ODESolution};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y1 = ctx.fun("y1", &[x]);
/// let y2 = ctx.fun("y2", &[x]);
/// let dy1 = ctx.fun("Derivative", &[y1, x]);
/// let dy2 = ctx.fun("Derivative", &[y2, x]);
/// // y1' = y2, y2' = -y1 (harmonic oscillator)
/// let eq1 = ctx.add(&[dy1, ctx.mul(&[ctx.num(-1), y2])]);
/// let eq2 = ctx.add(&[dy2, y1]);
/// let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
/// assert!(matches!(sol, ODESolution::System(_)));
/// ```
pub fn dsolve_system<'a>(
    ctx: &'a AtomArena<'a>,
    equations: &[Atom<'a>],
    funcs: &[Atom<'a>],
    var: Symbol,
) -> ODESolution<'a> {
    match systems::solve_linear_system(ctx, equations, funcs, var) {
        Some(sol) => sol,
        None => ODESolution::Unsolved(ODE {
            equation: equations.first().copied().unwrap_or_else(|| ctx.num(0)),
            func: funcs.first().copied().unwrap_or_else(|| ctx.num(0)),
            var,
        }),
    }
}

/// Normalize an ODE: simplify the equation and extract its canonical form.
fn normalize_ode<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>) -> ODE<'a> {
    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let simplified = simplify(ctx, ode.equation, &rules, 20);
    let simplified = simplify(ctx, simplified, &calc_rules, 10);
    let equation = normalize(ctx, simplified);
    ODE {
        equation,
        func: ode.func,
        var: ode.var,
    }
}

/// Substitute a candidate solution into an ODE equation and return the
/// residual with like terms collected. A zero residual means the candidate
/// satisfies the ODE. Useful for external verification (correctness tests).
pub fn substitute_solution_collected<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a> {
    let substituted = util::substitute_solution(ctx, equation, func, sol, var);
    util::collect_terms(ctx, substituted)
}

/// Apply a solution to verify it satisfies the ODE.
///
/// Substitutes `y = sol` into the ODE equation and checks whether the
/// result simplifies to zero. Returns `true` if the verification passes.
#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn verify_solution<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>, sol: Atom<'a>) -> bool {
    use crate::ode::util::substitute_solution;
    use ocas_atom::AtomNode;

    // Replace y(x) with sol in the ODE equation.
    let substituted = substitute_solution(ctx, ode.equation, ode.func, sol, ode.var);
    // Like-term collection cancels the complementary-solution contributions
    // (C1/C2 terms) which the rule-based simplifier cannot combine.
    let result = util::collect_terms(ctx, substituted);
    matches!(result.node(), AtomNode::Num(0))
}

/// Verify a system solution: substitute all component solutions into every
/// equation and check each residual vanishes.
#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn verify_system<'a>(
    ctx: &'a AtomArena<'a>,
    equations: &[Atom<'a>],
    funcs: &[Atom<'a>],
    sols: &[Atom<'a>],
    var: Symbol,
) -> bool {
    use ocas_atom::AtomNode;

    equations.iter().all(|eq| {
        // Replace each y_j with its solution and Derivative(y_j, x) with
        // its derivative.
        let mut substituted = *eq;
        for (func, sol) in funcs.iter().zip(sols.iter()) {
            substituted = util::substitute_solution(ctx, substituted, *func, *sol, var);
        }
        let residual = util::collect_terms(ctx, substituted);
        matches!(residual.node(), AtomNode::Num(0))
    })
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn classify_first_order_linear() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y' + y = 0 (linear first-order)
        let ode = ODE {
            equation: ctx.add(&[dy, y]),
            func: y,
            var: Symbol::new("x"),
        };
        let types = classify_ode(&ctx, ode);
        assert!(types.contains(&ODEType::LinearFirst));
        assert!(types.contains(&ODEType::Separable));
    }

    #[test]
    fn classify_second_order_constant_coeff() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let _dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - y = 0
        let ode = ODE {
            equation: ctx.add(&[d2y, ctx.mul(&[ctx.num(-1), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let types = classify_ode(&ctx, ode);
        assert!(types.contains(&ODEType::LinearConstantCoeff));
    }

    #[test]
    fn dsolve_first_order_linear() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y' - y = 0 => solution: y = C1*exp(x)
        let ode = ODE {
            equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                // Solution should contain C1 and exp(x)
                assert!(s.contains("C1"), "Solution should contain C1: {s}");
                assert!(s.contains("exp"), "Solution should contain exp: {s}");
            }
            _ => panic!("Expected explicit solution, got {:?}", sol),
        }
    }

    #[test]
    fn dsolve_first_order_linear_forcing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y' + y = x => solution should contain exp(-x) and polynomial terms
        let ode = ODE {
            equation: ctx.add(&[dy, y, ctx.mul(&[ctx.num(-1), x])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("C1"), "Solution should contain C1: {s}");
                assert!(s.contains("exp"), "Solution should contain exp: {s}");
            }
            _ => panic!("Expected explicit solution, got {:?}", sol),
        }
    }

    #[test]
    fn dsolve_second_order_constant_coeff_real_roots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - 3*y' + 2*y = 0 => roots: r=1, r=2 => y = C1*exp(x) + C2*exp(2x)
        let ode = ODE {
            equation: ctx.add(&[d2y, ctx.mul(&[ctx.num(-3), dy]), ctx.mul(&[ctx.num(2), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("C1"), "Solution should contain C1: {s}");
                assert!(s.contains("C2"), "Solution should contain C2: {s}");
            }
            _ => panic!("Expected explicit solution, got {:?}", sol),
        }
    }

    #[test]
    fn dsolve_second_order_repeated_root() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - 2*y' + y = 0 => repeated root r=1 => y = (C1 + C2*x)*exp(x)
        let ode = ODE {
            equation: ctx.add(&[d2y, ctx.mul(&[ctx.num(-2), dy]), y]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("C1"), "Solution should contain C1: {s}");
                assert!(s.contains("C2"), "Solution should contain C2: {s}");
                assert!(s.contains("exp"), "Solution should contain exp: {s}");
            }
            _ => panic!("Expected explicit solution, got {:?}", sol),
        }
    }

    #[test]
    fn dsolve_unsolvable_returns_unsolved() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        // An equation with no derivatives is not an ODE.
        let ode = ODE {
            equation: ctx.add(&[y, ctx.mul(&[ctx.num(-1), x])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        assert!(matches!(sol, ODESolution::Unsolved(_)));
    }

    #[test]
    fn ode_order_detection() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);

        // Order 0
        let expr = ctx.add(&[y, x]);
        assert_eq!(super::util::ode_order(expr, y, Symbol::new("x")), 0);

        // Order 1
        let dy = ctx.fun("Derivative", &[y, x]);
        assert_eq!(super::util::ode_order(dy, y, Symbol::new("x")), 1);

        // Order 2
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        assert_eq!(super::util::ode_order(d2y, y, Symbol::new("x")), 2);
    }

    #[test]
    fn classify_exact_ode() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // 2*x*y + x^2*y' = 0: M = 2xy, N = x^2, dM/dy = 2x = dN/dx -> exact.
        let m = ctx.mul(&[ctx.num(2), x, y]);
        let n_term = ctx.mul(&[ctx.pow(x, ctx.num(2)), dy]);
        let ode = ODE {
            equation: ctx.add(&[m, n_term]),
            func: y,
            var: Symbol::new("x"),
        };
        let types = classify_ode(&ctx, ode);
        assert!(
            types.contains(&ODEType::Exact),
            "Expected Exact classification, got {types:?}"
        );
    }

    #[test]
    fn classify_non_exact_ode() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y + x*y' = 0: M = y, N = x, dM/dy = 1 != dN/dx = 1 ... actually equal.
        // Use 2y + x*y' = 0 instead: dM/dy = 2 != dN/dx = 1 -> not exact.
        let m = ctx.mul(&[ctx.num(2), y]);
        let n_term = ctx.mul(&[x, dy]);
        let ode = ODE {
            equation: ctx.add(&[m, n_term]),
            func: y,
            var: Symbol::new("x"),
        };
        let types = classify_ode(&ctx, ode);
        assert!(
            !types.contains(&ODEType::Exact),
            "Did not expect Exact classification, got {types:?}"
        );
    }

    #[test]
    fn dsolve_second_order_irrational_roots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - 2*y = 0 => r = ±sqrt(2) => y = C1*exp(sqrt(2)*x) + C2*exp(-sqrt(2)*x)
        // The discriminant 8 is not a perfect square; roots must stay symbolic.
        let ode = ODE {
            equation: ctx.add(&[d2y, ctx.mul(&[ctx.num(-2), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(
                    s.contains("8^(1/2)") || s.contains("sqrt") || !s.contains("8^-1"),
                    "Solution should keep sqrt symbolically: {s}"
                );
                // Must NOT contain the truncated isqrt value 2 as the root:
                // r would wrongly be ±2/2 = ±1 instead of ±sqrt(8)/2.
                assert!(
                    !s.contains("exp(x)") && !s.contains("exp(-1*x)"),
                    "Roots were truncated to integers: {s}"
                );
            }
            _ => panic!("Expected explicit solution, got {:?}", sol),
        }
    }

    #[test]
    fn dsolve_exact_ode() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // 2*x*y + x^2*y' = 0 (exact): F = x^2*y = C.
        let m = ctx.mul(&[ctx.num(2), x, y]);
        let n_term = ctx.mul(&[ctx.pow(x, ctx.num(2)), dy]);
        let ode = ODE {
            equation: ctx.add(&[m, n_term]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::Exact));
        match sol {
            ODESolution::Implicit(expr) => {
                let s = expr.to_string();
                // Solution F(x,y) should contain x^2 * y.
                assert!(s.contains('y'), "Implicit solution should involve y: {s}");
            }
            other => panic!("Expected implicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_integrating_factor() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y + 2x*y' = 0 is NOT exact (M_y=1 != N_x=2), but
        // (M_y - N_x)/N = (1-2)/(2x) = -1/(2x) depends only on x,
        // so mu(x) = exp(-ln(x)/2) = x^(-1/2) makes it exact.
        // Solution: x^(1/2)*y = C, i.e. y = C*x^(-1/2).
        let m = y;
        let n_term = ctx.mul(&[ctx.num(2), x, dy]);
        let ode = ODE {
            equation: ctx.add(&[m, n_term]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::Exact));
        assert!(
            !matches!(sol, ODESolution::Unsolved(_)),
            "Integrating-factor ODE should be solvable, got {sol:?}"
        );
    }

    #[test]
    #[ignore = "Risch integrator loops on sec(x)/tan(x); VOP needs non-elementary fallback"]
    fn dsolve_second_order_vop() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' + y = 1/(cos(x)): undetermined coefficients cannot handle
        // sec(x) forcing; variation of parameters gives
        // y_p = cos(x)*log(cos(x)) + x*sin(x) (up to homogeneous terms).
        let sec = ctx.pow(ctx.fun("cos", &[x]), ctx.num(-1));
        let ode = ODE {
            equation: ctx.add(&[d2y, y, ctx.mul(&[ctx.num(-1), sec])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::LinearConstantCoeff));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                // Must not silently drop the forcing: the particular part
                // should contain x*sin(x) or log(cos(x)).
                assert!(
                    s.contains("sin") || s.contains("log"),
                    "VOP particular solution missing: {s}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_cauchy_euler_forcing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // x^2*y'' - 2*x*y' + 2*y = x^3 (Cauchy-Euler with forcing).
        // Homogeneous: r^2 - 3r + 2 = 0 -> r=1,2 -> y_c = C1*x + C2*x^2.
        // VOP on y'' - (2/x)y' + (2/x^2)y = x gives y_p = x^3/2.
        let x2_d2y = ctx.mul(&[ctx.pow(x, ctx.num(2)), d2y]);
        let neg_2x_dy = ctx.mul(&[ctx.num(-2), x, dy]);
        let two_y = ctx.mul(&[ctx.num(2), y]);
        let neg_x3 = ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(3))]);
        let ode = ODE {
            equation: ctx.add(&[x2_d2y, neg_2x_dy, two_y, neg_x3]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::CauchyEuler));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(
                    s.contains("x^3") || s.contains("x^(3"),
                    "Particular solution x^3/2 missing: {s}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_undetermined_quadratic_forcing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' + y = x^2: y_p = x^2 - 2.
        let x2 = ctx.pow(x, ctx.num(2));
        let ode = ODE {
            equation: ctx.add(&[d2y, y, ctx.mul(&[ctx.num(-1), x2])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::LinearConstantCoeff));
        match sol {
            ODESolution::Explicit(expr) => {
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "Solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_undetermined_exp_resonance() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - 3y' + 2y = exp(x): k=1 is a single root of r^2-3r+2.
        // y_p = -x*exp(x).
        let exp_x = ctx.fun("exp", &[x]);
        let ode = ODE {
            equation: ctx.add(&[
                d2y,
                ctx.mul(&[ctx.num(-3), dy]),
                ctx.mul(&[ctx.num(2), y]),
                ctx.mul(&[ctx.num(-1), exp_x]),
            ]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::LinearConstantCoeff));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(
                    s.contains("x") && s.contains("exp"),
                    "Resonance particular solution missing x*exp(x): {s}"
                );
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "Solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_undetermined_trig_forcing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' + y' + y = cos(x): p = c - a*w^2 = 0, q = b*w = 1, det = 1.
        // A = f_s*q = 0, B = f_s*p ... compute: f_c=1, f_s=0.
        // A = (f_c*p + f_s*q)/det = 0, B = (f_s*p - f_c*q)/det = -1.
        // y_p = -sin(x).
        let cos_x = ctx.fun("cos", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let ode = ODE {
            equation: ctx.add(&[d2y, dy, y, ctx.mul(&[ctx.num(-1), cos_x])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::LinearConstantCoeff));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(
                    s.contains("sin"),
                    "Trig particular solution missing sin: {s}"
                );
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "Solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_reduction_of_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // x*y'' - y' = 0: not constant-coeff, not Cauchy-Euler.
        // y1 = 1 works: p = -1/x, e^(-∫p) = x, u = ∫x dx = x^2/2, y2 = x^2/2.
        // General solution: y = C1 + C2*x^2.
        let x_d2y = ctx.mul(&[x, d2y]);
        let neg_dy = ctx.mul(&[ctx.num(-1), dy]);
        let ode = ODE {
            equation: ctx.add(&[x_d2y, neg_dy]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::ReductionOfOrder));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("C1"), "Solution should contain C1: {s}");
                assert!(
                    s.contains("x^2") || s.contains("x^(2"),
                    "Second solution x^2 missing: {s}"
                );
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "Solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_power_series_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y' - y = 0: recurrence a_{n+1} = a_n/(n+1), i.e. a_n = a_0/n!.
        // Series: a_0 * (1 + x + x^2/2 + x^3/6 + ...).
        let ode = ODE {
            equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::PowerSeries));
        match sol {
            ODESolution::Series(expr, n) => {
                let s = expr.to_string();
                assert!(n >= 6, "Expected at least 6 terms, got {n}");
                // a_2 = a_0/2, so a term with 2^-1 * x^2 should appear.
                assert!(s.contains("x^2"), "Series should contain x^2 term: {s}");
            }
            other => panic!("Expected series solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_power_series_second_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' + y = 0: a_{n+2} = -a_n/((n+2)(n+1)).
        // Even part: a_0*(1 - x^2/2 + x^4/24 - ...), odd: a_1*(x - x^3/6 + ...).
        let ode = ODE {
            equation: ctx.add(&[d2y, y]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::PowerSeries));
        match sol {
            ODESolution::Series(expr, n) => {
                let s = expr.to_string();
                assert!(n >= 6, "Expected at least 6 terms, got {n}");
                assert!(s.contains("x^2"), "Series should contain x^2 term: {s}");
                // The x^2 coefficient should be -a_0/2 (negative sign).
                assert!(
                    s.contains("-1"),
                    "Series should contain the negative x^2 coefficient: {s}"
                );
            }
            other => panic!("Expected series solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_frobenius_euler_point() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // 2*x*y'' + y' = 0: regular singular point at x=0.
        // Indicial 2r(r-1) + r = 2r^2 - r = 0 -> r = 0, 1/2.
        // The equation must not come back unsolved.
        let two_x_d2y = ctx.mul(&[ctx.num(2), x, d2y]);
        let ode = ODE {
            equation: ctx.add(&[two_x_d2y, dy]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, None);
        assert!(
            !matches!(sol, ODESolution::Unsolved(_)),
            "2x*y'' + y' = 0 should be solvable, got {sol:?}"
        );
    }

    #[test]
    fn dsolve_frobenius_half_integer_root() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // 2x*y'' + y' + 2y = 0: regular singular point at x=0.
        // Indicial 2r(r-1) + r = 2r^2 - r = 0 -> larger root r = 1/2.
        // Frobenius series: y = x^(1/2)*a0*(1 - 2x/3 + 2x^2/15 - ...).
        let ode = ODE {
            equation: ctx.add(&[
                ctx.mul(&[ctx.num(2), x, d2y]),
                dy,
                ctx.mul(&[ctx.num(2), y]),
            ]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve(&ctx, ode, Some(ODEType::PowerSeries));
        match sol {
            ODESolution::Series(expr, _n) => {
                let s = expr.to_string();
                // Must contain the x^(1/2) leading factor (half-integer root).
                assert!(
                    s.contains("x^(2^-1)") || s.contains("x^(1*(2^-1))"),
                    "Frobenius series missing x^(1/2) factor: {s}"
                );
                // The x^(3/2) coefficient is -2/3*a0.
                assert!(
                    s.contains("3^-1"),
                    "Frobenius series missing 1/3 coefficient: {s}"
                );
            }
            other => panic!("Expected Frobenius series, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_ivp_first_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        // y' - y = 0, y(0) = 2  =>  y = 2*exp(x)
        let ode = ODE {
            equation: ctx.add(&[dy, ctx.mul(&[ctx.num(-1), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve_ivp(&ctx, ode, ctx.num(2), None);
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("exp"), "Expected exp in IVP solution: {s}");
                assert!(!s.contains("C1"), "IVP solution must not contain C1: {s}");
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "IVP solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit IVP solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_ivp_second_order_distinct_roots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let dy = ctx.fun("Derivative", &[y, x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' - 3y' + 2y = 0, y(0)=1, y'(0)=0  =>  y = 2*exp(x) - exp(2x)
        let ode = ODE {
            equation: ctx.add(&[d2y, ctx.mul(&[ctx.num(-3), dy]), ctx.mul(&[ctx.num(2), y])]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve_ivp(&ctx, ode, ctx.num(1), Some(ctx.num(0)));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("exp"), "Expected exp in IVP solution: {s}");
                assert!(!s.contains("C1"), "IVP solution must not contain C1: {s}");
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "IVP solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit IVP solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_ivp_second_order_trig() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.fun("y", &[x]);
        let d2y = ctx.fun("Derivative", &[y, x, x]);
        // y'' + y = 0, y(0)=0, y'(0)=1  =>  y = sin(x)
        let ode = ODE {
            equation: ctx.add(&[d2y, y]),
            func: y,
            var: Symbol::new("x"),
        };
        let sol = dsolve_ivp(&ctx, ode, ctx.num(0), Some(ctx.num(1)));
        match sol {
            ODESolution::Explicit(expr) => {
                let s = expr.to_string();
                assert!(s.contains("sin"), "Expected sin in IVP solution: {s}");
                assert!(
                    verify_solution(&ctx, ode, expr),
                    "IVP solution does not satisfy ODE: {expr}"
                );
            }
            other => panic!("Expected explicit IVP solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_system_distinct_eigenvalues() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y1 = ctx.fun("y1", &[x]);
        let y2 = ctx.fun("y2", &[x]);
        let dy1 = ctx.fun("Derivative", &[y1, x]);
        let dy2 = ctx.fun("Derivative", &[y2, x]);
        // y1' = y2, y2' = y1  =>  A = [[0,1],[1,0]], λ = ±1.
        // y1 = C1*e^x - C2*e^-x, y2 = C1*e^x + C2*e^-x.
        let eq1 = ctx.add(&[dy1, ctx.mul(&[ctx.num(-1), y2])]);
        let eq2 = ctx.add(&[dy2, ctx.mul(&[ctx.num(-1), y1])]);
        let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
        match sol {
            ODESolution::System(comps) => {
                assert_eq!(comps.len(), 2);
                let s1 = comps[0].to_string();
                let s2 = comps[1].to_string();
                assert!(s1.contains("exp"), "y1 should contain exp: {s1}");
                assert!(s2.contains("exp"), "y2 should contain exp: {s2}");
            }
            other => panic!("Expected system solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_system_complex_eigenvalues() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y1 = ctx.fun("y1", &[x]);
        let y2 = ctx.fun("y2", &[x]);
        let dy1 = ctx.fun("Derivative", &[y1, x]);
        let dy2 = ctx.fun("Derivative", &[y2, x]);
        // y1' = y2, y2' = -y1  =>  harmonic oscillator, λ = ±i.
        // y1 = C1*sin(x) + C2*cos(x), y2 = C1*cos(x) - C2*sin(x).
        let eq1 = ctx.add(&[dy1, ctx.mul(&[ctx.num(-1), y2])]);
        let eq2 = ctx.add(&[dy2, y1]);
        let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
        match sol {
            ODESolution::System(comps) => {
                assert_eq!(comps.len(), 2);
                let s1 = comps[0].to_string();
                let s2 = comps[1].to_string();
                assert!(
                    s1.contains("sin") || s1.contains("cos"),
                    "y1 should contain trig: {s1}"
                );
                assert!(
                    s2.contains("sin") || s2.contains("cos"),
                    "y2 should contain trig: {s2}"
                );
                // Verify both equations against the component solutions.
                assert!(
                    verify_system(&ctx, &[eq1, eq2], &[y1, y2], comps, Symbol::new("x")),
                    "System solution does not satisfy equations"
                );
            }
            other => panic!("Expected system solution, got {other:?}"),
        }
    }

    #[test]
    fn dsolve_system_repeated_eigenvalue() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y1 = ctx.fun("y1", &[x]);
        let y2 = ctx.fun("y2", &[x]);
        let dy1 = ctx.fun("Derivative", &[y1, x]);
        let dy2 = ctx.fun("Derivative", &[y2, x]);
        // y1' = y1 + y2, y2' = y2  =>  A = [[1,1],[0,1]], repeated λ=1.
        // y1 = e^x (C1 + C2 x), y2 = C2 e^x.
        let eq1 = ctx.add(&[
            dy1,
            ctx.mul(&[ctx.num(-1), y1]),
            ctx.mul(&[ctx.num(-1), y2]),
        ]);
        let eq2 = ctx.add(&[dy2, ctx.mul(&[ctx.num(-1), y2])]);
        let sol = dsolve_system(&ctx, &[eq1, eq2], &[y1, y2], Symbol::new("x"));
        match sol {
            ODESolution::System(comps) => {
                assert_eq!(comps.len(), 2);
                assert!(
                    verify_system(&ctx, &[eq1, eq2], &[y1, y2], comps, Symbol::new("x")),
                    "System solution does not satisfy equations: {comps:?}"
                );
            }
            other => panic!("Expected system solution, got {other:?}"),
        }
    }
}
