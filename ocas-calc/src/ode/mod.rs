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
            ODEType::PowerSeries => series::solve_power_series(ctx, ode, ctx.num(0), 8),
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
    // Differentiate sol as needed and substitute.
    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let result = simplify(ctx, substituted, &rules, 20);
    let result = normalize(ctx, result);
    matches!(result.node(), AtomNode::Num(0))
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
}
