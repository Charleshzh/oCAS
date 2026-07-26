//! ODE classification engine.
//!
//! Analyzes an ODE equation and determines which solving methods may apply.
//! The classifier returns a list of [`ODEType`] candidates in priority order.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::ODE;
use super::util::{contains_func, is_linear_in, ode_order};

/// The type of an ODE, determining which solver to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ODEType {
    /// Separable: $f(x)\,g(y)\,dy = h(x)\,k(y)\,dx$.
    Separable,
    /// First-order linear: $y' + p(x)\,y = q(x)$.
    LinearFirst,
    /// Bernoulli: $y' + p(x)\,y = q(x)\,y^n$.
    Bernoulli,
    /// Exact: $M(x,y) + N(x,y)\,y' = 0$ with $\partial M/\partial y = \partial N/\partial x$.
    Exact,
    /// Homogeneous: $y' = f(y/x)$.
    Homogeneous,
    /// Second-order (or higher) linear with constant coefficients.
    LinearConstantCoeff,
    /// Cauchy-Euler: $a\,x^2\,y'' + b\,x\,y' + c\,y = f(x)$.
    CauchyEuler,
    /// Power series solution around an ordinary point.
    PowerSeries,
}

/// Classify an ODE and return all applicable solving methods in priority order.
///
/// The priority order is chosen to try the most specific methods first:
/// 1. Separable
/// 2. First-order linear
/// 3. Bernoulli
/// 4. Exact
/// 5. Homogeneous
/// 6. Second-order constant coefficient
/// 7. Cauchy-Euler
/// 8. Power series (fallback for linear ODEs)
pub fn classify_ode<'a>(ctx: &'a AtomArena<'a>, ode: ODE<'a>) -> Vec<ODEType> {
    let ODE {
        equation,
        func,
        var,
    } = ode;

    let order = ode_order(equation, func, var);
    let mut types = Vec::new();

    if order == 0 {
        // Not actually a differential equation.
        return types;
    }

    // Normalize equation form for analysis: collect terms on one side.
    // The equation is already in the form expr = 0.

    if order == 1 {
        // --- First-order classification ---
        let linear = is_first_order_linear(ctx, equation, func, var);

        if linear {
            types.push(ODEType::LinearFirst);
        }

        // Check Bernoulli: y' + p(x)*y = q(x)*y^n for n != 0, 1.
        if is_bernoulli(ctx, equation, func, var) {
            types.push(ODEType::Bernoulli);
        }

        // Check separable: the equation can be written as g(y)*y' = f(x).
        if is_separable(ctx, equation, func, var) {
            types.push(ODEType::Separable);
        }

        // Check exact: M + N*y' = 0 with dM/dy = dN/dx.
        if is_exact(ctx, equation, func, var) {
            types.push(ODEType::Exact);
        }

        // Check homogeneous: y' = F(y/x).
        if is_homogeneous(ctx, equation, func, var) {
            types.push(ODEType::Homogeneous);
        }
    } else {
        // --- Second-order (and higher) classification ---
        if is_linear_in(equation, func, var) {
            if is_constant_coeff_linear(ctx, equation, func, var) {
                types.push(ODEType::LinearConstantCoeff);
            }
            if is_cauchy_euler(ctx, equation, func, var) {
                types.push(ODEType::CauchyEuler);
            }
            // If linear but neither constant-coeff nor Cauchy-Euler,
            // we can still try constant_coeff as a fallback (it will fail gracefully).
            if types.is_empty() {
                types.push(ODEType::LinearConstantCoeff);
            }
        }
    }

    // Power series is always a fallback for linear ODEs.
    if is_linear_in(equation, func, var) && ode_order(equation, func, var) >= 1 {
        types.push(ODEType::PowerSeries);
    }

    types
}

/// Check if the ODE is first-order linear: y' + p(x)*y = q(x).
///
/// This is true when:
/// - The equation has order 1
/// - y and y' appear linearly (degree 1)
/// - No nonlinear terms in y (no y^2, sin(y), etc.)
fn is_first_order_linear<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    // The equation must be linear in func and its first derivative.
    if !is_linear_in(equation, func, var) {
        return false;
    }

    // Must contain y (otherwise it's just 0 = f(x)).
    if !contains_func(equation, func, var) {
        return false;
    }

    // Order must be exactly 1.
    ode_order(equation, func, var) == 1
}

/// Check if the ODE is Bernoulli: y' + p(x)*y = q(x)*y^n.
///
/// We detect this by checking:
/// - Linear in y' (order 1)
/// - Contains a nonlinear y^n term (n != 0, 1) alongside linear y terms.
fn is_bernoulli<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    if ode_order(equation, func, var) != 1 {
        return false;
    }

    // For Bernoulli detection, we look for a pattern where the equation
    // has both linear y terms and a y^n term. This is a heuristic check.
    has_nonlinear_func_term(equation, func, var)
}

/// Check if a nonlinear func term (func^n for n >= 2) appears in expr.
fn has_nonlinear_func_term<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) => args.iter().any(|a| has_nonlinear_func_term(*a, func, var)),
        AtomNode::Mul(args) => args.iter().any(|a| has_nonlinear_func_term(*a, func, var)),
        AtomNode::Pow(base, exp) => {
            if contains_func(*base, func, var) {
                if let AtomNode::Num(n) = exp.node() {
                    *n >= 2
                } else {
                    true // Non-constant exponent with func base is nonlinear.
                }
            } else {
                has_nonlinear_func_term(*exp, func, var)
            }
        }
        AtomNode::Fun(_, args) => args.iter().any(|a| has_nonlinear_func_term(*a, func, var)),
    }
}

/// Check if the ODE is separable.
///
/// A first-order ODE is separable if it can be written as g(y)*y' = f(x),
/// meaning all terms with y can be grouped with y', and all other terms
/// are functions of x only.
fn is_separable<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    if ode_order(equation, func, var) != 1 {
        return false;
    }

    // Get the y' symbol representation.
    let _dy = derivative_atom(func, var);

    // Try to express the equation as A(x,y)*dy + B(x,y) = 0.
    // For separability: A must factor as g(y)*a(x) and B must factor as f(x)*b(y).
    // This is a heuristic: we check if the equation is a product of a function
    // of x and a function involving y and y'.

    // Simple heuristic: check if the equation is already in the form
    // h(y)*y' - f(x) = 0 or equivalently M(x)*N(y) + P(x)*Q(y)*y' = 0.
    // We do a basic check: are there additive terms that are free of func?
    let has_free_terms = has_func_free_additive_terms(equation, func, var);
    let has_func_terms = contains_func(equation, func, var);

    // A separable ODE typically has both kinds of terms.
    has_free_terms && has_func_terms && ode_order(equation, func, var) == 1
}

/// Check if the ODE is exact: M(x,y) + N(x,y)*y' = 0 with dM/dy = dN/dx.
fn is_exact<'a>(_ctx: &'a AtomArena<'a>, equation: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    if ode_order(equation, func, var) != 1 {
        return false;
    }

    // For an exact ODE, we need to identify M and N from the equation
    // written as M + N*y' = 0. Then check dM/dy = dN/dx.
    // This is a deeper check that requires coefficient extraction.
    // We do a simplified version here.

    // For now, any first-order linear ODE is also exact if the
    // integrating factor condition holds, but we use a simpler heuristic:
    // Check if the equation has the structure where M and N can be identified.
    let y_sym = func_symbol(func);
    if y_sym.is_none() {
        return false;
    }

    // Simplified check: if the equation is linear, it's already covered by
    // LinearFirst. Exact is mainly for nonlinear first-order equations.
    if is_linear_in(equation, func, var) {
        return false;
    }

    // For nonlinear first-order, we attempt to identify M and N.
    // This is a placeholder: full exactness checking requires symbolic
    // partial differentiation with respect to the dependent variable.
    false
}

/// Check if the ODE is homogeneous: y' = F(y/x).
///
/// For a first-order ODE in the form dy/dx = f(x,y), it's homogeneous
/// if f(tx, ty) = f(x,y) (degree 0 homogeneity).
fn is_homogeneous<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    if ode_order(equation, func, var) != 1 {
        return false;
    }
    if !is_linear_in(equation, func, var) {
        return false;
    }

    // Check if all terms have the same total degree in (x, y).
    // For linear first-order ODEs, this means every term is degree 1.
    // E.g., x*y' - y = 0 is homogeneous (each term is degree 1 in x,y).
    // y' + y/x = 0 is also homogeneous.

    // Heuristic: check if the equation, when written in standard form
    // y' + p(x)*y = q(x), has q(x) = 0 (homogeneous linear case)
    // or if p(x) and q(x) are such that the equation is scale-invariant.

    // For now, check if the equation is linear and the "forcing" term (free of y)
    // is zero or a homogeneous function.
    let has_free_terms = has_func_free_additive_terms(equation, func, var);

    // True homogeneous: all terms have the same degree.
    // Quick check: does every additive term involve func or func' to degree 1?
    all_terms_homogeneous_degree(equation, func, var) && has_free_terms
}

/// Check if all additive terms in expr have the same total degree in (func, x).
fn all_terms_homogeneous_degree<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    // Collect the "degree" of each additive term.
    // A term's degree = total power of x + total power of y (func).
    match expr.node() {
        AtomNode::Add(args) => {
            let degrees: Vec<i64> = args.iter().map(|a| term_degree(*a, func, var)).collect();
            if degrees.is_empty() {
                return true;
            }
            let first = degrees[0];
            first > 0 && degrees.iter().all(|&d| d == first)
        }
        _ => {
            let d = term_degree(expr, func, var);
            d > 0
        }
    }
}

/// Compute the total degree of a multiplicative term in x and func.
fn term_degree<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> i64 {
    match expr.node() {
        AtomNode::Num(_) => 0,
        AtomNode::Var(v) => {
            if *v == var {
                1
            } else {
                0
            }
        }
        AtomNode::Mul(args) => args.iter().map(|a| term_degree(*a, func, var)).sum(),
        AtomNode::Pow(base, exp) => {
            if let AtomNode::Num(n) = exp.node() {
                term_degree(*base, func, var) * *n
            } else {
                1 // conservative
            }
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                // Derivative of func counts as degree 1 in func.
                if args[0].to_string() == func.to_string() {
                    1
                } else {
                    0
                }
            } else {
                // func itself (y(x))
                if expr.to_string() == func.to_string() {
                    1
                } else {
                    0
                }
            }
        }
        AtomNode::Add(_) => 1, // treat sub-expressions as degree 1 (conservative)
    }
}

/// Check if expr has additive terms that are free of func.
fn has_func_free_additive_terms<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Add(args) => args.iter().any(|a| !contains_func(*a, func, var)),
        _ => !contains_func(expr, func, var),
    }
}

/// Check if a second-order (or higher) linear ODE has constant coefficients.
fn is_constant_coeff_linear<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    if !is_linear_in(equation, func, var) {
        return false;
    }
    // Check that coefficients of y, y', y'' etc. are free of var.
    // Heuristic: collect all terms, check that the factor multiplying
    // each Derivative(func, var, ...) or func itself is free of var.
    coefficients_free_of_var(equation, func, var)
}

/// Check if the coefficient of each func/func' term is free of the independent variable.
fn coefficients_free_of_var<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Add(args) => args.iter().all(|a| coefficient_free_of_var(*a, func, var)),
        AtomNode::Mul(args) => {
            // One factor should be func/derivative, rest should be free of var.
            let func_factor_count = args
                .iter()
                .filter(|a| is_func_or_derivative(**a, func, var))
                .count();
            if func_factor_count == 1 {
                args.iter().all(|a| {
                    if is_func_or_derivative(*a, func, var) {
                        true
                    } else {
                        !contains_var(*a, var)
                    }
                })
            } else if func_factor_count == 0 {
                // No func factor — this is a constant term, fine.
                true
            } else {
                false
            }
        }
        // Bare func or derivative (implicit coefficient 1).
        AtomNode::Fun(name, args) if *name == Symbol::new("Derivative") && args.len() >= 2 => {
            args[0].to_string() == func.to_string()
        }
        _ => true,
    }
}

fn coefficient_free_of_var<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Mul(args) => {
            let func_factor = args.iter().find(|a| is_func_or_derivative(**a, func, var));
            if func_factor.is_some() {
                args.iter().all(|a| {
                    if is_func_or_derivative(*a, func, var) {
                        true
                    } else {
                        !contains_var(*a, var)
                    }
                })
            } else {
                // No func factor — pure coefficient term, must be free of var.
                !contains_var(expr, var)
            }
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string()
            } else if expr.to_string() == func.to_string() {
                true
            } else {
                !contains_var(expr, var)
            }
        }
        _ => !contains_var(expr, var),
    }
}

/// Check if the ODE is Cauchy-Euler: a*x^2*y'' + b*x*y' + c*y = f(x).
fn is_cauchy_euler<'a>(
    _ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> bool {
    if !is_linear_in(equation, func, var) {
        return false;
    }
    let order = ode_order(equation, func, var);
    if order < 2 {
        return false;
    }

    // Check the pattern: coefficient of y^(k) is a_k * x^k.
    is_cauchy_euler_pattern(equation, func, var)
}

/// Heuristic check for Cauchy-Euler pattern.
fn is_cauchy_euler_pattern<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Add(args) => args.iter().all(|a| is_ce_term(*a, func, var)),
        _ => is_ce_term(expr, func, var),
    }
}

/// Check if a single term matches the Cauchy-Euler pattern: c_k * x^k * y^(k).
fn is_ce_term<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Mul(args) => {
            // Find the func/derivative factor.
            let func_factor = args.iter().find(|a| is_func_or_derivative(**a, func, var));
            if let Some(ff) = func_factor {
                let order = func_derivative_order(*ff, func, var);
                // The remaining factors should form c_k * x^order.
                let coeff_factors: Vec<_> = args
                    .iter()
                    .filter(|a| !is_func_or_derivative(**a, func, var))
                    .copied()
                    .collect();
                match coeff_factors.len() {
                    0 => order == 0, // bare y with no coefficient
                    1 => {
                        // Should be c * x^order or just c (for order 0).
                        if order == 0 {
                            !contains_var(coeff_factors[0], var)
                        } else {
                            is_x_power(coeff_factors[0], var, order as i64)
                        }
                    }
                    _ => {
                        // Multiple coefficient factors: one should be x^order,
                        // the rest constants.
                        let has_x_power = coeff_factors
                            .iter()
                            .any(|c| is_x_power(*c, var, order as i64));
                        let rest_const = coeff_factors
                            .iter()
                            .filter(|c| !is_x_power(**c, var, order as i64))
                            .all(|c| !contains_var(*c, var));
                        has_x_power && rest_const
                    }
                }
            } else {
                // Pure coefficient term (no func), must be free of var.
                !contains_var(expr, var)
            }
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string()
            } else {
                expr.to_string() == func.to_string()
            }
        }
        _ => !contains_var(expr, var),
    }
}

/// Check if expr is x^n for a specific n.
fn is_x_power<'a>(expr: Atom<'a>, var: Symbol, n: i64) -> bool {
    match expr.node() {
        AtomNode::Var(v) => *v == var && n == 1,
        AtomNode::Pow(base, exp) => {
            if let AtomNode::Num(e) = exp.node() {
                if let AtomNode::Var(v) = base.node() {
                    *v == var && *e == n
                } else {
                    false
                }
            } else {
                false
            }
        }
        AtomNode::Num(1) => n == 0,
        _ => false,
    }
}

/// Check if expr is func or a derivative of func.
fn is_func_or_derivative<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str()
            } else {
                expr.to_string() == func.to_string()
            }
        }
        _ => expr.to_string() == func.to_string(),
    }
}

/// Get the derivative order of a func/derivative expression. 0 for bare func.
fn func_derivative_order<'a>(expr: Atom<'a>, _func: Atom<'a>, _var: Symbol) -> usize {
    match expr.node() {
        AtomNode::Fun(name, args) if *name == Symbol::new("Derivative") && args.len() >= 2 => {
            args.len() - 1
        }
        _ => 0,
    }
}

/// Check if expr contains the variable `var`.
fn contains_var<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => *v == var,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_var(*a, var)),
        AtomNode::Pow(base, exp) => contains_var(*base, var) || contains_var(*exp, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_var(*a, var)),
    }
}

/// Helper: build a `Derivative(func, var)` Atom from parts.
fn derivative_atom<'a>(func: Atom<'a>, var: Symbol) -> String {
    format!("Derivative({}, {})", func, var.as_str())
}

/// Extract the function name symbol from a `y(x)` Atom.
fn func_symbol<'a>(func: Atom<'a>) -> Option<Symbol> {
    match func.node() {
        AtomNode::Fun(name, _) => Some(*name),
        AtomNode::Var(v) => Some(*v),
        _ => None,
    }
}
