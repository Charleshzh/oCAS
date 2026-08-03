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
    /// Reduction of order for second-order linear ODEs (tries simple
    /// candidate solutions, then builds the second solution).
    ReductionOfOrder,
    /// Power series solution around an ordinary point.
    PowerSeries,
}

/// Classify an ODE and return all applicable solving methods in priority order.
///
/// The priority order is chosen to try the most specific methods first.
/// First-order ODEs are checked in this order (matching the push order in
/// `classify_ode`):
/// 1. First-order linear
/// 2. Bernoulli
/// 3. Separable
/// 4. Exact
/// 5. Homogeneous
/// Second-order ODEs additionally check:
/// 6. Constant coefficient
/// 7. Cauchy-Euler
/// 8. Reduction of order
/// Finally, every linear ODE of order ≥ 1 gets the power-series fallback.
///
/// # Example
///
/// ```no_run
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_core::arena::Arena;
/// use ocas_calc::ode::{classify_ode, ODE, ODEType};
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.fun("y", &[x]);
/// let eq = ctx.add(&[y, ctx.mul(&[ctx.num(-1), x])]);
/// let ode = ODE { equation: eq, func: ctx.fun("y", &[x]), var: Symbol::new("x") };
/// let methods = classify_ode(&ctx, ode);
/// assert!(methods.contains(&ODEType::LinearFirst));
/// ```
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
            // Reduction of order works for any second-order linear ODE
            // when a simple candidate solution can be found.
            if ode_order(equation, func, var) == 2 {
                types.push(ODEType::ReductionOfOrder);
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
fn is_exact<'a>(ctx: &'a AtomArena<'a>, equation: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    use ocas_rewrite::rules::default_rules;
    use ocas_rewrite::simplify::simplify;

    if ode_order(equation, func, var) != 1 {
        return false;
    }

    let Some(y_sym) = func_symbol(func) else {
        return false;
    };

    // Split equation into M (terms without y') and N (coefficient of y').
    // y' must appear linearly; terms containing y' inside powers or
    // nonlinear functions disqualify the exact form.
    let Some((m, n)) = split_mn(ctx, equation, func, var) else {
        return false;
    };

    // Replace y(x) with a plain symbol y so that M and N are bivariate
    // expressions in x and y; then compare partial derivatives.
    let y_var = ctx.var(y_sym.as_str());
    let m_sub = replace_atom(ctx, m, func, y_var);
    let n_sub = replace_atom(ctx, n, func, y_var);

    let dm_dy = crate::derivative::diff(ctx, m_sub, y_sym);
    let dn_dx = crate::derivative::diff(ctx, n_sub, var);

    // Exactness: dM/dy == dN/dx. Compare normalized forms structurally
    // (the simplifier does not combine like terms such as 2x - 2x, so a
    // difference-is-zero check alone is too weak).
    let dm_norm = ocas_atom::normalize::normalize(ctx, dm_dy);
    let dn_norm = ocas_atom::normalize::normalize(ctx, dn_dx);
    if dm_norm.to_string() == dn_norm.to_string() {
        return true;
    }

    // Fallback: try the difference after rule-based simplification.
    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let difference = simplify(
        ctx,
        ctx.add(&[dm_dy, ctx.mul(&[ctx.num(-1), dn_dx])]),
        &rules,
        20,
    );
    let difference = ocas_atom::normalize::normalize(ctx, difference);
    matches!(difference.node(), AtomNode::Num(0))
}

/// Split a first-order equation `M + N*y' = 0` into (M, N).
///
/// Returns `None` when y' does not appear linearly (e.g. inside a power or
/// another function), or when either side is empty.
pub(crate) fn split_mn<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let dy_str = derivative_atom(func, var);

    let mut m_terms = Vec::new();
    let mut n_terms = Vec::new();

    let terms: Vec<Atom<'a>> = match equation.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![equation],
    };

    for term in terms {
        let s = term.to_string();
        if s == dy_str {
            // Bare y' with coefficient 1.
            n_terms.push(ctx.num(1));
            continue;
        }
        match term.node() {
            AtomNode::Mul(args) => {
                let dy_count = args.iter().filter(|a| a.to_string() == dy_str).count();
                if dy_count == 1 {
                    // Coefficient of y' is the product of the other factors.
                    let rest: Vec<_> = args
                        .iter()
                        .filter(|a| a.to_string() != dy_str)
                        .copied()
                        .collect();
                    // The remaining factors must not contain y' nonlinearly.
                    if rest.iter().any(|a| contains_derivative_str(*a, &dy_str)) {
                        return None;
                    }
                    n_terms.push(if rest.is_empty() {
                        ctx.num(1)
                    } else {
                        ctx.mul(&rest)
                    });
                } else if dy_count == 0 {
                    if contains_derivative_str(term, &dy_str) {
                        // y' inside a power or function — nonlinear.
                        return None;
                    }
                    m_terms.push(term);
                } else {
                    // (y')^2 or higher — nonlinear.
                    return None;
                }
            }
            _ => {
                if contains_derivative_str(term, &dy_str) {
                    // y' inside Pow/Fun — nonlinear.
                    return None;
                }
                m_terms.push(term);
            }
        }
    }

    if m_terms.is_empty() || n_terms.is_empty() {
        return None;
    }

    let m = if m_terms.len() == 1 {
        m_terms[0]
    } else {
        ctx.add(&m_terms)
    };
    let n = if n_terms.len() == 1 {
        n_terms[0]
    } else {
        ctx.add(&n_terms)
    };
    Some((m, n))
}

/// Check whether `expr` contains the derivative string anywhere.
pub(crate) fn contains_derivative_str<'a>(expr: Atom<'a>, dy_str: &str) -> bool {
    if expr.to_string() == dy_str {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_derivative_str(*a, dy_str))
        }
        AtomNode::Pow(base, exp) => {
            contains_derivative_str(*base, dy_str) || contains_derivative_str(*exp, dy_str)
        }
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_derivative_str(*a, dy_str)),
    }
}

/// Replace every occurrence of `target` inside `expr` with `replacement`.
pub(crate) fn replace_atom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    target: Atom<'a>,
    replacement: Atom<'a>,
) -> Atom<'a> {
    if expr.to_string() == target.to_string() {
        return replacement;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| replace_atom(ctx, *a, target, replacement))
                .collect();
            ctx.add(&mapped)
        }
        AtomNode::Mul(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| replace_atom(ctx, *a, target, replacement))
                .collect();
            ctx.mul(&mapped)
        }
        AtomNode::Pow(base, exp) => {
            let b = replace_atom(ctx, *base, target, replacement);
            let e = replace_atom(ctx, *exp, target, replacement);
            ctx.pow(b, e)
        }
        AtomNode::Fun(name, args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| replace_atom(ctx, *a, target, replacement))
                .collect();
            ctx.fun(name.as_str(), &mapped)
        }
    }
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
