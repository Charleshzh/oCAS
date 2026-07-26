//! ODE utility functions: coefficient extraction, linearity checks, order detection.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use crate::derivative::diff;

/// Return the order (highest derivative degree) of the ODE for unknown `func`.
///
/// Scans for `Derivative(func, var)`, `Derivative(func, var, var)`, etc.
/// Returns 0 if no derivatives of `func` appear.
pub(crate) fn ode_order<'a>(equation: Atom<'a>, func: Atom<'a>, var: Symbol) -> usize {
    ode_order_inner(equation, func, var)
}

fn ode_order_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> usize {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => 0,
        AtomNode::Add(args) | AtomNode::Mul(args) => args
            .iter()
            .map(|a| ode_order_inner(*a, func, var))
            .max()
            .unwrap_or(0),
        AtomNode::Pow(base, exp) => {
            ode_order_inner(*base, func, var).max(ode_order_inner(*exp, func, var))
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                // Check if the first arg is `func` and second arg is `var`.
                if args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str() {
                    // Order = number of extra `var` args beyond the first two.
                    // Derivative(y(x), x) -> order 1
                    // Derivative(y(x), x, x) -> order 2
                    args.len() - 1
                } else {
                    0
                }
            } else {
                args.iter()
                    .map(|a| ode_order_inner(*a, func, var))
                    .max()
                    .unwrap_or(0)
            }
        }
    }
}

/// Check if `expr` is linear in `func` and its derivatives.
///
/// Linear means: `func`, `Derivative(func, var)`, etc. appear only as
/// first-degree terms multiplied by expressions free of `func`.
pub(crate) fn is_linear_in<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    is_linear_inner(expr, func, var)
}

fn is_linear_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Add(args) => args.iter().all(|a| is_linear_inner(*a, func, var)),
        AtomNode::Mul(args) => {
            // At most one factor may contain `func`.
            let func_dependent_count = args
                .iter()
                .filter(|a| contains_func(**a, func, var))
                .count();
            if func_dependent_count > 1 {
                return false;
            }
            // The factor containing func must itself be linear (first degree).
            args.iter().all(|a| {
                if contains_func(*a, func, var) {
                    is_func_first_degree(*a, func, var)
                } else {
                    true
                }
            })
        }
        AtomNode::Pow(base, exp) => {
            // If base contains func, exponent must be a positive constant == 1.
            if contains_func(*base, func, var) {
                if let AtomNode::Num(n) = exp.node() {
                    *n == 1
                } else {
                    false
                }
            } else {
                // exp may contain func only if it's first degree and base is constant.
                !contains_func(*exp, func, var)
                    || (is_func_first_degree(*exp, func, var) && !contains_func(*base, func, var))
            }
        }
        AtomNode::Fun(_, args) => {
            // A function call containing func is nonlinear (e.g. sin(y)).
            !contains_func(expr, func, var)
                || args.iter().all(|a| {
                    if contains_func(*a, func, var) {
                        is_func_first_degree(*a, func, var)
                    } else {
                        true
                    }
                })
        }
    }
}

/// Check if `expr` contains `func` or its derivatives.
pub(crate) fn contains_func<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    contains_func_inner(expr, func, var)
}

fn contains_func_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_func_inner(*a, func, var))
        }
        AtomNode::Pow(base, exp) => {
            contains_func_inner(*base, func, var) || contains_func_inner(*exp, func, var)
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str()
            } else {
                args.iter().any(|a| contains_func_inner(*a, func, var))
            }
        }
    }
}

/// Check if `func` or its derivatives appear only to the first power.
fn is_func_first_degree<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str()
            } else {
                false
            }
        }
        AtomNode::Pow(base, exp) => {
            if contains_func(*base, func, var) {
                if let AtomNode::Num(n) = exp.node() {
                    *n == 1
                } else {
                    false
                }
            } else {
                !contains_func(*exp, func, var)
            }
        }
        AtomNode::Mul(args) => {
            let func_count: usize = args
                .iter()
                .filter(|a| contains_func(**a, func, var))
                .count();
            func_count <= 1
        }
        AtomNode::Add(args) => args.iter().all(|a| is_func_first_degree(*a, func, var)),
        _ => !contains_func(expr, func, var),
    }
}

/// Substitute `y(x) → sol_expr` in `equation`, also replacing
/// `Derivative(y(x), x)` with `diff(sol_expr, x)`, and higher derivatives
/// accordingly.
pub(crate) fn substitute_solution<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a> {
    substitute_inner(ctx, equation, func, sol, var)
}

fn substitute_inner<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_inner(ctx, *a, func, sol, var))
                .collect();
            ctx.add(&mapped)
        }
        AtomNode::Mul(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_inner(ctx, *a, func, sol, var))
                .collect();
            ctx.mul(&mapped)
        }
        AtomNode::Pow(base, exp) => {
            let b = substitute_inner(ctx, *base, func, sol, var);
            let e = substitute_inner(ctx, *exp, func, sol, var);
            ctx.pow(b, e)
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                if args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str() {
                    // Derivative(y(x), x, x, ...) -> differentiate sol n times.
                    let order = args.len() - 1;
                    let mut result = sol;
                    for _ in 0..order {
                        result = diff(ctx, result, var);
                    }
                    result
                } else {
                    let mapped: Vec<_> = args
                        .iter()
                        .map(|a| substitute_inner(ctx, *a, func, sol, var))
                        .collect();
                    ctx.fun(name.as_str(), &mapped)
                }
            } else if args.iter().any(|a| contains_func_inner(*a, func, var)) {
                // func appears inside a regular function call (e.g. sin(y))
                let mapped: Vec<_> = args
                    .iter()
                    .map(|a| substitute_inner(ctx, *a, func, sol, var))
                    .collect();
                ctx.fun(name.as_str(), &mapped)
            } else {
                expr
            }
        }
    }
}
