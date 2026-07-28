//! First-order ODE solvers.
//!
//! Implements five classical methods for solving first-order ODEs:
//! - Separable equations
//! - First-order linear equations
//! - Bernoulli equations
//! - Exact equations
//! - Homogeneous equations

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use super::ODESolution;
use super::util::{contains_func, substitute_solution};
use crate::derivative::diff;
use crate::integral::integrate;

/// Check if an Atom is the number zero.
fn is_atom_zero(expr: Atom<'_>) -> bool {
    matches!(expr.node(), AtomNode::Num(0))
}

/// Attempt to solve a separable ODE.
///
/// A separable ODE can be written as $g(y)\,y' = f(x)$, so that
/// $\int g(y)\,dy = \int f(x)\,dx + C$.
pub(crate) fn solve_separable<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    // Extract y(x) name: the function symbol.
    let func_name = match func.node() {
        AtomNode::Fun(name, _) => *name,
        _ => return None,
    };

    let y = ctx.var(func_name.as_str());
    let _dy = ctx.fun("Derivative", &[func, ctx.var(var.as_str())]);

    // The equation is: expr = 0. We need to identify f(x) and g(y) such
    // that the equation can be written as g(y)*dy - f(x) = 0.
    //
    // Strategy: collect terms containing dy or func on one side, terms free
    // of func on the other.
    let (y_terms, x_terms) = separate_by_func(ctx, equation, func, var)?;

    if is_atom_zero(x_terms) && is_atom_zero(y_terms) {
        return None;
    }

    // Integrate both sides.
    // The x-side: integrate f(x) dx
    let x_integral = integrate(ctx, x_terms, var);
    // The y-side: integrate g(y) dy (replace x with y for integration)
    let y_in_x = substitute_var(ctx, y_terms, var, y);
    let y_integral = integrate(ctx, y_in_x, func_name);

    // Build implicit solution: y_integral - x_integral = C
    // For explicit, we'd need to solve for y, which is generally hard.
    // Return as implicit.
    let implicit = ctx.add(&[y_integral, ctx.mul(&[ctx.num(-1), x_integral])]);
    Some(ODESolution::Implicit(implicit))
}

/// Attempt to solve a first-order linear ODE: y' + p(x)*y = q(x).
///
/// Uses the integrating factor method:
/// $\mu(x) = e^{\int p(x)\,dx}$
/// $y = \frac{1}{\mu(x)} \int \mu(x)\,q(x)\,dx + \frac{C}{\mu(x)}$
pub(crate) fn solve_linear_first<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    let _x = ctx.var(var.as_str());
    let _y_sym = match func.node() {
        AtomNode::Fun(name, _) => ctx.var(name.as_str()),
        _ => return None,
    };

    // Extract p(x) and q(x) from y' + p(x)*y - q(x) = 0.
    let (p, q) = extract_linear_coeffs(ctx, equation, func, var)?;

    // Integrating factor: mu = exp(integral(p, x))
    let p_integral = integrate(ctx, p, var);
    let mu = ctx.fun("exp", &[p_integral]);

    // General solution: y = (1/mu) * (integral(mu * q, x) + C)
    let mu_q = ctx.mul(&[mu, q]);
    let int_mu_q = integrate(ctx, mu_q, var);

    // y = int_mu_q / mu + C/mu
    // We express this as: y = exp(-P) * (integral(exp(P)*q, x) + C)
    let neg_p_integral = ctx.mul(&[ctx.num(-1), p_integral]);
    let exp_neg_p = ctx.fun("exp", &[neg_p_integral]);

    // The particular solution (without C term).
    let particular = ctx.mul(&[exp_neg_p, int_mu_q]);

    // General solution includes C*exp(-P).
    // Return as explicit with C1 constant.
    let c1 = ctx.var("C1");
    let homogeneous = ctx.mul(&[c1, exp_neg_p]);
    let solution = ctx.add(&[particular, homogeneous]);

    Some(ODESolution::Explicit(solution))
}

/// Attempt to solve a Bernoulli ODE: y' + p(x)*y = q(x)*y^n.
///
/// Substitutes v = y^(1-n) to reduce to a linear ODE in v.
pub(crate) fn solve_bernoulli<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    let _x = ctx.var(var.as_str());
    let y_name = match func.node() {
        AtomNode::Fun(name, _) => name,
        _ => return None,
    };
    let _y = ctx.var(y_name.as_str());

    // Find the Bernoulli power n.
    let n = find_bernoulli_power(equation, func, var)?;

    if n == 0 || n == 1 {
        // n=0 is linear, n=1 is linear. These should be handled by LinearFirst.
        return None;
    }

    // Extract p(x) and the q(x) coefficient from the nonlinear term.
    let (p, q) = extract_linear_coeffs(ctx, equation, func, var)?;

    // Substitution: v = y^(1-n)
    // v' = (1-n) * y^(-n) * y'
    // The ODE becomes: v' + (1-n)*p(x)*v = (1-n)*q(x)
    let one_minus_n = 1 - n;
    let new_p = ctx.mul(&[ctx.num(one_minus_n), p]);
    let new_q = ctx.mul(&[ctx.num(one_minus_n), q]);

    // Solve the linear ODE in v: v' + new_p * v = new_q
    let new_p_integral = integrate(ctx, new_p, var);
    let mu = ctx.fun("exp", &[new_p_integral]);

    let mu_new_q = ctx.mul(&[mu, new_q]);
    let int_mu_new_q = integrate(ctx, mu_new_q, var);

    let neg_new_p_int = ctx.mul(&[ctx.num(-1), new_p_integral]);
    let exp_neg = ctx.fun("exp", &[neg_new_p_int]);

    let v_particular = ctx.mul(&[exp_neg, int_mu_new_q]);
    let c1 = ctx.var("C1");
    let v_homogeneous = ctx.mul(&[c1, exp_neg]);
    let v_solution = ctx.add(&[v_particular, v_homogeneous]);

    // Back-substitute: y = v^(1/(1-n))
    let power = if one_minus_n == 1 {
        // v = y, so y = v directly.
        v_solution
    } else {
        // y = v^(1/(1-n))
        // We represent 1/(1-n) as a rational power.
        let exponent = ctx.pow(ctx.num(one_minus_n), ctx.num(-1));
        ctx.pow(v_solution, exponent)
    };

    Some(ODESolution::Explicit(power))
}

/// Attempt to solve an exact ODE: M(x,y) + N(x,y)*y' = 0
/// where dM/dy = dN/dx.
///
/// The solution is F(x,y) = C where dF/dx = M and dF/dy = N.
///
/// When the equation is not exact as given, an integrating factor is
/// attempted: if $(M_y - N_x)/N$ depends only on $x$, then
/// $\mu(x) = \exp(\int (M_y - N_x)/N\,dx)$ makes the equation exact;
/// symmetrically, if $(N_x - M_y)/M$ depends only on $y$, then
/// $\mu(y) = \exp(\int (N_x - M_y)/M\,dy)$ works.
pub(crate) fn solve_exact<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    // Split into M + N*y' = 0 (y' must appear linearly).
    let (m, n) = super::classify::split_mn(ctx, equation, func, var)?;

    let y_name = match func.node() {
        AtomNode::Fun(name, _) => *name,
        _ => return None,
    };
    let y_var = ctx.var(y_name.as_str());

    // Work with bivariate expressions: replace y(x) by a plain symbol y.
    let m_sub = super::classify::replace_atom(ctx, m, func, y_var);
    let n_sub = super::classify::replace_atom(ctx, n, func, y_var);

    // Exactness check: dM/dy == dN/dx.
    let (m_eff, n_eff) = if partials_equal(ctx, m_sub, n_sub, y_name, var) {
        (m_sub, n_sub)
    } else {
        // Try integrating factors.
        find_integrating_factor(ctx, m_sub, n_sub, y_name, var)?
    };

    // Integrate M with respect to x (treating y as constant).
    let f_partial = integrate(ctx, m_eff, var);

    // Compute dF/dy of the partial integral.
    let df_dy = diff(ctx, f_partial, y_name);

    // The correction: g(y) = N - dF/dy. Use like-term collection so that
    // equivalent terms with different coefficient representations cancel.
    let correction =
        super::util::collect_terms(ctx, ctx.add(&[n_eff, ctx.mul(&[ctx.num(-1), df_dy])]));

    // If correction depends only on y, integrate it.
    if !contains_x(correction, var) {
        // Correction is a function of y only (or constant).
        let g_y = integrate(ctx, correction, y_name);
        let solution = ctx.add(&[f_partial, g_y]);
        // Substitute the plain symbol y back to y(x) for presentation.
        let solution = super::classify::replace_atom(ctx, solution, y_var, func);
        return Some(ODESolution::Implicit(solution));
    }

    None
}

/// Check whether dM/dy == dN/dx (after normalization).
fn partials_equal<'a>(
    ctx: &'a AtomArena<'a>,
    m: Atom<'a>,
    n: Atom<'a>,
    y_sym: Symbol,
    var: Symbol,
) -> bool {
    let dm_dy = diff(ctx, m, y_sym);
    let dn_dx = diff(ctx, n, var);
    let dm_norm = ocas_atom::normalize::normalize(ctx, dm_dy);
    let dn_norm = ocas_atom::normalize::normalize(ctx, dn_dx);
    if dm_norm.to_string() == dn_norm.to_string() {
        return true;
    }
    let difference =
        super::util::collect_terms(ctx, ctx.add(&[dm_dy, ctx.mul(&[ctx.num(-1), dn_dx])]));
    matches!(difference.node(), AtomNode::Num(0))
}

/// Attempt to find an integrating factor mu(x) or mu(y) that makes
/// M + N*y' = 0 exact. Returns the multiplied (mu*M, mu*N) on success.
fn find_integrating_factor<'a>(
    ctx: &'a AtomArena<'a>,
    m: Atom<'a>,
    n: Atom<'a>,
    y_sym: Symbol,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let dm_dy = diff(ctx, m, y_sym);
    let dn_dx = diff(ctx, n, var);

    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);

    // Candidate 1: (M_y - N_x)/N depends only on x.
    // mu(x) = exp(integral((M_y - N_x)/N dx)).
    let diff1 = simplify(
        ctx,
        ctx.add(&[dm_dy, ctx.mul(&[ctx.num(-1), dn_dx])]),
        &rules,
        20,
    );
    let ratio1 = ocas_atom::normalize::normalize(
        ctx,
        simplify(ctx, ctx.mul(&[diff1, ctx.pow(n, ctx.num(-1))]), &rules, 20),
    );
    if !ratio1.to_string().contains(y_sym.as_str()) && !contains_fun_named(ratio1, y_sym) {
        let exponent = integrate(ctx, ratio1, var);
        let mu = exp_simplify(ctx, exponent);
        let m_new = ctx.mul(&[mu, m]);
        let n_new = ctx.mul(&[mu, n]);
        // Verify exactness of the multiplied pair.
        if partials_equal(ctx, m_new, n_new, y_sym, var) {
            return Some((m_new, n_new));
        }
    }

    // Candidate 2: (N_x - M_y)/M depends only on y.
    // mu(y) = exp(integral((N_x - M_y)/M dy)).
    let diff2 = simplify(
        ctx,
        ctx.add(&[dn_dx, ctx.mul(&[ctx.num(-1), dm_dy])]),
        &rules,
        20,
    );
    let ratio2 = ocas_atom::normalize::normalize(
        ctx,
        simplify(ctx, ctx.mul(&[diff2, ctx.pow(m, ctx.num(-1))]), &rules, 20),
    );
    if !contains_x(ratio2, var) {
        let exponent = integrate(ctx, ratio2, y_sym);
        let mu = exp_simplify(ctx, exponent);
        let m_new = ctx.mul(&[mu, m]);
        let n_new = ctx.mul(&[mu, n]);
        if partials_equal(ctx, m_new, n_new, y_sym, var) {
            return Some((m_new, n_new));
        }
    }

    None
}

/// Simplify `exp(exponent)` when the exponent is a constant times a single
/// logarithm: `exp(k * log(u)) = u^k`. Falls back to the literal `exp` form.
///
/// Numeric constant factors inside the logarithm's argument are dropped:
/// they only scale the integrating factor by a constant, which is harmless.
fn exp_simplify<'a>(ctx: &'a AtomArena<'a>, exponent: Atom<'a>) -> Atom<'a> {
    super::util::exp_simplify(ctx, exponent)
}

/// Check whether expr contains a Fun node with the given symbol name.
fn contains_fun_named<'a>(expr: Atom<'a>, name: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_fun_named(*a, name))
        }
        AtomNode::Pow(base, exp) => {
            contains_fun_named(*base, name) || contains_fun_named(*exp, name)
        }
        AtomNode::Fun(n, args) => *n == name || args.iter().any(|a| contains_fun_named(*a, name)),
    }
}

/// Attempt to solve a homogeneous ODE: y' = f(y/x).
///
/// Substitutes v = y/x (so y = vx, y' = v + xv') to get a separable ODE in v and x.
pub(crate) fn solve_homogeneous<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    let x = ctx.var(var.as_str());
    let y_name = match func.node() {
        AtomNode::Fun(name, _) => name,
        _ => return None,
    };
    let y = ctx.var(y_name.as_str());
    let v = ctx.var("v");

    // Substitute y = v*x into the equation.
    let y_replacement = ctx.mul(&[v, x]);
    let substituted = substitute_solution(ctx, equation, func, y_replacement, var);

    // Also substitute y' = v + x*v'.
    // We need to express the substituted equation in terms of v and v' = Derivative(v, x).
    // This requires expanding and simplifying.

    // For a homogeneous ODE of the form y' = F(y/x):
    // After substitution: v + x*v' = F(v)
    // So: x*v' = F(v) - v
    // Separable: dv/(F(v)-v) = dx/x

    // Heuristic: try to collect the equation into the form x*v' = G(v).
    // The v' atom is only needed conceptually here; actual derivative
    // handling goes through `separate_by_var` below.
    let _dv_symbol = ctx.var("dv");

    // Simplify the substituted equation.
    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let simplified = simplify(ctx, substituted, &rules, 20);

    // Try to extract F(v) - v from the simplified equation.
    // This is a best-effort heuristic.
    let (v_part, x_part) = separate_by_var(ctx, simplified, var, v)?;

    // Integrate both sides.
    let x_integral = integrate(ctx, x_part, var);
    let v_sym = Symbol::new("v");
    let v_integral = integrate(ctx, v_part, v_sym);

    // Back-substitute v = y/x.
    let v_in_yx = ctx.mul(&[y, ctx.pow(x, ctx.num(-1))]);
    let v_sol = substitute_var(ctx, v_integral, v_sym, v_in_yx);

    let implicit = ctx.add(&[v_sol, ctx.mul(&[ctx.num(-1), x_integral])]);
    Some(ODESolution::Implicit(implicit))
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Separate an equation into terms containing func and terms free of func.
/// Returns (func_terms, free_terms) such that equation = func_terms + free_terms.
fn separate_by_func<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match equation.node() {
        AtomNode::Add(args) => {
            let mut y_terms = Vec::new();
            let mut x_terms = Vec::new();
            for a in args.iter() {
                if contains_func(*a, func, var) {
                    y_terms.push(*a);
                } else {
                    x_terms.push(*a);
                }
            }
            if y_terms.is_empty() || x_terms.is_empty() {
                return None;
            }
            let y_sum = ctx.add(&y_terms);
            let x_sum = ctx.add(&x_terms);
            Some((y_sum, ctx.mul(&[ctx.num(-1), x_sum])))
        }
        _ => {
            if contains_func(equation, func, var) {
                Some((equation, ctx.num(0)))
            } else {
                None
            }
        }
    }
}

/// Extract p(x) and q(x) from a first-order linear ODE written as
/// `Derivative(y(x), x) + p(x)*y(x) + rest = 0`.
///
/// Returns `(p(x), -rest)` where rest is everything not involving y or y'.
fn extract_linear_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let x = ctx.var(var.as_str());
    let dy = ctx.fun("Derivative", &[func, x]);

    // Decompose equation into sum of terms.
    let terms = flatten_add(equation);

    let mut p_coeff: Option<Atom<'a>> = None; // coefficient of y
    let mut q_neg: Vec<Atom<'a>> = Vec::new(); // terms not involving y or y'

    for term in &terms {
        let s = term.to_string();
        // Check if this term is y' (Derivative(y(x), x)).
        if s == dy.to_string() {
            // y' coefficient is 1, which is the standard form.
            continue;
        }
        // Bare y(x): coefficient is 1.
        if s == func.to_string() {
            p_coeff = Some(ctx.num(1));
            continue;
        }
        // Check if this term contains y(x) as a factor.
        if contains_func(*term, func, var) && !is_derivative(*term, func, var) {
            // This should be p(x)*y(x).
            // Extract p(x) by dividing out y(x).
            if let AtomNode::Mul(args) = term.node() {
                let factors: Vec<_> = args
                    .iter()
                    .filter(|a| a.to_string() != func.to_string())
                    .copied()
                    .collect();
                if factors.is_empty() {
                    p_coeff = Some(ctx.num(1));
                } else if factors.len() == 1 {
                    p_coeff = Some(factors[0]);
                } else {
                    p_coeff = Some(ctx.mul(&factors));
                }
            } else if term.to_string() == func.to_string() {
                p_coeff = Some(ctx.num(1));
            }
        } else if !contains_func(*term, func, var) {
            // This is a forcing term (free of y).
            q_neg.push(*term);
        }
    }

    let p = p_coeff.unwrap_or(ctx.num(0));
    let q = if q_neg.is_empty() {
        ctx.num(0)
    } else if q_neg.len() == 1 {
        ctx.mul(&[ctx.num(-1), q_neg[0]])
    } else {
        ctx.mul(&[ctx.num(-1), ctx.add(&q_neg)])
    };

    Some((p, q))
}

/// Find the Bernoulli power n from equation containing y^n alongside linear y terms.
fn find_bernoulli_power<'a>(equation: Atom<'a>, func: Atom<'a>, var: Symbol) -> Option<i64> {
    find_power_inner(equation, func, var)
}

fn find_power_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> Option<i64> {
    match expr.node() {
        AtomNode::Add(args) => args.iter().find_map(|a| find_power_inner(*a, func, var)),
        AtomNode::Mul(args) => args.iter().find_map(|a| find_power_inner(*a, func, var)),
        AtomNode::Pow(base, exp) => {
            if contains_func(*base, func, var)
                && let AtomNode::Num(n) = exp.node()
                && *n >= 2
            {
                return Some(*n);
            }
            find_power_inner(*base, func, var).or_else(|| find_power_inner(*exp, func, var))
        }
        _ => None,
    }
}

/// Separate terms into those involving `dep_var` and those free of it.
fn separate_by_var<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    _x_var: Symbol,
    dep_var: Atom<'a>,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match expr.node() {
        AtomNode::Add(args) => {
            let mut dep_terms = Vec::new();
            let mut free_terms = Vec::new();
            for a in args.iter() {
                if contains_var_atom(*a, dep_var) {
                    dep_terms.push(*a);
                } else {
                    free_terms.push(*a);
                }
            }
            if dep_terms.is_empty() || free_terms.is_empty() {
                return None;
            }
            Some((ctx.add(&dep_terms), ctx.add(&free_terms)))
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

/// Flatten an addition tree into a list of terms.
fn flatten_add<'a>(expr: Atom<'a>) -> Vec<Atom<'a>> {
    match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    }
}

/// Check if expr is `Derivative(func, var)`.
fn is_derivative<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            *name == Symbol::new("Derivative")
                && args.len() >= 2
                && args[0].to_string() == func.to_string()
                && args[1].to_string() == var.as_str()
        }
        _ => false,
    }
}

/// Check if expr contains a specific variable.
fn contains_x<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => *v == var,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_x(*a, var)),
        AtomNode::Pow(base, exp) => contains_x(*base, var) || contains_x(*exp, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_x(*a, var)),
    }
}

/// Check if expr contains a specific Atom (by string comparison).
fn contains_var_atom<'a>(expr: Atom<'a>, target: Atom<'a>) -> bool {
    let target_str = target.to_string();
    contains_str(expr, &target_str)
}

fn contains_str<'a>(expr: Atom<'a>, target: &str) -> bool {
    if expr.to_string() == target {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_str(*a, target)),
        AtomNode::Pow(base, exp) => contains_str(*base, target) || contains_str(*exp, target),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_str(*a, target)),
    }
}

/// Replace all occurrences of `var_sym` with `replacement` in `expr`.
fn substitute_var<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var_sym: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    crate::series::substitute(ctx, expr, var_sym, replacement)
}
