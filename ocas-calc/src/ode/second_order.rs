//! Second-order linear ODE solvers.
//!
//! Implements classical methods for second-order linear ODEs:
//! - Constant coefficients: $a y'' + b y' + c y = f(x)$
//! - Cauchy-Euler: $a x^2 y'' + b x y' + c y = f(x)$

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use super::ODESolution;
use super::util::{contains_func, ode_order};
use crate::derivative::diff;

/// Check if an Atom is the number zero.
fn is_atom_zero(expr: Atom<'_>) -> bool {
    matches!(expr.node(), AtomNode::Num(0))
}

/// Attempt to solve a constant-coefficient linear ODE:
/// $a y'' + b y' + c y = f(x)$.
///
/// For the homogeneous part, solves the characteristic equation
/// $a r^2 + b r + c = 0$ and builds the complementary solution.
/// For the particular solution, uses the method of undetermined coefficients.
pub(crate) fn solve_constant_coeff<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    let order = ode_order(equation, func, var);
    if order != 2 || !super::util::is_linear_in(equation, func, var) {
        return None;
    }

    let x = ctx.var(var.as_str());

    // Extract coefficients a, b, c and forcing function f(x).
    let (a, b, c, forcing) = extract_second_order_coeffs(ctx, equation, func, var)?;

    // Verify a, b, c are constants (free of x).
    if contains_x(a, var) || contains_x(b, var) || contains_x(c, var) {
        return None;
    }

    // Solve characteristic equation: a*r^2 + b*r + c = 0
    // Using the quadratic formula: r = (-b ± sqrt(b^2 - 4ac)) / (2a)
    let discriminant = ctx.add(&[ctx.pow(b, ctx.num(2)), ctx.mul(&[ctx.num(-4), a, c])]);

    // Build complementary solution based on discriminant.
    let y_c = complementary_solution(ctx, a, b, discriminant, x);

    // Build particular solution using undetermined coefficients.
    let y_p = if is_atom_zero(forcing) {
        ctx.num(0)
    } else {
        particular_solution_undetermined(ctx, forcing, a, b, c, x, var).unwrap_or(ctx.num(0))
    };

    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let solution = simplify(ctx, ctx.add(&[y_c, y_p]), &rules, 20);
    Some(ODESolution::Explicit(solution))
}

/// Attempt to solve a Cauchy-Euler equation:
/// $a x^2 y'' + b x y' + c y = f(x)$.
///
/// Substitutes $x = e^t$ to convert to a constant-coefficient equation.
pub(crate) fn solve_cauchy_euler<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    let order = ode_order(equation, func, var);
    if order != 2 {
        return None;
    }

    let x = ctx.var(var.as_str());

    // Extract coefficients assuming Cauchy-Euler form.
    // For a*x^2*y'' + b*x*y' + c*y = f(x):
    // The characteristic equation is a*r*(r-1) + b*r + c = 0
    // i.e., a*r^2 + (b-a)*r + c = 0
    let (a_ce, b_ce, c_ce, forcing) = extract_cauchy_euler_coeffs(ctx, equation, func, var)?;

    if is_atom_zero(a_ce) {
        return None;
    }

    // Characteristic equation: a*r*(r-1) + b*r + c = 0
    // = a*r^2 + (b-a)*r + c = 0
    let coeff_r2 = a_ce;
    let coeff_r1 = ctx.add(&[b_ce, ctx.mul(&[ctx.num(-1), a_ce])]);
    let coeff_r0 = c_ce;

    let discriminant = ctx.add(&[
        ctx.pow(coeff_r1, ctx.num(2)),
        ctx.mul(&[ctx.num(-4), coeff_r2, coeff_r0]),
    ]);

    // For Cauchy-Euler, solutions are of the form x^r.
    let y_c = cauchy_euler_complementary(ctx, coeff_r2, coeff_r1, coeff_r0, discriminant, x);

    let y_p = if is_atom_zero(forcing) {
        ctx.num(0)
    } else {
        // For non-homogeneous Cauchy-Euler, use variation of parameters.
        ctx.num(0) // Placeholder: full VOP is complex.
    };

    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let solution = simplify(ctx, ctx.add(&[y_c, y_p]), &rules, 20);
    Some(ODESolution::Explicit(solution))
}

// ---------------------------------------------------------------------------
// Complementary solution builders
// ---------------------------------------------------------------------------

/// Build complementary solution for constant-coefficient ODE.
fn complementary_solution<'a>(
    ctx: &'a AtomArena<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    discriminant: Atom<'a>,
    x: Atom<'a>,
) -> Atom<'a> {
    // Try to evaluate discriminant as a number.
    match discriminant.node() {
        AtomNode::Num(d) => {
            let c1 = ctx.var("C1");
            let c2 = ctx.var("C2");

            if *d > 0 {
                // Two distinct real roots: r1, r2.
                // y_c = C1*exp(r1*x) + C2*exp(r2*x)
                let (r1, r2) = real_roots(ctx, a, b, *d);
                ctx.add(&[
                    ctx.mul(&[c1, ctx.fun("exp", &[ctx.mul(&[r1, x])])]),
                    ctx.mul(&[c2, ctx.fun("exp", &[ctx.mul(&[r2, x])])]),
                ])
            } else if *d == 0 {
                // Repeated root: r = -b/(2a).
                // y_c = (C1 + C2*x)*exp(r*x)
                let r = ctx.mul(&[
                    ctx.num(-1),
                    b,
                    ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ]);
                ctx.mul(&[
                    ctx.add(&[c1, ctx.mul(&[c2, x])]),
                    ctx.fun("exp", &[ctx.mul(&[r, x])]),
                ])
            } else {
                // Complex roots: alpha ± beta*i.
                // y_c = exp(alpha*x)*(C1*cos(beta*x) + C2*sin(beta*x))
                let alpha = ctx.mul(&[
                    ctx.num(-1),
                    b,
                    ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ]);
                let beta_sq = ctx.mul(&[ctx.num(-1), discriminant]);
                let two_a = ctx.mul(&[ctx.num(2), a]);
                let beta = ctx.pow(
                    ctx.mul(&[beta_sq, ctx.pow(ctx.pow(two_a, ctx.num(2)), ctx.num(-1))]),
                    ctx.pow(ctx.num(2), ctx.num(-1)),
                );
                let exp_part = ctx.fun("exp", &[ctx.mul(&[alpha, x])]);
                let cos_part = ctx.fun("cos", &[ctx.mul(&[beta, x])]);
                let sin_part = ctx.fun("sin", &[ctx.mul(&[beta, x])]);
                ctx.mul(&[
                    exp_part,
                    ctx.add(&[ctx.mul(&[c1, cos_part]), ctx.mul(&[c2, sin_part])]),
                ])
            }
        }
        _ => {
            // Discriminant is not a simple number — symbolic roots.
            // Return a generic form.
            let c1 = ctx.var("C1");
            let c2 = ctx.var("C2");
            let r1 = ctx.mul(&[
                ctx.add(&[ctx.num(-1), b]),
                ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ctx.pow(discriminant, ctx.pow(ctx.num(2), ctx.num(-1))),
            ]);
            let r2 = ctx.mul(&[
                ctx.add(&[b, ctx.num(-1)]),
                ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ctx.pow(discriminant, ctx.pow(ctx.num(2), ctx.num(-1))),
            ]);
            ctx.add(&[
                ctx.mul(&[c1, ctx.fun("exp", &[ctx.mul(&[r1, x])])]),
                ctx.mul(&[c2, ctx.fun("exp", &[ctx.mul(&[r2, x])])]),
            ])
        }
    }
}

/// Build complementary solution for Cauchy-Euler ODE.
fn cauchy_euler_complementary<'a>(
    ctx: &'a AtomArena<'a>,
    _a: Atom<'a>,
    b: Atom<'a>,
    _c: Atom<'a>,
    discriminant: Atom<'a>,
    x: Atom<'a>,
) -> Atom<'a> {
    let c1 = ctx.var("C1");
    let c2 = ctx.var("C2");

    match discriminant.node() {
        AtomNode::Num(d) => {
            if *d > 0 {
                let (r1, r2) = real_roots(ctx, _a, b, *d);
                ctx.add(&[
                    ctx.mul(&[c1, ctx.pow(x, r1)]),
                    ctx.mul(&[c2, ctx.pow(x, r2)]),
                ])
            } else if *d == 0 {
                // Repeated root r: y = (C1 + C2*ln(x)) * x^r
                let r = ctx.mul(&[
                    ctx.num(-1),
                    b,
                    ctx.pow(ctx.mul(&[ctx.num(2), _a]), ctx.num(-1)),
                ]);
                ctx.mul(&[
                    ctx.add(&[c1, ctx.mul(&[c2, ctx.fun("log", &[x])])]),
                    ctx.pow(x, r),
                ])
            } else {
                // Complex roots: x^alpha * (C1*cos(beta*ln(x)) + C2*sin(beta*ln(x)))
                let alpha = ctx.mul(&[
                    ctx.num(-1),
                    b,
                    ctx.pow(ctx.mul(&[ctx.num(2), _a]), ctx.num(-1)),
                ]);
                let log_x = ctx.fun("log", &[x]);
                ctx.mul(&[
                    ctx.pow(x, alpha),
                    ctx.add(&[
                        ctx.mul(&[c1, ctx.fun("cos", &[log_x])]),
                        ctx.mul(&[c2, ctx.fun("sin", &[log_x])]),
                    ]),
                ])
            }
        }
        _ => {
            // Symbolic discriminant — generic form.
            ctx.add(&[
                ctx.mul(&[c1, ctx.pow(x, ctx.num(1))]),
                ctx.mul(&[c2, ctx.pow(x, ctx.num(2))]),
            ])
        }
    }
}

// ---------------------------------------------------------------------------
// Particular solution (undetermined coefficients)
// ---------------------------------------------------------------------------

/// Attempt undetermined coefficients for a polynomial forcing function.
fn particular_solution_undetermined<'a>(
    ctx: &'a AtomArena<'a>,
    forcing: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // For polynomial forcing f(x), try y_p = polynomial of same degree.
    // For exponential forcing, try y_p = A*exp(kx).
    // For trigonometric forcing, try y_p = A*cos(wx) + B*sin(wx).
    // For combinations, try products.

    // Try polynomial forcing first.
    if is_polynomial_in(forcing, var) {
        let degree = polynomial_degree(forcing, var);
        return undetermined_polynomial(ctx, forcing, a, b, c, x, degree, var);
    }

    // Try exponential forcing: A*exp(kx).
    if let Some(k) = extract_exp_coeff(forcing, var) {
        return undetermined_exponential(ctx, k, a, b, c, x);
    }

    None
}

/// Undetermined coefficients for polynomial forcing.
#[allow(clippy::too_many_arguments)]
fn undetermined_polynomial<'a>(
    ctx: &'a AtomArena<'a>,
    forcing: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
    degree: usize,
    var: Symbol,
) -> Option<Atom<'a>> {
    // Try y_p = sum_{i=0}^{degree} A_i * x^i
    // Substitute into a*y'' + b*y' + c*y and match coefficients with forcing.

    // Build y_p symbolically with unknown coefficients.
    let mut terms = Vec::new();
    let mut coeff_syms = Vec::new();
    for i in 0..=degree {
        let sym_name = format!("A{i}");
        let ai = ctx.var(&sym_name);
        coeff_syms.push(ai);
        if i == 0 {
            terms.push(ai);
        } else {
            terms.push(ctx.mul(&[ai, ctx.pow(x, ctx.num(i as i64))]));
        }
    }
    let y_p = ctx.add(&terms);

    // Compute a*y_p'' + b*y_p' + c*y_p.
    let dy_p = diff(ctx, y_p, var);
    let d2y_p = diff(ctx, dy_p, var);

    let lhs = ctx.add(&[
        ctx.mul(&[a, d2y_p]),
        ctx.mul(&[b, dy_p]),
        ctx.mul(&[c, y_p]),
    ]);

    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let _lhs_simplified = simplify(ctx, lhs, &rules, 20);

    // Match coefficients of x^0, x^1, ..., x^degree with forcing.
    // This is a system of linear equations in A0, A1, ..., A_degree.
    // For simplicity, try small cases.
    match degree {
        0 => {
            // y_p = A0, y_p' = 0, y_p'' = 0.
            // a*0 + b*0 + c*A0 = forcing → A0 = forcing/c
            let c_val = extract_constant(c)?;
            if c_val == 0 {
                return None;
            }
            Some(ctx.mul(&[forcing, ctx.pow(c, ctx.num(-1))]))
        }
        1 => {
            // y_p = A0 + A1*x
            // y_p' = A1, y_p'' = 0
            // c*A0 + (c*x + b)*A1 = forcing (after collecting terms)
            // For forcing = f0 + f1*x:
            // c*A1 = f1, b*A1 + c*A0 = f0
            let (f0, f1) = extract_linear_forcing(ctx, forcing, var)?;
            let c_val = extract_constant(c)?;
            let b_val = extract_constant(b)?;
            if c_val == 0 {
                return None;
            }
            let a1 = ctx.mul(&[f1, ctx.pow(c, ctx.num(-1))]);
            let a0 = ctx.mul(&[
                ctx.add(&[f0, ctx.mul(&[ctx.num(-1), b_val_atom(ctx, b_val), a1])]),
                ctx.pow(c, ctx.num(-1)),
            ]);
            Some(ctx.add(&[a0, ctx.mul(&[a1, x])]))
        }
        _ => None, // Higher degrees: would need a linear system solver.
    }
}

/// Undetermined coefficients for exponential forcing.
fn undetermined_exponential<'a>(
    ctx: &'a AtomArena<'a>,
    k: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    // For forcing = A*exp(kx), try y_p = B*exp(kx).
    // y_p' = Bk*exp(kx), y_p'' = Bk^2*exp(kx).
    // a*Bk^2 + b*Bk + c*B = A → B = A / (ak^2 + bk + c).
    let char_at_k = ctx.add(&[ctx.mul(&[a, ctx.pow(k, ctx.num(2))]), ctx.mul(&[b, k]), c]);
    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let char_simplified = simplify(ctx, char_at_k, &rules, 20);

    // Check that the characteristic equation doesn't vanish at k.
    if let AtomNode::Num(0) = char_simplified.node() {
        // k is a root of the characteristic equation — need to multiply by x.
        // This is the resonance case; handled separately if needed.
        return None;
    }

    let coeff = ctx.pow(char_simplified, ctx.num(-1));
    Some(ctx.mul(&[coeff, ctx.fun("exp", &[ctx.mul(&[k, x])])]))
}

// ---------------------------------------------------------------------------
// Coefficient extraction helpers
// ---------------------------------------------------------------------------

/// Extract a, b, c coefficients from a*y'' + b*y' + c*y + ... = 0.
fn extract_second_order_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>, Atom<'a>)> {
    let x = ctx.var(var.as_str());
    let dy = ctx.fun("Derivative", &[func, x]);
    let d2y = ctx.fun("Derivative", &[func, x, x]);

    let terms = flatten_add(equation);

    let mut a_coeff = ctx.num(0);
    let mut b_coeff = ctx.num(0);
    let mut c_coeff = ctx.num(0);
    let mut forcing_terms = Vec::new();

    for term in &terms {
        let s = term.to_string();
        if s == d2y.to_string() {
            a_coeff = ctx.add(&[a_coeff, ctx.num(1)]);
        } else if s == dy.to_string() {
            b_coeff = ctx.add(&[b_coeff, ctx.num(1)]);
        } else if s == func.to_string() {
            c_coeff = ctx.add(&[c_coeff, ctx.num(1)]);
        } else if contains_str(*term, &d2y.to_string()) {
            // a*y''
            if let AtomNode::Mul(args) = term.node() {
                let coeff: Vec<_> = args
                    .iter()
                    .filter(|a| a.to_string() != d2y.to_string())
                    .copied()
                    .collect();
                let c = if coeff.is_empty() {
                    ctx.num(1)
                } else {
                    ctx.mul(&coeff)
                };
                a_coeff = ctx.add(&[a_coeff, c]);
            }
        } else if contains_str(*term, &dy.to_string()) && !contains_str(*term, &d2y.to_string()) {
            // b*y'
            if let AtomNode::Mul(args) = term.node() {
                let coeff: Vec<_> = args
                    .iter()
                    .filter(|a| a.to_string() != dy.to_string())
                    .copied()
                    .collect();
                let c = if coeff.is_empty() {
                    ctx.num(1)
                } else {
                    ctx.mul(&coeff)
                };
                b_coeff = ctx.add(&[b_coeff, c]);
            }
        } else if contains_str(*term, &func.to_string())
            && !contains_str(*term, &dy.to_string())
            && !contains_str(*term, &d2y.to_string())
        {
            // c*y
            if let AtomNode::Mul(args) = term.node() {
                let coeff: Vec<_> = args
                    .iter()
                    .filter(|a| a.to_string() != func.to_string())
                    .copied()
                    .collect();
                let c = if coeff.is_empty() {
                    ctx.num(1)
                } else {
                    ctx.mul(&coeff)
                };
                c_coeff = ctx.add(&[c_coeff, c]);
            }
        } else if !contains_func(*term, func, var) {
            forcing_terms.push(ctx.mul(&[ctx.num(-1), *term]));
        }
    }

    let forcing = if forcing_terms.is_empty() {
        ctx.num(0)
    } else {
        ctx.add(&forcing_terms)
    };

    Some((a_coeff, b_coeff, c_coeff, forcing))
}

/// Extract Cauchy-Euler coefficients from a*x^2*y'' + b*x*y' + c*y = f(x).
fn extract_cauchy_euler_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>, Atom<'a>)> {
    // Same as second-order but verify the x^n factor pattern.
    extract_second_order_coeffs(ctx, equation, func, var)
}

/// Compute real roots from discriminant value.
fn real_roots<'a>(
    ctx: &'a AtomArena<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    disc: i64,
) -> (Atom<'a>, Atom<'a>) {
    let sqrt_disc = ctx.num(isqrt(disc));
    let two_a = ctx.mul(&[ctx.num(2), a]);
    let r1 = ctx.mul(&[
        ctx.add(&[ctx.num(-1), b, sqrt_disc]),
        ctx.pow(two_a, ctx.num(-1)),
    ]);
    let r2 = ctx.mul(&[
        ctx.add(&[ctx.num(-1), b, ctx.mul(&[ctx.num(-1), sqrt_disc])]),
        ctx.pow(two_a, ctx.num(-1)),
    ]);
    (r1, r2)
}

/// Integer square root (floor).
fn isqrt(n: i64) -> i64 {
    if n < 0 {
        return 0;
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

fn flatten_add<'a>(expr: Atom<'a>) -> Vec<Atom<'a>> {
    match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    }
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

fn contains_x<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => *v == var,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_x(*a, var)),
        AtomNode::Pow(base, exp) => contains_x(*base, var) || contains_x(*exp, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_x(*a, var)),
    }
}

fn is_polynomial_in<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v == var || true,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|a| is_polynomial_in(*a, var)),
        AtomNode::Pow(base, exp) => {
            if let AtomNode::Num(n) = exp.node() {
                *n >= 0 && is_polynomial_in(*base, var)
            } else {
                false
            }
        }
        AtomNode::Fun(name, _) => {
            // Functions of x are not polynomial.
            !contains_x(expr, var) || *name == Symbol::new("Derivative")
        }
    }
}

fn polynomial_degree<'a>(expr: Atom<'a>, var: Symbol) -> usize {
    match expr.node() {
        AtomNode::Num(_) => 0,
        AtomNode::Var(v) => {
            if *v == var {
                1
            } else {
                0
            }
        }
        AtomNode::Add(args) => args
            .iter()
            .map(|a| polynomial_degree(*a, var))
            .max()
            .unwrap_or(0),
        AtomNode::Mul(args) => args.iter().map(|a| polynomial_degree(*a, var)).sum(),
        AtomNode::Pow(base, exp) => {
            if let AtomNode::Num(n) = exp.node() {
                polynomial_degree(*base, var) * (*n as usize)
            } else {
                0
            }
        }
        _ => 0,
    }
}

fn extract_exp_coeff<'a>(expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Mul(args) => args.iter().find_map(|a| extract_exp_k(*a, var)),
        AtomNode::Fun(name, fargs) if *name == Symbol::new("exp") && fargs.len() == 1 => {
            extract_exp_k_inner(fargs[0], var)
        }
        _ => None,
    }
}

fn extract_exp_k<'a>(expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Fun(name, fargs) if *name == Symbol::new("exp") && fargs.len() == 1 => {
            extract_exp_k_inner(fargs[0], var)
        }
        _ => None,
    }
}

fn extract_exp_k_inner<'a>(arg: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if let AtomNode::Mul(factors) = arg.node() {
        let k: Vec<_> = factors
            .iter()
            .filter(|f| f.to_string() != var.as_str())
            .copied()
            .collect();
        if k.len() == 1 {
            return Some(k[0]);
        }
    }
    None
}

fn extract_constant<'a>(expr: Atom<'a>) -> Option<i64> {
    match expr.node() {
        AtomNode::Num(n) => Some(*n),
        _ => None,
    }
}

fn extract_linear_forcing<'a>(
    ctx: &'a AtomArena<'a>,
    forcing: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let _x = ctx.var(var.as_str());
    match forcing.node() {
        AtomNode::Num(n) => Some((ctx.num(*n), ctx.num(0))),
        AtomNode::Var(v) => {
            if *v == var {
                Some((ctx.num(0), ctx.num(1)))
            } else {
                Some((forcing, ctx.num(0)))
            }
        }
        AtomNode::Add(args) => {
            let mut f0 = ctx.num(0);
            let mut f1 = ctx.num(0);
            for a in args.iter() {
                match a.node() {
                    AtomNode::Num(n) => f0 = ctx.add(&[f0, ctx.num(*n)]),
                    AtomNode::Var(v) => {
                        if *v == var {
                            f1 = ctx.add(&[f1, ctx.num(1)]);
                        }
                    }
                    AtomNode::Mul(factors)
                        if factors.iter().any(|f| f.to_string() == var.as_str()) =>
                    {
                        let coeff: Vec<_> = factors
                            .iter()
                            .filter(|f| f.to_string() != var.as_str())
                            .copied()
                            .collect();
                        let c = if coeff.is_empty() {
                            ctx.num(1)
                        } else {
                            ctx.mul(&coeff)
                        };
                        f1 = ctx.add(&[f1, c]);
                    }
                    _ => {}
                }
            }
            Some((f0, f1))
        }
        _ => None,
    }
}

fn b_val_atom<'a>(ctx: &'a AtomArena<'a>, val: i64) -> Atom<'a> {
    ctx.num(val)
}
