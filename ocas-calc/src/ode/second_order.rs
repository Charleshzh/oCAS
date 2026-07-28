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
use crate::integral::integrate;

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
    // Fold numeric factors so the discriminant reduces to a Num when possible.
    let discriminant = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.pow(b, ctx.num(2)), ctx.mul(&[ctx.num(-4), a, c])]),
    );

    // Build complementary solution based on discriminant.
    let y_c = complementary_solution(ctx, a, b, discriminant, x);

    // Build particular solution: try undetermined coefficients first,
    // then fall back to variation of parameters.
    let y_p = if is_atom_zero(forcing) {
        ctx.num(0)
    } else {
        particular_solution_undetermined(ctx, forcing, a, b, c, x, var)
            .or_else(|| {
                // Standard form: g(x) = f(x)/a.
                let g =
                    super::util::collect_terms(ctx, ctx.mul(&[forcing, ctx.pow(a, ctx.num(-1))]));
                let (y1, y2) = constant_coeff_basis(ctx, a, b, discriminant, x)?;
                variation_of_parameters(ctx, y1, y2, g, var)
            })
            .unwrap_or(ctx.num(0))
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

    // The extracted coefficients carry the Cauchy-Euler x-power factors:
    // a_ce = a*x^2, b_ce = b*x, c_ce = c. Divide them out and require the
    // reduced coefficients to be constant (free of x).
    let a_const = super::util::collect_terms(
        ctx,
        ctx.mul(&[a_ce, ctx.pow(ctx.pow(x, ctx.num(2)), ctx.num(-1))]),
    );
    let b_const = super::util::collect_terms(ctx, ctx.mul(&[b_ce, ctx.pow(x, ctx.num(-1))]));
    let c_const = c_ce;
    if contains_x(a_const, var) || contains_x(b_const, var) || contains_x(c_const, var) {
        return None;
    }

    // Characteristic equation: a*r*(r-1) + b*r + c = 0
    // = a*r^2 + (b-a)*r + c = 0
    let coeff_r2 = a_const;
    let coeff_r1 = ctx.add(&[b_const, ctx.mul(&[ctx.num(-1), a_const])]);
    let coeff_r0 = c_const;

    let discriminant = super::util::collect_terms(
        ctx,
        ctx.add(&[
            ctx.pow(coeff_r1, ctx.num(2)),
            ctx.mul(&[ctx.num(-4), coeff_r2, coeff_r0]),
        ]),
    );

    // For Cauchy-Euler, solutions are of the form x^r.
    let y_c = cauchy_euler_complementary(ctx, coeff_r2, coeff_r1, coeff_r0, discriminant, x);

    let y_p = if is_atom_zero(forcing) {
        ctx.num(0)
    } else {
        // For non-homogeneous Cauchy-Euler, use variation of parameters
        // on the standard-form equation y'' + p(x)y' + q(x)y = g(x)
        // where g(x) = f(x) / (a*x^2).
        let a_x2 = ctx.mul(&[coeff_r2, ctx.pow(x, ctx.num(2))]);
        let g = super::util::collect_terms(ctx, ctx.mul(&[forcing, ctx.pow(a_x2, ctx.num(-1))]));
        cauchy_euler_basis(ctx, coeff_r2, coeff_r1, discriminant, x)
            .and_then(|(y1, y2)| variation_of_parameters(ctx, y1, y2, g, var))
            .unwrap_or(ctx.num(0))
    };

    let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let solution = simplify(ctx, ctx.add(&[y_c, y_p]), &rules, 20);
    Some(ODESolution::Explicit(solution))
}

// ---------------------------------------------------------------------------
// Complementary solution builders
// ---------------------------------------------------------------------------

/// Build the two fundamental solutions (y1, y2) of the homogeneous
/// constant-coefficient equation a*y'' + b*y' + c*y = 0.
fn constant_coeff_basis<'a>(
    ctx: &'a AtomArena<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    discriminant: Atom<'a>,
    x: Atom<'a>,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match discriminant.node() {
        AtomNode::Num(d) => {
            if *d > 0 {
                let (r1, r2) = real_roots(ctx, a, b, *d);
                Some((
                    ctx.fun("exp", &[ctx.mul(&[r1, x])]),
                    ctx.fun("exp", &[ctx.mul(&[r2, x])]),
                ))
            } else if *d == 0 {
                let r = ctx.mul(&[
                    ctx.num(-1),
                    b,
                    ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ]);
                let e = ctx.fun("exp", &[ctx.mul(&[r, x])]);
                Some((e, ctx.mul(&[x, e])))
            } else {
                let alpha = super::util::collect_terms(
                    ctx,
                    ctx.mul(&[
                        ctx.num(-1),
                        b,
                        ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                    ]),
                );
                // beta = sqrt(-d) / (2a). When a is numeric, evaluate the
                // rational under the square root exactly; keep the symbolic
                // radical only when the numerator is not a perfect square.
                let a_folded = super::util::collect_terms(ctx, a);
                let beta = if let AtomNode::Num(a_num) = a_folded.node() {
                    let num = -*d; // > 0
                    let den = 2 * a_num.abs();
                    let sn = isqrt(num);
                    if sn * sn == num {
                        // beta = sn/den exactly (rational).
                        ctx.mul(&[ctx.num(sn), ctx.pow(ctx.num(den), ctx.num(-1))])
                    } else {
                        // beta = sqrt(num)/den with symbolic numerator.
                        ctx.mul(&[
                            ctx.pow(ctx.num(num), ctx.pow(ctx.num(2), ctx.num(-1))),
                            ctx.pow(ctx.num(den), ctx.num(-1)),
                        ])
                    }
                } else {
                    let beta_sq = ctx.mul(&[ctx.num(-1), discriminant]);
                    let two_a = ctx.mul(&[ctx.num(2), a]);
                    super::util::collect_terms(
                        ctx,
                        ctx.pow(
                            ctx.mul(&[beta_sq, ctx.pow(ctx.pow(two_a, ctx.num(2)), ctx.num(-1))]),
                            ctx.pow(ctx.num(2), ctx.num(-1)),
                        ),
                    )
                };
                let exp_part = if matches!(alpha.node(), AtomNode::Num(0)) {
                    ctx.num(1)
                } else {
                    ctx.fun("exp", &[ctx.mul(&[alpha, x])])
                };
                Some((
                    ctx.mul(&[exp_part, ctx.fun("cos", &[ctx.mul(&[beta, x])])]),
                    ctx.mul(&[exp_part, ctx.fun("sin", &[ctx.mul(&[beta, x])])]),
                ))
            }
        }
        _ => None,
    }
}

/// Build the two fundamental solutions (y1, y2) of the homogeneous
/// Cauchy-Euler equation a*x^2*y'' + b*x*y' + c*y = 0.
///
/// `discriminant` is the discriminant of the indicial equation
/// a*r^2 + (b-a)*r + c = 0, and `b_ind` is the coefficient (b-a).
fn cauchy_euler_basis<'a>(
    ctx: &'a AtomArena<'a>,
    a: Atom<'a>,
    b_ind: Atom<'a>,
    discriminant: Atom<'a>,
    x: Atom<'a>,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match discriminant.node() {
        AtomNode::Num(d) => {
            if *d > 0 {
                let (r1, r2) = real_roots(ctx, a, b_ind, *d);
                Some((ctx.pow(x, r1), ctx.pow(x, r2)))
            } else if *d == 0 {
                let r = ctx.mul(&[
                    ctx.num(-1),
                    b_ind,
                    ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ]);
                let xr = ctx.pow(x, r);
                Some((xr, ctx.mul(&[xr, ctx.fun("log", &[x])])))
            } else {
                // Complex roots alpha ± beta*i:
                // y = x^alpha * cos/sin(beta * log(x))
                let alpha = ctx.mul(&[
                    ctx.num(-1),
                    b_ind,
                    ctx.pow(ctx.mul(&[ctx.num(2), a]), ctx.num(-1)),
                ]);
                let beta_sq = ctx.mul(&[ctx.num(-1), discriminant]);
                let two_a = ctx.mul(&[ctx.num(2), a]);
                let beta = ctx.pow(
                    ctx.mul(&[beta_sq, ctx.pow(ctx.pow(two_a, ctx.num(2)), ctx.num(-1))]),
                    ctx.pow(ctx.num(2), ctx.num(-1)),
                );
                let x_alpha = ctx.pow(x, alpha);
                let beta_logx = ctx.mul(&[beta, ctx.fun("log", &[x])]);
                Some((
                    ctx.mul(&[x_alpha, ctx.fun("cos", &[beta_logx])]),
                    ctx.mul(&[x_alpha, ctx.fun("sin", &[beta_logx])]),
                ))
            }
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Reduction of order
// ---------------------------------------------------------------------------

/// Attempt reduction of order for a second-order linear ODE
/// $a(x)\,y'' + b(x)\,y' + c(x)\,y = f(x)$.
///
/// Tries simple candidate solutions ($1$, $x$, $x^2$, $e^x$, $e^{-x}$,
/// $e^{2x}$) for the homogeneous equation; when one is found, the second
/// solution is built as
/// $y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$ with $p = b/a$.
/// A particular solution for nonzero $f$ is obtained via variation of
/// parameters.
pub(crate) fn solve_reduction_of_order<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;

    if ode_order(equation, func, var) != 2 || !super::util::is_linear_in(equation, func, var) {
        return None;
    }

    let x = ctx.var(var.as_str());
    let (a, b, _c, forcing) = extract_second_order_coeffs(ctx, equation, func, var)?;
    if is_atom_zero(a) {
        return None;
    }

    // Standard form: p = b/a (needed for the reduction-of-order integral).
    let p = super::util::collect_terms(ctx, ctx.mul(&[b, ctx.pow(a, ctx.num(-1))]));

    // Candidate first solutions for the homogeneous equation.
    let candidates = [
        ctx.num(1),
        x,
        ctx.pow(x, ctx.num(2)),
        ctx.fun("exp", &[x]),
        ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), x])]),
        ctx.fun("exp", &[ctx.mul(&[ctx.num(2), x])]),
    ];

    for y1 in candidates {
        // Check that a*y1'' + b*y1' + c*y1 == 0 using the extracted coeffs.
        if !satisfies_extracted(ctx, y1, a, b, _c, var) {
            continue;
        }

        // y2 = y1 * ∫( e^(-∫p) / y1^2 ) dx.
        let p_int = integrate(ctx, p, var);
        let exp_neg_p = super::util::exp_simplify(ctx, ctx.mul(&[ctx.num(-1), p_int]));
        let integrand = super::util::collect_terms(
            ctx,
            ctx.mul(&[exp_neg_p, ctx.pow(ctx.pow(y1, ctx.num(2)), ctx.num(-1))]),
        );
        let u = integrate(ctx, integrand, var);
        if is_integral_fallback(u) {
            continue;
        }
        let y2 = super::util::collect_terms(ctx, ctx.mul(&[y1, u]));
        if y2.to_string() == y1.to_string() {
            continue; // degenerate: need an independent second solution
        }

        let c1 = ctx.var("C1");
        let c2 = ctx.var("C2");
        let y_c = ctx.add(&[ctx.mul(&[c1, y1]), ctx.mul(&[c2, y2])]);

        // Particular solution via VOP when forcing is present.
        let y_p = if is_atom_zero(forcing) {
            ctx.num(0)
        } else {
            let g = super::util::collect_terms(ctx, ctx.mul(&[forcing, ctx.pow(a, ctx.num(-1))]));
            variation_of_parameters(ctx, y1, y2, g, var).unwrap_or(ctx.num(0))
        };

        let rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
        let solution = simplify(ctx, ctx.add(&[y_c, y_p]), &rules, 20);
        return Some(ODESolution::Explicit(solution));
    }

    None
}

/// Check whether `y1` satisfies a*y'' + b*y' + c*y = 0.
fn satisfies_extracted<'a>(
    ctx: &'a AtomArena<'a>,
    y1: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    var: Symbol,
) -> bool {
    let dy1 = diff(ctx, y1, var);
    let d2y1 = diff(ctx, dy1, var);
    let residual = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.mul(&[a, d2y1]), ctx.mul(&[b, dy1]), ctx.mul(&[c, y1])]),
    );
    matches!(residual.node(), AtomNode::Num(0))
}

// ---------------------------------------------------------------------------
// Variation of parameters
// ---------------------------------------------------------------------------

/// Variation of parameters for the standard-form equation
/// $y'' + p(x)\,y' + q(x)\,y = g(x)$ given two fundamental solutions
/// of the homogeneous equation.
///
/// Returns $y_p = -y_1 \int \frac{y_2\,g}{W}\,dx + y_2 \int \frac{y_1\,g}{W}\,dx$
/// where $W = y_1 y_2' - y_1' y_2$ is the Wronskian, or `None` if the
/// integrals cannot be evaluated in closed form.
fn variation_of_parameters<'a>(
    ctx: &'a AtomArena<'a>,
    y1: Atom<'a>,
    y2: Atom<'a>,
    g: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let dy1 = diff(ctx, y1, var);
    let dy2 = diff(ctx, y2, var);

    // Wronskian W = y1*y2' - y1'*y2.
    let w = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.mul(&[y1, dy2]), ctx.mul(&[ctx.num(-1), dy1, y2])]),
    );
    if matches!(w.node(), AtomNode::Num(0)) {
        return None;
    }
    let w_inv = ctx.pow(w, ctx.num(-1));

    // u1' = -y2*g/W, u2' = y1*g/W.
    let u1_prime = super::util::collect_terms(ctx, ctx.mul(&[ctx.num(-1), y2, g, w_inv]));
    let u2_prime = super::util::collect_terms(ctx, ctx.mul(&[y1, g, w_inv]));

    let u1 = integrate(ctx, u1_prime, var);
    if is_integral_fallback(u1) {
        return None;
    }
    let u2 = integrate(ctx, u2_prime, var);
    if is_integral_fallback(u2) {
        return None;
    }

    let y_p = ctx.add(&[ctx.mul(&[y1, u1]), ctx.mul(&[y2, u2])]);
    Some(super::util::collect_terms(ctx, y_p))
}

/// Check if an integration result is the unevaluated `Integral(...)` form.
fn is_integral_fallback<'a>(expr: Atom<'a>) -> bool {
    matches!(expr.node(), AtomNode::Fun(name, _) if name.as_str() == "Integral")
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
    let c1 = ctx.var("C1");
    let c2 = ctx.var("C2");
    if let Some((y1, y2)) = constant_coeff_basis(ctx, a, b, discriminant, x) {
        return ctx.add(&[ctx.mul(&[c1, y1]), ctx.mul(&[c2, y2])]);
    }

    // Discriminant is not a simple number — symbolic roots.
    // Return a generic form.
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
    // Flatten nested products and fold numeric signs first so that factor
    // extraction sees a canonical forcing term.
    let forcing = super::util::collect_terms(ctx, forcing);

    // Superposition: for a sum forcing, solve each term separately.
    if let AtomNode::Add(args) = forcing.node()
        && args.len() > 1
    {
        let mut parts = Vec::with_capacity(args.len());
        for term in args.iter() {
            parts.push(particular_solution_undetermined(
                ctx, *term, a, b, c, x, var,
            )?);
        }
        return Some(super::util::collect_terms(ctx, ctx.add(&parts)));
    }

    // Try polynomial forcing first.
    if is_polynomial_in(forcing, var) {
        let degree = polynomial_degree(forcing, var);
        return undetermined_polynomial(ctx, forcing, a, b, c, x, degree, var);
    }

    // Try exponential forcing: F*exp(kx).
    if let Some((f_coeff, k)) = extract_exp_forcing(ctx, forcing, var) {
        return undetermined_exponential(ctx, f_coeff, k, a, b, c, x);
    }

    // Try trigonometric forcing: f_c*cos(wx) + f_s*sin(wx).
    if let Some((f_c, f_s, w)) = extract_trig_forcing(ctx, forcing, var) {
        return undetermined_trig(ctx, f_c, f_s, w, a, b, c, x);
    }

    None
}

/// Undetermined coefficients for polynomial forcing of any degree.
///
/// Solves the triangular coefficient system by back-substitution. Handles
/// resonance (zero as a characteristic root with multiplicity s) by
/// multiplying the ansatz by x^s.
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
    // Extract forcing coefficients f_0..f_degree.
    let f = polynomial_coeffs(ctx, forcing, var, degree)?;

    // Resonance shift: multiplicity of 0 as a characteristic root.
    let c_zero = is_zero_atom(ctx, c);
    let b_zero = is_zero_atom(ctx, b);
    let s: i64 = if c_zero {
        if b_zero { 2 } else { 1 }
    } else {
        0
    };

    // Back-substitution. y_p = sum_{k=0}^{degree} A_k x^{k+s} with
    // A_{d+1} = A_{d+2} = 0.
    let mut coeffs: Vec<Atom<'a>> = vec![ctx.num(0); degree + 2];
    coeffs.push(ctx.num(0)); // A_{d+2} slot
    for k in (0..=degree).rev() {
        let k1 = ctx.num(k as i64 + 1);
        let k2 = ctx.num(k as i64 + 2);
        let a_next = coeffs[k + 1];
        let a_next2 = coeffs[k + 2];
        let numerator = match s {
            0 => {
                // f_k - b*(k+1)*A_{k+1} - a*(k+2)*(k+1)*A_{k+2}
                ctx.add(&[
                    f[k],
                    ctx.mul(&[ctx.num(-1), b, k1, a_next]),
                    ctx.mul(&[ctx.num(-1), a, k2, k1, a_next2]),
                ])
            }
            1 => {
                // f_k - a*(k+2)*(k+1)*A_{k+1}
                ctx.add(&[f[k], ctx.mul(&[ctx.num(-1), a, k2, k1, a_next])])
            }
            _ => f[k],
        };
        let denominator = match s {
            0 => c,
            1 => ctx.mul(&[b, k1]),
            _ => ctx.mul(&[a, k2, k1]),
        };
        let ak = super::util::collect_terms(
            ctx,
            ctx.mul(&[numerator, ctx.pow(denominator, ctx.num(-1))]),
        );
        coeffs[k] = ak;
    }

    // y_p = sum A_k x^{k+s}
    let mut terms = Vec::with_capacity(degree + 1);
    for (k, ak) in coeffs.iter().take(degree + 1).enumerate() {
        let power = ctx.pow(x, ctx.num(k as i64 + s));
        terms.push(ctx.mul(&[*ak, power]));
    }
    Some(super::util::collect_terms(ctx, ctx.add(&terms)))
}

/// Extract coefficients f_0..=f_degree of a polynomial in `var`.
///
/// Coefficients may be arbitrary x-free atoms (symbolic constants).
fn polynomial_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    poly: Atom<'a>,
    var: Symbol,
    degree: usize,
) -> Option<Vec<Atom<'a>>> {
    let collected = super::util::collect_terms(ctx, poly);
    let terms: Vec<Atom<'a>> = match collected.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![collected],
    };
    let mut coeffs = vec![ctx.num(0); degree + 1];
    for term in terms {
        // Split each term into coefficient (x-free) and x-power.
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let mut deg: usize = 0;
        let mut coeff_factors: Vec<Atom<'a>> = Vec::new();
        let mut ok = true;
        for f in factors {
            match f.node() {
                AtomNode::Var(v) if *v == var => deg += 1,
                AtomNode::Pow(base, exp) => {
                    if let (AtomNode::Var(v), AtomNode::Num(n)) = (base.node(), exp.node())
                        && *v == var
                        && *n >= 0
                    {
                        deg += *n as usize;
                    } else if contains_x(f, var) {
                        ok = false;
                        break;
                    } else {
                        coeff_factors.push(f);
                    }
                }
                _ => {
                    if contains_x(f, var) {
                        ok = false;
                        break;
                    }
                    coeff_factors.push(f);
                }
            }
        }
        if !ok || deg > degree {
            return None;
        }
        let coeff = if coeff_factors.is_empty() {
            ctx.num(1)
        } else {
            ctx.mul(&coeff_factors)
        };
        coeffs[deg] = super::util::collect_terms(ctx, ctx.add(&[coeffs[deg], coeff]));
    }
    Some(coeffs)
}

/// Check whether an atom simplifies to zero.
fn is_zero_atom<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> bool {
    matches!(
        super::util::collect_terms(ctx, expr).node(),
        AtomNode::Num(0)
    )
}

/// Undetermined coefficients for exponential forcing F*exp(kx),
/// including the resonance cases (k a single or double characteristic root).
fn undetermined_exponential<'a>(
    ctx: &'a AtomArena<'a>,
    f_coeff: Atom<'a>,
    k: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let exp_kx = ctx.fun("exp", &[ctx.mul(&[k, x])]);

    // char(k) = a*k^2 + b*k + c
    let char_at_k = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.mul(&[a, ctx.pow(k, ctx.num(2))]), ctx.mul(&[b, k]), c]),
    );
    if !matches!(char_at_k.node(), AtomNode::Num(0)) {
        // No resonance: B = F / char(k).
        let coeff =
            super::util::collect_terms(ctx, ctx.mul(&[f_coeff, ctx.pow(char_at_k, ctx.num(-1))]));
        return Some(ctx.mul(&[coeff, exp_kx]));
    }

    // char'(k) = 2*a*k + b: single root if nonzero, double root if zero.
    let char_prime = super::util::collect_terms(ctx, ctx.add(&[ctx.mul(&[ctx.num(2), a, k]), b]));
    if !matches!(char_prime.node(), AtomNode::Num(0)) {
        // Single root: y_p = F*x*exp(kx) / char'(k).
        let coeff =
            super::util::collect_terms(ctx, ctx.mul(&[f_coeff, ctx.pow(char_prime, ctx.num(-1))]));
        return Some(ctx.mul(&[coeff, x, exp_kx]));
    }

    // Double root: y_p = F*x^2*exp(kx) / (2a).
    let two_a = ctx.mul(&[ctx.num(2), a]);
    let coeff = super::util::collect_terms(ctx, ctx.mul(&[f_coeff, ctx.pow(two_a, ctx.num(-1))]));
    Some(ctx.mul(&[coeff, ctx.pow(x, ctx.num(2)), exp_kx]))
}

/// Undetermined coefficients for trigonometric forcing
/// f_c*cos(wx) + f_s*sin(wx).
#[allow(clippy::too_many_arguments)]
fn undetermined_trig<'a>(
    ctx: &'a AtomArena<'a>,
    f_c: Atom<'a>,
    f_s: Atom<'a>,
    w: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    c: Atom<'a>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let wx = ctx.mul(&[w, x]);
    let cos_wx = ctx.fun("cos", &[wx]);
    let sin_wx = ctx.fun("sin", &[wx]);

    // System matrix entries: p = c - a*w^2, q = b*w.
    let w2 = ctx.pow(w, ctx.num(2));
    let p = super::util::collect_terms(ctx, ctx.add(&[c, ctx.mul(&[ctx.num(-1), a, w2])]));
    let q = super::util::collect_terms(ctx, ctx.mul(&[b, w]));

    let det = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.pow(p, ctx.num(2)), ctx.pow(q, ctx.num(2))]),
    );

    if matches!(det.node(), AtomNode::Num(0)) {
        // Resonance: p = q = 0, i.e. b=0 and c = a*w^2.
        // y_p = x/(2w) * (f_c/a * sin(wx) - f_s/a * cos(wx)).
        if !matches!(q.node(), AtomNode::Num(0)) || !matches!(p.node(), AtomNode::Num(0)) {
            return None;
        }
        let two_w_a = ctx.mul(&[ctx.num(2), w, a]);
        let scale = super::util::collect_terms(ctx, ctx.pow(two_w_a, ctx.num(-1)));
        return Some(super::util::collect_terms(
            ctx,
            ctx.mul(&[
                scale,
                x,
                ctx.add(&[
                    ctx.mul(&[f_c, sin_wx]),
                    ctx.mul(&[ctx.num(-1), f_s, cos_wx]),
                ]),
            ]),
        ));
    }

    // Cramer for p*A + q*B = f_c, -q*A + p*B = f_s with det = p^2 + q^2:
    // A = (f_c*p - f_s*q)/det, B = (f_c*q + f_s*p)/det.
    let det_inv = ctx.pow(det, ctx.num(-1));
    let big_a = super::util::collect_terms(
        ctx,
        ctx.mul(&[
            ctx.add(&[ctx.mul(&[f_c, p]), ctx.mul(&[ctx.num(-1), f_s, q])]),
            det_inv,
        ]),
    );
    let big_b = super::util::collect_terms(
        ctx,
        ctx.mul(&[ctx.add(&[ctx.mul(&[f_c, q]), ctx.mul(&[f_s, p])]), det_inv]),
    );
    Some(super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.mul(&[big_a, cos_wx]), ctx.mul(&[big_b, sin_wx])]),
    ))
}

/// Extract (f_c, f_s, w) from forcing of the form
/// f_c*cos(wx) + f_s*sin(wx), where f_c, f_s, w are x-free atoms.
/// Either f_c or f_s may be zero (checked by the caller via Num(0)).
fn extract_trig_forcing<'a>(
    ctx: &'a AtomArena<'a>,
    forcing: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    let collected = super::util::collect_terms(ctx, forcing);
    let terms: Vec<Atom<'a>> = match collected.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![collected],
    };
    let mut f_c = ctx.num(0);
    let mut f_s = ctx.num(0);
    let mut w: Option<Atom<'a>> = None;
    for term in terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let mut trig: Option<(&str, Atom<'a>)> = None;
        let mut coeff_factors: Vec<Atom<'a>> = Vec::new();
        for f in factors {
            match f.node() {
                AtomNode::Fun(name, fargs)
                    if (name.as_str() == "cos" || name.as_str() == "sin") && fargs.len() == 1 =>
                {
                    if trig.is_some() {
                        return None; // product of two trig factors
                    }
                    trig = Some((name.as_str(), fargs[0]));
                }
                _ => {
                    if contains_x(f, var) {
                        return None;
                    }
                    coeff_factors.push(f);
                }
            }
        }
        let (kind, arg) = trig?;
        // The trig argument must be w*x with w x-free.
        let w_term = match arg.node() {
            AtomNode::Var(v) if *v == var => ctx.num(1),
            AtomNode::Mul(margs) => {
                let non_x: Vec<_> = margs
                    .iter()
                    .filter(|f| !matches!(f.node(), AtomNode::Var(v) if *v == var))
                    .copied()
                    .collect();
                let x_count = margs
                    .iter()
                    .filter(|f| matches!(f.node(), AtomNode::Var(v) if *v == var))
                    .count();
                if x_count != 1 || non_x.iter().any(|f| contains_x(*f, var)) {
                    return None;
                }
                if non_x.is_empty() {
                    ctx.num(1)
                } else {
                    ctx.mul(&non_x)
                }
            }
            _ => return None,
        };
        match w {
            None => w = Some(w_term),
            Some(existing) if existing.to_string() == w_term.to_string() => {}
            _ => return None, // different frequencies
        }
        let coeff = if coeff_factors.is_empty() {
            ctx.num(1)
        } else {
            ctx.mul(&coeff_factors)
        };
        if kind == "cos" {
            f_c = coeff;
        } else {
            f_s = coeff;
        }
    }
    w.map(|w| (f_c, f_s, w))
}

/// Extract (F, k) from forcing of the form F*exp(k*x), F x-free.
fn extract_exp_forcing<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut k: Option<Atom<'a>> = None;
    let mut coeff_factors: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        match f.node() {
            AtomNode::Fun(name, fargs) if name.as_str() == "exp" && fargs.len() == 1 => {
                if k.is_some() {
                    return None;
                }
                k = Some(match extract_exp_k_inner(fargs[0], var)? {
                    Some(kk) => kk,
                    None => ctx.num(1),
                });
            }
            _ => {
                if contains_x(f, var) {
                    return None;
                }
                coeff_factors.push(f);
            }
        }
    }
    let k = k?;
    let f_coeff = if coeff_factors.is_empty() {
        ctx.num(1)
    } else {
        ctx.mul(&coeff_factors)
    };
    Some((f_coeff, k))
}

// ---------------------------------------------------------------------------
// Coefficient extraction helpers
// ---------------------------------------------------------------------------

/// Extract a, b, c coefficients from a*y'' + b*y' + c*y + ... = 0.
pub(crate) fn extract_second_order_coeffs<'a>(
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
    // If the discriminant is a perfect square, use the exact integer root;
    // otherwise keep sqrt(disc) in symbolic form so that, e.g., disc=2
    // yields (-b ± sqrt(2))/(2a) instead of the truncated isqrt value.
    let s = isqrt(disc);
    let sqrt_disc = if s * s == disc {
        ctx.num(s)
    } else {
        ctx.pow(ctx.num(disc), ctx.pow(ctx.num(2), ctx.num(-1)))
    };
    let two_a = ctx.mul(&[ctx.num(2), a]);
    let neg_b = ctx.mul(&[ctx.num(-1), b]);
    let r1 = super::util::collect_terms(
        ctx,
        ctx.mul(&[ctx.add(&[neg_b, sqrt_disc]), ctx.pow(two_a, ctx.num(-1))]),
    );
    let r2 = super::util::collect_terms(
        ctx,
        ctx.mul(&[
            ctx.add(&[neg_b, ctx.mul(&[ctx.num(-1), sqrt_disc])]),
            ctx.pow(two_a, ctx.num(-1)),
        ]),
    );
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

/// Extract k from the argument of exp(...) where the argument is k*x.
/// Returns k=1 implicitly via `Some(None)` when the argument is bare x.
fn extract_exp_k_inner<'a>(arg: Atom<'a>, var: Symbol) -> Option<Option<Atom<'a>>> {
    match arg.node() {
        AtomNode::Var(v) if *v == var => Some(None),
        AtomNode::Mul(factors) => {
            let k: Vec<_> = factors
                .iter()
                .filter(|f| !matches!(f.node(), AtomNode::Var(v) if *v == var))
                .copied()
                .collect();
            let x_count = factors
                .iter()
                .filter(|f| matches!(f.node(), AtomNode::Var(v) if *v == var))
                .count();
            if k.len() == 1 && x_count == 1 && !contains_x(k[0], var) {
                return Some(Some(k[0]));
            }
            None
        }
        _ => None,
    }
}
