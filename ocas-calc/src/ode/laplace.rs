//! Laplace transform methods for linear ODE initial value problems.
//!
//! Provides [`solve_laplace`] for first- and second-order linear
//! constant-coefficient ODEs with given initial conditions.
//!
//! The method transforms the ODE into an algebraic equation in $Y(s)$,
//! solves for $Y(s)$, and inverts the transform using a table of standard
//! pairs: exponentials, sines, cosines (optionally damped), powers of $x$,
//! and sums thereof via linearity and simple partial fractions.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::ODESolution;
use super::second_order::extract_second_order_coeffs;
use super::util::{collect_terms, is_linear_in, ode_order};

/// Solve a linear constant-coefficient IVP by the Laplace transform method.
///
/// Handles $a y'' + b y' + c y = f(x)$ and $a y' + b y = f(x)$ with
/// initial conditions $y(0) = y_0$ (and $y'(0) = y_1$ for second order),
/// when the forcing $f$ has a table Laplace transform (polynomial,
/// exponential, sine, cosine, and exponential times sine/cosine).
///
/// The result is an explicit solution containing no free constants.
pub(crate) fn solve_laplace<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
    y0: Atom<'a>,
    y1: Option<Atom<'a>>,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;
    let order = ode_order(equation, func, var);
    if !(1..=2).contains(&order) || !is_linear_in(equation, func, var) {
        return None;
    }

    let (a_num, b_num, c_num, forcing) = if order == 1 {
        // p*y' + q*y + g(x) = 0  =>  forcing = -g(x)
        let (p, q, g) = extract_first_order_pqg(ctx, equation, func, var)?;
        let (Some(pn), Some(qn)) = (
            const_i64(collect_terms(ctx, p)),
            const_i64(collect_terms(ctx, q)),
        ) else {
            return None;
        };
        (pn, qn, 0i64, ctx.mul(&[ctx.num(-1), g]))
    } else {
        let (a, b, c, forcing) = extract_second_order_coeffs(ctx, equation, func, var)?;
        let (Some(an), Some(bn), Some(cn)) = (
            const_i64(collect_terms(ctx, a)),
            const_i64(collect_terms(ctx, b)),
            const_i64(collect_terms(ctx, c)),
        ) else {
            return None;
        };
        (an, bn, cn, forcing)
    };

    let f_s = laplace_transform(ctx, forcing, var)?;

    // Build Y(s) from the transformed equation.
    let s = ctx.var("s");
    let y1_val = y1.unwrap_or_else(|| ctx.num(0));
    let (num_s, den_s) = if order == 1 {
        // a*y' + b*y = f  =>  Y = (F + a*y0) / (a*s + b)
        let num = ctx.add(&[f_s, ctx.mul(&[ctx.num(a_num), y0])]);
        let den = ctx.add(&[ctx.mul(&[ctx.num(a_num), s]), ctx.num(b_num)]);
        (num, den)
    } else {
        // a*y'' + b*y' + c*y = f
        // Y = (F + a*(s*y0 + y1) + b*y0) / (a*s^2 + b*s + c)
        let num = ctx.add(&[
            f_s,
            ctx.mul(&[ctx.num(a_num), ctx.add(&[ctx.mul(&[s, y0]), y1_val])]),
            ctx.mul(&[ctx.num(b_num), y0]),
        ]);
        let den = ctx.add(&[
            ctx.mul(&[ctx.num(a_num), ctx.pow(s, ctx.num(2))]),
            ctx.mul(&[ctx.num(b_num), s]),
            ctx.num(c_num),
        ]);
        (num, den)
    };

    let y_s = collect_terms(ctx, ctx.mul(&[num_s, ctx.pow(den_s, ctx.num(-1))]));
    let y_x = inverse_laplace(ctx, y_s, Symbol::new("s"), var)?;

    Some(ODESolution::Explicit(collect_terms(ctx, y_x)))
}

// ---------------------------------------------------------------------------
// Forward Laplace transform
// ---------------------------------------------------------------------------

/// Laplace transform of `expr` with respect to `var` (x -> s).
fn laplace_transform<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let expr = collect_terms(ctx, expr);
    let terms: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut parts = Vec::with_capacity(terms.len());
    for term in terms {
        parts.push(laplace_term(ctx, term, var)?);
    }
    Some(collect_terms(ctx, ctx.add(&parts)))
}

/// Transform a single term (x-free constant times one kernel).
fn laplace_term<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let (coeff, kernel) = split_const_factors(ctx, term, var);
    let transformed = laplace_kernel(ctx, kernel, var)?;
    Some(collect_terms(ctx, ctx.mul(&[coeff, transformed])))
}

/// Transform the x-dependent kernel.
fn laplace_kernel<'a>(ctx: &'a AtomArena<'a>, kernel: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let s = ctx.var("s");
    match kernel.node() {
        AtomNode::Num(_) => Some(ctx.pow(s, ctx.num(-1))),
        AtomNode::Var(v) if *v == var => Some(ctx.pow(ctx.pow(s, ctx.num(2)), ctx.num(-1))),
        AtomNode::Pow(base, exp) => {
            if let (AtomNode::Var(v), AtomNode::Num(n)) = (base.node(), exp.node())
                && *v == var
                && *n >= 0
                && *n <= 12
            {
                let fact = (1..=*n).product::<i64>();
                return Some(ctx.mul(&[
                    ctx.num(fact),
                    ctx.pow(ctx.pow(s, ctx.num(*n + 1)), ctx.num(-1)),
                ]));
            }
            None
        }
        AtomNode::Fun(name, args) if name.as_str() == "exp" && args.len() == 1 => {
            let k = linear_coeff_of(ctx, args[0], var)?;
            Some(ctx.pow(ctx.add(&[s, ctx.mul(&[ctx.num(-1), k])]), ctx.num(-1)))
        }
        AtomNode::Fun(name, args)
            if (name.as_str() == "sin" || name.as_str() == "cos") && args.len() == 1 =>
        {
            let w = linear_coeff_of(ctx, args[0], var)?;
            let denom = ctx.add(&[ctx.pow(s, ctx.num(2)), ctx.pow(w, ctx.num(2))]);
            let numerator = if name.as_str() == "sin" { w } else { s };
            Some(ctx.mul(&[numerator, ctx.pow(denom, ctx.num(-1))]))
        }
        AtomNode::Mul(factors) => {
            // exp(kx) * sin(wx) or exp(kx) * cos(wx)
            let mut k: Option<Atom<'a>> = None;
            let mut trig: Option<(&str, Atom<'a>)> = None;
            for f in factors.iter() {
                match f.node() {
                    AtomNode::Fun(name, args) if name.as_str() == "exp" && args.len() == 1 => {
                        k = Some(linear_coeff_of(ctx, args[0], var)?);
                    }
                    AtomNode::Fun(name, args)
                        if (name.as_str() == "sin" || name.as_str() == "cos")
                            && args.len() == 1 =>
                    {
                        trig = Some((name.as_str(), linear_coeff_of(ctx, args[0], var)?));
                    }
                    _ => return None,
                }
            }
            let (k, (kind, w)) = (k?, trig?);
            let s_minus_k = ctx.add(&[s, ctx.mul(&[ctx.num(-1), k])]);
            let denom = ctx.add(&[ctx.pow(s_minus_k, ctx.num(2)), ctx.pow(w, ctx.num(2))]);
            let numerator = if kind == "sin" { w } else { s_minus_k };
            Some(ctx.mul(&[numerator, ctx.pow(denom, ctx.num(-1))]))
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Inverse Laplace transform
// ---------------------------------------------------------------------------

/// Inverse Laplace transform of `expr` (s -> x) via linearity + standard pairs.
fn inverse_laplace<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    s_sym: Symbol,
    var: Symbol,
) -> Option<Atom<'a>> {
    let expr = collect_terms(ctx, expr);
    let terms: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut parts = Vec::with_capacity(terms.len());
    for term in terms {
        parts.push(inverse_laplace_term(ctx, term, s_sym, var)?);
    }
    Some(collect_terms(ctx, ctx.add(&parts)))
}

/// Inverse-transform a single rational term in s.
fn inverse_laplace_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    s_sym: Symbol,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (num_s, den_s) = split_fraction(ctx, term, s_sym)?;
    let (n0, n1) = linear_coeffs_i64(ctx, num_s, s_sym)?;
    let (ca, cb, cc) = quadratic_coeffs_i64(ctx, den_s, s_sym)?;

    if ca == 0 {
        // First-order denominator b*s + c with constant numerator:
        // n0/(b*s + c) = (n0/b) / (s + c/b) -> (n0/b) * e^{-(c/b) x}
        if cb == 0 || n1 != 0 || n0 % cb != 0 || cc % cb != 0 {
            return None;
        }
        let scale = n0 / cb;
        let k = -(cc / cb);
        let x = ctx.var(var.as_str());
        return Some(collect_terms(
            ctx,
            ctx.mul(&[ctx.num(scale), ctx.fun("exp", &[ctx.mul(&[ctx.num(k), x])])]),
        ));
    }

    // Quadratic denominator a*s^2 + b*s + c.
    let disc = cb * cb - 4 * ca * cc;
    let den2 = 2 * ca;
    if disc > 0 {
        let sd = isqrt_i64(disc);
        if sd * sd == disc && (-cb + sd) % den2 == 0 && (-cb - sd) % den2 == 0 {
            let r1 = (-cb + sd) / den2;
            let r2 = (-cb - sd) / den2;
            if r1 != r2 {
                return Some(inverse_distinct_roots(ctx, var, r1, r2, n0, n1));
            }
        }
        return None;
    }
    if disc == 0 {
        // Repeated root r = -b/(2a).
        if (-cb) % den2 == 0 {
            let r = -cb / den2;
            return Some(inverse_repeated_root(ctx, var, r, n0, n1));
        }
        return None;
    }
    // Complex roots k ± i*w with k = -b/(2a), w = sqrt(-disc)/(2a).
    let disc_neg = -disc;
    let sw = isqrt_i64(disc_neg);
    if (-cb) % den2 == 0 && sw * sw == disc_neg && sw % den2 == 0 {
        let k = -cb / den2;
        let w = sw / den2;
        return inverse_complex_roots(ctx, var, k, w, n0, n1);
    }
    None
}

/// Inverse for distinct real roots r1 != r2:
/// (n1 s + n0)/((s-r1)(s-r2)) = A/(s-r1) + B/(s-r2),
/// A = (n1 r1 + n0)/(r1 - r2), B = (n1 r2 + n0)/(r2 - r1).
fn inverse_distinct_roots<'a>(
    ctx: &'a AtomArena<'a>,
    var: Symbol,
    r1: i64,
    r2: i64,
    n0: i64,
    n1: i64,
) -> Atom<'a> {
    let x = ctx.var(var.as_str());
    let a_num = n1 * r1 + n0;
    let a_den = r1 - r2;
    let b_num = n1 * r2 + n0;
    let b_den = r2 - r1;
    let term1 = ctx.mul(&[
        rat_atom(ctx, a_num, a_den),
        ctx.fun("exp", &[ctx.mul(&[ctx.num(r1), x])]),
    ]);
    let term2 = ctx.mul(&[
        rat_atom(ctx, b_num, b_den),
        ctx.fun("exp", &[ctx.mul(&[ctx.num(r2), x])]),
    ]);
    collect_terms(ctx, ctx.add(&[term1, term2]))
}

/// Inverse for repeated root r:
/// (n1 s + n0)/(s-r)^2 = n1/(s-r) + (n1 r + n0)/(s-r)^2
/// -> n1 e^{rx} + (n1 r + n0) x e^{rx}.
fn inverse_repeated_root<'a>(
    ctx: &'a AtomArena<'a>,
    var: Symbol,
    r: i64,
    n0: i64,
    n1: i64,
) -> Atom<'a> {
    let x = ctx.var(var.as_str());
    let e_rx = ctx.fun("exp", &[ctx.mul(&[ctx.num(r), x])]);
    let term1 = ctx.mul(&[ctx.num(n1), e_rx]);
    let term2 = ctx.mul(&[ctx.num(n1 * r + n0), x, e_rx]);
    collect_terms(ctx, ctx.add(&[term1, term2]))
}

/// Inverse for complex roots k ± i*w:
/// (n1 s + n0)/((s-k)^2 + w^2) = e^{kx} [ n1 cos(wx) + ((n0 + k n1)/w) sin(wx) ].
fn inverse_complex_roots<'a>(
    ctx: &'a AtomArena<'a>,
    var: Symbol,
    k: i64,
    w: i64,
    n0: i64,
    n1: i64,
) -> Option<Atom<'a>> {
    if w == 0 || (n0 + k * n1) % w != 0 {
        return None;
    }
    let x = ctx.var(var.as_str());
    let cos_wx = ctx.fun("cos", &[ctx.mul(&[ctx.num(w), x])]);
    let sin_wx = ctx.fun("sin", &[ctx.mul(&[ctx.num(w), x])]);
    let cos_part = ctx.mul(&[ctx.num(n1), cos_wx]);
    let sin_part = ctx.mul(&[ctx.num((n0 + k * n1) / w), sin_wx]);
    let bracket = ctx.add(&[cos_part, sin_part]);
    // k = 0 means no damping: drop the exp factor entirely.
    if k == 0 {
        return Some(collect_terms(ctx, bracket));
    }
    let e_kx = ctx.fun("exp", &[ctx.mul(&[ctx.num(k), x])]);
    Some(collect_terms(ctx, ctx.mul(&[e_kx, bracket])))
}

// ---------------------------------------------------------------------------
// Term decomposition helpers
// ---------------------------------------------------------------------------

/// Split a term into (numerator_in_s, denominator_in_s).
/// Non-numeric s-free factors are rejected (keeps integer arithmetic).
fn split_fraction<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    s_sym: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut num_factors: Vec<Atom<'a>> = Vec::new();
    let mut den: Option<Atom<'a>> = None;
    for f in factors {
        if let AtomNode::Pow(base, exp) = f.node()
            && let AtomNode::Num(n) = exp.node()
            && *n < 0
            && contains_var(*base, s_sym)
        {
            if den.is_some() {
                return None;
            }
            den = Some(ctx.pow(*base, ctx.num(-*n)));
            continue;
        }
        if !contains_var(f, s_sym) && !matches!(f.node(), AtomNode::Num(_)) {
            // Non-numeric s-free factor (symbolic constant): reject.
            return None;
        }
        num_factors.push(f);
    }
    let num = if num_factors.is_empty() {
        ctx.num(1)
    } else {
        ctx.mul(&num_factors)
    };
    let den = den?;
    Some((num, den))
}

/// Extract integer (A, B, C) from A*s^2 + B*s + C (any subset of terms).
fn quadratic_coeffs_i64<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    s_sym: Symbol,
) -> Option<(i64, i64, i64)> {
    let expr = collect_terms(ctx, expr);
    let terms: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    };
    let (mut ca, mut cb, mut cc) = (0i64, 0i64, 0i64);
    for t in terms {
        let (deg, num) = poly_term_deg(t, s_sym)?;
        match deg {
            2 => ca += num,
            1 => cb += num,
            0 => cc += num,
            _ => return None,
        }
    }
    Some((ca, cb, cc))
}

/// Extract (n0, n1) from a constant-or-linear expression n0 + n1*s.
fn linear_coeffs_i64<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    s_sym: Symbol,
) -> Option<(i64, i64)> {
    let (ca, cb, cc) = quadratic_coeffs_i64(ctx, expr, s_sym)?;
    if ca != 0 {
        return None;
    }
    Some((cc, cb))
}

/// Determine (degree, numeric_coeff) of a single polynomial term in s.
fn poly_term_deg<'a>(term: Atom<'a>, s_sym: Symbol) -> Option<(i64, i64)> {
    match term.node() {
        AtomNode::Num(n) => Some((0, *n)),
        AtomNode::Var(v) if *v == s_sym => Some((1, 1)),
        AtomNode::Pow(base, exp) => {
            if let (AtomNode::Var(v), AtomNode::Num(k)) = (base.node(), exp.node())
                && *v == s_sym
            {
                Some((*k, 1))
            } else {
                None
            }
        }
        AtomNode::Mul(factors) => {
            let mut num: i64 = 1;
            let mut deg: i64 = 0;
            for f in factors.iter() {
                match f.node() {
                    AtomNode::Num(n) => num *= n,
                    AtomNode::Var(v) if *v == s_sym => deg += 1,
                    AtomNode::Pow(b, e) => {
                        if let (AtomNode::Var(v), AtomNode::Num(k)) = (b.node(), e.node())
                            && *v == s_sym
                        {
                            deg += k;
                        } else {
                            return None;
                        }
                    }
                    _ => return None,
                }
            }
            Some((deg, num))
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

fn const_i64<'a>(expr: Atom<'a>) -> Option<i64> {
    match expr.node() {
        AtomNode::Num(n) => Some(*n),
        _ => None,
    }
}

/// Extract (p, q, g) from a first-order linear equation
/// p*y' + q*y + g(x) = 0, where p, q are y-free and g contains no y/y'.
fn extract_first_order_pqg<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    let x = ctx.var(var.as_str());
    let dy = ctx.fun("Derivative", &[func, x]);
    let dy_str = dy.to_string();
    let func_str = func.to_string();

    let terms: Vec<Atom<'a>> = match equation.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![equation],
    };

    let mut p_terms: Vec<Atom<'a>> = Vec::new();
    let mut q_terms: Vec<Atom<'a>> = Vec::new();
    let mut g_terms: Vec<Atom<'a>> = Vec::new();

    for term in terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let has_dy = factors.iter().any(|f| f.to_string() == dy_str);
        let has_y = factors.iter().any(|f| f.to_string() == func_str);
        if has_dy {
            let rest: Vec<_> = factors
                .iter()
                .filter(|f| f.to_string() != dy_str)
                .copied()
                .collect();
            p_terms.push(if rest.is_empty() {
                ctx.num(1)
            } else {
                ctx.mul(&rest)
            });
        } else if has_y {
            let rest: Vec<_> = factors
                .iter()
                .filter(|f| f.to_string() != func_str)
                .copied()
                .collect();
            q_terms.push(if rest.is_empty() {
                ctx.num(1)
            } else {
                ctx.mul(&rest)
            });
        } else if term.to_string() == dy_str {
            p_terms.push(ctx.num(1));
        } else if term.to_string() == func_str {
            q_terms.push(ctx.num(1));
        } else if !super::util::contains_func(term, func, var) {
            g_terms.push(term);
        } else {
            return None;
        }
    }

    let p = if p_terms.is_empty() {
        ctx.num(0)
    } else {
        collect_terms(ctx, ctx.add(&p_terms))
    };
    let q = if q_terms.is_empty() {
        ctx.num(0)
    } else {
        collect_terms(ctx, ctx.add(&q_terms))
    };
    let g = if g_terms.is_empty() {
        ctx.num(0)
    } else {
        collect_terms(ctx, ctx.add(&g_terms))
    };
    Some((p, q, g))
}

fn rat_atom<'a>(ctx: &'a AtomArena<'a>, n: i64, d: i64) -> Atom<'a> {
    if d == 0 {
        return ctx.num(0);
    }
    let (n, d) = if d < 0 { (-n, -d) } else { (n, d) };
    if d == 1 {
        ctx.num(n)
    } else {
        ctx.mul(&[ctx.num(n), ctx.pow(ctx.num(d), ctx.num(-1))])
    }
}

fn isqrt_i64(n: i64) -> i64 {
    if n <= 0 {
        return 0;
    }
    let mut x = n;
    let mut y = x / 2 + 1;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

/// Split a term into (constant_coeff, kernel): kernel holds x-dependent factors.
fn split_const_factors<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
) -> (Atom<'a>, Atom<'a>) {
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut kernel: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        if contains_var(f, var) {
            kernel.push(f);
        } else {
            coeff.push(f);
        }
    }
    let coeff = if coeff.is_empty() {
        ctx.num(1)
    } else {
        ctx.mul(&coeff)
    };
    let kernel = if kernel.is_empty() {
        ctx.num(1)
    } else if kernel.len() == 1 {
        kernel[0]
    } else {
        ctx.mul(&kernel)
    };
    (coeff, kernel)
}

/// Extract k from k*x (linear in var, x-free k). Bare x yields k = 1.
fn linear_coeff_of<'a>(ctx: &'a AtomArena<'a>, arg: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match arg.node() {
        AtomNode::Var(v) if *v == var => Some(ctx.num(1)),
        AtomNode::Mul(factors) => {
            let mut x_count = 0;
            let mut rest: Vec<Atom<'a>> = Vec::new();
            for f in factors.iter() {
                if matches!(f.node(), AtomNode::Var(v) if *v == var) {
                    x_count += 1;
                } else if contains_var(*f, var) {
                    return None;
                } else {
                    rest.push(*f);
                }
            }
            if x_count == 1 && !rest.is_empty() {
                Some(if rest.len() == 1 {
                    rest[0]
                } else {
                    ctx.mul(&rest)
                })
            } else {
                None
            }
        }
        _ => None,
    }
}

fn contains_var<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => false,
        AtomNode::Var(v) => *v == var,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_var(*a, var)),
        AtomNode::Pow(base, exp) => contains_var(*base, var) || contains_var(*exp, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_var(*a, var)),
    }
}
