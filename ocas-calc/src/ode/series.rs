//! Power series and Frobenius methods for ODE solutions.
//!
//! Provides [`solve_power_series`] for computing Taylor-series solutions of
//! linear ODEs around ordinary points, and [`solve_frobenius`] for solutions
//! around regular singular points.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::ODESolution;
use super::util::ode_order;
use crate::derivative::diff;

/// Compute a power series solution of a linear ODE around `x0`.
///
/// Assumes `x0` is an **ordinary point** (all coefficients are analytic at `x0`).
/// Sets $y = \sum_{n=0}^{N} a_n (x - x_0)^n$, substitutes into the ODE, and
/// solves for the coefficients $a_n$ recursively.
///
/// Returns `ODESolution::Series(truncated_expr, n_terms)`.
///
/// # Limitations
///
/// - Only works for linear ODEs with polynomial or analytic coefficients.
/// - The recursion is symbolic; coefficients are left as rational expressions
///   of the initial conditions $a_0, a_1, \ldots, a_{k-1}$ (where $k$ is the
///   ODE order).
pub(crate) fn solve_power_series<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
    x0: Atom<'a>,
    n_terms: usize,
) -> Option<ODESolution<'a>> {
    let super::ODE {
        equation,
        func,
        var,
    } = ode;
    let order = ode_order(equation, func, var);
    if order == 0 {
        return None;
    }

    let x = ctx.var(var.as_str());
    let h = ctx.add(&[x, ctx.mul(&[ctx.num(-1), x0])]);

    // Build symbolic coefficients a0, a1, ... as plain variables so that
    // `diff` on the residual treats them as constants.
    let max_coeff = n_terms;
    let mut a_syms: Vec<Atom<'a>> = Vec::with_capacity(max_coeff);
    for i in 0..max_coeff {
        a_syms.push(ctx.var(&format!("a{i}")));
    }

    // Build y = sum a_n (x-x0)^n and its derivatives, treating a_n as
    // constants (using `diff` here would differentiate the a_n symbols as
    // unknown functions of x).
    let (y_series, y1, y2) = build_series_triple(ctx, &a_syms, h);

    // Substitute y, y', y'' into the ODE equation: residual R(x) must be
    // identically zero, so R(x0) = R'(x0) = R''(x0) = ... = 0. The k-th
    // condition is linear in the coefficients and determines the highest
    // remaining coefficient, giving a recurrence.
    let mut residual = substitute_series(ctx, equation, func, var, y_series, y1, y2);
    residual = super::util::collect_terms(ctx, residual);

    // Solved coefficient values (a_0 .. a_{order-1} stay free parameters).
    let mut solved: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();

    for k in 0..n_terms.saturating_sub(order) {
        // Evaluate the k-th derivative of the residual at x = x0.
        if k > 0 {
            residual = super::util::collect_terms(ctx, diff(ctx, residual, var));
        }
        let mut cond = super::classify::replace_atom(ctx, residual, x, x0);
        cond = super::util::collect_terms(ctx, cond);

        // Substitute already-solved coefficients.
        for (sym, val) in &solved {
            cond = super::classify::replace_atom(ctx, cond, *sym, *val);
        }
        cond = super::util::collect_terms(ctx, cond);

        // Solve for the highest-index coefficient that is not yet solved
        // and not a free parameter (a_0 .. a_{order-1}).
        let target = a_syms.iter().skip(order).rev().find(|s| {
            cond.to_string().contains(&s.to_string())
                && !solved
                    .iter()
                    .any(|(sym, _)| sym.to_string() == s.to_string())
        });
        match target {
            Some(&target) => {
                if let Some(value) = solve_linear_coeff(ctx, cond, target) {
                    solved.push((target, value));
                }
            }
            None => {
                // No new coefficient in this condition. A nontrivial
                // constraint on the free parameters means x0 is not an
                // ordinary point (e.g. a regular singular point), so the
                // Taylor ansatz is inconsistent — decline and let the
                // Frobenius method handle it.
                if !matches!(cond.node(), AtomNode::Num(0)) {
                    return None;
                }
            }
        }
    }

    if solved.is_empty() {
        return None;
    }

    // Rebuild the series with solved coefficients substituted.
    let mut series_expr = y_series;
    for (sym, val) in &solved {
        series_expr = super::classify::replace_atom(ctx, series_expr, *sym, *val);
    }
    series_expr = super::util::collect_terms(ctx, series_expr);

    Some(ODESolution::Series(series_expr, n_terms))
}

/// Solve `cond` of the form `coeff*target + rest = 0` for `target`.
///
/// Returns `-rest/coeff` with like terms collected, or `None` when `target`
/// does not appear linearly with a nonzero coefficient.
fn solve_linear_coeff<'a>(
    ctx: &'a AtomArena<'a>,
    cond: Atom<'a>,
    target: Atom<'a>,
) -> Option<Atom<'a>> {
    let terms: Vec<Atom<'a>> = match cond.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![cond],
    };

    let target_str = target.to_string();
    let mut coeff_terms: Vec<Atom<'a>> = Vec::new();
    let mut rest_terms: Vec<Atom<'a>> = Vec::new();

    for term in terms {
        // Does this term contain `target` as a linear factor?
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let target_count = factors
            .iter()
            .filter(|f| f.to_string() == target_str)
            .count();
        if target_count == 1 {
            let rest: Vec<_> = factors
                .iter()
                .filter(|f| f.to_string() != target_str)
                .copied()
                .collect();
            // The remaining factors must not contain the target nonlinearly.
            if rest.iter().any(|f| f.to_string().contains(&target_str)) {
                return None;
            }
            coeff_terms.push(if rest.is_empty() {
                ctx.num(1)
            } else {
                ctx.mul(&rest)
            });
        } else if target_count == 0 {
            if term.to_string().contains(&target_str) {
                return None; // target inside a sub-expression: not linear
            }
            rest_terms.push(term);
        } else {
            return None; // target^2 or higher: not linear
        }
    }

    if coeff_terms.is_empty() {
        return None;
    }

    let coeff = super::util::collect_terms(ctx, ctx.add(&coeff_terms));
    if matches!(coeff.node(), AtomNode::Num(0)) {
        return None;
    }
    let rest = if rest_terms.is_empty() {
        ctx.num(0)
    } else {
        super::util::collect_terms(ctx, ctx.add(&rest_terms))
    };

    let value = super::util::collect_terms(
        ctx,
        ctx.mul(&[ctx.num(-1), rest, ctx.pow(coeff, ctx.num(-1))]),
    );
    Some(value)
}

/// Solve a second-order linear ODE with regular singular point at `x0 = 0`
/// using the Frobenius method.
///
/// Sets $y = x^r \sum_{n=0}^{N-1} a_n x^n$, substitutes into the ODE, and
/// groups the residual by powers of $x$. The lowest-power group yields the
/// indicial equation $A r^2 + B r + C = 0$; subsequent groups determine
/// $a_n$ recursively. Only real rational roots of the indicial equation are
/// handled; the series for the larger root is returned.
pub(crate) fn solve_frobenius<'a>(
    ctx: &'a AtomArena<'a>,
    ode: super::ODE<'a>,
    x0: Atom<'a>,
    n_terms: usize,
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
    // Only the singular point x0 = 0 is supported.
    if !matches!(x0.node(), AtomNode::Num(0)) {
        return None;
    }

    let x = ctx.var(var.as_str());
    let (a, b, c, forcing) =
        super::second_order::extract_second_order_coeffs(ctx, equation, func, var)?;
    // Only homogeneous equations.
    if !matches!(
        super::util::collect_terms(ctx, forcing).node(),
        AtomNode::Num(0)
    ) {
        return None;
    }

    // Placeholder atoms: u = x^r (opaque factor), r = indicial root symbol.
    let u = ctx.var("XR");
    let r = ctx.var("r");

    // S = sum a_n x^n, with derivatives built manually so that the a_n
    // coefficient symbols are not differentiated by `diff` (which would
    // treat them as unknown functions of x).
    let a_syms: Vec<Atom<'a>> = (0..n_terms).map(|i| ctx.var(&format!("a{i}"))).collect();
    let (s, s1, s2) = build_series_triple(ctx, &a_syms, x);

    // y = u*S, y' = r*x^-1*u*S + u*S', y'' = r(r-1)*x^-2*u*S + 2r*x^-1*u*S' + u*S''.
    let y = ctx.mul(&[u, s]);
    let x_inv = ctx.pow(x, ctx.num(-1));
    let x_inv2 = ctx.pow(x, ctx.num(-2));
    let y1 = ctx.add(&[ctx.mul(&[r, x_inv, u, s]), ctx.mul(&[u, s1])]);
    let y2 = ctx.add(&[
        ctx.mul(&[r, ctx.add(&[r, ctx.num(-1)]), x_inv2, u, s]),
        ctx.mul(&[ctx.num(2), r, x_inv, u, s1]),
        ctx.mul(&[u, s2]),
    ]);

    // Residual R = a*y2 + b*y1 + c*y, with the common factor u removed.
    let residual = super::util::collect_terms(
        ctx,
        ctx.add(&[ctx.mul(&[a, y2]), ctx.mul(&[b, y1]), ctx.mul(&[c, y])]),
    );

    // Group residual terms by the integer power of x (after removing u).
    let terms: Vec<Atom<'a>> = match residual.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![residual],
    };
    let mut groups: Vec<(i64, Vec<Atom<'a>>)> = Vec::new();
    for term in terms {
        let (power, stripped) = strip_x_and_u(ctx, term, x, u)?;
        if let Some(g) = groups.iter_mut().find(|(p, _)| *p == power) {
            g.1.push(stripped);
        } else {
            groups.push((power, vec![stripped]));
        }
    }
    groups.sort_by_key(|(p, _)| *p);

    // The lowest group gives the indicial equation A*r^2 + B*r + C = 0
    // (its terms must be a_0 times a quadratic in r).
    let (_, lowest_terms) = groups.first()?.clone();
    let (ca, cb, cc) = indicial_coeffs(ctx, &lowest_terms, r, a_syms[0])?;
    // Solve A*r^2 + B*r + C = 0 for rational roots.
    let disc = cb * cb - 4 * ca * cc;
    if disc < 0 || ca == 0 {
        return None;
    }
    let sd = isqrt_i64(disc);
    if sd * sd != disc {
        return None; // irrational roots: not handled
    }
    // Larger root r1 = (-B + sd) / (2A).
    let r1_num = -cb + sd;
    let r1_den = 2 * ca;
    if r1_den == 0 {
        return None;
    }

    // r as an atom (reduced fraction).
    let g = gcd_i64(r1_num.unsigned_abs() as i64, r1_den.unsigned_abs() as i64).max(1);
    let (rn, rd) = (r1_num / g, r1_den / g);
    let r_val = if rd == 1 {
        ctx.num(rn)
    } else {
        ctx.mul(&[ctx.num(rn), ctx.pow(ctx.num(rd), ctx.num(-1))])
    };

    // Substitute r = r_val into the groups and solve for a_n recursively.
    let mut solved: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    for (idx, (_, g_terms)) in groups.iter().enumerate() {
        let mut cond = ctx.add(g_terms);
        cond = super::classify::replace_atom(ctx, cond, r, r_val);
        cond = super::util::collect_terms(ctx, cond);
        // Substitute already-solved coefficients.
        for (sym, val) in &solved {
            cond = super::classify::replace_atom(ctx, cond, *sym, *val);
        }
        cond = super::util::collect_terms(ctx, cond);

        if idx == 0 {
            // Indicial group: must vanish identically for the chosen root.
            if !matches!(cond.node(), AtomNode::Num(0)) {
                return None;
            }
            continue;
        }

        // Solve for the highest-index a_n not yet solved, excluding the
        // free parameter a_0. Groups whose coefficients are all determined
        // are truncation artifacts (the dropped tail a_N, a_{N+1}, ...
        // would contribute here) and are simply skipped.
        let target = a_syms.iter().skip(1).rev().find(|s| {
            cond.to_string().contains(&s.to_string())
                && !solved
                    .iter()
                    .any(|(sym, _)| sym.to_string() == s.to_string())
        });
        let Some(&target) = target else {
            continue;
        };
        let value = solve_linear_coeff(ctx, cond, target)?;
        solved.push((target, value));
    }

    if solved.is_empty() {
        return None;
    }

    // Rebuild y = x^r * S with solved coefficients.
    let mut series = s;
    for (sym, val) in &solved {
        series = super::classify::replace_atom(ctx, series, *sym, *val);
    }
    series = super::util::collect_terms(ctx, series);
    let y_sol = super::util::collect_terms(ctx, ctx.mul(&[ctx.pow(x, r_val), series]));

    Some(ODESolution::Series(y_sol, n_terms))
}

/// Build (S, S', S'') for S = sum_{n} a_n x^n, treating a_n as constants.
fn build_series_triple<'a>(
    ctx: &'a AtomArena<'a>,
    a_syms: &[Atom<'a>],
    x: Atom<'a>,
) -> (Atom<'a>, Atom<'a>, Atom<'a>) {
    let mut s_terms = Vec::with_capacity(a_syms.len());
    let mut s1_terms = Vec::new();
    let mut s2_terms = Vec::new();
    for (n, an) in a_syms.iter().enumerate() {
        if n == 0 {
            s_terms.push(*an);
        } else {
            s_terms.push(ctx.mul(&[*an, ctx.pow(x, ctx.num(n as i64))]));
            // S' term: n * a_n * x^(n-1)
            let d1 = if n == 1 {
                ctx.mul(&[ctx.num(n as i64), *an])
            } else {
                ctx.mul(&[ctx.num(n as i64), *an, ctx.pow(x, ctx.num(n as i64 - 1))])
            };
            s1_terms.push(d1);
            // S'' term: n*(n-1) * a_n * x^(n-2)
            if n >= 2 {
                let d2 = if n == 2 {
                    ctx.mul(&[ctx.num((n * (n - 1)) as i64), *an])
                } else {
                    ctx.mul(&[
                        ctx.num((n * (n - 1)) as i64),
                        *an,
                        ctx.pow(x, ctx.num(n as i64 - 2)),
                    ])
                };
                s2_terms.push(d2);
            }
        }
    }
    let s = ctx.add(&s_terms);
    let s1 = if s1_terms.is_empty() {
        ctx.num(0)
    } else {
        ctx.add(&s1_terms)
    };
    let s2 = if s2_terms.is_empty() {
        ctx.num(0)
    } else {
        ctx.add(&s2_terms)
    };
    (s, s1, s2)
}

/// Remove the opaque `u` factor and the `x`-power from a residual term,
/// returning (power_of_x, remaining_expression).
fn strip_x_and_u<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    x: Atom<'a>,
    u: Atom<'a>,
) -> Option<(i64, Atom<'a>)> {
    let factors: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut power: i64 = 0;
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        if f.to_string() == u.to_string() {
            continue; // drop u
        }
        match f.node() {
            AtomNode::Var(v) if *v == Symbol::new("x") && f.to_string() == x.to_string() => {
                power += 1;
            }
            AtomNode::Pow(base, exp) => {
                if base.to_string() == x.to_string() {
                    if let AtomNode::Num(n) = exp.node() {
                        power += n;
                    } else {
                        return None;
                    }
                } else {
                    rest.push(f);
                }
            }
            _ => rest.push(f),
        }
    }
    let stripped = if rest.is_empty() {
        ctx.num(1)
    } else {
        ctx.mul(&rest)
    };
    Some((power, stripped))
}

/// Extract (A, B, C) of the indicial equation A*r^2 + B*r + C = 0 from the
/// lowest-power group terms. The terms must be `a_0` times a quadratic in r
/// with integer coefficients.
fn indicial_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    terms: &[Atom<'a>],
    r: Atom<'a>,
    a0: Atom<'a>,
) -> Option<(i64, i64, i64)> {
    let r_str = r.to_string();
    let a0_str = a0.to_string();
    let mut ca = 0i64;
    let mut cb = 0i64;
    let mut cc = 0i64;
    for term in terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![*term],
        };
        let mut r_pow = 0i64;
        let mut num: i64 = 1;
        for f in factors {
            match f.node() {
                AtomNode::Num(n) => num *= n,
                AtomNode::Var(_) if f.to_string() == r_str => r_pow += 1,
                AtomNode::Pow(base, exp) => {
                    if base.to_string() == r_str {
                        if let AtomNode::Num(n) = exp.node() {
                            r_pow += n;
                        } else {
                            return None;
                        }
                    } else {
                        return None;
                    }
                }
                AtomNode::Add(args) => {
                    // (r + -1) style: expand manually for r(r-1).
                    // Support Add of [r, Num(-1)] only.
                    let mut has_r = false;
                    let mut has_m1 = false;
                    for aa in args.iter() {
                        if aa.to_string() == r_str {
                            has_r = true;
                        } else if matches!(aa.node(), AtomNode::Num(-1)) {
                            has_m1 = true;
                        }
                    }
                    if has_r && has_m1 && args.len() == 2 {
                        // This factor is (r - 1): combined with an existing
                        // r factor it forms r(r-1) = r^2 - r. Handle by
                        // returning a marker: we expand below.
                        return expand_indicial(ctx, terms, r, a0);
                    }
                    return None;
                }
                _ => {
                    if f.to_string() != a0_str {
                        return None;
                    }
                }
            }
        }
        match r_pow {
            2 => ca += num,
            1 => cb += num,
            0 => cc += num,
            _ => return None,
        }
    }
    Some((ca, cb, cc))
}

/// Expand products in the indicial group (handles (r-1) factors) and retry
/// coefficient extraction.
fn expand_indicial<'a>(
    ctx: &'a AtomArena<'a>,
    terms: &[Atom<'a>],
    r: Atom<'a>,
    a0: Atom<'a>,
) -> Option<(i64, i64, i64)> {
    // Multiply out each term containing (r-1) factors, then re-extract.
    let mut expanded: Vec<Atom<'a>> = Vec::new();
    for term in terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![*term],
        };
        let mut current: Vec<Atom<'a>> = vec![ctx.num(1)];
        for f in factors {
            if let AtomNode::Add(args) = f.node() {
                let mut next: Vec<Atom<'a>> = Vec::new();
                for c in &current {
                    for aa in args.iter() {
                        next.push(ctx.mul(&[*c, *aa]));
                    }
                }
                current = next;
            } else {
                for c in &mut current {
                    *c = ctx.mul(&[*c, f]);
                }
            }
        }
        expanded.extend(current);
    }
    let expanded_sum = super::util::collect_terms(ctx, ctx.add(&expanded));
    let new_terms: Vec<Atom<'a>> = match expanded_sum.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expanded_sum],
    };
    // Retry extraction on the expanded terms (no Add factors remain).
    let r_str = r.to_string();
    let a0_str = a0.to_string();
    let mut ca = 0i64;
    let mut cb = 0i64;
    let mut cc = 0i64;
    for term in new_terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let mut r_pow = 0i64;
        let mut num: i64 = 1;
        for f in factors {
            match f.node() {
                AtomNode::Num(n) => num *= n,
                AtomNode::Var(_) if f.to_string() == r_str => r_pow += 1,
                _ => {
                    if f.to_string() != a0_str {
                        return None;
                    }
                }
            }
        }
        match r_pow {
            2 => ca += num,
            1 => cb += num,
            0 => cc += num,
            _ => return None,
        }
    }
    Some((ca, cb, cc))
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

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

/// Substitute series expressions for y, y', y'' into the ODE equation.
fn substitute_series<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
    y_val: Atom<'a>,
    y1_val: Atom<'a>,
    y2_val: Atom<'a>,
) -> Atom<'a> {
    substitute_series_inner(ctx, equation, func, var, y_val, y1_val, y2_val)
}

fn substitute_series_inner<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
    y_val: Atom<'a>,
    y1_val: Atom<'a>,
    y2_val: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_series_inner(ctx, *a, func, var, y_val, y1_val, y2_val))
                .collect();
            ctx.add(&mapped)
        }
        AtomNode::Mul(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_series_inner(ctx, *a, func, var, y_val, y1_val, y2_val))
                .collect();
            ctx.mul(&mapped)
        }
        AtomNode::Pow(base, exp) => {
            let b = substitute_series_inner(ctx, *base, func, var, y_val, y1_val, y2_val);
            let e = substitute_series_inner(ctx, *exp, func, var, y_val, y1_val, y2_val);
            ctx.pow(b, e)
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                let is_func =
                    args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str();
                if is_func {
                    let deriv_order = args.len() - 1;
                    return match deriv_order {
                        1 => y1_val,
                        2 => y2_val,
                        _ => {
                            // Higher derivatives: compute from series.
                            let mut result = y_val;
                            for _ in 0..deriv_order {
                                result = diff(ctx, result, var);
                            }
                            result
                        }
                    };
                }
            }
            // Check if this is func itself (not a derivative).
            if expr.to_string() == func.to_string() {
                return y_val;
            }
            // Otherwise, recurse into arguments.
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_series_inner(ctx, *a, func, var, y_val, y1_val, y2_val))
                .collect();
            ctx.fun(name.as_str(), &mapped)
        }
    }
}
