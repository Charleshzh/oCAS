//! Linear ODE system solvers.
//!
//! Provides [`solve_linear_system`] for constant-coefficient systems of the
//! form $\mathbf{Y}' = A\mathbf{Y} + \mathbf{g}(x)$.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::ODESolution;
use super::util::collect_terms;

/// Solve a constant-coefficient linear ODE system.
///
/// Given a system $\mathbf{Y}' = A\mathbf{Y}$ where $A$ is a constant
/// matrix, attempts to find the general solution. Currently supports
/// 2×2 systems whose characteristic polynomial has closed-form roots
/// (integer or rational eigenvalues, or complex conjugate pairs).
///
/// Each equation must be in the form `Derivative(y_i, x) - (A_i·Y) = 0`.
///
/// Returns `Some(ODESolution::System)` with the component solutions
/// containing free constants C1, C2, or `None` when the system is not
/// supported.
pub(crate) fn solve_linear_system<'a>(
    ctx: &'a AtomArena<'a>,
    equations: &[Atom<'a>],
    funcs: &[Atom<'a>],
    var: Symbol,
) -> Option<ODESolution<'a>> {
    if equations.len() != 2 || funcs.len() != 2 {
        return None;
    }
    let x = ctx.var(var.as_str());

    // Extract the 2×2 coefficient matrix A such that y_i' = A_i1*y1 + A_i2*y2.
    let a11 = extract_coeff(ctx, equations[0], funcs[0], var)?;
    let a12 = extract_coeff(ctx, equations[0], funcs[1], var)?;
    let a21 = extract_coeff(ctx, equations[1], funcs[0], var)?;
    let a22 = extract_coeff(ctx, equations[1], funcs[1], var)?;

    // Characteristic polynomial: λ^2 - (tr)λ + det = 0.
    let trace = a11 + a22;
    let det = a11 * a22 - a12 * a21;
    let disc = trace * trace - 4 * det;
    if disc < 0 {
        // Complex eigenvalues: α ± βi with α = trace/2, β = sqrt(-disc)/2.
        // Only even trace and perfect-square -disc are supported.
        if trace % 2 != 0 {
            return None;
        }
        let beta_sq = -disc;
        let sb = isqrt_i64(beta_sq);
        if sb * sb != beta_sq || sb % 2 != 0 {
            return None;
        }
        let alpha = trace / 2;
        let beta = sb / 2;
        return Some(solve_complex_2x2(ctx, x, a11, a12, a21, alpha, beta));
    }

    let sd = isqrt_i64(disc);
    if sd * sd != disc || (trace + sd) % 2 != 0 {
        return None; // irrational eigenvalues not supported
    }
    let l1 = (trace + sd) / 2;
    let l2 = (trace - sd) / 2;

    let c1 = ctx.var("C1");
    let c2 = ctx.var("C2");
    let e1 = ctx.fun("exp", &[ctx.mul(&[ctx.num(l1), x])]);

    if l1 == l2 {
        // Repeated eigenvalue: y = e^{λx} (C1 v + C2 (t v + w)).
        // v = eigenvector, w = generalized eigenvector solving (A-λI)w = v.
        // Only the case A ≠ λI (defective) or A = λI (complete) is handled.
        let e_lx = e1;
        if a11 == l1 && a12 == 0 && a21 == 0 && a22 == l1 {
            // A = λI: two independent eigenvectors e1, e2.
            let y1 = ctx.mul(&[c1, e_lx]);
            let y2 = ctx.mul(&[c2, e_lx]);
            return Some(ODESolution::System(ctx.slice(&[y1, y2])));
        }
        // Defective: v from (A-λI)v = 0; w from (A-λI)w = v.
        let (v1, v2) = eigenvector(a11, a12, a21, a22, l1)?;
        let (w1, w2) = generalized_eigenvector(a11, a12, a21, a22, l1, v1, v2)?;
        // y1 = e^{λx}(C1 v1 + C2 (x v1 + w1)), y2 = e^{λx}(C1 v2 + C2 (x v2 + w2))
        let xv1_w1 = collect_terms(ctx, ctx.add(&[ctx.mul(&[x, ctx.num(v1)]), ctx.num(w1)]));
        let xv2_w2 = collect_terms(ctx, ctx.add(&[ctx.mul(&[x, ctx.num(v2)]), ctx.num(w2)]));
        let y1 = collect_terms(
            ctx,
            ctx.mul(&[
                e_lx,
                ctx.add(&[ctx.mul(&[c1, ctx.num(v1)]), ctx.mul(&[c2, xv1_w1])]),
            ]),
        );
        let y2 = collect_terms(
            ctx,
            ctx.mul(&[
                e_lx,
                ctx.add(&[ctx.mul(&[c1, ctx.num(v2)]), ctx.mul(&[c2, xv2_w2])]),
            ]),
        );
        return Some(ODESolution::System(ctx.slice(&[y1, y2])));
    }

    // Distinct real eigenvalues: y = C1 v1 e^{λ1 x} + C2 v2 e^{λ2 x}.
    let e2 = ctx.fun("exp", &[ctx.mul(&[ctx.num(l2), x])]);
    let (v11, v12) = eigenvector(a11, a12, a21, a22, l1)?;
    let (v21, v22) = eigenvector(a11, a12, a21, a22, l2)?;
    let y1 = collect_terms(
        ctx,
        ctx.add(&[
            ctx.mul(&[c1, ctx.num(v11), e1]),
            ctx.mul(&[c2, ctx.num(v21), e2]),
        ]),
    );
    let y2 = collect_terms(
        ctx,
        ctx.add(&[
            ctx.mul(&[c1, ctx.num(v12), e1]),
            ctx.mul(&[c2, ctx.num(v22), e2]),
        ]),
    );
    Some(ODESolution::System(ctx.slice(&[y1, y2])))
}

/// Complex eigenvalues α ± βi: real fundamental solutions.
///
/// With (v + iw) a complex eigenvector for α + βi, real solutions are
/// e^{αx}(v cos βx - w sin βx) and e^{αx}(v sin βx + w cos βx).
/// For real 2×2 A with eigenvector components we use the standard real form:
/// y = e^{αx}[ C1 (p cos βx - q sin βx) + C2 (p sin βx + q cos βx) ]
/// where p, q come from the real/imaginary parts of the eigenvector.
fn solve_complex_2x2<'a>(
    ctx: &'a AtomArena<'a>,
    x: Atom<'a>,
    a11: i64,
    a12: i64,
    a21: i64,
    alpha: i64,
    beta: i64,
) -> ODESolution<'a> {
    // For the eigenvalue α + βi, the eigenvector equation (A - λI)v = 0
    // gives (a11 - α - iβ) v1 + a12 v2 = 0. Choose v1 = a12, then
    // v2 = -(a11 - α - iβ) = (α - a11) + iβ.
    // Real part p = (a12, α - a11), imaginary part q = (0, β).
    let p1 = a12;
    let p2 = alpha - a11;
    let q1 = 0i64;
    let q2 = beta;
    let _ = a21; // symmetric form uses only the first row

    let c1 = ctx.var("C1");
    let c2 = ctx.var("C2");
    let bx = ctx.mul(&[ctx.num(beta), x]);
    let cos_bx = ctx.fun("cos", &[bx]);
    let sin_bx = ctx.fun("sin", &[bx]);
    let e_ax = if alpha == 0 {
        ctx.num(1)
    } else {
        ctx.fun("exp", &[ctx.mul(&[ctx.num(alpha), x])])
    };

    // y1 = e^{αx}[ C1 (p1 cos - q1 sin) + C2 (p1 sin + q1 cos) ]
    let y1_inner = collect_terms(
        ctx,
        ctx.add(&[
            ctx.mul(&[
                c1,
                ctx.add(&[
                    ctx.mul(&[ctx.num(p1), cos_bx]),
                    ctx.mul(&[ctx.num(-q1), sin_bx]),
                ]),
            ]),
            ctx.mul(&[
                c2,
                ctx.add(&[
                    ctx.mul(&[ctx.num(p1), sin_bx]),
                    ctx.mul(&[ctx.num(q1), cos_bx]),
                ]),
            ]),
        ]),
    );
    // y2 = e^{αx}[ C1 (p2 cos - q2 sin) + C2 (p2 sin + q2 cos) ]
    let y2_inner = collect_terms(
        ctx,
        ctx.add(&[
            ctx.mul(&[
                c1,
                ctx.add(&[
                    ctx.mul(&[ctx.num(p2), cos_bx]),
                    ctx.mul(&[ctx.num(-q2), sin_bx]),
                ]),
            ]),
            ctx.mul(&[
                c2,
                ctx.add(&[
                    ctx.mul(&[ctx.num(p2), sin_bx]),
                    ctx.mul(&[ctx.num(q2), cos_bx]),
                ]),
            ]),
        ]),
    );

    let y1 = if alpha == 0 {
        y1_inner
    } else {
        collect_terms(ctx, ctx.mul(&[e_ax, y1_inner]))
    };
    let y2 = if alpha == 0 {
        y2_inner
    } else {
        collect_terms(ctx, ctx.mul(&[e_ax, y2_inner]))
    };
    ODESolution::System(ctx.slice(&[y1, y2]))
}

/// Extract the numeric coefficient of `func` in `eq` (which contains
/// `Derivative(func_i, x)` and linear `func` terms).
fn extract_coeff<'a>(
    ctx: &'a AtomArena<'a>,
    eq: Atom<'a>,
    func: Atom<'a>,
    var: Symbol,
) -> Option<i64> {
    let func_str = func.to_string();
    let dy_str = ctx
        .fun("Derivative", &[func, ctx.var(var.as_str())])
        .to_string();
    let _ = dy_str;

    let terms: Vec<Atom<'a>> = match eq.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![eq],
    };
    let mut total: i64 = 0;
    for term in terms {
        let factors: Vec<Atom<'a>> = match term.node() {
            AtomNode::Mul(args) => args.to_vec(),
            _ => vec![term],
        };
        let has_func = factors.iter().any(|f| f.to_string() == func_str);
        if has_func {
            let coeff: i64 = factors
                .iter()
                .filter(|f| f.to_string() != func_str)
                .map(|f| match f.node() {
                    AtomNode::Num(n) => *n,
                    _ => 0,
                })
                .product::<i64>()
                .max(i64::MIN + 1); // empty product = 1
            let coeff = if factors.len() == 1 { 1 } else { coeff };
            // The equation is y' - (A Y) = 0, so the coefficient of y in
            // the equation is -A_ij: negate to recover A_ij.
            total += -coeff;
        }
    }
    Some(total)
}

/// Find an eigenvector (v1, v2) with integer entries for eigenvalue λ.
fn eigenvector(a11: i64, a12: i64, a21: i64, a22: i64, lambda: i64) -> Option<(i64, i64)> {
    // (A - λI)v = 0: (a11-λ) v1 + a12 v2 = 0, a21 v1 + (a22-λ) v2 = 0.
    let m11 = a11 - lambda;
    let m22 = a22 - lambda;
    if a12 != 0 {
        // v = (a12, -m11) = (a12, λ - a11)
        return Some((a12, -m11));
    }
    if a21 != 0 {
        // v = (-m22, a21) = (λ - a22, a21)
        return Some((-m22, a21));
    }
    // Diagonal case: distinct handled by caller; here both rows zero.
    None
}

/// Find a generalized eigenvector w with (A - λI) w = v.
fn generalized_eigenvector(
    a11: i64,
    a12: i64,
    a21: i64,
    a22: i64,
    lambda: i64,
    v1: i64,
    v2: i64,
) -> Option<(i64, i64)> {
    // Solve (a11-λ) w1 + a12 w2 = v1, a21 w1 + (a22-λ) w2 = v2.
    // One equation is a multiple of the other; use the first nonzero row.
    let m11 = a11 - lambda;
    let m22 = a22 - lambda;
    if a12 != 0 {
        // w2 = (v1 - m11 w1)/a12 with w1 = 0.
        if v1 % a12 != 0 {
            // Try w1 = 1 if direct fails.
            let w1 = 1;
            let num = v1 - m11 * w1;
            if num % a12 == 0 {
                return Some((w1, num / a12));
            }
            return None;
        }
        return Some((0, v1 / a12));
    }
    if a21 != 0 {
        if v2 % a21 == 0 {
            return Some((v2 / a21, 0));
        }
        let w2 = 1;
        let num = v2 - m22 * w2;
        if num % a21 == 0 {
            return Some((num / a21, w2));
        }
        return None;
    }
    None
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
