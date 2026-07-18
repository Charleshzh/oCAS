//! Special-function antiderivatives (the Meijer-G fallback's endpoint).
//!
//! When the Risch algorithm proves an integral has no elementary
//! antiderivative, many common cases still have closed forms in terms of
//! special functions. Rather than routing through the full Meijer
//! G-function machinery (which requires hypergeometric series, gamma
//! functions, and Slater expansions oCAS does not yet provide), the
//! antiderivatives of the standard non-elementary integrals are encoded
//! directly:
//!
//! - `exp(±x²)` → `erf` / `erfi`
//! - `exp(x)/x` → `Ei`
//! - `sin(x)/x` → `Si`,  `cos(x)/x` → `Ci`,  `cosh(x)/x` → `Chi`, `sinh(x)/x` → `Shi`
//! - `sin(x²)`, `cos(x²)` → Fresnel `S` / `C`
//!
//! The entries below follow the same definitions as SymPy, so results
//! compare equal against `sympy.integrate`.

use ocas_atom::{Atom, AtomArena, AtomNode};

/// Try to integrate `expr` in terms of special functions.
///
/// Only the patterns listed in the module docs are recognized; returns
/// `None` otherwise (the caller then emits the unevaluated form).
pub(crate) fn special_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    // All table entries need the integrand normalized as a product; we
    // inspect it as a flat list of factors.
    let factors: Vec<Atom> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    erf_family(ctx, &factors, x)
        .or_else(|| ei_family(ctx, &factors, x))
        .or_else(|| trig_integral_family(ctx, &factors, x))
        .or_else(|| fresnel_family(ctx, &factors, x))
}

// ------------------------------------------------------------------
//  Pattern helpers
// ------------------------------------------------------------------

/// Match `exp(u)`; return `u`.
fn as_exp(f: Atom) -> Option<Atom> {
    if let AtomNode::Fun(name, args) = f.node() {
        if name.as_str() == "exp" && args.len() == 1 {
            return Some(args[0]);
        }
    }
    None
}

/// Match `x^-1`.
fn is_x_inv<'a>(f: Atom<'a>, x: Atom<'a>) -> bool {
    matches!(f.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(-1)))
}

/// Match `c·x^2` (with optional sign): returns the coefficient atom `c`.
fn as_quadratic<'a>(u: Atom<'a>, x: Atom<'a>) -> Option<Atom<'a>> {
    // x^2
    if matches!(u.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(2))) {
        return None; // coefficient 1 handled by caller
    }
    // c·x^2 or -x^2
    if let AtomNode::Mul(factors) = u.node() {
        if factors.len() == 2 {
            if matches!(factors[1].node(), AtomNode::Pow(b, e)
                if *b == x && matches!(e.node(), AtomNode::Num(2)))
            {
                return Some(factors[0]);
            }
        }
    }
    None
}

/// Whether `u` is exactly `x²`, or `(-x)²` which normalizes to it.
fn is_x_squared<'a>(u: Atom<'a>, x: Atom<'a>) -> bool {
    if matches!(u.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(2))) {
        return true;
    }
    // (-x)^2 or (-1·x)^2 → x²
    if let AtomNode::Pow(b, e) = u.node() {
        if matches!(e.node(), AtomNode::Num(2)) {
            if let AtomNode::Mul(factors) = b.node() {
                let all_neg_one_or_x = factors
                    .iter()
                    .all(|f| matches!(f.node(), AtomNode::Num(-1)) || *f == x);
                let has_x = factors.contains(&x);
                return all_neg_one_or_x && has_x;
            }
        }
    }
    false
}

/// Whether `u` is exactly `x`.
fn is_x<'a>(u: Atom<'a>, x: Atom<'a>) -> bool {
    u == x
}

// ------------------------------------------------------------------
//  erf family: exp(-x²), exp(c·x²)
// ------------------------------------------------------------------

fn erf_family<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], x: Atom<'a>) -> Option<Atom<'a>> {
    if factors.len() != 1 {
        return None;
    }
    let u = as_exp(factors[0])?;
    // exp(-x²) → (√π/2)·erf(x)
    if is_x_squared(u, x) {
        // exp(+x²) → (√π/2)·erfi(x)
        let sqrt_pi = ctx.fun("sqrt", &[ctx.var("pi")]);
        let erfi = ctx.fun("erfi", &[x]);
        return Some(ctx.mul(&[sqrt_pi, ctx.pow(ctx.num(2), ctx.num(-1)), erfi]));
    }
    // exp(c·x²) with negative c: √π/(2√(-c))·erf(√(-c)·x)
    if let Some(c) = as_quadratic(u, x) {
        let neg_c = ctx.mul(&[ctx.num(-1), c]);
        let root = ctx.fun("sqrt", &[neg_c]);
        let sqrt_pi = ctx.fun("sqrt", &[ctx.var("pi")]);
        let erf = ctx.fun("erf", &[ctx.mul(&[root, x])]);
        return Some(ctx.mul(&[
            sqrt_pi,
            ctx.pow(ctx.mul(&[ctx.num(2), root]), ctx.num(-1)),
            erf,
        ]));
    }
    None
}

// ------------------------------------------------------------------
//  Ei family: exp(x)/x, exp(c·x)/x
// ------------------------------------------------------------------

fn ei_family<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], x: Atom<'a>) -> Option<Atom<'a>> {
    if factors.len() != 2 {
        return None;
    }
    let (exp_f, inv_f) = if as_exp(factors[0]).is_some() {
        (factors[0], factors[1])
    } else if as_exp(factors[1]).is_some() {
        (factors[1], factors[0])
    } else {
        return None;
    };
    if !is_x_inv(inv_f, x) {
        return None;
    }
    let u = as_exp(exp_f)?;
    // exp(x)/x → Ei(x)
    if is_x(u, x) {
        return Some(ctx.fun("Ei", &[x]));
    }
    // exp(c·x)/x → Ei(c·x) for constant c.
    if let AtomNode::Mul(mf) = u.node() {
        if mf.len() == 2 && mf[1] == x && matches!(mf[0].node(), AtomNode::Num(_)) {
            return Some(ctx.fun("Ei", &[u]));
        }
    }
    None
}

// ------------------------------------------------------------------
//  Si/Ci/Shi/Chi family: sin(x)/x, cos(x)/x, sinh(x)/x, cosh(x)/x
// ------------------------------------------------------------------

fn trig_integral_family<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    if factors.len() != 2 {
        return None;
    }
    let (fun_f, inv_f) = if is_x_inv(factors[1], x) {
        (factors[0], factors[1])
    } else if is_x_inv(factors[0], x) {
        (factors[1], factors[0])
    } else {
        return None;
    };
    let _ = inv_f;
    let AtomNode::Fun(name, args) = fun_f.node() else {
        return None;
    };
    if args.len() != 1 || !is_x(args[0], x) {
        return None;
    }
    let target = match name.as_str() {
        "sin" => "Si",
        "cos" => "Ci",
        "sinh" => "Shi",
        "cosh" => "Chi",
        _ => return None,
    };
    Some(ctx.fun(target, &[x]))
}

// ------------------------------------------------------------------
//  Fresnel family: sin(x²), cos(x²)
// ------------------------------------------------------------------

fn fresnel_family<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    if factors.len() != 1 {
        return None;
    }
    let AtomNode::Fun(name, args) = factors[0].node() else {
        return None;
    };
    if args.len() != 1 || !is_x_squared(args[0], x) {
        return None;
    }
    // ∫ sin(x²) dx = √(π/2)·S(√(2/π)·x); same prefactor for cos → C.
    let target = match name.as_str() {
        "sin" => "fresnels",
        "cos" => "fresnelc",
        _ => return None,
    };
    let pi = ctx.var("pi");
    let two = ctx.num(2);
    // √(π/2) = sqrt(pi·2⁻¹)
    let prefactor = ctx.fun("sqrt", &[ctx.mul(&[pi, ctx.pow(two, ctx.num(-1))])]);
    // √(2/π)·x = sqrt(2·pi⁻¹)·x
    let inner = ctx.mul(&[
        ctx.fun("sqrt", &[ctx.mul(&[two, ctx.pow(pi, ctx.num(-1))])]),
        x,
    ]);
    Some(ctx.mul(&[prefactor, ctx.fun(target, &[inner])]))
}

#[cfg(test)]
mod tests {
    use ocas_atom::{AtomArena, Symbol};
    use ocas_core::arena::Arena;

    use super::*;

    fn int<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Option<Atom<'a>> {
        special_integrate(ctx, expr, ctx.var("x"))
    }

    #[test]
    fn exp_neg_x_squared_gives_erf() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let r = int(&ctx, expr).expect("erf form");
        assert!(r.to_string().contains("erf"), "got {r}");
    }

    #[test]
    fn exp_x_over_x_gives_ei() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("exp", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Ei form");
        assert_eq!(r.to_string(), "Ei(x)");
    }

    #[test]
    fn sin_x_over_x_gives_si() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("sin", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Si form");
        assert_eq!(r.to_string(), "Si(x)");
    }

    #[test]
    fn cos_x_over_x_gives_ci() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("cos", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Ci form");
        assert_eq!(r.to_string(), "Ci(x)");
    }

    #[test]
    fn sin_x_squared_gives_fresnel() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[ctx.pow(x, ctx.num(2))]);
        let r = int(&ctx, expr).expect("Fresnel form");
        assert!(r.to_string().contains("fresnels"), "got {r}");
    }

    #[test]
    fn unmatched_returns_none() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x) itself is elementary (handled by Risch).
        assert!(int(&ctx, ctx.fun("exp", &[x])).is_none());
        let _ = Symbol::new("x");
    }
}
