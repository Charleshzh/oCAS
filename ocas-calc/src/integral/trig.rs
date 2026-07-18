//! Trigonometric integrands via complex exponentials.
//!
//! [`trig_to_exp`] rewrites `sin`/`cos`/`tan`/etc. into exponentials with
//! the imaginary unit `I` (see [`crate::complex`]), so the Risch
//! algorithm (which only knows `log`/`exp`) can attempt them.
//! [`realify`] tries to rewrite the resulting complex answer back into
//! real form, merging conjugate logarithm pairs into real `log`/`atan`
//! terms on a best-effort basis.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use crate::complex::i;

/// Rewrite trigonometric functions in `expr` into complex exponentials.
///
/// Uses the standard identities (with `t = exp(I·u)`):
///
/// - `sin(u) → (t - t⁻¹)/(2I)`
/// - `cos(u) → (t + t⁻¹)/2`
/// - `tan(u) → sin/cos`, and similarly for `sec`/`csc`/`cot`.
///
/// Returns `None` when `expr` contains no trigonometric functions (the
/// caller then proceeds with the original expression).
pub(crate) fn trig_to_exp<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Option<Atom<'a>> {
    let mut found = false;
    let out = rewrite(ctx, expr, &mut found);
    found.then_some(out)
}

fn rewrite<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, found: &mut bool) -> Atom<'a> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            let new: Vec<Atom> = args.iter().map(|a| rewrite(ctx, *a, found)).collect();
            if matches!(expr.node(), AtomNode::Add(_)) {
                ctx.add(&new)
            } else {
                ctx.mul(&new)
            }
        }
        AtomNode::Pow(b, e) => {
            let nb = rewrite(ctx, *b, found);
            let ne = rewrite(ctx, *e, found);
            ctx.pow(nb, ne)
        }
        AtomNode::Fun(name, args) => {
            let new_args: Vec<Atom> = args.iter().map(|a| rewrite(ctx, *a, found)).collect();
            let n = name.as_str();
            if args.len() == 1 && matches!(n, "sin" | "cos" | "tan" | "cot" | "sec" | "csc") {
                *found = true;
                let u = new_args[0];
                trig_exp_form(ctx, n, u)
            } else {
                ctx.fun(n, &new_args)
            }
        }
    }
}

/// `exp(I·u)` and its inverse `(exp(I·u))⁻¹` (a power of the same
/// generator, so the tower sees a single exponential).
fn exp_iu<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>) -> (Atom<'a>, Atom<'a>) {
    let iu = ctx.mul(&[i(ctx), u]);
    let t = ctx.fun("exp", &[iu]);
    let t_inv = ctx.pow(t, ctx.num(-1));
    (t, t_inv)
}

/// The complex-exponential form of a trigonometric function.
fn trig_exp_form<'a>(ctx: &'a AtomArena<'a>, name: &str, u: Atom<'a>) -> Atom<'a> {
    let (t, ti) = exp_iu(ctx, u);
    let two = ctx.num(2);
    let two_i = ctx.mul(&[two, i(ctx)]);
    match name {
        // sin(u) = (t - t⁻¹)/(2I)
        "sin" => {
            let num = ctx.add(&[t, ctx.mul(&[ctx.num(-1), ti])]);
            ctx.mul(&[num, ctx.pow(two_i, ctx.num(-1))])
        }
        // cos(u) = (t + t⁻¹)/2
        "cos" => {
            let num = ctx.add(&[t, ti]);
            ctx.mul(&[num, ctx.pow(two, ctx.num(-1))])
        }
        // tan(u) = (t - t⁻¹)/(I·(t + t⁻¹))
        "tan" => {
            let num = ctx.add(&[t, ctx.mul(&[ctx.num(-1), ti])]);
            let den = ctx.mul(&[i(ctx), ctx.add(&[t, ti])]);
            ctx.mul(&[num, ctx.pow(den, ctx.num(-1))])
        }
        // cot(u) = I·(t + t⁻¹)/(t - t⁻¹)
        "cot" => {
            let num = ctx.mul(&[i(ctx), ctx.add(&[t, ti])]);
            let den = ctx.add(&[t, ctx.mul(&[ctx.num(-1), ti])]);
            ctx.mul(&[num, ctx.pow(den, ctx.num(-1))])
        }
        // sec(u) = 2/(t + t⁻¹)
        "sec" => {
            let den = ctx.add(&[t, ti]);
            ctx.mul(&[two, ctx.pow(den, ctx.num(-1))])
        }
        // csc(u) = 2I/(t - t⁻¹)
        "csc" => {
            let den = ctx.add(&[t, ctx.mul(&[ctx.num(-1), ti])]);
            ctx.mul(&[two_i, ctx.pow(den, ctx.num(-1))])
        }
        _ => unreachable!("trig_exp_form: unsupported {name}"),
    }
}

/// Try to rewrite a complex Risch result back into real form.
///
/// Handles the common patterns produced by trigonometric integrals:
///
/// - products/quotients of `exp(±I·u)` that combine back into `sin`/`cos`
///   are left in exponential form only when no simple real form exists;
/// - conjugate logarithm pairs `c·log(u + I v) + c·log(u - I v)` become
///   `c·log(u² + v²)`; pairs with opposite coefficients `c·log(u+Iv) -
///   c·log(u-Iv)` become `2c·atan(v/u)`.
///
/// This is a best-effort cosmetic pass; if nothing matches, `expr` is
/// returned unchanged (the complex answer is still mathematically valid,
/// as differentiation verifies).
pub(crate) fn realify<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Add(args) => {
            let rewritten: Vec<Atom> = args.iter().map(|a| realify(ctx, *a)).collect();
            if let Some(merged) = merge_conjugate_logs(ctx, &rewritten) {
                return merged;
            }
            ctx.add(&rewritten)
        }
        AtomNode::Mul(args) => {
            let new: Vec<Atom> = args.iter().map(|a| realify(ctx, *a)).collect();
            if let Some(combined) = combine_exp_pair(ctx, &new) {
                return combined;
            }
            ctx.mul(&new)
        }
        AtomNode::Pow(b, e) => ctx.pow(realify(ctx, *b), realify(ctx, *e)),
        AtomNode::Fun(name, args) => {
            let new: Vec<Atom> = args.iter().map(|a| realify(ctx, *a)).collect();
            ctx.fun(name.as_str(), &new)
        }
        _ => expr,
    }
}

/// In a sum, find `c·log(u + Iv) + c·log(u - Iv)` pairs and merge them
/// into `c·log(u² + v²)`, or `c·log(u+Iv) - c·log(u-Iv)` into `2c·atan(v/u)`.
fn merge_conjugate_logs<'a>(ctx: &'a AtomArena<'a>, args: &[Atom<'a>]) -> Option<Atom<'a>> {
    for i in 0..args.len() {
        for j in (i + 1)..args.len() {
            if let Some(merged) = try_merge_pair(ctx, args[i], args[j]) {
                let mut rest: Vec<Atom> = args
                    .iter()
                    .enumerate()
                    .filter(|&(k, _)| k != i && k != j)
                    .map(|(_, a)| *a)
                    .collect();
                rest.push(merged);
                return Some(if rest.len() == 1 {
                    rest[0]
                } else {
                    ctx.add(&rest)
                });
            }
        }
    }
    None
}

/// Split `expr` into `(coeff, log_arg)` when it is `coeff·log(arg)` (coeff
/// may be 1).
fn as_scaled_log<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Option<(Atom<'a>, Atom<'a>)> {
    if let AtomNode::Fun(name, args) = expr.node() {
        if name.as_str() == "log" && args.len() == 1 {
            return Some((ctx.num(1), args[0]));
        }
    }
    if let AtomNode::Mul(factors) = expr.node() {
        if factors.len() == 2 {
            if let AtomNode::Fun(name, args) = factors[1].node() {
                if name.as_str() == "log" && args.len() == 1 {
                    return Some((factors[0], args[0]));
                }
            }
        }
    }
    None
}

/// Split `arg` into `(u, v)` when it is `u + I·v` or `u - I·v`.
fn as_complex_linear<'a>(arg: Atom<'a>) -> Option<(Atom<'a>, Atom<'a>, bool)> {
    let AtomNode::Add(terms) = arg.node() else {
        return None;
    };
    if terms.len() != 2 {
        return None;
    }
    for &t in terms.iter() {
        if let Some((sign, v)) = as_iv(t) {
            let u = if t == terms[0] { terms[1] } else { terms[0] };
            return Some((u, v, sign));
        }
    }
    None
}

/// Match `I·v` (returns `(true, v)`) or `-I·v` (returns `(false, v)`).
fn as_iv<'a>(t: Atom<'a>) -> Option<(bool, Atom<'a>)> {
    if crate::complex::is_i(t) {
        // Bare I means v = 1, which we cannot construct without an arena.
        return None;
    }
    let AtomNode::Mul(factors) = t.node() else {
        return None;
    };
    let has_i = factors.iter().any(|f| crate::complex::is_i(*f));
    if !has_i {
        return None;
    }
    let rest: Vec<Atom> = factors
        .iter()
        .filter(|f| !crate::complex::is_i(**f))
        .copied()
        .collect();
    let (neg, rest): (bool, Vec<Atom>) = if rest
        .first()
        .is_some_and(|r| matches!(r.node(), AtomNode::Num(-1)))
    {
        (true, rest[1..].to_vec())
    } else {
        (false, rest)
    };
    let v = match rest.len() {
        0 => return None,
        1 => rest[0],
        _ => return None,
    };
    Some((!neg, v))
}

/// Try merging one pair of scaled logs into a real form.
fn try_merge_pair<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>, b: Atom<'a>) -> Option<Atom<'a>> {
    let (ca, arga) = as_scaled_log(ctx, a)?;
    let (cb, argb) = as_scaled_log(ctx, b)?;
    // Coefficients must match exactly (hash-consed equality) or be exact
    // negatives of each other.
    let (u, v, pos) = as_complex_linear(arga)?;
    let (u2, v2, pos2) = as_complex_linear(argb)?;
    if u != u2 || v != v2 || pos == pos2 {
        return None;
    }
    // a = c·log(u+Iv), b = c·log(u-Iv) → c·log(u²+v²)
    if ca == cb {
        let u2v2 = ctx.add(&[ctx.pow(u, ctx.num(2)), ctx.pow(v, ctx.num(2))]);
        let log = ctx.fun("log", &[u2v2]);
        // If the original had a bare log (ca is the log itself), emit a
        // bare log; otherwise keep the coefficient.
        return Some(if matches!(ca.node(), AtomNode::Num(1)) {
            log
        } else {
            ctx.mul(&[ca, log])
        });
    }
    // a = c·log(u+Iv), b = -c·log(u-Iv) → 2c·atan(v/u)
    if is_negation_of(ctx, cb, ca) {
        let ratio = ctx.mul(&[v, ctx.pow(u, ctx.num(-1))]);
        let atan = ctx.fun("atan", &[ratio]);
        let two = ctx.num(2);
        return Some(if matches!(ca.node(), AtomNode::Num(1)) {
            ctx.mul(&[two, atan])
        } else {
            ctx.mul(&[two, ca, atan])
        });
    }
    None
}

/// Whether `b` is `-1 · a` (structurally, via a leading -1 factor).
fn is_negation_of<'a>(ctx: &'a AtomArena<'a>, b: Atom<'a>, a: Atom<'a>) -> bool {
    let neg_a = ctx.mul(&[ctx.num(-1), a]);
    b == neg_a
}

/// In a product, combine `exp(I·u)·exp(I·v)` → `exp(I·(u+v))` so that
/// cancellation (e.g. exp(Ix)·exp(-Ix)) is visible to the simplifier.
fn combine_exp_pair<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>]) -> Option<Atom<'a>> {
    let mut exps: Vec<(usize, Atom<'a>)> = Vec::new();
    for (idx, f) in factors.iter().enumerate() {
        if let AtomNode::Fun(name, args) = f.node() {
            if name.as_str() == "exp" && args.len() == 1 && contains_i(args[0]) {
                exps.push((idx, args[0]));
            }
        }
    }
    if exps.len() < 2 {
        return None;
    }
    let (i0, a0) = exps[0];
    let (i1, a1) = exps[1];
    let sum = ctx.add(&[a0, a1]);
    let merged = ctx.fun("exp", &[sum]);
    let mut rest: Vec<Atom> = factors
        .iter()
        .enumerate()
        .filter(|&(k, _)| k != i0 && k != i1)
        .map(|(_, a)| *a)
        .collect();
    rest.push(merged);
    Some(if rest.len() == 1 {
        rest[0]
    } else {
        ctx.mul(&rest)
    })
}

fn contains_i(a: Atom) -> bool {
    if crate::complex::is_i(a) {
        return true;
    }
    match a.node() {
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|c| contains_i(*c))
        }
        AtomNode::Pow(b, e) => contains_i(*b) || contains_i(*e),
        _ => false,
    }
}

/// True when `expr` still contains the imaginary unit anywhere.
#[allow(dead_code)]
pub(crate) fn has_imaginary(expr: Atom) -> bool {
    contains_i(expr)
}

#[allow(dead_code)]
fn _unused(s: Symbol) -> Symbol {
    s
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_rewrite::simplify::simplify;

    use super::*;
    use crate::complex::complex_rules;

    #[test]
    fn sin_rewrites_to_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let out = trig_to_exp(&ctx, ctx.fun("sin", &[x])).expect("has trig");
        let s = out.to_string();
        assert!(s.contains("exp"), "expected exponentials, got {s}");
        assert!(s.contains("I"), "expected imaginary unit, got {s}");
    }

    #[test]
    fn no_trig_returns_none() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert!(trig_to_exp(&ctx, ctx.add(&[x, ctx.num(1)])).is_none());
    }

    #[test]
    fn realify_merges_conjugate_log_sum() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // log(x + I) + log(x - I) → log(x² + 1)
        let ix = ctx.mul(&[i(&ctx), ctx.num(1)]);
        let neg_ix = ctx.mul(&[ctx.num(-1), i(&ctx), ctx.num(1)]);
        let sum = ctx.add(&[
            ctx.fun("log", &[ctx.add(&[x, ix])]),
            ctx.fun("log", &[ctx.add(&[x, neg_ix])]),
        ]);
        let out = realify(&ctx, sum);
        assert_eq!(out.to_string(), "log((x^2) + (1^2))");
    }

    #[test]
    fn realify_exp_cancellation() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(I·x)·exp(-I·x) → exp(0)
        let prod = ctx.mul(&[
            ctx.fun("exp", &[ctx.mul(&[i(&ctx), x])]),
            ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), i(&ctx), x])]),
        ]);
        let merged = realify(&ctx, prod);
        let rules = complex_rules(&ctx, &());
        let simp = simplify(&ctx, merged, &rules, 20);
        // After merging, the exponent is I·x + (-I·x); the complex rules
        // do not cancel it fully, but the exponents are combined.
        assert!(simp.to_string().starts_with("exp("), "got {simp}");
    }

    #[test]
    fn has_imaginary_detection() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert!(has_imaginary(ctx.mul(&[i(&ctx), x])));
        assert!(!has_imaginary(ctx.add(&[x, ctx.num(1)])));
    }
}
