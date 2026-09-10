//! Normalization for [`Atom`] expression trees.
//!
//! The normalizer puts expressions into a deterministic canonical form:
//! nested additions and multiplications are flattened, arguments are sorted,
//! and numeric coefficients are merged.

use crate::{Atom, AtomArena, AtomNode};

/// Function heads whose argument order is semantic and must survive
/// normalization. Every other head has its arguments sorted into canonical
/// order.
///
/// - `Derivative` / `Integral` carry the variable of differentiation or
///   integration, which is positional.
/// - `EllipticF` / `EllipticE` / `EllipticPi` carry an amplitude and a
///   parameter `m = k²` (plus a characteristic for `EllipticPi`) that are not
///   interchangeable; sorting them would make `EllipticF(u, m)` and
///   `EllipticF(m, u)` the same atom.
/// - `Ei` is the exponential integral; its two-argument form is Rubi's
///   `Ei(n, z)` (the incomplete gamma `Eₙ`), whose argument order is semantic.
pub fn preserves_argument_order(name: &str) -> bool {
    matches!(
        name,
        "Derivative" | "Integral" | "EllipticF" | "EllipticE" | "EllipticPi" | "Ei"
    )
}

/// Normalize an atom into canonical form.
///
/// The result is allocated in the same arena as the input via `ctx`.
///
/// # Example
///
/// ```
/// use ocas_atom::normalize::normalize;
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let y = ctx.var("y");
/// let z = ctx.var("z");
/// let inner = ctx.add(&[x, y]);
/// let outer = ctx.add(&[inner, z, ctx.num(2), ctx.num(3)]);
/// let result = normalize(&ctx, outer);
/// assert_eq!(result.to_string(), "5 + x + y + z");
/// ```
pub fn normalize<'a>(ctx: &AtomArena<'a>, atom: Atom<'a>) -> Atom<'a> {
    match atom.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => atom,
        AtomNode::Fun(name, args) => {
            let mut normalized: Vec<Atom<'a>> = args.iter().map(|a| normalize(ctx, *a)).collect();
            // Preserve argument order for forms where order is semantic.
            if !preserves_argument_order(name.as_str()) {
                normalized.sort();
            }
            ctx.fun(name.as_str(), &normalized)
        }
        AtomNode::Add(args) => {
            // Normalize children FIRST, then flatten — this ensures any child
            // that normalizes into an Add node gets flattened, guaranteeing
            // idempotency (normalize(normalize(x)) == normalize(x)).
            let normalized_children: Vec<Atom<'a>> =
                args.iter().map(|a| normalize(ctx, *a)).collect();
            let mut flat = Vec::new();
            collect_add(&normalized_children, &mut flat);
            let mut normalized = flat;
            // Drop explicit zero terms first (covers the common `x + 0` case),
            // then sort and merge numeric literals. Merging can itself produce
            // a new zero (e.g. `93 + -93`), so drop zeros AGAIN after merging.
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(0)));
            normalized.sort();
            merge_numbers(ctx, &mut normalized, true);
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(0)));
            if normalized.is_empty() {
                ctx.num(0)
            } else if normalized.len() == 1 {
                normalized[0]
            } else {
                ctx.add(&normalized)
            }
        }
        AtomNode::Mul(args) => {
            // Normalize children FIRST, then flatten — same reasoning as Add.
            let normalized_children: Vec<Atom<'a>> =
                args.iter().map(|a| normalize(ctx, *a)).collect();
            let mut flat = Vec::new();
            collect_mul(&normalized_children, &mut flat);
            let mut normalized = flat;
            if normalized
                .iter()
                .any(|a| matches!(a.node(), AtomNode::Num(0)))
            {
                return ctx.num(0);
            }
            // Drop explicit unit terms first, then sort and merge numeric
            // literals. Merging can produce a new unit (e.g. `-1 * -1 = 1`),
            // so drop units AGAIN after merging — mirrors the Add branch.
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(1)));
            normalized.sort();
            merge_numbers(ctx, &mut normalized, false);
            normalized.retain(|a| !matches!(a.node(), AtomNode::Num(1)));
            if normalized.is_empty() {
                ctx.num(1)
            } else if normalized.len() == 1 {
                normalized[0]
            } else {
                ctx.mul(&normalized)
            }
        }
        AtomNode::Pow(base, exp) => {
            let base = normalize(ctx, *base);
            let exp = normalize(ctx, *exp);
            // Exact numeric folds: `u^0 → 1`, `u^1 → u`, `0^n → 0` (n > 0),
            // `1^_ → 1`, `(-1)^-1 → -1`, and exact integer powers `b^e`.
            if let AtomNode::Num(e) = exp.node() {
                if *e == 0 {
                    return ctx.num(1);
                }
                if *e == 1 {
                    return base;
                }
                if let AtomNode::Num(b) = base.node() {
                    if *b == 0 {
                        if *e > 0 {
                            return ctx.num(0);
                        }
                    } else if *b == 1 {
                        return ctx.num(1);
                    } else if *e > 0 {
                        if let Ok(e32) = u32::try_from(*e)
                            && let Some(v) = b.checked_pow(e32)
                        {
                            return ctx.num(v);
                        }
                    } else if *e == -1 && *b == -1 {
                        return ctx.num(-1);
                    }
                }
            }
            // Fold (u^r)^n → u^(r·n) when the outer exponent n is an
            // integer (r rational or atom): sound formal power-of-a-power
            // folding; non-integer outer exponents stay unfolded
            // ((x²)^(1/2) ≠ x).
            if let AtomNode::Num(n) = exp.node()
                && let AtomNode::Pow(inner_base, inner_exp) = base.node()
            {
                let merged = normalize(ctx, ctx.mul(&[*inner_exp, ctx.num(*n)]));
                return ctx.pow(*inner_base, merged);
            }
            ctx.pow(base, exp)
        }
    }
}

fn collect_add<'a>(args: &[Atom<'a>], out: &mut Vec<Atom<'a>>) {
    for &arg in args {
        match arg.node() {
            AtomNode::Add(inner) => collect_add(inner, out),
            _ => out.push(arg),
        }
    }
}

fn collect_mul<'a>(args: &[Atom<'a>], out: &mut Vec<Atom<'a>>) {
    for &arg in args {
        match arg.node() {
            AtomNode::Mul(inner) => collect_mul(inner, out),
            _ => out.push(arg),
        }
    }
}

fn merge_numbers<'a>(ctx: &AtomArena<'a>, args: &mut Vec<Atom<'a>>, is_add: bool) {
    let count = args
        .iter()
        .take_while(|a| matches!(a.node(), AtomNode::Num(_)))
        .count();

    if count >= 2 {
        let nums: Vec<i64> = args[0..count]
            .iter()
            .map(|a| match a.node() {
                AtomNode::Num(n) => *n,
                _ => unreachable!(),
            })
            .collect();
        // Use wrapping arithmetic to avoid panics on overflow in debug mode.
        // This matches Rust's release-mode behavior for i64 arithmetic.
        let merged = if is_add {
            nums.into_iter().fold(0i64, |acc, n| acc.wrapping_add(n))
        } else {
            nums.into_iter().fold(1i64, |acc, n| acc.wrapping_mul(n))
        };
        args.drain(0..count);
        args.insert(0, ctx.num(merged));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn normalize_leaves_atom_unchanged() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert_eq!(normalize(&ctx, x).to_string(), "x");
    }

    #[test]
    fn normalize_flattens_nested_add() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let inner = ctx.add(&[x, y]);
        let outer = ctx.add(&[inner, z]);
        assert_eq!(normalize(&ctx, outer).to_string(), "x + y + z");
    }

    #[test]
    fn normalize_drops_zero_from_opposite_numerics() {
        // `93 + (-93) + sin(x)` must collapse to `sin(x)`: the two numerics
        // merge to 0, which must then be dropped (regression for the
        // retain-before-merge ordering bug found by proptest).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sinx = ctx.fun("sin", &[x]);
        let a1 = ctx.add(&[ctx.num(93)]);
        let a2 = ctx.add(&[ctx.num(-93)]);
        let atom = ctx.add(&[a1, a2, sinx]);
        assert_eq!(normalize(&ctx, atom).to_string(), "sin(x)");
    }

    #[test]
    fn normalize_drops_unit_from_opposite_numerics() {
        // `(-1) * ((-1) * x)` must collapse to `x`: the two units merge to 1,
        // which must then be dropped.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let neg1 = ctx.num(-1);
        let inner = ctx.mul(&[x, neg1]);
        let atom = ctx.mul(&[neg1, inner]);
        assert_eq!(normalize(&ctx, atom).to_string(), "x");
    }

    #[test]
    fn normalize_sorts_arguments() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let z = ctx.var("z");
        let expr = ctx.add(&[z, x, y]);
        assert_eq!(normalize(&ctx, expr).to_string(), "x + y + z");
    }

    #[test]
    fn normalize_merges_numeric_literals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let one = ctx.num(1);
        let two = ctx.num(2);
        let x = ctx.var("x");
        let expr = ctx.add(&[one, x, two]);
        assert_eq!(normalize(&ctx, expr).to_string(), "3 + x");
    }

    #[test]
    fn normalize_pow() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two = ctx.num(2);
        let pow = ctx.pow(x, two);
        assert_eq!(normalize(&ctx, pow).to_string(), "x^2");
    }

    #[test]
    fn normalize_sorts_fun_arguments() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let f = ctx.fun("f", &[y, x]);
        assert_eq!(normalize(&ctx, f).to_string(), "f(x, y)");
    }

    #[test]
    fn normalize_preserves_elliptic_argument_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        // Amplitude first, parameter second: the order is semantic.
        let f = ctx.fun("EllipticF", &[x, y]);
        assert_eq!(normalize(&ctx, f).to_string(), "EllipticF(x, y)");
        let g = ctx.fun("EllipticF", &[y, x]);
        assert_eq!(normalize(&ctx, g).to_string(), "EllipticF(y, x)");
        // Round-trip is idempotent.
        assert_eq!(
            normalize(&ctx, normalize(&ctx, g)).to_string(),
            "EllipticF(y, x)"
        );
        // Three-argument form keeps all three positions.
        let p = ctx.fun("EllipticPi", &[y, x, ctx.num(1)]);
        assert_eq!(normalize(&ctx, p).to_string(), "EllipticPi(y, x, 1)");
        // Ei's two-argument (En) form is order-sensitive too.
        let ei = ctx.fun("Ei", &[ctx.num(1), x]);
        assert_eq!(normalize(&ctx, ei).to_string(), "Ei(1, x)");
        // Ordinary heads keep sorting.
        let h = ctx.fun("sin", &[y]);
        assert_eq!(normalize(&ctx, h).to_string(), "sin(y)");
    }
}
