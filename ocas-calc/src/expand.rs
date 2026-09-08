//! Bounded distributive expansion for the integration fallback path.
//!
//! Distributes multiplication over addition (`a*(b + c)` → `a*b + a*c`) and
//! expands small integer powers of sums, so that a failed product integral
//! can be retried termwise. The expansion is budgeted: when the number of
//! produced sum terms would exceed [`MAX_TERMS`], the attempt is abandoned
//! (`None`) so pathological inputs cannot blow up the integration chain.
//!
//! Unlike `normalize`, this pass can make expressions strictly larger; it is
//! therefore only used as a last-resort retry after the direct integration
//! methods (rational/Risch/rules/heuristics) have declined the input.

use ocas_atom::{Atom, AtomArena, AtomNode};

/// Maximum number of sum terms the expansion may produce in total.
const MAX_TERMS: usize = 64;

/// Expand `expr` distributively, returning `None` when the budget is
/// exceeded or the expression is already fully distributed (no change).
pub(crate) fn expand_bounded<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Option<Atom<'a>> {
    let mut terms = 0usize;
    let out = expand_rec(ctx, expr, &mut terms)?;
    if out == expr { None } else { Some(out) }
}

fn expand_rec<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, terms: &mut usize) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Add(args) => {
            let mut mapped = Vec::with_capacity(args.len());
            for a in args.iter() {
                mapped.push(expand_rec(ctx, *a, terms)?);
            }
            Some(ctx.add(&mapped))
        }
        AtomNode::Mul(args) => {
            let mut factors = Vec::with_capacity(args.len());
            for a in args.iter() {
                factors.push(expand_rec(ctx, *a, terms)?);
            }
            // No sum factor: the product is already fully distributed.
            let Some(pos) = factors
                .iter()
                .position(|f| matches!(f.node(), AtomNode::Add(_)))
            else {
                return Some(ctx.mul(&factors));
            };
            let AtomNode::Add(add_args) = factors[pos].node() else {
                unreachable!()
            };
            *terms += add_args.len();
            if *terms > MAX_TERMS {
                return None;
            }
            let mut result_terms = Vec::with_capacity(add_args.len());
            for term in add_args.iter() {
                let mut new_factors = factors.clone();
                new_factors[pos] = *term;
                result_terms.push(expand_rec(ctx, ctx.mul(&new_factors), terms)?);
            }
            Some(ctx.add(&result_terms))
        }
        AtomNode::Pow(base, exp) => {
            // (a + b + ...)^n with small integer n expands by repeated
            // multiplication; pre-check the worst-case term count.
            if let AtomNode::Num(n) = exp.node()
                && (2..=8).contains(n)
                && let AtomNode::Add(add_args) = base.node()
            {
                let worst = add_args.len().checked_pow(*n as u32)?;
                *terms += worst;
                if *terms > MAX_TERMS {
                    return None;
                }
                let mut acc = *base;
                for _ in 1..*n {
                    acc = expand_rec(ctx, ctx.mul(&[acc, *base]), terms)?;
                }
                return Some(acc);
            }
            Some(expr)
        }
        _ => Some(expr),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn distributes_product_over_sum() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // x*(x + 1) -> x + x*x (normalize does not fold like factors)
        let expr = ctx.mul(&[x, ctx.add(&[x, ctx.num(1)])]);
        let expanded = expand_bounded(&ctx, expr).expect("expansion");
        let normalized = ocas_atom::normalize::normalize(&ctx, expanded);
        assert_eq!(normalized.to_string(), "x + (x*x)");
    }

    #[test]
    fn expands_square_of_sum() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // (x + 1)^2 -> 1 + x + x + x*x (unfolded like terms)
        let expr = ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(2));
        let expanded = expand_bounded(&ctx, expr).expect("expansion");
        let normalized = ocas_atom::normalize::normalize(&ctx, expanded);
        assert_eq!(normalized.to_string(), "1 + x + x + (x*x)");
    }

    #[test]
    fn unchanged_without_sum_factor() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[x, ctx.fun("sin", &[x])]);
        assert!(expand_bounded(&ctx, expr).is_none());
    }

    #[test]
    fn budget_rejects_huge_powers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // 4-term sum to the 8th power: 4^8 = 65536 terms, over budget.
        let sum = ctx.add(&[x, ctx.num(1), ctx.num(2), ctx.num(3)]);
        let expr = ctx.pow(sum, ctx.num(8));
        assert!(expand_bounded(&ctx, expr).is_none());
    }

    #[test]
    fn does_not_expand_function_arguments() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[ctx.mul(&[ctx.num(2), ctx.add(&[x, ctx.num(1)])])]);
        assert!(expand_bounded(&ctx, expr).is_none());
    }
}
