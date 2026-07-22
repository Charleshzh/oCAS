//! Cross-domain property test: forward-mode AD via [`HyperDual`] matches
//! symbolic differentiation via [`ocas_calc::diff`].
//!
//! We check the headline 0.18.0 target: full first-order partials of a
//! three-variable product. The symbolic side is evaluated at the sample point
//! with a tiny inline evaluator (Num/Var/Add/Mul/Pow only) so the test does
//! not depend on `ocas-eval`'s floating-point trait bounds.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_calc::diff;
use ocas_core::arena::Arena;
use ocas_domain::dual::{HyperDual, new_first_order};
use ocas_domain::{Domain, Rational, RationalDomain};
use proptest::prelude::*;

/// Hand-rolled evaluator over the subset of atoms produced by `diff` for a
/// polynomial in `x`, `y`, `z`: Num, Var, Add, Mul, Pow (with integer exp).
fn eval_rational(atom: &Atom<'_>, x: &Rational, y: &Rational, z: &Rational) -> Rational {
    match atom.node() {
        AtomNode::Num(n) => Rational::new(*n, 1),
        AtomNode::Var(s) => match s.as_str() {
            "x" => x.clone(),
            "y" => y.clone(),
            "z" => z.clone(),
            other => panic!("unexpected variable in diff output: {other}"),
        },
        AtomNode::Add(children) => {
            let dom = RationalDomain;
            children.iter().fold(dom.zero(), |acc, c| {
                dom.add(&acc, &eval_rational(c, x, y, z))
            })
        }
        AtomNode::Mul(children) => {
            let dom = RationalDomain;
            children.iter().fold(dom.one(), |acc, c| {
                dom.mul(&acc, &eval_rational(c, x, y, z))
            })
        }
        AtomNode::Pow(base, exp) => {
            let b = eval_rational(base, x, y, z);
            match exp.node() {
                AtomNode::Num(e) => {
                    let dom = RationalDomain;
                    if *e >= 0 {
                        dom.pow(&b, *e as u64)
                    } else {
                        // b^(-n) = 1 / b^n
                        let pn = dom.pow(&b, (-*e) as u64);
                        dom.div(&dom.one(), &pn).expect("nonzero base")
                    }
                }
                _ => panic!("non-integer exponent in diff output"),
            }
        }
        AtomNode::Fun(_, _) => panic!("unexpected function node in polynomial diff output"),
    }
}

/// Build `f = x * y * z` in the supplied arena.
fn build_xyz<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    let x = ctx.var("x");
    let y = ctx.var("y");
    let z = ctx.var("z");
    ctx.mul(&[x, y, z])
}

fn small_rational() -> impl Strategy<Value = Rational> {
    (1i64..10).prop_map(|i| Rational::new(i, 1))
}

proptest! {
    /// ∂(xyz)/∂x at (x0,y0,z0): HyperDual AD must equal symbolic diff = y·z.
    #[test]
    fn dual_matches_diff_xyz_partial_x(
        x0 in small_rational(),
        y0 in small_rational(),
        z0 in small_rational(),
    ) {
        // --- HyperDual path ---
        let shape = new_first_order::<Rational>(3);
        let x = HyperDual::variable(&shape, 0, x0.clone());
        let y = HyperDual::variable(&shape, 1, y0.clone());
        let z = HyperDual::variable(&shape, 2, z0.clone());
        let f = x * y * z;
        let ad = f.deriv(0).expect("x-derivative slot present").clone();

        // --- Symbolic path ---
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = build_xyz(&ctx);
        let d = diff(&ctx, expr, Symbol::new("x"));
        let sym = eval_rational(&d, &x0, &y0, &z0);

        prop_assert_eq!(ad, sym);
    }

    /// ∂(xyz)/∂y must equal x·z on both paths.
    #[test]
    fn dual_matches_diff_xyz_partial_y(
        x0 in small_rational(),
        y0 in small_rational(),
        z0 in small_rational(),
    ) {
        let shape = new_first_order::<Rational>(3);
        let x = HyperDual::variable(&shape, 0, x0.clone());
        let y = HyperDual::variable(&shape, 1, y0.clone());
        let z = HyperDual::variable(&shape, 2, z0.clone());
        let f = x * y * z;
        let ad = f.deriv(1).expect("y-derivative slot present").clone();

        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = build_xyz(&ctx);
        let d = diff(&ctx, expr, Symbol::new("y"));
        let sym = eval_rational(&d, &x0, &y0, &z0);

        prop_assert_eq!(ad, sym);
    }

    /// ∂(xyz)/∂z must equal x·y on both paths.
    #[test]
    fn dual_matches_diff_xyz_partial_z(
        x0 in small_rational(),
        y0 in small_rational(),
        z0 in small_rational(),
    ) {
        let shape = new_first_order::<Rational>(3);
        let x = HyperDual::variable(&shape, 0, x0.clone());
        let y = HyperDual::variable(&shape, 1, y0.clone());
        let z = HyperDual::variable(&shape, 2, z0.clone());
        let f = x * y * z;
        let ad = f.deriv(2).expect("z-derivative slot present").clone();

        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = build_xyz(&ctx);
        let d = diff(&ctx, expr, Symbol::new("z"));
        let sym = eval_rational(&d, &x0, &y0, &z0);

        prop_assert_eq!(ad, sym);
    }
}
