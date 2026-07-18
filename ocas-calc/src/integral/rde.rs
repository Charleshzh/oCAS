//! The Risch differential equation `D q + f·q = g` (Bronstein ch. 6).
//!
//! Only polynomial solutions are sought: at each tower level the unknown
//! `q` must be a polynomial in the top generator, and at the base level a
//! polynomial in `x`. Rational-function solutions (requiring denominator
//! bounds) are outside this fragment and cause a fallback to the caller's
//! unevaluated form.
//!
//! The solver recurses down the tower: at a primitive level the
//! coefficients of `q` are eliminated top-down, at a hyperexponential
//! level the degree layers decouple, and at the base `ℚ[x]` level a
//! degree-bound + coefficient elimination closes the recursion.

use ocas_atom::AtomArena;
use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

use super::rational::poly_integrate;
use crate::tower::build::{GenKind, Tower};
use crate::tower::elem::{KElem, KPoly};

type DPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Solve `D q + f·q = g` for `q ∈ k_ℓ`, where `k_ℓ` is the field at tower
/// `level` and `f` does not involve the top generator. Returns `None`
/// when no solution exists in the implemented polynomial fragment.
pub(crate) fn rde_solve<'a>(
    _ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    f: &KElem,
    g: &KElem,
) -> Option<KElem> {
    if g.is_zero() {
        return Some(KElem::zero(f.n_vars()));
    }
    if level == 0 {
        return base_rde_kelem(f, g);
    }

    let top = level;
    let tgen = &tower.gens[level - 1];

    // f must be free of the top generator.
    if f.num.degree_in(top) > 0 || f.den.degree_in(top) > 0 {
        return None;
    }
    // g must reduce to a polynomial in the top generator. Rather than
    // checking the raw denominator (which may carry a cancellable factor
    // such as x/x), multiply through and let the KPoly view decide.
    let n = f.n_vars();
    let g_den_inv = KElem::from_poly(g.den.clone()).inv()?;
    let g_poly = KPoly::from_sparse(&g.num, top).mul_kelem(&g_den_inv);
    // After reduction the KPoly view may still have non-polynomial
    // coefficients; that is caught later by the coefficient extractions.
    let m = g_poly.degree()?;

    let mut q_coeffs = vec![KElem::zero(n); m + 1];
    match tgen.kind {
        GenKind::Constant | GenKind::Log => {
            // Primitive: D(q_j t^j) = Dq_j·t^j + j·q_j·t^{j-1}·Dt.
            // Layer j: Dq_j + f·q_j = g_j - (j+1)·q_{j+1}·Dt.
            for j in (0..=m).rev() {
                let mut rhs = g_poly.coeff_at(j);
                if j < m {
                    let shift = q_coeffs[j + 1]
                        .mul_rational(&Rational::new((j + 1) as i64, 1))
                        .mul(&tgen.dt);
                    rhs = rhs.sub(&shift);
                }
                q_coeffs[j] = rde_solve(_ctx, tower, level - 1, f, &rhs)?;
            }
        }
        GenKind::Exp => {
            // Hyperexponential (Dt = Du·t): layers decouple.
            // Layer j: Dq_j + (f + j·Du)·q_j = g_j.
            let t = KElem::var(top, n);
            let du = tgen.dt.div(&t)?;
            for j in (0..=m).rev() {
                let fj = f.add(&du.mul_rational(&Rational::new(j as i64, 1)));
                q_coeffs[j] = rde_solve(_ctx, tower, level - 1, &fj, &g_poly.coeff_at(j))?;
            }
        }
    }

    let q = KPoly {
        top,
        coeffs: q_coeffs,
        n_vars: n,
    };
    Some(q.kelem())
}

/// Base level `k₀ = ℚ(x)`: convert to dense univariate and solve.
fn base_rde_kelem(f: &KElem, g: &KElem) -> Option<KElem> {
    let n = f.n_vars();
    let fd = kelem_to_dpoly(f)?;
    let gd = kelem_to_dpoly(g)?;
    let q = base_rde(&fd, &gd)?;
    Some(embed(&dpoly_to_kelem(&q), n))
}

/// Embed a field element into a larger polynomial ring (unused trailing
/// variables), keeping exponent vectors valid.
fn embed(e: &KElem, n: usize) -> KElem {
    if e.n_vars() == n {
        return e.clone();
    }
    let embed_poly = |p: &crate::tower::elem::Sparse| -> crate::tower::elem::Sparse {
        crate::tower::elem::Sparse::from_terms(
            RationalDomain,
            n,
            p.terms_ref()
                .iter()
                .map(|(exp, c)| (exp.to_vec(), c.clone()))
                .collect(),
        )
    };
    KElem::new(embed_poly(&e.num), embed_poly(&e.den))
}

/// Field element ↔ dense polynomial over `ℚ[x]` (constant denominators).
fn kelem_to_dpoly(e: &KElem) -> Option<DPoly> {
    let dom = RationalDomain;
    // Case A: e is a rational constant (possibly unreduced, e.g. t/t).
    if let Some(c) = e.as_rational() {
        return Some(DPoly::from_coeffs(RationalDomain, vec![c]));
    }
    // Case B: constant scalar denominator and numerator in x only.
    let den_is_const = e
        .den
        .terms_ref()
        .keys()
        .all(|ex| ex.iter().all(|&k| k == 0));
    if !den_is_const {
        return None;
    }
    let dc = e.den.coeff(&vec![0; e.n_vars()]);
    let dc_inv = dom.inv(&dc)?;
    if e.num
        .terms_ref()
        .keys()
        .any(|ex| ex.iter().skip(1).any(|&k| k != 0))
    {
        return None;
    }
    let deg = e.num.degree_in(0);
    let mut coeffs = vec![dom.zero(); deg + 1];
    for (exp, c) in e.num.terms_ref() {
        coeffs[exp[0]] = dom.mul(c, &dc_inv);
    }
    Some(DPoly::from_coeffs(RationalDomain, coeffs))
}

fn dpoly_to_kelem(p: &DPoly) -> KElem {
    let terms = p
        .coeffs()
        .iter()
        .enumerate()
        .filter(|&(_, c)| !RationalDomain.is_zero(c))
        .map(|(i, c)| (vec![i], c.clone()))
        .collect();
    KElem::from_poly(crate::tower::elem::Sparse::from_terms(
        RationalDomain,
        1,
        terms,
    ))
}

/// Polynomial RDE over `ℚ[x]`: solve `q' + f·q = g` for `q ∈ ℚ[x]`.
fn base_rde(f: &DPoly, g: &DPoly) -> Option<DPoly> {
    let dom = RationalDomain;
    if g.is_zero() {
        return Some(DPoly::from_coeffs(RationalDomain, vec![]));
    }
    if f.is_zero() {
        // q' = g: termwise integration (constant of integration = 0).
        return Some(poly_integrate(g));
    }
    let mf = f.degree()?;
    let mg = g.degree()?;
    let mut r = g.clone();

    if mf == 0 {
        // f = c ≠ 0: deg q = deg g; eliminate top-down via q_j = r_j / c.
        let c = f.lcoeff();
        let mut q_coeffs = vec![dom.zero(); mg + 1];
        for j in (0..=mg).rev() {
            let rc = coeff_at(&r, j);
            if dom.is_zero(&rc) {
                continue;
            }
            let qj = dom.div(&rc, &c)?;
            q_coeffs[j] = qj.clone();
            let term = monomial(qj, j);
            r = r.sub(&term.derivative().add(&term.mul_scalar(&c)));
        }
        if r.is_zero() {
            Some(DPoly::from_coeffs(RationalDomain, q_coeffs))
        } else {
            None
        }
    } else {
        // deg(f·q) = deg q + mf dominates deg(q') = deg q - 1, so
        // deg q = mg - mf uniquely (no nonzero homogeneous solutions).
        let m = mg.checked_sub(mf)?;
        let flc = f.lcoeff();
        let mut q_coeffs = vec![dom.zero(); m + 1];
        for j in (0..=m).rev() {
            let rc = coeff_at(&r, j + mf);
            if dom.is_zero(&rc) {
                continue;
            }
            let qj = dom.div(&rc, &flc)?;
            q_coeffs[j] = qj.clone();
            let term = monomial(qj, j);
            r = r.sub(&term.derivative().add(&term.mul(f)));
        }
        if r.is_zero() {
            Some(DPoly::from_coeffs(RationalDomain, q_coeffs))
        } else {
            None
        }
    }
}

fn coeff_at(p: &DPoly, i: usize) -> Rational {
    p.coeffs()
        .get(i)
        .cloned()
        .unwrap_or_else(|| RationalDomain.zero())
}

fn monomial(c: Rational, k: usize) -> DPoly {
    let mut coeffs = vec![RationalDomain.zero(); k];
    coeffs.push(c);
    DPoly::from_coeffs(RationalDomain, coeffs)
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;
    use crate::tower::build::build_tower;
    use crate::tower::convert::atom_to_rational_extended;
    use ocas_atom::Symbol;

    fn rat(p: i64, q: i64) -> Rational {
        Rational::new(p, q)
    }

    fn dpoly(coeffs: &[(i64, i64)]) -> DPoly {
        DPoly::from_coeffs(
            RationalDomain,
            coeffs.iter().map(|&(p, q)| rat(p, q)).collect(),
        )
    }

    #[test]
    fn base_rde_constant_f() {
        // q' + q = x → q = x - 1
        let f = dpoly(&[(1, 1)]);
        let g = dpoly(&[(0, 1), (1, 1)]);
        let q = base_rde(&f, &g).expect("solution");
        assert_eq!(q, dpoly(&[(-1, 1), (1, 1)]));
    }

    #[test]
    fn base_rde_polynomial_f() {
        // q' + x·q = x^2 + 1 → q = x
        let f = dpoly(&[(0, 1), (1, 1)]);
        let g = dpoly(&[(1, 1), (0, 1), (1, 1)]);
        let q = base_rde(&f, &g).expect("solution");
        assert_eq!(q, dpoly(&[(0, 1), (1, 1)]));
    }

    #[test]
    fn base_rde_no_solution() {
        // q' + x·q = 1: deg q = 0 - 1 < 0 → no polynomial solution.
        let f = dpoly(&[(0, 1), (1, 1)]);
        let g = dpoly(&[(1, 1)]);
        assert!(base_rde(&f, &g).is_none());
    }

    #[test]
    fn base_rde_zero_f_integrates() {
        // q' = 2x → q = x^2
        let f = DPoly::from_coeffs(RationalDomain, vec![]);
        let g = dpoly(&[(0, 1), (2, 1)]);
        let q = base_rde(&f, &g).expect("solution");
        assert_eq!(q, dpoly(&[(0, 1), (0, 1), (1, 1)]));
    }

    #[test]
    fn rde_hyperexponential_layer() {
        // Tower [x, t = exp(x)]: solve Dq + q = x in k = ℚ(x, t) — the
        // answer is x - 1 (an element of ℚ(x) ⊂ k).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, ctx.fun("exp", &[x]), Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let one = KElem::one(n);
        let g = KElem::var(0, n);
        let q = rde_solve(&ctx, &tower, 1, &one, &g).expect("solution");
        let expect = KElem::var(0, n).sub(&KElem::one(n));
        assert!(q.eq_cross(&expect));
    }

    #[test]
    fn rde_primitive_layer() {
        // Tower [x, t = log(x)]: solve Dq + q = t in k = ℚ(x, t).
        // Try q = a·t + b: Dq + q = a'·t + a/x + b' + a·t + b
        //   = t·(a' + a) + (a/x + b' + b) = t → a' + a = 1 → a = 1;
        //   1/x + b' + b = 0 → b' + b = -1/x: no polynomial sol in ℚ(x)…
        // Use instead: Dq + q = t·(1) + (1/x) has q = t·? — simpler check:
        // Dq + q = 1 in k: q = 1 works.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, ctx.fun("log", &[x]), Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let one = KElem::one(n);
        let q = rde_solve(&ctx, &tower, 1, &one, &one).expect("solution");
        assert!(q.eq_cross(&KElem::one(n)));
    }

    #[test]
    fn rde_primitive_top_down() {
        // Tower [x, t = log(x)]: Dq + q = (x+1)·t + 1 has q = x·t:
        // D(x·t) = t + x·(1/x) = t + 1, so Dq + q = t + 1 + x·t. ✓
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, ctx.fun("log", &[x]), Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let f = KElem::one(n);
        let g = KElem::var(0, n)
            .add(&KElem::one(n))
            .mul(&KElem::var(1, n))
            .add(&KElem::one(n));
        let q = rde_solve(&ctx, &tower, 1, &f, &g).expect("solution");
        // Forward check: Dq + q == g.
        let dq = crate::tower::build::tower_diff(&q, &tower.gens);
        assert!(dq.add(&q).eq_cross(&g));
        // And q == x·t.
        assert!(q.eq_cross(&KElem::var(0, n).mul(&KElem::var(1, n))));
    }

    #[test]
    fn kelem_dpoly_roundtrip() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(3)]);
        let rf = atom_to_rational_extended(expr, &[x], 1).unwrap();
        let e = KElem::new(rf.numerator, rf.denominator);
        let d = kelem_to_dpoly(&e).unwrap();
        assert_eq!(d, dpoly(&[(3, 1), (0, 1), (1, 1)]));
        let back = dpoly_to_kelem(&d);
        assert!(back.eq_cross(&e));
    }
}
