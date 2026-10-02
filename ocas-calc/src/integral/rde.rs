//! The Risch differential equation `D q + f·q = g` (Bronstein ch. 6) — the
//! full rational fragment (0.29.0).
//!
//! 0.28.0 and earlier only sought **polynomial** solutions; 0.29.0 ports the
//! complete pipeline (cross-checked function by function against SymPy
//! 1.14's `integrals/rde.py`, which mirrors Bronstein's theorems):
//!
//! 1. **Weak normalization** ([`weak_normalizer`], Thm 6.1.1): normalize `f`
//!    so no residue at a normal irreducible is a positive integer. The
//!    substitution relation used here is `z = q_wn·y`,
//!    `Dz + (f − D q_wn/q_wn)·z = q_wn·g` (derived from the product rule).
//! 2. **Normal denominator bound** ([`normal_denom`], Thm 6.1.2): compute
//!    `hn` with `q = y·hn` polynomial at normal places.
//! 3. **Special denominator bound** ([`special_denom`], Thm 6.2.1): the
//!    `t`-adic part at hyperexponential levels, including the
//!    `parametric_log_deriv` refinement (implemented at the base field
//!    only; higher levels conservatively skip it).
//! 4. **Degree bound** ([`bound_degree`], §6.3): the `limited_integrate`
//!    refinements of the cancellation cases are *not* ported — skipping
//!    them only lowers the bound (missed solutions become honest declines,
//!    never wrong answers).
//! 5. **SPDE** ([`spde`], Rothstein's algorithm): reduce `a·Dq + b·q = c`
//!    to `a` constant.
//! 6. **Polynomial RDE dispatch** ([`solve_poly_rde`]): `no_cancel_b_large`
//!    at every level; cancellation via [`cancel_exp`] / [`cancel_primitive`]
//!    recurses into [`rischde`] one level down (this is where rational
//!    solutions at lower levels get used). The `D q = C` case at
//!    hyperexponential levels is handled layerwise; primitive `b = 0` needs
//!    `is_deriv_in_field` (an integrator) and declines honestly.
//!
//! Everything outside the fragment returns `None` (the caller falls back);
//! the tower check in `risch.rs` verifies every emitted answer.
//!
//! Substitution bookkeeping: `z = q_wn·y`, `q = z·hn`, `r = q·hs`,
//! `r = alpha·y + beta`, so `y = (alpha·y_poly + beta)/(q_wn·hn·hs)`.

use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

use super::logpart::{clean_kpoly, rational_root_candidates, rt_resultant_polynomial};
use super::rational::poly_integrate;
use super::risch::{kpoly_monomial, t_adic_valuation};
use crate::tower::build::{GenKind, Tower, tower_diff, tower_diff_kpoly};
use crate::tower::elem::{KElem, KPoly};

type DPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Defensive cap on the SPDE round count (the mathematical termination
/// argument is degree-driven; the cap only guards against a porting bug).
const MAX_SPDE_ROUNDS: usize = 64;

/// Defensive cap on the polynomial degree bound.
const MAX_RDE_DEGREE_BOUND: i64 = 64;

/// Solve `D q + f·q = g` for `q ∈ k_ℓ`, where `k_ℓ` is the field at tower
/// `level`. Returns `None` when no solution exists in the implemented
/// fragment (module docs).
pub(crate) fn rde_solve(tower: &Tower, level: usize, f: &KElem, g: &KElem) -> Option<KElem> {
    rischde(tower, level, f, g)
}

/// The full rational RDE pipeline; see the module docs.
fn rischde(tower: &Tower, level: usize, f: &KElem, g: &KElem) -> Option<KElem> {
    if g.is_zero() {
        return Some(KElem::zero(f.n_vars()));
    }
    // 1. Weak normalization: `z = q_wn·y`, `Dz + f'·z = q_wn·g`.
    let (q_wn, f2) = weak_normalizer(tower, level, f)?;
    let g2 = g.mul(&q_wn);
    // 2. Normal denominator bound: `q2 = z·hn` satisfies `a·Dq2 + b·q2 = c`.
    let (a, b, c, hn) = normal_denom(tower, level, &f2, &g2)?;
    // 3. Special denominator bound: `r = q2·hs ∈ k[t]` satisfies
    //    `A·Dr + B·r = C`.
    let (a2, b2, c2, hs) = special_denom(tower, level, &a, &b, &c)?;
    // The coefficient side of the reduced equation must be polynomial in
    // `t` for the degree bound and SPDE (the denominator-bound theorems
    // guarantee it when a solution exists); decline otherwise.
    let b2p = as_kpoly(&b2, level)?;
    let c2p = as_kpoly(&c2, level)?;
    // 4. Degree bound.
    let n = bound_degree(tower, level, &a2, &b2p, &c2p)?;
    // 5. SPDE reduction.
    let (b3, c3, m, alpha, beta) = spde(tower, level, &a2, &b2p, &c2p, n)?;
    // 6. Polynomial solve.
    let y = if c3.is_zero() {
        KElem::zero(f.n_vars())
    } else {
        solve_poly_rde(tower, level, &b3, &c3, m)?
    };
    // 7. Reassemble: `r = alpha·y + beta`; `y_orig = r/(q_wn·hn·hs)`.
    let r = alpha.kelem().mul(&y).add(&beta.kelem());
    let hn_e = hn.kelem();
    r.div(&hs)?.div(&hn_e)?.div(&q_wn)
}

// ------------------------------------------------------------------
//  Level-specific helpers
// ------------------------------------------------------------------

/// The case tag of a level: `exp` levels have a special polynomial (`t`),
/// everything else in our towers is primitive-like (all irreducibles
/// normal; the base field included).
#[derive(Clone, Copy, PartialEq, Eq)]
enum LevelCase {
    Base,
    Primitive,
    Exp,
}

fn level_case(tower: &Tower, level: usize) -> LevelCase {
    if level == 0 {
        LevelCase::Base
    } else if tower.gens[level - 1].kind == GenKind::Exp {
        LevelCase::Exp
    } else {
        LevelCase::Primitive
    }
}

/// Full tower derivation of a `k[t]` polynomial, canonicalized back to
/// coefficients in `k` (`tower_diff_kpoly` alone leaves `t` inside
/// coefficients at hyperexponential levels — see `logpart::clean_kpoly`).
fn derivation(tower: &Tower, level: usize, p: &KPoly) -> Option<KPoly> {
    if level == 0 {
        // `D x = 1`: the formal derivative is already canonical.
        return Some(p.derivative_dt());
    }
    clean_kpoly(&tower_diff_kpoly(
        p,
        &tower.gens[..level - 1],
        &tower.gens[level - 1].dt,
    ))
}

/// Split a denominator into its normal and special parts. Only
/// hyperexponential levels have a nontrivial special polynomial, and there
/// it is exactly the monomial `t^ν` (the tower merge keeps generators
/// independent).
fn splitfactor(d: &KPoly, tower: &Tower, level: usize) -> (KPoly, KPoly) {
    if level_case(tower, level) == LevelCase::Exp
        && let Some(k) = t_adic_valuation(d)
        && k > 0
    {
        let normal = KPoly {
            top: d.top,
            coeffs: d.coeffs[k..].to_vec(),
            n_vars: d.n_vars,
        };
        return (normal, kpoly_monomial(d.top, d.n_vars, k));
    }
    (d.clone(), KPoly::one(d.top, d.n_vars))
}

/// The `t`-adic order of a field element (`+∞`-ish for zero), as
/// `valuation(numerator) − valuation(denominator)`.
fn order_t(e: &KElem, top: usize) -> i64 {
    fn val(p: &crate::tower::elem::Sparse, top: usize) -> i64 {
        if p.is_zero() {
            return i64::MAX / 4;
        }
        p.terms_ref()
            .keys()
            .map(|ex| ex[top] as i64)
            .min()
            .unwrap_or(0)
    }
    val(&e.num, top) - val(&e.den, top)
}

/// Evaluate a field element at `t = 0`; `None` when the element has a pole
/// there.
fn eval_t0(e: &KElem, top: usize) -> Option<KElem> {
    let keep = |p: &crate::tower::elem::Sparse| {
        crate::tower::elem::Sparse::from_terms(
            RationalDomain,
            p.n_vars(),
            p.terms_ref()
                .iter()
                .filter(|(ex, _)| ex[top] == 0)
                .map(|(ex, c)| (ex.to_vec(), c.clone()))
                .collect(),
        )
    };
    let den = keep(&e.den);
    if den.is_zero() {
        return None;
    }
    Some(KElem::new(keep(&e.num), den))
}

/// A field element as a `k[t]` polynomial, or `None` when it is not one.
fn as_kpoly(e: &KElem, top: usize) -> Option<KPoly> {
    if e.den.degree_in(top) > 0 {
        return None;
    }
    let den_inv = KElem::from_poly(e.den.clone()).inv()?;
    Some(KPoly::from_sparse(&e.num, top).mul_kelem(&den_inv))
}

// ------------------------------------------------------------------
//  Step 1: weak normalization (Bronstein Thm 6.1.1)
// ------------------------------------------------------------------

/// Compute `q_wn` such that `f − D q_wn/q_wn` is weakly normalized (no
/// positive-integer residue at a normal irreducible). Returns
/// `(q_wn, f − D q_wn/q_wn)`.
fn weak_normalizer(tower: &Tower, level: usize, f: &KElem) -> Option<(KElem, KElem)> {
    let top = level;
    let n = f.n_vars();
    let fa = KPoly::from_sparse(&f.num, top);
    let fd = KPoly::from_sparse(&f.den, top);
    let (dn, _) = splitfactor(&fd, tower, level);
    // `d1` = the product of the simple roots of `dn` (squarefree kernel
    // of the squarefree part).
    let g = dn.gcd(&dn.derivative_dt());
    let (d_sqf, _) = dn.div_rem(&g);
    let g2 = d_sqf.gcd(&g);
    let (d1, _) = d_sqf.div_rem(&g2);
    if d1.degree().unwrap_or(0) == 0 {
        return Some((KElem::one(n), f.clone()));
    }
    let dd1 = derivation(tower, level, &d1)?;
    // Residues of `f` at the roots of `d1` are the roots of
    // `R(z) = resultant_t(fa − z·D d1, d1)` (the Rothstein–Trager
    // construction, reused from `logpart`).
    let r_poly = rt_resultant_polynomial(&fa, &d1, &dd1)?;
    let candidates = rational_root_candidates(&r_poly)?;
    let mut q = KPoly::one(top, n);
    let mut any = false;
    for c in &candidates {
        // Positive integers only.
        if c.denom().to_i64() != Some(1) {
            continue;
        }
        let Some(numer) = c.numer().to_i64() else {
            continue;
        };
        if numer <= 0 {
            continue;
        }
        // Confirm the root exactly (sieve candidates are not pre-verified).
        let mut acc = r_poly.last()?.clone();
        for i in (0..r_poly.len() - 1).rev() {
            acc = acc.mul_rational(c).add(&r_poly[i]);
        }
        if !acc.is_zero() {
            continue;
        }
        let v = d1.gcd(&fa.sub(&dd1.mul_kelem(&KElem::from_rational(c, n))));
        q = q.mul(&v);
        any = true;
    }
    if !any {
        return Some((KElem::one(n), f.clone()));
    }
    let qe = q.kelem();
    let dq = tower_diff(&qe, &tower.gens[..level]);
    let fw = f.sub(&dq.div(&qe)?);
    Some((qe, fw))
}

// ------------------------------------------------------------------
//  Step 2: normal denominator bound (Bronstein Thm 6.1.2)
// ------------------------------------------------------------------

/// Given `D y + f·y = g` with `f` weakly normalized, compute
/// `(a, b, c, hn)` such that every solution `y` makes `q = y·hn` satisfy
/// `a·Dq + b·q = c`, with `a ∈ k[t]` and `q` polynomial at normal places.
/// `None` when a divisibility side-condition fails (the equation then has
/// no solution at all — Bronstein makes this `NonElementary`; we decline).
fn normal_denom(
    tower: &Tower,
    level: usize,
    f: &KElem,
    g: &KElem,
) -> Option<(KPoly, KElem, KElem, KPoly)> {
    let top = level;
    let fd = KPoly::from_sparse(&f.den, top);
    let gd = KPoly::from_sparse(&g.den, top);
    let (dn, _) = splitfactor(&fd, tower, level);
    let (en, _) = splitfactor(&gd, tower, level);
    let p = dn.gcd(&en);
    let g1 = en.gcd(&en.derivative_dt());
    let g2 = p.gcd(&p.derivative_dt());
    let (h, rem) = g1.div_rem(&g2);
    if !rem.is_zero() {
        return None;
    }
    let a = dn.mul(&h);
    // `c = a·h` must be divisible by `en` (else no solution).
    let c_poly = a.mul(&h);
    let (_, rem) = c_poly.div_rem(&en);
    if !rem.is_zero() {
        return None;
    }
    // c = a·h·g (as a field element).
    let c = c_poly.kelem().mul(g);
    // b = a·f − dn·D h.
    let dh = derivation(tower, level, &h)?;
    let b = a.kelem().mul(f).sub(&dn.kelem().mul(&dh.kelem()));
    Some((a, b, c, h))
}

// ------------------------------------------------------------------
//  Step 3: special denominator bound (Bronstein Thm 6.2.1)
// ------------------------------------------------------------------

/// The `t`-adic part of the denominator bound at a hyperexponential level:
/// with `ν = min(0, ν_c − min(0, ν_b))` and `N = max(0, −ν_b, ν − ν_c)`,
/// `r = q·t^ν ∈ k[t]` satisfies `A·Dr + B·r = C` with
/// `A = a·t^N`, `B = b + ν·a·(Dt/t)`, `C = c·t^(N−ν)`.
/// Primitive/base levels are the identity. Returns `(A, B, C, t^−ν)`.
fn special_denom(
    tower: &Tower,
    level: usize,
    a: &KPoly,
    b: &KElem,
    c: &KElem,
) -> Option<(KPoly, KElem, KElem, KElem)> {
    let top = level;
    let n = a.n_vars;
    if level_case(tower, level) != LevelCase::Exp {
        return Some((a.clone(), b.clone(), c.clone(), KElem::one(n)));
    }
    let tgen = &tower.gens[level - 1];
    let du = tgen.dt.div(&KElem::var(top, n))?;
    let nb = order_t(b, top);
    let nc = order_t(c, top);
    let mut nu = 0.min(nc - 0.min(nb));
    if nb == 0 {
        // Refinement (Bronstein's cancellation case at `t`): when
        // `−b(0)/a(0) = m·(Dt/t) + Dv/v` exactly, `ν = min(ν, m)`.
        let a0 = eval_t0(&a.kelem(), top);
        let b0 = eval_t0(b, top);
        if let (Some(a0), Some(b0)) = (a0, b0)
            && let Some(alpha) = b0.div(&a0).map(|q| q.neg())
            && let Some((1, m)) = parametric_log_deriv(tower, level - 1, &alpha, &du)
        {
            nu = nu.min(m);
        }
    }
    let big_n = 0.max(-nb).max(nu - nc);
    debug_assert!(big_n >= nu, "N ≥ ν always holds since ν ≤ 0");
    let big_a = a.mul(&kpoly_monomial(top, n, big_n as usize));
    // B = b + ν·a·(Dt/t)
    let big_b = b.add(&a.kelem().mul(&du).mul_rational(&Rational::new(nu, 1)));
    // C = c·t^(N−ν) — the exponent is non-negative because N ≥ ν.
    let big_c = c.mul(&KElem::var(top, n).pow((big_n - nu) as u64));
    let h = KElem::var(top, n).pow((-nu) as u64);
    Some((big_a, big_b, big_c, h))
}

// ------------------------------------------------------------------
//  Step 4: degree bound (Bronstein §6.3)
// ------------------------------------------------------------------

/// A bound `n` with `deg q ≤ n` for every polynomial solution of
/// `a·Dq + b·q = c`. The `limited_integrate` refinements of Bronstein's
/// cancellation cases are not ported (skipping them only lowers the bound:
/// missed solutions become declines, never wrong answers).
fn bound_degree(tower: &Tower, level: usize, a: &KPoly, b: &KPoly, c: &KPoly) -> Option<i64> {
    let da = a.degree().unwrap_or(0) as i64;
    // deg(0) conventions: −∞.
    let db = b.degree().map(|d| d as i64).unwrap_or(-1);
    let dc = c.degree().map(|d| d as i64).unwrap_or(-1);
    let n = match level_case(tower, level) {
        LevelCase::Base => {
            let mut n = 0.max(dc - db.max(da - 1));
            if db == da - 1 {
                // alpha = −lc(b)/lc(a) integer ⇒ cancellation possible.
                let alpha = b.lc().div(&a.lc())?.neg();
                if let Some(r) = alpha.as_rational()
                    && r.denom().to_i64() == Some(1)
                    && let Some(k) = r.numer().to_i64()
                {
                    n = n.max(k).max(dc - db);
                }
            }
            n
        }
        LevelCase::Primitive => {
            if db > da {
                0.max(dc - db)
            } else {
                0.max(dc - da + 1)
            }
        }
        LevelCase::Exp => 0.max(dc - db.max(da)),
    };
    if n > MAX_RDE_DEGREE_BOUND {
        return None;
    }
    Some(n)
}

// ------------------------------------------------------------------
//  Step 5: SPDE (Rothstein's special polynomial differential equation)
// ------------------------------------------------------------------

/// Reduce `a·Dq + b·q = c` (with the degree bound `n`) to an equation with
/// a constant `a`. Returns `(B, C, m, alpha, beta)` such that every
/// solution `q` of degree `≤ n` has the form `q = alpha·h + beta` with
/// `D h + B·h = C`, `deg h ≤ m`. `None` = no solution (or a defensive
/// budget tripped).
fn spde(
    tower: &Tower,
    level: usize,
    a0: &KPoly,
    b0: &KPoly,
    c0: &KPoly,
    n0: i64,
) -> Option<(KElem, KElem, i64, KPoly, KPoly)> {
    let top = a0.top;
    let nv = a0.n_vars;
    let (mut a, mut b, mut c) = (a0.clone(), b0.clone(), c0.clone());
    let mut alpha = KPoly::one(top, nv);
    let mut beta = KPoly::zero(top, nv);
    let mut n = n0;
    for _round in 0..MAX_SPDE_ROUNDS {
        if c.is_zero() {
            return Some((
                KElem::zero(nv),
                KElem::zero(nv),
                0,
                KPoly::zero(top, nv),
                beta,
            ));
        }
        if n < 0 {
            return None;
        }
        let g = a.gcd(&b);
        let (c_q, c_r) = c.div_rem(&g);
        if !c_r.is_zero() {
            // g ∤ c: no solution.
            return None;
        }
        a = a.div_rem(&g).0;
        b = b.div_rem(&g).0;
        c = c_q;
        if a.degree() == Some(0) {
            let lc_inv = a.lc().inv()?;
            return Some((
                b.mul_kelem(&lc_inv).kelem(),
                c.mul_kelem(&lc_inv).kelem(),
                n,
                alpha,
                beta,
            ));
        }
        // Bézout: `s·b + t·a = 1` (coprime after the division);
        // `r = c·s mod a`, and `c − r·b` is exactly divisible by `a` with
        // quotient `z`.
        let (g0, s, _t) = b.eea(&a);
        if !g0.is_one() {
            return None;
        }
        let (_, r) = c.mul(&s).div_rem(&a);
        let (z, zr) = c.sub(&r.mul(&b)).div_rem(&a);
        if !zr.is_zero() {
            return None;
        }
        let da = derivation(tower, level, &a)?;
        let dr = derivation(tower, level, &r)?;
        b = b.add(&da);
        c = z.sub(&dr);
        n -= a.degree()? as i64;
        beta = beta.add(&alpha.mul(&r));
        alpha = alpha.mul(&a);
    }
    None
}

// ------------------------------------------------------------------
//  Step 6: polynomial RDE dispatch
// ------------------------------------------------------------------

/// Solve `D q + b·q = c` for `q ∈ k[t]` with `deg q ≤ n`.
fn solve_poly_rde(tower: &Tower, level: usize, b: &KElem, c: &KElem, n: i64) -> Option<KElem> {
    let top = level;
    if level == 0 {
        // Base field ℚ(x): the dense univariate solver (covers `b = 0` by
        // termwise integration and `b ≠ 0` by the no-cancellation peel).
        let fd = kelem_to_dpoly(b)?;
        let gd = kelem_to_dpoly(c)?;
        let q = base_rde(&fd, &gd)?;
        return Some(embed(&dpoly_to_kelem(&q), b.n_vars()));
    }
    let bp = as_kpoly(b, top)?;
    let cp = as_kpoly(c, top)?;
    let case = level_case(tower, level);
    let deg_dt: i64 = if case == LevelCase::Exp { 1 } else { 0 };
    if !bp.is_zero() && bp.degree().unwrap_or(0) as i64 > 0.max(deg_dt - 1) {
        no_cancel_b_large(tower, level, &bp, &cp, n)
    } else if bp.is_zero() {
        match case {
            // `D q = c` at a hyperexponential level: layers decouple via
            // `D(r_m·t^m) = (D r_m + m·(Dt/t)·r_m)·t^m`.
            LevelCase::Exp => integrate_level_dr(tower, level, &cp, n),
            // Primitive `b = 0` needs `is_deriv_in_field` (an integrator):
            // honest decline.
            _ => None,
        }
    } else {
        // Cancellation: `deg b ≤ max(0, deg Dt − 1)`, i.e. `b ∈ k`.
        if bp.degree().unwrap_or(0) > 0 {
            return None;
        }
        match case {
            LevelCase::Exp => cancel_exp(tower, level, &bp, &cp, n),
            LevelCase::Primitive => cancel_primitive(tower, level, &bp, &cp, n),
            LevelCase::Base => None,
        }
    }
}

/// `D q + b·q = c` with `deg b` large enough that no cancellation occurs:
/// peel the leading term per round (`m = deg c − deg b`).
fn no_cancel_b_large(tower: &Tower, level: usize, b: &KPoly, c: &KPoly, n: i64) -> Option<KElem> {
    let top = b.top;
    let mut q = KPoly::zero(top, b.n_vars);
    let mut c = c.clone();
    let mut n = n;
    while !c.is_zero() {
        let m = c.degree()? as i64 - b.degree()? as i64;
        if !(0..=n).contains(&m) {
            return None;
        }
        let lc_ratio = c.lc().div(&b.lc())?;
        let p = KPoly::from_kelem(lc_ratio, top).monomial_shift(&KElem::one(b.n_vars), m as usize);
        q = q.add(&p);
        n = m - 1;
        let dp = derivation(tower, level, &p)?;
        c = c.sub(&dp).sub(&b.mul(&p));
    }
    Some(q.kelem())
}

/// `D q = c` at a hyperexponential level: `D(r_m·t^m) = (D r_m +
/// m·(Dt/t)·r_m)·t^m`, so the layers decouple into lower-level RDEs. The
/// `t⁰` layer needs an antiderivative in `k` — only the trivial case
/// `c₀ = 0` is supported (anything else is an integration problem, not an
/// RDE).
fn integrate_level_dr(tower: &Tower, level: usize, c: &KPoly, n: i64) -> Option<KElem> {
    let top = level;
    let tgen = &tower.gens[level - 1];
    let du = tgen.dt.div(&KElem::var(top, tower.n_vars()))?;
    let mut q = KPoly::zero(top, tower.n_vars());
    for (i, coeff) in c.coeffs.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        if i == 0 {
            // Would need `∫ c₀` in the field below: not an RDE. Decline.
            return None;
        }
        if n >= 0 && i as i64 > n {
            return None;
        }
        let f_m = du.mul_rational(&Rational::new(i as i64, 1));
        let s = rischde(tower, level - 1, &f_m, coeff)?;
        q = q.add(&KPoly::from_kelem(s, top).monomial_shift(&KElem::one(tower.n_vars()), i));
    }
    Some(q.kelem())
}

/// Cancellation at a hyperexponential level (`b ∈ k` nonzero): peel
/// `m = deg c` per round with `s` solving `D s + (b + m·Dt/t)·s = lc(c)`
/// one level down.
fn cancel_exp(tower: &Tower, level: usize, b: &KPoly, c: &KPoly, n: i64) -> Option<KElem> {
    let top = level;
    let tgen = &tower.gens[level - 1];
    let du = tgen.dt.div(&KElem::var(top, tower.n_vars()))?;
    let bv = b.lc();
    let mut q = KPoly::zero(top, tower.n_vars());
    let mut c = c.clone();
    let mut n = n;
    while !c.is_zero() {
        let m = c.degree()? as i64;
        if n < m {
            return None;
        }
        let f_lower = bv.add(&du.mul_rational(&Rational::new(m, 1)));
        let s = rischde(tower, level - 1, &f_lower, &c.lc())?;
        let stm = KPoly::from_kelem(s, top).monomial_shift(&KElem::one(tower.n_vars()), m as usize);
        q = q.add(&stm);
        n = m - 1;
        let dstm = derivation(tower, level, &stm)?;
        c = c.sub(&stm.mul_kelem(&bv)).sub(&dstm);
    }
    Some(q.kelem())
}

/// Cancellation at a primitive level (`b ∈ k` nonzero, `D t ∈ k`): peel
/// `m = deg c` per round with `s` solving `D s + b·s = lc(c)` one level
/// down.
fn cancel_primitive(tower: &Tower, level: usize, b: &KPoly, c: &KPoly, n: i64) -> Option<KElem> {
    let top = level;
    let bv = b.lc();
    let mut q = KPoly::zero(top, tower.n_vars());
    let mut c = c.clone();
    let mut n = n;
    while !c.is_zero() {
        let m = c.degree()? as i64;
        if n < m {
            return None;
        }
        let s = rischde(tower, level - 1, &bv, &c.lc())?;
        let stm = KPoly::from_kelem(s, top).monomial_shift(&KElem::one(tower.n_vars()), m as usize);
        q = q.add(&stm);
        n = m - 1;
        let dstm = derivation(tower, level, &stm)?;
        c = c.sub(&stm.mul_kelem(&bv)).sub(&dstm);
    }
    Some(q.kelem())
}

// ------------------------------------------------------------------
//  Parametric logarithmic derivative (base field only)
// ------------------------------------------------------------------

/// Decide whether `n·alpha = Dv/v + m·eta` for integers `n > 0`, `m` and `v`
/// in a radical extension of the field; returns `(n, m)`.
///
/// Only the base field `ℚ(x)` with constant `eta` is implemented (that is
/// the case the `special_denom` refinement needs at level 1); everything
/// else declines, which just skips the refinement.
///
/// Method: the polynomial part of `alpha` must be the constant `m·eta/n`;
/// the proper remainder must be `Dv/v`, i.e. its (reduced) denominator is
/// squarefree and every residue `e_p = (A·(Dp·B/p)⁻¹) mod p` modulo the
/// irreducible denominator factors is a rational constant. The reassembly
/// `alpha₂ == Σ e_p·Dp·(B/p) / B` is verified exactly before returning.
fn parametric_log_deriv(
    _tower: &Tower,
    level: usize,
    alpha: &KElem,
    eta: &KElem,
) -> Option<(i64, i64)> {
    if level != 0 {
        return None;
    }
    let dom = RationalDomain;
    let e = eta.as_rational()?;
    let an = sparse_to_dpoly_x(&alpha.num)?;
    let ad = sparse_to_dpoly_x(&alpha.den)?;
    // Polynomial part of `alpha` must be a constant `rho`.
    let (pp, _rem) = an.div_rem(&ad)?;
    let rho = if pp.is_zero() {
        RationalDomain.zero()
    } else {
        if pp.degree()? != 0 {
            return None;
        }
        pp.coeffs()[0].clone()
    };
    // `n·rho = m·e` in ℚ.
    let (n_i, m_i) = if dom.is_zero(&e) {
        if dom.is_zero(&rho) {
            (1, 0)
        } else {
            return None;
        }
    } else {
        let r = dom.div(&rho, &e)?;
        let n = r.denom().to_i64()?;
        let m = r.numer().to_i64()?;
        if n <= 0 { (-n, -m) } else { (n, m) }
    };
    // `alpha₂ = n·alpha − m·eta` must be `Dv/v`: zero polynomial part,
    // squarefree denominator, constant residues.
    let alpha2 = alpha
        .mul_rational(&Rational::new(n_i, 1))
        .sub(&eta.mul_rational(&Rational::new(m_i, 1)));
    if alpha2.is_zero() {
        return Some((n_i, m_i));
    }
    let a2n = sparse_to_dpoly_x(&alpha2.num)?;
    let a2d = sparse_to_dpoly_x(&alpha2.den)?;
    // Reduce A/B.
    let (g0, _, _) = a2n.extended_gcd_poly(&a2d);
    let (a2n, a2d) = if g0.is_zero() || g0.degree() == Some(0) {
        (a2n, a2d)
    } else {
        (a2n.div_rem(&g0)?.0, a2d.div_rem(&g0)?.0)
    };
    if a2d.is_zero() || a2d.degree()? == 0 {
        // A polynomial `alpha₂` is `Dv/v` only when it vanishes (handled).
        return None;
    }
    if a2n.degree().unwrap_or(0) >= a2d.degree()? {
        return None;
    }
    // Squarefree check.
    let d_a2d = a2d.derivative();
    let (sqg, _, _) = a2d.extended_gcd_poly(&d_a2d);
    if !sqg.is_zero() && sqg.degree()? > 0 {
        return None;
    }
    // Factor B over ℤ and compute the residue modulo each factor.
    let mut lcm: i64 = 1;
    for cf in a2d.coeffs() {
        lcm = num_integer::lcm(lcm, cf.denom().to_i64()?);
    }
    if lcm > 1_000_000 {
        return None;
    }
    let zcoeffs: Option<Vec<ocas_domain::Integer>> = a2d
        .coeffs()
        .iter()
        .map(|cf| {
            let scaled = dom.mul(cf, &Rational::new(lcm, 1));
            scaled.numer().to_i64().map(ocas_domain::Integer::from)
        })
        .collect();
    let zpoly = DenseUnivariatePolynomial::from_coeffs(ocas_domain::IntegerDomain, zcoeffs?);
    let factors = zpoly.primitive_part().factor();
    let mut sum_num = DPoly::from_coeffs(RationalDomain, vec![]);
    for (fac_i64, _mult) in &factors {
        if fac_i64.degree().unwrap_or(0) == 0 {
            continue;
        }
        // To ℚ, monic.
        let pcoeffs: Option<Vec<Rational>> = fac_i64
            .coeffs()
            .iter()
            .map(|cf| cf.to_i64().map(|v| Rational::new(v, 1)))
            .collect();
        let mut p = DPoly::from_coeffs(RationalDomain, pcoeffs?);
        let lc_inv = dom.inv(&p.lcoeff())?;
        p = p.mul_scalar(&lc_inv);
        let dp = p.derivative();
        let (bp, bp_rem) = a2d.div_rem(&p)?;
        if !bp_rem.is_zero() {
            return None;
        }
        // Inverse of `dp·bp` mod `p` via the extended gcd.
        let (_, w) = dp.mul(&bp).div_rem(&p)?;
        let (gw, sw, _) = w.extended_gcd_poly(&p);
        if gw.is_zero() || gw.degree()? > 0 {
            return None;
        }
        let gw_inv = dom.inv(&gw.lcoeff())?;
        let (_, e_poly) = a2n.mul(&sw.mul_scalar(&gw_inv)).div_rem(&p)?;
        if !e_poly.is_zero() && e_poly.degree()? != 0 {
            return None; // residue not a constant
        }
        let e_p = if e_poly.is_zero() {
            RationalDomain.zero()
        } else {
            e_poly.coeffs()[0].clone()
        };
        sum_num = sum_num.add(&dp.mul(&bp).mul_scalar(&e_p));
    }
    // Verify the reassembly: `alpha₂ == sum_num / a2d` — both sides share
    // the reduced denominator `a2d`, so compare numerators exactly.
    if sum_num != a2n {
        return None;
    }
    Some((n_i, m_i))
}

/// Dense univariate view of a sparse polynomial that uses only variable 0.
fn sparse_to_dpoly_x(p: &crate::tower::elem::Sparse) -> Option<DPoly> {
    if p.terms_ref()
        .keys()
        .any(|ex| ex.iter().skip(1).any(|&k| k != 0))
    {
        return None;
    }
    let dom = RationalDomain;
    let deg = p.degree_in(0);
    let mut coeffs = vec![dom.zero(); deg + 1];
    for (exp, c) in p.terms_ref() {
        coeffs[exp[0]] = c.clone();
    }
    Some(DPoly::from_coeffs(RationalDomain, coeffs))
}

// ------------------------------------------------------------------
//  Base level `k₀ = ℚ(x)` (dense univariate, from the 0.27.x fragment)
// ------------------------------------------------------------------

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
        let q = rde_solve(&tower, 1, &one, &g).expect("solution");
        let expect = KElem::var(0, n).sub(&KElem::one(n));
        assert!(q.eq_cross(&expect));
    }

    #[test]
    fn rde_primitive_layer() {
        // Tower [x, t = log(x)]: solve Dq + q = 1 in k: q = 1 works.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, ctx.fun("log", &[x]), Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let one = KElem::one(n);
        let q = rde_solve(&tower, 1, &one, &one).expect("solution");
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
        let q = rde_solve(&tower, 1, &f, &g).expect("solution");
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

    // --------------------------------------------------------------
    //  0.29.0: the rational fragment
    // --------------------------------------------------------------

    /// Verify `Dq + f·q == g` in the tower field.
    fn assert_rde_solution(tower: &Tower, level: usize, f: &KElem, g: &KElem, q: &KElem) {
        let dq = tower_diff(q, &tower.gens[..level]);
        let lhs = dq.add(&f.mul(q));
        assert!(lhs.eq_cross(g), "Dq + f·q != g: level {level}");
    }

    #[test]
    fn rde_rational_solution_at_exp_level() {
        // Tower [x, t = exp(x)]: Dq + (t+1)·q = 1 has the rational solution
        // q = 1/t (D(1/t) = −1/t). This is rubi-00992's inner equation —
        // the polynomial fragment cannot represent it.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, ctx.fun("exp", &[x]), Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let t = KElem::var(1, n);
        let f = t.add(&KElem::one(n));
        let g = KElem::one(n);
        let q = rde_solve(&tower, 1, &f, &g).expect("rational solution");
        assert_rde_solution(&tower, 1, &f, &g, &q);
        // q is exactly 1/t.
        let expect = KElem::one(n).div(&KElem::var(1, n)).unwrap();
        assert!(q.eq_cross(&expect));
    }

    #[test]
    fn rde_rational_solution_at_base_level() {
        // At ℚ(x): Dq + (2/x)·q = x → q = x²/4. `f` has a pole at x,
        // exercising weak normalization + the normal denominator bound.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, x, Symbol::new("x")).unwrap();
        let n = tower.n_vars();
        let xe = KElem::var(0, n);
        let f = KElem::from_rational(&rat(2, 1), n).div(&xe).unwrap();
        let g = xe.clone();
        let q = rde_solve(&tower, 0, &f, &g).expect("base rational RDE");
        assert_rde_solution(&tower, 0, &f, &g, &q);
        let expect = xe.mul(&xe).mul_rational(&rat(1, 4));
        assert!(q.eq_cross(&expect), "q = {q:?}, expected x²/4");
    }
}
