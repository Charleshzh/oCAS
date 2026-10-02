//! The general logarithmic part at a tower level: Rothstein–Trager with
//! **rational roots only** (0.29.0; algebraic roots are 0.30.0's algebraic
//! extension wave).
//!
//! After Hermite reduction, the remainder `a1/d1` (with `d1` squarefree,
//! monic, coprime to `a1`) is integrated as `Σᵢ cᵢ·log(vᵢ)` where the `cᵢ`
//! are the roots of the Rothstein–Trager resultant
//! `R(z) = resultant_t(a1 − z·D d1, d1)` and `vᵢ = gcd(a1 − cᵢ·D d1, d1)`.
//! This replaces the 0.28.0 restriction to the logarithmic-derivative
//! identity `a1 == c·D d1` (which stays as a fast path at the call site —
//! it is exactly the single-root case `v = d1`).
//!
//! Only roots in `ℚ` are accepted, and only when they are *complete*
//! (deflating `R` by every rational linear factor found must leave a
//! constant): a leftover factor means irrational/complex constant roots
//! (0.30.0) or a non-constant root (the integral is then non-elementary in
//! `k(t)`), which 0.29.0 cannot distinguish — both decline honestly.
//!
//! Formulas are cross-checked against SymPy 1.14 `risch.py`
//! (`residue_reduce`) and Bronstein ch. 5; see
//! `docs/planning/LOGPART_RT_SPEC_CN.md` for the full specification with
//! line references.

use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

use crate::tower::build::{Tower, tower_diff_kpoly};
use crate::tower::elem::{KElem, KPoly};

type DPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Degree cap on `d1` for the resultant computation (deterministic budget:
/// the interpolation below costs `deg d1 + 1` polynomial resultants, and
/// each resultant is exact field arithmetic with no gcd normalization —
/// sizes grow quickly).
const MAX_RT_DEGREE: usize = 8;

/// Cap on rational-root candidates sieved from the integer contents.
const MAX_RT_CANDIDATES: usize = 64;

/// Absolute-value cap for trial-division factoring of the integer content.
/// Larger contents would only enlarge the candidate set; declining is the
/// honest outcome (the completeness check turns a missed root into a
/// refusal, never into a wrong answer).
const MAX_CONTENT: i64 = 1_000_000_000;

/// Integrate the simple fraction `a1/d1` at `level` as `Σ cᵢ·log(vᵢ)` plus
/// a remainder in the coefficient field `k`.
///
/// Returns `None` (the caller emits the unevaluated `Integral` fallback)
/// when `R(z)` has no complete set of rational roots, when any internal
/// budget trips, or when the self-verification fails. On success the
/// remainder is guaranteed `t`-free; it is zero except at hyperexponential
/// levels, where `D log v` contributes a `k`-valued rational part.
pub(crate) fn rothstein_trager(
    tower: &Tower,
    level: usize,
    a1: &KPoly,
    d1: &KPoly,
) -> Option<(Vec<(Rational, KElem)>, KElem)> {
    let top = level;
    let n = d1.n_vars;
    let deg_d = d1.degree()?;
    if deg_d == 0 || deg_d > MAX_RT_DEGREE {
        return None;
    }
    // Work with a monic denominator: `a1/d1 = (a1/lc)/(d1/lc)`. The tower
    // derivative of the monic form is computed directly below (`lc ∈ k` is
    // not necessarily `D`-constant, so `D(d1/lc) ≠ (D d1)/lc` in general).
    let lc = d1.lc();
    let (a1m, d1m) = if lc.eq_cross(&KElem::one(n)) {
        (a1.clone(), d1.clone())
    } else {
        let lc_inv = lc.inv()?;
        (a1.mul_kelem(&lc_inv), d1.monic())
    };
    // The simple-fraction precondition; the caller's `KRat` reduction
    // already ensures it, so a failure here means an upstream invariant
    // broke — decline rather than guess.
    if !d1m.gcd(&a1m).is_one() {
        return None;
    }
    let dd1m = clean_kpoly(&tower_diff_kpoly(
        &d1m,
        &tower.gens[..level - 1],
        &tower.gens[level - 1].dt,
    ))?;

    // R(z) by value interpolation: `deg_z R ≤ deg d1` (see the spec §1.2),
    // so `deg d1 + 1` rational nodes suffice.
    let r_poly = rt_resultant_polynomial(&a1m, &d1m, &dd1m)?;

    // Rational roots, complete or decline.
    let roots = rational_roots_complete(&r_poly)?;

    // One gcd per *distinct* root; `vᵢ` monic.
    let mut v_polys: Vec<KPoly> = Vec::with_capacity(roots.len());
    for c in &roots {
        let shifted = a1m.sub(&dd1m.mul_kelem(&KElem::from_rational(c, n)));
        let v = d1m.gcd(&shifted).monic();
        // `R(c) = 0` ⟺ `deg gcd ≥ 1`; a degree-0 gcd means an internal
        // inconsistency — decline honestly.
        if v.degree()? == 0 {
            return None;
        }
        v_polys.push(v);
    }

    // Self-verification (spec §1.4):
    // (i) `Π vᵢ == d1` exactly (distinct roots, `d1` squarefree).
    let mut prod = KPoly::one(top, n);
    for v in &v_polys {
        prod = prod.mul(v);
    }
    if !prod.sub(&d1m).is_zero() {
        return None;
    }
    // (ii) the remainder `a1/d1 − Σ cᵢ·D vᵢ/vᵢ` must be `t`-free.
    let mut rem = a1m.kelem().div(&d1m.kelem())?;
    let mut logs_out: Vec<(Rational, KElem)> = Vec::with_capacity(v_polys.len());
    for (c, v) in roots.iter().zip(v_polys.iter()) {
        let dv = tower_diff_kpoly(v, &tower.gens[..level - 1], &tower.gens[level - 1].dt);
        let term = dv.kelem().div(&v.kelem())?.mul_rational(c);
        rem = rem.sub(&term);
        logs_out.push((c.clone(), v.kelem()));
    }
    if rem.num.degree_in(top) > 0 || rem.den.degree_in(top) > 0 {
        return None;
    }

    Some((logs_out, rem))
}

/// Canonicalize a `KPoly` whose coefficients mention the top variable back
/// to coefficients in `k` — or `None` when the field element is not a
/// polynomial in `t`.
///
/// `tower_diff_kpoly` produces the polluted form at hyperexponential
/// levels: `D t = Du·t` puts the field element `t` into the coefficient
/// slots (`Σ (cᵢ·t)·tⁱ` instead of `Σ cᵢ·t^{i+1}`). The field value is
/// correct either way, but degree/`gcd`/`div_rem` are structural — they
/// need the canonical coefficients-in-`k` form.
pub(crate) fn clean_kpoly(p: &KPoly) -> Option<KPoly> {
    let top = p.top;
    let e = p.kelem();
    if e.den.degree_in(top) > 0 {
        return None;
    }
    let out = KPoly::from_sparse(&e.num, top);
    let den_inv = KElem::from_poly(e.den.clone()).inv()?;
    Some(out.mul_kelem(&den_inv))
}

/// `R(z) = resultant_t(a1 − z·D d1, d1)` as a dense coefficient vector in
/// `z` (ascending), built by value interpolation.
///
/// The resultant is taken with the **formal degree template**
/// `m = max(deg a1, deg D d1)` for the first argument — i.e. over `ℚ(z)`,
/// where the `z`-carrying leading coefficient is generically nonzero. At a
/// node `j` where the leading coefficient of `a1 − j·D d1` cancels (at
/// most one node: the leading-coefficient ratio), the naive resultant's
/// Sylvester matrix shrinks and specialization *fails* (measured: for
/// `a1 = 1, d1 = 1 + t, D d1 = t`, node `0` gives `Res(1, 1+t) = 1` while
/// the formal-template polynomial `R(z) = −(1+z)` requires `R(0) = −1`).
/// Such nodes are skipped; `deg d1 + 1` good nodes always exist and suffice
/// because `deg_z R ≤ deg d1`.
pub(crate) fn rt_resultant_polynomial(a1: &KPoly, d1: &KPoly, dd1: &KPoly) -> Option<Vec<KElem>> {
    let n = d1.n_vars;
    let deg_d = d1.degree()?;
    let template_deg = a1.degree().unwrap_or(0).max(dd1.degree().unwrap_or(0));
    let mut nodes: Vec<i64> = Vec::with_capacity(deg_d + 1);
    let mut values: Vec<KElem> = Vec::with_capacity(deg_d + 1);
    let mut j: i64 = 0;
    while values.len() < deg_d + 1 {
        if j > 4 * (deg_d as i64 + 2) {
            // Defensive: at most one node can drop the leading term, so
            // this is unreachable; decline rather than loop forever.
            return None;
        }
        let u = a1.sub(&dd1.mul_kelem(&KElem::from_rational(&Rational::new(j, 1), n)));
        j += 1;
        if u.degree() != Some(template_deg) {
            continue; // leading-coefficient cancellation: specialization fails here
        }
        nodes.push(j - 1);
        values.push(u.resultant(d1));
    }
    // Lagrange interpolation over the (distinct, integer) nodes:
    // `[z^k] R = Σ_j y_j · [z^k] Π_{i≠j} (z − z_i) / Π_{i≠j} (z_j − z_i)`.
    let mut out = vec![KElem::zero(n); deg_d + 1];
    for (idx, y) in values.iter().enumerate() {
        let zj = nodes[idx];
        let mut basis = vec![Rational::new(1, 1)];
        let mut denom: i64 = 1;
        for (i, &zi) in nodes.iter().enumerate() {
            if i == idx {
                continue;
            }
            // basis *= (z − z_i)
            let mut next = vec![RationalDomain.zero(); basis.len() + 1];
            for (k, c) in basis.iter().enumerate() {
                let term = RationalDomain.mul(c, &Rational::new(zi, 1));
                next[k] = RationalDomain.sub(&next[k], &term);
                next[k + 1] = RationalDomain.add(&next[k + 1], c);
            }
            basis = next;
            denom = denom.checked_mul(zj - zi)?;
        }
        for (k, c) in basis.iter().enumerate() {
            let scale = RationalDomain.div(c, &Rational::new(denom, 1))?;
            out[k] = out[k].add(&y.mul_rational(&scale));
        }
    }
    // Trim trailing zeros.
    while out.last().is_some_and(|c| c.is_zero()) {
        out.pop();
    }
    if out.is_empty() {
        // R ≡ 0 would mean `a1 − z·D d1` and `d1` share a factor for every
        // z — impossible for a simple fraction with `deg d1 ≥ 1`; decline.
        return None;
    }
    Some(out)
}

/// Rational-root candidates for `R(z)` (coefficients in the tower field).
///
/// On the fast path (`R ∈ ℚ[z]`) these are the exact rational roots from
/// integer factorization; otherwise a finite sieve from the rational-root
/// theorem over `ℤ[x, t₁, …][z]` (Gauss: `p | content(a₀)`,
/// `q | content(aₘ)`). `None` means the sieve's budgets tripped. Callers
/// must still evaluate each candidate exactly.
pub(crate) fn rational_root_candidates(r: &[KElem]) -> Option<Vec<Rational>> {
    // Fast path: `R ∈ ℚ[z]` — the exact factorization gives the roots
    // themselves (whether or not they are complete).
    if let Some(qcoeffs) = r
        .iter()
        .map(|c| c.as_rational())
        .collect::<Option<Vec<_>>>()
    {
        let (roots, _fully_split) =
            super::rational::rational_roots(&DPoly::from_coeffs(RationalDomain, qcoeffs))?;
        return Some(roots);
    }

    // General path: rational-root theorem over `ℤ[x, t₁, …][z]`. Clearing
    // denominators gives `A(z) = Σ aᵢ·zⁱ` with `aᵢ ∈ ℤ[vars]`; a reduced
    // rational root `p/q` satisfies `p | content(a₀)` and `q | content(aₘ)`
    // (Gauss), so the candidate set is finite.
    let coeffs: Vec<(&crate::tower::elem::Sparse, &crate::tower::elem::Sparse)> =
        r.iter().map(|e| (&e.num, &e.den)).collect();
    // L = lcm of every ℚ-denominator appearing in any numerator/denominator.
    let mut lcm: i64 = 1;
    for (num, den) in &coeffs {
        for c in num.terms_ref().values().chain(den.terms_ref().values()) {
            let d = c.denom().to_i64()?;
            lcm = num_integer::lcm(lcm, d);
            if lcm > MAX_CONTENT {
                return None;
            }
        }
    }
    let content = |p: &crate::tower::elem::Sparse| -> Option<i64> {
        // Integer content of `p·L` (L clears all denominators by construction).
        let mut g: i64 = 0;
        for c in p.terms_ref().values() {
            let v = RationalDomain.mul(c, &Rational::new(lcm, 1));
            let v = v.numer().to_i64()?;
            g = num_integer::gcd(g, v.checked_abs()?);
        }
        Some(g)
    };
    let m = r.len() - 1;
    // content(a'₀) = content(n₀·L)·Π_{j≥1} content(dⱼ·L); content(a'ₘ)
    // symmetric (Gauss: content of a product is the product of contents).
    let mut g0 = content(coeffs[0].0)?;
    let mut gm = content(coeffs[m].0)?;
    for (i, (_, den)) in coeffs.iter().enumerate() {
        let cd = content(den)?;
        if i != 0 {
            g0 = g0.checked_mul(cd)?;
        }
        if i != m {
            gm = gm.checked_mul(cd)?;
        }
    }
    if g0 == 0 || gm == 0 || g0 > MAX_CONTENT || gm > MAX_CONTENT {
        return None;
    }
    // Candidates: ±p/q with p | g0, q | gm (deduplicated; the set is small).
    let mut candidates: Vec<Rational> = Vec::new();
    let divs = |g: i64| -> Vec<i64> {
        let mut out = Vec::new();
        let mut d = 1i64;
        while d <= g / d {
            if g % d == 0 {
                out.push(d);
                if d != g / d {
                    out.push(g / d);
                }
            }
            d += 1;
        }
        out
    };
    for p in divs(g0) {
        for q in divs(gm) {
            for sign in [1, -1] {
                let c = Rational::new(sign * p, q);
                if !candidates.contains(&c) {
                    candidates.push(c);
                }
            }
        }
    }
    if candidates.len() > MAX_RT_CANDIDATES {
        return None;
    }
    Some(candidates)
}

/// All distinct rational roots of `R(z)` (coefficients in the tower field),
/// or `None` when the root set is not provably a complete set of rationals.
pub(crate) fn rational_roots_complete(r: &[KElem]) -> Option<Vec<Rational>> {
    // Deflate factors of z first: a zero root contributes
    // `gcd(a1 − 0·D d1, d1) = gcd(a1, d1) = 1` (simple-fraction
    // precondition), so it never produces a log term.
    let mut r: Vec<KElem> = r.to_vec();
    while r.first().is_some_and(|c| c.is_zero()) {
        r.remove(0);
    }
    if r.is_empty() {
        return Some(Vec::new());
    }

    let candidates = rational_root_candidates(&r)?;

    // Evaluate exactly and deflate per found root (with multiplicity), then
    // require a constant remainder (completeness).
    let mut cur = r.clone();
    let mut roots: Vec<Rational> = Vec::new();
    for c in candidates {
        loop {
            // Horner evaluation of `cur` at `c`.
            let mut acc = cur.last()?.clone();
            for i in (0..cur.len() - 1).rev() {
                acc = acc.mul_rational(&c).add(&cur[i]);
            }
            if !acc.is_zero() {
                break;
            }
            // Deflate by (z − c): synthetic division on ascending coeffs.
            let mut quo = vec![KElem::zero(cur[0].n_vars()); cur.len() - 1];
            let mut carry = cur.last()?.clone();
            for i in (1..cur.len()).rev() {
                quo[i - 1] = carry.clone();
                carry = cur[i - 1].add(&carry.mul_rational(&c));
            }
            cur = quo;
            if !roots.contains(&c) {
                roots.push(c.clone());
            }
            if cur.is_empty() {
                break;
            }
        }
    }
    // Completeness: the deflated polynomial must be a nonzero constant.
    if cur.len() != 1 || cur[0].is_zero() {
        return None;
    }
    Some(roots)
}

#[cfg(test)]
mod tests {
    use ocas_atom::{Atom, AtomArena, Symbol};
    use ocas_core::arena::Arena;

    use super::*;
    use crate::tower::build::build_tower;

    /// Build the tower for `expr`, convert it to a `KRat` at the top level,
    /// and split off the Hermite remainder `(a1, d1)` as in `integrate_level`.
    fn simple_remainder<'a>(
        ctx: &'a AtomArena<'a>,
        expr: Atom<'a>,
    ) -> (Tower<'a>, usize, KPoly, KPoly) {
        let tower = build_tower(ctx, expr, Symbol::new("x")).expect("tower");
        let level = tower.gens.len();
        let rf = crate::tower::convert::atom_to_rational_extended(
            tower.expr,
            &tower.gen_atoms(),
            tower.n_vars(),
        )
        .expect("field element");
        let f = crate::tower::elem::KRat::new(
            KPoly::from_sparse(&rf.numerator, level),
            KPoly::from_sparse(&rf.denominator, level),
        );
        let (_p, r) = f.num.div_rem(&f.den);
        // The remainders in the test cases below are already
        // squarefree-denominator; mirror `integrate_level`'s Hermite call.
        let (_g, a1, d1) =
            crate::integral::risch::hermite_tower(&tower, level, &r, &f.den).expect("hermite");
        (tower, level, a1, d1)
    }

    #[test]
    fn rt_exp_level_single_root() {
        // ∫ dx/(1 + exp(x)): a1 = 1, d1 = 1 + t, D d1 = t → root −1, v = 1 + t.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let f = ctx.pow(ctx.add(&[ctx.num(1), ctx.fun("exp", &[x])]), ctx.num(-1));
        let (tower, level, a1, d1) = simple_remainder(&ctx, f);
        let (logs, rem) = rothstein_trager(&tower, level, &a1, &d1).expect("RT");
        assert_eq!(logs.len(), 1);
        assert_eq!(logs[0].0, Rational::new(-1, 1));
        // Remainder must be t-free (it is exactly 1).
        assert!(rem.num.degree_in(level) == 0 && rem.den.degree_in(level) == 0);
    }

    #[test]
    fn rt_exp_level_double_root_once() {
        // ∫ dx/(exp(2x) − 1): d1 = t² − 1 (t = exp(2x), D t = 2t),
        // R = (1 − 2z)² → single root 1/2 counted once.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two_x = ctx.mul(&[ctx.num(2), x]);
        let f = ctx.pow(
            ctx.add(&[ctx.fun("exp", &[two_x]), ctx.num(-1)]),
            ctx.num(-1),
        );
        let (tower, level, a1, d1) = simple_remainder(&ctx, f);
        let (logs, rem) = rothstein_trager(&tower, level, &a1, &d1).expect("RT");
        assert_eq!(logs.len(), 1);
        assert_eq!(logs[0].0, Rational::new(1, 2));
        assert!(rem.num.degree_in(level) == 0 && rem.den.degree_in(level) == 0);
    }

    #[test]
    fn rt_log_level_two_roots() {
        // ∫ dx/(x·log(x)·(log(x)+1)): R = (1 − z²)/x² — coefficients carry
        // the lower generator x, exercising the general-candidate path.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let log_x = ctx.fun("log", &[x]);
        let den = ctx.mul(&[x, log_x, ctx.add(&[log_x, ctx.num(1)])]);
        let f = ctx.pow(den, ctx.num(-1));
        let (tower, level, a1, d1) = simple_remainder(&ctx, f);
        let (logs, rem) = rothstein_trager(&tower, level, &a1, &d1).expect("RT");
        assert_eq!(logs.len(), 2);
        let cs: Vec<Rational> = logs.iter().map(|(c, _)| c.clone()).collect();
        assert!(cs.contains(&Rational::new(-1, 1)) && cs.contains(&Rational::new(1, 1)));
        // At a log level the remainder must vanish exactly.
        assert!(rem.is_zero());
    }

    #[test]
    fn rt_declines_a_non_constant_root() {
        // ∫ dx/(exp(x) + x): R(z) = z·(1−x) − 1, root 1/(1−x) ∉ ℚ — the
        // integral is non-elementary; decline honestly.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let f = ctx.pow(ctx.add(&[ctx.fun("exp", &[x]), x]), ctx.num(-1));
        let (tower, level, a1, d1) = simple_remainder(&ctx, f);
        assert!(rothstein_trager(&tower, level, &a1, &d1).is_none());
    }
}
