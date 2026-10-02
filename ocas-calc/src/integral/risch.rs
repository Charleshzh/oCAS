//! The Risch algorithm for elementary transcendental integrals.
//!
//! [`risch_integrate`] builds an elementary extension tower for the
//! integrand (see [`crate::tower::build`]) and integrates recursively
//! down the tower (Bronstein, *Symbolic Integration I*, ch. 5):
//!
//! - at each level, the rational part is split off by Hermite reduction
//!   in `k(t)`;
//! - the logarithmic part is handled when it matches the
//!   logarithmic-derivative identity `c·D d/d`;
//! - the polynomial part is integrated by undetermined coefficients at
//!   primitive (`log`) levels and by the Risch differential equation at
//!   hyperexponential (`exp`) levels (see [`super::rde`]);
//! - the base level `ℚ(x)` delegates to the rational-function integrator
//!   (see [`super::rational`]).
//!
//! Scope limits (all documented in the book chapter): only polynomial
//! solutions of the Risch differential equation are sought, and only the
//! logarithmic-derivative identity is used for logarithmic parts. Inputs
//! outside this fragment yield `None` (the caller falls back), except
//! individual terms that may be returned in the unevaluated form
//! `Integral(term, var)`.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, Symbol};
use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use super::rational::integrate_dpoly_structured;
use super::rde::rde_solve;
use crate::rules::calculus_rules;
use crate::tower::build::{GenKind, Tower, build_tower, tower_diff_kpoly};
use crate::tower::convert::{
    GeneratorField, atom_to_rational_extended, rational_const_to_atom, rational_to_atom,
};
use crate::tower::elem::{KElem, KPoly, KRat};

type DPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Maximum recursion depth for `risch_integrate` / `integrate_level`.
/// Prevents infinite loops when the Risch algorithm cannot handle an
/// integrand but keeps retrying (e.g. sec(x) forcing in VOP).
const MAX_RISCH_DEPTH: usize = 16;

// clippy 1.97 `missing_const_for_thread_local` false-positive: the
// initializer already is the canonical `const {}` form (same convention
// as ocas-atom/src/workspace.rs and ocas-poly/src/factor/eez.rs).
thread_local! {
    #[allow(clippy::missing_const_for_thread_local)]
    static RISCH_DEPTH: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

struct RischDepthGuard;

impl Drop for RischDepthGuard {
    fn drop(&mut self) {
        RISCH_DEPTH.set(RISCH_DEPTH.get().saturating_sub(1));
    }
}

/// Integrate `expr` (an elementary function of `var` built from rational
/// functions, `log`, and `exp`) by the Risch algorithm.
///
/// Returns `None` when the expression is outside the implemented fragment
/// (see module docs); the caller should fall back to other integrators.
pub(crate) fn risch_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let depth = RISCH_DEPTH.get();
    if depth >= MAX_RISCH_DEPTH {
        return None;
    }
    RISCH_DEPTH.set(depth + 1);
    let _guard = RischDepthGuard;
    let tower = build_tower(ctx, expr, var)?;
    let level = tower.gens.len();
    // The tower may have merged dependent generators and rewritten the
    // integrand over them; the field element must come from that form.
    let rf = atom_to_rational_extended(tower.expr, &tower.gen_atoms(), tower.n_vars())?;
    let f = KRat::new(
        KPoly::from_sparse(&rf.numerator, level),
        KPoly::from_sparse(&rf.denominator, level),
    );
    let result = integrate_level(ctx, &tower, level, &f)?;
    let raw = level_result_to_atom(ctx, &tower, level, &result)?;
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let after_default = simplify(ctx, raw, &default_rules, 20);
    let after_calc = simplify(ctx, after_default, &calc_rules, 10);
    let answer = normalize(ctx, after_calc);
    // 0.28.0: never emit a result the exact checker can falsify.
    //
    // The check is done **inside the tower's own field**: the generators are
    // independent by construction (the dependency merge in `tower::build`
    // makes sure of it), so `D(F) − f` converting to a nonzero field element
    // is a real wrong answer, not a prover gap. Measured: the exp-level
    // rational part mis-scales negative powers, and
    // `∫ (exp(x)²+1)³/(8·exp(x)⁴) dx` came out wrong by a factor of 2 and
    // 4/3 on its `t⁻²`/`t⁻⁴` terms — a latent solver bug that the hyperbolic
    // front-end made reachable from the corpus. Declining keeps the engine's
    // hard invariant ("no uncertified answer") while the solver bug is fixed
    // separately. When the difference cannot even be converted (no such
    // shape is known for a tower result), the answer is accepted as before.
    let verified = tower_difference_is_zero(ctx, answer, tower.expr, var, &tower);
    if verified == Some(false) {
        return None;
    }
    Some(answer)
}

/// Whether `D(answer) − integrand` is exactly zero in the tower field.
///
/// `None` when the difference is not a rational function of the tower
/// generators (the caller then accepts the answer, as before).
fn tower_difference_is_zero<'a>(
    ctx: &'a AtomArena<'a>,
    answer: Atom<'a>,
    integrand: Atom<'a>,
    var: Symbol,
    tower: &Tower<'a>,
) -> Option<bool> {
    let derived = crate::derivative::diff(ctx, answer, var);
    let minus_one = ctx.num(-1);
    let difference = ctx.add(&[derived, ctx.mul(&[minus_one, integrand])]);
    let gens = tower.gen_atoms();
    let element = atom_to_rational_extended(difference, &gens, tower.n_vars())?;
    Some(element.numerator.n_terms() == 0)
}

// ------------------------------------------------------------------
//  Structured per-level result
// ------------------------------------------------------------------

/// The integral computed at one tower level.
pub(crate) struct LevelResult<'a> {
    /// Rational part in `k(t)` (from Hermite reduction).
    pub rational: KRat,
    /// Polynomial part in `k[t]`.
    pub poly: KPoly,
    /// Constant part in `k`.
    pub constant: KElem,
    /// New logarithmic terms `c·log(v)` with `v` in the current field.
    pub logs: Vec<(Rational, KElem)>,
    /// Other terms (atan, unevaluated integrals) as atoms.
    pub extras: Vec<Atom<'a>>,
}

impl<'a> LevelResult<'a> {
    fn empty(top: usize, n: usize) -> Self {
        Self {
            rational: KRat::new(KPoly::zero(top, n), KPoly::one(top, n)),
            poly: KPoly::zero(top, n),
            constant: KElem::zero(n),
            logs: Vec::new(),
            extras: Vec::new(),
        }
    }

    /// Fold the element parts (`rational + poly + constant`) into a
    /// single coefficient-field element.
    fn elem_part(&self) -> KElem {
        self.rational
            .kelem()
            .add(&self.poly.kelem())
            .add(&self.constant)
    }
}

// ------------------------------------------------------------------
//  Level integrator
// ------------------------------------------------------------------

/// Integrate `f` (an element of `k_level(t_level)`) at the given level.
fn integrate_level<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    f: &KRat,
) -> Option<LevelResult<'a>> {
    if level == 0 {
        return integrate_level_base(ctx, tower, f);
    }

    let top = level;
    let n = tower.n_vars();
    let tgen = &tower.gens[level - 1];

    // Split off the polynomial part.
    let (mut p, r) = f.num.div_rem(&f.den);

    let mut out = LevelResult::empty(top, n);

    // 0.29.0: at a hyperexponential level `t` divides its own derivative
    // (`D t = Du·t`), so a `t`-factor of the denominator is *special*:
    // Hermite reduction requires `gcd(v, Dv) = 1` and silently mis-scales
    // the answer when that fails (the latent bug the 0.28.0 tower check
    // caught: `t⁻²`/`t⁻⁴` coefficients off by 2× and 4/3× in
    // `∫ (exp(x)²+1)³/(8·exp(x)⁴) dx`). Split `den = t^k·d′` with
    // `gcd(d′, t) = 1`; the `t^k` part integrates as a Laurent polynomial
    // (one RDE per negative layer, exactly like the positive layers in
    // [`integrate_kpoly_hyperexp`]), and only `d′` reaches Hermite.
    let (r_eff, d_eff) = if tgen.kind == GenKind::Exp
        && let Some(k) = t_adic_valuation(&f.den)
        && k > 0
    {
        let t = KElem::var(top, n);
        let du = tgen.dt.div(&t)?;
        let tk = kpoly_monomial(top, n, k);
        let d_prime = shift_down(&f.den, k);
        // Partial-fraction split over the coprime factors: with
        // `s·d′ + w·t^k = 1`, `r/(t^k·d′) = r·s/t^k + r·w/d′`.
        let (g0, s, w) = d_prime.eea(&tk);
        if !g0.is_one() {
            return None;
        }
        let (q_poly, low) = r.mul(&s).div_rem(&tk);
        p = p.add(&q_poly);
        // Laurent layers: `low/t^k = Σ_j low_j·t^(j−k)`; each layer `m = j−k`
        // is solved by `D q + m·Du·q = low_j` in the field below.
        for (j, c) in low.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let m = j as i64 - k as i64;
            let fj = du.mul_rational(&Rational::new(m, 1));
            let qj = rde_solve(tower, level - 1, &fj, c)?;
            let den_m = kpoly_monomial(top, n, k - j);
            out.rational = krat_add(&out.rational, &KRat::new(KPoly::from_kelem(qj, top), den_m));
        }
        let (p2, r2) = r.mul(&w).div_rem(&d_prime);
        p = p.add(&p2);
        (r2, d_prime)
    } else {
        (r, f.den.clone())
    };

    // Hermite reduction on the proper part (coprime to `t` at exp levels).
    let (g, a1, d1) = hermite_tower(tower, level, &r_eff, &d_eff)?;

    out.rational = krat_add(&out.rational, &g);

    // Logarithmic part: the logarithmic-derivative identity `a1 == c·D d1`
    // stays as the fast path (it is exactly the single-root case `v = d1`);
    // everything else goes through Rothstein–Trager (0.29.0, rational roots
    // only — anything else declines honestly).
    if !a1.is_zero() {
        let dd1 = tower_diff_kpoly(&d1, &tower.gens[..level - 1], &tgen.dt);
        if let Some(c) = kpoly_scalar_multiple(&a1, &dd1) {
            out.logs.push((c, d1.kelem()));
        } else if let Some((logs, rem)) = super::logpart::rothstein_trager(tower, level, &a1, &d1) {
            out.logs.extend(logs);
            if !rem.is_zero() {
                // Hyperexponential levels: `D log v` leaves a `k`-valued
                // remainder that integrates one level down.
                let res = integrate_kelem_or_fallback(ctx, tower, level - 1, rem)?;
                out.constant = out.constant.add(&res.elem_part());
                out.logs.extend(res.logs);
                out.extras.extend(res.extras);
            }
        } else {
            out.extras
                .push(integral_fallback_atom(ctx, tower, level, &a1, &d1)?);
        }
    }

    // Polynomial part.
    if !p.is_zero() {
        let (poly_ans, const_ans, logs_p, extras_p) = integrate_kpoly(ctx, tower, level, &p)?;
        out.poly = poly_ans;
        out.constant = const_ans;
        out.logs.extend(logs_p);
        out.extras.extend(extras_p);
    }

    Some(out)
}

/// The `t`-adic valuation of a dense polynomial: the count of vanishing
/// low-order coefficients (`None` for the zero polynomial).
pub(crate) fn t_adic_valuation(p: &KPoly) -> Option<usize> {
    if p.is_zero() {
        return None;
    }
    Some(p.coeffs.iter().position(|c| !c.is_zero()).unwrap_or(0))
}

/// The monomial `t^k` as a polynomial over the coefficient field.
pub(crate) fn kpoly_monomial(top: usize, n: usize, k: usize) -> KPoly {
    let mut coeffs = vec![KElem::zero(n); k];
    coeffs.push(KElem::one(n));
    KPoly {
        top,
        coeffs,
        n_vars: n,
    }
}

/// `p / t^k` when `t^k` divides `p`: drop the vanishing low-order
/// coefficients (the result is already trimmed since `p` was).
fn shift_down(p: &KPoly, k: usize) -> KPoly {
    debug_assert!(k <= p.coeffs.len());
    debug_assert!(p.coeffs[..k].iter().all(|c| c.is_zero()));
    KPoly {
        top: p.top,
        coeffs: p.coeffs[k..].to_vec(),
        n_vars: p.n_vars,
    }
}

/// Hermite reduction in `k(t)`: returns `(g, a1, d1)` with
/// `a/d = D g + a1/d1` and `d1` squarefree.
///
/// Returns `None` when a denominator factor violates the reduction's
/// precondition `gcd(u·Dv, v) = 1`. That can only happen for *special*
/// polynomials; at a hyperexponential level the special factors are the
/// powers of `t`, which [`integrate_level`] splits off before calling, so a
/// failure here means an unmodelled dependency slipped past the tower merge
/// — declining is honest, proceeding would silently mis-scale the answer
/// (the 0.28.0 latent bug's mechanism, previously hidden behind a
/// `debug_assert`).
pub(crate) fn hermite_tower(
    tower: &Tower,
    level: usize,
    a: &KPoly,
    d: &KPoly,
) -> Option<(KRat, KPoly, KPoly)> {
    let top = level;
    let n = d.n_vars;
    let factors = d.square_free();
    let m = factors.iter().map(|&(_, k)| k).max().unwrap_or(1);
    if m <= 1 {
        return Some((
            KRat::new(KPoly::zero(top, n), KPoly::one(top, n)),
            a.clone(),
            d.clone(),
        ));
    }
    let (v, _) = factors
        .iter()
        .rev()
        .find(|&(_, k)| *k == m)
        .expect("max multiplicity factor exists");
    let vm = kpoly_pow(v, m as u64);
    let (u, rem) = d.div_rem(&vm);
    debug_assert!(rem.is_zero());

    // Solve s·B + t·v = a with B = u·D v (full tower derivative).
    let tgen = &tower.gens[level - 1];
    let dv = tower_diff_kpoly(v, &tower.gens[..level - 1], &tgen.dt);
    let b = u.mul(&dv);
    let (g0, s0, _t0) = b.eea(v);
    if !g0.is_one() {
        return None;
    }
    let (_, s) = a.mul(&s0).div_rem(v);
    let (t, rem) = a.sub(&s.mul(&b)).div_rem(v);
    debug_assert!(rem.is_zero());

    // ∫ a/(u v^m) = -(s/(m-1))/v^(m-1) + ∫ (t + u·D(s/(m-1))) / (u v^(m-1))
    let inv_m1 = Rational::new(1, (m - 1) as i64);
    let s_scaled = s.mul_kelem(&KElem::from_rational(&inv_m1, n));
    let vm1 = kpoly_pow(v, (m - 1) as u64);
    let ds = tower_diff_kpoly(&s_scaled, &tower.gens[..level - 1], &tgen.dt);
    let new_a = t.add(&u.mul(&ds));
    let new_d = u.mul(&vm1);
    let (g2, a1, d1) = hermite_tower(tower, level, &new_a, &new_d)?;

    // Combine: g = -s_scaled/vm1 + g2.
    let g_term = KRat::new(s_scaled.neg(), vm1);
    Some((krat_add(&g_term, &g2), a1, d1))
}

fn kpoly_pow(p: &KPoly, mut k: u64) -> KPoly {
    let mut result = KPoly::one(p.top, p.n_vars);
    let mut base = p.clone();
    while k > 0 {
        if k & 1 == 1 {
            result = result.mul(&base);
        }
        base = base.mul(&base);
        k >>= 1;
    }
    result
}

fn krat_add(a: &KRat, b: &KRat) -> KRat {
    if a.num.is_zero() {
        return b.clone();
    }
    if b.num.is_zero() {
        return a.clone();
    }
    KRat::new(a.num.mul(&b.den).add(&b.num.mul(&a.den)), a.den.mul(&b.den))
}

/// Whether `p == c·q` for a rational constant `c`; returns the constant.
fn kpoly_scalar_multiple(p: &KPoly, q: &KPoly) -> Option<Rational> {
    if q.is_zero() {
        return None;
    }
    if p.degree() != q.degree() {
        return None;
    }
    let dom = RationalDomain;
    // Candidate from leading coefficients: c = lc(p)/lc(q) must be rational.
    let c = kelem_to_rational(&p.lc().div(&q.lc())?)?;
    // Verify across all coefficients.
    for i in 0..=p.degree()? {
        let pc = p.coeff_at(i);
        let qc = q.coeff_at(i).mul_rational(&c);
        if !pc.eq_cross(&qc) {
            return None;
        }
    }
    let _ = dom;
    Some(c)
}

/// Extract a rational number from a constant field element.
fn kelem_to_rational(e: &KElem) -> Option<Rational> {
    e.as_rational()
}

// ------------------------------------------------------------------
//  Polynomial part
// ------------------------------------------------------------------

/// The tuple returned by the polynomial-part integrators:
/// `(poly, constant, logs, extras)`.
type KPolyIntegral<'a> = (KPoly, KElem, Vec<(Rational, KElem)>, Vec<Atom<'a>>);

/// Integrate `p ∈ k[t]` at the given level. Returns
/// `(poly, constant, logs, extras)`.
fn integrate_kpoly<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    p: &KPoly,
) -> Option<KPolyIntegral<'a>> {
    let top = level;
    let n = tower.n_vars();
    if p.is_zero() {
        return Some((KPoly::zero(top, n), KElem::zero(n), Vec::new(), Vec::new()));
    }
    if level == 0 {
        // Base ℚ(x): termwise integration.
        let dp = kpoly_to_dpoly(p)?;
        let int = super::rational::poly_integrate(&dp);
        return Some((
            dpoly_to_kpoly(&int, top, n),
            KElem::zero(n),
            Vec::new(),
            Vec::new(),
        ));
    }

    let tgen = &tower.gens[level - 1];
    let m = p.degree()?;
    match tgen.kind {
        GenKind::Constant | GenKind::Log => integrate_kpoly_primitive(ctx, tower, level, p, m),
        GenKind::Exp => integrate_kpoly_hyperexp(ctx, tower, level, p, m),
    }
}

/// Primitive level (`t = log(u)`): undetermined coefficients with the
/// top constant determined by a logarithm constraint.
fn integrate_kpoly_primitive<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    p: &KPoly,
    m: usize,
) -> Option<KPolyIntegral<'a>> {
    let top = level;
    let n = tower.n_vars();
    let tgen = &tower.gens[level - 1];
    let dt = &tgen.dt;
    // The argument `u` as a field element (for recognizing log(u) = t).
    let prefix = tower.gen_atoms()[..level].to_vec();
    let u_rf = atom_to_rational_extended(tgen.arg, &prefix, n)?;
    let u_k = KElem::new(u_rf.numerator, u_rf.denominator);

    let mut logs: Vec<(Rational, KElem)> = Vec::new();
    let mut extras: Vec<Atom> = Vec::new();
    let mut q = vec![KElem::zero(n); m + 2];

    // Layer m: free constant c = q_{m+1}, fixed by the log(u) constraint.
    let int_pm = integrate_kelem_or_fallback(ctx, tower, level - 1, p.coeff_at(m))?;
    // ∫Dt must be exactly log(u) = t.
    let int_dt = integrate_kelem_or_fallback(ctx, tower, level - 1, dt.clone())?;
    let (w_elem, w_logs, w_extras) = split_result(&int_dt);
    if !w_elem.is_zero()
        || !w_extras.is_empty()
        || w_logs.len() != 1
        || !RationalDomain.is_one(&w_logs[0].0)
        || !w_logs[0].1.eq_cross(&u_k)
    {
        return None;
    }
    let (a_elem, a_logs, a_extras) = split_result(&int_pm);
    // `q_m = ∫(p_m − (m+1)·q_{m+1}·Dt)` must be an element of the
    // coefficient field `k`. A logarithm `log(v)` with `v ≠ u` or an
    // unevaluated special term (`atan`, `Integral`, …) is *not* in `k`, so
    // the layer-`m` coefficient cannot be formed: the recurrence has no
    // solution in `k[t]` and we must decline.
    //
    // Only for `m = 0` may those terms be passed through: there the
    // coefficient *is* the whole answer term (`t⁰ = 1`), so accumulating
    // them at the top level is sound.
    if m >= 1 && !a_extras.is_empty() {
        return None;
    }
    let mut c = RationalDomain.zero();
    let mut found = false;
    for (d, v) in a_logs {
        if v.eq_cross(&u_k) {
            if found {
                return None;
            }
            found = true;
            c = RationalDomain.div(&d, &Rational::new((m + 1) as i64, 1))?;
        } else if m >= 1 {
            return None;
        } else {
            logs.push((d, v));
        }
    }
    extras.extend(a_extras);
    q[m + 1] = KElem::from_rational(&c, n);
    q[m] = a_elem;

    // Lower layers: no freedom left; a log(u) term here is fatal.
    for i in (0..m).rev() {
        let h = p.coeff_at(i).sub(
            &q[i + 1]
                .mul_rational(&Rational::new((i + 1) as i64, 1))
                .mul(dt),
        );
        let int_h = integrate_kelem_or_fallback(ctx, tower, level - 1, h)?;
        let (e, ls, es) = split_result(&int_h);
        // Same field membership requirement as above, for every layer but
        // the `t⁰` one: `q_i = ∫h` must lie in `k`, so a surviving
        // logarithm or special term makes the recurrence unsolvable.
        if i >= 1 && (!ls.is_empty() || !es.is_empty()) {
            return None;
        }
        q[i] = e;
        for (d, v) in ls {
            if v.eq_cross(&u_k) {
                return None;
            }
            logs.push((d, v));
        }
        extras.extend(es);
    }

    Some((
        KPoly {
            top,
            coeffs: q,
            n_vars: n,
        },
        KElem::zero(n),
        logs,
        extras,
    ))
}

/// Hyperexponential level (`t = exp(u)`): each exponential layer is an
/// independent Risch differential equation.
fn integrate_kpoly_hyperexp<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    p: &KPoly,
    m: usize,
) -> Option<KPolyIntegral<'a>> {
    let top = level;
    let n = tower.n_vars();
    let tgen = &tower.gens[level - 1];
    let t = KElem::var(top, n);
    let du = tgen.dt.div(&t)?;

    let mut q = vec![KElem::zero(n); m + 1];
    for i in (1..=m).rev() {
        let f_i = du.mul_rational(&Rational::new(i as i64, 1));
        q[i] = rde_solve(tower, level - 1, &f_i, &p.coeff_at(i))?;
    }

    // The t⁰ coefficient is integrated recursively in the field below.
    let mut constant = KElem::zero(n);
    let mut logs = Vec::new();
    let mut extras = Vec::new();
    let p0 = p.coeff_at(0);
    if !p0.is_zero() {
        let int_0 = integrate_kelem_or_fallback(ctx, tower, level - 1, p0)?;
        let (e, ls, es) = split_result(&int_0);
        constant = e;
        logs = ls;
        extras = es;
    }

    Some((
        KPoly {
            top,
            coeffs: q,
            n_vars: n,
        },
        constant,
        logs,
        extras,
    ))
}

// ------------------------------------------------------------------
//  Base level (ℚ(x)) adapter
// ------------------------------------------------------------------

/// Integrate a single coefficient-field element `h` at `level`, falling
/// back to an unevaluated `Integral` atom when the Risch recursion cannot
/// handle it (e.g. `∫(1/x) = log(x)` needs a fresh generator the current
/// tower does not contain).
fn integrate_kelem_or_fallback<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    h: KElem,
) -> Option<LevelResult<'a>> {
    let top = level;
    let n = tower.n_vars();
    if let Some(res) = integrate_level(ctx, tower, level, &krat_const(h.clone(), top, n)) {
        return Some(res);
    }
    let mut out = LevelResult::empty(top, n);
    let gens = tower.gen_atoms();
    let h_atom = kelem_to_atom(ctx, &h, &gens)?;
    out.extras.push(ctx.fun("Integral", &[h_atom, tower.x]));
    Some(out)
}

fn integrate_level_base<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    f: &KRat,
) -> Option<LevelResult<'a>> {
    let n = tower.n_vars();
    // At the base level the element may not be a polynomial in x with
    // constant coefficients (e.g. 1/x) — collapse the whole rational
    // function first, then view num/den as dense univariate polynomials.
    let elem = f.kelem();
    let num = sparse_to_dpoly(&elem.num)?;
    let den = sparse_to_dpoly(&elem.den)?;
    let s = integrate_dpoly_structured(ctx, &num, &den, tower.x)?;
    let mut out = LevelResult::empty(0, n);
    out.poly = dpoly_to_kpoly(&s.poly, 0, n);
    if let Some((g_num, g_den)) = &s.rational {
        out.rational = KRat::new(dpoly_to_kpoly(g_num, 0, n), dpoly_to_kpoly(g_den, 0, n));
    }
    for (c, v) in &s.log_terms {
        out.logs.push((c.clone(), dpoly_to_kelem(v, n)));
    }
    out.extras = s.extra_atoms;
    Some(out)
}

// ------------------------------------------------------------------
//  Conversions
// ------------------------------------------------------------------

fn krat_const(e: KElem, top: usize, n: usize) -> KRat {
    KRat::new(KPoly::from_kelem(e, top), KPoly::one(top, n))
}

/// Split a level result into `(elem, logs, extras)`.
fn split_result<'a>(res: &LevelResult<'a>) -> (KElem, Vec<(Rational, KElem)>, Vec<Atom<'a>>) {
    (res.elem_part(), res.logs.clone(), res.extras.clone())
}

/// Dense univariate view of a `KPoly` whose coefficients are constants.
fn kpoly_to_dpoly(p: &KPoly) -> Option<DPoly> {
    let mut coeffs = Vec::with_capacity(p.coeffs.len());
    for c in &p.coeffs {
        coeffs.push(kelem_to_rational(c)?);
    }
    Some(DPoly::from_coeffs(RationalDomain, coeffs))
}

/// Dense univariate view of a sparse polynomial that uses only variable 0.
fn sparse_to_dpoly(p: &crate::tower::elem::Sparse) -> Option<DPoly> {
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

fn dpoly_to_kelem(p: &DPoly, n: usize) -> KElem {
    let terms = p
        .coeffs()
        .iter()
        .enumerate()
        .filter(|&(_, c)| !RationalDomain.is_zero(c))
        .map(|(i, c)| (vec_i(i, n), c.clone()))
        .collect();
    KElem::from_poly(crate::tower::elem::Sparse::from_terms(
        RationalDomain,
        n,
        terms,
    ))
}

fn dpoly_to_kpoly(p: &DPoly, top: usize, n: usize) -> KPoly {
    let coeffs = p
        .coeffs()
        .iter()
        .map(|c| KElem::from_rational(c, n))
        .collect();
    let mut r = KPoly {
        top,
        coeffs,
        n_vars: n,
    };
    // Trim trailing zeros (KPoly invariant).
    while let Some(last) = r.coeffs.last() {
        if last.is_zero() {
            r.coeffs.pop();
        } else {
            break;
        }
    }
    r
}

fn vec_i(i: usize, n: usize) -> Vec<usize> {
    let mut v = vec![0; n];
    v[0] = i;
    v
}

fn integral_fallback_atom<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    a: &KPoly,
    d: &KPoly,
) -> Option<Atom<'a>> {
    let gens = tower.gen_atoms();
    let num = kelem_to_atom(ctx, &a.kelem(), &gens)?;
    let den = kelem_to_atom(ctx, &d.kelem(), &gens)?;
    let expr = ctx.mul(&[num, ctx.pow(den, ctx.num(-1))]);
    let _ = level;
    Some(ctx.fun("Integral", &[expr, tower.x]))
}

fn kelem_to_atom<'a>(ctx: &'a AtomArena<'a>, e: &KElem, gens: &[Atom<'a>]) -> Option<Atom<'a>> {
    rational_to_atom(
        ctx,
        &GeneratorField::from_num_den(e.num.clone(), e.den.clone()),
        gens,
    )
}

/// Assemble a level result into an atom over the generator atoms.
fn level_result_to_atom<'a>(
    ctx: &'a AtomArena<'a>,
    tower: &Tower<'a>,
    level: usize,
    res: &LevelResult<'a>,
) -> Option<Atom<'a>> {
    let gens = tower.gen_atoms();
    let t_atom = gens[level];
    let n = tower.n_vars();
    let mut parts: Vec<Atom> = Vec::new();

    if !res.rational.num.is_zero() {
        parts.push(kelem_to_atom(ctx, &res.rational.kelem(), &gens)?);
    }
    // Polynomial part: Σ cᵢ·tⁱ.
    let mut poly_terms: Vec<Atom> = Vec::new();
    for (i, c) in res.poly.coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let c_atom = kelem_to_atom(ctx, c, &gens)?;
        let tp = match i {
            0 => c_atom,
            1 => ctx.mul(&[c_atom, t_atom]),
            _ => ctx.mul(&[c_atom, ctx.pow(t_atom, ctx.num(i as i64))]),
        };
        poly_terms.push(tp);
    }
    match poly_terms.len() {
        0 => {}
        1 => parts.push(poly_terms[0]),
        _ => parts.push(ctx.add(&poly_terms)),
    }
    if !res.constant.is_zero() {
        parts.push(kelem_to_atom(ctx, &res.constant, &gens)?);
    }
    for (c, v) in &res.logs {
        let v_atom = kelem_to_atom(ctx, v, &gens)?;
        let log = ctx.fun("log", &[v_atom]);
        parts.push(if RationalDomain.is_one(c) {
            log
        } else {
            ctx.mul(&[rational_const_to_atom(ctx, c)?, log])
        });
    }
    parts.extend(res.extras.iter().copied());

    let _ = n;
    Some(match parts.len() {
        0 => ctx.num(0),
        1 => parts[0],
        _ => ctx.add(&parts),
    })
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;
    use crate::derivative::diff;
    use crate::tower::convert::atom_to_rational;

    /// Exact check: d/dx(result) == integrand as rational functions over
    /// the tower generators.
    fn assert_risch_antiderivative<'a>(ctx: &'a AtomArena<'a>, integrand: Atom<'a>, var: Symbol) {
        let result = risch_integrate(ctx, integrand, var).expect("risch result");
        let tower = build_tower(ctx, integrand, var).expect("tower");
        let gens = tower.gen_atoms();
        let d = normalize(ctx, diff(ctx, result, var));
        let lhs = atom_to_rational(d, &gens).expect("derivative converts");
        // The tower may rewrite the integrand over merged generators
        // (0.28.0), so the comparison must use the rewritten form.
        let rhs = atom_to_rational(normalize(ctx, tower.expr), &gens).expect("integrand converts");
        // Compare semantically: cross-multiply since the two sides may be
        // different representatives of the same rational function.
        let lhs_e = KElem::new(lhs.numerator, lhs.denominator);
        let rhs_e = KElem::new(rhs.numerator, rhs.denominator);
        assert!(lhs_e.eq_cross(&rhs_e), "d/dx(result) != integrand");
    }

    /// The 0.28.0 latent bug (caught by the tower check there, fixed in
    /// 0.29.0): the exp-level rational part mis-scaled negative powers of
    /// `t` because Hermite reduction was applied to the special factor
    /// `t^k` (`gcd(t, Dt) = t ≠ 1`). The Laurent split now routes those
    /// layers through per-layer RDEs.
    #[test]
    fn integrate_exp_negative_power_layers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // Minimal repro from the 0.28.0 post-mortem:
        // ∫ (exp(x)²+1)³/(8·exp(x)⁴) dx
        //   = t²/16 + 3x/8 − 3/(16·t²) − 1/(32·t⁴).
        let t = ctx.fun("exp", &[x]);
        let t2 = ctx.pow(t, ctx.num(2));
        let num = ctx.pow(ctx.add(&[t2, ctx.num(1)]), ctx.num(3));
        let den = ctx.mul(&[ctx.num(8), ctx.pow(t, ctx.num(4))]);
        let f = ctx.mul(&[num, ctx.pow(den, ctx.num(-1))]);
        assert_risch_antiderivative(&ctx, f, Symbol::new("x"));
        // ∫ (1 + exp(x))/exp(x) dx = x − exp(−x): a single negative layer.
        let g = ctx.mul(&[ctx.add(&[ctx.num(1), t]), ctx.pow(t, ctx.num(-1))]);
        assert_risch_antiderivative(&ctx, g, Symbol::new("x"));
    }

    /// rubi-01524's shape: dependent `exp(4x³) = exp(x³)⁴` merges (0.29.0),
    /// leaving a plain polynomial in `t` at the hyperexponential level.
    #[test]
    fn integrate_exp_integer_multiple_layers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x3 = ctx.pow(ctx.var("x"), ctx.num(3));
        let t = ctx.fun("exp", &[x3]);
        let t4 = ctx.fun("exp", &[ctx.mul(&[ctx.num(4), x3])]);
        let inner = ctx.pow(
            ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), t4])]),
            ctx.num(2),
        );
        let f = ctx.mul(&[t, inner, ctx.pow(ctx.var("x"), ctx.num(2))]);
        assert_risch_antiderivative(&ctx, f, Symbol::new("x"));
    }

    #[test]
    fn integrate_log_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ log(x) dx = x·log(x) - x
        assert_risch_antiderivative(&ctx, ctx.fun("log", &[x]), Symbol::new("x"));
    }

    #[test]
    fn integrate_log_x_over_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ log(x)/x dx = log(x)²/2
        let expr = ctx.mul(&[ctx.fun("log", &[x]), ctx.pow(x, ctx.num(-1))]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_x_log_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x·log(x) dx
        let expr = ctx.mul(&[x, ctx.fun("log", &[x])]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_exp_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert_risch_antiderivative(&ctx, ctx.fun("exp", &[x]), Symbol::new("x"));
    }

    #[test]
    fn integrate_x_exp_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x·exp(x) dx = (x-1)·exp(x)
        let expr = ctx.mul(&[x, ctx.fun("exp", &[x])]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_x2_exp_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x²·exp(x) dx = (x² - 2x + 2)·exp(x)
        let expr = ctx.mul(&[ctx.pow(x, ctx.num(2)), ctx.fun("exp", &[x])]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_exp_x_over_x_fails() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ exp(x)/x dx is not elementary (Ei) → None.
        let expr = ctx.mul(&[ctx.fun("exp", &[x]), ctx.pow(x, ctx.num(-1))]);
        assert!(risch_integrate(&ctx, expr, Symbol::new("x")).is_none());
    }

    #[test]
    fn integrate_log_x_plus_1() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ log(x+1) dx = (x+1)·log(x+1) - x. The primitive-layer
        // undetermined-constant choice (q_{m+1}) must make lower layers
        // integrable; our current fragment only fixes it via a log(u)
        // term in ∫p_m, which misses this case. Documented limitation.
        let expr = ctx.fun("log", &[ctx.add(&[x, ctx.num(1)])]);
        let _ = risch_integrate(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_exp_x2_times_x() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x·exp(x²) dx = exp(x²)/2
        let expr = ctx.mul(&[x, ctx.fun("exp", &[ctx.pow(x, ctx.num(2))])]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_one_over_x_log_x_sq() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x·log(x)²) dx = -1/(x·log(x))·x … = -1/log(x) + C…
        // Precisely: Hermite in k(t) yields -1/log(x) + ∫(-1/(x²·log(x)))
        // whose remainder is non-elementary here → expect None overall or
        // a partial answer; require at least no panic and None handling.
        let expr = ctx.pow(
            ctx.mul(&[x, ctx.pow(ctx.fun("log", &[x]), ctx.num(2))]),
            ctx.num(-1),
        );
        let _ = risch_integrate(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_log_squared() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ log(x)² dx = x·log(x)² - 2x·log(x) + 2x
        let expr = ctx.pow(ctx.fun("log", &[x]), ctx.num(2));
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    #[test]
    fn integrate_exp_plus_poly() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ (x² + exp(x)) dx
        let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.fun("exp", &[x])]);
        assert_risch_antiderivative(&ctx, expr, Symbol::new("x"));
    }

    /// Regression for the 0.27.1 wrong answers `rubi-01482`, `rubi-01662`
    /// and `rubi-00531`.
    ///
    /// The primitive-level recurrence requires every coefficient `q_i` to be
    /// an element of the coefficient field `k`: `q_i = ∫h_i` must therefore
    /// contribute no logarithm and no "extra" (a special function or an
    /// unevaluated `Integral`). The old code split those terms off, kept only
    /// the element part, and re-attached them at the top level — sound for the
    /// `t⁰` layer (where the coefficient *is* the answer term), unsound for
    /// every layer above it, because there the term would have to be
    /// multiplied by `t^i`.
    ///
    /// The two shapes below are the exact residuals the `heuristic` parts step
    /// fed back into risch for those corpus cases; both used to come back with
    /// a bogus closed form (`log(1+x²)/(1+x²) → atan(x)`), dropping the `log`
    /// factor. The true antiderivatives are non-elementary (dilogarithm), so
    /// the only correct answer is `None`.
    #[test]
    fn primitive_level_declines_when_a_coefficient_integral_leaves_the_field() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for input in [
            // Minimal reproducer: `∫dx/(1+x²) = atan(x)` is an extra, not a
            // field element, so the layer-1 coefficient cannot be formed.
            "log(1 + x^2)*(1 + x^2)^-1",
            // The same residual with the polynomial part that `∫ x^m/(1+x²)^n`
            // contributes; risch used to return `x - 2*atan(x) - …`.
            "(x^2 + (-1)*log(1 + x^2) + (-1*(2^-1))*((1 + x^2)^-1))*((1 + x^2)^-1)",
        ] {
            let expr = ocas_parse::parse(&ctx, input).expect("parse");
            assert!(
                risch_integrate(&ctx, expr, Symbol::new("x")).is_none(),
                "risch must decline `{input}`: its primitive-level coefficient \
                 integral is not an element of the coefficient field"
            );
        }
    }

    /// The pipeline must fall back to an honest `Integral(...)` residue for
    /// those shapes rather than emit a closed form risch cannot justify.
    #[test]
    fn out_of_field_residuals_fall_back_to_an_unevaluated_integral() {
        for input in [
            "log(1 + x^2)*(1 + x^2)^-1",
            "x^5*atan(x)/(1 + x^2)^2",
            "x^3*atan(x)/(1 + x^2)^2",
            "atanh(1 + x)/(2 + 2*x)",
        ] {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = ocas_parse::parse(&ctx, input).expect("parse");
            let r = crate::integrate(&ctx, expr, Symbol::new("x")).to_string();
            assert!(
                r.contains("Integral("),
                "`{input}` must be left unevaluated, got a closed form: {r}"
            );
        }
    }

    /// The restriction must not cost the primitive level its working cases:
    /// logarithms produced by the `t⁰` layer are still accumulated, and the
    /// polynomial recurrence still solves exact shapes.
    #[test]
    fn primitive_level_still_solves_its_fragment() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x·log(x) dx and ∫ log(x)² dx exercise the layer-1/layer-2
        // recurrence end to end.
        let xl = ctx.mul(&[x, ctx.fun("log", &[x])]);
        assert_risch_antiderivative(&ctx, xl, Symbol::new("x"));
        let l2 = ctx.pow(ctx.fun("log", &[x]), ctx.num(2));
        assert_risch_antiderivative(&ctx, l2, Symbol::new("x"));
        // ∫ x²·log(x)² dx: the `t⁰` layer's own `log(x)` terms stay additive.
        let x2l2 = ctx.mul(&[ctx.pow(x, ctx.num(2)), l2]);
        assert_risch_antiderivative(&ctx, x2l2, Symbol::new("x"));
    }
}
