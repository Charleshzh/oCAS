//! Algebraic Rothstein–Trager roots at the base level ℚ(x) (0.30.0, Wave 1).
//!
//! The 0.29.0 logarithmic part accepts only rational roots of the
//! Rothstein–Trager resultant `R(z)`. This module extends the base level
//! ([`super::rational`]) to **algebraic constant roots** in a bounded
//! fragment (see `docs/planning/ALGEBRAIC_RT_SPEC_CN.md`):
//!
//! - irreducible **quadratic** factors of `R(z)`: the two roots live in one
//!   quadratic field `ℚ(√Δ)`; several quadratic factors are accepted as
//!   long as they share one squarefree discriminant signature;
//! - one irreducible **biquadratic quartic** factor `A·z⁴ + B·z² + C` with
//!   `B² − 4AC < 0` (the complex 4-cycle: the `1/(1+x⁴)` family).
//!
//! For every accepted root `c` the factor `v = gcd(a − c·d′, d)` is computed
//! over the algebraic number field (generic Euclidean gcd; the same
//! `DenseUnivariatePolynomial<AlgebraicNumberField>` machinery the 0.17.0
//! Trager factorization uses), and the answer is verified **inside the
//! extension field** before emission (`Πv == d` plus the cross-multiplied
//! residue identity). Complex roots are emitted as real `log(u²+w²)` /
//! `atan(w/u)` conjugate-pair forms; the formal identity behind the
//! realification is
//!
//! ```text
//! c·log(u + i·w) + c̄·log(u − i·w) = p·log(u² + w²) − 2q·atan(w/u),
//! c = p + q·i
//! ```
//!
//! (the `atan(w/u)` branch differs from other real forms by a per-interval
//! constant; the derivative is exact — the certificate layer is the final
//! gate). Factors outside the fragment (degree 3, non-biquadratic degree 4,
//! degree ≥ 5, mixed extension fields, real-rooted biquadratics) decline
//! honestly, exactly like the 0.29.0 rational-root path.

use ocas_atom::{Atom, AtomArena};
use ocas_domain::number_theory::factor::factor_integer;
use ocas_domain::{
    AlgebraicElement, AlgebraicNumberField, Domain, Integer, IntegerDomain, Rational,
    RationalDomain,
};
use ocas_poly::DenseUnivariatePolynomial;

use super::rational::{
    DLogOrAtom, Sqrt, charge_rational_work, poly_atom, scale_atom, sqrt_positive_rational,
};
use crate::tower::convert::rational_const_to_atom;

type DPoly = DenseUnivariatePolynomial<RationalDomain>;
type ZPoly = DenseUnivariatePolynomial<IntegerDomain>;
type APoly = DenseUnivariatePolynomial<AlgebraicNumberField>;
type AElem = AlgebraicElement<Rational>;

/// Integrate the logarithmic part `a/d` when the RT resultant `rt` has
/// non-rational constant roots inside the supported fragment.
///
/// `linear_roots` are the rational roots already found by the caller; their
/// log terms are included in the answer and in the in-field verification.
/// Returns `None` outside the fragment or when any internal check fails —
/// the caller then emits the usual unevaluated `Integral` term.
pub(crate) fn rothstein_trager_algebraic<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    rt: &DPoly,
    linear_roots: &[Rational],
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let zpoly = integer_primitive(rt)?;
    if charge_rational_work(1 + zpoly.degree().unwrap_or(0) as u64 * 64) {
        return None;
    }
    let factors = zpoly.factor();

    // Classify the distinct irreducible factors.
    let mut quadratics: Vec<QuadFactor> = Vec::new();
    let mut quartics: Vec<QuarticFactor> = Vec::new();
    for (fac, _mult) in &factors {
        match fac.degree()? {
            0 | 1 => {} // rational roots are the caller's input
            2 => quadratics.push(QuadFactor::from_zpoly(fac)),
            4 => quartics.push(QuarticFactor::from_zpoly(fac)?),
            _ => return None,
        }
    }
    if quadratics.is_empty() && quartics.is_empty() {
        return None;
    }
    if !quartics.is_empty() && (!quadratics.is_empty() || quartics.len() > 1) {
        // Mixed extension fields are outside the 0.30.0 fragment.
        return None;
    }

    if let Some(q4) = quartics.first() {
        quartic_log_part(ctx, a, d, linear_roots, q4, x)
    } else {
        quadratic_log_part(ctx, a, d, linear_roots, &quadratics, x)
    }
}

// ---------------------------------------------------------------------
//  Factor classification
// ---------------------------------------------------------------------

/// `A·z² + B·z + C`, integer coefficients, irreducible over ℚ.
struct QuadFactor {
    a: Integer,
    b: Integer,
    c: Integer,
}

impl QuadFactor {
    fn from_zpoly(f: &ZPoly) -> Self {
        let c = |i: usize| {
            f.coeffs()
                .get(i)
                .cloned()
                .unwrap_or_else(|| Integer::from(0))
        };
        Self {
            a: c(2),
            b: c(1),
            c: c(0),
        }
    }

    /// Discriminant `B² − 4AC` (nonzero: the factor is irreducible).
    fn discriminant(&self) -> Integer {
        let b2 = &self.b * &self.b;
        let four_ac = Integer::from(4) * &self.a * &self.c;
        b2 - four_ac
    }
}

/// `A·z⁴ + B·z² + C`, integer coefficients, irreducible over ℚ.
struct QuarticFactor {
    a: Integer,
    b: Integer,
    c: Integer,
}

impl QuarticFactor {
    /// Biquadratic only: the `z³` and `z` coefficients must vanish.
    fn from_zpoly(f: &ZPoly) -> Option<Self> {
        let c = |i: usize| {
            f.coeffs()
                .get(i)
                .cloned()
                .unwrap_or_else(|| Integer::from(0))
        };
        if !c(3).is_zero() || !c(1).is_zero() {
            return None;
        }
        Some(Self {
            a: c(4),
            b: c(2),
            c: c(0),
        })
    }
}

/// Squarefree signature of a nonzero integer: `(is_negative, squarefree
/// part, square root of the square part)` with `n = sign · sqf · k²`.
/// Also used by the tower-level algebraic roots path (`logpart`).
pub(crate) fn squarefree_signature(n: &Integer) -> Option<(bool, Integer, Integer)> {
    let neg = n.is_negative();
    let abs = n.abs();
    let mut sqf = Integer::from(1);
    let mut k = Integer::from(1);
    for (p, e) in factor_integer(&abs) {
        if e % 2 == 1 {
            sqf *= &p;
        }
        for _ in 0..(e / 2) {
            k *= &p;
        }
    }
    // Reassemble check: |n| == sqf · k².
    if &sqf * &k * &k != abs {
        return None;
    }
    Some((neg, sqf, k))
}

// ---------------------------------------------------------------------
//  Field helpers
// ---------------------------------------------------------------------

/// Embed a ℚ-polynomial into `K[x]`.
fn embed_poly(field: &AlgebraicNumberField, p: &DPoly) -> APoly {
    APoly::from_coeffs(
        field.clone(),
        p.coeffs()
            .iter()
            .map(|c| field.from_base(c.clone()))
            .collect(),
    )
}

/// `a − c·d′` in `K[x]`.
fn sub_scaled(field: &AlgebraicNumberField, a: &DPoly, c: &AElem, dp: &DPoly) -> APoly {
    let deg = a.coeffs().len().max(dp.coeffs().len());
    let zero = RationalDomain.zero();
    let coeffs = (0..deg)
        .map(|i| {
            let ai = a.coeffs().get(i).unwrap_or(&zero);
            let di = dp.coeffs().get(i).unwrap_or(&zero);
            field.sub(
                &field.from_base(ai.clone()),
                &field.mul(c, &field.from_base(di.clone())),
            )
        })
        .collect();
    APoly::from_coeffs(field.clone(), coeffs)
}

/// Monic gcd of `a − c·d′` and `d` over `K`, or `None` when it is trivial
/// (a root of `R` must give a nontrivial gcd).
fn residue_factor(field: &AlgebraicNumberField, a: &DPoly, d: &DPoly, c: &AElem) -> Option<APoly> {
    let dp = d.derivative();
    let shifted = sub_scaled(field, a, c, &dp);
    let d_k = embed_poly(field, d);
    let (g, _, _) = shifted.extended_gcd_poly(&d_k);
    if g.degree().unwrap_or(0) == 0 {
        return None;
    }
    Some(g)
}

/// In-field verification (spec §2.6): `Πv == d` and the cross-multiplied
/// residue identity `a·Πv == d·Σ c·Dv·(Π/v)` over `K[x]`.
fn verify_in_field(
    field: &AlgebraicNumberField,
    a: &DPoly,
    d: &DPoly,
    pairs: &[(AElem, APoly)],
) -> bool {
    let d_k = embed_poly(field, d);
    let a_k = embed_poly(field, a);
    let mut prod = APoly::from_coeffs(field.clone(), vec![field.one()]);
    for (_, v) in pairs {
        prod = prod.mul(v);
    }
    if prod != d_k {
        return false;
    }
    // Residue identity: a·Πv == d·Σ cᵢ·Dvᵢ·(Π/vᵢ). Division is exact by (i).
    let mut acc = APoly::from_coeffs(field.clone(), vec![field.zero()]);
    for (c, v) in pairs {
        let dv = v.derivative();
        let (quot, rem) = prod.div_rem(v).expect("product divisible by its factor");
        debug_assert!(rem.is_zero());
        acc = acc.add(&dv.mul(&quot).mul_scalar(c));
    }
    let lhs = a_k.mul(&prod);
    let rhs = d_k.mul(&acc);
    lhs == rhs
}

/// Evaluate a coefficient (a polynomial in α of degree < deg m) at a root
/// element ρ ∈ K.
fn eval_coeff(field: &AlgebraicNumberField, coeff: &AElem, rho: &AElem) -> AElem {
    let mut acc = field.zero();
    for c in coeff.coeffs().iter().rev() {
        acc = field.add(&field.mul(&acc, rho), &field.from_base(c.clone()));
    }
    acc
}

/// Substitute `α := ρ` into every coefficient of `v` (a polynomial in x
/// whose coefficients are polynomials in α).
fn substitute_root(field: &AlgebraicNumberField, v: &APoly, rho: &AElem) -> APoly {
    APoly::from_coeffs(
        field.clone(),
        v.coeffs()
            .iter()
            .map(|c| eval_coeff(field, c, rho))
            .collect(),
    )
}

// ---------------------------------------------------------------------
//  Quadratic factors
// ---------------------------------------------------------------------

fn quadratic_log_part<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    linear_roots: &[Rational],
    quads: &[QuadFactor],
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let dom = RationalDomain;
    // All factors must share one squarefree discriminant signature.
    let first = quads[0].discriminant();
    let (neg, sqf, _k0) = squarefree_signature(&first)?;
    let sqf_rat = int_to_rational(&sqf)?;
    let mut per_factor = Vec::with_capacity(quads.len());
    for qf in quads {
        let delta = qf.discriminant();
        let (neg_i, sqf_i, k_i) = squarefree_signature(&delta)?;
        if neg_i != neg || sqf_i != sqf {
            // A second, independent quadratic field: outside the fragment.
            return None;
        }
        per_factor.push((qf, k_i));
    }

    // K = ℚ(α) with α² = ±sqf (sign = discriminant sign).
    let min_poly = if neg {
        // α² + sqf = 0
        vec![sqf_rat.clone(), dom.zero(), dom.one()]
    } else {
        vec![dom.neg(&sqf_rat), dom.zero(), dom.one()]
    };
    let field = AlgebraicNumberField::new(RationalDomain, min_poly);
    let alpha = field.alpha();

    // Linear roots contribute their (rational) v's, embedded into K.
    let mut pairs: Vec<(AElem, APoly)> = Vec::new();
    let mut out: Vec<DLogOrAtom<'a>> = Vec::new();
    let dp = d.derivative();
    for c in linear_roots {
        let v = d.gcd(&a.sub(&dp.mul_scalar(c)));
        if v.degree().unwrap_or(0) == 0 {
            return None;
        }
        pairs.push((field.from_base(c.clone()), embed_poly(&field, &v)));
        out.push(DLogOrAtom::Log(c.clone(), v));
    }

    // √sqf as an atom (Δ and sqf share their square root up to the rational
    // factor k, which the coefficients below carry explicitly).
    let rad_atom = match sqrt_positive_rational(ctx, &sqf_rat)? {
        Sqrt::Rat(_) => return None, // a perfect square would have split over ℚ
        Sqrt::Rad(atom) => atom,
    };

    for (qf, k) in &per_factor {
        let big_a = int_to_rational(&qf.a)?;
        let big_b = int_to_rational(&qf.b)?;
        let two_a = dom.add(&big_a, &big_a);
        let c0 = dom.div(&dom.neg(&big_b), &two_a)?; // −B/(2A)
        let k_rat = int_to_rational(k)?;
        let c1 = dom.div(&k_rat, &two_a)?; // k/(2A): roots c0 ± c1·α

        if charge_rational_work(64) {
            return None;
        }
        if neg {
            // Complex pair: c = c0 ± i·c1·√sqf. One gcd at α = generator,
            // one at its conjugate (both are needed for the in-field check).
            let c_elem = field.element(vec![c0.clone(), c1.clone()]);
            let v = residue_factor(&field, a, d, &c_elem)?;
            let c_conj = field.element(vec![c0.clone(), dom.neg(&c1)]);
            let v_conj = residue_factor(&field, a, d, &c_conj)?;
            pairs.push((c_elem.clone(), v.clone()));
            pairs.push((c_conj, v_conj));
            out.push(complex_pair_term(ctx, &v, &c0, &c1, &sqf_rat, rad_atom, x)?);
        } else {
            // Two real roots: emit c·log(v) per root with √ literals.
            for sign in [1i64, -1] {
                let s = Rational::new(sign, 1);
                let c_elem = field.element(vec![c0.clone(), dom.mul(&s, &c1)]);
                let v = residue_factor(&field, a, d, &c_elem)?;
                pairs.push((c_elem, v.clone()));
                out.push(real_root_term(ctx, &v, &c0, &c1, sign, rad_atom, x)?);
            }
        }
    }
    let _ = alpha;

    if !verify_in_field(&field, a, d, &pairs) {
        return None;
    }
    Some(out)
}

/// `c·log(v)` for a real root `c = c0 + sign·c1·√sqf`; the coefficients of
/// `v` are `e₀ + e₁·α` with `α = √sqf` the (fixed) generator of K — the
/// sign chose *which* root's `v` was computed, the coefficient substitution
/// itself is sign-free.
fn real_root_term<'a>(
    ctx: &'a AtomArena<'a>,
    v: &APoly,
    c0: &Rational,
    c1: &Rational,
    sign: i64,
    rad: Atom<'a>,
    x: Atom<'a>,
) -> Option<DLogOrAtom<'a>> {
    let dom = RationalDomain;
    let s = Rational::new(sign, 1);
    // e0 + e1·α ↦ e0 + e1·√sqf
    let coeff_atom = |e: &AElem| -> Option<Atom<'a>> {
        let e0 = e.coeffs().first().cloned().unwrap_or_else(|| dom.zero());
        let e1 = e.coeffs().get(1).cloned().unwrap_or_else(|| dom.zero());
        let rat_part = rational_const_to_atom(ctx, &e0)?;
        if dom.is_zero(&e1) {
            return Some(rat_part);
        }
        let rad_part = scale_atom(ctx, &e1, rad)?;
        if dom.is_zero(&e0) {
            return Some(rad_part);
        }
        Some(ctx.add(&[rat_part, rad_part]))
    };
    let v_atom = apoly_atom(ctx, v, &coeff_atom, x)?;
    // c = c0 + sign·c1·√sqf
    let c_lin = dom.mul(c1, &s);
    let c_atom = if dom.is_zero(c0) {
        scale_atom(ctx, &c_lin, rad)?
    } else {
        ctx.add(&[
            rational_const_to_atom(ctx, c0)?,
            scale_atom(ctx, &c_lin, rad)?,
        ])
    };
    Some(DLogOrAtom::Atom(
        ctx.mul(&[c_atom, ctx.fun("log", &[v_atom])]),
    ))
}

/// The realified conjugate pair `p·log(u²+w²) − 2q·atan(w/u)` for the root
/// `c = c0 + i·c1·√sqf` (`p = c0`, `q = c1·√sqf`). `v`'s coefficients are
/// `e₀ + e₁·α` with `α = i·√sqf`: `u` collects the `e₀` parts and
/// `w = √sqf · w̃` the `e₁` parts.
fn complex_pair_term<'a>(
    ctx: &'a AtomArena<'a>,
    v: &APoly,
    c0: &Rational,
    c1: &Rational,
    sqf: &Rational,
    rad: Atom<'a>,
    x: Atom<'a>,
) -> Option<DLogOrAtom<'a>> {
    let dom = RationalDomain;
    let deg = v.coeffs().len();
    let mut u_coeffs = Vec::with_capacity(deg);
    let mut w_coeffs = Vec::with_capacity(deg); // w̃: w = √sqf·w̃
    for e in v.coeffs() {
        u_coeffs.push(e.coeffs().first().cloned().unwrap_or_else(|| dom.zero()));
        w_coeffs.push(e.coeffs().get(1).cloned().unwrap_or_else(|| dom.zero()));
    }
    let u = DPoly::from_coeffs(RationalDomain, u_coeffs);
    let w_tilde = DPoly::from_coeffs(RationalDomain, w_coeffs);
    let u_atom = poly_atom(ctx, &u, x)?;
    let w_atom = poly_atom(ctx, &w_tilde, x)?;

    let mut out = Vec::new();
    // p·log(u² + sqf·w̃²) — the log argument is rational (sqf·w̃² has no
    // radical left).
    if !dom.is_zero(c0) {
        let u2 = u.mul(&u);
        let w2 = w_tilde.mul(&w_tilde).mul_scalar(sqf);
        let log_arg = poly_atom(ctx, &u2.add(&w2), x)?;
        out.push(scale_atom(ctx, c0, ctx.fun("log", &[log_arg]))?);
    }
    // +2·c1·√sqf·atan(u/(√sqf·w̃)): the reciprocal argument of the textbook
    // pair formula — `atan(w/u)` and `−atan(u/w)` differ by a per-interval
    // constant, so the derivative is identical, but the `atan(u/w)` shape is
    // what the canonical answers (and the downstream parts machinery) use.
    if !dom.is_zero(c1) {
        let atan_arg = ctx.mul(&[u_atom, ctx.pow(ctx.mul(&[rad, w_atom]), ctx.num(-1))]);
        let atan = ctx.fun("atan", &[atan_arg]);
        let coeff = dom.mul(&Rational::new(2, 1), c1);
        out.push(ctx.mul(&[rational_const_to_atom(ctx, &coeff)?, rad, atan]));
    }
    Some(DLogOrAtom::Atom(match out.len() {
        0 => ctx.num(0),
        1 => out[0],
        _ => ctx.add(&out),
    }))
}

/// Build the atom `Σ coeffᵢ·xⁱ` where each coefficient is produced by the
/// `coeff_atom` callback from its field element.
fn apoly_atom<'a>(
    ctx: &'a AtomArena<'a>,
    v: &APoly,
    coeff_atom: &dyn Fn(&AElem) -> Option<Atom<'a>>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let mut terms = Vec::new();
    for (i, e) in v.coeffs().iter().enumerate() {
        if e.coeffs().is_empty() {
            continue;
        }
        let c = coeff_atom(e)?;
        let term = match i {
            0 => c,
            1 => ctx.mul(&[c, x]),
            _ => ctx.mul(&[c, ctx.pow(x, ctx.num(i as i64))]),
        };
        terms.push(term);
    }
    Some(match terms.len() {
        0 => ctx.num(0),
        1 => terms[0],
        _ => ctx.add(&terms),
    })
}

// ---------------------------------------------------------------------
//  Biquadratic quartic factor (complex 4-cycle)
// ---------------------------------------------------------------------

fn quartic_log_part<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    linear_roots: &[Rational],
    q4: &QuarticFactor,
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let dom = RationalDomain;
    let big_a = int_to_rational(&q4.a)?;
    let big_b = int_to_rational(&q4.b)?;
    let big_c = int_to_rational(&q4.c)?;
    // Monic minimal polynomial z⁴ + (B/A)z² + C/A.
    let a_inv = dom.inv(&big_a)?;
    let b_m = dom.mul(&big_b, &a_inv);
    let c_m = dom.mul(&big_c, &a_inv);
    let field = AlgebraicNumberField::new(
        RationalDomain,
        vec![c_m.clone(), dom.zero(), b_m.clone(), dom.zero(), dom.one()],
    );

    // Monic discriminant Δ₂ = (B/A)² − 4·(C/A); only the complex 4-cycle
    // (Δ₂ < 0) is in the 0.30.0 fragment.
    let delta2 = dom.sub(&dom.mul(&b_m, &b_m), &dom.mul(&Rational::new(4, 1), &c_m));
    if dom.is_zero(&delta2) || !rat_is_neg(&delta2) {
        return None;
    }
    // s± = (sᵣ ± i·sᵢ) with sᵣ = −(B/A)/2, sᵢ = √|Δ₂|/2; |s| = √(C/A) needs
    // C/A > 0 (it equals sᵣ² + sᵢ²).
    if rat_is_neg(&c_m) || dom.is_zero(&c_m) {
        return None;
    }
    let s_re = dom.div(&dom.neg(&b_m), &Rational::new(2, 1))?;

    let alpha = field.alpha();
    if charge_rational_work(256) {
        return None;
    }
    let v_alpha = residue_factor(&field, a, d, &alpha)?;

    // Verification needs all four roots as K-elements. For a biquadratic
    // z⁴ + bz² + c the roots are ±α, ±β with β = γ/α, γ² = c: β ∈ K iff c
    // is a square in ℚ (K is then the splitting field — the cyclotomic-like
    // biquadratics of the acceptance families all qualify). This avoids a
    // full Trager factorization of m over K (measured: it dominated the
    // wave's integration time).
    let gamma = sqrt_rational(&c_m)?;
    let beta = field.div(&field.from_base(gamma), &alpha)?;
    let roots: Vec<AElem> = vec![
        alpha.clone(),
        field.neg(&alpha),
        beta.clone(),
        field.neg(&beta),
    ];
    let mut pairs: Vec<(AElem, APoly)> = Vec::new();
    let dp = d.derivative();
    for c in linear_roots {
        let v = d.gcd(&a.sub(&dp.mul_scalar(c)));
        if v.degree().unwrap_or(0) == 0 {
            return None;
        }
        pairs.push((field.from_base(c.clone()), embed_poly(&field, &v)));
    }
    for rho in &roots {
        pairs.push((rho.clone(), substitute_root(&field, &v_alpha, rho)));
    }
    if !verify_in_field(&field, a, d, &pairs) {
        return None;
    }

    // Emission: roots ±u ± i·v with u = √((|s|+sᵣ)/2), v = √((|s|−sᵣ)/2).
    let abs_s = match sqrt_positive_rational(ctx, &c_m)? {
        Sqrt::Rat(r) => rational_const_to_atom(ctx, &r)?,
        Sqrt::Rad(atom) => atom,
    };
    let abs_d2 = dom.neg(&delta2);
    let s_im = {
        let rad = match sqrt_positive_rational(ctx, &abs_d2)? {
            Sqrt::Rat(r) => rational_const_to_atom(ctx, &r)?,
            Sqrt::Rad(atom) => atom,
        };
        ctx.mul(&[rad, ctx.pow(ctx.num(2), ctx.num(-1))])
    };
    let s_re_atom = rational_const_to_atom(ctx, &s_re)?;
    let half = || ctx.pow(ctx.num(2), ctx.num(-1));
    let u_atom = ctx.pow(ctx.mul(&[ctx.add(&[abs_s, s_re_atom]), half()]), half());
    let v_atom = ctx.pow(
        ctx.mul(&[
            ctx.add(&[abs_s, ctx.mul(&[ctx.num(-1), s_re_atom])]),
            half(),
        ]),
        half(),
    );

    let mut out: Vec<DLogOrAtom<'a>> = Vec::new();
    for c in linear_roots {
        let v = d.gcd(&a.sub(&dp.mul_scalar(c)));
        out.push(DLogOrAtom::Log(c.clone(), v));
    }
    for sigma in [1i64, -1] {
        out.push(quartic_complex_pair(
            ctx, &v_alpha, s_re_atom, s_im, u_atom, v_atom, sigma, x,
        )?);
    }
    Some(out)
}

/// The realified pair for the quartic complex 4-cycle, one per sign of the
/// real part: the pair representative is the root `σ·u + i·v`, so the
/// contribution is `(σ·u)·log(U²+W²) − 2v·atan(W/U)` with `U + iW` the
/// factor `v` evaluated at that root.
#[allow(clippy::too_many_arguments)]
fn quartic_complex_pair<'a>(
    ctx: &'a AtomArena<'a>,
    v_alpha: &APoly,
    s_re: Atom<'a>,
    s_im: Atom<'a>,
    u_atom: Atom<'a>,
    v_atom: Atom<'a>,
    sigma: i64,
    x: Atom<'a>,
) -> Option<DLogOrAtom<'a>> {
    let dom = RationalDomain;
    let re = ctx.mul(&[ctx.num(sigma), u_atom]);
    let im = v_atom;
    // α⁰=1, α¹ = re + i·im, α² = s = s_re + i·(σ·s_im) — the square of
    // σ·u + i·v is (u²−v²) + 2σuv·i = s_re + i·σ·s_im — and α³ = α·α².
    let re2 = s_re;
    let im2 = ctx.mul(&[ctx.num(sigma), s_im]);
    let re3 = ctx.add(&[ctx.mul(&[re, re2]), ctx.mul(&[ctx.num(-1), im, im2])]);
    let im3 = ctx.add(&[ctx.mul(&[re, im2]), ctx.mul(&[im, re2])]);
    let re_powers = [ctx.num(1), re, re2, re3];
    let im_powers = [ctx.num(0), im, im2, im3];

    let mut u_terms = Vec::new();
    let mut w_terms = Vec::new();
    for (k, e) in v_alpha.coeffs().iter().enumerate() {
        // e = Σⱼ eⱼ·αʲ with rational eⱼ; split into real/imaginary parts.
        let mut re_acc = ctx.num(0);
        let mut im_acc = ctx.num(0);
        for (j, ej) in e.coeffs().iter().enumerate() {
            if dom.is_zero(ej) {
                continue;
            }
            let ec = rational_const_to_atom(ctx, ej)?;
            re_acc = ctx.add(&[re_acc, ctx.mul(&[ec, re_powers[j]])]);
            im_acc = ctx.add(&[im_acc, ctx.mul(&[ec, im_powers[j]])]);
        }
        let xk = match k {
            0 => ctx.num(1),
            1 => x,
            _ => ctx.pow(x, ctx.num(k as i64)),
        };
        u_terms.push(ctx.mul(&[re_acc, xk]));
        w_terms.push(ctx.mul(&[im_acc, xk]));
    }
    let u_poly = ctx.add(&u_terms);
    let w_poly = ctx.add(&w_terms);

    let sq_sum = ctx.add(&[ctx.pow(u_poly, ctx.num(2)), ctx.pow(w_poly, ctx.num(2))]);
    let log_term = ctx.mul(&[re, ctx.fun("log", &[sq_sum])]);
    // Reciprocal atan argument (canonical shape; formally equivalent).
    let atan = ctx.fun("atan", &[ctx.mul(&[u_poly, ctx.pow(w_poly, ctx.num(-1))])]);
    let atan_term = ctx.mul(&[ctx.num(2), im, atan]);
    Some(DLogOrAtom::Atom(ctx.add(&[log_term, atan_term])))
}

// ---------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------

fn int_to_rational(n: &Integer) -> Option<Rational> {
    Some(Rational::new(n.to_i64()?, 1))
}

fn rat_is_neg(r: &Rational) -> bool {
    use num_traits::Signed;
    r.inner().is_negative()
}

/// The exact square root of a positive rational, when it is rational.
fn sqrt_rational(r: &Rational) -> Option<Rational> {
    let p = r.numer().to_i64()?;
    let q = r.denom().to_i64()?;
    if p <= 0 || q <= 0 {
        return None;
    }
    let (sp, sq) = (p.checked_isqrt()?, q.checked_isqrt()?);
    if sp * sp == p && sq * sq == q {
        Some(Rational::new(sp, sq))
    } else {
        None
    }
}

/// Integer primitive part of a ℚ-polynomial (shared preamble with
/// `rational::rational_roots`). Also used by the tower-level algebraic
/// roots path (`logpart`).
pub(crate) fn integer_primitive(f: &DPoly) -> Option<ZPoly> {
    let dom = RationalDomain;
    if f.is_zero() || f.degree()? == 0 {
        return None;
    }
    let mut lcm: i64 = 1;
    for c in f.coeffs() {
        lcm = num_integer::lcm(lcm, c.denom().to_i64()?);
    }
    let zcoeffs: Option<Vec<Integer>> = f
        .coeffs()
        .iter()
        .map(|c| {
            let scaled = dom.mul(c, &Rational::new(lcm, 1));
            scaled.numer().to_i64().map(Integer::from)
        })
        .collect();
    let zpoly = ZPoly::from_coeffs(IntegerDomain, zcoeffs?);
    Some(zpoly.primitive_part())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::{AtomArena, AtomNode};
    use ocas_core::arena::Arena;
    use ocas_parse::parse;

    use crate::integral::rational::integrate_rational;
    use ocas_atom::Symbol;

    fn integrate_str<'a>(ctx: &'a AtomArena<'a>, src: &str) -> Option<Atom<'a>> {
        let expr = parse(ctx, src).expect("parse");
        integrate_rational(ctx, expr, Symbol::new("x"))
    }

    /// Minimal f64 evaluator for numeric derivative checks.
    fn eval_f64(atom: Atom, x: f64) -> f64 {
        match atom.node() {
            AtomNode::Num(n) => *n as f64,
            AtomNode::Var(_) => x,
            AtomNode::Add(args) => args.iter().map(|a| eval_f64(*a, x)).sum(),
            AtomNode::Mul(args) => args.iter().map(|a| eval_f64(*a, x)).product(),
            AtomNode::Pow(b, e) => eval_f64(*b, x).powf(eval_f64(*e, x)),
            AtomNode::Fun(name, args) => match name.as_str() {
                "log" => eval_f64(args[0], x).abs().ln(),
                "atan" => eval_f64(args[0], x).atan(),
                other => panic!("unsupported function in test eval: {other}"),
            },
        }
    }

    fn assert_numeric_antiderivative(result: Atom, integrand: Atom) {
        let h = 1e-6;
        for px in [0.7f64, 1.3, 2.1] {
            let approx = (eval_f64(result, px + h) - eval_f64(result, px - h)) / (2.0 * h);
            let exact = eval_f64(integrand, px);
            assert!(
                (approx - exact).abs() < 1e-4,
                "numeric derivative mismatch at x={px}: {approx} != {exact}"
            );
        }
    }

    #[test]
    fn quartic_biquadratic_one_over_one_plus_x4() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // The acceptance family: R(z) = 256z⁴+1, complex 4-cycle.
        let result = integrate_str(&ctx, "1/(1+x^4)").expect("solved");
        let s = result.to_string();
        assert!(!s.contains("Integral"), "still a residue: {s}");
        assert!(s.contains("atan"), "expected the realified atan form: {s}");
        let integrand = parse(&ctx, "1/(1+x^4)").unwrap();
        assert_numeric_antiderivative(result, integrand);
    }

    #[test]
    fn quadratic_real_roots_shared_field() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // R(z) = (20z²−10z+1)(20z²+10z+1): two quadratic factors sharing
        // the ℚ(√5) field; four real roots.
        let result = integrate_str(&ctx, "1/(1-3*x^2+x^4)").expect("solved");
        let s = result.to_string();
        assert!(!s.contains("Integral"), "still a residue: {s}");
        assert!(s.contains("log"), "expected log terms: {s}");
        let integrand = parse(&ctx, "1/(1-3*x^2+x^4)").unwrap();
        assert_numeric_antiderivative(result, integrand);
    }

    #[test]
    fn quadratic_complex_roots_pure_atan() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // R(z) = 4(24z²+1)²: repeated complex quadratic; the log part
        // vanishes (p = 0) and the answer is a pure atan form.
        let result = integrate_str(&ctx, "(1+x^2)/(1+4*x^2+x^4)").expect("solved");
        let s = result.to_string();
        assert!(!s.contains("Integral"), "still a residue: {s}");
        assert!(s.contains("atan"), "expected atan: {s}");
        let integrand = parse(&ctx, "(1+x^2)/(1+4*x^2+x^4)").unwrap();
        assert_numeric_antiderivative(result, integrand);
    }

    #[test]
    fn repeated_quartic_factor_realified() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // R(z) = (2304z⁴+1)²: repeated biquadratic; log+atan form with √6.
        let result = integrate_str(&ctx, "1/(1-x^4+x^8)").expect("solved");
        let s = result.to_string();
        assert!(!s.contains("Integral"), "still a residue: {s}");
        let integrand = parse(&ctx, "1/(1-x^4+x^8)").unwrap();
        assert_numeric_antiderivative(result, integrand);
    }

    #[test]
    fn hermite_residue_then_algebraic_logpart() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Hermite strips the squared factor; the residue 3/4·1/(1+x⁴) needs
        // the algebraic logarithmic part.
        let result = integrate_str(&ctx, "1/(x^4+1)^2").expect("solved");
        let s = result.to_string();
        assert!(!s.contains("Integral"), "still a residue: {s}");
        let integrand = parse(&ctx, "1/(x^4+1)^2").unwrap();
        assert_numeric_antiderivative(result, integrand);
    }

    #[test]
    fn out_of_fragment_degrees_decline() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // 1 + x + x³ is irreducible over ℚ: R(z) keeps a cubic factor,
        // which is outside the fragment — honest residue.
        let result = integrate_str(&ctx, "1/(1+x+x^3)").expect("partial");
        assert!(result.to_string().contains("Integral"));
    }
}
