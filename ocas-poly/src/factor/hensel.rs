//! Hensel lifting and Zassenhaus factor combination over $\mathbb{Z}$.
//!
//! This module lifts factorizations found modulo a prime back to the integers.
//! The finite-field factorization in [`super::finite_field`] is the foundation
//! it builds upon.
//!
//! Pipeline:
//! 1. Pick a prime $p$ not dividing the leading coefficient, with $f \bmod p$
//!    square-free.
//! 2. Factor $f \bmod p$ into monic irreducibles.
//! 3. Compute the Mignotte bound $B$; lift $p \to p^k$ with $p^k > 2B$ via
//!    linear Hensel lifting.
//! 4. Recombine the lifted factors by trial division (Zassenhaus).
//!
//! References: Zassenhaus (1969); Geddes, Czapor, Labahn, *Algorithms for
//! Computer Algebra*; Knuth, TAOCP vol. 2 §4.6.2.

use ocas_domain::{Domain, EuclideanDomain, FiniteField, Integer, IntegerDomain, number_theory};

use crate::dense::DenseUnivariatePolynomial;
use crate::factor::finite_field::{self, FpPoly};

/// Univariate polynomial over the integers.
pub type ZPoly = DenseUnivariatePolynomial<IntegerDomain>;

/// Convert a $\mathbb{Z}[x]$ polynomial to $\mathbb{F}_p[x]$ by reducing each
/// coefficient modulo $p$.
fn to_finite_field(f: &ZPoly, p: &Integer) -> FpPoly {
    let field = FiniteField::new(p.to_bigint());
    let coeffs = f
        .coeffs()
        .iter()
        .map(|c| field.element(c.to_bigint()))
        .collect();
    FpPoly::from_coeffs(field, coeffs)
}

/// Convert an $\mathbb{F}_p[x]$ polynomial back to $\mathbb{Z}[x]$ using the
/// symmetric representative of each coefficient (range $(-p/2, p/2]$).
fn from_finite_field_symmetric(g: &FpPoly) -> ZPoly {
    let domain = IntegerDomain;
    let p = g.domain().prime();
    let p_int = Integer::from(p.clone());
    let coeffs = g
        .coeffs()
        .iter()
        .map(|c| {
            let c_int = Integer::from(c.value().clone());
            number_theory::symmetric_mod(&c_int, &p_int)
        })
        .collect();
    ZPoly::from_coeffs(domain, coeffs)
}

/// Normalize an $\mathbb{F}_p[x]$ polynomial to monic form.
fn monic_fp(f: &FpPoly) -> FpPoly {
    if f.is_zero() {
        return f.zero();
    }
    let lc = f.leading_coeff().unwrap();
    if f.domain().is_one(lc) {
        return f.clone();
    }
    let inv = f.domain().inv(lc).expect("nonzero lc");
    f.mul_scalar(&inv)
}

/// Landau–Mignotte bound: an upper bound on the absolute value of any
/// coefficient of a factor of $f$ in $\mathbb{Z}[x]$.
///
/// For a degree-$n$ polynomial $f$ with coefficient 2-norm $\|f\|_2$, every
/// factor $g$ satisfies $\|g\|_\infty \le 2^n \|f\|_2$.
pub(crate) fn mignotte_bound(f: &ZPoly) -> Integer {
    let n = f.degree().unwrap_or(0);
    let mut sum_sq = Integer::from(0);
    for c in f.coeffs() {
        let v = c.abs();
        sum_sq += &(&v * &v);
    }
    let norm = sum_sq.sqrt() + &Integer::from(1);
    &Integer::from(2).pow_u32(n as u32) * &norm
}

/// Bézout coefficients over $\mathbb{F}_p$ for coprime `g`, `h`: returns
/// `(s, t)` with $s g + t h \equiv 1 \pmod p$.
///
/// Over a field the GCD is only defined up to a unit, so the raw extended
/// Euclid yields `s·g + t·h = c` for some nonzero constant `c`; we rescale by
/// `1/c` to obtain the normalized identity.
fn bezout_mod_p(g: &FpPoly, h: &FpPoly) -> (FpPoly, FpPoly) {
    let field = g.domain().clone();
    let one = FpPoly::from_coeffs(field.clone(), vec![field.one()]);
    let zero = FpPoly::from_coeffs(field.clone(), vec![field.zero()]);
    let (mut old_r, mut r) = (g.clone(), h.clone());
    let (mut old_s, mut s) = (one.clone(), zero.clone());
    let (mut old_t, mut t) = (zero, one);
    while !r.is_zero() {
        let (q, rem) = old_r.div_rem(&r).unwrap_or_else(|| (r.zero(), r.zero()));
        old_r = r;
        r = rem;
        let new_s = old_s.sub(&q.mul(&s));
        let new_t = old_t.sub(&q.mul(&t));
        old_s = s;
        s = new_s;
        old_t = t;
        t = new_t;
    }
    // Normalize so that s·g + t·h = 1 (old_r may be any nonzero constant c).
    if let Some(c) = old_r.leading_coeff()
        && !field.is_one(c)
    {
        let c_inv = field.inv(c).expect("gcd constant is nonzero");
        old_s = old_s.mul_scalar(&c_inv);
        old_t = old_t.mul_scalar(&c_inv);
    }
    (old_s, old_t)
}

/// Lift a two-factor factorization $f \equiv g \cdot h \pmod p$ (with
/// $\gcd(g, h) = 1$ over $\mathbb{F}_p$) towards $f = g \cdot h$ over
/// $\mathbb{Z}$, for monic $f$ with monic true factors.
///
/// Lifts modulo $p^k$ until $p^k > \text{bound}$, then returns the lifted
/// pair. Returns `None` if the lift is inconsistent.
fn hensel_lift_pair(
    f: &ZPoly,
    g0: &FpPoly,
    h0: &FpPoly,
    p: &Integer,
    bound: &Integer,
) -> Option<(ZPoly, ZPoly)> {
    // Bézout coefficient t for h0 in 1 = s·g0 + t·h0 (mod p); only t is used.
    let (s, t) = bezout_mod_p(g0, h0);
    debug_assert!(
        {
            let one = g0.one();
            let id = s.mul(g0).add(&t.mul(h0));
            id.sub(&one).is_zero()
        },
        "Bézout identity s·g0 + t·h0 = 1 failed"
    );
    let mut g = from_finite_field_symmetric(g0);
    let mut h = from_finite_field_symmetric(h0);
    let mut m = p.clone();

    loop {
        let gh = g.mul(&h);
        let e = f.sub(&gh);
        if e.is_zero() {
            return Some((g, h));
        }
        // e must be divisible by m coefficientwise.
        let mut e_over_m = Vec::new();
        for c in e.coeffs() {
            let (q, r) = c.div_rem(&m);
            if r.is_zero() {
                e_over_m.push(q);
            } else {
                return None;
            }
        }
        let e_over_m = ZPoly::from_coeffs(IntegerDomain, e_over_m);
        if e_over_m.is_zero() {
            return Some((g, h));
        }
        let e_bar = to_finite_field(&e_over_m, p);
        // Δg = (t·ē) mod g0  (degree < deg g0).
        let (_q1, dg) = t.mul(&e_bar).div_rem(g0)?;
        // Δh = (ē − Δg·h0) / g0  (exact over the field, degree < deg h0).
        // This direct division avoids the large-intermediate cancellation
        // that the q·h0 + s·ē form suffers when deg(h0) ≫ deg(g0).
        let dividend = e_bar.sub(&dg.mul(h0));
        let (dh, dh_rem) = dividend.div_rem(g0)?;
        debug_assert!(
            dh_rem.is_zero(),
            "Δh division not exact; Bézout identity may be broken"
        );
        let dg_z = from_finite_field_symmetric(&dg);
        let dh_z = from_finite_field_symmetric(&dh);
        let m_int = m.clone();
        g = g.add(&dg_z.mul_scalar(&m_int));
        h = h.add(&dh_z.mul_scalar(&m_int));
        m *= p;
        if &m > bound {
            return Some((g, h));
        }
    }
}

/// Lift many monic mod-$p$ factors back to $\mathbb{Z}$ by peeling off one
/// factor at a time (each step is a pairwise Hensel lift of `g0` against the
/// product of the remaining factors). `bound` is a power of $p$ exceeding
/// $2\,\mathrm{Mignotte}(f)$; the returned factors are reduced into the
/// symmetric range of `bound` for subsequent Zassenhaus recombination.
fn hensel_lift_multi(
    f: &ZPoly,
    factors: &[FpPoly],
    p: &Integer,
    bound: &Integer,
) -> Option<Vec<ZPoly>> {
    if factors.len() <= 1 {
        return Some(vec![f.clone()]);
    }
    let mut lifted: Vec<ZPoly> = Vec::new();
    let mut work = factors.to_vec();
    let mut f_current = f.clone();
    while work.len() > 1 {
        let g0 = monic_fp(&work[0].clone());
        let h0 = monic_fp(&work[1..].iter().cloned().reduce(|a, b| a.mul(&b)).unwrap());
        let (g, h) = hensel_lift_pair(&f_current, &g0, &h0, p, bound)?;
        // f_current must be kept mod-p faithful, so reduce only the emitted
        // factor g; carry h forward unreduced.
        lifted.push(reduce_symmetric(&g, bound));
        f_current = h;
        work = work[1..].to_vec();
    }
    lifted.push(reduce_symmetric(&f_current, bound));
    Some(lifted)
}

/// Reduce each coefficient of a $\mathbb{Z}[x]$ polynomial into the symmetric
/// range $(-M/2, M/2]$ modulo $M$. Used after Hensel lifting to recover the
/// true integer factors from their modular images.
fn reduce_symmetric(f: &ZPoly, modulus: &Integer) -> ZPoly {
    let coeffs = f
        .coeffs()
        .iter()
        .map(|c| number_theory::symmetric_mod(c, modulus))
        .collect();
    ZPoly::from_coeffs(IntegerDomain, coeffs)
}

/// Generate all index combinations of `k` elements from `0..n`.
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    combos(0, n, k, &mut cur, &mut out);
    out
}

fn combos(start: usize, n: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if cur.len() == k {
        out.push(cur.clone());
        return;
    }
    for i in start..n {
        cur.push(i);
        combos(i + 1, n, k, cur, out);
        cur.pop();
    }
}

/// Zassenhaus factor combination: enumerate subsets of the lifted factors in
/// order of increasing size and accept the first that divides $f$ in
/// $\mathbb{Z}[x]$. Each candidate's coefficients are reduced into the symmetric
/// range of the lifting modulus before trial-division, since a true
/// $\mathbb{Z}$-factor $h$ (with $\|h\|_\infty < B/2$) equals the subset
/// product modulo $B$. If `f` is not monic, an integer divisor of the leading
/// coefficient is attached to the candidate before testing. Returns the monic
/// irreducible factors.
fn zassenhaus_combine(f: &ZPoly, lifted: &[ZPoly], modulus: &Integer) -> Vec<ZPoly> {
    if lifted.is_empty() {
        return Vec::new();
    }
    let one = f.one();
    let mut remaining: Vec<ZPoly> = lifted.to_vec();
    let mut result = Vec::new();
    let mut rest = f.clone();
    let mut size = 1usize;
    while size <= remaining.len() && !remaining.is_empty() {
        let n = remaining.len();
        let mut found = false;
        for combo in combinations(n, size) {
            let mut prod = one.clone();
            for &idx in &combo {
                prod = prod.mul(&remaining[idx]);
            }
            // Scale by the leading coefficient of the current cofactor, then
            // take the primitive part: a true integer factor g satisfies
            // lc(rest)·∏(subset) = c·g for some content c (Zassenhaus).
            let scaled = prod.mul_scalar(
                &rest
                    .leading_coeff()
                    .cloned()
                    .unwrap_or_else(|| Integer::from(1)),
            );
            let candidate = reduce_symmetric(&scaled, modulus);
            let content = candidate
                .coeffs()
                .iter()
                .fold(Integer::from(0), |acc, c| IntegerDomain.gcd(&acc, &c.abs()));
            let content = if content.is_zero() {
                Integer::from(1)
            } else {
                content.abs()
            };
            let primitive = if content == Integer::from(1) {
                candidate.clone()
            } else {
                let coeffs = candidate
                    .coeffs()
                    .iter()
                    .map(|c| IntegerDomain.div(c, &content).unwrap_or_else(|| c.clone()))
                    .collect();
                ZPoly::from_coeffs(IntegerDomain, coeffs)
            };
            if primitive.is_one() || primitive.is_zero() {
                continue;
            }
            if let Some((q, r)) = rest.div_rem(&primitive)
                && r.is_zero()
            {
                result.push(primitive);
                let mut nr = Vec::new();
                for (i, fac) in remaining.iter().enumerate() {
                    if !combo.contains(&i) {
                        nr.push(fac.clone());
                    }
                }
                remaining = nr;
                rest = q;
                found = true;
                size = 1;
                break;
            }
        }
        if !found {
            size += 1;
        }
    }
    // The final cofactor is the last factor (primitive part).
    if !rest.is_one() && !rest.is_zero() {
        let content = rest
            .coeffs()
            .iter()
            .fold(Integer::from(0), |acc, c| IntegerDomain.gcd(&acc, &c.abs()));
        let content = if content.is_zero() {
            Integer::from(1)
        } else {
            content.abs()
        };
        let last = if content == Integer::from(1) {
            rest.clone()
        } else {
            let coeffs = rest
                .coeffs()
                .iter()
                .map(|c| IntegerDomain.div(c, &content).unwrap_or_else(|| c.clone()))
                .collect();
            ZPoly::from_coeffs(IntegerDomain, coeffs)
        };
        if !last.is_one() {
            result.push(last);
        }
    }
    result
}

/// Factor a monic square-free polynomial over $\mathbb{Z}$ into monic
/// irreducible factors. The input must be monic in its leading variable;
/// non-monic polynomials are returned unfactored (the multivariate driver
/// works over a finite field where monic normalization is always possible).
pub fn factor_square_free_monic(f: &ZPoly) -> Vec<ZPoly> {
    if f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    if f.degree() == Some(1) {
        return vec![f.clone()];
    }
    let lc_abs = f.leading_coeff().unwrap().abs();
    if lc_abs != Integer::from(1) {
        // Hensel lifting below assumes a monic leading coefficient.
        return vec![f.clone()];
    }
    let bound = mignotte_bound(f);
    let two_bound = &Integer::from(2) * &bound;
    let lc = f.leading_coeff().unwrap().abs();

    let mut prime_count = 0usize;
    for p in number_theory::primes_from(&Integer::from(2)) {
        prime_count += 1;
        if prime_count > 30 {
            break;
        }
        if (&lc % &p).is_zero() {
            continue;
        }
        let fp = to_finite_field(f, &p);
        let fpp = fp.derivative();
        if fpp.is_zero() || fp.gcd(&fpp).degree().unwrap_or(0) > 0 {
            continue; // not square-free mod p
        }
        let factors_p = finite_field::cantor_zassenhaus(&monic_fp(&fp));
        if factors_p.len() <= 1 {
            return vec![f.clone()]; // irreducible over Z
        }
        let mut lift_mod = p.clone();
        while lift_mod <= two_bound {
            lift_mod *= &p;
        }
        if let Some(lifted) = hensel_lift_multi(f, &factors_p, &p, &lift_mod) {
            let irreducibles = zassenhaus_combine(f, &lifted, &lift_mod);
            if !irreducibles.is_empty() {
                return irreducibles;
            }
        }
    }
    vec![f.clone()]
}

/// Factor a square-free (not necessarily monic) primitive polynomial over
/// $\mathbb{Z}$ into irreducible factors.
///
/// Non-monic inputs are handled by the leading-coefficient transformation:
/// for $f$ of degree $d$ with leading coefficient $a$, the polynomial
/// $a^{d-1} f(x/a)$ is monic and is factored via [`factor_square_free_monic`];
/// each monic factor $G$ maps back to the primitive part of $G(a x)$.
pub fn factor_square_free(f: &ZPoly) -> Vec<ZPoly> {
    let d = f.degree().unwrap_or(0);
    if d == 0 {
        return Vec::new();
    }
    if d == 1 {
        return vec![f.clone()];
    }
    let a = f.leading_coeff().cloned().unwrap();
    if a.abs() == Integer::from(1) {
        return factor_square_free_monic(f);
    }
    // g(x) = a^{d-1} f(x/a) = Σ_k c_k·a^{d-1-k}·x^k. For k = d this is
    // c_d·a^{-1} = a/a = 1, so g is monic; for k < d, g_k = c_k·a^{d-1-k}.
    let mut g_coeffs = vec![Integer::from(0); d + 1];
    for (k, c) in f.coeffs().iter().enumerate().take(d + 1) {
        if IntegerDomain.is_zero(c) {
            continue;
        }
        g_coeffs[k] = if k == d {
            Integer::from(1)
        } else {
            IntegerDomain.mul(c, &a.pow_u32((d - 1 - k) as u32))
        };
    }
    let g = ZPoly::from_coeffs(IntegerDomain, g_coeffs);
    let monic_factors = factor_square_free_monic(&g);
    if monic_factors.len() <= 1 && g != *f {
        // Transformation produced no split; try direct (monic-only) as a
        // fallback for degenerate cases.
        return vec![f.clone()];
    }
    // Map back: H(x) = primitive part of G(a x).
    let mut out = Vec::new();
    for gm in &monic_factors {
        let mut h = ZPoly::from_coeffs(IntegerDomain, vec![Integer::from(0)]);
        for (k, c) in gm.coeffs().iter().enumerate() {
            if IntegerDomain.is_zero(c) {
                continue;
            }
            let a_pow = a.pow_u32(k as u32);
            let coeff = IntegerDomain.mul(c, &a_pow);
            let term = ZPoly::from_coeffs(IntegerDomain, {
                let mut v = vec![Integer::from(0); k + 1];
                v[k] = coeff;
                v
            });
            h = h.add(&term);
        }
        let hp = h.primitive_part();
        // Normalize sign so the leading coefficient matches f's sign pattern.
        if hp.leading_coeff().is_some_and(|l| l.is_negative()) {
            out.push(hp.mul_scalar(&Integer::from(-1)));
        } else {
            out.push(hp);
        }
    }
    // Verify reconstruction; fall back to unfactored on any mismatch.
    let mut prod = f.one();
    for h in &out {
        prod = prod.mul(h);
    }
    let (q, r) = f.div_rem(&prod).unwrap_or((f.zero(), f.clone()));
    if r.is_zero() && q.degree() == Some(0) {
        out
    } else {
        vec![f.clone()]
    }
}

/// Factor a primitive polynomial over $\mathbb{Z}$ into irreducible factors
/// with multiplicities, using square-free factorization (characteristic 0,
/// so the generic Yun algorithm applies) followed by
/// [`factor_square_free_monic`] on each square-free component.
pub fn factor_primitive(f: &ZPoly) -> Vec<(ZPoly, usize)> {
    if f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    // Square-free factorization over Z (generic Yun works in characteristic 0).
    let sqfree = f.square_free_factorization();
    let mut result = Vec::new();
    for (g, mult) in sqfree {
        // Normalize sign so the factor is monic-ish: if leading coeff is
        // negative, negate the polynomial (absorbed into the content/sign).
        let lc = g.leading_coeff().cloned().unwrap();
        let sign = if lc.is_negative() {
            Integer::from(-1i64)
        } else {
            Integer::from(1i64)
        };
        let g_pos = g.mul_scalar(&sign);
        for irr in factor_square_free(&g_pos.primitive_part()) {
            result.push((irr, mult));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zpoly(coeffs: &[i64]) -> ZPoly {
        ZPoly::from_coeffs(
            IntegerDomain,
            coeffs.iter().map(|&c| Integer::from(c)).collect(),
        )
    }

    fn reconstruct(f: &ZPoly, factors: &[ZPoly]) -> bool {
        let mut acc = f.one();
        for g in factors {
            acc = acc.mul(g);
        }
        let (q, r) = f.div_rem(&acc).unwrap_or((f.zero(), f.clone()));
        r.is_zero() && q.degree() == Some(0)
    }

    #[test]
    fn factor_x100_minus_1_over_z() {
        let mut coeffs = vec![Integer::from(-1i64)];
        coeffs.resize(101, Integer::from(0));
        coeffs[100] = Integer::from(1);
        let f = ZPoly::from_coeffs(IntegerDomain, coeffs);
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        // 9 cyclotomic factors for d | 100.
        assert_eq!(factors.len(), 9);
    }

    #[test]
    fn factor_quadratic_split() {
        let f = zpoly(&[6, 5, 1]); // (x+2)(x+3)
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        assert_eq!(factors.len(), 2);
    }

    #[test]
    fn factor_irreducible_quadratic() {
        let f = zpoly(&[1, 0, 1]); // x^2+1
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        assert_eq!(factors.len(), 1);
    }

    #[test]
    fn factor_three_linear() {
        let f = zpoly(&[-6, 11, -6, 1]); // (x-1)(x-2)(x-3)
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        assert_eq!(factors.len(), 3);
    }

    #[test]
    fn factor_three_mixed() {
        // (x^2+1)(x^2+x+1)(x+1)
        let a = zpoly(&[1, 0, 1]);
        let b = zpoly(&[1, 1, 1]);
        let c = zpoly(&[1, 1]);
        let f = a.mul(&b).mul(&c);
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        assert_eq!(factors.len(), 3);
    }

    #[test]
    fn mignotte_bound_sanity() {
        let f = zpoly(&[1, 0, 1]); // x^2 + 1, ||f||_2 = sqrt(2)
        let b = mignotte_bound(&f);
        // True bound = 2^2 * sqrt(2) ≈ 5.66; conservative integer version
        // rounds up and is therefore >= 6.
        assert!(
            b >= Integer::from(6) && b <= Integer::from(10),
            "mignotte(x^2+1) = {b}, expected ~6-10"
        );
    }

    #[test]
    fn factor_cyclotomic_matches_sympy_over_z() {
        // Ground-truth degree histograms from SymPy 1.14
        // `Poly(x^n-1, x, domain='ZZ').factor_list()` for the cyclotomic
        // decomposition of x^n - 1 over Z.
        let cases: &[(usize, &[(usize, usize)])] = &[
            (12, &[(1, 2), (2, 3), (4, 1)]),
            (30, &[(1, 2), (2, 2), (4, 2), (8, 2)]),
            (60, &[(1, 2), (2, 3), (4, 3), (8, 3), (16, 1)]),
            (100, &[(1, 2), (2, 1), (4, 2), (8, 1), (20, 2), (40, 1)]),
        ];
        for &(n, expected) in cases {
            let mut coeffs = vec![Integer::from(-1i64)];
            coeffs.resize(n + 1, Integer::from(0));
            coeffs[n] = Integer::from(1);
            let f = ZPoly::from_coeffs(IntegerDomain, coeffs);
            let factors = factor_square_free_monic(&f);
            assert!(
                reconstruct(&f, &factors),
                "x^{n}-1: factors do not reconstruct"
            );
            let mut obs: std::collections::BTreeMap<usize, usize> = Default::default();
            for g in &factors {
                *obs.entry(g.degree().unwrap()).or_insert(0) += 1;
            }
            let observed: Vec<(usize, usize)> = obs.into_iter().collect();
            assert_eq!(
                observed.as_slice(),
                expected,
                "x^{n}-1 over Z: degree histogram mismatch"
            );
        }
    }

    #[test]
    fn factor_x30_minus_1_over_z() {
        // x^30 - 1 has 8 cyclotomic factors (d | 30: 1,2,3,5,6,10,15,30).
        let mut coeffs = vec![Integer::from(-1i64)];
        coeffs.resize(31, Integer::from(0));
        coeffs[30] = Integer::from(1);
        let f = ZPoly::from_coeffs(IntegerDomain, coeffs);
        let factors = factor_square_free_monic(&f);
        assert!(reconstruct(&f, &factors));
        assert_eq!(factors.len(), 8, "expected 8 cyclotomic factors");
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn any_int_poly(max_degree: usize) -> impl Strategy<Value = ZPoly> {
        (0..=max_degree)
            .prop_flat_map(move |deg| {
                let n = deg + 1;
                prop::collection::vec(i64_range(), n)
            })
            .prop_map(|coeffs| {
                let c: Vec<Integer> = coeffs.into_iter().map(Integer::from).collect();
                ZPoly::from_coeffs(IntegerDomain, c)
            })
    }

    fn i64_range() -> impl Strategy<Value = i64> {
        -100i64..=100i64
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(500))]

        #[test]
        fn factor_then_multiply_roundtrip(p in any_int_poly(6)) {
            // The input may have content > 1; we factor the primitive part.
            let f = p.primitive_part();
            if f.degree().unwrap_or(0) == 0 {
                return Ok(());
            }
            let factors = factor_primitive(&f);
            let mut acc = f.one();
            for (g, m) in &factors {
                for _ in 0..*m {
                    acc = acc.mul(g);
                }
            }
            // acc and f should be equal up to a constant factor.
            if let Some((q, r)) = f.div_rem(&acc) {
                prop_assert!(r.is_zero());
                prop_assert_eq!(q.degree(), Some(0));
            }
        }
    }
}
