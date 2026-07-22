//! Factorization over algebraic number fields $\mathbb{Q}(\alpha)$
//! (Trager's algorithm).
//!
//! The norm of $f \in K[x]$ down to $\mathbb{Q}[x]$ is computed by
//! evaluation–interpolation of the scalar resultant
//! $\operatorname{Res}_\alpha(m(\alpha), f(x, \alpha))$; Trager's shift
//! finds $s \ge 0$ such that the norm of $f(x - s\alpha)$ is square-free.
//! Factoring the norm over $\mathbb{Q}$ and taking GCDs over $K$ recovers
//! the irreducible factors of $f$.

use ocas_domain::number_theory::{crt, primes_from};
use ocas_domain::{
    AlgebraicElement, AlgebraicExtension, AlgebraicNumberField, Domain, EuclideanDomain,
    FiniteField, Integer, IntegerDomain, Rational, RationalDomain,
};

use crate::dense::DenseUnivariatePolynomial;
use crate::factor::finite_field::factor_over_finite_field;
use crate::factor::hensel::factor_primitive;
use crate::rational_reconstruction::rational_reconstruction;

/// Dense univariate polynomial over domain `D`.
type UP<D> = DenseUnivariatePolynomial<D>;
/// An element of an algebraic number field.
type AnfElement = AlgebraicElement<Rational>;

/// Compose `f` with the linear polynomial `x + shift` (Horner's scheme).
pub(crate) fn compose_linear<D: Domain>(f: &UP<D>, shift: &D::Element) -> UP<D> {
    let domain = f.domain().clone();
    let linear = UP::from_coeffs(domain.clone(), vec![shift.clone(), domain.one()]);
    let mut acc = UP::from_coeffs(domain.clone(), Vec::new());
    for c in f.coeffs().iter().rev() {
        acc = acc
            .mul(&linear)
            .add(&UP::from_coeffs(domain.clone(), vec![c.clone()]));
    }
    acc
}

/// Factor a non-constant **square-free** univariate polynomial over
/// $\mathbb{Q}$ into irreducible factors: clear denominators, take the
/// primitive part, and factor over $\mathbb{Z}$ (Hensel path — no
/// pseudo-remainder sequences).
///
/// The rational unit (content over $\mathbb{Q}$) is dropped, since units do
/// not affect factorization over a field.
pub(crate) fn factor_square_free_rationals(f: &UP<RationalDomain>) -> Vec<UP<RationalDomain>> {
    if f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    let primitive = clear_denominators(f).primitive_part();
    crate::factor::hensel::factor_square_free(&primitive)
        .into_iter()
        .map(|g| {
            UP::<RationalDomain>::from_coeffs(
                RationalDomain,
                g.coeffs()
                    .iter()
                    .map(|c| Rational::from_integer(c.clone()))
                    .collect(),
            )
        })
        .collect()
}

/// Interpolate the unique polynomial of degree `< xs.len()` through the
/// points `(xs[i], ys[i])` via Newton's divided differences.
fn interpolate_rational(xs: &[Rational], ys: &[Rational]) -> UP<RationalDomain> {
    let q = RationalDomain;
    debug_assert_eq!(xs.len(), ys.len());
    let n = xs.len();
    let mut coef: Vec<Rational> = ys.to_vec();
    for k in 1..n {
        for i in (k..n).rev() {
            let num = q.sub(&coef[i], &coef[i - 1]);
            let den = q.sub(&xs[i], &xs[i - k]);
            coef[i] = q.div(&num, &den).expect("distinct interpolation points");
        }
    }
    // Expand the Newton form into the standard basis.
    let mut poly = UP::<RationalDomain>::from_coeffs(q, vec![coef[n - 1].clone()]);
    for k in (0..n - 1).rev() {
        let factor = UP::<RationalDomain>::from_coeffs(q, vec![q.neg(&xs[k]), q.one()]);
        poly = poly
            .mul(&factor)
            .add(&UP::<RationalDomain>::from_coeffs(q, vec![coef[k].clone()]));
    }
    poly
}

/// Norm of `g` from $K[x]$ down to $\mathbb{Q}[x]$ by
/// evaluation–interpolation: at $\deg_x(g) \cdot [K:\mathbb{Q}] + 1$
/// rational points $x_j$ compute the scalar resultant
/// $\operatorname{Res}_\alpha(m, g(x_j, \alpha))$ (the norm of the value),
/// then interpolate. Exact because the norm has degree
/// $\deg_x(g) \cdot [K:\mathbb{Q}]$.
fn norm_eval_interp(
    field: &AlgebraicNumberField,
    g: &UP<AlgebraicNumberField>,
) -> UP<RationalDomain> {
    let deg_x = g.degree().unwrap_or(0);
    let d = field.extension_degree();
    let n_points = d * deg_x + 1;
    let m_poly = UP::<RationalDomain>::from_coeffs(RationalDomain, field.min_poly().to_vec());
    let mut xs = Vec::with_capacity(n_points);
    let mut ys = Vec::with_capacity(n_points);
    for j in 0..n_points {
        // Symmetric points 0, 1, −1, 2, −2, … keep the arithmetic small.
        let v = ((j as i64 + 1) / 2) * if j % 2 == 1 { -1 } else { 1 };
        let x = Rational::new(v, 1);
        let x_elem = field.from_base(x.clone());
        // Horner evaluation of g at x inside the extension field.
        let mut acc = field.zero();
        for c in g.coeffs().iter().rev() {
            acc = field.add(&field.mul(&acc, &x_elem), c);
        }
        let acc_poly = UP::<RationalDomain>::from_coeffs(RationalDomain, acc.coeffs().to_vec());
        ys.push(m_poly.resultant(&acc_poly));
        xs.push(x);
    }
    interpolate_rational(&xs, &ys)
}

/// Clear the denominators of a rational polynomial: returns the integer
/// polynomial `lcm(denoms) · f`.
fn clear_denominators(f: &UP<RationalDomain>) -> UP<IntegerDomain> {
    let zd = IntegerDomain;
    let mut scale = Integer::from(1);
    for c in f.coeffs() {
        let den = c.denom();
        let g = zd.gcd(&scale, &den);
        let prod = zd.mul(&scale, &den);
        scale = zd.div_rem(&prod, &g).expect("gcd divides product").0;
    }
    let scale_q = Rational::from_integer(scale);
    UP::<IntegerDomain>::from_coeffs(
        zd,
        f.coeffs()
            .iter()
            .map(|c| RationalDomain.mul(c, &scale_q).numer())
            .collect(),
    )
}

/// Square-freeness of a rational polynomial via modular checks: if `f mod
/// p` (degree-preserving) is square-free over $\mathbb{F}_p$ for some prime
/// `p`, then `f` is square-free over $\mathbb{Q}$ — an exact accept.
/// Rejection after all primes is heuristic (a square-free `f` is only
/// rejected when every listed prime divides its discriminant); the Trager
/// shift loop simply advances to the next shift in that case.
///
/// This avoids the pseudo-remainder GCD over $\mathbb{Z}$/$\mathbb{Q}$,
/// whose coefficients explode already at moderate degrees.
fn is_square_free_rational(f: &UP<RationalDomain>) -> bool {
    if f.degree().unwrap_or(0) <= 1 {
        return true;
    }
    const PRIMES: [u64; 8] = [101, 103, 107, 109, 113, 127, 131, 137];
    for p in PRIMES {
        let fp = FiniteField::new(num_bigint::BigInt::from(p));
        let mut coeffs = Vec::with_capacity(f.coeffs().len());
        let mut ok = true;
        for c in f.coeffs() {
            match map_rational_fp(c, &fp) {
                Some(v) => coeffs.push(v),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue; // p divides a denominator
        }
        let f_p = UP::<FiniteField>::from_coeffs(fp, coeffs);
        if f_p.degree() != f.degree() {
            continue; // p divides the leading coefficient
        }
        let g = f_p.gcd(&f_p.derivative());
        if g.degree() == Some(0) {
            return true;
        }
    }
    false
}

/// Maximum number of shifts tried in [`norm_with_shift`] before giving up.
const MAX_TRAGER_SHIFTS: u64 = 100;

/// Maximum number of primes tried in [`gcd_anf`] before falling back to the
/// dense Euclidean GCD over the number field.
const MAX_ANF_GCD_PRIMES: usize = 1000;

/// Galois field $\mathrm{GF}(p^d)$ used for modular images.
type GaloisField = AlgebraicExtension<FiniteField>;

/// Map a rational coefficient into $\mathbb{F}_p$: numerator times the
/// inverse of the denominator. Returns `None` when the denominator
/// vanishes mod `p`.
fn map_rational_fp(c: &Rational, fp: &FiniteField) -> Option<ocas_domain::FiniteFieldElement> {
    let num = fp.element(c.numer().to_bigint());
    let den = fp.element(c.denom().to_bigint());
    fp.div(&num, &den)
}

/// Map an ANF element (α-polynomial with rational coefficients) into
/// $\mathrm{GF}(p^d)$. Returns `None` when a denominator vanishes mod `p`.
fn map_element_gf(
    e: &AnfElement,
    gf: &GaloisField,
) -> Option<AlgebraicElement<ocas_domain::FiniteFieldElement>> {
    let fp = gf.base_domain();
    let mut coeffs = Vec::with_capacity(e.coeffs().len());
    for c in e.coeffs() {
        coeffs.push(map_rational_fp(c, fp)?);
    }
    Some(gf.element(coeffs))
}

/// Map a univariate ANF polynomial into $\mathrm{GF}(p^d)[x]$.
fn map_poly_gf(f: &UP<AlgebraicNumberField>, gf: &GaloisField) -> Option<UP<GaloisField>> {
    let mut coeffs = Vec::with_capacity(f.coeffs().len());
    for c in f.coeffs() {
        coeffs.push(map_element_gf(c, gf)?);
    }
    Some(UP::<GaloisField>::from_coeffs(gf.clone(), coeffs))
}

/// Check whether the minimal polynomial of `field` is irreducible mod `p`
/// (and all its coefficient denominators are units mod `p`), so that the
/// quotient $\mathbb{F}_p[\alpha]/(m)$ is indeed $\mathrm{GF}(p^d)$.
fn min_poly_irreducible_mod(field: &AlgebraicNumberField, fp: &FiniteField) -> bool {
    let mut coeffs = Vec::with_capacity(field.min_poly().len());
    for c in field.min_poly() {
        match map_rational_fp(c, fp) {
            Some(v) => coeffs.push(v),
            None => return false,
        }
    }
    let m_fp = UP::<FiniteField>::from_coeffs(fp.clone(), coeffs);
    let d = field.extension_degree();
    if m_fp.degree() != Some(d) {
        return false;
    }
    let factors = factor_over_finite_field(&m_fp);
    factors.len() == 1 && factors[0].1 == 1 && factors[0].0.degree() == Some(d)
}

/// Make a polynomial over a field monic by scaling with the inverse of the
/// leading coefficient.
fn monic_over<D: EuclideanDomain>(f: &UP<D>) -> UP<D> {
    match f.leading_coeff() {
        Some(lc) if !f.domain().is_one(lc) => {
            let inv = f.domain().inv(lc).expect("nonzero leading coefficient");
            f.mul_scalar(&inv)
        }
        _ => f.clone(),
    }
}

/// GCD of two univariate polynomials over the number field $K =
/// \mathbb{Q}(\alpha)$ via the modular method (Encarnación): map to
/// $\mathrm{GF}(p^d)$ for primes with $m$ irreducible mod $p$, combine
/// monic modular GCDs by CRT, rational-reconstruct the coefficients, and
/// verify by trial division over $K$.
///
/// Primes are unlucky when the modular GCD degree exceeds the true degree;
/// they are detected by comparing degrees across primes and discarded.
/// After [`MAX_ANF_GCD_PRIMES`] primes the dense Euclidean GCD over $K$ is
/// used as a correctness fallback.
///
/// Returns the monic GCD.
pub(crate) fn gcd_anf(
    field: &AlgebraicNumberField,
    a: &UP<AlgebraicNumberField>,
    b: &UP<AlgebraicNumberField>,
) -> UP<AlgebraicNumberField> {
    if a.is_zero() {
        return monic_over(b);
    }
    if b.is_zero() {
        return monic_over(a);
    }
    if a.degree() == Some(0) || b.degree() == Some(0) {
        return UP::<AlgebraicNumberField>::from_coeffs(field.clone(), vec![field.one()]);
    }
    let d = field.extension_degree();
    // CRT state: residues of the monic modular GCD coefficients,
    // [x-degree][α-degree], together with the modulus.
    let mut residues: Option<Vec<Vec<Integer>>> = None;
    let mut modulus = Integer::from(1);
    let mut primes = primes_from(&Integer::from(1));
    for _ in 0..MAX_ANF_GCD_PRIMES {
        let p_int = primes.next().expect("infinite prime iterator");
        let fp = FiniteField::new(p_int.to_bigint());
        // Skip primes where a denominator or the minimal polynomial's
        // irreducibility breaks down.
        if !min_poly_irreducible_mod(field, &fp) {
            continue;
        }
        let gf = GaloisField::new(
            fp.clone(),
            field
                .min_poly()
                .iter()
                .map(|c| map_rational_fp(c, &fp).expect("denominators checked"))
                .collect(),
        );
        // Skip primes dividing a leading coefficient.
        let lc_nonzero = [a.leading_coeff(), b.leading_coeff()]
            .into_iter()
            .flatten()
            .all(|e| map_element_gf(e, &gf).is_some_and(|lc| !gf.is_zero(&lc)));
        if !lc_nonzero {
            continue;
        }
        let (Some(a_fp), Some(b_fp)) = (map_poly_gf(a, &gf), map_poly_gf(b, &gf)) else {
            continue; // p divides a coefficient denominator
        };
        let g_fp = monic_over(&a_fp.gcd(&b_fp));
        let deg_g = g_fp.degree().unwrap_or(0);
        // Residue table of this modular GCD, padded to [deg+1][d].
        let mut table = vec![vec![Integer::from(0); d]; deg_g + 1];
        for (i, e) in g_fp.coeffs().iter().enumerate() {
            for (j, c) in e.coeffs().iter().enumerate() {
                table[i][j] = Integer::from(c.value().clone());
            }
        }
        let p_as_int = p_int;
        let dominated = match &residues {
            None => true,
            Some(prev) => deg_g + 1 < prev.len(),
        };
        let superseded = matches!(&residues, Some(prev) if deg_g + 1 > prev.len());
        if superseded {
            continue; // unlucky prime: modular degree too large
        }
        if dominated {
            residues = Some(table);
            modulus = p_as_int;
        } else {
            let prev = residues.take().expect("CRT state present");
            let mut merged = prev;
            for (i, row) in merged.iter_mut().enumerate() {
                for (j, r) in row.iter_mut().enumerate() {
                    // All coefficients merge against the same old modulus.
                    let (res, _) =
                        crt(r, &modulus, &table[i][j], &p_as_int).expect("coprime moduli");
                    *r = res;
                }
            }
            modulus = IntegerDomain.mul(&modulus, &p_as_int);
            residues = Some(merged);
        }
        // Try rational reconstruction + trial division.
        let table = residues.clone().expect("CRT state present");
        let mut coeffs: Vec<AnfElement> = Vec::with_capacity(table.len());
        let mut ok = true;
        for row in &table {
            let mut alpha_coeffs = Vec::with_capacity(d);
            for r in row {
                match rational_reconstruction(r, &modulus) {
                    Some((n, dd)) => {
                        alpha_coeffs.push(Rational::from_bigints(n.to_bigint(), dd.to_bigint()));
                    }
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                break;
            }
            coeffs.push(field.element(alpha_coeffs));
        }
        if !ok {
            continue;
        }
        let candidate = UP::<AlgebraicNumberField>::from_coeffs(field.clone(), coeffs);
        if candidate.degree() != Some(deg_g) {
            continue; // reconstruction produced a spurious leading zero
        }
        let divides_a = a.div_rem(&candidate).is_some_and(|(_, r)| r.is_zero());
        let divides_b = b.div_rem(&candidate).is_some_and(|(_, r)| r.is_zero());
        if divides_a && divides_b {
            return candidate;
        }
    }
    // Correctness fallback (never reached in practice): dense Euclidean
    // GCD over the number field, made monic.
    monic_over(&a.gcd(b))
}

/// Embed a rational polynomial into $K[x]$ coefficient-wise.
fn embed_rational_poly(
    field: &AlgebraicNumberField,
    f: &UP<RationalDomain>,
) -> UP<AlgebraicNumberField> {
    UP::<AlgebraicNumberField>::from_coeffs(
        field.clone(),
        f.coeffs()
            .iter()
            .map(|c| field.from_base(c.clone()))
            .collect(),
    )
}

/// Trager factorization of a square-free polynomial over $K$:
/// norm → factor over $\mathbb{Q}$ → GCD back over $K$.
///
/// Returns the monic irreducible factors of `f` (each with multiplicity 1).
/// Falls back to returning `f` itself if no Trager shift is found within
/// [`MAX_TRAGER_SHIFTS`] attempts (theoretically impossible).
fn factor_square_free_anf(
    field: &AlgebraicNumberField,
    f: &UP<AlgebraicNumberField>,
) -> Vec<UP<AlgebraicNumberField>> {
    if f.degree().unwrap_or(0) <= 1 {
        return vec![monic_over(f)];
    }
    let Some((s, g, r)) = norm_with_shift(field, f) else {
        return vec![monic_over(f)];
    };
    let rational_factors = factor_square_free_rationals(&r);
    if rational_factors.len() <= 1 {
        // The norm is irreducible over ℚ, so f is irreducible over K.
        return vec![monic_over(f)];
    }
    // Undoing the shift: the norm used x ↦ x − s·α, so factors of
    // g(x) = f(x − sα) are substituted back with x ↦ x + s·α.
    let back_shift = field.mul(&field.alpha(), &field.from_base(Rational::new(s as i64, 1)));
    let mut remaining = g;
    let mut out = Vec::new();
    for n_i in &rational_factors {
        if remaining.degree().unwrap_or(0) == 0 {
            break;
        }
        let n_k = embed_rational_poly(field, n_i);
        let h = gcd_anf(field, &n_k, &remaining);
        if h.degree().unwrap_or(0) == 0 {
            continue; // no common factor (should not happen for good shifts)
        }
        remaining = remaining
            .div_rem(&h)
            .map(|(q, _)| q)
            .unwrap_or_else(|| remaining.clone());
        out.push(monic_over(&compose_linear(&h, &back_shift)));
    }
    out
}

/// Yun square-free factorization over the number field, using the modular
/// [`gcd_anf`] instead of the generic pseudo-remainder GCD (whose
/// coefficients explode over `BigRational`-valued field elements).
fn square_free_anf(
    field: &AlgebraicNumberField,
    f: &UP<AlgebraicNumberField>,
) -> Vec<(UP<AlgebraicNumberField>, usize)> {
    let mut out = Vec::new();
    if f.is_zero() {
        return out;
    }
    let f = monic_over(f);
    let df = f.derivative();
    if df.is_zero() {
        out.push((f, 1));
        return out;
    }
    let mut g = gcd_anf(field, &f, &df);
    let mut w = match f.div_rem(&g) {
        Some((q, _)) => q,
        None => return out,
    };
    let mut k = 1;
    while !w.is_one() {
        let h = gcd_anf(field, &w, &g);
        if let Some((z, _)) = w.div_rem(&h)
            && z.degree().unwrap_or(0) > 0
        {
            out.push((monic_over(&z), k));
        }
        match g.div_rem(&h) {
            Some((q, _)) => g = q,
            None => break,
        }
        w = h;
        k += 1;
    }
    out
}

/// Factor a univariate polynomial over the algebraic number field $K =
/// \mathbb{Q}(\alpha)$ into monic irreducible factors with multiplicities
/// (Trager's algorithm).
///
/// The product of `factor^multiplicity` equals `self` up to a unit of $K$
/// (the leading coefficient of `self`).
fn factor_anf(
    field: &AlgebraicNumberField,
    f: &UP<AlgebraicNumberField>,
) -> Vec<(UP<AlgebraicNumberField>, usize)> {
    if f.is_zero() || f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    // Fast path: all coefficients are rational constants. Factor over ℚ
    // first (much cheaper than a norm of degree deg(f)·[K:ℚ], and a
    // rational input's norm is always a perfect power, forcing a shift),
    // then split each ℚ-irreducible factor over K individually.
    if f.coeffs().iter().all(|c| c.coeffs().len() <= 1) {
        let f_q = UP::<RationalDomain>::from_coeffs(
            RationalDomain,
            f.coeffs()
                .iter()
                .map(|c| {
                    c.coeffs()
                        .first()
                        .cloned()
                        .unwrap_or_else(|| Rational::new(0, 1))
                })
                .collect(),
        );
        let primitive = clear_denominators(&f_q).primitive_part();
        let mut out = Vec::new();
        for (h, mult) in factor_primitive(&primitive) {
            let h_k = UP::<AlgebraicNumberField>::from_coeffs(
                field.clone(),
                h.coeffs()
                    .iter()
                    .map(|c| field.from_base(Rational::from_integer(c.clone())))
                    .collect(),
            );
            for g in factor_square_free_anf(field, &h_k) {
                out.push((g, mult));
            }
        }
        return out;
    }
    let mut out = Vec::new();
    for (g, mult) in square_free_anf(field, f) {
        for h in factor_square_free_anf(field, &g) {
            out.push((h, mult));
        }
    }
    out
}

impl DenseUnivariatePolynomial<AlgebraicNumberField> {
    /// Factor this polynomial over its algebraic number field
    /// $\mathbb{Q}(\alpha)$ into monic irreducible factors with
    /// multiplicities, using Trager's algorithm: square-free
    /// factorization, then for each component a norm down to
    /// $\mathbb{Q}[x]$ (evaluation–interpolation of the resultant), a
    /// rational factorization, and GCDs back over the number field
    /// (modular method over $\mathrm{GF}(p^d)$ with CRT and rational
    /// reconstruction).
    ///
    /// The product of `factor^multiplicity` equals `self` up to a unit
    /// (the leading coefficient).
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// // ℚ(√2): minimal polynomial α² − 2.
    /// let field = AlgebraicNumberField::new(
    ///     RationalDomain,
    ///     vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
    /// );
    /// // x² − 2 = (x − √2)(x + √2) splits over ℚ(√2).
    /// let f = DenseUnivariatePolynomial::from_coeffs(
    ///     field.clone(),
    ///     vec![
    ///         field.from_base(Rational::new(-2, 1)),
    ///         field.zero(),
    ///         field.one(),
    ///     ],
    /// );
    /// let factors = f.factor();
    /// assert_eq!(factors.len(), 2);
    /// assert!(factors.iter().all(|(g, m)| *m == 1 && g.degree() == Some(1)));
    /// ```
    pub fn factor(&self) -> crate::factor::Factors<AlgebraicNumberField> {
        factor_anf(self.domain(), self)
    }
}

/// Trager's shift: find $s \ge 0$ such that the norm of
/// $g(x) = f(x - s\alpha)$ is square-free over $\mathbb{Q}$.
///
/// Returns `(s, g, norm)`. Such a shift always exists (the bad shifts are
/// finite in number); `None` is returned only after [`MAX_TRAGER_SHIFTS`]
/// attempts.
pub(crate) fn norm_with_shift(
    field: &AlgebraicNumberField,
    f: &UP<AlgebraicNumberField>,
) -> Option<(u64, UP<AlgebraicNumberField>, UP<RationalDomain>)> {
    let alpha = field.alpha();
    for s in 0..MAX_TRAGER_SHIFTS {
        let shift = field.mul(&alpha, &field.from_base(Rational::new(-(s as i64), 1)));
        let g = compose_linear(f, &shift);
        let r = norm_eval_interp(field, &g);
        if is_square_free_rational(&r) {
            return Some((s, g, r));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn q(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    /// ℚ(√2): minimal polynomial α² − 2.
    pub(crate) fn sqrt2_field() -> AlgebraicNumberField {
        AlgebraicNumberField::new(RationalDomain, vec![q(-2, 1), q(0, 1), q(1, 1)])
    }

    /// Build an ANF polynomial from per-coefficient α-polynomials.
    pub(crate) fn anf_poly(
        field: &AlgebraicNumberField,
        coeffs: Vec<Vec<Rational>>,
    ) -> UP<AlgebraicNumberField> {
        UP::<AlgebraicNumberField>::from_coeffs(
            field.clone(),
            coeffs.into_iter().map(|c| field.element(c)).collect(),
        )
    }

    /// Multiply by the inverse of the leading coefficient.
    fn monic(p: &UP<RationalDomain>) -> UP<RationalDomain> {
        let lc = p.leading_coeff().cloned().expect("nonzero polynomial");
        p.mul_scalar(&RationalDomain.inv(&lc).expect("nonzero lc"))
    }

    #[test]
    fn compose_linear_over_integers() {
        let f = UP::<IntegerDomain>::from_coeffs(
            IntegerDomain,
            vec![Integer::from(0), Integer::from(0), Integer::from(1)],
        ); // x²
        let g = compose_linear(&f, &Integer::from(1)); // (x + 1)²
        assert_eq!(
            g.coeffs(),
            &[Integer::from(1), Integer::from(2), Integer::from(1)]
        );
    }

    #[test]
    fn compose_linear_over_anf() {
        let field = sqrt2_field();
        let alpha = field.alpha();
        let f = anf_poly(&field, vec![vec![], vec![], vec![q(1, 1)]]); // x²
        let shift = field.neg(&alpha);
        let g = compose_linear(&f, &shift); // (x − α)² = x² − 2αx + 2
        let expected = anf_poly(
            &field,
            vec![vec![q(2, 1)], vec![q(0, 1), q(-2, 1)], vec![q(1, 1)]],
        );
        assert_eq!(g, expected);
    }

    #[test]
    fn factor_over_rationals_basic() {
        // (1/2)x² − 2 = (1/2)(x − 2)(x + 2)
        let f = UP::<RationalDomain>::from_coeffs(RationalDomain, vec![q(-2, 1), q(0, 1), q(1, 2)]);
        let factors = factor_square_free_rationals(&f);
        assert_eq!(factors.len(), 2);
        let monics: Vec<UP<RationalDomain>> = factors.iter().map(monic).collect();
        let x_minus_2 = UP::<RationalDomain>::from_coeffs(RationalDomain, vec![q(-2, 1), q(1, 1)]);
        let x_plus_2 = UP::<RationalDomain>::from_coeffs(RationalDomain, vec![q(2, 1), q(1, 1)]);
        assert!(
            monics.contains(&x_minus_2),
            "x − 2 expected, got {monics:?}"
        );
        assert!(monics.contains(&x_plus_2), "x + 2 expected, got {monics:?}");
    }

    #[test]
    fn norm_of_x_minus_alpha() {
        let field = sqrt2_field();
        // f = x − α; norm = (x − √2)(x + √2) = x² − 2, square-free at s = 0.
        let f = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(1, 1)]]);
        let (s, _g, r) = norm_with_shift(&field, &f).expect("shift exists");
        assert_eq!(s, 0);
        assert_eq!(monic(&r).coeffs(), &[q(-2, 1), q(0, 1), q(1, 1)]);
    }

    #[test]
    fn norm_shift_for_x2_minus_2() {
        let field = sqrt2_field();
        // f = x² − 2: norms at s = 0 ((x²−2)²) and s = 1 (x²(x²−8)) are not
        // square-free; at s = 2 the norm is x⁴ − 20x² + 36.
        let f = anf_poly(&field, vec![vec![q(-2, 1)], vec![], vec![q(1, 1)]]);
        let (s, _g, r) = norm_with_shift(&field, &f).expect("shift exists");
        assert_eq!(s, 2);
        assert_eq!(
            monic(&r).coeffs(),
            &[q(36, 1), q(0, 1), q(-20, 1), q(0, 1), q(1, 1)]
        );
    }

    #[test]
    fn gcd_anf_sqrt2_common_factor() {
        let field = sqrt2_field();
        // a = (x − α)(x + 1) = x² + (1 − α)x − α
        let a = anf_poly(
            &field,
            vec![
                vec![q(0, 1), q(-1, 1)],
                vec![q(1, 1), q(-1, 1)],
                vec![q(1, 1)],
            ],
        );
        // b = (x − α)(x − 1) = x² − (1 + α)x + α
        let b = anf_poly(
            &field,
            vec![
                vec![q(0, 1), q(1, 1)],
                vec![q(-1, 1), q(-1, 1)],
                vec![q(1, 1)],
            ],
        );
        let g = gcd_anf(&field, &a, &b);
        let expected = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(1, 1)]]);
        assert_eq!(g, expected, "gcd must be x − α");
    }

    #[test]
    fn gcd_anf_coprime_and_zero() {
        let field = sqrt2_field();
        let a = anf_poly(&field, vec![vec![q(1, 1)], vec![q(1, 1)]]); // x + 1
        let b = anf_poly(&field, vec![vec![q(2, 1)], vec![q(1, 1)]]); // x + 2
        let g = gcd_anf(&field, &a, &b);
        assert!(g.degree() == Some(0) && g.is_one(), "coprime → gcd 1");
        // gcd(0, b) = monic(b).
        let zero = UP::<AlgebraicNumberField>::from_coeffs(field.clone(), vec![]);
        let g = gcd_anf(&field, &zero, &b);
        assert_eq!(g, b);
    }

    #[test]
    fn gcd_anf_number_field_symbolica_mirror() {
        // Mirror of Symbolica's `gcd_number_field` test:
        // min_poly a³ + 3a² − 46a + 1.
        let field =
            AlgebraicNumberField::new(RationalDomain, vec![q(1, 1), q(-46, 1), q(3, 1), q(1, 1)]);
        // a = x³ − 2x² + (−2a² + 8a + 2)x − a² + 11a − 1
        let a = anf_poly(
            &field,
            vec![
                vec![q(-1, 1), q(11, 1), q(-1, 1)],
                vec![q(2, 1), q(8, 1), q(-2, 1)],
                vec![q(-2, 1)],
                vec![q(1, 1)],
            ],
        );
        // b = x³ − 2x² − x + 1
        let b = anf_poly(
            &field,
            vec![vec![q(1, 1)], vec![q(-1, 1)], vec![q(-2, 1)], vec![q(1, 1)]],
        );
        let g = gcd_anf(&field, &a, &b);
        // Expected: x − 50/91 − (23/91)a − (1/91)a²
        let expected = anf_poly(
            &field,
            vec![vec![q(-50, 91), q(-23, 91), q(-1, 91)], vec![q(1, 1)]],
        );
        assert_eq!(g, expected);
    }

    /// Check `f == lc(f) · ∏ factor^mult` (factors are monic).
    pub(crate) fn reconstructs(
        field: &AlgebraicNumberField,
        f: &UP<AlgebraicNumberField>,
        factors: &[(UP<AlgebraicNumberField>, usize)],
    ) -> bool {
        let mut acc = UP::<AlgebraicNumberField>::from_coeffs(field.clone(), vec![field.one()]);
        for (h, e) in factors {
            for _ in 0..*e {
                acc = acc.mul(h);
            }
        }
        let lc = f.leading_coeff().cloned().expect("nonzero input");
        acc.mul_scalar(&lc) == *f
    }

    #[test]
    fn factor_x2_minus_2_over_sqrt2() {
        let field = sqrt2_field();
        let f = anf_poly(&field, vec![vec![q(-2, 1)], vec![], vec![q(1, 1)]]);
        let factors = f.factor();
        assert_eq!(factors.len(), 2);
        assert!(
            factors
                .iter()
                .all(|(g, m)| *m == 1 && g.degree() == Some(1))
        );
        let x_minus_a = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(1, 1)]]);
        let x_plus_a = anf_poly(&field, vec![vec![q(0, 1), q(1, 1)], vec![q(1, 1)]]);
        let monics: Vec<_> = factors.iter().map(|(g, _)| g.clone()).collect();
        assert!(monics.contains(&x_minus_a));
        assert!(monics.contains(&x_plus_a));
        assert!(reconstructs(&field, &f, &factors));
    }

    #[test]
    fn factor_x2_plus_1_over_gaussian() {
        // ℚ(i): α² + 1; x² + 1 = (x − i)(x + i).
        let field = AlgebraicNumberField::new(RationalDomain, vec![q(1, 1), q(0, 1), q(1, 1)]);
        let f = anf_poly(&field, vec![vec![q(1, 1)], vec![], vec![q(1, 1)]]);
        let factors = f.factor();
        assert_eq!(factors.len(), 2);
        assert!(reconstructs(&field, &f, &factors));
    }

    #[test]
    fn factor_x3_minus_2_over_cbrt2() {
        // ℚ(∛2): α³ − 2; x³ − 2 = (x − α)(x² + αx + α²).
        let field =
            AlgebraicNumberField::new(RationalDomain, vec![q(-2, 1), q(0, 1), q(0, 1), q(1, 1)]);
        let f = anf_poly(&field, vec![vec![q(-2, 1)], vec![], vec![], vec![q(1, 1)]]);
        let factors = f.factor();
        assert_eq!(factors.len(), 2);
        let linear = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(1, 1)]]);
        let quadratic = anf_poly(
            &field,
            vec![
                vec![q(0, 1), q(0, 1), q(1, 1)],
                vec![q(0, 1), q(1, 1)],
                vec![q(1, 1)],
            ],
        );
        let monics: Vec<_> = factors.iter().map(|(g, _)| g.clone()).collect();
        assert!(monics.contains(&linear), "x − α expected");
        assert!(monics.contains(&quadratic), "x² + αx + α² expected");
        assert!(reconstructs(&field, &f, &factors));
    }

    #[test]
    fn factor_repeated_over_sqrt2() {
        let field = sqrt2_field();
        // (x − α)²(x + α) — exercises the Yun square-free stage.
        let x_minus_a = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(1, 1)]]);
        let x_plus_a = anf_poly(&field, vec![vec![q(0, 1), q(1, 1)], vec![q(1, 1)]]);
        let f = x_minus_a.mul(&x_minus_a).mul(&x_plus_a);
        let factors = f.factor();
        assert_eq!(factors.len(), 2);
        let mult_of =
            |g: &UP<AlgebraicNumberField>| factors.iter().find(|(h, _)| h == g).map(|(_, m)| *m);
        assert_eq!(mult_of(&x_minus_a), Some(2));
        assert_eq!(mult_of(&x_plus_a), Some(1));
        assert!(reconstructs(&field, &f, &factors));
    }

    #[test]
    fn factor_symbolica_quartic_mirror() {
        // Mirror of Symbolica's `algebraic_extension` test:
        // ℚ(⁴√3) with m = a⁴ − 3, and
        // f = z⁴ + z³ + (2 + a − a²)z² + (1 + a² − 2a³)z − 2
        //   = (z² + (1 − a)z + (1 − a²))(z² + az + (1 + a²)).
        let field = AlgebraicNumberField::new(
            RationalDomain,
            vec![q(-3, 1), q(0, 1), q(0, 1), q(0, 1), q(1, 1)],
        );
        let f = anf_poly(
            &field,
            vec![
                vec![q(-2, 1)],
                vec![q(1, 1), q(0, 1), q(1, 1), q(-2, 1)],
                vec![q(2, 1), q(1, 1), q(-1, 1)],
                vec![q(1, 1)],
                vec![q(1, 1)],
            ],
        );
        let factors = f.factor();
        assert_eq!(factors.len(), 2);
        let f1 = anf_poly(
            &field,
            vec![
                vec![q(1, 1), q(0, 1), q(-1, 1)],
                vec![q(1, 1), q(-1, 1)],
                vec![q(1, 1)],
            ],
        );
        let f2 = anf_poly(
            &field,
            vec![
                vec![q(1, 1), q(0, 1), q(1, 1)],
                vec![q(0, 1), q(1, 1)],
                vec![q(1, 1)],
            ],
        );
        let monics: Vec<_> = factors.iter().map(|(g, _)| g.clone()).collect();
        assert!(monics.contains(&f1), "z² + (1−a)z + (1−a²) expected");
        assert!(monics.contains(&f2), "z² + az + (1+a²) expected");
        assert!(reconstructs(&field, &f, &factors));
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;
    use tests::{anf_poly, q, reconstructs, sqrt2_field};

    /// A random small ANF element c₀ + c₁·α over ℚ(√2).
    fn any_sqrt2_element() -> impl Strategy<Value = Vec<Rational>> {
        (-3i64..=3, -3i64..=3).prop_map(|(c0, c1)| vec![q(c0, 1), q(c1, 1)])
    }

    /// A random monic linear or quadratic polynomial over ℚ(√2).
    fn any_small_factor(
        field: &AlgebraicNumberField,
    ) -> impl Strategy<Value = UP<AlgebraicNumberField>> {
        let linear = any_sqrt2_element().prop_map(|c| vec![c, vec![q(1, 1)]]);
        let quadratic = (any_sqrt2_element(), any_sqrt2_element())
            .prop_map(|(c0, c1)| vec![c0, c1, vec![q(1, 1)]]);
        let field = field.clone();
        prop::strategy::Union::new([linear.boxed(), quadratic.boxed()])
            .prop_map(move |coeffs| anf_poly(&field, coeffs))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        /// Factoring a product of 2–3 small monic factors over ℚ(√2) must
        /// reconstruct it.
        ///
        /// Marked `ignore` because ANF factorization is slow enough that a
        /// property-based sweep does not fit the unit-test budget; run
        /// manually or via the audit report.
        #[test]
        #[ignore = "slow ANF factorization proptest: run manually"]
        fn factor_product_roundtrip_sqrt2(
            factors in prop::collection::vec(
                any_small_factor(&sqrt2_field()),
                2..=3,
            ),
        ) {
            let field = sqrt2_field();
            let mut f = UP::<AlgebraicNumberField>::from_coeffs(
                field.clone(),
                vec![field.one()],
            );
            for g in &factors {
                f = f.mul(g);
            }
            let factored = f.factor();
            prop_assert!(
                reconstructs(&field, &f, &factored),
                "factorization does not reconstruct input: {:?}",
                factored
            );
        }
    }
}
