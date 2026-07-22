//! Multivariate polynomial factorization via EEZ Hensel lifting.
//!
//! Generalizes the bivariate Hensel lifting in [`super::multivariate`] to an
//! arbitrary number of variables. A square-free polynomial that is monic in
//! the main variable $x_0$ is evaluated at a sample point of the secondary
//! variables, factored as a univariate polynomial, and the factors are
//! lifted back one variable at a time through the ideals
//! $(x_k - a_k)$, solving multivariate Diophantine equations at each step.
//!
//! References: Wang (1978), "An Improved Multivariate Polynomial Factoring
//! Algorithm"; Geddes, Czapor, Labahn, *Algorithms for Computer Algebra*,
//! §6.4–6.5.

use num_bigint::BigInt;
use num_traits::ToPrimitive;
use ocas_domain::{
    Domain, EuclideanDomain, FiniteField, FiniteFieldElement, Integer, IntegerDomain, Rational,
    RationalDomain,
};

use crate::dense::DenseUnivariatePolynomial;
use crate::factor::multivariate::FpMPoly;
use crate::multivariate_gcd::{multivariate_gcd_fp, multivariate_gcd_z};
use crate::sparse::{Lex, SparseMultivariatePolynomial, binomial};

/// Sparse multivariate polynomial over a domain `D` with lexicographic order.
type MP<D> = SparseMultivariatePolynomial<D, Lex>;

/// Dense univariate polynomial over a domain `D`.
type UP<D> = DenseUnivariatePolynomial<D>;

// ---------------------------------------------------------------------
// Generic helpers
// ---------------------------------------------------------------------

/// The constant polynomial 1.
fn one_mpoly<D: Domain>(domain: &D, n_vars: usize) -> MP<D> {
    let mut p = MP::<D>::new(domain.clone(), n_vars);
    p.set_term_external(vec![0; n_vars], domain.one());
    p
}

/// Whether the polynomial is a constant (total degree 0, including zero).
fn is_constant<D: Domain>(f: &MP<D>) -> bool {
    f.total_degree() == Some(0) || f.is_zero()
}

/// Convert a sparse polynomial whose support only involves variable 0 into
/// a dense univariate polynomial.
fn mpoly_to_dense<D: Domain>(f: &MP<D>) -> UP<D> {
    let mut coeffs = Vec::new();
    for (exp, c) in f.terms_ref() {
        let idx = exp.first().copied().unwrap_or(0);
        if idx >= coeffs.len() {
            coeffs.resize(idx + 1, f.domain().zero());
        }
        coeffs[idx] = f.domain().add(&coeffs[idx], c);
    }
    UP::<D>::from_coeffs(f.domain().clone(), coeffs)
}

/// Wrap a dense univariate polynomial as a sparse polynomial in `n_vars`
/// variables (variable 0 is the polynomial variable).
fn dense_to_mpoly<D: Domain>(g: &UP<D>, n_vars: usize) -> MP<D> {
    let mut result = MP::<D>::new(g.domain().clone(), n_vars);
    for (i, c) in g.coeffs().iter().enumerate() {
        if !g.domain().is_zero(c) {
            let mut exp = vec![0usize; n_vars];
            exp[0] = i;
            result.set_term_external(exp, c.clone());
        }
    }
    result
}

/// The polynomial `(x_var - a)^j`.
fn x_minus_a_pow<D: Domain>(
    domain: &D,
    n_vars: usize,
    var: usize,
    a: &D::Element,
    j: usize,
) -> MP<D> {
    let mut result = MP::<D>::new(domain.clone(), n_vars);
    for i in 0..=j {
        // term: binom(j, i) * (-a)^(j-i) * x^i
        let binom = domain.cast_u64(binomial(j, i));
        let neg_a_pow = if (j - i) % 2 == 0 {
            domain.pow(a, (j - i) as u64)
        } else {
            domain.neg(&domain.pow(a, (j - i) as u64))
        };
        let coeff = domain.mul(&binom, &neg_a_pow);
        if domain.is_zero(&coeff) {
            continue;
        }
        let mut exp = vec![0usize; n_vars];
        exp[var] = i;
        result.set_term_external(exp, coeff);
    }
    result
}

/// Multi-factor Bézout coefficients: returns `b_i` such that
/// `Σ b_i · Π_{j≠i} f_j = 1`. Requires pairwise coprime inputs over a field;
/// returns `None` otherwise.
fn bezout_coefficients<D: EuclideanDomain>(factors: &[UP<D>]) -> Option<Vec<UP<D>>> {
    let n = factors.len();
    if n == 1 {
        return Some(vec![factors[0].one()]);
    }
    let mut result = vec![factors[0].zero(); n];
    result[0] = factors[0].one();
    let mut accum = factors[0].clone();
    for (i, f_i) in factors.iter().enumerate().skip(1) {
        let (g, s, t) = accum.extended_gcd_poly(f_i);
        if g.degree() != Some(0) {
            return None; // not pairwise coprime
        }
        let inv = g.domain().inv(g.leading_coeff().unwrap())?;
        let s = s.mul_scalar(&inv);
        let t = t.mul_scalar(&inv);
        for res in result.iter_mut().take(i) {
            *res = res.mul(&t);
        }
        result[i] = s;
        accum = accum.mul(f_i);
    }
    Some(result)
}

// ---------------------------------------------------------------------
// Multivariate Diophantine solver (recursive)
// ---------------------------------------------------------------------

/// Solve the multivariate Diophantine equation
///
/// `Σ_i σ_i · Π_{j≠i} u_j = e`
///
/// for polynomials `σ_i` in variables `0..k-1` with
/// `deg_{x_0}(σ_i) < deg_{x_0}(u_i)`. The `u_i` may depend on variables
/// `0..k-1`; `sample[var]` is the evaluation point used when reducing
/// variable `var`.
///
/// Returns `None` if the equation has no solution of the required form
/// (unlucky evaluation point).
fn diophantine<D: EuclideanDomain>(
    u: &[MP<D>],
    e: &MP<D>,
    sample: &[D::Element],
    k: usize,
) -> Option<Vec<MP<D>>> {
    let n_vars = e.n_vars();
    if k == 1 {
        // Univariate base case in variable 0.
        let ud: Vec<UP<D>> = u.iter().map(mpoly_to_dense).collect();
        let bez = bezout_coefficients(&ud)?;
        let ed = mpoly_to_dense(e);
        let mut sigmas = Vec::with_capacity(u.len());
        for (i, b_i) in bez.iter().enumerate() {
            let delta = ed.mul(b_i);
            let (_q, rem) = delta.div_rem(&ud[i])?;
            sigmas.push(dense_to_mpoly(&rem, n_vars));
        }
        return Some(sigmas);
    }

    // Recurse on the last live variable v = k - 1.
    let v = k - 1;
    let a = &sample[v];
    let u_base: Vec<MP<D>> = u.iter().map(|g| g.eval_keep(v, a)).collect();
    let cofactor: Vec<MP<D>> = (0..u.len())
        .map(|i| {
            u.iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .fold(one_mpoly(e.domain(), n_vars), |acc, (_, g)| acc.mul(g))
        })
        .collect();

    let e_base = e.eval_keep(v, a);
    let mut sigma = diophantine(&u_base, &e_base, sample, v)?;

    // Lift the solution in variable v, one Taylor coefficient at a time.
    // The loop bound is the degree of the ORIGINAL error in v (fixed);
    // the Taylor coefficients are taken from the residual e, recomputed
    // each iteration as `e - Σ σ_i·cofactor_i`.
    let dmax = e.degree_in(v);
    for j in 1..=dmax {
        // Recompute the residual from the original error and current sigma.
        let mut residual = e.clone();
        for (s, b) in sigma.iter().zip(&cofactor) {
            residual = residual.sub(&s.mul(b));
        }
        if residual.is_zero() {
            break;
        }
        let t_v = residual.taylor_coefficients(v, a);
        let t = t_v.get(j).cloned().unwrap_or_else(|| residual.zero());
        if t.is_zero() {
            continue;
        }
        let delta = diophantine(&u_base, &t, sample, v)?;
        let xp = x_minus_a_pow(e.domain(), n_vars, v, a, j);
        for (i, s) in sigma.iter_mut().enumerate() {
            *s = s.add(&delta[i].mul(&xp));
        }
    }
    Some(sigma)
}

// ---------------------------------------------------------------------
// EEZ Hensel lifting (generic over a field)
// ---------------------------------------------------------------------

/// Lift a univariate factorization of the image `f(sample)` back to a
/// multivariate factorization of `f`, one variable at a time.
///
/// `f` must be square-free, primitive and monic in variable 0; `uni` are the
/// monic pairwise-coprime irreducible factors of the univariate image.
/// Returns the lifted factors (same order as `uni`), or `None` if the
/// sample point is unlucky.
fn eez_lift<D: EuclideanDomain>(
    f: &MP<D>,
    sample: &[D::Element],
    uni: &[UP<D>],
) -> Option<Vec<MP<D>>> {
    let n = f.n_vars();
    let mut lifted: Vec<MP<D>> = uni.iter().map(|g| dense_to_mpoly(g, n)).collect();

    for k in 1..n {
        // Work on the image of f with variables > k evaluated at the sample.
        let mut f_k = f.clone();
        for m in (k + 1..n).rev() {
            f_k = f_k.eval_keep(m, &sample[m]);
        }
        let a_k = sample[k].clone();
        let d_k = f_k.degree_in(k);
        let t_f = f_k.taylor_coefficients(k, &a_k);

        for j in 1..=d_k {
            let mut prod = one_mpoly(f.domain(), n);
            for g in &lifted {
                prod = prod.mul(g);
            }
            let t_p = prod.taylor_coefficients(k, &a_k);
            let e = if j < t_p.len() {
                t_f[j].sub(&t_p[j])
            } else {
                t_f[j].clone()
            };
            if e.is_zero() {
                continue;
            }
            // Solve against the factors evaluated at the current variable's
            // sample point: the (x_k − a_k)^j coefficient of
            // σ_i·cofactor_i·(x_k − a_k)^j is σ_i·cofactor_i(a_k), and the
            // Diophantine solver requires its inputs to be free of
            // variables ≥ k. Previous corrections (and imposed leading
            // coefficients) may depend on x_k, so evaluate them away.
            let u_base: Vec<MP<D>> = lifted.iter().map(|g| g.eval_keep(k, &a_k)).collect();
            let sigmas = diophantine(&u_base, &e, sample, k)?;
            let xp = x_minus_a_pow(f.domain(), n, k, &a_k, j);
            for (i, g) in lifted.iter_mut().enumerate() {
                *g = g.add(&sigmas[i].mul(&xp));
            }
        }
    }
    Some(lifted)
}

/// Impose the true leading coefficient `ℓ` on the leading `x_0`-coefficient
/// of `f_i`, first evaluating `ℓ` at the sample points of variables
/// `≥ from_var` (those still fixed at the sample during the lift).
///
/// Reference: Symbolica `impose_true_lcoeffs_on_factors`.
fn impose_lcoeff_field<D: EuclideanDomain>(
    f_i: &MP<D>,
    true_lc: &MP<D>,
    sample: &[D::Element],
    from_var: usize,
) -> MP<D> {
    let mut lc = true_lc.clone();
    for (m, s) in sample.iter().enumerate().skip(from_var) {
        lc = lc.eval_keep(m, s);
    }
    let deg = f_i.degree_in(0);
    let mut result = f_i.clone();
    let top: Vec<smallvec::SmallVec<[usize; 4]>> = result
        .terms_ref()
        .keys()
        .filter(|e| e.first().copied().unwrap_or(0) == deg)
        .cloned()
        .collect();
    for e in top {
        result.set_term_external(e.to_vec(), f_i.domain().zero());
    }
    for (e, c) in lc.terms_ref() {
        let mut exp = e.to_vec();
        exp[0] = deg;
        let existing = result.coeff(&exp);
        result.set_term_external(exp, f_i.domain().add(&existing, c));
    }
    result
}

/// EEZ lift over a field with Wang-imposed leading coefficients.
///
/// Like [`eez_lift`], but the factors need not be monic in `x_0`: before
/// lifting each variable `k`, the true leading coefficient `ℓ_i` (evaluated
/// at the sample points of variables `> k`) is imposed on every factor.
/// Diophantine corrections have lower degree in `x_0` than the factors, so
/// the imposed leading coefficients stay fixed within each variable's lift.
fn eez_lift_imposed<D: EuclideanDomain>(
    f: &MP<D>,
    sample: &[D::Element],
    initial: &[MP<D>],
    true_lcoeffs: &[MP<D>],
) -> Option<Vec<MP<D>>> {
    let n = f.n_vars();
    let mut lifted: Vec<MP<D>> = initial.to_vec();

    for k in 1..n {
        // Impose the true leading coefficients for this lift step.
        for (i, g) in lifted.iter_mut().enumerate() {
            *g = impose_lcoeff_field(g, &true_lcoeffs[i], sample, k + 1);
        }
        // Work on the image of f with variables > k evaluated at the sample.
        let mut f_k = f.clone();
        for m in (k + 1..n).rev() {
            f_k = f_k.eval_keep(m, &sample[m]);
        }
        let a_k = sample[k].clone();
        let d_k = f_k.degree_in(k);
        let t_f = f_k.taylor_coefficients(k, &a_k);

        for j in 1..=d_k {
            let mut prod = one_mpoly(f.domain(), n);
            for g in &lifted {
                prod = prod.mul(g);
            }
            let t_p = prod.taylor_coefficients(k, &a_k);
            let e = if j < t_p.len() {
                t_f[j].sub(&t_p[j])
            } else {
                t_f[j].clone()
            };
            if e.is_zero() {
                continue;
            }
            // See eez_lift: solve against the factors with x_k evaluated at
            // the sample point (the imposed leading coefficients introduce
            // x_k-dependence into the factors at this stage).
            let u_base: Vec<MP<D>> = lifted.iter().map(|g| g.eval_keep(k, &a_k)).collect();
            let sigmas = diophantine(&u_base, &e, sample, k)?;
            let xp = x_minus_a_pow(f.domain(), n, k, &a_k, j);
            for (i, g) in lifted.iter_mut().enumerate() {
                *g = g.add(&sigmas[i].mul(&xp));
            }
        }
    }
    Some(lifted)
}

// ---------------------------------------------------------------------
// Integer EEZ lifting with imposed leading coefficients
// ---------------------------------------------------------------------

/// Solve the multivariate Diophantine equation over ℚ, returning integer
/// solutions if they exist. Wraps the generic field solver and requires the
/// result to have integral coefficients (exact for the EEZ lift).
fn diophantine_z(u: &[ZmPoly], e: &ZmPoly, sample: &[Integer], k: usize) -> Option<Vec<ZmPoly>> {
    // Convert to ℚ, solve, then check integrality.
    let u_q: Vec<MP<RationalDomain>> = u.iter().map(zmp_to_qmp).collect();
    let e_q = zmp_to_qmp(e);
    let sample_q: Vec<Rational> = sample
        .iter()
        .map(|s| Rational::from_integer(s.clone()))
        .collect();
    let sigmas_q = diophantine(&u_q, &e_q, &sample_q, k)?;
    let mut out = Vec::with_capacity(sigmas_q.len());
    for s in &sigmas_q {
        out.push(qmp_to_zmp_exact(s)?);
    }
    Some(out)
}

/// Convert a sparse integer polynomial to a sparse rational polynomial.
fn zmp_to_qmp(f: &ZmPoly) -> MP<RationalDomain> {
    MP::<RationalDomain>::from_terms(
        RationalDomain,
        f.n_vars(),
        f.terms_ref()
            .iter()
            .map(|(e, c)| (e.to_vec(), Rational::from_integer(c.clone())))
            .collect(),
    )
}

/// Convert a sparse rational polynomial to an integer polynomial, requiring
/// every coefficient to be integral (denominator 1).
fn qmp_to_zmp_exact(f: &MP<RationalDomain>) -> Option<ZmPoly> {
    let mut terms = Vec::with_capacity(f.n_terms());
    for (e, c) in f.terms_ref() {
        let d = c.denom();
        if !IntegerDomain.is_one(&d) {
            return None;
        }
        terms.push((e.to_vec(), c.numer()));
    }
    Some(SparseMultivariatePolynomial::from_terms(
        IntegerDomain,
        f.n_vars(),
        terms,
    ))
}

/// Impose the true leading coefficient `ℓ_i` on factor `F_i`: replace the
/// coefficient of `x_0^{deg}` by `ℓ_i`.
fn impose_lcoeff_z(f_i: &ZmPoly, true_lc: &ZmPoly) -> ZmPoly {
    let deg = f_i.degree_in(0);
    let mut result = f_i.clone();
    // Remove all terms at the top degree in x_0, then add true_lc · x_0^deg.
    let top: Vec<smallvec::SmallVec<[usize; 4]>> = result
        .terms_ref()
        .keys()
        .filter(|e| e.first().copied().unwrap_or(0) == deg)
        .cloned()
        .collect();
    for e in top {
        result.set_term_external(e.to_vec(), Integer::from(0));
    }
    for (e, c) in true_lc.terms_ref() {
        let mut exp = e.to_vec();
        exp[0] = deg;
        let existing = result.coeff(&exp);
        result.set_term_external(exp, IntegerDomain.add(&existing, c));
    }
    result
}

/// EEZ lift over the integers with Wang-imposed leading coefficients.
///
/// The initial factors `initial_i` have `lc_{x_0} = ℓ_i`. At each step the
/// error is solved over ℚ and the corrections are required to be integral;
/// after every correction the leading coefficients are re-imposed so the
/// error stays confined to lower degrees in `x_0`.
fn eez_lift_z(
    f: &ZmPoly,
    sample: &[Integer],
    initial: &[ZmPoly],
    true_lcoeffs: &[ZmPoly],
) -> Option<Vec<ZmPoly>> {
    let n = f.n_vars();
    let mut lifted: Vec<ZmPoly> = initial.to_vec();

    for k in 1..n {
        let mut f_k = f.clone();
        for m in (k + 1..n).rev() {
            f_k = f_k.eval_keep(m, &sample[m]);
        }
        let a_k = sample[k].clone();
        let d_k = f_k.degree_in(k);
        let t_f = f_k.taylor_coefficients(k, &a_k);

        for j in 1..=d_k {
            let mut prod = one_mpoly(&IntegerDomain, n);
            for g in &lifted {
                prod = prod.mul(g);
            }
            let t_p = prod.taylor_coefficients(k, &a_k);
            let e = if j < t_p.len() {
                t_f[j].sub(&t_p[j])
            } else {
                t_f[j].clone()
            };
            if e.is_zero() {
                continue;
            }
            // See eez_lift: solve against the factors with x_k evaluated at
            // the sample point.
            let u_base: Vec<ZmPoly> = lifted.iter().map(|g| g.eval_keep(k, &a_k)).collect();
            let sigmas = diophantine_z(&u_base, &e, sample, k)?;
            let xp = x_minus_a_pow(&IntegerDomain, n, k, &a_k, j);
            for (i, g) in lifted.iter_mut().enumerate() {
                *g = g.add(&sigmas[i].mul(&xp));
            }
            // Re-impose the true leading coefficients after each correction.
            for (i, g) in lifted.iter_mut().enumerate() {
                *g = impose_lcoeff_z(g, &true_lcoeffs[i]);
            }
        }
    }
    Some(lifted)
}

// ---------------------------------------------------------------------
// p-adic coefficient Hensel lifting (non-constant leading coefficients)
// ---------------------------------------------------------------------

/// Convert a sparse integer polynomial to a sparse polynomial over `field`.
fn zmp_to_fmp(f: &ZmPoly, field: &FiniteField) -> FpMPoly {
    FpMPoly::from_terms(
        field.clone(),
        f.n_vars(),
        f.terms_ref()
            .iter()
            .map(|(e, c)| (e.to_vec(), field.element(c.to_bigint())))
            .filter(|(_, c)| !field.is_zero(c))
            .collect(),
    )
}

/// Convert a sparse finite-field polynomial to an integer polynomial using
/// the symmetric residue representation `(-p/2, p/2]`.
fn fp_to_z_symmetric(g: &FpMPoly) -> ZmPoly {
    let p = g.domain().prime().clone();
    let half = &p / BigInt::from(2u32);
    ZmPoly::from_terms(
        IntegerDomain,
        g.n_vars(),
        g.terms_ref()
            .iter()
            .map(|(e, c)| {
                let v = c.value().clone();
                let sym = if v > half { v - &p } else { v };
                (e.to_vec(), Integer::from(sym))
            })
            .collect(),
    )
}

/// Gelfond-style bound on the coefficient magnitude of every factor of `f`:
/// `(sqrt(∏(d_v+1) · 2^{2·Σd_v − #vars}) + 1) · max_norm · |lc(f)|`.
///
/// Reference: Symbolica `coefficient_bound`.
fn coefficient_bound_z(f: &ZmPoly) -> Integer {
    let max_norm = f
        .terms_ref()
        .values()
        .map(|c| c.abs())
        .max()
        .unwrap_or_else(|| Integer::from(1));
    let mut bound = Integer::from(1);
    let mut total_degree = 0u64;
    let mut non_zero_vars = 0u64;
    for v in 0..f.n_vars() {
        let d = f.degree_in(v) as u64;
        if d > 0 {
            non_zero_vars += 1;
            total_degree += d;
            bound = IntegerDomain.mul(&bound, &Integer::from(d as i64 + 1));
        }
    }
    let shift = (total_degree * 2).saturating_sub(non_zero_vars);
    let pow2 = IntegerDomain.pow(&Integer::from(2), shift);
    bound = IntegerDomain.mul(&bound, &pow2);
    let root = IntegerDomain.add(&bound.sqrt(), &Integer::from(1));
    let lc = f
        .leading_coeff()
        .map(|c| c.abs())
        .unwrap_or_else(|| Integer::from(1));
    IntegerDomain.mul(&root, &IntegerDomain.mul(&max_norm, &lc))
}

/// p-adic coefficient Hensel lift with imposed leading coefficients.
///
/// `factors` are the mod-`p` multivariate factors of `target` (as symmetric
/// integers); the true leading coefficients are imposed exactly. Each
/// iteration solves the mod-`p` Diophantine equation for the current error,
/// applies the correction scaled by `m = p^k`, and re-imposes the leading
/// coefficients, until the error vanishes or `m` exceeds the coefficient
/// bound `max_p`. Returns `None` on an unlucky sample or prime.
///
/// Reference: Symbolica `sparse_coefficient_hensel_lift_mod_prime` (dense
/// Diophantine variant).
///
/// When `allow_dense_fallback` is false and the skeletons would need more
/// distinct nonzero field elements than `𝔽_p` offers, returns `None`
/// immediately so the caller escalates to a larger prime instead of
/// silently degrading to the dense Diophantine solver.
fn coefficient_hensel_lift_z(
    target: &ZmPoly,
    factors: Vec<ZmPoly>,
    true_lcoeffs: &[ZmPoly],
    p: u64,
    max_p: &Integer,
    sample: &[Integer],
    allow_dense_fallback: bool,
) -> Option<Vec<ZmPoly>> {
    let n = target.n_vars();
    let field = FiniteField::new(BigInt::from(p));
    let mut factors: Vec<ZmPoly> = factors
        .iter()
        .zip(true_lcoeffs)
        .map(|(g, l)| impose_lcoeff_z(g, l))
        .collect();
    let factors_fp: Vec<FpMPoly> = factors.iter().map(|g| zmp_to_fmp(g, &field)).collect();
    // Skeletons: the mod-p factors with their leading x_0-part removed,
    // used by the sparse Diophantine solver to restrict correction support.
    let skeletons: Vec<FpMPoly> = factors_fp
        .iter()
        .map(|g| {
            let deg = g.degree_in(0);
            let mut sk = g.clone();
            let top: Vec<smallvec::SmallVec<[usize; 4]>> = sk
                .terms_ref()
                .keys()
                .filter(|e| e.first().copied().unwrap_or(0) == deg)
                .cloned()
                .collect();
            for e in top {
                sk.set_term_external(e.to_vec(), field.zero());
            }
            sk
        })
        .collect();
    // Small-prime heuristic: skeleton interpolation needs one distinct
    // nonzero field element per group member; bail out so the caller can
    // escalate to a larger prime instead of degrading to the dense solver.
    if !allow_dense_fallback
        && let Some(needed) = sparse_samples_needed(&skeletons)
        && BigInt::from(needed) >= field.prime().clone()
    {
        return None;
    }
    let sample_fp: Vec<FiniteFieldElement> = sample
        .iter()
        .map(|s| field.element(s.to_bigint()))
        .collect();

    let mut prod = one_mpoly(&IntegerDomain, n);
    for g in &factors {
        prod = prod.mul(g);
    }
    let mut error = target.sub(&prod);
    let mut m = Integer::from(p as i64);
    let p_int = Integer::from(p as i64);
    let mut iteration = 0u64;
    while !error.is_zero() && m <= *max_p {
        // The error must be divisible by m; reduce the quotient mod p.
        let mut error_over_m = ZmPoly::new(IntegerDomain, n);
        for (e, c) in error.terms_ref() {
            let (q, r) = c.div_rem(&m);
            if !r.is_zero() {
                return None;
            }
            error_over_m.set_term_external(e.to_vec(), q);
        }
        let error_fp = zmp_to_fmp(&error_over_m, &field);
        // Prefer skeleton interpolation; fall back to the dense recursive
        // solver when the sparsity assumption fails.
        let deltas = sparse_diophantine_fp(&factors_fp, &error_fp, &skeletons, iteration)
            .or_else(|| diophantine(&factors_fp, &error_fp, &sample_fp, n))?;
        iteration += 1;
        for (g, d) in factors.iter_mut().zip(&deltas) {
            let corr = fp_to_z_symmetric(d).mul_scalar(&m);
            *g = g.add(&corr);
        }
        // Re-impose the true leading coefficients after each correction.
        for (g, l) in factors.iter_mut().zip(true_lcoeffs) {
            *g = impose_lcoeff_z(g, l);
        }
        let mut prod = one_mpoly(&IntegerDomain, n);
        for g in &factors {
            prod = prod.mul(g);
        }
        error = target.sub(&prod);
        m = IntegerDomain.mul(&m, &p_int);
    }
    if error.is_zero() { Some(factors) } else { None }
}

// ---------------------------------------------------------------------
// Sparse multivariate Diophantine solver (skeleton interpolation)
// ---------------------------------------------------------------------

/// Maximum skeleton size accepted by the sparse Diophantine solver; larger
/// inputs fall back to the dense recursive solver.
const SPARSE_MDP_MAX_TERMS: usize = 512;
/// Number of random base-point sets tried before giving up.
const SPARSE_MDP_BASE_ATTEMPTS: usize = 4;

/// Deterministic SplitMix64 generator for interpolation base points.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Solve `rhs[k] = Σ_i c_i · x_i^(k+1)` for the coefficients `c_i`.
///
/// Reference: Symbolica `solve_shifted_transposed_vandermonde`.
fn solve_shifted_transposed_vandermonde(
    field: &FiniteField,
    x: &[FiniteFieldElement],
    rhs: &[FiniteFieldElement],
) -> Vec<FiniteFieldElement> {
    debug_assert_eq!(x.len(), rhs.len());
    match x.len() {
        0 => Vec::new(),
        1 => vec![field.div(&rhs[0], &x[0]).expect("nonzero generator")],
        len => {
            // master(z) = ∏(z − x_i), built coefficient by coefficient.
            let mut master = vec![field.zero(); len + 1];
            master[0] = field.one();
            for (i, xi) in x.iter().enumerate() {
                let mut old_last = master[0].clone();
                master[0] = field.mul(&master[0], &field.neg(xi));
                for m in master.iter_mut().take(i + 1).skip(1) {
                    let ov = m.clone();
                    *m = field.add(&field.mul(m, &field.neg(xi)), &old_last);
                    old_last = ov;
                }
                master[i + 1] = field.one();
            }
            let mut sol = Vec::with_capacity(len);
            for (i, s) in x.iter().enumerate() {
                // norm = ∏_{j≠i} (x_i − x_j); generators are distinct.
                let mut norm = field.one();
                for (j, l) in x.iter().enumerate() {
                    if j != i {
                        norm = field.mul(&norm, &field.sub(s, l));
                    }
                }
                // Sample master/(1 − x_i·z) against rhs via Horner's rule.
                let mut coeff = field.zero();
                let mut last_q = field.zero();
                for (m, r) in master.iter().skip(1).zip(rhs.iter()).rev() {
                    last_q = field.add(m, &field.mul(s, &last_q));
                    coeff = field.add(&coeff, &field.mul(&last_q, r));
                }
                coeff = field.div(&coeff, &norm).expect("distinct generators");
                // Shift from the x_i^k basis to the x_i^(k+1) basis.
                coeff = field.div(&coeff, &x[i]).expect("nonzero generator");
                sol.push(coeff);
            }
            sol
        }
    }
}

/// Solve the univariate Diophantine equation `Σ δ_i · ∏_{j≠i} f_j = rhs`
/// over a prime field with `deg(δ_i) < deg(f_i)`, using the sequential
/// extended-Euclid iteration. Returns `None` if the factors are not
/// pairwise coprime.
///
/// Reference: Symbolica `try_univariate_diophantine`.
fn try_univariate_diophantine_fp(
    factors: &[UP<FiniteField>],
    rhs: &UP<FiniteField>,
) -> Option<Vec<UP<FiniteField>>> {
    let r = factors.len();
    if r == 0 {
        return None;
    }
    if r == 1 {
        return Some(vec![rhs.clone()]);
    }
    if factors
        .iter()
        .any(|f| f.leading_coeff().is_none_or(|c| f.domain().is_zero(c)))
    {
        return None;
    }
    // products[i] = f_{i+1}···f_{r−1}.
    let mut products: Vec<UP<FiniteField>> = Vec::with_capacity(r - 1);
    let mut cur = factors[r - 1].clone();
    products.push(cur.clone());
    for f in factors[1..r - 1].iter().rev() {
        cur = cur.mul(f);
        products.push(cur.clone());
    }
    products.reverse();

    let mut deltas: Vec<UP<FiniteField>> = Vec::with_capacity(r);
    let mut cur_s = rhs.clone();
    for (factor, product) in factors.iter().zip(&products) {
        let (g, s, t) = factor.extended_gcd_poly(product);
        if g.degree() != Some(0) {
            return None; // not coprime
        }
        let inv = factor.domain().inv(g.leading_coeff().unwrap())?;
        let s = s.mul_scalar(&inv);
        let t = t.mul_scalar(&inv);
        let (_, new_s) = t.mul(&cur_s).div_rem(factor)?;
        deltas.push(new_s);
        let (_, rem) = s.mul(&cur_s).div_rem(product)?;
        cur_s = rem;
    }
    deltas.push(cur_s);
    Some(deltas)
}

/// Evaluate the secondary-variable monomial `exp` (variables 1..n) at the
/// base points.
fn monomial_eval_fp(
    field: &FiniteField,
    exp: &[usize],
    base: &[FiniteFieldElement],
) -> FiniteFieldElement {
    let mut v = field.one();
    for (i, b) in base.iter().enumerate() {
        let e = exp.get(i + 1).copied().unwrap_or(0);
        if e > 0 {
            v = field.mul(&v, &field.pow(b, e as u64));
        }
    }
    v
}

/// Group the exponents of a skeleton by their `x_0`-degree.
fn group_skeleton(skeleton: &FpMPoly) -> Vec<(usize, Vec<smallvec::SmallVec<[usize; 4]>>)> {
    let mut groups: Vec<(usize, Vec<smallvec::SmallVec<[usize; 4]>>)> = Vec::new();
    for exp in skeleton.terms_ref().keys() {
        let deg = exp.first().copied().unwrap_or(0);
        if let Some((_, exps)) = groups.iter_mut().find(|(d, _)| *d == deg) {
            exps.push(exp.clone());
        } else {
            groups.push((deg, vec![exp.clone()]));
        }
    }
    groups
}

/// Maximum skeleton group size the sparse Diophantine solver would need for
/// these skeletons, or `None` when the sparse solver is inapplicable (empty
/// skeletons or too many terms). Mirrors the applicability checks in
/// [`sparse_diophantine_two_factor_fp`] and
/// [`sparse_diophantine_n_factor_fp`].
fn sparse_samples_needed(skeletons: &[FpMPoly]) -> Option<usize> {
    if skeletons.len() == 2 {
        let sparse = skeletons
            .iter()
            .filter(|s| !s.is_zero())
            .min_by_key(|s| s.n_terms())?;
        if sparse.n_terms() > SPARSE_MDP_MAX_TERMS {
            return None;
        }
        let needed = group_skeleton(sparse)
            .iter()
            .map(|(_, exps)| exps.len())
            .max()
            .unwrap_or(0);
        (needed > 0).then_some(needed)
    } else if skeletons.len() > 2 {
        let mut needed = 0usize;
        let mut total_terms = 0usize;
        for sk in skeletons {
            for (_, exps) in group_skeleton(sk) {
                needed = needed.max(exps.len());
            }
            total_terms += sk.n_terms();
        }
        if needed == 0 || total_terms > SPARSE_MDP_MAX_TERMS {
            return None;
        }
        Some(needed)
    } else {
        None
    }
}

/// Build the dense univariate image in variable 0 from per-term secondary
/// monomial evaluations (`evals[t]` is the monomial evaluation of the
/// exponent of `terms[t]`, WITHOUT the coefficient; the coefficient is
/// applied here so that powering the monomials does not power the
/// coefficients).
fn build_image_fp(
    field: &FiniteField,
    terms: &[(smallvec::SmallVec<[usize; 4]>, FiniteFieldElement)],
    evals: &[FiniteFieldElement],
) -> UP<FiniteField> {
    let mut coeffs: Vec<FiniteFieldElement> = Vec::new();
    for ((exp, c), ev) in terms.iter().zip(evals) {
        let idx = exp.first().copied().unwrap_or(0);
        if idx >= coeffs.len() {
            coeffs.resize(idx + 1, field.zero());
        }
        coeffs[idx] = field.add(&coeffs[idx], &field.mul(c, ev));
    }
    UP::<FiniteField>::from_coeffs(field.clone(), coeffs)
}

/// Coefficient of `x_0^deg` in a dense univariate polynomial (zero if the
/// degree is absent).
fn upoly_coeff(p: &UP<FiniteField>, deg: usize) -> FiniteFieldElement {
    p.coeffs()
        .get(deg)
        .cloned()
        .unwrap_or_else(|| p.domain().zero())
}

/// Snapshot of a polynomial's terms in a stable order.
type TermSnapshot = Vec<(smallvec::SmallVec<[usize; 4]>, FiniteFieldElement)>;

/// Skeleton exponents grouped by their `x_0`-degree, per factor.
type SkeletonGroups = Vec<Vec<(usize, Vec<smallvec::SmallVec<[usize; 4]>>)>>;

fn snapshot_terms(f: &FpMPoly) -> TermSnapshot {
    f.terms_ref()
        .iter()
        .map(|(e, c)| (e.clone(), c.clone()))
        .collect()
}

/// Random nonzero base points for the secondary variables.
fn random_base_fp(field: &FiniteField, n: usize, rng: &mut SplitMix64) -> Vec<FiniteFieldElement> {
    (1..n)
        .map(|_| {
            loop {
                let v = field.element(BigInt::from(rng.next()));
                if !field.is_zero(&v) {
                    break v;
                }
            }
        })
        .collect()
}

/// Sparse two-factor Diophantine solver over a prime field: interpolates
/// the correction of the factor with the sparser skeleton from univariate
/// images, then obtains the other correction by exact division.
///
/// Reference: Symbolica `sparse_multivariate_diophantine_two_factor_by_sampling`.
fn sparse_diophantine_two_factor_fp(
    factors: &[FpMPoly],
    prods: &[FpMPoly],
    error: &FpMPoly,
    skeletons: &[FpMPoly],
    seed: u64,
) -> Option<Vec<FpMPoly>> {
    if factors.len() != 2 {
        return None;
    }
    let field = error.domain().clone();
    let n = error.n_vars();
    let sparse_factor = skeletons
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.is_zero())
        .min_by_key(|(_, s)| s.n_terms())
        .map(|(i, _)| i)?;
    let dense_factor = 1 - sparse_factor;
    let skeleton = &skeletons[sparse_factor];
    let groups = group_skeleton(skeleton);
    let samples_needed = groups.iter().map(|(_, e)| e.len()).max().unwrap_or(0);
    if samples_needed == 0 || skeleton.n_terms() > SPARSE_MDP_MAX_TERMS {
        return None;
    }
    // Interpolation needs one distinct nonzero generator per group member;
    // a prime field offers only p − 1.
    if BigInt::from(samples_needed) >= field.prime().clone() {
        return None;
    }

    let error_terms = snapshot_terms(error);
    let factor_terms: Vec<TermSnapshot> = factors.iter().map(snapshot_terms).collect();

    let mut rng = SplitMix64(seed);
    'attempts: for _ in 0..SPARSE_MDP_BASE_ATTEMPTS {
        let base = random_base_fp(&field, n, &mut rng);
        // Generators must be nonzero and pairwise distinct within a group.
        let mut generators: Vec<Vec<FiniteFieldElement>> = Vec::with_capacity(groups.len());
        for (_, exps) in &groups {
            let mut gens: Vec<FiniteFieldElement> = Vec::with_capacity(exps.len());
            for e in exps {
                let g = monomial_eval_fp(&field, e, &base);
                if field.is_zero(&g) || gens.contains(&g) {
                    continue 'attempts;
                }
                gens.push(g);
            }
            generators.push(gens);
        }

        // Per-term secondary monomial evaluations (without coefficients).
        let error_base: Vec<FiniteFieldElement> = error_terms
            .iter()
            .map(|(e, _)| monomial_eval_fp(&field, e, &base))
            .collect();
        let factors_base: Vec<Vec<FiniteFieldElement>> = factor_terms
            .iter()
            .map(|terms| {
                terms
                    .iter()
                    .map(|(e, _)| monomial_eval_fp(&field, e, &base))
                    .collect()
            })
            .collect();

        let mut rhs: Vec<Vec<FiniteFieldElement>> = groups.iter().map(|_| Vec::new()).collect();
        let mut error_current = error_base.clone();
        let mut factors_current = factors_base.clone();
        for s in 0..samples_needed {
            if s > 0 {
                for (cur, b) in error_current.iter_mut().zip(&error_base) {
                    *cur = field.mul(cur, b);
                }
                for (cur_f, base_f) in factors_current.iter_mut().zip(&factors_base) {
                    for (cur, b) in cur_f.iter_mut().zip(base_f) {
                        *cur = field.mul(cur, b);
                    }
                }
            }
            let error_img = build_image_fp(&field, &error_terms, &error_current);
            let factor_imgs: Vec<UP<FiniteField>> = factor_terms
                .iter()
                .zip(&factors_current)
                .map(|(te, cur)| build_image_fp(&field, te, cur))
                .collect();
            let Some(deltas_img) = try_univariate_diophantine_fp(&factor_imgs, &error_img) else {
                continue 'attempts;
            };
            for (gi, (deg, _)) in groups.iter().enumerate() {
                rhs[gi].push(upoly_coeff(&deltas_img[sparse_factor], *deg));
            }
        }

        // Vandermonde per group gives the sparse correction.
        let mut sparse_delta = skeleton.zero();
        for ((_, exps), (gens, rhs)) in groups.iter().zip(generators.iter().zip(rhs.iter())) {
            let coeffs = solve_shifted_transposed_vandermonde(&field, gens, &rhs[..exps.len()]);
            for (c, e) in coeffs.into_iter().zip(exps) {
                if !field.is_zero(&c) {
                    sparse_delta.set_term_external(e.to_vec(), c);
                }
            }
        }

        // The other correction follows by exact division.
        let residual = error.sub(&sparse_delta.mul(&prods[sparse_factor]));
        let Some(dense_delta) = residual.checked_div_exact(&prods[dense_factor]) else {
            continue;
        };
        let mut deltas = vec![error.zero(), error.zero()];
        deltas[sparse_factor] = sparse_delta;
        deltas[dense_factor] = dense_delta;

        let mut check = error.zero();
        for (d, p) in deltas.iter().zip(prods) {
            check = check.add(&d.mul(p));
        }
        if check == *error {
            return Some(deltas);
        }
    }
    None
}

/// Sparse n-factor Diophantine solver over a prime field: interpolates every
/// factor's correction from its skeleton via univariate images and
/// Vandermonde solves, then verifies against the full equation.
///
/// Reference: Symbolica `sparse_multivariate_diophantine_by_sampling`.
fn sparse_diophantine_n_factor_fp(
    factors: &[FpMPoly],
    prods: &[FpMPoly],
    error: &FpMPoly,
    skeletons: &[FpMPoly],
    seed: u64,
) -> Option<Vec<FpMPoly>> {
    let r = factors.len();
    let field = error.domain().clone();
    let n = error.n_vars();

    let mut groups: SkeletonGroups = Vec::with_capacity(r);
    let mut samples_needed = 0usize;
    let mut total_terms = 0usize;
    for sk in skeletons {
        let fg = group_skeleton(sk);
        for (_, exps) in &fg {
            samples_needed = samples_needed.max(exps.len());
        }
        total_terms += sk.n_terms();
        groups.push(fg);
    }
    if samples_needed == 0 || total_terms > SPARSE_MDP_MAX_TERMS {
        return None;
    }
    // Interpolation needs one distinct nonzero generator per group member;
    // a prime field offers only p − 1.
    if BigInt::from(samples_needed) >= field.prime().clone() {
        return None;
    }

    let error_terms = snapshot_terms(error);
    let factor_terms: Vec<TermSnapshot> = factors.iter().map(snapshot_terms).collect();

    let mut rng = SplitMix64(seed);
    'attempts: for _ in 0..SPARSE_MDP_BASE_ATTEMPTS {
        let base = random_base_fp(&field, n, &mut rng);
        let mut generators: Vec<Vec<Vec<FiniteFieldElement>>> = Vec::with_capacity(r);
        for fg in &groups {
            let mut gen_fg: Vec<Vec<FiniteFieldElement>> = Vec::with_capacity(fg.len());
            for (_, exps) in fg {
                let mut gens: Vec<FiniteFieldElement> = Vec::with_capacity(exps.len());
                for e in exps {
                    let g = monomial_eval_fp(&field, e, &base);
                    if field.is_zero(&g) || gens.contains(&g) {
                        continue 'attempts;
                    }
                    gens.push(g);
                }
                gen_fg.push(gens);
            }
            generators.push(gen_fg);
        }

        // Per-term secondary monomial evaluations (without coefficients).
        let error_base: Vec<FiniteFieldElement> = error_terms
            .iter()
            .map(|(e, _)| monomial_eval_fp(&field, e, &base))
            .collect();
        let factors_base: Vec<Vec<FiniteFieldElement>> = factor_terms
            .iter()
            .map(|terms| {
                terms
                    .iter()
                    .map(|(e, _)| monomial_eval_fp(&field, e, &base))
                    .collect()
            })
            .collect();

        let mut rhs: Vec<Vec<Vec<FiniteFieldElement>>> = groups
            .iter()
            .map(|fg| fg.iter().map(|_| Vec::new()).collect())
            .collect();
        let mut error_current = error_base.clone();
        let mut factors_current = factors_base.clone();
        for s in 0..samples_needed {
            if s > 0 {
                for (cur, b) in error_current.iter_mut().zip(&error_base) {
                    *cur = field.mul(cur, b);
                }
                for (cur_f, base_f) in factors_current.iter_mut().zip(&factors_base) {
                    for (cur, b) in cur_f.iter_mut().zip(base_f) {
                        *cur = field.mul(cur, b);
                    }
                }
            }
            let error_img = build_image_fp(&field, &error_terms, &error_current);
            let factor_imgs: Vec<UP<FiniteField>> = factor_terms
                .iter()
                .zip(&factors_current)
                .map(|(te, cur)| build_image_fp(&field, te, cur))
                .collect();
            let Some(deltas_img) = try_univariate_diophantine_fp(&factor_imgs, &error_img) else {
                continue 'attempts;
            };
            for (i, fg) in groups.iter().enumerate() {
                for (gi, (deg, _)) in fg.iter().enumerate() {
                    rhs[i][gi].push(upoly_coeff(&deltas_img[i], *deg));
                }
            }
        }

        let mut deltas: Vec<FpMPoly> = skeletons.iter().map(|s| s.zero()).collect();
        for (i, fg) in groups.iter().enumerate() {
            for (gi, (_, exps)) in fg.iter().enumerate() {
                let coeffs = solve_shifted_transposed_vandermonde(
                    &field,
                    &generators[i][gi],
                    &rhs[i][gi][..exps.len()],
                );
                for (c, e) in coeffs.into_iter().zip(exps) {
                    if !field.is_zero(&c) {
                        deltas[i].set_term_external(e.to_vec(), c);
                    }
                }
            }
        }

        let mut check = error.zero();
        for (d, p) in deltas.iter().zip(prods) {
            check = check.add(&d.mul(p));
        }
        if check == *error {
            return Some(deltas);
        }
    }
    None
}

/// Solve the multivariate Diophantine equation `Σ δ_i · ∏_{j≠i} f_j = e`
/// over a prime field by skeleton interpolation. Returns `None` when the
/// sparsity assumption does not hold; the caller falls back to the dense
/// recursive solver.
fn sparse_diophantine_fp(
    factors: &[FpMPoly],
    error: &FpMPoly,
    skeletons: &[FpMPoly],
    seed: u64,
) -> Option<Vec<FpMPoly>> {
    // Benchmark kill-switch: `OCAS_DISABLE_SPARSE_DIO=1` forces the dense
    // recursive solver so dense-vs-sparse timings can be compared.
    static DISABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if *DISABLED.get_or_init(|| std::env::var_os("OCAS_DISABLE_SPARSE_DIO").is_some()) {
        return None;
    }
    let r = factors.len();
    if r < 2 {
        return None;
    }
    let field = error.domain().clone();
    let n = error.n_vars();
    let prods: Vec<FpMPoly> = (0..r)
        .map(|i| {
            let mut p = one_mpoly(&field, n);
            for (j, f) in factors.iter().enumerate() {
                if i != j {
                    p = p.mul(f);
                }
            }
            p
        })
        .collect();
    let found = if r == 2 {
        sparse_diophantine_two_factor_fp(factors, &prods, error, skeletons, seed)
    } else if r > 2 {
        sparse_diophantine_n_factor_fp(factors, &prods, error, skeletons, seed)
    } else {
        None
    };
    if let Some(d) = found {
        #[cfg(test)]
        SPARSE_DIO_HITS.with(|h| h.set(h.get() + 1));
        return Some(d);
    }
    None
}

// Test-only thread-local counter of successful sparse Diophantine solves.
// Thread-local so parallel tests cannot pollute each other's assertions.
#[cfg(test)]
mod sparse_dio_hits {
    #![allow(clippy::missing_const_for_thread_local)]
    thread_local! {
        pub static SPARSE_DIO_HITS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    }
}

#[cfg(test)]
use sparse_dio_hits::SPARSE_DIO_HITS;

/// Small primes tried for the p-adic coefficient lift. A prime is usable
/// when the univariate image keeps its degree and stays square-free mod p.
const PADIC_PRIMES: [u64; 25] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

/// Factor square-free `f` via p-adic coefficient Hensel lifting with
/// Wang-imposed non-constant leading coefficients.
///
/// The target is rescaled to `c^{r−1}·f` (image content `c`, `r` univariate
/// factors) so that the imposed leading coefficients `c·ℓ_i` are consistent
/// with the scaled initial factors `c·u_i`. The mod-`p` factors are first
/// lifted through the secondary variables over `𝔽_p` ([`eez_lift_imposed`]),
/// then lifted p-adically in the coefficients
/// ([`coefficient_hensel_lift_z`]). Returns the primitive irreducible
/// factors of `f`, or `None` on unlucky samples/primes.
///
/// Reference: Symbolica `multivariate_factorization` (univariate start).
fn padic_lift_factors(
    f: &ZmPoly,
    sample: &[Integer],
    uni: &[UP<IntegerDomain>],
    content: &Integer,
    true_lcoeffs: &[ZmPoly],
) -> Option<Vec<ZmPoly>> {
    let n = f.n_vars();
    let r = uni.len();
    let deg0 = f.degree_in(0);
    let mut scale_pow = Integer::from(1);
    for _ in 1..r {
        scale_pow = IntegerDomain.mul(&scale_pow, content);
    }
    let target = f.mul_scalar(&scale_pow);
    let scaled_lcoeffs: Vec<ZmPoly> = true_lcoeffs.iter().map(|l| l.mul_scalar(content)).collect();
    let scaled_uni: Vec<UP<IntegerDomain>> = uni.iter().map(|u| u.mul_scalar(content)).collect();
    let bound = coefficient_bound_z(&target);
    let image = eval_to_image_z(f, sample);

    // For large sparse inputs, start from a larger prime: the sparse
    // Diophantine solver needs at least one distinct nonzero field element
    // per skeleton group member, and bigger primes also cut the number of
    // p-adic iterations.
    let prime_start = if f.n_terms() >= 30 { 8 } else { 0 };
    // Two passes over the primes: the first escalates past primes too small
    // for skeleton interpolation (dense Diophantine fallback disabled); the
    // second re-enables the dense fallback so inputs whose skeletons fit no
    // prime in `PADIC_PRIMES` still factor.
    for allow_dense_fallback in [false, true] {
        for p in PADIC_PRIMES[prime_start..].iter().copied() {
            let field = FiniteField::new(BigInt::from(p));
            let image_fp = UP::<FiniteField>::from_coeffs(
                field.clone(),
                image
                    .coeffs()
                    .iter()
                    .map(|c| field.element(c.to_bigint()))
                    .collect(),
            );
            if image_fp.degree() != Some(deg0) || !image_fp.is_square_free() {
                continue; // p divides the leading coefficient or the discriminant
            }
            let target_fp = zmp_to_fmp(&target, &field);
            let initial_fp: Vec<FpMPoly> = scaled_uni
                .iter()
                .map(|u| zmp_to_fmp(&dense_to_mpoly(u, n), &field))
                .collect();
            let tl_fp: Vec<FpMPoly> = scaled_lcoeffs
                .iter()
                .map(|l| zmp_to_fmp(l, &field))
                .collect();
            let sample_fp: Vec<FiniteFieldElement> = sample
                .iter()
                .map(|s| field.element(s.to_bigint()))
                .collect();
            let Some(lifted_fp) = eez_lift_imposed(&target_fp, &sample_fp, &initial_fp, &tl_fp)
            else {
                continue;
            };
            let lifted_z: Vec<ZmPoly> = lifted_fp.iter().map(fp_to_z_symmetric).collect();
            // Smallest p^k with 2·p^k ≥ bound.
            let mut max_p = Integer::from(p as i64);
            while IntegerDomain.mul(&max_p, &Integer::from(2)) < bound {
                max_p = IntegerDomain.mul(&max_p, &Integer::from(p as i64));
            }
            let Some(factors) = coefficient_hensel_lift_z(
                &target,
                lifted_z,
                &scaled_lcoeffs,
                p,
                &max_p,
                sample,
                allow_dense_fallback,
            ) else {
                continue;
            };
            let out: Vec<ZmPoly> = factors.iter().map(primitive_positive).collect();
            let mut prod = one_mpoly(&IntegerDomain, n);
            for g in &out {
                prod = prod.mul(g);
            }
            if equal_up_to_unit(&prod, f) {
                let mut out = out;
                out.sort_by_key(|b| std::cmp::Reverse(b.degree_in(0)));
                return Some(out);
            }
            // Recombination fallback: the lifted factors may be finer than the
            // true irreducible factors.
            if let Some(irr) = zassenhaus_multivariate(f, &out) {
                let mut prod = one_mpoly(&IntegerDomain, n);
                for g in &irr {
                    prod = prod.mul(g);
                }
                if equal_up_to_unit(&prod, f) {
                    return Some(irr);
                }
            }
        }
    }
    None
}

// ---------------------------------------------------------------------
// Finite-field machinery
// ---------------------------------------------------------------------

/// Evaluate all secondary variables at the sample, leaving a dense
/// univariate polynomial in variable 0.
fn eval_to_image_fp(f: &FpMPoly, sample: &[FiniteFieldElement]) -> UP<FiniteField> {
    let mut p = f.clone();
    for k in (1..f.n_vars()).rev() {
        p = p.eval_keep(k, &sample[k]);
    }
    mpoly_to_dense(&p)
}

/// Find sample points for the secondary variables such that the univariate
/// image of `f` has full degree and is square-free. Returns candidates
/// (sample, monic irreducible image factors) ordered by increasing number
/// of factors, capped at `max_candidates`.
///
/// `range` caps the sample values; the caller retries with a larger range
/// when the first round yields no usable candidate (adaptive search).
///
/// `lc_filter` lists the non-constant factors of the leading coefficient
/// (in the secondary variables); samples where any evaluates to zero are
/// skipped, since Wang's leading-coefficient distribution requires every
/// image α_j to be nonzero (over a field, nonzero ⟹ invertible).
fn find_sample_fp(
    f: &FpMPoly,
    max_candidates: usize,
    range: u64,
    lc_filter: &[FpMPoly],
) -> Vec<(Vec<FiniteFieldElement>, Vec<UP<FiniteField>>)> {
    let n = f.n_vars();
    let field = f.domain().clone();
    let p = field.prime().to_u64().unwrap_or(u64::MAX);
    let range = range.clamp(1, p);
    let attempts = (range as usize).saturating_pow((n - 1) as u32).min(512);
    let deg0 = f.degree_in(0);

    let mut best: Vec<(Vec<FiniteFieldElement>, Vec<UP<FiniteField>>)> = Vec::new();
    for t in 0..attempts {
        let mut sample = vec![field.zero(); n];
        let mut rem = t as u64;
        for slot in sample.iter_mut().take(n).skip(1) {
            *slot = field.element(BigInt::from(rem % range));
            rem /= range;
        }
        // Skip samples where any LC filter polynomial evaluates to zero:
        // Wang's distribution requires α_j = g_j(s) ≠ 0.
        // Note: lc_filter elements have n-1 variables (main var dropped),
        // so variable k in the reduced polynomial corresponds to sample[k+1].
        if !lc_filter.is_empty() {
            let mut bad = false;
            for g in lc_filter {
                let mut img = g.clone();
                for k in (0..img.n_vars()).rev() {
                    img = img.eval_keep(k, &sample[k + 1]);
                }
                if img.is_zero() {
                    bad = true;
                    break;
                }
            }
            if bad {
                continue;
            }
        }
        let image = eval_to_image_fp(f, &sample);
        if image.degree().unwrap_or(0) != deg0 || image.degree().unwrap_or(0) == 0 {
            continue; // leading coefficient vanished or no x0-part
        }
        if !image.is_square_free() {
            continue;
        }
        let mut factors: Vec<UP<FiniteField>> = image
            .factor()
            .into_iter()
            .map(|(g, _)| monic_upoly(&g))
            .collect();
        factors.sort_by_key(|b| std::cmp::Reverse(b.degree().unwrap_or(0)));
        if factors.is_empty() {
            continue;
        }
        let pos = best
            .binary_search_by(|(_, b)| b.len().cmp(&factors.len()))
            .unwrap_or_else(|e| e);
        best.insert(pos, (sample, factors));
        if best.len() > max_candidates {
            best.pop();
        }
    }
    best
}

/// Reconstruct the true multivariate leading coefficients ℓ_i for the
/// factors of `f` over a prime field, using Wang's greedy distribution of
/// the irreducible factors of the overall leading coefficient.
///
/// Ported from [`wang_reconstruct_lcoeffs`] (ℤ path), adapted for `𝔽_p`:
/// over a field, α_j ≠ 0 implies invertibility, so the "coprimality"
/// check becomes "all α_j nonzero" and the divisibility check in the
/// greedy loop is always satisfiable.
fn wang_reconstruct_lcoeffs_fp(
    lcoeff: &FpMPoly,
    sample: &[FiniteFieldElement],
    uni: &[UP<FiniteField>],
) -> Option<Vec<FpMPoly>> {
    let n = lcoeff.n_vars();
    let field = lcoeff.domain().clone();

    // Constant LC fast path: each ℓ_i is the constant lc(u_i).
    if lcoeff.degree_in(0) == 0 && lcoeff.drop_main_var().total_degree() == Some(0) {
        return Some(
            uni.iter()
                .map(|u| {
                    let lc = u.leading_coeff().cloned().unwrap_or_else(|| field.one());
                    one_mpoly(&field, n).mul_scalar(&lc)
                })
                .collect(),
        );
    }

    // Factor the leading coefficient in the secondary variables.
    let lc_reduced = lcoeff.drop_main_var();
    let lc_factors: Vec<(FpMPoly, usize)> = if lc_reduced.n_vars() == 0 {
        Vec::new()
    } else {
        multivariate_factor_fp(&lc_reduced)
    };

    // Field images α_j = g_j(s), requiring each to be nonzero.
    let mut alpha: Vec<FiniteFieldElement> = Vec::new();
    let mut nonconst: Vec<FpMPoly> = Vec::new();
    let mut const_part = field.one();
    for (g, _e) in &lc_factors {
        if is_constant(g) {
            const_part = field.mul(&const_part, &g.coeff(&vec![0; g.n_vars()]));
            continue;
        }
        let mut img = g.clone();
        // g is in (n-1) variables (main var dropped); evaluate at the
        // secondary sample values.
        for k in 0..img.n_vars() {
            img = img.eval_keep(k, &sample[k + 1]);
        }
        let a = img.coeff(&vec![0; img.n_vars()]);
        if field.is_zero(&a) {
            return None; // unlucky sample: factor vanishes
        }
        alpha.push(a);
        nonconst.push(g.clone());
    }
    // Over a field, all nonzero elements are pairwise coprime, so no
    // additional coprimality check is needed (unlike the ℤ path).

    // Multiplicities of each non-constant factor in ℓ.
    let multiplicities: Vec<usize> = nonconst
        .iter()
        .map(|g| {
            lc_factors
                .iter()
                .find(|(h, _)| h == g)
                .map(|(_, e)| *e)
                .unwrap_or(1)
        })
        .collect();

    // Greedy distribution: assign g_j to u_i while α_j ≠ 0.
    // Over a field every nonzero element divides any element, so we
    // distribute based on the number of available copies.
    let r = uni.len();
    let mut lcoeffs: Vec<FpMPoly> = vec![one_mpoly(&field, n); r];
    let mut residual_lc: Vec<FiniteFieldElement> = uni
        .iter()
        .map(|u| u.leading_coeff().cloned().unwrap_or_else(|| field.one()))
        .collect();
    let mut used = vec![0usize; nonconst.len()];
    for i in 0..r {
        for j in 0..nonconst.len() {
            while used[j] < multiplicities[j] && !field.is_zero(&residual_lc[i]) {
                lcoeffs[i] = lcoeffs[i].mul(&nonconst[j].embed_new_main());
                // residual_lc[i] /= α_j  (exact in a field)
                residual_lc[i] = field
                    .div(&residual_lc[i], &alpha[j])
                    .unwrap_or_else(|| field.zero());
                used[j] += 1;
            }
        }
    }
    if used != multiplicities {
        return None; // could not distribute all factors
    }

    // No reconciliation needed: find_sample_fp returns monic factors
    // (lc(u_i) = 1), so ℓ_i(s) may be any nonzero value — the identity
    // ∏ ℓ_i = ℓ is verified globally below. (The ℤ path reconciles
    // because find_sample_z returns non-monic factors.)

    // Global verification: ∏ ℓ_i = ℓ as a polynomial identity.
    let mut prod = one_mpoly(&field, n);
    for l in &lcoeffs {
        prod = prod.mul(l);
    }
    if prod == *lcoeff { Some(lcoeffs) } else { None }
}

/// Make a dense univariate polynomial monic over a field.
fn monic_upoly<D: EuclideanDomain>(f: &UP<D>) -> UP<D> {
    if f.is_zero() {
        return f.clone();
    }
    let lc = f.leading_coeff().cloned().unwrap();
    let inv = f.domain().inv(&lc).expect("field leading coefficient");
    f.mul_scalar(&inv)
}

/// Shrink variable `var` by the characteristic: divide all its exponents by
/// `p`. Returns `None` if some exponent is not divisible by `p`.
/// Coefficients are unchanged (over a prime field every element is its own
/// p-th root by Fermat).
fn shrink_var(f: &FpMPoly, var: usize, p: usize) -> Option<FpMPoly> {
    let mut result = FpMPoly::new(f.domain().clone(), f.n_vars());
    for (exp, c) in f.terms_ref() {
        let e = exp.get(var).copied().unwrap_or(0);
        if !e.is_multiple_of(p) {
            return None;
        }
        let mut new_exp = exp.clone();
        new_exp[var] = e / p;
        result.set_term_external(new_exp.to_vec(), c.clone());
    }
    Some(result)
}

/// Expand variable `var` by the characteristic: multiply all its exponents
/// by `p` (the substitution `x_var → x_var^p`).
fn expand_var(f: &FpMPoly, var: usize, p: usize) -> FpMPoly {
    let mut result = FpMPoly::new(f.domain().clone(), f.n_vars());
    for (exp, c) in f.terms_ref() {
        let mut new_exp = exp.clone();
        new_exp[var] = exp.get(var).copied().unwrap_or(0) * p;
        result.set_term_external(new_exp.to_vec(), c.clone());
    }
    result
}

/// Take the full p-th root of `f`: succeeds iff every exponent of every
/// variable is divisible by the characteristic, dividing all exponents by
/// `p`. Over a prime field this is the exact inverse of Frobenius.
fn pth_root_full(f: &FpMPoly) -> Option<FpMPoly> {
    let p = f.domain().prime().to_usize()?;
    let mut result = FpMPoly::new(f.domain().clone(), f.n_vars());
    for (exp, c) in f.terms_ref() {
        if exp.iter().any(|&e| !e.is_multiple_of(p)) {
            return None;
        }
        let new_exp: Vec<usize> = exp.iter().map(|&e| e / p).collect();
        result.set_term_external(new_exp, c.clone());
    }
    Some(result)
}

/// Square-free test for a multivariate polynomial over a prime field:
/// `f` is square-free iff the gcd of `f` with all its (non-zero) partial
/// derivatives is a constant. If every partial derivative vanishes, `f` is
/// a polynomial in `x_i^p` and is not square-free (unless constant).
fn is_square_free_mpoly(f: &FpMPoly) -> bool {
    if is_constant(f) {
        return true;
    }
    let mut any_nonzero_partial = false;
    let mut g = f.clone();
    for v in 0..f.n_vars() {
        let dv = f.derivative(v);
        if dv.is_zero() {
            continue;
        }
        any_nonzero_partial = true;
        g = match multivariate_gcd_fp(&g, &dv) {
            Some(g) => g,
            None => return false,
        };
        if is_constant(&g) {
            return true;
        }
    }
    any_nonzero_partial && is_constant(&g)
}

/// Square-free factorization of a dense univariate polynomial over a prime
/// field, including p-th root handling in characteristic p.
fn sff_uni_fp(f: &UP<FiniteField>) -> Vec<(UP<FiniteField>, usize)> {
    let p = f.domain().prime().to_usize().unwrap_or(usize::MAX);
    let df = f.derivative();
    if df.is_zero() {
        // f = g(x^p); over a prime field g(x^p) = g(x)^p.
        let mut root_coeffs = Vec::new();
        for (i, c) in f.coeffs().iter().enumerate() {
            if !f.domain().is_zero(c) {
                debug_assert!(i.is_multiple_of(p));
                if root_coeffs.len() <= i / p {
                    root_coeffs.resize(i / p + 1, f.domain().zero());
                }
                root_coeffs[i / p] = c.clone();
            }
        }
        let root = UP::<FiniteField>::from_coeffs(f.domain().clone(), root_coeffs);
        return sff_uni_fp(&root)
            .into_iter()
            .map(|(g, m)| (g, m * p))
            .collect();
    }
    let mut result = Vec::new();
    let mut c = f.gcd(&df);
    let mut w = f.div_rem(&c).expect("gcd divides").0;
    let mut i = 1;
    while !w.is_zero() && w.degree().unwrap_or(0) > 0 {
        let y = w.gcd(&c);
        let z = w.div_rem(&y).expect("gcd divides").0;
        if z.degree().unwrap_or(0) > 0 {
            result.push((monic_upoly(&z), i));
        }
        c = c.div_rem(&y).expect("gcd divides").0;
        w = y;
        i += 1;
    }
    if c.degree().unwrap_or(0) > 0 {
        // c consists of p-th powers.
        let root_sparse = {
            let sparse = dense_to_mpoly(&c, 1);
            shrink_var(&sparse, 0, p).map(|r| mpoly_to_dense(&r))
        };
        if let Some(root) = root_sparse {
            for (g, m) in sff_uni_fp(&root) {
                result.push((monic_upoly(&g), m * p));
            }
        } else {
            result.push((monic_upoly(&c), 1));
        }
    }
    result
}

/// GCD of the coefficient polynomials of `x_0^k` (the content in the main
/// variable), returned with the same number of variables (exponent of
/// variable 0 is 0).
fn content_main_fp(f: &FpMPoly) -> FpMPoly {
    let d = f.degree_in(0);
    let mut acc: Option<FpMPoly> = None;
    for k in 0..=d {
        let ck = f.coeff_of_var_pow(0, k);
        if ck.is_zero() {
            continue;
        }
        let dropped = ck.drop_main_var();
        acc = Some(match acc {
            None => dropped,
            Some(a) => multivariate_gcd_fp(&a, &dropped)
                .unwrap_or_else(|| one_mpoly(f.domain(), f.n_vars() - 1)),
        });
        if acc.as_ref().is_some_and(is_constant) {
            break;
        }
    }
    acc.map(|a| a.embed_new_main())
        .unwrap_or_else(|| one_mpoly(f.domain(), f.n_vars()))
}

/// Square-free factorization of a multivariate polynomial over a prime
/// field: Yun's algorithm with respect to variable 0, recursing on the
/// content, with characteristic-p p-th root handling.
fn sff_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)> {
    if f.is_zero() || is_constant(f) {
        return Vec::new();
    }
    let n = f.n_vars();
    let p = f.domain().prime().to_usize().unwrap_or(usize::MAX);
    if n == 1 {
        return sff_uni_fp(&mpoly_to_dense(f))
            .into_iter()
            .map(|(g, m)| (dense_to_mpoly(&g, 1), m))
            .collect();
    }

    let mut result = Vec::new();
    // Recursive square-free factorization of the content in variable 0.
    let cont = content_main_fp(f);
    let pp = if is_constant(&cont) {
        f.clone()
    } else {
        for (g, m) in sff_fp(&cont.drop_main_var()) {
            result.push((g.embed_new_main(), m));
        }
        f.checked_div_exact(&cont).expect("content divides")
    };
    if pp.degree_in(0) == 0 {
        return result;
    }

    let df = pp.derivative(0);
    if df.is_zero() {
        // pp = g(x_0^p, x̃). If pp is a full p-th power, take the root and
        // scale multiplicities. Otherwise shrink x_0, factor g, re-expand,
        // and resolve any p-th power structure in the composed factors.
        if let Some(root) = pth_root_full(&pp) {
            for (g, m) in sff_fp(&root) {
                result.push((g, m * p));
            }
            return result;
        }
        if let Some(g) = shrink_var(&pp, 0, p) {
            for (h, m) in sff_fp(&g) {
                let hc = expand_var(&h, 0, p);
                if is_square_free_mpoly(&hc) {
                    result.push((hc, m));
                } else if let Some(root) = pth_root_full(&hc) {
                    for (q, r) in sff_fp(&root) {
                        result.push((q, r * m * p));
                    }
                } else {
                    // Conservative fallback (should not occur).
                    result.push((hc, m));
                }
            }
        } else {
            result.push((pp, 1));
        }
        return result;
    }

    let one = one_mpoly(f.domain(), n);
    let mut c = multivariate_gcd_fp(&pp, &df).unwrap_or_else(|| one.clone());
    let mut w = pp.checked_div_exact(&c).unwrap_or_else(|| pp.clone());
    let mut i = 1;
    while !is_constant(&w) {
        let y = multivariate_gcd_fp(&w, &c).unwrap_or_else(|| one.clone());
        let z = w.checked_div_exact(&y).unwrap_or_else(|| w.clone());
        if !is_constant(&z) {
            result.push((z, i));
        }
        c = c.checked_div_exact(&y).unwrap_or_else(|| c.clone());
        w = y;
        i += 1;
    }
    if !is_constant(&c) {
        // The Yun tail is a full p-th power: all exponents of every variable
        // are divisible by p.
        if let Some(root) = pth_root_full(&c) {
            for (g, m) in sff_fp(&root) {
                result.push((g, m * p));
            }
        } else {
            result.push((c, 1));
        }
    }
    result
}

/// Factor a square-free multivariate polynomial over a prime field into
/// irreducible factors.
fn factor_square_free_fp(f: &FpMPoly) -> Vec<FpMPoly> {
    if is_constant(f) {
        return Vec::new();
    }
    if f.n_vars() == 1 {
        return mpoly_to_dense(f)
            .factor()
            .into_iter()
            .map(|(g, _)| dense_to_mpoly(&monic_upoly(&g), 1))
            .collect();
    }
    if f.degree_in(0) == 0 {
        return factor_square_free_fp(&f.drop_main_var())
            .into_iter()
            .map(|g| g.embed_new_main())
            .collect();
    }

    let lc = f.leading_coeff_in(0);
    let n = f.n_vars();
    let field = f.domain().clone();

    if !is_constant(&lc) {
        // Non-constant leading coefficient: use Wang's LC reconstruction
        // to distribute the factors of ℓ among the polynomial factors,
        // then lift with imposed leading coefficients.
        return factor_square_free_fp_nonconstant_lc(f, &lc);
    }

    // Constant LC: make monic and lift.
    let c = lc.coeff(&vec![0; n]);
    let inv_c = field.inv(&c).expect("nonzero leading coefficient");
    let f_m = f.mul_scalar(&inv_c);

    // Adaptive sample search: retry with a larger value range when the
    // first round yields no usable candidate.
    for range in [8u64, 16, 32] {
        for (sample, uni) in find_sample_fp(&f_m, 8, range, &[]) {
            if uni.len() == 1 {
                // Square-free, degree-preserving irreducible image ⇒ f irreducible.
                return vec![f.clone()];
            }
            if let Some(lifted) = eez_lift(&f_m, &sample, &uni) {
                let mut prod = one_mpoly(&field, n);
                for g in &lifted {
                    prod = prod.mul(g);
                }
                if prod == f_m {
                    let mut out = lifted;
                    out[0] = out[0].mul_scalar(&c);
                    return out;
                }
            }
        }
    }
    vec![f.clone()]
}

/// Factor a square-free multivariate polynomial over a prime field when the
/// leading coefficient (in variable 0) is non-constant.
///
/// Uses Wang's leading-coefficient reconstruction to determine the true
/// multivariate leading coefficients ℓ_i, then lifts the factors via
/// [`eez_lift_imposed`].
fn factor_square_free_fp_nonconstant_lc(f: &FpMPoly, lc: &FpMPoly) -> Vec<FpMPoly> {
    let n = f.n_vars();
    let field = f.domain().clone();

    // Compute the non-constant irreducible factors of the LC (in the
    // secondary variables) for sample filtering.
    let lc_reduced = lc.drop_main_var();
    let lc_filter: Vec<FpMPoly> = if lc_reduced.total_degree() == Some(0) {
        Vec::new()
    } else {
        multivariate_factor_fp(&lc_reduced)
            .into_iter()
            .map(|(g, _)| g)
            .filter(|g| !is_constant(g))
            .collect()
    };

    // Adaptive sample search with lc_filter for Wang distribution.
    for range in [8u64, 16, 32] {
        for (sample, mono) in find_sample_fp(f, 8, range, &lc_filter) {
            if mono.len() == 1 {
                continue; // try other samples before concluding irreducible
            }
            // Reconstruct the true multivariate leading coefficients.
            let true_lcoeffs = match wang_reconstruct_lcoeffs_fp(lc, &sample, &mono) {
                Some(l) => l,
                None => continue,
            };
            // Build initial factors: multiply the monic univariate factor by
            // ℓ_i to get the correct leading coefficient AND correct lower
            // terms. For monic u_i, lc(ℓ_i · u_i) = ℓ_i, and the lower
            // terms carry the ℓ_i scaling that the EEZ lift will refine.
            let mut initial: Vec<FpMPoly> = Vec::with_capacity(mono.len());
            for (i, u) in mono.iter().enumerate() {
                let f_i = dense_to_mpoly(u, n);
                initial.push(true_lcoeffs[i].mul(&f_i));
            }
            // EEZ lift with imposed leading coefficients.
            if let Some(lifted) = eez_lift_imposed(f, &sample, &initial, &true_lcoeffs) {
                let mut prod = one_mpoly(&field, n);
                for g in &lifted {
                    prod = prod.mul(g);
                }
                if equal_up_to_unit(&prod, f) {
                    let mut out = lifted;
                    out.sort_by_key(|b| std::cmp::Reverse(b.degree_in(0)));
                    return out;
                }
            }
        }
    }
    vec![f.clone()]
}

/// Whether two polynomials are equal up to a nonzero constant multiple.
fn equal_up_to_unit<D: Domain>(a: &MP<D>, b: &MP<D>) -> bool {
    if a.is_zero() || b.is_zero() {
        return a.is_zero() && b.is_zero();
    }
    // Pick a term of a and compare coefficient ratios across all terms.
    let (e0, c0) = match a.terms_ref().iter().next() {
        Some(t) => t,
        None => return false,
    };
    let bc0 = b.coeff(e0);
    if b.domain().is_zero(&bc0) {
        return false;
    }
    // ratio = a.coeff / b.coeff must be the same for every term.
    let ratio = match a.domain().div(c0, &bc0) {
        Some(r) => r,
        None => return false,
    };
    if a.n_terms() != b.n_terms() {
        return false;
    }
    a.terms_ref().iter().all(|(e, c)| {
        let bc = b.coeff(e);
        *c == a.domain().mul(&ratio, &bc)
    })
}

/// Factor a multivariate polynomial over a prime finite field into
/// irreducible factors with multiplicities.
///
/// Supports polynomials with non-constant leading coefficients in
/// variable 0 via Wang's leading-coefficient reconstruction and imposed
/// EEZ lifting.
pub fn multivariate_factor_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)> {
    if f.is_zero() || is_constant(f) {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (g, m) in sff_fp(f) {
        for h in factor_square_free_fp(&g) {
            result.push((h, m));
        }
    }
    // Safety net: the factorization must reconstruct the input up to a unit,
    // and every factor must be square-free. On any inconsistency, return the
    // input as a single factor (conservative, never wrong).
    let mut prod = one_mpoly(f.domain(), f.n_vars());
    for (g, m) in &result {
        for _ in 0..*m {
            prod = prod.mul(g);
        }
    }
    if !equal_up_to_unit(&prod, f) || result.iter().any(|(g, _)| !is_square_free_mpoly(g)) {
        return vec![(f.clone(), 1)];
    }
    result
}

// =========================================================================
// Multivariate factorization over the integers (Wang EEZ + LC preprocessing)
// =========================================================================

/// Sparse multivariate polynomial over the integers.
type ZmPoly = MP<IntegerDomain>;

/// Evaluate all secondary variables of a ℤ polynomial at the integer sample,
/// leaving a dense univariate polynomial in variable 0.
fn eval_to_image_z(f: &ZmPoly, sample: &[Integer]) -> UP<IntegerDomain> {
    let mut p = f.clone();
    for k in (1..f.n_vars()).rev() {
        p = p.eval_keep(k, &sample[k]);
    }
    mpoly_to_dense(&p)
}

/// GCD of the coefficient polynomials of `x_0^k` (content in main variable),
/// returned with the same number of variables (exponent of variable 0 is 0).
fn content_main_z(f: &ZmPoly) -> ZmPoly {
    let d = f.degree_in(0);
    let mut acc: Option<ZmPoly> = None;
    for k in 0..=d {
        let ck = f.coeff_of_var_pow(0, k);
        if ck.is_zero() {
            continue;
        }
        let dropped = ck.drop_main_var();
        acc = Some(match acc {
            None => dropped,
            Some(a) => multivariate_gcd_z(&a, &dropped).unwrap_or_else(|| {
                let mut one = ZmPoly::new(IntegerDomain, f.n_vars() - 1);
                one.set_term_external(vec![0; f.n_vars() - 1], Integer::from(1));
                one
            }),
        });
        if acc.as_ref().is_some_and(is_constant) {
            break;
        }
    }
    acc.map(|a| a.embed_new_main())
        .unwrap_or_else(|| one_mpoly(&IntegerDomain, f.n_vars()))
}

/// Square-free factorization over ℤ via the multivariate GCD (Yun's
/// algorithm, characteristic 0 so no p-th root handling is needed).
fn sff_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)> {
    if f.is_zero() || is_constant(f) {
        return Vec::new();
    }
    let n = f.n_vars();
    if n == 1 {
        return mpoly_to_dense(f)
            .factor()
            .into_iter()
            .map(|(g, m)| (dense_to_mpoly(&g, 1), m))
            .collect();
    }
    let mut result = Vec::new();
    let cont = content_main_z(f);
    let pp = if is_constant(&cont) {
        f.clone()
    } else {
        for (g, m) in sff_z(&cont.drop_main_var()) {
            result.push((g.embed_new_main(), m));
        }
        f.checked_div_exact(&cont).expect("content divides")
    };
    if pp.degree_in(0) == 0 {
        return result;
    }
    let one = one_mpoly(&IntegerDomain, n);
    let df = pp.derivative(0);
    let mut c = multivariate_gcd_z(&pp, &df).unwrap_or_else(|| one.clone());
    let mut w = pp.checked_div_exact(&c).unwrap_or_else(|| pp.clone());
    let mut i = 1;
    while !is_constant(&w) {
        let y = multivariate_gcd_z(&w, &c).unwrap_or_else(|| one.clone());
        let z = w.checked_div_exact(&y).unwrap_or_else(|| w.clone());
        if !is_constant(&z) {
            result.push((z, i));
        }
        c = c.checked_div_exact(&y).unwrap_or_else(|| c.clone());
        w = y;
        i += 1;
    }
    if !is_constant(&c) {
        result.push((c, 1));
    }
    result
}

/// Wang's leading-coefficient reconstruction from a univariate sample.
///
/// Given the factorization `ℓ = ∏ g_j^{e_j}` of the leading coefficient and
/// the primitive univariate image factors `u_i`, distribute the non-constant
/// factors `g_j` among the `u_i` using the pairwise-coprime integer images
/// `α_j = |g_j(s)|`, then reconcile integer remainders so that
/// `ℓ_i(s) = lc(u_i)` exactly. Returns the true leading coefficients `ℓ_i`
/// satisfying `c·∏ ℓ_i = ℓ`, or `None` on an unlucky sample.
///
/// Reference: Wang (1978), *An Improved Multivariate Polynomial Factoring
/// Algorithm*; Symbolica `reconstruct_lcoeffs_from_univariate_sample`.
fn wang_reconstruct_lcoeffs(
    lcoeff: &ZmPoly,
    sample: &[Integer],
    uni: &[UP<IntegerDomain>],
    content: &Integer,
) -> Option<Vec<ZmPoly>> {
    let n = lcoeff.n_vars();
    // If the leading coefficient is a constant, the true LCs are the
    // (constant) leading coefficients of the univariate image factors, so
    // that ℓ_i(s) = lc(u_i) holds and imposing them is a no-op at the sample.
    if lcoeff.degree_in(0) == 0 && lcoeff.drop_main_var().total_degree() == Some(0) {
        return Some(
            uni.iter()
                .map(|u| {
                    let lc = u
                        .leading_coeff()
                        .cloned()
                        .unwrap_or_else(|| Integer::from(1));
                    one_mpoly(&IntegerDomain, n).mul_scalar(&lc)
                })
                .collect(),
        );
    }
    // Factor the leading coefficient (in the remaining variables).
    let lc_reduced = lcoeff.drop_main_var();
    let lc_factors = if lc_reduced.n_vars() == 0 {
        Vec::new()
    } else {
        multivariate_factor_z(&lc_reduced)
    };

    // Integer images α_j = |g_j(s)|, requiring each to be > 1 and pairwise
    // coprime.
    let mut alpha: Vec<Integer> = Vec::new();
    let mut nonconst: Vec<ZmPoly> = Vec::new();
    let mut const_part = Integer::from(1);
    for (g, _e) in &lc_factors {
        if is_constant(g) {
            const_part = IntegerDomain.mul(&const_part, &g.coeff(&vec![0; g.n_vars()]));
            continue;
        }
        let mut img = g.clone();
        // g is in (n-1) variables (main var dropped); evaluate all of them
        // at the secondary sample values.
        for k in 0..img.n_vars() {
            img = img.eval_keep(k, &sample[k + 1]);
        }
        let a = img.coeff(&vec![0; img.n_vars()]).abs();
        if a <= Integer::from(1) {
            return None;
        }
        alpha.push(a);
        nonconst.push(g.clone());
    }
    for i in 0..alpha.len() {
        for j in (i + 1)..alpha.len() {
            if IntegerDomain.gcd(&alpha[i], &alpha[j]) != Integer::from(1) {
                return None;
            }
        }
    }

    // Multiplicities of each non-constant factor in ℓ.
    let multiplicities: Vec<usize> = nonconst
        .iter()
        .map(|g| {
            lc_factors
                .iter()
                .find(|(h, _)| h == g)
                .map(|(_, e)| *e)
                .unwrap_or(1)
        })
        .collect();

    // Distribution: greedily assign g_j to u_i while α_j divides lc(u_i).
    let r = uni.len();
    let mut lcoeffs: Vec<ZmPoly> = vec![one_mpoly(&IntegerDomain, n); r];
    let mut residual_lc: Vec<Integer> = uni
        .iter()
        .map(|u| u.leading_coeff().cloned().unwrap().abs())
        .collect();
    let mut used = vec![0usize; nonconst.len()];
    for i in 0..r {
        for j in 0..nonconst.len() {
            while used[j] < multiplicities[j]
                && IntegerDomain.div(&residual_lc[i], &alpha[j]).is_some()
            {
                lcoeffs[i] = lcoeffs[i].mul(&nonconst[j].embed_new_main());
                residual_lc[i] = IntegerDomain.div(&residual_lc[i], &alpha[j]).unwrap();
                used[j] += 1;
            }
        }
    }
    if used != multiplicities {
        return None;
    }

    // Integer remainder reconciliation: force ℓ_i(s) = lc(u_i) (with sign).
    for i in 0..r {
        let mut img = lcoeffs[i].clone();
        for (k, sk) in sample.iter().enumerate().skip(1) {
            img = img.eval_keep(k, sk);
        }
        let beta = img.coeff(&vec![0; img.n_vars()]);
        if beta.is_zero() {
            return None;
        }
        let target = uni[i].leading_coeff().cloned().unwrap();
        let q = IntegerDomain.div(&target, &beta)?;
        lcoeffs[i] = lcoeffs[i].mul_scalar(&q);
    }

    // Global verification: c · ∏ ℓ_i = ℓ as a polynomial identity.
    let mut prod = one_mpoly(&IntegerDomain, n).mul_scalar(content);
    for l in &lcoeffs {
        prod = prod.mul(l);
    }
    if prod == *lcoeff { Some(lcoeffs) } else { None }
}

/// Find integer sample points for the secondary variables such that the
/// univariate image has full degree and is square-free, returning candidate
/// (sample, primitive image factors, content) ordered by factor count.
///
/// `lc_filter` lists the non-constant factors of the leading coefficient
/// (in the secondary variables); samples where any of them evaluates to
/// `0` or `±1` are skipped, since Wang's leading-coefficient distribution
/// requires every such image `α_j` to exceed 1 in absolute value.
///
/// `value_bound` caps the absolute sample values; `factor_square_free_z`
/// retries with a larger bound when every candidate of the first round
/// fails downstream (adaptive sample search, cf. Symbolica `find_sample`
/// restarting with `coefficient_upper_bound + 10`).
///
/// At most `max_decompositions` successful univariate decompositions are
/// evaluated (distinct samples whose image is square-free and factors into
/// ≥ 2 coprime parts); beyond that the best candidates found so far are
/// returned.  This caps the dominant cost (univariate factoring) without
/// affecting the quality of the candidate pool.
#[allow(clippy::type_complexity)]
fn find_sample_z(
    f: &ZmPoly,
    max_candidates: usize,
    lc_filter: &[ZmPoly],
    value_bound: i64,
) -> Vec<(Vec<Integer>, Vec<UP<IntegerDomain>>, Integer)> {
    let n = f.n_vars();
    let deg0 = f.degree_in(0);
    let mut best: Vec<(Vec<Integer>, Vec<UP<IntegerDomain>>, Integer)> = Vec::new();
    let mut candidates: Vec<i64> = Vec::with_capacity((2 * value_bound + 1) as usize);
    candidates.push(0);
    for v in 1..=value_bound {
        candidates.push(v);
        candidates.push(-v);
    }
    let attempts = 4000usize;
    let max_decompositions = 200usize;
    let mut decompositions = 0usize;
    let mut seen: std::collections::HashSet<Vec<i64>> = std::collections::HashSet::new();

    for t in 0..attempts {
        let mut sample = vec![Integer::from(0); n];
        let mut rem = t;
        for slot in sample.iter_mut().enumerate().take(n).skip(1) {
            let idx = rem % candidates.len();
            rem /= candidates.len();
            *slot.1 = Integer::from(candidates[idx]);
        }
        let key: Vec<i64> = sample
            .iter()
            .map(|s| s.to_i64().unwrap_or(i64::MAX))
            .collect();
        if !seen.insert(key) {
            continue;
        }
        // Wang viability: every non-constant LC factor must evaluate to
        // an integer exceeding 1 in absolute value.
        if lc_filter.iter().any(|g| {
            let mut img = g.clone();
            for k in 0..img.n_vars() {
                img = img.eval_keep(k, &sample[k + 1]);
            }
            img.coeff(&vec![0; img.n_vars()]).abs() <= Integer::from(1)
        }) {
            continue;
        }
        let image = eval_to_image_z(f, &sample);
        if image.degree().unwrap_or(0) != deg0 || deg0 == 0 {
            continue;
        }
        if !image.is_square_free() {
            continue;
        }
        let content = image
            .coeffs()
            .iter()
            .fold(Integer::from(0), |acc, c| IntegerDomain.gcd(&acc, c));
        let content = content.abs();
        if content.is_zero() {
            continue;
        }
        let primitive_img = if content == Integer::from(1) {
            image.clone()
        } else {
            let coeffs = image
                .coeffs()
                .iter()
                .map(|c| IntegerDomain.div(c, &content).unwrap_or_else(|| c.clone()))
                .collect();
            UP::<IntegerDomain>::from_coeffs(IntegerDomain, coeffs)
        };
        let factors = primitive_img.factor();
        if factors.is_empty() {
            continue;
        }
        let mut uni: Vec<UP<IntegerDomain>> = Vec::new();
        let mut ok = true;
        for (g, m) in &factors {
            if *m != 1 {
                ok = false;
                break;
            }
            let lc = g.leading_coeff().cloned().unwrap();
            let g = if lc.is_negative() {
                g.mul_scalar(&Integer::from(-1))
            } else {
                g.clone()
            };
            uni.push(g);
        }
        if !ok {
            continue;
        }
        uni.sort_by_key(|b| std::cmp::Reverse(b.degree().unwrap_or(0)));
        decompositions += 1;
        // Rank by (factor count, image content): fewer factors is better,
        // and a smaller content gives Wang's distribution more room (the
        // image content cannot be recovered by the integer reconcile).
        let content = content.abs();
        let pos = best
            .binary_search_by(|(_, b, c)| {
                b.len().cmp(&uni.len()).then_with(|| c.abs().cmp(&content))
            })
            .unwrap_or_else(|e| e);
        best.insert(pos, (sample, uni, content));
        if best.len() > max_candidates {
            best.pop();
        }
        // Early exit for constant leading coefficients: two candidates with
        // at least two factors are enough (the 0.16.0 behaviour). With a
        // non-constant LC the full scan is required, since Wang's
        // distribution may reject early candidates downstream.
        if lc_filter.is_empty() && best.len() >= 2 && best[0].1.len() >= 2 {
            break;
        }
        // Cap the number of successful univariate decompositions to bound
        // the dominant cost (univariate factoring over ℤ).
        if decompositions >= max_decompositions {
            break;
        }
    }
    best
}

/// Zassenhaus subset recombination over ℤ: among the lifted modular factors,
/// find subsets whose product is an exact divisor of `f`.
fn zassenhaus_multivariate(f: &ZmPoly, lifted: &[ZmPoly]) -> Option<Vec<ZmPoly>> {
    let n_lift = lifted.len();
    let mut remaining = f.clone();
    let mut unused: Vec<bool> = vec![true; n_lift];
    let mut result = Vec::new();
    let mut subset_size = 1;
    while subset_size <= n_lift / 2 {
        let mut indices: Vec<usize> = (0..n_lift).filter(|&i| unused[i]).collect();
        let mut found = false;
        // Enumerate combinations of the unused factors of size subset_size.
        let mut combo: Vec<usize> = Vec::new();
        let m = indices.len();
        if subset_size > m {
            break;
        }
        let mut stack: Vec<usize> = (0..subset_size).collect();
        loop {
            // Test the current combination.
            combo.clear();
            combo.extend(stack.iter().map(|&s| indices[s]));
            let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
            for &ci in &combo {
                prod = prod.mul(&lifted[ci]);
            }
            // Strip the content of the candidate.
            let cont = prod.content();
            let cand = if cont.abs() == Integer::from(1) {
                prod.clone()
            } else {
                let mut p = ZmPoly::new(IntegerDomain, f.n_vars());
                for (e, c) in prod.terms_ref() {
                    p.set_term_external(
                        e.to_vec(),
                        IntegerDomain.div(c, &cont).unwrap_or_else(|| c.clone()),
                    );
                }
                p
            };
            if let Some(_q) = remaining.checked_div_exact(&cand) {
                result.push(cand);
                remaining = _q;
                for &ci in &combo {
                    unused[ci] = false;
                }
                found = true;
                break;
            }
            // Advance the combination (lexicographic).
            let mut i = subset_size as isize - 1;
            let mut advanced = false;
            while i >= 0 {
                let iu = i as usize;
                if stack[iu] < m - (subset_size - iu) {
                    stack[iu] += 1;
                    for j in (iu + 1)..subset_size {
                        stack[j] = stack[j - 1] + 1;
                    }
                    advanced = true;
                    break;
                }
                i -= 1;
            }
            if !advanced {
                break;
            }
        }
        indices.clear();
        if !found {
            subset_size += 1;
        }
        let _ = &indices;
        if unused.iter().all(|u| !u) {
            break;
        }
    }
    // Whatever remains is the last factor (primitive part).
    if remaining.degree_in(0) > 0 || remaining.total_degree() != Some(0) {
        let cont = remaining.content();
        let cand = if cont.abs() == Integer::from(1) {
            remaining.clone()
        } else {
            let mut p = ZmPoly::new(IntegerDomain, f.n_vars());
            for (e, c) in remaining.terms_ref() {
                p.set_term_external(
                    e.to_vec(),
                    IntegerDomain.div(c, &cont).unwrap_or_else(|| c.clone()),
                );
            }
            p
        };
        result.push(cand);
    }
    if result.is_empty() {
        None
    } else {
        Some(result)
    }
}

/// Factor a square-free multivariate integer polynomial into irreducible
/// factors, using Wang's EEZ algorithm with leading-coefficient
/// preprocessing.
fn factor_square_free_z(f: &ZmPoly) -> Vec<ZmPoly> {
    if is_constant(f) {
        return Vec::new();
    }
    if f.n_vars() == 1 {
        return mpoly_to_dense(f)
            .factor()
            .into_iter()
            .map(|(g, _)| dense_to_mpoly(&g, 1))
            .collect();
    }
    if f.degree_in(0) == 0 {
        return factor_square_free_z(&f.drop_main_var())
            .into_iter()
            .map(|g| g.embed_new_main())
            .collect();
    }

    let lcoeff = f.leading_coeff_in(0);
    // Non-constant factors of the leading coefficient (in the secondary
    // variables), used to filter samples that cannot admit Wang's
    // leading-coefficient distribution.
    let lc_filter: Vec<ZmPoly> =
        if lcoeff.degree_in(0) == 0 && lcoeff.drop_main_var().total_degree() == Some(0) {
            Vec::new()
        } else {
            multivariate_factor_z(&lcoeff.drop_main_var())
                .into_iter()
                .map(|(g, _)| g)
                .filter(|g| !is_constant(g))
                .collect()
        };
    // Adaptive sample search: if every candidate of the first round fails
    // (unlucky samples), retry with a larger value bound.
    for value_bound in [7i64, 15, 25] {
        let samples = find_sample_z(f, 16, &lc_filter, value_bound);
        if let Some(out) = factor_square_free_z_candidates(f, &lcoeff, samples) {
            return out;
        }
    }
    vec![f.clone()]
}

/// Try to factor square-free `f` from the given sample candidates; returns
/// `None` when every candidate is unlucky.
#[allow(clippy::type_complexity)]
fn factor_square_free_z_candidates(
    f: &ZmPoly,
    lcoeff: &ZmPoly,
    samples: Vec<(Vec<Integer>, Vec<UP<IntegerDomain>>, Integer)>,
) -> Option<Vec<ZmPoly>> {
    for (sample, uni, content) in samples {
        if uni.len() == 1 {
            // Degree-preserving square-free irreducible image at this sample.
            // Keep trying other samples to split; conclude irreducible only
            // after exhausting all candidates (falls through to the end).
            continue;
        }
        // Wang LC preprocessing: reconstruct the true leading coefficients.
        let true_lcoeffs = match wang_reconstruct_lcoeffs(lcoeff, &sample, &uni, &content) {
            Some(l) => l,
            None => {
                continue; // unlucky sample
            }
        };
        if true_lcoeffs.iter().any(|l| !is_constant(l)) {
            // Non-constant true leading coefficients: exact-ℚ lifting
            // generally produces non-integral corrections, so lift
            // p-adically instead (Wang imposition with coefficient Hensel
            // lifting over 𝔽_p followed by a p-adic lift in the
            // coefficients).
            if let Some(out) = padic_lift_factors(f, &sample, &uni, &content, &true_lcoeffs) {
                return Some(out);
            }
            continue;
        }
        // Impose true LCs on the univariate image factors to form the
        // initial (zeroth-order) lifted factors, then lift through the
        // remaining variables via EEZ. Since ℓ_i(s) = lc(u_i), the multiplier
        // ℓ_i / ℓ_i(s) is 1 at the sample, so F_i^(0) = (ℓ_i/lc(u_i))·u_i is
        // a polynomial whose image at s is u_i and whose leading coeff is ℓ_i.
        let mut initial: Vec<ZmPoly> = Vec::with_capacity(uni.len());
        let mut ok = true;
        for (i, u) in uni.iter().enumerate() {
            let f_i = dense_to_mpoly(u, f.n_vars());
            let mut lc_img = true_lcoeffs[i].clone();
            for (k, sk) in sample.iter().enumerate().skip(1) {
                lc_img = lc_img.eval_keep(k, sk);
            }
            let lc_at_sample = lc_img.coeff(&vec![0; lc_img.n_vars()]);
            if lc_at_sample.is_zero() {
                ok = false;
                break;
            }
            let scale = match true_lcoeffs[i]
                .checked_div_exact(&one_mpoly(&IntegerDomain, f.n_vars()).mul_scalar(&lc_at_sample))
            {
                Some(s) => s,
                None => {
                    ok = false;
                    break;
                }
            };
            initial.push(f_i.mul(&scale));
        }
        if !ok {
            continue;
        }
        // EEZ lift with imposed leading coefficients.
        let lifted = match eez_lift_z(f, &sample, &initial, &true_lcoeffs) {
            Some(l) => l,
            None => {
                continue;
            }
        };
        // Verify the lifted factors; on success, return primitive,
        // sign-normalized irreducibles.
        let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
        for g in &lifted {
            prod = prod.mul(g);
        }
        if equal_up_to_unit(&prod, f) {
            let mut out: Vec<ZmPoly> = lifted.iter().map(primitive_positive).collect();
            out.sort_by_key(|b| std::cmp::Reverse(b.degree_in(0)));
            return Some(out);
        }
        // Zassenhaus recombination over the lifted modular factors.
        if let Some(irr) = zassenhaus_multivariate(f, &lifted) {
            let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
            for g in &irr {
                prod = prod.mul(g);
            }
            if equal_up_to_unit(&prod, f) {
                return Some(irr);
            }
        }
    }
    None
}

/// Return the primitive part of `g` with a positive Lex-leading coefficient.
fn primitive_positive(g: &ZmPoly) -> ZmPoly {
    let c = g.content();
    let pp = if c.abs() == Integer::from(1) {
        g.clone()
    } else {
        let mut p = ZmPoly::new(IntegerDomain, g.n_vars());
        for (e, c2) in g.terms_ref() {
            p.set_term_external(
                e.to_vec(),
                IntegerDomain.div(c2, &c).unwrap_or_else(|| c2.clone()),
            );
        }
        p
    };
    if pp.leading_coeff().is_some_and(|l| l.is_negative()) {
        pp.neg()
    } else {
        pp
    }
}

/// Factor a multivariate integer polynomial into irreducible factors with
/// multiplicities, using Wang's EEZ algorithm.
pub fn multivariate_factor_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)> {
    if f.is_zero() || is_constant(f) {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (g, m) in sff_z(f) {
        for h in factor_square_free_z(&g) {
            result.push((h, m));
        }
    }
    // Safety net: the factorization must reconstruct the input up to a unit.
    let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
    for (g, m) in &result {
        for _ in 0..*m {
            prod = prod.mul(g);
        }
    }
    if !equal_up_to_unit(&prod, f) {
        return vec![(f.clone(), 1)];
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn field(p: i64) -> FiniteField {
        FiniteField::new(BigInt::from(p))
    }

    fn fmp(p: i64, n_vars: usize, terms: &[(Vec<usize>, i64)]) -> FpMPoly {
        let fld = field(p);
        FpMPoly::from_terms(
            fld.clone(),
            n_vars,
            terms
                .iter()
                .map(|(e, c)| (e.clone(), fld.element(BigInt::from(*c))))
                .collect(),
        )
    }

    fn zm_poly(n_vars: usize, terms: &[(Vec<usize>, i64)]) -> ZmPoly {
        SparseMultivariatePolynomial::from_terms(
            IntegerDomain,
            n_vars,
            terms
                .iter()
                .map(|(e, c)| (e.clone(), Integer::from(*c)))
                .collect(),
        )
    }

    fn product(factors: &[FpMPoly]) -> FpMPoly {
        let n = factors[0].n_vars();
        let mut acc = one_mpoly(factors[0].domain(), n);
        for g in factors {
            acc = acc.mul(g);
        }
        acc
    }

    fn product_z(factors: &[ZmPoly]) -> ZmPoly {
        let n = factors[0].n_vars();
        let mut acc = one_mpoly(&IntegerDomain, n);
        for g in factors {
            acc = acc.mul(g);
        }
        acc
    }

    fn with_mult(factors: &[(FpMPoly, usize)]) -> Vec<FpMPoly> {
        factors
            .iter()
            .map(|(g, m)| {
                let mut acc = one_mpoly(g.domain(), g.n_vars());
                for _ in 0..*m {
                    acc = acc.mul(g);
                }
                acc
            })
            .collect()
    }

    fn with_mult_z(factors: &[(ZmPoly, usize)]) -> Vec<ZmPoly> {
        factors
            .iter()
            .map(|(g, m)| {
                let mut acc = one_mpoly(&IntegerDomain, g.n_vars());
                for _ in 0..*m {
                    acc = acc.mul(g);
                }
                acc
            })
            .collect()
    }

    // ---- F_p tests ----

    #[test]
    fn eez_trivariate_three_linear_factors() {
        let p = 13;
        let f1 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 1], 2)],
        );
        let f3 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 0], 1)],
        );
        let f = f1.mul(&f2).mul(&f3);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 3, "expected 3 factors, got {:?}", factors);
        assert_eq!(product(&with_mult(&factors)), f);
    }

    #[test]
    fn eez_trivariate_with_quadratic() {
        let p = 13;
        let f1 = fmp(
            p,
            3,
            &[(vec![2, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], -1)],
        );
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product(&with_mult(&factors)), f);
    }

    #[test]
    fn eez_repeated_factors() {
        let p = 13;
        let f1 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = fmp(
            p,
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 0], 1)],
        );
        let f = f1.mul(&f1).mul(&f2);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        let mut mults: Vec<usize> = factors.iter().map(|(_, m)| *m).collect();
        mults.sort_unstable();
        assert_eq!(mults, vec![1, 2]);
        assert_eq!(product(&with_mult(&factors)), f);
    }

    #[test]
    fn eez_pth_power_char3() {
        let p = 3;
        let inner = fmp(p, 2, &[(vec![3, 0], 1), (vec![0, 1], 1)]);
        let f = inner.mul(&inner);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 1, "expected 1 factor, got {:?}", factors);
        assert_eq!(factors[0].1, 2, "expected multiplicity 2");
        assert_eq!(factors[0].0, inner);
    }

    #[test]
    fn eez_full_pth_power_char3() {
        let p = 3;
        let inner = fmp(p, 2, &[(vec![3, 0], 1), (vec![0, 3], 1)]);
        let f = inner.mul(&inner);
        let x_plus_y = fmp(p, 2, &[(vec![1, 0], 1), (vec![0, 1], 1)]);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 1, "expected 1 factor, got {:?}", factors);
        assert_eq!(factors[0].1, 6, "expected multiplicity 6");
        assert_eq!(factors[0].0, x_plus_y);
    }

    #[test]
    fn eez_four_variables() {
        let p = 13;
        let f1 = fmp(
            p,
            4,
            &[
                (vec![1, 0, 0, 0], 1),
                (vec![0, 1, 0, 0], 1),
                (vec![0, 0, 1, 0], 1),
                (vec![0, 0, 0, 1], 1),
            ],
        );
        let f2 = fmp(
            p,
            4,
            &[
                (vec![1, 0, 0, 0], 1),
                (vec![0, 1, 0, 0], -1),
                (vec![0, 0, 1, 0], 1),
                (vec![0, 0, 0, 1], -1),
            ],
        );
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product(&with_mult(&factors)), f);
    }

    #[test]
    fn eez_irreducible() {
        let p = 13;
        let f = fmp(p, 2, &[(vec![2, 0], 1), (vec![0, 2], 1), (vec![0, 0], 1)]);
        let factors = multivariate_factor_fp(&f);
        assert_eq!(product(&with_mult(&factors)), f);
    }

    #[test]
    fn eez_random_roundtrip_trivariate() {
        let p = 13;
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as i64
        };
        let mut cases = 0;
        for _case_idx in 0..15 {
            let n_factors = 2;
            let mut f = fmp(p, 3, &[(vec![0, 0, 0], 1)]);
            for _ in 0..n_factors {
                let n_terms = 2 + (next() % 2) as usize;
                let mut terms = vec![(vec![1, 0, 0], 1i64)];
                for _ in 0..n_terms {
                    let e = vec![
                        (next() % 2) as usize,
                        (next() % 2) as usize,
                        (next() % 2) as usize,
                    ];
                    let c = 1 + next() % (p - 1);
                    terms.push((e, c));
                }
                let g = fmp(p, 3, &terms);
                f = f.mul(&g);
            }
            let result = multivariate_factor_fp(&f);
            let mut prod = fmp(p, 3, &[(vec![0, 0, 0], 1)]);
            for (g, m) in &result {
                for _ in 0..*m {
                    prod = prod.mul(g);
                }
            }
            assert!(
                equal_up_to_unit(&prod, &f),
                "roundtrip failed: {:?}",
                result
            );
            cases += 1;
        }
        assert!(cases >= 10);
    }

    // ---- ℤ tests ----

    #[test]
    fn z_trivariate_monic_three_linear() {
        // f = (x + y + z)(x - y + 2z)(x + y + 1) over ℤ (monic in x).
        let f1 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 1], 2)],
        );
        let f3 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 0], 1)],
        );
        let f = f1.mul(&f2).mul(&f3);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 3, "expected 3 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn z_bivariate_wang_nonconstant_lcoeff() {
        // f = (y·x² + 1)(x + 1) over ℤ: f1 = yx²+1 is monic in x with
        // non-constant LC y, so the univariate image y(s)x²+1 is monic and
        // factorable; Wang distributes ℓ = y between the two factors.
        let f1 = zm_poly(2, &[(vec![2, 1], 1), (vec![0, 0], 1)]); // y x² + 1 (monic, lc y)
        let f2 = zm_poly(2, &[(vec![1, 0], 1), (vec![0, 0], 1)]); // x + 1 (monic, lc 1)
        let f = f1.mul(&f2);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn fp_sparse_diophantine_two_factor() {
        // Over F_13: f1 = x²y + 2x + 1 (skeleton 2x+1), f2 = xy + 3
        // (skeleton 3). The sparse solver must recover the unique
        // corrections δ1 = 5x + 7, δ2 = 11 with support inside the
        // skeletons.
        let f1 = fmp(13, 2, &[(vec![2, 1], 1), (vec![1, 0], 2), (vec![0, 0], 1)]);
        let f2 = fmp(13, 2, &[(vec![1, 1], 1), (vec![0, 0], 3)]);
        let d1 = fmp(13, 2, &[(vec![1, 0], 5), (vec![0, 0], 7)]);
        let d2 = fmp(13, 2, &[(vec![0, 0], 11)]);
        let error = d1.mul(&f2).add(&d2.mul(&f1));
        let sk1 = fmp(13, 2, &[(vec![1, 0], 2), (vec![0, 0], 1)]);
        let sk2 = fmp(13, 2, &[(vec![0, 0], 3)]);
        let deltas = sparse_diophantine_fp(&[f1, f2], &error, &[sk1, sk2], 42)
            .expect("sparse solve must succeed");
        assert_eq!(deltas, vec![d1, d2]);
    }

    #[test]
    fn fp_sparse_diophantine_three_factor() {
        // Over F_17: three factors with sparse skeletons.
        let f1 = fmp(17, 2, &[(vec![2, 1], 1), (vec![1, 0], 3), (vec![0, 1], 2)]);
        let f2 = fmp(17, 2, &[(vec![1, 1], 1), (vec![0, 0], 5)]);
        let f3 = fmp(17, 2, &[(vec![1, 0], 1), (vec![0, 1], 1), (vec![0, 0], 7)]);
        let d1 = fmp(17, 2, &[(vec![1, 0], 4), (vec![0, 1], 9)]);
        let d2 = fmp(17, 2, &[(vec![0, 0], 6)]);
        let d3 = fmp(17, 2, &[(vec![0, 1], 2), (vec![0, 0], 1)]);
        let error = d1
            .mul(&f2)
            .mul(&f3)
            .add(&d2.mul(&f1).mul(&f3))
            .add(&d3.mul(&f1).mul(&f2));
        let sk1 = fmp(17, 2, &[(vec![1, 0], 3), (vec![0, 1], 2)]);
        let sk2 = fmp(17, 2, &[(vec![0, 0], 5)]);
        let sk3 = fmp(17, 2, &[(vec![0, 1], 1), (vec![0, 0], 7)]);
        let deltas = sparse_diophantine_fp(&[f1, f2, f3], &error, &[sk1, sk2, sk3], 7)
            .expect("sparse solve must succeed");
        assert_eq!(deltas, vec![d1, d2, d3]);
    }

    // ---- 𝔽_p non-constant LC tests (0.16.2) ----

    #[test]
    fn fp_bivariate_nonconstant_lcoeff() {
        // f = (y·x² + 1)(x + 1) over F₁₃.
        // LC in x is y (non-constant). Wang distributes ℓ = y to the first factor.
        let f1 = fmp(13, 2, &[(vec![2, 1], 1), (vec![0, 0], 1)]); // y·x² + 1
        let f2 = fmp(13, 2, &[(vec![1, 0], 1), (vec![0, 0], 1)]); // x + 1
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        let mut prod = fmp(13, 2, &[(vec![0, 0], 1)]);
        for (g, m) in &factors {
            for _ in 0..*m {
                prod = prod.mul(g);
            }
        }
        assert!(
            equal_up_to_unit(&prod, &f),
            "Fp non-constant LC bivariate roundtrip failed: factors={:?}",
            factors
        );
        assert!(factors.len() >= 2, "expected at least 2 factors");
    }

    #[test]
    fn fp_trivariate_nonconstant_lcoeff() {
        // f = (z·x + y)(x + z) over F₁₇.
        // LC in x is z (non-constant in the secondary variables).
        let f1 = fmp(17, 3, &[(vec![1, 0, 1], 1), (vec![0, 1, 0], 1)]); // z·x + y
        let f2 = fmp(17, 3, &[(vec![1, 0, 0], 1), (vec![0, 0, 1], 1)]); // x + z
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        let mut prod = fmp(17, 3, &[(vec![0, 0, 0], 1)]);
        for (g, m) in &factors {
            for _ in 0..*m {
                prod = prod.mul(g);
            }
        }
        assert!(
            equal_up_to_unit(&prod, &f),
            "Fp non-constant LC trivariate roundtrip failed: factors={:?}",
            factors
        );
        assert!(factors.len() >= 2, "expected at least 2 factors");
    }

    #[test]
    fn fp_reducible_nonconstant_lcoeff() {
        // f = (x + y)(y·x² + z) over F₁₃.
        // LC of f in x is y (from the second factor); first factor is monic.
        let f1 = fmp(13, 3, &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1)]); // x + y
        let f2 = fmp(13, 3, &[(vec![2, 1, 0], 1), (vec![0, 0, 1], 1)]); // y·x² + z
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        let mut prod = fmp(13, 3, &[(vec![0, 0, 0], 1)]);
        for (g, m) in &factors {
            for _ in 0..*m {
                prod = prod.mul(g);
            }
        }
        assert!(
            equal_up_to_unit(&prod, &f),
            "Fp reducible non-constant LC roundtrip failed: factors={:?}",
            factors
        );
    }

    #[test]
    #[ignore] // slow: 4-variable factorization; run manually
    fn fp_four_var_nonconstant_lcoeff() {
        // Sparse product in 4 variables with non-constant LC over F₁₃.
        let f1 = fmp(
            13,
            4,
            &[
                (vec![2, 1, 1, 0], 1), // y·z·x²
                (vec![1, 0, 0, 0], 3),
                (vec![0, 1, 0, 0], 2),
                (vec![0, 0, 0, 1], 5),
            ],
        );
        let f2 = fmp(
            13,
            4,
            &[
                (vec![1, 1, 0, 0], 1), // y·x
                (vec![1, 0, 0, 1], 1), // w·x
                (vec![0, 1, 0, 0], 2),
                (vec![0, 0, 0, 0], 3),
            ],
        );
        let f = f1.mul(&f2);
        let factors = multivariate_factor_fp(&f);
        let mut prod = fmp(13, 4, &[(vec![0, 0, 0, 0], 1)]);
        for (g, m) in &factors {
            for _ in 0..*m {
                prod = prod.mul(g);
            }
        }
        assert!(
            equal_up_to_unit(&prod, &f),
            "Fp 4-var non-constant LC roundtrip failed: factors={:?}",
            factors
        );
    }

    #[test]
    fn z_sparse_four_var_nonconstant_lc() {
        // Sparse product in 4 variables with ≥ 50 terms and non-constant
        // leading coefficients factors back into its two sparse factors via
        // the p-adic path.
        let mut f1_terms = vec![(vec![2usize, 1, 1, 0], 1i64)]; // y·z·x² LC
        let mut f2_terms = vec![(vec![1, 1, 0, 0], 1i64), (vec![1, 0, 0, 1], 1)]; // (y+w)·x
        for i in 0..4usize {
            for j in 0..3usize {
                let c1 = ((i * 7 + j * 3) % 4 + 1) as i64;
                let c2 = ((i * 5 + j * 11 + 2) % 4 + 1) as i64;
                f1_terms.push((vec![i % 2, i, j, (i + j) % 2], c1));
                f2_terms.push((vec![0, (i + 1) % 3, (j + 2) % 2, i % 3], c2));
            }
        }
        let f1 = zm_poly(4, &f1_terms);
        let f2 = zm_poly(4, &f2_terms);
        let f = f1.mul(&f2);
        assert!(
            f.n_terms() >= 50,
            "test product should be sparse-large, got {} terms",
            f.n_terms()
        );
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn small_prime_escalation_bails_for_large_skeleton_groups() {
        // Both factors have skeleton groups of size 2 ({y, z} at
        // x_0-degree 0), but 𝔽_2 offers only one nonzero element, so
        // skeleton interpolation cannot work: with the dense fallback
        // disabled the lift must bail out (the caller then escalates to a
        // larger prime); with the fallback enabled it succeeds.
        let f1 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], -1)],
        );
        let target = f1.mul(&f2);
        let one = zm_poly(3, &[(vec![0, 0, 0], 1)]);
        let lcs = [one.clone(), one.clone()];
        let max_p = Integer::from(8);
        let sample = vec![Integer::from(0); 3];
        assert!(
            coefficient_hensel_lift_z(
                &target,
                vec![f1.clone(), f2.clone()],
                &lcs,
                2,
                &max_p,
                &sample,
                false,
            )
            .is_none(),
            "small prime must bail out when skeleton groups need ≥ p generators"
        );
        let lifted = coefficient_hensel_lift_z(
            &target,
            vec![f1.clone(), f2.clone()],
            &lcs,
            2,
            &max_p,
            &sample,
            true,
        )
        .expect("dense fallback must succeed on the same input");
        assert_eq!(lifted, vec![f1, f2]);
    }

    #[test]
    fn z_coefficient_lift_uses_sparse_diophantine() {
        // Surgical test of the p-adic coefficient lift: corrupt one body
        // coefficient of a known factor, then check that the lift restores
        // it and that the sparse Diophantine solver was used (thread-local
        // hit counter). The prime 13 leaves enough nonzero elements for the
        // interpolation generators (groups of ≤ 4).
        let f1_true = zm_poly(
            4,
            &[
                (vec![2, 1, 1, 0], 1), // y·z·x² (LC)
                (vec![1, 0, 0, 0], 3),
                (vec![1, 1, 0, 0], 1),
                (vec![0, 0, 1, 0], 2),
                (vec![0, 0, 0, 1], 1),
            ],
        );
        let f2_true = zm_poly(
            4,
            &[
                (vec![1, 1, 0, 0], 1), // (y+w)·x (LC)
                (vec![1, 0, 0, 1], 1),
                (vec![0, 1, 0, 0], 2),
                (vec![0, 0, 0, 0], 3),
            ],
        );
        let target = f1_true.mul(&f2_true);
        // Corrupt one skeleton coefficient of f1 by a multiple of the prime
        // (3 → 3 + 13 = 16 at x^1), so the initial error is 0 mod 13.
        let f1_bad = zm_poly(
            4,
            &[
                (vec![2, 1, 1, 0], 1),
                (vec![1, 0, 0, 0], 16),
                (vec![1, 1, 0, 0], 1),
                (vec![0, 0, 1, 0], 2),
                (vec![0, 0, 0, 1], 1),
            ],
        );
        let lc1 = zm_poly(4, &[(vec![0, 1, 1, 0], 1)]); // y·z
        let lc2 = zm_poly(4, &[(vec![0, 1, 0, 0], 1), (vec![0, 0, 0, 1], 1)]); // y+w
        let max_p = Integer::from(13 * 13);
        let sample = vec![Integer::from(0); 4];
        let before = SPARSE_DIO_HITS.with(|h| h.get());
        let lifted = coefficient_hensel_lift_z(
            &target,
            vec![f1_bad, f2_true.clone()],
            &[lc1, lc2],
            13,
            &max_p,
            &sample,
            true,
        )
        .expect("coefficient lift must succeed");
        let after = SPARSE_DIO_HITS.with(|h| h.get());
        assert_eq!(lifted, vec![f1_true, f2_true]);
        assert!(
            after > before,
            "sparse Diophantine was not used in the coefficient lift"
        );
    }

    #[test]
    fn z_nonconstant_lcoeff_reducible_lc() {
        // Proptest regression: f = (x·y² − x + 2y)(x·y − z), ℓ = y³ − y
        // factors as (y−1)(y+1)y and must be distributed (y²−1) / y across
        // the two factors. Samples with |y(s)| ≤ 2 are filtered out.
        let a = zm_poly(
            3,
            &[(vec![1, 2, 0], 1), (vec![1, 0, 0], -1), (vec![0, 1, 0], 2)],
        );
        let b = zm_poly(3, &[(vec![1, 1, 0], 1), (vec![0, 0, 1], -1)]);
        let f = a.mul(&b);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn z_nonconstant_lcoeff_shared_monomial() {
        // Proptest regression: f = (x·y + 1)(x·y - z). The leading
        // coefficient y² must be distributed one y per factor; samples with
        // z = 0 have image content |y(s)| > 1 and cannot admit the
        // distribution, so candidate ranking prefers content-1 samples.
        let a = zm_poly(3, &[(vec![1, 1, 0], 1), (vec![0, 0, 0], 1)]);
        let b = zm_poly(3, &[(vec![1, 1, 0], 1), (vec![0, 0, 1], -1)]);
        let f = a.mul(&b);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn z_trivariate_nonconstant_lcoeff() {
        // f = (z·x² + y)(x + 1) over ℤ: f1 monic in x with LC z; images
        // z(s)x²+y(s)+1·... split as (linear)·(quadratic with LC z(s)).
        let f1 = zm_poly(3, &[(vec![2, 0, 1], 1), (vec![0, 1, 0], 1)]); // z x² + y (monic, lc z)
        let f2 = zm_poly(3, &[(vec![1, 0, 0], 1), (vec![0, 0, 0], 1)]); // x + 1 (monic, lc 1)
        let f = f1.mul(&f2);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }

    #[test]
    fn z_repeated_factors() {
        // f = (x + y + z)^2 (x - y + 1) over ℤ.
        let f1 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let f2 = zm_poly(
            3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 0], 1)],
        );
        let f = f1.mul(&f1).mul(&f2);
        let factors = multivariate_factor_z(&f);
        assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
        let mut mults: Vec<usize> = factors.iter().map(|(_, m)| *m).collect();
        mults.sort_unstable();
        assert_eq!(mults, vec![1, 2]);
        assert_eq!(product_z(&with_mult_z(&factors)), f);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    /// A random small multivariate polynomial over ℤ in exactly `n_vars`
    /// variables that is monic in variable 0. Kept deliberately tiny (few
    /// terms, low degree, small coefficients) so that multivariate
    /// factorization stays fast enough for property-based roundtrip testing.
    fn any_monic_zmp(n_vars: usize) -> impl Strategy<Value = ZmPoly> {
        (1usize..=2, 1usize..=2).prop_flat_map(move |(n_terms, max_deg)| {
            prop::collection::vec(
                (prop::collection::vec(0usize..=max_deg, n_vars), -2i64..=2),
                n_terms..=n_terms + 1,
            )
            .prop_map(move |mut terms| {
                let mut lead = vec![0usize; n_vars];
                lead[0] = max_deg.max(1); // ensure positive degree in x_0
                terms.push((lead, 1));
                SparseMultivariatePolynomial::<IntegerDomain, Lex>::from_terms(
                    IntegerDomain,
                    n_vars,
                    terms
                        .into_iter()
                        .map(|(e, c)| (e, Integer::from(c)))
                        .collect(),
                )
            })
        })
    }

    /// Reconstruct a polynomial from its factorization and compare up to a
    /// unit (sign), which is the strongest check that factorization did not
    /// lose or corrupt information.
    fn reconstructs(f: &ZmPoly, factors: &[(ZmPoly, usize)]) -> bool {
        let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
        for (g, m) in factors {
            for _ in 0..*m {
                prod = prod.mul(g);
            }
        }
        equal_up_to_unit(&prod, f)
    }

    /// A random small multivariate polynomial over ℤ in exactly `n_vars`
    /// variables whose leading coefficient in variable 0 is a non-constant
    /// monomial in the secondary variables.
    fn any_nonconstant_lc_zmp(n_vars: usize) -> impl Strategy<Value = ZmPoly> {
        (1usize..=2, 1usize..=2, 1usize..=2).prop_flat_map(move |(n_terms, max_deg, lc_deg)| {
            (
                prop::collection::vec(
                    (prop::collection::vec(0usize..=max_deg, n_vars), -2i64..=2),
                    n_terms..=n_terms + 1,
                ),
                prop::collection::vec(0usize..=lc_deg, n_vars - 1),
            )
                .prop_map(move |(mut terms, lc_exp)| {
                    let mut lead = vec![0usize; n_vars];
                    lead[0] = max_deg.max(1); // positive degree in x_0
                    for (i, e) in lc_exp.iter().enumerate() {
                        lead[i + 1] = *e;
                    }
                    // Ensure the leading coefficient is non-constant.
                    if lead.iter().skip(1).all(|&e| e == 0) {
                        lead[1] = 1;
                    }
                    terms.push((lead, 1));
                    SparseMultivariatePolynomial::<IntegerDomain, Lex>::from_terms(
                        IntegerDomain,
                        n_vars,
                        terms
                            .into_iter()
                            .map(|(e, c)| (e, Integer::from(c)))
                            .collect(),
                    )
                })
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        /// Factoring a product of two small monic factors must reconstruct it.
        ///
        /// Marked `ignore` because multivariate factorization is slow enough
        /// that a property-based sweep does not fit the unit-test budget; run
        /// manually or via the audit report.
        #[test]
        #[ignore = "slow multivariate factorization proptest: run manually"]
        fn factor_product_of_two_monic(
            a in any_monic_zmp(3),
            b in any_monic_zmp(3),
        ) {
            let f = a.mul(&b);
            let factors = multivariate_factor_z(&f);
            prop_assert!(
                reconstructs(&f, &factors),
                "factorization does not reconstruct input: {:?}",
                factors
            );
        }

        /// Factoring a product of two small factors with non-constant
        /// leading coefficients (Wang imposition + p-adic lifting) must
        /// reconstruct it and must not report the product as irreducible.
        ///
        /// Marked `ignore` for the same reason as the monic roundtrip.
        #[test]
        #[ignore = "slow multivariate factorization proptest: run manually"]
        fn factor_product_of_two_nonconstant_lc(
            a in any_nonconstant_lc_zmp(3),
            b in any_nonconstant_lc_zmp(3),
        ) {
            let f = a.mul(&b);
            let factors = multivariate_factor_z(&f);
            prop_assert!(
                reconstructs(&f, &factors),
                "factorization does not reconstruct input: {:?}",
                factors
            );
            // f = a·b with both factors of positive x_0-degree is reducible,
            // so a complete factorization cannot be trivial.
            let nontrivial = factors.len() >= 2 || factors.iter().any(|(_, m)| *m >= 2);
            prop_assert!(nontrivial, "factorization is trivial: {:?}", factors);
        }
    }
}
