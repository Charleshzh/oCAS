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
            let sigmas = diophantine(&lifted, &e, sample, k)?;
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
            let sigmas = diophantine_z(&lifted, &e, sample, k)?;
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
fn find_sample_fp(
    f: &FpMPoly,
    max_candidates: usize,
) -> Vec<(Vec<FiniteFieldElement>, Vec<UP<FiniteField>>)> {
    let n = f.n_vars();
    let field = f.domain().clone();
    let p = field.prime().to_u64().unwrap_or(u64::MAX);
    let range = p.clamp(1, 8);
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
    if !is_constant(&lc) {
        // Wang leading-coefficient preprocessing lands in a later phase;
        // conservatively report the input as unfactored.
        return vec![f.clone()];
    }
    let c = lc.coeff(&vec![0; f.n_vars()]);
    let inv_c = f.domain().inv(&c).expect("nonzero leading coefficient");
    let f_m = f.mul_scalar(&inv_c);

    for (sample, uni) in find_sample_fp(&f_m, 8) {
        if uni.len() == 1 {
            // Square-free, degree-preserving irreducible image ⇒ f irreducible.
            return vec![f.clone()];
        }
        if let Some(lifted) = eez_lift(&f_m, &sample, &uni) {
            let mut prod = one_mpoly(f.domain(), f.n_vars());
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
/// Currently the leading coefficient in variable 0 must be a nonzero field
/// constant; otherwise the (square-free part of the) input is returned
/// unfactored.
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
#[allow(clippy::type_complexity)]
fn find_sample_z(
    f: &ZmPoly,
    max_candidates: usize,
) -> Vec<(Vec<Integer>, Vec<UP<IntegerDomain>>, Integer)> {
    let n = f.n_vars();
    let deg0 = f.degree_in(0);
    let mut best: Vec<(Vec<Integer>, Vec<UP<IntegerDomain>>, Integer)> = Vec::new();
    let candidates: [i64; 15] = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7];
    let attempts = 4000usize;

    for t in 0..attempts {
        let mut sample = vec![Integer::from(0); n];
        let mut rem = t;
        for (k, slot) in sample.iter_mut().enumerate().take(n).skip(1) {
            let idx = rem % candidates.len();
            rem /= candidates.len();
            *slot = Integer::from(candidates[idx]);
            let _ = k;
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
        let pos = best
            .binary_search_by(|(_, b, _)| b.len().cmp(&uni.len()))
            .unwrap_or_else(|e| e);
        best.insert(pos, (sample, uni, content));
        if best.len() > max_candidates {
            best.pop();
        }
        // Early exit: a sample with ≥2 factors is good enough.
        if !best.is_empty() && best[0].1.len() >= 2 && best.len() >= 2 {
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
    let samples = find_sample_z(f, 8);
    for (sample, uni, content) in samples {
        if uni.len() == 1 {
            // Degree-preserving square-free irreducible image at this sample.
            // Keep trying other samples to split; conclude irreducible only
            // after exhausting all candidates (falls through to the end).
            continue;
        }
        // Wang LC preprocessing: reconstruct the true leading coefficients.
        let true_lcoeffs = match wang_reconstruct_lcoeffs(&lcoeff, &sample, &uni, &content) {
            Some(l) => l,
            None => {
                continue; // unlucky sample
            }
        };
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
            return out;
        }
        // Zassenhaus recombination over the lifted modular factors.
        if let Some(irr) = zassenhaus_multivariate(f, &lifted) {
            let mut prod = one_mpoly(&IntegerDomain, f.n_vars());
            for g in &irr {
                prod = prod.mul(g);
            }
            if equal_up_to_unit(&prod, f) {
                return irr;
            }
        }
    }
    vec![f.clone()]
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
    #[ignore = "non-constant LC imposition requires mod-p Hensel lifting (0.16.1)"]
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
    #[ignore = "non-constant LC imposition requires mod-p Hensel lifting (0.16.1)"]
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
