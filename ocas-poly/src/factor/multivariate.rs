//! Multivariate polynomial factorization.
//!
//! Currently implements bivariate factorization over $\mathbb{Z}$ and
//! $\mathbb{F}_p$ via evaluation and Hensel lifting (Wang's algorithm).
//!
//! The bivariate polynomial is treated as a univariate polynomial in the main
//! variable $x$ with coefficients in $\mathbb{Z}[y]$ (or $\mathbb{F}_p[y]$).
//! It is evaluated at $y = \alpha$ to obtain a univariate image over the base
//! domain, factored there, and the factors are lifted back to bivariate
//! polynomials by linear Hensel lifting in the ideal $(y - \alpha)$.
//!
//! References: Wang (1978), "An Improved Multivariate Polynomial Factoring
//! Algorithm"; Geddes, Czapor, Labahn, *Algorithms for Computer Algebra*.

use num_bigint::BigInt;
use num_traits::{One, Signed};
use ocas_domain::{
    Domain, FiniteField, FiniteFieldElement, Integer, IntegerDomain, Rational, RationalDomain,
};

use crate::dense::DenseUnivariatePolynomial;
use crate::factor::hensel;
use crate::sparse::{Lex, MonomialOrder, SparseMultivariatePolynomial};

/// Bivariate polynomial over the integers with lexicographic order.
pub type ZMPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

/// Univariate polynomial over the rationals.
pub type QPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Univariate polynomial over the integers.
pub type ZPoly = DenseUnivariatePolynomial<IntegerDomain>;

/// Univariate polynomial over a prime finite field.
pub type FpPoly = DenseUnivariatePolynomial<FiniteField>;

/// Bivariate polynomial over a prime finite field with lexicographic order.
pub type FpMPoly = SparseMultivariatePolynomial<FiniteField, Lex>;

/// Return the maximum degree of `poly` in variable `var_index`, or `0` for the
/// zero polynomial.
fn degree_in_var<D: Domain, O: MonomialOrder>(
    poly: &SparseMultivariatePolynomial<D, O>,
    var_index: usize,
) -> usize {
    poly.terms_ref()
        .keys()
        .map(|e| e.get(var_index).copied().unwrap_or(0))
        .max()
        .unwrap_or(0)
}

/// Evaluate `poly` at `y = value` and interpret the result as a univariate
/// polynomial in the main variable `x` (variable 0).
fn eval_to_univariate(poly: &ZMPoly, y_var: usize, value: &Integer) -> ZPoly {
    let evaluated = poly.eval(y_var, value);
    let mut coeffs = Vec::new();
    for (exp, c) in evaluated.terms_ref() {
        let idx = exp.first().copied().unwrap_or(0);
        if idx >= coeffs.len() {
            coeffs.resize(idx + 1, IntegerDomain.zero());
        }
        coeffs[idx] = c.clone();
    }
    ZPoly::from_coeffs(IntegerDomain, coeffs)
}

/// Lift a univariate polynomial in `x` back to a bivariate polynomial with
/// no dependence on the secondary variable `y`.
fn univariate_to_bivariate(g: &ZPoly, n_vars: usize, x_var: usize) -> ZMPoly {
    let mut terms = Vec::new();
    for (i, c) in g.coeffs().iter().enumerate() {
        if !IntegerDomain.is_zero(c) {
            let mut exp = vec![0usize; n_vars];
            exp[x_var] = i;
            terms.push((exp, c.clone()));
        }
    }
    ZMPoly::from_terms(IntegerDomain, n_vars, terms)
}

/// Multiply a univariate polynomial in `x` by $(y - \alpha)^k$ and return the
/// bivariate result. This is used to add a Hensel correction term.
fn univariate_times_y_minus_alpha_k(
    g: &ZPoly,
    k: usize,
    alpha: &Integer,
    n_vars: usize,
    x_var: usize,
    y_var: usize,
) -> ZMPoly {
    let mut terms = Vec::new();
    for (i, c) in g.coeffs().iter().enumerate() {
        if IntegerDomain.is_zero(c) {
            continue;
        }
        for j in 0..=k {
            let mut exp = vec![0usize; n_vars];
            exp[x_var] = i;
            exp[y_var] = j;
            let sign = if (k - j) % 2 == 0 { 1 } else { -1 };
            let coeff = Integer::from(
                BigInt::from(binomial(k, j)) * sign * alpha.inner().pow((k - j) as u32),
            );
            let prod = IntegerDomain.mul(c, &coeff);
            terms.push((exp, prod));
        }
    }
    ZMPoly::from_terms(IntegerDomain, n_vars, terms)
}

/// Binomial coefficient $\binom{n}{k}$.
fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    let k = k.min(n - k);
    let mut num = 1u64;
    let mut den = 1u64;
    for i in 0..k {
        num *= (n - i) as u64;
        den *= (i + 1) as u64;
    }
    num / den
}

/// Take the partial derivative of `poly` with respect to `var_index`.
fn derivative_in_var(poly: &ZMPoly, var_index: usize) -> ZMPoly {
    let mut result = ZMPoly::new(IntegerDomain, poly.n_vars());
    for (exp, coeff) in poly.terms_ref() {
        let power = exp.get(var_index).copied().unwrap_or(0);
        if power == 0 {
            continue;
        }
        let mut new_exp = exp.to_vec();
        new_exp[var_index] = power - 1;
        let scalar = IntegerDomain.cast_u64(power as u64);
        let new_coeff = IntegerDomain.mul(coeff, &scalar);
        result.set_term_external(new_exp, new_coeff);
    }
    result
}

/// Compute the Taylor coefficients of `poly` viewed as a polynomial in
/// $(y - \alpha)$, up to degree `max_k`. The coefficients are univariate
/// polynomials in `x`.
fn taylor_coeffs_in_y(poly: &ZMPoly, y_var: usize, alpha: &Integer, max_k: usize) -> Vec<ZPoly> {
    let mut coeffs = Vec::with_capacity(max_k + 1);
    let mut current = poly.clone();
    for k in 0..=max_k {
        let value = eval_to_univariate(&current, y_var, alpha);
        coeffs.push(divide_by_k_factorial(value, k));
        current = derivative_in_var(&current, y_var);
    }
    coeffs
}

/// Divide every coefficient of a univariate integer polynomial by `k!`.
fn divide_by_k_factorial(poly: ZPoly, k: usize) -> ZPoly {
    let mut fact = BigInt::one();
    for i in 1..=k {
        fact *= BigInt::from(i);
    }
    let fact_int = Integer::from(fact);
    let coeffs = poly
        .coeffs()
        .iter()
        .map(|c| IntegerDomain.div(c, &fact_int).unwrap_or_else(|| c.clone()))
        .collect();
    ZPoly::from_coeffs(IntegerDomain, coeffs)
}

/// Normalize a univariate integer polynomial to be monic by adjusting the sign
/// if necessary. The input is assumed to be primitive.
fn monic_zpoly(f: &ZPoly) -> ZPoly {
    if f.is_zero() {
        return f.clone();
    }
    let lc = f.leading_coeff().cloned().unwrap();
    if lc.inner().is_negative() {
        f.mul_scalar(&Integer::from(-1))
    } else {
        f.clone()
    }
}

/// Factor a primitive univariate integer polynomial into monic irreducible
/// factors.
fn factor_univariate_z(f: &ZPoly) -> Vec<ZPoly> {
    hensel::factor_primitive(f)
        .into_iter()
        .map(|(g, _)| monic_zpoly(&g))
        .collect()
}

/// Lift a square-free bivariate integer polynomial $f$ from its univariate
/// factorization at $y = \alpha$ back to a bivariate factorization, assuming
/// $f$ is monic in $x$ (its leading coefficient in $x$ is an integer constant).
fn hensel_lift_bivariate(
    f: &ZMPoly,
    alpha: &Integer,
    univariate_factors: &[ZPoly],
    x_var: usize,
    y_var: usize,
) -> Option<Vec<ZMPoly>> {
    let n_vars = f.n_vars();
    let d_y = degree_in_var(f, y_var);

    let c_f = taylor_coeffs_in_y(f, y_var, alpha, d_y);
    let q_factors: Vec<QPoly> = univariate_factors.iter().map(zpoly_to_qpoly).collect();
    let bezout_q = bezout_coefficients_q(&q_factors);

    let mut lifted: Vec<ZMPoly> = univariate_factors
        .iter()
        .map(|g| univariate_to_bivariate(g, n_vars, x_var))
        .collect();

    for k in 1..=d_y {
        let mut product = ZMPoly::from_terms(
            IntegerDomain,
            n_vars,
            vec![(vec![0; n_vars], Integer::from(1))],
        );
        for g in &lifted {
            product = product.mul(g);
        }

        let c_product = taylor_coeffs_in_y(&product, y_var, alpha, d_y);
        let error = c_f[k].sub(&c_product[k]);
        let error_q = zpoly_to_qpoly(&error);

        for i in 0..lifted.len() {
            let delta_q = error_q.mul(&bezout_q[i]);
            let (_q, remainder_q) = delta_q.div_rem(&q_factors[i]).unwrap();
            let remainder_z = qpoly_to_zpoly(&remainder_q)?;
            let correction =
                univariate_times_y_minus_alpha_k(&remainder_z, k, alpha, n_vars, x_var, y_var);
            lifted[i] = lifted[i].add(&correction);
        }
    }

    Some(lifted)
}

/// Choose evaluation points $y = \alpha$ such that the univariate image
/// $f(x, \alpha)$ is square-free and has the fewest irreducible factors.
/// A "lucky" point for Wang's Hensel lifting must preserve the bivariate
/// factorization pattern, so fewer factors are preferred over more.
///
/// Returns all candidates with at least two factors, ordered from fewest to
/// most factors, so the caller can retry on unlucky points.
fn choose_evaluation_points(f: &ZMPoly, y_var: usize) -> Vec<(Integer, Vec<ZPoly>)> {
    let candidates: [i64; 11] = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5];
    let mut best: Vec<(Integer, Vec<ZPoly>)> = Vec::new();
    for alpha in candidates {
        let alpha_int = Integer::from(alpha);
        let image = eval_to_univariate(f, y_var, &alpha_int);
        if image.degree().unwrap_or(0) < 1 || !image.is_square_free() {
            continue;
        }
        let factors = factor_univariate_z(&image);
        if factors.len() < 2 {
            continue;
        }
        // Keep candidates ordered by increasing number of factors.
        let insert_pos = best
            .binary_search_by(|(_, b)| b.len().cmp(&factors.len()))
            .unwrap_or_else(|e| e);
        best.insert(insert_pos, (alpha_int, factors));
    }
    best
}

/// Check whether a bivariate polynomial is the constant 1.
fn is_one_mpoly(f: &ZMPoly) -> bool {
    f.terms_ref().len() == 1
        && f.terms_ref()
            .iter()
            .next()
            .map(|(e, c)| e.iter().all(|&p| p == 0) && IntegerDomain.is_one(c))
            .unwrap_or(false)
}

/// Factor a square-free bivariate integer polynomial that is monic in $x$.
fn bivariate_factor_square_free_monic(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<ZMPoly> {
    if degree_in_var(f, x_var) == 0 || degree_in_var(f, y_var) == 0 {
        return vec![f.clone()];
    }

    let candidates = choose_evaluation_points(f, y_var);
    if candidates.is_empty() {
        return vec![f.clone()];
    }

    for (alpha, mut univariate_factors) in candidates {
        if univariate_factors.len() <= 1 {
            continue;
        }

        univariate_factors.sort_by_key(|b| std::cmp::Reverse(b.degree().unwrap_or(0)));

        let lifted = match hensel_lift_bivariate(f, &alpha, &univariate_factors, x_var, y_var) {
            Some(v) => v,
            None => continue,
        };

        let mut product = ZMPoly::from_terms(
            IntegerDomain,
            f.n_vars(),
            vec![(vec![0; f.n_vars()], Integer::from(1))],
        );
        for g in &lifted {
            product = product.mul(g);
        }
        if product == f.clone() || product == f.neg() {
            return lifted;
        }
    }

    vec![f.clone()]
}

/// Check whether the leading coefficient of `f` in variable `x_var` is a
/// nonzero integer constant (i.e., independent of all other variables).
fn lc_x_is_constant(f: &ZMPoly, x_var: usize) -> bool {
    let deg_x = degree_in_var(f, x_var);
    if deg_x == 0 {
        return true;
    }
    for exp in f.terms_ref().keys() {
        if exp.get(x_var).copied().unwrap_or(0) == deg_x {
            for (i, &e) in exp.iter().enumerate() {
                if i != x_var && e != 0 {
                    return false;
                }
            }
        }
    }
    true
}

/// Factor a primitive, square-free bivariate integer polynomial into irreducible
/// factors. Currently requires the leading coefficient in `x` to be a constant.
fn bivariate_factor_square_free(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<ZMPoly> {
    if !lc_x_is_constant(f, x_var) {
        return vec![f.clone()];
    }
    bivariate_factor_square_free_monic(f, x_var, y_var)
}

/// Factor a bivariate polynomial over the integers into irreducible factors with
/// multiplicities.
///
/// The input is treated as a polynomial in variable `x` (index `x_var`) with
/// coefficients in $\mathbb{Z}[y]$ (variable index `y_var`). The current
/// implementation handles the case where the leading coefficient in $x$ is a
/// nonzero integer constant.
pub fn bivariate_factor_z(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<(ZMPoly, usize)> {
    if f.is_zero() || f.total_degree() == Some(0) {
        return Vec::new();
    }

    let content = f.content();
    let mut result = Vec::new();
    if !IntegerDomain.is_one(&content) {
        result.push((
            ZMPoly::from_terms(
                IntegerDomain,
                f.n_vars(),
                vec![(vec![0; f.n_vars()], content)],
            ),
            1,
        ));
    }

    let primitive = f.primitive_part();
    if primitive.total_degree() == Some(0) {
        return result;
    }

    let sqfree = square_free_factorization_bivariate(&primitive, x_var, y_var);
    for (g, m) in sqfree {
        if is_one_mpoly(&g) {
            continue;
        }
        for irr in bivariate_factor_square_free(&g, x_var, y_var) {
            result.push((irr, m));
        }
    }

    result
}

/// Square-free factorization of a bivariate integer polynomial using the
/// heuristic bivariate GCD.
fn square_free_factorization_bivariate(
    f: &ZMPoly,
    x_var: usize,
    y_var: usize,
) -> Vec<(ZMPoly, usize)> {
    let f_deriv = derivative_in_var(f, x_var);
    let mut g = crate::multivariate_gcd::bivariate_gcd(f, &f_deriv)
        .unwrap_or_else(|| one_mpoly(f.n_vars()));
    let mut w = divide_bivariate_by_gcd(f, &g, x_var, y_var);

    let mut result = Vec::new();
    let mut k = 1usize;
    while !is_one_mpoly(&w) && w.total_degree() != Some(0) {
        let h =
            crate::multivariate_gcd::bivariate_gcd(&w, &g).unwrap_or_else(|| one_mpoly(f.n_vars()));
        let z = divide_bivariate_by_gcd(&w, &h, x_var, y_var);
        if !is_one_mpoly(&z) && z.total_degree() != Some(0) {
            result.push((z, k));
        }
        w = h;
        g = divide_bivariate_by_gcd(&g, &w, x_var, y_var);
        k += 1;
    }
    result
}

/// The constant polynomial 1 in `n_vars` variables.
fn one_mpoly(n_vars: usize) -> ZMPoly {
    ZMPoly::from_terms(
        IntegerDomain,
        n_vars,
        vec![(vec![0; n_vars], Integer::from(1))],
    )
}

/// Divide bivariate polynomial `a` by `b` assuming `b` divides `a` exactly in
/// $(\mathbb{Z}[y])[x]$. Returns `a` if the division fails.
fn divide_bivariate_by_gcd(a: &ZMPoly, b: &ZMPoly, x_var: usize, y_var: usize) -> ZMPoly {
    if b.is_zero() || is_one_mpoly(b) {
        return a.clone();
    }
    if a.is_zero() {
        return a.clone();
    }
    let deg_y_a = degree_in_var(a, y_var);
    let deg_y_b = degree_in_var(b, y_var);
    let n_points = deg_y_a.max(deg_y_b) + 2;

    let mut images: Vec<(Integer, ZPoly)> = Vec::new();
    let mut eval_point = Integer::from(0);
    for _ in 0..n_points + 10 {
        if images.len() >= n_points {
            break;
        }
        let a_eval = eval_to_univariate(a, y_var, &eval_point);
        let b_eval = eval_to_univariate(b, y_var, &eval_point);
        if b_eval.is_zero() || a_eval.is_zero() {
            eval_point = Integer::from(eval_point.inner() + 1);
            continue;
        }
        let (q, r) = a_eval.div_rem(&b_eval).unwrap();
        if !r.is_zero() {
            eval_point = Integer::from(eval_point.inner() + 1);
            continue;
        }
        images.push((eval_point.clone(), q));
        eval_point = Integer::from(eval_point.inner() + 1);
    }

    if images.len() < n_points {
        return a.clone();
    }

    interpolate_bivariate_quotient(&images, a.n_vars(), x_var, y_var)
}

/// Interpolate a bivariate quotient from univariate images at distinct y-values.
fn interpolate_bivariate_quotient(
    images: &[(Integer, ZPoly)],
    n_vars: usize,
    x_var: usize,
    y_var: usize,
) -> ZMPoly {
    let mut result = ZMPoly::new(IntegerDomain, n_vars);
    if images.is_empty() {
        return result;
    }
    let max_x_deg = images
        .iter()
        .map(|(_, g)| g.degree().unwrap_or(0))
        .max()
        .unwrap_or(0);
    for x_pow in 0..=max_x_deg {
        let mut y_points: Vec<(Integer, Integer)> = Vec::new();
        for (y_val, g) in images {
            if let Some(c) = g.coeff(x_pow) {
                y_points.push((y_val.clone(), c.clone()));
            }
        }
        if y_points.len() < 2 {
            continue;
        }
        let y_poly = lagrange_interpolate(&y_points);
        for (y_pow, c) in y_poly.coeffs().iter().enumerate() {
            if !IntegerDomain.is_zero(c) {
                let mut exp = vec![0; n_vars];
                exp[x_var] = x_pow;
                exp[y_var] = y_pow;
                result.set_term_external(exp, c.clone());
            }
        }
    }
    result
}

/// Lagrange interpolation of an integer polynomial from points
/// $(y_i, v_i)$. Returns a dense univariate polynomial in $y$.
fn lagrange_interpolate(points: &[(Integer, Integer)]) -> ZPoly {
    let n = points.len();
    let mut result = ZPoly::from_coeffs(IntegerDomain, Vec::new());
    for i in 0..n {
        let (y_i, v_i) = &points[i];
        let mut numerator = ZPoly::from_coeffs(IntegerDomain, vec![Integer::from(1)]);
        let mut denom = BigInt::one();
        for (j, (y_j, _v_j)) in points.iter().enumerate().take(n) {
            if i == j {
                continue;
            }
            let factor = ZPoly::from_coeffs(
                IntegerDomain,
                vec![Integer::from(-y_j.inner().clone()), Integer::from(1)],
            );
            numerator = numerator.mul(&factor);
            denom *= y_i.inner() - y_j.inner();
        }
        let q = v_i.inner() / &denom;
        debug_assert!(
            &q * &denom == *v_i.inner(),
            "lagrange_interpolate: non-exact division; \
             v_i={v_i}, denom={denom}, quotient={q}"
        );
        result = result.add(&numerator.mul_scalar(&Integer::from(q)));
    }
    result
}

fn zpoly_to_qpoly(f: &ZPoly) -> QPoly {
    QPoly::from_coeffs(
        RationalDomain,
        f.coeffs()
            .iter()
            .map(|c| Rational::from_bigints(c.inner().clone(), BigInt::one()))
            .collect(),
    )
}

fn qpoly_to_zpoly(f: &QPoly) -> Option<ZPoly> {
    let mut coeffs = Vec::new();
    for r in f.coeffs() {
        let numer = r.inner().numer();
        let denom = r.inner().denom();
        if !denom.is_one() {
            return None;
        }
        coeffs.push(Integer::from(numer.clone()));
    }
    Some(ZPoly::from_coeffs(IntegerDomain, coeffs))
}

fn monic_qpoly(f: &QPoly) -> QPoly {
    if f.is_zero() {
        return f.clone();
    }
    let lc = f.leading_coeff().unwrap();
    let inv = RationalDomain.inv(lc).unwrap();
    f.mul_scalar(&inv)
}

fn extended_gcd_qpoly(a: &QPoly, b: &QPoly) -> (QPoly, QPoly, QPoly) {
    if b.is_zero() {
        let monic_a = monic_qpoly(a);
        let lc = a.leading_coeff().unwrap();
        let inv = RationalDomain.inv(lc).unwrap();
        let s = QPoly::from_coeffs(RationalDomain, vec![inv]);
        return (monic_a, s, QPoly::new(RationalDomain));
    }
    if a.degree().unwrap_or(0) < b.degree().unwrap_or(0) {
        let (g, s, t) = extended_gcd_qpoly(b, a);
        return (g, t, s);
    }
    let (q, r) = a.div_rem(b).expect("Q is a field");
    let (g, s1, t1) = extended_gcd_qpoly(b, &r);
    let s = t1.clone();
    let t = s1.sub(&q.mul(&t1));
    (g, s, t)
}

fn bezout_coefficients_q(factors: &[QPoly]) -> Vec<QPoly> {
    let n = factors.len();
    if n == 1 {
        return vec![factors[0].one()];
    }
    let mut result = vec![factors[0].zero(); n];
    result[0] = factors[0].one();
    let mut accum = factors[0].clone();
    for i in 1..n {
        let (_g, s, t) = extended_gcd_qpoly(&accum, &factors[i]);
        for res in result.iter_mut().take(i) {
            *res = res.mul(&t);
        }
        result[i] = s;
        accum = accum.mul(&factors[i]);
    }
    result
}

fn eval_to_univariate_fp(poly: &FpMPoly, y_var: usize, value: &FiniteFieldElement) -> FpPoly {
    let evaluated = poly.eval(y_var, value);
    let mut coeffs = Vec::new();
    for (exp, c) in evaluated.terms_ref() {
        let idx = exp.first().copied().unwrap_or(0);
        if idx >= coeffs.len() {
            coeffs.resize(idx + 1, poly.domain().zero());
        }
        coeffs[idx] = c.clone();
    }
    FpPoly::from_coeffs(poly.domain().clone(), coeffs)
}

fn univariate_to_bivariate_fp(g: &FpPoly, n_vars: usize, x_var: usize) -> FpMPoly {
    let mut terms = Vec::new();
    for (i, c) in g.coeffs().iter().enumerate() {
        if !g.domain().is_zero(c) {
            let mut exp = vec![0usize; n_vars];
            exp[x_var] = i;
            terms.push((exp, c.clone()));
        }
    }
    FpMPoly::from_terms(g.domain().clone(), n_vars, terms)
}

fn derivative_in_var_fp(poly: &FpMPoly, var_index: usize) -> FpMPoly {
    let mut result = FpMPoly::new(poly.domain().clone(), poly.n_vars());
    for (exp, coeff) in poly.terms_ref() {
        let power = exp.get(var_index).copied().unwrap_or(0);
        if power == 0 {
            continue;
        }
        let mut new_exp = exp.to_vec();
        new_exp[var_index] = power - 1;
        let scalar = poly.domain().cast_u64(power as u64);
        let new_coeff = poly.domain().mul(coeff, &scalar);
        result.set_term_external(new_exp, new_coeff);
    }
    result
}

fn divide_by_k_factorial_fp(poly: FpPoly, k: usize) -> FpPoly {
    let domain = poly.domain().clone();
    let mut fact = domain.one();
    for i in 1..=k {
        fact = domain.mul(&fact, &domain.cast_u64(i as u64));
    }
    let fact_inv = domain.inv(&fact).expect("k! must be invertible mod p");
    let coeffs = poly
        .coeffs()
        .iter()
        .map(|c| domain.mul(c, &fact_inv))
        .collect();
    FpPoly::from_coeffs(domain, coeffs)
}

fn taylor_coeffs_in_y_fp(
    poly: &FpMPoly,
    y_var: usize,
    alpha: &FiniteFieldElement,
    max_k: usize,
) -> Vec<FpPoly> {
    let mut coeffs = Vec::with_capacity(max_k + 1);
    let mut current = poly.clone();
    for k in 0..=max_k {
        let value = eval_to_univariate_fp(&current, y_var, alpha);
        coeffs.push(divide_by_k_factorial_fp(value, k));
        current = derivative_in_var_fp(&current, y_var);
    }
    coeffs
}

fn monic_fppoly(f: &FpPoly) -> FpPoly {
    if f.is_zero() {
        return f.clone();
    }
    let lc = f.leading_coeff().cloned().unwrap();
    let inv = f.domain().inv(&lc).expect("nonzero leading coefficient");
    f.mul_scalar(&inv)
}

fn extended_gcd_fppoly(a: &FpPoly, b: &FpPoly) -> (FpPoly, FpPoly, FpPoly) {
    if b.is_zero() {
        return (a.clone(), a.one(), a.zero());
    }
    if a.degree().unwrap_or(0) < b.degree().unwrap_or(0) {
        let (g, s, t) = extended_gcd_fppoly(b, a);
        return (g, t, s);
    }
    let (q, r) = a.div_rem(b).expect("field division");
    let (g, s1, t1) = extended_gcd_fppoly(b, &r);
    let s = t1.clone();
    let t = s1.sub(&q.mul(&t1));
    (g, s, t)
}

fn bezout_coefficients_fp(factors: &[FpPoly]) -> Vec<FpPoly> {
    let n = factors.len();
    if n == 1 {
        return vec![factors[0].one()];
    }
    let mut result = vec![factors[0].zero(); n];
    result[0] = factors[0].one();
    let mut accum = factors[0].clone();
    for i in 1..n {
        let (_g, s, t) = extended_gcd_fppoly(&accum, &factors[i]);
        for res in result.iter_mut().take(i) {
            *res = res.mul(&t);
        }
        result[i] = s;
        accum = accum.mul(&factors[i]);
    }
    result
}

fn hensel_lift_bivariate_fp(
    f: &FpMPoly,
    alpha: &FiniteFieldElement,
    univariate_factors: &[FpPoly],
    x_var: usize,
    y_var: usize,
) -> Vec<FpMPoly> {
    let n_vars = f.n_vars();
    let d_y = degree_in_var(f, y_var);

    let c_f = taylor_coeffs_in_y_fp(f, y_var, alpha, d_y);
    let bezout = bezout_coefficients_fp(univariate_factors);

    let mut lifted: Vec<FpMPoly> = univariate_factors
        .iter()
        .map(|g| univariate_to_bivariate_fp(g, n_vars, x_var))
        .collect();

    for k in 1..=d_y {
        let mut product = FpMPoly::from_terms(
            f.domain().clone(),
            n_vars,
            vec![(vec![0; n_vars], f.domain().one())],
        );
        for g in &lifted {
            product = product.mul(g);
        }

        let c_product = taylor_coeffs_in_y_fp(&product, y_var, alpha, d_y);
        let error = c_f[k].sub(&c_product[k]);

        for i in 0..lifted.len() {
            let delta = error.mul(&bezout[i]);
            let (_q, remainder) = delta.div_rem(&univariate_factors[i]).unwrap();
            let correction =
                univariate_times_y_minus_alpha_k_fp(&remainder, k, alpha, n_vars, x_var, y_var);
            lifted[i] = lifted[i].add(&correction);
        }
    }

    lifted
}

fn univariate_times_y_minus_alpha_k_fp(
    g: &FpPoly,
    k: usize,
    alpha: &FiniteFieldElement,
    n_vars: usize,
    x_var: usize,
    y_var: usize,
) -> FpMPoly {
    let domain = g.domain().clone();
    let mut terms = Vec::new();
    for (i, c) in g.coeffs().iter().enumerate() {
        if domain.is_zero(c) {
            continue;
        }
        for j in 0..=k {
            let mut exp = vec![0usize; n_vars];
            exp[x_var] = i;
            exp[y_var] = j;
            let alpha_pow = domain.pow(alpha, (k - j) as u64);
            let binom = domain.cast_u64(binomial(k, j));
            let sign = if (k - j) % 2 == 0 {
                domain.one()
            } else {
                domain.neg(&domain.one())
            };
            let coeff = domain.mul(c, &domain.mul(&binom, &domain.mul(&sign, &alpha_pow)));
            terms.push((exp, coeff));
        }
    }
    FpMPoly::from_terms(domain, n_vars, terms)
}

fn choose_evaluation_point_fp(
    f: &FpMPoly,
    y_var: usize,
) -> Option<(FiniteFieldElement, Vec<FpPoly>)> {
    let domain = f.domain().clone();
    let p = domain.prime().clone();
    let mut best: Option<(FiniteFieldElement, Vec<FpPoly>)> = None;
    for a in 0i64..20 {
        if BigInt::from(a) >= p {
            break;
        }
        let alpha = domain.element(a);
        let image = eval_to_univariate_fp(f, y_var, &alpha);
        if image.degree().unwrap_or(0) < 1 || !image.is_square_free() {
            continue;
        }
        let mut factors = image.factor();
        factors.sort_by(|a, b| b.0.degree().unwrap_or(0).cmp(&a.0.degree().unwrap_or(0)));
        let factors: Vec<FpPoly> = factors.into_iter().map(|(g, _)| monic_fppoly(&g)).collect();
        if factors.len() < 2 {
            continue;
        }
        match &best {
            None => best = Some((alpha.clone(), factors)),
            Some((_, best_factors)) => {
                if factors.len() < best_factors.len() {
                    best = Some((alpha.clone(), factors));
                }
            }
        }
    }
    best
}

fn lc_x_is_constant_fp(f: &FpMPoly, x_var: usize) -> bool {
    let deg_x = degree_in_var(f, x_var);
    if deg_x == 0 {
        return true;
    }
    for exp in f.terms_ref().keys() {
        if exp.get(x_var).copied().unwrap_or(0) == deg_x {
            for (i, &e) in exp.iter().enumerate() {
                if i != x_var && e != 0 {
                    return false;
                }
            }
        }
    }
    true
}

fn bivariate_factor_square_free_monic_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<FpMPoly> {
    if degree_in_var(f, x_var) == 0 || degree_in_var(f, y_var) == 0 {
        return vec![f.clone()];
    }

    let (alpha, univariate_factors) = match choose_evaluation_point_fp(f, y_var) {
        Some(v) => v,
        None => return vec![f.clone()],
    };

    if univariate_factors.len() <= 1 {
        return vec![f.clone()];
    }

    let lifted = hensel_lift_bivariate_fp(f, &alpha, &univariate_factors, x_var, y_var);

    let mut product = FpMPoly::from_terms(
        f.domain().clone(),
        f.n_vars(),
        vec![(vec![0; f.n_vars()], f.domain().one())],
    );
    for g in &lifted {
        product = product.mul(g);
    }
    if product == f.clone() {
        return lifted;
    }

    vec![f.clone()]
}

/// Factor a square-free bivariate polynomial over a prime finite field into
/// irreducible factors. The current implementation requires the leading
/// coefficient in `x` to be a nonzero field constant.
fn bivariate_factor_square_free_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<FpMPoly> {
    if !lc_x_is_constant_fp(f, x_var) {
        return vec![f.clone()];
    }
    bivariate_factor_square_free_monic_fp(f, x_var, y_var)
}

/// Factor a bivariate polynomial over a prime finite field into irreducible
/// factors with multiplicities.
///
/// The current implementation handles square-free polynomials whose leading
/// coefficient in `x` is a field constant. Non-square-free inputs are returned
/// as a single factor (a conservative fallback).
pub fn bivariate_factor_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<(FpMPoly, usize)> {
    if f.is_zero() || f.total_degree() == Some(0) {
        return Vec::new();
    }

    let content = f.content();
    let mut result = Vec::new();
    if !f.domain().is_one(&content) {
        result.push((
            FpMPoly::from_terms(
                f.domain().clone(),
                f.n_vars(),
                vec![(vec![0; f.n_vars()], content)],
            ),
            1,
        ));
    }

    let primitive = f.primitive_part();
    if primitive.total_degree() == Some(0) {
        return result;
    }

    // Conservative square-free check: derivative in x must be non-zero.
    let deriv = derivative_in_var_fp(&primitive, x_var);
    if !deriv.is_zero() {
        for irr in bivariate_factor_square_free_fp(&primitive, x_var, y_var) {
            result.push((irr, 1));
        }
    } else {
        result.push((primitive, 1));
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use ocas_domain::{FiniteField, Integer};

    fn mpoly_from_str(coeffs: &[((usize, usize), i64)]) -> ZMPoly {
        let terms: Vec<(Vec<usize>, Integer)> = coeffs
            .iter()
            .map(|((x, y), c)| (vec![*x, *y], Integer::from(*c)))
            .collect();
        ZMPoly::from_terms(IntegerDomain, 2, terms)
    }

    fn fpoly_from_str(coeffs: &[((usize, usize), i64)], p: i64) -> FpMPoly {
        let domain = FiniteField::new(BigInt::from(p));
        let terms: Vec<(Vec<usize>, FiniteFieldElement)> = coeffs
            .iter()
            .map(|((x, y), c)| (vec![*x, *y], domain.element(*c)))
            .collect();
        FpMPoly::from_terms(domain, 2, terms)
    }

    fn one_mpoly_fp(n_vars: usize, domain: &FiniteField) -> FpMPoly {
        FpMPoly::from_terms(
            domain.clone(),
            n_vars,
            vec![(vec![0; n_vars], domain.one())],
        )
    }

    #[test]
    fn factor_monic_bivariate() {
        // (x^2 + y + 1)(x + y + 2)
        // = x^3 + x^2*y + 2*x^2 + x*y + y^2 + 3*y + x + 2
        let f = mpoly_from_str(&[
            ((3, 0), 1),
            ((2, 1), 1),
            ((2, 0), 2),
            ((1, 1), 1),
            ((1, 0), 1),
            ((0, 2), 1),
            ((0, 1), 3),
            ((0, 0), 2),
        ]);
        let factors = bivariate_factor_z(&f, 0, 1);
        let mut product = one_mpoly(2);
        for (g, m) in &factors {
            for _ in 0..*m {
                product = product.mul(g);
            }
        }
        assert!(
            product == f || product == f.neg(),
            "product did not reconstruct f"
        );
        assert!(factors.len() >= 2, "expected at least two factors");
    }

    #[test]
    #[ignore = "non-monic leading coefficient requires Wang LC handling"]
    fn factor_textbook_bivariate_non_monic() {
        // (x^2 + y + x + 1)(3x + y^2 + 4)
        // Leading coefficient in x is 3, so this case requires the full
        // Wang leading-coefficient algorithm (currently not implemented).
        let f = mpoly_from_str(&[
            ((3, 0), 3),
            ((2, 2), 1),
            ((2, 0), 7),
            ((1, 2), 1),
            ((1, 1), 3),
            ((1, 0), 7),
            ((0, 3), 1),
            ((0, 2), 1),
            ((0, 1), 4),
            ((0, 0), 4),
        ]);
        let factors = bivariate_factor_z(&f, 0, 1);
        let mut product = one_mpoly(2);
        for (g, m) in &factors {
            for _ in 0..*m {
                product = product.mul(g);
            }
        }
        assert!(
            product == f || product == f.neg(),
            "product did not reconstruct f"
        );
        assert!(factors.len() >= 2, "expected at least two factors");
    }

    #[test]
    fn factor_monic_bivariate_over_finite_field() {
        // (x^2 + y + 1)(x + y + 2) over F_5
        // = x^3 + x^2*y + 2*x^2 + x*y + y^2 + 3*y + x + 2  (mod 5)
        let f = fpoly_from_str(
            &[
                ((3, 0), 1),
                ((2, 1), 1),
                ((2, 0), 2),
                ((1, 1), 1),
                ((1, 0), 1),
                ((0, 2), 1),
                ((0, 1), 3),
                ((0, 0), 2),
            ],
            5,
        );
        let factors = bivariate_factor_fp(&f, 0, 1);
        let mut product = one_mpoly_fp(2, f.domain());
        for (g, m) in &factors {
            for _ in 0..*m {
                product = product.mul(g);
            }
        }
        assert_eq!(product, f, "product did not reconstruct f");
        assert!(factors.len() >= 2, "expected at least two factors");
    }
}
