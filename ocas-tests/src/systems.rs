//! Polynomial-system generators (cyclic-n, Katsura-n) over ℚ and ℤ_p.
//!
//! These are the standard test systems used by Singular, msolve and the
//! Symbolica benchmark suite.

use num_bigint::BigInt;
use ocas_domain::{FiniteField, Rational, RationalDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::{Grevlex, Lex, MonomialOrder};

/// Build a single rational term.
pub fn term(exps: Vec<usize>, num: i64, den: i64) -> (Vec<usize>, Rational) {
    (exps, Rational::new(num, den))
}

/// The cyclic-n ideal generators over ℚ under an arbitrary monomial order.
///
///   f_k = Σ_{start=0}^{n-1} x_{start} * x_{start+1} * ... * x_{start+k-1}   (indices mod n)
///   f_n = x_0 * x_1 * ... * x_{n-1} - 1
///
/// The terms themselves are order-independent (each monomial appears
/// exactly once with coefficient 1); only the order type parameter `O`
/// differs between instantiations.
pub fn cyclic_q_with_order<O: MonomialOrder + Default>(
    n: usize,
) -> Vec<SparseMultivariatePolynomial<RationalDomain, O>> {
    let d = RationalDomain;
    let mut gens = Vec::with_capacity(n);

    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut exps = vec![0usize; n];
            for j in 0..k {
                exps[(start + j) % n] = 1;
            }
            terms.push(term(exps, 1, 1));
        }
        gens.push(SparseMultivariatePolynomial::from_terms(d, n, terms));
    }

    let full_exps = vec![1usize; n];
    gens.push(SparseMultivariatePolynomial::from_terms(
        d,
        n,
        vec![term(full_exps, 1, 1), term(vec![0usize; n], -1, 1)],
    ));

    gens
}

/// The cyclic-n ideal generators over ℚ (Lex order; legacy signature).
pub fn cyclic_q(n: usize) -> Vec<SparseMultivariatePolynomial<RationalDomain, Lex>> {
    cyclic_q_with_order::<Lex>(n)
}

/// The cyclic-n ideal generators over ℚ (Grevlex order).
pub fn cyclic_q_grevlex(n: usize) -> Vec<SparseMultivariatePolynomial<RationalDomain, Grevlex>> {
    cyclic_q_with_order::<Grevlex>(n)
}

/// The cyclic-n ideal generators over ℤ_p under an arbitrary monomial order.
pub fn cyclic_fp_with_order<O: MonomialOrder + Default>(
    n: usize,
    p: u32,
) -> Vec<SparseMultivariatePolynomial<FiniteField, O>> {
    let field = FiniteField::new(BigInt::from(p));
    let mut gens = Vec::with_capacity(n);

    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut exps = vec![0usize; n];
            for j in 0..k {
                exps[(start + j) % n] = 1;
            }
            terms.push((exps, field.element(1)));
        }
        gens.push(SparseMultivariatePolynomial::from_terms(
            field.clone(),
            n,
            terms,
        ));
    }

    let full_exps = vec![1usize; n];
    gens.push(SparseMultivariatePolynomial::from_terms(
        field.clone(),
        n,
        vec![
            (full_exps, field.element(1)),
            (vec![0usize; n], field.element(p - 1)), // -1 mod p
        ],
    ));

    gens
}

/// The cyclic-n ideal generators over ℤ_p (Lex order; legacy signature).
pub fn cyclic_fp(n: usize, p: u32) -> Vec<SparseMultivariatePolynomial<FiniteField, Lex>> {
    cyclic_fp_with_order::<Lex>(n, p)
}

/// The cyclic-n ideal generators over ℤ_p (Grevlex order).
pub fn cyclic_fp_grevlex(
    n: usize,
    p: u32,
) -> Vec<SparseMultivariatePolynomial<FiniteField, Grevlex>> {
    cyclic_fp_with_order::<Grevlex>(n, p)
}

/// The Katsura-n ideal generators over ℚ.
///
/// Variables are x_0..x_n (n+1 variables); write u_i = x_{|i|} for
/// −n ≤ i ≤ n. The system is:
///
///   eq_k = Σ_{i=−n..n} u_i·u_{k−i} − u_k,   k = 0..n−1
///   eq_lin = u_0 + 2·Σ_{i=1..n} u_i − 1
///
/// giving n+1 equations. katsura-6 (n = 6, 7 variables) is the standard
/// benchmark size.
pub fn katsura_q(n: usize) -> Vec<SparseMultivariatePolynomial<RationalDomain, Lex>> {
    let n_vars = n + 1;
    let d = RationalDomain;
    let mut gens = Vec::with_capacity(n + 1);

    for k in 0..n {
        let mut terms = Vec::new();
        for i in -(n as i64)..=(n as i64) {
            let j = k as i64 - i;
            if !(-(n as i64)..=(n as i64)).contains(&j) {
                continue;
            }
            let mut exps = vec![0usize; n_vars];
            exps[i.unsigned_abs() as usize] += 1;
            exps[j.unsigned_abs() as usize] += 1;
            terms.push(term(exps, 1, 1));
        }
        let mut exps = vec![0usize; n_vars];
        exps[k] = 1;
        terms.push(term(exps, -1, 1));
        gens.push(SparseMultivariatePolynomial::from_terms(d, n_vars, terms));
    }

    let mut terms = vec![term(vec![0usize; n_vars], -1, 1)];
    let mut e0 = vec![0usize; n_vars];
    e0[0] = 1;
    terms.push(term(e0, 1, 1));
    for i in 1..=n {
        let mut exps = vec![0usize; n_vars];
        exps[i] = 1;
        terms.push(term(exps, 2, 1));
    }
    gens.push(SparseMultivariatePolynomial::from_terms(d, n_vars, terms));

    gens
}

/// The Katsura-n ideal generators over ℤ_p.
pub fn katsura_fp(n: usize, p: u32) -> Vec<SparseMultivariatePolynomial<FiniteField, Lex>> {
    let n_vars = n + 1;
    let field = FiniteField::new(BigInt::from(p));
    let mut gens = Vec::with_capacity(n + 1);

    for k in 0..n {
        let mut terms = Vec::new();
        for i in -(n as i64)..=(n as i64) {
            let j = k as i64 - i;
            if !(-(n as i64)..=(n as i64)).contains(&j) {
                continue;
            }
            let mut exps = vec![0usize; n_vars];
            exps[i.unsigned_abs() as usize] += 1;
            exps[j.unsigned_abs() as usize] += 1;
            terms.push((exps, field.element(1)));
        }
        let mut exps = vec![0usize; n_vars];
        exps[k] = 1;
        terms.push((exps, field.element(p - 1))); // -1 mod p
        gens.push(SparseMultivariatePolynomial::from_terms(
            field.clone(),
            n_vars,
            terms,
        ));
    }

    let mut terms = vec![(vec![0usize; n_vars], field.element(p - 1))];
    let mut e0 = vec![0usize; n_vars];
    e0[0] = 1;
    terms.push((e0, field.element(1)));
    for i in 1..=n {
        let mut exps = vec![0usize; n_vars];
        exps[i] = 1;
        terms.push((exps, field.element(2)));
    }
    gens.push(SparseMultivariatePolynomial::from_terms(
        field.clone(),
        n_vars,
        terms,
    ));

    gens
}
