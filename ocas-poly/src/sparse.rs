//! Sparse multivariate polynomial implementation.
//!
//! A [`SparseMultivariatePolynomial`] stores only non-zero terms as a map from
//! exponent vectors to coefficients. The exponent vector `vec![e1, e2, ...]`
//! represents the monomial `x1^e1 * x2^e2 * ...`. Monomial ordering is
//! controlled by the [`MonomialOrder`] type parameter.

use std::collections::HashMap;
use std::marker::PhantomData;

use ocas_domain::Domain;
use smallvec::SmallVec;

/// A monomial ordering determines how terms are sorted and compared.
///
/// Orderings are implemented as zero-sized types with an associated method
/// that compares two exponent vectors.
///
/// # Example
///
/// ```
/// use ocas_poly::sparse::{Grevlex, Lex, MonomialOrder};
///
/// let a = [2, 1];
/// let b = [1, 1];
/// assert_eq!(Lex::cmp(&a, &b), std::cmp::Ordering::Greater);
/// assert_eq!(Grevlex::cmp(&a, &b), std::cmp::Ordering::Less);
/// ```
pub trait MonomialOrder: Clone + Copy + PartialEq + Eq + std::fmt::Debug {
    /// Compare two exponent vectors.
    ///
    /// Returns `std::cmp::Ordering::Less` if `lhs` should appear before `rhs`
    /// in the ordering.
    fn cmp(lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering;
}

/// Lexicographic ordering: compare exponents left-to-right.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Lex;

impl MonomialOrder for Lex {
    fn cmp(lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        lhs.cmp(rhs)
    }
}

/// Graded reverse lexicographic ordering: first by total degree descending,
/// then reverse lexicographic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Grevlex;

impl MonomialOrder for Grevlex {
    fn cmp(lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        let deg_lhs: usize = lhs.iter().sum();
        let deg_rhs: usize = rhs.iter().sum();
        deg_rhs
            .cmp(&deg_lhs)
            .then_with(|| rhs.iter().rev().cmp(lhs.iter().rev()))
    }
}

/// A sparse multivariate polynomial with coefficients in a domain `D` and
/// monomial ordering `O`.
///
/// # Example
///
/// ```
/// use ocas_domain::{IntegerDomain, Integer};
/// use ocas_poly::sparse::Grevlex;
/// use ocas_poly::SparseMultivariatePolynomial;
///
/// let domain = IntegerDomain;
/// let p = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
///     domain,
///     2,
///     vec![(vec![1, 0], Integer::from(2)), (vec![0, 1], Integer::from(3))],
/// );
/// let q = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
///     domain,
///     2,
///     vec![(vec![1, 0], Integer::from(1)), (vec![0, 0], Integer::from(1))],
/// );
/// let r = p.mul(&q);
/// assert_eq!(r.coeff(&[1, 0]), Integer::from(2));
/// assert_eq!(r.coeff(&[0, 1]), Integer::from(3));
/// assert_eq!(r.coeff(&[2, 0]), Integer::from(2));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseMultivariatePolynomial<D: Domain, O: MonomialOrder = Grevlex> {
    /// Non-zero terms indexed by exponent vector.
    terms: HashMap<SmallVec<[usize; 4]>, D::Element>,
    /// The coefficient domain.
    domain: D,
    /// Number of variables. Exponent vectors are padded/trimmed to this length.
    n_vars: usize,
    _marker: PhantomData<O>,
}

impl<D: Domain, O: MonomialOrder> SparseMultivariatePolynomial<D, O> {
    /// Create the zero polynomial in `n_vars` variables over `domain`.
    pub fn new(domain: D, n_vars: usize) -> Self {
        Self {
            terms: HashMap::new(),
            domain,
            n_vars,
            _marker: PhantomData,
        }
    }

    /// Create a polynomial from a list of (exponent vector, coefficient) pairs.
    ///
    /// Zero coefficients and empty terms are dropped automatically.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::sparse::Grevlex;
    /// use ocas_poly::SparseMultivariatePolynomial;
    ///
    /// let domain = IntegerDomain;
    /// let p = SparseMultivariatePolynomial::<IntegerDomain, Grevlex>::from_terms(
    ///     domain,
    ///     2,
    ///     vec![(vec![1, 0], Integer::from(2)), (vec![0, 1], Integer::from(3))],
    /// );
    /// assert_eq!(p.n_terms(), 2);
    /// assert_eq!(p.coeff(&[1, 0]), Integer::from(2));
    /// ```
    pub fn from_terms(domain: D, n_vars: usize, terms: Vec<(Vec<usize>, D::Element)>) -> Self {
        let mut poly = Self::new(domain, n_vars);
        for (exp, coeff) in terms {
            poly.set_term(exp, coeff);
        }
        poly
    }

    /// Return a reference to the coefficient domain.
    pub fn domain(&self) -> &D {
        &self.domain
    }

    /// Return the number of variables.
    pub fn n_vars(&self) -> usize {
        self.n_vars
    }

    /// Return the number of non-zero terms.
    pub fn n_terms(&self) -> usize {
        self.terms.len()
    }

    /// Return whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Return the total degree, or `None` for the zero polynomial.
    pub fn total_degree(&self) -> Option<usize> {
        self.terms.keys().map(|e| e.iter().sum::<usize>()).max()
    }

    /// Return the coefficient of the given monomial, or zero if absent.
    pub fn coeff(&self, exp: &[usize]) -> D::Element {
        let key = Self::normalize_exp(exp, self.n_vars);
        self.terms
            .get(&key)
            .cloned()
            .unwrap_or_else(|| self.domain.zero())
    }

    /// Set the coefficient of a monomial. Zero coefficients remove the term.
    fn set_term(&mut self, exp: Vec<usize>, coeff: D::Element) {
        let key = Self::normalize_exp(&exp, self.n_vars);
        if self.domain.is_zero(&coeff) {
            self.terms.remove(&key);
        } else {
            self.terms.insert(key, coeff);
        }
    }

    fn normalize_exp(exp: &[usize], n_vars: usize) -> SmallVec<[usize; 4]> {
        let mut v = SmallVec::with_capacity(n_vars);
        for i in 0..n_vars {
            v.push(*exp.get(i).unwrap_or(&0));
        }
        v
    }

    /// Return the zero polynomial with the same shape.
    pub fn zero(&self) -> Self {
        Self::new(self.domain.clone(), self.n_vars)
    }

    /// Return the constant polynomial `1` over the same shape.
    pub fn one(&self) -> Self {
        let mut poly = Self::new(self.domain.clone(), self.n_vars);
        let mut exp = SmallVec::with_capacity(self.n_vars);
        exp.resize(self.n_vars, 0);
        poly.terms.insert(exp, self.domain.one());
        poly
    }

    /// Return the negation of this polynomial.
    pub fn neg(&self) -> Self {
        let mut poly = self.zero();
        for (exp, coeff) in &self.terms {
            poly.terms.insert(exp.clone(), self.domain.neg(coeff));
        }
        poly
    }

    /// Add another polynomial.
    ///
    /// Panics if the polynomials have different numbers of variables.
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(
            self.n_vars, other.n_vars,
            "polynomials must have the same number of variables"
        );
        let mut poly = self.clone();
        for (exp, coeff) in &other.terms {
            let existing = poly
                .terms
                .get(exp)
                .cloned()
                .unwrap_or_else(|| poly.domain.zero());
            let sum = poly.domain.add(&existing, coeff);
            if poly.domain.is_zero(&sum) {
                poly.terms.remove(exp);
            } else {
                poly.terms.insert(exp.clone(), sum);
            }
        }
        poly
    }

    /// Subtract another polynomial.
    ///
    /// Panics if the polynomials have different numbers of variables.
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Multiply by a scalar coefficient.
    pub fn mul_scalar(&self, scalar: &D::Element) -> Self {
        if self.domain.is_zero(scalar) {
            return self.zero();
        }
        let mut poly = self.zero();
        for (exp, coeff) in &self.terms {
            poly.terms
                .insert(exp.clone(), self.domain.mul(coeff, scalar));
        }
        poly
    }

    /// Multiply two polynomials.
    ///
    /// Panics if the polynomials have different numbers of variables.
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(
            self.n_vars, other.n_vars,
            "polynomials must have the same number of variables"
        );
        if self.is_zero() || other.is_zero() {
            return self.zero();
        }
        let mut poly = self.zero();
        for (e1, c1) in &self.terms {
            for (e2, c2) in &other.terms {
                let mut exp = SmallVec::with_capacity(self.n_vars);
                for i in 0..self.n_vars {
                    exp.push(e1[i] + e2[i]);
                }
                let prod = self.domain.mul(c1, c2);
                let existing = poly
                    .terms
                    .get(&exp)
                    .cloned()
                    .unwrap_or_else(|| poly.domain.zero());
                let sum = poly.domain.add(&existing, &prod);
                if poly.domain.is_zero(&sum) {
                    poly.terms.remove(&exp);
                } else {
                    poly.terms.insert(exp, sum);
                }
            }
        }
        poly
    }

    /// Return the terms sorted according to the monomial ordering.
    pub fn sorted_terms(&self) -> Vec<(&SmallVec<[usize; 4]>, &D::Element)> {
        let mut terms: Vec<_> = self.terms.iter().collect();
        terms.sort_by(|(a, _), (b, _)| O::cmp(a, b));
        terms
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain, Rational, RationalDomain};

    #[test]
    fn sparse_create_and_coeff() {
        let domain = IntegerDomain;
        let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Integer::from(2)),
                (vec![0, 1], Integer::from(3)),
            ],
        );
        assert_eq!(p.coeff(&[1, 0]), Integer::from(2));
        assert_eq!(p.coeff(&[0, 1]), Integer::from(3));
        assert_eq!(p.coeff(&[0, 0]), Integer::from(0));
    }

    #[test]
    fn sparse_total_degree() {
        let domain = IntegerDomain;
        let p = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(
            domain,
            2,
            vec![
                (vec![2, 1], Integer::from(1)),
                (vec![1, 0], Integer::from(1)),
            ],
        );
        assert_eq!(p.total_degree(), Some(3));
    }

    #[test]
    fn sparse_add_and_sub() {
        let domain = IntegerDomain;
        let a = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Integer::from(1)),
                (vec![0, 1], Integer::from(2)),
            ],
        );
        let b = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Integer::from(3)),
                (vec![0, 0], Integer::from(4)),
            ],
        );
        let sum = a.add(&b);
        assert_eq!(sum.coeff(&[1, 0]), Integer::from(4));
        assert_eq!(sum.coeff(&[0, 1]), Integer::from(2));
        assert_eq!(sum.coeff(&[0, 0]), Integer::from(4));

        let diff = b.sub(&a);
        assert_eq!(diff.coeff(&[1, 0]), Integer::from(2));
        assert_eq!(diff.coeff(&[0, 1]), Integer::from(-2));
        assert_eq!(diff.coeff(&[0, 0]), Integer::from(4));
    }

    #[test]
    fn sparse_multiplication() {
        let domain = RationalDomain;
        // (x + 2y) * (3x + y) = 3x^2 + 7xy + 2y^2
        let a = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Rational::new(1, 1)),
                (vec![0, 1], Rational::new(2, 1)),
            ],
        );
        let b = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Rational::new(3, 1)),
                (vec![0, 1], Rational::new(1, 1)),
            ],
        );
        let prod = a.mul(&b);
        assert_eq!(prod.coeff(&[2, 0]), Rational::new(3, 1));
        assert_eq!(prod.coeff(&[1, 1]), Rational::new(7, 1));
        assert_eq!(prod.coeff(&[0, 2]), Rational::new(2, 1));
    }

    #[test]
    fn sparse_sorted_terms_grevlex() {
        let domain = IntegerDomain;
        let p = SparseMultivariatePolynomial::<_, Grevlex>::from_terms(
            domain,
            2,
            vec![
                (vec![1, 0], Integer::from(1)),
                (vec![2, 0], Integer::from(1)),
                (vec![0, 1], Integer::from(1)),
            ],
        );
        let sorted = p.sorted_terms();
        let exps: Vec<_> = sorted.into_iter().map(|(e, _)| e.to_vec()).collect();
        // Grevlex order for these terms: x^2 (degree 2), x (degree 1), y (degree 1).
        // Among degree-1 terms, reverse lex compares the last non-zero exponent:
        // y = [0,1] comes before x = [1,0] because 1 > 0 in the last position.
        assert_eq!(exps, vec![vec![2, 0], vec![0, 1], vec![1, 0]]);
    }
}
