//! Sparse multivariate polynomial implementation.
//!
//! A [`SparseMultivariatePolynomial`] stores only non-zero terms as a map from
//! exponent vectors to coefficients. The exponent vector `vec![e1, e2, ...]`
//! represents the monomial `x1^e1 * x2^e2 * ...`. Monomial ordering is
//! controlled by the [`MonomialOrder`] type parameter.

use ocas_core::FastHashMap as HashMap;
use ocas_domain::{Domain, EuclideanDomain, FiniteField, IntegerDomain};
use smallvec::SmallVec;

use crate::factor::multivariate::{bivariate_factor_fp, bivariate_factor_z};

/// A monomial ordering determines how terms are sorted and compared.
///
/// Simple orderings (Lex, Grevlex, Grlex) are zero-sized types with no
/// runtime data. Parameterized orderings (WeightOrder, BlockOrder) carry
/// configuration at runtime.
///
/// # Example
///
/// ```
/// use ocas_poly::sparse::{Grevlex, Lex, MonomialOrder};
///
/// let a = [2, 1];
/// let b = [1, 1];
/// assert_eq!(Lex.cmp(&a, &b), std::cmp::Ordering::Greater);
/// assert_eq!(Grevlex.cmp(&a, &b), std::cmp::Ordering::Less);
/// ```
pub trait MonomialOrder: Clone + PartialEq + Eq + std::fmt::Debug + Default {
    /// Compare two exponent vectors.
    ///
    /// Returns `std::cmp::Ordering::Less` if `lhs` should appear before `rhs`
    /// in the ordering.
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering;
}

/// Lexicographic ordering: compare exponents left-to-right.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Lex;

impl MonomialOrder for Lex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        lhs.cmp(rhs)
    }
}

/// Graded reverse lexicographic ordering: first by total degree descending,
/// then reverse lexicographic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grevlex;

impl MonomialOrder for Grevlex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        let deg_lhs: usize = lhs.iter().sum();
        let deg_rhs: usize = rhs.iter().sum();
        deg_rhs
            .cmp(&deg_lhs)
            .then_with(|| rhs.iter().rev().cmp(lhs.iter().rev()))
    }
}

/// Graded lexicographic ordering: first by total degree descending,
/// then lexicographic.
///
/// Grlex is sometimes preferred over grevlex in Gröbner basis computations
/// because it can lead to smaller intermediate matrices in the F4 algorithm.
///
/// # Example
///
/// ```
/// use ocas_poly::sparse::{Grlex, MonomialOrder};
///
/// let a = [2, 0]; // x^2, degree 2
/// let b = [1, 1]; // x*y, degree 2
/// let c = [0, 3]; // y^3, degree 3
/// // c has highest degree, so it comes first
/// assert_eq!(Grlex.cmp(&c, &a), std::cmp::Ordering::Less);
/// // a and b have same degree; a > b lexicographically
/// assert_eq!(Grlex.cmp(&a, &b), std::cmp::Ordering::Greater);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Grlex;

impl MonomialOrder for Grlex {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        let deg_lhs: usize = lhs.iter().sum();
        let deg_rhs: usize = rhs.iter().sum();
        deg_rhs.cmp(&deg_lhs).then_with(|| lhs.cmp(rhs))
    }
}

/// Weighted ordering: compare by $\sum_i w_i \cdot e_i$ descending.
///
/// The weight vector is stored at construction time, enabling arbitrary
/// elimination orderings that cannot be expressed as zero-sized types.
///
/// # Example
///
/// ```
/// use ocas_poly::sparse::{MonomialOrder, WeightOrder};
/// use smallvec::smallvec;
///
/// let ord = WeightOrder::new(smallvec![2, 1]);
/// // [1,0] → weight 2, [0,1] → weight 1 → [1,0] is "larger"
/// assert_eq!(ord.cmp(&[1, 0], &[0, 1]), std::cmp::Ordering::Less);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightOrder {
    weights: SmallVec<[i64; 4]>,
}

impl WeightOrder {
    /// Create a new weighted ordering with the given weight vector.
    pub fn new(weights: SmallVec<[i64; 4]>) -> Self {
        Self { weights }
    }

    /// Create a weighted ordering from a slice.
    pub fn from_slice(weights: &[i64]) -> Self {
        Self {
            weights: SmallVec::from_slice(weights),
        }
    }
}

impl Default for WeightOrder {
    /// Default: all-ones weights (total degree ordering).
    fn default() -> Self {
        Self {
            weights: smallvec::smallvec![1; 4],
        }
    }
}

impl MonomialOrder for WeightOrder {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        let w_lhs: i64 = lhs
            .iter()
            .zip(self.weights.iter())
            .map(|(&e, &w)| w * e as i64)
            .sum();
        let w_rhs: i64 = rhs
            .iter()
            .zip(self.weights.iter())
            .map(|(&e, &w)| w * e as i64)
            .sum();
        // Higher weight first (descending).
        w_rhs.cmp(&w_lhs)
    }
}

/// A sub-ordering used inside [`BlockOrder`] for each variable block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubOrder {
    /// Lexicographic within the block.
    Lex,
    /// Graded reverse lexicographic within the block.
    Grevlex,
    /// Graded lexicographic within the block.
    Grlex,
}

impl SubOrder {
    fn cmp_block(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        match self {
            SubOrder::Lex => lhs.cmp(rhs),
            SubOrder::Grevlex => {
                let deg_l: usize = lhs.iter().sum();
                let deg_r: usize = rhs.iter().sum();
                deg_r
                    .cmp(&deg_l)
                    .then_with(|| rhs.iter().rev().cmp(lhs.iter().rev()))
            }
            SubOrder::Grlex => {
                let deg_l: usize = lhs.iter().sum();
                let deg_r: usize = rhs.iter().sum();
                deg_r.cmp(&deg_l).then_with(|| lhs.cmp(rhs))
            }
        }
    }
}

/// Block ordering: partition variables into contiguous blocks, each compared
/// under its own sub-ordering.
///
/// Blocks are defined by `boundaries`: a sorted list of split points
/// (exclusive upper bounds, *not* including `n_vars`). For example,
/// `boundaries = [2]` with `orders = [Lex, Grevlex]` on a 4-variable
/// polynomial means: compare variables 0–1 under Lex first; if equal,
/// compare variables 2–3 under Grevlex.
///
/// # Example
///
/// ```
/// use ocas_poly::sparse::{BlockOrder, MonomialOrder, SubOrder};
/// use smallvec::smallvec;
///
/// let ord = BlockOrder::new(smallvec![2], smallvec![SubOrder::Lex, SubOrder::Grevlex]);
/// // First compare variables 0–1 lex, then variables 2–3 grevlex.
/// let a = [1, 0, 0, 0]; // x₀
/// let b = [0, 1, 0, 0]; // x₁
/// // Lex: [1,0] > [0,1], so a is "greater" (comes first in ordering)
/// assert_eq!(ord.cmp(&a, &b), std::cmp::Ordering::Greater);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockOrder {
    /// Sorted split points (exclusive upper bounds, excluding n_vars).
    boundaries: SmallVec<[usize; 4]>,
    /// One sub-ordering per block (len = boundaries.len() + 1).
    orders: SmallVec<[SubOrder; 4]>,
}

impl BlockOrder {
    /// Create a new block ordering.
    ///
    /// `boundaries` must be sorted in ascending order and not include `n_vars`.
    /// `orders.len()` must equal `boundaries.len() + 1`.
    pub fn new(boundaries: SmallVec<[usize; 4]>, orders: SmallVec<[SubOrder; 4]>) -> Self {
        debug_assert_eq!(orders.len(), boundaries.len() + 1);
        Self { boundaries, orders }
    }
}

impl Default for BlockOrder {
    /// Default: single block with Grevlex.
    fn default() -> Self {
        Self {
            boundaries: SmallVec::new(),
            orders: smallvec::smallvec![SubOrder::Grevlex],
        }
    }
}

impl MonomialOrder for BlockOrder {
    fn cmp(&self, lhs: &[usize], rhs: &[usize]) -> std::cmp::Ordering {
        let mut start = 0;
        for (i, &end) in self.boundaries.iter().enumerate() {
            match self.orders[i].cmp_block(&lhs[start..end], &rhs[start..end]) {
                std::cmp::Ordering::Equal => {}
                ord => return ord,
            }
            start = end;
        }
        // Last block: from `start` to end of slice.
        self.orders[self.boundaries.len()].cmp_block(&lhs[start..], &rhs[start..])
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
    /// The monomial ordering used for leading-term and sorting operations.
    pub order: O,
}

impl<D: Domain, O: MonomialOrder> SparseMultivariatePolynomial<D, O> {
    /// Create the zero polynomial in `n_vars` variables over `domain`
    /// with the default monomial ordering.
    pub fn new(domain: D, n_vars: usize) -> Self {
        Self {
            terms: HashMap::default(),
            domain,
            n_vars,
            order: O::default(),
        }
    }

    /// Create the zero polynomial with an explicit monomial ordering.
    pub fn new_with_order(domain: D, n_vars: usize, order: O) -> Self {
        Self {
            terms: HashMap::default(),
            domain,
            n_vars,
            order,
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

    /// Return a reference to the internal term map (exponent → coefficient).
    pub fn terms_ref(&self) -> &HashMap<SmallVec<[usize; 4]>, D::Element> {
        &self.terms
    }

    /// Set the coefficient of a monomial (public version of `set_term`).
    /// Zero coefficients remove the term.
    pub fn set_term_external(&mut self, exp: Vec<usize>, coeff: D::Element) {
        self.set_term(exp, coeff);
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
        terms.sort_by(|(a, _), (b, _)| self.order.cmp(a, b));
        terms
    }

    // ------------------------------------------------------------------
    //  Gröbner-basis support
    // ------------------------------------------------------------------

    /// Return the leading term `(exponent_vector, coefficient)` or `None`
    /// for the zero polynomial.
    ///
    /// This scans the HashMap in O(n) without allocating — faster than
    /// `sorted_terms()` for repeated calls during reduction.
    pub fn leading_term(&self) -> Option<(&SmallVec<[usize; 4]>, &D::Element)> {
        self.terms
            .iter()
            .max_by(|(a, _), (b, _)| self.order.cmp(a, b))
    }

    /// Return the leading monomial (exponent vector) or `None`.
    pub fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.terms.keys().max_by(|a, b| self.order.cmp(a, b))
    }

    /// Return the leading coefficient or `None`.
    pub fn leading_coeff(&self) -> Option<&D::Element> {
        let lm = self.leading_monomial()?;
        self.terms.get(lm)
    }

    /// Multiply every term's exponent vector by `exp` element-wise.
    ///
    /// Panics if `exp.len() != self.n_vars`.
    pub fn mul_monomial(&self, exp: &[usize]) -> Self {
        assert_eq!(
            exp.len(),
            self.n_vars,
            "exponent vector must have length {}",
            self.n_vars
        );
        let mut poly = self.zero();
        for (e, c) in &self.terms {
            let mut new_exp = SmallVec::with_capacity(self.n_vars);
            for i in 0..self.n_vars {
                new_exp.push(e[i] + exp[i]);
            }
            poly.terms.insert(new_exp, c.clone());
        }
        poly
    }

    /// Reduce `self` by the given basis (a list of polynomials).
    ///
    /// Implements multivariate polynomial division: repeatedly look for a
    /// basis element whose leading monomial divides the current leading
    /// monomial, subtract the appropriate multiple, or else move the leading
    /// term into the remainder.  Requires that `div` on the domain succeeds
    /// (i.e. the domain is effectively a field).
    pub fn reduce(&self, basis: &[Self]) -> Self {
        let mut remainder = self.clone();
        let mut result = self.zero();

        // Cache each basis element's leading term.
        let basis_lts: Vec<_> = basis
            .iter()
            .filter_map(|g| g.leading_term().map(|(e, c)| (g, e.clone(), c.clone())))
            .collect();

        let max_iter = 10000;

        for _ in 0..max_iter {
            if remainder.is_zero() {
                break;
            }
            let (rm, rc) = match remainder.leading_term() {
                Some((e, c)) => (e.clone(), c.clone()),
                None => break,
            };

            let mut reduced = false;
            for (g, lm, lc) in &basis_lts {
                if monomial_divides(&rm, lm) {
                    let qm: SmallVec<[usize; 4]> =
                        rm.iter().zip(lm.iter()).map(|(a, b)| a - b).collect();
                    let qc = match self.domain.div(&rc, lc) {
                        Some(q) => q,
                        None => break,
                    };
                    let sub = g.mul_monomial(&qm).mul_scalar(&qc);
                    remainder = remainder.sub(&sub);
                    reduced = true;
                    break;
                }
            }

            if !reduced {
                let key = rm;
                let val = rc;
                result.terms.insert(key.clone(), val);
                remainder.terms.remove(&key);
            }
        }

        result
    }

    /// Compute the S-polynomial of `self` and `other`:
    ///
    /// S(f, g) = f·lc(g)·x^(lcm-lm(f)) - g·lc(f)·x^(lcm-lm(g))
    pub fn spoly(&self, other: &Self) -> Self {
        let (lm_f, lc_f) = match self.leading_term() {
            Some(t) => (t.0.clone(), t.1.clone()),
            None => return self.zero(),
        };
        let (lm_g, lc_g) = match other.leading_term() {
            Some(t) => (t.0.clone(), t.1.clone()),
            None => return self.zero(),
        };

        let lcm = monomial_lcm(&lm_f, &lm_g);

        let m_f: SmallVec<[usize; 4]> = lcm.iter().zip(lm_f.iter()).map(|(a, b)| a - b).collect();
        let m_g: SmallVec<[usize; 4]> = lcm.iter().zip(lm_g.iter()).map(|(a, b)| a - b).collect();

        let term1 = self.mul_monomial(&m_f).mul_scalar(&lc_g);
        let term2 = other.mul_monomial(&m_g).mul_scalar(&lc_f);

        term1.sub(&term2)
    }

    // ------------------------------------------------------------------
    //  Multivariate GCD support
    // ------------------------------------------------------------------

    /// Compute the content: the GCD of all coefficients.
    ///
    /// For the zero polynomial the content is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::SparseMultivariatePolynomial;
    /// use ocas_poly::Lex;
    ///
    /// let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    ///     IntegerDomain, 1,
    ///     vec![(vec![2], Integer::from(6)), (vec![1], Integer::from(9)), (vec![0], Integer::from(3))],
    /// );
    /// assert_eq!(p.content(), Integer::from(3));
    /// ```
    pub fn content(&self) -> D::Element
    where
        D: EuclideanDomain,
    {
        if self.is_zero() {
            return self.domain.zero();
        }
        let mut g = self.domain.zero();
        for c in self.terms.values() {
            g = self.domain.gcd(&g, c);
        }
        g
    }

    /// Return the primitive part: `self / content`.
    ///
    /// The result has content 1 (or is zero).
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::SparseMultivariatePolynomial;
    /// use ocas_poly::Lex;
    ///
    /// let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    ///     IntegerDomain, 1,
    ///     vec![(vec![2], Integer::from(6)), (vec![0], Integer::from(3))],
    /// );
    /// let pp = p.primitive_part();
    /// // After dividing by content=3: 2*x^2 + 1
    /// assert_eq!(pp.coeff(&[2]), Integer::from(2));
    /// assert_eq!(pp.coeff(&[0]), Integer::from(1));
    /// ```
    pub fn primitive_part(&self) -> Self
    where
        D: EuclideanDomain,
    {
        if self.is_zero() {
            return self.clone();
        }
        let content = self.content();
        if self.domain.is_one(&content) {
            return self.clone();
        }
        let mut result = self.zero();
        for (exp, c) in &self.terms {
            let q = self.domain.div(c, &content).unwrap_or_else(|| c.clone());
            result.terms.insert(exp.clone(), q);
        }
        result
    }

    /// Divide this polynomial by another, assuming the division is exact
    /// (no remainder).
    ///
    /// Each term of `self` is divided by the corresponding factor from
    /// `divisor`. This is used in rational-function canonicalization where
    /// the GCD is known to divide both numerator and denominator.
    ///
    /// # Panics
    ///
    /// Panics if the division is not exact.
    pub fn div_exact(&self, divisor: &Self) -> Self {
        if divisor.n_terms() <= 1 {
            // Check if divisor is constant 1 (or zero).
            let const_val = divisor.coeff(&vec![0; divisor.n_vars]);
            if self.domain.is_one(&const_val) {
                return self.clone();
            }
        }
        let (quot, rem) = self.div_rem_sparse(divisor);
        debug_assert!(rem.is_zero(), "div_exact: division had non-zero remainder");
        quot
    }

    /// Sparse polynomial long division returning (quotient, remainder).
    fn div_rem_sparse(&self, divisor: &Self) -> (Self, Self) {
        if divisor.is_zero() {
            panic!("division by zero polynomial");
        }
        let (_, div_lm) = match divisor.leading_term() {
            Some(t) => (t.0.clone(), t.1.clone()),
            None => return (self.zero(), self.clone()),
        };
        let div_lc = div_lm;
        let div_exp = divisor.leading_monomial().unwrap().clone();

        let mut remainder = self.clone();
        let mut quotient = self.zero();

        while !remainder.is_zero() {
            let (rem_exp, rem_lc) = match remainder.leading_term() {
                Some(t) => (t.0.clone(), t.1.clone()),
                None => break,
            };
            // Check if leading monomial of divisor divides leading monomial of remainder.
            if !monomial_divides(&rem_exp, &div_exp) {
                break;
            }
            let q_coeff = match self.domain.div(&rem_lc, &div_lc) {
                Some(q) => q,
                None => break,
            };
            let q_exp: SmallVec<[usize; 4]> = rem_exp
                .iter()
                .zip(div_exp.iter())
                .map(|(a, b)| a - b)
                .collect();
            // quotient += q_coeff * x^q_exp
            let existing = quotient
                .terms
                .get(&q_exp)
                .cloned()
                .unwrap_or_else(|| self.domain.zero());
            let sum = self.domain.add(&existing, &q_coeff);
            if self.domain.is_zero(&sum) {
                quotient.terms.remove(&q_exp);
            } else {
                quotient.terms.insert(q_exp, sum);
            }
            // remainder -= q_coeff * x^q_exp * divisor
            let scaled = divisor.mul_monomial(
                &remainder
                    .leading_monomial()
                    .unwrap()
                    .iter()
                    .zip(div_exp.iter())
                    .map(|(a, b)| a - b)
                    .collect::<SmallVec<[usize; 4]>>(),
            );
            let scaled = scaled.mul_scalar(&q_coeff);
            remainder = remainder.sub(&scaled);
        }
        (quotient, remainder)
    }

    /// Return the degree of this polynomial in the given variable.
    ///
    /// Returns 0 for the zero polynomial (by convention) or if the variable
    /// does not appear.
    pub fn degree_in(&self, var_index: usize) -> usize {
        self.terms
            .keys()
            .map(|e| e.get(var_index).copied().unwrap_or(0))
            .max()
            .unwrap_or(0)
    }

    // ------------------------------------------------------------------
    //  Multivariate factorization support
    // ------------------------------------------------------------------

    /// Return the coefficient polynomial of `x_var^pow`: the sum of all terms
    /// whose exponent in `var_index` equals `pow`, with that exponent zeroed
    /// out. The result has the same number of variables and does not depend
    /// on `var_index`.
    pub fn coeff_of_var_pow(&self, var_index: usize, pow: usize) -> Self {
        let mut result = Self::new(self.domain.clone(), self.n_vars);
        for (exp, coeff) in &self.terms {
            if exp.get(var_index).copied().unwrap_or(0) == pow {
                let mut new_exp = exp.clone();
                if var_index < new_exp.len() {
                    new_exp[var_index] = 0;
                }
                result.terms.insert(new_exp, coeff.clone());
            }
        }
        result
    }

    /// Return the leading coefficient when this polynomial is viewed as a
    /// univariate polynomial in `var_index`. The result is a polynomial in
    /// the remaining variables (same `n_vars`, exponent of `var_index` is 0).
    pub fn leading_coeff_in(&self, var_index: usize) -> Self {
        self.coeff_of_var_pow(var_index, self.degree_in(var_index))
    }

    /// Compute the formal partial derivative with respect to `var_index`.
    pub fn derivative(&self, var_index: usize) -> Self {
        let mut result = Self::new(self.domain.clone(), self.n_vars);
        for (exp, coeff) in &self.terms {
            let power = exp.get(var_index).copied().unwrap_or(0);
            if power == 0 {
                continue;
            }
            let mut new_exp = exp.clone();
            new_exp[var_index] = power - 1;
            let scalar = self.domain.cast_u64(power as u64);
            let new_coeff = self.domain.mul(coeff, &scalar);
            let existing = result
                .terms
                .get(&new_exp)
                .cloned()
                .unwrap_or_else(|| self.domain.zero());
            let sum = self.domain.add(&existing, &new_coeff);
            if self.domain.is_zero(&sum) {
                result.terms.remove(&new_exp);
            } else {
                result.terms.insert(new_exp, sum);
            }
        }
        result
    }

    /// Compute the Taylor coefficients in variable `var_index` around `a`:
    /// `f = Σ_j t_j · (x_var - a)^j` where each `t_j` does not depend on
    /// `var_index` (its exponent is zeroed).
    ///
    /// Returns `t_0, t_1, ..., t_d` with `d = degree_in(var_index)`.
    pub fn taylor_coefficients(&self, var_index: usize, a: &D::Element) -> Vec<Self> {
        let d = self.degree_in(var_index);
        let mut coeffs = vec![Self::new(self.domain.clone(), self.n_vars); d + 1];
        for (exp, coeff) in &self.terms {
            let e = exp.get(var_index).copied().unwrap_or(0);
            let mut base_exp = exp.clone();
            if var_index < base_exp.len() {
                base_exp[var_index] = 0;
            }
            // x_v^e = Σ_j binom(e, j) · a^(e-j) · (x_v - a)^j
            for (j, t_j) in coeffs.iter_mut().enumerate().take(e + 1) {
                let binom = self.domain.cast_u64(binomial(e, j));
                let a_pow = self.domain.pow(a, (e - j) as u64);
                let contrib = self.domain.mul(coeff, &self.domain.mul(&binom, &a_pow));
                let existing = t_j
                    .terms
                    .get(&base_exp)
                    .cloned()
                    .unwrap_or_else(|| self.domain.zero());
                let sum = self.domain.add(&existing, &contrib);
                if self.domain.is_zero(&sum) {
                    t_j.terms.remove(&base_exp);
                } else {
                    t_j.terms.insert(base_exp.clone(), sum);
                }
            }
        }
        coeffs
    }

    /// Drop variable 0, which must not occur in any term. Returns a
    /// polynomial in `n_vars - 1` variables with indices shifted down.
    ///
    /// Panics in debug builds if variable 0 occurs with a non-zero exponent.
    pub fn drop_main_var(&self) -> Self {
        debug_assert!(
            self.terms_ref()
                .keys()
                .all(|e| e.first().copied().unwrap_or(0) == 0),
            "drop_main_var: variable 0 must not occur"
        );
        let new_n_vars = self.n_vars.saturating_sub(1);
        let mut result = Self::new(self.domain.clone(), new_n_vars);
        for (exp, coeff) in &self.terms {
            if exp.first().copied().unwrap_or(0) != 0 {
                continue;
            }
            let new_exp: SmallVec<[usize; 4]> = exp.iter().skip(1).copied().collect();
            result.terms.insert(new_exp, coeff.clone());
        }
        result
    }

    /// Embed into one more variable by inserting a new variable 0 with
    /// exponent 0 (all existing variable indices shift up by one).
    pub fn embed_new_main(&self) -> Self {
        let new_n_vars = self.n_vars + 1;
        let mut result = Self::new(self.domain.clone(), new_n_vars);
        for (exp, coeff) in &self.terms {
            let mut new_exp = SmallVec::with_capacity(new_n_vars);
            new_exp.push(0);
            new_exp.extend(exp.iter().copied());
            result.terms.insert(new_exp, coeff.clone());
        }
        result
    }

    /// Permute variables: the result's exponent at position `i` is the old
    /// exponent at position `perm[i]`. `perm` must be a permutation of
    /// `0..n_vars`.
    pub fn permute_variables(&self, perm: &[usize]) -> Self {
        assert_eq!(perm.len(), self.n_vars, "perm must be a permutation");
        let mut result = Self::new(self.domain.clone(), self.n_vars);
        for (exp, coeff) in &self.terms {
            let mut new_exp = SmallVec::with_capacity(self.n_vars);
            for &p in perm {
                new_exp.push(exp.get(p).copied().unwrap_or(0));
            }
            result.terms.insert(new_exp, coeff.clone());
        }
        result
    }

    /// Divide this polynomial by `divisor`, returning the quotient only if
    /// the division is exact (zero remainder).
    pub fn checked_div_exact(&self, divisor: &Self) -> Option<Self> {
        if divisor.is_zero() {
            return None;
        }
        let (quot, rem) = self.div_rem_sparse(divisor);
        if rem.is_zero() { Some(quot) } else { None }
    }

    /// Evaluate variable `var_index` at `value` while keeping the total
    /// number of variables unchanged (the variable disappears from the
    /// support but all indices are preserved).
    ///
    /// This is the substitution used by multivariate Hensel lifting, where
    /// variable positions must stay fixed across recursion levels.
    pub fn eval_keep(&self, var_index: usize, value: &D::Element) -> Self {
        let mut result = Self::new(self.domain.clone(), self.n_vars);
        for (exp, coeff) in &self.terms {
            let power = self.domain.pow(value, exp[var_index] as u64);
            let new_coeff = self.domain.mul(coeff, &power);
            if self.domain.is_zero(&new_coeff) {
                continue;
            }
            let mut new_exp = exp.clone();
            new_exp[var_index] = 0;
            let existing = result
                .terms
                .get(&new_exp)
                .cloned()
                .unwrap_or_else(|| self.domain.zero());
            let sum = self.domain.add(&existing, &new_coeff);
            if self.domain.is_zero(&sum) {
                result.terms.remove(&new_exp);
            } else {
                result.terms.insert(new_exp, sum);
            }
        }
        result
    }

    // ------------------------------------------------------------------
    //  F4 / Gröbner support helpers
    // ------------------------------------------------------------------

    /// Return the exponent vector of the leading monomial, or `None` for zero.
    ///
    /// This is an alias for [`leading_monomial`](Self::leading_monomial) that
    /// matches the Symbolica naming convention used in the F4 algorithm.
    #[inline]
    pub fn max_exp(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.leading_monomial()
    }

    /// Return the leading coefficient, or `None` for zero.
    ///
    /// This is an alias for [`leading_coeff`](Self::leading_coeff) that
    /// matches the Symbolica naming convention used in the F4 algorithm.
    #[inline]
    pub fn max_coeff(&self) -> Option<&D::Element> {
        self.leading_coeff()
    }

    /// Iterate over all exponent vectors in sorted order (descending by
    /// the monomial ordering).
    ///
    /// The F4 algorithm uses this to enumerate every monomial in a
    /// polynomial for symbolic preprocessing.
    pub fn exponents_iter(&self) -> impl Iterator<Item = &SmallVec<[usize; 4]>> {
        let mut sorted: Vec<_> = self.terms.keys().collect();
        sorted.sort_by(|a, b| self.order.cmp(a, b));
        sorted.into_iter()
    }

    /// Divide every term by the leading coefficient, making the polynomial
    /// monic. Returns `false` if the polynomial is zero or the leading
    /// coefficient has no inverse.
    pub fn make_monic_inplace(&mut self) -> bool {
        if self.is_zero() {
            return false;
        }
        let lc = self.leading_coeff().cloned().unwrap();
        match self.domain.inv(&lc) {
            Some(inv_lc) => {
                for coeff in self.terms.values_mut() {
                    *coeff = self.domain.mul(coeff, &inv_lc);
                }
                true
            }
            None => false,
        }
    }

    /// Create a zero polynomial with the same domain and variable count.
    ///
    /// This is identical to [`zero`](Self::zero) but named to match the
    /// Symbolica convention used in F4 code.
    #[inline]
    pub fn zero_with_capacity(&self, _cap: usize) -> Self {
        self.zero()
    }

    /// Append a single monomial term `coeff * x^exp`.
    ///
    /// If the monomial already exists, the coefficients are summed.
    /// Zero coefficients remove the term.
    pub fn append_monomial(&mut self, coeff: D::Element, exp: &[usize]) {
        let key = Self::normalize_exp(exp, self.n_vars);
        let existing = self
            .terms
            .get(&key)
            .cloned()
            .unwrap_or_else(|| self.domain.zero());
        let sum = self.domain.add(&existing, &coeff);
        if self.domain.is_zero(&sum) {
            self.terms.remove(&key);
        } else {
            self.terms.insert(key, sum);
        }
    }

    /// Evaluate the polynomial by substituting `value` for variable `var_index`.
    ///
    /// Returns a polynomial in one fewer variable (all remaining variables
    /// keep their relative order). If `var_index` is the only variable, the
    /// result is a zero-variable (constant) polynomial.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::SparseMultivariatePolynomial;
    /// use ocas_poly::Lex;
    ///
    /// let p = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    ///     IntegerDomain, 2,
    ///     vec![
    ///         (vec![1, 1], Integer::from(1)), // x*y
    ///         (vec![0, 1], Integer::from(2)), // 2*y
    ///     ],
    /// );
    /// // Substitute x=3: result = 3*y + 2*y = 5*y
    /// let r = p.eval(0, &Integer::from(3));
    /// assert_eq!(r.coeff(&[1]), Integer::from(5));
    /// ```
    pub fn eval(&self, var_index: usize, value: &D::Element) -> Self {
        let new_n_vars = self.n_vars.saturating_sub(1);
        let mut result = Self::new(self.domain.clone(), new_n_vars);
        for (exp, coeff) in &self.terms {
            // Compute coefficient * value^exp[var_index].
            let power = self.domain.pow(value, exp[var_index] as u64);
            let new_coeff = self.domain.mul(coeff, &power);
            if self.domain.is_zero(&new_coeff) {
                continue;
            }
            // Build new exponent vector without var_index.
            let mut new_exp = SmallVec::with_capacity(new_n_vars);
            for i in 0..self.n_vars {
                if i != var_index {
                    new_exp.push(exp[i]);
                }
            }
            let existing = result
                .terms
                .get(&new_exp)
                .cloned()
                .unwrap_or_else(|| self.domain.zero());
            let sum = self.domain.add(&existing, &new_coeff);
            if self.domain.is_zero(&sum) {
                result.terms.remove(&new_exp);
            } else {
                result.terms.insert(new_exp, sum);
            }
        }
        result
    }
}

// ------------------------------------------------------------------
//  Factorization entry points for sparse multivariate polynomials
// ------------------------------------------------------------------

impl SparseMultivariatePolynomial<IntegerDomain, Lex> {
    /// Factor this bivariate integer polynomial into irreducible factors with
    /// multiplicities.
    ///
    /// With a constant leading coefficient in $x$ the polynomial is treated
    /// as univariate in $x$ with coefficients in $\mathbb{Z}[y]$ and factored
    /// via Wang's bivariate Hensel-lifting algorithm. With a non-constant
    /// leading coefficient the general EEZ path with Wang leading-coefficient
    /// imposition (p-adic coefficient Hensel lifting) is used instead.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{Integer, IntegerDomain};
    /// use ocas_poly::SparseMultivariatePolynomial;
    /// use ocas_poly::Lex;
    ///
    /// // (x^2 + y + 1)(x + y + 2)
    /// let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
    ///     IntegerDomain, 2,
    ///     vec![
    ///         (vec![3, 0], Integer::from(1)),
    ///         (vec![2, 1], Integer::from(1)),
    ///         (vec![2, 0], Integer::from(2)),
    ///         (vec![1, 1], Integer::from(1)),
    ///         (vec![1, 0], Integer::from(1)),
    ///         (vec![0, 2], Integer::from(1)),
    ///         (vec![0, 1], Integer::from(3)),
    ///         (vec![0, 0], Integer::from(2)),
    ///     ],
    /// );
    /// let factors = f.factor();
    /// assert!(factors.len() >= 2);
    /// ```
    pub fn factor(&self) -> Vec<(Self, usize)> {
        if self.n_vars() >= 3 {
            crate::factor::eez::multivariate_factor_z(self)
        } else if self
            .leading_coeff_in(0)
            .terms_ref()
            .keys()
            .any(|e| e.iter().skip(1).any(|&d| d > 0))
        {
            // Non-constant leading coefficient in x: the bivariate path
            // requires a constant LC, so use the EEZ path with Wang
            // leading-coefficient imposition.
            crate::factor::eez::multivariate_factor_z(self)
        } else {
            bivariate_factor_z(self, 0, 1)
        }
    }
}

impl SparseMultivariatePolynomial<FiniteField, Lex> {
    /// Factor this multivariate polynomial over a prime finite field into
    /// irreducible factors with multiplicities.
    ///
    /// Bivariate polynomials use the evaluation–Hensel path; polynomials in
    /// three or more variables use EEZ Hensel lifting. Both currently require
    /// the leading coefficient in the main variable to be a nonzero field
    /// constant.
    pub fn factor(&self) -> Vec<(Self, usize)> {
        if self.n_vars() >= 3 {
            crate::factor::eez::multivariate_factor_fp(self)
        } else {
            bivariate_factor_fp(self, 0, 1)
        }
    }
}

// ------------------------------------------------------------------
//  Monomial utilities
// ------------------------------------------------------------------

/// Check whether monomial `a` divides monomial `b`: `a[i] >= b[i]` for all i.
pub fn monomial_divides(a: &[usize], b: &[usize]) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| x >= y)
}

/// Compute the least common multiple of two monomials: element-wise max.
pub fn monomial_lcm(a: &[usize], b: &[usize]) -> SmallVec<[usize; 4]> {
    a.iter().zip(b.iter()).map(|(x, y)| *x.max(y)).collect()
}

/// Return true if the two monomials are coprime (no variable appears in both).
pub fn monomial_are_coprime(a: &[usize], b: &[usize]) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| *x == 0 || *y == 0)
}

/// Binomial coefficient `n choose k`.
pub(crate) fn binomial(n: usize, k: usize) -> u64 {
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
