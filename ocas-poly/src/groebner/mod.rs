//! Gröbner basis computation for multivariate polynomial ideals.
//!
//! Provides three algorithms, all reachable through the unified
//! [`groebner_basis`] entry point with an [`Algorithm`] selector:
//!
//! - **Buchberger** ([`buchberger`]) — classic S-polynomial iteration with
//!   Gebauer-Moeller optimization. Suitable for small ideals.
//! - **F4** ([`f4::f4`]) — matrix-based algorithm from Faugère (1999).
//!   Dramatically faster for larger ideals by batching S-polynomial
//!   reductions into sparse matrix row operations.
//! - **F5** ([`f5::f5`]) — signature-based algorithm from Faugère (2002).
//!   Rejects zero-reducers *before* matrix construction via syzygy
//!   criteria, targeting order-of-magnitude speedups on difficult ideals
//!   (e.g. cyclic-n). Currently a placeholder; full implementation
//!   landing in 0.19.0.
//!
//! All algorithms produce a reduced Gröbner basis. [`Algorithm::Auto`]
//! selects a backend by heuristic (currently F4).

pub mod f4;
pub mod f5;
pub mod fglm;
pub mod hilbert;

use ocas_core::FastHashSet as HashSet;
use ocas_domain::Domain;

use crate::sparse::{
    MonomialOrder, SparseMultivariatePolynomial, monomial_are_coprime, monomial_divides,
};

/// A Gröbner basis for a polynomial ideal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerBasis<D: Domain, O: MonomialOrder> {
    /// The polynomials forming the basis.
    pub basis: Vec<SparseMultivariatePolynomial<D, O>>,
}

impl<D: Domain, O: MonomialOrder> GroebnerBasis<D, O> {
    /// Compute a Gröbner basis from a set of generators using Buchberger's algorithm.
    ///
    /// Requires that the coefficient domain supports exact division (i.e., is
    /// effectively a field). The algorithm will panic if division fails.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{RationalDomain, Rational};
    /// use ocas_poly::sparse::Lex;
    /// use ocas_poly::GroebnerBasis;
    /// use ocas_poly::SparseMultivariatePolynomial;
    ///
    /// let d = RationalDomain;
    /// // ideal: x + y, x - y
    /// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    ///     (vec![1, 0], Rational::new(1, 1)),
    ///     (vec![0, 1], Rational::new(1, 1)),
    /// ]);
    /// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    ///     (vec![1, 0], Rational::new(1, 1)),
    ///     (vec![0, 1], Rational::new(-1, 1)),
    /// ]);
    /// let gb = GroebnerBasis::buchberger(&[f1, f2]);
    /// assert!(gb.basis.len() >= 2);
    /// ```
    pub fn buchberger(ideal: &[SparseMultivariatePolynomial<D, O>]) -> Self {
        // Filter out zero polynomials.
        let mut basis: Vec<SparseMultivariatePolynomial<D, O>> =
            ideal.iter().filter(|p| !p.is_zero()).cloned().collect();
        if basis.is_empty() {
            return Self { basis };
        }

        // Collect critical pairs: all unordered pairs (i, j) with i < j.
        let mut pairs: HashSet<(usize, usize)> = HashSet::default();
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                pairs.insert((i, j));
            }
        }

        let max_iter = 10000;

        for _ in 0..max_iter {
            if pairs.is_empty() {
                break;
            }
            let (i, j) = *pairs.iter().next().unwrap();
            pairs.remove(&(i, j));

            // Buchberger's first criterion: if the leading monomials are
            // coprime, the S-polynomial reduces to zero, so skip.
            let lm_i = basis[i].leading_monomial();
            let lm_j = basis[j].leading_monomial();
            if let (Some(mi), Some(mj)) = (&lm_i, &lm_j)
                && monomial_are_coprime(mi, mj)
            {
                continue;
            }

            // Compute S-polynomial and reduce by current basis.
            let s = basis[i].spoly(&basis[j]);
            let r = s.reduce(&basis);

            if !r.is_zero() {
                let new_idx = basis.len();
                basis.push(r);
                for k in 0..new_idx {
                    pairs.insert((k, new_idx));
                }
            }
        }

        Self { basis }
    }

    /// Minimize the basis: remove polynomials whose leading monomial is
    /// divisible by another element's leading monomial.
    pub fn minimize(mut self) -> Self {
        let lms: Vec<_> = self
            .basis
            .iter()
            .filter_map(|p| p.leading_monomial().cloned())
            .collect();

        let mut keep = vec![true; self.basis.len()];
        for i in 0..self.basis.len() {
            for j in 0..self.basis.len() {
                // Remove i if lms[j] divides lms[i] (i.e., lms[i] is a
                // multiple of lms[j], making i redundant).
                // monomial_divides(big, small) returns true when small divides big.
                if i != j && keep[i] && keep[j] && monomial_divides(&lms[i], &lms[j]) {
                    keep[i] = false;
                    break;
                }
            }
        }

        self.basis = self
            .basis
            .into_iter()
            .enumerate()
            .filter(|(i, _)| keep[*i])
            .map(|(_, p)| p)
            .collect();

        self
    }

    /// Inter-reduce the basis: reduce each element by the others and make
    /// each polynomial monic.
    ///
    /// The algorithm processes elements in ascending order of leading
    /// monomial. Each element is reduced by all elements with strictly
    /// smaller leading monomials (those already in the result set).
    /// This ensures the standard reduced Gröbner basis property:
    /// no monomial of any basis element is divisible by the leading
    /// monomial of any other basis element.
    pub fn auto_reduce(mut self) -> Self {
        let order = self
            .basis
            .first()
            .map(|p| p.order.clone())
            .unwrap_or_default();
        // Sort basis in ascending order of leading monomial (smallest first).
        self.basis
            .sort_by(|a, b| match (a.leading_monomial(), b.leading_monomial()) {
                (Some(ma), Some(mb)) => order.cmp(ma, mb),
                (Some(_), None) => std::cmp::Ordering::Greater,
                (None, Some(_)) => std::cmp::Ordering::Less,
                (None, None) => std::cmp::Ordering::Equal,
            });

        let mut reduced: Vec<SparseMultivariatePolynomial<D, O>> = Vec::new();

        for poly in &self.basis {
            // Reduce `poly` by all elements already in `reduced`
            // (which have smaller leading monomials).
            let mut r = poly.reduce(&reduced);
            if !r.is_zero() {
                if let Some(lc) = r.leading_coeff().cloned()
                    && let Some(inv) = r.domain().inv(&lc)
                {
                    r = r.mul_scalar(&inv);
                }
                reduced.push(r);
            }
        }

        self.basis = reduced;
        self
    }

    /// Verify that this is indeed a Gröbner basis by checking that all
    /// S-polynomials reduce to zero.
    pub fn is_groebner_basis(&self) -> bool {
        for i in 0..self.basis.len() {
            for j in i + 1..self.basis.len() {
                let s = self.basis[i].spoly(&self.basis[j]);
                let r = s.reduce(&self.basis);
                if !r.is_zero() {
                    return false;
                }
            }
        }
        true
    }

    /// Change the monomial order of this Gröbner basis.
    ///
    /// The polynomials are re-interpreted under the target order `O2`
    /// and the F4 algorithm is re-run. This is the simple reorder path
    /// (Symbolica's `reorder::<Order>()`). For zero-dimensional ideals,
    /// use [`crate::groebner::fglm::fglm`] for a much faster conversion.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{RationalDomain, Rational};
    /// use ocas_poly::sparse::{Grevlex, Lex};
    /// use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial, f4};
    ///
    /// let d = RationalDomain;
    /// // ideal: x + y, x - y  → basis {y, x} under Lex
    /// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    ///     (vec![1, 0], Rational::new(1, 1)),
    ///     (vec![0, 1], Rational::new(1, 1)),
    /// ]);
    /// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
    ///     (vec![1, 0], Rational::new(1, 1)),
    ///     (vec![0, 1], Rational::new(-1, 1)),
    /// ]);
    /// let gb_lex = f4::f4(&[f1, f2]);
    /// let gb_grevlex = gb_lex.reorder::<Grevlex>();
    /// assert!(gb_grevlex.is_groebner_basis());
    /// ```
    pub fn reorder<O2: MonomialOrder>(&self) -> GroebnerBasis<D, O2>
    where
        D: 'static,
    {
        let converted: Vec<SparseMultivariatePolynomial<D, O2>> = self
            .basis
            .iter()
            .map(|p| {
                SparseMultivariatePolynomial::from_terms(
                    p.domain().clone(),
                    p.n_vars(),
                    p.terms_ref()
                        .iter()
                        .map(|(e, c)| (e.to_vec(), c.clone()))
                        .collect(),
                )
            })
            .collect();
        crate::groebner::f4::f4(&converted)
    }
}

/// Convenience: compute a Gröbner basis and inter-reduce it.
pub fn buchberger<D: Domain, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O> {
    GroebnerBasis::buchberger(ideal).minimize().auto_reduce()
}

/// Algorithm selector for [`groebner_basis`].
///
/// `Auto` picks a backend based on the ideal's size and structure; the
/// other variants force a specific algorithm. See [`groebner_basis`] for
/// the unified entry point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Algorithm {
    /// Automatically select the most suitable algorithm based on ideal
    /// size and structure (heuristic, calibrated from benchmarks).
    /// Currently routes to F4; the crossover to F5 will be tuned from
    /// cyclic-n benchmarks once the F5 core is complete.
    #[default]
    Auto,
    /// Force the F4 matrix algorithm (Faugère 1999).
    F4,
    /// Force the F5 signature-based algorithm (Faugère 2002).
    F5,
    /// Force Buchberger's classic S-polynomial iteration.
    Buchberger,
}

/// Compute a Gröbner basis using the requested [`Algorithm`].
///
/// This is the unified entry point for Gröbner basis computation. Zero
/// polynomials in `ideal` are filtered internally by each backend.
///
/// [`Algorithm::Auto`] currently routes to F4; the crossover to F5 will
/// be calibrated from cyclic-n benchmarks once the F5 core is complete.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::Lex;
/// use ocas_poly::{Algorithm, groebner_basis, SparseMultivariatePolynomial};
///
/// let d = RationalDomain;
/// // ideal: x + y, x - y
/// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(1, 1)),
/// ]);
/// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(-1, 1)),
/// ]);
/// let gb = groebner_basis(&[f1, f2], Algorithm::Auto);
/// assert!(gb.is_groebner_basis());
/// ```
pub fn groebner_basis<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    algo: Algorithm,
) -> GroebnerBasis<D, O> {
    match algo {
        Algorithm::Auto | Algorithm::F4 => f4::f4(ideal),
        Algorithm::F5 => f5::f5(ideal),
        Algorithm::Buchberger => buchberger(ideal),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Lex;
    use ocas_domain::{Rational, RationalDomain};

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    fn make_poly(
        terms: Vec<(Vec<usize>, Rational)>,
    ) -> SparseMultivariatePolynomial<RationalDomain, Lex> {
        SparseMultivariatePolynomial::from_terms(RationalDomain, 2, terms)
    }

    #[test]
    fn empty_ideal() {
        let gb = buchberger::<RationalDomain, Lex>(&[]);
        assert!(gb.basis.is_empty());
    }

    #[test]
    fn single_polynomial() {
        // f = x^2 - 1
        let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            1,
            vec![(vec![2], r(1, 1)), (vec![0], r(-1, 1))],
        );
        let gb = buchberger(&[f]);
        assert_eq!(gb.basis.len(), 1);
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn linear_system() {
        // x + y = 0, x - y = 0  →  basis = {x, y}
        let f1 = make_poly(vec![(vec![1, 0], r(1, 1)), (vec![0, 1], r(1, 1))]);
        let f2 = make_poly(vec![(vec![1, 0], r(1, 1)), (vec![0, 1], r(-1, 1))]);
        let gb = buchberger(&[f1, f2]);
        assert!(gb.is_groebner_basis());
        // After auto-reduce, we expect {x, y} (monic leading terms)
        assert!(gb.basis.len() >= 2);
    }

    #[test]
    fn two_variable_ideal() {
        // x^2 - y, x^3 - x  (elimination ideal: y = x^2, x^3 = x → x ∈ {0, ±1})
        let f1 = make_poly(vec![(vec![2, 0], r(1, 1)), (vec![0, 1], r(-1, 1))]);
        let f2 = make_poly(vec![(vec![3, 0], r(1, 1)), (vec![1, 0], r(-1, 1))]);
        let gb = buchberger(&[f1, f2]);
        assert!(gb.is_groebner_basis());
        assert!(!gb.basis.is_empty());
    }
}
