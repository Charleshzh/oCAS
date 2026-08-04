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
//!   (e.g. cyclic-n). Production-grade since 0.19.0.
//! - **MultiModular** ([`multi_modular`]) — multi-prime strategy for ℚ
//!   ideals since 0.25.0: parallel F5 over lucky primes, CRT + rational
//!   reconstruction, exact ℚ verification, and a p-adic Hensel-lift
//!   shortcut.
//!
//! All algorithms produce a reduced Gröbner basis. [`Algorithm::Auto`]
//! routes ℚ ideals through the multi-modular pipeline and other domains
//! through F4.

pub mod f4;
pub mod f5;
pub mod fglm;
pub mod hilbert;
pub mod multi_modular;
pub(crate) mod packed;

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
    /// Automatically select the most suitable algorithm: the multi-modular
    /// pipeline for ℚ ideals (fast path since 0.25.0), F4 otherwise.
    #[default]
    Auto,
    /// Force the F4 matrix algorithm (Faugère 1999).
    F4,
    /// Force the F5 signature-based algorithm (Faugère 2002).
    F5,
    /// Force Buchberger's classic S-polynomial iteration.
    Buchberger,
    /// Force the multi-modular pipeline ([`crate::groebner::multi_modular`]):
    /// parallel F5 over lucky primes + CRT/rational reconstruction + exact
    /// ℚ verification, with a Hensel-lift shortcut. Only applies to ℚ
    /// coefficients; other domains fall back to F4.
    MultiModular,
}

/// Compute a Gröbner basis using the requested [`Algorithm`].
///
/// This is the unified entry point for Gröbner basis computation. Zero
/// polynomials in `ideal` are filtered internally by each backend.
///
/// [`Algorithm::Auto`] routes ℚ ideals through the multi-modular pipeline
/// (the fast path for rational coefficients since 0.25.0) and other
/// domains through F4.
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
pub fn groebner_basis<D: Domain + 'static, O: MonomialOrder + Send + Sync>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    algo: Algorithm,
) -> GroebnerBasis<D, O> {
    match algo {
        // Auto: multi-modular for ℚ ideals (the internal Any check returns
        // None for other domains, which then take the F4 path).
        Algorithm::Auto => match multi_modular::groebner_basis_mm(ideal) {
            Some(gb) => gb,
            None => f4::f4(ideal),
        },
        Algorithm::F4 => f4::f4(ideal),
        Algorithm::F5 => f5::f5(ideal),
        Algorithm::Buchberger => buchberger(ideal),
        Algorithm::MultiModular => match multi_modular::groebner_basis_mm(ideal) {
            Some(gb) => gb,
            None => f4::f4(ideal),
        },
    }
}

/// Eliminate variables from an ideal.
///
/// Returns the Gröbner basis of `I ∩ k[x_{elim_vars}, ..., x_{n-1}]`, i.e.,
/// the polynomials in the basis that do not involve the first `elim_vars`
/// variables. Uses Lex ordering which is a natural elimination order:
/// under Lex, the reduced GB of an ideal automatically contains the
/// elimination ideal's generators.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::Lex;
/// use ocas_poly::{SparseMultivariatePolynomial, eliminate, Algorithm};
///
/// let d = RationalDomain;
/// // Ideal: x + y + z, x*y + x*z in k[x,y,z]; eliminate x.
/// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
///     (vec![1, 0, 0], Rational::new(1, 1)),
///     (vec![0, 1, 0], Rational::new(1, 1)),
///     (vec![0, 0, 1], Rational::new(1, 1)),
/// ]);
/// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 3, vec![
///     (vec![1, 1, 0], Rational::new(1, 1)),
///     (vec![1, 0, 1], Rational::new(1, 1)),
/// ]);
/// let elim = eliminate(&[f1, f2], 1, Algorithm::Auto);
/// // Result should be in k[y,z]
/// for p in &elim.basis {
///     assert!(p.degree_in(0) == 0, "eliminated variable x should not appear");
/// }
/// ```
pub fn eliminate<D: Domain + 'static>(
    ideal: &[SparseMultivariatePolynomial<D, crate::sparse::Lex>],
    elim_vars: usize,
    algo: Algorithm,
) -> GroebnerBasis<D, crate::sparse::Lex> {
    let n_vars = ideal.first().map(|p| p.n_vars()).unwrap_or(0);
    assert!(
        elim_vars <= n_vars,
        "elim_vars ({elim_vars}) must be <= n_vars ({n_vars})"
    );
    if ideal.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Compute Gröbner basis under Lex ordering.
    // Lex is a natural elimination order: polynomials in the GB that
    // don't involve x_0,...,x_{s-1} form a GB of the elimination ideal.
    let gb = groebner_basis(ideal, algo);

    // Filter: keep only polynomials that don't involve the eliminated variables.
    let filtered: Vec<SparseMultivariatePolynomial<D, crate::sparse::Lex>> = gb
        .basis
        .into_iter()
        .filter(|p| {
            p.terms_ref()
                .keys()
                .all(|exp| exp.iter().take(elim_vars).all(|&e| e == 0))
        })
        .collect();

    GroebnerBasis { basis: filtered }
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

    // --- Step 1a: Lex order verification ---

    #[test]
    fn lex_cyclic_3() {
        // Cyclic-3: x+y+z, xy+yz+zx, xyz-1 under Lex.
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![
                (vec![1, 0, 0], r(1, 1)),
                (vec![0, 1, 0], r(1, 1)),
                (vec![0, 0, 1], r(1, 1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![
                (vec![1, 1, 0], r(1, 1)),
                (vec![0, 1, 1], r(1, 1)),
                (vec![1, 0, 1], r(1, 1)),
            ],
        );
        let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![(vec![1, 1, 1], r(1, 1)), (vec![0, 0, 0], r(-1, 1))],
        );
        let gb = groebner_basis(&[f1, f2, f3], Algorithm::F4);
        assert!(gb.is_groebner_basis());
        // Print basis for debugging.
        for (i, p) in gb.basis.iter().enumerate() {
            eprintln!("lex_cyclic_3 gb[{i}]: {p:?}");
        }
        // Under Lex, the GB should be triangular (each poly introduces
        // one fewer variable). The smallest variable (z) should appear
        // in a univariate polynomial.
        let has_univariate_in_z = gb
            .basis
            .iter()
            .any(|p| p.terms_ref().keys().all(|e| e[0] == 0 && e[1] == 0));
        assert!(
            has_univariate_in_z,
            "Lex GB should contain a univariate poly in z"
        );
    }

    #[test]
    fn lex_two_variable_elimination() {
        // Ideal: x^2 - y, x^3 - x in k[x,y] under Lex.
        // Lex GB should eliminate x: expect y^2 - y, xy - x, x^2 - y.
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![2, 0], r(1, 1)), (vec![0, 1], r(-1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![3, 0], r(1, 1)), (vec![1, 0], r(-1, 1))],
        );
        let gb = groebner_basis(&[f1, f2], Algorithm::F4);
        assert!(gb.is_groebner_basis());
        // The GB should be triangular: first poly in y only, then xy, then x^2.
        // Find the univariate poly in y.
        let y_poly = gb
            .basis
            .iter()
            .find(|p| p.terms_ref().keys().all(|e| e[0] == 0));
        assert!(
            y_poly.is_some(),
            "Lex GB should contain a univariate poly in y"
        );
    }

    // --- Step 1c: eliminate() tests ---

    #[test]
    fn eliminate_simple() {
        // Eliminate x from {x + y, x - y} in k[x,y].
        // x + y = 0 and x - y = 0 ⟹ x = 0, y = 0.
        // Eliminating x should give {y} (or just y = 0).
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 0], r(1, 1)), (vec![0, 1], r(1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 0], r(1, 1)), (vec![0, 1], r(-1, 1))],
        );
        let elim = eliminate(&[f1, f2], 1, Algorithm::F4);
        assert!(!elim.basis.is_empty());
        // All result polynomials should be in y only.
        for p in &elim.basis {
            assert_eq!(p.degree_in(0), 0, "eliminated var x should not appear");
        }
    }

    #[test]
    fn eliminate_cox_little_oshea() {
        // Cox-Little-O'Shea §3.1: eliminate x from {x+y+z-1, xy+xz, xyz}.
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![
                (vec![1, 0, 0], r(1, 1)),
                (vec![0, 1, 0], r(1, 1)),
                (vec![0, 0, 1], r(1, 1)),
                (vec![0, 0, 0], r(-1, 1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![(vec![1, 1, 0], r(1, 1)), (vec![1, 0, 1], r(1, 1))],
        );
        let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![(vec![1, 1, 1], r(1, 1))],
        );
        let elim = eliminate(&[f1, f2, f3], 1, Algorithm::F4);
        // All result polynomials should be in y, z only.
        for p in &elim.basis {
            assert_eq!(p.degree_in(0), 0, "x should be eliminated");
        }
        // Should contain y^2 + z^2 - y - z and yz + z^2 - z (or equivalent).
        assert!(
            elim.basis.len() >= 2,
            "expected at least 2 generators, got {}",
            elim.basis.len()
        );
    }
}
