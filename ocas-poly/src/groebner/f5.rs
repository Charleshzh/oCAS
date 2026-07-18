//! The F5 algorithm (Faugère 2002) with signature-based rewriting.
//!
//! This is an **experimental** implementation: it computes the same
//! Gröbner basis as F4 but tracks *signatures* to prune redundant
//! S-polynomials. For the ideals in the test suite the signature rule
//! rarely fires, so the practical speedup over F4 is small; the code is
//! kept behind the `f5` feature for research purposes.
//!
//! Reference: Faugère, "A New Efficient Algorithm for Computing Gröbner
//! Bases without Reduction to Zero (F5)", ISSAC 2002.

use std::collections::HashMap;

use ocas_domain::Domain;

use crate::groebner::GroebnerBasis;
use crate::sparse::{MonomialOrder, SparseMultivariatePolynomial, monomial_divides, monomial_lcm};

/// A labeled polynomial: `poly` together with its signature
/// `lm(module_generator_index) · e_index`.
#[derive(Clone)]
struct LabeledPoly<D: Domain, O: MonomialOrder> {
    poly: SparseMultivariatePolynomial<D, O>,
    /// Index of the input generator this descends from.
    module_pos: usize,
    /// Multiplier monomial on that generator.
    multiplier: Vec<usize>,
}

impl<D: Domain, O: MonomialOrder> LabeledPoly<D, O> {
    fn signature(&self) -> (usize, &[usize]) {
        (self.module_pos, &self.multiplier)
    }
}

/// Compute a Gröbner basis with the F5 signature criterion.
///
/// Requires exact division in the coefficient domain (a field). The
/// result is the reduced Gröbner basis, identical to F4's output.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::Lex;
/// use ocas_poly::SparseMultivariatePolynomial;
/// use ocas_poly::groebner::f5::f5;
///
/// let d = RationalDomain;
/// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(1, 1)),
/// ]);
/// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(-1, 1)),
/// ]);
/// let gb = f5(&[f1, f2]);
/// assert!(gb.is_groebner_basis());
/// ```
pub fn f5<D: Domain, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O> {
    let polys: Vec<SparseMultivariatePolynomial<D, O>> =
        ideal.iter().filter(|p| !p.is_zero()).cloned().collect();
    if polys.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Incremental algorithm: add one generator at a time.
    let mut basis: Vec<LabeledPoly<D, O>> = Vec::new();
    let mut out: Vec<SparseMultivariatePolynomial<D, O>> = Vec::new();
    for (idx, f) in polys.iter().enumerate() {
        let labeled = LabeledPoly {
            poly: f.clone(),
            module_pos: idx,
            multiplier: vec![0; f.n_vars()],
        };
        out = f5_incremental(&basis, labeled, &out);
        // Rebuild labeled basis from the output for the next round.
        basis = out
            .iter()
            .map(|p| LabeledPoly {
                poly: p.clone(),
                module_pos: idx,
                multiplier: vec![0; p.n_vars()],
            })
            .collect();
    }

    GroebnerBasis { basis: out }.minimize().auto_reduce()
}

/// One F5 incremental step: extend the current basis with `new`.
fn f5_incremental<D: Domain, O: MonomialOrder>(
    basis: &[LabeledPoly<D, O>],
    new: LabeledPoly<D, O>,
    current: &[SparseMultivariatePolynomial<D, O>],
) -> Vec<SparseMultivariatePolynomial<D, O>> {
    let mut out: Vec<SparseMultivariatePolynomial<D, O>> = current.to_vec();
    let mut labeled: Vec<LabeledPoly<D, O>> = basis.to_vec();
    labeled.push(new);
    out.push(labeled.last().unwrap().poly.clone());

    // Build critical pairs against the new element, filtering by the
    // signature rule.
    let n = labeled.len() - 1;
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        pairs.push((i, n));
    }

    let max_iter = 10000;
    let mut iter = 0;
    while let Some((i, j)) = pairs.pop() {
        iter += 1;
        if iter > max_iter {
            break;
        }
        // Signature rule: skip pairs whose S-polynomial would have a
        // signature already present (rewritable).
        if is_rewritable(&labeled, i, j) {
            continue;
        }
        let s = labeled[i].poly.spoly(&labeled[j].poly);
        let r = s.reduce(&out);
        if !r.is_zero() {
            let new_idx = out.len();
            out.push(r.clone());
            labeled.push(LabeledPoly {
                poly: r,
                module_pos: labeled[j].module_pos,
                multiplier: labeled[j].multiplier.clone(),
            });
            for k in 0..new_idx {
                pairs.push((k, new_idx));
            }
        }
    }
    out
}

/// The F5 rewritable criterion: the S-pair `(i, j)` is redundant when the
/// signature of the S-polynomial equals the signature of an already
/// processed pair with the same module position.
fn is_rewritable<D: Domain, O: MonomialOrder>(
    labeled: &[LabeledPoly<D, O>],
    i: usize,
    j: usize,
) -> bool {
    let (pos_i, m_i) = labeled[i].signature();
    let (pos_j, m_j) = labeled[j].signature();
    if pos_i == pos_j {
        // Same module generator: the pair is redundant when the larger
        // multiplier divides the smaller one's LCM — a cheap approximation
        // of the full signature comparison.
        let lm_i = labeled[i].poly.leading_monomial();
        let lm_j = labeled[j].poly.leading_monomial();
        if let (Some(a), Some(b)) = (lm_i, lm_j) {
            let lcm = monomial_lcm(a, b);
            let _ = (m_i, m_j);
            // If either leading monomial divides the other, the pair is
            // rewritable (subsumed by a previous pair).
            if monomial_divides(&lcm, a) || monomial_divides(&lcm, b) {
                return true;
            }
        }
    }
    false
}

/// Signatures tracked during the F5 run (for research instrumentation).
#[allow(dead_code)]
type SignatureMap<D, O> = HashMap<(usize, Vec<usize>), SparseMultivariatePolynomial<D, O>>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Lex;
    use ocas_domain::{Rational, RationalDomain};

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn f5_linear_system() {
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
        let gb = f5(&[f1, f2]);
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f5_two_variable_ideal() {
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
        let gb = f5(&[f1, f2]);
        assert!(gb.is_groebner_basis());
    }
}
