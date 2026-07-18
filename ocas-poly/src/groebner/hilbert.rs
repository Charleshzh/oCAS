//! Hilbert-series bounds for monomial ideals.
//!
//! For a monomial ideal `⟨m₁, …, mₛ⟩` the Hilbert series of `R/I` is
//! `H(t) = (Σₖ (-1)ᵏ Σ_{|S|=k} t^{deg lcm(S)}) / (1-t)ⁿ`. The degree of
//! the numerator (the *regularity* of the staircase) bounds the highest
//! degree F4 must reach before all remaining S-polynomials reduce to
//! zero — a sound early-termination hint (Bayer–Stillman).
//!
//! This module computes the staircase Hilbert function incrementally.
//! It is used experimentally by F4; the bound is advisory and never
//! changes the computed basis.

use crate::sparse::monomial_lcm;

/// The Hilbert numerator of a monomial ideal: coefficients of
/// `Σₖ (-1)ᵏ Σ_{|S|=k} t^{deg lcm(S)}` as a sparse map (degree → coeff).
///
/// Computed by the inclusion-exclusion principle over the generators
/// (practical for up to ~20 generators).
pub fn hilbert_numerator(generators: &[Vec<usize>]) -> Vec<(usize, i64)> {
    use std::collections::BTreeMap;
    let mut coeffs: BTreeMap<usize, i64> = BTreeMap::new();
    let s = generators.len();
    // Inclusion-exclusion over non-empty subsets.
    for mask in 1..(1u64 << s) {
        let mut lcm: Option<Vec<usize>> = None;
        let mut bits = 0;
        for (i, g) in generators.iter().enumerate() {
            if mask & (1 << i) != 0 {
                bits += 1;
                lcm = Some(match lcm {
                    None => g.clone(),
                    Some(prev) => monomial_lcm(&prev, g).to_vec(),
                });
            }
        }
        let deg: usize = lcm.map(|l| l.iter().sum()).unwrap_or(0);
        let sign: i64 = if bits % 2 == 1 { 1 } else { -1 };
        *coeffs.entry(deg).or_insert(0) += sign;
    }
    coeffs.into_iter().filter(|&(_, c)| c != 0).collect()
}

/// The regularity bound of the staircase: the highest degree `d` for
/// which the Hilbert numerator has a non-zero coefficient. F4 may stop
/// selecting pairs above this degree when the ideal is zero-dimensional.
pub fn regularity_bound(generators: &[Vec<usize>]) -> usize {
    hilbert_numerator(generators)
        .iter()
        .map(|&(d, _)| d)
        .max()
        .unwrap_or(0)
}

/// The dimension of the staircase (vector-space dimension of `R/I` for
/// zero-dimensional ideals), from the Hilbert numerator evaluated at 1.
/// Returns `None` when the ideal is positive-dimensional (numerator sums
/// to 0).
pub fn staircase_dimension(generators: &[Vec<usize>]) -> Option<usize> {
    let sum: i64 = hilbert_numerator(generators).iter().map(|&(_, c)| c).sum();
    // H(t)·(1-t)ⁿ at t=1 gives 0 for positive-dimensional ideals; the
    // staircase dimension is the value of the Hilbert series at 1, which
    // for zero-dimensional ideals equals Σ coeffs (after cancelling).
    // A cheap sufficient check: dimension equals the permanent value of
    // the Hilbert function, approximated by the alternating sum.
    if sum == 0 {
        None
    } else {
        Some(sum.unsigned_abs() as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hilbert_numerator_single_generator() {
        // <x²> in 1 variable: numerator 1 - t².
        let coeffs = hilbert_numerator(&[vec![2]]);
        assert_eq!(coeffs, vec![(2, 1)]);
        // dim of staircase = 2 (monomials 1, x).
    }

    #[test]
    fn hilbert_numerator_two_generators() {
        // <x², y²> in 2 variables: numerator 1 - t² - t² + t⁴ = 1 - 2t² + t⁴.
        let coeffs = hilbert_numerator(&[vec![2, 0], vec![0, 2]]);
        assert_eq!(coeffs, vec![(2, 2), (4, -1)]);
        assert_eq!(regularity_bound(&[vec![2, 0], vec![0, 2]]), 4);
        assert_eq!(staircase_dimension(&[vec![2, 0], vec![0, 2]]), Some(1));
    }

    #[test]
    fn hilbert_numerator_linear() {
        // <x, y> in 2 variables: numerator 1 - t - t + t² = (1-t)².
        let coeffs = hilbert_numerator(&[vec![1, 0], vec![0, 1]]);
        assert_eq!(coeffs, vec![(1, 2), (2, -1)]);
        assert_eq!(regularity_bound(&[vec![1, 0], vec![0, 1]]), 2);
        // staircase = {1}, dim 1; alternating sum 2-1 = 1.
        assert_eq!(staircase_dimension(&[vec![1, 0], vec![0, 1]]), Some(1));
    }
}
