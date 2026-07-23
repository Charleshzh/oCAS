//! FGLM Gröbner-basis conversion for zero-dimensional ideals.
//!
//! The FGLM algorithm (Faugère, Gianni, Lazard, Mora 1993) converts a
//! Gröbner basis of a **zero-dimensional** ideal from one monomial order
//! to another in `O(n·D³)` field operations, where `D` is the vector-space
//! dimension of `R/I`. This is dramatically faster than re-running F4
//! under the target order for large zero-dimensional ideals.
//!
//! Reference: Faugère et al., "Efficient Computation of Zero-dimensional
//! Gröbner Bases by Change of Ordering", JSC 1993.

use ocas_core::FastHashMap as HashMap;
use ocas_domain::Domain;

use crate::groebner::GroebnerBasis;
use crate::sparse::{MonomialOrder, SparseMultivariatePolynomial};

/// Convert a zero-dimensional Gröbner basis to the target order `O2`.
///
/// Returns `None` when the ideal is not zero-dimensional (infinitely many
/// monomials under the staircase). The input basis must be reduced.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::{Grevlex, Lex};
/// use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial, f4};
/// use ocas_poly::groebner::fglm::fglm;
///
/// let d = RationalDomain;
/// // ideal: x + y, x - y  →  zero-dimensional
/// let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(1, 1)),
/// ]);
/// let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(d, 2, vec![
///     (vec![1, 0], Rational::new(1, 1)),
///     (vec![0, 1], Rational::new(-1, 1)),
/// ]);
/// let gb_lex = f4::f4(&[f1, f2]);
/// let gb_grevlex = fglm::<_, Grevlex>(&gb_lex).expect("zero-dimensional");
/// assert!(gb_grevlex.is_groebner_basis());
/// ```
pub fn fglm<D: Domain, O2: MonomialOrder>(
    gb: &GroebnerBasis<D, impl MonomialOrder>,
) -> Option<GroebnerBasis<D, O2>> {
    let n_vars = gb.basis.first()?.n_vars();
    let domain = gb.basis.first()?.domain().clone();

    // The staircase: monomials not divisible by any leading monomial.
    // Zero-dimensional ⟺ the staircase is finite.
    let lms: Vec<Vec<usize>> = gb
        .basis
        .iter()
        .filter_map(|p| p.leading_monomial().map(|m| m.to_vec()))
        .collect();
    let staircase = compute_staircase(&lms, n_vars)?;
    let dim = staircase.len();

    // Normal forms of staircase monomials are themselves (already reduced).
    // Multiplication matrices M_i: for each variable x_i and each staircase
    // monomial m, nf(x_i·m) expressed in the staircase basis.
    let mut mult_matrices: Vec<Vec<Vec<D::Element>>> = vec![Vec::new(); n_vars];
    for var in 0..n_vars {
        let mut mat = vec![vec![domain.zero(); dim]; dim];
        for (col, m) in staircase.iter().enumerate() {
            let mut xm = m.clone();
            xm[var] += 1;
            let nf = normal_form_monomial(&xm, gb, &staircase, &domain);
            for (row, coeff) in nf.into_iter().enumerate() {
                mat[row][col] = coeff;
            }
        }
        mult_matrices[var] = mat;
    }

    // FGLM main loop: walk monomials of the target order in increasing
    // order, computing their normal forms; when a linear dependency with
    // previous normal forms is found, emit a new basis polynomial.
    let mut new_basis: Vec<SparseMultivariatePolynomial<D, O2>> = Vec::new();
    let mut seen_nfs: Vec<Vec<D::Element>> = Vec::new(); // normal forms seen
    let mut seen_mons: Vec<Vec<usize>> = Vec::new(); // corresponding monomials
    let mut boundary: Vec<Vec<usize>> = vec![vec![0; n_vars]]; // B, start at 1
    let mut visited: HashMap<Vec<usize>, bool> = HashMap::default();

    let max_steps = dim * (dim + 1) + n_vars * 4;
    let mut steps = 0;

    while !boundary.is_empty() && steps < max_steps {
        steps += 1;
        // Take the smallest monomial under O2 from the boundary.
        let (pos, m) = boundary
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| O2::default().cmp(a, b))
            .map(|(i, a)| (i, a.clone()))?;
        boundary.remove(pos);
        if visited.contains_key(&m) {
            continue;
        }
        visited.insert(m.clone(), true);

        // Normal form of m (in staircase coordinates).
        let nf = normal_form_monomial(&m, gb, &staircase, &domain);

        // Is nf linearly dependent on the previously seen normal forms?
        if let Some(relation) = find_relation(&seen_nfs, &nf, &domain) {
            // Emit polynomial: m + Σ relation_i · seen_mons_i.
            let mut terms: Vec<(Vec<usize>, D::Element)> = vec![(m.clone(), domain.one())];
            for (coeff, mon) in relation.into_iter().zip(seen_mons.iter()) {
                if !domain.is_zero(&coeff) {
                    terms.push((mon.clone(), coeff));
                }
            }
            new_basis.push(SparseMultivariatePolynomial::from_terms(
                domain.clone(),
                n_vars,
                terms,
            ));
            // Multiples of m need not be considered (they are in the ideal
            // generated by the new leading monomials).
            mark_multiples(&mut visited, &m, n_vars, dim * 2);
        } else {
            // No dependency: m is a new staircase monomial of the target
            // order; add its neighbours to the boundary.
            seen_nfs.push(nf);
            seen_mons.push(m.clone());
            for var in 0..n_vars {
                let mut next = m.clone();
                next[var] += 1;
                if !visited.contains_key(&next) {
                    boundary.push(next);
                }
            }
        }
    }

    if new_basis.is_empty() {
        return None;
    }
    let out = GroebnerBasis { basis: new_basis };
    Some(out.minimize().auto_reduce())
}

/// Monomials not divisible by any leading monomial (the staircase).
/// Returns `None` when the set is infinite (positive-dimensional ideal).
fn compute_staircase(lms: &[Vec<usize>], n_vars: usize) -> Option<Vec<Vec<usize>>> {
    let mut staircase = Vec::new();
    // BFS from the unit monomial; a monomial is in the staircase iff no LM
    // divides it. Zero-dimensional ⟺ BFS terminates.
    let mut queue = vec![vec![0usize; n_vars]];
    let mut seen: HashMap<Vec<usize>, ()> = HashMap::default();
    let limit = 100_000; // safety bound
    while let Some(m) = queue.pop() {
        if seen.contains_key(&m) {
            continue;
        }
        seen.insert(m.clone(), ());
        if seen.len() > limit {
            return None; // almost surely positive-dimensional
        }
        if lms.iter().any(|lm| monomial_divides_big(lm, &m)) {
            continue; // divisible by some LM → not in staircase
        }
        for var in 0..n_vars {
            let mut next = m.clone();
            next[var] += 1;
            queue.push(next);
        }
        staircase.push(m);
    }
    Some(staircase)
}

/// Whether `big` is divisible by `lm` (i.e. `lm` divides `big`).
fn monomial_divides_big(lm: &[usize], big: &[usize]) -> bool {
    lm.iter().zip(big.iter()).all(|(a, b)| a <= b)
}

/// Normal form of a monomial modulo the Gröbner basis, expressed in
/// staircase coordinates (a vector of length `dim`).
fn normal_form_monomial<D: Domain>(
    m: &[usize],
    gb: &GroebnerBasis<D, impl MonomialOrder>,
    staircase: &[Vec<usize>],
    domain: &D,
) -> Vec<D::Element> {
    // Reduce the monomial x^m by the basis; the result is a linear
    // combination of staircase monomials.
    let poly = SparseMultivariatePolynomial::from_terms(
        domain.clone(),
        m.len(),
        vec![(m.to_vec(), domain.one())],
    );
    let nf = poly.reduce(&gb.basis);
    let mut coords = vec![domain.zero(); staircase.len()];
    for (exp, coeff) in nf.terms_ref() {
        if let Some(pos) = staircase
            .iter()
            .position(|s| s.as_slice() == exp.as_slice())
        {
            coords[pos] = coeff.clone();
        }
    }
    coords
}

/// Find coefficients `c_i` with `nf = Σ c_i · seen_i`, or `None` when
/// linearly independent (Gaussian elimination over the field).
fn find_relation<D: Domain>(
    seen: &[Vec<D::Element>],
    nf: &[D::Element],
    domain: &D,
) -> Option<Vec<D::Element>> {
    if seen.is_empty() {
        return None;
    }
    let rows = seen.len();
    let cols = nf.len();
    // Solve [seen | nf] for the augmented column: seen·c = nf.
    let mut mat: Vec<Vec<D::Element>> = (0..cols)
        .map(|r| {
            let mut row: Vec<D::Element> = seen.iter().map(|s| s[r].clone()).collect();
            row.push(nf[r].clone());
            row
        })
        .collect();

    // Forward elimination.
    let mut pivot_cols = Vec::new();
    let mut r = 0;
    for c in 0..rows {
        // Find pivot.
        let mut piv = None;
        for (rr, row) in mat.iter().enumerate().skip(r) {
            if !domain.is_zero(&row[c]) {
                piv = Some(rr);
                break;
            }
        }
        let Some(piv) = piv else { continue };
        mat.swap(r, piv);
        let inv = domain.inv(&mat[r][c].clone())?;
        for elt in mat[r].iter_mut().take(rows + 1).skip(c) {
            *elt = domain.mul(elt, &inv);
        }
        for rr in 0..cols {
            if rr != r && !domain.is_zero(&mat[rr][c]) {
                let factor = mat[rr][c].clone();
                #[allow(clippy::needless_range_loop)]
                for cc in c..=rows {
                    let sub = domain.mul(&factor, &mat[r][cc]);
                    mat[rr][cc] = domain.sub(&mat[rr][cc], &sub);
                }
            }
        }
        pivot_cols.push(c);
        r += 1;
        if r == cols {
            break;
        }
    }

    // Check consistency: any row with all-zero seen-part but nonzero
    // augmented part means no solution (nf is independent).
    for row in &mat {
        let seen_zero = (0..rows).all(|c| domain.is_zero(&row[c]));
        if seen_zero && !domain.is_zero(&row[rows]) {
            return None;
        }
    }

    // Back-substitute.
    let mut c_vec = vec![domain.zero(); rows];
    for (i, &pc) in pivot_cols.iter().enumerate() {
        c_vec[pc] = mat[i][rows].clone();
    }
    Some(c_vec)
}

/// Mark all multiples of `m` (up to a total-degree bound) as visited so
/// they are never inserted into the boundary again.
fn mark_multiples(
    visited: &mut HashMap<Vec<usize>, bool>,
    m: &[usize],
    n_vars: usize,
    max_deg: usize,
) {
    let mut queue = vec![m.to_vec()];
    while let Some(cur) = queue.pop() {
        for var in 0..n_vars {
            let mut next = cur.clone();
            next[var] += 1;
            if next.iter().sum::<usize>() > max_deg {
                continue;
            }
            if visited.insert(next.clone(), true).is_none() {
                queue.push(next);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::{Grevlex, Lex};
    use ocas_domain::{Rational, RationalDomain};

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn fglm_linear_ideal() {
        // ideal: x + y, x - y over ℚ — zero-dimensional, dim 1.
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
        let gb_lex = crate::groebner::f4::f4(&[f1, f2]);
        let gb_grevlex = fglm::<_, Grevlex>(&gb_lex).expect("zero-dimensional");
        assert!(gb_grevlex.is_groebner_basis());
    }

    #[test]
    fn fglm_zero_dim_quadratic() {
        // ideal: x² - 1, y - x over ℚ — zero-dimensional, dim 2.
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![2, 0], r(1, 1)), (vec![0, 0], r(-1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![0, 1], r(1, 1)), (vec![1, 0], r(-1, 1))],
        );
        let gb_lex = crate::groebner::f4::f4(&[f1, f2]);
        let gb_grevlex = fglm::<_, Grevlex>(&gb_lex).expect("zero-dimensional");
        assert!(gb_grevlex.is_groebner_basis());
    }
}
