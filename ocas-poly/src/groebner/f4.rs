//! F4 algorithm for computing Gröbner bases.
//!
//! Implements the matrix-based F4 algorithm from Faugère's 1999 paper.
//!
//! The key idea: replace sequential S-polynomial reductions with batched
//! sparse-matrix row operations (Gaussian elimination over the coefficient field).

use std::collections::HashMap;

use smallvec::SmallVec;

use ocas_domain::{Domain, FiniteField};

use super::GroebnerBasis;
use crate::sparse::{MonomialOrder, SparseMultivariatePolynomial, monomial_divides, monomial_lcm};

/// A critical pair with precomputed lcm metadata for Gebauer-Moeller filtering.
#[derive(Debug, Clone)]
struct CriticalPair {
    idx1: usize,
    idx2: usize,
    lcm: SmallVec<[usize; 4]>,
    /// Total degree of the lcm (primary selection key).
    degree: usize,
}

impl CriticalPair {
    fn new<D: Domain, O: MonomialOrder>(
        basis: &[SparseMultivariatePolynomial<D, O>],
        i: usize,
        j: usize,
    ) -> Option<Self> {
        let lm_i = basis[i].leading_monomial()?;
        let lm_j = basis[j].leading_monomial()?;
        let lcm = monomial_lcm(lm_i, lm_j);
        let degree: usize = lcm.iter().sum();
        Some(Self {
            idx1: i,
            idx2: j,
            lcm,
            degree,
        })
    }
}

/// Cache entry: (exponent_diff, cached_polynomial).
type SimpCache<D, O> = Vec<(SmallVec<[usize; 4]>, SparseMultivariatePolynomial<D, O>)>;

/// Tracks a monomial's state during symbolic preprocessing.
#[derive(Debug, Clone)]
struct MonomialData {
    /// Whether this monomial has been processed (reducer found or confirmed absent).
    present: bool,
    /// Column index in the matrix (assigned during column construction).
    column: usize,
}

// =========================================================================
//  Public entry point
// =========================================================================

/// Compute a Gröbner basis using the F4 algorithm.
///
/// Operates over any field domain. For `FiniteField`, a specialized ℤ_p
/// fast path is used for row echelon form computation.
///
/// # Example
///
/// ```
/// use ocas_domain::{RationalDomain, Rational};
/// use ocas_poly::sparse::Lex;
/// use ocas_poly::groebner::f4::f4;
/// use ocas_poly::SparseMultivariatePolynomial;
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
/// let gb = f4(&[f1, f2]);
/// assert!(!gb.basis.is_empty());
/// ```
pub fn f4<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O> {
    if ideal.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Filter zeros and make monic.
    let mut basis: Vec<SparseMultivariatePolynomial<D, O>> =
        ideal.iter().filter(|p| !p.is_zero()).cloned().collect();
    for p in &mut basis {
        make_monic(p);
    }
    if basis.is_empty() {
        return GroebnerBasis { basis };
    }

    // Build initial critical pairs with Gebauer-Moeller filtering.
    let mut pairs: Vec<CriticalPair> = Vec::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            if let Some(cp) = CriticalPair::new(&basis, i, j) {
                // First criterion: skip coprime pairs.
                let lm_i = basis[i].leading_monomial().unwrap();
                let lm_j = basis[j].leading_monomial().unwrap();
                let is_coprime = lm_i
                    .iter()
                    .zip(lm_j.iter())
                    .all(|(a, b)| *a == 0 || *b == 0);
                if !is_coprime {
                    pairs.push(cp);
                }
            }
        }
    }

    // Check if we can use ℤ_p fast path.
    let use_fp = std::any::TypeId::of::<D>() == std::any::TypeId::of::<FiniteField>();
    let prime_i64 = if use_fp {
        let domain_ptr = basis[0].domain() as *const D;
        let ff = unsafe { &*domain_ptr.cast::<FiniteField>() };
        ff.prime_u64() as i64
    } else {
        0
    };

    // Reusable buffers.
    // MonomialData tracks whether a monomial has been processed during
    // symbolic preprocessing (present=true means "already handled").
    let mut all_monomials: HashMap<SmallVec<[usize; 4]>, MonomialData> = HashMap::new();
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut matrix: Vec<Vec<(D::Element, usize)>> = Vec::new();
    let mut fp_matrix: Vec<Vec<(i64, usize)>> = Vec::new();
    let mut pivots: Vec<Option<usize>> = Vec::new();
    let mut fp_buffer: Vec<i64> = Vec::new();

    // Simplification cache: for each basis element, a list of
    // (exponent_diff, cached_poly) to avoid recomputing products.
    let mut simplifications: Vec<SimpCache<D, O>> = basis
        .iter()
        .map(|p| {
            let zero_exp = SmallVec::from_elem(0, p.n_vars());
            vec![(zero_exp, p.clone())]
        })
        .collect();

    while !pairs.is_empty() {
        // --- Selection: find minimum lcm degree ---
        let min_deg = pairs.iter().map(|cp| cp.degree).min().unwrap();

        // Collect all pairs with minimum degree.
        let selected: Vec<CriticalPair> = pairs
            .iter()
            .filter(|cp| cp.degree == min_deg)
            .cloned()
            .collect();

        // Remove selected pairs from the pool.
        let sel_set: std::collections::HashSet<(usize, usize)> =
            selected.iter().map(|cp| (cp.idx1, cp.idx2)).collect();
        pairs.retain(|cp| !sel_set.contains(&(cp.idx1, cp.idx2)));

        if selected.is_empty() {
            continue;
        }

        // --- Build matrix rows from selected pairs ---
        all_monomials.clear();
        monomial_list.clear();
        matrix.clear();

        for cp in &selected {
            let i = cp.idx1;
            let j = cp.idx2;
            let lm_i = basis[i].leading_monomial().unwrap();
            let lm_j = basis[j].leading_monomial().unwrap();
            let lcm_exp = &cp.lcm;

            // Multiply f_i by x^(lcm/lm_i) and f_j by x^(lcm/lm_j).
            let diff_i: SmallVec<[usize; 4]> = lcm_exp
                .iter()
                .zip(lm_i.iter())
                .map(|(&a, b)| a - b)
                .collect();
            let diff_j: SmallVec<[usize; 4]> = lcm_exp
                .iter()
                .zip(lm_j.iter())
                .map(|(&a, b)| a - b)
                .collect();

            let fi_mult = basis[i].mul_monomial(&diff_i);
            let fj_mult = basis[j].mul_monomial(&diff_j);

            // S-polynomial (monic basis ⇒ lc = 1).
            let s_poly = fi_mult.sub(&fj_mult);

            if !s_poly.is_zero() {
                add_poly_to_matrix(&s_poly, &mut matrix, &mut all_monomials, &mut monomial_list);
            }
        }

        if matrix.is_empty() {
            continue;
        }

        // --- Iterative symbolic preprocessing (Faugère F4 key innovation) ---
        // Scan each matrix row's monomials. For each unseen monomial, search
        // for a reducer in the basis and append it as a new row. Repeat until
        // no new rows are added. This ensures the matrix contains all
        // necessary reduction information for the echelon step.
        //
        // Reference: Symbolica groebner.rs L262-288.
        let mut i = 0;
        while i < matrix.len() {
            // Collect monomials from this row by scanning the column indices.
            // We need the actual exponent vectors, so we look them up in monomial_list.
            // But monomial_list isn't populated yet at this point — we use all_monomials.
            // Instead, scan the polynomial that produced this row.
            // Since we store rows as (coeff, col_idx), we need a different approach:
            // iterate all_monomials entries that are not yet marked present.

            // Mark all monomials in current matrix rows as present.
            // For the first pass, the S-polynomial rows already registered their monomials
            // in add_poly_to_matrix. We need to find reducers for monomials that are NOT
            // leading monomials of the matrix rows.

            // Collect all monomials from the matrix that haven't been processed yet.
            let mut new_monomials: Vec<SmallVec<[usize; 4]>> = Vec::new();
            for (exp, md) in all_monomials.iter() {
                if !md.present {
                    new_monomials.push(exp.clone());
                }
            }

            if new_monomials.is_empty() {
                break;
            }

            // Mark them as present.
            for exp in &new_monomials {
                if let Some(md) = all_monomials.get_mut(exp) {
                    md.present = true;
                }
            }

            // For each new monomial, find a reducer in the basis.
            for exp in &new_monomials {
                let mut best: Option<usize> = None;
                for (bi, bp) in basis.iter().enumerate() {
                    if let Some(blm) = bp.leading_monomial()
                        && monomial_divides(exp, blm)
                    {
                        match best {
                            Some(b) if basis[b].n_terms() <= bp.n_terms() => {}
                            _ => best = Some(bi),
                        }
                    }
                }
                if let Some(bi) = best {
                    let blm = basis[bi].leading_monomial().unwrap();
                    let diff: SmallVec<[usize; 4]> =
                        exp.iter().zip(blm.iter()).map(|(a, b)| a - b).collect();
                    let reducer = get_simplified(&simplifications[bi], &diff, &basis[bi]);
                    // Add reducer row to matrix, registering any new monomials.
                    add_poly_to_matrix(
                        &reducer,
                        &mut matrix,
                        &mut all_monomials,
                        &mut monomial_list,
                    );
                }
            }

            i += 1;
        }

        if matrix.is_empty() || monomial_list.is_empty() {
            continue;
        }

        let ncols = monomial_list.len();

        // --- Sort columns: descending monomial order ---
        let mut col_order: Vec<usize> = (0..ncols).collect();
        col_order.sort_unstable_by(|&a, &b| O::cmp(&monomial_list[a], &monomial_list[b]));

        // Build inverse column map: old_col → new_col.
        let mut col_inv = vec![0usize; ncols];
        for (new_col, &old_col) in col_order.iter().enumerate() {
            col_inv[old_col] = new_col;
        }

        // Remap column indices in matrix.
        for row in &mut matrix {
            for (_, col) in row.iter_mut() {
                *col = col_inv[*col];
            }
        }

        // Build sorted monomial list.
        let mut sorted_monomials: Vec<SmallVec<[usize; 4]>> = vec![SmallVec::new(); ncols];
        for (new_col, &old_col) in col_order.iter().enumerate() {
            sorted_monomials[new_col] = monomial_list[old_col].clone();
        }

        // --- Row echelon form ---
        if use_fp {
            fp_matrix.clear();
            fp_matrix.resize(matrix.len(), vec![]);
            for (row_idx, row) in matrix.iter().enumerate() {
                for (coeff, col) in row {
                    let c = domain_to_i64_fp::<D>(coeff, prime_i64);
                    if c != 0 {
                        fp_matrix[row_idx].push((c, *col));
                    }
                }
            }
            echelonize_fp(
                &mut fp_matrix,
                ncols,
                prime_i64,
                &mut pivots,
                &mut fp_buffer,
            );
            matrix.clear();
            matrix.resize(fp_matrix.len(), vec![]);
            for (row_idx, row) in fp_matrix.iter().enumerate() {
                for &(c, col) in row {
                    matrix[row_idx].push((i64_to_domain_fp(basis[0].domain(), c, prime_i64), col));
                }
            }
        } else {
            echelonize_generic(&mut matrix, ncols, basis[0].domain(), &mut pivots);
        }

        // --- Extract new polynomials from reduced rows ---
        for row in &matrix {
            if row.is_empty() {
                continue;
            }
            let mut poly = basis[0].zero();
            for (coeff, col) in row.iter().rev() {
                poly.append_monomial(coeff.clone(), &sorted_monomials[*col]);
            }
            make_monic(&mut poly);
            if poly.is_zero() {
                continue;
            }

            // Reduce by existing basis to get a minimal representative.
            let reduced = poly.reduce(&basis);
            if !reduced.is_zero() {
                poly = reduced;
                make_monic(&mut poly);
            } else {
                continue;
            }

            let new_lm = poly.leading_monomial().unwrap().clone();

            // Skip if a polynomial with this leading monomial already exists.
            if basis.iter().any(|bp| {
                bp.leading_monomial()
                    .is_some_and(|blm| blm.as_slice() == new_lm.as_slice())
            }) {
                continue;
            }

            // Add to basis with Gebauer-Moeller pair filtering.
            update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
        }
    }

    // Post-processing: minimize and inter-reduce.
    GroebnerBasis { basis }.minimize().auto_reduce()
}

// =========================================================================
//  Gebauer-Moeller pair management
// =========================================================================

/// Add a new polynomial to the basis and update critical pairs using
/// the Gebauer-Moeller criteria.
///
/// Reference: Symbolica groebner.rs L475-545; Becker-Weispfenning
/// "A Computational Approach to Commutative Algebra".
fn update_pairs<D: Domain, O: MonomialOrder>(
    basis: &mut Vec<SparseMultivariatePolynomial<D, O>>,
    pairs: &mut Vec<CriticalPair>,
    simplifications: &mut Vec<SimpCache<D, O>>,
    new_poly: SparseMultivariatePolynomial<D, O>,
) {
    let new_lm = match new_poly.leading_monomial() {
        Some(m) => m.clone(),
        None => {
            basis.push(new_poly);
            return;
        }
    };
    let new_idx = basis.len();
    basis.push(new_poly);
    // Initialize simplification cache for the new basis element.
    simplifications.push(vec![(
        SmallVec::from_elem(0, basis[new_idx].n_vars()),
        basis[new_idx].clone(),
    )]);

    // Generate new pairs with existing basis elements.
    let mut new_pairs: Vec<(CriticalPair, bool)> = Vec::new();
    for i in 0..new_idx {
        if let Some(cp) = CriticalPair::new(basis, i, new_idx) {
            // First criterion: skip coprime pairs (disjoint LMs).
            let lm_b = basis[i].leading_monomial().unwrap();
            let is_coprime = lm_b
                .iter()
                .zip(new_lm.iter())
                .all(|(a, b)| *a == 0 || *b == 0);
            if !is_coprime {
                new_pairs.push((cp, true));
            }
        }
    }

    // Second criterion (Gebauer-Moeller): among new pairs, keep only those
    // whose lcm is minimal — no other new pair has a strictly smaller lcm.
    for i in 0..new_pairs.len() {
        if !new_pairs[i].1 {
            continue;
        }
        let dominated = new_pairs.iter().enumerate().any(|(j, (pj, kj))| {
            if !*kj || i == j {
                return false;
            }
            // Check if lcm[j] strictly divides lcm[i].
            new_pairs[i]
                .0
                .lcm
                .iter()
                .zip(pj.lcm.iter())
                .all(|(a, b)| a >= b)
                && new_pairs[i]
                    .0
                    .lcm
                    .iter()
                    .zip(pj.lcm.iter())
                    .any(|(a, b)| a > b)
        });
        if dominated {
            new_pairs[i].1 = false;
        }
    }

    // Clean existing pairs: remove pairs made redundant by the new polynomial.
    // A pair is redundant if the new polynomial's LM divides the pair's lcm
    // in a way that makes the old pair no longer needed.
    pairs.retain(|cp| {
        // Keep pair if new_lm does NOT strictly divide its lcm.
        let dominated = cp.lcm.iter().zip(new_lm.iter()).all(|(a, b)| a >= b)
            && cp.lcm.iter().zip(new_lm.iter()).any(|(a, b)| a > b);
        !dominated
    });

    // Add non-redundant new pairs.
    for (cp, keep) in new_pairs {
        if keep {
            pairs.push(cp);
        }
    }

    // NOTE: Do NOT remove basis elements here — it invalidates pair indices.
    // The minimize() post-processing step handles LM-divisible removal.
}

// =========================================================================
//  Helper: add a polynomial as a matrix row
// =========================================================================

fn add_poly_to_matrix<D: Domain, O: MonomialOrder>(
    poly: &SparseMultivariatePolynomial<D, O>,
    matrix: &mut Vec<Vec<(D::Element, usize)>>,
    monomial_map: &mut HashMap<SmallVec<[usize; 4]>, MonomialData>,
    monomial_list: &mut Vec<SmallVec<[usize; 4]>>,
) {
    let mut row: Vec<(D::Element, usize)> = Vec::new();
    for (exp, coeff) in poly.sorted_terms().iter().rev() {
        if poly.domain().is_zero(coeff) {
            continue;
        }
        monomial_map.entry((*exp).clone()).or_insert_with(|| {
            let idx = monomial_list.len();
            monomial_list.push((*exp).clone());
            MonomialData {
                present: false,
                column: idx,
            }
        });
        let md = monomial_map.get(*exp).unwrap();
        row.push(((*coeff).clone(), md.column));
    }
    if !row.is_empty() {
        matrix.push(row);
    }
}

// =========================================================================
//  Simplification cache
// =========================================================================

/// Look up a cached polynomial for the given exponent diff, or compute
/// `basis_poly * x^diff` and cache it.
///
/// Reference: Symbolica groebner.rs L167-185.
fn get_simplified<D: Domain, O: MonomialOrder>(
    cache: &SimpCache<D, O>,
    diff: &[usize],
    basis_poly: &SparseMultivariatePolynomial<D, O>,
) -> SparseMultivariatePolynomial<D, O> {
    // Check exact match first.
    for (cached_diff, cached_poly) in cache.iter().rev() {
        if cached_diff.as_slice() == diff {
            return cached_poly.clone();
        }
    }
    // Check if any cached diff divides the requested diff.
    for (cached_diff, cached_poly) in cache.iter().rev() {
        if diff.iter().zip(cached_diff.iter()).all(|(d, c)| d >= c) {
            let remaining: SmallVec<[usize; 4]> = diff
                .iter()
                .zip(cached_diff.iter())
                .map(|(d, c)| d - c)
                .collect();
            return cached_poly.mul_monomial(&remaining);
        }
    }
    // Fallback: compute directly.
    basis_poly.mul_monomial(diff)
}

// =========================================================================
//  Row echelon form: ℤ_p fast path
// =========================================================================

#[allow(clippy::needless_range_loop)]
fn echelonize_fp(
    matrix: &mut Vec<Vec<(i64, usize)>>,
    ncols: usize,
    prime: i64,
    pivots: &mut Vec<Option<usize>>,
    buffer: &mut Vec<i64>,
) {
    let p = prime;
    let p2 = p * p;

    sort_rows(matrix);

    pivots.clear();
    pivots.resize(ncols, None);

    // Identify initial pivots.
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }
        let col = matrix[r][0].1;
        if pivots[col].is_none() {
            pivots[col] = Some(r);
            if matrix[r][0].0 != 1 {
                let inv = mod_inv(matrix[r][0].0, p);
                for (c, _) in &mut matrix[r] {
                    *c = (*c * inv) % p;
                }
            }
        }
    }

    // Reduce rows.
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }

        let first_col = matrix[r][0].1;
        if pivots[first_col].is_none() {
            pivots[first_col] = Some(r);
            if matrix[r][0].0 != 1 {
                let inv = mod_inv(matrix[r][0].0, p);
                for (c, _) in &mut matrix[r] {
                    *c = (*c * inv) % p;
                }
            }
        }

        if pivots[first_col] == Some(r) {
            continue;
        }

        // Dense buffer elimination.
        buffer.clear();
        buffer.resize(ncols, 0);
        for &(c, col) in &matrix[r] {
            buffer[col] = c;
        }

        for i in 0..ncols {
            buffer[i] %= p;
            if buffer[i] == 0 {
                continue;
            }
            let pi = match pivots[i] {
                Some(pi) => pi,
                None => {
                    // New pivot.
                    pivots[i] = Some(r);
                    let inv = mod_inv(buffer[i], p);
                    buffer[i] = 1;
                    for j in (i + 1)..ncols {
                        buffer[j] = (buffer[j] * inv) % p;
                    }
                    matrix[r].clear();
                    for (col, val) in buffer.iter_mut().enumerate() {
                        let v = *val % p;
                        if v != 0 {
                            matrix[r].push((v, col));
                            *val = 0;
                        }
                    }
                    continue;
                }
            };

            let c = buffer[i];
            buffer[i] = 0;
            for &(pc, pcol) in &matrix[pi] {
                if pcol <= i {
                    continue;
                }
                let m = pc * c;
                let t = buffer[pcol];
                buffer[pcol] = if t >= m { t - m } else { t + p2 - m };
            }
        }

        // Write back if not already written as new pivot.
        if matrix[r].is_empty() || matrix[r][0].1 != first_col {
            matrix[r].clear();
            for (col, val) in buffer.iter_mut().enumerate() {
                let v = *val % p;
                if v != 0 {
                    matrix[r].push((v, col));
                    *val = 0;
                }
            }
        }
    }

    matrix.retain(|r| !r.is_empty());
}

// =========================================================================
//  Row echelon form: generic domain path
// =========================================================================

#[allow(clippy::needless_range_loop, clippy::collapsible_if)]
fn echelonize_generic<D: Domain>(
    matrix: &mut Vec<Vec<(D::Element, usize)>>,
    ncols: usize,
    domain: &D,
    pivots: &mut Vec<Option<usize>>,
) {
    sort_rows(matrix);

    pivots.clear();
    pivots.resize(ncols, None);

    let mut buffer: Vec<D::Element> = vec![domain.zero(); ncols];

    // Identify initial pivots.
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }
        let col = matrix[r][0].1;
        if pivots[col].is_none() {
            pivots[col] = Some(r);
            let lc = matrix[r][0].0.clone();
            if !domain.is_one(&lc) {
                if let Some(inv) = domain.inv(&lc) {
                    for (c, _) in &mut matrix[r] {
                        *c = domain.mul(c, &inv);
                    }
                }
            }
        }
    }

    // Reduce rows.
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }

        let first_col = matrix[r][0].1;
        if pivots[first_col].is_none() {
            pivots[first_col] = Some(r);
            let lc = matrix[r][0].0.clone();
            if !domain.is_one(&lc) {
                if let Some(inv) = domain.inv(&lc) {
                    for (c, _) in &mut matrix[r] {
                        *c = domain.mul(c, &inv);
                    }
                }
            }
            continue;
        }

        if pivots[first_col] == Some(r) {
            continue;
        }

        // Copy to dense buffer.
        for b in buffer.iter_mut() {
            *b = domain.zero();
        }
        for (c, col) in &matrix[r] {
            buffer[*col] = c.clone();
        }

        // Eliminate.
        for i in 0..ncols {
            if domain.is_zero(&buffer[i]) {
                continue;
            }
            let pi = match pivots[i] {
                Some(pi) => pi,
                None => {
                    pivots[i] = Some(r);
                    if let Some(inv) = domain.inv(&buffer[i]) {
                        buffer[i] = domain.one();
                        for j in (i + 1)..ncols {
                            buffer[j] = domain.mul(&buffer[j], &inv);
                        }
                    }
                    matrix[r].clear();
                    for (col, val) in buffer.iter_mut().enumerate() {
                        if !domain.is_zero(val) {
                            matrix[r].push((val.clone(), col));
                            *val = domain.zero();
                        }
                    }
                    continue;
                }
            };

            let c = buffer[i].clone();
            buffer[i] = domain.zero();
            for (pc, pcol) in &matrix[pi] {
                if *pcol <= i {
                    continue;
                }
                let product = domain.mul(pc, &c);
                buffer[*pcol] = domain.sub(&buffer[*pcol], &product);
            }
        }

        // Write back if not already written as new pivot.
        if matrix[r].is_empty() || matrix[r][0].1 != first_col {
            matrix[r].clear();
            for (col, val) in buffer.iter_mut().enumerate() {
                if !domain.is_zero(val) {
                    matrix[r].push((val.clone(), col));
                    *val = domain.zero();
                }
            }
        }
    }

    matrix.retain(|r| !r.is_empty());
}

// =========================================================================
//  Utilities
// =========================================================================

fn sort_rows<T>(matrix: &mut [Vec<(T, usize)>]) {
    matrix.sort_unstable_by(|a, b| match (a.first(), b.first()) {
        (Some((_, ca)), Some((_, cb))) => ca.cmp(cb).then(a.len().cmp(&b.len())),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
}

fn mod_inv(a: i64, p: i64) -> i64 {
    let a = ((a % p) + p) % p;
    if a == 0 {
        return 0;
    }
    let (mut old_r, mut r) = (a, p);
    let (mut old_s, mut s) = (1i64, 0i64);
    while r != 0 {
        let q = old_r / r;
        let tmp = r;
        r = old_r - q * r;
        old_r = tmp;
        let tmp = s;
        s = old_s - q * s;
        old_s = tmp;
    }
    ((old_s % p) + p) % p
}

#[allow(clippy::collapsible_if)]
fn make_monic<D: Domain, O: MonomialOrder>(p: &mut SparseMultivariatePolynomial<D, O>) {
    if p.is_zero() {
        return;
    }
    if let Some(lc) = p.leading_coeff().cloned()
        && let Some(inv_lc) = p.domain().inv(&lc)
    {
        let terms: Vec<(Vec<usize>, D::Element)> = p
            .terms_ref()
            .iter()
            .map(|(exp, coeff)| (exp.to_vec(), p.domain().mul(coeff, &inv_lc)))
            .collect();
        *p = SparseMultivariatePolynomial::from_terms(p.domain().clone(), p.n_vars(), terms);
    }
}

fn domain_to_i64_fp<D: Domain + 'static>(elem: &D::Element, prime: i64) -> i64 {
    if std::any::TypeId::of::<D>() == std::any::TypeId::of::<FiniteField>() {
        let ff_elem =
            unsafe { &*(elem as *const D::Element as *const <FiniteField as Domain>::Element) };
        let val = ff_elem.value();
        let (_, digits) = val.to_u64_digits();
        if digits.is_empty() {
            0
        } else {
            (digits[0] as i64) % prime
        }
    } else {
        0
    }
}

fn i64_to_domain_fp<D: Domain + 'static>(domain: &D, val: i64, prime: i64) -> D::Element {
    if std::any::TypeId::of::<D>() == std::any::TypeId::of::<FiniteField>() {
        let ff_domain = unsafe { &*(domain as *const D as *const FiniteField) };
        let v = ((val % prime) + prime) % prime;
        let elem = ff_domain.element(num_bigint::BigInt::from(v));
        unsafe {
            (&*(&elem as *const <FiniteField as Domain>::Element as *const D::Element)).clone()
        }
    } else {
        domain.zero()
    }
}

// =========================================================================
//  Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Lex;
    use num_bigint::BigInt;
    use ocas_domain::{FiniteField, Rational, RationalDomain};

    fn rat(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn f4_empty_ideal() {
        let gb = f4::<RationalDomain, Lex>(&[]);
        assert!(gb.basis.is_empty());
    }

    #[test]
    fn f4_single_polynomial() {
        let f = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            2,
            vec![(vec![2, 0], rat(1, 1)), (vec![0, 0], rat(-1, 1))],
        );
        let gb = f4(&[f]);
        assert_eq!(gb.basis.len(), 1);
    }

    #[test]
    fn f4_linear_system() {
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            2,
            vec![(vec![1, 0], rat(1, 1)), (vec![0, 1], rat(1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            2,
            vec![(vec![1, 0], rat(1, 1)), (vec![0, 1], rat(-1, 1))],
        );
        let gb = f4(&[f1, f2]);
        assert!(gb.basis.len() >= 2, "expected >= 2, got {}", gb.basis.len());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f4_two_variable_ideal() {
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            2,
            vec![(vec![2, 0], rat(1, 1)), (vec![0, 1], rat(-1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            RationalDomain,
            2,
            vec![(vec![3, 0], rat(1, 1)), (vec![1, 0], rat(-1, 1))],
        );
        let gb = f4(&[f1, f2]);
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f4_cyclic_3_zp() {
        // cyclic-3 over ℚ — requires complete F4 symbolic preprocessing
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![
                (vec![1, 0, 0], rat(1, 1)),
                (vec![0, 1, 0], rat(1, 1)),
                (vec![0, 0, 1], rat(1, 1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![
                (vec![1, 1, 0], rat(1, 1)),
                (vec![0, 1, 1], rat(1, 1)),
                (vec![1, 0, 1], rat(1, 1)),
            ],
        );
        let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            3,
            vec![(vec![1, 1, 1], rat(1, 1)), (vec![0, 0, 0], rat(-1, 1))],
        );
        let gb = f4(&[f1, f2, f3]);
        assert!(!gb.basis.is_empty());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f4_cyclic_4_zp() {
        let field = FiniteField::new(BigInt::from(13u32));
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            4,
            vec![
                (vec![1, 0, 0, 0], field.element(1)),
                (vec![0, 1, 0, 0], field.element(1)),
                (vec![0, 0, 1, 0], field.element(1)),
                (vec![0, 0, 0, 1], field.element(1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            4,
            vec![
                (vec![1, 1, 0, 0], field.element(1)),
                (vec![0, 1, 1, 0], field.element(1)),
                (vec![0, 0, 1, 1], field.element(1)),
                (vec![1, 0, 0, 1], field.element(1)),
            ],
        );
        let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            4,
            vec![
                (vec![1, 1, 1, 0], field.element(1)),
                (vec![0, 1, 1, 1], field.element(1)),
                (vec![1, 0, 1, 1], field.element(1)),
                (vec![1, 1, 0, 1], field.element(1)),
            ],
        );
        let f4_poly = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            4,
            vec![
                (vec![1, 1, 1, 1], field.element(1)),
                (vec![0, 0, 0, 0], field.element(12)),
            ],
        );
        let gb = f4(&[f1, f2, f3, f4_poly]);
        assert!(!gb.basis.is_empty());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    #[ignore = "timing test: ~55s per run"]
    fn f4_cyclic_5_fp13_timing() {
        let field = FiniteField::new(BigInt::from(13u32));
        let n = 5;
        let mut gens = Vec::new();
        for k in 1..n {
            let mut terms = Vec::new();
            for start in 0..n {
                let mut exps = vec![0usize; n];
                for j in 0..k {
                    exps[(start + j) % n] = 1;
                }
                terms.push((exps, field.element(1)));
            }
            gens.push(SparseMultivariatePolynomial::<_, Lex>::from_terms(
                field.clone(),
                n,
                terms,
            ));
        }
        let full_exps = vec![1usize; n];
        gens.push(SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            n,
            vec![
                (full_exps, field.element(1)),
                (vec![0usize; n], field.element(12)),
            ],
        ));
        let start = std::time::Instant::now();
        let gb = f4(&gens);
        let elapsed = start.elapsed();
        eprintln!("cyclic-5 Fp13: {:.2?}, basis={}", elapsed, gb.basis.len());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn mod_inv_basic() {
        assert_eq!(mod_inv(3, 7), 5);
        assert_eq!(mod_inv(2, 7), 4);
        assert_eq!(mod_inv(1, 13), 1);
    }

    #[test]
    fn grlex_ordering() {
        use crate::sparse::Grlex;
        assert_eq!(Grlex::cmp(&[2, 0], &[1, 1]), std::cmp::Ordering::Greater);
        assert_eq!(Grlex::cmp(&[1, 1], &[0, 2]), std::cmp::Ordering::Greater);
        assert_eq!(Grlex::cmp(&[0, 2], &[1, 0]), std::cmp::Ordering::Less);
        assert_eq!(Grlex::cmp(&[1, 0], &[0, 2]), std::cmp::Ordering::Greater);
    }
}
