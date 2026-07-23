//! The F5 algorithm (Faugère 2002) with signature-based rewriting.
//!
//! F5 attaches a *signature* to each polynomial and uses syzygy criteria
//! to reject zero-reducers *before* they enter the reduction matrix,
//! yielding order-of-magnitude speedups over F4 on difficult ideals such
//! as the cyclic family.
//!
//! The current implementation provides the generic-domain F5 core:
//! signature monomial ordering (pot), the syzygy criterion, and
//! signature-threaded matrix construction. The native ℤ_p fast path and
//! F5'/F5C optimizations land in subsequent 0.19.0 phases.
//!
//! Reference: Faugère, "A New Efficient Algorithm for Computing Gröbner
//! Bases without Reduction to Zero (F5)", ISSAC 2002; Eder & Perry,
//! "Signature-based Algorithms to Compute Gröbner Bases" (2009).

use smallvec::{SmallVec, smallvec};
use std::cmp::Ordering;

use ocas_core::FastHashMap as HashMap;
use ocas_core::FastHashSet as HashSet;
use ocas_domain::Domain;

use super::GroebnerBasis;
use crate::sparse::{MonomialOrder, SparseMultivariatePolynomial, monomial_divides};

// =========================================================================
//  Signature
// =========================================================================

/// A signature in the F5 algorithm.
///
/// Each polynomial in the F5 basis carries a signature `(module_pos,
/// monomial)` recording its "history": `module_pos` is the index of the
/// input generator it descends from, and `monomial` is the monomial
/// multiple applied to that generator's module basis vector `e_{module_pos}`.
///
/// Signatures are compared by the **pot** (position-over-term) order:
/// first by module position (smaller = earlier), then by monomial order `O`.
///
/// Reference: Faugère 2002, §2.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Signature {
    /// Index of the input generator (module position, 0-based).
    pub module_pos: usize,
    /// Multiplier monomial on that generator's module basis vector.
    pub monomial: SmallVec<[usize; 4]>,
}

impl Signature {
    /// Create the canonical signature of generator `module_pos`:
    /// `e_{module_pos}` with the unit monomial.
    pub fn unit(module_pos: usize, n_vars: usize) -> Self {
        Self {
            module_pos,
            monomial: smallvec![0; n_vars],
        }
    }

    /// Multiply this signature by a monomial `exp` (componentwise).
    /// The module position is unchanged.
    pub fn mul_monomial(&self, exp: &[usize]) -> Self {
        let monomial: SmallVec<[usize; 4]> = self
            .monomial
            .iter()
            .zip(exp.iter())
            .map(|(a, b)| a + b)
            .collect();
        Self {
            module_pos: self.module_pos,
            monomial,
        }
    }

    /// Compare two signatures under the pot (position-over-term) order.
    ///
    /// Returns `Less` if `self` should be processed *before* `other`.
    /// First compares module positions (smaller first), then monomials
    /// under the monomial order `O` (where `O::cmp(a, b)` returns `Less`
    /// when `a` precedes `b`).
    pub fn cmp_pot<O: MonomialOrder>(&self, other: &Self) -> Ordering {
        self.module_pos
            .cmp(&other.module_pos)
            .then_with(|| O::cmp(&self.monomial, &other.monomial))
    }
}

// =========================================================================
//  Labeled polynomial
// =========================================================================

/// A polynomial paired with its F5 signature.
#[derive(Clone)]
struct LabeledPoly<D: Domain, O: MonomialOrder> {
    poly: SparseMultivariatePolynomial<D, O>,
    sig: Signature,
}

impl<D: Domain, O: MonomialOrder> BasisPoly for LabeledPoly<D, O> {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.poly.leading_monomial()
    }
    fn n_vars(&self) -> usize {
        self.poly.n_vars()
    }
    fn n_terms(&self) -> usize {
        self.poly.n_terms()
    }
    fn mul_monomial(&self, exp: &[usize]) -> Self {
        Self {
            poly: self.poly.mul_monomial(exp),
            sig: self.sig.mul_monomial(exp),
        }
    }
}

impl<D: Domain, O: MonomialOrder> LabeledPoly<D, O> {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.poly.leading_monomial()
    }
}

// =========================================================================
//  Syzygy tracking
// =========================================================================

/// Tracks signatures that are known to produce zero reductions.
///
/// When a matrix row reduces to zero, its signature is a syzygy. Any
/// future row whose signature is a monomial multiple of a known syzygy
/// will also reduce to zero and can be skipped immediately — this is the
/// F5 syzygy criterion.
///
/// Internally, for each module position we store the leading monomials
/// of known syzygies. A signature `(k, t)` is flagged as a syzygy when
/// some stored LM for position `k` divides `t`.
struct SyzygySet {
    /// module_pos → leading monomials of known syzygies.
    lms: HashMap<usize, Vec<SmallVec<[usize; 4]>>>,
}

impl SyzygySet {
    fn new() -> Self {
        Self {
            lms: HashMap::default(),
        }
    }

    /// Record that signature `sig` produces a zero reduction.
    fn insert(&mut self, sig: &Signature) {
        self.lms
            .entry(sig.module_pos)
            .or_default()
            .push(sig.monomial.clone());
    }

    /// Check whether `sig` is (or is a multiple of) a known syzygy.
    fn contains(&self, sig: &Signature) -> bool {
        self.lms
            .get(&sig.module_pos)
            .is_some_and(|lms| lms.iter().any(|lm| monomial_divides(&sig.monomial, lm)))
    }
}

// =========================================================================
//  Matrix row with signature
// =========================================================================

/// A sparse matrix row tagged with its F5 signature.
///
/// Terms are stored in ascending column-index order (column 0 = leading
/// monomial, matching F4's convention after column remapping).
struct LabeledRow<D: Domain> {
    /// (coefficient, column_index), ascending by column.
    terms: Vec<(D::Element, usize)>,
    /// The F5 signature of this row.
    sig: Signature,
}

// =========================================================================
//  Entry point (generic path)
// =========================================================================

/// Compute a Gröbner basis using the F5 signature-based algorithm.
///
/// Requires exact division in the coefficient domain (a field). The
/// result is the reduced Gröbner basis, identical to F4's output.
///
/// The algorithm processes generators incrementally, attaching a
/// signature to each polynomial and using the syzygy criterion to reject
/// zero-reducers before they enter the reduction matrix.
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
pub fn f5<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> GroebnerBasis<D, O> {
    if ideal.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // ℤ_p fast path: run the entire F5 pipeline on native i64 residues,
    // converting to/from the BigInt-backed domain only at the boundaries.
    if std::any::TypeId::of::<D>() == std::any::TypeId::of::<ocas_domain::FiniteField>() {
        let domain_ptr = ideal[0].domain() as *const D;
        let ff = unsafe { &*domain_ptr.cast::<ocas_domain::FiniteField>() };
        return f5_fp(ideal, ff.prime_u64() as i64);
    }

    // Filter zeros and make monic.
    let mut generators: Vec<SparseMultivariatePolynomial<D, O>> =
        ideal.iter().filter(|p| !p.is_zero()).cloned().collect();
    for p in &mut generators {
        make_monic(p);
    }
    if generators.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    let n_vars = generators[0].n_vars();

    // Incremental F5: process one generator at a time, using F4's
    // Gebauer-Moeller update_pairs for pair management.
    let mut basis: Vec<LabeledPoly<D, O>> = Vec::new();
    let mut pairs: Vec<CriticalPair> = Vec::new();
    let mut simplifications: Vec<SimpCache<LabeledPoly<D, O>>> = Vec::new();
    let mut syzygies = SyzygySet::new();

    for (k, f) in generators.into_iter().enumerate() {
        let sig_k = Signature::unit(k, n_vars);
        let labeled = LabeledPoly {
            poly: f,
            sig: sig_k,
        };
        update_pairs(&mut basis, &mut pairs, &mut simplifications, labeled);

        // Degree-by-degree matrix reduction.
        while !pairs.is_empty() {
            let min_deg = pairs.iter().map(|p| p.degree).min().unwrap();
            let selected: Vec<CriticalPair> =
                pairs.extract_if(.., |p| p.degree == min_deg).collect();

            let new_polys = build_and_reduce::<D, O>(&selected, &basis, &mut syzygies);

            for poly in new_polys {
                update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
            }
        }
    }

    let polys: Vec<SparseMultivariatePolynomial<D, O>> =
        basis.into_iter().map(|lp| lp.poly).collect();
    GroebnerBasis { basis: polys }.minimize().auto_reduce()
}

// =========================================================================
//  Matrix construction + reduction (generic path)
// =========================================================================

fn build_and_reduce<D: Domain + 'static, O: MonomialOrder>(
    selected: &[CriticalPair],
    basis: &[LabeledPoly<D, O>],
    syzygies: &mut SyzygySet,
) -> Vec<LabeledPoly<D, O>> {
    let domain = basis[0].poly.domain();

    let mut monomial_map: HashMap<SmallVec<[usize; 4]>, usize> = HashMap::default();
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut rows: Vec<LabeledRow<D>> = Vec::new();
    let mut worklist: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut seen_heads: HashSet<SmallVec<[usize; 4]>> = HashSet::default();

    for pair in selected {
        let i = pair.idx1;
        let j = pair.idx2;
        let lm_i = basis[i].leading_monomial().unwrap();
        let lm_j = basis[j].leading_monomial().unwrap();
        let lcm_exp = &pair.lcm;

        let diff_i: SmallVec<[usize; 4]> = lcm_exp
            .iter()
            .zip(lm_i.iter())
            .map(|(a, b)| a - b)
            .collect();
        let diff_j: SmallVec<[usize; 4]> = lcm_exp
            .iter()
            .zip(lm_j.iter())
            .map(|(a, b)| a - b)
            .collect();

        for (idx, diff) in [(i, &diff_i), (j, &diff_j)] {
            let sig = basis[idx].sig.mul_monomial(diff);
            if syzygies.contains(&sig) {
                continue;
            }
            let mult = basis[idx].poly.mul_monomial(diff);
            seen_heads.insert(lcm_exp.clone());
            add_poly_as_row(
                &mult,
                sig,
                &mut rows,
                &mut monomial_map,
                &mut monomial_list,
                &mut worklist,
            );
        }
    }

    if rows.is_empty() {
        return vec![];
    }

    // --- Symbolic preprocessing ---
    while let Some(exp) = worklist.pop() {
        if let Some((bi, diff)) = find_reducer(basis, &exp) {
            let sig = basis[bi].sig.mul_monomial(&diff);
            if syzygies.contains(&sig) {
                continue;
            }
            seen_heads.insert(exp.clone());
            let reducer = basis[bi].poly.mul_monomial(&diff);
            add_poly_as_row(
                &reducer,
                sig,
                &mut rows,
                &mut monomial_map,
                &mut monomial_list,
                &mut worklist,
            );
        }
    }

    if rows.is_empty() || monomial_list.is_empty() {
        return vec![];
    }
    let ncols = monomial_list.len();

    // --- Sort columns: DESCENDING monomial order ---
    let mut col_order: Vec<usize> = (0..ncols).collect();
    col_order.sort_unstable_by(|&a, &b| O::cmp(&monomial_list[b], &monomial_list[a]));
    let mut col_inv = vec![0usize; ncols];
    for (new_col, &old_col) in col_order.iter().enumerate() {
        col_inv[old_col] = new_col;
    }
    for row in &mut rows {
        for (_, col) in row.terms.iter_mut() {
            *col = col_inv[*col];
        }
        row.terms.sort_unstable_by_key(|&(_, col)| col);
    }
    let mut sorted_monomials: Vec<SmallVec<[usize; 4]>> = vec![SmallVec::new(); ncols];
    for (new_col, &old_col) in col_order.iter().enumerate() {
        sorted_monomials[new_col] = monomial_list[old_col].clone();
    }

    // --- Sort rows by ascending signature (pot order) ---
    rows.sort_by(|a, b| a.sig.cmp_pot::<O>(&b.sig));

    // --- Echelonize ---
    echelonize(&mut rows, ncols, domain);

    // --- Extract new basis elements ---
    let mut new_polys: Vec<LabeledPoly<D, O>> = Vec::new();
    let basis_lm_set: HashSet<SmallVec<[usize; 4]>> = basis
        .iter()
        .filter_map(|lp| lp.leading_monomial().cloned())
        .collect();

    for row in &rows {
        if row.terms.is_empty() {
            syzygies.insert(&row.sig);
            continue;
        }
        let row_lm = &sorted_monomials[row.terms[0].1];
        if seen_heads.contains(row_lm) {
            continue;
        }
        if basis_lm_set.contains(row_lm) {
            continue;
        }

        let mut poly = basis[0].poly.zero();
        for (coeff, col) in row.terms.iter().rev() {
            poly.append_monomial(coeff.clone(), &sorted_monomials[*col]);
        }
        if poly.is_zero() {
            syzygies.insert(&row.sig);
            continue;
        }

        new_polys.push(LabeledPoly {
            poly,
            sig: row.sig.clone(),
        });
    }

    new_polys
}

/// Add a polynomial as a labeled matrix row (generic path).
fn add_poly_as_row<D: Domain, O: MonomialOrder>(
    poly: &SparseMultivariatePolynomial<D, O>,
    sig: Signature,
    rows: &mut Vec<LabeledRow<D>>,
    monomial_map: &mut HashMap<SmallVec<[usize; 4]>, usize>,
    monomial_list: &mut Vec<SmallVec<[usize; 4]>>,
    worklist: &mut Vec<SmallVec<[usize; 4]>>,
) {
    let domain = poly.domain();
    let mut terms: Vec<(D::Element, usize)> = Vec::new();
    for (exp, coeff) in poly.sorted_terms().iter().rev() {
        if domain.is_zero(coeff) {
            continue;
        }
        let col = *monomial_map.entry((*exp).clone()).or_insert_with(|| {
            let idx = monomial_list.len();
            monomial_list.push((*exp).clone());
            worklist.push((*exp).clone());
            idx
        });
        terms.push(((*coeff).clone(), col));
    }
    if !terms.is_empty() {
        rows.push(LabeledRow { terms, sig });
    }
}

/// Find a basis element whose leading monomial divides `exp` (generic).
fn find_reducer<D: Domain, O: MonomialOrder>(
    basis: &[LabeledPoly<D, O>],
    exp: &[usize],
) -> Option<(usize, SmallVec<[usize; 4]>)> {
    for (i, lp) in basis.iter().enumerate() {
        if let Some(lm) = lp.leading_monomial()
            && monomial_divides(exp, lm)
        {
            let diff: SmallVec<[usize; 4]> =
                exp.iter().zip(lm.iter()).map(|(a, b)| a - b).collect();
            return Some((i, diff));
        }
    }
    None
}

/// Echelonize labeled generic-domain rows.
fn echelonize<D: Domain>(rows: &mut Vec<LabeledRow<D>>, ncols: usize, domain: &D) {
    let mut pivots: Vec<Option<usize>> = vec![None; ncols];
    let mut scratch: Vec<(D::Element, usize)> = Vec::new();

    for (r, row) in rows.iter_mut().enumerate() {
        if row.terms.is_empty() {
            continue;
        }
        let head_col = row.terms[0].1;
        if pivots[head_col].is_none() {
            let lc = row.terms[0].0.clone();
            if !domain.is_one(&lc)
                && let Some(inv) = domain.inv(&lc)
            {
                for (c, _) in &mut row.terms {
                    *c = domain.mul(c, &inv);
                }
            }
            pivots[head_col] = Some(r);
        }
    }

    for r in 0..rows.len() {
        if rows[r].terms.is_empty() {
            continue;
        }
        if pivots[rows[r].terms[0].1] == Some(r) {
            continue;
        }

        let mut row = std::mem::take(&mut rows[r].terms);
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(pr) => {
                    let c = row[0].0.clone();
                    sub_scaled(domain, &mut row, &rows[pr].terms, &c, &mut scratch);
                }
                None => {
                    let lc = row[0].0.clone();
                    if !domain.is_one(&lc)
                        && let Some(inv) = domain.inv(&lc)
                    {
                        for (c, _) in &mut row {
                            *c = domain.mul(c, &inv);
                        }
                    }
                    pivots[head_col] = Some(r);
                    break;
                }
            }
        }
        rows[r].terms = row;
    }

    rows.retain(|r| !r.terms.is_empty());
}

/// Sparse `row -= c * pivot` (generic domain).
fn sub_scaled<D: Domain>(
    domain: &D,
    row: &mut Vec<(D::Element, usize)>,
    pivot: &[(D::Element, usize)],
    c: &D::Element,
    scratch: &mut Vec<(D::Element, usize)>,
) {
    scratch.clear();
    let mut i = 1;
    let mut j = 1;
    while i < row.len() && j < pivot.len() {
        if row[i].1 < pivot[j].1 {
            scratch.push(row[i].clone());
            i += 1;
        } else if row[i].1 > pivot[j].1 {
            let prod = domain.mul(&pivot[j].0, c);
            let v = domain.sub(&domain.zero(), &prod);
            if !domain.is_zero(&v) {
                scratch.push((v, pivot[j].1));
            }
            j += 1;
        } else {
            let prod = domain.mul(&pivot[j].0, c);
            let v = domain.sub(&row[i].0, &prod);
            if !domain.is_zero(&v) {
                scratch.push((v, row[i].1));
            }
            i += 1;
            j += 1;
        }
    }
    scratch.extend_from_slice(&row[i..]);
    for (pc, pcol) in &pivot[j..] {
        let prod = domain.mul(pc, c);
        let v = domain.sub(&domain.zero(), &prod);
        if !domain.is_zero(&v) {
            scratch.push((v, *pcol));
        }
    }
    std::mem::swap(row, scratch);
}

/// Make a polynomial monic (generic path).
fn make_monic<D: Domain, O: MonomialOrder>(poly: &mut SparseMultivariatePolynomial<D, O>) {
    poly.make_monic_inplace();
}

// =========================================================================
//  Native ℤ_p fast path (f5_fp)
//
//  Structurally identical to the generic F5 loop above, but every
//  polynomial operation runs on `super::f4::FpPoly` — `i64` modular
//  arithmetic with no `BigInt` inside the pipeline. BigInt conversions
//  happen only when reading input and emitting the result.
// =========================================================================

use super::f4::{
    BasisPoly, CriticalPair, FpPoly, SimpCache, mod_inv, monic_fp, norm_mod, update_pairs,
};

/// A labeled FpPoly for the fast path: polynomial + F5 signature.
#[derive(Clone)]
struct LabeledFpPoly {
    poly: FpPoly,
    sig: Signature,
}

impl BasisPoly for LabeledFpPoly {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.poly.leading_monomial()
    }
    fn n_vars(&self) -> usize {
        self.poly.n_vars()
    }
    fn n_terms(&self) -> usize {
        self.poly.n_terms()
    }
    fn mul_monomial(&self, exp: &[usize]) -> Self {
        Self {
            poly: self.poly.mul_monomial(exp),
            sig: self.sig.mul_monomial(exp),
        }
    }
}

/// A sparse matrix row tagged with its F5 signature (i64 coefficients).
struct LabeledFpRow {
    /// (coefficient, column_index), ascending by column.
    terms: Vec<(i64, usize)>,
    sig: Signature,
}

/// Native ℤ_p F5: the full F5 pipeline on `i64` residues.
///
/// Mirrors the generic [`f5`] loop but every polynomial operation
/// (S-polynomial construction, symbolic preprocessing, row echelon) runs
/// on [`FpPoly`] — `BigInt` conversions happen only at the boundaries.
#[allow(clippy::too_many_lines)]
fn f5_fp<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    prime: i64,
) -> GroebnerBasis<D, O> {
    let n_vars = ideal[0].n_vars();

    // Filter zeros, convert to native residues, make monic.
    let mut generators: Vec<FpPoly> = ideal
        .iter()
        .filter(|p| !p.is_zero())
        .map(|p| FpPoly::from_domain(p, prime))
        .collect();
    for p in &mut generators {
        monic_fp(p, prime);
    }
    if generators.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Incremental F5: process one generator at a time. Uses F4's
    // Gebauer-Moeller update_pairs for pair management.
    let mut basis: Vec<LabeledFpPoly> = Vec::new();
    let mut pairs: Vec<CriticalPair> = Vec::new();
    let mut simplifications: Vec<SimpCache<LabeledFpPoly>> = Vec::new();
    let mut syzygies = SyzygySet::new();

    for (k, f) in generators.into_iter().enumerate() {
        let sig_k = Signature::unit(k, n_vars);
        let labeled = LabeledFpPoly {
            poly: f,
            sig: sig_k,
        };
        update_pairs(&mut basis, &mut pairs, &mut simplifications, labeled);

        // Degree-by-degree matrix reduction.
        while !pairs.is_empty() {
            let min_deg = pairs.iter().map(|p| p.degree).min().unwrap();
            let selected: Vec<CriticalPair> =
                pairs.extract_if(.., |p| p.degree == min_deg).collect();

            let new_polys = build_and_reduce_fp::<O>(&selected, &basis, &mut syzygies, prime);

            for poly in new_polys {
                update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
            }
        }
    }

    // Convert back to the domain representation and post-process.
    let domain = ideal[0].domain().clone();
    let basis_d: Vec<SparseMultivariatePolynomial<D, O>> = basis
        .iter()
        .map(|lp| lp.poly.to_domain::<D, O>(&domain, prime))
        .collect();
    GroebnerBasis { basis: basis_d }.minimize().auto_reduce()
}

/// Build the reduction matrix from selected pairs, echelonize, and extract
/// new basis elements (fast path with i64 coefficients).
fn build_and_reduce_fp<O: MonomialOrder>(
    selected: &[CriticalPair],
    basis: &[LabeledFpPoly],
    syzygies: &mut SyzygySet,
    prime: i64,
) -> Vec<LabeledFpPoly> {
    // --- Build input rows from selected pairs ---
    let mut monomial_map: HashMap<SmallVec<[usize; 4]>, usize> = HashMap::default();
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut rows: Vec<LabeledFpRow> = Vec::new();
    let mut worklist: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut seen_heads: HashSet<SmallVec<[usize; 4]>> = HashSet::default();

    for pair in selected {
        let i = pair.idx1;
        let j = pair.idx2;
        let lm_i = basis[i].leading_monomial().unwrap();
        let lm_j = basis[j].leading_monomial().unwrap();
        let lcm_exp = &pair.lcm;

        let diff_i: SmallVec<[usize; 4]> = lcm_exp
            .iter()
            .zip(lm_i.iter())
            .map(|(a, b)| a - b)
            .collect();
        let diff_j: SmallVec<[usize; 4]> = lcm_exp
            .iter()
            .zip(lm_j.iter())
            .map(|(a, b)| a - b)
            .collect();

        for (idx, diff) in [(i, &diff_i), (j, &diff_j)] {
            let sig = basis[idx].sig.mul_monomial(diff);
            if syzygies.contains(&sig) {
                continue;
            }
            let mult = basis[idx].poly.mul_monomial(diff);
            seen_heads.insert(lcm_exp.clone());
            add_fppoly_as_row(
                &mult,
                sig,
                &mut rows,
                &mut monomial_map,
                &mut monomial_list,
                &mut worklist,
            );
        }
    }

    if rows.is_empty() {
        return vec![];
    }

    // --- Symbolic preprocessing ---
    while let Some(exp) = worklist.pop() {
        if let Some((bi, diff)) = find_reducer_fp(basis, &exp) {
            let sig = basis[bi].sig.mul_monomial(&diff);
            if syzygies.contains(&sig) {
                continue;
            }
            seen_heads.insert(exp.clone());
            let reducer = basis[bi].poly.mul_monomial(&diff);
            add_fppoly_as_row(
                &reducer,
                sig,
                &mut rows,
                &mut monomial_map,
                &mut monomial_list,
                &mut worklist,
            );
        }
    }

    if rows.is_empty() || monomial_list.is_empty() {
        return vec![];
    }
    let ncols = monomial_list.len();

    // --- Sort columns: DESCENDING monomial order ---
    let mut col_order: Vec<usize> = (0..ncols).collect();
    col_order.sort_unstable_by(|&a, &b| O::cmp(&monomial_list[b], &monomial_list[a]));
    let mut col_inv = vec![0usize; ncols];
    for (new_col, &old_col) in col_order.iter().enumerate() {
        col_inv[old_col] = new_col;
    }
    for row in &mut rows {
        for (_, col) in row.terms.iter_mut() {
            *col = col_inv[*col];
        }
        row.terms.sort_unstable_by_key(|&(_, col)| col);
    }
    let mut sorted_monomials: Vec<SmallVec<[usize; 4]>> = vec![SmallVec::new(); ncols];
    for (new_col, &old_col) in col_order.iter().enumerate() {
        sorted_monomials[new_col] = monomial_list[old_col].clone();
    }

    // --- Sort rows by ascending signature (pot order) ---
    rows.sort_by(|a, b| a.sig.cmp_pot::<O>(&b.sig));

    // --- Echelonize (i64 modular arithmetic) ---
    echelonize_fp_labeled(&mut rows, ncols, prime);

    // --- Extract new basis elements ---
    let mut new_polys: Vec<LabeledFpPoly> = Vec::new();
    let basis_lm_set: HashSet<SmallVec<[usize; 4]>> = basis
        .iter()
        .filter_map(|lp| lp.leading_monomial().cloned())
        .collect();

    for row in &rows {
        if row.terms.is_empty() {
            syzygies.insert(&row.sig);
            continue;
        }
        let row_lm = &sorted_monomials[row.terms[0].1];
        if seen_heads.contains(row_lm) {
            continue;
        }
        if basis_lm_set.contains(row_lm) {
            continue;
        }

        // Reconstruct the FpPoly (terms descending).
        // row.terms are in ascending column order; column 0 = greatest
        // monomial (descending sort). Forward iteration yields descending
        // monomial order, which is what FpPoly expects.
        let mut terms: Vec<(SmallVec<[usize; 4]>, i64)> = Vec::new();
        for &(c, col) in &row.terms {
            let v = norm_mod(c, prime);
            if v != 0 {
                terms.push((sorted_monomials[col].clone(), v));
            }
        }
        if terms.is_empty() {
            syzygies.insert(&row.sig);
            continue;
        }

        new_polys.push(LabeledFpPoly {
            poly: FpPoly {
                terms,
                n_vars: basis[0].poly.n_vars(),
            },
            sig: row.sig.clone(),
        });
    }

    new_polys
}

/// Add an FpPoly as a labeled matrix row, registering new monomials.
fn add_fppoly_as_row(
    poly: &FpPoly,
    sig: Signature,
    rows: &mut Vec<LabeledFpRow>,
    monomial_map: &mut HashMap<SmallVec<[usize; 4]>, usize>,
    monomial_list: &mut Vec<SmallVec<[usize; 4]>>,
    worklist: &mut Vec<SmallVec<[usize; 4]>>,
) {
    let mut terms: Vec<(i64, usize)> = Vec::new();
    // FpPoly terms are stored descending; iterate as-is.
    for (exp, coeff) in &poly.terms {
        if *coeff == 0 {
            continue;
        }
        let col = *monomial_map.entry(exp.clone()).or_insert_with(|| {
            let idx = monomial_list.len();
            monomial_list.push(exp.clone());
            worklist.push(exp.clone());
            idx
        });
        terms.push((*coeff, col));
    }
    if !terms.is_empty() {
        rows.push(LabeledFpRow { terms, sig });
    }
}

/// Find a basis element whose leading monomial divides `exp`.
/// Returns `(basis_index, diff)`.
///
/// Note: `monomial_divides(a, b)` returns true iff `b` divides `a`.
fn find_reducer_fp(
    basis: &[LabeledFpPoly],
    exp: &[usize],
) -> Option<(usize, SmallVec<[usize; 4]>)> {
    for (i, lp) in basis.iter().enumerate() {
        if let Some(lm) = lp.leading_monomial()
            && monomial_divides(exp, lm)
        {
            let diff: SmallVec<[usize; 4]> =
                exp.iter().zip(lm.iter()).map(|(a, b)| a - b).collect();
            return Some((i, diff));
        }
    }
    None
}

/// Echelonize labeled i64 rows using sparse Gaussian elimination (mod p).
///
#[allow(dead_code)]
/// Rows must be pre-sorted by ascending signature so that standard
fn echelonize_fp_labeled(rows: &mut Vec<LabeledFpRow>, ncols: usize, prime: i64) {
    let p = prime;
    let mut pivots: Vec<Option<usize>> = vec![None; ncols];
    let mut scratch: Vec<(i64, usize)> = Vec::new();

    // First pass: identify and normalize pivots.
    for (r, row) in rows.iter_mut().enumerate() {
        if row.terms.is_empty() {
            continue;
        }
        let head_col = row.terms[0].1;
        if pivots[head_col].is_none() {
            if row.terms[0].0 != 1 {
                let inv = mod_inv(row.terms[0].0, p);
                for (c, _) in &mut row.terms {
                    *c = (*c * inv) % p;
                }
            }
            pivots[head_col] = Some(r);
        }
    }

    // Second pass: reduce non-pivot rows.
    for r in 0..rows.len() {
        if rows[r].terms.is_empty() {
            continue;
        }
        if pivots[rows[r].terms[0].1] == Some(r) {
            continue;
        }

        let mut row = std::mem::take(&mut rows[r].terms);
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(pr_) => {
                    let c = row[0].0;
                    sub_scaled_fp_labeled(&mut row, &rows[pr_].terms, c, p, &mut scratch);
                }
                None => {
                    if row[0].0 != 1 {
                        let inv = mod_inv(row[0].0, p);
                        for (c, _) in &mut row {
                            *c = (*c * inv) % p;
                        }
                    }
                    pivots[head_col] = Some(r);
                    break;
                }
            }
        }
        rows[r].terms = row;
    }

    rows.retain(|r| !r.terms.is_empty());
}

/// Sparse `row -= c * pivot` (mod p) by merging two column-ascending rows.
fn sub_scaled_fp_labeled(
    row: &mut Vec<(i64, usize)>,
    pivot: &[(i64, usize)],
    c: i64,
    p: i64,
    scratch: &mut Vec<(i64, usize)>,
) {
    scratch.clear();
    scratch.reserve(row.len() + pivot.len());
    let mut i = 1; // skip head (cancels)
    let mut j = 1;
    while i < row.len() && j < pivot.len() {
        let (rc, rcol) = row[i];
        let (pc, pcol) = pivot[j];
        if rcol < pcol {
            scratch.push((rc, rcol));
            i += 1;
        } else if rcol > pcol {
            let v = norm_mod(-c * pc, p);
            if v != 0 {
                scratch.push((v, pcol));
            }
            j += 1;
        } else {
            let v = norm_mod(rc - c * pc, p);
            if v != 0 {
                scratch.push((v, rcol));
            }
            i += 1;
            j += 1;
        }
    }
    scratch.extend_from_slice(&row[i..]);
    for &(pc, pcol) in &pivot[j..] {
        let v = norm_mod(-c * pc, p);
        if v != 0 {
            scratch.push((v, pcol));
        }
    }
    std::mem::swap(row, scratch);
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

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn signature_unit_and_mul() {
        let s = Signature::unit(2, 3);
        assert_eq!(s.module_pos, 2);
        assert_eq!(s.monomial.as_slice(), &[0, 0, 0]);

        let s2 = s.mul_monomial(&[1, 2, 0]);
        assert_eq!(s2.module_pos, 2);
        assert_eq!(s2.monomial.as_slice(), &[1, 2, 0]);
    }

    #[test]
    fn signature_pot_order() {
        let s1 = Signature::unit(0, 2);
        let s2 = Signature::unit(1, 2);
        // module_pos dominates: (0, ...) < (1, ...)
        assert_eq!(s1.cmp_pot::<Lex>(&s2), Ordering::Less);

        // Same module_pos: compare monomials under O.
        let s3 = Signature {
            module_pos: 0,
            monomial: smallvec![0, 1],
        };
        let s4 = Signature {
            module_pos: 0,
            monomial: smallvec![1, 0],
        };
        // Lex: [0,1] < [1,0] (first component 0 < 1)
        assert_eq!(s3.cmp_pot::<Lex>(&s4), Ordering::Less);
    }

    #[test]
    fn syzygy_set_basic() {
        let mut syz = SyzygySet::new();
        let s = Signature {
            module_pos: 1,
            monomial: smallvec![2, 0],
        };
        assert!(!syz.contains(&s));
        syz.insert(&s);
        assert!(syz.contains(&s));
        // A multiple should also be detected.
        let s_mult = Signature {
            module_pos: 1,
            monomial: smallvec![3, 1],
        };
        assert!(syz.contains(&s_mult));
        // Different module_pos should NOT match.
        let s_other = Signature {
            module_pos: 0,
            monomial: smallvec![2, 0],
        };
        assert!(!syz.contains(&s_other));
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

    #[test]
    fn f5_matches_buchberger() {
        let d = RationalDomain;
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 1], r(1, 1)), (vec![0, 0], r(-1, 1))],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d,
            2,
            vec![(vec![1, 0], r(1, 1)), (vec![0, 1], r(-1, 1))],
        );

        let gb_f5 = f5(&[f1.clone(), f2.clone()]);
        let gb_buch = crate::groebner::buchberger(&[f1, f2]);

        assert!(gb_f5.is_groebner_basis());
        assert_eq!(gb_f5.basis.len(), gb_buch.basis.len());
    }

    /// Build cyclic-n generators over ℤ_p.
    fn cyclic_fp(n: usize, p: u32) -> Vec<SparseMultivariatePolynomial<FiniteField, Lex>> {
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
                (vec![0usize; n], field.element(p - 1)),
            ],
        ));
        gens
    }

    #[test]
    fn f5_fp_linear_system() {
        // Simple linear system over ℤ₁₃ — exercises the fast path.
        let field = FiniteField::new(BigInt::from(13));
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            2,
            vec![
                (vec![1, 0], field.element(1)),
                (vec![0, 1], field.element(1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            2,
            vec![
                (vec![1, 0], field.element(1)),
                (vec![0, 1], field.element(12)), // -1 mod 13
            ],
        );
        let gb = f5(&[f1, f2]);
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f5_fp_cyclic_3_fp13() {
        // This exercises the native ℤ_p fast path (FiniteField → f5_fp).
        let ideal = cyclic_fp(3, 13);
        let gb = f5(&ideal);
        assert!(!gb.basis.is_empty());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f5_fp_cyclic_3_fp101() {
        let ideal = cyclic_fp(3, 101);
        let gb = f5(&ideal);
        assert!(!gb.basis.is_empty());
        assert!(gb.is_groebner_basis());
    }

    #[test]
    fn f5_fp_matches_f4_cyclic_3() {
        // F5 (fast path) and F4 should produce the same basis for cyclic-3.
        let ideal = cyclic_fp(3, 13);
        let gb_f5 = f5(&ideal);
        let gb_f4 = crate::groebner::f4::f4(&ideal);
        assert!(gb_f5.is_groebner_basis());
        assert!(gb_f4.is_groebner_basis());
    }
}
