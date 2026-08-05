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
use rayon::prelude::*;

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
    /// under the given monomial `order` (where `order.cmp(a, b)` returns
    /// `Less` when `a` precedes `b`).
    pub fn cmp_pot<O: MonomialOrder>(&self, other: &Self, order: &O) -> Ordering {
        self.module_pos
            .cmp(&other.module_pos)
            .then_with(|| order.cmp(&self.monomial, &other.monomial))
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
///
/// The stored monomials are bucketed by support mask (like F4's
/// `DivisorIndex`): a divisor of `t` must have its support inside
/// `support(t)`, so only the submasks of `support(t)` are examined. This
/// replaces the O(syzygies) linear scan that dominated F5 preprocessing
/// on large ideals (e.g. cyclic-6).
struct SyzygySet {
    /// module_pos → support-mask-bucketed leading monomials of known syzygies.
    lms: HashMap<usize, MonomialBucketSet>,
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
            .insert(sig.monomial.clone());
    }

    /// Check whether `sig` is (or is a multiple of) a known syzygy.
    fn contains(&self, sig: &Signature) -> bool {
        self.lms
            .get(&sig.module_pos)
            .is_some_and(|lms| lms.any_divisor_of(&sig.monomial))
    }
}

/// A set of monomials supporting fast "some stored monomial divides `exp`"
/// queries, bucketed by exact support mask.
struct MonomialBucketSet {
    /// Exact support mask of a monomial → monomials with that mask.
    buckets: HashMap<u64, Vec<SmallVec<[usize; 4]>>>,
}

impl Default for MonomialBucketSet {
    fn default() -> Self {
        Self::new()
    }
}

impl MonomialBucketSet {
    fn new() -> Self {
        Self {
            buckets: HashMap::default(),
        }
    }

    fn insert(&mut self, m: SmallVec<[usize; 4]>) {
        self.buckets.entry(support_mask(&m)).or_default().push(m);
    }

    /// True iff some stored monomial divides `exp`.
    fn any_divisor_of(&self, exp: &[usize]) -> bool {
        let mask = support_mask(exp);
        let mut sub = mask;
        loop {
            if let Some(ms) = self.buckets.get(&sub)
                && ms.iter().any(|m| monomial_divides(exp, m))
            {
                return true;
            }
            if sub == 0 {
                break;
            }
            sub = (sub - 1) & mask;
        }
        false
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
    // Requires p < 2^31 so that products of two residues fit in i64;
    // larger primes fall back to the generic domain path.
    if let Some(ff) =
        (ideal[0].domain() as &dyn std::any::Any).downcast_ref::<ocas_domain::FiniteField>()
        && ff.prime_u64() < (1u64 << 31)
    {
        let prime = ff.prime_u64() as i64;
        // Packed fast path: u128 SWAR monomials when the ideal fits the
        // 8-variable / 15-bit exponent contract; otherwise the i64 path.
        if packed_eligible(ideal) {
            return f5_fp_packed(ideal, prime);
        }
        return f5_fp(ideal, prime);
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
    let order = basis[0].poly.order.clone();

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
    col_order.sort_unstable_by(|&a, &b| order.cmp(&monomial_list[b], &monomial_list[a]));
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
    rows.sort_by(|a, b| a.sig.cmp_pot::<O>(&b.sig, &order));

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
    BasisPoly, CriticalPair, DivisorIndex, FpPoly, SimpCache, domain_to_i64_fp, i64_to_domain_fp,
    mod_inv, monic_fp, norm_mod, support_mask, update_pairs,
};
use super::packed::{PackedMono, PackedMonomialBucketSet, PackedSig};

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

/// A sparse matrix row tagged with its F5 signature (i32 coefficients).
///
/// Coefficients are residues in `[0, p)` with `p < 2^31`, so they fit in
/// `i32`; products are widened to `i64` inside the echelon arithmetic.
/// Generic over the signature type so the packed ℤ_p path can carry a
/// [`PackedSig`] while the SmallVec path carries a [`Signature`]; the
/// echelon machinery only touches `terms`.
struct LabeledFpRow<S = Signature> {
    /// (coefficient, column_index), ascending by column.
    terms: Vec<(i32, usize)>,
    sig: S,
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
    let order = ideal[0].order.clone();
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
    // Reducer divisor index over basis leading monomials; the basis only
    // grows, so the index is extended in lockstep after each push.
    let mut div_index = DivisorIndex::new();

    for (k, f) in generators.into_iter().enumerate() {
        let sig_k = Signature::unit(k, n_vars);
        let labeled = LabeledFpPoly {
            poly: f,
            sig: sig_k,
        };
        update_pairs(&mut basis, &mut pairs, &mut simplifications, labeled);
        if let Some(lm) = basis.last().and_then(|lp| lp.leading_monomial()) {
            div_index.push(lm, basis.len() - 1);
        }

        // Degree-by-degree matrix reduction.
        while !pairs.is_empty() {
            let min_deg = pairs.iter().map(|p| p.degree).min().unwrap();
            let selected: Vec<CriticalPair> =
                pairs.extract_if(.., |p| p.degree == min_deg).collect();
            let new_polys = build_and_reduce_fp::<O>(
                &selected,
                &basis,
                &mut syzygies,
                &div_index,
                prime,
                &order,
            );

            for poly in new_polys {
                update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
                if let Some(lm) = basis.last().and_then(|lp| lp.leading_monomial()) {
                    div_index.push(lm, basis.len() - 1);
                }
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
    div_index: &DivisorIndex,
    prime: i64,
    order: &O,
) -> Vec<LabeledFpPoly> {
    // --- Build input rows from selected pairs ---
    // Capacity preallocation: each pair contributes up to 2 rows and the
    // column table typically grows a few monomials per pair.
    let map_cap = selected.len() * 4;
    let mut monomial_map: HashMap<SmallVec<[usize; 4]>, usize> =
        HashMap::with_capacity_and_hasher(map_cap, Default::default());
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::with_capacity(map_cap);
    let mut rows: Vec<LabeledFpRow> = Vec::with_capacity(selected.len() * 2);
    let mut worklist: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut seen_heads: HashSet<SmallVec<[usize; 4]>> = HashSet::default();

    // Row content per pair is independent: compute the raw (monomial,
    // coefficient) rows in parallel, then register monomials into the
    // shared column table in the original pair order (identical result,
    // deterministic).
    type RawPairRows = Vec<(Signature, Vec<(SmallVec<[usize; 4]>, i64)>)>;
    let raw_rows: Vec<RawPairRows> = selected
        .par_iter()
        .map(|pair| {
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

            let mut out: RawPairRows = Vec::new();
            for (idx, diff) in [(i, &diff_i), (j, &diff_j)] {
                let sig = basis[idx].sig.mul_monomial(diff);
                if syzygies.contains(&sig) {
                    continue;
                }
                let mult: FpPoly = basis[idx].poly.mul_monomial(diff);
                let terms: Vec<(SmallVec<[usize; 4]>, i64)> = mult
                    .terms
                    .iter()
                    .filter(|t| t.1 != 0)
                    .map(|t| (t.0.clone(), t.1))
                    .collect();
                if !terms.is_empty() {
                    out.push((sig, terms));
                }
            }
            out
        })
        .collect();

    for (pair, raw) in selected.iter().zip(raw_rows) {
        let lcm_exp = &pair.lcm;
        for (sig, terms) in raw {
            seen_heads.insert(lcm_exp.clone());
            let mut mapped: Vec<(i32, usize)> = Vec::with_capacity(terms.len());
            for (exp, coeff) in terms {
                let col = match monomial_map.get(&exp) {
                    Some(&c) => c,
                    None => {
                        let idx = monomial_list.len();
                        monomial_list.push(exp.clone());
                        worklist.push(exp.clone());
                        monomial_map.insert(exp, idx);
                        idx
                    }
                };
                mapped.push((coeff as i32, col));
            }
            rows.push(LabeledFpRow { terms: mapped, sig });
        }
    }

    if rows.is_empty() {
        return vec![];
    }

    // --- Symbolic preprocessing ---
    while let Some(exp) = worklist.pop() {
        if let Some((bi, diff)) = find_reducer_fp(div_index, basis, &exp) {
            let sig = basis[bi].sig.mul_monomial(&diff);
            if syzygies.contains(&sig) {
                continue;
            }
            seen_heads.insert(exp.clone());
            add_scaled_fppoly_as_row(
                &basis[bi].poly,
                &diff,
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
    col_order.sort_unstable_by(|&a, &b| order.cmp(&monomial_list[b], &monomial_list[a]));
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
    rows.sort_by(|a, b| a.sig.cmp_pot::<O>(&b.sig, order));

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
            let v = norm_mod(c as i64, prime);
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

/// Register `poly * x^diff` as a labeled matrix row without materializing
/// the product: each term's exponent is shifted by `diff` and registered
/// directly into the column table (get-before-insert avoids cloning the
/// key on hits, which dominate). Monomials are registered in the same
/// order as `add_fppoly_as_row` on the materialized product, so the
/// column assignment is identical.
fn add_scaled_fppoly_as_row(
    poly: &FpPoly,
    diff: &[usize],
    sig: Signature,
    rows: &mut Vec<LabeledFpRow>,
    monomial_map: &mut HashMap<SmallVec<[usize; 4]>, usize>,
    monomial_list: &mut Vec<SmallVec<[usize; 4]>>,
    worklist: &mut Vec<SmallVec<[usize; 4]>>,
) {
    let mut terms: Vec<(i32, usize)> = Vec::new();
    // FpPoly terms are stored descending; iterate as-is.
    // Stack buffer: avoids a heap allocation per term for n_vars <= 16
    // (SmallVec<[usize; 4]> spills to the heap beyond 4 variables, and
    // this loop is the F5 hot path). The map lookup borrows the buffer as
    // `&[usize]` (SmallVec: Borrow<[T]>); the owned key is built only on
    // the (minority) miss path.
    let mut buf: SmallVec<[usize; 16]> = SmallVec::new();
    for (exp, coeff) in &poly.terms {
        if *coeff == 0 {
            continue;
        }
        buf.clear();
        for (v, &dv) in diff.iter().enumerate() {
            buf.push(dv + exp.get(v).copied().unwrap_or(0));
        }
        let col = match monomial_map.get(&buf[..]) {
            Some(&c) => c,
            None => {
                let key: SmallVec<[usize; 4]> = SmallVec::from_slice(&buf);
                let idx = monomial_list.len();
                monomial_list.push(key.clone());
                worklist.push(key.clone());
                monomial_map.insert(key, idx);
                idx
            }
        };
        terms.push((*coeff as i32, col));
    }
    if !terms.is_empty() {
        rows.push(LabeledFpRow { terms, sig });
    }
}

/// Find a basis element whose leading monomial divides `exp`, via the
/// [`DivisorIndex`] over basis leading monomials.
/// Returns `(basis_index, diff)`.
///
/// Selection semantics match the original linear scan: the lowest basis
/// index whose leading monomial divides `exp` wins.
///
/// Note: `monomial_divides(a, b)` returns true iff `b` divides `a`.
fn find_reducer_fp(
    index: &DivisorIndex,
    basis: &[LabeledFpPoly],
    exp: &[usize],
) -> Option<(usize, SmallVec<[usize; 4]>)> {
    let mask = support_mask(exp);
    let mut best: Option<usize> = None;
    // Enumerate all submasks of `mask`, including `mask` itself and 0.
    let mut sub = mask;
    loop {
        if let Some(ids) = index.buckets.get(&sub) {
            for &bi in ids {
                if let Some(lm) = basis[bi].leading_monomial()
                    && monomial_divides(exp, lm)
                {
                    match best {
                        Some(b) if b <= bi => {}
                        _ => best = Some(bi),
                    }
                }
            }
        }
        if sub == 0 {
            break;
        }
        sub = (sub - 1) & mask;
    }
    best.map(|bi| {
        let lm = basis[bi].leading_monomial().unwrap();
        let diff: SmallVec<[usize; 4]> = exp.iter().zip(lm.iter()).map(|(a, b)| a - b).collect();
        (bi, diff)
    })
}

/// Where a pivot row lives after the pass-1 pivot extraction: in the
/// `pivot_store` (a pass-1 pivot, moved out of `rows` so Phase A workers
/// never alias) or still in `rows` (a pivot claimed during Phase B).
#[derive(Clone, Copy)]
enum PivotLoc {
    Store(usize),
    Row(usize),
}

/// Echelonize labeled i32 rows using sparse Gaussian elimination (mod p).
///
/// Rows must be pre-sorted by ascending signature so that standard
/// signature-rewriting pivot claims are reproduced exactly.
///
/// Two-phase, no-clone: pass-1 pivot rows are *moved* into `pivot_store`
/// (zero cloning) and Phase A reduces every non-pivot row in parallel
/// against that read-only store; Phase B resumes interrupted rows in the
/// original order and claims new pivots in place. The echelon form is
/// bit-identical to a fully sequential run, and the final row order is
/// the original order minus zero rows.
fn echelonize_fp_labeled<S: Send + Sync>(
    rows: &mut Vec<LabeledFpRow<S>>,
    ncols: usize,
    prime: i64,
) {
    let p = prime;
    let mut pivots: Vec<Option<PivotLoc>> = vec![None; ncols];

    // First pass: identify and normalize pivots.
    for (r, row) in rows.iter_mut().enumerate() {
        if row.terms.is_empty() {
            continue;
        }
        let head_col = row.terms[0].1;
        if pivots[head_col].is_none() {
            if row.terms[0].0 != 1 {
                let inv = mod_inv(row.terms[0].0 as i64, p);
                for (c, _) in &mut row.terms {
                    *c = ((*c as i64 * inv) % p) as i32;
                }
            }
            pivots[head_col] = Some(PivotLoc::Row(r));
        }
    }

    // Move pass-1 pivot rows out of `rows` into the store, recording their
    // original row numbers (for Phase B pivot lookups and final reinsertion).
    // `rows` slots become empty; Phase A/B skip them, and the extraction
    // order is unchanged once the terms are moved back.
    let mut pivot_store: Vec<Vec<(i32, usize)>> = Vec::new();
    let mut pivot_orig_row: Vec<usize> = Vec::new();
    for slot in pivots.iter_mut() {
        if let Some(PivotLoc::Row(r)) = *slot {
            pivot_store.push(std::mem::take(&mut rows[r].terms));
            pivot_orig_row.push(r);
            *slot = Some(PivotLoc::Store(pivot_store.len() - 1));
        }
    }

    // Phase A (parallel): reduce every non-pivot row against the pass-1
    // pivots (read-only store, zero cloning, no aliasing). A row whose head
    // reaches a column without a pass-1 pivot stops there; the reduction is
    // resumed in Phase B in the original row order, so the pivot-claim
    // sequence — and therefore the echelon form — is bit-identical to a
    // fully sequential run.
    rows.par_iter_mut().for_each(|labeled| {
        if labeled.terms.is_empty() {
            return;
        }
        let mut row = std::mem::take(&mut labeled.terms);
        let mut scratch: Vec<(i32, usize)> = Vec::new();
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(PivotLoc::Store(si)) => {
                    let c = row[0].0;
                    sub_scaled_fp_labeled(&mut row, &pivot_store[si], c, p, &mut scratch);
                }
                _ => break, // no pass-1 pivot: defer to Phase B
            }
        }
        labeled.terms = row;
    });

    // Phase B (sequential): continue interrupted rows in original order;
    // new pivots are claimed exactly as in a sequential run.
    let mut scratch: Vec<(i32, usize)> = Vec::new();
    for r in 0..rows.len() {
        if rows[r].terms.is_empty() {
            continue;
        }
        if matches!(pivots[rows[r].terms[0].1], Some(PivotLoc::Row(pr)) if pr == r) {
            continue;
        }

        let mut row = std::mem::take(&mut rows[r].terms);
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(PivotLoc::Store(si)) => {
                    let c = row[0].0;
                    sub_scaled_fp_labeled(&mut row, &pivot_store[si], c, p, &mut scratch);
                }
                Some(PivotLoc::Row(pr)) => {
                    let c = row[0].0;
                    sub_scaled_fp_labeled(&mut row, &rows[pr].terms, c, p, &mut scratch);
                }
                None => {
                    if row[0].0 != 1 {
                        let inv = mod_inv(row[0].0 as i64, p);
                        for (c, _) in &mut row {
                            *c = ((*c as i64 * inv) % p) as i32;
                        }
                    }
                    pivots[head_col] = Some(PivotLoc::Row(r));
                    break;
                }
            }
        }
        rows[r].terms = row;
    }

    // Move pass-1 pivot rows back to their original positions, then drop
    // zero rows: the surviving row order matches the sequential run.
    for (si, &orig) in pivot_orig_row.iter().enumerate() {
        rows[orig].terms = std::mem::take(&mut pivot_store[si]);
    }
    rows.retain(|r| !r.terms.is_empty());
}

/// Sparse `row -= c * pivot` (mod p) by merging two column-ascending rows.
/// Coefficients are i32 residues; products are computed in i64.
fn sub_scaled_fp_labeled(
    row: &mut Vec<(i32, usize)>,
    pivot: &[(i32, usize)],
    c: i32,
    p: i64,
    scratch: &mut Vec<(i32, usize)>,
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
            let v = norm_mod(-(c as i64) * (pc as i64), p) as i32;
            if v != 0 {
                scratch.push((v, pcol));
            }
            j += 1;
        } else {
            let v = norm_mod(rc as i64 - (c as i64) * (pc as i64), p) as i32;
            if v != 0 {
                scratch.push((v, rcol));
            }
            i += 1;
            j += 1;
        }
    }
    scratch.extend_from_slice(&row[i..]);
    for &(pc, pcol) in &pivot[j..] {
        let v = norm_mod(-(c as i64) * (pc as i64), p) as i32;
        if v != 0 {
            scratch.push((v, pcol));
        }
    }
    std::mem::swap(row, scratch);
}

// =========================================================================
//  Packed ℤ_p fast path (f5_fp_packed)
//
//  Same F5 pipeline as f5_fp, but monomials are u128 SWAR-packed
//  (super::packed::PackedMono): exponent addition, divisibility, lcm and
//  support masks become single u128 operations, and HashMap keys are
//  Copy values — no SmallVec cloning and no heap allocation per term.
//  Eligible when n_vars ≤ 8 and every input exponent < 2^15 (the SWAR
//  formulas' exactness contract; see packed.rs).
// =========================================================================

/// Whether the packed fast path may handle `ideal`: at most 8 variables
/// and every input exponent below 2^15 (the packed SWAR formulas are the
/// standard 15-bit-field tricks — exact only for fields < 2^15; all F5
/// matrix monomials are bounded fieldwise by the input exponents, so the
/// input bound covers every intermediate monomial).
pub(crate) fn packed_eligible<D: Domain, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
) -> bool {
    ideal.iter().all(|p| p.n_vars() <= 8)
        && ideal
            .iter()
            .all(|p| p.terms_ref().keys().all(|e| e.iter().all(|&x| x < 1 << 15)))
}

/// A native ℤ_p polynomial with packed monomials.
///
/// Terms are stored as a `Vec` sorted by **descending** monomial order
/// (leading term first) with coefficients in `[0, p)`; monomials are
/// [`PackedMono`] values. `lm_sv` is a SmallVec sidecar of the leading
/// monomial, computed once at construction, kept only for
/// [`BasisPoly`]/[`DivisorIndex`] compatibility.
#[derive(Debug, Clone)]
struct PackedFpPoly {
    /// Terms sorted descending by monomial order; coeffs in `[0, p)`.
    terms: Vec<(PackedMono, i64)>,
    n_vars: usize,
    lm_sv: Option<SmallVec<[usize; 4]>>,
}

/// Sidecar of the leading monomial of `terms` (smallest allocation-free
/// form compatible with `DivisorIndex::push` and `BasisPoly`).
fn lm_smallvec(terms: &[(PackedMono, i64)], n_vars: usize) -> Option<SmallVec<[usize; 4]>> {
    terms
        .first()
        .map(|t| SmallVec::from_slice(&t.0.unpack_sv(n_vars)))
}

impl PackedFpPoly {
    /// Leading monomial as a packed value (the authoritative one; the
    /// `lm_sv` sidecar is only for `BasisPoly` compatibility).
    fn leading_monomial_packed(&self) -> Option<PackedMono> {
        self.terms.first().map(|t| t.0)
    }

    /// Convert a domain polynomial to native residues (one-time cost at
    /// the pipeline boundary). Mirrors `FpPoly::from_domain`.
    fn from_domain<D: Domain + 'static, O: MonomialOrder>(
        p: &SparseMultivariatePolynomial<D, O>,
        prime: i64,
    ) -> Self {
        let mut terms: Vec<(PackedMono, i64)> = Vec::with_capacity(p.n_terms());
        // sorted_terms ascending → rev gives descending (leading first).
        for (exp, coeff) in p.sorted_terms().iter().rev() {
            let c = norm_mod(domain_to_i64_fp::<D>(coeff, prime), prime);
            if c != 0 {
                terms.push((PackedMono::pack(exp), c));
            }
        }
        let n_vars = p.n_vars();
        let lm_sv = lm_smallvec(&terms, n_vars);
        Self {
            terms,
            n_vars,
            lm_sv,
        }
    }

    /// Convert back to the BigInt-backed domain (one-time cost per
    /// surviving basis element and at the end of the pipeline). Mirrors
    /// `FpPoly::to_domain`.
    fn to_domain<D: Domain + 'static, O: MonomialOrder>(
        &self,
        domain: &D,
        prime: i64,
    ) -> SparseMultivariatePolynomial<D, O> {
        let mut poly = SparseMultivariatePolynomial::new(domain.clone(), self.n_vars);
        for (exp, c) in &self.terms {
            poly.append_monomial(
                i64_to_domain_fp::<D>(domain, *c, prime),
                &exp.unpack_sv(self.n_vars),
            );
        }
        poly
    }

    /// Multiply by a packed monomial (single u128 add per term).
    fn mul_monomial_packed(&self, diff: PackedMono) -> Self {
        let terms: Vec<(PackedMono, i64)> =
            self.terms.iter().map(|(e, c)| (e.add(diff), *c)).collect();
        let n_vars = self.n_vars;
        let lm_sv = lm_smallvec(&terms, n_vars);
        Self {
            terms,
            n_vars,
            lm_sv,
        }
    }
}

impl BasisPoly for PackedFpPoly {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.lm_sv.as_ref()
    }
    fn n_vars(&self) -> usize {
        self.n_vars
    }
    fn n_terms(&self) -> usize {
        self.terms.len()
    }
    fn mul_monomial(&self, exp: &[usize]) -> Self {
        self.mul_monomial_packed(PackedMono::pack(exp))
    }
}

/// Make a `PackedFpPoly` monic (scale by the modular inverse of the
/// leading coefficient). Mirrors `monic_fp`; only coefficients change.
fn monic_packed(p: &mut PackedFpPoly, prime: i64) {
    if let Some(&(_, lc)) = p.terms.first()
        && lc != 1
    {
        let inv = mod_inv(lc, prime);
        for (_, c) in &mut p.terms {
            *c = norm_mod(*c * inv, prime);
        }
    }
}

/// A labeled packed polynomial for the fast path: polynomial + packed
/// F5 signature.
#[derive(Clone)]
struct LabeledPackedFpPoly {
    poly: PackedFpPoly,
    sig: PackedSig,
}

impl BasisPoly for LabeledPackedFpPoly {
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
            sig: self.sig.mul_monomial(PackedMono::pack(exp)),
        }
    }
}

/// Syzygy tracking for packed signatures: `PackedMonomialBucketSet` per
/// module position (the `PackedSig` analogue of [`SyzygySet`]).
struct PackedSyzygySet {
    /// module_pos → support-mask-bucketed leading monomials of known syzygies.
    lms: HashMap<usize, PackedMonomialBucketSet>,
}

impl PackedSyzygySet {
    fn new() -> Self {
        Self {
            lms: HashMap::default(),
        }
    }

    /// Record that signature `sig` produces a zero reduction.
    fn insert(&mut self, sig: &PackedSig) {
        self.lms.entry(sig.pos).or_default().insert(sig.mono);
    }

    /// Check whether `sig` is (or is a multiple of) a known syzygy.
    fn contains(&self, sig: &PackedSig) -> bool {
        self.lms
            .get(&sig.pos)
            .is_some_and(|lms| lms.any_divisor_of(sig.mono))
    }
}

/// Native ℤ_p F5 on packed monomials: mirrors [`f5_fp`] exactly, with
/// `PackedFpPoly`/`PackedSig` in place of `FpPoly`/`Signature`. BigInt
/// conversions happen only at the boundaries.
#[allow(clippy::too_many_lines)]
fn f5_fp_packed<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    prime: i64,
) -> GroebnerBasis<D, O> {
    let order = ideal[0].order.clone();
    let mut generators: Vec<PackedFpPoly> = ideal
        .iter()
        .filter(|p| !p.is_zero())
        .map(|p| PackedFpPoly::from_domain(p, prime))
        .collect();
    for p in &mut generators {
        monic_packed(p, prime);
    }
    if generators.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Incremental F5: process one generator at a time. Uses F4's
    // Gebauer-Moeller update_pairs for pair management.
    let mut basis: Vec<LabeledPackedFpPoly> = Vec::new();
    let mut pairs: Vec<CriticalPair> = Vec::new();
    let mut simplifications: Vec<SimpCache<LabeledPackedFpPoly>> = Vec::new();
    let mut syzygies = PackedSyzygySet::new();
    // Reducer divisor index over basis leading monomials; the basis only
    // grows, so the index is extended in lockstep after each push.
    let mut div_index = DivisorIndex::new();

    for (k, f) in generators.into_iter().enumerate() {
        let sig_k = PackedSig::unit(k);
        let labeled = LabeledPackedFpPoly {
            poly: f,
            sig: sig_k,
        };
        update_pairs(&mut basis, &mut pairs, &mut simplifications, labeled);
        if let Some(lm) = basis.last().and_then(|lp| lp.leading_monomial()) {
            div_index.push(lm, basis.len() - 1);
        }

        // Degree-by-degree matrix reduction.
        while !pairs.is_empty() {
            let min_deg = pairs.iter().map(|p| p.degree).min().unwrap();
            let selected: Vec<CriticalPair> =
                pairs.extract_if(.., |p| p.degree == min_deg).collect();
            let new_polys = build_and_reduce_fp_packed::<O>(
                &selected,
                &basis,
                &mut syzygies,
                &div_index,
                prime,
                &order,
            );

            for poly in new_polys {
                update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
                if let Some(lm) = basis.last().and_then(|lp| lp.leading_monomial()) {
                    div_index.push(lm, basis.len() - 1);
                }
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
/// new basis elements (packed-monomial fast path). Row construction and
/// extraction mirror [`build_and_reduce_fp`] with `PackedMono` keys.
fn build_and_reduce_fp_packed<O: MonomialOrder>(
    selected: &[CriticalPair],
    basis: &[LabeledPackedFpPoly],
    syzygies: &mut PackedSyzygySet,
    div_index: &DivisorIndex,
    prime: i64,
    order: &O,
) -> Vec<LabeledPackedFpPoly> {
    let n_vars = basis[0].poly.n_vars();

    // --- Build input rows from selected pairs ---
    // Capacity preallocation: each pair contributes up to 2 rows and the
    // column table typically grows a few monomials per pair.
    let map_cap = selected.len() * 4;
    let mut monomial_map: HashMap<PackedMono, usize> =
        HashMap::with_capacity_and_hasher(map_cap, Default::default());
    let mut monomial_list: Vec<PackedMono> = Vec::with_capacity(map_cap);
    let mut rows: Vec<LabeledFpRow<PackedSig>> = Vec::with_capacity(selected.len() * 2);
    let mut worklist: Vec<PackedMono> = Vec::new();
    let mut seen_heads: HashSet<PackedMono> = HashSet::default();

    // Row content per pair is independent: compute the raw (monomial,
    // coefficient) rows in parallel, then register monomials into the
    // shared column table in the original pair order (identical result,
    // deterministic).
    type RawPairRows = Vec<(PackedSig, Vec<(PackedMono, i64)>)>;
    let raw_rows: Vec<RawPairRows> = selected
        .par_iter()
        .map(|pair| {
            let i = pair.idx1;
            let j = pair.idx2;
            let lm_i = basis[i].poly.leading_monomial_packed().unwrap();
            let lm_j = basis[j].poly.leading_monomial_packed().unwrap();
            let lcm_exp = PackedMono::pack(&pair.lcm);

            let diff_i = lcm_exp.sub(lm_i);
            let diff_j = lcm_exp.sub(lm_j);

            let mut out: RawPairRows = Vec::new();
            for (idx, diff) in [(i, diff_i), (j, diff_j)] {
                let sig = basis[idx].sig.mul_monomial(diff);
                if syzygies.contains(&sig) {
                    continue;
                }
                let mult: PackedFpPoly = basis[idx].poly.mul_monomial_packed(diff);
                let terms: Vec<(PackedMono, i64)> = mult
                    .terms
                    .iter()
                    .filter(|t| t.1 != 0)
                    .map(|t| (t.0, t.1))
                    .collect();
                if !terms.is_empty() {
                    out.push((sig, terms));
                }
            }
            out
        })
        .collect();

    for (pair, raw) in selected.iter().zip(raw_rows) {
        let lcm_packed = PackedMono::pack(&pair.lcm);
        for (sig, terms) in raw {
            seen_heads.insert(lcm_packed);
            let mut mapped: Vec<(i32, usize)> = Vec::with_capacity(terms.len());
            for (exp, coeff) in terms {
                let col = match monomial_map.get(&exp) {
                    Some(&c) => c,
                    None => {
                        let idx = monomial_list.len();
                        monomial_list.push(exp);
                        worklist.push(exp);
                        monomial_map.insert(exp, idx);
                        idx
                    }
                };
                mapped.push((coeff as i32, col));
            }
            rows.push(LabeledFpRow { terms: mapped, sig });
        }
    }

    if rows.is_empty() {
        return vec![];
    }

    // --- Symbolic preprocessing ---
    while let Some(exp) = worklist.pop() {
        if let Some((bi, diff)) = find_reducer_packed(div_index, basis, exp) {
            let sig = basis[bi].sig.mul_monomial(diff);
            if syzygies.contains(&sig) {
                continue;
            }
            seen_heads.insert(exp);
            add_scaled_packed_as_row(
                &basis[bi].poly,
                diff,
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
    col_order.sort_unstable_by(|&a, &b| {
        order.cmp(
            &monomial_list[b].unpack_sv(n_vars),
            &monomial_list[a].unpack_sv(n_vars),
        )
    });
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
    let mut sorted_monomials: Vec<PackedMono> = vec![PackedMono(0); ncols];
    for (new_col, &old_col) in col_order.iter().enumerate() {
        sorted_monomials[new_col] = monomial_list[old_col];
    }

    // --- Sort rows by ascending signature (pot order) ---
    rows.sort_by(|a, b| a.sig.cmp_pot::<O>(&b.sig, order, n_vars));

    // --- Echelonize (i64 modular arithmetic; signature type is opaque) ---
    echelonize_fp_labeled(&mut rows, ncols, prime);

    // --- Extract new basis elements ---
    let mut new_polys: Vec<LabeledPackedFpPoly> = Vec::new();
    let basis_lm_set: HashSet<PackedMono> = basis
        .iter()
        .filter_map(|lp| lp.poly.leading_monomial_packed())
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

        // Reconstruct the PackedFpPoly (terms descending).
        // row.terms are in ascending column order; column 0 = greatest
        // monomial (descending sort). Forward iteration yields descending
        // monomial order, which is what PackedFpPoly expects.
        let mut terms: Vec<(PackedMono, i64)> = Vec::new();
        for &(c, col) in &row.terms {
            let v = norm_mod(c as i64, prime);
            if v != 0 {
                terms.push((sorted_monomials[col], v));
            }
        }
        if terms.is_empty() {
            syzygies.insert(&row.sig);
            continue;
        }

        new_polys.push(LabeledPackedFpPoly {
            poly: PackedFpPoly {
                n_vars,
                lm_sv: lm_smallvec(&terms, n_vars),
                terms,
            },
            sig: row.sig.clone(),
        });
    }

    new_polys
}

/// Register `poly * x^diff` as a labeled matrix row without materializing
/// the product: each term's packed exponent is shifted by `diff` (one
/// u128 add) and registered directly into the column table. Keys are
/// `Copy` [`PackedMono`] values, so the get-before-insert hit path has
/// zero cloning and zero allocation.
fn add_scaled_packed_as_row(
    poly: &PackedFpPoly,
    diff: PackedMono,
    sig: PackedSig,
    rows: &mut Vec<LabeledFpRow<PackedSig>>,
    monomial_map: &mut HashMap<PackedMono, usize>,
    monomial_list: &mut Vec<PackedMono>,
    worklist: &mut Vec<PackedMono>,
) {
    let mut terms: Vec<(i32, usize)> = Vec::new();
    // PackedFpPoly terms are stored descending; iterate as-is.
    for (exp, coeff) in &poly.terms {
        if *coeff == 0 {
            continue;
        }
        let key = exp.add(diff);
        let col = match monomial_map.get(&key) {
            Some(&c) => c,
            None => {
                let idx = monomial_list.len();
                monomial_list.push(key);
                worklist.push(key);
                monomial_map.insert(key, idx);
                idx
            }
        };
        terms.push((*coeff as i32, col));
    }
    if !terms.is_empty() {
        rows.push(LabeledFpRow { terms, sig });
    }
}

/// Find a basis element whose leading monomial divides `exp`, via the
/// [`DivisorIndex`] over basis leading monomials, with SWAR divisibility
/// on the packed leading monomials. Returns `(basis_index, diff)`.
///
/// Selection semantics match [`find_reducer_fp`]: the lowest basis index
/// whose leading monomial divides `exp` wins.
fn find_reducer_packed(
    index: &DivisorIndex,
    basis: &[LabeledPackedFpPoly],
    exp: PackedMono,
) -> Option<(usize, PackedMono)> {
    let mask = exp.support_mask();
    let mut best: Option<usize> = None;
    // Enumerate all submasks of `mask`, including `mask` itself and 0.
    let mut sub = mask;
    loop {
        if let Some(ids) = index.buckets.get(&sub) {
            for &bi in ids {
                if let Some(lm) = basis[bi].poly.leading_monomial_packed()
                    && exp.divides(lm)
                {
                    match best {
                        Some(b) if b <= bi => {}
                        _ => best = Some(bi),
                    }
                }
            }
        }
        if sub == 0 {
            break;
        }
        sub = (sub - 1) & mask;
    }
    best.map(|bi| {
        let lm = basis[bi].poly.leading_monomial_packed().unwrap();
        (bi, exp.sub(lm))
    })
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
        assert_eq!(s1.cmp_pot::<Lex>(&s2, &Lex), Ordering::Less);

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
        assert_eq!(s3.cmp_pot::<Lex>(&s4, &Lex), Ordering::Less);
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

    #[test]
    fn packed_eligible_checks() {
        // cyclic-6 ℤ₁₃ fits the packed contract (6 vars, exponents ≤ 1).
        let ideal6 = cyclic_fp(6, 13);
        assert!(packed_eligible(&ideal6));

        // A 9-variable ideal does not fit.
        let field = FiniteField::new(BigInt::from(13));
        let mut terms = Vec::new();
        for i in 0..9 {
            let mut e = vec![0usize; 9];
            e[i] = 1;
            terms.push((e, field.element(1)));
        }
        let big = SparseMultivariatePolynomial::<_, Lex>::from_terms(field.clone(), 9, terms);
        let p2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            9,
            vec![
                (vec![1, 1, 0, 0, 0, 0, 0, 0, 0], field.element(1)),
                (vec![0usize; 9], field.element(12)),
            ],
        );
        assert!(!packed_eligible(&[big, p2]));

        // An exponent at or beyond 2^15 (the SWAR exactness bound) does
        // not fit, even with ≤ 8 variables.
        let f = FiniteField::new(BigInt::from(13));
        let h1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            f.clone(),
            2,
            vec![
                (vec![70000, 0], f.element(1)),
                (vec![0usize; 2], f.element(12)),
            ],
        );
        let h2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            f.clone(),
            2,
            vec![(vec![0, 1], f.element(1)), (vec![0usize; 2], f.element(12))],
        );
        assert!(!packed_eligible(&[h1, h2]));
        // Boundary: exactly 2^15 - 1 still fits.
        let b1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            f.clone(),
            2,
            vec![
                (vec![(1 << 15) - 1, 0], f.element(1)),
                (vec![0usize; 2], f.element(12)),
            ],
        );
        let b2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            f.clone(),
            2,
            vec![(vec![0, 1], f.element(1)), (vec![0usize; 2], f.element(12))],
        );
        assert!(packed_eligible(&[b1, b2]));
    }
}
