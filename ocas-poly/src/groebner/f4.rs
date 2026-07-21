//! F4 algorithm for computing Gröbner bases.
//!
//! Implements the matrix-based F4 algorithm from Faugère's 1999 paper.
//!
//! The key idea: replace sequential S-polynomial reductions with batched
//! sparse-matrix row operations (Gaussian elimination over the coefficient field).

use smallvec::SmallVec;

use ocas_core::FastHashMap as HashMap;
use ocas_core::FastHashSet as HashSet;
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

/// Minimal interface the F4 driver needs from a basis polynomial.
///
/// Implemented by both [`SparseMultivariatePolynomial`] (generic path)
/// and [`FpPoly`] (native ℤ_p fast path), so pair management and the
/// simplification cache are shared between the two pipelines.
trait BasisPoly: Clone {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>>;
    fn n_vars(&self) -> usize;
    fn n_terms(&self) -> usize;
    fn mul_monomial(&self, exp: &[usize]) -> Self;
}

impl<D: Domain, O: MonomialOrder> BasisPoly for SparseMultivariatePolynomial<D, O> {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        SparseMultivariatePolynomial::leading_monomial(self)
    }
    fn n_vars(&self) -> usize {
        SparseMultivariatePolynomial::n_vars(self)
    }
    fn n_terms(&self) -> usize {
        SparseMultivariatePolynomial::n_terms(self)
    }
    fn mul_monomial(&self, exp: &[usize]) -> Self {
        SparseMultivariatePolynomial::mul_monomial(self, exp)
    }
}

impl CriticalPair {
    fn new<P: BasisPoly>(basis: &[P], i: usize, j: usize) -> Option<Self> {
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
type SimpCache<P> = Vec<(SmallVec<[usize; 4]>, P)>;

/// Tracks a monomial's state during symbolic preprocessing.
#[derive(Debug, Clone)]
struct MonomialData {
    /// Column index in the matrix (assigned during column construction).
    column: usize,
}

// =========================================================================
//  Leading-monomial divisor index
// =========================================================================

/// Support bitmask of a monomial: bit `v` is set iff `exp[v] > 0`.
///
/// Variables beyond index 63 are ignored in the mask; this only weakens
/// the filter's selectivity, never its correctness (a false-positive
/// bucket match is filtered out by the exact divisibility check).
fn support_mask(exp: &[usize]) -> u64 {
    let mut mask = 0u64;
    for (v, &e) in exp.iter().enumerate() {
        if e > 0 && v < 64 {
            mask |= 1 << v;
        }
    }
    mask
}

/// Incremental index over basis leading monomials for fast reducer
/// (divisor) queries, replacing the O(monomials × basis) linear scan.
///
/// Basis elements are bucketed by the exact support mask of their leading
/// monomial. A query for `exp` enumerates the submasks of `support(exp)`
/// (a basis LM dividing `exp` must have its support inside `exp`'s) and
/// checks exact divisibility only within those buckets.
struct DivisorIndex {
    /// Exact support mask of a basis LM → basis indices with that mask.
    buckets: HashMap<u64, Vec<usize>>,
}

impl DivisorIndex {
    fn new() -> Self {
        Self {
            buckets: HashMap::default(),
        }
    }

    /// Register the leading monomial of basis element `idx`.
    fn push(&mut self, lm: &[usize], idx: usize) {
        self.buckets.entry(support_mask(lm)).or_default().push(idx);
    }
}

/// Find a reducer for `exp` in the basis: a basis element whose leading
/// monomial divides `exp`, chosen with the smallest number of terms
/// (ties keep the lowest basis index, matching the old linear scan).
fn find_reducer<P: BasisPoly>(index: &DivisorIndex, basis: &[P], exp: &[usize]) -> Option<usize> {
    let mask = support_mask(exp);
    let mut best: Option<usize> = None;
    // Enumerate all submasks of `mask`, including `mask` itself and 0.
    let mut sub = mask;
    loop {
        if let Some(ids) = index.buckets.get(&sub) {
            for &bi in ids {
                if let Some(blm) = basis[bi].leading_monomial()
                    && monomial_divides(exp, blm)
                {
                    match best {
                        Some(b) if basis[b].n_terms() <= basis[bi].n_terms() => {}
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
    best
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

    // ℤ_p fast path: run the entire F4 pipeline on native i64 residues,
    // converting to/from the BigInt-backed domain only at the boundaries.
    if std::any::TypeId::of::<D>() == std::any::TypeId::of::<FiniteField>() {
        let domain_ptr = ideal[0].domain() as *const D;
        let ff = unsafe { &*domain_ptr.cast::<FiniteField>() };
        return f4_fp(ideal, ff.prime_u64() as i64);
    }

    // Filter zeros and make monic.
    let mut initial: Vec<SparseMultivariatePolynomial<D, O>> =
        ideal.iter().filter(|p| !p.is_zero()).cloned().collect();
    for p in &mut initial {
        make_monic(p);
    }
    if initial.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Feed generators one at a time so that initial critical pairs go
    // through the same Gebauer–Moeller filtering (as Symbolica does).
    let mut basis: Vec<SparseMultivariatePolynomial<D, O>> = Vec::new();
    let mut pairs: Vec<CriticalPair> = Vec::new();
    let mut simplifications: Vec<SimpCache<SparseMultivariatePolynomial<D, O>>> = Vec::new();
    // Reducer divisor index and LM set, maintained in lockstep with `basis`.
    let mut div_index = DivisorIndex::new();
    let mut basis_lm_set: HashSet<SmallVec<[usize; 4]>> = HashSet::default();
    for p in initial {
        update_pairs(&mut basis, &mut pairs, &mut simplifications, p);
        let idx = basis.len() - 1;
        if let Some(lm) = basis[idx].leading_monomial() {
            div_index.push(lm, idx);
            basis_lm_set.insert(lm.clone());
        }
    }

    // Reusable buffers.
    // MonomialData tracks the column index assigned to each monomial;
    // newly registered monomials are pushed onto `worklist` and processed
    // exactly once by symbolic preprocessing.
    let mut all_monomials: HashMap<SmallVec<[usize; 4]>, MonomialData> = HashMap::default();
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut matrix: Vec<Vec<(D::Element, usize)>> = Vec::new();
    let mut pivots: Vec<Option<usize>> = Vec::new();
    // Head monomials of every input row of the current matrix (basis
    // multiples). The extraction step adds a reduced row to the basis
    // only when its leading monomial is NOT in this set.
    let mut input_heads: HashSet<SmallVec<[usize; 4]>> = HashSet::default();
    // Deduplication of (basis index, exponent diff) pairs within one
    // matrix construction.
    let mut seen_rows: HashSet<(usize, SmallVec<[usize; 4]>)> = HashSet::default();
    // Monomials awaiting a reducer search (filled by add_poly_to_matrix).
    let mut worklist: Vec<SmallVec<[usize; 4]>> = Vec::new();

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
        input_heads.clear();
        seen_rows.clear();
        worklist.clear();

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

            // Classic F4: add BOTH multiples as separate rows (never the
            // precomputed difference). Every input row is then a basis
            // multiple, so its head monomial lies in the basis LM ideal —
            // the invariant the extraction criterion relies on.
            // Reference: Symbolica groebner.rs, pair selection loop.
            input_heads.insert(lcm_exp.clone());
            for (idx, diff) in [(i, &diff_i), (j, &diff_j)] {
                if seen_rows.insert((idx, diff.clone())) {
                    let mult = get_simplified(&simplifications[idx], diff, &basis[idx]);
                    add_poly_to_matrix(
                        &mult,
                        &mut matrix,
                        &mut all_monomials,
                        &mut monomial_list,
                        &mut worklist,
                    );
                }
            }
        }

        if matrix.is_empty() {
            continue;
        }

        // --- Iterative symbolic preprocessing (Faugère F4 key innovation) ---
        // Every monomial registered in the matrix is pushed onto the
        // worklist exactly once. Pop a monomial, search for a reducer in
        // the basis via the LM divisor index, and append it as a new row
        // (whose own monomials may in turn extend the worklist). Repeat
        // until no unprocessed monomials remain. This ensures the matrix
        // contains all necessary reduction information for the echelon
        // step.
        //
        // Reference: Symbolica groebner.rs L262-288.
        while let Some(exp) = worklist.pop() {
            if let Some(bi) = find_reducer(&div_index, &basis, &exp) {
                let blm = basis[bi].leading_monomial().unwrap();
                let diff: SmallVec<[usize; 4]> =
                    exp.iter().zip(blm.iter()).map(|(a, b)| a - b).collect();
                let reducer = get_simplified(&simplifications[bi], &diff, &basis[bi]);
                // The reducer's head is exactly `exp` — record it as an
                // input head for the extraction criterion.
                input_heads.insert(exp.clone());
                // Add reducer row to matrix, registering any new monomials.
                add_poly_to_matrix(
                    &reducer,
                    &mut matrix,
                    &mut all_monomials,
                    &mut monomial_list,
                    &mut worklist,
                );
            }
        }

        if matrix.is_empty() || monomial_list.is_empty() {
            continue;
        }

        let ncols = monomial_list.len();

        // --- Sort columns: DESCENDING monomial order ---
        //
        // Column 0 must be the leading (greatest) monomial: rows store the
        // leading term first, `sort_rows` orders rows by first column, and
        // the elimination scan processes columns 0..ncols. Sorting columns
        // ascending instead would put pivots on TRAILING terms and break
        // the entire echelon step (this was the root cause of the
        // extraction blowup: the echelon form was decorative and all real
        // work fell back to polynomial division).
        let mut col_order: Vec<usize> = (0..ncols).collect();
        col_order.sort_unstable_by(|&a, &b| O::cmp(&monomial_list[b], &monomial_list[a]));

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
        echelonize_generic(&mut matrix, ncols, basis[0].domain(), &mut pivots);

        // --- Extract new polynomials from reduced rows ---
        //
        // Symbolica/Faugère criterion: a reduced row joins the basis only
        // when its leading monomial differs from every input row's head.
        // Every input row is a basis multiple, so "row LM == some input
        // head" means the LM already lies in the basis LM ideal and the
        // row carries no new information. Crucially, NO reduction of the
        // extracted polynomial against the basis is needed — rows of the
        // echelonized matrix are already fully inter-reduced, and tail
        // reduction is deferred to the final auto_reduce pass. This
        // eliminates the dominant cost of the old pipeline (extraction was
        // 99.98% of runtime on cyclic-5).
        for row in &matrix {
            if row.is_empty() {
                continue;
            }
            let row_lm = &sorted_monomials[row[0].1];
            if input_heads.contains(row_lm) {
                continue;
            }
            let mut poly = basis[0].zero();
            for (coeff, col) in row.iter().rev() {
                poly.append_monomial(coeff.clone(), &sorted_monomials[*col]);
            }
            if poly.is_zero() {
                continue;
            }

            let new_lm = poly.leading_monomial().unwrap().clone();

            // Skip if a polynomial with this leading monomial already exists.
            if basis_lm_set.contains(&new_lm) {
                continue;
            }

            // Add to basis with Gebauer-Moeller pair filtering.
            update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
            let idx = basis.len() - 1;
            if let Some(lm) = basis[idx].leading_monomial() {
                div_index.push(lm, idx);
                basis_lm_set.insert(lm.clone());
            }
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
fn update_pairs<P: BasisPoly>(
    basis: &mut Vec<P>,
    pairs: &mut Vec<CriticalPair>,
    simplifications: &mut Vec<SimpCache<P>>,
    new_poly: P,
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

    let is_disjoint = |cp: &CriticalPair| {
        let a = basis[cp.idx1].leading_monomial().unwrap();
        let b = basis[cp.idx2].leading_monomial().unwrap();
        a.iter().zip(b.iter()).all(|(x, y)| *x == 0 || *y == 0)
    };

    // Generate ALL new pairs. Coprime pairs participate in the dominance
    // computation below (their lcm can dominate other pairs) but are
    // never kept themselves (product criterion).
    let mut new_pairs: Vec<(CriticalPair, bool)> = (0..new_idx)
        .filter_map(|i| CriticalPair::new(basis, i, new_idx))
        .map(|cp| {
            let disjoint = is_disjoint(&cp);
            (cp, disjoint)
        })
        .collect();

    // Gebauer–Moeller chain criterion (sequential elimination): pair i is
    // redundant when another still-live pair j has lcm(j) dividing lcm(i)
    // (equality included — the later duplicate is dropped). Self-exclusion
    // works by clearing the flag before evaluation.
    // Reference: Symbolica groebner.rs `update`; Becker–Weispfenning.
    for i in 0..new_pairs.len() {
        new_pairs[i].1 = false;
        let disjoint = is_disjoint(&new_pairs[i].0);
        let survive = disjoint
            || new_pairs.iter().all(|p2| {
                !p2.1
                    || new_pairs[i]
                        .0
                        .lcm
                        .iter()
                        .zip(p2.0.lcm.iter())
                        .any(|(a, b)| a < b)
            });
        new_pairs[i].1 = survive;
    }
    let kept: Vec<CriticalPair> = new_pairs
        .into_iter()
        .filter(|(cp, k)| *k && !is_disjoint(cp))
        .map(|(cp, _)| cp)
        .collect();

    // Gebauer–Moeller update criterion for existing pairs: drop {i,j} only
    // when new_lm divides lcm(i,j) AND both lcm(i,new) and lcm(j,new) are
    // strictly smaller than lcm(i,j). The earlier implementation checked
    // only strict divisibility, which incorrectly removed pairs whose
    // lcm is reproduced by (i, new) — dropping S-polynomials that are
    // required for completeness (cyclic-5 failed is_groebner_basis).
    pairs.retain(|cp| {
        let new_divides = cp.lcm.iter().zip(new_lm.iter()).all(|(a, b)| a >= b);
        if !new_divides {
            return true;
        }
        let lm1 = basis[cp.idx1].leading_monomial().unwrap();
        let lm2 = basis[cp.idx2].leading_monomial().unwrap();
        let same1 = lm1
            .iter()
            .zip(new_lm.iter())
            .zip(cp.lcm.iter())
            .all(|((a, b), c)| (*a).max(*b) == *c);
        let same2 = lm2
            .iter()
            .zip(new_lm.iter())
            .zip(cp.lcm.iter())
            .all(|((a, b), c)| (*a).max(*b) == *c);
        same1 || same2
    });

    pairs.extend(kept);

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
    worklist: &mut Vec<SmallVec<[usize; 4]>>,
) {
    let mut row: Vec<(D::Element, usize)> = Vec::new();
    for (exp, coeff) in poly.sorted_terms().iter().rev() {
        if poly.domain().is_zero(coeff) {
            continue;
        }
        let mut new_col = None;
        let md = monomial_map.entry((*exp).clone()).or_insert_with(|| {
            let idx = monomial_list.len();
            monomial_list.push((*exp).clone());
            new_col = Some(idx);
            MonomialData { column: idx }
        });
        if new_col.is_some() {
            worklist.push((*exp).clone());
        }
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
fn get_simplified<P: BasisPoly>(cache: &SimpCache<P>, diff: &[usize], basis_poly: &P) -> P {
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
//  Native ℤ_p pipeline (FpPoly)
// =========================================================================

/// Reduce `a` into `[0, p)`.
#[inline]
fn norm_mod(a: i64, p: i64) -> i64 {
    let r = a % p;
    if r < 0 { r + p } else { r }
}

/// A native ℤ_p polynomial for the F4 fast path.
///
/// Terms are stored as a `Vec` sorted by **descending** monomial order
/// (leading term first) with coefficients in `[0, p)`. All arithmetic is
/// `i64` modular arithmetic — no `BigInt` is touched inside the pipeline.
/// Requires `p < 2^31` so that products of two residues fit in `i64`
/// (the same assumption `echelonize_fp` already makes with its `p²`
/// slack).
#[derive(Debug, Clone)]
struct FpPoly {
    /// Terms sorted descending by monomial order; coeffs in `[0, p)`.
    terms: Vec<(SmallVec<[usize; 4]>, i64)>,
    n_vars: usize,
}

impl FpPoly {
    fn zero(n_vars: usize) -> Self {
        Self {
            terms: Vec::new(),
            n_vars,
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn n_terms(&self) -> usize {
        self.terms.len()
    }

    /// Convert a domain polynomial to native residues (one-time cost at
    /// the pipeline boundary).
    fn from_domain<D: Domain + 'static, O: MonomialOrder>(
        p: &SparseMultivariatePolynomial<D, O>,
        prime: i64,
    ) -> Self {
        let mut terms: Vec<(SmallVec<[usize; 4]>, i64)> = Vec::with_capacity(p.n_terms());
        // sorted_terms ascending → rev gives descending (leading first).
        for (exp, coeff) in p.sorted_terms().iter().rev() {
            let c = norm_mod(domain_to_i64_fp::<D>(coeff, prime), prime);
            if c != 0 {
                terms.push(((*exp).clone(), c));
            }
        }
        Self {
            terms,
            n_vars: p.n_vars(),
        }
    }

    /// Convert back to the BigInt-backed domain (one-time cost per
    /// surviving basis element and at the end of the pipeline).
    fn to_domain<D: Domain + 'static, O: MonomialOrder>(
        &self,
        domain: &D,
        prime: i64,
    ) -> SparseMultivariatePolynomial<D, O> {
        let mut poly = SparseMultivariatePolynomial::new(domain.clone(), self.n_vars);
        for (exp, c) in &self.terms {
            poly.append_monomial(i64_to_domain_fp::<D>(domain, *c, prime), exp);
        }
        poly
    }
}

impl BasisPoly for FpPoly {
    fn leading_monomial(&self) -> Option<&SmallVec<[usize; 4]>> {
        self.terms.first().map(|t| &t.0)
    }
    fn n_vars(&self) -> usize {
        self.n_vars
    }
    fn n_terms(&self) -> usize {
        FpPoly::n_terms(self)
    }
    fn mul_monomial(&self, exp: &[usize]) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(e, c)| {
                    (
                        e.iter()
                            .zip(exp.iter())
                            .map(|(a, b)| a + b)
                            .collect::<SmallVec<[usize; 4]>>(),
                        *c,
                    )
                })
                .collect(),
            n_vars: self.n_vars,
        }
    }
}

/// Make an `FpPoly` monic (scale by the modular inverse of the leading
/// coefficient).
fn monic_fp(p: &mut FpPoly, prime: i64) {
    if let Some(&(_, lc)) = p.terms.first()
        && lc != 1
    {
        let inv = mod_inv(lc, prime);
        for (_, c) in &mut p.terms {
            *c = norm_mod(*c * inv, prime);
        }
    }
}

/// Fetch (or build and cache) the polynomial content of the basis
/// multiple `basis[basis_idx] * x^diff`. The row content is cached across
/// rounds so a repeated reducer/S-multiple is built only once; monomial
/// registration into the per-round column table happens at materialize
/// time. Returns the row-store index, or `None` for a zero polynomial.
fn register_row_fp(
    basis_idx: usize,
    diff: &SmallVec<[usize; 4]>,
    basis: &[FpPoly],
    simplifications: &[SimpCache<FpPoly>],
    row_store: &mut Vec<FpPoly>,
    row_cache: &mut HashMap<(usize, SmallVec<[usize; 4]>), usize>,
) -> Option<usize> {
    let key = (basis_idx, diff.clone());
    if let Some(&rs) = row_cache.get(&key) {
        return Some(rs);
    }
    let poly = get_simplified(&simplifications[basis_idx], diff, &basis[basis_idx]);
    if poly.is_zero() {
        return None;
    }
    let rs = row_store.len();
    row_store.push(poly);
    row_cache.insert(key, rs);
    Some(rs)
}

/// Native ℤ_p F4: the full F4 pipeline on `i64` residues.
///
/// Structurally identical to the generic [`f4`] loop, but every
/// polynomial operation (S-polynomials, symbolic preprocessing, row
/// echelon, tail reduction) runs on [`FpPoly`] — `BigInt` conversions
/// happen only when reading the input and emitting surviving basis
/// elements.
#[allow(clippy::too_many_lines)]
fn f4_fp<D: Domain + 'static, O: MonomialOrder>(
    ideal: &[SparseMultivariatePolynomial<D, O>],
    prime: i64,
) -> GroebnerBasis<D, O> {
    let n_vars = ideal[0].n_vars();

    // Filter zeros, convert to native residues, and make monic.
    let mut initial: Vec<FpPoly> = ideal
        .iter()
        .filter(|p| !p.is_zero())
        .map(|p| FpPoly::from_domain(p, prime))
        .collect();
    for p in &mut initial {
        monic_fp(p, prime);
    }
    if initial.is_empty() {
        return GroebnerBasis { basis: vec![] };
    }

    // Feed generators one at a time (same GM filtering as the generic path).
    let mut basis: Vec<FpPoly> = Vec::new();
    let mut pairs: Vec<CriticalPair> = Vec::new();
    let mut simplifications: Vec<SimpCache<FpPoly>> = Vec::new();
    // Reducer divisor index and LM set, maintained in lockstep with `basis`.
    let mut div_index = DivisorIndex::new();
    let mut basis_lm_set: HashSet<SmallVec<[usize; 4]>> = HashSet::default();
    for p in initial {
        update_pairs(&mut basis, &mut pairs, &mut simplifications, p);
        let idx = basis.len() - 1;
        if let Some(lm) = basis[idx].leading_monomial() {
            div_index.push(lm, idx);
            basis_lm_set.insert(lm.clone());
        }
    }

    // Monomial table: rebuilt each round (clearing keeps the column space
    // per-round compact — a persistent table accumulated every historical
    // monomial and exploded the echelon's pivots array and column sort).
    let mut all_monomials: HashMap<SmallVec<[usize; 4]>, MonomialData> = HashMap::default();
    let mut monomial_list: Vec<SmallVec<[usize; 4]>> = Vec::new();
    let mut matrix: Vec<Vec<(i64, usize)>> = Vec::new();
    let mut pivots: Vec<Option<usize>> = Vec::new();
    // Head monomials of every input row of the current matrix (basis
    // multiples). See the generic path for the extraction invariant.
    let mut input_heads: HashSet<SmallVec<[usize; 4]>> = HashSet::default();
    let mut seen_rows: HashSet<(usize, SmallVec<[usize; 4]>)> = HashSet::default();
    // Monomials awaiting a reducer search this round.
    let mut worklist: Vec<SmallVec<[usize; 4]>> = Vec::new();
    // Per-round dedup for the worklist monomial scan.
    let mut seen_monomials: HashSet<SmallVec<[usize; 4]>> = HashSet::default();

    // Row template cache (persistent across rounds): a basis multiple is
    // built once and its polynomial content stored; subsequent rounds
    // reuse the stored content instead of re-running the simplification
    // cache + polynomial walk. Basis elements are only ever appended, so
    // `(basis_idx, diff)` keys stay valid.
    let mut row_store: Vec<FpPoly> = Vec::new();
    let mut row_cache: HashMap<(usize, SmallVec<[usize; 4]>), usize> = HashMap::default();
    // Row-store indices participating in the current round's matrix.
    let mut round_rows: Vec<usize> = Vec::new();

    // Optional section timing for performance diagnosis (OCAS_F4_STATS=1).
    let stats = std::env::var("OCAS_F4_STATS").is_ok();
    let mut rounds = 0usize;
    let mut added = 0usize;
    let mut t_build = std::time::Duration::ZERO;
    let mut t_pre = std::time::Duration::ZERO;
    let mut t_ech = std::time::Duration::ZERO;
    let mut t_ext = std::time::Duration::ZERO;
    // 每轮矩阵规模直方（OCAS_F4_ROUND_STATS=1 逐行打印）。
    let round_stats = std::env::var("OCAS_F4_ROUND_STATS").is_ok();

    while !pairs.is_empty() {
        rounds += 1;
        let t0 = std::time::Instant::now();
        // --- Selection: find minimum lcm degree ---
        let min_deg = pairs.iter().map(|cp| cp.degree).min().unwrap();

        let selected: Vec<CriticalPair> = pairs
            .iter()
            .filter(|cp| cp.degree == min_deg)
            .cloned()
            .collect();

        let sel_set: std::collections::HashSet<(usize, usize)> =
            selected.iter().map(|cp| (cp.idx1, cp.idx2)).collect();
        pairs.retain(|cp| !sel_set.contains(&(cp.idx1, cp.idx2)));

        if selected.is_empty() {
            continue;
        }

        // --- Build matrix rows from selected pairs ---
        // Per-round column table; row content comes from the persistent
        // row cache.
        all_monomials.clear();
        monomial_list.clear();
        matrix.clear();
        input_heads.clear();
        seen_rows.clear();
        worklist.clear();
        seen_monomials.clear();
        round_rows.clear();

        for cp in &selected {
            let i = cp.idx1;
            let j = cp.idx2;
            let lm_i = basis[i].leading_monomial().unwrap();
            let lm_j = basis[j].leading_monomial().unwrap();
            let lcm_exp = &cp.lcm;

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

            // Classic F4: add BOTH multiples as separate rows (see the
            // generic path for why the difference is never precomputed).
            input_heads.insert(lcm_exp.clone());
            for (idx, diff) in [(i, &diff_i), (j, &diff_j)] {
                if seen_rows.insert((idx, diff.clone()))
                    && let Some(rs) = register_row_fp(
                        idx,
                        diff,
                        &basis,
                        &simplifications,
                        &mut row_store,
                        &mut row_cache,
                    )
                {
                    round_rows.push(rs);
                }
            }
        }

        if round_rows.is_empty() {
            continue;
        }

        t_build += t0.elapsed();
        let t1 = std::time::Instant::now();

        // --- Iterative symbolic preprocessing ---
        // Collect every monomial appearing in this round's rows into the
        // worklist (deduplicated per round), then process LIFO: search a
        // reducer via the LM divisor index and append it as a new row
        // (whose own monomials extend the worklist). The row cache keeps
        // repeated reducers free; the per-round worklist scan still
        // reduces every monomial against the CURRENT basis.
        for &rs in &round_rows {
            for (exp, _) in &row_store[rs].terms {
                if seen_monomials.insert(exp.clone()) {
                    worklist.push(exp.clone());
                }
            }
        }
        while let Some(exp) = worklist.pop() {
            if let Some(bi) = find_reducer(&div_index, &basis, &exp) {
                let blm = basis[bi].leading_monomial().unwrap();
                let diff: SmallVec<[usize; 4]> =
                    exp.iter().zip(blm.iter()).map(|(a, b)| a - b).collect();
                input_heads.insert(exp);
                if let Some(rs) = register_row_fp(
                    bi,
                    &diff,
                    &basis,
                    &simplifications,
                    &mut row_store,
                    &mut row_cache,
                ) {
                    round_rows.push(rs);
                    for (mexp, _) in &row_store[rs].terms {
                        if seen_monomials.insert(mexp.clone()) {
                            worklist.push(mexp.clone());
                        }
                    }
                }
            }
        }

        if round_rows.is_empty() {
            continue;
        }

        t_pre += t1.elapsed();
        let t2 = std::time::Instant::now();

        // Materialize this round's matrix: register each stored row's
        // monomials into the per-round column table (in stored order —
        // rows store the leading term first, so column 0 of each row is
        // its head before the global sort).
        for &rs in &round_rows {
            let poly = &row_store[rs];
            let mut row: Vec<(i64, usize)> = Vec::with_capacity(poly.terms.len());
            for (exp, coeff) in &poly.terms {
                let col = all_monomials.entry(exp.clone()).or_insert_with(|| {
                    let idx = monomial_list.len();
                    monomial_list.push(exp.clone());
                    MonomialData { column: idx }
                });
                row.push((*coeff, col.column));
            }
            if !row.is_empty() {
                matrix.push(row);
            }
        }

        let ncols = monomial_list.len();

        // --- Sort columns: DESCENDING monomial order ---
        // (See the generic path for why descending order is essential.)
        let mut col_order: Vec<usize> = (0..ncols).collect();
        col_order.sort_unstable_by(|&a, &b| O::cmp(&monomial_list[b], &monomial_list[a]));

        let mut col_inv = vec![0usize; ncols];
        for (new_col, &old_col) in col_order.iter().enumerate() {
            col_inv[old_col] = new_col;
        }

        for row in &mut matrix {
            for (_, col) in row.iter_mut() {
                *col = col_inv[*col];
            }
        }

        let mut sorted_monomials: Vec<SmallVec<[usize; 4]>> = vec![SmallVec::new(); ncols];
        for (new_col, &old_col) in col_order.iter().enumerate() {
            sorted_monomials[new_col] = monomial_list[old_col].clone();
        }

        // --- Row echelon form (native i64) ---
        echelonize_fp(&mut matrix, ncols, prime, &mut pivots);

        if round_stats {
            let nnz: usize = matrix.iter().map(Vec::len).sum();
            let max_len = matrix.iter().map(Vec::len).max().unwrap_or(0);
            eprintln!(
                "  round {rounds}: rows={} cols={ncols} nnz={nnz} maxlen={max_len} sel={}",
                matrix.len(),
                selected.len()
            );
        }

        t_ech += t2.elapsed();
        let t3 = std::time::Instant::now();

        // --- Extract new polynomials from reduced rows ---
        //
        // Same Symbolica/Faugère criterion as the generic path: add a row
        // only when its leading monomial differs from every input row's
        // head; NO reduction against the basis is performed here.
        for row in &matrix {
            if row.is_empty() {
                continue;
            }
            let row_lm = &sorted_monomials[row[0].1];
            if input_heads.contains(row_lm) {
                continue;
            }
            let mut poly = FpPoly::zero(n_vars);
            for &(c, col) in row {
                let v = norm_mod(c, prime);
                if v != 0 {
                    poly.terms.push((sorted_monomials[col].clone(), v));
                }
            }
            if poly.is_zero() {
                continue;
            }

            let new_lm = poly.leading_monomial().unwrap().clone();

            // Skip if a polynomial with this leading monomial already exists.
            if basis_lm_set.contains(&new_lm) {
                continue;
            }

            update_pairs(&mut basis, &mut pairs, &mut simplifications, poly);
            let idx = basis.len() - 1;
            if let Some(lm) = basis[idx].leading_monomial() {
                div_index.push(lm, idx);
                basis_lm_set.insert(lm.clone());
            }
            added += 1;
        }

        t_ext += t3.elapsed();
    }

    if stats {
        eprintln!(
            "f4_fp stats: rounds={rounds} added={added} | build={t_build:.2?} pre={t_pre:.2?} echelon={t_ech:.2?} extract={t_ext:.2?}"
        );
    }

    // Convert back to the domain representation and post-process with the
    // shared minimize/inter-reduce pipeline.
    let domain = ideal[0].domain().clone();
    let basis_d: Vec<SparseMultivariatePolynomial<D, O>> = basis
        .iter()
        .map(|p| p.to_domain::<D, O>(&domain, prime))
        .collect();
    GroebnerBasis { basis: basis_d }.minimize().auto_reduce()
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
) {
    let p = prime;

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

    // Scratch row reused by every merge step (avoids per-subtraction
    // allocation churn).
    let mut scratch: Vec<(i64, usize)> = Vec::new();

    // Reduce rows. Each row repeatedly cancels its head against the pivot
    // row registered at the head column (sparse two-pointer merge), until
    // the row vanishes or its head column has no pivot yet — in which case
    // the row is normalized to monic and registered as the new pivot.
    // This is semantically identical to the old dense-buffer scan but
    // costs O(nnz) per cancellation instead of O(ncols) per row.
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }
        if pivots[matrix[r][0].1] == Some(r) {
            continue;
        }

        let mut row = std::mem::take(&mut matrix[r]);
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(pi) if pi != r => {
                    // Pivot rows are monic, so the heads cancel exactly.
                    let c = row[0].0;
                    sub_scaled_fp(&mut row, &matrix[pi], c, p, &mut scratch);
                }
                Some(_) => break,
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
        matrix[r] = row;
    }

    matrix.retain(|r| !r.is_empty());
}

/// Sparse `row -= c * pivot` (mod `p`) by merging the two
/// column-ascending rows. The head columns coincide and cancel (pivot is
/// monic), so both heads are skipped. All emitted coefficients are
/// normalized into `[0, p)`; zeros are dropped.
fn sub_scaled_fp(
    row: &mut Vec<(i64, usize)>,
    pivot: &[(i64, usize)],
    c: i64,
    p: i64,
    scratch: &mut Vec<(i64, usize)>,
) {
    scratch.clear();
    scratch.reserve(row.len() + pivot.len());
    let mut i = 1;
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
//  Row echelon form: generic domain path
// =========================================================================

#[allow(clippy::needless_range_loop)]
fn echelonize_generic<D: Domain>(
    matrix: &mut Vec<Vec<(D::Element, usize)>>,
    ncols: usize,
    domain: &D,
    pivots: &mut Vec<Option<usize>>,
) {
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
            let lc = matrix[r][0].0.clone();
            if !domain.is_one(&lc)
                && let Some(inv) = domain.inv(&lc)
            {
                for (c, _) in &mut matrix[r] {
                    *c = domain.mul(c, &inv);
                }
            }
        }
    }

    let mut scratch: Vec<(D::Element, usize)> = Vec::new();

    // Reduce rows with sparse merge cancellation (same structure as the
    // ℤ_p fast path — see `echelonize_fp`).
    for r in 0..matrix.len() {
        if matrix[r].is_empty() {
            continue;
        }
        if pivots[matrix[r][0].1] == Some(r) {
            continue;
        }

        let mut row = std::mem::take(&mut matrix[r]);
        loop {
            if row.is_empty() {
                break;
            }
            let head_col = row[0].1;
            match pivots[head_col] {
                Some(pi) if pi != r => {
                    let c = row[0].0.clone();
                    sub_scaled_generic(domain, &mut row, &matrix[pi], &c, &mut scratch);
                }
                Some(_) => break,
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
        matrix[r] = row;
    }

    matrix.retain(|r| !r.is_empty());
}

/// Sparse `row -= c * pivot` over a generic domain by merging the two
/// column-ascending rows. The head columns coincide and cancel (pivot is
/// monic), so both heads are skipped; zero coefficients are dropped.
fn sub_scaled_generic<D: Domain>(
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
    fn f4_cyclic_3_fp13_matches_q() {
        let field = FiniteField::new(BigInt::from(13u32));
        let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            3,
            vec![
                (vec![1, 0, 0], field.element(1)),
                (vec![0, 1, 0], field.element(1)),
                (vec![0, 0, 1], field.element(1)),
            ],
        );
        let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            3,
            vec![
                (vec![1, 1, 0], field.element(1)),
                (vec![0, 1, 1], field.element(1)),
                (vec![1, 0, 1], field.element(1)),
            ],
        );
        let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
            field.clone(),
            3,
            vec![
                (vec![1, 1, 1], field.element(1)),
                (vec![0, 0, 0], field.element(12)),
            ],
        );
        let gb = f4(&[f1, f2, f3]);
        assert!(!gb.basis.is_empty());
        assert!(gb.is_groebner_basis());
        // cyclic-3 (zero-dim, degree 6) has a reduced Gröbner basis with
        // exactly 3 elements over any field of characteristic ≠ 2, 3.
        // This is the definitive regression test for the extraction LM
        // pre-skip: if it collapses the ideal, the basis shrinks.
        assert_eq!(gb.basis.len(), 3);
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
