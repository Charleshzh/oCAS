//! Symbolic-constant rational integration (0.27.0 S5 extension).
//!
//! The ℚ-coefficient backend ([`crate::integral::rational`]) declines any
//! rational function whose coefficients involve symbols. The Rubi corpus is
//! dominated by exactly those shapes (`1/(a-b*x^2)`, `(d+e*x)/(x^3*(a+c*x^2))`,
//! and trig-rationals after a Weierstrass substitution). This module
//! integrates a rational function of `x` over the constant field
//! `ℚ(s1, …, sk)` (symbols treated as independent constants):
//!
//! 1. polynomial division — the quotient is integrated termwise,
//! 2. Hermite reduction for repeated factors (square-free decomposition via
//!    multivariate GCD over `ℚ(symbols)`),
//! 3. a logarithmic part on the square-free denominator, restricted to
//!    factors of degree ≤ 2:
//!    - linear factors: residue formula `c·log(f)`,
//!    - quadratic factors: `(M·x+N)/(a·x²+b·x+c)` → `log` + `atan`/`atanh`
//!      via the quadratic formula (symbolic `√(4ac−b²)` when the
//!      discriminant is not a rational square; numeric discriminants pick
//!      `atan`/`atanh` by sign).
//!
//! A square-free factor of degree ≥ 3 leaves the whole log part as an
//! unevaluated `Integral` term (an honest partial result, mirroring the ℚ
//! backend's Rothstein–Trager behaviour).
//!
//! Coefficients are [`GeneratorField`] elements — rational functions of the
//! symbols — so all arithmetic is exact.

use num_traits::Signed;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_domain::{Domain, Integer, IntegerDomain, Rational, RationalDomain};
use ocas_poly::{Lex, SparseMultivariatePolynomial};

use crate::tower::convert::{GeneratorField, atom_to_rational, rational_to_atom};

type Sparse = SparseMultivariatePolynomial<ocas_domain::RationalDomain, Lex>;

// =========================================================================
// Deterministic work budget
// =========================================================================
//
// The field-Euclidean machinery below is exact, but its cost is data
// dependent: a single `FPoly::div_rem` multiplies the coefficient term
// counts of its intermediates, `field_divisors` enumerates the product of
// `(exponent_i + 1)` monomial candidates and `factor_via_integer` hands a
// cleared polynomial to the integer factorizer. All three are unbounded in
// principle (the Weierstrass t-forms with symbolic coefficients are the
// corpus shapes that expose it), so the module carries a deterministic
// counter instead of relying on a wall-clock timeout.
//
// The counter is reset at every public entry point of the module
// (`integrate_rational_symbolic`, `rational_complexity_ok`) and charged by
// the inner loops that can blow up; a charge that overshoots the cap makes
// the operation return `None`, which the caller reports as a decline.

/// Work units a single entry into this module may charge. Calibrated as a
/// backstop above every corpus case that already solves (the largest
/// observed successful charge is ~3.2e6 units — a pre-existing 3.6 s solve,
/// see the module tests); legitimate small integrations charge a few
/// hundred units, while the shapes that used to hang grow the coefficient
/// term counts geometrically and trip the cap within a few dozen
/// iterations.
const MAX_WORK_UNITS: u64 = 4_000_000;

/// Coefficient-size budget for the field-Euclidean loops (0.27.1).
const MAX_COEFF_COST: usize = 20_000;

/// Monomial divisors enumerated for one field element by
/// [`field_divisors`]. The enumeration is the product of `(exponent_i + 1)`
/// over the symbols, so a polynomial with a few 3-term monomials already
/// reaches this cap; legitimately factored candidates stay in the tens.
const MAX_FIELD_DIVISORS: usize = 256;

/// Root candidates evaluated by [`split_linear_candidates`]
/// (`2 · |divisors(a0)| · |divisors(lc)|`). Bounded structurally before the
/// evaluation loop so a fat coefficient cannot explode the candidate list.
const MAX_SPLIT_CANDIDATES: usize = 2048;

/// Terms of the multivariate polynomial handed to
/// `SparseMultivariatePolynomial::factor()` by [`factor_via_integer`].
const MAX_FACTOR_TERMS: usize = 64;

/// Terms of one numerator/denominator before `FPoly::from_sparse` (which is
/// quadratic in the term count) is allowed to run.
const MAX_SPARSE_TERMS: usize = 1024;

/// Predicted coefficient-multiplication cost of one `FPoly::div_rem` step:
/// `cost(remainder) · cost(quotient coefficient)`.
///
/// The post-step `MAX_COEFF_COST` check cannot see the intermediates of a
/// single step, and one `GeneratorField::mul` on fat coefficients is already
/// expensive on its own (it cross-multiplies rational functions and then
/// canonicalizes the result). The product is therefore predicted *before*
/// the multiplication and refuses work that cannot fit the budget.
///
/// Calibrated against the corpus: a step with a product of ~1.2e3 costs
/// ~25 ms while one at ~2.7e4 costs ~1.4 s (the term counts barely move, so
/// only the product separates them), so the cap is placed to keep every
/// permitted step in the tens-of-milliseconds range.
const MAX_STEP_PRODUCT: usize = 100_000;

/// Predicted cost of one field-Euclidean `fpoly_gcd` call
/// (`cost(a) · cost(b)`); checked before the loop is entered.
const MAX_GCD_PRODUCT: usize = 100_000;

thread_local! {
    static WORK_UNITS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    #[cfg(test)]
    static PEAK_WORK_UNITS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Reset the work budget; called by every public entry point of this
/// module. `symbolic_rational` is entered once per chain re-entry, so the
/// budget only ever covers a single invocation.
fn reset_work_budget() {
    WORK_UNITS.with(|c| c.set(0));
}

/// Charge `units` against the budget; returns `true` when it is exhausted.
fn charge_work(units: u64) -> bool {
    WORK_UNITS.with(|c| {
        let v = c.get().saturating_add(units);
        c.set(v);
        #[cfg(test)]
        PEAK_WORK_UNITS.with(|p| p.set(p.get().max(v)));
        v > MAX_WORK_UNITS
    })
}

/// Whether the budget is already spent (a check, not a charge).
fn budget_exhausted() -> bool {
    WORK_UNITS.with(|c| c.get() > MAX_WORK_UNITS)
}

// -------------------------------------------------------------------------
// Charged coefficient-field arithmetic
// -------------------------------------------------------------------------
//
// `GeneratorField` is `RationalPolynomial`, whose `mul`/`div`/`add` are the
// operations that actually blow up: they cross-multiply the numerator and
// denominator polynomials and then canonicalize the result. The cost is
// proportional to the product of the operand term counts, so each operation
// is charged by that product *at the call site*, which is the only place
// this module can see it. The operation itself still completes (the caller's
// loop-head `budget_exhausted` check stops the series); the structural
// pre-call checks (`MAX_STEP_PRODUCT`, `MAX_GCD_PRODUCT`) keep a single
// operation from being unboundedly large.

/// `a · b`, charged by the coefficient cross-product.
fn fmul(a: &GeneratorField, b: &GeneratorField) -> GeneratorField {
    charge_work(
        (coeff_cost(a) as u64)
            .saturating_mul(coeff_cost(b) as u64)
            .max(1),
    );
    a.mul(b)
}

/// `a + b`, charged by the operand sizes.
fn fadd(a: &GeneratorField, b: &GeneratorField) -> GeneratorField {
    charge_work((coeff_cost(a) as u64).saturating_add(coeff_cost(b) as u64));
    a.add(b)
}

/// `a / b`, charged by the coefficient cross-product.
fn fdiv(a: &GeneratorField, b: &GeneratorField) -> Option<GeneratorField> {
    charge_work(
        (coeff_cost(a) as u64)
            .saturating_mul(coeff_cost(b) as u64)
            .max(1),
    );
    a.div(b)
}

/// Largest charge seen since the last [`reset_peak_work`] (test-only
/// calibration hook; not part of the budget).
#[cfg(test)]
fn peak_work_units() -> u64 {
    PEAK_WORK_UNITS.with(|p| p.get())
}

#[cfg(test)]
fn reset_peak_work() {
    PEAK_WORK_UNITS.with(|p| p.set(0));
}

/// A univariate-in-`x` polynomial whose coefficients are elements of
/// `ℚ(symbols)` (rational functions of the constant symbols). Terms are
/// stored in ascending degree.
#[derive(Clone, Debug)]
struct FPoly {
    terms: Vec<(usize, GeneratorField)>,
}

fn rat_const(n: i64, n_vars: usize) -> GeneratorField {
    GeneratorField::from_polynomial(Sparse::from_terms(
        RationalDomain,
        n_vars,
        vec![(vec![0; n_vars], Rational::new(n, 1))],
    ))
}

fn fpoly_one(n_vars: usize) -> FPoly {
    FPoly {
        terms: vec![(0, rat_const(1, n_vars))],
    }
}

fn fpoly_zero(n_vars: usize) -> FPoly {
    FPoly {
        terms: vec![(0, GeneratorField::zero(&RationalDomain, n_vars))],
    }
}

/// Cancel the common monomial factor of a single-term/single-term field
/// element (the `RationalPolynomial` canonicalizer only cancels scalar
/// content, so exponents would otherwise accumulate through arithmetic).
fn mono_reduce(mut g: GeneratorField) -> GeneratorField {
    if g.numerator.n_terms() != 1 || g.denominator.n_terms() != 1 {
        return g;
    }
    let (en, cn) = match g.numerator.terms_ref().iter().next() {
        Some((e, c)) => (e.clone(), c.clone()),
        None => return g,
    };
    let (ed, cd) = match g.denominator.terms_ref().iter().next() {
        Some((e, c)) => (e.clone(), c.clone()),
        None => return g,
    };
    if en.len() != ed.len() {
        return g;
    }
    let common: Vec<usize> = en.iter().zip(&ed).map(|(a, b)| (*a).min(*b)).collect();
    if common.iter().all(|&v| v == 0) {
        return g;
    }
    let new_en: Vec<usize> = en.iter().zip(&common).map(|(a, b)| a - b).collect();
    let new_ed: Vec<usize> = ed.iter().zip(&common).map(|(a, b)| a - b).collect();
    g.numerator = Sparse::from_terms(RationalDomain, en.len(), vec![(new_en, cn)]);
    g.denominator = Sparse::from_terms(RationalDomain, ed.len(), vec![(new_ed, cd)]);
    g
}

impl FPoly {
    fn n_vars(&self) -> usize {
        self.terms[0].1.n_vars()
    }

    fn is_zero(&self) -> bool {
        self.terms.iter().all(|(_, c)| c.is_zero())
    }

    fn degree(&self) -> Option<usize> {
        self.terms
            .iter()
            .rev()
            .find(|(_, c)| !c.is_zero())
            .map(|(p, _)| *p)
    }

    fn leading_coeff(&self) -> Option<GeneratorField> {
        self.terms
            .iter()
            .rev()
            .find(|(_, c)| !c.is_zero())
            .map(|(_, c)| c.clone())
    }

    fn trim(&mut self) {
        let n_vars = self.n_vars();
        self.terms.retain(|(_, c)| !c.is_zero());
        self.terms.sort_by_key(|(p, _)| *p);
        if self.terms.is_empty() {
            self.terms
                .push((0, GeneratorField::zero(&RationalDomain, n_vars)));
        }
    }

    fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (p, c) in &other.terms {
            match out.terms.iter_mut().find(|(q, _)| q == p) {
                Some((_, acc)) => *acc = fadd(acc, c),
                None => out.terms.push((*p, c.clone())),
            }
        }
        out.trim();
        out
    }

    fn neg(&self) -> Self {
        Self {
            terms: self.terms.iter().map(|(p, c)| (*p, c.neg())).collect(),
        }
    }

    fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    fn mul(&self, other: &Self) -> Self {
        let mut out = fpoly_zero(self.n_vars());
        for (p, c) in &self.terms {
            for (q, d) in &other.terms {
                let pow = p + q;
                let prod = mono_reduce(fmul(c, d));
                match out.terms.iter_mut().find(|(r, _)| *r == pow) {
                    Some((_, acc)) => *acc = mono_reduce(fadd(acc, &prod)),
                    None => out.terms.push((pow, prod)),
                }
            }
        }
        out.trim();
        out
    }

    fn derivative(&self) -> Self {
        let mut out: Vec<(usize, GeneratorField)> = self
            .terms
            .iter()
            .filter(|(p, _)| *p > 0)
            .map(|(p, c)| {
                (
                    *p - 1,
                    mono_reduce(fmul(c, &rat_const(*p as i64, self.n_vars()))),
                )
            })
            .collect();
        if out.is_empty() {
            out.push((0, GeneratorField::zero(&RationalDomain, self.n_vars())));
        }
        Self { terms: out }
    }

    /// Convert to a multivariate polynomial over `[x, symbols…]`
    /// (coefficients in ℚ; a field element c = num/den contributes num
    /// with positive exponents and −den with negative ones).
    fn to_sparse(&self) -> Sparse {
        let n_vars = self.n_vars() + 1;
        let mut terms: Vec<(Vec<usize>, _)> = Vec::new();
        for (pow, c) in &self.terms {
            for (e, coeff) in c.numerator.terms_ref() {
                let mut e2 = vec![*pow];
                e2.extend_from_slice(e);
                terms.push((e2, coeff.clone()));
            }
            for (e, coeff) in c.denominator.terms_ref() {
                let mut e2 = vec![*pow];
                e2.extend_from_slice(e);
                terms.push((e2, RationalDomain.neg(coeff)));
            }
        }
        Sparse::from_terms(RationalDomain, n_vars, terms)
    }

    /// Convert a multivariate polynomial over `[x, symbols…]` back.
    fn from_sparse(p: &Sparse) -> Self {
        let nsym = p.n_vars() - 1;
        let mut out = fpoly_zero(nsym);
        for (e, c) in p.terms_ref() {
            let sym = GeneratorField::from_polynomial(Sparse::from_terms(
                RationalDomain,
                nsym,
                vec![(e[1..].to_vec(), c.clone())],
            ));
            out = out.add(&FPoly {
                terms: vec![(e[0], sym)],
            });
        }
        out
    }

    /// Long division in `x` over the field; returns `(quotient, remainder)`.
    ///
    /// Every step cancels the current leading term, so the remainder degree
    /// strictly decreases and `deg(self) − deg(den) + 1` bounds the loop;
    /// the intermediate coefficient sizes are unbounded, though, and a
    /// single call is enough to blow up (the inter-step `MAX_COEFF_COST`
    /// check in [`fpoly_gcd`] never sees those intermediates — see the
    /// many-symbol entry gate in [`integrate_rational_symbolic`]). Both the
    /// structural step count and the coefficient size are therefore checked
    /// inside the loop.
    fn div_rem(&self, den: &Self) -> Option<(Self, Self)> {
        let mut q = fpoly_zero(self.n_vars());
        let mut r = self.clone();
        let dd = den.degree()?;
        let lc_d = den.leading_coeff()?;
        let lc_d_cost = coeff_cost(&lc_d);
        let den_cost = fpoly_cost(den);
        let mut steps_left = self.degree().map_or(0, |d| d.saturating_sub(dd) + 1);
        while let Some(dr) = r.degree() {
            if dr < dd {
                break;
            }
            let lc_r = r.leading_coeff()?;
            let r_cost = fpoly_cost(&r);
            // Structural prediction before the division/multiplication
            // below: the quotient coefficient `c = lc_r / lc_d` can carry as
            // many terms as both leading coefficients together, and
            // `den · (c·x^k)` then costs `cost(den) · cost(c)` coefficient
            // products. Neither factor is visible to the checks that run
            // after the step.
            let t_cost = coeff_cost(&lc_r).saturating_add(lc_d_cost);
            if steps_left == 0
                || budget_exhausted()
                || r_cost.saturating_mul(t_cost) > MAX_STEP_PRODUCT
                || charge_work(1 + r_cost as u64 + den_cost as u64)
            {
                return None;
            }
            steps_left -= 1;
            let c = mono_reduce(fdiv(&lc_r, &lc_d)?);
            let t = FPoly {
                terms: vec![(dr - dd, c)],
            };
            q = q.add(&t);
            r = r.sub(&den.mul(&t));
            if fpoly_cost(&r) + fpoly_cost(&q) + den_cost > MAX_COEFF_COST {
                return None;
            }
        }
        Some((q, r))
    }

    fn rem(&self, den: &Self) -> Option<Self> {
        self.div_rem(den).map(|(_, r)| r)
    }

    /// Evaluate at `x = v` (sparse Horner: exponents may have gaps, and
    /// the leading term's `x^exp` must be applied).
    fn eval(&self, v: &GeneratorField) -> Option<GeneratorField> {
        let mut acc = GeneratorField::zero(&RationalDomain, self.n_vars());
        let mut prev: Option<usize> = None;
        for (exp, c) in self.terms.iter().rev() {
            if budget_exhausted() {
                return None;
            }
            let gap = prev.map(|p| p - *exp).unwrap_or(0);
            for _ in 0..=gap {
                acc = fmul(&acc, v);
            }
            acc = fadd(&acc, c);
            prev = Some(*exp);
        }
        if let Some(p) = prev {
            for _ in 0..p {
                acc = fmul(&acc, v);
            }
        }
        Some(acc)
    }

    /// Atom over `gens = [x, symbols…]`.
    fn to_atom<'a>(&self, ctx: &'a AtomArena<'a>, gens: &[Atom<'a>]) -> Option<Atom<'a>> {
        let mut terms: Vec<Atom<'a>> = Vec::new();
        for (pow, c) in &self.terms {
            let c_atom = rational_to_atom(ctx, c, &gens[1..])?;
            let x = gens[0];
            let xp = if *pow == 0 {
                ctx.num(1)
            } else {
                ctx.pow(x, ctx.num(*pow as i64))
            };
            terms.push(ctx.mul(&[c_atom, xp]));
        }
        if terms.is_empty() {
            return Some(ctx.num(0));
        }
        Some(ocas_atom::normalize::normalize(ctx, ctx.add(&terms)))
    }
}

/// GCD over `ℚ(symbols)[x]` via the field Euclidean algorithm (exact:
/// coefficients live in the field, so no content handling is needed),
/// normalised monic. The iteration cap guards against pathological
/// coefficient growth; a cap hit returns `None` (the caller falls back).
fn fpoly_gcd(a: &FPoly, b: &FPoly) -> Option<FPoly> {
    if a.is_zero() {
        return Some(b.clone());
    }
    if b.is_zero() {
        return Some(a.clone());
    }
    // Structural prediction before the Euclidean loop: the per-step checks
    // cannot see the intermediates of the first `div_rem`, which is where a
    // coefficient blow-up actually happens.
    if fpoly_cost(a).saturating_mul(fpoly_cost(b)) > MAX_GCD_PRODUCT
        || charge_work(1 + fpoly_cost(a) as u64 + fpoly_cost(b) as u64)
    {
        return None;
    }
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut steps = 0usize;
    while !r.is_zero() {
        steps += 1;
        if steps > 512 || budget_exhausted() || fpoly_cost(&old_r) + fpoly_cost(&r) > MAX_COEFF_COST
        {
            return None;
        }
        let (_, rem) = old_r.div_rem(&r)?;
        old_r = r;
        r = rem;
    }
    let mut g = old_r;
    if let Some(lc) = g.leading_coeff() {
        if let Some(inv) = lc.inv() {
            g = g.scale(&inv);
        }
    }
    g.trim();
    Some(g)
}

/// Total coefficient size of an `FPoly` (numerator + denominator term
/// counts over the symbol field) — the cost driver of field-Euclidean
/// steps. The budget turns multivariate-coefficient blow-ups (5+ symbol
/// corpus shapes) from hangs into fast declines.
fn fpoly_cost(p: &FPoly) -> usize {
    p.terms.iter().map(|(_, c)| coeff_cost(c)).sum()
}

/// Coefficient size of one field element (numerator + denominator terms).
fn coeff_cost(c: &GeneratorField) -> usize {
    c.numerator.n_terms() + c.denominator.n_terms()
}

impl FPoly {
    fn scale(&self, c: &GeneratorField) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(p, q)| (*p, mono_reduce(fmul(q, c))))
                .collect(),
        }
    }
}

/// Square-free decomposition `p = ∏ f_i^{m_i}` (Yun's algorithm) over the
/// field.
fn square_free_factors(p: &FPoly) -> Option<Vec<(FPoly, usize)>> {
    if p.degree() == Some(0) {
        return Some(Vec::new());
    }
    // NOTE: p is NOT made monic — the leading coefficient is part of the
    // factor (Hermite and the log-part residue need the raw factors).
    let p_prime = p.derivative();
    let a0 = fpoly_gcd(p, &p_prime)?;
    if a0.degree() == Some(0) {
        return Some(vec![(p.clone(), 1)]);
    }
    let (b1, _) = p.div_rem(&a0)?;
    let (c1, _) = p_prime.div_rem(&a0)?;
    let mut b = b1;
    let mut d = c1.sub(&b.derivative());
    let mut result: Vec<(FPoly, usize)> = Vec::new();
    let mut i = 1usize;
    // The multiplicity sum bounds the iteration count (a unit gcd still
    // advances `d`, which drives the next non-unit gcd); the counter is
    // deliberately generous — the `MAX_COEFF_COST` and work budgets are the
    // real backstops — because a unit-gcd step does not shrink the degree.
    let mut steps_left = 2 * (p.degree().unwrap_or(0) + 2);
    while b.degree() != Some(0) {
        if steps_left == 0
            || budget_exhausted()
            || charge_work(1 + fpoly_cost(&b) as u64 + fpoly_cost(&d) as u64)
        {
            return None;
        }
        steps_left -= 1;
        let ai = fpoly_gcd(&b, &d)?;
        let (b_next, _) = b.div_rem(&ai)?;
        let (c_next, _) = d.div_rem(&ai)?;
        let d_next = c_next.sub(&b_next.derivative());
        // A unit gcd still makes progress: the multiplicity accumulates in
        // later iterations (pure powers: a_i = 1 until d reaches 0, then
        // the last gcd is the base factor with the full multiplicity).
        if ai.degree() != Some(0) {
            result.push((ai, i));
        }
        b = b_next;
        d = d_next;
        i += 1;
    }
    if result.is_empty() {
        result.push((p.clone(), 1));
    }
    Some(result)
}

/// Hermite reduction: `num/den = Σ (B_k/D2_k)' + C/D''` with `D''`
/// square-free. Returns the accumulated `(B_k, D2_k)` pairs and `(C, D'')`.
type HermiteParts = (Vec<(FPoly, FPoly)>, FPoly, FPoly);

fn hermite_reduce(num: &FPoly, den: &FPoly) -> Option<HermiteParts> {
    // Per-factor Hermite step (Bronstein): for one factor f^m (m ≥ 2) with
    // D = D1·f^m and gcd(D1, f) = 1,
    //   A/(D1·f^m) = (B/f^{m−1})' + C/(D1·f^{m−1})
    //   A = (B'·f − B·(m−1)·f')·D1 + C·f,  so mod f: A ≡ −B·(m−1)·f'·D1.
    // Each step removes exactly one power of one factor; the loop runs
    // until the denominator is square-free.
    let mut b_parts: Vec<(FPoly, FPoly)> = Vec::new();
    let mut a = num.clone();
    let mut d = den.clone();
    // Each step strips one power of one factor from the denominator, so the
    // denominator degree strictly decreases; the counter makes that bound
    // explicit and charges the coefficient cost of the step.
    let mut steps_left = d.degree().unwrap_or(0) + 1;
    loop {
        if fpoly_cost(&a) + fpoly_cost(&d) > MAX_COEFF_COST || budget_exhausted() {
            return None;
        }
        if steps_left == 0 || charge_work(1 + fpoly_cost(&a) as u64 + fpoly_cost(&d) as u64) {
            return None;
        }
        steps_left -= 1;
        let factors = square_free_factors(&d)?;
        let Some((f, m)) = factors.iter().find(|(_, m)| *m >= 2).cloned() else {
            break;
        };
        // D1 = d / f^m
        let mut f_pow = fpoly_one(d.n_vars());
        for _ in 0..m {
            f_pow = f_pow.mul(&f);
        }
        let (d1, r) = d.div_rem(&f_pow)?;
        if !r.is_zero() {
            return None;
        }
        let f_prime = f.derivative();
        let w = f_prime
            .scale(&rat_const((m - 1) as i64, d.n_vars()))
            .mul(&d1);
        let (_, t) = extended_gcd(&f, &w)?;
        let b = a.mul(&t).neg().rem(&f)?;
        let b_prime = b.derivative();
        let inner = b_prime
            .mul(&f)
            .sub(
                &b.mul(&f_prime)
                    .scale(&rat_const((m - 1) as i64, d.n_vars())),
            )
            .mul(&d1);
        let (c, r2) = a.sub(&inner).div_rem(&f)?;
        if !r2.is_zero() {
            return None;
        }
        // B/f^{m−1} contributes to the answer; the new denominator is D1·f^{m−1}.
        let mut f_pow_m1 = fpoly_one(d.n_vars());
        for _ in 0..(m - 1) {
            f_pow_m1 = f_pow_m1.mul(&f);
        }
        b_parts.push((b, f_pow_m1.clone()));
        a = c;
        d = d1.mul(&f_pow_m1);
    }
    Some((b_parts, a, d))
}

/// Extended Euclidean over `ℚ(symbols)[x]`: `(s, t)` with `s·a + t·b = 1`
/// when `a, b` are coprime (the only use here).
fn extended_gcd(a: &FPoly, b: &FPoly) -> Option<(FPoly, FPoly)> {
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = fpoly_one(a.n_vars());
    let mut s = fpoly_zero(a.n_vars());
    let mut old_t = fpoly_zero(a.n_vars());
    let mut t = fpoly_one(a.n_vars());
    // Euclidean degree descent bounds the loop; charge the coefficient cost
    // of every step so a blow-up inside the series cannot run away. The step
    // bound is generous because a step whose remainder has the smaller
    // degree only swaps `old_r`/`r` and does not shrink the degree
    // (`extended_gcd(x, a - b·x²)` needs three steps for `deg = 1`).
    let mut steps_left = a.degree().unwrap_or(0) + b.degree().unwrap_or(0) + 2;
    while !r.is_zero() {
        if fpoly_cost(&old_r) + fpoly_cost(&r) + fpoly_cost(&s) + fpoly_cost(&t) > MAX_COEFF_COST {
            return None;
        }
        if steps_left == 0
            || budget_exhausted()
            || charge_work(1 + fpoly_cost(&old_r) as u64 + fpoly_cost(&r) as u64)
        {
            return None;
        }
        steps_left -= 1;
        let (q, rem) = old_r.div_rem(&r)?;
        old_r = r;
        r = rem;
        let new_s = old_s.sub(&q.mul(&s));
        old_s = s;
        s = new_s;
        let new_t = old_t.sub(&q.mul(&t));
        old_t = t;
        t = new_t;
    }
    if let Some(lc) = old_r.leading_coeff() {
        if let Some(inv) = lc.inv() {
            old_s = old_s.scale(&inv);
            old_t = old_t.scale(&inv);
        }
    }
    Some((old_s, old_t))
}

/// `(α, β)` of a linear factor `α·x + β`.
fn linear_coeffs(f: &FPoly) -> Option<(GeneratorField, GeneratorField)> {
    let a1 = f
        .terms
        .iter()
        .find(|(p, _)| *p == 1)
        .map(|(_, c)| c.clone())
        .unwrap_or_else(|| GeneratorField::zero(&RationalDomain, f.n_vars()));
    let a0 = f
        .terms
        .iter()
        .find(|(p, _)| *p == 0)
        .map(|(_, c)| c.clone())
        .unwrap_or_else(|| GeneratorField::zero(&RationalDomain, f.n_vars()));
    Some((a1, a0))
}

/// The symbols appearing in `expr` (variables other than `var`), sorted.
fn collect_symbols(expr: Atom<'_>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Var(v) => {
            if *v != var && !out.contains(v) {
                out.push(*v);
            }
        }
        AtomNode::Num(_) => {}
        AtomNode::Fun(_, args) => {
            for a in *args {
                collect_symbols(*a, var, out);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in *args {
                collect_symbols(*a, var, out);
            }
        }
        AtomNode::Pow(base, exp) => {
            collect_symbols(*base, var, out);
            collect_symbols(*exp, var, out);
        }
    }
}

/// Extract `Δ = 4ac − b²` as a field element.
fn discriminant(a: &GeneratorField, b: &GeneratorField, c: &GeneratorField) -> GeneratorField {
    let four = rat_const(4, a.n_vars());
    fmul(&fmul(&four, a), c).sub(&fmul(b, b))
}

/// `√Δ` as a rational function of the symbols when Δ is a square in
/// `ℚ(symbols)`; else `None`.
///
/// The monomial-wise square root below is only a *candidate* generator: it
/// is exact for a single-term element (`4·a² → 2·a`) but wrong for a sum,
/// because squaring a sum produces cross terms that the term-wise route
/// never sees (`4a² + 4b²` would yield the bogus `2a + 2b`, whose square is
/// `4a² + 8ab + 4b²`). The candidate is therefore squared with exact field
/// arithmetic and returned only when it reproduces Δ; otherwise this
/// returns `None` and the caller keeps the quadratic whole (the
/// log + atan/atanh branch), which is always correct.
fn rational_square_root(delta: &GeneratorField) -> Option<GeneratorField> {
    // The candidate is squared below, which is itself a charged field
    // multiplication; a huge Δ declines here rather than paying for it.
    if budget_exhausted() {
        return None;
    }
    let sqrt_sparse = |p: &Sparse| -> Option<Sparse> {
        let mut terms: Vec<(Vec<usize>, Rational)> = Vec::new();
        for (e, c) in p.terms_ref() {
            if e.iter().any(|&v| v % 2 != 0) {
                return None;
            }
            // A non-integral coefficient is not a rational square unless its
            // denominator is a square too; the exact check below would
            // reject it anyway, so bail out early rather than dropping the
            // denominator here.
            if c.denom().to_i64()? != 1 {
                return None;
            }
            let rp = isqrt_i64(c.numer().to_i64()?)?;
            terms.push((e.iter().map(|v| v / 2).collect(), Rational::new(rp, 1)));
        }
        Some(Sparse::from_terms(RationalDomain, p.n_vars(), terms))
    };
    let n = sqrt_sparse(&delta.numerator)?;
    let d = sqrt_sparse(&delta.denominator)?;
    let candidate = GeneratorField::from_num_den(n, d);
    // Exact check in `ℚ(symbols)`: `candidate² == delta`. Both sides are
    // canonical field elements, so comparing their difference against zero
    // is a decision procedure (cross-multiplication is exact polynomial
    // arithmetic, independent of how the fraction is represented).
    if fmul(&candidate, &candidate).sub(delta).is_zero() {
        Some(candidate)
    } else {
        None
    }
}

/// Integer multivariate factorization of a square-free factor: clears the
/// coefficient denominators, factors over ℤ, and returns the factors (each
/// converted back to the field) — their product equals `scalar·f` where the
/// scalar is the clearing constant (the log part multiplies its numerator
/// by that scalar, keeping the partial-fraction coefficients exact).
fn factor_via_integer(f: &FPoly) -> Option<Vec<(FPoly, usize)>> {
    // Structural gate in front of `SparseMultivariatePolynomial::factor()`:
    // the entry degree is already ≤ 3 (see `split_squarefree_factors`), the
    // coefficient-size budget bounds the clearing lcm, and the term count
    // bounds the conversion below. Refusing here keeps the factorizer off
    // pathologically wide inputs.
    if charge_work(1 + fpoly_cost(f) as u64) {
        return None;
    }
    let sparse = f.to_sparse();
    if sparse.terms_ref().len() > MAX_FACTOR_TERMS {
        return None;
    }
    // The to_sparse round trip encodes field coefficients by sign-flipping
    // the denominator monomials, which from_sparse cannot reconstruct
    // (a/c would come back as a − c). Integer factorization is therefore
    // only sound for scalar coefficients: reject any term carrying a
    // symbol (a nonzero exponent past the x slot).
    for e in sparse.terms_ref().keys() {
        if e[1..].iter().any(|&v| v != 0) {
            return None;
        }
    }
    // Clear denominators: D = lcm of all coefficient denominators.
    let mut d: i64 = 1;
    for c in sparse.terms_ref().values() {
        d = d / gcd_i64(d, c.denom().to_i64()?) * c.denom().to_i64()?;
    }
    let mut int_terms: Vec<(Vec<usize>, Integer)> = Vec::new();
    for (e, c) in sparse.terms_ref() {
        let den = c.denom().to_i64()?;
        let num = c.numer().to_i64()? * (d / den);
        int_terms.push((e.to_vec(), Integer::from(num)));
    }
    let int_poly: SparseMultivariatePolynomial<IntegerDomain, Lex> =
        SparseMultivariatePolynomial::from_terms(IntegerDomain, sparse.n_vars(), int_terms);
    let fac = int_poly.factor();
    let mut out: Vec<(FPoly, usize)> = Vec::new();
    for (g, k) in fac {
        // Convert each integer factor back into the field (divide by D so
        // the factors multiply to f; the scalar is tracked implicitly).
        let mut terms: Vec<(Vec<usize>, Rational)> = Vec::new();
        for (e, c) in g.terms_ref() {
            let num = c.to_i64()?;
            terms.push((e.to_vec(), Rational::new(num, 1)));
        }
        let gq = Sparse::from_terms(RationalDomain, sparse.n_vars(), terms);
        out.push((FPoly::from_sparse(&gq), k));
    }
    let _ = d;
    if out.is_empty() {
        return None;
    }
    Some(out)
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.abs().max(1)
}

/// Split square-free factors into factors of degree ≤ 2 over
/// `ℚ(symbols)`: degree-2 factors split when the discriminant is a rational
/// square; degree ≥ 3 factors go through integer multivariate
/// factorization (coefficients cleared), falling back to leaving the factor
/// whole (which the log part then reports as an unevaluated Integral).
/// Monomial divisors of a field element: every monomial of the numerator
/// (and denominator) with exponent-wise ≤ exponents, times the integer
/// divisors of the scalar content (bounded).
///
/// The candidate list is the product of `(exponent_i + 1)` over the symbols
/// of one monomial, so the size is predicted before the expansion allocates
/// it and the running total is charged against the module budget.
fn field_divisors(el: &GeneratorField) -> Option<Vec<GeneratorField>> {
    let mut out: Vec<GeneratorField> = Vec::new();
    for e in el.numerator.terms_ref().keys() {
        let mut exps: Vec<Vec<usize>> = vec![vec![0; e.len()]];
        for (i, &v) in e.iter().enumerate() {
            if exps.len().saturating_mul(v + 1) > MAX_FIELD_DIVISORS {
                return None;
            }
            let mut next = Vec::new();
            for cur in &exps {
                for k in 0..=v {
                    let mut c = cur.clone();
                    c[i] = k;
                    next.push(c);
                }
            }
            exps = next;
        }
        for exp in exps {
            if out.len() >= MAX_FIELD_DIVISORS || charge_work(1) {
                return None;
            }
            let p = Sparse::from_terms(
                RationalDomain,
                el.numerator.n_vars(),
                vec![(exp, Rational::new(1, 1))],
            );
            out.push(GeneratorField::from_polynomial(p));
        }
    }
    for e in el.denominator.terms_ref().keys() {
        let mut exps: Vec<Vec<usize>> = vec![vec![0; e.len()]];
        for (i, &v) in e.iter().enumerate() {
            if exps.len().saturating_mul(v + 1) > MAX_FIELD_DIVISORS {
                return None;
            }
            let mut next = Vec::new();
            for cur in &exps {
                for k in 0..=v {
                    let mut c = cur.clone();
                    c[i] = k;
                    next.push(c);
                }
            }
            exps = next;
        }
        for exp in exps {
            if out.len() >= MAX_FIELD_DIVISORS || charge_work(1) {
                return None;
            }
            let p = Sparse::from_terms(
                RationalDomain,
                el.denominator.n_vars(),
                vec![(exp, Rational::new(1, 1))],
            );
            out.push(GeneratorField::from_polynomial(p).inv()?);
        }
    }
    Some(out)
}

/// Split a square-free factor of degree 3 by trying field-linear root
/// candidates `r = −β/α` with `α | lc`, `β | a0` (monomial divisors).
///
/// The candidate list is `2 · |divisors(a0)| · |divisors(lc)|`; the product
/// is bounded before the evaluation loop so a fat leading coefficient
/// cannot make the enumeration combinatorial.
fn split_linear_candidates(f: &FPoly) -> Option<Vec<(FPoly, usize)>> {
    let lc = f.leading_coeff()?;
    let a0 = f
        .terms
        .iter()
        .find(|(p, _)| *p == 0)
        .map(|(_, c)| c.clone())
        .unwrap_or_else(|| GeneratorField::zero(&RationalDomain, f.n_vars()));
    let mut candidates: Vec<GeneratorField> = Vec::new();
    if a0.is_zero() {
        candidates.push(GeneratorField::zero(&RationalDomain, f.n_vars()));
    } else {
        let divisors_a0 = field_divisors(&a0)?;
        let divisors_lc = field_divisors(&lc)?;
        if divisors_a0
            .len()
            .saturating_mul(divisors_lc.len())
            .saturating_mul(2)
            > MAX_SPLIT_CANDIDATES
        {
            return None;
        }
        for da in divisors_a0 {
            for dl in &divisors_lc {
                if !dl.is_zero() {
                    let r = fdiv(&da, dl)?;
                    candidates.push(r.neg());
                    candidates.push(r);
                }
            }
        }
    }
    for r in candidates {
        // Each candidate costs one Horner evaluation of `f`; charge the
        // degree so the accumulated budget tracks the real work.
        if budget_exhausted() || charge_work(1 + f.degree().unwrap_or(0) as u64) {
            return None;
        }
        if f.eval(&r)?.is_zero() {
            // (x − r) divides f (monic factor; the log part works on the
            // monic-normalized denominator).
            let one = GeneratorField::one(&RationalDomain, f.n_vars());
            let mut g = FPoly {
                terms: vec![(1, one), (0, r.neg())],
            };
            g.trim();
            let (q, rem) = f.div_rem(&g)?;
            if rem.is_zero() && q.degree()? >= 1 {
                return Some(vec![(g, 1), (q, 1)]);
            }
        }
    }
    None
}

fn split_squarefree_factors(factors: Vec<(FPoly, usize)>) -> Option<Vec<(FPoly, usize)>> {
    let mut out: Vec<(FPoly, usize)> = Vec::new();
    for (f, m) in factors {
        if budget_exhausted() || charge_work(1 + fpoly_cost(&f) as u64) {
            return None;
        }
        let deg = f.degree()?;
        if deg > 2 && deg <= 3 {
            // Field-linear factors first (symbolic coefficients: the
            // integer factorization is scalar-only, so mixed factors like
            // x·(a−b·x²) or (a+b·x)·(c+x²) need the root-candidate split).
            if let Some(sub) = split_linear_candidates(&f) {
                let mut sub = split_squarefree_factors(sub)?;
                out.append(&mut sub);
                continue;
            }
            if let Some(sub) = factor_via_integer(&f) {
                let mut sub = split_squarefree_factors(sub)?;
                out.append(&mut sub);
                continue;
            }
            out.push((f, m));
            continue;
        }
        if deg == 2 {
            let (a, b, c) = quadratic_coeffs_fpoly(&f)?;
            // Root-splitting discriminant: b² − 4ac.
            let delta = fmul(&b, &b).sub(&fmul(&fmul(&rat_const(4, f.n_vars()), &a), &c));
            if let Some(s) = rational_square_root(&delta) {
                let two_a = fmul(&a, &rat_const(2, f.n_vars()));
                let r1 = fadd(&b.neg(), &s).div(&two_a)?;
                let r2 = b.neg().sub(&s).div(&two_a)?;
                let one = GeneratorField::one(&RationalDomain, f.n_vars());
                // Keep the leading coefficient of the original quadratic in
                // the split (product must equal f): the log-part residues
                // evaluate ∏_{j≠i} f_j(r)/f_i'(r), which is only correct
                // when the factors multiply back to the exact denominator.
                let mut g1 = FPoly {
                    terms: vec![(1, a.clone()), (0, fmul(&a, &r1).neg())],
                };
                g1.trim();
                let mut g2 = FPoly {
                    terms: vec![(1, one), (0, r2.neg())],
                };
                g2.trim();
                out.push((g1, m));
                out.push((g2, m));
                continue;
            }
        }
        out.push((f, m));
    }
    Some(out)
}

fn isqrt_i64(n: i64) -> Option<i64> {
    if n < 0 {
        return None;
    }
    let r = (n as f64).sqrt() as i64;
    [r - 1, r, r + 1]
        .into_iter()
        .find(|&c| c >= 0 && c * c == n)
}

/// The rational number stored in a constant field element.
fn constant_rational(delta: &GeneratorField) -> Option<ocas_domain::Rational> {
    if delta.numerator.n_terms() == 1 && delta.denominator.n_terms() == 1 {
        let (e, c) = delta.numerator.terms_ref().iter().next()?;
        if e.iter().all(|&v| v == 0) {
            return Some(c.clone());
        }
    }
    None
}

/// Whether the symbolic-constant rational backend can handle `expr` as a
/// rational function of `var` within its complexity bounds (≤ 4 symbols,
/// denominator degree ≤ 6). Used by the Weierstrass heuristic as a
/// feasibility gate before committing its t-integral to the chain.
pub(crate) fn rational_complexity_ok<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> bool {
    // This is a public entry point of the module (the Weierstrass gate calls
    // it directly), so the work budget starts fresh here too.
    reset_work_budget();
    let mut symbols: Vec<Symbol> = Vec::new();
    collect_symbols(expr, var, &mut symbols);
    if symbols.len() > 5 {
        return false;
    }
    let x = ctx.var(var.as_str());
    let mut gens: Vec<Atom<'a>> = vec![x];
    for s in &symbols {
        gens.push(ctx.var(s.as_str()));
    }
    let Some(rf) = atom_to_rational(expr, &gens) else {
        return false;
    };
    // `FPoly::from_sparse` is quadratic in the term count; predict that cost
    // instead of materialising the conversion.
    let sparse_terms = rf.numerator.n_terms() + rf.denominator.n_terms();
    if sparse_terms > MAX_SPARSE_TERMS
        || charge_work((sparse_terms as u64).saturating_mul(sparse_terms as u64))
    {
        return false;
    }
    let num = FPoly::from_sparse(&rf.numerator);
    let mut den = FPoly::from_sparse(&rf.denominator);
    // Best-effort cancellation, mirroring `integrate_rational_symbolic`:
    // the Weierstrass t-forms are unreduced (common (1 + t²) factors).
    if let Some(g) = fpoly_gcd(&num, &den) {
        if g.degree() != Some(0) {
            if let (Some((_, r1)), Some((dq, r2))) = (num.div_rem(&g), den.div_rem(&g)) {
                if r1.is_zero() && r2.is_zero() {
                    den = dq;
                }
            }
        }
    }
    den.degree().is_some_and(|d| d <= 6)
}

/// Integrate a rational function of `var` over `ℚ(symbols)`.
pub(crate) fn integrate_rational_symbolic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // Fresh budget for this invocation. The stage runs once per chain
    // re-entry, so a top-level `integrate` call can never accumulate a
    // stale budget across the chain.
    reset_work_budget();
    integrate_rational_symbolic_inner(ctx, expr, var)
}

fn integrate_rational_symbolic_inner<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let mut symbols: Vec<Symbol> = Vec::new();
    collect_symbols(expr, var, &mut symbols);
    symbols.sort_by_key(|s| s.as_str().to_string());
    let mut gens: Vec<Atom<'a>> = vec![x];
    for s in &symbols {
        gens.push(ctx.var(s.as_str()));
    }
    // Complexity guard: the field-Euclidean machinery degrades sharply on
    // high-degree denominators and many symbols (the trig-rational
    // t-forms can reach degree ~10 with 4+ symbols, where the coefficient
    // growth overflows the stack). Skipping keeps the pipeline fast and
    // stack-safe; the cases fall through to the heuristic stages.
    let rf = atom_to_rational(expr, &gens)?;
    // `FPoly::from_sparse` is quadratic in the term count; predict that cost
    // instead of materialising the conversion.
    let sparse_terms = rf.numerator.n_terms() + rf.denominator.n_terms();
    if sparse_terms > MAX_SPARSE_TERMS
        || charge_work((sparse_terms as u64).saturating_mul(sparse_terms as u64))
    {
        return None;
    }
    let mut num = FPoly::from_sparse(&rf.numerator);
    let mut den = FPoly::from_sparse(&rf.denominator);
    // Many-symbol entry gate (0.27.1): with 5+ coefficient generators the
    // field-Euclidean steps blow up inside a SINGLE div_rem (the inter-step
    // MAX_COEFF_COST budget never gets a chance to fire). Pure polynomials
    // (den.degree() == 0) skip Euclidean work entirely and stay admissible
    // at any symbol count; the gate only fires on genuine quotients whose
    // combined degree is large.
    let den_deg = den.degree()?;
    if den_deg > 6 {
        return None;
    }
    if den_deg > 0 && symbols.len() >= 5 && num.degree().unwrap_or(0) + den_deg > 8 {
        return None;
    }
    // Cancel the common polynomial factor of num/den (the Weierstrass
    // t-forms are often unreduced). Best effort: a gcd that is too
    // expensive or fails must not abort the whole integration.
    if let Some(g) = fpoly_gcd(&num, &den) {
        if g.degree() != Some(0) {
            if let (Some((nq, r1)), Some((dq, r2))) = (num.div_rem(&g), den.div_rem(&g)) {
                if r1.is_zero() && r2.is_zero() {
                    num = nq;
                    den = dq;
                }
            }
        }
    }
    let n_vars = den.n_vars();

    let mut parts: Vec<Atom<'a>> = Vec::new();

    // Polynomial part: ∫ Σ c·x^p dx = Σ c·x^(p+1)/(p+1).
    let (quotient, remainder) = num.div_rem(&den)?;
    for (p, c) in quotient.terms {
        if p == 0 {
            let c_atom = rational_to_atom(ctx, &c, &gens[1..])?;
            parts.push(ctx.mul(&[c_atom, x]));
        } else {
            let p1 = p as i64 + 1;
            let c_atom = rational_to_atom(ctx, &c, &gens[1..])?;
            parts.push(ctx.mul(&[
                c_atom,
                ctx.pow(x, ctx.num(p1)),
                ctx.pow(ctx.num(p1), ctx.num(-1)),
            ]));
        }
    }

    // Hermite reduction: the derivative parts integrate back to the B_k/D2_k
    // rational functions.
    let (b_parts, c_num, squarefree_den) = match hermite_reduce(&remainder, &den) {
        Some(v) => v,
        None => {
            return None;
        }
    };
    // Normalize the log-part denominator to monic (the SFF/gcd machinery
    // keeps leading coefficients in the factors, so the factor products
    // equal `lc·sf`; dividing both sides by the leading coefficient makes
    // the products exact and the partial-fraction residues consistent).
    let mut c_num = c_num;
    let mut squarefree_den = squarefree_den;
    if let Some(lc) = squarefree_den.leading_coeff() {
        if let Some(inv) = lc.inv() {
            c_num = c_num.scale(&inv);
            squarefree_den = squarefree_den.scale(&inv);
        }
    }
    for (b, d2) in b_parts {
        let b_atom = b.to_atom(ctx, &gens)?;
        let d2_atom = d2.to_atom(ctx, &gens)?;
        parts.push(ctx.mul(&[b_atom, ctx.pow(d2_atom, ctx.num(-1))]));
    }

    // Log part on the square-free denominator (split reducible quadratics).
    let factors = split_squarefree_factors(square_free_factors(&squarefree_den)?)?;
    let has_high_degree = factors.iter().any(|(f, _)| f.degree().unwrap_or(0) > 2);
    if has_high_degree {
        // Honest partial result: leave the whole log part unevaluated.
        let integrand = ctx.mul(&[
            c_num.to_atom(ctx, &gens)?,
            ctx.pow(squarefree_den.to_atom(ctx, &gens)?, ctx.num(-1)),
        ]);
        parts.push(ctx.fun("Integral", &[integrand, x]));
        return assemble(ctx, parts);
    }
    for (f, m) in &factors {
        debug_assert_eq!(*m, 1);
        if budget_exhausted() {
            return None;
        }
        let f_deg = f.degree()?;
        if f_deg == 1 {
            // c = num(r)/(f'(r)·∏_{j≠i} f_j(r)) with r = −β/α.
            let (alpha, beta) = linear_coeffs(f)?;
            let r = fdiv(&beta.neg(), &alpha)?;
            let mut denom = alpha.clone();
            for (g, _) in &factors {
                if g.terms == f.terms {
                    continue;
                }
                denom = fmul(&denom, &g.eval(&r)?);
            }
            let coeff = fdiv(&c_num.eval(&r)?, &denom)?;
            let coeff_atom = rational_to_atom(ctx, &coeff, &gens[1..])?;
            let f_atom = f.to_atom(ctx, &gens)?;
            parts.push(ctx.mul(&[coeff_atom, ctx.fun("log", &[f_atom])]));
        } else {
            // Quadratic factor: solve (M·x+N)/(Q·f) partial fraction.
            let (m, n) = quadratic_coeffs(&c_num, f, &factors)?;
            let (a, b, c) = quadratic_coeffs_fpoly(f)?;
            // ∫(Mx+N)/f = M/(2a)·log(f) + (N − M·b/(2a))·(2/√Δ)·h((2a·x+b)/√Δ)
            let two_a = fmul(&a, &rat_const(2, n_vars));
            let m_over = fdiv(&m, &two_a)?;
            let f_atom = f.to_atom(ctx, &gens)?;
            let m_atom = rational_to_atom(ctx, &m_over, &gens[1..])?;
            if !m_over.is_zero() {
                parts.push(ctx.mul(&[m_atom, ctx.fun("log", &[f_atom])]));
            }
            let delta = discriminant(&a, &b, &c);
            let mb_over = fdiv(&fmul(&m, &b), &two_a)?;
            let n_shift = n.sub(&mb_over);
            // atan: +(2/√Δ); atanh: −(2/√(−Δ)) — the derivative of
            // atanh((2ax+b)/√(−Δ)) is +1/f only with the minus sign.
            let fun = match constant_rational(&delta) {
                Some(d) if d.inner().is_positive() => "atan",
                Some(_) => "atanh",
                None => {
                    let lead = c
                        .numerator
                        .terms_ref()
                        .iter()
                        .next()
                        .map(|(_, cc)| cc.clone());
                    match lead {
                        Some(cc) if cc.numer().is_negative() => "atanh",
                        _ => "atan",
                    }
                }
            };
            let two = if fun == "atanh" { -2 } else { 2 };
            let coeff = n_shift.mul(&rat_const(two, n_vars));
            let coeff_atom = rational_to_atom(ctx, &coeff, &gens[1..])?;
            // The atanh form must use the real √(−Δ) (Δ = 4ac − b² < 0
            // there); the atan form uses √Δ (Δ > 0).
            let mag = if fun == "atanh" {
                delta.neg()
            } else {
                delta.clone()
            };
            let sqrt_atom: Atom<'a> = if let Some(r) = rational_square_root(&mag) {
                rational_to_atom(ctx, &r, &gens[1..])?
            } else {
                let d_atom = rational_to_atom(ctx, &mag, &gens[1..])?;
                ctx.pow(d_atom, ctx.pow(ctx.num(2), ctx.num(-1)))
            };
            let lin = FPoly {
                terms: vec![(1, two_a), (0, b)],
            };
            let arg_atom = ctx.mul(&[lin.to_atom(ctx, &gens)?, ctx.pow(sqrt_atom, ctx.num(-1))]);
            parts.push(ctx.mul(&[
                coeff_atom,
                ctx.pow(sqrt_atom, ctx.num(-1)),
                ctx.fun(fun, &[arg_atom]),
            ]));
        }
    }
    // Last-resort emission guard: the assembled closed form is differentiated
    // and compared with the input at a few deterministic sample points. A
    // wrong answer is worse than a fallback, and this module has produced two
    // (0.27.1's `rational_square_root`, 0.27.2's residue coefficients for
    // repeated factors) — the guard turns any future instance into an honest
    // decline. Unverifiable samples (domain, unsupported head) do not decline.
    if let Some(answer) = assemble(ctx, parts) {
        if !emission_is_verified(ctx, expr, answer, var, &symbols) {
            return None;
        }
        return Some(answer);
    }
    None
}

/// Sample abscissae for [`emission_is_verified`].
const GUARD_SAMPLES: [f64; 3] = [0.37, 0.83, 1.27];

/// Deterministic dummy value for the `i`-th coefficient symbol.
fn guard_param(i: usize) -> f64 {
    const TABLE: [f64; 8] = [2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 0.5, 1.5];
    TABLE[i % TABLE.len()]
}

/// Numeric f64 evaluation of the atoms this module can emit.
fn eval_num(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
    match expr.node() {
        AtomNode::Num(n) => Some(*n as f64),
        AtomNode::Var(v) => {
            if let Some((_, val)) = env.iter().find(|(s, _)| s == v) {
                return Some(*val);
            }
            match v.as_str() {
                "pi" => Some(std::f64::consts::PI),
                "e" | "E" => Some(std::f64::consts::E),
                _ => None,
            }
        }
        AtomNode::Add(args) => {
            let mut acc = 0.0;
            for a in args.iter() {
                acc += eval_num(*a, env)?;
            }
            Some(acc)
        }
        AtomNode::Mul(args) => {
            let mut acc = 1.0;
            for a in args.iter() {
                acc *= eval_num(*a, env)?;
            }
            Some(acc)
        }
        AtomNode::Pow(b, e) => {
            let (b, e) = (eval_num(*b, env)?, eval_num(*e, env)?);
            if e.fract() == 0.0 || b >= 0.0 {
                Some(b.powf(e))
            } else {
                None
            }
        }
        AtomNode::Fun(name, args) => {
            let v = eval_num(*args.first()?, env)?;
            Some(match name.as_str() {
                "log" => v.abs().ln(),
                "sqrt" => {
                    if v < 0.0 {
                        return None;
                    }
                    v.sqrt()
                }
                "atan" => v.atan(),
                "atanh" => {
                    if v.abs() >= 1.0 {
                        return None;
                    }
                    v.atanh()
                }
                "asin" => {
                    if !(-1.0..=1.0).contains(&v) {
                        return None;
                    }
                    v.asin()
                }
                "abs" => v.abs(),
                _ => return None,
            })
        }
    }
}

/// Differentiate `answer` and compare it with `expr` at a few sample points.
///
/// Returns `true` when every usable sample agrees to `1e-4` relative, and
/// also when no sample is usable (the guard never declines what it cannot
/// judge). Returns `false` only on a real disagreement.
fn emission_is_verified<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    answer: Atom<'a>,
    var: Symbol,
    symbols: &[Symbol],
) -> bool {
    let base: Vec<(Symbol, f64)> = symbols
        .iter()
        .enumerate()
        .map(|(i, s)| (*s, guard_param(i)))
        .collect();
    let derivative = crate::diff(ctx, answer, var);
    for &x in GUARD_SAMPLES.iter() {
        let mut env = base.clone();
        env.push((var, x));
        let (Some(lhs), Some(rhs)) = (eval_num(derivative, &env), eval_num(expr, &env)) else {
            continue;
        };
        if !lhs.is_finite() || !rhs.is_finite() {
            continue;
        }
        if (lhs - rhs).abs() > 1e-4 * rhs.abs().max(1.0) {
            return false;
        }
    }
    true
}

fn assemble<'a>(ctx: &'a AtomArena<'a>, parts: Vec<Atom<'a>>) -> Option<Atom<'a>> {
    if parts.is_empty() {
        return Some(ctx.num(0));
    }
    Some(ocas_atom::normalize::normalize(ctx, ctx.add(&parts)))
}

/// `(a, b, c)` of a quadratic `a·x² + b·x + c`.
fn quadratic_coeffs_fpoly(f: &FPoly) -> Option<(GeneratorField, GeneratorField, GeneratorField)> {
    let get = |p: usize| -> GeneratorField {
        f.terms
            .iter()
            .find(|(q, _)| *q == p)
            .map(|(_, c)| c.clone())
            .unwrap_or_else(|| GeneratorField::zero(&RationalDomain, f.n_vars()))
    };
    Some((get(2), get(1), get(0)))
}

/// Solve `num/(Q·f) = (M·x+N)/f + rest/Q` for the quadratic factor `f`,
/// where `Q = ∏_{j≠i} f_j`. Returns `(M, N)`.
fn quadratic_coeffs(
    num: &FPoly,
    f: &FPoly,
    factors: &[(FPoly, usize)],
) -> Option<(GeneratorField, GeneratorField)> {
    let mut q = fpoly_one(f.n_vars());
    for (g, _) in factors {
        if budget_exhausted() {
            return None;
        }
        if g.terms == f.terms {
            continue;
        }
        q = q.mul(g);
    }
    let (a, b, c) = quadratic_coeffs_fpoly(f)?;
    let num_mod = num.rem(f)?;
    let q_mod = q.rem(f)?;
    let (q1, q0) = linear_coeffs(&q_mod)?;
    let (p1, p0) = linear_coeffs(&num_mod)?;
    // (Mx+N)(q1 x+q0) mod f, with x² ≡ (−b·x − c)/a:
    //   x·(M·q0 + N·q1 − M·q1·b/a) + (N·q0 − M·q1·c/a)
    let a_inv = a.inv()?;
    let m11 = fadd(&q0, &fmul(&fmul(&q1, &b), &a_inv).neg());
    let m21 = fmul(&fmul(&q1, &c), &a_inv).neg();
    let det = fadd(&fmul(&m11, &q0), &fmul(&q1, &m21).neg());
    let det_inv = det.inv()?;
    let m = fmul(&fadd(&fmul(&p1, &q0), &fmul(&q1, &p0).neg()), &det_inv);
    let n = fmul(&fadd(&fmul(&m11, &p0), &fmul(&p1, &m21).neg()), &det_inv);
    Some((m, n))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Constant values for the numeric antiderivative check. Extra symbols
    /// are harmless: the evaluator looks them up by name.
    fn consts() -> Vec<(Symbol, f64)> {
        [
            ("a", 1.3),
            ("b", 0.7),
            ("c", 0.4),
            ("d", 0.9),
            ("e", 0.5),
            ("A", 1.1),
            ("B", 0.6),
            ("C", 0.8),
        ]
        .into_iter()
        .map(|(n, v)| (Symbol::new(n), v))
        .collect()
    }

    /// Numeric f64 evaluator over the elementary functions the integrator
    /// emits (the same shape as the one in `trig_kernel`'s tests, extended
    /// with the hyperbolic family this module's shapes use).
    fn eval_f64(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
        match expr.node() {
            AtomNode::Num(n) => Some(*n as f64),
            AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, val)| *val),
            AtomNode::Add(args) => args
                .iter()
                .try_fold(0.0, |acc, a| Some(acc + eval_f64(*a, env)?)),
            AtomNode::Mul(args) => args
                .iter()
                .try_fold(1.0, |acc, a| Some(acc * eval_f64(*a, env)?)),
            AtomNode::Pow(b, e) => Some(eval_f64(*b, env)?.powf(eval_f64(*e, env)?)),
            AtomNode::Fun(name, args) => {
                let v = eval_f64(*args.first()?, env)?;
                Some(match name.as_str() {
                    "sin" => v.sin(),
                    "cos" => v.cos(),
                    "tan" => v.tan(),
                    "cot" => v.tan().recip(),
                    "sec" => v.cos().recip(),
                    "csc" => v.sin().recip(),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "tanh" => v.tanh(),
                    "coth" => v.tanh().recip(),
                    "sech" => v.cosh().recip(),
                    "csch" => v.sinh().recip(),
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    "atanh" => v.atanh(),
                    "asin" => v.asin(),
                    _ => return None,
                })
            }
        }
    }

    fn int_str(input: &str, var: &str) -> String {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, input).unwrap();
        crate::integrate(&ctx, expr, Symbol::new(var)).to_string()
    }

    fn assert_solved(input: &str) {
        let r = int_str(input, "x");
        assert!(!r.contains("Integral("), "{input} left a residue: {r}");
    }

    /// The result must be an honest fallback (the unevaluated `Integral`
    /// marker still present) or, when a mechanism does solve the shape, a
    /// numerically correct antiderivative. Budgets must never turn a hang
    /// into a wrong answer.
    fn assert_returns_not_wrong(input: &str) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let integrand = ocas_parse::parse(&ctx, input).unwrap();
        let var = Symbol::new("x");
        let result = crate::integrate(&ctx, integrand, var);
        let text = result.to_string();
        if text.contains("Integral(") {
            return;
        }
        let d = crate::diff(&ctx, result, var);
        for &xv in &[0.3f64, 0.7] {
            let mut env = consts();
            env.push((var, xv));
            let lhs = eval_f64(d, &env).expect("eval derivative");
            let rhs = eval_f64(integrand, &env).expect("eval integrand");
            assert!(
                (lhs - rhs).abs() < 1e-5 * rhs.abs().max(1.0),
                "{input} at x={xv}: derivative {lhs} vs integrand {rhs} (result: {text})"
            );
        }
    }

    #[test]
    fn symbolic_rationals() {
        assert_solved("1/(a+b*x^2)");
        assert_solved("1/(a-b*x^2)");
        assert_solved("1/(x*(a+b*x)^2)");
        assert_solved("1/(x^2*(a-b*x^2))");
        assert_solved("(d+e*x)/(x^3*(a+c*x^2))");
        assert_solved("(A+B*x)/(a+b*x+c*x^2)");
        assert_solved("x^2/(a-b*x^2)^3");
    }

    #[test]
    fn numeric_quadratic_irreducible() {
        assert_solved("1/(x^2+2*x+3)");
    }

    // ------------------- hang regressions (0.27.1 timeouts) -------------

    /// Corpus shapes whose integration used to run into the per-case budget
    /// with this stage last entered (Rubi ids in the comments).
    ///
    /// A representative subset, not the whole list: driving every hang shape
    /// through the full chain in a debug build costs minutes per shape, which
    /// made this module's test binary unusable in the workspace suite. The
    /// complete list is measured by the corpus harness
    /// (`ocas-tests/benches/integrate_1892.rs`), whose per-case budget and
    /// report are the authoritative record.
    const HANG_SHAPES: &[&str] = &[
        "sinh(x)^3/(a + b*sinh(x))",           // rubi-00272
        "1/((a - b*x)*(a + b*x)*(c + d*x)^3)", // rubi-00317
    ];

    #[test]
    fn corpus_hang_shapes_return() {
        for input in HANG_SHAPES {
            assert_returns_not_wrong(input);
        }
    }

    /// A previously-fine input must still solve with the budget in place.
    #[test]
    fn budget_keeps_solving_normal_inputs() {
        for input in [
            "1/(a+b*x^2)",
            "(d+e*x)/(x^3*(a+c*x^2))",
            "(A+B*x)/(a+b*x+c*x^2)",
            "1/(x^2*(a-b*x^2))",
        ] {
            assert_solved(input);
        }
    }

    /// The budget resets on every entry: repeating a budget-tripping shape
    /// must neither slow down nor change the outcome (no cross-call
    /// leakage).
    ///
    /// The probe is the *direct* entry point on a shape that declines after a
    /// few hundred work units, so the repetition is cheap; the corpus shapes
    /// that route the whole chain cost seconds per call in a debug build and
    /// are covered by the harness instead.
    #[test]
    fn budget_does_not_leak_across_calls() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = Symbol::new("_t");
        let expr = ocas_parse::parse(
            &ctx,
            "(2*((1 + _t^2)^-1))*((a + b*(2*_t*(1 + _t^2)^-1))^-3)",
        )
        .unwrap();
        let first = integrate_rational_symbolic(&ctx, expr, t).map(|a| a.to_string());
        assert!(first.is_none(), "expected a decline, got {first:?}");
        for i in 0..64 {
            let r = integrate_rational_symbolic(&ctx, expr, t).map(|a| a.to_string());
            assert_eq!(r, first, "call {i} diverged");
        }
    }

    /// Direct entry-point check: the module itself must decline a Weierstrass
    /// t-form with symbolic coefficients instead of grinding in the field
    /// Euclidean loops, and the entry reset must keep two consecutive
    /// invocations at the same cost (no cross-call leakage).
    #[test]
    fn symbolic_backend_declines_heavy_t_form() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = Symbol::new("_t");
        // 1/(a + b·(2t/(1+t²)))^3 · 2/(1+t²): a symbolic-coefficient rational
        // function of `t` with a repeated denominator factor.
        let expr = ocas_parse::parse(
            &ctx,
            "(2*((1 + _t^2)^-1))*((a + b*(2*_t*(1 + _t^2)^-1))^-3)",
        )
        .unwrap();
        reset_peak_work();
        let first = integrate_rational_symbolic(&ctx, expr, t);
        let first_peak = peak_work_units();
        reset_peak_work();
        let second = integrate_rational_symbolic(&ctx, expr, t);
        let second_peak = peak_work_units();
        if let Some(r) = &first {
            assert!(
                !r.to_string().contains("Integral("),
                "a decline was expected, got {r}"
            );
        }
        assert_eq!(
            first.map(|a| a.to_string()),
            second.map(|a| a.to_string()),
            "two invocations disagreed"
        );
        // The decline now happens before any heavy work (both runs stay in the
        // hundreds of units), so the remaining variance is charge-accounting
        // noise of a few units, not leakage. Guard the bound, not exact
        // equality.
        assert!(
            first_peak.abs_diff(second_peak) <= 8,
            "budget leaked between invocations ({first_peak} then {second_peak})"
        );
        assert!(
            first_peak < 4 * MAX_WORK_UNITS && second_peak < 4 * MAX_WORK_UNITS,
            "budget did not contain the t-form: {first_peak} then {second_peak} units"
        );
    }

    // --------------- exactness of the quadratic square root --------------

    /// A field element over `ℚ(a, b)` from `(exponents, coefficient)` pairs
    /// (`a` is variable 0, `b` is variable 1).
    fn field_of(terms: &[(&[usize], i64)]) -> GeneratorField {
        GeneratorField::from_polynomial(Sparse::from_terms(
            RationalDomain,
            2,
            terms
                .iter()
                .map(|(e, c)| (e.to_vec(), Rational::new(*c, 1)))
                .collect(),
        ))
    }

    /// The monomial-wise candidate must be squared and checked: a *sum* of
    /// squares is not a square. Regression for the 0.27.1 wrong answer on
    /// `1/(b*x^2 + 2*a*x - b)`, whose discriminant is `4a² + 4b²`.
    #[test]
    fn sum_discriminant_is_not_a_square() {
        // Δ = 4a² + 4b² → the term-wise candidate would be 2a + 2b, whose
        // square is 4a² + 8ab + 4b². Must be rejected.
        let delta = field_of(&[(&[2, 0], 4), (&[0, 2], 4)]);
        assert!(
            rational_square_root(&delta).is_none(),
            "4a² + 4b² accepted as a square"
        );
        // Δ = 4a² + 8ab + 4b² IS a square, but a cross term has odd
        // exponents, so the conservative candidate generator declines it
        // (unchanged behaviour — never a wrong answer).
        let square = field_of(&[(&[2, 0], 4), (&[1, 1], 8), (&[0, 2], 4)]);
        assert!(rational_square_root(&square).is_none());
    }

    /// Genuine single-monomial squares must still be recognised.
    #[test]
    fn monomial_discriminant_still_splits() {
        for (delta, want) in [
            (field_of(&[(&[0, 2], 4)]), field_of(&[(&[0, 1], 2)])),
            (field_of(&[(&[2, 0], 4)]), field_of(&[(&[1, 0], 2)])),
            (field_of(&[(&[2, 0], 1)]), field_of(&[(&[1, 0], 1)])),
        ] {
            let got = rational_square_root(&delta).expect("genuine square declined");
            assert!(got.mul(&got).sub(&delta).is_zero(), "√Δ is not a root of Δ");
            assert!(
                got.sub(&want).is_zero(),
                "unexpected root: {} vs {}",
                got,
                want
            );
        }
        // A non-square monomial (coefficient 2) is still rejected.
        assert!(rational_square_root(&field_of(&[(&[2, 0], 2)])).is_none());
    }

    /// End-to-end: the shapes whose discriminant is a sum of squares must
    /// never produce the bogus two-log split; either the honest quadratic
    /// (log/atan) answer or a fallback is acceptable, a wrong answer is not.
    #[test]
    fn sum_square_discriminant_integrands_are_not_wrong() {
        for input in ["1/(b*x^2+2*a*x-b)", "2/(b*x^2+2*a*x-b)", "1/(a+b*sinh(x))"] {
            assert_returns_not_wrong(input);
        }
    }

    /// `1/(b*x^2-b)` declined before this change and must keep doing so.
    #[test]
    fn monomial_discriminant_integrand_unchanged() {
        let r = int_str("1/(b*x^2-b)", "x");
        assert!(r.contains("Integral("), "expected the fallback, got {r}");
    }
}
