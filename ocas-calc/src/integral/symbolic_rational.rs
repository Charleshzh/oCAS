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
                Some((_, acc)) => *acc = acc.add(c),
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
                let prod = mono_reduce(c.mul(d));
                match out.terms.iter_mut().find(|(r, _)| *r == pow) {
                    Some((_, acc)) => *acc = mono_reduce(acc.add(&prod)),
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
                    mono_reduce(c.mul(&rat_const(*p as i64, self.n_vars()))),
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
    fn div_rem(&self, den: &Self) -> Option<(Self, Self)> {
        let mut q = fpoly_zero(self.n_vars());
        let mut r = self.clone();
        let dd = den.degree()?;
        let lc_d = den.leading_coeff()?;
        while let Some(dr) = r.degree() {
            if dr < dd {
                break;
            }
            let lc_r = r.leading_coeff()?;
            let c = mono_reduce(lc_r.div(&lc_d)?);
            let t = FPoly {
                terms: vec![(dr - dd, c)],
            };
            q = q.add(&t);
            r = r.sub(&den.mul(&t));
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
            let gap = prev.map(|p| p - *exp).unwrap_or(0);
            for _ in 0..=gap {
                acc = acc.mul(v);
            }
            acc = acc.add(c);
            prev = Some(*exp);
        }
        if let Some(p) = prev {
            for _ in 0..p {
                acc = acc.mul(v);
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
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut steps = 0usize;
    while !r.is_zero() {
        steps += 1;
        if steps > 512 {
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

impl FPoly {
    fn scale(&self, c: &GeneratorField) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(p, q)| (*p, mono_reduce(q.mul(c))))
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
    while b.degree() != Some(0) {
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
    loop {
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
    while !r.is_zero() {
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
    four.mul(a).mul(c).sub(&b.mul(b))
}

/// `√Δ` as a rational function of the symbols when Δ is a square in
/// `ℚ(symbols)` (all monomial exponents even, constant coefficient a
/// rational square); else `None`.
fn rational_square_root(delta: &GeneratorField) -> Option<GeneratorField> {
    let sqrt_sparse = |p: &Sparse| -> Option<Sparse> {
        let mut terms: Vec<(Vec<usize>, Rational)> = Vec::new();
        for (e, c) in p.terms_ref() {
            if e.iter().any(|&v| v % 2 != 0) {
                return None;
            }
            let rp = isqrt_i64(c.numer().to_i64()?)?;
            terms.push((e.iter().map(|v| v / 2).collect(), Rational::new(rp, 1)));
        }
        Some(Sparse::from_terms(RationalDomain, p.n_vars(), terms))
    };
    let n = sqrt_sparse(&delta.numerator)?;
    let d = sqrt_sparse(&delta.denominator)?;
    Some(GeneratorField::from_num_den(n, d))
}

/// Integer multivariate factorization of a square-free factor: clears the
/// coefficient denominators, factors over ℤ, and returns the factors (each
/// converted back to the field) — their product equals `scalar·f` where the
/// scalar is the clearing constant (the log part multiplies its numerator
/// by that scalar, keeping the partial-fraction coefficients exact).
fn factor_via_integer(f: &FPoly) -> Option<Vec<(FPoly, usize)>> {
    let sparse = f.to_sparse();
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
fn field_divisors(el: &GeneratorField) -> Option<Vec<GeneratorField>> {
    let mut out: Vec<GeneratorField> = Vec::new();
    for e in el.numerator.terms_ref().keys() {
        let mut exps: Vec<Vec<usize>> = vec![vec![0; e.len()]];
        for (i, &v) in e.iter().enumerate() {
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
        for da in field_divisors(&a0)? {
            for dl in field_divisors(&lc)? {
                if !dl.is_zero() {
                    let r = da.div(&dl)?;
                    candidates.push(r.neg());
                    candidates.push(r);
                }
            }
        }
    }
    for r in candidates {
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
            let delta = b.mul(&b).sub(&rat_const(4, f.n_vars()).mul(&a).mul(&c));
            if let Some(s) = rational_square_root(&delta) {
                let two_a = a.mul(&rat_const(2, f.n_vars()));
                let r1 = b.neg().add(&s).div(&two_a)?;
                let r2 = b.neg().sub(&s).div(&two_a)?;
                let one = GeneratorField::one(&RationalDomain, f.n_vars());
                // Keep the leading coefficient of the original quadratic in
                // the split (product must equal f): the log-part residues
                // evaluate ∏_{j≠i} f_j(r)/f_i'(r), which is only correct
                // when the factors multiply back to the exact denominator.
                let mut g1 = FPoly {
                    terms: vec![(1, a.clone()), (0, a.mul(&r1).neg())],
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
    let mut num = FPoly::from_sparse(&rf.numerator);
    let mut den = FPoly::from_sparse(&rf.denominator);
    if den.degree()? > 6 {
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
        let f_deg = f.degree()?;
        if f_deg == 1 {
            // c = num(r)/(f'(r)·∏_{j≠i} f_j(r)) with r = −β/α.
            let (alpha, beta) = linear_coeffs(f)?;
            let r = beta.neg().div(&alpha)?;
            let mut denom = alpha.clone();
            for (g, _) in &factors {
                if g.terms == f.terms {
                    continue;
                }
                denom = denom.mul(&g.eval(&r)?);
            }
            let coeff = c_num.eval(&r)?.div(&denom)?;
            let coeff_atom = rational_to_atom(ctx, &coeff, &gens[1..])?;
            let f_atom = f.to_atom(ctx, &gens)?;
            parts.push(ctx.mul(&[coeff_atom, ctx.fun("log", &[f_atom])]));
        } else {
            // Quadratic factor: solve (M·x+N)/(Q·f) partial fraction.
            let (m, n) = quadratic_coeffs(&c_num, f, &factors)?;
            let (a, b, c) = quadratic_coeffs_fpoly(f)?;
            // ∫(Mx+N)/f = M/(2a)·log(f) + (N − M·b/(2a))·(2/√Δ)·h((2a·x+b)/√Δ)
            let two_a = a.mul(&rat_const(2, n_vars));
            let m_over = m.div(&two_a)?;
            let f_atom = f.to_atom(ctx, &gens)?;
            let m_atom = rational_to_atom(ctx, &m_over, &gens[1..])?;
            if !m_over.is_zero() {
                parts.push(ctx.mul(&[m_atom, ctx.fun("log", &[f_atom])]));
            }
            let delta = discriminant(&a, &b, &c);
            let mb_over = m.mul(&b).div(&two_a)?;
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
    assemble(ctx, parts)
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
    let m11 = q0.sub(&q1.mul(&b).mul(&a_inv));
    let m21 = q1.mul(&c).mul(&a_inv).neg();
    let det = m11.mul(&q0).sub(&q1.mul(&m21));
    let det_inv = det.inv()?;
    let m = p1.mul(&q0).sub(&q1.mul(&p0)).mul(&det_inv);
    let n = m11.mul(&p0).sub(&p1.mul(&m21)).mul(&det_inv);
    Some((m, n))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

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
}
