//! Elements of the differential fields in an extension tower.
//!
//! A tower field `ℚ(x, t₁, …, tₙ)` is represented in flat multivariate
//! form: an element of the field is a [`KElem`] (a numerator/denominator
//! pair of sparse polynomials over `ℚ`), and a univariate polynomial over
//! the top generator is a [`KPoly`] (dense coefficient vector in the top
//! variable). [`KRat`] is a rational function in the top generator.
//!
//! No multivariate GCD reduction is performed on [`KElem`]; only zero
//! detection (via the numerator) and cross-multiplied equality are
//! reliable. [`KPoly`] carries a genuine Euclidean algorithm because the
//! coefficient field is exact.

use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::{Lex, SparseMultivariatePolynomial};

pub(crate) type Sparse = SparseMultivariatePolynomial<RationalDomain, Lex>;

fn sparse_one(n: usize) -> Sparse {
    Sparse::new(RationalDomain, n).one()
}

/// If `num = c·den` as polynomials for a rational `c`, return `c`.
///
/// Works on unreduced fractions by comparing supports and per-term
/// coefficient ratios.
pub(crate) fn ratio_constant(num: &Sparse, den: &Sparse) -> Option<Rational> {
    let dom = RationalDomain;
    if num.is_zero() {
        return Some(dom.zero());
    }
    let nt = num.terms_ref();
    let dt = den.terms_ref();
    if nt.len() != dt.len() {
        return None;
    }
    let mut out = None;
    for (exp, nc) in nt {
        let dc = dt.get(exp)?;
        let ratio = dom.div(nc, dc)?;
        match &out {
            None => out = Some(ratio),
            Some(c) if *c == ratio => {}
            _ => return None,
        }
    }
    out
}

/// Cancel a shared monomial factor between `num` and `den`
/// (`a·tᵏ / (b·tᵏ) → a/b` for monomial-free leftovers).
///
/// Only handles the case where `den` is a single term; returns
/// `(num/den, one)` when every num term dominates the den exponent.
pub(crate) fn cancel_monomial(num: &Sparse, den: &Sparse) -> Option<(Sparse, Sparse)> {
    let dom = RationalDomain;
    if den.terms_ref().len() != 1 {
        return None;
    }
    let (d_exp, d_coeff) = den.terms_ref().iter().next().unwrap();
    let mut out_terms = Vec::new();
    for (n_exp, n_coeff) in num.terms_ref() {
        if n_exp.iter().zip(d_exp.iter()).any(|(a, b)| a < b) {
            return None;
        }
        let new_exp: Vec<usize> = n_exp.iter().zip(d_exp.iter()).map(|(a, b)| a - b).collect();
        out_terms.push((new_exp, dom.div(n_coeff, d_coeff)?));
    }
    let n = Sparse::from_terms(RationalDomain, num.n_vars(), out_terms);
    let one = num.one();
    Some((n, one))
}

/// Formal partial derivative of a sparse polynomial w.r.t. variable `var`.
pub(crate) fn sparse_partial_deriv(p: &Sparse, var: usize) -> Sparse {
    let dom = RationalDomain;
    let mut terms = Vec::new();
    for (exp, c) in p.terms_ref() {
        let e = exp[var];
        if e > 0 {
            let mut ne = exp.to_vec();
            ne[var] -= 1;
            terms.push((ne, dom.mul(c, &Rational::new(e as i64, 1))));
        }
    }
    Sparse::from_terms(RationalDomain, p.n_vars(), terms)
}

// ------------------------------------------------------------------
//  KElem: element of k = ℚ(x, t₁, …, t_{n-1})
// ------------------------------------------------------------------

/// An element of the coefficient field as a multivariate rational function.
#[derive(Clone, Debug)]
pub(crate) struct KElem {
    /// Numerator.
    pub num: Sparse,
    /// Denominator (never zero).
    pub den: Sparse,
}

impl KElem {
    /// Create without reduction. `den` must be non-zero.
    pub fn new(num: Sparse, den: Sparse) -> Self {
        debug_assert!(!den.is_zero());
        Self { num, den }
    }

    /// Create and reduce by cancelling a constant-support ratio
    /// (`num = c·den`) or a shared monomial factor (`a·tᵏ / b·tᵏ`).
    ///
    /// This is a best-effort normalizer: without a multivariate GCD the
    /// fraction is not fully reduced, but the common cases produced by
    /// tower arithmetic (quotients like `t/t`, cross-multiplied sums)
    /// collapse to their minimal form.
    pub fn reduced(num: Sparse, den: Sparse) -> Self {
        debug_assert!(!den.is_zero());
        if let Some(c) = ratio_constant(&num, &den) {
            let one = den.one();
            let const_poly = Sparse::from_terms(
                RationalDomain,
                den.n_vars(),
                vec![(vec![0; den.n_vars()], c)],
            );
            return Self {
                num: const_poly,
                den: one,
            };
        }
        if let Some((n2, d2)) = cancel_monomial(&num, &den) {
            return Self { num: n2, den: d2 };
        }
        Self { num, den }
    }

    pub fn zero(n: usize) -> Self {
        Self::new(Sparse::new(RationalDomain, n), sparse_one(n))
    }

    pub fn one(n: usize) -> Self {
        Self::new(sparse_one(n), sparse_one(n))
    }

    pub fn from_poly(p: Sparse) -> Self {
        let one = p.one();
        Self::new(p, one)
    }

    pub fn from_rational(r: &Rational, n: usize) -> Self {
        Self::from_poly(Sparse::from_terms(
            RationalDomain,
            n,
            vec![(vec![0; n], r.clone())],
        ))
    }

    /// The variable with the given index as a field element.
    pub fn var(idx: usize, n: usize) -> Self {
        let mut e = vec![0; n];
        e[idx] = 1;
        Self::from_poly(Sparse::from_terms(
            RationalDomain,
            n,
            vec![(e, RationalDomain.one())],
        ))
    }

    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    pub fn n_vars(&self) -> usize {
        self.num.n_vars()
    }

    /// Cross-multiplied equality (reliable without GCD reduction).
    pub fn eq_cross(&self, other: &Self) -> bool {
        self.num.mul(&other.den) == other.num.mul(&self.den)
    }

    /// If this element equals a rational constant `c` (i.e. `num = c·den`
    /// as polynomials), return it. Works on unreduced fractions by
    /// comparing supports and per-term ratios.
    pub fn as_rational(&self) -> Option<Rational> {
        ratio_constant(&self.num, &self.den)
    }

    pub fn neg(&self) -> Self {
        Self::new(self.num.neg(), self.den.clone())
    }

    pub fn add(&self, o: &Self) -> Self {
        if self.is_zero() {
            return o.clone();
        }
        if o.is_zero() {
            return self.clone();
        }
        if self.den == o.den {
            return Self::reduced(self.num.add(&o.num), self.den.clone());
        }
        Self::reduced(
            self.num.mul(&o.den).add(&o.num.mul(&self.den)),
            self.den.mul(&o.den),
        )
    }

    pub fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }

    pub fn mul(&self, o: &Self) -> Self {
        if self.is_zero() || o.is_zero() {
            return Self::zero(self.n_vars());
        }
        Self::reduced(self.num.mul(&o.num), self.den.mul(&o.den))
    }

    pub fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self::new(self.den.clone(), self.num.clone()))
        }
    }

    pub fn div(&self, o: &Self) -> Option<Self> {
        Some(self.mul(&o.inv()?))
    }

    pub fn mul_rational(&self, r: &Rational) -> Self {
        Self::new(self.num.mul_scalar(r), self.den.clone())
    }

    #[allow(dead_code)] // used by unit tests
    pub fn pow(&self, mut k: u64) -> Self {
        let mut result = Self::one(self.n_vars());
        let mut base = self.clone();
        while k > 0 {
            if k & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            k >>= 1;
        }
        result
    }

    /// Formal partial derivative w.r.t. variable `var` (quotient rule).
    pub fn partial_deriv(&self, var: usize) -> Self {
        let a = &self.num;
        let b = &self.den;
        let av = sparse_partial_deriv(a, var);
        let bv = sparse_partial_deriv(b, var);
        Self::new(av.mul(b).add(&a.mul(&bv).neg()), b.mul(b))
    }
}

// ------------------------------------------------------------------
//  KPoly: univariate polynomial over k in the top variable
// ------------------------------------------------------------------

/// A dense univariate polynomial over the coefficient field `k`.
#[derive(Clone, Debug)]
pub(crate) struct KPoly {
    /// Index of the polynomial variable (the top generator of the level).
    pub top: usize,
    /// Coefficients by ascending degree; trailing zeros are trimmed.
    pub coeffs: Vec<KElem>,
    /// Total number of variables of the underlying multivariate ring.
    pub n_vars: usize,
}

impl KPoly {
    fn trim(&mut self) {
        while let Some(last) = self.coeffs.last() {
            if last.is_zero() {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }

    pub fn zero(top: usize, n: usize) -> Self {
        Self {
            top,
            coeffs: Vec::new(),
            n_vars: n,
        }
    }

    pub fn one(top: usize, n: usize) -> Self {
        Self {
            top,
            coeffs: vec![KElem::one(n)],
            n_vars: n,
        }
    }

    /// A constant polynomial (degree 0 in the top variable).
    pub fn from_kelem(c: KElem, top: usize) -> Self {
        let n = c.n_vars();
        if c.is_zero() {
            return Self::zero(top, n);
        }
        Self {
            top,
            coeffs: vec![c],
            n_vars: n,
        }
    }

    /// View a sparse polynomial as a univariate polynomial in `top`.
    pub fn from_sparse(p: &Sparse, top: usize) -> Self {
        let n = p.n_vars();
        let deg = p.degree_in(top);
        let mut buckets: Vec<Vec<(Vec<usize>, Rational)>> = vec![Vec::new(); deg + 1];
        for (exp, c) in p.terms_ref() {
            let mut e = exp.to_vec();
            let d = e[top];
            e[top] = 0;
            buckets[d].push((e, c.clone()));
        }
        let coeffs = buckets
            .into_iter()
            .map(|terms| KElem::from_poly(Sparse::from_terms(RationalDomain, n, terms)))
            .collect();
        let mut r = Self {
            top,
            coeffs,
            n_vars: n,
        };
        r.trim();
        r
    }

    /// Collapse to a single field element via Horner-style summation.
    pub fn kelem(&self) -> KElem {
        let mut acc = KElem::zero(self.n_vars);
        let t = KElem::var(self.top, self.n_vars);
        let mut tp = KElem::one(self.n_vars);
        for c in &self.coeffs {
            acc = acc.add(&c.mul(&tp));
            tp = tp.mul(&t);
        }
        acc
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    pub fn lc(&self) -> KElem {
        self.coeffs
            .last()
            .cloned()
            .unwrap_or_else(|| KElem::zero(self.n_vars))
    }

    pub fn is_one(&self) -> bool {
        self.degree() == Some(0) && self.coeffs[0].eq_cross(&KElem::one(self.n_vars))
    }

    /// Coefficient of `t^i`, or zero.
    pub fn coeff_at(&self, i: usize) -> KElem {
        self.coeffs
            .get(i)
            .cloned()
            .unwrap_or_else(|| KElem::zero(self.n_vars))
    }

    pub fn neg(&self) -> Self {
        Self {
            top: self.top,
            coeffs: self.coeffs.iter().map(|c| c.neg()).collect(),
            n_vars: self.n_vars,
        }
    }

    pub fn add(&self, o: &Self) -> Self {
        debug_assert_eq!(self.top, o.top);
        let len = self.coeffs.len().max(o.coeffs.len());
        let mut coeffs = Vec::with_capacity(len);
        for i in 0..len {
            coeffs.push(self.coeff_at(i).add(&o.coeff_at(i)));
        }
        let mut r = Self {
            top: self.top,
            coeffs,
            n_vars: self.n_vars,
        };
        r.trim();
        r
    }

    pub fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }

    pub fn mul(&self, o: &Self) -> Self {
        debug_assert_eq!(self.top, o.top);
        if self.is_zero() || o.is_zero() {
            return Self::zero(self.top, self.n_vars);
        }
        let mut coeffs = vec![KElem::zero(self.n_vars); self.coeffs.len() + o.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in o.coeffs.iter().enumerate() {
                coeffs[i + j] = coeffs[i + j].add(&a.mul(b));
            }
        }
        let mut r = Self {
            top: self.top,
            coeffs,
            n_vars: self.n_vars,
        };
        r.trim();
        r
    }

    pub fn mul_kelem(&self, c: &KElem) -> Self {
        let mut r = Self {
            top: self.top,
            coeffs: self.coeffs.iter().map(|a| a.mul(c)).collect(),
            n_vars: self.n_vars,
        };
        r.trim();
        r
    }

    /// `c · t^k · self`.
    pub fn monomial_shift(&self, c: &KElem, k: usize) -> Self {
        if self.is_zero() || c.is_zero() {
            return Self::zero(self.top, self.n_vars);
        }
        let mut coeffs = vec![KElem::zero(self.n_vars); k];
        coeffs.extend(self.coeffs.iter().map(|a| a.mul(c)));
        Self {
            top: self.top,
            coeffs,
            n_vars: self.n_vars,
        }
    }

    /// Formal derivative w.r.t. the polynomial variable `t`.
    pub fn derivative_dt(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Self::zero(self.top, self.n_vars);
        }
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c.mul_rational(&Rational::new(i as i64, 1)))
            .collect();
        Self {
            top: self.top,
            coeffs,
            n_vars: self.n_vars,
        }
    }

    /// Exact division over the field `k`: returns `(quotient, remainder)`.
    pub fn div_rem(&self, d: &Self) -> (Self, Self) {
        assert!(!d.is_zero(), "KPoly::div_rem: division by zero");
        let mut q = Self::zero(self.top, self.n_vars);
        let mut r = self.clone();
        let d_lc = d.lc();
        let d_deg = match d.degree() {
            Some(v) => v,
            None => unreachable!("non-zero polynomial has a degree"),
        };
        while let Some(r_deg) = r.degree() {
            if r_deg < d_deg {
                break;
            }
            let k = r_deg - d_deg;
            let c = r.lc().div(&d_lc).expect("field division");
            // q += c·t^k
            if q.coeffs.len() <= k {
                q.coeffs.resize(k + 1, KElem::zero(self.n_vars));
            }
            q.coeffs[k] = q.coeffs[k].add(&c);
            r = r.sub(&d.monomial_shift(&c, k));
        }
        q.trim();
        (q, r)
    }

    /// Monic GCD over `k[t]` (Euclidean algorithm).
    pub fn gcd(&self, o: &Self) -> Self {
        let mut a = self.clone();
        let mut b = o.clone();
        while !b.is_zero() {
            let (_, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
        if a.is_zero() {
            return a;
        }
        a.monic()
    }

    /// Divide by the leading coefficient (no-op for zero).
    pub fn monic(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        let inv = self.lc().inv().expect("non-zero leading coefficient");
        self.mul_kelem(&inv)
    }

    /// Extended GCD: returns `(g, s, t)` with `s·self + t·other = g` monic.
    pub fn eea(&self, other: &Self) -> (Self, Self, Self) {
        let mut old_r = self.clone();
        let mut r = other.clone();
        let mut old_s = Self::one(self.top, self.n_vars);
        let mut s = Self::zero(self.top, self.n_vars);
        let mut old_t = Self::zero(self.top, self.n_vars);
        let mut t = Self::one(self.top, self.n_vars);
        while !r.is_zero() {
            let (q, rem) = old_r.div_rem(&r);
            old_r = r;
            r = rem;
            let new_s = old_s.sub(&q.mul(&s));
            old_s = s;
            s = new_s;
            let new_t = old_t.sub(&q.mul(&t));
            old_t = t;
            t = new_t;
        }
        if !old_r.is_zero() {
            let inv = old_r.lc().inv().expect("non-zero leading coefficient");
            old_r = old_r.mul_kelem(&inv);
            old_s = old_s.mul_kelem(&inv);
            old_t = old_t.mul_kelem(&inv);
        }
        (old_r, old_s, old_t)
    }

    /// Yun's square-free factorization over `k[t]`.
    pub fn square_free(&self) -> Vec<(Self, usize)> {
        let mut out = Vec::new();
        if self.is_zero() {
            return out;
        }
        let f = self.monic();
        let df = f.derivative_dt();
        let mut g = f.gcd(&df);
        if g.is_zero() {
            return out;
        }
        let (mut w, _) = f.div_rem(&g);
        let mut k = 1;
        while !w.is_one() && !w.is_zero() {
            let h = w.gcd(&g);
            let (z, _) = w.div_rem(&h);
            if !z.is_one() && !z.is_zero() {
                out.push((z, k));
            }
            let (new_g, _) = g.div_rem(&h);
            g = new_g;
            w = h;
            k += 1;
        }
        out
    }
}

// ------------------------------------------------------------------
//  KRat: rational function in the top variable
// ------------------------------------------------------------------

/// A rational function `num / den` in the top generator, kept coprime with
/// a monic denominator.
#[derive(Clone, Debug)]
pub(crate) struct KRat {
    /// Numerator.
    pub num: KPoly,
    /// Denominator (monic, non-zero).
    pub den: KPoly,
}

impl KRat {
    /// Create and reduce to coprime numerator/denominator.
    pub fn new(num: KPoly, den: KPoly) -> Self {
        assert!(!den.is_zero(), "KRat: zero denominator");
        if num.is_zero() {
            let one = KPoly::one(den.top, den.n_vars);
            return Self { num, den: one };
        }
        let g = num.gcd(&den);
        let (n, _) = num.div_rem(&g);
        let (d, _) = den.div_rem(&g);
        let lc_inv = d.lc().inv().expect("non-zero denominator lc");
        Self {
            num: n.mul_kelem(&lc_inv),
            den: d.mul_kelem(&lc_inv),
        }
    }

    /// Collapse to a single coefficient-field element.
    pub fn kelem(&self) -> KElem {
        self.num
            .kelem()
            .div(&self.den.kelem())
            .expect("non-zero denominator")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rat(p: i64, q: i64) -> Rational {
        Rational::new(p, q)
    }

    /// (t - c) as a KPoly with top variable `top` and constant c.
    fn linear(top: usize, n: usize, c: i64) -> KPoly {
        KPoly {
            top,
            coeffs: vec![KElem::from_rational(&rat(-c, 1), n), KElem::one(n)],
            n_vars: n,
        }
    }

    #[test]
    fn sparse_mul_exp_add() {
        // (x + 1)(x + 2) = x^2 + 3x + 2 in n_vars=2 (y unused).
        let dom = RationalDomain;
        let p = Sparse::from_terms(
            dom,
            2,
            vec![(vec![1, 0], rat(1, 1)), (vec![0, 0], rat(1, 1))],
        );
        let q = Sparse::from_terms(
            dom,
            2,
            vec![(vec![1, 0], rat(1, 1)), (vec![0, 0], rat(2, 1))],
        );
        let r = p.mul(&q);
        assert_eq!(r.n_terms(), 3);
        assert_eq!(r.coeff(&[2, 0]), rat(1, 1));
        assert_eq!(r.coeff(&[1, 0]), rat(3, 1));
        assert_eq!(r.coeff(&[0, 0]), rat(2, 1));
    }

    #[test]
    fn kelem_arithmetic() {
        let n = 2;
        let x = KElem::var(0, n);
        let t = KElem::var(1, n);
        // (x/t) + (t/x) = (x^2 + t^2)/(x·t)
        let sum = x.div(&t).unwrap().add(&t.div(&x).unwrap());
        let expect_num = x.mul(&x).add(&t.mul(&t));
        let expect = expect_num.div(&x.mul(&t)).unwrap();
        assert!(sum.eq_cross(&expect));
        // inv
        assert!(
            x.div(&t)
                .unwrap()
                .mul(&t.div(&x).unwrap())
                .eq_cross(&KElem::one(n))
        );
    }

    #[test]
    fn kelem_partial_deriv() {
        let n = 2;
        // f = t^2 / x; df/dx = -t^2/x^2; df/dt = 2t/x
        let x = KElem::var(0, n);
        let t = KElem::var(1, n);
        let f = t.pow(2).div(&x).unwrap();
        let dfdx = f.partial_deriv(0);
        assert!(dfdx.eq_cross(&t.pow(2).div(&x.pow(2)).unwrap().neg()));
        let dfdt = f.partial_deriv(1);
        let two_t = t.mul_rational(&rat(2, 1));
        assert!(dfdt.eq_cross(&two_t.div(&x).unwrap()));
    }

    #[test]
    fn kpoly_gcd_clean_exponents() {
        // gcd(t + x, 1) over n_vars=2, top=1: monic 1, no garbage exponents.
        let n = 2;
        let a = KPoly {
            top: 1,
            coeffs: vec![KElem::var(0, n), KElem::one(n)],
            n_vars: n,
        };
        let b = KPoly::one(1, n);
        let g = a.gcd(&b);
        assert!(g.is_one());
        for c in &g.coeffs {
            for e in c.num.terms_ref().keys() {
                assert!(e.len() <= n, "garbage exponent {e:?} in num");
            }
        }
    }

    #[test]
    fn kpoly_div_rem() {
        // (t^2 - 1) / (t - 1) = (t + 1), remainder 0
        let n = 2;
        let a = linear(1, n, 1).mul(&linear(1, n, -1));
        let b = linear(1, n, 1);
        let (q, r) = a.div_rem(&b);
        assert!(r.is_zero());
        assert!(q.kelem().eq_cross(&linear(1, n, -1).kelem()));
    }

    #[test]
    fn kpoly_gcd() {
        // gcd(t^2 - 1, t^2 - 2t + 1) = t - 1
        let n = 2;
        let a = linear(1, n, 1).mul(&linear(1, n, -1));
        let b = linear(1, n, 1).mul(&linear(1, n, 1));
        let g = a.gcd(&b);
        assert!(g.kelem().eq_cross(&linear(1, n, 1).kelem()));
    }

    #[test]
    fn kpoly_square_free() {
        // (t-1)^2 (t-2) → [(t-1, 2), (t-2, 1)]
        let n = 2;
        let f = linear(1, n, 1).mul(&linear(1, n, 1)).mul(&linear(1, n, 2));
        let sf = f.square_free();
        assert_eq!(sf.len(), 2);
        let mut by_mult: Vec<(usize, KElem)> = sf.iter().map(|(p, k)| (*k, p.kelem())).collect();
        by_mult.sort_by_key(|(k, _)| *k);
        assert_eq!(by_mult[0].0, 1);
        assert!(by_mult[0].1.eq_cross(&linear(1, n, 2).kelem()));
        assert_eq!(by_mult[1].0, 2);
        assert!(by_mult[1].1.eq_cross(&linear(1, n, 1).kelem()));
    }

    #[test]
    fn kpoly_eea() {
        // s·(t-1) + t·(t+1) = 1
        let n = 2;
        let a = linear(1, n, 1);
        let b = linear(1, n, -1);
        let (g, s, tt) = a.eea(&b);
        assert!(g.is_one());
        let check = s.mul(&a).add(&tt.mul(&b));
        assert!(check.is_one());
    }

    #[test]
    fn from_sparse_roundtrip() {
        // p = x·t + t^2 over n=2, top=1
        let n = 2;
        let dom = RationalDomain;
        let p = Sparse::from_terms(
            dom,
            n,
            vec![(vec![1, 1], rat(1, 1)), (vec![0, 2], rat(3, 1))],
        );
        let kp = KPoly::from_sparse(&p, 1);
        assert_eq!(kp.degree(), Some(2));
        // coeffs: [0, x, 3]
        assert!(kp.coeffs[1].eq_cross(&KElem::var(0, n)));
        assert!(kp.coeffs[2].eq_cross(&KElem::from_rational(&rat(3, 1), n)));
    }
}
