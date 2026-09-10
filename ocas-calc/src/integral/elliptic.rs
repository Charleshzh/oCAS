//! Elliptic-integral reduction (0.27.2 Wave D).
//!
//! Integrates `∫ R(u, y) du` with `y² = Q(u)`, `deg Q = 4` and `Q` squarefree,
//! by Legendre reduction into the symmetric forms
//!
//! - `EllipticF(φ, m)` — incomplete elliptic integral of the first kind,
//! - `EllipticE(φ, m)` — second kind,
//! - `EllipticPi(n, φ, m)` — third kind,
//!
//! with the SymPy parameter convention `m = k²` and `φ` the amplitude.
//! Argument order is semantic: `ocas_atom::normalize::preserves_argument_order`
//! keeps these heads (and `Ei`) in the order they are built, so the emitted
//! `EllipticF(φ, m)` never has its amplitude and parameter swapped.
//!
//! # Reduction outline
//!
//! 1. The integrand is normalised to `W(u)/√(Q(u))` with `W` a rational
//!    function and `Q` a polynomial. A half-power `b^{p/2}` with `p` odd
//!    contributes `b^{(p−1)/2}` to `W` and `b` to the radicand.
//! 2. `W` is converted to the Hermite form `N(u)/Q(u)^{3/2}` — possible exactly
//!    when the denominator of `W` divides `Q`, which holds for every integrand
//!    whose only pole comes from the radical itself.
//! 3. `N` is divided by `Q`; the quotient `S` is integrated directly and the
//!    remainder `R` (`deg R < 4`) is Hermite-reduced to `A(u)/√Q + B(u)/√Q`
//!    through the coefficient identity `R = A'Q − A Q'/2 + B Q`.
//!    [`hermite_even`] solves it in closed form; because `Q` is even the
//!    system is triangular (see the derivation there).
//! 4. `S + B` is integrated against `1/√Q`. `Q` must already be in **Legendre
//!    normal form** `c(1 − u²)(1 − M u²)`, i.e. `q₀ + q₂ + q₄ ≡ 0`; then
//!    `∫du/√Q = F(asin u, M)/√c` and `∫u²du/√Q = (F − E)/(M√c)`.
//!    Odd powers of `u` integrate to elementary functions and are declined.
//!
//! # Symbolic arithmetic
//!
//! All polynomial arithmetic is done over [`Coeff`], a canonical form for
//! symbolic coefficients: a rational linear combination of monomials over
//! opaque atom bases. Canonicalisation makes `(a+b)²`, `a²+2ab+b²` and
//! `(a+b)·M` vs `Ma+Mb` the same object and makes `X/X = 1`, so the exact
//! polynomial divisions above really are exact for symbolic coefficients.
//! Without it, coefficient comparisons built from raw atoms are unreliable and
//! the reduction would silently mis-decline or mis-solve.
//!
//! # Amplitude / branch caveat
//!
//! The amplitude is emitted as `asin(u)` with `u` the radical's variable. When
//! [`super::halfpower`] reached that variable through `u = sin θ`, the caller
//! applies the corresponding branch correction (see that module); the engine
//! itself only claims the principal branch `|θ| ≤ π/2`.
//!
//! # Known gaps (declined, never guessed)
//!
//! - `deg Q = 3` (cubic radicands such as `1/√(x³+1)`): needs the cubic
//!   `F`/`E` transform or the `u = 1/(v − c)` reduction; not implemented.
//! - Radicands not already in Legendre normal form: the general quartic needs
//!   a Möbius normalisation through the cross-ratio of its four roots. Only
//!   the normal form is handled.
//! - Odd numerators against an even `Q` (elementary `asinh`/`log` results).
//! - Third kind: `1/((u − c)√Q)` with `c` not a branch point — the `Π`
//!   reduction is not implemented (`EllipticPi` is emitted by no engine here,
//!   though the numeric oracle supports the head).
//! - `∫ S(u)√Q du` with `deg S > 0` (a non-constant polynomial times the
//!   radical).
//!
//! Every emitted form is checked by the module's numeric tests: the
//! antiderivative is differentiated symbolically (a local differentiator that
//! knows the elliptic derivatives) and compared against the integrand at fixed
//! sample points, with the elliptic values produced by the Carlson symmetric
//! forms `RF`/`RD`/`RJ` implemented in [`testnum`].

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::is_constant;

// =========================================================================
// Atom helpers
// =========================================================================

/// `2^{-1}`, the denominator of a half-integer exponent.
fn half<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    ctx.pow(ctx.num(2), ctx.num(-1))
}

/// The atom `p/2`.
fn half_exp<'a>(ctx: &'a AtomArena<'a>, p: i64) -> Atom<'a> {
    ctx.mul(&[ctx.num(p), half(ctx)])
}

/// Parse a rational exponent atom `p/q` (small integers, `q > 0`).
///
/// Accepts the shapes produced by the parser for `p/q`: `Num(p)`, `q^{-1}`,
/// and `p · q^{-1}` in either factor order (with `1` already folded out by
/// `normalize`).
fn exp_fraction<'a>(exp: Atom<'a>) -> Option<(i64, i64)> {
    match exp.node() {
        AtomNode::Num(n) => Some((*n, 1)),
        AtomNode::Pow(b, e) => {
            if let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
                && *ee == -1
                && *bb > 0
            {
                return Some((1, *bb));
            }
            None
        }
        AtomNode::Mul(args) => {
            let mut num: Option<i64> = None;
            let mut den: Option<i64> = None;
            for a in args.iter() {
                match a.node() {
                    AtomNode::Num(n) => {
                        if num.is_some() {
                            return None;
                        }
                        num = Some(*n);
                    }
                    AtomNode::Pow(b, e) => {
                        if let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
                            && *ee == -1
                            && *bb > 0
                        {
                            if den.is_some() {
                                return None;
                            }
                            den = Some(*bb);
                        } else {
                            return None;
                        }
                    }
                    _ => return None,
                }
            }
            match (num, den) {
                (Some(p), Some(q)) => Some((p, q)),
                (Some(p), None) => Some((p, 1)),
                (None, Some(q)) => Some((1, q)),
                (None, None) => None,
            }
        }
        _ => None,
    }
}

// =========================================================================
// Canonical symbolic coefficients
// =========================================================================

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

/// Reduce `(n, d)` to lowest terms with `d > 0`.
fn reduce_rat(n: i128, d: i128) -> (i128, i128) {
    if d == 0 {
        return (n, d);
    }
    let (n, d) = if d < 0 { (-n, -d) } else { (n, d) };
    let g = gcd_i128(n, d);
    (n / g, d / g)
}

/// One monomial: a rational coefficient times a product of `base^exp` factors.
///
/// Bases are *opaque* atoms — a variable, a function application, or a power
/// that cannot be expanded (a fractional power) — so that every sum is a sum
/// of monomials and equality is structural.
#[derive(Clone, PartialEq)]
struct Term<'a> {
    key: String,
    c: (i128, i128),
    f: Vec<(Atom<'a>, i64)>,
}

/// A canonical symbolic coefficient: a rational linear combination of
/// monomials. Two mathematically equal coefficients over the same bases have
/// identical normal forms.
#[derive(Clone, PartialEq, Default)]
struct Coeff<'a> {
    t: Vec<Term<'a>>,
}

fn factor_key<'a>(f: &[(Atom<'a>, i64)]) -> String {
    let mut parts: Vec<String> = f.iter().map(|(b, e)| format!("{b}^{e}")).collect();
    parts.sort();
    parts.join("*")
}

impl<'a> Coeff<'a> {
    fn zero() -> Self {
        Self { t: Vec::new() }
    }

    fn one() -> Self {
        Self {
            t: vec![Term {
                key: String::new(),
                c: (1, 1),
                f: Vec::new(),
            }],
        }
    }

    fn num(n: i64) -> Self {
        if n == 0 {
            return Self::zero();
        }
        Self {
            t: vec![Term {
                key: String::new(),
                c: (n as i128, 1),
                f: Vec::new(),
            }],
        }
    }

    /// A single opaque base with exponent one.
    fn base(a: Atom<'a>) -> Self {
        let f = vec![(a, 1i64)];
        Self {
            t: vec![Term {
                key: factor_key(&f),
                c: (1, 1),
                f,
            }],
        }
    }

    fn is_zero(&self) -> bool {
        self.t.is_empty()
    }

    fn neg(&self) -> Self {
        Self {
            t: self
                .t
                .iter()
                .map(|m| Term {
                    key: m.key.clone(),
                    c: (-m.c.0, m.c.1),
                    f: m.f.clone(),
                })
                .collect(),
        }
    }

    fn add(&self, other: &Self) -> Self {
        let mut out: Vec<Term<'a>> = Vec::with_capacity(self.t.len() + other.t.len());
        let (mut i, mut j) = (0usize, 0usize);
        while i < self.t.len() || j < other.t.len() {
            let take_left = match (self.t.get(i), other.t.get(j)) {
                (Some(a), Some(b)) => a.key <= b.key,
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (None, None) => break,
            };
            if take_left {
                let a = &self.t[i];
                if let Some(b) = other.t.get(j)
                    && b.key == a.key
                {
                    let c = reduce_rat(a.c.0 * b.c.1 + b.c.0 * a.c.1, a.c.1 * b.c.1);
                    if c.0 != 0 {
                        out.push(Term {
                            key: a.key.clone(),
                            c,
                            f: a.f.clone(),
                        });
                    }
                    i += 1;
                    j += 1;
                } else {
                    out.push(a.clone());
                    i += 1;
                }
            } else {
                out.push(other.t[j].clone());
                j += 1;
            }
        }
        Self { t: out }
    }

    fn mul(&self, other: &Self) -> Option<Self> {
        let mut out = Self::zero();
        for a in &self.t {
            for b in &other.t {
                let c = reduce_rat(a.c.0 * b.c.0, a.c.1 * b.c.1);
                if c.0 == 0 {
                    continue;
                }
                let mut f = a.f.clone();
                for (base, e) in &b.f {
                    if let Some(slot) = f.iter_mut().find(|(x, _)| x == base) {
                        slot.1 += e;
                    } else {
                        f.push((*base, *e));
                    }
                }
                f.retain(|(_, e)| *e != 0);
                f.sort_by_key(|(b, _)| b.to_string());
                let term = Term {
                    key: factor_key(&f),
                    c,
                    f,
                };
                out = out.add(&Self { t: vec![term] });
                if out.t.len() > COEFF_TERM_LIMIT {
                    return None;
                }
            }
        }
        Some(out)
    }

    fn pow(&self, n: i64) -> Option<Self> {
        if n.unsigned_abs() > 8 || self.t.len() > COEFF_TERM_LIMIT {
            return None;
        }
        let mut acc = Self::one();
        for _ in 0..n.unsigned_abs() {
            acc = acc.mul(self)?;
        }
        Some(if n < 0 { acc.inv()? } else { acc })
    }

    /// Reciprocal (fails on zero and on non-monomial values).
    fn inv(&self) -> Option<Self> {
        if self.t.is_empty() {
            return None;
        }
        let mut out = Self::zero();
        for m in &self.t {
            let f: Vec<(Atom<'a>, i64)> = m.f.iter().map(|(b, e)| (*b, -e)).collect();
            let mut key_f = f.clone();
            key_f.sort_by_key(|(b, _)| b.to_string());
            let term = Term {
                key: factor_key(&key_f),
                c: reduce_rat(m.c.1, m.c.0),
                f,
            };
            out = out.add(&Self { t: vec![term] });
        }
        Some(out)
    }
}

/// Cap on the number of monomials; past it the canonical form degrades to an
/// opaque base, which is still a sound (if weaker) representation.
const COEFF_TERM_LIMIT: usize = 512;

/// Canonicalise an atom as a coefficient. Never fails: anything that would
/// blow the budget is kept as an opaque base.
fn coeff_of<'a>(_ctx: &'a AtomArena<'a>, a: Atom<'a>, budget: &mut usize) -> Coeff<'a> {
    if *budget == 0 {
        return Coeff::base(a);
    }
    *budget -= 1;
    match a.node() {
        AtomNode::Num(n) => Coeff::num(*n),
        AtomNode::Var(_) | AtomNode::Fun(_, _) => Coeff::base(a),
        AtomNode::Add(args) => {
            let mut acc = Coeff::zero();
            for x in args.iter() {
                let c = coeff_of(_ctx, *x, budget);
                acc = acc.add(&c);
            }
            acc
        }
        AtomNode::Mul(args) => {
            let mut acc = Coeff::one();
            for x in args.iter() {
                let c = coeff_of(_ctx, *x, budget);
                match acc.mul(&c) {
                    Some(m) => acc = m,
                    None => return Coeff::base(a),
                }
            }
            acc
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Num(n) = e.node()
                && *n != 0
                && n.unsigned_abs() <= 8
            {
                let inner = coeff_of(_ctx, *b, budget);
                // A negative power of a *sum* cannot be expanded into
                // monomials; it falls through to the formal-reciprocal branch,
                // which is exact and still merges with other powers of the
                // same base.
                if *n > 0 || inner.t.len() == 1 {
                    if let Some(p) = inner.pow(*n) {
                        return p;
                    }
                } else {
                    let f = vec![(*b, *n)];
                    return Coeff {
                        t: vec![Term {
                            key: factor_key(&f),
                            c: (1, 1),
                            f,
                        }],
                    };
                }
            }
            Coeff::base(a)
        }
    }
}

/// Reciprocal of a coefficient: exact for a single monomial, and a formal
/// reciprocal (keyed by the coefficient's own atom) otherwise. Never wrong,
/// merely not always simplified.
fn coeff_recip<'a>(ctx: &'a AtomArena<'a>, c: &Coeff<'a>) -> Option<Coeff<'a>> {
    if c.is_zero() {
        return None;
    }
    if c.t.len() == 1 {
        return c.inv();
    }
    let atom = coeff_to_atom(ctx, c);
    let f = vec![(atom, -1i64)];
    Some(Coeff {
        t: vec![Term {
            key: factor_key(&f),
            c: (1, 1),
            f,
        }],
    })
}

fn coeff_is_zero(c: &Coeff<'_>) -> bool {
    c.is_zero()
}

/// Rebuild an atom from a canonical coefficient.
fn coeff_to_atom<'a>(ctx: &'a AtomArena<'a>, c: &Coeff<'a>) -> Atom<'a> {
    let mut terms: Vec<Atom<'a>> = Vec::with_capacity(c.t.len());
    for m in &c.t {
        let mut factors: Vec<Atom<'a>> = Vec::with_capacity(m.f.len() + 2);
        if m.c.0 != 1 || m.f.is_empty() {
            factors.push(ctx.num(i64::try_from(m.c.0).unwrap_or(i64::MAX)));
        }
        if m.c.1 != 1 {
            factors.push(ctx.pow(
                ctx.num(i64::try_from(m.c.1).unwrap_or(i64::MAX)),
                ctx.num(-1),
            ));
        }
        for (base, e) in &m.f {
            factors.push(if *e == 1 {
                *base
            } else {
                ctx.pow(*base, ctx.num(*e))
            });
        }
        terms.push(match factors.len() {
            0 => ctx.num(1),
            1 => factors[0],
            _ => ctx.mul(&factors),
        });
    }
    match terms.len() {
        0 => ctx.num(0),
        1 => terms[0],
        _ => ctx.add(&terms),
    }
}

// =========================================================================
// Polynomials with symbolic coefficients
// =========================================================================

/// Dense univariate polynomial: `p[i]` is the coefficient of `u^i`.
type Poly<'a> = Vec<Coeff<'a>>;

fn poly_trim<'a>(mut p: Poly<'a>) -> Poly<'a> {
    while let Some(last) = p.last() {
        if coeff_is_zero(last) {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn poly_deg(p: &[Coeff<'_>]) -> Option<usize> {
    p.len().checked_sub(1)
}

fn poly_add<'a>(a: &[Coeff<'a>], b: &[Coeff<'a>]) -> Poly<'a> {
    let mut out = a.to_vec();
    if out.len() < b.len() {
        out.resize(b.len(), Coeff::zero());
    }
    for (i, c) in b.iter().enumerate() {
        out[i] = out[i].add(c);
    }
    poly_trim(out)
}

fn poly_mul<'a>(a: &[Coeff<'a>], b: &[Coeff<'a>]) -> Option<Poly<'a>> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let mut out: Poly<'a> = vec![Coeff::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if coeff_is_zero(x) {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            let term = x.mul(y)?;
            out[i + j] = out[i + j].add(&term);
        }
    }
    Some(poly_trim(out))
}

/// Exact long division with symbolic coefficients.
///
/// The leading coefficient of the running remainder is removed by
/// construction (the quotient term is `lc(r)/lc(d)`, computed in the canonical
/// coefficient algebra), so the remainder returned is empty exactly when the
/// division is exact.
fn poly_divmod<'a>(num: &[Coeff<'a>], den: &[Coeff<'a>]) -> Option<(Poly<'a>, Poly<'a>)> {
    let den = poly_trim(den.to_vec());
    let dd = poly_deg(&den)?;
    if coeff_is_zero(&den[dd]) {
        return None;
    }
    let mut r = poly_trim(num.to_vec());
    let mut q: Poly<'a> = Vec::new();
    let mut guard = 0usize;
    while r.len() > dd {
        guard += 1;
        if guard > 64 {
            return None;
        }
        let shift = r.len() - 1 - dd;
        let factor = r[r.len() - 1].mul(&den[dd].inv()?)?;
        poly_add_term_dummy(&mut q, shift, &factor);
        let mut nr: Poly<'a> = r[..r.len() - 1].to_vec();
        for (i, d) in den.iter().enumerate().take(dd) {
            let sub = factor.mul(d)?;
            poly_add_term_dummy(&mut nr, shift + i, &sub.neg());
        }
        r = poly_trim(nr);
    }
    Some((poly_trim(q), r))
}

/// `poly_add_term` without the arena (the coefficient algebra is self
/// contained); resizing uses canonical zeros.
fn poly_add_term_dummy<'a>(p: &mut Poly<'a>, deg: usize, c: &Coeff<'a>) {
    if coeff_is_zero(c) {
        return;
    }
    if p.len() <= deg {
        p.resize(deg + 1, Coeff::zero());
    }
    p[deg] = p[deg].add(c);
}

fn poly_to_atom<'a>(ctx: &'a AtomArena<'a>, p: &[Coeff<'a>], u: Symbol) -> Atom<'a> {
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for (i, c) in p.iter().enumerate() {
        if coeff_is_zero(c) {
            continue;
        }
        let c_atom = coeff_to_atom(ctx, c);
        let t = match i {
            0 => c_atom,
            1 => ctx.mul(&[c_atom, ctx.var(u.as_str())]),
            _ => ctx.mul(&[c_atom, ctx.pow(ctx.var(u.as_str()), ctx.num(i as i64))]),
        };
        terms.push(t);
    }
    match terms.len() {
        0 => ctx.num(0),
        1 => terms[0],
        _ => ctx.add(&terms),
    }
}

// =========================================================================
// Rational functions with symbolic coefficients
// =========================================================================

#[derive(Clone)]
struct Rat<'a> {
    num: Poly<'a>,
    den: Poly<'a>,
}

impl<'a> Rat<'a> {
    fn one<'b>(_ctx: &'b AtomArena<'b>) -> Rat<'b> {
        Rat {
            num: vec![Coeff::one()],
            den: vec![Coeff::one()],
        }
    }

    fn constant(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Self {
        let mut budget = 4096usize;
        Rat {
            num: poly_trim(vec![coeff_of(ctx, a, &mut budget)]),
            den: vec![Coeff::one()],
        }
    }

    fn prod(&self, other: &Self) -> Option<Self> {
        Some(Rat {
            num: poly_mul(&self.num, &other.num)?,
            den: poly_mul(&self.den, &other.den)?,
        })
    }

    fn sum(_ctx: &'a AtomArena<'a>, a: &Self, b: &Self) -> Option<Self> {
        let n1 = poly_mul(&a.num, &b.den)?;
        let n2 = poly_mul(&b.num, &a.den)?;
        Some(Rat {
            num: poly_add(&n1, &n2),
            den: poly_mul(&a.den, &b.den)?,
        })
    }

    fn powi(&self, n: i64) -> Option<Self> {
        if n.unsigned_abs() > 24 {
            return None;
        }
        let one = vec![Coeff::one()];
        let mut num = one.clone();
        let mut den = one;
        for _ in 0..n.unsigned_abs() {
            num = poly_mul(&num, &self.num)?;
            den = poly_mul(&den, &self.den)?;
        }
        Some(if n >= 0 {
            Rat { num, den }
        } else {
            Rat { num: den, den: num }
        })
    }
}

/// Divide `num` by `den` when the division is exact, so that rationals carried
/// through the decomposition never accumulate spurious common factors (which
/// would inflate `deg Q`).
fn rat_reduce<'a>(num: &[Coeff<'a>], den: &[Coeff<'a>]) -> Option<(Poly<'a>, Poly<'a>)> {
    poly_deg(den)?;
    let (q, r) = poly_divmod(num, den)?;
    if r.is_empty() {
        Some((q, vec![Coeff::one()]))
    } else {
        Some((num.to_vec(), den.to_vec()))
    }
}

// =========================================================================
// Fractional-power decomposition
// =========================================================================

/// `expr = rat · √rad` when `has_rad`, and `expr = rat` otherwise.
struct Dec<'a> {
    rat: Rat<'a>,
    rad: Rat<'a>,
    has_rad: bool,
}

fn dec_const<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Dec<'a> {
    Dec {
        rat: Rat::constant(ctx, a),
        rad: Rat::one(ctx),
        has_rad: false,
    }
}

/// Apply an exponent `p/q` to an already decomposed base.
///
/// Only `q ∈ {1, 2}` is handled: integer powers, and half-powers of a rational
/// base. Everything else (cube roots, radicals of radicals) is declined.
fn apply_pow<'a>(_ctx: &'a AtomArena<'a>, inner: Dec<'a>, p: i64, q: i64) -> Option<Dec<'a>> {
    if q == 1 {
        return pow_int(inner, p);
    }
    if q != 2 || p.unsigned_abs() > 63 {
        return None;
    }
    if p % 2 == 0 {
        return pow_int(inner, p / 2);
    }
    if inner.has_rad {
        // A radical of a radical: `(r√s)^{p/2}` carries quarter powers.
        return None;
    }
    let m = (p - 1) / 2;
    let rat = inner.rat.powi(m)?;
    Some(Dec {
        rat,
        rad: inner.rat,
        has_rad: true,
    })
}

/// `base^n` for an integer `n`, where `base = rat·√rad`.
fn pow_int<'a>(d: Dec<'a>, n: i64) -> Option<Dec<'a>> {
    if n.unsigned_abs() > 24 {
        return None;
    }
    let k = n.div_euclid(2);
    let eps = n.rem_euclid(2);
    let rat = d.rat.powi(n)?.prod(&d.rad.powi(k)?)?;
    let has_rad = d.has_rad && eps == 1;
    Some(Dec {
        rat,
        rad: d.rad,
        has_rad,
    })
}

fn decompose<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    u: Symbol,
    budget: &mut usize,
) -> Option<Dec<'a>> {
    if *budget == 0 {
        return None;
    }
    *budget -= 1;
    match expr.node() {
        AtomNode::Num(_) => Some(dec_const(ctx, expr)),
        AtomNode::Var(v) => {
            if *v == u {
                // The polynomial variable itself: coefficient 1 at degree 1.
                Some(Dec {
                    rat: Rat {
                        num: vec![Coeff::zero(), Coeff::one()],
                        den: vec![Coeff::one()],
                    },
                    rad: Rat::one(ctx),
                    has_rad: false,
                })
            } else {
                Some(dec_const(ctx, expr))
            }
        }
        AtomNode::Add(args) => {
            let mut acc: Option<Rat<'a>> = None;
            for a in args.iter() {
                let d = decompose(ctx, *a, u, budget)?;
                if d.has_rad {
                    return None;
                }
                acc = Some(match acc {
                    None => d.rat,
                    Some(r) => Rat::sum(ctx, &r, &d.rat)?,
                });
            }
            let rat = acc?;
            Some(Dec {
                rat,
                rad: Rat::one(ctx),
                has_rad: false,
            })
        }
        AtomNode::Mul(args) => {
            let mut acc = Dec {
                rat: Rat::one(ctx),
                rad: Rat::one(ctx),
                has_rad: false,
            };
            for a in args.iter() {
                let d = decompose(ctx, *a, u, budget)?;
                let rad = acc.rad.prod(&d.rad)?;
                let rat = acc.rat.prod(&d.rat)?;
                acc = Dec {
                    rat,
                    rad,
                    has_rad: acc.has_rad || d.has_rad,
                };
            }
            Some(acc)
        }
        AtomNode::Pow(b, e) => {
            let inner = decompose(ctx, *b, u, budget)?;
            match exp_fraction(*e) {
                Some((p, q)) => apply_pow(ctx, inner, p, q),
                None => {
                    if is_constant(expr, u) {
                        Some(dec_const(ctx, expr))
                    } else {
                        None
                    }
                }
            }
        }
        AtomNode::Fun(name, args) => {
            if name.as_str() == "sqrt" && args.len() == 1 {
                let inner = decompose(ctx, args[0], u, budget)?;
                return apply_pow(ctx, inner, 1, 2);
            }
            if is_constant(expr, u) {
                Some(dec_const(ctx, expr))
            } else {
                None
            }
        }
    }
}

// =========================================================================
// Hermite reduction of `R/Q^{3/2}` for even `Q`
// =========================================================================

/// Reduce `R(u)/Q(u)^{3/2}`, `deg R < 4`, `Q = q₄u⁴ + q₂u² + q₀`, to
/// `A(u)/√Q + B(u)/√Q`.
///
/// `A = a₃u³ + a₂u² + a₁u + a₀` and `B = b₂u² + b₀` solve
/// `R = A'Q − A Q'/2 + B Q`; matching coefficients of `u⁰…u⁴` gives the
/// triangular system (with `Δ = q₂² − 4q₄q₀`):
///
/// ```text
///   a₃ = q₄(q₂n₀ − 2q₀n₂)/(q₀Δ)      b₂ = q₄(2q₀n₂ − q₂n₀)/(q₀Δ)
///   a₂ = (q₂n₃ − 2q₄n₁)/Δ            b₀ = (q₂n₂ − 2q₄n₀)/Δ
///   a₁ = (n₀(q₂² − 2q₄q₀) − q₂q₀n₂)/(q₀Δ)
///   a₀ = (2q₀n₃ − q₂n₁)/Δ
/// ```
///
/// The odd coefficients of `A` and `b₁ = 0` are what let an odd part of `R` be
/// reduced too. `Δ = (q₀ − q₄)²` for the Legendre normal form, which is used
/// directly to keep the emitted atoms small.
fn hermite_even<'a>(
    ctx: &'a AtomArena<'a>,
    q: &[Coeff<'a>],
    r: &[Coeff<'a>],
) -> Option<(Poly<'a>, Poly<'a>)> {
    let q0 = q.first().cloned().unwrap_or_else(Coeff::zero);
    let q2 = q.get(2).cloned().unwrap_or_else(Coeff::zero);
    let q4 = q.get(4).cloned().unwrap_or_else(Coeff::zero);
    let n0 = r.first().cloned().unwrap_or_else(Coeff::zero);
    let n1 = r.get(1).cloned().unwrap_or_else(Coeff::zero);
    let n2 = r.get(2).cloned().unwrap_or_else(Coeff::zero);
    let n3 = r.get(3).cloned().unwrap_or_else(Coeff::zero);
    if coeff_is_zero(&n0) && coeff_is_zero(&n1) && coeff_is_zero(&n2) && coeff_is_zero(&n3) {
        return Some((Vec::new(), Vec::new()));
    }
    let two = Coeff::num(2);
    let q0_q4 = q0.add(&q4.neg());
    let delta = q0_q4.mul(&q0_q4)?;
    let inv_delta = coeff_recip(ctx, &delta)?;
    let q0_delta = q0.mul(&delta)?;
    let inv_q0_delta = coeff_recip(ctx, &q0_delta)?;
    let a3 = q4
        .mul(&q2.mul(&n0)?.add(&two.mul(&q0)?.mul(&n2)?.neg()))?
        .mul(&inv_q0_delta)?;
    let a2 = q2
        .mul(&n3)?
        .add(&two.mul(&q4)?.mul(&n1)?.neg())
        .mul(&inv_delta)?;
    let inner1 = q2.mul(&q2)?.add(&two.mul(&q4)?.mul(&q0)?.neg());
    let a1 = n0
        .mul(&inner1)?
        .add(&q2.mul(&q0)?.mul(&n2)?.neg())
        .mul(&inv_q0_delta)?;
    let a0 = two
        .mul(&q0)?
        .mul(&n3)?
        .add(&q2.mul(&n1)?.neg())
        .mul(&inv_delta)?;
    let b2 = q4
        .mul(&two.mul(&q0)?.mul(&n2)?.add(&q2.mul(&n0)?.neg()))?
        .mul(&inv_q0_delta)?;
    let b0 = q2
        .mul(&n2)?
        .add(&two.mul(&q4)?.mul(&n0)?.neg())
        .mul(&inv_delta)?;
    let a = poly_trim(vec![a0, a1, a2, a3]);
    let b = poly_trim(vec![b0, Coeff::zero(), b2]);
    Some((a, b))
}

// =========================================================================
// Entry point
// =========================================================================

/// Integrate `expr` by elliptic reduction.
///
/// Returns `None` outside the squarefree quartic, Legendre-normal-form family
/// documented at the module level. Every accepted shape is verified by the
/// module tests against a Carlson-form numeric oracle.
pub(crate) fn integrate_elliptic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if super::node_count(expr) > 400 {
        return None;
    }
    let mut budget = 4_000usize;
    let dec = decompose(ctx, expr, var, &mut budget)?;
    if !dec.has_rad {
        return None;
    }
    // Reduce the rationals first: unreduced common factors would inflate the
    // degree of Q and defeat the normal-form test below.
    let (rnum, rden) = rat_reduce(&dec.rat.num, &dec.rat.den)?;
    let (qnum, qden) = rat_reduce(&dec.rad.num, &dec.rad.den)?;
    // The integrand is `W/√Q` with `Q = num(rad)·den(rad)` and
    // `W = rat·num(rad)`; see the module docs.
    let q = poly_mul(&qnum, &qden)?;
    if poly_deg(&q) != Some(4) {
        return None;
    }
    let q0 = q[0].clone();
    let q1 = q[1].clone();
    let q2 = q[2].clone();
    let q3 = q[3].clone();
    let q4 = q[4].clone();
    if !coeff_is_zero(&q1) || !coeff_is_zero(&q3) {
        return None;
    }
    if coeff_is_zero(&q0) || coeff_is_zero(&q4) {
        return None;
    }
    // Legendre normal form: q₀ + q₂ + q₄ ≡ 0.
    if !coeff_is_zero(&q0.add(&q2).add(&q4)) {
        return None;
    }
    let c = q0.clone();
    let m = q4.mul(&coeff_recip(ctx, &q0)?)?;
    let m_inv = q0.mul(&coeff_recip(ctx, &q4)?)?;
    // N = W·Q = rat·num(rad)·Q / den(rat) must be a polynomial.
    let wq = poly_mul(&poly_mul(&rnum, &qnum)?, &q)?;
    let (nprime, rem) = poly_divmod(&wq, &rden)?;
    if !rem.is_empty() {
        return None;
    }
    let (s_quot, r) = poly_divmod(&nprime, &q)?;
    if s_quot.len() > 4 {
        return None;
    }
    let (a_h, b_h) = hermite_even(ctx, &q, &r)?;
    let t = poly_add(&s_quot, &b_h);
    if t.len() > 4 {
        return None;
    }
    let t0 = t.first().cloned().unwrap_or_else(Coeff::zero);
    let t1 = t.get(1).cloned().unwrap_or_else(Coeff::zero);
    let t2 = t.get(2).cloned().unwrap_or_else(Coeff::zero);
    let t3 = t.get(3).cloned().unwrap_or_else(Coeff::zero);
    if !coeff_is_zero(&t1) || !coeff_is_zero(&t3) {
        return None;
    }
    let t2_over_m = t2.mul(&m_inv)?;
    let c_f_inner = t0.add(&t2_over_m);
    let c_e_inner = t2_over_m.neg();
    // `1/√c` is a fractional power, so it stays an atom factor. `c = 1` is
    // folded away so the Legendre-normalised cases read cleanly.
    let c_is_one = c.t.len() == 1 && c.t[0].f.is_empty() && c.t[0].c == (1, 1);
    let inv_sqrt_c = if c_is_one {
        ctx.num(1)
    } else {
        ctx.pow(coeff_to_atom(ctx, &c), half_exp(ctx, -1))
    };

    let u_atom = ctx.var(var.as_str());
    let phi = ctx.fun("asin", &[u_atom]);
    let m_atom = coeff_to_atom(ctx, &m);
    let mut terms: Vec<Atom<'a>> = Vec::new();
    if !a_h.is_empty() {
        let ah = poly_to_atom(ctx, &a_h, var);
        let qa = poly_to_atom(ctx, &q, var);
        terms.push(ctx.mul(&[ah, ctx.pow(qa, half_exp(ctx, -1))]));
    }
    if !coeff_is_zero(&c_f_inner) {
        terms.push(ctx.mul(&[
            inv_sqrt_c,
            coeff_to_atom(ctx, &c_f_inner),
            ctx.fun("EllipticF", &[phi, m_atom]),
        ]));
    }
    if !coeff_is_zero(&c_e_inner) {
        terms.push(ctx.mul(&[
            inv_sqrt_c,
            coeff_to_atom(ctx, &c_e_inner),
            ctx.fun("EllipticE", &[phi, m_atom]),
        ]));
    }
    if terms.is_empty() {
        return None;
    }
    let out = if terms.len() == 1 {
        terms[0]
    } else {
        ctx.add(&terms)
    };
    Some(ocas_atom::normalize::normalize(ctx, out))
}

// =========================================================================
// Numeric oracle (test-only): Carlson symmetric forms and a differentiator
// =========================================================================

#[cfg(test)]
pub(crate) mod testnum {
    use std::f64::consts::FRAC_PI_2;

    use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

    const ERRTOL: f64 = 0.0025;
    const MAXIT: usize = 60;
    /// Coarse sanity bound on the parameters of the elliptic heads. The real
    /// domain guard is the integrand itself: `1 − m sin²θ > 0` (and
    /// `1 − n sin²θ > 0`) is checked explicitly below, so values `m > 1` that
    /// keep the amplitude on the real branch — e.g. `F(asin x, 4)` for
    /// `|x| < 1/2` — are evaluated rather than declined.
    const PARAM_LIMIT: f64 = 1.0e4;

    /// Carlson's `R_F(x, y, z)` (duplication algorithm, Numerical Recipes
    /// §6.11). Returns `None` outside the real domain or on non-convergence.
    pub(crate) fn rf(mut x: f64, mut y: f64, mut z: f64) -> Option<f64> {
        if x < 0.0 || y < 0.0 || z < 0.0 {
            return None;
        }
        let mut ave;
        let mut delx;
        let mut dely;
        let mut delz;
        for _ in 0..MAXIT {
            let sx = x.sqrt();
            let sy = y.sqrt();
            let sz = z.sqrt();
            let alamb = sx * (sy + sz) + sy * sz;
            x = 0.25 * (x + alamb);
            y = 0.25 * (y + alamb);
            z = 0.25 * (z + alamb);
            ave = (x + y + z) / 3.0;
            delx = (ave - x) / ave;
            dely = (ave - y) / ave;
            delz = (ave - z) / ave;
            if delx.abs().max(dely.abs()).max(delz.abs()) < ERRTOL {
                let e2 = delx * dely - delz * delz;
                let e3 = delx * dely * delz;
                return Some(
                    (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0)
                        / ave.sqrt(),
                );
            }
        }
        None
    }

    /// `R_C(x, y) = R_F(x, y, y)`; only the principal real branch with `y > 0`
    /// is needed by the elliptic heads here.
    pub(crate) fn rc(x: f64, y: f64) -> Option<f64> {
        rf(x, y, y)
    }

    /// Carlson's `R_D(x, y, z)`.
    pub(crate) fn rd(mut x: f64, mut y: f64, mut z: f64) -> Option<f64> {
        const C1: f64 = 3.0 / 14.0;
        const C2: f64 = 1.0 / 6.0;
        const C3: f64 = 3.0 / 22.0;
        const C4: f64 = 3.0 / 26.0;
        const C5: f64 = 0.25 * 3.0 / 22.0;
        const C6: f64 = 1.5 * 3.0 / 26.0;
        if x < 0.0 || y < 0.0 || z <= 0.0 {
            return None;
        }
        let mut sum = 0.0;
        let mut fac = 1.0;
        for _ in 0..MAXIT {
            let sx = x.sqrt();
            let sy = y.sqrt();
            let sz = z.sqrt();
            let alamb = sx * (sy + sz) + sy * sz;
            sum += fac / (sz * (z + alamb));
            fac *= 0.25;
            x = 0.25 * (x + alamb);
            y = 0.25 * (y + alamb);
            z = 0.25 * (z + alamb);
            let ave = 0.2 * (x + y + 3.0 * z);
            let delx = (ave - x) / ave;
            let dely = (ave - y) / ave;
            let delz = (ave - z) / ave;
            if delx.abs().max(dely.abs()).max(delz.abs()) < ERRTOL {
                let ea = delx * dely;
                let eb = delz * delz;
                let ec = ea - eb;
                let ed = ea - 6.0 * eb;
                let ee = ed + ec + ec;
                return Some(
                    3.0 * sum
                        + fac
                            * (1.0
                                + ed * (-C1 + C5 * ed - C6 * delz * ee)
                                + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea)))
                            / (ave * ave.sqrt()),
                );
            }
        }
        None
    }

    /// Carlson's `R_J(x, y, z, p)` for positive arguments.
    pub(crate) fn rj(mut x: f64, mut y: f64, mut z: f64, mut p: f64) -> Option<f64> {
        const C1: f64 = 3.0 / 14.0;
        const C2: f64 = 1.0 / 3.0;
        const C3: f64 = 3.0 / 22.0;
        const C4: f64 = 3.0 / 26.0;
        const C5: f64 = 0.75 * 3.0 / 22.0;
        const C6: f64 = 1.5 * 3.0 / 26.0;
        const C7: f64 = 0.5 * 1.0 / 3.0;
        const C8: f64 = 3.0 / 22.0 + 3.0 / 22.0;
        if x < 0.0 || y < 0.0 || z < 0.0 || p <= 0.0 {
            return None;
        }
        let mut sum = 0.0;
        let mut fac = 1.0;
        for _ in 0..MAXIT {
            let sx = x.sqrt();
            let sy = y.sqrt();
            let sz = z.sqrt();
            let alamb = sx * (sy + sz) + sy * sz;
            let alpha = (p * (sx + sy + sz) + sx * sy * sz).powi(2);
            let beta = p * (p + alamb).powi(2);
            sum += fac * rc(alpha, beta)?;
            fac *= 0.25;
            x = 0.25 * (x + alamb);
            y = 0.25 * (y + alamb);
            z = 0.25 * (z + alamb);
            p = 0.25 * (p + alamb);
            let ave = 0.2 * (x + y + z + 2.0 * p);
            let delx = (ave - x) / ave;
            let dely = (ave - y) / ave;
            let delz = (ave - z) / ave;
            let delp = (ave - p) / ave;
            if delx.abs().max(dely.abs()).max(delz.abs()).max(delp.abs()) < ERRTOL {
                let ea = delx * (dely + delz) + dely * delz;
                let eb = delx * dely * delz;
                let ec = delp * delp;
                let ed = ea - 3.0 * ec;
                let ee = eb + 2.0 * delp * (ea - ec);
                return Some(
                    3.0 * sum
                        + fac
                            * (1.0
                                + ed * (-C1 + C5 * ed - C6 * ee)
                                + eb * (C7 + delp * (-C8 + delp * C4))
                                + delp * ea * (C2 - delp * C3)
                                - C2 * delp * ec)
                            / (ave * ave.sqrt()),
                );
            }
        }
        None
    }

    /// `EllipticF(φ, m) = sin φ · R_F(cos²φ, 1 − m sin²φ, 1)`.
    pub(crate) fn ellip_f(phi: f64, m: f64) -> Option<f64> {
        if phi.abs() > FRAC_PI_2 + 1e-12 || m.abs() > PARAM_LIMIT {
            return None;
        }
        let s = phi.sin();
        let q = 1.0 - m * s * s;
        if q <= 0.0 {
            return None;
        }
        Some(s * rf(phi.cos() * phi.cos(), q, 1.0)?)
    }

    /// `EllipticE(φ, m) = sin φ R_F − (m/3) sin³φ R_D`.
    pub(crate) fn ellip_e(phi: f64, m: f64) -> Option<f64> {
        if phi.abs() > FRAC_PI_2 + 1e-12 || m.abs() > PARAM_LIMIT {
            return None;
        }
        let s = phi.sin();
        let q = 1.0 - m * s * s;
        if q <= 0.0 {
            return None;
        }
        let c2 = phi.cos() * phi.cos();
        Some(s * rf(c2, q, 1.0)? - (m / 3.0) * s.powi(3) * rd(c2, q, 1.0)?)
    }

    /// `EllipticPi(n, φ, m) = sin φ R_F + (n/3) sin³φ R_J`.
    pub(crate) fn ellip_pi(n: f64, phi: f64, m: f64) -> Option<f64> {
        if phi.abs() > FRAC_PI_2 + 1e-12 || m.abs() > PARAM_LIMIT || n.abs() > PARAM_LIMIT {
            return None;
        }
        let s = phi.sin();
        let q = 1.0 - m * s * s;
        let pp = 1.0 - n * s * s;
        if q <= 0.0 || pp <= 0.0 {
            return None;
        }
        let c2 = phi.cos() * phi.cos();
        Some(s * rf(c2, q, 1.0)? + (n / 3.0) * s.powi(3) * rj(c2, q, 1.0, pp)?)
    }

    /// Evaluate `expr` numerically; `None` marks a domain/unsupported failure.
    pub(crate) fn eval_f64(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
        match expr.node() {
            AtomNode::Num(n) => Some(*n as f64),
            AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, val)| *val),
            AtomNode::Add(args) => {
                let mut acc = 0.0;
                for a in args.iter() {
                    acc += eval_f64(*a, env)?;
                }
                Some(acc)
            }
            AtomNode::Mul(args) => {
                let mut acc = 1.0;
                for a in args.iter() {
                    acc *= eval_f64(*a, env)?;
                }
                Some(acc)
            }
            AtomNode::Pow(b, e) => {
                let base = eval_f64(*b, env)?;
                let ex = eval_f64(*e, env)?;
                if base < 0.0 && ex.fract() != 0.0 {
                    return None;
                }
                Some(base.powf(ex))
            }
            AtomNode::Fun(name, args) => {
                let out = match name.as_str() {
                    "EllipticF" | "EllipticE" => {
                        let phi = eval_f64(*args.first()?, env)?;
                        let m = eval_f64(*args.get(1)?, env)?;
                        if name.as_str() == "EllipticF" {
                            ellip_f(phi, m)?
                        } else {
                            ellip_e(phi, m)?
                        }
                    }
                    "EllipticPi" => {
                        let n = eval_f64(*args.first()?, env)?;
                        let phi = eval_f64(*args.get(1)?, env)?;
                        let m = eval_f64(*args.get(2)?, env)?;
                        ellip_pi(n, phi, m)?
                    }
                    _ => {
                        let v = eval_f64(*args.first()?, env)?;
                        match name.as_str() {
                            "sin" => v.sin(),
                            "cos" => v.cos(),
                            "tan" => v.tan(),
                            "cot" => 1.0 / v.tan(),
                            "sec" => 1.0 / v.cos(),
                            "csc" => 1.0 / v.sin(),
                            "exp" => v.exp(),
                            "log" => v.abs().ln(),
                            "sqrt" => v.sqrt(),
                            "asin" => {
                                if !(-1.0..=1.0).contains(&v) {
                                    return None;
                                }
                                v.asin()
                            }
                            "acos" => {
                                if !(-1.0..=1.0).contains(&v) {
                                    return None;
                                }
                                v.acos()
                            }
                            "atan" => v.atan(),
                            "sinh" => v.sinh(),
                            "cosh" => v.cosh(),
                            "tanh" => v.tanh(),
                            "asinh" => v.asinh(),
                            "acosh" => {
                                if v < 1.0 {
                                    return None;
                                }
                                v.acosh()
                            }
                            _ => return None,
                        }
                    }
                };
                if out.is_finite() { Some(out) } else { None }
            }
        }
    }

    /// Symbolic differentiation supporting the three elliptic heads.
    ///
    /// `m` (and `n`) are parameters of the amplitude, not of the integration
    /// variable, so only the amplitude derivative is needed for the emitted
    /// forms; their dependence on `var` is asserted away so a form that *did*
    /// depend on them could not be silently mis-verified.
    pub(crate) fn diff_local<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
        match expr.node() {
            AtomNode::Num(_) => ctx.num(0),
            AtomNode::Var(v) => {
                if *v == var {
                    ctx.num(1)
                } else {
                    ctx.num(0)
                }
            }
            AtomNode::Add(args) => {
                let d: Vec<Atom<'a>> = args.iter().map(|a| diff_local(ctx, *a, var)).collect();
                ctx.add(&d)
            }
            AtomNode::Mul(args) => {
                let mut terms = Vec::with_capacity(args.len());
                for i in 0..args.len() {
                    let mut factors = Vec::with_capacity(args.len());
                    for (j, a) in args.iter().enumerate() {
                        if i == j {
                            factors.push(diff_local(ctx, *a, var));
                        } else {
                            factors.push(*a);
                        }
                    }
                    terms.push(ctx.mul(&factors));
                }
                ctx.add(&terms)
            }
            AtomNode::Pow(b, e) => {
                let base = *b;
                let ex = *e;
                let db = diff_local(ctx, base, var);
                let de = diff_local(ctx, ex, var);
                let t1 = ctx.mul(&[ex, ctx.pow(base, ctx.add(&[ex, ctx.num(-1)])), db]);
                let t2 = ctx.mul(&[ctx.pow(base, ex), ctx.fun("log", &[base]), de]);
                ctx.add(&[t1, t2])
            }
            AtomNode::Fun(name, args) => {
                let a0 = *args.first().unwrap();
                let da0 = diff_local(ctx, a0, var);
                let d0 = match name.as_str() {
                    "sin" => ctx.fun("cos", &[a0]),
                    "cos" => ctx.mul(&[ctx.num(-1), ctx.fun("sin", &[a0])]),
                    "exp" => ctx.fun("exp", &[a0]),
                    "log" => ctx.pow(a0, ctx.num(-1)),
                    "sqrt" => ctx.pow(ctx.mul(&[ctx.num(2), ctx.fun("sqrt", &[a0])]), ctx.num(-1)),
                    "tan" => ctx.pow(ctx.fun("sec", &[a0]), ctx.num(2)),
                    "asin" => ctx.pow(
                        ctx.fun(
                            "sqrt",
                            &[ctx.add(&[
                                ctx.num(1),
                                ctx.mul(&[ctx.num(-1), ctx.pow(a0, ctx.num(2))]),
                            ])],
                        ),
                        ctx.num(-1),
                    ),
                    "atan" => ctx.pow(ctx.add(&[ctx.num(1), ctx.pow(a0, ctx.num(2))]), ctx.num(-1)),
                    "sinh" => ctx.fun("cosh", &[a0]),
                    "cosh" => ctx.fun("sinh", &[a0]),
                    "tanh" => ctx.pow(ctx.fun("sech", &[a0]), ctx.num(2)),
                    "asinh" => ctx.pow(
                        ctx.fun("sqrt", &[ctx.add(&[ctx.pow(a0, ctx.num(2)), ctx.num(1)])]),
                        ctx.num(-1),
                    ),
                    "acosh" => ctx.pow(
                        ctx.fun("sqrt", &[ctx.add(&[ctx.pow(a0, ctx.num(2)), ctx.num(-1)])]),
                        ctx.num(-1),
                    ),
                    "EllipticF" | "EllipticE" => {
                        let m = *args.get(1).unwrap();
                        let s = ctx.fun("sin", &[a0]);
                        let root = ctx.fun(
                            "sqrt",
                            &[ctx.add(&[
                                ctx.num(1),
                                ctx.mul(&[ctx.num(-1), m, ctx.pow(s, ctx.num(2))]),
                            ])],
                        );
                        let dphi = if name.as_str() == "EllipticF" {
                            ctx.pow(root, ctx.num(-1))
                        } else {
                            root
                        };
                        let dm = diff_local(ctx, m, var);
                        assert_eq!(
                            ocas_atom::normalize::normalize(ctx, dm).to_string(),
                            "0",
                            "modulus depends on {var:?}: {m}"
                        );
                        ctx.mul(&[dphi, da0])
                    }
                    "EllipticPi" => {
                        let n = *args.first().unwrap();
                        let phi = *args.get(1).unwrap();
                        let m = *args.get(2).unwrap();
                        assert_eq!(
                            ocas_atom::normalize::normalize(ctx, diff_local(ctx, m, var))
                                .to_string(),
                            "0",
                            "modulus depends on {var:?}: {m}"
                        );
                        assert_eq!(
                            ocas_atom::normalize::normalize(ctx, diff_local(ctx, n, var))
                                .to_string(),
                            "0",
                            "characteristic depends on {var:?}: {n}"
                        );
                        let s = ctx.fun("sin", &[phi]);
                        let root = ctx.fun(
                            "sqrt",
                            &[ctx.add(&[
                                ctx.num(1),
                                ctx.mul(&[ctx.num(-1), m, ctx.pow(s, ctx.num(2))]),
                            ])],
                        );
                        let den = ctx.add(&[
                            ctx.num(1),
                            ctx.mul(&[ctx.num(-1), n, ctx.pow(s, ctx.num(2))]),
                        ]);
                        let dphi = ctx.pow(ctx.mul(&[den, root]), ctx.num(-1));
                        ctx.mul(&[dphi, diff_local(ctx, phi, var)])
                    }
                    _ => ctx.fun(
                        "Derivative",
                        &[ctx.fun(name.as_str(), args), ctx.var(var.as_str())],
                    ),
                };
                if matches!(name.as_str(), "EllipticF" | "EllipticE" | "EllipticPi") {
                    d0
                } else {
                    ctx.mul(&[d0, da0])
                }
            }
        }
    }

    /// Reference values by adaptive Simpson quadrature of the defining
    /// integral — an independent cross-check of the Carlson forms.
    pub(crate) fn ellip_f_simpson(phi: f64, m: f64) -> Option<f64> {
        let f = |t: f64| 1.0 / (1.0 - m * t.sin() * t.sin()).sqrt();
        simpson(&f, 0.0, phi, 1e-13, 40)
    }

    pub(crate) fn ellip_e_simpson(phi: f64, m: f64) -> Option<f64> {
        let f = |t: f64| (1.0 - m * t.sin() * t.sin()).sqrt();
        simpson(&f, 0.0, phi, 1e-13, 40)
    }

    pub(crate) fn ellip_pi_simpson(n: f64, phi: f64, m: f64) -> Option<f64> {
        let f = |t: f64| {
            let s = t.sin();
            1.0 / ((1.0 - n * s * s) * (1.0 - m * s * s).sqrt())
        };
        simpson(&f, 0.0, phi, 1e-13, 40)
    }

    fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64, depth: u32) -> Option<f64> {
        fn s(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
            let m = 0.5 * (a + b);
            (b - a) / 6.0 * (f(a) + 4.0 * f(m) + f(b))
        }
        fn rec(f: &dyn Fn(f64) -> f64, a: f64, b: f64, whole: f64, tol: f64, depth: u32) -> f64 {
            let m = 0.5 * (a + b);
            let left = s(f, a, m);
            let right = s(f, m, b);
            let delta = left + right - whole;
            if depth == 0 || delta.abs() <= 15.0 * tol {
                return left + right + delta / 15.0;
            }
            rec(f, a, m, left, 0.5 * tol, depth - 1) + rec(f, m, b, right, 0.5 * tol, depth - 1)
        }
        if a == b {
            return Some(0.0);
        }
        let whole = s(f, a, b);
        let v = rec(f, a, b, whole, tol, depth);
        if v.is_finite() { Some(v) } else { None }
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::testnum::*;
    use super::*;
    use ocas_core::arena::Arena;

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, s).unwrap_or_else(|_| panic!("parse failed: {s}"))
    }

    /// Differentiate the emitted antiderivative and compare against the
    /// integrand at fixed sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_elliptic(ctx, integrand, var).expect("mechanism declined");
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        let d = diff_local(ctx, result, var);
        let mut checked = 0usize;
        let mut worst = 0.0f64;
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = match eval_f64(d, &env) {
                Some(v) => v,
                None => continue,
            };
            let rhs = eval_f64(integrand, &env).expect("eval integrand");
            let tol = 1e-9 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at {var:?}={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
            checked += 1;
            worst = worst.max((lhs - rhs).abs() / rhs.abs().max(1.0));
        }
        assert!(checked >= 2, "only {checked} usable samples for {result}");
        println!("  checked={checked} worst_rel={worst:e}  ->  {result}");
    }

    fn declines<'a>(ctx: &'a AtomArena<'a>, s: &str) {
        let e = parse(ctx, s);
        assert!(
            integrate_elliptic(ctx, e, Symbol::new("x")).is_none(),
            "expected decline for {s}"
        );
    }

    // ---------------- Carlson forms against quadrature ----------------

    #[test]
    fn carlson_matches_quadrature() {
        for &(phi, m) in &[
            (0.3, 0.5),
            (1.0, -0.5),
            (std::f64::consts::FRAC_PI_2, 0.5),
            (0.7, 0.9),
            (-0.9, 0.25),
            (std::f64::consts::FRAC_PI_2, -0.99),
        ] {
            let f = ellip_f(phi, m).unwrap();
            let fs = ellip_f_simpson(phi, m).unwrap();
            assert!((f - fs).abs() < 1e-10, "F: {f} vs {fs} at ({phi},{m})");
            let e = ellip_e(phi, m).unwrap();
            let es = ellip_e_simpson(phi, m).unwrap();
            assert!((e - es).abs() < 1e-10, "E: {e} vs {es} at ({phi},{m})");
        }
        for &(n, phi, m) in &[(0.5, 0.8, 0.4), (-0.5, 1.2, 0.6), (0.3, 1.5, -0.4)] {
            let p = ellip_pi(n, phi, m).unwrap();
            let ps = ellip_pi_simpson(n, phi, m).unwrap();
            assert!((p - ps).abs() < 1e-10, "Pi: {p} vs {ps}");
        }
        // K(1/2) and E(π/2, 1/2) reference values.
        use std::f64::consts::FRAC_PI_2;
        assert!((ellip_f(FRAC_PI_2, 0.5).unwrap() - 1.854_074_677_301_372).abs() < 1e-12);
        assert!((ellip_e(FRAC_PI_2, 0.5).unwrap() - 1.350_643_881_047_675_5).abs() < 1e-12);
        // Out-of-domain arguments decline: amplitude past π/2, and a
        // parameter that pushes `1 − m sin²φ` off the real branch.
        assert!(ellip_f(2.0, 0.5).is_none());
        assert!(ellip_f(0.5, 100.0).is_none());
        // `m > 1` is *inside* the domain while `1 − m sin²φ > 0`.
        assert!(ellip_f(0.3, 4.0).is_some());
    }

    // ---------------- coefficient algebra ----------------

    #[test]
    fn canonical_coefficients_cancel() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for src in [
            "(a+b)*(a+b) - a^2 - 2*a*b - b^2",
            "(a+b)*M - (M*a + M*b)",
            "b*(a+b)^-1*(a+b)^-2 - b*(a+b)^-3",
        ] {
            let mut b = 4096usize;
            let e = parse(&ctx, src);
            assert!(
                coeff_is_zero(&coeff_of(&ctx, e, &mut b)),
                "not canonical zero: {src}"
            );
        }
    }

    #[test]
    fn symbolic_polynomial_division_is_exact() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Q = ((a+b) - b z²)(1 - z²) = (a+b) - (a+b+b)z² + b z⁴, even and in
        // Legendre normal form (q₀+q₂+q₄ = a+b-(a+2b)+b = 0).
        let q_atom = parse(&ctx, "((a+b) - b*z^2)*(1 - z^2)");
        let mut b = 4096usize;
        // Coefficients of Q² by exact polynomial multiplication of Q with
        // itself: the division Q²/Q must have empty remainder and quotient Q.
        let q_poly: Poly<'_> = vec![
            coeff_of(&ctx, parse(&ctx, "a+b"), &mut b),
            Coeff::zero(),
            coeff_of(&ctx, parse(&ctx, "-(a+b) - b"), &mut b),
            Coeff::zero(),
            coeff_of(&ctx, parse(&ctx, "b"), &mut b),
        ];
        let _ = q_atom;
        let sq = poly_mul(&q_poly, &q_poly).expect("mul");
        let (quot, rem) = poly_divmod(&sq, &q_poly).expect("division");
        assert!(rem.is_empty(), "remainder not empty");
        assert!(quot == q_poly, "quotient mismatch");
        // q₀ + q₂ + q₄ is exactly zero in the canonical algebra.
        let sum = q_poly[0].add(&q_poly[2]).add(&q_poly[4]);
        assert!(coeff_is_zero(&sum));
    }

    // ---------------- correctness ----------------

    #[test]
    fn legendre_quartic_first_kind() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "1/sqrt((1-x^2)*(1-k^2*x^2))");
        let env = [(Symbol::new("k"), 0.5)];
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &env, &[-0.8, -0.4, 0.2, 0.6]);
    }

    #[test]
    fn one_over_sqrt_one_minus_x4() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "1/sqrt(1-x^4)");
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &[], &[-0.8, -0.3, 0.25, 0.7]);
    }

    #[test]
    fn hermite_third_power_radical() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // ∫dx/(1-x⁴)^{3/2} = x/(2√(1-x⁴)) + F(asin x, -1)/2
        let e = parse(&ctx, "1/(1-x^4)^(3/2)");
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &[], &[-0.7, -0.2, 0.35, 0.8]);
    }

    #[test]
    fn hermite_with_even_u2_numerator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // ∫x²dx/(1-x⁴)^{3/2} = x³/(2√(1-x⁴)) + (F - E)/2
        let e = parse(&ctx, "x^2/(1-x^4)^(3/2)");
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &[], &[-0.7, -0.25, 0.3, 0.85]);
    }

    #[test]
    fn legendre_quartic_with_u2_numerator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "x^2/sqrt((1-x^2)*(1-(4^-1)*x^2))");
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &[], &[-0.8, -0.3, 0.4, 0.9]);
    }

    #[test]
    fn legendre_quartic_with_constant_numerator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "(1+x^2)/sqrt((1-x^2)*(1-(2^-1)*x^2))");
        assert_antiderivative_num(&ctx, e, Symbol::new("x"), &[], &[-0.8, -0.4, 0.3, 0.9]);
    }

    #[test]
    fn symbolic_modulus_and_parameter_families() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let s = Symbol::new("x");
        for (src, env, samples) in [
            (
                "1/sqrt((1-x^2)*(1-m*x^2))",
                vec![(Symbol::new("m"), 0.3)],
                vec![-0.8, -0.2, 0.4, 0.85],
            ),
            (
                "1/((1-x^2)*(1-m*x^2))^(3/2)",
                vec![(Symbol::new("m"), 0.3)],
                vec![-0.7, -0.3, 0.35, 0.8],
            ),
            (
                "x^2/((1-x^2)*(1-m*x^2))^(3/2)",
                vec![(Symbol::new("m"), 0.4)],
                vec![-0.75, -0.2, 0.5, 0.9],
            ),
            (
                "(1+2*x^2)/sqrt((1-x^2)*(1-m*x^2))",
                vec![(Symbol::new("m"), 0.5)],
                vec![-0.8, -0.3, 0.45, 0.9],
            ),
            (
                "1/sqrt((1-x^2)*(1+2*x^2))",
                vec![],
                vec![-0.8, -0.4, 0.3, 0.85],
            ),
            (
                "1/sqrt((1-x^2)*((2^-1)+(4^-1)*x^2))",
                vec![],
                vec![-0.8, -0.4, 0.3, 0.85],
            ),
            (
                "1/sqrt((1-x^2)*(1-(4*5^-1)*x^2))",
                vec![],
                vec![-0.8, -0.3, 0.4, 0.9],
            ),
            (
                "3/sqrt((1 - 4*x^2)*(1-x^2))",
                vec![],
                vec![-0.4, -0.2, 0.15, 0.45],
            ),
        ] {
            let e = parse(&ctx, src);
            assert_antiderivative_num(&ctx, e, s, &env, &samples);
        }
    }

    #[test]
    fn emitted_heads_keep_argument_order() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "1/sqrt((1-x^2)*(1-k^2*x^2))");
        let r = integrate_elliptic(&ctx, e, Symbol::new("x")).expect("solved");
        let s = r.to_string();
        assert!(s.contains("EllipticF"), "got {s}");
        let norm = ocas_atom::normalize::normalize(&ctx, r);
        assert_eq!(norm.to_string(), s, "normalize changed the result");
    }

    /// The rewritten form produced by `halfpower` for the quadratic
    /// `cos²` branch: `∫ (A − b z²)^{p/2} (1 − z²)^{-1/2} dz`.
    #[test]
    fn halfpower_rewritten_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [(Symbol::new("a"), 2.0), (Symbol::new("b"), 1.0)];
        let s = Symbol::new("z");
        for (src, samples) in [
            (
                "(a+b-b*z^2)^(-1/2)*(1-z^2)^(-1/2)",
                vec![-0.8, -0.4, 0.3, 0.9],
            ),
            (
                "(a+b-b*z^2)^(-3/2)*(1-z^2)^(-1/2)",
                vec![-0.75, -0.3, 0.4, 0.85],
            ),
            (
                "2*(a+b - 2*b*z^2)^(-1/2)*(1-z^2)^(-1/2)",
                vec![-0.8, -0.3, 0.5, 0.9],
            ),
            (
                "2*(a+b - 2*b*z^2)^(-3/2)*(1-z^2)^(-1/2)",
                vec![-0.7, -0.35, 0.45, 0.9],
            ),
            (
                "2*(a+b - 2*b*z^2)^(1/2)*(1-z^2)^(-1/2)",
                vec![-0.8, -0.4, 0.3, 0.85],
            ),
        ] {
            let e = parse(&ctx, src);
            assert_antiderivative_num(&ctx, e, s, &env, &samples);
        }
    }

    // ---------------- declines ----------------

    #[test]
    fn declines_unsupported_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Cubic radicand: needs the cubic F/E transform (not implemented).
        declines(&ctx, "1/sqrt(x^3+1)");
        declines(&ctx, "1/sqrt((x - 1)*(x + 1)*(x + 2))");
        // Quintic radicand.
        declines(&ctx, "1/sqrt(x^5+1)");
        // Odd numerator: elementary, not elliptic.
        declines(&ctx, "x/sqrt(1-x^4)");
        declines(&ctx, "x^3/sqrt(1-x^4)");
        // Third kind: simple pole away from the branch points.
        declines(&ctx, "1/((x-c)*sqrt(1-x^4))");
        // Non-Legendre even quartic: (1+x²)(1+2x²) needs a Möbius
        // normalisation.
        declines(&ctx, "1/sqrt((1+x^2)*(1+2*x^2))");
        declines(&ctx, "1/sqrt(1+x^4)");
        // Rational (non-radical) and transcendental integrands.
        declines(&ctx, "1/(1+x^2)");
        declines(&ctx, "exp(x)/sqrt(1-x^4)");
        declines(&ctx, "sin(x)");
        declines(&ctx, "1/sqrt(1-x^2)");
        declines(&ctx, "x^(1/3)");
        declines(&ctx, "1/sqrt(sqrt(1-x^4))");
    }

    #[test]
    fn accepts_direct_identity_variable_shapes() {
        // The engine is also reachable on plain algebraic inputs (var = x),
        // where the "substitution" is the identity.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = parse(&ctx, "1/sqrt((1-x^2)*(1-(5^-1)*x^2))");
        let r = integrate_elliptic(&ctx, e, Symbol::new("x")).expect("solved");
        assert!(r.to_string().contains("EllipticF"), "got {r}");
    }

    // ---------------- stress ----------------

    #[test]
    fn stress_stable_outcomes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let solid = parse(&ctx, "1/sqrt((1-x^2)*(1-(4^-1)*x^2))");
        let mut first: Option<String> = None;
        for i in 0..300 {
            let r = integrate_elliptic(&ctx, solid, Symbol::new("x"));
            let s = r.map(|a| a.to_string());
            if i == 0 {
                first = s.clone();
            }
            assert_eq!(s, first, "non-deterministic outcome at iteration {i}");
        }
        assert!(first.is_some(), "stable outcome must be a solve");
    }
}
