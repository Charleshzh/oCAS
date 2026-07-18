//! Rational function integration (Bronstein, *Symbolic Integration I*, ch. 2).
//!
//! [`integrate_rational`] computes the integral of a rational function of
//! the integration variable over `ℚ`:
//!
//! 1. the polynomial part by termwise integration,
//! 2. the rational part by Hermite reduction,
//! 3. the logarithmic part by the logarithmic-derivative identity
//!    (`c·f'/f → c·log(f)`), completing the square (degree-2 denominators,
//!    yielding `log` or `atan`), or the Rothstein–Trager resultant for
//!    general squarefree denominators.
//!
//! Terms whose Rothstein–Trager resultant does not split into rational
//! roots cannot be expressed without `RootSum`-style notation; those terms
//! are returned in the unevaluated form `Integral(term, var)`.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, Symbol};
use ocas_domain::{Domain, Integer, Rational, RationalDomain};
use ocas_poly::{DenseUnivariatePolynomial, Lex, SparseMultivariatePolynomial};
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;

use crate::rules::calculus_rules;
use crate::tower::convert::{
    GeneratorField, atom_to_rational, rational_const_to_atom, rational_to_atom,
};

type DPoly = DenseUnivariatePolynomial<RationalDomain>;
type Sparse = SparseMultivariatePolynomial<RationalDomain, Lex>;

/// Integrate a rational function of `var` over `ℚ`.
///
/// Returns `None` when `expr` is not a rational function of `var` alone
/// (other variables or function applications are present). Terms that
/// would require algebraic numbers beyond `√` are returned in the
/// unevaluated form `Integral(term, var)` inside the sum.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::integral::rational::integrate_rational;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// // ∫ x^-1 dx = log(x)
/// let result = integrate_rational(&ctx, ctx.pow(x, ctx.num(-1)), Symbol::new("x"));
/// assert_eq!(result.unwrap().to_string(), "log(x)");
/// ```
pub fn integrate_rational<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let gens = [x];
    let rf = atom_to_rational(expr, &gens)?;
    let num = sparse_to_dense(&rf.numerator);
    let den = sparse_to_dense(&rf.denominator);
    let raw = integrate_dpoly(ctx, &num, &den, x)?;
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let after_default = simplify(ctx, raw, &default_rules, 20);
    let after_calc = simplify(ctx, after_default, &calc_rules, 10);
    Some(normalize(ctx, after_calc))
}

/// Core integrator for `num / den` over `ℚ[x]`, building atoms over `x`.
pub(crate) fn integrate_dpoly<'a>(
    ctx: &'a AtomArena<'a>,
    num: &DPoly,
    den: &DPoly,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let s = integrate_dpoly_structured(ctx, num, den, x)?;
    let mut parts: Vec<Atom> = Vec::new();
    if !s.poly.is_zero() {
        parts.push(poly_atom(ctx, &s.poly, x)?);
    }
    if let Some((g_num, g_den)) = &s.rational {
        let g = GeneratorField::from_num_den(dense_to_sparse(g_num), dense_to_sparse(g_den));
        parts.push(rational_to_atom(ctx, &g, &[x])?);
    }
    for (c, v) in &s.log_terms {
        let v_atom = poly_atom(ctx, v, x)?;
        let log = ctx.fun("log", &[v_atom]);
        parts.push(scale_atom(ctx, c, log)?);
    }
    parts.extend(s.extra_atoms.iter().copied());
    Some(match parts.len() {
        0 => ctx.num(0),
        1 => parts[0],
        _ => ctx.add(&parts),
    })
}

/// A structured rational-function integral (used by the Risch tower's
/// base level, which needs the pieces rather than a assembled atom).
pub(crate) struct DPolyIntegral<'a> {
    /// Integral of the polynomial part.
    pub poly: DPoly,
    /// Hermite rational part `(num, den)`, if any.
    pub rational: Option<(DPoly, DPoly)>,
    /// Logarithmic terms `c·log(v)` with `v ∈ ℚ[x]`.
    pub log_terms: Vec<(Rational, DPoly)>,
    /// Other terms (atan, radical logs, unevaluated `Integral`s).
    pub extra_atoms: Vec<Atom<'a>>,
}

/// Structured core of [`integrate_dpoly`].
pub(crate) fn integrate_dpoly_structured<'a>(
    ctx: &'a AtomArena<'a>,
    num: &DPoly,
    den: &DPoly,
    x: Atom<'a>,
) -> Option<DPolyIntegral<'a>> {
    if den.is_zero() {
        return None;
    }
    let (p, r) = num.div_rem(den)?;
    let mut out = DPolyIntegral {
        poly: DPoly::from_coeffs(RationalDomain, vec![]),
        rational: None,
        log_terms: Vec::new(),
        extra_atoms: Vec::new(),
    };
    if !p.is_zero() {
        out.poly = poly_integrate(&p);
    }
    if !r.is_zero() {
        let (g_num, g_den, a1, d1) = hermite_reduce(&r, den);
        if !g_num.is_zero() {
            out.rational = Some((g_num, g_den));
        }
        if !a1.is_zero() {
            for term in log_part(ctx, &a1, &d1, x)? {
                match term {
                    DLogOrAtom::Log(c, v) => out.log_terms.push((c, v)),
                    DLogOrAtom::Atom(a) => out.extra_atoms.push(a),
                }
            }
        }
    }
    Some(out)
}

// ------------------------------------------------------------------
//  Polynomial part
// ------------------------------------------------------------------

/// Termwise integral of a polynomial: `∫ Σ cₖ xᵏ = Σ cₖ/(k+1) xᵏ⁺¹`.
pub(crate) fn poly_integrate(p: &DPoly) -> DPoly {
    let dom = RationalDomain;
    let mut coeffs = vec![dom.zero()];
    for (k, c) in p.coeffs().iter().enumerate() {
        let scale = Rational::new(1, (k + 1) as i64);
        coeffs.push(dom.mul(c, &scale));
    }
    DPoly::from_coeffs(RationalDomain, coeffs)
}

// ------------------------------------------------------------------
//  Hermite reduction
// ------------------------------------------------------------------

/// Hermite reduction of a proper fraction `a / d` with `d` monic.
///
/// Returns `(g_num, g_den, a1, d1)` such that `a/d = (g_num/g_den)' + a1/d1`
/// with `d1` squarefree.
fn hermite_reduce(a: &DPoly, d: &DPoly) -> (DPoly, DPoly, DPoly, DPoly) {
    let factors = d.square_free_factorization();
    let m = factors.iter().map(|&(_, k)| k).max().unwrap_or(1);
    if m <= 1 {
        return (a.zero(), a.one(), a.clone(), d.clone());
    }
    // Pick a factor v of maximal multiplicity m; write d = u * v^m.
    let (v, _) = factors
        .iter()
        .rev()
        .find(|&(_, k)| *k == m)
        .expect("max multiplicity factor exists");
    let vm = v.pow(m as u32);
    let (u, rem) = d.div_rem(&vm).expect("v^m divides d");
    debug_assert!(rem.is_zero());

    // Solve s*B + t*v = a with B = u*v' and deg s < deg v (gcd(B, v) = 1).
    let b = u.mul(&v.derivative());
    let (g, s0, _t0) = b.extended_gcd_poly(v);
    debug_assert!(g.is_one());
    let (_, s) = a.mul(&s0).div_rem(v).expect("field division");
    let (t, rem) = a.sub(&s.mul(&b)).div_rem(v).expect("exact quotient");
    debug_assert!(rem.is_zero());

    // ∫ a/(u v^m) = -(s/(m-1))/v^(m-1) + ∫ (t + u*(s/(m-1))') / (u v^(m-1))
    let inv_m1 = Rational::new(1, (m - 1) as i64);
    let s_scaled = s.mul_scalar(&inv_m1);
    let vm1 = v.pow((m - 1) as u32);
    let new_a = t.add(&u.mul(&s_scaled.derivative()));
    let new_d = u.mul(&vm1);
    let (g2n, g2d, a1, d1) = hermite_reduce(&new_a, &new_d);

    // Combine rational parts: -s_scaled/vm1 + g2n/g2d.
    let g_num = s_scaled.neg().mul(&g2d).add(&g2n.mul(&vm1));
    let g_den = vm1.mul(&g2d);
    (g_num, g_den, a1, d1)
}

// ------------------------------------------------------------------
//  Logarithmic part
// ------------------------------------------------------------------

/// A logarithmic-part term: either a structured `c·log(v)` or an opaque
/// atom (atan, radical forms, unevaluated integrals).
pub(crate) enum DLogOrAtom<'a> {
    /// `c·log(v)` with `v ∈ ℚ[x]`.
    Log(Rational, DPoly),
    /// Any other term as an atom.
    Atom(Atom<'a>),
}

/// Integrate `a / d` with `d` monic and squarefree. Returns the summands;
/// may include an unevaluated `Integral` term.
fn log_part<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let dom = RationalDomain;
    // Case 1: a == c·d' → c·log(d) (covers deg d == 1).
    let dp = d.derivative();
    if !dp.is_zero() && a.degree() == dp.degree() {
        let c = dom.div(&a.lcoeff(), &dp.lcoeff())?;
        if *a == dp.mul_scalar(&c) {
            return Some(vec![DLogOrAtom::Log(c, d.clone())]);
        }
    }
    // Case 2: degree-2 denominator → completing the square.
    if d.degree() == Some(2) {
        return complete_square(ctx, a, d, x);
    }
    // Case 3: Rothstein–Trager for general squarefree denominators.
    rothstein_trager(ctx, a, d, x)
}

/// Integrate `(A·x + B) / (x² + b·x + c)` with `x² + b·x + c` squarefree.
fn complete_square<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let dom = RationalDomain;
    let c0 = d.coeffs().first().cloned().unwrap_or_else(|| dom.zero());
    let c1 = d.coeffs().get(1).cloned().unwrap_or_else(|| dom.zero());
    let big_a = a.coeffs().get(1).cloned().unwrap_or_else(|| dom.zero());
    let big_b = a.coeffs().first().cloned().unwrap_or_else(|| dom.zero());
    let half = Rational::new(1, 2);

    let mut out = Vec::new();
    // (A/2)·log(d)
    if !dom.is_zero(&big_a) {
        let coeff = dom.mul(&big_a, &half);
        out.push(DLogOrAtom::Log(coeff, d.clone()));
    }
    // Remaining constant numerator: (B - A·b/2)·∫dx/(x²+b·x+c).
    let rem = dom.sub(&big_b, &dom.mul(&big_a, &dom.mul(&c1, &half)));
    if dom.is_zero(&rem) {
        return Some(out);
    }
    // Discriminant Δ = b² - 4c (nonzero: d is squarefree).
    let delta = dom.sub(&dom.mul(&c1, &c1), &dom.mul(&Rational::new(4, 1), &c0));
    let lin = linear_atom(ctx, &Rational::new(2, 1), &c1, x)?; // 2x + b
    if rat_is_negative(&delta) {
        // atan branch: ∫ = rem·(2/s)·atan((2x+b)/s), s = √(4c-b²).
        let s2 = dom.neg(&delta);
        match sqrt_positive_rational(ctx, &s2)? {
            Sqrt::Rat(sr) => {
                let sr_inv = dom.inv(&sr)?;
                let arg = linear_atom(
                    ctx,
                    &dom.mul(&Rational::new(2, 1), &sr_inv),
                    &dom.mul(&c1, &sr_inv),
                    x,
                )?;
                let coeff = dom.mul(&rem, &dom.mul(&Rational::new(2, 1), &sr_inv));
                out.push(DLogOrAtom::Atom(scale_atom(
                    ctx,
                    &coeff,
                    ctx.fun("atan", &[arg]),
                )?));
            }
            Sqrt::Rad(s) => {
                let arg = ctx.mul(&[lin, ctx.pow(s, ctx.num(-1))]);
                let atan = ctx.fun("atan", &[arg]);
                let two_rem = dom.mul(&rem, &Rational::new(2, 1));
                out.push(DLogOrAtom::Atom(ctx.mul(&[
                    rational_const_to_atom(ctx, &two_rem)?,
                    ctx.pow(s, ctx.num(-1)),
                    atan,
                ])));
            }
        }
    } else {
        // log branch: ∫ = rem·(1/s)·log((2x+b-s)/(2x+b+s)), s = √Δ.
        match sqrt_positive_rational(ctx, &delta)? {
            Sqrt::Rat(sr) => {
                let num_atom = linear_atom(ctx, &Rational::new(2, 1), &dom.sub(&c1, &sr), x)?;
                let den_atom = linear_atom(ctx, &Rational::new(2, 1), &dom.add(&c1, &sr), x)?;
                let ratio = ctx.mul(&[num_atom, ctx.pow(den_atom, ctx.num(-1))]);
                let log = ctx.fun("log", &[ratio]);
                let coeff = dom.mul(&rem, &dom.inv(&sr)?);
                out.push(DLogOrAtom::Atom(scale_atom(ctx, &coeff, log)?));
            }
            Sqrt::Rad(s) => {
                let num_atom = ctx.add(&[lin, ctx.mul(&[ctx.num(-1), s])]);
                let den_atom = ctx.add(&[lin, s]);
                let ratio = ctx.mul(&[num_atom, ctx.pow(den_atom, ctx.num(-1))]);
                let log = ctx.fun("log", &[ratio]);
                out.push(DLogOrAtom::Atom(ctx.mul(&[
                    rational_const_to_atom(ctx, &rem)?,
                    ctx.pow(s, ctx.num(-1)),
                    log,
                ])));
            }
        }
    }
    Some(out)
}

/// Rothstein–Trager: `∫ a/d = Σᵢ cᵢ·log(gcd(d, a - cᵢ·d'))` where the `cᵢ`
/// are the roots of `R(t) = resultantₓ(d, a - t·d')`.
///
/// `R(t)` is recovered by interpolation (resultants are univariate over
/// `ℚ`, so `t` cannot appear symbolically in the coefficients). Falls back
/// to an unevaluated `Integral` when `R` does not split over `ℚ`.
fn rothstein_trager<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    x: Atom<'a>,
) -> Option<Vec<DLogOrAtom<'a>>> {
    let n = d.degree()?;
    let dp = d.derivative();

    // Interpolate R(t) = resultant(d, a - t·d') at t = 0, …, n.
    let mut points = Vec::with_capacity(n + 1);
    for j in 0..=(n as i64) {
        let tj = Rational::new(j, 1);
        let shifted = a.sub(&dp.mul_scalar(&tj));
        let val = d.resultant(&shifted);
        points.push((tj, val));
    }
    let rt = lagrange_interpolate(&points);

    let (roots, fully_split) = rational_roots(&rt)?;
    if !fully_split {
        return Some(vec![DLogOrAtom::Atom(integral_fallback(ctx, a, d, x)?)]);
    }

    let mut out = Vec::new();
    for c in &roots {
        let v = d.gcd(&a.sub(&dp.mul_scalar(c)));
        if v.degree().unwrap_or(0) == 0 {
            continue;
        }
        out.push(DLogOrAtom::Log(c.clone(), v));
    }
    if out.is_empty() {
        return Some(vec![DLogOrAtom::Atom(integral_fallback(ctx, a, d, x)?)]);
    }
    Some(out)
}

// ------------------------------------------------------------------
//  Helpers
// ------------------------------------------------------------------

/// Lagrange interpolation over `ℚ` from point values.
fn lagrange_interpolate(points: &[(Rational, Rational)]) -> DPoly {
    let dom = RationalDomain;
    let mut result = DPoly::from_coeffs(RationalDomain, vec![]);
    for (j, (xj, yj)) in points.iter().enumerate() {
        let mut basis = DPoly::from_coeffs(RationalDomain, vec![dom.one()]);
        let mut denom = dom.one();
        for (k, (xk, _)) in points.iter().enumerate() {
            if k == j {
                continue;
            }
            basis = basis.mul(&DPoly::from_coeffs(
                RationalDomain,
                vec![dom.neg(xk), dom.one()],
            ));
            denom = dom.mul(&denom, &dom.sub(xj, xk));
        }
        let scale = dom.div(yj, &denom).expect("distinct interpolation nodes");
        result = result.add(&basis.mul_scalar(&scale));
    }
    result
}

/// Distinct rational roots of a polynomial over `ℚ`, plus whether the
/// polynomial splits completely into linear factors over `ℚ`.
fn rational_roots(f: &DPoly) -> Option<(Vec<Rational>, bool)> {
    let dom = RationalDomain;
    if f.is_zero() || f.degree()? == 0 {
        return Some((Vec::new(), true));
    }
    // Clear denominators to get an integer polynomial.
    let mut lcm: i64 = 1;
    for c in f.coeffs() {
        lcm = num_integer::lcm(lcm, c.denom().to_i64()?);
    }
    let zcoeffs: Option<Vec<Integer>> = f
        .coeffs()
        .iter()
        .map(|c| {
            let scaled = dom.mul(c, &Rational::new(lcm, 1));
            scaled.numer().to_i64().map(Integer::from)
        })
        .collect();
    let zpoly = DenseUnivariatePolynomial::from_coeffs(ocas_domain::IntegerDomain, zcoeffs?);
    let primitive = zpoly.primitive_part();
    let factors = primitive.factor();

    let mut roots = Vec::new();
    let mut split = true;
    for (fac, _mult) in &factors {
        if fac.degree() == Some(1) {
            let c0 = fac.coeffs().first()?.to_i64()?;
            let c1 = fac.coeffs().get(1)?.to_i64()?;
            roots.push(Rational::new(-c0, c1));
        } else {
            split = false;
        }
    }
    Some((roots, split))
}

/// Square root of a positive rational: exact when the squarefree numerator
/// is a perfect square, otherwise a radical atom `√(p·q) / q`.
enum Sqrt<'a> {
    Rat(Rational),
    Rad(Atom<'a>),
}

fn sqrt_positive_rational<'a>(ctx: &'a AtomArena<'a>, r: &Rational) -> Option<Sqrt<'a>> {
    let p = r.numer().to_i64()?;
    let q = r.denom().to_i64()?;
    let n = (p as i128).checked_mul(q as i128)?;
    let m = isqrt_i128(n);
    if m * m == n {
        let m = i64::try_from(m).ok()?;
        Some(Sqrt::Rat(Rational::new(m, q)))
    } else {
        let n = i64::try_from(n).ok()?;
        Some(Sqrt::Rad(
            ctx.pow(ctx.num(n), ctx.pow(ctx.num(2), ctx.num(-1))),
        ))
    }
}

fn isqrt_i128(n: i128) -> i128 {
    if n <= 0 {
        return 0;
    }
    let mut x = (n as f64).sqrt() as i128 + 2;
    while x * x > n {
        x = (x + n / x) / 2;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

fn rat_is_negative(r: &Rational) -> bool {
    use num_traits::Signed;
    r.inner().is_negative()
}

/// `coeff · atom`, eliding a unit coefficient.
fn scale_atom<'a>(ctx: &'a AtomArena<'a>, coeff: &Rational, atom: Atom<'a>) -> Option<Atom<'a>> {
    if RationalDomain.is_one(coeff) {
        return Some(atom);
    }
    Some(ctx.mul(&[rational_const_to_atom(ctx, coeff)?, atom]))
}

/// `a·x + b` with `a ≠ 0`, eliding a zero constant term.
fn linear_atom<'a>(
    ctx: &'a AtomArena<'a>,
    a: &Rational,
    b: &Rational,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let dom = RationalDomain;
    let ax = if dom.is_one(a) {
        x
    } else {
        ctx.mul(&[rational_const_to_atom(ctx, a)?, x])
    };
    if dom.is_zero(b) {
        return Some(ax);
    }
    Some(ctx.add(&[ax, rational_const_to_atom(ctx, b)?]))
}

fn poly_atom<'a>(ctx: &'a AtomArena<'a>, p: &DPoly, x: Atom<'a>) -> Option<Atom<'a>> {
    rational_to_atom(
        ctx,
        &GeneratorField::from_polynomial(dense_to_sparse(p)),
        &[x],
    )
}

fn integral_fallback<'a>(
    ctx: &'a AtomArena<'a>,
    a: &DPoly,
    d: &DPoly,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let rf = GeneratorField::from_num_den(dense_to_sparse(a), dense_to_sparse(d));
    let expr = rational_to_atom(ctx, &rf, &[x])?;
    Some(ctx.fun("Integral", &[expr, x]))
}

fn sparse_to_dense(p: &Sparse) -> DPoly {
    debug_assert_eq!(p.n_vars(), 1);
    let deg = p.degree_in(0);
    let mut coeffs = vec![RationalDomain.zero(); deg + 1];
    for (exp, coeff) in p.terms_ref() {
        coeffs[exp[0]] = coeff.clone();
    }
    DPoly::from_coeffs(RationalDomain, coeffs)
}

fn dense_to_sparse(p: &DPoly) -> Sparse {
    let terms = p
        .coeffs()
        .iter()
        .enumerate()
        .filter(|&(_, c)| !RationalDomain.is_zero(c))
        .map(|(i, c)| (vec![i], c.clone()))
        .collect();
    Sparse::from_terms(RationalDomain, 1, terms)
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;
    use crate::derivative::diff;
    use ocas_atom::AtomNode;

    /// Verify that `diff(result) == integrand` as canonical rational functions.
    fn assert_antiderivative<'a>(
        ctx: &'a AtomArena<'a>,
        result: Atom<'a>,
        integrand: Atom<'a>,
        x: Atom<'a>,
        var: Symbol,
    ) {
        let d = diff(ctx, result, var);
        let lhs = atom_to_rational(d, &[x]).expect("derivative is rational in x");
        let rhs = atom_to_rational(integrand, &[x]).expect("integrand is rational in x");
        assert_eq!(lhs, rhs, "d/dx(result) != integrand");
    }

    /// Minimal f64 evaluator for numeric checks on radical-containing results.
    fn eval_f64(atom: Atom, x: f64) -> f64 {
        match atom.node() {
            AtomNode::Num(n) => *n as f64,
            AtomNode::Var(_) => x,
            AtomNode::Add(args) => args.iter().map(|a| eval_f64(*a, x)).sum(),
            AtomNode::Mul(args) => args.iter().map(|a| eval_f64(*a, x)).product(),
            AtomNode::Pow(b, e) => eval_f64(*b, x).powf(eval_f64(*e, x)),
            AtomNode::Fun(name, args) => match name.as_str() {
                // |u| keeps the real evaluation finite; log|u| and log(u)
                // differ by a constant on each interval, so derivatives match.
                "log" => eval_f64(args[0], x).abs().ln(),
                "atan" => eval_f64(args[0], x).atan(),
                "sin" => eval_f64(args[0], x).sin(),
                "cos" => eval_f64(args[0], x).cos(),
                "exp" => eval_f64(args[0], x).exp(),
                other => panic!("unsupported function in test eval: {other}"),
            },
        }
    }

    /// Numeric finite-difference check that `result'` matches `integrand`.
    fn assert_numeric_antiderivative(result: Atom, integrand: Atom) {
        let h = 1e-6;
        for px in [0.7f64, 1.3, 2.1] {
            let approx = (eval_f64(result, px + h) - eval_f64(result, px - h)) / (2.0 * h);
            let exact = eval_f64(integrand, px);
            assert!(
                (approx - exact).abs() < 1e-4,
                "numeric derivative mismatch at x={px}: {approx} != {exact}"
            );
        }
    }

    fn integrate_str<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
        integrate_rational(ctx, expr, Symbol::new("x")).expect("rational integral")
    }

    #[test]
    fn inverse_gives_log() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let result = integrate_str(&ctx, ctx.pow(x, ctx.num(-1)));
        assert_eq!(result.to_string(), "log(x)");
    }

    #[test]
    fn polynomial_part() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ (x^2 + 3x + 2) dx
        let expr = ctx.add(&[
            ctx.pow(x, ctx.num(2)),
            ctx.mul(&[ctx.num(3), x]),
            ctx.num(2),
        ]);
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn atan_over_quadratic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x^2 + 1) dx = atan(x)
        let expr = ctx.pow(ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(1)]), ctx.num(-1));
        let result = integrate_str(&ctx, expr);
        assert_eq!(result.to_string(), "atan(x)");
    }

    #[test]
    fn log_over_quadratic_irrational_discriminant() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x^2 - 2) dx = (1/(2√2))·log((x-√2)/(x+√2))
        let expr = ctx.pow(ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(-2)]), ctx.num(-1));
        let result = integrate_str(&ctx, expr);
        assert!(result.to_string().contains("log"));
        assert_numeric_antiderivative(result, expr);
    }

    #[test]
    fn log_over_quadratic_rational_roots() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x^2 - 1) dx = (1/2)·log((x-1)/(x+1))
        let expr = ctx.pow(ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(-1)]), ctx.num(-1));
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn logarithmic_derivative_identity() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ (2x + 3)/(x^2 + 3x + 5) dx = log(x^2 + 3x + 5)
        let num = ctx.add(&[ctx.mul(&[ctx.num(2), x]), ctx.num(3)]);
        let den = ctx.add(&[
            ctx.pow(x, ctx.num(2)),
            ctx.mul(&[ctx.num(3), x]),
            ctx.num(5),
        ]);
        let expr = ctx.mul(&[num, ctx.pow(den, ctx.num(-1))]);
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn hermite_repeated_factor() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x+1)^2 dx = -(x+1)^-1
        let expr = ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(-2));
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn hermite_mixed_multiplicities() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x/( (x+1)^2 (x+2) ) dx — repeated factor times simple factor.
        let den = ctx.mul(&[
            ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(2)),
            ctx.add(&[x, ctx.num(2)]),
        ]);
        let expr = ctx.mul(&[x, ctx.pow(den, ctx.num(-1))]);
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn rothstein_trager_cubic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ 1/(x^3 + x) dx: squarefree cubic, roots 0, ±i — splits over ℚ
        // only partially; RT resultant has one rational root (0)... The
        // integrand = 1/x - x/(x^2+1), so a fully rational-form answer
        // exists via partial fractions with rational terms only after
        // combining conjugates. Verify by differentiation.
        let expr = ctx.pow(ctx.add(&[ctx.pow(x, ctx.num(3)), x]), ctx.num(-1));
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }

    #[test]
    fn non_rational_returns_none() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // sin(x) is not a rational function of x.
        let expr = ctx.fun("sin", &[x]);
        assert!(integrate_rational(&ctx, expr, Symbol::new("x")).is_none());
    }

    #[test]
    fn higher_multiplicity_hermite() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ (x^2 + 1)/(x-1)^3 dx — multiplicity 3 denominator.
        let num = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(1)]);
        let den = ctx.pow(ctx.add(&[x, ctx.num(-1)]), ctx.num(3));
        let expr = ctx.mul(&[num, ctx.pow(den, ctx.num(-1))]);
        let result = integrate_str(&ctx, expr);
        assert_antiderivative(&ctx, result, expr, x, Symbol::new("x"));
    }
}
