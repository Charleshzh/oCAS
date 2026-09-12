//! Special-function antiderivatives (the Meijer-G fallback's endpoint).
//!
//! When the Risch algorithm proves an integral has no elementary
//! antiderivative, many common cases still have closed forms in terms of
//! special functions. Rather than routing through the full Meijer
//! G-function machinery (which requires hypergeometric series, gamma
//! functions, and Slater expansions oCAS does not yet provide), the
//! antiderivatives of the standard non-elementary integrals are encoded
//! directly:
//!
//! - `exp(±x²)` → `erf` / `erfi`
//! - `exp(x)/x` → `Ei`
//! - `sin(x)/x` → `Si`,  `cos(x)/x` → `Ci`,  `cosh(x)/x` → `Chi`, `sinh(x)/x` → `Shi`
//! - `sin(x²)`, `cos(x²)` → Fresnel `S` / `C`
//!
//! 0.27.3 Phase B adds the *reduction* families that turn an arbitrary
//! polynomial (or a monomial denominator power) multiplying a special head
//! into a closed form, using only identities that are checked in this file:
//!
//! - **Polynomial × `F(u)`** with `u = a + b·x` and `F'` elementary
//!   (`erf`, `erfc`, `erfi`, `Si`, `Ci`, `Shi`, `Chi`, `Ei`): substitute
//!   `u = a + b·x`, expand `(u − a)^k`, then repeatedly use
//!   `∫u^j F(u)du = u^(j+1)F(u)/(j+1) − (1/(j+1))∫u^(j+1)F'(u)du`, whose
//!   residual is elementary for every head above.
//! - **Polynomial × `Ei(n, u)`** (integer order `n`): the inverse recurrence
//!   `E_n(u) = (e^−u − n·E_(n+1)(u))/u` reduces the polynomial degree one step
//!   at a time and terminates on the elementary form of `E_0`; the base case
//!   `∫E_n(u)du = −E_(n+1)(u)` (all orders) covers a bare kernel.
//! - **`F(b·x)/x^m`** (`m ≥ 2`, `F = Ei` or `Ei(n, ·)`): the parts descent
//!   `∫F(u)u^−m du = −F(u)u^(1−m)/(m−1) + (1/(m−1))∫F'(u)u^(1−m)du` walks
//!   down to `∫e^(±u)u^−j du`, which terminates in `Ei(±u)` plus
//!   elementary terms.
//! - **Polynomial × `F(u)²`** (`erf`, `erfi`, `Si`, `Ci`, `Shi`, `Chi`):
//!   parts with `∫u^j F² du = u^(j+1)F²/(j+1) − (2/(j+1))∫u^(j+1)F·F' du`,
//!   and the remaining `F·F'` integrals are reduced by the Gaussian
//!   moment recursion (for `erf`/`erfi`) or the coupled `∫u^j F·cos/sin`
//!   recursion (for the integral functions), with all base cases verified
//!   by differentiation in the unit tests.
//!
//! The entries follow the same definitions as SymPy, so results compare
//! equal against `sympy.integrate`. Every family returns `None` unless the
//! integrand really contains one of the special heads it handles, and every
//! reduction chain carries a hard step budget, so the stage can never hang.
//! It runs *before* the rule table, so breadth is deliberately traded for
//! precision: a shape that is not recognized is declined, never guessed.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::{is_constant, linear_form};

/// Try to integrate `expr` in terms of special functions.
///
/// Only the patterns listed in the module docs are recognized; returns
/// `None` otherwise (the caller then emits the unevaluated form).
pub(crate) fn special_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    // All table entries need the integrand normalized as a product; we
    // inspect it as a flat list of factors. The flattening mirrors
    // `normalize`'s `Mul` handling so the stage also works on the residual
    // shapes other reduction formulas hand it directly.
    let mut factors: Vec<Atom> = Vec::new();
    flat_factors(expr, &mut factors);
    erf_family(ctx, &factors, x)
        .or_else(|| ei_family(ctx, &factors, x))
        .or_else(|| trig_integral_family(ctx, &factors, x))
        .or_else(|| fresnel_family(ctx, &factors, x))
        .or_else(|| reduction_families(ctx, &factors, x))
}

/// Collect the flat factor list of a product, descending into nested `Mul`s.
fn flat_factors<'a>(expr: Atom<'a>, out: &mut Vec<Atom<'a>>) {
    match expr.node() {
        AtomNode::Mul(args) => {
            for a in args.iter() {
                flat_factors(*a, out);
            }
        }
        _ => out.push(expr),
    }
}

// ------------------------------------------------------------------
//  Pattern helpers
// ------------------------------------------------------------------

/// Match `exp(u)`; return `u`.
fn as_exp(f: Atom) -> Option<Atom> {
    if let AtomNode::Fun(name, args) = f.node() {
        if name.as_str() == "exp" && args.len() == 1 {
            return Some(args[0]);
        }
    }
    None
}

/// Match `x^-1`.
fn is_x_inv<'a>(f: Atom<'a>, x: Atom<'a>) -> bool {
    matches!(f.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(-1)))
}

/// Match `c·x^2` (with optional sign): returns the coefficient atom `c`.
fn as_quadratic<'a>(u: Atom<'a>, x: Atom<'a>) -> Option<Atom<'a>> {
    // x^2
    if matches!(u.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(2))) {
        return None; // coefficient 1 handled by caller
    }
    // c·x^2 or -x^2
    if let AtomNode::Mul(factors) = u.node() {
        if factors.len() == 2 {
            if matches!(factors[1].node(), AtomNode::Pow(b, e)
                if *b == x && matches!(e.node(), AtomNode::Num(2)))
            {
                return Some(factors[0]);
            }
        }
    }
    None
}

/// Whether `u` is exactly `x²`, or `(-x)²` which normalizes to it.
fn is_x_squared<'a>(u: Atom<'a>, x: Atom<'a>) -> bool {
    if matches!(u.node(), AtomNode::Pow(b, e) if *b == x && matches!(e.node(), AtomNode::Num(2))) {
        return true;
    }
    // (-x)^2 or (-1·x)^2 → x²
    if let AtomNode::Pow(b, e) = u.node() {
        if matches!(e.node(), AtomNode::Num(2)) {
            if let AtomNode::Mul(factors) = b.node() {
                let all_neg_one_or_x = factors
                    .iter()
                    .all(|f| matches!(f.node(), AtomNode::Num(-1)) || *f == x);
                let has_x = factors.contains(&x);
                return all_neg_one_or_x && has_x;
            }
        }
    }
    false
}

/// Whether `u` is exactly `x`.
fn is_x<'a>(u: Atom<'a>, x: Atom<'a>) -> bool {
    u == x
}

// ------------------------------------------------------------------
//  erf family: exp(-x²), exp(c·x²)
// ------------------------------------------------------------------

fn erf_family<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], x: Atom<'a>) -> Option<Atom<'a>> {
    if factors.len() != 1 {
        return None;
    }
    let u = as_exp(factors[0])?;
    // exp(-x²) → (√π/2)·erf(x)
    if is_x_squared(u, x) {
        // exp(+x²) → (√π/2)·erfi(x)
        let sqrt_pi = ctx.fun("sqrt", &[ctx.var("pi")]);
        let erfi = ctx.fun("erfi", &[x]);
        return Some(ctx.mul(&[sqrt_pi, ctx.pow(ctx.num(2), ctx.num(-1)), erfi]));
    }
    // exp(c·x²) with negative c: √π/(2√(-c))·erf(√(-c)·x)
    if let Some(c) = as_quadratic(u, x) {
        let neg_c = ctx.mul(&[ctx.num(-1), c]);
        let root = ctx.fun("sqrt", &[neg_c]);
        let sqrt_pi = ctx.fun("sqrt", &[ctx.var("pi")]);
        let erf = ctx.fun("erf", &[ctx.mul(&[root, x])]);
        return Some(ctx.mul(&[
            sqrt_pi,
            ctx.pow(ctx.mul(&[ctx.num(2), root]), ctx.num(-1)),
            erf,
        ]));
    }
    None
}

// ------------------------------------------------------------------
//  Ei family: exp(x)/x, exp(c·x)/x
// ------------------------------------------------------------------

fn ei_family<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], x: Atom<'a>) -> Option<Atom<'a>> {
    if factors.len() != 2 {
        return None;
    }
    let (exp_f, inv_f) = if as_exp(factors[0]).is_some() {
        (factors[0], factors[1])
    } else if as_exp(factors[1]).is_some() {
        (factors[1], factors[0])
    } else {
        return None;
    };
    if !is_x_inv(inv_f, x) {
        return None;
    }
    let u = as_exp(exp_f)?;
    // exp(x)/x → Ei(x)
    if is_x(u, x) {
        return Some(ctx.fun("Ei", &[x]));
    }
    // exp(c·x)/x → Ei(c·x) for constant c.
    if let AtomNode::Mul(mf) = u.node() {
        if mf.len() == 2 && mf[1] == x && matches!(mf[0].node(), AtomNode::Num(_)) {
            return Some(ctx.fun("Ei", &[u]));
        }
    }
    None
}

// ------------------------------------------------------------------
//  Si/Ci/Shi/Chi family: sin(x)/x, cos(x)/x, sinh(x)/x, cosh(x)/x
// ------------------------------------------------------------------

fn trig_integral_family<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    if factors.len() != 2 {
        return None;
    }
    let (fun_f, inv_f) = if is_x_inv(factors[1], x) {
        (factors[0], factors[1])
    } else if is_x_inv(factors[0], x) {
        (factors[1], factors[0])
    } else {
        return None;
    };
    let _ = inv_f;
    let AtomNode::Fun(name, args) = fun_f.node() else {
        return None;
    };
    if args.len() != 1 || !is_x(args[0], x) {
        return None;
    }
    let target = match name.as_str() {
        "sin" => "Si",
        "cos" => "Ci",
        "sinh" => "Shi",
        "cosh" => "Chi",
        _ => return None,
    };
    Some(ctx.fun(target, &[x]))
}

// ------------------------------------------------------------------
//  Fresnel family: sin(x²), cos(x²)
// ------------------------------------------------------------------

fn fresnel_family<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    if factors.len() != 1 {
        return None;
    }
    let AtomNode::Fun(name, args) = factors[0].node() else {
        return None;
    };
    if args.len() != 1 || !is_x_squared(args[0], x) {
        return None;
    }
    // ∫ sin(x²) dx = √(π/2)·S(√(2/π)·x); same prefactor for cos → C.
    let target = match name.as_str() {
        "sin" => "fresnels",
        "cos" => "fresnelc",
        _ => return None,
    };
    let pi = ctx.var("pi");
    let two = ctx.num(2);
    // √(π/2) = sqrt(pi·2⁻¹)
    let prefactor = ctx.fun("sqrt", &[ctx.mul(&[pi, ctx.pow(two, ctx.num(-1))])]);
    // √(2/π)·x = sqrt(2·pi⁻¹)·x
    let inner = ctx.mul(&[
        ctx.fun("sqrt", &[ctx.mul(&[two, ctx.pow(pi, ctx.num(-1))])]),
        x,
    ]);
    Some(ctx.mul(&[prefactor, ctx.fun(target, &[inner])]))
}

// ==================================================================
//  0.27.3 Phase B — generic special-function reduction families
// ==================================================================

/// Maximum polynomial degree expanded by the Phase-B families.
const MAX_SPECIAL_DEG: i64 = 8;
/// Hard step budget for every Phase-B reduction chain.
const MAX_SPECIAL_STEPS: usize = 16;

/// The special-function heads the reduction families know.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Head {
    Erf,
    Erfc,
    Erfi,
    Si,
    Ci,
    Shi,
    Chi,
    /// `Ei(u)`, the one-argument exponential integral (`Ei' = e^u/u`).
    Ei,
}

impl Head {
    fn from_name(name: &str) -> Option<Head> {
        Some(match name {
            "erf" => Head::Erf,
            "erfc" => Head::Erfc,
            "erfi" => Head::Erfi,
            "Si" => Head::Si,
            "Ci" => Head::Ci,
            "Shi" => Head::Shi,
            "Chi" => Head::Chi,
            "Ei" => Head::Ei,
            _ => return None,
        })
    }

    fn name(self) -> &'static str {
        match self {
            Head::Erf => "erf",
            Head::Erfc => "erfc",
            Head::Erfi => "erfi",
            Head::Si => "Si",
            Head::Ci => "Ci",
            Head::Shi => "Shi",
            Head::Chi => "Chi",
            Head::Ei => "Ei",
        }
    }

    /// Whether `∫F(u)² du` is implemented for this head.
    fn square_capable(self) -> bool {
        !matches!(self, Head::Ei | Head::Erfc)
    }
}

/// The elementary companions of a head's derivative.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Sin,
    Cos,
    Sinh,
    Cosh,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Sin => "sin",
            Kind::Cos => "cos",
            Kind::Sinh => "sinh",
            Kind::Cosh => "cosh",
        }
    }

    fn atom<'a>(self, ctx: &'a AtomArena<'a>, t: Atom<'a>) -> Atom<'a> {
        ctx.fun(self.name(), &[t])
    }
}

/// One recognized special-function kernel, extracted from the integrand.
#[derive(Clone, Copy, Debug)]
enum Kernel<'a> {
    /// `F(u)` with `u` linear in the integration variable.
    Plain(Head, Atom<'a>),
    /// `F(u)^2` with `u` linear and `F` square-capable.
    Square(Head, Atom<'a>),
    /// `Ei(n, u)`, the order-`n` exponential integral, `n` a small integer.
    Order(i64, Atom<'a>),
}

/// The linear kernel argument `u = a + b·x`, with the pieces the reduction
/// formulas need.
#[derive(Clone, Copy)]
struct Lin<'a> {
    u: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    x: Atom<'a>,
    /// `a` is literally zero, so `x = u/b` is a monomial.
    a_zero: bool,
}

// ------------------------------------------------------------------
//  Small atom builders
// ------------------------------------------------------------------

fn is_zero(a: Atom<'_>) -> bool {
    matches!(a.node(), AtomNode::Num(0))
}

/// Product, folding the trivial unit/zero factors so emitted forms stay tidy.
fn muln<'a>(ctx: &'a AtomArena<'a>, args: &[Atom<'a>]) -> Atom<'a> {
    let mut out = Vec::with_capacity(args.len());
    for a in args {
        match a.node() {
            AtomNode::Num(0) => return ctx.num(0),
            AtomNode::Num(1) => {}
            _ => out.push(*a),
        }
    }
    match out.len() {
        0 => ctx.num(1),
        1 => out[0],
        _ => ctx.mul(&out),
    }
}

/// Sum, dropping zero terms.
fn addn<'a>(ctx: &'a AtomArena<'a>, args: &[Atom<'a>]) -> Atom<'a> {
    let mut out = Vec::with_capacity(args.len());
    for a in args {
        if !is_zero(*a) {
            out.push(*a);
        }
    }
    match out.len() {
        0 => ctx.num(0),
        1 => out[0],
        _ => ctx.add(&out),
    }
}

fn neg<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    match a.node() {
        AtomNode::Num(v) => ctx.num(-v),
        _ => ctx.mul(&[ctx.num(-1), a]),
    }
}

/// `base^e` for an integer exponent, folding the degenerate cases.
fn powi_clean<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, e: i64) -> Atom<'a> {
    if e == 0 {
        return ctx.num(1);
    }
    if let AtomNode::Num(v) = base.node() {
        match *v {
            1 => return ctx.num(1),
            0 => {
                return if e > 0 {
                    ctx.num(0)
                } else {
                    ctx.pow(base, ctx.num(e))
                };
            }
            -1 => return ctx.num(if e % 2 == 0 { 1 } else { -1 }),
            _ => {
                if e > 0
                    && let Ok(e32) = u32::try_from(e)
                    && let Some(p) = v.checked_pow(e32)
                {
                    return ctx.num(p);
                }
            }
        }
    }
    if e == 1 {
        return base;
    }
    ctx.pow(base, ctx.num(e))
}

/// The atom `p/q` for small integers, in the tidiest canonical shape.
fn q_frac<'a>(ctx: &'a AtomArena<'a>, p: i64, q: i64) -> Atom<'a> {
    debug_assert!(q != 0, "q_frac with zero denominator");
    if p % q == 0 {
        return ctx.num(p / q);
    }
    if p == 1 {
        return powi_clean(ctx, ctx.num(q), -1);
    }
    if p == -1 {
        return muln(ctx, &[ctx.num(-1), powi_clean(ctx, ctx.num(q), -1)]);
    }
    muln(ctx, &[ctx.num(p), powi_clean(ctx, ctx.num(q), -1)])
}

fn exp_atom<'a>(ctx: &'a AtomArena<'a>, arg: Atom<'a>) -> Atom<'a> {
    ctx.fun("exp", &[arg])
}

fn sqrt_pi<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    ctx.fun("sqrt", &[ctx.var("pi")])
}

fn inv_sqrt_pi<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    powi_clean(ctx, sqrt_pi(ctx), -1)
}

/// `C(n, k)` for the small degrees the expansions need.
fn binom_k(n: i64, k: i64) -> i64 {
    if k < 0 || k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut acc = 1i64;
    for i in 0..k {
        acc = acc * (n - i) / (i + 1);
    }
    acc
}

// ------------------------------------------------------------------
//  Kernel recognition
// ------------------------------------------------------------------

/// Whether `u` is a genuine linear form `a + b·x` with `b ≠ 0`.
fn lin_of<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, x: Symbol) -> Option<Lin<'a>> {
    let (b, a) = linear_form(ctx, u, x)?;
    if !is_constant(a, x) || !is_constant(b, x) || is_zero(b) {
        return None;
    }
    let b = normalize(ctx, b);
    let a = normalize(ctx, a);
    Some(Lin {
        u,
        a_zero: is_zero(a),
        a,
        b,
        x: ctx.var(x.as_str()),
    })
}

/// Match a single special head applied to a linear argument.
fn match_plain<'a>(ctx: &'a AtomArena<'a>, f: Atom<'a>, x: Symbol) -> Option<Kernel<'a>> {
    let AtomNode::Fun(name, args) = f.node() else {
        return None;
    };
    if name.as_str() == "Ei" {
        if args.len() == 2 {
            if let AtomNode::Num(n) = args[0].node()
                && (-MAX_SPECIAL_DEG..=MAX_SPECIAL_DEG).contains(n)
                && lin_of(ctx, args[1], x).is_some()
            {
                return Some(Kernel::Order(*n, args[1]));
            }
        }
        if args.len() == 1 && lin_of(ctx, args[0], x).is_some() {
            return Some(Kernel::Plain(Head::Ei, args[0]));
        }
        return None;
    }
    let h = Head::from_name(name.as_str())?;
    if args.len() != 1 {
        return None;
    }
    lin_of(ctx, args[0], x)?;
    Some(Kernel::Plain(h, args[0]))
}

/// Match `F(u)^2` for a square-capable head.
fn match_square<'a>(ctx: &'a AtomArena<'a>, f: Atom<'a>, x: Symbol) -> Option<Kernel<'a>> {
    let AtomNode::Pow(base, exp) = f.node() else {
        return None;
    };
    if !matches!(exp.node(), AtomNode::Num(2)) {
        return None;
    }
    let AtomNode::Fun(name, args) = base.node() else {
        return None;
    };
    let h = Head::from_name(name.as_str()).filter(|h| h.square_capable())?;
    if args.len() != 1 {
        return None;
    }
    lin_of(ctx, args[0], x)?;
    Some(Kernel::Square(h, args[0]))
}

/// Split the integrand factors into the unique special kernel and the rest.
///
/// Two kernels in one product are declined: every reduction in this file is
/// written for a single special factor.
fn split_kernel<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Symbol,
) -> Option<(Kernel<'a>, Vec<Atom<'a>>)> {
    let mut found: Option<(usize, Kernel<'a>)> = None;
    for (i, f) in factors.iter().enumerate() {
        if let Some(k) = match_plain(ctx, *f, x).or_else(|| match_square(ctx, *f, x)) {
            if found.is_some() {
                return None;
            }
            found = Some((i, k));
        }
    }
    let (idx, kernel) = found?;
    let rest = factors
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != idx)
        .map(|(_, a)| *a)
        .collect();
    Some((kernel, rest))
}

/// Entry point of the Phase-B families.
fn reduction_families<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let xs = match x.node() {
        AtomNode::Var(v) => *v,
        _ => return None,
    };
    let (kernel, rest) = split_kernel(ctx, factors, xs)?;
    let rest = muln(ctx, &rest);
    match kernel {
        Kernel::Order(n, u) => order_family(ctx, n, u, rest, xs),
        Kernel::Plain(h, u) => plain_family(ctx, h, u, rest, xs),
        Kernel::Square(h, u) => square_family(ctx, h, u, rest, xs),
    }
}

// ------------------------------------------------------------------
//  Polynomial arithmetic (coefficient lists in x)
// ------------------------------------------------------------------

fn poly_add<'a>(ctx: &'a AtomArena<'a>, a: &[Atom<'a>], b: &[Atom<'a>]) -> Option<Vec<Atom<'a>>> {
    let len = a.len().max(b.len());
    if len as i64 > MAX_SPECIAL_DEG + 1 {
        return None;
    }
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let x = a.get(i).copied().unwrap_or_else(|| ctx.num(0));
        let y = b.get(i).copied().unwrap_or_else(|| ctx.num(0));
        out.push(addn(ctx, &[x, y]));
    }
    Some(out)
}

fn poly_mul<'a>(ctx: &'a AtomArena<'a>, a: &[Atom<'a>], b: &[Atom<'a>]) -> Option<Vec<Atom<'a>>> {
    if a.len() + b.len() - 1 > (MAX_SPECIAL_DEG + 1) as usize {
        return None;
    }
    let mut out = vec![ctx.num(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if is_zero(*x) {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if is_zero(*y) {
                continue;
            }
            out[i + j] = addn(ctx, &[out[i + j], muln(ctx, &[*x, *y])]);
        }
    }
    Some(out)
}

/// Expand `e` as a polynomial in `x` (ascending coefficients), declining
/// anything that is not a polynomial of degree `≤ MAX_SPECIAL_DEG`.
fn poly_coeffs<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>, x: Symbol) -> Option<Vec<Atom<'a>>> {
    if is_constant(e, x) {
        return Some(vec![e]);
    }
    match e.node() {
        AtomNode::Var(v) if *v == x => Some(vec![ctx.num(0), ctx.num(1)]),
        AtomNode::Add(args) => {
            let mut acc = vec![ctx.num(0)];
            for a in args.iter() {
                let c = poly_coeffs(ctx, *a, x)?;
                acc = poly_add(ctx, &acc, &c)?;
            }
            Some(acc)
        }
        AtomNode::Mul(args) => {
            let mut acc = vec![ctx.num(1)];
            for a in args.iter() {
                let c = poly_coeffs(ctx, *a, x)?;
                acc = poly_mul(ctx, &acc, &c)?;
            }
            Some(acc)
        }
        AtomNode::Pow(b, ex) => {
            let AtomNode::Num(n) = ex.node() else {
                return None;
            };
            if *n < 0 || *n > MAX_SPECIAL_DEG {
                return None;
            }
            let base = poly_coeffs(ctx, *b, x)?;
            let mut acc = vec![ctx.num(1)];
            for _ in 0..*n {
                acc = poly_mul(ctx, &acc, &base)?;
            }
            Some(acc)
        }
        _ => None,
    }
}

/// `x^k` for `k ≥ 1`, used by the monomial-denominator matcher.
fn x_power(e: Atom<'_>, x: Symbol) -> Option<i64> {
    match e.node() {
        AtomNode::Var(v) if *v == x => Some(1),
        AtomNode::Pow(b, ex) => {
            let AtomNode::Num(k) = ex.node() else {
                return None;
            };
            if *k >= 1 && x_power(*b, x) == Some(1) {
                Some(*k)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Match `c·x^(−m)` with constant `c` and `m ≥ 1`.
fn neg_power_monomial<'a>(
    ctx: &'a AtomArena<'a>,
    e: Atom<'a>,
    x: Symbol,
) -> Option<(Atom<'a>, i64)> {
    let neg_pow = |f: Atom<'a>| -> Option<i64> {
        let AtomNode::Pow(b, ex) = f.node() else {
            return None;
        };
        let AtomNode::Num(m) = ex.node() else {
            return None;
        };
        if *m >= 0 {
            return None;
        }
        let k = x_power(*b, x)?;
        let eff = k.checked_mul(-*m)?;
        if eff >= 1 { Some(eff) } else { None }
    };
    match e.node() {
        AtomNode::Pow(_, _) => Some((ctx.num(1), neg_pow(e)?)),
        AtomNode::Mul(args) => {
            let mut coeff: Vec<Atom<'a>> = Vec::new();
            let mut power: Option<i64> = None;
            for a in args.iter() {
                if is_constant(*a, x) {
                    coeff.push(*a);
                    continue;
                }
                if power.is_some() {
                    return None;
                }
                power = Some(neg_pow(*a)?);
            }
            Some((muln(ctx, &coeff), power?))
        }
        _ => None,
    }
}

// ------------------------------------------------------------------
//  Elementary building blocks
// ------------------------------------------------------------------

/// `∫ t^m e^(c·t) dt`, elementary by repeated parts.
fn int_pow_exp<'a>(
    ctx: &'a AtomArena<'a>,
    t: Atom<'a>,
    c: Atom<'a>,
    m: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || m < 0 {
        return None;
    }
    let cinv = powi_clean(ctx, c, -1);
    if m == 0 {
        return Some(muln(ctx, &[exp_atom(ctx, muln(ctx, &[c, t])), cinv]));
    }
    let rec = int_pow_exp(ctx, t, c, m - 1, budget - 1)?;
    Some(addn(
        ctx,
        &[
            muln(
                ctx,
                &[
                    powi_clean(ctx, t, m),
                    exp_atom(ctx, muln(ctx, &[c, t])),
                    cinv,
                ],
            ),
            neg(ctx, muln(ctx, &[ctx.num(m), cinv, rec])),
        ],
    ))
}

/// `∫ t^m e^(α·t²) dt` for `α = ±1, ±2`; elementary, ending in
/// `erf(√|α|·t)` / `erfi(√|α|·t)` for even `m`.
fn int_pow_exp_quad<'a>(
    ctx: &'a AtomArena<'a>,
    t: Atom<'a>,
    m: i64,
    alpha: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || m < 0 || alpha == 0 {
        return None;
    }
    let t2 = powi_clean(ctx, t, 2);
    let ex = exp_atom(ctx, muln(ctx, &[ctx.num(alpha), t2]));
    if m == 0 {
        // √π/(2√|α|) · erf/erfi(√|α|·t)
        let root = match alpha.abs() {
            1 => None,
            2 => Some(ctx.fun("sqrt", &[ctx.num(2)])),
            _ => return None,
        };
        let arg = match root {
            Some(r) => muln(ctx, &[r, t]),
            None => t,
        };
        let f = if alpha > 0 {
            ctx.fun("erfi", &[arg])
        } else {
            ctx.fun("erf", &[arg])
        };
        let coef = match root {
            Some(r) => muln(
                ctx,
                &[q_frac(ctx, 1, 2), sqrt_pi(ctx), powi_clean(ctx, r, -1)],
            ),
            None => muln(ctx, &[q_frac(ctx, 1, 2), sqrt_pi(ctx)]),
        };
        return Some(muln(ctx, &[coef, f]));
    }
    if m % 2 == 1 {
        let half = (m - 1) / 2;
        let ie = int_pow_exp(ctx, t2, ctx.num(alpha), half, budget - 1)?;
        return Some(muln(ctx, &[q_frac(ctx, 1, 2), ie]));
    }
    // [t^(m−1)e^(αt²) − (m−1)∫t^(m−2)e^(αt²)dt] / (2α)
    let rec = int_pow_exp_quad(ctx, t, m - 2, alpha, budget - 1)?;
    let inner = addn(
        ctx,
        &[
            muln(ctx, &[powi_clean(ctx, t, m - 1), ex]),
            neg(ctx, muln(ctx, &[ctx.num(m - 1), rec])),
        ],
    );
    Some(muln(ctx, &[q_frac(ctx, 1, 2 * alpha), inner]))
}

/// `∫ t^m sin/cos/sinh/cosh(t) dt`, elementary by repeated parts.
fn int_pow_trig<'a>(
    ctx: &'a AtomArena<'a>,
    t: Atom<'a>,
    m: i64,
    kind: Kind,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || m < 0 {
        return None;
    }
    if m == 0 {
        return Some(match kind {
            Kind::Sin => neg(ctx, Kind::Cos.atom(ctx, t)),
            Kind::Cos => Kind::Sin.atom(ctx, t),
            Kind::Sinh => Kind::Cosh.atom(ctx, t),
            Kind::Cosh => Kind::Sinh.atom(ctx, t),
        });
    }
    let (next, sign) = match kind {
        Kind::Sin => (Kind::Cos, 1),
        Kind::Cos => (Kind::Sin, -1),
        Kind::Sinh => (Kind::Cosh, -1),
        Kind::Cosh => (Kind::Sinh, -1),
    };
    let rec = int_pow_trig(ctx, t, m - 1, next, budget - 1)?;
    let prim = match kind {
        Kind::Sin => neg(ctx, Kind::Cos.atom(ctx, t)),
        Kind::Cos => Kind::Sin.atom(ctx, t),
        Kind::Sinh => Kind::Cosh.atom(ctx, t),
        Kind::Cosh => Kind::Sinh.atom(ctx, t),
    };
    let lead = muln(ctx, &[powi_clean(ctx, t, m), prim]);
    let tail = muln(ctx, &[ctx.num(m), rec]);
    Some(addn(
        ctx,
        &[lead, if sign > 0 { tail } else { neg(ctx, tail) }],
    ))
}

/// `∫ e^(c·u) u^(−j) du` for `j ≥ 1`; terminates in `Ei(c·u)`.
///
/// `∫e^(cu)u^{−j}du = −e^(cu)u^{1−j}/(j−1) + (c/(j−1))∫e^(cu)u^{1−j}du`,
/// and `d/du Ei(cu) = e^(cu)/u` (no `1/c` factor: the chain rule cancels it).
fn int_exp_over_pow<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    c: Atom<'a>,
    j: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || j < 1 {
        return None;
    }
    let ex = exp_atom(ctx, muln(ctx, &[c, u]));
    if j == 1 {
        return Some(ctx.fun("Ei", &[muln(ctx, &[c, u])]));
    }
    let rec = int_exp_over_pow(ctx, u, c, j - 1, budget - 1)?;
    Some(addn(
        ctx,
        &[
            neg(
                ctx,
                muln(ctx, &[ex, powi_clean(ctx, u, 1 - j), q_frac(ctx, 1, j - 1)]),
            ),
            muln(ctx, &[c, q_frac(ctx, 1, j - 1), rec]),
        ],
    ))
}

// ------------------------------------------------------------------
//  Family 1: polynomial × F(a + b·x)
// ------------------------------------------------------------------

fn plain_family<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    rest: Atom<'a>,
    x: Symbol,
) -> Option<Atom<'a>> {
    let lin = lin_of(ctx, u, x)?;
    if let Some(coeffs) = poly_coeffs(ctx, rest, x) {
        return poly_plain(ctx, h, lin, &coeffs);
    }
    // `F(b·x)/x^m`: only `Ei` has an elementary-logarithm descent.
    if h == Head::Ei && lin.a_zero {
        let (c, m) = neg_power_monomial(ctx, rest, x)?;
        let d = d_ei(ctx, lin.u, m)?;
        return Some(muln(ctx, &[c, powi_clean(ctx, lin.b, m - 1), d]));
    }
    None
}

/// `∫ P(x)·F(a + b·x) dx` via `x = (u − a)/b` and the parts recursion.
fn poly_plain<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    lin: Lin<'a>,
    coeffs: &[Atom<'a>],
) -> Option<Atom<'a>> {
    let mut terms = Vec::new();
    for (i, c) in coeffs.iter().enumerate() {
        if is_zero(*c) {
            continue;
        }
        let i = i as i64;
        let mut inner = Vec::new();
        for j in 0..=i {
            let s = s_plain(ctx, h, lin.u, j, MAX_SPECIAL_STEPS)?;
            let w = muln(
                ctx,
                &[
                    ctx.num(binom_k(i, j)),
                    powi_clean(ctx, neg(ctx, lin.a), i - j),
                    s,
                ],
            );
            inner.push(w);
        }
        let scale = muln(ctx, &[*c, powi_clean(ctx, lin.b, -(i + 1))]);
        terms.push(muln(ctx, &[scale, addn(ctx, &inner)]));
    }
    Some(addn(ctx, &terms))
}

/// `∫ u^j F(u) du = u^(j+1)F(u)/(j+1) − (1/(j+1))∫u^(j+1)F'(u)du`.
fn s_plain<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    j: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || j < 0 {
        return None;
    }
    let jp = j + 1;
    let lead = muln(
        ctx,
        &[
            powi_clean(ctx, u, jp),
            ctx.fun(h.name(), &[u]),
            q_frac(ctx, 1, jp),
        ],
    );
    let r = r_plain(ctx, h, u, jp, budget - 1)?;
    Some(addn(
        ctx,
        &[lead, neg(ctx, muln(ctx, &[q_frac(ctx, 1, jp), r]))],
    ))
}

/// `∫ u^k F'(u) du` for the heads with an elementary derivative.
fn r_plain<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    k: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || k < 1 {
        return None;
    }
    let two_over_sqrt_pi = muln(ctx, &[ctx.num(2), inv_sqrt_pi(ctx)]);
    Some(match h {
        Head::Erf => muln(
            ctx,
            &[
                two_over_sqrt_pi,
                int_pow_exp_quad(ctx, u, k, -1, budget - 1)?,
            ],
        ),
        Head::Erfc => neg(
            ctx,
            muln(
                ctx,
                &[
                    two_over_sqrt_pi,
                    int_pow_exp_quad(ctx, u, k, -1, budget - 1)?,
                ],
            ),
        ),
        Head::Erfi => muln(
            ctx,
            &[
                two_over_sqrt_pi,
                int_pow_exp_quad(ctx, u, k, 1, budget - 1)?,
            ],
        ),
        Head::Si => int_pow_trig(ctx, u, k - 1, Kind::Sin, budget - 1)?,
        Head::Ci => int_pow_trig(ctx, u, k - 1, Kind::Cos, budget - 1)?,
        Head::Shi => int_pow_trig(ctx, u, k - 1, Kind::Sinh, budget - 1)?,
        Head::Chi => int_pow_trig(ctx, u, k - 1, Kind::Cosh, budget - 1)?,
        Head::Ei => int_pow_exp(ctx, u, ctx.num(1), k - 1, budget - 1)?,
    })
}

/// `∫ Ei(u) u^(−m) du` for `m ≥ 2` (`Ei' = e^u/u`).
///
/// `∫F(u)u^(−m)du = −F(u)u^(1−m)/(m−1) + (1/(m−1))∫F'(u)u^(1−m)du`.
fn d_ei<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, m: i64) -> Option<Atom<'a>> {
    if m < 2 {
        return None;
    }
    let lead = neg(
        ctx,
        muln(
            ctx,
            &[
                ctx.fun("Ei", &[u]),
                powi_clean(ctx, u, 1 - m),
                q_frac(ctx, 1, m - 1),
            ],
        ),
    );
    let tail = muln(
        ctx,
        &[
            q_frac(ctx, 1, m - 1),
            int_exp_over_pow(ctx, u, ctx.num(1), m, MAX_SPECIAL_STEPS)?,
        ],
    );
    Some(addn(ctx, &[lead, tail]))
}

// ------------------------------------------------------------------
//  Family 2: polynomial × Ei(n, a + b·x)
// ------------------------------------------------------------------

fn order_family<'a>(
    ctx: &'a AtomArena<'a>,
    n: i64,
    u: Atom<'a>,
    rest: Atom<'a>,
    x: Symbol,
) -> Option<Atom<'a>> {
    let lin = lin_of(ctx, u, x)?;
    if let Some(coeffs) = poly_coeffs(ctx, rest, x) {
        // Degrees ≥ 1 need `a = 0` (otherwise `x/u` is a genuine rational
        // factor, not the constant `1/b`).
        if coeffs.len() == 1 || lin.a_zero {
            let mut terms = Vec::new();
            for (k, c) in coeffs.iter().enumerate() {
                if is_zero(*c) {
                    continue;
                }
                let j = j_ei(ctx, &lin, k as i64, n, MAX_SPECIAL_STEPS)?;
                terms.push(muln(ctx, &[*c, j]));
            }
            return Some(addn(ctx, &terms));
        }
        return None;
    }
    if lin.a_zero {
        let (c, m) = neg_power_monomial(ctx, rest, x)?;
        let d = d_order(ctx, lin.u, m, n, MAX_SPECIAL_STEPS)?;
        return Some(muln(ctx, &[c, powi_clean(ctx, lin.b, m - 1), d]));
    }
    None
}

/// `∫ x^k E_n(a + b·x) dx`.
///
/// `k = 0`: `−E_(n+1)(u)/b` for every order. For `a = 0` the inverse
/// recurrence `E_n(u) = (e^−u − n·E_(n+1)(u))/u` gives
/// `∫x^kE_n dx = (1/b)∫x^(k−1)e^−u dx − (n/b)∫x^(k−1)E_(n+1)dx`, which
/// lowers the degree by one per step and ends in elementary moments.
fn j_ei<'a>(
    ctx: &'a AtomArena<'a>,
    lin: &Lin<'a>,
    k: i64,
    n: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || k < 0 {
        return None;
    }
    if k == 0 {
        let head = ctx.fun("Ei", &[ctx.num(n + 1), lin.u]);
        return Some(neg(ctx, muln(ctx, &[head, powi_clean(ctx, lin.b, -1)])));
    }
    if !lin.a_zero {
        return None;
    }
    let cc = neg(ctx, lin.b);
    let binv = powi_clean(ctx, lin.b, -1);
    let base = int_pow_exp(ctx, lin.x, cc, k - 1, MAX_SPECIAL_STEPS)?;
    if n == 0 {
        return Some(muln(ctx, &[binv, base]));
    }
    let rec = j_ei(ctx, lin, k - 1, n + 1, budget - 1)?;
    Some(addn(
        ctx,
        &[
            muln(ctx, &[binv, base]),
            neg(ctx, muln(ctx, &[ctx.num(n), binv, rec])),
        ],
    ))
}

/// `∫ E_n(u) u^(−m) du`.
///
/// `n > 0`: the parts descent `∫E_n u^−m = −E_n u^(1−m)/(m−1) − (1/(m−1))∫E_(n−1)u^(1−m)`
/// (using `E_n' = −E_(n−1)`), which needs `m > n` so the chain ends at
/// `E_0(u) = e^−u/u`.
///
/// `n < 0`: the recurrence form `E_n = (e^−u − n·E_(n+1))/u` gives
/// `∫E_n u^−m = ∫e^−u u^(−m−1) − n∫E_(n+1)u^(−m−1)`, which raises the order
/// toward zero and terminates at `n = 0`.
fn d_order<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    m: i64,
    n: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || m < 1 {
        return None;
    }
    if n == 0 {
        return int_exp_over_pow(ctx, u, ctx.num(-1), m + 1, MAX_SPECIAL_STEPS);
    }
    if n > 0 {
        if m <= n {
            return None;
        }
        let lead = neg(
            ctx,
            muln(
                ctx,
                &[
                    ctx.fun("Ei", &[ctx.num(n), u]),
                    powi_clean(ctx, u, 1 - m),
                    q_frac(ctx, 1, m - 1),
                ],
            ),
        );
        let rec = d_order(ctx, u, m - 1, n - 1, budget - 1)?;
        return Some(addn(
            ctx,
            &[lead, neg(ctx, muln(ctx, &[q_frac(ctx, 1, m - 1), rec]))],
        ));
    }
    let lead = int_exp_over_pow(ctx, u, ctx.num(-1), m + 1, MAX_SPECIAL_STEPS)?;
    let rec = d_order(ctx, u, m + 1, n + 1, budget - 1)?;
    Some(addn(ctx, &[lead, neg(ctx, muln(ctx, &[ctx.num(n), rec]))]))
}

// ------------------------------------------------------------------
//  Family 3: polynomial × F(a + b·x)²
// ------------------------------------------------------------------

fn square_family<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    rest: Atom<'a>,
    x: Symbol,
) -> Option<Atom<'a>> {
    let lin = lin_of(ctx, u, x)?;
    let coeffs = poly_coeffs(ctx, rest, x)?;
    let mut terms = Vec::new();
    for (i, c) in coeffs.iter().enumerate() {
        if is_zero(*c) {
            continue;
        }
        let i = i as i64;
        let mut inner = Vec::new();
        for j in 0..=i {
            let t = t_square(ctx, h, lin.u, j, MAX_SPECIAL_STEPS)?;
            let w = muln(
                ctx,
                &[
                    ctx.num(binom_k(i, j)),
                    powi_clean(ctx, neg(ctx, lin.a), i - j),
                    t,
                ],
            );
            inner.push(w);
        }
        let scale = muln(ctx, &[*c, powi_clean(ctx, lin.b, -(i + 1))]);
        terms.push(muln(ctx, &[scale, addn(ctx, &inner)]));
    }
    Some(addn(ctx, &terms))
}

/// `∫ u^j F(u)² du = u^(j+1)F²/(j+1) − (2/(j+1))∫u^(j+1)F·F' du`.
fn t_square<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    j: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || j < 0 {
        return None;
    }
    let jp = j + 1;
    let sq = ctx.pow(ctx.fun(h.name(), &[u]), ctx.num(2));
    let lead = muln(ctx, &[powi_clean(ctx, u, jp), sq, q_frac(ctx, 1, jp)]);
    let us = u_square(ctx, h, u, jp, budget - 1)?;
    Some(addn(
        ctx,
        &[
            lead,
            neg(ctx, muln(ctx, &[ctx.num(2), q_frac(ctx, 1, jp), us])),
        ],
    ))
}

/// `∫ u^m F(u)F'(u) du` for `m ≥ 1`.
fn u_square<'a>(
    ctx: &'a AtomArena<'a>,
    h: Head,
    u: Atom<'a>,
    m: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || m < 1 {
        return None;
    }
    let two_over_sqrt_pi = muln(ctx, &[ctx.num(2), inv_sqrt_pi(ctx)]);
    Some(match h {
        Head::Erf => muln(
            ctx,
            &[two_over_sqrt_pi, m_gauss(ctx, u, -1, m, budget - 1)?],
        ),
        Head::Erfi => muln(ctx, &[two_over_sqrt_pi, m_gauss(ctx, u, 1, m, budget - 1)?]),
        Head::Ci => ab_pair(ctx, u, false, true, m - 1, budget - 1)?.0,
        Head::Si => ab_pair(ctx, u, false, false, m - 1, budget - 1)?.1,
        Head::Chi => ab_pair(ctx, u, true, true, m - 1, budget - 1)?.0,
        Head::Shi => ab_pair(ctx, u, true, false, m - 1, budget - 1)?.1,
        _ => return None,
    })
}

/// `M(σ, j) = ∫u^j F_σ(u)e^(σu²)du` with `F_−1 = erf`, `F_+1 = erfi`.
///
/// From `d/du[F_σ e^(σu²)] = (2/√π)e^(2σu²) + 2σu F_σ e^(σu²)`:
/// `M(j) = [u^(j−1)F_σ e^(σu²) − (j−1)M(j−2) − (2/√π)G(j−1)]/(2σ)`,
/// where `G(k) = ∫u^k e^(2σu²)du` is elementary.
fn m_gauss<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    sigma: i64,
    j: i64,
    budget: usize,
) -> Option<Atom<'a>> {
    if budget == 0 || j < 0 || (sigma != 1 && sigma != -1) {
        return None;
    }
    let f = if sigma > 0 {
        ctx.fun("erfi", &[u])
    } else {
        ctx.fun("erf", &[u])
    };
    let ex = exp_atom(ctx, muln(ctx, &[ctx.num(sigma), powi_clean(ctx, u, 2)]));
    if j == 0 {
        // (√π/4)·F_σ²
        return Some(muln(
            ctx,
            &[q_frac(ctx, 1, 4), sqrt_pi(ctx), ctx.pow(f, ctx.num(2))],
        ));
    }
    if j == 1 {
        // F e^(σu²)/(2σ) − σ·G(0)/√π
        let t1 = muln(ctx, &[f, ex, q_frac(ctx, 1, 2 * sigma)]);
        let g = int_pow_exp_quad(ctx, u, 0, 2 * sigma, budget - 1)?;
        let t2 = muln(ctx, &[ctx.num(sigma), inv_sqrt_pi(ctx), g]);
        return Some(addn(ctx, &[t1, neg(ctx, t2)]));
    }
    let lead = muln(ctx, &[powi_clean(ctx, u, j - 1), f, ex]);
    let rec = m_gauss(ctx, u, sigma, j - 2, budget - 1)?;
    let g = int_pow_exp_quad(ctx, u, j - 1, 2 * sigma, budget - 1)?;
    let two_over_sqrt_pi = muln(ctx, &[ctx.num(2), inv_sqrt_pi(ctx)]);
    let inner = addn(
        ctx,
        &[
            lead,
            neg(ctx, muln(ctx, &[ctx.num(j - 1), rec])),
            neg(ctx, muln(ctx, &[two_over_sqrt_pi, g])),
        ],
    );
    Some(muln(ctx, &[q_frac(ctx, 1, 2 * sigma), inner]))
}

/// `W_s(p) = ∫u^p s(2u) du` and `W_c(p) = ∫u^p C(2u) du` for the
/// trigonometric (`s = sin`, `C = cos`) or hyperbolic pair.
///
/// `W_s(p) = u^p·σC(2u)/2 − (pσ/2)W_c(p−1)`,
/// `W_c(p) = u^p·s(2u)/2 − (p/2)W_s(p−1)` with `σ = −1` (trig) or `+1`
/// (hyperbolic).
fn w_pair<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, hyp: bool, p: i64) -> (Atom<'a>, Atom<'a>) {
    let sigma = if hyp { 1 } else { -1 };
    let two_u = muln(ctx, &[ctx.num(2), u]);
    let s2 = if hyp {
        Kind::Sinh.atom(ctx, two_u)
    } else {
        Kind::Sin.atom(ctx, two_u)
    };
    let c2 = if hyp {
        Kind::Cosh.atom(ctx, two_u)
    } else {
        Kind::Cos.atom(ctx, two_u)
    };
    let mut ws = muln(ctx, &[ctx.num(sigma), q_frac(ctx, 1, 2), c2]);
    let mut wc = muln(ctx, &[q_frac(ctx, 1, 2), s2]);
    for q in 1..=p.max(0) {
        let nws = addn(
            ctx,
            &[
                muln(
                    ctx,
                    &[powi_clean(ctx, u, q), ctx.num(sigma), q_frac(ctx, 1, 2), c2],
                ),
                neg(
                    ctx,
                    muln(ctx, &[ctx.num(q), ctx.num(sigma), q_frac(ctx, 1, 2), wc]),
                ),
            ],
        );
        let nwc = addn(
            ctx,
            &[
                muln(ctx, &[powi_clean(ctx, u, q), q_frac(ctx, 1, 2), s2]),
                neg(ctx, muln(ctx, &[ctx.num(q), q_frac(ctx, 1, 2), ws])),
            ],
        );
        ws = nws;
        wc = nwc;
    }
    (ws, wc)
}

/// `A(j) = ∫u^j F(u)c(u) du`, `B(j) = ∫u^j F(u)s(u) du` for the integral
/// heads `Ci`/`Si` (trig) and `Chi`/`Shi` (hyperbolic), `cos_type` marking
/// the heads with `F' = c/u`.
///
/// Parts gives the coupled recursion
/// `A(j) = u^j F s − jB(j−1) − ∫u^jF's du`,
/// `B(j) = σ[u^j F c − jA(j−1) − ∫u^jF'c du]`,
/// whose residual integrals are elementary (listed below), so the pair
/// terminates at the verified `j = 0` bases.
fn ab_pair<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    hyp: bool,
    cos_type: bool,
    j: i64,
    budget: usize,
) -> Option<(Atom<'a>, Atom<'a>)> {
    if budget == 0 || j < 0 {
        return None;
    }
    let sigma = if hyp { 1 } else { -1 };
    let (s_kind, c_kind) = if hyp {
        (Kind::Sinh, Kind::Cosh)
    } else {
        (Kind::Sin, Kind::Cos)
    };
    let two_u = muln(ctx, &[ctx.num(2), u]);
    let f_name = match (hyp, cos_type) {
        (false, true) => "Ci",
        (false, false) => "Si",
        (true, true) => "Chi",
        (true, false) => "Shi",
    };
    let f = ctx.fun(f_name, &[u]);
    if j == 0 {
        // Bases, all verified by differentiation in the unit tests.
        let a0 = match (hyp, cos_type) {
            // Ci
            (false, true) => addn(
                ctx,
                &[
                    muln(ctx, &[f, s_kind.atom(ctx, u)]),
                    neg(
                        ctx,
                        muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Si", &[two_u])]),
                    ),
                ],
            ),
            // Si
            (false, false) => addn(
                ctx,
                &[
                    muln(ctx, &[f, s_kind.atom(ctx, u)]),
                    neg(ctx, muln(ctx, &[q_frac(ctx, 1, 2), log_u(ctx, u)])),
                    muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Ci", &[two_u])]),
                ],
            ),
            // Chi
            (true, true) => addn(
                ctx,
                &[
                    muln(ctx, &[f, s_kind.atom(ctx, u)]),
                    neg(
                        ctx,
                        muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Shi", &[two_u])]),
                    ),
                ],
            ),
            // Shi
            (true, false) => addn(
                ctx,
                &[
                    muln(ctx, &[f, s_kind.atom(ctx, u)]),
                    neg(
                        ctx,
                        muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Chi", &[two_u])]),
                    ),
                    muln(ctx, &[q_frac(ctx, 1, 2), log_u(ctx, u)]),
                ],
            ),
        };
        let b0 = match (hyp, cos_type) {
            // Ci
            (false, true) => addn(
                ctx,
                &[
                    neg(ctx, muln(ctx, &[f, c_kind.atom(ctx, u)])),
                    muln(ctx, &[q_frac(ctx, 1, 2), log_u(ctx, u)]),
                    muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Ci", &[two_u])]),
                ],
            ),
            // Si
            (false, false) => addn(
                ctx,
                &[
                    neg(ctx, muln(ctx, &[f, c_kind.atom(ctx, u)])),
                    muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Si", &[two_u])]),
                ],
            ),
            // Chi
            (true, true) => addn(
                ctx,
                &[
                    muln(ctx, &[f, c_kind.atom(ctx, u)]),
                    neg(ctx, muln(ctx, &[q_frac(ctx, 1, 2), log_u(ctx, u)])),
                    neg(
                        ctx,
                        muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Chi", &[two_u])]),
                    ),
                ],
            ),
            // Shi
            (true, false) => addn(
                ctx,
                &[
                    muln(ctx, &[f, c_kind.atom(ctx, u)]),
                    neg(
                        ctx,
                        muln(ctx, &[q_frac(ctx, 1, 2), ctx.fun("Shi", &[two_u])]),
                    ),
                ],
            ),
        };
        return Some((a0, b0));
    }

    let (ap, bp) = ab_pair(ctx, u, hyp, cos_type, j - 1, budget - 1)?;
    let (ws, wc) = w_pair(ctx, u, hyp, j - 1);
    let uj = powi_clean(ctx, u, j);
    let uj_over_j = muln(ctx, &[uj, q_frac(ctx, 1, j)]);
    // ∫u^(j−1) s(2u) du / 2  (the `sc` cross term for both types),
    let sc = muln(ctx, &[q_frac(ctx, 1, 2), ws]);
    // ∫u^(j−1)c² du = ½(u^j/j + W_c) for cos² and cosh² alike,
    let c_sq = muln(ctx, &[q_frac(ctx, 1, 2), addn(ctx, &[uj_over_j, wc])]);
    // ∫u^(j−1)s² du: (1 − cos2u)/2 for sine, but (cosh2u − 1)/2 for sinh.
    let s_sq = if hyp {
        muln(
            ctx,
            &[q_frac(ctx, 1, 2), addn(ctx, &[wc, neg(ctx, uj_over_j)])],
        )
    } else {
        muln(
            ctx,
            &[q_frac(ctx, 1, 2), addn(ctx, &[uj_over_j, neg(ctx, wc)])],
        )
    };
    // ∫u^j F' s du and ∫u^j F' c du, by head type.
    let (ta, tb) = if cos_type { (sc, c_sq) } else { (s_sq, sc) };
    let a = addn(
        ctx,
        &[
            muln(ctx, &[uj, f, s_kind.atom(ctx, u)]),
            neg(ctx, muln(ctx, &[ctx.num(j), bp])),
            neg(ctx, ta),
        ],
    );
    let b_inner = addn(
        ctx,
        &[
            muln(ctx, &[uj, f, c_kind.atom(ctx, u)]),
            neg(ctx, muln(ctx, &[ctx.num(j), ap])),
            neg(ctx, tb),
        ],
    );
    let b = muln(ctx, &[ctx.num(sigma), b_inner]);
    Some((a, b))
}

/// `log(u)` for the logarithm terms of the integral-function bases.
fn log_u<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>) -> Atom<'a> {
    ctx.fun("log", &[u])
}

// ------------------------------------------------------------------
//  Tests
// ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use ocas_core::arena::Arena;

    use super::*;

    fn int<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Option<Atom<'a>> {
        special_integrate(ctx, expr, ctx.var("x"))
    }

    // ---------------- existing table entries ----------------

    #[test]
    fn exp_neg_x_squared_gives_erf() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let r = int(&ctx, expr).expect("erf form");
        assert!(r.to_string().contains("erf"), "got {r}");
    }

    #[test]
    fn exp_x_over_x_gives_ei() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("exp", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Ei form");
        assert_eq!(r.to_string(), "Ei(x)");
    }

    #[test]
    fn sin_x_over_x_gives_si() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("sin", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Si form");
        assert_eq!(r.to_string(), "Si(x)");
    }

    #[test]
    fn cos_x_over_x_gives_ci() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("cos", &[x]), ctx.pow(x, ctx.num(-1))]);
        let r = int(&ctx, expr).expect("Ci form");
        assert_eq!(r.to_string(), "Ci(x)");
    }

    #[test]
    fn sin_x_squared_gives_fresnel() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[ctx.pow(x, ctx.num(2))]);
        let r = int(&ctx, expr).expect("Fresnel form");
        assert!(r.to_string().contains("fresnels"), "got {r}");
    }

    #[test]
    fn unmatched_returns_none() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x) itself is elementary (handled by Risch).
        assert!(int(&ctx, ctx.fun("exp", &[x])).is_none());
        let _ = Symbol::new("x");
    }

    // ---------------- numeric oracle (test-only) ----------------

    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

    /// Composite Simpson rule, enough for the defining integrals below.
    fn quad(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
        let n = if n % 2 == 1 { n + 1 } else { n };
        let h = (b - a) / n as f64;
        let mut s = f(a) + f(b);
        for i in 1..n {
            s += if i % 2 == 1 { 4.0 } else { 2.0 } * f(a + i as f64 * h);
        }
        s * h / 3.0
    }

    fn sin_over(t: f64) -> f64 {
        if t == 0.0 { 1.0 } else { t.sin() / t }
    }
    fn sinh_over(t: f64) -> f64 {
        if t == 0.0 { 1.0 } else { t.sinh() / t }
    }
    fn expm1_over(t: f64) -> f64 {
        if t == 0.0 { 1.0 } else { (t.exp() - 1.0) / t }
    }
    fn cosm1_over(t: f64) -> f64 {
        if t == 0.0 { 0.0 } else { (t.cos() - 1.0) / t }
    }
    fn coshm1_over(t: f64) -> f64 {
        if t == 0.0 { 0.0 } else { (t.cosh() - 1.0) / t }
    }

    fn erf_o(v: f64) -> f64 {
        2.0 / std::f64::consts::PI.sqrt() * quad(&|t| (-t * t).exp(), 0.0, v, 800)
    }
    fn erfi_o(v: f64) -> f64 {
        2.0 / std::f64::consts::PI.sqrt() * quad(&|t| (t * t).exp(), 0.0, v, 800)
    }
    fn si_o(v: f64) -> f64 {
        quad(&sin_over, 0.0, v, 800)
    }
    fn ci_o(v: f64) -> f64 {
        EULER_GAMMA + v.abs().ln() + quad(&cosm1_over, 0.0, v.abs(), 800)
    }
    fn shi_o(v: f64) -> f64 {
        quad(&sinh_over, 0.0, v, 800)
    }
    fn chi_o(v: f64) -> f64 {
        EULER_GAMMA + v.abs().ln() + quad(&coshm1_over, 0.0, v.abs(), 800)
    }
    fn ei_o(v: f64) -> f64 {
        EULER_GAMMA + v.abs().ln() + quad(&expm1_over, 0.0, v, 800)
    }
    fn factorial_o(n: u32) -> f64 {
        (1..=n).map(f64::from).product()
    }

    /// `Ei(n, z)` = `E_n(z)` (Rubi/Mathematica convention).
    fn ei_order_o(n: i64, z: f64) -> Option<f64> {
        if n == 0 {
            return if z == 0.0 { None } else { Some((-z).exp() / z) };
        }
        if n < 0 {
            let m = (-n) as u32;
            let mut sum = 0.0;
            for k in 0..=m {
                sum += z.powi(k as i32 - m as i32 - 1) / factorial_o(k);
            }
            return Some((-z).exp() * factorial_o(m) * sum);
        }
        if z <= 0.0 {
            return None;
        }
        // E_1(z) = −Ei(−z); upward recurrence for higher orders.
        let mut e = -ei_o(-z);
        for k in 1..n {
            e = ((-z).exp() - z * e) / k as f64;
        }
        Some(e)
    }

    fn eval(e: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
        match e.node() {
            AtomNode::Num(n) => Some(*n as f64),
            AtomNode::Var(v) => {
                if v.as_str() == "pi" {
                    return Some(std::f64::consts::PI);
                }
                env.iter().find(|(s, _)| s == v).map(|(_, val)| *val)
            }
            AtomNode::Add(args) => {
                let mut acc = 0.0;
                for a in args.iter() {
                    acc += eval(*a, env)?;
                }
                Some(acc)
            }
            AtomNode::Mul(args) => {
                let mut acc = 1.0;
                for a in args.iter() {
                    acc *= eval(*a, env)?;
                }
                Some(acc)
            }
            AtomNode::Pow(b, ex) => {
                let base = eval(*b, env)?;
                let e = eval(*ex, env)?;
                if base == 0.0 && e < 0.0 {
                    return None;
                }
                let v = base.powf(e);
                if v.is_finite() { Some(v) } else { None }
            }
            AtomNode::Fun(name, args) => {
                if name.as_str() == "Ei" && args.len() == 2 {
                    let n = eval(args[0], env)?;
                    let z = eval(args[1], env)?;
                    return ei_order_o(n as i64, z);
                }
                let v = eval(*args.first()?, env)?;
                let out = match name.as_str() {
                    "sin" => v.sin(),
                    "cos" => v.cos(),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "exp" => v.exp(),
                    "log" => v.abs().ln(),
                    "sqrt" => {
                        if v < 0.0 {
                            return None;
                        }
                        v.sqrt()
                    }
                    "erf" => erf_o(v),
                    "erfi" => erfi_o(v),
                    "Si" => si_o(v),
                    "Ci" => {
                        if v == 0.0 {
                            return None;
                        }
                        ci_o(v)
                    }
                    "Shi" => {
                        if v == 0.0 {
                            return None;
                        }
                        shi_o(v)
                    }
                    "Chi" => {
                        if v == 0.0 {
                            return None;
                        }
                        chi_o(v)
                    }
                    "Ei" if args.len() == 1 => {
                        if v == 0.0 {
                            return None;
                        }
                        ei_o(v)
                    }
                    _ => return None,
                };
                if out.is_finite() { Some(out) } else { None }
            }
        }
    }

    /// Symbolic differentiation of the emitted forms (independent of the
    /// reduction algebra, but sharing the head derivative identities).
    fn diff_local<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>, var: Symbol) -> Atom<'a> {
        match e.node() {
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
            AtomNode::Pow(b, ex) => {
                let (base, exp) = (*b, *ex);
                let db = diff_local(ctx, base, var);
                let de = diff_local(ctx, exp, var);
                let t1 = ctx.mul(&[exp, ctx.pow(base, ctx.add(&[exp, ctx.num(-1)])), db]);
                let t2 = ctx.mul(&[ctx.pow(base, exp), ctx.fun("log", &[base]), de]);
                ctx.add(&[t1, t2])
            }
            AtomNode::Fun(name, args) => {
                // `Ei(n, u)` carries the order first and the argument second.
                let a0 = if name.as_str() == "Ei" && args.len() == 2 {
                    args[1]
                } else {
                    *args.first().unwrap()
                };
                let d0 = diff_local(ctx, a0, var);
                let inner = match name.as_str() {
                    "sin" => ctx.fun("cos", &[a0]),
                    "cos" => ctx.mul(&[ctx.num(-1), ctx.fun("sin", &[a0])]),
                    "sinh" => ctx.fun("cosh", &[a0]),
                    "cosh" => ctx.fun("sinh", &[a0]),
                    "exp" => ctx.fun("exp", &[a0]),
                    "log" => ctx.pow(a0, ctx.num(-1)),
                    "sqrt" => ctx.pow(ctx.mul(&[ctx.num(2), ctx.fun("sqrt", &[a0])]), ctx.num(-1)),
                    "erf" => ctx.mul(&[
                        ctx.num(2),
                        inv_sqrt_pi(ctx),
                        ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(a0, ctx.num(2))])]),
                    ]),
                    "erfi" => ctx.mul(&[
                        ctx.num(2),
                        inv_sqrt_pi(ctx),
                        ctx.fun("exp", &[ctx.pow(a0, ctx.num(2))]),
                    ]),
                    "Si" => ctx.mul(&[ctx.fun("sin", &[a0]), ctx.pow(a0, ctx.num(-1))]),
                    "Ci" => ctx.mul(&[ctx.fun("cos", &[a0]), ctx.pow(a0, ctx.num(-1))]),
                    "Shi" => ctx.mul(&[ctx.fun("sinh", &[a0]), ctx.pow(a0, ctx.num(-1))]),
                    "Chi" => ctx.mul(&[ctx.fun("cosh", &[a0]), ctx.pow(a0, ctx.num(-1))]),
                    "Ei" if args.len() == 1 => {
                        ctx.mul(&[ctx.fun("exp", &[a0]), ctx.pow(a0, ctx.num(-1))])
                    }
                    "Ei" if args.len() == 2 => {
                        let AtomNode::Num(n) = args[0].node() else {
                            panic!("non-numeric Ei order in emitted form");
                        };
                        if *n == 0 {
                            ctx.mul(&[
                                ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), a0])]),
                                ctx.pow(a0, ctx.num(-1)),
                            ])
                        } else {
                            ctx.mul(&[ctx.num(-1), ctx.fun("Ei", &[ctx.num(n - 1), a0])])
                        }
                    }
                    _ => panic!("diff_local: unsupported head {}", name.as_str()),
                };
                ctx.mul(&[inner, d0])
            }
        }
    }

    /// Assert that `special_integrate` solves `s` and that its output really
    /// differentiates back to the integrand at the given samples.
    fn check(s: &str, consts: &[(&str, f64)], samples: &[f64]) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = ocas_parse::parse(&ctx, s).unwrap_or_else(|_| panic!("parse: {s}"));
        let x = ctx.var("x");
        let r = special_integrate(&ctx, e, x).unwrap_or_else(|| panic!("declined: {s}"));
        assert!(!r.to_string().contains("Integral("), "residue for {s}: {r}");
        let d = diff_local(&ctx, r, Symbol::new("x"));
        let base: Vec<(Symbol, f64)> = consts.iter().map(|(n, v)| (Symbol::new(n), *v)).collect();
        let mut checked = 0usize;
        let mut worst = 0.0f64;
        for &xv in samples {
            let mut env = base.clone();
            env.push((Symbol::new("x"), xv));
            let lhs = match eval(d, &env) {
                Some(v) => v,
                None => continue,
            };
            let rhs = eval(e, &env).expect("integrand evaluated");
            let tol = 1e-5 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{s} at x={xv}: d/dx F = {lhs}, integrand = {rhs}\n  F = {r}"
            );
            checked += 1;
            worst = worst.max((lhs - rhs).abs() / rhs.abs().max(1.0));
        }
        assert!(checked >= 2, "only {checked} usable samples for {s} -> {r}");
        println!("  checked={checked} worst_rel={worst:e}  {s}  ->  {r}");
    }

    fn declines(s: &str) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = ocas_parse::parse(&ctx, s).unwrap_or_else(|_| panic!("parse: {s}"));
        assert!(
            special_integrate(&ctx, e, ctx.var("x")).is_none(),
            "expected decline for {s}"
        );
    }

    // ---------------- polynomial × F(u) ----------------

    #[test]
    fn plain_family_polynomial_cases() {
        // 01052: Ci(b*x)
        check("Ci(b*x)", &[("b", 1.5)], &[0.6, 1.1, 1.9]);
        // 01112: x*Chi(b*x)
        check("x*Chi(b*x)", &[("b", 1.3)], &[0.5, 1.0, 1.7]);
        // 00587: x^3*Shi(a + b*x)
        check(
            "x^3*Shi(a + b*x)",
            &[("a", 0.4), ("b", 1.2)],
            &[0.4, 0.9, 1.6],
        );
        // 01763: x^4*erf(b*x)
        check("x^4*erf(b*x)", &[("b", 1.1)], &[0.5, 1.0, 1.6]);
        // Extra heads/shapes covered by the same recursion.
        check("x^2*Si(b*x)", &[("b", 1.4)], &[0.5, 1.1, 1.8]);
        check("x*erfi(b*x)", &[("b", 0.9)], &[0.5, 1.0, 1.5]);
        check("x^2*Ei(b*x)", &[("b", 1.2)], &[0.6, 1.1, 1.7]);
        check(
            "x^2*Shi(a + b*x)",
            &[("a", 0.3), ("b", 1.1)],
            &[0.5, 1.0, 1.5],
        );
    }

    // ---------------- polynomial × Ei(n, u) ----------------

    #[test]
    fn order_family_cases() {
        // 00276: x^4*Ei(-2, b*x)
        check("x^4*Ei(-2, b*x)", &[("b", 1.2)], &[0.6, 1.1, 1.7]);
        // 00953: x^3*Ei(-1, b*x)
        check("x^3*Ei(-1, b*x)", &[("b", 1.3)], &[0.6, 1.1, 1.7]);
        // 01172: Ei(1, a + b*x)
        check(
            "Ei(1, a + b*x)",
            &[("a", 0.5), ("b", 1.2)],
            &[0.4, 0.9, 1.6],
        );
        // Higher positive/negative orders and mixed offsets.
        check("Ei(2, b*x)", &[("b", 1.2)], &[0.6, 1.1, 1.7]);
        check("x^2*Ei(3, b*x)", &[("b", 1.1)], &[0.6, 1.1, 1.7]);
        check(
            "Ei(-2, a + b*x)",
            &[("a", 0.4), ("b", 1.2)],
            &[0.4, 0.9, 1.5],
        );
    }

    // ---------------- F(b*x)/x^m ----------------

    #[test]
    fn power_descent_cases() {
        // 01009: Ei(b*x)/x^4
        check("Ei(b*x)/x^4", &[("b", 1.2)], &[0.7, 1.2, 1.9]);
        // 01267: Ei(-3, b*x)/x^3
        check("Ei(-3, b*x)/x^3", &[("b", 1.2)], &[0.7, 1.2, 1.9]);
        // 01633: Ei(1, b*x)/x^3
        check("Ei(1, b*x)/x^3", &[("b", 1.2)], &[0.7, 1.2, 1.9]);
        // Deeper/shallower denominators and other orders.
        check("Ei(b*x)/x^2", &[("b", 1.3)], &[0.7, 1.2, 1.9]);
        check("Ei(2, b*x)/x^3", &[("b", 1.1)], &[0.7, 1.2, 1.9]);
        check("Ei(-1, b*x)/x^4", &[("b", 1.1)], &[0.7, 1.2, 1.9]);
        check("3*Ei(b*x)/x^3", &[("b", 1.1)], &[0.7, 1.2, 1.9]);
    }

    // ---------------- polynomial × F(u)² ----------------

    #[test]
    fn square_family_cases() {
        // 00123: erfi(b*x)^2
        check("erfi(b*x)^2", &[("b", 1.2)], &[0.5, 1.0, 1.5]);
        // 00575: x^2*Ci(b*x)^2
        check("x^2*Ci(b*x)^2", &[("b", 1.2)], &[0.7, 1.2, 1.8]);
        // 01516: Chi(a + b*x)^2
        check(
            "Chi(a + b*x)^2",
            &[("a", 0.5), ("b", 1.2)],
            &[0.4, 0.9, 1.5],
        );
        // 01729: (c + d*x)*erf(a + b*x)^2
        check(
            "(c + d*x)*erf(a + b*x)^2",
            &[("a", 0.4), ("b", 1.2), ("c", 0.7), ("d", 1.3)],
            &[0.4, 0.9, 1.5],
        );
        // 01799: x*Si(a + b*x)^2
        check(
            "x*Si(a + b*x)^2",
            &[("a", 0.5), ("b", 1.2)],
            &[0.4, 0.9, 1.5],
        );
        // Extra square shapes over the same recursions.
        check("erf(b*x)^2", &[("b", 1.1)], &[0.5, 1.0, 1.5]);
        check("x*Ci(b*x)^2", &[("b", 1.2)], &[0.7, 1.2, 1.8]);
        check(
            "Shi(a + b*x)^2",
            &[("a", 0.5), ("b", 1.1)],
            &[0.4, 0.9, 1.5],
        );
        check("x^3*Shi(b*x)^2", &[("b", 1.1)], &[0.6, 1.1, 1.6]);
        // Non-zero offsets with degree ≥ 2 exercise the full binomial
        // expansion of `(u − a)^k`.
        check(
            "x^2*Ci(a + b*x)^2",
            &[("a", 0.4), ("b", 1.2)],
            &[0.5, 1.0, 1.5],
        );
        check(
            "x^3*Chi(a + b*x)^2",
            &[("a", 0.3), ("b", 1.1)],
            &[0.5, 1.0, 1.5],
        );
        check(
            "x^3*erf(a + b*x)^2",
            &[("a", 0.4), ("b", 1.2)],
            &[0.4, 0.9, 1.5],
        );
        check(
            "x^2*Si(a + b*x)^2",
            &[("a", 0.4), ("b", 1.2)],
            &[0.4, 0.9, 1.5],
        );
        check("x^2*erfi(b*x)^2", &[("b", 1.1)], &[0.5, 1.0, 1.5]);
    }

    // ---------------- declines ----------------

    #[test]
    fn unscoped_shapes_decline() {
        // Rubi itself returns Unintegrable for these two.
        declines("erf(a + b*x)^2/(c + d*x)");
        declines("Ei(n, a + b*x)/(c + d*x)^2");
        // Not implemented (Fresnel denominator chains).
        declines("fresnels(b*x)/x^10");
        declines("fresnels(b*x)/x^8");
        // Two special factors in one product.
        declines("Ci(b*x)*Si(b*x)");
        // Non-linear kernel arguments.
        declines("Ci(b*x^2)");
        // Positive-power denominators stay with the elementary stages.
        declines("Ci(b*x)/x");
        declines("Ei(b*x)/x");
        declines("exp(x)");
        // `sin(x)/x` belongs to the pre-existing Si table entry, so it is
        // deliberately *not* in this list; assert the owner instead.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = ocas_parse::parse(&ctx, "sin(x)/x").expect("parse");
        assert_eq!(
            special_integrate(&ctx, e, ctx.var("x"))
                .expect("Si entry")
                .to_string(),
            "Si(x)"
        );
    }

    // ---------------- budgets and builders ----------------

    #[test]
    fn budget_and_builder_edges() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let u = ctx.var("u");
        // Exhausted budgets decline instead of recursing forever.
        assert!(int_pow_exp(&ctx, u, ctx.num(1), 3, 0).is_none());
        assert!(int_pow_exp_quad(&ctx, u, 3, -1, 0).is_none());
        assert!(int_pow_trig(&ctx, u, 3, Kind::Sin, 0).is_none());
        assert!(int_exp_over_pow(&ctx, u, ctx.num(1), 2, 0).is_none());
        assert!(d_ei(&ctx, u, 1).is_none());
        assert!(d_order(&ctx, u, 2, 3, MAX_SPECIAL_STEPS).is_none());
        assert!(m_gauss(&ctx, u, 0, 2, MAX_SPECIAL_STEPS).is_none());
        // `0^0 = 1`, `0^k = 0` (k > 0), `(-1)^k` folded.
        assert_eq!(powi_clean(&ctx, ctx.num(0), 0).to_string(), "1");
        assert_eq!(powi_clean(&ctx, ctx.num(0), 2).to_string(), "0");
        assert_eq!(powi_clean(&ctx, ctx.num(-1), 3).to_string(), "-1");
        assert_eq!(powi_clean(&ctx, ctx.num(-1), 4).to_string(), "1");
        assert_eq!(muln(&ctx, &[ctx.num(1), ctx.num(1)]).to_string(), "1");
        assert_eq!(muln(&ctx, &[ctx.num(0), ctx.var("x")]).to_string(), "0");
        assert_eq!(addn(&ctx, &[ctx.num(0), ctx.num(0)]).to_string(), "0");
        // Polynomial extraction declines non-polynomials and deep degrees.
        let xv = ctx.var("x");
        assert!(poly_coeffs(&ctx, ctx.pow(xv, ctx.num(-1)), Symbol::new("x")).is_none());
        assert!(poly_coeffs(&ctx, ctx.pow(xv, ctx.num(9)), Symbol::new("x")).is_none());
        assert_eq!(binom_k(5, 2), 10);
        assert_eq!(binom_k(0, 0), 1);
        assert_eq!(binom_k(3, 4), 0);
    }
}
