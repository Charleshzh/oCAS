//! Rational-derivative kernel substitution (0.27.2 Wave C1).
//!
//! Kernels whose derivative is a *rational* function of the kernel itself
//! (`t = tan u`, `cot u`, `tanh u`, `coth u`; `dT/du = ±(1 ± T²)`) let a
//! trig/hyperbolic integrand be rewritten as an algebraic one in `t` without
//! introducing a second radical.
//!
//! For an integrand built from one of these kernels at a shared linear
//! argument `u = a·x + b` the substitution `t = T(u)` maps the integral to
//!
//! ```text
//! ∫ f(T(u)) dx = (1/a) · ∫ g(t) dt / (1 ± t²)
//! ```
//!
//! (`+` for the `tan`/`cot` family with the `cot` sign folded into the
//! constant factors, `−` for `tanh`/`coth`). The t-form is pure algebra: it
//! is split into its constant (in `t`) factors and a kernel-free core, the
//! core is integrated by the existing algebraic engines
//! ([`binomial`](super::binomial), then
//! [`sqrt_quadratic`](super::sqrt_quadratic), then
//! [`rational`](super::rational)), and the answer is back-substituted with
//! `t = T(u)`.
//!
//! Two local clean-ups make the algebraic core match the engines' shapes:
//! constant factors are pulled out of powers (`(d·t)^(1/2) → d^(1/2)·t^(1/2)`),
//! equal-base factors with rational exponents are merged
//! (`√(1 + t²)·(1 + t²)⁻¹ → (1 + t²)^(−1/2)`), integer exponents written as
//! unreduced fraction products (`2·2⁻¹`) are folded, and a bounded
//! distributive expansion retries the engines term-by-term.
//!
//! Shapes whose core is a single fractional-power binomial times a rational
//! factor (which `binomial` rationalizes internally) are handled by
//! `binomial` itself. The recurring half-integer family that remains —
//! `t^(k/2)·R(t²)` — is rationalized locally by `t = w²` and the resulting
//! rational form is completed by [`integrate_rational_complete`], a bounded
//! partial-fraction completion over `ℚ` (the ℚ engine declines every
//! square-free denominator of degree ≥ 3 and every irreducible quartic, and
//! the module must not route algebraic t-forms into `symbolic_rational`,
//! which is not bounded on them — it overflows the stack on `1/((2+3·t)(1+t²))`).
//!
//! **Certainty gate.** Every candidate answer is checked numerically before it
//! is emitted: the engine order (`binomial` → `sqrt_quadratic` → `rational`)
//! is arbitrated by a derivative check in the substitution variable, and the
//! back-substituted answer is checked once more against the integrand as it
//! was requested, over the same abscissae the project's numerical oracle uses.
//! A candidate that cannot be confirmed — including a branch-limited engine
//! form that only holds for `t > 0`, and anything the evaluator cannot model —
//! is declined, so a later stage keeps first claim rather than being masked by
//! an unverified answer. This costs nothing on the shapes the gate confirms
//! and is what keeps the stage from being a source of wrong answers.
//!
//! This module also owns the Pythagorean rewrites used by the integrator
//! (`1 − tanh² → sech²`, `1 − sin² → cos²`, `sec² − 1 → tan²`,
//! `1 + cot² → csc²`, `cosh² − sinh² → 1`, `1 + tan² → sec²`). They are
//! deliberately local to the integrator rather than `ocas_atom::normalize`
//! to keep the canonical form blast radius small. Every rewrite is an
//! identity, so an antiderivative of the rewritten integrand is one of the
//! original. A rewrite that would leave the rational-derivative family
//! (`1 − tanh² → sech²`) is only used when the family survives it; otherwise
//! the original in-family form is substituted directly.
//!
//! The family is exactly the four rational-derivative kernels: `sec`, `csc`,
//! `sech`, `csch` and `sin`/`cos`/`sinh`/`cosh` are declined (their
//! derivatives are irrational in the kernel, so substituting them would
//! introduce a second radical and would intercept shapes the rule table
//! already solves correctly).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_domain::{Domain, Integer, IntegerDomain, Rational, RationalDomain};
use ocas_poly::{DenseUnivariatePolynomial, Lex, SparseMultivariatePolynomial};

use crate::tower::convert::{GeneratorField, atom_to_rational, rational_to_atom};

use super::rules::rat_of;
use super::{
    binomial, contains_integral, gcd_i64, int_pow, inv, is_constant, lcm_i64, linear_form,
    node_count, pick_subst_symbol, rat_atom, rational, replace_symbol, sqrt_quadratic,
};

/// Dense univariate polynomial over `ℚ` (the local completion's workhorse).
type DPoly = DenseUnivariatePolynomial<RationalDomain>;
/// Sparse multivariate polynomial over `ℚ` (the atom converter's format).
type Sparse = SparseMultivariatePolynomial<RationalDomain, Lex>;

/// Node budget for the integrand and for every derived t-form.
const MAX_NODES: usize = 300;
/// Cap on the number of kernel sites in one integrand.
const MAX_SITES: usize = 24;
/// Cap on the exponent denominator of a kernel power (`t^(p/q)`).
const MAX_EXP_DEN: i64 = 6;
/// Fixpoint cap for the Pythagorean pre-pass.
const MAX_REWRITE_ROUNDS: usize = 12;
/// Recursion cap for the bounded-expansion retry.
const MAX_EXPAND_DEPTH: usize = 2;
/// Cap on the number of merged exponent groups in [`combine_like_powers`].
const MAX_POW_GROUPS: usize = 32;
/// Cap on the magnitude of a merged rational exponent.
const MAX_POW_EXP: i64 = 64;
/// Cap on the single-radical rationalization degree (lcm of the exponents'
/// denominators) in [`try_radical_rationalize`].
const MAX_RADICAL_DEG: i64 = 4;
/// Cap on the number of fractional-power sites accepted in a rationalization.
const MAX_RADICAL_SITES: usize = 8;
/// Cap on the denominator degree accepted by [`integrate_rational_complete`].
const MAX_RAT_DEG: usize = 8;
/// Cap on the numerator degree accepted by [`integrate_rational_complete`].
const MAX_RAT_NUM_DEG: usize = 12;
/// Cap on the number of irreducible factors and partial-fraction pieces.
const MAX_FACTORS: usize = 8;
const MAX_PIECES: usize = 16;

/// Reduced rational exponent `(p, q)` with `q > 0` for an exponent atom, or
/// `None` when the atom is not a rational number.
///
/// `rat_of` returns the raw `p`/`q` pair of the atom shape, which is *not*
/// reduced: the canonicalizer writes `2·2⁻¹` for the integer `1`, so `(2, 2)`
/// must be reduced before any `q == 1` / `q > 1` test.
fn rat_exponent(atom: Atom<'_>) -> Option<(i64, i64)> {
    let (p, q) = rat_of(atom)?;
    if q == 0 {
        return None;
    }
    let (p, q) = if q < 0 {
        (p.checked_neg()?, q.checked_neg()?)
    } else {
        (p, q)
    };
    let g = gcd_i64(p, q).max(1);
    Some((p / g, q / g))
}

/// Integrate `expr` via a rational-derivative kernel substitution.
///
/// Returns `None` when no such substitution applies (the caller continues
/// down the pipeline).
pub(crate) fn integrate_kernel_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    // Top-level sums are distributed over terms by the caller chain.
    if matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    // The Pythagorean identities are exact, so integrating the rewritten
    // integrand integrates the original one. A rewrite that would leave the
    // rational-derivative family (`1 − tanh² → sech²` introduces `sech`) is
    // only used when the family survives it: otherwise the original, already
    // in-family form is the better input (`(1 − tanh²)^(3/2)` substitutes
    // directly to `(1 − t²)^(3/2)`).
    let original = expr;
    let rewritten = pythagorean_rewrite(ctx, expr, var);
    let (expr, info) = match rewritten {
        Some(r) => match match_family(ctx, r, var) {
            Some(info) => (r, info),
            None => (expr, match_family(ctx, expr, var)?),
        },
        None => (expr, match_family(ctx, expr, var)?),
    };
    for primary in info.candidates() {
        if let Some(r) = try_primary(ctx, original, expr, var, &info, primary) {
            return Some(r);
        }
    }
    None
}

// =========================================================================
// Pythagorean rewrites
// =========================================================================

/// Apply the Pythagorean rewrites to a fixpoint (bounded). Returns `None`
/// when no rewrite applies.
fn pythagorean_rewrite<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    // The pipeline hands over normalized integrands; normalizing here too
    // keeps the pre-pass usable on raw parses (nested `Mul`/`Add` nodes would
    // otherwise hide the two-term patterns).
    let mut current = normalize(ctx, expr);
    let mut changed = false;
    for _ in 0..MAX_REWRITE_ROUNDS {
        let (next, hit) = rewrite_pass(ctx, current, var);
        if !hit {
            break;
        }
        changed = true;
        current = next;
    }
    if changed { Some(current) } else { None }
}

/// One bottom-up rewrite pass. The flag reports whether any identity fired.
fn rewrite_pass<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> (Atom<'a>, bool) {
    match expr.node() {
        AtomNode::Add(args) => {
            let mut kids = Vec::with_capacity(args.len());
            let mut changed = false;
            for a in args.iter() {
                let (na, hit) = rewrite_pass(ctx, *a, var);
                changed |= hit;
                kids.push(na);
            }
            if let Some(r) = rewrite_add(ctx, &kids, var) {
                return (r, true);
            }
            if !changed {
                return (expr, false);
            }
            (ctx.add(&kids), true)
        }
        AtomNode::Mul(args) => {
            let mut kids = Vec::with_capacity(args.len());
            let mut changed = false;
            for a in args.iter() {
                let (na, hit) = rewrite_pass(ctx, *a, var);
                changed |= hit;
                kids.push(na);
            }
            if !changed {
                return (expr, false);
            }
            (ctx.mul(&kids), true)
        }
        AtomNode::Pow(b, e) => {
            let (nb, c1) = rewrite_pass(ctx, *b, var);
            let (ne, c2) = rewrite_pass(ctx, *e, var);
            if !c1 && !c2 {
                return (expr, false);
            }
            (ctx.pow(nb, ne), true)
        }
        AtomNode::Num(_) | AtomNode::Var(_) | AtomNode::Fun(_, _) => (expr, false),
    }
}

/// The six identities on a two-term sum, in either argument order.
fn rewrite_add<'a>(ctx: &'a AtomArena<'a>, kids: &[Atom<'a>], var: Symbol) -> Option<Atom<'a>> {
    if kids.len() != 2 {
        return None;
    }
    let (x, y) = (kids[0], kids[1]);
    // `c ± c·T(u)²` against a constant term: the constant lead carries the
    // sign of the identity (`1 − tanh²` versus `1 + cot²`).
    if let Some(r) = rewrite_const_pair(ctx, x, y, var) {
        return Some(r);
    }
    if let Some(r) = rewrite_const_pair(ctx, y, x, var) {
        return Some(r);
    }
    // `c·sec(u)² − c` and `c·cosh(u)² − c·sinh(u)²`.
    if let Some(r) = rewrite_square_pair(ctx, x, y, var) {
        return Some(r);
    }
    rewrite_square_pair(ctx, y, x, var)
}

/// `c <op> c'·T(u)²` where `c` is constant: the `1 ± T²` identities.
///
/// `sign` is `+1` for the `1 + T² = S²` family (cot→csc, tan→sec) and `−1`
/// for the `1 − T² = S²` family (tanh→sech, sin→cos). The two coefficients
/// must be equal up to that sign.
fn rewrite_const_pair<'a>(
    ctx: &'a AtomArena<'a>,
    constant_term: Atom<'a>,
    square_term: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let c = constant_term_of(ctx, constant_term, var)?;
    let (c2, name, arg) = square_term_of(ctx, square_term, var)?;
    let neg_c = normalize(ctx, ctx.mul(&[ctx.num(-1), c]));
    let target = match name {
        "tanh" | "sin" => ("sech", "cos", neg_c),
        "cot" | "tan" => ("csc", "sec", c),
        _ => return None,
    };
    if c2 != target.2 {
        return None;
    }
    let result_name = if name == "tanh" {
        target.0
    } else if name == "sin" {
        target.1
    } else if name == "cot" {
        target.0
    } else {
        target.1
    };
    let sq = int_pow(ctx, ctx.fun(result_name, &[arg]), 2);
    Some(normalize(ctx, ctx.mul(&[c, sq])))
}

/// `c·T(u)² <op> c'·S(u)²` or `c·T(u)² <op> c'` (the `sec² − 1` and
/// `cosh² − sinh²` identities).
fn rewrite_square_pair<'a>(
    ctx: &'a AtomArena<'a>,
    first: Atom<'a>,
    second: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (c1, name1, arg1) = square_term_of(ctx, first, var)?;
    // `sec(u)² − c` → `c·tan(u)²`.
    if name1 == "sec"
        && let Some(c2) = constant_term_of(ctx, second, var)
        && c2 == normalize(ctx, ctx.mul(&[ctx.num(-1), c1]))
    {
        let sq = int_pow(ctx, ctx.fun("tan", &[arg1]), 2);
        return Some(normalize(ctx, ctx.mul(&[c1, sq])));
    }
    // `cosh(u)² − sinh(u)²` → `1` (scaled by the common coefficient).
    if name1 == "cosh"
        && let Some((c2, name2, arg2)) = square_term_of(ctx, second, var)
        && name2 == "sinh"
        && arg2 == arg1
        && c2 == normalize(ctx, ctx.mul(&[ctx.num(-1), c1]))
    {
        return Some(c1);
    }
    None
}

/// Match a `var`-constant term; returns it in normalized form.
fn constant_term_of<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if is_constant(term, var) {
        Some(normalize(ctx, term))
    } else {
        None
    }
}

/// Match `c·T(u)²` with `T` one of the Pythagorean kernels; returns
/// `(c, name, u)`.
fn square_term_of<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, &'static str, Atom<'a>)> {
    let (coeff, rest) = split_coeff(ctx, term, var);
    let AtomNode::Pow(base, exp) = rest.node() else {
        return None;
    };
    if !matches!(exp.node(), AtomNode::Num(2)) {
        return None;
    }
    let AtomNode::Fun(name, args) = base.node() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let name = match name.as_str() {
        "tan" => "tan",
        "cot" => "cot",
        "tanh" => "tanh",
        "sin" => "sin",
        "cos" => "cos",
        "sinh" => "sinh",
        "cosh" => "cosh",
        "sec" => "sec",
        "csc" => "csc",
        _ => return None,
    };
    Some((coeff, name, args[0]))
}

/// Split a product into its `var`-constant factors and the remaining factor
/// product (`c·rest`).
fn split_coeff<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>, var: Symbol) -> (Atom<'a>, Atom<'a>) {
    let AtomNode::Mul(args) = term.node() else {
        return (ctx.num(1), term);
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for a in args.iter() {
        if is_constant(*a, var) {
            coeff.push(*a);
        } else {
            rest.push(*a);
        }
    }
    let c = match coeff.len() {
        0 => ctx.num(1),
        1 => normalize(ctx, coeff[0]),
        _ => normalize(ctx, ctx.mul(&coeff)),
    };
    let r = match rest.len() {
        0 => ctx.num(1),
        1 => rest[0],
        _ => normalize(ctx, ctx.mul(&rest)),
    };
    (c, r)
}

// =========================================================================
// Kernel family detection
// =========================================================================

/// The kernel family shared by every kernel site of an integrand.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    /// `tan` / `cot` (`dT/du = ±(1 + T²)`).
    Trig,
    /// `tanh` / `coth` (`dT/du = 1 − T²`).
    Hyper,
}

/// The substitution kernel `t = T(u)`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Tan,
    Cot,
    Tanh,
    Coth,
}

impl Kernel {
    fn name(self) -> &'static str {
        match self {
            Kernel::Tan => "tan",
            Kernel::Cot => "cot",
            Kernel::Tanh => "tanh",
            Kernel::Coth => "coth",
        }
    }

    fn family(self) -> Family {
        match self {
            Kernel::Tan | Kernel::Cot => Family::Trig,
            Kernel::Tanh | Kernel::Coth => Family::Hyper,
        }
    }

    /// The Jacobian denominator `1 ± t²` of `dt/du` (sign folded by the
    /// caller for `cot`).
    fn quad<'a>(self, ctx: &'a AtomArena<'a>, t: Atom<'a>) -> Atom<'a> {
        let t2 = int_pow(ctx, t, 2);
        match self {
            Kernel::Tan | Kernel::Cot => ctx.add(&[ctx.num(1), t2]),
            Kernel::Tanh | Kernel::Coth => {
                ctx.add(&[ctx.num(1), normalize(ctx, ctx.mul(&[ctx.num(-1), t2]))])
            }
        }
    }
}

/// Shared kernel argument `u` and its slope `a` in `u = a·x + b`.
struct KernelInfo<'a> {
    u: Atom<'a>,
    slope: Atom<'a>,
    family: Family,
    /// Kernel names present in the integrand, in discovery order.
    seen: Vec<&'static str>,
}

impl<'a> KernelInfo<'a> {
    /// Substitution kernels to try: the ones literally present first (so the
    /// common case needs no rewriting of the other kernel), then the
    /// reciprocal of the family.
    fn candidates(&self) -> Vec<Kernel> {
        let all: [Kernel; 2] = match self.family {
            Family::Trig => [Kernel::Tan, Kernel::Cot],
            Family::Hyper => [Kernel::Tanh, Kernel::Coth],
        };
        let mut out: Vec<Kernel> = Vec::with_capacity(2);
        for k in all {
            if self.seen.iter().any(|n| *n == k.name()) {
                out.push(k);
            }
        }
        for k in all {
            if !out.contains(&k) {
                out.push(k);
            }
        }
        out
    }
}

/// Scan for a single rational-derivative kernel family at one shared linear
/// argument. Declines on any other function of `var`, on mixed
/// trig/hyperbolic kernels, on disagreeing arguments, and on nonlinear
/// arguments (`tan(c + d·x^(1/3))`).
fn match_family<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<KernelInfo<'a>> {
    let mut u: Option<Atom<'a>> = None;
    let mut family: Option<Family> = None;
    let mut seen: Vec<&'static str> = Vec::new();
    let mut sites = 0usize;
    scan_family(ctx, expr, var, &mut u, &mut family, &mut seen, &mut sites)?;
    let u = u?;
    let family = family?;
    let (slope, _b) = linear_form(ctx, u, var)?;
    if matches!(normalize(ctx, slope).node(), AtomNode::Num(0)) {
        return None;
    }
    Some(KernelInfo {
        u: normalize(ctx, u),
        slope: normalize(ctx, slope),
        family,
        seen,
    })
}

/// Recursive worker of [`match_family`].
fn scan_family<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    u: &mut Option<Atom<'a>>,
    family: &mut Option<Family>,
    seen: &mut Vec<&'static str>,
    sites: &mut usize,
) -> Option<()> {
    if is_constant(expr, var) {
        return Some(());
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => None,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                scan_family(ctx, *a, var, u, family, seen, sites)?;
            }
            Some(())
        }
        AtomNode::Pow(b, e) => {
            scan_family(ctx, *b, var, u, family, seen, sites)?;
            scan_family(ctx, *e, var, u, family, seen, sites)
        }
        AtomNode::Fun(name, args) => {
            // `sqrt` is a power, not a function of the argument: recurse into
            // its radicand so `sqrt(tan(u))` stays in family.
            if name.as_str() == "sqrt" && args.len() == 1 {
                return scan_family(ctx, args[0], var, u, family, seen, sites);
            }
            let (k, fam) = match name.as_str() {
                "tan" => (Kernel::Tan, Family::Trig),
                "cot" => (Kernel::Cot, Family::Trig),
                "tanh" => (Kernel::Tanh, Family::Hyper),
                "coth" => (Kernel::Coth, Family::Hyper),
                // Every other head — including `sec`/`csc`/`sech`/`csch`, whose
                // derivatives are *irrational* in the kernel, and
                // `sin`/`cos`/`sinh`/`cosh` — is outside the
                // rational-derivative family. Declining them here keeps the
                // stages that already solve them (rule table, trig
                // reductions) in first claim.
                _ => return None,
            };
            if args.len() != 1 {
                return None;
            }
            match *family {
                Some(f) if f != fam => return None,
                None => *family = Some(fam),
                _ => {}
            }
            let arg = normalize(ctx, args[0]);
            match *u {
                Some(u0) if u0 != arg => return None,
                None => *u = Some(arg),
                _ => {}
            }
            *sites += 1;
            if *sites > MAX_SITES {
                return None;
            }
            if !seen.contains(&k.name()) {
                seen.push(k.name());
            }
            Some(())
        }
    }
}

// =========================================================================
// Substitution to the kernel variable
// =========================================================================

/// Build the algebraic t-integrand for one choice of substitution kernel.
///
/// `original` is the integrand exactly as requested (before the Pythagorean
/// pre-pass); the candidate antiderivative is verified against it.
fn try_primary<'a>(
    ctx: &'a AtomArena<'a>,
    original: Atom<'a>,
    expr: Atom<'a>,
    var: Symbol,
    info: &KernelInfo<'a>,
    primary: Kernel,
) -> Option<Atom<'a>> {
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let substituted = subst_expr(ctx, expr, var, info, primary, t)?;
    if !is_admissible_t(substituted, var, t_sym) {
        return None;
    }
    let mut factors: Vec<Atom<'a>> = vec![substituted];
    // dx = du/a.
    let slope = info.slope;
    if !matches!(slope.node(), AtomNode::Num(1)) {
        factors.push(inv(ctx, slope));
    }
    // dx = dt/(±(1 ± t²)).
    factors.push(inv(ctx, primary.quad(ctx, t)));
    if primary == Kernel::Cot {
        factors.push(ctx.num(-1));
    }
    let t_form = normalize(ctx, ctx.mul(&factors));
    let t_form = combine_like_powers(ctx, t_form, t_sym);
    if node_count(t_form) > MAX_NODES || !is_admissible_t(t_form, var, t_sym) {
        return None;
    }
    let (consts, core) = split_constants(ctx, t_form, t_sym);
    // A `t`-constant core means the integrand is exactly that constant times
    // the kernel derivative: `∫ c·T'(u) du = c·T(u)`.
    let anti = if is_constant(core, t_sym) {
        if matches!(normalize(ctx, core).node(), AtomNode::Num(0)) {
            return None;
        }
        normalize(ctx, ctx.mul(&[core, t]))
    } else {
        let anti = integrate_algebraic(ctx, core, t_sym, 0)?;
        if contains_integral(anti) {
            return None;
        }
        anti
    };
    let back = replace_symbol(ctx, anti, t_sym, ctx.fun(primary.name(), &[info.u]));
    let back = normalize(ctx, back);
    let candidate = if consts.is_empty() {
        back
    } else {
        let mut out: Vec<Atom<'a>> = consts;
        out.push(back);
        normalize(ctx, ctx.mul(&out))
    };
    // Certainty gate: the candidate must differentiate back to the integrand
    // at a deterministic set of sample points (with deterministic dummy
    // values for the free parameters). An engine answer that cannot be
    // confirmed is declined rather than emitted.
    if verify_candidate(ctx, original, candidate, var) {
        Some(candidate)
    } else {
        None
    }
}

// =========================================================================
// Certainty gate
// =========================================================================

/// Deterministic dummy values for the free parameters of the integrand,
/// chosen positive and non-round so a formal mismatch cannot cancel
/// accidentally at the sample point.
const VERIFY_PARAMS: [f64; 10] = [0.7, 1.3, 0.37, 0.91, 1.73, 0.53, 1.11, 0.83, 1.47, 0.61];
/// Sample points for the final (integrand-level) gate. The spread matches the
/// project's numerical oracle (`ocas-tests::integral_eval`), so an answer this
/// gate accepts is not later flagged there: it deliberately includes negative
/// and larger abscissae, which is what exposes a branch-limited engine form
/// (`(x^{−n}+b)^{1/s}` answers that only hold for `x > 0`).
const VERIFY_SAMPLES: [f64; 8] = [-1.7, -0.9, -0.37, 0.31, 0.77, 1.23, 1.91, 2.53];
/// Sample points for the substitution-variable gate (the t-form). Signed and
/// away from the kernel poles.
const VERIFY_T_SAMPLES: [f64; 5] = [-0.77, -0.41, 0.31, 0.77, 1.23];
/// Relative tolerance of the derivative check (the project's numerical
/// oracle rejects above `1e-4` and accepts below `1e-5`; the gate rejects
/// above `1e-4` so it never fails an answer the oracle would accept).
const VERIFY_TOL: f64 = 1e-4;
/// Minimum number of usable sample points for an answer to be considered
/// verified (the rest may lie outside a radical's real domain).
const VERIFY_MIN_SAMPLES: usize = 2;

/// Confirm `candidate' == form` numerically before emitting it.
///
/// The derivative is taken with the same 5-point central stencil the project's
/// numerical oracle uses (`h = 1e-4·max(1, |x|)`) rather than symbolically:
/// `diff` on a large candidate is itself a deep recursion, and this gate sits
/// inside the integration chain, where stack depth is at a premium.
///
/// The gate is deliberately conservative in one direction only: a candidate
/// is accepted when at least [`VERIFY_MIN_SAMPLES`] sample points give finite
/// real values on both sides and every one of them agrees within
/// [`VERIFY_TOL`]; anything else — a mismatch, or too few usable points —
/// declines. It cannot prove an antiderivative, but it does catch a wrong
/// engine answer, which would otherwise mask a correct later stage.
fn verify_derivative<'a>(
    _ctx: &'a AtomArena<'a>,
    form: Atom<'a>,
    candidate: Atom<'a>,
    var: Symbol,
    samples: &[f64],
) -> bool {
    let mut symbols: Vec<Symbol> = Vec::new();
    collect_symbols(form, var, &mut symbols);
    collect_symbols(candidate, var, &mut symbols);
    symbols.sort_by_key(|s| s.as_str().to_string());
    symbols.dedup();
    if symbols.len() > VERIFY_PARAMS.len() {
        return false;
    }
    let mut usable = 0usize;
    for &xv in samples.iter() {
        let mut env: Vec<(Symbol, f64)> = symbols
            .iter()
            .enumerate()
            .map(|(i, s)| (*s, VERIFY_PARAMS[i]))
            .collect();
        env.push((var, xv));
        let Some(rhs) = eval_real(form, &env) else {
            continue;
        };
        if !rhs.is_finite() {
            continue;
        }
        let h = 1e-4 * xv.abs().max(1.0);
        let mut stencil = [0.0f64; 4];
        let mut ok = true;
        for (slot, factor) in stencil.iter_mut().zip([-2.0f64, -1.0, 1.0, 2.0]) {
            if let Some((_, last)) = env.last_mut() {
                *last = xv + factor * h;
            }
            let Some(v) = eval_real(candidate, &env) else {
                ok = false;
                break;
            };
            if !v.is_finite() {
                ok = false;
                break;
            }
            *slot = v;
        }
        if !ok {
            continue;
        }
        // (f(x−2h) − 8f(x−h) + 8f(x+h) − f(x+2h)) / (12h).
        let lhs = (stencil[0] - 8.0 * stencil[1] + 8.0 * stencil[2] - stencil[3]) / (12.0 * h);
        usable += 1;
        if (lhs - rhs).abs() > VERIFY_TOL * rhs.abs().max(1.0) {
            if std::env::var_os("OCAS_KERNEL_VERIFY_DEBUG").is_some() {
                eprintln!(
                    "[kernel_subst verify] mismatch on {form} at {var:?}={xv}:\n  \
                     diff={lhs} form={rhs}\n  candidate: {candidate}"
                );
            }
            return false;
        }
    }
    if usable < VERIFY_MIN_SAMPLES && std::env::var_os("OCAS_KERNEL_VERIFY_DEBUG").is_some() {
        eprintln!("[kernel_subst verify] only {usable} usable samples for {form}\n  {candidate}");
    }
    usable >= VERIFY_MIN_SAMPLES
}

/// Gate an answer for the t-form `form` (substitution-variable samples).
fn verify_local<'a>(
    ctx: &'a AtomArena<'a>,
    form: Atom<'a>,
    candidate: Atom<'a>,
    var: Symbol,
) -> bool {
    verify_derivative(ctx, form, candidate, var, &VERIFY_T_SAMPLES)
}

/// Final gate: the back-substituted answer must differentiate back to the
/// integrand as it was requested.
fn verify_candidate<'a>(
    ctx: &'a AtomArena<'a>,
    integrand: Atom<'a>,
    candidate: Atom<'a>,
    var: Symbol,
) -> bool {
    verify_derivative(ctx, integrand, candidate, var, &VERIFY_SAMPLES)
}

/// Collect the symbols occurring in `expr`, excluding the integration
/// variable `var`.
fn collect_symbols<'a>(expr: Atom<'a>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Num(_) => {}
        AtomNode::Var(v) => {
            if *v != var {
                out.push(*v);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_symbols(*a, var, out);
            }
        }
        AtomNode::Pow(b, e) => {
            collect_symbols(*b, var, out);
            collect_symbols(*e, var, out);
        }
    }
}

/// `f64` evaluator for the certainty gate. Returns `None` for a head it does
/// not model, so an exotic answer is declined rather than mis-checked.
fn eval_real(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
    match expr.node() {
        AtomNode::Num(n) => Some(*n as f64),
        AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, val)| *val),
        AtomNode::Add(args) => args
            .iter()
            .try_fold(0.0, |acc, a| Some(acc + eval_real(*a, env)?)),
        AtomNode::Mul(args) => args
            .iter()
            .try_fold(1.0, |acc, a| Some(acc * eval_real(*a, env)?)),
        AtomNode::Pow(b, e) => Some(eval_real(*b, env)?.powf(eval_real(*e, env)?)),
        AtomNode::Fun(name, args) => {
            let v = eval_real(*args.first()?, env)?;
            Some(match name.as_str() {
                "sin" => v.sin(),
                "cos" => v.cos(),
                "tan" => v.tan(),
                "cot" => 1.0 / v.tan(),
                "sec" => 1.0 / v.cos(),
                "csc" => 1.0 / v.sin(),
                "sinh" => v.sinh(),
                "cosh" => v.cosh(),
                "tanh" => v.tanh(),
                "coth" => 1.0 / v.tanh(),
                "sech" => 1.0 / v.cosh(),
                "csch" => 1.0 / v.sinh(),
                "exp" => v.exp(),
                // |·| keeps the real branch finite; log|u| and log(u) differ
                // by a constant on each interval, so the derivative matches.
                "log" => v.abs().ln(),
                "sqrt" => v.sqrt(),
                "atan" => v.atan(),
                "asin" => v.asin(),
                "acos" => v.acos(),
                "atanh" => v.atanh(),
                "asinh" => v.asinh(),
                "acosh" => v.acosh(),
                _ => return None,
            })
        }
    }
}

/// Map every family kernel of `expr` to its algebraic image in `t`, leaving
/// constants untouched. Returns `None` when `var` survives outside a kernel
/// argument.
fn subst_expr<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    info: &KernelInfo<'a>,
    primary: Kernel,
    t: Atom<'a>,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) {
        return Some(expr);
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => None,
        AtomNode::Add(args) => {
            let mut kids = Vec::with_capacity(args.len());
            for a in args.iter() {
                kids.push(subst_expr(ctx, *a, var, info, primary, t)?);
            }
            Some(ctx.add(&kids))
        }
        AtomNode::Mul(args) => {
            let mut kids = Vec::with_capacity(args.len());
            for a in args.iter() {
                kids.push(subst_expr(ctx, *a, var, info, primary, t)?);
            }
            Some(ctx.mul(&kids))
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, args) = b.node()
                && args.len() == 1
                && let Some(img) = kernel_image(ctx, name.as_str(), primary, t)
            {
                if normalize(ctx, args[0]) != info.u {
                    return None;
                }
                let ne = subst_expr(ctx, *e, var, info, primary, t)?;
                return Some(normalize(ctx, ctx.pow(img, ne)));
            }
            let nb = subst_expr(ctx, *b, var, info, primary, t)?;
            let ne = subst_expr(ctx, *e, var, info, primary, t)?;
            Some(normalize(ctx, ctx.pow(nb, ne)))
        }
        AtomNode::Fun(name, args) => {
            if args.len() != 1 {
                return None;
            }
            // `sqrt(g)` is `g^(1/2)`: represent it as a power so the
            // algebraic dispatch sees a single power grammar.
            if name.as_str() == "sqrt" {
                let inner = subst_expr(ctx, args[0], var, info, primary, t)?;
                return Some(normalize(ctx, ctx.pow(inner, rat_atom(ctx, 1, 2))));
            }
            let img = kernel_image(ctx, name.as_str(), primary, t)?;
            if normalize(ctx, args[0]) != info.u {
                return None;
            }
            Some(img)
        }
    }
}

/// The image of a family kernel under `t = primary(u)`: the kernel itself or
/// its reciprocal. Only the four rational-derivative kernels are in the
/// family, so no radical is ever introduced by the substitution itself.
fn kernel_image<'a>(
    ctx: &'a AtomArena<'a>,
    name: &str,
    primary: Kernel,
    t: Atom<'a>,
) -> Option<Atom<'a>> {
    let in_family = matches!(name, "tan" | "cot" | "tanh" | "coth");
    let same_group = match primary.family() {
        Family::Trig => matches!(name, "tan" | "cot"),
        Family::Hyper => matches!(name, "tanh" | "coth"),
    };
    if !in_family || !same_group {
        return None;
    }
    if name == primary.name() {
        return Some(t);
    }
    Some(inv(ctx, t))
}

/// Structural gate on a candidate t-form: no leftover integration variable,
/// no functions of `t`, and no symbolic exponent on a `t`-dependent base
/// (`(b·t²)^n` is not an algebraic rational form).
fn is_admissible_t<'a>(expr: Atom<'a>, var: Symbol, t_sym: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v != var,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| is_admissible_t(*a, var, t_sym))
        }
        AtomNode::Pow(b, e) => {
            if !is_admissible_t(*e, var, t_sym) {
                return false;
            }
            if is_constant(*b, t_sym) {
                return true;
            }
            match rat_exponent(*e) {
                Some((_p, q)) => q <= MAX_EXP_DEN && is_admissible_t(*b, var, t_sym),
                None => false,
            }
        }
        AtomNode::Fun(_, _) => is_constant(expr, t_sym),
    }
}

// =========================================================================
// t-form clean-up
// =========================================================================

/// Distribute a rational power over its constant (in `t`) factors, so
/// `(d·t)^(1/2)` becomes `d^(1/2)·t^(1/2)` and `d^(1/2)` can be pulled out
/// as a constant multiplier. The `t`-dependent part is left untouched, so
/// no branch of a `t`-power is altered.
fn split_constant_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, t_sym: Symbol) -> Atom<'a> {
    match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_constant_powers(ctx, *a, t_sym))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_constant_powers(ctx, *a, t_sym))
                .collect();
            ctx.mul(&kids)
        }
        AtomNode::Pow(b, e) => {
            let nb = split_constant_powers(ctx, *b, t_sym);
            let ne = split_constant_powers(ctx, *e, t_sym);
            if let Some((_p, q)) = rat_exponent(ne)
                && q > 1
                && let AtomNode::Mul(factors) = nb.node()
            {
                let mut consts: Vec<Atom<'a>> = Vec::new();
                let mut rest: Vec<Atom<'a>> = Vec::new();
                for f in factors.iter() {
                    if is_constant(*f, t_sym) {
                        consts.push(*f);
                    } else {
                        rest.push(*f);
                    }
                }
                if !consts.is_empty() && !rest.is_empty() {
                    let c = normalize(ctx, ctx.pow(normalize(ctx, ctx.mul(&consts)), ne));
                    let r = normalize(ctx, ctx.pow(normalize(ctx, ctx.mul(&rest)), ne));
                    return normalize(ctx, ctx.mul(&[c, r]));
                }
            }
            normalize(ctx, ctx.pow(nb, ne))
        }
        AtomNode::Num(_) | AtomNode::Var(_) | AtomNode::Fun(_, _) => expr,
    }
}

/// Merge equal-base factors with numeric rational exponents, bottom-up.
fn combine_like_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, t_sym: Symbol) -> Atom<'a> {
    let expr = split_constant_powers(ctx, expr, t_sym);
    match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| combine_like_powers(ctx, *a, t_sym))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| combine_like_powers(ctx, *a, t_sym))
                .collect();
            combine_mul(ctx, &kids, t_sym)
        }
        AtomNode::Pow(b, e) => {
            let nb = combine_like_powers(ctx, *b, t_sym);
            let ne = combine_like_powers(ctx, *e, t_sym);
            // Fold exponents the canonicalizer leaves written as fraction
            // products (`2·2⁻¹` reads as the integer 1): an integer exponent
            // is exact for every base, so this is branch-safe.
            match rat_exponent(ne) {
                Some((p, 1)) => int_pow(ctx, nb, p),
                Some((p, q)) => {
                    // `(u^(2k))^(p/q) = u^(2kp/q)` is also exact when the
                    // folded exponent is an *even* integer: the radicand is a
                    // square, so the principal root drops the sign. An odd
                    // folded exponent would need `|u|` and is left alone.
                    if let AtomNode::Pow(inner, inner_exp) = nb.node()
                        && let Some((m, 1)) = rat_exponent(*inner_exp)
                        && m % 2 == 0
                        && let Some(n) = m.checked_mul(p)
                        && n % q == 0
                        && (n / q) % 2 == 0
                    {
                        return int_pow(ctx, *inner, n / q);
                    }
                    ctx.pow(nb, ne)
                }
                None => ctx.pow(nb, ne),
            }
        }
        AtomNode::Num(_) | AtomNode::Var(_) | AtomNode::Fun(_, _) => expr,
    }
}

/// Merge the factors of one product (see [`combine_like_powers`]).
fn combine_mul<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], t_sym: Symbol) -> Atom<'a> {
    let mut groups: Vec<(Atom<'a>, i64, i64)> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        let Some((base, p, q)) = groupable_factor(*f, t_sym) else {
            rest.push(*f);
            continue;
        };
        if let Some(slot) = groups.iter_mut().find(|(b, _, _)| *b == base) {
            let num = slot
                .1
                .checked_mul(q)
                .and_then(|v| v.checked_add(p * slot.2));
            let den = slot.2.checked_mul(q);
            let (Some(num), Some(den)) = (num, den) else {
                rest.push(*f);
                continue;
            };
            let g = gcd_i64(num, den).max(1);
            slot.1 = num / g;
            slot.2 = den / g;
        } else if groups.len() < MAX_POW_GROUPS {
            groups.push((base, p, q));
        } else {
            rest.push(*f);
        }
    }
    let mut out: Vec<Atom<'a>> = Vec::with_capacity(rest.len() + groups.len());
    for (base, p, q) in groups {
        if p == 0 {
            continue;
        }
        if p.abs() > MAX_POW_EXP || q.abs() > MAX_POW_EXP {
            out.push(base);
            continue;
        }
        if q == 1 {
            out.push(int_pow(ctx, base, p));
        } else {
            out.push(ctx.pow(base, rat_atom(ctx, p, q)));
        }
    }
    out.extend(rest);
    if out.is_empty() {
        // Every factor cancelled (`(1 − t²)/(1 − t²)`): the product is 1.
        return ctx.num(1);
    }
    normalize(ctx, ctx.mul(&out))
}

/// A `t`-dependent factor that carries a numeric rational exponent, as
/// `base^(p/q)`; `q = 1` for a bare factor.
fn groupable_factor<'a>(f: Atom<'a>, t_sym: Symbol) -> Option<(Atom<'a>, i64, i64)> {
    match f.node() {
        AtomNode::Pow(b, e) if !is_constant(*b, t_sym) => {
            let (p, q) = rat_exponent(*e)?;
            Some((*b, p, q))
        }
        _ if !is_constant(f, t_sym) => Some((f, 1, 1)),
        _ => None,
    }
}

/// Partition a product into its `t`-constant factors and the remaining
/// core (which is `1` when the whole product is constant).
fn split_constants<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    t_sym: Symbol,
) -> (Vec<Atom<'a>>, Atom<'a>) {
    let AtomNode::Mul(args) = expr.node() else {
        if is_constant(expr, t_sym) {
            return (vec![expr], ctx.num(1));
        }
        return (Vec::new(), expr);
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for a in args.iter() {
        if is_constant(*a, t_sym) {
            consts.push(*a);
        } else {
            rest.push(*a);
        }
    }
    let core = match rest.len() {
        0 => ctx.num(1),
        1 => rest[0],
        _ => normalize(ctx, ctx.mul(&rest)),
    };
    (consts, core)
}

// =========================================================================
// Algebraic dispatch
// =========================================================================

/// Integrate a kernel-free algebraic form in `var` via the algebraic engines.
fn integrate_algebraic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    depth: usize,
) -> Option<Atom<'a>> {
    if node_count(expr) > MAX_NODES {
        return None;
    }
    let (consts, core) = split_constants(ctx, expr, var);
    if is_constant(core, var) {
        return None;
    }
    if let Some(r) = try_engines(ctx, core, var) {
        return Some(scale_by(ctx, &consts, r));
    }
    if depth < MAX_EXPAND_DEPTH
        && let Some(expanded) = crate::expand::expand_bounded(ctx, core)
    {
        let folded = crate::ode::util::collect_terms(ctx, expanded);
        if let Some(r) = integrate_terms(ctx, folded, var, depth + 1)
            && verify_local(ctx, core, r, var)
        {
            return Some(scale_by(ctx, &consts, r));
        }
    }
    if let Some(r) = try_radical_rationalize(ctx, core, var)
        && verify_local(ctx, core, r, var)
    {
        return Some(scale_by(ctx, &consts, r));
    }
    if let Some(r) = integrate_rational_complete(ctx, core, var)
        && verify_local(ctx, core, r, var)
    {
        return Some(scale_by(ctx, &consts, r));
    }
    None
}

/// `consts · r`, folded.
fn scale_by<'a>(ctx: &'a AtomArena<'a>, consts: &[Atom<'a>], r: Atom<'a>) -> Atom<'a> {
    if consts.is_empty() {
        return r;
    }
    let mut factors: Vec<Atom<'a>> = consts.to_vec();
    factors.push(r);
    normalize(ctx, ctx.mul(&factors))
}

/// Integrate every additive term of an expanded form.
fn integrate_terms<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    depth: usize,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Add(args) => {
            let mut out: Vec<Atom<'a>> = Vec::with_capacity(args.len());
            for a in args.iter() {
                out.push(integrate_algebraic(ctx, *a, var, depth)?);
            }
            Some(normalize(ctx, ctx.add(&out)))
        }
        _ => integrate_algebraic(ctx, expr, var, depth),
    }
}

/// The three algebraic engines, in order of preference. An engine answer that
/// still contains an `Integral` residue, or that fails the local derivative
/// gate, is skipped in favour of the next engine: the gate is what stops a
/// branch-limited closed form (`(x^{−n}+b)^{1/s}` answers that hold only for
/// `x > 0`) from being emitted when a globally valid alternative exists.
fn try_engines<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if let Some(r) = binomial::integrate_binomial(ctx, expr, var)
        && !contains_integral(r)
        && verify_local(ctx, expr, r, var)
    {
        return Some(r);
    }
    if let Some(r) = sqrt_quadratic::integrate_sqrt_quadratic(ctx, expr, var)
        && !contains_integral(r)
        && verify_local(ctx, expr, r, var)
    {
        return Some(r);
    }
    if let Some(r) = rational::integrate_rational(ctx, expr, var)
        && !contains_integral(r)
        && verify_local(ctx, expr, r, var)
    {
        return Some(r);
    }
    None
}

// =========================================================================
// Single-radical rationalization + local rational completion
// =========================================================================

/// Rationalize the single fractional-power site `g^(k/s)` of `expr` (base `g`
/// linear in `var`) by `var = (w^s − α)/β`, leaving a rational `w`-form that
/// [`integrate_rational_complete`] can finish. Answers are back-substituted
/// `w = g^(1/s)`.
fn try_radical_rationalize<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let mut base: Option<Atom<'a>> = None;
    let mut dens: Vec<i64> = Vec::new();
    let mut sites = 0usize;
    collect_radical_sites(ctx, expr, var, &mut base, &mut dens, &mut sites)?;
    let g = base?;
    let (beta, alpha) = linear_form(ctx, g, var)?;
    if matches!(normalize(ctx, beta).node(), AtomNode::Num(0)) {
        return None;
    }
    let l = dens.iter().try_fold(1i64, |acc, s| lcm_i64(acc, *s))?;
    if !(2..=MAX_RADICAL_DEG).contains(&l) {
        return None;
    }
    let w_sym = pick_subst_symbol(expr, var)?;
    let w = ctx.var(w_sym.as_str());
    // var = (w^l − α)/β.
    let var_w = normalize(
        ctx,
        ctx.mul(&[
            normalize(
                ctx,
                ctx.add(&[
                    int_pow(ctx, w, l),
                    normalize(ctx, ctx.mul(&[ctx.num(-1), alpha])),
                ]),
            ),
            inv(ctx, beta),
        ]),
    );
    let subbed = subst_radical(ctx, expr, var, g, l, w, var_w)?;
    // dvar/dw = l·w^(l−1)/β.
    let dvar = normalize(
        ctx,
        ctx.mul(&[ctx.num(l), inv(ctx, beta), int_pow(ctx, w, l - 1)]),
    );
    let w_form = normalize(ctx, ctx.mul(&[subbed, dvar]));
    if node_count(w_form) > MAX_NODES || !is_rational_form(w_form, var) {
        return None;
    }
    let anti = integrate_rational_complete(ctx, w_form, w_sym)?;
    let back = replace_symbol(ctx, anti, w_sym, ctx.pow(g, rat_atom(ctx, 1, l)));
    Some(normalize(ctx, back))
}

/// Collect the fractional-power sites of `expr`: `Pow(base, k/s)` with `s > 1`
/// and `Fun("sqrt", [base])`, plus `Pow(Mul(c, base), k/s)` whose constant part
/// has already been split off. All sites must share one (normalized) base.
fn collect_radical_sites<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    base: &mut Option<Atom<'a>>,
    dens: &mut Vec<i64>,
    sites: &mut usize,
) -> Option<()> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if let Some((_, s)) = rat_exponent(*e)
                && s > 1
                && !is_constant(*b, var)
            {
                if s > MAX_RADICAL_DEG {
                    return None;
                }
                let nb = normalize(ctx, *b);
                match *base {
                    Some(g0) if nb != g0 => return None,
                    None => *base = Some(nb),
                    _ => {}
                }
                dens.push(s);
                *sites += 1;
                if *sites > MAX_RADICAL_SITES {
                    return None;
                }
            }
            collect_radical_sites(ctx, *b, var, base, dens, sites)?;
            collect_radical_sites(ctx, *e, var, base, dens, sites)
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            if !is_constant(args[0], var) {
                let nb = normalize(ctx, args[0]);
                match *base {
                    Some(g0) if nb != g0 => return None,
                    None => *base = Some(nb),
                    _ => {}
                }
                dens.push(2);
                *sites += 1;
                if *sites > MAX_RADICAL_SITES {
                    return None;
                }
            }
            collect_radical_sites(ctx, args[0], var, base, dens, sites)
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_radical_sites(ctx, *a, var, base, dens, sites)?;
            }
            Some(())
        }
        AtomNode::Num(_) | AtomNode::Var(_) => Some(()),
    }
}

/// Substitute `var → var_w` and every site `g^(k/s) → w^(k·l/s)` in `expr`
/// (integer exponents, since `s | l`), converting `sqrt` to a power.
fn subst_radical<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    g: Atom<'a>,
    l: i64,
    w: Atom<'a>,
    var_w: Atom<'a>,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        AtomNode::Var(v) => {
            if *v == var {
                Some(var_w)
            } else {
                Some(expr)
            }
        }
        AtomNode::Add(args) => {
            let mut kids = Vec::with_capacity(args.len());
            for a in args.iter() {
                kids.push(subst_radical(ctx, *a, var, g, l, w, var_w)?);
            }
            Some(ctx.add(&kids))
        }
        AtomNode::Mul(args) => {
            let mut kids = Vec::with_capacity(args.len());
            for a in args.iter() {
                kids.push(subst_radical(ctx, *a, var, g, l, w, var_w)?);
            }
            Some(ctx.mul(&kids))
        }
        AtomNode::Pow(b, e) => {
            if normalize(ctx, *b) == g
                && let Some((k, s)) = rat_exponent(*e)
                && s > 1
            {
                return Some(int_pow(ctx, w, k.checked_mul(l / s)?));
            }
            let nb = subst_radical(ctx, *b, var, g, l, w, var_w)?;
            let ne = subst_radical(ctx, *e, var, g, l, w, var_w)?;
            Some(normalize(ctx, ctx.pow(nb, ne)))
        }
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            if normalize(ctx, args[0]) == g {
                return Some(int_pow(ctx, w, l / 2));
            }
            let na = subst_radical(ctx, args[0], var, g, l, w, var_w)?;
            Some(normalize(ctx, ctx.pow(na, rat_atom(ctx, 1, 2))))
        }
        AtomNode::Fun(name, args) => {
            let mut kids = Vec::with_capacity(args.len());
            for a in args.iter() {
                kids.push(subst_radical(ctx, *a, var, g, l, w, var_w)?);
            }
            Some(ctx.fun(name.as_str(), &kids))
        }
    }
}

/// Structural check that `expr` is a rational form in `var`: only integer
/// powers, no functions, and no leftover occurrence of the original
/// integration variable's kernel form.
fn is_rational_form<'a>(expr: Atom<'a>, orig_var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v != orig_var,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| is_rational_form(*a, orig_var))
        }
        AtomNode::Pow(b, e) => {
            matches!(e.node(), AtomNode::Num(_)) && is_rational_form(*b, orig_var)
        }
        AtomNode::Fun(_, _) => false,
    }
}

/// Integrate a rational function of `var` over `ℚ` by an exact
/// partial-fraction completion, delegating every piece to the ℚ rational
/// engine. Covers what the engine itself declines as a whole: square-free
/// denominators mixing linear and quadratic factors, repeated quadratic
/// factors, and even quartics (via the `ℚ(√d)` split below).
fn integrate_rational_complete<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let rf = atom_to_rational(expr, &[x])?;
    let num = sparse_to_dense(&rf.numerator)?;
    let den = sparse_to_dense(&rf.denominator)?;
    if den.is_zero() || den.degree()? > MAX_RAT_DEG || num.degree().unwrap_or(0) > MAX_RAT_NUM_DEG {
        return None;
    }
    let (quotient, remainder) = num.div_rem(&den)?;
    let mut parts: Vec<Atom<'a>> = Vec::new();
    if !quotient.is_zero() {
        parts.push(poly_atom(ctx, &rational::poly_integrate(&quotient), x)?);
    }
    if !remainder.is_zero() {
        let lc = den.lcoeff();
        let inv_lc = RationalDomain.inv(&lc)?;
        let den = den.mul_scalar(&inv_lc);
        let remainder = remainder.mul_scalar(&inv_lc);
        let (num, den) = reduce_pair(&remainder, &den)?;
        if !num.is_zero() {
            let even = if den.degree() == Some(4) && is_even_poly(&den) && is_even_poly(&num) {
                even_quartic_atom(ctx, &num, &den, x)
            } else {
                None
            };
            let piece = match even {
                Some(p) => p,
                None => partial_fraction_atom(ctx, &num, &den, x, var)?,
            };
            parts.push(piece);
        }
    }
    match parts.len() {
        0 => Some(ctx.num(0)),
        1 => Some(parts.remove(0)),
        _ => Some(normalize(ctx, ctx.add(&parts))),
    }
}

/// Cancel the polynomial gcd of `num`/`den` (best effort).
fn reduce_pair(num: &DPoly, den: &DPoly) -> Option<(DPoly, DPoly)> {
    let cloned = || (num.clone(), den.clone());
    if num.is_zero() {
        return Some(cloned());
    }
    let (g, _, _) = num.extended_gcd_poly(den);
    if g.is_zero() || g.degree() == Some(0) {
        return Some(cloned());
    }
    let (nq, nr) = num.div_rem(&g)?;
    let (dq, dr) = den.div_rem(&g)?;
    if !nr.is_zero() || !dr.is_zero() {
        return Some(cloned());
    }
    Some((nq, dq))
}

/// True when only even powers of the variable occur.
fn is_even_poly(p: &DPoly) -> bool {
    p.coeffs()
        .iter()
        .enumerate()
        .all(|(i, c)| i % 2 == 0 || RationalDomain.is_zero(c))
}

/// Full partial-fraction decomposition of `num/den` (denominator monic with
/// irreducible factors of degree ≤ 2), each piece integrated by the ℚ engine.
fn partial_fraction_atom<'a>(
    ctx: &'a AtomArena<'a>,
    num: &DPoly,
    den: &DPoly,
    x: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let factors = factor_den(den)?;
    if factors.is_empty() || factors.len() > MAX_FACTORS {
        return None;
    }
    let (poly, pieces) = peel(num, den, &factors)?;
    if pieces.len() > MAX_PIECES {
        return None;
    }
    let mut parts: Vec<Atom<'a>> = Vec::new();
    if !poly.is_zero() {
        parts.push(poly_atom(ctx, &rational::poly_integrate(&poly), x)?);
    }
    for (n, d) in pieces {
        let piece = ratio_atom(ctx, &n, &d, x)?;
        let integrated = rational::integrate_rational(ctx, piece, var)?;
        if contains_integral(integrated) {
            return None;
        }
        parts.push(integrated);
    }
    match parts.len() {
        0 => Some(ctx.num(0)),
        1 => Some(parts.remove(0)),
        _ => Some(normalize(ctx, ctx.add(&parts))),
    }
}

/// Factor a monic denominator into `(irreducible factor, multiplicity)` pairs
/// with every factor of degree ≤ 2. Declines (rather than guessing) on a
/// higher-degree irreducible factor.
fn factor_den(den: &DPoly) -> Option<Vec<(DPoly, usize)>> {
    let mut out: Vec<(DPoly, usize)> = Vec::new();
    for (squarefree, mult) in den.square_free_factorization() {
        if squarefree.degree() == Some(0) {
            continue;
        }
        let irreducibles = factor_squarefree(&squarefree)?;
        if out.len() + irreducibles.len() > MAX_FACTORS {
            return None;
        }
        for f in irreducibles {
            out.push((f, mult));
        }
    }
    Some(out)
}

/// Factor a monic square-free polynomial over `ℚ` into monic irreducibles,
/// declining when any factor has degree > 2.
fn factor_squarefree(f: &DPoly) -> Option<Vec<DPoly>> {
    // Clear denominators and go through the integer factorization.
    let mut lcm: i64 = 1;
    for c in f.coeffs() {
        lcm = num_integer::lcm(lcm, c.denom().to_i64()?);
    }
    let mut zcoeffs: Vec<Integer> = Vec::with_capacity(f.coeffs().len());
    for c in f.coeffs() {
        let scaled = RationalDomain.mul(c, &Rational::new(lcm, 1));
        zcoeffs.push(Integer::from(scaled.numer().to_i64()?));
    }
    let zpoly = DenseUnivariatePolynomial::from_coeffs(IntegerDomain, zcoeffs);
    let factors = zpoly.primitive_part().factor();
    let mut out: Vec<DPoly> = Vec::with_capacity(factors.len());
    for (fac, _mult) in &factors {
        let lc = fac.coeffs().last()?.to_i64()?;
        if lc == 0 {
            return None;
        }
        let coeffs: Option<Vec<Rational>> = fac
            .coeffs()
            .iter()
            .map(|c| Some(Rational::new(c.to_i64()?, lc)))
            .collect();
        let dense = DPoly::from_coeffs(RationalDomain, coeffs?);
        if dense.degree()? > 2 {
            return None;
        }
        out.push(dense);
    }
    Some(out)
}

/// Peel `num/den` into a polynomial part and the partial fractions
/// `A/f^j` for every irreducible factor `f^j || den` (deg `f` ≤ 2).
///
/// The peeling uses the Bézout identity `s·f^m + t·g = 1` for `den = f^m·g`:
/// `num/den = num·t/f^m + num·s/g`, so the `f^m` part drops out and the
/// cofactor part is processed recursively.
fn peel(
    num: &DPoly,
    den: &DPoly,
    factors: &[(DPoly, usize)],
) -> Option<(DPoly, Vec<(DPoly, DPoly)>)> {
    let zero = || DPoly::from_coeffs(RationalDomain, vec![]);
    let mut remaining = den.clone();
    let mut current = num.clone();
    let mut poly = zero();
    let mut pieces: Vec<(DPoly, DPoly)> = Vec::new();
    for (f, mult) in factors {
        let fm = f.pow(*mult as u32);
        let (cofactor, rem) = remaining.div_rem(&fm)?;
        if !rem.is_zero() {
            return None;
        }
        if cofactor.degree() == Some(0) {
            // Last factor: the unit cofactor only scales the numerator.
            let c = cofactor.coeffs().first().cloned()?;
            let inv_c = RationalDomain.inv(&c)?;
            current = current.mul_scalar(&inv_c);
            let (quotient, proper) = current.div_rem(&fm)?;
            poly = poly.add(&quotient);
            push_power_pieces(f, *mult, &proper, &mut pieces)?;
            current = zero();
            continue;
        }
        // Bézout: s·f^m + t·cofactor = 1.
        let (gcd, s, t) = fm.extended_gcd_poly(&cofactor);
        if gcd.degree() != Some(0) {
            return None;
        }
        let c = gcd.coeffs().first().cloned()?;
        let inv_c = RationalDomain.inv(&c)?;
        let s = s.mul_scalar(&inv_c);
        let t = t.mul_scalar(&inv_c);
        // num/den = num·t/f^m + num·s/cofactor.
        let nt = current.mul(&t);
        let (quotient, proper) = nt.div_rem(&fm)?;
        poly = poly.add(&quotient);
        push_power_pieces(f, *mult, &proper, &mut pieces)?;
        current = current.mul(&s);
        remaining = cofactor;
    }
    if !current.is_zero() {
        return None;
    }
    Some((poly, pieces))
}

/// Split `proper/f^mult` (with `deg proper < mult·deg f`) into `A_j/f^j`.
fn push_power_pieces(
    f: &DPoly,
    mult: usize,
    proper: &DPoly,
    pieces: &mut Vec<(DPoly, DPoly)>,
) -> Option<()> {
    let mut rest = proper.clone();
    for j in (1..=mult).rev() {
        let (quotient, rem) = rest.div_rem(f)?;
        if !rem.is_zero() {
            pieces.push((rem, f.pow(j as u32)));
        }
        rest = quotient;
    }
    if !rest.is_zero() {
        return None;
    }
    Some(())
}

/// ∫ (n₁·w² + n₀)/(w⁴ + p·w² + q) dw — the even-quartic gap left by the ℚ
/// engine (which cannot factor `w⁴ + 1` and has no algebraic extensions).
///
/// - Case A (rational roots `r₁,₂`, discriminant a perfect square):
///   `A/(w² + r₁) + B/(w² + r₂)`, each an `atan`/`log` form.
/// - Case B (no real roots, `p² < 4q`): the classical
///   `ℚ(√q, √(2√q − p))` split
///   `w⁴ + p·w² + q = (w² − e·w + f)(w² + e·w + f)` with `f = √q`,
///   `e = √(2√q − p)`, whose answer needs only `log` and `atan` forms over
///   the radical atoms `f`, `e` and `√Δ` (`Δ = p + 2√q`).
fn even_quartic_atom<'a>(
    ctx: &'a AtomArena<'a>,
    num: &DPoly,
    den: &DPoly,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let zero = RationalDomain.zero();
    if num.degree().unwrap_or(0) > 2 || den.degree() != Some(4) {
        return None;
    }
    // Only the monic normal form is handled (the caller makes the
    // denominator monic; a scaled form would need the coefficients scaled
    // back before the closed form below).
    let lc = den.coeffs().get(4).cloned().unwrap_or_else(|| zero.clone());
    if lc != Rational::new(1, 1) {
        return None;
    }
    let p = den.coeffs().get(2).cloned().unwrap_or_else(|| zero.clone());
    let q = den
        .coeffs()
        .first()
        .cloned()
        .unwrap_or_else(|| zero.clone());
    let n0 = num
        .coeffs()
        .first()
        .cloned()
        .unwrap_or_else(|| zero.clone());
    let n1 = num.coeffs().get(2).cloned().unwrap_or_else(|| zero.clone());
    if RationalDomain.is_zero(&n1) && RationalDomain.is_zero(&n0) {
        return None;
    }
    let four = Rational::new(4, 1);
    let disc = RationalDomain.sub(&RationalDomain.mul(&p, &p), &RationalDomain.mul(&four, &q));
    if let Some(root) = sqrt_exact(&disc) {
        // Case A: r = (−p ± √disc)/2 with rational roots.
        let two = Rational::new(2, 1);
        let neg_p = RationalDomain.neg(&p);
        let r1 = RationalDomain.div(&RationalDomain.add(&neg_p, &root), &two)?;
        let r2 = RationalDomain.div(&RationalDomain.sub(&neg_p, &root), &two)?;
        if r1 == r2 {
            return None;
        }
        // A/(w²+r₁) + B/(w²+r₂) with A = (n₀ − n₁·r₁)/(r₂ − r₁), B = n₁ − A.
        let a_coef = RationalDomain.div(
            &RationalDomain.sub(&n0, &RationalDomain.mul(&n1, &r1)),
            &RationalDomain.sub(&r2, &r1),
        )?;
        let b_coef = RationalDomain.sub(&n1, &a_coef);
        let ta = quadratic_reciprocal_atom(ctx, &a_coef, &r1, x)?;
        let tb = quadratic_reciprocal_atom(ctx, &b_coef, &r2, x)?;
        return Some(normalize(ctx, ctx.add(&[ta, tb])));
    }
    // Case B needs disc < 0 (no real roots), which also forces q > 0.
    if !rat_is_neg(&disc) {
        return None;
    }
    let two = ctx.num(2);
    let f = sqrt_pos_atom(ctx, &q)?;
    let p_atom = const_atom(ctx, &p)?;
    let two_f = normalize(ctx, ctx.mul(&[two, f]));
    let e_sq = normalize(
        ctx,
        ctx.add(&[two_f, normalize(ctx, ctx.mul(&[ctx.num(-1), p_atom]))]),
    );
    let e = ctx.pow(e_sq, rat_atom(ctx, 1, 2));
    let delta = normalize(ctx, ctx.add(&[p_atom, two_f]));
    let sq_delta = ctx.pow(delta, rat_atom(ctx, 1, 2));
    let n0_atom = const_atom(ctx, &n0)?;
    let n1_atom = const_atom(ctx, &n1)?;
    // log term: (A/2)·log((w² − e·w + f)/(w² + e·w + f)), A = (n₁·f − n₀)/(2·e·f).
    let a_num = normalize(
        ctx,
        ctx.add(&[
            ctx.mul(&[n1_atom, f]),
            normalize(ctx, ctx.mul(&[ctx.num(-1), n0_atom])),
        ]),
    );
    let a_den = normalize(ctx, ctx.mul(&[two, e, f]));
    let a_coef = normalize(ctx, ctx.mul(&[a_num, inv(ctx, a_den)]));
    let w2 = ctx.pow(x, ctx.num(2));
    let ew = normalize(ctx, ctx.mul(&[e, x]));
    let top = normalize(
        ctx,
        ctx.add(&[w2, normalize(ctx, ctx.mul(&[ctx.num(-1), ew])), f]),
    );
    let bottom = normalize(ctx, ctx.add(&[w2, ew, f]));
    let ratio = normalize(ctx, ctx.mul(&[top, inv(ctx, bottom)]));
    let log_term = normalize(
        ctx,
        ctx.mul(&[a_coef, rat_atom(ctx, 1, 2), ctx.fun("log", &[ratio])]),
    );
    // atan term: ((n₀ + n₁·f)/(2·f·√Δ))·[atan((2w − e)/√Δ) + atan((2w + e)/√Δ)].
    let c_num = normalize(ctx, ctx.add(&[n0_atom, ctx.mul(&[n1_atom, f])]));
    let c_den = normalize(ctx, ctx.mul(&[two, f, sq_delta]));
    let c_coef = normalize(ctx, ctx.mul(&[c_num, inv(ctx, c_den)]));
    let two_w = normalize(ctx, ctx.mul(&[two, x]));
    let shift = normalize(ctx, ctx.mul(&[ctx.num(-1), e]));
    let arg1 = normalize(
        ctx,
        ctx.mul(&[normalize(ctx, ctx.add(&[two_w, shift])), inv(ctx, sq_delta)]),
    );
    let arg2 = normalize(
        ctx,
        ctx.mul(&[normalize(ctx, ctx.add(&[two_w, e])), inv(ctx, sq_delta)]),
    );
    let atan_term = normalize(
        ctx,
        ctx.mul(&[
            c_coef,
            ctx.add(&[ctx.fun("atan", &[arg1]), ctx.fun("atan", &[arg2])]),
        ]),
    );
    Some(normalize(ctx, ctx.add(&[log_term, atan_term])))
}

/// ∫ c/(w² + r) dw for rational `c` and `r` (`r = 0` gives `−c/w`).
fn quadratic_reciprocal_atom<'a>(
    ctx: &'a AtomArena<'a>,
    c: &Rational,
    r: &Rational,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    if RationalDomain.is_zero(c) {
        return Some(ctx.num(0));
    }
    let c_atom = const_atom(ctx, c)?;
    if RationalDomain.is_zero(r) {
        return Some(normalize(ctx, ctx.mul(&[ctx.num(-1), c_atom, inv(ctx, x)])));
    }
    if rat_is_pos(r) {
        // c/√r·atan(w/√r).
        let root = sqrt_pos_atom(ctx, r)?;
        let arg = normalize(ctx, ctx.mul(&[x, inv(ctx, root)]));
        return Some(normalize(
            ctx,
            ctx.mul(&[c_atom, inv(ctx, root), ctx.fun("atan", &[arg])]),
        ));
    }
    // c/(2√(−r))·log((w − √(−r))/(w + √(−r))).
    let minus_r = RationalDomain.neg(r);
    let root = sqrt_pos_atom(ctx, &minus_r)?;
    let top = normalize(
        ctx,
        ctx.add(&[x, normalize(ctx, ctx.mul(&[ctx.num(-1), root]))]),
    );
    let bottom = normalize(ctx, ctx.add(&[x, root]));
    let ratio = normalize(ctx, ctx.mul(&[top, inv(ctx, bottom)]));
    let log = ctx.fun("log", &[ratio]);
    Some(normalize(
        ctx,
        ctx.mul(&[c_atom, inv(ctx, ctx.num(2)), inv(ctx, root), log]),
    ))
}

/// The atom for a rational constant.
fn const_atom<'a>(ctx: &'a AtomArena<'a>, r: &Rational) -> Option<Atom<'a>> {
    Some(rat_atom(ctx, r.numer().to_i64()?, r.denom().to_i64()?))
}

/// Exact rational square root, `None` when `r < 0` or irrational.
fn sqrt_exact(r: &Rational) -> Option<Rational> {
    if rat_is_neg(r) {
        return None;
    }
    let p = r.numer().to_i64()?;
    let q = r.denom().to_i64()?;
    let n = (p as i128).checked_mul(q as i128)?;
    let m = n.checked_isqrt()?;
    if m * m != n {
        return None;
    }
    Some(Rational::new(i64::try_from(m).ok()?, q))
}

/// `√r` for `r ≥ 0` as an atom: exact when the root is rational, otherwise
/// the radical `√(p·q)/q`.
fn sqrt_pos_atom<'a>(ctx: &'a AtomArena<'a>, r: &Rational) -> Option<Atom<'a>> {
    if let Some(exact) = sqrt_exact(r) {
        return const_atom(ctx, &exact);
    }
    if rat_is_neg(r) {
        return None;
    }
    let p = r.numer().to_i64()?;
    let q = r.denom().to_i64()?;
    if p == 0 {
        return Some(ctx.num(0));
    }
    let n = p.checked_mul(q)?;
    let root = ctx.pow(ctx.num(n), rat_atom(ctx, 1, 2));
    Some(normalize(ctx, ctx.mul(&[root, inv(ctx, ctx.num(q))])))
}

fn rat_is_pos(r: &Rational) -> bool {
    use num_traits::Signed;
    r.inner().is_positive()
}

fn rat_is_neg(r: &Rational) -> bool {
    use num_traits::Signed;
    r.inner().is_negative()
}

/// The atom for the dense polynomial `p` in the variable `x`.
fn poly_atom<'a>(ctx: &'a AtomArena<'a>, p: &DPoly, x: Atom<'a>) -> Option<Atom<'a>> {
    if p.is_zero() {
        return Some(ctx.num(0));
    }
    let field = GeneratorField::from_polynomial(dense_to_sparse(p));
    rational_to_atom(ctx, &field, &[x])
}

/// The atom for the rational function `num/den` in the variable `x`.
fn ratio_atom<'a>(
    ctx: &'a AtomArena<'a>,
    num: &DPoly,
    den: &DPoly,
    x: Atom<'a>,
) -> Option<Atom<'a>> {
    let field = GeneratorField::from_num_den(dense_to_sparse(num), dense_to_sparse(den));
    rational_to_atom(ctx, &field, &[x])
}

fn sparse_to_dense(p: &Sparse) -> Option<DPoly> {
    if p.n_vars() != 1 {
        return None;
    }
    let deg = p.degree_in(0);
    if deg > MAX_RAT_NUM_DEG {
        return None;
    }
    let mut coeffs = vec![RationalDomain.zero(); deg + 1];
    for (exp, coeff) in p.terms_ref() {
        let index = *exp.first()?;
        *coeffs.get_mut(index)? = coeff.clone();
    }
    Some(DPoly::from_coeffs(RationalDomain, coeffs))
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
    use super::*;
    use ocas_core::arena::Arena;

    /// Sample points where every kernel family is positive and smooth
    /// (`tan x`, `cot x`, `sinh x`, `cosh x` all `> 0`).
    const SAMPLE_POS: [f64; 3] = [0.35, 0.75, 1.15];
    /// Sample points that keep a symbolic linear argument `c + d·x` inside
    /// `(0, π/2)` for the parameter values used by the symbolic tests.
    const SAMPLE_SMALL: [f64; 3] = [0.4, 0.8, 1.2];

    /// Numeric f64 evaluator for test verification; `log` goes through the
    /// absolute value so partial-fraction log signs do not matter.
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
                    "cot" => 1.0 / v.tan(),
                    "sec" => 1.0 / v.cos(),
                    "csc" => 1.0 / v.sin(),
                    "exp" => v.exp(),
                    "log" => v.abs().ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    "asin" => v.asin(),
                    "acos" => v.acos(),
                    "sinh" => v.sinh(),
                    "cosh" => v.cosh(),
                    "tanh" => v.tanh(),
                    "coth" => 1.0 / v.tanh(),
                    "sech" => 1.0 / v.cosh(),
                    "csch" => 1.0 / v.sinh(),
                    "asinh" => v.asinh(),
                    "acosh" => v.acosh(),
                    "atanh" => v.atanh(),
                    _ => return None,
                })
            }
        }
    }

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, s).unwrap()
    }

    /// Run the mechanism directly, require `Some` without residue, and check
    /// `diff(result) == integrand` numerically at every sample point.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_kernel_subst(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined: {integrand}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64(d, &env).unwrap_or_else(|| {
                panic!(
                    "eval diff failed\n  integrand: {integrand}\n  result: {result}\n  diff: {d}"
                )
            });
            let rhs = eval_f64(integrand, &env)
                .unwrap_or_else(|| panic!("eval integrand failed: {integrand}"));
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "integrand: {integrand}\nat x={xv}: diff={lhs} integrand={rhs}\nresult: {result}"
            );
        }
    }

    /// The Pythagorean rewrite fires on `1 − tanh²`, and the family-preserving
    /// selection then falls back to the original in-family form, which
    /// substitutes to `(1 − t²)^(3/2)` and integrates via the binomial engine.
    #[test]
    fn rewritten_integrand_stays_integrated() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "(1 - tanh(x)^2)^(3/2)"),
            var,
            &[],
            &SAMPLE_POS,
        );
    }

    #[test]
    fn numeric_tan_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for s in [
            "sqrt(tan(x))",
            "1/(2 + 3*tan(x))",
            "1/(5 + 3*tan(x))^3",
            "tan(x)^3/(1 + tan(x)^2)^(3/2)",
            "tan(x)^2/sqrt(1 + tan(x)^2)",
            "tan(x)^4",
            "1/(1 + tan(x)^2)^2",
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, s), var, &[], &SAMPLE_POS);
        }
    }

    #[test]
    fn numeric_cot_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for s in [
            "sqrt(cot(x))",
            "1/(2 + 3*cot(x))",
            "(1 + cot(x)^2)^(3/2)",
            "(1 + cot(x)^2)^(5/2)",
            "cot(x)^2*(1 + cot(x)^2)^(3/2)",
            "cot(x)^4/(1 + cot(x)^2)",
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, s), var, &[], &SAMPLE_POS);
        }
    }

    #[test]
    fn numeric_tanh_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for s in [
            "(1 - tanh(x)^2)^(3/2)",
            "(1 - tanh(x)^2)^(5/2)",
            "tanh(x)^3",
            "1/(2 + 3*tanh(x))",
            "tanh(x)^2/(1 + tanh(x)^2)",
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, s), var, &[], &SAMPLE_POS);
        }
    }

    #[test]
    fn numeric_coth_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for s in [
            "(1 + coth(x))^(7/2)",
            "sqrt(coth(x))",
            "1/(2 + 3*coth(x))",
            "coth(x)/(1 + coth(x)^2)",
            "coth(x)^3/(1 + coth(x)^2)",
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, s), var, &[], &SAMPLE_POS);
        }
    }

    // ---------------------------------------------------------------------
    // (b) symbolic-coefficient corpus shapes
    // ---------------------------------------------------------------------

    #[test]
    fn symbolic_coefficient_corpus_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let env = [
            (Symbol::new("a"), 1.3),
            (Symbol::new("b"), 0.7),
            (Symbol::new("c"), 0.15),
            (Symbol::new("d"), 0.8),
            (Symbol::new("e"), 0.1),
            (Symbol::new("f"), 0.7),
            // The imaginary unit is carried as a free symbol; evaluating it
            // at 1 still checks the differentiation identity.
            (Symbol::new("i"), 1.0),
        ];
        for s in [
            "coth(c + d*x)^4*(a + b*tanh(c + d*x)^2)^2",
            "(1 + coth(c + d*x))^(7/2)",
            "1/(5 + 3*tan(c + d*x))^3",
            "(a + i*a*tan(e + f*x))^3/(d*tan(e + f*x))^(7/2)",
            "sqrt(d*tan(e + f*x))*(a + i*a*tan(e + f*x))^2",
            "a*tan(b*x)^3",
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, s), var, &env, &SAMPLE_SMALL);
        }
    }

    // ---------------------------------------------------------------------
    // (c) honest declines
    // ---------------------------------------------------------------------

    #[test]
    fn declines_outside_the_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for s in [
            // Not a rational-derivative kernel.
            "sin(x)",
            "sqrt(sin(x))",
            "cos(tanh(a + b*x))",
            // `sec`/`csc` have irrational derivatives in the kernel
            // (`sec' = sec·tan`), so substitution would introduce a second
            // radical; they are left to the rule table and trig reductions.
            "1/(a + b*sec(x))",
            "csc(3*x + 1)^3",
            "sec(2*x + 1)^4",
            "sech(x)^3",
            // Rational-derivative kernels mixed with an out-of-family head.
            "(a*sin(e + f*x))^(3/2)*sqrt(b*tan(e + f*x))",
            // Kernels paired with an unrelated exponential.
            "tan(x)*exp(x)",
            // Nonlinear kernel argument (`c + d·x^(1/3)`).
            "1/(a + b*tan(c + d*x^(1/3)))^2",
            // Symbolic exponent on a kernel-dependent base.
            "(b*tan(c + d*x)^2)^n",
            // Mixed trig and hyperbolic kernels.
            "tan(x)*tanh(x)",
            // A kernel-free algebraic integrand.
            "1/((2+3*t)*(1+t^2))",
            // Shapes whose only engine answers are branch-limited (they hold
            // on `t > 0` but not on the whole real line) and for which no
            // globally valid alternative exists are declined by the
            // certainty gate rather than emitted as a wrong answer.
            "coth(x)^2/(1 + coth(x)^2)^(3/2)",
        ] {
            assert!(
                integrate_kernel_subst(&ctx, parse(&ctx, s), var).is_none(),
                "expected a decline: {s}"
            );
        }
    }

    // ---------------------------------------------------------------------
    // Pythagorean rewrites
    // ---------------------------------------------------------------------

    #[test]
    fn pythagorean_rewrites_are_exact_and_idempotent() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for (input, expected) in [
            ("1 - tanh(x)^2", "sech(x)^2"),
            ("1 - sin(x)^2", "cos(x)^2"),
            ("sec(x)^2 - 1", "tan(x)^2"),
            ("1 + cot(x)^2", "csc(x)^2"),
            ("cosh(x)^2 - sinh(x)^2", "1"),
            ("1 + tan(x)^2", "sec(x)^2"),
            ("a - a*tanh(x)^2", "a*sech(x)^2"),
            ("a - a*sin(x)^2", "a*cos(x)^2"),
        ] {
            let expr = parse(&ctx, input);
            let out = pythagorean_rewrite(&ctx, expr, var)
                .unwrap_or_else(|| panic!("no rewrite for {input}"));
            let want = normalize(&ctx, parse(&ctx, expected));
            assert_eq!(out, want, "{input} → {out}, want {want}");
            // Idempotent: a second pass finds nothing.
            assert!(
                pythagorean_rewrite(&ctx, out, var).is_none(),
                "not idempotent: {input}"
            );
        }
        // Nothing to rewrite: no change is reported.
        for s in [
            "1 + 2*tanh(x)^2",
            "tan(x)",
            "1 + x",
            "2 - tanh(x)^2",
            "1 - tanh(x)",
            "1 + sin(x)^2 - 2*sin(x)^2",
        ] {
            assert!(
                pythagorean_rewrite(&ctx, parse(&ctx, s), var).is_none(),
                "unexpected rewrite: {s}"
            );
        }
    }

    // ---------------------------------------------------------------------
    // (d) stress: repeated runs stay stable and bounded
    // ---------------------------------------------------------------------

    #[test]
    fn stress_repeated_runs_are_stable_and_bounded() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let cases = [
            "sqrt(tan(x))",
            "1/(2 + 3*tan(x))",
            "(1 + coth(x))^(7/2)",
            "(1 - tanh(x)^2)^(3/2)",
            "tanh(x)^3",
            "coth(c + d*x)^4*(a + b*tanh(c + d*x)^2)^2",
            "(a + i*a*tan(e + f*x))^3/(d*tan(e + f*x))^(7/2)",
            "sin(x)",
            "tan(x)*exp(x)",
            "1/(a + b*sec(x))",
            "csc(3*x + 1)^3",
            "(1 + cot(x)^2)^(5/2)",
        ];
        let mut baseline: Vec<bool> = vec![false; cases.len()];
        // 25 rounds × 12 cases = 300 mechanism runs.
        for round in 0..25 {
            for (index, s) in cases.iter().enumerate() {
                let expr = parse(&ctx, s);
                let outcome = integrate_kernel_subst(&ctx, expr, var);
                let solved = outcome.is_some();
                if let Some(r) = outcome {
                    assert!(r.to_string().len() <= 8000, "answer grew: {s}");
                    assert!(
                        node_count(r) <= 1200,
                        "answer node count grew: {s} ({})",
                        node_count(r)
                    );
                    // Stable and self-consistent: the emitted answer always
                    // passes the certainty gate it was accepted under.
                    assert!(
                        verify_candidate(&ctx, expr, r, var),
                        "unverified answer emitted for {s}"
                    );
                }
                if round == 0 {
                    baseline[index] = solved;
                } else {
                    assert_eq!(baseline[index], solved, "unstable outcome for {s}");
                }
            }
        }
        // The mixed list must contain both solved and declined entries, so
        // the stability check is not vacuous.
        assert!(baseline.iter().any(|s| *s));
        assert!(baseline.iter().any(|s| !*s));
    }
}
