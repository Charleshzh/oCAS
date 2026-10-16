//! Single-radical atomic-layer reduction (0.30.0 Wave 6a) — **shipped
//! unwired**.
//!
//! # Status: experiment asset, not in the pipeline
//!
//! The 0.30.0 Wave 6a experiment implemented and measured this stage; the
//! wiring was reverted because the full-corpus measurement broke hard gates,
//! and the code is kept (the `coupled.rs` precedent: infrastructure with its
//! own tests, deliberately not called) as the starting point for 0.31.
//! Measured on 1892 with the stage wired after `multi_radical`:
//!
//! - **+29 solved** (radical +18, power-binomial +9, inverse-trig-hyper +2),
//!   0 regressions — the mechanism itself works;
//! - **1 wrong answer**, `rubi-00401` (`sqrt(-1 + 1/x^2)/(x*(-1 + x^2)^2)`):
//!   the `x = 1/t` substitution's sign/branch handling is incomplete for a
//!   quotient radicand that carries a branch factor — the emitted answer is
//!   off by a global sign (oracle rel err 2.0). **Fix this first in 0.31**
//!   (the refusal gate below is meant to cover exactly this shape and does
//!   not yet);
//! - **timeouts 12 → 14** (two shapes pushed past the 10 s per-case budget);
//! - **wall clock 489 → 523 s**, over the 498 s wave gate (the
//!   `(a+b·x²)^(4/3)/x` decline spends ~1.07 s in downstream `binomial`;
//!   newly accepted shapes cost 6–33 ms each plus a few ~1 s declines).
//!
//! Re-wiring requires the branch fix plus fresh timeout/wall budgets.
//!
//! Scope: integrands built from **one** radical base `Q^(p/q)` whose radicand
//! is affine or quadratic, in spellings the existing single-base engines
//! (`rules`, [`sqrt_quadratic`](super::sqrt_quadratic),
//! [`binomial`](super::binomial), [`kernel_subst`](super::kernel_subst))
//! decline. No algebraic tower is built (`crate::tower` is untouched) and no
//! new radical relation is introduced: every mechanism is a *change of
//! variable* whose substituted form is handed back to the chain, plus the
//! Wave-4a sheet factor `σ = w·(w²)^(−1/2)` where the branch demands one.
//!
//! - **W6a-1, affine-square slice.** A half-integer-power site whose radicand
//!   is the square of an affine `w` — written `w²` syntactically or as the
//!   expanded trinomial `p² + 2pq + q²` ([`super::affine_square_base`]) — is
//!   rewritten `(w²)^(k/2) → w^k·σ` (`σ² = 1`, `Dσ = 0`; the
//!   [`multi_radical`](super::multi_radical) sheet precedent, MR spec §4.2).
//!   When that consumes *every* radical site the remainder is a rational
//!   function and the case is an ordinary rational integral.
//! - **W6a-2, even substitution.** An integrand that is odd in the variable,
//!   `f(x) = x·G(x²)`, becomes `G(v)/2` with `v = x²` (`x·dx = dv/2`, an
//!   exact identity for every `x ≠ 0`). Every radicand that was a polynomial
//!   in `x²` is a polynomial in `v`, so the substituted form is again a
//!   single-radical shape for the existing engines.
//! - **W6a-3, inverse-radicand substitution.** When a radicand is a
//!   polynomial in `x^(−1)` (e.g. `a + b/x`), `x = 1/t` turns it into a
//!   polynomial in `t`; the substitution introduces **no** branch factor,
//!   because the radicand is evaluated after the change of variable rather
//!   than rewritten through a quotient.
//!
//! Everything is refused outside its fragment: a substituted form that still
//! carries a radical site (W6a-1), a fractional power of the variable alone
//! (W6a-2/W6a-3, which would let `normalize` fold `(x²)^(3/2)` into `x³` —
//! `|x|³ ≠ x³`, the branch error `fold_linear_squares` documents), a
//! non-polynomial radicand, a radicand that is provably not real anywhere in
//! the substituted domain (`a + c·v` with numeric `a, c ≤ 0`: there would be
//! no real abscissa to verify an antiderivative at, and no complex branch is
//! modelled here), a residue, or a budget overflow.
//!
//! **Measured.** The stage runs after every dedicated single-radical engine
//! (`rules`, [`sqrt_quadratic`](super::sqrt_quadratic),
//! [`binomial`](super::binomial),
//! [`multi_radical`](super::multi_radical)) and before the transcendental
//! stages and the general heuristic, so it can only pick up shapes they
//! declined. On the 83 single-radical cases of
//! `ocas-tests/data/radical_attribution_030.csv` it takes the solved count
//! from 4 to 17 (13 unlocks, no regression, no new timeout) at 6–33 ms per
//! accepted case in debug; nine of the thirteen carry a `radical`
//! certificate, and four are verified numerically because the certificate
//! layer does not model their radical (see the tests). Declines cost 0.3 µs
//! to 63 ms except `(a + b·x²)^(4/3)/x`, whose substituted form keeps the
//! downstream `binomial` completion busy for ~1.07 s before it gives up.
//!
//! It declines outright inside a `chain::parts_active()` scope (the Wave-4a
//! parts discipline) and charges a deterministic work budget so a decline
//! stays cheap.

// Deliberately unwired (see the status note above): the only callers are the
// unit tests, so the plain-lib target reports dead code.
#![allow(dead_code)]

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    affine_square_base, chain, contains_integral, contains_symbol, gcd_i64, int_pow, integrate_raw,
    inv, is_constant, is_zero_expr, linear_form, node_count, rat_atom, replace_symbol,
};

/// Node budget for the input integrand and for every derived substituted
/// form. The chain re-entry costs far more than the substitution itself, so
/// the bound is deliberately the same order as `kernel_subst`'s.
const MAX_NODES: usize = 200;
/// Node budget for a substituted form handed back to the chain.
const MAX_REENTRY_NODES: usize = 300;
/// Work units a single entry into this module may charge.
const MAX_SR_WORK: u64 = 100_000;

thread_local! {
    static SR_WORK: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Reset the work budget; called by the public stage entry point.
fn reset_sr_budget() {
    SR_WORK.with(|c| c.set(0));
}

/// Charge `units` against the module budget; `true` when exhausted.
fn charge_sr_work(units: u64) -> bool {
    SR_WORK.with(|c| {
        let v = c.get().saturating_add(units);
        c.set(v);
        v > MAX_SR_WORK
    })
}

/// Reduced rational exponent `(p, q)`, `q > 0`, or `None` for a non-rational
/// exponent atom.
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

/// A radical site `base^(p/q)`, `q ≥ 2`, with a base that depends on `var`.
///
/// `sqrt(u)` is the same site as `u^(1/2)`; a constant base is not a site
/// (it is an ordinary constant of the integrand).
fn site_of(expr: Atom<'_>, var: Symbol) -> Option<(Atom<'_>, i64, i64)> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            let base = args[0];
            if is_constant(base, var) {
                None
            } else {
                Some((base, 1, 2))
            }
        }
        AtomNode::Pow(b, e) => {
            if is_constant(*b, var) {
                return None;
            }
            let (p, q) = rat_exponent(*e)?;
            if q >= 2 { Some((*b, p, q)) } else { None }
        }
        _ => None,
    }
}

/// Collect the distinct (normalized) site bases of `expr`.
fn collect_site_bases<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    out: &mut Vec<Atom<'a>>,
) {
    if let Some((b, _, _)) = site_of(expr, var) {
        let nb = normalize(ctx, b);
        if !out.contains(&nb) {
            out.push(nb);
        }
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => {}
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_site_bases(ctx, *a, var, out);
            }
        }
        AtomNode::Pow(b, e) => {
            collect_site_bases(ctx, *b, var, out);
            collect_site_bases(ctx, *e, var, out);
        }
        AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_site_bases(ctx, *a, var, out);
            }
        }
    }
}

/// Whether `expr` still contains a radical site.
fn has_site(expr: Atom<'_>, var: Symbol) -> bool {
    if site_of(expr, var).is_some() {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| has_site(*a, var)),
        AtomNode::Pow(b, e) => has_site(*b, var) || has_site(*e, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| has_site(*a, var)),
    }
}

/// Whether `expr` is a rational form in `var`: integer powers only, function
/// heads only over constants, and the variable itself only as a bare factor.
fn is_rational_in(expr: Atom<'_>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|a| is_rational_in(*a, var)),
        AtomNode::Pow(b, e) => {
            matches!(e.node(), AtomNode::Num(_))
                && is_rational_in(*b, var)
                && is_rational_in(*e, var)
        }
        AtomNode::Fun(_, args) => args.iter().all(|a| is_constant(*a, var)),
    }
}

/// `w` with `base ≡ w²` — the base written as `w²` directly, or as the
/// expanded trinomial of an affine `w` — where `w` must be affine in `var`.
fn square_root_form<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let w = match base.node() {
        AtomNode::Pow(b, e) if matches!(e.node(), AtomNode::Num(2)) => *b,
        AtomNode::Add(_) => affine_square_base(ctx, base, var)?,
        _ => return None,
    };
    let (slope, _) = linear_form(ctx, w, var)?;
    if is_zero_expr(ctx, slope) {
        return None;
    }
    Some(w)
}

// =========================================================================
// W6a-1: affine-square slice
// =========================================================================

/// Rewrite every site whose base is an affine square as `w^k·σ` (`k` odd) or
/// `w^k` (`k` even), pushing each `σ` onto `sigmas`.
fn slice_squares<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    targets: &[(Atom<'a>, Atom<'a>)],
    sigmas: &mut Vec<Atom<'a>>,
) -> Option<Atom<'a>> {
    if let Some((b, p, q)) = site_of(expr, var) {
        let nb = normalize(ctx, b);
        if q == 2
            && let Some((_, w)) = targets.iter().find(|(base, _)| *base == nb)
        {
            // `(w²)^(p/2) = w^p` for even `p`; odd `p` carries the sheet.
            if p % 2 == 0 {
                return Some(int_pow(ctx, *w, p));
            }
            sigmas.push(ctx.mul(&[*w, ctx.pow(nb, rat_atom(ctx, -1, 2))]));
            return Some(int_pow(ctx, *w, p));
        }
        // A site this pass does not own stays untouched; the "no residue"
        // gate in the caller then declines the whole mechanism.
        return Some(expr);
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(expr),
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_squares(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_squares(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => Some(ctx.pow(slice_squares(ctx, *b, var, targets, sigmas)?, *e)),
        AtomNode::Fun(name, args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_squares(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.fun(name.as_str(), &kids?))
        }
    }
}

/// W6a-1: `√(w²) = w·σ` for an affine `w`, which turns a perfect-square
/// radicand into an ordinary rational integral (times the formally constant
/// sheet factor).
fn try_square_slice<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let mut bases = Vec::new();
    collect_site_bases(ctx, expr, var, &mut bases);
    let targets: Vec<(Atom<'a>, Atom<'a>)> = bases
        .into_iter()
        .filter_map(|b| square_root_form(ctx, b, var).map(|w| (b, w)))
        .collect();
    if targets.is_empty() {
        return None;
    }
    if charge_sr_work(node_count(expr) as u64) {
        return None;
    }
    let mut sigmas: Vec<Atom<'a>> = Vec::new();
    let rest = slice_squares(ctx, expr, var, &targets, &mut sigmas)?;
    if sigmas.is_empty() || has_site(rest, var) || !is_rational_in(rest, var) {
        return None;
    }
    let rest = normalize(ctx, rest);
    if node_count(rest) > MAX_REENTRY_NODES {
        return None;
    }
    let anti = integrate_raw(ctx, rest, var, 0, true, 0, 0);
    if contains_integral(anti) {
        return None;
    }
    sigmas.push(anti);
    Some(normalize(ctx, ctx.mul(&sigmas)))
}

// =========================================================================
// W6a-2: even substitution `v = x²`
// =========================================================================

/// `g` with `f = var·g`, i.e. the integrand with one factor of the variable
/// divided out. `None` unless a decomposition exists on this branch; the
/// caller validates that the remainder really is even in `var`.
fn strip_one_var<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Var(v) if *v == var => Some(ctx.num(1)),
        AtomNode::Pow(b, e) => match (b.node(), e.node()) {
            (AtomNode::Var(v), AtomNode::Num(k)) if *v == var => {
                if *k == 1 {
                    Some(ctx.num(1))
                } else {
                    Some(ctx.pow(*b, ctx.num(*k - 1)))
                }
            }
            // An integer power of a product distributes over its factors
            // (`(x³·√u)^(−1) = x^(−3)·(√u)^(−1)`), which is what exposes the
            // odd power of the variable in the corpus's spelling.
            (AtomNode::Mul(inner), AtomNode::Num(_)) => {
                let kids: Vec<Atom<'a>> = inner.iter().map(|a| ctx.pow(*a, *e)).collect();
                strip_one_var(ctx, ctx.mul(&kids), var)
            }
            _ => None,
        },
        // Exactly one factor is divided; which one is decided by the caller's
        // evenness check, so every candidate index is offered in turn.
        AtomNode::Mul(_) => None,
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> =
                args.iter().map(|a| strip_one_var(ctx, *a, var)).collect();
            Some(ctx.add(&kids?))
        }
        _ => None,
    }
}

/// Rewrite an *even* form `f(var)` as `g(v)` with `v = var²`. `None` when
/// `f` is not even in `var` (an odd power / odd factor survives).
fn to_even_form<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    v: Atom<'a>,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => {
            if matches!(expr.node(), AtomNode::Var(s) if *s == var) {
                return None;
            }
            Some(expr)
        }
        AtomNode::Pow(b, e) => {
            if let (AtomNode::Var(s), AtomNode::Num(k)) = (b.node(), e.node())
                && *s == var
            {
                if k % 2 != 0 {
                    return None;
                }
                return Some(if *k == 0 {
                    ctx.num(1)
                } else {
                    int_pow(ctx, v, *k / 2)
                });
            }
            if contains_symbol(*e, var) {
                return None;
            }
            Some(ctx.pow(to_even_form(ctx, *b, var, v)?, *e))
        }
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> =
                args.iter().map(|a| to_even_form(ctx, *a, var, v)).collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> =
                args.iter().map(|a| to_even_form(ctx, *a, var, v)).collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Fun(name, args) => {
            if name.as_str() == "sqrt" && args.len() == 1 {
                return Some(ctx.fun("sqrt", &[to_even_form(ctx, args[0], var, v)?]));
            }
            if args.iter().all(|a| is_constant(*a, var)) {
                Some(expr)
            } else {
                None
            }
        }
    }
}

/// `g` with `f = var·g` and `g` even in `var`, or `None`.
fn odd_split<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    v: Atom<'a>,
) -> Option<Atom<'a>> {
    if let AtomNode::Mul(args) = expr.node() {
        for i in 0..args.len() {
            let Some(stripped) = strip_one_var(ctx, args[i], var) else {
                continue;
            };
            let mut kids: Vec<Atom<'a>> = args.to_vec();
            kids[i] = stripped;
            let candidate = ctx.mul(&kids);
            if let Some(even) = to_even_form(ctx, candidate, var, v) {
                return Some(even);
            }
        }
        return None;
    }
    let stripped = strip_one_var(ctx, expr, var)?;
    to_even_form(ctx, stripped, var, v)
}

/// Whether the radicand of any radical site of `expr` contains a negative
/// power of `var` — the substituted radicand must be a polynomial for the
/// existing single-radical engines to see their own shape.
fn radicand_has_neg_power(expr: Atom<'_>, var: Symbol) -> bool {
    if let Some((b, _, _)) = site_of(expr, var)
        && has_neg_power(b, var)
    {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| radicand_has_neg_power(*a, var))
        }
        AtomNode::Pow(b, e) => radicand_has_neg_power(*b, var) || radicand_has_neg_power(*e, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| radicand_has_neg_power(*a, var)),
    }
}

/// Whether `expr` contains `var` raised to a negative integer power.
fn has_neg_power(expr: Atom<'_>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if let (AtomNode::Var(s), AtomNode::Num(k)) = (b.node(), e.node())
                && *s == var
                && *k < 0
            {
                return true;
            }
            has_neg_power(*b, var) || has_neg_power(*e, var)
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| has_neg_power(*a, var)),
        AtomNode::Fun(_, args) => args.iter().any(|a| has_neg_power(*a, var)),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Every fractional power of a non-constant base in `expr` has a base that
/// cannot fold into a bare power of `var` after back-substitution.
///
/// This is the branch gate of W6a-2/W6a-3: `normalize` is free to fold
/// `(x²)^(3/2)` into `x³`, which is `|x|³` only for `x ≥ 0`, so such an answer
/// must never be emitted. A base that is affine in `var` with a **non-zero**
/// intercept (`a + c·var`) — or any other non-monomial shape — is safe.
fn fractional_bases_safe<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            let fractional = match rat_exponent(*e) {
                Some((p, q)) => q >= 2 || p.abs() > 1_000,
                None => true,
            };
            if fractional && !is_constant(*b, var) && !safe_power_base(ctx, *b, var) {
                return false;
            }
            fractional_bases_safe(ctx, *b, var) && fractional_bases_safe(ctx, *e, var)
        }
        AtomNode::Fun(name, args) => {
            if name.as_str() == "sqrt"
                && let Some(b) = args.first()
                && !is_constant(*b, var)
                && !safe_power_base(ctx, *b, var)
            {
                return false;
            }
            args.iter().all(|a| fractional_bases_safe(ctx, *a, var))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| fractional_bases_safe(ctx, *a, var))
        }
        AtomNode::Num(_) | AtomNode::Var(_) => true,
    }
}

/// A fractional power's base that cannot collapse to a power of `var`: a bare
/// monomial (`var^k`, `c·var^k`) is unsafe; anything else is safe.
fn safe_power_base<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, var: Symbol) -> bool {
    let bare = match base.node() {
        AtomNode::Var(s) => *s == var,
        AtomNode::Mul(args) => {
            let mut varcount = 0;
            let mut others_ok = true;
            for a in args.iter() {
                match a.node() {
                    AtomNode::Var(s) if *s == var => varcount += 1,
                    AtomNode::Pow(b, e) => {
                        if let (AtomNode::Var(s), AtomNode::Num(_)) = (b.node(), e.node())
                            && *s == var
                        {
                            varcount += 1;
                        } else if !is_constant(*a, var) {
                            others_ok = false;
                        }
                    }
                    _ if is_constant(*a, var) => {}
                    _ => others_ok = false,
                }
            }
            varcount == 1 && others_ok
        }
        _ => false,
    };
    if bare {
        return false;
    }
    if let Some((_, intercept)) = linear_form(ctx, base, var)
        && is_zero_expr(ctx, intercept)
    {
        return false;
    }
    true
}

/// Push an integer power into a product: `(u·v)^k → u^k·v^k` for `k ∈ ℤ`, an
/// exact identity. The corpus spells `1/(x³·√u)` as `(x³·√u)^(−1)`, and the
/// even substitution needs the odd power of the variable as a *factor* to
/// strip.
fn split_int_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            if let (AtomNode::Mul(inner), AtomNode::Num(_)) = (b.node(), e.node()) {
                let kids: Vec<Atom<'a>> = inner
                    .iter()
                    .map(|a| split_int_powers(ctx, ctx.pow(*a, *e)))
                    .collect();
                return ctx.mul(&kids);
            }
            ctx.pow(split_int_powers(ctx, *b), split_int_powers(ctx, *e))
        }
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| split_int_powers(ctx, *a)).collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| split_int_powers(ctx, *a)).collect();
            ctx.mul(&kids)
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| split_int_powers(ctx, *a)).collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// W6a-2: `f(x) = x·G(x²)` → `G(v)/2` with `v = x²`.
fn try_even_substitution<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let v_sym = aux_symbol(expr, var, &["v", "w"])?;
    let v = ctx.var(v_sym.as_str());
    if charge_sr_work(node_count(expr) as u64) {
        return None;
    }
    // `normalize` last: it flattens the nested product the integer-power
    // split introduces, which is what lets the factor search see the factors.
    let expr = normalize(ctx, fold_pow_pow(ctx, split_int_powers(ctx, expr)));
    let even = odd_split(ctx, expr, var, v)?;
    let v_form = normalize(ctx, ctx.mul(&[rat_atom(ctx, 1, 2), even]));
    if !has_site(v_form, v_sym)
        || radicand_has_neg_power(v_form, v_sym)
        || provably_nonpositive_radicand(ctx, v_form, v_sym)
    {
        return None;
    }
    if node_count(v_form) > MAX_REENTRY_NODES {
        return None;
    }
    let anti = integrate_raw(ctx, v_form, v_sym, 0, true, 0, 0);
    if contains_integral(anti) || !fractional_bases_safe(ctx, anti, v_sym) {
        return None;
    }
    let back = replace_symbol(ctx, anti, v_sym, ctx.pow(ctx.var(var.as_str()), ctx.num(2)));
    let back = normalize(ctx, back);
    if node_count(back) > MAX_REENTRY_NODES {
        return None;
    }
    Some(back)
}

// =========================================================================
// W6a-3: inverse-radicand substitution `x = 1/t`
// =========================================================================

/// Whether some radical site's base contains a negative power of `var`.
fn has_inverse_radicand(expr: Atom<'_>, var: Symbol) -> bool {
    if let Some((b, _, _)) = site_of(expr, var)
        && has_neg_power(b, var)
    {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| has_inverse_radicand(*a, var))
        }
        AtomNode::Pow(b, e) => has_inverse_radicand(*b, var) || has_inverse_radicand(*e, var),
        AtomNode::Fun(_, args) => args.iter().any(|a| has_inverse_radicand(*a, var)),
    }
}

/// Fold `(u^m)^n` into `u^(m·n)` for rational `m`, `n` — `normalize` does not
/// collapse a power of a power, so the reciprocal substitution's
/// `(t^(−1))^(−1)` would otherwise survive as a non-polynomial radicand.
fn fold_pow_pow<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Pow(b, e) => {
            let fb = fold_pow_pow(ctx, *b);
            let fe = fold_pow_pow(ctx, *e);
            if let AtomNode::Pow(inner, ie) = fb.node()
                && let (Some((p1, q1)), Some((p2, q2))) = (rat_exponent(*ie), rat_exponent(fe))
                && let (Some(p), Some(q)) = (p1.checked_mul(p2), q1.checked_mul(q2))
            {
                return ctx.pow(*inner, rat_atom(ctx, p, q));
            }
            ctx.pow(fb, fe)
        }
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_pow_pow(ctx, *a)).collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_pow_pow(ctx, *a)).collect();
            ctx.mul(&kids)
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_pow_pow(ctx, *a)).collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// W6a-3: `x = 1/t` when every radicand is a polynomial in `x^(−1)`.
fn try_reciprocal<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    if !has_inverse_radicand(expr, var) {
        return None;
    }
    let t_sym = aux_symbol(expr, var, &["t", "u"])?;
    let t = ctx.var(t_sym.as_str());
    if charge_sr_work(node_count(expr) as u64) {
        return None;
    }
    let subbed = fold_pow_pow(ctx, replace_symbol(ctx, expr, var, inv(ctx, t)));
    // dx = −t^(−2)·dt. Integer powers are pushed into their products so the
    // chain sees plain factors (`t^(−4)·(c+d·t²)^(−1/2)`) rather than a power
    // of a product.
    let t_form = normalize(
        ctx,
        fold_pow_pow(
            ctx,
            split_int_powers(
                ctx,
                ctx.mul(&[subbed, ctx.num(-1), ctx.pow(t, ctx.num(-2))]),
            ),
        ),
    );
    if radicand_has_neg_power(t_form, t_sym)
        || !has_site(t_form, t_sym)
        || provably_nonpositive_radicand(ctx, t_form, t_sym)
    {
        return None;
    }
    if node_count(t_form) > MAX_REENTRY_NODES {
        return None;
    }
    let anti = integrate_raw(ctx, t_form, t_sym, 0, true, 0, 0);
    if contains_integral(anti) {
        return None;
    }
    let back = fold_pow_pow(
        ctx,
        replace_symbol(ctx, anti, t_sym, inv(ctx, ctx.var(var.as_str()))),
    );
    let back = normalize(ctx, back);
    if !fractional_bases_safe(ctx, back, var) || node_count(back) > MAX_REENTRY_NODES {
        return None;
    }
    Some(back)
}

/// Whether some radical site of `expr` has an affine radicand `a + c·var`
/// with **numeric** `a ≤ 0` and `c ≤ 0` — i.e. the radical is provably not
/// real anywhere in the substituted domain (`var ≥ 0`, which is what the even
/// and inverse substitutions produce).
///
/// Such a shape has no real abscissa to verify an antiderivative at, and this
/// fragment deliberately models no complex branch, so it is declined rather
/// than emitted unverifiable. Symbolic coefficients are undecidable here and
/// are left to the downstream engines.
fn provably_nonpositive_radicand<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> bool {
    if let Some((b, _, _)) = site_of(expr, var)
        && let Some((slope, intercept)) = linear_form(ctx, b, var)
        // `linear_form` accumulates its coefficients with `ctx.add`/`ctx.mul`,
        // so they are not necessarily folded (`0 + (−9)`); `rat_of` reads the
        // canonical rational shape, which needs the fold first.
        && let (Some(p), Some(q)) = (
            rat_of(normalize(ctx, intercept)),
            rat_of(normalize(ctx, slope)),
        )
        && p.1 > 0
        && q.1 > 0
        && p.0 <= 0
        && q.0 <= 0
        && (p.0 < 0 || q.0 < 0)
    {
        return true;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => args
            .iter()
            .any(|a| provably_nonpositive_radicand(ctx, *a, var)),
        AtomNode::Pow(b, e) => {
            provably_nonpositive_radicand(ctx, *b, var)
                || provably_nonpositive_radicand(ctx, *e, var)
        }
        AtomNode::Fun(_, args) => args
            .iter()
            .any(|a| provably_nonpositive_radicand(ctx, *a, var)),
    }
}

/// Pick a substitution symbol from `candidates` that collides with neither the
/// integration variable nor any symbol already present in `expr`.
fn aux_symbol(expr: Atom<'_>, var: Symbol, candidates: &[&str]) -> Option<Symbol> {
    for name in candidates {
        let s = Symbol::new(name);
        if s != var && !contains_symbol(expr, s) {
            return Some(s);
        }
    }
    None
}

// =========================================================================
// Entry point
// =========================================================================

/// Integrate `expr` via the single-radical atomic-layer reductions (see the
/// module docs). Returns `None` outside the fragment, inside a
/// `chain::parts_active()` scope, or when a budget is exhausted.
pub(crate) fn integrate_single_radical<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if chain::parts_active() {
        return None;
    }
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    // Cheap structural gate: no radical site at all, no work.
    if !has_site(expr, var) {
        return None;
    }
    reset_sr_budget();
    if charge_sr_work(node_count(expr) as u64) {
        return None;
    }
    if let Some(r) = try_square_slice(ctx, expr, var) {
        return Some(r);
    }
    if let Some(r) = try_even_substitution(ctx, expr, var) {
        return Some(r);
    }
    try_reciprocal(ctx, expr, var)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;
    use std::time::Instant;

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, s).unwrap_or_else(|e| panic!("parse {s}: {e:?}"))
    }

    /// Numeric evaluator for test verification.
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
                    "log" => v.abs().ln(),
                    "sqrt" => v.sqrt(),
                    "atan" => v.atan(),
                    "asin" => v.asin(),
                    "atanh" => v.atanh(),
                    "asinh" => v.asinh(),
                    _ => return None,
                })
            }
        }
    }

    /// Run the mechanism directly, require `Some` without residue, and check
    /// `D(result) == integrand` numerically at every sample point.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_single_radical(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined: {integrand}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64(d, &env)
                .unwrap_or_else(|| panic!("eval diff failed\n  result: {result}\n  diff: {d}"));
            let rhs = eval_f64(integrand, &env)
                .unwrap_or_else(|| panic!("eval integrand failed: {integrand}"));
            let tol = 1e-5 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "integrand: {integrand}\nat x={xv}: diff={lhs} integrand={rhs}\nresult: {result}"
            );
        }
    }

    /// The mechanism must decline `integrand` without an answer.
    ///
    /// The bound is deliberately loose: the per-case wall-clock budget is
    /// measured by the corpus probe (one child process per case, the
    /// `integrate_1892` convention), not by a unit test that shares the
    /// machine with the rest of the suite. Measured in isolation in debug, the
    /// whole decline family costs 0.3 µs–63 ms except
    /// `(a + b·x²)^(4/3)/x`, whose substituted form makes the downstream
    /// `binomial` completion work for 1.07 s before it gives up. The
    /// assertion only catches a pathological hang.
    fn assert_declined<'a>(ctx: &'a AtomArena<'a>, integrand: Atom<'a>, var: Symbol) {
        let start = Instant::now();
        let out = integrate_single_radical(ctx, integrand, var);
        let elapsed = start.elapsed();
        assert!(out.is_none(), "expected a decline, got {:?}", out);
        assert!(
            elapsed.as_secs() < 60,
            "decline took {elapsed:?} (hang): {integrand}"
        );
    }

    /// An accepted family must also carry a machine-checkable certificate for
    /// the shape the stage emitted.
    fn assert_certified<'a>(ctx: &'a AtomArena<'a>, integrand: Atom<'a>, var: Symbol) {
        let result = integrate_single_radical(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined: {integrand}"));
        let normalized = normalize(ctx, integrand);
        if let Err(e) = super::super::certify::certify(ctx, normalized, result, var) {
            panic!("certify declined ({e:?}) for {integrand}\n  result: {result}");
        }
    }

    // -- W6a-1: affine-square slice -------------------------------------

    /// A perfect-square radicand (written expanded, as the corpus writes it)
    /// is exactly a rational integrand times the sheet factor.
    #[test]
    fn expanded_square_radicand_slices() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for src in [
            "(d + e*x)^4*sqrt(a^2 + 2*a*b*x + b^2*x^2)",
            "x^3*(A + B*x)*(a^2 + 2*a*b*x + b^2*x^2)^(3/2)",
            "(A + B*x)*(a^2 + 2*a*b*x + b^2*x^2)^(5/2)/x",
            "(a + b*x)/((d + e*x)^5*sqrt(a^2 + 2*a*b*x + b^2*x^2))",
        ] {
            let integrand = parse(&ctx, src);
            let result = integrate_single_radical(&ctx, integrand, var)
                .unwrap_or_else(|| panic!("declined: {src}"));
            assert!(
                !result.to_string().contains("Integral"),
                "residue: {src} -> {result}"
            );
            assert_certified(&ctx, integrand, var);
        }
    }

    /// The slice answer is numerically the integrand's antiderivative: the
    /// sheet factor `σ = w·(w²)^(−1/2)` evaluates to `sign(w)`, which is
    /// exactly the `√(w²) = |w|` branch of the original radicand.
    #[test]
    fn expanded_square_radicand_is_numerically_correct() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let a = Symbol::new("a");
        let b = Symbol::new("b");
        let d = Symbol::new("d");
        let e = Symbol::new("e");
        // `a + b·x > 0` on the sample points, so `σ = +1` and the check is on
        // the positive branch of the original radicand.
        for (src, extra) in [
            ("(d + e*x)^4*sqrt(a^2 + 2*a*b*x + b^2*x^2)", Vec::new()),
            (
                "x^3*(A + B*x)*(a^2 + 2*a*b*x + b^2*x^2)^(3/2)",
                vec![(Symbol::new("A"), 1.0), (Symbol::new("B"), 1.0)],
            ),
        ] {
            let integrand = parse(&ctx, src);
            let mut env: Vec<(Symbol, f64)> = vec![(a, 2.0), (b, 3.0), (d, 1.0), (e, 1.0)];
            env.extend(extra);
            let result = integrate_single_radical(&ctx, integrand, var)
                .unwrap_or_else(|| panic!("declined: {src}"));
            let dv = crate::diff(&ctx, result, var);
            // Samples on **both** sides of `w = a + b·x = 0`: the sheet factor
            // `σ = w·(w²)^(−1/2)` must keep its `sign(w)`, otherwise the answer
            // is the `√(w²) = w` branch and is wrong for `w < 0`.
            for &xv in &[-1.5_f64, -0.8, 0.2, 0.5, 1.1] {
                let mut env = env.clone();
                env.push((var, xv));
                let lhs = eval_f64(dv, &env).expect("eval diff");
                let rhs = eval_f64(integrand, &env).expect("eval integrand");
                let tol = 1e-5 * rhs.abs().max(1.0);
                assert!(
                    (lhs - rhs).abs() < tol,
                    "{src} at x={xv}: diff={lhs} rhs={rhs}\nresult: {result}"
                );
            }
        }
    }

    // -- W6a-2: even substitution ---------------------------------------

    /// `x/((5 − x²)·√(3 − x²))`: odd in `x` with an even quadratic radicand,
    /// so `v = x²` makes it a single-radical linear form an existing engine
    /// owns.
    #[test]
    fn even_substitution_handles_odd_quadratic_radicals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "x/((5 - x^2)*sqrt(3 - x^2))"),
            var,
            &[],
            &[0.1, 0.4, 0.9],
        );
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "x*sqrt(2 - 3*x^2)/(1 + x^2)"),
            var,
            &[],
            &[0.1, 0.4, 0.7],
        );
    }

    /// The certificate layer's limits are pre-existing and independent of this
    /// module: `∫√(a + b/x) dx = (2/(3b))·(a + b/x)^(3/2)` is elementary and
    /// still declines (`NotInField` — a radicand with a negative power of the
    /// variable is a "deeper composite" the radical layer does not model), and
    /// `asin(rational)` answers decline with `NonZero`. The Wave-6a accepted
    /// cases in those families are therefore verified numerically in the tests
    /// above, which the module's contract allows.
    #[test]
    fn certificate_limitations_are_pre_existing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let a = ctx.var("a");
        let b = ctx.var("b");
        let base = ctx.add(&[a, ctx.mul(&[b, ctx.pow(ctx.var("x"), ctx.num(-1))])]);
        let ig = ctx.fun("sqrt", &[base]);
        let ans = ctx.mul(&[
            rat_atom(&ctx, 2, 3),
            inv(&ctx, b),
            ctx.pow(base, rat_atom(&ctx, 3, 2)),
        ]);
        assert!(
            super::super::certify::certify(&ctx, ig, ans, var).is_err(),
            "the certificate layer unexpectedly models x^(−1) radicands"
        );
    }

    /// Fractional-power radicands that are polynomials in `x²` need no
    /// rationalizing substitution at all once `v = x²` has exposed the shape.
    #[test]
    fn even_substitution_handles_fractional_powers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "x^3*(1 + x^2)^(9/14)"),
            var,
            &[],
            &[0.4, 0.9, 1.6],
        );
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "x*(2 + 3*x^2)^(1/3)"),
            var,
            &[],
            &[0.4, 0.9, 1.6],
        );
        assert_certified(&ctx, parse(&ctx, "x*(2 + 3*x^2)^(1/3)"), var);
        // `(1+x²)^(37/14)` is outside the certificate's radical layer (it
        // models square-root generators), so this case is verified numerically
        // above instead — the module's contract is "certificate *or* numeric
        // derivative re-verification".
    }

    // -- W6a-3: inverse radicand ---------------------------------------

    /// `a + b/x` and `a + b/x²` radicands are polynomials in `x^(−1)`, so
    /// `x = 1/t` is an exact change of variable into a single-radical shape.
    #[test]
    fn inverse_radicand_substitution() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let a = Symbol::new("a");
        let b = Symbol::new("b");
        let consts = [(a, 2.0), (b, 3.0)];
        for src in ["1/(x^5*sqrt(a + b/x))", "sqrt(a + b/x)/x^3"] {
            assert_antiderivative_num(&ctx, parse(&ctx, src), var, &consts, &[0.4, 0.9, 1.7]);
        }
        // Both answers keep the radical in the integrand's own base, but the
        // certificate layer does not model an `x^(−1)` radicand at all (see
        // `certificate_limitations_are_pre_existing`), so they are verified
        // numerically here.
    }

    // -- declines --------------------------------------------------------

    #[test]
    fn declines_out_of_fragment_shapes() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        // Cube-root quadratic radicand: no rationalizing substitution.
        assert_declined(&ctx, parse(&ctx, "(b*x + c*x^2)^(2/3)"), var);
        // Mixed parity: the numerator carries an odd linear factor, so the
        // even substitution does not apply.
        assert_declined(
            &ctx,
            parse(&ctx, "(5 - x)*(2 + 3*x^2)^(3/2)/(3 + 2*x)^5"),
            var,
        );
        // Symbolic quadratic radicand over a linear denominator: the sign
        // evidence the existing tables need is undecidable.
        assert_declined(&ctx, parse(&ctx, "sqrt(a + c*x^2)/(d + e*x)"), var);
        // The substituted form needs the field `ℚ(a^(1/3))`, i.e. a tower
        // element this fragment deliberately does not build.
        assert_declined(&ctx, parse(&ctx, "(a + b*x^2)^(4/3)/x"), var);
        // The Chebyshev completion reaches `ℚ(√d)` (the `(w²−d)` partial
        // fractions of `x³·√(a+b/x)`), again outside the fragment.
        assert_declined(&ctx, parse(&ctx, "x^3*sqrt(a + b/x)"), var);
        assert_declined(&ctx, parse(&ctx, "(a + b/x^2)/(x^4*sqrt(c + d/x^2))"), var);
        // A negative power of the variable hidden inside an *integer* power of
        // a rational factor (`(−1+x⁻²)^(−2)`) is not cleared by this fragment:
        // the chain never sees the `t³/(t²−1)^(3/2)` shape the substitution
        // would produce.
        assert_declined(&ctx, parse(&ctx, "sqrt(-1 + 1/x^2)/(x*(-1 + x^2)^2)"), var);
        // `1/(x³·√(−9 − 4x²))` has a radicand that is negative on the whole
        // real axis: no real abscissa to verify an antiderivative at, and the
        // fragment models no complex branch.
        assert_declined(&ctx, parse(&ctx, "1/(x^3*sqrt(-9 - 4*x^2))"), var);
        // Not a radical shape at all.
        assert_declined(&ctx, parse(&ctx, "1/(1 + x^2)"), var);
    }

    #[test]
    fn declines_inside_parts_scope() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "x/((5 - x^2)*sqrt(3 - x^2))");
        let out = super::super::chain::with_parts_scope(|| {
            integrate_single_radical(&ctx, integrand, var)
        });
        assert!(out.is_none(), "must decline inside a parts scope");
    }

    /// The stage's own answers are idempotent shapes: re-entering the stage on
    /// its output does not produce a new answer, so the chain cannot ping-pong.
    #[test]
    fn stage_output_does_not_re_enter() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "x/((5 - x^2)*sqrt(3 - x^2))");
        let first = integrate_single_radical(&ctx, integrand, var).expect("first pass");
        let second = integrate_single_radical(&ctx, first, var);
        assert!(
            second.is_none(),
            "stage re-entered on its own output: {second:?}"
        );
    }
}
