//! Multi-radical-base reduction (0.30.0 Wave 4a).
//!
//! Implements `docs/planning/MULTI_RADICAL_SPEC_CN.md`: integrands with
//! **two or more distinct radical bases** `Sᵢ^(pᵢ/qᵢ)`, which every
//! single-base engine declines (`sqrt_quadratic` accepts one `√q`,
//! `binomial` one fractional-power base, `halfpower` one trig half-power —
//! see its module docs). The module only performs *atomic-layer* reductions
//! to a rational or single-base substitution form, then hands that form back
//! to the chain through [`integrate_raw`]; a wrong answer is never emitted —
//! every failure mode (precondition, non-rational t-form, `Integral`
//! residue, budget) returns `None` so the chain continues.
//!
//! # Reductions
//!
//! - **MR-L (Möbius rationalization of two affine bases).** For
//!   `u = a₀ + a₁·x`, `v = c₀ + c₁·x` with `Δ = a₀c₁ − a₁c₀ ≢ 0` and all
//!   exponent denominators dividing `L ∈ {2, 4}`, the substitution
//!   `t^L = u/v` gives
//!   `x = (a₀ − c₀·t^L)/(c₁·t^L − a₁)`, `v = Δ/(c₁·t^L − a₁)`,
//!   `dx/dt = −L·Δ·t^(L−1)/(c₁·t^L − a₁)²`, and `u = t^L·v`, so every site
//!   `u^(p/q) → t^(pL/q)·v^(p/q)` and the integrand becomes a rational
//!   function of `t` exactly when the `v`-powers combine to integer
//!   exponents (checked structurally after a like-power merge). The rational
//!   t-form is reintegrated by the chain and the answer back-substituted
//!   with `t = (u/v)^(1/L)`. `Δ ≢ 0` is decided after a `collect_terms`
//!   fold (symbolic coefficients); an undecidable or vanishing `Δ` declines
//!   (proportional bases are a single-base shape in disguise).
//! - **Same-radicand merge** (`√u·√v = √(u·v)`; spec §4.1). Only for
//!   *product* shapes: both bases appear as direct factors of one product
//!   with the *same* exponent `p/q`, `u·v` expands to a polynomial of
//!   degree ≤ 2 in `x` (the `sqrt_quadratic` acceptance surface), and a
//!   numeric evidence gate confirms `√u·√v = √(u·v)` in principal value at
//!   ≥ 2 sample points of the integrand's real domain. Quotients
//!   (`√u/√v`) are **never** merged — the Möbius `t = √(u/v)` carries them.
//!   With ≥ 3 bases, mergeable pairs are collapsed one per round
//!   (≤ `MAX_MR_ROUNDS` rounds, same-slope pairs first), reducing to the
//!   2-base cases; a round with no mergeable pair declines.
//!
//! # Fragment boundary (what the collector rejects)
//!
//! The integrand must be rational in `x` and in the radical sites: any
//! `var`-dependent function head other than `sqrt` (`asin`/`log`/`exp`/
//! `sin`/…) — whether as a factor or **inside a radical base** — declines
//! the whole module (the W5 interface contract, INV_SUBST_SPEC §2.7), as do
//! symbolic exponents on `var`-dependent bases, more than `MAX_MR_BASES`
//! distinct bases, exponent denominators above `MAX_MR_Q`, and bases nested
//! inside other bases beyond the denesting depth (slice MR-I). Constant
//! radical coefficients (`sqrt(3)`-style, no `var`) are inert and allowed.
//!
//! # Branch semantics (spec §4)
//!
//! All radicals are principal-value. The Möbius substitution
//! `t = (u/v)^(1/L)` *defines* the branch it integrates on: the answer is a
//! valid formal antiderivative on every interval with `u/v > 0`, and the
//! symbolic certificate (`certify.rs` radical layer, one generator
//! `rᵢ` per distinct base with `rᵢ^qᵢ = Sᵢ`) certifies the emitted form
//! independently of intervals. Merges additionally carry the numeric
//! evidence gate above; the certificate never sees the merge step (it
//! certifies `D(F) − f ≡ 0` in the original bases) — the two-gate split
//! (numeric evidence, symbolic certificate) of `kernel_subst` (0.27.2).
//! No `I` is ever introduced (`tower/build.rs:124` invariant).
//!
//! # Chain placement and parts discipline
//!
//! Wired between `binomial` and `exp_log` (spec §3.1): single-base engines
//! keep first claim, and the stage must precede `heuristic`, whose parts
//! recursion would burn the chain budget on these products. Inside a
//! `chain::parts_active()` scope the module declines outright — its
//! reductions either succeed completely or not at all (spec §1.3-7,
//! the `rubi-01646` lesson).
//!
//! # Budgets (deterministic; spec §3.2)
//!
//! `MAX_MR_NODES` (input nodes), `MAX_MR_BASES`, `MAX_MR_Q`,
//! `MAX_MR_SUBST_NODES` (substituted forms), `MAX_MR_NEST` (denesting),
//! `MAX_MR_ROUNDS` (pair-merging rounds) and `MAX_MR_WORK` (stage total,
//! same accounting model as `sqrt_quadratic`'s `MAX_SQRT_WORK`). The
//! Möbius t-power guard `MAX_MR_T_EXP` keeps a large literal exponent
//! (`x^1000`-style, one atom but a huge polynomial) inside the substituted
//! node budget's intent; the same role is played for site exponents by
//! `MAX_MR_SITE_P`.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    chain, contains_integral, gcd_i64, int_pow, integrate_raw, inv, is_constant, lcm_i64,
    linear_form, node_count, pick_subst_symbol, rat_atom, replace_symbol,
};

/// Node budget for the input integrand.
const MAX_MR_NODES: usize = 200;
/// Cap on the number of distinct radical bases.
const MAX_MR_BASES: usize = 4;
/// Cap on exponent denominators (`L = lcm ≤ 4` restricts the Möbius lift).
const MAX_MR_Q: i64 = 4;
/// Node budget for substituted t-forms (aligned with `sqrt_quadratic`'s
/// `MAX_REENTRY_NODES`).
const MAX_MR_SUBST_NODES: usize = 400;
/// Cap on nested-radical denesting rounds.
const MAX_MR_NEST: usize = 2;
/// Cap on pair-merging rounds for ≥ 3 bases.
const MAX_MR_ROUNDS: usize = 3;
/// Stage work budget (same accounting model as `sqrt_quadratic`'s
/// `MAX_SQRT_WORK`): every full-tree scan/substitution charges its node
/// count.
const MAX_MR_WORK: u64 = 100_000;
/// Cap on `|p·L/q|`, the t-power a Möbius site maps to (keeps the t-form's
/// *value* size inside the node budget's intent — a one-atom `t^1000` is
/// otherwise invisible to it).
const MAX_MR_T_EXP: i64 = 16;
/// Cap on the numerator magnitude of a site's fractional exponent.
const MAX_MR_SITE_P: i64 = 16;
/// Cap on the exponent produced by the like-power merge pass.
const MAX_MR_MERGE_EXP: i64 = 64;

thread_local! {
    static MR_WORK: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Reset the stage work budget; called once per module entry.
fn reset_mr_budget() {
    MR_WORK.with(|c| c.set(0));
}

/// Charge `units` against the stage budget; `true` when exhausted.
fn charge_mr_work(units: u64) -> bool {
    MR_WORK.with(|c| {
        let v = c.get().saturating_add(units);
        c.set(v);
        v > MAX_MR_WORK
    })
}

/// `(−1)·a`, normalized so numeric atoms fold.
fn neg<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    normalize(ctx, ctx.mul(&[ctx.num(-1), a]))
}

/// `mul` of a slice with the empty/singleton cases folded; normalized.
/// Nested `Mul` factors are spliced in (`ctx.mul` does not flatten, and a
/// nested product is invisible to the like-power merge's grouping).
fn mul_atoms<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>]) -> Atom<'a> {
    let mut flat: Vec<Atom<'a>> = Vec::with_capacity(factors.len());
    for f in factors {
        if let AtomNode::Mul(args) = f.node() {
            flat.extend_from_slice(args);
        } else {
            flat.push(*f);
        }
    }
    match flat.len() {
        0 => ctx.num(1),
        1 => flat[0],
        _ => normalize(ctx, ctx.mul(&flat)),
    }
}

/// Whether `e` folds to the literal zero under `collect_terms` (the
/// symbolic-coefficient zero test of spec §8-5).
fn identically_zero<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> bool {
    matches!(
        crate::ode::util::collect_terms(ctx, e).node(),
        AtomNode::Num(0)
    )
}

/// A reduced rational exponent `(p, q)` with `q > 0`, or `None` when the
/// atom is not a rational number. `rat_of` does not reduce (`2·2⁻¹` reads
/// as `(2, 2)`); reduction happens here.
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

/// Whether `expr` contains a radical site (a non-constant `sqrt` argument or
/// fractional power).
fn contains_radical(expr: Atom<'_>, var: Symbol) -> bool {
    site_of(expr, var).is_some()
        || match expr.node() {
            AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
                args.iter().any(|a| contains_radical(*a, var))
            }
            AtomNode::Pow(b, e) => contains_radical(*b, var) || contains_radical(*e, var),
            AtomNode::Num(_) | AtomNode::Var(_) => false,
        }
}

/// Read `atom` as a genuine radical site `base^(p/q)` with `q ≥ 2`,
/// `base` non-constant in `var`. Handles the three spellings: `sqrt(S)`,
/// `S^(p/q)`, and `sqrt(S)^n` with odd `n` (even `n` is exactly `S^(n/2)`,
/// no radical). Returns the *unnormalized* base; callers normalize.
fn site_of(atom: Atom<'_>, var: Symbol) -> Option<(Atom<'_>, i64, i64)> {
    match atom.node() {
        AtomNode::Fun(name, args)
            if name.as_str() == "sqrt" && args.len() == 1 && !is_constant(args[0], var) =>
        {
            Some((args[0], 1, 2))
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && !is_constant(args[0], var)
            {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                return if n % 2 != 0 {
                    Some((args[0], *n, 2))
                } else {
                    None
                };
            }
            let (p, q) = rat_exponent(*e)?;
            if q < 2 || is_constant(*b, var) {
                return None;
            }
            Some((*b, p, q))
        }
        _ => None,
    }
}

/// Read a product factor as `base^(p/q)` for the like-power merge: any
/// non-constant base, integer exponents included (`q = 1`); a bare
/// non-constant factor reads as `(f, 1, 1)`. The base is normalized so
/// equal bases hash to the same atom.
fn radical_factor<'a>(
    ctx: &'a AtomArena<'a>,
    f: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, i64, i64)> {
    if is_constant(f, var) {
        return None;
    }
    match f.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            Some((normalize(ctx, args[0]), 1, 2))
        }
        AtomNode::Pow(b, e) => {
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && !is_constant(args[0], var)
                && let AtomNode::Num(n) = e.node()
            {
                // Even powers of a square root are exact integer powers of
                // the radicand.
                return if n % 2 == 0 {
                    Some((normalize(ctx, args[0]), n / 2, 1))
                } else {
                    Some((normalize(ctx, args[0]), *n, 2))
                };
            }
            let (p, q) = rat_exponent(*e)?;
            Some((normalize(ctx, *b), p, q))
        }
        // A bare non-constant factor is `f^1`; normalize the key so it
        // groups with `Pow(f, …)` spellings (bases built by
        // `collect_terms` are not necessarily `normalize`-canonical).
        _ => Some((normalize(ctx, f), 1, 1)),
    }
}

// =========================================================================
// P1: like-power merge (same-base fractional exponents are summed)
// =========================================================================

/// `base^k` for integer `k`, folding a numeric inner power
/// (`(g^m)^k → g^(m·k)`, exact for integers): every occurrence of a base
/// must stay the *same atom* — the like-power merge keys on it.
fn pow_i<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, k: i64) -> Atom<'a> {
    if let AtomNode::Pow(g, e) = base.node()
        && let AtomNode::Num(m) = e.node()
        && let Some(mk) = m.checked_mul(k)
    {
        return int_pow(ctx, *g, mk);
    }
    int_pow(ctx, base, k)
}

/// Distribute integer powers over product bases (`(N·D⁻¹)² → N²·D⁻²`),
/// bottom-up. Exact for integer exponents; this flattens the nested
/// quotient powers a substitution builds (`x_t = N·D⁻¹` sitting under the
/// integrand's `x^k`) into the plain numerator/denominator powers the
/// rational backends expect — measured: the nested spelling stalls the
/// symbolic-coefficient rational stage on the Möbius t-form of
/// rubi-00293 (`x²·√(a+b·x)/(c+d·x)^(3/2)`).
fn distribute_int_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| distribute_int_powers(ctx, *a))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| distribute_int_powers(ctx, *a))
                .collect();
            ctx.mul(&kids)
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| distribute_int_powers(ctx, *a))
                .collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Pow(b, e) => {
            let nb = distribute_int_powers(ctx, *b);
            if let AtomNode::Num(k) = e.node()
                && k.abs() <= MAX_MR_T_EXP
            {
                // Fold nested powers first: `(g^m)^k → g^(m·k)`.
                if let AtomNode::Pow(g, e2) = nb.node()
                    && let AtomNode::Num(m) = e2.node()
                    && let Some(mk) = m.checked_mul(*k)
                {
                    return int_pow(ctx, *g, mk);
                }
                if let AtomNode::Mul(args) = nb.node() {
                    let parts: Vec<Atom<'a>> = args.iter().map(|a| pow_i(ctx, *a, *k)).collect();
                    return mul_atoms(ctx, &parts);
                }
            }
            if nb == *b { expr } else { ctx.pow(nb, *e) }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Merge equal-base power factors throughout `expr` (`u^(1/2)·u^(3/2) → u²`,
/// `sqrt(u)·u^(-1/2) → 1`), bottom-up. Sums of powers of one base are exact
/// in the radical field, so no gate is needed. Unmerged factors keep their
/// original spelling.
fn merge_radical_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| merge_radical_powers(ctx, *a, var))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| merge_radical_powers(ctx, *a, var))
                .collect();
            merge_mul(ctx, &kids, var)
        }
        // Recurse into power bases and function arguments: a like-power pair
        // under a denominator (`(t·V^(1/2)·V^(1/2))^-1`) must merge too.
        // Nested integer powers fold (`(g^m)^n → g^(m·n)`, exact), keeping
        // every occurrence of a base the same atom.
        AtomNode::Pow(b, e) => {
            let nb = merge_radical_powers(ctx, *b, var);
            let ne = merge_radical_powers(ctx, *e, var);
            if let AtomNode::Pow(g, e2) = nb.node()
                && let (AtomNode::Num(m), AtomNode::Num(n)) = (e2.node(), ne.node())
                && let Some(mn) = m.checked_mul(*n)
            {
                return int_pow(ctx, *g, mn);
            }
            if nb == *b && ne == *e {
                expr
            } else {
                ctx.pow(nb, ne)
            }
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| merge_radical_powers(ctx, *a, var))
                .collect();
            let rebuilt = ctx.fun(name.as_str(), &kids);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Merge the factors of one product (see [`merge_radical_powers`]).
/// Nested `Mul` factors are spliced in first (`ctx.mul` does not flatten).
fn merge_mul<'a>(ctx: &'a AtomArena<'a>, factors: &[Atom<'a>], var: Symbol) -> Atom<'a> {
    // (base, p, q, original factors of the group)
    let mut groups: Vec<(Atom<'a>, i64, i64, Vec<Atom<'a>>)> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut flat: Vec<Atom<'a>> = Vec::with_capacity(factors.len());
    for f in factors {
        if let AtomNode::Mul(args) = f.node() {
            flat.extend_from_slice(args);
            continue;
        }
        // `(A·B)^n` with integer `n` (normalize's spelling of a negative
        // power of a product is `(A³·B³)⁻¹`): split into per-factor powers
        // — exact for integers — so the factors join their power groups.
        if let AtomNode::Pow(b, e) = f.node()
            && let AtomNode::Mul(args) = b.node()
            && let AtomNode::Num(n) = e.node()
            && n.abs() <= MAX_MR_MERGE_EXP
        {
            for a in args.iter() {
                flat.push(pow_i(ctx, *a, *n));
            }
            continue;
        }
        flat.push(*f);
    }
    for f in &flat {
        let Some((base, p, q)) = radical_factor(ctx, *f, var) else {
            rest.push(*f);
            continue;
        };
        if let Some(slot) = groups.iter_mut().find(|(b, _, _, _)| *b == base) {
            let num = slot
                .1
                .checked_mul(q)
                .and_then(|v| p.checked_mul(slot.2).and_then(|w| v.checked_add(w)));
            let den = slot.2.checked_mul(q);
            if let (Some(num), Some(den)) = (num, den) {
                let g = gcd_i64(num, den).max(1);
                slot.1 = num / g;
                slot.2 = den / g;
                slot.3.push(*f);
            } else {
                // Arithmetic overflow: keep the factor unmerged.
                rest.push(*f);
            }
            continue;
        }
        groups.push((base, p, q, vec![*f]));
    }
    let mut out: Vec<Atom<'a>> = Vec::with_capacity(rest.len() + groups.len());
    for (base, p, q, orig) in groups {
        if p == 0 {
            continue;
        }
        if p.abs() > MAX_MR_MERGE_EXP || q > MAX_MR_MERGE_EXP {
            out.extend(orig);
            continue;
        }
        if orig.len() == 1 {
            // Unchanged single factor: keep the original spelling.
            out.push(orig[0]);
        } else if q == 1 {
            out.push(int_pow(ctx, base, p));
        } else if (p, q) == (1, 2) {
            out.push(ctx.fun("sqrt", &[base]));
        } else {
            out.push(ctx.pow(base, rat_atom(ctx, p, q)));
        }
    }
    out.extend(rest);
    mul_atoms(ctx, &out)
}

// =========================================================================
// Base collection (P4) and the fragment boundary
// =========================================================================

/// One radical site: `base^(p/q)` with `q ≥ 2`; `base` is normalized.
#[derive(Clone, Copy)]
struct Site<'a> {
    base: Atom<'a>,
    p: i64,
    q: i64,
}

/// The collected radical sites and the distinct (normalized) bases, in
/// first-occurrence order (deterministic).
struct Collect<'a> {
    sites: Vec<Site<'a>>,
    bases: Vec<Atom<'a>>,
}

/// Collect the radical sites of `expr`; `None` when the expression leaves
/// the fragment (see the module docs): a `var`-dependent non-`sqrt` function
/// head anywhere (including inside a base), a symbolic exponent on a
/// `var`-dependent base, or a budget overflow.
fn collect<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Collect<'a>> {
    let mut sites = Vec::new();
    collect_walk(ctx, expr, var, &mut sites)?;
    let mut bases: Vec<Atom<'a>> = Vec::new();
    for s in &sites {
        if !bases.contains(&s.base) {
            bases.push(s.base);
        }
    }
    if bases.len() > MAX_MR_BASES
        || sites
            .iter()
            .any(|s| s.q > MAX_MR_Q || s.p.abs() > MAX_MR_SITE_P)
    {
        return None;
    }
    Some(Collect { sites, bases })
}

fn collect_walk<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sites: &mut Vec<Site<'a>>,
) -> Option<()> {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(()),
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_walk(ctx, *a, var, sites)?;
            }
            Some(())
        }
        AtomNode::Fun(name, args) => {
            if name.as_str() == "sqrt" && args.len() == 1 && !is_constant(args[0], var) {
                let base = normalize(ctx, args[0]);
                if !base_heads_ok(base, var) {
                    return None;
                }
                sites.push(Site { base, p: 1, q: 2 });
                // The base may itself contain radicals (MR-I nesting).
                collect_walk(ctx, args[0], var, sites)
            } else if is_constant(expr, var) {
                // A constant head (`sqrt(3)`, `sin(2)`, …) is an inert
                // coefficient.
                Some(())
            } else {
                // A transcendental kernel of the variable: another stage's
                // shape (the W5 interface contract).
                None
            }
        }
        AtomNode::Pow(b, e) => {
            // `sqrt(S)^n`: odd `n` is a genuine site `(S, n, 2)`; even `n`
            // is exactly `S^(n/2)` (no radical), and either way nested
            // radicals inside `S` are collected from `S` directly.
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && !is_constant(args[0], var)
            {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                if n % 2 != 0 {
                    let base = normalize(ctx, args[0]);
                    if !base_heads_ok(base, var) {
                        return None;
                    }
                    sites.push(Site { base, p: *n, q: 2 });
                }
                return collect_walk(ctx, args[0], var, sites);
            }
            if let Some((base_raw, p, q)) = site_of(expr, var) {
                let base = normalize(ctx, base_raw);
                if !base_heads_ok(base, var) {
                    return None;
                }
                sites.push(Site { base, p, q });
                return collect_walk(ctx, base_raw, var, sites);
            }
            match e.node() {
                AtomNode::Num(_) => collect_walk(ctx, *b, var, sites),
                _ => {
                    if is_constant(*b, var) && is_constant(*e, var) {
                        Some(())
                    } else {
                        // Symbolic exponent on a `var`-dependent base, or an
                        // exponential: outside the fragment.
                        None
                    }
                }
            }
        }
    }
}

/// A radical base may not contain any `var`-dependent function head other
/// than (nested) `sqrt` — the W5/INV_SUBST contract. Constant heads are
/// inert coefficients.
fn base_heads_ok(base: Atom<'_>, var: Symbol) -> bool {
    match base.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|a| base_heads_ok(*a, var)),
        AtomNode::Pow(b, e) => base_heads_ok(*b, var) && base_heads_ok(*e, var),
        AtomNode::Fun(name, args) => {
            if is_constant(base, var) {
                return true;
            }
            name.as_str() == "sqrt" && args.len() == 1 && base_heads_ok(args[0], var)
        }
    }
}

// =========================================================================
// Numeric evidence gate (merges only; spec §4.1-2)
// =========================================================================

/// Sample abscissae (the same spread as the project's numerical oracle).
const GATE_SAMPLES: [f64; 8] = [-1.7, -0.9, -0.37, 0.31, 0.77, 1.23, 1.91, 2.53];
/// Deterministic dummy values for the free parameters.
const GATE_PARAMS: [f64; 10] = [0.7, 1.3, 0.37, 0.91, 1.73, 0.53, 1.11, 0.83, 1.47, 0.61];
/// Minimum number of usable sample points for a merge to be accepted.
const GATE_MIN_SAMPLES: usize = 2;

/// `f64` evaluator for the evidence gate; `None` for a head it does not
/// model (the gate then declines).
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
                "sqrt" => v.sqrt(),
                _ => return None,
            })
        }
    }
}

/// Collect the symbols of `expr` other than the integration variable.
fn collect_free_symbols(expr: Atom<'_>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Num(_) => {}
        AtomNode::Var(v) => {
            if *v != var {
                out.push(*v);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                collect_free_symbols(*a, var, out);
            }
        }
        AtomNode::Pow(b, e) => {
            collect_free_symbols(*b, var, out);
            collect_free_symbols(*e, var, out);
        }
    }
}

/// Numeric evidence that `rewritten` equals `original` on the integrand's
/// real domain (spec §4.1-2): at every sample point where **both** sides
/// evaluate to finite reals, they must agree; at least
/// [`GATE_MIN_SAMPLES`] such points are required. A point where the original
/// is not real (a negative radicand under a principal `sqrt`) but the
/// rewrite is — the merged radicand `u·v > 0` with `u, v < 0` — is outside
/// the original's real domain and is skipped, not counted.
fn gate_equivalent(original: Atom<'_>, rewritten: Atom<'_>, var: Symbol) -> bool {
    let mut symbols = Vec::new();
    collect_free_symbols(original, var, &mut symbols);
    collect_free_symbols(rewritten, var, &mut symbols);
    symbols.sort_by_key(|s| s.as_str().to_string());
    symbols.dedup();
    if symbols.len() > GATE_PARAMS.len() {
        return false;
    }
    let mut usable = 0usize;
    for &xv in &GATE_SAMPLES {
        let mut env: Vec<(Symbol, f64)> = symbols
            .iter()
            .enumerate()
            .map(|(i, s)| (*s, GATE_PARAMS[i]))
            .collect();
        env.push((var, xv));
        let (Some(o), Some(r)) = (eval_f64(original, &env), eval_f64(rewritten, &env)) else {
            continue;
        };
        if !o.is_finite() || !r.is_finite() {
            continue;
        }
        usable += 1;
        if (o - r).abs() > 1e-6 * o.abs().max(1.0) {
            return false;
        }
    }
    usable >= GATE_MIN_SAMPLES
}

// =========================================================================
// Constant-radical splitting (certificate surface)
// =========================================================================

/// Split a fractional power of a *constant* product over its factors:
/// `(4·b·d⁻¹)^(−1/2) → 4^(−1/2)·b^(−1/2)·d^(1/2)`. Exact in the radical
/// field (both sides square to the same element), and needed because the
/// certificate's radical layer requires each relation base to be a
/// polynomial — a quotient or product base like `4·b·d⁻¹` is outside its
/// fragment, while the atomic factors (`4`, `b`, `d`) are not.
fn split_const_power<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, p: i64, q: i64) -> Atom<'a> {
    match base.node() {
        AtomNode::Mul(args) => {
            let parts: Vec<Atom<'a>> = args
                .iter()
                .map(|f| split_const_power(ctx, *f, p, q))
                .collect();
            mul_atoms(ctx, &parts)
        }
        AtomNode::Pow(g, e2) => {
            if let Some((p2, q2)) = rat_exponent(*e2)
                && let (Some(n), Some(d)) = (p.checked_mul(p2), q.checked_mul(q2))
            {
                let g_ = gcd_i64(n, d).max(1);
                return ctx.pow(*g, rat_atom(ctx, n / g_, d / g_));
            }
            ctx.pow(base, rat_atom(ctx, p, q))
        }
        _ => ctx.pow(base, rat_atom(ctx, p, q)),
    }
}

/// Apply [`split_const_power`] to every constant-base fractional power in
/// `expr` (bottom-up). Numeric constants are folded atoms already, so only
/// symbolic-constant answers change.
fn split_constant_radicals<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let rebuilt = match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_constant_radicals(ctx, *a, var))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_constant_radicals(ctx, *a, var))
                .collect();
            ctx.mul(&kids)
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_constant_radicals(ctx, *a, var))
                .collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Pow(b, e) => {
            let nb = split_constant_radicals(ctx, *b, var);
            if let Some((p, q)) = rat_exponent(*e)
                && q >= 2
                && is_constant(nb, var)
                && matches!(nb.node(), AtomNode::Mul(_) | AtomNode::Pow(_, _))
            {
                return split_const_power(ctx, nb, p, q);
            }
            ctx.pow(nb, *e)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    };
    normalize(ctx, rebuilt)
}

// =========================================================================
// Möbius rationalization (MR-L)
// =========================================================================

/// The substitution data for `t^L = u/v` over affine bases
/// `u = a₀ + a₁·x`, `v = c₀ + c₁·x`.
struct Mobius<'a> {
    u: Atom<'a>,
    v: Atom<'a>,
    l: i64,
    t: Atom<'a>,
    /// `v(x(t)) = Δ/(c₁·t^L − a₁)`.
    v_t: Atom<'a>,
    /// `x(t) = (a₀ − c₀·t^L)/(c₁·t^L − a₁)`.
    x_t: Atom<'a>,
}

/// The Möbius substitution `t^L = u/v` (spec §2.2). Returns `None` unless
/// both bases are affine with `Δ ≢ 0`, `L = lcm(qᵢ) ∈ {2, 4}`, and the
/// substituted t-form is genuinely rational in `t` (a leftover radical —
/// e.g. a third base or an uncombinable half power — declines).
fn try_mobius<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let u = coll.bases[0];
    let v = coll.bases[1];
    let (a1, a0) = linear_form(ctx, u, var)?;
    let (c1, c0) = linear_form(ctx, v, var)?;
    if identically_zero(ctx, a1) || identically_zero(ctx, c1) {
        return None;
    }
    let mut l = 1i64;
    for s in &coll.sites {
        l = lcm_i64(l, s.q)?;
    }
    if l != 2 && l != 4 {
        return None;
    }
    for s in &coll.sites {
        // `q | l` holds by construction; the guard is on the t-power size.
        let te = s.p.checked_mul(l / s.q)?;
        if te.abs() > MAX_MR_T_EXP {
            return None;
        }
    }
    // Δ = a₀·c₁ − a₁·c₀; a vanishing determinant means proportional bases —
    // a single-base shape in disguise.
    let delta = crate::ode::util::collect_terms(
        ctx,
        normalize(
            ctx,
            ctx.add(&[ctx.mul(&[a0, c1]), neg(ctx, ctx.mul(&[a1, c0]))]),
        ),
    );
    if matches!(delta.node(), AtomNode::Num(0)) {
        return None;
    }
    let t_sym = pick_subst_symbol(expr, var)?;
    let t = ctx.var(t_sym.as_str());
    let tl = int_pow(ctx, t, l);
    // D = c₁·t^L − a₁.
    let d_t = normalize(ctx, ctx.add(&[ctx.mul(&[c1, tl]), neg(ctx, a1)]));
    let v_t = normalize(ctx, ctx.mul(&[delta, inv(ctx, d_t)]));
    let x_t = normalize(
        ctx,
        ctx.mul(&[
            normalize(ctx, ctx.add(&[a0, neg(ctx, ctx.mul(&[c0, tl]))])),
            inv(ctx, d_t),
        ]),
    );
    // dx/dt = −L·Δ·t^(L−1)/D².
    let dx_dt = normalize(
        ctx,
        ctx.mul(&[
            neg(ctx, delta),
            ctx.num(l),
            int_pow(ctx, t, l - 1),
            int_pow(ctx, d_t, -2),
        ]),
    );
    let pair = Mobius {
        u,
        v,
        l,
        t,
        v_t,
        x_t,
    };
    let subbed = mobius_subst(ctx, expr, var, &pair)?;
    // merge → distribute → merge: the first merge combines the `V` powers
    // (`V^(1/2)·V^(−3/2) → V⁻¹`), the distribution then flattens the
    // resulting product powers (`V⁻¹ = (Δ·D⁻¹)⁻¹ → Δ⁻¹·D`), and the second
    // merge sums the `D` powers into one denominator.
    let t_form = merge_radical_powers(
        ctx,
        distribute_int_powers(
            ctx,
            merge_radical_powers(ctx, normalize(ctx, ctx.mul(&[subbed, dx_dt])), t_sym),
        ),
        t_sym,
    );
    if node_count(t_form) > MAX_MR_SUBST_NODES || charge_mr_work(node_count(t_form) as u64) {
        return None;
    }
    if !is_rational_tform(t_form, t_sym, var) {
        return None;
    }
    let result_t = integrate_raw(ctx, t_form, t_sym, 0, true, 0, 0);
    if contains_integral(result_t) {
        return None;
    }
    // Back-substitute `t = u^(1/L)·v^(−1/L)` — per-base powers, *not* the
    // quotient radical `(u/v)^(1/L)`: the certificate's radical layer needs
    // each generator's relation base to be a polynomial in `x`, and a
    // quotient base has a nontrivial denominator (spec §3.3's `√(u/v)`
    // shape is the same function on the `u, v > 0` sheet).
    let t_back = normalize(
        ctx,
        ctx.mul(&[
            ctx.pow(u, rat_atom(ctx, 1, l)),
            ctx.pow(v, rat_atom(ctx, -1, l)),
        ]),
    );
    let back = replace_symbol(ctx, result_t, t_sym, t_back);
    Some(normalize(ctx, back))
}

/// Rewrite `expr` under the Möbius substitution: sites of `u` become
/// `t^(pL/q)·V^(p/q)`, sites of `v` become `V^(p/q)`, and `x` becomes
/// `x_t`. Any other site (a third base) or non-rational structure returns
/// `None` — the structural rationality check downstream is the final word.
fn mobius_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    mp: &Mobius<'a>,
) -> Option<Atom<'a>> {
    if let Some((base_raw, p, q)) = site_of(expr, var) {
        let base = normalize(ctx, base_raw);
        if base == mp.u {
            let te = p.checked_mul(mp.l / q)?;
            return Some(ctx.mul(&[int_pow(ctx, mp.t, te), ctx.pow(mp.v_t, rat_atom(ctx, p, q))]));
        }
        if base == mp.v {
            return Some(ctx.pow(mp.v_t, rat_atom(ctx, p, q)));
        }
        return None;
    }
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        AtomNode::Var(vv) => {
            if *vv == var {
                Some(mp.x_t)
            } else {
                Some(expr)
            }
        }
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| mobius_subst(ctx, *a, var, mp))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| mobius_subst(ctx, *a, var, mp))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => {
            if is_constant(expr, var) {
                return Some(expr);
            }
            // Only integer exponents can remain (fractional ones are sites).
            let nb = mobius_subst(ctx, *b, var, mp)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Fun(_, _) => {
            if is_constant(expr, var) {
                Some(expr)
            } else {
                // The collector rejected transcendental kernels, so this is
                // unreachable in practice; decline rather than guess.
                None
            }
        }
    }
}

/// Structural check: `expr` is a rational function of the substitution
/// variable `t_sym`, modulo `t_sym`-constant radical coefficients (a
/// `sqrt(3)`-style factor is a coefficient, not a radical of `t`). No
/// leftover occurrence of the original integration variable is allowed.
fn is_rational_tform(expr: Atom<'_>, t_sym: Symbol, orig_var: Symbol) -> bool {
    if super::contains_symbol(expr, orig_var) {
        return false;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().all(|a| is_rational_tform(*a, t_sym, orig_var))
        }
        AtomNode::Pow(b, e) => {
            if is_constant(*b, t_sym) {
                return true;
            }
            matches!(e.node(), AtomNode::Num(_)) && is_rational_tform(*b, t_sym, orig_var)
        }
        AtomNode::Fun(_, _) => is_constant(expr, t_sym),
    }
}

// =========================================================================
// Same-radicand merge (spec §4.1)
// =========================================================================

/// Per-base exponents of a product (`(base, p, q)`; `q = 1` for integer
/// powers) plus the remaining site-free factors.
type PowerMap<'a> = (Vec<(Atom<'a>, i64, i64)>, Vec<Atom<'a>>);

/// Decompose the top-level product of `expr` into per-base exponents (every
/// radical site must be a *direct factor* of the product — a site nested in
/// a sum or power declines) and the remaining site-free factors.
/// Integer-power factors of a base are folded into that base's exponent.
fn product_power_map<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<PowerMap<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut map: Vec<(Atom<'a>, i64, i64)> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    // `1/(√u·√v)` parses as `(√u·√v)^-1`: split integer powers of products
    // into per-factor powers (exact for integer exponents) so the sites are
    // seen as direct factors.
    let mut queue = factors;
    while let Some(f) = queue.pop() {
        if let AtomNode::Pow(b, e) = f.node()
            && let AtomNode::Mul(args) = b.node()
            && let AtomNode::Num(n) = e.node()
        {
            for a in args.iter() {
                queue.push(int_pow(ctx, *a, *n));
            }
            continue;
        }
        if let Some((base, p, q)) = radical_factor(ctx, f, var) {
            if let Some(slot) = map.iter_mut().find(|(b, _, _)| *b == base) {
                let num = slot
                    .1
                    .checked_mul(q)
                    .and_then(|v| v.checked_add(p.checked_mul(slot.2)?))?;
                let den = slot.2.checked_mul(q)?;
                let g = gcd_i64(num, den).max(1);
                slot.1 = num / g;
                slot.2 = den / g;
            } else {
                map.push((base, p, q));
            }
        } else {
            if contains_radical(f, var) {
                return None;
            }
            rest.push(f);
        }
    }
    // Every collected base must appear as a direct factor.
    if map.iter().filter(|(_, _, q)| *q >= 2).count() != coll.bases.len() {
        return None;
    }
    Some((map, rest))
}

/// The degree of a polynomial in `var`; `None` when not a polynomial
/// (negative/fractional powers, other structure).
fn poly_degree(expr: Atom<'_>, var: Symbol) -> Option<usize> {
    match expr.node() {
        AtomNode::Num(_) => Some(0),
        AtomNode::Var(v) => Some(usize::from(*v == var)),
        AtomNode::Add(args) => args
            .iter()
            .try_fold(0usize, |acc, a| Some(acc.max(poly_degree(*a, var)?))),
        AtomNode::Mul(args) => args
            .iter()
            .try_fold(0usize, |acc, a| acc.checked_add(poly_degree(*a, var)?)),
        AtomNode::Pow(b, e) => {
            let AtomNode::Num(n) = e.node() else {
                return None;
            };
            if *n < 0 {
                return None;
            }
            poly_degree(*b, var)?.checked_mul(usize::try_from(*n).ok()?)
        }
        AtomNode::Fun(_, _) => None,
    }
}

/// A merge performed during the reduction: the merged radicand (expanded,
/// normalized — the spelling `sqrt_quadratic`'s answers carry) and the two
/// original bases. The final answer is expressed in the **original** bases
/// again (spec §4.1-4: the certificate certifies `D(F) − f` with the
/// radicals in their original bases; the merge only steers the solving
/// path), by splitting every `(u·v)^(p/q)` site back into
/// `u^(p/q)·v^(p/q)` — exact on the gate-validated real domain.
struct MergeRecord<'a> {
    merged: Atom<'a>,
    b1: Atom<'a>,
    b2: Atom<'a>,
}

/// Split every radical site whose base is `merged` into the per-base powers
/// of `b1` and `b2` (see [`MergeRecord`]). Integer powers of the merged
/// base are polynomial and stay.
fn split_base_in<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    merged: Atom<'a>,
    b1: Atom<'a>,
    b2: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if name.as_str() == "sqrt" && args.len() == 1 && normalize(ctx, args[0]) == merged {
                return mul_atoms(ctx, &[ctx.fun("sqrt", &[b1]), ctx.fun("sqrt", &[b2])]);
            }
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_base_in(ctx, *a, merged, b1, b2))
                .collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Pow(b, e) => {
            // `sqrt(u·v)^n` spelling.
            if let AtomNode::Fun(name, args) = b.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && normalize(ctx, args[0]) == merged
                && let AtomNode::Num(n) = e.node()
            {
                if n % 2 == 0 {
                    return int_pow(ctx, merged, n / 2);
                }
                return mul_atoms(
                    ctx,
                    &[
                        int_pow(ctx, ctx.fun("sqrt", &[b1]), *n),
                        int_pow(ctx, ctx.fun("sqrt", &[b2]), *n),
                    ],
                );
            }
            if let Some((p, q)) = rat_exponent(*e)
                && q >= 2
                && normalize(ctx, *b) == merged
            {
                return mul_atoms(
                    ctx,
                    &[
                        ctx.pow(b1, rat_atom(ctx, p, q)),
                        ctx.pow(b2, rat_atom(ctx, p, q)),
                    ],
                );
            }
            let nb = split_base_in(ctx, *b, merged, b1, b2);
            let ne = split_base_in(ctx, *e, merged, b1, b2);
            if nb == *b && ne == *e {
                expr
            } else {
                ctx.pow(nb, ne)
            }
        }
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_base_in(ctx, *a, merged, b1, b2))
                .collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args
                .iter()
                .map(|a| split_base_in(ctx, *a, merged, b1, b2))
                .collect();
            ctx.mul(&kids)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Express `answer` in the original bases of every recorded merge.
fn split_merged_bases<'a>(
    ctx: &'a AtomArena<'a>,
    answer: Atom<'a>,
    records: &[MergeRecord<'a>],
) -> Atom<'a> {
    let mut out = answer;
    for rec in records {
        out = split_base_in(ctx, out, rec.merged, rec.b1, rec.b2);
    }
    normalize(ctx, out)
}

/// Read `term` as `coeff·var^k` with `k ∈ {0, 1, 2}` (the coefficient may
/// be any product of constants). `None` outside that shape.
fn term_coeff_pow<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, usize)> {
    match term.node() {
        AtomNode::Var(v) if *v == var => Some((ctx.num(1), 1)),
        AtomNode::Pow(b, e)
            if matches!(b.node(), AtomNode::Var(v) if *v == var)
                && matches!(e.node(), AtomNode::Num(1 | 2)) =>
        {
            let AtomNode::Num(n) = e.node() else {
                return None;
            };
            Some((ctx.num(1), usize::try_from(*n).ok()?))
        }
        AtomNode::Mul(args) => {
            let mut coeff: Vec<Atom<'a>> = Vec::new();
            let mut k = 0usize;
            for f in args.iter() {
                match f.node() {
                    AtomNode::Var(v) if *v == var => k += 1,
                    AtomNode::Pow(b, e)
                        if matches!(b.node(), AtomNode::Var(v) if *v == var)
                            && matches!(e.node(), AtomNode::Num(1 | 2)) =>
                    {
                        let AtomNode::Num(n) = e.node() else {
                            return None;
                        };
                        k += usize::try_from(*n).ok()?;
                    }
                    _ if is_constant(*f, var) => coeff.push(*f),
                    _ => return None,
                }
            }
            if k > 2 {
                return None;
            }
            Some((mul_atoms(ctx, &coeff), k))
        }
        _ if is_constant(term, var) => Some((term, 0)),
        _ => None,
    }
}

/// Coefficients `[c0, c1, c2]` of a polynomial of degree ≤ 2 in `var`;
/// `None` outside that shape.
fn poly_coeffs_le2<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<[Atom<'a>; 3]> {
    let terms: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut out = [ctx.num(0), ctx.num(0), ctx.num(0)];
    for t in terms {
        let (c, k) = term_coeff_pow(ctx, t, var)?;
        out[k] = normalize(ctx, ctx.add(&[out[k], c]));
    }
    Some(out)
}

/// The product of two bases, expanded, folded, and respelled as
/// `c0 + c1·x + c2·x²` with the zero coefficients dropped. The respelling
/// matters twice: `quadratic_coeffs` (the parser shared by
/// `sqrt_quadratic` and rule family G) only reads two-factor terms, while
/// `collect_terms` spells `8·b·x` as a three-factor product; and
/// `sqrt_quadratic` builds its answer's radicand in exactly this canonical
/// spelling, which is the merge record's key.
fn base_product_expanded<'a>(
    ctx: &'a AtomArena<'a>,
    b1: Atom<'a>,
    b2: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let prod = normalize(ctx, ctx.mul(&[b1, b2]));
    let prod = crate::expand::expand_bounded(ctx, prod)?;
    let prod = crate::ode::util::collect_terms(ctx, prod);
    let [c0, c1, c2] = poly_coeffs_le2(ctx, prod, var)?;
    let x = ctx.var(var.as_str());
    let mut terms = vec![c0];
    if !matches!(c1.node(), AtomNode::Num(0)) {
        terms.push(normalize(ctx, ctx.mul(&[c1, x])));
    }
    if !matches!(c2.node(), AtomNode::Num(0)) {
        terms.push(normalize(ctx, ctx.mul(&[c2, int_pow(ctx, x, 2)])));
    }
    Some(normalize(ctx, ctx.add(&terms)))
}

/// Merge a product-shape pair of radical bases into one base
/// (`rest·u^e·v^e → rest·(u·v)^e`), gated by the degree bound (≤ 2, the
/// `sqrt_quadratic` acceptance surface) and the numeric evidence gate.
/// Returns the merged *integrand* (with the radicand in expanded form).
fn merge_pair<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    map: &[(Atom<'a>, i64, i64)],
    rest: &[Atom<'a>],
    i: usize,
    j: usize,
) -> Option<(Atom<'a>, MergeRecord<'a>)> {
    let (b1, p1, q1) = map[i];
    let (b2, p2, q2) = map[j];
    if q1 < 2 || q2 < 2 || (p1, q1) != (p2, q2) {
        return None;
    }
    // Degenerate merged radicands decline (constant or empty product).
    if is_constant(b1, var) || is_constant(b2, var) {
        return None;
    }
    let uv = base_product_expanded(ctx, b1, b2, var)?;
    let deg = poly_degree(uv, var)?;
    if deg == 0 || deg > 2 {
        return None;
    }
    // The merged site is emitted in `Pow` spelling: `binomial`'s Chebyshev
    // matcher does not accept the `sqrt` head (measured: the G2 cross term
    // `x³·√(1−x²)` needs it).
    let merged_site = ctx.pow(uv, rat_atom(ctx, p1, q1));
    let mut factors: Vec<Atom<'a>> = rest.to_vec();
    for (k, (b, p, q)) in map.iter().enumerate() {
        if k == i || k == j {
            continue;
        }
        factors.push(if *q == 1 {
            int_pow(ctx, *b, *p)
        } else {
            ctx.pow(*b, rat_atom(ctx, *p, *q))
        });
    }
    factors.push(merged_site);
    let merged = mul_atoms(ctx, &factors);
    if !gate_equivalent(expr, merged, var) {
        return None;
    }
    Some((merged, MergeRecord { merged: uv, b1, b2 }))
}

/// Integrate the merged single-base form with the single-base engines
/// directly (the `kernel_subst::try_engines` precedent): `binomial`'s
/// Chebyshev cases first (its t-forms are the small answers), then
/// `sqrt_quadratic`'s direct forms, then the rational engine. The whole
/// chain is *not* re-entered: its stage order would let `sqrt_quadratic`'s
/// Euler III claim Chebyshev shapes with a degree-10 t-form (measured on
/// the G2 cross term `x³·√(1−x²)`), producing answers too large to
/// certify. A residue declines.
fn integrate_single_base<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if let Some(r) = super::binomial::integrate_binomial(ctx, expr, var)
        && !contains_integral(r)
    {
        return Some(r);
    }
    if let Some(r) = super::sqrt_quadratic::integrate_sqrt_quadratic(ctx, expr, var)
        && !contains_integral(r)
    {
        return Some(r);
    }
    if let Some(r) = super::rational::integrate_rational(ctx, expr, var)
        && !contains_integral(r)
    {
        return Some(r);
    }
    None
}

/// Two-base merge (`√u·√v → √(u·v)` and its reciprocal form): spec §2.2's
/// "同根合并子情形", which takes priority over the Möbius substitution for
/// product shapes. The merged single-base integrand is handed back to the
/// chain; a residue there declines the whole merge (spec §1.3-3).
fn try_merge_two<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
    records: &mut Vec<MergeRecord<'a>>,
) -> Option<Atom<'a>> {
    let (map, rest) = product_power_map(ctx, expr, var, coll)?;
    if map.iter().filter(|(_, _, q)| *q >= 2).count() != 2 {
        return None;
    }
    let idx: Vec<usize> = map
        .iter()
        .enumerate()
        .filter(|(_, (_, _, q))| *q >= 2)
        .map(|(k, _)| k)
        .collect();
    let (merged, rec) = merge_pair(ctx, expr, var, &map, &rest, idx[0], idx[1])?;
    if charge_mr_work(node_count(merged) as u64) {
        return None;
    }
    let r = integrate_single_base(ctx, merged, var)?;
    records.push(rec);
    Some(normalize(ctx, r))
}

/// With ≥ 3 bases, collapse one mergeable pair per round (spec §2.2's
/// pairing loop): same-slope affine pairs first (the `01384` shape), then
/// any pair whose product stays within degree ≤ 2. Returns the merged
/// integrand for the next round.
fn try_merge_round<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
    records: &mut Vec<MergeRecord<'a>>,
) -> Option<Atom<'a>> {
    let (map, rest) = product_power_map(ctx, expr, var, coll)?;
    let idx: Vec<usize> = map
        .iter()
        .enumerate()
        .filter(|(_, (_, _, q))| *q >= 2)
        .map(|(k, _)| k)
        .collect();
    if idx.len() < 3 {
        return None;
    }
    // Same-slope affine pairs first.
    for &i in &idx {
        for &j in &idx {
            if i >= j {
                continue;
            }
            let (Some((s1, _)), Some((s2, _))) = (
                linear_form(ctx, map[i].0, var),
                linear_form(ctx, map[j].0, var),
            ) else {
                continue;
            };
            if !identically_zero(ctx, normalize(ctx, ctx.add(&[s1, neg(ctx, s2)]))) {
                continue;
            }
            if let Some((m, rec)) = merge_pair(ctx, expr, var, &map, &rest, i, j) {
                records.push(rec);
                return Some(m);
            }
        }
    }
    // Then any mergeable pair, in collection order.
    for &i in &idx {
        for &j in &idx {
            if i >= j {
                continue;
            }
            if let Some((m, rec)) = merge_pair(ctx, expr, var, &map, &rest, i, j) {
                records.push(rec);
                return Some(m);
            }
        }
    }
    None
}

/// Fold integer powers of `sqrt` heads: `(√S)^n → S^(n/2)` (even `n`) or
/// `S^((n−1)/2)·√S` (odd `n`) — exact identities, and the reason a denested
/// answer's radicands spell themselves with the integrand's bases (the
/// certificate's radical layer keys generators by base).
fn fold_sqrt_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    let rebuilt = match expr.node() {
        AtomNode::Add(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_sqrt_powers(ctx, *a)).collect();
            ctx.add(&kids)
        }
        AtomNode::Mul(args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_sqrt_powers(ctx, *a)).collect();
            ctx.mul(&kids)
        }
        AtomNode::Fun(name, args) => {
            let kids: Vec<Atom<'a>> = args.iter().map(|a| fold_sqrt_powers(ctx, *a)).collect();
            ctx.fun(name.as_str(), &kids)
        }
        AtomNode::Pow(b, e) => {
            let nb = fold_sqrt_powers(ctx, *b);
            if let AtomNode::Fun(name, args) = nb.node()
                && name.as_str() == "sqrt"
                && args.len() == 1
                && let AtomNode::Num(n) = e.node()
                && n.abs() >= 2
            {
                if n % 2 == 0 {
                    return int_pow(ctx, args[0], n / 2);
                }
                return mul_atoms(ctx, &[int_pow(ctx, args[0], (n - 1) / 2), nb]);
            }
            ctx.pow(nb, *e)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    };
    normalize(ctx, rebuilt)
}

// =========================================================================
// MR-G: conjugate rationalization of (√u ± √v)^k factors (spec §2.5)
// =========================================================================

/// Read `term` as `±sqrt(w)`: returns `(w, sign)`.
fn signed_sqrt_term(term: Atom<'_>) -> Option<(Atom<'_>, i64)> {
    match term.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => {
            Some((args[0], 1))
        }
        AtomNode::Mul(args) => {
            let mut sign = 1i64;
            let mut rad: Option<Atom<'_>> = None;
            for f in args.iter() {
                match f.node() {
                    AtomNode::Num(-1) => sign = -sign,
                    AtomNode::Fun(name, args2)
                        if name.as_str() == "sqrt" && args2.len() == 1 && rad.is_none() =>
                    {
                        rad = Some(args2[0]);
                    }
                    _ => return None,
                }
            }
            rad.map(|w| (w, sign))
        }
        _ => None,
    }
}

/// Find a factor `(±√u ± √v)^k` with `k ∈ {−2, −1, 2}` and `u, v` affine in
/// `var`; returns `(u, v, s1, s2, k)`.
fn radical_sum_factor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, i64, i64, i64)> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    for f in factors {
        let AtomNode::Pow(b, e) = f.node() else {
            continue;
        };
        let AtomNode::Num(k) = e.node() else {
            continue;
        };
        if !matches!(k, -2 | -1 | 2) {
            continue;
        }
        let AtomNode::Add(args) = b.node() else {
            continue;
        };
        if args.len() != 2 {
            continue;
        }
        let (Some((u, s1)), Some((v, s2))) = (signed_sqrt_term(args[0]), signed_sqrt_term(args[1]))
        else {
            continue;
        };
        // Affine bases keep `u − v` polynomial (constant when the slopes
        // agree — the corpus shape).
        if linear_form(ctx, u, var).is_none() || linear_form(ctx, v, var).is_none() {
            continue;
        }
        return Some((normalize(ctx, u), normalize(ctx, v), s1, s2, *k));
    }
    None
}

/// Conjugate rationalization, exact identities:
/// `(s1·√u + s2·√v)² = u + v + 2·s1·s2·√u·√v`,
/// `(s1·√u + s2·√v)^(−1) = (s1·√u − s2·√v)/(u − v)`,
/// `(s1·√u + s2·√v)^(−2) = (u + v − 2·s1·s2·√u·√v)/(u − v)²`.
/// The rewritten integrand is handed back to the chain: the `√u·√v` cross
/// term is exactly the merge path's product shape, and the two single-base
/// halves distribute. `k = 1` needs no conjugate (the chain distributes the
/// sum); `|k| > 2` exceeds the denominator-power budget (spec §2.5).
fn try_conjugate<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    _coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let (u, v, s1, s2, k) = radical_sum_factor(ctx, expr, var)?;
    let su = ctx.fun("sqrt", &[u]);
    let sv = ctx.fun("sqrt", &[v]);
    let s = s1 * s2;
    let two_s_uv = ctx.mul(&[ctx.num(2 * s), su, sv]);
    let u_plus_v = normalize(ctx, ctx.add(&[u, v]));
    let replacement = match k {
        2 => normalize(ctx, ctx.add(&[u_plus_v, two_s_uv])),
        -1 | -2 => {
            let diff =
                crate::ode::util::collect_terms(ctx, normalize(ctx, ctx.add(&[u, neg(ctx, v)])));
            if matches!(diff.node(), AtomNode::Num(0)) {
                return None;
            }
            let numerator = if k == -1 {
                // s1·√u − s2·√v
                normalize(
                    ctx,
                    ctx.add(&[ctx.mul(&[ctx.num(s1), su]), ctx.mul(&[ctx.num(-s2), sv])]),
                )
            } else {
                normalize(ctx, ctx.add(&[u_plus_v, neg(ctx, two_s_uv)]))
            };
            normalize(ctx, ctx.mul(&[numerator, int_pow(ctx, diff, k)]))
        }
        _ => return None,
    };
    // Rebuild the integrand with the sum-power factor replaced. The factor
    // is identified structurally: `Pow(Add of two ±√ terms, k)`.
    let rebuilt = replace_factor(ctx, expr, |f| {
        if let AtomNode::Pow(b, e) = f.node()
            && let AtomNode::Num(n) = e.node()
            && *n == k
            && let AtomNode::Add(args) = b.node()
            && args.len() == 2
            && signed_sqrt_term(args[0]).is_some()
            && signed_sqrt_term(args[1]).is_some()
        {
            return Some(replacement);
        }
        None
    })?;
    if charge_mr_work(node_count(rebuilt) as u64) {
        return None;
    }
    // Distribute first: the `√u·√v` cross term must reach the merge path as
    // its own (single-product) integral — integrating the undistributed
    // form lets a Möbius re-entry swallow the whole sum (measured: a
    // sextic-denominator t-form, > 60 s and an uncertifiable answer shape).
    let distributed = match crate::expand::expand_bounded(ctx, rebuilt) {
        Some(e) => crate::ode::util::collect_terms(ctx, e),
        None => rebuilt,
    };
    let r = integrate_raw(ctx, distributed, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    Some(normalize(ctx, r))
}

/// Replace the first top-level product factor for which `f` returns
/// `Some(_)`; `None` when no factor matches.
fn replace_factor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    f: impl Fn(Atom<'a>) -> Option<Atom<'a>>,
) -> Option<Atom<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut replaced = false;
    let mut out = Vec::with_capacity(factors.len());
    for fac in factors {
        if !replaced && let Some(r) = f(fac) {
            replaced = true;
            out.push(r);
            continue;
        }
        out.push(fac);
    }
    if !replaced {
        return None;
    }
    Some(mul_atoms(ctx, &out))
}

// =========================================================================
// MR-R: reciprocal-power bases `a + b/x` (spec §2.8)
// =========================================================================

/// Read `term` as `c·x⁻¹` with `c` constant; returns `c`.
fn reciprocal_term<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    match term.node() {
        AtomNode::Pow(b, e)
            if matches!(b.node(), AtomNode::Var(v) if *v == var)
                && matches!(e.node(), AtomNode::Num(-1)) =>
        {
            Some(ctx.num(1))
        }
        AtomNode::Mul(args) => {
            let mut coeff: Vec<Atom<'a>> = Vec::new();
            let mut found = false;
            for f in args.iter() {
                if matches!(f.node(), AtomNode::Pow(b, e)
                    if matches!(b.node(), AtomNode::Var(v) if *v == var)
                        && matches!(e.node(), AtomNode::Num(-1)))
                {
                    if found {
                        return None;
                    }
                    found = true;
                } else if is_constant(*f, var) {
                    coeff.push(*f);
                } else {
                    return None;
                }
            }
            if !found {
                return None;
            }
            Some(mul_atoms(ctx, &coeff))
        }
        _ => None,
    }
}

/// Read `base` as `c0 + c1/x`; returns `(c0, c1)`.
fn reciprocal_affine<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match base.node() {
        AtomNode::Add(args) => {
            let mut c0 = ctx.num(0);
            let mut c1 = ctx.num(0);
            for t in args.iter() {
                if is_constant(*t, var) {
                    c0 = ctx.add(&[c0, *t]);
                } else {
                    c1 = ctx.add(&[c1, reciprocal_term(ctx, *t, var)?]);
                }
            }
            Some((normalize(ctx, c0), normalize(ctx, c1)))
        }
        _ => {
            let c1 = reciprocal_term(ctx, base, var)?;
            Some((ctx.num(0), c1))
        }
    }
}

/// Rewrite `expr` under `u = 1/x`: `x → u⁻¹` and `x^e → u^(−e)` (so a
/// reciprocal-affine base `c0 + c1/x` becomes the affine base `c0 + c1·u`).
fn reciprocal_subst<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    u: Atom<'a>,
) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Num(_) => Some(expr),
        AtomNode::Var(v) => {
            if *v == var {
                Some(inv(ctx, u))
            } else {
                Some(expr)
            }
        }
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| reciprocal_subst(ctx, *a, var, u))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| reciprocal_subst(ctx, *a, var, u))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => {
            if matches!(b.node(), AtomNode::Var(v) if *v == var) {
                let (p, q) = rat_exponent(*e)?;
                return Some(if q == 1 {
                    int_pow(ctx, u, p.checked_neg()?)
                } else {
                    ctx.pow(u, rat_atom(ctx, p.checked_neg()?, q))
                });
            }
            let nb = reciprocal_subst(ctx, *b, var, u)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Fun(name, args) => {
            if is_constant(expr, var) {
                return Some(expr);
            }
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| reciprocal_subst(ctx, *a, var, u))
                .collect();
            Some(ctx.fun(name.as_str(), &kids?))
        }
    }
}

/// MR-R: every radical base is reciprocal-affine (`c0 + c1/x`) or `x`
/// itself. `u = 1/x` (with `dx = −du/u²`) makes the bases affine in `u`;
/// the chain (MR-L's Möbius) takes it from there. Back-substitution is
/// `u = x⁻¹`.
fn try_reciprocal<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let mut any_reciprocal = false;
    for &base in &coll.bases {
        if base == x {
            continue;
        }
        let (_c0, c1) = reciprocal_affine(ctx, base, var)?;
        if identically_zero(ctx, c1) {
            return None;
        }
        any_reciprocal = true;
    }
    if !any_reciprocal {
        return None;
    }
    let u_sym = pick_subst_symbol(expr, var)?;
    let u = ctx.var(u_sym.as_str());
    let subbed = reciprocal_subst(ctx, expr, var, u)?;
    // dx = −u⁻² du.
    let u_form = merge_radical_powers(
        ctx,
        normalize(ctx, ctx.mul(&[subbed, neg(ctx, int_pow(ctx, u, -2))])),
        u_sym,
    );
    if node_count(u_form) > MAX_MR_SUBST_NODES || charge_mr_work(node_count(u_form) as u64) {
        return None;
    }
    let r = integrate_raw(ctx, u_form, u_sym, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    let back = replace_symbol(ctx, r, u_sym, inv(ctx, x));
    Some(normalize(ctx, back))
}

// =========================================================================
// MR-I: nested-radical denesting via `t² = g` (spec §2.8)
// =========================================================================

/// Whether `base` contains a radical site whose own base is `target`
/// (normalized comparison).
fn base_contains_radical_of<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    target: Atom<'a>,
    var: Symbol,
) -> bool {
    if let Some((b, _, _)) = site_of(base, var) {
        return normalize(ctx, b) == target;
    }
    match base.node() {
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => args
            .iter()
            .any(|a| base_contains_radical_of(ctx, *a, target, var)),
        AtomNode::Pow(b, e) => {
            base_contains_radical_of(ctx, *b, target, var)
                || base_contains_radical_of(ctx, *e, target, var)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Replace every radical site with base `target` by `t^p` (the denesting
/// image of `g^(p/2)` under `t² = g`).
fn replace_sites_of_base<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    target: Atom<'a>,
    t: Atom<'a>,
) -> Option<Atom<'a>> {
    if let Some((b, p, q)) = site_of(expr, var) {
        if q == 2 && normalize(ctx, b) == target {
            return Some(int_pow(ctx, t, p));
        }
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(expr),
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| replace_sites_of_base(ctx, *a, var, target, t))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| replace_sites_of_base(ctx, *a, var, target, t))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => {
            let nb = replace_sites_of_base(ctx, *b, var, target, t)?;
            let ne = replace_sites_of_base(ctx, *e, var, target, t)?;
            Some(ctx.pow(nb, ne))
        }
        AtomNode::Fun(name, args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| replace_sites_of_base(ctx, *a, var, target, t))
                .collect();
            Some(ctx.fun(name.as_str(), &kids?))
        }
    }
}

/// MR-I: a radical base containing another radical (e.g.
/// `sqrt(1 − x − sqrt(x))`). Repeatedly (≤ [`MAX_MR_NEST`] rounds) peel the
/// innermost radical `g^(1/2)` with affine `g` by `t² = g` — the principal
/// sheet `t = √g`, valid on the `g ≥ 0` half of the domain. Each round
/// lowers the nesting depth; the denested form is integrated by the chain
/// and the substitutions are unwound in reverse order.
fn try_denest<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let mut cur = expr;
    let mut cur_var = var;
    let mut backs: Vec<(Symbol, Atom<'a>)> = Vec::new();
    for _ in 0..MAX_MR_NEST {
        let coll = collect(ctx, cur, cur_var)?;
        if !coll.sites.iter().any(|s| contains_radical(s.base, cur_var)) {
            break;
        }
        // The innermost site: a radical-free base occurring inside another
        // site's base.
        let inner = coll
            .sites
            .iter()
            .find(|s| {
                !contains_radical(s.base, cur_var)
                    && coll
                        .sites
                        .iter()
                        .any(|o| base_contains_radical_of(ctx, o.base, s.base, cur_var))
            })
            .copied()?;
        if inner.q != 2 {
            return None;
        }
        let (g1, g0) = linear_form(ctx, inner.base, cur_var)?;
        if identically_zero(ctx, g1) {
            return None;
        }
        let t_sym = pick_subst_symbol(cur, cur_var)?;
        let t = ctx.var(t_sym.as_str());
        // t² = g: cur_var = (t² − g₀)/g₁, d(cur_var) = 2t/g₁·dt.
        let replaced = replace_sites_of_base(ctx, cur, cur_var, inner.base, t)?;
        let t2 = int_pow(ctx, t, 2);
        let x_t = normalize(
            ctx,
            ctx.mul(&[normalize(ctx, ctx.add(&[t2, neg(ctx, g0)])), inv(ctx, g1)]),
        );
        let dvar = normalize(ctx, ctx.mul(&[ctx.num(2), t, inv(ctx, g1)]));
        let next = merge_radical_powers(
            ctx,
            distribute_int_powers(
                ctx,
                normalize(
                    ctx,
                    ctx.mul(&[replace_symbol(ctx, replaced, cur_var, x_t), dvar]),
                ),
            ),
            t_sym,
        );
        if node_count(next) > MAX_MR_SUBST_NODES || charge_mr_work(node_count(next) as u64) {
            return None;
        }
        backs.push((t_sym, ctx.fun("sqrt", &[inner.base])));
        cur = next;
        cur_var = t_sym;
    }
    // After the budgeted rounds the nesting must be gone.
    let coll = collect(ctx, cur, cur_var)?;
    if coll.sites.iter().any(|s| contains_radical(s.base, cur_var)) {
        return None;
    }
    let mut r = integrate_raw(ctx, cur, cur_var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    for (sym, back) in backs.iter().rev() {
        r = replace_symbol(ctx, r, *sym, *back);
    }
    Some(normalize(ctx, r))
}

// =========================================================================
// MR-P: perfect-square bases and the sheet factor σ (spec §2.7)
// =========================================================================

/// Rewrite every radical site whose base is the square of an affine `w`
/// (from `targets`) as `w^p · σ` with `σ = w·(w²)^(−1/2)`, pushing `σ` onto
/// `sigmas`. `(w²)^(p/2) = w^p·σ` holds for every odd `p` (formally
/// `σ² = 1` under the relation `r² = w²`); even `p` is the exact power
/// `w^p` with no sheet factor.
fn slice_rewrite<'a>(
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
            // Sites have q = 2 by construction of `targets`. Even `p` is
            // the exact power `w^p` with no sheet factor.
            if p % 2 == 0 {
                return Some(int_pow(ctx, *w, p));
            }
            let sigma = ctx.mul(&[*w, ctx.pow(nb, rat_atom(ctx, -1, 2))]);
            sigmas.push(sigma);
            return Some(int_pow(ctx, *w, p));
        }
        // A non-target site is left unchanged.
        return Some(expr);
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(expr),
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_rewrite(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_rewrite(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => {
            let nb = slice_rewrite(ctx, *b, var, targets, sigmas)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Fun(name, args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| slice_rewrite(ctx, *a, var, targets, sigmas))
                .collect();
            Some(ctx.fun(name.as_str(), &kids?))
        }
    }
}

/// True when `expr` carries a negative integer power (`k ≤ −2`) of a
/// `var`-dependent base whose coefficients mention a symbol other than
/// `var` (a symbolic coefficient). The sliced `rest` in MR-P combines such
/// a factor with a surviving fractional-power site: the single-base
/// rationalization (`binomial`'s `t = base^(1/q)`) then produces a rational
/// t-form with a *powered symbolic denominator* — the documented
/// `symbolic_rational` weak spot (0.32.0). Measured: rubi-01871
/// (`(d+ex)^(5/2)/((a+bx)²)^(3/2)` → rest `σ·(d+ex)^(5/2)·(a+bx)⁻³`)
/// burned the whole 256-entry chain budget (199 s) without ever solving.
/// The slice path declines this conjunction instead of feeding the chain.
fn has_symbolic_denom_power<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    factors.iter().any(|f| {
        let AtomNode::Pow(b, e) = f.node() else {
            return false;
        };
        let AtomNode::Num(k) = e.node() else {
            return false;
        };
        if *k > -2 || is_constant(*b, var) {
            return false;
        }
        // A free symbol other than the integration variable = symbolic
        // coefficient (numeric-coefficient bases make ℚ-rational t-forms,
        // which the ℚ engines finish quickly).
        fn has_other_var(e: Atom<'_>, var: Symbol) -> bool {
            match e.node() {
                AtomNode::Var(v) => *v != var,
                AtomNode::Add(args) | AtomNode::Mul(args) => {
                    args.iter().any(|a| has_other_var(*a, var))
                }
                AtomNode::Pow(b, x) => has_other_var(*b, var) || has_other_var(*x, var),
                AtomNode::Fun(_, args) => args.iter().any(|a| has_other_var(*a, var)),
                AtomNode::Num(_) => false,
            }
        }
        has_other_var(*b, var)
    })
}

/// MR-P: sites whose base is `(w)²` (syntactically) or the expanded square
/// of an affine `w` (via [`super::affine_square_base`]'s exact
/// re-expansion) rewrite to `w^p·σ` with the sheet factor
/// `σ = w·(w²)^(−1/2)` — formally constant (`Dσ = 0` in the differential
/// field), so the `σ` product factors out of the integral. The rest is
/// integrated by the chain; the answer carries `σ` (the `halfpower` sheet
/// precedent; spec §4.2). The certificate sees `σ` through the radical
/// layer's relation `r² = (w)²` on the original base.
fn try_slice_factor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let mut targets: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    for &base in &coll.bases {
        // Syntactic square `w²` with affine `w`, or an expanded trinomial.
        // Non-square bases are skipped, never abort the pass.
        let w = match base.node() {
            AtomNode::Pow(w, e) if matches!(e.node(), AtomNode::Num(2)) => {
                let Some((slope, _)) = linear_form(ctx, *w, var) else {
                    continue;
                };
                if identically_zero(ctx, slope) {
                    continue;
                }
                *w
            }
            AtomNode::Add(_) => {
                let Some(w) = super::affine_square_base(ctx, base, var) else {
                    continue;
                };
                w
            }
            _ => continue,
        };
        targets.push((base, w));
    }
    if targets.is_empty() {
        return None;
    }
    let mut sigmas: Vec<Atom<'a>> = Vec::new();
    let rest = slice_rewrite(ctx, expr, var, &targets, &mut sigmas)?;
    if sigmas.is_empty() {
        return None;
    }
    // Fast-decline gate (see `has_symbolic_denom_power`): a surviving
    // fractional-power site plus a powered symbolic denominator is the
    // measured 199 s chain burn (rubi-01871).
    if has_symbolic_denom_power(rest, var)
        && let Some(rc) = collect(ctx, rest, var)
        && !rc.sites.is_empty()
    {
        return None;
    }
    if charge_mr_work(node_count(rest) as u64) {
        return None;
    }
    let r = integrate_raw(ctx, rest, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    let mut factors = sigmas;
    factors.push(r);
    Some(normalize(ctx, mul_atoms(ctx, &factors)))
}

// =========================================================================
// P2/P3: collapsible quadratic bases (spec §2.1)
// =========================================================================

/// Re-express radical sites of `cofactor` in the answer through the
/// original base `base = var·cofactor` (`cofactor^(p/q) →
/// base^(p/q)·var^(−p/q)`) — the split's inverse, valid on the
/// gate-validated domain. The certificate certifies `D(F) − f` with the
/// radicals in their original bases (spec §4.1-4); an answer carrying the
/// cofactor base instead would need the relation `r_B = r_x·r_cofactor`,
/// which the radical layer deliberately cannot see (the wrong-merge
/// guard).
fn rejoin_cofactors<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    pairs: &[(Atom<'a>, Atom<'a>)],
) -> Atom<'a> {
    let mut out = expr;
    for &(cofactor, base) in pairs {
        let x = ctx.var(var.as_str());
        let result = rewrite_sites_of_base(ctx, out, var, cofactor, &|ctx, p, q| {
            mul_atoms(
                ctx,
                &[
                    ctx.pow(base, rat_atom(ctx, p, q)),
                    ctx.pow(x, rat_atom(ctx, -p, q)),
                ],
            )
        });
        if let Some(r) = result {
            out = r;
        }
    }
    normalize(ctx, out)
}

/// Rewrite every site whose normalized base is `target` as
/// `∏ factorsᵢ^(p/q)`; `repl` is called with the site's `(p, q)` and
/// returns the replacement.
fn rewrite_sites_of_base<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    target: Atom<'a>,
    repl: &impl Fn(&'a AtomArena<'a>, i64, i64) -> Atom<'a>,
) -> Option<Atom<'a>> {
    if let Some((b, p, q)) = site_of(expr, var)
        && normalize(ctx, b) == target
    {
        return Some(repl(ctx, p, q));
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => Some(expr),
        AtomNode::Add(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| rewrite_sites_of_base(ctx, *a, var, target, repl))
                .collect();
            Some(ctx.add(&kids?))
        }
        AtomNode::Mul(args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| rewrite_sites_of_base(ctx, *a, var, target, repl))
                .collect();
            Some(ctx.mul(&kids?))
        }
        AtomNode::Pow(b, e) => {
            let nb = rewrite_sites_of_base(ctx, *b, var, target, repl)?;
            Some(ctx.pow(nb, *e))
        }
        AtomNode::Fun(name, args) => {
            let kids: Option<Vec<Atom<'a>>> = args
                .iter()
                .map(|a| rewrite_sites_of_base(ctx, *a, var, target, repl))
                .collect();
            Some(ctx.fun(name.as_str(), &kids?))
        }
    }
}

/// P2 (spec §2.1): a quadratic base with zero constant term factors out
/// `x`: `(b·x + c·x²)^(p/q) = x^(p/q)·(b + c·x)^(p/q) · (gate)`. The
/// rewrite needs the numeric evidence gate (the split direction is
/// branch-sensitive); the rewritten integrand re-enters the chain, where
/// the now-affine bases go to MR-L/merge.
fn try_extract_x_factor<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let mut rewritten = expr;
    let mut pairs: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    for &base in &coll.bases {
        let (c2, c1, c0) = super::rules::quadratic_coeffs(ctx, ctx.fun("sqrt", &[base]), var)?;
        if !identically_zero(ctx, c0) || identically_zero(ctx, c1) || identically_zero(ctx, c2) {
            continue;
        }
        // cofactor = c1 + c2·x.
        let cofactor = normalize(ctx, ctx.add(&[c1, ctx.mul(&[c2, x])]));
        rewritten = rewrite_sites_of_base(ctx, rewritten, var, base, &|ctx, p, q| {
            mul_atoms(
                ctx,
                &[
                    ctx.pow(x, rat_atom(ctx, p, q)),
                    ctx.pow(cofactor, rat_atom(ctx, p, q)),
                ],
            )
        })?;
        pairs.push((cofactor, base));
    }
    if pairs.is_empty() {
        return None;
    }
    let rewritten = normalize(ctx, rewritten);
    let rewritten = unify_proportional_bases(ctx, rewritten, var).unwrap_or(rewritten);
    let rewritten = merge_radical_powers(ctx, rewritten, var);
    if !gate_equivalent(expr, rewritten, var) {
        return None;
    }
    if charge_mr_work(node_count(rewritten) as u64) {
        return None;
    }
    let r = integrate_raw(ctx, rewritten, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    // Re-express cofactor radicals in the original bases (spec §4.1-4) so
    // the certificate's radical layer sees only the integrand's bases.
    Some(rejoin_cofactors(ctx, r, var, &pairs))
}

/// P3 (spec §2.1): a quadratic base with integer coefficients and a
/// positive-square discriminant splits into two affine factors over ℚ:
/// `c2·x² + c1·x + c0 = c2·(x − r1)·(x − r2)`. Only this bounded shape (no
/// general factorization call — cost discipline). `c2 > 0` keeps the
/// split's constant factor real; `r1 = r2` is the perfect square (MR-P's
/// case). The rewrite is gated numerically and re-enters the chain.
fn try_factor_quadratic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    coll: &Collect<'a>,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let mut rewritten = expr;
    let mut changed = false;
    for &base in &coll.bases {
        let (c2, c1, c0) = super::rules::quadratic_coeffs(ctx, ctx.fun("sqrt", &[base]), var)?;
        let (AtomNode::Num(a), AtomNode::Num(b), AtomNode::Num(c)) =
            (c2.node(), c1.node(), c0.node())
        else {
            continue;
        };
        let (a, b, c) = (i128::from(*a), i128::from(*b), i128::from(*c));
        if a <= 0 {
            continue;
        }
        let disc = b
            .checked_mul(b)?
            .checked_sub(a.checked_mul(4)?.checked_mul(c)?)?;
        if disc <= 0 {
            continue; // Δ = 0 is MR-P's perfect square; Δ < 0 stays irreducible.
        }
        let s = disc.checked_isqrt()?;
        if s * s != disc {
            continue;
        }
        // r1,2 = (−b ± s)/(2a), distinct since s > 0.
        let den = a.checked_mul(2)?;
        let to_atom = |num: i128| -> Option<Atom<'a>> {
            let g = {
                let mut u = num.unsigned_abs();
                let mut v = den.unsigned_abs();
                while v != 0 {
                    let t = u % v;
                    u = v;
                    v = t;
                }
                u.max(1)
            };
            let (n, d) = (num / g as i128, den / g as i128);
            Some(rat_atom(
                ctx,
                i64::try_from(n).ok()?,
                i64::try_from(d).ok()?,
            ))
        };
        let r1 = to_atom((-b).checked_sub(s)?)?;
        let r2 = to_atom((-b).checked_add(s)?)?;
        let c2a = ctx.num(i64::try_from(a).ok()?);
        // factors f1 = x − r1, f2 = x − r2.
        let f1 = normalize(ctx, ctx.add(&[x, neg(ctx, r1)]));
        let f2 = normalize(ctx, ctx.add(&[x, neg(ctx, r2)]));
        rewritten = rewrite_sites_of_base(ctx, rewritten, var, base, &|ctx, p, q| {
            mul_atoms(
                ctx,
                &[
                    ctx.pow(c2a, rat_atom(ctx, p, q)),
                    ctx.pow(f1, rat_atom(ctx, p, q)),
                    ctx.pow(f2, rat_atom(ctx, p, q)),
                ],
            )
        })?;
        changed = true;
    }
    if !changed {
        return None;
    }
    let rewritten = normalize(ctx, rewritten);
    let rewritten = unify_proportional_bases(ctx, rewritten, var).unwrap_or(rewritten);
    let rewritten = merge_radical_powers(ctx, rewritten, var);
    if !gate_equivalent(expr, rewritten, var) {
        return None;
    }
    if charge_mr_work(node_count(rewritten) as u64) {
        return None;
    }
    let r = integrate_raw(ctx, rewritten, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    Some(normalize(ctx, r))
}

/// Unify proportional affine bases: when `b_i = c·b_j` with `c` a positive
/// rational (decided by a `collect_terms` determinant fold), every site of
/// `b_i` rewrites to `c^(p/q)·b_j^(p/q)` — exact for a positive constant
/// factor. This is what lets a P2/P3 split's monic factors (`x − r`) meet
/// the integrand's other (`c2·(x − r)`-spelled) bases.
fn unify_proportional_bases<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let coll = collect(ctx, expr, var)?;
    let mut rewritten = expr;
    let mut changed = false;
    for i in 0..coll.bases.len() {
        for j in 0..i {
            let (Some((ai, a0)), Some((bj, b0))) = (
                linear_form(ctx, coll.bases[i], var),
                linear_form(ctx, coll.bases[j], var),
            ) else {
                continue;
            };
            // Proportional iff a0·bj − ai·b0 ≡ 0.
            if !identically_zero(
                ctx,
                normalize(
                    ctx,
                    ctx.add(&[ctx.mul(&[a0, bj]), neg(ctx, ctx.mul(&[ai, b0]))]),
                ),
            ) {
                continue;
            }
            // c = ai/bj (the slope ratio); only a positive rational is a
            // branch-safe constant split.
            let c = normalize(ctx, ctx.mul(&[ai, inv(ctx, bj)]));
            let Some((cp, cq)) = rat_of(c) else {
                continue;
            };
            if cp <= 0 {
                continue;
            }
            let c_atom = rat_atom(ctx, cp, cq);
            let target = coll.bases[j];
            rewritten = rewrite_sites_of_base(ctx, rewritten, var, coll.bases[i], &|ctx, p, q| {
                mul_atoms(
                    ctx,
                    &[
                        ctx.pow(c_atom, rat_atom(ctx, p, q)),
                        ctx.pow(target, rat_atom(ctx, p, q)),
                    ],
                )
            })?;
            changed = true;
            break;
        }
    }
    if changed {
        Some(normalize(ctx, rewritten))
    } else {
        None
    }
}

// =========================================================================
// Dispatch
// =========================================================================

/// The base-count reduction loop (spec §2.2): exactly two bases try the
/// merge first, then the Möbius substitution; three or more try one merge
/// round at a time. Every round re-collects (a merge changes the base set).
fn mr_reduce<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    rounds: usize,
    records: &mut Vec<MergeRecord<'a>>,
) -> Option<Atom<'a>> {
    if charge_mr_work(node_count(expr) as u64) {
        return None;
    }
    let coll = collect(ctx, expr, var)?;
    match coll.bases.len() {
        0 | 1 => None,
        2 => try_merge_two(ctx, expr, var, &coll, records)
            .or_else(|| try_mobius(ctx, expr, var, &coll)),
        _ => {
            if rounds >= MAX_MR_ROUNDS {
                return None;
            }
            let merged = try_merge_round(ctx, expr, var, &coll, records)?;
            mr_reduce(ctx, merged, var, rounds + 1, records)
        }
    }
}

/// Integrate `expr` via multi-radical-base reduction (see the module docs).
/// Returns `None` outside the fragment, on any failed precondition, and
/// inside `chain::parts_active()` scopes (spec §3.1 parts discipline).
pub(crate) fn integrate_multi_radical<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if chain::parts_active() {
        return None;
    }
    if is_constant(expr, var) || matches!(expr.node(), AtomNode::Add(_)) {
        return None;
    }
    if node_count(expr) > MAX_MR_NODES {
        return None;
    }
    reset_mr_budget();
    if charge_mr_work(node_count(expr) as u64) {
        return None;
    }
    // P1: same-base powers are summed first (exact in the radical field).
    let merged_input = merge_radical_powers(ctx, expr, var);
    let coll = collect(ctx, merged_input, var)?;
    // Single-base shapes belong to the single-base engines (they ran
    // earlier in the chain); this module only fires on ≥ 2 distinct bases.
    // Exception: a *genuine radical cancellation* (the merge strictly
    // reduced the number of radical sites and < 2 bases remain) re-enters
    // the chain on the collapsed form. A mere respelling must NOT re-enter:
    // the collapse is a full chain sweep, and on non-radical grinders
    // (measured: `sinh(x)^3/(a + b*sinh(x))`, whose inv-hyp t-forms the
    // merge pass respells) it multiplied the chain's work — 130 s → >300 s.
    if coll.bases.len() < 2 {
        if merged_input != expr {
            let before = collect(ctx, expr, var).map(|c| c.sites.len());
            if before.is_some_and(|n| coll.sites.len() < n) {
                let r = integrate_raw(ctx, merged_input, var, 0, true, 0, 0);
                if !contains_integral(r) {
                    return Some(normalize(ctx, r));
                }
            }
        }
        return None;
    }
    let expr = merged_input;
    let mut records = Vec::new();
    // MR-I nested radicals, then the MR-G conjugate sums and MR-R
    // reciprocal bases (their preconditions are disjoint from the affine
    // pair machinery), then the base-count reduction loop (merge → Möbius).
    let r = if coll.sites.iter().any(|s| contains_radical(s.base, var)) {
        try_denest(ctx, expr, var)?
    } else {
        try_conjugate(ctx, expr, var, &coll)
            .or_else(|| try_reciprocal(ctx, expr, var, &coll))
            .or_else(|| try_slice_factor(ctx, expr, var, &coll))
            .or_else(|| try_extract_x_factor(ctx, expr, var, &coll))
            .or_else(|| try_factor_quadratic(ctx, expr, var, &coll))
            .or_else(|| mr_reduce(ctx, expr, var, 0, &mut records))?
    };
    // Express the answer in the original bases (spec §4.1-4), then split
    // constant-radical products for the certificate's polynomial-base
    // fragment.
    let r = split_merged_bases(ctx, r, &records);
    // `(√S)^n` spellings fold to the radicand's own powers (exact).
    let r = fold_sqrt_powers(ctx, r);
    Some(split_constant_radicals(ctx, r, var))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification (covers the heads the
    /// answers carry; `log` goes through the absolute value so sign-only
    /// differences inside log arguments do not matter).
    fn eval_f64_test(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Option<f64> {
        match expr.node() {
            AtomNode::Num(n) => Some(*n as f64),
            AtomNode::Var(v) => env.iter().find(|(s, _)| s == v).map(|(_, val)| *val),
            AtomNode::Add(args) => args
                .iter()
                .try_fold(0.0, |acc, a| Some(acc + eval_f64_test(*a, env)?)),
            AtomNode::Mul(args) => args
                .iter()
                .try_fold(1.0, |acc, a| Some(acc * eval_f64_test(*a, env)?)),
            AtomNode::Pow(b, e) => Some(eval_f64_test(*b, env)?.powf(eval_f64_test(*e, env)?)),
            AtomNode::Fun(name, args) => {
                let v = eval_f64_test(*args.first()?, env)?;
                Some(match name.as_str() {
                    "sin" => v.sin(),
                    "cos" => v.cos(),
                    "exp" => v.exp(),
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

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        normalize(ctx, ocas_parse::parse(ctx, s).expect("parse"))
    }

    /// Run the mechanism directly, require `Some` without residue, check the
    /// symbolic certificate, and verify `diff(result) == integrand`
    /// numerically at the sample points.
    fn assert_solved_certified<'a>(
        ctx: &'a AtomArena<'a>,
        src: &str,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) -> Atom<'a> {
        let var = Symbol::new("x");
        let integrand = parse(ctx, src);
        let result = integrate_multi_radical(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined: {src}"));
        assert!(
            !contains_integral(result),
            "residue in {src} result: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64_test(d, &env).expect("eval diff");
            let rhs = eval_f64_test(integrand, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{src} at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
        let cert = super::super::certify::certify(ctx, integrand, result, var);
        assert!(
            cert.is_ok(),
            "{src}: emitted answer must certify, got {cert:?} for {result}"
        );
        result
    }

    /// Spec §5 L1 (quotient of two affine half powers → Möbius).
    #[test]
    fn l1_quotient_mobius() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "sqrt(1+2*x)/sqrt(3+4*x)", &[], &[0.3, 0.8, 1.5]);
    }

    /// Spec §5 L2 (`x` is itself an affine base).
    #[test]
    fn l2_sqrt_x_times_linear() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "1/(sqrt(x)*sqrt(2+3*x))", &[], &[0.3, 0.9, 1.8]);
    }

    /// Spec §5 L3 (higher half powers; the t-form is a pure rational
    /// function, SymPy declines the original).
    #[test]
    fn l3_mobius_higher_powers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "x^2*sqrt(1+x)/(2+3*x)^(3/2)", &[], &[0.2, 0.9, 2.0]);
    }

    /// Spec §5 F1 (fourth-root pair, `L = 4`).
    #[test]
    fn f1_fourth_root_pair() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(
            &ctx,
            "(1+x)^(1/4)/(1-x)^(5/4)",
            &[],
            &[-0.7, -0.2, 0.4, 0.9],
        );
    }

    /// Symbolic coefficients (spec §8-5): `Δ = a·d − b·c` is decided after
    /// a `collect_terms` fold.
    #[test]
    fn mobius_symbolic_coefficients() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let env = [
            (Symbol::new("a"), 1.0),
            (Symbol::new("b"), 2.0),
            (Symbol::new("c"), 3.0),
            (Symbol::new("d"), 4.0),
        ];
        assert_solved_certified(&ctx, "sqrt(a+b*x)/sqrt(c+d*x)", &env, &[0.3, 0.8, 1.5]);
    }

    /// Proportional bases (`Δ ≡ 0`) are a single-base shape: decline.
    #[test]
    fn proportional_bases_decline() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sqrt(1+2*x)/sqrt(2+4*x)");
        assert!(integrate_multi_radical(&ctx, expr, Symbol::new("x")).is_none());
    }

    /// Spec §5 L4 (rubi-01384 numeric shape): equal-exponent product pair
    /// merges into a single quadratic radicand.
    #[test]
    fn l4_merge_same_slope_pair() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "1/(sqrt(2+3*x)*sqrt(6+3*x))", &[], &[0.1, 0.9, 2.0]);
    }

    /// Spec §5 Q1 (rubi-01373 numeric shape): opposite-translate pair
    /// merges to `√(4x²−1)`.
    #[test]
    fn q1_merge_quadratic_product() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "1/(sqrt(-1+2*x)*sqrt(1+2*x))", &[], &[0.8, 1.3, 2.2]);
    }

    /// Spec §5 Q2 (guard rail): two irreducible quadratic bases — the
    /// Möbius t-form keeps a radical (elliptic class), so the module must
    /// decline honestly.
    #[test]
    fn q2_quadratic_pair_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sqrt(-1+2*x^2)/sqrt(3+4*x^2)");
        assert!(integrate_multi_radical(&ctx, expr, Symbol::new("x")).is_none());
    }

    /// The certificate must reject the pseudo-answer produced by a *wrongly
    /// directed* merge: `D(F)` carries `√(u·v)` while the integrand carries
    /// `√u·√v` — independent generators in the radical field (spec §6-9's
    /// wrong-merge guard).
    #[test]
    fn wrong_merge_certificate_guard() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "sqrt(1+x)*sqrt(1-x)");
        // The true antiderivative of the *merged* radicand `1−x²` …
        let pseudo = ocas_parse::parse(&ctx, "x*sqrt(1-x^2)/2 + asin(x)/2").expect("parse");
        assert!(
            super::super::certify::certify(&ctx, integrand, pseudo, var).is_err(),
            "a wrong-direction merge pseudo-answer must not certify"
        );
    }

    /// Bases containing a transcendental kernel of the variable are outside
    /// the fragment (the W5 interface contract), as are transcendental
    /// factors.
    #[test]
    fn function_heads_decline() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for src in [
            "sqrt(sin(x))*sqrt(1+x)",
            "sqrt(1+x)*sin(x)",
            "sqrt(1+log(x))/sqrt(2+3*x)",
            "x^n*sqrt(1+x)*sqrt(2+x)",
        ] {
            let expr = parse(&ctx, src);
            assert!(
                integrate_multi_radical(&ctx, expr, Symbol::new("x")).is_none(),
                "expected a decline: {src}"
            );
        }
    }

    /// Parts-scope discipline (spec §3.1): inside an integration-by-parts
    /// sub-integral the module declines outright.
    #[test]
    fn parts_scope_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sqrt(1+2*x)/sqrt(3+4*x)");
        let r = chain::with_parts_scope(|| integrate_multi_radical(&ctx, expr, Symbol::new("x")));
        assert!(r.is_none(), "must decline inside a parts scope: {r:?}");
    }

    /// Chain wiring: the pipeline solves the Möbius shapes end to end, and
    /// the Chebyshev single-base shape E1 stays with `binomial` (the
    /// multi-base collector must not preempt it — spec §2.6 E1).
    #[test]
    fn chain_solves_and_does_not_preempt_chebyshev() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for src in ["sqrt(1+2*x)/sqrt(3+4*x)", "x^2*sqrt(1+x)/(2+3*x)^(3/2)"] {
            let expr = ocas_parse::parse(&ctx, src).expect("parse");
            let r = crate::integral::integrate(&ctx, expr, var);
            assert!(
                !r.to_string().contains("Integral"),
                "chain should solve {src}, got: {r}"
            );
        }
        // E1 (rubi-01033): Chebyshev case 3 — `binomial` owns it, and this
        // module declines (its bases are not affine).
        let expr = parse(&ctx, "(1+x^2)^(1/3)/x^(11/3)");
        assert!(integrate_multi_radical(&ctx, expr, var).is_none());
        let r = crate::integral::integrate(
            &ctx,
            ocas_parse::parse(&ctx, "(1+x^2)^(1/3)/x^(11/3)").unwrap(),
            var,
        );
        assert!(
            !r.to_string().contains("Integral"),
            "E1 regressed (binomial path preempted?): {r}"
        );
    }

    /// Run the mechanism directly, require `Some` without residue, and
    /// verify `diff(result) == integrand` numerically. Used where the
    /// certificate surface cannot see the shape (MR-R's reciprocal bases
    /// have rational — not polynomial — relation bases).
    fn assert_solved_numeric<'a>(
        ctx: &'a AtomArena<'a>,
        src: &str,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) -> Atom<'a> {
        let var = Symbol::new("x");
        let integrand = parse(ctx, src);
        let result = integrate_multi_radical(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined: {src}"));
        assert!(
            !contains_integral(result),
            "residue in {src} result: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64_test(d, &env).expect("eval diff");
            let rhs = eval_f64_test(integrand, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{src} at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
        result
    }

    /// The merge pass must see through normalize's `(A³·B³)⁻¹` spelling of
    /// a negative power of a product (the rubi-00287 t-form's shape):
    /// `t⁴·D⁻²·D⁻¹·(N³·D⁻³)⁻¹` merges the `D` powers into a single factor.
    #[test]
    fn merge_sees_through_negative_product_powers() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let t = ctx.var("t");
        let d = normalize(
            &ctx,
            ctx.add(&[
                ctx.mul(&[ctx.num(-1), ctx.var("d")]),
                ctx.mul(&[ctx.var("b"), ctx.pow(t, ctx.num(2))]),
            ]),
        );
        let n = normalize(
            &ctx,
            ctx.add(&[
                ctx.var("c"),
                ctx.mul(&[ctx.num(-1), ctx.var("a"), ctx.pow(t, ctx.num(2))]),
            ]),
        );
        let f = ctx.mul(&[
            int_pow(&ctx, t, 4),
            int_pow(&ctx, d, -2),
            inv(&ctx, d),
            ctx.pow(
                normalize(&ctx, ctx.mul(&[int_pow(&ctx, n, 3), int_pow(&ctx, d, -3)])),
                ctx.num(-1),
            ),
        ]);
        let m = merge_radical_powers(&ctx, f, Symbol::new("t"));
        // All D powers collapse: D⁻²·D⁻¹·(D⁻³)⁻¹ = 1 — only t⁴ and N⁻³
        // survive (plus the dropped unit power).
        assert_eq!(m.to_string(), "(t^4)*((c + (-1*a*(t^2)))^-3)");
    }

    /// Spec §5 G1 (rubi-00282 numeric shape): denominator conjugate, `u − v`
    /// constant.
    #[test]
    fn g1_conjugate_denominator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(
            &ctx,
            "x^2/(sqrt(1+2*x)+sqrt(3+2*x))^2",
            &[],
            &[0.3, 0.9, 1.8],
        );
    }

    /// Spec §5 G2 (rubi-00343): numerator square of a radical sum.
    #[test]
    fn g2_numerator_square() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "x^3*(sqrt(1-x)+sqrt(1+x))^2", &[], &[0.2, 0.5, 0.8]);
    }

    /// Spec §5 R1 (rubi-00870 numeric shape): reciprocal bases via
    /// `u = 1/x`, then MR-L. The certificate cannot see reciprocal bases
    /// (rational relation bases), so this test is numeric-only.
    #[test]
    fn r1_reciprocal_bases() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_numeric(&ctx, "sqrt(1+2/x)/(3+4/x)^(3/2)", &[], &[0.3, 0.9, 1.8]);
    }

    /// Spec §5 I1 (rubi-01039): one denesting round (`x = t²`) leaves a
    /// single quadratic base. Numeric-only: the answer's `asin` form
    /// differentiates to `sqrt((4/5)·S)`-spelled terms, a constant multiple
    /// of the other radical base — the certificate's radical layer cannot
    /// equate constant multiples of generators (the same reason
    /// `sqrt_quadratic`'s own asin-branch tests are numeric-only).
    #[test]
    fn i1_nested_sqrt_linear() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_numeric(&ctx, "sqrt(1-x-sqrt(x))", &[], &[0.05, 0.15, 0.3]);
    }

    /// Spec §5 I2 (rubi-01588): two denesting rounds leave a single linear
    /// base (`binomial` territory).
    #[test]
    fn i2_double_nested() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_numeric(&ctx, "1/sqrt(2+sqrt(1+sqrt(x)))", &[], &[0.3, 0.8, 1.5]);
    }

    /// Nested radicals with a non-affine inner base decline honestly.
    #[test]
    fn nested_non_affine_inner_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "sqrt(1-x^2-sqrt(1+x^2))");
        assert!(integrate_multi_radical(&ctx, expr, Symbol::new("x")).is_none());
    }

    /// Spec §5 P1 (rubi-01701 numeric shape): the sheet factor
    /// `σ = (1+x)·((1+x)²)^(−1/2)` factors out; samples on both sheets
    /// (`x > −1` and `x < −1`).
    #[test]
    fn p1_slice_factor() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(
            &ctx,
            "sqrt((1+x)^2)/(2+x)^(7/2)",
            &[],
            &[-1.7, -0.5, 0.3, 1.2],
        );
    }

    /// P1 with the square written as the expanded trinomial (the
    /// `affine_square_base` path).
    #[test]
    fn p1_expanded_trinomial() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(
            &ctx,
            "sqrt(1+2*x+x^2)/(2+x)^(7/2)",
            &[],
            &[-1.7, -0.5, 0.3, 1.2],
        );
    }

    /// The MR-P fast-decline gate: a surviving fractional-power site plus a
    /// powered symbolic denominator (the `symbolic_rational` weak spot) is
    /// declined instead of burning the chain budget (measured: rubi-01871,
    /// 199 s without the gate).
    #[test]
    fn slice_with_symbolic_denom_power_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, "(d + e*x)^(5/2)/(a^2 + 2*a*b*x + b^2*x^2)^(3/2)");
        assert!(integrate_multi_radical(&ctx, expr, Symbol::new("x")).is_none());
    }

    /// Spec §5 W1 (rubi-00537 numeric shape, P2): the zero-constant-term
    /// quadratic base splits off `x`, collapsing to a single-base shape.
    #[test]
    fn p2_shared_x_factor() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_certified(&ctx, "x^(3/2)/sqrt(2*x+3*x^2)", &[], &[0.3, 0.9, 1.8]);
    }

    /// P3 (rubi-00247 family): the factorable quadratic base
    /// `2+5x+3x² = (1+x)(2+3x)` splits; the affine factors cancel against
    /// the other bases, leaving a rational integrand. Numeric-only: the
    /// answer is rational while the integrand carries the factorization
    /// relation `r_c = r_a·r_b` between its radical generators, which the
    /// certificate deliberately cannot see (the wrong-merge guard) — the
    /// numeric evidence gate carries the split.
    #[test]
    fn p3_factorable_quadratic() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        assert_solved_numeric(
            &ctx,
            "sqrt(2+3*x)*sqrt(1+x)/(2+5*x+3*x^2)^(3/2)",
            &[],
            &[0.3, 0.9, 1.8],
        );
    }
}
