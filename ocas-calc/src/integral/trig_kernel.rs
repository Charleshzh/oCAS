//! Single-trig-kernel rational forms and tan/sec-family reductions (0.27.1).
//!
//! Two mechanisms over a single trigonometric kernel `T(u)` with
//! `u = e + f·x` linear (symbolic coefficients allowed):
//!
//! - **K1 [`integrate_trig_kernel`] (single-kernel rational)**:
//!   `P(T)/(a + b·T)^n` where `P` is a polynomial ( Laurent after
//!   `sec`/`csc` conversion) in `T ∈ {sin, cos}`. Positive-degree terms use
//!   the shift `S = a + b·T` (so `T = (S − a)/b`) and binomial expansion;
//!   negative-degree terms use the 2-D peel
//!   `1/(T^s Q^n) = (1/a)·1/(T^s Q^{n−1}) − (b/a)·1/(T^{s−1} Q^n)`, which
//!   terminates at `1/Q^j` (trig_reduction T1) and `1/T^i` (rules C13/C14).
//!   The reduced sum is re-integrated through the chain once and declined
//!   on any `Integral` residue.
//! - **K2 tan^m·sec^n / cot^m·csc^n reductions**: classic even-`n` peel of
//!   `sec²` / odd-`m` peel of `tan·sec` (and the cot/csc mirror), emitted
//!   as closed forms; the even-`m`/odd-`n` case is rewritten to a sum of
//!   pure `sec^j`/`csc^j` powers and re-entered through the chain.
//!
//! Fractional kernel powers (the elliptic family, e.g. `1/√(a+b·sec)`)
//! decline honestly. The chain-order guard: K1 declines pure
//! `const·Q^(−n)` (all numerator degrees zero) — that is T1's shape — so
//! reduced terms re-entering the chain never re-match here.

use std::collections::BTreeMap;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::{
    contains_integral, int_pow, integrate_raw, inv, is_constant, linear_form, node_count, rat_atom,
};

/// Cap on kernel powers, polynomial degrees and the denominator power.
const MAX_POW: i64 = 8;
/// Node budget for the input integrand.
const MAX_NODES: usize = 300;
/// Node budget for the reduced sums handed back to the chain.
const MAX_REENTRY_NODES: usize = 20_000;

/// Work units a single entry into this module may charge.
///
/// The K1 peel recurses as `peel(s, n) → peel(s, n−1) + peel(s−1, n)`,
/// which has `C(s+n, s)` leaves — already bounded by `2^(2·MAX_POW)` calls
/// for the accepted shapes, but each leaf builds atoms, so the enumeration
/// is charged too. The cap is a backstop far above the worst legitimate
/// enumeration (the `sec(u)^n/(a+b·sec(u))^m` family): tripping it would
/// make K1 decline and hand the shape to the K2 rewrite, which returns a
/// different — still correct, but much larger — closed form.
const MAX_KERNEL_WORK: u64 = 5_000_000;

thread_local! {
    static KERNEL_WORK: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Reset the work budget; called by the public stage entry point so a
/// single invocation cannot accumulate a stale budget.
fn reset_kernel_budget() {
    KERNEL_WORK.with(|c| c.set(0));
}

/// Charge `units` against the module budget; `true` when exhausted.
fn charge_kernel_work(units: u64) -> bool {
    KERNEL_WORK.with(|c| {
        let v = c.get().saturating_add(units);
        c.set(v);
        v > MAX_KERNEL_WORK
    })
}

/// Integrate `expr` via the K1/K2 mechanisms (see module docs).
pub(crate) fn integrate_trig_kernel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    reset_kernel_budget();
    if is_constant(expr, var) || node_count(expr) > MAX_NODES {
        return None;
    }
    try_single_kernel(ctx, expr, var)
        .or_else(|| try_tan_sec_family(ctx, expr, var, Family::TanSec))
        .or_else(|| try_tan_sec_family(ctx, expr, var, Family::CotCsc))
}

// =========================================================================
// Shared matching
// =========================================================================

/// Trig function conversion to sin/cos exponents: returns `None` for
/// non-trig names. `(sin_exp, cos_exp)` per unit power of the named fun.
fn trig_pair(name: &str) -> Option<(i64, i64)> {
    match name {
        "sin" => Some((1, 0)),
        "cos" => Some((0, 1)),
        "sec" => Some((0, -1)),
        "csc" => Some((-1, 0)),
        "tan" => Some((1, -1)),
        "cot" => Some((-1, 1)),
        _ => None,
    }
}

/// One trig-kernel occurrence: function name, integer exponent, argument.
struct TrigPow<'a> {
    name: &'static str,
    exp: i64,
    arg: Atom<'a>,
}

/// Match `f` as `T(u)^k` (bare `T(u)` reads as k = 1). Integer |k| ≤
/// MAX_POW only; fractional powers decline.
fn match_trig_pow<'a>(f: Atom<'a>) -> Option<TrigPow<'a>> {
    let (name, args, exp) = match f.node() {
        AtomNode::Fun(name, args) if args.len() == 1 => (name.as_str(), args, 1),
        AtomNode::Pow(b, e) => match (b.node(), e.node()) {
            (AtomNode::Fun(name, args), AtomNode::Num(k)) if args.len() == 1 => {
                (name.as_str(), args, *k)
            }
            _ => return None,
        },
        _ => return None,
    };
    if trig_pair(name).is_none() || exp.abs() > MAX_POW {
        return None;
    }
    // Leak the &str to 'static is impossible; re-map below.
    Some(TrigPow {
        name: match name {
            "sin" => "sin",
            "cos" => "cos",
            "tan" => "tan",
            "cot" => "cot",
            "sec" => "sec",
            "csc" => "csc",
            _ => return None,
        },
        exp,
        arg: args[0],
    })
}

/// Consistent linear argument for all kernel factors: returns the shared
/// `u` and its slope `f` in `u = e + f·x`.
fn shared_linear_arg<'a>(
    ctx: &'a AtomArena<'a>,
    u: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    let (f, _e) = linear_form(ctx, u, var)?;
    if is_zero(ctx, f) {
        return None;
    }
    Some((u, f))
}

/// Unify a kernel argument with the shared one (normalized comparison).
fn unify_u<'a>(
    ctx: &'a AtomArena<'a>,
    slot: &mut Option<Atom<'a>>,
    u: Atom<'a>,
) -> Option<Atom<'a>> {
    let nu = normalize(ctx, u);
    match *slot {
        None => {
            *slot = Some(nu);
            Some(nu)
        }
        Some(u0) => {
            if nu == u0 {
                Some(u0)
            } else {
                None
            }
        }
    }
}

fn is_zero<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> bool {
    matches!(normalize(ctx, a).node(), AtomNode::Num(0))
}

/// Small binomial coefficient (checked; declines on overflow via caller).
fn binom(n: i64, k: i64) -> i64 {
    let mut r = 1i64;
    for i in 0..k {
        r = r.saturating_mul(n - i) / (i + 1);
    }
    r
}

// =========================================================================
// K1: single-kernel rational normalization
// =========================================================================

/// Match result for `P(T(u)) / (a + b·T(u))^n` over one kernel.
struct K1Match<'a> {
    /// true for a sin kernel, false for cos.
    sin: bool,
    /// The shared linear argument.
    u: Atom<'a>,
    /// Constant factors.
    consts: Vec<Atom<'a>>,
    /// Numerator Laurent terms: coefficient (constant atoms) × T-degree.
    numer: Vec<(Atom<'a>, i64)>,
    /// Denominator `(a + b·T)^n`, n ≥ 1, b ≠ 0.
    da: Atom<'a>,
    db: Atom<'a>,
    dn: i64,
}

/// Analyze an `Add` as a Laurent binomial `T^d0·(a' + b'·T)` in the given
/// kernel: returns `(d0, a', b')`. Every term must be `const·T^d` with the
/// same kernel; the distinct degrees must be exactly `{d0, d0 + 1}`.
fn add_as_kernel_binomial<'a>(
    ctx: &'a AtomArena<'a>,
    args: &[Atom<'a>],
    sin: bool,
    u: Atom<'a>,
    var: Symbol,
) -> Option<(i64, Atom<'a>, Atom<'a>)> {
    let mut terms: BTreeMap<i64, Vec<Atom<'a>>> = BTreeMap::new();
    for t in args {
        let (coeff, deg) = term_as_kernel_power(ctx, *t, sin, u, var)?;
        terms.entry(deg).or_default().push(coeff);
    }
    let degrees: Vec<i64> = terms.keys().copied().collect();
    if degrees.len() != 2 || degrees[1] != degrees[0] + 1 {
        return None;
    }
    let d0 = degrees[0];
    let a = normalize(ctx, ctx.add(&terms[&d0]));
    let b = normalize(ctx, ctx.add(&terms[&(d0 + 1)]));
    if is_zero(ctx, a) || is_zero(ctx, b) {
        return None;
    }
    Some((d0, a, b))
}

/// Match a single `Add` term (or any factor) as `const·T(same u)^deg`.
fn term_as_kernel_power<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    sin: bool,
    u: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, i64)> {
    if is_constant(term, var) {
        return Some((term, 0));
    }
    let sub: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };
    let mut coeff: Vec<Atom<'a>> = Vec::new();
    let mut deg: Option<i64> = None;
    for f in sub {
        if is_constant(f, var) {
            coeff.push(f);
            continue;
        }
        let tp = match_trig_pow(f)?;
        if deg.is_some() {
            return None;
        }
        if normalize(ctx, tp.arg) != normalize(ctx, u) {
            return None;
        }
        let (se, ce) = trig_pair(tp.name)?;
        // The term must stay within the single kernel.
        let d = if sin { se } else { ce };
        let other = if sin { ce } else { se };
        if other != 0 {
            return None;
        }
        deg = Some(d.checked_mul(tp.exp)?);
    }
    let coeff = if coeff.is_empty() {
        ctx.num(1)
    } else {
        normalize(ctx, ctx.mul(&coeff))
    };
    Some((coeff, deg?))
}

fn try_single_kernel<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
    let m = match_k1(ctx, expr, var)?;
    // Guard: pure `const·Q^(−n)` (all numerator degrees zero) is the
    // trig_reduction T1 shape — declining it keeps K1 off its own re-entry
    // terminals.
    if m.numer.iter().all(|(_, d)| *d == 0) {
        return None;
    }
    let t_atom = ctx.fun(if m.sin { "sin" } else { "cos" }, &[m.u]);
    let q_atom = normalize(ctx, ctx.add(&[m.da, ctx.mul(&[m.db, t_atom])]));
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for (c, d) in &m.numer {
        if *d >= 0 {
            shift_expand(ctx, *c, *d, m.dn, m.da, m.db, t_atom, q_atom, &mut terms)?;
        } else {
            peel(ctx, *c, -*d, m.dn, m.da, m.db, t_atom, q_atom, &mut terms)?;
        }
    }
    if terms.is_empty() {
        return None;
    }
    let mut all = m.consts.clone();
    all.push(ctx.add(&terms));
    let sum = normalize(ctx, ctx.mul(&all));
    if sum == expr {
        return None;
    }
    // Size gate before the chain re-entry: the reduced sum can be much
    // larger than the integrand, and the re-entered chain pays for every
    // node.
    if node_count(sum) > MAX_REENTRY_NODES {
        return None;
    }
    let r = integrate_raw(ctx, sum, var, 0, true, 0, 0);
    if contains_integral(r) {
        return None;
    }
    Some(r)
}

/// Shift-expand `c·T^m/Q^n` with `m ≥ 0`: `T = (S − a)/b`, `S = Q`; each
/// `S^k·Q^(−n)` lands as `(a+bT)^(k−n)` (re-expanded in T) or `Q^(k−n)`.
#[allow(clippy::too_many_arguments)]
fn shift_expand<'a>(
    ctx: &'a AtomArena<'a>,
    c: Atom<'a>,
    m: i64,
    n: i64,
    a: Atom<'a>,
    b: Atom<'a>,
    t_atom: Atom<'a>,
    q_atom: Atom<'a>,
    terms: &mut Vec<Atom<'a>>,
) -> Option<()> {
    let b_inv_m = ctx.pow(b, ctx.num(-m));
    for k in 0..=m {
        if charge_kernel_work(1 + terms.len() as u64) {
            return None;
        }
        // C(m,k)·(−a)^(m−k)·b^(−m)·S^k·Q^(−n)
        let mut coeff = vec![c, ctx.num(binom(m, k)), int_pow(ctx, a, m - k), b_inv_m];
        if (m - k) % 2 == 1 {
            coeff.push(ctx.num(-1));
        }
        let coeff = normalize(ctx, ctx.mul(&coeff));
        let j = k - n;
        if j >= 0 {
            // (a + b·T)^j re-expanded in T.
            for l in 0..=j {
                let c3 = normalize(
                    ctx,
                    ctx.mul(&[
                        coeff,
                        ctx.num(binom(j, l)),
                        int_pow(ctx, a, j - l),
                        int_pow(ctx, b, l),
                    ]),
                );
                terms.push(normalize(ctx, ctx.mul(&[c3, int_pow(ctx, t_atom, l)])));
            }
        } else {
            terms.push(normalize(
                ctx,
                ctx.mul(&[coeff, ctx.pow(q_atom, ctx.num(j))]),
            ));
        }
    }
    Some(())
}

/// Peel `c·T^(−s)·Q^(−n)` (s, n ≥ 1) via
/// `1/(T^s Q^n) = (1/a)/(T^s Q^(n−1)) − (b/a)/(T^(s−1) Q^n)`, terminating
/// at `T^(−s)` (sec/csc powers) and `Q^(−n)` (T1 shapes).
///
/// `T^(−s)` is emitted as the NAMED reciprocal function (`sec^s` /
/// `csc^s`): the rules table owns those, while a bare `cos(u)^(−1)` power
/// has no rule coverage and would fall back.
#[allow(clippy::too_many_arguments)]
fn peel<'a>(
    ctx: &'a AtomArena<'a>,
    c: Atom<'a>,
    s: i64,
    n: i64,
    a: Atom<'a>,
    b: Atom<'a>,
    t_atom: Atom<'a>,
    q_atom: Atom<'a>,
    terms: &mut Vec<Atom<'a>>,
) -> Option<()> {
    if s + n > 2 * MAX_POW {
        return None;
    }
    // The recursion below branches into `peel(s, n−1)` and `peel(s−1, n)`,
    // so one call enumerates `C(s+n, s)` leaves. Charge every call (and the
    // emitted terms) so the enumeration cannot run away.
    if charge_kernel_work(1 + terms.len() as u64) {
        return None;
    }
    if n == 0 {
        // Reciprocal kernel power → named sec/csc (see fn docs).
        let (name, u) = match t_atom.node() {
            AtomNode::Fun(name, args) => (name.as_str(), args[0]),
            _ => return None,
        };
        let recip = ctx.fun(if name == "cos" { "sec" } else { "csc" }, &[u]);
        terms.push(normalize(ctx, ctx.mul(&[c, int_pow(ctx, recip, s)])));
        return Some(());
    }
    if s == 0 {
        terms.push(normalize(ctx, ctx.mul(&[c, ctx.pow(q_atom, ctx.num(-n))])));
        return Some(());
    }
    let a_inv = inv(ctx, a);
    peel(
        ctx,
        normalize(ctx, ctx.mul(&[c, a_inv])),
        s,
        n - 1,
        a,
        b,
        t_atom,
        q_atom,
        terms,
    )?;
    peel(
        ctx,
        normalize(ctx, ctx.mul(&[c, ctx.num(-1), b, a_inv])),
        s - 1,
        n,
        a,
        b,
        t_atom,
        q_atom,
        terms,
    )
}

/// Match the K1 shape; see [`K1Match`].
///
/// Numerator structure: pure kernel powers and the denominator's clearing
/// shift are MULTIPLICATIVE (their degrees add), while an `Add` numerator
/// or a positive binomial-power factor contributes ADDITIVE Laurent terms.
/// The two combine as `(Σ terms)·T^shift`.
fn match_k1<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<K1Match<'a>> {
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    // Additive numerator Laurent terms (at most one structural source).
    let mut numer_add: Vec<(Atom<'a>, i64)> = Vec::new();
    // Multiplicative kernel-degree shift (pure powers + denominator clear).
    let mut kernel_shift: i64 = 0;
    let mut denom: Option<(Atom<'a>, Atom<'a>, i64)> = None;
    // Kernel usage flags and the shared argument.
    let mut use_sin = false;
    let mut use_cos = false;
    let mut shared_u: Option<Atom<'a>> = None;

    for f in factors {
        if is_constant(f, var) {
            consts.push(f);
            continue;
        }
        match f.node() {
            // (Laurent binomial in one kernel)^n, n ≠ 0 integer.
            AtomNode::Pow(b, e) if matches!(b.node(), AtomNode::Add(_)) => {
                let AtomNode::Num(n) = e.node() else {
                    return None;
                };
                let n = *n;
                if n == 0 || n.abs() > MAX_POW {
                    return None;
                }
                // Determine the kernel of this Add from its non-constant
                // terms before conversion: probe with both kernels.
                let AtomNode::Add(add_args) = b.node() else {
                    return None;
                };
                // Pick the kernel by the first trig term's name.
                let mut probe_sin = None;
                for t in add_args.iter() {
                    if is_constant(*t, var) {
                        continue;
                    }
                    let sub: Vec<Atom<'a>> = match t.node() {
                        AtomNode::Mul(m) => m.to_vec(),
                        _ => vec![*t],
                    };
                    for sf in sub {
                        if is_constant(sf, var) {
                            continue;
                        }
                        let tp = match_trig_pow(sf)?;
                        let (se, ce) = trig_pair(tp.name)?;
                        if se != 0 && ce != 0 {
                            // tan/cot mix kernels — K1 declines.
                            return None;
                        }
                        let is_sin = se != 0;
                        if let Some(prev) = probe_sin {
                            if prev != is_sin {
                                return None;
                            }
                        } else {
                            probe_sin = Some(is_sin);
                        }
                        unify_u(ctx, &mut shared_u, tp.arg)?;
                        if is_sin {
                            use_sin = true;
                        } else {
                            use_cos = true;
                        }
                    }
                }
                let sin = probe_sin?;
                let u = shared_u?;
                let (d0, a, bb) = add_as_kernel_binomial(ctx, add_args, sin, u, var)?;
                if n < 0 {
                    if denom.is_some() {
                        return None;
                    }
                    denom = Some((a, bb, -n));
                    // F = (T^d0·Q)^(−dn) = T^(d0·n)·Q^(−dn) (n < 0): the
                    // factored T^d0 power multiplies the numerator.
                    kernel_shift = kernel_shift.checked_add(d0.checked_mul(n)?)?;
                } else {
                    // Numerator binomial power: expand (T^d0·(a + b·T))^n.
                    // Multiplicative with any Add numerator — decline the
                    // combination rather than convolving.
                    if !numer_add.is_empty() {
                        return None;
                    }
                    for k in 0..=n {
                        let coeff = normalize(
                            ctx,
                            ctx.mul(&[
                                ctx.num(binom(n, k)),
                                int_pow(ctx, a, n - k),
                                int_pow(ctx, bb, k),
                            ]),
                        );
                        let deg = d0.checked_mul(n)?.checked_add(k)?;
                        numer_add.push((coeff, deg));
                    }
                }
            }
            // Numerator polynomial (Add of const·T^k terms).
            AtomNode::Add(add_args) => {
                // Kernel determined by the first trig term; constants pass.
                let mut this_sin = None;
                for t in add_args.iter() {
                    if is_constant(*t, var) {
                        continue;
                    }
                    let sub: Vec<Atom<'a>> = match t.node() {
                        AtomNode::Mul(m) => m.to_vec(),
                        _ => vec![*t],
                    };
                    for sf in &sub {
                        if is_constant(*sf, var) {
                            continue;
                        }
                        let tp = match_trig_pow(*sf)?;
                        let (se, ce) = trig_pair(tp.name)?;
                        if se != 0 && ce != 0 {
                            return None;
                        }
                        let is_sin = se != 0;
                        if let Some(prev) = this_sin {
                            if prev != is_sin {
                                return None;
                            }
                        } else {
                            this_sin = Some(is_sin);
                        }
                        unify_u(ctx, &mut shared_u, tp.arg)?;
                    }
                }
                let sin = this_sin?;
                if sin {
                    use_sin = true;
                } else {
                    use_cos = true;
                }
                let u = shared_u?;
                if !numer_add.is_empty() {
                    // An Add numerator combined with a numerator binomial
                    // power would need a convolution — decline.
                    return None;
                }
                for t in add_args.iter() {
                    let (coeff, deg) = term_as_kernel_power(ctx, *t, sin, u, var)?;
                    numer_add.push((coeff, deg));
                }
            }
            _ => {
                // Pure kernel power factor (or decline): multiplicative.
                let tp = match_trig_pow(f)?;
                let (se, ce) = trig_pair(tp.name)?;
                if se != 0 && ce != 0 {
                    return None;
                }
                let is_sin = se != 0;
                if is_sin {
                    use_sin = true;
                } else {
                    use_cos = true;
                }
                unify_u(ctx, &mut shared_u, tp.arg)?;
                let deg = if is_sin { se } else { ce }.checked_mul(tp.exp)?;
                kernel_shift = kernel_shift.checked_add(deg)?;
            }
        }
    }
    if use_sin == use_cos {
        // Both used (mixed kernels) or none at all.
        return None;
    }
    let (da, db, dn) = denom?;
    let u = shared_u?;
    shared_linear_arg(ctx, u, var)?;
    // Apply the multiplicative shift to the additive numerator terms (a
    // bare power numerator is a single term of degree 0 before shifting).
    if numer_add.is_empty() {
        numer_add.push((ctx.num(1), 0));
    }
    let shifted: Vec<(Atom<'a>, i64)> = numer_add
        .into_iter()
        .map(|(c, d)| Some((c, d.checked_add(kernel_shift)?)))
        .collect::<Option<_>>()?;
    // Merge equal degrees and fold coefficients.
    let mut merged: BTreeMap<i64, Vec<Atom<'a>>> = BTreeMap::new();
    for (c, d) in shifted {
        merged.entry(d).or_default().push(c);
    }
    let numer: Vec<(Atom<'a>, i64)> = merged
        .into_iter()
        .map(|(d, cs)| (normalize(ctx, ctx.add(&cs)), d))
        .collect();
    Some(K1Match {
        sin: use_sin,
        u,
        consts,
        numer,
        da,
        db,
        dn,
    })
}

// =========================================================================
// K2: tan^m·sec^n / cot^m·csc^n
// =========================================================================

#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    TanSec,
    CotCsc,
}

impl Family {
    /// (tan-like fun, sec-like fun)
    fn names(self) -> (&'static str, &'static str) {
        match self {
            Family::TanSec => ("tan", "sec"),
            Family::CotCsc => ("cot", "csc"),
        }
    }
}

fn try_tan_sec_family<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fam: Family,
) -> Option<Atom<'a>> {
    let (tname, sname) = fam.names();
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut m = 0i64; // tan/cot exponent
    let mut n = 0i64; // sec/csc exponent
    let mut shared_u: Option<Atom<'a>> = None;
    for f in factors {
        if is_constant(f, var) {
            consts.push(f);
            continue;
        }
        let tp = match_trig_pow(f)?;
        if tp.name != tname && tp.name != sname {
            return None;
        }
        let nu = normalize(ctx, tp.arg);
        match shared_u {
            None => shared_u = Some(nu),
            Some(u0) => {
                if nu != u0 {
                    return None;
                }
            }
        }
        if tp.exp < 0 {
            // Negative powers are the K1/rules territory.
            return None;
        }
        if tp.name == tname {
            m = m.checked_add(tp.exp)?;
        } else {
            n = n.checked_add(tp.exp)?;
        }
    }
    let u = shared_u?;
    let (_, slope) = shared_linear_arg(ctx, u, var)?;
    if m + n < 2 || m > MAX_POW || n > MAX_POW {
        return None;
    }
    if m == 0 || n == 0 {
        // Pure tan^m / sec^n: the rules table owns these.
        return None;
    }
    let tan = ctx.fun(tname, &[u]);
    let sec = ctx.fun(sname, &[u]);
    let f_inv = inv(ctx, slope);
    // Sign of the antiderivative of the peeled derivative factor:
    // tan' = +sec², sec' = +sec·tan; cot' = −csc², csc' = −csc·cot.
    let sign: i64 = match fam {
        Family::TanSec => 1,
        Family::CotCsc => -1,
    };
    let mut terms: Vec<Atom<'a>> = Vec::new();
    if n >= 2 && n % 2 == 0 {
        // Peel sec²: sec^(n−2) = (1 + tan²)^((n−2)/2) (csc² = 1 + cot²).
        let h = (n - 2) / 2;
        for k in 0..=h {
            let j = m + 2 * k;
            // ± C(h,k)·tan^(j+1)/((j+1)·f)
            let c0 = normalize(
                ctx,
                ctx.mul(&[ctx.num(binom(h, k) * sign), f_inv, rat_atom(ctx, 1, j + 1)]),
            );
            terms.push(normalize(ctx, ctx.mul(&[c0, int_pow(ctx, tan, j + 1)])));
        }
    } else if m % 2 == 1 && n >= 1 {
        // Peel tan·sec: tan^(m−1) = (sec²−1)^((m−1)/2).
        let h = (m - 1) / 2;
        for k in 0..=h {
            let s2: i64 = if (h - k) % 2 == 0 { 1 } else { -1 };
            let j = n - 1 + 2 * k;
            // ∫ sec^j·(sec·tan) du = sec^(j+1)/(j+1) (j ≥ 1); j = 0 → sec.
            let body = if j == 0 {
                normalize(ctx, ctx.mul(&[f_inv, sec]))
            } else {
                normalize(
                    ctx,
                    ctx.mul(&[f_inv, rat_atom(ctx, 1, j + 1), int_pow(ctx, sec, j + 1)]),
                )
            };
            terms.push(normalize(
                ctx,
                ctx.mul(&[ctx.num(binom(h, k) * s2 * sign), body]),
            ));
        }
    } else if m >= 2 && m % 2 == 0 && n % 2 == 1 {
        // tan^m = (sec² − 1)^(m/2): rewrite to a sum of pure sec powers
        // and re-enter the chain (rules C13/C14 own sec^j).
        let h = m / 2;
        for k in 0..=h {
            let s2: i64 = if (h - k) % 2 == 0 { 1 } else { -1 };
            let j = n + 2 * k;
            terms.push(normalize(
                ctx,
                ctx.mul(&[ctx.num(binom(h, k) * s2), int_pow(ctx, sec, j)]),
            ));
        }
        let mut all = consts.clone();
        all.push(ctx.add(&terms));
        let sum = normalize(ctx, ctx.mul(&all));
        if sum == expr {
            return None;
        }
        // Size gate before the chain re-entry (see `try_single_kernel`).
        if node_count(sum) > MAX_REENTRY_NODES {
            return None;
        }
        let r = integrate_raw(ctx, sum, var, 0, true, 0, 0);
        if contains_integral(r) {
            return None;
        }
        return Some(r);
    } else {
        return None;
    }
    if terms.is_empty() {
        return None;
    }
    let mut all = consts;
    all.push(ctx.add(&terms));
    Some(normalize(ctx, ctx.mul(&all)))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integral::integrate;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification.
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
                    "sec" => v.cos().recip(),
                    "csc" => v.sin().recip(),
                    "cot" => v.tan().recip(),
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
                    "asinh" => v.asinh(),
                    _ => return None,
                })
            }
        }
    }

    fn parse_norm<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        let e = ocas_parse::parse(ctx, s).unwrap();
        normalize(ctx, e)
    }

    /// Module-entry numeric check (the chain may route elsewhere, so call
    /// the module entry directly).
    fn assert_module_solves(input: &str, consts: &[(Symbol, f64)], samples: &[f64]) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, input);
        let result = integrate_trig_kernel(&ctx, expr, Symbol::new("x"))
            .unwrap_or_else(|| panic!("declined: {input}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue for {input}: {result}"
        );
        let d = crate::diff(&ctx, result, Symbol::new("x"));
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((Symbol::new("x"), xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(expr, &env).expect("eval integrand");
            let tol = 1e-5 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "{input} at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
    }

    /// Structural check for shapes whose emitted antiderivative is
    /// complex-branched: `sec(u)^n/(a + b·sec(u))^m` comes with radicals of
    /// `b² − a²` and an `atanh` of an argument outside `(−1, 1)` for one sign
    /// of `b² − a²`, so a real-valued elementwise evaluator reads `NaN` at
    /// *every* constant choice. The mechanism must still fire and produce a
    /// definite (residue-free) closed form.
    fn assert_module_no_residue(input: &str) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, input);
        let result = integrate_trig_kernel(&ctx, expr, Symbol::new("x"))
            .unwrap_or_else(|| panic!("declined: {input}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue for {input}: {result}"
        );
    }

    fn assert_module_declined(input: &str) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, input);
        assert!(
            integrate_trig_kernel(&ctx, expr, Symbol::new("x")).is_none(),
            "expected None for {input}"
        );
    }

    // ------------------------- K1 -------------------------

    #[test]
    fn k1_sec_over_linear_sec_squared() {
        // Corpus shape: sec(u)/(a + a·sec(u))² → cos/(a·cos + a)² (T2 shape
        // after conversion).
        let env = [
            (Symbol::new("a"), 2.0),
            (Symbol::new("c"), 0.3),
            (Symbol::new("d"), 0.7),
        ];
        assert_module_solves(
            "sec(c + d*x)/(a + a*sec(c + d*x))^2",
            &env,
            &[0.2, 0.5, 0.9],
        );
    }

    #[test]
    fn k1_linear_numerator_over_sin() {
        // (3 − 3·sin(2x))²/(1 + sin(2x)) — numeric variant of the corpus
        // (c − c·sin)⁴/(a + a·sin) shape.
        assert_module_solves("(3 - 3*sin(2*x))^2/(1 + sin(2*x))", &[], &[0.2, 0.5, 0.9]);
    }

    #[test]
    fn k1_corpus_sin_fourth_power() {
        // The corpus shape itself (a, c symbolic).
        let env = [
            (Symbol::new("a"), 2.0),
            (Symbol::new("c"), 1.5),
            (Symbol::new("e"), 0.4),
            (Symbol::new("f"), 0.8),
        ];
        assert_module_solves(
            "(c - c*sin(e + f*x))^4/(a + a*sin(e + f*x))",
            &env,
            &[0.15, 0.35, 0.6],
        );
    }

    #[test]
    fn k1_sec_cubed_over_sec_squared_denom() {
        // Corpus shape: sec³(u)/(a + b·sec(u))² — Laurent after conversion.
        // Structural check only: see `assert_module_no_residue` for why the
        // emitted complex-branched form cannot be verified elementwise.
        assert_module_no_residue("sec(c + d*x)^3/(a + b*sec(c + d*x))^2");
    }

    #[test]
    fn k1_abc_numerator_poly() {
        // Corpus shape: (A + B·sec + C·sec²)/(a + b·sec)². Structural check
        // only: see `assert_module_no_residue`.
        assert_module_no_residue("(A + B*sec(c + d*x) + C*sec(c + d*x)^2)/(a + b*sec(c + d*x))^2");
    }

    // ------------------------- K2 -------------------------

    #[test]
    fn k2_tan2_sec2() {
        // ∫ tan²·sec² dx = tan³/3 (u-substitution level; numeric check).
        assert_module_solves("tan(x)^2*sec(x)^2", &[], &[0.3, 0.6, 1.0]);
    }

    #[test]
    fn k2_tan3_sec2_linear() {
        // Odd m, even n ≥ 2 → the n-even branch fires first.
        let env = [(Symbol::new("c"), 0.2), (Symbol::new("d"), 0.9)];
        assert_module_solves("tan(c + d*x)^3*sec(c + d*x)^2", &env, &[0.2, 0.4, 0.7]);
    }

    #[test]
    fn k2_tan3_sec3() {
        // Odd m, odd n ≥ 1 → the tan·sec peel branch.
        assert_module_solves("tan(x)^3*sec(x)^3", &[], &[0.3, 0.6, 1.0]);
    }

    #[test]
    fn k2_tan2_sec3_chain_reentry() {
        // Even m, odd n → rewritten to sec powers, chain re-entry.
        assert_module_solves("tan(x)^2*sec(x)^3", &[], &[0.3, 0.6, 1.0]);
    }

    #[test]
    fn k2_cot_csc_mirror() {
        assert_module_solves("cot(x)^2*csc(x)^2", &[], &[0.4, 0.8, 1.1]);
        assert_module_solves("cot(x)^3*csc(x)^3", &[], &[0.4, 0.8, 1.1]);
    }

    // ------------------------- declines -------------------------

    #[test]
    fn declines_out_of_scope() {
        // T1's shape (constant numerator over a binomial power).
        assert_module_declined("1/(2 + 3*cos(x))^2");
        // Fractional kernel power (elliptic family).
        assert_module_declined("1/(2*sec(x))^(3/2)");
        // Mixed sin/cos kernels.
        assert_module_declined("sin(x)^2*cos(x)/(1 + cos(x))");
        // Two different arguments.
        assert_module_declined("sec(x)/(1 + sec(2*x))");
        // Nonlinear argument.
        assert_module_declined("sec(x^2)/(1 + sec(x^2))");
    }

    #[test]
    fn chain_end_to_end() {
        // Through the public chain (wired by the main line).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "(3 - 3*sin(2*x))^2/(1 + sin(2*x))");
        let r = integrate(&ctx, expr, Symbol::new("x"));
        assert!(
            !r.to_string().contains("Integral"),
            "chain left residue: {r}"
        );
    }

    // ------------------------- deterministic budget ---------------------

    /// Corpus shapes that used to exhaust the per-case budget with this
    /// stage last entered (Rubi ids in the comments). They are hyperbolic,
    /// so the kernel matchers must decline them immediately instead of
    /// driving the peel/re-entry machinery.
    const HANG_SHAPES: &[&str] = &[
        "sech(x)/(a + b*sinh(x))",   // rubi-00262
        "sech(x)/(a + b*csch(x))",   // rubi-00918
        "tanh(x)^3/(a + b*csch(x))", // rubi-01100
    ];

    #[test]
    fn corpus_hang_shapes_decline() {
        for input in HANG_SHAPES {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse_norm(&ctx, input);
            // Direct module entry: must return (either a definite answer or
            // a decline) without driving the K1 peel.
            match integrate_trig_kernel(&ctx, expr, Symbol::new("x")) {
                None => {}
                Some(r) => assert!(
                    !r.to_string().contains("Integral"),
                    "{input} produced a residue: {r}"
                ),
            }
        }
    }

    /// The peel/expansion reductions that already solve must keep solving.
    #[test]
    fn budget_keeps_solving_normal_inputs() {
        assert_module_solves("tan(x)^2*sec(x)^2", &[], &[0.3, 0.6, 1.0]);
        assert_module_solves("tan(x)^3*sec(x)^3", &[], &[0.3, 0.6, 1.0]);
        // `sec^n/(a+b·sec)^m` emits a complex-branched form; structural check
        // only (see `assert_module_no_residue`).
        assert_module_no_residue("sec(c + d*x)^3/(a + b*sec(c + d*x))^2");
    }

    /// The work budget resets on every entry: repeating a hanging shape 256
    /// times must give identical results and stay fast.
    #[test]
    fn budget_does_not_leak_across_calls() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "sech(x)/(a + b*sinh(x))");
        let start = std::time::Instant::now();
        let first = integrate_trig_kernel(&ctx, expr, Symbol::new("x")).map(|a| a.to_string());
        for i in 0..256 {
            let r = integrate_trig_kernel(&ctx, expr, Symbol::new("x")).map(|a| a.to_string());
            assert_eq!(r, first, "call {i} diverged");
        }
        assert!(
            start.elapsed().as_secs() < 30,
            "256 calls took {:?}; the budget is not containing the shape",
            start.elapsed()
        );
    }

    /// A heavy peel enumeration (s + n at the cap) must still terminate and
    /// stay consistent across repeated calls.
    #[test]
    fn peel_enumeration_is_charged_and_stable() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse_norm(&ctx, "sec(c + d*x)^7/(a + b*sec(c + d*x))^7");
        let var = Symbol::new("x");
        let first = integrate_trig_kernel(&ctx, expr, var).map(|a| a.to_string());
        for i in 0..64 {
            let r = integrate_trig_kernel(&ctx, expr, var).map(|a| a.to_string());
            assert_eq!(r, first, "call {i} diverged");
        }
    }
}
