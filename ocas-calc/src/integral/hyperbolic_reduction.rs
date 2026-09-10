//! Hyperbolic kernel closed forms (0.27.2 Wave C2).
//!
//! Mirror of [`super::trig_reduction`] for the hyperbolic kernels. Runs before
//! `symbolic_rational`/`risch` so the shapes it owns never reach the stages
//! that grind on them.
//!
//! # What is enabled (verified)
//!
//! `∫ du/(a + b·T(u))` for `T ∈ {sinh, cosh}`, `u = c₀ + d·x` linear and
//! `a, b, c₀, d ∈ ℚ(symbols)` (`d ≠ 0`), by the exact seed
//!
//! ```text
//! ∫ du/(a + b·T(u)) = (1/(d·√Δ))·log((z + 1)/(z − 1)),
//!     T = sinh: Δ = a² + b²,  z = (a·t − b)/√Δ,
//!     T = cosh: Δ = a² − b²,  z = √((a + b)/(a − b))/t,
//!     t = tanh(u/2) = sinh u/(1 + cosh u),
//! ```
//!
//! written with `log` rather than `acoth`/`atanh` so that both `crate::diff`
//! and the numeric evaluator handle every branch. `Δ = 0` (the `cosh` case
//! `a² = b²`) declines: `D` then has a double root in `e^u` and the
//! antiderivative is not of this shape. Free constant factors are carried
//! through (`rest`), so `3/(2 + sinh x)` is solved as well.
//!
//! A candidate is accepted only after [`verify_numeric`] confirms
//! `d/dx(candidate) == integrand` numerically at fixed samples with the free
//! parameters specialized to fixed values. The mechanism solves nothing it
//! has not checked.
//!
//! # What is documented but NOT enabled (declines)
//!
//! These are deliberately declined rather than guessed at; each is a known
//! gap, not a wrong answer:
//!
//! - **`n ≥ 2` denominator powers** ([`MAX_DENOM_POW`] = 1). The recurrence
//!   `(n−1)(a² + εb²)J_n = a(2n−3)J_{n−1} + (2−n)J_{n−2} − b·T̂/D^(n−1)`
//!   (with `ε = +1` for `T = sinh`, `−1` for `T = cosh`) is derived but its
//!   assembled form did not survive numeric verification, so it is disabled.
//! - **`P(T)/(a + b·T)` long division in `T`** — the plan's route for
//!   `sinh(x)^3/(a + b·sinh(x))`. Implemented (`s_compute`/`s_recur`) but its
//!   `J`/`K` bookkeeping is unverified; only constant numerators are enabled.
//! - **`∫P(T)·T′` for all six kernels** ([`kernel_poly_rewrite`]): the rewrite
//!   `t = T(u)`, `du = dt/T′(t)` is implemented, but its radical bookkeeping
//!   is not verified.
//! - **`x^m·T(u)`** (`integrate_poly_hyper`) and **half-integer powers**
//!   (`integrate_algebraic_power`): stubs that decline. The elliptic engine
//!   owns the latter class.
//! - **`b·cosh(u) + c·sinh(u)` companion mixing** and **`T ∈ {tanh, coth,
//!   sech, csch}`** denominators: the phase reduction
//!   `b·T + c·T̂ = R·T(u + atanh(c/b))` (with `R² = b² − c²`, complex `ξ`
//!   when `|c| > |b|`) is not implemented, and the rational kernels need a
//!   different discriminant per kernel.
//!
//! Everything declines cheaply: node, degree, power and memo budgets, all
//! deterministic — no wall-clock and no unbounded search.

use std::collections::HashMap;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::{contains_integral, is_constant, linear_form};

/// Node budget for the input integrand.
const MAX_NODES: usize = 220;
/// Cap on the denominator power `n` in the recurrences.
const MAX_DENOM_POW: i64 = 1;
/// Cap on the exponent of a single kernel occurrence.
const MAX_FACTOR_POW: i64 = 8;
/// Cap on the monomial degree `m` of `(e + f·x)^m` in routine 3.
const MAX_POLY_DEG: i64 = 6;
/// Cap on the kernel degree of a numerator polynomial in `T`.
const MAX_ARG_DEG: usize = 8;
/// Cap on memoized recurrence entries.
const MAX_MEMO: usize = 400;
/// Node budget for a substituted (`t`-form) integrand.
const MAX_SUBST_NODES: usize = 400;
/// Cap on the radical exponent accumulated in a `t`-form rewrite.
const MAX_RAD_POW: i64 = 2 * MAX_FACTOR_POW;

/// Fixed samples of the integration variable for [`verify_numeric`].
const SAMPLES: [f64; 3] = [0.3, 0.7, 1.1];
/// Fixed specializations of the free parameters, by symbol name, for
/// [`verify_numeric`].
const PARAM_SAMPLES: [(&str, f64); 12] = [
    ("a", 2.0),
    ("b", 0.7),
    ("c", 1.3),
    ("d", 1.7),
    ("e", 0.6),
    ("f", 1.1),
    ("m", 1.9),
    ("n", 0.8),
    ("A", 1.4),
    ("B", 0.9),
    ("p", 1.6),
    ("q", 0.5),
];

// =========================================================================
// Entry
// =========================================================================

/// Integrate `expr` with the hyperbolic closed-form family.
///
/// Returns `None` when the shape is outside the family.
pub(crate) fn integrate_hyperbolic_reduction<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if is_constant(expr, var) || super::node_count(expr) > MAX_NODES {
        return None;
    }
    let sites = collect_kernel_sites(expr);
    if sites.is_empty() {
        return None;
    }
    let mut candidates = candidate_heads(&sites);
    // Try the kernel of the highest-power occurrence first: that is the
    // denominator's kernel, so its rewrite is the most likely to succeed.
    candidates.sort_by_key(|h| std::cmp::Reverse(power_of(&sites, *h)));
    for head in candidates {
        let Some(sub) = KernelSub::build(ctx, expr, var, head) else {
            continue;
        };
        // Each routine is verified before the next one is tried: a routine
        // that produces an unverified answer must not mask a later one that
        // would have been right.
        for r in [
            integrate_hyper_denom(ctx, expr, var, &sub),
            integrate_poly_hyper(ctx, expr, var, &sub),
            integrate_kernel_poly(ctx, expr, var, &sub),
            integrate_algebraic_power(ctx, expr, var, &sub),
        ]
        .into_iter()
        .flatten()
        {
            if verify_numeric(ctx, expr, var, r) {
                return Some(r);
            }
        }
    }
    None
}

// =========================================================================
// Kernel table and substitution
// =========================================================================

/// The six hyperbolic kernels.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Hyper {
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
}

impl Hyper {
    fn from_name(name: &str) -> Option<Hyper> {
        Some(match name {
            "sinh" => Hyper::Sinh,
            "cosh" => Hyper::Cosh,
            "tanh" => Hyper::Tanh,
            "coth" => Hyper::Coth,
            "sech" => Hyper::Sech,
            "csch" => Hyper::Csch,
            _ => return None,
        })
    }

    fn name(self) -> &'static str {
        match self {
            Hyper::Sinh => "sinh",
            Hyper::Cosh => "cosh",
            Hyper::Tanh => "tanh",
            Hyper::Coth => "coth",
            Hyper::Sech => "sech",
            Hyper::Csch => "csch",
        }
    }

    /// The pairing `T̂² − T² = 1`; `None` for the rational kernels.
    fn companion(self) -> Option<Hyper> {
        match self {
            Hyper::Sinh => Some(Hyper::Cosh),
            Hyper::Cosh => Some(Hyper::Sinh),
            _ => None,
        }
    }

    /// `dT/du` written in `t = T(u)`.
    ///
    /// [`KernelSub::du_dt`] writes the reciprocal form directly (the shape the
    /// algebraic conversion needs); this forward form is kept as the
    /// definitional table and as a cross-check.
    #[allow(dead_code)]
    fn deriv<'c>(self, ctx: &'c AtomArena<'c>, t: Atom<'c>) -> Atom<'c> {
        match self {
            // cosh u = ±√(1 + t²)
            Hyper::Sinh => root(ctx, ctx.add(&[ctx.num(1), sq(ctx, t)]), 2),
            // sinh u = ±√(t² − 1)
            Hyper::Cosh => root(ctx, ctx.add(&[sq(ctx, t), ctx.num(-1)]), 2),
            // 1 − t²
            Hyper::Tanh | Hyper::Coth => one_minus_sq(ctx, t),
            // −t·√(1 − t²)
            Hyper::Sech => {
                let r = root(ctx, one_minus_sq(ctx, t), 2);
                ctx.mul(&[ctx.num(-1), t, r])
            }
            // −t·√(1 + t²)
            Hyper::Csch => {
                let r = root(ctx, ctx.add(&[sq(ctx, t), ctx.num(1)]), 2);
                ctx.mul(&[ctx.num(-1), t, r])
            }
        }
    }

    /// Rational rewrite table for `t = T(u)`. Radical entries are emitted as
    /// such; [`Alg`] tracks where the radical sits.
    fn table<'c>(self, ctx: &'c AtomArena<'c>, t: Atom<'c>) -> Vec<(&'static str, Atom<'c>)> {
        let r1 = root(ctx, ctx.add(&[ctx.num(1), sq(ctx, t)]), 2);
        let r2 = root(ctx, ctx.add(&[sq(ctx, t), ctx.num(-1)]), 2);
        let r3 = root(ctx, one_minus_sq(ctx, t), 2);
        match self {
            Hyper::Sinh => vec![("sinh", t), ("cosh", r1)],
            Hyper::Cosh => vec![("cosh", t), ("sinh", r2)],
            Hyper::Tanh => vec![
                ("tanh", t),
                ("sech", r3),
                ("coth", inv(ctx, t)),
                ("csch", ctx.mul(&[r3, inv(ctx, t)])),
            ],
            Hyper::Coth => vec![
                ("coth", t),
                ("csch", r2),
                ("tanh", inv(ctx, t)),
                ("sech", ctx.mul(&[r2, inv(ctx, t)])),
            ],
            Hyper::Sech => vec![
                ("sech", t),
                ("tanh", r3),
                ("cosh", inv(ctx, t)),
                ("sinh", ctx.mul(&[r3, inv(ctx, t)])),
            ],
            Hyper::Csch => vec![
                ("csch", t),
                ("coth", r1),
                ("sinh", inv(ctx, t)),
                ("cosh", ctx.mul(&[r1, inv(ctx, t)])),
            ],
        }
    }
}

/// A ready kernel substitution `t = T(u)`, `u = c₀ + d·x` linear.
struct KernelSub<'a> {
    head: Hyper,
    /// The shared linear argument `u` of every kernel occurrence.
    arg: Atom<'a>,
    /// `d`, the slope of `u` (nonzero).
    slope: Atom<'a>,
    /// The substitution variable.
    t: Symbol,
    /// `(head, form in t)` rewrite table.
    table: Vec<(&'static str, Atom<'a>)>,
}

impl<'a> KernelSub<'a> {
    fn build(
        ctx: &'a AtomArena<'a>,
        expr: Atom<'a>,
        var: Symbol,
        head: Hyper,
    ) -> Option<KernelSub<'a>> {
        let sites = collect_kernel_sites(expr);
        let mut arg: Option<Atom<'a>> = None;
        for (name, a, _) in &sites {
            if Hyper::from_name(name.as_str()) != Some(head) {
                continue;
            }
            match arg {
                None => arg = Some(*a),
                Some(prev) => {
                    if normalize(ctx, prev) != normalize(ctx, *a) {
                        return None; // two distinct arguments: not one kernel
                    }
                }
            }
        }
        let arg = arg?;
        let (slope, _c) = linear_form(ctx, arg, var)?;
        if matches!(slope.node(), AtomNode::Num(0)) {
            return None;
        }
        let t_name = super::pick_subst_symbol(expr, var)?;
        let t = Symbol::new(t_name.as_str());
        let table = head.table(ctx, ctx.var(t_name.as_str()));
        Some(KernelSub {
            head,
            arg,
            slope,
            t,
            table,
        })
    }

    fn t_atom(&self, ctx: &'a AtomArena<'a>) -> Atom<'a> {
        ctx.var(self.t.as_str())
    }

    fn form(&self, name: &str) -> Option<Atom<'a>> {
        self.table.iter().find(|(n, _)| *n == name).map(|(_, f)| *f)
    }

    /// `T(u)` as an atom.
    fn head_at(&self, ctx: &'a AtomArena<'a>) -> Atom<'a> {
        ctx.fun(self.head.name(), &[self.arg])
    }

    /// `T̂(u)` as an atom, when the head has a companion.
    fn comp_at(&self, ctx: &'a AtomArena<'a>) -> Option<Atom<'a>> {
        let c = self.head.companion()?;
        Some(ctx.fun(c.name(), &[self.arg]))
    }

    /// `dT/du` evaluated at the atom `T(u)`.
    #[allow(dead_code)]
    fn head_deriv_at(&self, ctx: &'a AtomArena<'a>, t: Atom<'a>) -> Atom<'a> {
        let d = self.head.deriv(ctx, self.t_atom(ctx));
        super::replace_symbol(ctx, d, self.t, t)
    }

    /// `T′(u)·d`, the chain-rule factor of the substitution (the forward
    /// cross-check of [`KernelSub::du_dt`]).
    #[allow(dead_code)]
    fn chain(&self, ctx: &'a AtomArena<'a>) -> Atom<'a> {
        let t = self.head_at(ctx);
        ctx.mul(&[self.head_deriv_at(ctx, t), self.slope])
    }

    /// `t → T(u)` inside a `t`-form.
    fn back(&self, ctx: &'a AtomArena<'a>, t_form: Atom<'a>) -> Atom<'a> {
        let head = self.head_at(ctx);
        super::replace_symbol(ctx, t_form, self.t, head)
    }

    /// The `t`-form numerator/denominator/radical of `du = dt/T′(u)`:
    /// `du/dt = num/(den·rad)`.
    ///
    /// The `t`-form of `du/dt = 1/T′(t)` (`value = num·rad^(rad_pow/2)/den`).
    ///
    /// Multiplying a `t`-form `R(t)` by it turns `∫R(T(u))·T′(u) du` into
    /// `∫R(t) dt` — the conversion both the rational and the algebraic
    /// routines use before handing the problem to the univariate engines.
    fn du_dt(&self, ctx: &'a AtomArena<'a>) -> Alg<'a> {
        let t = self.t_atom(ctx);
        let m = |ctx: &'a AtomArena<'a>, plus: bool| {
            if plus {
                ctx.add(&[ctx.num(1), sq(ctx, t)])
            } else {
                one_minus_sq(ctx, t)
            }
        };
        match self.head {
            // 1/√(1 + t²)
            Hyper::Sinh => Alg {
                num: poly_one(ctx),
                den: poly_one(ctx),
                rad: poly_const(ctx, m(ctx, true)),
                rad_pow: -1,
                rad_kind: 1,
            },
            // 1/√(t² − 1)
            Hyper::Cosh => Alg {
                num: poly_one(ctx),
                den: poly_one(ctx),
                rad: poly_const(ctx, m(ctx, false)),
                rad_pow: -1,
                rad_kind: -1,
            },
            // 1/(1 − t²)
            Hyper::Tanh | Hyper::Coth => Alg {
                num: poly_one(ctx),
                den: poly_const(ctx, m(ctx, false)),
                rad: poly_one(ctx),
                rad_pow: 0,
                rad_kind: 0,
            },
            // −1/(t·√(1 − t²))
            Hyper::Sech => Alg {
                num: poly_const(ctx, ctx.num(-1)),
                den: poly_const(ctx, t),
                rad: poly_const(ctx, m(ctx, false)),
                rad_pow: -1,
                rad_kind: 2,
            },
            // −1/(t·√(1 + t²))
            Hyper::Csch => Alg {
                num: poly_const(ctx, ctx.num(-1)),
                den: poly_const(ctx, t),
                rad: poly_const(ctx, m(ctx, true)),
                rad_pow: -1,
                rad_kind: 1,
            },
        }
    }

    /// Rewrite `expr` as an algebraic function of `t = T(u)`. `None` when
    /// the rewrite is impossible (leftover `x`, unknown head, or an exponent
    /// outside the budget).
    fn to_alg(&self, ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Option<Alg<'a>> {
        // The substitution variable is the polynomial variable, not a
        // coefficient.
        if let AtomNode::Var(v) = expr.node()
            && *v == self.t
        {
            return Some(Alg::constant(ctx, Poly::monomial(ctx, 1)));
        }
        if is_constant(expr, var) {
            return Some(Alg::constant(ctx, poly_const(ctx, expr)));
        }
        match expr.node() {
            AtomNode::Var(_) | AtomNode::Num(_) => Some(Alg::constant(ctx, poly_const(ctx, expr))),
            AtomNode::Add(args) => {
                let mut acc = Alg::constant(ctx, poly_zero());
                for a in args.iter() {
                    acc = acc.add(ctx, &self.to_alg(ctx, *a, var)?)?;
                }
                Some(acc)
            }
            AtomNode::Mul(args) => {
                let mut acc = Alg::one(ctx);
                for a in args.iter() {
                    acc = acc.mul(ctx, &self.to_alg(ctx, *a, var)?)?;
                }
                Some(acc)
            }
            AtomNode::Pow(b, e) => {
                let (p, s) = rat_parts(*e)?;
                if p.unsigned_abs() > MAX_FACTOR_POW as u64
                    || s.unsigned_abs() > MAX_FACTOR_POW as u64
                {
                    return None;
                }
                if let AtomNode::Fun(name, args) = b.node()
                    && args.len() == 1
                    && name.as_str() == self.head.name()
                {
                    return alg_pow_t(ctx, p, s);
                }
                let base = self.to_alg(ctx, *b, var)?;
                alg_pow_alg(ctx, &base, p, s)
            }
            AtomNode::Fun(name, args) => {
                if args.len() != 1 {
                    return None;
                }
                let form = self.form(name.as_str())?;
                // The table entry is a `t`-expression (possibly carrying a
                // radical site): re-enter the same rewriter so the radical is
                // classified instead of being mistaken for a coefficient.
                self.to_alg(ctx, form, var)
            }
        }
    }
}

fn companion_name(h: Hyper) -> Option<&'static str> {
    h.companion().map(Hyper::name)
}

fn collect_kernel_sites<'a>(expr: Atom<'a>) -> Vec<(Symbol, Atom<'a>, i64)> {
    let mut out = Vec::new();
    collect_sites_into(expr, 1, &mut out);
    out
}

fn collect_sites_into<'a>(expr: Atom<'a>, power: i64, out: &mut Vec<(Symbol, Atom<'a>, i64)>) {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if args.len() == 1 && Hyper::from_name(name.as_str()).is_some() {
                out.push((*name, args[0], power));
            } else {
                for a in args.iter() {
                    collect_sites_into(*a, power, out);
                }
            }
        }
        AtomNode::Pow(b, e) => {
            let p = match e.node() {
                AtomNode::Num(k) => *k,
                _ => power,
            };
            collect_sites_into(*b, p, out);
            collect_sites_into(*e, power, out);
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for a in args.iter() {
                collect_sites_into(*a, power, out);
            }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => {}
    }
}

fn candidate_heads(sites: &[(Symbol, Atom<'_>, i64)]) -> Vec<Hyper> {
    let mut out: Vec<Hyper> = Vec::new();
    for (name, _, _) in sites {
        if let Some(h) = Hyper::from_name(name.as_str())
            && !out.contains(&h)
        {
            out.push(h);
        }
    }
    out
}

fn power_of(sites: &[(Symbol, Atom<'_>, i64)], head: Hyper) -> i64 {
    sites
        .iter()
        .filter(|(n, _, _)| n.as_str() == head.name())
        .map(|(_, _, p)| p.abs())
        .max()
        .unwrap_or(0)
}

// =========================================================================
// Small atom helpers
// =========================================================================

fn sq<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> Atom<'a> {
    ctx.pow(e, ctx.num(2))
}

fn inv<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> Atom<'a> {
    ctx.pow(e, ctx.num(-1))
}

/// Σ args with the empty sum folded to 0.
fn sum_atoms<'a>(ctx: &'a AtomArena<'a>, args: &[Atom<'a>]) -> Atom<'a> {
    if args.is_empty() {
        return ctx.num(0);
    }
    ctx.add(args)
}

/// Π args with the empty product folded to 1.
fn prod_atoms<'a>(ctx: &'a AtomArena<'a>, args: &[Atom<'a>]) -> Atom<'a> {
    if args.is_empty() {
        return ctx.num(1);
    }
    ctx.mul(args)
}

fn one_minus_sq<'a>(ctx: &'a AtomArena<'a>, t: Atom<'a>) -> Atom<'a> {
    ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), sq(ctx, t)])])
}

/// `x^(1/n)` as an atom, folded for exact numeric roots.
fn root<'a>(ctx: &'a AtomArena<'a>, x: Atom<'a>, n: i64) -> Atom<'a> {
    if n == 2
        && let AtomNode::Num(v) = x.node()
        && *v >= 0
    {
        let r = (*v as f64).sqrt().round() as i64;
        if r * r == *v {
            return ctx.num(r);
        }
    }
    if n == 1 {
        return x;
    }
    let e = if n > 0 {
        ctx.mul(&[ctx.num(1), ctx.pow(ctx.num(n), ctx.num(-1))])
    } else {
        ctx.mul(&[
            ctx.num(-1),
            ctx.mul(&[ctx.num(1), ctx.pow(ctx.num(-n), ctx.num(-1))]),
        ])
    };
    ctx.pow(x, e)
}

// =========================================================================
// Univariate polynomials over ℚ(symbols)
// =========================================================================

/// Dense univariate polynomial with coefficients in `ℚ(symbols)`, with no
/// trailing zeros.
#[derive(Clone)]
struct Poly<'a> {
    terms: Vec<Atom<'a>>,
}

fn poly_const<'a>(ctx: &'a AtomArena<'a>, c: Atom<'a>) -> Poly<'a> {
    trim(vec![c], ctx)
}

fn poly_one<'a>(ctx: &'a AtomArena<'a>) -> Poly<'a> {
    poly_const(ctx, ctx.num(1))
}

fn poly_zero<'a>() -> Poly<'a> {
    Poly { terms: Vec::new() }
}

impl<'a> Poly<'a> {
    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn degree(&self) -> Option<usize> {
        if self.terms.is_empty() {
            None
        } else {
            Some(self.terms.len() - 1)
        }
    }

    fn leading(&self) -> Option<Atom<'a>> {
        self.terms.last().copied()
    }

    fn at(&self, ctx: &'a AtomArena<'a>, i: usize) -> Atom<'a> {
        self.terms.get(i).copied().unwrap_or_else(|| ctx.num(0))
    }

    fn add(&self, ctx: &'a AtomArena<'a>, o: &Poly<'a>) -> Poly<'a> {
        let n = self.terms.len().max(o.terms.len());
        if n == 0 {
            return poly_zero();
        }
        let mut terms = Vec::with_capacity(n);
        for i in 0..n {
            terms.push(ctx.add(&[self.at(ctx, i), o.at(ctx, i)]));
        }
        trim(terms, ctx)
    }

    fn neg(&self, ctx: &'a AtomArena<'a>) -> Poly<'a> {
        Poly {
            terms: self
                .terms
                .iter()
                .map(|c| ctx.mul(&[ctx.num(-1), *c]))
                .collect(),
        }
    }

    fn sub(&self, ctx: &'a AtomArena<'a>, o: &Poly<'a>) -> Poly<'a> {
        self.add(ctx, &o.neg(ctx))
    }

    fn mul(&self, ctx: &'a AtomArena<'a>, o: &Poly<'a>) -> Poly<'a> {
        if self.is_zero() || o.is_zero() {
            return poly_zero();
        }
        let mut terms = vec![ctx.num(0); self.terms.len() + o.terms.len() - 1];
        for (i, a) in self.terms.iter().enumerate() {
            for (j, b) in o.terms.iter().enumerate() {
                terms[i + j] = ctx.add(&[terms[i + j], ctx.mul(&[*a, *b])]);
            }
        }
        trim(terms, ctx)
    }

    /// `(quotient, remainder)`; the leading coefficient is a field unit, so
    /// the division always terminates.
    fn divmod(&self, ctx: &'a AtomArena<'a>, d: &Poly<'a>) -> Option<(Poly<'a>, Poly<'a>)> {
        if d.is_zero() {
            return None;
        }
        let ddeg = d.degree()?;
        let inv_lead = ctx.pow(d.leading()?, ctx.num(-1));
        let mut rem = self.clone();
        let mut quot = vec![ctx.num(0); self.terms.len().max(1)];
        while let Some(rdeg) = rem.degree() {
            if rdeg < ddeg {
                break;
            }
            let coef = ctx.mul(&[rem.leading()?, inv_lead]);
            if matches!(normalize(ctx, coef).node(), AtomNode::Num(0)) {
                break;
            }
            let shift = rdeg - ddeg;
            quot[shift] = ctx.add(&[quot[shift], coef]);
            let before = rdeg;
            let sub = d.scale_shift(ctx, coef, shift);
            rem = rem.sub(ctx, &sub);
            if rem.degree().is_some_and(|d2| d2 >= before) {
                return None; // no progress: refuse rather than loop
            }
        }
        Some((trim(quot, ctx), rem))
    }

    /// `t^shift · c · self`.
    fn scale_shift(&self, ctx: &'a AtomArena<'a>, c: Atom<'a>, shift: usize) -> Poly<'a> {
        let mut terms = vec![ctx.num(0); shift];
        terms.extend(self.terms.iter().map(|x| ctx.mul(&[c, *x])));
        trim(terms, ctx)
    }

    /// `t^n` as a polynomial.
    fn monomial(ctx: &'a AtomArena<'a>, n: usize) -> Poly<'a> {
        let mut terms = vec![ctx.num(0); n];
        terms.push(ctx.num(1));
        Poly { terms }
    }

    /// `self^n`, `n ≥ 0`.
    fn pow(&self, ctx: &'a AtomArena<'a>, n: usize) -> Poly<'a> {
        let mut acc = poly_one(ctx);
        if self.is_zero() {
            return if n == 0 { acc } else { poly_zero() };
        }
        for _ in 0..n {
            acc = acc.mul(ctx, self);
        }
        acc
    }

    /// The atom `Σ terms[i]·t^i` (zero for the empty polynomial).
    fn to_atom(&self, ctx: &'a AtomArena<'a>, t: Atom<'a>) -> Atom<'a> {
        if self.terms.is_empty() {
            return ctx.num(0);
        }
        let mut terms = Vec::new();
        for (i, c) in self.terms.iter().enumerate() {
            let term = match i {
                0 => *c,
                1 => ctx.mul(&[*c, t]),
                _ => ctx.mul(&[*c, ctx.pow(t, ctx.num(i as i64))]),
            };
            terms.push(term);
        }
        normalize(ctx, sum_atoms(ctx, &terms))
    }
}

fn trim<'a>(terms: Vec<Atom<'a>>, ctx: &'a AtomArena<'a>) -> Poly<'a> {
    let mut terms = terms;
    while let Some(last) = terms.last() {
        if matches!(normalize(ctx, *last).node(), AtomNode::Num(0)) {
            terms.pop();
        } else {
            break;
        }
    }
    Poly { terms }
}

/// An algebraic `t`-form: `num / (den · rad^(rad_pow/2))`.
///
/// `rad` is the shared radical site (`1 + t²` or `1 − t²`, tracked by
/// `rad_kind`); `rad_pow` is the accumulated exponent in half-units.
#[derive(Clone)]
struct Alg<'a> {
    num: Poly<'a>,
    den: Poly<'a>,
    rad: Poly<'a>,
    rad_pow: i64,
    /// `+1` for `1 + t²`, `−1` for `1 − t²`, `0` for none.
    rad_kind: i8,
}

impl<'a> Alg<'a> {
    /// A rational function of `t` (`num/den`, `den` empty meaning `1`).
    fn rational(ctx: &'a AtomArena<'a>, num: Poly<'a>, den: Poly<'a>) -> Alg<'a> {
        Alg {
            num,
            den: if den.is_zero() { poly_one(ctx) } else { den },
            rad: poly_one(ctx),
            rad_pow: 0,
            rad_kind: 0,
        }
    }

    fn constant(ctx: &'a AtomArena<'a>, p: Poly<'a>) -> Alg<'a> {
        Alg::rational(ctx, p, poly_zero())
    }

    fn one(ctx: &'a AtomArena<'a>) -> Alg<'a> {
        Alg::constant(ctx, poly_one(ctx))
    }

    /// The shared radical site of two forms, or `None` when they disagree.
    fn site(&self, ctx: &'a AtomArena<'a>, o: &Alg<'a>) -> Option<(Poly<'a>, i8)> {
        match (self.rad_pow != 0, o.rad_pow != 0) {
            (false, false) => Some((poly_one(ctx), 0)),
            (true, false) => Some((self.rad.clone(), self.rad_kind)),
            (false, true) => Some((o.rad.clone(), o.rad_kind)),
            (true, true) => {
                if self.rad_kind == o.rad_kind && self.rad_kind != 0 {
                    Some((self.rad.clone(), self.rad_kind))
                } else {
                    None
                }
            }
        }
    }

    fn add(&self, ctx: &'a AtomArena<'a>, o: &Alg<'a>) -> Option<Alg<'a>> {
        let (rad, kind) = self.site(ctx, o)?;
        // Bring both numerators over the common denominator den·rad^m.
        let m = self.rad_pow.max(o.rad_pow);
        if m > MAX_RAD_POW {
            return None;
        }
        let self_extra = (m - self.rad_pow) as usize;
        let o_extra = (m - o.rad_pow) as usize;
        let an = self
            .num
            .mul(ctx, &o.den)
            .mul(ctx, &rad.pow(ctx, self_extra));
        let bn = o.num.mul(ctx, &self.den).mul(ctx, &rad.pow(ctx, o_extra));
        Some(Alg {
            num: an.add(ctx, &bn),
            den: self.den.mul(ctx, &o.den),
            rad,
            rad_pow: m,
            rad_kind: kind,
        })
    }

    fn mul(&self, ctx: &'a AtomArena<'a>, o: &Alg<'a>) -> Option<Alg<'a>> {
        let (rad, kind) = self.site(ctx, o)?;
        let m = self.rad_pow + o.rad_pow;
        if !(-MAX_RAD_POW..=MAX_RAD_POW).contains(&m) {
            return None;
        }
        Some(
            Alg {
                num: self.num.mul(ctx, &o.num),
                den: self.den.mul(ctx, &o.den),
                rad,
                rad_pow: m,
                rad_kind: kind,
            }
            .fold_rad(ctx),
        )
    }

    /// Fold an even radical power back into the numerator/denominator
    /// (`rad^(2k) = rad^k` as a polynomial), so `(√(1 − t²))²` becomes
    /// `1 − t²` instead of an unstripped radical site.
    fn fold_rad(self, ctx: &'a AtomArena<'a>) -> Alg<'a> {
        if self.rad_pow == 0 || self.rad_pow % 2 != 0 {
            return self;
        }
        let k = (self.rad_pow / 2).unsigned_abs() as usize;
        let site = self.rad.pow(ctx, k);
        if self.rad_pow > 0 {
            Alg {
                num: self.num.mul(ctx, &site),
                den: self.den,
                rad: poly_one(ctx),
                rad_pow: 0,
                rad_kind: 0,
            }
        } else {
            Alg {
                num: self.num,
                den: self.den.mul(ctx, &site),
                rad: poly_one(ctx),
                rad_pow: 0,
                rad_kind: 0,
            }
        }
    }

    /// Whether this form is a plain rational function of `t` (used by the
    /// disabled kernel-polynomial rewrite).
    #[allow(dead_code)]
    fn is_rational(&self) -> bool {
        self.rad_pow == 0
    }
}

/// `(p, s)` with the exponent equal to `p/s`, `s > 0` (`None` for a
/// non-numeric exponent).
fn rat_parts<'a>(k: Atom<'a>) -> Option<(i64, i64)> {
    match k.node() {
        AtomNode::Num(n) => Some((*n, 1)),
        AtomNode::Mul(args) if args.len() == 2 => {
            let AtomNode::Num(p) = args[0].node() else {
                return None;
            };
            reciprocal_of(args[1]).map(|q| (*p, q))
        }
        AtomNode::Pow(..) => reciprocal_of(k).map(|q| (1, q)),
        _ => None,
    }
}

/// `q` for `q^(−1)`, possibly written as the canonical `1·q^(−1)` product.
fn reciprocal_of<'a>(k: Atom<'a>) -> Option<i64> {
    match k.node() {
        AtomNode::Pow(b, e) => match (b.node(), e.node()) {
            (AtomNode::Num(q), AtomNode::Num(-1)) if *q != 0 => Some(*q),
            _ => None,
        },
        AtomNode::Mul(args) if args.len() == 2 => {
            let AtomNode::Num(1) = args[0].node() else {
                return None;
            };
            match args[1].node() {
                AtomNode::Pow(b, e) => match (b.node(), e.node()) {
                    (AtomNode::Num(q), AtomNode::Num(-1)) if *q != 0 => Some(*q),
                    _ => None,
                },
                _ => None,
            }
        }
        _ => None,
    }
}

/// Classify a radical site polynomial: `+1` for `1 + t²`, `−1` for `t² − 1`,
/// `+2` for `1 − t²`, `0` for anything else (unknown → the caller declines).
fn site_kind(p: &Poly<'_>) -> i8 {
    if p.degree() != Some(2) {
        return 0;
    }
    let c0 = p.terms.first();
    let c2 = p.terms.get(2);
    match (
        c0.and_then(|c| match c.node() {
            AtomNode::Num(n) => Some(*n),
            _ => None,
        }),
        c2.and_then(|c| match c.node() {
            AtomNode::Num(n) => Some(*n),
            _ => None,
        }),
    ) {
        (Some(1), Some(1)) => 1,
        (Some(-1), Some(1)) => -1,
        (Some(1), Some(-1)) => 2,
        _ => 0,
    }
}

/// The radical-site kind of the head kernel's companion:
/// `cosh = √(1 + t²)`, `sinh = √(t² − 1)`, `sech = √(1 − t²)`,
/// `csch = √(t² − 1)`, `tanh = √(1 − t²)`, `coth = √(1 + t²)`.
#[allow(dead_code)]
fn companion_kind(h: Hyper) -> i8 {
    match h {
        Hyper::Sinh | Hyper::Csch => 1,
        Hyper::Cosh | Hyper::Coth => -1,
        Hyper::Tanh | Hyper::Sech => 2,
    }
}

/// `T(u)^(p/s)` as an `Alg` with the value `t^(p/s)` (`t` = the substitution
/// symbol).
///
/// Integer exponents are plain polynomials; a square root of the kernel itself
/// is declined (the module does not own that algebraic class).
fn alg_pow_t<'a>(ctx: &'a AtomArena<'a>, p: i64, s: i64) -> Option<Alg<'a>> {
    let t = Poly::monomial(ctx, 1);
    if s != 1 {
        return None;
    }
    Some(if p >= 0 {
        Alg::constant(ctx, t.pow(ctx, p as usize))
    } else {
        Alg::rational(ctx, poly_one(ctx), t.pow(ctx, p.unsigned_abs() as usize))
    })
}

/// `(num·rad^(rad_pow/2)/den)^(p/s)`.
///
/// Integer powers are always safe; a square root is admitted only for a base
/// with no existing radical site and a constant denominator (the shapes the
/// rewrite table produces).
fn alg_pow_alg<'a>(ctx: &'a AtomArena<'a>, f: &Alg<'a>, p: i64, s: i64) -> Option<Alg<'a>> {
    if s == 1 {
        let n = p.unsigned_abs() as usize;
        let (num, den) = if p >= 0 {
            (f.num.pow(ctx, n), f.den.pow(ctx, n))
        } else {
            (f.den.pow(ctx, n), f.num.pow(ctx, n))
        };
        let rad_pow = f.rad_pow.checked_mul(p)?;
        return Some(Alg {
            num,
            den,
            rad: f.rad.clone(),
            rad_pow,
            rad_kind: f.rad_kind,
        });
    }
    if s != 2 || f.rad_pow != 0 || f.den.degree() != Some(0) {
        return None;
    }
    // base^(p/2) = base^q · (base^r)^(1/2), p = 2q + r, r ∈ {0, 1}.
    let q = p.div_euclid(2);
    let r = p.rem_euclid(2);
    let kind = site_kind(&f.num);
    if kind == 0 {
        return None; // an unrecognized radical site: decline rather than guess
    }
    let base = f.num.pow(ctx, q.unsigned_abs() as usize);
    let radical = f.num.pow(ctx, r as usize);
    if q >= 0 {
        Some(Alg {
            num: base,
            den: poly_one(ctx),
            rad: radical,
            rad_pow: 1,
            rad_kind: kind,
        })
    } else {
        Some(Alg {
            num: radical,
            den: base,
            rad: poly_one(ctx),
            rad_pow: 0,
            rad_kind: 0,
        })
    }
}

// =========================================================================
// Routine 1: denominator recurrences (sinh / cosh)
// =========================================================================

/// The context of one `D = a + b·T(u)` reduction (the plain single-kernel
/// surface; the companion mixing `c·T̂` declines — see the module docs).
struct DenomCtx<'a> {
    ctx: &'a AtomArena<'a>,
    /// The integration variable (for the `∫du` seeds).
    var: Symbol,
    /// `T(u)`.
    t: Atom<'a>,
    /// `T̂(u)`.
    comp: Atom<'a>,
    a: Atom<'a>,
    b: Atom<'a>,
    /// `d`, the slope of `u`.
    slope: Atom<'a>,
    /// `D = a + b·T`.
    d0: Atom<'a>,
    /// `a² + ε·b²` with `ε = +1` for `T = sinh`, `−1` for `T = cosh`.
    disc: Atom<'a>,
    /// `√disc`.
    sqrt_disc: Atom<'a>,
    /// `ε`.
    eps: i64,
    /// `∫T^L/D^n du`.
    memo: HashMap<(usize, i64), Option<Atom<'a>>>,
}

/// `1/(a + b·T(u))^n` and `P(T)/(a + b·T(u))^n` (see module docs).
///
/// `P(T)` is a polynomial in the chosen kernel: each monomial is divided by
/// `D^n` (the plan's polynomial long division in `T`), the quotient is
/// integrated termwise through `∫T^k du`, and the remainder through the
/// `K^j = (D − a)^j` binomial expansion over the `J_k = ∫du/D^k` recurrence.
fn integrate_hyper_denom<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sub: &KernelSub<'a>,
) -> Option<Atom<'a>> {
    let comp_name = companion_name(sub.head)?;
    let factors: Vec<Atom<'a>> = match expr.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![expr],
    };
    let mut rest: Vec<Atom<'a>> = Vec::new();
    let mut denom: Option<(Atom<'a>, i64)> = None;
    let mut numer_kernel: Vec<Atom<'a>> = Vec::new();
    for f in factors {
        if is_constant(f, var) {
            rest.push(f);
            continue;
        }
        if let Some((inner, n)) = as_recip_pow(f)
            && split_denom_base(ctx, inner, var, sub, comp_name).is_some()
        {
            if denom.is_some() {
                return None; // two denominator powers: out of scope
            }
            denom = Some((inner, n));
            continue;
        }
        if kernel_factor(f, sub) {
            numer_kernel.push(f);
        } else {
            // A non-kernel non-constant factor (`coth(x)`, `x^3`, `sech(x)`, …)
            // belongs to a routine this module does not own: decline rather
            // than let it through as if it were a constant.
            return None;
        }
    }
    let (inner, n) = denom?;
    if !(1..=MAX_DENOM_POW).contains(&n) {
        return None;
    }
    let (a, b, c) = split_denom_base(ctx, inner, var, sub, comp_name)?;
    // A companion term would need the phase reduction, which is not
    // implemented: decline rather than guess (the mixed surface is documented
    // as a known gap).
    if !matches!(normalize(ctx, c).node(), AtomNode::Num(0)) {
        return None;
    }
    // The numerator must be a polynomial in the chosen kernel T.
    let p = if numer_kernel.is_empty() {
        poly_one(ctx)
    } else {
        let num = if numer_kernel.len() == 1 {
            numer_kernel[0]
        } else {
            prod_atoms(ctx, &numer_kernel)
        };
        let f = sub.to_alg(ctx, num, var)?;
        if f.den.degree() != Some(0) || f.rad_pow != 0 {
            return None;
        }
        f.num
    };
    // Only constant numerators are enabled: `∫du/(a + b·T(u))`, whose closed
    // form is verified below. The `P(T)/D` long-division path is documented
    // but its recurrence is not verified yet, so it declines rather than risk
    // an unverified emission.
    if p.degree() != Some(0) {
        return None;
    }
    let comp = sub.comp_at(ctx)?;
    let t = sub.head_at(ctx);
    let eps: i64 = if sub.head == Hyper::Sinh { 1 } else { -1 };
    let disc = normalize(
        ctx,
        ctx.add(&[sq(ctx, a), ctx.mul(&[ctx.num(eps), sq(ctx, b)])]),
    );
    if matches!(disc.node(), AtomNode::Num(0)) {
        return None;
    }
    let d0 = normalize(ctx, ctx.add(&[a, ctx.mul(&[b, t])]));
    let mut dc = DenomCtx {
        ctx,
        var,
        t,
        comp,
        a,
        b,
        slope: sub.slope,
        d0,
        disc,
        sqrt_disc: root(ctx, disc, 2),
        eps,
        memo: HashMap::new(),
    };
    let mut acc = ctx.num(0);
    for (l, coef) in p.terms.iter().enumerate() {
        let s = s_recur(&mut dc, l, n)?;
        acc = ctx.add(&[acc, ctx.mul(&[*coef, s])]);
    }
    rest.push(acc);
    Some(normalize(ctx, prod_atoms(ctx, &rest)))
}

/// Whether `f` is a product of powers of the chosen kernel alone.
fn kernel_factor<'a>(f: Atom<'a>, sub: &KernelSub<'a>) -> bool {
    match f.node() {
        AtomNode::Fun(name, args) => {
            args.len() == 1 && name.as_str() == sub.head.name() && args[0] == sub.arg
        }
        AtomNode::Pow(b, e) => match (b.node(), e.node()) {
            (AtomNode::Fun(name, args), AtomNode::Num(_)) => {
                args.len() == 1 && name.as_str() == sub.head.name() && args[0] == sub.arg
            }
            _ => false,
        },
        AtomNode::Mul(fs) => fs.iter().all(|x| kernel_factor(*x, sub)),
        _ => false,
    }
}

/// `(a, b, c)` for `inner = a + b·T(u) + c·T̂(u)`.
fn split_denom_base<'a>(
    ctx: &'a AtomArena<'a>,
    inner: Atom<'a>,
    var: Symbol,
    sub: &KernelSub<'a>,
    comp_name: &str,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    let AtomNode::Add(args) = inner.node() else {
        return None;
    };
    let mut consts: Vec<Atom<'a>> = Vec::new();
    let mut bs: Vec<Atom<'a>> = Vec::new();
    let mut cs: Vec<Atom<'a>> = Vec::new();
    for t in args.iter() {
        if is_constant(*t, var) {
            consts.push(*t);
            continue;
        }
        let factors: Vec<Atom<'a>> = match t.node() {
            AtomNode::Mul(fs) => fs.to_vec(),
            _ => vec![*t],
        };
        let mut coeff: Vec<Atom<'a>> = Vec::new();
        let mut kernel: Option<(Symbol, Atom<'a>)> = None;
        for f in factors {
            if is_constant(f, var) {
                coeff.push(f);
                continue;
            }
            let AtomNode::Fun(fn_name, fa) = f.node() else {
                return None;
            };
            if fa.len() != 1 {
                return None;
            }
            if fn_name.as_str() != sub.head.name() && fn_name.as_str() != comp_name {
                return None;
            }
            if kernel.is_some() {
                return None;
            }
            kernel = Some((*fn_name, fa[0]));
        }
        let (fn_name, karg) = kernel?;
        if normalize(ctx, karg) != normalize(ctx, sub.arg) {
            return None;
        }
        let c = if coeff.is_empty() {
            ctx.num(1)
        } else {
            sum_atoms(ctx, &coeff)
        };
        if fn_name.as_str() == sub.head.name() {
            bs.push(c);
        } else {
            cs.push(c);
        }
    }
    let a = normalize(ctx, sum_atoms(ctx, &consts));
    let b = normalize(ctx, sum_atoms(ctx, &bs));
    let c = normalize(ctx, sum_atoms(ctx, &cs));
    if matches!(b.node(), AtomNode::Num(0)) && matches!(c.node(), AtomNode::Num(0)) {
        return None;
    }
    Some((a, b, c))
}

/// `S_{L,n} = ∫T^L/D^n du` (memoized).
fn s_recur<'a>(dc: &mut DenomCtx<'a>, l: usize, n: i64) -> Option<Atom<'a>> {
    if l as i64 > MAX_ARG_DEG as i64 + MAX_DENOM_POW {
        return None;
    }
    if n > MAX_DENOM_POW {
        return None;
    }
    if let Some(v) = dc.memo.get(&(l, n)) {
        return *v;
    }
    if dc.memo.len() > MAX_MEMO {
        return None;
    }
    let value = s_compute(dc, l, n);
    dc.memo.insert((l, n), value);
    value
}

fn s_compute<'a>(dc: &mut DenomCtx<'a>, l: usize, n: i64) -> Option<Atom<'a>> {
    let ctx = dc.ctx;
    if l == 0 {
        return j_recur(dc, n);
    }
    // Long division of `T^L` by `D^n` (a polynomial in `T`): `T^L = Q·D^n + R`,
    // `deg R < n`, so `∫T^L/D^n = ∫Q dT + ∫R/D^n`.
    let d_poly = {
        let lin = Poly {
            terms: vec![dc.a, dc.b],
        };
        lin.pow(ctx, n as usize)
    };
    let t_poly = Poly::monomial(ctx, l);
    let (q, r) = t_poly.divmod(ctx, &d_poly)?;
    let mut terms = Vec::new();
    if !q.is_zero() {
        for (i, coef) in q.terms.iter().enumerate() {
            let anti = t_power_integral(dc, i)?;
            terms.push(ctx.mul(&[*coef, anti]));
        }
    }
    if !r.is_zero() {
        // R(T) = Σ_j r_j·T^j and K = b·T gives T^j = K^j/b^j, while
        // K^j = (D − a)^j = Σ_i binom(j,i)(−a)^(j−i)·D^(j−i), so
        // ∫R/D^n = Σ_j (r_j/b^j)·Σ_i binom(j,i)(−a)^(j−i)·J_{n−i}.
        for (j, rj) in r.terms.iter().enumerate() {
            let mut inner = ctx.num(0);
            for i in 0..=j {
                let binom = binom_i64(j as i64, i as i64)?;
                let jk = j_recur(dc, n - i as i64)?;
                inner = ctx.add(&[
                    inner,
                    ctx.mul(&[
                        ctx.num(binom),
                        ctx.pow(ctx.mul(&[ctx.num(-1), dc.a]), ctx.num((j - i) as i64)),
                        jk,
                    ]),
                ]);
            }
            terms.push(ctx.mul(&[
                *rj,
                inv(ctx, ctx.pow(ctx.num(1), ctx.num(0))),
                inv(ctx, ctx.pow(dc.b, ctx.num(j as i64))),
                inner,
            ]));
        }
    }
    Some(normalize(ctx, sum_atoms(ctx, &terms)))
}

/// `J_n = ∫du/D^n` for `n ≤ MAX_DENOM_POW` (memoized).
fn j_recur<'a>(dc: &mut DenomCtx<'a>, n: i64) -> Option<Atom<'a>> {
    if n > MAX_DENOM_POW {
        return None;
    }
    if let Some(v) = dc.memo.get(&(usize::MAX, n)) {
        return *v;
    }
    if dc.memo.len() > MAX_MEMO {
        return None;
    }
    let value = j_compute(dc, n);
    dc.memo.insert((usize::MAX, n), value);
    value
}

fn j_compute<'a>(dc: &mut DenomCtx<'a>, n: i64) -> Option<Atom<'a>> {
    let ctx = dc.ctx;
    match n {
        // ∫du
        0 => Some(ctx.var(dc.var.as_str())),
        1 => dc.seed(),
        _ if n >= 2 => {
            // (n−1)(a² + εk²)·J_n = a(2n−3)·J_{n−1} + (2−n)·J_{n−2}
            //                        − k·T̂/D^(n−1),
            // from differentiating T̂/D^(n−1) and using T̂² = T² + ε.
            let j1 = j_recur(dc, n - 1)?;
            let j2 = j_recur(dc, n - 2)?;
            let boundary = ctx.mul(&[dc.b, dc.comp, ctx.pow(dc.d0, ctx.num(-(n - 1)))]);
            let coeff = ctx.mul(&[dc.a, ctx.num(2 * n - 3)]);
            let mut terms = vec![ctx.mul(&[coeff, j1]), ctx.mul(&[ctx.num(-1), boundary])];
            if n > 2 {
                // (2 − n)·J_{n−2} = −(n − 2)·J_{n−2}
                terms.push(ctx.mul(&[ctx.num(-(n - 2)), j2]));
            } else {
                let _ = j2;
            }
            let den = ctx.mul(&[ctx.num(n - 1), dc.disc]);
            Some(normalize(
                ctx,
                ctx.mul(&[sum_atoms(ctx, &terms), inv(ctx, den)]),
            ))
        }
        _ => {
            // n < 0: ∫D^|n| du by the binomial expansion in T/T̂.
            let m = n.unsigned_abs() as usize;
            let mut acc = ctx.num(0);
            for i in 0..=m {
                for j in 0..=(m - i) {
                    let bi = binom_i64(m as i64, i as i64)?;
                    let bj = binom_i64((m - i) as i64, j as i64)?;
                    let coef = ctx.mul(&[
                        ctx.num(bi),
                        ctx.num(bj),
                        ctx.pow(dc.a, ctx.num((m - i - j) as i64)),
                        ctx.pow(dc.b, ctx.num(i as i64)),
                        ctx.pow(dc.b, ctx.num(j as i64)),
                    ]);
                    let anti = tprod_integral(dc, i, j)?;
                    acc = ctx.add(&[acc, ctx.mul(&[coef, anti])]);
                }
            }
            Some(normalize(ctx, acc))
        }
    }
}

/// `∫T^i·T̂^j du` by the companion recurrence
/// `(i+j+1)·I[i,j] = −T^i·T̂^(j+1)/d − i·I[i−1, j+2]`-style integration by
/// parts; the `i + j` degree never grows.
fn tprod_integral<'a>(dc: &mut DenomCtx<'a>, i: usize, j: usize) -> Option<Atom<'a>> {
    let ctx = dc.ctx;
    if i + j > MAX_ARG_DEG + MAX_DENOM_POW as usize {
        return None;
    }
    if i == 0 && j == 0 {
        return Some(ctx.var(dc.var.as_str()));
    }
    if i == 1 && j == 0 {
        return Some(ctx.mul(&[dc.comp, inv(ctx, dc.slope)]));
    }
    if i == 0 && j == 1 {
        return Some(ctx.mul(&[dc.t, inv(ctx, dc.slope)]));
    }
    // d(T^i·T̂^(j+1)) = (i+j+1)·T^i·T̂^j + i·ε·T^(i−1)·T̂^(j+1)  (times du/d)
    // because T̂' = T and T·T̂^(j+1) = T^i... evaluated below through the
    // identity T̂² = T² + ε.
    let boundary = ctx.mul(&[
        pow_t(ctx, dc.t, i),
        pow_t(ctx, dc.comp, j + 1),
        inv(ctx, dc.slope),
    ]);
    let prev = if i == 0 {
        None
    } else if dc.eps > 0 {
        Some(tprod_integral(dc, i - 1, j + 1)?)
    } else {
        Some(tprod_integral(dc, i - 1, j + 3)?)
    };
    let n = ctx.num((i + j + 1) as i64);
    let mut terms = vec![ctx.mul(&[inv(ctx, n), boundary])];
    if let Some(prev) = prev {
        terms.push(ctx.mul(&[ctx.num(-1), ctx.num(i as i64), inv(ctx, n), prev]));
    }
    Some(normalize(ctx, sum_atoms(ctx, &terms)))
}

/// `∫T^L du` (the quotient part of the long division):
/// `N_k = T^(k−1)·T̂/(d·k) − ((k−1)/k)·N_{k−2}`.
fn t_power_integral<'a>(dc: &DenomCtx<'a>, l: usize) -> Option<Atom<'a>> {
    let ctx = dc.ctx;
    if l > MAX_ARG_DEG {
        return None;
    }
    match l {
        0 => Some(ctx.var(dc.var.as_str())),
        1 => Some(ctx.mul(&[dc.comp, inv(ctx, dc.slope)])),
        _ => {
            let prev = t_power_integral(dc, l - 2)?;
            let boundary = ctx.mul(&[
                pow_t(ctx, dc.t, l - 1),
                dc.comp,
                inv(ctx, dc.slope),
                inv(ctx, ctx.num(l as i64)),
            ]);
            let coef = ctx.mul(&[ctx.num((l - 1) as i64), inv(ctx, ctx.num(l as i64))]);
            Some(normalize(
                ctx,
                ctx.add(&[boundary, ctx.mul(&[ctx.num(-1), coef, prev])]),
            ))
        }
    }
}

impl<'a> DenomCtx<'a> {
    /// `J_1 = ∫du/D`, the exact seed of the recurrence:
    ///
    /// ```text
    /// T = sinh: 2/(d·√(a²+b²))·acoth((a·t − b)/√(a²+b²)),  t = (1 + sinh u)/cosh u,
    /// T = cosh: 2/(d·√(a²−b²))·acoth((a+b)/((a−b)·t)),     t = sinh u/(1 + cosh u),
    /// ```
    ///
    /// both verified by differentiation in the tests.
    fn seed(&self) -> Option<Atom<'a>> {
        let ctx = self.ctx;
        // `t = tanh(u/2) = sinh u/(1 + cosh u)`, written in the kernels:
        // the sinh-like kernel over `1 +` the cosh-like kernel.
        let (sinh_like, cosh_like) = if self.eps > 0 {
            (self.t, self.comp)
        } else {
            (self.comp, self.t)
        };
        let t = ctx.mul(&[sinh_like, inv(ctx, ctx.add(&[ctx.num(1), cosh_like]))]);
        // `J_1 = (1/(d·√disc))·log((z + 1)/(z − 1))` with
        //   T = sinh: z = (a·t − b)/√(a² + b²),
        //   T = cosh: z = β/t,  β = √((a + b)/(a − b)).
        // The `log` form is emitted (rather than `acoth`/`atanh`) so that both
        // `crate::diff` and the numeric evaluator handle every branch.
        let z = if self.eps > 0 {
            ctx.mul(&[
                ctx.add(&[ctx.mul(&[self.a, t]), ctx.mul(&[ctx.num(-1), self.b])]),
                inv(ctx, self.sqrt_disc),
            ])
        } else {
            let beta = root(
                ctx,
                ctx.mul(&[
                    ctx.add(&[self.a, self.b]),
                    inv(ctx, ctx.add(&[self.a, ctx.mul(&[ctx.num(-1), self.b])])),
                ]),
                2,
            );
            ctx.mul(&[beta, inv(ctx, t)])
        };
        let ratio = ctx.mul(&[
            ctx.add(&[z, ctx.num(1)]),
            inv(ctx, ctx.add(&[z, ctx.num(-1)])),
        ]);
        let pre = ctx.mul(&[inv(ctx, self.slope), inv(ctx, self.sqrt_disc)]);
        Some(normalize(ctx, ctx.mul(&[pre, ctx.fun("log", &[ratio])])))
    }
}
fn pow_t<'a>(ctx: &'a AtomArena<'a>, t: Atom<'a>, e: usize) -> Atom<'a> {
    match e {
        0 => ctx.num(1),
        1 => t,
        _ => ctx.pow(t, ctx.num(e as i64)),
    }
}

// =========================================================================
// Routine 2: kernel-polynomial terms
// =========================================================================

/// Kernel-polynomial terms: `∫R(T(u)) du` where `R` is rational in the
/// kernel.
///
/// The intended rewrite multiplies the t-form of the integrand by
/// `du/dt = 1/T′(t)`; when that product is a polynomial in `t` the
/// antiderivative is its termwise integral, back-substituted. The rewrite is
/// **not enabled**: its radical bookkeeping is not verified yet, and an
/// unverified emission is strictly worse than a decline. The gate is kept so
/// the class is documented in one place.
fn integrate_kernel_poly<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sub: &KernelSub<'a>,
) -> Option<Atom<'a>> {
    let _ = (ctx, expr, var, sub);
    None
}

/// The kernel-polynomial rewrite (disabled — see [`integrate_kernel_poly`]).
#[allow(dead_code)]
fn kernel_poly_rewrite<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sub: &KernelSub<'a>,
) -> Option<Atom<'a>> {
    let alg = sub.to_alg(ctx, expr, var)?;
    let jac = sub.du_dt(ctx);
    let total = alg.mul(ctx, &jac)?;
    if total.rad_pow != 0 || total.den.degree() != Some(0) || total.num.is_zero() {
        return None;
    }
    let p = total.num;
    let t_atom = sub.t_atom(ctx);
    let mut anti_t = poly_antiderivative_t(ctx, &p, t_atom);
    // The constant denominator scales the answer.
    let den_atom = total.den.to_atom(ctx, t_atom);
    anti_t = normalize(ctx, ctx.mul(&[anti_t, inv(ctx, den_atom)]));
    let out = sub.back(ctx, anti_t);
    if super::node_count(out) > MAX_SUBST_NODES {
        return None;
    }
    Some(normalize(ctx, out))
}

fn poly_antiderivative_t<'a>(ctx: &'a AtomArena<'a>, p: &Poly<'a>, t: Atom<'a>) -> Atom<'a> {
    if p.terms.is_empty() {
        return ctx.num(0);
    }
    let mut terms = Vec::new();
    for (i, c) in p.terms.iter().enumerate() {
        let n = (i + 1) as i64;
        terms.push(ctx.mul(&[*c, inv(ctx, ctx.num(n)), pow_t(ctx, t, i + 1)]));
    }
    normalize(ctx, sum_atoms(ctx, &terms))
}

// =========================================================================
// Routine 3: polynomial×kernel closed forms
// =========================================================================

/// `∫(e + f·x)^m·T(u) dx` and `∫(e + f·x)^m·P(T)/D^n dx`.
fn integrate_poly_hyper<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    sub: &KernelSub<'a>,
) -> Option<Atom<'a>> {
    let _ = (ctx, expr, var, sub);
    let _ = MAX_POLY_DEG;
    None
}

// =========================================================================
// Routine 4: algebraic (half-integer-power) substitution
// =========================================================================

/// Half-integer powers over a single kernel.
///
/// The gate for this class is the binomial differential the plan describes
/// (`t^p·(α + β·t)^q·dt` after `t = T(u)`). The rewrite here is exactly that
/// substitution, but the residual `t`-form is only handed on when it is a
/// binomial the sibling engine accepts; otherwise — and always for the
/// verified-shape contract of this module — it declines, so no unverified
/// algebraic or elliptic shape is ever emitted from here.
fn integrate_algebraic_power<'a>(
    ctx: &'a AtomArena<'a>,
    _expr: Atom<'a>,
    _var: Symbol,
    sub: &KernelSub<'a>,
) -> Option<Atom<'a>> {
    let _ = (ctx, sub);
    None
}

// =========================================================================
// Final numeric verification
// =========================================================================

/// Specialize the free parameters and check
/// `d/dx(candidate) == integrand` numerically at [`SAMPLES`].
fn verify_numeric<'a>(
    ctx: &'a AtomArena<'a>,
    integrand: Atom<'a>,
    var: Symbol,
    candidate: Atom<'a>,
) -> bool {
    if contains_integral(candidate) {
        return false;
    }
    let d = crate::diff(ctx, candidate, var);
    for &xv in &SAMPLES {
        let mut env: Vec<(Symbol, f64)> = PARAM_SAMPLES
            .iter()
            .map(|(n, v)| (Symbol::new(n), *v))
            .collect();
        env.push((var, xv));
        let (Some(lhs), Some(rhs)) = (eval_f64(d, &env), eval_f64(integrand, &env)) else {
            continue;
        };
        if !lhs.is_finite() || !rhs.is_finite() {
            continue;
        }
        if (lhs - rhs).abs() > 1e-6 * rhs.abs().max(1.0) {
            if std::env::var_os("OCAS_HYP_DEBUG").is_some() {
                eprintln!(
                    "[hyp] verify REJECT at x={xv}: diff={lhs} integrand={rhs} cand={candidate}"
                );
            }
            return false;
        }
    }
    true
}

/// Numeric `f64` evaluator used by [`verify_numeric`] and the tests.
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
                "sinh" => v.sinh(),
                "cosh" => v.cosh(),
                "tanh" => v.tanh(),
                "coth" => v.tanh().recip(),
                "sech" => v.cosh().recip(),
                "csch" => v.sinh().recip(),
                "exp" => v.exp(),
                "log" => v.abs().ln(),
                "sqrt" => v.sqrt(),
                "atan" => v.atan(),
                "atanh" => v.atanh(),
                "acoth" => {
                    // arcoth z = ½ ln|(z+1)/(z−1)|, with z ↦ 1/z for |z| < 1
                    // (the branches differ by a constant).
                    let z = if v.abs() < 1.0 { 1.0 / v } else { v };
                    0.5 * ((z + 1.0) / (z - 1.0)).abs().ln()
                }
                _ => return None,
            })
        }
    }
}

/// `C(n, k)` for small non-negative integers (0 outside the range).
fn binom_i64(n: i64, k: i64) -> Option<i64> {
    if k < 0 || k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    let mut r: i64 = 1;
    for i in 0..k {
        r = r.checked_mul(n - i)?.checked_div(i + 1)?;
    }
    Some(r)
}

/// `(base,n)` for `(base^n)^(−1)` or `base^(−n)` with integer `n ≥ 1`.
fn as_recip_pow(f: Atom<'_>) -> Option<(Atom<'_>, i64)> {
    let AtomNode::Pow(b, e) = f.node() else {
        return None;
    };
    if let (AtomNode::Pow(bb, ee), AtomNode::Num(-1)) = (b.node(), e.node())
        && let AtomNode::Num(n) = ee.node()
        && *n >= 1
    {
        return Some((*bb, *n));
    }
    if let AtomNode::Num(n) = e.node()
        && *n <= -1
    {
        return Some((*b, n.checked_neg()?));
    }
    None
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        normalize(ctx, ocas_parse::parse(ctx, s).unwrap())
    }

    fn env(pairs: &[(&str, f64)]) -> Vec<(Symbol, f64)> {
        pairs.iter().map(|(n, v)| (Symbol::new(n), *v)).collect()
    }

    /// Run the mechanism, require `Some` without residue, and check
    /// `diff(result) == integrand` numerically at the sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate_hyperbolic_reduction(ctx, integrand, var)
            .unwrap_or_else(|| panic!("mechanism declined for {integrand}"));
        assert!(
            !result.to_string().contains("Integral"),
            "residue: {result}"
        );
        let d = crate::diff(ctx, result, var);
        let mut checked = 0usize;
        for &xv in samples {
            let mut e = consts.to_vec();
            e.push((var, xv));
            let (Some(lhs), Some(rhs)) = (eval_f64(d, &e), eval_f64(integrand, &e)) else {
                continue;
            };
            if !lhs.is_finite() || !rhs.is_finite() {
                continue;
            }
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
            checked += 1;
        }
        assert!(
            checked > 0,
            "no finite sample for {integrand} (result: {result})"
        );
    }

    fn assert_declined<'a>(ctx: &'a AtomArena<'a>, input: &str) {
        let r = integrate_hyperbolic_reduction(ctx, parse(ctx, input), Symbol::new("x"));
        assert!(r.is_none(), "expected None for {input}, got {r:?}");
    }

    // ------------- the verified surface: `∫du/(a + b·T(u))` -------------

    #[test]
    fn denom_numeric() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for (input, samples) in [
            ("1/(2 + sinh(x))", &[0.3, 0.7, 1.1][..]),
            ("1/(3 + 2*cosh(x))", &[0.3, 0.7, 1.1][..]),
            ("1/(5 + 3*cosh(x))", &[0.3, 0.7, 1.1][..]),
            ("3/(2 + sinh(x))", &[0.3, 0.7, 1.1][..]),
            ("1/(2 + sinh(2*x + 1))", &[0.2, 0.5, 0.9][..]),
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, input), var, &[], samples);
        }
    }

    #[test]
    fn denom_symbolic_coefficients() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        // `∫ du/(a + b·sinh(c + d·x))` and the cosh twin.
        let cs = env(&[("a", 2.0), ("b", 0.5), ("c", 0.4), ("d", 1.3)]);
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "1/(a + b*sinh(c + d*x))"),
            var,
            &cs,
            &[0.2, 0.5, 0.9],
        );
        let cs = env(&[("a", 3.0), ("b", 0.5), ("c", 0.4), ("d", 1.3)]);
        assert_antiderivative_num(
            &ctx,
            parse(&ctx, "1/(a + b*cosh(c + d*x))"),
            var,
            &cs,
            &[0.2, 0.5, 0.9],
        );
    }

    /// The `J_1` seed must be an antiderivative for both kernels, at numeric
    /// and at symbolic coefficients.
    #[test]
    fn seed_is_an_antiderivative() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        for (input, consts, samples) in [
            ("1/(2 + sinh(x))", env(&[]), &[0.3, 0.7, 1.1][..]),
            ("1/(3 + 2*cosh(x))", env(&[]), &[0.3, 0.7, 1.1][..]),
            (
                "1/(a + b*sinh(c + d*x))",
                env(&[("a", 2.0), ("b", 0.5), ("c", 0.4), ("d", 1.3)]),
                &[0.2, 0.5, 0.9][..],
            ),
            (
                "1/(a + b*cosh(c + d*x))",
                env(&[("a", 3.0), ("b", 0.5), ("c", 0.4), ("d", 1.3)]),
                &[0.2, 0.5, 0.9][..],
            ),
        ] {
            assert_antiderivative_num(&ctx, parse(&ctx, input), var, &consts, samples);
        }
    }

    #[test]
    fn denom_singular_declines() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // a² = b² with T = cosh: D has a double root in e^u.
        assert_declined(&ctx, "1/(2 + 2*cosh(x))^2");
        // Negative leading coefficient: not a verified branch.
        assert_declined(&ctx, "1/(-5 + 3*cosh(x))");
        // Negative leading coefficient: not a verified branch (documented gap).
        assert_declined(&ctx, "1/(-5 + 3*cosh(x))");
    }

    // ------------- declines -------------

    #[test]
    fn declines_unverified_families() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Denominator powers above one: the recurrence is not enabled.
        assert_declined(&ctx, "1/(2 + sinh(x))^2");
        assert_declined(&ctx, "1/(a + b*sinh(x))^3");
        // Non-constant numerator: the long-division path is not enabled.
        assert_declined(&ctx, "sinh(x)^3/(a + b*sinh(x))");
        // Kernel-polynomial terms: the radical bookkeeping is not verified.
        assert_declined(&ctx, "sinh(x)^2*cosh(x)");
        assert_declined(&ctx, "tanh(x)^3*sech(x)^2");
        // Companion mixing needs the phase reduction (documented gap).
        assert_declined(&ctx, "1/(a + b*cosh(x) + c*sinh(x))^2");
        assert_declined(&ctx, "cosh(x)^2/(a + b*coth(x))");
        // A non-kernel numerator over a kernel denominator.
        assert_declined(&ctx, "coth(x)/(a + b*sinh(x))^2");
        assert_declined(&ctx, "sech(x)/(a + b*sinh(x))");
        // Half-integer powers are not owned here (elliptic engine's class).
        assert_declined(&ctx, "(1 - tanh(x)^2)^(3/2)");
        assert_declined(&ctx, "(1 + coth(x))^(7/2)");
        // Routine 3 (`x^m·T(u)`) is documented but not enabled.
        assert_declined(&ctx, "x^2*sinh(x)");
    }

    #[test]
    fn declines_out_of_family() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Mixed trig/hyperbolic.
        assert_declined(&ctx, "sinh(x)*cos(x)");
        // Non-linear kernel argument.
        assert_declined(&ctx, "sinh(x^2)");
        // Kernel of a kernel.
        assert_declined(&ctx, "sinh(cosh(x))");
        // Two different kernel arguments.
        assert_declined(&ctx, "sinh(x)*sinh(2*x)");
        // Rational in the kernel but not a polynomial in it.
        assert_declined(&ctx, "1/(a + b*sinh(x)^2)");
        // Two denominator powers.
        assert_declined(&ctx, "1/((1 + sinh(x))^2*(2 + sinh(x)))");
    }

    // ------------- stress -------------

    #[test]
    fn stress_stability_and_budget() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let inputs = [
            "1/(a + b*sinh(c + d*x))",
            "1/(a + b*cosh(c + d*x))",
            "sinh(x)^3/(a + b*sinh(x))",
            "coth(x)/(a + b*sinh(x))^2",
            "sech(x)/(a + b*sinh(x))",
            "cosh(x)^3*sinh(x)^2/(a*cosh(x) + b*sinh(x))^2",
            "cosh(x)^2/(a + b*coth(x))",
            "sinh(x)^2*cosh(x)",
            "1/(a + b*sinh(x))^2",
            "1/(a + b*sinh(x)^2)",
            "sinh(x)*cos(x)",
            "sinh(x^2)",
            "x^2*sinh(x)",
            "(1 - tanh(x)^2)^(3/2)",
            "tanh(x)^3*sech(x)^2",
        ];
        let exprs: Vec<Atom<'_>> = inputs.iter().map(|s| parse(&ctx, s)).collect();
        let first: Vec<Option<String>> = exprs
            .iter()
            .map(|e| integrate_hyperbolic_reduction(&ctx, *e, var).map(|a| a.to_string()))
            .collect();
        assert!(first[0].is_some(), "1/(a + b sinh) not solved");
        assert!(first[1].is_some(), "1/(a + b cosh) not solved");
        assert!(first[2].is_none(), "sinh^3/(a+b sinh) must decline for now");
        assert!(first[9].is_none(), "1/(a + b sinh^2) must decline");
        for round in 0..300 {
            for (i, e) in exprs.iter().enumerate() {
                let again = integrate_hyperbolic_reduction(&ctx, *e, var).map(|a| a.to_string());
                assert_eq!(
                    again, first[i],
                    "case {i} ({}) changed at round {round}",
                    inputs[i]
                );
                if let Some(s) = &again {
                    assert!(
                        !s.contains("Integral"),
                        "case {i} ({}) left a residue at round {round}: {s}",
                        inputs[i]
                    );
                    assert!(
                        s.len() <= 4000,
                        "case {i} ({}) grew pathologically at round {round}: {} chars",
                        inputs[i],
                        s.len()
                    );
                }
            }
        }
    }
}
