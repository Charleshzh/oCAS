//! Construction of elementary extension towers from expressions.
//!
//! [`build_tower`] inspects an integrand and builds the differential field
//! tower `ℚ(x, t₁, …, tₙ)` in which it lives, where each `tᵢ` is a
//! logarithm or exponential over the field below. It also computes the
//! derivative of each generator, which the Risch algorithm needs.
//!
//! Limitations (the caller falls back to other integrators):
//!
//! - only `log` / `exp` function applications are admitted (trigonometric
//!   integrands are rewritten into exponentials before this entry point);
//! - algebraic functions such as `√x` (non-integer exponents) are
//!   rejected;
//! - algebraically dependent generators (e.g. `log(x)` and `log(2x)`, or
//!   `exp(x)` and `exp(x+1)`) are rejected rather than merged.

use ocas_atom::walk::collect_funs;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::convert::atom_to_rational_extended;
use super::elem::{KElem, KPoly};

/// Kind of a tower generator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum GenKind {
    /// A constant symbol (e.g. the imaginary unit), `D t = 0`.
    Constant,
    Log,
    Exp,
}

/// A tower generator together with its derivative.
pub(crate) struct GenInfo<'a> {
    /// Whether this generator is a logarithm or an exponential.
    pub kind: GenKind,
    /// The application atom, e.g. `log(u)`.
    pub atom: Atom<'a>,
    /// The argument `u`.
    pub arg: Atom<'a>,
    /// `D(generator)` as an element of the full tower field.
    pub dt: KElem,
}

/// A differential extension tower `ℚ(x, t₁, …, tₙ)`.
pub(crate) struct Tower<'a> {
    /// Integration variable atom.
    pub x: Atom<'a>,
    /// Generators from bottom to top; `gens[i]` ↔ variable index `i + 1`.
    pub gens: Vec<GenInfo<'a>>,
}

impl<'a> Tower<'a> {
    /// Number of polynomial variables (`1 + #generators`).
    pub fn n_vars(&self) -> usize {
        1 + self.gens.len()
    }

    /// All generator atoms `[x, t₁, …, tₙ]` in variable order.
    pub fn gen_atoms(&self) -> Vec<Atom<'a>> {
        let mut v = Vec::with_capacity(self.n_vars());
        v.push(self.x);
        v.extend(self.gens.iter().map(|g| g.atom));
        v
    }
}

/// Full derivation of a field element w.r.t. the tower derivation.
///
/// `gens` must be the prefix of tower generators the element depends on
/// (`gens[i]` ↔ variable index `i + 1`); `D x = 1`.
pub(crate) fn tower_diff(e: &KElem, gens: &[GenInfo]) -> KElem {
    let mut acc = e.partial_deriv(0);
    for (i, g) in gens.iter().enumerate() {
        acc = acc.add(&e.partial_deriv(i + 1).mul(&g.dt));
    }
    acc
}

/// Full derivation of a `k[t]` polynomial: `D(Σ aᵢ tⁱ) = Σ D(aᵢ) tⁱ +
/// (dp/dt)·Dt`. `gens` is the prefix below the top; `dt_top` is `D t`.
pub(crate) fn tower_diff_kpoly(p: &KPoly, gens: &[GenInfo], dt_top: &KElem) -> KPoly {
    let derived = KPoly {
        top: p.top,
        coeffs: p.coeffs.iter().map(|c| tower_diff(c, gens)).collect(),
        n_vars: p.n_vars,
    };
    derived.add(&p.derivative_dt().mul_kelem(dt_top))
}

/// Build the extension tower for `expr` over the integration variable
/// `var`, or `None` when the expression is not elementary-admissible (see
/// module docs).
pub(crate) fn build_tower<'a>(
    _ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Tower<'a>> {
    let x = _ctx.var(var.as_str());
    if !only_integer_powers(expr) {
        return None;
    }

    // Pass 1: collect and validate generators (innermost first).
    let mut gens: Vec<GenInfo<'a>> = Vec::new();
    // The imaginary unit, when present, becomes a constant generator.
    if contains_var(expr, "I") {
        gens.push(GenInfo {
            kind: GenKind::Constant,
            atom: _ctx.var("I"),
            arg: _ctx.var("I"),
            dt: KElem::zero(0),
        });
    }
    for (name, app) in collect_funs(expr) {
        let kind = match name.as_str() {
            "log" => GenKind::Log,
            "exp" => GenKind::Exp,
            _ => return None,
        };
        let args = app.children();
        if args.len() != 1 {
            return None;
        }
        let arg = args[0];
        if is_rational_constant(arg) {
            // log/exp of a constant should be a plain number; not a tower.
            return None;
        }
        for g in &gens {
            if algebraically_dependent(kind, arg, g) {
                return None;
            }
        }
        gens.push(GenInfo {
            kind,
            atom: app,
            arg,
            // Filled in during pass 2; placeholder is never read before then.
            dt: KElem::zero(0),
        });
    }

    // Pass 2: derivatives with the final variable count.
    let n = 1 + gens.len();
    for i in 0..gens.len() {
        let (done, rest) = gens.split_at_mut(i);
        let g = &mut rest[0];
        let mut prefix_atoms = Vec::with_capacity(i + 1);
        prefix_atoms.push(x);
        prefix_atoms.extend(done.iter().map(|d| d.atom));
        let u_rf = atom_to_rational_extended(g.arg, &prefix_atoms, n)?;
        let u_k = KElem::new(u_rf.numerator, u_rf.denominator);
        let du = tower_diff(&u_k, done);
        g.dt = match g.kind {
            GenKind::Constant => KElem::zero(n),
            GenKind::Log => du.div(&u_k)?,
            GenKind::Exp => du.mul(&KElem::var(i + 1, n)),
        };
    }

    Some(Tower { x, gens })
}

/// Whether `atom` is a rational constant (contains no variables at all).
fn is_rational_constant(atom: Atom) -> bool {
    atom_to_rational_extended(atom, &[], 0).is_some()
}

/// Whether `atom` contains the variable with the given name anywhere.
fn contains_var(atom: Atom, name: &str) -> bool {
    match atom.node() {
        AtomNode::Var(s) => s.as_str() == name,
        AtomNode::Num(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|a| contains_var(*a, name))
        }
        AtomNode::Pow(b, e) => contains_var(*b, name) || contains_var(*e, name),
    }
}

/// Reject expressions containing non-integer powers of non-constant bases
/// (algebraic functions such as `√x`).
fn only_integer_powers(atom: Atom) -> bool {
    if let Some((base, exp)) = atom.binary_children() {
        let ok = matches!(exp.node(), AtomNode::Num(_))
            || (is_rational_constant(base) && is_rational_constant(exp));
        return ok && only_integer_powers(base) && only_integer_powers(exp);
    }
    atom.children().iter().all(|c| only_integer_powers(*c))
}

/// Conservative algebraic-dependence check between a candidate generator
/// `(kind, arg)` and an existing one.
fn algebraically_dependent(kind: GenKind, arg: Atom, existing: &GenInfo) -> bool {
    // log(exp(v)) or exp(log(v)) collapse to v.
    if arg == existing.atom {
        return true;
    }
    let u = arg;
    let v = existing.arg;
    match (kind, existing.kind) {
        (GenKind::Log, GenKind::Log) => {
            // log(u) ∓ log(v) constant  ⟺  u/v or u·v constant.
            let ratio = is_rational_constant_div(u, v);
            let product = is_rational_constant_mul(u, v);
            // u = v^k or v = u^k for small integer k (e.g. log(x^2)).
            let powers = [2, 3, -2, -3]
                .iter()
                .any(|&k| pow_int_eq(u, v, k) || pow_int_eq(v, u, k));
            ratio || product || powers
        }
        (GenKind::Exp, GenKind::Exp) => {
            // exp(u)/exp(v) constant  ⟺  u - v constant.
            // exp(u)·exp(v) constant  ⟺  u + v constant (reciprocal pair).
            is_rational_constant_sub(u, v) || is_rational_constant_sum(u, v)
        }
        _ => false,
    }
}

/// Whether `u + v` is a rational constant, i.e. `u = -v + c`.
fn is_rational_constant_sum(u: Atom, v: Atom) -> bool {
    // Direct negation: u == -v (encoded as Mul with a -1 factor).
    if let AtomNode::Mul(factors) = u.node()
        && factors.len() == 2
        && matches!(factors[0].node(), AtomNode::Num(-1))
        && factors[1] == v
    {
        return true;
    }
    if let AtomNode::Mul(factors) = v.node()
        && factors.len() == 2
        && matches!(factors[0].node(), AtomNode::Num(-1))
        && factors[1] == u
    {
        return true;
    }
    // u = -v + c via Add shape: u = (-v) + c or v = (-u) + c.
    match u.node() {
        AtomNode::Add(args) if args.len() == 2 => {
            let neg_v = matches!(args[0].node(), AtomNode::Mul(f)
                if f.len() == 2 && matches!(f[0].node(), AtomNode::Num(-1)) && f[1] == v);
            if (neg_v && is_rational_constant(args[1]))
                || (args[0] == v && is_rational_constant(args[1]))
            {
                return true;
            }
        }
        _ => {}
    }
    false
}

// The dependence helpers below work on raw atoms; constant detection uses
// the rational converter with zero generators.

fn is_rational_constant_div(u: Atom, v: Atom) -> bool {
    // u/v — build via raw node inspection instead of an arena: walk both.
    // We cannot construct new atoms here (no ctx), so compare structurally:
    // u/v constant ⇔ u = c·v for a rational c — check via Mul shape.
    structurally_proportional(u, v)
}

fn is_rational_constant_mul(u: Atom, v: Atom) -> bool {
    // u·v constant ⇔ u = c/v — only cheap to detect when v = 1/w and
    // u = c·w; covered by structural proportionality on inverses.
    matches!(v.node(), AtomNode::Pow(b, e) if matches!(e.node(), AtomNode::Num(-1)) && structurally_proportional(u, *b))
}

fn is_rational_constant_sub(u: Atom, v: Atom) -> bool {
    // u - v constant ⇔ u = v + c. Detect via Add shape: (v + c) or (c + v).
    match u.node() {
        AtomNode::Add(args) if args.len() == 2 => {
            (args[0] == v && is_rational_constant(args[1]))
                || (args[1] == v && is_rational_constant(args[0]))
        }
        _ => false,
    }
}

/// Whether `u == v^k` structurally for integer `k` (hash-consed equality).
fn pow_int_eq(u: Atom, v: Atom, k: i64) -> bool {
    matches!(u.node(), AtomNode::Pow(b, e) if *b == v && matches!(e.node(), AtomNode::Num(n) if *n == k))
}

/// Whether `u = c·v` for a rational constant `c` (allowing an extra
/// constant factor on either side).
fn structurally_proportional(u: Atom, v: Atom) -> bool {
    if u == v {
        return true;
    }
    strip_const_factor(u) == Some(v) || strip_const_factor(v) == Some(u)
}

/// If `atom` is `c·w` with `c` a rational constant, return `w`.
fn strip_const_factor(atom: Atom) -> Option<Atom> {
    if let AtomNode::Mul(args) = atom.node()
        && args.len() == 2
    {
        if is_rational_constant(args[0]) {
            return Some(args[1]);
        }
        if is_rational_constant(args[1]) {
            return Some(args[0]);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_domain::{Domain, Rational, RationalDomain};

    use super::*;

    fn sym(name: &str) -> Symbol {
        Symbol::new(name)
    }

    #[test]
    fn tower_single_log() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[ctx.fun("log", &[x]), ctx.num(1)]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 1);
        assert_eq!(tower.gens[0].kind, GenKind::Log);
        // D log(x) = 1/x
        let one_over_x = KElem::one(2).div(&KElem::var(0, 2)).unwrap();
        assert!(tower.gens[0].dt.eq_cross(&one_over_x));
    }

    #[test]
    fn tower_nested_exp_log() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x·log(x)): tower [x, log(x), exp(x·log(x))]
        let log_x = ctx.fun("log", &[x]);
        let expr = ctx.fun("exp", &[ctx.mul(&[x, log_x])]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 2);
        assert_eq!(tower.gens[0].kind, GenKind::Log);
        assert_eq!(tower.gens[1].kind, GenKind::Exp);
        // D exp(x·log(x)) = (log(x) + 1)·t₂
        let log_var = KElem::var(1, 3);
        let exp_var = KElem::var(2, 3);
        let expect = log_var.add(&KElem::one(3)).mul(&exp_var);
        assert!(tower.gens[1].dt.eq_cross(&expect));
    }

    #[test]
    fn tower_rejects_dependent_logs() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // log(x) + log(2x): algebraically dependent.
        let expr = ctx.add(&[
            ctx.fun("log", &[x]),
            ctx.fun("log", &[ctx.mul(&[ctx.num(2), x])]),
        ]);
        assert!(build_tower(&ctx, expr, sym("x")).is_none());
    }

    #[test]
    fn tower_allows_independent_logs() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // log(x) + log(x+1): independent generators.
        let expr = ctx.add(&[
            ctx.fun("log", &[x]),
            ctx.fun("log", &[ctx.add(&[x, ctx.num(1)])]),
        ]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 2);
    }

    #[test]
    fn tower_rejects_dependent_exps() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x)·exp(x+1): dependent (ratio e).
        let expr = ctx.mul(&[
            ctx.fun("exp", &[x]),
            ctx.fun("exp", &[ctx.add(&[x, ctx.num(1)])]),
        ]);
        assert!(build_tower(&ctx, expr, sym("x")).is_none());
    }

    #[test]
    fn tower_rejects_trig_and_sqrt_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        assert!(build_tower(&ctx, ctx.fun("sin", &[x]), sym("x")).is_none());
        // x^(1/2): algebraic function.
        let sqrt_x = ctx.pow(x, ctx.pow(ctx.num(2), ctx.num(-1)));
        assert!(build_tower(&ctx, sqrt_x, sym("x")).is_none());
    }

    #[test]
    fn tower_rejects_log_of_constant() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[ctx.fun("log", &[ctx.num(3)]), x]);
        assert!(build_tower(&ctx, expr, sym("x")).is_none());
    }

    #[test]
    fn tower_diff_polynomial() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[x]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        // D(x·t₁) where t₁ = exp(x): = t₁ + x·t₁
        let xt = KElem::var(0, 2).mul(&KElem::var(1, 2));
        let d = tower_diff(&xt, &tower.gens);
        let t = KElem::var(1, 2);
        let expect = t.add(&KElem::var(0, 2).mul(&t));
        assert!(d.eq_cross(&expect));
        let _ = RationalDomain.one();
        let _ = Rational::new(1, 1);
    }
}
