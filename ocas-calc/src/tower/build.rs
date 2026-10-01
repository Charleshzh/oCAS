//! Construction of elementary extension towers from expressions.
//!
//! [`build_tower`] inspects an integrand and builds the differential field
//! tower `ℚ(x, t₁, …, tₙ)` in which it lives, where each `tᵢ` is a
//! logarithm, an exponential, or a constant over the field below. It also
//! computes the derivative of each generator, which the Risch algorithm
//! needs, and returns the integrand rewritten over the merged generators.
//!
//! Limitations (the caller falls back to other integrators):
//!
//! - only `log` / `exp` function applications are admitted (trigonometric
//!   integrands are rewritten into exponentials before this entry point);
//! - algebraic functions such as `√x` (non-integer exponents) are
//!   rejected.
//!
//! Algebraically dependent generators (`log(x)` with `log(2x)`, `exp(x)`
//! with `exp(x+1)`, `exp(u)` with `exp(−u)`) are **merged** rather than
//! rejected since 0.28.0 — see [`super::merge`].

use ocas_atom::walk::collect_funs;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::convert::atom_to_rational_extended;
use super::elem::{KElem, KPoly};
use super::merge;

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
    /// The integrand rewritten over the merged generators.
    pub expr: Atom<'a>,
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
///
/// 0.28.0: algebraically dependent generators are **merged** rather than
/// rejected (see [`super::merge`]); the tower carries the rewritten
/// integrand in [`Tower::expr`], which is the form the caller must convert
/// into the field.
pub(crate) fn build_tower<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Tower<'a>> {
    let x = ctx.var(var.as_str());
    if !only_integer_powers(expr) {
        return None;
    }

    // Pass 1: collect and merge generators (innermost first).
    let mut gens: Vec<GenInfo<'a>> = Vec::new();
    // Expressions carrying the imaginary unit (the trigonometric rewrite's
    // exponential form) are still rejected: the constant generator `I` is
    // accepted by this pass but its tower then unlocks the trig-exp Risch
    // path, which returns `exp(I·…)` partial results for shapes the
    // trigonometric mechanisms solve in real form (measured 0.28.0: the
    // `rules`/`trig_kernel` real-form regression suites fail once it is
    // enabled). Enabling it needs a realification pass of its own and is
    // recorded as a follow-up, not part of this wave.
    if contains_var(expr, "I") {
        return None;
    }
    let mut current = expr;
    loop {
        let mut rewritten = false;
        for (name, app) in collect_funs(current) {
            // Already a generator, including the constant atoms a previous
            // merge introduced.
            if gens.iter().any(|g| g.atom == app) {
                continue;
            }
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
            if let Some(merge) = merge::merge_candidate(ctx, kind, arg, app, &gens, var) {
                if let Some(constant) = merge.constant
                    && !gens.iter().any(|g| g.atom == constant)
                {
                    gens.push(GenInfo {
                        kind: GenKind::Constant,
                        atom: constant,
                        arg: constant,
                        dt: KElem::zero(0),
                    });
                }
                current = merge::substitute_atom(ctx, current, app, merge.replacement);
                rewritten = true;
                // The rewrite may have introduced new function atoms;
                // rescan from the top.
                break;
            }
            // `log`/`exp` of a numeric constant should be a plain number.
            if is_rational_constant(arg) {
                return None;
            }
            gens.push(GenInfo {
                kind,
                atom: app,
                arg,
                // Filled in during pass 2; placeholder is never read before then.
                dt: KElem::zero(0),
            });
        }
        if !rewritten {
            break;
        }
    }

    // Pass 2: derivatives with the final variable count.
    let n = 1 + gens.len();
    for i in 0..gens.len() {
        let (done, rest) = gens.split_at_mut(i);
        let g = &mut rest[0];
        if g.kind == GenKind::Constant {
            // Constant generators (`I`, `log(2)`, `exp(1)`, …) are
            // `D t = 0` without any argument conversion: their argument may
            // not even live in the field below.
            g.dt = KElem::zero(n);
            continue;
        }
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

    Some(Tower {
        x,
        expr: current,
        gens,
    })
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
    fn tower_merges_dependent_logs() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // log(x) + log(2x): dependent. 0.28.0 merges log(2x) into
        // log(x) + log(2), with `log(2)` a new constant generator.
        let expr = ctx.add(&[
            ctx.fun("log", &[x]),
            ctx.fun("log", &[ctx.mul(&[ctx.num(2), x])]),
        ]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 2);
        assert_eq!(tower.gens[0].kind, GenKind::Log);
        assert_eq!(tower.gens[1].kind, GenKind::Constant);
        assert_eq!(tower.gens[1].atom.to_string(), "log(2)");
        assert!(tower.gens[1].dt.is_zero());
        // The integrand is rewritten over the merged generators.
        assert_eq!(tower.expr.to_string(), "(log(x)) + ((log(x)) + (log(2)))");
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
    fn tower_merges_reciprocal_exponentials() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x)·exp(−x) = 1: the hyperbolic rewrites produce this pair.
        let expr = ctx.mul(&[
            ctx.fun("exp", &[x]),
            ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), x])]),
        ]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 1);
        assert_eq!(tower.gens[0].kind, GenKind::Exp);
        assert_eq!(tower.expr.to_string(), "(exp(x))*((exp(x))^-1)");
    }

    #[test]
    fn tower_merges_shifted_exponentials() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(x)·exp(x+1): dependent (ratio e). 0.28.0 merges the shifted
        // generator with `exp(1)` as a constant generator.
        let expr = ctx.mul(&[
            ctx.fun("exp", &[x]),
            ctx.fun("exp", &[ctx.add(&[x, ctx.num(1)])]),
        ]);
        let tower = build_tower(&ctx, expr, sym("x")).expect("tower");
        assert_eq!(tower.gens.len(), 2);
        assert_eq!(tower.gens[1].kind, GenKind::Constant);
        assert_eq!(tower.gens[1].atom.to_string(), "exp(1)");
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
    fn tower_rejects_the_imaginary_unit_for_now() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // exp(I·x) — the trigonometric rewrite's exponential form. The
        // constant generator `I` is deliberately still rejected (0.27.3
        // behaviour): unlocking it also unlocks the trig-exp Risch path,
        // whose `exp(I·…)` answers regress the real-form suites. The merge
        // machinery below is covered by the exponential/log tests instead.
        let expr = ctx.fun("exp", &[ctx.mul(&[ctx.var("I"), x])]);
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
