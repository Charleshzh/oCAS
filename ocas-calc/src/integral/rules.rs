//! Rule-table integration engine (0.27.0).
//!
//! A table-driven extension of the integration pipeline: a library of
//! standard-calculus rules (own ruleset, structurally Rubi-inspired but with
//! original rule text) is applied when the rational/Risch/special/heuristic
//! stages all fail, extending the breadth of elementary antiderivatives the
//! engine can produce — e.g. `sin(x)^n` reductions, `x^m*exp(a*x)` reductions
//! and linear-argument forms.
//!
//! Rules are declared as [`RuleSpec`]s in a `static` table. Patterns use the
//! literal integration variable `x`; [`build_rule_table`] bakes the actual
//! variable in with word-boundary replacement and parses the patterns into a
//! head-indexed [`RuleTable`]. Templates may leave residual `Integral(g, x)`
//! terms, which [`integrate_rules`] resolves by recursively integrating `g`
//! with a fresh structural depth but an incremented rule-reduction budget
//! ([`MAX_RULE_DEPTH`]), so reduction formulas terminate.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::matcher::{Bindings, MatchValue, match_pattern};
use ocas_rewrite::pattern::Pattern;
use ocas_rewrite::rules::{HeadKey, Rule, RuleTable, head_of};

/// Reduction-formula recursion budget, independent of the structural depth
/// limit (`MAX_DEPTH` in mod.rs): `sin(x)^n`-style reductions must not be
/// truncated at structural depth 8, and the budget also stops self-loops.
pub(crate) const MAX_RULE_DEPTH: usize = 64;

/// Predicate over a successful match's bindings, with the integration
/// variable. Used to constrain free parameters (e.g. `n_` must be an integer
/// ≠ -1).
pub(crate) type Pred = for<'a> fn(&Bindings<'a>, Symbol) -> bool;

/// A single integration rule declaration.
///
/// - [`RuleSpec::Template`]: pattern/template strings (`Rule::from_template`).
/// - [`RuleSpec::Closure`]: hand-written replacement builder for shapes
///   templates cannot express (binomial expansion, odd-power peeling, ...).
pub(crate) enum RuleSpec {
    Template {
        pat: &'static str,
        tmpl: &'static str,
        cond: Option<Pred>,
    },
    Closure {
        pat: &'static str,
        f: for<'a> fn(&'a AtomArena<'a>, &Bindings<'a>, Symbol) -> Option<Atom<'a>>,
        cond: Option<Pred>,
    },
}

/// A compiled closure rule: pattern + optional predicate + builder.
struct ClosureRule<'a> {
    head: HeadKey,
    pattern: Pattern<'a>,
    cond: Option<Pred>,
    f: for<'b> fn(&'b AtomArena<'b>, &Bindings<'b>, Symbol) -> Option<Atom<'b>>,
}

/// Compiled rule table for one integration variable.
pub(crate) struct IntegralRuleTable<'a> {
    table: RuleTable<'a>,
    closures: Vec<ClosureRule<'a>>,
}

/// Unwrap the quadratic under a square root: handles `sqrt(q)`, `sqrt(q)^-1`,
/// `q^(1/2)` and `q^(-1/2)` (normalized forms).
fn sqrt_quadratic_base<'a>(expr: Atom<'a>) -> Option<Atom<'a>> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "sqrt" && args.len() == 1 => Some(args[0]),
        AtomNode::Pow(base, exp) => match (base.node(), exp.node()) {
            (AtomNode::Fun(name, args), AtomNode::Num(-1))
                if name.as_str() == "sqrt" && args.len() == 1 =>
            {
                Some(args[0])
            }
            // q^(2^-1) or q^(-1*2^-1)
            (_, AtomNode::Pow(e2, e3)) => {
                if matches!(e2.node(), AtomNode::Num(2)) && matches!(e3.node(), AtomNode::Num(-1)) {
                    Some(*base)
                } else {
                    None
                }
            }
            (_, AtomNode::Mul(args)) => {
                let has_neg_half = args.iter().any(|a| matches!(a.node(), AtomNode::Num(-1)))
                    && args.iter().any(|a| {
                        matches!(a.node(), AtomNode::Pow(b, e)
                        if matches!(b.node(), AtomNode::Num(2))
                            && matches!(e.node(), AtomNode::Num(-1)))
                    });
                if has_neg_half && args.len() == 2 {
                    Some(*base)
                } else {
                    None
                }
            }
            _ => None,
        },
        _ => None,
    }
}

/// Extract `(a, b, c)` from a quadratic square-root form
/// `√(a·x² + b·x + c)` (any of the four parse shapes), with implicit
/// coefficients made explicit (`x²` → `a = 1`, missing terms → `0`).
///
/// The square root may sit anywhere in `expr` (e.g. `1/(x·√(...))` is
/// `(x·√(...))^-1` after normalization); the first qualifying node is
/// returned.
///
/// Shared by the Euler substitution (heuristic) and rule family G.
pub(crate) fn quadratic_coeffs<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    if let Some(inner) = sqrt_quadratic_base(expr) {
        if let Some(coeffs) = quadratic_from_inner(ctx, inner, var) {
            return Some(coeffs);
        }
    }
    // Recurse into children.
    match expr.node() {
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                if let Some(coeffs) = quadratic_coeffs(ctx, *a, var) {
                    return Some(coeffs);
                }
            }
        }
        AtomNode::Pow(base, exp) => {
            if let Some(coeffs) = quadratic_coeffs(ctx, *base, var) {
                return Some(coeffs);
            }
            if let Some(coeffs) = quadratic_coeffs(ctx, *exp, var) {
                return Some(coeffs);
            }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => {}
    }
    None
}

/// Parse a quadratic `a·x² + b·x + c` (any of the four parse shapes), with
/// implicit coefficients made explicit.
fn quadratic_from_inner<'a>(
    ctx: &'a AtomArena<'a>,
    inner: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>, Atom<'a>)> {
    let (mut a, mut b, mut c) = (ctx.num(0), ctx.num(0), ctx.num(0));
    let mut found = false;
    match inner.node() {
        AtomNode::Add(args) => {
            for term in args.iter() {
                match term.node() {
                    AtomNode::Pow(base, exp) => {
                        if matches!(base.node(), AtomNode::Var(v) if *v == var)
                            && matches!(exp.node(), AtomNode::Num(2))
                        {
                            a = ctx.num(1); // implicit coefficient 1
                            found = true;
                        } else if crate::integral::is_constant(*term, var) {
                            c = *term;
                            found = true;
                        } else {
                            return None;
                        }
                    }
                    AtomNode::Var(v) if *v == var => {
                        b = ctx.num(1); // implicit coefficient 1
                        found = true;
                    }
                    AtomNode::Mul(args) if args.len() == 2 => {
                        let (coeff, rest) = (args[0], args[1]);
                        match rest.node() {
                            AtomNode::Pow(base, exp)
                                if matches!(base.node(), AtomNode::Var(v) if *v == var)
                                    && matches!(exp.node(), AtomNode::Num(2)) =>
                            {
                                if !crate::integral::is_constant(coeff, var) {
                                    return None;
                                }
                                a = coeff;
                                found = true;
                            }
                            AtomNode::Var(v) if *v == var => {
                                if !crate::integral::is_constant(coeff, var) {
                                    return None;
                                }
                                b = coeff;
                                found = true;
                            }
                            _ => {
                                if crate::integral::is_constant(*term, var) {
                                    c = *term;
                                    found = true;
                                } else {
                                    return None;
                                }
                            }
                        }
                    }
                    AtomNode::Num(_) => {
                        c = *term;
                        found = true;
                    }
                    _ => {
                        if crate::integral::is_constant(*term, var) {
                            c = *term;
                            found = true;
                        } else {
                            return None;
                        }
                    }
                }
            }
        }
        AtomNode::Pow(base, exp) => {
            if matches!(base.node(), AtomNode::Var(v) if *v == var)
                && matches!(exp.node(), AtomNode::Num(2))
            {
                a = ctx.num(1);
                found = true;
            }
        }
        AtomNode::Var(v) if *v == var => {
            b = ctx.num(1);
            found = true;
        }
        _ => {
            if crate::integral::is_constant(inner, var) {
                c = inner;
                found = true;
            }
        }
    }
    if !found {
        return None;
    }
    Some((a, b, c))
}

/// Evaluate a rational atom as `(p, q)` with `q > 0` (integers).
pub(crate) fn rat_of(atom: Atom<'_>) -> Option<(i64, i64)> {
    match atom.node() {
        AtomNode::Num(n) => Some((*n, 1)),
        _ => crate::integral::fraction_exponent(atom),
    }
}

/// `p/q` (in lowest terms) is the square of a rational `r/s`; returns
/// `(r, s)` with `s > 0` when it is, else `None`.
pub(crate) fn rational_sqrt(p: i64, q: i64) -> Option<(i64, i64)> {
    if q <= 0 || p < 0 {
        return None;
    }
    let g = gcd_i64(p, q);
    let (p, q) = (p / g, q / g);
    let r = isqrt(p)?;
    let s = isqrt(q)?;
    Some((r, s))
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.abs().max(1)
}

/// Integer square root if `n` is a perfect square.
fn isqrt(n: i64) -> Option<i64> {
    if n < 0 {
        return None;
    }
    let r = (n as f64).sqrt() as i64;
    [r - 1, r, r + 1]
        .into_iter()
        .find(|&c| c >= 0 && c * c == n)
}

/// Predicate helpers over match bindings (free-parameter constraints).
///
/// Bound atoms are read by name; `var` is the integration variable.
mod pred {
    use super::*;

    /// The atom bound to `name`, if it was bound as a single atom.
    pub(super) fn bound<'a>(bindings: &Bindings<'a>, name: &str) -> Option<Atom<'a>> {
        match bindings.get(Symbol::new(name))? {
            MatchValue::Single(a) => Some(*a),
            MatchValue::Sequence(_) => None,
        }
    }

    /// `free_q`: the atom does not contain `var`.
    pub(super) fn free_q(bindings: &Bindings<'_>, var: Symbol) -> bool {
        ["a", "b", "c", "d", "m", "n", "p"]
            .iter()
            .all(|n| match bound(bindings, n) {
                Some(a) => crate::integral::is_constant(a, var),
                None => true,
            })
    }

    /// The integer value bound to `name` (only valid for integer bindings).
    pub(super) fn int_val(bindings: &Bindings<'_>, name: &str) -> Option<i64> {
        match bound(bindings, name) {
            Some(a) => match a.node() {
                AtomNode::Num(n) => Some(*n),
                _ => None,
            },
            None => None,
        }
    }

    /// `pos_int_q`: integer literal ≥ 0.
    pub(super) fn pos_int_q(bindings: &Bindings<'_>, name: &str) -> bool {
        int_val(bindings, name).is_some_and(|n| n >= 0)
    }

    /// Integer literal ≥ 2.
    pub(super) fn int_ge_2(bindings: &Bindings<'_>, name: &str) -> bool {
        int_val(bindings, name).is_some_and(|n| n >= 2)
    }

    /// Integer literal ≠ -1; symbolic (non-integer) values pass.
    pub(super) fn not_minus_one(bindings: &Bindings<'_>, name: &str) -> bool {
        int_val(bindings, name).is_none_or(|n| n != -1)
    }

    /// The atom is not the literal 0.
    pub(super) fn nonzero(bindings: &Bindings<'_>, name: &str) -> bool {
        !matches!(bound(bindings, name), Some(a) if matches!(a.node(), AtomNode::Num(0)))
    }

    /// `n` is an integer in `[lo, hi]`.
    pub(super) fn int_in_range(bindings: &Bindings<'_>, name: &str, lo: i64, hi: i64) -> bool {
        int_val(bindings, name).is_some_and(|n| (lo..=hi).contains(&n))
    }

    /// `n` is a non-negative odd integer ≤ 9.
    pub(super) fn odd_le_9(bindings: &Bindings<'_>, name: &str) -> bool {
        int_val(bindings, name).is_some_and(|n| (1..=9).contains(&n) && n % 2 == 1)
    }

    /// Product-to-sum guard: not both `a` and `b` bound to the same literal.
    pub(super) fn a_neq_b(bindings: &Bindings<'_>, _var: Symbol) -> bool {
        match (int_val(bindings, "a"), int_val(bindings, "b")) {
            (Some(na), Some(nb)) => na != nb,
            _ => true,
        }
    }
}

/// The integration-rule library (0.27.0 rule families A–H).
///
/// Patterns use the literal variable `x` (baked per-call); free parameters
/// are `a_ b_ c_ d_ m_ n_ p_` constrained by predicates. Templates may
/// contain residual `Integral(g, x)` reduction terms.
fn rule_specs() -> &'static [RuleSpec] {
    use pred::*;
    &[
        // ------------------------------------------------------------------
        // A. Powers / binomials
        // ------------------------------------------------------------------
        // A1: x^n → x^(n+1)/(n+1), n a free integer ≠ -1 (or symbolic n).
        RuleSpec::Template {
            pat: "x^n_",
            tmpl: "x^(n_+1)*(n_+1)^-1",
            cond: Some(|b, var| free_q(b, var) && not_minus_one(b, "n")),
        },
        // A2: 1/(a+b*x) → log(a+b*x)/b (specific form, before A3).
        RuleSpec::Template {
            pat: "(a_+b_*x)^-1",
            tmpl: "log(a_+b_*x)*b_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        // A3: (a+b*x)^n → (a+b*x)^(n+1)/(b*(n+1)).
        RuleSpec::Template {
            pat: "(a_+b_*x)^n_",
            tmpl: "(a_+b_*x)^(n_+1)*(b_*(n_+1))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b") && not_minus_one(b, "n")),
        },
        // A4: x^m*(a+b*x)^n with non-negative integer m,n ≤ 8: expand and
        // integrate termwise. `c___` absorbs any constant coefficient.
        RuleSpec::Closure {
            pat: "c___*x^m_*(a_+b_*x)^n_",
            f: binomial_integrate,
            cond: Some(|b, _var| {
                int_in_range(b, "m", 0, 8)
                    && int_in_range(b, "n", 0, 8)
                    && bound(b, "a").is_some()
                    && bound(b, "b").is_some()
            }),
        },
        // ------------------------------------------------------------------
        // B. Exponentials / logarithms
        // ------------------------------------------------------------------
        // B1: exp(a*x+b) → exp(a*x+b)/a.
        RuleSpec::Template {
            pat: "exp(a_*x+b_)",
            tmpl: "exp(a_*x+b_)*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // B2: x^n*exp(a*x) reduction, 0 ≤ n ≤ 12 integer.
        RuleSpec::Template {
            pat: "x^n_*exp(a_*x)",
            tmpl: "x^n_*exp(a_*x)*a_^-1 + (-1)*n_*a_^-1*Integral(x^(n_ - 1)*exp(a_*x), x)",
            cond: Some(|b, var| free_q(b, var) && int_in_range(b, "n", 0, 12) && nonzero(b, "a")),
        },
        // B2b: same with a phase b in the exponent.
        RuleSpec::Template {
            pat: "x^n_*exp(a_*x+b_)",
            tmpl: "x^n_*exp(a_*x+b_)*a_^-1 + (-1)*n_*a_^-1*Integral(x^(n_ - 1)*exp(a_*x+b_), x)",
            cond: Some(|b, var| free_q(b, var) && int_in_range(b, "n", 0, 12) && nonzero(b, "a")),
        },
        // B3: exp(a*x)*sin(b*x) → e^(ax)(a sin(bx) − b cos(bx))/(a²+b²).
        RuleSpec::Template {
            pat: "exp(a_*x)*sin(b_*x)",
            tmpl: "exp(a_*x)*(a_*sin(b_*x) + (-1)*b_*cos(b_*x))*((a_^2+b_^2))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // B4: exp(a*x)*cos(b*x) → e^(ax)(a cos(bx) + b sin(bx))/(a²+b²).
        RuleSpec::Template {
            pat: "exp(a_*x)*cos(b_*x)",
            tmpl: "exp(a_*x)*(a_*cos(b_*x) + b_*sin(b_*x))*((a_^2+b_^2))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // B3b–B4b: bare-x variants (`exp(x)`/`sin(x)`/`cos(x)` args are not
        // `a*x` products, so the general patterns above do not match them).
        RuleSpec::Template {
            pat: "exp(x)*sin(b_*x)",
            tmpl: "exp(x)*(sin(b_*x) + (-1)*b_*cos(b_*x))*((1+b_^2))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        RuleSpec::Template {
            pat: "exp(x)*cos(b_*x)",
            tmpl: "exp(x)*(cos(b_*x) + b_*sin(b_*x))*((1+b_^2))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        RuleSpec::Template {
            pat: "exp(a_*x)*sin(x)",
            tmpl: "exp(a_*x)*(a_*sin(x) + (-1)*cos(x))*((a_^2+1))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "exp(a_*x)*cos(x)",
            tmpl: "exp(a_*x)*(a_*cos(x) + sin(x))*((a_^2+1))^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // B5: log(x)^n reduction, n ≥ 1.
        RuleSpec::Template {
            pat: "log(x)^n_",
            tmpl: "x*log(x)^n_ + (-1)*n_*Integral(log(x)^(n_ - 1), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // B6: x^m*log(x)^n reduction, m ≠ -1.
        RuleSpec::Template {
            pat: "x^m_*log(x)^n_",
            tmpl: "x^(m_+1)*log(x)^n_*(m_+1)^-1 + (-1)*n_*(m_+1)^-1*Integral(x^m_*log(x)^(n_ - 1), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n") && not_minus_one(b, "m")),
        },
        // B7: 1/(x*log(x)) → log(log(x)).
        RuleSpec::Template {
            pat: "x^-1*log(x)^-1",
            tmpl: "log(log(x))",
            cond: Some(free_q),
        },
        // B8: x*exp(a*x^2) → exp(a*x^2)/(2a) (substitution form).
        RuleSpec::Template {
            pat: "x*exp(a_*x^2)",
            tmpl: "exp(a_*x^2)*(2*a_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // ------------------------------------------------------------------
        // C. Trigonometry
        // ------------------------------------------------------------------
        // C1: sin(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "sin(x)^n_",
            tmpl: "(-1)*sin(x)^(n_ - 1)*cos(x)*n_^-1 + (n_ - 1)*n_^-1*Integral(sin(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C2: cos(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "cos(x)^n_",
            tmpl: "cos(x)^(n_ - 1)*sin(x)*n_^-1 + (n_ - 1)*n_^-1*Integral(cos(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C3: tan(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "tan(x)^n_",
            tmpl: "tan(x)^(n_ - 1)*(n_ - 1)^-1 + (-1)*Integral(tan(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C3b: tan(x) → −log(cos(x)); cot(x) → log(sin(x)) (bases for the
        // n-reductions' base cases).
        RuleSpec::Template {
            pat: "tan(x)",
            tmpl: "(-1)*log(cos(x))",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "cot(x)",
            tmpl: "log(sin(x))",
            cond: Some(free_q),
        },
        // C4: sec(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "sec(x)^n_",
            tmpl: "sec(x)^(n_ - 2)*tan(x)*(n_ - 1)^-1 + (n_ - 2)*(n_ - 1)^-1*Integral(sec(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C4b: sec(x) → log(sec(x)+tan(x)); csc(x) → log(tan(x/2)).
        RuleSpec::Template {
            pat: "sec(x)",
            tmpl: "log(sec(x)+tan(x))",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "csc(x)",
            tmpl: "log(tan(x*2^-1))",
            cond: Some(free_q),
        },
        // C13: linear-argument bases (the corpus is dominated by
        // `sin(c+d*x)`-style shapes; the bare-x rules above do not match).
        RuleSpec::Template {
            pat: "tan(a_*x+b_)",
            tmpl: "(-1)*log(cos(a_*x+b_))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "cot(a_*x+b_)",
            tmpl: "log(sin(a_*x+b_))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "sec(a_*x+b_)",
            tmpl: "log(sec(a_*x+b_)+tan(a_*x+b_))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "csc(a_*x+b_)",
            tmpl: "log(tan((a_*x+b_)*2^-1))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // C14: linear-argument power reductions. ∫ sin(u)^n dx with
        // u = a·x+b is the du-formula divided by a per level.
        RuleSpec::Template {
            pat: "sin(a_*x+b_)^n_",
            tmpl: "(-1)*sin(a_*x+b_)^(n_ - 1)*cos(a_*x+b_)*(n_*a_)^-1 + (n_ - 1)*(n_*a_)^-1*Integral(sin(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "cos(a_*x+b_)^n_",
            tmpl: "cos(a_*x+b_)^(n_ - 1)*sin(a_*x+b_)*(n_*a_)^-1 + (n_ - 1)*(n_*a_)^-1*Integral(cos(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "tan(a_*x+b_)^n_",
            tmpl: "tan(a_*x+b_)^(n_ - 1)*((n_ - 1)*a_)^-1 + (-1)*a_^-1*Integral(tan(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "cot(a_*x+b_)^n_",
            tmpl: "(-1)*cot(a_*x+b_)^(n_ - 1)*((n_ - 1)*a_)^-1 + (-1)*a_^-1*Integral(cot(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "sec(a_*x+b_)^n_",
            tmpl: "sec(a_*x+b_)^(n_ - 2)*tan(a_*x+b_)*((n_ - 1)*a_)^-1 + (n_ - 2)*((n_ - 1)*a_)^-1*Integral(sec(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "csc(a_*x+b_)^n_",
            tmpl: "(-1)*csc(a_*x+b_)^(n_ - 2)*cot(a_*x+b_)*((n_ - 1)*a_)^-1 + (n_ - 2)*((n_ - 1)*a_)^-1*Integral(csc(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        // C15: two-linear-argument product-to-sum (the dominant corpus
        // shape); subsumes one-side-bare variants.
        RuleSpec::Template {
            pat: "sin(a_*x+b_)*sin(c_*x+d_)",
            tmpl: "sin((a_+(-1)*c_)*x+(b_+(-1)*d_))*((a_+(-1)*c_)*2)^-1 + (-1)*sin((a_+c_)*x+(b_+d_))*((a_+c_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        RuleSpec::Template {
            pat: "cos(a_*x+b_)*cos(c_*x+d_)",
            tmpl: "sin((a_+(-1)*c_)*x+(b_+(-1)*d_))*((a_+(-1)*c_)*2)^-1 + sin((a_+c_)*x+(b_+d_))*((a_+c_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        RuleSpec::Template {
            pat: "sin(a_*x+b_)*cos(c_*x+d_)",
            tmpl: "(-1)*cos((a_+c_)*x+(b_+d_))*((a_+c_)*2)^-1 + (-1)*cos((a_+(-1)*c_)*x+(b_+(-1)*d_))*((a_+(-1)*c_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        // C5: cot(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "cot(x)^n_",
            tmpl: "(-1)*cot(x)^(n_ - 1)*(n_ - 1)^-1 + (-1)*Integral(cot(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C6: csc(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "csc(x)^n_",
            tmpl: "(-1)*csc(x)^(n_ - 2)*cot(x)*(n_ - 1)^-1 + (n_ - 2)*(n_ - 1)^-1*Integral(csc(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // C7: sin(a*x)*sin(b*x) product-to-sum.
        RuleSpec::Template {
            pat: "sin(a_*x)*sin(b_*x)",
            tmpl: "sin((a_+(-1)*b_)*x)*((a_+(-1)*b_)*2)^-1 + (-1)*sin((a_+b_)*x)*((a_+b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        // C8: cos(a*x)*cos(b*x) product-to-sum.
        RuleSpec::Template {
            pat: "cos(a_*x)*cos(b_*x)",
            tmpl: "sin((a_+(-1)*b_)*x)*((a_+(-1)*b_)*2)^-1 + sin((a_+b_)*x)*((a_+b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        // C9: sin(a*x)*cos(b*x) product-to-sum.
        RuleSpec::Template {
            pat: "sin(a_*x)*cos(b_*x)",
            tmpl: "(-1)*cos((a_+b_)*x)*((a_+b_)*2)^-1 + (-1)*cos((a_+(-1)*b_)*x)*((a_+(-1)*b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && a_neq_b(b, var)),
        },
        // C7b–C9b: bare-x product-to-sum variants (`sin(x)`/`cos(x)` args
        // are not `a*x` products, so the general patterns do not match).
        RuleSpec::Template {
            pat: "sin(x)*sin(b_*x)",
            tmpl: "sin((1+(-1)*b_)*x)*((1+(-1)*b_)*2)^-1 + (-1)*sin((1+b_)*x)*((1+b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "b") != Some(1)),
        },
        RuleSpec::Template {
            pat: "sin(a_*x)*sin(x)",
            tmpl: "sin((a_+(-1)*1)*x)*((a_+(-1)*1)*2)^-1 + (-1)*sin((a_+1)*x)*((a_+1)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "a") != Some(1)),
        },
        RuleSpec::Template {
            pat: "cos(x)*cos(b_*x)",
            tmpl: "sin((1+(-1)*b_)*x)*((1+(-1)*b_)*2)^-1 + sin((1+b_)*x)*((1+b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "b").is_none_or(|v| v != 1)),
        },
        RuleSpec::Template {
            pat: "cos(a_*x)*cos(x)",
            tmpl: "sin((a_+(-1)*1)*x)*((a_+(-1)*1)*2)^-1 + sin((a_+1)*x)*((a_+1)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "a").is_none_or(|v| v != 1)),
        },
        RuleSpec::Template {
            pat: "sin(x)*cos(b_*x)",
            tmpl: "(-1)*cos((1+b_)*x)*((1+b_)*2)^-1 + (-1)*cos((1+(-1)*b_)*x)*((1+(-1)*b_)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "b").is_none_or(|v| v != 1)),
        },
        RuleSpec::Template {
            pat: "sin(a_*x)*cos(x)",
            tmpl: "(-1)*cos((a_+1)*x)*((a_+1)*2)^-1 + (-1)*cos((a_+(-1)*1)*x)*((a_+(-1)*1)*2)^-1",
            cond: Some(|b, var| free_q(b, var) && int_val(b, "a").is_none_or(|v| v != 1)),
        },
        // C10: sin(x)*cos(x)^n → −cos(x)^(n+1)/(n+1), n ≥ 0.
        RuleSpec::Template {
            pat: "sin(x)*cos(x)^n_",
            tmpl: "(-1)*cos(x)^(n_+1)*(n_+1)^-1",
            cond: Some(|b, var| free_q(b, var) && pos_int_q(b, "n")),
        },
        // C11: sin(x)^m*cos(x) → sin(x)^(m+1)/(m+1), m ≥ 0.
        RuleSpec::Template {
            pat: "sin(x)^m_*cos(x)",
            tmpl: "sin(x)^(m_+1)*(m_+1)^-1",
            cond: Some(|b, var| free_q(b, var) && pos_int_q(b, "m")),
        },
        // C12: sin^m*cos^n with an odd exponent: peel the odd power and
        // expand the even remainder (closure).
        RuleSpec::Closure {
            pat: "sin(x)^m_*cos(x)^n_",
            f: trig_odd_power,
            cond: Some(|b, _var| {
                (odd_le_9(b, "m") && int_in_range(b, "n", 2, 9))
                    || (int_in_range(b, "m", 2, 9) && odd_le_9(b, "n"))
            }),
        },
        // C12b: same-argument linear trig products (u = a·x+b): templates
        // for the m=1/n=1 edges, closure for the odd-power peel.
        RuleSpec::Template {
            pat: "sin(a_*x+b_)*cos(a_*x+b_)^n_",
            tmpl: "(-1)*cos(a_*x+b_)^(n_ + 1)*((n_ + 1)*a_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && pos_int_q(b, "n")),
        },
        RuleSpec::Template {
            pat: "sin(a_*x+b_)^m_*cos(a_*x+b_)",
            tmpl: "sin(a_*x+b_)^(m_ + 1)*((m_ + 1)*a_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && pos_int_q(b, "m")),
        },
        RuleSpec::Closure {
            pat: "sin(a_*x+b_)^m_*cos(a_*x+b_)^n_",
            f: trig_odd_power_linear,
            cond: Some(|b, _var| {
                (odd_le_9(b, "m") && int_in_range(b, "n", 2, 9))
                    || (int_in_range(b, "m", 2, 9) && odd_le_9(b, "n"))
            }),
        },
        // ------------------------------------------------------------------
        // D. Hyperbolic functions
        // ------------------------------------------------------------------
        // D1: sinh(x)^n, n ≥ 2.
        RuleSpec::Template {
            pat: "sinh(x)^n_",
            tmpl: "sinh(x)^(n_ - 1)*cosh(x)*n_^-1 + (-1)*(n_ - 1)*n_^-1*Integral(sinh(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D2: cosh(x)^n, n ≥ 2.
        RuleSpec::Template {
            pat: "cosh(x)^n_",
            tmpl: "cosh(x)^(n_ - 1)*sinh(x)*n_^-1 + (n_ - 1)*n_^-1*Integral(cosh(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D3: tanh(x)^n reduction, n ≥ 2 (∫tanh^n = −tanh^(n−1)/(n−1) + ∫tanh^(n−2)).
        RuleSpec::Template {
            pat: "tanh(x)^n_",
            tmpl: "(-1)*tanh(x)^(n_ - 1)*(n_ - 1)^-1 + Integral(tanh(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D4: coth(x)^n reduction, n ≥ 2.
        RuleSpec::Template {
            pat: "coth(x)^n_",
            tmpl: "(-1)*coth(x)^(n_ - 1)*(n_ - 1)^-1 + Integral(coth(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D5: sech(x)^n, n ≥ 2.
        RuleSpec::Template {
            pat: "sech(x)^n_",
            tmpl: "sech(x)^(n_ - 2)*tanh(x)*(n_ - 1)^-1 + (n_ - 2)*(n_ - 1)^-1*Integral(sech(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D6: csch(x)^n, n ≥ 2.
        RuleSpec::Template {
            pat: "csch(x)^n_",
            tmpl: "(-1)*csch(x)^(n_ - 2)*coth(x)*(n_ - 1)^-1 + (-1)*(n_ - 2)*(n_ - 1)^-1*Integral(csch(x)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && int_ge_2(b, "n")),
        },
        // D1b: sinh(x) → cosh(x); cosh(x) → sinh(x) (base cases).
        RuleSpec::Template {
            pat: "sinh(x)",
            tmpl: "cosh(x)",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "cosh(x)",
            tmpl: "sinh(x)",
            cond: Some(free_q),
        },
        // D3b: tanh(x) → log(cosh(x)); coth(x) → log(sinh(x)) (base cases).
        RuleSpec::Template {
            pat: "tanh(x)",
            tmpl: "log(cosh(x))",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "coth(x)",
            tmpl: "log(sinh(x))",
            cond: Some(free_q),
        },
        // D5b: sech(x) → atan(sinh(x)); csch(x) → log(tanh(x/2)).
        RuleSpec::Template {
            pat: "sech(x)",
            tmpl: "atan(sinh(x))",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "csch(x)",
            tmpl: "log(tanh(x*2^-1))",
            cond: Some(free_q),
        },
        // D7b: linear-argument hyperbolic power reductions.
        RuleSpec::Template {
            pat: "sinh(a_*x+b_)^n_",
            tmpl: "sinh(a_*x+b_)^(n_ - 1)*cosh(a_*x+b_)*(n_*a_)^-1 + (-1)*(n_ - 1)*(n_*a_)^-1*Integral(sinh(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "cosh(a_*x+b_)^n_",
            tmpl: "cosh(a_*x+b_)^(n_ - 1)*sinh(a_*x+b_)*(n_*a_)^-1 + (n_ - 1)*(n_*a_)^-1*Integral(cosh(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "tanh(a_*x+b_)^n_",
            tmpl: "(-1)*tanh(a_*x+b_)^(n_ - 1)*((n_ - 1)*a_)^-1 + a_^-1*Integral(tanh(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        RuleSpec::Template {
            pat: "coth(a_*x+b_)^n_",
            tmpl: "(-1)*coth(a_*x+b_)^(n_ - 1)*((n_ - 1)*a_)^-1 + a_^-1*Integral(coth(a_*x+b_)^(n_ - 2), x)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a") && int_ge_2(b, "n")),
        },
        // D8: sinh(a*x+b) → cosh(a*x+b)/a (and mirrors).
        RuleSpec::Template {
            pat: "sinh(a_*x+b_)",
            tmpl: "cosh(a_*x+b_)*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "cosh(a_*x+b_)",
            tmpl: "sinh(a_*x+b_)*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "tanh(a_*x+b_)",
            tmpl: "log(cosh(a_*x+b_))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "coth(a_*x+b_)",
            tmpl: "log(sinh(a_*x+b_))*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // ------------------------------------------------------------------
        // E. Inverse trig / hyperbolic
        // ------------------------------------------------------------------
        RuleSpec::Template {
            pat: "asin(x)",
            tmpl: "x*asin(x) + sqrt(1-x^2)",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "acos(x)",
            tmpl: "x*acos(x) + (-1)*sqrt(1-x^2)",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "atan(x)",
            tmpl: "x*atan(x) + (-1)*log(x^2+1)*2^-1",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "asinh(x)",
            tmpl: "x*asinh(x) + (-1)*sqrt(x^2+1)",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "acosh(x)",
            tmpl: "x*acosh(x) + (-1)*sqrt(x^2 - 1)",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "atanh(x)",
            tmpl: "x*atanh(x) + log(1-x^2)*2^-1",
            cond: Some(free_q),
        },
        // E2: x·asin(x) and x·atan(x) (parts-based textbook forms).
        RuleSpec::Template {
            pat: "x*asin(x)",
            tmpl: "x^2*asin(x)*2^-1 + (-1)*asin(x)*4^-1 + x*(1-x^2)^(2^-1)*4^-1",
            cond: Some(free_q),
        },
        RuleSpec::Template {
            pat: "x*atan(x)",
            tmpl: "x^2*atan(x)*2^-1 + (-1)*x*2^-1 + atan(x)*2^-1",
            cond: Some(free_q),
        },
        // ------------------------------------------------------------------
        // F. Rational intercepts
        // ------------------------------------------------------------------
        // F1: 1/(a²+x²) → atan(x/a)/a.
        RuleSpec::Template {
            pat: "(a_^2+x^2)^-1",
            tmpl: "atan(x*a_^-1)*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // F2: 1/(a²−x²) → atanh(x/a)/a.
        RuleSpec::Template {
            pat: "(a_^2+(-1)*x^2)^-1",
            tmpl: "atanh(x*a_^-1)*a_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // F3: x/(a+b*x²) → log(a+b*x²)/(2b).
        RuleSpec::Template {
            pat: "x*(a_+b_*x^2)^-1",
            tmpl: "log(a_+b_*x^2)*(2*b_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        // ------------------------------------------------------------------
        // G. Radicals
        // ------------------------------------------------------------------
        // G1: sqrt(a+b*x) → 2(a+b*x)^(3/2)/(3b) (Fun and Pow forms).
        RuleSpec::Template {
            pat: "sqrt(a_+b_*x)",
            tmpl: "2*(a_+b_*x)^(3*2^-1)*(3*b_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        RuleSpec::Template {
            pat: "(a_+b_*x)^(2^-1)",
            tmpl: "2*(a_+b_*x)^(3*2^-1)*(3*b_)^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        // G2: 1/sqrt(a+b*x) → 2√(a+bx)/b (Fun and Pow forms).
        RuleSpec::Template {
            pat: "sqrt(a_+b_*x)^-1",
            tmpl: "2*sqrt(a_+b_*x)*b_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        RuleSpec::Template {
            pat: "(a_+b_*x)^(-2^-1)",
            tmpl: "2*(a_+b_*x)^(2^-1)*b_^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        // G3: x/sqrt(a+b*x) (closure, Fun and Pow forms).
        RuleSpec::Closure {
            pat: "x*sqrt(a_+b_*x)^-1",
            f: x_over_sqrt_linear,
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        RuleSpec::Closure {
            pat: "x*(a_+b_*x)^(-2^-1)",
            f: x_over_sqrt_linear,
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "b")),
        },
        // G4: sqrt(a²−x²) → x√(a²−x²)/2 + a² asin(x/a)/2.
        RuleSpec::Template {
            pat: "sqrt(a_^2+(-1)*x^2)",
            tmpl: "x*sqrt(a_^2+(-1)*x^2)*2^-1 + a_^2*asin(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "(a_^2+(-1)*x^2)^(2^-1)",
            tmpl: "x*(a_^2+(-1)*x^2)^(2^-1)*2^-1 + a_^2*asin(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // G5: 1/sqrt(a²−x²) → asin(x/a) (Pow and sqrt-Fun forms).
        RuleSpec::Template {
            pat: "(a_^2+(-1)*x^2)^(-2^-1)",
            tmpl: "asin(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "sqrt(a_^2+(-1)*x^2)^-1",
            tmpl: "asin(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // G6: 1/sqrt(x²+a²) → asinh(x/a) (Pow and sqrt-Fun forms).
        RuleSpec::Template {
            pat: "(x^2+a_^2)^(-2^-1)",
            tmpl: "asinh(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "sqrt(x^2+a_^2)^-1",
            tmpl: "asinh(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // G7: 1/sqrt(x²−a²) → acosh(x/a) (Pow and sqrt-Fun forms).
        RuleSpec::Template {
            pat: "(x^2+(-1)*a_^2)^(-2^-1)",
            tmpl: "acosh(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "sqrt(x^2+(-1)*a_^2)^-1",
            tmpl: "acosh(x*a_^-1)",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // G8: sqrt(x²+a²) → x√(x²+a²)/2 + a² asinh(x/a)/2.
        RuleSpec::Template {
            pat: "sqrt(x^2+a_^2)",
            tmpl: "x*sqrt(x^2+a_^2)*2^-1 + a_^2*asinh(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "(x^2+a_^2)^(2^-1)",
            tmpl: "x*(x^2+a_^2)^(2^-1)*2^-1 + a_^2*asinh(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // G9: sqrt(x²−a²) → x√(x²−a²)/2 − a² acosh(x/a)/2.
        RuleSpec::Template {
            pat: "sqrt(x^2+(-1)*a_^2)",
            tmpl: "x*sqrt(x^2+(-1)*a_^2)*2^-1 + (-1)*a_^2*acosh(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        RuleSpec::Template {
            pat: "(x^2+(-1)*a_^2)^(2^-1)",
            tmpl: "x*(x^2+(-1)*a_^2)^(2^-1)*2^-1 + (-1)*a_^2*acosh(x*a_^-1)*2^-1",
            cond: Some(|b, var| free_q(b, var) && nonzero(b, "a")),
        },
        // ------------------------------------------------------------------
        // H. Special-function extensions (shapes special.rs does not cover)
        // ------------------------------------------------------------------
        // H1: x*sin(x²) → −cos(x²)/2.
        RuleSpec::Template {
            pat: "x*sin(x^2)",
            tmpl: "(-1)*cos(x^2)*2^-1",
            cond: Some(free_q),
        },
        // H2: x*cos(x²) → sin(x²)/2.
        RuleSpec::Template {
            pat: "x*cos(x^2)",
            tmpl: "sin(x^2)*2^-1",
            cond: Some(free_q),
        },
        // H3: x*sinh(x²) → cosh(x²)/2.
        RuleSpec::Template {
            pat: "x*sinh(x^2)",
            tmpl: "cosh(x^2)*2^-1",
            cond: Some(free_q),
        },
        // H4: x*cosh(x²) → sinh(x²)/2.
        RuleSpec::Template {
            pat: "x*cosh(x^2)",
            tmpl: "sinh(x^2)*2^-1",
            cond: Some(free_q),
        },
    ]
}

/// Binomial expansion closure: ∫ x^m (a+bx)^n dx with non-negative integer
/// m, n ≤ 8, absorbing a constant coefficient bound to `c___`.
fn binomial_integrate<'a>(
    ctx: &'a AtomArena<'a>,
    bindings: &Bindings<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let m = pred::int_val(bindings, "m")?;
    let n = pred::int_val(bindings, "n")?;
    let a = pred::bound(bindings, "a")?;
    let b = pred::bound(bindings, "b")?;
    if !(0..=8).contains(&m) || !(0..=8).contains(&n) {
        return None;
    }
    // Constant coefficient: product of the `c___` sequence (1 when empty).
    let coeff: Atom<'a> = match bindings.get(Symbol::new("c")) {
        Some(MatchValue::Sequence(slice)) if !slice.is_empty() => ctx.mul(slice),
        _ => ctx.num(1),
    };
    let x = ctx.var(var.as_str());
    let mut terms: Vec<Atom<'a>> = Vec::new();
    for k in 0..=n {
        let binom = binomial_coeff(n, k);
        // C(n,k) * a^(n-k) * b^k * x^(m+k+1) / (m+k+1)
        let term = ctx.mul(&[
            ctx.num(binom),
            ctx.pow(a, ctx.num(n - k)),
            ctx.pow(b, ctx.num(k)),
            ctx.pow(x, ctx.num(m + k + 1)),
            ctx.pow(ctx.num(m + k + 1), ctx.num(-1)),
        ]);
        terms.push(term);
    }
    let sum = ctx.add(&terms);
    Some(if matches!(coeff.node(), AtomNode::Num(1)) {
        sum
    } else {
        ctx.mul(&[coeff, sum])
    })
}

/// Small binomial coefficient (n ≤ 8).
fn binomial_coeff(n: i64, k: i64) -> i64 {
    let mut r = 1i64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

/// sin^m cos^n with an odd exponent: peel the odd power, expand the even
/// remainder in terms of the other function, and integrate termwise.
fn trig_odd_power<'a>(
    ctx: &'a AtomArena<'a>,
    bindings: &Bindings<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let m = pred::int_val(bindings, "m")?;
    let n = pred::int_val(bindings, "n")?;
    let x = ctx.var(var.as_str());
    let sin = ctx.fun("sin", &[x]);
    let cos = ctx.fun("cos", &[x]);
    // Peel sin (m odd): ∫ sin^m cos^n dx = −∫ (1−cos²)^((m-1)/2) cos^n d(cos).
    if (1..=9).contains(&m) && m % 2 == 1 && (2..=9).contains(&n) {
        let h = (m - 1) / 2; // exponent of (1−cos²)
        let mut terms: Vec<Atom<'a>> = Vec::new();
        for k in 0..=h {
            // Term of −∫ C(h,k)(−1)^k cos^(n+2k) d(cos):
            // C(h,k)(−1)^(k+1) cos^(n+2k+1)/(n+2k+1)
            let sign = if k % 2 == 0 { -1i64 } else { 1i64 };
            let binom = binomial_coeff(h, k);
            let denom = n + 2 * k + 1;
            terms.push(ctx.mul(&[
                ctx.num(sign * binom),
                ctx.pow(cos, ctx.num(denom)),
                ctx.pow(ctx.num(denom), ctx.num(-1)),
            ]));
        }
        return Some(ctx.add(&terms));
    }
    // Peel cos (n odd): ∫ sin^m cos^n dx = ∫ (1−sin²)^((n-1)/2) sin^m d(sin).
    if (2..=9).contains(&m) && (1..=9).contains(&n) && n % 2 == 1 {
        let h = (n - 1) / 2;
        let mut terms: Vec<Atom<'a>> = Vec::new();
        for k in 0..=h {
            // ∫ C(h,k)(−1)^k sin^(m+2k) d(sin):
            // C(h,k)(−1)^k sin^(m+2k+1)/(m+2k+1)
            let sign = if k % 2 == 0 { 1i64 } else { -1i64 };
            let binom = binomial_coeff(h, k);
            let denom = m + 2 * k + 1;
            terms.push(ctx.mul(&[
                ctx.num(sign * binom),
                ctx.pow(sin, ctx.num(denom)),
                ctx.pow(ctx.num(denom), ctx.num(-1)),
            ]));
        }
        return Some(ctx.add(&terms));
    }
    None
}

/// sin^m cos^n with the same linear argument u = a·x+b and one odd
/// exponent: peel the odd power, expand the even remainder in terms of the
/// other function, and integrate termwise (du = a·dx, so every term carries
/// a 1/a).
fn trig_odd_power_linear<'a>(
    ctx: &'a AtomArena<'a>,
    bindings: &Bindings<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let m = pred::int_val(bindings, "m")?;
    let n = pred::int_val(bindings, "n")?;
    let a = pred::bound(bindings, "a")?;
    let b = pred::bound(bindings, "b")?;
    let x = ctx.var(var.as_str());
    let u = ctx.add(&[ctx.mul(&[a, x]), b]);
    let sin = ctx.fun("sin", &[u]);
    let cos = ctx.fun("cos", &[u]);
    let a_inv = ctx.pow(a, ctx.num(-1));
    // Peel sin (m odd): ∫ sin^m cos^n dx = −(1/a)∫(1−cos²)^h cos^n d(cos).
    if (1..=9).contains(&m) && m % 2 == 1 && (2..=9).contains(&n) {
        let h = (m - 1) / 2;
        let mut terms: Vec<Atom<'a>> = Vec::new();
        for k in 0..=h {
            let sign = if k % 2 == 0 { -1i64 } else { 1i64 };
            let binom = binomial_coeff(h, k);
            let denom = n + 2 * k + 1;
            terms.push(ctx.mul(&[
                ctx.num(sign * binom),
                a_inv,
                ctx.pow(cos, ctx.num(denom)),
                ctx.pow(ctx.num(denom), ctx.num(-1)),
            ]));
        }
        return Some(ctx.add(&terms));
    }
    // Peel cos (n odd): ∫ sin^m cos^n dx = (1/a)∫(1−sin²)^h sin^m d(sin).
    if (2..=9).contains(&m) && (1..=9).contains(&n) && n % 2 == 1 {
        let h = (n - 1) / 2;
        let mut terms: Vec<Atom<'a>> = Vec::new();
        for k in 0..=h {
            let sign = if k % 2 == 0 { 1i64 } else { -1i64 };
            let binom = binomial_coeff(h, k);
            let denom = m + 2 * k + 1;
            terms.push(ctx.mul(&[
                ctx.num(sign * binom),
                a_inv,
                ctx.pow(sin, ctx.num(denom)),
                ctx.pow(ctx.num(denom), ctx.num(-1)),
            ]));
        }
        return Some(ctx.add(&terms));
    }
    None
}

/// ∫ x/√(a+bx) dx = 2(bx − 2a)√(a+bx)/(3b²) (substitution closure).
fn x_over_sqrt_linear<'a>(
    ctx: &'a AtomArena<'a>,
    bindings: &Bindings<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let a = pred::bound(bindings, "a")?;
    let b = pred::bound(bindings, "b")?;
    let x = ctx.var(var.as_str());
    // sqrt(a+b*x) — keep the Fun form for printing consistency.
    let sqrt = ctx.fun("sqrt", &[ctx.add(&[a, ctx.mul(&[b, x])])]);
    let inner = ctx.add(&[ctx.mul(&[b, x]), ctx.mul(&[ctx.num(-2), a])]);
    let denom = ctx.mul(&[ctx.num(3), b, b]);
    Some(ctx.mul(&[ctx.num(2), inner, sqrt, ctx.pow(denom, ctx.num(-1))]))
}

/// True if `name` is a plain identifier `[A-Za-z][A-Za-z0-9]*`.
fn is_identifier(name: &str) -> bool {
    let mut chars = name.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric())
}

/// Replace whole-word `x` with `var` (word = `[A-Za-z0-9_]+` run), so
/// `exp(x)`/`x_`/`x` inside patterns bake the caller's variable while
/// wildcard names and function names like `exp` stay untouched.
fn bake_var(s: &str, var: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_' {
            let start = i;
            while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
                i += 1;
            }
            let word = &s[start..i];
            if word == "x" {
                out.push_str(var);
            } else {
                out.push_str(word);
            }
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out
}

/// Collect wildcard base names used in a pattern (for the variable-name
/// collision guard).
fn wildcard_names(pat: &Pattern<'_>, out: &mut Vec<String>) {
    match pat {
        Pattern::Wildcard(name, _) => out.push(name.as_str().to_string()),
        Pattern::Add(pats) | Pattern::Mul(pats) | Pattern::Fun(_, pats) => {
            for p in pats {
                wildcard_names(p, out);
            }
        }
        Pattern::Pow(p) => {
            wildcard_names(&p.0, out);
            wildcard_names(&p.1, out);
        }
        Pattern::Literal(_) => {}
    }
}

/// Build (and compile) the rule table for `var`.
///
/// Returns `None` when the variable name is not a plain identifier or
/// collides with a wildcard name used by the rules — the rules path is then
/// skipped entirely and the pipeline degrades to the 0.26 behaviour.
pub(crate) fn build_rule_table<'a>(
    ctx: &'a AtomArena<'a>,
    var: Symbol,
) -> Option<IntegralRuleTable<'a>> {
    let var_name = var.as_str();
    if !is_identifier(var_name) {
        return None;
    }
    let mut rules: Vec<Rule<'a>> = Vec::new();
    let mut closures: Vec<ClosureRule<'a>> = Vec::new();
    for spec in rule_specs() {
        let baked_pat = bake_var(spec_pat(spec), var_name);
        let parsed = ocas_parse::parse(ctx, &baked_pat).ok()?;
        let pattern = Pattern::from_atom(&crate::pattern_alloc::VecAlloc, parsed);
        // Variable-name collision guard: a free-parameter wildcard whose
        // base name equals the integration variable would make the rule
        // match the wrong things; skip the whole rules path.
        let mut names = Vec::new();
        wildcard_names(&pattern, &mut names);
        if names.iter().any(|n| n == var_name) {
            return None;
        }
        let head = pattern_head_key(&pattern);
        match spec {
            RuleSpec::Template { tmpl, cond, .. } => {
                let baked_tmpl = bake_var(tmpl, var_name);
                let mut rule = Rule::from_template(
                    ctx,
                    &crate::pattern_alloc::VecAlloc,
                    &baked_pat,
                    &baked_tmpl,
                );
                if let Some(pred) = cond {
                    rule = rule.with_condition(move |b: &Bindings<'a>| pred(b, var));
                }
                rules.push(rule);
            }
            RuleSpec::Closure { f, cond, .. } => {
                closures.push(ClosureRule {
                    head,
                    pattern,
                    cond: *cond,
                    f: *f,
                });
            }
        }
    }
    Some(IntegralRuleTable {
        table: RuleTable::from_rules(rules),
        closures,
    })
}

fn spec_pat(spec: &RuleSpec) -> &'static str {
    match spec {
        RuleSpec::Template { pat, .. } | RuleSpec::Closure { pat, .. } => pat,
    }
}

fn pattern_head_key(pattern: &Pattern<'_>) -> HeadKey {
    match pattern {
        Pattern::Fun(name, _) => HeadKey::Fun(*name),
        Pattern::Add(_) => HeadKey::Add,
        Pattern::Mul(_) => HeadKey::Mul,
        Pattern::Pow(_) => HeadKey::Pow,
        Pattern::Literal(_) | Pattern::Wildcard(_, _) => HeadKey::Any,
    }
}

impl<'a> IntegralRuleTable<'a> {
    fn apply(&self, ctx: &'a AtomArena<'a>, atom: Atom<'a>, var: Symbol) -> Option<Atom<'a>> {
        if let Some(r) = self.table.apply(ctx, atom) {
            return Some(r);
        }
        let key = head_of(atom);
        for rule in &self.closures {
            if rule.head != key && rule.head != HeadKey::Any {
                continue;
            }
            // A failed match on one closure must not abort the scan of the
            // remaining closures.
            let Ok(bindings) = match_pattern(rule.pattern.clone(), atom) else {
                continue;
            };
            if let Some(pred) = rule.cond
                && !pred(&bindings, var)
            {
                continue;
            }
            if let Some(r) = (rule.f)(ctx, &bindings, var) {
                return Some(r);
            }
        }
        None
    }
}

/// Fold `f(x)^1 -> f(x)` and `f(x)^0 -> 1` throughout an expression.
///
/// The rules engine sees expressions that have not been through the
/// simplify stage, so reduction templates instantiate to shapes like
/// `cot(x)^1` that no rule matches. Folding restores the canonical form
/// before matching and before resolving residuals.
fn fold_trivial_powers<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Pow(base, exp) => match exp.node() {
            AtomNode::Num(1) => fold_trivial_powers(ctx, *base),
            AtomNode::Num(0) => ctx.num(1),
            _ => {
                let b = fold_trivial_powers(ctx, *base);
                let e = fold_trivial_powers(ctx, *exp);
                ctx.pow(b, e)
            }
        },
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> =
                args.iter().map(|a| fold_trivial_powers(ctx, *a)).collect();
            let rebuilt = ctx.fun(name.as_str(), &rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> =
                args.iter().map(|a| fold_trivial_powers(ctx, *a)).collect();
            let rebuilt = ctx.add(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> =
                args.iter().map(|a| fold_trivial_powers(ctx, *a)).collect();
            let rebuilt = ctx.mul(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

/// Apply the rule table to `expr`, resolving residual `Integral(g, var)`
/// terms by recursive integration.
///
/// A residual that itself fails to integrate is a legitimate partial result
/// (reduction-formula semantics) and is returned as-is.
pub(crate) fn integrate_rules<'a>(
    ctx: &'a AtomArena<'a>,
    rules: &IntegralRuleTable<'a>,
    expr: Atom<'a>,
    var: Symbol,
    rule_depth: usize,
) -> Option<Atom<'a>> {
    if rule_depth >= MAX_RULE_DEPTH {
        return None;
    }
    let expr = fold_trivial_powers(ctx, expr);
    let applied = rules.apply(ctx, expr, var)?;
    let applied = fold_trivial_powers(ctx, applied);
    Some(resolve_integrals(ctx, applied, var, rule_depth))
}

/// Replace every `Integral(g, v)` with `v == var` by `integrate_raw(g, ...)`
/// at fresh structural depth and `rule_depth + 1` reduction budget.
fn resolve_integrals<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    rule_depth: usize,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Fun(name, args) if name.as_str() == "Integral" && args.len() == 2 => {
            let v = args[1];
            if matches!(v.node(), AtomNode::Var(s) if *s == var) {
                let g = fold_trivial_powers(ctx, args[0]);
                let resolved =
                    crate::integral::integrate_raw(ctx, g, var, 0, true, rule_depth + 1, 0);
                if resolved == expr {
                    // No progress (e.g. rule-depth budget exhausted): keep
                    // the residual rather than recursing forever.
                    return expr;
                }
                // The substitution result may itself carry nested
                // `Integral(.., var)` nodes (e.g. a heuristic partial
                // result); resolve those too.
                return resolve_integrals(ctx, resolved, var, rule_depth);
            }
            // A different integration variable: keep the residual.
            expr
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| resolve_integrals(ctx, *a, var, rule_depth))
                .collect();
            // Rebuild through the arena so hash-consing reuses the node when
            // nothing changed.
            let rebuilt = ctx.fun(name.as_str(), &rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| resolve_integrals(ctx, *a, var, rule_depth))
                .collect();
            let rebuilt = ctx.add(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| resolve_integrals(ctx, *a, var, rule_depth))
                .collect();
            let rebuilt = ctx.mul(&rebuilt);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Pow(base, exp) => {
            let b = resolve_integrals(ctx, *base, var, rule_depth);
            let e = resolve_integrals(ctx, *exp, var, rule_depth);
            let rebuilt = ctx.pow(b, e);
            if rebuilt == expr { expr } else { rebuilt }
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_core::arena::Arena;

    #[test]
    fn bake_var_replaces_whole_words_only() {
        assert_eq!(bake_var("sin(x)^n_ + exp(x)", "t"), "sin(t)^n_ + exp(t)");
        assert_eq!(bake_var("x_*exp(x)", "y"), "x_*exp(y)");
        assert_eq!(bake_var("exp(x)", "t"), "exp(t)");
        assert_eq!(bake_var("x", "t"), "t");
        // `x` inside a longer identifier is untouched.
        assert_eq!(bake_var("exp(x)+x_2", "t"), "exp(t)+x_2");
    }

    #[test]
    fn build_rule_table_guards_bad_variable_names() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // `1 x` is not an identifier (contains a space) — rules path skipped.
        assert!(build_rule_table(&ctx, Symbol::new("1 x")).is_none());
        // Plain identifiers are accepted.
        assert!(build_rule_table(&ctx, Symbol::new("x")).is_some());
        assert!(build_rule_table(&ctx, Symbol::new("t")).is_some());
        // Wildcard collision: `m` is a free parameter in the rules; a
        // variable named `m` must skip the rules path.
        assert!(build_rule_table(&ctx, Symbol::new("m")).is_none());
    }

    /// Integrate via the full public chain and return the result string.
    fn int_str(input: &str, var: &str) -> String {
        use ocas_core::arena::Arena;
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, input).unwrap();
        crate::integrate(&ctx, expr, Symbol::new(var)).to_string()
    }

    /// Assert the integral contains no unevaluated `Integral(...)` residue.
    fn assert_solved(input: &str) {
        let r = int_str(input, "x");
        assert!(
            !r.contains("Integral("),
            "integrate({input}) left a residue: {r}"
        );
    }

    #[test]
    fn family_a_powers() {
        assert_solved("x^n");
        assert_solved("(a+b*x)^n");
        assert_solved("x^2*(a+b*x)^3");
        assert_solved("2*x^3*(a+b*x)^2");
        assert_solved("x^5*(a+b*x)^4");
    }

    #[test]
    fn family_b_exp_log() {
        assert_solved("exp(3*x+1)");
        assert_solved("x^3*exp(2*x)");
        assert_solved("x*exp(-x^2)");
        assert_solved("exp(2*x)*sin(3*x)");
        assert_solved("log(x)^3");
        assert_solved("x^2*log(x)^2");
        assert_solved("1/(x*log(x))");
    }

    #[test]
    fn family_c_trig() {
        assert_solved("sin(x)^5");
        assert_solved("cos(x)^4");
        assert_solved("tan(x)^4");
        assert_solved("sec(x)^4");
        assert_solved("cot(x)^3");
        assert_solved("csc(x)^3");
        assert_solved("sin(x)^3*cos(x)^2");
        assert_solved("sin(x)^2*cos(x)^5");
        assert_solved("sin(x)*cos(x)^3");
        assert_solved("sin(x)^3*cos(x)");
        assert_solved("sin(2*x)*sin(3*x)");
        assert_solved("cos(2*x)*cos(3*x)");
        assert_solved("sin(2*x)*cos(3*x)");
        // Plan acceptance: tan(x)^4 → tan(x)^3/3 − tan(x) + x (equivalence
        // in oCAS print form: `(-1*((tan(x)) + (-1*x))) + ((3^-1)*((tan(x))^3))`).
        let r = int_str("tan(x)^4", "x");
        assert!(!r.contains("Integral("), "tan^4 left a residue: {r}");
        assert!(r.contains("tan(x)") && r.contains("x"), "tan^4 result: {r}");
        assert!(r.contains("^3"), "tan^4 result lacks the cubic term: {r}");
    }

    #[test]
    fn family_d_hyperbolic() {
        assert_solved("sinh(x)^4");
        assert_solved("cosh(x)^4");
        assert_solved("tanh(x)^3");
        assert_solved("coth(x)^3");
        assert_solved("sech(x)^3");
        assert_solved("csch(x)^3");
        assert_solved("sinh(2*x+1)");
        assert_solved("cosh(2*x+1)");
        assert_solved("tanh(2*x+1)");
        assert_solved("coth(2*x+1)");
    }

    #[test]
    fn family_e_inverse_trig() {
        assert_solved("asin(x)");
        assert_solved("acos(x)");
        assert_solved("atan(x)");
        assert_solved("asinh(x)");
        assert_solved("acosh(x)");
        assert_solved("atanh(x)");
    }

    #[test]
    fn family_f_rational_intercepts() {
        assert_solved("1/(a^2+x^2)");
        assert_solved("1/(a^2-x^2)");
        assert_solved("x/(a+b*x^2)");
    }

    #[test]
    fn family_g_radicals() {
        assert_solved("sqrt(a+b*x)");
        assert_solved("1/sqrt(a+b*x)");
        assert_solved("x/sqrt(a+b*x)");
        assert_solved("sqrt(a^2-x^2)");
        assert_solved("1/sqrt(a^2-x^2)");
        assert_solved("1/sqrt(x^2+a^2)");
        assert_solved("1/sqrt(x^2-a^2)");
        assert_solved("sqrt(x^2+a^2)");
        assert_solved("sqrt(x^2-a^2)");
    }

    #[test]
    fn family_h_special_forms() {
        assert_solved("x*sin(x^2)");
        assert_solved("x*cos(x^2)");
        assert_solved("x*sinh(x^2)");
        assert_solved("x*cosh(x^2)");
    }

    #[test]
    fn rules_off_returns_unevaluated() {
        use ocas_core::arena::Arena;
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "tan(x)^4").unwrap();
        let r = crate::integrate_with_options(
            &ctx,
            expr,
            Symbol::new("x"),
            crate::IntegrateOptions { rules: false },
        );
        assert!(r.to_string().contains("Integral("));
    }
}
