//! Symbolic differentiation for oCAS.
//!
//! This module provides [`diff`], which computes the symbolic derivative of an
//! expression with respect to a variable. Elementary function rules are
//! hard-coded in a small table, and the chain rule is applied automatically for
//! compound arguments.
//!
//! # Special-function heads (0.27.3)
//!
//! The table also covers the special functions the integration pipeline emits:
//!
//! - `erf`, `erfc`, `erfi` (with the `2/√π` prefactor),
//! - `Ei`, `Si`, `Ci`, `Shi`, `Chi`,
//! - `fresnels`, `fresnelc` (`d/dz S(z) = sin(π z²/2)`,
//!   `d/dz C(z) = cos(π z²/2)`).
//!
//! Three heads take several arguments and are differentiated per argument:
//!
//! - `Ei(n, z)` is the exponential integral `Eₙ(z)` (the Rubi/Mathematica
//!   convention, matching `sympy.expint(n, z)`), with `∂/∂z Eₙ(z) = −Eₙ₋₁(z)`
//!   and `E₀(z) = e⁻ᶻ/z` kept elementary,
//! - `EllipticF(φ, m)`, `EllipticE(φ, m)`, `EllipticPi(n, φ, m)` with the
//!   standard `∂/∂φ` partials.
//!
//! Only the partials the integration pipeline needs are implemented (the
//! order/modulus/characteristic slots are treated as symbolic parameters). If
//! an **unmapped** slot's argument actually depends on the differentiation
//! variable, the whole expression is returned as an unevaluated
//! `Derivative(f, x)` rather than silently dropping the term — this module
//! never guesses a partial derivative.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_rewrite::simplify::simplify;

use crate::rules::calculus_rules;

/// Differentiate `expr` with respect to `var`.
///
/// The function implements the standard recursive derivative rules for sums,
/// products, powers, and a built-in table of elementary functions. The result
/// is simplified using the default rewrite rules plus a few calculus-specific
/// identities.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::diff;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let sin_x = ctx.fun("sin", &[x]);
/// let result = diff(&ctx, sin_x, Symbol::new("x"));
/// assert_eq!(result.to_string(), "cos(x)");
/// ```
///
/// For unknown functions, the derivative is returned as an unevaluated
/// `Derivative(f, x)` form.
pub fn diff<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let raw = diff_raw(ctx, expr, var);
    let simplified = simplify(ctx, raw, &rules, 20);
    normalize(ctx, simplified)
}

fn diff_raw<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    match expr.node() {
        AtomNode::Num(_) => ctx.num(0),
        AtomNode::Var(v) => {
            if *v == var {
                ctx.num(1)
            } else {
                ctx.num(0)
            }
        }
        AtomNode::Add(args) => {
            let mut terms = Vec::with_capacity(args.len());
            for a in args.iter() {
                terms.push(diff_raw(ctx, *a, var));
            }
            ctx.add(&terms)
        }
        AtomNode::Mul(args) => {
            // Product rule: d/dx (a * b * c) = a' * b * c + a * b' * c + ...
            let mut sum_terms = Vec::with_capacity(args.len());
            for i in 0..args.len() {
                let mut factors = Vec::with_capacity(args.len());
                for (j, a) in args.iter().enumerate() {
                    if i == j {
                        factors.push(diff_raw(ctx, *a, var));
                    } else {
                        factors.push(*a);
                    }
                }
                sum_terms.push(ctx.mul(&factors));
            }
            ctx.add(&sum_terms)
        }
        AtomNode::Pow(base, exp) => {
            let base = *base;
            let exp = *exp;
            let d_base = diff_raw(ctx, base, var);
            let d_exp = diff_raw(ctx, exp, var);

            let is_exp_const = matches!(exp.node(), AtomNode::Num(_));
            let is_base_const = matches!(base.node(), AtomNode::Num(_));

            if is_exp_const {
                // d/dx (b^n) = n * b^(n-1) * b'
                let n = if let AtomNode::Num(n) = exp.node() {
                    *n
                } else {
                    0
                };
                ctx.mul(&[ctx.num(n), ctx.pow(base, ctx.num(n - 1)), d_base])
            } else if is_base_const {
                // d/dx (a^u) = a^u * log(a) * u'
                ctx.mul(&[ctx.pow(base, exp), ctx.fun("log", &[base]), d_exp])
            } else {
                // Generalized power rule:
                //   d/dx (b^e) = b^e * (log(b) * e' + e * b' / b)
                let log_b = ctx.fun("log", &[base]);
                let term1 = ctx.mul(&[log_b, d_exp]);
                let term2 = ctx.mul(&[exp, d_base, ctx.pow(base, ctx.num(-1))]);
                let factor = ctx.add(&[term1, term2]);
                ctx.mul(&[ctx.pow(base, exp), factor])
            }
        }
        AtomNode::Fun(name, args) => diff_function(ctx, *name, args, var),
    }
}

fn diff_function<'a>(
    ctx: &'a AtomArena<'a>,
    name: Symbol,
    args: &'a [Atom<'a>],
    var: Symbol,
) -> Atom<'a> {
    debug_assert!(
        !args.is_empty(),
        "diff_function should only be called with non-empty function arguments"
    );
    let name_str = name.as_str();
    // Multi-argument heads carry one partial per argument; handle them before
    // the single-argument table so `args[1]` is not silently ignored.
    if let Some(d) = multi_arg_partial(ctx, name_str, args, var) {
        return d;
    }
    let u = args[0];
    let du = diff_raw(ctx, u, var);

    let derivative_of_arg: Atom<'a> = match name_str {
        "sin" => ctx.fun("cos", &[u]),
        "cos" => ctx.mul(&[ctx.num(-1), ctx.fun("sin", &[u])]),
        "exp" => ctx.fun("exp", &[u]),
        "log" => ctx.pow(u, ctx.num(-1)),
        "sqrt" => ctx.pow(ctx.mul(&[ctx.num(2), ctx.fun("sqrt", &[u])]), ctx.num(-1)),
        "tan" => ctx.pow(ctx.fun("sec", &[u]), ctx.num(2)),
        "sec" => ctx.mul(&[ctx.fun("sec", &[u]), ctx.fun("tan", &[u])]),
        "cot" => ctx.mul(&[ctx.num(-1), ctx.pow(ctx.fun("csc", &[u]), ctx.num(2))]),
        "csc" => ctx.mul(&[ctx.num(-1), ctx.fun("csc", &[u]), ctx.fun("cot", &[u])]),
        "asin" => ctx.pow(
            ctx.fun(
                "sqrt",
                &[ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(u, ctx.num(2))])])],
            ),
            ctx.num(-1),
        ),
        "acos" => ctx.mul(&[
            ctx.num(-1),
            ctx.pow(
                ctx.fun(
                    "sqrt",
                    &[ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(u, ctx.num(2))])])],
                ),
                ctx.num(-1),
            ),
        ]),
        "atan" => ctx.pow(ctx.add(&[ctx.num(1), ctx.pow(u, ctx.num(2))]), ctx.num(-1)),
        "atanh" => ctx.pow(
            ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), ctx.pow(u, ctx.num(2))])]),
            ctx.num(-1),
        ),
        "sinh" => ctx.fun("cosh", &[u]),
        "cosh" => ctx.fun("sinh", &[u]),
        "tanh" => ctx.pow(ctx.fun("sech", &[u]), ctx.num(2)),
        "coth" => ctx.mul(&[ctx.num(-1), ctx.pow(ctx.fun("csch", &[u]), ctx.num(2))]),
        "sech" => ctx.mul(&[ctx.num(-1), ctx.fun("sech", &[u]), ctx.fun("tanh", &[u])]),
        "csch" => ctx.mul(&[ctx.num(-1), ctx.fun("csch", &[u]), ctx.fun("coth", &[u])]),
        "asinh" => ctx.pow(
            ctx.fun("sqrt", &[ctx.add(&[ctx.pow(u, ctx.num(2)), ctx.num(1)])]),
            ctx.num(-1),
        ),
        "acosh" => ctx.pow(
            ctx.fun("sqrt", &[ctx.add(&[ctx.pow(u, ctx.num(2)), ctx.num(-1)])]),
            ctx.num(-1),
        ),
        // Special-function heads (0.27.3): each entry is `f'(u)`; the chain
        // rule below multiplies by `u'`.
        "erf" => ctx.mul(&[
            ctx.num(2),
            inv_sqrt_pi(ctx),
            ctx.fun("exp", &[neg_square(ctx, u)]),
        ]),
        "erfc" => ctx.mul(&[
            ctx.num(-2),
            inv_sqrt_pi(ctx),
            ctx.fun("exp", &[neg_square(ctx, u)]),
        ]),
        "erfi" => ctx.mul(&[
            ctx.num(2),
            inv_sqrt_pi(ctx),
            ctx.fun("exp", &[ctx.pow(u, ctx.num(2))]),
        ]),
        // `Ei(u) = ∫ eᵗ/t dt`, so `d/du Ei(u) = eᵘ/u`. The two-argument form
        // `Ei(n, z)` is a different function (`Eₙ`) and is handled above.
        "Ei" => ctx.mul(&[ctx.fun("exp", &[u]), ctx.pow(u, ctx.num(-1))]),
        "Si" => ctx.mul(&[ctx.fun("sin", &[u]), ctx.pow(u, ctx.num(-1))]),
        "Ci" => ctx.mul(&[ctx.fun("cos", &[u]), ctx.pow(u, ctx.num(-1))]),
        "Shi" => ctx.mul(&[ctx.fun("sinh", &[u]), ctx.pow(u, ctx.num(-1))]),
        "Chi" => ctx.mul(&[ctx.fun("cosh", &[u]), ctx.pow(u, ctx.num(-1))]),
        // SymPy convention: `S(z) = ∫₀ᶻ sin(π t²/2) dt`.
        "fresnels" => ctx.fun("sin", &[pi_half_square(ctx, u)]),
        "fresnelc" => ctx.fun("cos", &[pi_half_square(ctx, u)]),
        _ => {
            // Unknown function: return an unevaluated Derivative form.
            return ctx.fun(
                "Derivative",
                &[ctx.fun(name_str, args), ctx.var(var.as_str())],
            );
        }
    };

    ctx.mul(&[derivative_of_arg, du])
}

/// `1/√π`, the prefactor shared by the `erf` family.
fn inv_sqrt_pi<'a>(ctx: &'a AtomArena<'a>) -> Atom<'a> {
    ctx.pow(ctx.fun("sqrt", &[ctx.var("pi")]), ctx.num(-1))
}

/// `-u²`.
fn neg_square<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>) -> Atom<'a> {
    ctx.mul(&[ctx.num(-1), ctx.pow(u, ctx.num(2))])
}

/// `π·u²/2`, the Fresnel argument.
fn pi_half_square<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>) -> Atom<'a> {
    ctx.mul(&[
        ctx.var("pi"),
        ctx.pow(u, ctx.num(2)),
        ctx.pow(ctx.num(2), ctx.num(-1)),
    ])
}

/// Whether a raw derivative is structurally zero (the argument does not
/// depend on the differentiation variable).
fn is_zero_derivative(d: Atom<'_>) -> bool {
    matches!(d.node(), AtomNode::Num(0))
}

/// The unevaluated `Derivative(f(args), var)` form.
fn unevaluated<'a>(
    ctx: &'a AtomArena<'a>,
    name: &str,
    args: &'a [Atom<'a>],
    var: Symbol,
) -> Atom<'a> {
    ctx.fun("Derivative", &[ctx.fun(name, args), ctx.var(var.as_str())])
}

/// Per-argument partials for the multi-argument special heads.
///
/// Returns `None` when `name`/arity is not one of them (the caller then uses
/// the single-argument table). Returns the unevaluated `Derivative` form when
/// an argument sits in a slot whose partial is not implemented **and** that
/// argument really depends on `var`: declining is the only honest answer, a
/// zero term would be a wrong answer.
fn multi_arg_partial<'a>(
    ctx: &'a AtomArena<'a>,
    name: &str,
    args: &'a [Atom<'a>],
    var: Symbol,
) -> Option<Atom<'a>> {
    // `Ei(n, z)` = `Eₙ(z)`; `∂/∂z = −Eₙ₋₁(z)`, `∂/∂n` is not implemented.
    if name == "Ei" && args.len() == 2 {
        let (n, z) = (args[0], args[1]);
        if !is_zero_derivative(diff_raw(ctx, n, var)) {
            return Some(unevaluated(ctx, name, args, var));
        }
        let dz = diff_raw(ctx, z, var);
        if is_zero_derivative(dz) {
            return Some(ctx.num(0));
        }
        let partial = match n.node() {
            // `E₁' = −E₀ = −e⁻ᶻ/z`: keep the elementary form instead of
            // emitting a degenerate zero-order head.
            AtomNode::Num(1) => ctx.mul(&[
                ctx.num(-1),
                ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), z])]),
                ctx.pow(z, ctx.num(-1)),
            ]),
            AtomNode::Num(k) => ctx.mul(&[ctx.num(-1), ctx.fun("Ei", &[ctx.num(k - 1), z])]),
            _ => ctx.mul(&[ctx.num(-1), ctx.fun("Ei", &[ctx.add(&[n, ctx.num(-1)]), z])]),
        };
        return Some(ctx.mul(&[partial, dz]));
    }
    // Elliptic integrals: `EllipticF(φ, m)`, `EllipticE(φ, m)`,
    // `EllipticPi(n, φ, m)`. Only `∂/∂φ` is implemented; the modulus and
    // characteristic slots are symbolic parameters.
    let phi_index = match name {
        "EllipticF" | "EllipticE" if args.len() == 2 => 0usize,
        "EllipticPi" if args.len() == 3 => 1usize,
        _ => return None,
    };
    let unmapped: &[usize] = match name {
        "EllipticPi" => &[0, 2],
        _ => &[1],
    };
    if unmapped
        .iter()
        .any(|&i| !is_zero_derivative(diff_raw(ctx, args[i], var)))
    {
        return Some(unevaluated(ctx, name, args, var));
    }
    let phi = args[phi_index];
    let dphi = diff_raw(ctx, phi, var);
    if is_zero_derivative(dphi) {
        return Some(ctx.num(0));
    }
    let sin2 = ctx.pow(ctx.fun("sin", &[phi]), ctx.num(2));
    let (m, n) = match name {
        "EllipticPi" => (args[2], Some(args[0])),
        _ => (args[1], None),
    };
    // `1 − m·sin²φ`; the definition of all three integrals is real exactly
    // where this is positive, and an out-of-domain sample is reported as such
    // by the numerical oracle rather than as a disagreement.
    let one_minus_m_sin2 = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), m, sin2])]);
    let root = ctx.fun("sqrt", &[one_minus_m_sin2]);
    let partial = match name {
        "EllipticF" => ctx.pow(root, ctx.num(-1)),
        "EllipticE" => root,
        _ => {
            let n = n.expect("EllipticPi has a characteristic argument");
            let one_minus_n_sin2 = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), n, sin2])]);
            ctx.pow(ctx.mul(&[one_minus_n_sin2, root]), ctx.num(-1))
        }
    };
    Some(ctx.mul(&[partial, dphi]))
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn diff_number() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let n = ctx.num(7);
        let result = diff(&ctx, n, Symbol::new("x"));
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn diff_variable_same() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let result = diff(&ctx, x, Symbol::new("x"));
        assert_eq!(result.to_string(), "1");
    }

    #[test]
    fn diff_variable_other() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let y = ctx.var("y");
        let result = diff(&ctx, y, Symbol::new("x"));
        assert_eq!(result.to_string(), "0");
    }

    #[test]
    fn diff_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(2));
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "2*x");
    }

    #[test]
    fn diff_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "cos(x)");
    }

    #[test]
    fn diff_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("cos", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "-1*(sin(x))");
    }

    #[test]
    fn diff_sqrt() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sqrt", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(2*(sqrt(x)))^-1");
    }

    #[test]
    fn diff_atan() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("atan", &[x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(1 + (x^2))^-1");
    }

    #[test]
    fn diff_exp_squared() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let x2 = ctx.pow(x, ctx.num(2));
        let expr = ctx.fun("exp", &[x2]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "2*x*(exp(x^2))");
    }

    #[test]
    fn diff_product() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let sin_x = ctx.fun("sin", &[x]);
        let expr = ctx.mul(&[x, sin_x]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(sin(x)) + (x*(cos(x)))");
    }

    #[test]
    fn diff_unknown_function() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let y = ctx.var("y");
        let expr = ctx.fun("f", &[x, y]);
        let result = diff(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Derivative(f(x, y), x)");
    }

    /// Render `diff(f(x), x)` for a parsed one-variable expression.
    fn d(expr: &str) -> String {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let e = ocas_parse::parse(&ctx, expr).expect("parse");
        diff(&ctx, e, Symbol::new("x")).to_string()
    }

    #[test]
    fn diff_special_function_heads() {
        // The canonical table entries as `(input, required substrings)`. The
        // rendering is normalized by `simplify`, so the assertions pin the
        // mathematical content (which factors must be present) rather than a
        // byte-exact form; the numerical oracle in `ocas-tests` is the
        // independent check that the content is right.
        let cases: &[(&str, &[&str])] = &[
            ("erf(x)", &["exp(-1*(x^2))", "sqrt(pi)"]),
            ("erfc(x)", &["-2", "exp(-1*(x^2))", "sqrt(pi)"]),
            ("erfi(x)", &["2", "exp(x^2)", "sqrt(pi)"]),
            ("Ei(x)", &["exp(x)", "x^-1"]),
            ("Si(x)", &["sin(x)", "x^-1"]),
            ("Ci(x)", &["cos(x)", "x^-1"]),
            ("Shi(x)", &["sinh(x)", "x^-1"]),
            ("Chi(x)", &["cosh(x)", "x^-1"]),
            ("fresnels(x)", &["sin(", "pi", "x^2"]),
            ("fresnelc(x)", &["cos(", "pi", "x^2"]),
        ];
        for (input, needles) in cases {
            let got = d(input);
            for needle in *needles {
                assert!(
                    got.contains(needle),
                    "diff({input}) = {got:?} does not contain {needle:?}"
                );
            }
            // The derivative must not be left unevaluated.
            assert!(
                !got.starts_with("Derivative("),
                "diff({input}) fell through to an unevaluated form: {got}"
            );
        }
    }

    #[test]
    fn diff_ei_order_argument() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∂/∂z E₃(z) = −E₂(z).
        let e3 = ctx.fun("Ei", &[ctx.num(3), x]);
        assert_eq!(
            diff(&ctx, e3, Symbol::new("x")).to_string(),
            "-1*(Ei(2, x))"
        );
        // ∂/∂z E₁(z) = −E₀(z) = −e⁻ᶻ/z, kept elementary.
        let e1 = ctx.fun("Ei", &[ctx.num(1), x]);
        let got = diff(&ctx, e1, Symbol::new("x")).to_string();
        assert!(
            got.contains("exp(-1*x)") && !got.contains("Ei(0"),
            "expected an elementary -e^-x/x form, got {got}"
        );
    }

    #[test]
    fn diff_ei_order_slot_declines_instead_of_guessing() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // The order slot depends on the differentiation variable and its
        // partial is not implemented: an unevaluated Derivative is the only
        // honest answer (a zero term would be a wrong answer).
        let e = ctx.fun("Ei", &[x, ctx.var("y")]);
        let got = diff(&ctx, e, Symbol::new("x")).to_string();
        assert!(got.starts_with("Derivative("), "got {got}");
    }

    #[test]
    fn diff_elliptic_integrals() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let m = ctx.var("m");
        let phi = ctx.var("x");
        let f = ctx.fun("EllipticF", &[phi, m]);
        let got = diff(&ctx, f, Symbol::new("x")).to_string();
        // ∂F/∂φ = 1/√(1 − m·sin²φ)
        assert!(
            got.contains("sqrt(1")
                && got.contains("m*")
                && got.contains("sin(x)")
                && got.ends_with("^-1"),
            "diff(EllipticF) = {got}"
        );
        // The modulus slot's partial is not implemented: differentiating with
        // respect to `m` must decline rather than drop the term.
        assert_eq!(
            diff(&ctx, f, Symbol::new("m")).to_string(),
            "Derivative(EllipticF(x, m), m)"
        );
        // An unrelated variable makes both slots independent, so the answer is
        // a genuine zero.
        assert_eq!(diff(&ctx, f, Symbol::new("y")).to_string(), "0");
        // ... and if the modulus really depends on `x`, decline honestly.
        let f_bad = ctx.fun("EllipticF", &[phi, x]);
        let got = diff(&ctx, f_bad, Symbol::new("x")).to_string();
        assert!(got.starts_with("Derivative("), "got {got}");
    }
}
