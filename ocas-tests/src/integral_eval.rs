//! Numerical verification oracle for symbolic integration results.
//!
//! The Rubi coverage harness (`benches/integrate_1892.rs`) counts a problem as
//! "solved" when the result carries no `Integral(...)` residue. That string
//! criterion alone cannot tell a correct antiderivative from a
//! plausible-looking wrong one — 0.27.1 shipped two wrong-answer classes that
//! the string criterion accepted. This module adds an independent numerical
//! criterion:
//!
//! 1. assign deterministic dummy values to every free parameter,
//! 2. evaluate the claimed antiderivative `F` at `x ± h` and `x ± 2h`,
//! 3. compare the 5-point central difference against the integrand `f` at `x`.
//!
//! A result that passes is *verified*; one that fails is a `Mismatch`; one
//! whose functions the oracle cannot evaluate (Ei, Fresnel, the imaginary
//! unit, domains it cannot sample) is `Indeterminate`. The harness reports the
//! three counts separately and never merges them.

use ocas_atom::{Atom, AtomNode, Symbol};

/// Outcome of evaluating an atom at a concrete environment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Eval {
    /// The expression evaluated to a finite real value.
    Value(f64),
    /// The point lies outside the function's real domain.
    Domain,
    /// The expression overflowed or divided by zero (a non-finite value):
    /// distinct from [`Eval::Domain`] because it usually means the *emitted
    /// form* is invalid rather than the sample being unusable.
    NonFinite,
    /// A function head the oracle deliberately does not implement.
    Unsupported,
}

/// Result of one numerical verification attempt.
#[derive(Debug, Clone)]
pub enum Verify {
    /// Every usable sample agreed with the integrand.
    Verified {
        /// Number of sample points that produced a comparison.
        checked: usize,
        /// Largest observed relative error.
        worst_rel: f64,
    },
    /// At least one sample disagreed: the claimed antiderivative is wrong.
    Mismatch {
        /// Number of sample points that produced a comparison.
        checked: usize,
        /// Human-readable description of the worst disagreement.
        detail: String,
    },
    /// The oracle could not decide (unsupported heads, domain, too few
    /// usable samples). Never counts as verified, never counts as wrong.
    Indeterminate {
        /// Why the check was abandoned.
        reason: String,
    },
}

/// Evaluate `expr` under `env`.
///
/// `log` is evaluated as `ln|x|` so the check is robust to the sign
/// convention of logarithmic partial fractions (the same convention the
/// module-level numeric tests use).
pub fn eval_f64(expr: Atom<'_>, env: &[(Symbol, f64)]) -> Eval {
    match expr.node() {
        AtomNode::Num(n) => Eval::Value(*n as f64),
        AtomNode::Var(v) => {
            if let Some((_, val)) = env.iter().find(|(s, _)| s == v) {
                return Eval::Value(*val);
            }
            // Constants the corpus writes as bare symbols.
            match v.as_str() {
                "pi" => Eval::Value(std::f64::consts::PI),
                "e" | "E" => Eval::Value(std::f64::consts::E),
                _ => Eval::Domain,
            }
        }
        AtomNode::Add(args) => {
            let mut acc = 0.0;
            for a in args.iter() {
                match eval_f64(*a, env) {
                    Eval::Value(v) => acc += v,
                    other => return other,
                }
            }
            finite(acc)
        }
        AtomNode::Mul(args) => {
            let mut acc = 1.0;
            for a in args.iter() {
                match eval_f64(*a, env) {
                    Eval::Value(v) => acc *= v,
                    other => return other,
                }
            }
            finite(acc)
        }
        AtomNode::Pow(b, e) => {
            let base = eval_f64(*b, env);
            let exp = eval_f64(*e, env);
            match (base, exp) {
                (Eval::Value(b), Eval::Value(e)) => {
                    // Integer exponents keep `(-1)^0.5` out of the real path.
                    if e.fract() == 0.0 {
                        finite(b.powf(e))
                    } else if b < 0.0 {
                        Eval::Domain
                    } else {
                        finite(b.powf(e))
                    }
                }
                (Eval::Domain, _) | (_, Eval::Domain) => Eval::Domain,
                (Eval::NonFinite, _) | (_, Eval::NonFinite) => Eval::NonFinite,
                _ => Eval::Unsupported,
            }
        }
        AtomNode::Fun(name, args) => eval_fun(name.as_str(), args, env),
    }
}

/// Report a non-finite result as a zero-division/overflow failure rather than
/// a value.
fn finite(v: f64) -> Eval {
    if v.is_finite() {
        Eval::Value(v)
    } else {
        Eval::NonFinite
    }
}

/// Single-argument function heads the oracle implements.
fn unary(name: &str, v: f64) -> Eval {
    let out = match name {
        "exp" => v.exp(),
        "log" => v.abs().ln(),
        "sqrt" => {
            if v < 0.0 {
                return Eval::Domain;
            }
            v.sqrt()
        }
        "abs" => v.abs(),
        "sin" => v.sin(),
        "cos" => v.cos(),
        "tan" => v.tan(),
        "cot" => 1.0 / v.tan(),
        "sec" => 1.0 / v.cos(),
        "csc" => 1.0 / v.sin(),
        "asin" => {
            if !(-1.0..=1.0).contains(&v) {
                return Eval::Domain;
            }
            v.asin()
        }
        "acos" => {
            if !(-1.0..=1.0).contains(&v) {
                return Eval::Domain;
            }
            v.acos()
        }
        "atan" => v.atan(),
        "acot" => std::f64::consts::FRAC_PI_2 - v.atan(),
        "asec" => {
            if v.abs() < 1.0 {
                return Eval::Domain;
            }
            (1.0 / v).acos()
        }
        "acsc" => {
            if v.abs() < 1.0 {
                return Eval::Domain;
            }
            (1.0 / v).asin()
        }
        "sinh" => v.sinh(),
        "cosh" => v.cosh(),
        "tanh" => v.tanh(),
        "coth" => 1.0 / v.tanh(),
        "sech" => 1.0 / v.cosh(),
        "csch" => 1.0 / v.sinh(),
        "asinh" => v.asinh(),
        "acosh" => {
            if v < 1.0 {
                return Eval::Domain;
            }
            v.acosh()
        }
        "atanh" => {
            if v.abs() >= 1.0 {
                return Eval::Domain;
            }
            v.atanh()
        }
        "acoth" => {
            if v.abs() <= 1.0 {
                return Eval::Domain;
            }
            (1.0 / v).atanh()
        }
        "asech" => {
            if v <= 0.0 || v > 1.0 {
                return Eval::Domain;
            }
            (1.0 / v).acosh()
        }
        "acsch" => {
            if v == 0.0 {
                return Eval::Domain;
            }
            (1.0 / v).asinh()
        }
        "erf" => {
            if v.abs() > 3.0 {
                return Eval::Domain;
            }
            erf_series(v)
        }
        "erfc" => {
            if v.abs() > 3.0 {
                return Eval::Domain;
            }
            1.0 - erf_series(v)
        }
        "erfi" => {
            if v.abs() > 5.0 {
                return Eval::Domain;
            }
            erfi_series(v)
        }
        // Heads the oracle deliberately leaves undecided: they appear in a
        // handful of corpus cases and implementing them accurately is not
        // worth the risk of a false verification.
        _ => return Eval::Unsupported,
    };
    finite(out)
}

fn eval_fun(name: &str, args: &[Atom<'_>], env: &[(Symbol, f64)]) -> Eval {
    match name {
        "EllipticF" | "EllipticE" | "EllipticPi" => elliptic(name, args, env),
        _ => {
            if args.len() != 1 {
                return Eval::Unsupported;
            }
            match eval_f64(args[0], env) {
                Eval::Value(v) => unary(name, v),
                other => other,
            }
        }
    }
}

// ------------------------------------------------------------------
//  Error function family
// ------------------------------------------------------------------

const TWO_OVER_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;

/// Maclaurin series `2/√π · Σ (−1)ⁿ x^(2n+1) / (n! (2n+1))`.
///
/// Accurate to ~1e-14 for `|x| ≤ 3`; the alternating series loses digits
/// beyond that (at `x = 4` the absolute error is ~4e-12, at `x = 6` it is
/// ~7e-3), so callers restrict the domain to `|x| ≤ 3`.
fn erf_series(x: f64) -> f64 {
    let mut term = x;
    let mut sum = x;
    let x2 = x * x;
    let mut n = 1.0;
    for _ in 0..120 {
        term *= -x2 / n;
        let add = term / (2.0 * n + 1.0);
        sum += add;
        if add.abs() < 1e-18 * sum.abs().max(1e-300) {
            break;
        }
        n += 1.0;
    }
    TWO_OVER_SQRT_PI * sum
}

/// `erfi(x) = 2/√π · Σ x^(2n+1) / (n! (2n+1))` (all-positive terms, so no
/// cancellation; restricted to `|x| ≤ 5` to stay finite).
fn erfi_series(x: f64) -> f64 {
    let mut term = x;
    let mut sum = x;
    let x2 = x * x;
    let mut n = 1.0;
    for _ in 0..120 {
        term *= x2 / n;
        let add = term / (2.0 * n + 1.0);
        sum += add;
        if add.abs() < 1e-18 * sum.abs() {
            break;
        }
        n += 1.0;
    }
    TWO_OVER_SQRT_PI * sum
}

// ------------------------------------------------------------------
//  Elliptic integrals (defining-integral quadrature)
// ------------------------------------------------------------------

/// Evaluate `EllipticF(φ, m)`, `EllipticE(φ, m)` or `EllipticPi(n, φ, m)`
/// by adaptive Simpson on the defining integral.
///
/// `m = k²` (SymPy parameter convention). Restricted to `|φ| ≤ π/2` and
/// `|m|, |n| ≤ 0.9`, where the integrand is smooth; outside that range the
/// result is `Domain` (an honest "cannot decide", never a guess).
fn elliptic(name: &str, args: &[Atom<'_>], env: &[(Symbol, f64)]) -> Eval {
    let (n_arg, phi_arg, m_arg) = match (name, args) {
        ("EllipticPi", [n, phi, m]) => (Some(*n), *phi, *m),
        ("EllipticF" | "EllipticE", [phi, m]) => (None, *phi, *m),
        _ => return Eval::Unsupported,
    };
    let n_val = match n_arg {
        Some(a) => match eval_f64(a, env) {
            Eval::Value(v) => v,
            other => return other,
        },
        None => 0.0,
    };
    let phi = match eval_f64(phi_arg, env) {
        Eval::Value(v) => v,
        other => return other,
    };
    let m = match eval_f64(m_arg, env) {
        Eval::Value(v) => v,
        other => return other,
    };
    if phi.abs() > std::f64::consts::FRAC_PI_2 + 1e-12 || m.abs() > 0.9 || n_val.abs() > 0.9 {
        return Eval::Domain;
    }
    let f = |theta: f64| -> f64 {
        let s = theta.sin();
        let root = (1.0 - m * s * s).sqrt();
        match name {
            "EllipticF" => 1.0 / root,
            "EllipticE" => root,
            _ => 1.0 / ((1.0 - n_val * s * s) * root),
        }
    };
    finite(adaptive_simpson(&f, 0.0, phi, 1e-13, 40))
}

/// Adaptive Simpson quadrature with a relative/absolute tolerance and a
/// recursion cap. Returns `f64::NAN` when the cap is hit without converging.
fn adaptive_simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64, depth: u32) -> f64 {
    fn simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
        let m = 0.5 * (a + b);
        (b - a) / 6.0 * (f(a) + 4.0 * f(m) + f(b))
    }
    fn rec(f: &dyn Fn(f64) -> f64, a: f64, b: f64, whole: f64, tol: f64, depth: u32) -> f64 {
        let m = 0.5 * (a + b);
        let left = simpson(f, a, m);
        let right = simpson(f, m, b);
        let delta = left + right - whole;
        if depth == 0 || delta.abs() <= 15.0 * tol {
            return left + right + delta / 15.0;
        }
        rec(f, a, m, left, 0.5 * tol, depth - 1) + rec(f, m, b, right, 0.5 * tol, depth - 1)
    }
    if a == b {
        return 0.0;
    }
    let whole = simpson(f, a, b);
    rec(f, a, b, whole, tol, depth)
}

// ------------------------------------------------------------------
//  Antiderivative verification
// ------------------------------------------------------------------

/// Sample abscissae used by [`verify_antiderivative`]. Chosen to stay away
/// from the common singularities of the corpus (0, ±1) while covering both
/// signs, both sides of 1, and the region past π/2 ≈ 1.57 where the
/// half-power/elliptic sign corrections matter. A wider spread matters for
/// the *domain* rather than the accuracy: inverse-hyperbolic and radical
/// heads are real only on part of the line, so more abscissae mean more
/// cases get a genuine comparison instead of "no usable sample".
const SAMPLES: [f64; 16] = [
    -2.7, -2.3, -1.9, -1.55, -1.2, -0.85, -0.55, -0.3, 0.3, 0.55, 0.85, 1.2, 1.55, 1.9, 2.3, 2.7,
];

/// Parameters whose real numeric treatment is either impossible or wrong.
fn is_special_symbol(name: &str) -> bool {
    matches!(name, "i" | "I" | "inf" | "oo" | "nan")
}

/// Collect free variables of `expr`, excluding `var` and the constants the
/// evaluator handles itself.
pub fn free_symbols(expr: Atom<'_>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Var(v) => {
            if *v != var && !matches!(v.as_str(), "pi" | "e" | "E") && !out.contains(v) {
                out.push(*v);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            for a in args.iter() {
                free_symbols(*a, var, out);
            }
        }
        AtomNode::Pow(b, e) => {
            free_symbols(*b, var, out);
            free_symbols(*e, var, out);
        }
        AtomNode::Num(_) => {}
    }
}

/// Deterministic dummy value for a parameter, keyed by its name, a regime and
/// a salt so repeated runs assign identical values.
///
/// The magnitude stays modest so that function arguments (e.g. `b*x`) do not
/// leave the oracle's domains. `regime` selects the sign policy: regime 0 is
/// mixed-sign, regimes 1 and 2 are all-positive with distinct values — the
/// elliptic and radical closed forms carry `1/sqrt(a + b)`-style factors that
/// are real only when the parameters keep the radicand positive, so a
/// positive round is what lets those emissions be *verified* rather than
/// merely undecided.
fn param_value(name: &str, regime: usize, salt: u64) -> f64 {
    const MIXED: [f64; 12] = [
        2.0,
        3.0,
        5.0,
        0.5,
        1.5,
        0.25,
        4.0,
        1.0 / 3.0,
        1.0,
        2.5,
        0.75,
        6.0,
    ];
    const POSITIVE_SMALL: [f64; 12] = [
        0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 0.75, 1.25, 0.25, 6.0,
    ];
    const POSITIVE_LARGE: [f64; 12] = [
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
    ];
    let table: &[f64; 12] = match regime {
        0 => &MIXED,
        1 => &POSITIVE_SMALL,
        _ => &POSITIVE_LARGE,
    };
    let mut h: u64 = 0xcbf2_9ce4_8422_2325 ^ salt ^ ((regime as u64) << 32);
    for b in name.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    let v = table[(h % 12) as usize];
    if regime == 0 && (h >> 8) & 1 == 1 {
        -v
    } else {
        v
    }
}

/// One verification round with a fixed parameter regime and salt.
fn verify_round(
    integrand: Atom<'_>,
    antiderivative: Atom<'_>,
    var: Symbol,
    regime: usize,
    salt: u64,
) -> Verify {
    let mut params: Vec<Symbol> = Vec::new();
    free_symbols(integrand, var, &mut params);
    free_symbols(antiderivative, var, &mut params);
    if params.iter().any(|s| is_special_symbol(s.as_str())) {
        return Verify::Indeterminate {
            reason: "non-real parameter (i/inf)".to_string(),
        };
    }
    let base_env: Vec<(Symbol, f64)> = params
        .iter()
        .map(|s| (*s, param_value(s.as_str(), regime, salt)))
        .collect();

    let mut checked = 0usize;
    let mut unstable = 0usize;
    let mut nonfinite_anti = 0usize;
    let mut worst_rel = 0.0f64;
    let mut worst_detail = String::new();

    for &x in SAMPLES.iter() {
        let h = 1e-4 * x.abs().max(1.0);
        let eval_at = |t: f64| -> Eval {
            let mut env = base_env.clone();
            env.push((var, t));
            eval_f64(antiderivative, &env)
        };
        let mut env = base_env.clone();
        env.push((var, x));
        let rhs = eval_f64(integrand, &env);
        // An unsupported head anywhere makes the whole round undecidable.
        if matches!(rhs, Eval::Unsupported) {
            return Verify::Indeterminate {
                reason: "integrand uses a head the oracle does not implement".to_string(),
            };
        }
        let pts = [
            eval_at(x - 4.0 * h),
            eval_at(x - 2.0 * h),
            eval_at(x - h),
            eval_at(x + h),
            eval_at(x + 2.0 * h),
            eval_at(x + 4.0 * h),
        ];
        if pts.iter().any(|p| matches!(p, Eval::Unsupported)) {
            return Verify::Indeterminate {
                reason: "antiderivative uses a head the oracle does not implement".to_string(),
            };
        }
        // A finite integrand with a non-finite antiderivative is a broken
        // closed form (typically a coefficient like `1/(d - d)`), not an
        // unusable sample: count it so the case is reported as wrong.
        if matches!(rhs, Eval::Value(_)) && pts.iter().any(|p| matches!(p, Eval::NonFinite)) {
            nonfinite_anti += 1;
            continue;
        }
        let Some(vals) = pts
            .iter()
            .map(|p| match p {
                Eval::Value(v) => Some(*v),
                _ => None,
            })
            .collect::<Option<Vec<f64>>>()
        else {
            continue;
        };
        let d_h = (-vals[4] + 8.0 * vals[3] - 8.0 * vals[2] + vals[1]) / (12.0 * h);
        let d_2h = (-vals[5] + 8.0 * vals[4] - 8.0 * vals[1] + vals[0]) / (24.0 * h);
        let Eval::Value(f) = rhs else { continue };
        if !d_h.is_finite() || !d_2h.is_finite() {
            continue;
        }
        // Self-consistency gate: when the antiderivative is a ratio of nearly
        // cancelling quantities (e.g. `(1-2x)^3/((2+3x)^7(3+5x)^2)` at
        // x=-0.9) the finite difference is dominated by round-off and is not a
        // valid check. Two step sizes agreeing with each other is the
        // precondition for using either of them.
        if (d_h - d_2h).abs() / f.abs().max(1.0) > 1e-5 {
            unstable += 1;
            continue;
        }
        checked += 1;
        let rel = (d_h - f).abs() / f.abs().max(1.0);
        if rel > worst_rel {
            worst_rel = rel;
            worst_detail = format!("x={x}: d/dx F = {d_h}, integrand = {f}, rel err = {rel:.3e}");
        }
    }

    if checked < 2 {
        if nonfinite_anti >= 3 && nonfinite_anti == SAMPLES.len() {
            return Verify::Mismatch {
                checked: 0,
                detail: format!(
                    "antiderivative is non-finite at all {nonfinite_anti} samples while the \
                     integrand is finite (a coefficient that divides by an identically zero \
                     parameter expression?)"
                ),
            };
        }
        return Verify::Indeterminate {
            reason: format!(
                "only {checked} usable sample point(s) ({unstable} unstable, \
                 {nonfinite_anti} non-finite)"
            ),
        };
    }
    if worst_rel > 1e-4 {
        return Verify::Mismatch {
            checked,
            detail: worst_detail,
        };
    }
    if worst_rel > 1e-5 {
        return Verify::Indeterminate {
            reason: format!("inconclusive: worst relative error {worst_rel:.3e}"),
        };
    }
    Verify::Verified { checked, worst_rel }
}

/// Verify that `antiderivative` differentiates back to `integrand`.
///
/// Uses a 5-point central difference with `h = 1e-4·max(1, |x|)`. Three
/// parameter rounds are attempted — mixed-sign, then two all-positive regimes
/// (radical and elliptic forms are only real in the latter). Every round is
/// evaluated and the best evidence wins: any round that verifies makes the
/// case `Verified`, otherwise any round that disagrees makes it a `Mismatch`,
/// otherwise it is `Indeterminate`.
///
/// Trying *all* rounds matters for correctness of the verdict: a valid
/// generic formula can have a pole at an isolated parameter value (Rubi's
/// `∫(d+e·x)^m dx = (d+e·x)^(m+1)/(e(1+m))` is undefined at `m = −1`, where
/// the antiderivative is a logarithm). Stopping at the first round would
/// report such a formula as wrong when it is merely specialised, while a
/// genuinely wrong formula fails every round.
pub fn verify_antiderivative(integrand: Atom<'_>, antiderivative: Atom<'_>, var: Symbol) -> Verify {
    let mut best: Option<Verify> = None;
    for regime in 0..3usize {
        let round = verify_round(integrand, antiderivative, var, regime, 0);
        match round {
            Verify::Verified { .. } => return round,
            Verify::Mismatch { .. } => {
                if best.is_none() {
                    best = Some(round);
                }
            }
            other @ Verify::Indeterminate { .. } => {
                if best.is_none() {
                    best = Some(other);
                }
            }
        }
    }
    best.unwrap_or(Verify::Indeterminate {
        reason: "no parameter round produced a comparison".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    fn parse<'a>(ctx: &'a AtomArena<'a>, s: &str) -> Atom<'a> {
        ocas::prelude::parse(ctx, s).expect("parse")
    }

    #[test]
    fn oracle_accepts_a_correct_antiderivative() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let f = parse(&ctx, "x^2 + sin(x)");
        let big_f = parse(&ctx, "(3^-1)*x^3 - cos(x)");
        match verify_antiderivative(f, big_f, Symbol::new("x")) {
            Verify::Verified { checked, worst_rel } => {
                assert!(checked >= 2);
                assert!(worst_rel < 1e-5);
            }
            other => panic!("expected Verified, got {other:?}"),
        }
    }

    #[test]
    fn oracle_rejects_a_wrong_antiderivative() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let f = parse(&ctx, "sin(x)");
        // +cos(x) is a sign error: must be reported, not accepted.
        let big_f = parse(&ctx, "cos(x)");
        assert!(matches!(
            verify_antiderivative(f, big_f, Symbol::new("x")),
            Verify::Mismatch { .. }
        ));
    }

    #[test]
    fn oracle_is_indeterminate_for_special_heads_and_constants() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // Ei is deliberately unimplemented: honest indeterminate.
        let f = parse(&ctx, "Ei(x)");
        let big_f = parse(&ctx, "x");
        assert!(matches!(
            verify_antiderivative(f, big_f, Symbol::new("x")),
            Verify::Indeterminate { .. }
        ));
        // The imaginary unit cannot be sampled on the real line.
        let g = parse(&ctx, "x^i");
        let g_anti = parse(&ctx, "x^(1+i)/(1+i)");
        assert!(matches!(
            verify_antiderivative(g, g_anti, Symbol::new("x")),
            Verify::Indeterminate { .. }
        ));
    }

    #[test]
    fn erf_family_matches_known_values() {
        // References evaluated with SymPy to 20 digits.
        assert!((erf_series(0.5) - 0.520_499_877_813_046_5).abs() < 1e-14);
        assert!((erf_series(2.0) - 0.995_322_265_018_952_7).abs() < 1e-14);
        assert!((erf_series(3.0) - 0.999_977_909_503_001_4).abs() < 1e-13);
        // erfc(x) = 1 − erf(x) is a difference of nearly equal numbers here:
        // absolute error ~3e-14 at x = 3 (relative ~1e-9).
        assert!((1.0 - erf_series(3.0) - 2.209_049_699_858_544e-5).abs() < 1e-13);
        // Odd symmetry.
        assert!((erf_series(0.5) + erf_series(-0.5)).abs() < 1e-15);
        // erfi is NOT e^{x²}·erf(x) (that is the Dawson identity); these are
        // the SymPy values.
        assert!((erfi_series(0.5) - 0.614_952_094_696_511).abs() < 1e-14);
        assert!((erfi_series(1.0) - 1.650_425_758_797_543).abs() < 1e-13);
    }

    #[test]
    fn elliptic_quadrature_matches_known_values() {
        // K(1/2) = F(π/2, 1/2) = 1.8540746773013719
        let k = adaptive_simpson(
            &|t: f64| 1.0 / (1.0 - 0.5 * t.sin() * t.sin()).sqrt(),
            0.0,
            std::f64::consts::FRAC_PI_2,
            1e-13,
            40,
        );
        assert!((k - 1.854_074_677_301_372).abs() < 1e-9);
        // E(π/2, 1/2) = 1.3506438810476755
        let e = adaptive_simpson(
            &|t: f64| (1.0 - 0.5 * t.sin() * t.sin()).sqrt(),
            0.0,
            std::f64::consts::FRAC_PI_2,
            1e-13,
            40,
        );
        assert!((e - 1.350_643_881_047_675_5).abs() < 1e-9);
    }

    #[test]
    fn elliptic_heads_evaluate_through_the_atom_path() {
        // References from SymPy (17 digits): elliptic_f(pi/3, 1/2),
        // elliptic_e(pi/3, 1/2), elliptic_pi(3/10, pi/3, 1/2). The harness
        // oracle dispatches the elliptic heads to the defining-integral
        // quadrature, so these pin the argument convention (m = k², φ the
        // amplitude) as well as the numerics.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let cases: [(&str, f64); 3] = [
            ("EllipticF(pi*(3^-1), (2^-1))", 1.142_429_058_045_777_3),
            ("EllipticE(pi*(3^-1), (2^-1))", 0.964_951_457_642_992_6),
            (
                "EllipticPi((3*(10^-1)), pi*(3^-1), (2^-1))",
                1.268_121_649_431_691_6,
            ),
        ];
        for (src, want) in cases {
            let expr = parse(&ctx, src);
            match eval_f64(expr, &[]) {
                Eval::Value(v) => assert!((v - want).abs() < 1e-9, "{src}: got {v}, want {want}"),
                other => panic!("{src}: expected a value, got {other:?}"),
            }
        }
        // Argument order is semantic: swapping amplitude and parameter must
        // change the value (EllipticF(x, m) ≠ EllipticF(m, x) in general).
        let f1 = parse(&ctx, "EllipticF(x, (2^-1))");
        let f2 = parse(&ctx, "EllipticF((2^-1), x)");
        assert_ne!(f1.to_string(), f2.to_string());
        let env = [(Symbol::new("x"), 0.7)];
        let (a, b) = (eval_f64(f1, &env), eval_f64(f2, &env));
        assert!(matches!((a, b), (Eval::Value(_), Eval::Value(_))));
        assert_ne!(a, b);
    }
}
