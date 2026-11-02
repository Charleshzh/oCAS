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
//! whose functions the oracle cannot evaluate (heads outside the table, the
//! imaginary unit, domains it cannot sample) is `Indeterminate`. The harness
//! reports the three counts separately and never merges them.
//!
//! The special-function table (0.27.3) covers `erf`/`erfc`/`erfi`, the
//! exponential-integral family (`Ei`, `Ei(n, z)` = `Eₙ(z)`), the trigonometric
//! integrals (`Si`, `Ci`, `Shi`, `Chi`) and the Fresnel integrals
//! (`fresnels`, `fresnelc`). 0.31a added the dilogarithm head
//! `polylog(2, z)` = `Li₂(z)`, real-valued on `z ≤ 1`; `z > 1` is outside the
//! real domain and is reported as `Domain` (the sample is skipped) instead of
//! `Unsupported` (which would abort the whole verification round). Every
//! algorithm here was checked against `mpmath`/SymPy at 40 digits before being
//! written down; the crossovers and the branch conventions are documented at
//! their definitions. Heads that are *not* implemented stay `Unsupported`,
//! which is an honest "cannot decide" rather than a guess.

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
        // Exponential-integral family. `Ei` (one argument) is the classical
        // real exponential integral, real on the whole line; `Ei(n, z)` is a
        // different function and is dispatched in `eval_fun`.
        "Ei" => {
            if v == 0.0 {
                return Eval::Domain;
            }
            return ei(v);
        }
        "Si" => si(v),
        "Ci" => {
            if v == 0.0 {
                return Eval::Domain;
            }
            ci(v)
        }
        "Shi" => {
            if v == 0.0 {
                return Eval::Domain;
            }
            shi(v)
        }
        "Chi" => {
            if v == 0.0 {
                return Eval::Domain;
            }
            chi(v)
        }
        "fresnels" => fresnel_s(v),
        "fresnelc" => fresnel_c(v),
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
        "Ei" if args.len() == 2 => ei_order(args, env),
        "polylog" => polylog(args, env),
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
//  Exponential integral / trigonometric integral / Fresnel family
// ------------------------------------------------------------------

/// Euler–Mascheroni constant `γ` (25 digits).
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// `Si`/`Ci` crossover between the Taylor series and the asymptotic auxiliary
/// series. Measured against `mpmath` at 40 digits: the series is ~4e-11 at
/// `x = 18` and degrades to 5e-10 at 20, while the asymptotic series reaches
/// 7.5e-10 at 20 and 4.6e-12 at 24, so 20 is where the two curves cross and
/// neither side is worse than ~1e-9 (a derivative error of ~1e-5, i.e. at the
/// oracle's own tolerance).
const SI_CI_SERIES_MAX: f64 = 20.0;

/// Fresnel crossover. The Taylor series is unusable past `x ≈ 4` (its terms
/// grow like `x^{4k}`, so it cancels catastrophically), while the asymptotic
/// auxiliary series is already accurate to 7e-13 at `x = 4` and 1e-16 above 6.
const FRESNEL_SERIES_MAX: f64 = 3.0;

/// `Ei(x)` for `x > 0`: `γ + ln x + Σ_{k≥1} x^k/(k·k!)`.
///
/// Every term is positive, so there is no cancellation; the series is good to
/// ~1e-16 relative over the whole range the oracle samples (`x ≲ 50`).
fn ei_positive(x: f64) -> f64 {
    let mut term = 1.0_f64;
    let mut sum = 0.0_f64;
    for k in 1..=900u32 {
        let kf = f64::from(k);
        term *= x / kf;
        let add = term / kf;
        sum += add;
        if kf > x + 25.0 && add.abs() < 1e-20 * sum.abs() {
            break;
        }
    }
    EULER_GAMMA + x.ln() + sum
}

/// `E₁(y)` for `y > 0`.
///
/// Two cheap, allocation-free branches instead of a quadrature — the oracle
/// evaluates this many thousands of times per corpus case, and an adaptive
/// Simpson rule with a 1e-16 tolerance made the *debug* test suite take minutes
/// (a single 14-integrand verification test went from 0.4 s to 376 s):
///
/// - `y ≤ 15`: the convergent series `E₁(y) = −γ − ln y + Σ_{k≥1} (−1)^{k+1}
///   y^k/(k·k!)`. Cancellation stays mild in this range (relative error ≈2e-7 at
///   the cutoff), and `E₁(y) ≤ 1e-5` there, so the absolute error is ≤1e-12 —
///   which is what matters, because the verification compares against
///   `max(|integrand|, 1)`.
/// - `y > 15`: the asymptotic series `e^{−y}/y · Σ (−1)^k k!/y^k`, truncated at
///   its smallest term. `E₁` is at most 1e-8 there, and it only ever appears
///   next to `Ei(y)` terms that dwarf it (`Shi`/`Chi`) or as the whole value of
///   `Ei(−y)` (≤1e-8 absolute), so the truncation error is far below the
///   oracle's tolerance either way.
///
/// Both branches are a few dozen flops; the previous quadrature was the single
/// hot spot of every debug-mode suite run.
fn e1_positive(y: f64) -> f64 {
    if y <= 15.0 {
        let mut term = 1.0_f64;
        let mut sum = 0.0_f64;
        for k in 1..=200u32 {
            let kf = f64::from(k);
            term *= -y / kf;
            let add = term / kf;
            sum += add;
            if add.abs() < 1e-22 * sum.abs().max(1e-300) {
                break;
            }
        }
        return -EULER_GAMMA - y.ln() - sum;
    }
    let mut term = 1.0_f64;
    let mut sum = 1.0_f64;
    let mut prev = f64::INFINITY;
    for k in 1..=200u32 {
        term *= -f64::from(k) / y;
        let magnitude = term.abs();
        if magnitude > prev {
            break;
        }
        sum += term;
        prev = magnitude;
        if magnitude < 1e-18 {
            break;
        }
    }
    (-y).exp() / y * sum
}

/// `Ei(x)` for any real `x ≠ 0`.
fn ei(x: f64) -> Eval {
    let v = if x > 0.0 {
        ei_positive(x)
    } else {
        // Ei(−y) = −E₁(y) is real and exponentially small.
        -e1_positive(-x)
    };
    if v.is_finite() {
        Eval::Value(v)
    } else {
        Eval::Domain
    }
}

/// `(f, g)` for the trigonometric integrals:
/// `Si(x) = π/2 − f cos x − g sin x`, `Ci(x) = f sin x − g cos x`, with the
/// asymptotic auxiliary series
/// `f ~ Σ (−1)^k (2k)!/x^{2k+1}`, `g ~ Σ (−1)^k (2k+1)!/x^{2k+2}`.
///
/// Truncated at the smallest term: the series is asymptotic (divergent), so
/// adding terms past that point only adds error.
fn trig_aux(x: f64) -> (f64, f64) {
    let inv2 = 1.0 / (x * x);
    let mut tf = 1.0 / x;
    let mut tg = inv2;
    let (mut f, mut g) = (0.0_f64, 0.0_f64);
    let (mut prev_f, mut prev_g) = (f64::INFINITY, f64::INFINITY);
    for k in 0..80u32 {
        if k > 0 {
            let kf = f64::from(k);
            tf = -tf * (2.0 * kf) * (2.0 * kf - 1.0) * inv2;
            tg = -tg * (2.0 * kf + 1.0) * (2.0 * kf) * inv2;
        }
        if tf.abs() > prev_f || tg.abs() > prev_g {
            break;
        }
        f += tf;
        g += tg;
        prev_f = tf.abs();
        prev_g = tg.abs();
    }
    (f, g)
}

/// `Si(x) = Σ_{k≥0} (−1)^k x^{2k+1}/((2k+1)!(2k+1))`.
fn si_series(x: f64) -> f64 {
    let mut p = x;
    let mut sum = 0.0_f64;
    for k in 0..900u32 {
        let kf = f64::from(k);
        let add = p / (2.0 * kf + 1.0);
        sum += if k % 2 == 0 { add } else { -add };
        p *= x * x / ((2.0 * kf + 2.0) * (2.0 * kf + 3.0));
        if (p / (2.0 * kf + 3.0)).abs() < 1e-24 * sum.abs().max(1e-300) {
            break;
        }
    }
    sum
}

/// `Ci(x) = γ + ln|x| + Σ_{k≥1} (−1)^k x^{2k}/((2k)!(2k))` (even in `x`).
fn ci_series(x: f64) -> f64 {
    let mut t = 1.0_f64;
    let mut sum = 0.0_f64;
    for k in 1..=900u32 {
        let kf = f64::from(k);
        t *= -x * x / ((2.0 * kf - 1.0) * (2.0 * kf));
        let add = t / (2.0 * kf);
        sum += add;
        if add.abs() < 1e-24 * sum.abs().max(1e-300) {
            break;
        }
    }
    EULER_GAMMA + x.abs().ln() + sum
}

/// `Si(x)`, odd, real on the whole line.
fn si(x: f64) -> f64 {
    let y = x.abs();
    let v = if y <= SI_CI_SERIES_MAX {
        si_series(y)
    } else {
        let (f, g) = trig_aux(y);
        std::f64::consts::FRAC_PI_2 - f * y.cos() - g * y.sin()
    };
    if x < 0.0 { -v } else { v }
}

/// `Ci(x)`, even, real on the whole line (the real convention
/// `Ci(x) = γ + ln|x| + ∫_0^x (cos t − 1)/t dt`; SymPy's principal branch
/// differs from it by the constant `iπ` on the negative axis, which cancels in
/// every derivative the oracle checks).
fn ci(x: f64) -> f64 {
    let y = x.abs();
    if y <= SI_CI_SERIES_MAX {
        ci_series(y)
    } else {
        let (f, g) = trig_aux(y);
        f * y.sin() - g * y.cos()
    }
}

/// `Shi(x) = (Ei(x) + E₁(|x|))/2`, odd.
fn shi(x: f64) -> f64 {
    let y = x.abs();
    let v = 0.5 * (ei_positive(y) + e1_positive(y));
    if x < 0.0 { -v } else { v }
}

/// `Chi(x) = (Ei(x) − E₁(|x|))/2`, even.
fn chi(x: f64) -> f64 {
    let y = x.abs();
    0.5 * (ei_positive(y) - e1_positive(y))
}

/// `(f, g)` for the Fresnel integrals, with the asymptotic auxiliary series
///
/// ```text
/// f(x) ~ (1/(πx)) Σ (−1)^k a_k/(π^{2k} x^{4k}),  a_k = Π_{j≤k} (4j−3)(4j−1)
/// g(x) ~ (1/(π²x³)) Σ (−1)^k b_k/(π^{2k} x^{4k}), b_k = Π_{j≤k} (4j−1)(4j+1)
/// ```
///
/// derived by repeated integration by parts of `∫_x^∞ sin(πt²/2) dt` (the
/// coefficients were derived symbolically here and confirmed against `mpmath`:
/// 1e-16 at `x = 12`, and better for larger `x`). Truncated at the smallest
/// term, like [`trig_aux`].
fn fresnel_aux(x: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    let step = 1.0 / (pi * pi * x * x * x * x);
    let mut tf = 1.0_f64;
    let mut tg = 1.0_f64;
    let (mut sf, mut sg) = (0.0_f64, 0.0_f64);
    let (mut prev_f, mut prev_g) = (f64::INFINITY, f64::INFINITY);
    for k in 0..80u32 {
        if k > 0 {
            let kf = f64::from(k);
            tf = -tf * (4.0 * kf - 3.0) * (4.0 * kf - 1.0) * step;
            tg = -tg * (4.0 * kf - 1.0) * (4.0 * kf + 1.0) * step;
        }
        if tf.abs() > prev_f || tg.abs() > prev_g {
            break;
        }
        sf += tf;
        sg += tg;
        prev_f = tf.abs();
        prev_g = tg.abs();
    }
    (sf / (pi * x), sg / (pi * pi * x * x * x))
}

/// `fresnels(x) = Σ (−1)^k (π/2)^{2k+1} x^{4k+3}/((2k+1)!(4k+3))` for small `x`.
fn fresnel_s_series(x: f64) -> f64 {
    let a = std::f64::consts::FRAC_PI_2;
    let x4 = x.powi(4);
    let mut t = a * x * x * x / 3.0;
    let mut sum = 0.0_f64;
    for k in 0..80u32 {
        if k > 0 {
            let kf = f64::from(k);
            t = -t * a * a * x4 / ((2.0 * kf) * (2.0 * kf + 1.0)) * (4.0 * kf - 1.0)
                / (4.0 * kf + 3.0);
        }
        sum += t;
        if t.abs() < 1e-24 * sum.abs().max(1e-300) {
            break;
        }
    }
    sum
}

/// `fresnelc(x) = Σ (−1)^k (π/2)^{2k} x^{4k+1}/((2k)!(4k+1))` for small `x`.
fn fresnel_c_series(x: f64) -> f64 {
    let a = std::f64::consts::FRAC_PI_2;
    let x4 = x.powi(4);
    let mut t = x;
    let mut sum = 0.0_f64;
    for k in 0..80u32 {
        if k > 0 {
            let kf = f64::from(k);
            t = -t * a * a * x4 / ((2.0 * kf - 1.0) * (2.0 * kf)) * (4.0 * kf - 3.0)
                / (4.0 * kf + 1.0);
        }
        sum += t;
        if t.abs() < 1e-24 * sum.abs().max(1e-300) {
            break;
        }
    }
    sum
}

/// `fresnels(x)`, odd (SymPy convention `S(z) = ∫_0^z sin(πt²/2) dt`).
fn fresnel_s(x: f64) -> f64 {
    let y = x.abs();
    let v = if y <= FRESNEL_SERIES_MAX {
        fresnel_s_series(y)
    } else {
        let th = std::f64::consts::FRAC_PI_2 * y * y;
        let (f, g) = fresnel_aux(y);
        0.5 - f * th.cos() - g * th.sin()
    };
    if x < 0.0 { -v } else { v }
}

/// `fresnelc(x)`, odd.
fn fresnel_c(x: f64) -> f64 {
    let y = x.abs();
    let v = if y <= FRESNEL_SERIES_MAX {
        fresnel_c_series(y)
    } else {
        let th = std::f64::consts::FRAC_PI_2 * y * y;
        let (f, g) = fresnel_aux(y);
        0.5 + f * th.sin() - g * th.cos()
    };
    if x < 0.0 { -v } else { v }
}

/// `n!` for the small orders the exponential-integral recurrences need.
fn factorial(n: u32) -> f64 {
    (1..=n).map(f64::from).product()
}

/// `Ei(n, z)` = `Eₙ(z)`, the exponential integral of order `n`
/// (the Rubi/Mathematica convention, matching `sympy.expint(n, z)`).
///
/// - `n = 0`: `E₀(z) = e^{−z}/z`.
/// - `n < 0`: the closed elementary form
///   `E_{−m}(z) = m! e^{−z} Σ_{k≤m} z^{k−m−1}/k!`, real for every `z ≠ 0`.
/// - `n ≥ 1`, `z > 0`: `E₁(z)` followed by the upward recurrence
///   `E_{n+1}(z) = (e^{−z} − z Eₙ(z))/n`.
/// - `n ≥ 1`, `z ≤ 0`: **declined**. The analytic continuation of `Eₙ` across
///   the negative axis is complex (it differs from the real expression by the
///   constant `iπ`), so a real-valued comparison there would be meaningless.
///   The all-positive verification regimes still have eight usable samples.
fn ei_order(args: &[Atom<'_>], env: &[(Symbol, f64)]) -> Eval {
    let order = match eval_f64(args[0], env) {
        Eval::Value(v) => v,
        other => return other,
    };
    let z = match eval_f64(args[1], env) {
        Eval::Value(v) => v,
        other => return other,
    };
    if order != order.trunc() || order.abs() > 12.0 {
        return Eval::Domain;
    }
    let order = order as i64;
    if order == 0 {
        if z == 0.0 {
            return Eval::Domain;
        }
        return finite((-z).exp() / z);
    }
    if order < 0 {
        let m = (-order) as u32;
        let mut sum = 0.0_f64;
        for k in 0..=m {
            sum += z.powi(k as i32 - m as i32 - 1) / factorial(k);
        }
        return finite((-z).exp() * factorial(m) * sum);
    }
    if z <= 0.0 {
        return Eval::Domain;
    }
    let mut e = e1_positive(z);
    for k in 1..order {
        e = ((-z).exp() - z * e) / k as f64;
    }
    if e.is_finite() {
        Eval::Value(e)
    } else {
        Eval::Domain
    }
}

// ------------------------------------------------------------------
//  Dilogarithm (`polylog(2, z)` = Li₂)
// ------------------------------------------------------------------

/// `ζ(2) = π²/6 = Li₂(1)`.
const LI2_AT_ONE: f64 = std::f64::consts::PI * std::f64::consts::PI / 6.0;

/// `Li₂(−1) = −π²/12`.
const LI2_AT_MINUS_ONE: f64 = -std::f64::consts::PI * std::f64::consts::PI / 12.0;

/// Maximum number of terms of the Maclaurin series used for `|z| ≤ 0.5`
/// (SPEC §5's head-evaluation cap of 4096 terms). At `|z| = 0.5` the series
/// already converges to ~1e-16 in 50 terms, so the cap is a safety net rather
/// than the working regime.
const LI2_SERIES_MAX_TERMS: u32 = 4096;

/// Maclaurin series `Li₂(z) = Σ_{k≥1} z^k/k²`.
///
/// Fast and cancellation-free for `|z| ≤ 0.5` (including the alternating
/// negative half-line: the worst grid point `z = −0.9` needs 289 terms and is
/// accurate to ~4e-16). The `k²` denominator makes the series converge only
/// like `1/k` as `|z| → 1`, hence the [`LI2_SERIES_MAX_TERMS`] cap on the
/// branches that can reach that region.
fn li2_series(z: f64) -> f64 {
    let mut term = z;
    let mut sum = z;
    for k in 2..=LI2_SERIES_MAX_TERMS {
        let kf = f64::from(k);
        term *= z;
        let add = term / (kf * kf);
        sum += add;
        // 1e-18 relative is far below the oracle's 1e-12 accuracy target.
        if add.abs() <= 1e-18 * sum.abs() {
            break;
        }
    }
    sum
}

/// `polylog(s, z)`: real-valued dilogarithm `Li₂(z)`, order `s = 2` only.
///
/// The order argument is evaluated (not pattern-matched on the atom) and must
/// be exactly 2; every other order — and every arity other than two — is an
/// honest [`Eval::Unsupported`], matching SPEC §5 ("order rejected at the
/// head").
fn polylog(args: &[Atom<'_>], env: &[(Symbol, f64)]) -> Eval {
    if args.len() != 2 {
        return Eval::Unsupported;
    }
    let order = match eval_f64(args[0], env) {
        Eval::Value(v) => v,
        other => return other,
    };
    if order != 2.0 {
        return Eval::Unsupported;
    }
    match eval_f64(args[1], env) {
        Eval::Value(z) => li2(z),
        other => other,
    }
}

/// Real-valued dilogarithm `Li₂(z)` on its real domain `z ≤ 1`.
///
/// Branch selection (each formula was checked against `mpmath` at 40 digits
/// before being written down, per SPEC §3):
///
/// - `|z| ≤ 0.5`: direct Maclaurin series ([`li2_series`]).
/// - `0.5 < z < 1`: Euler's reflection
///   `Li₂(z) = π²/6 − ln z·ln(1−z) − Li₂(1−z)`; the argument `1−z ∈ (0, 0.5)`
///   goes back to the series.
/// - `−1 < z ≤ −0.5`: the series itself (SPEC §3: "−1 < z < 0: the series
///   converges directly, no branch problem"). Its `1/k` tail limits accuracy
///   only in a sliver next to `z = −1` (≈6e-10 at `z = −0.999`, worse closer
///   in); exact `z = −1` is the constant branch below.
/// - `z < −1`: inversion `Li₂(z) = −π²/6 − ½·ln²(−z) − Li₂(1/z)`, where
///   `1/z ∈ (−1, 0)` is handled by the previous two bullets. One recursive
///   call, so the depth is at most two.
/// - `z = −1` is the exact constant `−π²/12` (the inversion identity is
///   self-referential there, and this is also required to terminate the
///   recursion).
/// - `z = 1` is the exact constant `π²/6`; `z > 1` is **outside the real
///   domain** (Li₂ acquires an imaginary part) and is reported as
///   [`Eval::Domain`], so [`verify_antiderivative`] skips that sample instead
///   of aborting the whole round the way [`Eval::Unsupported`] would.
///
/// Non-finite arguments are [`Eval::NonFinite`].
fn li2(z: f64) -> Eval {
    if !z.is_finite() {
        return Eval::NonFinite;
    }
    let v = if z == 1.0 {
        LI2_AT_ONE
    } else if z > 1.0 {
        return Eval::Domain;
    } else if z == -1.0 {
        LI2_AT_MINUS_ONE
    } else if z.abs() <= 0.5 {
        li2_series(z)
    } else if z > 0.5 {
        LI2_AT_ONE - z.ln() * (1.0 - z).ln() - li2_series(1.0 - z)
    } else if z < -1.0 {
        let inner = match li2(1.0 / z) {
            Eval::Value(v) => v,
            other => return other,
        };
        -LI2_AT_ONE - 0.5 * (-z).ln().powi(2) - inner
    } else {
        // −1 < z < −0.5: the series converges directly.
        li2_series(z)
    };
    finite(v)
}

// ------------------------------------------------------------------
//  Elliptic integrals (defining-integral quadrature)
// ------------------------------------------------------------------

/// Evaluate `EllipticF(φ, m)`, `EllipticE(φ, m)` or `EllipticPi(n, φ, m)`
/// by a fixed composite Simpson rule on the defining integral.
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
    finite(composite_simpson(&f, 0.0, phi))
}

/// Panels of the fixed composite Simpson rule used by [`elliptic`].
///
/// The elliptic integrand is smooth on `|m| ≤ 0.9`, `|φ| ≤ π/2`; 2000 panels
/// put the quadrature error below 1e-12 (checked against the `K(1/2)` and
/// `E(π/2, 1/2)` references in the tests), which is what the 5-point central
/// difference needs. An *adaptive* rule chasing that accuracy instead cost
/// ~1 s per evaluation in debug builds — a single 14-integrand verification
/// test spent 320 s inside the quadrature, because `verify_antiderivative`
/// calls the oracle hundreds of times per case.
const ELLIPTIC_PANELS: usize = 2000;

/// Composite Simpson rule with a fixed even panel count.
fn composite_simpson(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let n = ELLIPTIC_PANELS;
    if a == b {
        return 0.0;
    }
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(a + i as f64 * h);
    }
    sum * h / 3.0
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
    fn oracle_is_indeterminate_for_unsupported_heads_and_constants() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // A head outside the oracle's table: honest indeterminate.
        let f = parse(&ctx, "Zeta(x)");
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
    fn special_function_table_matches_known_values() {
        // References evaluated with `mpmath` at 25 digits (the values SymPy
        // prints for the same heads; `Ci`/`Chi` on the negative axis use the
        // real even convention, whose real part is SymPy's principal value).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let cases: &[(&str, f64)] = &[
            ("Ei((2^-1))", 0.454_219_904_863_173_6),
            ("Ei(2)", 4.954_234_356_001_89),
            ("Ei(-1)", -0.219_383_934_395_520_3),
            ("Si(1)", 0.946_083_070_367_183),
            ("Si(-1)", -0.946_083_070_367_183),
            ("Si(5)", 1.549_931_244_944_674),
            ("Ci(1)", 0.337_403_922_900_968_13),
            ("Ci(5)", -0.190_029_749_656_643_87),
            ("Shi(2)", 2.501_567_433_354_975_6),
            ("Chi(2)", 2.452_666_922_646_914_5),
            ("fresnels(1)", 0.438_259_147_390_354_8),
            ("fresnelc(1)", 0.779_893_400_376_822_8),
            ("fresnels(5)", 0.499_191_381_917_116_9),
            ("fresnelc(5)", 0.563_631_188_704_012_2),
            // `Ei(n, z)` = `E_n(z)` = `sympy.expint(n, z)`, including the
            // elementary negative orders.
            ("Ei(1, 1)", 0.219_383_934_395_520_27),
            ("Ei(2, 2)", 0.037_534_261_820_490_45),
            ("Ei(0, 2)", 0.067_667_641_618_306_35),
            ("Ei(-1, 2)", 0.101_501_462_427_459_52),
            ("Ei(-2, 2)", 0.169_169_104_045_765_86),
            ("Ei(-3, 2)", 0.321_421_297_686_955_14),
        ];
        // Relative tolerance: the Taylor/asymptotic crossovers are ~1e-9 in
        // the worst case (`Si`/`Ci` near x = 20), everything else is ≤1e-12.
        for (src, want) in cases {
            let expr = parse(&ctx, src);
            match eval_f64(expr, &[]) {
                Eval::Value(v) => assert!(
                    (v - want).abs() <= 1e-8 * want.abs().max(1.0),
                    "{src}: got {v}, want {want}"
                ),
                other => panic!("{src}: expected a value, got {other:?}"),
            }
        }
    }

    #[test]
    fn special_function_crossovers_stay_accurate() {
        // The series/asymptotic crossovers are where accuracy is worst; check
        // both sides of both thresholds against the same references.
        assert!((si(20.0) - 1.548_241_701_043_439_8).abs() < 1e-8);
        assert!((si(24.0) - 1.554_738_691_722_919_1).abs() < 1e-9);
        assert!((ci(20.0) - 0.044_419_820_845_353_32).abs() < 1e-8);
        assert!((ci(24.0) - -0.038_333_015_551_247_15).abs() < 1e-9);
        assert!((fresnel_s(3.0) - 0.496_312_998_967_375).abs() < 1e-8);
        assert!((fresnel_s(4.0) - 0.420_515_754_246_928_4).abs() < 1e-9);
        assert!((fresnel_c(4.0) - 0.498_426_033_038_177_6).abs() < 1e-9);
        // The exponentially small negative-axis branch is where a naive
        // `γ + ln y + Σ` evaluation loses every significant digit:
        // Ei(−25) = −5.348899755340216640325e-13.
        let e = match ei(-25.0) {
            Eval::Value(v) => v,
            other => panic!("Ei(-25) should evaluate, got {other:?}"),
        };
        assert!(
            (e + 5.348_899_755_340_217e-13).abs() < 1e-20,
            "Ei(-25) = {e}"
        );
    }

    #[test]
    fn exponential_integral_asymptotic_branch_is_accurate() {
        fn value_of(e: Eval) -> f64 {
            match e {
                Eval::Value(v) => v,
                other => panic!("expected a value, got {other:?}"),
            }
        }
        // `|y| > 15` uses the asymptotic series rather than the series; check
        // both signs and the order-`n` recurrence there. References from
        // `mpmath` at 30 digits.
        let cases: &[(f64, f64)] = &[
            (value_of(ei(-20.0)), -9.835_525_290_649_882e-11),
            (value_of(ei(-30.0)), -3.021_552_010_688_812_5e-15),
            (shi(20.0), 12_807_826.332_028_294),
            (chi(20.0), 12_807_826.332_028_294),
        ];
        for (got, want) in cases {
            assert!(
                (got - want).abs() <= 1e-9 * want.abs().max(1.0),
                "got {got}, want {want}"
            );
        }
        // `Ei(n, z)` for `z > 15` via the recurrence from `E₁`.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for (src, want) in [
            ("Ei(1, 20)", 9.835_525_290_649_882e-11),
            ("Ei(2, 20)", 9.404_856_430_858_149e-11),
            ("Ei(3, 25)", 4.977_909_748_135_229e-13),
        ] {
            let expr = parse(&ctx, src);
            match eval_f64(expr, &[]) {
                Eval::Value(v) => assert!(
                    (v - want).abs() <= 1e-9 * want.abs().max(1.0),
                    "{src}: got {v}, want {want}"
                ),
                other => panic!("{src}: expected a value, got {other:?}"),
            }
        }
        // `Ei(−y) = −E₁(y)` must stay negative and tiny far out.
        let far = value_of(ei(-40.0));
        assert!(far < 0.0 && far.abs() < 1e-17, "Ei(-40) = {far}");
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
        let k = composite_simpson(
            &|t: f64| 1.0 / (1.0 - 0.5 * t.sin() * t.sin()).sqrt(),
            0.0,
            std::f64::consts::FRAC_PI_2,
        );
        assert!((k - 1.854_074_677_301_372).abs() < 1e-9);
        // E(π/2, 1/2) = 1.3506438810476755
        let e = composite_simpson(
            &|t: f64| (1.0 - 0.5 * t.sin() * t.sin()).sqrt(),
            0.0,
            std::f64::consts::FRAC_PI_2,
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

    #[test]
    fn dilogarithm_matches_mpmath_reference_grid() {
        // 40-digit references generated with mpmath (2026-11-02) from the
        // project venv, one value per grid point, e.g.
        //   ocas-tests/.venv/Scripts/python.exe -c "import mpmath; \
        //     mpmath.mp.dps=40; \
        //     print(mpmath.nstr(mpmath.polylog(2, mpmath.mpf('-0.9')), 40))"
        // These literals are deliberately NOT produced by the crate's own
        // series: the oracle must stay an independent reference (SPEC §3.4).
        const GRID: [(f64, f64); 8] = [
            (-0.9, -0.752_163_179_217_261_6),
            (-0.5, -0.448_414_206_923_646_2),
            (-0.1, -0.097_605_235_229_321_59),
            (0.1, 0.102_617_791_099_391_13),
            (0.3, 0.326_129_510_075_476_1),
            (0.5, 0.582_240_526_465_012_5),
            (0.7, 0.889_377_624_286_038_6),
            (0.9, 1.299_714_723_004_958_8),
        ];
        for (z, want) in GRID {
            match li2(z) {
                Eval::Value(v) => {
                    let rel = (v - want).abs() / want.abs();
                    assert!(
                        rel <= 1e-12,
                        "Li2({z}) = {v}, want {want}, relative error {rel:.3e}"
                    );
                }
                other => panic!("Li2({z}): expected a value, got {other:?}"),
            }
        }
        // The inversion branch (`z < −1`) needs its own independent probes;
        // same mpmath generation command, values at 40 digits:
        //   Li2(−1.2) = −0.9574053085587812190...
        //   Li2(−2)   = −1.4367463668836809463...
        //   Li2(−3)   = −1.9393754207667089530...
        const INVERTED: [(f64, f64); 3] = [
            (-1.2, -0.957_405_308_558_781_3),
            (-2.0, -1.436_746_366_883_680_8),
            (-3.0, -1.939_375_420_766_708_9),
        ];
        for (z, want) in INVERTED {
            match li2(z) {
                Eval::Value(v) => {
                    let rel = (v - want).abs() / want.abs();
                    assert!(
                        rel <= 1e-12,
                        "Li2({z}) = {v}, want {want}, relative error {rel:.3e}"
                    );
                }
                other => panic!("Li2({z}): expected a value, got {other:?}"),
            }
        }
    }

    #[test]
    fn dilogarithm_special_values() {
        // Li2(1) = pi^2/6 = 1.644934066848226436472415166646025189219,
        // Li2(0) = 0 and Li2(-1) = -pi^2/12
        // = -0.8224670334241132182362075833230125946095 (mpmath, 40 digits).
        // Evaluated through the atom path, so the dispatch is exercised too.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for (src, want) in [
            ("polylog(2, 1)", LI2_AT_ONE),
            ("polylog(2, 0)", 0.0),
            ("polylog(2, -1)", LI2_AT_MINUS_ONE),
        ] {
            match eval_f64(parse(&ctx, src), &[]) {
                Eval::Value(v) => assert!(
                    (v - want).abs() <= 1e-15 * want.abs().max(1.0),
                    "{src}: got {v}, want {want}"
                ),
                other => panic!("{src}: expected a value, got {other:?}"),
            }
        }
        // z = 1 and z = -1 hit the exact-constant branches, not a series sum:
        // pin the returned bits so a branch regression cannot hide behind the
        // tolerance above.
        assert_eq!(li2(1.0), Eval::Value(LI2_AT_ONE));
        assert_eq!(li2(0.0), Eval::Value(0.0));
        assert_eq!(li2(-1.0), Eval::Value(LI2_AT_MINUS_ONE));
        // The constants themselves agree with the 40-digit references.
        assert!((LI2_AT_ONE - 1.644_934_066_848_226_4).abs() < 1e-15);
        assert!((LI2_AT_MINUS_ONE + 0.822_467_033_424_113_2).abs() < 1e-15);
    }

    #[test]
    fn dilogarithm_domain_and_non_finite_inputs() {
        // z > 1 is outside the real domain (Li2 acquires an imaginary part).
        // It must be Domain, not Unsupported: Unsupported aborts the whole
        // verification round, while Domain only skips this sample.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // The parser has no decimal literals, so the out-of-domain probe is
        // written as the integer fraction 3/2 = 1.5.
        for src in ["polylog(2, 2)", "polylog(2, 3*(2^-1))"] {
            assert_eq!(
                eval_f64(parse(&ctx, src), &[]),
                Eval::Domain,
                "{src} should be out of the real domain"
            );
        }
        assert_eq!(li2(2.0), Eval::Domain);
        assert_eq!(li2(1.5), Eval::Domain);
        // Non-finite arguments are reported as non-finite values, not as a
        // domain error.
        for z in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(li2(z), Eval::NonFinite, "Li2({z}) should be NonFinite");
        }
    }

    #[test]
    fn dilogarithm_rejects_other_orders_and_arities() {
        // Only the real-valued order s = 2 is implemented in 0.31a; every
        // other order (and every arity except two) is an honest Unsupported
        // (SPEC §5: s ∉ {1, 2} is rejected at the head).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        for src in [
            "polylog(1, (2^-1))",
            "polylog(3, (2^-1))",
            "polylog((2^-1))",
            "polylog(2)",
            "polylog(2, (2^-1), 1)",
        ] {
            assert_eq!(
                eval_f64(parse(&ctx, src), &[]),
                Eval::Unsupported,
                "{src} should be Unsupported"
            );
        }
        // The order argument is evaluated rather than pattern-matched on the
        // atom: a symbolic order that evaluates to exactly 2 is accepted, one
        // that evaluates to any other value is not.
        let ok = parse(&ctx, "polylog(s, (2^-1))");
        assert!(matches!(
            eval_f64(ok, &[(Symbol::new("s"), 2.0)]),
            Eval::Value(_)
        ));
        assert_eq!(eval_f64(ok, &[(Symbol::new("s"), 3.0)]), Eval::Unsupported);
    }

    #[test]
    fn dilogarithm_derivative_matches_closed_form() {
        // d/dz Li2(z) = -ln(1-z)/z (SPEC §3.2). Central difference with
        // h = 1e-5: the truncation error is O(h^2), while the oracle's own
        // 5-point rule runs at h = 1e-4, too coarse for a 1e-9 target. z = 0
        // is excluded because the closed form divides by z; the grid's
        // closest point to z = 1 is 0.9, where the difference is still well
        // conditioned.
        const H: f64 = 1e-5;
        const GRID: [f64; 8] = [-0.9, -0.5, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9];
        for z in GRID {
            let at = |t: f64| match li2(t) {
                Eval::Value(v) => v,
                other => panic!("Li2({t}) did not evaluate: {other:?}"),
            };
            let numeric = (at(z + H) - at(z - H)) / (2.0 * H);
            let exact = -(1.0 - z).ln() / z;
            let rel = (numeric - exact).abs() / exact.abs();
            assert!(
                rel <= 1e-9,
                "d/dz Li2({z}) = {numeric}, want {exact}, relative error {rel:.3e}"
            );
        }
    }
}
