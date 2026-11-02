//! Polylogarithm heads and the dilogarithm (Li₂) reduction (0.31.0).
//!
//! `docs/planning/POLYLOG_SPEC_CN.md` is the specification of this wave; this
//! module is the 0.31a implementation of its two declared classes:
//!
//! - **A class** ([`integrate_polylog_head`]): an integrand that already
//!   carries `polylog(2, u)` with `u` affine in the integration variable is
//!   integrated by parts with the derivative table's
//!   `d/du Li₂(u) = −log(1−u)/u` (the same identity the reduction inverts,
//!   never a second copy of it). The residual is a log-rational integral and
//!   is handed to the B class.
//! - **B class** ([`integrate_log_rational`]): an elementary integrand
//!   `R(x)·log(A(x))` whose antiderivative needs Li₂. `R` is a rational
//!   weight with at most two linear denominator factors (and a numerator of
//!   lower degree), and `A` is a product/quotient/power of affine factors, so
//!   `log A` splits into a sum of `log(affine)` terms plus a constant. The
//!   substitution `u = A(x)` turns each term into one of the table entries
//!   below.
//!
//! # Identities (the whole reduction is a table, nothing is guessed)
//!
//! With `u` the (affine) logarithm argument and `γ`, `δ` constants:
//!
//! | integral | closed form |
//! |---|---|
//! | `∫ log(u)/u du` | `½·log²(u)` |
//! | `∫ log(u)/(u−β) du` | `log(u)·log(1−u/β) + Li₂(u/β)` |
//! | `∫ u^m·log(u) du`, `m ≥ 0` | `u^(m+1)·log(u)/(m+1) − u^(m+1)/(m+1)²` |
//! | `∫ 1/(γ+δu) du` | `log(γ+δu)/δ` |
//! | `1/((α₁+β₁x)(α₂+β₂x))` | `A/(α₁+β₁x) + B/(α₂+β₂x)` (exact two-factor split) |
//!
//! The pole entry is the SPEC's "pole term" (Landen form); the two-factor
//! split is the elementary identity
//! `N = A·l₂ + B·l₁` solved exactly. Sources were checked against SymPy
//! 1.14 (`sympy/integrals/rationaltools.py`, `prde.py`) and Lewin; the
//! per-entry record lives in the SPEC's implementation-deviation section.
//!
//! # Discipline
//!
//! - **Emission gate.** No candidate is returned before a local f64 check:
//!   the numerical derivative of the candidate is compared with the integrand
//!   at deterministic samples, over up to three parameter regimes. A shape
//!   with fewer than two usable samples, or with a relative disagreement
//!   above `GATE_TOL`, declines. The test oracle in `ocas-tests` is the
//!   independent check; this gate only stops unchecked emissions.
//! - **Refusals.** Non-affine logarithm arguments, products of logarithms,
//!   symbolic exponents, more than two denominator factors, improper
//!   fractions, complex coefficients, nested Li₂ and anything outside the
//!   fragment decline. The stage never emits an unevaluated residue and never
//!   commits to a path whose numerical evidence failed.
//! - **Budgets** are deterministic (`MAX_*` below); there is no wall-clock
//!   budget, so the stage cannot trade time for solves (POLYLOG SPEC §5).
//! - It declines inside a `chain::parts_active()` scope, like the other
//!   substitution engines (the Wave-4a parts discipline).

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use super::rules::rat_of;
use super::{
    chain, contains_integral, contains_symbol, is_constant, is_zero_expr, linear_form, node_count,
};

/// Maximum recursion depth of the term reduction (SPEC §5: 8).
const MAX_POLYLOG_STEPS: u32 = 8;
/// Maximum number of additive terms the folded integrand may have.
const MAX_TERMS: usize = 16;
/// Maximum number of logarithm factors one `log(…)` decomposes into.
const MAX_LOG_FACTORS: usize = 4;
/// Node ceiling for the input and for every derived form.
const MAX_POLYLOG_NODES: usize = 2000;
/// Node ceiling for an input that is handed to the (unbounded) product
/// expansion in `collect_terms`.
///
/// The B-class stage runs on every fallback of the pipeline, so a structural
/// product-over-sum expansion must be bounded before it is attempted; the
/// audited 0.31 targets are two orders of magnitude below this limit.
const MAX_FOLD_NODES: usize = 400;
/// Series term ceiling of the emission gate's Li₂ (SPEC §5: 4096).
const MAX_LI2_TERMS: usize = 4096;

/// Deterministic samples of the emission gate, in the integration variable.
///
/// The set is wider than the test oracle's because a Li₂ argument has to stay
/// on the real branch (`z ≤ 1`) for a sample to be usable: a correct answer
/// whose arguments leave the real domain on most samples still needs two
/// usable points to be accepted.
const GATE_SAMPLES: [f64; 11] = [
    0.15, 0.25, 0.35, 0.5, 0.6, 0.85, 1.2, 2.0, -0.2, -0.45, -0.8,
];
/// Relative tolerance of the emission gate.
///
/// The gate is a correctness guard, not a precision test: a wrong identity or
/// a lost sign is off by `O(1)`, while a correct one is at the finite-
/// difference noise floor (`~1e-10` for the shapes here). The bound is loose
/// enough to survive ill-conditioned samples and far below any real error.
const GATE_TOL: f64 = 1e-7;
/// Relative disagreement above which a candidate is rejected outright (as
/// opposed to "inconclusive, try another parameter regime").
const GROSS_TOL: f64 = 1e-3;

/// Diagnostic switch (`OCAS_POLYLOG_DEBUG=1`): print why the stage declined.
///
/// Same convention as `OCAS_INTEGRATE_TRACE`/`OCAS_CERT_TIMING`: a debugging
/// aid, not an API.
fn trace(message: &str) {
    if std::env::var_os("OCAS_POLYLOG_DEBUG").is_some() {
        eprintln!("[polylog] {message}");
    }
}

// =========================================================================
//  A/B class entry points
// =========================================================================

/// A class: `∫ R(x)·polylog(2, u(x)) dx`, `u` affine in `var`.
///
/// Called from the `special` stage (the polylog head makes the shape special
/// by definition, and no other stage claims it). Returns `None` outside the
/// fragment; never emits a residue.
pub(crate) fn integrate_polylog_head<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if chain::parts_active() || node_count(expr) > MAX_FOLD_NODES || !contains_polylog(expr) {
        return None;
    }
    let folded = crate::ode::util::collect_terms(ctx, expr);
    let terms = add_terms(folded);
    if terms.len() != 1 {
        return None;
    }
    let mut factors: Vec<Atom<'a>> = Vec::new();
    flat_factors(terms[0], &mut factors);
    let mut log_argument: Option<Atom<'a>> = None;
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for factor in factors {
        if let AtomNode::Fun(name, args) = factor.node() {
            if name.as_str() == "polylog" && args.len() == 2 {
                if log_argument.is_some() || !matches!(args[0].node(), AtomNode::Num(2)) {
                    return None;
                }
                log_argument = Some(args[1]);
                continue;
            }
            // Any other head (including a higher-order polylog) declines.
            return None;
        }
        rest.push(factor);
    }
    let raw_argument = match log_argument {
        Some(argument) => argument,
        None => {
            trace("A: no polylog(2, ·) factor");
            return None;
        }
    };
    // `linear_form` needs the expanded affine form: `c·(a+b·x)` is a product
    // until the product-over-sum expansion runs, and the folded argument is
    // the canonical one the residual and the emitted Li₂ must use.
    let u_arg = crate::ode::util::collect_terms(ctx, raw_argument);
    // `u` must really depend on the variable. `linear_form` returns the slope
    // and intercept; a zero slope is not a substitution.
    let Some((slope, intercept)) = linear_form(ctx, u_arg, var) else {
        trace("A: log argument is not affine");
        return None;
    };
    if is_zero_expr(ctx, slope) {
        trace("A: zero slope");
        return None;
    }
    let Some((weight_coef, weight_power)) = parse_monomial_weight(ctx, &rest, var) else {
        trace("A: weight outside the monomial fragment");
        return None;
    };
    let Some(primitive) = monomial_primitive(ctx, weight_coef, weight_power, var) else {
        trace("A: primal weight is not a supported monomial");
        return None;
    };
    let one_minus_u = ctx.add(&[ctx.num(1), ctx.mul(&[ctx.num(-1), u_arg])]);
    let logarithm = ctx.fun("log", &[one_minus_u]);
    // Residual of the parts split: `+∫ v·u'·log(1−u)/u dx` with `u' = slope`
    // and `v` the primitive of `R`.
    //
    // When `u` is a monomial (`u = u'·x`, zero intercept) the ratio `u'/u`
    // is exactly `1/x`; building it as `u'·u⁻¹` would leave an uncancelled
    // `c·x·(c·x)⁻¹` that the simplifier does not reduce and that the weight
    // parser would then misread as an improper fraction.
    let residual = if is_zero_expr(ctx, intercept) {
        let x = ctx.var(var.as_str());
        cleanup(
            ctx,
            normalize(ctx, ctx.mul(&[primitive, inv(ctx, x), logarithm])),
        )
    } else {
        let inverse_u = inv(ctx, u_arg);
        cleanup(
            ctx,
            normalize(ctx, ctx.mul(&[primitive, slope, inverse_u, logarithm])),
        )
    };
    trace(&format!("A: residual = {residual}"));
    let residual_anti = match integrate_log_rational(ctx, residual, var) {
        Some(anti) => anti,
        None => {
            trace("A: residual declined by the B class");
            return None;
        }
    };
    let li2 = ctx.fun("polylog", &[ctx.num(2), u_arg]);
    let value = cleanup(
        ctx,
        normalize(ctx, ctx.add(&[ctx.mul(&[li2, primitive]), residual_anti])),
    );
    if contains_integral(value) || node_count(value) > MAX_POLYLOG_NODES {
        trace("A: answer carries a residue or exceeds the node budget");
        return None;
    }
    if !gate_holds(ctx, expr, value, var) {
        trace(&format!("A: gate rejected {value}"));
        return None;
    }
    Some(value)
}

/// B class: log-rational integrals whose antiderivative needs Li₂.
///
/// Wired as the `log_rational` stage just before the heuristic stage, so it
/// only ever converts shapes every earlier stage declined. Returns `None`
/// outside the fragment; the emission gate must pass before anything is
/// returned.
pub(crate) fn integrate_log_rational<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    if chain::parts_active() || node_count(expr) > MAX_FOLD_NODES {
        return None;
    }
    if !has_log_head(expr) || !only_log_heads(expr) {
        trace("B: not a log-rational integrand");
        return None;
    }
    // Bound the product-over-sum expansion first (`collect_terms` expands
    // without a budget). `expand_bounded` returns `None` both for an
    // already-distributed input and for a budget abort; only the latter may
    // not be handed to the unbounded expansion.
    let expanded = match crate::expand::expand_bounded(ctx, expr) {
        Some(expanded) => expanded,
        None if needs_expansion(expr) => {
            trace("B: expansion budget exceeded");
            return None;
        }
        None => expr,
    };
    let folded = crate::ode::util::collect_terms(ctx, expanded);
    let aux = aux_symbol(expr, var)?;
    integrate_folded(ctx, expr, folded, var, aux, 0)
}

/// Integrate a folded (expanded, like-terms-combined) expression term by term.
fn integrate_folded<'a>(
    ctx: &'a AtomArena<'a>,
    original: Atom<'a>,
    folded: Atom<'a>,
    var: Symbol,
    aux: Symbol,
    step: u32,
) -> Option<Atom<'a>> {
    if step > MAX_POLYLOG_STEPS {
        return None;
    }
    let terms = add_terms(folded);
    if terms.is_empty() || terms.len() > MAX_TERMS {
        return None;
    }
    let mut total = ctx.num(0);
    for term in terms {
        total = ctx.add(&[total, integrate_term(ctx, term, var, aux, step)?]);
    }
    let total = cleanup(ctx, normalize(ctx, total));
    if contains_integral(total) || node_count(total) > MAX_POLYLOG_NODES {
        trace("B: answer carries a residue or exceeds the node budget");
        return None;
    }
    if !gate_holds(ctx, original, total, var) {
        trace(&format!("B: gate rejected {total}"));
        return None;
    }
    Some(total)
}

/// One additive term of the B class: `coef · weight(x) · log(A(x))`.
fn integrate_term<'a>(
    ctx: &'a AtomArena<'a>,
    term: Atom<'a>,
    var: Symbol,
    aux: Symbol,
    step: u32,
) -> Option<Atom<'a>> {
    if step > MAX_POLYLOG_STEPS {
        return None;
    }
    let mut factors: Vec<Atom<'a>> = Vec::new();
    flat_factors(term, &mut factors);
    let mut logs: Vec<Atom<'a>> = Vec::new();
    let mut rest: Vec<Atom<'a>> = Vec::new();
    for factor in factors {
        match factor.node() {
            AtomNode::Fun(name, args) if name.as_str() == "log" && args.len() == 1 => {
                logs.push(factor);
            }
            AtomNode::Fun(_, _) => return None,
            _ => rest.push(factor),
        }
    }
    if logs.len() > 1 {
        // Products of logarithms are outside the 0.31a fragment.
        return None;
    }
    let weight = parse_weight(ctx, &rest, var)?;
    if logs.is_empty() {
        return integrate_weight(ctx, &weight, var);
    }
    let AtomNode::Fun(_, args) = logs[0].node() else {
        return None;
    };
    let mut pairs: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    let mut constants: Vec<Atom<'a>> = Vec::new();
    accumulate_log(ctx, args[0], ctx.num(1), var, &mut pairs, &mut constants)?;
    if pairs.is_empty() && constants.is_empty() {
        return Some(ctx.num(0));
    }
    let mut total = ctx.num(0);
    for (coef, argument) in pairs {
        total = ctx.add(&[
            total,
            integrate_log_times_weight(ctx, coef, &weight, argument, var, aux, step + 1)?,
        ]);
    }
    if !constants.is_empty() {
        let mut factor = ctx.num(1);
        for c in constants {
            factor = ctx.mul(&[factor, c]);
        }
        let weight_anti = integrate_weight(ctx, &weight, var)?;
        total = ctx.add(&[total, ctx.mul(&[factor, weight_anti])]);
    }
    Some(normalize(ctx, total))
}

// =========================================================================
//  Log-argument decomposition
// =========================================================================

/// Split `log(arg)` into `Σ coef·log(affine)` pairs plus constant `log(·)`
/// factors. `log` of a product is a sum, `log` of an integer power is a
/// multiple, and a variable-free factor contributes a constant (the pipeline's
/// logarithm is the real `ln|·|`, so the split only differs by a piecewise
/// constant that the antiderivative's additive constant absorbs).
fn accumulate_log<'a>(
    ctx: &'a AtomArena<'a>,
    arg: Atom<'a>,
    coef: Atom<'a>,
    var: Symbol,
    pairs: &mut Vec<(Atom<'a>, Atom<'a>)>,
    constants: &mut Vec<Atom<'a>>,
) -> Option<()> {
    if pairs.len() > MAX_LOG_FACTORS || constants.len() > MAX_LOG_FACTORS {
        return None;
    }
    if is_constant(arg, var) {
        constants.push(ctx.mul(&[coef, ctx.fun("log", &[arg])]));
        return Some(());
    }
    match arg.node() {
        AtomNode::Mul(factors) => {
            for factor in factors.iter() {
                accumulate_log(ctx, *factor, coef, var, pairs, constants)?;
            }
            Some(())
        }
        AtomNode::Pow(base, exponent) => {
            let power = integer_exponent(*exponent)?;
            accumulate_log(
                ctx,
                *base,
                ctx.mul(&[coef, ctx.num(power)]),
                var,
                pairs,
                constants,
            )
        }
        _ => {
            // The argument must be affine in the variable: the substitution
            // `u = a + b·x` is what makes the table entries applicable. The
            // per-term expansion has not necessarily distributed products
            // inside a function argument (`1 + (−1)·(a·c + b·c·x)`), so the
            // candidate argument is folded first.
            let expanded = crate::ode::util::collect_terms(ctx, arg);
            linear_form(ctx, expanded, var)?;
            pairs.push((coef, expanded));
            Some(())
        }
    }
}

// =========================================================================
//  Weight parsing (rational factor structure)
// =========================================================================

/// `alpha + beta·x` raised to `power`; `beta` is not identically zero.
#[derive(Clone, Copy)]
struct LinearFactor<'a> {
    alpha: Atom<'a>,
    beta: Atom<'a>,
    power: i64,
}

/// `coef · Π (alpha_i + beta_i·x)^(power_i)`.
struct Weight<'a> {
    coef: Atom<'a>,
    factors: Vec<LinearFactor<'a>>,
}

/// Parse a list of multiplicative factors as a rational weight.
///
/// Fragment limits (v1): at most two denominator factors, each with exponent
/// `−1`; the numerator degree must be lower than the denominator degree. That
/// is exactly what the 0.31a targets need, and it keeps the reduction free of
/// polynomial division.
fn parse_weight<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    var: Symbol,
) -> Option<Weight<'a>> {
    let mut coef = ctx.num(1);
    let mut linear: Vec<LinearFactor<'a>> = Vec::new();
    for &factor in factors {
        if is_constant(factor, var) {
            coef = ctx.mul(&[coef, factor]);
            continue;
        }
        linear.push(parse_linear_factor(ctx, factor, var)?);
    }
    let mut numerator_degree = 0i64;
    let mut denominator_factors = 0i64;
    for factor in &linear {
        if factor.power == 0 {
            continue;
        }
        if factor.power < 0 {
            if factor.power != -1 {
                // Powers of a pole need a by-parts recursion; out of scope.
                return None;
            }
            denominator_factors += 1;
        } else {
            numerator_degree += factor.power;
        }
    }
    if denominator_factors > 2 {
        return None;
    }
    if denominator_factors == 0 {
        if numerator_degree > 2 {
            return None;
        }
    } else if numerator_degree >= denominator_factors {
        // An improper fraction needs polynomial division first.
        return None;
    }
    Some(Weight {
        coef,
        factors: linear,
    })
}

/// Parse one factor as `(alpha + beta·x)^power`.
fn parse_linear_factor<'a>(
    ctx: &'a AtomArena<'a>,
    factor: Atom<'a>,
    var: Symbol,
) -> Option<LinearFactor<'a>> {
    if let AtomNode::Pow(base, exponent) = factor.node() {
        let power = integer_exponent(*exponent)?;
        let (beta, alpha) = linear_form(ctx, *base, var)?;
        if is_zero_expr(ctx, beta) {
            return None;
        }
        return Some(LinearFactor { alpha, beta, power });
    }
    let (beta, alpha) = linear_form(ctx, factor, var)?;
    if is_zero_expr(ctx, beta) {
        return None;
    }
    Some(LinearFactor {
        alpha,
        beta,
        power: 1,
    })
}

/// `∫ weight(x) dx` for the parsed fragment (Pareto: rational part only).
fn integrate_weight<'a>(
    ctx: &'a AtomArena<'a>,
    weight: &Weight<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    let mut numerator: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    let mut denominator: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    for factor in &weight.factors {
        if factor.power == 0 {
            continue;
        }
        if factor.power < 0 {
            denominator.push((factor.alpha, factor.beta));
        } else {
            for _ in 0..factor.power {
                numerator.push((factor.alpha, factor.beta));
            }
        }
    }
    let mut polynomial: Vec<Atom<'a>> = vec![weight.coef];
    for (alpha, beta) in &numerator {
        polynomial = mul_by_linear(ctx, &polynomial, *alpha, *beta, x)?;
    }
    match denominator.len() {
        0 => {
            let mut total = ctx.num(0);
            for (j, coefficient) in polynomial.iter().enumerate() {
                let power = ctx.num(j as i64 + 1);
                total = ctx.add(&[
                    total,
                    ctx.mul(&[*coefficient, ctx.pow(x, power), ctx.pow(power, ctx.num(-1))]),
                ]);
            }
            Some(total)
        }
        1 => {
            if polynomial.len() != 1 {
                return None;
            }
            let (alpha, beta) = denominator[0];
            let logarithm = ctx.fun("log", &[ctx.add(&[alpha, ctx.mul(&[beta, x])])]);
            Some(ctx.mul(&[polynomial[0], ctx.pow(beta, ctx.num(-1)), logarithm]))
        }
        _ => {
            if polynomial.len() > 2 {
                return None;
            }
            let n0 = polynomial[0];
            let n1 = polynomial.get(1).copied().unwrap_or_else(|| ctx.num(0));
            let (a, b) = split_two_linear(ctx, n0, n1, denominator[0], denominator[1])?;
            let l1 = ctx.fun(
                "log",
                &[ctx.add(&[denominator[0].0, ctx.mul(&[denominator[0].1, x])])],
            );
            let l2 = ctx.fun(
                "log",
                &[ctx.add(&[denominator[1].0, ctx.mul(&[denominator[1].1, x])])],
            );
            Some(ctx.add(&[
                ctx.mul(&[a, ctx.pow(denominator[0].1, ctx.num(-1)), l1]),
                ctx.mul(&[b, ctx.pow(denominator[1].1, ctx.num(-1)), l2]),
            ]))
        }
    }
}

/// `∫ coef·weight(x)·log(arg) dx` for one affine `arg`.
fn integrate_log_times_weight<'a>(
    ctx: &'a AtomArena<'a>,
    coef: Atom<'a>,
    weight: &Weight<'a>,
    arg: Atom<'a>,
    var: Symbol,
    aux: Symbol,
    step: u32,
) -> Option<Atom<'a>> {
    if step > MAX_POLYLOG_STEPS {
        return None;
    }
    let (slope, intercept) = linear_form(ctx, arg, var)?;
    if is_zero_expr(ctx, slope) {
        return None;
    }
    let slope_inverse = ctx.pow(slope, ctx.num(-1));
    let u = ctx.var(aux.as_str());
    // x = (u − intercept)/slope. Every x-linear factor becomes u-linear.
    let mut denominator: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    let mut numerator: Vec<(Atom<'a>, Atom<'a>)> = Vec::new();
    for factor in &weight.factors {
        if factor.power == 0 {
            continue;
        }
        let gamma = ctx.add(&[
            factor.alpha,
            ctx.mul(&[ctx.num(-1), factor.beta, intercept, slope_inverse]),
        ]);
        let delta = ctx.mul(&[factor.beta, slope_inverse]);
        if is_zero_expr(ctx, delta) {
            return None;
        }
        if factor.power < 0 {
            denominator.push((gamma, delta));
        } else {
            for _ in 0..factor.power {
                numerator.push((gamma, delta));
            }
        }
    }
    // `weight(x)·dx = coef·Π factors · du/slope`: the 1/slope enters here.
    let mut polynomial: Vec<Atom<'a>> = vec![ctx.mul(&[coef, weight.coef, slope_inverse])];
    for (gamma, delta) in &numerator {
        polynomial = mul_by_linear(ctx, &polynomial, *gamma, *delta, u)?;
    }
    let in_u = match denominator.len() {
        0 => poly_log_integral(ctx, &polynomial, u)?,
        1 => {
            if polynomial.len() != 1 {
                return None;
            }
            pole_integral(ctx, polynomial[0], denominator[0].0, denominator[0].1, u)?
        }
        _ => {
            if polynomial.len() > 2 {
                return None;
            }
            let n0 = polynomial[0];
            let n1 = polynomial.get(1).copied().unwrap_or_else(|| ctx.num(0));
            let (a, b) = split_two_linear(ctx, n0, n1, denominator[0], denominator[1])?;
            ctx.add(&[
                pole_integral(ctx, a, denominator[0].0, denominator[0].1, u)?,
                pole_integral(ctx, b, denominator[1].0, denominator[1].1, u)?,
            ])
        }
    };
    // Back-substitute `u = arg` and give the caller the x-expression.
    let substituted = super::replace_symbol(ctx, in_u, aux, arg);
    Some(normalize(ctx, substituted))
}

/// `∫ P(u)·log(u) du` for a polynomial of degree at most two.
fn poly_log_integral<'a>(
    ctx: &'a AtomArena<'a>,
    polynomial: &[Atom<'a>],
    u: Atom<'a>,
) -> Option<Atom<'a>> {
    if polynomial.len() > 3 {
        return None;
    }
    let logarithm = ctx.fun("log", &[u]);
    let mut total = ctx.num(0);
    for (j, coefficient) in polynomial.iter().enumerate() {
        let power = j as i64 + 1;
        let upow = ctx.pow(u, ctx.num(power));
        let first = ctx.mul(&[
            *coefficient,
            upow,
            ctx.pow(ctx.num(power), ctx.num(-1)),
            logarithm,
        ]);
        let second = ctx.mul(&[
            *coefficient,
            upow,
            ctx.pow(ctx.num(power * power), ctx.num(-1)),
        ]);
        total = ctx.add(&[total, first, ctx.mul(&[ctx.num(-1), second])]);
    }
    Some(normalize(ctx, total))
}

/// `∫ coef·log(u)/(γ + δ·u) du`.
fn pole_integral<'a>(
    ctx: &'a AtomArena<'a>,
    coef: Atom<'a>,
    gamma: Atom<'a>,
    delta: Atom<'a>,
    u: Atom<'a>,
) -> Option<Atom<'a>> {
    let log_u = ctx.fun("log", &[u]);
    if is_zero_expr(ctx, gamma) {
        // The pole sits at the origin: `∫ log(u)/u du = ½·log²(u)`.
        return Some(ctx.mul(&[coef, rat_atom(ctx, 1, 2), log_u, log_u]));
    }
    let ratio = ctx.mul(&[delta, u, ctx.pow(gamma, ctx.num(-1))]);
    let log_one_plus = ctx.fun("log", &[ctx.add(&[ctx.num(1), ratio])]);
    let li2 = ctx.fun("polylog", &[ctx.num(2), ctx.mul(&[ctx.num(-1), ratio])]);
    let core = ctx.add(&[ctx.mul(&[log_u, log_one_plus]), li2]);
    Some(normalize(
        ctx,
        ctx.mul(&[coef, ctx.pow(delta, ctx.num(-1)), core]),
    ))
}

/// Split `(n0 + n1·t)/(l1·l2)` into `A/l1 + B/l2`, both `li = αi + βi·t`.
fn split_two_linear<'a>(
    ctx: &'a AtomArena<'a>,
    n0: Atom<'a>,
    n1: Atom<'a>,
    l1: (Atom<'a>, Atom<'a>),
    l2: (Atom<'a>, Atom<'a>),
) -> Option<(Atom<'a>, Atom<'a>)> {
    // N = A·l2 + B·l1  ⇒  A = (n1·α1 − n0·β1)/Δ, B = (β2·n0 − α2·n1)/Δ,
    // Δ = β2·α1 − β1·α2.
    let det = ctx.add(&[ctx.mul(&[l2.1, l1.0]), ctx.mul(&[ctx.num(-1), l1.1, l2.0])]);
    if is_zero_expr(ctx, det) {
        return None;
    }
    let det_inverse = ctx.pow(det, ctx.num(-1));
    let a_numerator = ctx.add(&[ctx.mul(&[n1, l1.0]), ctx.mul(&[ctx.num(-1), n0, l1.1])]);
    let b_numerator = ctx.add(&[ctx.mul(&[l2.1, n0]), ctx.mul(&[ctx.num(-1), l2.0, n1])]);
    Some((
        ctx.mul(&[a_numerator, det_inverse]),
        ctx.mul(&[b_numerator, det_inverse]),
    ))
}

/// Multiply a polynomial (in some variable `t`) by `c0 + c1·t`.
fn mul_by_linear<'a>(
    ctx: &'a AtomArena<'a>,
    polynomial: &[Atom<'a>],
    c0: Atom<'a>,
    c1: Atom<'a>,
    _t: Atom<'a>,
) -> Option<Vec<Atom<'a>>> {
    if polynomial.len() >= 3 {
        return None;
    }
    let mut out: Vec<Atom<'a>> = vec![ctx.num(0); polynomial.len() + 1];
    for (j, coefficient) in polynomial.iter().enumerate() {
        out[j] = ctx.add(&[out[j], ctx.mul(&[*coefficient, c0])]);
        out[j + 1] = ctx.add(&[out[j + 1], ctx.mul(&[*coefficient, c1])]);
    }
    Some(out)
}

/// The primitive of a monomial weight `c·x^k` (A class): `k ∈ {−2, 0, 1}`.
///
/// The by-parts residual must stay inside the B class, so the primitive is
/// restricted to monomials: a logarithmic primitive would turn the residual
/// into a product of two logarithms (out of the 0.31a fragment).
fn monomial_primitive<'a>(
    ctx: &'a AtomArena<'a>,
    coef: Atom<'a>,
    power: i64,
    var: Symbol,
) -> Option<Atom<'a>> {
    let x = ctx.var(var.as_str());
    match power {
        0 => Some(ctx.mul(&[coef, x])),
        1 => Some(ctx.mul(&[coef, rat_atom(ctx, 1, 2), ctx.pow(x, ctx.num(2))])),
        -2 => Some(ctx.mul(&[coef, ctx.num(-1), ctx.pow(x, ctx.num(-1))])),
        _ => None,
    }
}

/// Parse the A-class's rational weight as a monomial `coef·x^k`.
///
/// The B-class [`parse_weight`] enforces the proper-fraction fragment; the A
/// class only needs the weight to have a monomial primitive, which allows
/// `x^k` with `k ∈ {−2, 0, 1}` (the by-parts residual must stay reducible).
fn parse_monomial_weight<'a>(
    ctx: &'a AtomArena<'a>,
    factors: &[Atom<'a>],
    var: Symbol,
) -> Option<(Atom<'a>, i64)> {
    let mut coef = ctx.num(1);
    let mut power = 0i64;
    for &factor in factors {
        if is_constant(factor, var) {
            coef = ctx.mul(&[coef, factor]);
            continue;
        }
        let linear = parse_linear_factor(ctx, factor, var)?;
        if !is_zero_expr(ctx, linear.alpha) {
            // Not a monomial in the variable.
            return None;
        }
        power = power.checked_add(linear.power)?;
    }
    Some((coef, power))
}

// =========================================================================
//  Numeric emission gate
// =========================================================================

/// Check the candidate by numerically differentiating it against the
/// integrand at deterministic samples, over up to three parameter regimes.
///
/// Returns `true` only when at least two samples are usable in some regime and
/// the worst relative disagreement there is within [`GATE_TOL`].
fn gate_holds<'a>(
    ctx: &'a AtomArena<'a>,
    integrand: Atom<'a>,
    candidate: Atom<'a>,
    var: Symbol,
) -> bool {
    let mut parameters: Vec<Symbol> = Vec::new();
    collect_parameters(integrand, var, &mut parameters);
    collect_parameters(candidate, var, &mut parameters);
    let regimes = if parameters.is_empty() { 1 } else { 3 };
    let mut polylog_arguments: Vec<Atom<'a>> = Vec::new();
    collect_polylog_arguments(candidate, &mut polylog_arguments);
    for regime in 0..regimes {
        let environment: Vec<(Symbol, f64)> = parameters
            .iter()
            .enumerate()
            .map(|(index, symbol)| (*symbol, gate_parameter(index, regime)))
            .collect();
        // Fixed samples plus samples that place each `polylog(2, affine)`
        // argument inside the real branch (`u ≈ ½`). A Li₂ antiderivative
        // usually carries both `Li₂(u)` and `Li₂(1−u)`, so its real domain is
        // the narrow `u ∈ [0, 1]` window; a fixed grid can easily miss it.
        let mut samples: Vec<f64> = GATE_SAMPLES.to_vec();
        for argument in &polylog_arguments {
            let Some((slope, intercept)) = linear_form(ctx, *argument, var) else {
                continue;
            };
            let Some(slope_value) = eval_real(slope, var, 0.0, &environment) else {
                continue;
            };
            let Some(intercept_value) = eval_real(intercept, var, 0.0, &environment) else {
                continue;
            };
            if slope_value == 0.0 {
                continue;
            }
            for target in [0.5, -0.5] {
                let x = (target - intercept_value) / slope_value;
                if x.is_finite() && x.abs() <= 50.0 && x.abs() >= 0.05 {
                    samples.push(x);
                }
            }
        }
        let mut checked = 0usize;
        let mut worst = 0.0f64;
        for &x in samples.iter() {
            let target = eval_real(integrand, var, x, &environment);
            let Some(target_value) = target else {
                if std::env::var_os("OCAS_POLYLOG_DEBUG").is_some() {
                    eprintln!("[polylog] gate regime {regime} x={x}: integrand not evaluable");
                }
                continue;
            };
            if !target_value.is_finite() {
                continue;
            }
            let h = 1e-5 * x.abs().max(1.0);
            let mut values = [0.0f64; 4];
            let mut usable = true;
            for (index, offset) in [-2.0f64, -1.0, 1.0, 2.0].iter().enumerate() {
                match eval_real(candidate, var, x + offset * h, &environment) {
                    Some(value) if value.is_finite() => values[index] = value,
                    _ => {
                        if std::env::var_os("OCAS_POLYLOG_DEBUG").is_some() {
                            eprintln!(
                                "[polylog] gate regime {regime} x={x}: candidate not evaluable at offset {offset}"
                            );
                        }
                        usable = false;
                        break;
                    }
                }
            }
            if !usable {
                continue;
            }
            let derivative =
                (-values[3] + 8.0 * values[2] - 8.0 * values[1] + values[0]) / (12.0 * h);
            if !derivative.is_finite() {
                continue;
            }
            checked += 1;
            let relative = (derivative - target_value).abs() / target_value.abs().max(1.0);
            if relative > worst {
                worst = relative;
            }
            if std::env::var_os("OCAS_POLYLOG_DEBUG").is_some() {
                eprintln!(
                    "[polylog] gate regime {regime} x={x}: d={derivative} f={target_value} rel={relative:e}"
                );
            }
        }
        if checked < 2 {
            // Not enough real-domain samples to decide in this regime; the
            // next parameter assignment may be better conditioned.
            continue;
        }
        if worst > GROSS_TOL {
            // A usable regime that disagrees by orders of magnitude is
            // decisive evidence against the candidate (a wrong identity, a
            // lost sign or a lost constant); parameter values cannot repair
            // it.
            return false;
        }
        if worst <= GATE_TOL {
            return true;
        }
        // Inconclusive in this regime (finite-difference conditioning): try
        // the next one rather than rejecting a candidate that only failed a
        // precision check.
    }
    false
}

/// Deterministic parameter values for the gate's three regimes.
fn gate_parameter(index: usize, regime: usize) -> f64 {
    const MIXED: [f64; 8] = [2.0, -0.75, 3.0, 0.5, -1.25, 4.0, 0.25, 1.5];
    const POSITIVE_A: [f64; 8] = [2.0, 0.75, 3.0, 0.5, 1.25, 4.0, 0.25, 1.5];
    const POSITIVE_B: [f64; 8] = [1.0, 2.5, 0.5, 3.0, 4.0, 0.8, 1.5, 2.0];
    match regime {
        0 => MIXED[index % MIXED.len()],
        1 => POSITIVE_A[index % POSITIVE_A.len()],
        _ => POSITIVE_B[index % POSITIVE_B.len()],
    }
}

/// Free parameters of `expr`, in first-seen order (deterministic).
///
/// Only `pi` is excluded: the pipeline emits it as the mathematical constant
/// (`polylog(2,1) → π²/6`, `atan` coefficients), so sampling it as a free
/// parameter would scale a correct answer and wrongly reject it. `e`/`E` are
/// **parameters**: the Rubi corpus uses `e` as a free symbol and the parser
/// never turns it into Euler's number, so sampling it exercises the identity
/// at several values instead of fixing it to one constant. (The test oracle
/// keeps its historical constant-name convention; the gate is an independent
/// check and may sample more.)
fn collect_parameters(expr: Atom<'_>, var: Symbol, out: &mut Vec<Symbol>) {
    match expr.node() {
        AtomNode::Num(_) => {}
        AtomNode::Var(symbol) => {
            if *symbol != var && symbol.as_str() != "pi" && !out.contains(symbol) {
                out.push(*symbol);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for arg in args.iter() {
                collect_parameters(*arg, var, out);
            }
        }
        AtomNode::Pow(base, exponent) => {
            collect_parameters(*base, var, out);
            collect_parameters(*exponent, var, out);
        }
        AtomNode::Fun(_, args) => {
            for arg in args.iter() {
                collect_parameters(*arg, var, out);
            }
        }
    }
}

/// Collect the argument of every `polylog(2, ·)` atom in `expr`.
///
/// Used by the emission gate to place extra samples where the dilogarithm
/// argument is inside the real branch: an A-class answer carries `Li₂(u)` and
/// `Li₂(1−u)`, whose common real domain is the narrow `u ∈ [0, 1]` window.
fn collect_polylog_arguments<'a>(expr: Atom<'a>, out: &mut Vec<Atom<'a>>) {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if name.as_str() == "polylog"
                && args.len() == 2
                && matches!(args[0].node(), AtomNode::Num(2))
            {
                out.push(args[1]);
            }
            for arg in args.iter() {
                collect_polylog_arguments(*arg, out);
            }
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            for arg in args.iter() {
                collect_polylog_arguments(*arg, out);
            }
        }
        AtomNode::Pow(base, exponent) => {
            collect_polylog_arguments(*base, out);
            collect_polylog_arguments(*exponent, out);
        }
        AtomNode::Num(_) | AtomNode::Var(_) => {}
    }
}

/// Real-valued f64 evaluation used only by the emission gate.
///
/// `log` is `ln|v|` (the pipeline's real convention) and `polylog(2, z)` is
/// real only for `z ≤ 1`; anything outside the real domain returns `None` and
/// the sample is skipped.
fn eval_real(expr: Atom<'_>, var: Symbol, x: f64, environment: &[(Symbol, f64)]) -> Option<f64> {
    match expr.node() {
        AtomNode::Num(n) => Some(*n as f64),
        AtomNode::Var(symbol) => {
            if *symbol == var {
                Some(x)
            } else if let Some((_, value)) = environment.iter().find(|(name, _)| name == symbol) {
                // A parameter named `e` (very common in the Rubi corpus) must
                // shadow Euler's number, exactly like the test oracle's
                // `eval_f64` resolves its environment first.
                Some(*value)
            } else {
                match symbol.as_str() {
                    "pi" => Some(std::f64::consts::PI),
                    "e" | "E" => Some(std::f64::consts::E),
                    _ => None,
                }
            }
        }
        AtomNode::Add(args) => {
            let mut total = 0.0;
            for arg in args.iter() {
                total += eval_real(*arg, var, x, environment)?;
            }
            Some(total)
        }
        AtomNode::Mul(args) => {
            let mut total = 1.0;
            for arg in args.iter() {
                total *= eval_real(*arg, var, x, environment)?;
            }
            Some(total)
        }
        AtomNode::Pow(base, exponent) => {
            let base = eval_real(*base, var, x, environment)?;
            let exponent = eval_real(*exponent, var, x, environment)?;
            if base < 0.0 && exponent.fract() != 0.0 {
                return None;
            }
            Some(base.powf(exponent))
        }
        AtomNode::Fun(name, args) => match name.as_str() {
            "log" if args.len() == 1 => {
                let value = eval_real(args[0], var, x, environment)?;
                if value == 0.0 {
                    None
                } else {
                    Some(value.abs().ln())
                }
            }
            "sqrt" if args.len() == 1 => {
                let value = eval_real(args[0], var, x, environment)?;
                if value < 0.0 {
                    None
                } else {
                    Some(value.sqrt())
                }
            }
            "polylog" if args.len() == 2 => {
                let order = eval_real(args[0], var, x, environment)?;
                if order != 2.0 {
                    return None;
                }
                li2_real(eval_real(args[1], var, x, environment)?)
            }
            _ => None,
        },
    }
}

/// Real dilogarithm for the gate (series / reflection / inversion).
fn li2_real(z: f64) -> Option<f64> {
    if !z.is_finite() {
        return None;
    }
    if z > 1.0 {
        return None;
    }
    if z == 1.0 {
        return Some(std::f64::consts::PI * std::f64::consts::PI / 6.0);
    }
    if z.abs() <= 0.5 {
        return Some(li2_series(z));
    }
    if z > 0.5 {
        // Li₂(z) = π²/6 − log(z)·log(1−z) − Li₂(1−z), with 1−z < 0.5.
        let reflected = 1.0 - z;
        let pi2_over_six = std::f64::consts::PI * std::f64::consts::PI / 6.0;
        return Some(pi2_over_six - z.ln() * reflected.ln() - li2_series(reflected));
    }
    if z >= -2.0 {
        // −2 ≤ z < −0.5: Landen's identity
        // Li₂(z) = −Li₂(z/(z−1)) − ½·log²(1−z), whose argument lies in
        // (1/3, 2/3) — the Maclaurin series is accurate there for both the
        // value and its derivative. (The plain series at z ≈ −1 converges too
        // slowly to differentiate through, which is what the gate does.)
        let mapped = z / (z - 1.0);
        return Some(-li2_series(mapped) - 0.5 * (1.0 - z).ln().powi(2));
    }
    // z < −2: Li₂(z) = −π²/6 − ½·log²(−z) − Li₂(1/z), with 1/z ∈ (−0.5, 0).
    let inverse = 1.0 / z;
    let pi2_over_six = std::f64::consts::PI * std::f64::consts::PI / 6.0;
    Some(-pi2_over_six - 0.5 * (-z).ln().powi(2) - li2_series(inverse))
}

/// Maclaurin series `Σ z^k/k²`, used for `|z| ≤ 0.5`.
fn li2_series(z: f64) -> f64 {
    let mut total = 0.0;
    let mut power = 1.0;
    for k in 1..=MAX_LI2_TERMS {
        power *= z;
        let term = power / (k as f64 * k as f64);
        total += term;
        if term.abs() < 1e-18 * total.abs().max(1e-300) {
            break;
        }
    }
    total
}

// =========================================================================
//  Small helpers
// =========================================================================

/// Whether the expression contains a product over a sum (or a small power of a
/// sum) that a distributive expansion would rewrite.
///
/// Mirrors the rewrite rules of [`crate::expand::expand_bounded`] exactly —
/// including its treatment of function applications as opaque atoms, since
/// the bounded expansion does not descend into function arguments. The
/// polynomial parts of a logarithm argument are expanded separately by
/// [`accumulate_log`], after the argument has been isolated.
fn needs_expansion(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            if args.iter().any(|arg| needs_expansion(*arg)) {
                return true;
            }
            matches!(expr.node(), AtomNode::Mul(_))
                && args
                    .iter()
                    .any(|arg| matches!(arg.node(), AtomNode::Add(_)))
        }
        AtomNode::Pow(base, exponent) => {
            needs_expansion(*base)
                || needs_expansion(*exponent)
                || (matches!(base.node(), AtomNode::Add(_))
                    && matches!(exponent.node(), AtomNode::Num(n) if (2..=8).contains(n)))
        }
        AtomNode::Fun(_, _) | AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Collect the flat factor list of a product, descending into nested `Mul`s.
fn flat_factors<'a>(expr: Atom<'a>, out: &mut Vec<Atom<'a>>) {
    match expr.node() {
        AtomNode::Mul(args) => {
            for arg in args.iter() {
                flat_factors(*arg, out);
            }
        }
        _ => out.push(expr),
    }
}

/// The additive terms of an expression.
fn add_terms<'a>(expr: Atom<'a>) -> Vec<Atom<'a>> {
    match expr.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![expr],
    }
}

/// Whether the expression contains a `log` head anywhere.
fn has_log_head(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            (name.as_str() == "log" && args.len() == 1) || args.iter().any(|arg| has_log_head(*arg))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|arg| has_log_head(*arg)),
        AtomNode::Pow(base, exponent) => has_log_head(*base) || has_log_head(*exponent),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Whether every function head of the expression is `log` (no `polylog`,
/// `sqrt`, `exp`, …): the B class is a log-rational fragment only.
fn only_log_heads(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            name.as_str() == "log" && args.len() == 1 && args.iter().all(|arg| only_log_heads(*arg))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().all(|arg| only_log_heads(*arg)),
        AtomNode::Pow(base, exponent) => only_log_heads(*base) && only_log_heads(*exponent),
        AtomNode::Num(_) | AtomNode::Var(_) => true,
    }
}

/// Whether the expression contains a `polylog` head.
fn contains_polylog(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            name.as_str() == "polylog" || args.iter().any(|arg| contains_polylog(*arg))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|arg| contains_polylog(*arg)),
        AtomNode::Pow(base, exponent) => contains_polylog(*base) || contains_polylog(*exponent),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// An integer exponent atom, or `None` for a non-integer one.
fn integer_exponent(exponent: Atom<'_>) -> Option<i64> {
    let (numerator, denominator) = rat_of(exponent)?;
    if denominator == 1 {
        Some(numerator)
    } else {
        None
    }
}

/// `1/a`, keeping integer reciprocals exact.
fn inv<'a>(ctx: &'a AtomArena<'a>, atom: Atom<'a>) -> Atom<'a> {
    match atom.node() {
        AtomNode::Num(1) => ctx.num(1),
        AtomNode::Num(-1) => ctx.num(-1),
        _ => ctx.pow(atom, ctx.num(-1)),
    }
}

/// `p/q` as an atom.
fn rat_atom<'a>(ctx: &'a AtomArena<'a>, numerator: i64, denominator: i64) -> Atom<'a> {
    if numerator == 0 {
        return ctx.num(0);
    }
    if denominator == 1 {
        return ctx.num(numerator);
    }
    ctx.mul(&[
        ctx.num(numerator),
        ctx.pow(ctx.num(denominator), ctx.num(-1)),
    ])
}

/// Clean a derived form with the same rule sets the pipeline uses for every
/// candidate (`default_rules` then `calculus_rules`).
///
/// The A class's by-parts residual is built as `v·u'·u⁻¹·log(1−u)`, and
/// `normalize` does not cancel the reciprocal against the factors of `u`
/// (e.g. `c·x·(c·x)⁻¹` stays unevaluated). Without this cleanup the residual
/// would misparse as an improper fraction and the reduction would decline a
/// shape it can actually handle.
fn cleanup<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    let algebraic = ocas_rewrite::rules::default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let after_algebraic = ocas_rewrite::simplify::simplify(ctx, expr, &algebraic, 20);
    let calculus = crate::rules::calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    ocas_rewrite::simplify::simplify(ctx, after_algebraic, &calculus, 10)
}

/// Pick a substitution symbol that collides with neither `var` nor any symbol
/// already present in `expr`.
fn aux_symbol(expr: Atom<'_>, var: Symbol) -> Option<Symbol> {
    for name in ["u", "t", "w", "v"] {
        let symbol = Symbol::new(name);
        if symbol != var && !contains_symbol(expr, symbol) {
            return Some(symbol);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use ocas_core::arena::Arena;

    use super::*;

    fn parse<'a>(ctx: &'a AtomArena<'a>, source: &str) -> Atom<'a> {
        ocas_parse::parse(ctx, source).unwrap_or_else(|error| panic!("parse {source}: {error:?}"))
    }

    fn environment(symbols: &[Symbol], regime: usize) -> Vec<(Symbol, f64)> {
        symbols
            .iter()
            .enumerate()
            .map(|(index, symbol)| (*symbol, gate_parameter(index, regime)))
            .collect()
    }

    /// Numeric check of `d/dx(antiderivative) = integrand` at `samples`,
    /// with the module's own real evaluator (the test oracle in `ocas-tests`
    /// is the independent check).
    fn check_numeric<'a>(
        _ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        antiderivative: Atom<'a>,
        var: Symbol,
        _symbols: &[Symbol],
        samples: &[f64],
    ) {
        // Derive the parameter list exactly like the emission gate does: the
        // assignment of values is positional, so a different order would test
        // a different parameter assignment (and could push every Li₂ argument
        // outside the real domain).
        let mut parameters: Vec<Symbol> = Vec::new();
        collect_parameters(integrand, var, &mut parameters);
        collect_parameters(antiderivative, var, &mut parameters);
        // Try the same parameter regimes the emission gate uses: a parameter
        // assignment can push a Li₂ argument outside the real domain, and the
        // test must not mistake that for a wrong answer.
        let mut best: Option<(usize, f64)> = None;
        for regime in 0..3 {
            let ready = environment(&parameters, regime);
            let mut checked = 0usize;
            let mut worst = 0.0f64;
            for &x in samples {
                let Some(target) = eval_real(integrand, var, x, &ready) else {
                    continue;
                };
                if !target.is_finite() {
                    continue;
                }
                let h = 1e-5 * x.abs().max(1.0);
                let values: Option<Vec<f64>> = [-2.0f64, -1.0, 1.0, 2.0]
                    .iter()
                    .map(|offset| eval_real(antiderivative, var, x + offset * h, &ready))
                    .collect();
                let Some(values) = values else {
                    continue;
                };
                if !values.iter().all(|value| value.is_finite()) {
                    continue;
                }
                let derivative =
                    (-values[3] + 8.0 * values[2] - 8.0 * values[1] + values[0]) / (12.0 * h);
                checked += 1;
                let relative = (derivative - target).abs() / target.abs().max(1.0);
                if relative > worst {
                    worst = relative;
                }
            }
            if checked >= 2 && worst < 1e-6 {
                return;
            }
            if best.is_none() || checked > best.map(|b| b.0).unwrap_or(0) {
                best = Some((checked, worst));
            }
        }
        let (checked, worst) = best.unwrap_or((0, 0.0));
        panic!("no regime verified: best {checked} usable sample(s), worst {worst:e}");
    }

    #[test]
    fn b_class_core_identity() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "log(1 - x)/x");
        let answer =
            integrate_log_rational(&ctx, integrand, var).expect("∫ log(1−x)/x must be solved");
        assert!(!contains_integral(answer), "answer carries a residue");
        assert!(
            answer.to_string().contains("polylog(2"),
            "expected a dilogarithm, got {answer}"
        );
        check_numeric(
            &ctx,
            integrand,
            answer,
            var,
            &[],
            &[-0.9, -0.5, -0.2, 0.3, 0.6, 0.85],
        );
    }

    #[test]
    fn b_class_pole_scaling() {
        // `∫ (a + b·log(−2 + e·x))/x dx`: the log argument is affine with a
        // non-unit slope and a non-zero pole (Rubi's `rubi-00746` shape).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "(a + b*log(-2 + e*x))/x");
        let answer = integrate_log_rational(&ctx, integrand, var).expect("pole scaling");
        assert!(!contains_integral(answer));
        assert!(answer.to_string().contains("polylog(2"));
        let symbols = [Symbol::new("a"), Symbol::new("b"), Symbol::new("e")];
        check_numeric(
            &ctx,
            integrand,
            answer,
            var,
            &symbols,
            &[0.3, 0.6, 0.85, 1.4, -0.4],
        );
    }

    #[test]
    fn b_class_mobius_quotient() {
        // `log(e·(a+bx)/(c+dx))/(f+gx)`: `log` of a quotient splits, and the
        // resulting affine-log-over-linear terms are the pole entry
        // (Rubi's `rubi-01004` shape).
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "(A + B*log(e*(a + b*x)/(c + d*x)))/(f + g*x)");
        let answer = integrate_log_rational(&ctx, integrand, var).expect("Möbius quotient");
        assert!(!contains_integral(answer));
        assert!(answer.to_string().contains("polylog(2"));
        let symbols = [
            Symbol::new("A"),
            Symbol::new("B"),
            Symbol::new("a"),
            Symbol::new("b"),
            Symbol::new("c"),
            Symbol::new("d"),
            Symbol::new("e"),
            Symbol::new("f"),
            Symbol::new("g"),
        ];
        check_numeric(&ctx, integrand, answer, var, &symbols, &GATE_SAMPLES);
    }

    #[test]
    fn b_class_polynomial_weight() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "x*log(1 - x)");
        let answer = integrate_log_rational(&ctx, integrand, var).expect("polynomial weight");
        check_numeric(&ctx, integrand, answer, var, &[], &[-0.6, -0.2, 0.3, 0.7]);
    }

    #[test]
    fn a_class_basic_forms() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let c = Symbol::new("c");
        // Bare dilogarithm: `∫ Li₂(c·x) dx`.
        let integrand = parse(&ctx, "polylog(2, c*x)");
        let answer = integrate_polylog_head(&ctx, integrand, var).expect("∫ Li₂(cx)");
        check_numeric(&ctx, integrand, answer, var, &[c], &[0.2, 0.4, -0.3, -0.6]);
        // Monomial weight: `∫ x·Li₂(c·x) dx`.
        let integrand = parse(&ctx, "x*polylog(2, c*x)");
        let answer = integrate_polylog_head(&ctx, integrand, var).expect("∫ x·Li₂(cx)");
        check_numeric(&ctx, integrand, answer, var, &[c], &[0.2, 0.4, -0.3, -0.6]);
        // Inverse-square weight (Rubi's `rubi-00050` shape):
        // `∫ Li₂(c·(a+b·x))/x² dx`.
        let integrand = parse(&ctx, "polylog(2, c*(a + b*x))/x^2");
        let answer = integrate_polylog_head(&ctx, integrand, var).expect("rubi-00050 shape");
        let symbols = [c, Symbol::new("a"), Symbol::new("b")];
        // The answer carries `Li₂(u)` and `Li₂(1−u)`, so the usable window is
        // `u ∈ [0, 1]`: with the first regime (`c=2, a=−0.75, b=3`) that is
        // `x ∈ (0.25, 0.4167)`. Sample inside it; the adaptive gate does the
        // same automatically.
        check_numeric(
            &ctx,
            integrand,
            answer,
            var,
            &symbols,
            &[0.28, 0.33, 0.38, 0.41],
        );
    }

    #[test]
    fn refusals_are_honest() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        // Non-affine logarithm argument.
        assert!(
            integrate_log_rational(&ctx, parse(&ctx, "log(d + e*x^2)/(1 - x^2)"), var).is_none()
        );
        // Product of logarithms.
        assert!(integrate_log_rational(&ctx, parse(&ctx, "log(x)*log(1 - x)"), var).is_none());
        // Second-order pole needs a by-parts recursion outside the fragment.
        assert!(integrate_log_rational(&ctx, parse(&ctx, "log(1 - x)/x^2"), var).is_none());
        // Radical weight.
        assert!(integrate_log_rational(&ctx, parse(&ctx, "sqrt(x)*log(1 - x)"), var).is_none());
        // Other function heads belong to other stages.
        assert!(integrate_log_rational(&ctx, parse(&ctx, "exp(x)*log(1 - x)"), var).is_none());
        // No logarithm at all.
        assert!(integrate_log_rational(&ctx, parse(&ctx, "1/(1 + x^2)"), var).is_none());
        // Nested / higher-order polylogarithm heads.
        assert!(integrate_polylog_head(&ctx, parse(&ctx, "polylog(3, x)"), var).is_none());
        assert!(integrate_polylog_head(&ctx, parse(&ctx, "polylog(2, x^2)"), var).is_none());
        assert!(
            integrate_polylog_head(&ctx, parse(&ctx, "polylog(2, x)*polylog(2, x)"), var).is_none()
        );
    }

    #[test]
    fn declines_inside_parts_scope() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "log(1 - x)/x");
        let outside = integrate_log_rational(&ctx, integrand, var);
        assert!(outside.is_some());
        let inside = chain::with_parts_scope(|| integrate_log_rational(&ctx, integrand, var));
        assert!(inside.is_none(), "must decline inside a parts scope");
    }

    #[test]
    fn stage_output_does_not_re_enter() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "log(1 - x)/x");
        let answer = integrate_log_rational(&ctx, integrand, var).expect("core identity");
        // The emitted answer carries the dilogarithm head, so the A/B matchers
        // decline it; re-entering the stage cannot ping-pong.
        assert!(integrate_log_rational(&ctx, answer, var).is_none());
        assert!(integrate_polylog_head(&ctx, answer, var).is_none());
    }

    /// The emission gate must *reject* wrong candidates, not just accept the
    /// ones the stage happens to produce. This is the negative half of the
    /// "never emit an unchecked answer" contract.
    #[test]
    fn gate_rejects_wrong_candidates() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let integrand = parse(&ctx, "log(1 - x)/x");
        // A global sign error (the W6a failure mode, in the dilogarithm form).
        let wrong_sign = parse(&ctx, "-1*log(1 - x)/x");
        assert!(
            !gate_holds(&ctx, integrand, wrong_sign, var),
            "the gate accepted a candidate with the wrong global sign"
        );
        // A constant-scale error.
        let wrong_scale = parse(&ctx, "2*log(1 - x)/x");
        assert!(
            !gate_holds(&ctx, integrand, wrong_scale, var),
            "the gate accepted a candidate scaled by 2"
        );
        // A candidate that is not even an antiderivative.
        assert!(
            !gate_holds(&ctx, integrand, ctx.num(0), var),
            "the gate accepted the zero candidate"
        );
        // Sanity: the real answer passes (the rejections above are not
        // universal).
        let good = parse(&ctx, "log(x)*log(1 - x) + polylog(2, 1 - x)");
        assert!(
            gate_holds(&ctx, integrand, good, var),
            "the gate rejected the correct antiderivative"
        );
    }

    /// Every emitted B-class answer must also hold at abscissae the emission
    /// gate does **not** sample: the gate is a filter, not the acceptance
    /// criterion, and the independent oracle runs on its own grid.
    #[test]
    fn emitted_answers_hold_off_gate_samples() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        // None of these is in `GATE_SAMPLES` (nor an adaptive `u = ±½` point).
        const EXTRA: [f64; 6] = [0.05, 0.12, 0.22, 0.47, 0.62, 0.95];
        for source in ["log(1 - x)/x", "x*log(1 - x)"] {
            let integrand = parse(&ctx, source);
            let answer = integrate_log_rational(&ctx, integrand, var)
                .unwrap_or_else(|| panic!("{source} must be solved"));
            check_numeric(&ctx, integrand, answer, var, &[], &EXTRA);
        }
        // A parameterized pole form: extra abscissae, all three regimes tried
        // by `check_numeric` (the first usable regime verifies).
        let integrand = parse(&ctx, "(a + b*log(-2 + e*x))/x");
        let answer = integrate_log_rational(&ctx, integrand, var).expect("pole form");
        check_numeric(&ctx, integrand, answer, var, &[], &EXTRA);
    }

    /// The gate samples a symbol named `e` as a free parameter (the Rubi
    /// corpus uses `e` as a free symbol) while `pi` stays the mathematical
    /// constant. Both conventions are pinned here, because getting either
    /// wrong either fixes the corpus `e` to Euler's number or wrongly scales a
    /// correct answer carrying π.
    #[test]
    fn gate_treats_e_as_a_parameter_and_pi_as_a_constant() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let mut parameters: Vec<Symbol> = Vec::new();
        let integrand = parse(&ctx, "(a + b*log(-2 + e*x))/x");
        collect_parameters(integrand, var, &mut parameters);
        assert!(
            parameters.iter().any(|symbol| symbol.as_str() == "e"),
            "the gate must sample the corpus parameter `e`; got {parameters:?}"
        );
        let mut parameters: Vec<Symbol> = Vec::new();
        let pi_answer = parse(&ctx, "pi^2*6^-1 + log(x)");
        collect_parameters(pi_answer, var, &mut parameters);
        assert!(
            !parameters.iter().any(|symbol| symbol.as_str() == "pi"),
            "`pi` must stay a constant in the gate; got {parameters:?}"
        );
    }

    /// Diagnostic (not a gate): time the B-class stage over the 1892
    /// fallbacks, so a wall-clock movement between two corpus runs can be
    /// attributed to this stage instead of to the machine.
    ///
    /// ```text
    /// cargo test --release -p ocas-calc --lib corpus_log_stage_timing -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "corpus timing diagnostic; needs ocas-tests/data"]
    fn corpus_log_stage_timing() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../ocas-tests/data/integrate_1892_failures.jsonl");
        let Ok(text) = std::fs::read_to_string(&path) else {
            eprintln!("corpus failures file missing: {}", path.display());
            return;
        };
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let var = Symbol::new("x");
        let mut total = std::time::Duration::ZERO;
        let mut count = 0usize;
        let mut slow: Vec<(u128, String)> = Vec::new();
        for line in text.lines() {
            let Some(start) = line.find("\"integrand\": \"") else {
                continue;
            };
            let rest = &line[start + 14..];
            let Some(end) = rest.find("\"}") else {
                continue;
            };
            let source = &rest[..end];
            let Ok(expression) = ocas_parse::parse(&ctx, source) else {
                continue;
            };
            let before = std::time::Instant::now();
            let _ = integrate_log_rational(&ctx, expression, var);
            let elapsed = before.elapsed();
            total += elapsed;
            count += 1;
            if elapsed.as_millis() > 20 {
                slow.push((elapsed.as_millis(), source.to_string()));
            }
        }
        slow.sort_by_key(|entry| std::cmp::Reverse(entry.0));
        eprintln!("log_rational stage: {count} fallbacks, total {total:?}");
        for (ms, source) in slow.iter().take(10) {
            eprintln!("  {ms}ms {source}");
        }
    }
}
