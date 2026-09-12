//! Symbolic integration for oCAS.
//!
//! This module provides [`integrate`], a heuristic integrator for expressions
//! involving polynomials and elementary functions. It uses a lookup table for
//! common antiderivatives and supports simple linear substitutions.
//!
//! Integrals that cannot be expressed with the built-in table are returned as
//! the unevaluated form `Integral(expr, var)`.

#![allow(clippy::collapsible_if)]
// The chain-budget thread_local initializer is already const; on targets
// whose std thread_local implementation does not const-mark the init fn,
// the lint fires anyway (clippy#12276 acknowledges the backend dependence).
#![allow(clippy::missing_const_for_thread_local)]

pub(crate) mod binomial;
pub(crate) mod elliptic;
pub(crate) mod exp_log;
pub(crate) mod halfpower;
pub(crate) mod heuristic;
pub(crate) mod hyperbolic_reduction;
pub(crate) mod inverse_trig;
pub(crate) mod kernel_subst;
pub(crate) mod quad_power;
pub mod rational;
pub(crate) mod rde;
pub(crate) mod risch;
pub(crate) mod rules;
pub(crate) mod rules_ext;
pub(crate) mod special;
pub(crate) mod sqrt_quadratic;
pub(crate) mod symbolic_rational;
pub(crate) mod trig;
pub(crate) mod trig_kernel;
pub(crate) mod trig_reduce;
pub(crate) mod trig_reduction;

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_core::error::Result;
use ocas_core::fuel::Fuel;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::{simplify, simplify_with_fuel};

use crate::rules::calculus_rules;

/// Maximum recursion depth for `integrate_raw`, preventing infinite loops
/// on patterns such as nested linear substitutions if the table is misapplied.
const MAX_DEPTH: usize = 8;

/// Maximum number of `try_risch_or_fallback` chain entries per top-level
/// `integrate` call. The per-stage budgets (structural depth, rule depth,
/// parts depth) reset at substitution boundaries (Weierstrass, rule
/// residuals, expansion retries), so a cyclic interaction between stages
/// — observed in the wild: parts ↔ Weierstrass ping-pong on t-forms
/// carrying `atan(_t)` factors — can otherwise loop until the stack
/// overflows. Legitimate integrations use far fewer entries (typically
/// < 50), so tripping the budget degrades to the unevaluated form.
const MAX_CHAIN_ENTRIES: u32 = 256;

thread_local! {
    static CHAIN_ENTRIES: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

/// Reset the chain-entry budget; called by every public entry point.
fn reset_chain_budget() {
    CHAIN_ENTRIES.with(|c| c.set(0));
}

/// Consume one chain entry; returns true when the budget is exhausted.
fn chain_budget_exhausted() -> bool {
    CHAIN_ENTRIES.with(|c| {
        let v = c.get().saturating_add(1);
        c.set(v);
        v > MAX_CHAIN_ENTRIES
    })
}

/// Whether stage tracing is enabled (`OCAS_INTEGRATE_TRACE=1`).
///
/// Diagnostic only: a hang leaves the last `enter <stage>` line on stderr,
/// which attributes the case to a pipeline stage without guesswork.
fn trace_enabled() -> bool {
    thread_local! {
        static TRACE: bool = std::env::var_os("OCAS_INTEGRATE_TRACE")
            .is_some_and(|v| v != "0");
    }
    TRACE.with(|t| *t)
}

/// Run an `Option`-returning pipeline stage under the trace switch.
fn traced_stage<'a, T>(stage: &str, expr: Atom<'a>, f: impl FnOnce() -> Option<T>) -> Option<T> {
    if !trace_enabled() {
        return f();
    }
    eprintln!("[trace] enter  {stage} :: {expr}");
    let out = f();
    if out.is_none() {
        eprintln!("[trace] decline {stage}");
    }
    out
}

/// Trace a stage whose control flow cannot be expressed as `Option` (blocks
/// that re-enter the chain).
fn trace_enter(stage: &str, expr: Atom<'_>) {
    if trace_enabled() {
        eprintln!("[trace] enter  {stage} :: {expr}");
    }
}

thread_local! {
    /// Re-entry guard for the bounded-expansion pre-pass.
    static EXPAND_PREPASS_ACTIVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether `expr` contains any function head (transcendental node).
fn has_function_head(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Fun(_, _) => true,
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| has_function_head(*a)),
        AtomNode::Pow(b, e) => has_function_head(*b) || has_function_head(*e),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Bounded-distribution pre-pass for unexpanded products.
///
/// An unexpanded product of sums hangs the symbolic-rational backend on
/// shapes whose expansion is a small Laurent polynomial (verified:
/// `(A + B*x^2)*(b*x^2 + c*x^4)/x^6` never returns, while its expanded form
/// solves in milliseconds). Distributing first — under the same term budget
/// the late retry uses — turns those into ordinary termwise integrals.
///
/// The pass is deliberately restricted to **purely algebraic** expansions
/// (no `Fun` head anywhere): transcendental products are owned by the later
/// rule-table / trig / hyperbolic mechanisms, and expanding them here would
/// steal those cases from the stages built for them (measured: termwise
/// integration of `sec*cot^3`-style products routed work into paths that
/// returned wrong answers, while the unexpanded shape had a correct owner
/// downstream).
///
/// Returns `None` when the input is not a product of two or more
/// non-constant factors, when the expansion contains a function head, when
/// the expansion is already done, or when the expanded form does not
/// integrate. The re-entry guard makes the pass idempotent: a nested call on
/// the same shape declines immediately.
fn expand_prepass<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
) -> Option<Atom<'a>> {
    let AtomNode::Mul(args) = expr.node() else {
        return None;
    };
    if args.iter().filter(|a| !is_constant(**a, var)).count() < 2 {
        return None;
    }
    // Cheap pre-filter: a product with no sum factor cannot distribute.
    if !args.iter().any(|a| matches!(a.node(), AtomNode::Add(_))) {
        return None;
    }
    if has_function_head(expr) || node_count(expr) > 64 {
        return None;
    }
    if EXPAND_PREPASS_ACTIVE.with(|c| c.replace(true)) {
        return None;
    }
    let out = (|| {
        let expanded = crate::expand::expand_bounded(ctx, expr)?;
        if has_function_head(expanded) {
            return None;
        }
        let folded = crate::ode::util::collect_terms(ctx, expanded);
        if node_count(folded) > 256 {
            return None;
        }
        let r = integrate_raw(ctx, folded, var, 0, rules_enabled, rule_depth, parts_depth);
        if contains_integral(r) { None } else { Some(r) }
    })();
    EXPAND_PREPASS_ACTIVE.with(|c| c.set(false));
    out
}

/// Read-only, allocation-free pre-check for [`fold_linear_squares`]: does
/// `expr` contain a three-term sum at all? Only such a sum can be a squared
/// linear form, so this lets the common case skip the rewriting walk entirely.
/// Short-circuits on the first match.
fn contains_three_term_sum(expr: Atom<'_>) -> bool {
    match expr.node() {
        AtomNode::Add(args) => args.len() == 3 || args.iter().any(|a| contains_three_term_sum(*a)),
        AtomNode::Mul(args) => args.iter().any(|a| contains_three_term_sum(*a)),
        AtomNode::Pow(b, e) => contains_three_term_sum(*b) || contains_three_term_sum(*e),
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_three_term_sum(*a)),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Exact structural square root of `term`, when it is *syntactically* a square:
/// `u²`, a product of such squares, or a non-negative perfect-square integer.
///
/// This is only a candidate generator for [`try_fold_linear_square`]; the
/// rewrite is accepted by exact re-expansion, never by this test, so an
/// imperfect recogniser can only miss a fold, not produce a wrong one.
fn structural_sqrt<'a>(ctx: &'a AtomArena<'a>, term: Atom<'a>) -> Option<Atom<'a>> {
    match term.node() {
        AtomNode::Pow(b, e) if matches!(e.node(), AtomNode::Num(2)) => Some(*b),
        AtomNode::Num(n) if *n >= 0 => {
            let r = (*n as f64).sqrt() as i64;
            if r * r == *n { Some(ctx.num(r)) } else { None }
        }
        AtomNode::Mul(factors) => {
            let mut roots = Vec::with_capacity(factors.len());
            for f in factors.iter() {
                roots.push(structural_sqrt(ctx, *f)?);
            }
            Some(ctx.mul(&roots))
        }
        _ => None,
    }
}

/// Fold `p² + 2·p·q + q²` into `(p + q)²` throughout `expr`.
///
/// This is an exact algebraic identity, and the fold is only accepted when
/// expanding the candidate square reproduces the original sum term for term
/// (see [`try_fold_linear_square`]), so a false positive is impossible.
///
/// It exists because the corpus writes a *linear form squared* as an expanded
/// trinomial — `a² + 2·a·b·x + b²·x²` — which the symbolic-rational backend and
/// the radical engines grind on (measured: `rubi-00854` never returns) while
/// the folded form is an ordinary partial-fraction or polynomial problem.
///
/// `fold_here` is `false` for the base of a non-integer power (including the
/// `sqrt` head, which is a `Fun` rather than a `Pow`): `(T)^{k/2}` with
/// `T = (p+q)²` is the same function as the original radicand, but the radical
/// and elliptic engines match on the *expanded* quadratic, and
/// `((p+q)²)^{1/2} = |p+q| ≠ p+q` — folding there would either lose those
/// engines or introduce a branch error. Integer powers are folded, where
/// `((p+q)²)^n = (p+q)^{2n}` is exact.
///
/// The fold also requires `p + q` to be **affine in the integration variable**
/// (see [`try_fold_linear_square`]); that is what keeps it from re-forming the
/// trigonometric squares the product-to-sum stage needs expanded.
fn fold_linear_squares<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fold_here: bool,
) -> Atom<'a> {
    let rebuilt = match expr.node() {
        AtomNode::Add(args) => {
            let args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| fold_linear_squares(ctx, *a, var, true))
                .collect();
            ctx.add(&args)
        }
        AtomNode::Mul(args) => {
            let args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| fold_linear_squares(ctx, *a, var, true))
                .collect();
            ctx.mul(&args)
        }
        AtomNode::Pow(b, e) => {
            let integer_power = matches!(e.node(), AtomNode::Num(_));
            let b = fold_linear_squares(ctx, *b, var, integer_power);
            let e = fold_linear_squares(ctx, *e, var, true);
            ctx.pow(b, e)
        }
        AtomNode::Fun(name, args) => {
            // `sqrt(u)` is a `Fun` head, not a `Pow`, but it is still a
            // non-integer power: an argument that is a perfect-square
            // trinomial must not be folded under it (`|p+q| ≠ p+q`).
            let fold_args = name.as_str() != "sqrt";
            let args: Vec<Atom<'a>> = args
                .iter()
                .map(|a| fold_linear_squares(ctx, *a, var, fold_args))
                .collect();
            ctx.fun(name.as_str(), &args)
        }
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
    };
    if fold_here && let AtomNode::Add(args) = rebuilt.node() {
        return try_fold_linear_square(ctx, args, var).unwrap_or(rebuilt);
    }
    rebuilt
}

/// `p² + 2·p·q + q² → (p+q)²`, confirmed by exact re-expansion.
///
/// Each pair of terms with a structural square root is tried; the candidate is
/// accepted only when `expand((p+q)²)` reproduces the input sum exactly, so the
/// rewrite can never change the value of the expression. Both sides go through
/// [`crate::ode::util::collect_terms`] first: the corpus writes the cross term
/// as `((2·a)·b)·x`, which is the same monomial as `2·a·b·x` but not the same
/// atom, and the comparison must be about the mathematics, not the spelling.
///
/// `p + q` must additionally be **affine in `var`**. Without that restriction
/// the fold also re-forms trigonometric squares — `a²cos²u + 2ab·cos u·sin u +
/// b²sin²u → (a·cos u + b·sin u)²` — undoing the expansion the product-to-sum
/// stage (`trig_reduce`) needs to see, which turned `∫(a·cos u + b·sin u)² dx`
/// (corpus `rubi-00334`) from a solved case into a per-case timeout. The
/// module's purpose is the expanded *linear* square the corpus writes, so the
/// affine test both fixes that regression and states the intent.
fn try_fold_linear_square<'a>(
    ctx: &'a AtomArena<'a>,
    terms: &[Atom<'a>],
    var: Symbol,
) -> Option<Atom<'a>> {
    // The parser builds left-nested sums (`(a² + 2abx) + b²x²`) and only
    // `normalize` flattens them; flatten here so the fold also works on
    // expressions the pipeline assembled itself.
    let mut flat: Vec<Atom<'a>> = Vec::with_capacity(terms.len());
    flatten_add(terms, &mut flat);
    let terms: &[Atom<'a>] = &flat;
    if terms.len() != 3 {
        return None;
    }
    let original = ctx.add(terms);
    if node_count(original) > 64 {
        return None;
    }
    let roots: Vec<Option<Atom<'a>>> = terms.iter().map(|t| structural_sqrt(ctx, *t)).collect();
    if roots.iter().filter(|r| r.is_some()).count() < 2 {
        return None;
    }
    // `rhs` (and the expand-and-collect fallback below) are only needed when
    // the cheap normalized comparison fails; computing `rhs` lazily keeps the
    // common "three-term sum that is not a square" case cheap, which matters
    // because this runs on every pipeline entry.
    let mut rhs: Option<Atom<'a>> = None;
    for (i, pi) in roots.iter().enumerate() {
        for (j, qj) in roots.iter().enumerate() {
            if i == j {
                continue;
            }
            let (Some(p), Some(q)) = (pi, qj) else {
                continue;
            };
            let base = ctx.add(&[*p, *q]);
            // Only an affine base: see the doc comment.
            let (slope, _intercept) = linear_form(ctx, base, var)?;
            if matches!(slope.node(), AtomNode::Num(0)) {
                continue;
            }
            let candidate = ctx.pow(base, ctx.num(2));
            // Fast, sound acceptance: the two square terms plus a third term
            // that is already `2·p·q` after `normalize` is a complete proof of
            // the identity. `normalize` flattens products, so the corpus's
            // `((2·a)·b)·x` matches `2·a·(b·x)` here.
            let k = 3 - i - j;
            let cross = normalize(ctx, ctx.mul(&[ctx.num(2), *p, *q]));
            if normalize(ctx, terms[k]) == cross {
                return Some(candidate);
            }
            // Fallback for spellings `normalize` does not canonicalise: prove
            // the identity by exact re-expansion.
            let rhs = match rhs {
                Some(r) => r,
                None => {
                    let r = crate::ode::util::collect_terms(ctx, original);
                    rhs = Some(r);
                    r
                }
            };
            let Some(expanded) = crate::expand::expand_bounded(ctx, candidate) else {
                continue;
            };
            let lhs = crate::ode::util::collect_terms(ctx, expanded);
            if lhs == rhs {
                return Some(candidate);
            }
        }
    }
    None
}

/// Flatten a possibly nested sum into its top-level terms.
fn flatten_add<'a>(terms: &[Atom<'a>], out: &mut Vec<Atom<'a>>) {
    for t in terms {
        match t.node() {
            AtomNode::Add(inner) => flatten_add(inner, out),
            _ => out.push(*t),
        }
    }
}

/// Options controlling the integration pipeline.
///
/// `rules` enables the rule-table engine (default `true`); set it to `false`
/// to restore the pre-0.27 behaviour for comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntegrateOptions {
    /// Enable the rule-table integration engine (default: enabled).
    pub rules: bool,
}
impl Default for IntegrateOptions {
    fn default() -> Self {
        Self { rules: true }
    }
}

/// Integrate `expr` with respect to `var` using the given options.
pub fn integrate_with_options<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    options: IntegrateOptions,
) -> Atom<'a> {
    let normalized = normalize(ctx, expr);
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    reset_chain_budget();
    let raw = integrate_raw(ctx, normalized, var, 0, options.rules, 0, 0);
    // Combine default algebraic simplification with calculus-specific rules,
    // then normalize to a canonical form (removing *1, +0, sorting, etc.).
    let after_default = simplify(ctx, raw, &default_rules, 20);
    let after_calc = simplify(ctx, after_default, &calc_rules, 10);
    normalize(ctx, after_calc)
}

/// Integrate `expr` with respect to `var`.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::integrate;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.pow(x, ctx.num(2));
/// let result = integrate(&ctx, expr, Symbol::new("x"));
/// assert_eq!(result.to_string(), "(3^-1)*(x^3)");
/// ```
///
/// For integrals not covered by the heuristic table, the result is returned
/// as `Integral(expr, var)`.
pub fn integrate<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    integrate_with_options(ctx, expr, var, IntegrateOptions::default())
}

/// Integrate with a [`Fuel`] budget bounding the post-integration
/// simplification passes.
///
/// The integration traversal itself uses the internal depth limit; this entry
/// point threads `fuel` through the two simplification stages so a pathological
/// result that would otherwise spin the rewriter can be cut off determin-
/// istically. Returns `Err` only when fuel was exhausted mid-simplification.
///
/// # Example
///
/// ```no_run
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_core::arena::Arena;
/// use ocas_core::fuel::Fuel;
/// use ocas_calc::integral::integrate_with_fuel;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let expr = ctx.var("x");
/// let fuel = Fuel::new(500);
/// let result = integrate_with_fuel(&ctx, expr, Symbol::new("x"), &fuel);
/// match result {
///     Ok(r) => println!("{}", r),
///     Err(_) => println!("fuel exhausted during simplification"),
/// }
/// ```
pub fn integrate_with_fuel<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    fuel: &Fuel,
) -> Result<Atom<'a>> {
    let normalized = normalize(ctx, expr);
    let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
    let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
    reset_chain_budget();
    let raw = integrate_raw(ctx, normalized, var, 0, true, 0, 0);
    let after_default = simplify_with_fuel(ctx, raw, &default_rules, 20, fuel)?;
    let after_calc = simplify_with_fuel(ctx, after_default, &calc_rules, 10, fuel)?;
    Ok(normalize(ctx, after_calc))
}

/// Try heuristic integration techniques (parts, trig substitution,
/// Weierstrass, Euler) on `expr`.
///
/// Returns the integrated expression if a heuristic succeeds, or the
/// unevaluated `Integral(expr, var)` form if none do.
pub fn integrate_heuristic<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    let normalized = normalize(ctx, expr);
    reset_chain_budget();
    if let Some(r) = heuristic::heuristic_integrate(ctx, normalized, var, 0) {
        let calc_rules = calculus_rules(ctx, &crate::pattern_alloc::VecAlloc);
        let default_rules = default_rules(ctx, &crate::pattern_alloc::VecAlloc);
        let after_default = simplify(ctx, r, &default_rules, 20);
        let after_calc = simplify(ctx, after_default, &calc_rules, 10);
        return normalize(ctx, after_calc);
    }
    fallback(ctx, expr, var)
}

pub(crate) fn integrate_raw<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    depth: usize,
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
) -> Atom<'a> {
    if depth > MAX_DEPTH {
        return fallback(ctx, expr, var);
    }

    // Exact algebraic fold: `p² + 2·p·q + q² → (p+q)²` (0.27.3), for a base
    // affine in `var`. Applied at the pipeline entry so every downstream stage
    // — and every re-entered residual — sees the factored shape. Idempotent:
    // the folded form has no matching three-term sum, so a nested call is a
    // no-op.
    let expr = if contains_three_term_sum(expr) {
        fold_linear_squares(ctx, expr, var, true)
    } else {
        expr
    };

    match expr.node() {
        AtomNode::Num(_) => {
            // ∫ c dx = c * x
            let x = ctx.var(var.as_str());
            ctx.mul(&[expr, x])
        }
        AtomNode::Var(v) => {
            if *v == var {
                ctx.mul(&[
                    ctx.pow(ctx.var(var.as_str()), ctx.num(2)),
                    ctx.pow(ctx.num(2), ctx.num(-1)),
                ])
            } else {
                ctx.mul(&[expr, ctx.var(var.as_str())])
            }
        }
        AtomNode::Add(args) => {
            let mut terms = Vec::with_capacity(args.len());
            for a in args.iter() {
                terms.push(integrate_raw(
                    ctx,
                    *a,
                    var,
                    depth,
                    rules_enabled,
                    rule_depth,
                    parts_depth,
                ));
            }
            ctx.add(&terms)
        }
        AtomNode::Mul(args) => {
            let r = integrate_product(
                ctx,
                args,
                var,
                depth,
                rules_enabled,
                rule_depth,
                parts_depth,
            );
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var, rules_enabled, rule_depth, parts_depth)
            } else {
                r
            }
        }
        AtomNode::Pow(base, exp) => {
            let r = integrate_power(ctx, *base, *exp, var, depth);
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var, rules_enabled, rule_depth, parts_depth)
            } else {
                r
            }
        }
        AtomNode::Fun(name, args) => {
            let r = integrate_function(ctx, *name, args, var, depth);
            if is_fallback(&r) {
                try_risch_or_fallback(ctx, expr, var, rules_enabled, rule_depth, parts_depth)
            } else {
                r
            }
        }
    }
}

/// Try the rational-function integrator, then the Risch algorithm, then the
/// rule table, before giving up with the unevaluated `Integral` form.
fn try_risch_or_fallback<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
) -> Atom<'a> {
    // Backstop against cyclic stage interactions (see MAX_CHAIN_ENTRIES).
    if chain_budget_exhausted() {
        return fallback(ctx, expr, var);
    }
    if let Some(r) = traced_stage("rational", expr, || {
        rational::integrate_rational(ctx, expr, var)
    }) {
        return r;
    }
    // Symbolic-constant rationals (coefficients in ℚ(symbols)): the ℚ
    // backend declines these; the symbolic backend also powers the
    // trig-rational class through Weierstrass t-rationals.
    // Quadratic/linear denominator power reductions (0.27.1 Phase 2):
    // closed-form recurrences for `P(x)/q^n` — these shapes stall the
    // symbolic rational backend's multivariate coefficient gcd, so the
    // recurrence must preempt it.
    if let Some(r) = traced_stage("quad_power", expr, || {
        quad_power::integrate_quad_power(ctx, expr, var)
    }) {
        return r;
    }
    // Closed-form kernel families (0.27.2): they must preempt
    // `symbolic_rational`/`risch`, whose field-Euclidean steps grind
    // unboundedly on these shapes with symbolic coefficients (verified:
    // `sinh(x)^3/(a+b*sinh(x))` and the Weierstrass t-forms of
    // `1/(a+b*cos(u)+c*sin(u))^2` never return).
    //
    // The bounded-expansion pre-pass runs first: unexpanded products whose
    // expansion is a small Laurent polynomial must not reach the symbolic
    // backend at all (0.27.2 A1).
    if let Some(r) = expand_prepass(ctx, expr, var, rules_enabled, rule_depth, parts_depth) {
        trace_enter("expand_prepass", expr);
        return r;
    }
    if let Some(r) = traced_stage("kernel_subst", expr, || {
        kernel_subst::integrate_kernel_subst(ctx, expr, var)
    }) {
        return r;
    }
    if let Some(r) = traced_stage("hyperbolic_reduction", expr, || {
        hyperbolic_reduction::integrate_hyperbolic_reduction(ctx, expr, var)
    }) {
        return r;
    }
    if let Some(r) = traced_stage("symbolic_rational", expr, || {
        symbolic_rational::integrate_rational_symbolic(ctx, expr, var)
    }) {
        return r;
    }
    if let Some(r) = traced_stage("risch", expr, || risch::risch_integrate(ctx, expr, var)) {
        return r;
    }
    // Trigonometric integrands: rewrite into complex exponentials, run
    // Risch, then try to bring the answer back to real form. The tower
    // grinds on symbolic linear arguments (`cos(c + d*x)`), so this stage
    // only runs when every sin/cos argument has numeric coefficients.
    if trig::trig_args_numeric(ctx, expr, var)
        && let Some(exp_form) = trig::trig_to_exp(ctx, expr)
        && let Some(complex_ans) = traced_stage("risch(trig-exp)", expr, || {
            risch::risch_integrate(ctx, exp_form, var)
        })
    {
        return trig::realify(ctx, complex_ans);
    }
    // Non-elementary integrals with special-function closed forms
    // (erf, Ei, Si, Ci, Fresnel, …).
    if let Some(r) = traced_stage("special", expr, || {
        special::special_integrate(ctx, expr, ctx.var(var.as_str()))
    }) {
        return r;
    }
    // Rule-table engine: standard-calculus breadth rules with residual
    // `Integral(g, x)` reduction formulas. The table is only built when the
    // earlier stages failed, so Risch-solvable problems pay no parse cost.
    //
    // NOTE: this runs BEFORE the heuristic stage (the 0.27 plan originally
    // placed it after). Running it first is required for reachability:
    // trig products like `sin(x)^2*cos(x)^5` otherwise fall into the
    // Weierstrass substitution, whose t-rational blows up the rational
    // backend and never returns, so the rules never fire on the exact
    // shapes the rule library is built for. Rules only fire where the
    // rational/Risch/special stages failed, and the heuristic still runs
    // for everything the rules do not cover.
    if rules_enabled
        && let Some(table) = rules::build_rule_table(ctx, var)
        && let Some(r) = traced_stage("rules", expr, || {
            rules::integrate_rules(ctx, &table, expr, var, rule_depth)
        })
    {
        return r;
    }
    // General quadratic-radical engine (0.27.1): direct forms for
    // √(a+b·x+c·x²) composites, reciprocal forms, Euler III. Runs BEFORE
    // binomial's Chebyshev cases: both accept `q^±1/2`-style radicands,
    // and this engine's asin/log direct forms are the canonical answers
    // (Chebyshev's t-form back-substitution produces uglier atan shapes).
    if let Some(r) = traced_stage("sqrt_quadratic", expr, || {
        sqrt_quadratic::integrate_sqrt_quadratic(ctx, expr, var)
    }) {
        return r;
    }
    // Chebyshev binomial differentials and fractional-power
    // rationalization: substitute to a rational t-form, reintegrate,
    // back-substitute. Declines (None) unless an exact integrability
    // condition holds.
    if let Some(r) = traced_stage("binomial", expr, || {
        binomial::integrate_binomial(ctx, expr, var)
    }) {
        return r;
    }
    // exp/log-kernel substitutions (0.27.1): rational-in-e^(ax) and
    // hyperbolic-rational forms, f(log x)/x, log-power gaps of rule B6.
    if let Some(r) = traced_stage("exp_log", expr, || {
        exp_log::integrate_exp_log(ctx, expr, var)
    }) {
        return r;
    }

    // Inverse-trig/hyperbolic mechanisms (0.27.1): kernel-derivative power
    // rule, inv-hyp substitution to hyperbolic t-forms, bare linear args.
    if let Some(r) = traced_stage("inverse_trig", expr, || {
        inverse_trig::integrate_inverse_trig(ctx, expr, var)
    }) {
        return r;
    }
    // Trig product-to-sum reduction: products of sin/cos at linear
    // arguments become a sum of single trig terms, then distribute and
    // integrate termwise. Runs before the heuristic stage: the reduction
    // yields clean multiple-angle forms where Weierstrass would return
    // tan(u/2) shapes (or grind on the t-rational).
    if let Some(reduced) = trig_reduce::trig_reduce_products(ctx, expr, var) {
        trace_enter("trig_reduce", reduced);
        let candidate = crate::expand::expand_bounded(ctx, reduced).unwrap_or(reduced);
        let folded = crate::ode::util::collect_terms(ctx, candidate);
        let r = integrate_raw(ctx, folded, var, 0, rules_enabled, rule_depth, parts_depth);
        if !is_fallback(&r) {
            return r;
        }
    }
    // Trig-denominator power reductions, linear-numerator decomposition and
    // polynomial×trig closed forms (0.27.1). Intercepts `1/(a+b·T(u))^n`
    // before Weierstrass blows the t-rational up, and `x^m·T(ax+b)` shapes
    // that parts cannot finish within budget.
    if let Some(r) = traced_stage("trig_reduction", expr, || {
        trig_reduction::integrate_trig_reduction(ctx, expr, var)
    }) {
        return r;
    }
    // Single-trig-kernel rational forms and tan/sec-family reductions
    // (0.27.1): after trig_reduction so plain `1/(a+b·T)^n` stays there.
    if let Some(r) = traced_stage("trig_kernel", expr, || {
        trig_kernel::integrate_trig_kernel(ctx, expr, var)
    }) {
        return r;
    }
    // Bounded distributive expansion: distribute products over sums and
    // integrate termwise. Runs BEFORE the heuristic stage: parts recursion
    // on multi-factor products can consume the whole chain budget, which
    // would starve this retry's per-term re-entries. The expansion is
    // budgeted and idempotent, so the re-entered chain cannot loop back
    // here on the same shape. Like terms are folded first so factors like
    // `x*x` reach the integrator as `x^2`.
    if let Some(expanded) = crate::expand::expand_bounded(ctx, expr) {
        trace_enter("expand_retry", expanded);
        let folded = crate::ode::util::collect_terms(ctx, expanded);
        let r = integrate_raw(ctx, folded, var, 0, rules_enabled, rule_depth, parts_depth);
        if !is_fallback(&r) {
            return r;
        }
    }
    // Half-power front-end then elliptic reduction (0.27.2 D): placed after
    // the bounded-expansion retry so the elementary engines and the expanded
    // single-term shapes keep first claim, and before the heuristic stage so
    // Weierstrass cannot route these radicals into the t-rational backend.
    if let Some(r) = traced_stage("halfpower", expr, || {
        halfpower::integrate_half_power(ctx, expr, var)
    }) {
        return r;
    }
    if let Some(r) = traced_stage("elliptic", expr, || {
        elliptic::integrate_elliptic(ctx, expr, var)
    }) {
        return r;
    }
    // Heuristic techniques: parts, trig sub, Weierstrass, Euler.
    if let Some(r) = traced_stage("heuristic", expr, || {
        heuristic::heuristic_integrate(ctx, expr, var, parts_depth)
    }) {
        return r;
    }
    fallback(ctx, expr, var)
}

pub(crate) fn fallback<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Atom<'a> {
    ctx.fun("Integral", &[expr, ctx.var(var.as_str())])
}

/// True if `expr` does not contain `var`.
pub(crate) fn is_constant<'a>(expr: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) => true,
        AtomNode::Var(v) => *v != var,
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().all(|a| is_constant(*a, var))
        }
        AtomNode::Pow(base, exp) => is_constant(*base, var) && is_constant(*exp, var),
    }
}

fn integrate_product<'a>(
    ctx: &'a AtomArena<'a>,
    args: &'a [Atom<'a>],
    var: Symbol,
    depth: usize,
    rules_enabled: bool,
    rule_depth: usize,
    parts_depth: usize,
) -> Atom<'a> {
    // Split into constant factors and the remaining factor.
    let mut constants: Vec<Atom<'a>> = Vec::new();
    let mut non_constant: Vec<Atom<'a>> = Vec::new();

    for a in args.iter() {
        if is_constant(*a, var) {
            constants.push(*a);
        } else {
            non_constant.push(*a);
        }
    }

    if non_constant.is_empty() {
        // All factors are constant: ∫ c dx = c * x
        return ctx.mul(&[ctx.mul(args), ctx.var(var.as_str())]);
    }

    let core = if non_constant.len() == 1 {
        non_constant[0]
    } else {
        ctx.mul(&non_constant)
    };

    let integrated_core = integrate_raw(
        ctx,
        core,
        var,
        depth + 1,
        rules_enabled,
        rule_depth,
        parts_depth,
    );

    // If integration failed, wrap the whole product.
    if is_fallback(&integrated_core) {
        return fallback(ctx, ctx.mul(args), var);
    }

    let mut result_factors = constants;
    result_factors.push(integrated_core);
    ctx.mul(&result_factors)
}

pub(crate) fn is_fallback<'a>(atom: &Atom<'a>) -> bool {
    matches!(atom.node(), AtomNode::Fun(name, _) if name.as_str() == "Integral")
}

/// Deep residue check: any `Integral(...)` node anywhere in the tree
/// (unlike `is_fallback`, which only inspects the head).
pub(crate) fn contains_integral<'a>(atom: Atom<'a>) -> bool {
    match atom.node() {
        AtomNode::Fun(name, args) => {
            name.as_str() == "Integral" || args.iter().any(|a| contains_integral(*a))
        }
        AtomNode::Add(args) | AtomNode::Mul(args) => args.iter().any(|a| contains_integral(*a)),
        AtomNode::Pow(b, e) => contains_integral(*b) || contains_integral(*e),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

// =========================================================================
// Shared substitution-mechanism plumbing (0.27.1: binomial / exp_log / …)
// =========================================================================

pub(crate) fn gcd_i64(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.unsigned_abs(), b.unsigned_abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    i64::try_from(a).unwrap_or(i64::MAX)
}

pub(crate) fn lcm_i64(a: i64, b: i64) -> Option<i64> {
    let g = gcd_i64(a, b);
    (a / g).checked_mul(b)
}

/// Total node count of the expression tree (saturating).
pub(crate) fn node_count(expr: Atom<'_>) -> usize {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => 1,
        AtomNode::Pow(b, e) => node_count(*b)
            .saturating_add(node_count(*e))
            .saturating_add(1),
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => args
            .iter()
            .fold(1usize, |acc, a| acc.saturating_add(node_count(*a))),
    }
}

/// True when `sym` occurs as a `Var` anywhere in `expr`.
pub(crate) fn contains_symbol(expr: Atom<'_>, sym: Symbol) -> bool {
    match expr.node() {
        AtomNode::Var(v) => *v == sym,
        AtomNode::Num(_) => false,
        AtomNode::Pow(b, e) => contains_symbol(*b, sym) || contains_symbol(*e, sym),
        AtomNode::Add(args) | AtomNode::Mul(args) | AtomNode::Fun(_, args) => {
            args.iter().any(|a| contains_symbol(*a, sym))
        }
    }
}

/// Pick a substitution variable that does not collide with `var` or with
/// any symbol already present in `expr`.
pub(crate) fn pick_subst_symbol(expr: Atom<'_>, var: Symbol) -> Option<Symbol> {
    for name in ["t", "u"] {
        let s = Symbol::new(name);
        if s != var && !contains_symbol(expr, s) {
            return Some(s);
        }
    }
    None
}

/// Substitute every `Var(sym)` occurrence in `expr` with `replacement`.
pub(crate) fn replace_symbol<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    sym: Symbol,
    replacement: Atom<'a>,
) -> Atom<'a> {
    match expr.node() {
        AtomNode::Var(v) => {
            if *v == sym {
                replacement
            } else {
                expr
            }
        }
        AtomNode::Num(_) => expr,
        AtomNode::Add(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_symbol(ctx, *a, sym, replacement))
                .collect();
            ctx.add(&rebuilt)
        }
        AtomNode::Mul(args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_symbol(ctx, *a, sym, replacement))
                .collect();
            ctx.mul(&rebuilt)
        }
        AtomNode::Pow(b, e) => {
            let nb = replace_symbol(ctx, *b, sym, replacement);
            let ne = replace_symbol(ctx, *e, sym, replacement);
            ctx.pow(nb, ne)
        }
        AtomNode::Fun(name, args) => {
            let rebuilt: Vec<Atom<'a>> = args
                .iter()
                .map(|a| replace_symbol(ctx, *a, sym, replacement))
                .collect();
            ctx.fun(name.as_str(), &rebuilt)
        }
    }
}

/// `base^e` for integer `e`, folding the degenerate exponents.
pub(crate) fn int_pow<'a>(ctx: &'a AtomArena<'a>, base: Atom<'a>, e: i64) -> Atom<'a> {
    match e {
        0 => ctx.num(1),
        1 => base,
        _ => ctx.pow(base, ctx.num(e)),
    }
}

/// `a^(−1)`, folded for numeric `a`.
pub(crate) fn inv<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>) -> Atom<'a> {
    match a.node() {
        AtomNode::Num(1) => ctx.num(1),
        AtomNode::Num(n) => rat_atom(ctx, 1, *n),
        _ => ctx.pow(a, ctx.num(-1)),
    }
}

/// Build the atom `p/q` (reduced, `q > 0`).
pub(crate) fn rat_atom<'a>(ctx: &'a AtomArena<'a>, p: i64, q: i64) -> Atom<'a> {
    if p == 0 {
        return ctx.num(0);
    }
    let (p, q) = if q < 0 { (-p, -q) } else { (p, q) };
    let g = gcd_i64(p, q);
    let (p, q) = (p / g, q / g);
    if q == 1 {
        return ctx.num(p);
    }
    ctx.mul(&[ctx.num(p), ctx.pow(ctx.num(q), ctx.num(-1))])
}

fn integrate_power<'a>(
    ctx: &'a AtomArena<'a>,
    base: Atom<'a>,
    exp: Atom<'a>,
    var: Symbol,
    _depth: usize,
) -> Atom<'a> {
    // ∫ f(x)^0 dx = ∫ 1 dx = x (any base; reaches this table through the
    // rules engine's reduction formulas, e.g. tan(x)^0).
    if matches!(exp.node(), AtomNode::Num(0)) {
        return ctx.var(var.as_str());
    }

    // Detect x^n where n is a constant integer.
    if let AtomNode::Var(v) = base.node()
        && *v == var
    {
        if let AtomNode::Num(n) = exp.node() {
            if *n == -1 {
                // ∫ x^(-1) dx = log(x)
                return ctx.fun("log", &[base]);
            }
            // ∫ x^n dx = x^(n+1) / (n+1)
            let new_exp = ctx.num(n + 1);
            let denom = ctx.num(n + 1);
            return ctx.mul(&[ctx.pow(base, new_exp), ctx.pow(denom, ctx.num(-1))]);
        }
        // Fractional exponents p/q: ∫ x^(p/q) dx = x^(p/q+1) / (p/q + 1).
        if let Some((p, q)) = fraction_exponent(exp) {
            if p != -q {
                // new exponent = (p+q)/q; coefficient = q/(p+q)
                let new_exp = ctx.mul(&[ctx.num(p + q), ctx.pow(ctx.num(q), ctx.num(-1))]);
                let denom = ctx.mul(&[ctx.num(p + q), ctx.pow(ctx.num(q), ctx.num(-1))]);
                return ctx.mul(&[ctx.pow(base, new_exp), ctx.pow(denom, ctx.num(-1))]);
            }
        }
    }

    // Detect linear substitution: (a*x + b)^n where n is constant integer.
    if let AtomNode::Num(n) = exp.node()
        && let Some((a, _b)) = linear_form(ctx, base, var)
    {
        if *n == -1 {
            // ∫ (a*x + b)^(-1) dx = log(a*x + b) / a
            return ctx.mul(&[ctx.fun("log", &[base]), ctx.pow(a, ctx.num(-1))]);
        }
        // ∫ (a*x + b)^n dx = (a*x + b)^(n+1) / (a * (n+1))
        let new_exp = ctx.num(n + 1);
        let denom = ctx.mul(&[a, ctx.num(n + 1)]);
        return ctx.mul(&[ctx.pow(base, new_exp), ctx.pow(denom, ctx.num(-1))]);
    }

    // Fractional exponent on a linear form: (a*x + b)^(p/q).
    if let Some((p, q)) = fraction_exponent(exp)
        && p != -q
        && let Some((a, _b)) = linear_form(ctx, base, var)
    {
        // ∫ (a*x+b)^(p/q) dx = (a*x+b)^((p+q)/q) * q / (a*(p+q))
        let new_exp = ctx.mul(&[ctx.num(p + q), ctx.pow(ctx.num(q), ctx.num(-1))]);
        let coeff_num = ctx.num(q);
        let coeff_den = ctx.mul(&[a, ctx.num(p + q)]);
        return ctx.mul(&[
            ctx.pow(base, new_exp),
            coeff_num,
            ctx.pow(coeff_den, ctx.num(-1)),
        ]);
    }

    fallback(ctx, ctx.pow(base, exp), var)
}

/// Parse an exponent atom as a fraction p/q (small integers).
///
/// Accepts `p * q^-1`, `p * (q^-1)` with integer p, q (q > 0), as produced by
/// rational arithmetic in the ODE solvers, and the normalized bare
/// reciprocal `q^-1` (= 1/q).
fn fraction_exponent<'a>(exp: Atom<'a>) -> Option<(i64, i64)> {
    if let AtomNode::Pow(b, e) = exp.node()
        && let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
        && *ee == -1
        && *bb > 0
    {
        return Some((1, *bb));
    }
    if let AtomNode::Mul(args) = exp.node() {
        let mut num: Option<i64> = None;
        let mut den: Option<i64> = None;
        for a in args.iter() {
            match a.node() {
                AtomNode::Num(n) => {
                    if num.is_some() {
                        return None;
                    }
                    num = Some(*n);
                }
                AtomNode::Pow(b, e) => {
                    if let (AtomNode::Num(bb), AtomNode::Num(ee)) = (b.node(), e.node())
                        && *ee == -1
                    {
                        if den.is_some() {
                            return None;
                        }
                        den = Some(*bb);
                    } else {
                        return None;
                    }
                }
                _ => return None,
            }
        }
        if let (Some(p), Some(q)) = (num, den)
            && q > 0
        {
            return Some((p, q));
        }
    }
    None
}

/// If `expr` is of the form `a*x + b` (with `a` and `b` constant w.r.t. `var`),
/// return `(a, b)`. Otherwise return None.
pub(crate) fn linear_form<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Atom<'a>, Atom<'a>)> {
    match expr.node() {
        AtomNode::Var(v) if *v == var => Some((ctx.num(1), ctx.num(0))),
        AtomNode::Mul(args) => {
            let mut coeff = ctx.num(1);
            let mut has_var = false;
            for a in args.iter() {
                if let AtomNode::Var(v) = a.node()
                    && *v == var
                {
                    has_var = true;
                    continue;
                }
                if is_constant(*a, var) {
                    coeff = ctx.mul(&[coeff, *a]);
                } else {
                    return None;
                }
            }
            if has_var {
                Some((coeff, ctx.num(0)))
            } else {
                None
            }
        }
        AtomNode::Add(args) => {
            let mut a_part = ctx.num(0);
            let mut b_part = ctx.num(0);
            for arg in args.iter() {
                if let Some((ca, _cb)) = linear_form(ctx, *arg, var) {
                    a_part = ctx.add(&[a_part, ca]);
                } else if is_constant(*arg, var) {
                    b_part = ctx.add(&[b_part, *arg]);
                } else {
                    return None;
                }
            }
            Some((a_part, b_part))
        }
        _ => None,
    }
}

fn integrate_function<'a>(
    ctx: &'a AtomArena<'a>,
    name: Symbol,
    args: &'a [Atom<'a>],
    var: Symbol,
    _depth: usize,
) -> Atom<'a> {
    if args.is_empty() {
        return fallback(ctx, ctx.fun(name.as_str(), args), var);
    }
    let u = args[0];

    // Simple linear substitution forms: f(a*x + b)
    if let Some((a, _b)) = linear_form(ctx, u, var)
        && is_constant(a, var)
        && !is_one(a)
    {
        let inner_integral = match name.as_str() {
            "sin" => ctx.mul(&[ctx.num(-1), ctx.fun("cos", &[u])]),
            "cos" => ctx.fun("sin", &[u]),
            "exp" => ctx.fun("exp", &[u]),
            _ => return fallback(ctx, ctx.fun(name.as_str(), args), var),
        };
        return ctx.mul(&[ctx.pow(a, ctx.num(-1)), inner_integral]);
    }

    // Direct table for f(x) where u == x.
    if let AtomNode::Var(v) = u.node()
        && *v == var
    {
        let antiderivative: Option<Atom<'a>> = match name.as_str() {
            "sin" => Some(ctx.mul(&[ctx.num(-1), ctx.fun("cos", &[u])])),
            "cos" => Some(ctx.fun("sin", &[u])),
            "exp" => Some(ctx.fun("exp", &[u])),
            "log" => Some(ctx.mul(&[u, ctx.add(&[ctx.fun("log", &[u]), ctx.num(-1)])])),
            _ => None,
        };
        if let Some(anti) = antiderivative {
            return anti;
        }
    }

    fallback(ctx, ctx.fun(name.as_str(), args), var)
}

fn is_one<'a>(expr: Atom<'a>) -> bool {
    matches!(expr.node(), AtomNode::Num(1))
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn csc_product_expansion_not_starved_by_parts() {
        // Regression (0.27.1): the parts heuristic consumed the whole
        // chain-entry budget on this three-factor product, so the expansion
        // retry's per-term re-entries tripped the budget and fell back.
        // The expansion must run before the heuristic stage.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(
            &ctx,
            "csc(c + d*x)^3*(a - a*csc(c + d*x))*(A - A*csc(c + d*x))",
        )
        .unwrap();
        let r = integrate(&ctx, expr, Symbol::new("x"));
        assert!(
            !r.to_string().contains("Integral("),
            "csc product left a residue: {r}"
        );
    }

    #[test]
    fn integrate_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(3^-1)*(x^3)");
    }

    #[test]
    fn integrate_inverse() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(x, ctx.num(-1));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "log(x)");
    }

    #[test]
    fn integrate_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("sin", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "-1*(cos(x))");
    }

    #[test]
    fn integrate_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("cos", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "sin(x)");
    }

    #[test]
    fn integrate_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("exp", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "exp(x)");
    }

    #[test]
    fn integrate_linear_substitution() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let two_x_plus_one = ctx.add(&[ctx.mul(&[ctx.num(2), x]), ctx.num(1)]);
        let expr = ctx.pow(two_x_plus_one, ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "(6^-1)*((1 + (2*x))^3)");
    }

    #[test]
    fn integrate_unknown() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.fun("f", &[x]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Integral(f(x), x)");
    }

    #[test]
    fn integrate_sin_times_cos_via_trig_path() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ sin(x)·cos(x) dx — solved by the product-to-sum reduction
        // (0.27.x): sin(x)cos(x) = ½ sin(2x).
        let expr = ctx.mul(&[ctx.fun("sin", &[x]), ctx.fun("cos", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn integrate_cos_squared_via_trig_path() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ cos(x)² dx = x/2 + sin(2x)/4 — via power reduction (0.27.x).
        let expr = ctx.pow(ctx.fun("cos", &[x]), ctx.num(2));
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
    }

    #[test]
    fn integrate_exp_neg_x_squared_gives_erf() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ exp(-x²) dx = (√π/2)·erf(x) — the 0.11.0 known gap, now closed
        // by the special-function table.
        let expr = ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(result.to_string().contains("erf"), "got {result}");
        assert!(!result.to_string().starts_with("Integral"), "got {result}");
    }

    #[test]
    fn integrate_expand_product_retry() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ x*(x+1)^2 dx — the direct chain declines the product; the
        // bounded-expansion retry distributes and integrates termwise.
        let expr = ctx.mul(&[x, ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(2))]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
        // Verify by differentiation: diff(result) - integrand folds to 0.
        let d = crate::diff(&ctx, result, Symbol::new("x"));
        let residual = ctx.add(&[d, ctx.mul(&[ctx.num(-1), expr])]);
        let folded = crate::ode::util::collect_terms(&ctx, residual);
        assert_eq!(folded.to_string(), "0", "residual: {folded}");
    }

    #[test]
    fn integrate_expand_deep_product() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ (x+1)^2 * (x+2)^2 dx — double distribution, 9 expanded terms.
        let s1 = ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(2));
        let s2 = ctx.pow(ctx.add(&[x, ctx.num(2)]), ctx.num(2));
        let expr = ctx.mul(&[s1, s2]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!result.to_string().contains("Integral"), "got {result}");
        let d = crate::diff(&ctx, result, Symbol::new("x"));
        let residual = ctx.add(&[d, ctx.mul(&[ctx.num(-1), expr])]);
        let folded = crate::ode::util::collect_terms(&ctx, residual);
        assert_eq!(folded.to_string(), "0", "residual: {folded}");
    }

    #[test]
    fn integrate_expand_budget_keeps_fallback() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // (x+1)^8 * f(x): the power alone fits the budget, but the unknown
        // function factor still cannot integrate — the original fallback
        // form must be preserved (not the expanded one).
        let s = ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(8));
        let expr = ctx.mul(&[s, ctx.fun("f", &[x])]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Integral((f(x))*((1 + x)^8), x)");
    }

    #[test]
    fn integrate_weierstrass_cubed_denominator_terminates() {
        // 1/(a + b*cos(c + d*x))^3 — previously hung: the Weierstrass
        // t-integrand re-entered the ℚ rational backend, whose dense gcd
        // (naive pseudo-remainder) exploded coefficients at degree ~30.
        // With the subresultant gcd the case terminates (~0.1 s in debug);
        // a partial result with an honest Integral residue is acceptable.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "1/(-5 + 3*cos(c + d*x))^3").expect("parse");
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let _ = result; // termination is the assertion
    }

    #[test]
    fn integrate_parts_weierstrass_cycle_terminates() {
        // (c + d*x)^2/(a + a*sin(e + f*x)) — the parts ↔ Weierstrass cycle
        // (parts produces atan(_t)·T forms whose v' = T re-enters
        // Weierstrass with a fresh parts budget) looped until stack
        // overflow. The chain-entry budget caps it with an honest partial
        // result.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "(c + d*x)^2/(a + a*sin(e + f*x))").expect("parse");
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let _ = result; // termination is the assertion
    }

    #[test]
    fn integrate_exp_x_over_x_gives_ei() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // ∫ exp(x)/x dx = Ei(x) — non-elementary, special-function table.
        let expr = ctx.mul(&[ctx.fun("exp", &[x]), ctx.pow(x, ctx.num(-1))]);
        let result = integrate(&ctx, expr, Symbol::new("x"));
        assert_eq!(result.to_string(), "Ei(x)");
    }

    #[test]
    fn expanded_linear_square_is_folded_into_a_power() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ocas_parse::parse(&ctx, "1/(a^2 + 2*a*b*x + b^2*x^2)").unwrap();
        let folded = fold_linear_squares(&ctx, expr, Symbol::new("x"), true);
        assert!(
            folded.to_string().contains("(a + (b*x))^2"),
            "the trinomial should become a squared linear form: {folded}"
        );
        // The rewrite is accepted by exact re-expansion inside
        // `try_fold_linear_square`; confirm that canonicalisation agrees with
        // the parsed trinomial, which is what makes the acceptance test work.
        let lhs = crate::ode::util::collect_terms(
            &ctx,
            ctx.pow(
                ctx.add(&[ctx.var("a"), ctx.mul(&[ctx.var("b"), ctx.var("x")])]),
                ctx.num(2),
            ),
        );
        let rhs = crate::ode::util::collect_terms(
            &ctx,
            ocas_parse::parse(&ctx, "a^2 + 2*a*b*x + b^2*x^2").unwrap(),
        );
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn perfect_square_corpus_hang_is_solved() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // `rubi-00854`: the expanded trinomial leaves the symbolic-rational
        // backend grinding (a baseline 10 s timeout in 0.27.2); the exact fold
        // reduces it to a two-linear-factor partial fraction.
        let expr =
            ocas_parse::parse(&ctx, "(a + b*x)/((d + e*x)^4*(a^2 + 2*a*b*x + b^2*x^2))").unwrap();
        let r = integrate(&ctx, expr, Symbol::new("x"));
        assert!(!contains_integral(r), "not solved: {r}");
    }

    #[test]
    fn linear_square_fold_keeps_radicands_and_non_squares() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // A non-integer power over the trinomial is left alone: the radical and
        // elliptic engines match the expanded quadratic, and
        // `((a+b*x)^2)^(1/2) = |a+b*x|` is not `a+b*x`.
        for src in [
            "sqrt(a^2 + 2*a*b*x + b^2*x^2)",
            "(a^2 + 2*a*b*x + b^2*x^2)^(1/2)",
        ] {
            let radicand = ocas_parse::parse(&ctx, src).unwrap();
            assert_eq!(
                fold_linear_squares(&ctx, radicand, Symbol::new("x"), true),
                radicand,
                "{src}"
            );
        }
        // A three-term sum that is not a perfect square is untouched.
        let not_square = ocas_parse::parse(&ctx, "a^2 + 2*a*b*x + c*x^2").unwrap();
        assert_eq!(
            fold_linear_squares(&ctx, not_square, Symbol::new("x"), true),
            not_square
        );
    }
}
