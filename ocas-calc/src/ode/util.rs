//! ODE utility functions: coefficient extraction, linearity checks, order detection.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};

use crate::derivative::diff;

/// Return the order (highest derivative degree) of the ODE for unknown `func`.
///
/// Scans for `Derivative(func, var)`, `Derivative(func, var, var)`, etc.
/// Returns 0 if no derivatives of `func` appear.
pub(crate) fn ode_order<'a>(equation: Atom<'a>, func: Atom<'a>, var: Symbol) -> usize {
    ode_order_inner(equation, func, var)
}

fn ode_order_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> usize {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => 0,
        AtomNode::Add(args) | AtomNode::Mul(args) => args
            .iter()
            .map(|a| ode_order_inner(*a, func, var))
            .max()
            .unwrap_or(0),
        AtomNode::Pow(base, exp) => {
            ode_order_inner(*base, func, var).max(ode_order_inner(*exp, func, var))
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                // Check if the first arg is `func` and second arg is `var`.
                if args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str() {
                    // Order = number of extra `var` args beyond the first two.
                    // Derivative(y(x), x) -> order 1
                    // Derivative(y(x), x, x) -> order 2
                    args.len() - 1
                } else {
                    0
                }
            } else {
                args.iter()
                    .map(|a| ode_order_inner(*a, func, var))
                    .max()
                    .unwrap_or(0)
            }
        }
    }
}

/// Check if `expr` is linear in `func` and its derivatives.
///
/// Linear means: `func`, `Derivative(func, var)`, etc. appear only as
/// first-degree terms multiplied by expressions free of `func`.
pub(crate) fn is_linear_in<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    is_linear_inner(expr, func, var)
}

fn is_linear_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => true,
        AtomNode::Add(args) => args.iter().all(|a| is_linear_inner(*a, func, var)),
        AtomNode::Mul(args) => {
            // At most one factor may contain `func`.
            let func_dependent_count = args
                .iter()
                .filter(|a| contains_func(**a, func, var))
                .count();
            if func_dependent_count > 1 {
                return false;
            }
            // The factor containing func must itself be linear (first degree).
            args.iter().all(|a| {
                if contains_func(*a, func, var) {
                    is_func_first_degree(*a, func, var)
                } else {
                    true
                }
            })
        }
        AtomNode::Pow(base, exp) => {
            // If base contains func, exponent must be a positive constant == 1.
            if contains_func(*base, func, var) {
                if let AtomNode::Num(n) = exp.node() {
                    *n == 1
                } else {
                    false
                }
            } else {
                // exp may contain func only if it's first degree and base is constant.
                !contains_func(*exp, func, var)
                    || (is_func_first_degree(*exp, func, var) && !contains_func(*base, func, var))
            }
        }
        AtomNode::Fun(_, args) => {
            // A function call containing func is nonlinear (e.g. sin(y)).
            !contains_func(expr, func, var)
                || args.iter().all(|a| {
                    if contains_func(*a, func, var) {
                        is_func_first_degree(*a, func, var)
                    } else {
                        true
                    }
                })
        }
    }
}

/// Check if `expr` contains `func` or its derivatives.
pub(crate) fn contains_func<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    contains_func_inner(expr, func, var)
}

fn contains_func_inner<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => false,
        AtomNode::Add(args) | AtomNode::Mul(args) => {
            args.iter().any(|a| contains_func_inner(*a, func, var))
        }
        AtomNode::Pow(base, exp) => {
            contains_func_inner(*base, func, var) || contains_func_inner(*exp, func, var)
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str()
            } else {
                args.iter().any(|a| contains_func_inner(*a, func, var))
            }
        }
    }
}

/// Check if `func` or its derivatives appear only to the first power.
fn is_func_first_degree<'a>(expr: Atom<'a>, func: Atom<'a>, var: Symbol) -> bool {
    match expr.node() {
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str()
            } else {
                false
            }
        }
        AtomNode::Pow(base, exp) => {
            if contains_func(*base, func, var) {
                if let AtomNode::Num(n) = exp.node() {
                    *n == 1
                } else {
                    false
                }
            } else {
                !contains_func(*exp, func, var)
            }
        }
        AtomNode::Mul(args) => {
            let func_count: usize = args
                .iter()
                .filter(|a| contains_func(**a, func, var))
                .count();
            func_count <= 1
        }
        AtomNode::Add(args) => args.iter().all(|a| is_func_first_degree(*a, func, var)),
        _ => !contains_func(expr, func, var),
    }
}

/// Substitute `y(x) → sol_expr` in `equation`, also replacing
/// `Derivative(y(x), x)` with `diff(sol_expr, x)`, and higher derivatives
/// accordingly.
pub(crate) fn substitute_solution<'a>(
    ctx: &'a AtomArena<'a>,
    equation: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a> {
    substitute_inner(ctx, equation, func, sol, var)
}

fn substitute_inner<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    func: Atom<'a>,
    sol: Atom<'a>,
    var: Symbol,
) -> Atom<'a> {
    // The unknown function itself (e.g. bare y(x)) maps to the solution.
    if expr.to_string() == func.to_string() {
        return sol;
    }
    match expr.node() {
        AtomNode::Num(_) | AtomNode::Var(_) => expr,
        AtomNode::Add(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_inner(ctx, *a, func, sol, var))
                .collect();
            ctx.add(&mapped)
        }
        AtomNode::Mul(args) => {
            let mapped: Vec<_> = args
                .iter()
                .map(|a| substitute_inner(ctx, *a, func, sol, var))
                .collect();
            ctx.mul(&mapped)
        }
        AtomNode::Pow(base, exp) => {
            let b = substitute_inner(ctx, *base, func, sol, var);
            let e = substitute_inner(ctx, *exp, func, sol, var);
            ctx.pow(b, e)
        }
        AtomNode::Fun(name, args) => {
            if *name == Symbol::new("Derivative") && args.len() >= 2 {
                if args[0].to_string() == func.to_string() && args[1].to_string() == var.as_str() {
                    // Derivative(y(x), x, x, ...) -> differentiate sol n times.
                    let order = args.len() - 1;
                    let mut result = sol;
                    for _ in 0..order {
                        result = diff(ctx, result, var);
                    }
                    result
                } else {
                    let mapped: Vec<_> = args
                        .iter()
                        .map(|a| substitute_inner(ctx, *a, func, sol, var))
                        .collect();
                    ctx.fun(name.as_str(), &mapped)
                }
            } else if args.iter().any(|a| contains_func_inner(*a, func, var)) {
                // func appears inside a regular function call (e.g. sin(y))
                let mapped: Vec<_> = args
                    .iter()
                    .map(|a| substitute_inner(ctx, *a, func, sol, var))
                    .collect();
                ctx.fun(name.as_str(), &mapped)
            } else {
                expr
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Like-term collection
// ---------------------------------------------------------------------------

/// A reduced rational coefficient/exponent `num/den` with `den > 0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Rat {
    pub(crate) num: i128,
    pub(crate) den: i128,
}

impl Rat {
    fn new(num: i128, den: i128) -> Self {
        debug_assert!(den != 0);
        let (num, den) = if den < 0 { (-num, -den) } else { (num, den) };
        let g = gcd_i128(num.unsigned_abs(), den as u128);
        Rat {
            num: num / g as i128,
            den: den / g as i128,
        }
    }

    fn one() -> Self {
        Rat { num: 1, den: 1 }
    }

    fn is_zero(&self) -> bool {
        self.num == 0
    }

    fn add(self, other: Rat) -> Rat {
        Rat::new(
            self.num.saturating_mul(other.den) + other.num.saturating_mul(self.den),
            self.den.saturating_mul(other.den),
        )
    }

    fn mul(self, other: Rat) -> Rat {
        Rat::new(
            self.num.saturating_mul(other.num),
            self.den.saturating_mul(other.den),
        )
    }
}

fn gcd_i128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.max(1)
}

/// Decompose a multiplicative term into a rational coefficient and a list of
/// `(base, exponent)` pairs with rational exponents.
///
/// Handles `Num`, `Pow(Num, Num)` (folded into the coefficient), `Var`/`Fun`
/// bases (exponent 1), `Pow(base, Num)`, and rational exponent atoms such as
/// `Pow(base, Mul[Num(-1), Pow(Num(2), Num(-1))])` (i.e. $x^{-1/2}$).
/// Any other factor is kept as an opaque base with exponent 1.
fn decompose_term<'a>(term: Atom<'a>) -> (Rat, Vec<(Atom<'a>, Rat)>) {
    decompose_term_impl(term)
}

fn decompose_term_impl<'a>(term: Atom<'a>) -> (Rat, Vec<(Atom<'a>, Rat)>) {
    let mut coeff = Rat::one();
    let mut factors: Vec<(Atom<'a>, Rat)> = Vec::new();

    let mut push_factor = |base: Atom<'a>, e: Rat| {
        if e.is_zero() {
            return;
        }
        if let Some(entry) = factors
            .iter_mut()
            .find(|(b, _)| b.to_string() == base.to_string())
        {
            entry.1 = entry.1.add(e);
        } else {
            factors.push((base, e));
        }
    };

    let factors_src: Vec<Atom<'a>> = match term.node() {
        AtomNode::Mul(args) => args.to_vec(),
        _ => vec![term],
    };

    for f in factors_src {
        match f.node() {
            AtomNode::Num(n) => {
                coeff = coeff.mul(Rat::new(*n as i128, 1));
            }
            AtomNode::Pow(base, exp) => {
                if let (AtomNode::Num(b), AtomNode::Num(e)) = (base.node(), exp.node()) {
                    // Pure numeric power folds into the coefficient.
                    if *e >= 0 {
                        coeff = coeff.mul(Rat::new((*b as i128).saturating_pow(*e as u32), 1));
                    } else {
                        coeff = coeff.mul(Rat::new(1, (*b as i128).saturating_pow((-*e) as u32)));
                    }
                } else if let Some(r) = numeric_pow_rat(*base, *exp) {
                    // Nested numeric power such as (2^-1)^-1 = 2.
                    coeff = coeff.mul(r);
                } else if let Some(e_rat) = exponent_rat(*exp) {
                    // Flatten nested powers with a variable base:
                    // (u^p)^q = u^(p*q) when p is rational.
                    if let AtomNode::Pow(inner_base, inner_exp) = base.node()
                        && let Some(inner_rat) = exponent_rat(*inner_exp)
                    {
                        push_factor(*inner_base, inner_rat.mul(e_rat));
                    } else {
                        push_factor(*base, e_rat);
                    }
                } else {
                    push_factor(f, Rat::one());
                }
            }
            _ => push_factor(f, Rat::one()),
        }
    }

    // Drop factors whose accumulated exponent cancelled to zero, and fold
    // numeric bases with exactly-evaluable exponents (e.g. 3^1, 4^(1/2))
    // back into the coefficient so that, e.g., sqrt(3)*sqrt(3) = 3 becomes
    // a pure coefficient and matches like terms from other paths.
    let mut folded: Vec<(Atom<'a>, Rat)> = Vec::new();
    for (base, e) in factors {
        if e.is_zero() {
            continue;
        }
        if let AtomNode::Num(b) = base.node()
            && let Some(r) = int_pow_rat(*b, e)
        {
            coeff = coeff.mul(r);
            continue;
        }
        folded.push((base, e));
    }
    (coeff, folded)
}

/// Evaluate `b^e` exactly for integer b and rational e when the result is
/// rational: e integer, or e = p/q with b a perfect q-th power.
fn int_pow_rat(b: i64, e: Rat) -> Option<Rat> {
    if e.den == 1 {
        let k = e.num;
        return if k >= 0 {
            let ku = u32::try_from(k).ok()?;
            Some(Rat::new((b as i128).saturating_pow(ku), 1))
        } else {
            let ku = u32::try_from(k.checked_neg()?).ok()?;
            Some(Rat::new(1, (b as i128).saturating_pow(ku)))
        };
    }
    // e = p/q with q > 1: need b to be a perfect q-th power, then take
    // the p-th power of the root. Only q = 2 is handled (covers sqrt).
    if e.den == 2 && b >= 0 {
        let r = isqrt_i128(b as i128);
        if r * r == b as i128 {
            return int_pow_rat(r as i64, Rat::new(e.num, 1));
        }
    }
    None
}

/// Integer square root for i128 (floor), n >= 0.
fn isqrt_i128(n: i128) -> i128 {
    if n <= 0 {
        return 0;
    }
    let mut x = n;
    let mut y = x.div_euclid(2) + 1;
    while y < x {
        x = y;
        y = (x + n / x).div_euclid(2);
    }
    x
}

/// Parse an exponent atom into a rational number.
///
/// Accepts `Num(n)`, `Pow(Num(d), Num(-1))` (i.e. $1/d$), and flat `Mul`
/// combinations of those (e.g. $-1/2$ as `-1 * 2^-1`).
fn exponent_rat<'a>(exp: Atom<'a>) -> Option<Rat> {
    match exp.node() {
        AtomNode::Num(n) => Some(Rat::new(*n as i128, 1)),
        AtomNode::Pow(base, e) => {
            if let (AtomNode::Num(b), AtomNode::Num(k)) = (base.node(), e.node()) {
                if *k < 0 {
                    Some(Rat::new(1, (*b as i128).saturating_pow((-*k) as u32)))
                } else {
                    Some(Rat::new((*b as i128).saturating_pow(*k as u32), 1))
                }
            } else {
                None
            }
        }
        AtomNode::Mul(args) => {
            let mut acc = Rat::one();
            for a in args.iter() {
                acc = acc.mul(exponent_rat(*a)?);
            }
            Some(acc)
        }
        _ => None,
    }
}

/// Evaluate `base^exp` as a rational when both are rational-representable
/// atoms with an exact rational value (e.g. `(2^-1)^-1 = 2`,
/// `(2^2)^-1 = 1/4`). Returns `None` for non-exact cases like `2^(1/2)`.
fn numeric_pow_rat<'a>(base: Atom<'a>, exp: Atom<'a>) -> Option<Rat> {
    let b = exponent_rat(base)?;
    let e = exponent_rat(exp)?;
    // e must be an integer for an exact rational result in general.
    if e.den != 1 {
        return None;
    }
    let k = e.num;
    if k >= 0 {
        let ku = u32::try_from(k).ok()?;
        Some(Rat::new(b.num.saturating_pow(ku), b.den.saturating_pow(ku)))
    } else {
        let ku = u32::try_from(k.checked_neg()?).ok()?;
        Some(Rat::new(b.den.saturating_pow(ku), b.num.saturating_pow(ku)))
    }
}

/// Collect like terms of an additive expression.
///
/// Groups terms by their non-numeric factor signature (with rational
/// exponents combined, so $x \cdot x^{-1/2} = x^{1/2}$), sums the rational
/// coefficients of each group, and rebuilds the expression. Terms whose
/// coefficients cancel to zero are dropped; an empty result yields `Num(0)`.
///
/// This compensates for the general simplifier's lack of like-term
/// collection (e.g. `x^2 - 2*(2^-1)*x^2` becomes `0`).
pub(crate) fn collect_terms<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    let normalized = ocas_atom::normalize::normalize(ctx, expr);
    let expanded = expand(ctx, normalized);
    let normalized = ocas_atom::normalize::normalize(ctx, expanded);

    let terms: Vec<Atom<'a>> = match normalized.node() {
        AtomNode::Add(args) => args.to_vec(),
        _ => vec![normalized],
    };

    // Group coefficients by factor signature.
    type FactorGroup<'a> = (String, Vec<(Atom<'a>, Rat)>, Rat);
    let mut groups: Vec<FactorGroup<'a>> = Vec::new();
    for term in terms {
        let (coeff, factors) = decompose_term(term);
        let mut key_parts: Vec<String> = factors
            .iter()
            .map(|(b, e)| format!("{}^({}/{})", b, e.num, e.den))
            .collect();
        key_parts.sort();
        let key = key_parts.join("*");
        if let Some(group) = groups.iter_mut().find(|(k, _, _)| *k == key) {
            group.2 = group.2.add(coeff);
        } else {
            groups.push((key, factors, coeff));
        }
    }

    // Rebuild terms, dropping zero-coefficient groups.
    let mut rebuilt: Vec<Atom<'a>> = Vec::new();
    for (_, factors, coeff) in groups {
        if coeff.is_zero() {
            continue;
        }
        rebuilt.push(build_term(ctx, coeff, &factors));
    }

    let combined = if rebuilt.is_empty() {
        ctx.num(0)
    } else if rebuilt.len() == 1 {
        rebuilt[0]
    } else {
        ctx.add(&rebuilt)
    };
    ocas_atom::normalize::normalize(ctx, combined)
}

/// Simplify `exp(exponent)` when the exponent is a constant times a single
/// logarithm: `exp(k * log(u)) = u^k`. Falls back to the literal `exp` form.
///
/// Numeric constant factors inside the logarithm's argument are dropped:
/// they only scale the result by a constant, which is harmless in ODE
/// contexts (integrating factors, reduction of order).
pub(crate) fn exp_simplify<'a>(ctx: &'a AtomArena<'a>, exponent: Atom<'a>) -> Atom<'a> {
    let exponent = ocas_atom::normalize::normalize(ctx, exponent);
    let (coeff, factors) = decompose_term(exponent);
    if factors.len() == 1
        && let AtomNode::Fun(name, args) = factors[0].0.node()
        && name.as_str() == "log"
        && args.len() == 1
        && factors[0].1.num == 1
        && factors[0].1.den == 1
    {
        // Strip numeric constant factors from the log argument.
        let arg = args[0];
        let stripped = match arg.node() {
            AtomNode::Mul(fargs) => {
                let non_const: Vec<_> = fargs
                    .iter()
                    .filter(|a| !matches!(a.node(), AtomNode::Num(_)))
                    .copied()
                    .collect();
                match non_const.len() {
                    0 => arg,
                    1 => non_const[0],
                    _ => ctx.mul(&non_const),
                }
            }
            _ => arg,
        };
        return ctx.pow(stripped, rat_to_atom(ctx, coeff));
    }
    ctx.fun("exp", &[exponent])
}

/// Rebuild a term from a rational coefficient and `(base, exponent)` factors.
fn build_term<'a>(ctx: &'a AtomArena<'a>, coeff: Rat, factors: &[(Atom<'a>, Rat)]) -> Atom<'a> {
    let mut parts: Vec<Atom<'a>> = Vec::new();

    if coeff.num != 1 || factors.is_empty() {
        parts.push(ctx.num(coeff.num as i64));
    }
    if coeff.den != 1 {
        parts.push(ctx.pow(ctx.num(coeff.den as i64), ctx.num(-1)));
    }
    for (base, e) in factors {
        if e.num == 1 && e.den == 1 {
            parts.push(*base);
        } else {
            parts.push(ctx.pow(*base, rat_to_atom(ctx, *e)));
        }
    }

    match parts.len() {
        0 => ctx.num(1),
        1 => parts[0],
        _ => ctx.mul(&parts),
    }
}

/// Build an atom for a rational number `num/den` (as exponent).
fn rat_to_atom<'a>(ctx: &'a AtomArena<'a>, r: Rat) -> Atom<'a> {
    rat_to_atom_impl(ctx, r)
}

/// Distribute multiplication over addition, expanding `a*(b + c)` into
/// `a*b + a*c`. Only expands one `Mul` level deep into `Add` children —
/// powers of sums and nested products are handled recursively.
fn expand<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>) -> Atom<'a> {
    match expr.node() {
        AtomNode::Add(args) => {
            let mapped: Vec<_> = args.iter().map(|a| expand(ctx, *a)).collect();
            ctx.add(&mapped)
        }
        AtomNode::Mul(args) => {
            // Expand each factor first.
            let factors: Vec<Atom<'a>> = args.iter().map(|a| expand(ctx, *a)).collect();
            // Find the first Add factor to distribute over.
            if let Some(pos) = factors
                .iter()
                .position(|f| matches!(f.node(), AtomNode::Add(_)))
            {
                let add_args = match factors[pos].node() {
                    AtomNode::Add(a) => a,
                    _ => unreachable!(),
                };
                let mut result_terms = Vec::with_capacity(add_args.len());
                for term in add_args.iter() {
                    let mut new_factors: Vec<Atom<'a>> = factors.clone();
                    new_factors[pos] = *term;
                    result_terms.push(expand(ctx, ctx.mul(&new_factors)));
                }
                ctx.add(&result_terms)
            } else {
                ctx.mul(&factors)
            }
        }
        AtomNode::Pow(base, exp) => {
            // (a + b)^n with small positive integer n expands by repeated
            // multiplication.
            if let AtomNode::Num(n) = exp.node()
                && *n >= 2
                && *n <= 8
                && matches!(base.node(), AtomNode::Add(_))
            {
                let mut acc = *base;
                for _ in 1..*n {
                    acc = expand(ctx, ctx.mul(&[acc, *base]));
                }
                return acc;
            }
            expr
        }
        _ => expr,
    }
}

fn rat_to_atom_impl<'a>(ctx: &'a AtomArena<'a>, r: Rat) -> Atom<'a> {
    if r.den == 1 {
        ctx.num(r.num as i64)
    } else if r.num == 1 {
        ctx.pow(ctx.num(r.den as i64), ctx.num(-1))
    } else {
        ctx.mul(&[
            ctx.num(r.num as i64),
            ctx.pow(ctx.num(r.den as i64), ctx.num(-1)),
        ])
    }
}
