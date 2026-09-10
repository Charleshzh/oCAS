//! Product-to-sum reduction for products of `sin`/`cos` with linear
//! arguments.
//!
//! The rule table matches whole-expression head patterns, so composite
//! products like `(c+d*x)^2*cos(a+b*x)*sin(a+b*x)^2` fall through every
//! stage. This pass rewrites the trigonometric factor product into a sum
//! of single `sin`/`cos` terms at multiple angles (power reduction and
//! product-to-sum identities), which the chain then integrates termwise
//! after distribution. All arguments must be linear in the integration
//! variable (symbolic coefficients allowed), so a multiple angle `k*u` is
//! still a linear argument.

use ocas_atom::Symbol;
use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, AtomNode};

use super::linear_form;

/// Cap on the total number of trig factors (powers expanded): the reduction
/// tree branches by two per factor, so this bounds the output at 2^k leaves.
const MAX_TRIG_FACTORS: usize = 8;
/// Cap on produced sum terms; the expansion retry downstream has its own
/// budget on top of this.
const MAX_OUTPUT_TERMS: usize = 64;

/// A `sin(a*x + b)` / `cos(a*x + b)` factor with the linear argument kept
/// decomposed, so multiple-angle arithmetic is exact atom arithmetic.
#[derive(Clone, Copy)]
struct TrigLin<'a> {
    sin: bool,
    a: Atom<'a>,
    b: Atom<'a>,
}

/// Rewrite trig factor products in `expr` into sum form. Returns `None`
/// when no reducible product is present (fewer than two trig factors,
/// non-linear arguments, or unsupported functions), or when the budgets
/// are exceeded.
pub(crate) fn trig_reduce_products<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<Atom<'a>> {
    let (rest, trig) = extract_trig_factors(ctx, expr, var)?;
    if trig.len() < 2 || trig.len() > MAX_TRIG_FACTORS {
        return None;
    }
    let mut leaves: Vec<(i64, u32, Option<TrigLin<'a>>)> = Vec::new();
    reduce(ctx, &trig, 1, 0, &mut leaves)?;
    if leaves.len() > MAX_OUTPUT_TERMS {
        return None;
    }

    let mut terms = Vec::with_capacity(leaves.len());
    for (sign, two_pow, factor) in leaves {
        let mut parts: Vec<Atom<'a>> = rest.clone();
        if sign < 0 {
            parts.push(ctx.num(-1));
        }
        if two_pow > 0 {
            parts.push(ctx.pow(ctx.num(2), ctx.num(-(two_pow as i64))));
        }
        if let Some(f) = factor {
            parts.push(build_trig(ctx, f, var));
        }
        terms.push(ctx.mul(&parts));
    }
    let sum = ctx.add(&terms);
    if sum == expr { None } else { Some(sum) }
}

/// Classify one factor: `sin(u)^k` / `cos(u)^k` push `k` copies onto `trig`
/// (with `u` required linear in `var`); anything else lands in `rest`.
fn classify_factor<'a>(
    ctx: &'a AtomArena<'a>,
    f: Atom<'a>,
    var: Symbol,
    rest: &mut Vec<Atom<'a>>,
    trig: &mut Vec<TrigLin<'a>>,
) -> Option<()> {
    let (name, arg, copies) = match f.node() {
        AtomNode::Fun(name, args) if args.len() == 1 => (name, args[0], 1usize),
        AtomNode::Pow(b, e) => match (b.node(), e.node()) {
            (AtomNode::Fun(name, args), AtomNode::Num(k))
                if args.len() == 1 && (1..=8i64).contains(k) =>
            {
                (name, args[0], *k as usize)
            }
            _ => {
                rest.push(f);
                return Some(());
            }
        },
        _ => {
            rest.push(f);
            return Some(());
        }
    };
    if name.as_str() != "sin" && name.as_str() != "cos" {
        rest.push(f);
        return Some(());
    }
    let (a, b) = linear_form(ctx, arg, var)?;
    for _ in 0..copies {
        trig.push(TrigLin {
            sin: name.as_str() == "sin",
            a,
            b,
        });
    }
    Some(())
}

/// Split `expr` into non-trig factors and trig factors; `sin(u)^k` /
/// `cos(u)^k` contribute `k` copies. Returns `None` when any trig argument
/// is not linear in `var`.
fn extract_trig_factors<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
) -> Option<(Vec<Atom<'a>>, Vec<TrigLin<'a>>)> {
    let mut rest = Vec::new();
    let mut trig = Vec::new();
    match expr.node() {
        AtomNode::Mul(args) => {
            for a in args.iter() {
                classify_factor(ctx, *a, var, &mut rest, &mut trig)?;
            }
        }
        _ => classify_factor(ctx, expr, var, &mut rest, &mut trig)?,
    }
    if trig.is_empty() {
        return None;
    }
    Some((rest, trig))
}

/// Reduce the factor list by applying a product-to-sum identity to the
/// first two factors, recursing on both branches. Accumulates the sign and
/// the power-of-two denominator; leaves hold at most one factor.
fn reduce<'a>(
    ctx: &'a AtomArena<'a>,
    fs: &[TrigLin<'a>],
    sign: i64,
    two_pow: u32,
    out: &mut Vec<(i64, u32, Option<TrigLin<'a>>)>,
) -> Option<()> {
    if out.len() > MAX_OUTPUT_TERMS {
        return None;
    }
    if fs.len() <= 1 {
        out.push((sign, two_pow, fs.first().copied()));
        return Some(());
    }
    let (f0, f1) = (fs[0], fs[1]);
    let rest = &fs[2..];
    // Identities (branch factor, branch sign):
    //   sin u sin v = ½cos(u−v) − ½cos(u+v)
    //   cos u cos v = ½cos(u−v) + ½cos(u+v)
    //   sin u cos v = ½sin(u+v) + ½sin(u−v)
    let (u, v) = if !f0.sin && f1.sin {
        (f1, f0)
    } else {
        (f0, f1)
    };
    let branches: [(TrigLin<'a>, i64); 2] = if u.sin && v.sin {
        [
            (
                TrigLin {
                    sin: false,
                    a: sub(ctx, u.a, v.a),
                    b: sub(ctx, u.b, v.b),
                },
                1,
            ),
            (
                TrigLin {
                    sin: false,
                    a: add(ctx, u.a, v.a),
                    b: add(ctx, u.b, v.b),
                },
                -1,
            ),
        ]
    } else if !u.sin && !v.sin {
        [
            (
                TrigLin {
                    sin: false,
                    a: sub(ctx, u.a, v.a),
                    b: sub(ctx, u.b, v.b),
                },
                1,
            ),
            (
                TrigLin {
                    sin: false,
                    a: add(ctx, u.a, v.a),
                    b: add(ctx, u.b, v.b),
                },
                1,
            ),
        ]
    } else {
        // sin u cos v (after the swap above `u` is the sine factor).
        [
            (
                TrigLin {
                    sin: true,
                    a: add(ctx, u.a, v.a),
                    b: add(ctx, u.b, v.b),
                },
                1,
            ),
            (
                TrigLin {
                    sin: true,
                    a: sub(ctx, u.a, v.a),
                    b: sub(ctx, u.b, v.b),
                },
                1,
            ),
        ]
    };
    for (g, branch_sign) in branches {
        // sin(0) = 0 drops the branch; cos(0) = 1 drops the factor.
        let zero_angle = is_zero(g.a) && is_zero(g.b);
        if zero_angle && g.sin {
            continue;
        }
        let mut next = Vec::with_capacity(rest.len() + 1);
        next.extend_from_slice(rest);
        if !zero_angle {
            next.push(g);
        }
        reduce(ctx, &next, sign * branch_sign, two_pow + 1, out)?;
    }
    Some(())
}

/// Rebuild `sin(a*x + b)` / `cos(a*x + b)` from the decomposed argument.
/// Difference angles can produce a numerically negative leading
/// coefficient; canonicalize `cos(-u) → cos(u)` and `sin(-u) → -sin(u)` so
/// equal angles merge downstream.
fn build_trig<'a>(ctx: &'a AtomArena<'a>, f: TrigLin<'a>, var: Symbol) -> Atom<'a> {
    let arg = build_linear(ctx, f.a, f.b, var);
    if let AtomNode::Num(n) = f.a.node()
        && *n < 0
    {
        let pos_arg = build_linear(ctx, ctx.num(-n), ctx.mul(&[ctx.num(-1), f.b]), var);
        let fun = ctx.fun(if f.sin { "sin" } else { "cos" }, &[pos_arg]);
        return if f.sin {
            ctx.mul(&[ctx.num(-1), fun])
        } else {
            fun
        };
    }
    ctx.fun(if f.sin { "sin" } else { "cos" }, &[arg])
}

/// Rebuild `a*x + b`, dropping zero parts.
fn build_linear<'a>(ctx: &'a AtomArena<'a>, a: Atom<'a>, b: Atom<'a>, var: Symbol) -> Atom<'a> {
    let x = ctx.var(var.as_str());
    let a = fold(ctx, a);
    let b = fold(ctx, b);
    let ax = if is_one(a) {
        x
    } else if is_zero(a) {
        ctx.num(0)
    } else {
        ctx.mul(&[a, x])
    };
    let sum = if is_zero(b) { ax } else { ctx.add(&[b, ax]) };
    fold(ctx, sum)
}

fn add<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, v: Atom<'a>) -> Atom<'a> {
    fold(ctx, ctx.add(&[u, v]))
}

/// `u − v` with like-term folding; structurally identical atoms cancel to 0.
///
/// The folding step is essential, not cosmetic: `normalize` alone does not
/// collect `d + (−1)·d`, so two *mathematically* equal slopes coming from
/// different spellings would produce a non-zero-looking slope atom. A later
/// stage would then integrate `cos(slope·x + …)` as `sin(…)/slope` and emit a
/// coefficient dividing by an identically zero expression — the resonant
/// product-to-sum defect (`(a·cos(c+d·x) + b·sin(c+d·x))^2` and friends).
fn sub<'a>(ctx: &'a AtomArena<'a>, u: Atom<'a>, v: Atom<'a>) -> Atom<'a> {
    if u == v {
        return ctx.num(0);
    }
    let folded = fold(ctx, ctx.add(&[u, ctx.mul(&[ctx.num(-1), v])]));
    if is_zero(folded) { ctx.num(0) } else { folded }
}

/// `normalize` plus like-term collection, so equal slopes/phases fold to `0`
/// instead of staying as `d + (−1)·d`.
fn fold<'a>(ctx: &'a AtomArena<'a>, e: Atom<'a>) -> Atom<'a> {
    normalize(ctx, crate::ode::util::collect_terms(ctx, e))
}

fn is_zero(e: Atom<'_>) -> bool {
    matches!(e.node(), AtomNode::Num(0))
}

fn is_one(e: Atom<'_>) -> bool {
    matches!(e.node(), AtomNode::Num(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integral::integrate;
    use ocas_core::arena::Arena;

    /// Numeric f64 evaluator for test verification (handles the operators
    /// the reduction tests produce).
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
                    "sec" => v.cos().recip(),
                    "csc" => v.sin().recip(),
                    "cot" => v.tan().recip(),
                    "exp" => v.exp(),
                    "log" => v.ln(),
                    "sqrt" => v.sqrt(),
                    _ => return None,
                })
            }
        }
    }

    /// Integrate, require no residue, and check `diff(result) == integrand`
    /// numerically at the given sample points.
    fn assert_antiderivative_num<'a>(
        ctx: &'a AtomArena<'a>,
        integrand: Atom<'a>,
        var: Symbol,
        consts: &[(Symbol, f64)],
        samples: &[f64],
    ) {
        let result = integrate(ctx, integrand, var);
        assert!(
            !result.to_string().contains("Integral"),
            "fallback: {result}"
        );
        let d = crate::diff(ctx, result, var);
        for &xv in samples {
            let mut env = consts.to_vec();
            env.push((var, xv));
            let lhs = eval_f64(d, &env).expect("eval diff");
            let rhs = eval_f64(integrand, &env).expect("eval integrand");
            let tol = 1e-6 * rhs.abs().max(1.0);
            assert!(
                (lhs - rhs).abs() < tol,
                "at x={xv}: diff={lhs} integrand={rhs} (result: {result})"
            );
        }
    }

    #[test]
    fn sin_times_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.mul(&[ctx.fun("sin", &[x]), ctx.fun("cos", &[x])]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn cos_squared() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.pow(ctx.fun("cos", &[x]), ctx.num(2));
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn sin_squared_times_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let s2 = ctx.pow(ctx.fun("sin", &[x]), ctx.num(2));
        let expr = ctx.mul(&[s2, ctx.fun("cos", &[x])]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn poly_times_trig_product_symbolic_linear_arg() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let a = ctx.var("a");
        let b = ctx.var("b");
        // cos(a + b*x) * sin(a + b*x)^2 — symbolic linear argument.
        let u = ctx.add(&[a, ctx.mul(&[b, x])]);
        let s2 = ctx.pow(ctx.fun("sin", &[u]), ctx.num(2));
        let expr = ctx.mul(&[ctx.fun("cos", &[u]), s2]);
        let env = [(Symbol::new("a"), 0.4), (Symbol::new("b"), 1.3)];
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn poly_factor_times_reduced_trig() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // (2*x + 1)*cos(x)*sin(x)^2 — corpus shape class.
        let poly = ctx.add(&[ctx.mul(&[ctx.num(2), x]), ctx.num(1)]);
        let s2 = ctx.pow(ctx.fun("sin", &[x]), ctx.num(2));
        let expr = ctx.mul(&[poly, ctx.fun("cos", &[x]), s2]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn different_linear_args_product_to_sum() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // sin(x) * cos(2*x) = ½ sin(3x) − ½ sin(x)
        let expr = ctx.mul(&[
            ctx.fun("sin", &[x]),
            ctx.fun("cos", &[ctx.mul(&[ctx.num(2), x])]),
        ]);
        assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &[], &[0.3, 0.7, 1.1]);
    }

    #[test]
    fn declines_nonlinear_argument() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let x2 = ctx.pow(x, ctx.num(2));
        let expr = ctx.mul(&[ctx.fun("sin", &[x2]), ctx.fun("cos", &[x2])]);
        assert!(trig_reduce_products(&ctx, expr, Symbol::new("x")).is_none());
    }

    /// Resonance: both factors share the *same* linear argument, so the
    /// difference angle is identically zero. Before the folding fix the
    /// reduction left a slope atom `d + (−1)·d` behind and the chain emitted
    /// `1/(2·(d + (−1)·d))` — infinite at every parameter value (0.27.2
    /// corpus class: `(a·cos(c+d·x) + b·sin(c+d·x))^2`, `cos(u)^2·…`,
    /// `sin(u)·cos(u)·…`).
    #[test]
    fn resonant_equal_slopes_stay_finite_and_correct() {
        let inputs = [
            "sin(c + d*x)*cos(c + d*x)",
            "cos(c + d*x)^2",
            "sin(c + d*x)^2",
            "sin(c + d*x)^2*(a + b*sin(c + d*x)^2)",
            "cos(c + d*x)^3*(a + a*sin(c + d*x))",
            "(a*cos(c + d*x) + b*sin(c + d*x))^2",
        ];
        for input in inputs {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = ocas_parse::parse(&ctx, input).expect("parse");
            let result = integrate(&ctx, expr, Symbol::new("x"));
            let text = result.to_string();
            assert!(!text.contains("Integral"), "{input}: fallback {text}");
            // No coefficient may divide by an identically zero expression.
            assert!(
                !text.contains("(-1*d)") || !text.contains("d + (-1*d)"),
                "{input}: resonant denominator left behind: {text}"
            );
            let env = [
                (Symbol::new("a"), 1.5),
                (Symbol::new("b"), 2.0),
                (Symbol::new("c"), 0.5),
                (Symbol::new("d"), 1.5),
            ];
            assert_antiderivative_num(&ctx, expr, Symbol::new("x"), &env, &[0.3, 0.7, 1.9]);
        }
    }
}
