//! Numerical regression guard for integration wrong answers.
//!
//! `integrate` results are counted as "solved" by the corpus harness when no
//! `Integral(...)` residue remains. That string criterion accepted two
//! wrong-answer classes in 0.27.0 that were only found by per-case
//! inspection and fixed in 0.27.1:
//!
//! - the C14/D7b linear-argument power reductions divided the residual
//!   `Integral(g, x)` coefficient by the slope a second time
//!   (`sin(2*x+1)^3` and eight siblings),
//! - `rational::sqrt_positive_rational` dropped a `1/q` factor from
//!   `sqrt(p/q)`, silently returning wrong answers for shapes such as
//!   `1/(3*x^2+4*x+3)`.
//!
//! This module re-verifies both classes (plus a broader sample of 0.27.x
//! mechanisms) through the independent numerical oracle in
//! [`ocas_tests::integral_eval`], so a regression fails here rather than
//! silently inflating the coverage number.

use ocas::prelude::*;
use ocas_core::arena::Arena;
use ocas_tests::integral_eval::{Verify, verify_antiderivative};

/// Integrate `input` and numerically verify the result.
///
/// Returns the oracle verdict. A `Mismatch` means the engine returned a wrong
/// antiderivative; `Indeterminate` means the oracle could not decide (which is
/// reported, never treated as success).
fn integrate_and_verify(input: &str) -> (String, Verify) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, input).expect("parse integrand");
    let result = integrate(&ctx, expr, Symbol::new("x"));
    let text = result.to_string();
    let verdict = verify_antiderivative(expr, result, Symbol::new("x"));
    (text, verdict)
}

/// Every listed integrand must either verify or be honestly undecided — but
/// never be a wrong answer.
#[test]
fn historical_wrong_answer_classes_are_not_reintroduced() {
    // C14/D7b: the nine linear-argument power reductions whose residual
    // coefficient was divided by the slope twice in 0.27.0.
    let linear_arg_powers = [
        "sin(2*x+1)^3",
        "cos(3*x+1)^4",
        "tan(2*x+1)^3",
        "sec(2*x+1)^4",
        "csc(3*x+1)^3",
        "cot(2*x+1)^3",
        "sinh(2*x+1)^3",
        "cosh(3*x+1)^4",
        "tanh(2*x+1)^3",
    ];
    // sqrt(p/q) handling: the rational partial-fraction path.
    let rational_quadratics = [
        "1/(3*x^2+4*x+3)",
        "1/(2*x^2+3*x+5)",
        "1/(5*x^2+x+7)",
        "(2*x+1)/(3*x^2+4*x+3)",
    ];

    let mut wrong = Vec::new();
    let mut verified = 0usize;
    for input in linear_arg_powers.iter().chain(rational_quadratics.iter()) {
        let (text, verdict) = integrate_and_verify(input);
        match verdict {
            Verify::Verified { .. } => verified += 1,
            Verify::Mismatch { detail, .. } => {
                wrong.push(format!("{input} -> {text}\n    {detail}"))
            }
            Verify::Indeterminate { reason } => {
                // An undecided case is not a failure, but the result must at
                // least not be an unevaluated residue for the shapes above.
                assert!(
                    !text.contains("Integral("),
                    "{input}: engine fell back ({reason}) where 0.27.1 solved it: {text}"
                );
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "wrong antiderivatives recovered for {} case(s):\n{}",
        wrong.len(),
        wrong.join("\n")
    );
    // The oracle must actually have checked something.
    assert!(
        verified >= 8,
        "only {verified} of {} cases were numerically verified",
        linear_arg_powers.len() + rational_quadratics.len()
    );
}

/// A broad sample of 0.27.x mechanism families must not regress into wrong
/// answers. Cases the engine declines are fine (the corpus harness counts
/// them as fallback); cases it claims to solve must be right.
#[test]
fn mechanism_sample_has_no_wrong_answers() {
    let sample = [
        // rational / symbolic rational
        "1/(a+b*x)",
        "1/((a+b*x)^2)",
        "(A+B*x)/x^2",
        // exp / log kernels
        "exp(x)/(1-exp(2*x))",
        "sinh(x)/(2+3*cosh(x))",
        // trig kernels
        "sin(x)^2*cos(x)^5",
        "1/(2+cos(x))",
        "tan(x)^3",
        // inverse functions
        "atan(x)",
        "x*log(x)",
        // radicals
        "sqrt(1-x^2)",
        "1/sqrt(1-x^2)",
        // hyperbolic
        "1/(a+b*sinh(x))",
    ];
    let mut wrong = Vec::new();
    for input in sample {
        let (text, verdict) = integrate_and_verify(input);
        if let Verify::Mismatch { detail, .. } = verdict {
            wrong.push(format!("{input} -> {text}\n    {detail}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "wrong antiderivatives in mechanism sample:\n{}",
        wrong.join("\n")
    );
}
