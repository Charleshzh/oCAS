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

/// `rubi-00401` branch-shape pin for the 0.30.0 W6a single-radical stage.
///
/// The 0.30.0 Wave 6a single-radical atomic-layer reduction
/// (`ocas-calc/src/integral/single_radical.rs`) shipped **unwired**: with the
/// stage in the pipeline, `rubi-00401` got a *wrong* answer — the `x = 1/t`
/// substitution mishandled the sign/branch factor of a quotient radicand, so
/// the emitted antiderivative was off by a global sign (oracle relative
/// error 2.0). The wiring was reverted and the branch bug is the first item
/// on the 0.31 W6a list.
///
/// 0.31.0 fixed it with the branch-safe even substitution plus a two-sheet
/// numeric gate on every emission, and wired the stage back. The shape is now
/// **solved and verified on both sheets** (`x = ±0.6` agree exactly), so this
/// pin is the hard gate the 0.30.0 note asked for: a fallback or an
/// unverified answer is a regression, and a wrong closed form is the failure
/// mode the numerical oracle exists to catch. Do not delete or loosen this
/// assertion to make a wiring experiment green.
#[test]
fn wave_030_w6a_branch_shape_never_wrong() {
    let input = "sqrt(-1 + 1/x^2)/(x*(-1 + x^2)^2)";
    let (text, verdict) = integrate_and_verify(input);
    match verdict {
        Verify::Verified { checked, .. } => {
            eprintln!("rubi-00401: verified by {checked} sample point(s): {text}");
        }
        Verify::Mismatch { detail, .. } => {
            panic!("rubi-00401 wrong-antiderivative regression: {input} -> {text}\n    {detail}")
        }
        Verify::Indeterminate { reason } => {
            panic!(
                "rubi-00401 not verified after the 0.31.0 W6a wiring ({reason}): {input} -> {text}"
            )
        }
    }
}

/// The 0.27.3 mechanisms must not introduce wrong answers. Each entry is a
/// shape a new mechanism claims:
///
/// - the exact algebraic fold `p² + 2·p·q + q² → (p+q)²` for an affine base
///   (`rubi-00854` ground truth: the unfolded trinomial is a corpus hang),
/// - the special-function families (`erf`, `Ci`, `Chi`, `Shi`, `Ei`,
///   `Ei(n, z)`, `erfi`),
/// - the affine-argument half-power front-end (`cos(c + d·x)`, previously a
///   documented gap).
///
/// Declining a case is fine (the corpus harness counts it as a fallback); what
/// this test forbids is a *wrong* closed form. `Indeterminate` is reported by
/// the harness but not treated as failure here, so the assertion is
/// deliberately limited to `Mismatch`.
///
/// The log-of-a-linear-fraction shapes (`rubi-01014`, `rubi-01672`) are
/// deliberately *not* sampled: they are the two 0.27.2 per-case timeouts, and
/// integrating them costs minutes even in debug — a guard should not carry a
/// known hang.
#[test]
fn wave_0273_mechanisms_have_no_wrong_answers() {
    let sample = [
        // perfect-square trinomial with an affine base
        "1/(a^2 + 2*a*b*x + b^2*x^2)",
        "(a + b*x)/((d + e*x)^4*(a^2 + 2*a*b*x + b^2*x^2))",
        "(a + b*x)/(a^2 + 2*a*b*x + b^2*x^2)^2",
        // special-function families
        "x^4*erf(b*x)",
        "Ci(b*x)",
        "x*Chi(b*x)",
        "x^3*Shi(a + b*x)",
        "Ei(b*x)/x^4",
        "Ei(1, a + b*x)",
        "erfi(b*x)^2",
        // affine-argument half powers
        "1/cos(a + b*x)^(1/2)",
        "1/(b*cos(c + d*x))^(7/2)",
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
        "wrong antiderivatives in 0.27.3 mechanisms:\n{}",
        wrong.join("\n")
    );
}

/// 0.31.0 polylog/Li₂ wave (`POLYLOG_SPEC_CN.md`).
///
/// The two exp-log Rubi targets (`rubi-00746`, `rubi-01004`) are the
/// quantified gate evidence for the `exp-log ≥ 10/83` line: both must be
/// **solved** (no `Integral(` residue). Every shape of the wave must either
/// verify numerically or be honestly undecided, never a wrong closed form —
/// the A-class answers carry `Li₂(u)` *and* `Li₂(1−u)`, so their real domain
/// can be narrow and an `Indeterminate` verdict is a legitimate outcome there.
#[test]
fn wave_031_polylog_shapes_are_solved_and_never_wrong() {
    let shapes = [
        // B class: log-rational integrals whose antiderivative needs Li₂.
        "log(1 - x)/x",
        "(a + b*log(-2 + e*x))/x",                      // rubi-00746
        "(A + B*log(e*(a + b*x)/(c + d*x)))/(f + g*x)", // rubi-01004
        // A class: the integrand already carries a dilogarithm.
        "polylog(2, c*x)",
        "x*polylog(2, c*x)",
        "polylog(2, c*(a + b*x))/x^2", // rubi-00050
    ];
    let mut wrong = Vec::new();
    let mut unsolved = Vec::new();
    let mut verified = 0usize;
    for input in shapes {
        let (text, verdict) = integrate_and_verify(input);
        if text.contains("Integral(") {
            unsolved.push(format!("{input} -> {text}"));
        }
        match verdict {
            Verify::Verified { .. } => verified += 1,
            Verify::Mismatch { detail, .. } => {
                wrong.push(format!("{input} -> {text}\n    {detail}"))
            }
            Verify::Indeterminate { .. } => {}
        }
    }
    assert!(
        wrong.is_empty(),
        "wrong antiderivatives in the 0.31 polylog wave:\n{}",
        wrong.join("\n")
    );
    assert!(
        unsolved.is_empty(),
        "0.31 polylog shapes that fell back:\n{}",
        unsolved.join("\n")
    );
    assert!(
        verified >= 2,
        "the oracle verified only {verified} of the 0.31 polylog shapes"
    );
}

/// The A/B target shapes from the audited 1892 corpus must stay solved: this
/// is the `exp-log` gate's per-case evidence, kept next to the mechanism.
#[test]
fn wave_031_polylog_corpus_targets_are_solved() {
    for input in [
        "(a + b*log(-2 + e*x))/x",
        "(A + B*log(e*(a + b*x)/(c + d*x)))/(f + g*x)",
        "polylog(2, c*(a + b*x))/x^2",
    ] {
        let (text, verdict) = integrate_and_verify(input);
        assert!(
            !text.contains("Integral("),
            "0.31 corpus target unsolved: {input} -> {text}"
        );
        assert!(
            !matches!(verdict, Verify::Mismatch { .. }),
            "0.31 corpus target is a wrong answer: {input} -> {text}"
        );
    }
}
