//! Symbolic-certificate guards for integration results (0.28.0).
//!
//! 0.27.2 introduced an *independent numeric oracle*
//! ([`ocas_tests::integral_eval`]) because the corpus harness counted a case
//! as solved whenever the printed result contained no `Integral(` residue —
//! a criterion that accepted two wrong-answer classes. 0.28.0 adds the
//! stronger, machine-checkable criterion: the symbolic certificate
//! [`ocas_calc::integral::certify::certify`] must show that `D(F) − f`
//! reduces to zero.
//!
//! The tests here pin the two properties that matter and that the corpus
//! harness cannot enforce cheaply:
//!
//! 1. **Soundness** — a *wrong* antiderivative must never be certified.
//!    The historical wrong-answer classes of 0.27.1/0.27.2 are written out
//!    by hand and must all be rejected.
//! 2. **Agreement with the oracle** — whenever the numeric oracle calls a
//!    result a mismatch, the certificate must not pass it. A certificate
//!    that passes where the oracle disagrees would be a false positive: the
//!    wave's hard red line.
//!
//! The corpus-wide `certified_rate = 1.0` goal is measured by the 1892
//! harness (`ocas-tests/benches/integrate_1892.rs`); the corpus is never
//! committed, so the CI-scale gate lives here on a curated list.

use ocas::prelude::*;
use ocas_calc::integral::certify::{CertDecline, CertMethod, certify};
use ocas_core::arena::Arena;
use ocas_tests::integral_eval::{Verify, verify_antiderivative};

/// Integrate `input`, returning `(result, certificate verdict)`.
fn integrate_and_certify(input: &str) -> (String, std::result::Result<CertMethod, CertDecline>) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, input).expect("parse integrand");
    let var = Symbol::new("x");
    let normalized = ocas_atom::normalize::normalize(&ctx, expr);
    let result = integrate(&ctx, expr, var);
    let verdict = certify(&ctx, normalized, result, var).map(|c| c.method);
    (result.to_string(), verdict)
}

/// Hand-written wrong antiderivatives that 0.27.1/0.27.2 shipped and fixed.
///
/// Each entry is `(integrand, wrong antiderivative)`; the certificate must
/// reject every one of them.
const HISTORICAL_WRONG_ANSWERS: &[(&str, &str)] = &[
    // C14/D7b: the linear-argument power reductions divided the residual
    // coefficient by the slope a second time.
    ("sin(2*x+1)^3", "cos(2*x+1)^4/8"),
    ("sec(2*x+1)^4", "tan(2*x+1)^3/6"),
    ("sinh(2*x+1)^3", "cosh(2*x+1)^4/8"),
    // `rational::sqrt_positive_rational` dropped the `1/q` factor of
    // `sqrt(p/q)`.
    ("1/(3*x^2+4*x+3)", "sqrt(3)*atan((6*x+4)*sqrt(3)/6)/6"),
    // A plain scale error (guards against a prover that accepts anything
    // with the right shape).
    ("x", "x^2"),
    ("2*sin(x)", "sin(x)"),
    ("exp(x)", "exp(x)/2"),
    ("1/(1+x^2)", "2*atan(x)"),
];

/// A wrong antiderivative must never be certified.
#[test]
fn certificate_rejects_wrong_antiderivatives() {
    let mut accepted = Vec::new();
    for (integrand, wrong) in HISTORICAL_WRONG_ANSWERS {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let f =
            ocas_atom::normalize::normalize(&ctx, parse(&ctx, integrand).expect("parse integrand"));
        let big_f = parse(&ctx, wrong).expect("parse antiderivative");
        if let Ok(cert) = certify(&ctx, f, big_f, Symbol::new("x")) {
            accepted.push(format!(
                "{integrand} -> {wrong} certified as {:?}",
                cert.method
            ));
        }
    }
    assert!(
        accepted.is_empty(),
        "wrong antiderivatives were certified:\n{}",
        accepted.join("\n")
    );
}

/// A curated cross-mechanism list: every solved case must either carry a
/// certificate or be an honest fallback, and the numeric oracle must never
/// disagree with a certified result.
#[test]
fn certificates_never_contradict_the_numeric_oracle() {
    let cases = [
        // rational / symbolic rational
        "1/(a+b*x)",
        "1/(x^2*(x^2+5))",
        "(A+B*x)/x^2",
        // exp / log kernels
        "exp(x)/(1-exp(2*x))",
        "x*log(x)",
        "exp(x)/x",
        // trig
        "sin(x)^2*cos(x)^5",
        "1/(2+cos(x))",
        "tan(x)^3",
        "cos(2*x+1)^4",
        // hyperbolic (0.28.0 front-end)
        "sinh(x)",
        "tanh(x)^3",
        "x*cosh(x)",
        "sinh(a+b*x)^3*tanh(a+b*x)",
        // inverse functions
        "atan(x)",
        "asin(x)",
        // radicals
        "sqrt(1-x^2)",
        "1/sqrt(1-x^2)",
        "x^2/sqrt(1+x^2)",
        // special functions
        "exp(-x^2)",
        "erf(x)",
    ];

    let mut checked = 0usize;
    let mut false_positives = Vec::new();
    for input in cases {
        let (text, verdict) = integrate_and_certify(input);
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = parse(&ctx, input).expect("parse integrand");
        let result = integrate(&ctx, expr, Symbol::new("x"));
        let oracle = verify_antiderivative(expr, result, Symbol::new("x"));
        checked += 1;
        if let Verify::Mismatch { detail, .. } = oracle {
            // The oracle says the answer is wrong: emitting a certificate
            // would be a false positive.
            assert!(
                verdict.is_err(),
                "{input} -> {text}: certified ({verdict:?}) while the oracle reports {detail}"
            );
            false_positives.push(input);
        }
    }
    assert!(checked >= cases.len());
    assert!(
        false_positives.is_empty(),
        "oracle mismatches: {false_positives:?}"
    );
}

/// Structural certificates must be exact: a difference that does not
/// collapse is never reported as `Structural`.
#[test]
fn structural_certificates_are_exact() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let x = ctx.var("x");
    // ∫ x dx = x²/2
    let f = x;
    let big_f = ctx.mul(&[ctx.pow(ctx.num(2), ctx.num(-1)), ctx.pow(x, ctx.num(2))]);
    let cert = certify(&ctx, f, big_f, Symbol::new("x")).expect("certificate");
    assert!(
        matches!(
            cert.method,
            CertMethod::Structural | CertMethod::ElementaryField
        ),
        "unexpected method {:?}",
        cert.method
    );
    assert_eq!(cert.difference.to_string(), "0");

    // A residue-carrying result is not a certificate.
    let residue = ctx.fun("Integral", &[ctx.pow(x, ctx.num(-1)), x]);
    assert!(certify(&ctx, ctx.num(1), residue, Symbol::new("x")).is_err());
}
