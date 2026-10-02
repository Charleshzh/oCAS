//! Property tests over random differential towers (0.29.0, Wave F).
//!
//! Two generators plus an adversarial family:
//!
//! 1. **Integrate-by-construction**: pick a random tower (1–3 `log`/`exp`
//!    generators over small `x`-polynomials) and a random element `F` over
//!    it, differentiate `f = D F` exactly, and integrate `f` through
//!    [`crate::integral::outcome::integrate_outcome`]. `Found` carries a
//!    certificate by API construction; `Unknown` is honest. Any panic or
//!    `ProvedNonElementary` (no producer exists yet) fails.
//! 2. **Random integrands**: random rational expressions of the tower
//!    generators go through the *legacy* [`crate::integrate`] path; every
//!    residue-free answer it emits must carry a valid certificate — the
//!    `certified_rate = 1.0` gate on this family (a decline is a hard
//!    failure: these expressions are small enough that the certificate
//!    budgets never trip legitimately).
//! 3. **Adversarial fixed cases**: dependent generators
//!    (`exp(u)·exp(−u)`, `log(x) + log(3x)`, `exp(x) + exp(2x)`), which
//!    exercise the merge machinery end to end through the public pipeline.

#[cfg(test)]
mod tests {
    use ocas_atom::normalize::normalize;
    use ocas_atom::{AtomArena, Symbol};
    use ocas_core::arena::Arena;
    use proptest::prelude::*;

    use crate::integral::certify::certify;
    use crate::integral::outcome::{Outcome, integrate_outcome};

    /// A random expression (as a string) over the given leaf atoms: a
    /// depth-≤2 tree of `+`/`·`/`^k` with `k ∈ {2, 3, −1, −2}`.
    fn expr_strategy(leaves: Vec<String>) -> impl Strategy<Value = String> {
        prop::sample::select(leaves).prop_recursive(2, 8, 2, |inner| {
            prop_oneof![
                (inner.clone(), inner.clone()).prop_map(|(a, b)| format!("({a}) + ({b})")),
                (inner.clone(), inner.clone()).prop_map(|(a, b)| format!("({a}) * ({b})")),
                (inner, prop::sample::select(vec![2i64, 3, -1, -2]))
                    .prop_map(|(a, k)| format!("({a})^({k})")),
            ]
        })
    }

    /// Tower shapes: 1–3 generators of log/exp over small `x`-polynomials.
    fn tower_args() -> impl Strategy<Value = Vec<(&'static str, String)>> {
        prop::collection::vec(
            (
                prop::sample::select(vec!["log", "exp"]),
                prop::sample::select(vec![
                    "x".to_string(),
                    "x^2".to_string(),
                    "x + 1".to_string(),
                    "x^2 + x".to_string(),
                    "2*x".to_string(),
                ]),
            ),
            1..=3,
        )
    }

    /// The leaf atoms for a tower: `x`, small constants, and the
    /// generator applications themselves.
    fn leaves_of(gens: &[(&'static str, String)]) -> Vec<String> {
        let mut leaves = vec!["x".to_string(), "2".to_string(), "-1".to_string()];
        leaves.extend(gens.iter().map(|(kind, arg)| format!("{kind}({arg})")));
        leaves
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        /// Integrate-by-construction over random towers.
        #[test]
        fn tower_constructed_elements_integrate_safely(
            (gens, f_src) in tower_args().prop_flat_map(|gens| {
                let leaves = leaves_of(&gens);
                expr_strategy(leaves).prop_map(move |f| (gens.clone(), f))
            }),
        ) {
            let _ = &gens; // the leaves already embed the generators
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x = Symbol::new("x");
            let big_f = ocas_parse::parse(&ctx, &f_src).expect("parse F");
            let f = normalize(&ctx, crate::diff(&ctx, big_f, x));
            // Size guard: differentiation can blow the expression up; the
            // certificate budgets are sized for corpus-scale inputs.
            if crate::integral::node_count(f) > 400 {
                return Ok(());
            }
            match integrate_outcome(&ctx, f, x) {
                Outcome::Found { value, .. } => {
                    prop_assert!(!value.to_string().contains("Integral("));
                }
                Outcome::Unknown { .. } => {}
                Outcome::ProvedNonElementary { .. } => {
                    return Err(TestCaseError::fail(
                        "no ProvedNonElementary producer exists in 0.29.0",
                    ));
                }
            }
        }

        /// Random tower integrands through the legacy path: residue-free
        /// answers must certify.
        #[test]
        fn tower_random_integrands_never_emit_uncertified(
            f_src in tower_args().prop_flat_map(|gens| {
                expr_strategy(leaves_of(&gens))
            }),
        ) {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let x = Symbol::new("x");
            let expr = ocas_parse::parse(&ctx, &f_src).expect("parse f");
            if crate::integral::node_count(expr) > 200 {
                return Ok(());
            }
            let normalized = normalize(&ctx, expr);
            let result = crate::integrate(&ctx, expr, x);
            if result.to_string().contains("Integral(") {
                return Ok(()); // honest fallback / partial
            }
            let cert = certify(&ctx, normalized, result, x);
            prop_assert!(
                cert.is_ok(),
                "legacy integrate emitted an uncertifiable answer for {f_src}: \
                 {result} ({:?})",
                cert.err(),
            );
        }
    }

    /// Dependent-generator adversarial family: the merge machinery must
    /// keep these certifiable end to end.
    #[test]
    fn adversarial_dependent_generators_certify() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = Symbol::new("x");
        for src in [
            // exp(u)·exp(−u) = 1 lives in the integrand itself.
            "exp(x)*exp(-x) + x",
            // log(x) + log(3x): merged with a constant generator.
            "log(x) + log(3*x)",
            // exp(2x) = exp(x)² (0.29.0 merge extension).
            "exp(x) + exp(2*x)",
        ] {
            let expr = ocas_parse::parse(&ctx, src).expect("parse");
            let result = crate::integrate(&ctx, expr, x);
            assert!(
                !result.to_string().contains("Integral("),
                "dependent-generator family must integrate: {src} -> {result}"
            );
            let cert = certify(&ctx, normalize(&ctx, expr), result, x);
            assert!(cert.is_ok(), "certificate failed for {src}: {result}");
        }
    }
}
