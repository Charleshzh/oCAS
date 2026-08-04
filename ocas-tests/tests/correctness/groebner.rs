use num_bigint::BigInt;
use ocas_domain::{FiniteField, Rational, RationalDomain};
use ocas_poly::groebner::f4::f4;
use ocas_poly::groebner::f5::f5;
use ocas_poly::sparse::Lex;
use ocas_poly::{Algorithm, GroebnerBasis, SparseMultivariatePolynomial, groebner_basis};
use ocas_tests::systems::*;

fn rat(n: i64, d: i64) -> Rational {
    Rational::new(n, d)
}

// =========================================================================
//  Buchberger tests (existing, with #[ignore] removed)
// =========================================================================

#[test]
fn groebner_simple_linear_system() {
    let d = RationalDomain;
    let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        2,
        vec![
            (vec![0, 0], rat(-1, 1)),
            (vec![1, 0], rat(1, 1)),
            (vec![0, 1], rat(1, 1)),
        ],
    );
    let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        2,
        vec![
            (vec![0, 0], rat(-3, 1)),
            (vec![1, 0], rat(1, 1)),
            (vec![0, 1], rat(-1, 1)),
        ],
    );
    let gb = GroebnerBasis::buchberger(&[f1, f2]);
    assert!(!gb.basis.is_empty());
}

#[test]
fn groebner_buchberger_cyclic_3() {
    let ideal = cyclic_q(3);
    let gb = GroebnerBasis::buchberger(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
#[ignore = "Buchberger cyclic-4 is very slow; use F4 instead"]
fn groebner_buchberger_cyclic_4() {
    let ideal = cyclic_q(4);
    let gb = GroebnerBasis::buchberger(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

// =========================================================================
//  F4 tests over ℚ
// =========================================================================

#[test]
fn groebner_f4_cyclic_3_q() {
    let ideal = cyclic_q(3);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f4_cyclic_4_q() {
    let ideal = cyclic_q(4);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

// =========================================================================
//  F4 tests over ℤ_p
// =========================================================================

#[test]
fn groebner_f4_cyclic_3_fp13() {
    let ideal = cyclic_fp(3, 13);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f4_cyclic_4_fp13() {
    let ideal = cyclic_fp(4, 13);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f4_cyclic_3_fp101() {
    let ideal = cyclic_fp(3, 101);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f4_cyclic_4_fp101() {
    let ideal = cyclic_fp(4, 101);
    let gb = f4(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

// =========================================================================
//  F4 correctness: Buchberger vs F4 agreement
// =========================================================================

#[test]
fn groebner_f4_vs_buchberger_cyclic_3() {
    let ideal = cyclic_q(3);
    let gb_buch = GroebnerBasis::buchberger(&ideal);
    let gb_f4 = f4(&ideal);
    // Both should produce valid Gröbner bases.
    assert!(gb_buch.is_groebner_basis());
    assert!(gb_f4.is_groebner_basis());
    // Both should generate the same ideal (verified by is_groebner_basis).
}

// =========================================================================
//  F5 tests
// =========================================================================

#[test]
fn groebner_f5_cyclic_3_q() {
    let ideal = cyclic_q(3);
    let gb = f5(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f5_cyclic_3_fp13() {
    let ideal = cyclic_fp(3, 13);
    let gb = f5(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f5_vs_f4_cyclic_3() {
    let ideal = cyclic_q(3);
    let gb_f5 = f5(&ideal);
    let gb_f4 = f4(&ideal);
    assert!(gb_f5.is_groebner_basis());
    assert!(gb_f4.is_groebner_basis());
}

#[test]
fn groebner_f5_cyclic_4_q() {
    let ideal = cyclic_q(4);
    let gb = f5(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f5_cyclic_4_fp13() {
    let ideal = cyclic_fp(4, 13);
    let gb = f5(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
fn groebner_f5_cyclic_5_fp13() {
    let ideal = cyclic_fp(5, 13);
    let gb = f5(&ideal);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

#[test]
#[ignore = "timing test; run with --release --ignored"]
fn groebner_f5_cyclic_6_fp13() {
    // Roadmap acceptance target: cyclic-6 ℤ₁₃ in < 0.5 s median (criterion
    // f5_cyclic_fp13/cyclic_6 records the median; this test keeps a 2 s
    // wall-clock bound for CI noise). Baseline (F4, 0.15.2): 3670s.
    let ideal = cyclic_fp(6, 13);
    let start = std::time::Instant::now();
    let gb = f5(&ideal);
    let elapsed = start.elapsed();
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
    eprintln!("cyclic-6 ℤ₁₃ F5: {elapsed:.2?}");
    assert!(
        elapsed.as_secs() < 2,
        "cyclic-6 ℤ₁₃ took {elapsed:.2?}, expected < 2s"
    );
}

#[test]
#[ignore = "reference benchmark; run manually with --ignored --release"]
fn groebner_f5_cyclic_7_fp13() {
    // Reference benchmark (no hard time bound): cyclic-7 ℤ₁₃ must be
    // solvable and produce a valid Gröbner basis.
    let ideal = cyclic_fp(7, 13);
    let start = std::time::Instant::now();
    let gb = f5(&ideal);
    let elapsed = start.elapsed();
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
    eprintln!("cyclic-7 ℤ₁₃ F5: {elapsed:.2?}");
}

// =========================================================================
//  Multi-modular tests
// =========================================================================

/// The MultiModular algorithm variant must route through the multi-modular
/// pipeline and produce the same reduced basis as F5 (ℚ coefficients).
#[test]
fn algorithm_multi_modular_variant_routes() {
    let ideal = cyclic_q(3);
    let gb_mm = groebner_basis(&ideal, Algorithm::MultiModular);
    let gb_f5 = groebner_basis(&ideal, Algorithm::F5);
    assert!(gb_mm.is_groebner_basis());
    assert_eq!(gb_mm, gb_f5, "MultiModular must agree with F5 over ℚ");
}

/// The Auto variant must route ℚ ideals through the multi-modular pipeline
/// with identical results, and stay correct for ℤ_p (F4 path).
#[test]
fn auto_routes_multi_modular_for_q() {
    let ideal = cyclic_q(3);
    let gb_auto = groebner_basis(&ideal, Algorithm::Auto);
    let gb_f4 = groebner_basis(&ideal, Algorithm::F4);
    assert!(gb_auto.is_groebner_basis());
    assert_eq!(gb_auto, gb_f4, "Auto must agree with F4 over ℚ");

    // ℤ_p is not a ℚ ideal: Auto must fall back to F4 and stay correct.
    let ideal_fp = cyclic_fp(3, 13);
    let gb_auto = groebner_basis(&ideal_fp, Algorithm::Auto);
    let gb_f4 = groebner_basis(&ideal_fp, Algorithm::F4);
    assert!(gb_auto.is_groebner_basis());
    assert_eq!(gb_auto, gb_f4, "Auto must agree with F4 over ℤ_p");
}

/// Multi-modular must agree with F5 on 100 deterministic random ℚ ideals.
///
/// Random ideals are 3-variable, 3-generator, total degree ≤ 3, integer
/// coefficients in [−10, 10]. The reduced bases must be term-for-term
/// identical.
#[test]
fn multi_modular_matches_f5_random() {
    use rand_xoshiro::Xoshiro256PlusPlus;
    use rand_xoshiro::rand_core::{RngCore, SeedableRng};

    fn rand_q_ideal(rng: &mut Xoshiro256PlusPlus) -> Vec<SparseMultivariatePolynomial<RationalDomain, Lex>> {
        let n_vars = 3;
        let mut gens = Vec::with_capacity(3);
        for _ in 0..3 {
            let n_terms = 1 + (rng.next_u64() % 5) as usize;
            let mut terms = Vec::with_capacity(n_terms);
            for _ in 0..n_terms {
                let mut exps = vec![0usize; n_vars];
                // Total degree ≤ 3, distributed among the variables.
                let deg = 1 + (rng.next_u64() % 3) as usize;
                for _ in 0..deg {
                    exps[(rng.next_u64() % n_vars as u64) as usize] += 1;
                }
                let c = (rng.next_u64() % 21) as i64 - 10;
                if c != 0 {
                    terms.push((exps, rat(c, 1)));
                }
            }
            if terms.is_empty() {
                terms.push((vec![0usize; n_vars], rat(1, 1)));
            }
            gens.push(SparseMultivariatePolynomial::from_terms(
                RationalDomain,
                n_vars,
                terms,
            ));
        }
        gens
    }

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xC0FFEE);
    for trial in 0..100 {
        let ideal = rand_q_ideal(&mut rng);
        let gb_mm = groebner_basis(&ideal, Algorithm::MultiModular);
        let gb_f5 = groebner_basis(&ideal, Algorithm::F5);
        assert!(
            gb_mm == gb_f5,
            "multi-modular and F5 disagree on random ideal #{trial}: {ideal:?}"
        );
    }
}

// =========================================================================
//  Packed ℤ_p fast path (0.26.0)
// =========================================================================

/// The packed path requires n_vars ≤ 8; a 9-variable ideal must fall back
/// to the i64 fast path and still produce a valid Gröbner basis.
#[test]
fn f5_fp_packed_fallback_nvars_gt_8() {
    let field = FiniteField::new(BigInt::from(13));
    let mut terms = Vec::new();
    for i in 0..9 {
        let mut e = vec![0usize; 9];
        e[i] = 1;
        terms.push((e, field.element(1)));
    }
    let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(field.clone(), 9, terms);
    let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        field.clone(),
        9,
        vec![
            (vec![1, 1, 0, 0, 0, 0, 0, 0, 0], field.element(1)),
            (vec![0usize; 9], field.element(12)),
        ],
    );
    let gb = f5(&[f1, f2]);
    assert!(!gb.basis.is_empty());
    assert!(gb.is_groebner_basis());
}

/// The packed path requires input exponents < 2^15; x^70000 must fall back
/// to the i64 fast path and produce the expected reduced basis.
#[test]
fn f5_fp_packed_fallback_exp_overflow() {
    let field = FiniteField::new(BigInt::from(13));
    let h1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        field.clone(),
        2,
        vec![
            (vec![70000, 0], field.element(1)),
            (vec![0usize; 2], field.element(12)),
        ],
    );
    let h2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        field.clone(),
        2,
        vec![
            (vec![0, 1], field.element(1)),
            (vec![0usize; 2], field.element(12)),
        ],
    );
    let gb = f5(&[h1, h2]);
    assert!(gb.is_groebner_basis());
    // Reduced basis of (x^70000 - 1, y - 1) is exactly {y - 1, x^70000 - 1}.
    assert_eq!(gb.basis.len(), 2);
    let mut exps: Vec<Vec<usize>> = gb.basis.iter().map(|p| p.leading_monomial().unwrap().to_vec()).collect();
    exps.sort();
    assert_eq!(exps, vec![vec![0, 1], vec![70000, 0]]);
}

/// The packed pipeline (parallel row construction + two-phase echelon)
/// is deterministic: three runs of cyclic-4 ℤ₁₃ produce term-for-term
/// identical reduced bases.
#[test]
fn f5_fp_packed_deterministic() {
    let ideal = cyclic_fp(4, 13);
    let gb1 = f5(&ideal);
    let gb2 = f5(&ideal);
    let gb3 = f5(&ideal);
    assert_eq!(gb1, gb2);
    assert_eq!(gb1, gb3);
    assert!(gb1.is_groebner_basis());
}
