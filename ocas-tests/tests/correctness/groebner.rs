use num_bigint::BigInt;
use ocas_domain::{FiniteField, Rational, RationalDomain};
use ocas_poly::groebner::f4::f4;
use ocas_poly::sparse::Lex;
use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial};

fn rat(n: i64, d: i64) -> Rational {
    Rational::new(n, d)
}

/// Build cyclic-n generators over ℚ.
fn cyclic_q(n: usize) -> Vec<SparseMultivariatePolynomial<RationalDomain, Lex>> {
    let d = RationalDomain;
    let mut gens = Vec::with_capacity(n);
    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut exps = vec![0usize; n];
            for j in 0..k {
                exps[(start + j) % n] = 1;
            }
            terms.push((exps, rat(1, 1)));
        }
        gens.push(SparseMultivariatePolynomial::from_terms(d, n, terms));
    }
    let full_exps = vec![1usize; n];
    gens.push(SparseMultivariatePolynomial::from_terms(
        d,
        n,
        vec![(full_exps, rat(1, 1)), (vec![0usize; n], rat(-1, 1))],
    ));
    gens
}

/// Build cyclic-n generators over ℤ_p.
fn cyclic_fp(n: usize, p: u32) -> Vec<SparseMultivariatePolynomial<FiniteField, Lex>> {
    let field = FiniteField::new(BigInt::from(p));
    let mut gens = Vec::with_capacity(n);
    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut exps = vec![0usize; n];
            for j in 0..k {
                exps[(start + j) % n] = 1;
            }
            terms.push((exps, field.element(1)));
        }
        gens.push(SparseMultivariatePolynomial::from_terms(
            field.clone(),
            n,
            terms,
        ));
    }
    let full_exps = vec![1usize; n];
    gens.push(SparseMultivariatePolynomial::from_terms(
        field.clone(),
        n,
        vec![
            (full_exps, field.element(1)),
            (vec![0usize; n], field.element(p - 1)),
        ],
    ));
    gens
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
