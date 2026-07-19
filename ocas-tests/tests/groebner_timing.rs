//! Single-iteration Gröbner timing for large cyclic-n ideals.
//!
//! Criterion requires ≥10 samples, which is impractical when one
//! iteration takes minutes. These `#[ignore]`d tests time a single F4
//! run each and are run manually in release mode:
//!
//! ```text
//! cargo test -p ocas-tests --release --test groebner_timing -- --ignored --nocapture
//! ```

use std::time::Instant;

use num_bigint::BigInt;
use ocas_domain::FiniteField;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::groebner::f4::f4;
use ocas_poly::sparse::Lex;

/// The cyclic-n ideal generators over ℤ_p (mirror of benches/groebner.rs).
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
            (vec![0usize; n], field.element(p as i64 - 1)), // -1 mod p
        ],
    ));

    gens
}

fn time_f4(n: usize, p: u32) {
    let ideal = cyclic_fp(n, p);
    let start = Instant::now();
    let gb = f4(&ideal);
    let elapsed = start.elapsed();
    println!(
        "f4 cyclic_{n} over Z_{p}: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_5_fp13() {
    time_f4(5, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_6_fp13() {
    time_f4(6, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_7_fp13() {
    time_f4(7, 13);
}
