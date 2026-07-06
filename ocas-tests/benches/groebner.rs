//! Benchmark: Gröbner basis computation.
//!
//! Compares Buchberger's algorithm vs F4 on cyclic-n ideals over ℚ and ℤ_p.
//!
//! Reference: Symbolica's `groebner_basis.rs` example, Faugère F4 (1999).

use criterion::{Criterion, criterion_group, criterion_main};
use num_bigint::BigInt;
use ocas_domain::{FiniteField, Rational, RationalDomain};
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::buchberger;
use ocas_poly::groebner::f4::f4;
use ocas_poly::sparse::Lex;
use std::hint::black_box;

// =========================================================================
//  Helpers
// =========================================================================

/// Build a single rational term.
fn term(exps: Vec<usize>, num: i64, den: i64) -> (Vec<usize>, Rational) {
    (exps, Rational::new(num, den))
}

/// The cyclic-n ideal generators over ℚ.
///
///   f_k = Σ_{start=0}^{n-1} x_{start} * x_{start+1} * ... * x_{start+k-1}   (indices mod n)
///   f_n = x_0 * x_1 * ... * x_{n-1} - 1
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
            terms.push(term(exps, 1, 1));
        }
        gens.push(SparseMultivariatePolynomial::from_terms(d, n, terms));
    }

    let full_exps = vec![1usize; n];
    gens.push(SparseMultivariatePolynomial::from_terms(
        d,
        n,
        vec![term(full_exps, 1, 1), term(vec![0usize; n], -1, 1)],
    ));

    gens
}

/// The cyclic-n ideal generators over ℤ_p.
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
            (vec![0usize; n], field.element(p - 1)), // -1 mod p
        ],
    ));

    gens
}

// =========================================================================
//  Benchmarks
// =========================================================================

/// Buchberger on cyclic-n over ℚ (only cyclic-3; cyclic-4 is too slow).
fn bench_buchberger_cyclic_q(c: &mut Criterion) {
    let mut group = c.benchmark_group("buchberger_cyclic_q");
    for n in [3] {
        let ideal = cyclic_q(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = buchberger(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℚ (generic path).
fn bench_f4_cyclic_q(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_q");
    for n in [3, 4] {
        let ideal = cyclic_q(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℤ_13 (fast i64 path).
fn bench_f4_cyclic_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_fp13");
    group.sample_size(10);
    for n in [3, 4, 5] {
        let ideal = cyclic_fp(n, 13);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℤ_101 (larger prime).
fn bench_f4_cyclic_fp101(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_fp101");
    group.sample_size(10);
    for n in [3, 4, 5] {
        let ideal = cyclic_fp(n, 101);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

// =========================================================================
//  Criterion harness
// =========================================================================

criterion_group!(
    benches,
    bench_buchberger_cyclic_q,
    bench_f4_cyclic_q,
    bench_f4_cyclic_fp,
    bench_f4_cyclic_fp101,
);
criterion_main!(benches);
