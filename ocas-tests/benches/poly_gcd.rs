//! Benchmark: polynomial GCD over the integers.
//!
//! Measures [`DenseUnivariatePolynomial::gcd`] on structured inputs where
//! the GCD is known, mirroring the `polynomial_gcd.rs` example in Symbolica.

use criterion::{Criterion, criterion_group, criterion_main};
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::multivariate_gcd::{bivariate_gcd, gcd_modular};
use ocas_poly::sparse::Lex;
use std::hint::black_box;

fn build_poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<IntegerDomain> {
    let coeffs: Vec<Integer> = coeffs.iter().map(|&i| Integer::from(i)).collect();
    DenseUnivariatePolynomial::from_coeffs(IntegerDomain, coeffs)
}

/// `x^n - 1` as a coefficient vector (constant term first).
fn x_pow_n_minus_1(n: usize) -> Vec<i64> {
    let mut c = vec![0i64; n + 1];
    c[0] = -1;
    c[n] = 1;
    c
}

fn poly_gcd(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_gcd");
    // gcd(x^n - 1, x - 1) = x - 1. Both share the factor (x - 1); this
    // exercises the Euclidean algorithm on a high-degree dividend.
    for degree in [10, 50, 100, 500] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        let b = build_poly(&[-1, 1]);
        group.bench_with_input(
            format!("gcd_x_pow_{degree}_minus_1_with_x_minus_1"),
            &degree,
            |bench, _| {
                bench.iter(|| {
                    let g = black_box(&a).gcd(black_box(&b));
                    black_box(g);
                });
            },
        );
    }

    // gcd of two random-ish coprime-ish polys of equal degree.
    for degree in [20, 100] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        // (x^n - 1) shifted: replace constant with 2 -> coprime to (x^n-1)
        let mut b_coeffs = x_pow_n_minus_1(degree);
        b_coeffs[0] = 2;
        let b = build_poly(&b_coeffs);
        group.bench_with_input(
            format!("gcd_coprime_degree_{degree}"),
            &degree,
            |bench, _| {
                bench.iter(|| {
                    let g = black_box(&a).gcd(black_box(&b));
                    black_box(g);
                });
            },
        );
    }
    group.finish();
}

/// Benchmark bivariate GCD over ℤ: heuristic vs modular approach.
fn bivariate_gcd_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("bivariate_gcd");

    // (x+y)(x+1) and (x+y)(x+2) — shared factor x+y
    let domain = IntegerDomain;
    type ZMPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

    let a: ZMPoly = SparseMultivariatePolynomial::from_terms(
        domain,
        2,
        vec![
            (vec![2, 0], Integer::from(1)),
            (vec![1, 1], Integer::from(1)),
            (vec![1, 0], Integer::from(1)),
            (vec![0, 1], Integer::from(1)),
        ],
    );
    let b: ZMPoly = SparseMultivariatePolynomial::from_terms(
        domain,
        2,
        vec![
            (vec![2, 0], Integer::from(1)),
            (vec![1, 1], Integer::from(1)),
            (vec![1, 0], Integer::from(2)),
            (vec![0, 1], Integer::from(2)),
        ],
    );

    group.bench_function("heuristic_bivariate", |bench| {
        bench.iter(|| {
            let g = bivariate_gcd(black_box(&a), black_box(&b));
            black_box(g);
        });
    });

    group.bench_function("modular_bivariate", |bench| {
        bench.iter(|| {
            let g = gcd_modular(black_box(&a), black_box(&b));
            black_box(g);
        });
    });

    group.finish();
}

criterion_group!(benches, poly_gcd, bivariate_gcd_bench);
criterion_main!(benches);
