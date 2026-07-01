//! Benchmark: polynomial GCD over the integers.
//!
//! Measures [`DenseUnivariatePolynomial::gcd`] on structured inputs where
//! the GCD is known, mirroring the `polynomial_gcd.rs` example in Symbolica.

use criterion::{Criterion, criterion_group, criterion_main};
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;
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

criterion_group!(benches, poly_gcd);
criterion_main!(benches);
