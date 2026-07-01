//! Benchmark: square-free factorization over the integers.
//!
//! Measures [`DenseUnivariatePolynomial::square_free_factorization`] on
//! cyclotomic-style inputs `x^n - 1`, which factor into many square-free
//! components (the cyclotomic polynomials). This mirrors Symbolica's
//! `factorization.rs` example.

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

fn poly_square_free(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_factor_square_free");

    // x^n - 1 factors into the product of cyclotomic polynomials
    // Φ_d(x) for d | n; all distinct, so square-free.
    for degree in [12, 30, 60, 100] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        group.bench_with_input(format!("x_pow_{degree}_minus_1"), &degree, |bench, _| {
            bench.iter(|| {
                let f = black_box(&a).square_free_factorization();
                black_box(f);
            });
        });
    }

    // (x - 1)^2 * (x + 1) = x^3 - x^2 - x + 1 ... use repeated factors.
    // p = (x^2 - 1)^2 = x^4 - 2x^2 + 1 -> factor = (x^2-1)^2
    let repeated = build_poly(&[1, 0, -2, 0, 1]); // (x^2-1)^2
    group.bench_function("x2_minus_1_squared", |b| {
        b.iter(|| {
            let f = black_box(&repeated).square_free_factorization();
            black_box(f);
        });
    });

    group.finish();
}

fn poly_is_square_free(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_is_square_free");
    for degree in [12, 60, 100] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        group.bench_with_input(format!("x_pow_{degree}_minus_1"), &degree, |bench, _| {
            bench.iter(|| {
                let r = black_box(&a).is_square_free();
                black_box(r);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, poly_square_free, poly_is_square_free);
criterion_main!(benches);
