//! Benchmark: polynomial factorization.
//!
//! Measures square-free and full factorization over ℤ and finite fields,
//! mirroring Symbolica's `factorization.rs` example.

use criterion::{Criterion, criterion_group, criterion_main};
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField, Integer, IntegerDomain};
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

// ── full factorization over ℤ ──────────────────────────────────────

fn poly_factor_z(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_factor_z");

    for degree in [12, 30, 60, 100] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        group.bench_with_input(format!("x_pow_{degree}_minus_1"), &degree, |bench, _| {
            bench.iter(|| {
                let f = black_box(&a).factor();
                black_box(f);
            });
        });
    }

    // (x^2+1)(x^2+x+1)(x+1) = x^5 + 2x^4 + 3x^3 + 3x^2 + 2x + 1
    let mixed = build_poly(&[1, 2, 3, 3, 2, 1]);
    group.bench_function("mixed_deg5", |b| {
        b.iter(|| {
            let f = black_box(&mixed).factor();
            black_box(f);
        });
    });

    group.finish();
}

// ── full factorization over F_p ────────────────────────────────────

fn poly_factor_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_factor_fp");

    for p in [5u64, 7, 17] {
        let field = FiniteField::new(BigInt::from(p));
        let fp = DenseUnivariatePolynomial::from_coeffs(
            field.clone(),
            (0..=100)
                .map(|i| {
                    if i == 0 {
                        field.element(-1i64)
                    } else if i == 100 {
                        field.element(1)
                    } else {
                        field.zero()
                    }
                })
                .collect(),
        );
        group.bench_with_input(format!("x100-1_over_F{p}"), &p, |bench, _| {
            bench.iter(|| {
                let f = black_box(&fp).factor();
                black_box(f);
            });
        });
    }

    group.finish();
}

// ── square-free benchmarks (existing) ──────────────────────────────

fn poly_square_free(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_factor_square_free");

    for degree in [12, 30, 60, 100] {
        let a = build_poly(&x_pow_n_minus_1(degree));
        group.bench_with_input(format!("x_pow_{degree}_minus_1"), &degree, |bench, _| {
            bench.iter(|| {
                let f = black_box(&a).square_free_factorization();
                black_box(f);
            });
        });
    }

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

criterion_group!(
    benches,
    poly_square_free,
    poly_is_square_free,
    poly_factor_z,
    poly_factor_fp
);
criterion_main!(benches);
