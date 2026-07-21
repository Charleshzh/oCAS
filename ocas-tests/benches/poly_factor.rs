//! Benchmark: polynomial factorization.
//!
//! Measures square-free and full factorization over ℤ and finite fields,
//! mirroring Symbolica's `factorization.rs` example.

use criterion::{Criterion, criterion_group, criterion_main};
use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField, Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};
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

// ── multivariate factorization over ℤ (Wang EEZ) ───────────────────

type ZmPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

fn zm_poly(n_vars: usize, terms: &[(Vec<usize>, i64)]) -> ZmPoly {
    SparseMultivariatePolynomial::from_terms(
        IntegerDomain,
        n_vars,
        terms
            .iter()
            .map(|(e, c)| (e.clone(), Integer::from(*c)))
            .collect(),
    )
}

fn poly_factor_multivariate_z(c: &mut Criterion) {
    let mut group = c.benchmark_group("poly_factor_multivariate_z");

    // (x + y + z)(x - y + 2z)(x + y + 1): 3-var, 3 linear factors.
    let f1 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
    );
    let f2 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 1], 2)],
    );
    let f3 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 0], 1)],
    );
    let tri = f1.mul(&f2).mul(&f3);
    group.bench_function("trivariate_3_linear", |b| {
        b.iter(|| {
            let f = black_box(&tri).factor();
            black_box(f);
        });
    });

    // (x^2 + y + z)(x + y - z): 3-var, quadratic + linear.
    let g1 = zm_poly(
        3,
        &[(vec![2, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
    );
    let g2 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], -1)],
    );
    let quad = g1.mul(&g2);
    group.bench_function("trivariate_quad_linear", |b| {
        b.iter(|| {
            let f = black_box(&quad).factor();
            black_box(f);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    poly_square_free,
    poly_is_square_free,
    poly_factor_z,
    poly_factor_fp,
    poly_factor_multivariate_z
);
criterion_main!(benches);
