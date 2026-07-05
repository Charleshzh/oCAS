//! Benchmark: Estrin vs Horner polynomial evaluation.
//!
//! Compares `fast_polynomial::poly` (Estrin's scheme) against classic
//! Horner's method for evaluating polynomials of various degrees.
//!
//! Run with: `cargo bench --bench eval_estrin --features fast-poly`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;

/// Horner's method (baseline).
fn horner(coeffs: &[f64], x: f64) -> f64 {
    let mut result = 0.0f64;
    for &c in coeffs.iter().rev() {
        result = result * x + c;
    }
    result
}

/// Estrin's scheme via `fast_polynomial`.
fn estrin(coeffs: &[f64], x: f64) -> f64 {
    ocas_eval::poly_eval::eval_estrin(coeffs, x)
}

fn build_coeffs(degree: usize) -> Vec<f64> {
    (0..=degree).map(|i| (i as f64) * 0.1 + 1.0).collect()
}

fn bench_estrin_vs_horner(c: &mut Criterion) {
    let mut group = c.benchmark_group("estrin_vs_horner");
    let x = 1.5_f64;

    for degree in [4, 8, 16, 32, 64, 128, 256] {
        let coeffs = build_coeffs(degree);

        group.bench_with_input(BenchmarkId::new("horner", degree), &degree, |bench, _| {
            bench.iter(|| {
                let result = horner(black_box(&coeffs), black_box(x));
                black_box(result);
            });
        });

        group.bench_with_input(BenchmarkId::new("estrin", degree), &degree, |bench, _| {
            bench.iter(|| {
                let result = estrin(black_box(&coeffs), black_box(x));
                black_box(result);
            });
        });
    }
    group.finish();
}

fn bench_eval_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("eval_batch");
    let coeffs = build_coeffs(32);
    let xs: Vec<f64> = (0..10_000).map(|i| i as f64 * 0.001).collect();

    group.bench_function("horner_batch_10k", |bench| {
        bench.iter(|| {
            let results: Vec<f64> = black_box(&xs)
                .iter()
                .map(|&x| horner(black_box(&coeffs), x))
                .collect();
            black_box(results);
        });
    });

    group.bench_function("estrin_batch_10k", |bench| {
        bench.iter(|| {
            let results =
                ocas_eval::poly_eval::eval_estrin_batch(black_box(&coeffs), black_box(&xs));
            black_box(results);
        });
    });

    group.finish();
}

fn bench_eval_single(c: &mut Criterion) {
    // Single evaluation at various degrees — the hot path for symbolic computation
    let mut group = c.benchmark_group("eval_single");

    for degree in [8, 32, 64, 128] {
        let coeffs = build_coeffs(degree);
        let x = 2.0_f64;

        group.bench_with_input(BenchmarkId::new("horner", degree), &degree, |bench, _| {
            bench.iter(|| {
                black_box(horner(black_box(&coeffs), black_box(x)));
            });
        });

        group.bench_with_input(BenchmarkId::new("estrin", degree), &degree, |bench, _| {
            bench.iter(|| {
                black_box(estrin(black_box(&coeffs), black_box(x)));
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_estrin_vs_horner,
    bench_eval_batch,
    bench_eval_single
);
criterion_main!(benches);
