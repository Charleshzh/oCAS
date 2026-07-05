//! Benchmark: pulp SIMD vs scalar expression evaluation.
//!
//! Compares the `VectorEvaluator` (pulp SIMD) against the scalar
//! `ExpressionEvaluator` on batch inputs of various sizes.
//!
//! Run with: `cargo bench --bench eval_pulp --features simd`

use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_eval::VectorEvaluator;
use std::hint::black_box;

fn build_evaluators(expr_str: &str) -> (ExpressionEvaluator<f64>, VectorEvaluator) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, expr_str).unwrap();
    let scalar: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
    let vector = scalar.compile_vector_evaluator().unwrap();
    (scalar, vector)
}

fn bench_poly_batch(c: &mut Criterion) {
    let (scalar, vector) = build_evaluators("x^4 + 3*x^3 + 2*x^2 + x + 5");

    for batch_size in [64, 256, 1000, 4000, 10_000] {
        let inputs: Vec<f64> = (0..batch_size).map(|i| i as f64 * 0.01).collect();
        let params = vec![inputs];

        let mut group = c.benchmark_group(format!("poly_batch_{batch_size}"));

        group.bench_function("scalar", |bench| {
            bench.iter(|| {
                let mut sum = 0.0_f64;
                for &x in black_box(&params[0]) {
                    let result = scalar.evaluate(black_box(&[x])).unwrap();
                    sum += result[0];
                }
                black_box(sum);
            });
        });

        group.bench_function("simd_pulp", |bench| {
            bench.iter(|| {
                let result = vector.evaluate(black_box(&params)).unwrap();
                black_box(result);
            });
        });

        group.finish();
    }
}

fn bench_trig_batch(c: &mut Criterion) {
    let (scalar, vector) = build_evaluators("sin(x)^2 + cos(x)^2 + exp(x)");

    for batch_size in [64, 256, 1000, 4000] {
        let inputs: Vec<f64> = (0..batch_size).map(|i| i as f64 * 0.01).collect();
        let params = vec![inputs];

        let mut group = c.benchmark_group(format!("trig_batch_{batch_size}"));

        group.bench_function("scalar", |bench| {
            bench.iter(|| {
                let mut sum = 0.0_f64;
                for &x in black_box(&params[0]) {
                    let result = scalar.evaluate(black_box(&[x])).unwrap();
                    sum += result[0];
                }
                black_box(sum);
            });
        });

        group.bench_function("simd_pulp", |bench| {
            bench.iter(|| {
                let result = vector.evaluate(black_box(&params)).unwrap();
                black_box(result);
            });
        });

        group.finish();
    }
}

fn bench_high_degree_batch(c: &mut Criterion) {
    // Higher-degree polynomial: x^8 + x^6 + x^4 + x^2 + 1
    let (scalar, vector) = build_evaluators("x^8 + x^6 + x^4 + x^2 + 1");

    let batch_size = 4000;
    let inputs: Vec<f64> = (0..batch_size).map(|i| i as f64 * 0.001).collect();
    let params = vec![inputs];

    let mut group = c.benchmark_group("high_degree_batch_4k");

    group.bench_function("scalar", |bench| {
        bench.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&params[0]) {
                let result = scalar.evaluate(black_box(&[x])).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });

    group.bench_function("simd_pulp", |bench| {
        bench.iter(|| {
            let result = vector.evaluate(black_box(&params)).unwrap();
            black_box(result);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_poly_batch,
    bench_trig_batch,
    bench_high_degree_batch
);
criterion_main!(benches);
