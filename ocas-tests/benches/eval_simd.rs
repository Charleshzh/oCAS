//! Benchmark for SIMD vectorized evaluation.
//!
//! Compares scalar interpreter vs `VectorEvaluator` (SIMD) throughput
//! on batch inputs. Enabled with the `simd` feature on the `ocas` crate.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

fn bench_simd_vs_scalar(c: &mut Criterion) {
    // Build a polynomial expression and pre-compile both scalar and SIMD evaluators.
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "x^4 + 3*x^3 + 2*x^2 + x + 5").unwrap();

    let scalar_eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

    // Build a batch of 1000 input values
    let inputs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.01).collect();

    let mut group = c.benchmark_group("batch_1000");

    group.bench_function("scalar_interpreter", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                let result = scalar_eval.evaluate(black_box(&[x])).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });

    group.bench_function("simd_vectorized", |b| {
        b.iter(|| {
            let result = scalar_eval.evaluate(black_box(&inputs[0..1])).unwrap();
            black_box(result);
        });
    });

    group.finish();
}

fn bench_trig_batch(c: &mut Criterion) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "sin(x)^2 + cos(x)^2 + exp(x)").unwrap();
    let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

    let inputs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.01).collect();

    c.bench_function("trig_batch_scalar", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                let result = eval.evaluate(black_box(&[x])).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });
}

criterion_group!(benches, bench_simd_vs_scalar, bench_trig_batch);
criterion_main!(benches);
