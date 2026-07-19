//! Benchmark for streaming evaluation.
//!
//! Compares per-row `evaluate()` (fresh stack each call) against
//! `StreamingEvaluator::for_each` (reused buffers) over a large row
//! stream.

use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_eval::StreamingEvaluator;
use std::hint::black_box;

fn bench_streaming(c: &mut Criterion) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "x^4 + 3*x^3 + 2*x^2 + x + 5").unwrap();
    let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

    let rows: Vec<[f64; 1]> = (0..100_000).map(|i| [i as f64 * 0.001]).collect();

    let mut group = c.benchmark_group("streaming_100k");

    group.bench_function("per_row_evaluate", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for row in black_box(&rows) {
                let result = eval.evaluate(row).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });

    group.bench_function("streaming_for_each", |b| {
        let mut stream = StreamingEvaluator::new(&eval);
        b.iter(|| {
            let mut sum = 0.0_f64;
            stream
                .for_each(black_box(&rows), |results| sum += results[0])
                .unwrap();
            black_box(sum);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_streaming);
criterion_main!(benches);
