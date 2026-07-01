//! Benchmark for the expression interpreter.
//!
//! Measures compilation time and evaluation throughput for the
//! `ExpressionEvaluator<f64>` stack VM.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

fn bench_compile(c: &mut Criterion) {
    c.bench_function("compile_polynomial", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x^4 + 3*x^3 + 2*x^2 + x + 5")).unwrap();
            let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
            black_box(eval);
        });
    });

    c.bench_function("compile_trig", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("sin(x)^2 + cos(x)^2 + exp(x)")).unwrap();
            let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
            black_box(eval);
        });
    });
}

fn bench_evaluate(c: &mut Criterion) {
    // Pre-compile expressions, then benchmark just the evaluation

    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let poly = parse(&ctx, "x^4 + 3*x^3 + 2*x^2 + x + 5").unwrap();
    let poly_eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(poly).unwrap();

    let arena2 = Arena::new();
    let ctx2 = AtomArena::new(&arena2);
    let trig = parse(&ctx2, "sin(x)^2 + cos(x)^2 + exp(x)").unwrap();
    let trig_eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(trig).unwrap();

    c.bench_function("eval_polynomial", |b| {
        b.iter(|| {
            let result = poly_eval.evaluate(black_box(&[2.0])).unwrap();
            black_box(result);
        });
    });

    c.bench_function("eval_trig_expression", |b| {
        b.iter(|| {
            let result = trig_eval.evaluate(black_box(&[1.5])).unwrap();
            black_box(result);
        });
    });

    // Batch evaluation: many evaluations of the same expression
    c.bench_function("eval_polynomial_batch_1000", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for i in 0..1000 {
                let x = i as f64 * 0.01;
                let result = poly_eval.evaluate(black_box(&[x])).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });
}

criterion_group!(benches, bench_compile, bench_evaluate);
criterion_main!(benches);
