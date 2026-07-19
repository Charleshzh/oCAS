//! Benchmark for the Cranelift JIT backend.
//!
//! Measures JIT compilation latency and compares evaluation throughput
//! against the scalar interpreter. Enabled with the `jit` feature.

use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;
use ocas_eval::instruction::Instr;
use ocas_eval::jit::JitEngine;
use std::hint::black_box;

fn bench_jit_compile(c: &mut Criterion) {
    // Build a simple instruction sequence: x + y
    let instructions = vec![Instr::Add {
        dst: 2,
        srcs: vec![0, 1],
    }];

    c.bench_function("jit_compile_add", |b| {
        b.iter(|| {
            let func = JitEngine::compile(black_box(&instructions), 2, &[], &[2]);
            let _ = black_box(func);
        });
    });
}

fn bench_jit_exec_single_output(c: &mut Criterion) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "x^4 + 3*x^3 + 2*x^2 + x + 5").unwrap();
    let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
    let jit = eval.compile_jit().unwrap();

    let inputs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.01).collect();

    let mut group = c.benchmark_group("jit_exec_poly_1000");

    group.bench_function("interpreter", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                let result = eval.evaluate(black_box(&[x])).unwrap();
                sum += result[0];
            }
            black_box(sum);
        });
    });

    group.bench_function("jit", |b| {
        let mut out = [0.0f64; 1];
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                jit.call_into(black_box(&[x]), &mut out);
                sum += out[0];
            }
            black_box(sum);
        });
    });

    group.finish();
}

fn bench_jit_exec_multi_output(c: &mut Criterion) {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let sin_x = parse(&ctx, "sin(x)").unwrap();
    let cos_x = parse(&ctx, "cos(x)").unwrap();
    let tanh_like = parse(&ctx, "sin(x) / cos(x)").unwrap();
    let eval: ExpressionEvaluator<f64> =
        ExpressionEvaluator::compile_multi(&[sin_x, cos_x, tanh_like]).unwrap();
    let jit = eval.compile_jit().unwrap();

    let inputs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.01).collect();

    let mut group = c.benchmark_group("jit_exec_multi3_1000");

    group.bench_function("interpreter", |b| {
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                let result = eval.evaluate(black_box(&[x])).unwrap();
                sum += result[0] + result[1] + result[2];
            }
            black_box(sum);
        });
    });

    group.bench_function("jit", |b| {
        let mut out = [0.0f64; 3];
        b.iter(|| {
            let mut sum = 0.0_f64;
            for &x in black_box(&inputs) {
                jit.call_into(black_box(&[x]), &mut out);
                sum += out[0] + out[1] + out[2];
            }
            black_box(sum);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_jit_compile,
    bench_jit_exec_single_output,
    bench_jit_exec_multi_output,
);
criterion_main!(benches);
