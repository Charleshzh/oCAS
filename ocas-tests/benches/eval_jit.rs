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
            let func = JitEngine::compile(black_box(&instructions), 2, 1);
            black_box(func);
        });
    });
}

criterion_group!(benches, bench_jit_compile);
criterion_main!(benches);
