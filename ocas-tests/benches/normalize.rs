use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_atom::normalize::normalize;
use ocas_core::arena::Arena;
use std::hint::black_box;

fn normalize_small(c: &mut Criterion) {
    c.bench_function("normalize_small", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x + 0 + 1*x + 0")).unwrap();
            let result = normalize(&ctx, expr);
            black_box(result);
        });
    });
}

fn normalize_medium(c: &mut Criterion) {
    let input = "(a + b + c + d + e + f + g + h + i + j)^3";
    c.bench_function("normalize_medium", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box(input)).unwrap();
            let result = normalize(&ctx, expr);
            black_box(result);
        });
    });
}

fn normalize_large(c: &mut Criterion) {
    let terms: Vec<String> = (0..100).map(|i| format!("x{} + {}", i, i * 2)).collect();
    let input = format!("({})^2", terms.join(" + "));
    c.bench_function("normalize_large", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box(&input)).unwrap();
            let result = normalize(&ctx, expr);
            black_box(result);
        });
    });
}

criterion_group!(benches, normalize_small, normalize_medium, normalize_large);
criterion_main!(benches);
