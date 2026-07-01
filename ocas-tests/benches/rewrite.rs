use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_core::arena::Arena;
use std::hint::black_box;

fn bench_simplify(c: &mut Criterion) {
    c.bench_function("simplify_small", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let rules = ocas_rewrite::rules::default_rules(&ctx, &());
            let expr = parse(&ctx, black_box("x + x + 0 + 1*x")).unwrap();
            let result = simplify(&ctx, expr, &rules, 20);
            black_box(result);
        });
    });

    c.bench_function("simplify_nested", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let rules = ocas_rewrite::rules::default_rules(&ctx, &());
            let expr = parse(&ctx, black_box("(x + 0) * (y * 1) + (x * 0) + (z^1)")).unwrap();
            let result = simplify(&ctx, expr, &rules, 20);
            black_box(result);
        });
    });
}

criterion_group!(benches, bench_simplify);
criterion_main!(benches);
