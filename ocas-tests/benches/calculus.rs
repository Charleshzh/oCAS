use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_core::arena::Arena;

fn bench_diff(c: &mut Criterion) {
    c.bench_function("diff_polynomial", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x^4 + 3*x^3 + 2*x^2 + x + 5")).unwrap();
            let result = diff(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });

    c.bench_function("diff_nested", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("sin(x^2 + 1) * exp(3*x)")).unwrap();
            let result = diff(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });
}

fn bench_taylor(c: &mut Criterion) {
    c.bench_function("taylor_exp_order_5", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("exp(x)")).unwrap();
            let result = taylor(&ctx, expr, Symbol::new("x"), ctx.num(0), 5);
            black_box(result);
        });
    });
}

criterion_group!(benches, bench_diff, bench_taylor);
criterion_main!(benches);
