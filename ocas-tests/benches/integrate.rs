use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_core::arena::Arena;
use std::hint::black_box;

fn bench_integrate_rational(c: &mut Criterion) {
    c.bench_function("integrate_rational_log", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("(2*x + 3)/(x^2 + 3*x + 5)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });

    c.bench_function("integrate_rational_atan", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("1/(x^2 + 1)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });
}

fn bench_integrate_risch(c: &mut Criterion) {
    c.bench_function("integrate_log_x", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("log(x)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });

    c.bench_function("integrate_x_exp_x", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x*exp(x)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });

    c.bench_function("integrate_x2_exp_x", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x^2*exp(x)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });
}

fn bench_integrate_special(c: &mut Criterion) {
    c.bench_function("integrate_exp_neg_x2_erf", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("exp(-x^2)")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });

    c.bench_function("integrate_exp_over_x_ei", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("exp(x)/x")).unwrap();
            let result = integrate(&ctx, expr, Symbol::new("x"));
            black_box(result);
        });
    });
}

criterion_group!(
    benches,
    bench_integrate_rational,
    bench_integrate_risch,
    bench_integrate_special
);
criterion_main!(benches);
