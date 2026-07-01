use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_core::arena::Arena;

fn parse_small(c: &mut Criterion) {
    c.bench_function("parse_small", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box("x^2 + 2*x + 1")).unwrap();
            black_box(expr);
        });
    });
}

fn parse_medium(c: &mut Criterion) {
    let input = "(a + b + c + d + e)^2 * (x + y + z) - 3*x*y*z";
    c.bench_function("parse_medium", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box(input)).unwrap();
            black_box(expr);
        });
    });
}

fn parse_large(c: &mut Criterion) {
    let terms: Vec<String> = (0..50)
        .map(|i| format!("x{}^{} + {}", i, i % 5 + 1, i))
        .collect();
    let input = terms.join(" + ");
    c.bench_function("parse_large", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let expr = parse(&ctx, black_box(&input)).unwrap();
            black_box(expr);
        });
    });
}

fn parse_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_throughput");
    for depth in [2, 4, 6].iter() {
        let input = build_nested_expr(*depth, 3);
        group.bench_with_input(format!("depth_{}", depth), depth, |b, _| {
            b.iter(|| {
                let arena = Arena::new();
                let ctx = AtomArena::new(&arena);
                let expr = parse(&ctx, black_box(&input)).unwrap();
                black_box(expr);
            });
        });
    }
    group.finish();
}

fn build_nested_expr(depth: usize, width: usize) -> String {
    if depth == 0 {
        "x".to_string()
    } else {
        let inner = build_nested_expr(depth - 1, width);
        let parts: Vec<String> = (0..width).map(|i| format!("{}^{}", inner, i + 1)).collect();
        format!("({})", parts.join(" + "))
    }
}

criterion_group!(
    benches,
    parse_small,
    parse_medium,
    parse_large,
    parse_throughput
);
criterion_main!(benches);
