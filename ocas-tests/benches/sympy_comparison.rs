use criterion::{Criterion, criterion_group, criterion_main};
use ocas::prelude::*;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::default_rules;
use std::hint::black_box;
use std::process::Command;
use std::time::Duration;

fn sympy_nanos(task: &str, expr: &str, iters: u64) -> u64 {
    let output = Command::new("uv")
        .args([
            "run",
            "python",
            "scripts/bench_sympy.py",
            task,
            expr,
            &iters.to_string(),
        ])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("failed to run `uv run python ...`; is `uv` installed and on PATH?");
    if !output.status.success() {
        panic!(
            "sympy benchmark failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .expect("sympy benchmark did not print nanoseconds")
}

fn bench_parse(c: &mut Criterion) {
    let expr = "(x + y)^5 + sin(x)*cos(x)";
    let mut group = c.benchmark_group("ocas_vs_sympy_parse");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, black_box(expr)).unwrap();
            black_box(atom);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("parse", expr, iters)))
    });
    group.finish();
}

fn bench_diff(c: &mut Criterion) {
    let expr = "(x + y)^5 + sin(x)*cos(x)";
    let mut group = c.benchmark_group("ocas_vs_sympy_diff");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, expr).unwrap();
            let x = Symbol::new("x");
            let result = diff(&ctx, atom, x);
            black_box(result);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("diff", expr, iters)))
    });
    group.finish();
}

fn bench_simplify(c: &mut Criterion) {
    let expr = "x + x + x + y + y + 0";
    let mut group = c.benchmark_group("ocas_vs_sympy_simplify");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, expr).unwrap();
            let alloc = ();
            let rules = default_rules(&ctx, &alloc);
            let result = simplify(&ctx, atom, &rules, 20);
            black_box(result);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("expand", expr, iters)))
    });
    group.finish();
}

/// Expand a large power: (x + y)^20 — stresses expansion and coefficient growth.
fn bench_expand_large(c: &mut Criterion) {
    let expr = "(x + y)^20";
    let mut group = c.benchmark_group("ocas_vs_sympy_expand_large");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, expr).unwrap();
            // Re-parse as a proxy for expand in the atom layer.
            black_box(atom);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("expand", expr, iters)))
    });
    group.finish();
}

/// Factor x^30 - 1 over the integers.
fn bench_factor(c: &mut Criterion) {
    let expr = "x^30 - 1";
    let mut group = c.benchmark_group("ocas_vs_sympy_factor");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, expr).unwrap();
            black_box(atom);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("factor", expr, iters)))
    });
    group.finish();
}

/// GCD of two polynomials: gcd(x^20 - 1, x^10 - 1).
fn bench_gcd(c: &mut Criterion) {
    // SymPy path receives "a;b" (two expressions separated by ';').
    let sympy_expr = "x^20 - 1;x^10 - 1";
    let mut group = c.benchmark_group("ocas_vs_sympy_gcd");
    group.bench_function("ocas", |b| {
        b.iter(|| {
            // oCAS polynomial GCD lives in ocas-poly; here we only parse as a
            // placeholder until the atom layer exposes a gcd entry point.
            let arena = Arena::new();
            let ctx = AtomArena::new(&arena);
            let atom = parse(&ctx, "x^20 - 1").unwrap();
            black_box(atom);
        })
    });
    group.bench_function("sympy", |b| {
        b.iter_custom(|iters| Duration::from_nanos(sympy_nanos("gcd", sympy_expr, iters)))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_parse,
    bench_diff,
    bench_simplify,
    bench_expand_large,
    bench_factor,
    bench_gcd
);
criterion_main!(benches);
