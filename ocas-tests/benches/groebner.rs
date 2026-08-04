//! Benchmark: Gröbner basis computation.
//!
//! Compares Buchberger, F4, F5 and the multi-modular pipeline on cyclic-n
//! and Katsura-n ideals over ℚ and ℤ_p.
//!
//! Reference: Symbolica's `groebner_basis.rs` example, Faugère F4 (1999).

use criterion::{Criterion, criterion_group, criterion_main};
use ocas_poly::buchberger;
use ocas_poly::groebner::f4::f4;
use ocas_poly::groebner::f5::f5;
use ocas_poly::groebner::multi_modular::groebner_basis_multi_modular;
use ocas_tests::systems::*;
use std::hint::black_box;

// =========================================================================
//  Benchmarks
// =========================================================================

/// Buchberger on cyclic-n over ℚ (only cyclic-3; cyclic-4 is too slow).
fn bench_buchberger_cyclic_q(c: &mut Criterion) {
    let mut group = c.benchmark_group("buchberger_cyclic_q");
    {
        let n = 3;
        let ideal = cyclic_q(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = buchberger(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℚ (generic path).
fn bench_f4_cyclic_q(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_q");
    for n in [3, 4] {
        let ideal = cyclic_q(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℤ_13 (fast i64 path).
fn bench_f4_cyclic_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_fp13");
    group.sample_size(10);
    for n in [3, 4, 5] {
        let ideal = cyclic_fp(n, 13);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F4 on cyclic-n over ℤ_101 (larger prime).
fn bench_f4_cyclic_fp101(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_cyclic_fp101");
    group.sample_size(10);
    for n in [3, 4, 5] {
        let ideal = cyclic_fp(n, 101);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f4(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F5 on cyclic-n over ℤ_13 (fast i64 path). cyclic-6 is the 0.25.0
/// acceptance gate (< 0.5 s median in release).
fn bench_f5_cyclic_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("f5_cyclic_fp13");
    group.sample_size(10);
    for n in [3, 4, 5, 6] {
        let ideal = cyclic_fp(n, 13);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f5(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F5 on cyclic-n over ℤ_13 under Grevlex. msolve's cyclic-n reference
/// numbers (cyclic-6 0.04 s) are measured under grevlex, so this group is
/// the msolve-aligned benchmark since 0.26.0; the Lex group remains for
/// legacy comparison.
fn bench_f5_cyclic_fp_grevlex(c: &mut Criterion) {
    let mut group = c.benchmark_group("f5_cyclic_fp13_grevlex");
    group.sample_size(10);
    for n in [5, 6] {
        let ideal = cyclic_fp_grevlex(n, 13);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f5(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// Multi-modular pipeline on cyclic-n over ℚ.
fn bench_multi_modular_cyclic_q(c: &mut Criterion) {
    let mut group = c.benchmark_group("multi_modular_cyclic_q");
    group.sample_size(10);
    for n in [3, 4, 5] {
        let ideal = cyclic_q(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = groebner_basis_multi_modular(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

/// F5 on Katsura-n over ℤ_13 (reference benchmark; no hard time bound).
fn bench_katsura_fp(c: &mut Criterion) {
    let mut group = c.benchmark_group("katsura_fp13");
    group.sample_size(10);
    for n in [6, 7] {
        let ideal = katsura_fp(n, 13);
        group.bench_with_input(format!("katsura_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = f5(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

// =========================================================================
//  Criterion harness
// =========================================================================

criterion_group!(
    benches,
    bench_buchberger_cyclic_q,
    bench_f4_cyclic_q,
    bench_f4_cyclic_fp,
    bench_f4_cyclic_fp101,
    bench_f5_cyclic_fp,
    bench_f5_cyclic_fp_grevlex,
    bench_multi_modular_cyclic_q,
    bench_katsura_fp,
);
criterion_main!(benches);
