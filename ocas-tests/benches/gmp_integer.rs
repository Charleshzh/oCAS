use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

#[cfg(feature = "gmp")]
use ocas_core::gmp::GmpInteger;

#[cfg(feature = "gmp")]
use ocas_domain::{Domain, Integer, IntegerDomain};

// ---------------------------------------------------------------------------
// Low-level GmpInteger benchmarks (ocas-core wrapper)
// ---------------------------------------------------------------------------

#[cfg(feature = "gmp")]
fn gmp_add(c: &mut Criterion) {
    let a = GmpInteger::from_i64(1234567890123456789);
    let b = GmpInteger::from_i64(9876543210987654321);
    c.bench_function("gmp_add", |bench| {
        bench.iter(|| {
            let result = a.add(&b);
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn gmp_mul(c: &mut Criterion) {
    let a = GmpInteger::from_i64(1234567890123456789);
    let b = GmpInteger::from_i64(9876543210987654321);
    c.bench_function("gmp_mul", |bench| {
        bench.iter(|| {
            let result = a.mul(&b);
            black_box(result);
        });
    });
}

// ---------------------------------------------------------------------------
// SOO Integer benchmarks (ocas-domain Integer with Small/Large paths)
// ---------------------------------------------------------------------------

#[cfg(feature = "gmp")]
fn soo_small_add(c: &mut Criterion) {
    let domain = IntegerDomain;
    let a = Integer::from(42i64);
    let b = Integer::from(17i64);
    c.bench_function("soo_small_add", |bench| {
        bench.iter(|| {
            let result = domain.add(black_box(&a), black_box(&b));
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_small_mul(c: &mut Criterion) {
    let domain = IntegerDomain;
    let a = Integer::from(42i64);
    let b = Integer::from(17i64);
    c.bench_function("soo_small_mul", |bench| {
        bench.iter(|| {
            let result = domain.mul(black_box(&a), black_box(&b));
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_large_add(c: &mut Criterion) {
    let domain = IntegerDomain;
    let a = Integer::from(i64::MAX);
    let b = Integer::from(1i64);
    // First addition promotes to Large.
    let a_large = &a + &b;
    let b_large = Integer::from(i64::MAX);
    c.bench_function("soo_large_add", |bench| {
        bench.iter(|| {
            let result = domain.add(black_box(&a_large), black_box(&b_large));
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_large_mul(c: &mut Criterion) {
    let domain = IntegerDomain;
    let a = Integer::from(i64::MAX);
    let b = Integer::from(i64::MAX);
    c.bench_function("soo_large_mul", |bench| {
        bench.iter(|| {
            let result = domain.mul(black_box(&a), black_box(&b));
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_to_bigint_small(c: &mut Criterion) {
    let a = Integer::from(42i64);
    c.bench_function("soo_to_bigint_small", |bench| {
        bench.iter(|| {
            let result = black_box(&a).to_bigint();
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_to_bigint_large(c: &mut Criterion) {
    let a = Integer::from(i64::MAX) + Integer::from(1i64);
    c.bench_function("soo_to_bigint_large", |bench| {
        bench.iter(|| {
            let result = black_box(&a).to_bigint();
            black_box(result);
        });
    });
}

#[cfg(feature = "gmp")]
fn soo_is_zero(c: &mut Criterion) {
    let a = Integer::from(0i64);
    c.bench_function("soo_is_zero", |bench| {
        bench.iter(|| {
            let result = black_box(&a).is_zero();
            black_box(result);
        });
    });
}

// ---------------------------------------------------------------------------
// Grouping
// ---------------------------------------------------------------------------

#[cfg(feature = "gmp")]
criterion_group!(
    benches,
    gmp_add,
    gmp_mul,
    soo_small_add,
    soo_small_mul,
    soo_large_add,
    soo_large_mul,
    soo_to_bigint_small,
    soo_to_bigint_large,
    soo_is_zero,
);

#[cfg(not(feature = "gmp"))]
fn placeholder(c: &mut Criterion) {
    c.bench_function("gmp_disabled", |bench| {
        bench.iter(|| black_box(1 + 1));
    });
}

#[cfg(not(feature = "gmp"))]
criterion_group!(benches, placeholder);

criterion_main!(benches);
