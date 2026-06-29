use criterion::{Criterion, black_box, criterion_group, criterion_main};

#[cfg(feature = "gmp")]
use ocas_core::gmp::GmpInteger;

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

#[cfg(feature = "gmp")]
criterion_group!(benches, gmp_add, gmp_mul);

#[cfg(not(feature = "gmp"))]
fn placeholder(c: &mut Criterion) {
    c.bench_function("gmp_disabled", |bench| {
        bench.iter(|| black_box(1 + 1));
    });
}

#[cfg(not(feature = "gmp"))]
criterion_group!(benches, placeholder);

criterion_main!(benches);
