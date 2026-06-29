use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn polynomial_addition(c: &mut Criterion) {
    c.bench_function("polynomial_addition", |b| {
        b.iter(|| {
            let _ = black_box(1 + 1);
        });
    });
}

criterion_group!(benches, polynomial_addition);
criterion_main!(benches);
