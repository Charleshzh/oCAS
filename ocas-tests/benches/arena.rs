use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas_core::arena::Arena;

fn arena_alloc_small(c: &mut Criterion) {
    c.bench_function("arena_alloc_small", |b| {
        b.iter(|| {
            let arena = Arena::new();
            for i in 0..100 {
                let value = arena.allocate_with(|| black_box(i));
                black_box(value);
            }
        });
    });
}

fn arena_alloc_large(c: &mut Criterion) {
    c.bench_function("arena_alloc_large", |b| {
        b.iter(|| {
            let arena = Arena::with_capacity(256);
            let data = [0u8; 4096];
            let ptr = arena.allocate_with(|| black_box(data));
            black_box(ptr);
        });
    });
}

fn arena_alloc_many(c: &mut Criterion) {
    let mut group = c.benchmark_group("arena_alloc_many");
    for count in [10, 100, 1_000, 10_000].iter() {
        group.bench_with_input(format!("count_{}", count), count, |b, &count| {
            b.iter(|| {
                let arena = Arena::with_capacity(256);
                for i in 0..count {
                    let value = arena.allocate_with(|| black_box(i as i64));
                    black_box(value);
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    arena_alloc_small,
    arena_alloc_large,
    arena_alloc_many
);
criterion_main!(benches);
