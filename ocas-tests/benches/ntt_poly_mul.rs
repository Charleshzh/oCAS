//! Benchmark: NTT vs Karatsuba vs Schoolbook polynomial multiplication over ℤ_p.
//!
//! Compares the three multiplication strategies for `DenseUnivariatePolynomial<FiniteField>`
//! at various degrees. Uses the NTT-friendly prime p = 998244353.
//!
//! Run with: `cargo bench --bench ntt_poly_mul --features ntt`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use num_bigint::BigInt;
use ocas_domain::FiniteField;
use ocas_poly::DenseUnivariatePolynomial;
use std::hint::black_box;

/// NTT-friendly prime: 998244353 = 119 * 2^23 + 1
const P: u64 = 998244353;

fn build_fp_poly(degree: usize, prime: u64) -> DenseUnivariatePolynomial<FiniteField> {
    let field = FiniteField::new(BigInt::from(prime));
    let coeffs: Vec<_> = (0..=degree)
        .map(|i| field.element(BigInt::from((i as u64 * 7 + 13) % prime)))
        .collect();
    DenseUnivariatePolynomial::from_coeffs(field, coeffs)
}

fn bench_ntt_vs_schoolbook(c: &mut Criterion) {
    let mut group = c.benchmark_group("fp_mul_ntt_vs_schoolbook");

    for degree in [64, 128, 256, 512, 1024, 2048] {
        // Schoolbook vs NTT on raw u64 coefficients for a fair comparison
        let a_u64: Vec<u64> = (0..=degree).map(|i| (i as u64 * 7 + 13) % P).collect();
        let b_u64: Vec<u64> = (0..=degree).map(|i| (i as u64 * 11 + 5) % P).collect();

        // Schoolbook O(n^2)
        group.bench_with_input(
            BenchmarkId::new("schoolbook", degree),
            &degree,
            |bench, _| {
                bench.iter(|| {
                    let result_len = a_u64.len() + b_u64.len() - 1;
                    let mut result = vec![0u64; result_len];
                    for (i, &ai) in a_u64.iter().enumerate() {
                        for (j, &bj) in b_u64.iter().enumerate() {
                            result[i + j] = (result[i + j]
                                + ((ai as u128 * bj as u128) % P as u128) as u64)
                                % P;
                        }
                    }
                    black_box(result);
                });
            },
        );

        // NTT O(n log n)
        group.bench_with_input(BenchmarkId::new("ntt", degree), &degree, |bench, _| {
            bench.iter(|| {
                let result = ocas_poly::ntt::ntt_mul(black_box(&a_u64), black_box(&b_u64), P);
                black_box(result);
            });
        });
    }
    group.finish();
}

fn bench_ntt_forward(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_forward");

    for n in [256, 512, 1024, 2048, 4096] {
        let root = ocas_poly::ntt::find_primitive_root(P, n).unwrap();
        let data: Vec<u64> = (0..n).map(|i| (i as u64 * 31 + 17) % P).collect();

        group.bench_with_input(BenchmarkId::new("n", n), &n, |bench, _| {
            bench.iter(|| {
                let mut buf = data.clone();
                ocas_poly::ntt::ntt_forward(&mut buf, root, P);
                black_box(buf);
            });
        });
    }
    group.finish();
}

fn bench_fp_poly_mul_via_api(c: &mut Criterion) {
    let mut group = c.benchmark_group("fp_poly_mul_api");

    for degree in [64, 256, 512, 1024] {
        let a = build_fp_poly(degree, P);
        let b = build_fp_poly(degree, P);

        // Generic mul (Karatsuba/Schoolbook)
        group.bench_with_input(
            BenchmarkId::new("karatsuba_generic", degree),
            &degree,
            |bench, _| {
                bench.iter(|| {
                    let result = black_box(&a).mul(black_box(&b));
                    black_box(result);
                });
            },
        );

        // NTT-accelerated mul_ntt
        group.bench_with_input(
            BenchmarkId::new("ntt_mul_ntt", degree),
            &degree,
            |bench, _| {
                bench.iter(|| {
                    let mut buf = Vec::new();
                    black_box(&a).mul_ntt(black_box(&b), &mut buf);
                    black_box(buf);
                });
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_ntt_vs_schoolbook,
    bench_ntt_forward,
    bench_fp_poly_mul_via_api
);
criterion_main!(benches);
