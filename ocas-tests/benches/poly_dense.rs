use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn build_poly(degree: usize) -> DenseUnivariatePolynomial<IntegerDomain> {
    let domain = IntegerDomain;
    let coeffs: Vec<Integer> = (0..=degree).map(|i| Integer::from(i as i64)).collect();
    DenseUnivariatePolynomial::from_coeffs(domain, coeffs)
}

fn dense_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("dense_mul");
    for degree in [10, 50, 100, 500].iter() {
        let a = build_poly(*degree);
        let b = build_poly(*degree);
        group.bench_with_input(format!("degree_{}", degree), degree, |bench, _| {
            bench.iter(|| {
                let result = black_box(&a).mul(black_box(&b));
                black_box(result);
            });
        });
    }
    group.finish();
}

fn dense_eval(c: &mut Criterion) {
    let p = build_poly(100);
    let x = Integer::from(7);
    c.bench_function("dense_eval", |b| {
        b.iter(|| {
            let result = black_box(&p).eval(black_box(&x));
            black_box(result);
        });
    });
}

fn dense_div_rem(c: &mut Criterion) {
    let domain = IntegerDomain;
    let dividend = DenseUnivariatePolynomial::from_coeffs(
        domain,
        vec![
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(-1),
        ],
    );
    let divisor =
        DenseUnivariatePolynomial::from_coeffs(domain, vec![Integer::from(1), Integer::from(-1)]);
    c.bench_function("dense_div_rem", |b| {
        b.iter(|| {
            let result = black_box(&dividend).div_rem(black_box(&divisor));
            black_box(result);
        });
    });
}

criterion_group!(benches, dense_mul, dense_eval, dense_div_rem);
criterion_main!(benches);
