use criterion::{Criterion, black_box, criterion_group, criterion_main};
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::SparseMultivariatePolynomial;

fn build_sparse(
    degree: usize,
    n_terms: usize,
) -> SparseMultivariatePolynomial<IntegerDomain, ocas_poly::Grevlex> {
    let domain = IntegerDomain;
    let mut terms = Vec::new();
    for i in 0..n_terms {
        let mut exp = vec![0usize; 3];
        exp[0] = i % (degree + 1);
        exp[1] = (i / (degree + 1)) % (degree + 1);
        terms.push((exp, Integer::from((i as i64) + 1)));
    }
    SparseMultivariatePolynomial::from_terms(domain, 3, terms)
}

fn sparse_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_mul");
    for (degree, n_terms) in [(5, 10), (5, 50), (10, 50), (10, 100)].iter() {
        let a = build_sparse(*degree, *n_terms);
        let b = build_sparse(*degree, *n_terms);
        let label = format!("degree_{}_terms_{}", degree, n_terms);
        group.bench_with_input(&label, &(*degree, *n_terms), |bench, _| {
            bench.iter(|| {
                let result = black_box(&a).mul(black_box(&b));
                black_box(result);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, sparse_mul);
criterion_main!(benches);
