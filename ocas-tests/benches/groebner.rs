//! Benchmark: Gröbner basis computation via Buchberger's algorithm.
//!
//! Measures [`GroebnerBasis::buchberger`] on the cyclic-n ideals, the
//! standard benchmark family for Gröbner basis algorithms (cf. Symbolica's
//! `groebner_basis.rs` example).

use criterion::{Criterion, criterion_group, criterion_main};
use ocas_domain::{Rational, RationalDomain};
use ocas_poly::GroebnerBasis;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::Lex;
use std::hint::black_box;

type RatPoly = SparseMultivariatePolynomial<RationalDomain, Lex>;

/// Build a single term `c * x0^e0 * x1^e1 * ...`.
fn term(exps: Vec<usize>, num: i64, den: i64) -> (Vec<usize>, Rational) {
    (exps, Rational::new(num, den))
}

/// The cyclic-n ideal generators.
///
/// cyclic-n is defined by:
///   σ_k = Σ_{i1<...<ik} x_{i1}...x_{ik}            for k = 1..n
///   f_k = σ_k - σ_{k+1} * (product of all vars?)   (standard variant)
///
/// We use the common textbook form:
///   f_1 = x1 + x2 + ... + xn
///   f_2 = x1*x2 + x2*x3 + ... + x(n-1)*xn + xn*x1  (cyclic pairs)
///   f_3 = x1*x2*x3 + x2*x3*x4 + ...                (cyclic triples)
///   ...
///   f_n = x1*x2*...*xn - 1
fn cyclic(n: usize) -> Vec<RatPoly> {
    let d = RationalDomain;
    let mut gens = Vec::with_capacity(n - 1);

    // f_k: sum of all cyclic products of length k, for k = 1..=n-1.
    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut exps = vec![0usize; n];
            for j in 0..k {
                exps[(start + j) % n] = 1;
            }
            terms.push(term(exps, 1, 1));
        }
        gens.push(RatPoly::from_terms(d, n, terms));
    }

    // f_n: product of all variables minus 1.
    let full_exps = vec![1usize; n];
    gens.push(RatPoly::from_terms(
        d,
        n,
        vec![term(full_exps, 1, 1), term(vec![0usize; n], -1, 1)],
    ));

    gens
}

fn groebner_cyclic(c: &mut Criterion) {
    let mut group = c.benchmark_group("groebner_cyclic");
    // cyclic-3 and cyclic-4 are the standard small benchmarks. cyclic-4 is
    // already non-trivial for a naive Buchberger implementation.
    for n in [3, 4] {
        let ideal = cyclic(n);
        group.bench_with_input(format!("cyclic_{n}"), &n, |bench, _| {
            bench.iter(|| {
                let gb = GroebnerBasis::buchberger(black_box(&ideal));
                black_box(gb);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, groebner_cyclic);
criterion_main!(benches);
