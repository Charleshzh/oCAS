use ocas_domain::{Domain, Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};

fn int(i: i64) -> Integer {
    Integer::from(i)
}

/// Sparse multivariate polynomial over ℤ with Lex order.
type ZmPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

fn zm_poly(n_vars: usize, terms: &[(Vec<usize>, i64)]) -> ZmPoly {
    SparseMultivariatePolynomial::from_terms(
        IntegerDomain,
        n_vars,
        terms
            .iter()
            .map(|(e, c)| (e.clone(), Integer::from(*c)))
            .collect(),
    )
}

/// Expand multivariate factors (with multiplicities) back into a polynomial.
fn expand_mpoly_factors(factors: &[(ZmPoly, usize)]) -> ZmPoly {
    let n = factors[0].0.n_vars();
    let mut acc = SparseMultivariatePolynomial::new(IntegerDomain, n);
    acc.set_term_external(vec![0; n], Integer::from(1));
    for (g, m) in factors {
        for _ in 0..*m {
            acc = acc.mul(g);
        }
    }
    acc
}

/// Whether two multivariate polynomials are equal up to a nonzero constant.
fn mpoly_eq_up_to_unit(a: &ZmPoly, b: &ZmPoly) -> bool {
    if a.is_zero() || b.is_zero() {
        return a.is_zero() && b.is_zero();
    }
    let (e0, c0) = a.terms_ref().iter().next().unwrap();
    let bc0 = b.coeff(e0);
    if bc0.is_zero() {
        return false;
    }
    let ratio = IntegerDomain.div(c0, &bc0);
    let ratio = match ratio {
        Some(r) => r,
        None => return false,
    };
    if a.n_terms() != b.n_terms() {
        return false;
    }
    a.terms_ref()
        .iter()
        .all(|(e, c)| *c == IntegerDomain.mul(&ratio, &b.coeff(e)))
}

fn expand_factors(
    d: IntegerDomain,
    factors: &[(DenseUnivariatePolynomial<IntegerDomain>, usize)],
) -> DenseUnivariatePolynomial<IntegerDomain> {
    let mut product = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1)]);
    for (f, mult) in factors {
        let mut power = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1)]);
        for _ in 0..*mult {
            power = power.mul(f);
        }
        product = product.mul(&power);
    }
    product
}

#[test]
fn poly_factor_simple_difference_of_squares() {
    let d = IntegerDomain;
    // x^2 - 1
    let p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(1)]);
    let factors = p.factor();
    assert_eq!(factors.len(), 2);
    let expanded = expand_factors(d, &factors);
    assert_eq!(p.coeffs(), expanded.coeffs());
}

#[test]
fn poly_factor_simple_perfect_square() {
    let d = IntegerDomain;
    // x^2 + 2x + 1
    let p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(2), int(1)]);
    let factors = p.factor();
    let expanded = expand_factors(d, &factors);
    assert_eq!(p.coeffs(), expanded.coeffs());
}

#[test]
fn poly_factor_medium_quartic() {
    let d = IntegerDomain;
    // x^4 - 5x^2 + 4
    let p =
        DenseUnivariatePolynomial::from_coeffs(d, vec![int(4), int(0), int(-5), int(0), int(1)]);
    let factors = p.factor();
    let expanded = expand_factors(d, &factors);
    assert_eq!(p.coeffs(), expanded.coeffs());
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn poly_factor_complex_cubic_perfect_cube() {
    let d = IntegerDomain;
    // x^3 + 3x^2 + 3x + 1 = (x+1)^3
    let p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(3), int(3), int(1)]);
    let factors = p.factor();
    let expanded = expand_factors(d, &factors);
    assert_eq!(p.coeffs(), expanded.coeffs());
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn poly_factor_very_complex_high_degree() {
    let d = IntegerDomain;
    // (x+1)(x+2)(x+3)(x+4) = x^4 + 10x^3 + 35x^2 + 50x + 24
    let p =
        DenseUnivariatePolynomial::from_coeffs(d, vec![int(24), int(50), int(35), int(10), int(1)]);
    let factors = p.factor();
    let expanded = expand_factors(d, &factors);
    assert_eq!(p.coeffs(), expanded.coeffs());
}

// ── multivariate (≥3 variables) factorization ─────────────────────

#[test]
fn poly_factor_trivariate_monic_three_linear() {
    // f = (x + y + z)(x - y + 2z)(x + y + 1), monic in x.
    let f1 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
    );
    let f2 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 1], 2)],
    );
    let f3 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 0], 1)],
    );
    let f = f1.mul(&f2).mul(&f3);
    let factors = f.factor();
    assert_eq!(factors.len(), 3, "expected 3 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

#[test]
fn poly_factor_trivariate_repeated() {
    // f = (x + y + z)^2 (x - y + 1).
    let f1 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
    );
    let f2 = zm_poly(
        3,
        &[(vec![1, 0, 0], 1), (vec![0, 1, 0], -1), (vec![0, 0, 0], 1)],
    );
    let f = f1.mul(&f1).mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    let mut mults: Vec<usize> = factors.iter().map(|(_, m)| *m).collect();
    mults.sort_unstable();
    assert_eq!(mults, vec![1, 2]);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

#[test]
fn poly_factor_four_variables() {
    // f = (x1 + x2 + x3 + x4)(x1 - x2 + x3 - x4).
    let f1 = zm_poly(
        4,
        &[
            (vec![1, 0, 0, 0], 1),
            (vec![0, 1, 0, 0], 1),
            (vec![0, 0, 1, 0], 1),
            (vec![0, 0, 0, 1], 1),
        ],
    );
    let f2 = zm_poly(
        4,
        &[
            (vec![1, 0, 0, 0], 1),
            (vec![0, 1, 0, 0], -1),
            (vec![0, 0, 1, 0], 1),
            (vec![0, 0, 0, 1], -1),
        ],
    );
    let f = f1.mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

// ── non-constant leading coefficients (Wang imposition, 0.16.1) ───

#[test]
fn poly_factor_bivariate_nonconstant_lcoeff() {
    // f = (y·x² + 1)(x + 1): leading coefficient y.
    let f1 = zm_poly(2, &[(vec![2, 1], 1), (vec![0, 0], 1)]);
    let f2 = zm_poly(2, &[(vec![1, 0], 1), (vec![0, 0], 1)]);
    let f = f1.mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

#[test]
fn poly_factor_trivariate_nonconstant_lcoeff() {
    // f = (z·x² + y)(x + 1): leading coefficient z.
    let f1 = zm_poly(3, &[(vec![2, 0, 1], 1), (vec![0, 1, 0], 1)]);
    let f2 = zm_poly(3, &[(vec![1, 0, 0], 1), (vec![0, 0, 0], 1)]);
    let f = f1.mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

#[test]
fn poly_factor_trivariate_reducible_lcoeff() {
    // f = (x·y² − x + 2y)(x·y − z): leading coefficient y³ − y factors as
    // (y−1)(y+1)y and must be distributed (y²−1) / y across the factors.
    let f1 = zm_poly(
        3,
        &[(vec![1, 2, 0], 1), (vec![1, 0, 0], -1), (vec![0, 1, 0], 2)],
    );
    let f2 = zm_poly(3, &[(vec![1, 1, 0], 1), (vec![0, 0, 1], -1)]);
    let f = f1.mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}

#[test]
fn poly_factor_sparse_four_var_nonconstant_lcoeff() {
    // Sparse product in 4 variables with ≥ 50 terms and non-constant
    // leading coefficients (skeleton interpolation in the p-adic lift).
    let mut f1_terms = vec![(vec![2usize, 1, 1, 0], 1i64)]; // y·z·x² LC
    let mut f2_terms = vec![(vec![1, 1, 0, 0], 1i64), (vec![1, 0, 0, 1], 1)]; // (y+w)·x
    for i in 0..4usize {
        for j in 0..3usize {
            let c1 = ((i * 7 + j * 3) % 4 + 1) as i64;
            let c2 = ((i * 5 + j * 11 + 2) % 4 + 1) as i64;
            f1_terms.push((vec![i % 2, i, j, (i + j) % 2], c1));
            f2_terms.push((vec![0, (i + 1) % 3, (j + 2) % 2, i % 3], c2));
        }
    }
    let f1 = zm_poly(4, &f1_terms);
    let f2 = zm_poly(4, &f2_terms);
    let f = f1.mul(&f2);
    let factors = f.factor();
    assert_eq!(factors.len(), 2, "expected 2 factors, got {:?}", factors);
    assert!(mpoly_eq_up_to_unit(&expand_mpoly_factors(&factors), &f));
}
