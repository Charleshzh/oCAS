use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn int(i: i64) -> Integer {
    Integer::from(i)
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
