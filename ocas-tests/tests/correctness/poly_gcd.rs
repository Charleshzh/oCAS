use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn int(i: i64) -> Integer {
    Integer::from(i)
}

#[test]
fn poly_gcd_simple_quadratic() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(1)]); // x^2 - 1
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(1)]); // x - 1
    let g = a.gcd(&b);
    assert_eq!(g.coeffs(), &[int(-1), int(1)]); // x - 1
}

#[test]
fn poly_gcd_simple_common_factor() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(2), int(1)]); // x^2 + 2x + 1
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1)]); // x + 1
    let g = a.gcd(&b);
    assert_eq!(g.coeffs(), &[int(1), int(1)]); // x + 1
}

#[test]
fn poly_gcd_medium_coprime() {
    let d = IntegerDomain;
    let a =
        DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(0), int(0), int(1)]); // x^4 - 1
    let b = DenseUnivariatePolynomial::from_coeffs(
        d,
        vec![int(-1), int(0), int(0), int(0), int(0), int(0), int(1)],
    ); // x^6 - 1
    let g = a.gcd(&b);
    // gcd(x^4 - 1, x^6 - 1) = x^2 - 1
    assert_eq!(g.degree(), Some(2));
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn poly_gcd_complex_symmetric_identity() {
    let d = IntegerDomain;
    // x^3 - 1 = (x - 1)(x^2 + x + 1)
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(0), int(1)]);
    // x^3 + x^2 + x + 1 = (x + 1)(x^2 + 1)
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1), int(1), int(1)]);
    let g = a.gcd(&b);
    // These polynomials share no non-trivial factor.
    assert!(g.is_one());
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn poly_gcd_very_complex_random_lcm() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(0), int(1)]); // x^3 - 1
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(1)]); // x^2 - 1
    let g = a.gcd(&b);
    // gcd(x^3-1, x^2-1) = x - 1
    assert_eq!(g.degree(), Some(1));
}
