use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn int(i: i64) -> Integer {
    Integer::from(i)
}

#[test]
fn poly_arithmetic_simple_add() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1)]); // x + 1
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(2), int(1)]); // x + 2
    let c = a.add(&b);
    assert_eq!(c.coeffs(), &[int(3), int(2)]); // 2*x + 3
}

#[test]
fn poly_arithmetic_simple_mul() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1)]); // x + 1
    let c = a.mul(&a);
    assert_eq!(c.coeffs(), &[int(1), int(2), int(1)]); // x^2 + 2x + 1
}

#[test]
fn poly_arithmetic_simple_div() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(0), int(1)]); // x^2 - 1
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-1), int(1)]); // x - 1
    let (q, r) = a.div_rem(&b).unwrap();
    assert_eq!(q.coeffs(), &[int(1), int(1)]); // x + 1
    assert!(r.is_zero());
}

#[test]
fn poly_arithmetic_medium_cubic_times_quadratic() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(4), int(3), int(2), int(1)]); // x^3 + 2x^2 + 3x + 4
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(6), int(5), int(1)]); // x^2 + 5x + 6
    let c = a.mul(&b);
    assert_eq!(c.degree(), Some(5));
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn poly_arithmetic_complex_random_roundtrip() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(0), int(-1)]);
    let b = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1)]);
    let (q, r) = a.div_rem(&b).unwrap();
    // a = b*q + r
    let reconstructed = b.mul(&q).add(&r);
    assert_eq!(a.coeffs(), reconstructed.coeffs());
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn poly_arithmetic_very_complex_high_degree_multiplication() {
    let d = IntegerDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1), int(1)]); // x+1
    let mut power = a.clone();
    for _ in 0..9 {
        power = power.mul(&a);
    }
    // (x+1)^10 has degree 10
    assert_eq!(power.degree(), Some(10));
    let last = power.leading_coeff().cloned().unwrap();
    assert_eq!(last, int(1));
}
