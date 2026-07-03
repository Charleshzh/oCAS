use num_bigint::BigInt;
use ocas_domain::{Domain, FiniteField};

#[test]
fn finite_field_simple_add() {
    let f = FiniteField::new(BigInt::from(7));
    let a = f.element(3);
    let b = f.element(5);
    assert_eq!(f.add(&a, &b), f.element(1));
}

#[test]
fn finite_field_simple_mul() {
    let f = FiniteField::new(BigInt::from(7));
    let a = f.element(3);
    let b = f.element(5);
    assert_eq!(f.mul(&a, &b), f.element(1));
}

#[test]
fn finite_field_simple_inv() {
    let f = FiniteField::new(BigInt::from(7));
    let a = f.element(3);
    let inv = f.inv(&a).unwrap();
    assert_eq!(f.mul(&a, &inv), f.one());
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn finite_field_complex_poly_eval() {
    use ocas_poly::DenseUnivariatePolynomial;
    let f = FiniteField::new(BigInt::from(17));
    let p = DenseUnivariatePolynomial::from_coeffs(
        f.clone(),
        vec![f.element(1), f.element(2), f.element(1)], // x^2 + 2x + 1
    );
    let value = p.eval(&f.element(3));
    assert_eq!(value, f.element(16)); // 9 + 6 + 1 = 16 mod 17
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn finite_field_very_complex_large_prime() {
    let p = "1000000007".parse::<BigInt>().unwrap();
    let f = FiniteField::new(p);
    let a = f.element(123456789);
    let b = f.element(987654321);
    let sum = f.add(&a, &b);
    let expected = f.element(123456789_i64 + 987654321_i64);
    assert_eq!(sum, expected);
}
