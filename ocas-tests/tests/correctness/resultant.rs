use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn int(i: i64) -> Integer {
    Integer::from(i)
}

fn poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<IntegerDomain> {
    DenseUnivariatePolynomial::from_coeffs(IntegerDomain, coeffs.iter().map(|&c| int(c)).collect())
}

#[test]
fn resultant_different_linear_roots() {
    // Res(x - 1, x - 2) = 1 - 2 = -1
    let a = poly(&[-1, 1]); // x - 1
    let b = poly(&[-2, 1]); // x - 2
    assert_eq!(a.resultant(&b), int(-1));
}

#[test]
fn resultant_common_root() {
    // Res(x^2 - 1, x - 1) = 0 (share root x=1)
    let a = poly(&[-1, 0, 1]); // x^2 - 1
    let b = poly(&[-1, 1]); // x - 1
    assert_eq!(a.resultant(&b), int(0));
}

#[test]
fn resultant_no_common_root_quadratic() {
    // Res(x^2 + 1, (x+1)^2) = 4
    let a = poly(&[1, 0, 1]); // x^2 + 1
    let b = poly(&[1, 2, 1]); // x^2 + 2x + 1
    assert_eq!(a.resultant(&b), int(4));
}

#[test]
fn resultant_shared_factor() {
    // Res((x-1)(x-2), (x-1)(x-3)) = 0
    let a = poly(&[2, -3, 1]); // x^2 - 3x + 2
    let b = poly(&[3, -4, 1]); // x^2 - 4x + 3
    assert_eq!(a.resultant(&b), int(0));
}

#[test]
fn resultant_constant_poly() {
    // Res(x^2 + 1, 3) = 3^2 = 9
    let a = poly(&[1, 0, 1]);
    let b = poly(&[3]);
    assert_eq!(a.resultant(&b), int(9));
}

#[test]
fn resultant_symmetric_up_to_sign() {
    // Res(a, b) = (-1)^(deg_a * deg_b) * Res(b, a)
    let a = poly(&[-1, 0, 1]); // x^2 - 1, deg=2
    let b = poly(&[-2, 1]); // x - 2, deg=1
    let r1 = a.resultant(&b);
    let r2 = b.resultant(&a);
    // deg_a * deg_b = 2, even, so Res(a,b) = Res(b,a)
    assert_eq!(r1, r2);
}

#[test]
fn resultant_zero_polynomial() {
    let a = poly(&[0]); // zero polynomial
    let b = poly(&[1, 1]); // x + 1
    assert_eq!(a.resultant(&b), int(0));
}

#[test]
fn resultant_cubic_no_common_root() {
    // Res(x^3 - 1, x^3 + 1) = (-1)^3 * (1^3 - (-1)^3) = ... actually
    // x^3 - 1 = (x-1)(x^2+x+1), x^3 + 1 = (x+1)(x^2-x+1)
    // No common roots, so resultant != 0
    let a = poly(&[-1, 0, 0, 1]); // x^3 - 1
    let b = poly(&[1, 0, 0, 1]); // x^3 + 1
    let r = a.resultant(&b);
    assert_ne!(r, int(0));
}
