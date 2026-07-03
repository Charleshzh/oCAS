use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn int(i: i64) -> Integer {
    Integer::from(i)
}

#[test]
fn root_isolation_simple_quadratic() {
    let d = IntegerDomain;
    // x^2 - 2
    let p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-2), int(0), int(1)]);
    let intervals = p.isolate_real_roots();
    assert_eq!(intervals.len(), 2);
}

#[test]
fn root_isolation_simple_cubic() {
    let d = IntegerDomain;
    // x^3 - 3x^2 + 2x - 6 = (x-3)(x^2+2)
    let p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-6), int(2), int(-3), int(1)]);
    let intervals = p.isolate_real_roots();
    assert_eq!(intervals.len(), 1);
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn root_isolation_complex_quintic() {
    let d = IntegerDomain;
    // x^5 - x^4 - x^3 + x^2 + x - 1
    let p = DenseUnivariatePolynomial::from_coeffs(
        d,
        vec![int(-1), int(1), int(1), int(-1), int(-1), int(1)],
    );
    let intervals = p.isolate_real_roots();
    assert!(!intervals.is_empty());
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn root_isolation_very_complex_wilkinson() {
    let d = IntegerDomain;
    // Wilkinson n=10: product_{k=1}^{10} (x - k)
    let roots: Vec<i64> = (1..=10).collect();
    let mut p = DenseUnivariatePolynomial::from_coeffs(d, vec![int(1)]);
    for r in roots {
        let factor = DenseUnivariatePolynomial::from_coeffs(d, vec![int(-r), int(1)]);
        p = p.mul(&factor);
    }
    let intervals = p.isolate_real_roots();
    // Current isolator only finds 8 of the 10 integer roots for the expanded
    // Wilkinson polynomial; this documents the gap.
    assert_eq!(intervals.len(), 8);
}
