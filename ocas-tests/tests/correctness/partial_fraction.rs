use ocas_domain::{Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;
use ocas_calc::partial_fraction::{apart, together, PartialFractionTerm};

fn rat_poly(coeffs: &[i64]) -> DenseUnivariatePolynomial<RationalDomain> {
    DenseUnivariatePolynomial::from_coeffs(
        RationalDomain,
        coeffs.iter().map(|&c| Rational::new(c, 1)).collect(),
    )
}

fn rat_poly_r(coeffs: &[(i64, i64)]) -> DenseUnivariatePolynomial<RationalDomain> {
    DenseUnivariatePolynomial::from_coeffs(
        RationalDomain,
        coeffs.iter().map(|&(n, d)| Rational::new(n, d)).collect(),
    )
}

#[test]
fn apart_proper_fraction_single_factor() {
    // 1 / (x^2 - 1) — square-free, single factor from square_free_factorization
    let num = rat_poly(&[1]);
    let den = rat_poly(&[-1, 0, 1]); // x^2 - 1
    let (poly_part, terms) = apart(&num, &den);
    assert!(poly_part.is_none());
    assert_eq!(terms.len(), 1);
    assert_eq!(terms[0].exp, 1);
}

#[test]
fn apart_with_polynomial_part() {
    // (x^2 + 1) / (x - 1) = (x + 1) + 2/(x-1)
    let num = rat_poly(&[1, 0, 1]); // x^2 + 1
    let den = rat_poly(&[-1, 1]); // x - 1
    let (poly_part, terms) = apart(&num, &den);
    assert!(poly_part.is_some());
    let pp = poly_part.unwrap();
    assert_eq!(pp.degree(), Some(1));
    assert_eq!(pp.coeff(0), Some(&Rational::new(1, 1)));
    assert_eq!(pp.coeff(1), Some(&Rational::new(1, 1)));
    assert_eq!(terms.len(), 1);
}

#[test]
fn apart_repeated_factor() {
    // 1 / (x-1)^2 — should give p-adic expansion terms
    let num = rat_poly(&[1]);
    let den = rat_poly(&[1, -2, 1]); // (x-1)^2
    let (poly_part, terms) = apart(&num, &den);
    assert!(poly_part.is_none());
    assert!(!terms.is_empty());
}

#[test]
fn apart_trivial_polynomial() {
    // x / x = 1 (no remainder)
    let num = rat_poly(&[0, 1]);
    let den = rat_poly(&[0, 1]);
    let (poly_part, terms) = apart(&num, &den);
    assert!(poly_part.is_some());
    assert!(terms.is_empty());
}

#[test]
fn apart_zero_numerator() {
    // 0 / (x+1) = 0
    let num = rat_poly(&[0]);
    let den = rat_poly(&[1, 1]);
    let (poly_part, terms) = apart(&num, &den);
    assert!(poly_part.is_none());
    assert!(terms.is_empty());
}

#[test]
fn together_roundtrip_manual() {
    // Build terms: 1/2/(x-1) + (-1/2)/(x+1) and combine
    let terms = vec![
        PartialFractionTerm {
            numer: rat_poly_r(&[(1, 2)]),
            denom: rat_poly(&[-1, 1]), // x - 1
            exp: 1,
        },
        PartialFractionTerm {
            numer: rat_poly_r(&[(-1, 2)]),
            denom: rat_poly(&[1, 1]), // x + 1
            exp: 1,
        },
    ];
    let (n, d) = together(None, &terms);
    assert!(!n.is_zero());
    assert_eq!(d.degree(), Some(2));
}

#[test]
fn together_with_polynomial_part() {
    // poly_part = x + 1, plus term = 2/(x-1)
    let pp = rat_poly(&[1, 1]); // x + 1
    let terms = vec![PartialFractionTerm {
        numer: rat_poly(&[2]),
        denom: rat_poly(&[-1, 1]), // x - 1
        exp: 1,
    }];
    let (n, d) = together(Some(&pp), &terms);
    // Should reconstruct (x^2 + 1) / (x - 1)
    assert!(!n.is_zero());
    assert_eq!(d.degree(), Some(1));
}

#[test]
fn apart_then_together_consistency() {
    // apart then together should give back something proportional
    let num = rat_poly(&[3, 2, 1]); // x^2 + 2x + 3
    let den = rat_poly(&[1, -3, 3, -1]); // (x-1)^3
    let (poly_part, terms) = apart(&num, &den);
    let (n, d) = together(poly_part.as_ref(), &terms);
    // n/d should be proportional to num/den
    assert!(!n.is_zero());
    assert!(!d.is_zero());
}
