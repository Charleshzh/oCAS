//! Integration tests for oCAS.
//!
//! These tests exercise end-to-end workflows across multiple crates:
//! parsing an expression string, normalizing it, and converting between the
//! symbolic `Atom` representation and polynomial representations.

use ocas::prelude::*;
use ocas_atom::normalize::normalize;
use ocas_core::arena::Arena;

fn parse_normalized(input: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).expect("parse should succeed");
    let norm = normalize(&ctx, atom);
    norm.to_string()
}

#[test]
fn parse_and_normalize_number() {
    assert_eq!(parse_normalized("42"), "42");
}

#[test]
fn parse_and_normalize_simple_sum() {
    assert_eq!(parse_normalized("1 + 2"), "3");
}

#[test]
fn parse_and_normalize_symbolic_sum() {
    // The normalizer currently flattens and sorts but does not combine like
    // terms; verify the deterministic canonical shape instead.
    assert_eq!(parse_normalized("x + 2*x"), "x + (2*x)");
}

#[test]
fn parse_and_normalize_polynomial_expression() {
    let s = parse_normalized("x^2 + 2*x + 1");
    assert_eq!(s, "1 + (2*x) + (x^2)");
}

#[test]
fn dense_polynomial_from_rational_domain() {
    let domain = RationalDomain;
    let a = DenseUnivariatePolynomial::from_coeffs(
        domain,
        vec![
            Rational::new(1, 1),
            Rational::new(2, 1),
            Rational::new(1, 1),
        ],
    );
    let b = DenseUnivariatePolynomial::from_coeffs(
        domain,
        vec![Rational::new(1, 1), Rational::new(1, 1)],
    );
    let c = a.mul(&b);
    assert_eq!(c.degree(), Some(3));
    assert_eq!(c.coeff(0).unwrap().inner().numer(), &1.into());
    assert_eq!(c.coeff(1).unwrap().inner().numer(), &3.into());
    assert_eq!(c.coeff(2).unwrap().inner().numer(), &3.into());
    assert_eq!(c.coeff(3).unwrap().inner().numer(), &1.into());
}
