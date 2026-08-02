//! Integration tests for oCAS.
//!
//! These tests exercise end-to-end workflows across multiple crates:
//! parsing an expression string, normalizing it, and converting between the
//! symbolic `Atom` representation and polynomial representations.

use ocas::prelude::*;
use ocas_atom::Symbol;
use ocas_atom::normalize::normalize;
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
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
    assert_eq!(c.coeff(0).unwrap().numer(), Integer::from(1));
    assert_eq!(c.coeff(1).unwrap().numer(), Integer::from(3));
    assert_eq!(c.coeff(2).unwrap().numer(), Integer::from(3));
    assert_eq!(c.coeff(3).unwrap().numer(), Integer::from(1));
}

#[test]
fn canon_symmetric_tensor_via_parser() {
    // Reproduce the Python binding path: parse from string, then canonicalise.
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let mut reg = TensorRegistry::new();
    reg.register(Symbol::new("g"), SymmetrySpec::fully_symmetric(2));

    let g_ab = parse(&ctx, "g(a,b)").unwrap();
    let g_ba = parse(&ctx, "g(b,a)").unwrap();

    let ct1 = canonicalize_tensors(&ctx, g_ab, &reg).unwrap();
    let ct2 = canonicalize_tensors(&ctx, g_ba, &reg).unwrap();
    eprintln!("ct1 = {}", ct1.canonical_form);
    eprintln!("ct2 = {}", ct2.canonical_form);
    assert_eq!(
        ct1.canonical_form.to_string(),
        ct2.canonical_form.to_string(),
        "symmetric slots via parser should canonicalise consistently"
    );
}
