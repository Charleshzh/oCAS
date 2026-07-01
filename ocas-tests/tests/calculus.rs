//! End-to-end calculus tests for oCAS.
//!
//! These tests exercise the full pipeline: parsing an expression string,
//! normalizing it, applying calculus operations, and simplifying the result.

use ocas::prelude::*;
use ocas_atom::normalize::normalize;
use ocas_core::arena::Arena;

fn diff_string(input: &str, var: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).expect("parse should succeed");
    let result = diff(&ctx, atom, Symbol::new(var));
    normalize(&ctx, result).to_string()
}

fn integrate_string(input: &str, var: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).expect("parse should succeed");
    let result = integrate(&ctx, atom, Symbol::new(var));
    normalize(&ctx, result).to_string()
}

fn taylor_string(input: &str, var: &str, point: i64, order: usize) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).expect("parse should succeed");
    let result = taylor(&ctx, atom, Symbol::new(var), ctx.num(point), order);
    normalize(&ctx, result).to_string()
}

#[test]
fn differentiate_polynomial() {
    assert_eq!(diff_string("x^3", "x"), "3*(x^2)");
}

#[test]
fn differentiate_trigonometric() {
    let s = diff_string("sin(x)", "x");
    assert_eq!(s, "cos(x)");
}

#[test]
fn differentiate_product_rule() {
    let s = diff_string("x*sin(x)", "x");
    assert!(s.contains("sin(x)"));
    assert!(s.contains("x*(cos(x))"));
}

#[test]
fn differentiate_chain_rule() {
    let s = diff_string("exp(x^2)", "x");
    assert!(s.contains("exp(x^2)"));
    assert!(s.contains("2*x"));
}

#[test]
fn integrate_power() {
    let s = integrate_string("x^2", "x");
    assert!(s.contains("x^3"));
}

#[test]
fn integrate_inverse() {
    assert_eq!(integrate_string("1/x", "x"), "log(x)");
}

#[test]
fn integrate_trigonometric() {
    assert_eq!(integrate_string("cos(x)", "x"), "sin(x)");
}

#[test]
fn taylor_exponential() {
    let s = taylor_string("exp(x)", "x", 0, 3);
    assert!(s.contains("1"));
    assert!(s.contains("x"));
    assert!(s.contains("x^2"));
    assert!(s.contains("x^3"));
}

#[test]
fn taylor_sine() {
    let s = taylor_string("sin(x)", "x", 0, 5);
    assert!(s.contains("x"));
    assert!(s.contains("x^3"));
    assert!(s.contains("x^5"));
}

#[test]
fn unknown_function_derivative() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, "f(x)").unwrap();
    let result = diff(&ctx, atom, Symbol::new("x"));
    assert_eq!(result.to_string(), "Derivative(f(x), x)");
}

#[test]
fn unknown_integral() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, "1/f(x)").unwrap();
    let result = integrate(&ctx, atom, Symbol::new("x"));
    assert!(result.to_string().starts_with("Integral("));
    assert!(result.to_string().contains("f(x)"));
}
