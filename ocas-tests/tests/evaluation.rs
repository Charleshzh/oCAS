//! End-to-end evaluation tests for oCAS.
//!
//! These tests exercise the full pipeline: parsing an expression string,
//! compiling it to an evaluator, and evaluating it with numeric inputs.

use ocas::prelude::*;
use ocas_core::arena::Arena;

fn eval_f64(input: &str, params: &[f64]) -> f64 {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).expect("parse should succeed");
    let eval: ExpressionEvaluator<f64> =
        ExpressionEvaluator::compile(atom).expect("compile should succeed");
    let results = eval.evaluate(params).expect("evaluate should succeed");
    results[0]
}

#[test]
fn eval_constant() {
    let result = eval_f64("42", &[]);
    assert!((result - 42.0).abs() < 1e-10);
}

#[test]
fn eval_variable() {
    let result = eval_f64("x", &[7.0]);
    assert!((result - 7.0).abs() < 1e-10);
}

#[test]
fn eval_add() {
    let result = eval_f64("x + 3", &[2.0]);
    assert!((result - 5.0).abs() < 1e-10);
}

#[test]
fn eval_mul() {
    let result = eval_f64("x * y", &[3.0, 4.0]);
    assert!((result - 12.0).abs() < 1e-10);
}

#[test]
fn eval_polynomial() {
    // x^2 + 2*x + 1 at x=3 → 9 + 6 + 1 = 16
    let result = eval_f64("x^2 + 2*x + 1", &[3.0]);
    assert!((result - 16.0).abs() < 1e-10);
}

#[test]
fn eval_nested_arithmetic() {
    // (x + 1) * (x - 1) = x^2 - 1 at x=5 → 24
    let result = eval_f64("(x + 1) * (x - 1)", &[5.0]);
    assert!((result - 24.0).abs() < 1e-10);
}

#[test]
fn eval_sin() {
    let result = eval_f64("sin(x)", &[std::f64::consts::FRAC_PI_2]);
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
fn eval_cos() {
    let result = eval_f64("cos(x)", &[std::f64::consts::PI]);
    assert!((result + 1.0).abs() < 1e-10);
}

#[test]
fn eval_exp_log_roundtrip() {
    // log(exp(x)) = x
    let result = eval_f64("log(exp(x))", &[2.5]);
    assert!((result - 2.5).abs() < 1e-10);
}

#[test]
fn eval_sqrt() {
    let result = eval_f64("sqrt(x)", &[16.0]);
    assert!((result - 4.0).abs() < 1e-10);
}

#[test]
fn eval_pythagorean_identity() {
    // sin(x)^2 + cos(x)^2 = 1
    let x = 1.234;
    let result = eval_f64("sin(x)^2 + cos(x)^2", &[x]);
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
fn eval_complex_expression() {
    // exp(sin(x)) + log(x^2 + 1) at x=1.0
    let x = 1.0_f64;
    let expected = x.sin().exp() + (x * x + 1.0).ln();
    let result = eval_f64("exp(sin(x)) + log(x^2 + 1)", &[x]);
    assert!((result - expected).abs() < 1e-10);
}

#[test]
fn eval_power_integer() {
    let result = eval_f64("x^5", &[2.0]);
    assert!((result - 32.0).abs() < 1e-10);
}

#[test]
fn eval_power_negative_exp() {
    let result = eval_f64("x^(-1)", &[4.0]);
    assert!((result - 0.25).abs() < 1e-10);
}

#[test]
fn eval_wrong_param_count() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, "x + y").unwrap();
    let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(atom).unwrap();
    // 2 params expected, only 1 provided
    assert!(eval.evaluate(&[1.0]).is_err());
}

#[test]
fn eval_with_function_map() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, "square(x) + 1").unwrap();

    let mut map = FunctionMap::<f64>::new();
    map.register("square", 1, Box::new(|args| args[0] * args[0]));

    let eval = ExpressionEvaluator::compile_with(atom, map).unwrap();
    let result = eval.evaluate(&[3.0]).unwrap();
    assert!((result[0] - 10.0).abs() < 1e-10);
}

#[test]
fn eval_zero_params() {
    let result = eval_f64("sin(1) + cos(0)", &[]);
    assert!((result - (1.0_f64.sin() + 1.0)).abs() < 1e-10);
}

#[test]
fn eval_param_count_property() {
    // Multiple uses of the same variable should only need one param
    let result = eval_f64("x + x + x", &[5.0]);
    assert!((result - 15.0).abs() < 1e-10);
}
