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
fn evaluation_simple_constant() {
    let result = eval_f64("42", &[]);
    assert!((result - 42.0).abs() < 1e-10);
}

#[test]
fn evaluation_simple_variable() {
    let result = eval_f64("x", &[7.0]);
    assert!((result - 7.0).abs() < 1e-10);
}

#[test]
fn evaluation_medium_polynomial() {
    let result = eval_f64("x^2 + 2*x + 1", &[3.0]);
    assert!((result - 16.0).abs() < 1e-10);
}

#[test]
fn evaluation_medium_trigonometric_identity() {
    let x = 1.234;
    let result = eval_f64("sin(x)^2 + cos(x)^2", &[x]);
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn evaluation_complex_multi_variable() {
    let result = eval_f64("x*y + x^2 + y^2", &[2.0, 3.0]);
    let expected = 2.0 * 3.0 + 2.0_f64.powi(2) + 3.0_f64.powi(2);
    assert!((result - expected).abs() < 1e-10);
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn evaluation_very_complex_nested_functions() {
    let x = 1.0_f64;
    let expected = x.sin().exp() + (x * x + 1.0).ln();
    let result = eval_f64("exp(sin(x)) + log(x^2 + 1)", &[x]);
    assert!((result - expected).abs() < 1e-10);
}
