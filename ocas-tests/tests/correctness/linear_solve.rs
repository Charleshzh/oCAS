use ocas::prelude::*;

#[test]
fn linear_solve_simple_rational_2x2() {
    let a = vec![vec![2, 1], vec![1, -1]];
    let b = vec![5, 1];
    let solution = solve_linear_rational(&a, &b).unwrap();
    assert_eq!(solution, vec![(2, 1), (1, 1)]); // x=2, y=1
}

#[test]
fn linear_solve_simple_integer_2x2() {
    let a = vec![vec![1, 1], vec![1, -1]];
    let b = vec![3, 1];
    let solution = solve_linear_integer(&a, &b).unwrap();
    assert_eq!(solution, vec![2, 1]); // x=2, y=1
}

#[test]
fn linear_solve_medium_rational_3x3() {
    let a = vec![vec![1, 1, 1], vec![1, -1, 0], vec![0, 1, -1]];
    let b = vec![6, 1, 1];
    let solution = solve_linear_rational(&a, &b).unwrap();
    let x: Vec<f64> = solution
        .iter()
        .map(|(n, d)| *n as f64 / *d as f64)
        .collect();
    assert_eq!(x, vec![3.0, 2.0, 1.0]);
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn linear_solve_complex_diophantine() {
    let sol = solve_diophantine(3, 5, 1).unwrap();
    assert_eq!(sol.particular, (2, -1));
    assert_eq!(sol.general, (5, -3));
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn linear_solve_very_complex_diophantine_verify() {
    let sol = solve_diophantine(17, 31, 100).unwrap();
    let (x0, y0) = sol.particular;
    let (tx, ty) = sol.general;
    // Verify a particular solution.
    assert_eq!(17 * x0 + 31 * y0, 100);
    // Verify homogeneous part.
    assert_eq!(17 * tx + 31 * ty, 0);
}
