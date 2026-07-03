use ocas::prelude::*;
use ocas_poly::Matrix;

fn rat(n: i64, d: i64) -> Rational {
    Rational::new(n, d)
}

#[test]
fn matrix_simple_2x2_determinant() {
    let d = RationalDomain;
    let m = Matrix::from_rows(
        vec![vec![rat(1, 1), rat(2, 1)], vec![rat(3, 1), rat(4, 1)]],
        d,
    );
    let det = m.determinant().unwrap();
    assert_eq!(det, rat(-2, 1));
}

#[test]
fn matrix_simple_inverse() {
    let d = RationalDomain;
    let m = Matrix::from_rows(
        vec![vec![rat(1, 1), rat(2, 1)], vec![rat(3, 1), rat(4, 1)]],
        d,
    );
    let inv = m.inverse().unwrap();
    let product = m.matmul(&inv).unwrap();
    assert_eq!(product, Matrix::identity(2, d));
}

#[test]
fn matrix_medium_3x3_determinant() {
    let d = RationalDomain;
    let m = Matrix::from_rows(
        vec![
            vec![rat(1, 1), rat(2, 1), rat(3, 1)],
            vec![rat(4, 1), rat(5, 1), rat(6, 1)],
            vec![rat(7, 1), rat(8, 1), rat(10, 1)],
        ],
        d,
    );
    let det = m.determinant().unwrap();
    assert_eq!(det, rat(-3, 1));
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn matrix_complex_solve_linear() {
    let d = RationalDomain;
    let a = Matrix::from_rows(
        vec![vec![rat(2, 1), rat(1, 1)], vec![rat(1, 1), rat(-1, 1)]],
        d,
    );
    let b = vec![rat(3, 1), rat(0, 1)];
    let x = a.solve(&b).unwrap();
    assert_eq!(x[0], rat(1, 1));
    assert_eq!(x[1], rat(1, 1));
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn matrix_very_complex_hilbert() {
    let d = RationalDomain;
    let n = 5;
    let rows: Vec<Vec<Rational>> = (0..n)
        .map(|i| (0..n).map(|j| rat(1, (i + j + 1) as i64)).collect())
        .collect();
    let m = Matrix::from_rows(rows, d);
    let det = m.determinant().unwrap();
    assert!(!m.domain().is_zero(&det));
}
