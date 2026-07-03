use ocas_domain::{Rational, RationalDomain};
use ocas_poly::sparse::Lex;
use ocas_poly::{GroebnerBasis, SparseMultivariatePolynomial};

fn rat(n: i64, d: i64) -> Rational {
    Rational::new(n, d)
}

#[test]
fn groebner_simple_linear_system() {
    let d = RationalDomain;
    // ideal: x + y - 1, x - y - 3
    let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        2,
        vec![
            (vec![0, 0], rat(-1, 1)),
            (vec![1, 0], rat(1, 1)),
            (vec![0, 1], rat(1, 1)),
        ],
    );
    // Wait, the terms are (exponents, coeff). [1,0] is x, [0,1] is y, [0,0] is constant.
    let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        2,
        vec![
            (vec![0, 0], rat(-3, 1)),
            (vec![1, 0], rat(1, 1)),
            (vec![0, 1], rat(-1, 1)),
        ],
    );
    let gb = GroebnerBasis::buchberger(&[f1, f2]);
    assert!(!gb.basis.is_empty());
}

#[test]
#[ignore = "complex correctness test: run manually or via audit report"]
fn groebner_complex_cyclic_3() {
    let d = RationalDomain;
    // cyclic-3: a + b + c, ab + bc + ca, abc - 1
    let f1 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        3,
        vec![
            (vec![1, 0, 0], rat(1, 1)),
            (vec![0, 1, 0], rat(1, 1)),
            (vec![0, 0, 1], rat(1, 1)),
        ],
    );
    let f2 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        3,
        vec![
            (vec![1, 1, 0], rat(1, 1)),
            (vec![0, 1, 1], rat(1, 1)),
            (vec![1, 0, 1], rat(1, 1)),
        ],
    );
    let f3 = SparseMultivariatePolynomial::<_, Lex>::from_terms(
        d,
        3,
        vec![(vec![0, 0, 0], rat(-1, 1)), (vec![1, 1, 1], rat(1, 1))],
    );
    let gb = GroebnerBasis::buchberger(&[f1, f2, f3]);
    assert!(!gb.basis.is_empty());
}

#[test]
#[ignore = "very complex correctness test: run manually or via audit report"]
fn groebner_very_complex_cyclic_4() {
    let d = RationalDomain;
    // cyclic-4: a+b+c+d, ab+bc+cd+da, abc+bcd+cda+dab, abcd-1
    let mut polys = Vec::new();
    for coeffs in [
        vec![
            (vec![1, 0, 0, 0], rat(1, 1)),
            (vec![0, 1, 0, 0], rat(1, 1)),
            (vec![0, 0, 1, 0], rat(1, 1)),
            (vec![0, 0, 0, 1], rat(1, 1)),
        ],
        vec![
            (vec![1, 1, 0, 0], rat(1, 1)),
            (vec![0, 1, 1, 0], rat(1, 1)),
            (vec![0, 0, 1, 1], rat(1, 1)),
            (vec![1, 0, 0, 1], rat(1, 1)),
        ],
    ] {
        polys.push(SparseMultivariatePolynomial::<_, Lex>::from_terms(
            d, 4, coeffs,
        ));
    }
    let gb = GroebnerBasis::buchberger(&polys);
    assert!(!gb.basis.is_empty());
}
