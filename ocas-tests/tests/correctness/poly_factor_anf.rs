//! Correctness tests for factorization over algebraic number fields
//! (Trager's algorithm, 0.17.0).
//!
//! Expected factor counts and structures were cross-validated with SymPy's
//! `factor(expr, extension=...)`; each case also verifies that the product
//! of the monic factors (with multiplicities) equals the input up to its
//! leading coefficient.

use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(n, d)
}

type UP = DenseUnivariatePolynomial<AlgebraicNumberField>;

/// ℚ(√2): minimal polynomial α² − 2.
fn sqrt2_field() -> AlgebraicNumberField {
    AlgebraicNumberField::new(RationalDomain, vec![q(-2, 1), q(0, 1), q(1, 1)])
}

/// ℚ(∛2): minimal polynomial α³ − 2.
fn cbrt2_field() -> AlgebraicNumberField {
    AlgebraicNumberField::new(RationalDomain, vec![q(-2, 1), q(0, 1), q(0, 1), q(1, 1)])
}

/// ℚ(i): minimal polynomial α² + 1.
fn gaussian_field() -> AlgebraicNumberField {
    AlgebraicNumberField::new(RationalDomain, vec![q(1, 1), q(0, 1), q(1, 1)])
}

/// ℚ(ζ₅): minimal polynomial α⁴ + α³ + α² + α + 1.
fn zeta5_field() -> AlgebraicNumberField {
    AlgebraicNumberField::new(
        RationalDomain,
        vec![q(1, 1), q(1, 1), q(1, 1), q(1, 1), q(1, 1)],
    )
}

/// ℚ(⁴√3): minimal polynomial α⁴ − 3.
fn fourthroot3_field() -> AlgebraicNumberField {
    AlgebraicNumberField::new(
        RationalDomain,
        vec![q(-3, 1), q(0, 1), q(0, 1), q(0, 1), q(1, 1)],
    )
}

/// Build an ANF polynomial from per-coefficient α-polynomials.
fn anf_poly(field: &AlgebraicNumberField, coeffs: Vec<Vec<Rational>>) -> UP {
    UP::from_coeffs(
        field.clone(),
        coeffs.into_iter().map(|c| field.element(c)).collect(),
    )
}

/// Build an ANF polynomial with rational (constant) coefficients.
fn rational_poly(field: &AlgebraicNumberField, coeffs: &[i64]) -> UP {
    anf_poly(field, coeffs.iter().map(|&c| vec![q(c, 1)]).collect())
}

/// Linear ANF polynomial x − c where c = Σ c_i α^i.
fn x_minus(field: &AlgebraicNumberField, c: &[i64]) -> UP {
    anf_poly(
        field,
        vec![c.iter().map(|&v| q(-v, 1)).collect(), vec![q(1, 1)]],
    )
}

/// Check `f == lc(f) · ∏ factor^mult` (factors are monic).
fn reconstructs(field: &AlgebraicNumberField, f: &UP, factors: &[(UP, usize)]) -> bool {
    let mut acc = UP::from_coeffs(field.clone(), vec![field.one()]);
    for (h, e) in factors {
        for _ in 0..*e {
            acc = acc.mul(h);
        }
    }
    let lc = f.leading_coeff().cloned().expect("nonzero input");
    acc.mul_scalar(&lc) == *f
}

/// Assert the factorization has the expected (monic factor, multiplicity)
/// multiset and reconstructs the input.
fn check_factorization(field: &AlgebraicNumberField, f: &UP, expected: &[(UP, usize)]) {
    let factors = f.factor();
    assert_eq!(
        factors.len(),
        expected.len(),
        "factor count mismatch: got {factors:?}"
    );
    for (g, m) in expected {
        assert!(
            factors.iter().any(|(h, e)| h == g && e == m),
            "missing expected factor {g:?} (mult {m}); got {factors:?}"
        );
    }
    assert!(
        reconstructs(field, f, &factors),
        "factors do not reconstruct the input"
    );
}

/// Assert `f` is irreducible over the field (single factor, multiplicity 1,
/// equal to `f` up to a unit).
fn check_irreducible(field: &AlgebraicNumberField, f: &UP) {
    let factors = f.factor();
    assert_eq!(factors.len(), 1, "expected irreducible, got {factors:?}");
    assert_eq!(factors[0].1, 1);
    assert!(reconstructs(field, f, &factors));
}

// ------------------------------------------------------------------
// ℚ(√2)
// ------------------------------------------------------------------

#[test]
fn anf_sqrt2_x2_minus_2_splits() {
    // SymPy: factor(x^2 - 2, extension=sqrt(2)) == (x - sqrt(2))(x + sqrt(2))
    let field = sqrt2_field();
    let f = rational_poly(&field, &[-2, 0, 1]);
    let expected = [
        (x_minus(&field, &[0, 1]), 1),
        (x_minus(&field, &[0, -1]), 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_x2_plus_1_irreducible() {
    // SymPy: factor(x^2 + 1, extension=sqrt(2)) == x^2 + 1
    let field = sqrt2_field();
    let f = rational_poly(&field, &[1, 0, 1]);
    check_irreducible(&field, &f);
}

#[test]
fn anf_sqrt2_x2_minus_3_irreducible() {
    // SymPy: factor(x^2 - 3, extension=sqrt(2)) == x^2 - 3
    let field = sqrt2_field();
    let f = rational_poly(&field, &[-3, 0, 1]);
    check_irreducible(&field, &f);
}

#[test]
fn anf_sqrt2_x2_plus_x_plus_1_irreducible() {
    // Discriminant −3 ∉ ℚ(√2).
    let field = sqrt2_field();
    let f = rational_poly(&field, &[1, 1, 1]);
    check_irreducible(&field, &f);
}

#[test]
fn anf_sqrt2_mixed_split_and_irreducible() {
    // (x² − 2)(x² − 3) = x⁴ − 5x² + 6: two linears and one irreducible
    // quadratic over ℚ(√2).
    let field = sqrt2_field();
    let f = rational_poly(&field, &[6, 0, -5, 0, 1]);
    let x2_minus_3 = rational_poly(&field, &[-3, 0, 1]);
    let expected = [
        (x_minus(&field, &[0, 1]), 1),
        (x_minus(&field, &[0, -1]), 1),
        (x2_minus_3, 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_linear_times_quadratic() {
    // (x − α)(x² + 1) = x³ − αx² + x − α.
    let field = sqrt2_field();
    let f = anf_poly(
        &field,
        vec![
            vec![q(0, 1), q(-1, 1)],
            vec![q(1, 1)],
            vec![q(0, 1), q(-1, 1)],
            vec![q(1, 1)],
        ],
    );
    let x2_plus_1 = rational_poly(&field, &[1, 0, 1]);
    let expected = [(x_minus(&field, &[0, 1]), 1), (x2_plus_1, 1)];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_repeated_factors() {
    // (x − α)²(x + α): multiplicities 2 and 1.
    let field = sqrt2_field();
    let f = x_minus(&field, &[0, 1])
        .mul(&x_minus(&field, &[0, 1]))
        .mul(&x_minus(&field, &[0, -1]));
    let expected = [
        (x_minus(&field, &[0, 1]), 2),
        (x_minus(&field, &[0, -1]), 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_repeated_split_quadratic() {
    // (x² − 2)²: both linear factors have multiplicity 2.
    let field = sqrt2_field();
    let base = rational_poly(&field, &[-2, 0, 1]);
    let f = base.mul(&base);
    let expected = [
        (x_minus(&field, &[0, 1]), 2),
        (x_minus(&field, &[0, -1]), 2),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_nonmonic() {
    // (2x − α)(3x + α): leading coefficient 6 must be recovered.
    let field = sqrt2_field();
    let a = anf_poly(&field, vec![vec![q(0, 1), q(-1, 1)], vec![q(2, 1)]]);
    let b = anf_poly(&field, vec![vec![q(0, 1), q(1, 1)], vec![q(3, 1)]]);
    let f = a.mul(&b);
    let factors = f.factor();
    assert_eq!(factors.len(), 2);
    assert!(reconstructs(&field, &f, &factors));
}

#[test]
fn anf_sqrt2_three_linears() {
    // (x − 1)(x − α)(x − 2α).
    let field = sqrt2_field();
    let f = x_minus(&field, &[1])
        .mul(&x_minus(&field, &[0, 1]))
        .mul(&x_minus(&field, &[0, 2]));
    let expected = [
        (x_minus(&field, &[1]), 1),
        (x_minus(&field, &[0, 1]), 1),
        (x_minus(&field, &[0, 2]), 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_sqrt2_perfect_fourth_power() {
    // (x + 1)⁴: a single factor with multiplicity 4.
    let field = sqrt2_field();
    let base = rational_poly(&field, &[1, 1]);
    let f = base.mul(&base).mul(&base).mul(&base);
    let expected = [(rational_poly(&field, &[1, 1]), 4)];
    check_factorization(&field, &f, &expected);
}

// ------------------------------------------------------------------
// ℚ(∛2)
// ------------------------------------------------------------------

#[test]
fn anf_cbrt2_x3_minus_2() {
    // SymPy: factor(x^3 - 2, extension=2**(1/3))
    //   == (x - a)(x^2 + a*x + a^2), a = 2**(1/3)
    let field = cbrt2_field();
    let f = rational_poly(&field, &[-2, 0, 0, 1]);
    let quadratic = anf_poly(
        &field,
        vec![
            vec![q(0, 1), q(0, 1), q(1, 1)],
            vec![q(0, 1), q(1, 1)],
            vec![q(1, 1)],
        ],
    );
    let expected = [(x_minus(&field, &[0, 1, 0]), 1), (quadratic, 1)];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_cbrt2_x3_plus_1() {
    // SymPy: factor(x^3 + 1, extension=2**(1/3)) == (x + 1)(x^2 - x + 1)
    let field = cbrt2_field();
    let f = rational_poly(&field, &[1, 0, 0, 1]);
    let quadratic = rational_poly(&field, &[1, -1, 1]);
    let expected = [(rational_poly(&field, &[1, 1]), 1), (quadratic, 1)];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_cbrt2_alpha_and_alpha_squared_roots() {
    // (x − α)(x − α²) = x² − (α + α²)x + 2.
    let field = cbrt2_field();
    let f = anf_poly(
        &field,
        vec![
            vec![q(2, 1)],
            vec![q(0, 1), q(-1, 1), q(-1, 1)],
            vec![q(1, 1)],
        ],
    );
    let expected = [
        (x_minus(&field, &[0, 1, 0]), 1),
        (x_minus(&field, &[0, 0, 1]), 1),
    ];
    check_factorization(&field, &f, &expected);
}

// ------------------------------------------------------------------
// ℚ(i)
// ------------------------------------------------------------------

#[test]
fn anf_gaussian_x4_minus_1() {
    // SymPy: factor(x^4 - 1, extension=I) == (x-1)(x+1)(x-I)(x+I)
    let field = gaussian_field();
    let f = rational_poly(&field, &[-1, 0, 0, 0, 1]);
    let expected = [
        (x_minus(&field, &[1]), 1),
        (x_minus(&field, &[-1]), 1),
        (x_minus(&field, &[0, 1]), 1),
        (x_minus(&field, &[0, -1]), 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_gaussian_x2_plus_2_irreducible() {
    // ±i√2 ∉ ℚ(i).
    let field = gaussian_field();
    let f = rational_poly(&field, &[2, 0, 1]);
    check_irreducible(&field, &f);
}

// ------------------------------------------------------------------
// ℚ(ζ₅)
// ------------------------------------------------------------------

#[test]
fn anf_zeta5_cyclotomic_splits() {
    // SymPy: factor(x^4 + x^3 + x^2 + x + 1, extension=exp(2*pi*I/5))
    // splits into four linear factors (x − ζ₅^k), k = 1..4.
    let field = zeta5_field();
    let f = rational_poly(&field, &[1, 1, 1, 1, 1]);
    // α⁴ = −1 − α − α² − α³, so (x − α⁴) = x + 1 + α + α² + α³.
    let x_minus_a4 = anf_poly(
        &field,
        vec![vec![q(1, 1), q(1, 1), q(1, 1), q(1, 1)], vec![q(1, 1)]],
    );
    let expected = [
        (x_minus(&field, &[0, 1, 0, 0]), 1),
        (x_minus(&field, &[0, 0, 1, 0]), 1),
        (x_minus(&field, &[0, 0, 0, 1]), 1),
        (x_minus_a4, 1),
    ];
    check_factorization(&field, &f, &expected);
}

#[test]
fn anf_zeta5_x5_minus_1_five_linears() {
    // x⁵ − 1 = (x − 1)(x⁴ + x³ + x² + x + 1) → five linear factors.
    let field = zeta5_field();
    let f = rational_poly(&field, &[-1, 0, 0, 0, 0, 1]);
    let factors = f.factor();
    assert_eq!(factors.len(), 5);
    assert!(
        factors
            .iter()
            .all(|(g, m)| *m == 1 && g.degree() == Some(1))
    );
    assert!(reconstructs(&field, &f, &factors));
}

// ------------------------------------------------------------------
// ℚ(⁴√3)
// ------------------------------------------------------------------

#[test]
fn anf_fourthroot3_symbolica_quartic() {
    // Symbolica `algebraic_extension` test:
    // z⁴ + z³ + (2 + a − a²)z² + (1 + a² − 2a³)z − 2
    //   = (z² + (1 − a)z + (1 − a²))(z² + az + (1 + a²)).
    let field = fourthroot3_field();
    let f = anf_poly(
        &field,
        vec![
            vec![q(-2, 1)],
            vec![q(1, 1), q(0, 1), q(1, 1), q(-2, 1)],
            vec![q(2, 1), q(1, 1), q(-1, 1)],
            vec![q(1, 1)],
            vec![q(1, 1)],
        ],
    );
    let f1 = anf_poly(
        &field,
        vec![
            vec![q(1, 1), q(0, 1), q(-1, 1)],
            vec![q(1, 1), q(-1, 1)],
            vec![q(1, 1)],
        ],
    );
    let f2 = anf_poly(
        &field,
        vec![
            vec![q(1, 1), q(0, 1), q(1, 1)],
            vec![q(0, 1), q(1, 1)],
            vec![q(1, 1)],
        ],
    );
    let expected = [(f1, 1), (f2, 1)];
    check_factorization(&field, &f, &expected);
}

// ------------------------------------------------------------------
// Larger degree (performance-target scale: degree ≤ 12)
// ------------------------------------------------------------------

#[test]
fn anf_sqrt2_degree_12_product() {
    // (x² − 2)²(x² + 1)²(x − 1)(x + 1)(x² − x + 1) — degree 12 with
    // repeated factors; the x² − 2 part splits over ℚ(√2).
    let field = sqrt2_field();
    let split = rational_poly(&field, &[-2, 0, 1]); // x² − 2
    let quad1 = rational_poly(&field, &[1, 0, 1]); // x² + 1
    let lin1 = rational_poly(&field, &[-1, 1]); // x − 1
    let lin2 = rational_poly(&field, &[1, 1]); // x + 1
    let quad2 = rational_poly(&field, &[1, -1, 1]); // x² − x + 1
    let f = split
        .mul(&split)
        .mul(&quad1.mul(&quad1))
        .mul(&lin1)
        .mul(&lin2)
        .mul(&quad2);
    assert_eq!(f.degree(), Some(12));
    let factors = f.factor();
    assert!(reconstructs(&field, &f, &factors));
    let mult_of = |g: &UP| factors.iter().find(|(h, _)| h == g).map(|(_, m)| *m);
    assert_eq!(mult_of(&x_minus(&field, &[0, 1])), Some(2));
    assert_eq!(mult_of(&x_minus(&field, &[0, -1])), Some(2));
    assert_eq!(mult_of(&quad1), Some(2));
    assert_eq!(mult_of(&lin1), Some(1));
    assert_eq!(mult_of(&lin2), Some(1));
    assert_eq!(mult_of(&quad2), Some(1));
}

#[test]
fn anf_cbrt2_degree_9_product() {
    // (x³ − 2)(x³ + 1)(x³ − 1) over ℚ(∛2): x³ − 2 → linear + quadratic,
    // x³ ± 1 → linear + irreducible quadratic each.
    let field = cbrt2_field();
    let a = rational_poly(&field, &[-2, 0, 0, 1]);
    let b = rational_poly(&field, &[1, 0, 0, 1]);
    let c = rational_poly(&field, &[-1, 0, 0, 1]);
    let f = a.mul(&b).mul(&c);
    assert_eq!(f.degree(), Some(9));
    let factors = f.factor();
    assert!(reconstructs(&field, &f, &factors));
    // Three linear factors: x − α, x + 1, x − 1.
    let linears = factors
        .iter()
        .filter(|(g, _)| g.degree() == Some(1))
        .count();
    assert_eq!(linears, 3, "expected three linear factors, got {factors:?}");
}
