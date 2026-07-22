"""Tests for the Python algebraic-number-field bindings."""

import pytest

import ocas


def _sqrt2_field():
    # Minimal polynomial α² - 2 (ascending coefficients).
    return ocas.AlgebraicExtension([-2, 0, 1])


def test_algebraic_extension_degree():
    f = _sqrt2_field()
    assert f.extension_degree() == 2


def test_algebraic_extension_cbrt2():
    # α³ - 2
    f = ocas.AlgebraicExtension([-2, 0, 0, 1])
    assert f.extension_degree() == 3


def test_algebraic_extension_rejects_non_monic():
    with pytest.raises(ValueError):
        ocas.AlgebraicExtension([-2, 0, 2])  # leading coefficient 2


def test_algebraic_extension_rejects_constant():
    with pytest.raises(ValueError):
        ocas.AlgebraicExtension([5])  # degree 0


def test_algebraic_alpha_generator():
    f = _sqrt2_field()
    alpha = f.alpha()
    # α is represented as 0 + 1·α.
    assert alpha.coeffs() == ["0", "1"]


def test_algebraic_from_base_embeds_rational():
    f = _sqrt2_field()
    two = f.from_base(2)
    assert two.coeffs() == ["2"]
    half = f.from_base((1, 2))
    assert half.coeffs() == ["1/2"]


def test_algebraic_element_from_alpha_coeffs():
    f = _sqrt2_field()
    # 3 + 5·α
    e = f.element([3, 5])
    assert e.coeffs() == ["3", "5"]


def test_algebraic_polynomial_degree_and_coeffs():
    f = _sqrt2_field()
    # x² − 2: constant term first, base-domain constants.
    p = ocas.AlgebraicPolynomial(f, [-2, 0, 1])
    assert p.degree() == 2
    assert p.len() == 3
    # Each coefficient is a list of α-polynomial rational strings. A zero
    # coefficient trims to an empty list (no α-terms).
    assert p.coeffs() == [["-2"], [], ["1"]]


def test_algebraic_polynomial_with_alpha_coefficient():
    f = _sqrt2_field()
    # x² − α: the x² coefficient is 1, the constant term is 0 + (−1)·α.
    p = ocas.AlgebraicPolynomial(f, [[0, -1], 0, 1])
    assert p.degree() == 2
    assert p.coeffs()[0] == ["0", "-1"]
    assert p.coeffs()[2] == ["1"]


def test_algebraic_polynomial_factor_sqrt2_splits():
    # Over ℚ(√2): x² − 2 = (x − α)(x + α).
    f = _sqrt2_field()
    p = ocas.AlgebraicPolynomial(f, [-2, 0, 1])
    factors = p.factor()
    assert len(factors) == 2
    for fac in factors:
        assert fac.factor.degree() == 1
        assert fac.multiplicity == 1


def test_algebraic_polynomial_factor_with_alpha_irreducible():
    # x² − α is irreducible over ℚ(√2): a root would be 2^(1/4) ∉ ℚ(√2).
    f = _sqrt2_field()
    # Constant term = −α = [0, −1], x² coefficient = 1.
    p = ocas.AlgebraicPolynomial(f, [[0, -1], 0, 1])
    factors = p.factor()
    assert len(factors) == 1
    assert factors[0].factor.degree() == 2


def test_algebraic_polynomial_factor_cbrt2():
    # Over ℚ(∛2): x³ − 2 = (x − α)(x² + αx + α²).
    f = ocas.AlgebraicExtension([-2, 0, 0, 1])
    p = ocas.AlgebraicPolynomial(f, [-2, 0, 0, 1])
    factors = p.factor()
    assert len(factors) == 2
    degrees = sorted(fac.factor.degree() for fac in factors)
    assert degrees == [1, 2]


def test_algebraic_polynomial_to_string_nonempty():
    f = _sqrt2_field()
    p = ocas.AlgebraicPolynomial(f, [-2, 0, 1])
    s = str(p)
    assert "x" in s
    assert "2" in s
