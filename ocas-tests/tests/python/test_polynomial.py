"""Tests for the Python `Polynomial` class."""

import pytest

import ocas


def test_polynomial_construct_integer():
    p = ocas.Polynomial([1, 2, 1])  # 1 + 2x + x^2
    assert p.degree() == 2
    assert p.len() == 3


def test_polynomial_coeffs():
    p = ocas.Polynomial([1, 2, 1])
    coeffs = p.coeffs()
    assert coeffs == ["1", "2", "1"]


def test_polynomial_eval_integer():
    p = ocas.Polynomial([1, 2, 1])  # x^2 + 2x + 1
    assert p.eval(2) == "9"  # 4 + 4 + 1


def test_polynomial_add():
    a = ocas.Polynomial([1, 1])  # 1 + x
    b = ocas.Polynomial([1, 1])
    c = a + b
    assert c.coeffs() == ["2", "2"]


def test_polynomial_sub():
    a = ocas.Polynomial([3, 3])
    b = ocas.Polynomial([1, 1])
    c = a - b
    assert c.coeffs() == ["2", "2"]


def test_polynomial_mul():
    a = ocas.Polynomial([1, 1])  # (1 + x)
    b = ocas.Polynomial([1, -1])  # (1 - x)
    c = a * b  # 1 - x^2
    assert c.coeffs() == ["1", "0", "-1"]


def test_polynomial_neg():
    p = ocas.Polynomial([1, -2, 3])
    n = -p
    assert n.coeffs() == ["-1", "2", "-3"]


def test_polynomial_eq():
    a = ocas.Polynomial([1, 2, 3])
    b = ocas.Polynomial([1, 2, 3])
    assert a == b
    c = ocas.Polynomial([1, 2])
    assert a != c


def test_polynomial_derivative():
    p = ocas.Polynomial([1, 2, 3])  # 1 + 2x + 3x^2
    d = p.derivative()  # 2 + 6x
    assert d.coeffs() == ["2", "6"]


def test_polynomial_integral_rational():
    p = ocas.Polynomial([1], domain="rational")  # constant 1
    integ = p.integral()  # x
    assert integ.coeffs() == ["0", "1"]


def test_polynomial_gcd():
    # gcd(x^2-1, x-1) = x-1  (up to unit)
    a = ocas.Polynomial([1, 0, 1])  # 1 + x^2 (use positive to keep simple)
    b = ocas.Polynomial([1, 0])  # 1 + x... let's use real gcd
    # gcd of (x^2-1) and (x+1) is (x+1)
    a2 = ocas.Polynomial([-1, 0, 1])  # x^2 - 1
    b2 = ocas.Polynomial([1, 1])  # x + 1
    g = a2.gcd(b2)
    assert g.degree() == 1


def test_polynomial_div_rem():
    a = ocas.Polynomial([1, 0, 0, 1])  # x^3 + 1
    b = ocas.Polynomial([1, 1])  # x + 1
    result = a.div_rem(b)
    assert result is not None
    q, r = result
    # x^3+1 = (x+1)(x^2-x+1) + 0
    assert r.is_zero()


def test_polynomial_square_free_factorization():
    # (x-1)^2 = x^2 - 2x + 1
    p = ocas.Polynomial([1, -2, 1])
    factors = p.square_free_factorization()
    assert len(factors) == 1
    assert factors[0].multiplicity == 2
    assert factors[0].factor.degree() == 1


def test_polynomial_is_square_free():
    p = ocas.Polynomial([1, 2, 1])  # not square-free: (x+1)^2
    assert not p.is_square_free()
    q = ocas.Polynomial([1, 1])  # square-free
    assert q.is_square_free()


def test_polynomial_finite_field():
    p = ocas.Polynomial([1, 2, 1], domain=ocas.FiniteField(5))
    assert p.degree() == 2
    # eval at 2 over GF(5): 1 + 2*2 + 1*4 = 9 = 4 (mod 5)
    assert p.eval(2) == "4"


def test_polynomial_domain_mismatch_add():
    a = ocas.Polynomial([1, 2])  # integer
    b = ocas.Polynomial([1, 2], domain="rational")
    with pytest.raises(TypeError):
        a + b


def test_polynomial_primitive_part():
    p = ocas.Polynomial([2, 4, 6])  # 2(1 + 2x + 3x^2)
    pp = p.primitive_part()
    assert pp.coeffs() == ["1", "2", "3"]


def test_polynomial_zero():
    p = ocas.Polynomial([])
    assert p.is_zero()
    assert p.degree() is None
