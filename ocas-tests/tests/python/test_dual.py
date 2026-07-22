"""Tests for the Python hyper-dual number (forward AD) bindings."""

import pytest

import ocas


def _shape2():
    return ocas.DualShape.first_order(2)


def test_dual_shape_first_order():
    s = ocas.DualShape.first_order(3)
    assert s.n_vars == 3
    assert s.n_components == 4  # value + 3 first-order derivatives
    assert "DualShape" in repr(s)


def test_dual_shape_rejects_zero():
    with pytest.raises(ValueError):
        ocas.DualShape.first_order(0)


def test_variable_value_and_deriv():
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 3)
    assert x.value() == "3"
    assert x.deriv(0) == "1"  # ∂x/∂x = 1
    assert x.deriv(1) == "0"  # ∂x/∂y = 0
    assert x.n_vars == 2


def test_variable_rational_coefficient():
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, (1, 2))
    assert x.value() == "1/2"


def test_constant_has_zero_derivatives():
    s = _shape2()
    c = ocas.HyperDual.constant(s, 7)
    assert c.value() == "7"
    assert c.deriv(0) == "0"
    assert c.deriv(1) == "0"


def test_product_of_two_variables():
    # f(x, y) = x * y at (3, 5); ∂f/∂x = y = 5, ∂f/∂y = x = 3.
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 3)
    y = ocas.HyperDual.variable(s, 1, 5)
    f = x * y
    assert f.value() == "15"
    assert f.deriv(0) == "5"
    assert f.deriv(1) == "3"


def test_sum_and_difference():
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 2)
    y = ocas.HyperDual.variable(s, 1, 7)
    s_op = x + y
    assert s_op.value() == "9"
    assert s_op.deriv(0) == "1"
    assert s_op.deriv(1) == "1"
    d_op = x - y
    assert d_op.value() == "-5"
    assert d_op.deriv(0) == "1"
    assert d_op.deriv(1) == "-1"


def test_quotient_chain_rule():
    # f(x, y) = x / y at (3, 5); ∂f/∂x = 1/y = 1/5, ∂f/∂y = -x/y² = -3/25.
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 3)
    y = ocas.HyperDual.variable(s, 1, 5)
    f = x / y
    assert f.value() == "3/5"
    assert f.deriv(0) == "1/5"
    assert f.deriv(1) == "-3/25"


def test_negation():
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 4)
    f = -x
    assert f.value() == "-4"
    assert f.deriv(0) == "-1"


def test_division_by_zero_raises():
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 3)
    z = ocas.HyperDual.constant(s, 0)
    with pytest.raises(ValueError):
        _ = x / z


def test_shape_mismatch_arithmetic_raises():
    s2 = ocas.DualShape.first_order(2)
    s3 = ocas.DualShape.first_order(3)
    x = ocas.HyperDual.variable(s2, 0, 1)
    y = ocas.HyperDual.variable(s3, 0, 1)
    with pytest.raises(ValueError):
        _ = x + y


def test_variable_index_out_of_range():
    s = _shape2()
    with pytest.raises(ValueError):
        ocas.HyperDual.variable(s, 5, 1)


def test_chain_rule_polynomial():
    # f(x) = x^3 at point 2; use repeated multiplication.
    # f'(x) = 3x² = 12, f = 8.
    s = ocas.DualShape.first_order(1)
    x = ocas.HyperDual.variable(s, 0, 2)
    f = x * x * x
    assert f.value() == "8"
    assert f.deriv(0) == "12"


def test_mixed_partial_via_repeated_mul():
    # f(x, y) = x^2 * y at (3, 4). ∂f/∂x = 2xy = 24, ∂f/∂y = x² = 9.
    s = _shape2()
    x = ocas.HyperDual.variable(s, 0, 3)
    y = ocas.HyperDual.variable(s, 1, 4)
    f = x * x * y
    assert f.value() == "36"
    assert f.deriv(0) == "24"
    assert f.deriv(1) == "9"
