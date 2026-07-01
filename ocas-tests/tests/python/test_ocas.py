"""End-to-end tests for the oCAS Python bindings.

Run with: `maturin develop` in ocas-py/, then `pytest ocas-tests/tests/python/`.
"""

import gc

import pytest

import ocas


# ---------------------------------------------------------------------------
# Version / import
# ---------------------------------------------------------------------------


def test_version_is_string():
    assert isinstance(ocas.__version__, str)
    assert ocas.__version__ != ""


def test_expression_equality():
    e1 = ocas.Expression("x^2 + 1")
    e2 = ocas.Expression("x^2 + 1")
    assert e1 == e2


def test_expression_inequality():
    e1 = ocas.Expression("x^2")
    e2 = ocas.Expression("x^3")
    assert e1 != e2


def test_expression_hash_consistent():
    e1 = ocas.Expression("x^2 + 1")
    e2 = ocas.Expression("x^2 + 1")
    assert hash(e1) == hash(e2)


def test_expression_hash_in_set():
    e1 = ocas.Expression("x")
    e2 = ocas.Expression("x")
    s = {e1, e2}
    assert len(s) == 1  # same value → same hash → deduplicated


# ---------------------------------------------------------------------------
# Expression parsing and string round-trip
# ---------------------------------------------------------------------------


def test_expression_str():
    e = ocas.Expression("x^2 + 2*x + 1")
    assert "x" in str(e)


def test_expression_repr():
    e = ocas.Expression("x")
    assert "Expression" in repr(e)


def test_parse_invalid_raises():
    with pytest.raises(ValueError):
        ocas.Expression("@@@invalid@@@")


# ---------------------------------------------------------------------------
# Calculus
# ---------------------------------------------------------------------------


def test_diff_basic():
    # d/dx(x^2) = 2*x
    e = ocas.Expression("x^2")
    d = e.diff("x")
    assert str(d) == "2*x"


def test_diff_constant():
    e = ocas.Expression("5")
    d = e.diff("x")
    assert str(d) == "0"


def test_integrate_basic():
    # ∫ 2*x dx — the integrator leaves the result as 2*(2^-1)*(x^2)
    e = ocas.Expression("2*x")
    result = e.integrate("x")
    s = str(result)
    assert "(x^2)" in s
    assert "2*(2^-1)" in s


def test_taylor_exp():
    e = ocas.Expression("exp(x)")
    series = e.taylor("x", ocas.Expression("0"), 3)
    s = str(series)
    assert "1 + x" in s
    assert "(x^2)" in s
    assert "(x^3)" in s


# ---------------------------------------------------------------------------
# Simplification
# ---------------------------------------------------------------------------


def test_simplify_mul_zero():
    e = ocas.Expression("x*0")
    assert str(e.simplify()) == "0"


def test_simplify_mul_one():
    e = ocas.Expression("x*1")
    assert str(e.simplify()) == "x"


# ---------------------------------------------------------------------------
# Substitution
# ---------------------------------------------------------------------------


def test_substitute():
    e = ocas.Expression("x^2 + 1")
    y = ocas.Expression("y")
    result = e.substitute("x", y)
    assert "(y^2)" in str(result)


def test_substitute_numeric():
    e = ocas.Expression("x^2")
    two = ocas.Expression("2")
    result = e.substitute("x", two)
    # The default rule set does not evaluate numeric powers, so 2^2 is left
    # as-is; the important thing is that x was replaced.
    assert str(result) == "2^2"


# ---------------------------------------------------------------------------
# Operators
# ---------------------------------------------------------------------------


def test_add_operator():
    x = ocas.Expression("x")
    y = ocas.Expression("y")
    assert str(x + y) == "x + y"


def test_mul_operator():
    x = ocas.Expression("x")
    y = ocas.Expression("y")
    assert str(x * y) == "x*y"


def test_sub_operator():
    x = ocas.Expression("x")
    y = ocas.Expression("y")
    assert "(-1*y)" in str(x - y)


def test_pow_operator():
    x = ocas.Expression("x")
    three = ocas.Expression("3")
    assert str(x ** three) == "x^3"


def test_neg_operator():
    x = ocas.Expression("x")
    assert str(-x) == "-1*x"


def test_clone():
    e = ocas.Expression("x + 1")
    c = e.clone()
    assert str(c) == str(e)


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------


def test_solve_linear_rational():
    # 2x + 3y = 8, x - y = 1 → x = 11/5, y = 6/5
    result = ocas.solve_linear_rational([[2, 3], [1, -1]], [8, 1])
    assert result == [(11, 5), (6, 5)]


def test_solve_linear_integer_simple():
    # x + y = 3, x - y = 1 → x = 2, y = 1
    result = ocas.solve_linear_integer([[1, 1], [1, -1]], [3, 1])
    assert result == [2, 1]


def test_solve_linear_integer_no_solution():
    # 2x + 3y = 8, x - y = 1 has non-integer solution
    with pytest.raises(ValueError):
        ocas.solve_linear_integer([[2, 3], [1, -1]], [8, 1])


def test_solve_diophantine():
    # 3x + 5y = 1 → particular (2, -1), general (5k, -3k)
    result = ocas.solve_diophantine(3, 5, 1)
    assert result is not None
    assert result.particular == (2, -1)
    assert result.general == (5, -3)


def test_solve_diophantine_no_solution():
    # 2x + 4y = 1 has no solution (gcd(2,4)=2 does not divide 1)
    assert ocas.solve_diophantine(2, 4, 1) is None


# ---------------------------------------------------------------------------
# Numeric evaluation
# ---------------------------------------------------------------------------


def test_evaluator_basic():
    ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
    assert ev.evaluate([3.0, 1.0]) == [10.0]
    assert ev.evaluate([2.0, 0.0]) == [4.0]


def test_evaluator_wrong_arity():
    ev = ocas.ExpressionEvaluator("x + 1", ["x"])
    with pytest.raises(ValueError):
        ev.evaluate([1.0, 2.0])


def test_evaluator_n_params():
    ev = ocas.ExpressionEvaluator("x + y + z", ["x", "y", "z"])
    assert ev.n_params == 3


def test_evaluator_sin():
    import math

    ev = ocas.ExpressionEvaluator("sin(x)", ["x"])
    [result] = ev.evaluate([0.0])
    assert abs(result) < 1e-12
    [result] = ev.evaluate([math.pi / 2])
    assert abs(result - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# Memory pressure tests
#
# Note: tracemalloc cannot detect leaks in Rust extension code (it only
# hooks Python's allocator). These tests instead verify that creating and
# dropping many objects does not crash or leave dangling references.
# For rigorous leak detection, run the C tests under ASan/valgrind.
# ---------------------------------------------------------------------------


def test_no_crash_parse_diff_cycle():
    """Parse and differentiate many expressions; assert no crash."""
    for _ in range(200):
        e = ocas.Expression("x^2 + 2*x + 1")
        d = e.diff("x")
        s = e.simplify()
        del e, d, s
    gc.collect()


def test_no_crash_evaluator_cycle():
    """Create and use many evaluators; assert no crash."""
    for _ in range(100):
        ev = ocas.ExpressionEvaluator("x^2 + y", ["x", "y"])
        ev.evaluate([3.0, 1.0])
        del ev
    gc.collect()
