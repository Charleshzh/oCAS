"""Tests for the Python numerical-integration (Vegas) bindings."""

import math

import pytest

import ocas


def test_integrate_1d_linear():
    # ∫₀¹ x dx = 1/2.
    r = ocas.integrate_1d(lambda x: x, 0.0, 1.0, n_samples=20000, iterations=8)
    assert abs(r.integral - 0.5) < 0.01
    assert r.error > 0.0


def test_integrate_1d_square():
    # ∫₀¹ x² dx = 1/3.
    r = ocas.integrate_1d(lambda x: x * x, 0.0, 1.0, n_samples=20000, iterations=8)
    assert abs(r.integral - 1.0 / 3.0) < 0.01


def test_integrate_1d_constant():
    # ∫₀¹ 7 dx = 7 with very small error.
    r = ocas.integrate_1d(lambda x: 7.0, 0.0, 1.0, n_samples=2000, iterations=4)
    assert abs(r.integral - 7.0) < 1e-6
    assert r.error < 1e-6


def test_integrate_1d_shifted_interval():
    # ∫₂⁵ x dx = (25 - 4) / 2 = 10.5.
    r = ocas.integrate_1d(lambda x: x, 2.0, 5.0, n_samples=20000, iterations=8)
    assert abs(r.integral - 10.5) < 0.1


def test_integrate_1d_gaussian_peak():
    # ∫₀¹ exp(-50 (x-0.5)²) dx ≈ sqrt(π/50) ≈ 0.2507.
    expected = math.sqrt(math.pi / 50.0)
    r = ocas.integrate_1d(
        lambda x: math.exp(-50.0 * (x - 0.5) ** 2),
        0.0,
        1.0,
        n_bins=128,
        n_samples=20000,
        iterations=12,
    )
    # Vegas should resolve the peak; allow 5% relative tolerance.
    assert abs(r.integral - expected) < 0.05 * expected


def test_integrate_result_fields_and_iteration():
    r = ocas.integrate_1d(lambda x: x, 0.0, 1.0)
    assert isinstance(r.integral, float)
    assert isinstance(r.error, float)
    # Indexing: r[0] is integral, r[1] is error.
    assert r[0] == r.integral
    assert r[1] == r.error
    assert len(r) == 2
    assert "IntegrateResult" in repr(r)


def test_integrate_1d_rejects_bad_bounds():
    with pytest.raises(ValueError):
        ocas.integrate_1d(lambda x: x, 1.0, 0.0)


def test_integrate_1d_rejects_zero_bins():
    with pytest.raises(ValueError):
        ocas.integrate_1d(lambda x: x, 0.0, 1.0, n_bins=0)


def test_integrate_1d_rejects_nonpositive_learning_rate():
    with pytest.raises(ValueError):
        ocas.integrate_1d(lambda x: x, 0.0, 1.0, learning_rate=0.0)


def test_integrate_1d_callable_exception_propagates():
    def raise_after_some_calls(x):
        if x > 0.5:
            raise RuntimeError("boom")
        return x

    with pytest.raises(ValueError):
        ocas.integrate_1d(raise_after_some_calls, 0.0, 1.0, n_samples=100, iterations=2)


def test_integrate_1d_callable_returning_non_float_propagates():
    with pytest.raises(ValueError):
        ocas.integrate_1d(lambda x: "not a float", 0.0, 1.0, n_samples=100, iterations=2)


def test_integrate_1d_seed_reproducible():
    r1 = ocas.integrate_1d(
        lambda x: x * x, 0.0, 1.0, n_samples=5000, iterations=4, seed=12345
    )
    r2 = ocas.integrate_1d(
        lambda x: x * x, 0.0, 1.0, n_samples=5000, iterations=4, seed=12345
    )
    assert r1.integral == r2.integral
    assert r1.error == r2.error


def test_vegas_class_2d():
    # ∫₀¹∫₀¹ x·y dx dy = 1/4.
    v = ocas.Vegas(2, n_samples=20000, iterations=8, seed=1)
    r = v.integrate(lambda xs: xs[0] * xs[1])
    assert abs(r.integral - 0.25) < 0.02
    assert v.iterations == 8
    latest = v.result
    assert latest.integral == r.integral
    assert "Vegas" in repr(v)


def test_vegas_class_rejects_zero_dims():
    with pytest.raises(ValueError):
        ocas.Vegas(0)
