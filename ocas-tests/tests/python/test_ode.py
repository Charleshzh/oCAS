"""Tests for the ODE solver bindings."""

import ocas
import pytest


def eq(s):
    return ocas.Expression(s)


def test_classify_first_order_linear():
    e = eq("Derivative(y(x), x) - y(x)")
    types = ocas.classify_ode(e, "y", "x")
    assert "LinearFirst" in types


def test_classify_second_order_constant_coeff():
    e = eq("Derivative(y(x), x, x) - y(x)")
    types = ocas.classify_ode(e, "y", "x")
    assert "LinearConstantCoeff" in types


def test_dsolve_first_order_linear():
    e = eq("Derivative(y(x), x) - y(x)")
    sol = ocas.dsolve(e, "y", "x")
    assert "exp" in sol
    assert "C1" in sol


def test_dsolve_separable():
    e = eq("Derivative(y(x), x) - x*y(x)")
    sol = ocas.dsolve(e, "y", "x")
    assert "unsolved" not in sol.lower()


def test_dsolve_hint_linear_first():
    e = eq("Derivative(y(x), x) - y(x)")
    sol = ocas.dsolve(e, "y", "x", hint="LinearFirst")
    assert "exp" in sol


def test_dsolve_ivp_first_order():
    e = eq("Derivative(y(x), x) - y(x)")
    sol = ocas.dsolve_ivp(e, "y", "x", "2")
    assert "exp" in sol
    assert "C1" not in sol


def test_dsolve_ivp_second_order_trig():
    e = eq("Derivative(y(x), x, x) + y(x)")
    sol = ocas.dsolve_ivp(e, "y", "x", "0", "1")
    assert "sin" in sol
    assert "C1" not in sol


def test_dsolve_ivp_distinct_roots():
    e = eq("Derivative(y(x), x, x) - 3*Derivative(y(x), x) + 2*y(x)")
    sol = ocas.dsolve_ivp(e, "y", "x", "1", "0")
    assert "exp" in sol
    assert "C1" not in sol
