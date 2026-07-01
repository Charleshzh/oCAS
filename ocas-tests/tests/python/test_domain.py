"""Tests for the Python `Domain` classes (IntegerDomain, RationalDomain, FiniteField)."""

import pytest

import ocas


def test_integer_domain_construct():
    d = ocas.IntegerDomain()
    assert repr(d) == "IntegerDomain()"


def test_rational_domain_construct():
    d = ocas.RationalDomain()
    assert repr(d) == "RationalDomain()"


def test_finite_field_construct():
    gf = ocas.FiniteField(7)
    assert "7" in repr(gf)


def test_finite_field_rejects_small_modulus():
    with pytest.raises(ValueError):
        ocas.FiniteField(1)
    with pytest.raises(ValueError):
        ocas.FiniteField(0)


def test_finite_field_modulus_getter():
    gf = ocas.FiniteField(13)
    assert gf.modulus == "13"


def test_all_domain_classes_exposed():
    assert hasattr(ocas, "IntegerDomain")
    assert hasattr(ocas, "RationalDomain")
    assert hasattr(ocas, "FiniteField")
