"""Tests for the Python tensor algebra bindings."""

import pytest

import ocas


def test_tensor_rank_and_slots():
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    assert t.name == "T"
    assert t.rank == 2
    assert t.symmetry == "none"
    assert t.slots() == [("i", "upper"), ("j", "lower")]


def test_tensor_position_aliases():
    # "up"/"down" and "contravariant"/"covariant" aliases work.
    t = ocas.Tensor("T", [("i", "up"), ("j", "down")])
    assert t.slots() == [("i", "upper"), ("j", "lower")]
    t2 = ocas.Tensor("T", [("i", "contravariant"), ("j", "covariant")])
    assert t2.slots() == [("i", "upper"), ("j", "lower")]


def test_tensor_symmetry_keywords():
    t = ocas.Tensor("S", [("a", "lower"), ("b", "lower")], symmetry="symmetric")
    assert t.symmetry == "symmetric"
    e = ocas.Tensor("eps", [("a", "lower"), ("b", "lower")], symmetry="antisymmetric")
    assert e.symmetry == "antisymmetric"
    n = ocas.Tensor("N", [("a", "lower")], symmetry="none")
    assert n.symmetry == "none"


def test_tensor_rejects_bad_position():
    with pytest.raises(ValueError):
        ocas.Tensor("T", [("i", "sideways")])


def test_tensor_rejects_bad_symmetry():
    with pytest.raises(ValueError):
        ocas.Tensor("T", [("i", "upper")], symmetry="weird")


def test_tensor_dummy_labels():
    t = ocas.Tensor("T", [("i", "upper"), ("i", "lower")])
    dummies = t.dummy_labels()
    assert len(dummies) == 1
    assert dummies[0] == "i"


def test_tensor_to_string_atom():
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    s = t.to_string_atom()
    # The atom is rendered as T(i, j) — exact format depends on the printer.
    assert "T" in s


def test_contract_partial_produces_product():
    # T^i_j · U^j_k = (TU)^i_k  (partial contraction over j).
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    u = ocas.Tensor("U", [("j", "upper"), ("k", "lower")])
    kind, payload = ocas.contract_tensors(t, u)
    assert kind == "product"
    assert isinstance(payload, list)
    assert len(payload) == 1
    result = payload[0]
    assert result.rank == 2  # surviving slots: i (upper) and k (lower)


def test_contract_no_dummy_is_plain_product():
    # T^i_j · U^k_l with no shared dummy → plain product (two factors).
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    u = ocas.Tensor("U", [("k", "upper"), ("l", "lower")])
    kind, payload = ocas.contract_tensors(t, u)
    assert kind == "product"
    assert len(payload) == 2


def test_contract_full_produces_scalar():
    # T^i_j · U^j_i — every index contracts → scalar.
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    u = ocas.Tensor("U", [("j", "upper"), ("i", "lower")])
    kind, payload = ocas.contract_tensors(t, u)
    assert kind == "scalar"
    assert isinstance(payload, str)
    assert "T" in payload
    assert "U" in payload


def test_symmetrise_sign_none_and_symmetric():
    t = ocas.Tensor("T", [("i", "upper"), ("j", "lower")])
    assert ocas.tensor_symmetrise_sign(t) == 1
    s = ocas.Tensor("S", [("a", "lower"), ("b", "lower")], symmetry="symmetric")
    assert ocas.tensor_symmetrise_sign(s) == 1


def test_symmetrise_sign_antisymmetric_is_pm1():
    e = ocas.Tensor("eps", [("a", "lower"), ("b", "lower")], symmetry="antisymmetric")
    sign = ocas.tensor_symmetrise_sign(e)
    assert sign in (1, -1)


def test_tensor_repr():
    t = ocas.Tensor("T", [("i", "upper")])
    r = repr(t)
    assert "Tensor" in r
    assert "T" in r
