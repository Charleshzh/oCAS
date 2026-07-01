"""Tests for the Python `Matrix` class."""

import pytest

import ocas


def test_matrix_construct():
    m = ocas.Matrix([[1, 2], [3, 4]])
    assert m.nrows == 2
    assert m.ncols == 2


def test_matrix_shape():
    m = ocas.Matrix([[1, 2, 3], [4, 5, 6]])
    assert m.shape() == (2, 3)


def test_matrix_getitem():
    m = ocas.Matrix([[1, 2], [3, 4]])
    assert m[0, 0] == "1"
    assert m[1, 1] == "4"


def test_matrix_getitem_out_of_bounds():
    m = ocas.Matrix([[1, 2]])
    with pytest.raises(ValueError):
        m[5, 5]


def test_matrix_transpose():
    m = ocas.Matrix([[1, 2, 3], [4, 5, 6]])
    t = m.transpose()
    assert t.shape() == (3, 2)
    assert t[0, 1] == "4"
    assert t[2, 1] == "6"


def test_matrix_trace():
    m = ocas.Matrix([[1, 2], [3, 4]])
    assert m.trace() == "5"


def test_matrix_determinant_2x2():
    m = ocas.Matrix([[1, 2], [3, 4]])
    assert m.determinant() == "-2"


def test_matrix_determinant_singular():
    m = ocas.Matrix([[1, 2], [2, 4]])
    assert m.determinant() == "0"


def test_matrix_rank_full():
    m = ocas.Matrix([[1, 0], [0, 1]])
    assert m.rank() == 2


def test_matrix_rank_deficient():
    m = ocas.Matrix([[1, 2], [2, 4]])
    assert m.rank() == 1


def test_matrix_matmul():
    a = ocas.Matrix([[1, 2], [3, 4]])
    b = ocas.Matrix([[5, 6], [7, 8]])
    c = a @ b
    assert c[0, 0] == "19"
    assert c[1, 1] == "50"


def test_matrix_matmul_identity():
    a = ocas.Matrix([[1, 2], [3, 4]])
    ident = ocas.Matrix([[1, 0], [0, 1]])
    c = a @ ident
    assert c[0, 0] == "1"
    assert c[0, 1] == "2"


def test_matrix_add():
    a = ocas.Matrix([[1, 2], [3, 4]])
    b = ocas.Matrix([[5, 6], [7, 8]])
    c = a + b
    assert c[0, 0] == "6"
    assert c[1, 1] == "12"


def test_matrix_sub():
    a = ocas.Matrix([[5, 6], [7, 8]])
    b = ocas.Matrix([[1, 2], [3, 4]])
    c = a - b
    assert c[0, 0] == "4"
    assert c[1, 1] == "4"


def test_matrix_inverse_unimodular():
    # det = -1, integer inverse exists
    a = ocas.Matrix([[1, 2], [3, 5]])
    inv = a.inverse()
    # Verify A @ A^{-1} = I
    prod = a @ inv
    assert prod[0, 0] == "1"
    assert prod[1, 1] == "1"
    assert prod[0, 1] == "0"
    assert prod[1, 0] == "0"


def test_matrix_inverse_singular_errors():
    m = ocas.Matrix([[1, 2], [2, 4]])
    with pytest.raises(ValueError):
        m.inverse()


def test_matrix_solve():
    a = ocas.Matrix([[2, 1], [1, 1]])  # 2x+y=4, x+y=3
    sol = a.solve([4, 3])
    # x=1, y=2
    assert sol == ["1", "2"]


def test_matrix_finite_field():
    m = ocas.Matrix([[1, 2], [3, 4]], domain=ocas.FiniteField(5))
    assert m.determinant() == "3"  # 1*4 - 2*3 = -2 = 3 mod 5


def test_matrix_rows():
    m = ocas.Matrix([[1, 2], [3, 4]])
    rows = m.rows()
    assert rows == [["1", "2"], ["3", "4"]]


def test_matrix_eq():
    a = ocas.Matrix([[1, 2], [3, 4]])
    b = ocas.Matrix([[1, 2], [3, 4]])
    assert a == b
    c = ocas.Matrix([[1, 2], [3, 5]])
    assert a != c
