"""Tests for the ocas number-theory bindings (ocas 0.21.0)."""

import pytest

import ocas


class TestFactorint:
    def test_small(self):
        assert ocas.factorint(12) == [("2", 2), ("3", 1)]

    def test_prime(self):
        assert ocas.factorint(97) == [("97", 1)]

    def test_prime_power(self):
        assert ocas.factorint(7**6) == [("7", 6)]

    def test_negative(self):
        assert ocas.factorint(-12) == [("-1", 1), ("2", 2), ("3", 1)]

    def test_semiprime(self):
        assert ocas.factorint(1000003 * 1000033) == [
            ("1000003", 1),
            ("1000033", 1),
        ]

    def test_one_and_zero(self):
        assert ocas.factorint(1) == []
        assert ocas.factorint(0) == []

    def test_string_input(self):
        assert ocas.factorint("360") == [("2", 3), ("3", 2), ("5", 1)]


class TestIsprime:
    @pytest.mark.parametrize("p", [2, 3, 5, 97, 7919, 1000003, 2**61 - 1])
    def test_primes(self, p):
        assert ocas.isprime(p)

    @pytest.mark.parametrize(
        "c",
        [0, 1, 4, 9, 561, 1105, 1729, 2047, 3277, 4033, 825265, 1000003 * 1000033],
    )
    def test_composites(self, c):
        assert not ocas.isprime(c)

    def test_large_prime(self):
        assert ocas.isprime(10**18 + 9)

    def test_negative(self):
        assert not ocas.isprime(-7)


class TestNextprime:
    def test_small(self):
        assert ocas.nextprime(10) == 11
        assert ocas.nextprime(13) == 17
        assert ocas.nextprime(0) == 2
        assert ocas.nextprime(2) == 3

    def test_large(self):
        assert ocas.nextprime(10**12) == 10**12 + 39


class TestDiscreteLog:
    def test_small_group(self):
        # 2 is a primitive root mod 11; 2^7 ≡ 7 (mod 11).
        x = ocas.discrete_log(11, 2, 7)
        assert pow(2, x, 11) == 7

    def test_smooth_prime(self):
        # p = 101, p − 1 = 2²·5².
        x = ocas.discrete_log(101, 2, pow(2, 83, 101))
        assert x == 83

    def test_composite_modulus_bsgs(self):
        # Falls back to BSGS for composite modulus.
        x = ocas.discrete_log(1009, 11, pow(11, 377, 1009))
        assert pow(11, x, 1009) == pow(11, 377, 1009)

    def test_no_solution(self):
        # 5 is not a power of 2 mod 23 (2 has order 11).
        with pytest.raises(ValueError):
            ocas.discrete_log(23, 2, 5)


class TestCrt:
    def test_sunzi(self):
        r, m = ocas.crt([3, 5, 7], [2, 3, 2])
        assert (r, m) == (23, 105)

    def test_non_coprime(self):
        r, m = ocas.crt([4, 6], [1, 3])
        assert m == 12
        assert r % 4 == 1
        assert r % 6 == 3

    def test_inconsistent(self):
        with pytest.raises(ValueError):
            ocas.crt([4, 4], [1, 2])

    def test_large(self):
        r, m = ocas.crt([1000003, 1000033, 1000037], [11, 222, 3333])
        assert m == 1000003 * 1000033 * 1000037
        assert r % 1000003 == 11
        assert r % 1000033 == 222
        assert r % 1000037 == 3333


class TestJacobiSymbol:
    def test_values(self):
        assert ocas.jacobi_symbol(2, 7) == 1
        assert ocas.jacobi_symbol(3, 7) == -1
        assert ocas.jacobi_symbol(0, 7) == 0
        assert ocas.jacobi_symbol(2, 15) == 1
        assert ocas.jacobi_symbol(7, 15) == -1

    def test_even_modulus_rejected(self):
        with pytest.raises(ValueError):
            ocas.jacobi_symbol(1, 4)


class TestNumberTheoryFunctions:
    def test_totient(self):
        assert ocas.totient(9) == 6
        assert ocas.totient(36) == 12
        assert ocas.totient(97) == 96
        assert ocas.totient(1000003 * 1000033) == 1000002 * 1000032

    def test_mobius(self):
        assert ocas.mobius(1) == 1
        assert ocas.mobius(6) == 1
        assert ocas.mobius(30) == -1
        assert ocas.mobius(12) == 0

    def test_divisor_count(self):
        assert ocas.divisor_count(12) == 6
        assert ocas.divisor_count(97) == 2
        assert ocas.divisor_count(3**7) == 8

    def test_divisor_sigma(self):
        assert ocas.divisor_sigma(12) == 28
        assert ocas.divisor_sigma(12, 2) == 210
        assert ocas.divisor_sigma(12, 0) == 6
        for n in (6, 28, 496, 8128):
            assert ocas.divisor_sigma(n) == 2 * n  # perfect numbers

    def test_liouville_lambda(self):
        assert ocas.liouville_lambda(12) == -1
        assert ocas.liouville_lambda(6) == 1
        assert ocas.liouville_lambda(1) == 1

    def test_isprime_u64(self):
        assert ocas.isprime_u64(97)
        assert not ocas.isprime_u64(561)
        assert ocas.isprime_u64(2**64 - 59)
        assert not ocas.isprime_u64(2**64 - 1)
