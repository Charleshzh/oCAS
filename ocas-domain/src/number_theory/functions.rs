//! Multiplicative number-theory functions built on [`factor_integer`]:
//! Euler's φ, the Möbius function μ, the divisor functions τ and σ_k, and
//! Liouville's λ.

use super::factor::factor_integer;
use crate::Integer;

/// Euler's totient `φ(n)`: the number of residues in `[1, |n|]` coprime to
/// `n`. Computed as `|n|·∏(1 − 1/p)` over the distinct prime factors.
/// `φ(0) = 0` and `φ(±1) = 1` by convention.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::functions::euler_phi;
///
/// assert_eq!(euler_phi(&Integer::from(9)), Integer::from(6));
/// assert_eq!(euler_phi(&Integer::from(36)), Integer::from(12));
/// assert_eq!(euler_phi(&Integer::from(97)), Integer::from(96)); // prime
/// ```
pub fn euler_phi(n: &Integer) -> Integer {
    let m = n.abs();
    if m <= Integer::from(1) {
        return m;
    }
    let mut phi = m.clone();
    for (p, _) in factor_integer(&m) {
        phi = (&phi / &p) * &(&p - &Integer::from(1));
    }
    phi
}

/// The Möbius function `μ(n)`: `1` for `n = ±1`, `0` when `n` has a squared
/// prime factor, otherwise `(−1)^k` for `k` distinct prime factors.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::functions::moebius_mu;
///
/// assert_eq!(moebius_mu(&Integer::from(1)), 1);
/// assert_eq!(moebius_mu(&Integer::from(6)), 1);   // 2 primes
/// assert_eq!(moebius_mu(&Integer::from(30)), -1); // 3 primes
/// assert_eq!(moebius_mu(&Integer::from(12)), 0);  // 2² divides 12
/// ```
pub fn moebius_mu(n: &Integer) -> i8 {
    if n.is_zero() {
        return 0;
    }
    let m = n.abs();
    if m.is_one() {
        return 1;
    }
    let factors = factor_integer(&m);
    if factors.iter().any(|(_, e)| *e >= 2) {
        return 0;
    }
    if factors.len().is_multiple_of(2) {
        1
    } else {
        -1
    }
}

/// The divisor function `τ(n)` (also `d(n)`): the number of positive
/// divisors of `|n|`, computed as `∏(e_i + 1)`. `τ(0) = 0` by convention.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::functions::divisor_tau;
///
/// assert_eq!(divisor_tau(&Integer::from(12)), Integer::from(6));
/// assert_eq!(divisor_tau(&Integer::from(97)), Integer::from(2)); // prime
/// ```
pub fn divisor_tau(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::from(0);
    }
    let mut tau = Integer::from(1);
    for (_, e) in factor_integer(&n.abs()) {
        tau *= &Integer::from((e + 1) as i64);
    }
    tau
}

/// The divisor function `σ_k(n)`: the sum of the `k`-th powers of the
/// positive divisors of `|n|`, computed as `∏ (p^{k(e+1)} − 1)/(p^k − 1)`.
/// `σ_0(n) = τ(n)`; `σ_k(0) = 0` by convention.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::functions::divisor_sigma;
///
/// assert_eq!(divisor_sigma(&Integer::from(12), 1), Integer::from(28));
/// assert_eq!(divisor_sigma(&Integer::from(12), 2), Integer::from(210));
/// assert_eq!(divisor_sigma(&Integer::from(12), 0), Integer::from(6));
/// ```
pub fn divisor_sigma(n: &Integer, k: u32) -> Integer {
    if n.is_zero() {
        return Integer::from(0);
    }
    if k == 0 {
        return divisor_tau(n);
    }
    let one = Integer::from(1);
    let mut sigma = one.clone();
    for (p, e) in factor_integer(&n.abs()) {
        let pk = p.pow_u32(k);
        let num = pk.pow_u32(e + 1) - &one;
        let den = &pk - &one;
        sigma *= &(num / den);
    }
    sigma
}

/// Liouville's function `λ(n) = (−1)^Ω(n)`, where `Ω(n)` counts prime
/// factors with multiplicity. `λ(0) = 0` and `λ(±1) = 1` by convention.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::functions::liouville_lambda;
///
/// assert_eq!(liouville_lambda(&Integer::from(12)), -1); // 12 = 2²·3, Ω = 3
/// assert_eq!(liouville_lambda(&Integer::from(6)), 1);   // Ω = 2
/// ```
pub fn liouville_lambda(n: &Integer) -> i8 {
    if n.is_zero() {
        return 0;
    }
    let omega: u32 = factor_integer(&n.abs()).iter().map(|&(_, e)| e).sum();
    if omega.is_multiple_of(2) { 1 } else { -1 }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn phi_values() {
        // φ(1..12) = 1,1,2,2,4,2,6,4,6,4,10,4.
        let expected = [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4];
        for (i, &want) in expected.iter().enumerate() {
            assert_eq!(euler_phi(&b(i as i64 + 1)), b(want), "φ({})", i + 1);
        }
        // Multiplicativity: φ(1000003·1000033) = (p−1)(q−1).
        let p = b(1_000_003);
        let q = b(1_000_033);
        assert_eq!(euler_phi(&(&p * &q)), (&p - &b(1)) * (&q - &b(1)));
    }

    #[test]
    fn phi_identity_sum_over_divisors() {
        // ∑_{d|n} φ(d) = n for n = 1..60.
        for n in 1..60i64 {
            let nb = b(n);
            let mut sum = b(0);
            for d in 1..=n {
                if n % d == 0 {
                    sum += &euler_phi(&b(d));
                }
            }
            assert_eq!(sum, nb, "∑φ(d) ≠ {n}");
        }
    }

    #[test]
    fn mu_values() {
        assert_eq!(moebius_mu(&b(1)), 1);
        assert_eq!(moebius_mu(&b(2)), -1);
        assert_eq!(moebius_mu(&b(6)), 1);
        assert_eq!(moebius_mu(&b(30)), -1);
        assert_eq!(moebius_mu(&b(210)), 1);
        assert_eq!(moebius_mu(&b(4)), 0);
        assert_eq!(moebius_mu(&b(12)), 0);
        assert_eq!(moebius_mu(&b(18)), 0);
        assert_eq!(moebius_mu(&b(0)), 0);
        assert_eq!(moebius_mu(&b(-6)), 1); // sign ignored
    }

    #[test]
    fn tau_values() {
        // τ(1..12) = 1,2,2,3,2,4,2,4,3,4,2,6.
        let expected = [1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6];
        for (i, &want) in expected.iter().enumerate() {
            assert_eq!(divisor_tau(&b(i as i64 + 1)), b(want), "τ({})", i + 1);
        }
        // Prime power: τ(p^k) = k+1.
        assert_eq!(divisor_tau(&b(3).pow_u32(7)), b(8));
    }

    #[test]
    fn sigma_values() {
        // σ(6) = 1+2+3+6 = 12; σ(12) = 28; σ_2(4) = 1+4+16 = 21.
        assert_eq!(divisor_sigma(&b(6), 1), b(12));
        assert_eq!(divisor_sigma(&b(12), 1), b(28));
        assert_eq!(divisor_sigma(&b(4), 2), b(21));
        assert_eq!(divisor_sigma(&b(12), 0), b(6));
        // Perfect numbers: σ(n) = 2n.
        for n in [6i64, 28, 496, 8128] {
            assert_eq!(divisor_sigma(&b(n), 1), b(2 * n));
        }
    }

    #[test]
    fn lambda_values() {
        assert_eq!(liouville_lambda(&b(1)), 1);
        assert_eq!(liouville_lambda(&b(2)), -1);
        assert_eq!(liouville_lambda(&b(4)), 1); // Ω = 2
        assert_eq!(liouville_lambda(&b(12)), -1); // Ω = 3
        assert_eq!(liouville_lambda(&b(6)), 1);
        assert_eq!(liouville_lambda(&b(30)), -1); // Ω = 3
        assert_eq!(liouville_lambda(&b(0)), 0);
    }
}
