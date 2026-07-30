//! BPSW probable-prime test and Lucas sequence machinery.
//!
//! The BPSW test combines a base-2 strong Miller–Rabin round with a strong
//! Lucas probable-prime test using Selfridge's parameter selection. No
//! composite is known to pass it. The Lucas chain code is shared with
//! Williams' p+1 factorization method.

use super::{is_prime, jacobi, mr_witness};
use crate::Integer;

/// Halve `x` modulo odd `n`: returns `y` in `[0, n)` with `2y ≡ x (mod n)`.
fn half_mod(x: &Integer, n: &Integer) -> Integer {
    let two = Integer::from(2);
    let x = x.mod_floor(n);
    if x.is_even() {
        x / &two
    } else {
        (&x + n) / &two
    }
}

/// Compute `(U_k, V_k, Q^k) mod n` for the Lucas sequence with parameters
/// `(P, Q)`, i.e. `U_0 = 0`, `U_1 = 1`, `V_0 = 2`, `V_1 = P` and the
/// recurrences `U_{m+1} = P·U_m − Q·U_{m−1}` (same for `V`).
///
/// Uses the binary ladder with the doubling identities
/// `U_{2m} = U_m·V_m`, `V_{2m} = V_m² − 2Q^m` and the increment identities
/// `2U_{m+1} = P·U_m + V_m`, `2V_{m+1} = D·U_m + P·V_m` where `D = P² − 4Q`.
/// Requires `k ≥ 0` and odd `n > 2` (for the halving step).
pub(crate) fn lucas_uv_mod(
    p: &Integer,
    q: &Integer,
    k: &Integer,
    n: &Integer,
) -> (Integer, Integer, Integer) {
    let two = Integer::from(2);
    if k.is_zero() {
        return (
            Integer::from(0),
            two.mod_floor(n),
            Integer::from(1).mod_floor(n),
        );
    }
    let d = (p * p) - &(q * &Integer::from(4));
    // Bits of k, least significant first.
    let mut bits = Vec::new();
    let mut kk = k.clone();
    while !kk.is_zero() {
        let (quot, rem) = kk.div_rem(&two);
        bits.push(rem.is_one());
        kk = quot;
    }
    let mut u = Integer::from(0);
    let mut v = two.mod_floor(n);
    let mut qk = Integer::from(1).mod_floor(n);
    for &bit in bits.iter().rev() {
        // Double the index.
        u = (&u * &v).mod_floor(n);
        v = (&v * &v - &(&qk * &two)).mod_floor(n);
        qk = (&qk * &qk).mod_floor(n);
        if bit {
            let u_new = half_mod(&(p * &u + &v), n);
            let v_new = half_mod(&(&d * &u + p * &v), n);
            u = u_new;
            v = v_new;
            qk = (&qk * q).mod_floor(n);
        }
    }
    (u, v, qk)
}

/// Strong Lucas probable-prime test with Selfridge's parameters: choose the
/// first `D` from `5, −7, 9, −11, …` with `jacobi(D, n) = −1`, set `P = 1`,
/// `Q = (1 − D)/4`, and write `n + 1 = d·2^r` with `d` odd. `n` passes if
/// `U_d ≡ 0` or `V_{d·2^i} ≡ 0 (mod n)` for some `0 ≤ i < r`.
fn strong_lucas_prp(n: &Integer) -> bool {
    // Selfridge parameter search.
    let mut d_abs: i64 = 5;
    let mut sign: i64 = 1;
    let q = loop {
        let d = Integer::from(sign * d_abs);
        match jacobi(&d, n) {
            -1 => break Integer::from((1 - sign * d_abs) / 4),
            // gcd(D, n) > 1: n is composite unless |D| == n itself (a composite
            // non-square n is always caught with |D| ≤ its smallest prime
            // factor < n; squares are excluded by the caller).
            0 => return Integer::from(d_abs) == *n,
            _ => {
                d_abs += 2;
                sign = -sign;
            }
        }
    };
    let p = Integer::from(1);
    // n + 1 = d·2^r with d odd.
    let one = Integer::from(1);
    let mut d = n + &one;
    let mut r = 0u64;
    while d.is_even() {
        d >>= 1;
        r += 1;
    }
    let (u, mut v, mut qk) = lucas_uv_mod(&p, &q, &d, n);
    if u.is_zero() || v.is_zero() {
        return true;
    }
    for _ in 1..r {
        v = (&v * &v - &(&qk * &Integer::from(2))).mod_floor(n);
        qk = (&qk * &qk).mod_floor(n);
        if v.is_zero() {
            return true;
        }
    }
    false
}

/// Test whether `n` is a BPSW probable prime: base-2 strong Miller–Rabin
/// followed by a strong Lucas test with Selfridge parameters.
///
/// No composite is known to pass BPSW. For `n < 3.317·10²⁴` the result is
/// provably correct (the Miller–Rabin part alone is deterministic there).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::primes::is_prime_bpsw;
///
/// assert!(is_prime_bpsw(&Integer::from(97)));
/// assert!(!is_prime_bpsw(&Integer::from(561)));    // Carmichael number
/// assert!(!is_prime_bpsw(&Integer::from(2047)));   // base-2 strong pseudoprime
/// ```
pub fn is_prime_bpsw(n: &Integer) -> bool {
    let two = Integer::from(2);
    if n < &two {
        return false;
    }
    for &p in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let pb = Integer::from(p as i64);
        if *n == pb {
            return true;
        }
        if n.mod_floor(&pb).is_zero() {
            return false;
        }
    }
    // BPSW excludes perfect squares up front (cheap isqrt check).
    let s = n.sqrt();
    if &s * &s == *n {
        return false;
    }
    // Base-2 strong Miller–Rabin round.
    let one = Integer::from(1);
    let mut d = n - &one;
    let mut r = 0u64;
    while d.is_even() {
        d >>= 1;
        r += 1;
    }
    if !mr_witness(n, &d, r, &two) {
        return false;
    }
    strong_lucas_prp(n)
}

/// Deterministic primality test for `n < 2^64`.
///
/// The fixed 12-witness Miller–Rabin set of [`is_prime`] is deterministic for
/// every `n < 3.317·10²⁴`, which covers the entire `u64` range.
///
/// # Example
///
/// ```
/// use ocas_domain::number_theory::primes::is_prime_u64;
///
/// assert!(is_prime_u64(97));
/// assert!(!is_prime_u64(561));
/// assert!(is_prime_u64(u64::MAX - 58)); // 2^64 − 59, the largest u64 prime
/// ```
pub fn is_prime_u64(n: u64) -> bool {
    is_prime(&Integer::new(n))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn lucas_fibonacci_identity() {
        // P = 1, Q = −1 gives the Fibonacci numbers U_k = F_k and the Lucas
        // numbers V_k = L_k: F_10 = 55, L_10 = 123, Q^10 = 1.
        let n = b(1009);
        let (u, v, qk) = lucas_uv_mod(&b(1), &b(-1), &b(10), &n);
        assert_eq!(u, b(55));
        assert_eq!(v, b(123));
        assert_eq!(qk, b(1));
    }

    #[test]
    fn lucas_matches_recurrence() {
        // Cross-check the ladder against the defining recurrence for P = 3,
        // Q = −2, k = 1..60, all modulo a large prime.
        let n = b(1_000_003);
        let p = b(3);
        let q = b(-2);
        let mut u_prev = b(0); // U_0
        let mut u_cur = b(1); // U_1
        let mut v_prev = b(2); // V_0
        let mut v_cur = p.clone(); // V_1
        for k in 1..60 {
            let kb = b(k);
            let (u, v, _) = lucas_uv_mod(&p, &q, &kb, &n);
            assert_eq!(u, u_cur.mod_floor(&n), "U_{k} mismatch");
            assert_eq!(v, v_cur.mod_floor(&n), "V_{k} mismatch");
            let u_next = (&p * &u_cur - &q * &u_prev).mod_floor(&n);
            let v_next = (&p * &v_cur - &q * &v_prev).mod_floor(&n);
            u_prev = u_cur;
            u_cur = u_next;
            v_prev = v_cur;
            v_cur = v_next;
        }
    }

    #[test]
    fn bpsw_accepts_primes() {
        let primes = [
            "97",
            "7919",
            "1000003",
            "2147483647",          // Mersenne M31
            "2305843009213693951", // Mersenne M61
            "1000000000000000009", // 10^18 + 9
        ];
        for s in primes {
            let n: Integer = s.parse::<num_bigint::BigInt>().unwrap().into();
            assert!(is_prime_bpsw(&n), "{s} should pass BPSW");
        }
    }

    #[test]
    fn bpsw_rejects_composites() {
        // Carmichael numbers and base-2 strong pseudoprimes must all fail.
        let composites = [
            561i64, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825265, 2047, 3277, 4033, 4681,
            8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665, 80581, 85489, 88357, 90751,
        ];
        for c in composites {
            assert!(!is_prime_bpsw(&b(c)), "{c} must be rejected by BPSW");
        }
        // Perfect squares of primes.
        for p in [41i64, 101, 1009] {
            assert!(!is_prime_bpsw(&b(p * p)));
        }
        // A large semiprime.
        let n = b(1_000_003) * b(1_000_033);
        assert!(!is_prime_bpsw(&n));
    }

    #[test]
    fn is_prime_u64_boundary() {
        assert!(is_prime_u64(2));
        assert!(is_prime_u64(u64::MAX - 58)); // 2^64 − 59
        assert!(!is_prime_u64(u64::MAX)); // 2^64 − 1 = 3·5·17·257·641·65537·6700417
        assert!(!is_prime_u64(0));
        assert!(!is_prime_u64(1));
    }

    #[cfg(test)]
    mod proptests {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// BPSW, plain Miller–Rabin and the u64 wrapper agree everywhere.
            #[test]
            fn primality_tests_agree(x in any::<u64>()) {
                let n = Integer::new(x);
                assert_eq!(is_prime(&n), is_prime_bpsw(&n), "is_prime vs BPSW at {x}");
                assert_eq!(is_prime(&n), is_prime_u64(x), "is_prime vs u64 at {x}");
            }
        }
    }
}
