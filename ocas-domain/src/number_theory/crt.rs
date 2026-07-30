//! Multi-modulus Chinese remainder accumulation.
//!
//! Builds on the pairwise [`crt`] primitive: a list of congruences is merged
//! left to right, so moduli need not be pairwise coprime and inconsistencies
//! are reported as `None`.

use super::crt;
use crate::Integer;

/// Combine congruences `x ≡ r_i (mod m_i)` into a single `x ≡ R (mod M)`.
///
/// Returns `(R, M)` with `M = lcm(m_1, …, m_k)` and `0 ≤ R < M`, or `None`
/// when the system is inconsistent (some pair violates the pairwise
/// solvability condition) or the list is empty.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::crt::crt_many;
///
/// // Sunzi's problem: x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7).
/// let cs = [
///     (Integer::from(2), Integer::from(3)),
///     (Integer::from(3), Integer::from(5)),
///     (Integer::from(2), Integer::from(7)),
/// ];
/// let (r, m) = crt_many(&cs).unwrap();
/// assert_eq!(r, Integer::from(23));
/// assert_eq!(m, Integer::from(105));
/// ```
pub fn crt_many(congruences: &[(Integer, Integer)]) -> Option<(Integer, Integer)> {
    let mut iter = congruences.iter();
    let first = iter.next()?;
    let mut r = first.0.mod_floor(&first.1);
    let mut m = first.1.clone();
    for (ri, mi) in iter {
        let (r2, m2) = crt(&r, &m, ri, mi)?;
        r = r2;
        m = m2;
    }
    Some((r, m))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn crt_many_basic() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) => x ≡ 23 (mod 105).
        let cs = [(b(2), b(3)), (b(3), b(5)), (b(2), b(7))];
        let (r, m) = crt_many(&cs).unwrap();
        assert_eq!(r, b(23));
        assert_eq!(m, b(105));
    }

    #[test]
    fn crt_many_single_and_empty() {
        assert!(crt_many(&[]).is_none());
        let (r, m) = crt_many(&[(b(7), b(11))]).unwrap();
        assert_eq!(r, b(7));
        assert_eq!(m, b(11));
        // Residue is normalized into [0, m).
        let (r, _) = crt_many(&[(b(-4), b(11))]).unwrap();
        assert_eq!(r, b(7));
    }

    #[test]
    fn crt_many_inconsistent() {
        // x ≡ 1 (mod 4) and x ≡ 2 (mod 4): inconsistent.
        assert!(crt_many(&[(b(1), b(4)), (b(2), b(4))]).is_none());
        // Inconsistency detected when merging a later congruence.
        let cs = [(b(1), b(2)), (b(2), b(4)), (b(3), b(5))];
        assert!(crt_many(&cs).is_none());
    }

    #[test]
    fn crt_many_non_coprime() {
        // m1 = 4, m2 = 6: lcm = 12; x ≡ 1 (mod 4), x ≡ 3 (mod 6) => x = 9.
        let (r, m) = crt_many(&[(b(1), b(4)), (b(3), b(6))]).unwrap();
        assert_eq!(m, b(12));
        assert_eq!(r, b(9));
    }

    #[test]
    fn crt_many_recovers_value() {
        // Pick x = 12345 and coprime moduli; residues must reconstruct x.
        let x = b(12345);
        let moduli = [b(97), b(101), b(103), b(107)];
        let cs: Vec<(Integer, Integer)> =
            moduli.iter().map(|m| (x.mod_floor(m), m.clone())).collect();
        let (r, m) = crt_many(&cs).unwrap();
        let expected_m: Integer = moduli.iter().fold(b(1), |a, m| &a * m);
        assert_eq!(m, expected_m);
        assert_eq!(r, x.mod_floor(&m));
    }

    #[test]
    fn crt_many_large_moduli() {
        // Reconstruct a 40-digit coefficient from 8 prime images.
        let x: Integer = "1234567890123456789012345678901234567890"
            .parse::<num_bigint::BigInt>()
            .unwrap()
            .into();
        let primes = [
            1000003u64, 1000033, 1000037, 1000039, 1000081, 1000099, 1000117, 1000121,
        ];
        let cs: Vec<(Integer, Integer)> = primes
            .iter()
            .map(|&p| {
                let pb = Integer::new(p);
                (x.mod_floor(&pb), pb)
            })
            .collect();
        let (r, m) = crt_many(&cs).unwrap();
        assert_eq!(r, x.mod_floor(&m));
        assert_eq!(r, x); // product of the 8 primes exceeds |x|
    }
}
