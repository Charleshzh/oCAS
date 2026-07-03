//! Number-theory primitives for algebraic algorithms.
//!
//! Provides primality testing and generation, modular inverses, the Chinese
//! remainder theorem, quadratic-residue symbols, and modular square roots.
//! All routines operate on [`crate::Integer`].
//!
//! These primitives underpin polynomial factorization (Berlekamp,
//! Cantor–Zassenhaus, Hensel lifting, Zassenhaus), rational reconstruction,
//! and modular GCD algorithms.

use crate::Integer;

/// Small primes used as deterministic Miller–Rabin witnesses.
///
/// Testing against this set is *deterministic* for every `n < 3.317·10²⁴`.
/// For larger `n` the same bases still form a strong probable-prime test; the
/// chance of a composite passing all of them is astronomically small.
const MR_WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Test whether `n` is a (strong) probable prime to a single base `a`.
///
/// Implements the round of Miller–Rabin: write `n - 1 = d · 2^r` with `d` odd,
/// and check that `a^d ≡ 1 (mod n)` or `a^(d·2^j) ≡ -1 (mod n)` for some
/// `j < r`. Returns `false` if `n` is a small prime divisor of `a` handled
/// by the caller; here `1 < a < n` is assumed.
fn mr_witness(n: &Integer, d: &Integer, r: u64, a: &Integer) -> bool {
    let mut x = a.modpow(d, n);
    let one = Integer::from(1);
    let n_minus_one = n - &one;
    if x == one || x == n_minus_one {
        return true;
    }
    for _ in 1..r {
        x = (&x * &x).mod_floor(n);
        if x == n_minus_one {
            return true;
        }
    }
    false
}

/// Test whether `n` is prime.
///
/// Uses the deterministic Miller–Rabin witness set for `n < 3.317·10²⁴` and a
/// strong probable-prime test beyond that. Handles the small/trivial cases
/// `n ≤ 3` and even `n` explicitly before the main loop.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::is_prime;
///
/// assert!(is_prime(&Integer::from(97)));
/// assert!(!is_prime(&Integer::from(561)));   // Carmichael number
/// assert!(is_prime(&Integer::from(2_147_483_647_i64))); // Mersenne prime M31
/// ```
pub fn is_prime(n: &Integer) -> bool {
    // Small / trivial cases.
    if n < &Integer::from(2) {
        return false;
    }
    if *n == Integer::from(2) || *n == Integer::from(3) {
        return true;
    }
    // Even numbers (and explicit small-prime divisibility) weed out composites
    // cheaply before the expensive modular exponentiations.
    if n.is_even() {
        return false;
    }
    for &p in &[3u64, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let pb = Integer::from(p as i64);
        if *n == pb {
            return true;
        }
        if n.mod_floor(&pb).is_zero() {
            return false;
        }
    }

    // Write n - 1 = d · 2^r with d odd.
    let one = Integer::from(1);
    let n_minus_one = n - &one;
    let mut d = n_minus_one.clone();
    let mut r = 0u64;
    while d.is_even() {
        d >>= 1;
        r += 1;
    }

    for &a in &MR_WITNESSES {
        let ab = Integer::from(a as i64);
        if ab >= *n {
            continue;
        }
        if !mr_witness(n, &d, r, &ab) {
            return false;
        }
    }
    true
}

/// Return the smallest prime strictly greater than `n`.
///
/// Starts at `n + 1` (or `2` if `n < 2`) and tests successive odd candidates.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::next_prime;
///
/// assert_eq!(next_prime(&Integer::from(10)), Integer::from(11));
/// assert_eq!(next_prime(&Integer::from(13)), Integer::from(17));
/// assert_eq!(next_prime(&Integer::from(0)), Integer::from(2));
/// ```
pub fn next_prime(n: &Integer) -> Integer {
    let two = Integer::from(2);
    let three = Integer::from(3);
    if n < &two {
        return two;
    }
    let mut candidate = n + &Integer::from(1);
    if candidate == three {
        return three;
    }
    // Make candidate odd.
    if candidate.is_even() {
        candidate += &Integer::from(1);
    }
    while !is_prime(&candidate) {
        candidate += &two;
    }
    candidate
}

/// An iterator yielding successive primes starting strictly after `n`.
///
/// Useful for scanning primes during Hensel lifting (find a prime not dividing
/// the leading coefficient and keeping `f mod p` square-free).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::primes_from;
///
/// let mut it = primes_from(&Integer::from(100));
/// assert_eq!(it.next().unwrap().to_string(), "101");
/// assert_eq!(it.next().unwrap().to_string(), "103");
/// ```
pub fn primes_from(n: &Integer) -> PrimesFrom {
    PrimesFrom {
        current: if n < &Integer::from(2) {
            Integer::from(2)
        } else {
            n.clone()
        },
    }
}

/// Iterator over successive primes (see [`primes_from`]).
pub struct PrimesFrom {
    current: Integer,
}

impl Iterator for PrimesFrom {
    type Item = Integer;

    fn next(&mut self) -> Option<Integer> {
        self.current = next_prime(&self.current);
        Some(self.current.clone())
    }
}
///
/// Compute the multiplicative inverse of `a` modulo `m`, i.e. the `x` with
/// `a·x ≡ 1 (mod m)`.
///
/// Returns `None` when `gcd(a, m) ≠ 1` or `m ≤ 1`. The result lies in
/// `[0, m)`.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::mod_inv;
///
/// assert_eq!(mod_inv(&Integer::from(3), &Integer::from(11)), Some(Integer::from(4)));
/// assert_eq!(mod_inv(&Integer::from(2), &Integer::from(4)), None);
/// ```
pub fn mod_inv(a: &Integer, m: &Integer) -> Option<Integer> {
    if m <= &Integer::from(1) {
        return None;
    }
    let (g, x, _) = extended_gcd(a, m);
    if !g.is_one() {
        return None;
    }
    // Normalize into [0, m).
    let mut r = x.mod_floor(m);
    if r.is_negative() {
        r += m;
    }
    Some(r)
}

/// Extended Euclidean algorithm: returns `(g, x, y)` with `g = a·x + b·y`,
/// `g = gcd(a, b)` (taken non-negative).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::extended_gcd;
///
/// let (g, x, y) = extended_gcd(&Integer::from(240), &Integer::from(46));
/// assert_eq!(g, Integer::from(2));
/// assert_eq!(&x * &Integer::from(240) + &y * &Integer::from(46), g);
/// ```
pub fn extended_gcd(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = Integer::from(1);
    let mut s = Integer::from(0);
    let mut old_t = Integer::from(0);
    let mut t = Integer::from(1);

    while !r.is_zero() {
        let (q, rem) = old_r.div_rem(&r);
        old_r = r;
        r = rem;

        let qs = &q * &s;
        let new_s = &old_s - &qs;
        old_s = s;
        s = new_s;

        let qt = &q * &t;
        let new_t = &old_t - &qt;
        old_t = t;
        t = new_t;
    }
    // Ensure gcd is reported non-negative.
    if old_r.is_negative() {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }
    (old_r, old_s, old_t)
}

/// Reduce `a` into the symmetric range `(-m/2, m/2]` modulo `m`.
///
/// Used by Hensel lifting when recovering integer coefficients from their
/// modular images.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::symmetric_mod;
///
/// // mod 7, range (-3.5, 3.5]: 3 stays 3, 5 wraps to -2, 6 wraps to -1.
/// assert_eq!(symmetric_mod(&Integer::from(3), &Integer::from(7)), Integer::from(3));
/// assert_eq!(symmetric_mod(&Integer::from(5), &Integer::from(7)), Integer::from(-2));
/// assert_eq!(symmetric_mod(&Integer::from(6), &Integer::from(7)), Integer::from(-1));
/// ```
pub fn symmetric_mod(a: &Integer, m: &Integer) -> Integer {
    let half = m / &Integer::from(2);
    let mut r = a.mod_floor(m);
    if r > half {
        r -= m;
    }
    r
}

/// Combine two congruences `x ≡ r1 (mod m1)`, `x ≡ r2 (mod m2)` via the Chinese
/// remainder theorem.
///
/// Returns `(r, m)` with `m = lcm(m1, m2)` and `r ≡ r1 (mod m1)`,
/// `r ≡ r2 (mod m2)`, `0 ≤ r < m`, or `None` when the system is inconsistent
/// (i.e. `r1 - r2` is not divisible by `gcd(m1, m2)`).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::crt;
///
/// // x ≡ 2 (mod 3), x ≡ 3 (mod 5)  =>  x ≡ 8 (mod 15).
/// let (r, m) = crt(&Integer::from(2), &Integer::from(3),
///                  &Integer::from(3), &Integer::from(5)).unwrap();
/// assert_eq!(r, Integer::from(8));
/// assert_eq!(m, Integer::from(15));
/// ```
pub fn crt(r1: &Integer, m1: &Integer, r2: &Integer, m2: &Integer) -> Option<(Integer, Integer)> {
    let (g, p, _q) = extended_gcd(m1, m2);
    let diff = r1 - r2;
    // Solvable iff gcd(m1, m2) divides (r1 - r2).
    if !diff.mod_floor(&g).is_zero() {
        return None;
    }
    // lcm(m1, m2) = m1 / g * m2.
    let lcm = (m1 / &g) * m2;
    // Using the Bézout identity s·m1 + t·m2 = g (here `p` is the s for m1):
    //   r = r1 + m1 · p · ((r2 - r1) / g)
    // satisfies r ≡ r1 (mod m1) and r ≡ r2 (mod m2).
    let step = (r2 - r1) / &g;
    let mut r = r1 + &(m1 * &p * &step);
    r = r.mod_floor(&lcm);
    Some((r, lcm))
}

/// The Legendre symbol `(a / p)` for an odd prime `p`.
///
/// Returns `1` if `a` is a quadratic residue mod `p`, `-1` if it is a
/// non-residue, and `0` if `p | a`. The primality of `p` is the caller's
/// responsibility; this is just `jacobi(a, p)` for prime `p`.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::legendre;
///
/// assert_eq!(legendre(&Integer::from(2), &Integer::from(7)), 1);  // 2 is a QR mod 7
/// assert_eq!(legendre(&Integer::from(3), &Integer::from(7)), -1); // 3 is a non-QR
/// ```
pub fn legendre(a: &Integer, p: &Integer) -> i8 {
    jacobi(a, p)
}

/// The Jacobi symbol `(a / n)` for odd positive `n`.
///
/// Computed by quadratic reciprocity. Returns `0`, `1`, or `-1`. If `n` is
/// even or non-positive the result is undefined (the function returns `0`).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::jacobi;
///
/// assert_eq!(jacobi(&Integer::from(2), &Integer::from(15)), 1);
/// assert_eq!(jacobi(&Integer::from(7), &Integer::from(15)), -1);
/// ```
pub fn jacobi(a: &Integer, n: &Integer) -> i8 {
    if n.is_zero() || n.is_negative() || n.is_even() {
        return 0;
    }
    let four = Integer::from(4);
    let eight = Integer::from(8);
    let mut a = a.mod_floor(n);
    let mut n = n.clone();
    let mut t: i8 = 1;
    while !a.is_zero() {
        // Remove factors of 2 from a.
        while a.is_even() {
            a >>= 1;
            let r = n.mod_floor(&eight);
            if r == Integer::from(3) || r == Integer::from(5) {
                t = -t;
            }
        }
        // Now a is odd: swap a, n.
        std::mem::swap(&mut a, &mut n);
        if a.mod_floor(&four) == Integer::from(3) && n.mod_floor(&four) == Integer::from(3) {
            t = -t;
        }
        a = a.mod_floor(&n);
    }
    if n.is_one() { t } else { 0 }
}

/// Compute a square root of `a` modulo the odd prime `p`, i.e. an `x` with
/// `x² ≡ a (mod p)`.
///
/// Returns `None` when `a` is a quadratic non-residue (`legendre(a,p) = -1`)
/// or when `p` is not an odd prime. Handles `p ≡ 3 (mod 4)` via the fast path
/// `x = a^((p+1)/4)` and falls back to the full Tonelli–Shanks algorithm.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::mod_sqrt;
///
/// // 2 is a QR mod 7: roots are 3 and 4 (9 ≡ 2, 16 ≡ 2 mod 7).
/// let r = mod_sqrt(&Integer::from(2), &Integer::from(7)).unwrap();
/// assert!(r == Integer::from(3) || r == Integer::from(4));
/// ```
pub fn mod_sqrt(a: &Integer, p: &Integer) -> Option<Integer> {
    if p <= &Integer::from(2) {
        return None;
    }
    let a = a.mod_floor(p);
    if a.is_zero() {
        return Some(Integer::from(0));
    }
    if legendre(&a, p) != 1 {
        return None;
    }
    // Fast path: p ≡ 3 (mod 4)  =>  x = a^((p+1)/4).
    if p.mod_floor(&Integer::from(4)) == Integer::from(3) {
        let exp = (p + &Integer::from(1)) / &Integer::from(4);
        let x = a.modpow(&exp, p);
        return Some(x);
    }
    // Full Tonelli–Shanks. Write p - 1 = q · 2^s with q odd.
    let one = Integer::from(1);
    let mut q = p - &one;
    let mut s = 0u64;
    while q.is_even() {
        q >>= 1;
        s += 1;
    }
    // Find a non-residue z.
    let mut z = Integer::from(2);
    while legendre(&z, p) != -1 {
        z += &one;
    }
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + &one) / &Integer::from(2)), p);
    while !t.is_one() {
        // Find the least i, 0 < i < m, with t^(2^i) ≡ 1.
        let mut i = 0u64;
        let mut t2i = t.clone();
        while !t2i.is_one() {
            t2i = (&t2i * &t2i).mod_floor(p);
            i += 1;
            if i >= m {
                return None;
            }
        }
        let mut b = c.clone();
        for _ in 0..(m - i - 1) {
            b = (&b * &b).mod_floor(p);
        }
        m = i;
        c = (&b * &b).mod_floor(p);
        t = (&t * &c).mod_floor(p);
        r = (&r * &b).mod_floor(p);
    }
    Some(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn primality_small() {
        let primes = [2, 3, 5, 7, 11, 13, 97, 101, 997, 7919];
        for p in primes {
            assert!(is_prime(&b(p)), "{p} should be prime");
        }
        let composites = [0, 1, 4, 6, 8, 9, 15, 21, 25, 100, 1001];
        for c in composites {
            assert!(!is_prime(&b(c)), "{c} should be composite");
        }
    }

    #[test]
    fn primality_carmichael() {
        // 561 = 3·11·17, the smallest Carmichael number (fools Fermat tests).
        assert!(!is_prime(&b(561)));
        assert!(!is_prime(&b(1105)));
        assert!(!is_prime(&b(1729)));
    }

    #[test]
    fn primality_large() {
        // 2^31 - 1 is the Mersenne prime M31.
        assert!(is_prime(&Integer::from(2_147_483_647_i64)));
        // 2^67 - 1 is composite (factors: 193707721 · 761838257287).
        let big = Integer::from(2).pow_u32(67) - &Integer::from(1);
        assert!(!is_prime(&big));
    }

    #[test]
    fn next_prime_works() {
        assert_eq!(next_prime(&b(0)), b(2));
        assert_eq!(next_prime(&b(1)), b(2));
        assert_eq!(next_prime(&b(2)), b(3));
        assert_eq!(next_prime(&b(10)), b(11));
        assert_eq!(next_prime(&b(13)), b(17));
        assert_eq!(next_prime(&b(100)), b(101));
    }

    #[test]
    fn primes_from_iterator() {
        let got: Vec<String> = primes_from(&b(10)).take(5).map(|x| x.to_string()).collect();
        assert_eq!(got, vec!["11", "13", "17", "19", "23"]);
    }

    #[test]
    fn modular_inverse() {
        assert_eq!(mod_inv(&b(3), &b(11)), Some(b(4))); // 3·4 = 12 ≡ 1
        assert_eq!(mod_inv(&b(7), &b(13)), Some(b(2))); // 7·2 = 14 ≡ 1
        assert_eq!(mod_inv(&b(2), &b(4)), None);
        assert_eq!(mod_inv(&b(6), &b(9)), None);
        // Inverse is always in [0, m). -3 ≡ 8 (mod 11); 8·7 = 56 ≡ 1, so inv = 7.
        let inv = mod_inv(&b(-3), &b(11)).unwrap();
        assert_eq!(inv, b(7));
    }

    #[test]
    fn symmetric_modulo_range() {
        for a in -14..=14 {
            let r = symmetric_mod(&b(a), &b(7));
            // Range is (-3.5, 3.5], i.e. r ∈ {-3,-2,-1,0,1,2,3}.
            assert!(
                r > b(-4) && r <= b(3),
                "symmetric_mod({a}, 7) = {r} out of range"
            );
            assert_eq!((&r * &r).mod_floor(&b(7)), (b(a) * b(a)).mod_floor(&b(7)));
        }
        assert_eq!(symmetric_mod(&b(6), &b(7)), b(-1));
        assert_eq!(symmetric_mod(&b(5), &b(7)), b(-2));
        assert_eq!(symmetric_mod(&b(3), &b(7)), b(3));
    }

    #[test]
    fn crt_basic() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5) => x ≡ 8 (mod 15).
        let (r, m) = crt(&b(2), &b(3), &b(3), &b(5)).unwrap();
        assert_eq!(r, b(8));
        assert_eq!(m, b(15));
        // Sunzi: x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) => x ≡ 23 (mod 105).
        let (r, m) = crt(&b(2), &b(3), &b(3), &b(5)).unwrap();
        let (r, m) = crt(&r, &m, &b(2), &b(7)).unwrap();
        assert_eq!(r, b(23));
        assert_eq!(m, b(105));
    }

    #[test]
    fn crt_inconsistent() {
        // x ≡ 1 (mod 4), x ≡ 2 (mod 4): inconsistent, gcd(4,4)=4 does not divide 1-2.
        assert!(crt(&b(1), &b(4), &b(2), &b(4)).is_none());
    }

    #[test]
    fn crt_non_coprime_compatible() {
        // m1=4, m2=6, gcd=2. x ≡ 1 (mod 4) and x ≡ 3 (mod 6): both hold for x=9, lcm=12.
        let (r, m) = crt(&b(1), &b(4), &b(3), &b(6)).unwrap();
        assert_eq!(m, b(12));
        assert_eq!(r.mod_floor(&b(4)), b(1));
        assert_eq!(r.mod_floor(&b(6)), b(3));
    }

    #[test]
    fn legendre_symbols() {
        // (2/7) = 1, (3/7) = -1, (5/7) = -1, (6/7) = -1.
        assert_eq!(legendre(&b(2), &b(7)), 1);
        assert_eq!(legendre(&b(3), &b(7)), -1);
        assert_eq!(legendre(&b(5), &b(7)), -1);
        assert_eq!(legendre(&b(0), &b(7)), 0);
    }

    #[test]
    fn jacobi_reciprocity() {
        // (2/15) = (2/3)(2/5) = (-1)(-1) = 1.
        assert_eq!(jacobi(&b(2), &b(15)), 1);
        // Cross-check the Jacobi symbol against Euler's criterion for prime
        // moduli: (a/p) ≡ a^((p-1)/2) (mod p), mapped to {-1, 0, 1}.
        let primes = [7u64, 11, 13, 17, 23, 41, 101, 1009, 9907];
        for &p in &primes {
            let pb = Integer::from(p as i64);
            assert!(is_prime(&pb), "{p} assumed prime");
            let exp = (&pb - &Integer::from(1)) / &Integer::from(2);
            for a in 0..p.min(60) {
                let ab = b(a as i64);
                let expected = match ab.modpow(&exp, &pb).to_string().as_str() {
                    "0" => 0,
                    "1" => 1,
                    _ => -1, // p-1
                };
                assert_eq!(jacobi(&ab, &pb), expected, "jacobi({a}/{p}) mismatch");
            }
        }
    }

    #[test]
    fn mod_sqrt_residue() {
        // mod 7: QRs are {0,1,2,4}. Roots of 2 are 3,4.
        for a in [b(0), b(1), b(2), b(4)] {
            let x = mod_sqrt(&a, &b(7)).unwrap();
            assert_eq!((&x * &x).mod_floor(&b(7)), a.mod_floor(&b(7)));
        }
        // non-residues
        assert!(mod_sqrt(&b(3), &b(7)).is_none());
        assert!(mod_sqrt(&b(5), &b(7)).is_none());
    }

    #[test]
    fn mod_sqrt_fast_path_p3_mod4() {
        // 11 ≡ 3 (mod 4): 3 is a QR, 5² = 25 ≡ 3 (mod 11).
        let x = mod_sqrt(&b(3), &b(11)).unwrap();
        assert_eq!((&x * &x).mod_floor(&b(11)), b(3));
    }

    #[test]
    fn mod_sqrt_general_prime() {
        // 17 ≡ 1 (mod 4): exercise the full Tonelli–Shanks path.
        // QRs mod 17: 1,2,4,8,9,13,15,16. Verify 2 has a square root.
        let x = mod_sqrt(&b(2), &b(17)).unwrap();
        assert_eq!((&x * &x).mod_floor(&b(17)), b(2));
        // Larger p ≡ 1 (mod 4): 41. 4·? Check several residues.
        for a in 0..41 {
            let ab = b(a);
            let r = (&ab * &ab).mod_floor(&b(41));
            let s = mod_sqrt(&r, &b(41)).unwrap();
            assert_eq!((&s * &s).mod_floor(&b(41)), r);
        }
    }

    #[test]
    fn extended_gcd_identity() {
        let (g, x, y) = extended_gcd(&b(240), &b(46));
        assert_eq!(g, b(2));
        assert_eq!(&x * &b(240) + &y * &b(46), b(2));
        // gcd of coprime numbers is 1.
        let (g, _, _) = extended_gcd(&b(17), &b(13));
        assert_eq!(g, b(1));
        // gcd is non-negative even for negative inputs.
        let (g, _, _) = extended_gcd(&b(-240), &b(46));
        assert_eq!(g, b(2));
    }
}
