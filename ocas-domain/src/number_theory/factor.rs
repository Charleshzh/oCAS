//! Integer factorization: trial division, Pollard's rho (Brent's variant),
//! Pollard's p−1, Williams' p+1, and Lenstra's elliptic-curve method (ECM).
//!
//! [`factor_integer`] is the high-level entry point. It peels off small
//! factors by trial division, then splits composite cofactors with an
//! escalating rho → p−1 → p+1 → ECM strategy with growing smoothness bounds.
//! Randomized methods take an explicit RNG so tests can be reproducible.

use rand::RngCore;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::primes::{is_prime_bpsw, lucas_uv_mod};
use super::{jacobi, mod_inv};
use crate::{EuclideanDomain, Integer, IntegerDomain};

/// GCD of two integers via the Euclidean domain (always non-negative).
fn gcd_int(a: &Integer, b: &Integer) -> Integer {
    IntegerDomain.gcd(&a.abs(), &b.abs())
}

/// All primes `≤ limit` via a sieve of Eratosthenes.
fn primes_up_to(limit: u64) -> Vec<u64> {
    if limit < 2 {
        return Vec::new();
    }
    let mut sieve = vec![true; (limit + 1) as usize];
    sieve[0] = false;
    sieve[1] = false;
    let mut p = 2u64;
    while p * p <= limit {
        if sieve[p as usize] {
            let mut multiple = p * p;
            while multiple <= limit {
                sieve[multiple as usize] = false;
                multiple += p;
            }
        }
        p += 1;
    }
    sieve
        .iter()
        .enumerate()
        .filter_map(|(i, &is)| is.then_some(i as u64))
        .collect()
}

/// Largest power `q^e ≤ bound` (with `q^e = q` when `q² > bound`).
fn prime_power_le(q: u64, bound: u64) -> u64 {
    let mut e = q;
    while e <= bound / q {
        e *= q;
    }
    e
}

/// Divide out all prime factors `≤ limit`, returning `(factors, cofactor)`.
///
/// The cofactor has no prime factor `≤ limit` (it need not be prime).
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::factor_trial;
///
/// let (factors, rest) = factor_trial(&Integer::from(2 * 2 * 3 * 7 * 1000003), 100);
/// assert_eq!(factors, vec![
///     (Integer::from(2), 2),
///     (Integer::from(3), 1),
///     (Integer::from(7), 1),
/// ]);
/// assert_eq!(rest, Integer::from(1000003));
/// ```
pub fn factor_trial(n: &Integer, limit: u64) -> (Vec<(Integer, u32)>, Integer) {
    let mut m = n.abs();
    let mut factors = Vec::new();
    for p in primes_up_to(limit) {
        if m.is_one() {
            break;
        }
        let pb = Integer::new(p);
        let mut e = 0u32;
        loop {
            let (quot, rem) = m.div_rem(&pb);
            if !rem.is_zero() {
                break;
            }
            m = quot;
            e += 1;
        }
        if e > 0 {
            factors.push((pb, e));
        }
    }
    (factors, m)
}

/// Brent's variant of Pollard's rho: search for a nontrivial factor of the
/// odd composite `n`. Retries with fresh random parameters a bounded number
/// of times; returns `None` only when repeatedly unlucky.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::pollard_rho_brent;
/// use rand::SeedableRng;
/// use rand_xoshiro::Xoshiro256PlusPlus;
///
/// let n = Integer::from(1000003) * Integer::from(1000033);
/// let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
/// let d = pollard_rho_brent(&n, &mut rng).unwrap();
/// assert!(d == Integer::from(1000003) || d == Integer::from(1000033));
/// ```
pub fn pollard_rho_brent(n: &Integer, rng: &mut Xoshiro256PlusPlus) -> Option<Integer> {
    if n.is_even() {
        return Some(Integer::from(2));
    }
    let one = Integer::from(1);
    for _restart in 0..2 {
        let c = Integer::from((((rng.next_u64() >> 1) % 1_000_003) + 1) as i64);
        let mut y = Integer::from(((rng.next_u64() >> 1) % 1_000_033) as i64);
        let m: u64 = 128;
        let mut g = one.clone();
        let mut r: u64 = 1;
        let mut x = Integer::from(0);
        let mut ys = Integer::from(0);
        let mut stalled = false;
        while g.is_one() {
            x = y.clone();
            for _ in 0..r {
                y = (&y * &y + &c).mod_floor(n);
            }
            let mut k = 0u64;
            while k < r && g.is_one() {
                ys = y.clone();
                let mut q = one.clone();
                for _ in 0..m.min(r - k) {
                    y = (&y * &y + &c).mod_floor(n);
                    let diff = (&x - &y).abs();
                    q = (&q * &diff).mod_floor(n);
                }
                g = gcd_int(&q, n);
                k += m;
            }
            r <<= 1;
            if r > (1u64 << 18) {
                stalled = true;
                break;
            }
        }
        if stalled {
            continue;
        }
        if g == *n {
            // Overshot: backtrack one step at a time from the last block.
            loop {
                ys = (&ys * &ys + &c).mod_floor(n);
                g = gcd_int(&(&x - &ys).abs(), n);
                if !g.is_one() {
                    break;
                }
            }
        }
        if g != *n {
            return Some(g);
        }
    }
    None
}

/// Pollard's p−1 method (stage 1): find a factor `p` of odd composite `n`
/// when `p − 1` is `b1`-powersmooth. Computes `a^M mod n` with
/// `M = ∏ q^e ≤ b1` and tests `gcd(a^M − 1, n)` periodically.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::pollard_pm1;
/// use rand::SeedableRng;
/// use rand_xoshiro::Xoshiro256PlusPlus;
///
/// // 65537 is prime and 65536 = 2^16 is 2^17-powersmooth.
/// let n = Integer::from(65537) * Integer::from(1000003);
/// let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
/// let d = pollard_pm1(&n, 1 << 17, &mut rng).unwrap();
/// assert_eq!(n.mod_floor(&d), Integer::from(0));
/// ```
pub fn pollard_pm1(n: &Integer, b1: u64, rng: &mut Xoshiro256PlusPlus) -> Option<Integer> {
    let one = Integer::from(1);
    let primes = primes_up_to(b1);
    for _attempt in 0..2 {
        let a0 = Integer::from((2 + ((rng.next_u64() >> 1) % 1_000_000)) as i64);
        let g0 = gcd_int(&a0, n);
        if g0 > one && g0 < *n {
            return Some(g0);
        }
        if g0 == *n {
            continue;
        }
        let mut a = a0;
        for (idx, &q) in primes.iter().enumerate() {
            a = a.modpow(&Integer::new(prime_power_le(q, b1)), n);
            if idx % 32 == 31 {
                let g = gcd_int(&(&a - &one), n);
                if g > one && g < *n {
                    return Some(g);
                }
                if g == *n {
                    break;
                }
            }
        }
        let g = gcd_int(&(&a - &one), n);
        if g > one && g < *n {
            return Some(g);
        }
        if g.is_one() {
            return None; // no B1-powersmooth p−1 among factors
        }
        // g == n: unlucky base (all factors smooth) — retry with another.
    }
    None
}

/// Williams' p+1 method (stage 1): find a factor `p` of odd composite `n`
/// when `p + 1` is `b1`-powersmooth. Uses the Lucas `V` sequence with
/// `Q = 1` and random `P` with `jacobi(P² − 4, n) = −1`.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::williams_pp1;
/// use rand::SeedableRng;
/// use rand_xoshiro::Xoshiro256PlusPlus;
///
/// // 31 is prime and 31 + 1 = 2^5 is 64-powersmooth.
/// let n = Integer::from(31) * Integer::from(1000003);
/// let mut rng = Xoshiro256PlusPlus::seed_from_u64(13);
/// let d = williams_pp1(&n, 64, &mut rng).unwrap();
/// assert_eq!(n.mod_floor(&d), Integer::from(0));
/// ```
pub fn williams_pp1(n: &Integer, b1: u64, rng: &mut Xoshiro256PlusPlus) -> Option<Integer> {
    let one = Integer::from(1);
    let two = Integer::from(2);
    let primes = primes_up_to(b1);
    for _attempt in 0..2 {
        let p0 = Integer::from((3 + ((rng.next_u64() >> 1) % 1_000_000)) as i64);
        let disc = &p0 * &p0 - &Integer::from(4);
        let gd = gcd_int(&disc, n);
        if gd > one && gd < *n {
            return Some(gd);
        }
        if gd == *n || jacobi(&disc, n) != -1 {
            continue; // need the p+1 twist
        }
        let mut v = p0;
        let mut unlucky = false;
        for (idx, &q) in primes.iter().enumerate() {
            let e = Integer::new(prime_power_le(q, b1));
            let (_, v_new, _) = lucas_uv_mod(&v, &one, &e, n);
            v = v_new;
            if idx % 32 == 31 {
                let g = gcd_int(&(&v - &two), n);
                if g > one && g < *n {
                    return Some(g);
                }
                if g == *n {
                    unlucky = true;
                    break;
                }
            }
        }
        if unlucky {
            continue;
        }
        let g = gcd_int(&(&v - &two), n);
        if g > one && g < *n {
            return Some(g);
        }
        if g.is_one() {
            return None; // no B1-powersmooth p+1 among factors
        }
    }
    None
}

/// Point in `(X : Z)` projective coordinates on a Montgomery curve.
struct ProjPoint {
    x: Integer,
    z: Integer,
}

/// Montgomery doubling: `[2](X : Z)` with curve constant `a24 = (A + 2)/4`.
fn ecm_double(p: &ProjPoint, a24: &Integer, n: &Integer) -> ProjPoint {
    let up = (&p.x + &p.z).mod_floor(n);
    let up = (&up * &up).mod_floor(n); // (X + Z)²
    let um = (&p.x - &p.z).mod_floor(n);
    let um = (&um * &um).mod_floor(n); // (X − Z)²
    let t = (&up - &um).mod_floor(n); // 4XZ
    let x2 = (&up * &um).mod_floor(n);
    let e = (&um + &(a24 * &t)).mod_floor(n);
    let z2 = (&t * &e).mod_floor(n);
    ProjPoint { x: x2, z: z2 }
}

/// Montgomery differential addition: `P + Q` given `P − Q = diff`.
fn ecm_add(p: &ProjPoint, q: &ProjPoint, diff: &ProjPoint, n: &Integer) -> ProjPoint {
    let a = (&p.x - &p.z).mod_floor(n);
    let b = (&p.x + &p.z).mod_floor(n);
    let c = (&q.x - &q.z).mod_floor(n);
    let d = (&q.x + &q.z).mod_floor(n);
    let da = (&d * &a).mod_floor(n);
    let cb = (&c * &b).mod_floor(n);
    let plus = (&da + &cb).mod_floor(n);
    let minus = (&da - &cb).mod_floor(n);
    let x = (&diff.z * &(&plus * &plus)).mod_floor(n);
    let z = (&diff.x * &(&minus * &minus)).mod_floor(n);
    ProjPoint { x, z }
}

/// Montgomery ladder: scalar multiplication `[k]P` for `k ≥ 1`, x-coordinate
/// only, with the base point as the fixed difference.
fn ecm_mul(p: &ProjPoint, k: u64, a24: &Integer, n: &Integer) -> ProjPoint {
    debug_assert!(k >= 1);
    let mut r0 = ProjPoint {
        x: p.x.clone(),
        z: p.z.clone(),
    };
    let mut r1 = ecm_double(p, a24, n);
    let top = 63 - k.leading_zeros();
    for i in (0..top).rev() {
        if (k >> i) & 1 == 0 {
            r1 = ecm_add(&r0, &r1, p, n);
            r0 = ecm_double(&r0, a24, n);
        } else {
            r0 = ecm_add(&r0, &r1, p, n);
            r1 = ecm_double(&r1, a24, n);
        }
    }
    r0
}

/// Outcome of the Suyama curve setup: either a curve, a lucky factor from
/// the denominator gcd, or a degenerate parameter to discard.
enum Suyama {
    Curve(ProjPoint, Integer),
    Factor(Integer),
    Retry,
}

/// Suyama's parametrization: from `σ ∉ {0, 1, 5}` build the Montgomery
/// curve with `a24 = (A + 2)/4`, `A = ((v−u)³(3u+v))/(4u³v) − 2`, and the
/// rational point `(u³ : v³)`, where `u = σ² − 5`, `v = 4σ`.
fn suyama_curve(sigma: u64, n: &Integer) -> Suyama {
    let one = Integer::from(1);
    let s = Integer::new(sigma);
    let u = (&s * &s - &Integer::from(5)).mod_floor(n);
    let v = (&Integer::from(4) * &s).mod_floor(n);
    let u3 = (&u * &u * &u).mod_floor(n);
    let v3 = (&v * &v * &v).mod_floor(n);
    let den = (&Integer::from(4) * &u3 * &v).mod_floor(n);
    let g = gcd_int(&den, n);
    if g == *n {
        return Suyama::Retry;
    }
    if g > one {
        return Suyama::Factor(g);
    }
    // a24 = (A + 2)/4 = numerator / (4 · denominator).
    let vu = (&v - &u).mod_floor(n);
    let num = (&vu * &vu * &vu * &(&Integer::from(3) * &u + &v)).mod_floor(n);
    let den4 = (&Integer::from(4) * &den).mod_floor(n);
    let inv = match mod_inv(&den4, n) {
        Some(i) => i,
        None => return Suyama::Retry,
    };
    let a24 = (&num * &inv).mod_floor(n);
    Suyama::Curve(ProjPoint { x: u3, z: v3 }, a24)
}

/// Lenstra's elliptic-curve method, stage 1: try up to `max_curves` Suyama
/// curves with smoothness bound `b1`. Finds a factor `p` of `n` when the
/// group order of some curve over `𝔽_p` is `b1`-powersmooth.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::ecm;
/// use rand::SeedableRng;
/// use rand_xoshiro::Xoshiro256PlusPlus;
///
/// let n = Integer::from(1000003) * Integer::from(1000033);
/// let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);
/// let d = ecm(&n, 2_000, 50, &mut rng).unwrap();
/// assert!(d == Integer::from(1000003) || d == Integer::from(1000033));
/// ```
pub fn ecm(n: &Integer, b1: u64, max_curves: u32, rng: &mut Xoshiro256PlusPlus) -> Option<Integer> {
    let one = Integer::from(1);
    let primes = primes_up_to(b1);
    for _curve in 0..max_curves {
        let sigma = 6 + (rng.next_u64() >> 1) % 1_000_000;
        let (mut pt, a24) = match suyama_curve(sigma, n) {
            Suyama::Factor(g) => return Some(g),
            Suyama::Retry => continue,
            Suyama::Curve(p, a) => (p, a),
        };
        let mut failed = false;
        for (idx, &q) in primes.iter().enumerate() {
            pt = ecm_mul(&pt, prime_power_le(q, b1), &a24, n);
            if idx % 64 == 63 {
                let g = gcd_int(&pt.z, n);
                if g > one && g < *n {
                    return Some(g);
                }
                if g == *n {
                    failed = true;
                    break;
                }
            }
        }
        if failed {
            continue;
        }
        let g = gcd_int(&pt.z, n);
        if g > one && g < *n {
            return Some(g);
        }
    }
    None
}

/// Record one prime factor occurrence (merging multiplicities).
fn record(factors: &mut Vec<(Integer, u32)>, p: Integer) {
    if let Some(entry) = factors.iter_mut().find(|(q, _)| *q == p) {
        entry.1 += 1;
    } else {
        factors.push((p, 1));
    }
}

/// Find one nontrivial factor of composite `n` by escalating methods with
/// growing smoothness bounds. Loops until successful.
fn find_factor(n: &Integer, rng: &mut Xoshiro256PlusPlus) -> Integer {
    // One quick rho attempt: cheap for small factors, pointless to repeat
    // once the ECM bounds grow (rho's expected cost is O(√p)).
    if let Some(d) = pollard_rho_brent(n, rng) {
        return d;
    }
    let mut b1: u64 = 2_000;
    loop {
        if let Some(d) = pollard_pm1(n, b1, rng) {
            return d;
        }
        if let Some(d) = williams_pp1(n, b1, rng) {
            return d;
        }
        // Curve budget scaled to the smoothness bound (≈ B1/550).
        let curves = (b1 / 550).clamp(10, 300) as u32;
        if let Some(d) = ecm(n, b1, curves, rng) {
            return d;
        }
        b1 = b1.saturating_mul(4);
    }
}

/// Fully factor `|n|` into primes, returning `(prime, exponent)` pairs in
/// ascending order. `n ∈ {0, ±1}` yields an empty list; the sign is ignored.
///
/// Strategy: trial division up to 1000, then each composite cofactor is
/// split by [`find_factor`] (rho → p−1 → p+1 → ECM with growing bounds) and
/// recursed until every leaf passes the BPSW primality test.
///
/// # Example
///
/// ```
/// use ocas_domain::Integer;
/// use ocas_domain::number_theory::factor::factor_integer;
///
/// let f = factor_integer(&Integer::from(2 * 2 * 3 * 5 * 101 * 1000003));
/// assert_eq!(f, vec![
///     (Integer::from(2), 2),
///     (Integer::from(3), 1),
///     (Integer::from(5), 1),
///     (Integer::from(101), 1),
///     (Integer::from(1000003), 1),
/// ]);
/// ```
pub fn factor_integer(n: &Integer) -> Vec<(Integer, u32)> {
    let mut rng = Xoshiro256PlusPlus::from_rng(&mut rand::rng());
    factor_integer_with_rng(n, &mut rng)
}

/// [`factor_integer`] with an explicit RNG for reproducible runs.
pub fn factor_integer_with_rng(n: &Integer, rng: &mut Xoshiro256PlusPlus) -> Vec<(Integer, u32)> {
    let mut factors: Vec<(Integer, u32)> = Vec::new();
    let m0 = n.abs();
    if m0 <= one() {
        return factors;
    }
    let (trial, rem) = factor_trial(&m0, 1_000);
    factors.extend(trial);
    let mut stack = Vec::new();
    if !rem.is_one() {
        stack.push(rem);
    }
    while let Some(m) = stack.pop() {
        if m.is_one() {
            continue;
        }
        if is_prime_bpsw(&m) {
            record(&mut factors, m);
            continue;
        }
        let d = find_factor(&m, rng);
        stack.push(d.clone());
        stack.push(&m / &d);
    }
    factors.sort_by(|a, b| a.0.cmp(&b.0));
    factors
}

fn one() -> Integer {
    Integer::from(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(n: i64) -> Integer {
        Integer::from(n)
    }

    fn check_factorization(n: &Integer, factors: &[(Integer, u32)]) {
        let mut product = Integer::from(1);
        for (p, e) in factors {
            assert!(is_prime_bpsw(p), "factor {p} of {n} is not prime");
            product *= &p.pow_u32(*e);
        }
        assert_eq!(&product, &n.abs(), "factors do not reconstruct {n}");
        // Sorted ascending.
        for w in factors.windows(2) {
            assert!(w[0].0 < w[1].0);
        }
    }

    #[test]
    fn trial_division_peels_small_primes() {
        let n = b(2).pow_u32(10) * b(3).pow_u32(5) * b(1_000_003);
        let (factors, rest) = factor_trial(&n, 100);
        assert_eq!(factors, vec![(b(2), 10), (b(3), 5)]);
        assert_eq!(rest, b(1_000_003));
        // Fully factored by trial division.
        let (factors, rest) = factor_trial(&b(2 * 3 * 5 * 7 * 11 * 13), 100);
        assert!(rest.is_one());
        assert_eq!(factors.len(), 6);
        // Nothing to peel.
        let (factors, rest) = factor_trial(&b(1_000_003), 100);
        assert!(factors.is_empty());
        assert_eq!(rest, b(1_000_003));
    }

    #[test]
    fn rho_finds_factors() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        for (p, q) in [
            (101i64, 103i64),
            (10_007, 10_009),
            (1_000_003, 1_000_033),
            (999_983, 1_000_003),
        ] {
            let n = b(p) * b(q);
            let d = pollard_rho_brent(&n, &mut rng).expect("rho should split");
            assert!(d > b(1) && d < n);
            assert!(n.mod_floor(&d).is_zero());
        }
    }

    #[test]
    fn pm1_smooth_factor() {
        // 65537 is prime; 65536 = 2^16 is powersmooth for any B1 ≥ 65536.
        let n = b(65_537) * b(1_000_003);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
        let d = pollard_pm1(&n, 1 << 17, &mut rng).expect("p−1 should split");
        assert_eq!(d, b(65_537));
    }

    #[test]
    fn pp1_smooth_factor() {
        // 31 is prime; 32 = 2^5 is powersmooth for B1 ≥ 32.
        let n = b(31) * b(1_000_003);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
        let d = williams_pp1(&n, 64, &mut rng).expect("p+1 should split");
        assert_eq!(d, b(31));
    }

    #[test]
    fn ecm_finds_factor() {
        // 9-digit factors: comfortable for stage-1 ECM with B1 = 2000.
        let p = super::super::next_prime(&b(100_000_000));
        let q = super::super::next_prime(&b(100_003_000));
        let n = &p * &q;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let d = ecm(&n, 2_000, 200, &mut rng).expect("ECM should split");
        assert!(d == p || d == q, "ECM returned non-factor {d}");
    }

    #[test]
    fn factor_integer_basic() {
        let f = factor_integer(&b(2 * 2 * 3 * 5 * 101 * 1_000_003));
        check_factorization(&b(2 * 2 * 3 * 5 * 101 * 1_000_003), &f);
        assert_eq!(
            f,
            vec![
                (b(2), 2),
                (b(3), 1),
                (b(5), 1),
                (b(101), 1),
                (b(1_000_003), 1)
            ]
        );
    }

    #[test]
    fn factor_integer_edge_cases() {
        assert!(factor_integer(&b(0)).is_empty());
        assert!(factor_integer(&b(1)).is_empty());
        assert!(factor_integer(&b(-1)).is_empty());
        // Sign is ignored.
        assert_eq!(factor_integer(&b(-12)), vec![(b(2), 2), (b(3), 1)]);
        // Prime input.
        assert_eq!(factor_integer(&b(97)), vec![(b(97), 1)]);
        // Prime power.
        assert_eq!(factor_integer(&b(7).pow_u32(6)), vec![(b(7), 6)]);
    }

    #[test]
    fn factor_integer_semiprimes() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(17);
        for (p, q) in [
            (1_000_003i64, 1_000_033i64),
            (999_983, 999_979),
            (10_000_019, 10_000_079),
        ] {
            let n = b(p) * b(q);
            let f = factor_integer_with_rng(&n, &mut rng);
            check_factorization(&n, &f);
        }
    }

    #[test]
    fn factor_integer_is_deterministic_with_seed() {
        let n = b(1_000_003) * b(1_000_033) * b(101);
        let mut r1 = Xoshiro256PlusPlus::seed_from_u64(123);
        let mut r2 = Xoshiro256PlusPlus::seed_from_u64(123);
        assert_eq!(
            factor_integer_with_rng(&n, &mut r1),
            factor_integer_with_rng(&n, &mut r2)
        );
    }

    #[test]
    fn factor_integer_stress_small_composites() {
        // 60 pseudorandom composites up to ~2^40: full correctness check.
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(2024);
        for _ in 0..60 {
            let p = super::super::next_prime(&b(((rng.next_u64() >> 1) % 1_000_000 + 2) as i64));
            let q = super::super::next_prime(&b(((rng.next_u64() >> 1) % 1_000_000 + 2) as i64));
            let r = super::super::next_prime(&b(((rng.next_u64() >> 1) % 10_000 + 2) as i64));
            let n = &p * &q * &r;
            let f = factor_integer_with_rng(&n, &mut rng);
            check_factorization(&n, &f);
        }
    }

    /// 0.21.0 acceptance: ECM/factor pipeline splits a 30-digit semiprime
    /// within 10 seconds (release mode).
    #[test]
    #[ignore = "performance acceptance: run with --release --ignored"]
    fn factor_integer_30_digit_semiprime_under_10s() {
        let p = super::super::next_prime(&b(123_456_789_012_349));
        let q = super::super::next_prime(&b(987_654_321_098_761));
        assert_eq!(p.to_string().len(), 15);
        assert_eq!(q.to_string().len(), 15);
        let n = &p * &q;
        assert!(n.to_string().len() >= 29, "n = {n}");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        let start = std::time::Instant::now();
        let f = factor_integer_with_rng(&n, &mut rng);
        let elapsed = start.elapsed();
        check_factorization(&n, &f);
        assert_eq!(f.len(), 2, "expected exactly two prime factors: {f:?}");
        assert!(
            elapsed.as_secs() < 10,
            "30-digit semiprime took {elapsed:?} (limit 10s)"
        );
    }
}
