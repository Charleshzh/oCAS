//! Number Theoretic Transform (NTT) for fast polynomial multiplication
//! over ℤ_p (prime finite fields).
//!
//! The NTT is the finite-field analogue of the FFT. For a prime $p$ such
//! that $p - 1$ is divisible by a sufficiently large power of 2, we can
//! perform convolution in $O(n \log n)$ instead of $O(n^2)$ (schoolbook)
//! or $O(n^{1.585})$ (Karatsuba).
//!
//! This module implements a radix-2 Cooley-Tukey NTT with bit-reversal
//! permutation. It only activates for primes where `p - 1` has a 2-power
//! factor ≥ the next power of 2 above the polynomial degree.

use num_bigint::BigInt;

/// Minimum polynomial degree to trigger NTT multiplication.
pub(crate) const NTT_THRESHOLD: usize = 256;

/// Returns `true` if `p - 1` is divisible by `n` (i.e. an `n`-th root of
/// unity exists in ℤ_p).  `n` must be a power of 2.
pub fn is_ntt_friendly(p: u64, n: usize) -> bool {
    if n == 0 {
        return true;
    }
    // p - 1 must be divisible by n
    let pm1 = p - 1;
    pm1 % (n as u64) == 0
}

/// Find the smallest primitive `n`-th root of unity in ℤ_p.
///
/// Returns `None` if no such root exists (i.e. `p - 1` is not divisible
/// by `n`).
pub fn find_primitive_root(p: u64, n: usize) -> Option<u64> {
    if !is_ntt_friendly(p, n) {
        return None;
    }
    let pm1 = p - 1;
    let exponent = pm1 / (n as u64);

    // Try small candidates g; return g^exponent mod p
    for g in 2u64.. {
        if g >= p {
            break;
        }
        let root = modpow(g, exponent, p);
        // root must be a primitive n-th root: root^n ≡ 1 (mod p) and
        // root^(n/q) ≢ 1 for every prime q | n. For power-of-2 n, we
        // only need root^(n/2) ≢ 1.
        if root == 1 {
            continue;
        }
        if n > 1 {
            let half = modpow(root, (n / 2) as u64, p);
            if half == 1 {
                continue;
            }
        }
        return Some(root);
    }
    None
}

/// Modular exponentiation: `base^exp mod p` using binary exponentiation.
///
/// Uses `u128` intermediate to avoid overflow when multiplying two `u64`
/// values whose product may approach $2^{128}$.
#[inline]
pub fn modpow(mut base: u64, mut exp: u64, p: u64) -> u64 {
    let p128 = p as u128;
    base %= p;
    let mut result: u64 = 1;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % p128) as u64;
        }
        base = ((base as u128 * base as u128) % p128) as u64;
        exp >>= 1;
    }
    result
}

/// Modular multiplication: `(a * b) mod p` using `u128` intermediate.
#[inline]
fn modmul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % (p as u128)) as u64
}

// ---------------------------------------------------------------------------
// Montgomery modular arithmetic
// ---------------------------------------------------------------------------

/// Precomputed context for Montgomery modular multiplication.
///
/// Montgomery multiplication replaces the expensive `u128 % p` operation
/// with a multiplication and bit-shift, which is significantly faster on
/// modern CPUs.
///
/// All values in the NTT are kept in Montgomery form ($x \cdot R \bmod p$,
/// where $R = 2^{64}$). Conversions happen once at entry and exit.
#[derive(Clone, Debug)]
pub struct MontgomeryContext {
    /// The prime modulus.
    p: u64,
    /// $p' = -p^{-1} \bmod 2^{64}$, used for reduction.
    p_inv: u64,
    /// $R^2 \bmod p = 2^{128} \bmod p$, used to convert to Montgomery form.
    r2: u64,
}

impl MontgomeryContext {
    /// Create a new Montgomery context for the given prime.
    ///
    /// Panics if `p` is 0 or even.
    pub fn new(p: u64) -> Self {
        assert!(p > 0 && p % 2 == 1, "p must be a positive odd integer");
        // Compute p_inv = -p^{-1} mod 2^64 via Newton's method.
        // Start with x0 = 1 (since p is odd, p * 1 ≡ 1 mod 2).
        let mut p_inv: u64 = 1;
        for _ in 0..6 {
            p_inv = p_inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(p_inv)));
        }
        // p_inv * p ≡ 1 (mod 2^64), but we want -p^{-1}:
        p_inv = p_inv.wrapping_neg();

        // Compute R^2 mod p = 2^128 mod p via repeated doubling.
        let mut r2: u128 = 1;
        for _ in 0..128 {
            r2 <<= 1;
            if r2 >= p as u128 {
                r2 -= p as u128;
            }
        }

        Self {
            p,
            p_inv,
            r2: r2 as u64,
        }
    }

    /// The prime modulus.
    #[inline]
    pub fn p(&self) -> u64 {
        self.p
    }

    /// Convert a value to Montgomery form: `a * R mod p`.
    ///
    /// Input must be in `[0, p)`.
    #[inline]
    pub fn to_montgomery(&self, a: u64) -> u64 {
        self.mul(a, self.r2)
    }

    /// Convert from Montgomery form back to standard: `a * R^{-1} mod p`.
    ///
    /// Input must be in Montgomery form.
    #[inline]
    pub fn from_montgomery(&self, a: u64) -> u64 {
        self.mont_reduce(a as u128)
    }

    /// Montgomery multiplication: `(a * b * R^{-1}) mod p`.
    ///
    /// Both `a` and `b` must be in Montgomery form. The result is also
    /// in Montgomery form.
    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let prod = a as u128 * b as u128;
        self.mont_reduce(prod)
    }

    /// Montgomery reduction: `t * R^{-1} mod p` where `t < p * R`.
    ///
    /// This is the core operation — replaces `u128 % p` with a
    /// multiplication and shift.
    #[inline]
    fn mont_reduce(&self, t: u128) -> u64 {
        let t_lo = t as u64;
        // m = t_lo * p_inv mod 2^64
        let m = t_lo.wrapping_mul(self.p_inv);
        // Compute (t + m * p) >> 64 in full 128-bit arithmetic.
        let mp = m as u128 * self.p as u128;
        let sum = t.wrapping_add(mp);
        // sum's low 64 bits are zero by construction
        let r = (sum >> 64) as u64;
        // r may be >= p, so conditional subtract
        if r >= self.p { r - self.p } else { r }
    }
}

// ---------------------------------------------------------------------------
// Bit-reverse permutation
// ---------------------------------------------------------------------------

/// Bit-reverse permutation of a slice (in-place).
fn bit_reverse_permute<T>(a: &mut [T]) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    let bits = (n as u64).trailing_zeros() as usize;
    for i in 0..n {
        let j = i.reverse_bits() >> (usize::BITS as usize - bits);
        if i < j {
            a.swap(i, j);
        }
    }
}

/// In-place forward NTT (Cooley-Tukey radix-2 DIT).
///
/// `a` must have a length that is a power of 2. `root` is a primitive
/// `n`-th root of unity in ℤ_p.
pub fn ntt_forward(a: &mut [u64], root: u64, p: u64) {
    let n = a.len();
    debug_assert!(n.is_power_of_two(), "NTT length must be a power of 2");

    bit_reverse_permute(a);

    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let w_step = modpow(root, (n / len) as u64, p);

        for start in (0..n).step_by(len) {
            let mut w: u64 = 1;
            for j in 0..half {
                let u = a[start + j];
                let v = modmul(a[start + j + half], w, p);
                a[start + j] = if u + v >= p { u + v - p } else { u + v };
                a[start + j + half] = if u >= v { u - v } else { u + p - v };
                w = modmul(w, w_step, p);
            }
        }
        len <<= 1;
    }
}

/// In-place inverse NTT.
///
/// `root_inv` is the inverse of the primitive `n`-th root of unity, and
/// `n_inv` is the modular inverse of `n` (i.e. $n^{-1} \bmod p$).
pub fn ntt_inverse(a: &mut [u64], root_inv: u64, p: u64, n_inv: u64) {
    ntt_forward(a, root_inv, p);
    for x in a.iter_mut() {
        *x = modmul(*x, n_inv, p);
    }
}

/// NTT butterfly operations in Montgomery form.
///
/// Assumes all elements of `a` are already in Montgomery form, and
/// `root` is the root of unity in **standard form**. Performs
/// bit-reversal + Cooley-Tukey butterflies. Output stays in
/// Montgomery form.
///
/// Precomputes all twiddle factors once to avoid repeated `modpow`.
fn ntt_butterfly_mont(a: &mut [u64], root: u64, ctx: &MontgomeryContext) {
    let n = a.len();
    let p = ctx.p();

    bit_reverse_permute(a);

    // Precompute w_step for each stage (log2(n) stages).
    let log_n = n.trailing_zeros() as usize;
    let mut stage_roots_m = vec![0u64; log_n];
    for (k, slot) in stage_roots_m.iter_mut().enumerate() {
        let len = 2usize << k;
        let w_step_raw = modpow(root, (n / len) as u64, p);
        *slot = ctx.to_montgomery(w_step_raw);
    }

    for (k, &w_step) in stage_roots_m.iter().enumerate() {
        let len = 2usize << k;
        let half = len / 2;

        for start in (0..n).step_by(len) {
            let mut w: u64 = ctx.to_montgomery(1);
            for j in 0..half {
                let u = a[start + j];
                let v = ctx.mul(a[start + j + half], w);
                a[start + j] = if u + v >= p { u + v - p } else { u + v };
                a[start + j + half] = if u >= v { u - v } else { u + p - v };
                w = ctx.mul(w, w_step);
            }
        }
    }
}

/// NTT-based polynomial multiplication over ℤ_p.
///
/// Given coefficient vectors `a` and `b` (constant term first), returns
/// the coefficient vector of their product. All arithmetic is performed
/// modulo `p`.
///
/// `p` must be a prime such that a suitable root of unity exists for the
/// required transform length (next power of 2 ≥ `a.len() + b.len() - 1`).
///
/// # Panics
///
/// Panics if `p` does not have a suitable root of unity for the required
/// transform length. Callers should check `is_ntt_friendly` first.
pub fn ntt_mul(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let result_len = a.len() + b.len() - 1;
    let n = result_len.next_power_of_two();

    let root = find_primitive_root(p, n)
        .expect("prime must be NTT-friendly for the required transform length");
    let root_inv = modpow(root, p - 2, p);
    let n_inv = modpow(n as u64, p - 2, p);

    let ctx = MontgomeryContext::new(p);

    let mut fa = vec![0u64; n];
    let mut fb = vec![0u64; n];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);

    // Convert to Montgomery form once
    for x in fa.iter_mut() {
        *x = ctx.to_montgomery(*x);
    }
    for x in fb.iter_mut() {
        *x = ctx.to_montgomery(*x);
    }

    // Forward NTT staying in Montgomery form
    ntt_butterfly_mont(&mut fa, root, &ctx);
    ntt_butterfly_mont(&mut fb, root, &ctx);

    // Pointwise multiplication (already in Montgomery form)
    for i in 0..n {
        fa[i] = ctx.mul(fa[i], fb[i]);
    }

    // Inverse NTT staying in Montgomery form
    ntt_butterfly_mont(&mut fa, root_inv, &ctx);

    // Scale by n_inv (in Montgomery form)
    let n_inv_m = ctx.to_montgomery(n_inv);
    for x in fa.iter_mut() {
        *x = ctx.mul(*x, n_inv_m);
    }

    // Convert back from Montgomery form
    for x in fa.iter_mut() {
        *x = ctx.from_montgomery(*x);
    }

    fa.truncate(result_len);
    fa
}

/// Attempt NTT-based multiplication for `DenseUnivariatePolynomial<FiniteField>`.
///
/// Returns `Some(product_coeffs)` if the prime is NTT-friendly and the
/// degree is above the threshold, otherwise `None` (caller should fall
/// back to Karatsuba/Schoolbook).
pub fn try_ntt_mul_fp(
    a_coeffs: &[BigInt],
    b_coeffs: &[BigInt],
    prime: &BigInt,
) -> Option<Vec<BigInt>> {
    // Only use NTT for sufficiently large polynomials
    if a_coeffs.len().min(b_coeffs.len()) < NTT_THRESHOLD {
        return None;
    }

    // Prime must fit in u64 for our NTT implementation
    let p = prime.to_u64_digits().1;
    if p.len() != 1 {
        return None; // prime too large for u64
    }
    let p64 = p[0];

    // Check if a suitable root of unity exists
    let result_len = a_coeffs.len() + b_coeffs.len() - 1;
    let n = result_len.next_power_of_two();
    if !is_ntt_friendly(p64, n) {
        return None;
    }

    // Convert BigInt coefficients to u64
    let a_u64: Vec<u64> = a_coeffs
        .iter()
        .map(|c| {
            let (_, digits) = c.to_u64_digits();
            if digits.is_empty() {
                0
            } else {
                digits[0] % p64
            }
        })
        .collect();
    let b_u64: Vec<u64> = b_coeffs
        .iter()
        .map(|c| {
            let (_, digits) = c.to_u64_digits();
            if digits.is_empty() {
                0
            } else {
                digits[0] % p64
            }
        })
        .collect();

    let result = ntt_mul(&a_u64, &b_u64, p64);

    // Convert back to BigInt
    Some(result.into_iter().map(BigInt::from).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modpow_basic() {
        assert_eq!(modpow(2, 10, 1000), 1024 % 1000);
        assert_eq!(modpow(3, 0, 7), 1);
        assert_eq!(modpow(3, 6, 7), 1); // Fermat: 3^6 ≡ 1 (mod 7)
    }

    #[test]
    fn ntt_friendly_primes() {
        // p = 998244353 = 119 * 2^23 + 1 — common NTT prime
        assert!(is_ntt_friendly(998244353, 1 << 23));
        assert!(!is_ntt_friendly(998244353, 1 << 24));
        // p = 7 — 7 - 1 = 6 = 2 * 3, so NTT works for n = 2
        assert!(is_ntt_friendly(7, 2));
        assert!(!is_ntt_friendly(7, 4));
    }

    #[test]
    fn primitive_root_p7() {
        // For p = 7, n = 2: need g^3 mod 7 where g is a generator
        // 3^3 = 27 ≡ 6 ≡ -1 (mod 7), so 6 is a primitive 2nd root of unity
        let root = find_primitive_root(7, 2).unwrap();
        assert_eq!(modpow(root, 2, 7), 1);
        assert_ne!(root, 1);
    }

    #[test]
    fn ntt_roundtrip_small() {
        // p = 998244353 has primitive 8-th root of unity
        let p: u64 = 998244353;
        let n = 8;
        let root = find_primitive_root(p, n).unwrap();
        let root_inv = modpow(root, p - 2, p);
        let n_inv = modpow(n as u64, p - 2, p);

        let original = vec![1u64, 2, 3, 4, 5, 6, 7, 8];
        let mut data = original.clone();

        ntt_forward(&mut data, root, p);
        ntt_inverse(&mut data, root_inv, p, n_inv);

        assert_eq!(data, original);
    }

    #[test]
    fn montgomery_roundtrip() {
        let p: u64 = 998244353;
        let ctx = MontgomeryContext::new(p);
        for x in [0u64, 1, 42, p / 2, p - 1] {
            let m = ctx.to_montgomery(x);
            let r = ctx.from_montgomery(m);
            assert_eq!(r, x, "roundtrip failed for {x}");
        }
    }

    #[test]
    fn montgomery_mul() {
        let p: u64 = 998244353;
        let ctx = MontgomeryContext::new(p);
        // 3 * 4 = 12 mod p
        let a_m = ctx.to_montgomery(3);
        let b_m = ctx.to_montgomery(4);
        let c_m = ctx.mul(a_m, b_m);
        let c = ctx.from_montgomery(c_m);
        assert_eq!(c, 12);
    }

    #[test]
    fn ntt_mul_trivial() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let p: u64 = 998244353;
        let a = vec![1u64, 2];
        let b = vec![3u64, 4];
        let result = ntt_mul(&a, &b, p);
        assert_eq!(result, vec![3, 10, 8]);
    }

    #[test]
    fn ntt_mul_schoolbook_cross_check() {
        // Cross-check NTT against schoolbook for random-looking polynomials
        let p: u64 = 998244353;
        let a: Vec<u64> = (0..32).map(|i| (i * 7 + 13) % p).collect();
        let b: Vec<u64> = (0..24).map(|i| (i * 11 + 5) % p).collect();

        // Schoolbook
        let result_len = a.len() + b.len() - 1;
        let mut schoolbook = vec![0u64; result_len];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                schoolbook[i + j] = (schoolbook[i + j] + modmul(ai, bj, p)) % p;
            }
        }

        let ntt_result = ntt_mul(&a, &b, p);
        assert_eq!(ntt_result, schoolbook);
    }

    #[test]
    fn ntt_mul_large() {
        // Larger polynomials to exercise the bit-reversal and butterfly stages
        let p: u64 = 998244353;
        let a: Vec<u64> = (0..128).map(|i| (i as u64 * 31 + 17) % p).collect();
        let b: Vec<u64> = (0..128).map(|i| (i as u64 * 47 + 23) % p).collect();

        // Schoolbook reference
        let result_len = a.len() + b.len() - 1;
        let mut schoolbook = vec![0u64; result_len];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                schoolbook[i + j] = (schoolbook[i + j] + modmul(ai, bj, p)) % p;
            }
        }

        let ntt_result = ntt_mul(&a, &b, p);
        assert_eq!(ntt_result.len(), schoolbook.len());
        assert_eq!(ntt_result, schoolbook);
    }

    #[test]
    fn try_ntt_mul_fp_small_returns_none() {
        // Polynomials below threshold should return None
        let a = vec![BigInt::from(1), BigInt::from(2)];
        let b = vec![BigInt::from(3), BigInt::from(4)];
        let prime = BigInt::from(998244353u64);
        assert!(try_ntt_mul_fp(&a, &b, &prime).is_none());
    }

    #[test]
    fn try_ntt_mul_fp_unfriendly_prime_returns_none() {
        // p = 1000000007 — (p-1) = 2 * 500000003, only 2^1 factor
        let n = NTT_THRESHOLD + 1;
        let a: Vec<BigInt> = (0..n).map(BigInt::from).collect();
        let b: Vec<BigInt> = (0..n).map(|i| BigInt::from(i % 10)).collect();
        let prime = BigInt::from(1_000_000_007u64);
        // Transform length would be next_power_of_two(2n-1) which needs > 2^1
        assert!(try_ntt_mul_fp(&a, &b, &prime).is_none());
    }
}
