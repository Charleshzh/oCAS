//! Modular (Brown) GCD for dense univariate polynomials over ℤ.
//!
//! The naive pseudo-remainder GCD in [`crate::gcd`] explodes coefficients
//! for degrees ≳ 16. Brown's algorithm instead computes monic GCDs modulo
//! several primes, reconstructs an integer multiple of the true primitive
//! GCD by CRT with symmetric representatives, and confirms it by exact
//! trial division. Primes where the modular GCD has a larger degree than
//! the true GCD ("unlucky" primes) are detected by degree comparison and
//! discarded.

use ocas_domain::number_theory::{crt::crt_many, primes_from, symmetric_mod};
use ocas_domain::{Domain, EuclideanDomain, FiniteField, Integer, IntegerDomain};
// The GMP backend's `Integer` is not `Sync`: the parallel batch loop is
// compiled out there, so the prelude (par_iter/into_par_iter) is unused.
#[cfg(not(feature = "gmp"))]
use rayon::prelude::*;

use crate::dense::DenseUnivariatePolynomial;
use crate::factor::finite_field::FpPoly;

/// Dense univariate polynomial over ℤ.
pub type ZPoly = DenseUnivariatePolynomial<IntegerDomain>;

/// Safety cap on the number of primes tried before falling back to the
/// pseudo-remainder GCD. In practice CRT succeeds within a few dozen
/// primes (coefficient bit-length divided by ~30).
const MAX_PRIMES: usize = 10_000;

/// Reduce a ℤ[x] polynomial modulo the prime of `field`.
fn reduce_mod_field(p: &ZPoly, field: &FiniteField) -> FpPoly {
    let coeffs = p
        .coeffs()
        .iter()
        .map(|c| field.element(c.to_bigint()))
        .collect();
    FpPoly::from_coeffs(field.clone(), coeffs)
}

/// Exact quotient `dividend / divisor` in ℤ[x], or `None` when the division
/// is not exact (some leading coefficient fails to divide, or a nonzero
/// remainder survives).
fn div_exact_z(dividend: &ZPoly, divisor: &ZPoly) -> Option<ZPoly> {
    if divisor.is_zero() {
        return None;
    }
    let dom = IntegerDomain;
    if dividend.is_zero() {
        return Some(ZPoly::new(dom));
    }
    let div_deg = divisor.degree()?;
    let div_lc = divisor.leading_coeff()?.clone();
    let mut remainder = dividend.clone();
    let mut qcoeffs: Vec<Integer> = Vec::new();
    while let Some(deg) = remainder.degree() {
        if deg < div_deg {
            break;
        }
        let lc = remainder.leading_coeff()?.clone();
        // Exact divisibility check on the leading coefficient.
        let q = dom.div(&lc, &div_lc)?;
        let t = deg - div_deg;
        if qcoeffs.len() <= t {
            qcoeffs.resize(t + 1, Integer::from(0));
        }
        qcoeffs[t] = q.clone();
        let mut sub_coeffs = vec![Integer::from(0); t];
        sub_coeffs.extend(divisor.coeffs().iter().map(|c| &q * c));
        remainder = remainder.sub(&ZPoly::from_coeffs(dom, sub_coeffs));
    }
    if remainder.is_zero() {
        Some(ZPoly::from_coeffs(dom, qcoeffs))
    } else {
        None
    }
}

/// CRT-reconstruct the coefficient vector of the scaled GCD from the
/// per-prime images, using symmetric representatives.
fn reconstruct(images: &[(Integer, FpPoly)], deg: usize) -> Option<ZPoly> {
    let mut coeffs = Vec::with_capacity(deg + 1);
    for i in 0..=deg {
        let cs: Vec<(Integer, Integer)> = images
            .iter()
            .map(|(p, g)| {
                let c = g
                    .coeff(i)
                    .map(|e| Integer::from(e.value().clone()))
                    .unwrap_or_else(|| Integer::from(0));
                (c, p.clone())
            })
            .collect();
        let (r, m) = crt_many(&cs)?;
        coeffs.push(symmetric_mod(&r, &m));
    }
    Some(ZPoly::from_coeffs(IntegerDomain, coeffs))
}

/// Compute one monic, γ-scaled modular GCD image at prime `p`.
///
/// Returns `None` when `p` divides γ or an input vanishes mod `p`
/// (unlucky prime). Extracted so the batch loop can run under either a
/// parallel or sequential iterator without duplicating the body.
fn modular_gcd_image(
    p: &Integer,
    gamma: &Integer,
    ap: &ZPoly,
    bp: &ZPoly,
) -> Option<(Integer, FpPoly, usize)> {
    if gamma.mod_floor(p).is_zero() {
        return None;
    }
    let field = FiniteField::new(p.to_bigint());
    let fa = reduce_mod_field(ap, &field);
    let fb = reduce_mod_field(bp, &field);
    let g = fa.gcd(&fb);
    let deg = g.degree()?; // one input vanished mod p: unlucky
    // Normalize: monic, then scaled by γ.
    let lc = g.leading_coeff()?.clone();
    let inv_lc = field.inv(&lc)?;
    let gamma_p = field.element(gamma.to_bigint());
    let scale = field.mul(&inv_lc, &gamma_p);
    Some((p.clone(), g.mul_scalar(&scale), deg))
}

/// Compute the primitive GCD of `a` and `b` in ℤ[x] by the modular Brown
/// algorithm: monic GCDs over `𝔽_p` are scaled by `γ = gcd(lc a, lc b)`,
/// combined across primes with CRT, and confirmed by exact trial division.
///
/// The result is primitive (like [`DenseUnivariatePolynomial::gcd`]); the
/// contents of the inputs are ignored. Falls back to the pseudo-remainder
/// GCD only if an implausible number of primes was exhausted.
///
/// # Example
///
/// ```
/// use ocas_domain::{IntegerDomain, Integer};
/// use ocas_poly::DenseUnivariatePolynomial;
/// use ocas_poly::gcd::modular::gcd_modular_z;
///
/// let d = IntegerDomain;
/// let i = |v: i64| Integer::from(v);
/// let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
/// let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);
/// let g = gcd_modular_z(&a, &b);
/// assert_eq!(g.coeffs(), &[i(1), i(1)]); // x + 1
/// ```
pub fn gcd_modular_z(a: &ZPoly, b: &ZPoly) -> ZPoly {
    let dom = IntegerDomain;
    if a.is_zero() {
        return b.primitive_part();
    }
    if b.is_zero() {
        return a.primitive_part();
    }
    let ap = a.primitive_part();
    let bp = b.primitive_part();
    // γ is a multiple of the true GCD's leading coefficient; scaling the
    // monic modular images by γ keeps the CRT targets integral.
    let gamma = dom.gcd(
        ap.leading_coeff().expect("nonzero polynomial"),
        bp.leading_coeff().expect("nonzero polynomial"),
    );

    let mut best_deg: Option<usize> = None;
    let mut images: Vec<(Integer, FpPoly)> = Vec::new();
    let mut prime_iter = primes_from(&Integer::from(1_073_741_824)); // > 2^30
    let batch_size = rayon::current_num_threads().max(1);
    let mut tried = 0usize;
    while tried < MAX_PRIMES {
        let batch: Vec<Integer> = prime_iter.by_ref().take(batch_size).collect();
        if batch.is_empty() {
            break;
        }
        tried += batch.len();
        // Monic scaled modular GCD images, computed in parallel (sequential
        // under the GMP backend, whose `Integer` is not `Sync`). Primes
        // dividing γ or where an input vanishes mod p contribute None.
        #[cfg(not(feature = "gmp"))]
        let computed: Vec<Option<(Integer, FpPoly, usize)>> = batch
            .par_iter()
            .map(|p| modular_gcd_image(p, &gamma, &ap, &bp))
            .collect();
        #[cfg(feature = "gmp")]
        let computed: Vec<Option<(Integer, FpPoly, usize)>> = batch
            .iter()
            .map(|p| modular_gcd_image(p, &gamma, &ap, &bp))
            .collect();
        for item in computed {
            let Some((p, g_scaled, deg)) = item else {
                continue;
            };
            match best_deg {
                None => {
                    best_deg = Some(deg);
                    images.push((p, g_scaled));
                }
                Some(bd) if deg < bd => {
                    // Earlier primes were unlucky; restart with the smaller GCD.
                    best_deg = Some(deg);
                    images.clear();
                    images.push((p, g_scaled));
                }
                Some(bd) if deg == bd => images.push((p, g_scaled)),
                _ => continue, // unlucky prime: modular GCD degree too large
            }
            let deg = best_deg.expect("set above");
            if deg == 0 {
                // GCD of the primitive parts is a constant.
                return ZPoly::from_coeffs(dom, vec![Integer::from(1)]);
            }
            // Trial reconstruction: accept only a common divisor of full degree.
            if let Some(candidate) = reconstruct(&images, deg) {
                let cand = candidate.primitive_part();
                if cand.degree() == Some(deg)
                    && div_exact_z(&ap, &cand).is_some()
                    && div_exact_z(&bp, &cand).is_some()
                {
                    return cand;
                }
            }
        }
    }
    // Unreachable in practice; keeps the function total.
    a.gcd(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn i(v: i64) -> Integer {
        Integer::from(v)
    }

    fn zpoly(coeffs: &[i64]) -> ZPoly {
        ZPoly::from_coeffs(IntegerDomain, coeffs.iter().map(|&v| i(v)).collect())
    }

    /// Deterministic pseudo-random polynomial with `deg` coefficients in
    /// `(-bound, bound)`.
    fn rand_poly(deg: usize, bound: i64, seed: &mut u64) -> ZPoly {
        let mut coeffs = Vec::with_capacity(deg + 1);
        for _ in 0..=deg {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = ((*seed >> 33) as i64) % (2 * bound) - bound;
            coeffs.push(i(v));
        }
        // Ensure exact degree.
        if coeffs[deg].is_zero() {
            coeffs[deg] = i(1);
        }
        ZPoly::from_coeffs(IntegerDomain, coeffs)
    }

    /// Normalize sign: make the leading coefficient positive.
    fn monic_sign(p: &ZPoly) -> ZPoly {
        if p.lcoeff().is_negative() {
            p.neg()
        } else {
            p.clone()
        }
    }

    #[test]
    fn small_cases_match_prs() {
        let a = zpoly(&[-1, 0, 1]);
        let b = zpoly(&[1, 2, 1]);
        let g = gcd_modular_z(&a, &b);
        assert_eq!(g.coeffs(), &[i(1), i(1)]);

        let g2 = gcd_modular_z(&zpoly(&[-1, 0, 1]), &zpoly(&[1, 1]));
        assert_eq!(g2.coeffs(), &[i(1), i(1)]);

        // Coprime.
        let g3 = gcd_modular_z(&zpoly(&[1, 1]), &zpoly(&[2, 1]));
        assert_eq!(g3.degree(), Some(0));

        // Zero handling.
        let g4 = gcd_modular_z(&a, &a.zero());
        assert_eq!(g4.coeffs(), a.primitive_part().coeffs());

        // Contents are ignored (primitive result).
        let g5 = gcd_modular_z(&zpoly(&[2, 2]), &zpoly(&[4, 4]));
        assert_eq!(g5.coeffs(), &[i(1), i(1)]);
    }

    #[test]
    fn constructed_common_factor() {
        // g = (3x + 5)(x² + 2) = 3x³ + 5x² + 6x + 10.
        let g = zpoly(&[10, 6, 5, 3]);
        let a = g.mul(&zpoly(&[-1, 2])); // (2x − 1)
        let b = g.mul(&zpoly(&[7, 1])); // (x + 7)
        let got = gcd_modular_z(&a, &b);
        assert_eq!(got, monic_sign(&g).primitive_part());
    }

    #[test]
    fn gcd_is_common_divisor_and_primitive() {
        let mut seed = 42u64;
        for _ in 0..20 {
            let g = rand_poly(4, 10, &mut seed).primitive_part();
            let a = g.mul(&rand_poly(3, 10, &mut seed));
            let b = g.mul(&rand_poly(5, 10, &mut seed));
            let got = gcd_modular_z(&a, &b);
            assert!(div_exact_z(&a, &got).is_some(), "gcd must divide a");
            assert!(div_exact_z(&b, &got).is_some(), "gcd must divide b");
            assert!(got.content().is_one(), "gcd must be primitive");
            assert!(div_exact_z(&got, &g).is_some(), "gcd must contain g");
        }
    }

    #[test]
    fn consistency_with_prs_on_small_polys() {
        let mut seed = 7u64;
        for _ in 0..30 {
            let a = rand_poly(6, 8, &mut seed);
            let b = rand_poly(5, 8, &mut seed);
            let got = monic_sign(&gcd_modular_z(&a, &b));
            let want = monic_sign(&a.gcd(&b));
            assert_eq!(got, want, "a={a:?} b={b:?}");
        }
    }

    #[test]
    fn big_coefficients_no_explosion() {
        // Degree 24 with ~50-digit coefficients: hopeless for naive PRS,
        // routine for the modular path.
        let mut seed = 99u64;
        let mut big = rand_poly(12, 1_000_000, &mut seed);
        // Square it twice to get large coefficients.
        big = big.mul(&big);
        let a = big.mul(&rand_poly(6, 100, &mut seed));
        let b = big.mul(&rand_poly(8, 100, &mut seed));
        let got = gcd_modular_z(&a, &b);
        assert!(div_exact_z(&a, &got).is_some());
        assert!(div_exact_z(&b, &got).is_some());
        assert_eq!(got.degree(), big.primitive_part().degree());
    }

    /// 0.21.0 acceptance: degree-50 polynomials with ~100-digit integer
    /// coefficients — the naive pseudo-remainder GCD explodes on these.
    #[test]
    #[ignore = "performance acceptance: run with --release --ignored"]
    fn modular_gcd_degree_50_100_digit_coeffs() {
        fn big_rand(digits: usize, seed: &mut u64) -> Integer {
            let mut s = String::from("9");
            for _ in 0..digits {
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                s.push((b'0' + ((*seed >> 33) % 10) as u8) as char);
            }
            Integer::from(s.parse::<num_bigint::BigInt>().unwrap())
        }
        let mut seed = 12345u64;
        let mk = |deg: usize, seed: &mut u64| {
            let coeffs: Vec<Integer> = (0..=deg).map(|_| big_rand(50, seed)).collect();
            ZPoly::from_coeffs(IntegerDomain, coeffs)
        };
        let g = mk(25, &mut seed).primitive_part();
        let r1 = mk(25, &mut seed);
        let r2 = mk(25, &mut seed);
        let a = g.mul(&r1);
        let b = g.mul(&r2);
        assert_eq!(a.degree(), Some(50));
        let start = std::time::Instant::now();
        let got = gcd_modular_z(&a, &b);
        let elapsed = start.elapsed();
        eprintln!("deg-50 / 100-digit modular gcd took {elapsed:?}");
        assert!(div_exact_z(&a, &got).is_some(), "gcd must divide a");
        assert!(div_exact_z(&b, &got).is_some(), "gcd must divide b");
        assert_eq!(got.degree(), g.degree());
        assert!(
            elapsed.as_secs() < 120,
            "modular gcd took {elapsed:?} (soft limit 120s)"
        );
    }
}
