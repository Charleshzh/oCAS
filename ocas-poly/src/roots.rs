//! Real root isolation and numerical root approximation.
//!
//! Uses Sturm sequences for exact real-root counting and isolation,
//! with bisection for refinement. For high-precision refinement, the
//! optional `mpfr` backend can be used.

use std::fmt::Display;

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};

use ocas_domain::EuclideanDomain;

use crate::dense::DenseUnivariatePolynomial;

/// A real interval known to contain exactly one root.
#[derive(Debug, Clone, PartialEq)]
pub struct RootInterval {
    /// Lower bound of the interval.
    pub low: f64,
    /// Upper bound of the interval.
    pub high: f64,
}

impl<D: EuclideanDomain> DenseUnivariatePolynomial<D>
where
    D::Element: Display,
{
    /// Compute the Sturm sequence for this polynomial.
    ///
    /// The Sturm sequence is: p0 = p, p1 = p', p_{i+1} = -rem(p_{i-1}, p_i).
    /// The number of sign changes at x gives the number of real roots > x.
    pub fn sturm_sequence(&self) -> Vec<Self> {
        let mut seq = Vec::new();
        if self.is_zero() {
            return seq;
        }

        seq.push(self.clone());
        let deriv = self.derivative();
        if deriv.is_zero() {
            return seq;
        }
        seq.push(deriv);

        loop {
            let a = &seq[seq.len() - 2];
            let b = &seq[seq.len() - 1];
            if b.is_zero() {
                break;
            }
            // Compute pseudo-remainder and negate.
            let rem = match a.pseudo_remainder(b) {
                Some(r) => r,
                None => break,
            };
            if rem.is_zero() {
                break;
            }
            // Negate: p_{i+1} = -rem(p_{i-1}, p_i)
            seq.push(rem.neg());
        }

        seq
    }

    /// Evaluate this polynomial at `x` as a floating-point value.
    ///
    /// Uses Horner's method with f64 arithmetic. For exact evaluation,
    /// use `eval()` with domain elements.
    pub fn eval_f64(&self, x: f64) -> f64 {
        let mut result = 0.0;
        for coeff in self.coeffs().iter().rev() {
            result = result * x + coeff_value(coeff);
        }
        result
    }

    /// Count the number of distinct real roots of this polynomial.
    ///
    /// Uses Sturm's theorem: count roots in (-∞, +∞).
    pub fn count_real_roots(&self) -> usize {
        let seq = self.sturm_sequence();
        if seq.len() < 2 {
            return 0;
        }
        let neg_inf = count_sign_changes_at_infinity(&seq, true);
        let pos_inf = count_sign_changes_at_infinity(&seq, false);
        neg_inf.saturating_sub(pos_inf)
    }

    /// Isolate real roots: return a list of intervals, each containing
    /// exactly one real root.
    ///
    /// Uses bisection with Sturm-based counting to find intervals.
    /// Sign evaluation is exact (dyadic-rational bisection points with
    /// big-integer arithmetic) whenever the coefficients parse as
    /// rationals; ill-conditioned polynomials like the expanded Wilkinson
    /// `∏(x−k)` lose roots under f64 sign evaluation. Falls back to the
    /// f64 path for coefficients without a rational representation.
    pub fn isolate_real_roots(&self) -> Vec<RootInterval> {
        let seq = self.sturm_sequence();
        if seq.len() < 2 {
            return vec![];
        }

        let total_roots = self.count_real_roots();
        if total_roots == 0 {
            return vec![];
        }

        // Find a bounding interval [-M, M] that contains all real roots.
        let m = root_bound(self);
        if let Some(exact) = ExactSturm::prepare(&seq) {
            return self.isolate_exact(&exact, m, total_roots);
        }
        self.isolate_f64(&seq, m, total_roots)
    }

    /// f64 bisection fallback (pre-0.27 behaviour), for Sturm sequences
    /// whose coefficients are not exactly representable as rationals.
    fn isolate_f64(
        &self,
        seq: &[DenseUnivariatePolynomial<D>],
        m: f64,
        total_roots: usize,
    ) -> Vec<RootInterval> {
        let mut intervals = Vec::new();
        let mut stack = vec![(-m, m)];

        while let Some((lo, hi)) = stack.pop() {
            if intervals.len() >= total_roots {
                break;
            }

            let lo_signs = count_sign_changes(seq, lo);
            let hi_signs = count_sign_changes(seq, hi);
            let count = lo_signs.saturating_sub(hi_signs);

            if count == 0 {
                continue;
            }
            if count == 1 && (hi - lo) < 1e-10 {
                intervals.push(RootInterval { low: lo, high: hi });
                continue;
            }
            if hi - lo < 1e-12 {
                if count == 1 {
                    intervals.push(RootInterval { low: lo, high: hi });
                }
                continue;
            }

            let mid = (lo + hi) / 2.0;
            stack.push((lo, mid));
            stack.push((mid, hi));
        }

        intervals
    }

    /// Exact bisection: endpoints are dyadics `m · 2⁻ᵏ`, sign counts are
    /// computed with big-integer arithmetic (no f64 rounding).
    fn isolate_exact(&self, seq: &ExactSturm, m: f64, total_roots: usize) -> Vec<RootInterval> {
        let bound = BigInt::from(m.ceil() as i64);
        let mut intervals = Vec::new();
        // Stack entries: (lo_m, hi_m, k) representing [lo_m·2⁻ᵏ, hi_m·2⁻ᵏ].
        let mut stack = vec![(-&bound, bound, 0u32)];

        while let Some((lo_m, hi_m, k)) = stack.pop() {
            if intervals.len() >= total_roots {
                break;
            }

            let lo_signs = seq.count_at(&lo_m, k);
            let hi_signs = seq.count_at(&hi_m, k);
            let count = lo_signs.saturating_sub(hi_signs);

            if count == 0 {
                continue;
            }
            // Width as f64 is only used for the termination thresholds,
            // never for sign decisions.
            let width = (&hi_m - &lo_m).to_f64().unwrap_or(f64::INFINITY) / 2.0f64.powi(k as i32);
            let low_f = dyadic_f64(&lo_m, k);
            let high_f = dyadic_f64(&hi_m, k);
            if count == 1 && width < 1e-10 {
                intervals.push(RootInterval {
                    low: low_f,
                    high: high_f,
                });
                continue;
            }
            if width < 1e-12 {
                if count == 1 {
                    intervals.push(RootInterval {
                        low: low_f,
                        high: high_f,
                    });
                }
                continue;
            }

            // mid = (lo + hi)/2; all endpoints move to scale k+1.
            let lo2 = &lo_m << 1usize;
            let hi2 = &hi_m << 1usize;
            let mid = &lo_m + &hi_m;
            stack.push((lo2, mid.clone(), k + 1));
            stack.push((mid, hi2, k + 1));
        }

        intervals
    }

    /// Refine a root interval using bisection to the given tolerance.
    pub fn refine_root(&self, interval: &RootInterval, tol: f64) -> RootInterval {
        let mut lo = interval.low;
        let mut hi = interval.high;
        let f_lo = self.eval_f64(lo);

        if f_lo.abs() < 1e-15 {
            return RootInterval { low: lo, high: lo };
        }

        while hi - lo > tol {
            let mid = (lo + hi) / 2.0;
            let f_mid = self.eval_f64(mid);
            if f_mid.abs() < 1e-15 {
                return RootInterval {
                    low: mid,
                    high: mid,
                };
            }
            if f_lo * f_mid < 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }

        RootInterval { low: lo, high: hi }
    }
}

/// A Sturm sequence pre-scaled to integer coefficients for exact dyadic
/// evaluation (per polynomial, coefficients are `cᵢ·lcm(denoms)`).
struct ExactSturm {
    polys: Vec<Vec<BigInt>>,
}

impl ExactSturm {
    /// Clear denominators of every Sturm polynomial. Returns `None` when a
    /// coefficient has no rational text form.
    fn prepare<D: EuclideanDomain>(seq: &[DenseUnivariatePolynomial<D>]) -> Option<Self>
    where
        D::Element: Display,
    {
        let mut polys = Vec::with_capacity(seq.len());
        for p in seq {
            let mut coeffs = Vec::with_capacity(p.coeffs().len());
            let mut lcm = BigInt::from(1);
            for c in p.coeffs() {
                let r = coeff_to_bigrational(c)?;
                lcm = bigint_lcm(&lcm, r.denom());
                coeffs.push(r);
            }
            let scale = BigRational::from_integer(lcm);
            polys.push(coeffs.iter().map(|c| (c * &scale).to_integer()).collect());
        }
        Some(Self { polys })
    }

    /// Exact sign of every Sturm polynomial at `m·2⁻ᵏ`: for a degree-`n`
    /// polynomial, `p(m·2⁻ᵏ) = 2⁻ᵏⁿ·L⁻¹·Σ aᵢ·mⁱ·2ᵏ⁽ⁿ⁻ⁱ⁾` and the scale
    /// factors are positive, so the sum's sign is the value's sign.
    fn count_at(&self, m: &BigInt, k: u32) -> usize {
        let mut count = 0;
        let mut prev: Option<bool> = None;
        for icoeffs in &self.polys {
            let sign = eval_sign_dyadic(icoeffs, m, k);
            if sign == 0 {
                continue;
            }
            let positive = sign > 0;
            if let Some(p) = prev
                && p != positive
            {
                count += 1;
            }
            prev = Some(positive);
        }
        count
    }
}

/// Sign of `Σ aᵢ·mⁱ·2ᵏ⁽ⁿ⁻ⁱ⁾` for integer coefficients `aᵢ` (ascending).
fn eval_sign_dyadic(coeffs: &[BigInt], m: &BigInt, k: u32) -> i8 {
    let n = coeffs.len().saturating_sub(1);
    let mut s = BigInt::zero();
    for (i, a) in coeffs.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        let mut t = a * m.pow(i as u32);
        t <<= k as usize * (n - i);
        s += t;
    }
    match s.sign() {
        Sign::Plus => 1,
        Sign::Minus => -1,
        Sign::NoSign => 0,
    }
}

/// Convert a dyadic `m·2⁻ᵏ` to f64 (display/intervals only).
fn dyadic_f64(m: &BigInt, k: u32) -> f64 {
    m.to_f64().unwrap_or(f64::NAN) / 2.0f64.powi(k as i32)
}

/// Parse a domain element's display text as an exact rational
/// ("n", "-n", or "n/d").
fn coeff_to_bigrational(elem: &(impl Display + ?Sized)) -> Option<BigRational> {
    elem.to_string().trim().parse::<BigRational>().ok()
}

/// lcm for positive big integers (Euclid).
fn bigint_lcm(a: &BigInt, b: &BigInt) -> BigInt {
    let mut x = a.clone();
    let mut y = b.clone();
    while !y.is_zero() {
        let r = x % &y;
        x = y;
        y = r;
    }
    if x.is_zero() {
        return BigInt::from(1);
    }
    (a / &x) * b
}

/// Count sign changes in the Sturm sequence at ±∞.
fn count_sign_changes_at_infinity<D: EuclideanDomain>(
    seq: &[DenseUnivariatePolynomial<D>],
    at_neg_inf: bool,
) -> usize
where
    D::Element: Display,
{
    let vals: Vec<f64> = seq
        .iter()
        .map(|p| {
            if p.is_zero() {
                return 0.0;
            }
            let deg = p.degree().unwrap_or(0);
            let lc = coeff_value(p.leading_coeff().unwrap());
            // At +∞: sign of leading coefficient
            // At -∞: sign depends on degree parity
            if at_neg_inf {
                if deg % 2 == 0 { lc } else { -lc }
            } else {
                lc
            }
        })
        .collect();
    count_sign_changes_in_vals(&vals)
}

/// Evaluate the Sturm sequence at `x` and count sign changes.
fn count_sign_changes<D: EuclideanDomain>(seq: &[DenseUnivariatePolynomial<D>], x: f64) -> usize
where
    D::Element: Display,
{
    let vals: Vec<f64> = seq.iter().map(|p| p.eval_f64(x)).collect();
    count_sign_changes_in_vals(&vals)
}

fn count_sign_changes_in_vals(vals: &[f64]) -> usize {
    let mut count = 0;
    let mut prev_sign: Option<bool> = None;
    for &v in vals {
        if v == 0.0 {
            continue;
        }
        let sign = v > 0.0;
        if let Some(p) = prev_sign
            && p != sign
        {
            count += 1;
        }
        prev_sign = Some(sign);
    }
    count
}

/// Compute a bound M such that all real roots lie in [-M, M].
fn root_bound<D: EuclideanDomain>(p: &DenseUnivariatePolynomial<D>) -> f64
where
    D::Element: Display,
{
    if p.is_zero() || p.degree().is_none() {
        return 1.0;
    }
    let coeffs = p.coeffs();
    let lc = coeff_value(coeffs.last().unwrap()).abs();
    let mut max_abs = 0.0f64;
    for c in &coeffs[..coeffs.len() - 1] {
        let v = coeff_value(c).abs();
        if v > max_abs {
            max_abs = v;
        }
    }
    1.0 + max_abs / lc.max(1e-10)
}

/// Convert a domain element to f64 for numerical evaluation.
fn coeff_value(elem: &(impl Display + ?Sized)) -> f64 {
    let s = elem.to_string();
    let trimmed = s.trim();
    // Try direct f64 parse first.
    if let Ok(v) = trimmed.parse::<f64>() {
        return v;
    }
    // Try integer parse.
    if let Ok(v) = trimmed.parse::<i64>() {
        return v as f64;
    }
    // Try rational format "n/d".
    if let Some((num_str, den_str)) = trimmed.split_once('/')
        && let (Ok(n), Ok(d)) = (num_str.trim().parse::<f64>(), den_str.trim().parse::<f64>())
        && d != 0.0
    {
        return n / d;
    }
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain};

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn count_roots_x2_minus_1() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
        assert_eq!(p.count_real_roots(), 2);
    }

    #[test]
    fn count_roots_x2_plus_1() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(0), i(1)]);
        assert_eq!(p.count_real_roots(), 0);
    }

    #[test]
    fn count_roots_perfect_square() {
        let d = IntegerDomain;
        // (x+1)^2 = x^2 + 2x + 1 has one distinct real root at x = -1
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);
        assert_eq!(p.count_real_roots(), 1);
    }

    #[test]
    fn isolate_roots_x2_minus_2() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-2), i(0), i(1)]);
        let intervals = p.isolate_real_roots();
        assert_eq!(intervals.len(), 2);
        // Refine one root to verify it's near sqrt(2) ≈ 1.414
        let refined = p.refine_root(&intervals[1], 1e-6);
        let approx = (refined.low + refined.high) / 2.0;
        assert!((approx.abs() - std::f64::consts::SQRT_2).abs() < 0.01);
    }

    #[test]
    fn sturm_sequence_length() {
        let d = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]);
        let seq = p.sturm_sequence();
        assert!(seq.len() >= 2);
    }
}
