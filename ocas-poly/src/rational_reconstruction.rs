//! Rational reconstruction algorithms.
//!
//! Given `a` and modulus `m`, find `(n, d)` such that
//! `a * d ≡ n (mod m)` with `gcd(n, d) = 1` and `2 * |n| * |d| < m`.
//!
//! Uses the extended Euclidean algorithm (continued fraction approach).

use ocas_domain::{Domain, EuclideanDomain, Integer, IntegerDomain};

/// Rational reconstruction using the extended Euclidean algorithm.
///
/// Given `a` and `m`, find `(n, d)` such that:
/// - `a * d ≡ n (mod m)`
/// - `gcd(n, d) = 1`
/// - `2 * |n| * |d| < m` (half-open condition)
///
/// Returns `None` if no such pair exists.
///
/// # Example
///
/// ```
/// use ocas_domain::{Integer, IntegerDomain};
/// use ocas_poly::rational_reconstruction::rational_reconstruction;
///
/// // 3/7 mod 100: 7^{-1} mod 100 = 43, so a = 3 * 43 mod 100 = 29
/// let m = Integer::from(100);
/// let a = Integer::from(29);
/// let result = rational_reconstruction(&a, &m);
/// assert!(result.is_some());
/// ```
pub fn rational_reconstruction(a: &Integer, m: &Integer) -> Option<(Integer, Integer)> {
    let d = IntegerDomain;

    if d.is_zero(m) {
        return None;
    }

    // Reduce a modulo m first: a_red = a mod m, ensuring 0 <= a_red < m.
    let (_, a_red_raw) = d.div_rem(a, m)?;
    let a_red = if is_negative(&a_red_raw) {
        d.add(&a_red_raw, m)
    } else {
        a_red_raw
    };

    if d.is_zero(&a_red) {
        return Some((d.zero(), d.one()));
    }

    // Run extended Euclidean on (m, a_red).
    // Track (r0, r1) and (t0, t1) where at each step:
    //   r_i = s_i * m + t_i * a_red
    let mut r0 = m.clone();
    let mut r1 = a_red.clone();
    let mut t0 = d.zero();
    let mut t1 = d.one();

    // Bound: we want |n| <= bound and |d| <= bound where bound = sqrt(m/2).
    let m_half = d.div(m, &Integer::from(2)).unwrap_or_else(|| m.clone());
    let bound = integer_sqrt(&m_half);

    loop {
        if d.is_zero(&r1) {
            return None;
        }

        // Check convergence: both |r1| and |t1| must be within bound.
        if abs_le(&r1, &bound) && abs_le(&t1, &bound) {
            // Found valid reconstruction: n = r1, d = t1.
            // Verify: a_red * d ≡ n (mod m).
            let check = d.sub(&d.mul(&a_red, &t1), &r1);
            let (_, rem) = d.div_rem(&check, m)?;
            if d.is_zero(&rem) && !d.is_zero(&t1) {
                // Normalize: make d positive.
                let (n, d_val) = if is_negative(&t1) {
                    (d.neg(&r1), d.neg(&t1))
                } else {
                    (r1, t1)
                };
                return Some((n, d_val));
            }
            return None;
        }

        let (q, r_new) = d.div_rem(&r0, &r1)?;
        let t_new = d.sub(&t0, &d.mul(&q, &t1));

        r0 = r1;
        r1 = r_new;
        t0 = t1;
        t1 = t_new;
    }
}

/// Check if `a < 0`.
fn is_negative(a: &Integer) -> bool {
    a < &Integer::from(0)
}

/// Check if `|a| <= b` (absolute value comparison).
fn abs_le(a: &Integer, b: &Integer) -> bool {
    let d = IntegerDomain;
    let abs_a = if is_negative(a) { d.neg(a) } else { a.clone() };
    &abs_a <= b
}

/// Integer square root using Newton's method.
///
/// Returns `floor(sqrt(n))` for `n >= 0`.
fn integer_sqrt(n: &Integer) -> Integer {
    let d = IntegerDomain;
    if d.is_zero(n) {
        return d.zero();
    }
    if n <= &Integer::from(1) {
        return n.clone();
    }

    let two = Integer::from(2);
    let mut x = n.clone();
    loop {
        // next = (x + n/x) / 2
        let n_div_x = d.div(n, &x).unwrap_or(d.zero());
        let sum = d.add(&x, &n_div_x);
        let next = d.div(&sum, &two).unwrap_or(x.clone());

        // Converged when next >= x (approaching from above).
        if next >= x {
            break;
        }
        x = next;
    }
    // Correction: Newton's method with integer division may overshoot.
    // If x^2 > n, decrement until x^2 <= n.
    let one = Integer::from(1);
    while d.mul(&x, &x) > *n && x > one {
        x = d.sub(&x, &one);
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{Domain, EuclideanDomain, Integer, IntegerDomain};

    fn int(i: i64) -> Integer {
        Integer::from(i)
    }

    #[test]
    fn rational_reconstruction_basic() {
        // n=3, d=7, m=101 (prime)
        // 7^{-1} mod 101 = 29 (since 7*29 = 203 ≡ 1 mod 101)
        // a = 3 * 29 mod 101 = 87 mod 101 = 87
        let m = int(101);
        let d_inv_7 = int(29); // 7 * 29 = 203 ≡ 1 (mod 101)
        let a = IntegerDomain
            .div_rem(&IntegerDomain.mul(&int(3), &d_inv_7), &m)
            .unwrap()
            .1;
        let result = rational_reconstruction(&a, &m);
        assert!(result.is_some());
        let (n, d) = result.unwrap();
        // Verify: a * d ≡ n (mod m)
        let check = IntegerDomain.sub(&IntegerDomain.mul(&a, &d), &n);
        let (_, rem) = IntegerDomain.div_rem(&check, &m).unwrap();
        assert!(IntegerDomain.is_zero(&rem));
        assert!(!IntegerDomain.is_zero(&d));
    }

    #[test]
    fn rational_reconstruction_zero() {
        let result = rational_reconstruction(&int(0), &int(100));
        assert_eq!(result, Some((int(0), int(1))));
    }

    #[test]
    fn rational_reconstruction_trivial() {
        // a = 5, m = 101. 5 = 5/1, so n=5, d=1
        let result = rational_reconstruction(&int(5), &int(101));
        assert!(result.is_some());
        let (n, d) = result.unwrap();
        assert_eq!(d, int(1));
        assert_eq!(n, int(5));
    }

    #[test]
    fn rational_reconstruction_failure() {
        // gcd(50, 100) = 50 != 1, so reconstruction fails
        let result = rational_reconstruction(&int(50), &int(100));
        assert!(result.is_none());
    }

    #[test]
    fn rational_reconstruction_zero_modulus() {
        let result = rational_reconstruction(&int(5), &int(0));
        assert!(result.is_none());
    }

    #[test]
    fn rational_reconstruction_one_half() {
        // n=1, d=2, m=101
        // 2^{-1} mod 101 = 51 (since 2*51=102 ≡ 1 mod 101)
        // a = 1 * 51 mod 101 = 51
        let m = int(101);
        let a = int(51);
        let result = rational_reconstruction(&a, &m);
        assert!(result.is_some());
        let (n, d) = result.unwrap();
        // Verify
        let check = IntegerDomain.sub(&IntegerDomain.mul(&a, &d), &n);
        let (_, rem) = IntegerDomain.div_rem(&check, &m).unwrap();
        assert!(IntegerDomain.is_zero(&rem));
        assert!(!IntegerDomain.is_zero(&d));
    }

    #[test]
    fn rational_reconstruction_two_thirds() {
        // n=2, d=3, m=101
        // 3^{-1} mod 101 = 34 (since 3*34=102 ≡ 1 mod 101)
        // a = 2 * 34 mod 101 = 68
        let m = int(101);
        let a = int(68);
        let result = rational_reconstruction(&a, &m);
        assert!(result.is_some());
        let (n, d) = result.unwrap();
        // Verify
        let check = IntegerDomain.sub(&IntegerDomain.mul(&a, &d), &n);
        let (_, rem) = IntegerDomain.div_rem(&check, &m).unwrap();
        assert!(IntegerDomain.is_zero(&rem));
        assert!(!IntegerDomain.is_zero(&d));
    }

    #[test]
    fn integer_sqrt_basic() {
        assert_eq!(integer_sqrt(&int(0)), int(0));
        assert_eq!(integer_sqrt(&int(1)), int(1));
        assert_eq!(integer_sqrt(&int(4)), int(2));
        assert_eq!(integer_sqrt(&int(9)), int(3));
        assert_eq!(integer_sqrt(&int(10)), int(3)); // floor
        assert_eq!(integer_sqrt(&int(50)), int(7)); // floor(sqrt(50)) = 7
    }
}
