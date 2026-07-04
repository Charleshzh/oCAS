//! Dense univariate polynomial implementation.
//!
//! A [`DenseUnivariatePolynomial`] stores all coefficients from the constant
//! term up to the leading coefficient in a contiguous vector. This is well
//! suited for univariate arithmetic with moderate degree.

use ocas_domain::{Domain, EuclideanDomain};

/// Threshold below which Karatsuba falls back to schoolbook multiplication.
const KARATSUBA_THRESHOLD: usize = 32;

/// A dense univariate polynomial with coefficients in a domain `D`.
///
/// # Example
///
/// ```
/// use ocas_domain::{IntegerDomain, Integer};
/// use ocas_poly::DenseUnivariatePolynomial;
///
/// let domain = IntegerDomain;
/// let p = DenseUnivariatePolynomial::from_coeffs(
///     domain,
///     vec![Integer::from(1), Integer::from(2), Integer::from(1)],
/// );
/// let q = DenseUnivariatePolynomial::from_coeffs(
///     domain,
///     vec![Integer::from(1), Integer::from(1)],
/// );
/// let r = p.mul(&q);
/// assert_eq!(r.coeffs(), &[
///     Integer::from(1),
///     Integer::from(3),
///     Integer::from(3),
///     Integer::from(1),
/// ]);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseUnivariatePolynomial<D: Domain> {
    /// Coefficients from constant term upward. Trailing zeros are removed so
    /// the zero polynomial is represented by an empty vector.
    coeffs: Vec<D::Element>,
    /// The coefficient domain. Stored in the polynomial so all operations can
    /// access it without passing it explicitly.
    domain: D,
}

impl<D: Domain> DenseUnivariatePolynomial<D> {
    /// Create the zero polynomial over `domain`.
    pub fn new(domain: D) -> Self {
        Self {
            coeffs: Vec::new(),
            domain,
        }
    }

    /// Create a polynomial from a vector of coefficients `[a0, a1, ..., an]`.
    ///
    /// Trailing zero coefficients are stripped automatically.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let domain = IntegerDomain;
    /// let p = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(0), Integer::from(2)],
    /// );
    /// assert_eq!(p.degree(), Some(2));
    /// assert_eq!(p.coeff(2), Some(&Integer::from(2)));
    /// ```
    pub fn from_coeffs(domain: D, coeffs: Vec<D::Element>) -> Self {
        let mut poly = Self { coeffs, domain };
        poly.trim_trailing_zeros();
        poly
    }

    /// Return a reference to the coefficient domain.
    pub fn domain(&self) -> &D {
        &self.domain
    }

    /// Return the coefficients from constant term upward.
    pub fn coeffs(&self) -> &[D::Element] {
        &self.coeffs
    }

    /// Return whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Return the degree of the polynomial, or `None` for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// Return the coefficient of `x^n`, or `None` if the term is absent.
    pub fn coeff(&self, n: usize) -> Option<&D::Element> {
        self.coeffs.get(n)
    }

    /// Return the leading coefficient, or `None` for the zero polynomial.
    pub fn leading_coeff(&self) -> Option<&D::Element> {
        self.coeffs.last()
    }

    /// Convenience alias: return the leading coefficient, or the domain's
    /// zero element for the zero polynomial.
    pub fn lcoeff(&self) -> D::Element {
        self.leading_coeff()
            .cloned()
            .unwrap_or_else(|| self.domain.zero())
    }

    /// Return the constant term (coefficient of $x^0$), or the domain's
    /// zero element for the zero polynomial.
    pub fn constant(&self) -> D::Element {
        self.coeff(0).cloned().unwrap_or_else(|| self.domain.zero())
    }

    /// Return the zero polynomial with the same domain.
    pub fn zero(&self) -> Self {
        Self::new(self.domain.clone())
    }

    /// Return the constant polynomial `1` over the same domain.
    pub fn one(&self) -> Self {
        Self::from_coeffs(self.domain.clone(), vec![self.domain.one()])
    }

    /// Return whether this is the constant polynomial 1.
    pub fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.domain.is_one(&self.coeffs[0])
    }

    /// Return the negation of this polynomial.
    pub fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter().map(|c| self.domain.neg(c)).collect();
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Add another polynomial.
    pub fn add(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(len);
        let zero = self.domain.zero();
        for i in 0..len {
            let a = self.coeffs.get(i).unwrap_or(&zero);
            let b = other.coeffs.get(i).unwrap_or(&zero);
            coeffs.push(self.domain.add(a, b));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Subtract another polynomial.
    pub fn sub(&self, other: &Self) -> Self {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(len);
        let zero = self.domain.zero();
        for i in 0..len {
            let a = self.coeffs.get(i).unwrap_or(&zero);
            let b = other.coeffs.get(i).unwrap_or(&zero);
            coeffs.push(self.domain.sub(a, b));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Multiply by a scalar coefficient.
    pub fn mul_scalar(&self, scalar: &D::Element) -> Self {
        if self.domain.is_zero(scalar) {
            return self.zero();
        }
        let coeffs = self
            .coeffs
            .iter()
            .map(|c| self.domain.mul(c, scalar))
            .collect();
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Multiply two polynomials.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let domain = IntegerDomain;
    /// let a = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(1)],
    /// );
    /// let b = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(-1)],
    /// );
    /// let c = a.mul(&b);
    /// assert_eq!(c.coeffs(), &[Integer::from(1), Integer::from(0), Integer::from(-1)]);
    /// ```
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return self.zero();
        }
        let mut buf = Vec::new();
        self.mul_into(other, &mut buf);
        Self::from_coeffs(self.domain.clone(), buf)
    }

    /// Multiply two polynomials, reusing the provided buffer for the result.
    ///
    /// The buffer is cleared and resized as needed. This avoids repeated
    /// heap allocation in hot loops (e.g. GCD, factorization).
    ///
    /// After the call, `buf` contains the coefficients of the product
    /// (constant term first). If either polynomial is zero, `buf` is cleared.
    pub fn mul_into(&self, other: &Self, buf: &mut Vec<D::Element>) {
        if self.is_zero() || other.is_zero() {
            buf.clear();
            return;
        }
        if self.coeffs.len().min(other.coeffs.len()) >= KARATSUBA_THRESHOLD {
            Self::karatsuba_mul_into(&self.coeffs, &other.coeffs, &self.domain, buf);
        } else {
            Self::schoolbook_mul_into(&self.coeffs, &other.coeffs, &self.domain, buf);
        }
    }

    /// Schoolbook O(n·m) polynomial multiplication into `buf`.
    fn schoolbook_mul_into(
        a: &[D::Element],
        b: &[D::Element],
        domain: &D,
        buf: &mut Vec<D::Element>,
    ) {
        let result_len = a.len() + b.len() - 1;
        buf.clear();
        buf.resize(result_len, domain.zero());
        for (i, ai) in a.iter().enumerate() {
            for (j, bj) in b.iter().enumerate() {
                let prod = domain.mul(ai, bj);
                buf[i + j] = domain.add(&buf[i + j], &prod);
            }
        }
    }

    /// Karatsuba fast multiplication into `buf`.
    ///
    /// Splits each polynomial at the midpoint `m = n/2` and computes the
    /// product using three half-size multiplications instead of four:
    ///
    /// ```text
    /// a = a0 + a1·x^m,  b = b0 + b1·x^m
    /// z0 = a0·b0,  z2 = a1·b1
    /// z1 = (a0+a1)·(b0+b1) − z0 − z2
    /// result = z0 + z1·x^m + z2·x^(2m)
    /// ```
    fn karatsuba_mul_into(
        a: &[D::Element],
        b: &[D::Element],
        domain: &D,
        buf: &mut Vec<D::Element>,
    ) {
        let n = a.len().max(b.len());
        if n < KARATSUBA_THRESHOLD {
            Self::schoolbook_mul_into(a, b, domain, buf);
            return;
        }

        let m = n / 2;
        let zero = domain.zero();

        // Split: a = a0 + a1·x^m, b = b0 + b1·x^m
        let (a0, a1) = if a.len() <= m {
            (a, &[][..])
        } else {
            (&a[..m], &a[m..])
        };
        let (b0, b1) = if b.len() <= m {
            (b, &[][..])
        } else {
            (&b[..m], &b[m..])
        };

        // z0 = a0 * b0
        let mut z0 = Vec::new();
        Self::karatsuba_mul_into(a0, b0, domain, &mut z0);

        // z2 = a1 * b1
        let mut z2 = Vec::new();
        Self::karatsuba_mul_into(a1, b1, domain, &mut z2);

        // a01 = a0 + a1,  b01 = b0 + b1
        let a01_len = a0.len().max(a1.len());
        let b01_len = b0.len().max(b1.len());
        let mut a01 = vec![zero.clone(); a01_len];
        let mut b01 = vec![zero.clone(); b01_len];
        for (i, a01_val) in a01.iter_mut().enumerate() {
            let ai = a0.get(i).unwrap_or(&zero);
            let aj = a1.get(i).unwrap_or(&zero);
            *a01_val = domain.add(ai, aj);
        }
        for (i, b01_val) in b01.iter_mut().enumerate() {
            let bi = b0.get(i).unwrap_or(&zero);
            let bj = b1.get(i).unwrap_or(&zero);
            *b01_val = domain.add(bi, bj);
        }

        // z1 = (a0+a1)*(b0+b1) - z0 - z2
        let mut z1 = Vec::new();
        Self::karatsuba_mul_into(&a01, &b01, domain, &mut z1);
        for (i, z1_val) in z1.iter_mut().enumerate() {
            let z0i = z0.get(i).unwrap_or(&zero);
            let z2i = z2.get(i).unwrap_or(&zero);
            *z1_val = domain.sub(z1_val, z0i);
            *z1_val = domain.sub(z1_val, z2i);
        }

        // Combine: result = z0 + z1·x^m + z2·x^(2m)
        let result_len = a.len() + b.len() - 1;
        buf.clear();
        buf.resize(result_len, zero);

        for (i, c) in z0.iter().enumerate() {
            buf[i] = domain.add(&buf[i], c);
        }
        for (i, c) in z1.iter().enumerate() {
            let idx = i + m;
            if idx < buf.len() {
                buf[idx] = domain.add(&buf[idx], c);
            }
        }
        for (i, c) in z2.iter().enumerate() {
            let idx = i + 2 * m;
            if idx < buf.len() {
                buf[idx] = domain.add(&buf[idx], c);
            }
        }
    }

    /// Evaluate the polynomial at `x` using Horner's method.
    ///
    /// The zero polynomial evaluates to the domain's zero element.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let domain = IntegerDomain;
    /// let p = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(2), Integer::from(3)],
    /// );
    /// let value = p.eval(&Integer::from(2));
    /// assert_eq!(value, Integer::from(17));
    /// ```
    pub fn eval(&self, x: &D::Element) -> D::Element {
        let mut result = self.domain.zero();
        for coeff in self.coeffs.iter().rev() {
            result = self.domain.mul(&result, x);
            result = self.domain.add(&result, coeff);
        }
        result
    }

    /// Return the formal derivative of this polynomial.
    ///
    /// For `p(x) = a_0 + a_1 x + a_2 x^2 + ...` the derivative is
    /// `p'(x) = a_1 + 2 a_2 x + 3 a_3 x^2 + ...`.
    pub fn derivative(&self) -> Self {
        if self.degree().is_none() {
            return self.zero();
        }
        let mut coeffs = Vec::with_capacity(self.coeffs.len().saturating_sub(1));
        for (i, c) in self.coeffs.iter().enumerate().skip(1) {
            let scalar = self.domain.cast_u64(i as u64);
            coeffs.push(self.domain.mul(c, &scalar));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }

    /// Return the formal integral of this polynomial, with constant term zero.
    ///
    /// For `p(x) = a_0 + a_1 x + a_2 x^2 + ...` the integral is
    /// `∫p(x) dx = 0 + a_0 x + (a_1/2) x^2 + (a_2/3) x^3 + ...`.
    pub fn integral(&self) -> Self {
        if self.is_zero() {
            return self.zero();
        }
        let mut coeffs = Vec::with_capacity(self.coeffs.len() + 1);
        coeffs.push(self.domain.zero());
        for (i, c) in self.coeffs.iter().enumerate() {
            let denom = self.domain.cast_u64((i + 1) as u64);
            let inv = self
                .domain
                .inv(&denom)
                .unwrap_or_else(|| self.domain.zero());
            coeffs.push(self.domain.mul(c, &inv));
        }
        Self::from_coeffs(self.domain.clone(), coeffs)
    }
}

impl<D: EuclideanDomain> DenseUnivariatePolynomial<D> {
    /// Multiply all coefficients by a constant.
    ///
    /// Equivalent to [`mul_scalar`](Self::mul_scalar) but restricted to
    /// [`EuclideanDomain`] for consistency with [`div_coeff`](Self::div_coeff).
    pub fn mul_coeff(&self, c: &D::Element) -> Self {
        self.mul_scalar(c)
    }

    /// Divide all coefficients by a constant (must divide exactly).
    ///
    /// Panics in debug mode if any coefficient is not divisible by `c`.
    pub fn div_coeff(&self, c: &D::Element) -> Self {
        let inv = self.domain.inv(c).expect("div_coeff: cannot invert zero");
        self.mul_scalar(&inv)
    }

    /// Divide this polynomial by another, returning `(quotient, remainder)`.
    ///
    /// Returns `None` if the divisor is the zero polynomial.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_domain::{IntegerDomain, Integer};
    /// use ocas_poly::DenseUnivariatePolynomial;
    ///
    /// let domain = IntegerDomain;
    /// let p = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(0), Integer::from(-1)],
    /// );
    /// let q = DenseUnivariatePolynomial::from_coeffs(
    ///     domain,
    ///     vec![Integer::from(1), Integer::from(1)],
    /// );
    /// let (quot, rem) = p.div_rem(&q).unwrap();
    /// assert_eq!(quot.coeffs(), &[Integer::from(1), Integer::from(-1)]);
    /// assert!(rem.is_zero());
    /// ```
    pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)> {
        if divisor.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some((self.zero(), self.zero()));
        }
        let mut remainder = self.clone();
        let mut quotient_coeffs: Vec<D::Element> = Vec::new();
        let domain = self.domain.clone();
        let divisor_degree = divisor.degree().unwrap_or(0);
        let divisor_lc = divisor.leading_coeff().unwrap().clone();

        while let Some(deg) = remainder.degree() {
            if deg < divisor_degree {
                break;
            }
            let lc = remainder.leading_coeff().unwrap().clone();
            let (q, _r) = domain.div_rem(&lc, &divisor_lc)?;
            let term_degree = deg - divisor_degree;

            // Ensure quotient_coeffs is long enough.
            if term_degree >= quotient_coeffs.len() {
                quotient_coeffs.resize(term_degree + 1, domain.zero());
            }
            quotient_coeffs[term_degree] = domain.add(&quotient_coeffs[term_degree], &q);

            // remainder -= q * x^term_degree * divisor
            let mut sub_coeffs = vec![domain.zero(); term_degree];
            sub_coeffs.extend(divisor.coeffs.iter().map(|c| domain.mul(c, &q)));
            let sub = Self::from_coeffs(domain.clone(), sub_coeffs);
            remainder = remainder.sub(&sub);

            // Stop if remainder did not shrink (defensive against non-exact division).
            if let Some(rem_deg) = remainder.degree() {
                if rem_deg >= deg {
                    break;
                }
            } else {
                break;
            }
        }

        let quotient = Self::from_coeffs(domain, quotient_coeffs);
        Some((quotient, remainder))
    }

    // ------------------------------------------------------------------
    //  Diophantine / CRT and p-adic expansion (for partial fractions)
    // ------------------------------------------------------------------

    /// Compute `self^n` by repeated squaring.
    pub fn pow(&self, n: u32) -> Self {
        if n == 0 {
            return self.one();
        }
        let mut base = self.clone();
        let mut exp = n;
        let mut result = self.one();
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exp >>= 1;
        }
        result
    }

    /// Compute the extended GCD of two polynomials: `(g, s, t)` such that
    /// `s * self + t * other = g` where `g = gcd(self, other)`.
    ///
    /// Uses the extended Euclidean algorithm.
    pub fn extended_gcd_poly(&self, other: &Self) -> (Self, Self, Self) {
        let d = self.domain();
        let mut old_r = self.clone();
        let mut r = other.clone();
        let mut old_s = self.one();
        let mut s = self.zero();
        let mut old_t = self.zero();
        let mut t = self.one();

        while !r.is_zero() {
            let (q, rem) = old_r.div_rem(&r).unwrap_or((old_r.zero(), old_r.clone()));
            old_r = r;
            r = rem;

            let new_s = old_s.sub(&q.mul(&s));
            old_s = s;
            s = new_s;

            let new_t = old_t.sub(&q.mul(&t));
            old_t = t;
            t = new_t;
        }

        // Normalize so that g is monic (or at least has positive leading coeff).
        if let Some(lc) = old_r.leading_coeff()
            && !d.is_one(lc)
            && let Some(lc_inv) = d.inv(lc)
        {
            old_r = old_r.mul_scalar(&lc_inv);
            old_s = old_s.mul_scalar(&lc_inv);
            old_t = old_t.mul_scalar(&lc_inv);
        }

        (old_r, old_s, old_t)
    }

    /// Polynomial CRT (diophantine solver).
    ///
    /// Given a list of pairwise coprime polynomials `polys` and a target `b`,
    /// returns `[s0, ..., sn]` such that:
    ///
    /// $$\sum_i s_i \cdot \prod_{j \neq i} p_j \equiv b \pmod{\prod_i p_i}$$
    ///
    /// Uses the extended Euclidean algorithm recursively.
    ///
    /// # Panics
    ///
    /// Panics if the polynomials are not pairwise coprime (i.e. the GCD is
    /// not a unit).
    pub fn diophantine(polys: &mut [Self], b: &Self) -> Vec<Self> {
        let n = polys.len();
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            let (_, r) = b
                .div_rem(&polys[0])
                .unwrap_or_else(|| (b.zero(), b.clone()));
            return vec![r];
        }

        // Compute suffix products: suffix[i] = Π_{j>i} polys[j]
        let mut suffix: Vec<Self> = Vec::with_capacity(n);
        let mut prod = polys[n - 1].one(); // empty product = 1
        for i in (0..n - 1).rev() {
            prod = prod.mul(&polys[i + 1]);
            suffix.push(prod.clone());
        }
        suffix.reverse();
        // suffix[i] = Π_{j>i} polys[j] for i in 0..n-1

        // Recursive EEA approach:
        // Start with cur = b, then for each i:
        //   (g, s, t) = extended_gcd(p[i], suffix[i])
        //   result[i] = (t * cur) mod p[i]
        //   cur = (s * cur) mod suffix[i]
        let mut cur = b.clone();
        let mut result = Vec::with_capacity(n);

        for i in 0..n {
            let (g, s, t) = polys[i].extended_gcd_poly(&suffix[i]);

            // g should be a unit (constant, ideally 1).
            // result[i] = (t * cur) / g  mod p[i]
            let ts = t.mul(&cur);
            let ts_div_g = if g.is_one() {
                ts
            } else {
                ts.div_rem(&g).map(|(q, _)| q).unwrap_or(ts)
            };
            let (_, ri) = ts_div_g
                .div_rem(&polys[i])
                .unwrap_or_else(|| (ts_div_g.zero(), ts_div_g));

            // Update: cur = (s * cur) / g  mod suffix[i]
            if i < n - 1 {
                let ss = s.mul(&cur);
                let ss_div_g = if g.is_one() {
                    ss
                } else {
                    ss.div_rem(&g).map(|(q, _)| q).unwrap_or(ss)
                };
                let (_, new_cur) = ss_div_g
                    .div_rem(&suffix[i])
                    .unwrap_or_else(|| (ss_div_g.zero(), ss_div_g));
                cur = new_cur;
            }

            result.push(ri);
        }

        result
    }

    /// p-adic expansion of `self` with respect to `p`.
    ///
    /// Returns `[a0, a1, a2, ...]` such that:
    ///
    /// $$\text{self} = a_0 + a_1 \cdot p + a_2 \cdot p^2 + \cdots$$
    ///
    /// where each $a_k$ has degree less than $\deg(p)$.
    ///
    /// This is computed by repeated polynomial division (like integer
    /// p-adic expansion).
    pub fn p_adic_expansion(&self, p: &Self) -> Vec<Self> {
        let mut result = Vec::new();
        let mut r = self.clone();
        while !r.is_zero() {
            let (q, rem) = match r.div_rem(p) {
                Some(v) => v,
                None => break,
            };
            result.push(rem);
            r = q;
        }
        result
    }
}

impl<D: Domain> DenseUnivariatePolynomial<D> {
    fn trim_trailing_zeros(&mut self) {
        while let Some(last) = self.coeffs.last() {
            if self.domain.is_zero(last) {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_domain::{FiniteField, Integer, IntegerDomain, Rational, RationalDomain};

    fn int(i: i64) -> Integer {
        Integer::from(i)
    }

    #[test]
    fn zero_polynomial_has_no_degree() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::new(domain);
        assert!(p.is_zero());
        assert_eq!(p.degree(), None);
    }

    #[test]
    fn degree_and_coeffs() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![3.into(), 0.into(), 2.into(), 0.into()],
        );
        assert_eq!(p.degree(), Some(2));
        assert_eq!(p.coeff(0).cloned(), Some(3.into()));
        assert_eq!(p.coeff(2).cloned(), Some(2.into()));
        assert_eq!(p.coeff(3), None);
        assert_eq!(p.leading_coeff().cloned(), Some(2.into()));
    }

    #[test]
    fn trailing_zeros_are_trimmed() {
        let domain = IntegerDomain;
        let p = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![1.into(), 2.into(), 0.into(), 0.into()],
        );
        assert_eq!(p.degree(), Some(1));
        assert_eq!(p.coeffs().len(), 2);
    }

    #[test]
    fn add_polynomials() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let sum = a.add(&b);
        assert_eq!(sum.coeffs().to_vec(), vec![4.into(), 2.into(), 4.into()]);
    }

    #[test]
    fn sub_polynomials() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let diff = a.sub(&b);
        assert_eq!(
            diff.coeffs().to_vec(),
            vec![(-2).into(), 2.into(), (-4).into()]
        );
    }

    #[test]
    fn mul_polynomials() {
        let domain = IntegerDomain;
        // (1 + 2x) * (3 + 4x^2) = 3 + 6x + 4x^2 + 8x^3
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let b = DenseUnivariatePolynomial::from_coeffs(domain, vec![3.into(), 0.into(), 4.into()]);
        let prod = a.mul(&b);
        assert_eq!(
            prod.coeffs().to_vec(),
            vec![3.into(), 6.into(), 4.into(), 8.into()]
        );
    }

    #[test]
    fn mul_by_zero_yields_zero() {
        let domain = IntegerDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(domain, vec![1.into(), 2.into()]);
        let zero = DenseUnivariatePolynomial::new(domain);
        let prod = a.mul(&zero);
        assert!(prod.is_zero());
    }

    #[test]
    fn rational_polynomial_multiplication() {
        let domain = RationalDomain;
        let a = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(1, 2), Rational::new(1, 1)],
        );
        let b = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(2, 1), Rational::new(1, 1)],
        );
        let prod = a.mul(&b);
        // (1/2 + x) * (2 + x) = 1 + (5/2)x + x^2
        assert_eq!(prod.coeff(0).cloned(), Some(Rational::new(1, 1)));
        assert_eq!(prod.coeff(1).cloned(), Some(Rational::new(5, 2)));
        assert_eq!(prod.coeff(2).cloned(), Some(Rational::new(1, 1)));
    }

    #[test]
    fn finite_field_polynomial_arithmetic() {
        let domain = FiniteField::new(num_bigint::BigInt::from(7));
        let a = DenseUnivariatePolynomial::from_coeffs(
            domain.clone(),
            vec![domain.element(3), domain.element(1)],
        );
        let b = DenseUnivariatePolynomial::from_coeffs(
            domain.clone(),
            vec![domain.element(2), domain.element(0), domain.element(1)],
        );
        let prod = a.mul(&b);
        // (3 + x) * (2 + x^2) = 6 + 2x + 3x^2 + x^3  (mod 7)
        assert_eq!(prod.coeff(0).cloned(), Some(domain.element(6)));
        assert_eq!(prod.coeff(1).cloned(), Some(domain.element(2)));
        assert_eq!(prod.coeff(2).cloned(), Some(domain.element(3)));
        assert_eq!(prod.coeff(3).cloned(), Some(domain.element(1)));
    }

    #[test]
    fn evaluate_polynomial() {
        let domain = IntegerDomain;
        // p(x) = 1 + 2x + 3x^2
        let p = DenseUnivariatePolynomial::from_coeffs(domain, vec![int(1), int(2), int(3)]);
        assert_eq!(p.eval(&int(0)), int(1));
        assert_eq!(p.eval(&int(1)), int(6));
        assert_eq!(p.eval(&int(2)), int(17));
    }

    #[test]
    fn polynomial_derivative() {
        let domain = IntegerDomain;
        // p(x) = 1 + 2x + 3x^2 + 4x^3 -> p'(x) = 2 + 6x + 12x^2
        let p =
            DenseUnivariatePolynomial::from_coeffs(domain, vec![int(1), int(2), int(3), int(4)]);
        let dp = p.derivative();
        assert_eq!(dp.coeffs().to_vec(), vec![int(2), int(6), int(12)]);
    }

    #[test]
    fn polynomial_integral() {
        let domain = RationalDomain;
        // p(x) = 1 + 2x -> int p = 0 + x + x^2
        let p = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(1, 1), Rational::new(2, 1)],
        );
        let ip = p.integral();
        assert_eq!(ip.coeff(0).cloned(), Some(Rational::new(0, 1)));
        assert_eq!(ip.coeff(1).cloned(), Some(Rational::new(1, 1)));
        assert_eq!(ip.coeff(2).cloned(), Some(Rational::new(1, 1)));
    }

    #[test]
    fn polynomial_division_with_remainder_over_integers() {
        let domain = IntegerDomain;
        // (x^2 + 1) / (x - 1) = x + 1 remainder 2
        let dividend = DenseUnivariatePolynomial::from_coeffs(domain, vec![int(1), int(0), int(1)]);
        let divisor = DenseUnivariatePolynomial::from_coeffs(domain, vec![int(-1), int(1)]);
        let (q, r) = dividend.div_rem(&divisor).unwrap();
        assert_eq!(q.coeffs().to_vec(), vec![int(1), int(1)]);
        assert_eq!(r.coeffs().to_vec(), vec![int(2)]);
    }

    #[test]
    fn polynomial_division_exact_over_rationals() {
        let domain = RationalDomain;
        // (x^2 - 1) / (x - 1) = x + 1, remainder 0
        let dividend = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![
                Rational::new(-1, 1),
                Rational::new(0, 1),
                Rational::new(1, 1),
            ],
        );
        let divisor = DenseUnivariatePolynomial::from_coeffs(
            domain,
            vec![Rational::new(-1, 1), Rational::new(1, 1)],
        );
        let (q, r) = dividend.div_rem(&divisor).unwrap();
        assert_eq!(q.degree(), Some(1));
        assert_eq!(q.coeff(0).cloned(), Some(Rational::new(1, 1)));
        assert_eq!(q.coeff(1).cloned(), Some(Rational::new(1, 1)));
        assert!(r.is_zero());
    }

    #[test]
    fn karatsuba_large_multiplication() {
        let d = IntegerDomain;
        // Create two degree-100 polynomials (exceeds KARATSUBA_THRESHOLD=32)
        let coeffs_a: Vec<Integer> = (0..=100).map(|i| Integer::from(i as i64)).collect();
        let coeffs_b: Vec<Integer> = (0..=100).map(|i| Integer::from((i + 1) as i64)).collect();
        let a = DenseUnivariatePolynomial::from_coeffs(d, coeffs_a);
        let b = DenseUnivariatePolynomial::from_coeffs(d, coeffs_b);

        let c = a.mul(&b);

        // Degree should be 200
        assert_eq!(c.degree(), Some(200));
        // Leading coefficient = 100 * 101 = 10100
        assert_eq!(c.leading_coeff(), Some(&Integer::from(10100)));
        // Constant term = 0 * 1 = 0
        assert_eq!(c.constant(), Integer::from(0));
        // Coefficient of x^1 = a0*b1 + a1*b0 = 0*2 + 1*1 = 1
        assert_eq!(c.coeff(1).cloned(), Some(Integer::from(1)));
    }

    #[test]
    fn karatsuba_cross_check_with_schoolbook() {
        use ocas_domain::RationalDomain;
        let d = RationalDomain;
        // Two degree-50 rational polynomials to force Karatsuba path
        let coeffs_a: Vec<Rational> = (0..=50)
            .map(|i| Rational::new(i as i64, (i + 1) as i64))
            .collect();
        let coeffs_b: Vec<Rational> = (0..=50)
            .map(|i| Rational::new((i + 1) as i64, (i + 2) as i64))
            .collect();
        let a = DenseUnivariatePolynomial::from_coeffs(d, coeffs_a.clone());
        let b = DenseUnivariatePolynomial::from_coeffs(d, coeffs_b.clone());

        let c_karat = a.mul(&b);

        // Compute with schoolbook manually for cross-check
        let result_len = coeffs_a.len() + coeffs_b.len() - 1;
        let mut expected = vec![Rational::new(0, 1); result_len];
        for (i, ai) in coeffs_a.iter().enumerate() {
            for (j, bj) in coeffs_b.iter().enumerate() {
                let prod = d.mul(ai, bj);
                expected[i + j] = d.add(&expected[i + j], &prod);
            }
        }

        assert_eq!(c_karat.coeffs(), &expected[..]);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use ocas_domain::{Integer, IntegerDomain};
    use proptest::prelude::*;

    fn any_int_poly(
        max_degree: usize,
    ) -> impl Strategy<Value = DenseUnivariatePolynomial<IntegerDomain>> {
        prop::collection::vec(any::<i64>(), 0..=max_degree).prop_map(|v| {
            let domain = IntegerDomain;
            DenseUnivariatePolynomial::from_coeffs(
                domain,
                v.into_iter().map(Integer::from).collect(),
            )
        })
    }

    proptest! {
        #[test]
        fn addition_is_commutative(a in any_int_poly(8), b in any_int_poly(8)) {
            assert_eq!(a.add(&b), b.add(&a));
        }

        #[test]
        fn multiplication_is_commutative(a in any_int_poly(5), b in any_int_poly(5)) {
            assert_eq!(a.mul(&b), b.mul(&a));
        }

        #[test]
        fn derivative_reduces_degree_or_zero(a in any_int_poly(6)) {
            let da = a.derivative();
            match a.degree() {
                None => assert!(da.is_zero()),
                Some(0) => assert!(da.is_zero()),
                Some(d) => assert!(da.degree().unwrap_or(0) < d),
            }
        }

        #[test]
        fn mul_then_div_exact_when_divisor_is_factor(
            a in any_int_poly(4),
            b in any_int_poly(4),
        ) {
            // Ensure b is not zero.
            prop_assume!(!b.is_zero());
            let prod = a.mul(&b);
            let (q, r) = prod.div_rem(&b).unwrap();
            assert!(r.is_zero());
            assert_eq!(q, a);
        }
    }
}
