//! Algebraic extension domain $K = D[\alpha]/(m(\alpha))$.
//!
//! An [`AlgebraicExtension`] is the quotient of a univariate polynomial ring
//! by a monic polynomial $m$. When the base domain $D$ is a field and $m$ is
//! irreducible, the quotient is a field:
//!
//! - `AlgebraicExtension<RationalDomain>` is an algebraic number field
//!   $\mathbb{Q}(\alpha)$;
//! - `AlgebraicExtension<FiniteField>` is a Galois field $\mathrm{GF}(p^d)$.
//!
//! Elements are residue classes represented by their unique polynomial
//! representative of degree less than $\deg(m)$. Inversion uses the extended
//! Euclidean algorithm over the base field (self-contained dense polynomial
//! arithmetic, so that `ocas-domain` does not depend on `ocas-poly`).
//!
//! Irreducibility of the minimal polynomial is **not** checked (mirroring
//! Symbolica); over a reducible modulus the ring has zero divisors and
//! [`Domain::inv`] returns `None` for non-units.
//!
//! # Example
//!
//! ```
//! use ocas_domain::{AlgebraicExtension, Domain, Rational, RationalDomain};
//!
//! // ℚ(√2): minimal polynomial α² − 2.
//! let two = Rational::new(2, 1);
//! let neg_two = RationalDomain.neg(&two);
//! let field = AlgebraicExtension::new(
//!     RationalDomain,
//!     vec![neg_two, Rational::new(0, 1), Rational::new(1, 1)],
//! );
//! let sqrt2 = field.alpha();
//! // √2·√2 = 2.
//! assert_eq!(field.mul(&sqrt2, &sqrt2), field.from_base(two));
//! ```

use crate::domain::{Domain, EuclideanDomain};
use crate::rational::RationalDomain;

/// Dense univariate polynomial coefficients over `D`, ascending degree
/// order with trailing zeros trimmed.
type PolyCoeffs<D> = Vec<<D as Domain>::Element>;

/// An element of an algebraic extension: the residue class of a polynomial
/// in $\alpha$ of degree less than the extension degree.
///
/// Coefficients are stored in ascending degree order with trailing zeros
/// trimmed, so the zero element has an empty coefficient vector and
/// equality of representatives is semantic equality.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AlgebraicElement<E> {
    coeffs: Vec<E>,
}

impl<E> AlgebraicElement<E> {
    /// Coefficients in ascending degree order (trailing zeros trimmed).
    pub fn coeffs(&self) -> &[E] {
        &self.coeffs
    }
}

impl<E: std::fmt::Display> std::fmt::Display for AlgebraicElement<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }
        for (i, c) in self.coeffs.iter().enumerate().rev() {
            if i < self.coeffs.len() - 1 {
                write!(f, " + ")?;
            }
            match i {
                0 => write!(f, "{c}")?,
                1 => write!(f, "({c})·α")?,
                _ => write!(f, "({c})·α^{i}")?,
            }
        }
        Ok(())
    }
}

/// An algebraic extension $D[\alpha]/(m(\alpha))$ over a base domain $D$.
///
/// The minimal polynomial is stored in ascending degree order and must be
/// monic of degree ≥ 1; both are checked with `debug_assert`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlgebraicExtension<D: Domain> {
    base: D,
    min_poly: Vec<D::Element>,
}

impl<D: Domain> AlgebraicExtension<D> {
    /// Create the extension $D[\alpha]/(m)$ from a monic minimal polynomial
    /// `min_poly` in ascending degree order.
    ///
    /// Irreducibility is not verified; over a reducible modulus some
    /// nonzero elements are zero divisors and [`Domain::inv`] returns
    /// `None` for them.
    pub fn new(base: D, min_poly: Vec<D::Element>) -> Self {
        debug_assert!(
            min_poly.len() >= 2,
            "minimal polynomial must have degree at least 1"
        );
        debug_assert!(
            min_poly.last() == Some(&base.one()),
            "minimal polynomial must be monic"
        );
        Self { base, min_poly }
    }

    /// The base domain $D$.
    pub fn base_domain(&self) -> &D {
        &self.base
    }

    /// The minimal polynomial in ascending degree order (monic).
    pub fn min_poly(&self) -> &[D::Element] {
        &self.min_poly
    }

    /// The extension degree $\deg(m)$.
    pub fn extension_degree(&self) -> usize {
        self.min_poly.len() - 1
    }

    /// Embed a base-domain constant into the extension.
    pub fn from_base(&self, c: D::Element) -> AlgebraicElement<D::Element> {
        let mut coeffs = vec![c];
        self.trim(&mut coeffs);
        AlgebraicElement { coeffs }
    }

    /// The generator $\alpha$ of the extension.
    pub fn alpha(&self) -> AlgebraicElement<D::Element> {
        AlgebraicElement {
            coeffs: vec![self.base.zero(), self.base.one()],
        }
    }

    /// Create an element from arbitrary coefficients (ascending order),
    /// reduced modulo the minimal polynomial.
    pub fn element(&self, mut coeffs: Vec<D::Element>) -> AlgebraicElement<D::Element> {
        self.reduce(&mut coeffs);
        AlgebraicElement { coeffs }
    }

    /// Trim trailing zero coefficients.
    fn trim(&self, v: &mut Vec<D::Element>) {
        while let Some(last) = v.last() {
            if self.base.is_zero(last) {
                v.pop();
            } else {
                break;
            }
        }
    }

    /// Reduce coefficients modulo the (monic) minimal polynomial in place.
    fn reduce(&self, v: &mut Vec<D::Element>) {
        let d = self.extension_degree();
        while v.len() > d {
            let c = v.pop().expect("nonempty while reducing");
            if self.base.is_zero(&c) {
                continue;
            }
            let offset = v.len() - d;
            for (i, mc) in self.min_poly.iter().take(d).enumerate() {
                let t = self.base.mul(&c, mc);
                let slot = &mut v[offset + i];
                *slot = self.base.sub(slot, &t);
            }
        }
        self.trim(v);
    }

    /// Polynomial addition (truncating) over the base domain.
    fn poly_add(&self, a: &[D::Element], b: &[D::Element]) -> Vec<D::Element> {
        let mut out = Vec::with_capacity(a.len().max(b.len()));
        for i in 0..a.len().max(b.len()) {
            let x = a.get(i);
            let y = b.get(i);
            let c = match (x, y) {
                (Some(x), Some(y)) => self.base.add(x, y),
                (Some(x), None) => x.clone(),
                (None, Some(y)) => y.clone(),
                (None, None) => unreachable!(),
            };
            out.push(c);
        }
        self.trim(&mut out);
        out
    }

    /// Polynomial subtraction (truncating) over the base domain.
    fn poly_sub(&self, a: &[D::Element], b: &[D::Element]) -> Vec<D::Element> {
        let mut out = Vec::with_capacity(a.len().max(b.len()));
        for i in 0..a.len().max(b.len()) {
            let x = a.get(i);
            let y = b.get(i);
            let c = match (x, y) {
                (Some(x), Some(y)) => self.base.sub(x, y),
                (Some(x), None) => x.clone(),
                (None, Some(y)) => self.base.neg(y),
                (None, None) => unreachable!(),
            };
            out.push(c);
        }
        self.trim(&mut out);
        out
    }

    /// Schoolbook polynomial multiplication over the base domain.
    fn poly_mul(&self, a: &[D::Element], b: &[D::Element]) -> Vec<D::Element> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut out = vec![self.base.zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            if self.base.is_zero(x) {
                continue;
            }
            for (j, y) in b.iter().enumerate() {
                let t = self.base.mul(x, y);
                let slot = &mut out[i + j];
                *slot = self.base.add(slot, &t);
            }
        }
        self.trim(&mut out);
        out
    }

    /// Polynomial division with remainder over the base domain: returns
    /// `(quotient, remainder)` with `a = quotient·b + remainder` and
    /// `deg(remainder) < deg(b)`, or `None` when `b` is zero or a leading
    /// coefficient division fails (i.e. the base is not a field).
    fn poly_quot_rem(
        &self,
        a: &[D::Element],
        b: &[D::Element],
    ) -> Option<(PolyCoeffs<D>, PolyCoeffs<D>)> {
        if b.is_empty() {
            return None;
        }
        let mut rem = a.to_vec();
        self.trim(&mut rem);
        let mut quot = vec![self.base.zero(); a.len().max(b.len()) - b.len() + 1];
        let deg_b = b.len() - 1;
        let lc_b = b.last().expect("nonempty divisor");
        while rem.len() > deg_b && !rem.is_empty() {
            let k = rem.len() - 1 - deg_b;
            let c = self
                .base
                .div(rem.last().expect("nonempty remainder"), lc_b)?;
            if !self.base.is_zero(&c) {
                quot[k] = self.base.add(&quot[k], &c);
                for (i, bc) in b.iter().enumerate() {
                    let t = self.base.mul(&c, bc);
                    let slot = &mut rem[k + i];
                    *slot = self.base.sub(slot, &t);
                }
            }
            self.trim(&mut rem);
        }
        self.trim(&mut quot);
        Some((quot, rem))
    }

    /// Extended Euclidean algorithm over the base field: returns
    /// `(g, s, t)` with `s·a + t·b = g` and `g` monic, or `None` when a
    /// leading coefficient inversion fails (the base is not a field).
    fn poly_extended_gcd(
        &self,
        a: &[D::Element],
        b: &[D::Element],
    ) -> Option<(PolyCoeffs<D>, PolyCoeffs<D>, PolyCoeffs<D>)> {
        let mut old_r = a.to_vec();
        let mut r = b.to_vec();
        let mut old_s = vec![self.base.one()];
        let mut s: Vec<D::Element> = Vec::new();
        while !r.is_empty() {
            let (q, rem) = self.poly_quot_rem(&old_r, &r)?;
            old_r = r;
            r = rem;
            let qs = self.poly_mul(&q, &s);
            let new_s = self.poly_sub(&old_s, &qs);
            old_s = s;
            s = new_s;
        }
        // Normalize g (and the Bezout coefficient) to be monic.
        let lc = old_r.last()?.clone();
        let lc_inv = self.base.inv(&lc)?;
        if !self.base.is_one(&lc) {
            for c in old_r.iter_mut().chain(old_s.iter_mut()) {
                *c = self.base.mul(c, &lc_inv);
            }
        }
        Some((old_r, old_s, Vec::new()))
    }
}

impl<D: Domain> Domain for AlgebraicExtension<D> {
    type Element = AlgebraicElement<D::Element>;

    fn zero(&self) -> Self::Element {
        AlgebraicElement { coeffs: Vec::new() }
    }

    fn one(&self) -> Self::Element {
        self.from_base(self.base.one())
    }

    fn add(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        AlgebraicElement {
            coeffs: self.poly_add(&a.coeffs, &b.coeffs),
        }
    }

    fn sub(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        AlgebraicElement {
            coeffs: self.poly_sub(&a.coeffs, &b.coeffs),
        }
    }

    fn neg(&self, a: &Self::Element) -> Self::Element {
        AlgebraicElement {
            coeffs: a.coeffs.iter().map(|c| self.base.neg(c)).collect(),
        }
    }

    fn mul(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        let mut coeffs = self.poly_mul(&a.coeffs, &b.coeffs);
        self.reduce(&mut coeffs);
        AlgebraicElement { coeffs }
    }

    fn div(&self, a: &Self::Element, b: &Self::Element) -> Option<Self::Element> {
        self.inv(b).map(|inv| self.mul(a, &inv))
    }

    fn inv(&self, a: &Self::Element) -> Option<Self::Element> {
        if self.is_zero(a) {
            return None;
        }
        // s·a + t·m = g over the base field; a is a unit iff deg(g) = 0,
        // in which case g = 1 (monic) and a⁻¹ = s mod m.
        let (g, s, _) = self.poly_extended_gcd(&a.coeffs, &self.min_poly)?;
        if !g.is_empty() && g.len() == 1 {
            let mut coeffs = s;
            self.reduce(&mut coeffs);
            Some(AlgebraicElement { coeffs })
        } else {
            None
        }
    }

    fn is_zero(&self, a: &Self::Element) -> bool {
        a.coeffs.is_empty()
    }

    fn cast_u64(&self, n: u64) -> Self::Element {
        self.from_base(self.base.cast_u64(n))
    }
}

impl<D: Domain> EuclideanDomain for AlgebraicExtension<D> {
    fn div_rem(
        &self,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Option<(Self::Element, Self::Element)> {
        // Over a field every nonzero element is a unit, so division is
        // exact and the remainder is always zero.
        self.div(a, b).map(|q| (q, self.zero()))
    }

    fn gcd(&self, a: &Self::Element, b: &Self::Element) -> Self::Element {
        // In a field the GCD is degenerate: 0 if both are zero, else 1.
        // This keeps `content()`/`primitive_part()` well-behaved for
        // polynomials over the extension (same convention as `FiniteField`).
        if self.is_zero(a) && self.is_zero(b) {
            self.zero()
        } else {
            self.one()
        }
    }
}

/// An algebraic number field $\mathbb{Q}(\alpha)$.
pub type AlgebraicNumberField = AlgebraicExtension<RationalDomain>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_field::FiniteField;
    use crate::rational::Rational;
    use num_bigint::BigInt;

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    /// ℚ(α) with min_poly α² − c.
    fn q_ext_sqrt(c: i64) -> AlgebraicNumberField {
        AlgebraicNumberField::new(RationalDomain, vec![r(-c, 1), r(0, 1), r(1, 1)])
    }

    #[test]
    fn sqrt2_arithmetic() {
        let field = q_ext_sqrt(2);
        let alpha = field.alpha();
        // α² = 2
        assert_eq!(field.mul(&alpha, &alpha), field.from_base(r(2, 1)));
        // (1 + α)(1 − α) = 1 − α² = −1
        let one = field.one();
        let a = field.add(&one, &alpha);
        let b = field.sub(&one, &alpha);
        assert_eq!(field.mul(&a, &b), field.from_base(r(-1, 1)));
    }

    #[test]
    fn sqrt2_inverse() {
        let field = q_ext_sqrt(2);
        let alpha = field.alpha();
        // α⁻¹ = α/2
        let inv = field.inv(&alpha).expect("α is a unit");
        assert_eq!(inv.coeffs(), &[r(0, 1), r(1, 2)], "1/√2 = √2/2");
        // (1 + 2α)·inv(1 + 2α) = 1
        let a = field.element(vec![r(1, 1), r(2, 1)]);
        let a_inv = field.inv(&a).expect("unit");
        assert_eq!(field.mul(&a, &a_inv), field.one());
    }

    #[test]
    fn gaussian_rationals() {
        // ℚ(i): α² + 1.
        let field = AlgebraicNumberField::new(RationalDomain, vec![r(1, 1), r(0, 1), r(1, 1)]);
        let i = field.alpha();
        // i² = −1
        assert_eq!(field.mul(&i, &i), field.from_base(r(-1, 1)));
        // (1 + i)² = 2i
        let one_plus_i = field.add(&field.one(), &i);
        let sq = field.mul(&one_plus_i, &one_plus_i);
        assert_eq!(sq.coeffs(), &[r(0, 1), r(2, 1)]);
        // (1 + i)⁻¹ = (1 − i)/2
        let inv = field.inv(&one_plus_i).expect("unit");
        assert_eq!(inv.coeffs(), &[r(1, 2), r(-1, 2)]);
    }

    #[test]
    fn cbrt2_inverse() {
        // ℚ(∛2): α³ − 2.
        let field =
            AlgebraicNumberField::new(RationalDomain, vec![r(-2, 1), r(0, 1), r(0, 1), r(1, 1)]);
        let alpha = field.alpha();
        // α⁻¹ = α²/2
        let inv = field.inv(&alpha).expect("unit");
        assert_eq!(inv.coeffs(), &[r(0, 1), r(0, 1), r(1, 2)]);
        // A full-degree element: (1 + α + α²)⁻¹ exists and inverts.
        let a = field.element(vec![r(1, 1), r(1, 1), r(1, 1)]);
        let a_inv = field.inv(&a).expect("unit");
        assert_eq!(field.mul(&a, &a_inv), field.one());
    }

    #[test]
    fn galois_field_gf9() {
        // GF(3²): α² + 1 is irreducible over 𝔽_3.
        let base = FiniteField::new(BigInt::from(3));
        let field = AlgebraicExtension::new(
            base.clone(),
            vec![base.element(1), base.element(0), base.element(1)],
        );
        let alpha = field.alpha();
        // α² = −1 = 2 (mod 3)
        assert_eq!(field.mul(&alpha, &alpha), field.from_base(base.element(2)));
        // (1 + α)⁻¹ · (1 + α) = 1; every nonzero element of GF(9) is a unit.
        let a = field.add(&field.one(), &alpha);
        let a_inv = field.inv(&a).expect("unit in GF(9)");
        assert_eq!(field.mul(&a, &a_inv), field.one());
        // The multiplicative group has order 8: (1+α)⁸ = 1.
        assert_eq!(field.pow(&a, 8), field.one());
    }

    #[test]
    fn reducible_modulus_has_zero_divisors() {
        // α² − 1 = (α − 1)(α + 1) over ℚ: α − 1 is a zero divisor.
        let field = q_ext_sqrt(1);
        let alpha = field.alpha();
        let a = field.sub(&alpha, &field.one());
        assert!(field.inv(&a).is_none(), "α − 1 is not a unit");
    }

    #[test]
    fn degenerate_gcd_and_div_rem() {
        let field = q_ext_sqrt(2);
        let a = field.alpha();
        let z = field.zero();
        assert_eq!(field.gcd(&a, &z), field.one());
        assert_eq!(field.gcd(&z, &z), field.zero());
        let (q, rem) = field.div_rem(&a, &a).expect("division by a unit");
        assert_eq!(q, field.one());
        assert_eq!(rem, field.zero());
        assert!(field.div_rem(&a, &z).is_none());
    }

    #[test]
    fn cast_and_element_reduction() {
        let field = q_ext_sqrt(2);
        assert_eq!(field.cast_u64(3), field.from_base(r(3, 1)));
        // Coefficients of degree ≥ deg(m) are reduced by `element`.
        let e = field.element(vec![r(0, 1), r(0, 1), r(1, 1)]); // α²
        assert_eq!(e, field.from_base(r(2, 1)));
        // Zero is the empty coefficient vector.
        assert!(field.is_zero(&field.element(vec![r(0, 1), r(0, 1)])));
    }
}
