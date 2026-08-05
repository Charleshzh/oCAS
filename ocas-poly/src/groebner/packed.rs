//! SWAR-packed monomial representation for the F5 ℤ_p fast path.
//!
//! A monomial in ≤ 8 variables is packed into a single `u128` as eight
//! 16-bit exponent fields, variable `i` occupying bits
//! `[112 - 16i, 127 - 16i]` (variable 0 in the highest field, so plain
//! integer comparison of packed values is exactly Lex order). All
//! arithmetic (exponent addition, divisibility, lcm, support masks) runs
//! as a handful of SIMD-within-a-register `u128` operations, eliminating
//! the per-monomial heap allocations and per-term `usize`-slice churn of
//! the generic path.
//!
//! Eligibility is enforced at the pipeline boundary ([`super::f5::packed_eligible`]):
//! `n_vars ≤ 8` and every input exponent `< 32768`. The SWAR formulas below
//! (divides/lcm/support_mask) are the standard 15-bit-field tricks: they
//! are exact only when all field values are `< 2^15`, so the eligibility
//! bound is 32768, not 65536. All matrix monomials in F5 are bounded
//! fieldwise by the input exponents (lcm closure), so the bound on inputs
//! covers every intermediate monomial; `debug_assert!`s additionally guard
//! the `add` invariant.

use std::cmp::Ordering;

use ocas_core::FastHashMap as HashMap;
use smallvec::SmallVec;

use crate::sparse::MonomialOrder;

/// High bit of every 16-bit field (bit 127, 111, …, 15).
const H: u128 = 0x8000_8000_8000_8000_8000_8000_8000_8000;
/// `1` in the lowest bit of every 16-bit field.
const ONE: u128 = 0x0001_0001_0001_0001_0001_0001_0001_0001;

/// SWAR-packed monomial: `u128` holding up to 8 exponent fields of 16 bits.
///
/// Variable `i` lives in bits `[112 - 16i, 127 - 16i]`; fields above
/// `n_vars` are zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct PackedMono(pub u128);

impl PackedMono {
    /// Pack an exponent vector. `exp.len()` must be ≤ 8 and every
    /// exponent < 2^16.
    #[inline]
    pub(crate) fn pack(exp: &[usize]) -> Self {
        debug_assert!(exp.len() <= 8, "packed monomials support ≤ 8 variables");
        let mut v = 0u128;
        for (i, &e) in exp.iter().enumerate() {
            debug_assert!(e < 1 << 16, "packed exponent field overflow: {e}");
            v |= (e as u128) << (112 - 16 * i);
        }
        Self(v)
    }

    /// Unpack the first `n_vars` fields into an inline SmallVec.
    #[inline]
    pub(crate) fn unpack_sv(&self, n_vars: usize) -> SmallVec<[usize; 8]> {
        debug_assert!(n_vars <= 8);
        let mut out: SmallVec<[usize; 8]> = SmallVec::with_capacity(n_vars);
        for i in 0..n_vars {
            out.push(((self.0 >> (112 - 16 * i)) & 0xFFFF) as usize);
        }
        out
    }

    /// Componentwise exponent addition (one `u128` add in release).
    #[inline]
    pub(crate) fn add(self, other: Self) -> Self {
        debug_assert!(
            {
                // Per-field no-overflow check (debug builds only).
                let mut ok = true;
                for i in 0..8 {
                    let a = (self.0 >> (112 - 16 * i)) & 0xFFFF;
                    let b = (other.0 >> (112 - 16 * i)) & 0xFFFF;
                    ok &= a + b < 0x10000;
                }
                ok
            },
            "packed exponent addition overflow"
        );
        Self(self.0.wrapping_add(other.0))
    }

    /// Componentwise exponent subtraction. The caller guarantees
    /// `other.divides(self)` (i.e. `self` is a multiple of `other`).
    #[inline]
    pub(crate) fn sub(self, other: Self) -> Self {
        Self(self.0.wrapping_sub(other.0))
    }

    /// True iff `other` divides `self` (componentwise `other ≤ self`);
    /// same semantics as `monomial_divides(self, other)`.
    #[inline]
    pub(crate) fn divides(self, other: Self) -> bool {
        (self.0 | H).wrapping_sub(other.0) & H == H
    }

    /// Componentwise maximum (the lcm of the two monomials).
    ///
    /// Currently unused by the pipeline (pair lcms are packed directly);
    /// kept as part of the packed primitive API and covered by unit tests.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn lcm(self, other: Self) -> Self {
        let m = (self.0 | H).wrapping_sub(other.0);
        let keep = m & H;
        let full = (keep >> 15).wrapping_mul(0xFFFF);
        Self((self.0 & full) | (other.0 & !full))
    }

    /// Bit `i` of the result is set iff field `i` is nonzero.
    ///
    /// Field `i` (variable `i`) has its high bit at `127 - 16i`; the fold
    /// shifts each high bit down to bit `i`, matching the SmallVec
    /// `support_mask` convention (bit `v` ⟺ variable `v`) so packed masks
    /// can be mixed with [`super::f4::DivisorIndex`] buckets.
    ///
    /// Exact for fields ≤ 0x8000 (the high-bit trick `(x | H) - 1` is the
    /// standard 15-bit-field SWAR form); the pipeline eligibility bound
    /// keeps every field < 2^15.
    #[inline]
    pub(crate) fn support_mask(self) -> u64 {
        let nz = (self.0 | H).wrapping_sub(ONE) & H;
        ((nz >> 127)
            | (nz >> 110)
            | (nz >> 93)
            | (nz >> 76)
            | (nz >> 59)
            | (nz >> 42)
            | (nz >> 25)
            | (nz >> 8)) as u64
            & 0xFF
    }
}

/// A signature in the F5 algorithm with a packed monomial: same semantics
/// as [`super::f5::Signature`], but multiplications by monomials are
/// single `u128` adds with no allocation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct PackedSig {
    /// Index of the input generator (module position, 0-based).
    pub pos: usize,
    /// Multiplier monomial on that generator's module basis vector.
    pub mono: PackedMono,
}

impl PackedSig {
    /// The canonical signature of generator `pos`: the unit monomial.
    #[inline]
    pub(crate) fn unit(pos: usize) -> Self {
        Self {
            pos,
            mono: PackedMono(0),
        }
    }

    /// Multiply this signature by a monomial (module position unchanged).
    #[inline]
    pub(crate) fn mul_monomial(&self, m: PackedMono) -> Self {
        Self {
            pos: self.pos,
            mono: self.mono.add(m),
        }
    }

    /// Compare under the pot (position-over-term) order.
    #[inline]
    pub(crate) fn cmp_pot<O: MonomialOrder>(
        &self,
        other: &Self,
        order: &O,
        n_vars: usize,
    ) -> Ordering {
        self.pos
            .cmp(&other.pos)
            .then_with(|| order.cmp(&self.mono.unpack_sv(n_vars), &other.mono.unpack_sv(n_vars)))
    }
}

/// A set of packed monomials supporting fast "some stored monomial divides
/// `exp`" queries, bucketed by exact support mask (the `PackedMono`
/// analogue of [`super::f5::MonomialBucketSet`]).
pub(crate) struct PackedMonomialBucketSet {
    /// Exact support mask of a monomial → monomials with that mask.
    buckets: HashMap<u64, Vec<PackedMono>>,
}

impl Default for PackedMonomialBucketSet {
    fn default() -> Self {
        Self::new()
    }
}

impl PackedMonomialBucketSet {
    pub(crate) fn new() -> Self {
        Self {
            buckets: HashMap::default(),
        }
    }

    pub(crate) fn insert(&mut self, m: PackedMono) {
        self.buckets.entry(m.support_mask()).or_default().push(m);
    }

    /// True iff some stored monomial divides `exp` (including `exp`
    /// itself and the unit monomial, matching the SmallVec-based set).
    pub(crate) fn any_divisor_of(&self, exp: PackedMono) -> bool {
        let mask = exp.support_mask();
        let mut sub = mask;
        loop {
            if let Some(ms) = self.buckets.get(&sub)
                && ms.iter().any(|m| exp.divides(*m))
            {
                return true;
            }
            if sub == 0 {
                break;
            }
            sub = (sub - 1) & mask;
        }
        false
    }
}

// =========================================================================
//  Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference support mask computed field-by-field on the unpacked slice.
    fn ref_support_mask(exp: &[usize]) -> u64 {
        let mut mask = 0u64;
        for (i, &e) in exp.iter().enumerate() {
            if e > 0 {
                mask |= 1 << i;
            }
        }
        mask
    }

    /// Reference lcm computed field-by-field.
    fn ref_lcm(a: &[usize], b: &[usize]) -> Vec<usize> {
        a.iter().zip(b.iter()).map(|(x, y)| (*x).max(*y)).collect()
    }

    /// Reference divisibility: `b` divides `a`.
    fn ref_divides(a: &[usize], b: &[usize]) -> bool {
        a.iter().zip(b.iter()).all(|(x, y)| x >= y)
    }

    /// All exponent vectors of length `n_vars` with each field in
    /// `[0, bound]`, for exhaustive small-space checks. Only usable with
    /// small `bound` (exhaustive 8-variable spaces grow as `(bound+1)^8`).
    fn exps(n_vars: usize, bound: usize) -> Vec<Vec<usize>> {
        let mut out = vec![vec![]];
        for _ in 0..n_vars {
            let mut next = Vec::new();
            for e in &out {
                for v in 0..=bound {
                    let mut c = e.clone();
                    c.push(v);
                    next.push(c);
                }
            }
            out = next;
        }
        out
    }

    /// Deterministic pseudo-random exponent vectors (LCG, no external
    /// deps) within the formulas' exactness contract (fields < 2^15):
    /// small values, the 32766/32767 top of the contract, and mid-range.
    fn random_exps(n_vars: usize, count: usize, seed: u64) -> Vec<Vec<usize>> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        (0..count)
            .map(|_| {
                (0..n_vars)
                    .map(|_| match next() % 10 {
                        0..=3 => (next() % 3) as usize,      // small values
                        4..=5 => 32766usize,                 // near top of contract
                        6..=7 => 32767usize,                 // top of contract
                        _ => 3 + ((next() as usize) % 1024), // mid-range
                    })
                    .collect()
            })
            .collect()
    }

    /// Random vectors with full 16-bit fields (65534/65535) — only valid
    /// where the operation is exact for all 16-bit values (pack/unpack).
    fn random_exps_full16(n_vars: usize, count: usize, seed: u64) -> Vec<Vec<usize>> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        (0..count)
            .map(|_| {
                (0..n_vars)
                    .map(|_| match next() % 10 {
                        0..=3 => (next() % 3) as usize,
                        4..=5 => 65534usize,
                        6..=7 => 65535usize,
                        _ => 3 + ((next() as usize) % 1024),
                    })
                    .collect()
            })
            .collect()
    }

    /// Test vectors for `add`: exhaustive small space plus random
    /// mid-range/edge values (fields ≤ 32767 keep sums < 2^16).
    fn add_vecs(n_vars: usize) -> Vec<Vec<usize>> {
        let mut out = exps(n_vars, 2);
        if n_vars <= 4 {
            return out;
        }
        out.extend(random_exps(n_vars, 120, 0xADD + n_vars as u64));
        out
    }

    #[test]
    fn pack_unpack_roundtrip() {
        for n_vars in 1..=8 {
            for exp in exps(n_vars, 2) {
                let m = PackedMono::pack(&exp);
                assert_eq!(m.unpack_sv(n_vars).as_slice(), exp.as_slice());
            }
            // Edge fields: 65534/65535 pack and roundtrip exactly
            // (pack/unpack are exact for all 16-bit values).
            for exp in random_exps_full16(n_vars, 60, 0x5EED + n_vars as u64) {
                let m = PackedMono::pack(&exp);
                assert_eq!(m.unpack_sv(n_vars).as_slice(), exp.as_slice());
            }
        }
    }

    #[test]
    fn pack_fields_in_expected_slots() {
        // Variable 0 lives in the highest field: x0^1 packs to bit 112.
        assert_eq!(PackedMono::pack(&[1]).0, 1u128 << 112);
        assert_eq!(PackedMono::pack(&[0, 1]).0, 1u128 << 96);
        assert_eq!(PackedMono::pack(&[0, 0, 0, 0, 0, 0, 0, 1]).0, 1u128);
        // Integer comparison is Lex order: [0,5] < [1,0] because the
        // first component (variable 0, highest field) decides.
        assert!(PackedMono::pack(&[0, 5]).0 < PackedMono::pack(&[1, 0]).0);
        assert!(PackedMono::pack(&[1, 0]).0 < PackedMono::pack(&[2, 0]).0);
    }

    #[test]
    fn add_sub_roundtrip() {
        for n_vars in 1..=8 {
            let vecs = add_vecs(n_vars);
            for a in &vecs {
                for b in &vecs {
                    let pa = PackedMono::pack(a);
                    let pb = PackedMono::pack(b);
                    let sum: Vec<usize> = a.iter().zip(b.iter()).map(|(x, y)| x + y).collect();
                    let ps = pa.add(pb);
                    assert_eq!(ps.unpack_sv(n_vars).as_slice(), sum.as_slice());
                    // a + b - b == a (exact subtraction).
                    assert_eq!(ps.sub(pb), pa);
                }
            }
        }
    }

    #[test]
    fn sub_borrow_chains() {
        // Exact subtraction across 0xFFFF fields exercises the borrow
        // propagation through lower fields.
        for n_vars in 1..=8 {
            for _ in 0..400 {
                let mut a = vec![0usize; n_vars];
                let mut b = vec![0usize; n_vars];
                for i in 0..n_vars {
                    a[i] = [0, 1, 65534, 65535][(i * 7 + n_vars) % 4];
                    b[i] = if a[i] > 0 { a[i] - 1 } else { 0 };
                }
                let pa = PackedMono::pack(&a);
                let pb = PackedMono::pack(&b);
                let diff: Vec<usize> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
                assert_eq!(pa.sub(pb).unpack_sv(n_vars).as_slice(), diff.as_slice());
            }
        }
    }

    #[test]
    fn divides_matches_reference() {
        for n_vars in 1..=8 {
            let small = exps(n_vars, 2);
            for a in &small {
                for b in &small {
                    let pa = PackedMono::pack(a);
                    let pb = PackedMono::pack(b);
                    assert_eq!(
                        pa.divides(pb),
                        ref_divides(a, b),
                        "divides mismatch: {a:?} vs {b:?}"
                    );
                }
            }
            let randoms = random_exps(n_vars, 150, 0xD17 + n_vars as u64);
            for a in &randoms {
                for b in &randoms {
                    let pa = PackedMono::pack(a);
                    let pb = PackedMono::pack(b);
                    assert_eq!(
                        pa.divides(pb),
                        ref_divides(a, b),
                        "divides mismatch (random): {a:?} vs {b:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn lcm_matches_reference() {
        for n_vars in 1..=8 {
            let small = exps(n_vars, 2);
            for a in &small {
                for b in &small {
                    let pa = PackedMono::pack(a);
                    let pb = PackedMono::pack(b);
                    let l = pa.lcm(pb);
                    assert_eq!(
                        l.unpack_sv(n_vars).as_slice(),
                        ref_lcm(a, b).as_slice(),
                        "lcm mismatch: {a:?} lcm {b:?}"
                    );
                }
            }
            let randoms = random_exps(n_vars, 150, 0x1CD + n_vars as u64);
            for a in &randoms {
                for b in &randoms {
                    let pa = PackedMono::pack(a);
                    let pb = PackedMono::pack(b);
                    let l = pa.lcm(pb);
                    assert_eq!(
                        l.unpack_sv(n_vars).as_slice(),
                        ref_lcm(a, b).as_slice(),
                        "lcm mismatch (random): {a:?} lcm {b:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn support_mask_matches_reference() {
        for n_vars in 1..=8 {
            for exp in exps(n_vars, 2) {
                let m = PackedMono::pack(&exp);
                assert_eq!(
                    m.support_mask(),
                    ref_support_mask(&exp),
                    "support_mask mismatch: {exp:?}"
                );
            }
            for exp in random_exps(n_vars, 80, 0x5EED + n_vars as u64) {
                let m = PackedMono::pack(&exp);
                assert_eq!(
                    m.support_mask(),
                    ref_support_mask(&exp),
                    "support_mask mismatch (random): {exp:?}"
                );
            }
        }
    }

    #[test]
    fn packed_sig_basic() {
        use crate::sparse::Lex;
        let s0 = PackedSig::unit(0);
        let s1 = PackedSig::unit(1);
        let order = Lex;
        // pot: module position dominates.
        assert_eq!(s0.cmp_pot::<Lex>(&s1, &order, 2), Ordering::Less);
        // Same position: compare monomials under O.
        let a = PackedSig {
            pos: 0,
            mono: PackedMono::pack(&[0, 1]),
        };
        let b = PackedSig {
            pos: 0,
            mono: PackedMono::pack(&[1, 0]),
        };
        assert_eq!(a.cmp_pot::<Lex>(&b, &order, 2), Ordering::Less);
        // mul_monomial keeps position, adds exponents.
        let m = PackedMono::pack(&[1, 1]);
        let prod = s1.mul_monomial(m);
        assert_eq!(prod.pos, 1);
        assert_eq!(prod.mono.unpack_sv(2).as_slice(), &[1, 1]);
    }

    #[test]
    fn bucket_set_divisor_queries() {
        let mut set = PackedMonomialBucketSet::new();
        set.insert(PackedMono::pack(&[2, 0]));
        assert!(set.any_divisor_of(PackedMono::pack(&[2, 0]))); // itself
        assert!(set.any_divisor_of(PackedMono::pack(&[3, 1]))); // multiple
        assert!(set.any_divisor_of(PackedMono::pack(&[2, 0, 5, 0, 0, 0, 0, 1])));
        assert!(!set.any_divisor_of(PackedMono::pack(&[1, 0]))); // smaller
        assert!(!set.any_divisor_of(PackedMono::pack(&[0, 2]))); // other var
        // Unit monomial divides everything.
        set.insert(PackedMono(0));
        assert!(set.any_divisor_of(PackedMono::pack(&[0, 0, 0, 0, 0, 0, 0, 0])));
    }
}
