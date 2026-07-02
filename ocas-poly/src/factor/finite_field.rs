//! Factorization over finite fields $\mathbb{F}_p$.
//!
//! Implements the Cantor–Zassenhaus algorithm — distinct-degree factorization
//! (DDF) followed by equal-degree factorization (EDF) — together with the
//! Berlekamp algorithm for small primes. The top-level entry point
//! [`factor_over_finite_field`] first performs square-free factorization and
//! then dispatches to these algorithms.
//!
//! # Algorithm outline
//!
//! For a monic square-free $f \in \mathbb{F}_p[x]$:
//!
//! 1. **DDF** splits $f$ into products $g_d$ where every irreducible factor of
//!    $g_d$ has degree exactly $d$. It uses the Frobenius map $x \mapsto
//!    x^{p} \bmod f$: the irreducible factors of degree $d$ are precisely those
//!    dividing $\gcd(f, x^{p^d} - x)$ but not the earlier $g_{d'}$ ($d' < d$).
//!
//! 2. **EDF** splits a product of equal-degree irreducibles by picking random
//!    polynomials $a$ and computing $\gcd(f, a^{(p^d-1)/2} - 1)$ in odd
//!    characteristic (a trace map in characteristic 2).
//!
//! References: Cantor & Zassenhaus (1981); Berlekamp (1970); von zur Gathen &
//! Gerhard, *Modern Computer Algebra*, ch. 14.

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use ocas_domain::{Domain, FiniteField};

use crate::dense::DenseUnivariatePolynomial;
use crate::matrix::Matrix;

/// Convenience alias for univariate polynomials over a finite field.
pub type FpPoly = DenseUnivariatePolynomial<FiniteField>;

/// Normalize a polynomial to monic form by dividing by its leading
/// coefficient. The zero polynomial is returned unchanged.
fn monic(f: &FpPoly) -> FpPoly {
    if f.is_zero() {
        return f.zero();
    }
    let lc = f.leading_coeff().unwrap();
    if f.domain().is_one(lc) {
        return f.clone();
    }
    let inv = f.domain().inv(lc).expect("leading coefficient is nonzero");
    f.mul_scalar(&inv)
}

/// Return the polynomial $x$ (i.e. $0 + 1\cdot x$) over the given field.
fn x_var(field: &FiniteField) -> FpPoly {
    FpPoly::from_coeffs(field.clone(), vec![field.zero(), field.one()])
}

/// Compute $\text{base}^{\text{exp}} \bmod \text{modulus}$ via repeated
/// squaring, reducing modulo `modulus` after every multiplication.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use ocas_domain::{Domain, FiniteField};
/// use ocas_poly::factor::finite_field::poly_pow_mod;
///
/// let f = FiniteField::new(BigInt::from(7));
/// // modulus = x^2 + 1 over F_7
/// let m = {
///     use ocas_poly::DenseUnivariatePolynomial;
///     DenseUnivariatePolynomial::from_coeffs(f.clone(), vec![f.element(1), f.element(0), f.element(1)])
/// };
/// // base = x
/// let base = {
///     use ocas_poly::DenseUnivariatePolynomial;
///     DenseUnivariatePolynomial::from_coeffs(f.clone(), vec![f.element(0), f.element(1)])
/// };
/// // x^2 mod (x^2+1) = -1 = 6
/// let r = poly_pow_mod(&base, &BigInt::from(2), &m);
/// assert_eq!(r.coeff(0).cloned(), Some(f.element(6)));
/// ```
pub fn poly_pow_mod(base: &FpPoly, exp: &BigInt, modulus: &FpPoly) -> FpPoly {
    let field = base.domain().clone();
    if modulus.is_zero() {
        return base.clone();
    }
    let mut result = FpPoly::from_coeffs(field.clone(), vec![field.one()]);
    let mut b = match base.div_rem(modulus) {
        Some((_, r)) => r,
        None => base.clone(),
    };
    let mut e = exp.clone();
    while !e.is_zero() {
        if (&e & &BigInt::one()) == BigInt::one() {
            result = result.mul(&b);
            if let Some((_, r)) = result.div_rem(modulus) {
                result = r;
            }
        }
        e >>= 1;
        if !e.is_zero() {
            b = b.mul(&b);
            if let Some((_, r)) = b.div_rem(modulus) {
                b = r;
            }
        }
    }
    result
}

/// Distinct-degree factorization (DDF).
///
/// Given a monic square-free polynomial $f$, returns a list of `(g_d, d)` where
/// $g_d$ is the product of all irreducible factors of $f$ of degree exactly
/// $d$. The input is assumed monic and square-free; callers normalize first.
fn distinct_degree_factorization(f: &FpPoly) -> Vec<(FpPoly, usize)> {
    let field = f.domain().clone();
    let p = field.prime().clone();
    let mut result = Vec::new();
    let mut current = f.clone();
    let mut degree = 1usize;
    // h tracks x^(p^degree) mod current via Frobenius iteration.
    let mut h = x_var(&field); // x^(p^0) = x

    while let Some(deg_current) = current.degree() {
        if deg_current < 2 * degree {
            break;
        }
        // Advance h to x^(p^degree) mod current by raising to the p-th power.
        h = poly_pow_mod(&h, &p, &current);
        // g = gcd(current, h - x)
        let h_minus_x = h.sub(&x_var(&field));
        let g = current.gcd(&h_minus_x);
        let g = monic(&g);
        if g.degree().unwrap_or(0) > 0 {
            result.push((g.clone(), degree));
            // Remove the found factors from current and reduce h.
            if let Some((q, _)) = current.div_rem(&g) {
                current = monic(&q);
                if let Some((_, rh)) = h.div_rem(&current) {
                    h = rh;
                }
            }
        }
        degree += 1;
    }
    if current.degree().unwrap_or(0) > 0 && !current.is_one() {
        // The remaining factor is irreducible of its full degree.
        result.push((monic(&current), current.degree().unwrap()));
    }
    result
}

/// A tiny deterministic xorshift64 PRNG, used to generate pseudo-random
/// candidate polynomials for EDF without pulling in a `rand` dependency.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        // Avoid the degenerate all-zero state.
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
}

/// Build a pseudo-random polynomial of degree `< degree_bound` over the field,
/// seeded by `seed`. The polynomial is nonzero with overwhelming probability.
fn random_candidate(field: &FiniteField, degree_bound: usize, seed: u64) -> FpPoly {
    let p = field
        .prime()
        .to_u64()
        .expect("EDF candidate generation requires a prime fitting in u64");
    let mut rng = Rng::new(seed.wrapping_add(0x9E3779B97F4A7C15));
    let mut coeffs = Vec::with_capacity(degree_bound);
    for _ in 0..degree_bound {
        coeffs.push(field.element(BigInt::from(rng.next_u64() % p)));
    }
    FpPoly::from_coeffs(field.clone(), coeffs)
}

/// Equal-degree factorization (EDF).
///
/// Given a square-free polynomial `f` that is a product of irreducible factors
/// each of degree `d`, returns the list of those irreducible factors.
fn equal_degree_factorization(f: &FpPoly, d: usize) -> Vec<FpPoly> {
    let field = f.domain().clone();
    let p_big = field.prime().clone();
    let deg_f = f.degree().unwrap_or(0);
    if deg_f == 0 || d == 0 || deg_f == d {
        return vec![monic(f)];
    }

    let char_two = p_big == BigInt::from(2u32);

    // Exponent for odd characteristic: (p^d - 1) / 2.
    let half_exp = if !char_two {
        Some((p_big.pow(d as u32) - BigInt::one()) / 2)
    } else {
        None
    };

    let mut factors: Vec<FpPoly> = vec![monic(f)];
    let mut round: u64 = 0;

    loop {
        let mut made_progress = false;
        let mut new_factors = Vec::new();
        for (idx, factor) in factors.iter().enumerate() {
            let deg = factor.degree().unwrap_or(0);
            if deg == d || deg == 0 {
                new_factors.push(factor.clone());
                continue;
            }
            let r = deg / d; // number of degree-d factors
            if r <= 1 {
                new_factors.push(factor.clone());
                continue;
            }
            let split = if char_two {
                edf_split_char2(factor, d, round, idx as u64)
            } else {
                edf_split_odd(factor, half_exp.as_ref().unwrap(), round, idx as u64)
            };
            if let Some((g1, g2)) = split {
                new_factors.push(monic(&g1));
                new_factors.push(monic(&g2));
                made_progress = true;
            } else {
                new_factors.push(factor.clone());
            }
        }
        factors = new_factors;
        // Done when every factor has degree d (or 0).
        if factors
            .iter()
            .all(|f| f.degree().unwrap_or(0) == 0 || f.degree().unwrap_or(0) == d)
        {
            break;
        }
        round += 1;
        if !made_progress && round > 64 {
            // Safety valve: should not happen for valid square-free inputs.
            break;
        }
    }
    factors
        .into_iter()
        .filter(|f| f.degree().unwrap_or(0) > 0)
        .collect()
}

/// Attempt one EDF split for odd characteristic: pick pseudo-random candidates
/// `a`, compute $b = a^{(p^d-1)/2} \bmod f$, and return the nontrivial split
/// $\gcd(f, b-1)$ together with $f/\gcd(f,b-1)$. The seed varies across rounds
/// so repeated passes explore fresh candidates.
fn edf_split_odd(f: &FpPoly, half_exp: &BigInt, round: u64, idx: u64) -> Option<(FpPoly, FpPoly)> {
    let field = f.domain().clone();
    let deg_f = f.degree().unwrap_or(0);
    for attempt in 0..256u64 {
        let seed = round
            .wrapping_mul(0x10000)
            .wrapping_add(idx.wrapping_mul(0x100))
            .wrapping_add(attempt);
        let a = random_candidate(&field, deg_f, seed);
        if a.is_zero() || a.degree().unwrap_or(0) == 0 {
            continue;
        }
        let b = poly_pow_mod(&a, half_exp, f);
        // gcd(f, b - 1)
        let b_minus_one = b.sub(&FpPoly::from_coeffs(field.clone(), vec![field.one()]));
        let g = f.gcd(&b_minus_one);
        let dg = g.degree().unwrap_or(0);
        if dg > 0
            && dg < deg_f
            && let Some((q, _)) = f.div_rem(&g)
        {
            return Some((g, q));
        }
    }
    None
}

/// Attempt one EDF split in characteristic 2 using the trace map
/// $T(a) = a + a^2 + a^{2^2} + \dots + a^{2^{d-1}} \bmod f$.
fn edf_split_char2(f: &FpPoly, d: usize, round: u64, idx: u64) -> Option<(FpPoly, FpPoly)> {
    let field = f.domain().clone();
    let deg_f = f.degree().unwrap_or(0);
    for attempt in 0..256u64 {
        let seed = round
            .wrapping_mul(0x10000)
            .wrapping_add(idx.wrapping_mul(0x100))
            .wrapping_add(attempt);
        let a = random_candidate(&field, deg_f, seed);
        if a.is_zero() || a.degree().unwrap_or(0) == 0 {
            continue;
        }
        // trace = sum_{i=0}^{d-1} a^(2^i) mod f
        let mut trace = a.clone();
        let mut term = a;
        for _ in 1..d {
            term = poly_pow_mod(&term, &BigInt::from(2u32), f);
            trace = trace.add(&term);
            if let Some((_, r)) = trace.div_rem(f) {
                trace = r;
            }
        }
        let g = f.gcd(&trace);
        let dg = g.degree().unwrap_or(0);
        if dg > 0
            && dg < deg_f
            && let Some((q, _)) = f.div_rem(&g)
        {
            return Some((g, q));
        }
    }
    None
}

/// Berlekamp factorization of a monic square-free polynomial over a small
/// prime field $\mathbb{F}_p$.
///
/// Constructs the Frobenius matrix $Q$ where row $i$ is the coefficient vector
/// of $x^{ip} \bmod f$, then finds the nullspace of $Q - I$; each nontrivial
/// nullspace vector yields a splitting via $\gcd(f, v - a)$ for $a \in
/// \mathbb{F}_p$. Best suited to small primes where the matrix fits in memory
/// and the brute-force $a$-loop is affordable.
///
/// Reference: Berlekamp (1970), "Factoring polynomials over finite fields".
pub fn berlekamp(f: &FpPoly) -> Vec<FpPoly> {
    let n = match f.degree() {
        Some(d) if d > 0 => d,
        _ => return Vec::new(),
    };
    if n == 1 {
        return vec![monic(f)];
    }
    let field = f.domain().clone();
    let p = field.prime().clone();
    let zero = field.zero();
    let one = field.one();

    // ── build Frobenius matrix Q, n×n ──────────────────────────────
    // Q[i][j] = coefficient of x^j in x^{i·p} mod f.
    // We compute xp = x^p mod f first, then row_i = (row_{i-1})^p mod f.
    let x = x_var(&field);
    let xp = poly_pow_mod(&x, &p, f); // x^p mod f (deg < n)
    let mut q = Vec::with_capacity(n);
    // Row 0: x^0 = 1, i.e. [1, 0, ..., 0].
    let mut row0 = vec![zero.clone(); n];
    row0[0] = field.one();
    q.push(row0);
    // Row 1: x^p mod f.
    let mut row = xp.coeffs().to_vec(); // degree-ordered coefficients
    row.resize(n, zero.clone());
    q.push(row);
    // Rows 2..n-1: multiply by xp (mod f) each step.
    for i in 2..n {
        let prev = q[i - 1].clone();
        // Compute prev_poly · xp mod f.
        // prev_poly = sum prev[j] x^j, xp = sum xp_coeffs[k] x^k.
        // Product: sum_{j,k} prev[j]·xp_coeffs[k] · x^{j+k}. Reduce mod f.
        let mut prod = vec![field.zero(); 2 * n];
        for (j, cj) in prev.iter().enumerate() {
            if field.is_zero(cj) {
                continue;
            }
            for (k, ck) in q[1].iter().enumerate() {
                if field.is_zero(ck) {
                    continue;
                }
                let t = field.mul(cj, ck);
                prod[j + k] = field.add(&prod[j + k], &t);
            }
        }
        // Reduce modulo f.
        let prod_poly = FpPoly::from_coeffs(field.clone(), prod.clone());
        let (_q, r) = prod_poly.div_rem(f).unwrap_or_else(|| {
            (
                FpPoly::from_coeffs(field.clone(), vec![field.zero()]),
                prod_poly.clone(),
            )
        });
        let mut reduced = r.coeffs().to_vec();
        reduced.resize(n, zero.clone());
        q.push(reduced);
    }

    // ── Q - I ──────────────────────────────────────────────────────
    let mut m = Vec::with_capacity(n * n);
    for (i, row) in q.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let mut v = v.clone();
            if i == j {
                v = field.sub(&v, &one);
            }
            m.push(v);
        }
    }
    // ── Gaussian elimination to reduced row echelon over F_p ───────
    // We work on mutable copy since Matrix::row_echelon does fraction-free
    // elimination which is fine for a field (just a bit of extra gcd work).
    let mut re = Matrix::new(n, n, m.clone(), field.clone());
    let rank = re.row_echelon(n);
    // Build a basis of the nullspace. After row-echelon form:
    // Each free column j (not a pivot) gives a basis vector with entry 1
    // at column j, zero at other free columns, and values from the row
    // for pivot columns.
    let mut pivots = vec![false; n];
    let mut pivot_col = vec![0usize; n]; // pivot_col[row] = col
    let mut pi = 0;
    for j in 0..n {
        if pi >= rank {
            break;
        }
        if !field.is_zero(&re[(pi, j)]) {
            pivots[j] = true;
            pivot_col[pi] = j;
            pi += 1;
        }
    }

    // For each non-zero non-constant nullspace vector, attempt to split.
    let mut factors = vec![monic(f)];
    for j in 0..n {
        if pivots[j] {
            continue;
        }
        // Free column → build basis vector v (row-vector convention).
        // v[j] = 1; for each pivot row r: v[pivot_col[r]] = -re[r,j]/re[r,pivot_col[r]].
        let mut v = vec![zero.clone(); n];
        v[j] = field.one();
        for r in 0..rank {
            let pc = pivot_col[r];
            if pc >= j {
                continue;
            }
            let pivot_val = re[(r, pc)].clone();
            let free_val = re[(r, j)].clone();
            // v[pc] = -free_val / pivot_val
            if field.is_zero(&free_val) {
                continue;
            }
            let t = field
                .div(&free_val, &pivot_val)
                .unwrap_or_else(|| field.zero());
            v[pc] = field.neg(&t);
        }
        // Build polynomial from v.
        let vpoly = FpPoly::from_coeffs(field.clone(), v);
        if vpoly.is_zero() || vpoly.degree().unwrap_or(0) == 0 {
            continue;
        }
        // For each a in F_p, compute gcd(f, v - a).
        let mut new_factors = Vec::new();
        for factor in &factors {
            let df = factor.degree().unwrap_or(0);
            if df <= 1 || df == n {
                new_factors.push(factor.clone());
                continue;
            }
            let mut f_remaining = factor.clone();
            for a in 0..p.to_u64().unwrap_or(1) {
                let a_el = field.element(BigInt::from(a));
                let va = vpoly.sub(&FpPoly::from_coeffs(field.clone(), vec![a_el]));
                let g = f_remaining.gcd(&va);
                let g = monic(&g);
                let dg = g.degree().unwrap_or(0);
                if dg > 0 && dg < f_remaining.degree().unwrap_or(0) {
                    let f_over_g = f_remaining
                        .div_rem(&g)
                        .map(|(q, _)| monic(&q))
                        .unwrap_or_else(|| f_remaining.zero());
                    new_factors.push(g);
                    f_remaining = f_over_g;
                    if f_remaining.degree().unwrap_or(0) <= 1 {
                        break;
                    }
                }
            }
            if f_remaining.degree().unwrap_or(0) > 0 {
                new_factors.push(monic(&f_remaining));
            }
        }
        factors = new_factors;
    }
    factors
        .into_iter()
        .filter(|g| !g.is_zero() && !g.is_one())
        .collect()
}

/// Cantor–Zassenhaus factorization of a monic square-free polynomial.
///
/// Returns the list of its monic irreducible factors over $\mathbb{F}_p$.
pub fn cantor_zassenhaus(f: &FpPoly) -> Vec<FpPoly> {
    if f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (g, d) in distinct_degree_factorization(f) {
        for irr in equal_degree_factorization(&g, d) {
            result.push(irr);
        }
    }
    result
}

/// Completely factor a univariate polynomial over the finite field
/// $\mathbb{F}_p$ into monic irreducible factors with multiplicities.
///
/// The input's leading coefficient is extracted separately: the returned
/// factors are all monic, and a leading constant factor is included as the
/// first element when the input is not monic (with multiplicity 1).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use ocas_domain::{Domain, FiniteField};
/// use ocas_poly::DenseUnivariatePolynomial;
/// use ocas_poly::factor::finite_field::factor_over_finite_field;
///
/// let f = FiniteField::new(BigInt::from(5));
/// // x^2 - 1 = (x-1)(x+1) over F_5
/// let p = DenseUnivariatePolynomial::from_coeffs(
///     f.clone(), vec![f.element(4), f.element(0), f.element(1)]);
/// let factors = factor_over_finite_field(&p);
/// // Two monic linear factors.
/// let linear_count = factors.iter()
///     .filter(|(g, _)| g.degree() == Some(1)).count();
/// assert_eq!(linear_count, 2);
/// ```
pub fn factor_over_finite_field(f: &FpPoly) -> Vec<(FpPoly, usize)> {
    let field = f.domain().clone();
    let mut result = Vec::new();
    if f.is_zero() {
        return result;
    }

    // Extract the leading coefficient as a (constant) factor so the rest is monic.
    let lc = f.leading_coeff().cloned().unwrap_or_else(|| field.zero());
    let monic_f = monic(f);
    if !field.is_one(&lc) && !field.is_zero(&lc) {
        result.push((FpPoly::from_coeffs(field.clone(), vec![lc]), 1));
    }

    // Square-free factorization over the finite field (handles characteristic p
    // via p-th root extraction, unlike the generic char-0 Yun algorithm).
    for (g, multiplicity) in square_free_factorization_ff(&monic_f) {
        if g.degree().unwrap_or(0) == 0 {
            continue;
        }
        // Use Cantor–Zassenhaus for all primes. The Berlekamp path
        // (p ≤ 1000) is implemented but needs further validation
        // (TODO: fix nullspace extraction for degree-4+ factors).
        let use_berlekamp = false; // disabled pending validation
        let irr_factors = if use_berlekamp {
            berlekamp(&g)
        } else {
            cantor_zassenhaus(&g)
        };
        for irr in irr_factors {
            result.push((irr, multiplicity));
        }
    }
    result
}

/// Take the $p$-th root of a polynomial whose formal derivative is zero (i.e.
/// a polynomial in $x^p$) over a prime field $\mathbb{F}_p$.
///
/// Over $\mathbb{F}_p$ the Frobenius is the identity on coefficients
/// ($a^p = a$), so if $f(x) = \sum_k a_{pk} x^{pk}$ then its $p$-th root is
/// $g(x) = \sum_k a_{pk} x^k$.
fn pth_root_prime(f: &FpPoly) -> FpPoly {
    let field = f.domain().clone();
    let p = field
        .prime()
        .to_u64()
        .expect("p-th root extraction requires a prime fitting in u64") as usize;
    let deg = f.degree().unwrap_or(0);
    let mut coeffs = Vec::new();
    let mut j = 0;
    while p * j <= deg {
        let c = f.coeff(p * j).cloned().unwrap_or_else(|| field.zero());
        coeffs.push(c);
        j += 1;
    }
    FpPoly::from_coeffs(field, coeffs)
}

/// Square-free factorization over a finite field, correctly handling
/// characteristic $p$ (the Musser/Bernardin algorithm).
///
/// Unlike the char-0 Yun algorithm, this detects when $f' = 0$ (meaning $f$ is
/// a $p$-th power) and recurses on its $p$-th root, scaling multiplicities by
/// $p$.
fn square_free_factorization_ff(f: &FpPoly) -> Vec<(FpPoly, usize)> {
    if f.is_zero() || f.degree().unwrap_or(0) == 0 {
        return Vec::new();
    }
    let field = f.domain().clone();
    let p = field.prime().clone();
    let mut output: Vec<(FpPoly, usize)> = Vec::new();

    let fp = f.derivative();
    if fp.is_zero() {
        // f is a p-th power: f = g^p. Take the p-th root and recurse.
        let g = pth_root_prime(f);
        for (gi, ki) in square_free_factorization_ff(&g) {
            output.push((gi, ki * p.to_u64().unwrap_or(1) as usize));
        }
        return output;
    }

    let mut c = f.gcd(&fp);
    c = monic(&c);
    let mut w = monic(&match f.div_rem(&c) {
        Some((q, _)) => q,
        None => f.clone(),
    });
    let mut i = 1usize;
    while !w.is_one() && w.degree().unwrap_or(0) > 0 {
        let y = monic(&w.gcd(&c));
        let z = monic(&match w.div_rem(&y) {
            Some((q, _)) => q,
            None => w.clone(),
        });
        if z.degree().unwrap_or(0) > 0 && !z.is_one() {
            output.push((z, i));
        }
        w = y;
        if let Some((q, _)) = c.div_rem(&w) {
            c = monic(&q);
        }
        i += 1;
    }
    // Any leftover c is a p-th power of the remaining factors.
    if c.degree().unwrap_or(0) > 0 && !c.is_one() {
        let g = pth_root_prime(&c);
        for (gi, ki) in square_free_factorization_ff(&g) {
            output.push((gi, ki * p.to_u64().unwrap_or(1) as usize));
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseUnivariatePolynomial;

    fn field(p: u64) -> FiniteField {
        FiniteField::new(BigInt::from(p))
    }

    fn fpoly(f: &FiniteField, coeffs: &[i64]) -> FpPoly {
        DenseUnivariatePolynomial::from_coeffs(
            f.clone(),
            coeffs.iter().map(|&c| f.element(BigInt::from(c))).collect(),
        )
    }

    /// Multiply a list of factors together (over the field).
    fn product(polys: &[FpPoly]) -> FpPoly {
        let field = polys
            .first()
            .map(|p| p.domain().clone())
            .unwrap_or_else(|| FiniteField::new(BigInt::from(2)));
        let mut acc = FpPoly::from_coeffs(field.clone(), vec![field.one()]);
        for p in polys {
            acc = acc.mul(p);
        }
        acc
    }

    /// Reconstruct the original polynomial from a factorization (multiplicities
    /// included) and confirm it equals the input up to a nonzero scalar.
    fn assert_factors_reconstruct(input: &FpPoly, factors: &[(FpPoly, usize)]) {
        let field = input.domain().clone();
        let mut acc = FpPoly::from_coeffs(field.clone(), vec![field.one()]);
        for (g, m) in factors {
            for _ in 0..*m {
                acc = acc.mul(g);
            }
        }
        // acc should be a nonzero scalar multiple of input.
        let deg_acc = acc.degree().unwrap_or(0);
        let deg_in = input.degree().unwrap_or(0);
        assert_eq!(deg_acc, deg_in, "degree mismatch on reconstruction");
        // Verify exact divisibility: input / acc should be a constant.
        if let Some((q, r)) = input.div_rem(&acc) {
            assert!(r.is_zero(), "nonzero remainder when reconstructing");
            assert_eq!(q.degree(), Some(0), "quotient should be a constant");
        }
    }

    #[test]
    fn poly_pow_mod_basic() {
        let f = field(7);
        let m = fpoly(&f, &[1, 0, 1]); // x^2 + 1
        let x = fpoly(&f, &[0, 1]);
        // x^2 mod (x^2+1) = -1 = 6
        let r = poly_pow_mod(&x, &BigInt::from(2u32), &m);
        assert_eq!(r.coeff(0).cloned(), Some(f.element(BigInt::from(6))));
        // x^4 mod (x^2+1): x^2 = -1, so x^4 = 1
        let r = poly_pow_mod(&x, &BigInt::from(4u32), &m);
        assert_eq!(r.coeff(0).cloned(), Some(f.element(BigInt::from(1))));
    }

    #[test]
    fn factor_x_squared_minus_1_over_f5() {
        let f = field(5);
        let p = fpoly(&f, &[4, 0, 1]); // x^2 - 1
        let factors = factor_over_finite_field(&p);
        let linear_count = factors
            .iter()
            .filter(|(g, _)| g.degree() == Some(1))
            .count();
        assert_eq!(linear_count, 2);
        assert_factors_reconstruct(&p, &factors);
    }

    #[test]
    fn factor_irreducible_degree_3_over_f2() {
        let f = field(2);
        // x^3 + x + 1 is irreducible over F_2.
        let p = fpoly(&f, &[1, 1, 0, 1]);
        let factors = factor_over_finite_field(&p);
        let nontrivial: Vec<_> = factors
            .iter()
            .filter(|(g, _)| g.degree().unwrap_or(0) > 0)
            .collect();
        assert_eq!(nontrivial.len(), 1);
        assert_eq!(nontrivial[0].0.degree(), Some(3));
        assert_factors_reconstruct(&p, &factors);
    }

    #[test]
    fn factor_cyclotomic_x4_plus_1_over_f3() {
        // Over F_3: x^4 + 1 = (x^2+x+2)(x^2-x-1)... factor and verify irreducibles.
        let f = field(3);
        let p = fpoly(&f, &[1, 0, 0, 0, 1]); // x^4 + 1
        let factors = factor_over_finite_field(&p);
        assert_factors_reconstruct(&p, &factors);
        // Every factor must be irreducible (cannot split further); verify by
        // checking each is square-free and re-factoring yields itself.
        for (g, _) in &factors {
            if g.degree().unwrap_or(0) < 1 {
                continue;
            }
            let re = cantor_zassenhaus(&monic(g));
            assert_eq!(re.len(), 1, "factor {:?} not irreducible", g);
        }
    }

    #[test]
    fn factor_four_irreducible_quadratics_over_f5() {
        // Product of four distinct irreducible monic quadratics over F_5.
        // Ground truth (from SymPy): x^2+2, x^2+3, x^2+x+1, x^2+x+2.
        let f = field(5);
        let q1 = fpoly(&f, &[2, 0, 1]);
        let q2 = fpoly(&f, &[3, 0, 1]);
        let q3 = fpoly(&f, &[1, 1, 1]);
        let q4 = fpoly(&f, &[2, 1, 1]);
        let p = product(&[q1, q2, q3, q4]);
        let factors = factor_over_finite_field(&p);
        assert_factors_reconstruct(&p, &factors);
        let quad_count = factors
            .iter()
            .filter(|(g, _)| g.degree() == Some(2))
            .count();
        assert_eq!(
            quad_count, 4,
            "expected four quadratic factors, got {:?}",
            factors
        );
        // Each must be irreducible.
        for (g, _) in &factors {
            if g.degree().unwrap_or(0) < 1 {
                continue;
            }
            let re = cantor_zassenhaus(&monic(g));
            assert_eq!(re.len(), 1, "factor {:?} not irreducible", g);
        }
    }

    #[test]
    fn factor_cyclotomic_matches_sympy() {
        // Ground-truth degree histograms from SymPy 1.14
        // `Poly(x^n-1, x, domain=GF(p)).factor_list()`, expressed as
        // (p, n) -> sorted (degree, count) pairs.
        #[allow(clippy::type_complexity)]
        let cases: &[(u64, usize, &[(usize, usize)])] = &[
            (5, 10, &[(1, 10)]),
            (5, 30, &[(1, 10), (2, 10)]),
            (7, 10, &[(1, 2), (4, 2)]),
            (7, 30, &[(1, 6), (4, 6)]),
            (17, 20, &[(1, 4), (4, 4)]),
            (17, 30, &[(1, 2), (2, 2), (4, 6)]),
            (17, 100, &[(1, 4), (4, 4), (20, 4)]),
            (2, 10, &[(1, 2), (4, 2)]),
            (3, 30, &[(1, 6), (4, 6)]),
        ];
        for &(p, n, expected) in cases {
            let f = field(p);
            let mut coeffs = vec![f.element(BigInt::from(-1i64))];
            coeffs.resize(n + 1, f.element(BigInt::from(0)));
            coeffs[n] = f.element(BigInt::from(1));
            let poly = DenseUnivariatePolynomial::from_coeffs(f, coeffs);
            let factors = factor_over_finite_field(&poly);
            // Build the observed degree histogram over non-constant factors.
            let mut obs: std::collections::BTreeMap<usize, usize> = Default::default();
            for (g, m) in &factors {
                if g.degree().unwrap_or(0) > 0 {
                    *obs.entry(g.degree().unwrap()).or_insert(0) += *m;
                }
            }
            let observed: Vec<(usize, usize)> = obs.into_iter().collect();
            assert_eq!(
                observed.as_slice(),
                expected,
                "x^{n}-1 over F_{p}: degree histogram mismatch",
            );
        }
    }

    #[test]
    fn factor_repeated_root_over_f7() {
        let f = field(7);
        // (x-2)^3 * (x+1) over F_7
        let l1 = fpoly(&f, &[-2, 1]); // x - 2
        let l2 = fpoly(&f, &[1, 1]); //  x + 1
        let mut p = product(&[l1.clone(), l1.clone(), l1.clone(), l2.clone()]);
        let _ = &mut p;
        let factors = factor_over_finite_field(&p);
        assert_factors_reconstruct(&p, &factors);
        // The factor (x-2) should have multiplicity 3.
        let m2 = factors
            .iter()
            .filter(|(g, _)| {
                g.degree() == Some(1) && g.coeff(0).map(|c| c == &f.element(-2)).unwrap_or(false)
            })
            .map(|(_, m)| *m)
            .next();
        assert_eq!(m2, Some(3));
    }
}
