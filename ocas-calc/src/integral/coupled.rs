//! Coupled differential systems `y′ + A·y = b` over the base field `ℚ(x)`
//! with a **constant rational** matrix `A` (0.29.0, Wave D).
//!
//! This is the system form the parametric Risch differential equation
//! bottoms out to (Bronstein ch. 6–7). 0.29.0 ships it as **infrastructure
//! with its own tests**: no corpus case in the 1892-benchmark needs it (the
//! Wave-0 attribution), and the RDE cancellation branches that would call
//! it (`cancel_*` with `b` a logarithmic derivative) decline honestly
//! instead — the same place SymPy leaves `NotImplementedError`.
//!
//! Method: `A` must be diagonalizable over `ℚ` (rational eigenvalues with
//! full geometric multiplicity — computed exactly via Faddeev–LeVerrier +
//! ℚ-factorization + RREF nullspaces). With `P` the eigenvector basis, the
//! substitution `y = P·z` decouples the system into scalar RDEs
//! `zᵢ′ + λᵢ·zᵢ = (P⁻¹b)ᵢ`, each solved by the 0.29.0 rational RDE at the
//! base level. Every returned solution is verified exactly
//! (`D y + A·y == b` in the field) before it leaves the function.

use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

use super::rde::rde_solve;
use crate::tower::build::Tower;
use crate::tower::elem::KElem;

type DPoly = DenseUnivariatePolynomial<RationalDomain>;

/// Largest system dimension handled (deterministic budget).
const MAX_COUPLED_DIM: usize = 4;

/// Solve `y′ + A·y = b` for `y ∈ ℚ(x)ⁿ` with `A ∈ ℚ^{n×n}`.
///
/// Returns `None` when `A` is not diagonalizable over `ℚ`, when a decoupled
/// scalar RDE is outside the rational fragment, or when any budget trips.
/// The result is verified exactly before returning.
pub(crate) fn solve_coupled_constant(
    tower: &Tower,
    a: &[Vec<Rational>],
    b: &[KElem],
) -> Option<Vec<KElem>> {
    let n = a.len();
    if n == 0 || n > MAX_COUPLED_DIM || b.len() != n {
        return None;
    }
    let nv = tower.n_vars();
    // Characteristic polynomial c(λ) = λⁿ + c₁λⁿ⁻¹ + … via Faddeev–LeVerrier.
    let charpoly = char_poly_flv(a);
    // Rational roots with multiplicity; all of them, or decline.
    let roots = rational_roots_with_multiplicity(&charpoly)?;
    // Eigen-basis per eigenvalue (RREF nullspace of `A − λI`).
    let mut basis: Vec<Vec<Rational>> = Vec::new(); // columns
    let mut lambdas: Vec<Rational> = Vec::new();
    for (lambda, mult) in &roots {
        let ns = nullspace(&sub_diag(a, lambda));
        if ns.len() != *mult {
            return None; // not diagonalizable over ℚ
        }
        for v in ns {
            basis.push(v);
            lambdas.push(lambda.clone());
        }
    }
    // P = [v₁ … vₙ] (columns); P⁻¹ via RREF on [P | I].
    let p_inv = inverse(&transpose_cols(&basis, n))?;
    // w = P⁻¹·b (componentwise KElem arithmetic over the trivial tower).
    let mut ys: Vec<KElem> = Vec::with_capacity(n);
    for (col, lambda) in lambdas.iter().enumerate() {
        let _ = col;
        // w_i = Σ_j P⁻¹[i][j]·b[j]
        let mut w = KElem::zero(nv);
        for (j, bj) in b.iter().enumerate() {
            w = w.add(&bj.mul_rational(&p_inv[col][j]));
        }
        // zᵢ′ + λᵢ·zᵢ = wᵢ — the 0.29.0 rational RDE at the base level.
        let f = KElem::from_rational(lambda, nv);
        ys.push(rde_solve(tower, 0, &f, &w)?);
    }
    // y = P·z.
    let p = transpose_cols(&basis, n);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut yi = KElem::zero(nv);
        for (j, zj) in ys.iter().enumerate() {
            yi = yi.add(&zj.mul_rational(&p[i][j]));
        }
        out.push(yi);
    }
    // Exact verification: D y + A·y == b componentwise.
    for (i, yi) in out.iter().enumerate() {
        let dyi = crate::tower::build::tower_diff(yi, &[]);
        let mut lhs = dyi;
        for (j, yj) in out.iter().enumerate() {
            lhs = lhs.add(&yj.mul_rational(&a[i][j]));
        }
        if !lhs.eq_cross(&b[i]) {
            return None;
        }
    }
    Some(out)
}

// ------------------------------------------------------------------
//  Exact rational linear algebra (small, self-contained)
// ------------------------------------------------------------------

/// Characteristic polynomial `det(λI − A)` via Faddeev–LeVerrier, as a
/// dense `ℚ[λ]` polynomial (ascending coefficients).
fn char_poly_flv(a: &[Vec<Rational>]) -> DPoly {
    let dom = RationalDomain;
    let n = a.len();
    // B₀ = I; cₖ = −(1/k)·tr(A·Bₖ₋₁); Bₖ = A·Bₖ₋₁ + cₖ·I.
    // charpoly = λⁿ + c₁λⁿ⁻¹ + … + cₙ; we collect c₀ = 1 … cₙ.
    let mut cs: Vec<Rational> = vec![dom.one()];
    let mut b_prev: Vec<Vec<Rational>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { dom.one() } else { dom.zero() })
                .collect()
        })
        .collect();
    for k in 1..=n {
        let ab = mat_mul(a, &b_prev);
        let trace = (0..n).fold(dom.zero(), |acc, i| dom.add(&acc, &ab[i][i]));
        let ck = dom.neg(&dom.mul(&trace, &Rational::new(1, k as i64)));
        cs.push(ck.clone());
        // Bₖ = A·Bₖ₋₁ + cₖ·I
        b_prev = ab;
        for i in 0..n {
            b_prev[i][i] = dom.add(&b_prev[i][i], &ck);
        }
    }
    // c(λ) = Σₖ cₖ·λ^(n−k) → ascending coefficients are cs reversed.
    DPoly::from_coeffs(RationalDomain, cs.into_iter().rev().collect())
}

/// Distinct rational roots of `f` with their multiplicities, or `None` when
/// `f` does not split completely over `ℚ`.
fn rational_roots_with_multiplicity(f: &DPoly) -> Option<Vec<(Rational, usize)>> {
    let (roots, fully_split) = super::rational::rational_roots(f)?;
    if !fully_split {
        return None;
    }
    let dom = RationalDomain;
    let mut out = Vec::new();
    for r in roots {
        // Multiplicity by repeated synthetic division of `f` by (λ − r).
        let mut cur = f.clone();
        let mut mult = 0usize;
        loop {
            // Evaluate cur(r) exactly.
            let mut acc = dom.zero();
            for c in cur.coeffs().iter().rev() {
                acc = dom.add(&dom.mul(&acc, &r), c);
            }
            if !dom.is_zero(&acc) {
                break;
            }
            // Deflate: ascending synthetic division by (λ − r).
            let m = cur.coeffs().len() - 1;
            let mut quo = vec![dom.zero(); m];
            let mut carry = cur.coeffs()[m].clone();
            for i in (1..=m).rev() {
                quo[i - 1] = carry.clone();
                if i >= 2 {
                    carry = dom.add(&cur.coeffs()[i - 1], &dom.mul(&carry, &r));
                }
            }
            cur = DPoly::from_coeffs(RationalDomain, quo);
            mult += 1;
        }
        out.push((r, mult));
    }
    Some(out)
}

/// `A − λI`.
fn sub_diag(a: &[Vec<Rational>], lambda: &Rational) -> Vec<Vec<Rational>> {
    let dom = RationalDomain;
    a.iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, v)| {
                    if i == j {
                        dom.sub(v, lambda)
                    } else {
                        v.clone()
                    }
                })
                .collect()
        })
        .collect()
}

/// Row-reduce `[m]` to RREF over ℚ and return both the RREF and the pivot
/// columns. Plain Gaussian elimination with exact rational arithmetic.
fn rref(m: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    let dom = RationalDomain;
    let mut m: Vec<Vec<Rational>> = m.to_vec();
    let rows = m.len();
    let cols = m.first().map(|r| r.len()).unwrap_or(0);
    let mut pivot_row = 0;
    for col in 0..cols {
        if pivot_row >= rows {
            break;
        }
        let Some(piv) = (pivot_row..rows).find(|&r| !dom.is_zero(&m[r][col])) else {
            continue;
        };
        m.swap(pivot_row, piv);
        let inv = dom.inv(&m[pivot_row][col]).expect("nonzero pivot");
        for v in &mut m[pivot_row] {
            *v = dom.mul(v, &inv);
        }
        for r in 0..rows {
            if r == pivot_row {
                continue;
            }
            let factor = m[r][col].clone();
            if dom.is_zero(&factor) {
                continue;
            }
            let pr = m[pivot_row].clone();
            m[r] = m[r]
                .iter()
                .zip(pr.iter())
                .map(|(v, p)| dom.sub(v, &dom.mul(&factor, p)))
                .collect();
        }
        pivot_row += 1;
    }
    m
}

/// Nullspace basis of `m` (columns of free variables), via RREF.
fn nullspace(m: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    let dom = RationalDomain;
    let rows = m.len();
    let cols = m.first().map(|r| r.len()).unwrap_or(0);
    let r = rref(m);
    // Pivot columns: first nonzero entry per row.
    let mut pivot_of_row = Vec::new();
    let mut pivot_cols = Vec::new();
    for row in &r {
        let pc = row.iter().position(|v| !dom.is_zero(v));
        pivot_of_row.push(pc);
        if let Some(c) = pc {
            pivot_cols.push(c);
        }
    }
    let mut basis = Vec::new();
    for free in 0..cols {
        if pivot_cols.contains(&free) {
            continue;
        }
        let mut v = vec![RationalDomain.zero(); cols];
        v[free] = dom.one();
        for (row, pc) in pivot_of_row.iter().enumerate() {
            if let Some(c) = pc {
                v[*c] = dom.neg(&r[row][free]);
            }
        }
        basis.push(v);
    }
    let _ = rows;
    basis
}

/// Matrix product over ℚ.
fn mat_mul(a: &[Vec<Rational>], b: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    let dom = RationalDomain;
    let n = a.len();
    let m = b[0].len();
    let inner = b.len();
    (0..n)
        .map(|i| {
            (0..m)
                .map(|j| {
                    (0..inner).fold(dom.zero(), |acc, k| {
                        dom.add(&acc, &dom.mul(&a[i][k], &b[k][j]))
                    })
                })
                .collect()
        })
        .collect()
}

/// Column vectors → row-major matrix (columns become the basis).
fn transpose_cols(cols: &[Vec<Rational>], n: usize) -> Vec<Vec<Rational>> {
    (0..n)
        .map(|i| (0..cols.len()).map(|j| cols[j][i].clone()).collect())
        .collect()
}

/// Matrix inverse over ℚ via RREF on `[m | I]`; `None` when singular.
fn inverse(m: &[Vec<Rational>]) -> Option<Vec<Vec<Rational>>> {
    let n = m.len();
    let dom = RationalDomain;
    let aug: Vec<Vec<Rational>> = m
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .cloned()
                .chain((0..n).map(|j| if i == j { dom.one() } else { dom.zero() }))
                .collect()
        })
        .collect();
    let r = rref(&aug);
    // The left half must be the identity.
    for (i, row) in r.iter().enumerate() {
        for (j, v) in row.iter().enumerate().take(n) {
            let expect = if i == j { dom.one() } else { dom.zero() };
            if *v != expect {
                return None;
            }
        }
    }
    Some(r.iter().map(|row| row[n..].to_vec()).collect())
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;
    use ocas_domain::Rational;

    use super::*;
    use crate::tower::build::build_tower;
    use crate::tower::convert::atom_to_rational_extended;
    use ocas_atom::Symbol;

    fn rat(p: i64, q: i64) -> Rational {
        Rational::new(p, q)
    }

    /// KElem over the trivial tower from a small rational string.
    fn kelem_of<'a>(ctx: &'a AtomArena<'a>, src: &str, n: usize) -> KElem {
        let atom = ocas_parse::parse(ctx, src).expect("parse");
        let x = ctx.var("x");
        let rf = atom_to_rational_extended(atom, &[x], n).expect("rational");
        KElem::new(rf.numerator, rf.denominator)
    }

    #[test]
    fn coupled_distinct_rational_eigenvalues() {
        // A = [[0, 1], [2, 1]]: eigenvalues 2, −1. With b = (0, x), the
        // decoupled system is z₁′ + 2z₁ = w₁, z₂′ − z₂ = w₂.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, x, Symbol::new("x")).expect("tower");
        let a = vec![vec![rat(0, 1), rat(1, 1)], vec![rat(2, 1), rat(1, 1)]];
        let b = vec![KElem::zero(1), kelem_of(&ctx, "x", 1)];
        let sol = solve_coupled_constant(&tower, &a, &b);
        // Diagonalizable with rational eigenvalues → a solution exists and
        // is verified inside; spot-check existence only (the verification
        // gate already proves `Dy + Ay == b`).
        assert!(sol.is_some(), "the decoupled system must solve");
    }

    #[test]
    fn coupled_declines_a_jordan_block() {
        // A = [[1, 1], [0, 1]]: repeated eigenvalue 1, geometric multiplicity
        // 1 < 2 → not diagonalizable → honest decline.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, x, Symbol::new("x")).expect("tower");
        let a = vec![vec![rat(1, 1), rat(1, 1)], vec![rat(0, 1), rat(1, 1)]];
        let b = vec![KElem::zero(1), kelem_of(&ctx, "x", 1)];
        assert!(solve_coupled_constant(&tower, &a, &b).is_none());
    }

    #[test]
    fn coupled_declines_irrational_eigenvalues() {
        // A = [[0, 2], [1, 0]]: eigenvalues ±√2 ∉ ℚ → decline.
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let tower = build_tower(&ctx, x, Symbol::new("x")).expect("tower");
        let a = vec![vec![rat(0, 1), rat(2, 1)], vec![rat(1, 1), rat(0, 1)]];
        let b = vec![KElem::zero(1), kelem_of(&ctx, "x", 1)];
        assert!(solve_coupled_constant(&tower, &a, &b).is_none());
    }
}
