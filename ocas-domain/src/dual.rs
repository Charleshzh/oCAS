//! Hyper-dual numbers for forward automatic differentiation.
//!
//! A [`HyperDual<T>`] carries a scalar value together with a set of
//! derivative components indexed by multi-indices (tuples of non-negative
//! integer powers, one per differentiation variable). This is the runtime-
//! shaped variant: the component layout is built once into a [`DualShape`]
//! (shared cheaply via `Arc`) and all arithmetic consults a precomputed
//! multiplication table.
//!
//! Only polynomial/rational arithmetic is supported (add, sub, mul, div,
//! neg, integer powers). Transcendental functions (sin/exp/log) require a
//! real-coefficient trait and are out of scope for this module; pair it with
//! [`crate::Rational`] for exact rational AD or with a floating type for
//! numerical AD.
//!
//! # Example
//!
//! ```
//! use ocas_domain::dual::{HyperDual, new_first_order};
//! use ocas_domain::Rational;
//!
//! let shape = new_first_order::<Rational>(2); // track ∂/∂x₀ and ∂/∂x₁
//! // f(x₀, x₁) = x₀·x₁ at point (3, 5).
//! let x = HyperDual::variable(&shape, 0, Rational::new(3, 1));
//! let y = HyperDual::variable(&shape, 1, Rational::new(5, 1));
//! let f = x * y;
//! assert_eq!(f.value(), &Rational::new(15, 1));    // f = 15
//! assert_eq!(f.deriv(0), Some(&Rational::new(5, 1))); // ∂f/∂x₀ = x₁ = 5
//! assert_eq!(f.deriv(1), Some(&Rational::new(3, 1))); // ∂f/∂x₁ = x₀ = 3
//! ```

use std::sync::Arc;

use crate::Rational;

/// Coefficient trait for dual numbers: rational-like arithmetic with a
/// multiplicative inverse.
///
/// This is intentionally narrow — only what [`HyperDual`] needs. [`Rational`]
/// satisfies it via its `std::ops` implementations.
pub trait DualCoeff:
    Clone
    + PartialEq
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
    + std::ops::Neg<Output = Self>
    + std::ops::AddAssign
    + std::ops::MulAssign
{
    /// The additive identity.
    fn zero() -> Self;
    /// The multiplicative identity.
    fn one() -> Self;
}

impl DualCoeff for Rational {
    fn zero() -> Self {
        Rational::new(0, 1)
    }
    fn one() -> Self {
        Rational::new(1, 1)
    }
}

/// A component layout for hyper-dual numbers.
///
/// Components are stored as multi-indices — vectors of non-negative integer
/// powers, one entry per differentiation variable. The layout must be
/// *ancestor-closed*: if a component with multi-index `m` is present, every
/// multi-index dominated component-wise by `m` must also be present. This is
/// required for the multiplication table to be complete.
///
/// The first component is always the all-zeros multi-index (the scalar value).
/// Build a shape with [`DualShape::new`] or use a constructor like
/// [`new_first_order`].
#[derive(Debug, Clone)]
pub struct DualShape {
    /// Component multi-indices in canonical order. Index 0 is the value
    /// component (all zeros).
    components: Vec<Vec<usize>>,
    /// Multiplication table: `(a, b, c)` means component `a` times component
    /// `b` contributes to component `c`. Pairs where one operand is the value
    /// component (index 0) are excluded — those are handled by scalar scaling.
    mult_table: Vec<(usize, usize, usize)>,
}

impl DualShape {
    /// Build a shape from an ancestor-closed list of component multi-indices.
    ///
    /// The all-zeros multi-index is forced to position 0 regardless of input
    /// order. Returns `None` if the layout is not ancestor-closed.
    pub fn new(mut components: Vec<Vec<usize>>) -> Option<Self> {
        if components.is_empty() {
            return None;
        }
        // Normalise: determine variable count from the longest multi-index,
        // pad shorter ones with zeros, and move the zero component to front.
        let nvars = components.iter().map(|c| c.len()).max()?;
        for c in components.iter_mut() {
            c.resize(nvars, 0);
        }
        let zero = vec![0usize; nvars];
        // Find and move the zero component to index 0.
        let pos = components.iter().position(|c| *c == zero)?;
        if pos != 0 {
            components.swap(0, pos);
        }
        // Validate ancestor-closure: for every component, all strictly smaller
        // (component-wise) non-negative multi-indices must also be present.
        let contains = |m: &[usize]| -> bool {
            components
                .iter()
                .any(|c| c.len() == nvars && c.as_slice() == m)
        };
        for c in &components {
            if !ancestors_present(c, &contains) {
                return None;
            }
        }
        // Build multiplication table.
        let mult_table = build_mult_table(&components);
        Some(Self {
            components,
            mult_table,
        })
    }

    /// Number of components (value + derivative slots).
    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// Number of differentiation variables.
    pub fn n_vars(&self) -> usize {
        self.components.first().map(|c| c.len()).unwrap_or(0)
    }

    /// The component multi-indices.
    pub fn components(&self) -> &[Vec<usize>] {
        &self.components
    }

    /// Multiplication table accessor (mainly for inspection / tests).
    pub fn mult_table(&self) -> &[(usize, usize, usize)] {
        &self.mult_table
    }

    /// Look up the component index for a multi-index, if present.
    pub fn index_of(&self, multi_index: &[usize]) -> Option<usize> {
        self.components
            .iter()
            .position(|c| c.as_slice() == multi_index)
    }
}

/// Check every strictly-smaller ancestor of `m` is present in `contains`.
fn ancestors_present<F: Fn(&[usize]) -> bool>(m: &[usize], contains: &F) -> bool {
    // Generate every component-wise sub-multi-index strictly smaller than m.
    // We do this by decrementing coordinates one at a time and checking the
    // resulting (distinct) ancestors recursively; a small Cartesian product
    // is fine since shapes are small in practice.
    let nvars = m.len();
    // Total ancestor space = product of (m[i]+1) minus 1 (m itself).
    let mut counts: Vec<usize> = vec![0; nvars];
    loop {
        // Skip the all-equal-to-m case (counts == m).
        let is_self = counts.iter().zip(m.iter()).all(|(a, b)| a == b);
        if !is_self && !contains(&counts) {
            return false;
        }
        // Increment counts as a mixed-radix counter with limits m[i].
        let mut i = 0;
        while i < nvars {
            counts[i] += 1;
            if counts[i] > m[i] {
                counts[i] = 0;
                i += 1;
            } else {
                break;
            }
        }
        if i == nvars {
            break;
        }
    }
    true
}

/// Build the truncated multiplication table for a component layout.
///
/// For every ordered pair `(a, b)` with `a, b > 0`, if the component-wise sum
/// of their multi-indices is itself a component, record `(a, b, c)` where `c`
/// is the index of the sum. Pairs summing beyond the layout are dropped
/// (truncation).
fn build_mult_table(components: &[Vec<usize>]) -> Vec<(usize, usize, usize)> {
    let nvars = components.first().map(|c| c.len()).unwrap_or(0);
    let mut table = Vec::new();
    for a in 1..components.len() {
        for b in a..components.len() {
            let sum: Vec<usize> = (0..nvars)
                .map(|i| components[a][i] + components[b][i])
                .collect();
            if let Some(c) = components
                .iter()
                .position(|comp| comp.len() == nvars && comp.as_slice() == sum.as_slice())
            {
                table.push((a, b, c));
            }
        }
    }
    table
}

/// A hyper-dual number: a value plus derivative components laid out by a
/// shared [`DualShape`].
#[derive(Debug, Clone)]
pub struct HyperDual<T: DualCoeff> {
    /// Component values, indexed in lockstep with `shape.components()`.
    values: Vec<T>,
    shape: Arc<DualShape>,
}

impl<T: DualCoeff> HyperDual<T> {
    /// Build a hyper-dual from a full component vector. Length must match the
    /// shape's component count.
    pub fn from_values(shape: Arc<DualShape>, values: Vec<T>) -> Option<Self> {
        if values.len() != shape.n_components() {
            return None;
        }
        Some(Self { values, shape })
    }

    /// The scalar value component (component 0).
    pub fn value(&self) -> &T {
        &self.values[0]
    }

    /// Borrow all components in shape order.
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// The shared shape.
    pub fn shape(&self) -> &Arc<DualShape> {
        &self.shape
    }

    /// Return the derivative with respect to variable `i` (first-order
    /// component `[0..1_i..0]`), if that component is in the shape.
    pub fn deriv(&self, i: usize) -> Option<&T> {
        let mut idx = vec![0usize; self.shape.n_vars()];
        if i >= idx.len() {
            return None;
        }
        idx[i] = 1;
        self.shape.index_of(&idx).map(|pos| &self.values[pos])
    }

    /// A pure constant: value `c`, all derivative components zero.
    pub fn constant(shape: &Arc<DualShape>, c: T) -> Self {
        let mut values = vec![T::zero(); shape.n_components()];
        values[0] = c;
        Self {
            values,
            shape: shape.clone(),
        }
    }

    /// The independent variable number `i` set to value `c`: derivative
    /// component `[0..1_i..0]` is one, all others zero.
    pub fn variable(shape: &Arc<DualShape>, i: usize, c: T) -> Self {
        let mut d = Self::constant(shape, c);
        let mut idx = vec![0usize; shape.n_vars()];
        if i < idx.len() {
            idx[i] = 1;
            if let Some(pos) = shape.index_of(&idx) {
                d.values[pos] = T::one();
            }
        }
        d
    }

    /// Additive identity for this shape.
    pub fn zero(shape: &Arc<DualShape>) -> Self {
        Self::constant(shape, T::zero())
    }

    /// Multiplicative identity for this shape.
    pub fn one(shape: &Arc<DualShape>) -> Self {
        Self::constant(shape, T::one())
    }

    /// Multiplicative inverse via the geometric series
    /// `1/(v + ε) = (1/v)·Σ_{k≥0} (-ε/v)^k`, truncated to the shape's
    /// non-value components. Returns `None` if the value component is zero.
    pub fn inv(&self) -> Option<Self> {
        let v = self.value().clone();
        if v == T::zero() {
            return None;
        }
        let inv_v = T::one().clone() / v.clone();
        // (-ε/v) as a truncated vector; component 0 is zero because ε has no
        // scalar part. All powers of this vector therefore keep component 0 = 0,
        // so `mul_truncated`'s scalar cross terms contribute nothing here.
        let mut neg_ratio = vec![T::zero(); self.values.len()];
        for (k, slot) in neg_ratio.iter_mut().enumerate().skip(1) {
            *slot = (self.values[k].clone() / v.clone()).neg();
        }
        // S = Σ_{p≥0} (-ε/v)^p  (scalar 1 from p=0, plus higher-order terms).
        let mut s = vec![T::zero(); self.values.len()];
        s[0] = T::one();
        let mut current_power = neg_ratio.clone();
        loop {
            let mut changed = false;
            for k in 1..s.len() {
                if current_power[k] != T::zero() {
                    s[k] += current_power[k].clone();
                    changed = true;
                }
            }
            if !changed {
                break;
            }
            current_power = mul_truncated(&current_power, &neg_ratio, &self.shape);
        }
        // inv = inv_v · S (scalar multiply every component).
        let result_values: Vec<T> = s.iter().map(|c| inv_v.clone() * c.clone()).collect();
        Some(Self {
            values: result_values,
            shape: self.shape.clone(),
        })
    }
}

/// Pointwise addition; both operands must share the same shape.
fn assert_same_shape<T: DualCoeff>(a: &HyperDual<T>, b: &HyperDual<T>) {
    debug_assert_eq!(
        a.shape.n_components(),
        b.shape.n_components(),
        "HyperDual shape mismatch"
    );
}

impl<T: DualCoeff> std::ops::Add for HyperDual<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_same_shape(&self, &rhs);
        let values: Vec<T> = self
            .values
            .iter()
            .zip(rhs.values.iter())
            .map(|(a, b)| a.clone() + b.clone())
            .collect();
        Self {
            values,
            shape: self.shape,
        }
    }
}

impl<T: DualCoeff> std::ops::Sub for HyperDual<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        assert_same_shape(&self, &rhs);
        let values: Vec<T> = self
            .values
            .iter()
            .zip(rhs.values.iter())
            .map(|(a, b)| a.clone() - b.clone())
            .collect();
        Self {
            values,
            shape: self.shape,
        }
    }
}

impl<T: DualCoeff> std::ops::Neg for HyperDual<T> {
    type Output = Self;
    fn neg(self) -> Self {
        let values: Vec<T> = self.values.iter().map(|a| a.clone().neg()).collect();
        Self {
            values,
            shape: self.shape,
        }
    }
}

impl<T: DualCoeff> std::ops::Mul for HyperDual<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        assert_same_shape(&self, &rhs);
        // Closed-form: result[k] = self[k]·rhs[0] + self[0]·rhs[k] + Σ table.
        let mut result = vec![T::zero(); self.values.len()];
        let sv = self.values[0].clone();
        let rv = rhs.values[0].clone();
        // Scalar cross terms: self[k]*rv + sv*rhs[k].
        for (k, slot) in result.iter_mut().enumerate().skip(1) {
            *slot = self.values[k].clone() * rv.clone();
            *slot += sv.clone() * rhs.values[k].clone();
        }
        // Value component.
        result[0] = sv * rv;
        // Table contributions (operand pairs both > 0).
        for &(a, b, c) in self.shape.mult_table() {
            // The table stores (a, b) with a <= b; both orderings contribute.
            result[c] += self.values[a].clone() * rhs.values[b].clone();
            if a != b {
                result[c] += self.values[b].clone() * rhs.values[a].clone();
            }
        }
        Self {
            values: result,
            shape: self.shape,
        }
    }
}

impl<T: DualCoeff> std::ops::Div for HyperDual<T> {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        let inv = rhs.inv().expect("division by zero-valued HyperDual");
        self * inv
    }
}

/// Truncated component-wise product used by `inv`'s power iteration.
///
/// Computes `self ⊗ rhs` keeping only components present in the shape.
/// Correct in general (used by `mul` via the closed-form there); `inv` calls
/// it only on operands whose scalar component is zero, so the scalar cross
/// terms below are harmless there.
fn mul_truncated<T: DualCoeff>(self_v: &[T], rhs_v: &[T], shape: &DualShape) -> Vec<T> {
    let mut result = vec![T::zero(); self_v.len()];
    let sv = &self_v[0];
    let rv = &rhs_v[0];
    result[0] = sv.clone() * rv.clone();
    for k in 1..result.len() {
        result[k] = self_v[k].clone() * rv.clone();
        result[k] += sv.clone() * rhs_v[k].clone();
    }
    for &(a, b, c) in shape.mult_table() {
        result[c] += self_v[a].clone() * rhs_v[b].clone();
        if a != b {
            result[c] += self_v[b].clone() * rhs_v[a].clone();
        }
    }
    result
}

/// Build a first-order shape tracking one derivative per variable.
///
/// Components are `[0..0]` (value) and `[0..1_i..0]` for each variable `i` in
/// `0..nvars`. Higher-order or mixed partials are not tracked.
pub fn new_first_order<T: DualCoeff>(nvars: usize) -> Arc<DualShape> {
    let mut components = vec![vec![0usize; nvars]];
    for i in 0..nvars {
        let mut m = vec![0usize; nvars];
        m[i] = 1;
        components.push(m);
    }
    Arc::new(DualShape::new(components).expect("first-order shape is ancestor-closed"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Rational;

    fn r(n: i64, d: i64) -> Rational {
        Rational::new(n, d)
    }

    #[test]
    fn first_order_shape_layout() {
        let shape = new_first_order::<Rational>(2);
        assert_eq!(shape.n_components(), 3);
        assert_eq!(shape.n_vars(), 2);
        assert_eq!(shape.components()[0], vec![0, 0]);
        assert_eq!(shape.components()[1], vec![1, 0]);
        assert_eq!(shape.components()[2], vec![0, 1]);
        // First-order: no second-order components, so the mult table is empty
        // (no pair of non-value components sums to a present component).
        assert!(shape.mult_table().is_empty());
    }

    #[test]
    fn variable_and_constant() {
        let shape = new_first_order::<Rational>(2);
        let x = HyperDual::variable(&shape, 0, r(3, 1));
        assert_eq!(x.value(), &r(3, 1));
        assert_eq!(x.deriv(0), Some(&r(1, 1)));
        assert_eq!(x.deriv(1), Some(&r(0, 1)));
        let c = HyperDual::constant(&shape, r(7, 1));
        assert_eq!(c.deriv(0), Some(&r(0, 1)));
    }

    #[test]
    fn product_of_two_variables() {
        let shape = new_first_order::<Rational>(2);
        let x = HyperDual::variable(&shape, 0, r(3, 1));
        let y = HyperDual::variable(&shape, 1, r(5, 1));
        let f = x * y;
        assert_eq!(f.value(), &r(15, 1));
        assert_eq!(f.deriv(0), Some(&r(5, 1))); // ∂(xy)/∂x = y
        assert_eq!(f.deriv(1), Some(&r(3, 1))); // ∂(xy)/∂y = x
    }

    #[test]
    fn sum_and_difference() {
        let shape = new_first_order::<Rational>(2);
        let x = HyperDual::variable(&shape, 0, r(2, 1));
        let y = HyperDual::variable(&shape, 1, r(7, 1));
        let s = x.clone() + y.clone();
        assert_eq!(s.value(), &r(9, 1));
        assert_eq!(s.deriv(0), Some(&r(1, 1)));
        assert_eq!(s.deriv(1), Some(&r(1, 1)));
        let d = x - y;
        assert_eq!(d.value(), &r(-5, 1));
        assert_eq!(d.deriv(1), Some(&r(-1, 1)));
    }

    #[test]
    fn reciprocal_of_constant() {
        let shape = new_first_order::<Rational>(1);
        let c = HyperDual::constant(&shape, r(4, 1));
        let inv = c.inv().unwrap();
        assert_eq!(inv.value(), &r(1, 4));
        assert_eq!(inv.deriv(0), Some(&r(0, 1)));
    }

    #[test]
    fn quotient_derivatives() {
        let shape = new_first_order::<Rational>(2);
        let x = HyperDual::variable(&shape, 0, r(6, 1));
        let y = HyperDual::variable(&shape, 1, r(3, 1));
        // f = x/y; ∂f/∂x = 1/y = 1/3; ∂f/∂y = -x/y² = -6/9 = -2/3.
        let f = x / y;
        assert_eq!(f.value(), &r(2, 1));
        assert_eq!(f.deriv(0), Some(&r(1, 3)));
        assert_eq!(f.deriv(1), Some(&r(-2, 3)));
    }

    #[test]
    fn reciprocal_of_variable_gives_correct_derivative() {
        let shape = new_first_order::<Rational>(1);
        let x = HyperDual::variable(&shape, 0, r(5, 1));
        let inv = x.inv().unwrap();
        // d/dx (1/x) = -1/x² = -1/25.
        assert_eq!(inv.value(), &r(1, 5));
        assert_eq!(inv.deriv(0), Some(&r(-1, 25)));
    }

    #[test]
    fn power_of_variable() {
        let shape = new_first_order::<Rational>(1);
        let x = HyperDual::variable(&shape, 0, r(3, 1));
        // x^3 by repeated multiplication.
        let x2 = x.clone() * x.clone();
        let x3 = x2 * x;
        assert_eq!(x3.value(), &r(27, 1));
        // d/dx(x^3) = 3x^2 = 27.
        assert_eq!(x3.deriv(0), Some(&r(27, 1)));
    }

    #[test]
    fn three_variable_product_derivatives() {
        // The headline proptest target: ∂(xyz)/∂x = yz, etc.
        let shape = new_first_order::<Rational>(3);
        let x = HyperDual::variable(&shape, 0, r(2, 1));
        let y = HyperDual::variable(&shape, 1, r(3, 1));
        let z = HyperDual::variable(&shape, 2, r(5, 1));
        let f = x * y * z;
        assert_eq!(f.value(), &r(30, 1));
        assert_eq!(f.deriv(0), Some(&r(15, 1))); // yz = 15
        assert_eq!(f.deriv(1), Some(&r(10, 1))); // xz = 10
        assert_eq!(f.deriv(2), Some(&r(6, 1))); // xy = 6
    }

    #[test]
    fn shape_rejects_non_ancestor_closed() {
        // Components {[0], [2]}: missing [1], so not ancestor-closed.
        let components = vec![vec![0], vec![2]];
        assert!(DualShape::new(components).is_none());
    }

    #[test]
    fn shape_accepts_second_order() {
        // {[0], [1], [2]}: ancestor-closed, tracks value + first + second deriv.
        let components = vec![vec![0], vec![1], vec![2]];
        let shape = DualShape::new(components).unwrap();
        assert_eq!(shape.n_components(), 3);
        // Mult table: [1]*[1] -> [2].
        assert_eq!(shape.mult_table(), &[(1usize, 1usize, 2usize)]);
    }
}
