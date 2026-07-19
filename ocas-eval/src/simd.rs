//! SIMD vectorized expression evaluation.
//!
//! The [`VectorEvaluator`] evaluates expressions on batches of input
//! vectors using `pulp` for portable SIMD computation. The SIMD width
//! is detected at runtime (SSE2/AVX2/AVX-512). [`VectorEvaluatorF32`]
//! provides the single-precision variant with twice the lane count.
//!
//! Enabled with the `simd` feature flag.

use pulp::Arch;

use crate::error::{EvaluationError, Result};
use crate::instruction::{BuiltinOp, Instr};

/// A scalar lane type for vectorized evaluation.
///
/// Implemented for `f64` and `f32`; abstracts the arithmetic needed by
/// the chunk engine so both precisions share one implementation.
pub(crate) trait Lane: Copy {
    const ZERO: Self;
    const ONE: Self;
    /// Maximum SIMD lane count for this type (AVX-512 width).
    const MAX_LANES: usize;

    fn add(self, other: Self) -> Self;
    fn mul(self, other: Self) -> Self;
    fn recip(self) -> Self;
    fn powi(self, exp: i32) -> Self;
    fn powf(self, exp: Self) -> Self;
    fn builtin(op: BuiltinOp, x: Self) -> Self;
}

impl Lane for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const MAX_LANES: usize = 8;

    fn add(self, other: Self) -> Self {
        self + other
    }
    fn mul(self, other: Self) -> Self {
        self * other
    }
    fn recip(self) -> Self {
        1.0 / self
    }
    fn powi(self, exp: i32) -> Self {
        self.powi(exp)
    }
    fn powf(self, exp: Self) -> Self {
        self.powf(exp)
    }
    fn builtin(op: BuiltinOp, x: Self) -> Self {
        match op {
            BuiltinOp::Sin => x.sin(),
            BuiltinOp::Cos => x.cos(),
            BuiltinOp::Tan => x.tan(),
            BuiltinOp::Sec => 1.0 / x.cos(),
            BuiltinOp::Csc => 1.0 / x.sin(),
            BuiltinOp::Cot => 1.0 / x.tan(),
            BuiltinOp::Exp => x.exp(),
            BuiltinOp::Log => x.ln(),
            BuiltinOp::Sqrt => x.sqrt(),
            BuiltinOp::Abs => x.abs(),
        }
    }
}

impl Lane for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const MAX_LANES: usize = 16;

    fn add(self, other: Self) -> Self {
        self + other
    }
    fn mul(self, other: Self) -> Self {
        self * other
    }
    fn recip(self) -> Self {
        1.0 / self
    }
    fn powi(self, exp: i32) -> Self {
        self.powi(exp)
    }
    fn powf(self, exp: Self) -> Self {
        self.powf(exp)
    }
    fn builtin(op: BuiltinOp, x: Self) -> Self {
        match op {
            BuiltinOp::Sin => x.sin(),
            BuiltinOp::Cos => x.cos(),
            BuiltinOp::Tan => x.tan(),
            BuiltinOp::Sec => 1.0 / x.cos(),
            BuiltinOp::Csc => 1.0 / x.sin(),
            BuiltinOp::Cot => 1.0 / x.tan(),
            BuiltinOp::Exp => x.exp(),
            BuiltinOp::Log => x.ln(),
            BuiltinOp::Sqrt => x.sqrt(),
            BuiltinOp::Abs => x.abs(),
        }
    }
}

/// Detect the SIMD lane count for lane type `L` at runtime.
fn detect_lanes<L: Lane>() -> usize {
    if cfg!(target_arch = "x86_64") {
        if std::is_x86_feature_detected!("avx512f") {
            L::MAX_LANES
        } else if std::is_x86_feature_detected!("avx2") {
            L::MAX_LANES / 2
        } else {
            L::MAX_LANES / 4 // SSE2 is baseline on x86_64
        }
    } else {
        1
    }
}

/// Generic batch evaluator over lane type `L`.
struct GenEvaluator<L: Lane> {
    instructions: Vec<Instr>,
    param_count: usize,
    #[allow(dead_code)]
    const_count: usize,
    stack_size: usize,
    result_indices: Vec<usize>,
    constants: Vec<L>,
}

/// Batch evaluator using SIMD for vectorized expression computation.
///
/// Processes input vectors in chunks whose width is detected at runtime,
/// falling back to scalar evaluation for any remainder.
pub struct VectorEvaluator {
    inner: GenEvaluator<f64>,
}

/// Single-precision batch evaluator (twice the lane count of
/// [`VectorEvaluator`] on the same hardware).
pub struct VectorEvaluatorF32 {
    inner: GenEvaluator<f32>,
}

impl VectorEvaluator {
    /// Create a vector evaluator from compiled instruction data.
    #[allow(dead_code)]
    pub(crate) fn new(
        instructions: Vec<Instr>,
        param_count: usize,
        const_count: usize,
        stack_size: usize,
        result_indices: Vec<usize>,
        constants: Vec<f64>,
    ) -> Self {
        Self {
            inner: GenEvaluator {
                instructions,
                param_count,
                const_count,
                stack_size,
                result_indices,
                constants,
            },
        }
    }

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.inner.param_count
    }

    /// Return the number of results produced per input row.
    pub fn result_count(&self) -> usize {
        self.inner.result_indices.len()
    }

    /// Evaluate the expression on a batch of inputs.
    ///
    /// Each element of `params` is a vector of input values for one
    /// parameter. All input vectors must have the same length.
    ///
    /// Returns one vector per result slot; each has the same length as
    /// the input vectors (row-major: `results[output][row]`).
    pub fn evaluate(&self, params: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
        self.inner.evaluate(params)
    }
}

impl VectorEvaluatorF32 {
    /// Create a single-precision vector evaluator from compiled
    /// instruction data. Constants are narrowed from f64 to f32.
    #[allow(dead_code)]
    pub(crate) fn new(
        instructions: Vec<Instr>,
        param_count: usize,
        const_count: usize,
        stack_size: usize,
        result_indices: Vec<usize>,
        constants: Vec<f32>,
    ) -> Self {
        Self {
            inner: GenEvaluator {
                instructions,
                param_count,
                const_count,
                stack_size,
                result_indices,
                constants,
            },
        }
    }

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.inner.param_count
    }

    /// Return the number of results produced per input row.
    pub fn result_count(&self) -> usize {
        self.inner.result_indices.len()
    }

    /// Evaluate the expression on a batch of single-precision inputs.
    /// Same contract as [`VectorEvaluator::evaluate`].
    pub fn evaluate(&self, params: &[Vec<f32>]) -> Result<Vec<Vec<f32>>> {
        self.inner.evaluate(params)
    }
}

impl<L: Lane> GenEvaluator<L> {
    /// Evaluate the expression on a batch of inputs.
    fn evaluate(&self, params: &[Vec<L>]) -> Result<Vec<Vec<L>>> {
        if params.is_empty() && self.param_count > 0 {
            return Err(EvaluationError::WrongArity {
                name: "batch".into(),
                expected: self.param_count,
                got: 0,
            });
        }
        if params.len() != self.param_count {
            return Err(EvaluationError::WrongArity {
                name: "batch".into(),
                expected: self.param_count,
                got: params.len(),
            });
        }

        let batch_size = if self.param_count > 0 {
            params[0].len()
        } else {
            1
        };

        for p in params.iter() {
            if p.len() != batch_size {
                return Err(EvaluationError::TypeMismatch {
                    expected: format!("vector of length {batch_size}"),
                    found: format!("vector of length {}", p.len()),
                });
            }
        }

        let arch = Arch::new();
        let mut results = vec![vec![L::ZERO; batch_size]; self.result_indices.len()];

        // Determine the SIMD lane count from CPU features and lane width.
        let lanes = detect_lanes::<L>();

        let chunks = batch_size / lanes;
        let remainder_start = chunks * lanes;

        if lanes > 1 {
            // Dispatch once for all SIMD chunks. The closure captures
            // everything by reference; `arch.dispatch` selects the best
            // available SIMD width at runtime.
            let instructions = &self.instructions;
            let constants = &self.constants;
            let param_count = self.param_count;
            let stack_size = self.stack_size;
            let result_indices = &self.result_indices;

            arch.dispatch(|| {
                eval_simd_chunks(
                    instructions,
                    constants,
                    param_count,
                    stack_size,
                    result_indices,
                    params,
                    &mut results,
                    chunks,
                    lanes,
                );
            });
        }

        // Scalar fallback for remainder elements (and full batch if lanes == 1).
        let scalar_start = if lanes > 1 { remainder_start } else { 0 };
        for i in scalar_start..batch_size {
            let scalar_params: Vec<L> = params.iter().map(|p| p[i]).collect();
            let scalar_results = self.evaluate_scalar(&scalar_params)?;
            for (ri, &res_idx) in self.result_indices.iter().enumerate() {
                results[ri][i] = scalar_results[res_idx];
            }
        }

        Ok(results)
    }

    /// Scalar evaluation for a single set of parameters.
    fn evaluate_scalar(&self, params: &[L]) -> Result<Vec<L>> {
        let mut stack: Vec<L> = vec![L::ZERO; self.stack_size];

        for (i, p) in params.iter().enumerate() {
            stack[i] = *p;
        }

        for (i, c) in self.constants.iter().enumerate() {
            stack[self.param_count + i] = *c;
        }

        for instr in &self.instructions {
            match instr {
                Instr::Add { dst, srcs } => {
                    let mut sum = stack[srcs[0]];
                    for idx in &srcs[1..] {
                        sum = sum.add(stack[*idx]);
                    }
                    stack[*dst] = sum;
                }
                Instr::Mul { dst, srcs } => {
                    let mut prod = stack[srcs[0]];
                    for idx in &srcs[1..] {
                        prod = prod.mul(stack[*idx]);
                    }
                    stack[*dst] = prod;
                }
                Instr::Pow { dst, base, exp } => {
                    stack[*dst] = stack[*base].powi(*exp as i32);
                }
                Instr::Powf { dst, base, exp } => {
                    stack[*dst] = stack[*base].powf(stack[*exp]);
                }
                Instr::BuiltinOp { dst, op, src } => {
                    stack[*dst] = L::builtin(*op, stack[*src]);
                }
                Instr::ExternalFun { .. } => {
                    return Err(EvaluationError::JitCompilationError {
                        message: "external functions not supported in SIMD".into(),
                    });
                }
                Instr::Copy { dst, src } => {
                    stack[*dst] = stack[*src];
                }
            }
        }

        Ok(stack)
    }
}

/// Evaluate SIMD chunks using a manually-unrolled fixed-width stack.
///
/// Each stack slot holds an `[L; 16]` array where only `lanes` elements
/// are active (depending on CPU features and lane width). Arithmetic is
/// performed lane-by-lane; transcendental functions fall back to scalar
/// per lane.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
#[inline]
fn eval_simd_chunks<L: Lane>(
    instructions: &[Instr],
    constants: &[L],
    param_count: usize,
    stack_size: usize,
    result_indices: &[usize],
    params: &[Vec<L>],
    out_results: &mut [Vec<L>],
    chunks: usize,
    lanes: usize,
) {
    // Pre-allocate stack buffer once, reuse across all chunks.
    let mut stack: Vec<[L; 16]> = vec![[L::ZERO; 16]; stack_size];

    for chunk in 0..chunks {
        let base = chunk * lanes;

        // Clear stack for this chunk (reuse allocation).
        stack.fill([L::ZERO; 16]);

        // Load parameters
        for (pi, p) in params.iter().enumerate() {
            stack[pi][..lanes].copy_from_slice(&p[base..(base + lanes)]);
        }

        // Load constants (broadcast to all lanes)
        for (ci, c) in constants.iter().enumerate() {
            for l in 0..lanes {
                stack[param_count + ci][l] = *c;
            }
        }

        // Execute instructions
        for instr in instructions {
            match instr {
                Instr::Add { dst, srcs } => {
                    let mut sum = stack[srcs[0]];
                    for &idx in &srcs[1..] {
                        let other = stack[idx];
                        for l in 0..lanes {
                            sum[l] = sum[l].add(other[l]);
                        }
                    }
                    stack[*dst] = sum;
                }
                Instr::Mul { dst, srcs } => {
                    let mut prod = stack[srcs[0]];
                    for &idx in &srcs[1..] {
                        let other = stack[idx];
                        for l in 0..lanes {
                            prod[l] = prod[l].mul(other[l]);
                        }
                    }
                    stack[*dst] = prod;
                }
                Instr::Pow {
                    dst,
                    base: base_idx,
                    exp,
                } => {
                    let base_val = stack[*base_idx];
                    let mut result = [L::ONE; 16];
                    let mut e = exp.abs();
                    let mut b = base_val;
                    while e > 0 {
                        if e & 1 == 1 {
                            for l in 0..lanes {
                                result[l] = result[l].mul(b[l]);
                            }
                        }
                        let mut new_b = [L::ZERO; 16];
                        for l in 0..lanes {
                            new_b[l] = b[l].mul(b[l]);
                        }
                        b = new_b;
                        e >>= 1;
                    }
                    if *exp < 0 {
                        for l in 0..lanes {
                            result[l] = result[l].recip();
                        }
                    }
                    stack[*dst] = result;
                }
                Instr::Powf {
                    dst,
                    base: base_idx,
                    exp: exp_idx,
                } => {
                    let base_val = stack[*base_idx];
                    let exp_val = stack[*exp_idx];
                    let mut result = [L::ZERO; 16];
                    for l in 0..lanes {
                        result[l] = base_val[l].powf(exp_val[l]);
                    }
                    stack[*dst] = result;
                }
                Instr::BuiltinOp { dst, op, src } => {
                    let arr = stack[*src];
                    let mut result = [L::ZERO; 16];
                    for l in 0..lanes {
                        result[l] = L::builtin(*op, arr[l]);
                    }
                    stack[*dst] = result;
                }
                Instr::ExternalFun { .. } => {
                    // Not supported in SIMD mode; caller should check.
                }
                Instr::Copy { dst, src } => {
                    stack[*dst] = stack[*src];
                }
            }
        }

        // Write results for every result index.
        for (ri, &res_idx) in result_indices.iter().enumerate() {
            out_results[ri][base..(base + lanes)].copy_from_slice(&stack[res_idx][..lanes]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_add_constant() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![1.0f64];
        let eval = VectorEvaluator::new(instructions, 1, 1, 3, vec![2], constants);

        let params = vec![vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]];
        let results = eval.evaluate(&params).unwrap();

        assert_eq!(results.len(), 1);
        for (i, &r) in results[0].iter().enumerate() {
            assert!((r - (i as f64 + 1.0)).abs() < 1e-10, "mismatch at {i}");
        }
    }

    #[test]
    fn simd_mul() {
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 0],
        }];
        let eval = VectorEvaluator::new(instructions, 1, 0, 3, vec![2], vec![]);

        let params = vec![vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]];
        let results = eval.evaluate(&params).unwrap();

        for (i, &r) in results[0].iter().enumerate() {
            let expected = (i as f64) * (i as f64);
            assert!((r - expected).abs() < 1e-10, "mismatch at {i}");
        }
    }

    #[test]
    fn simd_remainder() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![10.0f64];
        let eval = VectorEvaluator::new(instructions, 1, 1, 3, vec![2], constants);

        let params = vec![vec![1.0, 2.0, 3.0]];
        let results = eval.evaluate(&params).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 3);
        assert!((results[0][0] - 11.0).abs() < 1e-10);
        assert!((results[0][1] - 12.0).abs() < 1e-10);
        assert!((results[0][2] - 13.0).abs() < 1e-10);
    }

    #[test]
    fn simd_powi() {
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: 3,
        }];
        let eval = VectorEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);

        let params = vec![vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]];
        let results = eval.evaluate(&params).unwrap();

        for (i, &r) in results[0].iter().enumerate() {
            let expected = (i as f64).powi(3);
            assert!((r - expected).abs() < 1e-10, "mismatch at {i}");
        }
    }

    #[test]
    fn simd_multi_output() {
        // params: [x]; outputs: [x + 1, x * x]
        // stack: 0=x, 1=const(1.0), 2=add, 3=mul
        let instructions = vec![
            Instr::Add {
                dst: 2,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 3,
                srcs: vec![0, 0],
            },
        ];
        let eval = VectorEvaluator::new(instructions, 1, 1, 4, vec![2, 3], vec![1.0]);

        let params = vec![vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]];
        let results = eval.evaluate(&params).unwrap();

        assert_eq!(results.len(), 2);
        for i in 0..10 {
            let x = i as f64;
            assert!((results[0][i] - (x + 1.0)).abs() < 1e-10, "out0 at {i}");
            assert!((results[1][i] - (x * x)).abs() < 1e-10, "out1 at {i}");
        }
    }

    #[test]
    fn simd_f32_basic() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let eval = VectorEvaluatorF32::new(instructions, 1, 1, 3, vec![2], vec![1.0]);

        let params = vec![(0..20).map(|i| i as f32).collect::<Vec<_>>()];
        let results = eval.evaluate(&params).unwrap();

        assert_eq!(results.len(), 1);
        for (i, &r) in results[0].iter().enumerate() {
            assert!((r - (i as f32 + 1.0)).abs() < 1e-6, "mismatch at {i}");
        }
    }

    #[test]
    fn simd_f32_multi_output() {
        // outputs: [x^3, sin(x)]
        let instructions = vec![
            Instr::Pow {
                dst: 1,
                base: 0,
                exp: 3,
            },
            Instr::BuiltinOp {
                dst: 2,
                op: BuiltinOp::Sin,
                src: 0,
            },
        ];
        let eval = VectorEvaluatorF32::new(instructions, 1, 0, 3, vec![1, 2], vec![]);

        let params = vec![(0..20).map(|i| i as f32 * 0.1).collect::<Vec<_>>()];
        let results = eval.evaluate(&params).unwrap();

        assert_eq!(results.len(), 2);
        for (i, &x) in params[0].iter().enumerate() {
            assert!((results[0][i] - x.powi(3)).abs() < 1e-4, "pow at {i}");
            assert!((results[1][i] - x.sin()).abs() < 1e-5, "sin at {i}");
        }
    }

    #[test]
    fn simd_f32_matches_f64_within_precision() {
        // Same expression through f64 and f32 paths must agree to f32 precision
        let instructions = vec![
            Instr::Pow {
                dst: 2,
                base: 0,
                exp: 2,
            },
            Instr::Add {
                dst: 3,
                srcs: vec![2, 1],
            },
        ];
        let eval64 = VectorEvaluator::new(instructions.clone(), 1, 1, 4, vec![3], vec![1.0]);
        let eval32 = VectorEvaluatorF32::new(instructions, 1, 1, 4, vec![3], vec![1.0]);

        let inputs: Vec<f64> = (0..33).map(|i| i as f64 * 0.25).collect();
        let r64 = eval64.evaluate(std::slice::from_ref(&inputs)).unwrap();
        let inputs32: Vec<f32> = inputs.iter().map(|&x| x as f32).collect();
        let r32 = eval32.evaluate(&[inputs32]).unwrap();

        for i in 0..33 {
            let diff = (r32[0][i] as f64 - r64[0][i]).abs();
            let rel = diff / r64[0][i].abs().max(1.0);
            assert!(
                rel < 1e-5,
                "row {i}: f32 {} vs f64 {}",
                r32[0][i],
                r64[0][i]
            );
        }
    }
}
