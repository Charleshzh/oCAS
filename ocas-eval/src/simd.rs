//! SIMD vectorized expression evaluation.
//!
//! The [`VectorEvaluator`] evaluates expressions on batches of input
//! vectors using `pulp` for portable SIMD computation. The SIMD width
//! is detected at runtime (SSE2 → 2, AVX2 → 4, AVX-512 → 8 lanes).
//!
//! Enabled with the `simd` feature flag.

use pulp::Arch;

use crate::error::{EvaluationError, Result};
use crate::instruction::Instr;

/// Batch evaluator using SIMD for vectorized expression computation.
///
/// Processes input vectors in chunks whose width is detected at runtime,
/// falling back to scalar evaluation for any remainder.
pub struct VectorEvaluator {
    instructions: Vec<Instr>,
    param_count: usize,
    #[allow(dead_code)]
    const_count: usize,
    stack_size: usize,
    result_indices: Vec<usize>,
    constants: Vec<f64>,
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
            instructions,
            param_count,
            const_count,
            stack_size,
            result_indices,
            constants,
        }
    }

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.param_count
    }

    /// Evaluate the expression on a batch of inputs.
    ///
    /// Each element of `params` is a vector of input values for one
    /// parameter. All input vectors must have the same length.
    ///
    /// Returns a vector of result values, one per input row.
    pub fn evaluate(&self, params: &[Vec<f64>]) -> Result<Vec<f64>> {
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
        let mut results = vec![0.0f64; batch_size];

        // Determine the SIMD lane count. On AVX-512 this is 8, AVX2 → 4,
        // SSE2 → 2, scalar fallback → 1.
        let lanes = if cfg!(target_arch = "x86_64") {
            if std::is_x86_feature_detected!("avx512f") {
                8
            } else if std::is_x86_feature_detected!("avx2") {
                4
            } else {
                2 // SSE2 is baseline on x86_64
            }
        } else {
            1
        };

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
            let scalar_params: Vec<f64> = params.iter().map(|p| p[i]).collect();
            let scalar_results = self.evaluate_scalar(&scalar_params)?;
            for (ri, &res_idx) in self.result_indices.iter().enumerate() {
                if ri == 0 {
                    results[i] = scalar_results[res_idx];
                }
            }
        }

        Ok(results)
    }

    /// Scalar evaluation for a single set of parameters.
    fn evaluate_scalar(&self, params: &[f64]) -> Result<Vec<f64>> {
        let mut stack: Vec<f64> = vec![0.0f64; self.stack_size];

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
                        sum += stack[*idx];
                    }
                    stack[*dst] = sum;
                }
                Instr::Mul { dst, srcs } => {
                    let mut prod = stack[srcs[0]];
                    for idx in &srcs[1..] {
                        prod *= stack[*idx];
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
                    use crate::instruction::BuiltinOp;
                    let x = stack[*src];
                    let result = match op {
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
                    };
                    stack[*dst] = result;
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
/// Each stack slot holds an `[f64; 8]` array where only `lanes` elements
/// are active (2, 4, or 8 depending on CPU features). Arithmetic is
/// performed lane-by-lane; transcendental functions fall back to scalar
/// per lane.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
#[inline]
fn eval_simd_chunks(
    instructions: &[Instr],
    constants: &[f64],
    param_count: usize,
    stack_size: usize,
    result_indices: &[usize],
    params: &[Vec<f64>],
    out_results: &mut [f64],
    chunks: usize,
    lanes: usize,
) {
    // Pre-allocate stack buffer once, reuse across all chunks.
    let mut stack: Vec<[f64; 8]> = vec![[0.0; 8]; stack_size];

    for chunk in 0..chunks {
        let base = chunk * lanes;

        // Clear stack for this chunk (reuse allocation).
        stack.fill([0.0; 8]);

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
                            sum[l] += other[l];
                        }
                    }
                    stack[*dst] = sum;
                }
                Instr::Mul { dst, srcs } => {
                    let mut prod = stack[srcs[0]];
                    for &idx in &srcs[1..] {
                        let other = stack[idx];
                        for l in 0..lanes {
                            prod[l] *= other[l];
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
                    let mut result = [1.0f64; 8];
                    let mut e = exp.abs();
                    let mut b = base_val;
                    while e > 0 {
                        if e & 1 == 1 {
                            for l in 0..lanes {
                                result[l] *= b[l];
                            }
                        }
                        let mut new_b = [0.0f64; 8];
                        for l in 0..lanes {
                            new_b[l] = b[l] * b[l];
                        }
                        b = new_b;
                        e >>= 1;
                    }
                    if *exp < 0 {
                        for l in 0..lanes {
                            result[l] = 1.0 / result[l];
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
                    let mut result = [0.0f64; 8];
                    for l in 0..lanes {
                        result[l] = base_val[l].powf(exp_val[l]);
                    }
                    stack[*dst] = result;
                }
                Instr::BuiltinOp { dst, op, src } => {
                    use crate::instruction::BuiltinOp;
                    let arr = stack[*src];
                    let mut result = [0.0f64; 8];
                    for l in 0..lanes {
                        result[l] = match op {
                            BuiltinOp::Sin => arr[l].sin(),
                            BuiltinOp::Cos => arr[l].cos(),
                            BuiltinOp::Tan => arr[l].tan(),
                            BuiltinOp::Sec => 1.0 / arr[l].cos(),
                            BuiltinOp::Csc => 1.0 / arr[l].sin(),
                            BuiltinOp::Cot => 1.0 / arr[l].tan(),
                            BuiltinOp::Exp => arr[l].exp(),
                            BuiltinOp::Log => arr[l].ln(),
                            BuiltinOp::Sqrt => arr[l].sqrt(),
                            BuiltinOp::Abs => arr[l].abs(),
                        };
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

        // Write results for the first result index.
        for (ri, &res_idx) in result_indices.iter().enumerate() {
            if ri == 0 {
                out_results[base..(base + lanes)].copy_from_slice(&stack[res_idx][..lanes]);
            }
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

        for (i, &r) in results.iter().enumerate() {
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

        for (i, &r) in results.iter().enumerate() {
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
        assert_eq!(results.len(), 3);
        assert!((results[0] - 11.0).abs() < 1e-10);
        assert!((results[1] - 12.0).abs() < 1e-10);
        assert!((results[2] - 13.0).abs() < 1e-10);
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

        for (i, &r) in results.iter().enumerate() {
            let expected = (i as f64).powi(3);
            assert!((r - expected).abs() < 1e-10, "mismatch at {i}");
        }
    }
}
