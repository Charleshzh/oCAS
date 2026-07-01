//! SIMD vectorized expression evaluation.
//!
//! The [`VectorEvaluator`] evaluates expressions on batches of input
//! vectors using SIMD primitives for 4-wide parallel computation.
//!
//! Enabled with the `simd` feature flag.

use wide::f64x4;

use crate::domain::EvaluationDomain;
use crate::error::{EvaluationError, Result};
use crate::instruction::Instr;

/// Batch evaluator using SIMD for vectorized expression computation.
///
/// Processes input vectors in chunks of 4 lanes, falling back to scalar
/// evaluation for any remainder.
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
            // No params — evaluate once
            1
        };

        // Validate all param vectors have the same length
        for p in params.iter() {
            if p.len() != batch_size {
                return Err(EvaluationError::TypeMismatch {
                    expected: format!("vector of length {batch_size}"),
                    found: format!("vector of length {}", p.len()),
                });
            }
        }

        let mut results = vec![0.0f64; batch_size];
        let simd_chunks = batch_size / 4;
        let remainder_start = simd_chunks * 4;

        // SIMD path: process 4 elements at a time
        for chunk in 0..simd_chunks {
            let base = chunk * 4;
            let simd_params: Vec<f64x4> = params
                .iter()
                .map(|p| f64x4::new([p[base], p[base + 1], p[base + 2], p[base + 3]]))
                .collect();

            let simd_results = self.evaluate_simd(&simd_params)?;
            for (ri, &res_idx) in self.result_indices.iter().enumerate() {
                let arr = simd_results[res_idx].to_array();
                results[base + ri] = arr[0];
                // Store all 4 results at their positions
                if self.result_indices.len() == 1 {
                    for (j, &v) in arr.iter().enumerate() {
                        if base + j < batch_size {
                            results[base + j] = v;
                        }
                    }
                }
            }
        }

        // Scalar fallback for remainder
        for i in remainder_start..batch_size {
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

    /// SIMD evaluation for exactly 4 input values per parameter.
    fn evaluate_simd(&self, params: &[f64x4]) -> Result<Vec<f64x4>> {
        let mut stack: Vec<f64x4> = vec![f64x4::ZERO; self.stack_size];

        for (i, p) in params.iter().enumerate() {
            stack[i] = *p;
        }

        for (i, c) in self.constants.iter().enumerate() {
            stack[self.param_count + i] = f64x4::splat(*c);
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
                    let base_val = stack[*base];
                    // Manual powi for SIMD
                    let mut result = f64x4::splat(1.0);
                    let mut e = exp.abs();
                    let mut b = base_val;
                    while e > 0 {
                        if e & 1 == 1 {
                            result *= b;
                        }
                        b *= b;
                        e >>= 1;
                    }
                    if *exp < 0 {
                        result = f64x4::splat(1.0) / result;
                    }
                    stack[*dst] = result;
                }
                Instr::Powf { dst, base, exp } => {
                    // Use scalar fallback per lane for powf
                    let base_arr = stack[*base].to_array();
                    let exp_arr = stack[*exp].to_array();
                    let mut result = [0.0f64; 4];
                    for i in 0..4 {
                        result[i] = base_arr[i].powf(exp_arr[i]);
                    }
                    stack[*dst] = f64x4::new([result[0], result[1], result[2], result[3]]);
                }
                Instr::BuiltinFun { dst, name, src } => {
                    let arr = stack[*src].to_array();
                    let fn_name = name.as_str().to_lowercase();
                    let mut result = [0.0f64; 4];
                    for i in 0..4 {
                        result[i] = match fn_name.as_str() {
                            "sin" => arr[i].sin(),
                            "cos" => arr[i].cos(),
                            "tan" => arr[i].tan(),
                            "sec" => 1.0 / arr[i].cos(),
                            "csc" => 1.0 / arr[i].sin(),
                            "cot" => 1.0 / arr[i].tan(),
                            "exp" => arr[i].exp(),
                            "log" => {
                                if arr[i] <= 0.0 {
                                    return Err(EvaluationError::UnsupportedOperation {
                                        message: "log of non-positive number".into(),
                                    });
                                }
                                arr[i].ln()
                            }
                            "sqrt" => {
                                if arr[i] < 0.0 {
                                    return Err(EvaluationError::UnsupportedOperation {
                                        message: "sqrt of negative number".into(),
                                    });
                                }
                                arr[i].sqrt()
                            }
                            "abs" => arr[i].abs(),
                            _ => return Err(EvaluationError::FunctionNotFound { name: fn_name }),
                        };
                    }
                    stack[*dst] = f64x4::new([result[0], result[1], result[2], result[3]]);
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
                Instr::BuiltinFun { dst, name, src } => {
                    let result = f64::resolve_builtin(name.as_str(), &stack[*src])?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_add_constant() {
        // Evaluate: x + 1 for 4 values
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![1.0f64];
        let eval = VectorEvaluator::new(instructions, 1, 1, 3, vec![2], constants);

        let params = vec![vec![1.0, 2.0, 3.0, 4.0]];
        let results = eval.evaluate(&params).unwrap();
        assert_eq!(results, vec![2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn simd_mul_batch() {
        // Evaluate: x * 2 for 8 values (tests chunking)
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![2.0f64];
        let eval = VectorEvaluator::new(instructions, 1, 1, 3, vec![2], constants);

        let params = vec![vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]];
        let results = eval.evaluate(&params).unwrap();
        assert_eq!(results, vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]);
    }

    #[test]
    fn simd_sin() {
        let name = ocas_atom::Symbol::new("sin");
        let instructions = vec![Instr::BuiltinFun {
            dst: 1,
            name,
            src: 0,
        }];
        let eval = VectorEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);

        let pi_half = std::f64::consts::FRAC_PI_2;
        let params = vec![vec![pi_half, 0.0, pi_half, 0.0]];
        let results = eval.evaluate(&params).unwrap();
        assert!((results[0] - 1.0).abs() < 1e-10);
        assert!((results[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn simd_pow_integer() {
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: 3,
        }];
        let eval = VectorEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);

        let params = vec![vec![1.0, 2.0, 3.0, 4.0]];
        let results = eval.evaluate(&params).unwrap();
        assert_eq!(results, vec![1.0, 8.0, 27.0, 64.0]);
    }

    #[test]
    fn simd_two_params() {
        // Evaluate: x + y for pairs
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let eval = VectorEvaluator::new(instructions, 2, 0, 3, vec![2], vec![]);

        let params = vec![vec![1.0, 2.0, 3.0, 4.0], vec![5.0, 6.0, 7.0, 8.0]];
        let results = eval.evaluate(&params).unwrap();
        assert_eq!(results, vec![6.0, 8.0, 10.0, 12.0]);
    }
}
