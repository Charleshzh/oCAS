//! Stack-based expression evaluator.
//!
//! The [`ExpressionEvaluator`] executes a pre-compiled instruction sequence
//! on a flat stack, producing numeric results from input parameters.

use crate::domain::{EvaluationDomain, PowfExtension};
use crate::error::{EvaluationError, Result};
use crate::function_map::FunctionMap;
use crate::instruction::Instr;

/// A compiled expression ready for numeric evaluation.
///
/// The evaluator holds a sequence of [`Instr`]s and a pre-allocated stack.
/// Call [`evaluate`](ExpressionEvaluator::evaluate) with parameter values
/// to compute the result.
pub struct ExpressionEvaluator<T: EvaluationDomain> {
    /// The instruction sequence to execute.
    instructions: Vec<Instr>,
    /// Number of parameter slots at the start of the stack.
    param_count: usize,
    /// Number of constant slots (after params).
    #[allow(dead_code)]
    const_count: usize,
    /// Total stack size (params + constants + temporaries + outputs).
    stack_size: usize,
    /// Indices of result slots in the stack.
    result_indices: Vec<usize>,
    /// Pre-computed constant values.
    constants: Vec<T>,
    /// Optional user-defined function registry.
    function_map: Option<FunctionMap<T>>,
}

impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    /// Create an evaluator from compiled instruction data.
    ///
    /// This is a low-level constructor used by the compiler. Most users
    /// should use [`ExpressionEvaluator::compile`] instead.
    #[allow(dead_code)]
    pub(crate) fn new(
        instructions: Vec<Instr>,
        param_count: usize,
        const_count: usize,
        stack_size: usize,
        result_indices: Vec<usize>,
        constants: Vec<T>,
    ) -> Self {
        Self {
            instructions,
            param_count,
            const_count,
            stack_size,
            result_indices,
            constants,
            function_map: None,
        }
    }

    /// Create an evaluator with a function map for user-defined functions.
    #[allow(dead_code)]
    pub(crate) fn new_with_functions(
        instructions: Vec<Instr>,
        param_count: usize,
        const_count: usize,
        stack_size: usize,
        result_indices: Vec<usize>,
        constants: Vec<T>,
        function_map: FunctionMap<T>,
    ) -> Self {
        Self {
            instructions,
            param_count,
            const_count,
            stack_size,
            result_indices,
            constants,
            function_map: Some(function_map),
        }
    }

    /// Return the number of parameters expected by this evaluator.
    pub fn param_count(&self) -> usize {
        self.param_count
    }

    /// Return the number of results produced by this evaluator.
    pub fn result_count(&self) -> usize {
        self.result_indices.len()
    }

    /// Return the stack size required to evaluate this expression.
    ///
    /// Callers reusing a stack buffer across evaluations should
    /// pre-allocate a `Vec` with this capacity.
    pub fn stack_size(&self) -> usize {
        self.stack_size
    }

    /// Evaluate the expression with the given parameter values.
    ///
    /// Returns a vector of result values. The number of results equals
    /// `result_indices.len()`.
    ///
    /// # Errors
    ///
    /// Returns [`EvaluationError`] if the number of parameters does not
    /// match, or if an arithmetic error occurs (e.g. division by zero).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let result = evaluator.evaluate(&[2.0_f64])?;
    /// assert_eq!(result.len(), 1);
    /// ```
    pub fn evaluate(&self, params: &[T]) -> Result<Vec<T>> {
        let mut stack: Vec<T> = Vec::with_capacity(self.stack_size);
        let mut results: Vec<T> = Vec::with_capacity(self.result_indices.len());
        self.evaluate_with_stack(params, &mut stack, &mut results)?;
        Ok(results)
    }

    /// Evaluate the expression, reusing caller-provided buffers.
    ///
    /// `stack` must have capacity for at least
    /// [`stack_size`](ExpressionEvaluator::stack_size) elements (it is
    /// resized as needed); `results` is filled with the result values.
    /// Reusing buffers across calls avoids per-evaluation heap
    /// allocation, which matters for streaming and batch workloads.
    ///
    /// # Errors
    ///
    /// Returns [`EvaluationError`] if the number of parameters does not
    /// match, or if an arithmetic error occurs.
    pub fn evaluate_with_stack(
        &self,
        params: &[T],
        stack: &mut Vec<T>,
        results: &mut Vec<T>,
    ) -> Result<()> {
        if params.len() != self.param_count {
            return Err(EvaluationError::WrongArity {
                name: "<expr>".into(),
                expected: self.param_count,
                got: params.len(),
            });
        }

        stack.clear();
        stack.resize(self.stack_size, T::zero());

        // Fill parameters
        for (i, p) in params.iter().enumerate() {
            stack[i] = p.clone();
        }

        // Fill constants
        for (i, c) in self.constants.iter().enumerate() {
            stack[self.param_count + i] = c.clone();
        }

        // Execute instructions
        for instr in &self.instructions {
            match instr {
                Instr::Add { dst, srcs } => {
                    let mut sum = stack[srcs[0]].clone();
                    for idx in &srcs[1..] {
                        sum = sum.add_ref(&stack[*idx]);
                    }
                    stack[*dst] = sum;
                }
                Instr::Mul { dst, srcs } => {
                    let mut prod = stack[srcs[0]].clone();
                    for idx in &srcs[1..] {
                        prod = prod.mul_ref(&stack[*idx]);
                    }
                    stack[*dst] = prod;
                }
                Instr::Pow { dst, base, exp } => {
                    stack[*dst] = stack[*base].powi_ref(*exp);
                }
                Instr::Powf { dst, base, exp } => {
                    let result = stack[*base].powf_ref(&stack[*exp])?;
                    stack[*dst] = result;
                }
                Instr::BuiltinOp { dst, op, src } => {
                    let name = match op {
                        crate::instruction::BuiltinOp::Sin => "sin",
                        crate::instruction::BuiltinOp::Cos => "cos",
                        crate::instruction::BuiltinOp::Tan => "tan",
                        crate::instruction::BuiltinOp::Sec => "sec",
                        crate::instruction::BuiltinOp::Csc => "csc",
                        crate::instruction::BuiltinOp::Cot => "cot",
                        crate::instruction::BuiltinOp::Exp => "exp",
                        crate::instruction::BuiltinOp::Log => "log",
                        crate::instruction::BuiltinOp::Sqrt => "sqrt",
                        crate::instruction::BuiltinOp::Abs => "abs",
                    };
                    let result = T::resolve_builtin(name, &stack[*src])?;
                    stack[*dst] = result;
                }
                Instr::ExternalFun { dst, fn_idx, srcs } => {
                    let args: Vec<T> = srcs.iter().map(|&i| stack[i].clone()).collect();
                    let result = self
                        .function_map
                        .as_ref()
                        .and_then(|fm| fm.call_by_index(*fn_idx, &args))
                        .ok_or_else(|| EvaluationError::FunctionNotFound {
                            name: format!("external function at index {fn_idx}"),
                        })?;
                    stack[*dst] = result;
                }
                Instr::Copy { dst, src } => {
                    stack[*dst] = stack[*src].clone();
                }
            }
        }

        // Collect results
        results.clear();
        results.extend(self.result_indices.iter().map(|&i| stack[i].clone()));

        Ok(())
    }
}

/// Cranelift JIT compilation support.
#[cfg(feature = "jit")]
impl ExpressionEvaluator<f64> {
    /// Compile this evaluator's instruction sequence to native machine
    /// code via the Cranelift JIT backend.
    ///
    /// Constants are embedded as immediates and all result slots are
    /// written in order, so multi-output evaluators are fully supported.
    ///
    /// # Errors
    ///
    /// Returns [`EvaluationError::UnsupportedOperation`] if the expression
    /// contains instructions the JIT cannot lower (e.g. external
    /// functions or `sec`/`csc`/`cot`), or
    /// [`EvaluationError::JitCompilationError`] on backend failure.
    pub fn compile_jit(&self) -> Result<crate::jit::JitCompiledFunction> {
        let constants: Vec<f64> = self.constants.clone();
        crate::jit::JitEngine::compile(
            &self.instructions,
            self.param_count,
            &constants,
            &self.result_indices,
        )
    }

    /// Compile this evaluator's instruction sequence to single-precision
    /// native code. Constants are narrowed from f64 to f32; results have
    /// f32 precision.
    ///
    /// # Errors
    ///
    /// Same conditions as [`compile_jit`](ExpressionEvaluator::compile_jit).
    pub fn compile_jit_f32(&self) -> Result<crate::jit::JitCompiledF32> {
        let constants: Vec<f32> = self.constants.iter().map(|&c| c as f32).collect();
        crate::jit::JitEngine::compile_f32(
            &self.instructions,
            self.param_count,
            &constants,
            &self.result_indices,
        )
    }
}

/// SIMD batch evaluation support.
#[cfg(feature = "simd")]
impl ExpressionEvaluator<f64> {
    /// Compile this evaluator into a [`VectorEvaluator`] for batch SIMD evaluation.
    ///
    /// The resulting evaluator processes multiple input values simultaneously
    /// using the best available SIMD width (SSE2/AVX2/AVX-512).
    ///
    /// # Errors
    ///
    /// Returns [`EvaluationError::UnsupportedOperation`] if the expression
    /// contains external functions, which are not supported in SIMD mode.
    pub fn compile_vector_evaluator(&self) -> Result<crate::simd::VectorEvaluator> {
        // Check for unsupported instructions
        for instr in &self.instructions {
            if let Instr::ExternalFun { .. } = instr {
                return Err(EvaluationError::UnsupportedOperation {
                    message: "external functions not supported in SIMD mode".into(),
                });
            }
        }

        Ok(crate::simd::VectorEvaluator::new(
            self.instructions.clone(),
            self.param_count,
            self.const_count,
            self.stack_size,
            self.result_indices.clone(),
            self.constants.clone(),
        ))
    }

    /// Compile this evaluator into a single-precision
    /// [`VectorEvaluatorF32`](crate::simd::VectorEvaluatorF32) for batch
    /// SIMD evaluation. Constants are narrowed from f64 to f32; on the
    /// same hardware this doubles the SIMD lane count.
    ///
    /// # Errors
    ///
    /// Returns [`EvaluationError::UnsupportedOperation`] if the expression
    /// contains external functions, which are not supported in SIMD mode.
    pub fn compile_vector_evaluator_f32(&self) -> Result<crate::simd::VectorEvaluatorF32> {
        // Check for unsupported instructions
        for instr in &self.instructions {
            if let Instr::ExternalFun { .. } = instr {
                return Err(EvaluationError::UnsupportedOperation {
                    message: "external functions not supported in SIMD mode".into(),
                });
            }
        }

        let constants: Vec<f32> = self.constants.iter().map(|&c| c as f32).collect();
        Ok(crate::simd::VectorEvaluatorF32::new(
            self.instructions.clone(),
            self.param_count,
            self.const_count,
            self.stack_size,
            self.result_indices.clone(),
            constants,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_evaluator() -> ExpressionEvaluator<f64> {
        // Evaluate: x + 1
        // Stack layout: [param(0)=x] [const(0)=1.0] [temp(0)=result]
        // Instructions: Add(temp(0), param(0), const(0))
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![1.0f64];
        ExpressionEvaluator::new(instructions, 1, 1, 3, vec![2], constants)
    }

    #[test]
    fn simple_add() {
        let eval = make_simple_evaluator();
        assert_eq!(eval.param_count(), 1);
        let result = eval.evaluate(&[2.0]).unwrap();
        assert!((result[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn wrong_param_count() {
        let eval = make_simple_evaluator();
        assert!(eval.evaluate(&[1.0, 2.0]).is_err());
        assert!(eval.evaluate(&[]).is_err());
    }

    #[test]
    fn mul_expression() {
        // Evaluate: x * 2
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let constants = vec![2.0f64];
        let eval = ExpressionEvaluator::new(instructions, 1, 1, 3, vec![2], constants);
        assert!((eval.evaluate(&[3.0]).unwrap()[0] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn pow_expression() {
        // Evaluate: x^3
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: 3,
        }];
        let eval = ExpressionEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);
        assert!((eval.evaluate(&[2.0]).unwrap()[0] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn builtin_sin() {
        // Evaluate: sin(x)
        let instructions = vec![Instr::BuiltinOp {
            dst: 1,
            op: crate::instruction::BuiltinOp::Sin,
            src: 0,
        }];
        let eval = ExpressionEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);
        let result = eval.evaluate(&[std::f64::consts::FRAC_PI_2]).unwrap();
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn copy_instruction() {
        // Evaluate: x (identity)
        let instructions = vec![Instr::Copy { dst: 1, src: 0 }];
        let eval = ExpressionEvaluator::new(instructions, 1, 0, 2, vec![1], vec![]);
        assert!((eval.evaluate(&[42.0]).unwrap()[0] - 42.0).abs() < 1e-10);
    }

    #[test]
    fn evaluate_with_stack_reuses_buffers() {
        let eval = make_simple_evaluator();
        let mut stack: Vec<f64> = Vec::with_capacity(eval.stack_size());
        let mut results: Vec<f64> = Vec::with_capacity(eval.result_count());

        for x in 0..100 {
            eval.evaluate_with_stack(&[x as f64], &mut stack, &mut results)
                .unwrap();
            assert!((results[0] - (x as f64 + 1.0)).abs() < 1e-10);
        }
        // Buffers retain capacity across calls (no reallocation needed)
        assert!(stack.capacity() >= eval.stack_size());
    }

    #[test]
    fn evaluate_with_stack_wrong_arity() {
        let eval = make_simple_evaluator();
        let mut stack = Vec::new();
        let mut results = Vec::new();
        assert!(
            eval.evaluate_with_stack(&[1.0, 2.0], &mut stack, &mut results)
                .is_err()
        );
    }

    #[test]
    fn result_and_stack_getters() {
        let eval = make_simple_evaluator();
        assert_eq!(eval.result_count(), 1);
        assert_eq!(eval.stack_size(), 3);
    }
}
