//! Cranelift JIT backend for oCAS expression evaluation.
//!
//! Compiles [`Instr`](crate::Instr) sequences to native machine code via
//! Cranelift, enabling ≥10x faster evaluation compared to the interpreter.
//!
//! Enabled with the `jit` feature flag.

use cranelift_codegen::Context;
use cranelift_codegen::ir::types::{F32, F64, I64};
use cranelift_codegen::ir::{AbiParam, Function, InstBuilder, Signature, UserFuncName, Value};
use cranelift_frontend::{FunctionBuilder, FunctionBuilderContext};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{FuncId, Linkage, Module};
use ocas_core::FastHashMap as HashMap;

use crate::error::{EvaluationError, Result};
use crate::instruction::Instr;

/// Floating-point width for JIT code generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FloatWidth {
    /// 32-bit single precision (`f32` lane, libm `*f` variants).
    F32,
    /// 64-bit double precision (`f64` lane, standard libm).
    F64,
}

impl FloatWidth {
    fn ty(self) -> cranelift_codegen::ir::Type {
        match self {
            FloatWidth::F32 => F32,
            FloatWidth::F64 => F64,
        }
    }

    fn bytes(self) -> i64 {
        match self {
            FloatWidth::F32 => 4,
            FloatWidth::F64 => 8,
        }
    }

    /// libm symbol name for a builtin at this width.
    fn libm_name(self, base: &str) -> String {
        match self {
            FloatWidth::F32 => format!("{base}f"),
            FloatWidth::F64 => base.to_string(),
        }
    }
}

/// A JIT-compiled function ready for execution.
pub struct JitCompiledFunction {
    _module: JITModule,
    entry: extern "C" fn(*const f64, *mut f64),
    param_count: usize,
    result_count: usize,
}

impl JitCompiledFunction {
    /// Call the compiled function with the given f64 parameters.
    pub fn call(&self, params: &[f64]) -> Vec<f64> {
        let mut results = vec![0.0f64; self.result_count];
        (self.entry)(params.as_ptr(), results.as_mut_ptr());
        results
    }

    /// Call the compiled function, writing results into a
    /// caller-provided buffer (avoids per-call heap allocation in hot
    /// loops). `results.len()` must be at least
    /// [`result_count`](JitCompiledFunction::result_count).
    pub fn call_into(&self, params: &[f64], results: &mut [f64]) {
        debug_assert!(results.len() >= self.result_count);
        (self.entry)(params.as_ptr(), results.as_mut_ptr());
    }

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.param_count
    }

    /// Return the number of results produced.
    pub fn result_count(&self) -> usize {
        self.result_count
    }
}

/// A JIT-compiled single-precision function ready for execution.
pub struct JitCompiledF32 {
    _module: JITModule,
    entry: extern "C" fn(*const f32, *mut f32),
    param_count: usize,
    result_count: usize,
}

impl JitCompiledF32 {
    /// Call the compiled function with the given f32 parameters.
    pub fn call(&self, params: &[f32]) -> Vec<f32> {
        let mut results = vec![0.0f32; self.result_count];
        (self.entry)(params.as_ptr(), results.as_mut_ptr());
        results
    }

    /// Call the compiled function, writing results into a
    /// caller-provided buffer (avoids per-call heap allocation in hot
    /// loops). `results.len()` must be at least
    /// [`result_count`](JitCompiledF32::result_count).
    pub fn call_into(&self, params: &[f32], results: &mut [f32]) {
        debug_assert!(results.len() >= self.result_count);
        (self.entry)(params.as_ptr(), results.as_mut_ptr());
    }

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.param_count
    }

    /// Return the number of results produced.
    pub fn result_count(&self) -> usize {
        self.result_count
    }
}

/// JIT compilation engine for oCAS expressions.
pub struct JitEngine;

impl JitEngine {
    /// Compile an instruction sequence into a callable JIT function.
    ///
    /// The stack layout is `[params | constants | temps]`: stack indices
    /// `0..param_count` refer to parameters, indices
    /// `param_count..param_count + constants.len()` to the given constants
    /// (embedded as immediates), and higher indices to temporaries produced
    /// by the instruction sequence. `result_indices` lists the stack slots
    /// whose values are written to the output buffer, in order.
    pub fn compile(
        instructions: &[Instr],
        param_count: usize,
        constants: &[f64],
        result_indices: &[usize],
    ) -> Result<JitCompiledFunction> {
        let (module, entry_ptr, _, _) = compile_module(
            instructions,
            param_count,
            constants,
            result_indices,
            FloatWidth::F64,
        )?;
        // Safety: Cranelift guarantees the pointer is a valid function
        // with the signature we declared.
        let entry: extern "C" fn(*const f64, *mut f64) = unsafe { std::mem::transmute(entry_ptr) };
        Ok(JitCompiledFunction {
            _module: module,
            entry,
            param_count,
            result_count: result_indices.len(),
        })
    }

    /// Compile an instruction sequence into a single-precision JIT
    /// function. Same stack layout as [`compile`](JitEngine::compile);
    /// libm calls use the `*f` (f32) symbol variants.
    pub fn compile_f32(
        instructions: &[Instr],
        param_count: usize,
        constants: &[f32],
        result_indices: &[usize],
    ) -> Result<JitCompiledF32> {
        // f32 → f64 is lossless; the backend re-narrows when emitting.
        let constants64: Vec<f64> = constants.iter().map(|&c| c as f64).collect();
        let (module, entry_ptr, _, _) = compile_module(
            instructions,
            param_count,
            &constants64,
            result_indices,
            FloatWidth::F32,
        )?;
        // Safety: Cranelift guarantees the pointer is a valid function
        // with the signature we declared.
        let entry: extern "C" fn(*const f32, *mut f32) = unsafe { std::mem::transmute(entry_ptr) };
        Ok(JitCompiledF32 {
            _module: module,
            entry,
            param_count,
            result_count: result_indices.len(),
        })
    }
}

/// Shared codegen: build and finalize the JIT module, returning the
/// module (kept alive by the caller) and the finalized entry pointer.
fn compile_module(
    instructions: &[Instr],
    param_count: usize,
    constants: &[f64],
    result_indices: &[usize],
    width: FloatWidth,
) -> Result<(JITModule, *const u8, usize, usize)> {
    let builder = JITBuilder::with_flags(&[], cranelift_module::default_libcall_names())?;
    let mut module = JITModule::new(builder);

    let call_conv = module.isa().default_call_conv();

    let mut sig = Signature::new(call_conv);
    sig.params.push(AbiParam::new(I64));
    sig.params.push(AbiParam::new(I64));

    let mut fn_builder_ctx = FunctionBuilderContext::new();

    let mut ctx = Context::new();

    let func_id = build_ir(
        &mut module,
        &mut ctx,
        instructions,
        param_count,
        constants,
        result_indices,
        width,
        &sig,
        &mut fn_builder_ctx,
    )?;

    module.define_function(func_id, &mut ctx)?;
    module.finalize_definitions()?;

    let entry_ptr = module.get_finalized_function(func_id);
    Ok((module, entry_ptr, param_count, result_indices.len()))
}

#[allow(clippy::too_many_arguments)]
fn build_ir(
    module: &mut JITModule,
    ctx: &mut Context,
    instructions: &[Instr],
    param_count: usize,
    constants: &[f64],
    result_indices: &[usize],
    width: FloatWidth,
    sig: &Signature,
    fn_builder_ctx: &mut FunctionBuilderContext,
) -> Result<FuncId> {
    let func_id = module.declare_function("eval", Linkage::Export, sig)?;

    ctx.func = Function::with_name_signature(UserFuncName::user(0, 0), sig.clone());

    let mut builder = FunctionBuilder::new(&mut ctx.func, fn_builder_ctx);
    let block = builder.create_block();
    builder.append_block_params_for_function_params(block);
    builder.switch_to_block(block);

    let params_ptr = builder.block_params(block)[0];
    let results_ptr = builder.block_params(block)[1];

    let fty = width.ty();
    let stride = width.bytes();

    let mut slots: HashMap<usize, Value> = HashMap::default();

    // Emit a float constant of the active width.
    let emit_const = |builder: &mut FunctionBuilder, c: f64| match width {
        FloatWidth::F32 => builder.ins().f32const(c as f32),
        FloatWidth::F64 => builder.ins().f64const(c),
    };

    // Load parameters
    for i in 0..param_count {
        let offset = i as i64 * stride;
        let addr = builder.ins().iadd_imm(params_ptr, offset);
        let val = builder
            .ins()
            .load(fty, cranelift_codegen::ir::MemFlags::new(), addr, 0);
        slots.insert(i, val);
    }

    // Embed constants as immediates
    for (i, c) in constants.iter().enumerate() {
        let val = emit_const(&mut builder, *c);
        slots.insert(param_count + i, val);
    }

    // Resolve a stack slot, erroring on dangling references instead of
    // panicking.
    macro_rules! slot {
        ($idx:expr) => {
            slots
                .get($idx)
                .copied()
                .ok_or_else(|| EvaluationError::JitCompilationError {
                    message: format!("dangling stack slot reference {}", $idx),
                })?
        };
    }

    // Execute instructions
    for instr in instructions {
        match instr {
            Instr::Copy { dst, src } => {
                let val = slot!(src);
                slots.insert(*dst, val);
            }
            Instr::Add { dst, srcs } => {
                let mut acc = slot!(&srcs[0]);
                for src in &srcs[1..] {
                    acc = builder.ins().fadd(acc, slot!(src));
                }
                slots.insert(*dst, acc);
            }
            Instr::Mul { dst, srcs } => {
                let mut acc = slot!(&srcs[0]);
                for src in &srcs[1..] {
                    acc = builder.ins().fmul(acc, slot!(src));
                }
                slots.insert(*dst, acc);
            }
            Instr::Pow { dst, base, exp } => {
                let base_val = slot!(base);
                // Small integer exponents lower to multiplication chains
                // (exponentiation by squaring) — much faster than libm pow.
                let result = if exp.unsigned_abs() <= 16 {
                    powi_chain(&mut builder, base_val, *exp, width)
                } else {
                    let exp_val = emit_const(&mut builder, *exp as f64);
                    call_libm(
                        &mut builder,
                        module,
                        &width.libm_name("pow"),
                        &[base_val, exp_val],
                    )?
                };
                slots.insert(*dst, result);
            }
            Instr::Powf { dst, base, exp } => {
                let base_val = slot!(base);
                let exp_val = slot!(exp);
                let result = call_libm(
                    &mut builder,
                    module,
                    &width.libm_name("pow"),
                    &[base_val, exp_val],
                )?;
                slots.insert(*dst, result);
            }
            Instr::BuiltinOp { dst, op, src } => {
                use crate::instruction::BuiltinOp;
                let src_val = slot!(src);
                let fn_name = match op {
                    BuiltinOp::Sin => width.libm_name("sin"),
                    BuiltinOp::Cos => width.libm_name("cos"),
                    BuiltinOp::Tan => width.libm_name("tan"),
                    BuiltinOp::Exp => width.libm_name("exp"),
                    BuiltinOp::Log => width.libm_name("log"),
                    BuiltinOp::Sqrt => width.libm_name("sqrt"),
                    BuiltinOp::Abs => width.libm_name("fabs"),
                    _ => {
                        return Err(EvaluationError::UnsupportedOperation {
                            message: format!("JIT does not support builtin op {op:?}"),
                        });
                    }
                };
                let result = call_libm(&mut builder, module, &fn_name, &[src_val])?;
                slots.insert(*dst, result);
            }
            Instr::ExternalFun { .. } => {
                return Err(EvaluationError::JitCompilationError {
                    message: "external functions not supported in JIT".into(),
                });
            }
        }
    }

    // Store results in the deterministic order given by result_indices.
    for (i, slot_idx) in result_indices.iter().enumerate() {
        let val = slot!(slot_idx);
        let offset = i as i64 * stride;
        let addr = builder.ins().iadd_imm(results_ptr, offset);
        builder
            .ins()
            .store(cranelift_codegen::ir::MemFlags::new(), val, addr, 0);
    }

    builder.ins().return_(&[]);
    builder.seal_block(block);
    builder.finalize();

    Ok(func_id)
}

/// Emit `base^exp` for an integer exponent via exponentiation by
/// squaring (negative exponents produce `1 / base^|exp|`).
fn powi_chain(builder: &mut FunctionBuilder, base: Value, exp: i64, width: FloatWidth) -> Value {
    let one = match width {
        FloatWidth::F32 => builder.ins().f32const(1.0),
        FloatWidth::F64 => builder.ins().f64const(1.0),
    };
    let mut result = one;
    let mut b = base;
    let mut e = exp.unsigned_abs();
    while e > 0 {
        if e & 1 == 1 {
            result = builder.ins().fmul(result, b);
        }
        e >>= 1;
        if e > 0 {
            b = builder.ins().fmul(b, b);
        }
    }
    if exp < 0 {
        let one = match width {
            FloatWidth::F32 => builder.ins().f32const(1.0),
            FloatWidth::F64 => builder.ins().f64const(1.0),
        };
        result = builder.ins().fdiv(one, result);
    }
    result
}

/// Call a libm function (sin, cos, pow, etc.) and return the result Value.
fn call_libm(
    builder: &mut FunctionBuilder,
    module: &mut JITModule,
    name: &str,
    args: &[Value],
) -> Result<Value> {
    // Argument types match the function's float width; Cranelift
    // validates call signatures against declared argument Values.
    let arg_ty = builder.func.dfg.value_type(args[0]);
    let param_types: Vec<AbiParam> = args.iter().map(|_| AbiParam::new(arg_ty)).collect();
    let callee_sig = Signature {
        params: param_types,
        returns: vec![AbiParam::new(arg_ty)],
        call_conv: module.isa().default_call_conv(),
    };

    let callee_id = module.declare_function(name, Linkage::Import, &callee_sig)?;
    let local_callee = module.declare_func_in_func(callee_id, builder.func);
    let call = builder.ins().call(local_callee, args);
    Ok(builder.inst_results(call)[0])
}

impl From<cranelift_module::ModuleError> for EvaluationError {
    fn from(err: cranelift_module::ModuleError) -> Self {
        EvaluationError::JitCompilationError {
            message: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jit_constant_return() {
        // Simplest possible JIT: empty function, tests basic pipeline
        let func = JitEngine::compile(&[], 0, &[], &[]);
        assert!(func.is_ok());
    }

    #[test]
    fn jit_simple_add() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile(&instructions, 2, &[], &[2]).unwrap();
        let result = func.call(&[2.0, 3.0]);
        assert!((result[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn jit_simple_mul() {
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile(&instructions, 2, &[], &[2]).unwrap();
        let result = func.call(&[3.0, 4.0]);
        assert!((result[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn jit_builtin_sin() {
        let instructions = vec![Instr::BuiltinOp {
            dst: 1,
            op: crate::instruction::BuiltinOp::Sin,
            src: 0,
        }];
        let func = JitEngine::compile(&instructions, 1, &[], &[1]).unwrap();
        let result = func.call(&[std::f64::consts::FRAC_PI_2]);
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn jit_pow_integer() {
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: 3,
        }];
        let func = JitEngine::compile(&instructions, 1, &[], &[1]).unwrap();
        let result = func.call(&[2.0]);
        assert!((result[0] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn jit_pow_negative_exponent() {
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: -2,
        }];
        let func = JitEngine::compile(&instructions, 1, &[], &[1]).unwrap();
        let result = func.call(&[4.0]);
        assert!((result[0] - 0.0625).abs() < 1e-12);
    }

    #[test]
    fn jit_call_into_reuses_buffer() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile(&instructions, 2, &[], &[2]).unwrap();
        let mut out = [0.0f64; 1];
        for i in 0..10 {
            func.call_into(&[i as f64, 1.0], &mut out);
            assert!((out[0] - (i as f64 + 1.0)).abs() < 1e-10);
        }
    }

    #[test]
    fn jit_multi_output() {
        // params: [x, y]; outputs: [x + y, x * y]
        let instructions = vec![
            Instr::Add {
                dst: 2,
                srcs: vec![0, 1],
            },
            Instr::Mul {
                dst: 3,
                srcs: vec![0, 1],
            },
        ];
        let func = JitEngine::compile(&instructions, 2, &[], &[2, 3]).unwrap();
        let result = func.call(&[2.0, 3.0]);
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn jit_multi_output_shared_subexpression() {
        // params: [x]; outputs: [sin(x) + 1, sin(x) * 2]
        // stack: 0=x, 1..=2=constants (1.0, 2.0), 3=sin(x), 4=add, 5=mul
        let instructions = vec![
            Instr::BuiltinOp {
                dst: 3,
                op: crate::instruction::BuiltinOp::Sin,
                src: 0,
            },
            Instr::Add {
                dst: 4,
                srcs: vec![3, 1],
            },
            Instr::Mul {
                dst: 5,
                srcs: vec![3, 2],
            },
        ];
        let func = JitEngine::compile(&instructions, 1, &[1.0, 2.0], &[4, 5]).unwrap();
        let result = func.call(&[std::f64::consts::FRAC_PI_2]);
        assert!((result[0] - 2.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn jit_dangling_slot_errors() {
        // Referencing a slot that was never defined must be an error, not a panic.
        let instructions = vec![Instr::Copy { dst: 5, src: 9 }];
        let result = JitEngine::compile(&instructions, 1, &[], &[5]);
        assert!(result.is_err());
    }

    #[test]
    fn jit_from_evaluator_multi_output() {
        // End-to-end: parse-level atoms → multi-output evaluator → JIT
        use ocas_atom::AtomArena;
        use ocas_core::arena::Arena;

        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let sum = ctx.add(&[ctx.var("x"), ctx.var("y")]);
        let prod = ctx.mul(&[ctx.var("x"), ctx.var("y")]);
        let eval: crate::ExpressionEvaluator<f64> =
            crate::ExpressionEvaluator::compile_multi(&[sum, prod]).unwrap();
        let func = eval.compile_jit().unwrap();
        assert_eq!(func.param_count(), 2);
        assert_eq!(func.result_count(), 2);
        let result = func.call(&[2.0, 3.0]);
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn jit_from_evaluator_with_constants() {
        // x^4 + 3x^3 + 2x^2 + x + 5 via the full compile pipeline
        use ocas_atom::AtomArena;
        use ocas_core::arena::Arena;

        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[
            ctx.pow(x, ctx.num(4)),
            ctx.mul(&[ctx.num(3), ctx.pow(x, ctx.num(3))]),
            ctx.mul(&[ctx.num(2), ctx.pow(x, ctx.num(2))]),
            x,
            ctx.num(5),
        ]);
        let eval: crate::ExpressionEvaluator<f64> =
            crate::ExpressionEvaluator::compile(expr).unwrap();
        let func = eval.compile_jit().unwrap();
        let expected = eval.evaluate(&[1.5]).unwrap();
        let result = func.call(&[1.5]);
        assert!((result[0] - expected[0]).abs() < 1e-9);
    }

    #[test]
    fn jit_f32_simple() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile_f32(&instructions, 2, &[], &[2]).unwrap();
        let result = func.call(&[2.0f32, 3.0f32]);
        assert!((result[0] - 5.0f32).abs() < 1e-6);
    }

    #[test]
    fn jit_f32_multi_output_with_constants() {
        // params: [x]; outputs: [sin(x) + 1, x^3]
        // stack: 0=x, 1=const(1.0), 2=sin, 3=add, 4=pow
        let instructions = vec![
            Instr::BuiltinOp {
                dst: 2,
                op: crate::instruction::BuiltinOp::Sin,
                src: 0,
            },
            Instr::Add {
                dst: 3,
                srcs: vec![2, 1],
            },
            Instr::Pow {
                dst: 4,
                base: 0,
                exp: 3,
            },
        ];
        let func = JitEngine::compile_f32(&instructions, 1, &[1.0], &[3, 4]).unwrap();
        assert_eq!(func.result_count(), 2);
        let result = func.call(&[std::f32::consts::FRAC_PI_2]);
        assert!((result[0] - 2.0f32).abs() < 1e-5);
        let expect = std::f32::consts::FRAC_PI_2.powi(3);
        assert!((result[1] - expect).abs() < 1e-4);
    }

    #[test]
    fn jit_f32_matches_f64_within_precision() {
        // Compare f32 and f64 JIT results on the same expression
        use ocas_atom::AtomArena;
        use ocas_core::arena::Arena;

        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[
            ctx.pow(x, ctx.num(4)),
            ctx.mul(&[ctx.num(3), ctx.pow(x, ctx.num(3))]),
            ctx.mul(&[ctx.num(2), ctx.pow(x, ctx.num(2))]),
            x,
            ctx.num(5),
        ]);
        let eval: crate::ExpressionEvaluator<f64> =
            crate::ExpressionEvaluator::compile(expr).unwrap();
        let f64_func = eval.compile_jit().unwrap();
        let f32_func = eval.compile_jit_f32().unwrap();
        let r64 = f64_func.call(&[1.25]);
        let r32 = f32_func.call(&[1.25f32]);
        // f32 has ~7 decimal digits; require agreement within 1e-4 relative
        let rel = ((r32[0] as f64 - r64[0]) / r64[0]).abs();
        assert!(rel < 1e-4, "f32 {} vs f64 {}", r32[0], r64[0]);
    }

    #[test]
    fn jit_f32_call_into() {
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile_f32(&instructions, 2, &[], &[2]).unwrap();
        let mut out = [0.0f32; 1];
        func.call_into(&[3.0f32, 4.0f32], &mut out);
        assert!((out[0] - 12.0f32).abs() < 1e-6);
    }
}
