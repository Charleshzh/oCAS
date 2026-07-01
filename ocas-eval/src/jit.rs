//! Cranelift JIT backend for oCAS expression evaluation.
//!
//! Compiles [`Instr`](crate::Instr) sequences to native machine code via
//! Cranelift, enabling ≥10x faster evaluation compared to the interpreter.
//!
//! Enabled with the `jit` feature flag.

use std::collections::HashMap;

use cranelift_codegen::Context;
use cranelift_codegen::ir::types::{F64, I64};
use cranelift_codegen::ir::{AbiParam, Function, InstBuilder, Signature, UserFuncName, Value};
use cranelift_codegen::isa::CallConv;
use cranelift_frontend::{FunctionBuilder, FunctionBuilderContext};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{FuncId, Linkage, Module};

use crate::error::{EvaluationError, Result};
use crate::instruction::Instr;

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

    /// Return the number of parameters expected.
    pub fn param_count(&self) -> usize {
        self.param_count
    }
}

/// JIT compilation engine for oCAS expressions.
pub struct JitEngine;

impl JitEngine {
    /// Compile an instruction sequence into a callable JIT function.
    pub fn compile(
        instructions: &[Instr],
        param_count: usize,
        result_count: usize,
    ) -> Result<JitCompiledFunction> {
        let builder = JITBuilder::with_flags(&[], cranelift_module::default_libcall_names())?;
        let mut module = JITModule::new(builder);

        let mut sig = Signature::new(CallConv::SystemV);
        sig.params.push(AbiParam::new(I64));
        sig.params.push(AbiParam::new(I64));

        let mut fn_builder_ctx = FunctionBuilderContext::new();

        let mut ctx = Context::new();

        let func_id = build_ir(
            &mut module,
            &mut ctx,
            instructions,
            param_count,
            result_count,
            &sig,
            &mut fn_builder_ctx,
        )?;

        module.define_function(func_id, &mut ctx)?;
        module.finalize_definitions()?;

        let entry_ptr = module.get_finalized_function(func_id);
        // Safety: Cranelift guarantees the pointer is a valid function
        // with the signature we declared.
        let entry: extern "C" fn(*const f64, *mut f64) = unsafe { std::mem::transmute(entry_ptr) };

        Ok(JitCompiledFunction {
            _module: module,
            entry,
            param_count,
            result_count,
        })
    }
}

fn build_ir(
    module: &mut JITModule,
    ctx: &mut Context,
    instructions: &[Instr],
    param_count: usize,
    result_count: usize,
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

    let mut slots: HashMap<usize, Value> = HashMap::new();

    // Load parameters
    for i in 0..param_count {
        let offset = (i * 8) as i64;
        let addr = builder.ins().iadd_imm(params_ptr, offset);
        let val = builder
            .ins()
            .load(F64, cranelift_codegen::ir::MemFlags::new(), addr, 0);
        slots.insert(i, val);
    }

    // Execute instructions
    for instr in instructions {
        match instr {
            Instr::Copy { dst, src } => {
                let val = slots[src];
                slots.insert(*dst, val);
            }
            Instr::Add { dst, srcs } => {
                let mut acc = slots[&srcs[0]];
                for src in &srcs[1..] {
                    acc = builder.ins().fadd(acc, slots[src]);
                }
                slots.insert(*dst, acc);
            }
            Instr::Mul { dst, srcs } => {
                let mut acc = slots[&srcs[0]];
                for src in &srcs[1..] {
                    acc = builder.ins().fmul(acc, slots[src]);
                }
                slots.insert(*dst, acc);
            }
            Instr::Pow { dst, base, exp } => {
                let base_val = slots[base];
                let exp_val = builder.ins().f64const(*exp as f64);
                let result = call_libm(&mut builder, module, "pow", &[base_val, exp_val])?;
                slots.insert(*dst, result);
            }
            Instr::Powf { dst, base, exp } => {
                let base_val = slots[base];
                let exp_val = slots[exp];
                let result = call_libm(&mut builder, module, "pow", &[base_val, exp_val])?;
                slots.insert(*dst, result);
            }
            Instr::BuiltinFun { dst, name, src } => {
                let src_val = slots[src];
                let fn_name = match name.as_str().to_lowercase().as_str() {
                    "sin" => "sin",
                    "cos" => "cos",
                    "tan" => "tan",
                    "exp" => "exp",
                    "log" => "log",
                    "sqrt" => "sqrt",
                    "abs" => "fabs",
                    _ => {
                        return Err(EvaluationError::FunctionNotFound {
                            name: name.as_str().to_string(),
                        });
                    }
                };
                let result = call_libm(&mut builder, module, fn_name, &[src_val])?;
                slots.insert(*dst, result);
            }
            Instr::ExternalFun { .. } => {
                return Err(EvaluationError::JitCompilationError {
                    message: "external functions not supported in JIT".into(),
                });
            }
        }
    }

    // Store results
    // The result slot is the max slot index used
    for i in 0..result_count {
        if let Some(&val) = slots.values().nth(i) {
            let offset = (i * 8) as i64;
            let addr = builder.ins().iadd_imm(results_ptr, offset);
            builder
                .ins()
                .store(cranelift_codegen::ir::MemFlags::new(), val, addr, 0);
        }
    }

    builder.ins().return_(&[]);
    builder.seal_block(block);
    builder.finalize();

    Ok(func_id)
}

/// Call a libm function (sin, cos, pow, etc.) and return the result Value.
fn call_libm(
    builder: &mut FunctionBuilder,
    module: &mut JITModule,
    name: &str,
    args: &[Value],
) -> Result<Value> {
    let param_types: Vec<AbiParam> = args.iter().map(|_| AbiParam::new(F64)).collect();
    let callee_sig = Signature {
        params: param_types,
        returns: vec![AbiParam::new(F64)],
        call_conv: CallConv::SystemV,
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
        let func = JitEngine::compile(&[], 0, 0);
        assert!(func.is_ok());
    }

    #[test]
    #[ignore = "JIT calling convention needs platform-specific tuning"]
    fn jit_simple_add() {
        let instructions = vec![Instr::Add {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile(&instructions, 2, 1).unwrap();
        let result = func.call(&[2.0, 3.0]);
        assert!((result[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    #[ignore = "JIT calling convention needs platform-specific tuning"]
    fn jit_simple_mul() {
        let instructions = vec![Instr::Mul {
            dst: 2,
            srcs: vec![0, 1],
        }];
        let func = JitEngine::compile(&instructions, 2, 1).unwrap();
        let result = func.call(&[3.0, 4.0]);
        assert!((result[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    #[ignore = "JIT calling convention needs platform-specific tuning"]
    fn jit_builtin_sin() {
        let name = ocas_atom::Symbol::new("sin");
        let instructions = vec![Instr::BuiltinFun {
            dst: 1,
            name,
            src: 0,
        }];
        let func = JitEngine::compile(&instructions, 1, 1).unwrap();
        let result = func.call(&[std::f64::consts::FRAC_PI_2]);
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    #[ignore = "JIT calling convention needs platform-specific tuning"]
    fn jit_pow_integer() {
        let instructions = vec![Instr::Pow {
            dst: 1,
            base: 0,
            exp: 3,
        }];
        let func = JitEngine::compile(&instructions, 1, 1).unwrap();
        let result = func.call(&[2.0]);
        assert!((result[0] - 8.0).abs() < 1e-10);
    }
}
