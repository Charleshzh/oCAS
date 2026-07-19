//! AST-to-instruction compiler.
//!
//! Transforms an [`EvalTree`](crate::EvalTree) into an
//! [`ExpressionEvaluator`](crate::ExpressionEvaluator) by generating
//! a sequence of [`Instr`](crate::Instr)s.
//!
//! The compiler walks the tree in post-order, assigning stack slots
//! and emitting instructions for each node.

use std::collections::HashSet;

use ocas_atom::Atom;
use ocas_core::FastHashMap as HashMap;

use crate::domain::{EvaluationDomain, PowfExtension};
use crate::error::{EvaluationError, Result};
use crate::evaluator::ExpressionEvaluator;
use crate::function_map::FunctionMap;
use crate::instruction::Instr;
use crate::optimize;
use crate::tree::EvalTree;

/// Compile an [`Atom`] into an [`ExpressionEvaluator`].
pub fn compile_atom<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
) -> Result<ExpressionEvaluator<T>> {
    compile_atom_with(atom, None)
}

/// Compile an [`Atom`] into an [`ExpressionEvaluator`] with a function map.
pub fn compile_atom_with<T: EvaluationDomain + PowfExtension>(
    atom: Atom<'_>,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>> {
    let tree = EvalTree::from_atom(atom);
    compile_tree_with(&tree, function_map)
}

/// Compile an [`EvalTree`] into an [`ExpressionEvaluator`].
#[allow(dead_code)]
pub fn compile_tree<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
) -> Result<ExpressionEvaluator<T>> {
    compile_tree_with(tree, None)
}

/// Compile an [`EvalTree`] with an optional function map.
pub fn compile_tree_with<T: EvaluationDomain + PowfExtension>(
    tree: &EvalTree,
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>> {
    compile_trees(&[tree], function_map)
}

/// Compile multiple [`Atom`]s into a single multi-output
/// [`ExpressionEvaluator`], sharing constants, CSE, and stack slots
/// across all outputs.
pub fn compile_atoms_multi<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
) -> Result<ExpressionEvaluator<T>> {
    compile_atoms_multi_with(atoms, None)
}

/// Compile multiple [`Atom`]s with a function map into a single
/// multi-output [`ExpressionEvaluator`].
pub fn compile_atoms_multi_with<T: EvaluationDomain + PowfExtension>(
    atoms: &[Atom<'_>],
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>> {
    let trees: Vec<EvalTree> = atoms.iter().map(|a| EvalTree::from_atom(*a)).collect();
    let refs: Vec<&EvalTree> = trees.iter().collect();
    compile_trees(&refs, function_map)
}

/// Compile multiple [`EvalTree`]s into a single multi-output
/// [`ExpressionEvaluator`].
pub fn compile_trees_multi<T: EvaluationDomain + PowfExtension>(
    trees: &[&EvalTree],
) -> Result<ExpressionEvaluator<T>> {
    compile_trees(trees, None)
}

/// Shared compiler core: compile a set of trees into one evaluator
/// with multiple result slots.
fn compile_trees<T: EvaluationDomain + PowfExtension>(
    trees: &[&EvalTree],
    function_map: Option<FunctionMap<T>>,
) -> Result<ExpressionEvaluator<T>> {
    // Constant folding and algebraic simplification at the tree level.
    let folded: Vec<EvalTree> = trees.iter().map(|t| t.fold_constants()).collect();

    // Pass 1: collect all variable names and count constants
    let mut var_names = HashSet::new();
    let mut const_count = 0usize;
    for tree in &folded {
        scan_tree(tree, &mut var_names, &mut const_count);
    }
    let param_count = var_names.len();

    // Assign parameter slots: sort for deterministic ordering
    let mut sorted_vars: Vec<String> = var_names.into_iter().collect();
    sorted_vars.sort();
    let var_to_param: HashMap<String, usize> = sorted_vars
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();

    // Pass 2: compile all trees with one shared context so that CSE
    // deduplicates subexpressions across outputs.
    let temp_base = param_count + const_count;
    let (instructions, next_temp, constants, result_slots) = {
        let mut ctx =
            CompileContext::<T>::new(param_count, temp_base, var_to_param, function_map.as_ref());
        let mut result_slots = Vec::with_capacity(folded.len());
        for tree in &folded {
            result_slots.push(ctx.compile_node(tree)?);
        }
        (ctx.instructions, ctx.next_temp, ctx.constants, result_slots)
    };

    let actual_const_count = constants.len();
    let (instructions, temp_count, result_indices) =
        optimize::optimize(instructions, temp_base, next_temp, &result_slots);
    let stack_size = temp_base + temp_count;

    match function_map {
        Some(fm) => Ok(ExpressionEvaluator::new_with_functions(
            instructions,
            param_count,
            actual_const_count,
            stack_size,
            result_indices,
            constants,
            fm,
        )),
        None => Ok(ExpressionEvaluator::new(
            instructions,
            param_count,
            actual_const_count,
            stack_size,
            result_indices,
            constants,
        )),
    }
}

/// Pre-scan tree to count variables and constants.
fn scan_tree(tree: &EvalTree, vars: &mut HashSet<String>, const_count: &mut usize) {
    match tree {
        EvalTree::Num(_) => {
            *const_count += 1;
        }
        EvalTree::Var(name) => {
            vars.insert(name.clone());
        }
        EvalTree::Add(terms) | EvalTree::Mul(terms) => {
            for t in terms {
                scan_tree(t, vars, const_count);
            }
        }
        EvalTree::Pow(base, exp) => {
            scan_tree(base, vars, const_count);
            scan_tree(exp, vars, const_count);
        }
        EvalTree::Fun(_, args) => {
            for a in args {
                scan_tree(a, vars, const_count);
            }
        }
    }
}

struct CompileContext<'a, T: EvaluationDomain> {
    instructions: Vec<Instr>,
    /// Next available temp slot index. Temps start at `temp_base` in the actual stack.
    next_temp: usize,
    /// Parameter slots occupy stack[0..param_count].
    param_count: usize,
    /// Base index for temp slots in the actual stack (= param_count + estimated_const_count).
    temp_base: usize,
    /// variable name → stack slot index (0..param_count-1)
    variables: HashMap<String, usize>,
    /// constant values in order
    constants: Vec<T>,
    /// Optional function map for resolving external functions
    function_map: Option<&'a FunctionMap<T>>,
}

impl<'a, T: EvaluationDomain> CompileContext<'a, T> {
    fn new(
        param_count: usize,
        temp_base: usize,
        variables: HashMap<String, usize>,
        function_map: Option<&'a FunctionMap<T>>,
    ) -> Self {
        Self {
            instructions: Vec::new(),
            next_temp: 0,
            param_count,
            temp_base,
            variables,
            constants: Vec::new(),
            function_map,
        }
    }

    fn alloc_temp(&mut self) -> usize {
        let slot = self.next_temp + self.temp_base;
        self.next_temp += 1;
        slot
    }

    fn param_slot(&self, name: &str) -> usize {
        self.variables[name]
    }

    fn const_slot(&mut self, value: T) -> usize {
        let idx = self.constants.len();
        self.constants.push(value);
        self.param_count + idx
    }

    fn compile_node(&mut self, node: &EvalTree) -> Result<usize> {
        match node {
            EvalTree::Num(n) => {
                let dst = self.alloc_temp();
                let const_slot = self.const_slot(T::from_f64(*n));
                self.instructions.push(Instr::Copy {
                    dst,
                    src: const_slot,
                });
                Ok(dst)
            }
            EvalTree::Var(name) => {
                let dst = self.alloc_temp();
                let param_slot = self.param_slot(name);
                self.instructions.push(Instr::Copy {
                    dst,
                    src: param_slot,
                });
                Ok(dst)
            }
            EvalTree::Add(terms) => {
                let dst = self.alloc_temp();
                let mut srcs = Vec::with_capacity(terms.len());
                for term in terms {
                    srcs.push(self.compile_node(term)?);
                }
                self.instructions.push(Instr::Add { dst, srcs });
                Ok(dst)
            }
            EvalTree::Mul(factors) => {
                let dst = self.alloc_temp();
                let mut srcs = Vec::with_capacity(factors.len());
                for factor in factors {
                    srcs.push(self.compile_node(factor)?);
                }
                self.instructions.push(Instr::Mul { dst, srcs });
                Ok(dst)
            }
            EvalTree::Pow(base, exp) => {
                let base_slot = self.compile_node(base)?;
                let dst = self.alloc_temp();
                if let EvalTree::Num(n) = exp.as_ref()
                    && n.fract() == 0.0
                    && *n >= i64::MIN as f64
                    && *n <= i64::MAX as f64
                {
                    self.instructions.push(Instr::Pow {
                        dst,
                        base: base_slot,
                        exp: *n as i64,
                    });
                    return Ok(dst);
                }
                let exp_slot = self.compile_node(exp)?;
                self.instructions.push(Instr::Powf {
                    dst,
                    base: base_slot,
                    exp: exp_slot,
                });
                Ok(dst)
            }
            EvalTree::Fun(name, args) => {
                if is_builtin(name) && args.len() == 1 {
                    let arg_slot = self.compile_node(&args[0])?;
                    let dst = self.alloc_temp();
                    let op = crate::instruction::BuiltinOp::from_name(name)
                        .expect("is_builtin guarantees known name");
                    self.instructions.push(Instr::BuiltinOp {
                        dst,
                        op,
                        src: arg_slot,
                    });
                    Ok(dst)
                } else if let Some(fm) = self.function_map {
                    // Look up in function map
                    if let Some(_entry) = fm.resolve(name) {
                        let mut srcs = Vec::with_capacity(args.len());
                        for arg in args {
                            srcs.push(self.compile_node(arg)?);
                        }
                        let dst = self.alloc_temp();
                        // Find the function index in the map
                        let fn_idx =
                            fm.index_of(name)
                                .ok_or_else(|| EvaluationError::FunctionNotFound {
                                    name: name.clone(),
                                })?;
                        self.instructions
                            .push(Instr::ExternalFun { dst, fn_idx, srcs });
                        Ok(dst)
                    } else {
                        Err(EvaluationError::FunctionNotFound { name: name.clone() })
                    }
                } else {
                    Err(EvaluationError::FunctionNotFound { name: name.clone() })
                }
            }
        }
    }
}

fn is_builtin(name: &str) -> bool {
    matches!(
        name.to_lowercase().as_str(),
        "sin" | "cos" | "tan" | "sec" | "csc" | "cot" | "exp" | "log" | "sqrt" | "abs"
    )
}

// ---------------------------------------------------------------------------
// ExpressionEvaluator::compile
// ---------------------------------------------------------------------------

impl<T: EvaluationDomain + PowfExtension> ExpressionEvaluator<T> {
    /// Compile an [`Atom`] into an executable evaluator.
    pub fn compile(atom: Atom<'_>) -> Result<Self> {
        compile_atom(atom)
    }

    /// Compile an [`Atom`] with a [`FunctionMap`] for user-defined functions.
    pub fn compile_with(atom: Atom<'_>, map: FunctionMap<T>) -> Result<Self> {
        compile_atom_with(atom, Some(map))
    }

    /// Compile multiple [`Atom`]s into one multi-output evaluator,
    /// sharing common subexpressions across all outputs.
    pub fn compile_multi(atoms: &[Atom<'_>]) -> Result<Self> {
        compile_atoms_multi(atoms)
    }

    /// Compile multiple [`Atom`]s with a [`FunctionMap`] into one
    /// multi-output evaluator.
    pub fn compile_multi_with(atoms: &[Atom<'_>], map: FunctionMap<T>) -> Result<Self> {
        compile_atoms_multi_with(atoms, Some(map))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    #[test]
    fn compile_constant() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.num(42);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[]).unwrap();
        assert!((result[0] - 42.0).abs() < 1e-10);
    }

    #[test]
    fn compile_single_var() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.var("x");
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        assert_eq!(eval.param_count(), 1);
        let result = eval.evaluate(&[7.0]).unwrap();
        assert!((result[0] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn compile_add_two_vars() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.add(&[ctx.var("x"), ctx.var("y")]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[2.0, 3.0]).unwrap();
        assert!((result[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn compile_mul_var_const() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.mul(&[ctx.var("x"), ctx.num(3)]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[4.0]).unwrap();
        assert!((result[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn compile_pow_integer_exp() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.pow(ctx.var("x"), ctx.num(3));
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[2.0]).unwrap();
        assert!((result[0] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn compile_sin() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("sin", &[ctx.var("x")]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[std::f64::consts::FRAC_PI_2]).unwrap();
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn compile_cos() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("cos", &[ctx.var("x")]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[std::f64::consts::PI]).unwrap();
        assert!((result[0] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn compile_exp_log_roundtrip() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let exp_x = ctx.fun("exp", &[ctx.var("x")]);
        let expr = ctx.fun("log", &[exp_x]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[2.0]).unwrap();
        assert!((result[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn compile_nested_expression() {
        // (x + 1) * (x - 1) = x^2 - 1
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let x_plus_1 = ctx.add(&[x, ctx.num(1)]);
        let x_minus_1 = ctx.add(&[x, ctx.num(-1)]);
        let expr = ctx.mul(&[x_plus_1, x_minus_1]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[3.0]).unwrap();
        assert!((result[0] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn compile_sqrt() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("sqrt", &[ctx.num(16)]);
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        let result = eval.evaluate(&[]).unwrap();
        assert!((result[0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn compile_zero_params() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("sin", &[ctx.num(1)]); // sin(1 rad)
        let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
        assert_eq!(eval.param_count(), 0);
        let result = eval.evaluate(&[]).unwrap();
        assert!((result[0] - 1.0f64.sin()).abs() < 1e-10);
    }

    #[test]
    fn compile_with_external_function() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("square", &[ctx.var("x")]);

        let mut map = FunctionMap::<f64>::new();
        map.register("square", 1, Box::new(|args| args[0] * args[0]));

        let eval = ExpressionEvaluator::compile_with(expr, map).unwrap();
        let result = eval.evaluate(&[3.0]).unwrap();
        assert!((result[0] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn compile_external_function_not_registered() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("missing_fn", &[ctx.var("x")]);

        let result: Result<ExpressionEvaluator<f64>> = ExpressionEvaluator::compile(expr);
        assert!(result.is_err());
    }

    #[test]
    fn compile_with_case_insensitive_external() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = ctx.fun("Square", &[ctx.num(4)]);

        let mut map = FunctionMap::<f64>::new();
        map.register("square", 1, Box::new(|args| args[0] * args[0]));

        let eval = ExpressionEvaluator::compile_with(expr, map).unwrap();
        let result = eval.evaluate(&[]).unwrap();
        assert!((result[0] - 16.0).abs() < 1e-10);
    }

    #[test]
    fn compile_multi_two_outputs() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let sum = ctx.add(&[ctx.var("x"), ctx.var("y")]);
        let prod = ctx.mul(&[ctx.var("x"), ctx.var("y")]);
        let eval: ExpressionEvaluator<f64> =
            ExpressionEvaluator::compile_multi(&[sum, prod]).unwrap();
        assert_eq!(eval.result_count(), 2);
        assert_eq!(eval.param_count(), 2);
        let result = eval.evaluate(&[2.0, 3.0]).unwrap();
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn compile_multi_shared_subexpression() {
        // outputs: sin(x) + 1, sin(x) * 2 — sin(x) shared via CSE
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let sin_x = ctx.fun("sin", &[ctx.var("x")]);
        let out0 = ctx.add(&[sin_x, ctx.num(1)]);
        let out1 = ctx.mul(&[sin_x, ctx.num(2)]);
        let eval: ExpressionEvaluator<f64> =
            ExpressionEvaluator::compile_multi(&[out0, out1]).unwrap();
        let result = eval.evaluate(&[std::f64::consts::FRAC_PI_2]).unwrap();
        assert!((result[0] - 2.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn compile_multi_constant_folding() {
        // (2 + 3) * x and x^1 — folding reduces both
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let five_x = ctx.mul(&[ctx.add(&[ctx.num(2), ctx.num(3)]), ctx.var("x")]);
        let x_pow_1 = ctx.pow(ctx.var("x"), ctx.num(1));
        let eval: ExpressionEvaluator<f64> =
            ExpressionEvaluator::compile_multi(&[five_x, x_pow_1]).unwrap();
        let result = eval.evaluate(&[4.0]).unwrap();
        assert!((result[0] - 20.0).abs() < 1e-10);
        assert!((result[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn compile_multi_with_external_function() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let sq = ctx.fun("square", &[ctx.var("x")]);
        let cube_arg = ctx.fun("square", &[ctx.var("x")]);
        let out1 = ctx.mul(&[cube_arg, ctx.var("x")]);

        let mut map = FunctionMap::<f64>::new();
        map.register("square", 1, Box::new(|args| args[0] * args[0]));

        let eval = ExpressionEvaluator::compile_multi_with(&[sq, out1], map).unwrap();
        let result = eval.evaluate(&[3.0]).unwrap();
        assert!((result[0] - 9.0).abs() < 1e-10);
        assert!((result[1] - 27.0).abs() < 1e-10);
    }
}
