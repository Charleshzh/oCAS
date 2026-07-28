//! Python `ode` module — ordinary differential equation solvers.

use ocas_atom::Symbol;
use ocas_calc::ode::{
    ODE, ODESolution, ODEType, classify_ode as rs_classify, dsolve as rs_dsolve,
    dsolve_ivp as rs_dsolve_ivp,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::expression::Expression;

/// Classify an ODE and return the applicable method names.
///
/// `equation` is the ODE written as an expression equal to zero (e.g.
/// `"Derivative(y(x), x) - y(x)"`), `func` the unknown function name
/// (e.g. `"y"`), and `var` the independent variable (e.g. `"x"`).
///
/// Returns a list of method names such as `"LinearFirst"`, `"Separable"`,
/// `"LinearConstantCoeff"`, etc.
#[pyfunction]
#[pyo3(name = "classify_ode")]
pub fn py_classify_ode(equation: &Expression, func: &str, var: &str) -> PyResult<Vec<String>> {
    let (ctx, ode) = build_ode(equation, func, var)?;
    Ok(rs_classify(ctx, ode)
        .iter()
        .map(ode_type_name)
        .map(str::to_owned)
        .collect())
}

/// Solve an ODE symbolically.
///
/// - `equation`: expression equal to zero, e.g. `"Derivative(y(x), x) - y(x)"`.
/// - `func`: unknown function name, e.g. `"y"`.
/// - `var`: independent variable, e.g. `"x"`.
/// - `hint`: optional method name (one of the strings returned by
///   `classify_ode`) to force a specific solver.
///
/// Returns a string describing the solution: an explicit solution
/// `y = ...`, an implicit form, a truncated series, or the unevaluated ODE.
#[pyfunction]
#[pyo3(name = "dsolve", signature = (equation, func, var, hint=None))]
pub fn py_dsolve(
    equation: &Expression,
    func: &str,
    var: &str,
    hint: Option<&str>,
) -> PyResult<String> {
    let hint_type = hint.map(parse_ode_type).transpose()?;
    let (ctx, ode) = build_ode(equation, func, var)?;
    let sol = rs_dsolve(ctx, ode, hint_type);
    Ok(format_solution(&sol))
}

/// Solve a first- or second-order linear constant-coefficient IVP via the
/// Laplace transform.
///
/// - `y0`: value `y(0)` as a string expression (e.g. `"1"`).
/// - `y1`: value `y'(0)` (required for second-order problems).
///
/// Returns an explicit solution string with no free constants.
#[pyfunction]
#[pyo3(name = "dsolve_ivp", signature = (equation, func, var, y0, y1=None))]
pub fn py_dsolve_ivp(
    equation: &Expression,
    func: &str,
    var: &str,
    y0: &str,
    y1: Option<&str>,
) -> PyResult<String> {
    let (ctx, ode) = build_ode(equation, func, var)?;
    let y0_atom = parse_in(ctx, y0)?;
    let y1_atom = y1.map(|s| parse_in(ctx, s)).transpose()?;
    let sol = rs_dsolve_ivp(ctx, ode, y0_atom, y1_atom);
    Ok(format_solution(&sol))
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build an `ODE` inside the expression's arena, returning the arena
/// reference and the constructed ODE.
fn build_ode(
    equation: &Expression,
    func: &str,
    var: &str,
) -> PyResult<(&'static ocas_atom::AtomArena<'static>, ODE<'static>)> {
    let ctx = equation.ctx_ref();
    let x = ctx.var(var);
    let func_atom = ctx.fun(func, &[x]);
    let ode = ODE {
        equation: equation.atom(),
        func: func_atom,
        var: Symbol::new(var),
    };
    Ok((ctx, ode))
}

fn parse_in<'a>(ctx: &'a ocas_atom::AtomArena<'a>, input: &str) -> PyResult<ocas_atom::Atom<'a>> {
    let static_input = unsafe { crate::expression::extend_str_lifetime_pub(input) };
    ocas_parse::parse(ctx, static_input)
        .map_err(|e| PyValueError::new_err(format!("parse error: {e}")))
}

fn ode_type_name(t: &ODEType) -> &'static str {
    match t {
        ODEType::Separable => "Separable",
        ODEType::LinearFirst => "LinearFirst",
        ODEType::Bernoulli => "Bernoulli",
        ODEType::Exact => "Exact",
        ODEType::Homogeneous => "Homogeneous",
        ODEType::LinearConstantCoeff => "LinearConstantCoeff",
        ODEType::CauchyEuler => "CauchyEuler",
        ODEType::ReductionOfOrder => "ReductionOfOrder",
        ODEType::PowerSeries => "PowerSeries",
    }
}

fn parse_ode_type(name: &str) -> PyResult<ODEType> {
    match name {
        "Separable" => Ok(ODEType::Separable),
        "LinearFirst" => Ok(ODEType::LinearFirst),
        "Bernoulli" => Ok(ODEType::Bernoulli),
        "Exact" => Ok(ODEType::Exact),
        "Homogeneous" => Ok(ODEType::Homogeneous),
        "LinearConstantCoeff" => Ok(ODEType::LinearConstantCoeff),
        "CauchyEuler" => Ok(ODEType::CauchyEuler),
        "ReductionOfOrder" => Ok(ODEType::ReductionOfOrder),
        "PowerSeries" => Ok(ODEType::PowerSeries),
        other => Err(PyValueError::new_err(format!(
            "unknown ODE type hint: {other}"
        ))),
    }
}

fn format_solution(sol: &ODESolution<'_>) -> String {
    match sol {
        ODESolution::Explicit(e) => format!("y = {e}"),
        ODESolution::Implicit(e) => format!("{e} = C"),
        ODESolution::Parametric(a, b) => format!("x = {a}, y = {b}"),
        ODESolution::Series(e, n) => format!("series({n} terms): y = {e}"),
        ODESolution::System(comps) => {
            let joined: Vec<String> = comps
                .iter()
                .enumerate()
                .map(|(i, c)| format!("y{} = {c}", i + 1))
                .collect();
            joined.join(", ")
        }
        ODESolution::Unsolved(_) => "unsolved".to_string(),
    }
}
