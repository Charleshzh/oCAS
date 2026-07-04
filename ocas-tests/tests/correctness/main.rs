use std::io::Write;
use std::process::Command;
use std::sync::Once;

use ocas::prelude::*;
use ocas_atom::normalize::normalize;
use ocas_core::arena::Arena;
use ocas_rewrite::rules::default_rules;

static UV_CHECK: Once = Once::new();
static mut UV_AVAILABLE: bool = false;

fn ensure_uv_available() {
    UV_CHECK.call_once(|| {
        let output = Command::new("uv").arg("--version").output();
        let available = output.map(|o| o.status.success()).unwrap_or(false);
        // SAFETY: written only once under Once.
        unsafe {
            UV_AVAILABLE = available;
        }
    });
}

fn is_uv_available() -> bool {
    ensure_uv_available();
    // SAFETY: read after Once initialization.
    unsafe { UV_AVAILABLE }
}

/// Run the SymPy comparison script and return the normalized result string.
///
/// Returns `None` if `uv` is not available, so tests can be skipped gracefully.
pub fn sympy_result(task: &str, expr: &str) -> Option<String> {
    if !is_uv_available() {
        return None;
    }

    let output = Command::new("uv")
        .args([
            "run",
            "python",
            "scripts/compare_sympy.py",
            task,
            expr,
            "check",
        ])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("failed to run compare_sympy.py");

    if !output.status.success() {
        panic!(
            "compare_sympy.py failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
}

/// Run the Symbolica comparison script and return the normalized result string.
pub fn symbolica_result(task: &str, expr: &str) -> Option<String> {
    let output = Command::new("uv")
        .args(["run", "python", "scripts/compare_symbolica.py", task, expr])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("failed to run compare_symbolica.py");

    if !output.status.success() {
        panic!(
            "compare_symbolica.py failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
}

/// Parse an expression and return its raw printed form.
pub fn parse_to_string(input: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    parse(&ctx, input).unwrap().to_string()
}

/// Parse and normalize an expression.
pub fn normalize_to_string(input: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).unwrap();
    normalize(&ctx, atom).to_string()
}

/// Parse, normalize, and simplify with default rules.
pub fn simplify_to_string(input: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).unwrap();
    let normalized = normalize(&ctx, atom);
    let rules = default_rules(&ctx, &());
    simplify(&ctx, normalized, &rules, 20).to_string()
}

/// Parse an expression, compute a derivative, and return the printed form.
///
/// The result is already simplified by `diff` internally.
pub fn diff_to_string(input: &str, var: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).unwrap();
    let v = Symbol::new(var);
    diff(&ctx, atom, v).to_string()
}

/// Parse an expression, compute an integral, and return the printed form.
pub fn integrate_to_string(input: &str, var: &str) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).unwrap();
    let v = Symbol::new(var);
    integrate(&ctx, atom, v).to_string()
}

/// Parse an expression, compute a Taylor series, and return the printed form.
pub fn taylor_to_string(input: &str, var: &str, point: i64, order: usize) -> String {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let atom = parse(&ctx, input).unwrap();
    let v = Symbol::new(var);
    let p = ctx.num(point);
    taylor(&ctx, atom, v, p, order).to_string()
}

/// Normalize an oCAS result string so it can be compared with SymPy output.
///
/// This is intentionally shallow: it removes whitespace and rewrites `^` to `**`
/// to match SymPy's exponentiation operator. Deeper semantic equivalence is
/// handled by asking SymPy to expand/collect its own output before comparison.
pub fn normalize_ocas_output(s: &str) -> String {
    s.replace(' ', "").replace('^', "**")
}

/// Run the SymPy comparison script in `verify` mode and return whether the oCAS
/// result is equivalent to the SymPy reference.
///
/// Returns `None` if `uv` is not available, so tests can be skipped gracefully.
pub fn verify_sympy_result(task: &str, expr: &str, ocas_result: &str) -> Option<bool> {
    if !is_uv_available() {
        return None;
    }

    let task_escaped = task.replace('\\', "\\\\").replace('"', "\\\"");
    let expr_escaped = expr.replace('\\', "\\\\").replace('"', "\\\"");
    let ocas_escaped = ocas_result.replace('\\', "\\\\").replace('"', "\\\"");
    let payload = format!(
        "{{\"task\": \"{}\", \"expr\": \"{}\", \"ocas_result\": \"{}\"}}",
        task_escaped, expr_escaped, ocas_escaped
    );

    let mut child = Command::new("uv")
        .args([
            "run",
            "python",
            "scripts/compare_sympy.py",
            "verify",
            task,
            expr,
        ])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn compare_sympy.py verify");

    {
        let mut stdin = child.stdin.take().expect("failed to open stdin");
        stdin
            .write_all(payload.as_bytes())
            .expect("failed to write payload");
        // stdin is closed when the handle is dropped at the end of this block.
    }

    let output = child
        .wait_with_output()
        .expect("failed to read compare_sympy.py verify output");

    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();

    if !output.status.success() {
        if stdout == "false" && stderr.is_empty() {
            // The script determined that the oCAS result is not equivalent.
            return Some(false);
        }
        eprintln!(
            "compare_sympy.py verify failed: cmd=compare_sympy.py verify {} {} {:?}",
            task, expr, ocas_result
        );
        eprintln!("stderr: {}", stderr);
        panic!("compare_sympy.py verify failed: {}", stderr);
    }

    Some(stdout == "true")
}

/// Compare an oCAS result string with the SymPy reference for a task using
/// semantic equivalence.
///
/// Panics with a descriptive message if the results are not equivalent. If `uv`
/// is not available, the assertion is silently skipped.
pub fn assert_eq_sympy(ocas_result: &str, task: &str, expr: &str) {
    if let Some(ok) = verify_sympy_result(task, expr, ocas_result) {
        assert!(
            ok,
            "oCAS result is not equivalent to SymPy for {}({})\n  oCAS: {}",
            task, expr, ocas_result
        );
    }
}

mod calculus;
mod evaluation;
mod finite_field;
mod groebner;
mod linear_solve;
mod matrix;
mod normalize;
mod parse;
mod partial_fraction;
mod poly_arithmetic;
mod poly_factor;
mod poly_gcd;
mod resultant;
mod rewrite;
mod root_isolation;
