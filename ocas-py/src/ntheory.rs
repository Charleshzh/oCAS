//! Python `ntheory` module — number-theory functions (SymPy `ntheory` style).
//!
//! Integers of arbitrary size are accepted as Python `int` (or decimal
//! strings) and returned as Python `int` where the value fits the API shape.

use num_bigint::BigInt;
use ocas_domain::Integer;
use ocas_domain::number_theory::crt::crt_many;
use ocas_domain::number_theory::dlog::{dlog_bsgs, dlog_pohlig_hellman};
use ocas_domain::number_theory::factor::factor_integer;
use ocas_domain::number_theory::functions::{
    divisor_sigma, divisor_tau, euler_phi, liouville_lambda, moebius_mu,
};
use ocas_domain::number_theory::primes::{is_prime_bpsw, is_prime_u64};
use ocas_domain::number_theory::{jacobi, next_prime};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Parse a Python `int` or decimal string into an [`Integer`].
fn to_integer(v: &Bound<'_, PyAny>) -> PyResult<Integer> {
    let s = v.str()?;
    let s = s.to_str()?;
    s.parse::<BigInt>()
        .map(Integer::from)
        .map_err(|_| PyValueError::new_err(format!("expected an integer, got {s:?}")))
}

/// Convert an [`Integer`] to a Python `int` (arbitrary precision).
fn to_py_int(py: Python<'_>, n: &Integer) -> PyResult<Py<PyAny>> {
    let s = n.to_string();
    let builtins = PyModule::import(py, "builtins")?;
    let int_cls = builtins.getattr("int")?;
    Ok(int_cls.call1((s,))?.unbind())
}

/// Factor `|n|` into primes: returns a list of `(prime, exponent)` tuples in
/// ascending order. Like SymPy's `factorint`, negative input includes the
/// `(-1, 1)` marker pair first.
#[pyfunction]
#[pyo3(name = "factorint")]
pub fn py_factorint(n: &Bound<'_, PyAny>) -> PyResult<Vec<(String, u32)>> {
    let n = to_integer(n)?;
    let mut out: Vec<(String, u32)> = Vec::new();
    if n.is_negative() {
        out.push(("-1".to_string(), 1));
    }
    out.extend(
        factor_integer(&n)
            .into_iter()
            .map(|(p, e)| (p.to_string(), e)),
    );
    Ok(out)
}

/// BPSW probable-prime test (deterministic for `n < 2^64`; no known
/// composite passes at any size).
#[pyfunction]
#[pyo3(name = "isprime")]
pub fn py_isprime(n: &Bound<'_, PyAny>) -> PyResult<bool> {
    let n = to_integer(n)?;
    Ok(is_prime_bpsw(&n))
}

/// Smallest prime strictly greater than `n`.
#[pyfunction]
#[pyo3(name = "nextprime")]
pub fn py_nextprime(py: Python<'_>, n: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let n = to_integer(n)?;
    to_py_int(py, &next_prime(&n))
}

/// Solve `base^x ≡ target (mod p)`: Pohlig–Hellman for prime `p`, BSGS
/// otherwise. Raises `ValueError` when no logarithm exists.
#[pyfunction]
#[pyo3(name = "discrete_log")]
pub fn py_discrete_log(
    py: Python<'_>,
    p: &Bound<'_, PyAny>,
    base: &Bound<'_, PyAny>,
    target: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    let p = to_integer(p)?;
    let base = to_integer(base)?;
    let target = to_integer(target)?;
    let x = if is_prime_bpsw(&p) {
        dlog_pohlig_hellman(&base, &target, &p)
    } else {
        dlog_bsgs(&base, &target, &p)
    };
    match x {
        Some(x) => to_py_int(py, &x),
        None => Err(PyValueError::new_err(
            "no discrete logarithm exists for these inputs",
        )),
    }
}

/// Chinese remainder theorem: given `moduli` and `residues` lists, return
/// `(r, m)` with `r ≡ residues[i] (mod moduli[i])`. Raises `ValueError`
/// when the system is inconsistent.
#[pyfunction]
#[pyo3(name = "crt")]
pub fn py_crt(
    py: Python<'_>,
    moduli: &Bound<'_, PyAny>,
    residues: &Bound<'_, PyAny>,
) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
    let moduli = moduli
        .cast::<pyo3::types::PySequence>()
        .map_err(|_| PyValueError::new_err("moduli must be a sequence"))?;
    let residues = residues
        .cast::<pyo3::types::PySequence>()
        .map_err(|_| PyValueError::new_err("residues must be a sequence"))?;
    if moduli.len()? != residues.len()? {
        return Err(PyValueError::new_err(
            "moduli and residues must have the same length",
        ));
    }
    let mut cs = Vec::with_capacity(moduli.len()?);
    for i in 0..moduli.len()? {
        let m = moduli.get_item(i)?;
        let r = residues.get_item(i)?;
        cs.push((to_integer(&r)?, to_integer(&m)?));
    }
    let (r, m) = crt_many(&cs).ok_or_else(|| PyValueError::new_err("inconsistent CRT system"))?;
    Ok((to_py_int(py, &r)?, to_py_int(py, &m)?))
}

/// The Jacobi symbol `(a / n)` for odd positive `n`.
#[pyfunction]
#[pyo3(name = "jacobi_symbol")]
pub fn py_jacobi_symbol(a: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<i8> {
    let a = to_integer(a)?;
    let n = to_integer(n)?;
    if n.is_zero() || n.is_negative() || n.is_even() {
        return Err(PyValueError::new_err("n must be a positive odd integer"));
    }
    Ok(jacobi(&a, &n))
}

/// Euler's totient `φ(n)`.
#[pyfunction]
#[pyo3(name = "totient")]
pub fn py_totient(py: Python<'_>, n: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let n = to_integer(n)?;
    to_py_int(py, &euler_phi(&n))
}

/// The Möbius function `μ(n)`.
#[pyfunction]
#[pyo3(name = "mobius")]
pub fn py_mobius(n: &Bound<'_, PyAny>) -> PyResult<i8> {
    let n = to_integer(n)?;
    Ok(moebius_mu(&n))
}

/// Number of positive divisors `τ(n)`.
#[pyfunction]
#[pyo3(name = "divisor_count")]
pub fn py_divisor_count(py: Python<'_>, n: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let n = to_integer(n)?;
    to_py_int(py, &divisor_tau(&n))
}

/// Sum of `k`-th powers of the positive divisors `σ_k(n)`.
#[pyfunction]
#[pyo3(name = "divisor_sigma", signature = (n, k=1))]
pub fn py_divisor_sigma(py: Python<'_>, n: &Bound<'_, PyAny>, k: u32) -> PyResult<Py<PyAny>> {
    let n = to_integer(n)?;
    to_py_int(py, &divisor_sigma(&n, k))
}

/// Liouville's function `λ(n) = (−1)^Ω(n)`.
#[pyfunction]
#[pyo3(name = "liouville_lambda")]
pub fn py_liouville_lambda(n: &Bound<'_, PyAny>) -> PyResult<i8> {
    let n = to_integer(n)?;
    Ok(liouville_lambda(&n))
}

/// Deterministic primality for `n < 2^64` (u64 fast path).
#[pyfunction]
#[pyo3(name = "isprime_u64")]
pub fn py_isprime_u64(n: u64) -> bool {
    is_prime_u64(n)
}
