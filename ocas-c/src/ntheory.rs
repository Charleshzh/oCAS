//! C/C++ bindings for number-theory functions.
//!
//! Integers of arbitrary size are passed as decimal strings. String results
//! are heap-allocated and must be released with
//! [`ocas_string_free`](crate::expression::ocas_string_free).

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{CString, c_char, c_int};
use std::ptr;

use num_bigint::BigInt;
use ocas_domain::Integer;
use ocas_domain::number_theory::crt::crt_many;
use ocas_domain::number_theory::dlog::{dlog_bsgs, dlog_pohlig_hellman};
use ocas_domain::number_theory::factor::factor_integer;
use ocas_domain::number_theory::functions::{
    divisor_sigma, divisor_tau, euler_phi, liouville_lambda, moebius_mu,
};
use ocas_domain::number_theory::primes::is_prime_bpsw;
use ocas_domain::number_theory::{jacobi, next_prime};

use crate::error::{OCAS_ERROR_PARSE, OCAS_ERROR_RUNTIME, set};
use crate::expression::cstr_to_str_pub;

/// Parse a required decimal-integer argument.
fn parse_int_arg(s: Option<&str>) -> Option<Integer> {
    let s = s?;
    match s.parse::<BigInt>() {
        Ok(v) => Some(Integer::from(v)),
        Err(_) => {
            set(OCAS_ERROR_PARSE, "expected a decimal integer");
            None
        }
    }
}

macro_rules! int_arg {
    ($ptr:ident, $name:literal, $err_out:ident) => {
        match parse_int_arg(cstr_to_str_pub($ptr, $name)) {
            Some(v) => v,
            None => {
                crate::error::write_last_code($err_out);
                return ptr::null_mut();
            }
        }
    };
}

/// Factor `|n|` into primes, returned as `"p1:e1,p2:e2,..."` (ascending).
/// A leading `"-1:1"` pair marks negative input. Returns `NULL` on failure.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string. `err_out` if
/// non-null must point to writable memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_factorint(
    n: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let n = int_arg!(n, "n", err_out);
    let result = std::panic::catch_unwind(|| {
        let mut parts: Vec<String> = Vec::new();
        if n.is_negative() {
            parts.push("-1:1".to_string());
        }
        parts.extend(
            factor_integer(&n)
                .into_iter()
                .map(|(p, e)| format!("{p}:{e}")),
        );
        parts.join(",")
    });
    finish_string(result, err_out)
}

/// BPSW probable-prime test: returns 1 when `n` is (probably) prime, 0 when
/// composite, and −1 on error (with `err_out` set).
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_isprime(n: *const c_char, err_out: *mut c_int) -> c_int {
    crate::error::clear();
    let Some(n) = parse_int_arg(cstr_to_str_pub(n, "n")) else {
        crate::error::write_last_code(err_out);
        return -1;
    };
    let result = std::panic::catch_unwind(|| is_prime_bpsw(&n));
    match result {
        Ok(b) => {
            if !err_out.is_null() {
                unsafe { *err_out = crate::error::OCAS_OK };
            }
            i32::from(b)
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during primality test");
            crate::error::write_last_code(err_out);
            -1
        }
    }
}

/// Smallest prime strictly greater than `n`, as a decimal string.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_nextprime(
    n: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let n = int_arg!(n, "n", err_out);
    let result = std::panic::catch_unwind(|| next_prime(&n).to_string());
    finish_string(result, err_out)
}

/// Solve `base^x ≡ target (mod p)`: Pohlig–Hellman for prime `p`, BSGS
/// otherwise. Returns the logarithm as a decimal string, or `NULL` (with
/// `err_out = OCAS_ERROR_RUNTIME`) when no logarithm exists.
///
/// # Safety
///
/// All arguments must be valid null-terminated decimal strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_discrete_log(
    p: *const c_char,
    base: *const c_char,
    target: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let p = int_arg!(p, "p", err_out);
    let base = int_arg!(base, "base", err_out);
    let target = int_arg!(target, "target", err_out);
    let result = std::panic::catch_unwind(|| {
        let x = if is_prime_bpsw(&p) {
            dlog_pohlig_hellman(&base, &target, &p)
        } else {
            dlog_bsgs(&base, &target, &p)
        };
        match x {
            Some(x) => x.to_string(),
            None => String::new(),
        }
    });
    match result {
        Ok(s) if s.is_empty() => {
            set(OCAS_ERROR_RUNTIME, "no discrete logarithm exists");
            crate::error::write_last_code(err_out);
            ptr::null_mut()
        }
        other => finish_string(other, err_out),
    }
}

/// Chinese remainder theorem over comma-separated decimal lists:
/// `moduli = "m1,m2,..."`, `residues = "r1,r2,..."`. Returns `"r,m"` with
/// `r ≡ residues[i] (mod moduli[i])`, or `NULL` for inconsistent systems or
/// malformed input.
///
/// # Safety
///
/// `moduli` and `residues` must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_crt(
    moduli: *const c_char,
    residues: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let parse_list = |s: Option<&str>| -> Option<Vec<Integer>> {
        s?.split(',')
            .map(|t| t.trim().parse::<BigInt>().ok().map(Integer::from))
            .collect()
    };
    let (Some(ms), Some(rs)) = (
        parse_list(cstr_to_str_pub(moduli, "moduli")),
        parse_list(cstr_to_str_pub(residues, "residues")),
    ) else {
        set(
            OCAS_ERROR_PARSE,
            "expected comma-separated decimal integers",
        );
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    };
    if ms.len() != rs.len() {
        set(
            OCAS_ERROR_RUNTIME,
            "moduli and residues must have the same length",
        );
        crate::error::write_last_code(err_out);
        return ptr::null_mut();
    }
    let result = std::panic::catch_unwind(|| {
        let cs: Vec<(Integer, Integer)> = rs.into_iter().zip(ms).collect();
        crt_many(&cs).map(|(r, m)| format!("{r},{m}"))
    });
    match result {
        Ok(Some(s)) => finish_string(Ok(s), err_out),
        Ok(None) => {
            set(OCAS_ERROR_RUNTIME, "inconsistent CRT system");
            crate::error::write_last_code(err_out);
            ptr::null_mut()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during ntheory operation");
            crate::error::write_last_code(err_out);
            ptr::null_mut()
        }
    }
}

/// The Jacobi symbol `(a / n)` for odd positive `n`: returns −1, 0, or 1;
/// −2 with `err_out` set signals invalid input.
///
/// # Safety
///
/// `a` and `n` must be valid null-terminated decimal strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_jacobi(
    a: *const c_char,
    n: *const c_char,
    err_out: *mut c_int,
) -> c_int {
    crate::error::clear();
    let (Some(a), Some(n)) = (
        parse_int_arg(cstr_to_str_pub(a, "a")),
        parse_int_arg(cstr_to_str_pub(n, "n")),
    ) else {
        crate::error::write_last_code(err_out);
        return -2;
    };
    if n.is_zero() || n.is_negative() || n.is_even() {
        set(OCAS_ERROR_RUNTIME, "n must be a positive odd integer");
        crate::error::write_last_code(err_out);
        return -2;
    }
    let result = std::panic::catch_unwind(|| jacobi(&a, &n));
    match result {
        Ok(v) => {
            if !err_out.is_null() {
                unsafe { *err_out = crate::error::OCAS_OK };
            }
            c_int::from(v)
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during Jacobi symbol");
            crate::error::write_last_code(err_out);
            -2
        }
    }
}

/// Euler's totient `φ(n)` as a decimal string.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_totient(
    n: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let n = int_arg!(n, "n", err_out);
    let result = std::panic::catch_unwind(|| euler_phi(&n).to_string());
    finish_string(result, err_out)
}

/// The Möbius function `μ(n)`: returns −1, 0, or 1; −2 with `err_out` set
/// signals invalid input.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_mobius(n: *const c_char, err_out: *mut c_int) -> c_int {
    crate::error::clear();
    let Some(n) = parse_int_arg(cstr_to_str_pub(n, "n")) else {
        crate::error::write_last_code(err_out);
        return -2;
    };
    let result = std::panic::catch_unwind(|| moebius_mu(&n));
    match result {
        Ok(v) => {
            if !err_out.is_null() {
                unsafe { *err_out = crate::error::OCAS_OK };
            }
            c_int::from(v)
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during Möbius function");
            crate::error::write_last_code(err_out);
            -2
        }
    }
}

/// Number of positive divisors `τ(n)` as a decimal string.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_divisor_count(
    n: *const c_char,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let n = int_arg!(n, "n", err_out);
    let result = std::panic::catch_unwind(|| divisor_tau(&n).to_string());
    finish_string(result, err_out)
}

/// Sum of `k`-th powers of the positive divisors `σ_k(n)` as a decimal
/// string.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_divisor_sigma(
    n: *const c_char,
    k: u32,
    err_out: *mut c_int,
) -> *mut c_char {
    crate::error::clear();
    let n = int_arg!(n, "n", err_out);
    let result = std::panic::catch_unwind(|| divisor_sigma(&n, k).to_string());
    finish_string(result, err_out)
}

/// Liouville's function `λ(n) = (−1)^Ω(n)`: returns −1, 0, or 1; −2 with
/// `err_out` set signals invalid input.
///
/// # Safety
///
/// `n` must be a valid null-terminated decimal string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ocas_ntheory_liouville(n: *const c_char, err_out: *mut c_int) -> c_int {
    crate::error::clear();
    let Some(n) = parse_int_arg(cstr_to_str_pub(n, "n")) else {
        crate::error::write_last_code(err_out);
        return -2;
    };
    let result = std::panic::catch_unwind(|| liouville_lambda(&n));
    match result {
        Ok(v) => {
            if !err_out.is_null() {
                unsafe { *err_out = crate::error::OCAS_OK };
            }
            c_int::from(v)
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during Liouville function");
            crate::error::write_last_code(err_out);
            -2
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Convert a panic-catching result into a heap-allocated C string.
fn finish_string(result: std::thread::Result<String>, err_out: *mut c_int) -> *mut c_char {
    let s = match result {
        Ok(s) => s,
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "panic during ntheory operation");
            if !err_out.is_null() {
                unsafe { *err_out = OCAS_ERROR_RUNTIME };
            }
            return ptr::null_mut();
        }
    };
    match CString::new(s) {
        Ok(cs) => {
            if !err_out.is_null() {
                unsafe { *err_out = crate::error::OCAS_OK };
            }
            cs.into_raw()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "result string contains NUL");
            if !err_out.is_null() {
                unsafe { *err_out = OCAS_ERROR_RUNTIME };
            }
            ptr::null_mut()
        }
    }
}
