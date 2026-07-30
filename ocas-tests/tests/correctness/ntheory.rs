//! Number-theory correctness tests: `ocas_domain::number_theory` vs SymPy
//! `ntheory`. Each submodule is checked on ≥ 20 cases; SymPy comparison is
//! skipped silently when `uv` is unavailable (property/value assertions
//! still run).

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

use crate::sympy_result;

fn int(s: &str) -> Integer {
    s.parse::<BigInt>().unwrap().into()
}

fn factorint_string(n: &Integer) -> String {
    let mut parts: Vec<String> = Vec::new();
    if n.is_negative() {
        parts.push("-1:1".to_string());
    }
    parts.extend(
        factor_integer(n)
            .into_iter()
            .map(|(p, e)| format!("{p}:{e}")),
    );
    parts.join(",")
}

/// Compare an oCAS string result with the SymPy reference (when available)
/// and with the hardcoded expectation (always).
fn check_nt(task: &str, input: &str, got: String, expected: &str) {
    assert_eq!(got, expected, "oCAS value wrong for {task}({input})");
    if let Some(reference) = sympy_result(task, input) {
        assert_eq!(
            got, reference,
            "mismatch vs SymPy for {task}({input}): oCAS={got} sympy={reference}"
        );
    }
}

#[test]
fn factorint_vs_sympy() {
    let cases: &[(&str, &str)] = &[
        ("0", ""),
        ("1", ""),
        ("-1", "-1:1"),
        ("2", "2:1"),
        ("12", "2:2,3:1"),
        ("-12", "-1:1,2:2,3:1"),
        ("360", "2:3,3:2,5:1"),
        ("97", "97:1"),
        ("1024", "2:10"),
        ("6561", "3:8"),
        ("561", "3:1,11:1,17:1"),
        ("1729", "7:1,13:1,19:1"),
        ("41041", "7:1,11:1,13:1,41:1"),
        ("825265", "5:1,7:1,17:1,19:1,73:1"),
        ("1000003", "1000003:1"),
        ("1000006000009", "1000003:2"),
        ("1000036000099", "1000003:1,1000033:1"),
        ("2305843009213693951", "2305843009213693951:1"),
        ("9223372036854775808", "2:63"),
        ("1000000000000000009", "1000000000000000009:1"),
        (
            "2432902008176640000",
            "2:18,3:8,5:4,7:2,11:1,13:1,17:1,19:1",
        ),
        ("1000000016000000063", "1000000007:1,1000000009:1"),
    ];
    assert_eq!(cases.len(), 22);
    for &(n, expected) in cases {
        let got = factorint_string(&int(n));
        // Convention difference: SymPy reports factorint(0) = {0: 1}, oCAS
        // returns an empty factor list.
        if n == "0" {
            assert_eq!(got, expected);
            continue;
        }
        check_nt("nt_factorint", n, got, expected);
    }
    // Property check on a 16-digit semiprime: factors are prime and
    // reconstruct the input.
    let n = int("1000003") * int("10000000019");
    let f = factor_integer(&n);
    let mut product = Integer::from(1);
    for (p, e) in &f {
        assert!(is_prime_bpsw(p));
        product *= &p.pow_u32(*e);
    }
    assert_eq!(product, n);
}

#[test]
fn isprime_vs_sympy() {
    let cases: &[(&str, bool)] = &[
        ("0", false),
        ("1", false),
        ("2", true),
        ("3", true),
        ("4", false),
        ("97", true),
        ("561", false),
        ("1105", false),
        ("1729", false),
        ("2047", false),
        ("3277", false),
        ("4033", false),
        ("7919", true),
        ("1000003", true),
        ("1000015", false),
        ("2147483647", true),
        ("2305843009213693951", true),
        ("1000000000000000009", true),
        ("1000000000000000021", false),
        ("18446744073709551557", true),
        ("18446744073709551615", false),
        ("825265", false),
        ("3215031751", false),
    ];
    assert_eq!(cases.len(), 23);
    for &(n, expected) in cases {
        let got = is_prime_bpsw(&int(n));
        assert_eq!(got, expected, "is_prime_bpsw({n})");
        if let Some(reference) = sympy_result("nt_isprime", n) {
            assert_eq!(
                got.to_string(),
                reference.to_lowercase(),
                "vs SymPy isprime({n})"
            );
        }
        // u64 wrapper agrees on in-range values.
        if let Ok(v) = n.parse::<u64>() {
            assert_eq!(is_prime_u64(v), expected, "is_prime_u64({n})");
        }
    }
}

#[test]
fn nextprime_vs_sympy() {
    let cases: &[(&str, &str)] = &[
        ("0", "2"),
        ("1", "2"),
        ("2", "3"),
        ("3", "5"),
        ("4", "5"),
        ("10", "11"),
        ("13", "17"),
        ("90", "97"),
        ("97", "101"),
        ("100", "101"),
        ("1000", "1009"),
        ("7919", "7927"),
        ("10000", "10007"),
        ("100000", "100003"),
        ("1000000", "1000003"),
        ("1000000000000", "1000000000039"),
        ("-5", "2"),
        ("-100", "2"),
        ("999983", "1000003"),
        ("2147483647", "2147483659"),
    ];
    assert_eq!(cases.len(), 20);
    for &(n, expected) in cases {
        let got = next_prime(&int(n)).to_string();
        check_nt("nt_nextprime", n, got, expected);
    }
}

#[test]
fn totient_vs_sympy() {
    let expected: &[(&str, &str)] = &[
        ("1", "1"),
        ("2", "1"),
        ("9", "6"),
        ("10", "4"),
        ("12", "4"),
        ("36", "12"),
        ("97", "96"),
        ("100", "40"),
        ("360", "96"),
        ("561", "320"),
        ("1000", "400"),
        ("1024", "512"),
        ("6561", "4374"),
        ("1000003", "1000002"),
        ("1000006000009", "1000005000006"),
        ("12345", "6576"),
        ("999999", "466560"),
        ("20", "8"),
        ("50", "20"),
        ("210", "48"),
    ];
    assert_eq!(expected.len(), 20);
    for &(n, want) in expected {
        let got = euler_phi(&int(n)).to_string();
        check_nt("nt_totient", n, got, want);
    }
}

#[test]
fn mobius_vs_sympy() {
    let cases: &[(&str, i8)] = &[
        ("1", 1),
        ("2", -1),
        ("3", -1),
        ("4", 0),
        ("6", 1),
        ("10", 1),
        ("12", 0),
        ("18", 0),
        ("30", -1),
        ("210", 1),
        ("360", 0),
        ("561", -1),
        ("1024", 0),
        ("1000003", -1),
        ("97", -1),
        ("2310", -1),
        ("36", 0),
        ("48", 0),
        ("999999", 0),
        ("1000006000009", 0),
    ];
    assert_eq!(cases.len(), 20);
    for &(n, want) in cases {
        let got = moebius_mu(&int(n));
        assert_eq!(got, want, "moebius_mu({n})");
        if let Some(reference) = sympy_result("nt_mobius", n) {
            assert_eq!(got.to_string(), reference, "vs SymPy mobius({n})");
        }
    }
}

#[test]
fn divisor_count_vs_sympy() {
    let cases: &[(&str, &str)] = &[
        ("1", "1"),
        ("2", "2"),
        ("6", "4"),
        ("12", "6"),
        ("36", "9"),
        ("97", "2"),
        ("100", "9"),
        ("360", "24"),
        ("1024", "11"),
        ("561", "8"),
        ("1000", "16"),
        ("6561", "9"),
        ("2187", "8"),
        ("20", "6"),
        ("210", "16"),
        ("12345", "8"),
        ("999999", "64"),
        ("1000003", "2"),
        ("72", "12"),
        ("144", "15"),
    ];
    assert_eq!(cases.len(), 20);
    for &(n, want) in cases {
        let got = divisor_tau(&int(n)).to_string();
        check_nt("nt_divisor_count", n, got, want);
    }
}

#[test]
fn divisor_sigma_vs_sympy() {
    let cases: &[(&str, &str)] = &[
        ("1;1", "1"),
        ("6;1", "12"),
        ("12;1", "28"),
        ("28;1", "56"),
        ("12;2", "210"),
        ("12;0", "6"),
        ("4;2", "21"),
        ("10;3", "1134"),
        ("36;1", "91"),
        ("97;1", "98"),
        ("100;2", "13671"),
        ("360;1", "1170"),
        ("6;0", "4"),
        ("20;2", "546"),
        ("210;1", "576"),
        ("99;1", "156"),
        ("1024;1", "2047"),
        ("6561;1", "9841"),
        ("45;2", "2366"),
        ("1000;1", "2340"),
    ];
    assert_eq!(cases.len(), 20);
    for &(nk, want) in cases {
        let (n, k) = nk.split_once(';').unwrap();
        let got = divisor_sigma(&int(n), k.parse().unwrap()).to_string();
        check_nt("nt_divisor_sigma", nk, got, want);
    }
}

#[test]
fn liouville_vs_sympy() {
    let cases: &[(&str, i8)] = &[
        ("1", 1),
        ("2", -1),
        ("3", -1),
        ("4", 1),
        ("6", 1),
        ("10", 1),
        ("12", -1),
        ("18", -1),
        ("30", -1),
        ("36", 1),
        ("48", -1),
        ("210", 1),
        ("360", 1),
        ("561", -1),
        ("97", -1),
        ("100", 1),
        ("1024", 1),
        ("2310", -1),
        ("999999", -1),
        ("1000003", -1),
    ];
    assert_eq!(cases.len(), 20);
    for &(n, want) in cases {
        let got = liouville_lambda(&int(n));
        assert_eq!(got, want, "liouville_lambda({n})");
        if let Some(reference) = sympy_result("nt_liouville", n) {
            assert_eq!(got.to_string(), reference, "vs SymPy liouville({n})");
        }
    }
}

#[test]
fn crt_vs_sympy() {
    let cases: &[(&str, &str, &str)] = &[
        ("3,5,7", "2,3,2", "23,105"),
        ("4,6", "1,3", "9,12"),
        ("2,3,5", "1,2,4", "29,30"),
        ("11,13", "7,11", "128,143"),
        ("97,101", "5,7", "4855,9797"),
        ("3,5,7", "2,3,5", "68,105"),
        ("17", "5", "5,17"),
        ("12,18", "4,10", "28,36"),
        ("1000003,1000033", "11,222", "233334700013,1000036000099"),
        ("8,12", "4,8", "20,24"),
        ("5,9", "3,5", "23,45"),
        ("7,11,13", "1,2,3", "211,1001"),
        ("31,37", "15,16", "201,1147"),
        ("2,9", "1,5", "5,18"),
        ("10,15", "5,10", "25,30"),
        ("19,23", "2,3", "325,437"),
        ("29,31", "17,19", "887,899"),
        ("41,43", "33,35", "1755,1763"),
        ("13,17,19", "12,16,18", "4198,4199"),
        ("6,9,15", "2,5,2", "32,90"),
    ];
    assert_eq!(cases.len(), 20);
    for &(ms, rs, want) in cases {
        let cs: Vec<(Integer, Integer)> = rs
            .split(',')
            .zip(ms.split(','))
            .map(|(r, m)| (int(r), int(m)))
            .collect();
        let (r, m) = crt_many(&cs).expect("consistent CRT system");
        let got = format!("{r},{m}");
        check_nt("nt_crt", &format!("{ms}|{rs}"), got, want);
        // Property: r satisfies every congruence.
        for (ri, mi) in &cs {
            assert_eq!(r.mod_floor(mi), ri.mod_floor(mi));
        }
    }
}

#[test]
fn jacobi_vs_sympy() {
    let cases: &[(&str, &str)] = &[
        ("2;7", "1"),
        ("3;7", "-1"),
        ("0;7", "0"),
        ("5;7", "-1"),
        ("6;7", "-1"),
        ("2;15", "1"),
        ("7;15", "-1"),
        ("4;15", "1"),
        ("10;15", "0"),
        ("1;9", "1"),
        ("2;9", "1"),
        ("5;9", "1"),
        ("7;9", "1"),
        ("8;9", "1"),
        ("3;11", "1"),
        ("5;11", "1"),
        ("2;11", "-1"),
        ("6;11", "-1"),
        ("7;11", "-1"),
        ("10;11", "-1"),
        ("9;11", "1"),
        ("22;55", "0"),
        ("12;35", "1"),
    ];
    assert_eq!(cases.len(), 23);
    for &(an, want) in cases {
        let (a, n) = an.split_once(';').unwrap();
        let got = jacobi(&int(a), &int(n)).to_string();
        check_nt("nt_jacobi", an, got, want);
    }
}

#[test]
fn discrete_log_vs_sympy() {
    // Primitive-root cases: the logarithm is unique modulo the order, so a
    // direct value comparison with SymPy is valid.
    let primitive_cases: &[(&str, &str, &str, &str)] = &[
        ("11", "2", "7", "7"),
        ("11", "2", "5", "4"),
        ("11", "2", "1", "0"),
        ("11", "2", "9", "6"),
        ("11", "2", "4", "2"),
        ("23", "5", "10", "3"),
        ("23", "5", "3", "16"),
        ("23", "5", "1", "0"),
        ("101", "2", "66", "83"),
        ("101", "2", "27", "7"),
        ("101", "2", "1", "0"),
        ("101", "2", "6", "70"),
    ];
    for &(p, b, t, want) in primitive_cases {
        let pb = int(p);
        let bb = int(b);
        let tb = int(t);
        let x = dlog_pohlig_hellman(&bb, &tb, &pb).expect("log must exist");
        assert_eq!(bb.modpow(&x, &pb), tb.mod_floor(&pb));
        check_nt(
            "nt_discrete_log",
            &format!("{p};{b};{t}"),
            x.to_string(),
            want,
        );
    }
    // Non-primitive / composite cases: verify the defining equation instead.
    let verify_cases: &[(&str, &str, &str)] = &[
        ("23", "2", "13"),  // ord(2) = 11
        ("23", "4", "9"),   // ord(4) = 11
        ("101", "4", "36"), // ord(4) = 50
        ("105", "2", "4"),  // composite modulus, BSGS
        ("105", "2", "16"),
        ("1009", "11", "354"),
        ("97", "5", "35"),
        ("1000", "3", "27"),
    ];
    for &(p, b, t) in verify_cases {
        let pb = int(p);
        let bb = int(b);
        let tb = int(t);
        let x = if is_prime_bpsw(&pb) {
            dlog_pohlig_hellman(&bb, &tb, &pb)
        } else {
            dlog_bsgs(&bb, &tb, &pb)
        }
        .expect("log must exist");
        assert_eq!(bb.modpow(&x, &pb), tb.mod_floor(&pb), "dlog({p},{b},{t})");
        // Cross-check with SymPy only when the group is cyclic of prime order
        // modulo the base order; otherwise solutions may differ by order
        // multiples.
        if let Some(reference) = sympy_result("nt_discrete_log", &format!("{p};{b};{t}")) {
            let x_ref = int(&reference);
            assert_eq!(bb.modpow(&x_ref, &pb), tb.mod_floor(&pb));
        }
    }
    // Unsolvable case must be rejected.
    assert!(dlog_pohlig_hellman(&int("2"), &int("5"), &int("23")).is_none());
}
