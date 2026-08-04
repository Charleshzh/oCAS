//! Single-iteration Gröbner timing for large cyclic-n ideals.
//!
//! Criterion requires ≥10 samples, which is impractical when one
//! iteration takes minutes. These `#[ignore]`d tests time a single run
//! each and are run manually in release mode:
//!
//! ```text
//! cargo test -p ocas-tests --release --test groebner_timing -- --ignored --nocapture
//! ```

use std::time::Instant;

use ocas_poly::groebner::f4::f4;
use ocas_poly::groebner::f5::f5;
use ocas_poly::groebner::multi_modular::groebner_basis_multi_modular;
use ocas_tests::systems::*;

fn time_f4(n: usize, p: u32) {
    let ideal = cyclic_fp(n, p);
    let start = Instant::now();
    let gb = f4(&ideal);
    let elapsed = start.elapsed();
    println!(
        "f4 cyclic_{n} over Z_{p}: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

fn time_f5(n: usize, p: u32) {
    let ideal = cyclic_fp(n, p);
    let start = Instant::now();
    let gb = f5(&ideal);
    let elapsed = start.elapsed();
    println!(
        "f5 cyclic_{n} over Z_{p}: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

fn time_f5_grevlex(n: usize, p: u32) {
    let ideal = cyclic_fp_grevlex(n, p);
    let start = Instant::now();
    let gb = f5(&ideal);
    let elapsed = start.elapsed();
    println!(
        "f5 grevlex cyclic_{n} over Z_{p}: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

fn time_multi_modular_q(n: usize) {
    let ideal = cyclic_q(n);
    let start = Instant::now();
    let gb = groebner_basis_multi_modular(&ideal);
    let elapsed = start.elapsed();
    println!(
        "multi_modular cyclic_{n} over Q: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_5_fp13() {
    time_f4(5, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_6_fp13() {
    time_f4(6, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn cyclic_7_fp13() {
    time_f4(7, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn f5_cyclic_6_fp13() {
    time_f5(6, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn f5_cyclic_7_fp13() {
    time_f5(7, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn f5_cyclic_7_fp13_grevlex() {
    time_f5_grevlex(7, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn multi_modular_cyclic_6_q() {
    time_multi_modular_q(6);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn multi_modular_cyclic_7_q() {
    time_multi_modular_q(7);
}

fn time_f5_katsura(n: usize, p: u32) {
    let ideal = katsura_fp(n, p);
    let start = Instant::now();
    let gb = f5(&ideal);
    let elapsed = start.elapsed();
    println!(
        "f5 katsura_{n} over Z_{p}: {:.3} s ({} basis elements)",
        elapsed.as_secs_f64(),
        gb.basis.len()
    );
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn f5_katsura_6_fp13() {
    time_f5_katsura(6, 13);
}

#[test]
#[ignore = "single-iteration timing; run manually in release mode"]
fn f5_katsura_7_fp13() {
    time_f5_katsura(7, 13);
}
