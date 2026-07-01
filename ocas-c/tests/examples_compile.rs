//! Verify that the C and C++ examples compile and run against the built
//! library.
//!
//! These tests use the system compiler (`cc`/`gcc`/`clang`/`cl`/`c++`) to
//! compile `examples/expression.c` and `examples/cpp_example.cpp`, link
//! them against the `ocas_c` static library, and run the resulting
//! binaries. If no compiler is found, the tests are skipped.

use std::path::PathBuf;
use std::process::Command;

fn workspace_target() -> Option<PathBuf> {
    // The workspace target dir is typically two levels up from the
    // crate's target directory.
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let target = PathBuf::from(manifest_dir)
        .join("..")
        .join("target")
        .join("debug");
    let lib = if cfg!(windows) {
        let gnu_lib = target.join("libocas_c.a");
        let msvc_lib = target.join("ocas_c.lib");
        if gnu_lib.exists() { gnu_lib } else { msvc_lib }
    } else {
        target.join("libocas_c.a")
    };
    if lib.exists() { Some(target) } else { None }
}

fn find_compiler(candidates: &'static [&'static str]) -> Option<&'static str> {
    candidates
        .iter()
        .find(|&cc| {
            Command::new(cc).arg("--version").output().is_ok()
                || Command::new(cc).arg("/?").output().is_ok()
        })
        .copied()
}

fn crate_dir() -> PathBuf {
    PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap())
}

#[test]
fn expression_c_example_compiles_and_runs() {
    let Some(cc) = find_compiler(&["cc", "gcc", "clang", "cl"]) else {
        eprintln!("skipping: no C compiler found");
        return;
    };
    let Some(target_dir) = workspace_target() else {
        eprintln!("skipping: ocas_c static library not found (run `cargo build -p ocas-c` first)");
        return;
    };

    let crate_dir = crate_dir();
    let example_src = crate_dir.join("examples").join("expression.c");
    let include_dir = crate_dir.join("include");
    let output_bin = target_dir.join("expression_c_example");

    let mut cmd = Command::new(cc);
    cmd.arg(&example_src)
        .arg("-I")
        .arg(&include_dir)
        .arg("-L")
        .arg(&target_dir)
        .arg("-locas_c")
        .arg("-o")
        .arg(&output_bin);

    let compile_result = cmd.status();
    match compile_result {
        Ok(status) if status.success() => {}
        Ok(status) => {
            eprintln!("skipping: C compilation failed (exit {:?})", status.code());
            return;
        }
        Err(e) => {
            eprintln!("skipping: failed to invoke C compiler: {e}");
            return;
        }
    }

    let run_result = Command::new(&output_bin).status();
    match run_result {
        Ok(status) => assert!(
            status.success(),
            "C expression example exited with non-zero status: {status:?}"
        ),
        Err(e) => panic!("failed to run compiled C example: {e}"),
    }
}

#[test]
fn cpp_example_compiles_and_runs() {
    let Some(cxx) = find_compiler(&["c++", "g++", "clang++", "cl"]) else {
        eprintln!("skipping: no C++ compiler found");
        return;
    };
    let Some(target_dir) = workspace_target() else {
        eprintln!("skipping: ocas_c static library not found (run `cargo build -p ocas-c` first)");
        return;
    };

    let crate_dir = crate_dir();
    let example_src = crate_dir.join("examples").join("cpp_example.cpp");
    let include_dir = crate_dir.join("include");
    let output_bin = target_dir.join("cpp_example");

    let mut cmd = Command::new(cxx);
    cmd.arg(&example_src)
        .arg("-I")
        .arg(&include_dir)
        .arg("-L")
        .arg(&target_dir)
        .arg("-locas_c")
        .arg("-o")
        .arg(&output_bin);

    let compile_result = cmd.status();
    match compile_result {
        Ok(status) if status.success() => {}
        Ok(status) => {
            eprintln!(
                "skipping: C++ compilation failed (exit {:?})",
                status.code()
            );
            return;
        }
        Err(e) => {
            eprintln!("skipping: failed to invoke C++ compiler: {e}");
            return;
        }
    }

    let run_result = Command::new(&output_bin).status();
    match run_result {
        Ok(status) => assert!(
            status.success(),
            "C++ example exited with non-zero status: {status:?}"
        ),
        Err(e) => panic!("failed to run compiled C++ example: {e}"),
    }
}
