use std::env;
use std::hint::black_box;
use std::process::ExitCode;
use std::time::Instant;

use symbolica::prelude::*;

/// Build a closure that executes `task` on `expr` once per call.
///
/// Correctness mode calls it once and prints the result; `time` mode calls
/// it once as a warmup and then `iters` times, printing the total
/// nanoseconds — the same contract as `compare_sympy.py time` so the two
/// sides can be divided directly.
fn build_op<'a>(task: &'a str, expr: &'a str) -> Result<Box<dyn Fn() -> String + 'a>, String> {
    let x = symbol!("x");

    match task {
        "parse" => {
            // Re-parse on every call (mirrors the SymPy `parse` task).
            Ok(Box::new(move || {
                let a = parse!(expr);
                format!("{}", a.expand())
            }))
        }
        "diff" => {
            let a = parse!(expr);
            Ok(Box::new(move || format!("{}", a.derivative(x))))
        }
        "expand" => {
            let a = parse!(expr);
            Ok(Box::new(move || format!("{}", a.expand())))
        }
        "simplify" => {
            // Symbolica does not expose a single `simplify`; expand is a reasonable proxy.
            let a = parse!(expr);
            Ok(Box::new(move || format!("{}", a.expand())))
        }
        "factor" => {
            let a = parse!(expr).expand();
            let poly: MultivariatePolynomial<_, u8> = a.to_polynomial(&Z, None);
            Ok(Box::new(move || {
                let factors = poly.factor();
                let mut parts = Vec::new();
                for (f, pow) in factors {
                    if pow == 1 {
                        parts.push(format!("{}", f));
                    } else {
                        parts.push(format!("({})^{}", f, pow));
                    }
                }
                parts.join(" * ")
            }))
        }
        "series" => {
            // Validate once so errors surface at build time; recompute per call.
            let e = expr.to_string();
            let probe = parse!(e.as_str());
            probe
                .series(x, 0, 10)
                .map_err(|err| format!("series error: {:?}", err))?;
            let x0 = symbol!("x");
            Ok(Box::new(move || {
                let b = parse!(e.as_str());
                let s = b.series(x0, 0, 10).expect("series validated");
                format!("{}", s)
            }))
        }
        _ => Err(format!("unknown task: {}", task)),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();

    if args.len() >= 2 && args[1] == "time" {
        if args.len() != 5 {
            eprintln!("Usage: symbolica_runner time <task> <expr> <iters>");
            return ExitCode::from(1);
        }
        let iters: u32 = match args[4].parse() {
            Ok(n) => n,
            Err(_) => {
                eprintln!("invalid iters: {}", args[4]);
                return ExitCode::from(1);
            }
        };
        return match build_op(&args[2], &args[3]) {
            Ok(op) => {
                op(); // warmup
                let start = Instant::now();
                for _ in 0..iters {
                    black_box(op());
                }
                println!("{}", start.elapsed().as_nanos());
                ExitCode::SUCCESS
            }
            Err(err) => {
                eprintln!("{}", err);
                ExitCode::from(1)
            }
        };
    }

    if args.len() != 3 {
        eprintln!("Usage: symbolica_runner <task> <expr>");
        return ExitCode::from(1);
    }

    match build_op(&args[1], &args[2]) {
        Ok(op) => {
            println!("{}", op());
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("{}", err);
            ExitCode::from(1)
        }
    }
}
