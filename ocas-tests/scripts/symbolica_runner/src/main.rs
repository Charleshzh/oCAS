use std::env;
use std::process::ExitCode;

use symbolica::prelude::*;

fn run(task: &str, expr: &str) -> Result<String, String> {
    let x = symbol!("x");

    match task {
        "parse" => {
            let a = parse!(expr);
            Ok(format!("{}", a.expand()))
        }
        "diff" => {
            let a = parse!(expr);
            Ok(format!("{}", a.derivative(x)))
        }
        "expand" => {
            let a = parse!(expr);
            Ok(format!("{}", a.expand()))
        }
        "simplify" => {
            // Symbolica does not expose a single `simplify`; expand is a reasonable proxy.
            let a = parse!(expr);
            Ok(format!("{}", a.expand()))
        }
        "factor" => {
            let a = parse!(expr).expand();
            let poly: MultivariatePolynomial<_, u8> = a.to_polynomial(&Z, None);
            let factors = poly.factor();
            let mut parts = Vec::new();
            for (f, pow) in factors {
                if pow == 1 {
                    parts.push(format!("{}", f));
                } else {
                    parts.push(format!("({})^{}", f, pow));
                }
            }
            Ok(parts.join(" * "))
        }
        "series" => {
            let a = parse!(expr);
            let s = a
                .series(x, 0, 10)
                .map_err(|e| format!("series error: {:?}", e))?;
            Ok(format!("{}", s))
        }
        _ => Err(format!("unknown task: {}", task)),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: symbolica_runner <task> <expr>");
        return ExitCode::from(1);
    }

    match run(&args[1], &args[2]) {
        Ok(out) => {
            println!("{}", out);
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("{}", err);
            ExitCode::from(1)
        }
    }
}
