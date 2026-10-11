//! Single-integral probe: parse, integrate, and print the result of one
//! integrand string. Used by the wave plans for per-case probing.
//!
//! ```text
//! cargo run --release -p ocas --example integrate_probe -- "1/(1+x^4)"
//! ```

use ocas::prelude::*;

fn main() {
    let src = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("usage: integrate_probe <integrand in x>");
        std::process::exit(2);
    });
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = match parse(&ctx, &src) {
        Ok(e) => e,
        Err(err) => {
            println!("PARSE_ERR: {err}");
            return;
        }
    };
    let result = integrate(&ctx, expr, Symbol::new("x"));
    println!("{result}");
    // With `outcome`, also print the three-valued outcome's certificate.
    if std::env::args().any(|a| a == "outcome") {
        match integrate_outcome(&ctx, expr, Symbol::new("x")) {
            Outcome::Found { certificate, .. } => {
                println!("cert: {}", certificate.method.as_str());
            }
            Outcome::ProvedNonElementary { .. } => println!("cert: non-elementary"),
            Outcome::Unknown { .. } => println!("cert: unknown"),
        }
    }
}
