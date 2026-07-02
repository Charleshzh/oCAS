# Architecture

oCAS is organized as a Cargo workspace of 12 crates with strict downward
layering. Each crate may only depend on lower-level crates; no reverse or
cyclic dependencies are permitted.

| Level | crate | Responsibility |
|---|---|---|
| 5 Bindings | `ocas`, `ocas-py`, `ocas-c` | Rust, Python, C/C++ public API |
| 4 Application | `ocas-calc`, `ocas-eval`, `ocas-parse` | Calculus, evaluation, parsing |
| 3 Symbol engine | `ocas-atom`, `ocas-rewrite` | Atom, converters, pattern matching, e-graph |
| 2 Algebra kernel | `ocas-domain`, `ocas-poly` | Domains, polynomials, number theory |
| 1 Numerical backend | `ocas-core` | GMP/FLINT encapsulation |
| 0 Runtime | `ocas-core` | Arena, errors, thread pool, FFI |

See [ARCHITECTURE_EN.md](https://github.com/charleshzh/ocas/blob/main/ARCHITECTURE_EN.md)
in the repository for the full design document.
