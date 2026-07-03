# Architecture

oCAS is organized as a Cargo workspace of 12 public crates (plus
`ocas-tests` for integration tests and benchmarks) with strict downward
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

## Notable modules

**`ocas-poly::factor`** — Polynomial factorization using Hensel lifting and
finite-field techniques. Supports square-free decomposition and full
factorization over the integers and rationals.

**`ocas-poly::groebner`** — Gröbner basis computation (Buchberger algorithm)
with configurable monomial orders (Lex, Grevlex).

**`ocas-poly::roots`** — Real root isolation for univariate polynomials.

**`ocas-tests::correctness`** — Cross-validation framework with 82 tests
across 16 mathematical modules, comparing oCAS results against SymPy,
SageMath, and Symbolica.

See [ARCHITECTURE_EN.md](https://github.com/charleshzh/ocas/blob/main/docs/planning/ARCHITECTURE_EN.md)
in the repository for the full design document.
