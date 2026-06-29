# oCAS Roadmap

This roadmap outlines the phased development of oCAS from a minimal Rust CAS core to a multi-language, high-performance computer algebra system.

Dates are approximate and will shift based on contributor availability.

---

## Phase 0: Foundation (Weeks 1–4)

**Goal**: establish the workspace, runtime layer, and numerical backend bindings.

### Deliverables

- [ ] Workspace structure with `ocas-core` crate
- [ ] CI pipeline: `cargo test`, `cargo clippy`, `cargo-deny`, formatting
- [ ] Unified error type `OcasError`
- [ ] Arena / bump allocator for expression nodes
- [ ] Thread pool wrapper around `rayon`
- [ ] FFI glue types and conventions
- [ ] GMP bindings via `rug` or thin `bindgen` wrappers
- [ ] FLINT 3 `bindgen` setup and minimal polynomial wrapper

### Success Criteria

- `cargo build` succeeds.
- Arena allocates and drops a small tree without leaks ( Miri / valgrind clean).
- GMP integer arithmetic is callable from Rust.

---

## Phase 1: Expression Core (Weeks 5–10)

**Goal**: build the symbolic expression tree, parser, printer, and normalizer.

### Deliverables

- [ ] `ocas-atom` crate
- [ ] `Atom` tagged-union design: integer, rational, float, symbol, function, add, mul, pow
- [ ] Arena-backed AST with lifetime-safe public API
- [ ] Hash consing for common subexpressions
- [ ] Lexer using `logos`
- [ ] Recursive-descent / Pratt parser with Mathematica-like and Python-like syntax
- [ ] Printer: ASCII, LaTeX, and compact forms
- [ ] Normalizer: flatten `Add`/`Mul`, sort terms, merge numeric coefficients
- [ ] Property tests for parsing round-trips and algebraic identities

### Success Criteria

- `parse("x^2 + 2*x + 1")` produces the expected AST.
- `to_string(parse(s)) == s` for a broad set of expressions.
- Normalization is deterministic and tested.

---

## Phase 2: Domain and Polynomial System (Weeks 11–20)

**Goal**: implement algebraic domains and polynomial algorithms.

### Deliverables

- [ ] `ocas-domain` crate
- [ ] Domains: `Integer`, `Rational`, `FiniteField`, `RealBall`, `Complex`
- [ ] Domain trait: `add`, `mul`, `inv`, `eq`, `characteristic`
- [ ] `ocas-poly` crate
- [ ] Dense univariate polynomial
- [ ] Sparse multivariate polynomial
- [ ] Polynomial addition, multiplication, division with remainder
- [ ] GCD: Euclidean and subresultant PRS
- [ ] Factorization: univariate via FLINT, multivariate via Wang/EEZ fallback
- [ ] Gröbner basis: Buchberger and F4 (basic implementation)
- [ ] Taylor/Laurent series expansion
- [ ] Matrix type backed by `faer` and optional `LinBox`

### Success Criteria

- Polynomial operations match SymPy/SageMath outputs on regression suite.
- Factorization works for standard benchmark polynomials.
- Gröbner basis computes a simple system correctly.

---

## Phase 3: Calculus and Solvers (Weeks 21–32)

**Goal**: implement differentiation, integration, and equation solving.

### Deliverables

- [ ] `ocas-calc` crate
- [ ] Symbolic differentiation using recursive rules + pattern matching
- [ ] Derivative table for elementary functions
- [ ] Partial integration with Risch-like heuristics
- [ ] Integration table + pattern matching
- [ ] Taylor series expansion
- [ ] Linear system solver (`faer` / `LinBox`)
- [ ] Polynomial system solver (Gröbner + root isolation)
- [ ] Single-variable polynomial root finding (Arb)
- [ ] Diophantine and algebraic number basics

### Success Criteria

- Differentiation passes a comprehensive test suite.
- Integration succeeds on standard calculus problems.
- Linear and polynomial solvers produce correct results.

---

## Phase 4: Evaluation and JIT (Weeks 33–44)

**Goal**: build fast numerical evaluators and code generation.

### Deliverables

- [ ] `ocas-eval` crate
- [ ] Tree interpreter for scalar and vector evaluation
- [ ] AST-to-instruction compiler with constant folding
- [ ] Function registry for user-defined scalar/vector functions
- [ ] Cranelift JIT backend for scalar expressions
- [ ] Optional LLVM/Inkwell backend for AOT
- [ ] SIMD vectorized evaluation for dense arrays
- [ ] Error-propagating float evaluation
- [ ] Benchmark suite comparing against Symbolica and SageMath

### Success Criteria

- JIT evaluates repeated expressions faster than the interpreter by an order of magnitude.
- Benchmarks demonstrate parity or better on polynomial and linear-algebra tasks.

---

## Phase 5: Language Bindings (Weeks 45–52)

**Goal**: expose oCAS to Python and C/C++ users.

### Deliverables

- [ ] `ocas-py` crate with PyO3
- [ ] Python classes: `Expression`, `Polynomial`, `Matrix`, `Domain`
- [ ] Python wheel build via Maturin
- [ ] `ocas-c` crate with cbindgen
- [ ] Stable C API for expression lifecycle, parsing, evaluation
- [ ] C++ RAII wrapper
- [ ] Memory management verified across all bindings with valgrind / ASan
- [ ] Documentation and examples for each language

### Success Criteria

- Python `pip install ocas` works on Linux/macOS/Windows.
- C example compiles and runs against the shared library.
- No memory leaks detected in binding tests.

---

## Phase 6: GPL Backends and Ecosystem (Ongoing)

**Goal**: integrate optional GPL backends and grow the ecosystem.

### Deliverables

- [ ] `ocas-gpl` crate isolating GPL-only dependencies
- [ ] NTL wrapper for finite-field factorization
- [ ] Singular interface for advanced Gröbner bases and primary decomposition
- [ ] PARI/GP interface for algebraic number theory
- [ ] SageMath interoperability layer
- [ ] Symbolic regression / optimization examples
- [ ] SciML-style code generation for differential equations
- [ ] Plugin API for user-defined domains and rewrite rules

### Success Criteria

- `ocas-gpl` builds only when explicitly enabled.
- Optional backends are documented with license implications.
- Community contributes at least one non-trivial plugin or backend.

---

## Milestones

| Milestone | Target Date | Key Deliverable |
|---|---|---|
| M0 | Month 1 | Workspace builds, GMP/FLINT callable |
| M1 | Month 3 | Parse/print/normalize works end-to-end |
| M2 | Month 6 | Polynomial system competitive with SymPy |
| M3 | Month 9 | Calculus and solvers functional |
| M4 | Month 12 | JIT + Python/C bindings released |
| M5 | Month 18+ | Optional GPL backends and ecosystem |

---

## How to Contribute to the Roadmap

- Open an issue to propose a new feature or milestone adjustment.
- Pick an unclaimed deliverable from the current phase and comment on the tracking issue.
- See [CONTRIBUTING.md](CONTRIBUTING.md) for coding conventions and PR process.
