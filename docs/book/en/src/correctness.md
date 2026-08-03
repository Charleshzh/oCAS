# Correctness

oCAS includes an automated correctness cross-validation framework that
compares results against three reference systems: SymPy, SageMath, and
Symbolica. This chapter describes the framework, how to run it, and its
current known limitations.

---

## Framework overview

The correctness suite lives in `ocas-tests/tests/correctness/` and contains
**201 `#[test]` tests across 19 mathematical modules** (as of 0.24.0; 35 of
them are marked `#[ignore]` to track known gaps). Each test:

1. Generates an input (expression, polynomial, system of equations, etc.)
2. Computes a result using oCAS
3. Computes the equivalent result using a reference system
4. Asserts that both results are semantically equal

The modules cover the full breadth of oCAS functionality:

| Module | Tests | Covers |
|---|---|---|
| `calculus` | 16 | Differentiation, Taylor, integration (SymPy-checked) |
| `evaluation` | 6 | Numeric evaluation |
| `finite_field` | 5 | Finite field arithmetic |
| `groebner` | 18 | Gröbner basis computation |
| `integral_risch` | 15 | Risch symbolic integration |
| `linear_solve` | 5 | Linear solvers |
| `matrix` | 5 | Matrix operations |
| `normalize` | 8 | Expression normalization |
| `ntheory` | 11 | Number theory (cross-checked against SymPy `ntheory`) |
| `ode` | 34 | ODE solving |
| `parse` | 6 | Expression parsing and printing |
| `partial_fraction` | 8 | Partial fraction decomposition |
| `poly_arithmetic` | 6 | Dense/sparse polynomial arithmetic |
| `poly_factor` | 12 | Square-free and full factorization |
| `poly_factor_anf` | 21 | Algebraic number field factorization |
| `poly_gcd` | 5 | Polynomial GCD |
| `resultant` | 8 | Resultants |
| `rewrite` | 8 | Rewriting and simplification |
| `root_isolation` | 4 | Real root isolation |

---

## Ignored tests (known gaps)

Some tests are marked `#[ignore]` (35 in total) to track known gaps — they
are expected to fail and are only run manually when reproducing or
advancing the issue (`cargo test -p ocas-tests --test correctness
-- --ignored`). For example, Wilkinson polynomial root-finding: only 8 of
the 10 real roots are found.

---

## Running the tests

```bash
# Run all correctness tests
cargo test -p ocas-tests --test correctness

# Run a specific module (e.g. ODE solving)
cargo test -p ocas-tests --test correctness ode

# Run the ignored tests (known gaps)
cargo test -p ocas-tests --test correctness -- --ignored

# Run with verbose output to inspect failures
cargo test -p ocas-tests --test correctness -- --nocapture
```

The tests require no external dependencies — all reference computations
use SymPy through a `uv`-managed Python subprocess that provisions itself
automatically.

---

## Comparison harnesses

Separate scripts provide manual cross-checking against SageMath and Symbolica
for deeper investigation:

```bash
# SageMath (requires `sage` installed locally)
cd ocas-tests
sage scripts/bench_sage.py factor "x^30 - 1" 100

# Symbolica (requires Symbolica source checkout)
cd ../symbolica
cargo run --release --example factorization
```

These harnesses were used to validate the automated test suite during
development and are maintained for manual regression testing.

---

## Audit report

`ocas-tests/scripts/generate_audit_report.py` runs the full suite (including
the `--ignored` tests) and generates an audit report:

```bash
cd ocas-tests
python scripts/generate_audit_report.py
```

The report is written to `docs/planning/correctness/audit-<date>.md` and
includes:
- Pass/fail/ignored count summaries for the regular and ignored test runs
- A list of failing tests
- A factorization-timing comparison against Symbolica

---

## Known limitations

| Issue | Module | Status |
|---|---|---|
| Wilkinson n=10: only 8 of 10 real roots found | `root_isolation` | Known gap (`#[ignore]`) |
| `sin(x)^2 + cos(x)^2 → 1` requires `egg` feature | `rewrite` | Simplifies with the `egg` feature |
| Bernoulli forcing y^n confuses linear coefficient extraction | `ode` | Known limitation (`#[ignore]`) |
| Integrator missing tan/sec table entries | `ode` | Planned (`#[ignore]`) |

---

## See also

- [Performance](./performance.md) — benchmark suite details
- [Rust API](./api/rust.md) — core types used in tests
- [Contributing](./contributing.md) — how to add new correctness tests
