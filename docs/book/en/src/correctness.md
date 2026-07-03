# Correctness

oCAS includes an automated correctness cross-validation framework that
compares results against three reference systems: SymPy, SageMath, and
Symbolica. This chapter describes the framework, how to run it, and its
current known limitations.

---

## Framework overview

The correctness suite lives in `ocas-tests/tests/correctness/` and contains
**82 tests across 16 mathematical modules**. Each test:

1. Generates an input (expression, polynomial, system of equations, etc.)
2. Computes a result using oCAS
3. Computes the equivalent result using a reference system
4. Asserts that both results are semantically equal

The modules cover the full breadth of oCAS functionality:

| Module | Tests | Covers |
|---|---|---|
| `algebra` | 8 | Simplification, identity laws |
| `calculus_diff` | 7 | Symbolic differentiation |
| `calculus_int` | 5 | Heuristic integration |
| `calculus_series` | 4 | Taylor series expansion |
| `domain_integer` | 5 | Integer arithmetic, GCD |
| `domain_rational` | 4 | Rational arithmetic |
| `domain_finite_field` | 3 | Finite field arithmetic |
| `evaluation` | 6 | Numeric evaluation |
| `groebner` | 4 | Gröbner basis computation |
| `linear_algebra` | 6 | Matrix operations, linear solving |
| `parsing` | 5 | Expression parsing and printing |
| `poly_dense` | 5 | Dense polynomial arithmetic |
| `poly_factor` | 4 | Square-free and full factorization |
| `poly_gcd` | 5 | Polynomial GCD |
| `poly_sparse` | 5 | Sparse multivariate arithmetic |
| `solvers` | 6 | Linear and Diophantine solvers |

---

## Difficulty tiers

Tests are classified by difficulty to help target debugging effort:

| Tier | Description | Count |
|---|---|---|
| Trivial | Basic sanity checks (e.g. `x + 0 = x`) | ~20 |
| Easy | Single-step operations (e.g. `d/dx x^3`) | ~30 |
| Medium | Multi-step or moderate complexity | ~20 |
| Hard | Large expressions, edge cases | ~8 |
| Extreme | Known to exercise limitations | ~4 |

Extreme-tier tests are expected to **fail** and track known gaps
(e.g., Wilkinson polynomial root-finding: 8 of 10 real roots found).

---

## Running the tests

```bash
# Run all correctness tests
cargo test -p ocas-tests --test correctness

# Run a specific module
cargo test -p ocas-tests --test correctness algebra

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

Running the full suite generates `correctness_report.md` with a per-module
summary of pass/fail counts and annotations for known limitations.

```bash
cargo test -p ocas-tests --test correctness -- --generate-report
```

The report includes:
- Pass/fail counts per module
- List of failing tests with expected vs actual results
- Difficulty-tier breakdown
- Annotated known gaps with links to tracking issues

---

## Known limitations

| Issue | Module | Status |
|---|---|---|
| Wilkinson n=10: only 8 of 10 real roots found | `poly_sparse` / roots | Under investigation |
| `sin(x)^2 + cos(x)^2 → 1` requires `egg` feature | `algebra` | Works with `egg` feature |
| Some heuristic integrals not yet implemented | `calculus_int` | Expanding in 0.12 |
| Polynomial factorization limited to Z[x] | `poly_factor` | Q[x] planned for 0.12 |

---

## See also

- [Performance](./performance.md) — benchmark suite details
- [Rust API](./rust-api.md) — core types used in tests
- [Contributing](./contributing.md) — how to add new correctness tests
