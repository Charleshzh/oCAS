# oCAS — open Computer Algebra System

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.97%2B-orange.svg)](https://www.rust-lang.org)
[![CI](https://github.com/charleshzh/ocas/actions/workflows/ci.yml/badge.svg)](https://github.com/charleshzh/ocas/actions/workflows/ci.yml)
[![docs.rs](https://docs.rs/ocas/badge.svg)](https://docs.rs/ocas)
[![PyPI](https://img.shields.io/pypi/v/ocas.svg)](https://pypi.org/project/ocas/)

> **Languages / 语言**: [English](README.md) · [中文](docs/README_CN.md)

**oCAS** is a modern, high-performance computer algebra system written in Rust.
It is designed to match or exceed the core performance of Symbolica and
SageMath while remaining free and open under the **LGPL-3.0-or-later** license.

> **Status**: Beta (0.13.2). The full symbolic engine, polynomial algebra,
> equation solvers, JIT/SIMD evaluation, and tri-language bindings (Rust,
> Python, C/C++) are feature-complete. `pip install ocas` is available on PyPI.
> See the
> [documentation](https://charleshzh.github.io/ocas/latest/en/)
> ([中文](https://charleshzh.github.io/ocas/latest/zh/)),
> [Rust API docs](https://docs.rs/ocas), and
> [roadmap](docs/planning/ROADMAP_EN.md).
>
> **Note on docs.rs**: The hosted Rust API docs are built with portable
> features only (no system GMP/MPFR/FLINT backends). To browse the full API
> including backend features, build the documentation locally with
> `cargo doc -p ocas --features gmp,mpfr,flint --no-deps`.

---

## Why oCAS?

| Capability | oCAS | Symbolica | SageMath |
|---|---|---|---|
| Core language | Rust | Rust | Python/Cython |
| License | **LGPL-3.0+** | Proprietary / source-available | GPL-3.0 |
| Python API | ✅ PyO3 | ✅ | ✅ Native |
| C/C++ API | ✅ cbindgen | ✅ | ⚠️ Limited |
| Rust API | ✅ Native | ✅ | ❌ |
| High-performance backend | GMP / FLINT 3 / Arb / LinBox | Self-hosted + optional | FLINT / Singular / PARI |
| JIT code generation | Cranelift / LLVM | Custom | Limited |
| E-graph simplification | ✅ `egg` (optional) | ❌ | Partial |

oCAS targets researchers, engineers, and developers who need a fast,
embeddable, and license-clean symbolic computation engine.

---

## Key Features

- **Symbolic expressions**: variables, functions, rational numbers,
  arbitrary-precision integers, floats, and complex numbers.
- **Fast parsing & printing**: Mathematica-like and Python-like syntax, LaTeX
  output.
- **Simplification & rewriting**: normalization, hash consing, and equality
  saturation via `egg`.
- **Polynomial arithmetic**: dense and sparse multivariate polynomials, GCD,
  univariate and bivariate factorization (Hensel lifting), Gröbner bases,
  resultants.
- **Calculus**: symbolic differentiation, partial integration, Taylor/Laurent
  series.
- **Equation solving**: linear systems, polynomial systems, algebraic numbers.
- **Numerical evaluation**: arbitrary-precision floating point, interval
  arithmetic, numerical integration.
- **JIT compilation**: compile repeated expressions to native code with
  Cranelift or LLVM.
- **Multi-language bindings**: Rust, Python, and C/C++ from a single codebase.
- **Optimized integer arithmetic**: SOO (small-object optimization) stores
  values fitting in `i64` on the stack with GMP fallback for large values.
- **Modular multivariate GCD**: compute GCD over $\mathbb{Z}[x,y]$ via
  reduction mod prime and reconstruction.

---

## Architecture

```text
┌─────────────────────────────────────────────────────────────┐
│ Layer 5: Language Bindings                                  │
│  Python (PyO3) │ C/C++ (cbindgen) │ Rust API                │
├─────────────────────────────────────────────────────────────┤
│ Layer 4: Applications                                       │
│  Parser │ Printer │ Solver │ Calculus │ Series              │
├─────────────────────────────────────────────────────────────┤
│ Layer 3: Symbolic Engine                                    │
│  Expression (Atom) │ Transformer │ Pattern Matcher │ Rewriting │ E-Graph │
├─────────────────────────────────────────────────────────────┤
│ Layer 2: Algebraic Kernel                                   │
│  Polynomial │ Matrix │ Domain │ Number Theory │ Factorize  │
├─────────────────────────────────────────────────────────────┤
│ Layer 1: Numeric Backend                                    │
│  GMP/MPFR/MPC │ FLINT 3 / Arb │ NTL │ LinBox │ faer        │
├─────────────────────────────────────────────────────────────┤
│ Layer 0: Runtime                                            │
│  Arena Allocator │ Error Handling │ Thread Pool │ FFI Glue  │
└─────────────────────────────────────────────────────────────┘
```

See [ARCHITECTURE_EN.md](docs/planning/ARCHITECTURE_EN.md) for the full design.

---

## Quick Start

### Rust

```toml
[dependencies]
ocas = "0.13"
```

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

fn main() -> Result<(), ParseError> {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    let expr = parse(&ctx, "x^2 + 2*x + 1")?;
    println!("{}", expr); // 1 + 2*x + x^2

    Ok(())
}
```

Rule-based simplification, pattern matching, and optional `egg` equality
saturation are available through `ocas::prelude` (see `ocas-rewrite` examples).

---

## Building from Source

### Requirements

- Rust 1.97 or later
- C compiler (MSVC, GCC, or Clang)
- GMP, MPFR, and FLINT 3 development libraries

On Ubuntu/Debian:

```bash
sudo apt-get install libgmp-dev libmpfr-dev libflint-dev
```

On macOS:

```bash
brew install gmp mpfr flint
```

On Windows, see [docs/build-windows.md](docs/build-windows.md).

### Build

```bash
git clone https://github.com/charleshzh/ocas.git
cd ocas
cargo build --release
```

### Run Tests

```bash
cargo test --workspace --exclude ocas-py
```

---

## Performance

oCAS aims for competitive performance with Symbolica and SageMath on core
workloads:

- Integer and rational arithmetic via **GMP**.
- Polynomial multiplication, GCD, and factorization via **FLINT 3**.
- Exact linear algebra via **faer** (default) and **LinBox** (optional).
- JIT evaluation via **Cranelift**.

Benchmarks are under active development. See
[ocas-tests/benches](ocas-tests/benches). For the gap analysis against
Symbolica, SageMath, and SymPy, see [GAP_ANALYSIS_EN.md](docs/planning/GAP_ANALYSIS_EN.md).

---

## Roadmap

See [ROADMAP_EN.md](docs/planning/ROADMAP_EN.md) for the versioned plan, and
[EVOLUTION_PLAN_EN.md](docs/planning/EVOLUTION_PLAN_EN.md) for the fine-grained per-version
evolution plan (Beta → 1.0).

---

## License

This project is licensed under the **GNU Lesser General Public License v3.0 or
later** (LGPL-3.0-or-later). See [LICENSE](LICENSE) for the full text.

### Optional GPL Backends

Some optional backends (e.g., NTL, Singular, PARI/GP, SageMath interfaces) are
licensed under GPL-compatible licenses. These are isolated in the `ocas-gpl`
crate and disabled by default. Enabling them will produce a combined work
under GPL terms.

---

## Contributing

We welcome contributions. Please read [CONTRIBUTING.md](docs/book/en/src/contributing.md) before
opening a pull request.

---

## Acknowledgments

oCAS builds on the shoulders of many excellent open-source projects:

- [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/), [MPC](https://www.multiprecision.org/mpc/)
- [FLINT](https://flintlib.org/) and [Arb](https://arblib.org/) (now part of FLINT)
- [LinBox](https://linalg.org/) and [Givaro](https://github.com/linbox-team/givaro)
- [SymEngine](https://github.com/symengine/symengine)
- [faer](https://codeberg.org/sarah-quinones/faer) (primary development moved from [GitHub mirror](https://github.com/sarah-quinones/faer-rs))
- [egg](https://github.com/egraphs-good/egg)
- [Cranelift](https://github.com/bytecodealliance/wasmtime/tree/main/cranelift)
- [PyO3](https://github.com/PyO3/pyo3)
