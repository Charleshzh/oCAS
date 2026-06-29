# oCAS — open Computer Algebra System

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.89%2B-orange.svg)](https://www.rust-lang.org)

**oCAS** is a modern, high-performance computer algebra system written in Rust. It is designed to match or exceed the core performance of Symbolica and SageMath while remaining free and open under the **LGPL-3.0-or-later** license.

> **Status**: Early development. APIs and crate boundaries are subject to change.

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
| E-graph simplification | ✅ `egg` | ❌ | Partial |

oCAS targets researchers, engineers, and developers who need a fast, embeddable, and license-clean symbolic computation engine.

## Key Features

- **Symbolic expressions**: variables, functions, rational numbers, arbitrary-precision integers, floats, and complex numbers.
- **Fast parsing & printing**: Mathematica-like and Python-like syntax, LaTeX output.
- **Simplification & rewriting**: normalization, hash consing, and equality saturation via `egg`.
- **Polynomial arithmetic**: dense and sparse multivariate polynomials, GCD, factorization, Gröbner bases, resultants.
- **Calculus**: symbolic differentiation, partial integration, Taylor/Laurent series.
- **Equation solving**: linear systems, polynomial systems, algebraic numbers.
- **Numerical evaluation**: arbitrary-precision floating point, interval arithmetic, numerical integration.
- **JIT compilation**: compile repeated expressions to native code with Cranelift or LLVM.
- **Multi-language bindings**: Rust, Python, and C/C++ from a single codebase.

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
│  Expression (Atom) │ Transformer │ Pattern Matcher │ E-Graph │
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

See [ARCHITECTURE.md](ARCHITECTURE.md) for the full design.

## Quick Start

### Rust

```toml
[dependencies]
ocas = "0.1"
```

```rust
use ocas::prelude::*;

fn main() -> Result<(), OcasError> {
    let x = symbol!("x");
    let y = symbol!("y");

    let expr = parse!("x^2 + 2*x*y + y^2")?;
    let factored = expr.factor()?;
    println!("{}", factored); // (x + y)^2

    let deriv = expr.derivative(&x)?;
    println!("d/dx = {}", deriv); // 2*x + 2*y

    Ok(())
}
```

### Python

```python
import ocas

x, y = ocas.symbols("x y")
expr = ocas.parse("x^2 + 2*x*y + y^2")
print(expr.factor())        # (x + y)^2
print(expr.diff(x))         # 2*x + 2*y
```

## Building from Source

### Requirements

- Rust 1.89 or later
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

On Windows, see [docs/build-windows.md](docs/build-windows.md) (TODO).

### Build

```bash
git clone https://github.com/ocas-dev/ocas.git
cd ocas
cargo build --release
```

### Run Tests

```bash
cargo test --workspace
```

## Performance

oCAS aims for competitive performance with Symbolica and SageMath on core workloads:

- Integer and rational arithmetic via **GMP**.
- Polynomial multiplication, GCD, and factorization via **FLINT 3**.
- Exact linear algebra via **faer** (default) and **LinBox** (optional).
- JIT evaluation via **Cranelift**.

Benchmarks are under active development. See [ocas-tests/benches](ocas-tests/benches).

## Roadmap

See [ROADMAP.md](ROADMAP.md) for the phased development plan.

High-level timeline:

1. **Foundation** — arena, error handling, GMP/FLINT bindings
2. **Expression core** — AST, parser, printer, normalization
3. **Domain & polynomials** — integers, rationals, finite fields, polynomial algorithms
4. **Calculus & solvers** — differentiation, integration, equation solving
5. **Evaluation & JIT** — AST compiler, Cranelift backend
6. **Bindings & ecosystem** — Python, C/C++, optional GPL backends

## License

This project is licensed under the **GNU Lesser General Public License v3.0 or later** (LGPL-3.0-or-later).

See [LICENSE](LICENSE) for the full text.

### Optional GPL Backends

Some optional backends (e.g., NTL, Singular, PARI/GP, SageMath interfaces) are licensed under GPL-compatible licenses. These are isolated in the `ocas-gpl` crate and disabled by default. Enabling them will produce a combined work under GPL terms.

## Contributing

We welcome contributions. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

## Acknowledgments

oCAS builds on the shoulders of many excellent open-source projects:

- [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/), [MPC](https://www.multiprecision.org/mpc/)
- [FLINT](https://flintlib.org/) and [Arb](https://arblib.org/)
- [LinBox](https://linalg.org/) and [Givaro](https://github.com/linbox-team/givaro)
- [SymEngine](https://github.com/symengine/symengine)
- [faer](https://github.com/sarah-quinones/faer-rs)
- [egg](https://github.com/egraphs-good/egg)
- [Cranelift](https://github.com/bytecodealliance/wasmtime/tree/main/cranelift)
- [PyO3](https://github.com/PyO3/pyo3)
