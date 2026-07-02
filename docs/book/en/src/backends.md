# Numerical Backends

oCAS ships a pure-Rust default build (using `num-bigint`/`num-rational`) and
optional high-performance backends behind feature flags. The default build
uses only LGPL-compatible dependencies; GPL-exclusive code is isolated in the
optional `ocas-gpl` crate.

## Feature flags

| Feature | Crate | Backend | License |
|---|---|---|---|
| `gmp` | `ocas-core`, `ocas-domain` | GMP integers/rationals via `rug` | LGPL |
| `mpfr` | `ocas-domain` | MPFR real balls via `rug` | LGPL |
| `flint` | `ocas-poly` | FLINT 3 (Linux/WSL only) | LGPL |
| `jit` | `ocas-eval` | Cranelift JIT compilation | Apache-2.0 |
| `simd` | `ocas-eval` | SIMD-vectorized evaluation | MIT |
| `egg` | `ocas-rewrite` | E-graph equality saturation | MIT |

## Enabling backends

```bash
# GMP + MPFR backends (recommended on Windows via MSYS2 MINGW64)
cargo build -p ocas --features gmp,mpfr

# FLINT backend (Linux/WSL only; requires system libflint)
cargo build -p ocas --features gmp,mpfr,flint

# JIT + SIMD evaluation
cargo build -p ocas --features jit,simd
```

See [docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)
for Windows-specific setup of GMP/MPFR via MSYS2.

## GPL backends

The `ocas-gpl` crate (licensed GPL-3.0-or-later) is opt-in and never included
in the default build. It hosts GPL-exclusive algorithms that cannot ship under
LGPL. Default oCAS builds remain LGPL-clean.
