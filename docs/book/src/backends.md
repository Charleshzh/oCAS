# Numerical Backends / 数值后端

**English**

oCAS ships a pure-Rust default build (using `num-bigint`/`num-rational`) and
optional high-performance backends behind feature flags. The default build
uses only LGPL-compatible dependencies; GPL-exclusive code is isolated in the
optional `ocas-gpl` crate.

**中文**

oCAS 默认提供纯 Rust 构建（使用 `num-bigint`/`num-rational`），并通过 feature flag 提供可选的高性能后端。默认构建仅使用 LGPL 兼容依赖；GPL 专属代码隔离在可选的 `ocas-gpl` crate 中。

## Feature flags / 功能开关

| Feature | Crate | Backend | License |
|---|---|---|---|
| `gmp` | `ocas-core`, `ocas-domain` | GMP integers/rationals via `rug` | LGPL |
| `mpfr` | `ocas-domain` | MPFR real balls via `rug` | LGPL |
| `flint` | `ocas-poly` | FLINT 3 (Linux/WSL only) | LGPL |
| `jit` | `ocas-eval` | Cranelift JIT compilation | Apache-2.0 |
| `simd` | `ocas-eval` | SIMD-vectorized evaluation | MIT |
| `egg` | `ocas-rewrite` | E-graph equality saturation | MIT |

## Enabling backends / 启用后端

```bash
# GMP + MPFR backends (recommended on Windows via MSYS2 MINGW64)
cargo build -p ocas --features gmp,mpfr

# FLINT backend (Linux/WSL only; requires system libflint)
cargo build -p ocas --features gmp,mpfr,flint

# JIT + SIMD evaluation
cargo build -p ocas --features jit,simd
```

**English**

See [docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)
for Windows-specific setup of GMP/MPFR via MSYS2.

**中文**

Windows 上通过 MSYS2 安装 GMP/MPFR 的步骤见
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)。

## GPL backends / GPL 后端

**English**

The `ocas-gpl` crate (licensed GPL-3.0-or-later) is opt-in and never included
in the default build. It hosts GPL-exclusive algorithms that cannot ship under
LGPL. Default oCAS builds remain LGPL-clean.

**中文**

`ocas-gpl` crate（GPL-3.0-or-later 许可）为可选，从不包含在默认构建中。它承载无法在 LGPL 下发布的 GPL 专属算法。默认 oCAS 构建保持 LGPL 纯净。
