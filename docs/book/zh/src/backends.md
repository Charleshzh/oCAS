# 数值后端

oCAS 默认提供纯 Rust 构建（使用 `num-bigint`/`num-rational`），并通过 feature flag 提供可选的高性能后端。默认构建仅使用 LGPL 兼容依赖；GPL 专属代码隔离在可选的 `ocas-gpl` crate 中。

## 功能开关

| Feature | Crate | 后端 | 许可证 |
|---|---|---|---|
| `gmp` | `ocas-core`、`ocas-domain` | 经 `rug` 的 GMP 整数/有理数 | LGPL |
| `mpfr` | `ocas-domain` | 经 `rug` 的 MPFR 实数球 | LGPL |
| `flint` | `ocas-poly` | FLINT 3（仅 Linux/WSL） | LGPL |
| `jit` | `ocas-eval` | Cranelift JIT 编译 | Apache-2.0 |
| `simd` | `ocas-eval` | SIMD 向量化求值 | MIT |
| `egg` | `ocas-rewrite` | E-graph 等式饱和 | MIT |

## 启用后端

```bash
# GMP + MPFR 后端（Windows 上推荐通过 MSYS2 MINGW64）
cargo build -p ocas --features gmp,mpfr

# FLINT 后端（仅 Linux/WSL；需要系统 libflint）
cargo build -p ocas --features gmp,mpfr,flint

# JIT + SIMD 求值
cargo build -p ocas --features jit,simd
```

Windows 上通过 MSYS2 安装 GMP/MPFR 的步骤见
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)。

## GPL 后端

`ocas-gpl` crate（GPL-3.0-or-later 许可）为可选，从不包含在默认构建中。它承载无法在 LGPL 下发布的 GPL 专属算法。默认 oCAS 构建保持 LGPL 纯净。
