# Building on Windows / 在 Windows 上构建

**English**

The recommended Windows toolchain is MSYS2 MINGW64 (GNU target), which
provides GMP/MPFR for the optional numerical backends. MSVC builds work for
the pure-Rust default (no backends) configuration.

**中文**

Windows 推荐工具链为 MSYS2 MINGW64（GNU 目标），可为可选数值后端提供 GMP/MPFR。MSVC 构建适用于纯 Rust 默认（无后端）配置。

## MSYS2 MINGW64 setup / MSYS2 MINGW64 配置

```bash
# In the MSYS2 MINGW64 shell
pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
```

Then build with backends:

```bash
cargo build -p ocas --features gmp,mpfr
```

## MSVC (no backends) / MSVC（无后端）

```powershell
cargo build --target x86_64-pc-windows-msvc --no-default-features
```

> The `flint` feature is not supported on Windows due to upstream
> `flint3-sys` limitations (POSIX-only types).

**English**

For the full step-by-step guide including PATH configuration and
troubleshooting, see
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)
in the repository.

**中文**

包含 PATH 配置与问题排查的完整分步指南，见仓库中的
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)。
