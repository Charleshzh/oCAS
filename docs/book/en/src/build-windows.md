# Building on Windows

The recommended Windows toolchain is MSYS2 MINGW64 (GNU target), which
provides GMP/MPFR for the optional numerical backends. MSVC builds work for
the pure-Rust default (no backends) configuration.

## MSYS2 MINGW64 setup

```bash
# In the MSYS2 MINGW64 shell
pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
```

Then build with backends:

```bash
cargo build -p ocas --features gmp,mpfr
```

## MSVC (no backends)

```powershell
cargo build --target x86_64-pc-windows-msvc --no-default-features
```

> The `flint` feature is not supported on Windows due to upstream
> `flint3-sys` limitations (POSIX-only types).

For the full step-by-step guide including PATH configuration and
troubleshooting, see
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)
in the repository.
