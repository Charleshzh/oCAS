# Building oCAS on Windows

This guide covers building oCAS on Windows. The default backend stack requires GMP, MPFR, and FLINT, which are easiest to obtain via MSYS2 or vcpkg.

## Option 1: MSYS2 / MINGW64 (Recommended)

1. Install [MSYS2](https://www.msys2.org/).
2. Open the MINGW64 terminal and install dependencies:
   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr mingw-w64-x86_64-flint
   ```
3. Install Rust for the MINGW64 target:
   ```bash
   rustup toolchain install stable-x86_64-pc-windows-gnu
   rustup default stable-x86_64-pc-windows-gnu
   ```
4. Build oCAS:
   ```bash
   cargo build --release
   ```

## Option 2: vcpkg

1. Install [vcpkg](https://vcpkg.io/).
2. Install libraries:
   ```bash
   vcpkg install gmp mpfr flint
   ```
3. Set environment variables so Cargo can find the libraries:
   ```powershell
   $env:VCPKG_ROOT = "C:\path\to\vcpkg"
   cargo build --release
   ```

## Option 3: No External Backends (Limited)

If you only want to build the pure-Rust components:

```powershell
cargo build --no-default-features
```

This disables GMP, MPFR, and FLINT. Performance will be significantly lower.

## Known Issues

- MSVC + GMP is not officially supported by GMP. Use MINGW64 or vcpkg.
- Make sure `pkg-config` can locate `.pc` files for GMP/MPFR/FLINT.

## Getting Help

If you encounter build issues, please open a GitHub issue with the full build log.
