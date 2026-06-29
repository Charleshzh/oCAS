# Building oCAS on Windows / 在 Windows 上构建 oCAS

**English**

This guide covers building oCAS on Windows. The default backend stack requires GMP, MPFR, and FLINT, which are easiest to obtain via MSYS2 or vcpkg.

**中文**

本指南介绍在 Windows 上构建 oCAS。默认后端栈需要 GMP、MPFR 与 FLINT，最简单的方式是通过 MSYS2 或 vcpkg 获取。

---

## Option 1: MSYS2 / MINGW64 (Recommended) / 方案 1：MSYS2 / MINGW64（推荐）

**English**

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

**中文**

1. 安装 [MSYS2](https://www.msys2.org/)。
2. 打开 MINGW64 终端并安装依赖：
   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr mingw-w64-x86_64-flint
   ```
3. 安装 MINGW64 目标的 Rust：
   ```bash
   rustup toolchain install stable-x86_64-pc-windows-gnu
   rustup default stable-x86_64-pc-windows-gnu
   ```
4. 构建 oCAS：
   ```bash
   cargo build --release
   ```

---

## Option 2: vcpkg / 方案 2：vcpkg

**English**

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

**中文**

1. 安装 [vcpkg](https://vcpkg.io/)。
2. 安装库：
   ```bash
   vcpkg install gmp mpfr flint
   ```
3. 设置环境变量以便 Cargo 找到库：
   ```powershell
   $env:VCPKG_ROOT = "C:\path\to\vcpkg"
   cargo build --release
   ```

---

## Option 3: No External Backends (Limited) / 方案 3：无外部后端（有限）

**English**

If you only want to build the pure-Rust components:

```powershell
cargo build --no-default-features
```

This disables GMP, MPFR, and FLINT. Performance will be significantly lower.

**中文**

如果您只想构建纯 Rust 组件：

```powershell
cargo build --no-default-features
```

这将禁用 GMP、MPFR 和 FLINT，性能会显著降低。

---

## Known Issues / 已知问题

**English**

- MSVC + GMP is not officially supported by GMP. Use MINGW64 or vcpkg.
- The `gmp` Cargo feature is disabled by default on Windows for this reason.
  To run tests or benchmarks with GMP enabled, use the MINGW64 target or a
  Linux/macOS environment.
- Make sure `pkg-config` can locate `.pc` files for GMP/MPFR/FLINT.

**中文**

- GMP 官方不支持 MSVC + GMP，请使用 MINGW64 或 vcpkg。
- 确保 `pkg-config` 可以找到 GMP/MPFR/FLINT 的 `.pc` 文件。

---

## Getting Help / 获取帮助

**English**

If you encounter build issues, please open a GitHub issue with the full build log.

**中文**

如果遇到构建问题，请创建 GitHub issue 并附上完整构建日志。
