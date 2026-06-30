# Building oCAS on Windows / 在 Windows 上构建 oCAS

**English**

oCAS defaults to the MSYS2 MINGW64 toolchain on Windows. This is controlled by
[`rust-toolchain.toml`](../rust-toolchain.toml) (host toolchain) and
[`.cargo/config.toml`](../.cargo/config.toml) (default target). If you prefer
MSVC, use `cargo build --target x86_64-pc-windows-msvc --no-default-features`.

**中文**

oCAS 在 Windows 上默认使用 MSYS2 MINGW64 工具链。默认行为由
[`rust-toolchain.toml`](../rust-toolchain.toml)（主机工具链）和
[`.cargo/config.toml`](../.cargo/config.toml)（默认目标）共同控制。如需使用
MSVC，请执行 `cargo build --target x86_64-pc-windows-msvc --no-default-features`。

---

## Option 1: MSYS2 / MINGW64 (Recommended) / 方案 1：MSYS2 / MINGW64（推荐）

**English**

1. Install [MSYS2](https://www.msys2.org/).
2. Open the MINGW64 terminal and install dependencies:
   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
   ```
3. Open the project from the MINGW64 terminal so that `sh`, `gcc`, and
   `pkg-config` are on `PATH`:
   ```bash
   cd /d/rust/oCAS/ocas
   cargo build --release
   ```
   The `rust-toolchain.toml` file selects the MINGW64 host toolchain
   automatically, so you do not need to run `rustup default` manually.

4. To use the GMP-backed `Integer`/`Rational` backend or the MPFR-backed
   `RealBall` backend, make sure the MSYS2 `pkg-config` is found before any
   other `pkg-config` on `PATH`:
   ```bash
   cargo test -p ocas-domain --features gmp
   cargo test -p ocas-domain --features mpfr
   cargo test -p ocas-poly --features gmp
   ```
   From PowerShell outside the MSYS2 shell, the most reliable way is to run
   cargo inside the MINGW64 subshell via `msys2_shell.cmd`:
   ```powershell
   C:\msys64\msys2_shell.cmd -defterm -here -no-start -mingw64 -c "cd /d/rust/oCAS/ocas && cargo test -p ocas-domain --features gmp"
   C:\msys64\msys2_shell.cmd -defterm -here -no-start -mingw64 -c "cd /d/rust/oCAS/ocas && cargo test -p ocas-domain --features mpfr"
   ```
   This initializes the full MSYS2 environment and avoids `pkg-config`
   conflicts with tools like Strawberry Perl or Git for Windows.

**中文**

1. 安装 [MSYS2](https://www.msys2.org/)。
2. 打开 MINGW64 终端并安装依赖：
   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
   ```
3. 从 MINGW64 终端打开项目，确保 `sh`、`gcc` 与 `pkg-config` 在 `PATH` 中：
   ```bash
   cd /d/rust/oCAS/ocas
   cargo build --release
   ```
   项目中的 `rust-toolchain.toml` 会自动选择 MINGW64 主机工具链，无需手动执行
   `rustup default`。

4. 若要使用 GMP 后端的 `Integer`/`Rational` 或 MPFR 后端的严格
   `RealBall`，请确保 MSYS2 的 `pkg-config` 在 `PATH` 中最先被找到：
   ```bash
   cargo test -p ocas-domain --features gmp
   cargo test -p ocas-domain --features mpfr
   cargo test -p ocas-poly --features gmp
   ```
   若在 MSYS2 外部使用 PowerShell，最可靠的方式是通过 `msys2_shell.cmd` 启动
   MINGW64 子 shell 执行 cargo 命令：
   ```powershell
   C:\msys64\msys2_shell.cmd -defterm -here -no-start -mingw64 -c "cd /d/rust/oCAS/ocas && cargo test -p ocas-domain --features gmp"
   C:\msys64\msys2_shell.cmd -defterm -here -no-start -mingw64 -c "cd /d/rust/oCAS/ocas && cargo test -p ocas-domain --features mpfr"
   ```
   这会初始化完整的 MSYS2 环境，避免 PowerShell 中混用 Strawberry Perl、
   Git for Windows 等自带的 `pkg-config` 导致构建失败。

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

If you only want to build the pure-Rust components (default workspace build):

```powershell
cargo check --workspace --exclude ocas-py
```

This skips the Python bindings and avoids GMP/MPFR/FLINT. MSVC works fine for
this configuration. Performance will be significantly lower if you later enable
backends.

**中文**

如果您只想构建纯 Rust 组件（默认工作区构建）：

```powershell
cargo check --workspace --exclude ocas-py
```

这会跳过 Python 绑定并避免 GMP/MPFR/FLINT。此配置下 MSVC 也能正常工作。
如果后续启用后端，性能会显著降低。

---

## Known Issues / 已知问题

**English**

- MSVC + GMP is not officially supported by GMP. Use MINGW64 or vcpkg.
- The `flint` feature is experimental and currently does not build on Windows
  because `flint3-sys` relies on POSIX tooling and `libc::pthread_mutex_t`,
  which are not available in the Windows `libc` crate. Use a Linux runner or
  WSL for FLINT-backed tests.
- When building with the `gmp`/`mpfr` features from outside the MSYS2 shell,
  make sure `C:/msys64/usr/bin` and `C:/msys64/mingw64/bin` are on `PATH` so
  that `sh`, `gcc`, and `pkg-config` are found.

**中文**

- GMP 官方不支持 MSVC + GMP，请使用 MINGW64 或 vcpkg。
- `flint` 特性为实验性，目前在 Windows 上无法构建，因为 `flint3-sys` 依赖
  POSIX 工具链以及 Windows 版 `libc` crate 未提供的 `libc::pthread_mutex_t`。
  FLINT 后端测试请在 Linux 或 WSL 下运行。
- 若不在 MSYS2 shell 内构建而启用 `gmp`/`mpfr` 特性，请确保
  `C:/msys64/usr/bin` 与 `C:/msys64/mingw64/bin` 在 `PATH` 中，以便找到
  `sh`、`gcc` 与 `pkg-config`。

---

## Getting Help / 获取帮助

**English**

If you encounter build issues, please open a GitHub issue with the full build log.

**中文**

如果遇到构建问题，请创建 GitHub issue 并附上完整构建日志。
