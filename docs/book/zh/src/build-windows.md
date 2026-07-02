# 在 Windows 上构建

Windows 推荐工具链为 MSYS2 MINGW64（GNU 目标），可为可选数值后端提供 GMP/MPFR。MSVC 构建适用于纯 Rust 默认（无后端）配置。

## MSYS2 MINGW64 配置

```bash
# 在 MSYS2 MINGW64 shell 中
pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
```

然后带后端构建：

```bash
cargo build -p ocas --features gmp,mpfr
```

## MSVC（无后端）

```powershell
cargo build --target x86_64-pc-windows-msvc --no-default-features
```

> 由于上游 `flint3-sys` 的限制（仅支持 POSIX 类型），`flint` feature 在 Windows 上不可用。

包含 PATH 配置与问题排查的完整分步指南，见仓库中的
[docs/build-windows.md](https://github.com/charleshzh/ocas/blob/main/docs/build-windows.md)。
