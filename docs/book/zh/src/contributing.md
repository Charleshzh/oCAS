# 贡献

欢迎贡献！oCAS 在 `main` 主干上进行开发，使用短生命周期特性分支与 Pull Request。

## 开发环境

- Rust 1.89+（Edition 2024）
- C 编译器（MSVC、GCC 或 Clang）
- 可选：用于后端特性的 GMP/MPFR/FLINT 系统库

## 质量门禁

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

## Python 绑定

```bash
cd ocas-py
uv venv .venv
uv pip install pytest maturin
uv run maturin develop
uv run python -m pytest ../ocas-tests/tests/python/ -v
```

## 编码规范

- 遵循 [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)。
- 提交前运行 `cargo fmt --all`。
- 表达式树优先使用 arena 分配；避免逐节点 `Box`/`Rc`。
- 错误类型使用 `thiserror`；用 `?` 传播。
- 后端代码隐藏在 feature flag 后；GPL 代码只能放在 `ocas-gpl`。
- 稳定版本的公共 API 必须包含 rustdoc 示例。

完整约定、分层依赖规则与发布流程见仓库的 `CLAUDE.md` 与 `ARCHITECTURE_CN.md`。
