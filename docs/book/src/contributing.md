# Contributing / 贡献

**English**

Contributions are welcome! oCAS follows trunk-based development on `main`
with short-lived feature branches and pull requests.

**中文**

欢迎贡献！oCAS 在 `main` 主干上进行开发，使用短生命周期特性分支与 Pull Request。

## Development setup / 开发环境

- Rust 1.89+ (Edition 2024)
- A C compiler (MSVC, GCC, or Clang)
- Optional: GMP/MPFR/FLINT system libraries for backend features

## Quality gates / 质量门

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

## Python bindings / Python 绑定

```bash
cd ocas-py
uv venv .venv
uv pip install pytest maturin
uv run maturin develop
uv run python -m pytest ../ocas-tests/tests/python/ -v
```

## Coding standards / 编码规范

- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
- Run `cargo fmt --all` before submitting.
- Prefer arena allocation for expression trees; avoid per-node `Box`/`Rc`.
- Use `thiserror` for error types; propagate with `?`.
- Backend code lives behind feature flags; GPL code only in `ocas-gpl`.
- Public APIs for stable releases must include rustdoc examples.

**English**

See the repository's `CLAUDE.md` and `ARCHITECTURE_EN.md` for the full
conventions, layered-dependency rules, and release process.

**中文**

完整约定、分层依赖规则与发布流程见仓库的 `CLAUDE.md` 与 `ARCHITECTURE_CN.md`。
