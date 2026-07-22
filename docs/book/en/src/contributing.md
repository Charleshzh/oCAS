# Contributing

Contributions are welcome! oCAS follows trunk-based development on `main`
with short-lived feature branches and pull requests.

## Development setup

- Rust 1.97+ (Edition 2024)
- A C compiler (MSVC, GCC, or Clang)
- Optional: GMP/MPFR/FLINT system libraries for backend features

## Quality gates

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

## Python bindings

```bash
cd ocas-py
uv venv .venv
uv pip install pytest maturin
uv run maturin develop
uv run python -m pytest ../ocas-tests/tests/python/ -v
```

## Coding standards

- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
- Run `cargo fmt --all` before submitting.
- Prefer arena allocation for expression trees; avoid per-node `Box`/`Rc`.
- Use `thiserror` for error types; propagate with `?`.
- Backend code lives behind feature flags; GPL code only in `ocas-gpl`.
- Public APIs for stable releases must include rustdoc examples.

See the repository's `CLAUDE.md` and `ARCHITECTURE_EN.md` for the full
conventions, layered-dependency rules, and release process.
