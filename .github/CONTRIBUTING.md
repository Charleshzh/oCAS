# Contributing to oCAS

Thank you for your interest in oCAS! This document outlines how to contribute effectively to the project.

## License

By contributing to oCAS, you agree that your contributions will be licensed under the **LGPL-3.0-or-later** license.

## Getting Started

1. Fork the repository.
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/ocas.git
   cd ocas
   ```
3. Install dependencies:
   - Rust 1.89 or later
   - GMP, MPFR, FLINT 3 development libraries
4. Build the workspace:
   ```bash
   cargo build
   ```
5. Run tests:
   ```bash
   cargo test --workspace
   ```

## Development Workflow

1. Create a new branch for your work:
   ```bash
   git checkout -b feature/my-feature
   ```
2. Make your changes.
3. Ensure tests pass and code is formatted:
   ```bash
   cargo fmt --check
   cargo clippy --workspace -- -D warnings
   cargo test --workspace
   ```
4. Commit with a clear message:
   ```
   feat(poly): add sparse multivariate polynomial GCD
   ```
5. Open a pull request against the `main` branch.

## Code Style

- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
- Use `cargo fmt` for formatting.
- Use meaningful variable names; avoid single-letter names except in mathematical contexts.
- Document public APIs with `rustdoc` examples.
- Keep `unsafe` blocks minimal and well-documented.

## Testing

- Add unit tests for new functions.
- Add property tests for algebraic invariants using `proptest` where appropriate.
- Include regression tests in `ocas-tests` for comparisons with SymPy/SageMath.
- Run the full test suite before submitting a PR.

## Backend and License Considerations

- Do not introduce GPL-only dependencies into the default build.
- GPL-compatible backends must be placed in the `ocas-gpl` crate and guarded by the `gpl` feature.
- When adding a new dependency, run `cargo-deny` to verify license compatibility.

## Reporting Issues

When reporting bugs, please include:

- A minimal reproducible example
- The output of `rustc --version` and `cargo --version`
- Your operating system and installed backend versions (GMP, FLINT, etc.)
- The full error message or unexpected behavior

## Communication

- Use GitHub Issues for bug reports and feature requests.
- Use GitHub Discussions for questions and design proposals.

## Acknowledgments

Contributors will be acknowledged in the project release notes.
