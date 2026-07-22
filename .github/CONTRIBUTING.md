# Contributing to oCAS / 参与 oCAS

**English**

Thank you for your interest in oCAS! This document outlines how to contribute effectively to the project.

**中文**

感谢您对 oCAS 的兴趣！本文档说明如何高效地为项目做出贡献。

---

## License / 许可证

**English**

By contributing to oCAS, you agree that your contributions will be licensed under the **LGPL-3.0-or-later** license.

**中文**

向 oCAS 贡献代码即表示您同意您的贡献采用 **LGPL-3.0-or-later** 许可证。

---

## Getting Started / 开始

**English**

1. Fork the repository.
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/ocas.git
   cd ocas
   ```
3. Install dependencies:
   - Rust 1.97 or later
   - GMP, MPFR, FLINT 3 development libraries
4. Build the workspace:
   ```bash
   cargo build
   ```
5. Run tests:
   ```bash
   cargo test --workspace --exclude ocas-py
   ```

**中文**

1. Fork 仓库。
2. 克隆您的 fork：
   ```bash
   git clone https://github.com/YOUR_USERNAME/ocas.git
   cd ocas
   ```
3. 安装依赖：
   - Rust 1.97 或更高版本
   - GMP、MPFR、FLINT 3 开发库
4. 构建 workspace：
   ```bash
   cargo build
   ```
5. 运行测试：
   ```bash
   cargo test --workspace
   ```

---

## Development Workflow / 开发流程

**English**

1. Create a new branch for your work:
   ```bash
   git checkout -b feature/my-feature
   ```
2. Make your changes.
3. Ensure tests pass and code is formatted:
   ```bash
   cargo fmt --check
   cargo clippy --workspace -- -D warnings
   cargo test --workspace --exclude ocas-py
   ```
4. Commit with a clear message:
   ```
   feat(poly): add sparse multivariate polynomial GCD
   ```
5. Open a pull request against the `main` branch.

**中文**

1. 为您的改动创建新分支：
   ```bash
   git checkout -b feature/my-feature
   ```
2. 进行修改。
3. 确保测试通过且代码已格式化：
   ```bash
   cargo fmt --check
   cargo clippy --workspace -- -D warnings
   cargo test --workspace --exclude ocas-py
   ```
4. 使用清晰的提交信息：
   ```
   feat(poly): add sparse multivariate polynomial GCD
   ```
5. 向 `main` 分支提交 Pull Request。

---

## Code Style / 代码风格

**English**

- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
- Use `cargo fmt` for formatting.
- Use meaningful variable names; avoid single-letter names except in mathematical contexts.
- Document public APIs with `rustdoc` examples.
- Keep `unsafe` blocks minimal and well-documented.

**中文**

- 遵循 [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)。
- 使用 `cargo fmt` 进行格式化。
- 使用有意义的变量名；除非在数学上下文中，避免单字母名称。
- 使用 `rustdoc` 示例文档化公共 API。
- 保持 `unsafe` 块最小化并充分注释。

---

## Testing / 测试

**English**

- Add unit tests for new functions.
- Add property tests for algebraic invariants using `proptest` where appropriate.
- Include regression tests in `ocas-tests` for comparisons with SymPy/SageMath.
- Run the full test suite before submitting a PR.

**中文**

- 为新函数添加单元测试。
- 在适当时使用 `proptest` 为代数不变量添加属性测试。
- 在 `ocas-tests` 中包含与 SymPy/SageMath 对比的回归测试。
- 提交 PR 前运行完整测试套件。

---

## Backend and License Considerations / 后端与许可证注意事项

**English**

- Do not introduce GPL-only dependencies into the default build.
- GPL-compatible backends must be placed in the `ocas-gpl` crate and guarded by the `gpl` feature.
- When adding a new dependency, run `cargo-deny` to verify license compatibility.

**中文**

- 不要向默认构建引入 GPL-only 依赖。
- GPL 兼容后端必须放在 `ocas-gpl` crate 中，并通过 `gpl` feature 保护。
- 添加新依赖时，运行 `cargo-deny` 验证许可证兼容性。

---

## Reporting Issues / 报告问题

**English**

When reporting bugs, please include:

- A minimal reproducible example
- The output of `rustc --version` and `cargo --version`
- Your operating system and installed backend versions (GMP, FLINT, etc.)
- The full error message or unexpected behavior

**中文**

报告 bug 时，请包含：

- 最小可复现示例
- `rustc --version` 和 `cargo --version` 的输出
- 您的操作系统与已安装后端版本（GMP、FLINT 等）
- 完整错误信息或异常行为

---

## Communication / 沟通

**English**

- Use GitHub Issues for bug reports and feature requests.
- Use GitHub Discussions for questions and design proposals.

**中文**

- 使用 GitHub Issues 报告 bug 和提出功能请求。
- 使用 GitHub Discussions 提问和提交设计提案。

---

## Acknowledgments / 致谢

**English**

Contributors will be acknowledged in the project release notes.

**中文**

贡献者将在项目发布说明中得到致谢。
