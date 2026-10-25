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
5. Run the fast quality gate (see [Quality Gates](#quality-gates--质量门禁) below):
   ```bash
   cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
   ```

The layered-architecture rules are documented in `docs/planning/ARCHITECTURE_EN.md` (Chinese: `docs/planning/ARCHITECTURE_CN.md`).

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
5. 运行 fast 质量门禁（见下文[质量门禁](#quality-gates--质量门禁)）：
   ```bash
   cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
   ```

分层架构规则见 `docs/planning/ARCHITECTURE_CN.md`（英文：`docs/planning/ARCHITECTURE_EN.md`）。

---

## Quality Gates / 质量门禁

**English**

Two tiers: the **fast** gate for every change, the **full** gate before opening a pull request (CI parity).

Fast gate:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
```

Full gate (adds the cross-system correctness suite in `ocas-tests`):

```bash
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

Tests marked `#[ignore]` (slow property tests, performance acceptance tests, work-in-progress fixes) are outside both gates; the scheduled long-test workflow (`.github/workflows/full-tests.yml`) runs them, and they can be invoked with `cargo test ... -- --ignored` (add `--release` for the performance acceptance tests).

**中文**

分为两级：每次改动运行 **fast（快速）** 门禁，提交 Pull Request 前运行 **full（完整）** 门禁（与 CI 一致）。

Fast 门禁：

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
```

Full 门禁（额外包含 `ocas-tests` 中的跨系统正确性测试套件）：

```bash
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

标记为 `#[ignore]` 的测试（慢速属性测试、性能验收测试、进行中的修复）不在两级门禁内；它们由定时长测工作流（`.github/workflows/full-tests.yml`）运行，可通过 `cargo test ... -- --ignored` 显式执行（性能验收测试加 `--release`）。

---

## Development Workflow / 开发流程

**English**

1. Create a new branch for your work:
   ```bash
   git checkout -b feature/my-feature
   ```
2. Make your changes.
3. Run the quality gates described above; at minimum the fast gate must pass:
   ```bash
   cargo fmt --all -- --check
   cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
   cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
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
3. 运行上述质量门禁；至少 fast 门禁必须通过：
   ```bash
   cargo fmt --all -- --check
   cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
   cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
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
- Run the full quality gate before submitting a PR.

**中文**

- 为新函数添加单元测试。
- 在适当时使用 `proptest` 为代数不变量添加属性测试。
- 在 `ocas-tests` 中包含与 SymPy/SageMath 对比的回归测试。
- 提交 PR 前运行 full 质量门禁。

---

## Backend and License Considerations / 后端与许可证注意事项

**English**

- Do not introduce GPL-only dependencies into the default build.
- `ocas-gpl` is an independent placeholder crate for future GPL-3.0-or-later backends. The top-level `ocas` facade does not expose a `gpl` feature, so a GPL backend has to be opted into by depending on `ocas-gpl` directly.
- When adding a new dependency, run `cargo-deny` to verify license compatibility.

**中文**

- 不要向默认构建引入 GPL-only 依赖。
- `ocas-gpl` 是用于未来 GPL-3.0-or-later 后端的独立占位 crate。顶层 `ocas` facade 不再暴露 `gpl` feature，因此 GPL 后端必须通过直接依赖 `ocas-gpl` 来显式启用。
- 添加新依赖时，运行 `cargo-deny` 验证许可证兼容性。

---

## Release Process / 发布流程

**English**

Releases are automated and CI-gated; the root `vX.Y.Z` tag is the release switch. There is only one changelog — the root `CHANGELOG.md` (`release-plz.toml` sets `changelog_update = false`, so no per-crate changelogs or release PRs are generated).

1. Bump `workspace.package.version` and `CHANGELOG.md` on `main` yourself (12 publishable crates; `ocas-tests` is `publish = false`), then push. If the version did not change, release-plz is a no-op.
2. The `release-plz` workflow waits for CI to be green on that exact commit (`preflight`), publishes the crates to crates.io in dependency order, and then creates the root `vX.Y.Z` tag. Tag creation needs the `RELEASE_PLZ_TOKEN` secret (a fine-grained PAT with `Contents: Read and write`); a tag pushed with the default `GITHUB_TOKEN` does not trigger the downstream workflows.
3. The tag triggers `release.yml` (validates the tag against the workspace version, dry-runs packaging and attaches the C artifacts to the project-level GitHub Release), `wheels.yml` (builds wheels + sdist with maturin and publishes to PyPI; pre-release tags go to TestPyPI — note that the `pypi` environment asks an owner to approve the publish job), and `docs.yml` (publishes `/latest/{en,zh}/` plus a frozen `/v<version>/{en,zh}/` snapshot).
4. Manual fallback for a commit that is on `main` and green in CI:
   ```bash
   git tag -a vX.Y.Z -m "oCAS vX.Y.Z"
   git push origin vX.Y.Z
   ```
   `v*` tags are protected by the `protect-release-tags` ruleset: creation is allowed, deletion/force-push is not.
5. Keep `Cargo.toml`, `CHANGELOG.md`, the READMEs, the book version tables and the `Current version` line of `AGENTS.md` in sync; `consistency.yml` enforces this on every push and pull request. Rotate `CARGO_REGISTRY_TOKEN` / `RELEASE_PLZ_TOKEN` before they expire.

**中文**

发布已自动化且有 CI 门禁；顶层 `vX.Y.Z` 标签就是发布开关。变更日志只有根 `CHANGELOG.md` 一份（`release-plz.toml` 设 `changelog_update = false`，不生成 per-crate CHANGELOG，也不创建发版 PR）。

1. 自行在 `main` 上提升 `workspace.package.version` 并更新 `CHANGELOG.md`（12 个可发布 crate；`ocas-tests` 为 `publish = false`），然后 push。版本未变化时 release-plz 自然 no-op。
2. `release-plz` workflow 先等待该提交的 CI 变绿（`preflight`），再按依赖序发布各 crate 到 crates.io，最后创建顶层 `vX.Y.Z` 标签。打标签需要 `RELEASE_PLZ_TOKEN` secret（具备 `Contents: Read and write` 权限的 fine-grained PAT）；用默认 `GITHUB_TOKEN` 推送的标签不会触发下游 workflow。
3. 标签触发 `release.yml`（校验标签与工作区版本一致、做打包试运行，并把 C 绑定产物附加到项目级 GitHub Release）、`wheels.yml`（用 maturin 构建 wheel + sdist 并发布到 PyPI；预发布标签发到 TestPyPI——注意 `pypi` environment 会要求 owner 批准发布 job）以及 `docs.yml`（发布 `/latest/{en,zh}/` 与冻结的 `/v<version>/{en,zh}/` 版本快照）。
4. 手动兜底（对已在 `main` 上且 CI 通过的提交）：
   ```bash
   git tag -a vX.Y.Z -m "oCAS vX.Y.Z"
   git push origin vX.Y.Z
   ```
   `v*` 受 `protect-release-tags` ruleset 保护：可创建，不可删除/强推。
5. 保持 `Cargo.toml`、`CHANGELOG.md`、两个 README、书籍版本表以及 `AGENTS.md` 的 `Current version` 行同步；`consistency.yml` 会在每次 push 与 Pull Request 上强制执行该检查。请在到期前轮换 `CARGO_REGISTRY_TOKEN` / `RELEASE_PLZ_TOKEN`。

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
