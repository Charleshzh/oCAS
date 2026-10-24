# 贡献

欢迎贡献！oCAS 在 `main` 主干上进行开发，使用短生命周期特性分支与 Pull Request。

## 开发环境

- Rust 1.97+（Edition 2024）
- C 编译器（MSVC、GCC 或 Clang）
- 可选：用于后端特性的 GMP/MPFR/FLINT 系统库

## 质量门禁

分为两级：每次改动都运行 **fast（快速）** 门禁，提交 Pull Request 前运行
**full（完整）** 门禁（与 CI 一致）。

Fast（默认反馈回路）：

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
```

Full（与 CI 一致——额外包含 `ocas-tests` 中的跨系统正确性测试套件）：

```bash
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

标记为 `#[ignore]` 的测试（慢速属性测试、性能验收测试、进行中的修复）刻意不在
两级门禁内，由定时长测任务（`.github/workflows/full-tests.yml`，每周）执行，
也可以显式运行：

```bash
cargo test --workspace --exclude ocas-py --features egg -- --ignored
cargo test --release -p ocas-domain --lib -- --ignored   # 性能验收
```

## 编码规范

- 遵循 [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)。
- 提交前运行 `cargo fmt --all`。
- 表达式树优先使用 arena 分配；避免逐节点 `Box`/`Rc`。
- 错误类型使用 `thiserror`；用 `?` 传播。
- 后端代码隐藏在 feature flag 后。`ocas-gpl` 是用于未来 GPL-3.0-or-later 后端的
  独立占位 crate；顶层 `ocas` facade 不再暴露 `gpl` feature，因此 GPL 后端必须
  通过直接依赖 `ocas-gpl` 来显式启用。
- 稳定版本的公共 API 必须包含 rustdoc 示例。

完整约定与分层依赖规则见 `docs/planning/ARCHITECTURE_CN.md`。

## Python 绑定

```bash
cd ocas-py
uv venv .venv
uv pip install pytest maturin
uv run maturin develop
uv run python -m pytest ../ocas-tests/tests/python/ -v
```

## 发布流程

发布已自动化；顶层 `vX.Y.Z` 标签就是发布开关。

1. **发版 PR 与打标签。** `main` 前进后，release-plz
   （`.github/workflows/release-plz.yml`）会创建（并随后更新）一个提升
   `workspace.package.version` 的发版 PR。该 PR 合并后，workflow 先等待合并提交
   的 CI 变绿，再按依赖序发布各 crate，并创建顶层 `vX.Y.Z` 标签。打标签需要
   `RELEASE_PLZ_TOKEN` secret（具备 `contents: write` 权限的 PAT）：用默认
   `GITHUB_TOKEN` 推送的标签不会触发下游 workflow。
2. **标签触发的下游流程。**
   - `release.yml` 校验标签与工作区版本一致，对每个 crate 做 `cargo publish`
     试运行，并把 C 绑定产物附加到项目级 GitHub Release。
   - `wheels.yml` 用 maturin 构建 Python wheel 与 sdist 并发布到 PyPI（预发布
     标签发到 TestPyPI）。
   - `docs.yml` 重新构建双语 mdBook，并把 `/latest/{en,zh}/` 与冻结的
     `/v<version>/{en,zh}/` 版本快照部署到 GitHub Pages。
3. **手动兜底。** 若自动化失败，请对已在 `main` 上且 CI 通过的提交打标签：

   ```bash
   git tag -a vX.Y.Z -m "oCAS vX.Y.Z"
   git push origin vX.Y.Z
   ```

   不要在本地执行发布步骤；发布由 CI 完成。

4. **版本一致性。** 保持 `Cargo.toml`、`CHANGELOG.md`、两个 README 与
   `docs/book/*/src/introduction.md` 中的版本表同步。`consistency.yml` 会在每次
   push 与 Pull Request 上通过 `ocas-tests/scripts/check_release_consistency.py`
   强制执行该检查。
