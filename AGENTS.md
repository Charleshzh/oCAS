# oCAS — Agent Guide / 代理指南

**English**

Cross-harness context for AI agents (and reviewers) working on oCAS. This file is
tracked in git; machine-local preferences live in an untracked `CLAUDE.md`.
Read this file before changing code, tests, CI, or docs.

## Project

oCAS (open Computer Algebra System) is a Rust computer algebra system under
**LGPL-3.0-or-later**. Default builds are portable (pure Rust); all GPL-only code
must stay in the `ocas-gpl` placeholder crate (currently no implementation).
- Edition 2024, MSRV 1.97 (pinned via `rust-toolchain.toml` to 1.97.1).
- Current version: **0.30.0** (see `Cargo.toml` / `CHANGELOG.md`).
- Repository: https://github.com/Charleshzh/oCAS · Docs: https://charleshzh.github.io/oCAS/latest/en/

## Layout and layers

Dependencies point strictly downward; no cycles.

| Layer | Crates |
|---|---|
| 5 bindings | `ocas`, `ocas-py` (PyO3, cdylib), `ocas-c` (cbindgen) |
| 4 applications | `ocas-calc` (calculus/integration), `ocas-eval` (interpreter/JIT), `ocas-parse` |
| 3 symbolic engine | `ocas-atom` (arena expression tree), `ocas-rewrite` (rules/e-graph) |
| 2 algebraic kernel | `ocas-domain`, `ocas-poly` |
| 0+1 runtime/backend | `ocas-core` (arena, errors, thread pool, FFI glue) |
| tests | `ocas-tests` (integration, correctness-vs-SymPy, benchmarks) — `publish = false` |

Accepted deviations are recorded in `docs/planning/ARCHITECTURE_CN.md` /
`ARCHITECTURE_EN.md`: `ocas-rewrite -> ocas-parse` runtime dependency,
`ocas-gpl` placeholder, Cranelift-only JIT (LLVM planned for 0.34).

## Build and quality gates

Useful commands (Windows GNU/MSYS2 or Linux; Windows specifics:
`docs/build-windows.md`, `.cargo/config.toml` is intentionally untracked):

```bash
cargo build --workspace
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
# fast gate (~4-5 min): unit tests of every crate
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
# full gate before a release: adds the ocas-tests integration suite
cargo test --workspace --exclude ocas-py --features egg
cargo check -p ocas-py
cargo deny check
python ocas-tests/scripts/check_release_consistency.py
```

CI (`.github/workflows/ci.yml`) is the source of truth and runs more: two clippy
feature sets, 3-OS test matrix, LGPL backend tests, maturin wheel + pytest + C
example, strict rustdoc, and Miri. Do not push code that cannot pass the local
fast gate; CI can take ~70-85 minutes.

Weekly `full-tests.yml` runs the `#[ignore]`d long suites in release mode with
per-test timeouts. Known gaps are quarantined explicitly (issue #55) so the
workflow turns red only on new failures; **never delete or weaken a test** to
make it green, and shrink the quarantine lists as fixes land.

## Testing conventions

- Bottom layers (`ocas-core`, `ocas-domain`, `ocas-poly`) have the highest test
  priority: correctness + stress tests, and `proptest` invariants where useful.
- Regression pins matter for every bug fix: add a deterministic test (or
  proptest/regression seed) next to the fix, ideally machine-checked (numeric
  oracle in `ocas-tests/src/integral_eval.rs`, certificate gate for integration).
- Known, unwired mechanisms stay `#[ignore]`d and are tracked in issue #55
  instead of being deleted or relaxed.
- Integration benchmark corpus (Rubi 1892) is gitignored in
  `ocas-tests/data/`; fetch it with
  `uv run python ocas-tests/scripts/fetch_rubi_corpus.py`. The bench exits 0
  when the corpus is missing — automation must inspect the report fields
  (`solved`, `fallback`, `certified_rate`, `verify_mismatches`, `wall time`).

## Release process

Automated, CI-gated; the root `vX.Y.Z` tag is the release switch. One root
`CHANGELOG.md` is the only changelog (`release-plz.toml` sets
`changelog_update = false`; no per-crate changelogs or release PRs).

1. Bump `workspace.package.version` and `CHANGELOG.md` on `main` (12 publishable
   crates; `ocas-tests` is `publish = false`), then push.
2. `release-plz.yml` waits for CI on that commit (`preflight`), then publishes to
   crates.io and creates the root `vX.Y.Z` tag with the `RELEASE_PLZ_TOKEN` PAT
   (a tag pushed with `GITHUB_TOKEN` would not trigger downstream workflows).
   The tag step is idempotent and skips when the tag already exists.
3. The tag triggers `release.yml` (C artifacts on the GitHub Release),
   `wheels.yml` (5-platform wheels + sdist to PyPI; the `pypi` environment has
   required reviewers, so `publish-pypi` waits for an owner approval) and
   `docs.yml` (`/latest/{en,zh}/` plus frozen `/v<version>/{en,zh}/` snapshots).
4. Manual fallback after a green CI commit on `main`:
   `git tag -a vX.Y.Z -m "oCAS vX.Y.Z" && git push origin vX.Y.Z`.
   `v*` tags are protected by the `protect-release-tags` ruleset (creation is
   allowed; deletion/force-push is not — temporarily disable the ruleset to roll
   back).
5. Keep the version surfaces in sync: `Cargo.toml`, `CHANGELOG.md`, both
   READMEs, book version tables, and this file's `Current version` line.
   `consistency.yml` enforces this on every push/PR.
6. Credentials expire: `CARGO_REGISTRY_TOKEN` (crates.io, default 90 days) and
   `RELEASE_PLZ_TOKEN` (fine-grained PAT) must be rotated before expiry.

## Documentation conventions

- Top-level Markdown is bilingual (English first, then Chinese, blank line
  between paragraphs); paired documents (`README.md` + `docs/README_CN.md`,
  book `en/` + `zh/`) must be updated together.
- Exception: some planning/SPEC documents are Chinese-only
  (`BENCHMARK_RESULTS_CN.md`, `BENCHMARK_SUITE_CN.md`, `COMPETITIVE_MATRIX_CN.md`,
  `ALGEBRAIC_RT_SPEC_CN.md`, `MULTI_RADICAL_SPEC_CN.md`, `INV_SUBST_SPEC_CN.md`,
  `LOGPART_RT_SPEC_CN.md`, `POLYLOG_SPEC_CN.md`).
- Do not silently change measured benchmark numbers; record measurement date,
  command, and environment, and keep "measurement date" vs "release date"
  distinguishable.

## Operating rules for agents

- Git writes (commit/merge/push/rebase) and any remote-side action require the
  user's explicit consent in this project. A user-approved release runbook may
  authorize the specific push/tag/gh operations it lists.
- Never commit secrets, credentials, local env files (`.cargo/config.toml`),
  `CLAUDE.md`, or generated benchmark data (`ocas-tests/data/`).
- Prefer evidence over assumption: check source, docs, or official references;
  cite file paths/line numbers when reporting findings.
- Current release state: 0.30.0 published on crates.io / PyPI / GitHub Releases /
  docs (2026-10-25); 0.28.0/0.29.0 backfilled. Next line: 0.31a (polylog/Li2
  head + oracle + Li2 skeleton + W6a branch fix) and 0.31b (SPEC-first Meijer G /
  ProvedNonElementary / algebraic tower, deferred to 0.32/0.33).
- Open issues: #55 = `full-tests` triage / known gaps; #21 = `paste` RUSTSEC
  advisory tracking.

---

**中文**

面向 AI 代理（与评审者）的跨 harness 上下文，随仓库跟踪；机器本地的偏好放在未跟踪
的 `CLAUDE.md`。修改代码、测试、CI 或文档前先读本文件。

## 项目

oCAS 是 Rust 编写的计算机代数系统，许可证 **LGPL-3.0-or-later**。默认构建可移植
（纯 Rust）；GPL 专属代码只能放在 `ocas-gpl` 占位 crate（当前无实现）。
- Edition 2024，MSRV 1.97（`rust-toolchain.toml` 固定 1.97.1）。
- 当前版本：**0.30.0**（见 `Cargo.toml` / `CHANGELOG.md`）。
- 仓库：https://github.com/Charleshzh/oCAS · 文档：https://charleshzh.github.io/oCAS/latest/zh/

## 结构与分层

依赖严格向下，禁止环。

| 层 | crate |
|---|---|
| 5 绑定 | `ocas`、`ocas-py`（PyO3，cdylib）、`ocas-c`（cbindgen） |
| 4 应用 | `ocas-calc`、`ocas-eval`、`ocas-parse` |
| 3 符号引擎 | `ocas-atom`、`ocas-rewrite` |
| 2 代数核 | `ocas-domain`、`ocas-poly` |
| 0+1 运行时/后端 | `ocas-core` |
| 测试 | `ocas-tests`（`publish = false`） |

已接受的例外记录在 `docs/planning/ARCHITECTURE_CN.md` / `ARCHITECTURE_EN.md`：
`ocas-rewrite -> ocas-parse` 运行期依赖、`ocas-gpl` 占位、JIT 仅 Cranelift
（LLVM 规划在 0.34）。

## 构建与质量门禁

```bash
cargo build --workspace
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
# fast gate（约 4–5 分钟）：全部 crate 单测
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
# 发布前 full gate：加 ocas-tests 集成套件
cargo test --workspace --exclude ocas-py --features egg
cargo check -p ocas-py
cargo deny check
python ocas-tests/scripts/check_release_consistency.py
```

CI（`.github/workflows/ci.yml`）是最终标准，还包含：双档 clippy、三平台测试、
LGPL 后端测试、maturin wheel + pytest + C 示例、严格 rustdoc、Miri；一次约
70–85 分钟。不得 push 无法通过本地 fast gate 的代码。

每周 `full-tests.yml` 以 release 模式跑 `#[ignore]` 长测试，逐测试超时；已知缺口
显式隔离并记录在 issue #55，该 workflow 只在出现新失败时变红。**不要删除或弱化
测试**来让它变绿；修好后收缩隔离列表。

## 测试约定

- 底层（`ocas-core`、`ocas-domain`、`ocas-poly`）优先级最高：正确性 + 压力测试，
  必要时加 `proptest` 不变量。
- 每个 bug 修复都要有确定性回归测试（或 proptest 回归种子），优先用机器可验证的
  oracle（积分数值 oracle、证书门）。
- 已知未接线机制保持 `#[ignore]` 并记录在 issue #55，而不是删除/放宽断言。
- Rubi 1892 语料在 gitignore 的 `ocas-tests/data/`；用
  `uv run python ocas-tests/scripts/fetch_rubi_corpus.py` 获取。缺语料时 bench
  exit 0，自动化必须检查 report 字段（`solved`、`fallback`、`certified_rate`、
  `verify_mismatches`、`wall time`）。

## 发布流程

自动化、CI 门禁；顶层 `vX.Y.Z` tag 是发布开关。唯一变更日志是根 `CHANGELOG.md`
（`release-plz.toml` 设 `changelog_update = false`，不生成 per-crate CHANGELOG/发版 PR）。

1. 在 `main` 上 bump `workspace.package.version` 与 `CHANGELOG.md`（12 个可发布
   crate；`ocas-tests` 为 `publish = false`），然后 push。
2. `release-plz.yml` 先等待该 commit 的 CI 绿（preflight），再发布 crates.io，并用
   `RELEASE_PLZ_TOKEN`（PAT）创建顶层 tag；用 `GITHUB_TOKEN` 推 tag 不会触发下游
   workflow。tag 步骤幂等（tag 已存在则跳过）。
3. tag 触发 `release.yml`（C 产物）、`wheels.yml`（5 平台 wheel + sdist 到 PyPI；
   `pypi` environment 有 required reviewers，`publish-pypi` 会等待 owner 批准）、
   `docs.yml`（`/latest/{en,zh}/` 与冻结的 `/v<version>/{en,zh}/`）。
4. 手工兜底（对已在 `main` 且 CI 通过的提交）：
   `git tag -a vX.Y.Z -m "oCAS vX.Y.Z" && git push origin vX.Y.Z`。`v*` 受
   `protect-release-tags` ruleset 保护（可创建、不可删除/强推；需回滚时临时禁用）。
5. 保持版本表面一致：`Cargo.toml`、`CHANGELOG.md`、两个 README、书籍版本表、本文件
   的 `Current version` 行；`consistency.yml` 在每次 push/PR 强制校验。
6. 凭据会过期：`CARGO_REGISTRY_TOKEN`（crates.io，默认 90 天）与
   `RELEASE_PLZ_TOKEN`（fine-grained PAT）需在到期前轮换。

## 文档约定

- 顶层 Markdown 双语（英文在前、中文在后，段间空行）；成对文档（`README.md` +
  `docs/README_CN.md`、book `en/` + `zh/`）必须同步修改。
- 例外：部分规划/规格文档目前仅中文（`BENCHMARK_RESULTS_CN.md`、
  `BENCHMARK_SUITE_CN.md`、`COMPETITIVE_MATRIX_CN.md`、`ALGEBRAIC_RT_SPEC_CN.md`、
  `MULTI_RADICAL_SPEC_CN.md`、`INV_SUBST_SPEC_CN.md`、`LOGPART_RT_SPEC_CN.md`、
  `POLYLOG_SPEC_CN.md`）。
- 不要静默改基准数字；记录测量日期、命令与环境，并区分「测量日」与「发布日」。

## 代理操作约定

- 本项目的 git 写操作（commit/merge/push/rebase）与任何远程操作需要用户明确同意；
  经用户批准的发布 runbook 可授权其中列出的 push/tag/gh 操作。
- 不要提交密钥、凭据、本地环境文件（`.cargo/config.toml`）、`CLAUDE.md` 或生成的
  基准数据（`ocas-tests/data/`）。
- 以证据说话：查源码、文档或官方资料；汇报时给文件路径/行号。
- 当前发布状态：0.30.0 已于 2026-10-25 全渠道发布（crates.io/PyPI/GitHub Releases/
  文档），0.28.0/0.29.0 已回填。下一条线：0.31a（polylog/Li₂ 头 + oracle + Li₂
  骨架 + W6a 分支修复）与 0.31b（SPEC 先行；Meijer G / ProvedNonElementary / 代数塔
  移到 0.32/0.33）。
- open issue：#55（`full-tests` 已知缺口跟踪）、#21（`paste` RUSTSEC 跟踪）。
