# Contributing

Contributions are welcome! oCAS follows trunk-based development on `main`
with short-lived feature branches and pull requests.

## Development setup

- Rust 1.97+ (Edition 2024)
- A C compiler (MSVC, GCC, or Clang)
- Optional: GMP/MPFR/FLINT system libraries for backend features

## Quality gates

Two tiers: run the **fast** gate on every change, and the **full** gate before
opening a pull request (CI parity).

Fast (default feedback loop):

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --features egg --exclude ocas-py -- -D warnings
cargo test --workspace --exclude ocas-py --exclude ocas-tests --features egg
```

Full (CI parity — adds the cross-system correctness suite in `ocas-tests`):

```bash
cargo test --workspace --exclude ocas-py --features egg
cargo deny check
```

Tests marked `#[ignore]` (slow property tests, performance acceptance tests,
work-in-progress fixes) are deliberately outside both gates; the scheduled
long-test sweep (`.github/workflows/full-tests.yml`, weekly) runs them, and they
can also be invoked explicitly:

```bash
cargo test --workspace --exclude ocas-py --features egg -- --ignored
cargo test --release -p ocas-domain --lib -- --ignored   # performance acceptance
```

## Coding standards

- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).
- Run `cargo fmt --all` before submitting.
- Prefer arena allocation for expression trees; avoid per-node `Box`/`Rc`.
- Use `thiserror` for error types; propagate with `?`.
- Backend code lives behind feature flags. `ocas-gpl` is an independent
  placeholder crate for future GPL-3.0-or-later backends; the top-level `ocas`
  facade does not expose a `gpl` feature, so a GPL backend has to be opted into
  by depending on `ocas-gpl` directly.
- Public APIs for stable releases must include rustdoc examples.

See `docs/planning/ARCHITECTURE_EN.md` for the full conventions and the
layered-dependency rules.

## Python bindings

```bash
cd ocas-py
uv venv .venv
uv pip install pytest maturin
uv run maturin develop
uv run python -m pytest ../ocas-tests/tests/python/ -v
```

## Release process

Releases are automated and CI-gated; the root `vX.Y.Z` tag is the release switch.
There is only one changelog — the root `CHANGELOG.md` (`release-plz.toml` sets
`changelog_update = false`; no per-crate changelogs or release PRs).

1. **Version bump and tag creation.** Bump `workspace.package.version` and
   `CHANGELOG.md` on `main` yourself (12 publishable crates; `ocas-tests` is
   `publish = false`), then push. The `release-plz` workflow waits for CI to be
   green on that exact commit, publishes the crates in dependency order, and
   creates the root `vX.Y.Z` tag. Tag creation needs the `RELEASE_PLZ_TOKEN`
   secret (a fine-grained PAT with `Contents: Read and write`): a tag pushed with
   the default `GITHUB_TOKEN` does not trigger the downstream workflows. If the
   version did not change, release-plz is a no-op.
2. **Tag fan-out.**
   - `release.yml` validates the tag against the workspace version, dry-runs
     `cargo publish` for every crate, and attaches the C artifacts to the
     project-level GitHub Release.
   - `wheels.yml` builds the Python wheels and sdist with maturin and publishes
     to PyPI (pre-release tags go to TestPyPI); the `pypi` environment asks an
     owner to approve the publish job.
   - `docs.yml` rebuilds the bilingual mdBook and publishes both
     `/latest/{en,zh}/` and a frozen `/v<version>/{en,zh}/` snapshot to GitHub
     Pages.
3. **Manual fallback.** If the automation fails, tag a commit that is on
   `main` and green in CI:

   ```bash
   git tag -a vX.Y.Z -m "oCAS vX.Y.Z"
   git push origin vX.Y.Z
   ```

   Do not run the publish steps locally; publishing is done by CI. `v*` tags are
   protected, so creation is allowed but deletion/force-push is not.

4. **Version surfaces.** Keep `Cargo.toml`, `CHANGELOG.md`, the READMEs, the
   version tables in `docs/book/*/src/introduction.md`, and the `Current version`
   line in `AGENTS.md` in sync. `consistency.yml` enforces this on every push and
   pull request via `ocas-tests/scripts/check_release_consistency.py`. Rotate
   `CARGO_REGISTRY_TOKEN` / `RELEASE_PLZ_TOKEN` before they expire.
