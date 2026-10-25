#!/usr/bin/env python3
"""Release consistency checker for oCAS.

Standard library only (Python 3.8+); no third-party imports, no network access,
no file writes. Safe to run locally and in CI on the current working tree.

Source of truth: `version` in `[workspace.package]` of the root `Cargo.toml`
(call it V, e.g. ``0.30.0``). Every user-facing version surface must agree with
V:

  README.md                         status badge ``Beta (V)`` and the
                                    ``ocas = "X.Y"`` dependency snippet
  docs/README_CN.md                 the same two claims (Chinese edition)
  docs/book/en/src/introduction.md  ``currently at version **V**``
  docs/book/zh/src/introduction.md  ``当前版本为 **V**``
  CHANGELOG.md                      first release heading
                                    ``## [V] - YYYY-MM-DD`` with a real date
  README.md / docs/README_CN.md     every ``https://charleshzh.github.io/...``
                                    documentation URL lives under ``/oCAS/``

Extra drift guard (same spirit as the checks above): the release date in the
version tables of both book introductions must equal the CHANGELOG date for V.

Output: exactly one ``PASS: ...`` line and exit status 0 when everything is
consistent; otherwise one ``file:line: 期望 ...，实际 ...`` line per violation
and exit status 1.

Usage:
    python3 ocas-tests/scripts/check_release_consistency.py [--root DIR]
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import sys
from pathlib import Path

# ``0.30.0`` plus optional pre-release suffix (``-rc.1`` etc.).
_VERSION_RE = re.compile(r"\d+\.\d+\.\d+(?:-[0-9A-Za-z.]+)?")
_RELEASE_HEADING_RE = re.compile(
    r"^##\s+\[(?P<version>\d+\.\d+\.\d+)\]\s+-\s+(?P<date>\d{4}-\d{2}-\d{2})\s*$"
)
_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_URL_RE = re.compile(r"https?://charleshzh\.github\.io/(?P<path>[^\s)）\],，、]*)")
# ``| **0.30.0** | 2026-10-24 | ...`` in the book version tables.
_BOOK_ROW_RE = re.compile(
    r"^\|\s*\*\*(?P<version>[0-9][^*]*)\*\*\s*\|\s*(?P<date>\d{4}-\d{2}-\d{2})\s*\|"
)


class Report:
    """Collects violations, formatted as ``file:line: 期望 X，实际 Y``."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.errors: list[str] = []

    def _display(self, path: Path) -> str:
        try:
            return path.relative_to(self.root).as_posix()
        except ValueError:  # pragma: no cover - only for --root misconfiguration
            return path.as_posix()

    def fail(self, path: Path, line: int, expected: str, actual: str) -> None:
        self.errors.append(f"{self._display(path)}:{line}: 期望 {expected}，实际 {actual}")

    @property
    def ok(self) -> bool:
        return not self.errors


def read_lines(path: Path, report: Report) -> list[str]:
    """Read a UTF-8 text file, recording a violation (instead of crashing)."""
    try:
        return path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        report.fail(path, 0, "存在且可读取的文件", f"读取失败（{exc.__class__.__name__}）")
        return []


def find_first(lines: list[str], pattern: re.Pattern[str]):
    """Return ``(line_number, line_text, match)`` for the first match, or ``(0, '', None)``."""
    for line_number, line in enumerate(lines, 1):
        match = pattern.search(line)
        if match is not None:
            return line_number, line.strip(), match
    return 0, "", None


def read_workspace_version(report: Report) -> str | None:
    """Return V from ``[workspace.package] version`` in the root Cargo.toml."""
    path = report.root / "Cargo.toml"
    lines = read_lines(path, report)
    if not lines:
        return None
    in_section = False
    for line_number, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("["):
            in_section = stripped == "[workspace.package]"
            continue
        if not in_section:
            continue
        match = re.match(r'^version\s*=\s*"([^"]+)"', stripped)
        if match is None:
            continue
        version = match.group(1)
        if _VERSION_RE.fullmatch(version) is None:
            report.fail(path, line_number, '合法版本号 "X.Y.Z"', f'"{version}"')
            return None
        return version
    report.fail(
        path,
        0,
        '[workspace.package] 下的 version = "X.Y.Z"',
        "未找到",
    )
    return None


def check_status_badge(report: Report, version: str) -> None:
    """README/README_CN status lines must read ``Beta (V)``."""
    expected = f"Beta ({version})"
    pattern = re.compile(r"Beta \((?P<version>[0-9][^)]*)\)")
    for rel_path in ("README.md", "docs/README_CN.md"):
        path = report.root / rel_path
        lines = read_lines(path, report)
        if not lines:
            continue
        line_number, _, match = find_first(lines, pattern)
        if match is None:
            report.fail(path, 0, expected, "未找到匹配行")
        elif match.group("version") != version:
            report.fail(path, line_number, expected, f"Beta ({match.group('version')})")


def check_install_snippet(report: Report, version: str) -> None:
    """README/README_CN quick-start snippets must read ``ocas = "X.Y"``."""
    expected_version = ".".join(version.split(".")[:2])
    expected = f'ocas = "{expected_version}"'
    pattern = re.compile(r'^(?P<dep>ocas)\s*=\s*"(?P<version>[^"]+)"\s*$')
    for rel_path in ("README.md", "docs/README_CN.md"):
        path = report.root / rel_path
        lines = read_lines(path, report)
        if not lines:
            continue
        line_number, _, match = find_first(lines, pattern)
        if match is None:
            report.fail(path, 0, expected, "未找到匹配行")
        elif match.group("version") != expected_version:
            report.fail(
                path,
                line_number,
                expected,
                f'ocas = "{match.group("version")}"',
            )


def check_book_version_statement(report: Report, version: str) -> None:
    """Book introductions must state the current version (en + zh)."""
    checks = (
        ("docs/book/en/src/introduction.md", re.compile(r"currently at version \*\*(?P<version>[0-9][^*]*)\*\*")),
        ("docs/book/zh/src/introduction.md", re.compile(r"当前版本为 \*\*(?P<version>[0-9][^*]*)\*\*")),
    )
    for rel_path, pattern in checks:
        path = report.root / rel_path
        lines = read_lines(path, report)
        if not lines:
            continue
        line_number, _, match = find_first(lines, pattern)
        expected = f"**{version}**"
        if match is None:
            report.fail(path, 0, expected, "未找到匹配行")
        elif match.group("version") != version:
            report.fail(path, line_number, expected, f"**{match.group('version')}**")


def check_agents_version(report: Report, version: str) -> None:
    """AGENTS.md's cross-harness version line must match the workspace version."""
    path = report.root / "AGENTS.md"
    lines = read_lines(path, report)
    if not lines:
        return
    pattern = re.compile(r"Current version:\s*\*\*(?P<version>[0-9][^*]*)\*\*")
    line_number, _, match = find_first(lines, pattern)
    expected = f"Current version: **{version}**"
    if match is None:
        report.fail(path, 0, expected, "未找到匹配行")
    elif match.group("version") != version:
        report.fail(path, line_number, expected, f"Current version: **{match.group('version')}**")

def check_changelog(report: Report, version: str) -> str | None:
    """Check the CHANGELOG's first release heading; return its date (or None)."""
    path = report.root / "CHANGELOG.md"
    lines = read_lines(path, report)
    if not lines:
        return None
    heading = None
    for line_number, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped.startswith("## ["):
            continue
        # A future ``## [Unreleased]`` heading must not hide the first release.
        if "unreleased" in stripped.lower():
            continue
        heading = (line_number, _RELEASE_HEADING_RE.match(stripped), stripped)
        break
    expected = f"## [{version}] - YYYY-MM-DD"
    if heading is None:
        report.fail(path, 0, expected, "未找到任何发布标题行")
        return None
    line_number, match, stripped = heading
    if match is None:
        report.fail(path, line_number, expected, stripped)
        return None
    if match.group("version") != version:
        report.fail(path, line_number, expected, stripped)
        return None
    date = match.group("date")
    if _DATE_RE.match(date) is None:
        report.fail(path, line_number, "合法日期 YYYY-MM-DD", date)
        return None
    try:
        _dt.date.fromisoformat(date)
    except ValueError:
        report.fail(path, line_number, "合法日期 YYYY-MM-DD", date)
        return None
    return date


def check_book_table_dates(report: Report, version: str, changelog_date: str | None) -> None:
    """Book version-table row for V must carry the CHANGELOG release date."""
    if changelog_date is None:
        return
    for rel_path in ("docs/book/en/src/introduction.md", "docs/book/zh/src/introduction.md"):
        path = report.root / rel_path
        lines = read_lines(path, report)
        if not lines:
            continue
        row = None
        for line_number, line in enumerate(lines, 1):
            match = _BOOK_ROW_RE.match(line.strip())
            if match is not None and match.group("version") == version:
                row = (line_number, match)
                break
        expected = f"| **{version}** | {changelog_date} |"
        if row is None:
            report.fail(path, 0, expected, f"未找到 **{version}** 版本行")
        elif row[1].group("date") != changelog_date:
            report.fail(
                path,
                row[0],
                expected,
                f"| **{row[1].group('version')}** | {row[1].group('date')} |",
            )


def check_docs_site_urls(report: Report) -> None:
    """Every GitHub Pages URL must use the canonical ``/oCAS/`` path segment."""
    for rel_path in ("README.md", "docs/README_CN.md"):
        path = report.root / rel_path
        lines = read_lines(path, report)
        for line_number, line in enumerate(lines, 1):
            for match in _URL_RE.finditer(line):
                url_path = match.group("path")
                if url_path.split("/", 1)[0] == "oCAS":
                    continue
                report.fail(
                    path,
                    line_number,
                    "文档站 URL 路径包含 /oCAS/",
                    f"https://charleshzh.github.io/{url_path}",
                )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check that oCAS version and documentation surfaces are consistent."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="repository root (default: two levels above this script)",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve() if args.root is not None else Path(__file__).resolve().parents[2]

    report = Report(root)
    version = read_workspace_version(report)
    if version is not None:
        check_status_badge(report, version)
        check_install_snippet(report, version)
        check_book_version_statement(report, version)
        check_agents_version(report, version)
        changelog_date = check_changelog(report, version)
        check_book_table_dates(report, version, changelog_date)
        check_docs_site_urls(report)

    if not report.ok:
        for error in report.errors:
            print(error, file=sys.stderr)
        print(f"FAILED: {len(report.errors)} consistency violation(s)", file=sys.stderr)
        return 1

    print(f"PASS: release consistency checks passed (workspace version {version})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
