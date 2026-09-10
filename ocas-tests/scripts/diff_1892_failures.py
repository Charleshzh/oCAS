#!/usr/bin/env python3
"""Diff two integrate_1892 failure dumps (JSONL) at bucket granularity.

Usage:
    uv run python ocas-tests/scripts/diff_1892_failures.py BEFORE.jsonl AFTER.jsonl
    uv run python ocas-tests/scripts/diff_1892_failures.py BEFORE.jsonl AFTER.jsonl \
        [integrate_1892_unverified.jsonl]

Each input line is a JSON record produced by benches/integrate_1892.rs:
    {"id", "bucket", "outcome", "ms", "var", "integrand"}
Cases absent from a dump were solved in that run. The report lists per-bucket
newly-solved / regressed counts, timeout and crash deltas, and the integrands
that changed state so regressions are caught before merge.

The optional third argument is the numerical-verification dump
(`integrate_1892_unverified.jsonl`), whose records are
    {"id", "bucket", "class", "detail", "integrand", "result"}
with `class` in {"mismatch", "indeterminate"}. When given, the report adds the
wrong-answer count (`mismatch`) and the inconclusive count; mismatches are
printed in full because each one is a wrong antiderivative. When omitted, the
script looks for the file next to the AFTER dump and for the sibling
`integrate_1892_report.json`, and degrades gracefully if either is missing.
"""

import json
import os
import sys
from collections import defaultdict

OUTCOMES = ("fallback", "timeout", "crashed", "parse_err")
VERIFY_CLASSES = ("mismatch", "indeterminate")


def load(path: str) -> dict[str, dict]:
    cases: dict[str, dict] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            cases[rec["id"]] = rec
    return cases


def load_json(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def resolve_aux(after_path: str, explicit: str | None) -> tuple[str, str]:
    """Return (unverified_dump, report_json) paths, preferring explicit args."""
    directory = os.path.dirname(os.path.abspath(after_path))
    unverified = explicit or os.path.join(directory, "integrate_1892_unverified.jsonl")
    report = os.path.join(directory, "integrate_1892_report.json")
    return unverified, report


def report_verification(after: dict[str, dict], unverified_path: str, report_path: str) -> int:
    """Print verification stats; return the mismatch count."""
    stats = load_json(report_path)
    if stats:
        print("verification (from report json):")
        for key in (
            "solved",
            "verified_solved",
            "unverified_solved",
            "verify_indeterminate",
            "verify_mismatches",
            "verify_worst_rel",
            "timed_out",
            "crashed",
            "total_ms",
            "coverage_pct",
        ):
            if key in stats:
                print(f"  {key:>20}: {stats[key]}")

    classes: dict[str, dict[str, dict]] = defaultdict(dict)
    if os.path.exists(unverified_path):
        for case_id, rec in load(unverified_path).items():
            classes[rec.get("class", "unknown")][case_id] = rec
    elif not stats:
        print("verification: no dump and no report json found (string coverage only)")
        return 0

    for cls in VERIFY_CLASSES:
        recs = classes.get(cls, {})
        print(f"  {cls:>20}: {len(recs)}")

    mismatches = classes.get("mismatch", {})
    if mismatches:
        print("\nWRONG ANSWERS (numerical verification failed):")
        for rec in sorted(mismatches.values(), key=lambda r: r["id"]):
            print(f"  [{rec['bucket']}] {rec['id']}: {rec['integrand']}")
            print(f"      result: {rec.get('result', '')}")
            print(f"      detail: {rec.get('detail', '')}")
    unsolved = len(after)
    if unsolved and not mismatches:
        print(f"  ({unsolved} cases remain unsolved; none of the solved set is wrong)")
    return len(mismatches)


def main() -> int:
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        return 2
    before = load(sys.argv[1])
    after = load(sys.argv[2])

    newly_solved = {i: before[i] for i in before.keys() - after.keys()}
    regressed = {i: after[i] for i in after.keys() - before.keys()}

    by_bucket: dict[str, dict[str, int]] = defaultdict(
        lambda: {"newly_solved": 0, "regressed": 0}
    )
    for rec in newly_solved.values():
        by_bucket[rec["bucket"]]["newly_solved"] += 1
    for rec in regressed.values():
        by_bucket[rec["bucket"]]["regressed"] += 1

    def outcome_counts(dump: dict[str, dict]) -> dict[str, int]:
        c: dict[str, int] = defaultdict(int)
        for rec in dump.values():
            c[rec["outcome"]] += 1
        return c

    oc_before, oc_after = outcome_counts(before), outcome_counts(after)

    print(f"newly solved: {len(newly_solved)}  regressed: {len(regressed)}")
    print(f"{'bucket':>20}  {'+solved':>7}  {'-regr':>5}")
    for bucket in sorted(by_bucket):
        d = by_bucket[bucket]
        print(f"{bucket:>20}  {d['newly_solved']:>7}  {d['regressed']:>5}")
    print("outcomes (before -> after):")
    for outcome in OUTCOMES:
        print(f"  {outcome:>10}: {oc_before.get(outcome, 0)} -> {oc_after.get(outcome, 0)}")

    unverified_path, report_path = resolve_aux(
        sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None
    )
    mismatches = report_verification(after, unverified_path, report_path)

    if regressed:
        print("\nREGRESSIONS (must be zero before merge):")
        for rec in sorted(regressed.values(), key=lambda r: r["id"]):
            print(f"  [{rec['bucket']}] {rec['id']}: {rec['integrand']} ({rec['outcome']})")
    print("\nnewly solved sample (up to 40):")
    for rec in sorted(newly_solved.values(), key=lambda r: r["id"])[:40]:
        print(f"  [{rec['bucket']}] {rec['id']}: {rec['integrand']}")

    return 1 if regressed or mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
