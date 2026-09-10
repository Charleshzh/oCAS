#!/usr/bin/env python3
"""Diff two integrate_1892 failure dumps (JSONL) at bucket granularity.

Usage:
    uv run python ocas-tests/scripts/diff_1892_failures.py BEFORE.jsonl AFTER.jsonl

Each input line is a JSON record produced by benches/integrate_1892.rs:
    {"id", "bucket", "outcome", "ms", "var", "integrand"}
Cases absent from a dump were solved in that run. The report lists per-bucket
newly-solved / regressed counts, timeout and crash deltas, and the integrands
that changed state so regressions are caught before merge.
"""

import json
import sys
from collections import defaultdict

OUTCOMES = ("fallback", "timeout", "crashed", "parse_err")


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


def main() -> int:
    if len(sys.argv) != 3:
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

    if regressed:
        print("\nREGRESSIONS (must be zero before merge):")
        for rec in sorted(regressed.values(), key=lambda r: r["id"]):
            print(f"  [{rec['bucket']}] {rec['id']}: {rec['integrand']} ({rec['outcome']})")
    print("\nnewly solved sample (up to 40):")
    for rec in sorted(newly_solved.values(), key=lambda r: r["id"])[:40]:
        print(f"  [{rec['bucket']}] {rec['id']}: {rec['integrand']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
