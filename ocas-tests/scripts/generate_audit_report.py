import argparse
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = ROOT / "docs" / "planning" / "correctness"


def run_cargo_tests(args: list[str]) -> tuple[int, str, str]:
    cmd = ["cargo", "test", "-p", "ocas-tests", "correctness", "--", "--nocapture", *args]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def parse_test_summary(output: str) -> dict:
    lines = output.splitlines()
    summary = {
        "passed": 0,
        "failed": 0,
        "ignored": 0,
        "filtered": 0,
        "measured": 0,
        "time": "",
    }
    for line in lines:
        if line.startswith("test result:"):
            parts = line.split(";")
            for part in parts:
                part = part.strip()
                if "passed" in part:
                    summary["passed"] = int(part.split()[0])
                elif "failed" in part:
                    summary["failed"] = int(part.split()[0])
                elif "ignored" in part:
                    summary["ignored"] = int(part.split()[0])
                elif "filtered out" in part:
                    summary["filtered"] = int(part.split()[0])
                elif "measured" in part:
                    summary["measured"] = int(part.split()[0])
                elif "finished" in part:
                    summary["time"] = part
    return summary


def collect_failures(output: str) -> list[str]:
    failures = []
    for line in output.splitlines():
        if "FAILED" in line or "failures:" in line:
            failures.append(line)
    return failures


_RUNNER_DIR = Path(__file__).resolve().parent / "symbolica_runner"


def run_symbolica_factor_time(expr: str) -> float | None:
    """Time Symbolica `factor` over the integers; returns ns/op or None."""
    args = [
        "cargo",
        "run",
        "--quiet",
        "--release",
        "--manifest-path",
        str(_RUNNER_DIR / "Cargo.toml"),
        "--",
        "factor_time",
        expr,
    ]
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return None
    try:
        return float(result.stdout.strip())
    except ValueError:
        return None


_CRITERION_TIME = __import__("re").compile(r"time:\s+\[[\d.]+ (\w+)\s+[\d.]+ \w+\s+[\d.]+ \w+\]")


def run_ocas_bench(filter_str: str) -> float | None:
    """Run one criterion benchmark and return the median estimate in ns."""
    args = [
        "cargo",
        "bench",
        "-p",
        "ocas-tests",
        "--bench",
        "poly_factor",
        "--",
        "--warm-up-time",
        "1",
        "--measurement-time",
        "4",
        filter_str,
    ]
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return None
    import re

    m = re.search(r"time:\s+\[[\d.]+ \w+ ([\d.]+) (\w+) [\d.]+ \w+\]", result.stdout)
    if not m:
        return None
    value, unit = float(m.group(1)), m.group(2)
    scale = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "µs": 1e3, "ms": 1e6, "s": 1e9}
    return value * scale.get(unit, 1.0)


def sparse_4var_expr() -> str:
    """Sparse 4-variable nonconstant-LC benchmark input (factored form)."""

    def dedup(terms: list[tuple[list[int], int]]) -> list[tuple[list[int], int]]:
        merged: dict[tuple[int, ...], int] = {}
        for exp, c in terms:
            merged[tuple(exp)] = c  # last wins, matching set_term semantics
        return [(list(e), c) for e, c in merged.items()]

    f1: list[tuple[list[int], int]] = [([2, 1, 1, 0], 1)]
    f2: list[tuple[list[int], int]] = [([1, 1, 0, 0], 1), ([1, 0, 0, 1], 1)]
    for i in range(4):
        for j in range(3):
            c1 = (i * 7 + j * 3) % 4 + 1
            c2 = (i * 5 + j * 11 + 2) % 4 + 1
            f1.append(([i % 2, i, j, (i + j) % 2], c1))
            f2.append(([0, (i + 1) % 3, (j + 2) % 2, i % 3], c2))

    def to_str(terms: list[tuple[list[int], int]]) -> str:
        parts = []
        for exp, c in dedup(terms):
            mon = [f"{v}^{e}" if e > 1 else v for v, e in zip("xyzw", exp) if e]
            parts.append("*".join([str(c), *mon]) if mon else str(c))
        return "(" + "+".join(parts) + ")"

    return f"{to_str(f1)}*{to_str(f2)}"


def factorization_comparison() -> list[tuple[str, float | None, float | None]]:
    """(case, oCAS ns/op, Symbolica ns/op) for same-scale factorization."""
    trivariate_expr = "((z*x^2+y)*(x+1))"
    sparse_expr = sparse_4var_expr()
    return [
        (
            "trivariate_nonconstant_lcoeff",
            run_ocas_bench("trivariate_nonconstant_lcoeff"),
            run_symbolica_factor_time(trivariate_expr),
        ),
        (
            "sparse_4var_nonconstant_lcoeff",
            run_ocas_bench("sparse_4var_nonconstant_lcoeff"),
            run_symbolica_factor_time(sparse_expr),
        ),
    ]


def generate_report(
    simple_summary: dict,
    complex_summary: dict,
    failures: list[str],
    comparison: list[tuple[str, float | None, float | None]],
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    report = f"""# oCAS Correctness Audit Report

Generated: {now}

## Summary

| Category | Passed | Failed | Ignored | Filtered |
|---|---:|---:|---:|---:|
| Simple + Medium | {simple_summary['passed']} | {simple_summary['failed']} | {simple_summary['ignored']} | {simple_summary['filtered']} |
| Complex + Very Complex | {complex_summary['passed']} | {complex_summary['failed']} | {complex_summary['ignored']} | {complex_summary['filtered']} |

## Failures

"""
    if failures:
        for failure in failures:
            report += f"- {failure}\n"
    else:
        report += "No failures detected.\n"

    report += """
## Symbolica Factorization Comparison

Same-scale multivariate factorization timings (ns/op; `—` = unavailable).

| Case | oCAS | Symbolica | oCAS / Symbolica |
|---|---:|---:|---:|
"""
    for case, ocas_ns, sym_ns in comparison:
        def fmt(ns: float | None) -> str:
            return f"{ns:,.0f}" if ns is not None else "—"

        ratio = f"{ocas_ns / sym_ns:.2f}" if ocas_ns and sym_ns else "—"
        report += f"| {case} | {fmt(ocas_ns)} | {fmt(sym_ns)} | {ratio} |\n"

    report += """
## Notes

- Simple and medium tests are expected to run in CI and pass automatically.
- Complex and very complex tests are marked `#[ignore]` and run only during
  manual audits or via this report generator.
- SymPy comparison tests are skipped if `uv` is not available or no venv with
  SymPy is configured.
- Symbolica timings come from the isolated `symbolica_runner` crate
  (`factor_time`, 20 iterations); oCAS timings from the `poly_factor`
  criterion group.
"""
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate oCAS correctness audit report")
    parser.add_argument(
        "--output",
        type=Path,
        help="Output path for the report (default: docs/planning/correctness/audit-<date>.md)",
    )
    args = parser.parse_args()

    print("Running simple + medium tests...")
    _, simple_out, simple_err = run_cargo_tests([])
    simple_summary = parse_test_summary(simple_out + simple_err)
    print(f"Simple + medium: {simple_summary}")

    print("Running complex + very complex tests...")
    _, complex_out, complex_err = run_cargo_tests(["--ignored"])
    complex_summary = parse_test_summary(complex_out + complex_err)
    print(f"Complex + very complex: {complex_summary}")

    failures = collect_failures(simple_out + simple_err + complex_out + complex_err)

    print("Running Symbolica factorization comparison...")
    comparison = factorization_comparison()
    print(f"Comparison: {comparison}")

    report = generate_report(simple_summary, complex_summary, failures, comparison)

    output_path = args.output
    if output_path is None:
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
        date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        output_path = REPORT_DIR / f"audit-{date}.md"

    output_path.write_text(report, encoding="utf-8")
    print(f"Report written to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
