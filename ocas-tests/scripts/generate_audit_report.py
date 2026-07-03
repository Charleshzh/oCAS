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


def generate_report(simple_summary: dict, complex_summary: dict, failures: list[str]) -> str:
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
## Notes

- Simple and medium tests are expected to run in CI and pass automatically.
- Complex and very complex tests are marked `#[ignore]` and run only during
  manual audits or via this report generator.
- SymPy comparison tests are skipped if `uv` is not available or no venv with
  SymPy is configured.
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

    report = generate_report(simple_summary, complex_summary, failures)

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
