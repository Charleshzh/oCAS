import argparse
import subprocess
import sys
from pathlib import Path


_RUNNER_DIR = Path(__file__).resolve().parent / "symbolica_runner"


def _strip_banner(stdout: str) -> str:
    """Remove Symbolica's restricted-license banner.

    Without a license key the Symbolica library prints a boxed banner to
    stdout at startup; it is not part of the task result.
    """
    start = stdout.find("┌")
    if start != -1:
        end = stdout.find("└", start)
        if end != -1:
            nl = stdout.find("\n", end)
            if nl != -1:
                return stdout[nl + 1 :]
    return stdout


def run_symbolica(task: str, expr: str) -> str:
    """Run the isolated Symbolica runner and return normalized stdout."""
    args = ["cargo", "run", "--quiet", "--manifest-path", str(_RUNNER_DIR / "Cargo.toml"), "--", task, expr]
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"symbolica runner failed: {result.stderr}")
    return _strip_banner(result.stdout).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare oCAS results against Symbolica.")
    parser.add_argument("task", help="Operation to compare (parse, diff, expand, factor, series, simplify)")
    parser.add_argument("expr", help="Expression string, using ^ for exponentiation")
    args = parser.parse_args()
    print(run_symbolica(args.task, args.expr))
    return 0


if __name__ == "__main__":
    sys.exit(main())
