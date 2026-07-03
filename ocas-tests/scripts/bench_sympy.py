#!/usr/bin/env python3
"""Thin wrapper for the existing benchmark that calls compare_sympy.py.

The new compare_sympy.py supports both timing and correctness modes. This file
preserves the original CLI used by benches/sympy_comparison.rs.
"""
import sys
import subprocess


def main() -> None:
    task = sys.argv[1]
    expr_str = sys.argv[2]
    iters = sys.argv[3]
    result = subprocess.run(
        ["uv", "run", "python", "scripts/compare_sympy.py", "time", task, expr_str, iters],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    print(result.stdout.strip())


if __name__ == "__main__":
    main()
