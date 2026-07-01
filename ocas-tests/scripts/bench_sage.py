#!/usr/bin/env sage
"""SageMath performance-comparison harness for oCAS.

Run manually (SageMath must be installed):

    sage scripts/bench_sage.py <task> <expr> <iters>

Tasks: ``parse``, ``diff``, ``expand``, ``factor``, ``gcd``. For ``gcd`` pass
two polynomials separated by ``;``. The script prints the total elapsed time
in nanoseconds (one warmup run, then `iters` measured runs), matching the
output contract of ``bench_sympy.py`` so the two can be compared directly.

Note: SageMath uses ``^`` for exponentiation, identical to oCAS, so no
translation is needed (unlike the SymPy harness).
"""

import sys
import timeit

import sage.all as sage  # noqa: E402  (sage runtime)


def make_stmt(task, expr_str):
    x, y, z = sage.var("x y z")
    locs = {"x": x, "y": y, "z": z}

    def parse(s):
        return sage.SR(s)

    if task == "parse":
        return lambda: parse(expr_str)
    if task == "diff":
        return lambda: sage.diff(parse(expr_str), x)
    if task == "expand":
        return lambda: sage.expand(parse(expr_str))
    if task == "factor":
        return lambda: sage.factor(parse(expr_str))
    if task == "gcd":
        parts = expr_str.split(";")
        a = sage.expand(parse(parts[0]))
        b = sage.expand(parse(parts[1])) if len(parts) > 1 else a + 1
        # SageMath's polynomial gcd needs a polynomial ring; SR gcd works for
        # symbolic expressions via .gcd().
        return lambda: a.gcd(b)
    raise ValueError(f"unknown task: {task}")


def main():
    if len(sys.argv) != 4:
        print("usage: sage bench_sage.py <task> <expr> <iters>", file=sys.stderr)
        sys.exit(2)
    task = sys.argv[1]
    expr_str = sys.argv[2]
    iters = int(sys.argv[3])

    stmt = make_stmt(task, expr_str)
    stmt()  # warmup

    total_seconds = timeit.timeit(stmt, number=iters)
    print(int(total_seconds * 1e9))


if __name__ == "__main__":
    main()
