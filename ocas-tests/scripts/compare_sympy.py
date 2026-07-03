import argparse
import json
import sys
import timeit
import sympy as sp
from typing import Callable


_X, _Y, _Z = sp.symbols("x y z")


def _to_sympy(expr_str: str) -> str:
    """Convert oCAS ^ exponentiation to SymPy **."""
    return expr_str.replace("^", "**")


def _parse_expr(expr_str: str):
    return sp.parse_expr(_to_sympy(expr_str))


def _normalize(expr) -> str:
    """Return a canonical string form for comparison.

    SymPy does not guarantee term order, so we expand and collect on all free
    symbols before printing. This is sufficient for comparing polynomial /
    rational / elementary expressions without building a full CAS equivalence
    checker.
    """
    expanded = sp.expand(expr)
    symbols = list(expanded.free_symbols)
    if symbols:
        expanded = sp.collect(expanded, symbols)
    return str(expanded)


def _make_task(task: str, expr_str: str) -> Callable:
    x, y, z = _X, _Y, _Z

    if task == "series":
        # expr_str may encode the order as "expr:order" (default 10).
        if ":" in expr_str:
            expr_part, order_part = expr_str.rsplit(":", 1)
            order = int(order_part)
        else:
            expr_part = expr_str
            order = 10
        sympy_expr = _parse_expr(expr_part)
        return lambda: sp.series(sympy_expr, x, 0, order + 1).removeO()

    sympy_expr = _parse_expr(expr_str)

    if task == "parse":
        return lambda: _parse_expr(expr_str)
    if task == "diff":
        return lambda: sp.diff(sympy_expr, x)
    if task == "expand":
        return lambda: sp.expand(sympy_expr)
    if task == "factor":
        return lambda: sp.factor(sympy_expr)
    if task == "gcd":
        parts = expr_str.split(";")
        a = sp.expand(_parse_expr(parts[0]))
        b = sp.expand(_parse_expr(parts[1])) if len(parts) > 1 else a + 1
        return lambda: sp.gcd(a, b)
    if task == "simplify":
        return lambda: sp.simplify(sympy_expr)
    if task == "integrate":
        return lambda: sp.integrate(sympy_expr, x)
    if task == "eval":
        # Evaluate numerically at the point encoded in expr_str as "expr @ x=1,y=2".
        if "@" not in expr_str:
            return lambda: sympy_expr.evalf()
        expr_part, subs_part = expr_str.split("@", 1)
        expr = _parse_expr(expr_part)
        subs = {}
        for assignment in subs_part.split(","):
            var, val = assignment.split("=")
            subs[sp.symbols(var.strip())] = sp.sympify(val.strip())
        return lambda: expr.evalf(subs=subs)
    if task == "solve_linear":
        # expr_str is a semicolon-separated list of equations in x,y,z.
        equations = [sp.sympify(_to_sympy(eq)) for eq in expr_str.split(";")]
        symbols = [s for s in [x, y, z] if any(s in eq.free_symbols for eq in equations)]
        return lambda: sp.linsolve(equations, symbols)
    if task == "roots":
        return lambda: sp.nroots(sympy_expr)
    raise ValueError(f"unknown task: {task}")


def time_task(task: str, expr_str: str, iters: int) -> int:
    """Run the task `iters` times and return total nanoseconds."""
    stmt = _make_task(task, expr_str)
    stmt()  # warmup
    total_seconds = timeit.timeit(stmt, number=iters)
    return int(total_seconds * 1e9)


def compute_task(task: str, expr_str: str) -> str:
    """Run the task once and return a normalized string result."""
    stmt = _make_task(task, expr_str)
    result = stmt()
    return _normalize(result)


def _check_equivalent(ocas_expr: str, sympy_expr: str, task: str, input_expr: str) -> bool:
    """Check whether the oCAS result is equivalent to the SymPy reference.

    The strategy depends on the operation: for integration, the constant of
    integration is irrelevant, so we differentiate the oCAS result and compare
    with the original integrand; for factorization, both sides are checked by
    expansion; for most other cases a simple difference test is sufficient.
    """
    a = _parse_expr(ocas_expr)
    b = _parse_expr(sympy_expr)

    if task == "integrate":
        # Differentiate the oCAS antiderivative and compare with the original integrand.
        # The constant of integration is irrelevant.
        integrand = _parse_expr(input_expr)
        return sp.simplify(sp.diff(a, _X) - integrand) == 0
    if task == "factor":
        return sp.expand(a - b) == 0
    if task == "roots":
        # Compare numeric roots (complex) sorted by real then imaginary parts.
        if not (hasattr(a, "__len__") and hasattr(b, "__len__")):
            return False
        try:
            a_roots = sorted([complex(r) for r in a], key=lambda c: (c.real, c.imag))
            b_roots = sorted([complex(r) for r in b], key=lambda c: (c.real, c.imag))
            if len(a_roots) != len(b_roots):
                return False
            return all(abs(ar - br) < 1e-10 for ar, br in zip(a_roots, b_roots))
        except Exception:
            return False

    return sp.simplify(a - b) == 0


def verify_task(task: str, input_expr: str, ocas_result: str) -> bool:
    """Verify that `ocas_result` matches the SymPy reference for `task(input_expr)`."""
    reference = compute_task(task, input_expr)
    return _check_equivalent(ocas_result, reference, task, input_expr)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare oCAS results against SymPy.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check_parser = subparsers.add_parser("check", help="Print normalized reference result")
    check_parser.add_argument("task", help="Operation to compare")
    check_parser.add_argument("expr", help="Expression string, using ^ for exponentiation")

    verify_parser = subparsers.add_parser(
        "verify", help="Verify an oCAS result against SymPy (reads JSON from stdin)"
    )
    verify_parser.add_argument("task", help="Operation to compare")
    verify_parser.add_argument("expr", help="Input expression string")
    verify_parser.add_argument(
        "ocas_result", nargs="?", help="oCAS result string (or read from stdin)"
    )

    time_parser = subparsers.add_parser("time", help="Time the operation")
    time_parser.add_argument("task", help="Operation to time")
    time_parser.add_argument("expr", help="Expression string")
    time_parser.add_argument("iters", type=int, help="Number of iterations")

    args = parser.parse_args()

    if args.command == "check":
        print(compute_task(args.task, args.expr))
    elif args.command == "verify":
        if args.ocas_result is not None:
            ocas_result = args.ocas_result
        else:
            data = json.load(sys.stdin)
            ocas_result = data["ocas_result"]
        reference = compute_task(args.task, args.expr)
        ok = _check_equivalent(ocas_result, reference, args.task, args.expr)
        print("true" if ok else "false")
        return 0 if ok else 1
    elif args.command == "time":
        print(time_task(args.task, args.expr, args.iters))
    return 0


if __name__ == "__main__":
    sys.exit(main())
