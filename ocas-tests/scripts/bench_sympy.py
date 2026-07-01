import sys
import timeit
import sympy as sp


def make_stmt(task: str, expr_str: str):
    x, y, z = sp.symbols("x y z")
    # oCAS uses ^ for exponentiation; SymPy uses **.
    sympy_expr = expr_str.replace("^", "**")
    if task == "parse":
        return lambda: sp.parse_expr(sympy_expr)
    if task == "diff":
        return lambda: sp.diff(sp.parse_expr(sympy_expr), x)
    if task == "expand":
        return lambda: sp.expand(sp.parse_expr(sympy_expr))
    raise ValueError(f"unknown task: {task}")


def main() -> None:
    task = sys.argv[1]
    expr_str = sys.argv[2]
    iters = int(sys.argv[3])

    stmt = make_stmt(task, expr_str)
    # single warmup run to discount import / compilation caches
    stmt()

    total_seconds = timeit.timeit(stmt, number=iters)
    print(int(total_seconds * 1e9))


if __name__ == "__main__":
    main()
