#!/usr/bin/env python3
"""0.31 attribution: which 0.30.0 unsolved Rubi cases need dilogarithms (Li2)?

Reads the deterministic answer-augmented corpus sample
(`ocas-tests/data/rubi_1892_answers.tsv`, produced by
`fetch_rubi_corpus.py --answers-out`) plus the 0.30.0 failure dump
(`ocas-tests/data/failures_final_030.jsonl`, or a newer `--failures` file)
and writes `ocas-tests/data/polylog_attribution_031.csv` for the unsolved
cases whose integrand or reference answer mentions `polylog`:

* class **A**: the *integrand* contains `polylog` (Li2 is already in the
  integrand; needs the polylog head + derivative + by-parts);
* class **B**: the integrand is elementary but the *reference answer*
  contains `polylog` (needs a reduction to dilogarithms);
* an unsolved case with neither signature is not written to the CSV at
  all (its class would be empty); it is only counted in the summary.

Columns (properly quoted with the `csv` module):

    id, bucket, integrand, answer_polylog_orders, class, fragment_hint, reason

* `bucket` replicates the nine-way classification of
  `ocas-tests/benches/integrate_1892.rs` exactly (and is cross-checked
  against the `bucket` field that the bench stored in the failures JSONL).
  The replication includes the bench's `contains_fractional_pow` quirk: a
  `Pow` node is only inspected for a non-`Num` exponent at the node itself,
  it does not recurse into the base or the exponent.  With that quirk the
  replica agrees with all 1502 stored bench buckets, so a drift here means
  the bench changed.
* `answer_polylog_orders`: sorted unique integer first arguments of the
  `polylog(<s>, ...)` calls in the answer, `;`-joined (`2`, `2;3`, ...);
  empty when the answer only carries a symbolic order (e.g.
  `polylog(k, ...)`).  Only the fourth top-level answer field parsed by the
  fetch script is inspected (the corpus' alternative fifth-field answers
  are not merged).
* `fragment_hint` (class B only, empty for class A): how close the case is
  to the 0.31a log-rational/dilogarithm reducer fragment
  (`docs/planning/POLYLOG_SPEC_CN.md` section 4):
  - `candidate`: the answer needs only order 2 and the integrand looks like
    `R(x)*log(u(x))` with a *reducible* log argument `u` (affine in the
    variable, a Mobius transform affine/affine, or a constant power of an
    affine argument) and a rational weight whose variable-carrying
    denominator factors are all *linear* (first-order poles);
  - `non-affine`: order 2 only, but the shape misses the fragment above
    (non-affine/Mobius log argument, multiple or no `log` kernels,
    higher-order denominators, ...); the `reason` column says which check
    failed;
  - `order3+`: the answer needs order 3 (or a mixture whose highest order
    is 3) — out of 0.31a, needs 0.31b;
  - `higher-order`: the answer needs an order of 4 or more.
* `reason`: short deterministic evidence for the classification (answer
  orders, log-argument shape, denominator shape, or the A-class note).

The heuristic is text/tree-level and explicitly *not* a proof: it is an
audit to pick targets, not an acceptance gate.

Inputs are gitignored benchmark data (`ocas-tests/data/`).  When the
answers file is missing the script prints a clear message and exits 0
(graceful skip, matching the bench's missing-corpus convention); with
`--failures` missing it does the same.  The script is deterministic,
idempotent (it fully recomputes the CSV) and Python 3.11 standard library
only.

Usage (all defaults are repo-relative):

    uv run python ocas-tests/scripts/attribute_031.py
    uv run python ocas-tests/scripts/attribute_031.py \
        ocas-tests/data/rubi_1892_answers.tsv \
        --failures ocas-tests/data/failures_final_030.jsonl
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ANSWERS = Path("ocas-tests/data/rubi_1892_answers.tsv")
DEFAULT_FAILURES = Path("ocas-tests/data/failures_final_030.jsonl")
DEFAULT_OUT = Path("ocas-tests/data/polylog_attribution_031.csv")

CSV_HEADER = [
    "id",
    "bucket",
    "integrand",
    "answer_polylog_orders",
    "class",
    "fragment_hint",
    "reason",
]

# ---------------------------------------------------------------------------
# Bucket classification: exact port of `bucket()` from
# `ocas-tests/benches/integrate_1892.rs`.
# ---------------------------------------------------------------------------

TRIG_HEADS = {"sin", "cos", "tan", "cot", "sec", "csc"}
HYPERBOLIC_HEADS = {"sinh", "cosh", "tanh", "coth", "sech", "csch"}
INVERSE_HEADS = {
    "asin",
    "acos",
    "atan",
    "acot",
    "asec",
    "acsc",
    "asinh",
    "acosh",
    "atanh",
    "acoth",
    "asech",
    "acsch",
}
SPECIAL_HEADS = {
    "erf",
    "erfc",
    "erfi",
    "Ei",
    "Si",
    "Ci",
    "Shi",
    "Chi",
    "fresnels",
    "fresnelc",
    "EllipticF",
    "EllipticE",
    "EllipticPi",
}


class _Node:
    """Mini expression tree, only as much as the bench classification needs."""

    __slots__ = ("kind", "a", "b", "name", "args")

    def __init__(self, kind, a=None, b=None, name=None, args=None):
        self.kind = kind
        self.a = a
        self.b = b
        self.name = name
        self.args = args or []


_TOKEN_RE = re.compile(r"\s*(?:(\d+)|([A-Za-z_][A-Za-z0-9_]*)|(.))", re.S)


def _tokenize(text: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    pos = 0
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if m is None:
            break
        pos = m.end()
        if m.group(1) is not None:
            tokens.append(("num", m.group(1)))
        elif m.group(2) is not None:
            tokens.append(("id", m.group(2)))
        else:
            tokens.append((m.group(3), m.group(3)))
    return tokens


class _Parser:
    """Recursive-descent parser for the subset of oCAS parse syntax the
    converted Rubi integrands use (integer literals, identifiers, calls,
    + - * / ^, parentheses).  Mirrors `ocas-parse/src/parser.rs` precedence."""

    def __init__(self, text: str):
        self.tokens = _tokenize(text)
        self.pos = 0

    def _cur(self) -> str | None:
        return self.tokens[self.pos][0] if self.pos < len(self.tokens) else None

    def _advance(self) -> None:
        self.pos += 1

    def parse(self) -> _Node:
        node = self._expr()
        if self._cur() is not None:
            raise ValueError(f"trailing token {self._cur()!r}")
        return node

    def _expr(self) -> _Node:
        left = self._term()
        while self._cur() in ("+", "-"):
            self._advance()
            right = self._term()
            left = _Node("add", args=[left, right])
        return left

    def _term(self) -> _Node:
        left = self._unary()
        while self._cur() in ("*", "/"):
            op = self._cur()
            self._advance()
            right = self._unary()
            if op == "*":
                left = _Node("mul", args=[left, right])
            else:
                left = _Node("mul", args=[left, _Node("pow", a=right, b=_Node("num", a=-1))])
        return left

    def _unary(self) -> _Node:
        if self._cur() == "-" and self.pos + 1 < len(self.tokens) and self.tokens[self.pos + 1][0] == "num":
            return self._power()
        if self._cur() == "-":
            self._advance()
            return _Node("mul", args=[_Node("num", a=-1), self._unary()])
        return self._power()

    def _power(self) -> _Node:
        base = self._signed_primary()
        if self._cur() == "^":
            self._advance()
            exponent = self._unary()
            return _Node("pow", a=base, b=exponent)
        return base

    def _signed_primary(self) -> _Node:
        if self._cur() == "-" and self.pos + 1 < len(self.tokens) and self.tokens[self.pos + 1][0] == "num":
            self._advance()
            value = -int(self.tokens[self.pos][1])
            self._advance()
            return _Node("num", a=value)
        return self._primary()

    def _primary(self) -> _Node:
        token = self._cur()
        if token == "num":
            value = int(self.tokens[self.pos][1])
            self._advance()
            return _Node("num", a=value)
        if token == "id":
            name = self.tokens[self.pos][1]
            self._advance()
            if self._cur() == "(":
                self._advance()
                args: list[_Node] = []
                if self._cur() != ")":
                    while True:
                        args.append(self._expr())
                        if self._cur() != ",":
                            break
                        self._advance()
                if self._cur() != ")":
                    raise ValueError("unbalanced function call")
                self._advance()
                return _Node("fun", name=name, args=args)
            return _Node("var", name=name)
        if token == "(":
            self._advance()
            inner = self._expr()
            if self._cur() != ")":
                raise ValueError("unbalanced parentheses")
            self._advance()
            return inner
        raise ValueError(f"unexpected token {token!r}")


def parse_expr(text: str) -> _Node | None:
    """Parse an integrand, or return None when it does not parse."""
    try:
        return _Parser(text).parse()
    except ValueError:
        return None


def collect_heads(node: _Node, out: list[str]) -> None:
    """Collect the distinct function heads appearing anywhere in the tree."""
    if node.kind == "fun":
        out.append(node.name)
        for arg in node.args:
            collect_heads(arg, out)
    elif node.kind in ("add", "mul"):
        for arg in node.args:
            collect_heads(arg, out)
    elif node.kind == "pow":
        collect_heads(node.a, out)
        collect_heads(node.b, out)


def _contains_fractional_pow(node: _Node) -> bool:
    """Exact port of the bench helper: a `Pow` is fractional iff its exponent
    node is not a `Num`.  Note the bench quirk: the base and exponent trees of a
    `Pow` are NOT descended into (so a fractional power buried inside another
    `Pow` does not count), which is why a text-level heuristic disagrees with
    the bench.  Reproduced verbatim on purpose; it matches all stored buckets.
    """
    if node.kind == "pow":
        return node.b.kind != "num"
    if node.kind in ("add", "mul"):
        return any(_contains_fractional_pow(arg) for arg in node.args)
    if node.kind == "fun":
        return any(_contains_fractional_pow(arg) for arg in node.args)
    return False


def bucket_of(node: _Node) -> str:
    """Port of the bench's `bucket()` (nine plan buckets)."""
    heads: list[str] = []
    collect_heads(node, heads)
    heads = sorted(set(heads))
    if not heads:
        if _contains_fractional_pow(node):
            return "radical"
        if node.kind in ("pow", "mul"):
            return "power-binomial"
        if node.kind in ("add", "num", "var"):
            return "rational"
        return "mixed-other"
    if len(heads) == 1:
        name = heads[0]
        if name in TRIG_HEADS:
            return "trig"
        if name in HYPERBOLIC_HEADS:
            return "hyperbolic"
        if name in INVERSE_HEADS:
            return "inverse-trig-hyper"
        if name in ("exp", "log"):
            return "exp-log"
        if name == "sqrt":
            return "radical"
        if name in SPECIAL_HEADS:
            return "special"
        return "mixed-other"
    return "mixed-other"


# ---------------------------------------------------------------------------
# polylog order extraction
# ---------------------------------------------------------------------------

_POLYLOG_CALL_RE = re.compile(r"(?<![A-Za-z0-9_])polylog\s*\(", re.IGNORECASE)


def _call_arguments(text: str, lparen: int) -> list[str]:
    """Split the arguments of the call whose '(' is at `lparen`."""
    depth = 1
    args: list[str] = []
    start = lparen + 1
    i = start
    while i < len(text) and depth:
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                args.append(text[start:i])
                return args
        elif ch == "," and depth == 1:
            args.append(text[start:i])
            start = i + 1
        i += 1
    return []  # unbalanced


_INT_RE = re.compile(r"[+-]?\d+")


def polylog_orders(text: str) -> list[int]:
    """Sorted unique integer `polylog(<s>, ...)` orders in `text`.

    `polylog(k, ...)` (symbolic order) contributes nothing; the caller's
    `reason` records that case separately.
    """
    orders: set[int] = set()
    for match in _POLYLOG_CALL_RE.finditer(text):
        args = _call_arguments(text, match.end() - 1)
        if not args:
            continue
        first = args[0].strip()
        if _INT_RE.fullmatch(first):
            orders.add(int(first))
    return sorted(orders)


def has_symbolic_polylog_order(text: str) -> bool:
    """True when some `polylog(...)` call has a non-integer order argument."""
    for match in _POLYLOG_CALL_RE.finditer(text):
        args = _call_arguments(text, match.end() - 1)
        if not args:
            continue
        if not _INT_RE.fullmatch(args[0].strip()):
            return True
    return False


# ---------------------------------------------------------------------------
# Fragment-shape heuristic (class B rows)
# ---------------------------------------------------------------------------

_LOG_CALL_RE = re.compile(r"(?<![A-Za-z0-9_])log\s*\(", re.IGNORECASE)


def _log_arguments(node: _Node, out: list[_Node]) -> None:
    if node.kind == "fun":
        if node.name.lower() == "log":
            out.extend(node.args)
        for arg in node.args:
            _log_arguments(arg, out)
    elif node.kind in ("add", "mul"):
        for arg in node.args:
            _log_arguments(arg, out)
    elif node.kind == "pow":
        _log_arguments(node.a, out)
        _log_arguments(node.b, out)


def contains_var(node: _Node, var: str) -> bool:
    if node.kind == "num":
        return False
    if node.kind == "var":
        return node.name == var
    if node.kind == "fun":
        return any(contains_var(arg, var) for arg in node.args)
    if node.kind in ("add", "mul"):
        return any(contains_var(arg, var) for arg in node.args)
    if node.kind == "pow":
        return contains_var(node.a, var) or contains_var(node.b, var)
    return False


def affine_in(node: _Node, var: str) -> bool:
    """True when `node` is affine in `var` (constants are affine/degree 0)."""
    if not contains_var(node, var):
        return True
    if node.kind == "var":
        return node.name == var
    if node.kind == "add":
        return all(affine_in(arg, var) for arg in node.args)
    if node.kind == "mul":
        varying = [arg for arg in node.args if contains_var(arg, var)]
        return len(varying) == 1 and affine_in(varying[0], var)
    if node.kind == "pow":
        base_var = contains_var(node.a, var)
        exp_var = contains_var(node.b, var)
        if not base_var and not exp_var:
            return True
        if base_var and not exp_var:
            return node.b.kind == "num" and node.b.a == 1 and affine_in(node.a, var)
        return False
    if node.kind == "fun":  # a function of the variable is not affine
        return False
    return False


def _flatten(node: _Node, kind: str) -> list[_Node]:
    if node.kind == kind:
        out: list[_Node] = []
        for arg in node.args:
            out.extend(_flatten(arg, kind))
        return out
    return [node]


def _node_key(node: _Node):
    """Structural key so two occurrences of the same denominator compare equal."""
    if node.kind == "num":
        return ("num", node.a)
    if node.kind == "var":
        return ("var", node.name)
    if node.kind in ("add", "mul"):
        return (node.kind, tuple(_node_key(arg) for arg in node.args))
    if node.kind == "pow":
        return ("pow", _node_key(node.a), _node_key(node.b))
    if node.kind == "fun":
        return ("fun", node.name, tuple(_node_key(arg) for arg in node.args))
    return (node.kind,)


def _mobius_term(term: _Node, var: str) -> tuple[int, tuple] | None:
    """Analyse one additive term as `c * A(x) / D(x)` with affine A, D.

    Returns `(numerator_var_factors, denominator_keys)` or None when the term
    is not of that form.
    """
    numerator: list[_Node] = []
    denominators: list[_Node] = []
    for factor in _flatten(term, "mul"):
        if not contains_var(factor, var):
            continue
        if affine_in(factor, var):
            numerator.append(factor)
        elif (
            factor.kind == "pow"
            and factor.b.kind == "num"
            and factor.b.a == -1
            and affine_in(factor.a, var)
        ):
            denominators.append(factor.a)
        else:
            return None
    return len(numerator), tuple(_node_key(den) for den in denominators)


def mobius_in(node: _Node, var: str) -> bool:
    """True when `node` is a Mobius transform of `var` (affine/affine).

    Also accepts a constant power of an affine expression (its logarithm is
    again a logarithm of an affine expression).  Used only for the heuristic
    fragment hint, so it is deliberately conservative.
    """
    if not contains_var(node, var) or affine_in(node, var):
        return False
    if node.kind == "pow" and not contains_var(node.b, var):
        # log(u^k) = k*log(u): reducible when u itself is affine/Mobius.
        return affine_in(node.a, var) or mobius_in(node.a, var)
    terms = _flatten(node, "add")
    analysed: list[tuple[int, tuple]] = []
    for term in terms:
        analysis = _mobius_term(term, var)
        if analysis is None:
            return False
        numerator_count, denominator_keys = analysis
        if numerator_count == 0 and not denominator_keys:
            continue  # constant term
        if numerator_count > 1:
            return False  # product of two variable factors: degree >= 2
        analysed.append(analysis)
    if not analysed:
        return False
    denominator_sets = {denominators for _num, denominators in analysed}
    if len(denominator_sets) != 1:
        return False  # mixed denominators: common-denominator numerator is not affine
    return len(next(iter(denominator_sets))) >= 1


def _is_constant_times_affine_power(node: _Node, var: str) -> bool:
    """True for `c * u(x)^k` with `u` affine and `k` not containing the
    variable, e.g. `log(c*x^n)` (whose logarithm derivative is rational)."""
    varying = [factor for factor in _flatten(node, "mul") if contains_var(factor, var)]
    if len(varying) != 1:
        return False
    factor = varying[0]
    return (
        factor.kind == "pow"
        and not contains_var(factor.b, var)
        and affine_in(factor.a, var)
    )


def _log_arg_shape(node: _Node, var: str) -> str:
    """Classify one log argument: affine / mobius / power / non-affine."""
    if not contains_var(node, var):
        return "constant"
    if affine_in(node, var):
        return "affine"
    if _is_constant_times_affine_power(node, var):
        return "power"
    if mobius_in(node, var):
        return "mobius"
    return "non-affine"


def _collect_weight_denominators(node: _Node, var: str, out: list[_Node]) -> None:
    """Collect `Pow(..., -k)` factors of the rational weight.

    Does not descend into function arguments: a division inside
    `log(u(x))` is part of the log argument, not of the rational weight.
    """
    if node.kind == "fun":
        return
    if node.kind == "pow":
        if node.b.kind == "num" and node.b.a < 0 and contains_var(node.a, var):
            out.append(node)
        _collect_weight_denominators(node.a, var, out)
        _collect_weight_denominators(node.b, var, out)
    elif node.kind in ("add", "mul"):
        for arg in node.args:
            _collect_weight_denominators(arg, var, out)


def _denominator_shape(node: _Node, var: str) -> str:
    """Classify the variable-carrying denominators of the rational weight.

    `linear` means: at least one variable-carrying denominator factor exists
    and every such factor is `affine^(-1)` (first-order poles only, the
    partial-fraction reducible shape).  `nonlinear` means a variable-carrying
    denominator is not first-order-linear; `none` means no variable-carrying
    denominator at all.
    """
    denominators: list[_Node] = []
    _collect_weight_denominators(node, var, denominators)
    if not denominators:
        return "none"
    for factor in denominators:
        if factor.b.a != -1 or not affine_in(factor.a, var):
            return "nonlinear"
    return "linear"


def fragment_hint(orders: list[int], integrand: _Node | None, var: str) -> tuple[str, str]:
    """Return `(hint, shape_detail)` for a class-B row."""
    order_set = set(orders)
    if not order_set:
        return "higher-order", "answer has no numeric polylog order"
    if max(order_set) == 3:
        return "order3+", f"highest answer order {max(order_set)} (needs 0.31b)"
    if max(order_set) >= 4:
        return "higher-order", f"highest answer order {max(order_set)} (needs 0.31b)"
    # orders == {2}
    if integrand is None:
        return "non-affine", "integrand did not parse; shape not inspected"
    log_args: list[_Node] = []
    _log_arguments(integrand, log_args)
    if not log_args:
        return "non-affine", "no log() kernel in the integrand"
    if len(log_args) > 1:
        return "non-affine", f"{len(log_args)} log() kernels (single-kernel fragment)"
    shape = _log_arg_shape(log_args[0], var)
    denominator = _denominator_shape(integrand, var)
    detail = f"log argument {shape}; denominator {denominator}"
    if shape in ("affine", "mobius", "power") and denominator == "linear":
        return "candidate", detail
    return "non-affine", detail


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------


def read_answers(path: Path) -> list[tuple[str, str, str, str]]:
    """Read `id\\tintegrand\\tvar\\tanswer` rows (answer optional)."""
    rows: list[tuple[str, str, str, str]] = []
    malformed = 0
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split("\t", 3)
        if len(parts) < 3:
            malformed += 1
            print(f"warning: {path}:{lineno}: expected >=3 TSV fields, skipped", file=sys.stderr)
            continue
        answer = parts[3] if len(parts) >= 4 else ""
        rows.append((parts[0], parts[1], parts[2], answer))
    if malformed:
        print(f"warning: {malformed} malformed row(s) skipped", file=sys.stderr)
    return rows


def read_failures(path: Path) -> dict[str, dict]:
    """Read the JSONL failure dump into `{id: record}` (first occurrence wins)."""
    records: dict[str, dict] = {}
    duplicates = 0
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            print(f"warning: {path}:{lineno}: bad JSON ({exc}); skipped", file=sys.stderr)
            continue
        if not isinstance(record, dict):
            continue
        cid = record.get("id")
        if not isinstance(cid, str):
            continue
        if cid in records:
            duplicates += 1
            continue
        records[cid] = record
    if duplicates:
        print(f"warning: {duplicates} duplicate failure id(s) ignored", file=sys.stderr)
    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _order_label(orders: list[int]) -> str:
    return ";".join(str(order) for order in orders)


def _histogram_lines(title: str, counter: Counter) -> list[str]:
    lines = [f"  {title}:"]
    if not counter:
        lines.append("    (empty)")
        return lines
    for key in sorted(counter):
        lines.append(f"    {key}: {counter[key]}")
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Attribute the 0.30.0 unsolved Rubi cases that need polylog/Li2.",
    )
    parser.add_argument(
        "answers",
        nargs="?",
        type=Path,
        default=DEFAULT_ANSWERS,
        help=f"answer-augmented corpus TSV (default: {DEFAULT_ANSWERS})",
    )
    parser.add_argument(
        "--failures",
        type=Path,
        default=DEFAULT_FAILURES,
        help=f"unsolved-case JSONL dump (default: {DEFAULT_FAILURES})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help=f"output CSV (default: {DEFAULT_OUT})",
    )
    args = parser.parse_args(argv)

    answers_path = args.answers if args.answers.is_absolute() else REPO_ROOT / args.answers
    failures_path = args.failures if args.failures.is_absolute() else REPO_ROOT / args.failures
    out_path = args.out if args.out.is_absolute() else REPO_ROOT / args.out

    if not answers_path.exists():
        print(
            f"answers file not found: {answers_path}\n"
            "This file is generated benchmark data (gitignored). Fetch it first:\n"
            "\n    uv run python ocas-tests/scripts/fetch_rubi_corpus.py \\\n"
            "        --answers-out ocas-tests/data/rubi_1892_answers.tsv\n",
            file=sys.stderr,
        )
        return 0
    if not failures_path.exists():
        print(
            f"failures file not found: {failures_path}\n"
            "This file is generated benchmark data (gitignored); run\n"
            "\n    cargo bench -p ocas-tests --bench integrate_1892\n"
            "\nto produce ocas-tests/data/integrate_1892_failures.jsonl, or pass "
            "an existing dump via --failures.\n",
            file=sys.stderr,
        )
        return 0

    answers = read_answers(answers_path)
    failures = read_failures(failures_path)
    answer_ids = {row[0] for row in answers}

    failures_missing_from_answers = sorted(set(failures) - answer_ids)
    attr_rows: list[dict] = []
    bucket_mismatches: list[tuple[str, str, str]] = []
    parse_failures = 0

    for cid, integrand, var, answer in answers:
        if cid not in failures:
            continue
        integrand_lower = integrand.lower()
        answer_lower = answer.lower()
        if "polylog" in integrand_lower:
            cls = "A"
        elif "polylog" in answer_lower:
            cls = "B"
        else:
            continue  # not a polylog-attribution case; only counted in summary
        orders = polylog_orders(answer)
        tree = parse_expr(integrand)
        if tree is None:
            parse_failures += 1
            bucket = "mixed-other"
        else:
            bucket = bucket_of(tree)
        stored_bucket = failures[cid].get("bucket")
        if isinstance(stored_bucket, str) and stored_bucket != bucket:
            bucket_mismatches.append((cid, stored_bucket, bucket))

        if cls == "A":
            hint = ""
            if orders:
                reason = (
                    f"integrand contains polylog (A-class); answer needs orders "
                    f"{_order_label(orders)}"
                )
            elif has_symbolic_polylog_order(answer):
                reason = "integrand contains polylog (A-class); answer has symbolic polylog order"
            else:
                reason = "integrand contains polylog (A-class); answer has no numeric order"
        else:
            hint, shape_detail = fragment_hint(orders, tree, var)
            if has_symbolic_polylog_order(answer):
                symbolic = "; answer also has a symbolic polylog order"
            else:
                symbolic = ""
            reason = f"elementary integrand needs Li2; {shape_detail}{symbolic}"

        attr_rows.append(
            {
                "id": cid,
                "bucket": bucket,
                "integrand": integrand,
                "answer_polylog_orders": _order_label(orders),
                "class": cls,
                "fragment_hint": hint,
                "reason": reason,
            }
        )

    attr_rows.sort(key=lambda row: row["id"])
    class_a = [row for row in attr_rows if row["class"] == "A"]
    class_b = [row for row in attr_rows if row["class"] == "B"]
    candidates = [row for row in class_b if row["fragment_hint"] == "candidate"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(attr_rows)

    # ---- deterministic summary -------------------------------------------
    print("attribute_031: Rubi 1892 answer-aware polylog attribution")
    print(f"  answers:   {answers_path} ({len(answers)} cases)")
    print(f"  failures:  {failures_path} ({len(failures)} unsolved ids)")
    if failures_missing_from_answers:
        preview = ", ".join(failures_missing_from_answers[:10])
        more = "" if len(failures_missing_from_answers) <= 10 else ", ..."
        print(
            f"  warning:   {len(failures_missing_from_answers)} unsolved id(s) not in "
            f"the answers TSV ({preview}{more})"
        )
    print(f"  unsolved total: {len(failures)}")
    print(
        f"  polylog attribution rows: {len(attr_rows)} "
        f"(class A {len(class_a)}, class B {len(class_b)}); "
        f"{len(failures) - len(attr_rows)} unsolved case(s) have no polylog signature "
        "and stay out of the CSV"
    )
    print(f"  class A ids ({len(class_a)}): {', '.join(row['id'] for row in class_a) or '(none)'}")
    print(f"  class B ids ({len(class_b)}): {', '.join(row['id'] for row in class_b) or '(none)'}")

    bucket_counts = Counter(row["bucket"] for row in class_b)
    print(f"  class B buckets ({len(class_b)}):")
    for bucket in sorted(bucket_counts):
        print(f"    {bucket:>18}: {bucket_counts[bucket]}")

    def _orders_of(row: dict) -> tuple[int, ...]:
        return tuple(int(part) for part in row["answer_polylog_orders"].split(";") if part)

    for line in _histogram_lines(
        "class A answer polylog-order histogram",
        Counter(_orders_of(row) for row in class_a),
    ):
        print(line)
    for line in _histogram_lines(
        "class B answer polylog-order histogram",
        Counter(_orders_of(row) for row in class_b),
    ):
        print(line)

    hint_counts = Counter(row["fragment_hint"] for row in class_b)
    print(f"  class B fragment hints ({len(class_b)}):")
    for hint in sorted(hint_counts):
        print(f"    {hint:>12}: {hint_counts[hint]}")

    print(
        f"  candidate ids ({len(candidates)}): "
        f"{', '.join(row['id'] for row in candidates) or '(none)'}"
    )
    candidate_buckets = Counter(row["bucket"] for row in candidates)
    for bucket in sorted(candidate_buckets):
        ids = [row["id"] for row in candidates if row["bucket"] == bucket]
        print(f"    {bucket:>18}: {', '.join(ids)}")

    if bucket_mismatches:
        print(
            f"  warning: {len(bucket_mismatches)} bucket replica mismatch(es) against the "
            f"stored bench buckets, e.g. {bucket_mismatches[:5]}"
        )
    else:
        print("  bucket replica cross-check: all attributed rows match the stored bench buckets")
    if parse_failures:
        print(f"  warning: {parse_failures} integrand(s) did not parse for bucketing")

    print(f"  wrote {out_path} ({len(attr_rows)} rows, {len(CSV_HEADER)} columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
