"""0.30.0 Wave 0: shape attribution for the radical and inverse-trig-hyper
failure buckets of the 1892-case corpus.

Reads the per-case failure dump (JSONL) and classifies each failure by
regex/scanner-level shape features. The classes are estimates (same caveat
as the 0.27.3 investigation: regex bucketing, magnitude reference only).

Usage:
    uv run python ocas-tests/scripts/attribute_030.py \
        ocas-tests/data/failures_final_029.jsonl
"""

import json
import re
import sys
from collections import Counter
from pathlib import Path

INV_HEADS = [
    "asinh", "acosh", "atanh", "acoth", "asech", "acsch",
    "asin", "acos", "atan", "acot", "asec", "acsc",
]
TRIG_HEADS = {"sin", "cos", "tan", "sec", "csc", "cot"}

# Numeric fractional power `base^(p/q)` (also `base^(-p/q)`).
POW_WITH_BASE = re.compile(r"(\([^()]*\)|[A-Za-z_][\w]*)\s*\^\(?\s*(-?\d+)\s*/\s*(\d+)\s*\)?")


def has_symbolic_exponent(integrand: str) -> bool:
    """True when some `^` exponent contains a letter (a symbolic exponent).

    Replaces the 0.30.0-Wave-0 regex `\\^\\(?[a-zA-Z(]`, whose character
    class contained `(` — any exponent written with nested parens (e.g.
    `^((1+x)/(1-x))` … but also numeric `(p)/(q)` shapes starting with
    `((`) matched it, inflating the symbolic-exponent bucket by ~143 cases
    that were really numeric fractional powers.
    """
    for m in re.finditer(r"\^", integrand):
        i = m.end()
        while i < len(integrand) and integrand[i] in " \t":
            i += 1
        if i < len(integrand) and integrand[i] == "(":
            depth = 1
            j = i + 1
            while j < len(integrand) and depth:
                if integrand[j] == "(":
                    depth += 1
                elif integrand[j] == ")":
                    depth -= 1
                j += 1
            if re.search(r"[a-zA-Z]", integrand[i:j]):
                return True
        else:
            tok = re.match(r"-?[\w./]+", integrand[i:])
            if tok and re.search(r"[a-zA-Z]", tok.group(0)):
                return True
    return False


def sqrt_bases(integrand: str):
    """Extract the argument of every `sqrt(...)` call (balanced parens)."""
    bases = []
    for m in re.finditer(r"sqrt\(", integrand):
        depth = 1
        i = m.end()
        while i < len(integrand) and depth:
            if integrand[i] == "(":
                depth += 1
            elif integrand[i] == ")":
                depth -= 1
            i += 1
        bases.append(integrand[m.end() : i - 1])
    return bases


def frac_pow_bases(integrand: str):
    bases = []
    for m in POW_WITH_BASE.finditer(integrand):
        base, p, q = m.group(1), int(m.group(2)), int(m.group(3))
        if q != 0 and p % q != 0:
            bases.append(base)
    return bases


def radicand_degree(base: str) -> str:
    """Rough degree of the radicand in x."""
    powers = [int(p) for p in re.findall(r"x\^(\d+)", base)]
    deg = max(powers, default=(1 if re.search(r"(^|[^a-z])x([^a-z0-9]|$)", base) else 0))
    if deg == 0:
        return "constant"
    if deg == 1:
        return "linear"
    if deg == 2:
        return "quadratic"
    return f"deg{deg}+"


def classify_radical(integrand: str) -> str:
    heads = set(re.findall(r"([a-z]+)\(", integrand))
    has_trig = bool(heads & TRIG_HEADS)
    has_inv = bool(heads & set(INV_HEADS))
    has_log_exp = bool(heads & {"log", "exp"})

    if has_symbolic_exponent(integrand):
        cls = "symbolic-exponent"
    else:
        bases = sqrt_bases(integrand) + frac_pow_bases(integrand)
        distinct = set(bases)
        if not bases:
            cls = "unclassified"
        elif len(distinct) >= 2:
            cls = "multi-radical-base"
        else:
            deg = radicand_degree(bases[0])
            cls = f"single-radical-{deg}"
            if bases[0] == "x":
                cls = "frac-power-of-x"

    suffix = []
    if has_trig:
        suffix.append("trig")
    if has_inv:
        suffix.append("inv")
    if has_log_exp:
        suffix.append("logexp")
    return cls + ("|" + "+".join(suffix) if suffix else "")


def classify_invtrig(integrand: str) -> str:
    heads = [h for h in INV_HEADS if re.search(rf"\b{h}\(", integrand)]
    family = "invtrig" if heads and heads[0] in {"asin", "acos", "atan", "acot", "asec", "acsc"} else "invhyp"

    feats = []
    # (a + b*F(u))^k kernel power, k != 1
    if re.search(r"\((?:[a-zA-Z][\w]*\s*\+\s*)?[a-zA-Z]?\*?(?:" + "|".join(INV_HEADS) + r")\([^()]*\)\)\s*\^", integrand):
        feats.append("kernel-power")
    # Bare F(u)^k
    if re.search(r"(?:" + "|".join(INV_HEADS) + r")\([^()]*\)\s*\^", integrand):
        feats.append("bare-power")
    # F in a denominator
    if re.search(r"/\s*\([^)]*(?:" + "|".join(INV_HEADS) + r")", integrand):
        feats.append("F-in-denom")
    # Radical alongside
    if sqrt_bases(integrand) or frac_pow_bases(integrand):
        feats.append("radical")
    # Argument shape of the first inverse head
    m = re.search(r"(?:" + "|".join(INV_HEADS) + r")\(([^()]*)\)", integrand)
    if m:
        arg = m.group(1).replace(" ", "")
        if re.search(r"x\^[2-9]", arg):
            feats.append("nonlinear-arg")
        elif re.fullmatch(r"[a-z]?\*?x", arg):
            feats.append("zero-int-arg")
        else:
            feats.append("affine-arg")
    return family + ":" + ("+".join(feats) if feats else "bare")


def main() -> None:
    path = Path(sys.argv[1])
    out_dir = path.parent
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]

    for bucket, classify, out_name in [
        ("radical", classify_radical, "radical_attribution_030.csv"),
        ("inverse-trig-hyper", classify_invtrig, "invtrig_attribution_030.csv"),
    ]:
        entries = [r for r in rows if r["bucket"] == bucket]
        counter = Counter()
        classified = []
        for r in entries:
            cls = classify(r["integrand"])
            counter[cls] += 1
            classified.append((r["id"], cls, r["integrand"]))
        out = out_dir / out_name
        with out.open("w", encoding="utf-8", newline="") as fh:
            fh.write("id,class,integrand\n")
            for rid, cls, integ in classified:
                fh.write(f'{rid},"{cls}","{integ}"\n')
        print(f"=== {bucket}: {len(entries)} failures -> {out.name}")
        for cls, n in counter.most_common(24):
            print(f"  {n:4d}  {cls}")
        rest = sum(n for cls, n in counter.most_common()[24:])
        if rest:
            print(f"  {rest:4d}  (tail, {len(counter) - 24} more classes)")
        print()


if __name__ == "__main__":
    main()
