#!/usr/bin/env python3
"""Fetch the Rubi Axiom test corpus and build a deterministic oCAS sample.

Source: Rubi's official Axiom-syntax test-suite archive

    https://rulebasedintegration.org/TestFiles/AxiomSyntaxFiles/AxiomSyntaxTestFiles.zip

The archive's SHA-256 is pinned (recorded in the symbolica-integrate
repository's `tests/data/RUBITestFiles-README.md`, downloaded 2026-07-19):

    d9debdd24ada7f6c6aa7c87bde75f61c9bbbe5083a127dfb75f1fc9469c66653

A digest mismatch aborts unless `--ignore-digest` is passed, so a run is
reproducible from the fixed source + fixed seed (default 1892).

Fallback source (same corpus, same format): the extracted `.input` files
mirrored in the MIT-licensed symbolica-integrate repository
(https://github.com/symbolica-dev/symbolica-integrate), fetched as the
repository tarball when the primary URL is unreachable.

Each corpus case is a single line beginning with `[`; its first top-level
comma field is the integrand and the second is the integration variable
(Axiom syntax). This script converts integrands to oCAS parse syntax
(`Sin[x]`-style Axiom forms do not occur here; the corpus uses lowercase
parenthesised calls such as `sin(x)`), skips cases it cannot convert, and
deterministically samples 1892 convertible cases using the fixed seed.

Output: `ocas-tests/data/rubi_1892.tsv` (one row per case:
`id\\tintegrand\\tvar`) plus `ocas-tests/data/rubi_1892.meta.json` with
provenance (source URL, digest, seed, sample size). The corpus itself is
cached under `ocas-tests/data/cache/` and is gitignored; it is never
committed to the repository.

Run via uv (or any Python 3.11+; standard library only):

    uv run python ocas-tests/scripts/fetch_rubi_corpus.py
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

CORPUS_URL = (
    "https://rulebasedintegration.org/TestFiles/AxiomSyntaxFiles/"
    "AxiomSyntaxTestFiles.zip"
)
PINNED_SHA256 = "d9debdd24ada7f6c6aa7c87bde75f61c9bbbe5083a127dfb75f1fc9469c66653"
GITHUB_TARBALL_URL = "https://codeload.github.com/symbolica-dev/symbolica-integrate/tar.gz/refs/heads/main"
RUBI_DIR_PREFIX = "tests/data/RUBITestFiles/"

IDENT_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
INT_RE = re.compile(r"[0-9]+")

# Axiom function names -> oCAS function names. Names absent from this map
# (e.g. `sqrt`, `sin`, `exp`, `erf`, `Ei`) are already oCAS-compatible and
# pass through unchanged; unknown names are kept verbatim (the oCAS parser
# accepts any identifier as a function head, so the integrand still parses
# and is classified as an unsolved case rather than a conversion failure).
FN_MAP = {
    "ProductLog": "productlog",
    "GAMMA": "gamma",
    "lnGAMMA": "lngamma",
    "Psi": "psi",
    "Zeta": "zeta",
    "Li": "li",
    "Factorial": "factorial",
    "FresnelS": "fresnels",
    "FresnelC": "fresnelc",
}


class ConvertError(ValueError):
    """The Axiom expression uses syntax the converter does not support."""


class _Lexer:
    """Tokenise an Axiom integrand into oCAS-parseable tokens.

    Raises ConvertError on any character that cannot be part of an oCAS
    expression (floats, `!`, `_` subscripts, brackets, ...).
    """

    def __init__(self, src: str):
        self.src = src
        self.pos = 0

    def peek(self) -> str | None:
        while self.pos < len(self.src) and self.src[self.pos].isspace():
            self.pos += 1
        if self.pos >= len(self.src):
            return None
        ch = self.src[self.pos]
        if ch == "%":
            m = re.match(r"%[A-Za-z]+", self.src[self.pos:])
            if not m:
                raise ConvertError(f"bad % constant near {self.src[self.pos:]!r}")
            return m.group(0)
        if ch.isalpha() or ch == "_":
            m = IDENT_RE.match(self.src, self.pos)
            assert m is not None
            return m.group(0)
        if ch.isdigit():
            m = INT_RE.match(self.src, self.pos)
            assert m is not None
            return m.group(0)
        if ch in "+-*/^(),":
            return ch
        raise ConvertError(f"unsupported character {ch!r} near {self.src[self.pos:]!r}")

    def next(self) -> str | None:
        tok = self.peek()
        if tok is not None:
            self.pos += len(tok)
        return tok


class _Converter:
    """Recursive-descent Axiom -> oCAS converter.

    Grammar mirrors the oCAS parser (expr/term/factor/primary with the same
    precedence), so every emitted string parses in oCAS. Axiom's implicit
    multiplication (e.g. `Derivative(1)(F)(x)`) is not part of the grammar
    and raises ConvertError, which skips the case.
    """

    def __init__(self, src: str):
        self.lex = _Lexer(src)

    def _peek(self) -> str | None:
        return self.lex.peek()

    def _next(self) -> str | None:
        return self.lex.next()

    def convert(self) -> str:
        out = self._expr()
        if self._peek() is not None:
            raise ConvertError(f"trailing tokens near {self._peek()!r}")
        return out

    def _expr(self) -> str:
        parts = [self._term()]
        while True:
            tok = self._peek()
            if tok == "+":
                self._next()
                parts.append(" + ")
                parts.append(self._term())
            elif tok == "-":
                self._next()
                # A space before binary `-` is required: the oCAS lexer merges
                # `-<digit>` into a single negative-integer token, so `x-1`
                # would not parse (tokens: x, -1).
                parts.append(" - ")
                parts.append(self._term())
            else:
                break
        return "".join(parts)

    def _term(self) -> str:
        parts = [self._factor()]
        while True:
            tok = self._peek()
            if tok == "*":
                self._next()
                parts.append("*")
                parts.append(self._factor())
            elif tok == "/":
                self._next()
                parts.append("/")
                parts.append(self._factor())
            else:
                break
        return "".join(parts)

    def _factor(self) -> str:
        tok = self._peek()
        if tok == "-":
            self._next()
            return "-" + self._factor()
        base = self._primary()
        if self._peek() == "^":
            self._next()
            exp = self._factor()
            return f"{base}^({exp})" if exp.startswith("-") or "+" in exp else f"{base}^{exp}"
        return base

    def _primary(self) -> str:
        tok = self._next()
        if tok is None:
            raise ConvertError("unexpected end of integrand")
        if tok.isdigit():
            return tok
        if tok.startswith("%"):
            if tok == "%e":
                # `%e^X` -> `exp(X)`; standalone `%e` -> symbol `e`.
                if self._peek() == "^":
                    self._next()
                    exp = self._factor()
                    return f"exp({exp})"
                return "e"
            if tok == "%pi":
                return "pi"
            if tok == "%i":
                return "i"
            raise ConvertError(f"unsupported constant {tok!r}")
        if tok == "(":
            inner = self._expr()
            if self._next() != ")":
                raise ConvertError("unbalanced parentheses")
            return f"({inner})"
        # identifier: function call or variable
        name = FN_MAP.get(tok, tok)
        if self._peek() == "(":
            self._next()
            args = []
            if self._peek() != ")":
                while True:
                    args.append(self._expr())
                    if self._peek() != ",":
                        break
                    self._next()
            if self._next() != ")":
                raise ConvertError("unbalanced parentheses in call")
            return f"{name}({','.join(args)})"
        return name


def split_top_level(s: str) -> list[str]:
    """Split `s` on top-level commas (parens nest)."""
    parts: list[str] = []
    depth = 0
    cur: list[str] = []
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append("".join(cur).strip())
    return parts


def parse_cases(text: str, source: str) -> list[tuple[str, str, str]]:
    """Yield (integrand, var, source) for every case line in a corpus file."""
    cases: list[tuple[str, str, str]] = []
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("["):
            continue
        body = s[1:]
        if s.endswith("]"):
            body = s[1:-1]
        fields = split_top_level(body)
        if len(fields) < 2:
            continue
        integrand, var = fields[0], fields[1]
        if not IDENT_RE.fullmatch(var):
            continue
        cases.append((integrand, var, source))
    return cases


def splitmix64(seed: int) -> int:
    """Deterministic PRNG (SplitMix64) so sampling is version-independent."""
    seed &= (1 << 64) - 1
    while True:
        seed = (seed + 0x9E3779B97F4A7C15) & (1 << 64) - 1
        z = seed
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9 & (1 << 64) - 1
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB & (1 << 64) - 1
        yield z ^ (z >> 31)


def seeded_shuffle(items: list, seed: int) -> None:
    """In-place Fisher-Yates with the deterministic SplitMix64 PRNG."""
    rng = splitmix64(seed)
    for i in range(len(items) - 1, 0, -1):
        j = next(rng) % (i + 1)
        items[i], items[j] = items[j], items[i]


def download(url: str, dest: Path) -> None:
    print(f"downloading {url}")
    try:
        urllib.request.urlretrieve(url, dest)
    except urllib.error.URLError as e:
        raise RuntimeError(f"download failed: {url}: {e}") from e


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_corpus_zip(path: Path) -> list[tuple[str, str, str]]:
    cases: list[tuple[str, str, str]] = []
    with zipfile.ZipFile(path) as z:
        for name in sorted(z.namelist()):
            if not name.endswith(".input"):
                continue
            text = z.read(name).decode("utf-8", errors="replace")
            cases.extend(parse_cases(text, name))
    return cases


def load_corpus_tarball(path: Path) -> list[tuple[str, str, str]]:
    cases: list[tuple[str, str, str]] = []
    with tarfile.open(path, "r:gz") as t:
        for member in sorted(t.getmembers(), key=lambda m: m.name):
            if not member.isfile() or RUBI_DIR_PREFIX not in member.name:
                continue
            if not member.name.endswith(".input"):
                continue
            f = t.extractfile(member)
            if f is None:
                continue
            text = f.read().decode("utf-8", errors="replace")
            cases.extend(parse_cases(text, member.name))
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1892)
    parser.add_argument("--sample-size", type=int, default=1892)
    parser.add_argument("--out", type=Path, default=Path("ocas-tests/data/rubi_1892.tsv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("ocas-tests/data/cache"))
    parser.add_argument("--ignore-digest", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    cache_dir = args.cache_dir if args.cache_dir.is_absolute() else root / args.cache_dir
    out_path = args.out if args.out.is_absolute() else root / args.out
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    zip_path = cache_dir / "AxiomSyntaxTestFiles.zip"
    digest: str | None = None
    source_url = CORPUS_URL
    cases: list[tuple[str, str, str]] = []
    if zip_path.exists():
        digest = sha256_of(zip_path)
        print(f"using cached corpus {zip_path} (sha256 {digest[:16]}...)")
        if digest != PINNED_SHA256 and not args.ignore_digest:
            print(
                f"error: cached corpus digest {digest} != pinned {PINNED_SHA256}; "
                "pass --ignore-digest to proceed (breaks reproducibility)",
                file=sys.stderr,
            )
            return 2
        cases = load_corpus_zip(zip_path)
    else:
        try:
            download(CORPUS_URL, zip_path)
        except RuntimeError as e:
            print(f"warning: primary source failed ({e}); trying GitHub mirror")
            tar_path = cache_dir / "symbolica-integrate.tar.gz"
            try:
                download(GITHUB_TARBALL_URL, tar_path)
            except RuntimeError as e2:
                print(f"error: fallback source also failed: {e2}", file=sys.stderr)
                return 1
            digest = sha256_of(tar_path)
            cases = load_corpus_tarball(tar_path)
            source_url = GITHUB_TARBALL_URL
        else:
            digest = sha256_of(zip_path)
            print(f"sha256 {digest}")
            if digest != PINNED_SHA256 and not args.ignore_digest:
                print(
                    f"error: downloaded digest {digest} != pinned {PINNED_SHA256}; "
                    "pass --ignore-digest to proceed (breaks reproducibility)",
                    file=sys.stderr,
                )
                return 2
            cases = load_corpus_zip(zip_path)

    if not cases:
        print("error: no case lines found in corpus", file=sys.stderr)
        return 1
    print(f"parsed {len(cases)} cases")

    # Deterministic sampling: shuffle the whole corpus with the fixed seed,
    # convert in order, skip failures, take the first `sample_size` cases.
    seeded_shuffle(cases, args.seed)
    sampled: list[tuple[str, str, str, str]] = []  # (id, integrand, var, source)
    convertible = 0
    skipped = 0
    for idx, (integrand, var, source) in enumerate(cases):
        try:
            converted = _Converter(integrand).convert()
        except ConvertError:
            skipped += 1
            continue
        convertible += 1
        if len(sampled) >= args.sample_size:
            continue
        sampled.append((f"rubi-{idx + 1:05d}", converted, var, source))

    print(f"convertible {convertible}, skipped {skipped}, sampled {len(sampled)}")
    if len(sampled) < args.sample_size:
        print(
            f"warning: only {len(sampled)} cases convertible (wanted {args.sample_size}); "
            "coverage will be computed over the actual N",
            file=sys.stderr,
        )
    if len(sampled) < 1500:
        print(
            "error: fewer than 1500 convertible cases; switch corpus source and rerun",
            file=sys.stderr,
        )
        return 1

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        for cid, integrand, var, _src in sampled:
            f.write(f"{cid}\t{integrand}\t{var}\n")

    meta = {
        "seed": args.seed,
        "requested": args.sample_size,
        "n": len(sampled),
        "source_url": source_url,
        "source_sha256": digest,
        "pinned_sha256": PINNED_SHA256,
        "digest_matched": digest == PINNED_SHA256,
        "total_cases": len(cases),
        "convertible": convertible,
        "skipped": skipped,
        "date": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
    }
    meta_path = out_path.with_suffix(".meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"wrote {out_path} ({len(sampled)} rows) and {meta_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
