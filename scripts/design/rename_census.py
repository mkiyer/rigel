#!/usr/bin/env python
"""Which names does a vocabulary ruling touch, and which of them carry two senses?

The vocabulary is REGION (a genomic interval; was node), BOUNDARY (a single genomic position between two
regions; was edge, cut, line, seam) and SPLICE JUNCTION (a connection between two boundaries; `sj` is
allowed, bare `junction` is not), with splice donor / acceptor banned as names for a genomically-ordered
pair. This reports every site the ruling touches, by kind (Python identifiers by class / function /
argument / attribute / variable from the AST, C++ identifiers, prose in docstrings, docs and the memory
directory, file names, and the payload wire schema, whose renaming refuses every cache), and it never
renames. The dangerous part of a bulk rename is not the volume but the tokens that carry two senses,
where a correct-looking global replace silently corrupts the other one (`line` as a 0-bp boundary and as a
line of text, `node` as an interval and as an `ast` node, `donor` as a splice donor and as the toy
harness's source condition, `cut` at two extents); those are listed with no automatic target and
`--sense` dumps every site of one with context for a per-site ruling. Names in `EXEMPT` are correct as
they stand and must survive the rename.

Usage::

    python scripts/design/rename_census.py              # the whole landscape
    python scripts/design/rename_census.py --sense cut  # every site of one ambiguous token, with context
"""

from __future__ import annotations

import argparse
import ast
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

_REPO = Path(__file__).resolve().parents[2]

#: The tokens the ruling touches, and what each becomes. A ``None`` target is a token that needs a
#: per-site decision, which is the point of this census.
TOKENS = {
    "node": "region",
    "edge": "boundary",
    "line": None,  # ambiguous: a 0-bp boundary vs an ordinary line of text
    "seam": "boundary",
    "cut": None,  # ambiguous: position-including-termini vs interior boundary
    "junction": None,  # ambiguous: bare `junction` is banned, `splice_junction`/`sj` are not
    "donor": None,  # ambiguous: splice donor vs the toy harness's source condition
    "acceptor": None,
    #: a retired message policy. A surviving site is a record, a doc, or a label that still names it
    #: for a live mechanism; only the last is a defect, and the per-site dump tells the three apart.
    "relay": None,
}

#: Names that are correct and must survive the rename, each with its reason.
EXEMPT = {
    "splice_donor_acceptor": "sim/splice_motif.py — takes the STRAND and returns the GT..AG "
    "dinucleotides. The one place the biology term is used correctly.",
    "donor_dir": "the toy harness's SOURCE CONDITION, not a splice donor",
    "donor_on": "ditto", "donor_off": "ditto", "donor_name": "ditto", "donor_qname": "ditto",
    "toy_donor": "ditto", "_donor_sim_params": "ditto",
}

CODE_DIRS = ("src", "tests", "scripts")
DOC_GLOBS = ("docs/**/*.md", "CLAUDE.md", "README.md")
MEMORY = Path.home() / ".claude/projects/-Users-mkiyer-proj-rigel/memory"


def _tokens_in(name: str) -> list[str]:
    low = name.lower()
    return [t for t in TOKENS if re.search(rf"(?:^|[^a-z]){t}", low)]


class _Collect(ast.NodeVisitor):
    """Real identifiers, by kind; a grep cannot tell a class from a comment."""

    def __init__(self):
        self.by_kind: dict[str, Counter] = defaultdict(Counter)

    def visit_ClassDef(self, n):
        self.by_kind["class"][n.name] += 1
        self.generic_visit(n)

    def _fn(self, n):
        self.by_kind["function"][n.name] += 1
        for a in [*n.args.args, *n.args.kwonlyargs, *n.args.posonlyargs]:
            self.by_kind["argument"][a.arg] += 1
        self.generic_visit(n)

    visit_FunctionDef = _fn
    visit_AsyncFunctionDef = _fn

    def visit_Attribute(self, n):
        self.by_kind["attribute"][n.attr] += 1
        self.generic_visit(n)

    def visit_Name(self, n):
        self.by_kind["variable"][n.id] += 1
        self.generic_visit(n)


def python_identifiers():
    out: dict[str, Counter] = defaultdict(Counter)
    files: dict[str, set] = defaultdict(set)
    for d in CODE_DIRS:
        for p in sorted((_REPO / d).rglob("*.py")):
            try:
                tree = ast.parse(p.read_text())
            except SyntaxError:
                continue
            c = _Collect()
            c.visit(tree)
            for kind, names in c.by_kind.items():
                for name, k in names.items():
                    if _tokens_in(name):
                        out[kind][name] += k
                        files[name].add(str(p.relative_to(_REPO)))
    return out, files


def cpp_identifiers():
    out, files = Counter(), defaultdict(set)
    pat = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")
    for p in sorted((_REPO / "src").rglob("*.h")) + sorted((_REPO / "src").rglob("*.cpp")):
        for name in pat.findall(p.read_text()):
            if _tokens_in(name):
                out[name] += 1
                files[name].add(str(p.relative_to(_REPO)))
    return out, files


def prose_counts():
    """Docstrings/comments in code, plus docs, CLAUDE.md and the MEMORY directory."""
    buckets = Counter()
    per_token: dict[str, Counter] = defaultdict(Counter)

    def scan(path, bucket):
        try:
            text = path.read_text()
        except (OSError, UnicodeDecodeError):
            return
        for t in TOKENS:
            n = len(re.findall(rf"(?:^|[^a-zA-Z]){t}", text, re.I))
            if n:
                buckets[bucket] += n
                per_token[t][bucket] += n

    for g in DOC_GLOBS:
        for p in sorted(_REPO.glob(g)):
            scan(p, "docs")
    if MEMORY.is_dir():
        for p in sorted(MEMORY.rglob("*.md")):
            scan(p, "memory")
    for d in CODE_DIRS:
        for p in sorted((_REPO / d).rglob("*.py")):
            try:
                tree = ast.parse(p.read_text())
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                doc = ast.get_docstring(node) if isinstance(
                    node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
                ) else None
                if doc:
                    for t in TOKENS:
                        n = len(re.findall(rf"(?:^|[^a-zA-Z]){t}", doc, re.I))
                        if n:
                            buckets["docstrings"] += n
                            per_token[t]["docstrings"] += n
    return buckets, per_token


def filenames():
    return [
        str(p.relative_to(_REPO))
        for d in CODE_DIRS
        for p in sorted((_REPO / d).rglob("*.py"))
        if _tokens_in(p.stem)
    ]


def payload_names():
    sys.path.insert(0, str(_REPO / "src"))
    from rigel.scan_cache import _schema_names

    return [n for n in _schema_names() if _tokens_in(n)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sense", choices=sorted(TOKENS), default=None,
                    help="dump every site of one token, with context, for a per-site ruling")
    args = ap.parse_args()

    if args.sense:
        pat = re.compile(rf"(?:^|[^a-zA-Z]){args.sense}", re.I)
        for d in CODE_DIRS:
            for p in sorted((_REPO / d).rglob("*.py")):
                for i, line in enumerate(p.read_text().splitlines(), 1):
                    if pat.search(line):
                        print(f"{p.relative_to(_REPO)}:{i}: {line.strip()[:120]}")
        return 0

    py, py_files = python_identifiers()
    cpp, cpp_files = cpp_identifiers()
    prose, per_token = prose_counts()

    print("=" * 100)
    print("  ⭐ PYTHON IDENTIFIERS, by kind — what a grep cannot tell apart")
    print("=" * 100)
    for kind in ("class", "function", "argument", "attribute", "variable"):
        names = py.get(kind, Counter())
        if not names:
            continue
        print(f"\n  {kind.upper()}  ({len(names)} distinct, {sum(names.values())} occurrences)")
        for name, n in names.most_common(12):
            tag = "  ⭐EXEMPT" if name in EXEMPT else ""
            print(f"     {n:5d}  {name:<44s} {len(py_files[name])} files{tag}")
        if len(names) > 12:
            print(f"     … and {len(names) - 12} more distinct names")

    print("\n" + "=" * 100)
    print(f"  ⭐ C++ IDENTIFIERS  ({len(cpp)} distinct, {sum(cpp.values())} occurrences)")
    print("=" * 100)
    for name, n in cpp.most_common(20):
        print(f"     {n:5d}  {name:<44s} {len(cpp_files[name])} files")

    print("\n" + "=" * 100)
    print("  ⭐ PROSE — the half a code-only rename leaves behind and rots")
    print("=" * 100)
    for bucket, n in prose.most_common():
        print(f"     {n:6d}  {bucket}")
    print(f"\n     {'token':<12s} " + "".join(f"{b:>12s}" for b in ("docstrings", "docs", "memory")))
    for t in sorted(TOKENS):
        row = per_token[t]
        if sum(row.values()):
            print(f"     {t:<12s} " + "".join(f"{row.get(b, 0):12d}"
                                              for b in ("docstrings", "docs", "memory")))

    print("\n" + "=" * 100)
    print("  ⛔ FILENAMES — an import-graph change, not a text change")
    print("=" * 100)
    for f in filenames():
        print(f"     {f}")

    print("\n" + "=" * 100)
    print("  ⛔⛔ THE PAYLOAD WIRE SCHEMA — renaming any of these REFUSES EVERY CACHE")
    print("=" * 100)
    for n in payload_names():
        print(f"     {n}")

    print("\n" + "=" * 100)
    print("  ⛔⛔⛔ AMBIGUOUS TOKENS — a global replace is WRONG for each of these")
    print("=" * 100)
    for t, target in TOKENS.items():
        if target is None:
            print(f"     {t:<10s} needs a per-site ruling   (--sense {t})")
    print("\n  ⭐ EXEMPT, each for a measured reason:")
    for name, why in EXEMPT.items():
        print(f"     {name:<24s} {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
