"""Survey 3: top-level defs/classes in src/rigel/*.py (excl. calibration) with no callers."""

import ast
import pathlib
import re
import subprocess

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
files = sorted(ROOT.glob("src/rigel/*.py"))

# build one big corpus of all python + notebook-ish text we search
corpus_files = (
    list(ROOT.glob("src/rigel/**/*.py"))
    + list(ROOT.glob("tests/**/*.py"))
    + list(ROOT.glob("scripts/**/*.py"))
    + list(ROOT.glob("src/rigel/native/**/*.cpp"))
    + list(ROOT.glob("src/rigel/native/**/*.h"))
)
texts = {}
for f in corpus_files:
    try:
        texts[f] = f.read_text()
    except Exception:
        pass

rows = []
for p in files:
    tree = ast.parse(p.read_text())
    src_lines = p.read_text().splitlines()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = node.name
            if name.startswith("__"):
                continue
            pat = re.compile(r"\b" + re.escape(name) + r"\b")
            hits = []
            for f, t in texts.items():
                for m in pat.finditer(t):
                    ln = t[: m.start()].count("\n") + 1
                    if f == p and node.lineno <= ln <= (node.end_lineno or node.lineno):
                        continue
                    if f == p and ln == node.lineno:
                        continue
                    hits.append(f"{f.relative_to(ROOT)}:{ln}")
            # also class methods? skip
            if len(hits) == 0:
                rows.append((str(p.relative_to(ROOT)), node.lineno, node.end_lineno, name, "NO REFS"))
            elif len(hits) <= 2:
                rows.append(
                    (str(p.relative_to(ROOT)), node.lineno, node.end_lineno, name, "; ".join(hits))
                )

for r in sorted(rows):
    n = (r[2] or r[1]) - r[1] + 1
    print(f"{r[0]}:{r[1]}-{r[2]} ({n}L) {r[3]}  <- {r[4]}")
