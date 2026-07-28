"""Survey 4: class methods + dataclass fields in src/rigel/*.py (excl. calibration) with no refs."""

import ast
import pathlib
import re

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
files = sorted(ROOT.glob("src/rigel/*.py"))
corpus_files = (
    list(ROOT.glob("src/rigel/**/*.py"))
    + list(ROOT.glob("tests/**/*.py"))
    + list(ROOT.glob("scripts/**/*.py"))
    + list(ROOT.glob("src/rigel/native/**/*.cpp"))
    + list(ROOT.glob("src/rigel/native/**/*.h"))
    + list(ROOT.glob("docs/**/*.md"))
)
texts = {}
for f in corpus_files:
    try:
        texts[f] = f.read_text()
    except Exception:
        pass

DUNDER_OK = {"__post_init__", "__init__", "__repr__", "__len__", "__iter__", "__enter__",
             "__exit__", "__eq__", "__hash__", "__getitem__", "__setitem__", "__contains__"}


def refs(name, owner, lo, hi):
    pat = re.compile(r"\b" + re.escape(name) + r"\b")
    hits = []
    for f, t in texts.items():
        for m in pat.finditer(t):
            ln = t[: m.start()].count("\n") + 1
            if f == owner and lo <= ln <= hi:
                continue
            hits.append(f"{f.relative_to(ROOT)}:{ln}")
    return hits


for p in files:
    tree = ast.parse(p.read_text())
    for cls in tree.body:
        if not isinstance(cls, ast.ClassDef):
            continue
        for node in cls.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if node.name in DUNDER_OK or node.name.startswith("__"):
                    continue
                h = refs(node.name, p, node.lineno, node.end_lineno or node.lineno)
                if not h:
                    n = (node.end_lineno or node.lineno) - node.lineno + 1
                    print(f"DEAD-METHOD {p.relative_to(ROOT)}:{node.lineno}-{node.end_lineno} "
                          f"({n}L) {cls.name}.{node.name}")
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                nm = node.target.id
                if nm.startswith("_"):
                    continue
                h = refs(nm, p, node.lineno, node.end_lineno or node.lineno)
                # exclude self-file constructor refs
                if not h:
                    print(f"DEAD-FIELD {p.relative_to(ROOT)}:{node.lineno} {cls.name}.{nm}")
