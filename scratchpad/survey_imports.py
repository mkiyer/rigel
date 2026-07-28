"""Survey: for every script under scripts/, check that its rigel imports resolve."""

import ast
import importlib
import pathlib
import sys

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
bad = {}
for p in sorted(ROOT.glob("scripts/**/*.py")):
    try:
        tree = ast.parse(p.read_text())
    except SyntaxError as e:
        bad[str(p)] = [f"SYNTAX: {e}"]
        continue
    probs = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if not mod.startswith("rigel"):
                continue
            try:
                m = importlib.import_module(mod)
            except Exception as e:
                probs.append(f"L{node.lineno}: import module {mod!r} -> {type(e).__name__}: {e}")
                continue
            for a in node.names:
                if a.name == "*":
                    continue
                if not hasattr(m, a.name):
                    try:
                        importlib.import_module(f"{mod}.{a.name}")
                    except Exception:
                        probs.append(f"L{node.lineno}: {mod}.{a.name} MISSING")
        elif isinstance(node, ast.Import):
            for a in node.names:
                if not a.name.startswith("rigel"):
                    continue
                try:
                    importlib.import_module(a.name)
                except Exception as e:
                    probs.append(f"L{node.lineno}: import {a.name} -> {type(e).__name__}: {e}")
    if probs:
        bad[str(p.relative_to(ROOT))] = probs

for k in sorted(bad):
    print(f"### {k}")
    for v in bad[k]:
        print(f"    {v}")
print(f"\n{len(bad)} of {len(list(ROOT.glob('scripts/**/*.py')))} scripts have unresolved rigel imports",
      file=sys.stderr)
