"""Survey 2: attribute access on imported rigel modules that no longer exists."""

import ast
import importlib
import pathlib

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
report = {}
for p in sorted(ROOT.glob("scripts/**/*.py")) + sorted(ROOT.glob("src/rigel/**/*.py")):
    try:
        tree = ast.parse(p.read_text())
    except SyntaxError:
        continue
    alias = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith("rigel"):
                    alias[a.asname or a.name.split(".")[0]] = a.name
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod.startswith("rigel"):
                for a in node.names:
                    full = f"{mod}.{a.name}"
                    try:
                        importlib.import_module(full)
                    except Exception:
                        continue
                    alias[a.asname or a.name] = full
    probs = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            base = node.value.id
            if base in alias:
                try:
                    m = importlib.import_module(alias[base])
                except Exception:
                    continue
                if not hasattr(m, node.attr):
                    probs.append(f"L{node.lineno}: {alias[base]}.{node.attr} MISSING")
    if probs:
        report[str(p.relative_to(ROOT))] = sorted(set(probs))

for k in sorted(report):
    print(f"### {k}")
    for v in report[k]:
        print(f"    {v}")
