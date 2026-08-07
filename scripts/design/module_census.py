#!/usr/bin/env python
"""⭐⭐⭐ **WHERE DOES A CHANGE GO?** — the calibration package's real shape, re-derived from the AST.

⛔ **THE PROBLEM THIS EXISTS FOR, stated as the owner stated it: 35 modules and no sense of how or where to
develop.** And the graph says why, which is not what it looks like from the outside: there are **no import
cycles**, and **18 of the 35 modules have exactly one importer**. It is not a knot. It is a FLAT PILE of 35
peers, and nothing in the tree names the layers that already exist in the edges.

⭐⭐ **The second finding is the one that actually blocks a reader, and this instrument is the only thing
that can see it: the module docstrings MISDESCRIBE THE GRAPH.** ``run_fill`` said it was "shared by
`density_model`, `strand_deconv`, `priors`, and the `sweep` chain geometry" and had **one** importer.
``strand_likelihood`` said it was "Used by the per-node strand module (`strand_deconv`)" and
``strand_deconv`` does not import it at all. A developer who trusts either sentence goes looking for code
that is not there — and neither sentence can rot *loudly*, because nothing checks prose against edges.
⛔ That is the same failure as a stale doc citation, one layer down: **a claim about the code, inside the
code, that nothing gates.** Measured 14 on 2026-08-07; **6 were genuinely stale and are fixed**, and the
rest are data-flow statements the instrument cannot tell apart — see below.

What it reports
---------------
1. **THE LAYERING**, from ``rigel.calibration._layers``, with every module's assigned layer and any import
   that points UPWARD (which is what a layering violation is). ``test_layering.py`` gates this; here it is
   printed so a reader can see the shape rather than only be told it holds.
2. **THE GRAPH** — per module: layer, lines, public surface, and every importer inside the package and
   outside it. ⭐ A module with no importer inside the package is an ENTRY POINT or it is DEAD, and the two
   are distinguished by whether anything outside the package calls it.
3. ⭐⭐ **STALE DOCSTRING CROSS-REFERENCES** — every sibling module a docstring names for which **no import
   edge exists in either direction**. Each one is a sentence pointing at code that is not connected.
4. **DEAD PUBLIC SURFACE** — exported names that nothing anywhere imports. ⚠ Not automatically a cut: a
   name may be an executable reference that a test gates, which is a legitimate second home. The report
   says which, by naming the importers it does have.

⚠ **It reports, it does not judge.** A verdict needs a human, because three of the categories above have a
legitimate member: an entry point looks dead, a reference implementation looks duplicated, and a
single-consumer module may be a named concept worth its own file.

Usage
-----
    python scripts/design/module_census.py                 # the calibration package
    python scripts/design/module_census.py --package rigel # any package under src/
    python scripts/design/module_census.py --stale-only    # just the lying docstrings
"""

from __future__ import annotations

import argparse
import ast
import collections
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

#: a sibling module named inside a docstring or comment, as ``name`` / `name` / :mod:`name`
_REF = re.compile(r"[`:](?:mod:`)?([a-z_][a-z0-9_]{2,})`")


def _modules(pkg: pathlib.Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for p in sorted(pkg.rglob("*.py")):
        rel = p.relative_to(pkg)
        name = str(rel.with_suffix("")).replace("/__init__", "").replace("__init__", "<pkg>")
        src = p.read_text()
        tree = ast.parse(src)
        imports: set[str] = set()
        typing_only: set[str] = set()
        # ⭐ An import inside ``if TYPE_CHECKING:`` is an ANNOTATION, not an edge — it cannot form a cycle
        # and it does not constrain the layering. It is still reported, because an annotation reaching
        # upward is a hint that a TYPE belongs lower, which is exactly what it was here.
        guarded: set[int] = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.If):
                test = ast.unparse(n.test)
                if "TYPE_CHECKING" in test:
                    for sub in ast.walk(n):
                        guarded.add(id(sub))
        for n in ast.walk(tree):
            tgt = typing_only if id(n) in guarded else imports
            if isinstance(n, ast.ImportFrom) and n.level:
                if n.module:
                    tgt.add(n.module.split(".")[-1])
                for a in n.names:
                    tgt.add(a.name)
            elif isinstance(n, ast.Import):
                for a in n.names:
                    tgt.add(a.name.split(".")[-1])
        out[name] = {
            "path": p,
            "loc": len(src.splitlines()),
            "imports": imports,
            "typing_only": typing_only,
            "public": [
                n.name
                for n in tree.body
                if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and not n.name.startswith("_")
            ],
            "doc": ast.get_docstring(tree) or "",
        }
    return out


def _edges(mods: dict[str, dict]) -> tuple[dict, dict]:
    """``(uses, used_by)`` over intra-package edges, keyed by the module's own short name."""
    short = {n.split("/")[-1]: n for n in mods}
    uses = collections.defaultdict(set)
    used_by = collections.defaultdict(set)
    for m, d in mods.items():
        for i in d["imports"]:
            t = short.get(i)
            if t and t != m:
                uses[m].add(t)
                used_by[t].add(m)
    return uses, used_by


def _outside(pkg: pathlib.Path, mods: dict[str, dict]) -> dict[str, set[str]]:
    """Who imports each module from OUTSIDE the package — what separates an entry point from dead code."""
    outside: dict[str, set[str]] = {m: set() for m in mods}
    pkg_rel = pkg.relative_to(ROOT / "src") if (ROOT / "src") in pkg.parents else pkg
    dotted = str(pkg_rel).replace("/", ".")
    for base in ("src", "tests", "scripts"):
        for p in (ROOT / base).rglob("*.py"):
            if pkg in p.parents:
                continue
            txt = p.read_text(errors="ignore")
            for m in mods:
                short = m.split("/")[-1]
                if short == "<pkg>":
                    continue
                if (
                    f"{dotted}.{short}" in txt
                    or f"{dotted} import {short}" in txt
                    or re.search(rf"\bimport\b[^\n]*\b{short}\b", txt)
                    and short in txt
                ):
                    # a name-only match is weak, so require the package or an explicit from-import
                    if f"{dotted}.{short}" in txt or f"{dotted} import" in txt:
                        outside[m].add(str(p.relative_to(ROOT)))
    return outside


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--package", default="rigel/calibration")
    ap.add_argument("--stale-only", action="store_true")
    args = ap.parse_args()

    pkg = ROOT / "src" / args.package
    if not pkg.is_dir():
        raise SystemExit(f"⛔ no such package: {pkg}")
    mods = _modules(pkg)
    uses, used_by = _edges(mods)
    outside = _outside(pkg, mods)
    short = {n.split("/")[-1]: n for n in mods}

    try:
        from rigel.calibration._layers import LAYERS, layer_of
    except ImportError:
        LAYERS, layer_of = None, lambda _m: None  # noqa: E731

    total = sum(m["loc"] for m in mods.values())
    print()
    print(f"   {args.package}: {len(mods)} modules, {total:,} lines")

    # ── 1. THE LAYERING ────────────────────────────────────────────────────────────────────────────────
    if LAYERS is not None and not args.stale_only:
        print()
        print("   ⭐⭐ THE LAYERING — where a change goes. An import may point DOWN or sideways, never UP.")
        viol = []
        for i, (num, title, members) in enumerate(LAYERS):
            present = [m for m in members if m in mods]
            missing = [m for m in members if m not in mods]
            loc = sum(mods[m]["loc"] for m in present)
            print(f"\n   {num}. {title:<34} {loc:>6,} lines")
            for m in present:
                up = sorted(t for t in uses[m] if (layer_of(t) or 0) > (layer_of(m) or 0))
                ann = sorted(
                    short[x]
                    for x in mods[m]["typing_only"]
                    if x in short and (layer_of(short[x]) or 0) > (layer_of(m) or 0)
                )
                mark = "  ⛔ imports UP: " + ", ".join(up) if up else ""
                if ann and not up:
                    mark = "  ⚠ annotates UP (TYPE_CHECKING, not an edge): " + ", ".join(ann)
                print(f"        {m:<24} {mods[m]['loc']:>5}{mark}")
                viol += [(m, t) for t in up]
            for m in missing:
                print(f"        ⚠ {m} is declared in the layering and does not exist")
        unplaced = [m for m in mods if layer_of(m) is None]
        if unplaced:
            print(f"\n   ⛔ NOT IN ANY LAYER — a new module nobody placed: {', '.join(sorted(unplaced))}")
        print()
        print(f"   {'✅ no upward imports' if not viol else f'⛔ {len(viol)} UPWARD IMPORTS'}")

    # ── 3. STALE DOCSTRING CROSS-REFERENCES ────────────────────────────────────────────────────────────
    print()
    print("   ⭐⭐ UNBACKED DOCSTRING CROSS-REFERENCES — a sibling named in prose with NO import edge either way")
    print("      ⚠ NOT automatically wrong, and this is the distinction the instrument CANNOT make:")
    print("        an IMPORT claim (\"shared by X\", \"used by X\") with no edge is STALE and must be fixed;")
    print("        a DATA-FLOW claim (\"fitted in X\", \"X consumes this\") is true of a value passed through")
    print("        a caller and is worth keeping. Read the sentence; the count is a worklist, not a verdict.")
    n_stale = 0
    for m in sorted(mods):
        named = {short[t] for t in _REF.findall(mods[m]["doc"]) if t in short and short[t] != m}
        stale = sorted(named - uses[m] - used_by[m])
        if stale:
            n_stale += len(stale)
            print(f"      {m:<24} names {', '.join(stale)}")
    print(f"      {'✅ none' if not n_stale else f'{n_stale} to read'}"
          f"{'' if not n_stale else '  (was 14 on 2026-08-07; 6 were genuinely stale and are fixed)'}")
    if args.stale_only:
        return 1 if n_stale else 0

    # ── 2. THE GRAPH ───────────────────────────────────────────────────────────────────────────────────
    print()
    print("   THE GRAPH — a module with no importer INSIDE is an entry point, or it is dead")
    print(f"   {'module':<24} {'L':>2} {'LOC':>5} {'pub':>4} {'in':>3} {'out':>4}  importers inside")
    print("   " + "-" * 108)
    for m in sorted(mods, key=lambda x: (layer_of(x) or 99, x)):
        ins, outs = sorted(used_by[m]), outside[m]
        flag = ""
        if not ins and not outs and m != "<pkg>":
            flag = "  ⛔ DEAD: nothing imports it, anywhere"
        elif not ins and m != "<pkg>":
            flag = "  ⭐ ENTRY POINT (used only from outside)"
        print(
            f"   {m:<24} {layer_of(m) or '?':>2} {mods[m]['loc']:>5} {len(mods[m]['public']):>4} "
            f"{len(ins):>3} {len(outs):>4}  {', '.join(ins)[:44]}{flag}"
        )

    # ── 4. DEAD PUBLIC SURFACE ─────────────────────────────────────────────────────────────────────────
    print()
    print("   DEAD PUBLIC SURFACE — exported and imported by nothing")
    print("      ⚠ NOT automatically a cut: an executable REFERENCE that a test gates is a legitimate")
    print("        second home. The count of test-only importers is what tells the two apart.")
    all_src = "\n".join(
        p.read_text(errors="ignore")
        for base in ("src", "tests", "scripts")
        for p in (ROOT / base).rglob("*.py")
    )
    n_dead = 0
    for m in sorted(mods):
        dead = [
            nm
            for nm in mods[m]["public"]
            if len(re.findall(rf"\b{re.escape(nm)}\b", all_src)) <= 1
        ]
        if dead:
            n_dead += len(dead)
            print(f"      {m:<24} {', '.join(dead)[:78]}")
    print(f"      {'✅ none' if not n_dead else f'{n_dead} names defined and never imported'}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
