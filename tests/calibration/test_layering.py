"""⭐⭐⭐ THE LAYERING IS ENFORCED, NOT DOCUMENTED — because a documented layering rots silently.

The calibration package had 35 modules and no stated shape. Measured from the AST it was never a knot —
**no import cycles, and 18 of 35 modules with exactly one importer** — it was a **FLAT PILE of peers**, and a
flat pile is the one structure that cannot tell you where to add anything.

`rigel.calibration._layers` names the layers that were already in the edges. This file makes them true:

* **every module has a declared home** — an unplaced module fails, because "nobody decided where this goes"
  is exactly the state being ended;
* **no import points UP a layer** — the rule that makes the ordering mean something;
* **the declaration matches the tree** — a layer naming a module that does not exist fails too.

⛔ **Why a test and not a docstring.** The census that produced this found **13 module docstrings naming a
sibling with no import edge in either direction** — prose about the code, inside the code, that nothing
gated. A layering written only in prose would join them within a release. `scripts/design/module_census.py`
prints the same graph for a human; this file is what stops it drifting.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

from rigel.calibration._layers import LAYERS, layer_of

PKG = pathlib.Path(__file__).resolve().parents[2] / "src" / "rigel" / "calibration"


def _name(p: pathlib.Path) -> str:
    rel = str(p.relative_to(PKG).with_suffix(""))
    return rel.replace("/__init__", "").replace("__init__", "<pkg>")


def _runtime_imports(p: pathlib.Path) -> set[str]:
    """Sibling modules imported at RUNTIME.

    ⭐ Imports inside ``if TYPE_CHECKING:`` are excluded deliberately: an annotation cannot form a cycle and
    does not constrain the layering. ⚠ They are not ignorance — the census reports them, and one of them is
    how ``capture_eff_length`` annotates a type from layer 7. That is a hint, not a violation.
    """
    tree = ast.parse(p.read_text())
    guarded: set[int] = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.If) and "TYPE_CHECKING" in ast.unparse(n.test):
            for sub in ast.walk(n):
                guarded.add(id(sub))
    out: set[str] = set()
    for n in ast.walk(tree):
        if id(n) in guarded or not isinstance(n, ast.ImportFrom) or not n.level:
            continue
        if n.module:
            out.add(n.module.split(".")[-1])
        out.update(a.name for a in n.names)
    return out


ALL_FILES = sorted(PKG.rglob("*.py"))
SHORT = {_name(p).split("/")[-1]: _name(p) for p in ALL_FILES}


def test_every_module_has_a_declared_home():
    """⛔ An unplaced module is a new file nobody decided the home of — which is the flat pile returning one
    file at a time. Adding a module means adding it to `_layers.LAYERS`, and that is the point: it forces
    the question "which layer is this?" at the moment it is answerable."""
    unplaced = sorted(n for n in (_name(p) for p in ALL_FILES) if layer_of(n) is None)
    assert not unplaced, (
        f"modules with no declared layer: {unplaced}. Add each to rigel.calibration._layers.LAYERS — "
        f"see its docstring for what each layer is for."
    )


def test_the_declaration_matches_the_tree():
    """A layer naming a module that does not exist is the same rot in the other direction."""
    present = {_name(p) for p in ALL_FILES}
    declared = {m for _n, _t, members in LAYERS for m in members}
    assert not (declared - present), f"declared but absent: {sorted(declared - present)}"


def test_no_module_is_declared_twice():
    declared = [m for _n, _t, members in LAYERS for m in members]
    dupes = sorted({m for m in declared if declared.count(m) > 1})
    assert not dupes, f"a module in two layers has no home at all: {dupes}"


def test_the_layers_are_numbered_in_order():
    nums = [n for n, _t, _m in LAYERS]
    assert nums == sorted(nums) == list(range(len(nums)))


@pytest.mark.parametrize("path", ALL_FILES, ids=_name)
def test_no_import_points_UP_a_layer(path):
    """⛔⛔ **THE RULE.** An import may point DOWN a layer or SIDEWAYS within one. Never UP.

    ⭐ It found two real violations on the tree as it stood, and both were the SAME defect wearing two
    costumes: a TYPE defined too high. ``NodeDeconv`` — one slot's deconvolution result, the pie
    ``(f_pos, f_neg, f_g)`` that is the tool's central datum — was defined in the STRAND family at layer 4
    and imported by `node_geometry` and `simplex_logodds` at layer 3 and by `sweep` at layer 6. Three layers
    reached upward for it. The repair is the one a layering violation always asks for: **the type belongs at
    the bottom, not with the code that happened to define it first.** It is now layer 0.
    """
    me = _name(path)
    mine = layer_of(me)
    if mine is None:
        pytest.skip("unplaced — test_every_module_has_a_declared_home owns that failure")
    up = sorted(
        f"{SHORT[i]} (layer {layer_of(SHORT[i])})"
        for i in _runtime_imports(path)
        if i in SHORT and SHORT[i] != me and (layer_of(SHORT[i]) or 0) > mine
    )
    assert not up, (
        f"{me} is layer {mine} and imports UP into {up}. Either the thing it needs belongs lower — a TYPE "
        f"almost always does — or {me} belongs higher. See rigel.calibration._layers."
    )


def test_the_layering_is_not_vacuous():
    """⛔ TRAPS: could-the-arm-have-fired applied to this file: a layering with everything in one layer, or with no edges
    between layers, would pass every test above and constrain nothing. So assert the ordering is doing
    work — several layers, and real downward edges crossing them."""
    assert len(LAYERS) >= 5, "a layering with too few layers cannot express a direction"
    crossing = 0
    for p in ALL_FILES:
        mine = layer_of(_name(p))
        if mine is None:
            continue
        for i in _runtime_imports(p):
            t = SHORT.get(i)
            if t and layer_of(t) is not None and layer_of(t) < mine:
                crossing += 1
    assert crossing >= 20, (
        f"only {crossing} downward edges cross a layer boundary — the layering would be satisfied by "
        f"almost any assignment, so it is not constraining anything."
    )
