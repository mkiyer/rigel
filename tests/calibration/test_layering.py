"""The calibration package's layering, enforced.

`rigel.calibration._layers` names the layers of the import graph, and this module holds three properties of them: every
module has a declared home; no import points up a layer, the rule that gives the ordering meaning; and the declaration
matches the tree, so a layer naming a missing module fails too. Without the check the package would have no cycles but
would not say where anything new belongs, and a layering written only in prose would drift.
"""

from __future__ import annotations

import ast
import pathlib

from rigel.calibration._layers import LAYERS, layer_of

PKG = pathlib.Path(__file__).resolve().parents[2] / "src" / "rigel" / "calibration"


def _name(p: pathlib.Path) -> str:
    rel = str(p.relative_to(PKG).with_suffix(""))
    return rel.replace("/__init__", "").replace("__init__", "<pkg>")


def _runtime_imports(p: pathlib.Path) -> set[str]:
    """Sibling modules imported at runtime. Imports inside ``if TYPE_CHECKING:`` are excluded: an annotation cannot form
    a cycle and does not constrain the layering (``capture_eff_length`` annotates a type from layer 7 this way).
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
    """Adding a module means adding it to `_layers.LAYERS`, which forces the question "which layer is this?" when it is
    answerable."""
    unplaced = sorted(n for n in (_name(p) for p in ALL_FILES) if layer_of(n) is None)
    assert not unplaced, (
        f"modules with no declared layer: {unplaced}. Add each to rigel.calibration._layers.LAYERS — "
        f"see its docstring for what each layer is for."
    )


def test_the_declaration_matches_the_tree():
    """A layer naming a module that does not exist fails too."""
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


def test_no_import_points_UP_a_layer():
    """An import may point down a layer or sideways within one, never up. An upward import almost always means a type is
    defined too high (``RegionDeconv``, one slot's deconvolution result ``(f_pos, f_neg, f_g)``, is the example, reached
    for by three layers below the strand family that first defined it); the type belongs at the bottom.
    """
    bad = {}
    for path in ALL_FILES:
        me = _name(path)
        mine = layer_of(me)
        if mine is None:
            continue  # test_every_module_has_a_declared_home owns that failure
        up = sorted(
            f"{SHORT[i]} (layer {layer_of(SHORT[i])})"
            for i in _runtime_imports(path)
            if i in SHORT and SHORT[i] != me and (layer_of(SHORT[i]) or 0) > mine
        )
        if up:
            bad[f"{me} (layer {mine})"] = up
    assert not bad, (
        f"upward imports: {bad}. Either the thing needed belongs lower — a TYPE almost always does — or "
        f"the importer belongs higher. See rigel.calibration._layers."
    )


def test_the_layering_is_not_vacuous():
    """A layering with everything in one layer, or with no imports crossing layers, would pass every test above and
    constrain nothing, so this asserts several layers and real downward imports crossing them."""
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
        f"only {crossing} downward boundaries cross a layer boundary — the layering would be satisfied by "
        f"almost any assignment, so it is not constraining anything."
    )
