"""The calibration package's layering is enforced here rather than merely written down.

`rigel.calibration._layers` names the layers that are already in the import graph, and this file holds
three properties true of them: every module has a declared home, because an unplaced module is one
nobody decided the place of; no import points UP a layer, which is the rule that makes the ordering
mean anything; and the declaration matches the tree, so a layer naming a module that does not exist
fails as well. Without the enforcement the package is a flat pile of peers rather than a knot
(TRAPS: a-flat-pile-is-not-a-knot) — it has no cycles, but neither does it tell you where to add
anything. A layering written only in prose would rot the way module docstrings naming a sibling they
never import already do, and this file is what stops it drifting.
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

    Imports inside ``if TYPE_CHECKING:`` are excluded deliberately: an annotation cannot form a cycle
    and does not constrain the layering. One of them is how ``capture_eff_length`` annotates a type from layer 7. That is a hint, not a
    violation.
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
    """An unplaced module is a new file nobody decided the home of — the flat pile returning one file at
    a time. Adding a module means adding it to `_layers.LAYERS`, and that is the point: it forces the
    question "which layer is this?" at the moment it is answerable."""
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
    """⛔ The rule: an import may point DOWN a layer or SIDEWAYS within one, never UP.

    An upward import is almost always a TYPE defined too high — ``RegionDeconv``, one slot's
    deconvolution result ``(f_pos, f_neg, f_g)``, is the example, reached for by three layers below the
    strand family that first defined it. The repair a layering violation asks for is the same every
    time: the type belongs at the bottom, not with the code that happened to define it first.
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
    """TRAPS: could-the-arm-have-fired applied to this file: a layering with everything in one layer, or
    with no imports crossing between layers, would pass every test above and constrain nothing. So
    assert the ordering is doing work — several layers, and real downward imports crossing them."""
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
