"""The instrument shelf's index is checked in both directions, because a hand-written one drifts.

The only index of `scripts/design/` is a table in `CLAUDE.md`, maintained by hand, and an index nobody
can trust is worse than none: a reader takes its silence as "that script does not exist" and rebuilds
the instrument. These gates make the table's two promises true — every instrument on disk has a row, so
nothing is invisible, and every row names a file that exists, so nothing points at a ghost — and they
hold the same two promises for `scripts/profiling/` against `scripts/README.md`. Every script in all
three trees must also import and carry a module docstring, since an instrument that raises on import is
still indexed and still documented. No count of instruments or of lines is carried here: re-derive one
(`TRAPS: re-record-the-baseline`). The gates deliberately do not judge whether a row's prose is
ACCURATE — a test cannot, and pretending otherwise would read as more coverage than it has
(`TRAPS: a-gate-that-reconstructs`).
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
CLAUDE = ROOT / "CLAUDE.md"
SCRIPTS = ROOT / "scripts"
DESIGN_DIR = SCRIPTS / "design"

#: rows look like ``| `design/policy_benchmark.py` | … |``
_ROW = re.compile(r"`design/([a-z0-9_]+\.py)`")

#: Instruments that predate the current campaign, were never run in it, and are recorded as drift rather
#: than silently tolerated. Each is a DECISION owed: promote it to the table or delete it. Adding to this
#: list is not a fix — it is a way of saying "not yet", and the list should only ever shrink.
UNDOCUMENTED_DEBT: frozenset[str] = frozenset()

#: the instruments — every design/ file except the package marker and the `_`-prefixed helper modules,
#: which answer no question and have no row (they are still import- and docstring-gated below).
ON_DISK = frozenset(
    p.name
    for p in DESIGN_DIR.glob("*.py")
    if p.name != "__init__.py" and not p.name.startswith("_")
)
IN_TABLE = frozenset(_ROW.findall(CLAUDE.read_text()))

SIM_DIR = SCRIPTS / "sim"

#: The third tree. `scripts/profiling/` was covered by no gate while `design/` and `sim/` were, and it
#: rotted twice with the one defect class the gates here already catch — a driver reading ``sys.argv[1]``
#: at import time, and an instrument importing two names its helper did not export, so even ``--help``
#: raised. Only the gate's REACH differed (`TRAPS: a-green-suite-hid-five-dead-instruments`), so the
#: remedy was to extend it rather than delete the tree.
PROFILING_DIR = SCRIPTS / "profiling"

#: Instruments that do not import, each with the reason and the decision owed. Same contract as
#: ``UNDOCUMENTED_DEBT``: adding to this is not a fix, it is a way of saying "not yet", and the list
#: should only ever shrink. A name here is still gated — the test asserts it fails for the RECORDED
#: reason, so a script that starts working, and one that breaks a NEW way, both fail loudly. Empty, and
#: it should stay that way: the gate refuses a stale entry as loudly as a broken script, because an
#: exemption that outlives its defect hides the next real break
#: (`TRAPS: a-green-suite-hid-five-dead-instruments`).
BROKEN_ON_IMPORT: dict[str, str] = {}


def _instruments(directory: pathlib.Path) -> list[pathlib.Path]:
    return [p for p in directory.glob("*.py") if p.name != "__init__.py"]


#: every script the import + docstring gates cover, all THREE trees.
ALL_SCRIPTS = sorted(
    _instruments(DESIGN_DIR) + _instruments(SIM_DIR) + _instruments(PROFILING_DIR),
    key=lambda p: (p.parent.name, p.name),
)


def _case_id(path: pathlib.Path) -> str:
    """Name the tree, never the basename alone: two trees can carry one basename, and under a basename
    id pytest would de-duplicate them to `name.py0` / `name.py1`, so a failure would name neither file
    and `-k` could not select one."""
    return f"{path.parent.name}/{path.name}"


@pytest.mark.parametrize("path", ALL_SCRIPTS, ids=_case_id)
def test_every_instrument_still_imports(path):
    """A `src/` deletion kills instruments silently, and nothing else here can see it.

    Being indexed and having a docstring are both true of a script that raises on line 1, so importing
    is the cheapest check that an instrument is still connected to the code it measures — and it is the
    one that rots, since a deletion in `src/` is green everywhere else
    (`TRAPS: a-green-suite-hid-five-dead-instruments`).

    Import only, never execution: an instrument's numbers need its substrate, so this is a connectivity
    gate rather than a claim that the script is CORRECT (`TRAPS: a-gate-that-reconstructs`).
    """
    import importlib.util
    import sys

    expected = BROKEN_ON_IMPORT.get(path.name)
    name = f"_gate_{path.stem}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)

    # Reproduce what `python scripts/design/x.py` does, or the gate tests a different program. Two
    # things the naive spelling gets wrong, both of which read as false failures on instruments that run
    # perfectly (`TRAPS: a-gate-that-reconstructs` — a gate that rebuilds its subject tests the rebuild):
    #   * the module must be in `sys.modules` BEFORE `exec_module`, or a dataclass resolving its own
    #     `__module__` gets `None` and raises `'NoneType' object has no attribute '__dict__'`;
    #   * the script's OWN directory is `sys.path[0]` under a real invocation, which is how the
    #     sibling-importing instruments (`_sibling`) resolve.
    sys.modules[name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    except SystemExit:
        pass  # a script that argues with argv at import is still connected
    except Exception as exc:  # noqa: BLE001
        if expected and expected in str(exc):
            pytest.xfail(
                f"{path.name}: known broken on `{expected}` — a DECISION is owed, see the dict"
            )
        raise AssertionError(
            f"{path.name} no longer imports: {type(exc).__name__}: {exc}\n"
            f"An instrument that cannot be imported cannot be run, and nothing else in this file "
            f"would have caught it."
        ) from exc
    else:
        if expected:
            raise AssertionError(
                f"{path.name} imports fine now, but is still listed in BROKEN_ON_IMPORT as broken on "
                f"`{expected}`. Remove the entry — a stale exemption hides the next real break."
            )
    finally:
        sys.path.remove(str(path.parent))
        sys.modules.pop(name, None)


def test_every_instrument_is_in_the_index_or_named_as_debt():
    """A script nobody indexed is a script the next session rebuilds."""
    missing = sorted(ON_DISK - IN_TABLE - UNDOCUMENTED_DEBT)
    assert not missing, (
        f"{len(missing)} instruments are not in CLAUDE.md's table and are not listed as known debt: "
        f"{missing}. Add a row describing WHAT QUESTION IT ANSWERS, or delete the file."
    )


def test_the_index_does_not_point_at_ghosts():
    """The same rot in the other direction — a row for a file that was deleted."""
    ghosts = sorted(IN_TABLE - ON_DISK)
    assert not ghosts, (
        f"CLAUDE.md's table has rows for files that do not exist: {ghosts}. A reader who trusts the table "
        f"goes looking for an instrument that was deleted."
    )


def test_the_profiling_tree_is_indexed_in_the_scripts_readme():
    """The third tree has its own index, and it is `scripts/README.md` rather than `CLAUDE.md`.

    `CLAUDE.md`'s table indexes `scripts/design/`; the profiling drivers are indexed in the README's
    `profiling/` row instead, so this checks the same two promises there — nothing on disk is invisible,
    and no row points at a ghost. One test rather than two parametrised ones, because the tree is small
    and the failure message can name both directions at once.
    """
    readme = (SCRIPTS / "README.md").read_text()
    # Tree-qualified, like `CLAUDE.md`'s `design/…` rows: a bare basename would let a README mention of
    # a `design/` file satisfy the row for a different instrument of the same name.
    listed = frozenset(re.findall(r"`profiling/([a-z0-9_]+\.py)`", readme))
    on_disk = frozenset(p.name for p in _instruments(PROFILING_DIR))
    assert on_disk, (
        "scripts/profiling/ is empty — delete the row and this gate, or restore the tree"
    )
    missing = sorted(on_disk - listed)
    ghosts = sorted(listed - on_disk)
    assert not missing, (
        f"profiling instruments not named in scripts/README.md: {missing}. A reader takes the README's "
        f"silence as 'that script does not exist' and rebuilds it."
    )
    assert not ghosts, f"scripts/README.md names profiling files that do not exist: {ghosts}"


def test_the_documented_debt_is_real_debt():
    """The debt list may only name files that EXIST and are NOT in the table. A stale entry there would
    let a real gap hide behind a name that has already been dealt with."""
    gone = sorted(f for f in UNDOCUMENTED_DEBT if f not in ON_DISK)
    assert not gone, f"debt entries for deleted files — remove them: {gone}"
    resolved = sorted(f for f in UNDOCUMENTED_DEBT if f in IN_TABLE)
    assert not resolved, (
        f"these are now IN the table, so they are no longer debt — remove them from UNDOCUMENTED_DEBT: "
        f"{resolved}"
    )


@pytest.mark.parametrize("path", ALL_SCRIPTS, ids=_case_id)
def test_every_instrument_says_what_it_is_for(path):
    """A module docstring is the only thing a reader has before running an instrument. One with none is
    unusable without reading its argument parser. All three trees, the same reach as the import gate.
    """
    doc = ast.get_docstring(ast.parse(path.read_text()))
    assert doc and len(doc.strip()) > 60, (
        f"scripts/{_case_id(path)} has no usable module docstring. Lead with the QUESTION it answers — "
        f"that is what makes the shelf navigable."
    )


def test_the_index_is_not_vacuous():
    """`TRAPS: could-the-arm-have-fired` applied here: if the row regex matched nothing, every test
    above would pass while checking nothing at all."""
    assert IN_TABLE, "no rows parsed from CLAUDE.md — the table format moved"
    assert ON_DISK, "no instruments found — DESIGN_DIR moved"
