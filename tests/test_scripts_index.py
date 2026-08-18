"""⭐⭐ THE INSTRUMENT SHELF HAS AN INDEX, AND THE INDEX IS CHECKED — because a hand-written one drifts.

`scripts/design/` is **59 instruments** (re-derived 2026-08-17, after the two dark `flgap` instruments
were deleted), and it is the project's debug loop rather than sprawl. The problem was never the count. It
was that the only index — a table in `CLAUDE.md` — is hand-maintained, and by the time it was measured it
had **drifted by eight entries**, which `CLAUDE.md` itself had to admit in prose.

⚠ Even that one figure is hand-carried prose, so re-derive rather than trust
(`TRAPS: re-record-the-baseline`): `ls scripts/design/*.py | wc -l`. It is at least *constrained* — the
two gates below force `rows(CLAUDE.md) + UNDOCUMENTED_DEBT == on disk`, so a wrong count here is visible
from inside this file. ⛔ **A LINE COUNT IS NOT, AND ONE WAS CARRIED HERE UNTIL IT WENT STALE TWICE.**
The predecessor read "49 instruments, 15,904 boundaries" and its successor "59 instruments, 25,180 lines",
which was already wrong by five within the session that measured it — every prose edit to any instrument
moves it and nothing consumes it, so it has been dropped rather than re-measured. ⚠ "boundaries" in the
first of those was a bulk rename landing on the word *lines* in its ordinary-English sense, which is the
two-senses hazard `rename_census.py` exists to flag.

⛔ **An index nobody can trust is worse than none**, because a reader takes its silence as "that script
does not exist" and rebuilds it. This file makes the table's two promises true:

* every instrument on disk appears in the table, so nothing is invisible;
* every row of the table names a file that exists, so nothing points at a ghost.

⭐ Same shape as `tests/calibration/test_layering.py`: state the structure once, then let a test hold it,
rather than a convention everyone means to keep.

⚠ **It deliberately does NOT check the description.** Whether a row's prose is accurate is a judgement a
test cannot make, and pretending otherwise would be a gate that reads as coverage and is not
(`TRAPS: a-gate-that-reconstructs`). What it checks is presence in both directions, which is exactly the
part that rots silently.
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

#: rows look like ``| `design/toy_panel.py` | … |``
_ROW = re.compile(r"`design/([a-z0-9_]+\.py)`")

#: ⛔ Instruments that predate the current campaign, were never run in it, and are documented as drift
#: rather than silently tolerated. ⚠ Each is a DECISION owed: promote it to the table or delete it. Adding
#: to this list is not a fix — it is a way of saying "not yet", and the list should only ever shrink.
UNDOCUMENTED_DEBT: frozenset[str] = frozenset(
    {
        "bam_spans.py",
        "implicit_splice_census.py",
        "inv_L_limits.py",
        "length_sieve.py",
        "partition_test.py",
        "reference_on_real_data.py",
        "sigma_inv_L.py",
        "spanning.py",
    }
)

ON_DISK = frozenset(p.name for p in DESIGN_DIR.glob("*.py") if p.name != "__init__.py")
IN_TABLE = frozenset(_ROW.findall(CLAUDE.read_text()))

SIM_DIR = SCRIPTS / "sim"

#: ⛔⛔ **THE THIRD TREE, GATED SINCE 2026-08-17 — AND THE REASON IS TWO ROTS, NOT A PREFERENCE.**
#: `scripts/profiling/` was covered by NO gate while `design/` and `sim/` were, and it rotted **twice**
#: with the one defect class the gate here already catches: `pyspy_driver.py` read ``sys.argv[1]`` at
#: import time (found 2026-08-11) and `scan_profile.py` imported two names `profiler.py` did not export,
#: so even ``--help`` raised (found 2026-08-17). Both times the defect was identical to one this file
#: had been catching for months in the next directory along, and only its REACH differed
#: (`TRAPS: a-green-suite-hid-five-dead-instruments`).
#: ⭐ The owner's decision, recorded in `scripts/README.md`: extend the gate rather than delete the tree,
#: because the performance work before 0.8.0 needs it.
PROFILING_DIR = SCRIPTS / "profiling"

#: ⛔⛔ **INSTRUMENTS THAT DO NOT IMPORT, EACH WITH THE REASON AND THE DECISION OWED.** Same contract as
#: ``UNDOCUMENTED_DEBT``: adding to this is not a fix, it is a way of saying "not yet", and the list
#: should only ever SHRINK. ⚠ A name here is still gated — the test asserts it fails for the RECORDED
#: reason, so a script that starts working, or breaks a NEW way, both fail loudly.
#: ⭐ **EMPTY, AND IT SHOULD STAY THAT WAY.** An entry here is a script the gate KNOWS is dead; it buys
#: time to repair one, and nothing else. ⛔ The gate refuses a STALE entry as loudly as a broken script,
#: because an exemption that outlives its defect hides the next real break — which is exactly how five
#: dead instruments once sat behind a green suite (`TRAPS: a-green-suite-hid-five-dead-instruments`).
#: ⚠ `prior_units_check.py` was the last entry, exempted on the deleted `_component_region_arrays`, and
#: it was repaired on 2026-08-17 — at which point THIS gate fired on the stale exemption, which is the
#: behaviour it was written for.
BROKEN_ON_IMPORT: dict[str, str] = {}

def _instruments(directory: pathlib.Path) -> list[pathlib.Path]:
    return [p for p in directory.glob("*.py") if p.name != "__init__.py"]


#: every script the import + docstring gates cover, all THREE trees.
ALL_SCRIPTS = sorted(
    _instruments(DESIGN_DIR) + _instruments(SIM_DIR) + _instruments(PROFILING_DIR),
    key=lambda p: (p.parent.name, p.name),
)


def _case_id(path: pathlib.Path) -> str:
    """⛔⛔ **NAME THE TREE, NEVER THE BASENAME ALONE — A BASENAME COLLIDES ACROSS TREES.**

    `design/scan_profile.py` (the accumulator's ns/fragment, regressed over several BAMs) and
    `profiling/scan_profile.py` (the scan's own wall time and RSS across thread budgets) are different
    instruments that share a filename. Under a basename id pytest would silently de-duplicate them to
    `scan_profile.py0` / `scan_profile.py1`, whose ORDER is an implementation detail — so a failure
    would name neither file, and `-k` could not select one. ⚠ The collision itself is not fixed here;
    renaming a file is an owner call, and `scripts/README.md` records the proposal.
    """
    return f"{path.parent.name}/{path.name}"


@pytest.mark.parametrize("path", ALL_SCRIPTS, ids=_case_id)
def test_every_instrument_still_imports(path):
    """⛔⛔ **A `src/` DELETION KILLS INSTRUMENTS SILENTLY, AND NOTHING HERE COULD SEE IT.**

    This file's own docstring claimed "measured 2026-08-07, **every one of them imports cleanly**" — a
    measurement taken once, by hand, and never gated. By 2026-08-11 it was false for **five** scripts and
    the suite was green the whole time: three died when the fixed-point layer (`INV_LENGTH_SCALE`,
    `inv_length_quantum`) went at `94d283c0`, one when `enrichment_frame` went at `0d9d422b`, one on
    `_component_region_arrays`. Two commits, five dead instruments, 3,235 passing tests.

    ⭐ The other tests here check that a script is INDEXED and has a DOCSTRING — both true of a script
    that raises on boundary 1. Importing is the cheapest possible check that it is still connected to the
    code it measures, and it is the one that rots.

    ⚠ Import only, never execution: an instrument's numbers need its substrate, and this is a
    connectivity gate, not a claim that the script is CORRECT. `TRAPS: a-gate-that-reconstructs` —
    a gate that reads as more coverage than it has is worse than none.
    """
    import importlib.util
    import sys

    expected = BROKEN_ON_IMPORT.get(path.name)
    name = f"_gate_{path.stem}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)

    # ⛔⛔ REPRODUCE WHAT `python scripts/design/x.py` DOES, OR THE GATE TESTS A DIFFERENT PROGRAM.
    # Two things the naive spelling gets wrong, and both were measured as FALSE FAILURES on six
    # instruments that run perfectly (`TRAPS: a-gate-that-reconstructs` — a gate that rebuilds its
    # subject tests the rebuild):
    #   * the module must be in `sys.modules` BEFORE `exec_module`, or a dataclass resolving its own
    #     `__module__` gets `None` and raises `'NoneType' object has no attribute '__dict__'`;
    #   * the script's OWN directory is `sys.path[0]` under a real invocation, which is how the
    #     sibling-importing instruments (`_sibling`, `from reference_on_real_data import ...`) resolve.
    sys.modules[name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    except SystemExit:
        pass  # a script that argues with argv at import is still connected
    except Exception as exc:  # noqa: BLE001
        if expected and expected in str(exc):
            pytest.xfail(f"{path.name}: known broken on `{expected}` — a DECISION is owed, see the dict")
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
    """⛔ A script nobody indexed is a script the next session rebuilds."""
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
    """⛔ **THE THIRD TREE HAS ITS OWN INDEX, AND IT IS `scripts/README.md` RATHER THAN `CLAUDE.md`.**

    `CLAUDE.md`'s table indexes `scripts/design/`; the profiling drivers are indexed in the README's
    `profiling/` row instead, so this checks the same two promises there — nothing on disk is invisible,
    and no row points at a ghost. ⚠ One test rather than two parametrised ones, because the tree is
    small and the failure message can name both directions at once.
    """
    readme = (SCRIPTS / "README.md").read_text()
    # ⚠ TREE-QUALIFIED, like `CLAUDE.md`'s `design/…` rows: a bare basename would let the README's
    # mention of `design/scan_profile.py` satisfy the row for a different instrument of the same name.
    listed = frozenset(re.findall(r"`profiling/([a-z0-9_]+\.py)`", readme))
    on_disk = frozenset(p.name for p in _instruments(PROFILING_DIR))
    assert on_disk, "scripts/profiling/ is empty — delete the row and this gate, or restore the tree"
    missing = sorted(on_disk - listed)
    ghosts = sorted(listed - on_disk)
    assert not missing, (
        f"profiling instruments not named in scripts/README.md: {missing}. A reader takes the README's "
        f"silence as 'that script does not exist' and rebuilds it."
    )
    assert not ghosts, f"scripts/README.md names profiling files that do not exist: {ghosts}"


def test_the_documented_debt_is_real_debt():
    """⚠ The debt list may only name files that EXIST and are NOT in the table. A stale entry there would
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
    """⛔ A module docstring is the only thing a reader has before running it. An instrument with none is
    unusable without reading its argument parser.

    ⚠ It covers all three trees, where it used to cover `design/` alone — the same reach the import gate
    was given on 2026-08-17, and for the same reason.
    """
    doc = ast.get_docstring(ast.parse(path.read_text()))
    assert doc and len(doc.strip()) > 60, (
        f"scripts/{_case_id(path)} has no usable module docstring. Lead with the QUESTION it answers — "
        f"that is what makes the shelf navigable."
    )


def test_the_index_is_not_vacuous():
    """⚠ `TRAPS: could-the-arm-have-fired` applied here: if the row regex matched nothing, every test above
    would pass while checking nothing at all."""
    assert len(IN_TABLE) >= 30, (
        f"only {len(IN_TABLE)} rows parsed from CLAUDE.md — the table format moved"
    )
    assert len(ON_DISK) >= 30, f"only {len(ON_DISK)} instruments found — DESIGN_DIR moved"
