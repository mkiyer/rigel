"""⭐⭐ THE INSTRUMENT SHELF HAS AN INDEX, AND THE INDEX IS CHECKED — because a hand-written one drifts.

`scripts/design/` is **49 instruments, 15,904 lines**, and it is the project's debug loop rather than
sprawl: measured 2026-08-07, **every one of them imports cleanly** and 35 of 49 lead with the question they
answer. The problem was never the count. It was that the only index — a table in `CLAUDE.md` — is
hand-maintained, and by the time it was measured it had **drifted by eight entries**, which `CLAUDE.md`
itself had to admit in prose.

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
DESIGN_DIR = ROOT / "scripts" / "design"

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

SIM_DIR = ROOT / "scripts" / "sim"

#: ⛔⛔ **INSTRUMENTS THAT DO NOT IMPORT, EACH WITH THE REASON AND THE DECISION OWED.** Same contract as
#: ``UNDOCUMENTED_DEBT``: adding to this is not a fix, it is a way of saying "not yet", and the list
#: should only ever SHRINK. ⚠ A name here is still gated — the test asserts it fails for the RECORDED
#: reason, so a script that starts working, or breaks a NEW way, both fail loudly.
BROKEN_ON_IMPORT: dict[str, str] = {
    # `_component_region_arrays` was deleted from `calibration/priors.py` and has no successor: the module
    # now exposes `_region_locus_shares` / `_edge_locus_shares` / `_project_regions_to_loci`. ⭐ The
    # DECISION is not a repair — this instrument compares the EM prior "in fragment units vs the OLD
    # incidence sum", and `ROADMAP.md` §0 records that migration as DONE, so its question is closed.
    # Delete it or re-point it at the surviving API; do not leave it here.
    "prior_units_check.py": "_component_region_arrays",
}

#: every script the import gate covers, both trees.
ALL_SCRIPTS = sorted(
    [p for p in DESIGN_DIR.glob("*.py") if p.name != "__init__.py"]
    + [p for p in SIM_DIR.glob("*.py") if p.name != "__init__.py"],
    key=lambda p: p.name,
)


@pytest.mark.parametrize("path", ALL_SCRIPTS, ids=lambda p: p.name)
def test_every_instrument_still_imports(path):
    """⛔⛔ **A `src/` DELETION KILLS INSTRUMENTS SILENTLY, AND NOTHING HERE COULD SEE IT.**

    This file's own docstring claimed "measured 2026-08-07, **every one of them imports cleanly**" — a
    measurement taken once, by hand, and never gated. By 2026-08-11 it was false for **five** scripts and
    the suite was green the whole time: three died when the fixed-point layer (`INV_LENGTH_SCALE`,
    `inv_length_quantum`) went at `94d283c0`, one when `enrichment_frame` went at `0d9d422b`, one on
    `_component_region_arrays`. Two commits, five dead instruments, 3,235 passing tests.

    ⭐ The other tests here check that a script is INDEXED and has a DOCSTRING — both true of a script
    that raises on line 1. Importing is the cheapest possible check that it is still connected to the
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


@pytest.mark.parametrize("name", sorted(ON_DISK))
def test_every_instrument_says_what_it_is_for(name):
    """⛔ A module docstring is the only thing a reader has before running it. An instrument with none is
    unusable without reading its argument parser."""
    doc = ast.get_docstring(ast.parse((DESIGN_DIR / name).read_text()))
    assert doc and len(doc.strip()) > 60, (
        f"scripts/design/{name} has no usable module docstring. Lead with the QUESTION it answers — "
        f"35 of the 49 already do, and that is what makes the shelf navigable."
    )


def test_the_index_is_not_vacuous():
    """⚠ `TRAPS: could-the-arm-have-fired` applied here: if the row regex matched nothing, every test above
    would pass while checking nothing at all."""
    assert len(IN_TABLE) >= 30, (
        f"only {len(IN_TABLE)} rows parsed from CLAUDE.md — the table format moved"
    )
    assert len(ON_DISK) >= 30, f"only {len(ON_DISK)} instruments found — DESIGN_DIR moved"
