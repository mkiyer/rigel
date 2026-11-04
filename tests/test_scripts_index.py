"""The hand-written indexes of `scripts/` are checked against the disk in both directions.

`scripts/design/` is indexed by a table in `CLAUDE.md` and `scripts/profiling/` by `scripts/README.md`. Every script
on disk must have a row and every row must name a file that exists; every script in the three trees must also import
and carry a module docstring. The gates do not judge whether a row's prose is accurate, which no test can.
"""

from __future__ import annotations

import ast
import importlib.util
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CLAUDE = ROOT / "CLAUDE.md"
SCRIPTS = ROOT / "scripts"
DESIGN_DIR = SCRIPTS / "design"

#: rows look like ``| `design/policy_benchmark.py` | … |``
_ROW = re.compile(r"`design/([a-z0-9_]+\.py)`")

#: the instruments — every design/ file except the package marker and the `_`-prefixed helper modules,
#: which answer no question and have no row (they are still import- and docstring-gated below).
ON_DISK = frozenset(
    p.name
    for p in DESIGN_DIR.glob("*.py")
    if p.name != "__init__.py" and not p.name.startswith("_")
)
IN_TABLE = frozenset(_ROW.findall(CLAUDE.read_text()))

SIM_DIR = SCRIPTS / "sim"

#: The third tree, held to the same import and docstring gates as `design/` and `sim/`.
PROFILING_DIR = SCRIPTS / "profiling"


def _instruments(directory: pathlib.Path) -> list[pathlib.Path]:
    return [p for p in directory.glob("*.py") if p.name != "__init__.py"]


#: every script the import + docstring gates cover, all three trees.
ALL_SCRIPTS = sorted(
    _instruments(DESIGN_DIR) + _instruments(SIM_DIR) + _instruments(PROFILING_DIR),
    key=lambda p: (p.parent.name, p.name),
)


def _rel(path: pathlib.Path) -> str:
    return f"{path.parent.name}/{path.name}"


def _import_error(path: pathlib.Path) -> str | None:
    """Import *path* the way ``python scripts/<tree>/x.py`` would, and return the error, if any.

    Two things a naive import gets wrong, each a false failure on an instrument that runs: the module must
    be in `sys.modules` before `exec_module`, or a dataclass resolving its own `__module__` gets `None`; and
    the script's own directory is `sys.path[0]` under a real invocation, which is how the sibling-importing
    instruments (`_sibling`) resolve.
    """
    name = f"_gate_{path.stem}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    except SystemExit:
        pass  # a script that argues with argv at import is still connected
    except Exception as exc:  # noqa: BLE001
        return f"{type(exc).__name__}: {exc}"
    finally:
        sys.path.remove(str(path.parent))
        sys.modules.pop(name, None)
    return None


def test_every_instrument_still_imports():
    """A `src/` deletion breaks instruments silently, since the rest of the suite stays green.

    Importing is the cheapest check that an instrument is still connected to the code it measures. Import
    only, never execution: an instrument's numbers need its substrate, so this is a connectivity gate, not
    a claim that the script is correct.
    """
    broken = {_rel(p): err for p in ALL_SCRIPTS if (err := _import_error(p))}
    assert not broken, (
        f"scripts that no longer import: {broken}. An instrument that cannot be imported cannot be run, "
        f"and nothing else in this file would have caught it."
    )


def test_every_instrument_is_in_the_index():
    """Every instrument on disk has a row in `CLAUDE.md`'s table."""
    missing = sorted(ON_DISK - IN_TABLE)
    assert not missing, (
        f"{len(missing)} instruments are not in CLAUDE.md's table: {missing}. Add a row describing WHAT "
        f"QUESTION IT ANSWERS, or delete the file."
    )


def test_the_index_does_not_point_at_ghosts():
    """Every row in `CLAUDE.md`'s table names a file that exists."""
    ghosts = sorted(IN_TABLE - ON_DISK)
    assert not ghosts, (
        f"CLAUDE.md's table has rows for files that do not exist: {ghosts}. A reader who trusts the table "
        f"goes looking for an instrument that was deleted."
    )


def test_the_profiling_tree_is_indexed_in_the_scripts_readme():
    """The profiling drivers are indexed in `scripts/README.md`'s `profiling/` row, in both directions:
    every driver on disk is named there, and every name there exists.
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


def test_every_instrument_says_what_it_is_for():
    """Every script in the three trees carries a module docstring, the only thing a reader has before
    running it.
    """
    bare = sorted(
        _rel(p)
        for p in ALL_SCRIPTS
        if len((ast.get_docstring(ast.parse(p.read_text())) or "").strip()) <= 60
    )
    assert not bare, (
        f"scripts with no usable module docstring: {bare}. Lead with the QUESTION each answers — that "
        f"is what makes the shelf navigable."
    )


def test_the_index_is_not_vacuous():
    """If the row regex matched nothing, every test above would pass while checking nothing."""
    assert IN_TABLE, "no rows parsed from CLAUDE.md — the table format moved"
    assert ON_DISK, "no instruments found — DESIGN_DIR moved"
