"""Helpers shared by the instruments under ``scripts/design/``; a helper module, not an instrument.

An instrument runs as ``python scripts/design/<name>.py``, so ``scripts/design/`` is not a package and a
sibling cannot be imported by name: ``sibling("_oracle_arms.py")`` loads the file beside this one once and
shares the module. The module also holds the data root and the panels' paths, ``set_field`` (the one
``SECTION.FIELD=VALUE`` parser behind every ``--set``), the override fields, the stratum readers, the pool
ledger and the subprocess sharder behind every ``--jobs``.

Usage::

    from _shared import sibling
    OA = sibling("_oracle_arms.py")
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import os
import re
import subprocess
import sys
import typing
from pathlib import Path

DESIGN = Path(__file__).resolve().parent


#: The data root every script reads by default — ``$RIGEL_RUNS``, else ``Downloads/rigel_runs`` in the home
#: directory — and the panel and its index under it: one spelling, here, for every script. The sim configs
#: name their paths under ``$RIGEL_RUNS``, so a script that reads a config exports this first.
RUNS = Path(os.environ.get("RIGEL_RUNS") or Path.home() / "Downloads" / "rigel_runs").expanduser()
DEFAULT_SUITE = RUNS / "suite" / "ladder"
DEFAULT_INDEX = RUNS / "suite" / "rigel_index"

#: the substrates, each: (index dir, the panel dir holding `scan_cache/` and `oracle_cache/`)
PANELS = {
    "test": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios"),
    # the same chromosome under the two other probe designs (junction-only and sparse probes)
    "test_junction": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_junction"),
    "test_sparse": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_sparse"),
    "ladder": (DEFAULT_INDEX, DEFAULT_SUITE),
}

#: The mass arrays ``OracleTruth.override_masses`` replaces. Named here so the ``noop`` gate can re-inject
#: exactly this set from the shipped result and demand byte-identity: an override of a field nothing reads
#: would never run.
OVERRIDE_FIELDS = (
    "count_gdna_region",
    "count_rna_region",
    "count_gdna_boundary",
    "count_rna_boundary",
    "count_rna_spliced_boundary",
    "count_rna_sj",
)


#: The strand axis's two halves, keyed by the strand specificity a condition name spells (``_ss_0.50_``). A
#: library between them is in neither half: the test panel's ss 0.70 rows are a stratum of their own, never
#: pooled into a half whose bar they would move.
STRAND_HALVES = {"0.50": "unstranded", "0.99": "stranded"}

#: The four strata of the two halves, in report order. ``strata`` appends any other a panel carries.
STRATA = (
    ("stranded", "capture OFF"),
    ("stranded", "capture ON"),
    ("unstranded", "capture OFF"),
    ("unstranded", "capture ON"),
)


def strandedness(cond: str) -> str:
    """``unstranded`` (ss 0.50), ``stranded`` (ss 0.99), or ``ss <value>`` for any other specificity,
    which is reported apart. A name that spells no specificity raises rather than landing in a half."""
    token = re.search(r"_ss_(\d+\.\d+)(?:_|$)", cond)
    if token is None:
        raise ValueError(f"{cond!r} names no strand specificity (`_ss_<value>_`)")
    return STRAND_HALVES.get(token.group(1), f"ss {token.group(1)}")


def stratum(cond: str) -> tuple[str, str]:
    """The panel's two axes: the strand half (or the specificity it is apart at) and capture."""
    return (strandedness(cond), "capture ON" if "capture_on" in cond else "capture OFF")


def strata(conds) -> list[tuple[str, str]]:
    """:data:`STRATA`, then every other stratum among ``conds``, sorted: the list a per-stratum report
    iterates, so a condition apart from both halves is printed on its own row rather than dropped."""
    return [*STRATA, *sorted({stratum(c) for c in conds} - set(STRATA))]


def is_zero_gdna(cond: str) -> bool:
    """``g00``, the zero-gDNA control. Truth is exactly 0, so every gDNA fragment in the prior there is a
    false positive with nothing to cancel it, and a relative change is unbounded. Read per condition."""
    return "_g00_" in cond


def pool_ledger(condition_dir: Path) -> dict:
    """The simulator's starting fragment count per origin pool (``gdna``, ``mrna``, ``nrna``).

    ``calibrate`` answers in gDNA and RNA only, so its RNA answer is scored against ``mrna + nrna``; the
    split between annotated and nascent RNA is the EM's, scored in ``quant_accuracy.py``'s pool table. Read
    from ``truth_summary.json``, never from the condition name; a missing pool reads 0, a missing file raises.
    """
    summary = json.loads((Path(condition_dir) / "truth_summary.json").read_text())
    counts = summary["origin_counts"]
    return {k: float(counts.get(k, 0.0)) for k in ("gdna", "mrna", "nrna")}


def run_shards(script: Path, names: list, jobs: int, scratch: Path, shard_args) -> list[Path]:
    """Run ``script`` on ``names`` split round-robin into ``jobs`` shards, each its own subprocess, and return
    the shards' output files in shard order; ``shard_args(out, conditions)`` is one shard's command line.

    Shards, not threads: conditions share nothing but a read-only index and cache, so this changes no number.
    Every shard is waited on before any failure is raised, and one failed shard refuses the whole panel,
    because a short output would read as a complete one.
    """
    shards = [s for s in (names[i::jobs] for i in range(jobs)) if s]
    scratch.mkdir(parents=True, exist_ok=True)
    outs = [scratch / f"shard{i}.json" for i in range(len(shards))]
    procs = [
        subprocess.Popen([sys.executable, str(script), *shard_args(out, shard)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for out, shard in zip(outs, shards)
    ]
    failed = False
    for i, proc in enumerate(procs):
        log, _ = proc.communicate()
        if proc.returncode != 0:
            failed = True
            print(f"  ⛔ shard {i} FAILED (rc={proc.returncode}):\n{log}", flush=True)
        else:
            print(f"  shard {i}: {len(shards[i])} conditions ok", flush=True)
    if failed:
        raise SystemExit("⛔ a shard failed; refusing to report a partial panel")
    return outs


def set_args(settings) -> list[str]:
    """``--set`` specs as the command-line arguments a shard is handed."""
    return [arg for spec in settings for arg in ("--set", spec)]


def sibling(name: str):
    """Load ``scripts/design/<name>`` once and return the module; ``name`` includes ``.py``."""
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, DESIGN / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


def set_field(cfg, spec: str):
    """``SECTION.FIELD=VALUE`` applied to a frozen ``PipelineConfig``: the one parser behind every
    instrument's ``--set``, so an arm is spelled the same way everywhere. The value is typed from the field's
    annotation (a ``Literal`` by membership; ``bool``, ``int``, ``float``, ``str``; ``none`` where the field
    admits ``None``). An unknown section or field, a malformed spec, or a value the type refuses raises
    ``SystemExit`` naming it; the config's own ``__post_init__`` validation still runs on the result."""
    dotted, sep, raw = spec.partition("=")
    section_name, dot, field_name = dotted.partition(".")
    if not (sep and dot and section_name and field_name):
        raise SystemExit(f"⛔ --set expects SECTION.FIELD=VALUE, got {spec!r}")
    section = getattr(cfg, section_name, None)
    if not dataclasses.is_dataclass(section):
        raise SystemExit(f"⛔ --set: {type(cfg).__name__} has no section {section_name!r}")
    if field_name not in {f.name for f in dataclasses.fields(section)}:
        raise SystemExit(f"⛔ --set: {section_name} has no field {field_name!r}")
    hint = typing.get_type_hints(type(section))[field_name]
    # a ``Literal[...]`` field's args are its values, not callables, so it is matched by membership
    literal = typing.get_origin(hint) is typing.Literal
    kinds = [t for t in (typing.get_args(hint) or (hint,)) if t is not type(None)]
    try:
        if literal:
            allowed = {str(v): v for v in typing.get_args(hint)}
            value = allowed[raw]
        elif len(kinds) < len(typing.get_args(hint) or (hint,)) and raw.lower() == "none":
            value = None
        elif bool in kinds:
            words = {"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False}
            value = words[raw.lower()]
        elif int in kinds:
            value = int(raw)
        elif float in kinds:
            value = float(raw)
        else:
            value = kinds[0](raw)
    except (KeyError, ValueError, TypeError) as e:
        raise SystemExit(f"⛔ --set: {dotted} takes {hint!r}, got {raw!r}") from e
    return dataclasses.replace(
        cfg, **{section_name: dataclasses.replace(section, **{field_name: value})}
    )
