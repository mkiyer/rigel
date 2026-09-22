"""The one loader the instruments use to import a sibling instrument by path.

An instrument runs as ``python scripts/design/<name>.py``, so ``scripts/design/`` is not a package and a
sibling cannot be imported by name. ``sibling("pass0_vs_oracle.py")`` loads the file beside this one,
registers it in ``sys.modules`` under its stem BEFORE executing it (a dataclass defined in the sibling
resolves its own module through ``sys.modules`` at class-creation time), and returns the cached module on
every later call, so two instruments loading the same sibling share one copy of it. This is a helper, not
an instrument: it answers no question and has no row in the instrument table. It also holds
``set_field``, the one ``SECTION.FIELD=VALUE`` parser behind every instrument's ``--set``, and the panel's
default paths, the six override fields and the two stratum readers the oracle instruments share.

Usage::

    from _shared import sibling
    P0 = sibling("pass0_vs_oracle.py")
"""

from __future__ import annotations

import dataclasses
import importlib.util
import sys
import typing
from pathlib import Path

DESIGN = Path(__file__).resolve().parent


#: The panel every instrument reads by default, and its index — one spelling, here, for every instrument.
RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = RUNS / "suite" / "ladder"
DEFAULT_INDEX = RUNS / "suite" / "rigel_index"

#: The mass arrays ``OracleTruth.override_masses`` replaces. Named here so the ``noop`` gate can
#: re-inject exactly this set from the SHIPPED result and demand byte-identity — an override applied
#: to a field nothing reads is an override that never ran (TRAPS: an-ablation-that-never-ran).
OVERRIDE_FIELDS = (
    "count_gdna_region",
    "count_rna_region",
    "count_gdna_boundary",
    "count_rna_boundary",
    "count_rna_spliced_boundary",
    "count_rna_sj",
)


def stratum(cond: str) -> tuple[str, str]:
    """The panel's two binary axes."""
    return (
        "stranded" if "ss_0.99" in cond else "unstranded",
        "capture ON" if "capture_on" in cond else "capture OFF",
    )


def is_zero_gdna(cond: str) -> bool:
    """``g00`` — the owner-required ZERO-gDNA control. Truth is exactly 0, so every gDNA fragment in
    the prior there is a false positive with nothing to cancel it, and a relative change is unbounded.
    Reported on its own row, never inside ALL."""
    return "_g00_" in cond


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
    """``SECTION.FIELD=VALUE`` applied to a frozen ``PipelineConfig`` — the one parser behind every
    instrument's ``--set``, so an arm is a config value spelled the same way on every instrument. The
    value is typed from the field's annotation (``bool``, ``int``, ``float``, ``str``; ``none`` where the
    field admits ``None``), so an ``int | None`` field such as ``sweep_block_slots`` takes an int. An unknown
    section or field, a malformed spec, or a value the field's type refuses raises ``SystemExit``
    naming it; the config's own validation (``__post_init__``) still runs on the result."""
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
    # A ``Literal[...]`` field admits exactly its listed values, by membership: its "args" are the values
    # themselves, not callables, so the kind-based coercion below cannot apply to it (it once raised on
    # every ``em.mode=map`` and ``em.warm_start=...``, leaving no instrument able to run a MAP arm).
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
