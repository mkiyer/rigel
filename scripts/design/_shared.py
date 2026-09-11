"""The one loader the instruments use to import a sibling instrument by path.

An instrument runs as ``python scripts/design/<name>.py``, so ``scripts/design/`` is not a package and a
sibling cannot be imported by name. ``sibling("pass0_vs_oracle.py")`` loads the file beside this one,
registers it in ``sys.modules`` under its stem BEFORE executing it (a dataclass defined in the sibling
resolves its own module through ``sys.modules`` at class-creation time), and returns the cached module on
every later call, so two instruments loading the same sibling share one copy of it. This is a helper, not
an instrument: it answers no question and has no row in the instrument table.

Usage::

    from _shared import sibling
    P0 = sibling("pass0_vs_oracle.py")
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

DESIGN = Path(__file__).resolve().parent


def sibling(name: str):
    """Load ``scripts/design/<name>`` once and return the module; ``name`` includes ``.py``."""
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, DESIGN / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]
