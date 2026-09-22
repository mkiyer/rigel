"""The instruments that carry their own falsification, run inside the suite.

A `scripts/design/` instrument that perturbs every one of its own comparators and requires each to fire,
with no I/O. A `--self-test` flag that nothing invokes is a gate that goes stale silently: the
perturbations pass on the day they are written and nobody runs them again. This file runs them
in-process on every suite run, and asserts the COUNT as well as the pass, so losing a perturbation is a
failure rather than a quieter green.

In-process, never a subprocess: a subprocess would need the panel's conda environment resolved from
inside pytest and would hide an import error as a non-zero exit code, and importing the module is also
what proves the shared sibling loader still finds the instruments' helpers.
"""

from __future__ import annotations

import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "scripts" / "design"

#: ``(stem, gates, the line the instrument prints when they all fire)``. Re-derive a count from the
#: instrument's own output rather than adjusting it here: a self-test that quietly lost a case reads
#: exactly like one that never had it.
INSTRUMENTS = [
    ("calibration_vs_oracle", 43, "{n}/{n} self-test gates fired"),
    ("ruler_vs_truth", 30, "{n}/{n} self-test gates fired"),
]


def _load(stem: str):
    path = DESIGN / f"{stem}.py"
    sys.path.insert(0, str(DESIGN))
    spec = importlib.util.spec_from_file_location(stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[stem] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("stem,gates,fired", INSTRUMENTS, ids=[i[0] for i in INSTRUMENTS])
def test_every_self_test_perturbation_still_fires(stem, gates, fired):
    module = _load(stem)
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = module.self_test()
    out = buf.getvalue()
    assert rc == 0, f"{stem} --self-test FAILED:\n{out}"
    assert fired.format(n=gates) in out, (
        f"{stem}'s self-test no longer reports {gates} gates. Re-derive the count from the output "
        f"below and change it here only with a perturbation added or removed on purpose:\n{out}"
    )
    assert "⛔" not in out, f"a {stem} self-test gate did not fire:\n{out}"
