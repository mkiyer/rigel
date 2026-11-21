"""The instruments that carry their own falsification, run inside the suite.

A `scripts/design/` instrument's `--self-test` perturbs every one of its own comparators, with no I/O, and requires
each to fire. A self-test nothing invokes goes stale silently, so this file runs each in-process on every suite run
and asserts the count of gates as well as the pass: losing a perturbation fails. In-process, not a subprocess, so an
import error surfaces as itself and the instruments' imports of their sibling helpers are proved to resolve.
"""

from __future__ import annotations

import importlib
import io
from contextlib import redirect_stdout

import pytest

#: ``(stem, gates, the line the instrument prints when they all fire)``. Re-derive a count from the
#: instrument's own output rather than adjusting it here: a self-test that quietly lost a case reads
#: exactly like one that never had it.
INSTRUMENTS = [
    ("calibration_vs_oracle", 26, "{n}/{n} self-test gates fired"),
    ("ruler_vs_truth", 50, "{n}/{n} self-test gates fired"),
    ("paired_score", 47, "{n}/{n} self-test gates fired"),
]


@pytest.mark.parametrize("stem,gates,fired", INSTRUMENTS, ids=[i[0] for i in INSTRUMENTS])
def test_every_self_test_perturbation_fires(stem, gates, fired):
    module = importlib.import_module(stem)
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
