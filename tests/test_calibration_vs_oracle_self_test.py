"""``calibration_vs_oracle.py --self-test`` is a REGRESSION GATE, not a thing to remember to run.

⭐ The instrument carries its own falsification — every comparator perturbed by one ULP and required to
FIRE, with no I/O — which is the project's standing rule that writing a gate is only half of it. But a
``--self-test`` flag nothing invokes is a gate that goes stale silently: the perturbations pass on the
day they are written and nobody runs them again. This runs them in-process on every suite run.

⛔ **In-process, never a subprocess.** A subprocess would need the panel's conda environment resolved
from inside pytest and would hide an import error as a non-zero exit code; importing the module is also
what proves the ``_sibling`` loader still finds ``pass0_vs_oracle`` and ``prior_vs_oracle``.

⚠ It asserts the COUNT as well as the pass, so deleting a perturbation is a failure rather than a
quieter green — the same reason the suite's own baseline is re-derived and never adjusted.
"""

from __future__ import annotations

import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "design" / "calibration_vs_oracle.py"

#: The number of perturbations the instrument ships. ⛔ Re-derive from its output rather than adjusting
#: this: a self-test that quietly lost a case reads exactly like one that never had it.
EXPECTED_GATES = 36


def _load():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("calibration_vs_oracle", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["calibration_vs_oracle"] = module
    spec.loader.exec_module(module)
    return module


def test_every_self_test_perturbation_still_fires():
    module = _load()
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = module.self_test()
    out = buf.getvalue()
    assert rc == 0, f"calibration_vs_oracle --self-test FAILED:\n{out}"
    assert f"{EXPECTED_GATES}/{EXPECTED_GATES} self-test gates fired" in out, (
        f"the self-test no longer reports {EXPECTED_GATES} gates. Re-derive the count from the output "
        f"below and change EXPECTED_GATES only with a perturbation added or removed on purpose:\n{out}"
    )
    assert "⛔" not in out, f"a self-test gate did not fire:\n{out}"
