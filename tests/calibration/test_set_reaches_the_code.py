"""A ``--set`` the instrument never hands to the code it measures is refused, before any I/O.

`calibration_vs_oracle.py`, `policy_benchmark.py` and `ruler_vs_truth.py` calibrate a cached scan and pass
``config.calibration`` alone, so an ``em.``, ``scan.`` or ``scoring.`` override was applied to a config nothing read
and the run reported the shipped numbers under the override's name ("on BOTH arms"). `ruler_vs_truth.py --scale`
runs the pipeline to the prior assembly, so it consumes the scan and the scoring too, and never the EM.
`policy_benchmark.py` sets the message policy itself, once per arm, so a ``--set`` of it would make every row one
policy.

Every run here points at an empty directory, so the unfixed instruments fail later and differently (no panel, no
index): only the refusal's own message passes.
"""

from __future__ import annotations

import sys

import pytest

import _shared
import calibration_vs_oracle as cvo
import policy_benchmark as pb
import ruler_vs_truth as rvt

#: the refusal names the sections it hands on, so a different early exit cannot pass for it
REFUSAL = "never read"


def _run(module, argv, monkeypatch):
    monkeypatch.setattr(sys, "argv", [module.__name__, *argv])
    with pytest.raises(SystemExit) as exc:
        module.main()
    return str(exc.value)


@pytest.fixture
def empty_panels(tmp_path, monkeypatch):
    panels = {name: (tmp_path / "idx", tmp_path / "panel") for name in _shared.PANELS}
    monkeypatch.setattr(pb, "PANELS", panels)
    monkeypatch.setattr(rvt, "PANELS", panels)
    return tmp_path


@pytest.mark.parametrize(
    "spec", ["em.mode=map", "scan.total_threads=1", "scoring.overhang_log_penalty=-1.0"]
)
def test_calibration_vs_oracle_refuses_a_section_it_does_not_consume(spec, tmp_path, monkeypatch):
    msg = _run(
        cvo, ["--suite", str(tmp_path), "--index", str(tmp_path), "--set", spec], monkeypatch
    )
    assert REFUSAL in msg and spec in msg


@pytest.mark.parametrize("spec", ["em.mode=map", "scan.total_threads=1"])
def test_policy_benchmark_refuses_a_section_it_does_not_consume(spec, empty_panels, monkeypatch):
    msg = _run(pb, ["--panel", "test", "--set", spec], monkeypatch)
    assert REFUSAL in msg and spec in msg


def test_policy_benchmark_refuses_a_set_of_the_policy_it_varies(empty_panels, monkeypatch):
    msg = _run(pb, ["--panel", "test", "--set", "calibration.message_policy=silent"], monkeypatch)
    assert "once per arm" in msg


@pytest.mark.parametrize("spec", ["em.mode=map", "scan.total_threads=1"])
def test_ruler_vs_truth_refuses_a_section_its_arms_do_not_consume(spec, empty_panels, monkeypatch):
    msg = _run(rvt, ["--panel", "test", "--set", spec], monkeypatch)
    assert REFUSAL in msg and spec in msg


def test_ruler_vs_truth_scale_consumes_the_scan_and_the_scoring_but_never_the_em(
    empty_panels, monkeypatch
):
    msg = _run(rvt, ["--panel", "test", "--scale", "--set", "em.mode=map"], monkeypatch)
    assert REFUSAL in msg
    assert rvt.consumed_sections(True) == {"scan", "scoring", "calibration"}


def test_a_calibration_set_is_accepted_everywhere():
    """The refusal is not a blanket one: every instrument still takes the section it hands on."""
    for consumed in (cvo.CONSUMED_SECTIONS, pb.CONSUMED_SECTIONS, rvt.consumed_sections(False)):
        _shared.refuse_unconsumed(["calibration.sweep_logodds_step=0.1"], consumed, "x")
