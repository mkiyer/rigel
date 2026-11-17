"""`quant_accuracy.py` stamps every row with the config it was scored under, and the report reads only the stamps.

The report used to print the shipped ``CalibrationConfig().message_policy`` at render time and the rows carried no
``--set``, so an arm run under ``--set calibration.message_policy=silent`` read "transfer", and nothing could refuse a
pair scored under two configs. Each row now carries its ``--set`` overrides (the effective config's differences from
the shipped one), its message policy whatever set it, and the code revision it ran on; a pair whose stamps differ in
anything but the arm is refused. The revision is provenance, never a refusal: an A/B of a ``src/`` change is two
revisions by design.
"""

from __future__ import annotations

import dataclasses
import io
import json
from contextlib import redirect_stdout

import pytest

from rigel.config import CalibrationConfig, PipelineConfig
import _shared
import quant_accuracy as qa

COND = "gdna_g50_ss_0.99_nrna_mid_capture_off"
PINNED = ("scan.total_threads=1", "em.assignment_mode=fractional")


def _stamp(*specs, arm="base"):
    return qa.row_stamp(arm, COND, _shared.configured(specs), qa.ORACLE_BAM)


def test_the_shipped_config_stamps_no_override_and_the_shipped_policy():
    s = _stamp()
    assert s["settings"] == {}
    assert s["message_policy"] == CalibrationConfig().message_policy
    assert s["revision"]


def test_every_set_is_stamped_by_its_effective_value():
    s = _stamp(*PINNED)
    assert s["settings"] == {"scan.total_threads": 1, "em.assignment_mode": "fractional"}
    assert s["assignment_mode"] == "fractional"


def test_a_silent_run_is_stamped_silent():
    s = _stamp("calibration.message_policy=silent")
    assert s["message_policy"] == "silent"
    assert s["settings"] == {"calibration.message_policy": "silent"}


def test_the_stamp_is_json_and_round_trips():
    s = _stamp(*PINNED, "calibration.message_policy=silent")
    assert json.loads(json.dumps(s)) == s


def _arm_file(tmp_path, name, stamp, axis="transcript"):
    row = {**stamp, "condition": COND, "axis": axis, "count_abs_err": 1.0}
    p = tmp_path / f"{name}.jsonl"
    p.write_text(json.dumps(row) + "\n")
    return p


def test_a_silent_row_beside_a_transfer_row_is_refused(tmp_path):
    a = _arm_file(tmp_path, "transfer", _stamp(*PINNED))
    b = _arm_file(tmp_path, "silent", _stamp(*PINNED, "calibration.message_policy=silent"))
    with pytest.raises(SystemExit, match="stamp"):
        qa.report([a, b])


def test_any_other_set_that_differs_is_refused(tmp_path):
    a = _arm_file(tmp_path, "a", _stamp(*PINNED))
    b = _arm_file(tmp_path, "b", _stamp(*PINNED, "calibration.sweep_logodds_step=0.1"))
    with pytest.raises(SystemExit, match="stamp"):
        qa.report([a, b])


def test_an_unstamped_row_beside_a_stamped_one_is_refused(tmp_path):
    a = _arm_file(tmp_path, "a", _stamp(*PINNED))
    stripped = {k: v for k, v in _stamp(*PINNED).items() if k not in ("settings", "message_policy")}
    b = _arm_file(tmp_path, "b", stripped)
    with pytest.raises(SystemExit, match="stamp"):
        qa.report([a, b])


def test_two_arms_that_differ_in_the_arm_and_the_revision_only_are_compared(tmp_path):
    a = _arm_file(tmp_path, "base", _stamp(*PINNED))
    other = dict(_stamp(*PINNED, arm="no_capture"), revision="0123456789ab")
    b = _arm_file(tmp_path, "no_capture", other)
    qa._check_comparable([(qa._load(a), "base"), (qa._load(b), "no_capture")])


def _text_report(tmp_path, stamp):
    p = _arm_file(tmp_path, "arm", stamp)
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            qa.report([p])
        except (
            KeyError,
            ValueError,
        ):  # a one-row file renders its header before any table needs more
            pass
    return buf.getvalue()


def test_the_text_report_prints_the_stamped_policy_not_the_shipped_one(tmp_path):
    out = _text_report(tmp_path, _stamp("calibration.message_policy=silent"))
    assert "message_policy='silent'" in out, out
    assert "message_policy='transfer'" not in out, out


def test_the_markdown_header_prints_the_stamped_policy_not_the_shipped_one(tmp_path):
    stamp = _stamp("calibration.message_policy=silent")
    lines = qa._md_header(
        {(COND, "library"): stamp}, "arm", [COND], _arm_file(tmp_path, "arm", stamp)
    )
    head = "\n".join(lines)
    assert "message_policy='silent'" in head, head
    assert "calibration.message_policy=silent" in head, head


def test_the_overrides_see_a_field_outside_every_section():
    """Not only what ``--set`` can reach: any field of the run's config that differs from the shipped one."""
    assert _shared.config_overrides(PipelineConfig()) == {}
    assert _shared.config_overrides(dataclasses.replace(PipelineConfig(), second_pass_seed=7)) == {
        "second_pass_seed": 7
    }
