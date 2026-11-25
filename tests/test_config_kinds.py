"""A setting of the wrong kind is refused, from a --config YAML as from the API: a quoted "false" for a switch (read
as a true one), a fractional count (truncated downstream), an infinite window, a string size. Each got through."""

from __future__ import annotations

import pytest

from rigel.cli import build_parser, build_pipeline_config, resolve_quant_args
from rigel.config import BamScanConfig, CalibrationConfig


def _from_yaml(tmp_path, body: str):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(body)
    args = build_parser().parse_args(["quant", "--config", str(cfg)])
    resolve_quant_args(args)
    return build_pipeline_config(args)


@pytest.mark.parametrize(
    "body, error",
    [
        ('keep_duplicates: "false"\n', TypeError),
        ("calib_refit_iters: 1.5\n", TypeError),
        ("intron_overhang_tolerance: 2.5\n", TypeError),
        ('scan_buffer_size: "2"\n', TypeError),
    ],
    ids=["a-quoted-false", "a-fractional-refit-count", "a-fractional-tolerance", "a-string-size"],
)
def test_a_yaml_setting_of_the_wrong_kind_is_refused(tmp_path, body, error):
    with pytest.raises(error):
        _from_yaml(tmp_path, body)


def test_a_yaml_boolean_switch_is_read_as_given(tmp_path):
    assert _from_yaml(tmp_path, "keep_duplicates: true\n").scan.skip_duplicates is False


@pytest.mark.parametrize(
    "make",
    [
        lambda: BamScanConfig(skip_duplicates="false"),
        lambda: BamScanConfig(max_frag_length=True),
        lambda: CalibrationConfig(refit_iters=1.5),
    ],
    ids=["a-string-switch", "a-bool-count", "a-fractional-count"],
)
def test_a_config_of_the_wrong_kind_is_refused(make):
    with pytest.raises(TypeError):
        make()
