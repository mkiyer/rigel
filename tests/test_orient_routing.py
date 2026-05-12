"""Tests for calibration orientation helper conventions."""

import logging
from types import SimpleNamespace

import pytest

from rigel.calibration._orient import (
    ORIENT_OPP,
    ORIENT_SAME,
    ORIENT_UNINF,
    StrandSummary,
    classify_orient,
)
from rigel.calibration.regions import RegionStrand
from rigel.pipeline import _calibration_strand_summary, _warn_if_calibration_strand_unidentifiable
from rigel.strand_model import StrandModel
from rigel.types import Strand


def test_orient_matrix():
    assert classify_orient(int(RegionStrand.POS), int(Strand.POS)) == ORIENT_SAME
    assert classify_orient(int(RegionStrand.NEG), int(Strand.NEG)) == ORIENT_SAME
    assert classify_orient(int(RegionStrand.POS), int(Strand.NEG)) == ORIENT_OPP
    assert classify_orient(int(RegionStrand.NEG), int(Strand.POS)) == ORIENT_OPP
    assert classify_orient(int(RegionStrand.NONE), int(Strand.POS)) == ORIENT_UNINF
    assert classify_orient(int(RegionStrand.AMBIG), int(Strand.POS)) == ORIENT_UNINF
    assert classify_orient(int(RegionStrand.POS), int(Strand.NONE)) == ORIENT_UNINF
    assert classify_orient(int(RegionStrand.POS), int(Strand.AMBIGUOUS)) == ORIENT_UNINF


def test_strand_summary_uses_readable_signed_contrast_name():
    summary = StrandSummary(p_r1_sense=0.9, n_observations=50)

    assert summary.signed_strand_contrast == pytest.approx(0.8)
    assert summary.strand_specificity == pytest.approx(0.9)
    assert summary.read1_sense
    assert not hasattr(summary, "c")


def test_strand_summary_from_model():
    model = StrandModel()
    model.pos_pos = 9
    model.pos_neg = 1

    summary = StrandSummary.from_model(model)

    assert summary.p_r1_sense == pytest.approx(0.9)
    assert summary.n_observations == 10


def test_strand_summary_validates_probability():
    with pytest.raises(ValueError, match="p_r1_sense"):
        StrandSummary(p_r1_sense=1.5)


def _strand_model(*, same: int, opposite: int) -> StrandModel:
    model = StrandModel()
    model.pos_pos = same
    model.pos_neg = opposite
    return model


def test_calibration_summary_uses_spliced_model_even_when_exonic_diagnostic_identifiable():
    primary = _strand_model(same=10, opposite=18)
    diagnostic = _strand_model(same=652, opposite=1168)
    models = SimpleNamespace(exonic_spliced=primary, exonic=diagnostic)

    summary = _calibration_strand_summary(models)

    assert summary.n_observations == primary.n_observations
    assert summary.p_r1_sense == pytest.approx(primary.p_r1_sense)


def test_calibration_summary_keeps_unstranded_primary_when_diagnostic_is_noise():
    primary = _strand_model(same=162_000, opposite=162_162)
    diagnostic = _strand_model(same=500_000, opposite=500_400)
    models = SimpleNamespace(exonic_spliced=primary, exonic=diagnostic)

    summary = _calibration_strand_summary(models)

    assert summary.n_observations == primary.n_observations
    assert summary.p_r1_sense == pytest.approx(primary.p_r1_sense)


def test_unidentifiable_spliced_model_warns_but_does_not_use_exonic_diagnostic(caplog):
    primary = _strand_model(same=10, opposite=18)
    diagnostic = _strand_model(same=652, opposite=1168)
    models = SimpleNamespace(exonic_spliced=primary, exonic=diagnostic)

    with caplog.at_level(logging.WARNING, logger="rigel.pipeline"):
        _warn_if_calibration_strand_unidentifiable(models)

    assert "Spliced strand model is not identifiable" in caplog.text
    assert "Exonic diagnostic strand signal is identifiable" in caplog.text
    assert "not used for calibration" in caplog.text
