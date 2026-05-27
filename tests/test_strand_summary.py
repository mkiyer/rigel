"""Tests for calibration StrandSummary posterior fields."""

from __future__ import annotations

import pytest

from rigel.calibration.strand_summary import StrandSummary
from rigel.strand_model import StrandModel
from rigel.types import Strand


def test_strand_summary_from_model_copies_minor_rate_posterior() -> None:
    model = StrandModel()
    for _ in range(95):
        model.observe(Strand.POS, Strand.POS)
    for _ in range(5):
        model.observe(Strand.POS, Strand.NEG)

    summary = StrandSummary.from_model(model)

    assert summary.n_same == 95
    assert summary.n_opposite == 5
    assert summary.minor_rate_alpha == pytest.approx(6.0)
    assert summary.minor_rate_beta == pytest.approx(96.0)
    assert summary.minor_rate_mean == pytest.approx(6.0 / 102.0)


def test_strand_summary_uninformative_has_beta_one_one_minor_rate() -> None:
    summary = StrandSummary.uninformative()

    assert summary.p_r1_sense == pytest.approx(0.5)
    assert summary.n_observations == 0
    assert summary.n_same == 0
    assert summary.n_opposite == 0
    assert summary.minor_rate_alpha == pytest.approx(1.0)
    assert summary.minor_rate_beta == pytest.approx(1.0)


def test_strand_summary_requires_explicit_counts() -> None:
    with pytest.raises(ValueError, match=r"n_same \+ n_opposite"):
        StrandSummary(p_r1_sense=0.95, n_observations=100)


def test_strand_summary_derives_minor_rate_from_explicit_counts() -> None:
    summary = StrandSummary(p_r1_sense=0.95, n_observations=100, n_same=95, n_opposite=5)

    assert summary.n_same == 95
    assert summary.n_opposite == 5
    assert summary.minor_rate_alpha == pytest.approx(6.0)
    assert summary.minor_rate_beta == pytest.approx(96.0)


def test_strand_summary_rejects_probability_count_mismatch() -> None:
    with pytest.raises(ValueError, match="p_r1_sense"):
        StrandSummary(p_r1_sense=0.96, n_observations=100, n_same=95, n_opposite=5)


def test_strand_summary_rejects_invalid_minor_rate_parameters() -> None:
    with pytest.raises(ValueError, match="minor_rate_alpha"):
        StrandSummary(minor_rate_alpha=0.0)
    with pytest.raises(ValueError, match="minor_rate_beta"):
        StrandSummary(minor_rate_beta=-1.0)


def test_strand_summary_rejects_mismatched_minor_rate_parameters() -> None:
    with pytest.raises(ValueError, match="minor_rate_alpha"):
        StrandSummary(
            p_r1_sense=0.9,
            n_observations=10,
            n_same=9,
            n_opposite=1,
            minor_rate_alpha=5.0,
            minor_rate_beta=10.0,
        )
