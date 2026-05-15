"""Tests for gDNA overlap effective-length geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from rigel.calibration._exposure import gdna_eff_len_for_loci
from rigel.locus import Locus


@dataclass(frozen=True)
class _FakeFL:
    pmf: np.ndarray

    @property
    def mean(self) -> float:
        ell = np.arange(self.pmf.size, dtype=np.float64)
        return float(np.dot(ell, self.pmf))


def _fl(length_probs: dict[int, float]) -> _FakeFL:
    max_len = max(length_probs)
    pmf = np.zeros(max_len + 1, dtype=np.float64)
    for length, prob in length_probs.items():
        pmf[length] = prob
    pmf /= pmf.sum()
    return _FakeFL(pmf=pmf)


def test_single_interval_fast_path_matches_closed_form():
    fl = _fl({3: 0.5, 5: 0.5})
    locus = Locus(ref="chr1", ref_id=0, start=100, end=200)

    eff = gdna_eff_len_for_loci((locus,), [1000], fl)

    assert eff == pytest.approx(locus.span + fl.mean - 1.0)


def test_single_interval_clips_at_left_contig_edge():
    fl = _fl({5: 1.0})
    locus = Locus(ref="chr1", ref_id=0, start=1, end=4)

    eff = gdna_eff_len_for_loci((locus,), [10], fl)

    # Start window is [1 - 5 + 1, 4) = [-3, 4), clipped to [0, 4).
    assert eff == pytest.approx(4.0)


def test_nearby_intervals_merge_expanded_start_windows():
    fl = _fl({10: 1.0})
    loci = (
        Locus(ref="chr1", ref_id=0, start=100, end=110),
        Locus(ref="chr1", ref_id=0, start=115, end=125),
    )

    eff = gdna_eff_len_for_loci(loci, [1000], fl)

    # [91, 110) and [106, 125) merge to [91, 125).
    assert eff == pytest.approx(34.0)


def test_multi_reference_loci_sum_by_reference():
    fl = _fl({5: 1.0})
    loci = (
        Locus(ref="chr1", ref_id=0, start=10, end=20),
        Locus(ref="chr2", ref_id=1, start=0, end=5),
    )

    eff = gdna_eff_len_for_loci(loci, [1000, 10], fl)

    # chr1: [6, 20) -> 14 starts. chr2: [-4, 5) clipped to [0, 5) -> 5.
    assert eff == pytest.approx(19.0)
