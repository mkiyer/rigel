"""Per-node strand/count handoff (`strand_deconv._deconv_per_node`).

Two concerns: (1) the **routing** — a node takes the strand posterior iff the library is
strand-identifiable AND the node is strand-observable, else the count fraction; (2) the strand
module's **LLR bias** (``gdna_strand_llr_bias``) — the FP-aversion knob that siphons unspliced mass
into gDNA monotonically, neutral at 0, full at the extremes.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.strand_deconv import _deconv_per_node


def _frac(
    *, strand_observable=True, sense, antisense, count_gdna_frac,
    llr=0.0, overdispersion=0.0, rna_overdispersion=0.0, rna_sense_frac=0.99, mass=100.0,
):
    r = _deconv_per_node(
        np.array([mass]),
        np.array([0.0]),
        np.array([float(sense)]),
        np.array([float(antisense)]),
        np.array([bool(strand_observable)]),
        np.array([float(count_gdna_frac)]),
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=overdispersion,
        rna_strand_overdispersion=rna_overdispersion,
        gdna_strand_llr_bias=llr,
        n_grid=400,
    )
    return float(r.gdna_frac[0])


# --------------------------------------------------------------------------- #
# precision-weighted strand→count deference: g = w·g_strand + (1−w)·g_count, w = (2κ−1)²
# --------------------------------------------------------------------------- #


def test_unstranded_defers_fully_to_count():
    # κ=½ ⇒ w=(2κ−1)²=0 ⇒ the node is count-only regardless of its sense split.
    frac = _frac(rna_sense_frac=0.5, sense=50, antisense=50, count_gdna_frac=0.3)
    assert frac == pytest.approx(0.3)


def test_strand_unobservable_is_count_only():
    # No defined sense (e.g. AMBIG) ⇒ count_gdna_frac, even in a strongly stranded library.
    frac = _frac(strand_observable=False, sense=50, antisense=50, count_gdna_frac=0.7,
                 rna_sense_frac=0.99)
    assert frac == pytest.approx(0.7)


def test_high_ss_is_mostly_strand():
    # κ=0.99 ⇒ w≈0.96. A symmetric node (gDNA's signature) reads ~all gDNA, nearly ignoring the
    # count fraction (0): frac ≈ 0.96·g_strand + 0.04·0.
    frac = _frac(rna_sense_frac=0.99, sense=50, antisense=50, count_gdna_frac=0.0)
    assert frac > 0.85  # strand dominates; the 4% count weight pulls it down only slightly


def test_high_ss_rna_node_reads_rna():
    frac = _frac(rna_sense_frac=0.99, sense=99, antisense=1, count_gdna_frac=1.0)
    assert frac < 0.1  # strand says RNA; the count fraction (1.0) gets only ~4% weight


def test_intermediate_ss_blends():
    # κ=0.75 ⇒ w=0.25: a symmetric (gDNA-signature) node blends mostly toward the count fraction (0).
    frac = _frac(rna_sense_frac=0.75, sense=50, antisense=50, count_gdna_frac=0.0)
    assert 0.1 < frac < 0.4  # ≈ 0.25·g_strand + 0.75·0 — between count (0) and strand (~1)


def test_weight_monotone_in_strand_specificity():
    # The same gDNA-signature node calls more gDNA as strand specificity rises (w grows).
    fracs = [_frac(rna_sense_frac=k, sense=50, antisense=50, count_gdna_frac=0.0)
             for k in (0.50, 0.60, 0.75, 0.90, 0.99)]
    assert all(b >= a - 1e-9 for a, b in zip(fracs, fracs[1:]))  # non-decreasing in κ
    assert fracs[0] == pytest.approx(0.0)  # κ=½ → pure count (0)


# --------------------------------------------------------------------------- #
# FP-aversion LLR bias (on the final blended fraction)
# --------------------------------------------------------------------------- #


def _gdna_frac(llr, sense, antisense, *, overdispersion=0.0):
    # count_gdna_frac=1.0 so the universal LLR bias can siphon the full node into gDNA at λ→+∞.
    return _frac(sense=sense, antisense=antisense, count_gdna_frac=1.0, llr=llr,
                 overdispersion=overdispersion, rna_sense_frac=0.99)


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_llr_bias_siphons_monotonically(overdispersion):
    fracs = [_gdna_frac(llr, 8, 2, overdispersion=overdispersion)
             for llr in (-2.0, 0.0, 1.0, 2.0, 4.0)]
    assert all(b >= a - 1e-9 for a, b in zip(fracs, fracs[1:]))  # non-decreasing in λ
    assert fracs[-1] > fracs[1] + 0.1  # materially more gDNA at high λ than at neutral


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_llr_bias_neutral_then_full_siphon(overdispersion):
    neutral = _gdna_frac(0.0, 95, 5, overdispersion=overdispersion)
    extreme = _gdna_frac(40.0, 95, 5, overdispersion=overdispersion)
    rna_lean = _gdna_frac(-40.0, 95, 5, overdispersion=overdispersion)
    assert neutral < 0.2  # at λ=0 a strand-specific node reads RNA
    assert extreme > 0.99  # λ→+∞ siphons it (and every unspliced node) into gDNA
    assert rna_lean < 0.01  # λ→−∞ is the symmetric RNA limit


def test_llr_bias_zero_is_unbiased_median():
    frac = _gdna_frac(0.0, 8, 2)
    assert 0.0 < frac < 1.0
