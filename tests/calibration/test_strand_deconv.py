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
    *, strand_identifiable, strand_observable, sense, antisense, count_gdna_frac,
    llr=0.0, overdispersion=0.0, rna_overdispersion=0.0, rna_sense_frac=0.99, mass=100.0,
):
    r = _deconv_per_node(
        np.array([mass]),
        np.array([0.0]),
        np.array([float(sense)]),
        np.array([float(antisense)]),
        np.array([bool(strand_observable)]),
        np.array([float(count_gdna_frac)]),
        strand_identifiable=strand_identifiable,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=overdispersion,
        rna_strand_overdispersion=rna_overdispersion,
        gdna_strand_llr_bias=llr,
        n_grid=400,
    )
    return float(r.gdna_frac[0])


# --------------------------------------------------------------------------- #
# routing — strand vs count handoff (disjoint, no product)
# --------------------------------------------------------------------------- #


def test_unstranded_routes_to_count():
    # strand_identifiable=False ⇒ the node uses count_gdna_frac regardless of its sense split.
    frac = _frac(strand_identifiable=False, strand_observable=True, sense=50, antisense=50,
                 count_gdna_frac=0.3)
    assert frac == pytest.approx(0.3)


def test_strand_unobservable_routes_to_count():
    # strand-unobservable node (e.g. AMBIG) ⇒ count_gdna_frac, even in a strand-identifiable library.
    frac = _frac(strand_identifiable=True, strand_observable=False, sense=50, antisense=50,
                 count_gdna_frac=0.7)
    assert frac == pytest.approx(0.7)


def test_strand_observable_uses_strand_not_count():
    # Identifiable + observable ⇒ the strand posterior governs, NOT count_gdna_frac. A symmetric
    # node (gDNA's signature) reads ~all gDNA even though count_gdna_frac says 0.
    frac = _frac(strand_identifiable=True, strand_observable=True, sense=50, antisense=50,
                 count_gdna_frac=0.0, rna_sense_frac=0.99)
    assert frac > 0.9  # strand says gDNA; the count fraction (0) is ignored


def test_strand_observable_rna_node_reads_rna():
    # A node at the RNA sense rate reads ~all RNA from strand alone (count fraction ignored).
    frac = _frac(strand_identifiable=True, strand_observable=True, sense=99, antisense=1,
                 count_gdna_frac=1.0, rna_sense_frac=0.99)
    assert frac < 0.1


# --------------------------------------------------------------------------- #
# strand module — LLR bias (FP aversion)
# --------------------------------------------------------------------------- #


def _gdna_frac(llr, sense, antisense, *, overdispersion=0.0):
    return _frac(strand_identifiable=True, strand_observable=True, sense=sense, antisense=antisense,
                 count_gdna_frac=0.0, llr=llr, overdispersion=overdispersion)


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
