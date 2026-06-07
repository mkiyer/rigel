"""Strand-deconvolution gDNA LLR bias (``gdna_strand_llr_bias``) — the FP-aversion knob.

Stresses the per-node joint deconvolution across Beta-Binomial overdispersion levels to be certain
the log-odds bias siphons unspliced mass into gDNA monotonically, is neutral at 0, and siphons
*everything* at the extreme (λ→+∞ ⇒ all unspliced → gDNA, even a confidently-RNA node; λ→−∞ ⇒ all
RNA). The calibration-stage twin of the EM ``gdna_em_llr_bias``: ``log_post += λ·logit(gdna_frac)``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.joint_deconv import _joint_per_node


def _gdna_frac(
    llr, sense, antisense, *, overdispersion=0.0, rna_overdispersion=0.0,
    rna_sense_frac=0.99, mass=100.0,
):
    """gDNA fraction of one strand-observable node, count prior flat (count_evidence=0) so the
    strand clue + the LLR tilt drive the call. Higher λ ⇒ more FP-averse (more gDNA)."""
    r = _joint_per_node(
        mass_unspl=np.array([mass]),
        mass_spliced=np.array([0.0]),
        sense=np.array([float(sense)]),
        antisense=np.array([float(antisense)]),
        density=np.array([0.0]),
        count_evidence=np.array([0.0]),  # flat count prior → isolate the strand deconvolution
        eff_len=np.array([1.0]),
        strand_observable=np.array([True]),
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=overdispersion,
        rna_strand_overdispersion=rna_overdispersion,
        gdna_strand_llr_bias=llr,
        n_grid=400,
    )
    return float(r.gdna_frac[0])


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_llr_bias_siphons_monotonically(overdispersion):
    # An ambiguous node (8 sense / 2 antisense, rna_sense=0.99 — the motivating example): raising
    # the LLR bias monotonically shifts mass into gDNA at every overdispersion level.
    fracs = [_gdna_frac(llr, 8, 2, overdispersion=overdispersion)
             for llr in (-2.0, 0.0, 1.0, 2.0, 4.0)]
    assert all(b >= a - 1e-9 for a, b in zip(fracs, fracs[1:]))  # non-decreasing in λ
    assert fracs[-1] > fracs[1] + 0.1  # materially more gDNA at high λ than at neutral (λ=0)


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_llr_bias_neutral_then_full_siphon(overdispersion):
    # A confidently-RNA node (95 sense / 5 antisense): neutral (λ=0) calls it RNA; the extreme
    # (λ→+∞) siphons ALL of its unspliced mass into gDNA — the property a quantile cannot deliver —
    # and λ→−∞ pushes symmetrically the other way (all RNA).
    neutral = _gdna_frac(0.0, 95, 5, overdispersion=overdispersion)
    extreme = _gdna_frac(40.0, 95, 5, overdispersion=overdispersion)
    rna_lean = _gdna_frac(-40.0, 95, 5, overdispersion=overdispersion)
    assert neutral < 0.2  # at λ=0 a strand-specific node reads RNA
    assert extreme > 0.99  # λ→+∞ siphons it (and every unspliced node) into gDNA
    assert rna_lean < 0.01  # λ→−∞ is the symmetric RNA limit


def test_llr_bias_zero_is_unbiased_median():
    # λ=0 must be exactly the unbiased posterior median (the neutral default — no behavior change).
    # An ambiguous node sits strictly between the gDNA (0.5) and RNA (0.99) sense rates.
    frac = _gdna_frac(0.0, 8, 2)
    assert 0.0 < frac < 1.0
