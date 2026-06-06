"""Strand-deconvolution gDNA confidence (``gdna_strand_confidence_z``) — the FP-aversion knob.

Stresses the per-node joint deconvolution across Beta-Binomial overdispersion levels to be
certain the z-score knob siphons unspliced mass into gDNA monotonically, is neutral at z=0, and
siphons *everything* at the extreme (z→∞ ⇒ all unspliced → gDNA, even a confidently-RNA node).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.joint_deconv import _joint_per_node, strand_confidence_z_to_lambda


def _gdna_frac(z, sense, antisense, *, overdispersion=0.0, rna_sense_frac=0.99, mass=100.0):
    """gDNA fraction of one strand-observable node, count prior flat (count_evidence=0) so the
    strand clue + the z-prior drive the call. Higher z ⇒ more FP-averse (more gDNA)."""
    one = np.array([1.0])
    r = _joint_per_node(
        mass_unspl=np.array([mass]),
        mass_spliced=np.array([0.0]),
        sense=np.array([float(sense)]),
        antisense=np.array([float(antisense)]),
        density=np.array([0.0]),
        count_evidence=np.array([0.0]),  # flat count prior → isolate the strand deconvolution
        eff_len=one,
        strand_observable=np.array([True]),
        rna_sense_frac=rna_sense_frac,
        strand_overdispersion=overdispersion,
        gdna_strand_lambda=strand_confidence_z_to_lambda(z),
        n_grid=400,
    )
    return float(r.gdna_frac[0])


# --- the z -> lambda mapping ------------------------------------------------


def test_z_to_lambda_neutral_and_monotone():
    assert strand_confidence_z_to_lambda(0.0) == 0.0
    zs = [-3.0, -1.0, -0.2, 0.0, 0.2, 1.0, 2.0, 3.0, 5.0]
    lams = [strand_confidence_z_to_lambda(z) for z in zs]
    assert all(b > a for a, b in zip(lams, lams[1:]))  # strictly increasing in z
    assert strand_confidence_z_to_lambda(2.0) > 0.0  # z>0 favors gDNA
    assert -1.0 < strand_confidence_z_to_lambda(-2.0) < 0.0  # z<0 leans RNA, bounded at -1
    assert strand_confidence_z_to_lambda(8.0) > 1e6  # large z → very strong gDNA prior


# --- the knob siphons monotonically, across overdispersion ------------------


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_confidence_z_siphons_monotonically(overdispersion):
    # An ambiguous node (8 sense / 2 antisense, rna_sense=0.99 — the motivating example): raising z
    # monotonically shifts mass into gDNA at every overdispersion level.
    fracs = [_gdna_frac(z, 8, 2, overdispersion=overdispersion) for z in (0.0, 1.0, 2.0, 3.0, 5.0)]
    assert all(b >= a - 1e-9 for a, b in zip(fracs, fracs[1:]))  # non-decreasing in z
    assert fracs[-1] > fracs[0] + 0.1  # and materially higher at high confidence


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_confidence_z_neutral_then_full_siphon(overdispersion):
    # A confidently-RNA node (95 sense / 5 antisense): neutral (z=0) calls it RNA; the extreme
    # (z→∞) siphons ALL of its unspliced mass into gDNA — the property a quantile cannot deliver.
    neutral = _gdna_frac(0.0, 95, 5, overdispersion=overdispersion)
    extreme = _gdna_frac(12.0, 95, 5, overdispersion=overdispersion)
    assert neutral < 0.2  # at z=0 a strand-specific node reads RNA
    assert extreme > 0.99  # z→∞ siphons it (and every unspliced node) into gDNA
