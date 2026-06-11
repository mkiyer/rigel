"""Per-node strand/count handoff (`strand_deconv._deconv_per_node`).

Two concerns: (1) the **routing** — a node takes the strand posterior iff the library is
strand-identifiable AND the node is strand-observable, else the count fraction; (2) the **FP-rate
quantile knob** (``deconv_quantile``) — reads ``g(q)=clip(center+Φ⁻¹(q)·σ)``: neutral at ½ (no shift),
shifting toward gDNA as q rises. Crucially it is **uncertainty-aware** (the shift scales with the
per-node σ), so it widens the estimate but cannot siphon a *confident* node — the deliberate contrast
with the retired fixed log-odds tilt (see ``docs/calibration/phase2_design.md``).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.strand_deconv import _deconv_per_node


def _frac(
    *, strand_observable=True, sense, antisense, count_gdna_frac, count_gdna_frac_var=0.0,
    quantile=0.5, overdispersion=0.0, rna_overdispersion=0.0, rna_sense_frac=0.99, mass=100.0,
):
    r = _deconv_per_node(
        np.array([mass]),
        np.array([0.0]),
        np.array([float(sense)]),
        np.array([float(antisense)]),
        np.array([bool(strand_observable)]),
        np.array([float(count_gdna_frac)]),
        np.array([float(count_gdna_frac_var)]),
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=overdispersion,
        rna_strand_overdispersion=rna_overdispersion,
        deconv_quantile=quantile,
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
# FP-rate quantile knob (uncertainty-aware shift on the final blended fraction)
# --------------------------------------------------------------------------- #


def test_quantile_half_is_noop_center():
    # q=½ ⇒ Φ⁻¹=0 ⇒ g = center (the blended point estimate), independent of σ: the variance is
    # consumed only when q≠½. A count variance present must NOT change the q=½ answer.
    no_var = _frac(sense=12, antisense=8, count_gdna_frac=0.5, rna_sense_frac=0.99, quantile=0.5)
    with_var = _frac(sense=12, antisense=8, count_gdna_frac=0.5, count_gdna_frac_var=0.05,
                     rna_sense_frac=0.99, quantile=0.5)
    assert no_var == pytest.approx(with_var)


@pytest.mark.parametrize("overdispersion", [0.0, 0.1, 0.5])
def test_quantile_monotone_toward_gdna(overdispersion):
    # g(q) is non-decreasing in q: q>½ is FP-averse (more gDNA), q<½ leans RNA. Use a wide-posterior
    # node (mid, near-balanced counts) so the quantile has room to move.
    fracs = [_frac(sense=6, antisense=4, count_gdna_frac=0.5, count_gdna_frac_var=0.02,
                   rna_sense_frac=0.99, overdispersion=overdispersion, quantile=q)
             for q in (0.05, 0.25, 0.5, 0.75, 0.95)]
    assert all(b >= a - 1e-9 for a, b in zip(fracs, fracs[1:]))  # non-decreasing in q
    assert fracs[-1] > fracs[2] + 0.01  # q=0.95 materially more gDNA than the q=0.5 center


def test_quantile_is_uncertainty_aware():
    # The shift scales with σ. On a count-routed node σ=σ_count, so a higher count variance ⇒ a
    # larger q=½→q=0.95 shift. This is the property the retired fixed log-odds tilt lacked.
    def shift(var):
        center = _frac(strand_observable=False, sense=0, antisense=0, count_gdna_frac=0.5,
                       count_gdna_frac_var=var, quantile=0.5)
        hi = _frac(strand_observable=False, sense=0, antisense=0, count_gdna_frac=0.5,
                   count_gdna_frac_var=var, quantile=0.95)
        return hi - center
    assert shift(0.04) > shift(0.0025) + 0.05  # wider count posterior ⇒ larger FP-quantile shift


def test_quantile_cannot_siphon_a_confident_node():
    # The decisive contrast with the retired log-odds tilt (which forced even a confident-RNA node to
    # gDNA at λ→+∞): a strand-specific RNA node has σ→0, so an extreme q barely moves it — by design.
    rna_center = _frac(sense=99, antisense=1, count_gdna_frac=0.0, rna_sense_frac=0.99, quantile=0.5)
    rna_extreme = _frac(sense=99, antisense=1, count_gdna_frac=0.0, rna_sense_frac=0.99, quantile=0.99)
    assert rna_center < 0.1  # confidently RNA at the center
    assert rna_extreme < 0.3  # NOT siphoned into gDNA (σ is tiny) — a quantile cannot move it


def test_quantile_clips_to_unit_interval():
    # g(q)=clip(center+Φ⁻¹(q)·σ) stays in [0, 1] at extreme quantiles + large σ.
    hi = _frac(strand_observable=False, sense=0, antisense=0, count_gdna_frac=0.9,
               count_gdna_frac_var=0.09, quantile=0.999)
    lo = _frac(strand_observable=False, sense=0, antisense=0, count_gdna_frac=0.1,
               count_gdna_frac_var=0.09, quantile=0.001)
    assert 0.0 <= lo <= hi <= 1.0
