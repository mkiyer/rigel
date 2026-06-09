"""Unit tests for the density model: strand-clean closed form + local imputation.

Phase 1 of the density rework (docs/calibration/density_phase1_local_imputation_design.md): the
global sweep is replaced by a LOCAL estimator — observable regions use their own contained density,
non-observable (exon) regions are anchored from their observable boundary sides (crossing count /
fl_mean, unbiased post the accumulator span fix), run interiors are carried inward, and a region
with no observable node anywhere takes density 0 (⇒ defer to strand).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.density_model import (
    node_gdna_density,
    strand_clean_gdna_frac,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_INTRON_POS,
)
from rigel.calibration.substrate import CalibrationSubstrate, SubstrateView

INTRON = BIT_INTRON_POS  # 0x8 — intron+, count-observable
EXON = BIT_EXON_POS  # 0x2 — exon+, NOT count-observable
INTERGENIC = 0  # no bits — count-observable


# --------------------------------------------------------------------------- #
# strand_clean_gdna_frac
# --------------------------------------------------------------------------- #


def test_clean_rna_rate_reads_zero_gdna():
    # A node exactly at the RNA sense rate κ is pure RNA → gdna_frac 0.
    gf = strand_clean_gdna_frac(sense=[80.0], total=[100.0], rna_sense_frac=0.8)
    assert gf[0] == pytest.approx(0.0)


def test_clean_symmetric_reads_full_gdna():
    # Symmetric (sense_frac ½) = gDNA's signature → gdna_frac 1.
    gf = strand_clean_gdna_frac(sense=[50.0], total=[100.0], rna_sense_frac=0.8)
    assert gf[0] == pytest.approx(1.0)


def test_clean_clips_to_unit_interval():
    # More-antisense-than-RNA and more-sense-than-gDNA both saturate (no negatives, no >1).
    gf = strand_clean_gdna_frac(sense=[20.0, 95.0], total=[100.0, 100.0], rna_sense_frac=0.8)
    assert gf[0] == pytest.approx(1.0)  # (0.2−0.8)/(0.5−0.8)=2 → clip 1
    assert gf[1] == pytest.approx(0.0)  # (0.95−0.8)/(0.5−0.8)=−0.5 → clip 0


def test_clean_unstranded_is_degenerate_one():
    # κ ≈ ½ ⇒ gDNA and RNA both symmetric ⇒ cannot clean ⇒ 1.0 (keep raw), never a divide-by-zero.
    gf = strand_clean_gdna_frac(sense=[10.0, 90.0], total=[100.0, 100.0], rna_sense_frac=0.5)
    assert np.allclose(gf, 1.0)


def test_clean_empty_node_is_finite():
    gf = strand_clean_gdna_frac(sense=[0.0], total=[0.0], rna_sense_frac=0.9)
    assert np.isfinite(gf[0])


# --------------------------------------------------------------------------- #
# node_gdna_density — helpers to build a minimal substrate + region geometry
# --------------------------------------------------------------------------- #


def _view(pos, neg) -> SubstrateView:
    pos = np.asarray(pos, dtype=np.int64)
    neg = np.asarray(neg, dtype=np.int64)
    z = np.zeros_like(pos)
    return SubstrateView(
        n_unspliced_pos=pos,
        n_unspliced_neg=neg,
        n_spliced_sense=z,
        n_spliced_antisense=z,
        mass_unspliced=(pos + neg).astype(np.float64),
        mass_spliced=z.astype(np.float64),
    )


def _region_arrays(signatures, ref_names=None) -> RegionArrays:
    n = len(signatures)
    ref_names = ref_names or ["chr1"] * n
    starts, ends, pos, cur = [], [], 0, None
    for rn in ref_names:
        if rn != cur:
            pos, cur = 0, rn
        starts.append(pos)
        ends.append(pos + 100)
        pos += 100
    df = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": pd.array(ref_names, dtype="string"),
            "start": np.asarray(starts, dtype=np.int64),
            "end": np.asarray(ends, dtype=np.int64),
            "length": np.full(n, 100, dtype=np.int64),
            "signature": np.asarray(signatures, dtype=np.uint8),
        }
    )
    refmap = {rn: i for i, rn in enumerate(dict.fromkeys(ref_names))}
    return RegionArrays.from_region_df(df, refmap)


def _substrate(n, contained, left, right) -> CalibrationSubstrate:
    return CalibrationSubstrate(
        n_regions=n,
        region_len=np.full(n, 100.0),
        strand_class=np.zeros(n, dtype=np.int8),  # node uses region_arrays.strand_class
        contained=_view(*contained),
        left=_view(*left),
        right=_view(*right),
    )


def _zeros(n):
    return (np.zeros(n), np.zeros(n))


# --------------------------------------------------------------------------- #
# node_gdna_density — geometry cases
# --------------------------------------------------------------------------- #


def test_exon_between_introns_recovers_uniform_density():
    # intron+ | exon+ | intron+ ; both exon boundaries are observable.
    ra = _region_arrays([INTRON, EXON, INTRON])
    # observable introns: symmetric (pure-gDNA) contained, count 200 each.
    contained = ([100, 0, 100], [100, 0, 100])
    # exon's two boundary sides: symmetric crossing flux, total 100 each.
    left = ([0, 50, 0], [0, 50, 0])
    right = ([0, 50, 0], [0, 50, 0])
    sub = _substrate(3, contained, left, right)
    nd = node_gdna_density(
        sub, ra, region_eff_len=np.array([100.0, 100.0, 100.0]), fl_mean=50.0, rna_sense_frac=0.99
    )
    # ρ = 200/100 = 2 (introns) and 100/50 = 2 (exon from each side) → uniform 2.
    assert np.allclose(nd.density, [2.0, 2.0, 2.0])


def test_tiny_observable_region_anchors_from_boundaries():
    # intergenic | tiny intron (eff_len 0) | intergenic — the tiny region cannot use contained.
    ra = _region_arrays([INTERGENIC, INTRON, INTERGENIC])
    contained = ([100, 0, 100], [100, 0, 100])
    left = ([0, 30, 0], [0, 30, 0])  # r1 left side total 60 → 60/50 = 1.2
    right = ([0, 40, 0], [0, 40, 0])  # r1 right side total 80 → 80/50 = 1.6
    sub = _substrate(3, contained, left, right)
    nd = node_gdna_density(
        sub, ra, region_eff_len=np.array([100.0, 0.0, 100.0]), fl_mean=50.0, rna_sense_frac=0.99
    )
    assert np.isfinite(nd.density[1])  # not inf from /eff_len=0
    assert nd.density[1] == pytest.approx(np.mean([1.2, 1.6]))


def test_run_interior_filled_from_anchored_edges():
    # intron+ | exon+ | exon+ | exon+ | intron+ : the middle exon has BOTH boundaries
    # non-observable (shared exon bit) and is reachable only by the inward carry.
    ra = _region_arrays([INTRON, EXON, EXON, EXON, INTRON])
    contained = ([100, 0, 0, 0, 100], [100, 0, 0, 0, 100])
    left = ([0, 50, 0, 0, 0], [0, 50, 0, 0, 0])  # r1 anchors from its (observable) left boundary
    right = ([0, 0, 0, 100, 0], [0, 0, 0, 100, 0])  # r3 anchors from its (observable) right boundary
    sub = _substrate(5, contained, left, right)
    nd = node_gdna_density(
        sub, ra, region_eff_len=np.full(5, 100.0), fl_mean=50.0, rna_sense_frac=0.99
    )
    assert nd.density[1] == pytest.approx(2.0)  # 100/50
    assert nd.density[3] == pytest.approx(4.0)  # 200/50
    assert nd.density[2] == pytest.approx(3.0)  # carry mean of the two flanks


def test_no_observable_node_defers_to_strand_with_zero_density():
    # A single exon-only reference: no observable region, no observable boundary → density 0
    # ⇒ count prior collapses to Jeffreys ⇒ strand governs (never a deflated global average).
    ra = _region_arrays([EXON])
    sub = _substrate(1, _zeros(1), _zeros(1), _zeros(1))
    nd = node_gdna_density(
        sub, ra, region_eff_len=np.array([100.0]), fl_mean=50.0, rna_sense_frac=0.99
    )
    assert nd.density[0] == 0.0
    assert nd.count_evidence[0] == 0.0


def test_density_does_not_cross_references():
    # Two single-region refs: an observable intron and a no-anchor exon — the exon must NOT
    # inherit the intron's density (the carry is per-reference).
    ra = _region_arrays([INTRON, EXON], ref_names=["chr1", "chr2"])
    contained = ([100, 0], [100, 0])
    sub = _substrate(2, contained, _zeros(2), _zeros(2))
    nd = node_gdna_density(
        sub, ra, region_eff_len=np.array([100.0, 100.0]), fl_mean=50.0, rna_sense_frac=0.99
    )
    assert nd.density[0] == pytest.approx(2.0)
    assert nd.density[1] == 0.0  # no cross-ref carry
