"""Unit tests for the count module (density_model): raw-count local imputation.

The count module (decoupled architecture) estimates gDNA density from RAW unspliced counts — no
strand cleaning. Observable regions use their own contained density; non-observable (exon/AMBIG)
regions are anchored from their observable boundary sides (crossing count / fl_mean); run interiors
are carried inward; a region with no observable node anywhere takes the global count-weighted-mean
observable density. The per-node gDNA fraction is ``count_gdna_frac = clip(density·eff / mass)``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
)
from rigel.calibration.substrate import CalibrationSubstrate


def SubstrateView(*_args, **_kwargs):
    """⛔ Deleted in S5.d — the 4-channel per-face view is gone.

    It held `[unspliced+, unspliced-, spliced_sense, spliced_antisense]`, i.e. TWO strand conventions in
    one array, and a `mass` that the accumulator no longer emits. Its successor is
    `substrate.PopulationView`: three integer sums, genome strand only.

    Bound so this module still IMPORTS: `node_gdna_density` is an S5.e consumer and its tests fail there
    with a message naming the step, instead of an ImportError taking out the whole file. S5.e deletes it.
    """
    raise NotImplementedError(
        "SubstrateView was deleted in S5.d; use substrate.PopulationView. node_gdna_density is S5.e. "
        "See docs/S5_DESIGN_LOG.md §2."
    )


INTRON = BIT_INTRON_POS  # 0x8 — intron+, count-observable
EXON = BIT_EXON_POS  # 0x2 — exon+, NOT count-observable, strand POS
AMBIG = BIT_EXON_POS | BIT_EXON_NEG  # 0x3 — exon+ over exon−, NOT observable, strand AMBIG
INTERGENIC = 0  # no bits — count-observable


# --------------------------------------------------------------------------- #
# helpers to build a minimal substrate + region geometry
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
    return RegionArrays.from_frame(df, refmap)


def _substrate(n, contained, left, right) -> CalibrationSubstrate:
    return CalibrationSubstrate(
        n_regions=n,
        strand_class=np.zeros(n, dtype=np.int8),  # node uses region_arrays.strand_class
        contained=_view(*contained),
        left=_view(*left),
        right=_view(*right),
    )


def _zeros(n):
    return (np.zeros(n), np.zeros(n))


def _density(sub, ra, eff, fl_mean=50.0):
    return node_gdna_density(
        sub, ra, region_eff_len=np.asarray(eff, dtype=np.float64), fl_mean=fl_mean
    )


# --------------------------------------------------------------------------- #
# node_gdna_density — geometry cases
# --------------------------------------------------------------------------- #


def test_exon_between_introns_recovers_uniform_density():
    # intron+ | exon+ | intron+ ; both exon boundaries are observable.
    ra = _region_arrays([INTRON, EXON, INTRON])
    contained = ([100, 0, 100], [100, 0, 100])  # observable introns: contained count 200 each
    left = ([0, 50, 0], [0, 50, 0])  # exon's boundary sides: crossing flux total 100 each
    right = ([0, 50, 0], [0, 50, 0])
    sub = _substrate(3, contained, left, right)
    nd = _density(sub, ra, [100.0, 100.0, 100.0])
    # ρ = 200/100 = 2 (introns) and 100/50 = 2 (exon from each side) → uniform 2.
    assert np.allclose(nd.density, [2.0, 2.0, 2.0])


def test_tiny_observable_region_anchors_from_boundaries():
    # intergenic | tiny intron (eff_len 0) | intergenic — the tiny region cannot use contained.
    ra = _region_arrays([INTERGENIC, INTRON, INTERGENIC])
    contained = ([100, 0, 100], [100, 0, 100])
    left = ([0, 30, 0], [0, 30, 0])  # r1 left side total 60 → 60/50 = 1.2
    right = ([0, 40, 0], [0, 40, 0])  # r1 right side total 80 → 80/50 = 1.6
    sub = _substrate(3, contained, left, right)
    nd = _density(sub, ra, [100.0, 0.0, 100.0])
    assert np.isfinite(nd.density[1])  # not inf from /eff_len=0
    assert nd.density[1] == pytest.approx(np.mean([1.2, 1.6]))


def test_run_interior_filled_from_anchored_edges():
    # intron+ | exon+ | exon+ | exon+ | intron+ : the middle exon has BOTH boundaries
    # non-observable (shared exon bit) and is reachable only by the inward carry.
    ra = _region_arrays([INTRON, EXON, EXON, EXON, INTRON])
    contained = ([100, 0, 0, 0, 100], [100, 0, 0, 0, 100])
    left = ([0, 50, 0, 0, 0], [0, 50, 0, 0, 0])  # r1 anchors from its (observable) left boundary
    right = (
        [0, 0, 0, 100, 0],
        [0, 0, 0, 100, 0],
    )  # r3 anchors from its (observable) right boundary
    sub = _substrate(5, contained, left, right)
    nd = _density(sub, ra, np.full(5, 100.0))
    assert nd.density[1] == pytest.approx(2.0)  # 100/50
    assert nd.density[3] == pytest.approx(4.0)  # 200/50
    assert nd.density[2] == pytest.approx(3.0)  # carry mean of the two flanks


def test_count_gdna_frac_is_density_ratio():
    # count_gdna_frac = clip(density·eff / mass). For a pure-gDNA observable intron at uniform
    # density the contained mass equals density·eff, so the fraction is 1 (all gDNA).
    ra = _region_arrays([INTRON, EXON, INTRON])
    contained = ([100, 0, 100], [100, 0, 100])
    left = ([0, 50, 0], [0, 50, 0])
    right = ([0, 50, 0], [0, 50, 0])
    sub = _substrate(3, contained, left, right)
    nd = _density(sub, ra, np.full(3, 100.0))
    assert nd.count_gdna_frac[0] == pytest.approx(1.0)  # intron: density·eff = 2·100 = 200 = mass
    assert nd.count_gdna_frac[2] == pytest.approx(1.0)
    # exon has no contained mass here ⇒ fraction 0 (no gDNA to attribute from its own count)
    assert nd.count_gdna_frac[1] == pytest.approx(0.0)


def test_no_observable_node_takes_zero_baseline():
    # A single exon-only reference: no observable region, no observable boundary → no anchor.
    # Density takes the global baseline, which is 0 (there is no observable region anywhere).
    ra = _region_arrays([EXON])
    sub = _substrate(1, _zeros(1), _zeros(1), _zeros(1))
    nd = _density(sub, ra, [100.0])
    assert nd.density[0] == 0.0
    assert nd.count_gdna_frac[0] == 0.0


def test_density_does_not_cross_references():
    # chr1 carries a high-density observable intron (density 4.0); chr2 is a lone no-anchor exon.
    # The chr2 exon must NOT inherit chr1's 4.0 via the run-fill carry (the carry is per-reference)
    # — it takes the GLOBAL baseline (the count-weighted-mean observable density, 4.0 here).
    ra = _region_arrays([INTRON, EXON], ref_names=["chr1", "chr2"])
    contained = ([200, 0], [200, 0])
    sub = _substrate(2, contained, _zeros(2), _zeros(2))
    nd = _density(sub, ra, [100.0, 100.0])
    assert nd.density[0] == pytest.approx(4.0)  # chr1 observable intron: 400 / 100
    assert nd.density[1] == pytest.approx(4.0)  # chr2 no-anchor: GLOBAL baseline, not a carry
