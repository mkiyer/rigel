"""Unit tests for the count module (density_model): raw-count local imputation.

The count module estimates gDNA density from RAW unspliced counts — no strand cleaning. Observable
regions use their own contained density; non-observable (exon/AMBIG) regions are anchored from their
observable flanking BOUNDARIES; run interiors are carried inward; a region with no observable neighbour
anywhere takes the global count-weighted-mean observable density.

⭐ **What S5.e changed here.** A region used to carry its own ``left``/``right`` per-side views of the
flanking boundary flux, and the crossing divisor arrived separately as a scalar ``fl_mean``. An boundary is
a 0-bp boundary with ONE count, shared by the two regions it separates, and its divisor is the same
``eff_gdna`` array the solver uses — so both the side views and the second length model are gone, and a
caller can no longer hand this a divisor that disagrees with the geometry's.

⚠ The geometry here is built BY HAND rather than through ``build_region_geometry``: these tests are about
the imputation logic, and the divisors have their own gate in ``test_region_geometry.py``. Stating the
numbers directly is what keeps the expected densities arithmetic rather than a second derivation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.density_model import region_gdna_density
from rigel.calibration.region_chain import REGION, build_region_chain
from rigel.calibration.region_geometry import RegionGeometry
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
)


INTRON = BIT_INTRON_POS  # 0x8 — intron+, count-observable
EXON = BIT_EXON_POS  # 0x2 — exon+, NOT count-observable, strand POS
AMBIG = BIT_EXON_POS | BIT_EXON_NEG  # 0x3 — exon+ over exon−, NOT observable, strand AMBIG
INTERGENIC = 0  # no bits — count-observable


# --------------------------------------------------------------------------- #
# helpers to build a minimal substrate + region geometry
# --------------------------------------------------------------------------- #


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
            "node_id": np.arange(n, dtype=np.int64),
            "ref_name": pd.array(ref_names, dtype="string"),
            "start": np.asarray(starts, dtype=np.int64),
            "end": np.asarray(ends, dtype=np.int64),
            "length": np.full(n, 100, dtype=np.int64),
            "signature": np.asarray(signatures, dtype=np.uint8),
        }
    )
    refmap = {rn: i for i, rn in enumerate(dict.fromkeys(ref_names))}
    return RegionArrays.from_frame(df, refmap)


def _parts(signatures, region_count, region_eff, boundary_count, boundary_eff, ref_names=None):
    """A chain + a hand-built :class:`RegionGeometry` over ``signatures``.

    ``region_count``/``region_eff`` are per REGION, ``boundary_count``/``boundary_eff`` per contiguous BOUNDARY — and a
    reference with ``k`` regions owns exactly ``k - 1`` boundaries, with no terminal slots.
    """
    ra = _region_arrays(signatures, ref_names)
    rno = np.asarray(ra.ref_offsets, dtype=np.int64)
    reo = np.zeros_like(rno)
    np.cumsum(np.maximum(np.diff(rno) - 1, 0), out=reo[1:])
    chain = build_region_chain(rno, reo)

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    n_slots = chain.n_slots
    count = np.zeros((n_slots, 2))
    eff = np.zeros(n_slots)
    count[is_region, 0] = np.asarray(region_count, float)[obj[is_region]]
    eff[is_region] = np.asarray(region_eff, float)[obj[is_region]]
    if np.any(~is_region):
        count[~is_region, 0] = np.asarray(boundary_count, float)[obj[~is_region]]
        eff[~is_region] = np.asarray(boundary_eff, float)[obj[~is_region]]
    geometry = RegionGeometry(
        n_slots=n_slots,
        unspliced_count=count,
        # the length channels are irrelevant to region_gdna_density (it reads count and eff only);
        # shape-correct zeros say "this fixture makes no statement about fragment lengths".
        spliced_count=np.zeros((n_slots, 2)),
        eff_gdna=eff,
        eff_rna=eff,
        junction_count=np.zeros((n_slots, 2)),
        eff_junction=np.zeros((n_slots, 2)),
        # the per-FLANK split of the same flux; this fixture has no junction at all, so all four are 0
        # and the two flank totals coincide (`region_total_density`).
        junction_count_lo=np.zeros((n_slots, 2)),
        junction_count_hi=np.zeros((n_slots, 2)),
        eff_junction_lo=np.zeros((n_slots, 2)),
        eff_junction_hi=np.zeros((n_slots, 2)),
    )
    return region_gdna_density(chain, geometry, ra)


# --------------------------------------------------------------------------- #
# region_gdna_density — geometry cases
# --------------------------------------------------------------------------- #


def test_exon_between_introns_recovers_uniform_density():
    # intron+ | exon+ | intron+ ; both of the exon's flanking boundaries are observable.
    nd = _parts(
        [INTRON, EXON, INTRON],
        region_count=[200, 0, 200],
        region_eff=[100.0, 100.0, 100.0],
        boundary_count=[100, 100],
        boundary_eff=[50.0, 50.0],
    )
    # rho = 200/100 = 2 (introns) and 100/50 = 2 (exon from each boundary) -> uniform 2.
    assert np.allclose(nd.density, [2.0, 2.0, 2.0])


def test_tiny_observable_region_anchors_from_boundaries():
    # intergenic | tiny intron (eff_len 0) | intergenic — the tiny region cannot use contained.
    # ⚠ eff 0 must make the region fall through to its flanks, NOT divide and floor (trap 23).
    nd = _parts(
        [INTERGENIC, INTRON, INTERGENIC],
        region_count=[200, 0, 200],
        region_eff=[100.0, 0.0, 100.0],
        boundary_count=[60, 80],
        boundary_eff=[50.0, 50.0],
    )
    assert np.isfinite(nd.density[1])  # not inf from /eff_len = 0
    assert nd.density[1] == pytest.approx(np.mean([1.2, 1.6]))


def test_run_interior_filled_from_anchored_boundaries():
    # intron+ | exon+ | exon+ | exon+ | intron+ : the middle exon has BOTH boundaries non-observable
    # (shared exon bit) and is reachable only by the inward carry.
    nd = _parts(
        [INTRON, EXON, EXON, EXON, INTRON],
        region_count=[200, 0, 0, 0, 200],
        region_eff=np.full(5, 100.0),
        boundary_count=[100, 0, 0, 200],
        boundary_eff=np.full(4, 50.0),
    )
    assert nd.density[1] == pytest.approx(2.0)  # 100/50, from its observable LEFT boundary
    assert nd.density[3] == pytest.approx(4.0)  # 200/50, from its observable RIGHT boundary
    assert nd.density[2] == pytest.approx(3.0)  # carry mean of the two flanks


def test_count_gdna_frac_is_density_ratio():
    # count_gdna_frac = clip(density·eff / contained count). For a pure-gDNA observable intron at
    # uniform density the contained count equals density·eff, so the fraction is 1 (all gDNA).
    nd = _parts(
        [INTRON, EXON, INTRON],
        region_count=[200, 0, 200],
        region_eff=np.full(3, 100.0),
        boundary_count=[100, 100],
        boundary_eff=[50.0, 50.0],
    )
    assert nd.count_gdna_frac[0] == pytest.approx(1.0)  # density·eff = 2·100 = 200 = the count
    assert nd.count_gdna_frac[2] == pytest.approx(1.0)
    # the exon has no contained count here => fraction 0 (no gDNA to attribute from its own count)
    assert nd.count_gdna_frac[1] == pytest.approx(0.0)


def test_no_observable_region_takes_zero_baseline():
    # A single exon-only reference: no observable region, and with one region there is no boundary at all.
    # Density takes the global baseline, which is 0 (there is no observable region anywhere).
    nd = _parts([EXON], region_count=[0], region_eff=[100.0], boundary_count=[], boundary_eff=[])
    assert nd.density[0] == 0.0
    assert nd.count_gdna_frac[0] == 0.0


def test_density_does_not_cross_references():
    # chr1 carries a high-density observable intron (density 4.0); chr2 is a lone no-anchor exon.
    # The chr2 exon must NOT inherit chr1's 4.0 via the run-fill carry (the carry is per-reference)
    # — it takes the GLOBAL baseline (the count-weighted-mean observable density, 4.0 here).
    # ⭐ Each reference owns ONE region and therefore ZERO boundaries, so nothing can leak across even in
    # principle: `chain.left`/`chain.right` are -1 at every reference terminal.
    nd = _parts(
        [INTRON, EXON],
        region_count=[400, 0],
        region_eff=[100.0, 100.0],
        boundary_count=[],
        boundary_eff=[],
        ref_names=["chr1", "chr2"],
    )
    assert nd.density[0] == pytest.approx(4.0)  # chr1 observable intron: 400 / 100
    assert nd.density[1] == pytest.approx(4.0)  # chr2 no-anchor: GLOBAL baseline, not a carry
