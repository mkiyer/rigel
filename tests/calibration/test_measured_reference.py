"""STAGE 2 OF THE FIRST-PASS REDESIGN — the MEASURED reference location, falsified clause by clause.

`density_deconv.measured_reference_location` gives ψ's reference its per-slot LOCATION from the
measured gDNA background at exactly the stage-0 substrate's single-stranded intron REGIONs:
``m_i = rho_bg·E_g,i / M_i``, with ``rho_bg`` pooled over intergenic REGIONs (the stage-1 anchor pool,
ratio of sums), FIRM-clipped into the lattice window ``[sigma(-L), sigma(L)]`` (owner ruling
2026-08-26: the location clips, the strength stays one pseudo-fragment). Everything else keeps the
base location untouched.

⛔ The scope is REGIONS ONLY, and that is a measured refusal rather than a preference: the
boundary-inclusive form regressed stranded × capture-ON on the test chromosome AND the ladder (the
damage entirely on the boundary axis with messages off), while regions-only improved every stratum on
both — so a silent re-widening to boundaries must fail a test here.

The fixture is constructed arrays (no index): a nine-slot chain whose statics carve one intergenic
REGION pool, one claimed intron, and exon/boundary slots that must pass through — plus a stub
geometry, because the function reads exactly ``unspliced_count`` and ``eff_gdna`` and nothing else.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.density_deconv import measured_reference_location
from rigel.calibration.region_chain import build_region_chain
from rigel.calibration.region_geometry import RegionStatics
from rigel.calibration.splice_graph import FLAG_DONOR_POS, FLAG_TSS_POS


def _sigma(x: float) -> float:
    return 1.0 / (1.0 + np.exp(-x))


def _fixture(counts, eff_gdna):
    """Slots 0..8: N(ig) B(edge) N(exon) B(donor) N(intron) B(acceptor) N(exon) B(edge) N(ig) —
    the same shape the audit's self-test uses, so slot 4 is the ONE claimed intron REGION and slots
    3/5 are its (claimed but out-of-scope) intron BOUNDARIEs."""
    chain = build_region_chain(np.array([0, 5]), np.array([0, 4]))
    fp = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0], bool)
    fn = np.zeros(9, bool)
    mp = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0], bool)
    bflags = np.zeros(9, np.uint16)
    bflags[1] = bflags[7] = FLAG_TSS_POS
    bflags[3] = bflags[5] = FLAG_DONOR_POS
    statics = RegionStatics(
        n_slots=9,
        free_pos=fp,
        free_neg=fn,
        mrna_active_pos=mp,
        mrna_active_neg=np.zeros(9, bool),
        boundary_flags=np.where(np.asarray(chain.kind) == 1, bflags, 0).astype(np.uint16),
    )
    geometry = SimpleNamespace(
        unspliced_count=np.stack([np.asarray(counts, np.float64), np.zeros(9)], axis=1),
        eff_gdna=np.asarray(eff_gdna, np.float64),
    )
    return chain, statics, geometry


BASE = np.full(9, 0.75)


def test_only_the_claimed_intron_region_moves():
    """Exons, boundaries (including the CLAIMED intron boundaries — the refused scope) and intergenic
    slots keep the base location bit-exactly."""
    chain, statics, geometry = _fixture(
        counts=[100, 5, 50, 8, 20, 8, 50, 5, 100], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    moved = np.flatnonzero(out != BASE)
    np.testing.assert_array_equal(moved, [4])


def test_the_location_is_the_background_over_the_slots_own_total():
    """Hand-derived: rho_bg pools the two intergenic REGIONs ONLY (slots 0 and 8) as a ratio of sums —
    never a boundary, never the intron itself — and the intron's location is rho_bg·E_g/M, un-clipped
    when it lands inside the window."""
    counts = np.array([5.0, 999, 50, 999, 20, 999, 50, 999, 5])  # loud boundary counts: the pool
    eff = np.full(9, 1000.0)  # must ignore them or the location moves
    chain, statics, geometry = _fixture(counts, eff)
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    rho_bg = (5.0 + 5.0) / (1000.0 + 1000.0)
    assert out[4] == pytest.approx(rho_bg * 1000.0 / 20.0)  # = 0.25, strictly inside the window


def test_a_collision_clips_to_the_lattice_cap():
    """An intron no denser than background (including an EMPTY one under a live background) claims
    all-gDNA at the CAP sigma(L), never above it — the FIRM ruling's valid-CAP use."""
    chain, statics, geometry = _fixture(
        counts=[100, 0, 50, 0, 2, 0, 50, 0, 100], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    assert out[4] == pytest.approx(_sigma(10.0))
    chain, statics, geometry = _fixture(
        counts=[100, 0, 50, 0, 0, 0, 50, 0, 100], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    assert out[4] == pytest.approx(_sigma(10.0))


def test_zero_background_under_live_mass_clips_to_the_floor():
    """rho_bg = 0 with fragments present: whatever is here is not background — the location claims
    all-RNA at sigma(-L), the zero-gDNA control's truth."""
    chain, statics, geometry = _fixture(
        counts=[0, 0, 50, 0, 20, 0, 50, 0, 0], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    assert out[4] == pytest.approx(_sigma(-10.0))


def test_no_measurement_at_all_keeps_the_base():
    """0/0 — zero background AND zero mass — licenses nothing; the base location survives."""
    chain, statics, geometry = _fixture(
        counts=[0, 0, 50, 0, 0, 0, 50, 0, 0], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    assert out[4] == BASE[4]


def test_the_clip_follows_the_window():
    """The cap is sigma(L) for the WINDOW THE SOLVE RUNS ON, not a constant: L = 5 and L = 10 clip a
    collision to different values (`TRAPS: a-clamp-at-the-closed-end-escapes-the-window`)."""
    chain, statics, geometry = _fixture(
        counts=[100, 0, 50, 0, 1, 0, 50, 0, 100], eff_gdna=np.full(9, 10.0)
    )
    at10 = measured_reference_location(chain, statics, geometry, 10.0, base=BASE.copy())
    at5 = measured_reference_location(chain, statics, geometry, 5.0, base=BASE.copy())
    assert at10[4] == pytest.approx(_sigma(10.0)) and at5[4] == pytest.approx(_sigma(5.0))
    assert at5[4] < at10[4]


def test_no_base_means_neutral_everywhere_else():
    """base=None starts from the neutral half, so the reference term is EXACTLY ZERO at every slot
    the measurement does not reach — the same contract as `structural_reference=False`."""
    chain, statics, geometry = _fixture(
        counts=[100, 5, 50, 8, 20, 8, 50, 5, 100], eff_gdna=np.full(9, 10.0)
    )
    out = measured_reference_location(chain, statics, geometry, 10.0, base=None)
    assert (np.delete(out, 4) == 0.5).all() and out[4] != 0.5
