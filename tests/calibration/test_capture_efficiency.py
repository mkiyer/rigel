"""``capture_efficiency``: a piece's capture efficiency is the posterior mean of its clipped gDNA
density against the reference, under the fitted landscape, from its own contained count and the
crossings within a fragment's reach — no floor, no constant.

The two falsification tests were written first and verified failing on the shipped ruler through its
own fixture (`ISSUES: ruler-multimapper-floor-caps-the-correction`): an unprobed exon holding no gDNA
fragment read the multimapper floor ``1/(C+1)``, 0.024 on 40 RNA fragments against a depleted level of
0.001 (+3.19 nat); a transcript of exons shorter than a fragment read 0 on its nil contained supports and
then the floor. Here the posterior reads the depleted level from the population and the tiny exon from
its edge crossings. The rest pins the mechanism's parts: the plug-in limit at high depth, the
apportionment of a crossing to the pieces within reach, a boundary's own efficiency, and calibrate's
wiring (efficiencies exactly 1 with no reference).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.capture_efficiency import capture_efficiencies
from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    crossing_base_shares,
    crossing_eff_length,
)
from rigel.calibration.landscape import fit_landscape
from rigel.calibration.region_arrays import RegionArrays

RHO_REF, RHO_DEP = 1.0, 1e-3


def _pmf(mean=200.0, sd=50.0, lo=100, hi=300) -> np.ndarray:
    w = np.arange(hi + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((w - mean) / sd) ** 2)
    p[(w < lo) | (w > hi)] = 0.0
    return p / p.sum()


PMF = _pmf()
S_E = float(crossing_eff_length(PMF, np.array([UNBOUNDED_REACH]), np.array([UNBOUNDED_REACH]))[0])


def _regions(lengths) -> RegionArrays:
    lengths = np.asarray(lengths, dtype=np.int64)
    b = np.concatenate([[0], np.cumsum(lengths)])
    frame = pd.DataFrame(
        {
            "region_id": np.arange(lengths.size, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * lengths.size, dtype="string"),
            "start": b[:-1],
            "end": b[1:],
            "length": lengths,
            "signature": np.zeros(lengths.size, np.uint8),
        }
    )
    return RegionArrays.from_frame(frame, {"chr1": 0})


def _support(lengths, pmf=PMF) -> np.ndarray:
    """The contained support ``E_f[(ℓ − w + 1)⁺]`` per piece — 0 for a piece shorter than every fragment."""
    w = np.arange(pmf.shape[0], dtype=np.float64)
    L = np.asarray(lengths, dtype=np.float64)[:, None]
    return (pmf[None, :] * np.maximum(L - w[None, :] + 1.0, 0.0)).sum(1)


@pytest.fixture(scope="module")
def landscape():
    """A bimodal gDNA landscape fitted on a population: 900 depleted pieces at 10^-3 and 300 captured
    pieces at the reference, each a Poisson count on a 1 kb support, plus 200 zero-count anchors."""
    rng = np.random.default_rng(3)
    eff = np.full(1400, 1000.0)
    count = np.concatenate(
        [
            rng.poisson(RHO_DEP * 1000.0, 900).astype(float),
            rng.poisson(RHO_REF * 1000.0, 300).astype(float),
            np.zeros(200),
        ]
    )
    anchor = np.r_[np.zeros(1200, bool), np.ones(200, bool)]
    ls = fit_landscape(count, np.maximum(count * 2.0, 1.0), eff, np.full(1400, 0.1), anchor=anchor)
    assert ls is not None
    return ls


def _field(lengths, density, pmf=PMF, crossing_density=None):
    """A deposition-faithful field: contained counts ``ρ_p S_p`` and crossing counts
    ``Σ_q A_eq ρ_q`` (the crossing fragments' bases at their pieces' densities), or an explicit crossing
    density per boundary."""
    ra = _regions(lengths)
    density = np.asarray(density, dtype=np.float64)
    S = _support(lengths, pmf)
    k_reg = density * S
    shares = crossing_base_shares(ra, pmf)
    E, Q, A = shares
    k_bnd = np.zeros(len(lengths) - 1)
    if crossing_density is None:
        np.add.at(k_bnd, E, A * density[Q])
    else:
        k_bnd[:] = np.asarray(crossing_density, dtype=np.float64) * S_E
    return ra, k_reg, S, k_bnd, shares


# --- the two falsification tests ---------------------------------------------------------------------


def test_an_unprobed_exon_with_no_gdna_fragment_reads_the_depleted_level_not_a_floor(landscape):
    """A 1 kb exon holding no gDNA at all among depleted neighbours: its posterior sits at the
    depleted level, within a factor of two of 0.001 and nowhere near the floor's 0.024 or 1."""
    lengths = [5000, 1000, 5000]
    ra, k_reg, S, k_bnd, shares = _field(lengths, [RHO_DEP, RHO_DEP, RHO_DEP])
    k_reg[1] = 0.0
    c, _reach = capture_efficiencies(
        landscape, RHO_REF, k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E), shares
    )
    assert abs(np.log(c[1] / RHO_DEP)) < np.log(2.0), c[1]
    assert c[1] < 0.005


def test_a_tiny_exon_reads_its_edge_crossings(landscape):
    """40 bp exons between 1 kb introns hold no contained fragment (their support is 0). With captured
    crossings at their edges — the fragments there carry the exon's bases at the reference and the
    intron's at the depleted level — each exon reads within a factor of two of fully captured; with
    depleted crossings it reads depleted. PERTURBATION: with the crossing terms dropped every tiny exon
    reads the population's clipped mean, the same number either way."""
    lengths = [5000, 1000, 40, 1000, 40, 1000, 40, 1000, 5000]
    exon = np.array([2, 4, 6])
    density = np.full(len(lengths), RHO_DEP)
    density[exon] = RHO_REF
    ra, k_reg, S, k_bnd, shares = _field(lengths, density)
    assert np.all(S[exon] == 0.0), "the fixture's exons must hold no contained fragment"
    c, _ = capture_efficiencies(
        landscape, RHO_REF, k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E), shares
    )
    assert np.all(c[exon] > 0.5), c[exon]
    ra, k_reg, S, k_bnd, shares = _field(lengths, np.full(len(lengths), RHO_DEP))
    c_dep, _ = capture_efficiencies(
        landscape, RHO_REF, k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E), shares
    )
    assert np.all(c_dep[exon] < 0.02), c_dep[exon]
    assert np.all(c[exon] > 20.0 * c_dep[exon])


# --- the parts ----------------------------------------------------------------------------------------


def test_a_well_measured_piece_reads_its_plug_in(landscape):
    """At high depth the posterior is the plug-in: 50 kb pieces (fifty depleted fragments, thousands of
    captured ones) at the reference read 1, at the depleted level 0.001, and at a third of the reference
    a third."""
    lengths = [50000, 50000, 50000, 50000]
    density = [RHO_REF, RHO_DEP, RHO_REF / 3.0, RHO_REF]
    ra, k_reg, S, k_bnd, shares = _field(lengths, density)
    c, _ = capture_efficiencies(
        landscape, RHO_REF, k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E), shares
    )
    assert c[0] > 0.95 and c[3] > 0.95
    assert abs(np.log(c[1] / RHO_DEP)) < 0.2
    assert abs(np.log(c[2] * 3.0)) < 0.2


def test_every_efficiency_lies_in_the_unit_interval(landscape):
    lengths = [5000, 40, 1000, 7, 9, 300, 5000]
    ra, k_reg, S, k_bnd, shares = _field(lengths, np.full(7, RHO_REF))
    c, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, np.full(6, S_E), shares)
    assert np.all(c >= 0.0) and np.all(c <= 1.0)
    assert cb.shape == (6,) and np.all(cb >= 0.0) and np.all(cb <= 1.0)


def test_a_boundarys_efficiency_is_its_own_crossing_counts_posterior(landscape):
    """A boundary reads its own crossing count on its crossing support: at the reference on both sides
    it reads 1; between two depleted pieces it reads the depleted level; at a probed exon's edge
    beside a depleted intron it reads the crossing fragments' mean — about half. The locus prior reads
    this beside the count the same crossing mass sits in."""
    lengths = [5000, 5000, 5000]
    ra, k_reg, S, k_bnd, shares = _field(lengths, [RHO_REF, RHO_REF, RHO_REF])
    _, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, np.full(2, S_E), shares)
    assert np.all(cb > 0.9)
    ra, k_reg, S, k_bnd, shares = _field(lengths, [RHO_DEP, RHO_DEP, RHO_REF])
    _, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, np.full(2, S_E), shares)
    assert cb[0] < 0.02
    assert 0.3 < cb[1] < 0.7, cb


def test_the_crossing_is_apportioned_to_the_pieces_that_can_explain_it(landscape):
    """Two tiny exons share nothing, but a tiny exon and its long depleted intron share every crossing
    at their edge: the intron's own count pins it depleted, so the captured crossing is the exon's.
    PERTURBATION: attributing the whole count to every piece within reach (no apportionment) lifts the
    depleted intron beside a captured tiny exon well above its own level."""
    lengths = [5000, 1000, 40, 1000, 5000]
    density = np.array([RHO_DEP, RHO_DEP, RHO_REF, RHO_DEP, RHO_DEP])
    ra, k_reg, S, k_bnd, shares = _field(lengths, density)
    c, _ = capture_efficiencies(
        landscape, RHO_REF, k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E), shares
    )
    assert c[2] > 0.5
    assert c[1] < 0.01 and c[3] < 0.01, (c[1], c[3])


def test_calibrate_publishes_the_efficiencies_and_exactly_one_without_a_reference():
    """The wiring: a solve with no refit has no landscape and no reference, and the result's
    efficiencies are exactly 1 on both axes."""
    import sys

    sys.path.insert(0, "tests/calibration")
    from _synthetic import (
        make_gdna_fl_pmf,
        make_strand_models,
        make_synthetic_payload,
        make_synthetic_sj,
    )

    from rigel.calibration import calibrate
    from rigel.config import CalibrationConfig

    payload, ra = make_synthetic_payload()
    pmf = make_gdna_fl_pmf()
    res = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=make_strand_models(0.95, 40),
        gdna_fl_pmf=pmf,
        rna_fl_pmf=pmf,
        config=CalibrationConfig(calib_refit_iters=0),
        sj=make_synthetic_sj(),
    )
    assert res.gdna_reference_density is None
    np.testing.assert_array_equal(res.gdna_capture_efficiency_region, 1.0)
    np.testing.assert_array_equal(res.gdna_capture_efficiency_boundary, 1.0)
