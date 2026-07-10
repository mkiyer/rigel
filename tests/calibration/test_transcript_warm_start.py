"""Unit tests for the Phase-B warm-start builder (`build_transcript_warm_start` + `_capture_correction`).

Two concerns, tested in isolation:

  * `_capture_correction` — the per-node capture lift `rho·(1 + w·(1/ε − 1))`: the no-reference identity,
    the depleted-mode ε floor, and the evidence shrinkage `w = ev/(ev+1)`.
  * `build_transcript_warm_start` — the bottleneck ROLE logic (contained / crossing / junction), strand
    selection, observability exclusion, junction strand-gated discrimination, and the NaN fallback. Driven
    through a controlled mock incidence with `rho_ref = None` (< 5 gDNA nodes) so the identity correction
    isolates the min/strand/observability behaviour from the correction math.
"""

from __future__ import annotations

import numpy as np
import pytest
from types import SimpleNamespace

import rigel.calibration.capture_eff_length as cel
from rigel.calibration.capture_eff_length import _capture_correction, build_transcript_warm_start
from rigel.calibration.result import RnaWarmStart
from rigel.types import STRAND_NEG, STRAND_POS


# ---------------------------------------------------------------------------
# _capture_correction
# ---------------------------------------------------------------------------

def test_capture_correction_identity_without_reference():
    correct = _capture_correction(None, None)
    rho = np.array([5.0, 0.0, 9.0])
    np.testing.assert_array_equal(correct(rho, np.array([1.0, 1.0, 1.0]), np.array([9.0, 9.0, 9.0])), rho)


def test_capture_correction_lift_floor_and_shrinkage():
    correct = _capture_correction(rho_ref=10.0, rho_depleted=2.0)  # ε_floor = 0.2
    rho = np.array([1.0, 1.0, 1.0, 1.0])
    gdna = np.array([10.0, 1.0, 1.0, 5.0])  # ε = 1.0 (enriched) / 0.1→floored 0.2 / 0.1→0.2 / 0.5
    ev = np.array([1e6, 1e6, 0.0, 1e6])     # w ≈ 1 / 1 / 0 / 1
    out = correct(rho, gdna, ev)
    assert out[0] == pytest.approx(1.0)          # fully enriched (ε=1) → no lift
    assert out[1] == pytest.approx(1.0 * 5.0)    # ε floored at 0.2 → 1/ε = 5, full evidence → ×5
    assert out[2] == pytest.approx(1.0)          # zero evidence → w=0 → no lift despite low ε
    assert out[3] == pytest.approx(1.0 * 2.0)    # ε=0.5 → 1/ε = 2, full evidence → ×2


# ---------------------------------------------------------------------------
# build_transcript_warm_start — fixtures
# ---------------------------------------------------------------------------

def _rna_warm_start() -> RnaWarmStart:
    i8 = lambda v: np.asarray(v, dtype=np.int8)  # noqa: E731
    return RnaWarmStart(
        rho_contained_pos=np.array([5.0, 2.0, 8.0, 0.01]),
        rho_contained_neg=np.array([1.0, 3.0, 1.0, 0.01]),
        rho_crossing_pos=np.array([3.0, 7.0, 0.0, 0.0]),
        rho_crossing_neg=np.array([9.0, 9.0, 0.0, 0.0]),
        rho_spliced_right=np.array([0.0, 6.0, 0.0, 0.0]),  # donor at region 1's right
        rho_spliced_left=np.array([0.0, 0.0, 4.0, 0.0]),   # acceptor at region 2's left
        spliced_strand_right=i8([0, STRAND_POS, 0, 0]),
        spliced_strand_left=i8([0, 0, STRAND_POS, 0]),
    )


def _calibration():
    """4 regions, one reference. gDNA present on regions 0-2 only ⇒ < 5 valid nodes ⇒ rho_ref = None ⇒
    identity correction. Region 3 has NO observed mass (contained_ev = 0) → excluded. Seams 0,1 carry gDNA
    (observable); seam 2 does not."""
    return SimpleNamespace(
        mass_gdna_contained=np.array([1.0, 1.0, 1.0, 0.0]),
        gdna_region_eff_len=np.ones(4),
        mass_rna_contained=np.array([4.0, 4.0, 4.0, 0.0]),   # contained_ev = [5, 5, 5, 0]
        mass_gdna_right=np.array([1.0, 1.0, 0.0, 0.0]),
        mass_gdna_left=np.array([0.0, 1.0, 1.0, 0.0]),        # seam gDNA = [2, 2, 0, 0]
        mass_rna_right=np.zeros(4),
        mass_rna_left=np.zeros(4),
        gdna_boundary_len=np.ones(4),
        rna_warm_start=_rna_warm_start(),
    )


def _run(monkeypatch, *, rt, rr, bt, br, jt, jl, jr, strand, eff, n_t):
    e = np.empty(0, dtype=np.int64)
    inc = tuple(np.asarray(a, dtype=np.int64) if len(a) else e for a in (rt, rr, bt, br, jt, jl, jr))
    monkeypatch.setattr(cel, "_transcript_node_incidence", lambda index, region_arrays: inc)
    index = SimpleNamespace(num_transcripts=n_t, t_to_strand_arr=np.asarray(strand))
    region_arrays = SimpleNamespace(ref_id=np.zeros(4, dtype=np.int64))
    return build_transcript_warm_start(
        _calibration(), region_arrays, index, np.asarray(eff, dtype=np.float64)
    )


# ---------------------------------------------------------------------------
# build_transcript_warm_start — role logic
# ---------------------------------------------------------------------------

def test_contained_bottleneck_strand_and_observability(monkeypatch):
    # t0 (+) over regions 0-3: observable {0,1,2} pos densities {5,2,8} → min 2 (region 3 excluded despite
    #   its 0.01 density — proves unobserved nodes never crush). t1 (−) over {0,1}: neg {1,3} → min 1.
    # t2 not in the incidence → NaN (fall back to the coverage seed).
    warm = _run(
        monkeypatch,
        rt=[0, 0, 0, 0, 1, 1], rr=[0, 1, 2, 3, 0, 1],
        bt=[], br=[], jt=[], jl=[], jr=[],
        strand=[STRAND_POS, STRAND_NEG, STRAND_POS], eff=[10.0, 20.0, 30.0], n_t=3,
    )
    np.testing.assert_allclose(warm[:2], [2.0 * 10.0, 1.0 * 20.0])
    assert np.isnan(warm[2])


def test_crossing_seam_bottleneck_and_observability(monkeypatch):
    # t0 (+) over seams {0,1,2}: crossing pos {3,7,0}; seam 2 unobservable (seam gDNA 0) → excluded → min 3
    #   (not 0). Proves the seam observability gate.
    warm = _run(
        monkeypatch,
        rt=[], rr=[], bt=[0, 0, 0], br=[0, 1, 2], jt=[], jl=[], jr=[],
        strand=[STRAND_POS], eff=[10.0], n_t=1,
    )
    assert warm[0] == pytest.approx(3.0 * 10.0)


def test_junction_strand_gated_discrimination(monkeypatch):
    # Junction flanks (jl=1 donor, jr=2 acceptor). t0 (+): donor 6 & acceptor 4 both strand-match → min 4.
    # t1 (−): both junction strands are +, mismatching − → 0 & 0 → min 0 → the isoform is CRUSHED (0, not
    # NaN: its flanks are observed). This is the exon-sharing-isoform discrimination.
    warm = _run(
        monkeypatch,
        rt=[], rr=[], bt=[], br=[], jt=[0, 1], jl=[1, 1], jr=[2, 2],
        strand=[STRAND_POS, STRAND_NEG], eff=[10.0, 10.0], n_t=2,
    )
    assert warm[0] == pytest.approx(4.0 * 10.0)
    assert warm[1] == pytest.approx(0.0)


def test_missing_warm_start_or_efflen_returns_all_nan(monkeypatch):
    monkeypatch.setattr(cel, "_transcript_node_incidence",
                        lambda index, region_arrays: (np.empty(0, np.int64),) * 7)
    index = SimpleNamespace(num_transcripts=3, t_to_strand_arr=np.array([STRAND_POS] * 3))
    ra = SimpleNamespace(ref_id=np.zeros(4, dtype=np.int64))
    cal = _calibration()
    # no effective_lengths_em → full fallback
    assert np.all(np.isnan(build_transcript_warm_start(cal, ra, index, None)))
    # no rna_warm_start → full fallback
    cal_no_ws = SimpleNamespace(**{**cal.__dict__, "rna_warm_start": None})
    assert np.all(np.isnan(build_transcript_warm_start(cal_no_ws, ra, index, np.ones(3))))
