from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._arrays import RegionArrays
from rigel.calibration.adaptive_prior import (
    LOCUS_ESS_CEIL,
    MAX_ESS,
    PRIOR_BIAS_APPLIED,
    PRIOR_ESS_CAPPED,
    PRIOR_NO_UNSPLICED_MASS,
    PRIOR_STRUCTURAL_GATED,
    compute_adaptive_prior,
)
from rigel.locus import Locus, MultiLocus


GDNA_STATE = [1.0, 0.0]
RNA_STATE = [0.0, 1.0]
UNIFORM_STATE = [0.5, 0.5]
INTERIOR_STATE = [0.5, 0.5]


def _region_arrays(regions: list[tuple[int, int] | tuple[int, int, int]]) -> RegionArrays:
    ref_ids = np.array([row[0] if len(row) == 3 else 0 for row in regions], dtype=np.int32)
    starts = np.array([row[-2] for row in regions], dtype=np.int64)
    ends = np.array([row[-1] for row in regions], dtype=np.int64)
    order = np.lexsort((starts, ref_ids))
    ref_ids = ref_ids[order]
    starts = starts[order]
    ends = ends[order]
    n_refs = int(ref_ids.max(initial=0)) + 1 if ref_ids.size else 1
    counts = np.bincount(ref_ids, minlength=n_refs).astype(np.int32, copy=False)
    ref_offsets = np.empty(n_refs + 1, dtype=np.int32)
    ref_offsets[0] = 0
    np.cumsum(counts, out=ref_offsets[1:])
    return RegionArrays(
        ref_id=ref_ids,
        start=starts,
        end=ends,
        signature=np.zeros(ref_ids.size, dtype=np.uint8),
        ts_class=np.zeros(ref_ids.size, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64, copy=False),
        ref_offsets=ref_offsets,
        order=order.astype(np.int64, copy=False),
        n_refs=n_refs,
    )


def _multi_locus(locus_id: int, start: int, end: int, *, ref_id: int = 0) -> MultiLocus:
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([locus_id], dtype=np.int32),
        unit_indices=np.array([locus_id], dtype=np.int32),
        gdna_span=max(end - start, 1),
        loci=(Locus(ref=f"chr{ref_id + 1}", ref_id=ref_id, start=start, end=end),),
    )


def _compute(
    *,
    regions: list[tuple[int, int] | tuple[int, int, int]],
    p_states: list[list[float]],
    unspliced: list[float],
    loci: list[MultiLocus],
    has_gdna: list[bool] | None = None,
    gdna_unspliced: list[float] | None = None,
    rna_unspliced: list[float] | None = None,
    rna_call_bias: float = 0.5,
    max_ess: float = MAX_ESS,
    locus_ess: list[float] | None = None,
):
    if has_gdna is None:
        has_gdna = [True] * len(loci)
    if gdna_unspliced is None or rna_unspliced is None:
        gdna_unspliced, rna_unspliced = _legacy_test_mass_split(p_states, unspliced)
    if locus_ess is None:
        # Default to a value >= LOCUS_ESS_CEIL so the ESS shrink ramp is
        # inactive (shrink = 1.0) and existing assertions hold unchanged.
        locus_ess = [LOCUS_ESS_CEIL] * len(loci)
    return compute_adaptive_prior(
        region_arrays=_region_arrays(regions),
        multi_loci=loci,
        p_states=np.asarray(p_states, dtype=np.float64),
        total_mass=np.asarray(unspliced, dtype=np.float64),
        gdna_mass=np.asarray(gdna_unspliced, dtype=np.float64),
        rna_mass=np.asarray(rna_unspliced, dtype=np.float64),
        has_gdna_candidate=np.asarray(has_gdna, dtype=bool),
        locus_ess=np.asarray(locus_ess, dtype=np.float64),
        rna_call_bias=rna_call_bias,
        max_ess=max_ess,
    )


def _legacy_test_mass_split(
    p_states: list[list[float]], unspliced: list[float]
) -> tuple[list[float], list[float]]:
    states = np.asarray(p_states, dtype=np.float64)
    total = np.asarray(unspliced, dtype=np.float64)
    state_gdna_share = states[:, 0]
    state_rna_share = states[:, 1]
    state_total = state_gdna_share + state_rna_share
    gdna_share = np.divide(
        state_gdna_share,
        state_total,
        out=np.zeros_like(state_gdna_share),
        where=state_total > 0.0,
    )
    gdna = total * gdna_share
    rna = total - gdna
    return gdna.tolist(), rna.tolist()


def test_unexpressed_soft_gate_boundaries_and_monotonicity() -> None:
    result = _compute(
        regions=[(0, 10), (10, 20), (20, 30), (30, 40)],
        p_states=[GDNA_STATE, [0.8, 0.2], UNIFORM_STATE, RNA_STATE],
        unspliced=[1.0, 1.0, 1.0, 1.0],
        loci=[],
        has_gdna=[],
    )

    assert result.region_weight[0] == pytest.approx(1.0)
    assert result.region_weight[1] == pytest.approx(0.8)
    assert result.region_weight[2] == pytest.approx(0.5)
    assert result.region_weight[3] == pytest.approx(0.0)


def test_regional_pseudocount_boundaries() -> None:
    result = _compute(
        regions=[(0, 10), (10, 20), (20, 30), (30, 40)],
        p_states=[GDNA_STATE, UNIFORM_STATE, RNA_STATE, GDNA_STATE],
        unspliced=[0.0, 10.0, 20.0, 30.0],
        loci=[
            _multi_locus(0, 0, 10),
            _multi_locus(1, 10, 20),
            _multi_locus(2, 20, 30),
            _multi_locus(3, 30, 40),
        ],
    )

    np.testing.assert_allclose(
        result.n_local,
        np.array([[0.0, 0.0], [2.5, 2.5], [0.0, 0.0], [30.0, 0.0]]),
    )


def test_explicit_mass_split_overrides_latent_state_rna_gdna_split() -> None:
    result = _compute(
        regions=[(0, 100)],
        p_states=[GDNA_STATE],
        unspliced=[100.0],
        gdna_unspliced=[60.0],
        rna_unspliced=[40.0],
        loci=[_multi_locus(0, 0, 100)],
    )

    np.testing.assert_allclose(result.n_local, np.array([[60.0, 40.0]]))
    np.testing.assert_allclose([result.alpha_gdna_add[0], result.alpha_rna_add[0]], [60.0, 40.0])


def test_geometry_disjoint_regions_allocate_to_matching_loci() -> None:
    result = _compute(
        regions=[(0, 100), (100, 200)],
        p_states=[GDNA_STATE, GDNA_STATE],
        unspliced=[10.0, 20.0],
        gdna_unspliced=[10.0, 0.0],
        rna_unspliced=[0.0, 20.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200)],
    )

    np.testing.assert_allclose(result.n_local, np.array([[10.0, 0.0], [0.0, 20.0]]))
    np.testing.assert_allclose(result.locus_unspliced, np.array([10.0, 20.0]))
    np.testing.assert_array_equal(result.n_regions_touched, np.array([1, 1], dtype=np.int32))


def test_geometry_split_region_allocates_proportionally() -> None:
    result = _compute(
        regions=[(0, 100)],
        p_states=[GDNA_STATE],
        unspliced=[10.0],
        loci=[_multi_locus(0, 0, 60), _multi_locus(1, 40, 100)],
    )

    np.testing.assert_allclose(result.n_local[:, 0], np.array([5.0, 5.0]))
    np.testing.assert_allclose(result.multi_locus_region_mass, np.array([5.0, 5.0]))


def test_global_minus_local_is_componentwise_nonnegative() -> None:
    result = _compute(
        regions=[(0, 100), (100, 200), (200, 300)],
        p_states=[GDNA_STATE, GDNA_STATE, INTERIOR_STATE],
        unspliced=[10.0, 20.0, 40.0],
        gdna_unspliced=[10.0, 0.0, 20.0],
        rna_unspliced=[0.0, 20.0, 20.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)],
    )

    assert np.all(result.n_other >= 0.0)


def test_confident_locus_ignores_global_pool() -> None:
    result = _compute(
        regions=[(0, 100), (100, 200)],
        p_states=[GDNA_STATE, GDNA_STATE],
        unspliced=[10.0, 20.0],
        gdna_unspliced=[0.0, 20.0],
        rna_unspliced=[10.0, 0.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200)],
    )

    assert result.locus_weight[0] == pytest.approx(1.0)
    np.testing.assert_allclose(result.alpha_gdna_add[0], 0.0)
    np.testing.assert_allclose(result.alpha_rna_add[0], 10.0)


def test_ambiguous_locus_uses_local_soft_gate_without_global_fallback() -> None:
    result = _compute(
        regions=[(0, 100), (100, 200)],
        p_states=[UNIFORM_STATE, GDNA_STATE],
        unspliced=[10.0, 5.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200)],
    )

    assert result.locus_weight[0] == pytest.approx(0.5)
    assert result.shrink_weight[0] == pytest.approx(0.0)
    np.testing.assert_allclose(result.alpha_gdna_add[0], 2.5)
    np.testing.assert_allclose(result.alpha_rna_add[0], 2.5)


def test_cap_rescales_to_locus_unspliced_and_preserves_direction() -> None:
    result = _compute(
        regions=[(0, 100), (100, 200), (200, 300)],
        p_states=[UNIFORM_STATE, GDNA_STATE, GDNA_STATE],
        unspliced=[10_000.0, 100.0, 100.0],
        gdna_unspliced=[5_000.0, 100.0, 0.0],
        rna_unspliced=[5_000.0, 0.0, 100.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)],
    )

    assert result.ess_final[0] == pytest.approx(MAX_ESS)
    np.testing.assert_allclose([result.alpha_gdna_add[0], result.alpha_rna_add[0]], [1500.0, 1500.0])
    assert int(result.flags[0] & PRIOR_ESS_CAPPED) != 0


def test_structural_gate_zeroes_final_alphas() -> None:
    result = _compute(
        regions=[(0, 100)],
        p_states=[GDNA_STATE],
        unspliced=[10.0],
        loci=[_multi_locus(0, 0, 100)],
        has_gdna=[False],
    )

    assert result.alpha_gdna_add[0] == pytest.approx(0.0)
    assert result.alpha_rna_add[0] == pytest.approx(0.0)
    assert int(result.flags[0] & PRIOR_STRUCTURAL_GATED) != 0


def test_single_locus_sample_uses_only_local_evidence() -> None:
    result = _compute(
        regions=[(0, 100)],
        p_states=[GDNA_STATE],
        unspliced=[25.0],
        loci=[_multi_locus(0, 0, 100)],
    )

    np.testing.assert_allclose(result.n_other, np.array([[0.0, 0.0]]))
    assert result.alpha_gdna_add[0] == pytest.approx(25.0)
    assert result.alpha_rna_add[0] == pytest.approx(0.0)


def test_locus_order_shuffle_preserves_id_indexed_outputs() -> None:
    loci = [_multi_locus(0, 0, 100), _multi_locus(1, 100, 200)]
    forward = _compute(
        regions=[(0, 100), (100, 200)],
        p_states=[GDNA_STATE, GDNA_STATE],
        unspliced=[10.0, 20.0],
        gdna_unspliced=[10.0, 0.0],
        rna_unspliced=[0.0, 20.0],
        loci=loci,
    )
    shuffled = _compute(
        regions=[(0, 100), (100, 200)],
        p_states=[GDNA_STATE, GDNA_STATE],
        unspliced=[10.0, 20.0],
        gdna_unspliced=[10.0, 0.0],
        rna_unspliced=[0.0, 20.0],
        loci=list(reversed(loci)),
    )

    np.testing.assert_allclose(shuffled.alpha_gdna_add, forward.alpha_gdna_add)
    np.testing.assert_allclose(shuffled.alpha_rna_add, forward.alpha_rna_add)
    np.testing.assert_allclose(shuffled.n_local, forward.n_local)


def test_no_unspliced_mass_sets_flag_and_zeroes_alphas() -> None:
    result = _compute(
        regions=[(0, 100)],
        p_states=[GDNA_STATE],
        unspliced=[0.0],
        loci=[_multi_locus(0, 0, 100)],
    )

    assert result.ess_final[0] == pytest.approx(0.0)
    assert int(result.flags[0] & PRIOR_NO_UNSPLICED_MASS) != 0


def test_default_rna_call_bias_is_identical_to_explicit_half() -> None:
    kwargs = dict(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
    )
    default = _compute(**kwargs)
    explicit = _compute(**kwargs, rna_call_bias=0.5)

    np.testing.assert_array_equal(default.alpha_gdna_add, explicit.alpha_gdna_add)
    np.testing.assert_array_equal(default.alpha_rna_add, explicit.alpha_rna_add)


def test_rna_call_bias_near_zero_shifts_interior_share_toward_gdna() -> None:
    neutral = _compute(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
    )
    shifted = _compute(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
        rna_call_bias=1.0e-4,
    )

    assert shifted.alpha_rna_add[0] < 0.01
    assert shifted.alpha_gdna_add[0] == pytest.approx(neutral.ess_final[0], rel=1.0e-3)
    assert int(shifted.flags[0] & PRIOR_BIAS_APPLIED) != 0


def test_rna_call_bias_near_one_shifts_interior_share_toward_rna() -> None:
    neutral = _compute(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
    )
    shifted = _compute(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
        rna_call_bias=1.0 - 1.0e-4,
    )

    assert shifted.alpha_gdna_add[0] < 0.01
    assert shifted.alpha_rna_add[0] == pytest.approx(neutral.ess_final[0], rel=1.0e-3)
    assert int(shifted.flags[0] & PRIOR_BIAS_APPLIED) != 0


def test_rna_call_bias_preserves_mass_and_cap() -> None:
    neutral = _compute(
        regions=[(0, 100), (100, 200), (200, 300)],
        p_states=[UNIFORM_STATE, GDNA_STATE, RNA_STATE],
        unspliced=[10.0, 100.0, 100.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)],
    )
    shifted = _compute(
        regions=[(0, 100), (100, 200), (200, 300)],
        p_states=[UNIFORM_STATE, GDNA_STATE, RNA_STATE],
        unspliced=[10.0, 100.0, 100.0],
        loci=[_multi_locus(0, 0, 100), _multi_locus(1, 100, 200), _multi_locus(2, 200, 300)],
        rna_call_bias=0.8,
    )

    np.testing.assert_allclose(shifted.ess_final, neutral.ess_final)
    assert shifted.ess_final[0] <= min(shifted.locus_unspliced[0], MAX_ESS)


def test_structural_gate_wins_over_rna_call_bias() -> None:
    shifted = _compute(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
        has_gdna=[False],
        rna_call_bias=0.9,
    )

    assert shifted.alpha_gdna_add[0] == pytest.approx(0.0)
    assert shifted.alpha_rna_add[0] == pytest.approx(0.0)
    assert int(shifted.flags[0] & PRIOR_BIAS_APPLIED) == 0


def test_rna_call_bias_monotonically_increases_final_rna_share() -> None:
    kwargs = dict(
        regions=[(0, 100)],
        p_states=[INTERIOR_STATE],
        unspliced=[100.0],
        loci=[_multi_locus(0, 0, 100)],
    )
    low = _compute(**kwargs, rna_call_bias=0.2)
    neutral = _compute(**kwargs, rna_call_bias=0.5)
    high = _compute(**kwargs, rna_call_bias=0.8)

    assert low.alpha_rna_add[0] < neutral.alpha_rna_add[0] < high.alpha_rna_add[0]


def test_invalid_rna_call_bias_is_rejected() -> None:
    with pytest.raises(ValueError, match="rna_call_bias"):
        _compute(
            regions=[(0, 100)],
            p_states=[INTERIOR_STATE],
            unspliced=[100.0],
            loci=[_multi_locus(0, 0, 100)],
            rna_call_bias=1.0,
        )


# ---------------------------------------------------------------------------
# PR 03 Test 17 — locus_ess shrink ramp attenuates prior in low-evidence loci.
# ---------------------------------------------------------------------------


def test_17_locus_ess_shrink_attenuates_low_evidence_locus() -> None:
    """Two loci with identical mass but different locus_ess get different priors.

    The shrink ramp ``(ess - floor) / (ceil - floor)`` clipped to ``[0, 1]``
    multiplies the additive concentration. A locus at ``ess == floor`` collapses
    to zero prior; a locus at ``ess == ceil`` keeps full strength.
    """
    regions = [(0, 100), (100, 200)]
    p_states = [INTERIOR_STATE, INTERIOR_STATE]
    unspliced = [80.0, 80.0]
    loci = [_multi_locus(0, 0, 100), _multi_locus(1, 100, 200)]

    base = _compute(
        regions=regions, p_states=p_states, unspliced=unspliced, loci=loci
    )
    # Baseline: identical priors on both loci.
    assert base.alpha_gdna_add[0] == pytest.approx(base.alpha_gdna_add[1])
    assert base.alpha_rna_add[0] == pytest.approx(base.alpha_rna_add[1])

    # Now supply asymmetric locus_ess: locus 0 has full evidence (ess=ceil),
    # locus 1 has minimum evidence (ess=floor).
    region_arrays = _region_arrays(regions)
    shrunk = compute_adaptive_prior(
        region_arrays=region_arrays,
        multi_loci=loci,
        p_states=np.asarray(p_states, dtype=np.float64),
        total_mass=np.asarray(unspliced, dtype=np.float64),
        gdna_mass=np.asarray(
            _legacy_test_mass_split(p_states, unspliced)[0], dtype=np.float64
        ),
        rna_mass=np.asarray(
            _legacy_test_mass_split(p_states, unspliced)[1], dtype=np.float64
        ),
        has_gdna_candidate=np.array([True, True], dtype=bool),
        locus_ess=np.array([50.0, 1.0], dtype=np.float64),
        locus_ess_floor=1.0,
        locus_ess_ceil=50.0,
    )

    # Locus 0 (ess=ceil) keeps full strength.
    assert shrunk.alpha_gdna_add[0] == pytest.approx(base.alpha_gdna_add[0])
    assert shrunk.alpha_rna_add[0] == pytest.approx(base.alpha_rna_add[0])
    # Locus 1 (ess=floor) collapses to zero.
    assert shrunk.alpha_gdna_add[1] == pytest.approx(0.0)
    assert shrunk.alpha_rna_add[1] == pytest.approx(0.0)
    # ess_final tracks the post-shrink concentration.
    assert shrunk.ess_final[1] == pytest.approx(0.0)


def test_17b_locus_ess_validation() -> None:
    """``locus_ess`` shape/value validation."""
    regions = [(0, 100)]
    p_states = [INTERIOR_STATE]
    unspliced = [10.0]
    loci = [_multi_locus(0, 0, 100)]
    gdna, rna = _legacy_test_mass_split(p_states, unspliced)
    region_arrays = _region_arrays(regions)

    with pytest.raises(ValueError, match="locus_ess must have shape"):
        compute_adaptive_prior(
            region_arrays=region_arrays,
            multi_loci=loci,
            p_states=np.asarray(p_states, dtype=np.float64),
            total_mass=np.asarray(unspliced, dtype=np.float64),
            gdna_mass=np.asarray(gdna, dtype=np.float64),
            rna_mass=np.asarray(rna, dtype=np.float64),
            has_gdna_candidate=np.array([True], dtype=bool),
            locus_ess=np.array([1.0, 2.0], dtype=np.float64),
        )

    with pytest.raises(ValueError, match="locus_ess must be finite"):
        compute_adaptive_prior(
            region_arrays=region_arrays,
            multi_loci=loci,
            p_states=np.asarray(p_states, dtype=np.float64),
            total_mass=np.asarray(unspliced, dtype=np.float64),
            gdna_mass=np.asarray(gdna, dtype=np.float64),
            rna_mass=np.asarray(rna, dtype=np.float64),
            has_gdna_candidate=np.array([True], dtype=bool),
            locus_ess=np.array([-1.0], dtype=np.float64),
        )

    with pytest.raises(ValueError, match="locus_ess_floor"):
        compute_adaptive_prior(
            region_arrays=region_arrays,
            multi_loci=loci,
            p_states=np.asarray(p_states, dtype=np.float64),
            total_mass=np.asarray(unspliced, dtype=np.float64),
            gdna_mass=np.asarray(gdna, dtype=np.float64),
            rna_mass=np.asarray(rna, dtype=np.float64),
            has_gdna_candidate=np.array([True], dtype=bool),
            locus_ess=np.array([1.0], dtype=np.float64),
            locus_ess_floor=10.0,
            locus_ess_ceil=5.0,
        )
