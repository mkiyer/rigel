"""The fitted strand protocol governs local RNA witnesses, independent of remote expression."""

from dataclasses import fields, replace

import numpy as np
import pytest

from rigel.calibration.messages import ChainView
from rigel.calibration.messages.transfer import TransferPolicy
from rigel.calibration.splice_graph import FLAG_TES_NEG
from _transfer_harness import BlockContext, _leaves, _native_passes, _prepared


def _local_context(both_strands=False):
    """An intron/boundary pair; the right region can admit both RNA strands.

    The both-strand topology reproduces an annotation-derived intron next to
    a negative-strand exon edge. All junction fluxes are zero in these examples.
    """
    pair = np.zeros((3, 2))
    return BlockContext(
        eff_gdna=np.array([2980.0, 221.0, 379.0]) if both_strands else np.full(3, 100.0),
        eff_rna=np.array([2972.0, 228.0, 372.0]) if both_strands else np.full(3, 100.0),
        unspliced_count=(
            np.array([[2.0, 3.0], [0.0, 1.0], [0.0, 1.0]])
            if both_strands
            else np.array([[90.0, 10.0], [60.0, 40.0], [90.0, 10.0]])
        ),
        spliced_count=pair.copy(),
        sj_count=pair.copy(),
        sj_count_lo=pair.copy(),
        sj_count_hi=pair.copy(),
        route_rate_lo=pair.copy(),
        route_rate_hi=pair.copy(),
        left=np.array([-1, 0, 1]),
        right=np.array([1, 2, -1]),
        is_boundary=np.array([False, True, False]),
        is_exon_region=np.array([False, False, both_strands]),
        free_pos=np.ones(3, bool),
        free_neg=np.array([False, False, both_strands]),
        exon_pos=np.zeros(3, bool),
        exon_neg=np.array([False, False, both_strands]),
        boundary_flags=np.array([0, FLAG_TES_NEG if both_strands else 0, 0], np.uint16),
        strand_live=True,
        has_own_composition=np.array([True, True, not both_strands]),
        n_grid=101,
        logodds_window=10.0,
    )


def _library_view(ctx, remote_count):
    """Three disconnected objects hold the two coordinates fixed.

    Pure DNA and a both-strand exon set positive origins. The last exon
    has zero RNA opportunity, so changing its count cannot move either origin.
    """
    data = {
        f.name: np.asarray(getattr(ctx, f.name)).copy()
        for f in fields(ChainView)
        if f.name != "strand_live"
    }
    tail = {
        name: np.zeros((3, *value.shape[1:]), dtype=value.dtype) for name, value in data.items()
    }
    tail["eff_gdna"][:] = 100.0
    tail["eff_rna"][1] = 100.0
    tail["unspliced_count"][:] = [[50.0, 50.0], [50.0, 50.0], [remote_count, 0.0]]
    tail["is_exon_region"][1:] = True
    tail["free_pos"][1:] = True
    tail["free_neg"][1] = True
    tail["exon_pos"][1:] = True
    tail["exon_neg"][1] = True
    tail["left"][:] = tail["right"][:] = -1
    return ChainView(
        **{name: np.concatenate([value, tail[name]]) for name, value in data.items()},
        strand_live=ctx.strand_live,
    )


@pytest.mark.parametrize("kappa", [0.95, 0.05])
def test_counted_introns_do_not_need_a_counted_exon_to_use_the_strand_protocol(kappa):
    view = _library_view(_local_context(), 0.0)
    assert TransferPolicy(kappa).library(view).split_live


@pytest.mark.parametrize("kappa, live", [(None, True), (None, False), (0.5, False), (0.95, False)])
def test_unavailable_strand_evidence_cannot_be_enabled_by_remote_expression(kappa, live):
    view = replace(_library_view(_local_context(), 100.0), strand_live=live)
    assert not TransferPolicy(kappa).library(view).split_live


@pytest.mark.parametrize("kappa", [0.95, 0.05])
@pytest.mark.parametrize("both_strands", [False, True])
def test_disconnected_expression_preserves_local_sources_messages_and_delivery(kappa, both_strands):
    ctx, policy = _local_context(both_strands), TransferPolicy(kappa)
    results = []
    for remote_count in (0.0, 1.0, 10000.0):
        library = policy.library(_library_view(ctx, remote_count))
        assert library.rho_gdna == 1.0
        assert library.rho_rna == (101.0 / 472.0 if both_strands else 1.0)
        prepared = _prepared(policy, ctx, library)
        received = _native_passes(prepared, ctx)
        results.append((prepared, received, prepared.solve(*received)))
    baseline, base_received, base_delivered = results[0]
    assert baseline.own.mask.any()
    assert any(part.level_rna_pos.present.any() for part in base_received)
    assert base_delivered.lam_rows is not None
    assert np.ptp(base_delivered.lam_rows) > 0.0
    if both_strands:
        assert base_delivered.cube_rows is not None
        assert base_delivered.cube_rows.has_pos.any()
    for prepared, received, delivered in results[1:]:
        np.testing.assert_array_equal(prepared.own.mask, baseline.own.mask)
        np.testing.assert_array_equal(
            prepared.own.rows[prepared.own.mask], baseline.own.rows[baseline.own.mask]
        )
        for name in prepared.lanes:
            a, b = baseline.lanes[name].own_level, prepared.lanes[name].own_level
            np.testing.assert_array_equal(a.mask, b.mask)
            np.testing.assert_array_equal(a.rows[a.mask], b.rows[b.mask])
        for a, b in zip(base_received, received):
            for name, x, y in _leaves(a, b):
                np.testing.assert_array_equal(x, y, err_msg=name)
        np.testing.assert_array_equal(base_delivered.lam_rows, delivered.lam_rows)
        if both_strands:
            a, b = base_delivered.cube_rows, delivered.cube_rows
            np.testing.assert_array_equal(a.slot, b.slot)
            for lane in ("pos", "neg"):
                mask = getattr(a, "has_" + lane)
                np.testing.assert_array_equal(mask, getattr(b, "has_" + lane))
                np.testing.assert_array_equal(
                    getattr(a, "profile_" + lane)[mask], getattr(b, "profile_" + lane)[mask]
                )
