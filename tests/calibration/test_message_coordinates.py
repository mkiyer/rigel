"""A missing numerical coordinate must not discard observed RNA evidence."""

from dataclasses import replace

import numpy as np
import pytest

from _transfer_harness import _hop, _prepared
from test_message_locality import _library_view, _local_context
from test_transfer_rna_lanes import _empty_piece_ctx
from rigel.calibration.messages.transfer import TransferPolicy
from rigel.calibration.splice_graph import FLAG_ACCEPTOR_NEG, FLAG_DONOR_NEG, FLAG_DONOR_POS


def _flux_only(strand="pos", count=40.0, rate=0.02, side="hi"):
    ctx = _empty_piece_ctx(flux=count, rate=rate)
    ctx = replace(
        ctx,
        unspliced_count=np.zeros_like(ctx.unspliced_count),
        has_own_composition=np.zeros(ctx.n_slots, bool),
    )
    if side == "lo":
        counts = np.zeros_like(ctx.sj_count)
        rates = np.zeros_like(ctx.route_rate_hi)
        counts[3, 0], rates[3, 0] = count, rate
        flags = np.zeros_like(ctx.boundary_flags)
        flags[3] = FLAG_DONOR_POS
        ctx = replace(
            ctx,
            sj_count=counts.copy(),
            sj_count_lo=counts,
            sj_count_hi=np.zeros_like(counts),
            route_rate_lo=rates,
            route_rate_hi=np.zeros_like(rates),
            boundary_flags=flags,
        )
    if strand == "neg":
        flags = ctx.boundary_flags.copy()
        flags[1 if side == "hi" else 3] = FLAG_ACCEPTOR_NEG if side == "hi" else FLAG_DONOR_NEG
        ctx = replace(
            ctx,
            free_pos=ctx.free_neg,
            free_neg=ctx.free_pos,
            exon_pos=ctx.exon_neg,
            exon_neg=ctx.exon_pos,
            sj_count=ctx.sj_count[:, ::-1].copy(),
            sj_count_lo=ctx.sj_count_lo[:, ::-1].copy(),
            sj_count_hi=ctx.sj_count_hi[:, ::-1].copy(),
            route_rate_lo=ctx.route_rate_lo[:, ::-1].copy(),
            route_rate_hi=ctx.route_rate_hi[:, ::-1].copy(),
            boundary_flags=flags,
        )
    return ctx


def _scaled(ctx, scale):
    return replace(
        ctx,
        eff_rna=ctx.eff_rna * scale,
        route_rate_lo=ctx.route_rate_lo / scale,
        route_rate_hi=ctx.route_rate_hi / scale,
    )


@pytest.mark.parametrize("empty_single", [False, True])
def test_local_sources_survive_when_the_exon_reduction_is_zero(empty_single):
    ctx = _local_context()
    view = _library_view(ctx, 0.0)
    if empty_single:
        opportunity = view.eff_rna.copy()
        opportunity[-1] = 100.0
        view = replace(view, eff_rna=opportunity)
    else:
        counts = view.unspliced_count.copy()
        counts[-2] = 0.0
        view = replace(view, unspliced_count=counts)
    policy = TransferPolicy(0.95)
    facts = policy.library(view)
    assert facts.rho_rna > 0.0
    prepared = _prepared(policy, ctx, facts)
    assert prepared.lanes["pos"].own_level.mask.all()
    assert _hop(prepared, 0, 1).level_rna_pos.present[1]


@pytest.mark.parametrize("strand", ["pos", "neg"])
@pytest.mark.parametrize("side", ["lo", "hi"])
@pytest.mark.parametrize("kappa", [0.99, 0.5, None])
def test_flux_survives_when_every_unspliced_object_is_empty(strand, side, kappa):
    ctx, policy = _flux_only(strand, side=side), TransferPolicy(kappa)
    ctx = replace(ctx, strand_live=kappa is not None and kappa != 0.5)
    facts = policy.library(ctx)
    assert facts.rho_rna > 0.0
    prepared = _prepared(policy, ctx, facts)
    lane = prepared.lanes[strand]
    assert lane.empty[2] and lane.own_level[2] is not None
    target = 3 if side == "hi" else 1
    received = getattr(_hop(prepared, 2, target), "level_rna_" + strand)
    assert received.present[target]
    assert received.count[target] == 40.0
    assert received.opportunity[target] == 40.0 / 0.02
    assert not received.has_witness[target]


@pytest.mark.parametrize("kind", ["observed", "flux"])
def test_rna_unit_change_preserves_the_source_profiles(kind):
    ctx = _local_context() if kind == "observed" else _flux_only()
    other, policy = _scaled(ctx, 7.0), TransferPolicy(0.95)
    a, b = policy.library(ctx), policy.library(other)
    assert a.rho_rna > 0.0
    np.testing.assert_allclose(b.rho_rna * 7.0, a.rho_rna, rtol=1e-14)
    pa, pb = _prepared(policy, ctx, a), _prepared(policy, other, b)
    for strand in ("pos", "neg"):
        mask = pa.lanes[strand].own_level.mask
        np.testing.assert_array_equal(mask, pb.lanes[strand].own_level.mask)
        np.testing.assert_allclose(
            pa.lanes[strand].own_level.rows[mask],
            pb.lanes[strand].own_level.rows[mask],
            atol=1e-10,
            rtol=1e-12,
        )


@pytest.mark.parametrize("kind", ["observed", "flux"])
def test_dna_opportunity_does_not_set_the_rna_coordinate(kind):
    ctx = _local_context() if kind == "observed" else _flux_only()
    policy = TransferPolicy(0.95)
    a, b = policy.library(ctx), policy.library(replace(ctx, eff_gdna=ctx.eff_gdna * 11.0))
    assert a.rho_rna > 0.0 and a.rho_rna == b.rho_rna


@pytest.mark.parametrize("sense", [0.95, 0.05])
def test_an_empty_protocol_aligned_column_does_not_remove_the_coordinate(sense):
    ctx = _local_context()
    counts = np.zeros_like(ctx.unspliced_count)
    counts[:, 1 if sense > 0.5 else 0] = 5.0
    ctx = replace(ctx, unspliced_count=counts)
    policy = TransferPolicy(sense)
    a = policy.library(ctx)
    reflected = replace(
        ctx, unspliced_count=counts[:, ::-1].copy(), free_pos=ctx.free_neg, free_neg=ctx.free_pos
    )
    b = policy.library(reflected)
    assert a.rho_rna > 0.0 and a.rho_rna == b.rho_rna


def test_existing_positive_coordinate_and_other_library_facts_are_preserved():
    ctx, policy = _local_context(), TransferPolicy(0.95)
    view = _library_view(ctx, 2.0)
    counts = view.unspliced_count.copy()
    counts[-2] = [60.0, 63.0]
    view = replace(view, unspliced_count=counts)
    a = policy.library(view)
    assert a.rho_rna == 1.23 and a.rho_gdna == 1.0 and a.split_live
    b = policy.library(ctx)
    assert b.rho_gdna == a.rho_gdna and b.split_live == a.split_live


@pytest.mark.parametrize("count,rate", [(0.0, 0.0), (40.0, 0.0), (0.0, 0.02)])
def test_absent_or_invalid_flux_cannot_create_a_coordinate_or_source(count, rate):
    ctx, policy = _flux_only(count=count, rate=rate), TransferPolicy(0.95)
    facts = policy.library(ctx)
    assert facts.rho_rna == 0.0
    prepared = _prepared(policy, ctx, facts)
    assert not any(prepared.lanes[s].own_level.mask.any() for s in ("pos", "neg"))


def test_intergenic_counts_alone_cannot_create_an_rna_source():
    ctx = _local_context()
    ctx = replace(
        ctx, free_pos=np.zeros(ctx.n_slots, bool), has_own_composition=np.zeros(ctx.n_slots, bool)
    )
    policy = TransferPolicy(0.95)
    facts = policy.library(ctx)
    assert facts.rho_rna == 0.0
    prepared = _prepared(policy, ctx, facts)
    assert not any(prepared.lanes[s].own_level.mask.any() for s in ("pos", "neg"))
