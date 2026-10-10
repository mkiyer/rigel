"""Observed strand messages depend on counts, annotation and protocol, not beliefs."""

from dataclasses import replace

import numpy as np
import pytest
import _transfer_harness as harness
from _transfer_harness import BlockContext, _prepared, _rule
from scipy.special import expit
from scipy.stats import binom

from rigel.calibration.messages.transfer import TransferPolicy, _Library
from rigel.calibration.blocks import SweepCapture
from rigel.calibration.sweep import solve_chain


def _context():
    """An exon, its boundary and an intron; no archived scan or fitted prior.

    Contained opportunities use region lengths 500/900 and DNA/RNA fragment
    lengths 200/100. The boundary has 199/99 crossing starts, respectively.
    There are no observed junction fragments in this small source fixture.
    """
    return BlockContext(
        eff_gdna=np.array([301.0, 199.0, 701.0]),
        eff_rna=np.array([401.0, 99.0, 801.0]),
        sj_count=np.zeros((3, 2)),
        sj_count_lo=np.zeros((3, 2)),
        sj_count_hi=np.zeros((3, 2)),
        route_rate_lo=np.zeros((3, 2)),
        route_rate_hi=np.zeros((3, 2)),
        unspliced_count=np.array([[9.0, 2.0], [7.0, 3.0], [5.0, 4.0]]),
        spliced_count=np.zeros((3, 2)),
        left=np.array([-1, 0, 1]),
        right=np.array([1, 2, -1]),
        is_boundary=np.array([False, True, False]),
        is_exon_region=np.array([True, False, False]),
        free_pos=np.ones(3, dtype=bool),
        free_neg=np.zeros(3, dtype=bool),
        exon_pos=np.array([True, False, False]),
        exon_neg=np.zeros(3, dtype=bool),
        boundary_flags=np.zeros(3, dtype=np.uint16),
        strand_live=True,
        has_own_composition=np.ones(3, dtype=bool),
        n_grid=65,
        logodds_window=8.0,
    )


def _prepare(context, kappa=0.99):
    return _prepared(TransferPolicy(kappa), context, _Library(0.1, 0.2, True))


@pytest.mark.parametrize("counts", [(9, 0), (5, 4), (0, 9)])
def test_intron_observations_supply_the_same_claim_as_exon_observations(counts):
    context = _context()
    context.unspliced_count[2] = counts
    prepared = _prepare(context)
    assert prepared.own.mask[2], "A factory-free intron lost its own strand observations"
    exon = context.is_exon_region.copy()
    exon[2] = True
    reference = _prepare(replace(context, is_exon_region=exon))
    np.testing.assert_array_equal(prepared.own[2], reference.own[2])
    delivered = _rule(prepared, 2, 1)
    assert delivered is not None and np.ptp(delivered) > 0


@pytest.mark.parametrize("kappa", [0.5, 0.7, 0.99])
@pytest.mark.parametrize("slot", [0, 1, 2])
@pytest.mark.parametrize("counts", [(0, 1), (1, 0), (5, 4)])
def test_claim_equals_the_conditional_two_poisson_observation_law(slot, counts, kappa):
    context = _context()
    context.unspliced_count[slot] = counts
    prepared = _prepare(context, kappa)
    assert prepared.own.mask[slot]
    fraction = expit(prepared.faces.lam)
    probability = 0.5 * fraction + kappa * (1 - fraction)
    reference = binom.logpmf(counts[0], sum(counts), probability)
    reference -= reference.max()
    np.testing.assert_allclose(prepared.own[slot], reference, atol=2e-13, rtol=2e-13)


@pytest.mark.parametrize("name", ["belief_fg", "od_g", "od_r"])
def test_message_preparation_rejects_posterior_inputs(monkeypatch, name):
    prepare = harness.transfer_prepare

    def with_posterior_input(**kwargs):
        kwargs[name] = np.full(3, 0.5) if name == "belief_fg" else 0.0
        return prepare(**kwargs)

    monkeypatch.setattr(harness, "transfer_prepare", with_posterior_input)
    with pytest.raises(TypeError, match="incompatible function arguments"):
        _prepare(_context())


def test_strand_reversal_preserves_the_claim():
    context = _context()
    original = _prepare(context)
    assert original.own.mask.all()
    reverse = _prepare(
        replace(
            context,
            free_pos=context.free_neg,
            free_neg=context.free_pos,
            exon_pos=context.exon_neg,
            exon_neg=context.exon_pos,
            unspliced_count=context.unspliced_count[:, ::-1].copy(),
        )
    )
    np.testing.assert_array_equal(original.own.mask, reverse.own.mask)
    np.testing.assert_allclose(
        original.own.rows[original.own.mask],
        reverse.own.rows[reverse.own.mask],
        rtol=2e-13,
        atol=2e-11,
    )


@pytest.mark.parametrize("missing", ["protocol", "single_strand", "own_evidence"])
def test_missing_strand_information_does_not_create_a_claim(missing):
    context = _context()
    kappa = 0.99
    if missing == "protocol":
        kappa = None
    elif missing == "single_strand":
        context = replace(context, free_neg=np.ones(3, dtype=bool))
    else:
        context = replace(context, has_own_composition=np.zeros(3, dtype=bool))
    assert not _prepare(context, kappa).own.mask.any()


def test_delivered_evidence_is_invariant_to_incoming_beliefs(sweep_inputs):
    """Changing a valid count belief cannot change the observations sent to neighbours.

    Exercise the complete native sweep, including its own-evidence presence masks,
    at interior and exact-vertex beliefs. The count posterior is deliberately not
    compared: its frozen strand variance can depend on the incoming belief.
    """
    chain, statics, geometry, belief, regions = sweep_inputs["args"]
    free_pos, free_neg = statics.free_pos, statics.free_neg
    n_rna = free_pos.astype(int) + free_neg.astype(int)
    policy, _, _ = harness._full_policy(sweep_inputs)

    def delivery(fraction):
        fg = np.where(n_rna > 0, fraction, 1.0)
        share = (1 - fg) / np.maximum(n_rna, 1)
        seeded = replace(belief, f_g=fg, f_pos=share * free_pos, f_neg=share * free_neg)
        capture = SweepCapture()
        solve_chain(
            chain,
            statics,
            geometry,
            seeded,
            regions,
            **sweep_inputs["kw"],
            policy=policy,
            _capture=capture,
        )
        rows = {"lam_rows": capture.lam_rows}
        for side in ("from_left", "from_right"):
            received = getattr(capture, side)
            for name in ("has_neighbour", "has_composition"):
                rows[f"{side}/{name}"] = received[name]
            rows[f"{side}/composition"] = received["composition"][received["has_composition"]]
            for lane in harness.LANES:
                values = received[lane]
                rows[f"{side}/{lane}/present"] = values["present"]
                for name, value in values.items():
                    if name != "present":
                        rows[f"{side}/{lane}/{name}"] = value[values["present"]]
        cube = capture.cube_rows
        if cube is not None:
            for name in ("slot", "has_pos", "has_neg", "u", "total", "opportunity", "rho_ref"):
                rows[f"cube/{name}"] = getattr(cube, name)
            rows["cube/pos"] = cube.profile_pos[cube.has_pos]
            rows["cube/neg"] = cube.profile_neg[cube.has_neg]
        assert capture.lam_rows is not None and np.ptp(capture.lam_rows, axis=1).max() > 0
        return rows

    reference = delivery(0.5)
    for fraction in (0.0, 0.01, 0.99, 1.0):
        changed = delivery(fraction)
        assert changed.keys() == reference.keys()
        for name in reference:
            np.testing.assert_array_equal(
                changed[name], reference[name], err_msg=f"{fraction}: {name}"
            )
