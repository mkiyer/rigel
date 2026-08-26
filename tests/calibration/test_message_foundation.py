"""Gates for the message-propagation FOUNDATION SPEC (`calibration.messages.foundation`).

The spec is the owner's architecture (2026-08-26) made executable: one message type with the
unspliced and spliced lanes SEPARATE (provenance is structural), one propagation-time propagate
rule, one solve-time solve interface. Variance models implement the spec by overriding its
narrow extension points; the base classes enforce the rules no implementation may break.

Written BEFORE the module, verified failing; every gate watched firing against a deliberately
broken build.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.messages import PsiMessage


def _F():
    from rigel.calibration.messages import foundation as F

    return F


# ── the message type — provenance is structural ─────────────────────────────────────────────────


def test_the_message_carries_exactly_the_axiom_populations():
    """AXIOM 0 made structural: THREE unspliced populations (gDNA, RNA+, RNA−) and TWO spliced
    ones (certified RNA per strand — gDNA cannot splice, so a spliced gDNA lane is inexpressible).
    A fourth population cannot be written down."""
    F = _F()
    import dataclasses

    fields = {f.name for f in dataclasses.fields(F.Message)}
    assert fields == {
        "unspliced_gdna",
        "unspliced_rna_pos",
        "unspliced_rna_neg",
        "spliced_rna_pos",
        "spliced_rna_neg",
    }


def test_silence_is_expressible_and_detectable():
    F = _F()
    assert F.Message.silent().is_silent
    assert F.Claim.silent().is_silent
    loud = F.Message.silent().with_lane("spliced_rna_pos", F.Claim(abundance=2.0, precision=5.0))
    assert not loud.is_silent


# ── the propagation-time propagate rule ───────────────────────────────────────────────────────────


def test_pass_through_when_the_node_has_nothing():
    """THE RULE (owner): a node with no own belief relays the incoming message unchanged."""
    F = _F()
    incoming = F.Message(
        unspliced_gdna=F.Claim(3.0, 10.0),
        unspliced_rna_pos=F.Claim(1.5, 4.0),
        unspliced_rna_neg=F.Claim.silent(),
        spliced_rna_pos=F.Claim(0.8, 25.0),
        spliced_rna_neg=F.Claim.silent(),
    )
    out = F.PassThroughPropagation().propagate(F.Message.silent(), incoming)
    for lane in F.Message.LANES:
        a, b = out.lane(lane), incoming.lane(lane)
        assert a.abundance == b.abundance and a.precision == b.precision, lane


def test_a_strong_own_belief_integrates():
    """A node with its own belief blends: the relayed abundance lies between the two claims and
    the relayed precision is their (independent-evidence) sum."""
    F = _F()
    own = F.Message.silent().with_lane("unspliced_gdna", F.Claim(4.0, 20.0))
    incoming = F.Message.silent().with_lane("unspliced_gdna", F.Claim(1.0, 5.0))
    out = F.PassThroughPropagation().propagate(own, incoming)
    g = out.lane("unspliced_gdna")
    assert 1.0 < g.abundance < 4.0
    assert np.isclose(g.precision, 25.0)
    assert g.abundance > 2.5, "the stronger claim must dominate the blend"


def test_lanes_do_not_mix_at_propagation():
    """Provenance is load-bearing: spliced information may never leak into an unspliced lane (or
    vice versa) during propagation — the lanes travel separately so the solve can treat them
    by provenance."""
    F = _F()
    incoming = F.Message.silent().with_lane("spliced_rna_pos", F.Claim(9.0, 100.0))
    out = F.PassThroughPropagation().propagate(F.Message.silent(), incoming)
    assert out.lane("spliced_rna_pos").precision == 100.0
    for lane in F.Message.LANES:
        if lane != "spliced_rna_pos":
            assert out.lane(lane).is_silent, f"{lane} must stay silent"


def test_the_propagation_refuses_an_amplifying_attenuation():
    """The single-source rule, made structural: an attenuation (the propagation variance model's
    override point) may only LOWER a claim's precision. A model that amplifies must be refused by
    the BASE class — no subclass can smuggle an amplifier past the spec."""
    F = _F()

    class Amplifier(F.PassThroughPropagation):
        def attenuate(self, claim, lane, hop=None):
            if claim.precision > 0:
                return F.Claim(claim.abundance, claim.precision * 2.0)
            return claim

    incoming = F.Message.silent().with_lane("unspliced_gdna", F.Claim(1.0, 5.0))
    with pytest.raises(ValueError, match="amplif"):
        Amplifier().propagate(F.Message.silent(), incoming)


def test_a_lowering_attenuation_is_accepted_and_applied():
    F = _F()

    class Halver(F.PassThroughPropagation):
        def attenuate(self, claim, lane, hop=None):
            return F.Claim(claim.abundance, claim.precision * 0.5)

    incoming = F.Message.silent().with_lane("unspliced_rna_pos", F.Claim(2.0, 8.0))
    out = Halver().propagate(F.Message.silent(), incoming)
    assert np.isclose(out.lane("unspliced_rna_pos").precision, 4.0)


# ── the solve-time solve interface ────────────────────────────────────────────────────────────────


def test_the_solve_requires_both_lanes_to_be_implemented():
    """The solve model must say what it does with the unspliced lanes AND with the spliced lanes
    (the certified-flux treatment is a spliced-lane solve). A model implementing only one is not
    a solve model."""
    F = _F()

    class HalfDone(F.SolveModel):
        def solve_unspliced(self, own, forward, backward):
            return PsiMessage.silent()

    with pytest.raises(TypeError):
        HalfDone()


def test_the_solve_assembles_the_two_lanes_into_one_psi_message():
    """The template: the unspliced solve supplies the Gaussian channels, the spliced solve
    supplies the row factor, and the base assembles them — the seam the shipped tree already
    has (deliver's channels + the certified-flux rows), formalized."""
    F = _F()

    class Trivial(F.SolveModel):
        def solve_unspliced(self, own, forward, backward):
            return PsiMessage(lam_mode=np.zeros(3), lam_prec=np.ones(3))

        def solve_spliced(self, own, forward, backward):
            return np.ones((3, 5))

    msg = Trivial().solve(F.Message.silent(), F.Message.silent(), F.Message.silent())
    assert isinstance(msg, PsiMessage)
    assert msg.lam_prec is not None
    assert msg.lam_rows is not None and msg.lam_rows.shape == (3, 5)
