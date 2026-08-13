"""The policy that sends nothing — and it is **THE DEFAULT**.

       Gate: ``tests/calibration/test_sweep_backbone.py``

⭐⭐⭐ **This file is five boundaries of behaviour, and that is the entire point.** A new session reads
``sweep.py`` plus this and holds the whole working system in their head. :mod:`~.head` is opt-in and is
clearly labelled as the legacy arm being dismantled.

⛔ **It is also a MEASURED floor, not a placeholder.** With no belief propagation whatsoever the panel's
deliverable is worth −50 % overall, and the sign is not uniform: it is a net **improvement** on three of
the four strata (stranded × capture-ON **−58.3 %, 16/16**; stranded × capture-OFF **−43.7 %, 16/16**;
unstranded × capture-OFF **−32.1 %**) and a large regression on exactly one (unstranded × capture-ON
**+154.8 %, 0/16**) — the stratum where ``kappa = 1/2`` makes the strand λ-term exactly 0, so a slot has no
own composition evidence at all and a message is the only source there is.

⚠ **So this is a diagnostic, not a proposal.** Shipping it takes the zero-gDNA control to +8,003 %. What it
establishes is that the message layer's whole value is ONE stratum, which is what scoped the backbone.
"""

from __future__ import annotations

from . import NeighbourState, PsiMessage, StepContext

__all__ = ["SilentPolicy"]


class _SilentRelay:
    def scan(self, *, backward: bool):
        return None  # nothing to relay, so the backbone runs no loop at all

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        return PsiMessage.silent()


class SilentPolicy:
    """Sends nothing on any channel. ψ then carries the slot's OWN evidence alone — its two strand counts,
    its spliced count, the derived reference, the fitted gDNA prior and the intron factory."""

    name = "silent"

    def prepare(self, ctx: StepContext) -> _SilentRelay:
        return _SilentRelay()
