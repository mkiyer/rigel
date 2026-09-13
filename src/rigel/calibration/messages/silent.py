"""The policy that sends nothing — the OFF state (``message_policy = "silent"``) and the measured floor.

       Gate: ``tests/calibration/test_sweep_backbone.py``

Five boundaries of behaviour, and that is the point: a reader of ``sweep.py`` plus this file holds
the whole working backbone. :mod:`~.transfer` is the shipped policy; this one is what every policy
is priced against.

It is a MEASURED floor, not a placeholder. On strand-specific data a sighted exon's own solve is
already excellent and a message can mostly only disturb it, so the bar a policy is held to is: WIN
on unstranded data, do minimal HARM against this floor on stranded data, the two halves never
pooled (`scripts/design/policy_benchmark.py` prints them apart).
"""

from __future__ import annotations

from . import ChainView, PsiMessage, BlockContext

__all__ = ["SilentPolicy"]


class _PreparedSilence:
    def propagate(self, received, *, backward: bool):
        return None  # sends nothing: every node with a neighbour holds silence from this side

    def solve(self, from_left, from_right) -> PsiMessage:
        return PsiMessage.silent()


class SilentPolicy:
    """Sends nothing on any channel. ψ then carries the slot's OWN evidence alone — its two strand counts,
    its spliced count, the fitted gDNA prior and the intron factory."""

    name = "silent"

    def library(self, view: ChainView) -> None:
        return None  # nothing to reduce: no message needs a library-wide fact

    def prepare(self, ctx: BlockContext, library) -> _PreparedSilence:
        return _PreparedSilence()
