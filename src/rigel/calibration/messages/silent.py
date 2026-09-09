"""The policy that sends nothing — the OFF state (``message_propagation = False``) and the measured floor.

       Gate: ``tests/calibration/test_sweep_backbone.py``

⭐⭐⭐ **This file is five boundaries of behaviour, and that is the entire point.** A new session reads
``sweep.py`` plus this and holds the whole working system in their head. :mod:`~.transfer` is the SHIPPED
policy (the default since 2026-09-09); this one is what every policy is priced against.

⛔ **It is a MEASURED floor, not a placeholder.** On strand-specific data a sighted exon's own solve is
excellent and a message can mostly only disturb it, so the bar a policy is held to is: WIN on unstranded
data, do minimal HARM against this floor on stranded data, the two halves never pooled
(`scripts/design/policy_benchmark.py` prints them apart). The shipped policy beats this floor on 7 of 8
ladder rows of each half, the other two within 1 % (2026-09-09).
"""

from __future__ import annotations

from . import PsiMessage, StepContext

__all__ = ["SilentPolicy"]


class _PreparedSilence:
    def propagate(self, *, backward: bool):
        return None  # sends nothing: every node holds SILENCE from this side

    def solve(self, from_left, from_right) -> PsiMessage:
        return PsiMessage.silent()


class SilentPolicy:
    """Sends nothing on any channel. ψ then carries the slot's OWN evidence alone — its two strand counts,
    its spliced count, the fitted gDNA prior and the intron factory."""

    name = "silent"

    def prepare(self, ctx: StepContext) -> _PreparedSilence:
        return _PreparedSilence()
