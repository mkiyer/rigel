"""The policy that sends nothing — and it is **THE DEFAULT**.

       Gate: ``tests/calibration/test_sweep_backbone.py``

⭐⭐⭐ **This file is five boundaries of behaviour, and that is the entire point.** A new session reads
``sweep.py`` plus this and holds the whole working system in their head. :mod:`~.head` is opt-in — one
config flag, ``CalibrationConfig.message_propagation``, with every operator inside it behind its own named
switch. ⛔ This line called ``head`` *"clearly labelled as the legacy arm being dismantled"* until
2026-08-17: it carries no such label, and it is not being dismantled. Its own docstring says the opposite —
the switches exist so ``ladder_arm_ab.py`` can price the operators ONE AT A TIME.

⛔ **It is also a MEASURED floor, not a placeholder.** With no belief propagation whatsoever the panel's
deliverable is worth −50 % overall, and the sign is not uniform: it is a net **improvement** on three of
the four strata (stranded × capture-ON **−58.3 %, 16/16**; stranded × capture-OFF **−43.7 %, 16/16**;
unstranded × capture-OFF **−32.1 %**) and a large regression on exactly one (unstranded × capture-ON
**+154.8 %, 0/16**) — the stratum where ``kappa = 1/2`` makes the strand λ-term exactly 0, so a slot has no
own composition evidence at all and a message is the only source there is.

⛔⛔ **AND IT IS SHIPPED, SO THAT PRICE IS BEING PAID — this paragraph read "so this is a diagnostic, not a
proposal" until 2026-08-17, which contradicted this file's own first line.** ``message_propagation``
defaults ``False`` (owner, 2026-08-07), so ``+8,003 %`` on the zero-gDNA control is the standing cost of
the shipped configuration, not a hypothetical one. ⚠ It stays down until the tool is optimised end to end
across all scenarios (owner, 2026-08-10) — a STUDY configuration reversed by one flag, not a retirement.
What the floor establishes is that the message layer's whole value is ONE stratum, which is what scoped
the backbone.

⛔ **NEEDS AN OWNER RULING — the ``−50 % overall`` above and ``CalibrationConfig.message_propagation``'s
``panel TOTAL is +99.9 % worse`` describe the same mute with opposite signs.** They may be two different
yardsticks (``abs_err_all_final`` vs ``abs_err_all``), which ``arm_score.py`` prints side by side for
exactly this reason, but neither number says which it is. ⚠ A second candidate, and the arithmetic
argues for it rather than any measurement doing so: ``0.50`` and ``1.999`` are reciprocals, so the two
may be ONE ratio read in opposite directions — ``muted/relayed`` against ``relayed/muted`` — in which
case the defect is a comparison direction, not a yardstick. Both are preserved verbatim rather than
reconciled by guesswork.
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
