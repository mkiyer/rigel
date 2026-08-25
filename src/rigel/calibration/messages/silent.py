"""The policy that sends nothing — the OFF state (``message_propagation = False``) and the measured floor.

       Gate: ``tests/calibration/test_sweep_backbone.py``

⭐⭐⭐ **This file is five boundaries of behaviour, and that is the entire point.** A new session reads
``sweep.py`` plus this and holds the whole working system in their head. :mod:`~.relay` is the SHIPPED
policy (``message_propagation`` defaults ``True`` since 2026-08-18), with every operator inside it behind
its own named switch — the switches exist so ``ladder_arm_ab.py`` can price the operators ONE AT A TIME.

⛔ **It is also a MEASURED floor, not a placeholder.** With no belief propagation whatsoever the panel's
deliverable is worth −50 % overall, and the sign is not uniform: it is a net **improvement** on three of
the four strata (stranded × capture-ON **−58.3 %, 16/16**; stranded × capture-OFF **−43.7 %, 16/16**;
unstranded × capture-OFF **−32.1 %**) and a large regression on exactly one (unstranded × capture-ON
**+154.8 %, 0/16**) — the stratum where ``kappa = 1/2`` makes the strand λ-term exactly 0, so a slot has no
own composition evidence at all and a message is the only source there is.

⭐ **AND IT IS NO LONGER SHIPPED — this paragraph asserted the opposite until 2026-08-23, the third
wrong-fact fossil of the same family** (`calibrate.py` and `zero_controls.py` carried the other two).
``message_propagation`` has defaulted ``True`` since 2026-08-18, installing :class:`~.relay.RelayPolicy`;
this policy is the OFF state and the measured floor the message layer is priced against. What the floor
establishes is unchanged and still load-bearing: the relay's whole value is ONE stratum (the deferred
unstranded × capture-ON), it is a net harm on the other three, and that is what scoped the backbone —
and the reference every message policy is measured against.

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
    its spliced count, the fitted gDNA prior and the intron factory."""

    name = "silent"

    def prepare(self, ctx: StepContext) -> _SilentRelay:
        return _SilentRelay()
