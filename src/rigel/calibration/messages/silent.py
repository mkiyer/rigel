"""The policy that sends nothing — the OFF state (``message_policy = "silent"``) and the measured floor.

       Gate: ``tests/calibration/test_sweep_backbone.py``

The kernel runs no layer for it (`native/solve_kernel.cpp`): every node holds silence from each side it
has a neighbour on, and ψ solves every slot on its own evidence and the prior alone. Every message policy
is judged against this floor — win on unstranded data, minimal harm on stranded data, never pooled.
"""

from __future__ import annotations

__all__ = ["SilentPolicy"]


class SilentPolicy:
    name = "silent"
    strand = None

    def library(self, view):
        return None
