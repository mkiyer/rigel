"""The message policy that sends nothing (``message_policy = "silent"``): the floor every policy is judged against.

The kernel runs no message layer for it (``native/solve_kernel.cpp``), so ψ solves every slot on its own evidence and
the prior alone.

Gate: ``tests/calibration/test_sweep_backbone.py``.
"""

from __future__ import annotations

__all__ = ["SilentPolicy"]


class SilentPolicy:
    name = "silent"

    def library(self, view):
        return None
