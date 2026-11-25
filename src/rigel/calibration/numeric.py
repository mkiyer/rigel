"""The numerical primitives several calibration layers share, with no calibration dependency (layer 0): bisection to
adjacent floats, and the mass a beta distribution puts on an interval."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.special import betainc, betaincc

__all__ = ["beta_span", "bisect_to_adjacent"]


def bisect_to_adjacent(
    root_above: Callable[[float], bool], lo: float, hi: float
) -> tuple[float, float]:
    """The bracket ``(lo, hi)`` once bisection has closed it to adjacent floats, ``root_above(mid)`` saying the root
    lies above ``mid``. At adjacent floats the midpoint rounds to one end, so no iteration budget is needed."""
    while True:
        mid = 0.5 * (lo + hi)
        if not lo < mid < hi:
            return lo, hi
        if root_above(mid):
            lo = mid
        else:
            hi = mid


def beta_span(a, b, x1, x2):
    """``I_x2(a, b) − I_x1(a, b)``, the mass ``Beta(a, b)`` puts on ``[x1, x2]``, never negative. The upper tail goes
    through the complements, which keep their digits where both are near 1."""
    lower = betainc(a, b, x2) - betainc(a, b, x1)
    upper = betaincc(a, b, x1) - betaincc(a, b, x2)
    return np.maximum(np.where(betainc(a, b, x1) > 0.5, upper, lower), 0.0)
