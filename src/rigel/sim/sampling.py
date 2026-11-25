"""The simulator's fragment-length sampler, a rounded normal truncated to a window, drawn by rejection where the window
holds enough of the law and exactly from its integer masses where it does not; its caller is
:class:`wgs_engine.WholeGenomeSimulator`."""

from __future__ import annotations

import numpy as np
from scipy.special import log_ndtr, ndtr

from ..config import CONSTANTS


def truncated_normal_frag_lengths(
    rng: np.random.Generator,
    n: int,
    mean: float,
    std: float,
    frag_min: int,
    frag_max: int,
) -> np.ndarray:
    """Sample ``n`` integer fragment lengths: a ``Normal(mean, std)`` draw rounded to the nearest integer,
    truncated to ``[frag_min, frag_max]`` (inclusive), by rejection against ``rng``. The draw must round, not
    truncate toward zero, so that the law drawn is the rounded normal's mass at each integer, truncated: the closed
    form the instruments integrate."""
    if frag_min > frag_max:
        raise ValueError(
            f"The fragment-length window [{frag_min}, {frag_max}] holds no length, so no fragment can be "
            "drawn: frag_min must not exceed frag_max (an RNA window's frag_max is capped at the longest "
            "transcript)."
        )
    if not (np.isfinite(mean) and np.isfinite(std) and std >= 0.0):
        raise ValueError(
            f"A fragment-length law needs a finite mean and a finite, nonnegative spread; got {mean}, {std}"
        )
    if std == 0.0:
        # a point mass: every length is the rounded mean, which the window holds or no fragment can be drawn
        lengths, _ = window_law(mean, std, frag_min, frag_max)
        return np.full(n, int(lengths[0]), dtype=int)
    in_window = float(ndtr((frag_max + 0.5 - mean) / std) - ndtr((frag_min - 0.5 - mean) / std))
    if in_window < CONSTANTS.simulator.min_rejection_acceptance:
        return _window_masses_drawn(rng, n, mean, std, frag_min, frag_max)
    result = np.empty(n, dtype=int)
    filled = 0
    while filled < n:
        needed = n - filled
        size = (
            int(needed * CONSTANTS.simulator.oversample_ratio)
            + CONSTANTS.simulator.oversample_extra
        )
        raw = np.rint(rng.normal(mean, std, size)).astype(int)
        valid = raw[(raw >= frag_min) & (raw <= frag_max)]
        nkeep = min(len(valid), needed)
        result[filled : filled + nkeep] = valid[:nkeep]
        filled += nkeep
    return result


def window_law(
    mean: float, std: float, frag_min: int, frag_max: int
) -> tuple[np.ndarray, np.ndarray]:
    """``(lengths, pmf)``: the rounded normal's mass at each integer of the window, renormalised, the law the sampler
    draws. Its half-integer bin edges are read in log space from the near tail, so a window far into a tail keeps its
    relative masses; a law with no spread is its rounded mean. Raises when the window holds none of the law."""
    if std == 0.0:
        point = int(np.rint(mean))
        if not frag_min <= point <= frag_max:
            raise ValueError(
                f"A fragment-length law with no spread puts every fragment at {point} bp, outside the window "
                f"[{frag_min}, {frag_max}]: no fragment can be drawn."
            )
        return np.array([point]), np.array([1.0])
    lengths = np.arange(frag_min, frag_max + 1)
    lo, hi = (lengths - 0.5 - mean) / std, (lengths + 0.5 - mean) / std
    if frag_min > mean:
        # the upper tail: Q(lo) − Q(hi), Q the survival function, read as Φ of the negated edges
        log_near, log_far = log_ndtr(-lo), log_ndtr(-hi)
    else:
        log_near, log_far = log_ndtr(hi), log_ndtr(lo)
    with np.errstate(divide="ignore"):
        log_mass = log_near + np.log1p(-np.exp(log_far - log_near))
    if not np.isfinite(log_mass).any():
        raise ValueError(
            f"The fragment-length law Normal({mean}, {std}) holds no mass in the window [{frag_min}, {frag_max}] at "
            "double precision: no fragment can be drawn."
        )
    mass = np.exp(log_mass - log_mass[np.isfinite(log_mass)].max())
    return lengths, mass / mass.sum()


def _window_masses_drawn(
    rng, n: int, mean: float, std: float, frag_min: int, frag_max: int
) -> np.ndarray:
    """``n`` lengths drawn exactly from the window's law (:func:`window_law`)."""
    lengths, pmf = window_law(mean, std, frag_min, frag_max)
    return rng.choice(lengths, size=n, p=pmf)
