"""Shared fragment-length sampling for the simulator's read engine.

The single shared fragment-length sampler: :func:`truncated_normal_frag_lengths` draws
fragment lengths from a truncated normal via a rejection loop. Its one caller is the
whole-genome read engine (:class:`wgs_engine.WholeGenomeSimulator`).
"""

from __future__ import annotations

import numpy as np

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
    truncated to ``[frag_min, frag_max]`` (inclusive), via rejection sampling against ``rng``. Rounding, never
    truncation toward zero, is what makes the law `whole_genome.fl_pmf` states the one drawn."""
    if frag_min > frag_max:
        # The rejection loop would never fill: no draw lands in an empty window.
        raise ValueError(
            f"The fragment-length window [{frag_min}, {frag_max}] holds no length, so no fragment can be "
            "drawn: frag_min must not exceed frag_max (an RNA window's frag_max is capped at the longest "
            "transcript)."
        )
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
