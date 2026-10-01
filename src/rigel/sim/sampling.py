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
    """Sample ``n`` integer fragment lengths from ``Normal(mean, std)`` truncated to
    ``[frag_min, frag_max]`` (inclusive), via rejection sampling against ``rng``."""
    result = np.empty(n, dtype=int)
    filled = 0
    while filled < n:
        needed = n - filled
        size = (
            int(needed * CONSTANTS.simulator.oversample_ratio)
            + CONSTANTS.simulator.oversample_extra
        )
        raw = rng.normal(mean, std, size).astype(int)
        valid = raw[(raw >= frag_min) & (raw <= frag_max)]
        nkeep = min(len(valid), needed)
        result[filled : filled + nkeep] = valid[:nkeep]
        filled += nkeep
    return result
