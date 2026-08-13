"""Shared region-geometry helpers: same-reference neighbour masks + bidirectional run-fill.

A stage that walks the region partition looks at each region's left/right neighbour *within the same
reference*, and fills unset regions by carrying the nearest anchored value inward from both ends of a run.
That bookkeeping lives here once.

⚠ **ONE consumer today — `density_model`.** This docstring used to claim four (`strand_deconv`, `priors`
and the chain sweep as well) and was measured wrong by `scripts/design/module_census.py`: the other three
had stopped calling it and nothing said so. A sentence about the code, inside the code, that nothing gates
is exactly how a reader is sent looking for something that is not there.

Both helpers treat the region array as a sequence of per-reference runs: a neighbour relation or a
carry never crosses a reference boundary (``ref_id`` change).
"""

from __future__ import annotations

import numpy as np

__all__ = ["same_ref_left_right", "runfill_bidirectional"]


def same_ref_left_right(ref_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(left_same, right_same)`` — does region ``i`` have a same-reference neighbour on each side?

    ``left_same[i]`` is True iff region ``i−1`` exists and shares ``i``'s reference (so ``i`` has a
    left neighbour); ``right_same[i]`` likewise for ``i+1``. The reference edges are False. Both are
    ``bool[R]``. This is the single source of the ``ref_id[1:] == ref_id[:-1]`` / ``[:-1] == [1:]``
    shifted-equality bookkeeping used across the package; the per-boundary form (length ``R−1``) is
    ``right_same[:-1]`` (equivalently ``left_same[1:]``).
    """
    ref = np.asarray(ref_id)
    r = ref.shape[0]
    left_same = np.zeros(r, dtype=bool)
    right_same = np.zeros(r, dtype=bool)
    if r > 1:
        eq = ref[:-1] == ref[1:]  # eq[k] = (ref[k] == ref[k+1])
        left_same[1:] = eq  # region i (>0) has a same-ref left neighbour iff ref[i] == ref[i-1]
        right_same[:-1] = (
            eq  # region i (<R-1) has a same-ref right neighbour iff ref[i] == ref[i+1]
        )
    return left_same, right_same


def runfill_bidirectional(values: np.ndarray, ref_id: np.ndarray) -> np.ndarray:
    """Fill unset (``nan``) regions by carrying the nearest anchored value inward from both run ends.

    Within each reference run, an unset region inherits the nearest set value reached from the left
    (forward pass) and from the right (reverse pass); where both reach it takes their **mean**, where
    one reaches it takes that one, where neither reaches (a run with no anchor anywhere) it stays
    ``nan`` — the caller supplies its own global fallback. Originally-set regions are returned
    unchanged. A carry never crosses a reference boundary.

    This is the exact bidirectional run-fill `density_model.region_gdna_density` uses (extracted verbatim);
    the BP chain sweep reuses it in fraction/state space. Cost is two ``O(R)`` passes.
    """
    v = np.asarray(values, dtype=np.float64)
    ref = np.asarray(ref_id)
    r = v.shape[0]
    fwd = v.copy()
    for i in range(1, r):
        if np.isnan(fwd[i]) and ref[i] == ref[i - 1]:
            fwd[i] = fwd[i - 1]
    rev = v.copy()
    for i in range(r - 2, -1, -1):
        if np.isnan(rev[i]) and ref[i] == ref[i + 1]:
            rev[i] = rev[i + 1]
    stack = np.vstack([fwd, rev])
    n_valid = (~np.isnan(stack)).sum(axis=0)
    carried = np.where(n_valid > 0, np.nansum(stack, axis=0) / np.maximum(n_valid, 1), np.nan)
    return np.where(np.isnan(v), carried, v)
