"""rigel.calibration._locus_n_obs — Anchor-transcript routing of EM units
to the constituent ``Locus`` intervals of a ``MultiLocus``.

Per ``docs/calibration/calibration_v6_plan.md`` §2.8:

* Fast path (single-``Locus`` ``MultiLocus``, ≈99% of components):
  return ``ml.unit_indices`` directly with zero allocation.
* Slow path: gather each unit's anchor transcript from
  ``em_data.locus_t_indices``, look up the local ``Locus`` index via
  the per-``MultiLocus`` ``t_to_local_locus`` map (indexed by the
  local position within ``ml.transcript_indices``), and bin units
  per ``Locus``.

**Anti-pattern explicitly forbidden:** materializing per-fragment
genomic coordinates from the buffer to bin units geometrically.  See
``/memories/repo/multilocus-partition-design-2026-04.md``.
"""

from __future__ import annotations

import numpy as np

from ..locus import MultiLocus
from ..scored_fragments import ScoredFragments


__all__ = ["partition_units_to_loci", "build_t_to_local_locus"]


def build_t_to_local_locus(ml: MultiLocus, t_starts: np.ndarray, t_ref: np.ndarray) -> np.ndarray:
    """Build a per-``MultiLocus`` map: local transcript position → local
    ``Locus`` index.

    ``t_starts`` and ``t_ref`` are global per-transcript arrays
    (``index.t_df["start"].values``, etc.).  The returned array has
    length ``len(ml.transcript_indices)`` and dtype ``int32``; each
    entry is the index into ``ml.loci`` containing that transcript's
    start coordinate.

    Returns an empty ``int32[]`` for single-``Locus`` MultiLoci (the
    fast path of :func:`partition_units_to_loci` does not consult it).
    """
    n_loci = len(ml.loci)
    if n_loci == 1:
        return np.zeros(0, dtype=np.int32)

    t_idx = ml.transcript_indices
    ts = t_starts[t_idx]
    tr = t_ref[t_idx]
    out = np.full(t_idx.size, -1, dtype=np.int32)
    for j, loc in enumerate(ml.loci):
        m = (tr == loc.ref_id) & (ts >= loc.start) & (ts < loc.end)
        out[m] = j
    return out


def partition_units_to_loci(
    ml: MultiLocus,
    em_data: ScoredFragments,
    t_to_local_locus: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """Route the units of ``ml`` to its constituent ``Locus`` intervals.

    Parameters
    ----------
    ml
        The multi-locus to partition.
    em_data
        Global ScoredFragments providing ``locus_t_indices`` (anchor
        transcript per unit).
    t_to_local_locus
        Output of :func:`build_t_to_local_locus` for ``ml``.  Ignored
        on the single-``Locus`` fast path.

    Returns
    -------
    tuple[np.ndarray, ...]
        One ``int32[]`` per ``Locus`` in ``ml.loci``, in the same order.
        Fast path returns ``(ml.unit_indices,)`` (the original array).
    """
    n_loci = len(ml.loci)
    if n_loci == 1:
        return (ml.unit_indices,)

    anchor_t = em_data.locus_t_indices[ml.unit_indices]
    # ml.transcript_indices is sorted ascending by build_multi_loci ⇒
    # searchsorted gives the local position of each anchor transcript.
    local_pos = np.searchsorted(ml.transcript_indices, anchor_t)
    if local_pos.size > 0:
        # Bounds + identity check: every anchor must be present in
        # ml.transcript_indices (locus invariant).
        if (local_pos >= ml.transcript_indices.size).any() or (
            ml.transcript_indices[np.minimum(local_pos, ml.transcript_indices.size - 1)]
            != anchor_t
        ).any():
            raise RuntimeError(
                f"partition_units_to_loci: anchor transcript not in "
                f"transcript_indices of MultiLocus(multi_locus_id="
                f"{ml.multi_locus_id}). build_multi_loci invariant violation."
            )

    bins = t_to_local_locus[local_pos]
    if (bins < 0).any():
        raise RuntimeError(
            f"partition_units_to_loci: anchor transcript not in any Locus "
            f"of MultiLocus(multi_locus_id={ml.multi_locus_id}). "
            f"build_multi_loci invariant violation."
        )

    unit_indices = ml.unit_indices
    return tuple(
        np.ascontiguousarray(unit_indices[bins == j], dtype=np.int32)
        for j in range(n_loci)
    )
