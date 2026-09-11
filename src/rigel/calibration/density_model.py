"""Signature-derived count-observability masks: which regions and which contiguous boundaries hold a count
that is a gDNA count, rather than one unspliced mature RNA can also contribute to.

A region is count-observable iff it carries no exon bit. A contiguous boundary is count-observable iff its
two flanks share no exon bit — a shared exon-strand means one exon continues across the boundary, so
unspliced mature RNA crosses it and the crossing count is not a gDNA count.

Both masks are pure functions of the signature array and the reference ids: no count, no density and no
solve enters here. They are the seed selector for the gDNA strand-overdispersion fit
(:mod:`rigel.calibration.gdna_strand`), which needs seeds whose strand split is gDNA's own.
"""

from __future__ import annotations


import numpy as np

from .region_arrays import boundary_region_indices
from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG


def count_observable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-observability, on the two axes it describes.

    Returns ``(region_count_observable[N], boundary_count_observable[E])``. Each mask is on its own
    object's axis — a contiguous boundary is a first-class object, so the same-reference test that pairs
    a boundary with its two flanks belongs to
    :func:`~rigel.calibration.region_arrays.boundary_region_indices` rather than to this function.
    """
    sig = np.asarray(signature).astype(np.int64)
    region_count_observable = (sig & _EXON_BITS) == 0
    lo, hi = boundary_region_indices(ref_id)
    # observable ⇔ no exon bit is SHARED across the boundary ⇒ no single exon-strand continues across it
    # ⇒ no unspliced mature RNA crosses, so the crossing count is gDNA (+ nascent).
    boundary_count_observable = (sig[lo] & sig[hi] & _EXON_BITS) == 0
    return region_count_observable, boundary_count_observable


__all__ = ["count_observable_masks"]
