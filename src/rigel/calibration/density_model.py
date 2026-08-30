"""Signature-derived COUNT-OBSERVABILITY masks — which regions and which contiguous boundaries hold a
count that is gDNA rather than unspliced mature RNA.

A region is count-observable iff it carries no exon bit; a contiguous boundary iff its two flanks SHARE no
exon bit (a shared exon-strand means one exon continues across it, so unspliced mature RNA crosses it and
its count is not a gDNA count). These two masks are the gDNA strand-overdispersion fit's SEED SELECTOR
(:mod:`rigel.calibration.gdna_strand`).

⭐ **This module used to also impute a per-region gDNA DENSITY** — read off count-observable regions, carried
inward from flanking boundaries, run-filled across unanchored runs (:mod:`run_fill`), with a global
count-weighted baseline for a reference with no anchor at all. That density existed to weight the strand
fit's seeds. The fit moved to the away-half moment, which is unbiased under any RNA content and therefore
needs NO weight, leaving the density with no production consumer; it and `run_fill` were deleted on
2026-08-30 (owner ruling) rather than kept alive by their own tests.
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

    Returns ``(region_count_observable[N], boundary_count_observable[E])``.

    ⭐ **The boundary mask is on the BOUNDARY axis now.** It used to be a ``bool[R]`` in which entry ``r``
    described the boundary ``(r, r+1)`` and the last entry of every reference was a padding ``False`` — a
    region-shaped array standing in for an boundary-shaped one, with the ``same``-reference test carried
    inside it. A contiguous boundary is a first-class object with its own axis, so the mask is
    ``bool[E]`` and the ``same`` test is :func:`~rigel.calibration.region_arrays.boundary_region_indices`'s
    job rather than this function's.
    """
    sig = np.asarray(signature).astype(np.int64)
    region_count_observable = (sig & _EXON_BITS) == 0
    lo, hi = boundary_region_indices(ref_id)
    # observable ⇔ no exon bit is SHARED across the boundary ⇒ no single exon-strand continues across it
    # ⇒ no unspliced mature RNA crosses, so the crossing count is gDNA (+ nascent).
    boundary_count_observable = (sig[lo] & sig[hi] & _EXON_BITS) == 0
    return region_count_observable, boundary_count_observable


__all__ = ["count_observable_masks"]
