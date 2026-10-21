"""Test-only helpers for fragment resolution; NOT part of the shipped ``rigel`` package.

Production fragment resolution runs entirely in C++ (``rigel._resolve_impl``'s
``FragmentResolver``, reached as ``index.resolver``). This module provides:

- ``make_fragment()`` — a minimal fragment object carrying the three attributes
  the C++ resolver reads
- ``resolve_fragment()`` — a thin driver that hands such a fragment to
  ``index.resolver``; it does no resolution of its own
- ``MergeOutcome`` — the resolver's merge-level codes (``native/constants.h``), named for the
  gates
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

from enum import IntEnum
from types import SimpleNamespace


# ---------------------------------------------------------------------------
# The resolver's merge-level codes (`native/constants.h`)
# ---------------------------------------------------------------------------


class MergeOutcome(IntEnum):
    """Which relaxation level succeeded in progressive set merging.

    During fragment resolution, transcript/gene index sets from
    individual exon blocks or splice junctions are merged via
    progressive relaxation:

    0. INTERSECTION — intersection of *all* sets (most specific)
    1. INTERSECTION_NONEMPTY — intersection of non-empty sets only
    2. UNION — union of all sets (most sensitive)
    3. EMPTY — no sets to merge (no hits of this type)
    """

    INTERSECTION = 0
    INTERSECTION_NONEMPTY = 1
    UNION = 2
    EMPTY = 3


# ---------------------------------------------------------------------------
# Fragment construction helper
# ---------------------------------------------------------------------------


def make_fragment(exons=(), introns=()):
    """Create a minimal fragment-like object for :func:`resolve_fragment`.

    Builds a lightweight ``SimpleNamespace`` with ``.exons``, ``.introns``,
    and ``.genomic_footprint`` — the three attributes the C++ resolve
    kernel reads via ``frag.attr()``.

    Parameters
    ----------
    exons : iterable of GenomicInterval
        Exon blocks (will be coerced to a tuple).
    introns : iterable of GenomicInterval
        Splice junctions (will be coerced to a tuple).

    Returns
    -------
    SimpleNamespace
        Object accepted by ``resolve_fragment`` and the C++ kernel.
    """
    exons = tuple(exons)
    introns = tuple(introns)
    footprint = exons[-1].end - exons[0].start if exons else -1
    return SimpleNamespace(
        exons=exons,
        introns=introns,
        genomic_footprint=footprint,
    )


# ---------------------------------------------------------------------------
# Fragment resolution — C++ native kernel (required)
# ---------------------------------------------------------------------------


def resolve_fragment(frag, index):
    """Resolve a fragment to its compatible transcript set with the C++ resolver.

    Calls ``index.resolver.resolve_fragment(frag)`` (``rigel._resolve_impl``'s
    ``FragmentResolver``) and returns what it returns, a ``ResolvedFragment`` or
    ``None``; a fragment with no exon blocks is ``None`` without the call.
    """
    if not frag.exons:
        return None
    return index.resolver.resolve_fragment(frag)
