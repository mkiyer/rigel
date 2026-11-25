"""Test-only helpers for fragment resolution; not part of the shipped ``rigel`` package.

Production fragment resolution runs entirely in C++ (``rigel._resolve_impl``'s
``FragmentResolver``, reached as ``index.resolver``). This module provides:

- ``make_fragment()`` — a minimal fragment object carrying the three attributes
  the C++ resolver reads
- ``resolve_fragment()`` — a thin driver that hands such a fragment to the
  resolver's test seam; it does no resolution of its own
- ``GenomicInterval`` — an aligned block or intron on a reference and strand, the shape a fragment's
  exons and introns take
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

from types import SimpleNamespace
from typing import NamedTuple

from rigel._resolve_impl import testing as resolve_testing
from rigel.config import BamScanConfig
from rigel.types import Strand


class GenomicInterval(NamedTuple):
    """An interval on a reference and strand: an aligned block or an observed intron of a fragment."""

    ref: str
    start: int
    end: int
    strand: int = Strand.NONE


# ---------------------------------------------------------------------------
# Fragment construction helper
# ---------------------------------------------------------------------------


def make_fragment(exons=(), introns=()):
    """Create a minimal fragment-like object for :func:`resolve_fragment`.

    Builds a lightweight ``SimpleNamespace`` with ``.exons`` and ``.introns``,
    the two attributes the C++ resolve kernel reads via ``frag.attr()``.

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
    return SimpleNamespace(exons=tuple(exons), introns=tuple(introns))


# ---------------------------------------------------------------------------
# Fragment resolution — C++ native kernel (required)
# ---------------------------------------------------------------------------


def resolve_fragment(frag, index, *, intron_overhang_tolerance: int = 0):
    """Resolve a fragment to its compatible transcript set with the C++ resolver.

    Calls ``rigel._resolve_impl.testing.resolve_fragment`` on ``index.resolver``, which resolves it as the scan
    resolves a hit, and returns what it returns, a ``ResolvedFragment`` or ``None``; a fragment with no exon blocks is
    ``None`` without the call.
    """
    if not frag.exons:
        return None
    return resolve_testing.resolve_fragment(
        index.resolver,
        frag,
        max_fragment_length=BamScanConfig().max_frag_length,
        intron_overhang_tolerance=intron_overhang_tolerance,
    )
