"""Test-only helpers for fragment resolution; NOT part of the shipped ``rigel`` package.

Production fragment resolution runs entirely in C++ (``rigel._resolve_impl``'s
``FragmentResolver``, reached as ``index.resolver``). This module provides:

- ``make_fragment()`` — a minimal fragment object carrying the three attributes
  the C++ resolver reads
- ``resolve_fragment()`` — a thin driver that hands such a fragment to
  ``index.resolver``; it does no resolution of its own
- ``_detect_intrachromosomal_chimera()`` — a Python twin of the intrachromosomal
  chimera rule (transcript-set disjointness, compatibility checked first), gated
  in ``test_resolution.py`` beside the C++ kernel
- ``MergeOutcome`` / ``ChimeraType`` — the resolver's merge-level and chimera codes
  (``native/constants.h``), named for the gates
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

from enum import IntEnum
from types import SimpleNamespace

from rigel.types import Strand


# ---------------------------------------------------------------------------
# The resolver's merge-level and chimera codes (`native/constants.h`)
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


class ChimeraType(IntEnum):
    """Classification of chimeric fragments.

    A chimeric fragment has exon blocks that map to disjoint transcript
    sets, indicating fusion or read-through events.

    Values
    ------
    NONE : int
        All exon blocks map to a connected set of transcripts.
    TRANS : int
        Exon blocks span multiple reference sequences (interchromosomal).
        Suggestive of trans-splicing or gene fusions.
    CIS_STRAND_SAME : int
        Intrachromosomal chimera where both disjoint exon-block
        clusters align to the same strand.  Suggestive of
        transcriptional read-through between adjacent genes.
    CIS_STRAND_DIFF : int
        Intrachromosomal chimera where the disjoint exon-block
        clusters align to different strands.  Suggestive of genomic
        rearrangement or trans-splicing.
    """

    NONE = 0
    TRANS = 1
    CIS_STRAND_SAME = 2
    CIS_STRAND_DIFF = 3


# ---------------------------------------------------------------------------
# Intrachromosomal chimera detection
# ---------------------------------------------------------------------------


def _detect_intrachromosomal_chimera(
    exon_blocks: tuple,
    exon_t_sets: list[frozenset[int]],
    max_fragment_length: int,
) -> ChimeraType | None:
    """Detect intrachromosomal chimeras via transcript-set disjointness.

    Two exon blocks are "connected" if their transcript sets share at
    least one transcript index.  If the non-empty sets form more than
    one connected component, the fragment is a chimera CANDIDATE.

    COMPATIBILITY IS CHECKED FIRST.  Disjoint
    transcript sets are not evidence of a rearrangement: a genomic
    molecule is contiguous and routinely spans two transcripts that
    share nothing.  A candidate is only a chimera if the mates are
    genomically INCOMPATIBLE — an orientation that is not facing inward
    (the two components carry different reference strands), or an
    implied fragment length beyond ``max_fragment_length``.  Otherwise
    it is an ordinary molecule and this returns ``None``.

    Parameters
    ----------
    exon_blocks : tuple of GenomicInterval
        The fragment's exon blocks (all on the same reference).
    exon_t_sets : list of frozenset[int]
        Per-block transcript index sets (parallel to *exon_blocks*).
    max_fragment_length : int
        The library's fragment-length limit (``BamScanConfig.max_frag_length``).
        A candidate whose implied fragment length — outermost start to
        outermost end, NOT the gap between blocks — is within this is a
        compatible molecule rather than a chimera.

    Returns
    -------
    ChimeraType or None
        The chimera type if chimeric, ``None`` otherwise.
    """
    # Filter to blocks with non-empty transcript sets
    items = [(block, tset) for block, tset in zip(exon_blocks, exon_t_sets) if tset]
    if len(items) <= 1:
        return None  # Cannot be chimeric with 0 or 1 annotated block

    n = len(items)

    # Union-find for connected components
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    # Connect blocks whose transcript sets intersect
    for i in range(n):
        for j in range(i + 1, n):
            if items[i][1] & items[j][1]:  # non-empty intersection
                union(i, j)

    # Group into connected components
    components: dict[int, list[int]] = {}
    for i in range(n):
        components.setdefault(find(i), []).append(i)

    if len(components) <= 1:
        return None  # All connected — not chimeric

    # --- Strand characterisation ---
    comp_strands: list[int] = []
    for members in components.values():
        strand = Strand.NONE
        for idx in members:
            strand |= items[idx][0].strand
        comp_strands.append(strand)

    # Compare strands across components
    unique_strands = set(comp_strands)
    if len(unique_strands) == 1:
        chimera_type = ChimeraType.CIS_STRAND_SAME
    else:
        chimera_type = ChimeraType.CIS_STRAND_DIFF

    # COMPATIBILITY BEFORE CHIMERA.  One reference is already given (the caller
    # gates on `is_interchromosomal`); facing inward is exactly `len(unique_strands) == 1`,
    # because `build_fragment` keys blocks by (ref, ref_strand) with R2's orientation
    # flipped, so an inward-facing pair lands both mates on ONE strand.  What remains is
    # whether a molecule of the implied length could exist in this library.
    if chimera_type == ChimeraType.CIS_STRAND_SAME:
        span = max(b.end for b, _ in items) - min(b.start for b, _ in items)
        if span <= max_fragment_length:
            return None

    return chimera_type


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
