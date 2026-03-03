"""
hulkrna.resolution — Fragment resolution and chimera detection.

Resolves aligned fragments to their compatible transcript/gene sets
by querying the reference index and applying progressive set merging.

Uses the C++ native kernel (``hulkrna._resolve_impl``) exclusively.

This module contains:
- ``make_fragment()`` — lightweight fragment constructor for ``resolve_fragment``
- ``resolve_fragment()`` — core fragment-to-transcript resolution with
  chimera detection (interchromosomal and intrachromosomal)
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

from types import SimpleNamespace

from .types import (
    ChimeraType,
    Strand,
)



# ---------------------------------------------------------------------------
# Intrachromosomal chimera detection
# ---------------------------------------------------------------------------

def _detect_intrachromosomal_chimera(
    exon_blocks: tuple,
    exon_t_sets: list[frozenset[int]],
) -> tuple[ChimeraType, int] | None:
    """Detect intrachromosomal chimeras via transcript-set disjointness.

    Two exon blocks are "connected" if their transcript sets share at
    least one transcript index.  If the non-empty sets form more than
    one connected component, the fragment is chimeric.

    Parameters
    ----------
    exon_blocks : tuple of GenomicInterval
        The fragment's exon blocks (all on the same reference).
    exon_t_sets : list of frozenset[int]
        Per-block transcript index sets (parallel to *exon_blocks*).

    Returns
    -------
    tuple[ChimeraType, int] or None
        ``(chimera_type, chimera_gap)`` if chimeric, ``None`` otherwise.
        ``chimera_gap`` is the minimum genomic distance (bp) between
        exon blocks in different connected components.
    """
    # Filter to blocks with non-empty transcript sets
    items = [
        (block, tset)
        for block, tset in zip(exon_blocks, exon_t_sets)
        if tset
    ]
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

    # --- Compute minimum gap between components ---
    comp_list = list(components.values())
    min_gap = float("inf")
    for ci in range(len(comp_list)):
        for cj in range(ci + 1, len(comp_list)):
            for bi in comp_list[ci]:
                for bj in comp_list[cj]:
                    blk_i = items[bi][0]
                    blk_j = items[bj][0]
                    if blk_i.end <= blk_j.start:
                        gap = blk_j.start - blk_i.end
                    elif blk_j.end <= blk_i.start:
                        gap = blk_i.start - blk_j.end
                    else:
                        gap = 0  # overlapping
                    min_gap = min(min_gap, gap)

    return chimera_type, int(min_gap)


# ---------------------------------------------------------------------------
# Fragment construction helper
# ---------------------------------------------------------------------------

def make_fragment(exons=(), introns=()):
    """Create a minimal fragment-like object for :func:`resolve_fragment`.

    Builds a lightweight ``SimpleNamespace`` with ``.exons``, ``.introns``,
    and ``.genomic_footprint`` — the three attributes the C++ resolve
    kernel reads via ``frag.attr()``.

    Replaces the former ``Fragment`` dataclass (now removed).

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
        exons=exons, introns=introns, genomic_footprint=footprint,
    )


# ---------------------------------------------------------------------------
# Fragment resolution — C++ native kernel (required)
# ---------------------------------------------------------------------------

def resolve_fragment(frag, index):
    """Resolve a fragment to its compatible transcript set.

    Uses the C++ native kernel (``hulkrna._resolve_impl``) via
    ``ResolveContext.resolve_fragment()``.  Returns a C++
    ``ResolvedResult`` object that exposes attributes for model
    training and buffering.
    """
    if not frag.exons:
        return None
    return index._resolve_ctx.resolve_fragment(frag)

