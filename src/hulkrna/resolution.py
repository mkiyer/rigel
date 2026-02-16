"""
hulkrna.resolution — Fragment resolution and insert-size computation.

Resolves aligned fragments to their compatible transcript/gene sets
by querying the reference index and applying progressive set merging.
Also computes insert sizes with gap correction.

This module contains:
- ``merge_sets_with_criteria()`` — progressive set merging
- ``resolve_fragment()`` — core fragment-to-transcript resolution with
  chimera detection (interchromosomal and intrachromosomal)
- ``compute_insert_size()`` — insert size with gap correction
- ``fragment_insert_size()`` — simple insert size (no transcript context)
- ``ResolvedFragment`` — compact cached resolution result
"""

import logging
from dataclasses import dataclass

from .types import (
    EMPTY_MERGE,
    ChimeraType,
    IntervalType,
    MergeCriteria,
    MergeResult,
    Strand,
)
from .categories import CountCategory
from .fragment import Fragment

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Set merging
# ---------------------------------------------------------------------------

def merge_sets_with_criteria(
    t_sets: list[frozenset[int]],
) -> MergeResult:
    """Merge transcript index sets with progressive relaxation.

    Strategy
    --------
    1. **Intersection** of all sets (most specific).
    2. **Intersection of non-empty sets** — skip empty sets that would
       collapse everything to the empty set.
    3. **Union** of all sets (most sensitive).

    Parameters
    ----------
    t_sets : list of frozenset[int]
        Per-segment transcript index sets.

    Returns
    -------
    MergeResult
    """
    if not t_sets:
        return EMPTY_MERGE

    # 1. Intersection of all sets (most specific)
    t_inds = frozenset.intersection(*t_sets) if t_sets else frozenset()
    if t_inds:
        return MergeResult(t_inds, MergeCriteria.INTERSECTION)

    # 2. Intersection of non-empty sets
    t_nonempty = [s for s in t_sets if s]
    t_inds = (
        frozenset.intersection(*t_nonempty) if t_nonempty else frozenset()
    )
    if t_inds:
        return MergeResult(t_inds, MergeCriteria.INTERSECTION_NONEMPTY)

    # 3. Union of all sets (most sensitive)
    t_inds = frozenset.union(*t_sets) if t_sets else frozenset()
    if t_inds:
        return MergeResult(t_inds, MergeCriteria.UNION)

    return EMPTY_MERGE


# ---------------------------------------------------------------------------
# Insert size computation
# ---------------------------------------------------------------------------

def fragment_insert_size(frag: Fragment) -> int:
    """Compute insert size as footprint minus observed introns.

    No gap correction (no transcript context).  Used for intergenic
    fragments where there are no compatible transcripts.

    Returns
    -------
    int
        Insert size in bp, or -1 if the fragment has no exon blocks.
    """
    if not frag.exons:
        return -1
    sorted_exons = sorted((e.start, e.end) for e in frag.exons)
    footprint = sorted_exons[-1][1] - sorted_exons[0][0]
    observed_intron_total = sum(
        intron.end - intron.start for intron in frag.introns
    )
    return footprint - observed_intron_total


def compute_insert_size(
    frag: Fragment,
    compatible_t_inds: frozenset[int],
    index,
) -> int:
    """Compute insert size with gap correction for compatible transcripts.

    Returns the unambiguous insert size if all compatible transcripts
    agree, or -1 if they disagree or the fragment is empty.

    **Algorithm**:

    1. Footprint = max(end) - min(start) of fragment exon blocks.
    2. Subtract observed intron sizes (CIGAR N-ops) -> upper bound.
    3. Identify gaps between consecutive exon blocks that are *not*
       observed introns.
    4. For each compatible transcript, find reference introns fully
       contained in those gaps -> subtract from upper bound.
    5. If all transcripts yield the same positive size, return it;
       otherwise return -1.

    Parameters
    ----------
    frag : Fragment
        The consolidated fragment (must be single-ref).
    compatible_t_inds : frozenset[int]
        Transcript indices surviving resolution.
    index : HulkIndex
        The loaded reference index.

    Returns
    -------
    int
        Unambiguous insert size, or -1.
    """
    if not frag.exons or not compatible_t_inds:
        return -1

    ref = frag.exons[0].ref
    sorted_exons = sorted((e.start, e.end) for e in frag.exons)
    footprint = sorted_exons[-1][1] - sorted_exons[0][0]

    # Observed introns (CIGAR N-ops) on the fragment's reference
    observed_introns = set(
        (intron.start, intron.end)
        for intron in frag.introns
        if intron.ref == ref
    )
    observed_intron_total = sum(e - s for s, e in observed_introns)
    upper_bound = footprint - observed_intron_total

    # Gaps between consecutive exon blocks that aren't observed introns
    gaps: list[tuple[int, int]] = []
    for i in range(len(sorted_exons) - 1):
        gap_start = sorted_exons[i][1]
        gap_end = sorted_exons[i + 1][0]
        if gap_start < gap_end and (gap_start, gap_end) not in observed_introns:
            gaps.append((gap_start, gap_end))

    if not gaps:
        return upper_bound

    # Per-transcript gap correction
    t_gap_intron_size: dict[int, int] = {}
    for gap_start, gap_end in gaps:
        for t_idx, _g_idx, _strand, sj_start, sj_end in index.query_gap_sjs(
            ref, gap_start, gap_end
        ):
            if t_idx in compatible_t_inds:
                t_gap_intron_size[t_idx] = (
                    t_gap_intron_size.get(t_idx, 0) + (sj_end - sj_start)
                )

    # Check if all compatible transcripts agree
    insert_values: set[int] = set()
    for t_idx in compatible_t_inds:
        size = upper_bound - t_gap_intron_size.get(t_idx, 0)
        if size > 0:
            insert_values.add(size)

    if len(insert_values) == 1:
        return insert_values.pop()
    return -1


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
        chimera_type = ChimeraType.STRAND_SAME
    else:
        chimera_type = ChimeraType.STRAND_DIFF

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
# ResolvedFragment — compact cached resolution result
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ResolvedFragment:
    """Cached result of resolving one fragment against the reference index.

    Stores everything needed for model training (Phase 1), unique
    counting (Phase 2), and Bayesian reassignment (Phase 3) so that
    the BAM does not need to be re-read.

    Gene information is derived from transcript indices via
    ``index.t_to_g_arr`` rather than stored redundantly.

    Attributes
    ----------
    t_inds : frozenset[int]
        Compatible transcript indices.
    n_genes : int
        Number of distinct genes derived from ``t_inds`` via
        ``t_to_g_arr``.
    count_cat : CountCategory
        INTRON / UNSPLICED / SPLICED_UNANNOT / SPLICED_ANNOT.
    exon_strand : Strand
        Combined alignment strand of the fragment's exon blocks
        (after R2 flip).
    sj_strand : Strand
        Combined splice-junction strand from XS tags.
    insert_size : int
        Computed insert size (-1 if ambiguous or unavailable).
    merge_criteria : MergeCriteria
        Which relaxation level produced the transcript sets.
    num_hits : int
        Number of alignment pairs for this read name (NH tag).
    chimera_type : ChimeraType
        Chimera classification (NOT_CHIMERIC / INTERCHROMOSOMAL /
        STRAND_SAME / STRAND_DIFF).
    chimera_gap : int
        Minimum genomic distance between disjoint exon-block clusters
        for intrachromosomal chimeras, -1 for interchromosomal or
        non-chimeric fragments.
    """
    t_inds: frozenset[int]
    n_genes: int
    count_cat: CountCategory
    exon_strand: Strand
    sj_strand: Strand
    insert_size: int
    merge_criteria: MergeCriteria
    num_hits: int
    chimera_type: ChimeraType = ChimeraType.NOT_CHIMERIC
    chimera_gap: int = -1

    @property
    def is_chimeric(self) -> bool:
        """True if this fragment is classified as chimeric."""
        return self.chimera_type != ChimeraType.NOT_CHIMERIC

    @property
    def is_unique_gene(self) -> bool:
        """True if exactly one gene index in the compatible set."""
        return self.n_genes == 1

    @property
    def is_strand_qualified(self) -> bool:
        """True if this fragment qualifies for strand model training.

        Requirements:
        1. Has annotated splice junction(s).
        2. Resolves to a unique gene.
        3. Exon strand is unambiguous (POS or NEG).
        4. SJ strand is unambiguous (POS or NEG).
        5. Not chimeric.
        """
        return (
            self.count_cat == CountCategory.SPLICED_ANNOT
            and self.is_unique_gene
            and self.exon_strand in (Strand.POS, Strand.NEG)
            and self.sj_strand in (Strand.POS, Strand.NEG)
            and not self.is_chimeric
        )


# ---------------------------------------------------------------------------
# Fragment resolution
# ---------------------------------------------------------------------------

def resolve_fragment(
    frag: Fragment,
    index,
) -> ResolvedFragment | None:
    """Resolve a fragment to its compatible transcript set.

    Queries the index for overlapping intervals and splice junctions,
    merges transcript sets via progressive relaxation, classifies the
    overlap category, and detects chimeric fragments.

    Chimera detection
    -----------------
    - **Interchromosomal**: exon blocks span multiple references.
    - **Intrachromosomal**: exon blocks on a single reference whose
      per-block transcript sets form disjoint connected components
      (no shared transcript between clusters).

    Intrachromosomal chimeras are further classified as
    ``STRAND_SAME`` or ``STRAND_DIFF`` based on the alignment strand
    of the disjoint exon-block clusters.

    Parameters
    ----------
    frag : Fragment
        The consolidated fragment from paired-end reads.
    index : HulkIndex
        The loaded reference index.

    Returns
    -------
    ResolvedFragment or None
        ``None`` for intergenic or empty fragments.
    """
    if not frag.exons:
        return None

    # --- Interchromosomal chimera detection ---
    chimera_type = ChimeraType.NOT_CHIMERIC
    chimera_gap = -1

    refs = set(e.ref for e in frag.exons)
    is_interchromosomal = len(refs) > 1
    if is_interchromosomal:
        chimera_type = ChimeraType.INTERCHROMOSOMAL

    # --- Query each exon block against the unified cgranges index ---
    exon_t_sets: list[frozenset[int]] = []
    intron_t_sets: list[frozenset[int]] = []

    exon_strand = Strand.NONE

    for exon_block in frag.exons:
        exon_strand |= exon_block.strand

        block_exon_t: set[int] = set()
        block_intron_t: set[int] = set()

        for t_idx, g_idx, itype in index.query_exon(exon_block):
            if itype == IntervalType.EXON:
                block_exon_t.add(t_idx)
            elif itype == IntervalType.INTRON:
                block_intron_t.add(t_idx)

        exon_t_sets.append(frozenset(block_exon_t))
        intron_t_sets.append(frozenset(block_intron_t))

    # --- Intrachromosomal chimera detection (single-ref only) ---
    if not is_interchromosomal:
        chimera_result = _detect_intrachromosomal_chimera(
            frag.exons, exon_t_sets,
        )
        if chimera_result is not None:
            chimera_type, chimera_gap = chimera_result

    # --- Query SJ matches ---
    sj_t_sets: list[frozenset[int]] = []
    has_annotated_sj = False
    has_unannotated_sj = False
    sj_strand = Strand.NONE

    for intron in frag.introns:
        key = (intron.ref, intron.start, intron.end, intron.strand)
        match = index.sj_map.get(key)
        if match is not None:
            t_set, g_set = match
            sj_t_sets.append(t_set)
            has_annotated_sj = True
            sj_strand |= intron.strand
        else:
            has_unannotated_sj = True

    # --- Resolution ---
    any_exon = any(s for s in exon_t_sets)

    if any_exon:
        # Merge EXON + SJ sets together for maximum specificity
        all_t_sets = exon_t_sets + sj_t_sets
        merge_result = merge_sets_with_criteria(all_t_sets)

        # Determine count category
        if has_annotated_sj:
            count_cat = CountCategory.SPLICED_ANNOT
        elif has_unannotated_sj:
            count_cat = CountCategory.SPLICED_UNANNOT
        else:
            count_cat = CountCategory.UNSPLICED
    else:
        # No EXON overlap → try INTRON fallback
        any_intron = any(s for s in intron_t_sets)
        if any_intron:
            merge_result = merge_sets_with_criteria(intron_t_sets)
            count_cat = CountCategory.INTRON
        else:
            # Pure intergenic — no compatible hits
            return None

    if merge_result.is_empty:
        # Defensive — shouldn't happen with valid input
        return None

    # --- Derive n_genes from transcript indices ---
    t_inds = merge_result.t_inds
    n_genes = len(set(int(index.t_to_g_arr[t]) for t in t_inds))

    # --- Compute insert size (skip for chimeras — meaningless) ---
    if chimera_type == ChimeraType.NOT_CHIMERIC:
        insert_size = compute_insert_size(frag, t_inds, index)
    else:
        insert_size = -1

    return ResolvedFragment(
        t_inds=t_inds,
        n_genes=n_genes,
        count_cat=count_cat,
        exon_strand=exon_strand,
        sj_strand=sj_strand,
        insert_size=insert_size,
        merge_criteria=merge_result.criteria,
        num_hits=1,  # Caller must set this from the pairs list length
        chimera_type=chimera_type,
        chimera_gap=chimera_gap,
    )
