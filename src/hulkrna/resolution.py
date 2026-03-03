"""
hulkrna.resolution — Fragment resolution and fragment-length computation.

Resolves aligned fragments to their compatible transcript/gene sets
by querying the reference index and applying progressive set merging.
Also computes fragment lengths with gap correction.

Uses the C++ native kernel (``hulkrna._resolve_impl``) exclusively.

This module contains:
- ``merge_sets_with_criteria()`` — progressive set merging
- ``resolve_fragment()`` — core fragment-to-transcript resolution with
  chimera detection (interchromosomal and intrachromosomal)
- ``compute_frag_lengths()`` — per-transcript fragment lengths with gap correction
- ``ResolvedFragment`` — compact cached resolution result
- ``pair_multimapper_reads()`` — secondary alignment pairing (annotate pass)
"""

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


from dataclasses import dataclass

from .types import (
    EMPTY_MERGE,
    ChimeraType,
    MergeCriteria,
    MergeResult,
    Strand,
)
from .categories import SpliceType
from .fragment import Fragment


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
# Fragment length computation
# ---------------------------------------------------------------------------


def compute_frag_lengths(
    frag: Fragment,
    compatible_t_inds: frozenset[int],
    index,
) -> dict[int, int]:
    """Compute per-transcript SJ-corrected fragment lengths.

    Each compatible transcript may have different annotated introns in
    the gaps between the fragment's exon blocks, yielding a different
    fragment length.  Returning a per-transcript dict (rather than a single
    consensus value) allows the EM to score each mRNA candidate with
    its own fragment-length probability.

    **Algorithm**:

    1. Footprint = max(end) - min(start) of fragment exon blocks.
    2. Subtract observed intron sizes (CIGAR N-ops) → upper bound.
    3. Identify gaps between consecutive exon blocks that are *not*
       observed introns.
    4. For each compatible transcript, find reference introns fully
       contained in those gaps → subtract from upper bound.
    5. Return ``{t_idx: corrected_size}`` for each transcript with
       a positive corrected size; omit transcripts whose corrected
       size is ≤ 0 (geometrically impossible).

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
    dict[int, int]
        Mapping ``t_idx → SJ-corrected fragment length`` for each
        transcript with a positive size.  Empty if no exon blocks
        or no compatible transcripts.
    """
    if not frag.exons or not compatible_t_inds:
        return {}

    # A7: Fast path for single exon block (no gaps, shared length)
    if len(frag.exons) == 1:
        fl = frag.exons[0].end - frag.exons[0].start
        if fl > 0:
            return {t_idx: fl for t_idx in compatible_t_inds}
        return {}

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
        # No unaccounted gaps → all transcripts share the same size
        if upper_bound > 0:
            return {t_idx: upper_bound for t_idx in compatible_t_inds}
        return {}

    # Per-transcript gap correction
    t_gap_intron_size: dict[int, int] = {}
    for gap_start, gap_end in gaps:
        for t_idx, _strand, sj_start, sj_end in index.query_gap_sjs(
            ref, gap_start, gap_end
        ):
            if t_idx in compatible_t_inds:
                t_gap_intron_size[t_idx] = (
                    t_gap_intron_size.get(t_idx, 0) + (sj_end - sj_start)
                )

    # Build per-transcript fragment lengths (omit non-positive)
    result: dict[int, int] = {}
    for t_idx in compatible_t_inds:
        size = upper_bound - t_gap_intron_size.get(t_idx, 0)
        if size > 0:
            result[t_idx] = size
    return result


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
    splice_type : SpliceType
        UNSPLICED / SPLICED_UNANNOT / SPLICED_ANNOT.
    exon_strand : Strand
        Combined alignment strand of the fragment's exon blocks
        (after R2 flip).
    sj_strand : Strand
        Combined splice-junction strand from XS tags.
    frag_lengths : dict[int, int]
        Per-transcript SJ-corrected fragment lengths.  Maps
        ``t_idx → corrected_frag_length`` for each compatible
        transcript with a positive computed size.  Different
        transcripts may yield different sizes because they have
        different introns in the gap between exon blocks.  Empty
        dict when unavailable (e.g. chimeric fragments).
    genomic_footprint : int
        Full genomic span from first exon start to last exon end.
        No intron subtraction — represents the fragment length on
        an unspliced molecule (nRNA or gDNA).  -1 if unavailable.
    genomic_start : int
        Leftmost genomic coordinate of the fragment (first exon start).
        Used to map fragments to transcript-relative coordinates for
        the coverage-weight model.  -1 if unavailable.
    merge_criteria : MergeCriteria
        Which relaxation level produced the transcript sets.
    num_hits : int
        Number of alignment pairs for this read name (NH tag).
    chimera_type : ChimeraType
        Chimera classification (NONE / TRANS /
        CIS_STRAND_SAME / CIS_STRAND_DIFF).
    chimera_gap : int
        Minimum genomic distance between disjoint exon-block clusters
        for intrachromosomal chimeras, -1 for interchromosomal or
        non-chimeric fragments.
    """
    t_inds: frozenset[int]
    n_genes: int
    splice_type: SpliceType
    exon_strand: Strand
    sj_strand: Strand
    frag_lengths: dict[int, int]
    genomic_footprint: int
    genomic_start: int
    merge_criteria: MergeCriteria
    num_hits: int
    overlap_bp: dict[int, tuple[int, int, int]] | None = None
    read_length: int = 0
    nm: int = 0
    chimera_type: ChimeraType = ChimeraType.NONE
    chimera_gap: int = -1

    @property
    def is_chimeric(self) -> bool:
        """True if this fragment is classified as chimeric."""
        return self.chimera_type != ChimeraType.NONE

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
            self.splice_type == SpliceType.SPLICED_ANNOT
            and self.is_unique_gene
            and self.exon_strand in (Strand.POS, Strand.NEG)
            and self.sj_strand in (Strand.POS, Strand.NEG)
            and not self.is_chimeric
        )


# ---------------------------------------------------------------------------
# Fragment resolution — C++ native kernel (required)
# ---------------------------------------------------------------------------

def resolve_fragment(
    frag: Fragment,
    index,
) -> "ResolvedFragment | None":
    """Resolve a fragment to its compatible transcript set.

    Uses the C++ native kernel (``hulkrna._resolve_impl``) via
    ``ResolveContext.resolve_fragment()``.  Returns a C++
    ``ResolvedResult`` object (duck-typed to ``ResolvedFragment``)
    that exposes the same attributes for model training and buffering.
    """
    if not frag.exons:
        return None
    return index._resolve_ctx.resolve_fragment(frag)


# ---------------------------------------------------------------------------
# Multimapper secondary pairing (used by annotate second BAM pass)
# ---------------------------------------------------------------------------

def pair_multimapper_reads(
    sec_r1_locs: list[list],
    sec_r2_locs: list[list],
    index,
    sj_strand_tag,
) -> list[tuple[list, list]]:
    """Pair secondary R1/R2 alignments using transcript-set intersection.

    This function is used by ``annotate.write_annotated_bam()`` to
    reconstruct the same fragment ordering as the C++ BAM scanner
    during the second BAM pass.

    Strategy:
    1. STRICT — pair by R1/R2 transcript-set intersection.
    2. FALLBACK — same-reference closest distance (greedy 1:1).
    3. CROSS-PAIR — remaining unmatched R1 × R2.
    4. SINGLETONS — leftover unpaired locations.
    """
    if not sec_r1_locs and not sec_r2_locs:
        return []

    # --- Resolve each R1 location individually ---
    r1_resolved: list[tuple[list, frozenset, int, int]] = []
    for r1_reads in sec_r1_locs:
        frag = Fragment.from_reads(r1_reads, [], sj_strand_tag=sj_strand_tag)
        t_inds: frozenset = frozenset()
        if frag.exons:
            result = resolve_fragment(frag, index)
            if result is not None:
                t_inds = result.t_inds
        ref_id = r1_reads[0].reference_id if r1_reads else -1
        ref_start = r1_reads[0].reference_start if r1_reads else -1
        r1_resolved.append((r1_reads, t_inds, ref_id, ref_start))

    # --- Resolve each R2 location individually ---
    r2_resolved: list[tuple[list, frozenset, int, int]] = []
    for r2_reads in sec_r2_locs:
        frag = Fragment.from_reads([], r2_reads, sj_strand_tag=sj_strand_tag)
        t_inds = frozenset()
        if frag.exons:
            result = resolve_fragment(frag, index)
            if result is not None:
                t_inds = result.t_inds
        ref_id = r2_reads[0].reference_id if r2_reads else -1
        ref_start = r2_reads[0].reference_start if r2_reads else -1
        r2_resolved.append((r2_reads, t_inds, ref_id, ref_start))

    # --- STRICT — pair by transcript-set intersection ---
    paired: list[tuple[list, list]] = []
    r1_paired: set[int] = set()
    r2_paired: set[int] = set()

    for i, (r1_reads, r1_t, _, _) in enumerate(r1_resolved):
        if not r1_t:
            continue
        for j, (r2_reads, r2_t, _, _) in enumerate(r2_resolved):
            if not r2_t:
                continue
            if r1_t & r2_t:
                paired.append((list(r1_reads), list(r2_reads)))
                r1_paired.add(i)
                r2_paired.add(j)

    # --- FALLBACK — same-reference closest distance ---
    unmatched_r1 = [i for i in range(len(r1_resolved)) if i not in r1_paired]
    unmatched_r2 = [j for j in range(len(r2_resolved)) if j not in r2_paired]

    if unmatched_r1 and unmatched_r2:
        candidates: list[tuple[int, int, int]] = []
        for i in unmatched_r1:
            _, _, r1_ref, r1_pos = r1_resolved[i]
            for j in unmatched_r2:
                _, _, r2_ref, r2_pos = r2_resolved[j]
                if r1_ref == r2_ref and r1_ref >= 0:
                    dist = abs(r1_pos - r2_pos)
                    candidates.append((dist, i, j))
        candidates.sort()
        for _dist, i, j in candidates:
            if i not in r1_paired and j not in r2_paired:
                paired.append((
                    list(r1_resolved[i][0]),
                    list(r2_resolved[j][0]),
                ))
                r1_paired.add(i)
                r2_paired.add(j)

    # --- CROSS-PAIR remaining unmatched R1 × R2 ---
    final_unmatched_r1 = [i for i in range(len(r1_resolved))
                          if i not in r1_paired]
    final_unmatched_r2 = [j for j in range(len(r2_resolved))
                          if j not in r2_paired]

    if final_unmatched_r1 and final_unmatched_r2:
        for i in final_unmatched_r1:
            for j in final_unmatched_r2:
                paired.append((
                    list(r1_resolved[i][0]),
                    list(r2_resolved[j][0]),
                ))
                r1_paired.add(i)
                r2_paired.add(j)

    # --- Remaining singletons ---
    for i in range(len(r1_resolved)):
        if i not in r1_paired:
            paired.append((list(r1_resolved[i][0]), []))
    for j in range(len(r2_resolved)):
        if j not in r2_paired:
            paired.append(([], list(r2_resolved[j][0])))

    return paired
