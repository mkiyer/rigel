"""
hulkrna.resolution — Fragment resolution and fragment-length computation.

Resolves aligned fragments to their compatible transcript/gene sets
by querying the reference index and applying progressive set merging.
Also computes fragment lengths with gap correction.

This module contains:
- ``merge_sets_with_criteria()`` — progressive set merging
- ``resolve_fragment()`` — core fragment-to-transcript resolution with
  chimera detection (interchromosomal and intrachromosomal)
- ``compute_frag_lengths()`` — per-transcript fragment lengths with gap correction
- ``ResolvedFragment`` — compact cached resolution result
- ``pair_multimapper_reads()`` — secondary alignment pairing
"""

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


from dataclasses import dataclass

from .types import (
    EMPTY_MERGE,
    ChimeraType,
    IntervalType,
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
# Fragment resolution — dispatch (C++ native or Python fallback)
# ---------------------------------------------------------------------------

def resolve_fragment(
    frag: Fragment,
    index,
) -> "ResolvedFragment | None":
    """Resolve a fragment to its compatible transcript set.

    Dispatches to the C++ native kernel (``hulkrna._resolve_impl``)
    when available, otherwise falls back to the pure-Python
    implementation.  Both paths produce identical results.

    When the native kernel is available, returns a C++
    ``ResolvedResult`` object (duck-typed to ``ResolvedFragment``)
    that exposes the same attributes for model training and buffering.
    """
    if not frag.exons:
        return None

    ctx = getattr(index, '_resolve_ctx', None)
    if ctx is not None:
        return ctx.resolve_fragment(frag)
    return _resolve_python(frag, index)


def _resolve_native(
    frag: Fragment,
    index,
    ctx,
) -> ResolvedFragment | None:
    """Resolve via C++ kernel — called from ``resolve_fragment``."""
    rtid = index._resolve_ref_to_id
    exon_ref_ids = [rtid.get(e.ref, -1) for e in frag.exons]
    exon_starts = [e.start for e in frag.exons]
    exon_ends = [e.end for e in frag.exons]
    exon_strands = [int(e.strand) for e in frag.exons]

    intron_ref_ids = [rtid.get(i.ref, -1) for i in frag.introns]
    intron_starts = [i.start for i in frag.introns]
    intron_ends = [i.end for i in frag.introns]
    intron_strands = [int(i.strand) for i in frag.introns]

    result = ctx.resolve(
        exon_ref_ids, exon_starts, exon_ends, exon_strands,
        intron_ref_ids, intron_starts, intron_ends, intron_strands,
        frag.genomic_footprint,
    )
    if result is None:
        return None

    return ResolvedFragment(
        t_inds=frozenset(result[0]),
        n_genes=result[1],
        splice_type=SpliceType(result[2]),
        exon_strand=Strand(result[3]),
        sj_strand=Strand(result[4]),
        frag_lengths=result[5],
        genomic_footprint=result[6],
        genomic_start=result[7],
        merge_criteria=MergeCriteria(result[8]),
        num_hits=1,
        overlap_bp=result[9],
        read_length=result[10],
        chimera_type=ChimeraType(result[11]),
        chimera_gap=result[12],
    )


# ---------------------------------------------------------------------------
# Fragment resolution — Python reference implementation
# ---------------------------------------------------------------------------

def _resolve_python(
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
    ``CIS_STRAND_SAME`` or ``CIS_STRAND_DIFF`` based on the alignment strand
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
    chimera_type = ChimeraType.NONE
    chimera_gap = -1

    refs = set(e.ref for e in frag.exons)
    is_interchromosomal = len(refs) > 1
    if is_interchromosomal:
        chimera_type = ChimeraType.TRANS

    # --- Query each exon block + compute overlap profile (single pass) ---
    # Uses the collapsed index.query() once per block.  Each hit returns
    # (h_start, h_end, itype, t_set) where t_set is already a frozenset.
    exon_t_sets: list[frozenset[int]] = []
    transcript_t_sets: list[frozenset[int]] = []

    exon_strand = Strand.NONE

    # Overlap profile accumulators (A2: array-based, indexed by t_idx)
    n_tx = len(index.t_to_g_arr)
    _t_exon_bp = [0] * n_tx
    _t_transcript_bp = [0] * n_tx
    _t_unambig_intron_bp = [0] * n_tx
    read_length = 0

    for exon_block in frag.exons:
        exon_strand |= exon_block.strand
        block_start = exon_block.start
        block_end = exon_block.end
        read_length += block_end - block_start

        block_exon_t: set[int] = set()
        block_transcript_t: set[int] = set()

        # Single cgranges query per block (collapsed)
        for h_start, h_end, itype, t_set in index.query(exon_block):
            if itype == IntervalType.EXON:
                block_exon_t.update(t_set)
                clipped_lo = max(block_start, h_start)
                clipped_hi = min(block_end, h_end)
                if clipped_hi > clipped_lo:
                    bp = clipped_hi - clipped_lo
                    for t_idx in t_set:
                        _t_exon_bp[t_idx] += bp
            elif itype == IntervalType.TRANSCRIPT:
                block_transcript_t.update(t_set)
                clipped_lo = max(block_start, h_start)
                clipped_hi = min(block_end, h_end)
                if clipped_hi > clipped_lo:
                    bp = clipped_hi - clipped_lo
                    for t_idx in t_set:
                        _t_transcript_bp[t_idx] += bp
            elif itype == IntervalType.UNAMBIG_INTRON:
                clipped_lo = max(block_start, h_start)
                clipped_hi = min(block_end, h_end)
                if clipped_hi > clipped_lo:
                    bp = clipped_hi - clipped_lo
                    for t_idx in t_set:
                        _t_unambig_intron_bp[t_idx] += bp

        exon_t_sets.append(frozenset(block_exon_t))
        transcript_t_sets.append(frozenset(block_transcript_t))

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
        if match is None and intron.strand not in (Strand.POS, Strand.NEG):
            # Some aligners/outputs omit splice-junction strand tags.
            # Fall back to strand-agnostic coordinate matching so that
            # annotated SJs still contribute transcript evidence.
            pos_match = index.sj_map.get(
                (intron.ref, intron.start, intron.end, Strand.POS)
            )
            neg_match = index.sj_map.get(
                (intron.ref, intron.start, intron.end, Strand.NEG)
            )
            if pos_match is not None or neg_match is not None:
                t_set: set[int] = set()
                if pos_match is not None:
                    t_set.update(pos_match)
                if neg_match is not None:
                    t_set.update(neg_match)
                match = frozenset(t_set)
        if match is not None:
            sj_t_sets.append(match)
            has_annotated_sj = True
            sj_strand |= intron.strand
        else:
            has_unannotated_sj = True

    # --- Resolution: Three-state categorization ---
    # Priority order:
    # 1. Annotated SJs → SPLICED_ANNOT (definitively mature RNA)
    # 2. Unannotated SJs → SPLICED_UNANNOT
    # 3. No SJs → UNSPLICED (overlap profile captures exonic/intronic)
    #
    # Intronic vs exonic overlap is NOT used for category assignment.
    # The per-candidate overlap profile (n_exon_bp, n_intron_bp) handles
    # the distinction between mature RNA, nascent RNA, and gDNA in the
    # 3-pool likelihood scoring.

    any_exon = any(s for s in exon_t_sets)
    any_transcript = any(s for s in transcript_t_sets)

    if has_annotated_sj and any_exon:
        # --- SPLICED_ANNOT: annotated SJs prove mature RNA origin ---
        exon_merge = merge_sets_with_criteria(exon_t_sets)
        sj_merge = merge_sets_with_criteria(sj_t_sets)
        if not sj_merge.is_empty:
            if not exon_merge.is_empty:
                t_inds = exon_merge.t_inds & sj_merge.t_inds
                if not t_inds:
                    t_inds = sj_merge.t_inds
            else:
                t_inds = sj_merge.t_inds
            merge_result = MergeResult(t_inds, sj_merge.criteria)
        else:
            all_t_sets = exon_t_sets + sj_t_sets
            merge_result = merge_sets_with_criteria(all_t_sets)
        splice_type = SpliceType.SPLICED_ANNOT

    elif any_exon or any_transcript:
        # --- Genic fragment (exon and/or transcript-span overlap) ---
        # Merge exon and transcript-span transcript sets INDEPENDENTLY,
        # then union the results.  This prevents transcript-span overlap
        # from *narrowing* the exon-based candidate set — critical for
        # overlapping antisense genes where gene A's intronic span
        # spatially overlaps gene B's exon.  The per-candidate overlap
        # profile (exon_bp, intron_bp) still carries the intronic
        # information downstream for probabilistic scoring.
        parts: list[frozenset[int]] = []
        best_criteria = MergeCriteria.UNION

        if any_exon:
            exon_merge = merge_sets_with_criteria(exon_t_sets)
            if not exon_merge.is_empty:
                parts.append(exon_merge.t_inds)
                best_criteria = exon_merge.criteria

        if any_transcript:
            tx_merge = merge_sets_with_criteria(transcript_t_sets)
            if not tx_merge.is_empty:
                parts.append(tx_merge.t_inds)

        if parts:
            t_inds = frozenset.union(*parts)
            merge_result = MergeResult(
                t_inds,
                best_criteria if len(parts) == 1
                else MergeCriteria.UNION,
            )
        else:
            merge_result = EMPTY_MERGE

        if has_unannotated_sj:
            splice_type = SpliceType.SPLICED_UNANNOT
        else:
            splice_type = SpliceType.UNSPLICED

    else:
        # Pure intergenic — no compatible hits
        return None

    if merge_result.is_empty:
        # Defensive — shouldn't happen with valid input
        return None

    # --- Derive final transcript set ---
    t_inds = merge_result.t_inds
    if not t_inds:
        return None

    # --- Overlap profiles for all genic fragments ---
    # Overlap BP counts were computed in the single-pass query above.
    # intron_bp is derived as transcript_bp - exon_bp (not stored directly).
    # unambig_intron_bp comes from UNAMBIG_INTRON cgranges hits.
    all_overlap_t: set[int] = set()
    for _ts in exon_t_sets:
        all_overlap_t.update(_ts)
    for _ts in transcript_t_sets:
        all_overlap_t.update(_ts)
    if all_overlap_t and read_length > 0:
        overlap_bp = {
            t_idx: (
                _t_exon_bp[t_idx],
                max(_t_transcript_bp[t_idx] - _t_exon_bp[t_idx], 0),
                _t_unambig_intron_bp[t_idx],
            )
            for t_idx in all_overlap_t
        }
    else:
        overlap_bp = {}
    if not overlap_bp:
        # Pure intronic fragments have no exonic overlap; use 0 exon_bp.
        # intron_bp = read_length to indicate fully genic.
        # unambig_intron_bp = read_length (all intronic, no exonic overlap).
        overlap_bp = {
            t_idx: (0, read_length, read_length) for t_idx in t_inds
        }

    # --- Derive n_genes from transcript indices ---
    n_genes = len(set(int(index.t_to_g_arr[t]) for t in t_inds))

    # --- Compute fragment lengths (skip for chimeras — meaningless) ---
    # frag_lengths: per-transcript SJ-corrected
    if chimera_type == ChimeraType.NONE:
        frag_lengths = compute_frag_lengths(frag, t_inds, index)
    else:
        frag_lengths = {}

    return ResolvedFragment(
        t_inds=t_inds,
        n_genes=n_genes,
        splice_type=splice_type,
        exon_strand=exon_strand,
        sj_strand=sj_strand,
        frag_lengths=frag_lengths,
        genomic_footprint=frag.genomic_footprint,
        genomic_start=min(e.start for e in frag.exons) if frag.exons else -1,
        merge_criteria=merge_result.criteria,
        num_hits=1,  # Caller must set this from the pairs list length
        overlap_bp=overlap_bp,
        read_length=read_length,
        chimera_type=chimera_type,
        chimera_gap=chimera_gap,
    )


# ---------------------------------------------------------------------------
# Multimapper secondary pairing (resolve-then-pair)
# ---------------------------------------------------------------------------


def pair_multimapper_reads(
    sec_r1_locs: list[list],
    sec_r2_locs: list[list],
    index,
    sj_strand_tag,
) -> list[tuple[list, list]]:
    """Pair secondary R1/R2 alignments using transcript-set intersection.

    Instead of blind cross-product pairing (which creates false chimeric
    pairs between paralogs on different chromosomes), each secondary R1
    and R2 location is resolved independently to obtain its compatible
    transcript set.  Pairs are then formed using a multi-tier strategy:

    1. **STRICT** — create all pairs where the R1 and R2 transcript
       sets have a non-empty intersection.  These are the most likely
       true pairs since both mates agree on at least one transcript.

    2. **FALLBACK** — for unmatched locations, try same-reference
       pairing by genomic proximity (closest distance, greedy 1:1).
       This handles intergenic secondaries that fall near each other.

    3. **CROSS-PAIR** — any remaining unmatched R1 and R2 locations
       are cross-paired (all remaining R1 × all remaining R2).  This
       preserves interchromosomal chimera evidence — these pairs
       enter the EM as ambiguous multimappers and are scored
       alongside real transcript candidates.

    4. **SINGLETONS** — if only unmatched R1 or only unmatched R2
       locations remain (no cross-pair partner), they become
       singletons.

    No hard distance filter is applied.  Implausible pairs (e.g.
    artificial "Frankenstein fragments" spanning distant paralogs)
    receive per-candidate fragment-length penalties in the EM, which
    naturally drives their posterior weight toward zero.

    Parameters
    ----------
    sec_r1_locs : list[list]
        Secondary R1 locations.  Each element is a list containing
        one ``pysam.AlignedSegment``.
    sec_r2_locs : list[list]
        Secondary R2 locations.
    index : HulkIndex
        Reference index for resolution.
    sj_strand_tag : str or tuple[str, ...]
        Splice-junction strand tag(s).

    Returns
    -------
    list[tuple[list, list]]
        Paired ``(r1_reads, r2_reads)`` tuples for secondary hits.
    """
    if not sec_r1_locs and not sec_r2_locs:
        return []

    # --- Step 1: Resolve each R1 location individually ---
    # Store (reads, t_inds, ref_id, ref_start)
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

    # --- Step 2: Resolve each R2 location individually ---
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

    # --- Step 3: STRICT — pair by transcript-set intersection ---
    paired: list[tuple[list, list]] = []
    r1_paired: set[int] = set()
    r2_paired: set[int] = set()

    for i, (r1_reads, r1_t, _, _) in enumerate(r1_resolved):
        if not r1_t:
            continue  # intergenic R1, skip strict
        for j, (r2_reads, r2_t, _, _) in enumerate(r2_resolved):
            if not r2_t:
                continue  # intergenic R2, skip strict
            if r1_t & r2_t:  # non-empty intersection
                paired.append((list(r1_reads), list(r2_reads)))
                r1_paired.add(i)
                r2_paired.add(j)

    # --- Step 4: FALLBACK — same-reference closest distance ---
    unmatched_r1 = [i for i in range(len(r1_resolved)) if i not in r1_paired]
    unmatched_r2 = [j for j in range(len(r2_resolved)) if j not in r2_paired]

    if unmatched_r1 and unmatched_r2:
        # Build all same-reference (i, j, distance) candidates
        candidates: list[tuple[int, int, int]] = []
        for i in unmatched_r1:
            _, _, r1_ref, r1_pos = r1_resolved[i]
            for j in unmatched_r2:
                _, _, r2_ref, r2_pos = r2_resolved[j]
                if r1_ref == r2_ref and r1_ref >= 0:
                    dist = abs(r1_pos - r2_pos)
                    candidates.append((dist, i, j))

        # Greedy 1:1 matching by shortest distance
        candidates.sort()
        for _dist, i, j in candidates:
            if i not in r1_paired and j not in r2_paired:
                paired.append((
                    list(r1_resolved[i][0]),
                    list(r2_resolved[j][0]),
                ))
                r1_paired.add(i)
                r2_paired.add(j)

    # --- Step 5: CROSS-PAIR remaining unmatched R1 × R2 ---
    # Preserves interchromosomal chimera evidence.  These pairs
    # enter the EM as ambiguous multimappers; implausible ones
    # receive heavy fragment-length penalties and get near-zero weight.
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

    # --- Step 6: Remaining singletons (only R1 or only R2 left) ---
    for i in range(len(r1_resolved)):
        if i not in r1_paired:
            paired.append((list(r1_resolved[i][0]), []))
    for j in range(len(r2_resolved)):
        if j not in r2_paired:
            paired.append(([], list(r2_resolved[j][0])))

    return paired
