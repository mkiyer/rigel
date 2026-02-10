"""
hulkrna.count — Two-pass BAM counting pipeline.

Pass 1 (``pass1_learn``) scans a name-sorted BAM, resolves fragments
to compatible transcripts, and trains ``StrandModel`` +
``InsertSizeModel`` from unambiguous fragments.  No intermediate files
are written.

Pass 2 (``pass2_count``) re-scans the same BAM, re-resolves each
fragment identically, and assigns fractional counts into in-memory
``ReadCounter`` arrays using the trained models.

Both passes share ``_resolve_fragment()`` for the core resolution
logic.

Resolution algorithm
--------------------
1. Query each fragment exon block against the unified cgranges index.
   Partition overlaps by ``IntervalType``: EXON, INTRON, INTERGENIC.
2. Query each fragment intron (CIGAR N-op) against the SJ exact-match
   map.  Track annotated vs unannotated status and SJ strand.
3. If the fragment has **any EXON overlap**, merge EXON + SJ transcript
   sets via ``merge_sets_with_criteria()`` for maximum specificity.
4. If **no EXON overlap**, fall back to INTRON overlap sets.
5. If neither EXON nor INTRON overlap exists → **intergenic** (returns
   ``None``).
6. Compute per-transcript insert sizes for compatible transcripts.
"""

import logging
from typing import NamedTuple

import numpy as np
import pandas as pd

from .bam import parse_bam_file
from .core import (
    CountCategory,
    CountStrand,
    CountType,
    IntervalType,
    MergeResult,
    Strand,
    merge_sets_with_criteria,
)
from .fragment import Fragment
from .index import HulkIndex
from .insert_model import InsertSizeModels
from .strand_model import StrandModel

logger = logging.getLogger(__name__)

NUM_COUNT_TYPES = len(CountType)


# ---------------------------------------------------------------------------
# Resolution result
# ---------------------------------------------------------------------------

class _ResolveResult(NamedTuple):
    """Result of resolving a fragment to compatible transcripts/genes."""
    merge_result: MergeResult
    count_cat: CountCategory
    exon_strand: Strand
    sj_strand: Strand


# ---------------------------------------------------------------------------
# Insert size computation
# ---------------------------------------------------------------------------

def _fragment_insert_size(frag: Fragment) -> int:
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


def _compute_insert_size(
    frag: Fragment,
    compatible_t_inds: frozenset[int],
    index: HulkIndex,
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
# Fragment resolution
# ---------------------------------------------------------------------------

def _resolve_fragment(
    frag: Fragment,
    index: HulkIndex,
) -> _ResolveResult | None:
    """Resolve a fragment to its compatible transcript/gene set.

    Queries the index for overlapping intervals and splice junctions,
    merges transcript/gene sets via progressive relaxation, and computes
    per-transcript insert sizes.

    **Chimeric check** must be done by the caller — this function
    assumes all exon blocks are on a single reference.

    Parameters
    ----------
    frag : Fragment
        The consolidated fragment from paired-end reads.
    index : HulkIndex
        The loaded reference index.

    Returns
    -------
    _ResolveResult or None
        ``None`` for intergenic or empty fragments.
    """
    # --- Query each exon block against the unified cgranges index ---
    exon_t_sets: list[frozenset[int]] = []
    exon_g_sets: list[frozenset[int]] = []
    intron_t_sets: list[frozenset[int]] = []
    intron_g_sets: list[frozenset[int]] = []

    exon_strand = Strand.NONE

    for exon_block in frag.exons:
        exon_strand |= exon_block.strand

        block_exon_t: set[int] = set()
        block_exon_g: set[int] = set()
        block_intron_t: set[int] = set()
        block_intron_g: set[int] = set()

        for t_idx, g_idx, itype in index.query_exon(exon_block):
            if itype == IntervalType.EXON:
                block_exon_t.add(t_idx)
                block_exon_g.add(g_idx)
            elif itype == IntervalType.INTRON:
                block_intron_t.add(t_idx)
                block_intron_g.add(g_idx)

        exon_t_sets.append(frozenset(block_exon_t))
        exon_g_sets.append(frozenset(block_exon_g))
        intron_t_sets.append(frozenset(block_intron_t))
        intron_g_sets.append(frozenset(block_intron_g))

    # --- Query SJ matches ---
    sj_t_sets: list[frozenset[int]] = []
    sj_g_sets: list[frozenset[int]] = []
    has_annotated_sj = False
    has_unannotated_sj = False
    sj_strand = Strand.NONE

    for intron in frag.introns:
        key = (intron.ref, intron.start, intron.end, intron.strand)
        match = index.sj_map.get(key)
        if match is not None:
            t_set, g_set = match
            sj_t_sets.append(t_set)
            sj_g_sets.append(g_set)
            has_annotated_sj = True
            sj_strand |= intron.strand
        else:
            has_unannotated_sj = True

    # --- Resolution ---
    any_exon = any(s for s in exon_t_sets)

    if any_exon:
        # Merge EXON + SJ sets together for maximum specificity
        all_t_sets = exon_t_sets + sj_t_sets
        all_g_sets = exon_g_sets + sj_g_sets
        merge_result = merge_sets_with_criteria(all_t_sets, all_g_sets)

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
            merge_result = merge_sets_with_criteria(
                intron_t_sets, intron_g_sets
            )
            count_cat = CountCategory.INTRON
        else:
            # Pure intergenic — no compatible hits
            return None

    if merge_result.is_empty:
        # Defensive — shouldn't happen with valid input
        return None

    return _ResolveResult(
        merge_result=merge_result,
        count_cat=count_cat,
        exon_strand=exon_strand,
        sj_strand=sj_strand,
    )


# ---------------------------------------------------------------------------
# ReadCounter
# ---------------------------------------------------------------------------

class ReadCounter:
    """Accumulates fractional read counts into transcript and gene arrays.

    Count arrays have shape ``(N, 12)`` where the 12 columns correspond
    to ``CountType`` values (4 categories × 3 strands).
    """

    def __init__(self, num_transcripts: int, num_genes: int):
        self.t_counts = np.zeros(
            (num_transcripts, NUM_COUNT_TYPES), dtype=np.float32
        )
        self.g_counts = np.zeros(
            (num_genes, NUM_COUNT_TYPES), dtype=np.float32
        )

    def assign(
        self,
        result: _ResolveResult,
        index: HulkIndex,
        num_hits: int = 1,
    ) -> None:
        """Assign counts for a resolved fragment.

        Determines sense/antisense strand by comparing the fragment's
        exon alignment strand to the reference annotation strand.

        - **Unique gene**: SENSE if fragment aligns in same direction
          as the gene, ANTISENSE otherwise.
        - **Multi-gene**: AMBIGUOUS (fractional 1/N assignment).

        For multimappers (``num_hits > 1``), each alignment pair's
        contribution is weighted by ``1 / num_hits``.

        Parameters
        ----------
        result : _ResolveResult
            Resolution output from ``_resolve_fragment()``.
        index : HulkIndex
            The loaded reference index.
        num_hits : int
            Number of alignment pairs for this read name (NH tag).
        """
        mr = result.merge_result
        t_inds = mr.t_inds
        g_inds = mr.g_inds

        if not t_inds:
            return

        # --- Determine strand ---
        if mr.is_unique_gene:
            g_idx = next(iter(g_inds))
            ref_strand = int(index.g_to_strand_arr[g_idx])
            frag_strand = int(result.exon_strand)
            if frag_strand in (Strand.POS, Strand.NEG):
                if frag_strand == ref_strand:
                    count_strand = CountStrand.SENSE
                else:
                    count_strand = CountStrand.ANTISENSE
            else:
                count_strand = CountStrand.AMBIGUOUS
        else:
            count_strand = CountStrand.AMBIGUOUS

        count_type = int(result.count_cat) * 3 + int(count_strand)

        # --- Assign transcript counts ---
        t_weight = 1.0 / (len(t_inds) * num_hits)
        for t_idx in t_inds:
            self.t_counts[t_idx, count_type] += t_weight

        # --- Assign gene counts ---
        g_weight = 1.0 / (len(g_inds) * num_hits)
        for g_idx in g_inds:
            self.g_counts[g_idx, count_type] += g_weight

    def get_t_counts_df(self) -> pd.DataFrame:
        """Return transcript counts as a DataFrame."""
        return pd.DataFrame(
            self.t_counts, columns=list(CountType.columns())
        )

    def get_g_counts_df(self) -> pd.DataFrame:
        """Return gene counts as a DataFrame."""
        return pd.DataFrame(
            self.g_counts, columns=list(CountType.columns())
        )


# ---------------------------------------------------------------------------
# Pass 1: learn models
# ---------------------------------------------------------------------------

def pass1_learn(
    bam_iter,
    index: HulkIndex,
    *,
    skip_duplicates: bool = True,
    log_every: int = 1_000_000,
    insert_size_max: int = 1000,
) -> tuple[dict, StrandModel, InsertSizeModels]:
    """Scan a BAM file, resolve fragments, and train models.

    Iterates over read-pair groups from a name-sorted / collated BAM,
    builds ``Fragment`` objects, resolves each to its compatible
    transcript/gene set, and trains ``StrandModel`` and
    ``InsertSizeModel`` from unambiguous fragments.

    Multimapping reads (NH > 1) are always excluded from Pass 1 since
    they do not contribute useful information for model training.

    No intermediate files are written.

    Parameters
    ----------
    bam_iter
        Iterator over ``pysam.AlignedSegment`` objects.
    index : HulkIndex
        The loaded reference index.
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    log_every : int
        Log progress every *log_every* fragments.
    insert_size_max : int
        Maximum insert size for the histogram model.

    Returns
    -------
    tuple[dict, StrandModel, InsertSizeModels]
        ``(stats, strand_model, insert_models)`` where *stats* is a
        flat dict with BAM parsing stats and Pass 1 resolution stats.
    """
    # BAM parsing stats (populated by parse_bam_file)
    stats: dict = {}

    # Models
    strand_model = StrandModel()
    insert_models = InsertSizeModels(max_size=insert_size_max)

    # Counters
    frag_id = 0
    n_chimeric = 0
    n_intergenic = 0
    n_with_exon = 0
    n_with_intron_fallback = 0
    n_with_annotated_sj = 0
    n_with_unannotated_sj = 0
    n_unique_gene = 0
    n_multi_gene = 0
    n_strand_skipped_no_sj = 0
    n_strand_skipped_multi_gene = 0
    n_strand_skipped_ambiguous = 0
    n_insert_unambiguous = 0
    n_insert_ambiguous = 0
    n_insert_intergenic = 0

    logger.info("[START] Pass 1: scanning BAM → learning models")

    for pairs in parse_bam_file(
        bam_iter, stats,
        skip_duplicates=skip_duplicates,
        include_multimap=False,  # Always exclude multimappers in Pass 1
    ):
        for r1, r2 in pairs:
            frag = Fragment.from_reads(r1, r2)

            # Check chimeric (exon blocks on multiple refs)
            refs = set(e.ref for e in frag.exons)
            if len(refs) > 1:
                n_chimeric += 1
                continue
            if not frag.exons:
                continue

            result = _resolve_fragment(frag, index)

            if result is None:
                n_intergenic += 1
                # Compute simple insert size for intergenic model
                insert_size = _fragment_insert_size(frag)
                if insert_size > 0:
                    insert_models.observe(insert_size, count_cat=None)
                    n_insert_intergenic += 1
                continue

            # --- Stats: overlap type ---
            if result.count_cat == CountCategory.INTRON:
                n_with_intron_fallback += 1
            else:
                n_with_exon += 1
            if result.count_cat == CountCategory.SPLICED_ANNOT:
                n_with_annotated_sj += 1
            elif result.count_cat == CountCategory.SPLICED_UNANNOT:
                n_with_unannotated_sj += 1

            # --- Stats: gene ambiguity ---
            if result.merge_result.is_unique_gene:
                n_unique_gene += 1
            else:
                n_multi_gene += 1

            # --- Strand model training ---
            if result.count_cat != CountCategory.SPLICED_ANNOT:
                n_strand_skipped_no_sj += 1
            elif not result.merge_result.is_unique_gene:
                n_strand_skipped_multi_gene += 1
            elif (
                result.exon_strand not in (Strand.POS, Strand.NEG)
                or result.sj_strand not in (Strand.POS, Strand.NEG)
            ):
                n_strand_skipped_ambiguous += 1
            else:
                strand_model.observe(
                    result.exon_strand, result.sj_strand
                )

            # --- Insert size model training ---
            insert_size = _compute_insert_size(
                frag, result.merge_result.t_inds, index
            )
            if insert_size > 0:
                insert_models.observe(insert_size, result.count_cat)
                n_insert_unambiguous += 1
            else:
                n_insert_ambiguous += 1

        frag_id += 1

        if frag_id % log_every == 0:
            logger.debug(
                f"  processed {frag_id:,} fragments, "
                f"{strand_model.n_observations:,} strand obs, "
                f"{insert_models.n_observations:,} insert obs"
            )

    logger.info(
        f"[DONE] Pass 1: {frag_id:,} fragments → "
        f"{n_unique_gene:,} unique-gene, "
        f"{n_multi_gene:,} multi-gene, "
        f"{n_intergenic:,} intergenic"
    )
    strand_model.log_summary()
    insert_models.log_summary()

    # Merge resolution stats into the BAM stats dict
    stats["n_fragments"] = frag_id
    stats["n_chimeric"] = n_chimeric
    stats["n_intergenic"] = n_intergenic
    stats["n_with_exon"] = n_with_exon
    stats["n_with_intron_fallback"] = n_with_intron_fallback
    stats["n_with_annotated_sj"] = n_with_annotated_sj
    stats["n_with_unannotated_sj"] = n_with_unannotated_sj
    stats["n_unique_gene"] = n_unique_gene
    stats["n_multi_gene"] = n_multi_gene
    stats["n_strand_skipped_no_sj"] = n_strand_skipped_no_sj
    stats["n_strand_skipped_multi_gene"] = n_strand_skipped_multi_gene
    stats["n_strand_skipped_ambiguous"] = n_strand_skipped_ambiguous
    stats["n_insert_unambiguous"] = n_insert_unambiguous
    stats["n_insert_ambiguous"] = n_insert_ambiguous
    stats["n_insert_intergenic"] = n_insert_intergenic

    return stats, strand_model, insert_models


# ---------------------------------------------------------------------------
# Pass 2: assign counts
# ---------------------------------------------------------------------------

def pass2_count(
    bam_iter,
    index: HulkIndex,
    strand_model: StrandModel,
    insert_models: InsertSizeModels,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    log_every: int = 1_000_000,
) -> tuple[dict, ReadCounter]:
    """Re-scan BAM and assign fractional counts using trained models.

    Repeats the same fragment resolution as Pass 1, then assigns
    counts to transcript/gene arrays via ``ReadCounter``.

    The ``strand_model`` and ``insert_model`` parameters are accepted
    for forward compatibility with future Bayesian counting but are
    not used in the current fractional (1/N) placeholder.

    Parameters
    ----------
    bam_iter
        Iterator over ``pysam.AlignedSegment`` objects.
    index : HulkIndex
        The loaded reference index.
    strand_model : StrandModel
        Trained strand model from Pass 1 (reserved for future use).
    insert_models : InsertSizeModels
        Trained insert size models from Pass 1 (reserved for future use).
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include multimapping reads (NH > 1, default False).
    log_every : int
        Log progress every *log_every* fragments.

    Returns
    -------
    tuple[dict, ReadCounter]
        ``(stats, counter)`` where *stats* is a flat dict with BAM
        parsing stats and resolution stats, and *counter* holds the
        accumulated transcript/gene count arrays.
    """
    # BAM parsing stats
    stats: dict = {}

    # Count arrays
    counter = ReadCounter(index.num_transcripts, index.num_genes)

    # Counters
    frag_id = 0
    n_chimeric = 0
    n_intergenic = 0
    n_with_exon = 0
    n_with_intron_fallback = 0
    n_unique_gene = 0
    n_multi_gene = 0

    logger.info("[START] Pass 2: scanning BAM → assigning counts")

    for pairs in parse_bam_file(
        bam_iter, stats,
        skip_duplicates=skip_duplicates,
        include_multimap=include_multimap,
    ):
        num_hits = len(pairs)
        for r1, r2 in pairs:
            frag = Fragment.from_reads(r1, r2)

            # Check chimeric (exon blocks on multiple refs)
            refs = set(e.ref for e in frag.exons)
            if len(refs) > 1:
                n_chimeric += 1
                continue
            if not frag.exons:
                continue

            result = _resolve_fragment(frag, index)

            if result is None:
                n_intergenic += 1
                continue

            # --- Stats ---
            if result.count_cat == CountCategory.INTRON:
                n_with_intron_fallback += 1
            else:
                n_with_exon += 1
            if result.merge_result.is_unique_gene:
                n_unique_gene += 1
            else:
                n_multi_gene += 1

            # --- Assign counts ---
            counter.assign(result, index, num_hits=num_hits)

        frag_id += 1

        if frag_id % log_every == 0:
            logger.debug(f"  processed {frag_id:,} fragments")

    logger.info(
        f"[DONE] Pass 2: {frag_id:,} fragments → "
        f"{n_unique_gene:,} unique-gene, "
        f"{n_multi_gene:,} multi-gene, "
        f"{n_intergenic:,} intergenic"
    )

    # Populate stats
    stats["n_fragments"] = frag_id
    stats["n_chimeric"] = n_chimeric
    stats["n_intergenic"] = n_intergenic
    stats["n_with_exon"] = n_with_exon
    stats["n_with_intron_fallback"] = n_with_intron_fallback
    stats["n_unique_gene"] = n_unique_gene
    stats["n_multi_gene"] = n_multi_gene

    return stats, counter
