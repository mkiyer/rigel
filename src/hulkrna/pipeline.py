"""
hulkrna.pipeline — Single-pass BAM counting pipeline with EM.

Architecture
------------
The pipeline reads the BAM file **once** and processes counts in
two logical stages:

**BAM Scan**: Parse all fragments, resolve each against the
reference index, train strand and insert-size models from
unique-mapper fragments, and buffer all resolved fragments into a
memory-efficient columnar buffer (``FragmentBuffer``).

**Count Assignment**: Iterate the buffer once to:

1. **Unique assignment** — deterministic count += 1 for unambiguous
   SPLICED_ANNOT fragments (annotated splice junctions definitively
   prove RNA origin).
2. **Pre-compute EM data** — for every other fragment, compute the
   data likelihood for each candidate transcript *and* a single gDNA
   pseudo-component, then store in CSR arrays.
3. **EM iteration** — run vectorized Expectation-Maximization over
   *all* ambiguous fragments simultaneously to estimate transcript
   *and gDNA* abundances.  The gDNA component competes with
   transcripts for probability mass.
4. **Stochastic assignment** — use the converged EM posteriors to
   make a single discrete stochastic choice per fragment (transcript
   or gDNA).

Each fragment contributes exactly **one discrete count** (1.0) to a
single component (transcript or gDNA).  The EM approach follows the
standard algorithm used by RSEM, Salmon, and Kallisto, extended
with a genomic DNA contamination model.

The columnar buffer stores fragment data at ~40-50 bytes per fragment
(vs. ~600 bytes for Python ``ResolvedFragment`` objects).  When
memory exceeds a configurable threshold, chunks are spilled to disk
as Arrow IPC (Feather v2) files with LZ4 compression.
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bam import parse_bam_file, detect_sj_strand_tag
from .buffer import (
    FragmentBuffer,
    FRAG_UNIQUE,
    FRAG_ISOFORM_AMBIG,
    FRAG_GENE_AMBIG,
    FRAG_MULTIMAPPER,
    FRAG_CHIMERIC,
)
from .categories import CountCategory
from .counter import ReadCounter, EMData
from .fragment import Fragment
from .types import ChimeraType, Strand
from .index import HulkIndex
from .insert_model import InsertSizeModels
from .resolution import (
    fragment_insert_size,
    resolve_fragment,
)
from .stats import PipelineStats
from .strand_model import StrandModel, StrandModels

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline result
# ---------------------------------------------------------------------------

@dataclass
class PipelineResult:
    """Complete output of the three-phase pipeline."""

    stats: PipelineStats
    strand_models: StrandModels
    insert_models: InsertSizeModels
    counter: ReadCounter

    @property
    def strand_model(self) -> StrandModel:
        """Backward-compatible alias for the exonic strand model."""
        return self.strand_models.exonic


# ---------------------------------------------------------------------------
# BAM Scan: resolve, train models, buffer
# ---------------------------------------------------------------------------

def scan_and_buffer(
    bam_iter,
    index: HulkIndex,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    log_every: int = 1_000_000,
    insert_size_max: int = 1000,
    chunk_size: int = 1_000_000,
    max_memory_bytes: int = 2 * 1024**3,
    spill_dir: Path | None = None,
    sj_strand_tag: str | tuple[str, ...] = "XS",
) -> tuple[PipelineStats, StrandModels, InsertSizeModels, FragmentBuffer]:
    """Single-pass BAM scan: resolve fragments, train models, buffer results.

    Reads the BAM once.  For each fragment:

    1. Build a ``Fragment`` from the read pair.
    2. Resolve against the reference index.
    3. Train strand and insert-size models (unique mappers only).
    4. Buffer all resolved fragments into the columnar buffer.
       Each fragment receives a ``frag_id`` that groups alignments
       of the same molecule for Phase 4 multimapper handling.

    Three strand models are trained:

    - **exonic**: from SPLICED_ANNOT fragments with unique gene and
      unambiguous strands (same as the original single model).
    - **intronic**: from INTRON fragments with unique gene assignment,
      using the exon strand vs gene strand relationship.
    - **intergenic**: from intergenic fragments (trained using alignment
      strand when splice junction strand is available from nearby context).

    Models are trained exclusively from unique-mapper fragments
    (NH = 1) regardless of the *include_multimap* setting.

    Parameters
    ----------
    bam_iter
        Iterator over ``pysam.AlignedSegment`` objects.
    index : HulkIndex
        The loaded reference index.
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include multimapping reads (NH > 1) for Phase 4 counting
        (default False).  Model training always excludes
        multimappers.
    log_every : int
        Log progress every *log_every* read-name groups.
    insert_size_max : int
        Maximum insert size for the histogram model.
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory for in-memory buffer chunks before disk spill
        (default 2 GiB).  Set to 0 to disable spilling.
    spill_dir : Path or None
        Directory for spilled chunk files (default: system temp dir).
    sj_strand_tag : str or tuple of str
        BAM tag(s) for splice-junction strand (``"XS"`` for STAR,
        ``"ts"`` for minimap2; default ``"XS"``).  A tuple is
        checked in order; the first present tag on each read wins.

    Returns
    -------
    tuple[PipelineStats, StrandModels, InsertSizeModels, FragmentBuffer]
        All resolved fragments (unique and multimapper) are stored in
        the buffer with ``frag_id`` for grouping multimapper hits.
    """
    stats = PipelineStats()
    bam_stats = stats.as_bam_stats_dict()

    strand_models = StrandModels()
    insert_models = InsertSizeModels(max_size=insert_size_max)
    buffer = FragmentBuffer(
        t_to_g_arr=index.t_to_g_arr,
        chunk_size=chunk_size,
        max_memory_bytes=max_memory_bytes,
        spill_dir=spill_dir,
    )

    logger.info("[START] Scanning BAM → resolve + train models + buffer")

    group_id = 0
    for nh, hits in parse_bam_file(
        bam_iter,
        bam_stats,
        skip_duplicates=skip_duplicates,
        include_multimap=include_multimap,
    ):
        num_hits = nh
        is_unique_mapper = num_hits == 1
        n_buffered_mm = 0

        for r1_reads, r2_reads in hits:
            frag = Fragment.from_reads(r1_reads, r2_reads, sj_strand_tag=sj_strand_tag)
            stats.n_fragments += 1

            if not frag.exons:
                continue

            result = resolve_fragment(frag, index)

            if result is None:
                # --- Intergenic: deterministic gDNA assignment ---
                if frag.introns:
                    stats.n_intergenic_spliced += 1
                else:
                    stats.n_intergenic_unspliced += 1
                # Train intergenic insert model (unique mappers only)
                if is_unique_mapper:
                    isize = fragment_insert_size(frag)
                    if isize > 0:
                        insert_models.observe(isize, count_cat=None)
                        stats.n_insert_intergenic += 1
                continue

            # Set num_hits from pairs list length
            result.num_hits = num_hits

            # --- Handle chimeric fragments ---
            if result.is_chimeric:
                ct = result.chimera_type
                stats.n_chimeric += 1
                if ct == ChimeraType.INTERCHROMOSOMAL:
                    stats.n_chimeric_interchrom += 1
                elif ct == ChimeraType.STRAND_SAME:
                    stats.n_chimeric_strand_same += 1
                elif ct == ChimeraType.STRAND_DIFF:
                    stats.n_chimeric_strand_diff += 1
                # Buffer chimeric fragments for reporting, but skip
                # model training since they are not countable.
                buffer.append(result, frag_id=group_id)
                continue

            # --- Stats: overlap type ---
            if result.count_cat == CountCategory.INTRON:
                stats.n_with_intron_fallback += 1
            else:
                stats.n_with_exon += 1
            if result.count_cat == CountCategory.SPLICED_ANNOT:
                stats.n_with_annotated_sj += 1
            elif result.count_cat == CountCategory.SPLICED_UNANNOT:
                stats.n_with_unannotated_sj += 1

            # --- Stats: gene ambiguity ---
            if result.is_unique_gene:
                stats.n_unique_gene += 1
            else:
                stats.n_multi_gene += 1

            if is_unique_mapper:
                # --- Train models (unique mappers only) ---
                # Exonic strand model (from SPLICED_ANNOT with known SJ strand)
                if result.is_strand_qualified:
                    strand_models.exonic.observe(result.exon_strand, result.sj_strand)
                    stats.n_strand_trained += 1
                elif result.count_cat != CountCategory.SPLICED_ANNOT:
                    stats.n_strand_skipped_no_sj += 1
                elif not result.is_unique_gene:
                    stats.n_strand_skipped_multi_gene += 1
                else:
                    stats.n_strand_skipped_ambiguous += 1

                # Intronic strand model: INTRON + unique gene + known strands
                if (
                    result.count_cat == CountCategory.INTRON
                    and result.is_unique_gene
                    and result.exon_strand in (Strand.POS, Strand.NEG)
                ):
                    # Use gene strand as the "truth" for intronic fragments
                    t_idx = int(next(iter(result.t_inds)))
                    g_idx = int(index.t_to_g_arr[t_idx])
                    g_strand = Strand(int(index.g_to_strand_arr[g_idx]))
                    if g_strand in (Strand.POS, Strand.NEG):
                        strand_models.intronic.observe(result.exon_strand, g_strand)

                # Insert size model
                if result.insert_size > 0:
                    insert_models.observe(
                        result.insert_size, result.count_cat
                    )
                    stats.n_insert_unambiguous += 1
                else:
                    stats.n_insert_ambiguous += 1

            # --- Buffer ALL resolved fragments with frag_id ---
            buffer.append(result, frag_id=group_id)

            if not is_unique_mapper:
                n_buffered_mm += 1

        # Track multimapper molecule stats
        if n_buffered_mm > 0:
            stats.n_multimapper_groups += 1
            stats.n_multimapper_alignments += n_buffered_mm

        group_id += 1
        if group_id % log_every == 0:
            logger.debug(
                f"  {group_id:,} read-name groups, "
                f"{stats.n_fragments:,} fragments, "
                f"{strand_models.n_observations:,} strand obs, "
                f"{buffer.total_fragments:,} buffered"
            )

    buffer.finalize()

    logger.info(
        f"[DONE] Scan: {stats.n_fragments:,} fragments → "
        f"{stats.n_unique_gene:,} unique-gene, "
        f"{stats.n_multi_gene:,} multi-gene, "
        f"{stats.n_intergenic:,} intergenic, "
        f"{buffer.total_fragments:,} buffered "
        f"({buffer.memory_bytes / 1024**2:.1f} MB in memory, "
        f"{buffer.n_spilled} chunks spilled), "
        f"{stats.n_multimapper_groups:,} multimapper molecules"
    )
    strand_models.log_summary()
    insert_models.log_summary()

    return stats, strand_models, insert_models, buffer


# ---------------------------------------------------------------------------
# EM scoring helper
# ---------------------------------------------------------------------------


def _score_candidate(
    exon_strand: int,
    t_strand: int,
    count_cat: int,
    insert_size: int,
    strand_model: StrandModel,
    insert_models: InsertSizeModels,
) -> tuple[float, int]:
    """Compute data log-likelihood and count type for one candidate.

    Returns
    -------
    log_lik : float
        ``log(P_strand) + log(P_insert)``.
    count_type : int
        Flat ``CountType`` column index (0-11).
    """
    p_strand = strand_model.strand_likelihood(
        exon_strand, Strand(t_strand)
    )

    cat_model = insert_models.category_models.get(
        count_cat, insert_models.global_model
    )
    log_p_insert = (
        cat_model.log_likelihood(insert_size)
        if insert_size > 0
        else 0.0
    )

    log_lik = np.log(max(p_strand, 1e-10)) + log_p_insert

    count_strand = ReadCounter.classify_strand(
        exon_strand, t_strand, strand_model
    )
    count_type = int(count_cat) * 3 + int(count_strand)

    return log_lik, count_type


# Default gDNA splice penalties per CountCategory.
# INTRON and UNSPLICED are fully compatible with gDNA (penalty 1.0).
# SPLICED_UNANNOT gets a small penalty (most are real SJs, but some
# are alignment artifacts from gDNA).
# SPLICED_ANNOT never competes with gDNA (handled upstream).
_GDNA_SPLICE_PENALTIES = {
    CountCategory.INTRON: 1.0,
    CountCategory.UNSPLICED: 1.0,
    CountCategory.SPLICED_UNANNOT: 0.01,
    CountCategory.SPLICED_ANNOT: 0.0,  # never used
}


def _score_gdna_candidate(
    count_cat: int,
    insert_size: int,
    insert_models: InsertSizeModels,
    gdna_splice_penalties: dict | None = None,
) -> float:
    """Compute data log-likelihood for the gDNA pseudo-component.

    gDNA is modelled as:
    - **Unstranded**: ``P_strand = 0.5`` (gDNA has no strand bias).
    - **Insert size**: uses the *intergenic* insert-size model
      (learned from intergenic fragments which are high-confidence gDNA).
    - **Splice penalty**: INTRON/UNSPLICED get 1.0 (fully compatible),
      SPLICED_UNANNOT gets a configurable penalty (default 0.01).

    Returns
    -------
    float
        ``log(P_strand) + log(P_insert) + log(splice_penalty)``.
    """
    penalties = gdna_splice_penalties or _GDNA_SPLICE_PENALTIES
    splice_pen = penalties.get(CountCategory(count_cat), 1.0)

    log_p_strand = np.log(0.5)

    intergenic_model = insert_models.intergenic
    log_p_insert = (
        intergenic_model.log_likelihood(insert_size)
        if insert_size > 0
        else 0.0
    )

    log_p_splice = np.log(max(splice_pen, 1e-10))

    return log_p_strand + log_p_insert + log_p_splice


# ---------------------------------------------------------------------------
# Buffer scan + EM data construction
# ---------------------------------------------------------------------------


def _scan_and_build_em_data(
    buffer: FragmentBuffer,
    index: HulkIndex,
    strand_model: StrandModel,
    insert_models: InsertSizeModels,
    counter: ReadCounter,
    stats: PipelineStats,
    log_every: int,
    gdna_splice_penalties: dict | None = None,
) -> EMData:
    """Single buffer pass: assign uniques + pre-compute EM data with per-gene gDNA.

    Iterates the buffer once.  Fragment routing depends on category:

    - **SPLICED_ANNOT + FRAG_UNIQUE**: deterministic assignment to the
      sole transcript (annotated SJs prove RNA origin).
    - **All other non-chimeric fragments**: enter the EM with their
      transcript candidates *plus* per-gene gDNA shadow candidate(s).
      Each fragment's shadow candidate uses the gene index of the
      highest-scoring transcript candidate.

    Multimapper alignments sharing a ``frag_id`` are grouped into a
    single EM unit (one molecule = one count).

    Parameters
    ----------
    buffer, index, strand_model, insert_models, counter, stats, log_every
        Pipeline components (counter and stats are updated in-place).
    gdna_splice_penalties : dict or None
        Override default gDNA splice penalties per CountCategory.

    Returns
    -------
    EMData
    """
    gdna_base = counter.gdna_base_index
    t_to_g = index.t_to_g_arr

    # Accumulation lists for CSR (converted to arrays at end)
    offsets: list[int] = [0]
    t_indices_list: list[int] = []
    log_liks_list: list[float] = []
    count_types_list: list[int] = []

    # Per-unit locus tracking (best transcript per EM unit)
    locus_t_list: list[int] = []
    locus_ct_list: list[int] = []

    # Multimapper group state -- persists across chunks because a
    # molecule's alignments may span a chunk boundary.
    mm_pending: list = []
    mm_fid: int = -1

    def _add_transcript_candidates(bf) -> tuple[float, int]:
        """Add all transcript candidates for a fragment.

        Returns (best_log_lik, best_count_type) among transcripts
        for locus tracking.  Returns (-inf, 0) if no candidates.
        """
        best_ll = -np.inf
        best_ct = 0
        for t_idx in bf.t_inds:
            t_idx_int = int(t_idx)
            t_strand = int(index.t_to_strand_arr[t_idx_int])
            log_lik, ct = _score_candidate(
                bf.exon_strand,
                t_strand,
                bf.count_cat,
                bf.insert_size,
                strand_model,
                insert_models,
            )
            t_indices_list.append(t_idx_int)
            log_liks_list.append(log_lik)
            count_types_list.append(ct)
            if log_lik > best_ll:
                best_ll = log_lik
                best_ct = ct
        return best_ll, best_ct

    def _add_gdna_candidate(bf, best_t_idx: int) -> None:
        """Add a per-gene gDNA shadow candidate for a fragment.

        The shadow index is ``gdna_base + gene_index`` where
        ``gene_index`` comes from the best transcript candidate.
        """
        if best_t_idx < 0:
            return  # No transcript candidate → no gene to shadow
        gdna_ll = _score_gdna_candidate(
            bf.count_cat, bf.insert_size, insert_models,
            gdna_splice_penalties,
        )
        g_idx = int(t_to_g[best_t_idx])
        shadow_idx = gdna_base + g_idx
        t_indices_list.append(shadow_idx)
        log_liks_list.append(gdna_ll)
        # gDNA candidates use count_type 0 (placeholder; not used for
        # gDNA since gDNA assignments go to gdna_locus_counts instead)
        count_types_list.append(0)

    def _flush_mm_group() -> None:
        """Flush accumulated multimapper alignments as one EM unit."""
        best_ll = -np.inf
        best_ct = 0
        best_t = -1
        any_spliced_annot = False
        for bf in mm_pending:
            if bf.count_cat == int(CountCategory.SPLICED_ANNOT):
                any_spliced_annot = True
            ll, ct = _add_transcript_candidates(bf)
            if ll > best_ll:
                best_ll = ll
                best_ct = ct
        # Find best transcript candidate for this unit (before adding shadow)
        unit_start_pre = offsets[-1]
        unit_end_pre = len(t_indices_list)
        _best_ll = -np.inf
        _best_t = -1
        _best_ct = 0
        for j in range(unit_start_pre, unit_end_pre):
            if log_liks_list[j] > _best_ll:
                _best_ll = log_liks_list[j]
                _best_t = t_indices_list[j]
                _best_ct = count_types_list[j]
        # Add per-gene gDNA shadow candidate unless all hits are SPLICED_ANNOT
        if not any_spliced_annot and mm_pending and _best_t >= 0:
            _add_gdna_candidate(mm_pending[0], _best_t)
        offsets.append(len(t_indices_list))
        locus_t_list.append(_best_t)
        locus_ct_list.append(_best_ct)

    n_processed = 0
    for chunk in buffer.iter_chunks():
        frag_classes = chunk.fragment_classes

        for i in range(chunk.size):
            fc = frag_classes[i]

            if fc == FRAG_CHIMERIC:
                continue  # Skip chimeric fragments during counting

            bf = chunk[i]
            is_spliced_annot = bf.count_cat == int(CountCategory.SPLICED_ANNOT)

            if fc == FRAG_UNIQUE and is_spliced_annot:
                # --- Deterministic: SPLICED_ANNOT + unique ---
                # Annotated SJs prove RNA origin; no gDNA competition.
                counter.assign_unique(bf, index, strand_model)
                stats.n_counted_truly_unique += 1

            elif fc == FRAG_MULTIMAPPER:
                fid = int(chunk.frag_id[i])
                if fid != mm_fid:
                    if mm_pending:
                        _flush_mm_group()
                        stats.n_counted_multimapper += 1
                    mm_fid = fid
                    mm_pending.clear()
                    mm_pending.append(bf)
                else:
                    mm_pending.append(bf)

            else:
                # --- EM unit: transcript candidates + per-gene gDNA shadow ---
                # This covers:
                # - FRAG_UNIQUE + non-SPLICED_ANNOT (competes with gDNA)
                # - FRAG_ISOFORM_AMBIG (any category)
                # - FRAG_GENE_AMBIG (any category)
                best_ll, best_ct = _add_transcript_candidates(bf)

                # Find best transcript candidate for shadow gene assignment
                unit_start = offsets[-1]
                unit_end = len(t_indices_list)
                _best_ll = -np.inf
                _best_t = -1
                _best_ct = 0
                for j in range(unit_start, unit_end):
                    if log_liks_list[j] > _best_ll:
                        _best_ll = log_liks_list[j]
                        _best_t = t_indices_list[j]
                        _best_ct = count_types_list[j]

                # Add per-gene gDNA shadow for non-SPLICED_ANNOT fragments
                if not is_spliced_annot and _best_t >= 0:
                    _add_gdna_candidate(bf, _best_t)

                if len(t_indices_list) > offsets[-1]:
                    offsets.append(len(t_indices_list))
                    locus_t_list.append(_best_t)
                    locus_ct_list.append(_best_ct)

                if fc == FRAG_UNIQUE:
                    stats.n_counted_truly_unique += 1
                elif fc == FRAG_ISOFORM_AMBIG:
                    stats.n_counted_isoform_ambig += 1
                else:
                    stats.n_counted_gene_ambig += 1

            n_processed += 1
            if n_processed % log_every == 0:
                logger.debug(
                    f"  Scan: {n_processed:,} / "
                    f"{buffer.total_fragments:,}"
                )

    # Flush last multimapper group
    if mm_pending:
        _flush_mm_group()
        stats.n_counted_multimapper += 1

    n_units = len(offsets) - 1
    n_candidates = len(t_indices_list)

    return EMData(
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(t_indices_list, dtype=np.int32),
        log_liks=np.array(log_liks_list, dtype=np.float64),
        count_types=np.array(count_types_list, dtype=np.uint8),
        locus_t_indices=np.array(locus_t_list, dtype=np.int32),
        locus_count_types=np.array(locus_ct_list, dtype=np.uint8),
        n_units=n_units,
        n_candidates=n_candidates,
        gdna_base_index=gdna_base,
    )


# ---------------------------------------------------------------------------
# Hierarchical per-gene shadow initialization (empirical Bayes)
# ---------------------------------------------------------------------------


def compute_shadow_init(
    unique_counts: np.ndarray,
    t_to_g_arr: np.ndarray,
    num_genes: int,
    strand_model: StrandModel,
    n_intergenic: int,
    *,
    kappa: float = 1.0,
) -> np.ndarray:
    """Compute per-gene gDNA shadow priors via empirical Bayes shrinkage.

    Uses the "Iceberg" model of genomic DNA contamination: gDNA
    distributes equally on sense and antisense strands.  Observable
    antisense reads in genic regions are the tip of the iceberg; an
    equal mass of sense gDNA is hidden among RNA reads.

    Global gDNA rate ("Iceberg" estimator)::

        θ_global = (n_intergenic + 2 × Σ antisense) / n_total

    - **Intergenic** reads are 100% gDNA (both strands visible).
    - **Antisense genic** reads are ~100% gDNA; the matching sense
      gDNA is hidden among RNA sense reads → 2× multiplier.

    Per-gene initialization (depth-scaled shrinkage)::

        shadow_init_g = 2 × antisense_g + κ × depth_g × θ_global

    where:
    - ``antisense_g`` — sum of antisense unique counts for gene *g*
    - ``depth_g`` — total unique counts for gene *g* (all categories
      and strands)
    - ``κ`` — shrinkage strength toward global estimate

    The depth-scaling ensures the global prior is proportional to
    each gene's sequencing depth: a deeply-sequenced gene gets a
    larger absolute gDNA floor than a lowly-sequenced gene.

    Parameters
    ----------
    unique_counts : np.ndarray
        (N_t, 12) unique counts array from ReadCounter.
    t_to_g_arr : np.ndarray
        Transcript-to-gene mapping.
    num_genes : int
        Number of genes.
    strand_model : StrandModel
        Trained exonic strand model.
    n_intergenic : int
        Number of intergenic fragments (deterministic gDNA).
    kappa : float
        Shrinkage strength (default 1.0).

    Returns
    -------
    np.ndarray
        float64[num_genes] — per-gene shadow initialization values.
    """
    from .categories import CountStrand

    n_total = float(unique_counts.sum()) + n_intergenic
    if n_total == 0:
        return np.zeros(num_genes, dtype=np.float64)

    # Per-transcript total unique counts → per-gene depth
    t_depth = unique_counts.sum(axis=1)  # (N_t,)
    gene_depth = np.zeros(num_genes, dtype=np.float64)
    np.add.at(gene_depth, t_to_g_arr, t_depth)

    # Per-gene antisense counts
    # Antisense columns: indices 2, 5, 8, 11 (stride 3, offset 2)
    antisense_cols = [int(CountStrand.ANTISENSE) + cat * 3 for cat in range(4)]
    t_antisense = unique_counts[:, antisense_cols].sum(axis=1)  # (N_t,)
    g_antisense = np.zeros(num_genes, dtype=np.float64)
    np.add.at(g_antisense, t_to_g_arr, t_antisense)

    total_antisense = g_antisense.sum()

    # Global gDNA rate — "Iceberg" estimator:
    # Intergenic = 100% gDNA (both strands visible, no doubling).
    # Antisense genic = gDNA tip; matching sense gDNA hidden → 2× multiplier.
    estimated_gdna = n_intergenic + 2.0 * total_antisense
    theta_global = estimated_gdna / n_total

    # Per-gene shadow initialization — depth-scaled shrinkage
    shadow_init = 2.0 * g_antisense + kappa * gene_depth * theta_global

    return shadow_init


# ---------------------------------------------------------------------------
# Count from buffer (unified EM)
# ---------------------------------------------------------------------------


def count_from_buffer(
    buffer: FragmentBuffer,
    index: HulkIndex,
    strand_model: StrandModel,
    insert_models: InsertSizeModels,
    stats: PipelineStats,
    *,
    seed: int | None = None,
    alpha: float = 1.0,
    em_iterations: int = 10,
    log_every: int = 1_000_000,
    gdna_splice_penalties: dict | None = None,
) -> ReadCounter:
    """Assign counts from buffered fragments via unified EM with per-gene gDNA.

    Every fragment contributes exactly **one discrete count** to a
    single component (transcript or per-gene gDNA shadow).

    1. **SPLICED_ANNOT + unique** fragments are assigned deterministically
       to their transcript (annotated SJs prove RNA origin).
    2. **All other** fragments enter the EM with transcript candidates
       *plus* a per-gene gDNA shadow, then are assigned stochastically
       using the converged posteriors.
    3. **Shadow initialization** uses empirical Bayes shrinkage
       (``compute_shadow_init``) from unique antisense counts and the
       global gDNA rate ("Iceberg" estimator: intergenic + 2× antisense).

    Parameters
    ----------
    buffer : FragmentBuffer
        Finalized buffer from :func:`scan_and_buffer`.
    index : HulkIndex
        The loaded reference index.
    strand_model : StrandModel
        Trained exonic strand model.
    insert_models : InsertSizeModels
        Trained insert size models.
    stats : PipelineStats
        Stats object (updated in-place with counting stats).
    seed : int or None
        Random seed for reproducibility.
    alpha : float
        Dirichlet pseudocount for prior smoothing (default 1.0).
    em_iterations : int
        Maximum number of EM iterations (default 10).
    log_every : int
        Log progress every *log_every* fragments.
    gdna_splice_penalties : dict or None
        Override default gDNA splice penalties per CountCategory.

    Returns
    -------
    ReadCounter
    """
    # First pass: assign SPLICED_ANNOT uniques + build EM data.
    # We create an initial counter to accumulate unique_counts, which
    # we then use to compute shadow_init via empirical Bayes.
    counter = ReadCounter(
        index.num_transcripts, index.num_genes,
        seed=seed, alpha=alpha,
    )

    logger.info(
        f"[START] Counting {buffer.total_fragments:,} buffered fragments "
        f"(per-gene gDNA shadows, {index.num_genes} genes)"
    )

    # Single buffer pass: assign uniques + build EM data
    em_data = _scan_and_build_em_data(
        buffer, index, strand_model, insert_models,
        counter, stats, log_every,
        gdna_splice_penalties=gdna_splice_penalties,
    )

    # Compute per-gene shadow initialization via empirical Bayes
    shadow_init = compute_shadow_init(
        counter.unique_counts,
        index.t_to_g_arr,
        index.num_genes,
        strand_model,
        stats.n_intergenic,
    )
    counter.shadow_init = shadow_init

    logger.info(
        f"[DONE] Scan: {stats.n_counted_truly_unique:,} unique, "
        f"{em_data.n_units:,} ambiguous units "
        f"({em_data.n_candidates:,} candidates), "
        f"shadow_init sum={shadow_init.sum():.1f}"
    )

    # EM iteration + stochastic assignment
    if em_data.n_units > 0:
        logger.info(
            f"[START] EM: up to {em_iterations} iterations over "
            f"{em_data.n_units:,} ambiguous units"
        )
        counter.run_em(em_data, em_iterations=em_iterations)

        logger.info("[START] Stochastic assignment of ambiguous fragments")
        counter.assign_ambiguous(em_data)
        logger.info(
            f"[DONE] Assigned {em_data.n_units:,} ambiguous fragments "
            f"(gDNA EM: {counter.gdna_em_count:.0f})"
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")
        counter.run_em(em_data, em_iterations=0)

    # --- Update stats with gDNA results ---
    stats.n_gdna_em = int(counter.gdna_em_count)

    logger.info(
        f"[gDNA] total={counter.gdna_total:.0f} "
        f"(intergenic={counter.gdna_unique_count:.0f}, "
        f"EM={counter.gdna_em_count:.0f}), "
        f"contamination rate={counter.gdna_contamination_rate:.2%}"
    )

    return counter


# ---------------------------------------------------------------------------
# Full pipeline orchestration
# ---------------------------------------------------------------------------

def run_pipeline(
    bam_path,
    index: HulkIndex,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    insert_size_max: int = 1000,
    log_every: int = 1_000_000,
    chunk_size: int = 1_000_000,
    max_memory_bytes: int = 2 * 1024**3,
    spill_dir: Path | str | None = None,
    sj_strand_tag: str | tuple[str, ...] = "auto",
    seed: int | None = None,
    alpha: float = 1.0,
    em_iterations: int = 10,
    gdna_splice_penalty_unannot: float = 0.01,
) -> PipelineResult:
    """Run the complete counting pipeline with unified EM and gDNA modeling.

    Opens the BAM **once**, scans all fragments, trains models,
    buffers resolved fragments, then counts via EM-based assignment.
    The EM includes a genomic DNA pseudo-component that competes
    with transcripts for non-SPLICED_ANNOT fragments.

    Parameters
    ----------
    bam_path : str or Path
        Path to the name-sorted / collated BAM file.
    index : HulkIndex
        The loaded reference index.
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include multimapping reads (default False).
    insert_size_max : int
        Maximum insert size for histogram models.
    log_every : int
        Log progress every *log_every* read-name groups.
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory for in-memory buffer chunks before disk spill
        (default 2 GiB).  Set to 0 to disable spilling.
    spill_dir : Path, str, or None
        Directory for disk-spilled buffer chunks
        (default: system temp dir).
    sj_strand_tag : str, tuple of str, or ``"auto"``
        BAM tag(s) for splice-junction strand.  ``"auto"``
        (default) scans the BAM for spliced reads and detects
        which tag is present.  Pass ``"XS"`` for STAR,
        ``"ts"`` for minimap2, or a tuple like
        ``("XS", "ts")`` to check multiple tags in order.
    seed : int or None
        Random seed for reproducibility (default None).
    alpha : float
        Dirichlet pseudocount for prior smoothing (default 1.0).
    em_iterations : int
        Maximum number of EM iterations (default 10).
    gdna_splice_penalty_unannot : float
        gDNA splice penalty for SPLICED_UNANNOT fragments (default 0.01).
        Higher values make gDNA more competitive for unannotated
        spliced fragments (treating them more like alignment artifacts).

    Returns
    -------
    PipelineResult
    """
    import pysam

    bam_path = str(bam_path)
    if spill_dir is not None:
        spill_dir = Path(spill_dir)

    # -- Resolve sj_strand_tag --
    if sj_strand_tag == "auto":
        detected = detect_sj_strand_tag(bam_path)
        if len(detected) == 1:
            sj_strand_tag = detected[0]
        elif len(detected) > 1:
            sj_strand_tag = detected
        else:
            # No tags found — all SJs will get Strand.NONE
            sj_strand_tag = ()

    # -- Single BAM pass: scan + train + buffer --
    bamfh = pysam.AlignmentFile(bam_path, "rb")
    try:
        stats, strand_models, insert_models, buffer = (
            scan_and_buffer(
                bamfh.fetch(until_eof=True),
                index,
                skip_duplicates=skip_duplicates,
                include_multimap=include_multimap,
                log_every=log_every,
                insert_size_max=insert_size_max,
                chunk_size=chunk_size,
                max_memory_bytes=max_memory_bytes,
                spill_dir=spill_dir,
                sj_strand_tag=sj_strand_tag,
            )
        )
    finally:
        bamfh.close()

    # -- Count assignment via unified EM --
    # Build gDNA splice penalty dict from the user-facing parameter
    gdna_splice_penalties = dict(_GDNA_SPLICE_PENALTIES)
    gdna_splice_penalties[CountCategory.SPLICED_UNANNOT] = (
        gdna_splice_penalty_unannot
    )

    try:
        counter = count_from_buffer(
            buffer,
            index,
            strand_models.exonic,
            insert_models,
            stats,
            seed=seed,
            alpha=alpha,
            em_iterations=em_iterations,
            log_every=log_every,
            gdna_splice_penalties=gdna_splice_penalties,
        )
    finally:
        buffer.cleanup()

    return PipelineResult(
        stats=stats,
        strand_models=strand_models,
        insert_models=insert_models,
        counter=counter,
    )
