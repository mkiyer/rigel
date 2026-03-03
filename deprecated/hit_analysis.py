"""
hulkrna.hit_analysis — Pass 2 hit resolution and model learning.

After Pass 1 (``query_bam``) produces an intermediate Feather file of
fragment-to-reference hits, this module streams through those hits to:

1. **Resolve fragment ambiguity** — Merge per-exon-block and per-SJ
   transcript/gene index sets via progressive relaxation, tracking
   which merge criteria succeeded (intersection, intersection-nonempty,
   or union).

2. **Learn the strand protocol** — Accumulate observations from
   fragments with annotated splice junction matches to infer
   FR/RF/unstranded protocol and strand specificity.

3. **Learn the fragment length distribution** — Train an
   ``FragmentLengthModel`` histogram from fragments where all
   candidate transcripts agree on the same fragment length.

Performance opportunities
-------------------------
- **Parallelism**: Fragment resolution is embarrassingly parallel.
  Partition frag_id ranges across ``multiprocessing.Pool`` workers
  or use ``concurrent.futures.ProcessPoolExecutor``.
- **Vectorized resolution**: Replace per-fragment Python loops with
  NumPy/pandas vectorized operations for hit partitioning and
  set merging.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pandas as pd
import pyarrow as pa
import pyarrow.ipc

from .core import (
    EMPTY_MERGE,
    IntervalType,
    MergeOutcome,
    MergeResult,
    Strand,
    merge_sets_with_criteria,
)
from .frag_length_model import FragmentLengthModel
from .strand_model import StrandModel

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Resolved fragment
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class ResolvedFragment:
    """Resolved hit information for a single fragment after merging.

    Produced by :func:`resolve_fragment_hits` from the raw hit rows
    in the intermediate Feather file.

    Attributes
    ----------
    frag_id : int
        Fragment identifier (groups rows in the hit file).
    exon_merge : MergeResult
        Transcript/gene merge from EXON-type interval hits only.
    sj_merge : MergeResult
        Transcript/gene merge from annotated SJ hits.
        ``EMPTY_MERGE`` if the fragment has no SJ-type rows.
    exon_strand : Strand
        Bitwise OR of all exon block alignment strands (after r2 flip).
    sj_strand : Strand
        Bitwise OR of all annotated SJ reference strands.
    has_exon_overlap : bool
        True if any hit has ``interval_type == EXON``.
    has_intron_overlap : bool
        True if any hit has ``interval_type == INTRON``.
    has_intergenic_overlap : bool
        True if any hit has ``interval_type == INTERGENIC``.
    """
    frag_id: int
    exon_merge: MergeResult
    sj_merge: MergeResult
    exon_strand: Strand
    sj_strand: Strand
    has_exon_overlap: bool
    has_intron_overlap: bool
    has_intergenic_overlap: bool
    has_unannotated_sj: bool

    @property
    def has_annotated_sj(self) -> bool:
        """True if the fragment matched at least one annotated SJ."""
        return not self.sj_merge.is_empty

    @property
    def is_spliced(self) -> bool:
        """True if the fragment has any SJ (annotated or unannotated)."""
        return self.has_annotated_sj or self.has_unannotated_sj

    @property
    def is_unique_gene(self) -> bool:
        """True if the fragment maps to exactly one gene.

        Uses the exon merge result if available, otherwise falls back
        to the SJ merge result.
        """
        if not self.exon_merge.is_empty:
            return self.exon_merge.is_unique_gene
        return self.sj_merge.is_unique_gene

    @property
    def is_high_confidence(self) -> bool:
        """True if resolved via full intersection to a single gene."""
        return (
            self.exon_merge.is_unique_gene
            and self.exon_merge.criteria == MergeOutcome.INTERSECTION
        )


# ---------------------------------------------------------------------------
# Streaming parser
# ---------------------------------------------------------------------------

def iter_fragment_hits(
    hit_path: Path | str,
) -> Iterator[tuple[int, pd.DataFrame]]:
    """Stream hits from the intermediate Feather file, grouped by frag_id.

    Reads the Arrow IPC file in RecordBatch-sized chunks, carrying over
    any partial fragment group that spans a batch boundary.  Peak memory
    is bounded by the batch size written by Pass 1 plus one fragment's
    worth of carry-over rows.

    Yields ``(frag_id, hits_df)`` pairs in ascending frag_id order.

    Parameters
    ----------
    hit_path : Path or str
        Path to the Feather file produced by ``query_bam()``.

    Yields
    ------
    tuple[int, pd.DataFrame]
        ``(frag_id, hits)`` where *hits* is the subset of rows for
        that fragment.
    """
    with pa.OSFile(str(hit_path), "rb") as f:
        reader = pa.ipc.open_file(f)

        # Carry-over buffer for a partial fragment at a batch boundary
        carry: pd.DataFrame | None = None

        for batch_idx in range(reader.num_record_batches):
            batch_df = reader.get_batch(batch_idx).to_pandas()

            # Prepend carry-over from previous batch
            if carry is not None:
                batch_df = pd.concat([carry, batch_df], ignore_index=True)
                carry = None

            if batch_df.empty:
                continue

            # The last frag_id in this batch might be split across batches
            last_frag_id = int(batch_df["frag_id"].iloc[-1])

            # Split: complete fragments vs trailing partial
            complete_mask = batch_df["frag_id"] != last_frag_id

            # If the entire batch is one frag_id and it's not the last
            # batch, carry it all forward
            if not complete_mask.any() and batch_idx < reader.num_record_batches - 1:
                carry = batch_df
                continue

            # Yield all complete fragments
            complete = batch_df[complete_mask]
            if not complete.empty:
                for frag_id, group in complete.groupby("frag_id", sort=True):
                    yield int(frag_id), group

            # The trailing fragment goes to carry (or is yielded if last batch)
            trailing = batch_df[~complete_mask]
            if batch_idx < reader.num_record_batches - 1:
                carry = trailing
            else:
                # Last batch — yield the trailing fragment
                if not trailing.empty:
                    yield int(last_frag_id), trailing

        # If there's still a carry-over after all batches (shouldn't happen
        # normally, but defensive)
        if carry is not None and not carry.empty:
            last_id = int(carry["frag_id"].iloc[0])
            yield last_id, carry


# ---------------------------------------------------------------------------
# Fragment resolution
# ---------------------------------------------------------------------------

def resolve_fragment_hits(
    frag_id: int,
    hits: pd.DataFrame,
) -> ResolvedFragment:
    """Resolve a fragment's raw hits into merged transcript/gene sets.

    Partitions the hit rows into exon-block hits and splice-junction
    hits, builds per-block and per-SJ transcript/gene index sets, and
    merges each via progressive relaxation with criteria tracking.

    Exon merging uses only ``interval_type == EXON`` rows. INTRON and
    INTERGENIC hits are noted via boolean flags but do not participate
    in set merging.

    Parameters
    ----------
    frag_id : int
        Fragment identifier.
    hits : pd.DataFrame
        All hit rows for this fragment (columns per ``HIT_SCHEMA``).

    Returns
    -------
    ResolvedFragment
    """
    ivtypes = hits["interval_type"].values

    # Partition: annotated SJ, unannotated SJ, and everything else
    sj_mask = ivtypes == IntervalType.SJ
    sj_unannot_mask = ivtypes == IntervalType.SJ_UNANNOT
    block_mask = ~(sj_mask | sj_unannot_mask)
    block_hits = hits[block_mask]
    sj_hits = hits[sj_mask]

    # --- Exon block merging (EXON interval_type only) ---
    exon_only = block_hits[block_hits["interval_type"] == IntervalType.EXON]

    if not exon_only.empty:
        exon_t_sets: list[frozenset[int]] = []
        exon_g_sets: list[frozenset[int]] = []
        # Group by query exon block coordinates to reconstruct per-block sets
        for _key, block in exon_only.groupby(
            ["ref", "start", "end", "strand"]
        ):
            t_set = frozenset(
                int(v) for v in block["t_index"].values if v >= 0
            )
            g_set = frozenset(
                int(v) for v in block["g_index"].values if v >= 0
            )
            exon_t_sets.append(t_set)
            exon_g_sets.append(g_set)
        exon_merge = merge_sets_with_criteria(exon_t_sets, exon_g_sets)
    else:
        exon_merge = EMPTY_MERGE

    # --- Splice junction merging ---
    if not sj_hits.empty:
        sj_t_sets: list[frozenset[int]] = []
        sj_g_sets: list[frozenset[int]] = []
        for _key, sj_group in sj_hits.groupby(
            ["ref", "start", "end", "strand"]
        ):
            t_set = frozenset(
                int(v) for v in sj_group["t_index"].values if v >= 0
            )
            g_set = frozenset(
                int(v) for v in sj_group["g_index"].values if v >= 0
            )
            sj_t_sets.append(t_set)
            sj_g_sets.append(g_set)
        sj_merge = merge_sets_with_criteria(sj_t_sets, sj_g_sets)
    else:
        sj_merge = EMPTY_MERGE

    # --- Strand computation ---
    exon_strand = Strand.NONE
    for s in block_hits["strand"].values:
        exon_strand |= Strand(int(s))

    sj_strand = Strand.NONE
    for s in sj_hits["strand"].values:
        sj_strand |= Strand(int(s))

    # --- Interval type coverage ---
    itype_set = set(int(v) for v in ivtypes)

    return ResolvedFragment(
        frag_id=frag_id,
        exon_merge=exon_merge,
        sj_merge=sj_merge,
        exon_strand=exon_strand,
        sj_strand=sj_strand,
        has_exon_overlap=IntervalType.EXON in itype_set,
        has_intron_overlap=IntervalType.INTRON in itype_set,
        has_intergenic_overlap=IntervalType.INTERGENIC in itype_set,
        has_unannotated_sj=IntervalType.SJ_UNANNOT in itype_set,
    )


# ---------------------------------------------------------------------------
# Analysis statistics
# ---------------------------------------------------------------------------

@dataclass
class HitAnalysisStats:
    """Aggregate statistics from Pass 2 hit analysis."""

    n_fragments: int = 0

    # Interval type coverage
    n_with_exon_overlap: int = 0
    n_with_intron_overlap: int = 0
    n_with_intergenic_overlap: int = 0
    n_with_annotated_sj: int = 0
    n_with_unannotated_sj: int = 0

    # Gene ambiguity (from exon merge)
    n_unique_gene: int = 0
    n_multi_gene: int = 0
    n_no_gene: int = 0  # no EXON overlap → no gene match

    # Exon merge criteria distribution
    n_exon_intersection: int = 0
    n_exon_intersection_nonempty: int = 0
    n_exon_union: int = 0

    # SJ merge criteria distribution
    n_sj_intersection: int = 0
    n_sj_intersection_nonempty: int = 0
    n_sj_union: int = 0

    # Strand learning qualification
    n_strand_skipped_no_sj: int = 0
    n_strand_skipped_multi_gene: int = 0
    n_strand_skipped_ambiguous: int = 0

    # Fragment length learning
    n_insert_computed: int = 0
    n_insert_unambiguous: int = 0
    n_insert_ambiguous: int = 0
    n_insert_skipped: int = 0

    def to_dict(self) -> dict:
        """Return all fields as a plain dict."""
        return {k: v for k, v in self.__dict__.items()}


# ---------------------------------------------------------------------------
# Top-level analysis function
# ---------------------------------------------------------------------------

def analyze_hits(
    hit_path: Path | str,
    *,
    max_frag_length: int = 1000,
    log_every: int = 1_000_000,
) -> tuple[HitAnalysisStats, StrandModel, FragmentLengthModel]:
    """Run Pass 2 hit analysis: resolve fragments and learn models.

    Streams through the intermediate Feather file produced by Pass 1,
    resolves each fragment's hits via progressive set merging, and
    feeds qualified observations into a :class:`StrandModel` and an
    :class:`FragmentLengthModel`.

    Fragment length training
    --------------------
    For each fragment, the unique set of ``frag_length`` values (from
    rows with ``t_index >= 0`` and ``frag_length >= 0``) is examined:

    - If all transcripts agree on a single fragment length, the
      fragment is *unambiguous* and that size is observed.
    - If transcripts disagree, the fragment is *ambiguous* and
      skipped for training (deferred to Bayesian counting).

    Parameters
    ----------
    hit_path : Path or str
        Path to the intermediate Feather file from ``query_bam()``.
    max_frag_length : int
        Maximum fragment length for the histogram model.
    log_every : int
        Log progress every *log_every* fragments at DEBUG level.

    Returns
    -------
    stats : HitAnalysisStats
        Fragment-level statistics including ambiguity and merge criteria.
    strand_model : StrandModel
        Learned strand model (Beta posterior from qualified observations).
    frag_length_model : FragmentLengthModel
        Learned fragment length distribution.
    """
    hit_path = Path(hit_path)
    stats = HitAnalysisStats()
    strand_model = StrandModel()
    frag_length_model = FragmentLengthModel(max_size=max_frag_length)

    logger.info(f"[START] Analyzing hits from {hit_path}")

    for frag_id, hits in iter_fragment_hits(hit_path):
        resolved = resolve_fragment_hits(frag_id, hits)
        stats.n_fragments += 1

        # --- Interval type coverage ---
        if resolved.has_exon_overlap:
            stats.n_with_exon_overlap += 1
        if resolved.has_intron_overlap:
            stats.n_with_intron_overlap += 1
        if resolved.has_intergenic_overlap:
            stats.n_with_intergenic_overlap += 1
        if resolved.has_annotated_sj:
            stats.n_with_annotated_sj += 1
        if resolved.has_unannotated_sj:
            stats.n_with_unannotated_sj += 1

        # --- Gene ambiguity (from exon merge) ---
        if resolved.exon_merge.is_empty:
            stats.n_no_gene += 1
        elif resolved.exon_merge.is_unique_gene:
            stats.n_unique_gene += 1
        else:
            stats.n_multi_gene += 1

        # --- Exon merge criteria ---
        if not resolved.exon_merge.is_empty:
            if resolved.exon_merge.criteria == MergeOutcome.INTERSECTION:
                stats.n_exon_intersection += 1
            elif resolved.exon_merge.criteria == MergeOutcome.INTERSECTION_NONEMPTY:
                stats.n_exon_intersection_nonempty += 1
            elif resolved.exon_merge.criteria == MergeOutcome.UNION:
                stats.n_exon_union += 1

        # --- SJ merge criteria ---
        if not resolved.sj_merge.is_empty:
            if resolved.sj_merge.criteria == MergeOutcome.INTERSECTION:
                stats.n_sj_intersection += 1
            elif resolved.sj_merge.criteria == MergeOutcome.INTERSECTION_NONEMPTY:
                stats.n_sj_intersection_nonempty += 1
            elif resolved.sj_merge.criteria == MergeOutcome.UNION:
                stats.n_sj_union += 1

        # --- Strand learning (qualification + observation) ---
        if not resolved.has_annotated_sj:
            stats.n_strand_skipped_no_sj += 1
        elif not resolved.sj_merge.is_unique_gene:
            stats.n_strand_skipped_multi_gene += 1
        elif (
            resolved.exon_strand not in (Strand.POS, Strand.NEG)
            or resolved.sj_strand not in (Strand.POS, Strand.NEG)
        ):
            stats.n_strand_skipped_ambiguous += 1
        else:
            strand_model.observe(resolved.exon_strand, resolved.sj_strand)

        # --- Fragment length learning ---
        if "frag_length" in hits.columns:
            valid_mask = (hits["t_index"] >= 0) & (hits["frag_length"] >= 0)
            valid_sizes = hits.loc[valid_mask, "frag_length"]
            if valid_sizes.empty:
                stats.n_insert_skipped += 1
            else:
                unique_sizes = set(int(v) for v in valid_sizes.unique())
                stats.n_insert_computed += 1
                if len(unique_sizes) == 1:
                    stats.n_insert_unambiguous += 1
                    frag_length_model.observe(unique_sizes.pop())
                else:
                    stats.n_insert_ambiguous += 1

        # --- Progress ---
        if stats.n_fragments % log_every == 0:
            logger.debug(
                f"  analyzed {stats.n_fragments:,} fragments, "
                f"{strand_model.n_observations:,} strand observations, "
                f"{frag_length_model.n_observations:,} fragment length observations"
            )

    logger.info(
        f"[DONE] Analyzed {stats.n_fragments:,} fragments: "
        f"{stats.n_unique_gene:,} unique-gene, "
        f"{stats.n_multi_gene:,} multi-gene, "
        f"{stats.n_no_gene:,} no-gene"
    )
    strand_model.log_summary()
    frag_length_model.log_summary()

    return stats, strand_model, frag_length_model
