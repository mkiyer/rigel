"""
rigel.annotate — Per-fragment annotated BAM output.

Produces a BAM file where every original record is stamped with EM-derived
assignment tags, enabling read-level introspection of the quantification pipeline.

Two-pass architecture
---------------------
Pass 1 (the normal pipeline) builds a lightweight **annotation table** —
parallel numpy arrays keyed by ``frag_id`` — recording the pipeline's
assignment decision for each fragment:

- ``best_tid`` — transcript index (-1 = intergenic / gDNA / unassigned)
- ``best_gid`` — gene index (-1 = intergenic / unassigned)
- ``tx_flags`` — assignment flags bitfield (is_resolved, is_gdna, is_nrna, is_synthetic)
- ``posterior``— posterior probability of the winning assignment
- ``frag_class`` — fragment class (unambig / ambig_same_strand / …)
- ``n_candidates`` — number of competing EM candidates
- ``splice_type`` — splice classification

Pass 2 re-reads the input BAM (same name-sorted order) and stamps each
record with standard BAM tags.  This pass is implemented entirely in
C++ (``_bam_impl.BamAnnotationWriter``) using the same htslib-based
infrastructure as pass 1.

BAM Tag Schema
--------------
.. list-table::
   :widths: 5 5 40
   :header-rows: 1

   * - Tag
     - Type
     - Description
   * - ZT
     - Z
     - Transcript ID (``"."`` for intergenic / gDNA-only)
   * - ZG
     - Z
     - Gene ID (``"."`` for intergenic)
   * - ZR
     - Z
     - Gene name / symbol (``"."`` for intergenic)
   * - ZI
     - i
     - Transcript index into rigel reference (``-1`` if unassigned)
   * - ZJ
     - i
     - Gene index into rigel reference (``-1`` if unassigned)
   * - ZF
     - i
     - Assignment flags bitfield:
       bit 0 = is_resolved (fragment scored and assigned by EM),
       bit 1 = is_gdna (EM gDNA component won),
       bit 2 = is_nrna (assigned transcript is single-exon),
       bit 3 = is_synthetic (assigned transcript is rigel-generated nRNA span)
   * - ZW
     - f
     - Posterior probability of assignment
   * - ZC
     - Z
     - Fragment class: ``unambig``, ``ambig_same_strand``,
       ``ambig_opp_strand``, ``multimapper``, ``chimeric``, ``intergenic``
   * - ZH
     - i
     - Primary hit flag: 1 = winning alignment, 0 = secondary
   * - ZN
     - i
     - Number of EM candidate components
   * - ZS
     - Z
     - Splice type (``spliced_annot``, ``spliced_unannot``,
       ``unspliced``, ``unknown``)
   * - ZL
     - i
     - Locus ID (``-1`` if no locus)
"""

import logging
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Assignment flag bits (written to BAM as the ZF:i tag)
# ---------------------------------------------------------------------------
AF_RESOLVED: int = 0x1   # bit 0: fragment was scored and assigned by EM
AF_GDNA: int = 0x2       # bit 1: EM gDNA component won
AF_NRNA: int = 0x4       # bit 2: assigned transcript is single-exon
AF_SYNTHETIC: int = 0x8  # bit 3: assigned transcript is rigel-generated nRNA span

# Pre-computed valid assignment flag values for common outcomes.
AF_UNRESOLVED: int = 0                                      # 0  — not modeled
AF_TRANSCRIPT: int = AF_RESOLVED                             # 1  — multi-exon transcript
AF_GDNA_RESOLVED: int = AF_RESOLVED | AF_GDNA               # 3  — gDNA component
AF_NRNA_RESOLVED: int = AF_RESOLVED | AF_NRNA               # 5  — single-exon annotated
AF_SYNTH_RESOLVED: int = AF_RESOLVED | AF_NRNA | AF_SYNTHETIC  # 13 — synthetic nRNA span

# Fragment-class labels for the ZC tag.
_FRAG_CLASS_LABELS = {
    0: "unambig",
    1: "ambig_same_strand",
    2: "ambig_opp_strand",
    3: "multimapper",
    4: "chimeric",
    -1: "intergenic",
}


@dataclass(slots=True)
class AnnotationTable:
    """Lightweight per-fragment annotation table.

    Backed by parallel numpy arrays, one entry per annotated fragment.
    Fragments are looked up by ``frag_id`` via :pyattr:`frag_id_to_row`.

    Attributes
    ----------
    capacity : int
        Pre-allocated array size (= total frag_ids).
    frag_ids : np.ndarray
        int64 — frag_id for each row (populated densely from 0).
    best_tid : np.ndarray
        int32 — assigned transcript index (-1 = none).
    best_gid : np.ndarray
        int32 — assigned gene index (-1 = none).
    tx_flags : np.ndarray
        uint8 — assignment flags bitfield (see AF_* constants).
    posterior : np.ndarray
        float32 — posterior probability of assignment.
    frag_class : np.ndarray
        int8 — fragment class code (matches buffer FRAG_* constants, -1=intergenic).
    n_candidates : np.ndarray
        int16 — number of EM candidate components.
    splice_type : np.ndarray
        uint8 — SpliceType enum value.
    _size : int
        Number of rows currently written.
    frag_id_to_row : dict[int, int]
        Maps frag_id → row index for O(1) lookup.
    """

    capacity: int
    frag_ids: np.ndarray = field(repr=False)
    best_tid: np.ndarray = field(repr=False)
    best_gid: np.ndarray = field(repr=False)
    tx_flags: np.ndarray = field(repr=False)
    posterior: np.ndarray = field(repr=False)
    frag_class: np.ndarray = field(repr=False)
    n_candidates: np.ndarray = field(repr=False)
    splice_type: np.ndarray = field(repr=False)
    locus_id: np.ndarray = field(repr=False)
    _size: int = 0
    frag_id_to_row: dict = field(default_factory=dict, repr=False)

    @classmethod
    def create(cls, capacity: int) -> "AnnotationTable":
        """Allocate an empty table with the given capacity."""
        return cls(
            capacity=capacity,
            frag_ids=np.full(capacity, -1, dtype=np.int64),
            best_tid=np.full(capacity, -1, dtype=np.int32),
            best_gid=np.full(capacity, -1, dtype=np.int32),
            tx_flags=np.zeros(capacity, dtype=np.uint8),
            posterior=np.zeros(capacity, dtype=np.float32),
            frag_class=np.full(capacity, -1, dtype=np.int8),
            n_candidates=np.zeros(capacity, dtype=np.int16),
            splice_type=np.zeros(capacity, dtype=np.uint8),
            locus_id=np.full(capacity, -1, dtype=np.int32),
        )

    def add(
        self,
        frag_id: int,
        best_tid: int = -1,
        best_gid: int = -1,
        tx_flags: int = AF_UNRESOLVED,
        posterior: float = 0.0,
        frag_class: int = -1,
        n_candidates: int = 0,
        splice_type: int = 0,
    ) -> None:
        """Append one annotation row."""
        idx = self._size
        if idx >= self.capacity:
            self._grow()
        self.frag_ids[idx] = frag_id
        self.best_tid[idx] = best_tid
        self.best_gid[idx] = best_gid
        self.tx_flags[idx] = tx_flags
        self.posterior[idx] = posterior
        self.frag_class[idx] = frag_class
        self.n_candidates[idx] = n_candidates
        self.splice_type[idx] = splice_type
        self.frag_id_to_row[frag_id] = idx
        self._size += 1

    def add_batch(
        self,
        frag_ids: np.ndarray,
        best_tids: np.ndarray,
        best_gids: np.ndarray,
        tx_flags: np.ndarray,
        posteriors: np.ndarray,
        frag_classes: np.ndarray,
        n_candidates: np.ndarray,
        splice_types: np.ndarray,
    ) -> None:
        """Append multiple annotation rows in batch."""
        n = len(frag_ids)
        if n == 0:
            return

        needed = self._size + n
        if needed > self.capacity:
            self._grow_to(max(self.capacity * 2, needed))

        start = self._size
        end = start + n

        self.frag_ids[start:end] = frag_ids
        self.best_tid[start:end] = best_tids
        self.best_gid[start:end] = best_gids
        self.tx_flags[start:end] = tx_flags
        self.posterior[start:end] = posteriors
        self.frag_class[start:end] = frag_classes
        self.n_candidates[start:end] = n_candidates
        self.splice_type[start:end] = splice_types

        self.frag_id_to_row.update(
            zip(frag_ids.tolist(), range(start, end))
        )

        self._size = end

    def _grow_to(self, new_cap: int) -> None:
        """Grow capacity to at least new_cap."""
        if new_cap <= self.capacity:
            return
        for attr in (
            "frag_ids", "best_tid", "best_gid", "tx_flags",
            "posterior", "frag_class", "n_candidates", "splice_type",
            "locus_id",
        ):
            old = getattr(self, attr)
            new = np.empty(new_cap, dtype=old.dtype)
            new[:self.capacity] = old
            if attr in ("frag_ids", "best_tid", "best_gid", "frag_class", "locus_id"):
                new[self.capacity:] = -1
            else:
                new[self.capacity:] = 0
            setattr(self, attr, new)
        self.capacity = new_cap

    def _grow(self) -> None:
        """Double capacity when full."""
        self._grow_to(max(self.capacity * 2, 1024))

    def get(self, frag_id: int):
        """Return annotation dict for a frag_id, or None if absent."""
        row = self.frag_id_to_row.get(frag_id)
        if row is None:
            return None
        return {
            "best_tid": int(self.best_tid[row]),
            "best_gid": int(self.best_gid[row]),
            "tx_flags": int(self.tx_flags[row]),
            "posterior": float(self.posterior[row]),
            "frag_class": int(self.frag_class[row]),
            "n_candidates": int(self.n_candidates[row]),
            "splice_type": int(self.splice_type[row]),
        }

    @property
    def size(self) -> int:
        return self._size


def _splice_type_label(code: int) -> str:
    """Convert SpliceType int to lowercase label for the ZS tag."""
    from .splice import SpliceType
    try:
        return SpliceType(code).name.lower()
    except (ValueError, KeyError):
        return "unknown"


def write_annotated_bam(
    bam_path: str,
    output_path: str,
    annotations: AnnotationTable,
    index,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    sj_strand_tag: str | tuple[str, ...] = "auto",
    locus_id_per_transcript: np.ndarray | None = None,
) -> dict:
    """Second BAM pass: stamp every record with assignment tags.

    Re-reads the input BAM in the same name-sorted order as pass 1.
    For each read-name group, looks up the annotation by ``frag_id``
    (maintained in lockstep with the original scan) and writes tags.

    Implemented entirely in C++ via ``_bam_impl.BamAnnotationWriter``
    using the same htslib infrastructure as pass 1 — no pysam needed.

    Parameters
    ----------
    bam_path : str
        Path to input name-sorted BAM.
    output_path : str
        Path to output annotated BAM.
    annotations : AnnotationTable
        The annotation table built during quantification.
    index : TranscriptIndex
        Reference index (for transcript / gene ID lookup).
    skip_duplicates : bool
        Must match the value used in Pass 1.
    include_multimap : bool
        Must match the value used in Pass 1.
    sj_strand_tag : str or tuple
        Must match the value used in Pass 1.
    locus_id_per_transcript : np.ndarray or None
        int32 array mapping transcript index → locus_id. If provided,
        used to populate the ZL (locus_id) BAM tag.

    Returns
    -------
    dict
        Summary statistics of the annotation pass.
    """
    from .native import BamAnnotationWriter

    # Convert sj_strand_tag to spec string
    if isinstance(sj_strand_tag, (list, tuple)):
        sj_spec = ",".join(sj_strand_tag) if sj_strand_tag else "none"
    else:
        sj_spec = sj_strand_tag if sj_strand_tag else "none"

    # Build writer using the same FragmentResolver as pass 1
    writer = BamAnnotationWriter(
        index.resolver,
        sj_spec,
        skip_duplicates,
        include_multimap,
    )

    # Slice annotation arrays to populated size
    n = annotations.size

    # Populate locus_id from best_tid → locus_id_per_transcript
    if locus_id_per_transcript is not None:
        best_tid = annotations.best_tid[:n]
        valid = best_tid >= 0
        annotations.locus_id[:n] = np.where(
            valid,
            locus_id_per_transcript[np.clip(best_tid, 0, None)],
            -1,
        )

    t_ids = list(index.t_df["t_id"].values)
    g_ids = list(index.g_df["g_id"].values)
    g_names = list(index.g_df["g_name"].values)

    summary = writer.write(
        bam_path,
        output_path,
        annotations.frag_ids[:n],
        annotations.best_tid[:n],
        annotations.best_gid[:n],
        annotations.tx_flags[:n],
        annotations.posterior[:n],
        annotations.frag_class[:n],
        annotations.n_candidates[:n],
        annotations.splice_type[:n],
        annotations.locus_id[:n],
        n,
        t_ids,
        g_ids,
        g_names,
    )

    logger.info(
        f"[ANNOTATE] Wrote {summary['n_records_written']:,} records "
        f"to {output_path} "
        f"({summary['n_annotated']:,} annotated, "
        f"{summary['n_intergenic']:,} intergenic, "
        f"{summary['n_filtered_passthrough']:,} filtered pass-through)"
    )

    return summary
