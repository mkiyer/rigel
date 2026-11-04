"""rigel.annotate — the annotated BAM: every input record stamped with its fragment's assignment.

The pipeline fills an :class:`AnnotationTable`, parallel arrays keyed by ``frag_id`` holding each fragment's
assignment. :func:`write_annotated_bam` then re-reads the input BAM in the same name-sorted order and stamps every
record, in C++ (``BamAnnotationWriter``), with these tags:

====  ====  ==============================================================================================
Tag   Type  Meaning
====  ====  ==============================================================================================
ZT    Z     transcript ID (``.`` for intergenic or gDNA-only)
ZG    Z     gene ID (``.`` for intergenic)
ZR    Z     gene name (``.`` for intergenic)
ZI    i     transcript index into the rigel reference (``-1`` if unassigned)
ZJ    i     gene index into the rigel reference (``-1`` if unassigned)
ZF    i     the outcome bitfield (the ``AF_*`` constants below)
ZW    f     posterior probability of the assignment
ZC    Z     input-ambiguity class: ``unambig``, ``ambig_same_strand``, ``ambig_opp_strand``, ``multimapper``,
            or ``.`` for a fragment that never entered the EM; ``unused_alignment`` on a placed read pair's
            record in none of its usable alignments (a rejected one's, one that pairs with nothing, a
            supplementary record); for a read pair withheld with no usable alignment (ZF 0), its reason:
            ``off_index``, ``chimeric``, ``unannotated_splice`` or ``unpaired``
ZH    i     1 on the winning placement's records, 0 otherwise
ZN    i     number of EM candidate components
ZS    Z     splice type of the record's own placement (``unspliced``, ``spliced_annot``,
            ``spliced_implicit``, ``splice_artifact``), or ``unknown`` on a record that is no usable alignment
ZL    i     locus ID (``-1`` if none)
ZB    i     junctions in this record's CIGAR that matched the splice-artifact blacklist and were scored as
            unspliced; per record, not per fragment, and 0 when no blacklist is loaded
====  ====  ==============================================================================================
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import CONSTANTS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ZF, the per-fragment outcome bitfield: one AND-mask per predicate answers "what happened to this fragment?"
#
#   bit 0 (0x01)  is_resolved             a candidate won, or the intergenic fast path took it
#   bit 1 (0x02)  is_annotated            the winner is a reference transcript, single-exon included
#   bit 2 (0x04)  is_gdna                 the winner is gDNA (unspliced intergenic + EM-gDNA)
#   bit 3 (0x08)  reserved                never set, so an older BAM's bit is never read with a new meaning
#   bit 4 (0x10)  is_nascent              the winner is a synthetic span, which models a multi-exon
#                                         transcript's unspliced RNA
#   bit 5 (0x20)  is_intergenic           no transcript holds it (the deterministic fast path), so it is
#                                         unspliced, and gDNA
#   bit 6 (0x40)  reserved                never set, as bit 3
#   bit 7 (0x80)  is_multimapper_dropped  unresolved because --no-include-multimap dropped it
#
# Bits 1, 2 and 4 record which kind of candidate won. Invariants (gate: tests/test_annotate.py):
#   1. is_resolved and is_multimapper_dropped are never both set; ZF = 0, unresolved, sets neither
#   2. exactly one winner-kind bit is set when is_resolved is, and none otherwise
#   3. is_intergenic implies is_resolved, and never comes with is_annotated or is_nascent
#   4. bits 3 and 6 are never set
# ---------------------------------------------------------------------------

# Primitive bits
AF_RESOLVED: int = 0x01
AF_ANNOTATED_BIT: int = 0x02
AF_GDNA_BIT: int = 0x04
AF_NASCENT_BIT: int = 0x10
AF_INTERGENIC_BIT: int = 0x20
AF_MULTIMAPPER_DROP_BIT: int = 0x80

# Canonical composed values (the only legitimate ZF outputs).
AF_UNRESOLVED: int = 0x00  # 0x00
AF_ANNOTATED: int = AF_RESOLVED | AF_ANNOTATED_BIT  # 0x03
AF_NASCENT: int = AF_RESOLVED | AF_NASCENT_BIT  # 0x11
AF_GDNA_EM: int = AF_RESOLVED | AF_GDNA_BIT  # 0x05
AF_GDNA_INTERGENIC: int = AF_GDNA_EM | AF_INTERGENIC_BIT  # 0x25
AF_MULTIMAPPER_DROP: int = AF_MULTIMAPPER_DROP_BIT  # 0x80


def winner_flag(is_synthetic: bool) -> int:
    """The ZF value of a fragment a transcript won: ``AF_NASCENT`` for a synthetic span, ``AF_ANNOTATED`` for every
    reference transcript. A gDNA winner is stamped ``AF_GDNA_EM`` at the call site instead."""
    return AF_NASCENT if is_synthetic else AF_ANNOTATED


def winner_flags(is_synthetic: np.ndarray) -> np.ndarray:
    """:func:`winner_flag` over the winning transcripts' ``is_synthetic`` flags: one ``uint8`` ZF value per winner."""
    is_synth_b = np.asarray(is_synthetic, dtype=bool)
    return np.where(is_synth_b, AF_NASCENT, AF_ANNOTATED).astype(np.uint8)


@dataclass(slots=True)
class AnnotationTable:
    """The per-fragment annotation table: parallel numpy arrays, one entry per annotated fragment.

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
    locus_id : np.ndarray
        int32 — locus id of the EM subproblem (-1 = none, or not yet filled).
    _size : int
        Number of rows currently written.
    """

    capacity: int
    frag_ids: np.ndarray = field(repr=False)
    best_tid: np.ndarray = field(repr=False)
    best_gid: np.ndarray = field(repr=False)
    tx_flags: np.ndarray = field(repr=False)
    posterior: np.ndarray = field(repr=False)
    frag_class: np.ndarray = field(repr=False)
    n_candidates: np.ndarray = field(repr=False)
    locus_id: np.ndarray = field(repr=False)
    _size: int = 0

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
        locus_ids: np.ndarray | None = None,
    ) -> None:
        """Append a batch of annotation rows.

        ``locus_ids`` is each fragment's locus from the EM partition. When omitted, ``locus_id`` stays -1 and
        :func:`write_annotated_bam` fills it from the winning transcript.
        """
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
        if locus_ids is not None:
            self.locus_id[start:end] = locus_ids

        self._size = end

    def _grow_to(self, new_cap: int) -> None:
        """Grow capacity to at least new_cap."""
        if new_cap <= self.capacity:
            return
        for attr in (
            "frag_ids",
            "best_tid",
            "best_gid",
            "tx_flags",
            "posterior",
            "frag_class",
            "n_candidates",
            "locus_id",
        ):
            old = getattr(self, attr)
            new = np.empty(new_cap, dtype=old.dtype)
            new[: self.capacity] = old
            if attr in ("frag_ids", "best_tid", "best_gid", "frag_class", "locus_id"):
                new[self.capacity :] = -1
            else:
                new[self.capacity :] = 0
            setattr(self, attr, new)
        self.capacity = new_cap

    def _grow(self) -> None:
        """Double capacity when full."""
        self._grow_to(max(self.capacity * 2, CONSTANTS.resources.annotation_table_min_growth))

    @property
    def size(self) -> int:
        return self._size


def _tag_strings(values) -> list[str]:
    """IDs and names as the BAM's string tags: a missing or empty one is ``.``."""
    return [v if isinstance(v, str) and v else "." for v in values]


#: The counts the annotation pass reports, each named and counted as the scan counts it, and the digest of its read
#: pairs: every one must equal the scan's, or the pass stamped what the scan did not decide.
SAME_AS_THE_SCAN = (
    "n_read_pairs",
    "n_excluded",
    "n_withheld_unpaired",
    "n_withheld_unannotated_splice",
    "n_withheld_off_index",
    "n_withheld_chimeric",
    "n_withheld_unpaired_missing_mate",
    "n_withheld_unpaired_reference",
    "n_withheld_unpaired_strand",
    "n_withheld_unpaired_outward",
    "n_withheld_unpaired_length",
    "n_excluded_by_option",
    "n_placed_unique",
    "n_placed_multimapping",
    "n_records_used",
    "n_records_unused",
    "n_records_unplaced",
    "filtered",
    "template_digest",
)


def write_annotated_bam(
    bam_path: str,
    output_path: str,
    annotations: AnnotationTable,
    index,
    *,
    stats,
    skip_duplicates: bool,
    include_multimap: bool,
    sj_spec: str,
    locus_id_per_transcript: np.ndarray | None = None,
) -> dict:
    """The annotated-BAM pass: stamp every record with assignment tags.

    Re-reads the input BAM in the scan's name-sorted order, looks up each read-name group's annotation by
    ``frag_id`` and writes the tags, in C++ (``BamAnnotationWriter``).

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
    stats : PipelineStats
        The scan's counts. The pass re-reads the BAM and re-judges every read pair, and must see what the scan saw:
        every count in :data:`SAME_AS_THE_SCAN`, the digest of its read pairs among them, and every record written. On
        any difference the annotated BAM is deleted and the run stops.
    skip_duplicates : bool
        Must match the value used in Pass 1.
    include_multimap : bool
        Must match the value used in Pass 1.
    sj_spec : str
        The splice-junction tag spec Pass 1's scanner read with, resolved (``"XS"``, ``"ts"``,
        ``"XS,ts"``, ``"ts,XS"`` or ``"none"``), never ``"auto"``.
    locus_id_per_transcript : np.ndarray or None
        int32 array mapping transcript index → locus_id. If provided,
        used to populate the ZL (locus_id) BAM tag.

    Returns
    -------
    dict
        Summary statistics of the annotation pass.
    """
    from .native import BamAnnotationWriter

    # the scan's own resolver
    writer = BamAnnotationWriter(
        index.resolver,
        sj_spec,
        skip_duplicates,
        include_multimap,
    )

    n = annotations.size

    # ZL: the fragment's own locus from the EM scatter, which carries in-locus gDNA winners. A row still -1, such as
    # one from the deterministic path that bypasses the EM, falls back to its winning transcript's locus; intergenic
    # gDNA stays -1.
    if locus_id_per_transcript is not None:
        best_tid = annotations.best_tid[:n]
        from_tid = np.where(
            best_tid >= 0,
            locus_id_per_transcript[np.clip(best_tid, 0, None)],
            -1,
        )
        current = annotations.locus_id[:n]
        annotations.locus_id[:n] = np.where(current >= 0, current, from_tid)

    # The writer takes strings; an unnamed gene (NaN in the index) is written `.`
    t_ids = _tag_strings(index.t_df["t_id"].values)
    g_ids = _tag_strings(index.g_df["g_id"].values)
    g_names = _tag_strings(index.g_df["g_name"].values)

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
        annotations.locus_id[:n],
        n,
        t_ids,
        g_ids,
        g_names,
    )

    differ = {
        key: (summary[key], getattr(stats, key))
        for key in SAME_AS_THE_SCAN
        if summary[key] != getattr(stats, key)
    }
    if summary["n_records_written"] != stats.total:
        differ["records written"] = (summary["n_records_written"], stats.total)
    if differ:
        Path(output_path).unlink(missing_ok=True)
        raise RuntimeError(
            "The annotation pass did not see what the scan saw, so its BAM is deleted: "
            + "; ".join(
                f"{key} {seen:,} against the scan's {scanned:,}"
                for key, (seen, scanned) in differ.items()
            )
            + ". The BAM must not change between the two passes."
        )
    logger.info(
        f"[ANNOTATE] Wrote {summary['n_records_written']:,} records to {output_path}: "
        f"{summary['n_read_pairs']:,} read pairs, as the scan judged them"
    )

    return summary
