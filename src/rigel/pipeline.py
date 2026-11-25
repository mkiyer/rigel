"""The quantification pipeline: one BAM scan, the second pass, calibration, then the per-locus EM.

The native scanner reads the BAM once. ``scan_and_buffer`` resolves every fragment against the index, trains the
strand models, tallies the calibration payload and buffers the resolved fragments (``FragmentBuffer``).
``run_pipeline`` drains the held fragments, fits the fragment-length models on the drained payload and calibrates.
``quant_from_buffer`` scores the buffer, builds the loci, assembles each locus's gDNA prior from the calibration
result (``calibration.priors.assemble_priors``) and runs the native locus EM. This module only orchestrates: the
scorer is ``scoring.py``, the router ``routing.py``, the loci ``locus.py`` and the EM ``em_solver.cpp``.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, fields
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from .annotate import (
    AF_GDNA_EM,
    AF_GDNA_INTERGENIC,
    AF_UNRESOLVED,
    winner_flags,
)
from .buffer import FragmentBuffer, FinalizedChunk
from .config import (
    CONSTANTS,
    EMConfig,
    PipelineConfig,
    BamScanConfig,
    FragmentScoringConfig,
    TranscriptGeometry,
)
from .estimator import AbundanceEstimator
from .frag_length_model import FragmentLengthModel
from .index import TranscriptIndex
from .native import EM_WINNER_GDNA, EM_WINNER_INTERGENIC
from .native import BamScanner as _NativeBamScanner
from .native import UnpairedRecordError
from .routing import FragmentRouter
from .scoring import build_scorer, sequence_log_penalties
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from collections.abc import Callable

    import pandas as pd

    from .annotate import AnnotationTable
    from .calibration import CalibrationResult
    from .calibration.fl import FLModels
    from .calibration.region_arrays import RegionArrays
    from .calibration.shared_length_law import LengthTables
    from .scan_payload import AccumulatorPayload
    from .scored_fragments import ScoredFragments

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline result
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """Complete output of the pipeline."""

    stats: PipelineStats
    strand_models: StrandModels
    estimator: AbundanceEstimator
    pipeline_config: PipelineConfig
    calibration: CalibrationResult
    # Library-wide FL distributions actually built + used for calibration (global / RNA / gDNA); the
    # empirical views feed the summary QC report.
    fl_models: FLModels
    # The EM's two length tables and the decision behind them (`calibration.shared_length_law`).
    length_tables: LengthTables
    # The level a junction's capture price is capped at: calibration's saturated level at its capture weight, raised as
    # far as RNA's own spliced fragments support (`calibration.capture_eff_length.fit_junction_cap`).
    junction_cap: float
    # Per-region genome-wide gDNA track (ref/start/end + gDNA mass/density/frac),
    # built from the calibration result; feeds the report's genome track + bedGraph.
    calibration_track: pd.DataFrame


def _log_sj_tags(stats: PipelineStats, scan: BamScanConfig) -> None:
    """Say what the junction-strand tag gave on the library's spliced records, and warn when it gave none a strand.

    The tag is read on every record, so the counts are the whole library's. No share of untagged records is a
    threshold: only a library none of whose spliced records carries a usable tag trains no junction strand.
    """
    spliced = stats.n_sj_tag_tagged + stats.n_sj_tag_untagged + stats.n_sj_tag_conflicting

    def share(n: int) -> float:
        return 100.0 * n / spliced if spliced else 0.0

    logger.info(
        "[SCAN] splice-junction strand tag %s: %d spliced records, %d with a strand (%.2f %%), %d without one "
        "(%.2f %%), %d whose XS and ts disagree, read as ambiguous (%.2f %%); %d unreadable tag values",
        scan.sj_tag_spec,
        spliced,
        stats.n_sj_tag_tagged,
        share(stats.n_sj_tag_tagged),
        stats.n_sj_tag_untagged,
        share(stats.n_sj_tag_untagged),
        stats.n_sj_tag_conflicting,
        share(stats.n_sj_tag_conflicting),
        stats.n_sj_tag_invalid,
    )
    if spliced > 0 and stats.n_sj_tag_tagged == 0:
        logger.warning(
            "No splice-junction strand tag (XS, or minimap2's ts) is read on this library's spliced reads, so no "
            "junction trains the strand model and the library is quantified as unstranded. If it is stranded, "
            "align with the tag (STAR: --outSAMstrandField intronMotif) or name it with --sj-strand-tag."
        )


def _log_record_tags(stats: PipelineStats) -> None:
    """Say what NM and NH were on the library's records, and warn of an NM that is absent, unreadable or disputed.

    A record without NM reads as 0 edits, so it can outscore an alignment whose mismatches are known. A record without
    NH is ordinary for an aligner that marks its other alignments secondary instead.
    """
    logger.info(
        "[SCAN] record tags: NM absent on %d records, unreadable on %d, disputed by %d alignments; NH absent on %d "
        "records, unreadable on %d",
        stats.n_nm_tag_missing,
        stats.n_nm_tag_invalid,
        stats.n_nm_conflicting,
        stats.n_nh_tag_missing,
        stats.n_nh_tag_invalid,
    )
    if stats.n_nm_tag_missing or stats.n_nm_tag_invalid:
        logger.warning(
            "%d records carry no NM and %d an NM that is not an integer; each is read as 0 edits, which can favor it "
            "over an alignment whose mismatches are known. STAR writes NM with --outSAMattributes NH HI AS nM NM.",
            stats.n_nm_tag_missing,
            stats.n_nm_tag_invalid,
        )
    if stats.n_nm_conflicting:
        logger.warning(
            "%d alignments are written in several records that disagree on NM; each is read with the largest.",
            stats.n_nm_conflicting,
        )


def _log_hit_verdicts(stats: PipelineStats) -> None:
    """Say how many pairings of an R1 alignment with an R2 alignment were tried and dropped, and why, and how many
    templates that left with no usable hit."""
    rejected = {
        "one mate only": stats.n_hits_rejected_missing_mate,
        "two references": stats.n_hits_rejected_reference,
        "one strand": stats.n_hits_rejected_strand,
        "facing outward": stats.n_hits_rejected_outward,
        "too long on every path": stats.n_hits_rejected_length,
        "junctions that contradict": stats.n_hits_rejected_contradictory,
        "a novel junction": stats.n_hits_rejected_novel_junction,
        "junctions no one transcript over both mates holds": stats.n_hits_rejected_unsupported,
        "off the index": stats.n_hits_rejected_off_index,
        "a mate aligned in pieces (chimeric)": stats.n_hits_rejected_chimeric,
    }
    if stats.n_nh_above_reported:
        logger.warning(
            "%d templates carry an NH above the records reported for a mate, so their multiplicity is "
            "read from the alignments present; STAR reports every alignment with --outSAMmultNmax -1",
            stats.n_nh_above_reported,
        )
    if not any(rejected.values()):
        return
    logger.info(
        "[SCAN] %d pairings dropped (%s); templates left with no usable hit, withheld: %d unpaired, %d unannotated "
        "splices, %d off the index, %d chimeric",
        sum(rejected.values()),
        ", ".join(f"{why} {n}" for why, n in rejected.items() if n),
        stats.n_withheld_unpaired,
        stats.n_withheld_unannotated_splice,
        stats.n_withheld_off_index,
        stats.n_withheld_chimeric,
    )


def check_scan_accounting(stats: PipelineStats) -> None:
    """Every read pair has one outcome and every record one status, so the scan's counts reconcile by construction; a
    difference is a defect in the scan, and the run stops rather than report counts that do not add up."""
    withheld = (
        stats.n_withheld_unpaired
        + stats.n_withheld_unannotated_splice
        + stats.n_withheld_off_index
        + stats.n_withheld_chimeric
    )
    unpaired = (
        stats.n_withheld_unpaired_missing_mate
        + stats.n_withheld_unpaired_reference
        + stats.n_withheld_unpaired_strand
        + stats.n_withheld_unpaired_outward
        + stats.n_withheld_unpaired_length
    )
    placed = stats.n_placed_unique + stats.n_placed_multimapping
    pairs = stats.n_excluded + withheld + stats.n_excluded_by_option + placed
    records = (
        stats.filtered + stats.n_records_used + stats.n_records_unused + stats.n_records_unplaced
    )
    broken = [
        what
        for what, holds in (
            (
                f"{stats.n_read_pairs:,} read pairs but {pairs:,} outcomes",
                pairs == stats.n_read_pairs,
            ),
            (f"{stats.total:,} records but {records:,} statuses", records == stats.total),
            (
                f"{stats.n_withheld_unpaired:,} unpaired but {unpaired:,} reasons",
                unpaired == stats.n_withheld_unpaired,
            ),
            (f"{stats.n_fragments:,} fragments but {placed:,} placed", placed == stats.n_fragments),
        )
        if not holds
    ]
    if broken:
        raise RuntimeError(f"The scan's counts do not reconcile: {'; '.join(broken)}.")
    logger.info(
        "[SCAN] %s read pairs: %s placed (%s unique, %s multimapping), %s withheld, %s excluded (every record "
        "unmapped, QC-failed or a duplicate)%s",
        f"{stats.n_read_pairs:,}",
        f"{placed:,}",
        f"{stats.n_placed_unique:,}",
        f"{stats.n_placed_multimapping:,}",
        f"{withheld:,}",
        f"{stats.n_excluded:,}",
        f", {stats.n_excluded_by_option:,} excluded by --no-include-multimap"
        if stats.n_excluded_by_option
        else "",
    )


def _apply_scan_stats(stats: PipelineStats, stats_dict: dict) -> None:
    """Copy the C++ scan's counters onto ``stats``, each key onto the field of its name.

    A key with no field is refused, so a counter the scan exports cannot be dropped on the way to ``summary.json``.
    """
    unknown = sorted(set(stats_dict) - {f.name for f in fields(stats)})
    if unknown:
        raise KeyError(f"scan counters with no PipelineStats field: {unknown}")
    for key, value in stats_dict.items():
        setattr(stats, key, value)


def warn_if_calibration_strand_unidentifiable(strand_models: StrandModels) -> None:
    """Warn when calibration's strand decision (`strand_balance.strand_discriminability` on the sense fraction
    `strand_balance.fit_strand_balance` fits) reads the library as unstranded. Silent with no spliced observation,
    which is warned about where it is counted."""
    from .calibration.strand_balance import fit_strand_balance

    balance = fit_strand_balance(strand_models)
    if balance.fallback_used or balance.live:
        return
    logger.warning(
        "[CAL] Calibration reads the library as unstranded (n_spliced_obs=%d, "
        "rna_sense_frac=%.6f), so its strand channel is off. If this library was expected to be "
        "stranded, inspect its splice evidence.",
        balance.n_observations,
        balance.rna_sense_frac,
    )


class BamInputError(ValueError):
    """A BAM the scan would read wrong, refused before its records are trusted; the message names the remedy."""


def check_annotated_bam_path(bam_path, annotated_bam_path) -> None:
    """Refuse an annotated BAM that is the input BAM, under its own name or another (a hard or symbolic link): the
    annotation pass opens its output after its input, so it would truncate the BAM it is reading, and its check against
    the scan would then delete it."""
    if annotated_bam_path is None:
        return
    target, bam = Path(annotated_bam_path), Path(bam_path)
    if target.exists() and bam.exists() and target.samefile(bam):
        raise BamInputError(
            f"The annotated BAM {annotated_bam_path} is the input BAM {bam_path}: writing it would destroy the "
            "input. Name another path."
        )


def _check_bam_header(bam_path: str, index: TranscriptIndex) -> None:
    """Refuse, from its header alone, a BAM the scan would read wrong.

    A coordinate-sorted BAM puts a pair's mates apart, and the scan groups adjacent records by read name, so every mate
    would be its own fragment. A BAM none of whose references the index holds was aligned elsewhere, and a reference
    it holds at another length is another assembly. References beyond the index (decoys, alternates, spike-ins) are
    allowed and named in a warning: the scan withholds the reads aligned there, counted in no pool.
    """
    import pysam

    try:
        with pysam.AlignmentFile(bam_path, "r", check_sq=False) as f:
            header = f.header.to_dict()
    # pysam refuses a truncated BAM (no BGZF end-of-file marker) with an OSError, and a file that is no SAM or BAM at all
    # with a ValueError
    except (OSError, ValueError) as exc:
        raise BamInputError(f"{bam_path} cannot be read as a BAM: {exc}") from exc
    if header.get("HD", {}).get("SO") == "coordinate":
        raise BamInputError(
            f"{bam_path} is sorted by coordinate (@HD SO:coordinate), so a pair's mates are not adjacent and each "
            "would be read as its own fragment; sort it by read name (samtools sort -n) or collate it "
            "(samtools collate) first."
        )
    bam_refs = {sq["SN"]: int(sq["LN"]) for sq in header.get("SQ", [])}
    shared = [name for name in bam_refs if name in index.ref_lengths]
    some = CONSTANTS.qc.examples
    if bam_refs and not shared:
        raise BamInputError(
            f"None of {bam_path}'s {len(bam_refs):,} references is in the index ({some(bam_refs)}, against "
            f"{some(index.ref_lengths)}): the BAM and the index name their references differently, or were built on "
            "different genomes."
        )
    relengthed = [name for name in shared if bam_refs[name] != int(index.ref_lengths[name])]
    if relengthed:
        raise BamInputError(
            f"{len(relengthed):,} of {bam_path}'s references have another length in the index "
            f"({some(f'{n}: {bam_refs[n]:,} against {int(index.ref_lengths[n]):,}' for n in relengthed)}): "
            "the BAM and the index were built on different assemblies."
        )
    beyond = [name for name in bam_refs if name not in index.ref_lengths]
    if beyond:
        logger.warning(
            "%s of the BAM's %s references are not in the index (%s); the reads aligned there are withheld, "
            "counted in no pool (read_pairs.withheld.off_index).",
            f"{len(beyond):,}",
            f"{len(bam_refs):,}",
            some(beyond),
        )


def scan_and_buffer(
    bam_path: str,
    index: TranscriptIndex,
    scan: "BamScanConfig",
) -> tuple[PipelineStats, StrandModels, FragmentBuffer, "AccumulatorPayload"]:
    """One native BAM pass: parse, resolve, train the strand models and buffer the fragments.

    Parameters
    ----------
    bam_path : str
        Path to the name-sorted / collated BAM file.
    index : TranscriptIndex
        The loaded reference index.
    scan : BamScanConfig
        BAM scanning and buffering configuration.

    Returns
    -------
    tuple
        ``(stats, strand_models, buffer, calibration_payload)``. A scan that fails cleans its buffer up, so it leaves
        no spilled chunk behind.

    No fragment-length model comes out of here: every length law is built from the payload by
    :func:`rigel.calibration.fl.build_fl_models`, so length is measured by one rule, the accumulator's ``L``.
    """
    _check_bam_header(bam_path, index)
    buffer = FragmentBuffer(
        max_memory_bytes=scan.buffer_size_bytes,
        spill_dir=None if scan.spill_dir is None else Path(scan.spill_dir),
    )
    try:
        stats, strand_models, payload = _scan(bam_path, index, scan, buffer)
    except BaseException:
        buffer.cleanup()
        raise
    return stats, strand_models, buffer, payload


def _warn_lacking_mates(stats: PipelineStats) -> None:
    """Warn, with the count, of read pairs lacking a mate their own records call mapped. Their mates sit elsewhere in
    the BAM, which then is not grouped by read name (sorted by coordinate, its header not saying so), or were
    removed record by record; each such read pair is read without its mate."""
    if stats.n_read_pairs_lacking_mate:
        logger.warning(
            "%s of %s read pairs lack a mate their records call mapped (FLAG 0x8 unset): the mates are elsewhere in "
            "the BAM or were removed, and each such pair is read without its mate. A BAM sorted by coordinate whose "
            "header does not say so reads this way; sort it by read name (samtools sort -n) or collate it "
            "(samtools collate).",
            f"{stats.n_read_pairs_lacking_mate:,}",
            f"{stats.n_read_pairs:,}",
        )


def _scan(
    bam_path: str, index: TranscriptIndex, scan: "BamScanConfig", buffer: FragmentBuffer
) -> tuple[PipelineStats, StrandModels, "AccumulatorPayload"]:
    """The scan itself, its fragments streamed into *buffer*: the counters, checked to reconcile, the strand models
    and the calibration payload."""
    from .scan_payload import AccumulatorPayload

    logger.info("[START] Native C++ BAM scan → resolve + train + buffer")
    scanner = _native_scanner(index, scan)
    n_scan, n_bgzf = scan.resolved_scan_threads()
    if scan.bgzf_threads is not None and n_bgzf != scan.bgzf_threads:
        logger.info(
            "[scan] Capped BGZF decompression threads from %d to %d to fit "
            "within total thread budget %d.",
            scan.bgzf_threads,
            n_bgzf,
            scan.resolved_total_threads(),
        )
    try:
        result = scanner.scan(
            bam_path,
            # each chunk a zero-copy dict from C++
            chunk_callback=lambda raw: buffer.inject_chunk(FinalizedChunk.from_raw(raw)),
            chunk_size=scan.fragments_per_chunk,
            n_workers=n_scan,
            n_decomp_threads=n_bgzf,
            qname_batch_size=scan.read_name_batch_size,
        )
    except UnpairedRecordError as exc:
        raise BamInputError(f"{bam_path}: {exc}") from exc

    stats = PipelineStats()
    _apply_scan_stats(stats, result["stats"])
    check_scan_accounting(stats)
    _warn_lacking_mates(stats)
    strand_models = _strand_models(result["strand_observations"], stats)
    logger.info(
        f"[DONE] Native scan: {stats.n_fragments:,} fragments → "
        f"{stats.n_same_strand:,} same-strand, "
        f"{stats.n_ambig_strand:,} ambig-strand, "
        f"{stats.n_intergenic:,} intergenic, "
        f"{buffer.total_fragments:,} buffered "
        f"({buffer.memory_bytes / 1024**2:.1f} MB in memory, "
        f"{buffer.n_spilled} chunks spilled), "
        f"{stats.n_multimapper_groups:,} multimapper molecules"
    )
    strand_models.log_summary()
    # The calibration payload, keyed by `graph_hash` (regions and boundaries), because the payload's sj axis is
    # meaningless against a different sj CSR
    return (
        stats,
        strand_models,
        AccumulatorPayload.from_scan_result(result, graph_hash=index.graph_hash),
    )


def resolver_policy(scan: "BamScanConfig") -> dict:
    """The settings a pass resolves under, the scan's own, given to each native pass rather than set on the index's
    shared resolver: the pair rule's fragment-length limit (no pair whose molecule is longer is one) and the
    intron-overhang tolerance K (which pairs are one molecule, the candidates' lengths, calibration's tally and which
    unspliced fragments become SPLICED_IMPLICIT)."""
    return {
        "max_fragment_length": int(scan.max_frag_length),
        "intron_overhang_tolerance": int(scan.intron_overhang_tolerance),
    }


def _native_scanner(index: TranscriptIndex, scan: "BamScanConfig") -> _NativeBamScanner:
    """The native scanner over the index's resolver, the splice graph installed. The index installs what it knows
    (the transcript strands among it); the scan sets what is a setting."""
    scanner = _NativeBamScanner(
        index.resolver,
        scan.sj_tag_spec,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
        **resolver_policy(scan),
    )
    # The splice graph's region partition, so the scan tallies per-region evidence
    _wire_calibration_regions(scanner, index, scan.max_frag_length)
    return scanner


def _strand_models(strand_observations: dict, stats: PipelineStats) -> StrandModels:
    """The strand models, built once; the spliced 2×2 is the marginal of the per-sj strand table. One strand-qualified
    fragment credits one sj, so the 2×2 is exactly the table's marginal and both halves of the RNA strand
    Beta-Binomial are fitted on one population. The C++ counts the fragments independently of the table, so this is a
    real cross-check."""
    strand_models = StrandModels.from_scan(strand_observations)
    n_credited = strand_models.sj_table.n_observations
    if n_credited != stats.n_strand_trained:
        raise RuntimeError(
            f"SJ strand table credited {n_credited:,} observations but the scanner "
            f"qualified {stats.n_strand_trained:,} fragments; one qualified fragment must "
            "credit exactly one sj."
        )
    return strand_models


def library_fl_models(payload, index: TranscriptIndex) -> "FLModels":
    """The fragment-length models from one payload, exactly as production builds them.

    All three histograms come from the payload, so there is one definition of length; a fragment enters a pool
    only when exactly one hypothesis survived. The RNA pool is read at its own sj opportunity (using an annotated sj
    is a length-dependent selection), and the gDNA pools, none of them pure, are deconvolved by the two-pool
    contrast over the partition the scanner deposits into (`calibration.fl`). Every gDNA opportunity is read over
    the index's observable genome (`gdna_opportunity.observable_gdna_opportunity`). RNA's law is its own counts';
    gDNA's is empirical-Bayes shrunk toward the all-fragment anchor.

    The one home for this call: the second pass scores on pass one's payload and calibration fits on the drained
    one, and an instrument must read the same models.
    """
    from .calibration.fl import build_fl_models
    from .calibration.gdna_opportunity import observable_gdna_opportunity
    from .calibration.region_arrays import RegionArrays
    from .calibration.sj_opportunity import crossing_probability_from_index

    regions = RegionArrays.from_index(index)
    max_length = int(payload.max_length)
    return build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_length),
        gdna_opportunity=observable_gdna_opportunity(regions, max_length),
        regions=regions,
    )


def length_tables(
    fl_models: "FLModels", calibration: "CalibrationResult", region_arrays: "RegionArrays"
) -> "LengthTables":
    """The EM's two length tables, exactly as production decides them: the two measured laws and their sampling
    variances, against RNA's captured fraction under calibration's capture efficiencies
    (`calibration.shared_length_law`, `calibration.capture_eff_length.rna_captured_fraction`). The scorer and RNA's
    normaliser read the result, and so must any instrument that stands in for either."""
    from .calibration.capture_eff_length import rna_captured_fraction
    from .calibration.shared_length_law import shared_length_law

    return shared_length_law(
        fl_models.rna_pmf,
        fl_models.gdna_pmf,
        fl_models.rna_pmf_var,
        fl_models.gdna_pmf_var,
        rna_captured_fraction(calibration, region_arrays),
        rna_measured=fl_models.n_rna > 0.0,
    )


def drain_side_buffer(
    payload, index: TranscriptIndex, strand_models, *, seed: int, lift: dict | None = None
):
    """The second pass: score the held fragments, draw one hypothesis each, and re-deposit them.

    Pass 1 holds two kinds of fragment whole: one whose unsequenced gap has more than one surviving explanation
    (mostly the long ones, since a longer gap admits more hypotheses), and a multimapping fragment more than one of
    whose placements admits a surviving (placement, hypothesis). Nothing is tallied for them until this runs.

    It runs between the scan and calibration: every input the score needs comes from pass 1 alone, so calibration
    runs once, on the complete tally. The scorer uses pass one's models and calibration the drained tally's. Fit
    once, score once, drain once: the confident set is biased short, so feeding the drained law back into the score
    would favour the shorter, more-spliced path, and the loop could run away.

    ``lift`` is an out-parameter for an origin-split oracle (``lift_choices`` in `scripts/design/_oracle.py`):
    the drain conditions on the whole tally, so the oracle replays the whole's drawn choices inside each partition,
    and needs the choices, the undrained payload and the two index-derived arrays `drain` takes, which exist only
    here. Like ``calibrate(diagnostics=)``, it is inert in production. An empty side buffer leaves ``lift`` untouched.
    """
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import build_sj_arrays, build_region_partition_arrays
    from .second_pass import choose_hypotheses, drain, score_held_fragments

    held = payload.deferred.n_fragments
    if held == 0:
        # Not an error and not a no-op to hide: a library with no annotated intron in any mate gap is a
        # real state. Left undrained, so `payload.drain is None` continues to mean "pass one only".
        logger.info("[SP2] nothing held; the side buffer is empty and the drain is skipped")
        return payload

    start = time.perf_counter()
    _, _, region_types = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)
    scores = score_held_fragments(
        payload,
        # Pass one's length models, de-tilted as calibration's are: the scorer weighs a candidate length by
        # `f(L)`, and the raw sj-tilted RNA pool would favour the longer hypothesis
        fl_models=library_fl_models(payload, index),
        # `P(align_strand agrees | RNA)`: near 0 on an R1-antisense (dUTP) library, where disagreement is likely
        rna_sense_frac=strand_models.p_r1_sense,
        regions=RegionArrays.from_index(index),
        sj=sj,
    )
    choices = choose_hypotheses(scores, payload, seed=seed)
    drained = drain(payload, choices, region_types=region_types, sj=sj)
    if lift is not None:
        # ``undrained`` is the payload as it entered here, not ``drained``: the drained bank is empty by
        # design ("after it nothing is held"), so it cannot supply `lift_choices`' key pool.
        lift.update(undrained=payload, choices=choices, region_types=region_types, sj=sj)
    report = drained.drain
    logger.info(
        "[SP2] drained %d held fragments (%d held on their placements) in %.1f s: %d deposited, "
        "%d dropped (%d chose a hypothesis that implies no gap intron, %d one that implies a gap "
        "intron); %d had candidates tied for the lead",
        report.offered,
        report.offered_multimapper,
        time.perf_counter() - start,
        report.deposited,
        report.dropped,
        report.chose_no_implied_intron,
        report.chose_implied_intron,
        scores.n_undecided,
    )
    return drained


def _wire_calibration_regions(
    scanner,
    index: TranscriptIndex,
    max_frag_length: int,
) -> None:
    """Install the index's splice graph into a native BamScanner: the region partition, then the sj.

    :func:`~rigel.calibration.splice_graph.build_region_partition_arrays` flattens the per-reference
    partition into the ``(region_bounds, ref_region_bound_offsets, n_refs, region_types, max_length)`` ABI, and
    :func:`~rigel.calibration.splice_graph.build_sj_arrays` builds the sj CSR keyed by the
    flat region_bound index. Both are derived from ``index.regions_df``/``index.edges_df`` and ordered to match
    ``index.ref_names``, which is the resolver's reference-id space.

    Both calls are required: ``set_regions`` refuses to run twice, and ``scan`` refuses to run without ``set_sj``,
    because a missing sj table would read every intron as unannotated and leave the spliced banks silently empty.
    """
    from .calibration.splice_graph import build_sj_arrays, build_region_partition_arrays

    region_bounds, ref_region_bound_offsets, region_types = build_region_partition_arrays(index)
    n_refs = len(index.ref_names)
    scanner.set_regions(
        np.ascontiguousarray(region_bounds, dtype=np.int64),
        np.ascontiguousarray(ref_region_bound_offsets, dtype=np.int64),
        n_refs,
        np.ascontiguousarray(region_types, dtype=np.uint8),
        int(max_frag_length),
    )
    # ``edge_row`` is not passed: it is a join key back to ``index.edges_df``, not the sj id (the CSR slot is), and
    # indexing a sj bank with it would write past the end of the array
    sj = build_sj_arrays(index)
    scanner.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int64),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int64),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )


# ---------------------------------------------------------------------------
# Quantify from buffer (the locus EM)
# ---------------------------------------------------------------------------


def length_models(tables: "LengthTables") -> tuple[FragmentLengthModel, FragmentLengthModel]:
    """``(rna_fl, gdna_fl)``: the scorer's two length models, on the EM's tables (:func:`length_tables`)."""
    max_size = int(tables.rna.size) - 1
    return (
        FragmentLengthModel.from_pmf(tables.rna, max_size),
        FragmentLengthModel.from_pmf(tables.gdna, max_size),
    )


def transcript_geometry(
    index: TranscriptIndex,
    rna_fl,
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    tables: "LengthTables",
    n_threads: int = 1,
) -> tuple[float, TranscriptGeometry]:
    """``(junction_cap, geometry)``: the cap on a junction's capture price and the transcripts' effective lengths, read
    from one frame of the index's pieces and cuts (``calibration.capture_eff_length.TranscriptFrame``). The frame is
    ten million rows on a whole genome and is dropped on return, before a fragment is scored.

    The cap is the saturated level calibration publishes, raised as far as RNA's spliced fragments support
    (``fit_junction_cap``). ``effective_lengths_em`` are capture-contracted by the rule the gDNA component's length
    takes too (each transcript's share of every region, boundary and junction its fragments deposit on, at that
    object's capture efficiency, a junction's capped at the cap), so every component competes on one footing under
    capture. ``effective_lengths``, TPM's, stay the full FL-marginal length over the observable starts. Where every
    capture efficiency is exactly 1 (no landscape, or a reading of 1 everywhere) the two are bit-identical, whatever
    the length tables (``calibration/capture_eff_length.py``). The frame's shares run on ``n_threads`` threads.
    """
    if rna_fl.n_observations <= 0:
        raise ValueError(
            "the library carries no RNA fragment-length observation, so no transcript has an "
            "effective length; nothing can be quantified."
        )
    from .calibration.capture_eff_length import (
        TranscriptFrame,
        fit_junction_cap,
        observable_fl_eff_lengths,
        transcript_capture_eff_lengths,
    )

    frame = TranscriptFrame(index, region_arrays, tables, n_threads)
    junction_cap = fit_junction_cap(calibration, frame)
    effective_lengths = observable_fl_eff_lengths(frame, rna_fl)
    effective_lengths_em = transcript_capture_eff_lengths(
        calibration, frame, effective_lengths, junction_cap
    )
    return junction_cap, TranscriptGeometry(
        effective_lengths=effective_lengths,
        effective_lengths_em=effective_lengths_em,
    )


def score_fragments(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    rna_fl,
    gdna_fl,
    stats: PipelineStats,
    estimator: AbundanceEstimator,
    scoring: "FragmentScoringConfig",
    log_every: int,
    annotations,
    place_efficiencies,
) -> "ScoredFragments":
    """Build the scorer, route the buffer, and return ScoredFragments. A hit's sequence is priced at the mismatch
    rate the scan read on the library's unique unspliced hits."""
    mismatch_log_penalty, clip_log_penalty = sequence_log_penalties(stats.mismatch_rate)
    scorer = build_scorer(
        strand_models,
        rna_fl,
        gdna_fl,
        index,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        clip_log_penalty=clip_log_penalty,
        pruning_min_posterior=scoring.pruning_min_posterior,
    )
    router = FragmentRouter(
        scorer,
        estimator,
        stats,
        index,
        annotations=annotations,
        place_efficiencies=place_efficiencies,
    )
    # routing consumes the buffer; the router and the scorer go with this frame
    return router.route(buffer, log_every)


def _assign_locus_ids(estimator: AbundanceEstimator, multi_loci: list) -> None:
    """Stamp ``multi_locus_id`` onto every transcript on the estimator — the ``locus_id`` column of the
    output tables and the annotated BAM's ZL tag."""
    for locus in multi_loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.multi_locus_id


def populate_em_annotations(
    batch_parts,
    out_winner_tid,
    out_winner_post,
    out_n_candidates,
    annotations,
    index,
):
    """Populate AnnotationTable from EM assignment output and partition metadata."""
    t_to_g = index.t_to_g_arr
    is_synth_arr = index.t_df["is_synthetic"].values.astype(bool)

    # Concatenate per-locus metadata in same order as C++ output
    frag_ids = np.concatenate([p.frag_ids for p in batch_parts])
    frag_class = np.concatenate([p.frag_class for p in batch_parts])
    # True per-fragment locus-subproblem id (so in-locus gDNA winners keep their
    # locus in the ZL tag, instead of collapsing to -1 via the winning transcript).
    locus_ids = np.concatenate(
        [np.full(p.frag_ids.shape[0], p.locus_id, dtype=np.int32) for p in batch_parts]
    )

    n = len(frag_ids)
    best_tid = np.asarray(out_winner_tid, dtype=np.int32)
    posteriors = np.asarray(out_winner_post, dtype=np.float32)
    n_cand = np.asarray(out_n_candidates, dtype=np.int32)

    # Gene index: map valid transcript winners to gene
    valid_t = best_tid >= 0
    best_gid = np.full(n, -1, dtype=np.int32)
    best_gid[valid_t] = t_to_g[best_tid[valid_t]]

    # ZF assignment flags bitfield
    tx_flags = np.full(n, AF_UNRESOLVED, dtype=np.uint8)
    tx_flags[best_tid == EM_WINNER_GDNA] = AF_GDNA_EM
    tx_flags[best_tid == EM_WINNER_INTERGENIC] = AF_GDNA_INTERGENIC
    if valid_t.any():
        winner_tids = best_tid[valid_t]
        tx_flags[valid_t] = winner_flags(is_synth_arr[winner_tids])

    # Clean sentinels: the locus's gDNA and the intergenic term -> -1
    best_tid_clean = best_tid.copy()
    best_tid_clean[(best_tid == EM_WINNER_GDNA) | (best_tid == EM_WINNER_INTERGENIC)] = -1

    annotations.add_batch(
        frag_ids=frag_ids,
        best_tids=best_tid_clean,
        best_gids=best_gid,
        tx_flags=tx_flags,
        posteriors=posteriors,
        frag_classes=frag_class.view(np.int8),
        n_candidates=n_cand,
        locus_ids=locus_ids,
    )


def em_pseudocounts(
    gdna_count: np.ndarray, fragment_count: np.ndarray, n_gdna_eligible: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The locus EM's two pseudocounts, ``(gdna_prior_count, rna_prior_count)``, per locus.

    The odds are calibration's gDNA share of the locus. The EM's theta is each component's share of every fragment
    in the locus (the E-step's units and the deterministic spliced fragments, which skip the E-step but enter every
    M-step), and the grouped prior update is exact EM for ``P_g log theta_g + P_R log(1 − theta_g)``, which leaves
    the EM's own answer in place iff ``P_g : P_R = G : (N − G)``. ``gdna_count`` and ``fragment_count`` are
    calibration's two counts of the locus (``LocusPriors``): its gDNA fragments, and every fragment it deposited on
    the same objects. Both hold a multimapper once, where the second pass drew it, so the share is read over one set
    of fragments and applied to the whole locus; a fragment calibration never saw carries the locus's share, not RNA.

    The strength is one pseudocount per unit that has a gDNA candidate; it is chosen, not derived. A locus with no
    gDNA candidate, or none of whose fragments calibration counted, gets no prior. A gDNA count above the fragment
    count, beyond rounding, is refused: the two are one tally, so it is an accounting error, and clipping it would
    hand the locus a share of one and RNA no pseudocount at all.
    """
    gdna_count = np.asarray(gdna_count, dtype=np.float64)
    fragment_count = np.asarray(fragment_count, dtype=np.float64)
    n_gdna_eligible = np.asarray(n_gdna_eligible, dtype=np.float64)
    seen = fragment_count > 0
    share = np.divide(gdna_count, fragment_count, out=np.zeros_like(gdna_count), where=seen)
    if np.any(share > 1.0 + CONSTANTS.calibration.fraction_tolerance):
        i = int(np.argmax(share))
        raise ValueError(
            f"a locus's gDNA count ({float(gdna_count[i])!r}) exceeds its count of every fragment "
            f"({float(fragment_count[i])!r}); the two must be read over the same fragments."
        )
    share = np.clip(share, 0.0, 1.0)
    strength = np.where(seen, n_gdna_eligible, 0.0)
    gdna_prior = strength * share
    return gdna_prior, strength - gdna_prior


def run_locus_em_partitioned(
    estimator: AbundanceEstimator,
    partitions: dict,
    multi_loci: list,
    index: TranscriptIndex,
    gdna_count: np.ndarray,
    fragment_count: np.ndarray,
    gdna_eff_len: np.ndarray,
    *,
    em_config: EMConfig,
    intergenic_rate: float = 0.0,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> None:
    """Run the per-locus batch EM and record lean per-locus results.

    Hands every ``LocusPartition`` to the native EM, which reads each array by its name, and runs the whole batch in
    one ``run_batch_locus_em_partitioned`` call (the solver runs its own threads), and appends a per-locus dict to
    ``estimator.locus_results`` for ``get_loci_df``. Calibration enters as its per-locus gDNA count, which
    :func:`em_pseudocounts` reads against its count of every fragment it deposited on the locus. The loci table's
    ``has_gdna`` is structural: some unspliced unit carries a finite gDNA log-likelihood.
    """
    if not multi_loci:
        return

    parts = [partitions[loc.multi_locus_id] for loc in multi_loci]
    ids = [loc.multi_locus_id for loc in multi_loci]
    n_gdna_eligible = np.array(
        [np.count_nonzero((p.gdna_excluded == 0) & np.isfinite(p.gdna_log_liks)) for p in parts],
        dtype=np.float64,
    )
    gdna_prior, rna_prior = em_pseudocounts(
        np.asarray(gdna_count, dtype=np.float64)[ids],
        np.asarray(fragment_count, dtype=np.float64)[ids],
        n_gdna_eligible,
    )
    gdna_prior = np.ascontiguousarray(gdna_prior)
    rna_prior = np.ascontiguousarray(rna_prior)

    em_result = estimator.run_batch_locus_em_partitioned(
        parts,
        [loc.transcript_indices for loc in multi_loci],
        gdna_prior,
        rna_prior_count=rna_prior,
        gdna_eff_len=np.ascontiguousarray(gdna_eff_len[ids], dtype=np.float64),
        intergenic_rate=intergenic_rate,
        emit_locus_stats=emit_locus_stats,
        emit_assignments=annotations is not None,
    )
    total_gdna_em, rna_arr, gdna_arr = em_result[0], em_result[1], em_result[2]
    if annotations is not None:
        populate_em_annotations(parts, em_result[3], em_result[4], em_result[5], annotations, index)
    _record_locus_results(
        estimator,
        multi_loci,
        index,
        rna_arr,
        gdna_arr,
        gdna_prior,
        rna_prior,
        has_gdna=n_gdna_eligible > 0,
        gdna_eff_len=gdna_eff_len,
    )

    estimator._gdna_em_total = total_gdna_em
    n_units = sum(len(loc.unit_indices) for loc in multi_loci)
    logger.info(
        "[DONE] Per-locus EM: %d loci, %d ambiguous fragments, gDNA EM=%.0f",
        len(multi_loci),
        n_units,
        total_gdna_em,
    )


def _record_locus_results(
    estimator: AbundanceEstimator,
    multi_loci: list,
    index: TranscriptIndex,
    rna_arr: np.ndarray,
    gdna_arr: np.ndarray,
    gdna_prior: np.ndarray,
    rna_prior: np.ndarray,
    *,
    has_gdna: np.ndarray,
    gdna_eff_len: np.ndarray,
) -> None:
    """Append each locus's lean result to ``estimator.locus_results``, the rows ``get_loci_df`` builds the loci table
    from; ``n_genes`` counts the annotated genes among the locus's transcripts."""
    t_to_g = index.t_to_g_arr
    is_synth_g = index.g_df["is_synthetic"].to_numpy()
    for i, loc in enumerate(multi_loci):
        lid = loc.multi_locus_id
        gene_set = {
            int(t_to_g[int(t)])
            for t in loc.transcript_indices
            if not is_synth_g[int(t_to_g[int(t)])]
        }
        estimator.locus_results.append(
            {
                "locus_id": lid,
                "locus_span_bp": loc.span_bp,
                "n_genes": len(gene_set),
                "n_em_fragments": len(loc.unit_indices),
                "rna_total": float(rna_arr[i]),
                "count_gdna": float(gdna_arr[i]),
                "count_gdna_prior": float(gdna_prior[i]),
                "count_rna_prior": float(rna_prior[i]),
                "has_gdna": bool(has_gdna[i]),
                "gdna_em_effective_length": float(gdna_eff_len[lid]),
            }
        )


def quant_from_buffer(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    tables: "LengthTables",
    region_arrays: "RegionArrays",
    stats: PipelineStats,
    calibration: "CalibrationResult",
    junction_cap: float,
    geometry: TranscriptGeometry,
    *,
    em_config: EMConfig,
    scoring: FragmentScoringConfig,
    log_every: int = BamScanConfig.log_every,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> AbundanceEstimator:
    """Quantify buffered fragments: calibration prior → per-locus EM.

    Scores the buffer with the EM's two length tables (:func:`length_tables`), prices every junction under one cap
    (``junction_cap``, :func:`transcript_geometry`) in a multimapper's place weights as the transcripts' lengths
    (``geometry``) were priced, builds connected-component loci, turns the per-region ``CalibrationResult``
    into each locus's gDNA count, fragment count and gDNA length (``assemble_priors``), partitions
    the global CSR, and runs the per-locus EM, whose pseudocounts read the gDNA count
    against the fragment count (``em_pseudocounts``). Returns the populated
    ``AbundanceEstimator``.
    """
    from .calibration.capture_eff_length import place_efficiencies
    from .calibration.priors import assemble_priors, intergenic_gdna_rate
    from .locus import build_multi_loci
    from .locus_partition import partition_and_free

    # The scorer's two tables: one shared law unless the data show a difference capture cannot make. The
    # transcripts' geometry and RNA's normaliser read the same tables; calibration and the priors keep `gdna_pmf`.
    rna_fl, gdna_fl = length_models(tables)

    estimator = AbundanceEstimator(
        index.num_transcripts,
        em_config=em_config,
        geometry=geometry,
        is_synthetic=index.t_df["is_synthetic"].values,
    )

    em_data = score_fragments(
        buffer,
        index,
        strand_models,
        rna_fl,
        gdna_fl,
        stats,
        estimator,
        scoring,
        log_every,
        annotations,
        # each multimapper's place weighs its candidates by its capture, a junction's at the lengths' own cap
        place_efficiencies(calibration, region_arrays, index, junction_cap),
    )

    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)

    if em_data.n_units == 0 or not multi_loci:
        return estimator

    priors = assemble_priors(calibration, region_arrays, multi_loci)
    partitions = partition_and_free(em_data, multi_loci)
    run_locus_em_partitioned(
        estimator,
        partitions,
        multi_loci,
        index,
        priors.gdna_count,
        priors.fragment_count,
        priors.gdna_eff_len,
        em_config=em_config,
        intergenic_rate=intergenic_gdna_rate(calibration, region_arrays),
        annotations=annotations,
        emit_locus_stats=emit_locus_stats,
    )
    del partitions
    return estimator


# ---------------------------------------------------------------------------
# Full pipeline orchestration
# ---------------------------------------------------------------------------


def run_pipeline(
    bam_path,
    index: TranscriptIndex,
    config: PipelineConfig | None = None,
    *,
    before_annotated_bam: "Callable[[PipelineResult], None] | None" = None,
) -> PipelineResult:
    """Run the complete quantification pipeline with locus-level EM.

    Reads the BAM once with the native scanner, drains the held fragments, calibrates, then quantifies with the
    per-locus EM. The fragment buffer is cleaned up on every path, a failure included.

    Parameters
    ----------
    bam_path : str or Path
        Path to the name-sorted / collated BAM file.
    index : TranscriptIndex
        The loaded reference index.
    config : PipelineConfig or None
        Complete pipeline configuration.  If ``None``, uses defaults.
    before_annotated_bam : callable or None
        Called with the finished result before the optional annotated-BAM pass, so a caller can persist the
        quantification before a pass that re-reads the whole BAM; the CLI writes its tables here.

    Returns
    -------
    PipelineResult
    """
    if config is None:
        config = PipelineConfig()
    bam_path = str(bam_path)
    check_annotated_bam_path(bam_path, config.annotated_bam_path)
    scan = config.scan
    stats, strand_models, buffer, calibration_payload = scan_and_buffer(bam_path, index, scan)
    with buffer:
        _log_sj_tags(stats, scan)
        _log_record_tags(stats)
        _log_hit_verdicts(stats)
        # The second pass drains the side buffer before calibration
        calibration_payload = drain_side_buffer(
            calibration_payload, index, strand_models, seed=config.second_pass_seed
        )
        calibration, region_arrays, fl_models = _calibrate(
            calibration_payload, index, strand_models, config.calibration
        )
        # The EM's two length tables: one shared law unless the data show a difference capture cannot make
        tables = length_tables(fl_models, calibration, region_arrays)
        _log_length_tables(tables)
        # The cap on a junction's capture price (the saturated level, raised as far as RNA's spliced fragments
        # support) and the transcripts' lengths priced under it
        junction_cap, geometry = transcript_geometry(
            index,
            length_models(tables)[0],
            calibration,
            region_arrays,
            tables,
            n_threads=config.calibration.n_threads,
        )
        _log_junction_cap(calibration.gdna_capture_saturation, junction_cap)
        annotations = None if config.annotated_bam_path is None else _annotation_table(buffer)
        estimator = quant_from_buffer(
            buffer,
            index,
            strand_models,
            tables,
            region_arrays,
            stats,
            calibration,
            junction_cap,
            geometry,
            em_config=config.em,
            scoring=config.scoring,
            log_every=scan.log_every,
            annotations=annotations,
            emit_locus_stats=config.emit_locus_stats,
        )

    # The genome-wide gDNA track for the QC report and the bedGraph: the per-region calibration solution
    from .calibration.track import build_gdna_track

    result = PipelineResult(
        stats=stats,
        strand_models=strand_models,
        estimator=estimator,
        pipeline_config=config,
        calibration=calibration,
        fl_models=fl_models,
        length_tables=tables,
        junction_cap=junction_cap,
        calibration_track=build_gdna_track(calibration, region_arrays, index.ref_names),
    )
    if before_annotated_bam is not None:
        before_annotated_bam(result)

    # The annotated-BAM pass: re-read the BAM and write it with assignment tags (opt-in)
    if annotations is not None:
        from .annotate import write_annotated_bam

        write_annotated_bam(
            bam_path,
            str(config.annotated_bam_path),
            annotations,
            index,
            stats=stats,
            scan=scan,
            locus_id_per_transcript=estimator.locus_id_per_transcript,
        )
    return result


def _calibrate(
    payload, index: TranscriptIndex, strand_models: StrandModels, calibration_config
) -> tuple["CalibrationResult", "RegionArrays", "FLModels"]:
    """Calibration on the drained tally: deconvolve every region and boundary into gDNA / RNA+ / RNA−. Returns the
    result with the region frame it was solved over and the fragment-length models it read, refit on that tally."""
    from .calibration import calibrate
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays

    warn_if_calibration_strand_unidentifiable(strand_models)
    confidence = CONSTANTS.qc.strand_ci_confidence
    logger.info(
        "[CAL] Strand trainer: n_spliced_obs=%d  ss_est=%.6f  ε_CI(%g%%)=%.4g",
        strand_models.n_observations,
        strand_models.strand_specificity,
        100.0 * confidence,
        strand_models.strand_specificity_ci_epsilon(),
    )

    # No alignment check here: calibrate() runs it (CalibrationSubstrate.from_payload)
    region_arrays = RegionArrays.from_index(index)
    boundary_flags = build_boundary_flags_array(index)
    # The sj axis, in the accumulator's sj slot order: where each sj attaches, its transcript strand and its exonic
    # reach either side. The calibrator places each sj as a factor on its two endpoint regions, never as a message
    # channel, because every sj closes a loop and the graph is not a tree.
    sj = build_sj_geometry_arrays(index)

    # The two component length models the calibrator's effective lengths need, refit on the drained tally
    fl_models = library_fl_models(payload, index)
    bins = np.arange(fl_models.gdna_pmf.size, dtype=np.float64)
    logger.info(
        "[CAL] FL models: n_gdna_pool=%.0f gdna_fl_mean=%.1f rna_fl_mean=%.1f global_fl_mean=%.1f",
        fl_models.n_gdna,
        float(np.dot(bins, fl_models.gdna_pmf)),
        float(np.dot(bins, fl_models.rna_pmf)),
        float(np.dot(bins, fl_models.global_pmf)),
    )

    calibration = calibrate(
        payload=payload,
        region_arrays=region_arrays,
        strand_model=strand_models,
        gdna_fl_pmf=fl_models.gdna_pmf,
        rna_fl_pmf=fl_models.rna_pmf,
        config=calibration_config,
        sj=sj,
        boundary_flags=boundary_flags,
    )
    # Fragments, not object incidences: each axis (regions, boundaries, sj) is converted by its own population's
    # mass per crossing before the sum, since raw banks count a fragment once per boundary it crosses and once
    # per sj it uses
    logger.info(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gDNA_fragments=%.0f RNA_fragments=%.0f (sj incidences %.0f)",
        calibration.n_regions,
        calibration.n_boundaries,
        calibration.n_sj,
        calibration.gdna_density_global,
        calibration.rna_sense_frac,
        calibration.library_gdna_fragments,
        calibration.library_rna_fragments,
        float(calibration.count_rna_sj.sum()),
    )
    return calibration, region_arrays, fl_models


def _log_length_tables(tables: "LengthTables") -> None:
    decision = tables.decision
    bins = np.arange(tables.rna.size, dtype=np.float64)
    logger.info(
        "[FL] length tables: rna_captured=%.4f captured=%.4f contrast=%.4f bound=%.4f noise=%.3g weight=%.3f "
        "rna_mean=%.1f gdna_table_mean=%.1f",
        decision.rna_captured_fraction,
        decision.captured,
        decision.contrast,
        decision.bound,
        decision.noise,
        decision.weight,
        float(np.dot(bins, tables.rna)),
        float(np.dot(bins, tables.gdna)),
    )


def _log_junction_cap(saturation: float, junction_cap: float) -> None:
    """Say where the junction cap sits, and whether RNA's spliced fragments raised it."""
    if junction_cap == saturation:
        logger.info(
            "junction cap: the saturated level at the capture weight, %.4g× the bulk (the spliced counts raise nothing)",
            saturation,
        )
    else:
        logger.info(
            "junction cap: %.4g× the bulk, raised by RNA's spliced fragments from the saturated level at the capture "
            "weight, %.4g×",
            junction_cap,
            saturation,
        )


def _annotation_table(buffer: FragmentBuffer) -> "AnnotationTable":
    """The empty per-fragment table the annotated-BAM pass reads, sized to the buffered fragments."""
    from .annotate import AnnotationTable

    return AnnotationTable.create(
        capacity=max(
            buffer.total_fragments + CONSTANTS.resources.annotation_table_padding,
            CONSTANTS.resources.annotation_table_min_capacity,
        )
    )
