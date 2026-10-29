"""
rigel.pipeline — Single-pass BAM quantification pipeline with locus-level EM.

Architecture
------------
The pipeline reads the BAM file **once** via the C++ native BAM scanner
(``_bam_impl``) and processes fragments in two logical stages:

**BAM Scan** (``scan_and_buffer``): Parse all fragments in C++ via htslib,
resolve each against the reference index, train the strand models from
unique-mapper fragments, tally the calibration payload, and buffer all
resolved fragments into a memory-efficient columnar buffer (``FragmentBuffer``).

**The second pass and calibration** (``run_pipeline``): drain the held
fragments, fit the fragment-length models on the drained payload, and
calibrate.

**Quantification** (``quant_from_buffer``): score buffered fragments,
construct MultiLoci, assemble the per-locus gDNA/RNA split prior from the
calibration result (``calibration.priors.assemble_priors``), partition the
CSR data, and run locus-level native EM.

``scoring.py`` assembles the native scorer's parameters.  Locus
construction lives in ``locus.py``; the EM's warm start is computed inside
the native solver (``em_solver.cpp``).  The router that scores the buffer
and builds the CSR lives in ``routing.py``.  This module is a thin orchestrator.
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
from .buffer import FRAG_MULTIMAPPER, FragmentBuffer, _FinalizedChunk
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
from .native import EM_WINNER_INTERGENIC
from .native import BamScanner as _NativeBamScanner
from .routing import FragmentRouter
from .scoring import build_scorer, sequence_log_penalties
from .splice import SpliceType
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from collections.abc import Callable

    from .annotate import AnnotationTable
    from .calibration import CalibrationResult
    from .calibration.fl import FLModels
    from .calibration.region_arrays import RegionArrays
    from .calibration.shared_length_law import LengthTables
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
    pipeline_config: "PipelineConfig" = None
    calibration: "CalibrationResult" = None
    # Library-wide FL distributions actually built + used for calibration (global / RNA / gDNA); the
    # empirical views feed the summary QC report.
    fl_models: "FLModels" = None
    # The EM's two length tables and the decision behind them (`calibration.shared_length_law`).
    length_tables: "LengthTables" = None
    # Per-region genome-wide gDNA track (ref/start/end + gDNA mass/density/frac),
    # built from the calibration result; feeds the report's genome track + bedGraph.
    calibration_track: "object" = None


def _sj_tag_to_spec(sj_strand_tag) -> str:
    """Convert BamScanConfig.sj_strand_tag to the string spec for BamScanner."""
    if isinstance(sj_strand_tag, str):
        return sj_strand_tag if sj_strand_tag else "none"
    if isinstance(sj_strand_tag, (list, tuple)):
        return ",".join(sj_strand_tag) if sj_strand_tag else "none"
    return "none"


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
        _sj_tag_to_spec(scan.sj_strand_tag),
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


def _check_scan_accounting(stats: PipelineStats) -> None:
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


def _warn_if_calibration_strand_unidentifiable(strand_models: StrandModels) -> None:
    """Warn when calibration's strand decision (`region_init.strand_discriminability`, on the spliced
    sense fraction `strand_balance.fit_strand_balance` fits) reads the library as unstranded. With no
    spliced observation it says nothing: the balance falls back to κ = ½ and calibration reads the library as
    unstranded, and the empty spliced census is warned about where it is counted."""
    from .calibration.region_init import strand_discriminability
    from .calibration.strand_balance import fit_strand_balance

    balance = fit_strand_balance(strand_models)
    if balance.fallback_used:
        return
    if strand_discriminability(balance.rna_sense_frac, balance.n_observations) > 0.0:
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
    except (
        OSError
    ) as exc:  # pysam refuses a BAM missing its BGZF end-of-file marker: a truncated file
        raise BamInputError(f"{bam_path} cannot be read: {exc}") from exc
    if header.get("HD", {}).get("SO") == "coordinate":
        raise BamInputError(
            f"{bam_path} is sorted by coordinate (@HD SO:coordinate), so a pair's mates are not adjacent and each "
            "would be read as its own fragment; sort it by read name (samtools sort -n) or collate it "
            "(samtools collate) first."
        )
    bam_refs = {sq["SN"]: int(sq["LN"]) for sq in header.get("SQ", [])}
    shared = [name for name in bam_refs if name in index.ref_lengths]

    def some(names):
        names = list(names)
        return ", ".join(names[:3]) + (f" and {len(names) - 3:,} more" if len(names) > 3 else "")

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
) -> tuple[
    PipelineStats,
    StrandModels,
    FragmentBuffer,
    "object",  # AccumulatorPayload
]:
    """Single-pass C++ BAM scan: resolve, train models, buffer — all in one pass.

    The entire BAM parsing, fragment construction, overlap resolution,
    model training and columnar buffering happens in C++ via htslib.

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
        ``(stats, strand_models, buffer, calibration_payload)``.

    No fragment-length model comes out of here. Every fragment-length distribution the tool uses is
    built from the payload by :func:`rigel.calibration.fl.build_fl_models`, so length is measured by
    one rule — the accumulator's ``L`` — over one stated population.
    """
    _check_bam_header(bam_path, index)
    stats = PipelineStats()
    buffer = FragmentBuffer(
        max_memory_bytes=scan.buffer_size_bytes,
        spill_dir=None if scan.spill_dir is None else Path(scan.spill_dir),
    )
    logger.info("[START] Native C++ BAM scan → resolve + train + buffer")

    # The index installs what it knows (the transcript strands among it); the scan sets what is a setting.
    resolve_ctx = index.resolver

    # Intron-overhang tolerance K (bp): how many bases an aligned read end may lie inside an annotated
    # intron, at the edge of an unsequenced gap or at a fragment's outer end. The pair rule applies it
    # along a witness transcript, and the resolver when it enumerates a fragment's gap hypotheses: an
    # intron in an unaligned gap is cut whole within the gap widened by K, and a transcript is
    # contradicted only by an intron holding more than K aligned bases. So K decides which pairs are
    # one molecule, the candidates' lengths, calibration's tally (the scan-cache key includes it) and
    # which unspliced fragments are relabelled SPLICED_IMPLICIT.
    resolve_ctx.set_intron_overhang_tolerance(int(scan.intron_overhang_tolerance))
    # The pair rule's fragment-length limit: no pair whose molecule is longer is one.
    resolve_ctx.set_max_fragment_length(int(scan.max_frag_length))

    # Create the native scanner
    sj_spec = _sj_tag_to_spec(scan.sj_strand_tag)
    scanner = _NativeBamScanner(
        resolve_ctx,
        sj_spec,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
    )

    # Calibration region wiring: install the splice graph's region partition into the scanner so
    # per-region evidence is collected during the scan. The loader refuses an index without the graph.
    _wire_calibration_regions(
        scanner,
        index,
        scan.max_frag_length,
    )

    # Streaming chunk callback — receives zero-copy dict from C++
    def _on_chunk(raw: dict) -> None:
        chunk = _FinalizedChunk.from_raw(raw)
        buffer.inject_chunk(chunk)

    # Execute the full BAM scan in C++ with streaming chunk output
    n_scan, n_bgzf = scan.resolved_scan_threads()
    if scan.bgzf_threads is not None and n_bgzf != scan.bgzf_threads:
        logger.info(
            "[scan] Capped BGZF decompression threads from %d to %d to fit "
            "within total thread budget %d.",
            scan.bgzf_threads,
            n_bgzf,
            scan.resolved_total_threads(),
        )
    result = scanner.scan(
        bam_path,
        chunk_callback=_on_chunk,
        chunk_size=scan.fragments_per_chunk,
        n_workers=n_scan,
        n_decomp_threads=n_bgzf,
        qname_batch_size=scan.read_name_batch_size,
    )

    # Replay stats
    _apply_scan_stats(stats, result["stats"])
    _check_scan_accounting(stats)

    # Build the strand models from the scanner's observations. Immutable and built
    # once — the spliced 2×2 is the marginal of its per-sj strand table.
    strand_models = StrandModels.from_scan(result["strand_observations"])
    # ONE strand-qualified fragment credits ONE sj. That is what makes the 2×2
    # exactly the table's marginal, so both halves of the RNA strand Beta-Binomial are
    # fitted on one population. The
    # C++ counts the fragments independently of the table it builds, so this is a real
    # cross-check, and it is the invariant to break if the crediting rule ever changes.
    n_credited = strand_models.sj_table.n_observations
    if n_credited != stats.n_strand_trained:
        raise RuntimeError(
            f"SJ strand table credited {n_credited:,} observations but the scanner "
            f"qualified {stats.n_strand_trained:,} fragments; one qualified fragment must "
            "credit exactly one sj."
        )

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

    # The calibration payload: the C++ scanner's accumulator tally, typed.
    from .scan_payload import AccumulatorPayload

    # The provenance covers regions AND boundaries. The payload is boundary-keyed by
    # construction — its sj axis is meaningless against a different sj CSR — whereas
    # `partition_hash` covers `regions.feather` only, deliberately.
    calibration_payload = AccumulatorPayload.from_scan_result(result, graph_hash=index.graph_hash)

    return stats, strand_models, buffer, calibration_payload


def library_fl_models(payload, index: TranscriptIndex) -> "FLModels":
    """The fragment-length models from one payload, exactly as production builds them.

    All three histograms come from the PAYLOAD — one frame, one definition of length. A fragment
    enters a pool when exactly one hypothesis survived, so its ``L`` is not in doubt. The RNA pool is read at
    its own sj opportunity ("used an annotated sj" is a length-dependent selection), and
    the gDNA pools — none of them pure — are deconvolved by the two-pool contrast, which needs the
    partition the scanner deposits into (`calibration.fl`). Every gDNA opportunity is read over the
    index's observable genome: the pools' and the regions' (`gdna_opportunity.observable_gdna_opportunity`,
    the frame's `region_arrays.observable_pieces`). RNA's law is its own counts', the unconditional anchor only
    the law an empty pool returns; gDNA's is smooth-EB shrunk toward that anchor.

    The one home for this call: the second pass scores on pass one's payload, calibration fits on the
    drained one, and an instrument that measured a different length model than production's would be
    measuring a tool that does not ship.
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


def _drain_side_buffer(
    payload, index: TranscriptIndex, strand_models, *, seed: int, _lift: dict | None = None
):
    """THE SECOND PASS: score the held fragments, draw one hypothesis each, re-deposit.

    Pass 1 holds two kinds of fragment WHOLE: one whose unsequenced gap has more than one surviving
    explanation — systematically the long ones, because a longer gap admits more hypotheses — and a
    multimapping fragment more than one of whose placements admits a surviving (placement,
    hypothesis) pair. Nothing is lost, but nothing is tallied either until this runs.

    It runs HERE, between the scan and calibration, and that is a structural decision. Every input the
    score needs — the densities at each object, the fragment-length models, the strand model — comes
    from pass 1 alone, so no calibration output is required. That is what lets calibration run exactly
    once, on the complete tally, instead of once before the drain and again after it.

    The models the SCORER uses are pass one's; the models CALIBRATION uses are the drained tally's. Fit
    once, score once, drain once, stop: the confident set is biased short, so feeding the drained anchor
    back into the score would prefer the shorter — that is, the more-spliced — path, and that loop can
    run away.

    ``_lift`` is the out-parameter an origin-split ORACLE needs. The drain conditions on the WHOLE
    tally, so partitions drained independently do not sum to the whole drained.
    `second_pass.lift_choices` repairs that by replaying the
    whole's already-drawn choices inside each partition — which needs the choices, the undrained whole
    they were drawn on, and the two index-derived arrays `drain` takes. All four exist only inside this
    function, so it publishes them into ``_lift`` rather than letting a caller re-derive them and drift.
    Same convention as ``calibrate(_debug=)`` /
    ``solve_chain(_capture=)``, and inert in production, where nobody passes it. An empty side buffer
    leaves ``_lift`` UNTOUCHED — the early return below is the "nothing was drained" signal on this
    path too.
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
        # Pass one's length models, de-tilted the way calibration's are: the scorer weighs a candidate
        # length by `f(L)`, so the raw sj-tilted RNA pool would make it prefer the longer hypothesis for
        # the same reason that pool is long. Calibration refits on the drained tally — fit once, score
        # once, drain once (the docstring above).
        fl_models=library_fl_models(payload, index),
        # `P(align_strand agrees | RNA)`. On an R1-antisense (dUTP) library it is near 0, so
        # DISAGREEMENT is the likely case.
        rna_sense_frac=strand_models.p_r1_sense,
        regions=RegionArrays.from_index(index),
        sj=sj,
    )
    choices = choose_hypotheses(scores, payload, seed=seed)
    drained = drain(payload, choices, region_types=region_types, sj=sj)
    if _lift is not None:
        # ``undrained`` is the payload as it entered here, NOT ``drained``: the drained bank is empty by
        # design ("after it nothing is held"), so it cannot supply `lift_choices`' key pool.
        _lift.update(undrained=payload, choices=choices, region_types=region_types, sj=sj)
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

    Two calls, and both are required. ``set_regions`` refuses to run twice, which is why the sj are
    separate; and ``scan`` refuses to run if the second call is missing, because a missing sj table is
    invisible — every observed intron would simply read as unannotated, so every sj boundary and both
    spliced banks would come back empty from a scan that looked perfectly well-formed.
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
    # ``edge_row`` is deliberately NOT passed. It is a join key back to ``index.edges_df``, not the
    # sj-boundary id — that IS the CSR slot — and using it to index a sj bank would write far past
    # the end of the array.
    sj = build_sj_arrays(index)
    scanner.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int64),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int64),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )


# ---------------------------------------------------------------------------
# Quantify from buffer (locus-level EM, no global EM)
# ---------------------------------------------------------------------------


def _setup_geometry_and_estimator(
    index: TranscriptIndex,
    rna_fl,
    em_config: EMConfig,
    calibration: "CalibrationResult",
    region_arrays: "RegionArrays",
    tables: "LengthTables",
) -> AbundanceEstimator:
    """Compute the transcripts' effective lengths and create the AbundanceEstimator.

    The EM effective lengths (``effective_lengths_em``) are capture-contracted by the rule the gDNA component's
    length takes too — each transcript's conserved share of every region, boundary and junction its fragments
    deposit on, at that object's capture efficiency carried to the EM's two length tables — so every component
    competes on one contracted footing under capture. The output ``effective_lengths`` (used for TPM) stay the full
    FL-marginal length over the template's observable starts (a fragment touching the index's unobservable mask is
    never observed). No reference (capture off, or no gDNA) ⇒ ``effective_lengths_em`` == ``effective_lengths``
    (bit-identical). See ``calibration/capture_eff_length.py``.
    """
    if rna_fl.n_observations <= 0:
        raise ValueError(
            "the library carries no RNA fragment-length observation, so no transcript has an "
            "effective length; nothing can be quantified."
        )
    from .calibration.capture_eff_length import (
        observable_fl_eff_lengths,
        transcript_capture_eff_lengths,
    )

    effective_lengths = observable_fl_eff_lengths(index, region_arrays, rna_fl)

    effective_lengths_em = transcript_capture_eff_lengths(
        calibration, region_arrays, index, effective_lengths, tables
    )

    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        effective_lengths_em=effective_lengths_em,
    )

    estimator = AbundanceEstimator(
        index.num_transcripts,
        em_config=em_config,
        geometry=geometry,
        is_synthetic=index.t_df["is_synthetic"].values,
    )
    return estimator


def _score_fragments(
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
    place_efficiencies=None,
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


def _populate_em_annotations(
    batch_parts,
    out_winner_tid,
    out_winner_post,
    out_n_candidates,
    annotations,
    index,
):
    """Populate AnnotationTable from EM assignment output and partition metadata."""
    if out_winner_tid is None:
        return

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
    n_cand = np.asarray(out_n_candidates, dtype=np.int16)

    # Gene index: map valid transcript winners to gene
    valid_t = best_tid >= 0
    best_gid = np.full(n, -1, dtype=np.int32)
    best_gid[valid_t] = t_to_g[best_tid[valid_t]]

    # ZF assignment flags bitfield
    tx_flags = np.full(n, AF_UNRESOLVED, dtype=np.uint8)
    tx_flags[best_tid == -2] = AF_GDNA_EM
    tx_flags[best_tid == EM_WINNER_INTERGENIC] = AF_GDNA_INTERGENIC
    if valid_t.any():
        winner_tids = best_tid[valid_t]
        tx_flags[valid_t] = winner_flags(is_synth_arr[winner_tids])

    # Clean sentinels: -2 (the locus's gDNA) and the intergenic term -> -1
    best_tid_clean = best_tid.copy()
    best_tid_clean[(best_tid == -2) | (best_tid == EM_WINNER_INTERGENIC)] = -1

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
    gdna_count: np.ndarray, n_calibrated: np.ndarray, n_gdna_eligible: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The locus EM's two pseudocounts, ``(gdna_prior_count, rna_prior_count)``, per locus.

    THE ODDS ARE CALIBRATION'S gDNA SHARE OF THE LOCUS. The EM's theta is each component's share of every
    fragment in the locus — the E-step's units and the deterministic spliced fragments, which skip the
    E-step but enter every M-step — and the grouped prior update is exact EM for
    ``P_g log theta_g + P_R log(1 − theta_g)``. It leaves the EM's own answer in place iff
    ``P_g : P_R = G : (N − G)``. Calibration counts ``gdna_count`` from its own sample of the locus, and
    ``n_calibrated`` counts that sample's unique fragments (no multimapper, no fragment whose every junction the
    blacklist rejected), so the share is read over that sample and applied to the whole locus: which fragments
    either stage calls spliced never enters, and a fragment calibration never saw carries the locus's share, not
    RNA. Calibration also deposits each multimapper at the placement its second pass draws, and
    ``n_calibrated`` does not count those, so in a locus holding multimappers the count covers more fragments
    than the share is read over.

    THE STRENGTH is one pseudocount per unit that has a gDNA candidate, the shipped total before the odds
    were corrected; it is not derived (``tests/test_em_pseudocounts.py``). A share above one is a share of
    one, and a locus with no gDNA candidate, or none of whose fragments calibration counted, gets no prior.
    """
    gdna_count = np.asarray(gdna_count, dtype=np.float64)
    n_calibrated = np.asarray(n_calibrated, dtype=np.float64)
    n_gdna_eligible = np.asarray(n_gdna_eligible, dtype=np.float64)
    seen = n_calibrated > 0
    share = np.clip(
        np.divide(gdna_count, n_calibrated, out=np.zeros_like(gdna_count), where=seen), 0.0, 1.0
    )
    strength = np.where(seen, n_gdna_eligible, 0.0)
    gdna_prior = strength * share
    return gdna_prior, strength - gdna_prior


def _run_locus_em_partitioned(
    estimator: AbundanceEstimator,
    partitions: dict,
    multi_loci: list,
    index: TranscriptIndex,
    gdna_count: np.ndarray,
    gdna_eff_len: np.ndarray,
    *,
    em_config: EMConfig,
    intergenic_rate: float = 0.0,
    rna_prior_weight: np.ndarray | None = None,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> None:
    """Run the per-locus batch EM and record lean per-locus results.

    Packs each ``LocusPartition`` into the C++ 7-tuple, runs the whole batch in
    one ``run_batch_locus_em_partitioned`` call (the solver runs its own
    ``std::thread`` workers), and appends a per-locus dict to ``estimator.locus_results`` for
    ``get_loci_df``. Calibration enters as its per-locus gDNA count, which :func:`em_pseudocounts` reads
    against the unique part of the locus calibration counts from: its non-multimapper units that the blacklist
    did not reject, plus the deterministic fragments of its transcripts (the multimappers calibration deposits
    at their drawn placements are not counted). The
    loci table's ``has_gdna`` reports the structural eligibility (any unspliced unit carrying a finite
    gDNA log-lik — the rule the C++ extractor applies).

    ``rna_prior_weight`` is the one PER-TRANSCRIPT array here and it is passed FLAT, while every
    other prior array is subscripted by ``ids`` into per-locus order. That asymmetry is the point: the
    RNA prior's TOTAL is a property of the locus, its ALLOCATION is a property of the transcripts, and
    the C++ remaps the latter itself when it builds each sub-problem. ⛔ Subscripting it by ``ids``
    would silently read transcript rows at locus positions.
    """
    if not multi_loci:
        return

    parts = [partitions[loc.multi_locus_id] for loc in multi_loci]
    partition_tuples = [
        (
            p.offsets,
            p.t_indices,
            p.log_liks,
            p.coverage_weights,
            p.count_cols,
            p.is_spliced,
            p.gdna_log_liks,
            p.intergenic_log_liks,
        )
        for p in parts
    ]
    locus_t_lists = [loc.transcript_indices for loc in multi_loci]
    ids = [loc.multi_locus_id for loc in multi_loci]
    g_eff = np.ascontiguousarray(gdna_eff_len[ids], dtype=np.float64)
    n_gdna_eligible = np.array(
        [np.count_nonzero((p.is_spliced == 0) & np.isfinite(p.gdna_log_liks)) for p in parts],
        dtype=np.float64,
    )
    # The unique part of calibration's sample: the deterministic fragments (unique by construction) and every
    # unit it deposits that is not a multimapper or a fragment whose every junction the blacklist rejected. The
    # multimappers it deposits at their drawn placements are not counted.
    deterministic = estimator.unambig_counts.sum(axis=1)
    n_calibrated = np.array(
        [
            np.count_nonzero(
                (p.frag_class != FRAG_MULTIMAPPER) & (p.splice_type != SpliceType.SPLICE_ARTIFACT)
            )
            + deterministic[np.asarray(loc.transcript_indices, np.int64)].sum()
            for p, loc in zip(parts, multi_loci)
        ],
        dtype=np.float64,
    )
    gdna_prior, rna_prior = em_pseudocounts(
        np.asarray(gdna_count, dtype=np.float64)[ids], n_calibrated, n_gdna_eligible
    )
    gdna_prior = np.ascontiguousarray(gdna_prior)
    rna_prior = np.ascontiguousarray(rna_prior)
    has_gdna = n_gdna_eligible > 0

    em_result = estimator.run_batch_locus_em_partitioned(
        partition_tuples,
        locus_t_lists,
        gdna_prior,
        rna_prior_count=rna_prior,
        # FLAT and per-TRANSCRIPT — deliberately not `[ids]`. If this boundary is dropped the
        # parameter is accepted and silently ignored, and every allocation, however extreme, produces
        # byte-identical output. A fire counter cannot see that — it counts nonzero weights, not
        # weights the solver read; only an arm that injects a maximally WRONG allocation and checks
        # that the result MOVED can.
        rna_prior_weight=rna_prior_weight,
        gdna_eff_len=g_eff,
        intergenic_rate=intergenic_rate,
        em_iterations=em_config.iterations,
        em_convergence_delta=em_config.convergence_delta,
        emit_locus_stats=emit_locus_stats,
        emit_assignments=annotations is not None,
    )
    total_gdna_em, rna_arr, gdna_arr = em_result[0], em_result[1], em_result[2]
    if annotations is not None and len(em_result) > 3:
        _populate_em_annotations(
            parts, em_result[3], em_result[4], em_result[5], annotations, index
        )

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

    estimator._gdna_em_total = total_gdna_em
    n_units = sum(len(loc.unit_indices) for loc in multi_loci)
    logger.info(
        "[DONE] Per-locus EM: %d loci, %d ambiguous fragments, gDNA EM=%.0f",
        len(multi_loci),
        n_units,
        total_gdna_em,
    )


def quant_from_buffer(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    tables: "LengthTables",
    region_arrays: "RegionArrays",
    stats: PipelineStats,
    calibration: "CalibrationResult",
    *,
    em_config: EMConfig | None = None,
    scoring: FragmentScoringConfig | None = None,
    log_every: int = BamScanConfig.log_every,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> AbundanceEstimator:
    """Quantify buffered fragments: calibration prior → per-locus EM.

    Scores the buffer with the EM's two length tables (:func:`length_tables`),
    builds connected-component loci, turns the per-region ``CalibrationResult``
    into each locus's gDNA count and gDNA length (``assemble_priors``), partitions
    the global CSR, and runs the per-locus EM, whose pseudocounts read that count
    against the locus's own fragments (``em_pseudocounts``). Returns the populated
    ``AbundanceEstimator``.
    """
    from .calibration.capture_eff_length import place_efficiencies
    from .calibration.priors import assemble_priors, intergenic_gdna_rate
    from .locus import build_multi_loci
    from .locus_partition import partition_and_free

    em_config = em_config or EMConfig()
    scoring_cfg = scoring or FragmentScoringConfig()

    # The scorer's two tables: one shared law unless the data show a difference capture cannot make. The
    # transcripts' geometry and RNA's normaliser read the same tables; calibration and the priors keep `gdna_pmf`.
    max_size = int(tables.rna.size) - 1
    rna_fl = FragmentLengthModel.from_pmf(tables.rna, max_size)
    gdna_fl = FragmentLengthModel.from_pmf(tables.gdna, max_size)

    estimator = _setup_geometry_and_estimator(
        index, rna_fl, em_config, calibration, region_arrays, tables
    )

    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        rna_fl,
        gdna_fl,
        stats,
        estimator,
        scoring_cfg,
        log_every,
        annotations,
        # each multimapper's place weighs its candidates by its capture
        place_efficiencies(calibration, region_arrays, index),
    )

    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)

    if em_data.n_units == 0 or not multi_loci:
        return estimator

    priors = assemble_priors(calibration, region_arrays, multi_loci)
    partitions = partition_and_free(em_data, multi_loci)
    _run_locus_em_partitioned(
        estimator,
        partitions,
        multi_loci,
        index,
        priors.gdna_count,
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

    Opens the BAM **once** via the C++ BAM scanner, scans all fragments,
    trains models, buffers resolved fragments, then quantifies via per-locus EM.

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

    # -- Single BAM pass (C++ native scanner) --
    scan = config.scan
    stats, strand_models, buffer, calibration_payload = scan_and_buffer(bam_path, index, scan)
    _log_sj_tags(stats, scan)
    _log_record_tags(stats)
    _log_hit_verdicts(stats)

    # -- The second pass: drain the side buffer, BEFORE calibration --
    calibration_payload = _drain_side_buffer(
        calibration_payload, index, strand_models, seed=config.second_pass_seed
    )

    # -- Calibration: deconvolve every region and boundary into gDNA / RNA+ / RNA− --
    from .calibration import calibrate
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays

    _warn_if_calibration_strand_unidentifiable(strand_models)
    confidence = CONSTANTS.qc.strand_ci_confidence
    strand_ci_eps = strand_models.strand_specificity_ci_epsilon(confidence=confidence)
    logger.info(
        "[CAL] Strand trainer: n_spliced_obs=%d  ss_est=%.6f  ε_CI(%g%%)=%.4g",
        strand_models.n_observations,
        strand_models.strand_specificity,
        100.0 * confidence,
        strand_ci_eps,
    )

    # No alignment check here: CalibrationSubstrate.from_payload runs the identical one inside
    # calibrate(), microseconds later, and two copies of an invariant is one too many.
    region_arrays = RegionArrays.from_index(index)
    boundary_flags = build_boundary_flags_array(index)
    # The sj axis, in the accumulator's own sj slot order: where each sj attaches,
    # its TRANSCRIPT strand, and its exonic reach either side. The calibrator places it as a FACTOR on
    # its two endpoint regions — never as a message channel, since every sj closes an undirected
    # loop and the graph is not a polytree.
    sj = build_sj_geometry_arrays(index)

    # The two COMPONENT fragment-length models the calibrator's effective lengths need, refit on the
    # drained tally.
    fl_models = library_fl_models(calibration_payload, index)
    gdna_fl_pmf = fl_models.gdna_pmf
    _bins = np.arange(gdna_fl_pmf.size, dtype=np.float64)
    logger.info(
        "[CAL] FL models: n_gdna_pool=%.0f gdna_fl_mean=%.1f rna_fl_mean=%.1f global_fl_mean=%.1f",
        fl_models.n_gdna,
        float(np.dot(_bins, fl_models.gdna_pmf)),
        float(np.dot(_bins, fl_models.rna_pmf)),
        float(np.dot(_bins, fl_models.global_pmf)),
    )

    # NOTE: the buffer is NOT freed here — quant_from_buffer routes it below.
    calibration = calibrate(
        payload=calibration_payload,
        region_arrays=region_arrays,
        strand_model=strand_models,
        gdna_fl_pmf=gdna_fl_pmf,
        rna_fl_pmf=fl_models.rna_pmf,
        config=config.calibration,
        sj=sj,
        boundary_flags=boundary_flags,
    )

    # FRAGMENTS, not object incidences. The sum still runs over all three axes — gDNA is contained
    # in a region or crosses a boundary, RNA also jumps, and at a donor boundary the sj flux IS the
    # RNA that spliced there — but each axis is converted by its own population's
    # mass-per-crossing first. Adding the raw banks would count one fragment once per boundary it
    # crossed AND once per sj it used, which reads `f_gdna` badly low on a contaminated library.
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

    # -- The EM's two length tables: one shared law unless the data show a difference capture cannot make --
    tables = length_tables(fl_models, calibration, region_arrays)
    decision = tables.decision
    logger.info(
        "[FL] length tables: kappa=%.4f captured=%.4f contrast=%.4f bound=%.4f noise=%.3g weight=%.3f "
        "rna_mean=%.1f gdna_table_mean=%.1f",
        decision.kappa,
        decision.captured,
        decision.contrast,
        decision.bound,
        decision.noise,
        decision.weight,
        float(np.dot(_bins, tables.rna)),
        float(np.dot(_bins, tables.gdna)),
    )

    # -- Quantification: calibration prior → per-locus EM --
    # Annotation table for the optional annotated-BAM pass (opt-in).
    annotations = None
    if config.annotated_bam_path is not None:
        from .annotate import AnnotationTable

        annotations = AnnotationTable.create(
            capacity=max(
                buffer.total_fragments + CONSTANTS.resources.annotation_table_padding,
                CONSTANTS.resources.annotation_table_min_capacity,
            )
        )

    try:
        estimator = quant_from_buffer(
            buffer,
            index,
            strand_models,
            tables,
            region_arrays,
            stats,
            calibration,
            em_config=config.em,
            scoring=config.scoring,
            log_every=scan.log_every,
            annotations=annotations,
            emit_locus_stats=config.emit_locus_stats,
        )
    finally:
        buffer.cleanup()

    # Genome-wide gDNA track for the QC report (+ bedGraph). Pure persistence of
    # the per-region calibration solution.
    from .calibration.track import build_gdna_track

    calibration_track = build_gdna_track(calibration, region_arrays, index.ref_names)

    result = PipelineResult(
        stats=stats,
        strand_models=strand_models,
        estimator=estimator,
        pipeline_config=config,
        calibration=calibration,
        fl_models=fl_models,
        length_tables=tables,
        calibration_track=calibration_track,
    )
    if before_annotated_bam is not None:
        before_annotated_bam(result)

    # -- The annotated-BAM pass: re-read the BAM and write it with assignment tags (opt-in) --
    if annotations is not None:
        from .annotate import write_annotated_bam

        write_annotated_bam(
            bam_path,
            str(config.annotated_bam_path),
            annotations,
            index,
            stats=stats,
            skip_duplicates=scan.skip_duplicates,
            include_multimap=scan.include_multimap,
            sj_spec=_sj_tag_to_spec(scan.sj_strand_tag),
            locus_id_per_transcript=estimator.locus_id_per_transcript,
        )
    return result
