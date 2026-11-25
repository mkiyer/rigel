"""What ``rigel quant`` writes: its tables, its calibration track and ``summary.json``, by the contract of
:mod:`rigel.outputs`. It takes the run's results and the command as the run resolved it; the CLI owns the arguments
(the command record, ``config.yaml``).
"""

from __future__ import annotations

import json
import logging
import time
from collections.abc import Callable
from pathlib import Path

from . import outputs

#: summary.json's floats (fractions, strand rates, quantiles) are rounded to this many decimals.
SUMMARY_DECIMALS = 6

#: summary.json's fragment counts are rounded to this many decimals: a hundredth of a fragment.
SUMMARY_COUNT_DECIMALS = 2

#: summary.json's lengths in bp (a length law's mean, median and standard deviation) are rounded to this many
#: decimals: a hundredth of a base.
SUMMARY_LENGTH_DECIMALS = 2

#: summary.json's strand posterior variance keeps this many significant digits: it is p(1 − p)/n, of order 1e-8 at a
#: million fragments, where fixed decimals print 0, and at most 0.25, so 8 significant digits are at least 8 decimals.
SUMMARY_VARIANCE_DIGITS = 8

#: ``fragment_lengths.feather``'s pure-pool categories, in ``FragmentPool`` order. Named here, not derived from the
#: enum, because they are an output contract: a rename in the enum must not rename a reported category.
_POOL_REPORT_NAMES = (
    "pool_gdna_intergenic",
    "pool_gdna_intronic",
    "pool_gdna_intron_exon",
    "pool_gdna_intergenic_exon",
    "pool_rna_spliced",
)

#: What a store-built index records of the store its unobservable mask was read from (`rigel.observable`).
_MASK_STORE_KEYS = (
    "aligner",
    "aligner_version",
    "read_length_bins",
    "frag_len_mean",
    "frag_len_sd",
    "error_rate",
    "references_from_fasta",
)


def write_quant_outputs(
    result,
    index,
    output_dir: Path,
    *,
    tsv: bool,
    command: dict,
    write_config: Callable[[Path], None],
) -> None:
    """Write everything ``rigel quant`` produces (:mod:`rigel.outputs`): the four quantification tables,
    ``locus_stats.feather`` when the run collected it, ``fragment_lengths.feather``, the calibration track,
    ``config.yaml`` (``write_config(path)``), and last ``summary.json``, the completion marker, which records
    ``command`` as the run was invoked. The previous run's marker and products go first, so an interrupted run leaves
    no directory a reader takes for a complete one. ``tsv`` adds each table's TSV mirror."""
    estimator = result.estimator
    outputs.begin_run(output_dir)
    for key, val in sorted(result.stats.to_dict().items()):
        if isinstance(val, (int, float)):
            logging.info(f"  {key}: {val:,}")

    tables = _write_tables(estimator, index, output_dir, tsv)
    fl_summary, fl_histogram_df = fragment_length_report(result.fl_models)

    # The raw fragment-length histograms, long form: one (category, length, count) row per 1-bp bin
    fragment_lengths_path = outputs.table_path(output_dir, outputs.FRAGMENT_LENGTHS)
    _write_table(fl_histogram_df, fragment_lengths_path, tsv)
    logging.info(f"[DONE] Wrote {fragment_lengths_path} ({len(fl_histogram_df)} rows)")

    _write_calibration_track(result.calibration_track, output_dir)
    outputs.publish(output_dir / outputs.CONFIG, write_config)

    summary = _summary(result, index, command, tables, fl_summary)
    summary_path = output_dir / outputs.SUMMARY

    def _dump(path: Path) -> None:
        with open(path, "w") as f:
            json.dump(summary, f, indent=2)

    outputs.publish(summary_path, _dump)
    logging.info(f"[DONE] Summary written to {summary_path}")


def _write_table(df, path: Path, tsv: bool) -> None:
    """One output table: Feather v2 with zstd, and its TSV mirror under ``--tsv``, each published whole."""
    outputs.publish(path, lambda p: df.to_feather(str(p), compression=outputs.TABLE_COMPRESSION))
    if tsv:
        outputs.publish(
            path.with_suffix(".tsv"), lambda p: df.to_csv(str(p), sep="\t", index=False)
        )


def _write_tables(estimator, index, output_dir: Path, tsv: bool) -> dict:
    """Write the four quantification tables, and the per-locus EM profile when the run collected one; return the
    tables by name."""
    import pandas as pd

    tables = {
        "quant": estimator.get_counts_df(index),
        "gene_quant": estimator.get_gene_counts_df(index),
        "nascent_quant": estimator.get_nascent_quant_df(index),
        "loci": estimator.get_loci_df(index),
    }
    paths = [outputs.table_path(output_dir, name) for name in outputs.QUANT_TABLES]
    for name, path in zip(outputs.QUANT_TABLES, paths):
        _write_table(tables[name], path, tsv)
    logging.info(f"[DONE] Wrote {', '.join(str(path) for path in paths)}")

    if estimator.locus_stats:
        locus_stats_path = outputs.table_path(output_dir, outputs.LOCUS_STATS)
        locus_stats_df = pd.DataFrame(estimator.locus_stats)
        outputs.publish(
            locus_stats_path,
            lambda p: locus_stats_df.to_feather(str(p), compression=outputs.TABLE_COMPRESSION),
        )
        logging.info(f"[DONE] Wrote {locus_stats_path} ({len(locus_stats_df)} loci)")
    return tables


def _write_calibration_track(track, output_dir: Path) -> None:
    """The genome-wide gDNA track, calibration's per-region gDNA level: a Feather (the report's input) and a bedGraph
    (for a genome browser)."""
    from .calibration.track import write_bedgraph

    track_path = outputs.table_path(output_dir, outputs.CALIBRATION_TRACK)
    bedgraph_path = output_dir / outputs.BEDGRAPH
    outputs.publish(
        track_path, lambda p: track.to_feather(str(p), compression=outputs.TABLE_COMPRESSION)
    )
    outputs.publish(bedgraph_path, lambda p: write_bedgraph(track, p))
    logging.info(f"[DONE] Wrote {track_path} + {bedgraph_path.name} ({len(track)} regions)")


def _summary(result, index, command: dict, tables: dict, fl_summary: dict) -> dict:
    """``summary.json``: the run's command and configuration, then what the scan saw, the strand and calibration
    answers, the length laws and the quantification's pools."""
    from . import __version__

    stats = result.stats
    return {
        "schema_version": outputs.SUMMARY_SCHEMA_VERSION,
        "rigel_version": __version__,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "command": command,
        "configuration": result.pipeline_config.to_dict(),
        # the observable genome every opportunity of the run was read over: the index's mask and where it came from
        "observable_genome": observable_genome_summary(index),
        "alignment_stats": _alignment_summary(stats),
        "read_pairs": _read_pair_summary(stats),
        "fragment_stats": _fragment_summary(stats),
        "strand_model": _strand_model_summary(result.strand_models),
        "calibration": _calibration_summary(result.calibration, result.junction_cap),
        "gdna_em_effective_length": _series_summary(tables["loci"], "gdna_em_effective_length"),
        "fragment_length": fl_summary,
        "length_tables": _length_tables_summary(result.length_tables),
        "quantification": _quantification_summary(
            stats, result.estimator, index, float(tables["quant"]["count"].sum())
        ),
    }


def observable_genome_summary(index) -> dict:
    """The observable genome the run read (`rigel.observable`): the index's record of what its unobservable mask is —
    ``rule``, and the store's ``read_length`` bin it reads (``None`` without a store) — the observable and unobservable
    bases over the genome and ``by_reference``, and the ``store`` the mask was read from (``None`` without one, when the
    mask is the reference's unknown bases)."""
    from .index import load_manifest

    record = (load_manifest(index.index_dir) or {}).get("observable") or {}
    df = index.unobservable_df
    masked = {ref: 0 for ref in index.ref_lengths}
    if df is not None and len(df):
        for ref, n in (df["end"] - df["start"]).groupby(df["ref"].astype(str)).sum().items():
            masked[str(ref)] += int(n)
    by_reference = {
        ref: {"observable_bases": int(n) - masked[ref], "unobservable_bases": masked[ref]}
        for ref, n in index.ref_lengths.items()
    }
    return {
        "rule": record.get("rule"),
        "read_length": record.get("read_length"),
        "observable_bases": sum(v["observable_bases"] for v in by_reference.values()),
        "unobservable_bases": sum(v["unobservable_bases"] for v in by_reference.values()),
        "store": {k: record[k] for k in _MASK_STORE_KEYS} if "aligner" in record else None,
        "by_reference": by_reference,
    }


def _alignment_summary(stats) -> dict:
    """The records: each flag counter counts every record carrying its flag (``proper_pairs`` the primary records
    flagged as a proper pair), and every record has one status."""
    return {
        "total_reads": stats.total,
        "mapped_reads": stats.total - stats.unmapped,
        "unmapped_reads": stats.unmapped,
        "secondary_reads": stats.secondary,
        "supplementary_reads": stats.supplementary,
        "proper_pairs": stats.proper_pair,
        "duplicate_reads": stats.duplicate,
        "qc_fail_reads": stats.qc_fail,
        # One status per record, summing to ``total_reads``: filtered (unmapped, QC-failed or a skipped duplicate),
        # used (a usable alignment's), unused (in a placed read pair, in no usable alignment), or unplaced (in a read
        # pair withheld or excluded by the multimap option)
        "records": {
            "filtered": stats.filtered,
            "used": stats.n_records_used,
            "unused": stats.n_records_unused,
            "unplaced": stats.n_records_unplaced,
        },
        # Each spliced record's junction-strand tag, read on that record: given a strand, given none, or ambiguous
        # because its XS and ts disagree; ``invalid`` counts the tag values that were neither '+' nor '-', each read
        # as absent.
        "sj_strand_tag": {
            "tagged": stats.n_sj_tag_tagged,
            "untagged": stats.n_sj_tag_untagged,
            "conflicting": stats.n_sj_tag_conflicting,
            "invalid": stats.n_sj_tag_invalid,
        },
        # NM and NH on each record: ``missing`` without the tag, ``invalid`` holding it as anything but an integer,
        # each read as absent (NM 0, NH 1); ``conflicting`` counts the alignments whose records disagree on NM, each
        # read with the larger
        "nm_tag": {
            "missing": stats.n_nm_tag_missing,
            "invalid": stats.n_nm_tag_invalid,
            "conflicting": stats.n_nm_conflicting,
        },
        "nh_tag": {"missing": stats.n_nh_tag_missing, "invalid": stats.n_nh_tag_invalid},
    }


def _read_pair_summary(stats) -> dict:
    """One outcome per read pair (a run of records sharing a name), summing to ``total``: excluded (every record
    filtered), withheld with no usable alignment by reason (an unpaired one also by the furthest check any of its
    alignments reached), excluded by the multimap option, or placed, unique or multimapping."""
    return {
        "total": stats.n_read_pairs,
        "excluded": stats.n_excluded,
        "withheld": {
            "unpaired": stats.n_withheld_unpaired,
            "unannotated_splice": stats.n_withheld_unannotated_splice,
            "off_index": stats.n_withheld_off_index,
            "chimeric": stats.n_withheld_chimeric,
        },
        "unpaired_by_reason": {
            "missing_mate": stats.n_withheld_unpaired_missing_mate,
            "reference": stats.n_withheld_unpaired_reference,
            "strand": stats.n_withheld_unpaired_strand,
            "outward": stats.n_withheld_unpaired_outward,
            "length": stats.n_withheld_unpaired_length,
        },
        "excluded_by_option": stats.n_excluded_by_option,
        "placed": {
            "unique": stats.n_placed_unique,
            "multimapping": stats.n_placed_multimapping,
        },
    }


def _fragment_summary(stats) -> dict:
    """The fragments: genic and intergenic, the hits rejected by reason, the splice census and the mismatch rate."""
    rate = stats.mismatch_rate
    return {
        "total": stats.n_fragments,
        "genic": stats.n_fragments - stats.n_intergenic,
        "intergenic": stats.n_intergenic,
        "with_annotated_sj": stats.n_with_annotated_sj,
        # Templates whose NH exceeds the alignments reported for them (a truncated output)
        "nh_above_reported": stats.n_nh_above_reported,
        "hits_rejected": {
            "missing_mate": stats.n_hits_rejected_missing_mate,
            "reference": stats.n_hits_rejected_reference,
            "strand": stats.n_hits_rejected_strand,
            "outward": stats.n_hits_rejected_outward,
            "length": stats.n_hits_rejected_length,
            "contradictory": stats.n_hits_rejected_contradictory,
            "novel_junction": stats.n_hits_rejected_novel_junction,
            "unsupported": stats.n_hits_rejected_unsupported,
            "off_index": stats.n_hits_rejected_off_index,
            "chimeric": stats.n_hits_rejected_chimeric,
        },
        "splice": _splice_summary(stats),
        # The library's mismatch rate, read on its unique unspliced hits (their NM over the read bases they align):
        # the price of a multimapper's hit's mismatches and soft-clipped bases; null with no such hit
        "mismatch_rate": {
            "rate": None if rate is None else round(rate, SUMMARY_DECIMALS),
            "mismatches": stats.n_unique_unspliced_mismatches,
            "bases": stats.n_unique_unspliced_bases,
        },
    }


def _splice_summary(stats) -> dict:
    """The splice-type breakdown: the scanner's per-fragment census, which counts every fragment it offers the
    accumulator (and closes against the accumulator's qc block). It reports the ``spliced_implicit`` class the
    ``with_*_sj`` counters omit."""
    from .splice import SpliceType, census_field

    census = {
        "unspliced": SpliceType.UNSPLICED,
        "spliced_annotated": SpliceType.SPLICED_ANNOT,
        "spliced_implicit": SpliceType.SPLICED_IMPLICIT,
    }
    return {label: int(getattr(stats, census_field(stype))) for label, stype in census.items()}


def _strand_model_summary(sm) -> dict:
    """The strand model: the protocol and its specificity from the spliced reads, the junction table behind it, the
    junction fit (:func:`strand_fit_summary`) and the all-exonic diagnostic."""
    spliced = sm.exonic_spliced
    ci_lo, ci_hi = spliced.wald_interval_95()
    return {
        "strand_specificity": round(sm.strand_specificity, SUMMARY_DECIMALS),
        "p_r1_sense": round(sm.p_r1_sense, SUMMARY_DECIMALS),
        "read1_sense": bool(sm.read1_sense),
        "n_training_fragments": sm.n_observations,
        "posterior_variance": _significant(spliced.posterior_variance(), SUMMARY_VARIANCE_DIGITS),
        "ci_95": [round(ci_lo, SUMMARY_DECIMALS), round(ci_hi, SUMMARY_DECIMALS)],
        # The per-sj strand table this 2×2 is the marginal of, and the population the junction fit reads κ from.
        "sj": sm.sj_table.to_dict(),
        **strand_fit_summary(sm),
        # The all-exonic (RNA + gDNA) diagnostic model, never used for scoring. A positive ``contamination_gap`` means
        # unstranded gDNA pulls the all-exonic estimate toward 0.5 relative to the spliced-only model: a sign of gDNA
        # contamination.
        "diagnostics": {
            "exonic_all_specificity": round(sm.exonic.strand_specificity, SUMMARY_DECIMALS),
            "exonic_all_p_r1_sense": round(sm.exonic.p_r1_sense, SUMMARY_DECIMALS),
            "exonic_all_n_fragments": sm.exonic.n_observations,
            "contamination_gap": round(
                spliced.strand_specificity - sm.exonic.strand_specificity, SUMMARY_DECIMALS
            ),
        },
    }


def _significant(x: float, digits: int) -> float:
    """``x`` rounded to ``digits`` significant digits."""
    return float(f"{x:.{digits}g}")


def strand_fit_summary(strand_models) -> dict:
    """``summary.json``'s view of the junction fit behind κ.

    κ, the spliced reads' sense fraction, pooled over every junction and over the junctions the three-class fit
    credits to genuine RNA; each class's share of the junctions and of their fragments (shallow splice artifacts
    are most of the junctions and few of the fragments); whether the strand channel is live, the protocol decision
    calibration's strand channel reads, made here by the same two functions on the same 2×2; and the strand
    overdispersion, 0 by policy for gDNA and RNA alike, so the channel is binomial. Where the channel is dead the
    classes are not identified and calibration reads the pooled counts (``StrandModel.genuine_n_same``), so the
    genuine κ is the pooled one and there are no class shares."""
    from .calibration.strand_balance import fit_strand_balance
    from .junction_fit import junction_classes

    table = strand_models.sj_table
    n_sense, depth = float(table.n_sense.sum()), float(table.depth.sum())
    pooled = round((n_sense + 1.0) / (depth + 2.0), SUMMARY_DECIMALS)
    names = ("genuine", "artifact", "reversed")

    def shares(values):
        return (
            None
            if values is None
            else {k: round(v, SUMMARY_DECIMALS) for k, v in zip(names, values)}
        )

    live = fit_strand_balance(strand_models).live
    if live:
        classes = junction_classes(table.n_sense, table.depth)
        genuine = round(classes.kappa, SUMMARY_DECIMALS)
        by_junction, by_fragment = shares(classes.by_junction), shares(classes.by_fragment)
    else:
        genuine, by_junction, by_fragment = pooled, None, None
    return {
        "kappa": {"pooled": pooled, "genuine": genuine},
        "junction_classes": {"by_junction": by_junction, "by_fragment": by_fragment},
        "strand_live": live,
        "strand_overdispersion": {"value": 0.0, "model": "binomial"},
    }


def _calibration_summary(cal, junction_cap: float) -> dict:
    """Calibration's library scalars and the junction cap fitted from them, never a per-region dict."""
    return {
        # full precision: a density is a few fragments per Mb, which a fixed count of decimals rounds to 0
        "gdna_density_global": float(cal.gdna_density_global),
        # calibration's capture answer, which the ruler reads: the gDNA bulk density the capture efficiencies are
        # read against (None = no gDNA landscape, nothing contracted), the weight the reading is raised to (the share
        # of the reads whose own gDNA would show; 1 keeps the reading, 0 claims no capture) and the landscape's
        # saturated level against that bulk at that weight (1 = no level above it)
        "gdna_bulk_density": (
            None if cal.gdna_bulk_density is None else float(cal.gdna_bulk_density)
        ),
        "gdna_capture_weight": round(float(cal.gdna_capture_weight), SUMMARY_DECIMALS),
        "gdna_capture_saturation": float(cal.gdna_capture_saturation),
        # the cap on a junction's capture price in the EM's lengths: the saturated level, raised as far as RNA's own
        # spliced fragments support (equal to it where they raise nothing)
        "junction_capture_cap": float(junction_cap),
        "fraction_rna_sense": round(float(cal.rna_sense_frac), SUMMARY_DECIMALS),
        "n_regions": int(cal.n_regions),
        "n_boundaries": int(cal.n_boundaries),
        "n_sj": int(cal.n_sj),
    }


def _series_summary(df, column: str) -> dict[str, float]:
    """A table column's min, median, 95th percentile and max; zeros when the table is empty."""
    if df.empty:
        return {"min": 0.0, "median": 0.0, "p95": 0.0, "max": 0.0}
    values = df[column]
    return {
        "min": round(float(values.min()), SUMMARY_DECIMALS),
        "median": round(float(values.median()), SUMMARY_DECIMALS),
        "p95": round(float(values.quantile(0.95)), SUMMARY_DECIMALS),
        "max": round(float(values.max()), SUMMARY_DECIMALS),
    }


def _length_tables_summary(tables) -> dict:
    """The EM's two length tables and the decision behind them (`calibration.shared_length_law`): ``weight`` 0 is one
    shared law for both origins, 1 gDNA's own law brought into RNA's captured frame, continuous between.
    ``rna_captured_fraction`` is RNA's captured fraction, ``captured`` the capture mixture's fraction, ``contrast``
    and ``bound`` the Hellinger distances the weight compares, ``noise`` the squared distance sampling noise alone
    adds; the three means are the two tables' and gDNA's uncaptured law's, in bp."""
    import numpy as np

    d = tables.decision
    L = np.arange(tables.rna.size, dtype=np.float64)
    return {
        "weight": round(float(d.weight), SUMMARY_DECIMALS),
        "rna_captured_fraction": round(float(d.rna_captured_fraction), SUMMARY_DECIMALS),
        "captured": round(float(d.captured), SUMMARY_DECIMALS),
        "contrast": round(float(d.contrast), SUMMARY_DECIMALS),
        "bound": round(float(d.bound), SUMMARY_DECIMALS),
        "noise": round(float(d.noise), SUMMARY_DECIMALS),
        "rna_table_mean": round(float(tables.rna @ L), SUMMARY_LENGTH_DECIMALS),
        "gdna_table_mean": round(float(tables.gdna @ L), SUMMARY_LENGTH_DECIMALS),
        "gdna_uncaptured_mean": round(float(tables.gdna_uncaptured @ L), SUMMARY_LENGTH_DECIMALS),
    }


def fragment_length_report(fl_models):
    """Split the fragment-length models into a lean summary + a tidy histogram table.

    Returns ``(summary_by_category, histogram_df)`` where

    * ``summary_by_category`` is ``{category: {n_observations, mean, std, median, mode, max_size, overflow_count,
      overflow_fraction}}`` (:func:`length_summary`) — small, human-readable, destined for ``summary.json`` (no raw
      per-bin arrays).
    * ``histogram_df`` is a long ``(category, length, count)`` DataFrame — the raw 1-bp histograms
      (:func:`length_histogram`), destined for ``fragment_lengths.feather``, which keeps them out of ``summary.json``.

    Every category measures one quantity, the accumulator's ``L``, over a stated population:

    ``global``
        the unconditional anchor: every deposited fragment, no purity condition.
    ``gdna`` / ``rna``
        the two modelled pools, the empirical views of what scoring and calibration actually use.
    the five ``pool_*`` rows
        the five pools unaggregated (none is pure). ``gdna`` is the fitted uniform-frame gDNA
        histogram and ``rna`` the de-tilted ``pool_rna_spliced``; the two ``*_exon`` crossing rows are
        the on-target gDNA, shown on their own so the off-target / on-target comparison is an output.

    There are no per-splice-type histograms: the scanner measures those by other rules than ``L``. The splice
    counts come from the scanner's census instead.
    """
    import pandas as pd

    models = {
        "global": fl_models.global_model(),
        "gdna": fl_models.gdna_model(),
        "rna": fl_models.rna_model(),
        **{
            name: fl_models._empirical(fl_models.pool_counts[i])
            for i, name in enumerate(_POOL_REPORT_NAMES)
        },
    }
    cats: list = []
    lengths: list = []
    counts: list = []
    for name, model in models.items():
        bins, values = length_histogram(model)
        cats += [name] * len(bins)
        lengths += bins.tolist()
        counts += values.tolist()
    histogram_df = pd.DataFrame(
        {
            "category": pd.Series(cats, dtype="category"),
            "length": pd.Series(lengths, dtype="int32"),
            "count": pd.Series(counts, dtype="float64"),
        }
    )
    return {name: length_summary(model) for name, model in models.items()}, histogram_df


def length_summary(model) -> dict:
    """A length law's statistics as ``summary.json`` reports them. ``mean``, ``std`` and ``mode`` describe the
    in-range lengths; the right-censored tail (lengths ``>= max_size``) is ``overflow_count`` and
    ``overflow_fraction``."""
    total = model.total_weight
    overflow = model.overflow_count
    return {
        "n_observations": int(model.n_observations),
        "mean": float(round(model.mean, SUMMARY_LENGTH_DECIMALS)),
        "std": float(round(model.std, SUMMARY_LENGTH_DECIMALS)),
        "median": float(round(model.median, SUMMARY_LENGTH_DECIMALS)),
        "mode": int(model.mode),
        "max_size": int(model.max_size),
        "overflow_count": float(round(overflow, SUMMARY_COUNT_DECIMALS)),
        "overflow_fraction": float(round(overflow / total, SUMMARY_DECIMALS)) if total > 0 else 0.0,
    }


def length_histogram(model):
    """``(lengths, counts)``: a length law's 1-bp histogram, trimmed to its first and last non-empty bins."""
    import numpy as np

    nonzero = np.flatnonzero(model.counts)
    if nonzero.size == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float64)
    lo, hi = int(nonzero[0]), int(nonzero[-1]) + 1
    return np.arange(lo, hi), model.counts[lo:hi]


def _quantification_summary(stats, estimator, index, count_annotated: float) -> dict:
    """The quantification: what was counted and assigned, then the pools (:func:`quantification_pools`)."""
    is_synthetic = index.t_df["is_synthetic"].to_numpy(dtype=bool)
    n_nascent = int(is_synthetic.sum())
    return {
        # the annotated transcripts, single-exon included; ``n_nascent`` counts the synthetic spans
        "n_transcripts": int(is_synthetic.size) - n_nascent,
        "n_nascent": n_nascent,
        "n_genes": index.num_annotated_genes,
        "n_loci": len(estimator.locus_results),
        "n_unambig_assigned": stats.deterministic_unambig_units,
        "n_em_assigned": (
            stats.em_routed_unambig_units
            + stats.em_routed_ambig_same_strand_units
            + stats.em_routed_ambig_opp_strand_units
            + stats.em_routed_multimapper_units
        ),
        # dropped by a hard scoring gate (an alpha of 0): every candidate failed it, so no pool holds them
        "n_gated_out": stats.n_gated_out,
        **quantification_pools(
            stats,
            count_annotated,
            float(estimator.synthetic_count),
            float(estimator.gdna_em_count),
            float(estimator.gdna_intergenic_em_count),
        ),
    }


def quantification_pools(
    stats,
    count_annotated: float,
    count_nascent: float,
    count_gdna_em: float,
    count_gdna_intergenic_em: float,
) -> dict:
    """``summary.json``'s pools and their fractions of every fragment counted.

    An intergenic fragment reaches no locus. It is unspliced (a spliced hit the scan keeps has a transcript holding
    its junctions), so gDNA by construction: ``count_gdna`` is the EM's genic gDNA (``count_gdna_em``) plus the
    intergenic fragments; leaving them out would halve the observed contamination of a library with substantial
    gDNA. The annotated-transcript (``fraction_annotated``), synthetic-span (``fraction_nascent``) and gDNA fractions
    sum to one; ``count_gdna_em`` and ``fraction_gdna_em`` are the EM's genic part alone, and ``count_intergenic`` and
    ``fraction_intergenic`` the intergenic part, so ``count_gdna = count_gdna_em + count_intergenic``."""
    # The templates that reach no locus (``PipelineStats.n_intergenic``), and the gDNA the EM put at a multimapper's
    # intergenic place: neither is in a locus.
    intergenic = stats.n_intergenic + count_gdna_intergenic_em
    total_all = count_annotated + count_nascent + count_gdna_em + intergenic
    count_gdna = count_gdna_em + intergenic

    def share(x: float) -> float:
        return round(x / total_all, SUMMARY_DECIMALS) if total_all > 0 else 0.0

    return {
        "count_annotated": round(count_annotated, SUMMARY_COUNT_DECIMALS),
        "count_nascent": round(count_nascent, SUMMARY_COUNT_DECIMALS),
        "count_gdna": round(count_gdna, SUMMARY_COUNT_DECIMALS),
        "count_gdna_em": round(count_gdna_em, SUMMARY_COUNT_DECIMALS),
        "fraction_annotated": share(count_annotated),
        "fraction_nascent": share(count_nascent),
        "fraction_gdna": share(count_gdna),
        "fraction_gdna_em": share(count_gdna_em),
        "count_intergenic": round(intergenic, SUMMARY_COUNT_DECIMALS),
        "fraction_intergenic": share(intergenic),
    }
