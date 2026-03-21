#!/usr/bin/env python3
"""Re-run Rigel quantification on existing benchmark BAMs.

After the multimapper counting bugfix in bam_scanner.cpp, this script
re-runs only the Rigel pipeline (default + no_prune configs) on
existing name-sorted BAMs from benchmark_output_v6_prune, then
merges updated Rigel results with existing salmon/kallisto data.

Usage:
    conda activate rigel
    python scripts/debug/rerun_rigel_v6.py
"""

import copy
import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark import (
    RigelConfig,
    parse_truth_from_fastq,
    run_rigel_tool,
    score_tool,
    aggregate_to_genes,
    load_transcripts,
    write_per_transcript_csv,
    load_sim_manifest,
    _rigel_tool_name,
    _get_peak_rss_mb,
    AlignerConfig,
)
from rigel.index import TranscriptIndex

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Configuration ────────────────────────────────────────────────
BENCH_DIR = Path("~/Downloads/rigel_runs/benchmark_output_v6_prune").expanduser()
SIM_DIR = Path("~/Downloads/rigel_runs/sim_ccle_hela_salmon").expanduser()
GENOME = Path("~/Downloads/rigel_runs/refs/human/genome_controls.fasta.bgz").expanduser()
GTF = Path("~/Downloads/rigel_runs/refs/human/genes_controls.gtf.gz").expanduser()

RIGEL_CONFIGS = [
    RigelConfig(name="default", params={"prune_threshold": 0.1}),
    RigelConfig(name="no_prune", params={"prune_threshold": -1.0}),
]

ALIGNERS = [
    AlignerConfig(name="oracle", type="oracle"),
    AlignerConfig(name="minimap2", type="minimap2"),
]

MULTI_RIGEL = len(RIGEL_CONFIGS) > 1


def main():
    t0_total = time.monotonic()

    # Load original summary
    orig_json = BENCH_DIR / "summary.json"
    with open(orig_json) as f:
        orig_results = json.load(f)
    log.info("Loaded %d entries from original summary.json", len(orig_results))

    # Build index for original → fast lookup
    orig_lookup = {}
    for entry in orig_results:
        key = (entry["dataset_name"], entry["aligner"])
        orig_lookup[key] = entry

    # Load transcripts for scoring
    transcripts = load_transcripts(GTF, transcript_filter="basic")

    # Load rigel index
    rigel_index = TranscriptIndex.load(BENCH_DIR / "rigel_index")

    # Discover datasets from sim manifest
    manifest, sim_datasets = load_sim_manifest(str(SIM_DIR))

    # Process each dataset × aligner
    updated_results = []
    for ds in sim_datasets:
        for ac in ALIGNERS:
            key = (ds.name, ac.name)
            orig_entry = orig_lookup.get(key)
            if orig_entry is None:
                log.warning("No original entry for %s, skipping", key)
                continue

            align_dir = BENCH_DIR / ds.name / f"align_{ac.name}"
            bam_ns = align_dir / "reads_namesort.bam"

            if not bam_ns.exists():
                log.warning("BAM not found: %s, skipping", bam_ns)
                continue

            # Parse truth
            mrna_truth, nrna_truth, n_gdna_truth = {}, {}, 0
            if ds.fastq_r1 and ds.fastq_r1.exists():
                mrna_truth, nrna_truth, n_gdna_truth = parse_truth_from_fastq(
                    ds.fastq_r1
                )

            mrna_truth_f = {k: float(v) for k, v in mrna_truth.items()}

            # Run rigel configs
            new_tx_metrics = {}
            new_gene_metrics = {}
            new_pool_counts = {}
            new_elapsed = {}
            new_peak_rss = {}
            new_throughput = {}
            new_tool_tx_counts = {}

            for hc in RIGEL_CONFIGS:
                hn = _rigel_tool_name(ac, hc, MULTI_RIGEL)
                log.info("  %s / %s / %s ...", ds.name, ac.name, hn)
                annotated_bam = align_dir / f"annotated_{hc.name}.bam"
                counts, elapsed, pools, n_frags = run_rigel_tool(
                    bam_ns, rigel_index, ds.pipeline_seed, hc,
                    annotated_bam_path=annotated_bam,
                )
                new_tool_tx_counts[hn] = counts
                new_pool_counts[hn] = pools
                new_elapsed[hn] = elapsed
                new_peak_rss[hn] = _get_peak_rss_mb()
                new_throughput[hn] = n_frags / elapsed if elapsed > 0 else 0.0

                if mrna_truth_f:
                    from dataclasses import asdict
                    tm = score_tool(hn, mrna_truth_f, counts, elapsed,
                                    new_peak_rss[hn], new_throughput[hn])
                    new_tx_metrics[hn] = asdict(tm)
                    gene_truth = aggregate_to_genes(mrna_truth_f, transcripts)
                    gene_obs = aggregate_to_genes(counts, transcripts)
                    gm = score_tool(hn, gene_truth, gene_obs, elapsed,
                                    new_peak_rss[hn], new_throughput[hn])
                    new_gene_metrics[hn] = asdict(gm)

                log.info("    done (%.1fs, %.0fMB RSS, %,.0f frags/s)",
                         elapsed, new_peak_rss[hn], new_throughput[hn])

            # Update per-transcript CSV
            if mrna_truth and new_tool_tx_counts:
                write_per_transcript_csv(
                    transcripts, mrna_truth, nrna_truth, new_tool_tx_counts,
                    BENCH_DIR / ds.name / f"per_transcript_counts_{ac.name}.csv",
                )

            # Merge into original entry
            updated = copy.deepcopy(orig_entry)
            updated["transcript_metrics"].update(new_tx_metrics)
            updated["gene_metrics"].update(new_gene_metrics)
            updated["pool_counts"].update(new_pool_counts)
            updated["elapsed"].update(new_elapsed)
            updated["peak_rss"].update(new_peak_rss)
            updated["throughput"].update(new_throughput)
            updated_results.append(updated)

    # Write updated summary
    out_json = BENCH_DIR / "summary.json"
    out_json_bak = BENCH_DIR / "summary.json.pre_fix"
    if not out_json_bak.exists():
        import shutil
        shutil.copy2(out_json, out_json_bak)
        log.info("Backed up original summary.json to %s", out_json_bak)

    with open(out_json, "w") as f:
        json.dump(updated_results, f, indent=2)
    log.info("Wrote updated summary.json (%d entries)", len(updated_results))

    total = time.monotonic() - t0_total
    log.info("Re-run complete in %.1fs", total)


if __name__ == "__main__":
    main()
