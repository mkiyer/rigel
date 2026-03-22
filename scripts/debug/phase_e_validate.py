#!/usr/bin/env python3
"""Phase E Validation: re-run rigel with Phase C code on existing benchmark BAMs.

Rebuilds the index (Phase C: no nrna.feather), runs rigel on the 3 hardest
benchmark conditions (oracle aligner), and compares per-transcript counts
and pool-level metrics against the pre-Phase-C baseline.

Usage:
    conda activate rigel
    python scripts/debug/phase_e_validate.py
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

# ── Paths ──────────────────────────────────────────────────────
BENCH_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")
SIM_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_salmon")
GENOME = Path("/Users/mkiyer/Downloads/rigel_runs/refs/human/genome_controls.fasta.bgz")
GTF = Path("/Users/mkiyer/Downloads/rigel_runs/refs/human/genes_controls.gtf.gz")
BASELINE_JSON = BENCH_DIR / "summary.json"
OUTPUT_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune/phase_e_results")

CONDITIONS = [
    "gdna_high_ss_1.00_nrna_default",
    "gdna_high_ss_0.90_nrna_default",
    "gdna_high_ss_0.50_nrna_default",
]


def parse_truth_from_fastq(
    r1_path: Path,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Parse ground-truth from FASTQ read names."""
    mrna_counts: Counter[str] = Counter()
    nrna_counts: Counter[str] = Counter()
    n_gdna = 0
    opener = gzip.open if str(r1_path).endswith(".gz") else open
    with opener(r1_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            qname = line[1:].strip()
            if qname.endswith("/1"):
                qname = qname[:-2]
            t_id = qname.split(":")[0]
            if t_id.startswith("gdna"):
                n_gdna += 1
            elif t_id.startswith("nrna_"):
                nrna_counts[t_id[5:]] += 1
            else:
                mrna_counts[t_id] += 1
    return dict(mrna_counts), dict(nrna_counts), n_gdna


def score_tool(
    truth: dict[str, float], observed: dict[str, float]
) -> dict[str, float]:
    """Compare observed vs truth → metrics dict."""
    all_tids = sorted(set(truth) | set(observed))
    n = len(all_tids)
    truth_arr = np.array([truth.get(t, 0.0) for t in all_tids])
    obs_arr = np.array([observed.get(t, 0.0) for t in all_tids])
    abs_err = np.abs(truth_arr - obs_arr)
    mae = float(abs_err.mean())
    rmse = float(np.sqrt(np.mean((truth_arr - obs_arr) ** 2)))
    if truth_arr.std() > 0 and obs_arr.std() > 0:
        pearson = float(np.corrcoef(truth_arr, obs_arr)[0, 1])
        from scipy.stats import spearmanr

        spearman = float(spearmanr(truth_arr, obs_arr).statistic)
    else:
        pearson = 0.0
        spearman = 0.0
    return {
        "n_transcripts": n,
        "total_truth": float(truth_arr.sum()),
        "total_observed": float(obs_arr.sum()),
        "total_abs_error": float(abs_err.sum()),
        "mean_abs_error": mae,
        "rmse": rmse,
        "pearson": pearson,
        "spearman": spearman,
    }


def load_baseline(conditions: list[str]) -> dict:
    """Load baseline metrics from summary.json."""
    with open(BASELINE_JSON) as f:
        data = json.load(f)
    baseline = {}
    for entry in data:
        ds = entry["dataset_name"]
        al = entry["aligner"]
        if ds in conditions and al == "oracle":
            baseline[ds] = entry
    return baseline


def build_index(index_dir: Path):
    """Build rigel index with Phase C code."""
    from rigel.index import TranscriptIndex

    print(f"Building index at {index_dir} ...", flush=True)
    required = ["transcripts.feather", "intervals.feather", "ref_lengths.feather"]
    if all((index_dir / f).exists() for f in required):
        print("  Index already exists, loading...", flush=True)
    else:
        t0 = time.monotonic()
        TranscriptIndex.build(
            fasta_file=str(GENOME),
            gtf_file=str(GTF),
            output_dir=str(index_dir),
            write_tsv=True,
            gtf_parse_mode="warn-skip",
        )
        elapsed = time.monotonic() - t0
        print(f"  Index built in {elapsed:.1f}s", flush=True)
    return TranscriptIndex.load(str(index_dir))


def run_rigel_on_bam(bam_path: Path, index):
    """Run rigel pipeline on a BAM, return (counts_dict, pool_counts, elapsed)."""
    from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
    from rigel.pipeline import run_pipeline

    cfg = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto", include_multimap=True),
        scoring=FragmentScoringConfig(),
    )
    print(f"  Running rigel on {bam_path.name} ...", flush=True)
    t0 = time.monotonic()
    pipe = run_pipeline(bam_path, index, config=cfg)
    elapsed = time.monotonic() - t0

    counts_df = pipe.estimator.get_counts_df(index)
    counts = {
        row.transcript_id: float(row.mrna)
        for row in counts_df.itertuples(index=False)
    }
    mature_pred = float(sum(counts.values()))
    nascent_pred = float(pipe.estimator.nrna_em_count)
    genomic_pred = float(pipe.stats.n_gdna_total)
    pool_counts = {
        "mature_rna": mature_pred,
        "nascent_rna": nascent_pred,
        "genomic_dna": genomic_pred,
    }
    n_frags = pipe.stats.total
    print(f"    Done in {elapsed:.1f}s ({n_frags} fragments)", flush=True)
    return counts, pool_counts, elapsed


def compare_metrics(label: str, baseline_m: dict, new_m: dict):
    """Print side-by-side comparison."""
    keys = ["mean_abs_error", "rmse", "pearson", "spearman", "total_truth", "total_observed", "total_abs_error"]
    print(f"\n  {label}:")
    print(f"    {'metric':<20s} {'baseline':>14s} {'new':>14s} {'delta':>14s} {'pct_change':>12s}")
    print("    " + "-" * 76)
    for k in keys:
        bv = baseline_m[k]
        nv = new_m[k]
        delta = nv - bv
        pct = (delta / bv * 100) if bv != 0 else 0.0
        flag = " !!!" if abs(pct) > 5.0 else ""
        print(f"    {k:<20s} {bv:>14.4f} {nv:>14.4f} {delta:>+14.4f} {pct:>+11.2f}%{flag}")


def compare_pools(baseline_pools: dict, new_pools: dict):
    """Print pool count comparison."""
    print("\n  Pool counts:")
    print(f"    {'pool':<15s} {'baseline':>14s} {'new':>14s} {'delta':>14s} {'pct_change':>12s}")
    print("    " + "-" * 67)
    for k in ["mature_rna", "nascent_rna", "genomic_dna"]:
        bv = baseline_pools[k]
        nv = new_pools[k]
        delta = nv - bv
        pct = (delta / bv * 100) if bv != 0 else 0.0
        flag = " !!!" if abs(pct) > 2.0 else ""
        print(f"    {k:<15s} {bv:>14.1f} {nv:>14.1f} {delta:>+14.1f} {pct:>+11.2f}%{flag}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load baseline
    print("Loading baseline results...", flush=True)
    baseline = load_baseline(CONDITIONS)
    for cond in CONDITIONS:
        if cond not in baseline:
            print(f"ERROR: {cond} not found in baseline summary.json", file=sys.stderr)
            sys.exit(1)
    print(f"  Loaded {len(baseline)} baseline conditions\n")

    # Build index
    index_dir = OUTPUT_DIR / "rigel_index"
    index = build_index(index_dir)
    print(f"  Index has {index.num_transcripts} transcripts")
    n_syn = int(index.t_df["is_synthetic_nrna"].sum())
    print(f"  Synthetic nRNA transcripts: {n_syn}\n")

    # Build basic-filtered transcript → gene map (matching baseline methodology)
    from rigel.transcript import Transcript as _Tx

    basic_transcripts = _Tx.read_gtf(str(GTF), parse_mode="warn-skip")
    basic_transcripts = [t for t in basic_transcripts if t.is_basic]
    basic_tid_to_gene = {t.t_id: t.g_id for t in basic_transcripts}
    print(f"  Basic transcript → gene map: {len(basic_tid_to_gene)} entries\n")

    # Run rigel on each condition and compare
    all_results = {}
    for cond in CONDITIONS:
        print(f"\n{'='*72}")
        print(f"Condition: {cond}")
        print(f"{'='*72}")

        bam_path = BENCH_DIR / cond / "align_oracle" / "reads_namesort.bam"
        if not bam_path.exists():
            print(f"  BAM not found: {bam_path}", file=sys.stderr)
            continue

        # Parse truth
        r1_path = SIM_DIR / cond / "sim_R1.fq.gz"
        print(f"  Parsing truth from {r1_path.name} ...", flush=True)
        mrna_truth, nrna_truth, n_gdna_truth = parse_truth_from_fastq(r1_path)
        truth = {tid: float(c) for tid, c in mrna_truth.items()}
        print(f"    mRNA truth: {sum(mrna_truth.values())} frags, "
              f"nRNA truth: {sum(nrna_truth.values())} frags, "
              f"gDNA truth: {n_gdna_truth} frags")

        # Run rigel
        counts, pool_counts, elapsed = run_rigel_on_bam(bam_path, index)

        # Filter out synthetic nRNA transcripts from counts — they have mrna=0
        # and should not participate in mRNA comparison
        syn_tids = set(
            index.t_df.loc[index.t_df["is_synthetic_nrna"], "t_id"].values
        )
        counts = {tid: c for tid, c in counts.items() if tid not in syn_tids}

        # Score: transcript-level
        new_transcript_metrics = score_tool(truth, counts)
        new_transcript_metrics["elapsed_sec"] = elapsed

        # Score: gene-level — use basic-filtered transcripts for gene mapping
        # to match the baseline benchmark methodology exactly
        truth_gene: dict[str, float] = {}
        for tid, c in truth.items():
            gid = basic_tid_to_gene.get(tid, tid)
            truth_gene[gid] = truth_gene.get(gid, 0.0) + c
        obs_gene: dict[str, float] = {}
        for tid, c in counts.items():
            gid = basic_tid_to_gene.get(tid, tid)
            obs_gene[gid] = obs_gene.get(gid, 0.0) + c
        new_gene_metrics = score_tool(truth_gene, obs_gene)

        # Compare against baseline
        bl = baseline[cond]
        bl_tm = bl["transcript_metrics"]["rigel_default_oracle"]
        bl_gm = bl["gene_metrics"]["rigel_default_oracle"]
        bl_pools = bl["pool_counts"]["rigel_default_oracle"]

        compare_metrics("TRANSCRIPT-LEVEL", bl_tm, new_transcript_metrics)
        compare_metrics("GENE-LEVEL", bl_gm, new_gene_metrics)
        compare_pools(bl_pools, pool_counts)

        all_results[cond] = {
            "transcript_metrics": new_transcript_metrics,
            "gene_metrics": new_gene_metrics,
            "pool_counts": pool_counts,
        }

    # Save results
    results_file = OUTPUT_DIR / "phase_e_results.json"
    with open(results_file, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n\nResults saved to {results_file}")

    # Summary verdict
    print("\n" + "=" * 72)
    print("VALIDATION SUMMARY")
    print("=" * 72)
    any_regression = False
    for cond in CONDITIONS:
        if cond not in all_results:
            continue
        bl = baseline[cond]
        bl_tm = bl["transcript_metrics"]["rigel_default_oracle"]
        new_tm = all_results[cond]["transcript_metrics"]
        mae_pct = abs(new_tm["mean_abs_error"] - bl_tm["mean_abs_error"]) / bl_tm["mean_abs_error"] * 100
        pearson_delta = abs(new_tm["pearson"] - bl_tm["pearson"])
        rmse_pct = abs(new_tm["rmse"] - bl_tm["rmse"]) / bl_tm["rmse"] * 100
        status = "PASS" if mae_pct < 5 and pearson_delta < 0.001 and rmse_pct < 5 else "FAIL"
        if status == "FAIL":
            any_regression = True
        print(f"  {cond}: {status} (MAE delta={mae_pct:.2f}%, pearson delta={pearson_delta:.6f}, RMSE delta={rmse_pct:.2f}%)")

    if any_regression:
        print("\n*** REGRESSIONS DETECTED — investigate above ***")
        sys.exit(1)
    else:
        print("\n*** ALL CONDITIONS PASS — Phase C has minimal impact ***")
        sys.exit(0)


if __name__ == "__main__":
    main()
