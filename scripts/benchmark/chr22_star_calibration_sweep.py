"""STAR-aligned calibration sweep on human chr22.

Counterpart to ``chr22_calibration_sweep.py`` but, instead of writing an
oracle BAM, generates paired-end FASTQ reads and aligns them with STAR
against a chr22-only STAR index.  This exercises the full realistic
alignment path (soft-clipping, secondary alignments, splice mapping,
unmapped reads, mate-inconsistency, etc.) that the oracle sweep bypasses.

Axes, fragment counts, and truth semantics mirror
``chr22_calibration_sweep.py`` so results are directly comparable.

Usage::

    python scripts/benchmark/chr22_star_calibration_sweep.py \\
        --out scripts/benchmark/results/calibration_star_chr22 \\
        --fasta   .../chr22_only/chr22.fasta.bgz \\
        --index   .../chr22_only/rigel_index \\
        --star-index .../chr22_only/star_index \\
        --star-bin /path/to/STAR
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimulator, SimConfig
from rigel.transcript import Transcript

# Reuse helpers from the oracle sweep so the two scripts stay in lock-step.
from chr22_calibration_sweep import (
    REF_NAME,
    _StringGenome,
    load_chr22_genome,
    select_transcripts,
    truth_lambda_gdna,
)

logger = logging.getLogger("chr22_star_calibration_sweep")


DEFAULT_FASTA = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/chr22_only/chr22.fasta.bgz"
)
DEFAULT_INDEX = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/chr22_only/rigel_index"
)
DEFAULT_STAR_INDEX = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/chr22_only/star_index"
)
DEFAULT_STAR_BIN = Path(
    "/home/mkiyer/sw/miniforge3/envs/alignable/bin/STAR"
)


# ---------------------------------------------------------------------------
# FASTQ with exact pool split
# ---------------------------------------------------------------------------

def write_fastq_fixed_split(
    sim: ReadSimulator,
    out_dir: Path,
    pool_split: tuple[int, int, int],
    *,
    prefix: str = "sim",
) -> tuple[Path, Path, dict]:
    """Emit paired FASTQ where (n_mrna, n_nrna, n_gdna) are exact.

    The stock ``ReadSimulator.simulate`` calls ``_compute_pool_split`` to
    derive the 3-way split from abundances.  Here we force the split by
    monkey-patching that method, then fall through to the existing
    generator pipeline so read-name conventions, quality scores, strand
    flipping, etc. are identical to the oracle BAM path.
    """
    n_mrna, n_nrna, n_gdna = pool_split
    n_total = n_mrna + n_nrna + n_gdna

    out_dir.mkdir(parents=True, exist_ok=True)
    r1_path = out_dir / f"{prefix}_R1.fq.gz"
    r2_path = out_dir / f"{prefix}_R2.fq.gz"

    # Force the pool split for this simulator instance only.
    sim._compute_pool_split = lambda _n: (n_mrna, n_nrna, n_gdna)  # type: ignore[method-assign]

    stats = {"n_written": 0}
    with gzip.open(r1_path, "wt") as r1_fh, gzip.open(r2_path, "wt") as r2_fh:
        for tup in sim.simulate(n_total):
            r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual = tup
            r1_fh.write(f"@{r1_name}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_fh.write(f"@{r2_name}\n{r2_seq}\n+\n{r2_qual}\n")
            stats["n_written"] += 1

    logger.info(
        "FASTQ: %d pairs written (target=%d, split=%s) → %s",
        stats["n_written"], n_total, pool_split, r1_path.name,
    )
    return r1_path, r2_path, stats


# ---------------------------------------------------------------------------
# STAR alignment
# ---------------------------------------------------------------------------

def align_star(
    r1: Path, r2: Path, out_dir: Path,
    *, star_bin: Path, star_index: Path, threads: int = 8,
) -> Path:
    """Run STAR against the chr22-only index, produce a name-sorted BAM.

    STAR options chosen to roughly match rigel's production runner:

      - ``--outSAMtype BAM Unsorted``            → raw BAM, we n-sort ourselves
      - ``--outFilterMultimapNmax 20``           → allow genuine multimappers
      - ``--outSAMattributes NH HI AS nM MD XS`` → NH required by rigel scanner
      - ``--readFilesCommand zcat``              → gzipped FASTQ
      - ``--outSAMstrandField intronMotif``      → XS tag for spliced reads
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "star_"
    cmd = [
        str(star_bin),
        "--runMode", "alignReads",
        "--genomeDir", str(star_index),
        "--readFilesIn", str(r1), str(r2),
        "--readFilesCommand", "zcat",
        "--runThreadN", str(threads),
        "--outSAMtype", "BAM", "Unsorted",
        "--outSAMunmapped", "Within",
        "--outSAMattributes", "NH", "HI", "AS", "nM", "MD", "XS",
        "--outSAMstrandField", "intronMotif",
        "--outFilterMultimapNmax", "20",
        "--outFilterType", "BySJout",
        "--alignSJoverhangMin", "8",
        "--alignSJDBoverhangMin", "1",
        "--outFileNamePrefix", str(prefix),
    ]
    t0 = time.perf_counter()
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logger.error("STAR failed (rc=%d):\n%s", res.returncode, res.stderr)
        raise RuntimeError("STAR alignment failed")
    dt = time.perf_counter() - t0
    logger.info("STAR done in %.1fs", dt)

    raw_bam = out_dir / "star_Aligned.out.bam"
    nsort_bam = out_dir / "aligned.nsort.bam"
    # pysam.sort -n is slow for large BAMs but deterministic; `samtools
    # sort -n -@ threads` would be faster, use it if available.
    pysam.sort("-n", "-@", str(threads), "-o", str(nsort_bam), str(raw_bam))
    raw_bam.unlink(missing_ok=True)
    return nsort_bam


# ---------------------------------------------------------------------------
# Runner (per sweep cell)
# ---------------------------------------------------------------------------

def run_one(
    out_dir: Path,
    *,
    genome: _StringGenome,
    transcripts: list[Transcript],
    index: TranscriptIndex,
    ss: float,
    gdna_fraction: float,
    seed: int,
    n_rna_fragments: int,
    rec_idx: int,
    star_bin: Path,
    star_index: Path,
    star_threads: int,
    keep_artifacts: bool = False,
) -> dict:
    label = f"ss{ss:.2f}_gdna{gdna_fraction:.2f}_seed{seed}"
    scen_dir = out_dir / "scenarios" / label
    scen_dir.mkdir(parents=True, exist_ok=True)

    sim_cfg = SimConfig(
        frag_mean=200, frag_std=40, frag_min=80, frag_max=500,
        read_length=100, strand_specificity=ss, seed=seed,
    )
    gdna_cfg = (
        GDNAConfig(abundance=10.0, frag_mean=350, frag_std=120,
                   frag_min=100, frag_max=1000)
        if gdna_fraction > 0 else None
    )

    sim = ReadSimulator(genome, transcripts, config=sim_cfg, gdna_config=gdna_cfg)

    n_mrna, n_nrna = sim.compute_rna_split(n_rna_fragments)
    n_gdna = int(round(n_rna_fragments * gdna_fraction)) if gdna_cfg else 0
    pool_split = (n_mrna, n_nrna, n_gdna)
    n_rna_truth = n_mrna + n_nrna

    # ---- simulate FASTQ --------------------------------------------------
    t0 = time.perf_counter()
    r1, r2, _ = write_fastq_fixed_split(sim, scen_dir, pool_split)
    t_sim = time.perf_counter() - t0

    # ---- align with STAR -------------------------------------------------
    t1 = time.perf_counter()
    bam_path = align_star(
        r1, r2, scen_dir,
        star_bin=star_bin, star_index=star_index, threads=star_threads,
    )
    t_star = time.perf_counter() - t1

    # Count mapped/unmapped pairs (fragment-level sanity).
    n_mapped_pairs = 0
    n_unmapped_pairs = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam:
            if rec.is_secondary or rec.is_supplementary:
                continue
            if not rec.is_read1:
                continue
            if rec.is_unmapped and rec.mate_is_unmapped:
                n_unmapped_pairs += 1
            else:
                n_mapped_pairs += 1

    # ---- rigel pipeline --------------------------------------------------
    lam_truth = truth_lambda_gdna(index, n_gdna)
    cfg = PipelineConfig(
        em=EMConfig(seed=seed),
        scan=BamScanConfig(sj_strand_tag="XS"),
    )
    t2 = time.perf_counter()
    try:
        pr = run_pipeline(str(bam_path), index, config=cfg)
        ok = True
        err = ""
    except Exception as exc:  # noqa: BLE001
        ok = False
        err = repr(exc)
        pr = None
    t_pipe = time.perf_counter() - t2

    row: dict = {
        "label": label,
        "ss": ss,
        "gdna_fraction": gdna_fraction,
        "seed": seed,
        "n_rna_truth": n_rna_truth,
        "n_gdna_truth": n_gdna,
        "n_total_truth_pairs": n_rna_truth + n_gdna,
        "n_mapped_pairs": n_mapped_pairs,
        "n_unmapped_pairs": n_unmapped_pairs,
        "mapping_rate": (
            n_mapped_pairs / max(1, n_mapped_pairs + n_unmapped_pairs)
        ),
        "lambda_truth": lam_truth,
        "t_sim_s": t_sim,
        "t_star_s": t_star,
        "t_pipe_s": t_pipe,
        "ok": ok,
        "err": err,
        "rec_idx": rec_idx,
    }

    if pr is not None:
        cal = pr.calibration
        est = pr.estimator
        row.update({
            "ss_estimated": None if cal is None else float(cal.strand_specificity),
            "lambda_pool": None if cal is None else cal.lambda_gdna,
            "lambda_density": None if cal is None else cal.lambda_gdna_density,
            "lambda_strand": None if cal is None else cal.lambda_gdna_strand,
            "phi_density": None if cal is None else cal.phi_density,
            "phi_strand": None if cal is None else cal.phi_strand,
            "fisher_density": None if cal is None else cal.fisher_info_density,
            "fisher_strand": None if cal is None else cal.fisher_info_strand,
            "n_eff_density": None if cal is None else cal.n_admitted_eff_density,
            "n_eff_strand": None if cal is None else cal.n_admitted_eff_strand,
            "consistency_chi2": None if cal is None else cal.consistency_chi2,
            "gdna_em_count": float(est.gdna_em_count),
            "nrna_em_count": float(est.nrna_em_count),
            "mrna_em_count": float(est.mrna_em_count) if hasattr(est, "mrna_em_count") else None,
            "gdna_contamination_rate": float(est.gdna_contamination_rate),
            "n_total_fragments": None if cal is None else float(cal.region_n_total.sum()),
        })
        import math
        lp = row["lambda_pool"]
        lt = row["lambda_truth"]
        if lp is not None and lt is not None and lt > 0:
            row["lambda_relative_error"] = (lp - lt) / lt
            row["lambda_log2_ratio"] = math.log2(
                max(lp, 1e-30) / max(lt, 1e-30)
            )
        else:
            row["lambda_relative_error"] = None
            row["lambda_log2_ratio"] = None
        if n_gdna > 0:
            row["gdna_em_relative_error"] = (
                (row["gdna_em_count"] - n_gdna) / n_gdna
            )
        else:
            row["gdna_em_relative_error"] = None

    if not keep_artifacts:
        # Keep the BAM path for debugging unless caller wants it gone;
        # FASTQ and STAR artefacts are always removed.
        for p in (r1, r2):
            p.unlink(missing_ok=True)
        for logf in scen_dir.glob("star_Log*"):
            logf.unlink(missing_ok=True)
        for sj in scen_dir.glob("star_SJ*"):
            sj.unlink(missing_ok=True)
        shutil.rmtree(scen_dir / "star__STARtmp", ignore_errors=True)
        bam_path.unlink(missing_ok=True)

    return row


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--star-index", type=Path, default=DEFAULT_STAR_INDEX)
    ap.add_argument("--star-bin", type=Path, default=DEFAULT_STAR_BIN)
    ap.add_argument("--star-threads", type=int, default=8)
    ap.add_argument("--n-rna-fragments", type=int, default=150_000)
    ap.add_argument("--max-transcripts", type=int, default=1500)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--keep-artifacts", action="store_true")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    args.out.mkdir(parents=True, exist_ok=True)

    if args.quick:
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.1, 0.5, 2.0]
        seeds = [1]
        n_rna = max(50_000, args.n_rna_fragments // 4)
    else:
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.05, 0.25, 1.0, 2.0]
        seeds = [1, 2]
        n_rna = args.n_rna_fragments

    n_cells = len(ss_grid) * len(gdna_grid) * len(seeds)
    logger.info("chr22 STAR calibration sweep — %d cells", n_cells)
    logger.info("  STAR index = %s", args.star_index)
    logger.info("  n_rna/cell = %d", n_rna)

    logger.info("Loading rigel index from %s ...", args.index)
    index = TranscriptIndex.load(args.index)
    genome = load_chr22_genome(args.fasta)

    panels: dict[int, list[Transcript]] = {}
    for s in seeds:
        panels[s] = select_transcripts(
            index, seed=s, max_transcripts=args.max_transcripts,
        )

    rows: list[dict] = []
    rec_idx = 0
    for ss in ss_grid:
        for gd in gdna_grid:
            for s in seeds:
                logger.info(
                    "[%d/%d] ss=%.2f gdna=%.2f seed=%d",
                    rec_idx + 1, n_cells, ss, gd, s,
                )
                row = run_one(
                    args.out,
                    genome=genome,
                    transcripts=panels[s],
                    index=index,
                    ss=ss, gdna_fraction=gd, seed=s,
                    n_rna_fragments=n_rna,
                    rec_idx=rec_idx,
                    star_bin=args.star_bin,
                    star_index=args.star_index,
                    star_threads=args.star_threads,
                    keep_artifacts=args.keep_artifacts,
                )
                rows.append(row)
                rec_idx += 1
                # Flush after each cell.
                pd.DataFrame(rows).to_csv(
                    args.out / "results.tsv", sep="\t", index=False,
                )

    df = pd.DataFrame(rows)
    df.to_csv(args.out / "results.tsv", sep="\t", index=False)
    (args.out / "results.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    logger.info("Wrote %d rows → %s", len(df), args.out / "results.tsv")


if __name__ == "__main__":
    main()
