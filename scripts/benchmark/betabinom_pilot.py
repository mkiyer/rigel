"""Beta-Binomial strand-LLR pilot on the mini-stress grid.

Runs the pipeline twice per (ss, gdna_fraction, seed) cell:

  1. strand_llr_mode = "binomial"  — current shipping behaviour.
  2. strand_llr_mode = "betabinom" — v3-style shared κ.

Writes a tidy TSV with one row per (cell, mode) for downstream
comparison.  Reuses the mini-genome generator from
``scripts/benchmark/calibration_sweep.py`` but side-steps its stale
diagnostic paths (they reference calibration fields that no longer
exist on the v4 CalibrationResult).

Usage::

    python scripts/benchmark/betabinom_pilot.py \
        --out scripts/benchmark/results/betabinom_pilot \
        [--quick]

See ``docs/calibration/calibration_v3_vs_v4_comparison.md``.
"""
from __future__ import annotations

import argparse
import dataclasses as _dc
import json
import logging
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Allow importing the sibling harness for the truth-lambda helper only.
_THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS))
from calibration_sweep import (  # type: ignore[import-not-found]
    GeneSpec,
    truth_lambda_gdna,
)

from rigel.config import (
    BamScanConfig,
    CalibrationConfig,
    EMConfig,
    PipelineConfig,
)
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

logger = logging.getLogger("betabinom_pilot")


# ---------------------------------------------------------------------------
# Realistic mini-genome
# ---------------------------------------------------------------------------

def build_realistic_mini_genome(
    n_genes: int = 25,
    gene_pitch: int = 200_000,
    exon_len: int = 2_000,
    intron_len: int = 1_000,
    n_exons: int = 4,
    seed: int = 0,
) -> tuple[int, list[GeneSpec]]:
    """Long-exon, large-intergenic genome for strand-gate validation.

    Compared to ``build_mini_genome``, this variant widens exons from
    150 bp → 2 kb and gene_pitch from 80 kb → 200 kb so that:

    * the vast majority of RNA fragments (frag_mean 200 bp) fall
      *within* a single exon → unspliced-annotated reads are plentiful,
      which is what the aggregate-strand-z gate pools over;
    * the exonic footprint is ≤ 10 % of the genome
      (default: 25 genes × 4 × 2 kb = 200 kb / 5 Mb = 4 %),
      so gDNA contamination does not overrun the gene-strand bucket.

    Gene span is ``n_exons·exon_len + (n_exons−1)·intron_len`` —
    e.g. 4 × 2 kb + 3 × 1 kb = 11 kb, comfortably inside a 200 kb block.
    """
    rng = np.random.default_rng(seed)
    genes: list[GeneSpec] = []
    stride = exon_len + intron_len
    for i in range(n_genes):
        block = i * gene_pitch
        strand = "+" if i % 2 == 0 else "-"
        start0 = block + 10_000  # 10 kb intergenic lead
        exons = [(start0 + k * stride, start0 + k * stride + exon_len)
                 for k in range(n_exons)]
        base = float(rng.integers(20, 200))
        if i % 2 == 0:
            # multi-isoform: full + a skip-exon variant
            skip_idx = [0, 2, n_exons - 1] if n_exons >= 3 else [0, n_exons - 1]
            txs = [
                {"t_id": f"g{i}_t1", "exons": exons, "abundance": base},
                {"t_id": f"g{i}_t2",
                 "exons": [exons[k] for k in skip_idx],
                 "abundance": base * 0.4},
            ]
        else:
            txs = [{"t_id": f"g{i}_t1", "exons": exons, "abundance": base}]
        genes.append(GeneSpec(g_id=f"g{i}", strand=strand, transcripts=txs))
    genome_length = n_genes * gene_pitch
    return genome_length, genes


def _make_config(seed: int, strand_llr_mode: str) -> PipelineConfig:
    return PipelineConfig(
        em=EMConfig(seed=seed),
        scan=BamScanConfig(sj_strand_tag="XS"),
        calibration=CalibrationConfig(strand_llr_mode=strand_llr_mode),
    )


def _run_one(
    *,
    out_dir: Path,
    ss: float,
    gdna_fraction: float,
    seed: int,
    n_rna_fragments: int,
    n_genes: int,
    mode: str,
    shared_scenario_root: Path,
) -> dict:
    """Simulate the BAM once (cached per cell), run the pipeline once per mode."""
    label = f"ss{ss:.2f}_gdna{gdna_fraction:.2f}_seed{seed}"

    # One scenario BAM per cell, reused across modes.
    work_dir = shared_scenario_root / label
    work_dir.mkdir(parents=True, exist_ok=True)

    genome_length, genes = build_realistic_mini_genome(n_genes=n_genes, seed=seed)
    sc = Scenario(label, genome_length=genome_length, seed=seed,
                  work_dir=work_dir)
    for g in genes:
        sc.add_gene(g.g_id, g.strand, g.transcripts)

    sim_cfg = SimConfig(
        frag_mean=200, frag_std=40, frag_min=80, frag_max=500,
        read_length=100, strand_specificity=ss, seed=seed,
    )
    gdna_cfg = (GDNAConfig(abundance=10.0, frag_mean=350, frag_std=120,
                           frag_min=100, frag_max=1000)
                if gdna_fraction > 0 else None)

    t0 = time.perf_counter()
    result = sc.build_oracle(
        n_rna_fragments=n_rna_fragments,
        gdna_fraction=gdna_fraction,
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
    )
    t_sim = time.perf_counter() - t0

    n_rna_truth = sum(result.ground_truth_auto().values())
    n_gdna_truth = result.ground_truth_gdna_auto()
    lam_truth = truth_lambda_gdna(result, n_gdna_truth)

    cfg = _make_config(seed=seed, strand_llr_mode=mode)
    t1 = time.perf_counter()
    try:
        pr = run_pipeline(result.bam_path, result.index, config=cfg)
        ok = True
        err = ""
    except Exception as exc:  # noqa: BLE001
        ok = False
        err = repr(exc)
        pr = None
    t_pipe = time.perf_counter() - t1

    row: dict = {
        "label": label,
        "mode": mode,
        "ss": ss,
        "gdna_fraction": gdna_fraction,
        "seed": seed,
        "n_genes": n_genes,
        "genome_length": genome_length,
        "n_rna_target": n_rna_fragments,
        "n_rna_truth": int(n_rna_truth),
        "n_gdna_truth": int(n_gdna_truth),
        "lambda_truth": float(lam_truth),
        "t_sim_s": t_sim,
        "t_pipe_s": t_pipe,
        "ok": ok,
        "error": err,
    }

    if pr is not None and pr.calibration is not None:
        cal = pr.calibration
        est = pr.estimator
        lam_hat = float(cal.lambda_gdna)
        row.update({
            "ss_estimated": float(cal.strand_specificity),
            "lambda_hat": lam_hat,
            "lambda_ratio": (lam_hat / lam_truth) if lam_truth > 0 else float("nan"),
            "lambda_rel_err": ((lam_hat - lam_truth) / lam_truth)
                if lam_truth > 0 else float("nan"),
            "lambda_log2_ratio": math.log2(max(lam_hat, 1e-30) / max(lam_truth, 1e-30))
                if lam_truth > 0 else float("nan"),
            "mu_R": float(cal.mu_R) if cal.mu_R is not None else float("nan"),
            "sigma_R": float(cal.sigma_R) if cal.sigma_R is not None else float("nan"),
            "pi": float(cal.mixing_pi) if cal.mixing_pi is not None else float("nan"),
            "pi_soft": float(cal.mixing_pi_soft)
                if cal.mixing_pi_soft is not None else float("nan"),
            "strand_used": bool(cal.strand_used),
            "strand_z": float(cal.strand_z),
            "strand_llr_mode": str(cal.strand_llr_mode),
            "kappa": float(cal.kappa),
            "em_n_iter": int(cal.em_n_iter),
            "em_converged": bool(cal.em_converged),
            "n_eligible": int(cal.n_eligible),
            "n_soft": int(cal.n_soft),
            "n_spliced_hard": int(cal.n_spliced_hard),
            "total_e_gdna": float(cal.region_e_gdna.sum()),
            "gdna_em_count": float(est.gdna_em_count),
            "nrna_em_count": float(est.nrna_em_count),
            "gdna_contamination_rate": float(est.gdna_contamination_rate),
            "n_total_fragments": float(cal.region_n_total.sum()),
        })
        if n_gdna_truth > 0:
            row["gdna_em_relative_error"] = (
                (row["gdna_em_count"] - n_gdna_truth) / n_gdna_truth
            )
        else:
            row["gdna_em_relative_error"] = None

    sc.cleanup()
    return row


def make_grid(quick: bool) -> list[tuple[float, float, int]]:
    if quick:
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.1, 0.5, 2.0]
        seeds = [1]
    else:
        ss_grid = [0.5, 0.9, 0.99, 1.0]
        gdna_grid = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0]
        seeds = [1, 2, 3]
    return [(s, g, d) for s in ss_grid for g in gdna_grid for d in seeds]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", type=Path, required=True)
    parser.add_argument("--n-rna-fragments", type=int, default=50_000)
    parser.add_argument("--n-genes", type=int, default=25)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--filter-ss", type=float, default=None)
    parser.add_argument("--filter-gdna", type=float, default=None)
    parser.add_argument("--modes", type=str, default="binomial,betabinom",
                        help="Comma-separated list of modes to run.")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    scenarios_root = out_dir / "scenarios"
    scenarios_root.mkdir(parents=True, exist_ok=True)

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    grid = make_grid(args.quick)
    if args.filter_ss is not None:
        grid = [c for c in grid if c[0] == args.filter_ss]
    if args.filter_gdna is not None:
        grid = [c for c in grid if c[1] == args.filter_gdna]

    logger.info(
        "Pilot grid: %d cells × %d modes = %d runs",
        len(grid), len(modes), len(grid) * len(modes),
    )
    rows: list[dict] = []
    total = len(grid) * len(modes)
    i = 0
    for ss, gdna_frac, seed in grid:
        for mode in modes:
            i += 1
            logger.info(
                "[%d/%d] ss=%.2f  gdna=%.2f  seed=%d  mode=%s",
                i, total, ss, gdna_frac, seed, mode,
            )
            try:
                row = _run_one(
                    out_dir=out_dir,
                    ss=ss,
                    gdna_fraction=gdna_frac,
                    seed=seed,
                    n_rna_fragments=args.n_rna_fragments,
                    n_genes=args.n_genes,
                    mode=mode,
                    shared_scenario_root=scenarios_root,
                )
            except Exception as exc:  # noqa: BLE001
                logger.exception("Cell failed")
                row = {
                    "ss": ss, "gdna_fraction": gdna_frac, "seed": seed,
                    "mode": mode, "ok": False, "error": repr(exc),
                }
            rows.append(row)

    df = pd.DataFrame(rows)
    tsv_path = out_dir / "results.tsv"
    json_path = out_dir / "results.json"
    df.to_csv(tsv_path, sep="\t", index=False)
    json_path.write_text(json.dumps(rows, default=str, indent=2))
    logger.info("Wrote %s (%d rows)", tsv_path, len(rows))


if __name__ == "__main__":
    main()
