"""Phase-2 calibration stress sweep (M6).

Generates a multi-gene synthetic genome, then sweeps:

  * strand specificity (SS)   ∈ {0.5, 0.9, 0.99, 1.0}
  * gDNA fraction-of-RNA      ∈ {0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0}
  * RNG seed                  ∈ {1, 2, 3}    (replicates)

For each (SS, gdna_fraction, seed) point:

  1. build an oracle BAM with controlled gDNA contamination,
  2. run the full rigel pipeline,
  3. extract the calibration estimates and the EM-allocated gDNA count,
  4. compute the *truth* λ_G from the simulated fragments.

Results are written to a TSV that the companion analyser
(``analyze_calibration_sweep.py``) summarises into the
``calibration_v4_baseline.md`` report.

Run::

    python scripts/benchmark/calibration_sweep.py \\
        --out scripts/benchmark/results/calibration_v4 \\
        [--quick]

``--quick`` runs only the fastest subset (1 seed, 4 gdna levels,
3 SS levels) for development checks.  Full sweep is the default.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.calibration import calibrate_gdna
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, Scenario, SimConfig

logger = logging.getLogger("calibration_sweep")


# ------------------------------------------------------------------ genome


@dataclass(frozen=True)
class GeneSpec:
    g_id: str
    strand: str
    transcripts: list[dict]


def build_mini_genome(n_genes: int = 25,
                      gene_pitch: int = 80_000,
                      seed: int = 0) -> tuple[int, list[GeneSpec]]:
    """A 2 Mb mini-genome with 25 genes / ~40 transcripts.

    Each gene block is 80 kb wide; the transcript spans the first
    ~6 kb of the block, leaving ~74 kb of intergenic space — plenty
    of headroom for the calibration to see pure-gDNA regions.

    Exons are intentionally short (150 bp) so that fragments
    (mean 200 bp) routinely span exon-exon junctions, giving the
    strand model enough spliced unique-mapper observations to
    train an SS estimate.

    Roughly half the genes are multi-isoform (two transcripts
    sharing some exons), the rest single-iso.  Strand alternates +/-.
    """
    rng = np.random.default_rng(seed)
    genes: list[GeneSpec] = []
    for i in range(n_genes):
        block = i * gene_pitch
        strand = "+" if i % 2 == 0 else "-"
        # Six short exons, ~150 bp each, ~700 bp introns.
        ex_starts = [block + 1_000 + k * 850 for k in range(6)]
        exons = [(s, s + 150) for s in ex_starts]
        # vary expression by 10x within a sweep so the algorithm
        # sees realistic per-region count heterogeneity.
        base = float(rng.integers(20, 200))
        if i % 2 == 0:
            # multi-isoform: full transcript + a 3-exon variant
            txs = [
                {"t_id": f"g{i}_t1",
                 "exons": exons,
                 "abundance": base},
                {"t_id": f"g{i}_t2",
                 "exons": [exons[0], exons[2], exons[5]],
                 "abundance": base * 0.4},
            ]
        else:
            txs = [
                {"t_id": f"g{i}_t1",
                 "exons": [exons[0], exons[2], exons[3], exons[5]],
                 "abundance": base},
            ]
        genes.append(GeneSpec(g_id=f"g{i}", strand=strand, transcripts=txs))

    genome_length = n_genes * gene_pitch
    return genome_length, genes


# ------------------------------------------------------------------ truth


def truth_lambda_gdna(scenario_result, n_gdna_truth: int) -> float:
    """Compute the *true* λ_G from the simulator output.

    ``λ_G = (# unspliced gDNA fragments overlapping admissible regions)
             / (Σ E_i over admissible regions)``.

    For the synthetic indexes (no mappability mask), the index region
    table covers the full genome with E_i = region length, so a clean
    upper-bound proxy is ``n_gdna_truth / sum(region_length)``.

    This is what the calibration's λ̂ aims to recover.
    """
    region_df = scenario_result.index.region_df
    if "mappable_effective_length" in region_df.columns:
        total_E = float(region_df["mappable_effective_length"].sum())
    else:
        total_E = float(region_df["length"].sum())
    if total_E <= 0:
        return float("nan")
    return float(n_gdna_truth) / total_E


# ------------------------------------------------------------------ runner


def run_one(out_dir: Path,
            ss: float,
            gdna_fraction: float,
            seed: int,
            n_rna_fragments: int,
            n_genes: int,
            rec_idx: int) -> dict:
    """Build one scenario, run the pipeline, return a result row."""
    label = f"ss{ss:.2f}_gdna{gdna_fraction:.2f}_seed{seed}"
    work_dir = out_dir / "scenarios" / label
    work_dir.mkdir(parents=True, exist_ok=True)

    genome_length, genes = build_mini_genome(n_genes=n_genes, seed=seed)
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

    # Truth bookkeeping
    n_rna_truth = sum(result.ground_truth_auto().values())
    n_gdna_truth = result.ground_truth_gdna_auto()
    lam_truth = truth_lambda_gdna(result, n_gdna_truth)

    # Pipeline
    # Oracle BAM writes the STAR-convention XS tag on spliced reads, so we
    # point the scanner at XS explicitly (the default "auto" would also work).
    cfg = PipelineConfig(
        em=EMConfig(seed=seed),
        scan=BamScanConfig(sj_strand_tag="XS"),
    )
    t1 = time.perf_counter()
    try:
        pr = run_pipeline(result.bam_path, result.index, config=cfg)
        ok = True
        err = ""
    except Exception as exc:  # noqa: BLE001 — capture any pipeline failure
        ok = False
        err = repr(exc)
        pr = None
    t_pipe = time.perf_counter() - t1

    row: dict = {
        "label": label,
        "rec_idx": rec_idx,
        "ss": ss,
        "gdna_fraction": gdna_fraction,
        "seed": seed,
        "n_genes": n_genes,
        "genome_length": genome_length,
        "n_rna_target": n_rna_fragments,
        "n_rna_truth": n_rna_truth,
        "n_gdna_truth": n_gdna_truth,
        "lambda_truth": lam_truth,
        "t_sim_s": t_sim,
        "t_pipe_s": t_pipe,
        "ok": ok,
        "error": err,
    }

    if pr is not None:
        cal = pr.calibration
        est = pr.estimator
        row.update({
            "ss_estimated": (None if cal is None
                             else float(cal.strand_specificity)),
            "lambda_pool": None if cal is None else cal.lambda_gdna,
            "lambda_density": (None if cal is None
                               else cal.lambda_gdna_density),
            "lambda_strand": (None if cal is None
                              else cal.lambda_gdna_strand),
            "phi_density": None if cal is None else cal.phi_density,
            "phi_strand": None if cal is None else cal.phi_strand,
            "fisher_density": (None if cal is None
                               else cal.fisher_info_density),
            "fisher_strand": (None if cal is None
                              else cal.fisher_info_strand),
            "n_eff_density": (None if cal is None
                              else cal.n_admitted_eff_density),
            "n_eff_strand": (None if cal is None
                             else cal.n_admitted_eff_strand),
            "consistency_chi2": (None if cal is None
                                 else cal.consistency_chi2),
            "gdna_em_count": float(est.gdna_em_count),
            "nrna_em_count": float(est.nrna_em_count),
            "gdna_contamination_rate": float(est.gdna_contamination_rate),
            "n_total_fragments": (None if cal is None
                                  else float(cal.region_n_total.sum())),
        })
        # Derived metrics on the pipeline-realistic pooled λ̂.
        lp = row["lambda_pool"]
        lt = row["lambda_truth"]
        if lp is not None and lt is not None and lt > 0:
            row["lambda_relative_error"] = (lp - lt) / lt
            row["lambda_log2_ratio"] = math.log2(max(lp, 1e-30)
                                                  / max(lt, 1e-30))
        else:
            row["lambda_relative_error"] = None
            row["lambda_log2_ratio"] = None
        if n_gdna_truth > 0:
            row["gdna_em_relative_error"] = (
                (row["gdna_em_count"] - n_gdna_truth) / n_gdna_truth
            )
        else:
            row["gdna_em_relative_error"] = None

        # ---- Truth-SS calibration --------------------------------
        # The synthetic oracle BAM does not emit XS/ts tags rich
        # enough to train an SS estimate from spliced reads, so the
        # pipeline-estimated SS collapses to 0.5 and the strand
        # pathway is masked out.  Run the calibration once more
        # *with the simulated SS injected as ground truth* to give
        # the strand pathway a chance to fire — this isolates strand-
        # pathway behaviour from the orthogonal SS-estimation issue.
        try:
            _, _, frag_models2, _, region_counts2, fl_table2 = scan_and_buffer(
                str(result.bam_path), result.index, cfg.scan,
            )
            frag_models2.build_scoring_models()
            frag_models2.finalize(prior_ess=cfg.calibration.fl_prior_ess)
            cal_truth = calibrate_gdna(
                region_counts2, fl_table2, result.index.region_df,
                strand_specificity=ss,
                mean_frag_len=frag_models2.global_model.mean,
                intergenic_fl_model=frag_models2.intergenic,
                fl_prior_ess=cfg.calibration.fl_prior_ess,
            )
            row.update({
                "truthss_lambda_pool": cal_truth.lambda_gdna,
                "truthss_lambda_density": cal_truth.lambda_gdna_density,
                "truthss_lambda_strand": cal_truth.lambda_gdna_strand,
                "truthss_phi_density": cal_truth.phi_density,
                "truthss_phi_strand": cal_truth.phi_strand,
                "truthss_fisher_density": cal_truth.fisher_info_density,
                "truthss_fisher_strand": cal_truth.fisher_info_strand,
                "truthss_n_eff_density": cal_truth.n_admitted_eff_density,
                "truthss_n_eff_strand": cal_truth.n_admitted_eff_strand,
                "truthss_consistency_chi2": cal_truth.consistency_chi2,
            })
            tlp = cal_truth.lambda_gdna
            if tlp is not None and lt > 0:
                row["truthss_lambda_relative_error"] = (tlp - lt) / lt
            tls = cal_truth.lambda_gdna_strand
            if tls is not None and lt > 0:
                row["truthss_strand_relative_error"] = (tls - lt) / lt
        except Exception as exc:  # noqa: BLE001
            row["truthss_error"] = repr(exc)

    sc.cleanup()
    return row


# ------------------------------------------------------------------ sweep


def make_grid(quick: bool) -> list[tuple[float, float, int]]:
    if quick:
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.1, 0.5, 2.0]
        seeds = [1]
    else:
        ss_grid = [0.5, 0.9, 0.99, 1.0]
        gdna_grid = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0]
        seeds = [1, 2, 3]
    return [(ss, g, s) for ss in ss_grid for g in gdna_grid for s in seeds]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", type=Path, required=True,
                        help="Output directory for TSV + per-run JSON.")
    parser.add_argument("--n-rna-fragments", type=int, default=50_000,
                        help="RNA fragments per scenario (default 50k).")
    parser.add_argument("--n-genes", type=int, default=25,
                        help="Genes in the mini-genome (default 25).")
    parser.add_argument("--quick", action="store_true",
                        help="Run a much smaller sweep for development.")
    parser.add_argument("--filter-ss", type=float, default=None,
                        help="Run only this SS value (debugging).")
    parser.add_argument("--filter-gdna", type=float, default=None,
                        help="Run only this gdna_fraction value.")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    args.out.mkdir(parents=True, exist_ok=True)
    grid = make_grid(args.quick)
    if args.filter_ss is not None:
        grid = [g for g in grid if abs(g[0] - args.filter_ss) < 1e-6]
    if args.filter_gdna is not None:
        grid = [g for g in grid if abs(g[1] - args.filter_gdna) < 1e-6]

    print(f"Running {len(grid)} configurations...", flush=True)
    rows = []
    for k, (ss, gf, seed) in enumerate(grid):
        t = time.perf_counter()
        row = run_one(args.out, ss=ss, gdna_fraction=gf, seed=seed,
                      n_rna_fragments=args.n_rna_fragments,
                      n_genes=args.n_genes, rec_idx=k)
        dt = time.perf_counter() - t
        lam_str = (f"λ̂={row.get('lambda_pool'):.3e} "
                   f"truth={row['lambda_truth']:.3e} "
                   f"err={row.get('lambda_relative_error') or 0:+.2%}"
                   if row.get("ok") and row.get("lambda_pool") is not None
                   else "λ̂ = n/a")
        print(f"  [{k+1:3d}/{len(grid)}] ss={ss:.2f} gdna={gf:.2f} "
              f"seed={seed} n_rna={row['n_rna_truth']} "
              f"n_gdna={row['n_gdna_truth']}  {lam_str}  "
              f"({dt:.1f}s)", flush=True)
        rows.append(row)
        # Save progressively so we can analyse partial results.
        df = pd.DataFrame(rows)
        df.to_csv(args.out / "results.tsv", sep="\t", index=False)
        with open(args.out / "results.json", "w") as f:
            json.dump(rows, f, indent=2, default=str)

    print(f"\nDone. Results: {args.out / 'results.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
