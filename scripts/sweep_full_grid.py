#!/usr/bin/env python
"""Full hyperparameter sweep: gDNA × nRNA × mRNA × SS.

20 kb genome, two two-exon transcripts:
  t1 (+) exons (2000,4000) (8000,10000)  abundance = mRNA level (swept)
  t2 (-) exons (14000,16000) (18000,19000) abundance = 0 (always negative control)

Grid: 6 × 6 × 6 × 4 = 864 cases, 10 000 fragments each.
Outputs a TSV to stdout (pipe to file) with per-case metrics.
"""

import sys
import time
import tempfile
import pathlib
import logging

import numpy as np

from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logging.basicConfig(level=logging.WARNING, stream=sys.stderr)
logger = logging.getLogger(__name__)

# ── Grid ────────────────────────────────────────────────────────────
GDNA_LEVELS  = [0, 20, 50, 100, 200, 500]
NRNA_LEVELS  = [0, 20, 50, 100, 200, 500]
MRNA_LEVELS  = [0, 20, 50, 100, 200, 500]
SS_LEVELS    = [0.5, 0.9, 0.95, 1.0]

N_FRAGMENTS  = 10_000
SIM_SEED     = 42
PIPELINE_SEED = 42
GENOME_LEN   = 20_000

# ── Header ──────────────────────────────────────────────────────────
COLUMNS = [
    "gdna", "nrna", "mrna", "ss",
    # Fragment counts (ground truth)
    "n_frags_pipeline", "n_intergenic", "n_chimeric",
    "mrna_gt", "nrna_gt", "gdna_gt",
    # Pipeline estimates
    "mrna_pipe", "nrna_pipe", "gdna_pipe",
    # t1 transcript detail
    "t1_expected", "t1_observed", "t1_abs_diff", "t1_rel_err",
    # t2 negative control
    "t2_expected", "t2_observed",
    # Derived totals
    "total_rna_gt", "total_rna_pipe", "total_rna_rel_err",
    # Absolute pool errors
    "mrna_abs_err", "nrna_abs_err", "gdna_abs_err",
    # Relative pool errors (signed: positive = over-estimate)
    "mrna_rel_err", "nrna_rel_err", "gdna_rel_err",
]


def make_scenario(work_dir, mrna_abundance):
    """Create a fresh Scenario with the given mRNA abundance for t1."""
    sc = Scenario(
        "sweep",
        genome_length=GENOME_LEN,
        seed=SIM_SEED,
        work_dir=work_dir,
    )
    sc.add_gene("g1", "+", [
        {
            "t_id": "t1",
            "exons": [(2000, 4000), (8000, 10000)],
            "abundance": mrna_abundance,
        },
    ])
    sc.add_gene("g2", "-", [
        {
            "t_id": "t2",
            "exons": [(14000, 16000), (18000, 19000)],
            "abundance": 0,
        },
    ])
    return sc


def signed_rel(observed, expected):
    """Signed relative error: positive = over-estimate."""
    if expected == 0:
        return float("nan") if observed == 0 else float("inf")
    return (observed - expected) / expected


def run_one(gdna, nrna, mrna, ss, tmpdir):
    """Run a single grid point and return a dict of metrics."""
    label = f"g{gdna}_n{nrna}_m{mrna}_s{int(ss*100)}"
    work_dir = tmpdir / label

    sc = make_scenario(work_dir, mrna_abundance=mrna)

    sim_cfg = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=ss, seed=SIM_SEED,
    )
    gdna_cfg = GDNAConfig(
        abundance=gdna, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    ) if gdna > 0 else None

    # Skip degenerate case: no signal at all
    if mrna == 0 and nrna == 0 and gdna == 0:
        sc.cleanup()
        return None

    result = sc.build_oracle(
        n_fragments=N_FRAGMENTS,
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        nrna_abundance=nrna if nrna > 0 else 0.0,
    )

    cfg = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    bench = run_benchmark(result, pr, scenario_name=label)

    # Per-transcript
    t1 = next((t for t in bench.transcripts if t.t_id == "t1"), None)
    t2 = next((t for t in bench.transcripts if t.t_id == "t2"), None)

    # Pool totals
    mrna_gt = bench.total_expected
    nrna_gt = bench.n_nrna_expected
    gdna_gt = bench.n_gdna_expected
    mrna_pipe = bench.total_observed
    nrna_pipe = bench.n_nrna_pipeline
    gdna_pipe = bench.n_gdna_pipeline

    total_rna_gt = mrna_gt + nrna_gt
    total_rna_pipe = mrna_pipe + nrna_pipe
    total_rna_rel = signed_rel(total_rna_pipe, total_rna_gt)

    row = {
        "gdna": gdna, "nrna": nrna, "mrna": mrna, "ss": ss,
        "n_frags_pipeline": bench.n_fragments,
        "n_intergenic": bench.n_intergenic,
        "n_chimeric": bench.n_chimeric,
        "mrna_gt": mrna_gt, "nrna_gt": nrna_gt, "gdna_gt": gdna_gt,
        "mrna_pipe": f"{mrna_pipe:.1f}",
        "nrna_pipe": f"{nrna_pipe:.1f}",
        "gdna_pipe": f"{gdna_pipe:.1f}",
        "t1_expected": t1.expected if t1 else 0,
        "t1_observed": f"{t1.observed:.1f}" if t1 else "0.0",
        "t1_abs_diff": f"{t1.abs_diff:.1f}" if t1 else "0.0",
        "t1_rel_err": f"{t1.rel_error:.4f}" if t1 else "nan",
        "t2_expected": t2.expected if t2 else 0,
        "t2_observed": f"{t2.observed:.1f}" if t2 else "0.0",
        "total_rna_gt": total_rna_gt,
        "total_rna_pipe": f"{total_rna_pipe:.1f}",
        "total_rna_rel_err": f"{total_rna_rel:.4f}",
        "mrna_abs_err": f"{mrna_pipe - mrna_gt:.1f}",
        "nrna_abs_err": f"{nrna_pipe - nrna_gt:.1f}",
        "gdna_abs_err": f"{gdna_pipe - gdna_gt:.1f}",
        "mrna_rel_err": f"{signed_rel(mrna_pipe, mrna_gt):.4f}",
        "nrna_rel_err": f"{signed_rel(nrna_pipe, nrna_gt):.4f}",
        "gdna_rel_err": f"{signed_rel(gdna_pipe, gdna_gt):.4f}",
    }

    sc.cleanup()
    return row


def main():
    total = len(GDNA_LEVELS) * len(NRNA_LEVELS) * len(MRNA_LEVELS) * len(SS_LEVELS)
    print("\t".join(COLUMNS), flush=True)

    t0 = time.time()
    done = 0
    skipped = 0

    with tempfile.TemporaryDirectory(prefix="hulkrna_sweep_") as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        for gdna in GDNA_LEVELS:
            for nrna in NRNA_LEVELS:
                for mrna in MRNA_LEVELS:
                    for ss in SS_LEVELS:
                        row = run_one(gdna, nrna, mrna, ss, tmpdir)
                        done += 1
                        if row is None:
                            skipped += 1
                            continue
                        vals = [str(row[c]) for c in COLUMNS]
                        print("\t".join(vals), flush=True)

                        if done % 50 == 0:
                            elapsed = time.time() - t0
                            rate = done / elapsed
                            eta = (total - done) / rate if rate > 0 else 0
                            print(
                                f"# Progress: {done}/{total} "
                                f"({done/total:.0%}) "
                                f"elapsed={elapsed:.0f}s "
                                f"ETA={eta:.0f}s",
                                file=sys.stderr, flush=True,
                            )

    elapsed = time.time() - t0
    print(
        f"# Done: {done} cases ({skipped} skipped) in {elapsed:.1f}s",
        file=sys.stderr, flush=True,
    )


if __name__ == "__main__":
    main()
