#!/usr/bin/env python3
"""Generate VCaP benchmark simulation data with three conditions.

Conditions (additive model):
  1. pristine    — 10M mRNA, 0 gDNA, 0 nRNA  (10M total)
  2. gdna        — 10M mRNA, 10M gDNA, 0 nRNA (20M total)
  3. gdna_nrna   — 10M mRNA, 10M gDNA, 5M nRNA (25M total)

mRNA abundances come from salmon quant.sf (VCaP prostate).
nRNA is spiked at the abundance level, but fragment counts for
the gdna_nrna condition are controlled independently via separate
mRNA/nRNA pool sampling.

Usage:
    conda activate rigel
    python scripts/sim/run_vcap_benchmark.py
"""
from __future__ import annotations

import copy
import json
import logging
import sys
import time
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.sim.sim import (
    GDNASimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    _spike_in_nrna,
    assign_file_abundances,
    load_transcripts,
    write_truth_abundances,
)

# ── Configuration ─────────────────────────────────────────────────

GENOME = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/genome_controls.fasta.bgz"
GTF = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/genes_controls.gtf.gz"
SALMON_TSV = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/ccle_vcap_prostate/salmon/quant.sf.gz"
OUTDIR = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/vcap_sim")

TRANSCRIPT_FILTER = "basic"
SEED = 42
STRAND_SPECIFICITY = 1.0  # perfectly strand-specific

SIM_PARAMS = SimulationParams(
    n_rna_fragments=10_000_000,  # placeholder, overridden per condition
    sim_seed=SEED,
    frag_mean=250.0,
    frag_std=50.0,
    frag_min=50,
    frag_max=1000,
    read_length=150,
    error_rate=0.0,
)

GDNA_CONFIG = GDNASimConfig(
    rates=[],  # not used directly
    frag_mean=350.0,
    frag_std=100.0,
    frag_min=100,
    frag_max=1000,
)

# Conditions: list of dicts for clarity
#   n_rna:  total RNA fragments (used in combined-pool mode)
#   n_mrna / n_nrna: if both set, pools are sampled independently
#   nrna_frac_range: abundance-level nRNA spike-in fraction
CONDITIONS = [
    {"name": "pristine", "n_rna": 10_000_000, "n_gdna": 0},
    {"name": "gdna",     "n_rna": 10_000_000, "n_gdna": 10_000_000},
    {
        "name": "gdna_nrna",
        "n_mrna": 10_000_000,
        "n_nrna": 5_000_000,
        "n_gdna": 10_000_000,
        "nrna_frac_range": (0.333, 0.333),
    },
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Load transcripts and assign salmon abundances ─────────────
    logger.info("Loading transcripts...")
    transcripts = load_transcripts(GTF, transcript_filter=TRANSCRIPT_FILTER)
    if not transcripts:
        print("ERROR: No transcripts loaded", file=sys.stderr)
        return 1

    logger.info("Assigning salmon abundances from %s", SALMON_TSV)
    assign_file_abundances(transcripts, SALMON_TSV)

    # Save base abundances (all mRNA, no nRNA yet)
    base_abundances = [(t.abundance or 0.0, 0.0) for t in transcripts]

    # ── Run each condition ────────────────────────────────────────
    manifest_conditions = []
    total_t0 = time.monotonic()

    for cond in CONDITIONS:
        cond_name = cond["name"]
        n_gdna = cond.get("n_gdna", 0)
        n_mrna = cond.get("n_mrna")
        n_nrna = cond.get("n_nrna")
        n_rna = cond.get("n_rna", (n_mrna or 0) + (n_nrna or 0))
        nrna_frac_range = cond.get("nrna_frac_range")
        cond_dir = OUTDIR / cond_name

        # Display fragment counts
        if n_mrna is not None and n_nrna is not None:
            disp_mrna, disp_nrna = n_mrna, n_nrna
            n_total = n_mrna + n_nrna + n_gdna
        else:
            disp_mrna, disp_nrna = n_rna, 0
            n_total = n_rna + n_gdna

        print(
            f"\n{'='*60}\n"
            f"Condition: {cond_name}\n"
            f"  mRNA:  {disp_mrna:,}\n"
            f"  nRNA:  {disp_nrna:,}\n"
            f"  gDNA:  {n_gdna:,}\n"
            f"  Total: ~{n_total:,}\n"
            f"{'='*60}",
            flush=True,
        )

        # Skip if already generated
        if (cond_dir / "sim_R1.fq.gz").exists():
            print(f"  Output exists at {cond_dir}, skipping", flush=True)
            manifest_conditions.append({
                "name": cond_name,
                "n_mrna": n_mrna, "n_nrna": n_nrna,
                "n_rna": n_rna, "n_gdna": n_gdna,
                "nrna_frac_range": list(nrna_frac_range) if nrna_frac_range else None,
                "strand_specificity": STRAND_SPECIFICITY,
            })
            continue

        # Deep-copy transcripts and restore base abundances
        cond_transcripts = copy.deepcopy(transcripts)
        for t, (base_mrna, _) in zip(cond_transcripts, base_abundances):
            t.abundance = base_mrna
            t.nrna_abundance = 0.0

        # Spike in nRNA if requested
        if nrna_frac_range is not None:
            _spike_in_nrna(cond_transcripts, nrna_frac_range, seed=SEED)

        # Write truth abundances
        truth_path = OUTDIR / f"truth_abundances_{cond_name}.tsv"
        write_truth_abundances(cond_transcripts, truth_path)

        # Build and run simulator
        t0 = time.monotonic()
        sim = WholeGenomeSimulator(
            GENOME,
            cond_transcripts,
            SIM_PARAMS,
            GDNA_CONFIG,
            strand_specificity=STRAND_SPECIFICITY,
            seed=SEED,
        )
        r1, r2, bam = sim.simulate_and_write(
            cond_dir, n_rna, n_gdna,
            n_mrna=n_mrna, n_nrna=n_nrna,
            oracle_bam=True, prefix="sim",
        )
        sim.close()
        elapsed = time.monotonic() - t0

        print(
            f"  Done in {elapsed:.1f}s\n"
            f"  R1: {r1}\n"
            f"  R2: {r2}\n"
            f"  BAM: {bam}",
            flush=True,
        )

        manifest_conditions.append({
            "name": cond_name,
            "n_mrna": n_mrna, "n_nrna": n_nrna,
            "n_rna": n_rna, "n_gdna": n_gdna,
            "nrna_frac_range": list(nrna_frac_range) if nrna_frac_range else None,
            "strand_specificity": STRAND_SPECIFICITY,
            "fastq_r1": str(r1),
            "fastq_r2": str(r2),
            "oracle_bam": str(bam) if bam else None,
            "truth_abundances": str(truth_path),
            "elapsed_seconds": round(elapsed, 1),
        })

    # ── Write manifest ────────────────────────────────────────────
    manifest = {
        "version": 1,
        "genome": GENOME,
        "gtf": GTF,
        "salmon_abundances": SALMON_TSV,
        "transcript_filter": TRANSCRIPT_FILTER,
        "simulation": asdict(SIM_PARAMS),
        "gdna": asdict(GDNA_CONFIG),
        "conditions": manifest_conditions,
    }
    manifest_path = OUTDIR / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    total_elapsed = time.monotonic() - total_t0
    print(
        f"\n{'='*60}\n"
        f"All conditions complete in {total_elapsed:.1f}s\n"
        f"Output: {OUTDIR}\n"
        f"Manifest: {manifest_path}\n"
        f"{'='*60}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
