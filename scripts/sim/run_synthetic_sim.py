#!/usr/bin/env python3
"""Run the full synthetic mini-genome simulation pipeline.

Steps:
  1. Generate synthetic genome (10 Mb, ~50 genes, ~250 transcripts)
  2. Simulate reads across 10 conditions (2 strand × 5 gDNA)
  3. Each condition produces oracle BAM + paired-end FASTQ

Usage:
  python scripts/sim/run_synthetic_sim.py \\
    --outdir /Users/mkiyer/Downloads/rigel_runs/sim_synthetic
"""
from __future__ import annotations

import argparse
import importlib
import importlib.machinery
import importlib.util
import logging
import shutil
import sys
import time
from pathlib import Path

# Ensure scripts/ and src/ are importable
_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = _SCRIPT_DIR.parent.parent
sys.path.insert(0, str(_ROOT / "src"))
sys.path.insert(0, str(_SCRIPT_DIR))

from generate_synthetic_genome import (
    generate_genes,
    generate_random_genome,
    inject_splice_sites,
    write_fasta,
    write_gtf,
    GENOME_LENGTH,
    N_GENES,
    REF_NAME,
    SEED,
)
import numpy as np

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Full synthetic mini-genome simulation pipeline"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic"),
        help="Base output directory"
    )
    parser.add_argument(
        "--genome-length", type=int, default=GENOME_LENGTH,
        help=f"Genome length in bp (default: {GENOME_LENGTH:,})"
    )
    parser.add_argument(
        "--n-genes", type=int, default=N_GENES,
        help=f"Number of genes (default: {N_GENES})"
    )
    parser.add_argument(
        "--seed", type=int, default=SEED,
        help=f"Random seed (default: {SEED})"
    )
    parser.add_argument(
        "--n-rna", type=int, default=1_000_000,
        help="Number of RNA fragments per condition (default: 1,000,000)"
    )
    parser.add_argument(
        "-j", "--n-workers", type=int, default=1,
        help="Worker processes for parallel read generation"
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip conditions that already have output files"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    ref_dir = outdir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()

    # ══════════════════════════════════════════════════════════════════════
    # Step 1: Generate synthetic genome
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("STEP 1: Generating synthetic genome")
    print("=" * 70)

    fasta_path = ref_dir / "genome.fa"
    gtf_path = ref_dir / "genes.gtf"

    if fasta_path.exists() and gtf_path.exists() and args.skip_existing:
        print(f"  Genome already exists at {ref_dir}, skipping generation")
    else:
        genes = generate_genes(args.genome_length, args.n_genes, args.seed)
        n_transcripts = sum(len(g.transcripts) for g in genes)
        n_antisense = sum(1 for g in genes if g.gene_name.endswith("as"))

        # Generate genome sequence
        genome_seq = generate_random_genome(args.genome_length, args.seed)

        # Inject splice sites
        all_transcripts = [tx for g in genes for tx in g.transcripts]
        seq_array = np.frombuffer(genome_seq.encode("ascii"), dtype=np.uint8).copy()
        inject_splice_sites(seq_array, all_transcripts)
        genome_seq = seq_array.tobytes().decode("ascii")

        # Write files
        fasta_path = write_fasta(genome_seq, REF_NAME, ref_dir)
        gtf_path = write_gtf(genes, REF_NAME, ref_dir)

        iso_counts = [len(g.transcripts) for g in genes]
        n_multi_exon = sum(1 for tx in all_transcripts if len(tx.exons) > 1)

        print(f"  Genome:       {fasta_path} ({args.genome_length:,} bp)")
        print(f"  Annotation:   {gtf_path}")
        print(f"  Genes:        {len(genes)} ({n_antisense} antisense overlapping)")
        print(f"  Transcripts:  {n_transcripts} ({n_multi_exon} multi-exon)")
        print(f"  Isoforms/gene: {min(iso_counts)}-{max(iso_counts)} (mean {np.mean(iso_counts):.1f})")

    # ══════════════════════════════════════════════════════════════════════
    # Step 2: Run simulation across conditions
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("STEP 2: Simulating reads across conditions")
    print("=" * 70)

    # Import functions from sim.py (via sys.path already set up)
    # We need to work around the name collision with the 'sim' directory.
    # Temporarily rename the module file for import.
    import importlib
    sys.path.insert(0, str(_SCRIPT_DIR))
    # sim.py sits in _SCRIPT_DIR; give it a unique module name
    _sim_path = _SCRIPT_DIR / "sim.py"
    _loader = importlib.machinery.SourceFileLoader("whole_genome_sim", str(_sim_path))
    _spec = importlib.util.spec_from_loader("whole_genome_sim", _loader)
    _sim_mod = importlib.util.module_from_spec(_spec)
    sys.modules["whole_genome_sim"] = _sim_mod
    _loader.exec_module(_sim_mod)

    SimConfig = _sim_mod.SimConfig
    SimulationParams = _sim_mod.SimulationParams
    AbundanceConfig = _sim_mod.AbundanceConfig
    GDNASimConfig = _sim_mod.GDNASimConfig
    NRNAConfig = _sim_mod.NRNAConfig
    load_transcripts = _sim_mod.load_transcripts
    assign_random_abundances = _sim_mod.assign_random_abundances
    WholeGenomeSimulator = _sim_mod.WholeGenomeSimulator
    write_truth_abundances = _sim_mod.write_truth_abundances
    condition_dir_name = _sim_mod.condition_dir_name
    gdna_label_for_rate = _sim_mod.gdna_label_for_rate
    write_manifest = _sim_mod.write_manifest
    import copy
    import json

    # Configuration
    gdna_rates = [0.0, 0.1, 0.5, 1.0, 2.0]
    gdna_labels = ["none", "low", "med", "equal", "high"]
    strand_specs = [0.99, 0.50]
    n_rna = args.n_rna

    # Build config object for manifest
    cfg = SimConfig(
        genome=str(fasta_path),
        gtf=str(gtf_path),
        outdir=str(outdir),
        transcript_filter="all",
        simulation=SimulationParams(
            n_rna_fragments=n_rna,
            sim_seed=args.seed,
            frag_mean=250,
            frag_std=50,
            frag_min=80,
            frag_max=800,
            read_length=150,
            error_rate=0.0,
            n_workers=args.n_workers,
        ),
        abundance=AbundanceConfig(
            mode="random",
            seed=args.seed,
            min=1.0,
            max=10000.0,
            frac_expressed=0.5,
        ),
        gdna=GDNASimConfig(
            rates=gdna_rates,
            rate_labels=gdna_labels,
            frag_mean=350,
            frag_std=100,
            frag_min=100,
            frag_max=1000,
        ),
        nrna=NRNAConfig(fracs=[(0.0, 0.0)], frac_labels=["none"]),
        strand_specificities=strand_specs,
        oracle_bam=True,
    )

    # Load transcripts and assign abundances
    print(f"\n  Loading transcripts from {gtf_path}...")
    transcripts = load_transcripts(gtf_path)
    assign_random_abundances(transcripts, cfg.abundance)

    # Write truth abundances
    truth_path = outdir / "truth_abundances.tsv"
    write_truth_abundances(transcripts, truth_path)
    print(f"  Truth abundances: {truth_path}")

    n_expressed = sum(1 for t in transcripts if (t.abundance or 0) > 0)
    print(f"  Expressed: {n_expressed}/{len(transcripts)} transcripts")

    # Run conditions
    total_conditions = len(gdna_rates) * len(strand_specs)
    conditions: list[dict] = []
    cond_num = 0

    for i, gdna_rate in enumerate(gdna_rates):
        gdna_label = gdna_labels[i]
        for strand_spec in strand_specs:
            cond_num += 1
            n_gdna = round(gdna_rate * n_rna)

            cond_name = condition_dir_name(gdna_label, strand_spec, "none")
            cond_dir = outdir / cond_name

            print(
                f"\n  [{cond_num}/{total_conditions}] {cond_name}: "
                f"RNA={n_rna:,} gDNA={n_gdna:,} SS={strand_spec:.2f}",
                flush=True,
            )

            cond_entry = {
                "name": cond_name,
                "gdna_label": gdna_label,
                "gdna_rate": gdna_rate,
                "strand_specificity": strand_spec,
                "nrna_label": "none",
                "nrna_frac_min": 0.0,
                "nrna_frac_max": 0.0,
                "n_rna": n_rna,
                "n_gdna": n_gdna,
                "truth_abundances": "truth_abundances.tsv",
                "fastq_r1": f"{cond_name}/sim_R1.fq.gz",
                "fastq_r2": f"{cond_name}/sim_R2.fq.gz",
            }

            # Check if already done
            if (cond_dir / "sim_R1.fq.gz").exists() and args.skip_existing:
                print("    Output exists, skipping", flush=True)
                cond_entry["oracle_bam"] = f"{cond_name}/sim_oracle.bam"
            else:
                print("    Simulating...", end="", flush=True)
                t_cond = time.monotonic()
                simulator = WholeGenomeSimulator(
                    fasta_path, transcripts, cfg.simulation, cfg.gdna,
                    strand_specificity=strand_spec,
                )
                _, _, bam_path = simulator.simulate_and_write(
                    cond_dir, n_rna, n_gdna,
                    oracle_bam=True, prefix="sim",
                    n_workers=args.n_workers,
                )
                simulator.close()
                cond_entry["oracle_bam"] = (
                    f"{cond_name}/sim_oracle.bam" if bam_path else None
                )
                elapsed_cond = time.monotonic() - t_cond
                print(f" done ({elapsed_cond:.1f}s)", flush=True)

            conditions.append(cond_entry)

    # Write manifest
    write_manifest(outdir, cfg, conditions)

    elapsed_total = time.monotonic() - t0
    print(f"\n{'=' * 70}")
    print(f"SIMULATION COMPLETE")
    print(f"  Output:     {outdir}")
    print(f"  Conditions: {total_conditions}")
    print(f"  Elapsed:    {elapsed_total:.1f}s")
    print(f"{'=' * 70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
