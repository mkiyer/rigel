#!/usr/bin/env python3
"""Generate synthetic mini-genome simulation suites.

Steps:
  1. Generate synthetic genome (10 Mb, ~50 genes, ~250 transcripts)
    2. Simulate reads across strand x gDNA x nRNA-ratio conditions
  3. Each condition produces oracle BAM + paired-end FASTQ

Usage:
    python scripts/sim/simulate_suite.py \\
    --outdir /Users/mkiyer/Downloads/rigel_runs/sim_synthetic
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import logging
import sys
import time
from pathlib import Path

import numpy as np

from rigel.sim.synthetic_genome import (
    GENOME_LENGTH,
    N_GENES,
    REF_NAME,
    SEED,
    generate_genes,
    generate_random_genome,
    inject_splice_sites,
    write_fasta,
    write_gtf,
)
from rigel.sim.abundance import (
    apply_nrna_ratio,
    assign_random_abundances,
    load_transcripts,
    write_truth_abundances,
)
from rigel.sim.config import (
    AbundanceConfig,
    GDNASimConfig,
    NRNAConfig,
    SimConfig,
    SimulationParams,
)
from rigel.sim.manifest import condition_dir_name
from rigel.sim.whole_genome import (
    WholeGenomeSimulator,
    write_manifest,
)

logger = logging.getLogger(__name__)

GDNA_RATES = [0.0, 0.1, 0.5, 1.0, 2.0]
GDNA_LABELS = ["none", "low", "med", "equal", "high"]
STRAND_SPECIFICITIES = [0.99, 0.50]
PROFILE_NRNA = {
    "smoke": ([0.0], ["none"]),
    "full": ([0.0, 0.1, 0.5, 1.0], ["none", "low", "med", "high"]),
}


def stable_seed(base_seed: int, *parts: object) -> int:
    """Derive a reproducible 32-bit seed from a base seed and string parts."""
    text = "\0".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") & 0xFFFF_FFFF


def _parse_csv_floats(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def _parse_csv_labels(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


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
        "--profile", choices=("smoke", "full"), default="smoke",
        help="Condition profile: smoke = 10 no-nRNA baseline conditions, full = 40 conditions"
    )
    parser.add_argument(
        "--n-rna", type=int, default=1_000_000,
        help="Number of mature RNA fragments per condition (default: 1,000,000)"
    )
    parser.add_argument(
        "--nrna-ratios", type=str, default=None,
        help="Comma-separated additive nRNA:mRNA fragment ratios (overrides --profile)"
    )
    parser.add_argument(
        "--nrna-labels", type=str, default=None,
        help="Comma-separated labels for --nrna-ratios (required when overriding ratios)"
    )
    parser.add_argument(
        "--conditions", nargs="*", default=None,
        help="Optional condition names to generate from the selected profile/grid"
    )
    parser.add_argument(
        "-j", "--n-workers", type=int, default=1,
        help="Worker processes for parallel read generation"
    )
    parser.add_argument(
        "--reference-only", action="store_true",
        help="Generate only the synthetic reference FASTA/GTF and exit"
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
    if args.nrna_ratios is None:
        nrna_ratios, nrna_labels = PROFILE_NRNA[args.profile]
    else:
        if args.nrna_labels is None:
            raise ValueError("--nrna-labels is required when --nrna-ratios is set")
        nrna_ratios = _parse_csv_floats(args.nrna_ratios)
        nrna_labels = _parse_csv_labels(args.nrna_labels)
    if len(nrna_labels) != len(nrna_ratios):
        raise ValueError("--nrna-labels must have the same length as --nrna-ratios")
    selected_conditions = set(args.conditions or [])

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
        if iso_counts:
            print(
                f"  Isoforms/gene: {min(iso_counts)}-{max(iso_counts)} "
                f"(mean {np.mean(iso_counts):.1f})"
            )
        else:
            print("  Isoforms/gene: n/a")

    if args.reference_only:
        print(f"\nReference generated at {ref_dir}")
        return 0

    # ══════════════════════════════════════════════════════════════════════
    # Step 2: Run simulation across conditions
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("STEP 2: Simulating reads across conditions")
    print("=" * 70)

    # Configuration
    gdna_rates = GDNA_RATES
    gdna_labels = GDNA_LABELS
    strand_specs = STRAND_SPECIFICITIES
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
            seed=stable_seed(args.seed, "abundance"),
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
        nrna=NRNAConfig(
            mode="additive_ratio",
            ratios=nrna_ratios,
            ratio_labels=nrna_labels,
        ),
        strand_specificities=strand_specs,
        oracle_bam=True,
    )

    # Load transcripts and assign abundances
    print(f"\n  Loading transcripts from {gtf_path}...")
    transcripts = load_transcripts(gtf_path)
    assign_random_abundances(transcripts, cfg.abundance)

    n_expressed = sum(1 for t in transcripts if (t.abundance or 0) > 0)
    print(f"  Expressed: {n_expressed}/{len(transcripts)} transcripts")

    # Run conditions
    all_condition_names = [
        condition_dir_name(gdna_label, strand_spec, nrna_label)
        for nrna_label in nrna_labels
        for gdna_label in gdna_labels
        for strand_spec in strand_specs
    ]
    unknown_conditions = selected_conditions.difference(all_condition_names)
    if unknown_conditions:
        unknown = ", ".join(sorted(unknown_conditions))
        raise ValueError(f"Unknown condition(s) for selected profile/grid: {unknown}")
    total_conditions = (
        len(selected_conditions) if selected_conditions else len(all_condition_names)
    )
    conditions: list[dict] = []
    cond_num = 0

    for nrna_ratio, nrna_label in zip(nrna_ratios, nrna_labels):
        cond_transcripts = copy.deepcopy(transcripts)
        apply_nrna_ratio(cond_transcripts, nrna_ratio)
        truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        truth_path = outdir / truth_name
        write_truth_abundances(cond_transcripts, truth_path)
        print(f"  Truth abundances ({nrna_label}): {truth_path}")

        for i, gdna_rate in enumerate(gdna_rates):
            gdna_label = gdna_labels[i]
            for strand_spec in strand_specs:
                n_mrna = n_rna
                n_nrna = round(n_mrna * nrna_ratio)
                n_rna_total = n_mrna + n_nrna
                n_gdna = round(gdna_rate * n_mrna)

                cond_name = condition_dir_name(gdna_label, strand_spec, nrna_label)
                if selected_conditions and cond_name not in selected_conditions:
                    continue
                cond_num += 1
                condition_seed = stable_seed(args.seed, cond_name)
                cond_dir = outdir / cond_name

                print(
                    f"\n  [{cond_num}/{total_conditions}] {cond_name}: "
                    f"mRNA={n_mrna:,} nRNA={n_nrna:,} gDNA={n_gdna:,} "
                    f"SS={strand_spec:.2f}",
                    flush=True,
                )

                cond_entry = {
                    "name": cond_name,
                    "gdna_label": gdna_label,
                    "gdna_rate": gdna_rate,
                    "strand_specificity": strand_spec,
                    "nrna_label": nrna_label,
                    "nrna_mode": "additive_ratio",
                    "nrna_ratio": nrna_ratio,
                    "n_mrna": n_mrna,
                    "n_nrna": n_nrna,
                    "n_rna": n_rna_total,
                    "n_gdna": n_gdna,
                    "n_total": n_rna_total + n_gdna,
                    "seed": condition_seed,
                    "truth_abundances": truth_name,
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
                    cond_sim = copy.deepcopy(cfg.simulation)
                    cond_sim.sim_seed = condition_seed
                    simulator = WholeGenomeSimulator(
                        fasta_path, cond_transcripts, cond_sim, cfg.gdna,
                        strand_specificity=strand_spec,
                    )
                    _, _, bam_path = simulator.simulate_and_write(
                        cond_dir, n_rna_total, n_gdna,
                        oracle_bam=True, prefix="sim",
                        n_mrna=n_mrna,
                        n_nrna=n_nrna,
                        n_workers=cond_sim.n_workers,
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
    print("SIMULATION COMPLETE")
    print(f"  Output:     {outdir}")
    print(f"  Conditions: {total_conditions}")
    print(f"  Elapsed:    {elapsed_total:.1f}s")
    print(f"{'=' * 70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
